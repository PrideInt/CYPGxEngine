import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pgx.core.engine import (
    find_recommendations,
    list_drugs,
    search_diplotypes,
    generate_report,
    lookup_gene,
    normalize_diplotype,
    open_database
)
from pgx.models import is_resolved
from tests.make_fixture import build

class TestNormalizeDiplotype(unittest.TestCase):
    def test_orders_alleles_numerically(self):
        self.assertEqual(normalize_diplotype("*2/*1"), "*1/*2")
        self.assertEqual(normalize_diplotype("*10/*2"), "*2/*10")

    def test_tolerates_whitespace(self):
        self.assertEqual(normalize_diplotype(" *1 / *2 "), "*1/*2")

    def test_preserves_duplication_notation(self):
        self.assertEqual(normalize_diplotype("*1x2/*4"), "*1x2/*4")

    def test_rejects_malformed(self):
        with self.assertRaises(ValueError):
            normalize_diplotype("*1")
        with self.assertRaises(ValueError):
            normalize_diplotype("*1/*2/*3")

class TestEngine(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.mkdtemp()
        cls.database_path = build(str(Path(cls.directory) / "fixture.db"))

    def setUp(self):
        self.connection = open_database(self.database_path)

    def tearDown(self):
        self.connection.close()

    # Phenotype lookup
    def test_resolves_phenotype(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*2")
        self.assertTrue(is_resolved(result))
        self.assertEqual(result["phenotype"], "Intermediate Metabolizer")

    def test_resolves_regardless_of_allele_order(self):
        result = lookup_gene(self.connection, "CYP2C19", "*2/*1")
        self.assertEqual(result["phenotype"], "Intermediate Metabolizer")

    def test_case_insensitive_gene(self):
        self.assertTrue(is_resolved(lookup_gene(self.connection, "cyp2c19", "*1/*2")))

    def test_unknown_diplotype_is_an_error_not_a_guess(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*999")
        self.assertFalse(is_resolved(result))
        self.assertIsNone(result["phenotype"])
        self.assertIn("not found", result["error"])

    def test_unknown_gene_reports_distinctly(self):
        result = lookup_gene(self.connection, "NUDT15", "*1/*3")
        self.assertIn("not present", result["error"])

    # CYP2C19 phenotype categories
    def test_all_eight_cyp2c19_phenotypes_present(self):
        rows = self.connection.execute(
            "SELECT DISTINCT phenotype FROM diplotype WHERE gene = 'CYP2C19'"
        ).fetchall()
        self.assertEqual(8, len({row["phenotype"] for row in rows}))

    def test_likely_intermediate_is_distinct_from_intermediate(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*9")
        self.assertEqual(result["phenotype"], "Likely Intermediate Metabolizer")

    def test_rapid_is_distinct_from_ultrarapid(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*17")
        self.assertEqual(result["phenotype"], "Rapid Metabolizer")

    # Indeterminate is a result, not an error
    def test_indeterminate_resolves(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*12")
        self.assertTrue(is_resolved(result))
        self.assertEqual(result["phenotype"], "Indeterminate")
        self.assertIsNone(result["error"])

    def test_indeterminate_yields_a_no_recommendation_row(self):
        results = [lookup_gene(self.connection, "CYP2C19", "*1/*12")]
        
        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["escitalopram"]
        )
        self.assertEqual(1, len(recommendations))
        self.assertEqual("No Recommendation", recommendations[0]["classification"])

    # CYP2D6 activity score banding
    def test_activity_score_carried(self):
        result = lookup_gene(self.connection, "CYP2D6", "*1/*4")
        self.assertEqual("1.0", result["activity_score"])

    def test_activity_score_one_is_intermediate_not_normal(self):
        result = lookup_gene(self.connection, "CYP2D6", "*1/*4")
        self.assertEqual("Intermediate Metabolizer", result["phenotype"])

    def test_normal_metabolizer_floor_boundary(self):
        result = lookup_gene(self.connection, "CYP2D6", "*1/*10")
        self.assertEqual("1.25", result["activity_score"])
        self.assertEqual("Normal Metabolizer", result["phenotype"])

    def test_cyp2d6_indeterminate_keys_on_na(self):
        result = lookup_gene(self.connection, "CYP2D6", "*1/*22")
        self.assertEqual({"CYP2D6": "n/a"}, result["lookup_key"])

    def test_unknown_copy_number_keys_on_a_bucket(self):
        result = lookup_gene(self.connection, "CYP2D6", "*1/*1xN")
        self.assertEqual({"CYP2D6": "\u22653.0"}, result["lookup_key"])
        self.assertEqual("Ultrarapid Metabolizer", result["phenotype"])

    def test_bucketed_score_still_matches_recommendations(self):
        results = [lookup_gene(self.connection, "CYP2D6", "*1/*1xN")]

        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["paroxetine"]
        )
        self.assertEqual(1, len(recommendations))
        self.assertEqual("Moderate", recommendations[0]["classification"])

    def test_cyp2d6_joins_on_activity_score(self):
        result = lookup_gene(self.connection, "CYP2D6", "*4/*4")
        self.assertEqual({"CYP2D6": "0.0"}, result["lookup_key"])

    # The CYP2D6 drugs must not be interchangeable
    def test_ultrarapid_differs_across_the_cyp2d6_drugs(self):
        results = [lookup_gene(self.connection, "CYP2D6", "*1/*1x2")]

        classifications = {
            recommendation["drug"]: recommendation["classification"]
            for recommendation in find_recommendations(self.connection, results)
        }
        self.assertEqual("Moderate", classifications["paroxetine"])
        self.assertEqual("No Recommendation", classifications["fluvoxamine"])
        self.assertEqual("No Recommendation", classifications["venlafaxine"])
        self.assertEqual("Optional", classifications["vortioxetine"])

    def test_poor_metabolizer_differs_across_the_cyp2d6_drugs(self):
        results = [lookup_gene(self.connection, "CYP2D6", "*4/*4")]

        texts = {
            recommendation["drug"]: recommendation["recommendation"]
            for recommendation in find_recommendations(self.connection, results)
        }
        self.assertEqual(4, len(set(texts.values())))

    # Specificity rule
    def test_single_gene_patient_gets_single_gene_recommendation(self):
        results = [lookup_gene(self.connection, "CYP2C19", "*1/*1")]

        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"]
        )
        self.assertEqual(1, len(recommendations))
        self.assertEqual(["CYP2C19"], recommendations[0]["genes"])

    def test_dual_typed_patient_gets_the_combined_cell(self):
        results = [
            lookup_gene(self.connection, "CYP2C19", "*1/*1"),
            lookup_gene(self.connection, "CYP2B6", "*6/*6")
        ]
        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"]
        )
        self.assertEqual(1, len(recommendations))
        self.assertEqual(["CYP2B6", "CYP2C19"], recommendations[0]["genes"])

    def test_combined_can_differ_from_either_single_gene(self):
        combined = find_recommendations(
            self.connection,
            [
                lookup_gene(self.connection, "CYP2C19", "*1/*1"),
                lookup_gene(self.connection, "CYP2B6", "*6/*6")
            ],
            drugs=["sertraline"]
        )[0]
        single = find_recommendations(
            self.connection,
            [lookup_gene(self.connection, "CYP2C19", "*1/*1")],
            drugs=["sertraline"]
        )[0]
        self.assertNotEqual(combined["recommendation"], single["recommendation"])

    def test_grid_cells_are_mutually_exclusive(self):
        results = [
            lookup_gene(self.connection, "CYP2C19", "*1/*1"),
            lookup_gene(self.connection, "CYP2B6", "*6/*6")
        ]
        every = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"],
            specific_only=False
        )
        self.assertEqual(1, len(every))

    def test_untyped_gene_matches_a_no_result_cell(self):
        results = [lookup_gene(self.connection, "CYP2C19", "*1/*1")]

        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"]
        )
        self.assertEqual(1, len(recommendations))
        self.assertEqual("Strong", recommendations[0]["classification"])

    def test_no_result_gene_is_omitted_from_the_reported_genes(self):
        results = [lookup_gene(self.connection, "CYP2C19", "*1/*1")]

        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"]
        )
        self.assertEqual(["CYP2C19"], recommendations[0]["genes"])

    def test_typing_the_second_gene_changes_the_answer(self):
        single = find_recommendations(
            self.connection,
            [lookup_gene(self.connection, "CYP2C19", "*1/*1")],
            drugs=["sertraline"]
        )[0]
        both = find_recommendations(
            self.connection,
            [
                lookup_gene(self.connection, "CYP2C19", "*1/*1"),
                lookup_gene(self.connection, "CYP2B6", "*6/*6")
            ],
            drugs=["sertraline"]
        )[0]
        self.assertEqual(["CYP2C19"], single["genes"])
        self.assertEqual(["CYP2B6", "CYP2C19"], both["genes"])
        self.assertNotEqual(single["recommendation"], both["recommendation"])

    def test_partial_multigene_resolves_through_the_no_result_cell(self):
        results = [lookup_gene(self.connection, "CYP2B6", "*6/*6")]

        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"]
        )
        self.assertEqual(1, len(recommendations))
        self.assertEqual(["CYP2B6"], recommendations[0]["genes"])

    def test_a_typed_gene_never_matches_a_no_result_cell(self):
        results = [
            lookup_gene(self.connection, "CYP2C19", "*1/*1"),
            lookup_gene(self.connection, "CYP2B6", "*6/*6")
        ]
        every = find_recommendations(
            self.connection,
            results,
            drugs=["sertraline"],
            specific_only=False
        )
        for recommendation in every:
            self.assertEqual(
                ["CYP2B6", "CYP2C19"],
                recommendation["genes"]
            )

    def test_specificity_does_not_affect_unrelated_drugs(self):
        results = [
            lookup_gene(self.connection, "CYP2C19", "*1/*1"),
            lookup_gene(self.connection, "CYP2B6", "*6/*6"),
            lookup_gene(self.connection, "CYP2D6", "*4/*4")
        ]
        drugs = {
            recommendation["drug"]
            for recommendation in find_recommendations(self.connection, results)
        }
        self.assertEqual(
            {
                "escitalopram",
                "sertraline",
                "paroxetine",
                "fluvoxamine",
                "venlafaxine",
                "vortioxetine"
            },
            drugs
        )

    def test_unresolved_gene_contributes_nothing(self):
        results = [lookup_gene(self.connection, "CYP2C19", "*1/*999")]
        self.assertEqual([], find_recommendations(self.connection, results))

    def test_cpic_level_is_for_the_recommendations_own_gene(self):
        results = [lookup_gene(self.connection, "CYP2D6", "*4/*4")]
        recommendations = find_recommendations(
            self.connection,
            results,
            drugs=["paroxetine"]
        )
        self.assertEqual("A", recommendations[0]["cpic_level"])

    # Report assembly
    # Fields CPIC publishes alongside the phenotype
    def test_ehr_priority_carried(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*1")
        self.assertEqual("Normal/Routine/Low Risk", result["ehr_priority"])

    def test_abnormal_priority_carried(self):
        result = lookup_gene(self.connection, "CYP2D6", "*4/*4")
        self.assertEqual("Abnormal/Priority/High Risk", result["ehr_priority"])

    def test_allele_functions_carried_in_order(self):
        result = lookup_gene(self.connection, "CYP2D6", "*1/*10")
        self.assertEqual(["Normal function", "Decreased function"], result["allele_functions"])

    def test_consultation_text_carried(self):
        result = lookup_gene(self.connection, "CYP2B6", "*6/*6")
        self.assertIsNotNone(result["consultation_text"])

    def test_absent_fields_are_empty_not_missing(self):
        """A diplotype CPIC publishes no banding for still reports the keys."""
        result = lookup_gene(self.connection, "CYP2C19", "*1/*17")
        self.assertIsNone(result["ehr_priority"])
        self.assertEqual([], result["allele_functions"])

    def test_unresolved_gene_still_reports_the_new_keys(self):
        result = lookup_gene(self.connection, "CYP2C19", "*1/*999")

        for key in ("ehr_priority", "consultation_text", "allele_functions"):
            self.assertIn(key, result)

    # A knowledge base built before this schema must be refused, not misread
    def test_rejects_an_older_knowledge_base(self):
        import sqlite3

        stale = str(Path(self.directory) / "stale.db")
        connection = sqlite3.connect(stale)
        connection.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
        connection.execute("INSERT INTO meta VALUES ('cpic_db', 'FIXTURE')")
        connection.commit()
        connection.close()

        with self.assertRaises(ValueError) as caught:
            open_database(stale)

        self.assertIn("pre-1.1", str(caught.exception))

    # Browsing the knowledge base
    def test_search_narrows_by_substring(self):
        result = search_diplotypes(self.connection, "CYP2C19", "*17")
        self.assertTrue(result["total"] >= 2)

        for entry in result["diplotypes"]:
            self.assertIn("*17", entry["diplotype"])

    def test_search_without_a_query_returns_the_gene(self):
        every = search_diplotypes(self.connection, "CYP2C19", limit=200)
        narrowed = search_diplotypes(self.connection, "CYP2C19", "*17", limit=200)
        self.assertGreater(every["total"], narrowed["total"])

    def test_search_reports_the_total_beyond_the_limit(self):
        result = search_diplotypes(self.connection, "CYP2D6", limit=2)
        self.assertEqual(2, len(result["diplotypes"]))
        self.assertGreater(result["total"], 2)

    def test_search_is_case_insensitive_on_the_gene(self):
        self.assertEqual(
            search_diplotypes(self.connection, "cyp2c19")["total"],
            search_diplotypes(self.connection, "CYP2C19")["total"]
        )

    def test_search_treats_wildcards_as_literal_text(self):
        """A bare % must not match everything."""
        self.assertEqual(0, search_diplotypes(self.connection, "CYP2C19", "%")["total"])

    def test_search_on_an_unknown_gene_is_empty_not_an_error(self):
        result = search_diplotypes(self.connection, "NUDT15")
        self.assertEqual(0, result["total"])
        self.assertEqual([], result["diplotypes"])

    def test_search_carries_the_priority_band(self):
        result = search_diplotypes(self.connection, "CYP2D6", "*4/*4")
        self.assertEqual("Abnormal/Priority/High Risk", result["diplotypes"][0]["ehr_priority"])

    # Listing the drugs this build can actually answer for
    def test_lists_only_answerable_drugs_by_default(self):
        names = {entry["drug"] for entry in list_drugs(self.connection)}
        self.assertIn("sertraline", names)
        self.assertNotIn("warfarin", names)

    def test_unanswerable_drugs_are_listed_when_asked_for(self):
        entries = {
            entry["drug"]: entry
            for entry in list_drugs(self.connection, answerable_only=False)
        }
        self.assertIn("warfarin", entries)
        self.assertFalse(entries["warfarin"]["answerable"])
        self.assertEqual(["CYP2C9"], entries["warfarin"]["genes"])

    def test_a_two_gene_drug_reports_both_genes(self):
        entries = {entry["drug"]: entry for entry in list_drugs(self.connection)}
        self.assertEqual(["CYP2B6", "CYP2C19"], entries["sertraline"]["genes"])

    def test_drug_listing_carries_the_evidence_level(self):
        entries = {entry["drug"]: entry for entry in list_drugs(self.connection)}
        self.assertEqual("A", entries["paroxetine"]["cpic_level"])

    def test_report_is_serializable(self):
        report = generate_report(
            self.connection,
            {"CYP2C19": "*2/*2", "CYP2D6": "*4/*4"}
        )
        self.assertEqual("1.1", report["schema_version"])
        self.assertEqual(2, len(report["genes"]))
        self.assertIn("Educational use only", report["disclaimer"])

    def test_report_round_trips_through_json(self):
        report = generate_report(self.connection, {"CYP2C19": "*1/*2"})
        self.assertEqual(report, json.loads(json.dumps(report)))

    def test_warnings_surface_bad_input(self):
        report = generate_report(self.connection, {"CYP2C19": "garbage"})
        self.assertTrue(report["warnings"])
        self.assertEqual([], report["drugs"])

    def test_warns_on_undeliverable_drug_request(self):
        report = generate_report(
            self.connection,
            {"CYP2C19": "*1/*1"},
            drugs=["clopidogrel"]
        )
        self.assertTrue(any("clopidogrel" in warning for warning in report["warnings"]))

    def test_provenance_present(self):
        report = generate_report(self.connection, {"CYP2C19": "*1/*1"})
        self.assertEqual("FIXTURE", report["sources"]["cpic_db"])

if __name__ == "__main__":
    unittest.main(verbosity=2)
