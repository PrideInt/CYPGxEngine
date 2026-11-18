import sqlite3
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pgx.frequency import (
    allele_prior,
    brier,
    candidates_for,
    gene_model,
    log_loss,
    marginal,
    pooled_prior,
    posterior
)
from pgx.panel import is_structural, mask
from pgx.query.cpic import SCHEMA

GENE = "TESTG"

NORMAL = "Normal Metabolizer"
INTERMEDIATE = "Intermediate Metabolizer"
POOR = "Poor Metabolizer"

# A toy gene small enough to do by hand. *1 normal, *2 and *3 no function,
# *4 present in diplotypes but never characterised, and one structural
# diplotype to confirm it is kept out of the allele pool.
DIPLOTYPES = [
    ("*1/*1", NORMAL),
    ("*1/*2", INTERMEDIATE),
    ("*1/*3", INTERMEDIATE),
    ("*1/*4", INTERMEDIATE),
    ("*2/*2", POOR),
    ("*2/*3", POOR),
    ("*2/*4", POOR),
    ("*3/*3", POOR),
    ("*3/*4", POOR),
    ("*4/*4", POOR),
    ("*1/*2x2", POOR)
]

# Test:  *1 0.6, *2 0.3, *3 0.1  — sums to 1, so no renormalisation to hide
#        arithmetic errors behind.
# Zeros: *3 recorded as an observed zero, which must be excluded from the
#        prior while still counting toward published mass.
FREQUENCIES = [
    ("Test", "*1", 0.6), ("Test", "*2", 0.3), ("Test", "*3", 0.1),
    ("Zeros", "*1", 0.7), ("Zeros", "*2", 0.3), ("Zeros", "*3", 0.0)
]

def toy_connection():
    connection = sqlite3.connect(":memory:")
    connection.executescript(SCHEMA)

    for diplotype, phenotype in DIPLOTYPES:
        connection.execute(
            "INSERT INTO diplotype (gene, diplotype, phenotype) VALUES (?,?,?)",
            (GENE, diplotype, phenotype)
        )

    for population, allele, frequency in FREQUENCIES:
        connection.execute(
            "INSERT INTO allele_frequency VALUES (?,?,?,?)",
            (GENE, allele, population, frequency)
        )

    connection.commit()
    return connection

class FrequencyTestCase(unittest.TestCase):

    def setUp(self):
        self.connection = toy_connection()
        self.model = gene_model(self.connection, GENE)
        self.prior = allele_prior(self.connection, GENE, "Test", self.model["alleles"])

    def tearDown(self):
        self.connection.close()

class TestGeneModel(FrequencyTestCase):
    def test_collects_alleles_from_diplotypes(self):
        self.assertEqual(self.model["alleles"], ["*1", "*2", "*3", "*4"])

    def test_excludes_structural_alleles(self):
        self.assertNotIn("*2x2", self.model["alleles"])

    def test_keeps_structural_diplotypes_resolvable(self):
        # Excluded from the allele pool, but still a row the engine can look up.
        self.assertEqual(self.model["phenotypes"]["*1/*2x2"], POOR)


class TestAllelePrior(FrequencyTestCase):
    def test_frequencies_are_a_distribution(self):
        self.assertAlmostEqual(sum(self.prior["frequencies"].values()), 1.0, places=9)

    def test_uncharacterised_allele_is_reported_missing(self):
        self.assertIn("*4", self.prior["missing"])
        self.assertNotIn("*4", self.prior["frequencies"])

    def test_observed_zero_is_excluded_but_counted(self):
        prior = allele_prior(self.connection, GENE, "Zeros", self.model["alleles"])

        self.assertNotIn("*3", prior["frequencies"])
        self.assertIn("*3", prior["missing"])
        self.assertAlmostEqual(prior["coverage"], 1.0, places=9)

    def test_coverage_is_the_retained_share(self):
        self.assertAlmostEqual(self.prior["coverage"], 1.0, places=9)

    def test_pooled_averages_across_populations(self):
        pooled = pooled_prior(self.connection, GENE, self.model["alleles"], ["Test", "Zeros"])

        # Test contributes 0.6/0.3/0.1, Zeros 0.7/0.3/— ; unweighted mean.
        self.assertAlmostEqual(pooled["frequencies"]["*1"], 0.65, places=9)
        self.assertAlmostEqual(pooled["frequencies"]["*2"], 0.30, places=9)
        self.assertAlmostEqual(pooled["frequencies"]["*3"], 0.05, places=9)

class TestMarginal(FrequencyTestCase):
    def test_hardy_weinberg_by_hand(self):
        # *1/*1 .36 normal; *1/*2 .36 and *1/*3 .12 intermediate;
        # *2/*2 .09, *2/*3 .06, *3/*3 .01 poor.
        distribution = marginal(self.model, self.prior)

        self.assertAlmostEqual(distribution[NORMAL], 0.36, places=9)
        self.assertAlmostEqual(distribution[INTERMEDIATE], 0.48, places=9)
        self.assertAlmostEqual(distribution[POOR], 0.16, places=9)

    def test_sums_to_one(self):
        self.assertAlmostEqual(
            sum(marginal(self.model, self.prior).values()), 1.0, places=9)


class TestPosterior(FrequencyTestCase):
    def test_empty_panel_equals_the_marginal(self):
        # With nothing interrogated every diplotype reports *1/*1, so the
        # posterior has learned nothing and must be the prior marginal.
        expected = marginal(self.model, self.prior)
        actual = posterior(self.model, self.prior, "*1/*1", [])

        for phenotype in set(expected) | set(actual):
            self.assertAlmostEqual(actual.get(phenotype, 0.0), expected.get(phenotype, 0.0), places=12)

    def test_full_panel_is_a_point_mass(self):
        full = ["*2", "*3"]

        self.assertEqual(posterior(self.model, self.prior, "*1/*1", full),
                         {NORMAL: 1.0})
        self.assertEqual(posterior(self.model, self.prior, "*2/*2", full),
                         {POOR: 1.0})

    def test_partial_panel_by_hand(self):
        # Panel carries *2 only. An observed *1/*1 could be *1/*1 (.36),
        # *1/*3 (.12) or *3/*3 (.01): total .49.
        distribution = posterior(self.model, self.prior, "*1/*1", ["*2"])

        self.assertAlmostEqual(distribution[NORMAL], 0.36 / 0.49, places=9)
        self.assertAlmostEqual(distribution[INTERMEDIATE], 0.12 / 0.49, places=9)
        self.assertAlmostEqual(distribution[POOR], 0.01 / 0.49, places=9)

    def test_identified_allele_narrows_one_side_only(self):
        # Observed *1/*2 with *2 carried: *2 is certain, the other side is
        # *1 (.36) or *3 (.06).
        distribution = posterior(self.model, self.prior, "*1/*2", ["*2"])

        self.assertAlmostEqual(distribution[INTERMEDIATE], 0.36 / 0.42, places=9)
        self.assertAlmostEqual(distribution[POOR], 0.06 / 0.42, places=9)

    def test_sums_to_one(self):
        for panel in ([], ["*2"], ["*2", "*3"]):
            for observed in ("*1/*1", "*1/*2", "*2/*2"):
                distribution = posterior(self.model, self.prior, observed, panel)

                if distribution:
                    self.assertAlmostEqual(sum(distribution.values()), 1.0, places=9)

    def test_allele_without_a_prior_yields_no_posterior(self):
        # *4 is on the panel and was reported, but was never characterised.
        # An empty distribution, not a confident wrong one.
        self.assertEqual(posterior(self.model, self.prior, "*4/*4", ["*4"]), {})


class TestCandidates(FrequencyTestCase):
    def test_carried_allele_names_itself(self):
        self.assertEqual(
            candidates_for("*2", ["*1", "*2", "*3"], ["*2"]), ["*2"])

    def test_reference_stands_for_everything_uninterrogated(self):
        self.assertEqual(
            candidates_for("*1", ["*1", "*2", "*3"], ["*2"]), ["*1", "*3"])

    def test_carried_allele_without_a_prior_names_nothing(self):
        self.assertEqual(candidates_for("*4", ["*1", "*2"], ["*4"]), [])

    def test_uncarried_non_reference_is_unreachable(self):
        self.assertEqual(candidates_for("*3", ["*1", "*2", "*3"], ["*2"]), [])


class TestMask(unittest.TestCase):
    def test_uncarried_alleles_become_reference(self):
        self.assertEqual(mask("*2/*3", ["*2"]), "*1/*2")

    def test_empty_panel_reports_reference(self):
        self.assertEqual(mask("*2/*3", []), "*1/*1")

    def test_carried_alleles_survive(self):
        self.assertEqual(mask("*2/*2", ["*2"]), "*2/*2")

    def test_result_is_normalised(self):
        self.assertEqual(mask("*3/*2", ["*2"]), "*1/*2")

    def test_structural_detection(self):
        self.assertTrue(is_structural("*2x2"))
        self.assertTrue(is_structural("*36+*10"))
        self.assertFalse(is_structural("*2"))


class TestScoring(unittest.TestCase):
    def test_brier_rewards_a_correct_point_call(self):
        self.assertAlmostEqual(brier({NORMAL: 1.0}, NORMAL), 0.0, places=9)

    def test_brier_caps_a_confidently_wrong_call(self):
        # The reason Brier is the headline and log loss is not: wrong and
        # certain costs 2.0, not infinity.
        self.assertAlmostEqual(brier({NORMAL: 1.0}, POOR), 2.0, places=9)

    def test_brier_prefers_hedging_when_wrong(self):
        hedged = brier({NORMAL: 0.5, POOR: 0.5}, POOR)

        self.assertLess(hedged, brier({NORMAL: 1.0}, POOR))
        self.assertAlmostEqual(hedged, 0.5, places=9)

    def test_log_loss_is_floored(self):
        self.assertAlmostEqual(log_loss({NORMAL: 1.0}, NORMAL), 0.0, places=9)
        self.assertLess(log_loss({NORMAL: 1.0}, POOR), 14.0)

if __name__ == "__main__":
    unittest.main()
