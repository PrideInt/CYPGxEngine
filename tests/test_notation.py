import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pgx.notation import (
    ALTERNATE_DIPLOTYPE,
    AMBIGUOUS_CALL,
    NOT_A_DIPLOTYPE,
    PLAIN,
    SUBALLELE_ANNOTATION,
    UNCERTAIN_ALLELE,
    UNKNOWN_COPY_NUMBER,
    expand_unknown_copies,
    normalize_getrm
)

class TestNormalizeGetRm(unittest.TestCase):
    # Values needing no interpretation
    def test_plain_value_passes_through(self):
        result = normalize_getrm("*1/*6")
        self.assertEqual(["*1/*6"], result["candidates"])
        self.assertEqual(PLAIN, result["note"])

    def test_tandem_arrangement_is_preserved(self):
        result = normalize_getrm("*36x2+*10/*52")
        self.assertEqual(["*36x2+*10/*52"], result["candidates"])
        self.assertEqual(PLAIN, result["note"])

    def test_lowercase_duplication_is_preserved(self):
        self.assertEqual(["*2x2/*22"], normalize_getrm("*2x2/*22")["candidates"])

    # Values that are not diplotypes at all
    def test_placeholder_text_yields_nothing(self):
        for raw in ("no consensus", "duplication", "", None, "   "):
            result = normalize_getrm(raw)
            self.assertEqual([], result["candidates"])
            self.assertEqual(NOT_A_DIPLOTYPE, result["note"])

    def test_single_allele_is_rejected(self):
        result = normalize_getrm("*1(*5) or *1(*22)")
        self.assertEqual([], result["candidates"])
        self.assertEqual(NOT_A_DIPLOTYPE, result["note"])

    # Parentheses carry three different meanings
    def test_suballele_after_an_allele_is_stripped(self):
        result = normalize_getrm("*2(*35)/*4")
        self.assertEqual(["*2/*4"], result["candidates"])
        self.assertEqual(SUBALLELE_ANNOTATION, result["note"])

    def test_suballele_with_a_space_is_stripped(self):
        self.assertEqual(["*2/*9"], normalize_getrm("*2 (*35)/*9")["candidates"])

    def test_suballele_on_both_alleles_is_stripped(self):
        self.assertEqual(["*1/*1"], normalize_getrm("*1(*5)/*1(*27)")["candidates"])

    def test_standalone_parenthesised_allele_is_unwrapped(self):
        result = normalize_getrm("*2/(*4)")
        self.assertEqual(["*2/*4"], result["candidates"])
        self.assertEqual(UNCERTAIN_ALLELE, result["note"])

    def test_parenthesised_alternate_diplotype_is_offered_second(self):
        result = normalize_getrm("*1/*1 (*1/*17)")
        self.assertEqual(["*1/*1", "*1/*17"], result["candidates"])
        self.assertEqual(ALTERNATE_DIPLOTYPE, result["note"])

    def test_trailing_commentary_is_discarded(self):
        self.assertEqual(
            ["*2x2/*4x2"],
            normalize_getrm("*2x2/*4x2 (+hybrid)")["candidates"]
        )

    # Values offering more than one reading
    def test_explicit_or_yields_both_readings(self):
        result = normalize_getrm("*1/*4 or *4/*17")
        self.assertEqual(["*1/*4", "*4/*17"], result["candidates"])
        self.assertEqual(AMBIGUOUS_CALL, result["note"])

    def test_semicolon_expands_in_place(self):
        result = normalize_getrm("*1/(*2;*10)")
        self.assertEqual(["*1/*2", "*1/*10"], result["candidates"])
        self.assertEqual(AMBIGUOUS_CALL, result["note"])

    def test_semicolon_expands_on_the_left(self):
        self.assertEqual(
            ["*2/*6", "*10/*6"],
            normalize_getrm("(*2;*10)/*6")["candidates"]
        )

    # Unknown copy number
    def test_unknown_copy_number_is_flagged(self):
        result = normalize_getrm("*1XN/*2")
        self.assertEqual(UNKNOWN_COPY_NUMBER, result["note"])

    def test_uppercase_duplication_is_lowercased(self):
        self.assertEqual(["*1x2/*2"], normalize_getrm("*1X2/*2")["candidates"])

    def test_duplicate_candidates_are_collapsed(self):
        result = normalize_getrm("*1/*1 (*1/*1)")
        self.assertEqual(["*1/*1"], result["candidates"])


class TestExpandUnknownCopies(unittest.TestCase):
    def test_unknown_copies_offer_both_cpic_forms(self):
        self.assertEqual(
            ["*1x2/*2", "*1x≥3/*2"],
            expand_unknown_copies("*1xN/*2")
        )

    def test_uppercase_unknown_copies_expand(self):
        self.assertEqual(
            ["*2/*2x2", "*2/*2x≥3"],
            expand_unknown_copies("*2/*2XN")
        )

    def test_known_diplotype_is_untouched(self):
        self.assertEqual(["*1/*2"], expand_unknown_copies("*1/*2"))

if __name__ == "__main__":
    unittest.main(verbosity=2)
