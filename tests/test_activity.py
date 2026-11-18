import sys
import unittest
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pgx.activity import (
    INDETERMINATE,
    INTERMEDIATE,
    NORMAL,
    POOR,
    ULTRARAPID,
    band_activity_score,
    parse_activity_score,
    phenotype_from_values,
    total_activity_score
)

class TestParseActivityScore(unittest.TestCase):
    def test_exact_score(self):
        self.assertEqual((Decimal("1.25"), False), parse_activity_score("1.25"))

    def test_bucketed_score_is_a_lower_bound(self):
        self.assertEqual((Decimal("3.0"), True), parse_activity_score("≥3.0"))

    def test_absent_score(self):
        for text in ("n/a", "N/A", "", "   ", None, "No Result"):
            self.assertEqual((None, False), parse_activity_score(text))

    def test_unparseable_score(self):
        self.assertEqual((None, False), parse_activity_score("Normal Metabolizer"))

class TestBandActivityScore(unittest.TestCase):
    # Band boundaries from the 2020 standardization
    def test_zero_is_poor(self):
        self.assertEqual(POOR, band_activity_score(Decimal("0")))

    def test_just_above_zero_is_intermediate(self):
        self.assertEqual(INTERMEDIATE, band_activity_score(Decimal("0.25")))

    def test_one_is_intermediate_not_normal(self):
        self.assertEqual(INTERMEDIATE, band_activity_score(Decimal("1.0")))

    def test_normal_floor_is_inclusive(self):
        self.assertEqual(NORMAL, band_activity_score(Decimal("1.25")))

    def test_normal_ceiling_is_inclusive(self):
        self.assertEqual(NORMAL, band_activity_score(Decimal("2.25")))

    def test_just_above_ceiling_is_ultrarapid(self):
        self.assertEqual(ULTRARAPID, band_activity_score(Decimal("2.5")))

    def test_missing_score_is_indeterminate(self):
        self.assertEqual(INDETERMINATE, band_activity_score(None))

    # A lower bound only decides the phenotype above the ultrarapid threshold
    def test_lower_bound_above_threshold_is_ultrarapid(self):
        self.assertEqual(
            ULTRARAPID,
            band_activity_score(Decimal("3.0"), at_least=True)
        )

    def test_lower_bound_below_threshold_cannot_be_decided(self):
        self.assertEqual(
            INDETERMINATE,
            band_activity_score(Decimal("1.0"), at_least=True)
        )

class TestTotalActivityScore(unittest.TestCase):
    def test_sums_two_values(self):
        self.assertEqual(Decimal("1.25"), total_activity_score("1.0", "0.25"))

    def test_sums_to_zero(self):
        self.assertEqual(Decimal("0"), total_activity_score("0", "0"))

    def test_missing_value_yields_nothing(self):
        self.assertIsNone(total_activity_score("1.0", None))
        self.assertIsNone(total_activity_score(None, None))
        self.assertIsNone(total_activity_score("1.0", ""))

    def test_decimal_arithmetic_is_exact(self):
        """0.1 + 0.2 in float is not 0.3, which would misband a boundary."""
        self.assertEqual(Decimal("0.3"), total_activity_score("0.1", "0.2"))

    # Worked examples from the guideline
    def test_normal_and_no_function_is_intermediate(self):
        self.assertEqual(INTERMEDIATE, phenotype_from_values("1.0", "0"))

    def test_normal_and_star_ten_is_the_normal_floor(self):
        self.assertEqual(NORMAL, phenotype_from_values("1.0", "0.25"))

    def test_two_normal_alleles_is_normal(self):
        self.assertEqual(NORMAL, phenotype_from_values("1.0", "1.0"))

    def test_two_no_function_alleles_is_poor(self):
        self.assertEqual(POOR, phenotype_from_values("0", "0"))

if __name__ == "__main__":
    unittest.main(verbosity=2)
