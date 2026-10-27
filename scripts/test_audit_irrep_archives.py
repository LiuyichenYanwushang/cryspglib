"""Exact affine rank controls; never require archive regeneration."""
import unittest
from fractions import Fraction
from audit_irrep_archives import rank


class AffineRankTests(unittest.TestCase):
    def test_fixed_and_missing_parameters(self):
        self.assertEqual(rank([None, (0, 0, 0), None]), 0)

    def test_dependent_parameters_are_one_line(self):
        self.assertEqual(rank([(1, 2, 3), (2, 4, 6), (-1, -2, -3)]), 1)

    def test_planes_and_volume(self):
        self.assertEqual(rank([(1, 0, 1), (0, 1, 1), (1, 1, 2)]), 2)
        self.assertEqual(rank([(1, 0, 1), (0, 1, 1), (1, 1, 3)]), 3)

    def test_small_nonzero_direction_is_not_rounded_away(self):
        self.assertEqual(rank([(1, 0, 0), (1, Fraction(1, 10**80), 0)]), 2)


if __name__ == "__main__":
    unittest.main()
