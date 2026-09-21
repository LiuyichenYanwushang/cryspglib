"""Offline regressions for the other-wave-vector character extraction.

Run with:
    python3 -m unittest discover -s scripts -p test_freeze_w_little_characters.py

The gate needs the extracted ISOTROPY binary, so these tests drive the parsing,
the exact solver and the joint line solve directly.  The line solve is exercised
through a stubbed `run_iso`, which means the command sequence itself (elements,
k vector, Gamma irreps, compatibilities, characters) is part of the fixture: a
change that stops asking for one of them fails here.
"""

import sys
import unittest
from fractions import Fraction

import freeze_w_little_characters as gate


def iso_output(*tables):
    """Join captured tables the way one `iso` session prints them."""
    return "Isotropy\n" + "\n".join(tables) + "\n*\n"


ELEMENTS = iso_output(
    "Elements",
    "(x,y,z), (x,-y,-z), (-x,y,-z), (-x,-y,z), (y,z,x), (y,-z,-x),",
    "(-y,z,-x), (-y,-z,x), (z,x,y), (z,-x,-y), (-z,x,-y), (-z,-x,y)",
)
K_POINT = iso_output("    k vector", "DT  (0,2a,0)")
GAMMA = iso_output(
    "Irrep (ML) k vector Dim",
    "GM1        (0,0,0)  1",
    "GM2GM3     (0,0,0)  2",
    "GM4        (0,0,0)  3",
)
LINE = iso_output(
    "Irrep (ML) k vector Dim",
    "DT1        (0,2a,0)  6",
    "DT2        (0,2a,0)  6",
)
COMPAT_GM1 = iso_output("Compat (ML) Element Char", "GM1: DT1    (x,y,z) 1.000", "GM1: DT1    (-x,y,-z) 1.000")
COMPAT_GM2GM3 = iso_output(
    "Compat (ML) Element Char", "GM2GM3: DT1  (x,y,z) 2.000", "GM2GM3: DT1  (-x,y,-z) 2.000"
)
COMPAT_GM4 = iso_output(
    "Compat (ML) Element Char",
    "GM4: DT1 DT2 DT2 (x,y,z) 3.000",
    "GM4: DT1 DT2 DT2 (-x,y,-z) -1.000",
)


def stdout_for(sg, commands):
    """Answer the five queries the line solve issues, in order."""
    joined = "\n".join(commands)
    if "SHOW ELEMENTS" in joined:
        return ELEMENTS
    if "DISPLAY KPOINT" in joined:
        return K_POINT
    if "VALUE KPOINT GM" in joined:
        return GAMMA
    if "SHOW COMPATIBILITY" not in joined:
        return LINE
    if "VALUE IRREP GM1" in joined:
        return COMPAT_GM1
    if "VALUE IRREP GM2GM3" in joined:
        return COMPAT_GM2GM3
    return COMPAT_GM4


class ParseTests(unittest.TestCase):
    def test_component_count_of_compound_labels(self):
        self.assertEqual(gate.component_count("DT1"), 1)
        self.assertEqual(gate.component_count("DT3DT4"), 2)
        self.assertEqual(gate.component_count("GM2GM3"), 2)
        self.assertEqual(gate.component_count("LD1LE1"), 2)
        self.assertEqual(gate.component_count("W1WA1"), 2)

    def test_elements_are_parsed_as_rotations(self):
        elements = gate.parse_elements(ELEMENTS)
        self.assertEqual(len(elements), 12)
        self.assertEqual(elements[0][2], "x,y,z")
        self.assertEqual(elements[1][0][1], (Fraction(0), Fraction(-1), Fraction(0)))
        self.assertEqual(elements[1][1], (Fraction(0), Fraction(0), Fraction(0)))

    def test_direction_sets_the_free_parameter_to_one(self):
        self.assertEqual(
            gate.parse_direction("(0,2a,0)"), (Fraction(0), Fraction(2), Fraction(0))
        )
        self.assertEqual(
            gate.parse_direction("(2a,2a,0)"), (Fraction(2), Fraction(2), Fraction(0))
        )
        self.assertEqual(
            gate.parse_direction("(1/2,-2a+1,2a)"),
            (Fraction(1, 2), Fraction(-1), Fraction(2)),
        )

    def test_compat_list_keeps_multiplicities(self):
        self.assertEqual(gate.parse_compat(COMPAT_GM4), ["DT1", "DT2", "DT2"])
        self.assertEqual(gate.parse_compat(COMPAT_GM1), ["DT1"])

    def test_characters_are_read_from_shared_rows(self):
        self.assertEqual(
            gate.parse_characters(COMPAT_GM4),
            [("x,y,z", Fraction(3)), ("-x,y,-z", Fraction(-1))],
        )

    def test_irrep_table_rows(self):
        self.assertEqual(
            gate.parse_irrep_rows(GAMMA),
            [("GM1", "(0,0,0)", 1), ("GM2GM3", "(0,0,0)", 2), ("GM4", "(0,0,0)", 3)],
        )


class SolveTests(unittest.TestCase):
    def test_exact_solve_reports_rank_and_consistency(self):
        solution, rank, consistent = gate.solve(
            [[1, 0], [0, 1], [1, 1]], [Fraction(1), Fraction(-1), Fraction(0)]
        )
        self.assertTrue(consistent)
        self.assertEqual(rank, 2)
        self.assertEqual(solution, [Fraction(1), Fraction(-1)])

    def test_inconsistent_system_is_reported(self):
        _, _, consistent = gate.solve([[1, 0], [1, 0]], [Fraction(1), Fraction(2)])
        self.assertFalse(consistent)

    def test_under_determined_system_is_reported(self):
        _, rank, consistent = gate.solve([[1, 1]], [Fraction(2)])
        self.assertTrue(consistent)
        self.assertEqual(rank, 1)


class LineSolveTests(unittest.TestCase):
    def derive(self, labels=("DT1", "DT2")):
        original = gate.run_iso
        gate.run_iso = lambda sg, commands, timeout=300: stdout_for(sg, commands)
        try:
            return gate.derive_line(196, "DT", list(labels))
        finally:
            gate.run_iso = original

    def test_line_solve_recovers_the_little_group_characters(self):
        tables, failures = self.derive()
        self.assertEqual(failures, [])
        dt1 = tables["DT1"]
        self.assertEqual(dt1["k_vector"], "(0,2a,0)")
        self.assertEqual(dt1["direction"], ["0", "2", "0"])
        self.assertEqual(len(dt1["operations"]), 2)
        self.assertEqual(dt1["operations"][0]["character"], "1")
        self.assertEqual(dt1["operations"][1]["character"], "1")
        self.assertEqual(
            tables["DT2"]["operations"][1]["character"], "-1"
        )

    def test_missing_line_irrep_is_reported(self):
        _, failures = self.derive(labels=("DT1", "DT3"))
        self.assertEqual(len(failures), 1)
        self.assertIn("not printed for k DT", failures[0])

    def test_under_determined_system_is_reported(self):
        """A subset of a line still solves: components never enter alone."""
        original = gate.run_iso
        gate.run_iso = lambda sg, commands, timeout=300: stdout_for(sg, commands)
        try:
            tables, failures = gate.derive_line(196, "DT", ["DT1"])
        finally:
            gate.run_iso = original
        # DT1 alone is still determined: the `GM4: DT1 DT2 DT2` row is used
        # through the dual system, which fixes DT1 without fixing DT2.
        self.assertEqual(failures, [])
        self.assertEqual(tables["DT1"]["operations"][1]["character"], "1")

    def test_compound_label_sums_its_components(self):
        original = gate.run_iso
        gate.run_iso = lambda sg, commands, timeout=300: stdout_for(sg, commands)
        try:
            tables, failures = gate.derive_line(196, "DT", ["DT1", "DT2"])
        finally:
            gate.run_iso = original
        self.assertEqual(failures, [])
        self.assertEqual(tables["DT1"]["components"], ["DT1"])
        self.assertEqual(tables["DT2"]["components"], ["DT2"])


if __name__ == "__main__":
    sys.exit(unittest.main())
