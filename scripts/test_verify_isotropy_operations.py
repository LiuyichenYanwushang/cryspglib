"""Offline regressions for the isotropy *operation* gate.

Run with: python3 -m unittest discover -s scripts -p test_verify_isotropy_operations.py

The oracle process and the archive readers are mocked; the real parsers, the
label legend matching, the group checks and the gate exit status are exercised.
The checked rows are the ones pinned in tests/data/isotropy/operations.json
(iso 9.6.1, SET I ALL OR 1).
"""

import contextlib
import io
import json
import os
import tempfile
import unittest
from fractions import Fraction
from types import SimpleNamespace
from unittest import mock

import verify_isotropy_operations as ops


# The cubic label block the printed rows below use, with the archive's own
# rotations (data_space.txt `point_op_label`/`ipoint_op`).  Kept local so this
# test runs without extracting iso.zip.
CUBIC = {
    "E": (1, 0, 0, 0, 1, 0, 0, 0, 1),
    "C2x": (1, 0, 0, 0, -1, 0, 0, 0, -1),
    "C2y": (-1, 0, 0, 0, 1, 0, 0, 0, -1),
    "C2z": (-1, 0, 0, 0, -1, 0, 0, 0, 1),
    "I": (-1, 0, 0, 0, -1, 0, 0, 0, -1),
    "SGx": (-1, 0, 0, 0, 1, 0, 0, 0, 1),
    "SGy": (1, 0, 0, 0, -1, 0, 0, 0, 1),
    "SGz": (1, 0, 0, 0, 1, 0, 0, 0, -1),
    "C4x+": (1, 0, 0, 0, 0, -1, 0, 1, 0),
    "C4x-": (1, 0, 0, 0, 0, 1, 0, -1, 0),
    "S4x-": (-1, 0, 0, 0, 0, 1, 0, -1, 0),
    "S4x+": (-1, 0, 0, 0, 0, -1, 0, 1, 0),
    # Diagonal axis names, as printed for #123/#12 rows.
    "C2a": (0, 1, 0, 1, 0, 0, 0, 0, -1),
    "SGda": (0, -1, 0, -1, 0, 0, 0, 0, 1),
}
LEGEND = {
    label: {"rotation": rotation, "stokes": f"op{index}", "index": index}
    for index, (label, rotation) in enumerate(CUBIC.items(), start=1)
}

# SG 221 GM4+ P1, verbatim from the oracle.
ROW = (
    "**Subgroup Size Dir Basis Vectors            Origin  Elements\n"
    "83 P4/m  1    P1  (0,0,1),(0,-1,0),(1,0,0) (0,0,0) "
    "(E|0,0,0), (C2x|0,0,0), (C4x+|0,0,0), (C4x-|0,0,0), (I|0,0,0), "
    "(SGx|0,0,0), (S4x-|0,0,0), (S4x+|0,0,0)\n*\n"
)
# SG 221 GM3+ P1: the Elements column wraps onto a second physical line.
WRAPPED_ROW = (
    "**Subgroup   Size Dir Basis Vectors           Origin  Elements\n"
    "123 P4/mmm 1    P1  (1,0,0),(0,1,0),(0,0,1) (0,0,0) "
    "(E|0,0,0), (C2x|0,0,0), (C2y|0,0,0), (C2z|0,0,0), (I|0,0,0), "
    "(SGx|0,0,0), (SGy|0,0,0), (SGz|0,0,0)\n"
    "                    (C2a|0,0,0), (SGda|0,0,0)\n"
    "*\n"
)
CASE = {"sg": 221, "ml": "GM4+", "direction": "P1", "subgroup": 83, "size": 1, "elements": 8}
MACHINE_ROW = {
    "subgroup": 83,
    "size": 1,
    # Parent 221 is P, so the stored primitive basis is the conventional one.
    # Row-major (0,0,1),(0,-1,0),(1,0,0), i.e. the printed basis itself.
    "basis": [[Fraction(0), Fraction(0), Fraction(1)],
              [Fraction(0), Fraction(-1), Fraction(0)],
              [Fraction(1), Fraction(0), Fraction(0)]],
    "origin": [Fraction(0)] * 3,
}


def synthetic_data_space(labels, stokes, rotations, orders=None, per_sg=None):
    """Build a minimal data_space.txt for the two archive readers."""
    lines = ["point_op_label"]
    for start in range(0, len(labels), 8):
        lines.append(" ".join(f'"{label} "' for label in labels[start:start + 8]))
    lines.append("point_op_label_stokes")
    for start in range(0, len(stokes), 8):
        lines.append(" ".join(f'"{token} "' for token in stokes[start:start + 8]))
    lines.append("ipoint_op")
    numbers = [str(value) for rotation in rotations for value in rotation]
    for start in range(0, len(numbers), 12):
        lines.append(" ".join(numbers[start:start + 12]))
    lines.append("ispace_point_group")
    lines.append(" ".join(str(value) for value in (per_sg or [1] * 230)))
    lines.append("ipoint_group_order")
    lines.append(" ".join(str(value) for value in (orders or [1] * 32)))
    lines.append("end")
    return "\n".join(lines) + "\n"


class LegendTests(unittest.TestCase):
    def write(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="latin-1")
        handle.write(text)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_valid_legend_is_read_with_its_axis_notation(self):
        labels = ["E", "I"] + [f"L{index}" for index in range(3, 73)]
        stokes = ["1", "-1"] + [f"op{index}" for index in range(3, 73)]
        rotations = [CUBIC["E"], CUBIC["I"]] + [CUBIC["E"]] * 70
        path = self.write(synthetic_data_space(labels, stokes, rotations))
        legend = ops.load_point_op_legend(path)
        self.assertEqual(len(legend), 72)
        self.assertEqual(legend["I"]["rotation"], CUBIC["I"])
        self.assertEqual(legend["I"]["stokes"], "-1")

    def test_duplicate_label_with_a_different_rotation_is_rejected(self):
        labels = ["E", "I"] + [f"L{index}" for index in range(3, 72)] + ["E"]
        stokes = ["1"] * 72
        rotations = [CUBIC["E"], CUBIC["I"]] + [CUBIC["E"]] * 69 + [CUBIC["C2x"]]
        path = self.write(synthetic_data_space(labels, stokes, rotations))
        with self.assertRaisesRegex(RuntimeError, "two different rotations"):
            ops.load_point_op_legend(path)

    def test_non_unimodular_rotation_is_rejected(self):
        labels = ["E", "I"] + [f"L{index}" for index in range(3, 73)]
        stokes = ["1"] * 72
        rotations = [CUBIC["E"], CUBIC["I"]] + [CUBIC["E"]] * 69 + [(2, 0, 0, 0, 1, 0, 0, 0, 1)]
        path = self.write(synthetic_data_space(labels, stokes, rotations))
        with self.assertRaisesRegex(RuntimeError, "not unimodular"):
            ops.load_point_op_legend(path)

    def test_missing_section_is_rejected(self):
        path = self.write("point_op_label\n\"E \"\n")
        with self.assertRaisesRegex(RuntimeError, "missing section"):
            ops.load_point_op_legend(path)

    def test_point_group_orders_are_read_per_space_group(self):
        orders = list(range(1, 33))
        per_sg = [(index % 32) + 1 for index in range(230)]
        path = self.write(synthetic_data_space(["E"], ["1"], [CUBIC["E"]], orders, per_sg))
        table = ops.load_point_group_orders(path)
        self.assertEqual(table[1], 1)
        self.assertEqual(table[221], orders[per_sg[220] - 1])
        per_sg[220] = 33
        path = self.write(synthetic_data_space(["E"], ["1"], [CUBIC["E"]], orders, per_sg))
        with self.assertRaisesRegex(RuntimeError, "out-of-range point group"):
            ops.load_point_group_orders(path)


class ElementTests(unittest.TestCase):
    def test_exact_translations_are_kept_verbatim(self):
        legend = dict(LEGEND, SGv1={"rotation": CUBIC["SGx"], "stokes": "-2[100]", "index": 73})
        elements = ops.parse_elements(
            "(E|0,0,0), (C2x|0,2,2), (I|5/2,5/2,5/2), (SGv1|-2/3,-1/3,1/6)", legend
        )
        self.assertEqual([element["translation"] for element in elements],
                         [["0", "0", "0"], ["0", "2", "2"], ["5/2", "5/2", "5/2"],
                          ["-2/3", "-1/3", "1/6"]])
        self.assertEqual(elements[3]["raw"], "(SGv1|-2/3,-1/3,1/6)")

    def test_unknown_label_is_not_matched_fuzzily(self):
        for text in ("(C4z+|0,0,0)", "(c4x+|0,0,0)", "(C4x|0,0,0)"):
            with self.subTest(text=text), \
                    self.assertRaisesRegex(RuntimeError, "unknown operation label"):
                ops.parse_elements(text, LEGEND)

    def test_duplicate_and_malformed_elements_fail(self):
        with self.assertRaisesRegex(RuntimeError, "duplicate element"):
            ops.parse_elements("(C2x|0,0,0), (C2x|0,0,0)", LEGEND)
        for text in ("(C2x|0,0)", "(C2x|0,0,0,0)", "(C2x|a,0,0)", "(C2x 0,0,0)"):
            with self.subTest(text=text), \
                    self.assertRaisesRegex(RuntimeError, "malformed|unparsable"):
                ops.parse_elements(text, LEGEND)

    def test_truncated_element_list_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "unparsable"):
            ops.parse_elements("(E|0,0,0), (C2x|0,0", LEGEND)
        with self.assertRaisesRegex(RuntimeError, "empty"):
            ops.parse_elements("", LEGEND)


class RowTests(unittest.TestCase):
    def test_single_row_is_parsed_exactly(self):
        row = ops.parse_row(ROW)
        self.assertEqual(row["subgroup"], 83)
        self.assertEqual(row["size"], 1)
        self.assertEqual(row["direction"], "P1")
        self.assertEqual(row["basis"], [["0", "0", "1"], ["0", "-1", "0"], ["1", "0", "0"]])
        self.assertEqual(row["origin"], ["0", "0", "0"])
        self.assertEqual(len(ops.parse_elements(row["elements_text"], LEGEND)), 8)

    def test_wrapped_elements_are_joined(self):
        row = ops.parse_row(WRAPPED_ROW)
        self.assertEqual(row["subgroup"], 123)
        labels = [element["label"] for element in ops.parse_elements(row["elements_text"], LEGEND)]
        self.assertEqual(labels[-2:], ["C2a", "SGda"])
        self.assertEqual(len(labels), 10)

    def test_two_rows_mean_the_direction_filter_was_lost(self):
        with self.assertRaisesRegex(RuntimeError, "more than one direction row"):
            ops.parse_row(
                ROW.replace(
                    "*\n",
                    "83 P4/m  1    P2  (1,-1,0),(1,1,0),(0,0,1) (0,0,0) (E|0,0,0)\n*\n",
                )
            )

    def test_syntax_error_paging_and_missing_header_fail(self):
        with self.assertRaisesRegex(RuntimeError, "syntax error"):
            ops.parse_row("***Syntax error: VALUE\n*\n")
        with self.assertRaisesRegex(RuntimeError, "paged"):
            ops.parse_row("**Subgroup Size Dir Basis Vectors Origin Elements\n--More--\n")
        with self.assertRaisesRegex(RuntimeError, "table header"):
            ops.parse_row("**Subgroup Size Dir Basis Vectors Origin\n83 P4/m 1 P1 ...\n")

    def test_missing_elements_column_is_not_silently_accepted(self):
        # `SHOW ELEMENTS` without `VALUE DIRECTION` prints an empty column.
        with self.assertRaisesRegex(RuntimeError, "unparsable|malformed|empty"):
            ops.parse_row(ROW.replace(
                "(E|0,0,0), (C2x|0,0,0), (C4x+|0,0,0), (C4x-|0,0,0), (I|0,0,0), "
                "(SGx|0,0,0), (S4x-|0,0,0), (S4x+|0,0,0)", ""))


class CaseTests(unittest.TestCase):
    def check(self, **overrides):
        case = dict(CASE, **{key: overrides[key] for key in overrides if key in CASE})
        row = ops.parse_row(ROW)
        return ops.check_case(case, row, LEGEND, {83: 8}, overrides.get("machine", MACHINE_ROW))

    def test_pinned_row_passes(self):
        row = self.check()
        self.assertEqual(len(row["elements"]), 8)

    def test_wrong_count_point_group_and_subgroup_fail(self):
        for overrides, message in (
            ({"elements": 4}, "printed 8 elements"),
            ({"subgroup": 12}, "subgroup 83"),
            ({"machine": dict(MACHINE_ROW, subgroup=12)}, "machine record"),
            ({"machine": dict(MACHINE_ROW, size=2)}, "machine record"),
            ({"machine": None}, "no such record"),
        ):
            with self.subTest(overrides=overrides), self.assertRaisesRegex(RuntimeError, message):
                self.check(**overrides)

    def test_a_pinned_count_that_disagrees_with_the_point_group_fails(self):
        with self.assertRaisesRegex(RuntimeError, "point group order"):
            ops.check_case(
                dict(CASE, elements=8), ops.parse_row(ROW), LEGEND, {83: 16}, MACHINE_ROW
            )

    def test_rotations_must_form_a_group(self):
        # Same count, but one element is swapped for a rotation of a different
        # subgroup: the eight rotations are then not closed.
        row = ops.parse_row(ROW.replace("(S4x+|0,0,0)", "(C2y|0,0,0)"))
        with self.assertRaisesRegex(RuntimeError, "closed"):
            ops.check_case(CASE, row, LEGEND, {83: 8}, MACHINE_ROW)

    def test_a_non_zero_identity_translation_fails(self):
        row = ops.parse_row(ROW.replace("(E|0,0,0)", "(E|1,0,0)"))
        with self.assertRaisesRegex(RuntimeError, "identity carries translation"):
            ops.check_case(CASE, row, LEGEND, {83: 8}, MACHINE_ROW)

    def test_origin_and_lattice_disagreements_fail(self):
        with self.assertRaisesRegex(RuntimeError, "origin"):
            self.check(machine=dict(MACHINE_ROW, origin=[Fraction(1, 2), 0, 0]))
        # A unimodular re-orientation spans the same lattice (the geometry gate
        # owns orientation); a non-unimodular one must fail here.
        doubled = [list(row) for row in MACHINE_ROW["basis"]]
        doubled[0][2] = Fraction(2)
        with self.assertRaisesRegex(RuntimeError, "different lattices|volume"):
            self.check(machine=dict(MACHINE_ROW, basis=doubled))


class MainTests(unittest.TestCase):
    def run_main(self, stored, fresh, write=False):
        handle = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(stored, handle)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        argv = ["--write"] if write else []
        with mock.patch.object(ops, "FIXTURE", handle.name), \
                mock.patch.object(ops, "build_fixture", return_value=(fresh, 1)), \
                contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            status = ops.main(argv)
        with open(handle.name, encoding="utf-8") as reader:
            return status, json.load(reader)

    def test_replay_of_an_identical_fixture_succeeds(self):
        fresh = {"provenance": {}, "cases": [dict(CASE)]}
        status, stored = self.run_main(fresh, fresh)
        self.assertEqual(status, 0)
        self.assertEqual(stored, fresh)

    def test_a_changed_oracle_result_fails_the_replay(self):
        stored = {"provenance": {}, "cases": [dict(CASE)]}
        changed = {"provenance": {}, "cases": [dict(CASE, subgroup=12)]}
        status, _ = self.run_main(stored, changed)
        self.assertEqual(status, 1)

    def test_write_mode_stores_the_fresh_fixture(self):
        changed = {"provenance": {}, "cases": [dict(CASE, subgroup=12)]}
        status, stored = self.run_main({"provenance": {}, "cases": []}, changed, write=True)
        self.assertEqual(status, 0)
        self.assertEqual(stored, changed)

    def test_process_failure_is_not_an_empty_success(self):
        response = SimpleNamespace(returncode=1, stdout="", stderr="failed")
        with mock.patch.object(ops.os.path, "isfile", return_value=True), \
                mock.patch.object(ops.subprocess, "run", return_value=response), \
                self.assertRaisesRegex(RuntimeError, "iso exited with 1"):
            ops.run_oracle(221, "GM4+", "P1")

    def test_a_wrong_setting_banner_is_rejected(self):
        response = SimpleNamespace(returncode=0, stdout="Isotropy\n", stderr="")
        with mock.patch.object(ops.os.path, "isfile", return_value=True), \
                mock.patch.object(ops.subprocess, "run", return_value=response), \
                self.assertRaisesRegex(RuntimeError, "pinned setting"):
            ops.run_oracle(221, "GM4+", "P1")


if __name__ == "__main__":
    unittest.main()
