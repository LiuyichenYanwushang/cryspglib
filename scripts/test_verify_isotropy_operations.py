"""Offline regressions for the isotropy *operation* gate.

Run with: python3 -m unittest discover -s scripts -p test_verify_isotropy_operations.py

The oracle process and the archive readers are mocked; the real parsers, the
label legend matching, the group checks, the rotation-action convention and the
gate exit status are exercised.  The checked rows are the ones pinned in
tests/data/isotropy/operations.json (iso 9.6.1, SET I ALL OR 1).
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


# Row-action matrices (`x' = x M`) exactly as stored in `data_space.txt`
# `ipoint_op`, keyed by the label the program prints.  Kept local so this test
# runs without extracting iso.zip; the anchor labels are the ones
# `verify_isotropy_operations.LEGEND_ANCHORS` checks against the Stokes column.
ROW_ACTION_LEGEND = {
    "E": (1, 0, 0, 0, 1, 0, 0, 0, 1),  # 1
    "I": (-1, 0, 0, 0, -1, 0, 0, 0, -1),  # -1
    "C2x": (1, 0, 0, 0, -1, 0, 0, 0, -1),  # 2[100]
    "C2y": (-1, 0, 0, 0, 1, 0, 0, 0, -1),  # 2[010]
    "C2z": (-1, 0, 0, 0, -1, 0, 0, 0, 1),  # 2[001]
    "C2a": (0, 1, 0, 1, 0, 0, 0, 0, -1),  # 2[110]
    "C2b": (0, -1, 0, -1, 0, 0, 0, 0, -1),  # 2[-110]
    "C21''": (1, 0, 0, -1, -1, 0, 0, 0, -1),  # 2[100] (hexagonal)
    "SGv1": (-1, 0, 0, 1, 1, 0, 0, 0, 1),  # -2[100] (hexagonal)
    "C4x+": (1, 0, 0, 0, 0, 1, 0, -1, 0),  # 4[100]
    "C4x-": (1, 0, 0, 0, 0, -1, 0, 1, 0),  # 4[-100]
    "C4z+": (0, 1, 0, -1, 0, 0, 0, 0, 1),  # 4[001]
    "C4z-": (0, -1, 0, 1, 0, 0, 0, 0, 1),  # 4[00-1]
    "S4x-": (-1, 0, 0, 0, 0, -1, 0, 1, 0),  # -4[100]
    "S4x+": (-1, 0, 0, 0, 0, 1, 0, -1, 0),  # -4[-100]
    "SGx": (-1, 0, 0, 0, 1, 0, 0, 0, 1),  # -2[100]
    "SGy": (1, 0, 0, 0, -1, 0, 0, 0, 1),  # -2[010]
    "SGz": (1, 0, 0, 0, 1, 0, 0, 0, -1),  # -2[001]
    "SGda": (0, -1, 0, -1, 0, 0, 0, 0, 1),  # -2[110]
    "SGdb": (0, 1, 0, 1, 0, 0, 0, 0, 1),  # -2[-110]
}
LEGEND = {
    label: {"rotation": rotation, "stokes": f"op{index}", "index": index}
    for index, (label, rotation) in enumerate(ROW_ACTION_LEGEND.items(), start=1)
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
# SG 167 GM3+ P1 -> #15 C2/c: the record that discriminates row from column
# action, since #15's rotation set is not closed under transposition.
ANOMALOUS_ROW = (
    "**Subgroup Size Dir Basis Vectors                          Origin         Elements\n"
    "15 C2/c  1    P1  (1/3,2/3,2/3),(-1,0,0),(-1/3,-2/3,1/3) (-1/6,1/6,1/6) "
    "(E|0,0,0), (C21''|0,0,1/2), (I|-1/3,1/3,1/3), (SGv1|-2/3,-1/3,1/6)\n*\n"
)
CASE = {"sg": 221, "ml": "GM4+", "direction": "P1", "subgroup": 83, "size": 1, "elements": 8}
C167 = {"sg": 167, "ml": "GM3+", "direction": "P1", "subgroup": 15, "size": 1, "elements": 4}
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
MACHINE_167 = {
    "subgroup": 15,
    "size": 1,
    # Stored values for this record (data_isotropy.txt): the basis is in the
    # parent *primitive* frame and the origin in the record setting, so the
    # conversion under test is a real one (unlike the P-parent 221 row).
    "basis": [[Fraction(0), Fraction(0), Fraction(1)],
              [Fraction(0), Fraction(1), Fraction(0)],
              [Fraction(-1), Fraction(0), Fraction(0)]],
    "origin": [Fraction(0), Fraction(1, 2), Fraction(0)],
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


def anchor_legend(overrides=None):
    """A complete 72-entry synthetic legend containing every anchor label.

    The anchors come first (so `load_point_op_legend` can check the action
    convention), the rest are unique identity entries.
    """
    labels = list(ROW_ACTION_LEGEND)
    rotations = list(ROW_ACTION_LEGEND.values())
    overrides = overrides or {}
    for label, rotation in overrides.items():
        rotations[labels.index(label)] = rotation
    while len(labels) < 72:
        index = len(labels)
        labels.append(f"L{index}")
        rotations.append(ROW_ACTION_LEGEND["E"])
        assert len(labels) == len(rotations)
    return labels, [f"op{index}" for index in range(len(labels))], rotations


class LegendTests(unittest.TestCase):
    def write(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="latin-1")
        handle.write(text)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_valid_legend_is_read_with_its_axis_notation(self):
        path = self.write(synthetic_data_space(*anchor_legend()))
        legend = ops.load_point_op_legend(path)
        self.assertEqual(len(legend), 72)
        self.assertEqual(legend["I"]["rotation"], ROW_ACTION_LEGEND["I"])
        self.assertEqual(
            legend["C4z+"]["stokes"], f"op{list(ROW_ACTION_LEGEND).index('C4z+')}"
        )

    def test_legend_anchors_reject_a_transposed_convention(self):
        # Column-action form of `4[001]`: it maps [100] to [0,-1,0] under the
        # row action, so the stored legend would contradict its Stokes column.
        transposed = (0, -1, 0, 1, 0, 0, 0, 0, 1)
        path = self.write(synthetic_data_space(*anchor_legend({"C4z+": transposed})))
        with self.assertRaisesRegex(RuntimeError, "assumed convention"):
            ops.load_point_op_legend(path)

    def test_duplicate_label_with_a_different_rotation_is_rejected(self):
        labels, stokes, rotations = anchor_legend()
        labels[-1] = "E"
        rotations[-1] = ROW_ACTION_LEGEND["C2x"]
        path = self.write(synthetic_data_space(labels, stokes, rotations))
        with self.assertRaisesRegex(RuntimeError, "two different rotations"):
            ops.load_point_op_legend(path)

    def test_non_unimodular_rotation_is_rejected(self):
        labels, stokes, rotations = anchor_legend()
        rotations[-1] = (2, 0, 0, 0, 1, 0, 0, 0, 1)
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
        path = self.write(synthetic_data_space(["E"], ["1"], [ROW_ACTION_LEGEND["E"]], orders, per_sg))
        table = ops.load_point_group_orders(path)
        self.assertEqual(table[1], 1)
        self.assertEqual(table[221], orders[per_sg[220] - 1])
        per_sg[220] = 33
        path = self.write(
            synthetic_data_space(["E"], ["1"], [ROW_ACTION_LEGEND["E"]], orders, per_sg)
        )
        with self.assertRaisesRegex(RuntimeError, "out-of-range point group"):
            ops.load_point_group_orders(path)


class ElementTests(unittest.TestCase):
    def test_exact_translations_are_kept_verbatim(self):
        legend = dict(LEGEND)
        elements = ops.parse_elements(
            "(E|0,0,0), (C2x|0,2,2), (I|5/2,5/2,5/2), (SGv1|-2/3,-1/3,1/6)", legend
        )
        self.assertEqual([element["translation"] for element in elements],
                         [["0", "0", "0"], ["0", "2", "2"], ["5/2", "5/2", "5/2"],
                          ["-2/3", "-1/3", "1/6"]])
        self.assertEqual(elements[3]["raw"], "(SGv1|-2/3,-1/3,1/6)")

    def test_decoded_rotations_are_column_action(self):
        # `C4z+` is `4[001]`: stored row action maps [100] -> [0,1,0], so the
        # decoded column-action matrix must map [100] -> [0,-1,0] instead.
        decoded = ops.parse_elements("(C4z+|0,0,0)", LEGEND)[0]
        self.assertEqual(decoded["rotation"], [0, -1, 0, 1, 0, 0, 0, 0, 1])
        column_action = [
            [Fraction(decoded["rotation"][3 * i + j]) for j in range(3)] for i in range(3)
        ]
        # Both matrices describe the same operation: [100] -> [010].
        # Row action: e_x . M is the first row of M.
        row_matrix = [
            [Fraction(ROW_ACTION_LEGEND["C4z+"][3 * i + j]) for j in range(3)]
            for i in range(3)
        ]
        self.assertEqual(row_matrix[0], [Fraction(0), Fraction(1), Fraction(0)])
        # Column action: R . e_x is the first column of the decoded matrix.
        self.assertEqual(
            [column_action[i][0] for i in range(3)],
            [Fraction(0), Fraction(1), Fraction(0)],
        )
        legend_matrix = [
            [Fraction(ROW_ACTION_LEGEND["C4z+"][3 * i + j]) for j in range(3)] for i in range(3)
        ]
        self.assertEqual(
            column_action,
            [[legend_matrix[j][i] for j in range(3)] for i in range(3)],
        )

    def test_unknown_label_is_not_matched_fuzzily(self):
        # `C4y+`, a lower-case spelling and a bare `C4x` are all absent from the
        # legend, so they must never be matched to `C4x+`/`C4z+`.
        for text in ("(C4y+|0,0,0)", "(c4x+|0,0,0)", "(C4x|0,0,0)"):
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


class LatticeConventionTests(unittest.TestCase):
    """The 167 GM3+ P1 row (#15 C2/c) is what pins the rotation convention."""

    CENTRING_C = [
        [Fraction(1, 2), Fraction(1, 2), Fraction(0)],
        [Fraction(-1, 2), Fraction(1, 2), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1)],
    ]
    BASIS_167 = [
        [Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)],
        [Fraction(-1), Fraction(0), Fraction(0)],
        [Fraction(-1, 3), Fraction(-2, 3), Fraction(1, 3)],
    ]

    def setUp(self):
        self.lattice = ops.child_lattice(self.BASIS_167, self.CENTRING_C)

    def test_child_lattice_is_the_parent_r_lattice(self):
        # Same lattice as the parent's primitive basis: 167 is R-centred and
        # this subgroup has Size = 1.
        parent = [
            [Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)],
            [Fraction(-1, 3), Fraction(1, 3), Fraction(1, 3)],
            [Fraction(-1, 3), Fraction(-2, 3), Fraction(1, 3)],
        ]
        self.assertEqual(ops.determinant3([v for row in self.lattice for v in row]),
                         ops.determinant3([v for row in parent for v in row]))
        self.assertTrue(ops.geometry_gate.same_lattice(self.lattice, parent))

    def test_only_the_transposed_rotation_preserves_the_lattice(self):
        row_action = [
            [Fraction(value) for value in (1, 0, 0)],
            [Fraction(value) for value in (-1, -1, 0)],
            [Fraction(value) for value in (0, 0, -1)],
        ]
        column_action = [[row_action[j][i] for j in range(3)] for i in range(3)]
        self.assertFalse(ops.preserves_lattice(self.lattice, row_action))
        self.assertTrue(ops.preserves_lattice(self.lattice, column_action))
        # The other three elements are insensitive to the convention here.
        for rotation in (
            [[Fraction(1), Fraction(0), Fraction(0)], [Fraction(0), Fraction(1), Fraction(0)],
             [Fraction(0), Fraction(0), Fraction(1)]],
            [[Fraction(-1), Fraction(0), Fraction(0)], [Fraction(0), Fraction(-1), Fraction(0)],
             [Fraction(0), Fraction(0), Fraction(-1)]],
        ):
            self.assertTrue(ops.preserves_lattice(self.lattice, rotation))

    def test_the_167_row_passes_the_full_check(self):
        row = ops.parse_row(ANOMALOUS_ROW)
        checked = ops.check_case(
            C167, row, LEGEND, {15: 4}, MACHINE_167
        )
        self.assertEqual(len(checked["elements"]), 4)
        labels = [element["label"] for element in checked["elements"]]
        self.assertEqual(labels, ["E", "C21''", "I", "SGv1"])

    def test_an_off_grid_subgroup_translation_is_rejected(self):
        row = ops.parse_row(ANOMALOUS_ROW)
        elements = ops.parse_elements(row["elements_text"], LEGEND)
        elements[-1] = dict(elements[-1], translation=["1/5", "0", "0"])

        def fake_parse(text, legend):
            return [dict(element) for element in elements]

        with mock.patch.object(ops, "parse_elements", side_effect=fake_parse):
            with self.assertRaisesRegex(RuntimeError, "source grid"):
                ops.check_case(C167, row, LEGEND, {15: 4}, MACHINE_167)

    def test_a_row_action_rotation_set_is_rejected_by_check_case(self):
        row = ops.parse_row(ANOMALOUS_ROW)
        elements = ops.parse_elements(row["elements_text"], LEGEND)
        for element in elements:
            matrix = [
                [element["rotation"][3 * i + j] for j in range(3)] for i in range(3)
            ]
            element["rotation"] = [value for row_ in
                                   [[matrix[j][i] for j in range(3)] for i in range(3)]
                                   for value in row_]

        def fake_parse(text, legend):
            return [dict(element) for element in elements]

        with mock.patch.object(ops, "parse_elements", side_effect=fake_parse):
            with self.assertRaisesRegex(RuntimeError, "does not preserve the subgroup lattice"):
                ops.check_case(C167, row, LEGEND, {15: 4}, MACHINE_167)


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
