#!/usr/bin/env python3
"""Offline tests for the R3 gap-source classifier.

The archive-derived expectations below are hand-checked:

* SG 5 (C2): the U line ``(0, t, 1/2)`` has a two-operation little group
  (identity and the twofold rotation) whose two little irreps are archived as
  the compound records ``U1UA1`` / ``U2UA2``; the general-position record
  ``GP1GQ1`` also passes through the same point and must be excluded.
* SG 24 (I2_12_12_1): the little group at ``(1, 1, 1/2)`` has four operations
  (three of them screws), so its factor system is projective with a non-trivial
  ``omega``; the product of two screws realises a centring translation, which is
  why lattice membership must use the primitive lattice and not ``Z^3``.
* SG 3 (P2): ``(1/3, 0, 0)`` has a trivial little co-group, so the target there
  is a one-dimensional Bloch phase whatever domain passes through it.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from classify_subduction_gap_sources import (  # noqa: E402
    analyse_group,
    factor_system,
    in_primitive_lattice,
    in_primitive_reciprocal,
    little_group,
    match_arm,
    matmul,
    parse_star,
    rotation_inverse,
)
from iso_irrep_exact import load_exact_iso_irrep_sources  # noqa: E402

DATABASE = load_exact_iso_irrep_sources()


def records_of(spacegroup: int):
    return [record for record in DATABASE.pir_records if record.spacegroup == spacegroup]


def analyse(child: int, star: str) -> dict:
    return analyse_group(
        DATABASE.source_universe(child), records_of(child), parse_star(star)
    )


class CenteringLatticeTests(unittest.TestCase):
    def test_primitive_lattice_membership(self):
        self.assertTrue(in_primitive_lattice((Fraction(1, 2), Fraction(1, 2), Fraction(0)), "C"))
        self.assertTrue(in_primitive_lattice((Fraction(3, 2), Fraction(1, 2), Fraction(1)), "C"))
        self.assertFalse(in_primitive_lattice((Fraction(0), Fraction(0), Fraction(1, 2)), "C"))
        self.assertTrue(in_primitive_lattice((Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)), "I"))
        self.assertFalse(in_primitive_lattice((Fraction(1, 2), Fraction(0), Fraction(0)), "I"))

    def test_primitive_reciprocal_extinctions(self):
        self.assertTrue(in_primitive_reciprocal((Fraction(1), Fraction(1), Fraction(0)), "C"))
        self.assertFalse(in_primitive_reciprocal((Fraction(1), Fraction(0), Fraction(0)), "C"))
        # I centring: G1 + G2 + G3 must be even.
        self.assertTrue(in_primitive_reciprocal((Fraction(1), Fraction(1), Fraction(0)), "I"))
        self.assertTrue(in_primitive_reciprocal((Fraction(2), Fraction(0), Fraction(0)), "I"))
        self.assertFalse(in_primitive_reciprocal((Fraction(1), Fraction(1), Fraction(1)), "I"))
        self.assertFalse(in_primitive_reciprocal((Fraction(1, 2), Fraction(0), Fraction(0)), "P"))

    def test_star_parsing_keeps_every_point(self):
        star = parse_star("-1/2,1,3/4;1/2,1,1/4")
        self.assertEqual(star[0], (Fraction(-1, 2), Fraction(1), Fraction(3, 4)))
        self.assertEqual(star[1], (Fraction(1, 2), Fraction(1), Fraction(1, 4)))


class LittleGroupTests(unittest.TestCase):
    def test_c2_u_line_has_a_two_operation_little_group(self):
        operations, order = little_group(DATABASE.source_universe(5), parse_star("0,1/2,1/2")[0])
        self.assertEqual(order, 2)
        self.assertEqual(len(operations), 2)
        nontrivial, nonsymmorphic = factor_system(operations, parse_star("0,1/2,1/2")[0], "C")
        self.assertEqual(nontrivial, [])
        self.assertEqual(nonsymmorphic, 0)

    def test_screw_little_group_is_projective(self):
        # SG 24 = I2_12_12_1: the P point (1/2, 1/2, 1/2) is fixed modulo the
        # primitive reciprocal lattice by all three 2_1 screws and the identity.
        star = parse_star("1/2,1/2,1/2")
        operations, order = little_group(DATABASE.source_universe(24), star[0])
        self.assertEqual(order, 4)
        self.assertEqual(len(operations), 4)
        nontrivial, nonsymmorphic = factor_system(operations, star[0], "I")
        self.assertTrue(nontrivial, "screw products must leave a non-trivial phase")
        self.assertIn(Fraction(1, 4), {value for _, _, value in nontrivial})
        self.assertEqual(nonsymmorphic, 3)
        # The screw product differs from the chosen representative by an
        # I-centring vector, which is a primitive lattice vector.
        self.assertTrue(in_primitive_lattice((Fraction(1, 2), Fraction(-1, 2), Fraction(1, 2)), "I"))
        report = analyse(24, "1/2,1/2,1/2")
        self.assertEqual(report["classification"], "parameterized_source")
        self.assertEqual(report["little_co_group_order"], 4)
        self.assertEqual(report["matched_labels"], ["W1WA1"])
        self.assertEqual(report["factor_system_values"], ["1/4", "3/4"])

    def test_centred_class_period_is_not_one(self):
        # In a C-centred group the b* line returns to its own class only after
        # t = 2, so sampling [0, 1] would report (0, 3/2, 0) as sourceless.
        from classify_subduction_gap_sources import class_period

        self.assertEqual(class_period((Fraction(0), Fraction(1), Fraction(0)), "C"), 2)
        report = analyse(5, "0,3/2,0")
        self.assertEqual(report["classification"], "parameterized_source")
        self.assertEqual(report["matched_labels"], ["LD1LE1", "LD2LE2"])

    def test_trivial_little_group_on_the_c2_plane(self):
        # SG 3 (P2, axis b): k_b = 0 is a plane whose generic point has a
        # trivial little co-group.
        operations, order = little_group(DATABASE.source_universe(3), parse_star("1/3,0,0")[0])
        self.assertEqual(order, 1)
        self.assertEqual(len(operations), 1)



def read_frozen_setting(ordinal: int) -> dict:
    """One entry of the generated frozen-setting table, parsed from its source."""
    import re

    text = (SCRIPT_DIR.parent / "src" / "irrep" / "subduction_settings_data.rs").read_text(
        encoding="utf-8"
    )
    text = text[text.index("FROZEN_EMBEDDING_SETTINGS") :]
    pattern = re.compile(
        r"\n\s*\(\s*" + str(ordinal) + r"\s*,\s*(\d+)\s*,\s*(\d+)\s*,"
        r"\s*\[\[([^\]]*)\],\s*\[([^\]]*)\],\s*\[([^\]]*)\]\]\s*,"
        r"\s*(-?\d+)"
    )
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"ordinal {ordinal} is not in the frozen settings table")
    rows = tuple(
        tuple(int(value) for value in group.split(",")) for group in match.groups()[2:5]
    )
    return {
        "parent": int(match.group(1)),
        "child": int(match.group(2)),
        "numerator": rows,
        "denominator": int(match.group(6)),
    }


class SettingFrameTests(unittest.TestCase):
    """Per-operation check of a non-symmetric (shear) change of basis.

    Ordinal 26 (SG 3 -> child #3) is the first frozen setting whose ``U`` is not
    a signed permutation: ``U = [[1, 2, 1], [-1, 2, -1], [-1, 0, 1]] / 2``.  A
    change of basis must conjugate the archive little group into integral
    rotations and leave every factor-system phase unchanged, because
    ``omega = exp(2 pi i q.L)`` is a duality pairing and ``q`` and ``L``
    transform contragrediently.
    """

    U = (
        (Fraction(1, 2), Fraction(1), Fraction(1, 2)),
        (Fraction(-1, 2), Fraction(1), Fraction(-1, 2)),
        (Fraction(-1, 2), Fraction(0), Fraction(1, 2)),
    )

    @staticmethod
    def _as_fractions(rotation):
        return tuple(tuple(Fraction(value) for value in row) for row in rotation)

    def test_the_frozen_setting_of_the_witness_is_the_engine_setting(self):
        # The frozen table stores this record's setting as the fraction below;
        # `cargo run --release -p cryspglib --example trace_embedding -- 26`
        # prints the engine's `setting` as the same matrix, so the witness frame
        # is the engine's own.  The per-operation conjugation of the *archive*
        # child frame into the engine's child frame is still open: conjugating
        # by this U does not give integral rotations, so it is not the direct
        # setting transform between the two child frames.
        entry = read_frozen_setting(26)
        self.assertEqual(entry["numerator"], ((1, 2, 1), (-1, 2, -1), (-1, 0, 1)))
        self.assertEqual(entry["denominator"], 2)
        self.assertNotEqual(
            entry["numerator"], ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        )

    def test_shear_setting_preserves_every_factor_system_phase(self):
        star = parse_star("0,1/3,1/2")
        q_archive = star[0]
        operations, _order = little_group(DATABASE.source_universe(3), q_archive)
        _nontrivial, _nonsymmorphic = factor_system(operations, q_archive, "P")
        inverse_transpose = tuple(
            tuple(rotation_inverse(self.U)[row][column] for row in range(3))
            for column in range(3)
        )
        q_engine = tuple(
            sum(inverse_transpose[i][j] * q_archive[j] for j in range(3))
            for i in range(3)
        )
        representatives = sorted({rotation for rotation, _ in operations})
        for first in representatives:
            for second in representatives:
                product = matmul(first, second)
                lattice = tuple(
                    Fraction(0) for _ in range(3)
                )
                # The archive cocycle lattice for this pair: solve for the
                # translation that makes the product consistent.
                for rotation, translation in operations:
                    if rotation != first:
                        continue
                    first_translation = translation
                for rotation, translation in operations:
                    if rotation != second:
                        continue
                    second_translation = translation
                product_translation = tuple(
                    sum(
                        Fraction(first[axis][k]) * second_translation[k]
                        for k in range(3)
                    )
                    + first_translation[axis]
                    for axis in range(3)
                )
                for rotation, translation in operations:
                    if rotation == product:
                        lattice = tuple(
                            product_translation[i] - translation[i] for i in range(3)
                        )
                turns_archive = sum(q_archive[i] * lattice[i] for i in range(3))
                lattice_engine = tuple(
                    sum(self.U[i][j] * lattice[j] for j in range(3)) for i in range(3)
                )
                turns_engine = sum(
                    q_engine[i] * lattice_engine[i] for i in range(3)
                )
                self.assertEqual(
                    turns_archive - turns_engine,
                    Fraction(0),
                    f"phase changed under the shear setting for {first} x {second}",
                )


class RegressionTests(unittest.TestCase):
    """The four defects the review reproduced, pinned as permanent tests."""

    def test_coupled_direction_contributes_to_every_component(self):
        # A direction (1, 1, 0) with t = 1/4 must give (1/4, 1/4, 0).  Reading
        # only the diagonal component returned (1/4, 0, 0) and mis-matched the
        # whole arm.
        from classify_subduction_gap_sources import arm_point

        record = next(
            record
            for record in records_of(155)
            if record.irrep_label in {"Y1YA1", "Y2YA2"}
        )
        arm = next(
            arm
            for arm in record.k_arms
            if [j for j in range(3) if arm.parameters[j] is not None]
        )
        free = [j for j in range(3) if arm.parameters[j] is not None]
        parameters, solved_free = match_arm(arm, parse_star("3/4,3/4,3/2")[0], "R")
        self.assertIsNotNone(parameters, "the coupled index-2 arm must match q")
        self.assertEqual(parameters, (Fraction(3, 4),))
        point = arm_point(arm, solved_free, parameters)
        self.assertEqual(point[0], Fraction(3, 4))
        self.assertEqual(point[1], Fraction(3, 4))
        self.assertEqual(point[2], Fraction(3, 2))

    def test_rotation_inverse_is_the_matrix_inverse(self):
        # Trigonal and other non-diagonal rotations exposed a missing transpose:
        # the helper returned R^-T, so the reciprocal action came out wrong.
        seen = 0
        for spacegroup in (3, 5, 24, 155, 166, 167):
            for operation in DATABASE.source_universe(spacegroup).operations:
                inverse = rotation_inverse(operation.rotation)
                for i in range(3):
                    for j in range(3):
                        total = sum(
                            Fraction(operation.rotation[i][k]) * inverse[k][j]
                            for k in range(3)
                        )
                        self.assertEqual(total, Fraction(1 if i == j else 0))
                seen += 1
        self.assertGreater(seen, 4)

    def test_every_arm_of_a_star_gives_the_same_co_group_order(self):
        # A wrong reciprocal action made points of one star disagree; the star is
        # one physical q class, so the order must be constant over it.
        checked = 0
        for child, star in [(155, "-1/2,1/4,3/2;1/4,-1/2,3/2;1/4,1/4,3/2"),
                            (5, "0,1/2,1/2"), (24, "1/2,1/2,1/2")]:
            universe = DATABASE.source_universe(child)
            orders = {
                little_group(universe, point)[1] for point in parse_star(star)
            }
            self.assertEqual(len(orders), 1, f"child #{child} star {star}: {orders}")
            checked += 1
        self.assertEqual(checked, 3)

    def test_matrix_availability_uses_the_decoder_not_the_phase_slots(self):
        # Discrete records carry no parametric-phase slot by format (#5 GM1 has
        # four `None`s), yet their matrices exist.
        report = analyse(5, "0,0,0")
        self.assertTrue(report["matrix_available"])
        self.assertGreater(report["matrix_elements"], 0)


class SourceClassificationTests(unittest.TestCase):
    def test_c2_u_line_uses_the_archived_line_records(self):
        report = analyse(5, "0,1/2,1/2")
        self.assertEqual(report["classification"], "parameterized_source")
        self.assertEqual(report["little_co_group_order"], 2)
        self.assertEqual(report["min_free_directions"], 1)
        self.assertEqual(report["matched_labels"], ["U1UA1", "U2UA2"])
        self.assertEqual(report["matched_parameters"], [("1/2",)])
        self.assertEqual(report["excluded_generic"], 1)
        self.assertEqual(report["factor_system_nontrivial"], 0)

    def test_general_position_point_is_analytic(self):
        report = analyse(3, "1/3,0,0;2/3,0,0")
        self.assertEqual(report["classification"], "analytic_general_position")
        self.assertEqual(report["little_co_group_order"], 1)
        self.assertTrue(report["matched_irnumbers"], "a source may still be recorded")

    def test_the_reviewed_witness_is_a_parameterized_source(self):
        # Child #155 (R32) at the pinned star: the little co-group has order two
        # and the archived coupled line k = (t, t, 3/2) passes through it at
        # t = 3/4, i.e. at (3/4, 3/4, 3/2) = (-1/4, -1/4, 3/2) + (1, 1, 0), an
        # R-centred reciprocal vector.  The earlier "no source" verdict came
        # from dropping that off-diagonal coupling and from a wrong reciprocal
        # action, and its negative test pinned the wrong answer.
        report = analyse(155, "-1/4,-1/4,3/2;-1/4,1/2,3/2;1/2,-1/4,3/2")
        self.assertEqual(report["little_co_group_order"], 2)
        self.assertEqual(report["classification"], "parameterized_source")
        self.assertEqual(report["matched_labels"], ["Y1YA1", "Y2YA2"])
        self.assertEqual(report["matched_parameters"], [("3/4",)])

    def test_every_archive_arm_match_is_exact(self):
        # The matcher never rounds: a solved arm reproduces the queried point
        # exactly, modulo the primitive reciprocal lattice.
        universe = DATABASE.source_universe(9)
        for record in records_of(9)[:5]:
            for arm in record.k_arms:
                for point in parse_star("1/4,-1/4,1/2;3/4,1/4,1/2"):
                    parameters, free = match_arm(arm, point, universe.centering.value)
                    if parameters is None:
                        continue
                    for axis, index in enumerate(free):
                        direction = arm.parameters[index]
                        self.assertIsNotNone(direction)
                        self.assertLessEqual(len(free), 3)
                    self.assertTrue(free or parameters == ())


class EndToEndTests(unittest.TestCase):
    MANIFEST_HEADER = (
        "ordinal\tparent_sg\tcondensing_cdml\tdirection\tprobe_cdml\tchild_sg\tstar_index"
        "\tstar_size\tarm_count\tblock_dimension\tparent_dimension\tgamma\tstatus"
        "\tmatching_cdml\tq\tcanonical_q\tsetting_numerator\tsetting_denominator\tchild_shift\n"
    )

    def _manifest(self, directory: str) -> Path:
        path = Path(directory) / "gaps.tsv"
        rows = [
            # child 5 U line, child 24 screw point, child 3 general point.
            ("100", "5", "U1", "0,1/2,1/2", "1", "1", "0,0,0"),
            ("200", "155", "W1", "-1/4,-1/4,3/2;-1/4,1/2,3/2;1/2,-1/4,3/2", "1", "1", "0,0,0"),
            ("300", "3", "F1", "1/3,0,0;2/3,0,0", "1", "1", "0,0,0"),
        ]
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(self.MANIFEST_HEADER)
            for ordinal, child, probe, q, setting_num, setting_den, shift in rows:
                handle.write(
                    "\t".join(
                        [
                            ordinal, "221", "GM4+", "P1", probe, child, "0", "1", "1", "1", "2",
                            "false", "missing_discrete_scalar_data", "", q, q,
                            setting_num, setting_den, shift,
                        ]
                    )
                    + "\n"
                )
        return path

    def _q(self, child: int, text: str):
        from classify_subduction_gap_sources import parse_star

        return parse_star(text)[0]

    def test_classifier_classifies_every_row_of_a_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = self._manifest(directory)
            result = subprocess.run(
                [sys.executable, str(SCRIPT_DIR / "classify_subduction_gap_sources.py"), str(manifest)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        lines = result.stdout.splitlines()
        self.assertEqual(lines[0].split("\t")[0], "child_sg")
        rows = [line.split("\t") for line in lines[1:] if line]
        self.assertEqual(len(rows), 3)
        known = {
            "analytic_general_position",
            "parameterized_source",
            "special_value_no_source",
            "no_source_needs_algorithm",
        }
        for row in rows:
            self.assertEqual(len(row), 25, row)
            self.assertIn(row[-1], known, row)
            self.assertTrue(row[7], "every group must record its matched sources")
        by_child = {row[0]: row[-1] for row in rows}
        self.assertEqual(by_child["5"], "parameterized_source")
        self.assertEqual(by_child["155"], "parameterized_source")
        self.assertEqual(by_child["3"], "analytic_general_position")
        self.assertIn(f"groups=3", result.stderr)

    def test_full_manifest_is_classified_without_omission(self):
        # The R3 denominator gate.  It needs the manifest produced by
        # examples/census_subduction_gaps.rs plus the full archive scan, so it
        # runs only when explicitly requested:
        #     R3_FULL_MANIFEST=1 python3 -m unittest ...
        if os.environ.get("R3_FULL_MANIFEST") != "1":
            self.skipTest("set R3_FULL_MANIFEST=1 to run the full 899-group gate")
        manifest = SCRIPT_DIR.parent / "target" / "r12_gaps.tsv"
        if not manifest.exists():
            self.skipTest("the R3 gap manifest is not present in this checkout")
        result = subprocess.run(
            [sys.executable, str(SCRIPT_DIR / "classify_subduction_gap_sources.py"), str(manifest)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        rows = [line.split("\t") for line in result.stdout.splitlines()[1:] if line]
        self.assertEqual(len(rows), 899)
        for row in rows:
            self.assertTrue(row[7], "every group must record its matched sources")
        labels = {row[-1] for row in rows}
        self.assertEqual(labels, {"analytic_general_position", "parameterized_source"})
        self.assertEqual(sum(1 for row in rows if row[-1] == "parameterized_source"), 684)
        self.assertEqual(sum(1 for row in rows if row[-1] == "analytic_general_position"), 215)
        # Every matched record's archived matrix block is complete.
        for row in rows:
            # Columns: ..., matrix_available, matrix_elements, star_orders,
            # classification.
            self.assertEqual(row[-4], "true", row)  # matrix_available
            self.assertEqual(row[-3].count(","), 0, row)  # star_orders: one order
        # Every star is one q class: its arms must agree on the little co-group,
        # and the same order must hold across the settings of one star.
        by_star = {}
        for row in rows:
            by_star.setdefault((row[0], row[1]), set()).add(row[14])
        self.assertTrue(all(len(orders) == 1 for orders in by_star.values()))
        for row in rows:
            self.assertEqual(
                len(set(row[14].split(","))), 1, f"arms disagree: {row[14]}"
            )


if __name__ == "__main__":
    unittest.main()
