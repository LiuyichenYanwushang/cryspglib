#!/usr/bin/env python3
"""Focused tests for the pure exact ISO-IR data-Hall derivation."""

from dataclasses import FrozenInstanceError, fields, is_dataclass
from fractions import Fraction
import hashlib
import io
from pathlib import Path
import re
import unittest
import zipfile
from types import SimpleNamespace

from . import derive_iso_irrep_data_hall as hall
from . import iso_irrep_exact as exact


def _independent_mod1(value):
    return value - (value.numerator // value.denominator)


def _independent_cosets(matrix, denominator):
    generators = tuple(
        tuple(Fraction(matrix[3 * column + row], denominator)
              for row in range(3))
        for column in range(3)
    )
    zero = (Fraction(0), Fraction(0), Fraction(0))
    result = {zero}
    frontier = [zero]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = tuple(
                _independent_mod1(current[index] + generator[index])
                for index in range(3)
            )
            if candidate not in result:
                result.add(candidate)
                frontier.append(candidate)
    return tuple(sorted(result))


def _fortran_centering(member_name, marker_name):
    path = exact.PIR_ARCHIVE_PATH if member_name == "PIR_data.f" else exact.CIR_ARCHIVE_PATH
    expected_size = exact.PIR_ARCHIVE_BYTES if member_name == "PIR_data.f" else exact.CIR_ARCHIVE_BYTES
    expected_sha = exact.PIR_ARCHIVE_SHA256 if member_name == "PIR_data.f" else exact.CIR_ARCHIVE_SHA256
    archive_bytes = path.read_bytes()
    if len(archive_bytes) != expected_size:
        raise AssertionError("pinned archive byte length changed")
    if hashlib.sha256(archive_bytes).hexdigest() != expected_sha:
        raise AssertionError("pinned archive SHA changed")
    with zipfile.ZipFile(io.BytesIO(archive_bytes)) as archive:
        names = [info.filename for info in archive.infolist()]
        if names.count(member_name) != 1:
            raise AssertionError("Fortran source member is not unique")
        text = archive.read(member_name).decode("ascii", errors="strict")
    lower = text.lower()
    start = lower.index(marker_name.lower())
    section = text[start:]
    match = re.search(
        r"data\s+centeringmatrix/(.*?)data\s+centeringmatrix_denom/([^/]*)/",
        section,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        raise AssertionError("official centering DATA statements not found")
    matrix = tuple(
        int(token)
        for token in re.findall(r"(?<![A-Za-z])[+-]?\d+(?![A-Za-z])", match.group(1))
    )
    denominators = tuple(
        int(token)
        for token in re.findall(r"(?<![A-Za-z])[+-]?\d+(?![A-Za-z])", match.group(2))
    )
    return matrix, denominators


class _FakeProvenance:
    def __init__(self, operations, *, parent_halls=(1,)):
        self.operations = tuple(operations)
        self.parent_halls = frozenset(parent_halls)

    def spacegroup_number_for_hall(self, hall):
        return 1 if hall in self.parent_halls else 2

    def spg_operations(self, hall):
        return self.operations


def _fake_operation(rotation, translation=(Fraction(0),) * 3):
    return SimpleNamespace(rotation=rotation, translation=translation)


def _fake_source(operations, *, symbol="P1", spacegroup=1):
    return SimpleNamespace(
        spacegroup=spacegroup,
        space_group_symbol=symbol,
        centering=SimpleNamespace(value=symbol[0]),
        operations=tuple(operations),
    )


class DataHallDerivationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # This is the only full derivation in this focused process.
        cls.result = hall.derive_data_hall_frames()
        cls.source_database = exact.load_exact_iso_irrep_sources()
        cls.provenance = hall.spglib_magnetic_provenance.load_committed_provenance()

    def test_full_census_and_ordered_frames(self):
        census = self.result.census
        self.assertEqual(len(self.result.frames), 230)
        self.assertEqual(census.pir_records, 10_294)
        self.assertEqual(census.cir_records, 11_202)
        self.assertEqual(census.source_representatives, 2_609)
        self.assertEqual((census.raw_unique, census.raw_ambiguous, census.raw_missing), (220, 10, 0))
        self.assertEqual(
            census.raw_ambiguous_spacegroups,
            (5, 8, 9, 12, 15, 21, 38, 39, 65, 67),
        )
        self.assertEqual(
            (census.filtered_unique, census.filtered_ambiguous, census.filtered_missing),
            (230, 0, 0),
        )
        self.assertEqual(census.selected_hall_operations, 4_425)
        self.assertEqual(census.source_to_hall, 2_609)
        self.assertEqual(census.source_to_hall_nonzero, 0)
        self.assertEqual(census.hall_to_source, 4_425)
        self.assertEqual(census.hall_to_source_nonzero, 1_816)
        self.assertEqual(census.expanded_normalization_nonzero, 410)
        self.assertEqual(
            census.centering_counts,
            (("P", 149), ("A", 4), ("B", 0), ("C", 16),
             ("F", 16), ("I", 38), ("R", 7)),
        )
        self.assertEqual(
            sum(count for shift, count in census.hall_to_source_shifts
                if shift != (0, 0, 0)),
            1_816,
        )
        self.assertEqual(
            sum(count for shift, count in census.expanded_normalization_shifts
                if shift != (0, 0, 0)),
            410,
        )
        self.assertEqual(
            census.hall_to_source_cosets,
            (((0, 0, 0), 2_609), ((0, 6, 6), 392), ((4, 8, 8), 51),
             ((6, 0, 6), 376), ((6, 6, 0), 458), ((6, 6, 6), 488),
             ((8, 4, 4), 51)),
        )

    def test_fixed_witnesses_and_lookup(self):
        pir_first = self.source_database.pir_records[0]
        cir_first = self.source_database.cir_records[0]
        for record in (pir_first, cir_first):
            self.assertEqual(record.irnumber, 1)
            self.assertEqual(record.spacegroup, 1)
            self.assertEqual(record.operations[0].rotation, hall.IDENTITY_ROTATION)
            self.assertEqual(record.operations[0].translation, (Fraction(0),) * 3)

        sg1 = self.result.source_frame(1)
        self.assertEqual(sg1.data_hall, 1)
        self.assertEqual(sg1.raw_candidate_halls, (1,))
        self.assertEqual(sg1.source_to_hall[0].shift_numerator, (0, 0, 0))
        sg5 = self.result.source_frame(5)
        self.assertEqual(sg5.raw_candidate_halls, (9, 10, 11))
        self.assertEqual(sg5.data_hall, 9)
        self.assertEqual(sg5.centering, "C")
        self.assertEqual(self.result.source_frame(146).centering, "R")
        self.assertTrue(self.result.source_frame(146).source_symbol.startswith("R"))
        self.assertEqual(self.result.source_frame(225).centering, "F")
        self.assertTrue(self.result.source_frame(225).source_symbol.startswith("F"))
        for bad in (0, 231, True, False, 1.0, "1"):
            with self.subTest(spacegroup=bad):
                with self.assertRaises(hall.DataHallLookupError):
                    self.result.source_frame(bad)
        self.assertIs(self.result.source_frame(1), self.result.frames[0])

    def test_fortran_centering_matrices_and_independent_cosets(self):
        parsed = []
        for member_name, marker in (
            ("PIR_data.f", "subroutine pir_data_test_multiplication_table"),
            ("CIR_data.f", "subroutine cir_data_test_multiplication_table"),
        ):
            parsed.append(_fortran_centering(member_name, marker))
        self.assertEqual(parsed[0], parsed[1])
        expected_matrix = tuple(item[2] for item in hall.CENTERING_MATRIX_DATA)
        expected_denominators = tuple(item[1] for item in hall.CENTERING_MATRIX_DATA)
        self.assertEqual(
            parsed[0],
            (tuple(value for matrix in expected_matrix for value in matrix),
             expected_denominators),
        )
        for (name, denominator, matrix), (_, cosets) in zip(
            hall.CENTERING_MATRIX_DATA, hall.CENTERING_COSETS
        ):
            independent = _independent_cosets(matrix, denominator)
            self.assertEqual(cosets, independent)
            self.assertEqual(
                tuple(
                    tuple(value * hall.TRANSLATION_DENOMINATOR
                          for value in vector)
                    for vector in independent
                ),
                tuple(
                    tuple(value * hall.TRANSLATION_DENOMINATOR
                          for value in vector)
                    for vector in dict(hall.CENTERING_COSETS)[name]
                ),
            )

    def test_candidate_operation_permutation_preserves_raw_candidate_set(self):
        source = self.source_database.source_universe(1)
        normal, _ = hall._raw_candidates(source, self.provenance)

        class Permuted:
            def spacegroup_number_for_hall(self, number):
                return self_outer.provenance.spacegroup_number_for_hall(number)

            def spg_operations(self, number):
                return tuple(reversed(self_outer.provenance.spg_operations(number)))

        self_outer = self
        permuted, _ = hall._raw_candidates(source, Permuted())
        self.assertEqual(set(normal), set(permuted))
        self.assertEqual(normal, permuted)

    def test_sg5_centering_filter_is_the_only_tie_break(self):
        source = self.source_database.source_universe(5)
        candidates, details = hall._raw_candidates(source, self.provenance)
        self.assertEqual(candidates, (9, 10, 11))
        self.assertEqual(
            hall._filtered_candidates(candidates, details, hall._centering_cosets_for("C")),
            (9,),
        )
        self.assertEqual(
            hall._filtered_candidates(candidates, details, hall._centering_cosets_for("A")),
            (10,),
        )
        self.assertEqual(
            hall._filtered_candidates(candidates, details, hall._centering_cosets_for("I")),
            (11,),
        )
        for bad in ("X", "", None, []):
            with self.subTest(centering=bad):
                with self.assertRaises((hall.DataHallInvariantError, hall.DataHallSchemaError)):
                    hall._centering_cosets_for(bad)
        source = self.source_database.source_universe(5)
        candidates, details = hall._raw_candidates(source, self.provenance)
        with self.assertRaises(hall.DataHallDerivationError):
            hall._filtered_candidates(candidates, details, ())

    def test_synthetic_missing_ambiguous_and_nonintegral_fail_closed(self):
        identity = _fake_operation(hall.IDENTITY_ROTATION)
        source = _fake_source((identity,))
        missing = _FakeProvenance((_fake_operation(((-1, 0, 0), (0, -1, 0), (0, 0, -1))),))
        with self.assertRaises(hall.DataHallDerivationError):
            hall._derive_one(source, missing)

        ambiguous = _FakeProvenance((identity, identity))
        with self.assertRaises(hall.DataHallDerivationError):
            hall._raw_candidates(source, ambiguous)

        with self.assertRaises(hall.DataHallDerivationError):
            hall._mapping_for_selected_hall(
                hall._source_operations(source),
                (_fake_operation(((-1, 0, 0), (0, -1, 0), (0, 0, -1))),),
                hall._centering_cosets_for("P"),
                "synthetic missing mapping",
            )

        with self.assertRaises(hall.DataHallDerivationError):
            hall._numerators_over_12(
                (Fraction(1, 5), Fraction(0), Fraction(0)),
                "synthetic nonintegral shift",
            )

    def test_synthetic_nonclosed_source_quotient_fails(self):
        shear = ((1, 1, 0), (0, 1, 0), (0, 0, 1))
        source = _fake_source((_fake_operation(hall.IDENTITY_ROTATION), _fake_operation(shear)))
        provenance = _FakeProvenance(tuple(source.operations))
        with self.assertRaises(hall.DataHallInvariantError):
            hall._derive_one(source, provenance)

    def test_public_graph_is_frozen_and_has_no_mutable_nested_values(self):
        def check(value):
            self.assertNotIsInstance(value, (list, dict, set))
            if is_dataclass(value):
                self.assertFalse(hasattr(value, "__dict__"))
                for field in fields(value):
                    check(getattr(value, field.name))
            elif isinstance(value, tuple):
                for item in value:
                    check(item)

        check(self.result)
        for value in (
            self.result,
            self.result.frames[0],
            self.result.frames[0].source_to_hall[0],
            self.result.frames[0].hall_to_source[0],
            self.result.census,
        ):
            self.assertTrue(is_dataclass(value))
            first = fields(value)[0].name
            with self.assertRaises(FrozenInstanceError):
                setattr(value, first, None)

    def test_module_is_pure_and_has_no_runtime_fallback(self):
        text = Path(hall.__file__).read_text(encoding="utf-8")
        self.assertNotIn("generated_data", text)
        self.assertNotIn("SG_DATA_HALL", text)
        self.assertNotIn("hall_operations.json", text)
        self.assertNotIn("read_bytes", text)
        self.assertNotIn("write_outputs", text)
        self.assertNotIn("open(", text)


if __name__ == "__main__":
    unittest.main()
