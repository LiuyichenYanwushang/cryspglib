"""Focused tests for the strict ISO-IR PIR/CIR source-frame loader."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from fractions import Fraction
import hashlib
import tempfile
import threading
import unittest
from pathlib import Path
import zipfile

from . import iso_irrep_exact as exact


def _synthetic_source(archive, *, special=True, symbol="P1", label="GM1", irnumber=1,
                      matrix=None, k_payload=None, operation=None,
                      irtranslation=None, extra_lines=()):
    """Build one official-shape source record for parser seam tests."""

    titles = exact._PIR_TITLES if archive is exact.SourceArchive.PIR else exact._CIR_TITLES
    if matrix is None:
        matrix = "1" if archive is exact.SourceArchive.PIR else "(1,0)"
    if k_payload is None:
        if special:
            k_payload = (0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
        else:
            k_payload = (0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1)
    if operation is None:
        operation = (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
    if irtranslation is None:
        irtranslation = (0, 0, 0, 1)
    header = (
        f'    {irnumber}   1 "{symbol:<10}" "{label:<8}"  1  1  1  1  1'
    )
    lines = list(titles) + [
        header,
        " ".join(str(value) for value in k_payload),
        " ".join(str(value) for value in operation),
    ]
    if not special:
        lines.append(" ".join(str(value) for value in irtranslation))
    lines.append(matrix)
    lines.extend(extra_lines)
    return "\n".join(lines) + "\n"


class SyntheticParserTests(unittest.TestCase):
    def test_pir_and_cir_exact_payloads(self):
        for archive in (exact.SourceArchive.PIR, exact.SourceArchive.CIR):
            record = exact.parse_exact_source_text(
                _synthetic_source(archive, special=False), archive
            )[0]
            self.assertEqual(record.k_arms[0].constant, (Fraction(0),) * 3)
            self.assertEqual(record.k_arms[0].parameters[0], (Fraction(1), Fraction(0), Fraction(0)))
            self.assertEqual(record.operations[0].rotation, ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
            self.assertEqual(record.operations[0].translation, (Fraction(0),) * 3)
            self.assertEqual(record.irtranslations[0].vector, (Fraction(0),) * 3)
            self.assertFalse(record.special)

    def test_special_has_no_translation_row(self):
        record = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.PIR), exact.SourceArchive.PIR
        )[0]
        self.assertEqual(record.irtranslations, (None,))

    def test_canonical_integer_gate(self):
        for token in ("+1", "01", "-0", "١"):
            source = _synthetic_source(exact.SourceArchive.PIR).replace(
                "  1  1  1  1  1", f"  {token}  1  1  1  1", 1
            )
            with self.assertRaises(exact.SourceSchemaError):
                exact.parse_exact_source_text(source, exact.SourceArchive.PIR)

    def test_payload_short_extra_and_blank_are_fatal(self):
        good = _synthetic_source(exact.SourceArchive.PIR)
        lines = good.splitlines()
        lines[4] = "0 0 0"  # k payload short, followed by operation row
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text("\n".join(lines), exact.SourceArchive.PIR)
        lines = good.splitlines()
        lines[4] += " 0"  # terminal k line overrun
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text("\n".join(lines), exact.SourceArchive.PIR)
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text(good.replace("\n1\n", "\n\n1\n"), exact.SourceArchive.PIR)

    def test_operation_and_parameter_invariants(self):
        source = _synthetic_source(
            exact.SourceArchive.PIR,
            operation=(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1),
        )
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text(source, exact.SourceArchive.PIR)
        bad_k = (0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1)
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text(
                _synthetic_source(exact.SourceArchive.PIR, special=False, k_payload=bad_k),
                exact.SourceArchive.PIR,
            )

    def test_unknown_matrix_and_centering(self):
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text(
                _synthetic_source(exact.SourceArchive.PIR, matrix="1.00000"),
                exact.SourceArchive.PIR,
            )
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text(
                _synthetic_source(exact.SourceArchive.PIR, symbol="X1"),
                exact.SourceArchive.PIR,
            )


class UniverseAndCacheTests(unittest.TestCase):
    def test_universe_folding_compares_cross_archive_operations(self):
        pir = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.PIR, symbol="C2"), exact.SourceArchive.PIR
        )
        cir = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.CIR, symbol="C2"), exact.SourceArchive.CIR
        )
        database = exact._assemble_database(pir, cir, validate_census=False)
        self.assertEqual(database.universes[1].centering, exact.Centering.C)
        self.assertEqual(database.universes[1].pir_irnumbers, (1,))
        self.assertEqual(database.universes[1].cir_irnumbers, (1,))
        altered = exact.parse_exact_source_text(
            _synthetic_source(
                exact.SourceArchive.CIR,
                symbol="C2",
                operation=(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1),
            ),
            exact.SourceArchive.CIR,
        )
        with self.assertRaises(exact.SourceInvariantError):
            exact._assemble_database(pir, altered, validate_census=False)

    def test_lookup_is_strict(self):
        database = exact.ExactIsoIrrepDatabase(( ), ( ), (None,) * 231)
        for value in (True, False, 1.0, "1", None):
            with self.assertRaises(exact.SourceLookupError):
                database.source_universe(value)
        for value in (0, -1, 231, 10**100):
            with self.assertRaises(exact.SourceLookupError):
                database.source_universe(value)

    def test_frozen_slots_objects(self):
        record = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.PIR), exact.SourceArchive.PIR
        )[0]
        self.assertFalse(hasattr(record, "__dict__"))
        with self.assertRaises(FrozenInstanceError):
            record.dimension = 2
        with self.assertRaises(TypeError):
            record.k_arms[0].raw_augmented[0] = 1

    def test_single_flight(self):
        old_database = exact._DATABASE
        calls = []
        original = exact._load_uncached
        marker = exact.ExactIsoIrrepDatabase(( ), ( ), (None,) * 231)
        exact._DATABASE = None
        exact._load_uncached = lambda: (calls.append(True) or marker)
        try:
            outputs = []
            threads = [threading.Thread(target=lambda: outputs.append(exact.load_exact_iso_irrep_sources())) for _ in range(12)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
            self.assertEqual(len(calls), 1)
            self.assertEqual(len(outputs), 12)
            self.assertTrue(all(item is marker for item in outputs))
        finally:
            exact._load_uncached = original
            exact._DATABASE = old_database


class ArchiveTests(unittest.TestCase):
    def test_pinned_archive_sizes_and_hashes(self):
        for archive, path, size, digest in (
            (exact.SourceArchive.PIR, exact.PIR_ARCHIVE_PATH, exact.PIR_ARCHIVE_BYTES, exact.PIR_ARCHIVE_SHA256),
            (exact.SourceArchive.CIR, exact.CIR_ARCHIVE_PATH, exact.CIR_ARCHIVE_BYTES, exact.CIR_ARCHIVE_SHA256),
        ):
            self.assertEqual(path.stat().st_size, size)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
            self.assertTrue(exact._read_verified_member(archive).startswith("ISO-IR:"))

    def test_member_ambiguity_is_fatal(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ambiguous.zip"
            with zipfile.ZipFile(str(path), "w") as archive:
                archive.writestr("PIR_data.txt", "x")
                archive.writestr("nested/PIR_data.txt", "x")
            with self.assertRaises(exact.ArchiveIntegrityError):
                exact._read_verified_member(
                    exact.SourceArchive.PIR,
                    path=path,
                    expected_bytes=path.stat().st_size,
                    expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                )


class FullSourceCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = exact.load_exact_iso_irrep_sources()

    def test_full_census_and_universes(self):
        database = self.database
        self.assertEqual(len(database.pir_records), 10_294)
        self.assertEqual(len(database.cir_records), 11_202)
        self.assertEqual(len(database.pir_records) + len(database.cir_records), 21_496)
        self.assertEqual(sum(record.special for record in database.pir_records + database.cir_records), 10_073)
        self.assertEqual(
            sum(not record.special for record in database.pir_records + database.cir_records),
            11_423,
        )
        self.assertEqual(
            sum(translation is not None for record in database.pir_records for translation in record.irtranslations),
            64_588,
        )
        self.assertEqual(
            sum(translation is not None for record in database.cir_records for translation in record.irtranslations),
            68_612,
        )
        self.assertEqual(sum(universe is not None for universe in database.universes), 230)
        self.assertEqual(sum(len(universe.operations) for universe in database.universes[1:] if universe), 2_609)
        self.assertEqual(
            {operation.raw_augmented[15] for record in database.pir_records for operation in record.operations},
            {1, 2, 3, 4, 6},
        )
        self.assertEqual(
            {translation.raw[3] for record in database.pir_records for translation in record.irtranslations if translation},
            {1, 2, 3, 4, 6},
        )
        self.assertEqual(
            {operation.raw_augmented[15] for record in database.cir_records for operation in record.operations},
            {1, 2, 3, 4, 6},
        )
        self.assertEqual(
            {translation.raw[3] for record in database.cir_records for translation in record.irtranslations if translation},
            {1, 2, 3, 4, 6},
        )
        self.assertEqual(
            {universe.centering.value for universe in database.universes[1:] if universe},
            {"P", "A", "C", "F", "I", "R"},
        )

    def test_fixed_witnesses_and_cache_identity(self):
        database = self.database
        self.assertIs(database, exact.load_exact_iso_irrep_sources())
        for records in (database.pir_records, database.cir_records):
            self.assertEqual(records[0].spacegroup, 1)
            self.assertEqual(records[0].operations[0].rotation, ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
            self.assertEqual(records[0].operations[0].translation, (Fraction(0),) * 3)
            self.assertFalse(records[8].special)
            self.assertTrue(all(translation is not None for translation in records[8].irtranslations))
        self.assertEqual(database.source_universe(5).space_group_symbol, "C2")
        self.assertEqual(database.source_universe(5).centering, exact.Centering.C)
        self.assertEqual(database.source_universe(146).centering, exact.Centering.R)
        self.assertEqual(database.source_universe(225).centering, exact.Centering.F)


if __name__ == "__main__":
    unittest.main()
