"""Focused tests for the strict ISO-IR PIR/CIR source-frame loader."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import hashlib
import io
import tempfile
import threading
import unittest
from unittest import mock
from pathlib import Path
import zipfile

from . import iso_irrep_exact as exact


def _synthetic_source(archive, *, special=True, symbol="P1", label="GM1", irnumber=1,
                      matrix=None, k_payload=None, operation=None, operations=None,
                      irtranslation=None, irtranslations=None, dimension=1,
                      kcount=1, pmkcount=1, opcount=1, extra_lines=()):
    """Build one official-shape source record for parser seam tests."""

    titles = exact._PIR_TITLES if archive is exact.SourceArchive.PIR else exact._CIR_TITLES
    if matrix is None:
        matrix_token = "1" if archive is exact.SourceArchive.PIR else "(1,0)"
        matrix = " ".join(matrix_token for _ in range(dimension * dimension))
    if k_payload is None:
        if special:
            k_payload = (0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
        else:
            k_payload = (0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1)
    if operation is None:
        operation = (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
    if operations is None:
        operations = (operation,) * opcount
    if irtranslation is None:
        irtranslation = (0, 0, 0, 1)
    if irtranslations is None:
        irtranslations = (irtranslation,) * opcount
    label_width = 8 if archive is exact.SourceArchive.PIR else 4
    header = (
        f'{irnumber:5d}{1:4d} "{symbol:<10}" "{label:<{label_width}}"'
        f'{dimension:3d}{1:3d}{kcount:3d}{pmkcount:3d}{opcount:3d}'
    )
    lines = list(titles) + [
        header,
        " ".join(str(value) for value in k_payload),
    ]
    for operation_row, irtranslation_row in zip(operations, irtranslations):
        lines.append(" ".join(str(value) for value in operation_row))
        if not special:
            lines.append(" ".join(str(value) for value in irtranslation_row))
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
            exact.parse_exact_source_text("\n".join(lines) + "\n", exact.SourceArchive.PIR)
        lines = good.splitlines()
        lines[4] += " 0"  # terminal k line overrun
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_text("\n".join(lines) + "\n", exact.SourceArchive.PIR)
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

    def test_official_header_widths_and_ascii_line_contract(self):
        good = _synthetic_source(exact.SourceArchive.PIR)
        self.assertEqual(len(good.splitlines()[3]), 48)
        self.assertEqual(
            len(_synthetic_source(exact.SourceArchive.CIR).splitlines()[3]), 44
        )
        for malformed in (
            good[:-1],
            good.replace("\n", "\r\n"),
            good.replace("\n", "\t\n", 1),
            good.replace("P1", "P\u00a01", 1),
            good.replace("P1", "P\x001", 1),
            good.replace("\n    1", "\n     1", 1),
            good.replace("\n    1", "\n   01", 1),
        ):
            with self.assertRaises(exact.IsoIrrepExactError):
                exact.parse_exact_source_text(malformed, exact.SourceArchive.PIR)
        header = good.splitlines()[3]
        for token in ("+1", "01", "-0", "１"):
            # Keep a five-character i5 field while changing only its token.
            changed = header[:0] + (" " * max(0, 5 - len(token))) + token + header[5:]
            with self.assertRaises(exact.IsoIrrepExactError):
                exact.parse_exact_source_text(
                    "\n".join((good.splitlines()[:3] + [changed] + good.splitlines()[4:])) + "\n",
                    exact.SourceArchive.PIR,
                )
        lines = good.split("\n")[:-1]
        lines[3] += "\n"
        with self.assertRaises(exact.SourceSchemaError):
            exact.parse_exact_source_lines(lines, exact.SourceArchive.PIR)
        header = good.splitlines()[3]
        leading_symbol = header[:11] + " " + header[12:]
        leading_label = header[:24] + " " + header[25:]
        for malformed_header in (leading_symbol, leading_label):
            malformed_lines = good.splitlines()
            malformed_lines[3] = malformed_header
            with self.assertRaises(exact.SourceSchemaError):
                exact.parse_exact_source_text(
                    "\n".join(malformed_lines) + "\n",
                    exact.SourceArchive.PIR,
                )

    def test_payload_integer_paths_reject_noncanonical_ascii(self):
        for line_index in (4, 5):
            for token in ("+1", "01", "-0", "１"):
                lines = _synthetic_source(exact.SourceArchive.PIR).splitlines()
                values = lines[line_index].split(" ")
                for position, value in enumerate(values):
                    if value:
                        values[position] = token
                        break
                lines[line_index] = " ".join(values)
                with self.assertRaises(exact.SourceSchemaError):
                    exact.parse_exact_source_text("\n".join(lines) + "\n", exact.SourceArchive.PIR)

    def test_crystallographic_operation_and_dimension_invariants(self):
        identity = (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        determinant_zero = (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
        with self.assertRaises(exact.SourceInvariantError):
            exact.parse_exact_source_text(
                _synthetic_source(exact.SourceArchive.PIR, operation=determinant_zero),
                exact.SourceArchive.PIR,
            )
        duplicated = (1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        with self.assertRaises(exact.SourceInvariantError):
            exact.parse_exact_source_text(
                _synthetic_source(
                    exact.SourceArchive.PIR,
                    special=True,
                    operations=(identity, duplicated),
                    opcount=2,
                ),
                exact.SourceArchive.PIR,
            )
        with self.assertRaises(exact.SourceInvariantError):
            exact.parse_exact_source_text(
                _synthetic_source(
                    exact.SourceArchive.PIR,
                    special=True,
                    operations=(duplicated, identity),
                    opcount=2,
                ),
                exact.SourceArchive.PIR,
            )
        shear = (1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        with self.assertRaises(exact.SourceInvariantError):
            exact.parse_exact_source_text(
                _synthetic_source(
                    exact.SourceArchive.PIR,
                    special=True,
                    operations=(identity, shear),
                    opcount=2,
                ),
                exact.SourceArchive.PIR,
            )
        k_arm = (0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
        with self.assertRaises(exact.SourceInvariantError):
            exact.parse_exact_source_text(
                _synthetic_source(
                    exact.SourceArchive.PIR,
                    k_payload=k_arm + k_arm,
                    dimension=3,
                    kcount=2,
                    pmkcount=2,
                    matrix=" ".join("1" for _ in range(9)),
                ),
                exact.SourceArchive.PIR,
            )

        # The closure memo is keyed by the complete ordered rotation
        # signature, not merely by SG number.  A second same-SG source with a
        # different non-closed set must not inherit the first record's result.
        first = _synthetic_source(exact.SourceArchive.PIR, label="GM1")
        second = _synthetic_source(
            exact.SourceArchive.PIR,
            irnumber=2,
            label="GM2",
            operations=(identity, shear),
            opcount=2,
        )
        combined = first + "\n".join(second.splitlines()[3:]) + "\n"
        with self.assertRaises(exact.SourceInvariantError):
            exact.parse_exact_source_text(combined, exact.SourceArchive.PIR)


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
        identity = (1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        inversion = (-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1)
        altered = exact.parse_exact_source_text(
            _synthetic_source(
                exact.SourceArchive.CIR,
                symbol="C2",
                operations=(identity, inversion),
                opcount=2,
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
        with self.assertRaises(TypeError):
            replace(record, operations=list(record.operations))
        universe = exact.ExactSpaceGroupUniverse(
            1,
            "P1",
            exact.Centering.P,
            record.operations,
            (record.irnumber,),
            (record.irnumber,),
        )
        with self.assertRaises(TypeError):
            exact.ExactIsoIrrepDatabase([record], (), (None,) * 231)
        with self.assertRaises(TypeError):
            exact.ExactIsoIrrepDatabase((), (), [None] * 231)
        with self.assertRaises(TypeError):
            replace(universe, operations=list(universe.operations))

    def test_public_graph_rejects_raw_and_decoded_mismatches(self):
        special = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.PIR), exact.SourceArchive.PIR
        )[0]
        nonspecial = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.PIR, special=False), exact.SourceArchive.PIR
        )[0]
        zero = (Fraction(0),) * 3
        identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

        # A raw zero augmented row has no positive denominator, even though
        # the decoded/public fields below look like the identity.
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactSeitz(identity, zero, (0,) * 16)
        operation_raw = list(special.operations[0].raw_augmented)
        operation_raw[3] = 1
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactSeitz(
                special.operations[0].rotation,
                special.operations[0].translation,
                tuple(operation_raw),
            )

        arm_raw = list(special.k_arms[0].raw_augmented)
        arm_raw[0] = 1
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactKArm(
                special.k_arms[0].constant,
                special.k_arms[0].parameters,
                tuple(arm_raw),
            )
        parameter_raw = list(special.k_arms[0].raw_augmented)
        parameter_raw[4] = 1
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactKArm(
                special.k_arms[0].constant,
                special.k_arms[0].parameters,
                tuple(parameter_raw),
            )

        translation = nonspecial.irtranslations[0]
        translation_raw = list(translation.raw)
        translation_raw[0] = 1
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactIrTranslation(translation.vector, tuple(translation_raw))

        with self.assertRaises(exact.SourceInvariantError):
            replace(special, spacegroup=0)
        with self.assertRaises(exact.SourceInvariantError):
            replace(special, centering=exact.Centering.C)
        with self.assertRaises(exact.SourceInvariantError):
            replace(special, irtranslations=(translation,))
        with self.assertRaises(exact.SourceInvariantError):
            replace(nonspecial, irtranslations=(None,))

        inversion_raw = (-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1)
        inversion = exact.ExactSeitz(
            ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            zero,
            inversion_raw,
        )
        with self.assertRaises(exact.SourceInvariantError):
            replace(special, operations=(inversion,), irtranslations=(None,))

    def test_public_universe_and_database_slots_are_semantic(self):
        record = exact.parse_exact_source_text(
            _synthetic_source(exact.SourceArchive.PIR), exact.SourceArchive.PIR
        )[0]
        universe = exact.ExactSpaceGroupUniverse(
            1,
            "P1",
            exact.Centering.P,
            record.operations,
            (1,),
            (),
        )
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactSpaceGroupUniverse(
                1,
                "P1",
                exact.Centering.P,
                record.operations,
                (),
                (),
            )
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactSpaceGroupUniverse(
                1,
                "P1",
                exact.Centering.P,
                record.operations,
                (1, 1),
                (),
            )
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactSpaceGroupUniverse(
                1,
                "P1",
                exact.Centering.C,
                record.operations,
                (1,),
                (),
            )
        slots = [None] * 231
        slots[2] = universe
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactIsoIrrepDatabase((), (), tuple(slots))

        # A one-sided graph is a valid synthetic partial database when its
        # universe and record references agree exactly.
        slots = [None] * 231
        slots[1] = universe
        database = exact.ExactIsoIrrepDatabase((record,), (), tuple(slots))
        self.assertIs(database.universes[1], universe)

        # A universe must not dangle without records, and records must not
        # dangle without a universe slot.
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactIsoIrrepDatabase((), (), tuple(slots))
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactIsoIrrepDatabase((record,), (), (None,) * 231)

        for altered_universe in (
            replace(universe, space_group_symbol="P2"),
            replace(universe, pir_irnumbers=(2,)),
        ):
            altered_slots = [None] * 231
            altered_slots[1] = altered_universe
            with self.assertRaises(exact.SourceInvariantError):
                exact.ExactIsoIrrepDatabase((record,), (), tuple(altered_slots))

        inversion = exact.ExactSeitz(
            ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            (Fraction(0),) * 3,
            (-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1),
        )
        altered_universe = replace(
            universe,
            operations=(record.operations[0], inversion),
        )
        altered_slots = [None] * 231
        altered_slots[1] = altered_universe
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactIsoIrrepDatabase((record,), (), tuple(altered_slots))

        # The crosslink pass also compares every same-SG record, rather than
        # accepting the first record as an unchecked baseline.
        altered_record = replace(record, irrep_label="GM2", space_group_symbol="P2")
        with self.assertRaises(exact.SourceInvariantError):
            exact.ExactIsoIrrepDatabase(
                (record, altered_record), (), tuple(slots)
            )

    def test_source_lines_snapshot_and_input_boundary(self):
        source_lines = _synthetic_source(exact.SourceArchive.PIR).splitlines()
        self.assertEqual(
            exact.parse_exact_source_lines(source_lines, exact.SourceArchive.PIR)[0].irnumber,
            1,
        )

        class MutatingList(list):
            def __iter__(self):
                self[3] = "mutated while iterating"
                return super().__iter__()

        with self.assertRaises(TypeError):
            exact.parse_exact_source_lines(
                MutatingList(source_lines), exact.SourceArchive.PIR
            )
        with self.assertRaises(TypeError):
            exact.parse_exact_source_lines(
                iter(source_lines), exact.SourceArchive.PIR
            )

        original_lines = list(source_lines)
        expected_snapshot = tuple(original_lines)
        captured = []

        def capture(snapshot, archive, *, validate_census, ascii_validated):
            captured.append(snapshot)
            original_lines[3] = "mutated after snapshot"
            return ()

        with mock.patch.object(
            exact, "_parse_records_from_lines", side_effect=capture
        ):
            self.assertEqual(
                exact.parse_exact_source_lines(
                    original_lines, exact.SourceArchive.PIR
                ),
                (),
            )
        self.assertEqual(captured, [expected_snapshot])
        self.assertIs(type(captured[0]), tuple)

    def test_single_flight(self):
        old_database = exact._DATABASE
        calls = []
        original = exact._load_uncached
        original_lock = exact._DATABASE_LOCK
        marker = exact.ExactIsoIrrepDatabase(( ), ( ), (None,) * 231)
        exact._DATABASE = None
        builder_started = threading.Event()
        second_attempt = threading.Event()

        class TrackingLock:
            def __init__(self):
                self.lock = threading.Lock()
                self.attempts = 0

            def __enter__(self):
                self.attempts += 1
                if self.attempts == 2:
                    second_attempt.set()
                self.lock.acquire()
                return self

            def __exit__(self, _kind, _value, _traceback):
                self.lock.release()

        tracking_lock = TrackingLock()
        exact._DATABASE_LOCK = tracking_lock

        def build():
            calls.append(True)
            builder_started.set()
            if not second_attempt.wait(5):
                raise AssertionError("second caller never attempted the cache lock")
            return marker

        exact._load_uncached = build
        try:
            outputs = []
            errors = []

            def worker():
                try:
                    outputs.append(exact.load_exact_iso_irrep_sources())
                except BaseException as error:  # make thread failures visible
                    errors.append(error)

            first = threading.Thread(target=worker)
            first.start()
            self.assertTrue(builder_started.wait(2))
            second = threading.Thread(target=worker)
            second.start()
            self.assertTrue(second_attempt.wait(2))
            first.join(5)
            second.join(5)
            self.assertFalse(first.is_alive())
            self.assertFalse(second.is_alive())
            self.assertFalse(errors)
            self.assertEqual(len(calls), 1)
            self.assertEqual(len(outputs), 2)
            self.assertTrue(all(item is marker for item in outputs))

            exact._DATABASE = None

            def fail_once():
                raise RuntimeError("synthetic builder failure")

            exact._load_uncached = fail_once
            with self.assertRaises(RuntimeError):
                exact.load_exact_iso_irrep_sources()
            self.assertIsNone(exact._DATABASE)
            exact._load_uncached = lambda: marker
            self.assertIs(exact.load_exact_iso_irrep_sources(), marker)
        finally:
            exact._load_uncached = original
            exact._DATABASE_LOCK = original_lock
            exact._DATABASE = old_database


class ArchiveTests(unittest.TestCase):
    @staticmethod
    def _zip_bytes(members, compression=zipfile.ZIP_STORED):
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, "w", compression=compression) as archive:
            for name, contents in members:
                archive.writestr(name, contents)
        return stream.getvalue()

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

    def test_verified_member_uses_one_immutable_archive_snapshot(self):
        good = self._zip_bytes(
            [("PIR_data.txt", _synthetic_source(exact.SourceArchive.PIR).encode("ascii"))]
        )
        evil = self._zip_bytes([("PIR_data.txt", b"evil\n")])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "PIR_data.zip"
            path.write_bytes(evil)
            with mock.patch.object(Path, "read_bytes", return_value=good) as read_bytes:
                text = exact._read_verified_member(
                    exact.SourceArchive.PIR,
                    path=path,
                    expected_bytes=len(good),
                    expected_sha256=hashlib.sha256(good).hexdigest(),
                )
            self.assertEqual(text, _synthetic_source(exact.SourceArchive.PIR))
            self.assertEqual(read_bytes.call_count, 1)

    def test_root_member_name_and_archive_errors_are_typed(self):
        bad_member_sets = (
            (("nested/PIR_data.txt", b"x"),),
            (("../PIR_data.txt", b"x"),),
            (("PIR_data.txt", b"x"), ("nested/PIR_data.txt", b"x")),
            (("PIR_data.txt", b"x"), ("PIR_data.txt", b"y")),
        )
        with tempfile.TemporaryDirectory() as directory:
            for ordinal, members in enumerate(bad_member_sets):
                path = Path(directory) / ("bad%d.zip" % ordinal)
                raw = self._zip_bytes(members)
                path.write_bytes(raw)
                with self.assertRaises(exact.ArchiveIntegrityError):
                    exact._read_verified_member(
                        exact.SourceArchive.PIR,
                        path=path,
                        expected_bytes=len(raw),
                        expected_sha256=hashlib.sha256(raw).hexdigest(),
                    )
            with self.assertRaises(exact.ArchiveIntegrityError):
                exact._read_verified_member(
                    exact.SourceArchive.PIR,
                    path=path,
                    expected_bytes=len(raw) + 1,
                    expected_sha256=hashlib.sha256(raw).hexdigest(),
                )
            with self.assertRaises(exact.ArchiveIntegrityError):
                exact._read_verified_member(
                    exact.SourceArchive.PIR,
                    path=path,
                    expected_bytes=len(raw),
                    expected_sha256="0" * 64,
                )
            utf8 = self._zip_bytes([("PIR_data.txt", b"\xff")])
            path.write_bytes(utf8)
            with self.assertRaises(exact.SourceSchemaError):
                exact._read_verified_member(
                    exact.SourceArchive.PIR,
                    path=path,
                    expected_bytes=len(utf8),
                    expected_sha256=hashlib.sha256(utf8).hexdigest(),
                )
            not_a_zip = b"not a zip"
            path.write_bytes(not_a_zip)
            with self.assertRaises(exact.ArchiveIntegrityError):
                exact._read_verified_member(
                    exact.SourceArchive.PIR,
                    path=path,
                    expected_bytes=len(not_a_zip),
                    expected_sha256=hashlib.sha256(not_a_zip).hexdigest(),
                )
            with mock.patch.object(Path, "read_bytes", side_effect=NotImplementedError):
                with self.assertRaises(exact.ArchiveIntegrityError):
                    exact._read_verified_member(exact.SourceArchive.PIR, path=path)


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
