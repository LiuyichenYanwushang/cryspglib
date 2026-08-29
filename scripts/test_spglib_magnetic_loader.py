#!/usr/bin/env python3
"""Focused tests for the typed committed magnetic provenance loader."""

from dataclasses import FrozenInstanceError, fields, is_dataclass
from fractions import Fraction
import hashlib
import os
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest import mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import extract_spglib_magnetic_provenance as extractor
import spglib_magnetic_provenance as loader


ARTIFACT = Path(__file__).parent / "data/spglib_magnetic_provenance_v1.json"
MANIFEST = Path(__file__).parent / "data/spglib_magnetic_provenance_v1.manifest.json"
GOLDEN_ARTIFACT_LENGTH = 1_537_875
GOLDEN_ARTIFACT_SHA256 = (
    "933a52a6696e7f6a1a2e426825ad92c377c6e96330e18c5c045d659798d740b9"
)
GOLDEN_MANIFEST_LENGTH = 570
GOLDEN_MANIFEST_SHA256 = (
    "6a9e1b64c190c30a556d63e51e5b896b967d33e8821714beb745ae699fab84bf"
)


class MagneticLoaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = loader.load_committed_provenance()

    @staticmethod
    def _synchronized_pair(mutator):
        artifact_bytes = ARTIFACT.read_bytes()
        manifest_bytes = MANIFEST.read_bytes()
        artifact = extractor._parse_json_bytes(artifact_bytes, "artifact")
        manifest = extractor._parse_json_bytes(manifest_bytes, "manifest")
        mutator(artifact)
        artifact_bytes = extractor.canonical_json(artifact)
        manifest["artifact"]["bytes"] = len(artifact_bytes)
        manifest["artifact"]["sha256"] = hashlib.sha256(artifact_bytes).hexdigest()
        manifest_bytes = extractor.canonical_json(manifest)
        return artifact_bytes, manifest_bytes

    def _restore_cached_database(self):
        with loader._CACHE_LOCK:
            loader._CACHED_DATABASE = self.database

    def test_fixed_trust_root_and_cached_immutable_result(self):
        artifact = ARTIFACT.read_bytes()
        manifest = MANIFEST.read_bytes()
        self.assertEqual(len(artifact), GOLDEN_ARTIFACT_LENGTH)
        self.assertEqual(hashlib.sha256(artifact).hexdigest(), GOLDEN_ARTIFACT_SHA256)
        self.assertEqual(len(manifest), GOLDEN_MANIFEST_LENGTH)
        self.assertEqual(hashlib.sha256(manifest).hexdigest(), GOLDEN_MANIFEST_SHA256)
        self.assertIs(self.database, loader.load_committed_provenance())

    def test_relative_import_freezes_absolute_data_path(self):
        repo = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as other_directory:
            probe = "\n".join((
                "import os, sys",
                "sys.path.insert(0, 'scripts')",
                "import spglib_magnetic_provenance as loader",
                "assert not os.path.isabs(loader.__file__), loader.__file__",
                "assert loader._DATA_DIR.is_absolute(), loader._DATA_DIR",
                "os.chdir(sys.argv[1])",
                "database = loader.load_committed_provenance()",
                "print(len(database.spg.operation_index), len(database.msg.metadata))",
            ))
            environment = os.environ.copy()
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            result = subprocess.run(
                [sys.executable, "-c", probe, other_directory],
                cwd=str(repo),
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=90,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "531 1652")

    def test_loader_reads_each_committed_file_once(self):
        artifact = ARTIFACT.read_bytes()
        manifest = MANIFEST.read_bytes()
        loader._reset_cache_for_test()
        try:
            sentinel = object()
            with mock.patch.object(
                loader.Path, "read_bytes", side_effect=[artifact, manifest]
            ) as read_bytes:
                with mock.patch.object(loader, "_from_bytes", return_value=sentinel):
                    self.assertIs(loader.load_committed_provenance(), sentinel)
            self.assertEqual(read_bytes.call_count, 2)
        finally:
            loader._reset_cache_for_test()
            self._restore_cached_database()

    def test_concurrent_cold_load_is_single_flight(self):
        artifact = ARTIFACT.read_bytes()
        manifest = MANIFEST.read_bytes()
        loader._reset_cache_for_test()
        start = threading.Barrier(3)
        results = [None, None]
        errors = [None, None]
        reads = {}
        payloads = {ARTIFACT.name: artifact, MANIFEST.name: manifest}
        sentinel = object()

        def read_bytes(path):
            reads[path.name] = reads.get(path.name, 0) + 1
            return payloads[path.name]

        def build_database(*args, **kwargs):
            return sentinel

        def invoke(index):
            try:
                start.wait(timeout=10)
                results[index] = loader.load_committed_provenance()
            except BaseException as error:  # report thread failures below
                errors[index] = error

        threads = [threading.Thread(target=invoke, args=(index,)) for index in range(2)]
        try:
            with mock.patch.object(
                loader.Path, "read_bytes", autospec=True, side_effect=read_bytes
            ) as read_mock:
                with mock.patch.object(
                    loader, "_from_bytes", side_effect=build_database
                ) as builder_mock:
                    for thread in threads:
                        thread.start()
                    start.wait(timeout=10)
                    for thread in threads:
                        thread.join(timeout=90)
                    self.assertTrue(all(not thread.is_alive() for thread in threads))
                    self.assertEqual(read_mock.call_count, 2)
                    self.assertEqual(builder_mock.call_count, 1)
            self.assertEqual(reads, {
                ARTIFACT.name: 1,
                MANIFEST.name: 1,
            })
            self.assertEqual(errors, [None, None])
            self.assertIs(results[0], sentinel)
            self.assertIs(results[1], sentinel)
        finally:
            for thread in threads:
                if thread.is_alive():
                    thread.join(timeout=1)
            loader._reset_cache_for_test()
            self._restore_cached_database()

    def test_decoder_rejects_sentinel_limit_and_bool_with_decode_error(self):
        # A pair-level invalid encoding is rejected by the authoritative
        # extractor schema before typed conversion; these exercise the typed
        # codec classification directly.
        for encoded in (0, loader.MAGNETIC_OPERATION_ENCODING_LIMIT, True):
            with self.subTest(encoded=encoded):
                with self.assertRaises(loader.MagneticProvenanceDecodeError):
                    loader._decode_operation(encoded)

    def test_test_only_pair_seam_rejects_typed_duplicate(self):
        def mutate(artifact):
            order, offset = artifact["spg"]["symmetry_operation_index"][2]
            self.assertEqual(order, 2)
            operations = artifact["spg"]["symmetry_operations"]
            operations[offset] = operations[offset + 1]

        artifact_bytes, manifest_bytes = self._synchronized_pair(mutate)
        parsed = extractor.parse_and_validate_committed_pair(
            artifact_bytes, manifest_bytes, ARTIFACT.name
        )
        self.assertEqual(
            parsed["spg"]["symmetry_operations"][2],
            parsed["spg"]["symmetry_operations"][3],
        )
        with self.assertRaises(loader.MagneticProvenanceInvariantError):
            loader._from_uncommitted_pair_for_test(
                artifact_bytes, manifest_bytes, ARTIFACT.name
            )

    def test_test_only_pair_seam_wraps_schema_corruption(self):
        def mutate(artifact):
            artifact["spg"]["symmetry_operations"].pop()

        artifact_bytes, manifest_bytes = self._synchronized_pair(mutate)
        with self.assertRaises(loader.MagneticProvenanceSchemaError):
            loader._from_uncommitted_pair_for_test(
                artifact_bytes, manifest_bytes, ARTIFACT.name
            )

    def test_dataclasses_are_frozen_slotted_and_nested_data_are_tuples(self):
        objects = (
            self.database,
            self.database.spg,
            self.database.msg,
            self.database.msg.metadata[7],
            self.database.spg.operation_index[1],
            self.database.msg.operation_index[7][0],
        )
        for value in objects:
            self.assertTrue(is_dataclass(value))
            self.assertFalse(hasattr(value, "__dict__"))
            first_field = fields(value)[0].name
            with self.assertRaises(FrozenInstanceError):
                setattr(value, first_field, None)

        self.assertIsInstance(self.database.spg.spacegroup_numbers, tuple)
        self.assertIsInstance(self.database.msg.uni_mapping, tuple)
        self.assertIsInstance(self.database.msg.operation_index[7], tuple)
        self.assertIsInstance(self.database.msg.operation_index[7][0], loader.OperationSpan)
        self.assertIsInstance(self.database.msg.operation_index[7][0].order, int)
        with self.assertRaises(TypeError):
            self.database.msg.uni_mapping[1][0] = 2

        def assert_no_mutable(value):
            self.assertNotIsInstance(value, (list, dict))
            if is_dataclass(value):
                for field in fields(value):
                    assert_no_mutable(getattr(value, field.name))
            elif isinstance(value, tuple):
                for item in value:
                    assert_no_mutable(item)

        assert_no_mutable(self.database)

    def test_census_and_cross_table_closure_inputs(self):
        spg = self.database.spg
        msg = self.database.msg
        self.assertEqual(len(spg.spacegroup_numbers), 531)
        self.assertEqual(len(spg.operation_index), 531)
        self.assertEqual(len(spg.raw_operation_codes), 8147)
        self.assertEqual(sum(span.order for span in spg.operation_index[1:]), 7388)
        self.assertEqual(len(msg.metadata), 1652)
        self.assertEqual(len(msg.uni_mapping), 1652)
        self.assertEqual(len(msg.derived_hall_mapping), 531)
        self.assertEqual(len(msg.operation_index), 1652)
        self.assertEqual(len(msg.raw_operation_codes), 76683)
        self.assertEqual(len(msg.decoded_operations), 76683)
        self.assertEqual(sum(
            span.order
            for row in msg.operation_index[1:]
            for span in row
            if span.order
        ), 76682)
        self.assertEqual(sum(
            len(codes)
            for row in msg.alternative_codes
            for codes in row
        ), 536)
        self.assertEqual(msg.derived_hall_mapping[0], (0, 0))
        self.assertIsNone(spg.decoded_operations[0])
        self.assertIsNone(msg.decoded_operations[0])

        type_counts = {kind: 0 for kind in loader.MagneticKind}
        for metadata in msg.metadata[1:]:
            type_counts[metadata.kind] += 1
        self.assertEqual(
            {int(kind): count for kind, count in type_counts.items()},
            {1: 230, 2: 230, 3: 674, 4: 517},
        )
        for hall in range(1, 531):
            unis = self.database.unis_for_hall(hall)
            self.assertEqual(unis, tuple(range(unis[0], unis[-1] + 1)))
            for uni in unis:
                self.assertIn(hall, self.database.halls_for_uni(uni))
                self.assertEqual(
                    self.database.spacegroup_number_for_hall(hall),
                    self.database.magnetic_metadata(uni).parent_spacegroup,
                )

    def test_query_slices_and_exact_operation_properties(self):
        def determinant(rotation):
            (a, b, c), (d, e, f), (g, h, i) = rotation
            return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

        def reencode(operation):
            rotation_payload = sum(
                (value + 1) * 3 ** (8 - index)
                for index, value in enumerate(
                    value for row in operation.rotation for value in row
                )
            )
            translation_payload = sum(
                value * 12 ** (2 - index)
                for index, value in enumerate(operation.translation_numerator)
            )
            return (
                int(operation.time_reversal) * loader.SPACE_OPERATION_SCALE
                + rotation_payload
                + loader.ROTATION_PAYLOAD * translation_payload
            )

        for index, (raw, operation) in enumerate(zip(
                self.database.spg.raw_operation_codes,
                self.database.spg.decoded_operations)):
            if index == 0:
                self.assertIsNone(operation)
                continue
            self.assertEqual(reencode(operation), raw)
            self.assertIn(determinant(operation.rotation), (-1, 1))
            self.assertTrue(all(0 <= value < 12 for value in operation.translation_numerator))

        for index, (raw, operation) in enumerate(zip(
                self.database.msg.raw_operation_codes,
                self.database.msg.decoded_operations)):
            if index == 0:
                self.assertIsNone(operation)
                continue
            self.assertEqual(reencode(operation), raw)
            self.assertIn(determinant(operation.rotation), (-1, 1))
            self.assertTrue(all(0 <= value < 12 for value in operation.translation_numerator))

        for hall in range(1, 531):
            span = self.database.spg_operation_span(hall)
            raw = self.database.spg.raw_operation_codes[span.offset:span.offset + span.order]
            self.assertEqual(tuple(operation.encoded for operation in self.database.spg_operations(hall)), raw)
        for uni in range(1, 1652):
            for hall in self.database.halls_for_uni(uni):
                span = self.database.magnetic_operation_span(uni, hall)
                raw = self.database.msg.raw_operation_codes[span.offset:span.offset + span.order]
                self.assertEqual(tuple(operation.encoded for operation in self.database.magnetic_operations(uni, hall)), raw)

        self.assertEqual(
            self.database.spg_operations(1)[0].translation,
            (Fraction(0), Fraction(0), Fraction(0)),
        )

    def test_uni7_and_uni9_witnesses(self):
        uni7 = self.database.magnetic_metadata(7)
        self.assertEqual(
            (uni7.uni, uni7.litvin, uni7.bns, uni7.og,
             uni7.parent_spacegroup, uni7.kind),
            (7, 7, "2.7", "2.4.7", 2, loader.MagneticKind.ANTI_TRANSLATION),
        )
        self.assertEqual(self.database.halls_for_uni(7), (2,))
        self.assertEqual(self.database.magnetic_operation_span(7, 2), loader.OperationSpan(4, 14))
        self.assertEqual(
            tuple(operation.encoded for operation in self.database.magnetic_operations(7, 2)),
            (16484, 3198, 34146806, 34133520),
        )
        self.assertEqual(
            tuple(operation.rotation for operation in self.database.magnetic_operations(7, 2)),
            (
                ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
                ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            ),
        )
        self.assertEqual(
            tuple(operation.translation for operation in self.database.magnetic_operations(7, 2)),
            (
                (Fraction(0), Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(1, 2)),
                (Fraction(0), Fraction(0), Fraction(1, 2)),
            ),
        )
        self.assertEqual(
            tuple(operation.time_reversal for operation in self.database.magnetic_operations(7, 2)),
            (loader.TimeReversal.UNITARY, loader.TimeReversal.UNITARY,
             loader.TimeReversal.ANTIUNITARY, loader.TimeReversal.ANTIUNITARY),
        )
        self.assertEqual(self.database.raw_alternative_codes(7, 2),
                         (30, 90, 111, 810, 2301, 6831))
        transformations = self.database.std_transformations(7, 2)
        self.assertEqual(len(transformations), 7)
        self.assertEqual(transformations[0], loader._IDENTITY_OPERATION)

        uni9 = self.database.magnetic_metadata(9)
        self.assertEqual(
            (uni9.uni, uni9.litvin, uni9.bns, uni9.og,
             uni9.parent_spacegroup, uni9.kind),
            (9, 9, "3.2", "3.2.9", 3, loader.MagneticKind.GREY),
        )
        self.assertEqual(self.database.halls_for_uni(9), (3, 4, 5))
        self.assertEqual(
            tuple(self.database.magnetic_operation_span(9, hall) for hall in (3, 4, 5)),
            (loader.OperationSpan(4, 20), loader.OperationSpan(4, 40),
             loader.OperationSpan(4, 60)),
        )
        for hall in (3, 4, 5):
            operations = self.database.magnetic_operations(9, hall)
            self.assertEqual(sum(operation.time_reversal is loader.TimeReversal.UNITARY
                                 for operation in operations), 2)
            self.assertEqual(sum(operation.time_reversal is loader.TimeReversal.ANTIUNITARY
                                 for operation in operations), 2)
            self.assertIn(34028708, tuple(operation.encoded for operation in operations))
            self.assertEqual(self.database.raw_alternative_codes(9, hall), ())
            self.assertEqual(self.database.std_transformations(9, hall),
                             (loader._IDENTITY_OPERATION,))

    def test_invalid_lookups_are_structured_and_one_based(self):
        invalid_halls = (0, -1, 531, True, 1.0, "1", None)
        for hall in invalid_halls:
            with self.subTest(hall=hall):
                with self.assertRaises(loader.ArtifactLookupError):
                    self.database.spacegroup_number_for_hall(hall)
                with self.assertRaises(loader.ArtifactLookupError):
                    self.database.spg_operation_span(hall)
                with self.assertRaises(loader.ArtifactLookupError):
                    self.database.unis_for_hall(hall)
        invalid_unis = (0, -1, 1652, True, 1.0, "1", None)
        for uni in invalid_unis:
            with self.subTest(uni=uni):
                with self.assertRaises(loader.ArtifactLookupError):
                    self.database.magnetic_metadata(uni)
                with self.assertRaises(loader.ArtifactLookupError):
                    self.database.halls_for_uni(uni)
        with self.assertRaises(loader.ArtifactLookupError):
            self.database.magnetic_operation_span(7, 1)
        with self.assertRaises(loader.ArtifactLookupError):
            self.database.raw_alternative_codes(9, 6)
        with self.assertRaises(loader.ArtifactLookupError):
            self.database.std_transformations(9, 6)

    def test_private_bytes_corruption_fails_as_integrity(self):
        artifact = ARTIFACT.read_bytes()
        manifest = MANIFEST.read_bytes()
        broken_artifact = bytearray(artifact)
        broken_artifact[-2] ^= 1
        with self.assertRaises(loader.MagneticProvenanceIntegrityError):
            loader._from_bytes(bytes(broken_artifact), manifest)
        broken_manifest = bytearray(manifest)
        broken_manifest[-2] ^= 1
        with self.assertRaises(loader.MagneticProvenanceIntegrityError):
            loader._from_bytes(artifact, bytes(broken_manifest))


if __name__ == "__main__":
    unittest.main()
