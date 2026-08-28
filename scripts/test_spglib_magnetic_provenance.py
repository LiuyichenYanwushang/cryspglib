#!/usr/bin/env python3
"""Focused tests for the pinned spglib magnetic provenance extractor."""

import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import extract_spglib_magnetic_provenance as extractor


UPSTREAM = Path("/tmp/spglib-v2.5.0")
ARTIFACT = Path(__file__).parent / "data/spglib_magnetic_provenance_v1.json"
MANIFEST = Path(__file__).parent / "data/spglib_magnetic_provenance_v1.manifest.json"


class MagneticProvenanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact, cls.details = extractor.extract(UPSTREAM)

    def test_pinned_source_hashes_and_census(self):
        self.assertEqual(self.details["msg_database.c"], extractor.EXPECTED_SOURCES["msg_database.c"])
        self.assertEqual(self.details["spg_database.c"], extractor.EXPECTED_SOURCES["spg_database.c"])
        self.assertEqual(len(self.artifact["spg"]["spacegroup_number"]), 531)
        self.assertEqual(len(self.artifact["spg"]["symmetry_operation_index"]), 531)
        self.assertEqual(len(self.artifact["msg"]["magnetic_spacegroup_types"]), 1652)
        self.assertEqual(len(self.artifact["msg"]["magnetic_symmetry_operations"]), 76683)
        self.assertEqual(
            {str(kind): sum(row["type"] == kind for row in self.artifact["msg"]["magnetic_spacegroup_types"][1:])
             for kind in (1, 2, 3, 4)},
            extractor.EXPECTED_TYPE_COUNTS,
        )
        self.assertEqual(self.artifact["msg"]["magnetic_symmetry_operations"][0], 0)

    def test_normalized_dimensions_and_ranges(self):
        msg = self.artifact["msg"]
        self.assertEqual(
            {len(row) for row in msg["magnetic_spacegroup_operation_index"]}, {18}
        )
        self.assertTrue(all(len(pair) == 2 for row in msg["magnetic_spacegroup_operation_index"] for pair in row))
        self.assertEqual(
            {len(row) for row in msg["alternative_transformations"]}, {18}
        )
        self.assertTrue(all(len(entry) == 7 for row in msg["alternative_transformations"] for entry in row))
        operations = msg["magnetic_symmetry_operations"]
        for uni, mapping in enumerate(msg["magnetic_spacegroup_uni_mapping"]):
            hall_count, first_hall = mapping
            if uni == 0:
                self.assertEqual(mapping, [0, 0])
                continue
            self.assertGreaterEqual(first_hall, 1)
            for order, offset in msg["magnetic_spacegroup_operation_index"][uni][:hall_count]:
                self.assertGreater(order, 0)
                self.assertGreater(offset, 0)
                self.assertLessEqual(offset + order, len(operations))

    def test_decoder_witnesses_exact(self):
        expected = {
            16484: ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 0),
            34146806: ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 6], 1),
            3198: ([-1, 0, 0, 0, -1, 0, 0, 0, -1], [0, 0, 0], 0),
            34133520: ([-1, 0, 0, 0, -1, 0, 0, 0, -1], [0, 0, 6], 1),
            3360: ([-1, 0, 0, 0, 1, 0, 0, 0, -1], [0, 0, 0], 0),
            34028708: ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 1),
            34015584: ([-1, 0, 0, 0, 1, 0, 0, 0, -1], [0, 0, 0], 1),
            3200: ([-1, 0, 0, 0, -1, 0, 0, 0, 1], [0, 0, 0], 0),
            34015424: ([-1, 0, 0, 0, -1, 0, 0, 0, 1], [0, 0, 0], 1),
            16320: ([1, 0, 0, 0, -1, 0, 0, 0, -1], [0, 0, 0], 0),
            34028544: ([1, 0, 0, 0, -1, 0, 0, 0, -1], [0, 0, 0], 1),
        }
        for raw, (rotation, translation, time_reversal) in expected.items():
            self.assertEqual(
                extractor._decode_magnetic_operation(raw),
                {"rotation": rotation, "translation_numerator": translation,
                 "time_reversal": time_reversal},
            )

        uni7 = [extractor._decode_magnetic_operation(raw) for raw in
                [16484, 3198, 34146806, 34133520]]
        self.assertEqual([item["rotation"] for item in uni7], [
            [1, 0, 0, 0, 1, 0, 0, 0, 1],
            [-1, 0, 0, 0, -1, 0, 0, 0, -1],
            [1, 0, 0, 0, 1, 0, 0, 0, 1],
            [-1, 0, 0, 0, -1, 0, 0, 0, -1],
        ])
        self.assertEqual([item["translation_numerator"] for item in uni7],
                         [[0, 0, 0], [0, 0, 0], [0, 0, 6], [0, 0, 6]])
        self.assertEqual([item["time_reversal"] for item in uni7], [0, 0, 1, 1])

        uni9_groups = [
            [16484, 3360, 34028708, 34015584],
            [16484, 3200, 34028708, 34015424],
            [16484, 16320, 34028708, 34028544],
        ]
        expected_groups = [
            [
                ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 0),
                ([-1, 0, 0, 0, 1, 0, 0, 0, -1], [0, 0, 0], 0),
                ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 1),
                ([-1, 0, 0, 0, 1, 0, 0, 0, -1], [0, 0, 0], 1),
            ],
            [
                ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 0),
                ([-1, 0, 0, 0, -1, 0, 0, 0, 1], [0, 0, 0], 0),
                ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 1),
                ([-1, 0, 0, 0, -1, 0, 0, 0, 1], [0, 0, 0], 1),
            ],
            [
                ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 0),
                ([1, 0, 0, 0, -1, 0, 0, 0, -1], [0, 0, 0], 0),
                ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0], 1),
                ([1, 0, 0, 0, -1, 0, 0, 0, -1], [0, 0, 0], 1),
            ],
        ]
        for raw_group, expected_group in zip(uni9_groups, expected_groups):
            decoded = [extractor._decode_magnetic_operation(raw) for raw in raw_group]
            self.assertEqual(
                [(item["rotation"], item["translation_numerator"], item["time_reversal"])
                 for item in decoded], expected_group,
            )

    def test_key_uni_operations_decode_antiunitary(self):
        msg = self.artifact["msg"]
        for uni in (7, 9):
            order, offset = msg["magnetic_spacegroup_operation_index"][uni][0]
            decoded = [extractor._decode_magnetic_operation(value)
                       for value in msg["magnetic_symmetry_operations"][offset:offset + order]]
            self.assertTrue(any(operation["time_reversal"] == 1 for operation in decoded), uni)
            self.assertTrue(all(len(operation["rotation"]) == 9 for operation in decoded))
            self.assertTrue(all(len(operation["translation_numerator"]) == 3 for operation in decoded))

    def test_canonical_artifact_and_manifest_hash(self):
        artifact_bytes = extractor.canonical_json(self.artifact)
        self.assertEqual(artifact_bytes[-1:], b"\n")
        self.assertNotIn(b".0", artifact_bytes)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "artifact.json"
            manifest = Path(directory) / "manifest.json"
            extractor.write_outputs(UPSTREAM, output, manifest)
            self.assertEqual(output.read_bytes(), ARTIFACT.read_bytes())
            value = json.loads(manifest.read_text())
            self.assertEqual(value["artifact"]["bytes"], output.stat().st_size)
            self.assertEqual(value["artifact"]["sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
            self.assertEqual(value["schema"], extractor.MANIFEST_SCHEMA)

    def test_two_runs_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.json"
            first_manifest = Path(directory) / "first.manifest.json"
            second = Path(directory) / "second.json"
            second_manifest = Path(directory) / "second.manifest.json"
            extractor.write_outputs(UPSTREAM, first, first_manifest)
            extractor.write_outputs(UPSTREAM, second, second_manifest)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(first_manifest.read_bytes().replace(b'first.json', b'second.json'), second_manifest.read_bytes())

    def test_initializer_grammar_rejects_missing_duplicate_trailing_and_bad_shape(self):
        with self.assertRaises(extractor.ExtractionError):
            extractor._initializer_text("int x[] = {1};", "missing")
        with self.assertRaises(extractor.ExtractionError):
            extractor._initializer_text("int x[] = {1}; int x[] = {2};", "x")
        with self.assertRaises(extractor.ExtractionError):
            extractor._parse_initializer("{1} trailing")
        with self.assertRaises(extractor.ExtractionError):
            extractor._normalize_rows([[1, 2, 3]], 1, 2, "fixture")

    def test_initializer_partial_rows_zero_fill_and_corrupt_tokens_reject(self):
        self.assertEqual(extractor._normalize_rows([[1]], 2, 2, "fixture"), [[1, 0], [0, 0]])
        self.assertEqual(extractor._normalize_3d([[[1]]], 2, 2, 2, "fixture"),
                         [[[1, 0], [0, 0]], [[0, 0], [0, 0]]])
        with self.assertRaises(extractor.ExtractionError):
            extractor._ints(extractor._parse_initializer("{1, nope}"), "fixture")
        with self.assertRaises(extractor.ExtractionError):
            extractor._decode_magnetic_operation(-1)
        with self.assertRaises(extractor.ExtractionError):
            extractor._decode_magnetic_operation(extractor.MSG_OPERATION_SCALE * 2)

    def test_manifest_source_hash_corruption_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "msg_database.c"
            path.write_bytes((UPSTREAM / "src/msg_database.c").read_bytes() + b" ")
            with self.assertRaises(extractor.ExtractionError):
                extractor._source(path, extractor.EXPECTED_SOURCES["msg_database.c"])


if __name__ == "__main__":
    unittest.main()
