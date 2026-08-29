#!/usr/bin/env python3
"""Focused tests for the pinned spglib magnetic provenance extractor."""

import copy
import hashlib
import os
import subprocess
import sys
import tempfile
from unittest import mock
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import extract_spglib_magnetic_provenance as extractor


UPSTREAM = (Path(os.environ["SPGLIB_V2_5_0_SOURCE"])
            if os.environ.get("SPGLIB_V2_5_0_SOURCE") else None)
ARTIFACT = Path(__file__).parent / "data/spglib_magnetic_provenance_v1.json"
MANIFEST = Path(__file__).parent / "data/spglib_magnetic_provenance_v1.manifest.json"
GOLDEN_ARTIFACT_BYTES = 1_537_875
GOLDEN_ARTIFACT_SHA256 = "933a52a6696e7f6a1a2e426825ad92c377c6e96330e18c5c045d659798d740b9"
REAL_DECLARATION_SIGNATURES = {
    "spacegroup_types": "static SpacegroupType const spacegroup_types[] =",
    "symmetry_operation_index": "static int const symmetry_operation_index[][2] =",
    "symmetry_operations": "static int const symmetry_operations[] =",
    "magnetic_spacegroup_types": "static const MagneticSpacegroupType magnetic_spacegroup_types[] =",
    "magnetic_spacegroup_uni_mapping": "static const int magnetic_spacegroup_uni_mapping[][2] =",
    "magnetic_spacegroup_hall_mapping": "static const int magnetic_spacegroup_hall_mapping[][2] =",
    "magnetic_spacegroup_operation_index": "static const int magnetic_spacegroup_operation_index[][18][2] =",
    "magnetic_symmetry_operations": "static const int magnetic_symmetry_operations[] =",
    "alternative_transformations": "static const int alternative_transformations[][18][7] =",
}


class MagneticProvenanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if UPSTREAM is None:
            cls.artifact = None
            cls.details = None
        else:
            cls.artifact, cls.details = extractor.extract(UPSTREAM)

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for regeneration tests")
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

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for regeneration tests")
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

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for regeneration tests")
    def test_key_uni_operations_decode_antiunitary(self):
        msg = self.artifact["msg"]
        for uni in (7, 9):
            order, offset = msg["magnetic_spacegroup_operation_index"][uni][0]
            decoded = [extractor._decode_magnetic_operation(value)
                       for value in msg["magnetic_symmetry_operations"][offset:offset + order]]
            self.assertTrue(any(operation["time_reversal"] == 1 for operation in decoded), uni)
            self.assertTrue(all(len(operation["rotation"]) == 9 for operation in decoded))
            self.assertTrue(all(len(operation["translation_numerator"]) == 3 for operation in decoded))

    def test_committed_artifact_and_manifest_integrity(self):
        artifact_bytes = ARTIFACT.read_bytes()
        self.assertEqual(len(artifact_bytes), GOLDEN_ARTIFACT_BYTES)
        self.assertEqual(hashlib.sha256(artifact_bytes).hexdigest(), GOLDEN_ARTIFACT_SHA256)
        artifact = extractor._parse_json_bytes(artifact_bytes, str(ARTIFACT))
        extractor.validate_artifact(artifact)
        self.assertEqual(extractor.canonical_json(artifact), artifact_bytes)
        manifest_bytes = MANIFEST.read_bytes()
        manifest = extractor._parse_json_bytes(manifest_bytes, str(MANIFEST))
        extractor.validate_manifest(manifest, artifact_bytes, ARTIFACT.name)
        self.assertEqual(extractor.canonical_json(manifest), manifest_bytes)

    def test_thin_committed_pair_parser_closes_canonical_bytes(self):
        artifact_bytes = ARTIFACT.read_bytes()
        manifest_bytes = MANIFEST.read_bytes()
        artifact = extractor.parse_and_validate_committed_pair(
            artifact_bytes, manifest_bytes, ARTIFACT.name
        )
        self.assertEqual(extractor.canonical_json(artifact), artifact_bytes)
        with self.assertRaises(extractor.ExtractionError):
            extractor.parse_and_validate_committed_pair(
                b"{}", manifest_bytes, ARTIFACT.name
            )
        with self.assertRaises(extractor.ExtractionError):
            extractor.parse_and_validate_committed_pair(
                artifact_bytes, manifest_bytes[:-1], ARTIFACT.name
            )

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for regeneration tests")
    def test_canonical_artifact_and_manifest_hash(self):
        artifact_bytes = extractor.canonical_json(self.artifact)
        self.assertEqual(artifact_bytes[-1:], b"\n")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "artifact.json"
            manifest = Path(directory) / "manifest.json"
            extractor.write_outputs(UPSTREAM, output, manifest)
            self.assertEqual(output.read_bytes(), ARTIFACT.read_bytes())
            value = extractor._load_json(manifest)
            extractor.validate_manifest(value, output.read_bytes(), output.name)

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for regeneration tests")
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
            extractor._initializer_text("int x[] = junk {1};", "x")
        with self.assertRaises(extractor.ExtractionError):
            extractor._initializer_text("int x[] = {1} junk;", "x")
        with self.assertRaises(extractor.ExtractionError):
            extractor._initializer_text("int x[] = {1}", "x")
        with self.assertRaises(extractor.ExtractionError):
            extractor._initializer_text("x = {1};", "x")
        with self.assertRaises(extractor.ExtractionError):
            extractor._initializer_text("x[] = {1};", "x")
        with self.assertRaises(extractor.ExtractionError):
            extractor._normalize_rows([[1, 2, 3]], 1, 2, "fixture")

    def test_initializer_partial_rows_zero_fill_and_corrupt_tokens_reject(self):
        self.assertEqual(extractor._normalize_rows([[1]], 2, 2, "fixture"),
                         [[1, 0], [0, 0]])
        self.assertEqual(extractor._normalize_3d([[[1]], []], 2, 2, 2, "fixture"),
                         [[[1, 0], [0, 0]], [[0, 0], [0, 0]]])
        with self.assertRaises(extractor.ExtractionError):
            extractor._normalize_3d([[[1]]], 2, 2, 2, "fixture")
        with self.assertRaises(extractor.ExtractionError):
            extractor._ints(extractor._parse_initializer("{1, nope}"), "fixture")
        with self.assertRaises(extractor.ExtractionError):
            extractor._decode_magnetic_operation(-1)
        with self.assertRaises(extractor.ExtractionError):
            extractor._decode_magnetic_operation(extractor.MSG_OPERATION_SCALE * 2)

    def test_artifact_validator_rejects_boolean_and_missing_outer_row(self):
        artifact = extractor._load_json(ARTIFACT)
        tampered = copy.deepcopy(artifact)
        tampered["msg"]["magnetic_symmetry_operations"][1] = True
        with self.assertRaises(extractor.ExtractionError):
            extractor.validate_artifact(tampered)
        tampered = copy.deepcopy(artifact)
        tampered["msg"]["magnetic_spacegroup_operation_index"] = (
            tampered["msg"]["magnetic_spacegroup_operation_index"][:-1]
        )
        with self.assertRaises(extractor.ExtractionError):
            extractor.validate_artifact(tampered)

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for mapping tests")
    def test_hall_mapping_and_parent_spacegroup_invariants(self):
        msg_types = [
            [row["uni"], row["litvin"], row["bns"], row["og"],
             row["parent_spacegroup"], row["type"]]
            for row in self.artifact["msg"]["magnetic_spacegroup_types"]
        ]
        hall_mapping = copy.deepcopy(self.details["magnetic_spacegroup_hall_mapping"])
        extractor._validate_hall_mapping(
            hall_mapping,
            self.artifact["msg"]["magnetic_spacegroup_uni_mapping"],
            msg_types,
            self.artifact["spg"]["spacegroup_number"],
        )
        self.assertEqual(self.artifact["msg"]["magnetic_spacegroup_uni_mapping"][1], [1, 1])
        tampered = copy.deepcopy(self.artifact["msg"]["magnetic_spacegroup_uni_mapping"])
        tampered[1] = [1, 2]
        with self.assertRaises(extractor.ExtractionError):
            extractor._validate_hall_mapping(
                hall_mapping, tampered, msg_types,
                self.artifact["spg"]["spacegroup_number"],
            )

    def test_strict_c_comments_numbers_and_strings(self):
        self.assertEqual(extractor._strip_comments('"http://example.invalid"'),
                         '"http://example.invalid"')
        self.assertEqual(extractor._parse_initializer('{"http://example.invalid"}'),
                         ["http://example.invalid"])
        self.assertEqual(extractor._parse_initializer(r'{"\x41\101\n"}'),
                         ["AA\n"])
        with self.assertRaises(extractor.ExtractionError):
            extractor._parse_initializer("{1/**/2}")
        with self.assertRaises(extractor.ExtractionError):
            extractor._parse_initializer("{1\\\n,2}")
        parsed = extractor._parse_initializer("{\"quoted\", bare_identifier}")
        self.assertEqual(parsed[0], "quoted")
        self.assertIs(type(parsed[0]), extractor.CString)
        self.assertEqual(parsed[1], "bare_identifier")
        self.assertIs(type(parsed[1]), extractor.Identifier)
        for spelling in ("09", "012", "+1", "1.5", "1e3"):
            with self.subTest(spelling=spelling):
                with self.assertRaises(extractor.ExtractionError):
                    extractor._parse_initializer("{" + spelling + "}")
        for spelling in (r'"\u002f"', r'"\/"', r'"\x110000"', r'"\xD800"'):
            with self.subTest(spelling=spelling):
                with self.assertRaises(extractor.ExtractionError):
                    extractor._parse_initializer("{" + spelling + "}")

    def test_json_no_float_is_recursive_and_strict(self):
        for text in (
            b'{"x":1.5}', b'{"x":1e3}', b'{"x":NaN}',
            b'{"x":Infinity}', b'{"x":-Infinity}',
            b'{"x":[{"y":1.5}]}',
        ):
            with self.subTest(text=text):
                with self.assertRaises(extractor.ExtractionError):
                    extractor._parse_json_bytes(text, "fixture")
        for value in ({"x": 1.5}, {"x": float("nan")}, {"x": float("inf")}):
            with self.subTest(value=value):
                with self.assertRaises(extractor.ExtractionError):
                    extractor.canonical_json(value)
        with self.assertRaises(extractor.ExtractionError):
            extractor.canonical_json({"x": "\ud800"})

    def test_committed_manifest_corruption_fails_closed(self):
        artifact_bytes = ARTIFACT.read_bytes()
        manifest_bytes = MANIFEST.read_bytes()
        broken_bytes = manifest_bytes.replace(
            extractor.EXPECTED_UPSTREAM_COMMIT.encode("ascii"), b"0" * 40
        )
        self.assertNotEqual(broken_bytes, manifest_bytes)
        with tempfile.TemporaryDirectory() as directory:
            broken = Path(directory) / "broken.json"
            broken.write_bytes(broken_bytes)
            broken_manifest = extractor._load_json(broken)
            with self.assertRaises(extractor.ExtractionError):
                extractor.validate_manifest(broken_manifest, artifact_bytes, ARTIFACT.name)

    def test_upstream_provenance_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(extractor.ExtractionError):
                extractor._verify_upstream_provenance(Path(directory))
        with mock.patch.object(
            extractor, "_git_output", side_effect=["true", "wrong-commit"]
        ):
            with self.assertRaises(extractor.ExtractionError):
                extractor._verify_upstream_provenance(Path("/pinned/source"))
        with mock.patch.object(
            extractor, "_git_output",
            side_effect=["true", extractor.EXPECTED_UPSTREAM_COMMIT, "other-tag"],
        ):
            with self.assertRaises(extractor.ExtractionError):
                extractor._verify_upstream_provenance(Path("/pinned/source"))

    def test_output_targets_reject_same_path_and_inode_aliases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            same = root / "same.json"
            with self.assertRaises(extractor.ExtractionError):
                extractor._validate_output_targets(same, same)
            target = root / "target.json"
            target.write_bytes(b"artifact")
            hardlink = root / "hardlink.json"
            os.link(target, hardlink)
            with self.assertRaises(extractor.ExtractionError):
                extractor._validate_output_targets(target, hardlink)
            with self.assertRaises(extractor.ExtractionError):
                extractor.write_outputs(root / "not-a-checkout", target, hardlink)
            symlink = root / "symlink.json"
            os.symlink(target, symlink)
            with self.assertRaises(extractor.ExtractionError):
                extractor._validate_output_targets(symlink, hardlink)
            with self.assertRaises(extractor.ExtractionError):
                extractor.write_outputs(root / "not-a-checkout", symlink, hardlink)
            with self.assertRaises(extractor.ExtractionError):
                extractor.write_outputs(root / "not-a-checkout", same, same)

    def test_atomic_write_wraps_parent_errors_and_cleans_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            parent_file = root / "parent"
            parent_file.write_bytes(b"old parent")
            with self.assertRaises(extractor.ExtractionError):
                extractor._atomic_write(parent_file / "artifact.json", b"new")
            self.assertEqual(parent_file.read_bytes(), b"old parent")
            self.assertEqual(list(root.iterdir()), [parent_file])

            target = root / "existing.json"
            target.write_bytes(b"old target")
            with mock.patch.object(extractor.os, "replace", side_effect=OSError("injected")):
                with self.assertRaises(extractor.ExtractionError):
                    extractor._atomic_write(target, b"new target")
            self.assertEqual(target.read_bytes(), b"old target")
            self.assertEqual(sorted(path.name for path in root.iterdir()),
                             ["existing.json", "parent"])

    def test_real_pinned_declaration_signatures_are_strict(self):
        self.assertEqual(set(REAL_DECLARATION_SIGNATURES),
                         set(extractor.PINNED_DECLARATIONS))
        for name, signature in REAL_DECLARATION_SIGNATURES.items():
            with self.subTest(name=name):
                valid = signature + " {1};"
                self.assertEqual(extractor._initializer_text(valid, name), "{1}")

                missing_type = valid[valid.index(name):]
                with self.assertRaises(extractor.ExtractionError):
                    extractor._initializer_text(missing_type, name)

                wrong_dimension = signature.replace("[]", "[1]", 1) + " {1};"
                with self.assertRaises(extractor.ExtractionError):
                    extractor._initializer_text(wrong_dimension, name)

                trailing = signature + " {1} junk;"
                with self.assertRaises(extractor.ExtractionError):
                    extractor._initializer_text(trailing, name)

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for checkout tests")
    def test_pinned_blob_ignores_dirty_and_crlf_worktree(self):
        expected_artifact = ARTIFACT.read_bytes()
        with tempfile.TemporaryDirectory() as directory:
            checkout = Path(directory) / "checkout"
            subprocess.run(
                ["git", "clone", "--no-hardlinks", "--local",
                 str(UPSTREAM), str(checkout)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            source_paths = [checkout / "src/msg_database.c",
                            checkout / "src/spg_database.c"]
            for path in source_paths:
                path.write_bytes(path.read_bytes() + b"\n/* dirty worktree */\n")
            self.assertTrue(extractor._git_output(checkout, ["status", "--short"]))
            artifact, _ = extractor.extract(checkout)
            self.assertEqual(extractor.canonical_json(artifact), expected_artifact)

            for path in source_paths:
                path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
                self.assertIn(b"\r\n", path.read_bytes())
            artifact, _ = extractor.extract(checkout)
            self.assertEqual(extractor.canonical_json(artifact), expected_artifact)

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for blob tests")
    def test_git_source_uses_pinned_blob_without_working_tree_read(self):
        path = UPSTREAM / "src/msg_database.c"
        source_bytes = path.read_bytes()
        expected_hash = hashlib.sha256(source_bytes).hexdigest()
        with mock.patch.object(extractor, "_git_blob", return_value=source_bytes) as git_blob:
            with mock.patch.object(Path, "read_bytes", side_effect=AssertionError("working tree read")):
                source, actual_hash = extractor._git_source(
                    UPSTREAM, "src/msg_database.c", expected_hash
                )
        self.assertTrue(source)
        self.assertEqual(actual_hash, expected_hash)
        git_blob.assert_called_once_with(UPSTREAM, "src/msg_database.c")

    def test_git_blob_reports_binary_subprocess_errors(self):
        error = subprocess.CalledProcessError(
            1, ["git"], stderr=b"fatal: missing blob\n"
        )
        with mock.patch.object(extractor.subprocess, "run", side_effect=error):
            with self.assertRaises(extractor.ExtractionError) as raised:
                extractor._git_blob(Path("/pinned/source"), "src/msg_database.c")
        self.assertIn("missing blob", str(raised.exception))

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for source tests")
    def test_source_bytes_are_read_once_and_strictly_decoded(self):
        path = UPSTREAM / "src/msg_database.c"
        source_bytes = path.read_bytes()
        expected_hash = hashlib.sha256(source_bytes).hexdigest()
        with mock.patch.object(Path, "read_text", side_effect=AssertionError("TOCTOU read")):
            with mock.patch.object(Path, "read_bytes", wraps=path.read_bytes) as read_bytes:
                source, actual_hash = extractor._source(path, expected_hash)
        self.assertEqual(actual_hash, expected_hash)
        self.assertTrue(source)
        self.assertEqual(read_bytes.call_count, 1)
        with tempfile.TemporaryDirectory() as directory:
            invalid = Path(directory) / "invalid.c"
            invalid_bytes = b"\xff"
            invalid.write_bytes(invalid_bytes)
            with self.assertRaises(extractor.ExtractionError):
                extractor._source(invalid, hashlib.sha256(invalid_bytes).hexdigest())

    @unittest.skipUnless(UPSTREAM is not None, "set SPGLIB_V2_5_0_SOURCE for source tests")
    def test_manifest_source_hash_corruption_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "msg_database.c"
            path.write_bytes((UPSTREAM / "src/msg_database.c").read_bytes() + b" ")
            with self.assertRaises(extractor.ExtractionError):
                extractor._source(path, extractor.EXPECTED_SOURCES["msg_database.c"])


if __name__ == "__main__":
    unittest.main()
