#!/usr/bin/env python3
"""Focused tests for the canonical exact ISO-IR data--Hall sidecar."""

from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from . import freeze_iso_irrep_data_hall as freeze


DATA_DIR = Path(__file__).parent / "data"
ARTIFACT = DATA_DIR / "iso_irrep_data_hall_v1.json"
MANIFEST = DATA_DIR / "iso_irrep_data_hall_v1.manifest.json"
GOLDEN_ARTIFACT_BYTES = 697_730
GOLDEN_ARTIFACT_SHA256 = (
    "35bcb00958021eb6fc5a330f8dbf85a80be78ccec324f441e6138cdba4b617e0"
)
GOLDEN_MANIFEST_BYTES = 869
GOLDEN_MANIFEST_SHA256 = (
    "bc6aa7a94d698f2193e7cb623b16dded2dd8e0307d502cf28b68c554f364d7e2"
)


def _canonical(value):
    return freeze.canonical_json(value)


class FreezeDataHallTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact_bytes = ARTIFACT.read_bytes()
        cls.manifest_bytes = MANIFEST.read_bytes()
        cls.artifact, cls.manifest = freeze.parse_and_validate_pair(
            cls.artifact_bytes, cls.manifest_bytes
        )

    def test_committed_pair_golden_and_canonical(self):
        self.assertEqual(len(self.artifact_bytes), GOLDEN_ARTIFACT_BYTES)
        self.assertEqual(
            hashlib.sha256(self.artifact_bytes).hexdigest(),
            GOLDEN_ARTIFACT_SHA256,
        )
        self.assertEqual(len(self.manifest_bytes), GOLDEN_MANIFEST_BYTES)
        self.assertEqual(
            hashlib.sha256(self.manifest_bytes).hexdigest(),
            GOLDEN_MANIFEST_SHA256,
        )
        self.assertEqual(_canonical(self.artifact), self.artifact_bytes)
        self.assertEqual(_canonical(self.manifest), self.manifest_bytes)
        self.assertEqual(
            self.manifest["artifact"],
            {
                "path": "iso_irrep_data_hall_v1.json",
                "bytes": GOLDEN_ARTIFACT_BYTES,
                "sha256": GOLDEN_ARTIFACT_SHA256,
            },
        )

    def test_schema_census_and_fixed_witnesses(self):
        self.assertEqual(
            set(self.artifact),
            {
                "schema", "translation_denominator", "frame_semantics",
                "mapping_semantics", "inputs", "census", "spacegroups",
            },
        )
        self.assertEqual(self.artifact["schema"], freeze.SCHEMA)
        self.assertEqual(self.artifact["translation_denominator"], 12)
        self.assertEqual(
            self.artifact["frame_semantics"],
            "direct-source-frame-P=I-p=0",
        )
        self.assertEqual(
            self.artifact["mapping_semantics"],
            {
                "source_to_hall": "source=hall+shift",
                "hall_to_source": "hall=source+shift",
            },
        )
        census = self.artifact["census"]
        self.assertEqual(
            census,
            {
                **census,
                "pir_records": 10_294,
                "cir_records": 11_202,
                "source_representatives": 2_609,
                "raw_unique": 220,
                "raw_ambiguous": 10,
                "raw_missing": 0,
                "filtered_unique": 230,
                "filtered_ambiguous": 0,
                "filtered_missing": 0,
                "selected_hall_operations": 4_425,
                "source_to_hall": 2_609,
                "source_to_hall_nonzero": 0,
                "hall_to_source": 4_425,
                "hall_to_source_nonzero": 1_816,
                "expanded_normalization_nonzero": 410,
                "raw_ambiguous_spacegroups": [5, 8, 9, 12, 15, 21, 38, 39, 65, 67],
                "centering_counts": [
                    ["P", 149], ["A", 4], ["B", 0], ["C", 16],
                    ["F", 16], ["I", 38], ["R", 7],
                ],
            },
        )
        spacegroups = self.artifact["spacegroups"]
        self.assertEqual(spacegroups[0]["spacegroup"], 1)
        self.assertEqual(spacegroups[0]["raw_candidate_halls"], [1])
        self.assertEqual(spacegroups[0]["data_hall"], 1)
        sg5 = spacegroups[4]
        self.assertEqual(sg5["raw_candidate_halls"], [9, 10, 11])
        self.assertEqual(sg5["data_hall"], 9)
        self.assertEqual(sg5["centering"], "C")
        self.assertEqual(spacegroups[145]["centering"], "R")
        self.assertTrue(spacegroups[145]["source_symbol"].startswith("R"))
        self.assertEqual(spacegroups[224]["centering"], "F")
        self.assertTrue(spacegroups[224]["source_symbol"].startswith("F"))

    def test_mapping_and_distribution_semantics_are_rechecked(self):
        for index, record in enumerate(self.artifact["spacegroups"], 1):
            self.assertEqual(record["spacegroup"], index)
            source_maps = record["source_to_hall"]
            hall_maps = record["hall_to_source"]
            self.assertEqual(len(source_maps), record["source_operation_count"])
            self.assertEqual(len(hall_maps), record["hall_operation_count"])
            for source_index, mapping in enumerate(source_maps):
                self.assertEqual(mapping["source_operation_index"], source_index)
                self.assertEqual(
                    mapping["lattice_shift_numerator"], [0, 0, 0]
                )
                inverse = hall_maps[mapping["hall_operation_index"]]
                self.assertEqual(inverse["source_operation_index"], source_index)
                self.assertEqual(
                    inverse["lattice_shift_numerator"],
                    [-value for value in mapping["lattice_shift_numerator"]],
                )
            residues = set()
            for hall_index, mapping in enumerate(hall_maps):
                self.assertEqual(
                    mapping["hall_operation_index"], hall_index
                )
                source_index = mapping["source_operation_index"]
                residues.add(tuple(
                    value % freeze.TRANSLATION_DENOMINATOR
                    for value in mapping["lattice_shift_numerator"]
                ))
                self.assertGreaterEqual(source_index, 0)
            self.assertTrue(residues)

        # The production validator recomputes the full weighted distributions;
        # this explicit check protects the two distinct nonzero metrics.
        census = self.artifact["census"]
        self.assertEqual(
            sum(row[1] for row in census["hall_to_source_shifts"]),
            census["hall_to_source"],
        )
        self.assertEqual(
            sum(row[1] for row in census["expanded_normalization_shifts"]),
            census["hall_to_source"],
        )

    def test_fresh_build_matches_committed_pair(self):
        # One full production build is intentional; two additional cold
        # processes were used for the first-run byte-for-byte determinism gate.
        artifact = freeze.build_artifact()
        artifact_bytes = _canonical(artifact)
        manifest_bytes = _canonical(freeze.build_manifest(artifact_bytes))
        self.assertEqual(artifact_bytes, self.artifact_bytes)
        self.assertEqual(manifest_bytes, self.manifest_bytes)

    def test_canonical_json_rejects_non_schema_values(self):
        class StringSubclass(str):
            pass

        for value in (
            {"value": None}, {"value": True}, {"value": 1.0}, {"value": (1,)},
            {"value": StringSubclass("x")}, {"value": "é"},
        ):
            with self.subTest(value=value):
                with self.assertRaises(freeze.FreezeSchemaError):
                    _canonical(value)

        nested = 0
        for _ in range(2_000):
            nested = [nested]
        with self.assertRaises(freeze.FreezeSchemaError):
            _canonical(nested)

        nested_bytes = b"[" * 2_000 + b"0" + b"]" * 2_000 + b"\n"
        with self.assertRaises(freeze.FreezeSchemaError):
            freeze._parse_canonical_json(nested_bytes, "deep")

    def test_pair_parser_rejects_encoding_and_schema_corruption(self):
        corruptions = {
            "missing final LF": self.artifact_bytes[:-1],
            "extra final LF": self.artifact_bytes + b"\n",
            "CRLF": self.artifact_bytes.replace(b"\n", b"\r\n", 1),
            "nonascii": self.artifact_bytes[:-1] + "é\n".encode("utf-8"),
        }
        for name, broken in corruptions.items():
            with self.subTest(corruption=name):
                with self.assertRaises(freeze.FreezeError):
                    freeze.parse_and_validate_pair(broken, self.manifest_bytes)

        decoded = json.loads(self.artifact_bytes.decode("ascii"))
        decoded["unexpected"] = 1
        broken = _canonical(decoded)
        with self.assertRaises(freeze.FreezeIntegrityError):
            freeze.parse_and_validate_pair(broken, self.manifest_bytes)
        with self.assertRaises(freeze.FreezeSchemaError):
            freeze._parse_and_validate_uncommitted_pair(
                broken, self.manifest_bytes
            )

        decoded = json.loads(self.artifact_bytes.decode("ascii"))
        del decoded["census"]
        with self.assertRaises(freeze.FreezeIntegrityError):
            freeze.parse_and_validate_pair(_canonical(decoded), self.manifest_bytes)
        with self.assertRaises(freeze.FreezeSchemaError):
            freeze._parse_and_validate_uncommitted_pair(
                _canonical(decoded), self.manifest_bytes
            )

        duplicate = b'{"schema":"x","schema":"y"}\n'
        with self.assertRaises(freeze.FreezeIntegrityError):
            freeze.parse_and_validate_pair(duplicate, self.manifest_bytes)
        with self.assertRaises(freeze.FreezeSchemaError):
            freeze._parse_and_validate_uncommitted_pair(
                duplicate, self.manifest_bytes
            )

    def test_pair_parser_rejects_synchronized_semantic_mutations(self):
        def resigned(mutator):
            decoded = deepcopy(self.artifact)
            mutator(decoded)
            broken_artifact = _canonical(decoded)
            broken_manifest = deepcopy(self.manifest)
            broken_manifest["artifact"]["bytes"] = len(broken_artifact)
            broken_manifest["artifact"]["sha256"] = hashlib.sha256(
                broken_artifact
            ).hexdigest()
            return broken_artifact, _canonical(broken_manifest)

        mutations = {
            "SG5 Hall 9 to 10": lambda value: value["spacegroups"][4].__setitem__(
                "data_hall", 10
            ),
            "raw candidate replacement": lambda value: value["spacegroups"][4].__setitem__(
                "raw_candidate_halls", [9, 10, 12]
            ),
            "anchor increment": lambda value: value["spacegroups"][4].__setitem__(
                "pir_anchor_irnumber",
                value["spacegroups"][4]["pir_anchor_irnumber"] + 1,
            ),
            "fake C symbol": lambda value: value["spacegroups"][4].__setitem__(
                "source_symbol", "C FAKE"
            ),
        }
        for name, mutator in mutations.items():
            with self.subTest(mutation=name):
                broken_artifact, broken_manifest = resigned(mutator)
                with self.assertRaises(freeze.FreezeIntegrityError):
                    freeze.parse_and_validate_pair(broken_artifact, broken_manifest)

        def permute_p_mapping(value):
            record = value["spacegroups"][1]
            record["source_to_hall"][0]["hall_operation_index"], record[
                "source_to_hall"
            ][1]["hall_operation_index"] = (
                record["source_to_hall"][1]["hall_operation_index"],
                record["source_to_hall"][0]["hall_operation_index"],
            )
            record["hall_to_source"][0]["source_operation_index"], record[
                "hall_to_source"
            ][1]["source_operation_index"] = (
                record["hall_to_source"][1]["source_operation_index"],
                record["hall_to_source"][0]["source_operation_index"],
            )

        broken_artifact, broken_manifest = resigned(permute_p_mapping)
        with self.assertRaises(freeze.FreezeIntegrityError):
            freeze.parse_and_validate_pair(broken_artifact, broken_manifest)

        def shift_hall_mapping(value):
            record = value["spacegroups"][4]
            record["hall_to_source"][0]["lattice_shift_numerator"][0] += 12
            record["source_to_hall"][0]["lattice_shift_numerator"][0] -= 12
            record["census"] = freeze._aggregate_census(
                value["spacegroups"],
                value["census"]["pir_records"],
                value["census"]["cir_records"],
            )

        broken_artifact, broken_manifest = resigned(shift_hall_mapping)
        with self.assertRaises(freeze.FreezeIntegrityError):
            freeze.parse_and_validate_pair(broken_artifact, broken_manifest)

        decoded = deepcopy(self.artifact)
        decoded["inputs"]["pir_zip"]["sha256"] = "0" * 64
        broken_artifact = _canonical(decoded)
        broken_manifest = deepcopy(self.manifest)
        broken_manifest["inputs"]["pir_zip"]["sha256"] = "0" * 64
        broken_manifest["artifact"]["bytes"] = len(broken_artifact)
        broken_manifest["artifact"]["sha256"] = hashlib.sha256(
            broken_artifact
        ).hexdigest()
        with self.assertRaises(freeze.FreezeIntegrityError):
            freeze.parse_and_validate_pair(
                broken_artifact, _canonical(broken_manifest)
            )

    def test_manifest_builder_closes_pair(self):
        manifest = freeze.build_manifest(self.artifact_bytes)
        self.assertEqual(_canonical(manifest), self.manifest_bytes)
        with self.assertRaises(freeze.FreezeSchemaError):
            freeze.build_manifest(self.artifact_bytes[:-1])

    def test_atomic_write_paths_and_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "artifact.json"
            manifest = root / "manifest.json"
            with mock.patch.object(freeze, "build_artifact", return_value=self.artifact):
                written_artifact, written_manifest = freeze.write_outputs(
                    output, manifest
                )
            self.assertEqual(output.read_bytes(), written_artifact)
            self.assertEqual(manifest.read_bytes(), written_manifest)
            freeze.parse_and_validate_pair(written_artifact, written_manifest)
            self.assertEqual(list(root.glob(".*.tmp")), [])

            with self.assertRaises(freeze.FreezeInvariantError):
                freeze.write_outputs(output, output)

            hardlink = root / "hardlink.json"
            os.link(output, hardlink)
            with self.assertRaises(freeze.FreezeInvariantError):
                freeze.write_outputs(output, hardlink)

    def test_replacement_identity_rejects_aliases_before_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            real = root / "real"
            real.mkdir()
            alias = root / "alias"
            os.symlink(real, alias, target_is_directory=True)
            output = real / "same.json"
            manifest = alias / "same.json"
            with mock.patch.object(
                freeze, "build_artifact", side_effect=AssertionError("built")
            ) as builder:
                with self.assertRaises(freeze.FreezeInvariantError):
                    freeze.write_outputs(output, manifest)
            self.assertEqual(builder.call_count, 0)

            sub = root / "sub"
            sub.mkdir()
            output = sub / ".." / "same2.json"
            manifest = root / "same2.json"
            with mock.patch.object(
                freeze, "build_artifact", side_effect=AssertionError("built")
            ) as builder:
                with self.assertRaises(freeze.FreezeInvariantError):
                    freeze.write_outputs(output, manifest)
            self.assertEqual(builder.call_count, 0)

    def test_replacement_target_is_bound_before_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first"
            second = root / "second"
            first.mkdir()
            second.mkdir()
            alias = root / "alias"
            os.symlink(second, alias, target_is_directory=True)
            output = first / "same.json"
            manifest = alias / "same.json"

            def retarget_then_build():
                alias.unlink()
                os.symlink(first, alias, target_is_directory=True)
                return self.artifact

            with mock.patch.object(
                freeze, "build_artifact", side_effect=retarget_then_build
            ):
                written_artifact, written_manifest = freeze.write_outputs(
                    output, manifest
                )
            self.assertEqual(output.read_bytes(), written_artifact)
            self.assertEqual((second / "same.json").read_bytes(), written_manifest)
            self.assertEqual((alias / "same.json").read_bytes(), written_artifact)

    def test_native_path_string_errors_precede_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            invalid_paths = (
                root / "bad\x00.json",
                root / "bad\x00parent" / "artifact.json",
                root / "\ud800.json",
            )
            for invalid in invalid_paths:
                with self.subTest(path=repr(str(invalid))):
                    with mock.patch.object(
                        freeze,
                        "build_artifact",
                        side_effect=AssertionError("built"),
                    ) as builder:
                        with self.assertRaises(freeze.FreezeSchemaError):
                            freeze.write_outputs(invalid, root / "manifest.json")
                    self.assertEqual(builder.call_count, 0)

    def test_pathlike_failures_are_schema_errors_before_build(self):
        class ExplodingPath:
            def __init__(self, error):
                self.error = error

            def __fspath__(self):
                raise self.error

        errors = (
            OSError("synthetic path failure"),
            RuntimeError("synthetic path failure"),
            ValueError("synthetic path failure"),
            UnicodeEncodeError("ascii", "\ud800", 0, 1, "synthetic path failure"),
        )
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "manifest.json"
            for error in errors:
                with self.subTest(error=type(error).__name__):
                    with mock.patch.object(
                        freeze,
                        "build_artifact",
                        side_effect=AssertionError("built"),
                    ) as builder:
                        with self.assertRaises(freeze.FreezeSchemaError):
                            freeze.write_outputs(ExplodingPath(error), manifest)
                    self.assertEqual(builder.call_count, 0)

    def test_atomic_write_failure_cleans_staged_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "artifact.json"
            manifest = root / "manifest.json"
            real_replace = freeze.os.replace
            calls = []

            def replace_then_fail(source, target):
                calls.append((source, target))
                if len(calls) == 1:
                    return real_replace(source, target)
                raise OSError("synthetic replace failure")

            with mock.patch.object(freeze, "build_artifact", return_value=self.artifact):
                with mock.patch.object(
                    freeze.os,
                    "replace",
                    side_effect=replace_then_fail,
                ):
                    with self.assertRaises(freeze.FreezeIntegrityError):
                        freeze.write_outputs(output, manifest)
            self.assertEqual(len(calls), 2)
            self.assertTrue(output.exists())
            self.assertEqual(output.read_bytes(), self.artifact_bytes)
            self.assertFalse(manifest.exists())
            self.assertEqual(list(root.glob(".*.tmp")), [])

    def test_source_has_no_forbidden_runtime_fallbacks(self):
        source = Path(freeze.__file__).read_text(encoding="utf-8")
        for forbidden in (
            "generated_data", "SG_DATA_HALL", "hall_operations.json",
        ):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
