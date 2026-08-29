#!/usr/bin/env python3
"""Focused tests for the independent frozen data--Hall provenance loader."""

import ast
from dataclasses import FrozenInstanceError, fields, is_dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import tempfile
import threading
import unittest
from unittest import mock

from . import iso_irrep_data_hall as loader


DATA_DIR = Path(__file__).parent / "data"
ARTIFACT = DATA_DIR / "iso_irrep_data_hall_v1.json"
MANIFEST = DATA_DIR / "iso_irrep_data_hall_v1.manifest.json"


def _canonical(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii") + b"\n"


def _read_pair():
    return ARTIFACT.read_bytes(), MANIFEST.read_bytes()


class DataHallLoaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact_bytes, cls.manifest_bytes = _read_pair()
        cls.database = loader.load_committed_data_hall_provenance()

    def tearDown(self):
        loader._reset_cache_for_test()

    def test_fixed_pair_census_and_witnesses(self):
        self.assertEqual(len(self.artifact_bytes), loader.ARTIFACT_BYTE_LENGTH)
        self.assertEqual(
            hashlib.sha256(self.artifact_bytes).hexdigest(), loader.ARTIFACT_SHA256
        )
        self.assertEqual(len(self.manifest_bytes), loader.MANIFEST_BYTE_LENGTH)
        self.assertEqual(
            hashlib.sha256(self.manifest_bytes).hexdigest(), loader.MANIFEST_SHA256
        )
        database = self.database
        self.assertIs(type(database), loader.DataHallProvenanceDatabase)
        self.assertIs(type(database.frames), tuple)
        self.assertEqual(len(database.frames), 230)
        self.assertEqual(
            tuple(field.name for field in fields(loader.DataHallCensus)),
            (
                "pir_records", "cir_records", "source_representatives",
                "raw_unique", "raw_ambiguous", "raw_missing",
                "raw_ambiguous_spacegroups", "filtered_unique",
                "filtered_ambiguous", "filtered_missing",
                "selected_hall_operations", "source_to_hall",
                "source_to_hall_nonzero", "hall_to_source",
                "hall_to_source_nonzero", "hall_to_source_shifts",
                "hall_to_source_cosets", "expanded_normalization_nonzero",
                "expanded_normalization_shifts", "centering_counts",
            ),
        )
        census = database.census
        self.assertEqual(census.pir_records, 10_294)
        self.assertEqual(census.cir_records, 11_202)
        self.assertEqual(census.source_representatives, 2_609)
        self.assertEqual(census.selected_hall_operations, 4_425)
        self.assertEqual(
            (census.raw_unique, census.raw_ambiguous, census.raw_missing),
            (220, 10, 0),
        )
        self.assertEqual(
            census.raw_ambiguous_spacegroups,
            (5, 8, 9, 12, 15, 21, 38, 39, 65, 67),
        )
        self.assertEqual(
            (census.filtered_unique, census.filtered_ambiguous, census.filtered_missing),
            (230, 0, 0),
        )
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
        for spacegroup in (1, 3, 5, 146, 225):
            frame = database.frames[spacegroup - 1]
            self.assertEqual(frame.spacegroup, spacegroup)
            self.assertEqual(loader.data_hall_frame(spacegroup), frame)
        self.assertEqual(database.frames[0].data_hall, 1)
        self.assertEqual(database.frames[2].data_hall, 3)
        self.assertEqual(database.frames[4].raw_candidate_halls, (9, 10, 11))
        self.assertEqual(database.frames[4].data_hall, 9)
        self.assertEqual(database.frames[4].hall_to_source[2].shift_numerator, (6, 6, 0))
        self.assertEqual(database.frames[145].centering, "R")
        self.assertEqual(database.frames[224].centering, "F")

    def test_public_apis_identity_and_exact_lookup(self):
        first = loader.load_committed_data_hall_provenance()
        second = loader.load_committed_data_hall_provenance()
        self.assertIs(first, second)
        self.assertIs(loader.data_hall_frame(5), first.frames[4])
        self.assertEqual(loader.data_hall_for_spacegroup(5), 9)
        class IntSubclass(int):
            pass

        for value in (True, 0, 231, 1.0, IntSubclass(1)):
            with self.subTest(value=value):
                with self.assertRaises(loader.DataHallLookupError):
                    loader.data_hall_frame(value)

    def test_graph_is_frozen_and_has_no_mutable_payload(self):
        seen = set()

        def visit(value):
            if id(value) in seen:
                return
            seen.add(id(value))
            if isinstance(value, (list, dict)):
                self.fail("mutable payload leaked into public graph")
            if is_dataclass(value) and not isinstance(value, type):
                for field in fields(value):
                    visit(getattr(value, field.name))
            elif type(value) is tuple:
                for child in value:
                    visit(child)

        visit(self.database)
        with self.assertRaises(FrozenInstanceError):
            self.database.census = self.database.census
        with self.assertRaises(FrozenInstanceError):
            self.database.frames[0].data_hall = 1
        with self.assertRaises(loader.DataHallSchemaError):
            loader.SourceToHallMapping(0, 0, [0, 0, 0])
        with self.assertRaises(loader.DataHallSchemaError):
            replace(self.database.frames[0], basis=[1] * 9)
        with self.assertRaises(loader.DataHallSchemaError):
            replace(self.database.census, centering_counts=[])
        with self.assertRaises(loader.DataHallSchemaError):
            loader.DataHallProvenanceDatabase(list(self.database.frames), self.database.census)

    def test_private_parser_rejects_encoding_and_semantic_corruption(self):
        for data in (
            self.artifact_bytes[:-1],
            self.artifact_bytes + b"\n",
            self.artifact_bytes.replace(b"\n", b"\r\n", 1),
            self.artifact_bytes[:-1] + "é\n".encode("utf-8"),
        ):
            with self.assertRaises(loader.DataHallProvenanceError):
                loader._parse_canonical_json(data, "synthetic")
        for data in (
            b"{\"x\":null}\n",
            b"{\"x\":true}\n",
            b"{\"x\":1.0}\n",
            b"{\"x\":NaN}\n",
            b"{\"x\":1,\"x\":2}\n",
        ):
            with self.assertRaises(loader.DataHallSchemaError):
                loader._parse_canonical_json(data, "synthetic")
        deep = b"[" * 2_000 + b"0" + b"]" * 2_000 + b"\n"
        with self.assertRaises(loader.DataHallSchemaError):
            loader._parse_canonical_json(deep, "deep")

        artifact = json.loads(self.artifact_bytes.decode("ascii"))
        artifact["spacegroups"][4]["source_to_hall"][0]["hall_operation_index"] = 99
        with self.assertRaises(loader.DataHallProvenanceError):
            loader._parse_and_validate_pair(_canonical(artifact), self.manifest_bytes)
        artifact = json.loads(self.artifact_bytes.decode("ascii"))
        artifact["census"]["source_to_hall_nonzero"] = 1
        with self.assertRaises(loader.DataHallProvenanceError):
            loader._parse_and_validate_pair(_canonical(artifact), self.manifest_bytes)

    def test_public_fixed_gate_rejects_corrupt_files(self):
        real_read_bytes = loader.Path.read_bytes

        def fake_read(path):
            if path == loader._ARTIFACT_PATH:
                return self.artifact_bytes[:-1]
            return real_read_bytes(path)

        loader._reset_cache_for_test()
        with mock.patch.object(loader.Path, "read_bytes", autospec=True, side_effect=fake_read):
            with self.assertRaises(loader.DataHallIntegrityError):
                loader.load_committed_data_hall_provenance()
        self.assertIsNone(loader._DATABASE)

        def swapped_read(path):
            if path == loader._ARTIFACT_PATH:
                return self.manifest_bytes
            if path == loader._MANIFEST_PATH:
                return self.artifact_bytes
            return real_read_bytes(path)

        with mock.patch.object(loader.Path, "read_bytes", autospec=True, side_effect=swapped_read):
            with self.assertRaises(loader.DataHallIntegrityError):
                loader.load_committed_data_hall_provenance()
        self.assertIsNone(loader._DATABASE)

    def test_chdir_does_not_change_anchored_paths(self):
        loader._reset_cache_for_test()
        original = Path.cwd()
        try:
            with tempfile.TemporaryDirectory() as directory:
                os.chdir(directory)
                database = loader.load_committed_data_hall_provenance()
        finally:
            os.chdir(original)
        self.assertEqual(database.frames[4].data_hall, 9)

    def test_single_flight_reads_once_with_forced_overlap(self):
        class TrackingLock:
            def __init__(self):
                self.lock = threading.Lock()
                self.held = False
                self.second_attempt = threading.Event()

            def __enter__(self):
                if self.held:
                    self.second_attempt.set()
                self.lock.acquire()
                self.held = True
                return self

            def __exit__(self, exc_type, exc, traceback):
                self.held = False
                self.lock.release()

        loader._reset_cache_for_test()
        tracking_lock = TrackingLock()
        first_read_started = threading.Event()
        counts = {loader._ARTIFACT_PATH: 0, loader._MANIFEST_PATH: 0}
        results = []
        failures = []
        real_read_bytes = loader.Path.read_bytes
        real_parse = loader._parse_and_validate_pair
        parse_count = [0]

        def read(path):
            if path in counts:
                counts[path] += 1
                if path == loader._ARTIFACT_PATH and counts[path] == 1:
                    first_read_started.set()
                    self.assertTrue(tracking_lock.second_attempt.wait(5))
                return real_read_bytes(path)
            return real_read_bytes(path)

        def parse(artifact_bytes, manifest_bytes):
            parse_count[0] += 1
            return real_parse(artifact_bytes, manifest_bytes)

        def worker():
            try:
                results.append(loader.load_committed_data_hall_provenance())
            except BaseException as error:  # report thread failures in main thread
                failures.append(error)

        old_lock = loader._CACHE_LOCK
        loader._CACHE_LOCK = tracking_lock
        try:
            with mock.patch.object(loader.Path, "read_bytes", autospec=True, side_effect=read):
                with mock.patch.object(loader, "_parse_and_validate_pair", side_effect=parse):
                    first = threading.Thread(target=worker)
                    second = threading.Thread(target=worker)
                    first.start()
                    self.assertTrue(first_read_started.wait(5))
                    second.start()
                    first.join(10)
                    second.join(10)
        finally:
            loader._CACHE_LOCK = old_lock
            loader._reset_cache_for_test()
        self.assertFalse(first.is_alive())
        self.assertFalse(second.is_alive())
        self.assertEqual(failures, [])
        self.assertEqual(len(results), 2)
        self.assertIs(results[0], results[1])
        self.assertEqual(counts[loader._ARTIFACT_PATH], 1)
        self.assertEqual(counts[loader._MANIFEST_PATH], 1)
        self.assertEqual(parse_count[0], 1)

    def test_failed_first_load_is_retryable(self):
        loader._reset_cache_for_test()
        real_read_bytes = loader.Path.read_bytes
        failed = [False]

        def fail_once(path):
            if path == loader._ARTIFACT_PATH and not failed[0]:
                failed[0] = True
                return b"not the committed artifact"
            return real_read_bytes(path)

        with mock.patch.object(loader.Path, "read_bytes", autospec=True, side_effect=fail_once):
            with self.assertRaises(loader.DataHallIntegrityError):
                loader.load_committed_data_hall_provenance()
            self.assertIsNone(loader._DATABASE)
            retried = loader.load_committed_data_hall_provenance()
        self.assertEqual(retried.frames[4].data_hall, 9)
        self.assertIs(loader.load_committed_data_hall_provenance(), retried)

    def test_loader_source_is_independent(self):
        source = Path(loader.__file__).read_text(encoding="utf-8")
        tree = ast.parse(source)
        imported = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
        self.assertFalse(any(name.startswith("scripts.") for name in imported))
        for forbidden in (
            "freeze_iso_irrep_data_hall", "derive_iso_irrep_data_hall",
            "iso_irrep_exact", "generated_data", "SG_DATA_HALL",
            "hall_operations.json",
        ):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
