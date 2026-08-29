#!/usr/bin/env python3
"""Tests for the repo-only Rust/Python magnetic database parity frame."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import struct
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
import spglib_magnetic_provenance as provenance
import spglib_magnetic_runtime_parity as parity


REPOSITORY = Path(__file__).resolve().parents[1]
EXPECTED_PAYLOAD_LENGTHS = (
    2124, 4248, 32588, 64474, 13216, 4248, 237888, 306732,
    832608, 162920, 90776, 1023740, 87054,
)


def _first_difference(left: bytes, right: bytes) -> int | None:
    limit = min(len(left), len(right))
    for index in range(limit):
        if left[index] != right[index]:
            return index
    if len(left) != len(right):
        return limit
    return None


def _operation_bytes(payload: bytes, offset: int, magnetic: bool):
    width = 13 if magnetic else 12
    end = offset + width
    if end > len(payload):
        raise AssertionError("operation record is truncated")
    rotation = tuple(
        tuple(struct.unpack_from("<b", payload, offset + 3 * row + column)[0]
              for column in range(3))
        for row in range(3)
    )
    translation = tuple(payload[offset + 9 + index] for index in range(3))
    if any(value >= provenance.TRANSLATION_DENOMINATOR for value in translation):
        raise AssertionError("operation translation numerator is out of range")
    time_reversal = payload[offset + 12] if magnetic else None
    if magnetic and time_reversal not in (0, 1):
        raise AssertionError("operation time reversal is out of range")
    return rotation, translation, time_reversal, end


def _parse_variable_operations(
    payload: bytes, record_count: int, magnetic: bool, triple_header: bool = True
):
    header_width = 6 if triple_header else 4
    operation_width = 13 if magnetic else 12
    records = []
    offset = 0
    for _ in range(record_count):
        if offset + header_width > len(payload):
            raise AssertionError("variable section record header is truncated")
        if triple_header:
            uni_or_hall, hall_or_count, count = struct.unpack_from(
                "<HHH", payload, offset
            )
        else:
            uni_or_hall, count = struct.unpack_from("<HH", payload, offset)
            hall_or_count = count
        offset += header_width
        operations = []
        for _ in range(count):
            operation = _operation_bytes(payload, offset, magnetic)
            operations.append(operation[:3])
            offset += operation_width
        records.append((uni_or_hall, hall_or_count, count, tuple(operations)))
    if offset != len(payload):
        raise AssertionError("variable section has trailing bytes")
    return tuple(records)


class RuntimeParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = provenance.load_committed_provenance()
        cls.expected = parity.build_expected_frame(cls.database)

    def test_python_frame_is_deterministic_and_strictly_structured(self):
        self.assertEqual(self.expected, parity.build_expected_frame(self.database))
        sections = parity.parse_frame(self.expected)
        self.assertEqual(len(sections), parity.SECTION_COUNT)
        self.assertEqual(tuple(row[0] for row in sections), parity.SECTION_TAGS)
        self.assertEqual(tuple(row[1] for row in sections), parity.SECTION_COUNTS)
        self.assertEqual(
            tuple(len(row[2]) for row in sections), EXPECTED_PAYLOAD_LENGTHS
        )
        self.assertEqual(
            tuple(hashlib.sha256(row[2]).hexdigest() for row in sections),
            parity.GOLDEN_SECTION_PAYLOAD_SHA256,
        )
        self.assertEqual(len(self.expected), parity.GOLDEN_FRAME_LENGTH)
        self.assertEqual(
            hashlib.sha256(self.expected).hexdigest(), parity.GOLDEN_FRAME_SHA256
        )
        self.assertEqual(self.expected[:8], parity.MAGIC)
        self.assertEqual(struct.unpack_from("<II", self.expected, 8), (1, 13))

    def test_fixed_and_variable_section_record_contract(self):
        sections = parity.parse_frame(self.expected)
        fixed_widths = {
            "SGNO": 4, "SGIX": 8, "SGRW": 4, "MUNI": 8,
            "MHLL": 8, "MIDX": 8, "MRAW": 4, "MALT": 28,
            "SDEC": 20,
        }
        for tag, _, payload in sections:
            if tag in fixed_widths:
                self.assertEqual(len(payload) % fixed_widths[tag], 0, tag)
        sapi = _parse_variable_operations(
            sections[10][2], sections[10][1], False, triple_header=False
        )
        mapi = _parse_variable_operations(sections[11][2], sections[11][1], True)
        tapi = _parse_variable_operations(sections[12][2], sections[12][1], False)
        self.assertEqual(len(sapi), 530)
        self.assertEqual(len(mapi), parity.MSG_ACTIVE_SPAN_COUNT)
        self.assertEqual(len(tapi), parity.MSG_ACTIVE_SPAN_COUNT)
        self.assertEqual(sum(row[2] for row in sapi), 7388)
        self.assertEqual(sum(row[2] for row in mapi), 76682)
        self.assertEqual(sum(row[2] for row in tapi), 5015)
        self.assertEqual(sum(row[2] - 1 for row in tapi), 536)
        self.assertEqual(tuple((row[0], row[1]) for row in mapi), tuple(
            (uni, hall)
            for uni in range(1, provenance.MSG_UNI_COUNT)
            for hall in self.database.halls_for_uni(uni)
        ))
        self.assertEqual(tuple((row[0], row[1]) for row in tapi), tuple(
            (uni, hall)
            for uni in range(1, provenance.MSG_UNI_COUNT)
            for hall in self.database.halls_for_uni(uni)
        ))

    def test_uni7_uni9_records_and_time_reversal_census(self):
        sections = parity.parse_frame(self.expected)
        mapi = _parse_variable_operations(sections[11][2], sections[11][1], True)
        by_key = {(row[0], row[1]): row for row in mapi}
        uni7 = by_key[(7, 2)][3]
        self.assertEqual(
            tuple((rotation, translation, time) for rotation, translation, time in uni7),
            (
                (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0), 0),
                (((-1, 0, 0), (0, -1, 0), (0, 0, -1)), (0, 0, 0), 0),
                (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 6), 1),
                (((-1, 0, 0), (0, -1, 0), (0, 0, -1)), (0, 0, 6), 1),
            ),
        )
        expected_uni9 = (
            (3, ((-1, 0, 0), (0, 1, 0), (0, 0, -1))),
            (4, ((-1, 0, 0), (0, -1, 0), (0, 0, 1))),
            (5, ((1, 0, 0), (0, -1, 0), (0, 0, -1))),
        )
        for hall, rotation in expected_uni9:
            operations = by_key[(9, hall)][3]
            self.assertEqual(len(operations), 4)
            self.assertEqual(operations[0][0], ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
            self.assertEqual(operations[1][0], rotation)
            self.assertEqual(operations[2][0], ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
            self.assertEqual(operations[3][0], rotation)
            self.assertTrue(all(translation == (0, 0, 0) for _, translation, _ in operations))
            self.assertEqual(tuple(row[2] for row in operations), (0, 0, 1, 1))
        unitary = sum(operation[2] == 0 for row in mapi for operation in row[3])
        antiunitary = sum(operation[2] == 1 for row in mapi for operation in row[3])
        self.assertEqual(unitary, 42035)
        self.assertEqual(antiunitary, 34647)
        self.assertEqual(unitary + antiunitary, 76682)

        tapi = _parse_variable_operations(sections[12][2], sections[12][1], False)
        tapi_by_key = {(row[0], row[1]): row for row in tapi}
        self.assertEqual(tapi_by_key[(7, 2)][2], 7)
        for hall in (3, 4, 5):
            self.assertEqual(tapi_by_key[(9, hall)][2], 1)
            self.assertEqual(tapi_by_key[(9, hall)][3][0][:2], (
                ((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0)
            ))

    def _run_rust_frame(self):
        environment = os.environ.copy()
        environment.pop("SPGLIB_DEBUG", None)
        environment.pop("SPGLIB_INFO", None)
        environment["SPGLIB_WARNING"] = "OFF"
        # The surrounding workspace may expose a read-only target directory;
        # an explicit temporary target keeps this repo-only gate runnable.
        environment.setdefault(
            "CARGO_TARGET_DIR", "/tmp/cryspglib-magnetic-parity-target"
        )
        result = subprocess.run(
            [
                "cargo", "run", "--release", "--quiet", "--example",
                "magnetic_database_parity_dump",
            ],
            cwd=str(REPOSITORY),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=300,
        )
        if result.returncode != 0:
            self.fail(
                "Rust parity example failed with status "
                f"{result.returncode}: {result.stderr.decode('utf-8', 'replace')}"
            )
        return result.stdout

    def _assert_frame_equal(self, expected: bytes, actual: bytes):
        expected_sections = parity.parse_frame(expected)
        try:
            actual_sections = parity.parse_frame(actual)
        except parity.FrameError as error:
            offset = _first_difference(expected, actual)
            self.fail(f"Rust frame parse failed near byte {offset}: {error}")
        if len(expected) != len(actual):
            self.fail(
                f"frame length mismatch: expected {len(expected)}, actual {len(actual)}"
            )
        for index, (left, right) in enumerate(zip(expected_sections, actual_sections)):
            if left[0] != right[0]:
                self.fail(f"section {index} tag mismatch: {left[0]} != {right[0]}")
            if left[1] != right[1]:
                self.fail(
                    f"section {left[0]} record count mismatch: {left[1]} != {right[1]}"
                )
            if len(left[2]) != len(right[2]):
                self.fail(
                    f"section {left[0]} payload length mismatch: "
                    f"{len(left[2])} != {len(right[2])}"
                )
            difference = _first_difference(left[2], right[2])
            if difference is not None:
                self.fail(f"section {left[0]} first payload difference at byte {difference}")
        self.assertEqual(expected, actual)

    @unittest.skipUnless(
        os.environ.get("CRYSPGLIB_RUN_MAGNETIC_DB_PARITY") == "1",
        "set CRYSPGLIB_RUN_MAGNETIC_DB_PARITY=1 for the Rust runtime gate",
    )
    def test_rust_runtime_matches_typed_expected_frame_twice(self):
        first = self._run_rust_frame()
        second = self._run_rust_frame()
        self._assert_frame_equal(self.expected, first)
        self._assert_frame_equal(self.expected, second)
        self.assertEqual(len(first), parity.GOLDEN_FRAME_LENGTH)
        self.assertEqual(len(second), parity.GOLDEN_FRAME_LENGTH)
        self.assertEqual(hashlib.sha256(first).hexdigest(), parity.GOLDEN_FRAME_SHA256)
        self.assertEqual(hashlib.sha256(second).hexdigest(), parity.GOLDEN_FRAME_SHA256)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
