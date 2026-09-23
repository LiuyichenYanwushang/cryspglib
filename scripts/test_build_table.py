#!/usr/bin/env python3
"""Regression tests for the frozen-table assembler.

`scripts/task9/build_table.py` is the only writer of
``src/irrep/subduction_settings_data.rs``.  These tests exercise the states that
used to pass silently: an empty or truncated census, duplicated records, a
`--check` mismatch, and -- the case the first round of tests missed -- a **new
output path**, where no existing module can supply the expected ordinal set.
They run entirely on temporary files, pin their own `--baseline` fixture, and
never touch the committed module.

Usage::

    python3 -m unittest discover -s scripts -p test_build_table.py
"""

from __future__ import annotations

import io
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, "task9"))

import build_table  # noqa: E402

MODULE_HEAD = "pub static FROZEN_EMBEDDING_SETTINGS: &[FrozenEmbeddingSetting] = &[\n"


def module(entries):
    """A minimal module with the fields `parse_committed` reads."""
    lines = [MODULE_HEAD]
    for ordinal, parent, child in entries:
        lines.append(f"    ({ordinal}, {parent}, {child}, IDENTITY, 1, NO_SHIFT),\n")
    lines.append("];\n")
    return "".join(lines)


def census_line(ordinal, parent, child, u=("1", "0", "0", "0", "1", "0", "0", "0", "1")):
    return json.dumps(
        {
            "type": "record",
            "ordinal": ordinal,
            "parent": parent,
            "child": child,
            "u": [list(u[0:3]), list(u[3:6]), list(u[6:9])],
            "status": "printed",
        }
    )


class AssemblerTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.census = os.path.join(self.directory.name, "census.jsonl")
        self.shifts = os.path.join(self.directory.name, "shifts.json")
        self.empty = os.path.join(self.directory.name, "empty.jsonl")
        self.out = os.path.join(self.directory.name, "table.rs")
        self.baseline = os.path.join(self.directory.name, "baseline.rs")
        self.expected = os.path.join(self.directory.name, "expected.json")
        with open(self.shifts, "w", encoding="utf-8") as handle:
            json.dump({"solved": {}}, handle)
        open(self.empty, "w", encoding="utf-8").close()

    def write_census(self, lines):
        with open(self.census, "w", encoding="utf-8") as handle:
            for line in lines:
                handle.write(line + "\n")

    def write_baseline(self, entries):
        with open(self.baseline, "w", encoding="utf-8") as handle:
            handle.write(module(entries))

    def run_assembler(self, *extra):
        argv = [
            "--census", self.census,
            "--shifts", self.shifts,
            "--empty", self.empty,
            "--out", self.out,
            "--baseline", self.baseline,
            *extra,
        ]
        stdout, stderr = io.StringIO(), io.StringIO()
        code = 0
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                build_table.main(argv)
            except SystemExit as error:
                code = error.code if isinstance(error.code, int) else 1
                if isinstance(error.code, str):
                    # An uncaught SystemExit would print this; keep it visible so
                    # the assertions can check the refusal reason.
                    stderr.write(error.code + "\n")
        return code, stdout.getvalue(), stderr.getvalue()

    def test_an_empty_input_cannot_write_a_table(self):
        self.write_census([])
        with open(self.out, "w", encoding="utf-8") as handle:
            handle.write(module([(0, 1, 1), (1, 1, 1)]))
        self.write_baseline([(0, 1, 1), (1, 1, 1)])
        with open(self.out, encoding="utf-8") as handle:
            before = handle.read()
        code, _, _ = self.run_assembler()
        self.assertNotEqual(code, 0, "an empty census must not assemble")
        with open(self.out, encoding="utf-8") as handle:
            self.assertEqual(
                handle.read(),
                before,
                "the output must be untouched when validation fails",
            )

    def test_a_truncated_census_cannot_shrink_the_table(self):
        self.write_census([census_line(0, 1, 1)])
        with open(self.out, "w", encoding="utf-8") as handle:
            handle.write(module([(0, 1, 1), (1, 1, 1)]))
        self.write_baseline([(0, 1, 1), (1, 1, 1)])
        code, _, _ = self.run_assembler()
        self.assertNotEqual(code, 0)
        # An explicit `--partial` is the documented escape hatch.
        code, _, _ = self.run_assembler("--partial")
        self.assertEqual(code, 0)
        with open(self.out, encoding="utf-8") as handle:
            self.assertIn("(0, 1, 1,", handle.read())

    def test_a_duplicated_record_is_a_conflict_not_an_overwrite(self):
        self.write_census([census_line(0, 1, 1), census_line(0, 1, 1)])
        with open(self.out, "w", encoding="utf-8") as handle:
            handle.write(module([(0, 1, 1)]))
        self.write_baseline([(0, 1, 1)])
        code, _, stderr = self.run_assembler()
        self.assertNotEqual(code, 0)
        self.assertIn("duplicated", stderr)

    def test_a_complete_input_writes_and_check_agrees(self):
        self.write_census([census_line(0, 1, 1), census_line(1, 1, 1)])
        with open(self.out, "w", encoding="utf-8") as handle:
            handle.write(module([(0, 1, 1), (1, 1, 1)]))
        self.write_baseline([(0, 1, 1), (1, 1, 1)])
        code, _, _ = self.run_assembler()
        self.assertEqual(code, 0)
        code, stdout, _ = self.run_assembler("--check")
        self.assertEqual(code, 0, stdout)
        self.assertIn("identical", stdout)
        # A changed input makes `--check` fail instead of rewriting silently.
        self.write_census([census_line(0, 1, 1), census_line(1, 2, 1)])
        code, _, _ = self.run_assembler("--check")
        self.assertNotEqual(code, 0)
        code, _, _ = self.run_assembler()
        self.assertEqual(code, 0)
        code, _, _ = self.run_assembler("--check")
        self.assertEqual(code, 0)

    def test_a_new_output_needs_a_known_expected_universe(self):
        # The case the first round missed: `--out` does not exist yet, so only
        # `--baseline` (or `--expected`) can tell a complete table from a
        # truncated input.  An empty census must not write a zero-entry "full"
        # table, and a rejected build must not create the output at all.
        self.write_census([])
        self.write_baseline([(0, 1, 1), (1, 1, 1)])
        self.assertFalse(os.path.exists(self.out))
        code, _, stderr = self.run_assembler()
        self.assertNotEqual(code, 0, "an unknown universe must not permit a write")
        self.assertFalse(os.path.exists(self.out), "a rejected build wrote its output")
        # A complete input for the same fresh output path is accepted.
        self.write_census([census_line(0, 1, 1), census_line(1, 1, 1)])
        code, _, _ = self.run_assembler()
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(self.out))
        # With no baseline either, the assembler refuses instead of guessing.
        os.unlink(self.out)
        missing = os.path.join(self.directory.name, "absent.rs")
        code, _, stderr = self.run_assembler("--baseline", missing)
        self.assertNotEqual(code, 0)
        self.assertIn("cannot determine the expected ordinal set", stderr)
        self.assertFalse(os.path.exists(self.out))
        # `--partial` is the explicit escape hatch for an experimental table.
        code, _, _ = self.run_assembler("--baseline", missing, "--partial")
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(self.out))

    def test_output_and_baseline_must_agree_on_the_ordinal_set(self):
        self.write_census([census_line(0, 1, 1)])
        with open(self.out, "w", encoding="utf-8") as handle:
            handle.write(module([(0, 1, 1)]))
        self.write_baseline([(0, 1, 1), (1, 1, 1)])
        code, _, stderr = self.run_assembler()
        self.assertNotEqual(code, 0)
        self.assertIn("disagree", stderr)
        # An explicit `--expected` states the intended universe and wins.
        with open(self.expected, "w", encoding="utf-8") as handle:
            json.dump([0], handle)
        code, _, _ = self.run_assembler("--expected", self.expected)
        self.assertEqual(code, 0)

    def test_load_shifts_accepts_every_pipeline_shape(self):
        for payload, expected in (
            ({"solved": {"7": [1, 0, 0, 2]}}, {7: [1, 0, 0, 2]}),
            ({"accepted": [7], "shifts": {"7": [1, 0, 0, 2]}}, {7: [1, 0, 0, 2]}),
            ({"accepted": [], "shifts": {"7": [1, 0, 0, 2]}}, {}),
            ({"7": {"shift": [1, 0, 0, 2]}}, {7: [1, 0, 0, 2]}),
        ):
            with open(self.shifts, "w", encoding="utf-8") as handle:
                json.dump(payload, handle)
            self.assertEqual(build_table.load_shifts(self.shifts), expected, payload)


if __name__ == "__main__":
    unittest.main()
