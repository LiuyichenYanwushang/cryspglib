#!/usr/bin/env python3
"""Offline tests for the task-9 per-ordinal embedding settings.

No oracle process is started here: the committed module is checked against the
pinned isotropy tables in ``isotropy_subgroup/iso.zip`` (read, never an
extracted checkout), and a small set of **frozen oracle anchors** pins the
derivations.  The anchors store the exact printed basis and origin of a few
records, as printed by the pinned official ``iso`` binary; the derivation code
has to reproduce the committed ``U``/``delta`` from those values plus the
pinned machine table.  ``python3 scripts/generate_subduction_settings.py
--check`` is the online half of the same gate (it re-queries the oracle).

Usage::

    python3 -m unittest discover -s scripts -p test_subduction_settings.py
"""

from __future__ import annotations

import os
import sys
import unittest
from unittest import mock
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir))

sys.path.insert(0, HERE)
import generate_subduction_settings as gen  # noqa: E402

OUTPUT = os.path.join(ROOT, "src", "irrep", "subduction_settings_data.rs")

# Exact ``DISPLAY ISOTROPY`` basis/origin of five records, transcribed from the
# pinned `iso` 9.6.1 oracle (`SET I ALL OR 1`).  They span: a cubic parent with
# an F-centred supercell, a rhombohedral parent with the cyclic setting, the
# shear-free identity case, and the one record whose isotropy setting and Hall
# row use different ITA origin choices (#126: OR 1 prints (1,1,1) where OR 2
# prints (1/4,1/4,1/4) for the same cell).
ANCHORS = {
    314: {
        "parent": 16,
        "ml": "R1",
        "label": "P1",
        "child": 22,
        "basis": [["2", "0", "0"], ["0", "2", "0"], ["0", "0", "2"]],
        "origins": {1: ["0", "0", "0"], 2: ["0", "0", "0"]},
        "setting": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "child_shift": [0, 0, 0, 1],
    },
    7651: {
        "parent": 167,
        "ml": "GM3+",
        "label": "P1",
        "child": 15,
        "basis": [["1/3", "2/3", "2/3"], ["-1", "0", "0"], ["-1/3", "-2/3", "1/3"]],
        "origins": {1: ["-1/6", "1/6", "1/6"], 2: ["-1/6", "1/6", "1/6"]},
        "setting": [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
        "child_shift": [0, 0, 0, 1],
    },
    12400: {
        "parent": 221,
        "ml": "GM4+",
        "label": "P1",
        "child": 83,
        "basis": [["0", "0", "1"], ["0", "-1", "0"], ["1", "0", "0"]],
        "origins": {1: ["0", "0", "0"], 2: ["0", "0", "0"]},
        "setting": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "child_shift": [0, 0, 0, 1],
    },
    13729: {
        "parent": 225,
        "ml": "W5",
        "label": "4D22",
        "child": 8,
        "basis": [["2", "0", "2"], ["-2", "0", "2"], ["0", "-1", "0"]],
        "origins": {1: ["0", "0", "0"], 2: ["0", "0", "0"]},
        "setting": [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
        "child_shift": [0, 0, 0, 1],
    },
    6213: {
        "parent": 139,
        "ml": "M1-",
        "label": "P1",
        "child": 126,
        "basis": [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]],
        "origins": {1: ["1", "1", "1"], 2: ["1/4", "1/4", "1/4"]},
        "setting": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "child_shift": [1, 1, 1, 4],
    },
}


def frozen_basis(anchor):
    return [[Fraction(value) for value in row] for row in anchor["basis"]]


def frozen_origin(anchor, choice=1):
    return [Fraction(value) for value in anchor["origins"][choice]]


class CommittedModuleTest(unittest.TestCase):
    """The committed metadata must be addressable and internally consistent."""

    @classmethod
    def setUpClass(cls):
        with open(OUTPUT, encoding="utf-8") as handle:
            cls.text = handle.read()
        cls.entries = gen.parse_committed(cls.text)
        cls.records = gen.load_machine_records()
        cls.selected = gen.select_records(cls.records)

    def test_committed_entries_address_the_pinned_records(self):
        self.assertEqual(len(self.entries), 75)
        self.assertEqual(len({entry["ordinal"] for entry in self.entries}), 75)
        for entry in self.entries:
            self.assertTrue(
                gen.entry_addresses_record(entry, self.records),
                f"ordinal {entry['ordinal']} does not address its own record",
            )
        # A different ordinal is a different (parent, subgroup) pair and must
        # not pass the same identity check; the same holds for a wrong child.
        other = dict(self.entries[0])
        other["ordinal"] = 370  # the pinned record is (19, 19), not (16, 22)
        self.assertFalse(gen.entry_addresses_record(other, self.records))
        swapped = dict(self.entries[0])
        swapped["child"] = 1
        self.assertFalse(gen.entry_addresses_record(swapped, self.records))
        # An out-of-range ordinal is rejected instead of read past the table.
        outside = dict(self.entries[0])
        outside["ordinal"] = len(self.records)
        self.assertFalse(gen.entry_addresses_record(outside, self.records))

    def test_comments_carry_the_source_record_identity(self):
        sizes = {}
        for record in self.selected:
            sizes[record["ordinal"]] = abs(gen.geometry.det3(record["basis"]))
        for entry in self.entries:
            record = self.records[entry["ordinal"]]
            expected = (
                f"{record['parent']} {record['ml']} {record['label']} -> "
                f"#{record['child']} (size {sizes[entry['ordinal']]})"
            )
            self.assertTrue(
                entry["comment"].startswith(expected),
                f"ordinal {entry['ordinal']}: comment {entry['comment']!r} does "
                f"not name {expected!r}",
            )

    def test_settings_are_unimodular_and_match_the_previous_rules(self):
        for entry in self.entries:
            self.assertTrue(gen.integers(entry["setting"]))
            self.assertEqual(abs(gen.geometry.det3(entry["setting"])), 1)
            extra = gen.EXTRA_RECORDS.get(entry["ordinal"])
            expected = extra[4] if extra is not None else gen.PREVIOUS_RULES[(entry["parent"], entry["child"])]
            self.assertEqual(
                entry["setting"],
                expected,
                f"ordinal {entry['ordinal']}",
            )

    def test_only_the_126_child_shift_is_non_zero_and_reduced(self):
        non_zero = [
            entry
            for entry in self.entries
            if list(entry["child_shift"]) != list(gen.ZERO_SHIFT)
        ]
        self.assertEqual(len(non_zero), 1)
        entry = non_zero[0]
        self.assertEqual((entry["parent"], entry["child"]), (139, 126))
        self.assertEqual(entry["child_shift"], [1, 1, 1, 4])
        self.assertIn("origin choice", entry["comment"])
        for entry in self.entries:
            x, y, z, denominator = entry["child_shift"]
            self.assertGreater(denominator, 0)
            divisor = denominator
            for value in (x, y, z):
                divisor = gen._gcd(divisor, value)
            self.assertEqual(divisor, 1, f"ordinal {entry['ordinal']} is unreduced")


class FrozenAnchorTest(unittest.TestCase):
    """Frozen oracle evidence must reproduce the committed derivations."""

    @classmethod
    def setUpClass(cls):
        with open(OUTPUT, encoding="utf-8") as handle:
            cls.entries = {
                entry["ordinal"]: entry
                for entry in gen.parse_committed(handle.read())
            }
        cls.records = gen.load_machine_records()

    def test_anchors_reproduce_U_and_delta(self):
        for ordinal, anchor in ANCHORS.items():
            record = self.records[ordinal]
            self.assertEqual(record["parent"], anchor["parent"])
            self.assertEqual(record["ml"], anchor["ml"])
            self.assertEqual(record["label"], anchor["label"])
            self.assertEqual(record["child"], anchor["child"])
            setting = gen.derive_setting(record, frozen_basis(anchor))
            self.assertEqual(setting, anchor["setting"], f"ordinal {ordinal}")
            self.assertEqual(
                setting,
                self.entries[ordinal]["setting"],
                f"ordinal {ordinal} disagrees with the committed module",
            )
            gen.check_origin(record, frozen_origin(anchor))
            # OR 1 always matches by construction, so the shift is exactly the
            # printed origin difference of the setting the Hall row belongs to.
            choices = [1]
            committed = self.entries[ordinal]["child_shift"]
            if list(committed) != list(gen.ZERO_SHIFT):
                choices = [2]
            delta = gen.setting_shift(
                record,
                frozen_basis(anchor),
                {1: frozen_origin(anchor, 1), 2: frozen_origin(anchor, 2)},
                choices,
            )
            self.assertEqual(
                list(gen.encode_shift(delta)),
                list(committed),
                f"ordinal {ordinal}",
            )

    def test_anchor_origin_disagreement_is_detected(self):
        record = self.records[6213]
        gen.check_origin(record, frozen_origin(ANCHORS[6213]))
        with self.assertRaises(gen.GenerationError):
            gen.check_origin(record, [Fraction(1, 4), Fraction(1, 4), Fraction(1, 4)])
        record_314 = self.records[314]
        with self.assertRaises(gen.GenerationError):
            gen.check_origin(record_314, [Fraction(0), Fraction(0), Fraction(1, 2)])


class DerivationGateTest(unittest.TestCase):
    def test_generated_hall_map_cannot_override_source_frame_provenance(self):
        import tempfile

        values = gen.load_data_hall()
        values[43] = 1
        with tempfile.NamedTemporaryFile(mode="w", suffix=".rs") as source:
            source.write("pub static SG_DATA_HALL: [u16; 231] = [" +
                         ",".join(map(str, values)) + "];")
            source.flush()
            with mock.patch.object(gen, "GENERATED_DATA", source.name):
                with self.assertRaisesRegex(gen.GenerationError, "frame provenance"):
                    gen.load_data_hall()

    def test_generator_ignores_untrusted_local_extraction(self):
        import tempfile
        import zipfile

        def check_pinned_directory(directory):
            with zipfile.ZipFile(gen.census.ARCHIVE) as archive:
                for name in gen.census.ORACLE_MEMBERS:
                    with open(os.path.join(directory, name), "rb") as source:
                        self.assertEqual(source.read(), archive.read(name))
            self.assertNotEqual(directory, untrusted)
            return "verified"

        with tempfile.TemporaryDirectory() as untrusted:
            with open(os.path.join(untrusted, "iso"), "wb") as source:
                source.write(b"untrusted local binary")
            with mock.patch.object(gen, "ISO_DIR", untrusted), mock.patch.object(
                gen, "_build_entries", side_effect=check_pinned_directory
            ):
                self.assertEqual(gen.build_entries(), "verified")

    """The setting/shift gates reject the inputs they must reject."""

    def test_shear_setting_is_accepted(self):
        shear = [[1, 0, 0], [1, 1, 0], [0, 0, 1]]
        record = {"ordinal": 0, "parent": 221, "child": 83, "basis": shear}
        setting = gen.derive_setting(record, gen.IDENTITY)
        self.assertEqual(setting, shear)
        self.assertEqual(abs(gen.geometry.det3(setting)), 1)
        # A mirror (det -1) is a legitimate setting too; the gate is |det| = 1.
        mirror = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]
        record = {"ordinal": 0, "parent": 221, "child": 83, "basis": mirror}
        self.assertEqual(gen.derive_setting(record, gen.IDENTITY), mirror)

    def test_non_integer_and_non_unimodular_settings_are_rejected(self):
        record = {"ordinal": 0, "parent": 221, "child": 83, "basis": gen.IDENTITY}
        with self.assertRaises(gen.GenerationError):
            gen.check_unimodular(record, [[Fraction(1, 2), 0, 0], [0, 1, 0], [0, 0, 1]])
        with self.assertRaises(gen.GenerationError):
            gen.check_unimodular(record, [[2, 0, 0], [0, 1, 0], [0, 0, 1]])
        with self.assertRaises(gen.GenerationError):
            gen.check_unimodular(record, [[0, 0, 0], [0, 1, 0], [0, 0, 1]])
        # Bases that do not span the same lattice are refused before the setting
        # is even formed, rather than silently rescaling the cell.
        with self.assertRaises(gen.GenerationError):
            gen.derive_setting(record, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    def test_setting_selection_needs_exactly_one_origin(self):
        record = {"ordinal": 0}
        identity = ((1, 0, 0, 0, 1, 0, 0, 0, 1), (Fraction(0), Fraction(0), Fraction(0)))
        translated = (
            (1, 0, 0, 0, 1, 0, 0, 0, 1),
            (Fraction(1, 2), Fraction(0), Fraction(0)),
        )
        expansions = {
            1: {(1, 0, 0, 0, 1, 0, 0, 0, 1): {gen.canonical_translation([Fraction(0)] * 3)}},
            2: {(1, 0, 0, 0, 1, 0, 0, 0, 1): {gen.canonical_translation([Fraction(1, 2), Fraction(0), Fraction(0)])}},
        }
        self.assertEqual(gen.select_setting(record, expansions, [identity]), [1])
        self.assertEqual(gen.select_setting(record, expansions, [translated]), [2])
        with self.assertRaises(gen.GenerationError):
            gen.select_setting(record, expansions, [])
        # Two matching settings with different printed origins are ambiguous.
        with self.assertRaises(gen.GenerationError):
            gen.setting_shift(
                record,
                gen.IDENTITY,
                {1: [Fraction(0)] * 3, 2: [Fraction(1, 4)] * 3},
                [1, 2],
            )


if __name__ == "__main__":
    unittest.main()
