"""Regressions for generated Rust tables before and after cargo fmt.

Run: python3 -m unittest discover -s scripts -p test_rustfmt_tables.py
"""
import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import freeze_direction_descriptors as descriptors
import freeze_w_little_characters as little
import generate_subduction_settings as settings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "task9"))
import build_table
import probe_table

COMPACT = build_table.HEADER + """    (26, 3, 3, [[1, 222222222, 1], [-1, 2, -111111111], [-1, 0, 1]], 2, NO_SHIFT),
    (27, 3, 3, IDENTITY, 1, [1, -1, 0, 4]),
];
"""


def rustfmt(text):
    return subprocess.run(
        ["rustfmt", "--edition", "2024"], input=text, text=True,
        capture_output=True, check=True,
    ).stdout


class RustfmtTableTests(unittest.TestCase):
    def test_all_setting_readers_preserve_multiline_rows(self):
        expected = settings.parse_committed(COMPACT)
        self.assertEqual(len(expected), 2)
        for text in (COMPACT, rustfmt(COMPACT)):
            with self.subTest(formatted=text != COMPACT):
                self.assertEqual(settings.parse_committed(text), expected)
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "table.rs"
                    path.write_text(text, encoding="utf-8")
                    self.assertEqual(build_table.parse_committed(path), {26: (3, 3), 27: (3, 3)})
                    parsed = probe_table.parse(path)
                    self.assertEqual([row[0] for row in parsed], [26, 27])
                    self.assertEqual(parsed[0][3:], ([[1, 222222222, 1], [-1, 2, -111111111], [-1, 0, 1]], 2, [0, 0, 0, 1]))

    def test_full_settings_table_roundtrips_without_data_changes(self):
        text = (ROOT / "src/irrep/subduction_settings_data.rs").read_text(encoding="utf-8")
        rows = settings.parse_committed(text)
        settings.check_ordinal_coverage(rows, 15239)
        entries = {
            row["ordinal"]: (row["parent"], row["child"], row["setting"],
                             row["setting_denominator"], row["child_shift"])
            for row in rows
        }
        self.assertEqual(build_table.render(entries), text)

    def test_malformed_setting_rows_are_not_silently_skipped(self):
        for body in (
            "garbage,",
            "(26, 3, 3, IDENTITY, 1),",
            "(26, 3, 3, UNKNOWN, 1, NO_SHIFT),",
            "(26, 3, 3, [[1]], 1, NO_SHIFT),",
            "(26, 3, 3, IDENTITY, 0, NO_SHIFT),",
            "(26, 3, 3, IDENTITY, 1, [0, 0, 0, 0]),",
            "(26, 3, 3, IDENTITY, 1, NO_SHIFT), broken,",
        ):
            with self.subTest(body=body), self.assertRaises(settings.GenerationError):
                settings.parse_committed(build_table.HEADER + body + "\n];\n")

    def test_assembler_keeps_duplicate_detection_after_formatting(self):
        text = rustfmt(COMPACT.replace("(27,", "(26,"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "table.rs"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(SystemExit, "appears twice"):
                build_table.parse_committed(path)

    def test_settings_emitter_is_rustfmt_stable(self):
        text = build_table.render({
            15238: (230, 227, [[1, 222222222, 1], [-1, 2, -111111111], [-1, 0, 1]], 2, [1, -1, 0, 4]),
        })
        self.assertEqual(text, rustfmt(text))
        self.assertEqual(settings.parse_committed(text)[0]["child_shift"], [1, -1, 0, 4])

    def test_descriptor_emitter_matches_formatted_table_and_rejects_tampering(self):
        json_path = ROOT / "scripts/data/direction_descriptors_v1.json"
        text = descriptors.render(str(json_path))
        self.assertEqual(text, rustfmt(text))
        self.assertEqual(text, (ROOT / "src/irrep/direction_descriptors_data.rs").read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "table.rs"
            for contents, expected in ((text, 0), (text.replace('"P1"', '"WRONG"', 1), 1)):
                path.write_text(contents, encoding="utf-8")
                argv = ["freeze_direction_descriptors.py", "--check", "--json", str(json_path), "--output", str(path)]
                with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(descriptors.main(), expected)
                self.assertEqual(path.read_text(encoding="utf-8"), contents, "--check must never rewrite the input")

    def test_little_character_emitter_is_rustfmt_stable(self):
        table = {
            "space_group": 196, "label": "DT1", "k_label": "DT",
            "direction": ["0", "2", "0"], "dimension": 1,
            "operations": [{"element": "x,y,z", "rotation": settings.IDENTITY,
                            "translation": ["0", "0", "0"], "character": ["1", "0"]}],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "table.rs"
            self.assertEqual(little.emit_rust([table], [], path), 1)
            text = path.read_text(encoding="utf-8")
            self.assertEqual(text, rustfmt(text))


if __name__ == "__main__":
    unittest.main()
