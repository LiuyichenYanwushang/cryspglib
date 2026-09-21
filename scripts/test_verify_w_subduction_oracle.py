"""Offline regressions for the other-wave-vector oracle gate.

Run with:
    python3 -m unittest discover -s scripts -p test_verify_w_subduction_oracle.py

`verify_w_subduction_oracle.py` needs the extracted ISOTROPY binary to talk to
the real oracle, so these tests drive the parsing and the comparison logic
directly: a stub archive stands in for the pinned tables and a captured
`DISPLAY ISOTROPY` table stands in for the program.  Every gate is exercised by
mutating that baseline, so a gate that stops firing fails here instead of
passing silently.
"""

import sys
import unittest

import verify_w_subduction_oracle as gate


def stub_archive(records, w_labels=("DT1", "DT2", "SM1")):
    """A `PinnedArchive` whose record lookup is a plain dict.

    `records` maps `(parent sg, compact label)` to
    `{(subgroup, direction label): (ordinal, [(source label, frequency)])}`.
    """
    archive = object.__new__(gate.PinnedArchive)
    archive.w_label = list(w_labels)
    archive.other_wave_vector = {
        ordinal: rows for table in records.values() for ordinal, rows in table.values()
    }
    archive.records_of_irrep = lambda sg, compact: dict(records.get((sg, compact), {}))
    return archive


def oracle_text(rows):
    """Build a captured `DISPLAY ISOTROPY` table.

    `rows` is a sequence of `(subgroup, direction label, frequency string)`.
    The trailing `*` is the program's prompt, which the parser must strip.
    """
    lines = ["Isotropy, Version 9.6.1, Jan 2022", "*Subgroup      Size Dir Frequency"]
    for subgroup, direction, frequencies in rows:
        lines.append(f"{subgroup} P2_12_12   16   {direction} {frequencies}")
    lines.append("*")
    return "\n".join(lines)


# Two rows of SG 196 `W1`, as the program prints them: the frequency list mixes
# same-space-group entries (`GM1`, `GM2GM3`, `X1`, `W1`) with the
# other-wave-vector sources (`DT1`, `DT2`, `SM1`).  The `C11` row is printed
# twice because the program prints one row per domain.
BASELINE_ORACLE = oracle_text(
    [
        (18, "C5", "1 GM1, 2 GM2GM3, 1 X1, 2 W1, 1 DT1, 1 DT2, 1 SM1"),
        (18, "C11", "1 GM1, 2 GM2GM3, 3 X4, 2 W1, 2 DT2, 1 SM1"),
        (18, "C11", "1 GM1, 2 GM2GM3, 3 X4, 2 W1, 2 DT2, 1 SM1"),
    ]
)

BASELINE_RECORDS = {
    (196, "W1"): {
        (18, "C5"): (10031, [("DT1", 1), ("DT2", 1), ("SM1", 1)]),
        (18, "C11"): (10033, [("DT2", 2), ("SM1", 1)]),
    }
}


class ParseTests(unittest.TestCase):
    def test_frequency_list_parses_counts_and_labels(self):
        self.assertEqual(
            gate.parse_frequency_list("1 GM1, 2 GM2GM3, 12 DT5"),
            [("GM1", 1), ("GM2GM3", 2), ("DT5", 12)],
        )

    def test_empty_frequency_list_is_empty(self):
        self.assertEqual(gate.parse_frequency_list("   "), [])

    def test_unparsable_entry_is_rejected(self):
        with self.assertRaises(ValueError):
            gate.parse_frequency_list("GM1")

    def test_rows_are_keyed_by_subgroup_and_direction(self):
        rows = gate.parse_isotropy_rows(BASELINE_ORACLE)
        self.assertEqual(sorted(rows), [(18, "C11"), (18, "C5")])
        self.assertEqual(len(rows[(18, "C11")]), 2)
        self.assertEqual(
            rows[(18, "C5")][0],
            [
                ("GM1", 1),
                ("GM2GM3", 2),
                ("X1", 1),
                ("W1", 2),
                ("DT1", 1),
                ("DT2", 1),
                ("SM1", 1),
            ],
        )

    def test_header_and_prompt_lines_are_ignored(self):
        rows = gate.parse_isotropy_rows(BASELINE_ORACLE)
        self.assertEqual(len(rows), 2)

    def test_wrapped_frequency_list_is_joined(self):
        """A long list wraps onto indented lines; the tail must not be dropped."""
        text = "\n".join(
            [
                "Subgroup      Size Dir Frequency",
                "18 P2_12_12   16   C5  1 GM1, 2 GM2GM3, 1 DT1,",
                "                        1 DT2, 1 SM1",
                "18 P2_12_12   16   C11 2 DT2, 1 SM1",
                "*",
            ]
        )
        rows = gate.parse_isotropy_rows(text)
        self.assertEqual(
            rows[(18, "C5")][0],
            [("GM1", 1), ("GM2GM3", 2), ("DT1", 1), ("DT2", 1), ("SM1", 1)],
        )
        self.assertEqual(rows[(18, "C11")][0], [("DT2", 2), ("SM1", 1)])


class CompareTests(unittest.TestCase):
    def compare(self, records=None, text=BASELINE_ORACLE, w_labels=("DT1", "DT2", "SM1")):
        archive = stub_archive(records or BASELINE_RECORDS, w_labels=w_labels)
        oracle = gate.parse_isotropy_rows(text)
        return gate.compare_group(archive, 196, "W1", oracle, set(archive.w_label))

    def test_matching_rows_pass(self):
        compared, failures = self.compare()
        self.assertEqual(failures, [])
        self.assertEqual(compared["oracle_rows"], 2)
        # Seven printed rows: `C5` once plus the two domain rows of `C11`.
        self.assertEqual(compared["oracle_w_rows"], 7)
        self.assertEqual(compared["pinned_w_rows"], 5)
        self.assertEqual(compared["records"], 2)

    def test_same_k_entries_do_not_enter_the_comparison(self):
        """Dropping the label filter would make the baseline fail."""
        compared, failures = self.compare()
        self.assertEqual(failures, [])
        self.assertEqual(compared["oracle_w_rows"], 7)

    def test_missing_pinned_row_fails(self):
        records = {
            (196, "W1"): {
                (18, "C5"): (10031, [("DT1", 1), ("DT2", 1)]),
                (18, "C11"): (10033, [("DT2", 2), ("SM1", 1)]),
            }
        }
        _, failures = self.compare(records=records)
        self.assertEqual(len(failures), 1)
        self.assertIn("pinned 1 DT1, 1 DT2 | oracle 1 DT1, 1 DT2, 1 SM1", failures[0])

    def test_wrong_frequency_fails(self):
        records = {
            (196, "W1"): {
                (18, "C5"): (10031, [("DT1", 2), ("DT2", 1), ("SM1", 1)]),
                (18, "C11"): (10033, [("DT2", 2), ("SM1", 1)]),
            }
        }
        _, failures = self.compare(records=records)
        self.assertEqual(len(failures), 1)
        self.assertIn("pinned 2 DT1, 1 DT2, 1 SM1", failures[0])

    def test_oracle_w_row_without_a_pinned_record_fails(self):
        records = {(196, "W1"): {(18, "C11"): (10033, [("DT2", 2), ("SM1", 1)])}}
        _, failures = self.compare(records=records)
        self.assertEqual(len(failures), 1)
        self.assertIn("the pinned table has no such record", failures[0])

    def test_pinned_row_without_an_oracle_row_fails(self):
        records = {
            (196, "W1"): {
                (18, "C5"): (10031, [("DT1", 1), ("DT2", 1), ("SM1", 1)]),
                (18, "C11"): (10033, [("DT2", 2), ("SM1", 1)]),
                (99, "P1"): (10099, [("DT1", 1)]),
            }
        }
        _, failures = self.compare(records=records)
        self.assertEqual(len(failures), 1)
        self.assertIn("have no oracle row", failures[0])

    def test_domain_duplicate_of_the_oracle_row_is_accepted(self):
        """The program prints one row per domain; the pinned rows are a set."""
        _, failures = self.compare()
        self.assertEqual(failures, [])


class CheckTests(unittest.TestCase):
    def test_check_uses_the_live_oracle_per_irrep(self):
        archive = stub_archive(BASELINE_RECORDS)
        archive.records_with_other_wave_vectors = lambda: [10031, 10033]
        archive.context = lambda ordinal: {
            "parent": 196,
            "compact_label": "W1",
            "irrep_label": "W1",
            "subgroup": 18,
            "direction_label": "C5" if ordinal == 10031 else "C11",
        }
        calls = []

        def fake_oracle(sg, compact_label):
            calls.append((sg, compact_label))
            return gate.parse_isotropy_rows(BASELINE_ORACLE)

        original = gate.run_oracle
        gate.run_oracle = fake_oracle
        try:
            totals, labels, failures = gate.check(archive)
        finally:
            gate.run_oracle = original
        self.assertEqual(calls, [(196, "W1")])
        self.assertEqual(failures, [])
        self.assertEqual(totals["groups"], 1)
        self.assertEqual(totals["records"], 2)
        self.assertEqual(labels, ["DT1", "DT2", "SM1"])


if __name__ == "__main__":
    sys.exit(unittest.main())
