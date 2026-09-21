"""Offline regressions for the task 9 oracle setting census collector.

Run with: python3 -m unittest discover -s scripts -p test_audit_subduction_settings.py

The pinned transcript below is real iso 9.6.1 output for SG 43 ``L1`` under
``SET I ALL OR 1``.  Only the fake-oracle test starts a subprocess; every other
test drives the parsers and the candidate classification directly.
"""

import os
import stat
import sys
import tempfile
import unittest
from fractions import Fraction
from unittest import mock

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import audit_subduction_settings as census  # noqa: E402

BANNER = (
    "Isotropy, Version 9.6.1, Jan 2022\n"
    "Harold T. Stokes, Dorian M. Hatch, and Branton J. Campbell\n"
    "Brigham Young University\n"
    "Current setting is International (new ed.) with conventional basis vectors."
)
HEADER = "***********Subgroup Size Dir Basis Vectors           Origin"
EMPTY = BANNER + "\n************"
SG43_L1 = "\n".join(
    [
        BANNER,
        HEADER,
        "1 P1     2    P1  (1/2,1,1/2),(0,-1/2,1/2),(1/2,-1/2,0)  (0,0,0)",
        "5 C2     4    P3  (1,1,0),(0,0,-2),(-1/2,1/2,0)          (0,0,0)",
        "9 Cc     4    P4  (0,1,1),(2,0,0),(0,1/2,-1/2)           (7/8,7/16,7/16)",
        "9 Cc     4    P5  (1,0,1),(0,2,0),(-1/2,0,1/2)            (7/16,7/8,7/16)",
        "1 P1     4    C1  (1/2,1/2,-1),(1/2,1/2,1),(1/2,-1/2,0)  (0,0,0)",
        "1 P1     4    C2  (1,1/2,1/2),(0,-1/2,1/2),(1,-1/2,-1/2) (0,0,0)",
        "1 P1     4    C3  (1/2,1,1/2),(1/2,-1,1/2),(1/2,0,-1/2)  (0,0,0)",
        "5 C2     8    C8  (0,2,0),(0,0,-2),(-1,1,0)              (0,0,0)",
        "1 P1     8    4D1 (0,1,1),(1,0,1),(1,1,0)                (0,0,0)",
        "*",
    ]
)
IDENTITY = [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]]


def machine_row(label="P1", child=1, size=1, basis=None, ordinal=1, parent=1,
                index=1, irrep="GM1"):
    return {
        "ordinal": ordinal,
        "parent": parent,
        "source_irrep_index": index,
        "irrep": irrep,
        "direction_label": label,
        "child": child,
        "size": size,
        "dim": 3,
        "direction": "(a,0,0)",
        "basis": basis or [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "origin": [Fraction(0), Fraction(0), Fraction(0)],
    }


def oracle_row(label="P1", child=1, size=1, basis=None, origin=None):
    return {
        "child": child,
        "symbol": "P1",
        "size": size,
        "label": label,
        "basis": [
            [Fraction(token) for token in vector]
            for vector in (basis or IDENTITY)
        ],
        "origin": [Fraction(token) for token in (origin or ["0", "0", "0"])],
        "raw": f"{child} P1 {size} {label} (test)",
    }


class ClassificationTests(unittest.TestCase):
    def test_actual_shear_is_accepted_as_candidate(self):
        row = machine_row(basis=[[1, 1, 0], [0, 1, 0], [0, 0, 1]])
        evidence = census.derive_record(row, oracle_row(), 1)
        self.assertEqual(evidence["status"], "candidate")
        self.assertEqual(evidence["detail"], "origin_exact_u_unimodular")
        self.assertTrue(evidence["u_integral"])
        self.assertEqual(evidence["u_det"], "1")
        self.assertNotEqual(evidence["u"], IDENTITY)  # a real shear, not identity
        self.assertEqual(evidence["u"], [["1", "1", "0"], ["0", "1", "0"], ["0", "0", "1"]])

    def test_wrong_origin_is_reported_and_not_fixed(self):
        oracle = oracle_row(origin=["1/2", "0", "0"])
        evidence = census.derive_record(machine_row(), oracle, 1)
        self.assertEqual(evidence["status"], "origin_mismatch")
        self.assertEqual(evidence["detail"], "u_integral_unimodular_origin_delta")
        self.assertTrue(evidence["u_integral"])  # only the origin is wrong
        self.assertFalse(evidence["origin_exact"])
        self.assertEqual(evidence["origin_delta"], ["1/2", "0", "0"])
        self.assertEqual(evidence["origin_printed"], ["1/2", "0", "0"])
        self.assertEqual(evidence["origin_converted"], ["0", "0", "0"])

    def test_nonintegral_u_is_nonunimodular(self):
        oracle = oracle_row(basis=[["2", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]])
        evidence = census.derive_record(machine_row(), oracle, 1)
        self.assertEqual(evidence["status"], "nonunimodular")
        self.assertEqual(evidence["detail"], "u_nonintegral")
        self.assertFalse(evidence["u_integral"])
        self.assertEqual(evidence["u"][0], ["1/2", "0", "0"])

    def test_wrong_determinant_is_nonunimodular(self):
        row = machine_row(basis=[[2, 0, 0], [0, 1, 0], [0, 0, 1]])
        evidence = census.derive_record(row, oracle_row(), 1)
        self.assertEqual(evidence["status"], "nonunimodular")
        self.assertEqual(evidence["detail"], "u_determinant")
        self.assertTrue(evidence["u_integral"])
        self.assertEqual(evidence["u_det"], "2")

    def test_query_accepts_a_real_shear_and_keeps_origin_evidence(self):
        rows = [machine_row()]
        status, results, error = census.classify_query(
            rows, [oracle_row(origin=["0", "0", "1"])], 1
        )
        self.assertEqual(status, "nonempty")
        self.assertIsNone(error)
        self.assertEqual(results[0]["status"], "origin_mismatch")
        self.assertEqual(results[0]["origin_delta"], ["0", "0", "1"])


class FailClosedTests(unittest.TestCase):
    def test_duplicate_machine_label_is_rejected(self):
        rows = [machine_row(label="P1"), machine_row(label="P1", child=5)]
        status, results, error = census.classify_query(rows, [oracle_row()], 1)
        self.assertEqual(status, "error")
        self.assertEqual(error["kind"], "duplicate_machine_label")
        self.assertTrue(all(result["status"] == "duplicate_machine_label" for result in results))

    def test_duplicate_oracle_rows_are_rejected(self):
        rows = [machine_row(label="P1"), machine_row(label="C1")]
        oracles = [oracle_row(label="P1"), oracle_row(label="P1")]
        status, _, error = census.classify_query(rows, oracles, 1)
        self.assertEqual(status, "error")
        self.assertEqual(error["kind"], "duplicate_oracle_label")

    def test_missing_and_extra_rows_are_rejected(self):
        rows = [machine_row(label="P1"), machine_row(label="C1")]
        status, _, error = census.classify_query(rows, [oracle_row(label="P1")], 1)
        self.assertEqual((status, error["kind"]), ("error", "row_count_mismatch"))
        oracles = [oracle_row(label="P1"), oracle_row(label="X1")]
        status, _, error = census.classify_query(rows, oracles, 1)
        self.assertEqual((status, error["kind"]), ("error", "unexpected_oracle_row"))
        oracles = [oracle_row(label="P1"), oracle_row(label="C1", child=7, size=4)]
        status, _, error = census.classify_query(rows, oracles, 1)
        self.assertEqual((status, error["kind"]), ("error", "child_mismatch"))
        oracles = [oracle_row(label="P1"), oracle_row(label="C1", size=4)]
        status, _, error = census.classify_query(rows, oracles, 1)
        self.assertEqual((status, error["kind"]), ("error", "size_mismatch"))

    def test_malformed_rows_and_stray_text_are_syntax_errors(self):
        bad_rows = [
            BANNER + "\n" + HEADER + "\n1 P1 1 P1 (1,0,0),(0,1,0) (0,0,0)\n*",
            BANNER + "\n" + HEADER + "\n1 P1 1 P1 (1,0,0),(0,1,0),(0,0,1) (0,0,0)\nstray\n*",
            BANNER + "\n" + HEADER + "\n1 P1 x P1 (1,0,0),(0,1,0),(0,0,1) (0,0,0)\n*",
            "Isotropy, Version 9.5.0, Jan 2020\n" + "\n".join(BANNER.splitlines()[1:]),
            BANNER,
            BANNER + "\n" + HEADER + "\n*",
        ]
        for text in bad_rows:
            with self.subTest(text=text[-60:]):
                with self.assertRaises(census.OracleParseError):
                    census.parse_oracle_table(text)

    def test_empty_official_table_is_recorded_separately(self):
        self.assertEqual(census.parse_oracle_table(EMPTY), [])
        rows = [machine_row(), machine_row(label="C1", child=5)]
        status, results, error = census.classify_query(rows, [], 1)
        self.assertEqual(status, "empty")
        self.assertIsNone(error)
        self.assertTrue(all(result["status"] == census.UNRETURNED for result in results))

    def test_summary_partitions_sum_exactly(self):
        records = []
        outcomes = {}
        records.append(dict(census.base_record(machine_row(), 1),
                            status="candidate", detail="d", u=None, u_det="1",
                            u_integral=True, origin_exact=True, origin_delta=None,
                            origin_converted=None, origin_printed=None,
                            basis_converted=None, basis_printed=None, oracle_row=None))
        records.append(dict(census.base_record(machine_row(label="C1", child=5), 1),
                            status=census.UNRETURNED, detail="d", u=None, u_det=None,
                            u_integral=None, origin_exact=None, origin_delta=None,
                            origin_converted=None, origin_printed=None,
                            basis_converted=None, basis_printed=None, oracle_row=None))
        outcomes[(1, 1)] = "nonempty"
        outcomes[(1, 2)] = "empty"
        summary = census.build_summary(
            [{"parent": 1, "index": 1, "irrep": "GM1"}, {"parent": 1, "index": 2, "irrep": "GM2"}],
            records, outcomes, [], 2,
        )
        self.assertTrue(all(summary["accounting"].values()))
        self.assertTrue(summary["collector_completed"])
        self.assertEqual(summary["classification"]["candidate"], 1)
        self.assertEqual(summary["records_unreturned_empty"], 1)
        self.assertEqual(summary["records_total"], 2)
        inconsistent = census.build_summary(
            [{"parent": 1, "index": 1, "irrep": "GM1"}], records, outcomes, [], 5
        )
        self.assertFalse(inconsistent["collector_completed"])
        self.assertFalse(inconsistent["accounting"]["every_machine_record_emitted"])
        self.assertEqual(inconsistent["errors"][0]["kind"], "internal_accounting_error")

    def test_pinned_setting_commands_are_fixed(self):
        self.assertEqual(census.ORACLE_COMMANDS.count(census.SETTING_COMMAND), 1)
        self.assertEqual(census.ORACLE_COMMANDS[0], "PAGE 1000")
        self.assertEqual(census.ORACLE_COMMANDS[1], "SC 250")
        self.assertLess(census.ORACLE_COMMANDS.index("PAGE 1000"),
                        census.ORACLE_COMMANDS.index(census.SETTING_COMMAND))
        self.assertLess(census.ORACLE_COMMANDS.index(census.SETTING_COMMAND),
                        census.ORACLE_COMMANDS.index("VALUE PARENT {parent}"))
        self.assertEqual(census.ORACLE_COMMANDS[-2:], ("DISPLAY ISOTROPY", "QUIT"))
        self.assertEqual(census.parse_parents(["3,43", "43"]), [3, 43])
        with self.assertRaises(census.CollectorError):
            census.parse_parents(["0"])
        with self.assertRaises(census.CollectorError):
            census.parse_parents(["231"])


class RealArchiveTests(unittest.TestCase):
    def test_ordinals_match_the_rust_record_index(self):
        records = census.load_machine_records(census.load_inventory())
        self.assertEqual([r["ordinal"] for r in records], list(range(15239)))
        row = records[13345]
        self.assertEqual(
            (row["ordinal"], row["parent"], row["irrep"], row["direction_label"], row["child"]),
            (13345, 225, "GM4-", "C2", 8),
        )

    def setUp(self):
        if not os.path.isfile(census.ARCHIVE):
            self.skipTest("pinned archive is not present")
        self.inventory = census.load_inventory()
        self.records = census.load_machine_records(self.inventory)

    def test_real_sg43_l1_rows_classify_as_candidates(self):
        query = next(entry for entry in self.inventory
                     if entry["parent"] == 43 and entry["irrep"] == "L1")
        rows = [record for record in self.records
                if record["parent"] == 43
                and record["source_irrep_index"] == query["index"]]
        self.assertEqual(len(rows), 9)
        status, results, error = census.classify_query(
            rows, census.parse_oracle_table(SG43_L1), 43
        )
        self.assertEqual(status, "nonempty")
        self.assertIsNone(error)
        self.assertTrue(all(result["status"] == "candidate" for result in results))
        by_label = {row["direction_label"]: result for row, result in zip(rows, results)}
        self.assertEqual(by_label["P1"]["u"], IDENTITY)
        # A real non-identity unimodular change of basis is accepted.
        self.assertEqual(by_label["P3"]["u"], [["0", "0", "1"], ["1", "0", "0"], ["0", "1", "0"]])
        self.assertEqual(by_label["P4"]["origin_printed"], ["7/8", "7/16", "7/16"])
        self.assertTrue(by_label["P4"]["origin_exact"])

    def test_machine_inventory_matches_pinned_counts(self):
        self.assertEqual(len(self.inventory), census.EXPECTED_IRREPS)
        self.assertEqual(len(self.records), census.EXPECTED_RECORDS)
        keys = {(record["parent"], record["source_irrep_index"],
                 record["direction_label"], record["child"]) for record in self.records}
        self.assertEqual(len(keys), len(self.records))


class FakeOracleRunnerTests(unittest.TestCase):
    def test_isolated_runner_uses_private_cwd_and_fixed_setting(self):
        with tempfile.TemporaryDirectory() as root:
            capture = os.path.join(root, "stdin.txt")
            binary = os.path.join(root, "fake_iso")
            with open(binary, "w", encoding="utf-8") as handle:
                handle.write(
                    "#!/usr/bin/env python3\n"
                    "import os, sys\n"
                    "data = sys.stdin.read()\n"
                    "open(os.environ['FAKE_CAPTURE'], 'w').write(data)\n"
                    f"sys.stdout.write({SG43_L1!r})\n"
                )
            os.chmod(binary, os.stat(binary).st_mode | stat.S_IXUSR)
            work_root = os.path.join(root, "work")
            os.makedirs(work_root)
            with mock.patch.dict(os.environ, {"FAKE_CAPTURE": capture}):
                outcome = census.run_query(
                    binary, root, work_root, {"parent": 43, "irrep": "L1"}
                )
            self.assertIsNone(outcome["error"])
            self.assertEqual(len(outcome["rows"]), 9)
            with open(capture, encoding="utf-8") as handle:
                sent = handle.read()
            self.assertIn(census.SETTING_COMMAND + "\n", sent)
            self.assertIn("VALUE PARENT 43\n", sent)
            self.assertIn("VALUE IRREP L1\n", sent)
            self.assertFalse(os.listdir(work_root))  # per-query cwd was removed

    def test_failed_query_still_emits_every_record(self):
        with tempfile.TemporaryDirectory() as root:
            binary = os.path.join(root, "bad_iso")
            with open(binary, "w", encoding="utf-8") as handle:
                handle.write("#!/usr/bin/env python3\nimport sys\nsys.exit(3)\n")
            os.chmod(binary, os.stat(binary).st_mode | stat.S_IXUSR)
            work_root = os.path.join(root, "work")
            os.makedirs(work_root)
            queries = [{"parent": 1, "index": 1, "irrep": "GM1"}]
            rows_by_query = {(1, 1): [machine_row()]}
            records, outcomes, errors = census.collect(
                queries, rows_by_query, binary, root, work_root, 1
            )
            self.assertEqual(len(records), 1)  # the record is preserved, not dropped
            self.assertEqual(outcomes[(1, 1)], "error")
            self.assertEqual(records[0]["status"], "oracle_error")
            self.assertEqual(errors[0]["kind"], "nonzero_exit")
            summary = census.build_summary(queries, records, outcomes, errors, 1)
            self.assertFalse(summary["collector_completed"])
            self.assertEqual(summary["records_failed"], 1)
            self.assertTrue(summary["accounting"]["records_partition"])


if __name__ == "__main__":
    unittest.main()
