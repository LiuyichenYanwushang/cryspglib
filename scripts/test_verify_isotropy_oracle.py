"""Offline regressions for the isotropy oracle gate.

Run with: python3 -m unittest discover -s scripts -p test_verify_isotropy_oracle.py
The fixed rows below are from iso 9.6.1, SG 221 GM4+ P1, SET I ALL OR 1.
Only the archive input and process boundary are mocked; the real parsers,
frame conversion, lattice checks and gate exit status are exercised.
"""

import contextlib
import io
import unittest
from types import SimpleNamespace
from unittest import mock

import verify_isotropy_oracle as oracle
from direction_map import direction_str


GEOMETRY = """***Subgroup Size Dir Basis Vectors Origin
83 P4/m 1 P1 (0,0,1),(0,-1,0),(1,0,0) (0,0,0)
*
"""
VECTORS = """***Subgroup Dir
83 P4/m P1 (a,0,0)
*
"""


class IsotropyOracleTests(unittest.TestCase):
    def setUp(self):
        self.records = {(221, "GM4+"): [{
            "subgroup": 83,
            "label": "P1",
            "direction": "(a,0,0)",
            "dim": 3,
            "size": 1,
            "basis": [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
            "origin": [0, 0, 0],
        }]}

    def run_gate(self, geometry=GEOMETRY, vectors=VECTORS):
        responses = [
            SimpleNamespace(returncode=0, stdout=output, stderr="")
            for output in (geometry, vectors)
        ]
        output = io.StringIO()
        with mock.patch.object(oracle, "CASES", [(221, "GM4+")]), \
                mock.patch.object(oracle, "machine_records", return_value=self.records), \
                mock.patch.object(oracle.os.path, "isfile", return_value=True), \
                mock.patch.object(oracle.subprocess, "run", side_effect=responses) as run, \
                contextlib.redirect_stdout(output):
            status = oracle.main()
        for call in run.call_args_list:
            self.assertIn("SET I ALL OR 1\n", call[1]["input"])
        self.assertIn("SHOW DIRECTION VECTOR\n", run.call_args_list[1][1]["input"])
        return status, output.getvalue()

    def test_exact_reference_passes(self):
        status, output = self.run_gate()
        self.assertEqual(status, 0, output)
        self.assertIn("descriptor strings checked: 1", output)
        self.assertIn("origin on 1: all exact", output)

    def test_origin_shift_fails_even_by_a_parent_lattice_vector(self):
        for origin in ("(1,0,0)", "(1/2,0,0)"):
            with self.subTest(origin=origin):
                status, output = self.run_gate(geometry=GEOMETRY.replace("(0,0,0)", origin))
                self.assertEqual(status, 1, output)
                self.assertIn("origin", output)
                self.assertIn("exactly", output)

    def test_wrong_descriptor_fails_with_correct_label_and_geometry(self):
        status, output = self.run_gate(vectors=VECTORS.replace("(a,0,0)", "(0,a,0)"))
        self.assertEqual(status, 1, output)
        self.assertIn("descriptor", output)

    def test_missing_or_unexpected_vector_labels_fail(self):
        for vectors in ("", VECTORS.replace("P1", "P2"), VECTORS + "83 P4/m P2 (a,a,0)\n"):
            with self.subTest(vectors=vectors):
                status, output = self.run_gate(vectors=vectors)
                self.assertEqual(status, 1, output)
                self.assertIn("direction vector labels", output)

    def test_duplicate_vector_labels_fail_instead_of_overwriting(self):
        with self.assertRaisesRegex(RuntimeError, "duplicate.*P1"):
            self.run_gate(vectors=VECTORS + "83 P4/m P1 (a,0,0)\n")

    def test_descriptor_separator_and_whitespace_follow_rust_matching(self):
        status, output = self.run_gate(vectors=VECTORS.replace("(a,0,0)", "( a ;\t0 ; 0 )"))
        self.assertEqual(status, 0, output)
        self.assertIn("descriptor strings checked: 1", output)

    def test_equal_volume_does_not_mask_a_different_lattice(self):
        geometry = GEOMETRY.replace(
            "(0,0,1),(0,-1,0),(1,0,0)", "(2,0,0),(0,1/2,0),(0,0,1)"
        )
        status, output = self.run_gate(geometry=geometry)
        self.assertEqual(status, 1, output)
        self.assertIn("different lattices", output)

    def test_oracle_process_failure_is_not_an_empty_success(self):
        response = SimpleNamespace(returncode=1, stdout="", stderr="failed")
        for runner in (oracle.run_oracle, oracle.run_oracle_direction_vectors):
            with self.subTest(runner=runner.__name__), \
                    mock.patch.object(oracle.os.path, "isfile", return_value=True), \
                    mock.patch.object(oracle.subprocess, "run", return_value=response), \
                    self.assertRaisesRegex(RuntimeError, "iso exited with 1"):
                runner(221, "GM4+")

    def test_c2_mapping_retains_parent_context(self):
        # SHOW DIRECTION VECTOR: 177 L1 C2 versus 225 GM4- C2.  The
        # component order differs even after accepting ';' as ','.
        for sg, expected in ((177, "(a;b;a)"), (225, "(a,a,b)"), (177, "(a;b;a)")):
            with self.subTest(sg=sg):
                self.assertEqual(direction_str(3, 2, "C2", sg), expected)


if __name__ == "__main__":
    unittest.main()
