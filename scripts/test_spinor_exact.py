"""Focused tests for the exact, generation-only spin source sidecar."""

import hashlib
import os
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction

sys.path.insert(0, os.path.dirname(__file__))
import spinor_exact as exact


class ExactSpinSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources = exact.parse_all_exact()

    def test_codebook_and_zeta_are_exact(self):
        self.assertEqual(len(exact.ZETA24_POWERS), 24)
        self.assertEqual(exact.ZETA24**24, exact.Complex24(exact.ONE, exact.ZERO))
        self.assertTrue(all(exact.ZETA24**n != exact.Complex24(exact.ONE, exact.ZERO)
                            for n in range(1, 24)))
        self.assertEqual(exact.Radical24(Fraction(1)) * exact.Radical24(b=Fraction(1)),
                         exact.Radical24(b=Fraction(1)))
        self.assertEqual(exact.Radical24(b=Fraction(1)) * exact.Radical24(b=Fraction(1)),
                         exact.Radical24(Fraction(2)))
        self.assertEqual(exact.Radical24(c=Fraction(1)) * exact.Radical24(c=Fraction(1)),
                         exact.Radical24(Fraction(3)))

    def test_source_census_and_all_golden_hashes(self):
        counts, dimensions, _domains, su2_pairs, character_pairs = exact.source_spelling_census(self.sources)
        self.assertEqual(dict(counts), {
            "files": 230, "operations": 2609, "kblocks": 1350, "rows": 3611,
            "real_rows": 1470, "polar_rows": 2141, "character_columns": 31491,
            "operation_indices": 9945, "su2_pairs": 10436, "character_pairs": 14702,
        })
        self.assertEqual(dimensions, Counter({"1": 1895, "2": 1609, "4": 93, "3": 14}))
        self.assertEqual(len(su2_pairs), 24)
        self.assertEqual(len(character_pairs), 39)
        hashes = exact.source_spelling_hashes(self.sources)
        self.assertEqual(hashes["combined"], exact.COMBINED_GOLDEN)
        self.assertEqual(hashes["su2_pair"], exact.SU2_PAIR_GOLDEN)
        self.assertEqual(hashes["character_pair"], exact.CHARACTER_PAIR_GOLDEN)
        for domain, golden in exact.DOMAIN_GOLDENS.items():
            self.assertEqual(hashes[domain], golden, domain)

    def test_exact_structure_and_lattice_products(self):
        self.assertEqual(exact.validate_exact_sources(self.sources), (52381, 32546, 19835))
        self.assertEqual(
            Counter(len(exact._translation_lattice_cosets(source.operations)) for source in self.sources),
            Counter({1: 209, 2: 15, 4: 6}),
        )

    def test_identity_is_searched_and_character_orthogonality_is_exact(self):
        identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        for source in self.sources:
            groups = defaultdict(list)
            for row in source.rows:
                groups[(row.raw_k, row.operation_indices)].append(row)
            for rows in groups.values():
                identity_columns = [
                    i for i in rows[0].operation_indices
                    if source.operations[i].rotation == identity
                    and source.operations[i].translation == (Fraction(0),) * 3
                ]
                self.assertEqual(len(identity_columns), 1)
                identity_column = identity_columns[0]
                for row in rows:
                    self.assertEqual(row.characters[identity_column].re,
                                     exact.Radical24(Fraction(row.dimension)))
                    self.assertEqual(row.characters[identity_column].im, exact.ZERO)
                for left in rows:
                    for right in rows:
                        inner = sum(
                            (a * b.conjugate() for a, b in zip(left.characters, right.characters)),
                            exact.Complex24(),
                        )
                        expected = len(left.operation_indices) if left == right else 0
                        self.assertEqual(inner, exact.Complex24(exact.Radical24(Fraction(expected)), exact.ZERO))

    def _parse_rejecting_fixture(self, mutate, row="1.0"):
        operation = "1 0 0 0 1 0 0 0 1 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0"
        text = "\n".join([
            "SG=3", " name=P2", " nsym= 1", " spinor=True", "symmetries=",
            mutate(operation), "", " kpoint  GM : 0.0 0.0 0.0  : 1",
            f"-GM1 1    {row}", "",
        ])
        with tempfile.NamedTemporaryFile(mode="w", suffix="-spin.dat") as stream:
            stream.write(text)
            stream.flush()
            with self.assertRaises(exact.ExactSpinSourceError):
                exact.parse_spinor_file_exact(stream.name)

    def test_unknown_alternate_signed_zero_nonfinite_and_operation_shapes_reject(self):
        for bad in ("+0.0", "0.000000", "-0.0", "nan", "inf"):
            self._parse_rejecting_fixture(lambda operation, bad=bad: operation.replace("0.0", bad, 1))
        self._parse_rejecting_fixture(lambda operation: " ".join(operation.split()[:-1]))
        self._parse_rejecting_fixture(lambda operation: operation + " 0.0")

    def test_illegal_pair_and_operation_index_reject(self):
        self._parse_rejecting_fixture(lambda operation: " ".join(
            operation.split()[:16] + ["-1.0"] + operation.split()[17:]
        ))
        operation = "1 0 0 0 1 0 0 0 1 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0"
        text = "\n".join([
            "SG=3", " nsym= 1", " spinor=True", "symmetries=", operation,
            "", " kpoint GM : 0.0 0.0 0.0  : 0", "-GM1 1 1.0", "",
        ])
        with tempfile.NamedTemporaryFile(mode="w", suffix="-spin.dat") as stream:
            stream.write(text); stream.flush()
            with self.assertRaises(exact.ExactSpinSourceError):
                exact.parse_spinor_file_exact(stream.name)

    def test_pinned_small_direct_zero_is_position_sensitive(self):
        self.assertEqual(exact.DIRECT_CHAR_TOKENS["1e-05"], exact.ZERO)
        self._parse_rejecting_fixture(lambda operation: operation, row="1e-05")
        with self.assertRaises(exact.ExactSpinSourceError):
            exact._parse_direct_character("1e-05", "fixture", 1, 1, 3, 0, 0)

    def test_source_ordinal_linkage_is_frozen(self):
        source = self.sources[2]
        broken = replace(source, rows=(replace(source.rows[0], source_row_ordinal=1),) + source.rows[1:])
        with self.assertRaisesRegex(exact.ExactSpinSourceError, "ordinals"):
            exact._validate_source_row_ordinals(broken)


if __name__ == "__main__":
    unittest.main()
