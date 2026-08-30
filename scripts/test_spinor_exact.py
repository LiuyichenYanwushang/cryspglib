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


class ExactSpinMaterializationTests(unittest.TestCase):
    def test_terminal_materialization_and_pauli_components(self):
        radical = exact.Radical24(
            Fraction(1, 3), Fraction(1, 2), Fraction(-1, 4), Fraction(1, 6)
        )
        self.assertEqual(
            radical.materialize(),
            float(Fraction(1, 3))
            + float(Fraction(1, 2)) * 2.0**0.5
            - float(Fraction(1, 4)) * 3.0**0.5
            + float(Fraction(1, 6)) * 6.0**0.5,
        )
        value = exact.Complex24(radical, -radical)
        self.assertEqual(
            value.materialize(), (radical.materialize(), -radical.materialize())
        )

        u0 = exact.Radical24(Fraction(1, 2))
        u1 = exact.Radical24(b=Fraction(1, 2))
        u2 = exact.Radical24(c=Fraction(1, 2))
        u3 = exact.Radical24(d=Fraction(1, 2))
        operation = exact.ExactSpinOperation(
            (),
            (),
            (
                exact.Complex24(u0, u3),
                exact.Complex24(u2, u1),
                exact.Complex24(-u2, u1),
                exact.Complex24(u0, -u3),
            ),
        )
        self.assertEqual(operation.pauli_components(), (u0, u1, u2, u3))


class ExactSpinSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources = exact.parse_all_exact()
        cls.validation = exact.validate_exact_sources(cls.sources)

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
        self.assertEqual(self.validation, (52381, 32546, 19835))
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
                    local for local, global_index in enumerate(rows[0].operation_indices)
                    if source.operations[global_index].rotation == identity
                    and source.operations[global_index].translation == (Fraction(0),) * 3
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

    def test_identity_global_index_is_mapped_to_reordered_local_column(self):
        source = next(source for source in self.sources if any(
            len(row.operation_indices) > 1 for row in source.rows
        ))
        row = next(row for row in source.rows if len(row.operation_indices) > 1)
        identity_global = next(
            index for index in row.operation_indices
            if source.operations[index].rotation == (1, 0, 0, 0, 1, 0, 0, 0, 1)
            and source.operations[index].translation == (Fraction(0),) * 3
        )
        reordered = replace(
            row,
            operation_indices=tuple(reversed(row.operation_indices)),
            characters=tuple(reversed(row.characters)),
        )
        self.assertNotEqual(identity_global, exact._identity_local_column(source, reordered))
        self.assertEqual(
            exact._identity_local_column(source, reordered),
            len(reordered.operation_indices) - 1 - row.operation_indices.index(identity_global),
        )

    def test_production_parse_is_gated_by_exact_validator(self):
        import parse_spinor_data
        original = exact.parse_all_exact
        try:
            exact.parse_all_exact = lambda *_args, **_kwargs: (_ for _ in ()).throw(
                exact.ExactSpinSourceError("synthetic exact gate failure")
            )
            with self.assertRaisesRegex(exact.ExactSpinSourceError, "synthetic exact gate"):
                parse_spinor_data.parse_all_spinor()
        finally:
            exact.parse_all_exact = original

    def _parse_rejecting_fixture(self, mutate, row="1.0"):
        operation = "1 0 0 0 1 0 0 0 1 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0"
        text = "\n".join([
            "SG=3", " name=P2", " nsym= 1", " spinor=True", "symmetries=",
            mutate(operation), "", " kpoint  GM : 0.0 0.0 0.0  : 1",
            f"-GM1 1    {row}", "",
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "irreps-SG=3-spin.dat")
            with open(path, "w") as stream:
                stream.write(text)
            with self.assertRaises(exact.ExactSpinSourceError):
                exact.parse_spinor_file_exact(path)

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
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "irreps-SG=3-spin.dat")
            with open(path, "w") as stream:
                stream.write(text)
            with self.assertRaises(exact.ExactSpinSourceError):
                exact.parse_spinor_file_exact(path)

    def test_headers_k_blocks_and_labels_are_strict(self):
        operation = "1 0 0 0 1 0 0 0 1 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0"

        def parse(lines):
            with tempfile.TemporaryDirectory() as directory:
                path = os.path.join(directory, "irreps-SG=3-spin.dat")
                with open(path, "w") as stream:
                    stream.write("\n".join(lines))
                exact.parse_spinor_file_exact(path)

        valid = [
            "SG=3", " name=P2", " nsym= 1", " spinor=True", "symmetries=",
            operation, "", " kpoint GM : 0.0 0.0 0.0  : 1", "-GM1 1 1.0", "",
        ]
        for mutation in (
            lambda lines: ["SG=03"] + lines[1:],
            lambda lines: lines[:2] + [" nsym= +1"] + lines[3:],
            lambda lines: lines[:3] + [" spinor=False"] + lines[4:],
            lambda lines: lines[:4] + ["SG=3"] + lines[4:],
            lambda lines: lines[:6] + ["unknown=1"] + lines[6:],
            lambda lines: lines[:7] + [lines[7].replace("GM", "WHAT", 1)] + lines[8:],
            lambda lines: lines[:7] + [" kpoint GM : 0.0 0.0 0.0  : 1"] + lines[7:],
            lambda lines: lines[:7] + ["-X1 1 1.0"] + lines[8:],
        ):
            with self.assertRaises(exact.ExactSpinSourceError):
                parse(mutation(valid))

        duplicate_block = valid[:-1] + [" kpoint GM : 0.0 0.0 0.0  : 1", "-GM1 1 1.0", ""]
        with self.assertRaises(exact.ExactSpinSourceError):
            parse(duplicate_block)

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
