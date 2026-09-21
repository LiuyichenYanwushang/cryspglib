"""The source-equivalence gate must see complex lattice-translation phases."""

import unittest
import zipfile

import generate_subduction_compound_fixtures as fixtures


class ConjugateCharacterGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with zipfile.ZipFile(fixtures.ROOT / "isotropy_subgroup/CIR_data.zip") as archive:
            cls.lines = archive.read("CIR_data.txt").decode("ascii").splitlines()

    def check_source(self, irnumber):
        for index, line in enumerate(self.lines):
            if '"' not in line:
                continue
            header = fixtures.source._parse_cir_header(line, index + 1)
            if header["irnumber"] != irnumber:
                continue
            arms, next_index = fixtures.source._read_exact_cir_block(
                self.lines, index + 1, 16 * header["kcount"], "test source",
                fixtures.source._parse_cir_integer,
            )
            return fixtures.check_self_conjugate_characters(header, arms, next_index, self.lines)
        self.fail(f"source {irnumber} missing")

    def test_quaternionic_source_is_self_conjugate(self):
        self.assertEqual(self.check_source(559), 8)  # SG19 R1

    def test_disjoint_stars_fail_even_when_unshifted_characters_are_real(self):
        # SG23 W1 has real characters at the four archived representatives,
        # but its I-centring phase is -i. Checking those four traces alone
        # would falsely identify the seed with its conjugate.
        with self.assertRaisesRegex(ValueError, "non-self-conjugate.*716"):
            self.check_source(716)


if __name__ == "__main__":
    unittest.main()
