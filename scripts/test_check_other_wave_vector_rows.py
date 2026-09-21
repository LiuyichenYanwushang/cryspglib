"""Offline regressions for the other-wave-vector row gate.

Run with: python3 -m unittest discover -s scripts -p test_check_other_wave_vector_rows.py

The synthetic archive below mirrors the pinned layout at a size the tests can
mutate: two sources in one space group, three condensate records, two of them
carrying rows.  Every gate is driven by a mutation of that baseline, so a gate
that stops firing fails here instead of passing silently.
"""

import copy
import sys
import unittest

import check_other_wave_vector_rows as gate

EXPECTED = gate.Expectations(sources=2, records=2, rows=3, condensates=3)


def baseline():
    """A synthetic archive: two sources, three condensates, two rows.

    The little table carries one line (free direction) and one point (no free
    direction) so a source pinned to a fixed k trips the parameterization gate.
    """
    little_k = []
    little_labels = []
    for lattice in range(gate.LATTICES):
        for slot in range(gate.K_SLOTS):
            if lattice == 0 and slot == 1:
                groups = [(0, 0, 0, 1), (1, 0, 1, 1), (0, 0, 0, 1), (0, 0, 0, 1)]
            elif lattice == 0 and slot == 2:
                groups = [(0, 0, 0, 1), (1, 1, 2, 1), (0, 0, 0, 1), (0, 0, 0, 1)]
            elif lattice == 0 and slot == 3:
                groups = [(1, 0, 0, 2), (0, 0, 0, 1), (0, 0, 0, 1), (0, 0, 0, 1)]
            else:
                groups = [(0, 0, 0, 1)] * gate.K_GROUPS
            little_k.extend(entry for group in groups for entry in group)
            little_labels.append(
                {1: "DT", 2: "SM", 3: "X"}.get(slot if lattice == 0 else -1, "P")
            )
    return {
        "little_k": little_k,
        "little_k_count": [4] + [0] * (gate.LATTICES - 1),
        "little_k_label": little_labels,
        "little_label": ["DT1", "SM1"],
        "little_sg": [196, 196],
        "little_k_index": [2, 3],
        "little_dim": [6, 12],
        "sg_lattice": [1] * 230,
        "irrep_sections": sorted(gate.W_SECTIONS),
        "labels": ["DT1", "SM1"],
        "source_sg": [196, 196],
        "dimension": [6, 12],
        "type": [1, 1],
        "count": [0, 2, 1],
        "pointer": [0, 1, 3],
        "irrep": [1, 2, 2],
        "frequency": [1, 3, 2],
        "parent": [196, 196, 196],
    }


class OtherWaveVectorTests(unittest.TestCase):
    def check(self, data):
        return gate.check(data, EXPECTED)

    def test_the_baseline_passes_every_gate(self):
        self.assertEqual(self.check(baseline()), [])

    def test_a_source_on_a_fixed_wave_vector_fails(self):
        # The second source sits on the point (1,0,0)/2, so its rows have a
        # numeric k: the gate must refuse to keep reporting them as unresolved.
        data = baseline()
        data["little_k_index"] = [4, 3]
        failures = self.check(data)
        self.assertTrue(any("with a fixed k" in message for message in failures))

    def test_a_source_missing_from_the_little_table_fails(self):
        data = baseline()
        data["little_label"] = ["DT1", "SM9"]
        failures = self.check(data)
        self.assertTrue(any("little-table records" in message for message in failures))

    def test_a_source_naming_an_out_of_range_k_slot_fails(self):
        data = baseline()
        data["little_k_index"] = [2, 9]
        failures = self.check(data)
        self.assertTrue(any("outside its lattice" in message for message in failures))

    def test_a_little_k_array_of_the_wrong_length_fails(self):
        data = baseline()
        data["little_k"] = data["little_k"][:-16]
        failures = self.check(data)
        self.assertTrue(any("little_k has" in message for message in failures))

    def test_a_frequency_above_its_dimension_fails(self):
        data = baseline()
        data["frequency"][1] = 13  # SM1 has dimension 12
        self.assertIn("frequencies exceed", self.check(data)[0])

    def test_a_zero_frequency_fails(self):
        data = baseline()
        data["frequency"][0] = 0
        failures = self.check(data)
        self.assertTrue(any("zero frequency" in message for message in failures))

    def test_a_record_paired_with_another_space_group_fails(self):
        data = baseline()
        data["parent"][1] = 225
        failures = self.check(data)
        self.assertTrue(any("another space group" in message for message in failures))

    def test_an_out_of_range_irrep_index_fails(self):
        data = baseline()
        data["irrep"][2] = 3
        failures = self.check(data)
        self.assertTrue(any("out of range" in message for message in failures))
        # The out-of-range entry must not be dereferenced as a source.
        self.assertFalse(any("another space group" in message for message in failures))

    def test_counts_that_do_not_sum_to_the_packed_table_fail(self):
        data = baseline()
        data["count"][1] = 1
        failures = self.check(data)
        self.assertTrue(any("do not match the counts' sum" in message for message in failures))

    def test_a_pointer_index_that_disagrees_with_the_counts_fails(self):
        data = baseline()
        data["pointer"] = [0, 1, 0]
        failures = self.check(data)
        self.assertTrue(any("pointer starts" in message for message in failures))

    def test_a_new_other_wave_vector_section_fails(self):
        data = baseline()
        data["irrep_sections"] = sorted(gate.W_SECTIONS | {"irrep_w_k_points"})
        failures = self.check(data)
        self.assertTrue(any("expected exactly" in message for message in failures))

    def test_a_short_source_table_fails(self):
        data = baseline()
        data["dimension"] = [6]
        failures = self.check(data)
        self.assertTrue(any("not all 2 long" in message for message in failures))

    def test_a_count_array_that_is_not_per_record_fails(self):
        data = baseline()
        data["count"] = [0, 2]
        failures = self.check(data)
        self.assertTrue(any("isotropy_w_subduce_count has 2 entries" in message
                            for message in failures))

    def test_every_mutation_is_detected_on_a_copy(self):
        # Guard against a gate that mutates its input: the baseline must stay
        # usable and still pass after a failing run.
        data = baseline()
        broken = copy.deepcopy(data)
        broken["frequency"][0] = 7  # DT1 has dimension 6
        self.assertNotEqual(self.check(broken), [])
        self.assertEqual(self.check(data), [])


    def test_every_row_has_a_frozen_character_table(self):
        """No source is left open, so no row is blocked."""
        self.assertEqual(gate.UNRESOLVED_SOURCES, frozenset())
        data = baseline()
        frozen, blocked, sources = gate.frozen_coverage(data)
        self.assertEqual((frozen, blocked), (3, 0))
        self.assertEqual(sources, [])
        # The counter still distinguishes a source that has no table: the
        # pinned set is the only thing that decides, so an artificial entry
        # moves every row into the blocked bucket.
        gate.UNRESOLVED_SOURCES = frozenset({(202, "DT3")})
        try:
            data["source_sg"] = [202, 202]
            data["labels"] = ["DT3", "DT3"]
            frozen, blocked, sources = gate.frozen_coverage(data)
        finally:
            gate.UNRESOLVED_SOURCES = frozenset()
        self.assertEqual((frozen, blocked), (0, 3))
        self.assertEqual(sources, [(202, "DT3")])


if __name__ == "__main__":
    unittest.main()
