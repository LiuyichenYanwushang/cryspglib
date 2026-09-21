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
    return {
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


if __name__ == "__main__":
    unittest.main()
