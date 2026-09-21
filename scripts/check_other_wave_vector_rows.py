#!/usr/bin/env python3
"""What the pinned archive does and does not pin about the other-wave-vector rows.

`data_isotropy.txt` carries the `isotropy_w_subduce_*` arrays: 5756 entries over
1006 isotropy records, naming parent irreps **at other wave vectors** -- the
`DT`/`SM` star labels of the eleven cubic space groups 196, 202, 203, 209, 210,
216, 219, 225, 226, 227 and 228 (73 distinct labels).

Task 9's audit reports those rows separately instead of computing them, and this
script makes that decision executable:

* the 73 irreps they name are stored in `data_irreps.txt` as `irrep_w_label`,
  `irrep_w_space_group`, `irrep_w_dimension` and `irrep_w_type` only -- there is
  **no** other-wave-vector k-vector table and **no** character/matrix table to
  compute a subduction from, so no engine can produce their frequencies from the
  pinned archive;
* everything the pinned archive *does* pin is checked here: array lengths, the
  per-record counts against the packed table, the sparse pointer index against
  the record start offsets, in-range irrep indices, the record's parent space
  group against the source irrep's own space group, non-zero frequencies, and
  the bound `frequency <= dimension` (a trivial-representation multiplicity can
  never exceed the dimension of the representation containing it).

Exit codes: 0 = every check passed, 1 = a check failed, 2 = pinned input missing.
"""

import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import generate_irrep_data as generator  # noqa: E402

# The only other-wave-vector sections the pinned irrep file may carry.  A new
# section (a k vector or a character table) would mean the rows became
# computable, which must fail this gate loudly instead of passing unnoticed.
W_SECTIONS = {
    "irrep_w_label",
    "irrep_w_space_group",
    "irrep_w_dimension",
    "irrep_w_type",
}


class Expectations:
    """Row counts the pinned archive is expected to hold.

    The offline tests drive `check` with a small synthetic archive, so the
    counts are an argument instead of module constants baked into the gates.
    """

    def __init__(self, sources, records, rows, condensates):
        self.sources = sources
        self.records = records
        self.rows = rows
        self.condensates = condensates


PINNED = Expectations(sources=73, records=1006, rows=5756, condensates=15239)


def parse_pinned():
    """Read both pinned data files out of the SHA256-verified archive."""
    irrep_lines = generator.read_file("data_irreps.txt")
    isotropy_lines = generator.read_file("data_isotropy.txt")
    irrep_sec = generator.get_sections(irrep_lines)
    isotropy_sec = generator.get_sections(isotropy_lines)
    return {
        "irrep_sections": sorted(
            name for name in irrep_sec if name.startswith("irrep_w_")
        ),
        "labels": generator.parse_labels(irrep_lines, irrep_sec, "irrep_w_label"),
        "source_sg": generator.parse_ints(
            irrep_lines, irrep_sec, "irrep_w_space_group"
        ),
        "dimension": generator.parse_ints(
            irrep_lines, irrep_sec, "irrep_w_dimension"
        ),
        "type": generator.parse_ints(irrep_lines, irrep_sec, "irrep_w_type"),
        "count": generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_count"
        ),
        "pointer": generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_pointer"
        ),
        "irrep": generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_irrep"
        ),
        "frequency": generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_frequency"
        ),
        "parent": generator.parse_ints(isotropy_lines, isotropy_sec, "isotropy_parent"),
    }


def record_offsets(count):
    """Packed-table start offset of every record, in record order."""
    offsets = []
    total = 0
    for value in count:
        offsets.append(total)
        total += value
    return offsets


def check(data, expected=PINNED):
    """Return one message per violated invariant (empty means every gate passed)."""
    failures = []

    def require(condition, message):
        if not condition:
            failures.append(message)

    sections = data["irrep_sections"]
    require(
        sorted(sections) == sorted(W_SECTIONS),
        "other-wave-vector sections are "
        f"{sections}, expected exactly {sorted(W_SECTIONS)}: a k-vector or "
        "character section would make these rows computable",
    )
    sizes = {name: len(data[name]) for name in ("labels", "source_sg", "dimension", "type")}
    require(
        set(sizes.values()) == {expected.sources},
        f"other-wave-vector irrep tables are not all {expected.sources} long: {sizes}",
    )
    if failures:
        return failures

    count, packed = data["count"], data["irrep"]
    require(
        len(count) == expected.condensates,
        f"isotropy_w_subduce_count has {len(count)} entries, "
        f"expected {expected.condensates}",
    )
    require(
        len(data["parent"]) == len(count),
        "isotropy_w_subduce_count is not parallel to isotropy_parent",
    )
    if failures:
        return failures
    require(
        sum(count) == len(packed) == len(data["frequency"]),
        f"packed rows {len(packed)}/{len(data['frequency'])} do not match "
        f"the counts' sum {sum(count)}",
    )
    records = [index for index, value in enumerate(count) if value > 0]
    require(
        len(records) == expected.records,
        f"{len(records)} records carry other-wave-vector rows, "
        f"expected {expected.records}",
    )
    require(
        len(packed) == expected.rows,
        f"{len(packed)} packed rows, expected {expected.rows}",
    )
    if failures:
        return failures

    # The sparse pointer index must name exactly the start offsets of the
    # records with rows; two parallel pinned arrays disagreeing here would mean
    # the packed rows belong to different records than the counts say.
    offsets = record_offsets(count)
    starts = {offsets[index] + 1 for index in records}
    pointers = {value for value in data["pointer"] if value > 0}
    require(
        starts == pointers,
        f"pointer starts {sorted(pointers)[:4]}... differ from record starts "
        f"{sorted(starts)[:4]}... ({len(pointers)} vs {len(starts)})",
    )

    dimensions = data["dimension"]
    source_sg = data["source_sg"]
    parents = data["parent"]
    over_dimension = 0
    sg_mismatch = 0
    zero_frequency = 0
    out_of_range = 0
    entry = 0
    for index in records:
        for _ in range(count[index]):
            source = packed[entry]
            frequency = data["frequency"][entry]
            if not 1 <= source <= expected.sources:
                out_of_range += 1
            else:
                if source_sg[source - 1] != parents[index]:
                    sg_mismatch += 1
                if frequency > dimensions[source - 1]:
                    over_dimension += 1
            if frequency <= 0:
                zero_frequency += 1
            entry += 1
    require(out_of_range == 0, f"{out_of_range} irrep indices are out of range")
    require(
        sg_mismatch == 0,
        f"{sg_mismatch} rows pair a record with an irrep of another space group",
    )
    require(zero_frequency == 0, f"{zero_frequency} rows carry a zero frequency")
    require(
        over_dimension == 0,
        f"{over_dimension} frequencies exceed their irrep's dimension",
    )
    return failures


def main():
    try:
        data = parse_pinned()
    except (FileNotFoundError, ValueError) as error:
        print(f"pinned input unavailable: {error}", file=sys.stderr)
        return 2
    failures = check(data)
    for message in failures:
        print(f"FAIL {message}", file=sys.stderr)
    rows = len(data["irrep"])
    print(
        "other_wave_vector: sources={} records={} rows={} "
        "k_vector_data=absent character_data=absent checks_failed={}".format(
            len(data["labels"]),
            sum(1 for value in data["count"] if value > 0),
            rows,
            len(failures),
        )
    )
    if failures:
        return 1
    print(
        "note: the pinned archive stores these irreps as label/space group/"
        "dimension/type only, so their subduction frequencies are not "
        "computable from it; the audit reports them separately "
        "(--require-w-complete)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
