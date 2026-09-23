#!/usr/bin/env python3
"""What the pinned archive does and does not pin about the other-wave-vector rows.

`data_isotropy.txt` carries the `isotropy_w_subduce_*` arrays: 5756 entries over
1006 isotropy records, naming parent irreps **at other wave vectors** -- the
`DT`/`SM` star labels of the eleven cubic space groups 196, 202, 203, 209, 210,
216, 219, 225, 226, 227 and 228 (73 distinct labels).

Task 9's audit reports those rows separately instead of computing them, and this
script makes that decision executable:

* the 73 irreps they name are stored in `data_irreps.txt` as `irrep_w_label`,
  `irrep_w_space_group`, `irrep_w_dimension` and `irrep_w_type` only -- that file
  has **no** other-wave-vector k-vector table and **no** character/matrix table,
  so no engine can produce their frequencies from it.  (Their little-group
  tables are archived separately in `data_little.txt`, where all 73 sources are
  present with matching dimensions; decoding those is a separate follow-up, and
  the first paragraph of the gate below keeps the question visible.);
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


LATTICES = 14
K_SLOTS = 27
K_GROUPS = 4


def parse_pinned():
    """Read the pinned data files out of the SHA256-verified archive."""
    irrep_lines = generator.read_file("data_irreps.txt")
    isotropy_lines = generator.read_file("data_isotropy.txt")
    little_lines = generator.read_file("data_little.txt")
    space_lines = generator.read_file("data_space.txt")
    irrep_sec = generator.get_sections(irrep_lines)
    isotropy_sec = generator.get_sections(isotropy_lines)
    little_sec = generator.get_sections(little_lines)
    space_sec = generator.get_sections(space_lines)
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
        # Little-group tables: the wave vector of every little irrep, and with
        # them the reason the `w` rows cannot be given a numeric k.
        "little_k": generator.parse_ints(little_lines, little_sec, "little_k"),
        "little_k_count": generator.parse_ints(little_lines, little_sec, "little_k_count"),
        "little_k_label": generator.parse_labels(little_lines, little_sec, "little_k_label"),
        "little_label": generator.parse_labels(
            little_lines, little_sec, "little_irr_full_label"
        ),
        "little_sg": generator.parse_ints(
            little_lines, little_sec, "little_irr_space_group"
        ),
        "little_k_index": generator.parse_ints(little_lines, little_sec, "little_irr_k"),
        "little_dim": generator.parse_ints(
            little_lines, little_sec, "little_irr_full_dim"
        ),
        "sg_lattice": generator.parse_ints(space_lines, space_sec, "ispace_lattice"),
    }


def wave_vector(little_k, lattice, slot):
    """The four `(x, y, z, d)` groups of one (lattice, k slot).

    The pinned little table stores each wave vector as

    * the base point `k0` (the first group) and
    * up to three free-parameter directions (the remaining groups),

    so a point has no direction, a line has one, and the general position has
    three.  Groups are `(x, y, z, d)` = `(x/d, y/d, z/d)`.
    """
    start = (lattice * K_SLOTS + slot) * K_GROUPS * 4
    entries = little_k[start:start + K_GROUPS * 4]
    return [tuple(entries[index * 4:index * 4 + 4]) for index in range(K_GROUPS)]


def format_rational(group):
    x, y, z, d = group
    if (x, y, z) == (0, 0, 0):
        return "0"
    return f"({x},{y},{z})/{d}"


def little_wave_vectors(data, expected_sources):
    """Resolve every other-wave-vector source to its little-group wave vector.

    Returns `(resolved, failures)`, where each resolution is
    `(sg, label, dimension, lattice, slot label, k0, free directions)`.
    """
    failures = []
    little_k = data["little_k"]
    if len(little_k) != LATTICES * K_SLOTS * K_GROUPS * 4:
        return [], [
            f"little_k has {len(little_k)} entries, expected "
            f"{LATTICES * K_SLOTS * K_GROUPS * 4} (14 lattices x 27 k slots x 4 "
            "rational groups)"
        ]
    if len(data["little_k_count"]) != LATTICES:
        return [], [
            f"little_k_count has {len(data['little_k_count'])} entries, "
            f"expected {LATTICES} (one per Bravais lattice)"
        ]
    index = {}
    for position, (label, sg) in enumerate(
        zip(data["little_label"], data["little_sg"])
    ):
        index.setdefault((sg, label), []).append(position)
    resolved = []
    for label, sg in expected_sources:
        found = index.get((sg, label), [])
        if len(found) != 1:
            failures.append(
                f"other-wave-vector source {label} of space group {sg} has "
                f"{len(found)} little-table records, expected exactly one"
            )
            continue
        position = found[0]
        slot = data["little_k_index"][position] - 1
        lattice = data["sg_lattice"][sg - 1] - 1
        if not 0 <= lattice < LATTICES:
            failures.append(f"space group {sg} has lattice index {lattice + 1}")
            continue
        if not 0 <= slot < data["little_k_count"][lattice]:
            failures.append(
                f"other-wave-vector source {label} of space group {sg} names "
                f"k slot {slot + 1}, outside its lattice's "
                f"{data['little_k_count'][lattice]} slots"
            )
            continue
        groups = wave_vector(little_k, lattice, slot)
        free = [group for group in groups[1:] if group[:3] != (0, 0, 0)]
        slot_label = data["little_k_label"][lattice * K_SLOTS + slot]
        resolved.append(
            (sg, label, data["little_dim"][position], lattice, slot_label,
             format_rational(groups[0]), [format_rational(group) for group in free])
        )
    return resolved, failures


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
    if failures:
        return failures

    # Why the rows cannot be computed from the pinned data: every source is a
    # *parameterized* wave vector, so it has no single numeric k to fold.
    sources = sorted({(label, sg) for sg, label in zip(source_sg, data["labels"])})
    resolved, wave_failures = little_wave_vectors(data, sources)
    failures.extend(wave_failures)
    if not failures:
        unparameterized = [
            f"{label} (SG {sg})" for sg, label, _, _, _, _, free in resolved
            if not free
        ]
        require(
            not unparameterized,
            "other-wave-vector sources with a fixed k (a decode or pinning "
            f"change; they would become computable): {unparameterized}",
        )
    return failures


# Sources without a frozen little-group character table.  The list is empty:
# every one of the 73 sources is covered by `src/irrep/w_little_characters_data.rs`
# (65 from the Gamma compatibility rows, the eight SG 202/203/209/210
# `DT3`/`DT4` from the little-cogroup pair route), and the same list is pinned
# there as `W_LITTLE_CHARACTERS_UNRESOLVED` (regression
# `tests/w_little_characters.rs`).
UNRESOLVED_SOURCES = frozenset()


def frozen_coverage(data):
    """Rows whose source has a frozen character table, and rows that do not."""
    sources = list(zip(data["source_sg"], data["labels"]))
    packed = data["irrep"]
    frozen = blocked = 0
    blocked_sources = set()
    for index in packed:
        source = sources[index - 1]
        if (source[0], source[1].strip()) in UNRESOLVED_SOURCES:
            blocked += 1
            blocked_sources.add((source[0], source[1].strip()))
        else:
            frozen += 1
    return frozen, blocked, sorted(blocked_sources)


def main():
    try:
        data = parse_pinned()
    except (FileNotFoundError, ValueError) as error:
        print(f"pinned input unavailable: {error}", file=sys.stderr)
        return 2
    failures = check(data)
    for message in failures:
        print(f"FAIL {message}", file=sys.stderr)
    sources = sorted({(label, sg) for sg, label in zip(data["source_sg"], data["labels"])})
    resolved, _ = little_wave_vectors(data, sources)
    free_counts = sorted({len(free) for *_, free in resolved})
    directions = sorted({direction for *_, free in resolved for direction in free})
    print(
        "w_wave_vectors: sources={} resolved={} free_parameters={} "
        "base_points={} directions={}".format(
            len(sources),
            len(resolved),
            free_counts,
            sorted({base for *_, base, _ in resolved}),
            directions,
        )
    )
    rows = len(data["irrep"])
    frozen, blocked, blocked_sources = frozen_coverage(data)
    print(
        "frozen_characters: rows_with_frozen_table={} rows_blocked={} "
        "unresolved_sources={}".format(frozen, blocked, blocked_sources)
    )
    print(
        "other_wave_vector: sources={} records={} rows={} "
        "k_vectors=parameterized_lines character_data=absent checks_failed={}".format(
            len(data["labels"]),
            sum(1 for value in data["count"] if value > 0),
            rows,
            len(failures),
        )
    )
    if failures:
        return 1
    print(
        "note: the pinned irrep table stores these irreps as label/space group/"
        "dimension/type only.  Their little-group wave vectors are lines "
        "k = Gamma + t*v (one free parameter for all 73 sources, from the pinned "
        "little table; the direction components, the pinned k vectors and the "
        "frozen little-group operations are all in the parent's conventional "
        "reciprocal frame -- measured in R6.1, see docs/subduction-conventions.md "
        "§16 -- so there is no single numeric k to fold and the "
        "frequency needs the induced representation of the line: the frozen "
        "little-group character tables of all 73 sources are in "
        "src/irrep/w_little_characters_data.rs; the engine turns them into "
        "the 5,756 frequencies at the program's parameter t = 1/4 "
        "(`subduce_line_at_parameter` + `trivial_content()`; the R5 Gamma-only "
        "route `line_trivial_content_via_blocks` is kept as a second reading), "
        "and the "
        "audit reports that track separately under `--require-w-complete` "
        "(both gates exit 0 with `engine_errors=0`).  Their values also have a "
        "live oracle: scripts/verify_w_subduction_oracle.py compares all 5,756 "
        "rows against the official program's SHOW FREQUENCY list per parent "
        "irrep."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
