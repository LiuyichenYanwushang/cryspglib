#!/usr/bin/env python3
"""Generate the per-ordinal frozen embedding settings of task 9.

The engine's full-subduction layer (`src/irrep/subduction.rs`) embeds one
isotropy subgroup record into its parent.  The embedding needs two conventions
that the stored isotropy tables alone do not fix:

* the setting transform ``U`` relating the stored basis ``W`` (parent primitive
  frame) to the subgroup's own ITA conventional cell ``B``,
  ``W . P_parent = U . (P_sub . B)``,
* the child-frame origin shift ``delta`` between the ITA setting the isotropy
  tables were recorded in and the shipped ``SG_DATA_HALL`` frame of the
  subgroup, applied as ``t' = t + delta - R delta``.

Task 8 froze these per ``(parent, subgroup)`` pair for the 13 pairs its fixtures
cover.  This generator replaces that table with one entry **per isotropy
record ordinal** for those 69 records and six task-9 label/shear witnesses,
derived from the official ISOTROPY program.

Every entry is derived, never searched:

* the machine tables are read from the pinned ``isotropy_subgroup/iso.zip``
  archive (SHA-256 gated; extracted directories are ignored),
* the official ``iso`` binary is run once per record in an isolated temporary
  directory with ``SET I ALL OR 1`` (the setting the pinned tables were recorded
  in) and ``VALUE PARENT/IRREP/DIRECTION`` + ``SHOW BASIS/ORIGIN/SIZE/ELEMENTS``
  / ``DISPLAY ISOTROPY``,
* ``U = W . P_parent . (P_sub . B_oracle)^-1`` must be an exact integer matrix
  with ``|det| = 1`` (no signed-permutation restriction),
* ``delta`` is the **origin difference** between the recorded ITA setting and
  the printed setting the shipped ``SG_DATA_HALL[subgroup]`` row belongs to,
  mapped into the child cell.  The setting is selected by exact equality of the
  full expanded ``SHOW ELEMENTS`` operation set (both ``SET I ALL OR 1`` and
  ``OR 2`` are queried), so ``delta`` is a geometric quantity read off the
  printed origins; a Hall row that matches no setting fails, and a shift is
  never fitted to make a single operation match.

The origin column is an independent check: the stored origin converted through
the parent primitive basis must equal the printed origin exactly.

Usage::

    python3 scripts/generate_subduction_settings.py            # --check
    python3 scripts/generate_subduction_settings.py --check --evidence /tmp/x.jsonl

Scope of this tool today: it derives the conventions of its **75 legacy records**
(69 task-8 records plus six label/shear witnesses) from the official oracle and
checks those rows of the committed module.  The shipped
``src/irrep/subduction_settings_data.rs`` is the **15,239-record full table**
assembled by ``scripts/task9/build_table.py`` (see ``scripts/task9/README.md``);
this tool is the library those task-9 scripts import, and it deliberately has no
writer anymore -- ``--write`` fails instead of replacing the full table with its
75-record subset.  ``--check`` (the default) re-queries the official oracle,
compares the 75 derived rows with the committed module and verifies that the
module covers every pinned ordinal.  It never writes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import traceback
import zipfile
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir))
ISO_DIR = os.path.join(ROOT, "isotropy_subgroup")
OUTPUT = os.path.join(ROOT, "src", "irrep", "subduction_settings_data.rs")
GENERATED_DATA = os.path.join(ROOT, "src", "irrep", "generated_data.rs")

sys.path.insert(0, HERE)
import generate_irrep_data as machine  # noqa: E402
import verify_isotropy_operations as operations  # noqa: E402
import verify_isotropy_oracle as geometry  # noqa: E402
import audit_subduction_settings as census  # noqa: E402

# The engine snaps Hall translations to this grid (`SOURCE_TRANSLATION_GRID` in
# `src/irrep/subduction.rs`); every Hall row is checked against it exactly.
SOURCE_TRANSLATION_GRID = 12

# The records this table covers, in the order the task-8 fixtures enumerate
# them: the nine original pair rules (59 records) plus the four scalar
# self-restrictions (10 records).  The counts are part of the gate: a changed
# pinned table must fail here instead of silently changing the table's scope.
COVERAGE = [
    ((221, 83), 1),
    ((221, 12), 16),
    ((221, 148), 4),
    ((221, 123), 11),
    ((221, 47), 9),
    ((225, 8), 6),
    ((16, 22), 4),
    ((167, 15), 7),
    ((139, 126), 1),
    ((19, 19), 1),
    ((23, 23), 1),
    ((45, 45), 1),
    ((83, 83), 7),
]

# The pair rules the previous `FROZEN_EMBEDDINGS` table used, as `U` rows.  The
# generated per-ordinal entries must reproduce them record by record; this is
# the "all 69 match the current rules" probe, kept in the generator so the
# comparison is reproducible without any ad-hoc evidence file.
PREVIOUS_RULES = {
    (221, 83): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (221, 12): [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    (221, 148): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (221, 123): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (221, 47): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (225, 8): [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    (16, 22): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (167, 15): [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    (139, 126): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (19, 19): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (23, 23): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (45, 45): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    (83, 83): [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
}

# The previously frozen child shifts.  Only #126 needs one: the isotropy record
# is printed in ITA origin choice 1 while the shipped Hall row uses origin
# choice 2, which is exactly the (1/4, 1/4, 1/4) shift in the child cell.
PREVIOUS_SHIFTS = {
    (139, 126): (1, 1, 1, 4),
}
ZERO_SHIFT = (0, 0, 0, 1)

IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
CYCLIC = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
ZERO_VECTOR = [Fraction(0), Fraction(0), Fraction(0)]

# Task-9 witnesses, selected by ordinal rather than extending any pair rule.
# SHOW BASIS and the full SHOW ELEMENTS set independently fix these settings.
EXTRA_RECORDS = {
    427: (21, "R1", "P1", 22, [[1, -1, 0], [1, 0, -1], [1, 0, 0]]),
    984: (43, "GM1", "P1", 43, IDENTITY),
    1942: (64, "GM2+", "P1", 14, [[2, 0, 1], [-1, 0, -1], [0, 1, 0]]),
    2102: (67, "GM2+", "P1", 13, [[2, 0, 1], [-1, 0, -1], [0, 1, 0]]),
    15125: (230, "GM4-", "P2", 43, IDENTITY),
    15131: (230, "GM5-", "P2", 43, IDENTITY),
}


class GenerationError(RuntimeError):
    """A derivation or oracle check failed; the table is not generated."""


# ── exact linear algebra (independent of the Rust engine) ────────────────────


def matmul(left, right):
    return [
        [sum(left[i][t] * right[t][j] for t in range(3)) for j in range(3)]
        for i in range(3)
    ]


def matvec(matrix, vector):
    return [sum(matrix[i][t] * vector[t] for t in range(3)) for i in range(3)]


def inverse(matrix):
    """Exact inverse of a 3x3 rational matrix, or ``None`` if singular."""
    determ = Fraction(geometry.det3(matrix))
    if determ == 0:
        return None
    cofactors = [[None] * 3 for _ in range(3)]
    for row in range(3):
        for column in range(3):
            rows = [index for index in range(3) if index != row]
            columns = [index for index in range(3) if index != column]
            minor = Fraction(
                matrix[rows[0]][columns[0]] * matrix[rows[1]][columns[1]]
                - matrix[rows[0]][columns[1]] * matrix[rows[1]][columns[0]]
            )
            sign = 1 if (row + column) % 2 == 0 else -1
            cofactors[column][row] = sign * minor / determ
    return cofactors


def solve(matrix, vector):
    """Exact solution of ``matrix . x = vector``, or ``None`` if singular."""
    rows = [
        [Fraction(matrix[i][j]) for j in range(3)] + [Fraction(vector[i])]
        for i in range(3)
    ]
    for column in range(3):
        pivot = next((row for row in range(column, 3) if rows[row][column] != 0), None)
        if pivot is None:
            return None
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(3):
            if row != column and rows[row][column] != 0:
                factor = rows[row][column]
                rows[row] = [rows[row][j] - factor * rows[column][j] for j in range(4)]
    return [rows[i][3] for i in range(3)]


def canonical_translation(vector):
    """Representative of ``vector`` modulo the conventional cell ``Z^3``."""
    return tuple(
        value - Fraction(value.numerator // value.denominator) for value in vector
    )


def integers(matrix):
    return all(value.denominator == 1 for row in matrix for value in row)


def int_matrix(matrix):
    return [[int(value) for value in row] for row in matrix]


def encode_shift(vector):
    """``(x, y, z, d)`` encoding of an exact translation, on a common denominator."""
    denominator = 1
    for value in vector:
        denominator = denominator * value.denominator // _gcd(denominator, value.denominator)
    return tuple(
        [int(value * denominator) for value in vector] + [denominator]
    )


def _gcd(left, right):
    while right:
        left, right = right, left % right
    return abs(left)


# ── pinned inputs ────────────────────────────────────────────────────────────


def load_machine_records():
    """The complete pinned isotropy record table, in raw ordinal order."""
    lines = machine.read_file("data_isotropy.txt")
    sections = machine.get_sections(lines)
    parent = machine.parse_ints(lines, sections, "isotropy_parent")
    irrep = machine.parse_ints(lines, sections, "isotropy_irrep")
    subgroup = machine.parse_ints(lines, sections, "isotropy_subgroup")
    basis = machine.parse_ints(lines, sections, "isotropy_basis")
    origin = machine.parse_ints(lines, sections, "isotropy_origin")
    label = machine.parse_labels(lines, sections, "isotropy_orderparam_label")
    irrep_lines = machine.read_file("data_irreps.txt")
    irrep_sections = machine.get_sections(irrep_lines)
    irrep_labels = machine.parse_labels(irrep_lines, irrep_sections, "irrep_label")
    machine.require_per_record(basis, 9, len(parent), "isotropy_basis")
    machine.require_per_record(origin, 4, len(parent), "isotropy_origin")
    records = []
    for index in range(len(parent)):
        records.append(
            {
                "ordinal": index,
                "parent": parent[index],
                "ml": irrep_labels[irrep[index] - 1],
                "label": label[index],
                "child": subgroup[index],
                "basis": [
                    [basis[index * 9 + row * 3 + column] for column in range(3)]
                    for row in range(3)
                ],
                "origin": origin[index * 4:index * 4 + 4],
            }
        )
    return records


def select_records(records):
    """The original 69 records plus six independently pinned witnesses."""
    wanted = dict(COVERAGE)
    selected = [
        record
        for record in records
        if (record["parent"], record["child"]) in wanted
    ]
    counts = {}
    for record in selected:
        key = (record["parent"], record["child"])
        counts[key] = counts.get(key, 0) + 1
    if counts != wanted:
        raise GenerationError(
            f"covered record counts changed: expected {wanted}, found {counts}"
        )
    for ordinal, (parent, ml, label, child, _) in EXTRA_RECORDS.items():
        record = records[ordinal]
        if (record["parent"], record["ml"], record["label"], record["child"]) != (parent, ml, label, child):
            raise GenerationError(f"ordinal {ordinal}: task-9 source identity changed")
        selected.append(record)
    selected.sort(key=lambda record: record["ordinal"])
    ordinals = [record["ordinal"] for record in selected]
    if len(set(ordinals)) != len(ordinals):
        raise GenerationError("duplicate isotropy ordinals in the covered selection")
    for record in selected:
        if records[record["ordinal"]] is not record:
            raise GenerationError(
                f"ordinal {record['ordinal']} does not address its own raw record"
            )
    return selected


def load_data_hall():
    """``SG_DATA_HALL`` from the committed generated table the engine reads."""
    with open(GENERATED_DATA, encoding="utf-8") as handle:
        text = handle.read()
    marker = "pub static SG_DATA_HALL"
    if marker not in text:
        raise GenerationError(f"{GENERATED_DATA}: SG_DATA_HALL not found")
    body = text[text.index(marker):].split("=", 1)[1]
    body = body[body.index("[") + 1:body.index("];")]
    values = []
    for line in body.splitlines():
        for token in line.split("//")[0].replace(",", " ").split():
            values.append(int(token))
    if len(values) != 231:
        raise GenerationError(f"SG_DATA_HALL has {len(values)} entries, expected 231")
    expected = [0] * 231
    for frame in machine.load_committed_data_hall_provenance().frames:
        expected[frame.spacegroup] = frame.data_hall
    if values != expected:
        raise GenerationError("SG_DATA_HALL disagrees with the pinned ISO--IR frame provenance")
    return values


def load_hall_operations():
    """The digest-pinned Hall operation table, as exact 1/12-grid rationals."""
    table = {}
    for _sg, halls in machine._load_hall_operations().items():
        for hall, rotations, translations in halls:
            if hall in table:
                raise GenerationError(f"hall_operations.json repeats Hall {hall}")
            exact_operations = []
            for rotation, translation in zip(rotations, translations):
                exact = []
                for value in translation:
                    scaled = value * SOURCE_TRANSLATION_GRID
                    rounded = round(scaled)
                    if abs(scaled - rounded) > 1e-9:
                        raise GenerationError(
                            f"Hall {hall} translation {value} is off the "
                            f"1/{SOURCE_TRANSLATION_GRID} grid"
                        )
                    exact.append(Fraction(rounded, SOURCE_TRANSLATION_GRID))
                exact_operations.append((tuple(int(v) for v in rotation), tuple(exact)))
            table[hall] = {"sg": _sg, "operations": exact_operations}
    return table


def load_point_op_legend():
    """Decode the program's operation legend from the pinned archive."""
    with zipfile.ZipFile(os.path.join(ISO_DIR, "iso.zip")) as archive:
        payload = archive.read("data_space.txt")
    with tempfile.TemporaryDirectory() as temporary:
        path = os.path.join(temporary, "data_space.txt")
        with open(path, "wb") as handle:
            handle.write(payload)
        return operations.load_point_op_legend(path)


# ── official oracle queries ──────────────────────────────────────────────────


def query_record(parent, ml, label, legend, origin_choice=1, *, oracle_dir):
    """One record's full row (basis/origin/size/elements) from the oracle.

    The operation gate's parser is reused. This runner parameterizes
    ``SET I ALL OR <n>`` so the same
    isolated-process protocol can also print the other ITA origin choice, which
    is what pins the child shift.  The row is normalized so both query paths
    compare exact rationals (`direction` is the printed direction label).
    """
    binary = os.path.join(oracle_dir, "iso")
    if not os.path.isfile(binary):
        raise GenerationError(f"the official oracle is missing: {binary}")
    commands = [
        "PAGE 1000",
        "SC 250",
        f"SET I ALL OR {origin_choice}",
        f"VALUE PARENT {parent}",
        f"VALUE IRREP {ml}",
        f"VALUE DIRECTION {label}",
        "SHOW SUBGROUP",
        "SHOW DIRECTION",
        "SHOW BASIS",
        "SHOW ORIGIN",
        "SHOW SIZE",
        "SHOW ELEMENTS",
        "DISPLAY ISOTROPY",
        "QUIT",
    ]
    try:
        with tempfile.TemporaryDirectory() as workdir:
            result = subprocess.run(
                [binary],
                input="\n".join(commands) + "\n",
                capture_output=True,
                text=True,
                cwd=workdir,
                env=dict(os.environ, ISODATA=oracle_dir + os.sep),
                timeout=300,
            )
        if result.returncode != 0:
            raise RuntimeError(f"iso exited with {result.returncode}: {result.stderr}")
        if operations.BANNER not in result.stdout:
            raise RuntimeError(
                f"iso did not report the pinned setting: {result.stdout[:200]!r}"
            )
        row = operations.parse_row(result.stdout)
    except Exception as error:  # noqa: BLE001 - reported with record context
        raise GenerationError(
            f"SG {parent} {ml} {label} OR{origin_choice}: the oracle failed: {error}"
        ) from error
    if row["direction"] != label:
        raise GenerationError(
            f"SG {parent} {ml}: the oracle answered for {row['direction']!r}, "
            f"not {label!r}"
        )
    try:
        elements = operations.parse_elements(row.pop("elements_text"), legend)
    except Exception as error:  # noqa: BLE001 - reported with record context
        raise GenerationError(
            f"SG {parent} {ml} {label}: malformed SHOW ELEMENTS: {error}"
        ) from error
    return {
        "subgroup": row["subgroup"],
        "size": row["size"],
        "direction": row["direction"],
        "basis": [[Fraction(value) for value in vector] for vector in row["basis"]],
        "origin": [Fraction(value) for value in row["origin"]],
        "elements": elements,
    }


# ── derivations ──────────────────────────────────────────────────────────────


def check_unimodular(record, setting):
    """Reject a setting transform that is not an exact unimodular integer matrix.

    The derivation below can only produce such a matrix when the two bases span
    the same lattice, so this gate is defensive; it is also the shared
    contract the committed table is validated against.
    """
    if not integers(setting):
        raise GenerationError(
            f"ordinal {record['ordinal']}: U = {setting} is not an integer matrix"
        )
    determ = geometry.det3(setting)
    if abs(determ) != 1:
        raise GenerationError(
            f"ordinal {record['ordinal']}: |det U| = {abs(determ)}, not 1"
        )
    return int_matrix(setting)


def derive_setting(record, oracle_basis):
    """``U = W . P_parent . (P_sub . B_oracle)^-1`` with exact integer checks."""
    parent_basis = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[record["parent"]]]
    child_basis = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[record["child"]]]
    stored = [[Fraction(value) for value in row] for row in record["basis"]]
    printed = [[Fraction(value) for value in row] for row in oracle_basis]
    left = matmul(stored, parent_basis)
    right = matmul(child_basis, printed)
    if not geometry.same_lattice(left, right):
        raise GenerationError(
            f"ordinal {record['ordinal']}: stored basis {left} and the printed "
            f"cell {printed} (times the subgroup centring) span different lattices"
        )
    inverse_right = inverse(right)
    if inverse_right is None:
        raise GenerationError(f"ordinal {record['ordinal']}: the printed basis is singular")
    return check_unimodular(record, matmul(left, inverse_right))


def check_origin(record, oracle_origin):
    """The stored origin converted through ``P_parent`` must match exactly."""
    parent_basis = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[record["parent"]]]
    x, y, z, denominator = record["origin"]
    if denominator <= 0:
        raise GenerationError(
            f"ordinal {record['ordinal']}: origin denominator {denominator}"
        )
    stored = [
        Fraction(x, denominator),
        Fraction(y, denominator),
        Fraction(z, denominator),
    ]
    converted = geometry.to_conventional(stored, parent_basis)
    printed = [Fraction(value) for value in oracle_origin]
    if converted != printed:
        raise GenerationError(
            f"ordinal {record['ordinal']}: stored origin {converted} != printed "
            f"{printed} exactly"
        )
    return converted


def centring_classes(subgroup_sg):
    """The finite group ``L_sub mod Z^3``: the child's centring translations.

    The official ``SHOW ELEMENTS`` column lists one representative per point
    operation *modulo the subgroup lattice*, while the shipped Hall row lists
    the full conventional cell (one entry per centring vector).  Expanding the
    printed rows by these classes is what makes the two comparable; leaving out
    the zero class silently drops the primitive representative.
    """
    centring = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[subgroup_sg]]
    # Rows are the primitive vectors in conventional coordinates, so the
    # conventional cell holds `1 / |det|` of them (4 for F, 3 for R, 2 for I/C).
    order = round(abs(Fraction(1) / geometry.det3(centring)))
    classes = set()
    for first in range(order):
        for second in range(order):
            for third in range(order):
                vector = [
                    first * centring[0][i] + second * centring[1][i] + third * centring[2][i]
                    for i in range(3)
                ]
                classes.add(canonical_translation(vector))
    if len(classes) != order:
        raise GenerationError(
            f"SG {subgroup_sg}: the centring classes do not form a group of "
            f"order {order}"
        )
    return sorted(classes)


def child_frame_operations(record, oracle_basis, oracle_origin, elements):
    """The official elements expressed in the subgroup's own conventional cell."""
    basis = [[Fraction(value) for value in row] for row in oracle_basis]
    origin = [Fraction(value) for value in oracle_origin]
    transform = [[basis[j][i] for j in range(3)] for i in range(3)]
    inverse_transform = inverse(transform)
    if inverse_transform is None:
        raise GenerationError(f"ordinal {record['ordinal']}: singular child basis")
    classes = centring_classes(record["child"])
    expanded = {}
    for element in elements:
        rotation = [
            [Fraction(element["rotation"][3 * i + j]) for j in range(3)]
            for i in range(3)
        ]
        translation = [Fraction(value) for value in element["translation"]]
        child_rotation = matmul(matmul(inverse_transform, rotation), transform)
        if not integers(child_rotation):
            raise GenerationError(
                f"ordinal {record['ordinal']}: {element['raw']} does not invert to "
                "an integer rotation in the child cell"
            )
        rotated_origin = matvec(rotation, origin)
        child_translation = matvec(
            inverse_transform,
            [translation[i] + rotated_origin[i] - origin[i] for i in range(3)],
        )
        key = tuple(int(value) for row in child_rotation for value in row)
        for vector in classes:
            shifted = [child_translation[i] + vector[i] for i in range(3)]
            expanded.setdefault(key, set()).add(canonical_translation(shifted))
    return expanded


def select_setting(record, expansions, hall_operations):
    """The printed ITA settings whose full operation set equals the Hall row.

    Equality is over the **whole expanded operation set** (all rotations and
    all centring classes), never a single operation, and a Hall row that
    matches no printed setting fails the generator.
    """
    hall = {}
    for rotation, translation in hall_operations:
        hall.setdefault(rotation, set()).add(canonical_translation(translation))
    matches = [
        choice for choice, expanded in expansions.items() if expanded == hall
    ]
    if not matches:
        raise GenerationError(
            f"ordinal {record['ordinal']}: the Hall row matches none of the "
            f"printed ITA settings ({sorted(expansions)})"
        )
    return sorted(matches)


def setting_shift(record, basis, origins, choices):
    """``delta = T^-1 (o_setting - o_recorded)`` for the matching settings.

    Settings that print the same operations are the same setting only if they
    also print the same origin; two different origins would make the shift
    ambiguous, so exactly one value must survive.
    """
    transform = [[basis[j][i] for j in range(3)] for i in range(3)]
    inverse_transform = inverse(transform)
    if inverse_transform is None:
        raise GenerationError(f"ordinal {record['ordinal']}: singular child basis")
    candidates = set()
    for choice in choices:
        difference = [origins[choice][i] - origins[1][i] for i in range(3)]
        candidates.add(canonical_translation(matvec(inverse_transform, difference)))
    if len(candidates) != 1:
        raise GenerationError(
            f"ordinal {record['ordinal']}: the Hall row matches settings "
            f"{choices} with {len(candidates)} different origin shifts"
        )
    return list(candidates)[0]


def origin_choice_shift(record, basis, expansions, origins, hall_operations):
    """The child origin shift between the recorded ITA setting and the Hall row.

    The isotropy tables (and therefore ``U``) are recorded under
    ``SET I ALL OR 1``, which the geometry gate pins.  Which ITA setting the
    shipped Hall row uses is decided by *exact equality of the full expanded
    operation set* against the printed settings, and the shift itself is then
    the printed **origin difference** mapped into the child cell
    (``delta = T^-1 (o_setting - o_recorded)``) -- a geometric quantity, not a
    translation fitted to make one operation match.
    """
    choices = select_setting(record, expansions, hall_operations)
    delta = setting_shift(record, basis, origins, choices)
    # The shift must reproduce the official operation set from the Hall row.
    hall = {}
    for rotation, translation in hall_operations:
        hall.setdefault(rotation, set()).add(canonical_translation(translation))
    recalculated = {}
    for rotation, translations in hall.items():
        matrix = [[rotation[3 * i + j] for j in range(3)] for i in range(3)]
        correction = matvec(
            [[(1 if i == j else 0) - matrix[i][j] for j in range(3)] for i in range(3)],
            delta,
        )
        recalculated[rotation] = {
            canonical_translation([translation[i] + correction[i] for i in range(3)])
            for translation in translations
        }
    if recalculated != expansions[1]:
        raise GenerationError(
            f"ordinal {record['ordinal']}: the derived child shift "
            f"{[str(value) for value in delta]} does not reproduce the recorded "
            "operation set"
        )
    return delta


# ── assembly ─────────────────────────────────────────────────────────────────


def build_entries():
    """Derive every covered entry from the pinned inputs and the official oracle."""
    machine._verify_pinned_archives()
    with tempfile.TemporaryDirectory(prefix="subduction-settings-") as temporary:
        oracle_dir, _ = census.extract_oracle(temporary)
        return _build_entries(oracle_dir)


def _build_entries(oracle_dir):
    records = load_machine_records()
    selected = select_records(records)
    data_hall = load_data_hall()
    hall_table = load_hall_operations()
    legend = load_point_op_legend()
    entries = []
    for record in selected:
        row = query_record(record["parent"], record["ml"], record["label"], legend,
                           oracle_dir=oracle_dir)
        if row["subgroup"] != record["child"]:
            raise GenerationError(
                f"ordinal {record['ordinal']}: the oracle printed subgroup "
                f"{row['subgroup']}, the pinned table says {record['child']}"
            )
        stored_size = abs(geometry.det3(record["basis"]))
        if row["size"] != stored_size:
            raise GenerationError(
                f"ordinal {record['ordinal']}: the oracle printed size {row['size']}, "
                f"the pinned basis has |det| = {stored_size}"
            )
        check_origin(record, row["origin"])
        setting = derive_setting(record, row["basis"])
        hall_number = data_hall[record["child"]]
        other_setting = query_record(
            record["parent"], record["ml"], record["label"], legend, origin_choice=2,
            oracle_dir=oracle_dir,
        )
        if any(other_setting[field] != row[field] for field in ("subgroup", "size", "basis")):
            raise GenerationError(
                f"ordinal {record['ordinal']}: the two ITA settings print "
                "different subgroups or bases; an origin-only correction is invalid"
            )
        expansions = {}
        for choice, chosen in ((1, row), (2, other_setting)):
            if len(chosen["basis"]) != 3:
                raise GenerationError(
                    f"ordinal {record['ordinal']}: malformed basis in OR{choice}"
                )
            expansions[choice] = child_frame_operations(
                record, chosen["basis"], chosen["origin"], chosen["elements"]
            )
        origins = {1: row["origin"], 2: other_setting["origin"]}
        child_shift = origin_choice_shift(
            record,
            row["basis"],
            expansions,
            origins,
            hall_table[hall_number]["operations"],
        )
        pair = (record["parent"], record["child"])
        extra = EXTRA_RECORDS.get(record["ordinal"])
        previous = extra[4] if extra is not None else PREVIOUS_RULES[pair]
        if setting != previous:
            raise GenerationError(
                f"ordinal {record['ordinal']}: derived U = {setting} differs from "
                f"the frozen pair rule {previous}"
            )
        expected_shift = PREVIOUS_SHIFTS.get(pair, ZERO_SHIFT)
        encoded_shift = encode_shift(child_shift)
        if encoded_shift != expected_shift:
            raise GenerationError(
                f"ordinal {record['ordinal']}: derived child shift "
                f"{encoded_shift} differs from the recorded {expected_shift}"
            )
        entries.append(
            {
                "ordinal": record["ordinal"],
                "parent": record["parent"],
                "ml": record["ml"],
                "label": record["label"],
                "child": record["child"],
                "size": row["size"],
                "setting": setting,
                "child_shift": list(encoded_shift),
                "origin_choice": (
                    "ITA origin choice 1 vs 2"
                    if encoded_shift != ZERO_SHIFT
                    else ""
                ),
            }
        )
    entries.sort(key=lambda entry: entry["ordinal"])
    return entries


# ── Rust emission ────────────────────────────────────────────────────────────


def parse_committed(text):
    """Parse the committed entries back, for semantic comparison and tests.

    The shipped module is the 15,239-record full table, so one entry is
    `(ordinal, parent, child, U numerator, U denominator, child shift)`; a
    caller that only wants its own legacy rows filters by ordinal.
    """
    marker = "pub static FROZEN_EMBEDDING_SETTINGS"
    if marker not in text:
        raise GenerationError(f"{OUTPUT}: FROZEN_EMBEDDING_SETTINGS not found")
    body = text[text.index(marker):].split("=", 1)[1]
    body = body[body.index("[") + 1:body.rindex("];")]
    constants = {"IDENTITY": IDENTITY, "CYCLIC": CYCLIC}
    shifts = {"NO_SHIFT": list(ZERO_SHIFT)}
    entries = []
    for raw_line in body.splitlines():
        line = raw_line.split("//")[0].strip()
        if not line.startswith("("):
            continue
        fields = _split_top_level(line.strip().rstrip(",").strip("()"))
        if len(fields) != 6:
            raise GenerationError(f"cannot parse entry: {line!r}")
        ordinal, parent, child = (int(fields[i]) for i in range(3))
        setting_text, denominator_text, shift_text = fields[3], fields[4], fields[5]
        if setting_text in constants:
            setting = constants[setting_text]
        else:
            setting = [
                [int(value) for value in row.split(",")]
                for row in setting_text.strip("[]").split("], [")
            ]
        denominator = int(denominator_text)
        if denominator <= 0:
            raise GenerationError(f"non-positive setting denominator: {line!r}")
        if shift_text in shifts:
            shift = list(shifts[shift_text])
        else:
            shift = [int(value) for value in shift_text.strip("[]").split(",")]
        entries.append(
            {
                "ordinal": ordinal,
                "parent": parent,
                "child": child,
                "setting": setting,
                "setting_denominator": denominator,
                "child_shift": shift,
            }
        )
    return entries


def entry_addresses_record(entry, records):
    """Whether a committed entry's ordinal really addresses its own record.

    The engine's fail-closed lookup relies on the same identity: an ordinal
    whose pinned record is a different ``(parent, subgroup)`` pair must be
    rejected, never reinterpreted.
    """
    ordinal = entry["ordinal"]
    if not 0 <= ordinal < len(records):
        return False
    record = records[ordinal]
    return (
        record["parent"] == entry["parent"] and record["child"] == entry["child"]
    )


def check_ordinal_coverage(committed, record_count):
    """The committed module must address every pinned record exactly once.

    Count equality is not coverage: a module whose ordinals are ``[0, 1, 1]``
    over three pinned records has the right length and each surviving row still
    names its own record, yet ordinal 2 is missing -- and `compare()` keys by
    ordinal, so the duplicate would collapse silently.  This is the single
    implementation the offline tests and the online ``--check`` both call, so a
    test can never be stricter than the gate it is meant to pin.
    """
    ordinals = [entry["ordinal"] for entry in committed]
    expected = list(range(record_count))
    if ordinals == expected:
        return
    duplicated = sorted({value for value in ordinals if ordinals.count(value) > 1})
    missing = sorted(set(expected) - set(ordinals))
    outside = sorted(set(ordinals) - set(expected))
    raise GenerationError(
        f"the committed module must address every pinned record exactly once and in "
        f"ordinal order, but {len(ordinals)} entries leave {len(missing)} ordinals "
        f"uncovered (first {missing[:5]}), repeat {len(duplicated)} ordinals "
        f"(first {duplicated[:5]}) and carry {len(outside)} ordinals outside the "
        f"pinned range (first {outside[:5]})"
    )


def _split_top_level(text):
    """Split on commas that are not inside brackets."""
    fields = []
    depth = 0
    current = ""
    for character in text:
        if character in "[":
            depth += 1
        elif character in "]":
            depth -= 1
        if character == "," and depth == 0:
            fields.append(current.strip())
            current = ""
        else:
            current += character
    fields.append(current.strip())
    return fields


def rational_setting(entry):
    """`U` of one parsed entry as exact fractions (`numerator / denominator`)."""
    denominator = entry.get("setting_denominator", 1)
    return [
        [Fraction(int(value), denominator) for value in row] for row in entry["setting"]
    ]


def compare(entries, committed):
    """Check the derived legacy rows against the committed full table."""
    by_ordinal = {entry["ordinal"]: entry for entry in committed}
    for derived in entries:
        stored = by_ordinal.get(derived["ordinal"])
        if stored is None:
            raise GenerationError(
                f"ordinal {derived['ordinal']}: the committed module has no entry"
            )
        for field in ("parent", "child"):
            if derived[field] != stored[field]:
                raise GenerationError(
                    f"ordinal {derived['ordinal']}: derived {field} {derived[field]} "
                    f"!= committed {stored[field]}"
                )
        if rational_setting(derived) != rational_setting(stored):
            raise GenerationError(
                f"ordinal {derived['ordinal']}: derived U {derived['setting']} "
                f"!= committed {stored['setting']}/{stored['setting_denominator']}"
            )
        if list(derived["child_shift"]) != list(stored["child_shift"]):
            raise GenerationError(
                f"ordinal {derived['ordinal']}: derived child shift "
                f"{derived['child_shift']} != committed {stored['child_shift']}"
            )


def archive_hashes():
    hashes = {}
    for name, path in (
        ("iso.zip", os.path.join(ISO_DIR, "iso.zip")),
        ("hall_operations.json", os.path.join(HERE, "hall_operations.json")),
    ):
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        hashes[name] = digest.hexdigest()
    return hashes


# ── entry point ──────────────────────────────────────────────────────────────


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--write",
        action="store_true",
        help="refused: the full table is assembled by scripts/task9/build_table.py",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="query the official oracle and compare with the committed module",
    )
    parser.add_argument(
        "--evidence",
        metavar="PATH",
        help="write the derived entries as JSON lines (reproducible probe)",
    )
    arguments = parser.parse_args(argv)
    if arguments.write:
        raise GenerationError(
            f"{os.path.relpath(OUTPUT, ROOT)} is the 15,239-record table assembled "
            "by scripts/task9/build_table.py (see scripts/task9/README.md); this "
            f"tool derives only its {len(select_records(load_machine_records()))} "
            "legacy records and must not replace it. Run "
            "`python3 scripts/task9/build_table.py --check` against the same "
            "evidence to verify the full table."
        )
    hashes = archive_hashes()
    entries = build_entries()
    if arguments.evidence:
        with open(arguments.evidence, "w", encoding="utf-8") as handle:
            for entry in entries:
                handle.write(json.dumps(entry, sort_keys=True) + "\n")
    with open(OUTPUT, encoding="utf-8") as handle:
        committed_text = handle.read()
    committed = parse_committed(committed_text)
    records = load_machine_records()
    # The shipped module is the full table: it must address every pinned record
    # **exactly once, in ordinal order** (`check_ordinal_coverage`, the same
    # function the offline tests call), and every row must name its own record.
    # This tool checks the coverage and its own 75 rows; the conventions of the
    # other rows are checked by `scripts/task9/build_table.py --check`.
    check_ordinal_coverage(committed, len(records))
    for entry in committed:
        if not entry_addresses_record(entry, records):
            raise GenerationError(
                f"ordinal {entry['ordinal']}: the committed entry does not address "
                f"its own pinned record (claims parent {entry['parent']}, subgroup "
                f"{entry['child']})"
            )
    compare(entries, committed)
    nonzero = [
        entry for entry in entries if list(entry["child_shift"]) != list(ZERO_SHIFT)
    ]
    pairs = {(entry["parent"], entry["child"]) for entry in entries}
    print(
        f"checked {len(entries)} legacy entries against the official oracle "
        f"({len(pairs)} pairs, {len(nonzero)} non-zero child shift); the committed "
        f"module covers all {len(records)} pinned records"
    )
    print(f"pinned archives: {hashes}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except GenerationError as error:
        traceback.print_exc()
        print(f"FAILED: {error}", file=sys.stderr)
        sys.exit(1)
