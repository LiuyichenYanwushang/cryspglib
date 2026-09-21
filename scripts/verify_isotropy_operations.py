#!/usr/bin/env python3
"""Pin the subgroup operation embeddings against the ISO `iso` program.

Task 2 of `docs/full-irrep-subduction-plan.md`.  The geometry gate
(`verify_isotropy_oracle.py`) pins *which* subgroup a condensing direction
selects and its lattice/origin; this gate pins the **operations** of that
subgroup in the parent frame, which is what the full-subduction engine
(`src/irrep/subduction.rs`, task 3+) has to reproduce from
`SG_DATA_HALL[subgroup]` plus an affine embedding.

Evidence per (space group, irrep, direction) case, all from the pinned
`isotropy_subgroup/iso.zip`, one process per query in its own temporary working
directory (the program writes `iso.log` next to itself)::

    PAGE 1000
    SC 250
    SET I ALL OR 1
    VALUE PARENT <sg>
    VALUE IRREP <ml>
    VALUE DIRECTION <label>
    SHOW BASIS
    SHOW ORIGIN
    SHOW SIZE
    SHOW ELEMENTS
    DISPLAY ISOTROPY

`SHOW ELEMENTS` and `SHOW SIZE` are only *printed* by `DISPLAY ISOTROPY`, and
`SHOW ELEMENTS` additionally requires a preceding `VALUE DIRECTION <label>`
(without it the column stays empty).  Both facts are pinned here so a future
refactor cannot silently turn this gate into a no-op.

Two independent reference tables come from the same archive:

* the operation labels (`C4x+`, `SGda`, `C21''`, `SGv1`, ...) are decoded with
  the program's own legend in `data_space.txt`: `point_op_label` (72 labels),
  `point_op_label_stokes` (the same 72 operations in axis notation such as
  `2[110]`, `-4[001]`) and `ipoint_op` (72 row-major integer rotation
  matrices).  No naming grammar is re-invented here, and a label that is not in
  the legend is a hard failure rather than a guess.
* the expected number of elements is the subgroup's **point group order**, read
  from `ispace_point_group` + `ipoint_group_order`.  For a centred subgroup that
  is strictly smaller than the number of operations listed in its conventional
  cell (e.g. #12 `C2/m` prints 4 coset representatives, not 8), which is exactly
  the distinction `docs/subduction-conventions.md` §4 warns about.

Translations are kept as exact rationals and stored **verbatim**: the program
prints unreduced representatives such as `(C2x|0,2,2)` and `(I|5/2,5/2,5/2)`.

Usage::

    python3 scripts/verify_isotropy_operations.py            # replay + verify
    python3 scripts/verify_isotropy_operations.py --write    # refresh fixture

The fixture (`tests/data/isotropy/operations.json`) is deliberately small: it
holds the pinned expectations the Rust tests of tasks 4-5 consume; it is not a
copy of the isotropy tables.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir))
ISO_DIR = os.path.join(ROOT, "isotropy_subgroup")
DATA_SPACE = os.path.join(ISO_DIR, "data_space.txt")
ARCHIVE = os.path.join(ISO_DIR, "iso.zip")
FIXTURE = os.path.join(ROOT, "tests", "data", "isotropy", "operations.json")
# The same data in a line format the Rust tests parse without a JSON
# dependency (the workspace lockfile lives outside this crate, so adding one is
# not an option here).
TEXT_FIXTURE = os.path.join(ROOT, "tests", "data", "isotropy", "operations.txt")

# The geometry gate already owns the parent-primitive-basis conversion and the
# machine-table reader, so the basis/origin half of each row is checked with
# those helpers instead of a second copy of the same arithmetic.  The import
# also pins that this gate cannot silently disagree with the geometry gate.
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import verify_isotropy_oracle as geometry_gate  # noqa: E402

BANNER = "Current setting is International (new ed.) with conventional basis vectors."
# Denominator of the shipped operation tables; the Rust layer pins the same
# value as `irrep::subduction::SOURCE_TRANSLATION_GRID`.
SOURCE_TRANSLATION_GRID = 12
SETTING_COMMANDS = ["PAGE 1000", "SC 250", "SET I ALL OR 1"]
LABEL_SECTIONS = ("point_op_label", "point_op_label_stokes", "ipoint_op")

# One entry per queried direction.  `elements` is the point group order of the
# subgroup; `size` is the isotropy-table `Size` (index of the subgroup lattice
# in the parent lattice) and is `> 1` for the supercell cases.  The batch spans
# P/C/F/I/R parents, several centrings of the subgroup, two records with the
# same stored `(W, origin)` but different point groups (#123 vs #47), and a
# non-zero origin (#139).
CASES = [
    {"sg": 221, "ml": "GM4+", "direction": "P1", "subgroup": 83, "size": 1, "elements": 8},
    {"sg": 221, "ml": "GM4+", "direction": "P2", "subgroup": 12, "size": 1, "elements": 4},
    {"sg": 221, "ml": "GM4+", "direction": "P3", "subgroup": 148, "size": 1, "elements": 6},
    {"sg": 221, "ml": "GM3+", "direction": "P1", "subgroup": 123, "size": 1, "elements": 16},
    {"sg": 221, "ml": "GM3+", "direction": "C1", "subgroup": 47, "size": 1, "elements": 8},
    {"sg": 225, "ml": "GM4-", "direction": "C1", "subgroup": 8, "size": 1, "elements": 2},
    {"sg": 225, "ml": "GM4-", "direction": "C2", "subgroup": 8, "size": 1, "elements": 2},
    {"sg": 16, "ml": "R1", "direction": "P1", "subgroup": 22, "size": 2, "elements": 4},
    {"sg": 167, "ml": "GM3+", "direction": "P1", "subgroup": 15, "size": 1, "elements": 4},
    {"sg": 139, "ml": "M1-", "direction": "P1", "subgroup": 126, "size": 2, "elements": 16},
    # Self-restrictions pin the source labels for compound full-star tests:
    # quaternionic one/two-arm seeds, disjoint k/-k stars, and distinct CIRs.
    {"sg": 19, "ml": "GM1", "direction": "P1", "subgroup": 19, "size": 1, "elements": 4},
    {"sg": 23, "ml": "GM1", "direction": "P1", "subgroup": 23, "size": 1, "elements": 4},
    {"sg": 45, "ml": "GM1", "direction": "P1", "subgroup": 45, "size": 1, "elements": 4},
    {"sg": 83, "ml": "GM1+", "direction": "P1", "subgroup": 83, "size": 1, "elements": 8},
]

ROW_RE = re.compile(
    r"^(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+"
    r"\(([^()]*)\),\(([^()]*)\),\(([^()]*)\)\s+\(([^()]*)\)\s+(.+)$"
)
NEW_ROW_RE = re.compile(r"^\d+\s+\S+\s+\d+\s+\S+\s+\(")
ELEMENT_TOKEN_RE = re.compile(
    r"^[A-Za-z0-9'+\-]+\|(?:-?\d+)(?:/-?\d+)?(?:,(?:-?\d+)(?:/-?\d+)?){2}$"
)
IDENTITY = (1, 0, 0, 0, 1, 0, 0, 0, 1)


def _sections(path):
    """Split an ISOTROPY data file into ``key -> list of raw lines``."""
    with open(path, encoding="latin-1") as handle:
        lines = handle.read().splitlines()
    sections = {}
    key = None
    for line in lines:
        if re.match(r"^[a-z_0-9]+$", line):
            key = line
            sections[key] = []
        elif key is not None:
            sections[key].append(line)
    return sections


def _quoted(section):
    # The archive pads the fixed-width label columns with spaces; the printed
    # `SHOW ELEMENTS` labels carry no padding, so strip before matching.
    return [token.strip() for token in re.findall(r'"([^"]*)"', "\n".join(section))]


def _flat_ints(section):
    return [int(token) for line in section for token in line.split()]


def determinant3(matrix):
    return (
        matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7])
        - matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6])
        + matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6])
    )


def multiply3(left, right):
    return tuple(
        sum(left[3 * row + t] * right[3 * t + col] for t in range(3))
        for row in range(3)
        for col in range(3)
    )


# The legend stores **row-action** matrices: a printed op `M` acts as
# `x' = x M` (the Stokes notation beside it confirms this, see
# `LEGEND_ANCHORS`).  The engine and `crate::SymmetryOps` use the column
# convention `x' = R x`, so decoding transposes the legend matrix.  Without
# that transpose the 167 GM3+ P1 row (#15 C2/c) is visibly wrong: its printed
# rotations do not preserve the subgroup lattice, while their transposes do.
LEGEND_ANCHORS = [
    # (label, axis to test, expected image of that axis under the row action)
    ("C2a", (1, 1, 0), (1, 1, 0)),      # 2[110]: the axis is fixed
    ("C2b", (1, -1, 0), (1, -1, 0)),    # 2[-110]
    ("C4z+", (1, 0, 0), (0, 1, 0)),     # 4[001]: right-handed +90 degrees
    ("C4z-", (1, 0, 0), (0, -1, 0)),    # 4[00-1]
    ("SGda", (1, 1, 0), (-1, -1, 0)),   # -2[110]: the normal is reversed
    ("I", (1, 2, 3), (-1, -2, -3)),
    ("C21''", (1, 0, 0), (1, 0, 0)),    # hexagonal 2[100]
    ("SGv1", (1, 0, 0), (-1, 0, 0)),    # hexagonal -2[100]
]


def row_action(rotation, axis):
    """``axis . M`` for a row-major 3x3 matrix."""
    return tuple(sum(axis[row] * rotation[3 * row + column] for row in range(3)) for column in range(3))


def load_point_op_legend(path):
    """The program's own operation-label legend, straight from the archive.

    Returns ``{label: {"rotation": (9 ints), "stokes": str, "index": 1-based}}``
    with the matrices exactly as stored (row action, see `LEGEND_ANCHORS`).
    A label that appears twice must carry the same rotation (the cubic and
    hexagonal blocks both start with `E`); a conflicting duplicate is an error,
    because matching it by name alone would then be ambiguous.
    """
    sections = _sections(path)
    for key in LABEL_SECTIONS:
        if key not in sections:
            raise RuntimeError(f"{path}: missing section {key!r}")
    labels = _quoted(sections["point_op_label"])
    stokes = _quoted(sections["point_op_label_stokes"])
    numbers = _flat_ints(sections["ipoint_op"])
    if len(labels) != 72 or len(stokes) != 72:
        raise RuntimeError(f"{path}: expected 72 labels, got {len(labels)}/{len(stokes)}")
    if len(numbers) != 9 * len(labels):
        raise RuntimeError(
            f"{path}: expected {9 * len(labels)} rotation entries, got {len(numbers)}"
        )
    legend = {}
    for index, label in enumerate(labels):
        rotation = tuple(numbers[9 * index:9 * index + 9])
        if abs(determinant3(rotation)) != 1:
            raise RuntimeError(f"{path}: label {label!r} is not unimodular")
        entry = {"rotation": rotation, "stokes": stokes[index], "index": index + 1}
        previous = legend.get(label)
        if previous is not None and previous["rotation"] != rotation:
            raise RuntimeError(
                f"{path}: label {label!r} has two different rotations "
                f"({previous['stokes']} and {stokes[index]})"
            )
        legend.setdefault(label, entry)
    if legend["E"]["rotation"] != IDENTITY:
        raise RuntimeError(f"{path}: label 'E' is not the identity")
    if determinant3(legend["I"]["rotation"]) != -1:
        raise RuntimeError(f"{path}: label 'I' is not an inversion")
    # Pin the action convention against the Stokes notation shipped next to the
    # labels, so "the legend is row action" is a checked claim rather than an
    # assumption baked into the decoder.
    for label, axis, expected in LEGEND_ANCHORS:
        entry = legend.get(label)
        if entry is None:
            raise RuntimeError(f"{path}: anchor label {label!r} is missing")
        image = row_action(entry["rotation"], axis)
        if image != expected:
            raise RuntimeError(
                f"{path}: label {label!r} ({entry['stokes']}) maps {axis} to "
                f"{image} under the row action, expected {expected}: the legend "
                "is not in the assumed convention"
            )
    return legend


def load_point_group_orders(path):
    """``sg -> point group order`` from ``ispace_point_group``/``ipoint_group_order``."""
    sections = _sections(path)
    per_space_group = _flat_ints(sections["ispace_point_group"])
    orders = _flat_ints(sections["ipoint_group_order"])
    if len(per_space_group) != 230:
        raise RuntimeError(f"{path}: expected 230 point group indices")
    if len(orders) != 32:
        raise RuntimeError(f"{path}: expected 32 point group orders")
    table = {}
    for sg, index in enumerate(per_space_group, start=1):
        if not 1 <= index <= len(orders):
            raise RuntimeError(f"{path}: SG {sg} has out-of-range point group {index}")
        table[sg] = orders[index - 1]
    return table


def parse_rational(token):
    if not re.match(r"^-?\d+(?:/-?\d+)?$", token):
        raise RuntimeError(f"not a rational: {token!r}")
    return Fraction(token)


def parse_elements(text, legend):
    """Decode a `SHOW ELEMENTS` cell into exact ``(label, rotation, translation)``.

    The cell is a comma-separated list of ``(label|t1,t2,t3)`` tokens and may be
    wrapped over several physical lines by the program's column layout; callers
    join those lines before calling this.  Every label must exist in the
    legend (`no fuzzy matching`) and no operation may repeat.
    """
    leftover = re.sub(r"\([^()]*\)", "", text)
    if leftover.strip(" ,\t"):
        raise RuntimeError(f"unparsable ELEMENTS cell: {text[:120]!r}")
    tokens = [token.strip() for token in re.findall(r"\(([^()]*)\)", text)]
    if not tokens:
        raise RuntimeError("empty ELEMENTS cell")
    elements = []
    seen = set()
    for token in tokens:
        if not ELEMENT_TOKEN_RE.match(token):
            raise RuntimeError(f"malformed element token: {token!r}")
        label, translations = token.split("|", 1)
        entry = legend.get(label)
        if entry is None:
            raise RuntimeError(f"unknown operation label: {label!r}")
        # Row action (printed) -> column action (engine), see LEGEND_ANCHORS.
        rotation = [
            [entry["rotation"][3 * column + row] for column in range(3)]
            for row in range(3)
        ]
        translation = tuple(parse_rational(part) for part in translations.split(","))
        if len(translation) != 3:
            raise RuntimeError(f"element {token!r} does not have three components")
        key = (label, translation)
        if key in seen:
            raise RuntimeError(f"duplicate element: {token!r}")
        seen.add(key)
        elements.append(
            {
                "label": label,
                "rotation": [value for row in rotation for value in row],
                "translation": [str(value) for value in translation],
                "raw": f"({token})",
            }
        )
    return elements


def inverse3(matrix):
    """Exact inverse of a row-major 3x3 matrix of Fractions."""
    determ = determinant3([value for row in matrix for value in row])
    if determ == 0:
        raise RuntimeError("singular matrix in the operation fixture")
    cofactors = [[None] * 3 for _ in range(3)]
    for row in range(3):
        for column in range(3):
            rows = [index for index in range(3) if index != row]
            columns = [index for index in range(3) if index != column]
            minor = (
                matrix[rows[0]][columns[0]] * matrix[rows[1]][columns[1]]
                - matrix[rows[0]][columns[1]] * matrix[rows[1]][columns[0]]
            )
            sign = 1 if (row + column) % 2 == 0 else -1
            cofactors[column][row] = sign * minor / determ
    return cofactors


def matmul3(left, right):
    return [
        [sum(left[i][t] * right[t][j] for t in range(3)) for j in range(3)]
        for i in range(3)
    ]


def matvec3(matrix, vector):
    return [sum(matrix[i][t] * vector[t] for t in range(3)) for i in range(3)]


def child_lattice(basis, centring):
    """Rows of the subgroup's primitive lattice, in parent conventional axes.

    ``centring`` holds the centring vectors in the subgroup's own conventional
    basis, so the primitive vectors are ``centring . basis`` (row ``j`` of the
    centring matrix holds the coefficients of the conventional rows).
    """
    return [
        [sum(centring[row][inner] * basis[inner][column] for inner in range(3)) for column in range(3)]
        for row in range(3)
    ]


def preserves_lattice(lattice, rotation):
    """Whether ``rotation`` maps the lattice (rows = basis vectors) to itself.

    A point is ``rows^T . coordinates``, so its coordinates are
    ``(rows^T)^-1 . point`` -- the *transposed* inverse.  Using ``rows^-1``
    instead agrees only for symmetric bases.
    """
    coordinate_matrix = [
        [inverse3(lattice)[row][column] for row in range(3)] for column in range(3)
    ]
    for row in lattice:
        image = [
            sum(rotation[i][t] * row[t] for t in range(3)) for i in range(3)
        ]
        coordinates = [
            sum(coordinate_matrix[i][t] * image[t] for t in range(3)) for i in range(3)
        ]
        if any(value.denominator != 1 for value in coordinates):
            return False
    return True


def parse_row(text):
    """Parse the single printed row of a filtered `DISPLAY ISOTROPY` table.

    Continuation lines (the `Elements` column wraps) are joined first.  A
    truncated or paged table cannot produce a well-formed row, so this either
    returns the complete row or raises.
    """
    body = text.split("***", 1)[-1]
    if "Syntax error" in body:
        raise RuntimeError(f"oracle reported a syntax error: {body.strip()[:120]!r}")
    for marker in ("Press", "--More--", "more)"):
        if marker in body:
            raise RuntimeError("oracle output looks paged/truncated")
    lines = [line.rstrip() for line in body.splitlines()]
    header = next(
        (
            i
            for i, line in enumerate(lines)
            if "Subgroup" in line and "Dir" in line and "Elements" in line
        ),
        None,
    )
    if header is None:
        raise RuntimeError(
            "no `Subgroup ... Dir ... Elements` table header in the oracle output; "
            "SHOW SUBGROUP/SHOW DIRECTION frame the columns"
        )
    collected = []
    for line in lines[header + 1:]:
        stripped = line.strip()
        if stripped.startswith("*") or stripped.startswith("**"):
            break
        if not stripped:
            continue
        if collected and NEW_ROW_RE.match(stripped):
            raise RuntimeError(
                "the oracle printed more than one direction row; "
                "VALUE DIRECTION did not filter the table"
            )
        collected.append(stripped)
    if not collected:
        raise RuntimeError("the oracle printed no direction row")
    row = " ".join(collected)
    match = ROW_RE.match(row)
    if match is None:
        raise RuntimeError(f"unparsable isotropy row: {row[:160]!r}")
    basis = [match.group(5), match.group(6), match.group(7)]
    return {
        "subgroup": int(match.group(1)),
        "symbol": match.group(2),
        "size": int(match.group(3)),
        "direction": match.group(4),
        "basis": [[str(parse_rational(v)) for v in vector.split(",")] for vector in basis],
        "origin": [str(parse_rational(v)) for v in match.group(8).split(",")],
        "elements_text": match.group(9),
    }


def run_oracle(sg, ml, direction):
    """Run one query in its own temporary working directory and parse the row."""
    binary = os.path.join(ISO_DIR, "iso")
    if not os.path.isfile(binary):
        raise FileNotFoundError(
            "the ISOTROPY binary is not extracted; run\n"
            f"  unzip -o {ARCHIVE} -d {ISO_DIR}\n"
            "before using this oracle"
        )
    commands = SETTING_COMMANDS + [
        f"VALUE PARENT {sg}",
        f"VALUE IRREP {ml}",
        f"VALUE DIRECTION {direction}",
        # The Subgroup/Dir columns only appear when these are requested, and
        # `SHOW ELEMENTS` stays empty without the preceding VALUE DIRECTION.
        "SHOW SUBGROUP",
        "SHOW DIRECTION",
        "SHOW BASIS",
        "SHOW ORIGIN",
        "SHOW SIZE",
        "SHOW ELEMENTS",
        "DISPLAY ISOTROPY",
        "QUIT",
    ]
    env = dict(os.environ, ISODATA=ISO_DIR + os.sep)
    with tempfile.TemporaryDirectory() as workdir:
        result = subprocess.run(
            [binary],
            input="\n".join(commands) + "\n",
            capture_output=True,
            text=True,
            cwd=workdir,
            env=env,
            timeout=300,
        )
    if result.returncode != 0:
        raise RuntimeError(f"iso exited with {result.returncode}: {result.stderr}")
    if BANNER not in result.stdout:
        raise RuntimeError(f"iso did not report the pinned setting: {result.stdout[:200]!r}")
    row = parse_row(result.stdout)
    if row["direction"] != direction:
        raise RuntimeError(
            f"oracle printed direction {row['direction']!r}, asked for {direction!r}"
        )
    return row


def check_case(case, row, legend, orders, machine_row):
    """Verify one row against the pinned expectation and the archive tables.

    Two independent halves: the *elements* are decoded with the program's own
    label legend and must form the subgroup's point group, while the
    *basis/origin/size* half is compared against the stored isotropy tables
    through the geometry gate's conversion helpers (which the geometry gate
    pins on 48 other rows).
    """
    where = f"SG {case['sg']} {case['ml']} {case['direction']}"
    if row["subgroup"] != case["subgroup"]:
        raise RuntimeError(f"{where}: subgroup {row['subgroup']} != {case['subgroup']}")
    if row["size"] != case["size"]:
        raise RuntimeError(f"{where}: size {row['size']} != {case['size']}")
    if machine_row is None:
        raise RuntimeError(f"{where}: the stored isotropy tables have no such record")
    if machine_row["subgroup"] != case["subgroup"] or machine_row["size"] != case["size"]:
        raise RuntimeError(
            f"{where}: machine record is #({machine_row['subgroup']}, "
            f"size {machine_row['size']}), pinned #({case['subgroup']}, {case['size']})"
        )

    parent_basis = geometry_gate.PRIMITIVE_BASIS[
        geometry_gate.CENTERING_LETTER[case["sg"]]
    ]
    oracle_basis = [[Fraction(value) for value in vector] for vector in row["basis"]]
    oracle_origin = [Fraction(value) for value in row["origin"]]
    converted_basis = geometry_gate.to_conventional_basis(machine_row["basis"], parent_basis)
    sub_centring = geometry_gate.PRIMITIVE_BASIS[
        geometry_gate.CENTERING_LETTER[case["subgroup"]]
    ]
    printed_lattice = geometry_gate.primitive_of_conventional_cell(oracle_basis, sub_centring)
    if not geometry_gate.same_lattice(converted_basis, printed_lattice):
        raise RuntimeError(
            f"{where}: stored basis {converted_basis} and printed cell "
            f"{oracle_basis} (times the subgroup centring) span different lattices"
        )
    z_parent = geometry_gate.CENTERING_Z[geometry_gate.CENTERING_LETTER[case["sg"]]]
    z_sub = geometry_gate.CENTERING_Z[geometry_gate.CENTERING_LETTER[case["subgroup"]]]
    expected_volume = Fraction(z_sub * machine_row["size"], z_parent)
    volume = abs(geometry_gate.det3(oracle_basis))
    if volume != expected_volume:
        raise RuntimeError(
            f"{where}: |det Basis| {volume} != Z_sub {z_sub} * size "
            f"{machine_row['size']} / Z_parent {z_parent}"
        )
    converted_origin = geometry_gate.to_conventional(machine_row["origin"], parent_basis)
    if converted_origin != oracle_origin:
        raise RuntimeError(
            f"{where}: origin {oracle_origin} != stored {converted_origin} exactly"
        )

    centring = geometry_gate.CENTERING_LETTER[case["subgroup"]]
    elements = parse_elements(row.pop("elements_text"), legend)
    if len(elements) != case["elements"]:
        raise RuntimeError(
            f"{where}: printed {len(elements)} elements, pinned {case['elements']}"
        )
    expected_order = orders[case["subgroup"]]
    if case["elements"] != expected_order:
        raise RuntimeError(
            f"{where}: {case['elements']} representatives for point group order "
            f"{expected_order} of SG {case['subgroup']}"
        )
    lattice = child_lattice(
        [[Fraction(value) for value in vector] for vector in row["basis"]],
        [[Fraction(value) for value in vector]
         for vector in geometry_gate.PRIMITIVE_BASIS[centring]],
    )
    for element in elements:
        rotation = [
            [Fraction(element["rotation"][3 * i + j]) for j in range(3)] for i in range(3)
        ]
        if not preserves_lattice(lattice, rotation):
            raise RuntimeError(
                f"{where}: {element['raw']} does not preserve the subgroup "
                "lattice -- the rotation is not in the printed basis' frame"
            )
    # The printed row must invert to a genuine operation of the subgroup in the
    # subgroup's own conventional basis: an integer rotation that preserves the
    # subgroup cell, and a translation on the shipped source grid.  This is the
    # strongest check of the whole convention chain (T = B^T, column action,
    # t_G = T t_H + o - R_G o) and it is what makes the row-action decode of the
    # 167 GM3+ P1 row verifiable rather than assumed.
    transform = [
        [Fraction(row["basis"][j][i]) for j in range(3)] for i in range(3)
    ]
    transform_inverse = inverse3(transform)
    origin_vector = [Fraction(value) for value in row["origin"]]
    child_cell = [
        [Fraction(value) for value in vector]
        for vector in geometry_gate.PRIMITIVE_BASIS[centring]
    ]
    for element in elements:
        rotation = [
            [Fraction(element["rotation"][3 * i + j]) for j in range(3)] for i in range(3)
        ]
        translation = [Fraction(value) for value in element["translation"]]
        child_rotation = matmul3(matmul3(transform_inverse, rotation), transform)
        if any(value.denominator != 1 for vector in child_rotation for value in vector):
            raise RuntimeError(
                f"{where}: {element['raw']} does not invert to an integer rotation "
                "in the subgroup basis"
            )
        if not preserves_lattice(child_cell, child_rotation):
            raise RuntimeError(
                f"{where}: {element['raw']} does not preserve the subgroup cell "
                "lattice in the subgroup basis"
            )
        rotated_origin = matvec3(rotation, origin_vector)
        child_translation = matvec3(
            transform_inverse,
            [
                translation[i] + rotated_origin[i] - origin_vector[i]
                for i in range(3)
            ],
        )
        for value in child_translation:
            if SOURCE_TRANSLATION_GRID % value.denominator != 0:
                raise RuntimeError(
                    f"{where}: {element['raw']} inverts to translation "
                    f"{[str(v) for v in child_translation]}, which is off the "
                    f"1/{SOURCE_TRANSLATION_GRID} source grid"
                )
    rotations = [tuple(element["rotation"]) for element in elements]
    if len(set(rotations)) != len(rotations):
        raise RuntimeError(f"{where}: the same rotation appears twice")
    unique = set(rotations)
    if IDENTITY not in unique:
        raise RuntimeError(f"{where}: no identity operation")
    for rotation in unique:
        if not _is_inverse_pair(unique, rotation):
            raise RuntimeError(f"{where}: {rotation} has no inverse in the set")
        for other in unique:
            if multiply3(rotation, other) not in unique:
                raise RuntimeError(
                    f"{where}: rotations are not closed ({rotation} * {other})"
                )
    identity = next(element for element in elements if element["label"] == "E")
    if identity["rotation"] != list(IDENTITY):
        raise RuntimeError(f"{where}: label 'E' is not the identity rotation")
    if any(Fraction(value) != 0 for value in identity["translation"]):
        raise RuntimeError(f"{where}: identity carries translation {identity['translation']}")
    row["elements"] = elements
    return row


def _is_inverse_pair(unique, rotation):
    return any(multiply3(rotation, other) == IDENTITY for other in unique)


def build_fixture():
    legend = load_point_op_legend(DATA_SPACE)
    orders = load_point_group_orders(DATA_SPACE)
    records = geometry_gate.machine_records()
    cases = []
    element_total = 0
    for case in CASES:
        by_label = {record["label"]: record for record in records.get((case["sg"], case["ml"]), [])}
        row = run_oracle(case["sg"], case["ml"], case["direction"])
        row = check_case(case, row, legend, orders, by_label.get(case["direction"]))
        element_total += len(row["elements"])
        centring = geometry_gate.CENTERING_LETTER[case["subgroup"]]
        order = orders[case["subgroup"]]
        conventional = order * geometry_gate.CENTERING_Z[centring]
        if centring != "P" and conventional == len(row["elements"]):
            raise RuntimeError(
                f"SG {case['subgroup']}: centred cell must list more operations "
                f"than the {order} point group cosets it prints"
            )
        row.update(
            {
                "sg": case["sg"],
                "ml": case["ml"],
                "centring": centring,
                "point_group_order": order,
                "conventional_cell_operations": conventional,
            }
        )
        cases.append(row)
    provenance = {
        "oracle": "isotropy_subgroup/iso.zip",
        "oracle_banner": "Isotropy, Version 9.6.1, Jan 2022",
        "setting_commands": SETTING_COMMANDS,
        "element_columns": ["SHOW BASIS", "SHOW ORIGIN", "SHOW SIZE", "SHOW ELEMENTS"],
        "element_display": "DISPLAY ISOTROPY (with VALUE DIRECTION <label>)",
        "rotation_convention": (
            "column action x' = R x, i.e. the engine's convention.  The program "
            "and data_space.txt store row-action matrices (x' = x M); decoding "
            "transposes them, which the lattice-preservation check verifies on "
            "every element (and which is what exposes 167 GM3+ P1 otherwise)"
        ),
        "label_legend": {
            "file": "data_space.txt",
            "sections": list(LABEL_SECTIONS),
            "slots": 72,
            "unique_labels": len(legend),
            "sha256": legend_digest(legend),
        },
        "point_group_orders": {
            "file": "data_space.txt",
            "sections": ["ispace_point_group", "ipoint_group_order"],
        },
        "convention_cross_check": (
            "every operation is mapped back into the subgroup basis and must "
            "give an integer rotation that preserves the subgroup cell and a "
            "translation on the 1/12 source grid"
        ),
        "geometry_cross_check": (
            "each basis/origin is compared against the stored isotropy tables "
            "through scripts/verify_isotropy_oracle.py (lattice equality, "
            "Z_sub*size/Z_parent volume, exact origin)"
        ),
        "archive_sha256": file_digest(ARCHIVE),
        "note": (
            "elements are coset representatives of the subgroup's translation "
            "subgroup in the parent frame; their number is the subgroup point "
            "group order, not the conventional-cell operation count "
            "(`conventional_cell_operations` = order * Z_subgroup)"
        ),
    }
    return {"provenance": provenance, "cases": cases}, element_total


def render_text_fixture(fixture):
    """The compact line format read by `tests/irrep_subduction.rs`."""
    lines = [
        "# generated by scripts/verify_isotropy_operations.py --write",
        "# case <sg> <ml> <direction> <subgroup> <size> <centring> <point-group-order>",
        "# basis <a;b;c>   origin <o>   element <label> <R> <t>",
    ]
    for case in fixture["cases"]:
        lines.append(
            "case {} {} {} {} {} {} {}".format(
                case["sg"],
                case["ml"],
                case["direction"],
                case["subgroup"],
                case["size"],
                case["centring"],
                case["point_group_order"],
            )
        )
        lines.append("basis " + ";".join(",".join(row) for row in case["basis"]))
        lines.append("origin " + ",".join(case["origin"]))
        for element in case["elements"]:
            rotation = ";".join(
                ",".join(str(value) for value in element["rotation"][3 * i:3 * i + 3])
                for i in range(3)
            )
            lines.append(
                "element {} {} {}".format(
                    element["label"], rotation, ",".join(element["translation"])
                )
            )
    return "\n".join(lines) + "\n"


def legend_digest(legend):
    payload = json.dumps(
        {label: entry["rotation"] for label, entry in sorted(legend.items())},
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_fixture(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--write", action="store_true", help="rewrite the fixture instead of replaying it"
    )
    args = parser.parse_args(argv)
    if not os.path.isfile(DATA_SPACE):
        print(
            f"missing {DATA_SPACE}; run\n  unzip -o {ARCHIVE} -d {ISO_DIR}",
            file=sys.stderr,
        )
        return 2
    fixture, element_total = build_fixture()
    if args.write:
        os.makedirs(os.path.dirname(FIXTURE), exist_ok=True)
        with open(FIXTURE, "w", encoding="utf-8") as handle:
            json.dump(fixture, handle, indent=2, sort_keys=False)
            handle.write("\n")
        with open(TEXT_FIXTURE, "w", encoding="utf-8") as handle:
            handle.write(render_text_fixture(fixture))
        print(
            f"wrote {os.path.relpath(FIXTURE, ROOT)} and "
            f"{os.path.relpath(TEXT_FIXTURE, ROOT)}"
        )
    else:
        stored = load_fixture(FIXTURE)
        if stored != fixture:
            print("the stored operation fixture differs from the oracle:", file=sys.stderr)
            _report_difference(stored, fixture)
            return 1
        with open(TEXT_FIXTURE, encoding="utf-8") as handle:
            text = handle.read()
        if text != render_text_fixture(fixture):
            print(
                f"{os.path.relpath(TEXT_FIXTURE, ROOT)} is out of date with the "
                "oracle; rerun with --write",
                file=sys.stderr,
            )
            return 1
    print(
        f"operation fixtures checked: {len(fixture['cases'])} cases / "
        f"{element_total} coset representatives, all exact "
        f"(settings pinned, origins and translations verbatim)"
    )
    return 0


def _report_difference(stored, fresh):
    stored_cases = {(c["sg"], c["ml"], c["direction"]): c for c in stored.get("cases", [])}
    fresh_cases = {(c["sg"], c["ml"], c["direction"]): c for c in fresh.get("cases", [])}
    for key in sorted(fresh_cases.keys() - stored_cases.keys()):
        print(f"  missing from fixture: {key}", file=sys.stderr)
    for key in sorted(stored_cases.keys() - fresh_cases.keys()):
        print(f"  unexpected in fixture: {key}", file=sys.stderr)
    for key in sorted(stored_cases.keys() & fresh_cases.keys()):
        if stored_cases[key] != fresh_cases[key]:
            print(f"  {key}: fixture {stored_cases[key]}", file=sys.stderr)
            print(f"  {key}: oracle  {fresh_cases[key]}", file=sys.stderr)
    if stored.get("provenance") != fresh.get("provenance"):
        print("  provenance block differs", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
