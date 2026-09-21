#!/usr/bin/env python3
"""Cross-check the vendored isotropy geometry against the ISO `iso` program.

`isotropy_subgroup/iso.zip` ships the ISOTROPY suite binaries together with the
pinned data files that `generate_irrep_data.py` parses.  Running `iso` therefore
gives an independent oracle for the semantics that cannot be inferred from the
raw arrays alone.

For every requested (space group, Miller-Love irrep) pair the oracle prints one
row per order-parameter direction::

    <subgroup#> <symbol> <Size> <Dir> <Basis> <Origin>

The stored arrays and the oracle use *different frames*, which is the main thing
this script pins down:

* ``isotropy_basis`` / ``isotropy_origin`` are expressed in the parent's
  **primitive-cell frame** (the frame the ISOTROPY program works in internally),
* the oracle prints the subgroup's ITA **conventional** basis and the origin in
  the parent's conventional basis.

Agreement is therefore checked after converting the stored values with the
parent's primitive basis ``P`` (rows = primitive vectors in conventional
coordinates):

* ``Size`` equals ``|det W|``,
* ``|det Basis_oracle|`` equals ``Z(subgroup) * |det W| / Z(parent)``,
* ``w_machine . P`` equals ``w_oracle`` exactly,
* subgroup number and direction label are identical,
* dimension-3 descriptors agree under the Rust API's whitespace/separator
  normalization; missing or duplicate vector labels fail the gate.

The **origin** needs one more qualifier.  The program prints its ``Origin``
column in whatever ITA setting is currently selected (``SET I``), and its
factory default is *origin choice 2 for every space group*, while the pinned
data tables were recorded in a mixed setting: origin choice 1 for most parents
and origin choice 2 for the space groups 227/228 class.  This script therefore
runs ``SET I ALL OR 1`` so the oracle prints the setting the tables were
recorded in; without it roughly 18.6% of the whole table appears to disagree
(for SG 139 ``M1-`` ``P1`` the default prints ``(1/4,1/4,1/4)`` while the stored
``(2,2,2)`` is ``(1,1,1)``, which is exactly what OR 1 prints).  A residual
class of rows whose *cell/axis* choice differs (monoclinic, trigonal) is still
not reachable through ``SET I``; see ``docs/isotropy-data-semantics.md`` §3.

Usage::

    python3 scripts/verify_isotropy_oracle.py
"""

import os
import re
import subprocess
import sys
from fractions import Fraction

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
ISO_DIR = os.path.join(REPO_DIR, "isotropy_subgroup")

sys.path.insert(0, SCRIPT_DIR)
from generate_irrep_data import (  # noqa: E402
    read_file,
    get_sections,
    parse_ints,
    parse_labels,
    direction_str,
)

# Parent primitive basis per centering, in parent conventional coordinates.
# Pinning: the P1 rows printed by the oracle (SG 225 F, 229 I, 5/8 C, 167 R)
# reproduce these matrices, and the origin conversion below holds for every
# sampled space group.
PRIMITIVE_BASIS = {
    "P": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "I": [
        [Fraction(-1, 2), Fraction(1, 2), Fraction(1, 2)],
        [Fraction(1, 2), Fraction(-1, 2), Fraction(1, 2)],
        [Fraction(1, 2), Fraction(1, 2), Fraction(-1, 2)],
    ],
    "F": [
        [0, Fraction(1, 2), Fraction(1, 2)],
        [Fraction(1, 2), 0, Fraction(1, 2)],
        [Fraction(1, 2), Fraction(1, 2), 0],
    ],
    "A": [
        [0, Fraction(1, 2), Fraction(1, 2)],
        [0, Fraction(-1, 2), Fraction(1, 2)],
        [1, 0, 0],
    ],
    # No standard ITA setting used by the pinned data is B-centred, so this
    # entry is never exercised; the Rust table keeps it for completeness.
    "B": [
        [Fraction(1, 2), 0, Fraction(1, 2)],
        [Fraction(-1, 2), 0, Fraction(1, 2)],
        [0, 1, 0],
    ],
    "C": [
        [Fraction(1, 2), Fraction(1, 2), 0],
        [Fraction(-1, 2), Fraction(1, 2), 0],
        [0, 0, 1],
    ],
    # Rhombohedral lattice in hexagonal axes (ISOTROPY's default setting).
    "R": [
        [Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)],
        [Fraction(-1, 3), Fraction(1, 3), Fraction(1, 3)],
        [Fraction(-1, 3), Fraction(-2, 3), Fraction(1, 3)],
    ],
}

# Centering letter per space group (1-230).
CENTERING_LETTER = {
    1: "P", 2: "P", 3: "P", 4: "P", 5: "C", 6: "P", 7: "P", 8: "C", 9: "C",
    10: "P", 11: "P", 12: "C", 13: "P", 14: "P", 15: "C", 16: "P", 17: "P",
    18: "P", 19: "P", 20: "C", 21: "C", 22: "F", 23: "I", 24: "I", 25: "P",
    26: "P", 27: "P", 28: "P", 29: "P", 30: "P", 31: "P", 32: "P", 33: "P",
    34: "P", 35: "C", 36: "C", 37: "C", 38: "A", 39: "A", 40: "A", 41: "A",
    42: "F", 43: "F", 44: "I", 45: "I", 46: "I", 47: "P", 48: "P", 49: "P",
    50: "P", 51: "P", 52: "P", 53: "P", 54: "P", 55: "P", 56: "P", 57: "P",
    58: "P", 59: "P", 60: "P", 61: "P", 62: "P", 63: "C", 64: "C", 65: "C",
    66: "C", 67: "C", 68: "C", 69: "F", 70: "F", 71: "I", 72: "I", 73: "I",
    74: "I", 75: "P", 76: "P", 77: "P", 78: "P", 79: "I", 80: "I", 81: "P",
    82: "I", 83: "P", 84: "P", 85: "P", 86: "P", 87: "I", 88: "I", 89: "P",
    90: "P", 91: "P", 92: "P", 93: "P", 94: "P", 95: "P", 96: "P", 97: "I",
    98: "I", 99: "P", 100: "P", 101: "P", 102: "P", 103: "P", 104: "P",
    105: "P", 106: "P", 107: "I", 108: "I", 109: "I", 110: "I", 111: "P",
    112: "P", 113: "P", 114: "P", 115: "P", 116: "P", 117: "P", 118: "P",
    119: "I", 120: "I", 121: "I", 122: "I", 123: "P", 124: "P", 125: "P",
    126: "P", 127: "P", 128: "P", 129: "P", 130: "P", 131: "P", 132: "P",
    133: "P", 134: "P", 135: "P", 136: "P", 137: "P", 138: "P", 139: "I",
    140: "I", 141: "I", 142: "I", 143: "P", 144: "P", 145: "P", 146: "R",
    147: "P", 148: "R", 149: "P", 150: "P", 151: "P", 152: "P", 153: "P",
    154: "P", 155: "R", 156: "P", 157: "P", 158: "P", 159: "P", 160: "R",
    161: "R", 162: "P", 163: "P", 164: "P", 165: "P", 166: "R", 167: "R",
    168: "P", 169: "P", 170: "P", 171: "P", 172: "P", 173: "P", 174: "P",
    175: "P", 176: "P", 177: "P", 178: "P", 179: "P", 180: "P", 181: "P",
    182: "P", 183: "P", 184: "P", 185: "P", 186: "P", 187: "P", 188: "P",
    189: "P", 190: "P", 191: "P", 192: "P", 193: "P", 194: "P", 195: "P",
    196: "F", 197: "I", 198: "P", 199: "I", 200: "P", 201: "P", 202: "F",
    203: "F", 204: "I", 205: "P", 206: "I", 207: "P", 208: "P", 209: "F",
    210: "F", 211: "I", 212: "P", 213: "P", 214: "I", 215: "P", 216: "F",
    217: "I", 218: "P", 219: "F", 220: "I", 221: "P", 222: "P", 223: "P",
    224: "P", 225: "F", 226: "F", 227: "F", 228: "F", 229: "I", 230: "I",
}

CENTERING_Z = {"P": 1, "A": 2, "B": 2, "C": 2, "I": 2, "F": 4, "R": 3}

# (space group, Miller-Love irrep) pairs spanning every centering and
# 1-, 2-, 3- and 6-dimensional irreps.
CASES = [
    (2, "GM1-"),
    (5, "GM2"),
    (8, "GM2"),
    (12, "GM1+"),
    (14, "GM1+"),
    (15, "GM1+"),
    (16, "R1"),
    (16, "R2"),
    (38, "GM1"),
    (39, "GM1"),
    (64, "GM1+"),
    (123, "GM4+"),
    (139, "GM4+"),
    (167, "GM3+"),
    (177, "L1"),  # Complex separators and the non-cubic C2 component order.
    (194, "GM6+"),
    (213, "X2"),  # C23 -> R3: distinguish primitive and hexagonal cells.
    (221, "GM3+"),
    (221, "GM4+"),
    (225, "GM3+"),
    (225, "GM4-"),
    (229, "GM4-"),
    (230, "GM5+"),
]


def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def machine_records():
    """Isotropy records keyed by (parent SG, Miller-Love label)."""
    iso_lines = read_file("data_isotropy.txt")
    iso_sec = get_sections(iso_lines)
    ir_lines = read_file("data_irreps.txt")
    ir_sec = get_sections(ir_lines)

    labels = parse_labels(ir_lines, ir_sec, "irrep_label")
    sg_numbers = parse_ints(ir_lines, ir_sec, "irrep_space_group")
    assert len(sg_numbers) == len(labels)

    par = parse_ints(iso_lines, iso_sec, "isotropy_parent")
    irr = parse_ints(iso_lines, iso_sec, "isotropy_irrep")
    sub = parse_ints(iso_lines, iso_sec, "isotropy_subgroup")
    direction = parse_ints(iso_lines, iso_sec, "isotropy_direction")
    dim = parse_ints(iso_lines, iso_sec, "isotropy_orderparam_dim")
    free = parse_ints(iso_lines, iso_sec, "isotropy_orderparam_freeparam")
    label = parse_labels(iso_lines, iso_sec, "isotropy_orderparam_label")
    basis = parse_ints(iso_lines, iso_sec, "isotropy_basis")
    origin = parse_ints(iso_lines, iso_sec, "isotropy_origin")
    dir_text = [
        direction_str(dim[k], free[k], label[k], par[k]) for k in range(len(par))
    ]

    records = {}
    for k in range(len(par)):
        ml = labels[irr[k] - 1]
        b = basis[k * 9:(k + 1) * 9]
        o = origin[k * 4:(k + 1) * 4]
        matrix = [b[0:3], b[3:6], b[6:9]]
        records.setdefault((par[k], ml), []).append(
            {
                "subgroup": sub[k],
                "direction": dir_text[k],
                "label": label[k],
                "size": abs(det3(matrix)),
                "dim": dim[k],
                "basis": [[Fraction(value) for value in row] for row in matrix],
                "origin": [
                    Fraction(o[0], o[3]),
                    Fraction(o[1], o[3]),
                    Fraction(o[2], o[3]),
                ],
            }
        )
    return records


def run_oracle(sg, ml):
    """Run `iso` for one (space group, irrep) and return its printed rows."""
    binary = os.path.join(ISO_DIR, "iso")
    if not os.path.isfile(binary):
        raise FileNotFoundError(
            "the ISOTROPY binary is not extracted; run\n"
            f"  unzip -o {os.path.join(ISO_DIR, 'iso.zip')} -d {ISO_DIR}\n"
            "before using this oracle"
        )
    commands = [
        "PAGE 1000",
        "SC 250",
        # Print the Origin column in the ITA setting the data tables were
        # recorded in (origin choice 1); the factory default is origin choice 2
        # for every space group and disagrees for 18.6% of the whole table.
        "SET I ALL OR 1",
        f"VALUE PARENT {sg}",
        f"VALUE IRREP {ml}",
        "SHOW SUBGROUP",
        "SHOW BASIS",
        "SHOW ORIGIN",
        "SHOW SIZE",
        "SHOW DIRECTION",
        "DISPLAY ISOTROPY",
        "QUIT",
    ]
    env = dict(os.environ, ISODATA=ISO_DIR + os.sep)
    result = subprocess.run(
        [binary],
        input="\n".join(commands) + "\n",
        capture_output=True,
        text=True,
        cwd=ISO_DIR,
        env=env,
        timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(f"iso exited with {result.returncode}: {result.stderr}")
    rows = []
    for line in result.stdout.splitlines():
        line = line.strip().rstrip("*").strip()
        m = re.match(r"^(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)$", line)
        if not m:
            continue
        rows.append(
            {
                "subgroup": int(m.group(1)),
                "size": int(m.group(3)),
                "label": m.group(4),
                "basis": [
                    [Fraction(token) for token in vector.split(",")]
                    for vector in m.group(5).strip("()").split("),(")
                ],
                "origin": [
                    Fraction(token) for token in m.group(6).strip("()").split(",")
                ],
            }
        )
    return rows


def run_oracle_direction_vectors(sg, ml):
    """The program's own order-parameter component strings, keyed by direction.

    `SHOW DIRECTION VECTOR` is what validates the descriptor strings the Rust
    API accepts from users: the row labels (`P1`, `C1`, …) only pin the *rows*,
    not the component strings that `IsotropyDirection::Descriptor` matches.
    """
    binary = os.path.join(ISO_DIR, "iso")
    commands = [
        "PAGE 1000",
        "SC 250",
        "SET I ALL OR 1",
        f"VALUE PARENT {sg}",
        f"VALUE IRREP {ml}",
        "SHOW SUBGROUP",
        "SHOW DIRECTION VECTOR",
        "DISPLAY ISOTROPY",
        "QUIT",
    ]
    result = subprocess.run(
        [binary],
        input="\n".join(commands) + "\n",
        capture_output=True,
        text=True,
        cwd=ISO_DIR,
        env=dict(os.environ, ISODATA=ISO_DIR + os.sep),
        timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(f"iso exited with {result.returncode}: {result.stderr}")
    vectors = {}
    for line in result.stdout.splitlines():
        line = line.strip().rstrip("*").strip()
        m = re.match(r"^(\d+)\s+(\S+)\s+(\S+)\s+(\(.*\))$", line)
        if m:
            label = m.group(3)
            if label in vectors:
                raise RuntimeError(f"SG {sg} {ml}: duplicate direction vector label {label}")
            vectors[label] = "".join(m.group(4).split())
    return vectors


def same_lattice(left, right):
    """Whether two generating sets span the same lattice (rows, exact rationals).

    Neither set has to be a basis: the program prints the subgroup's
    *conventional* cell, whose lattice is generated by that cell together with
    the subgroup's own centring translations.
    """
    return all(in_lattice(row, right) for row in left) and all(
        in_lattice(row, left) for row in right
    )


def in_lattice(vector, primitive_basis):
    """Whether ``vector`` is an integer combination of ``primitive_basis`` rows."""
    matrix = [
        [Fraction(primitive_basis[i][j]) for i in range(3)] + [Fraction(vector[j])]
        for j in range(3)
    ]
    for column in range(3):
        pivot = next(
            (row for row in range(column, 3) if matrix[row][column] != 0), None
        )
        if pivot is None:
            return False
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for row in range(3):
            if row != column and matrix[row][column] != 0:
                factor = matrix[row][column]
                matrix[row] = [
                    matrix[row][j] - factor * matrix[column][j] for j in range(4)
                ]
    return all(matrix[i][3].denominator == 1 for i in range(3))


def to_conventional(vector, primitive_basis):
    """Express primitive-frame coordinates in the parent conventional basis."""
    return [
        sum(vector[i] * Fraction(primitive_basis[i][j]) for i in range(3))
        for j in range(3)
    ]


def to_conventional_basis(basis, primitive_basis):
    """Express a row-major basis in the parent conventional basis."""
    return [to_conventional(row, primitive_basis) for row in basis]


def primitive_of_conventional_cell(cell, centring):
    """Primitive vectors of a centred conventional cell, in the same frame.

    ``cell`` rows are the conventional vectors and ``centring`` rows are that
    lattice's primitive vectors *in its own conventional basis*, so the result
    is ``centring . cell``: row ``j`` of ``centring`` holds the coefficients of
    the conventional rows.  (Multiplying the other way round only agrees for
    symmetric centring matrices such as the F one.)
    """
    return [
        [
            sum(Fraction(centring[j][t]) * Fraction(cell[t][i]) for t in range(3))
            for i in range(3)
        ]
        for j in range(3)
    ]


def main():
    records = machine_records()
    checked_rows = 0
    checked_descriptors = 0
    exact_origins = 0
    lattice_origins = 0
    failures = []

    for sg, ml in CASES:
        expected = records.get((sg, ml))
        if not expected:
            failures.append(f"SG {sg} {ml}: no machine records")
            continue
        rows = run_oracle(sg, ml)
        if len(rows) != len(expected):
            failures.append(
                f"SG {sg} {ml}: oracle printed {len(rows)} rows, "
                f"machine has {len(expected)}"
            )
            continue
        # The component strings the user actually types must be the program's:
        # comparing only `P1`/`C1` labels would miss a wrong `(a,a,0)` mapping.
        vectors = run_oracle_direction_vectors(sg, ml)
        # Rows are matched by ISOTROPY direction label: neither the oracle nor
        # the data file guarantees the same row order.
        by_label = {row["label"]: row for row in rows}
        if len(by_label) != len(rows):
            failures.append(f"SG {sg} {ml}: oracle labels are not unique")
            continue
        if vectors.keys() != by_label.keys():
            failures.append(
                f"SG {sg} {ml}: direction vector labels differ from geometry "
                f"(missing {sorted(by_label.keys() - vectors.keys())}, "
                f"unexpected {sorted(vectors.keys() - by_label.keys())})"
            )
            continue
        primitive_basis = PRIMITIVE_BASIS[CENTERING_LETTER[sg]]
        for machine_row in expected:
            label = machine_row["label"]
            where = f"SG {sg} {ml} {label} {machine_row['direction']}"
            oracle_row = by_label.get(label)
            if oracle_row is None:
                failures.append(f"{where}: oracle printed no row for this label")
                continue
            checked_rows += 1
            if machine_row["dim"] == 3:
                # Only `dim = 3` descriptors are claimed to be the program's own
                # strings.  `dim = 2` is cryspglib's internal notation: the
                # program's spelling varies per irrep (SG 5 `L1` prints
                # `(a;a)`, SG 91 `A1` prints `(a,0)`, SG 194 `GM6+` prints
                # `(a,0.577a)`), and `dim >= 4` prints full component lists.
                checked_descriptors += 1
                # Match the public Descriptor selector: separators and
                # whitespace are aliases, but component order is significant.
                descriptor = "".join(machine_row["direction"].split()).replace(";", ",")
                if descriptor != vectors[label].replace(";", ","):
                    failures.append(
                        f"{where}: descriptor {machine_row['direction']!r} != "
                        f"oracle {vectors[label]!r}"
                    )
            if oracle_row["subgroup"] != machine_row["subgroup"]:
                failures.append(
                    f"{where}: subgroup {oracle_row['subgroup']} != "
                    f"{machine_row['subgroup']}"
                )
            if oracle_row["size"] != machine_row["size"]:
                failures.append(
                    f"{where}: size {oracle_row['size']} != {machine_row['size']}"
                )
            # The oracle prints the subgroup's conventional cell in parent
            # conventional coordinates, so its volume in parent conventional
            # units is Z(subgroup) * size / Z(parent).
            z_parent = CENTERING_Z[CENTERING_LETTER[sg]]
            z_sub = CENTERING_Z[CENTERING_LETTER[machine_row["subgroup"]]]
            oracle_det = abs(det3(oracle_row["basis"]))
            expected_det = Fraction(z_sub * machine_row["size"], z_parent)
            if oracle_det != expected_det:
                failures.append(
                    f"{where}: |det Basis| {oracle_det} != Z_sub {z_sub} * "
                    f"size {machine_row['size']} / Z_parent {z_parent}"
                )
            # A volume check cannot see a wrongly oriented parent primitive
            # basis (a transposed C matrix has the same determinant but spans a
            # different lattice), so compare the *lattices*: the stored basis
            # converted with the parent's primitive basis must generate the
            # lattice of the printed conventional cell together with the
            # subgroup's own centring translations.
            converted_basis = to_conventional_basis(machine_row["basis"], primitive_basis)
            sub_primitive = PRIMITIVE_BASIS[
                CENTERING_LETTER[machine_row["subgroup"]]
            ]
            printed_lattice = primitive_of_conventional_cell(
                oracle_row["basis"], sub_primitive
            )
            if not same_lattice(converted_basis, printed_lattice):
                failures.append(
                    f"{where}: converted basis {converted_basis} and printed cell "
                    f"{oracle_row['basis']} (times the subgroup centring) span "
                    f"different lattices"
                )
            # The stored origin is in the parent primitive frame; with the
            # oracle pinned to the recorded setting the two must agree
            # **exactly**.  A difference of a parent lattice vector is not an
            # acceptable answer here: it would let a wrong conversion (or a
            # mutation of the oracle output) pass, so it fails the gate and is
            # only reported as a diagnostic counter.
            converted = to_conventional(machine_row["origin"], primitive_basis)
            delta = [converted[i] - oracle_row["origin"][i] for i in range(3)]
            if all(value == 0 for value in delta):
                exact_origins += 1
            else:
                if in_lattice(delta, primitive_basis):
                    lattice_origins += 1
                failures.append(
                    f"{where}: origin {oracle_row['origin']} != {converted} exactly "
                    f"(delta {delta})"
                )

    print(f"oracle rows checked: {checked_rows}")
    print(f"descriptor strings checked: {checked_descriptors}")
    if failures:
        print(f"FAILURES: {len(failures)}")
        for failure in failures[:25]:
            print("  " + failure)
        return 1
    print(
        f"all oracle rows agree with the vendored isotropy data "
        f"(subgroup/size/basis/label on {checked_rows} rows, descriptors on "
        f"{checked_descriptors}, origin on {exact_origins}: all exact; "
        f"no exemptions)"
    )
    if lattice_origins:
        print(
            f"  note: {lattice_origins} row(s) differ by a parent lattice vector "
            "and failed above"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
