#!/usr/bin/env python3
"""Derive the little-group characters of the 73 other-wave-vector irreps.

The 5,756 `isotropy_w_subduce_*` rows name irreps of the eleven cubic parents
196, 202, 203, 209, 210, 216, 219, 225, 226, 227 and 228 that live at
**parametric** k domains (all of them lines through Gamma, one free parameter).
`data_irreps.txt` stores only their label/space group/dimension/type, so the task
9 audit reports them separately: without a character there is no frequency.

The official program knows those characters indirectly.  For a parent space
group it prints, for every Gamma irrep, the compatibility relation

    GM4: DT1 DT2 DT2          (constant `SHOW COMPATIBILITY`, k' = DT)

i.e. the decomposition of the Gamma irrep into the line's little irreps, and
`SHOW CHARACTER` gives that Gamma irrep's character on any space-group element.
At k = Gamma the line irreps carry no Bloch phase, so their characters are pure
point characters ``D_i(R)`` and the compatibility data gives one linear equation
per Gamma irrep and per point operation ``R`` of the little group::

    sum_i  n_GM * mult(GM -> i) * D_i(R)  =  chi_GM(R)

with ``n_GM`` the number of Miller-Love components of a compound Gamma row
(`GM2GM3` counts twice).  Solving that small exact rational system for every
``R`` recovers the little-group character table of the line irreps.  The
equations are over-determined in practice, which is the gate: every solved table
is re-substituted into all equations and must reproduce them exactly.

``--json`` writes the solved tables (with the equations they were derived from)
for the engine work that turns them into the 5,756 frequencies.

Usage::

    python3 scripts/freeze_w_little_characters.py --limit 1 --verbose
    python3 scripts/freeze_w_little_characters.py --json target/task9/w_little_characters.json

Exit codes: 0 = every source solved and verified, 1 = a source is inconsistent
or under-determined, 2 = the pinned archive or the extracted binary is missing.
"""

import argparse
import cmath
import json
import math
import os
import re
import subprocess
import sys
from fractions import Fraction

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
ISO_DIR = os.path.join(REPO_DIR, "isotropy_subgroup")

sys.path.insert(0, SCRIPT_DIR)

import generate_irrep_data as generator  # noqa: E402

ELEMENT_RE = re.compile(r"\(([^)]*)\)")
IRREP_ROW_RE = re.compile(r"^(\S+)\s+(\(.*?\))\s+(\d+)\s*$")
CHARACTER_ROW_RE = re.compile(r"\(([^)]*)\)\s+(-?[\d.]+)\s*$")
COMPONENT_RE = re.compile(r"[A-Z]{1,3}\d+[+-]?")


class Source:
    """One other-wave-vector irrep: parent, k domain, compact little label."""

    def __init__(self, sg, label, k_label, slot_dims, component_count):
        self.sg = sg
        self.label = label
        self.k_label = k_label
        self.slot_dims = slot_dims
        self.component_count = component_count

    def __repr__(self):
        return f"Source({self.sg}, {self.label}, k={self.k_label})"


def component_count(label):
    """Miller-Love components of a (possibly compound) little irrep label."""
    return len(COMPONENT_RE.findall(label))


def pinned_sources():
    """The 73 other-wave-vector irreps and the k slot each one lives in."""
    irrep_lines = generator.read_file("data_irreps.txt")
    little_lines = generator.read_file("data_little.txt")
    irrep_sec = generator.get_sections(irrep_lines)
    little_sec = generator.get_sections(little_lines)

    labels = generator.parse_labels(irrep_lines, irrep_sec, "irrep_w_label")
    space_groups = generator.parse_ints(irrep_lines, irrep_sec, "irrep_w_space_group")
    dimensions = generator.parse_ints(irrep_lines, irrep_sec, "irrep_w_dimension")

    little_label = generator.parse_labels(little_lines, little_sec, "little_irr_full_label")
    little_sg = generator.parse_ints(little_lines, little_sec, "little_irr_space_group")
    little_k = generator.parse_ints(little_lines, little_sec, "little_irr_k")
    little_count = generator.parse_ints(little_lines, little_sec, "little_irr_count")
    little_dims = generator.parse_ints(little_lines, little_sec, "little_irr_dim")
    k_labels = generator.parse_labels(little_lines, little_sec, "little_k_label")
    space_lines = generator.read_file("data_space.txt")
    space_sec = generator.get_sections(space_lines)
    lattice = generator.parse_ints(space_lines, space_sec, "ispace_lattice")

    compact = {}
    for index, label in enumerate(little_label):
        compact[(little_sg[index], little_k[index], label.strip())] = index

    sources = []
    for source, label in enumerate(labels):
        sg = space_groups[source]
        slot = None
        for key, index in compact.items():
            if key[0] == sg and key[2] == label.strip():
                slot = key[1]
                break
        if slot is None:
            raise RuntimeError(f"SG {sg} {label}: no compact little irrep entry")
        lat = lattice[sg - 1] - 1
        k_label = k_labels[lat * 27 + slot - 1].strip()
        offset = 0
        for previous in range(sg - 1):
            pass
        # dims of this (space group, k slot) live in the 12-wide buffer of the
        # slot; the sum of the non-zero entries is the little dimension sum.
        sources.append(
            Source(
                sg,
                label.strip(),
                k_label,
                None,
                component_count(label.strip()),
            )
        )
    return sources, dimensions, little_count, little_dims


K_ROW_RE = re.compile(r"^(\S+)\s+(\(.*?\))\s")


def parse_k_list(text):
    """`DISPLAY KPOINT` -> [(label, representative k vector), ...]."""
    rows = []
    for line in text.splitlines():
        match = K_ROW_RE.match(line.strip().rstrip("*").strip())
        if match:
            rows.append((match.group(1), match.group(2)))
    return rows


def points_on_line(k_rows, direction):
    """Special points `alpha * direction` on the same line (alpha rational)."""
    points = []
    for label, vector in k_rows:
        if re.search(r"[abg]", vector):
            continue  # a domain (line/plane), not a point
        point = parse_direction(vector)
        alpha = None
        on_line = True
        for axis in range(3):
            if direction[axis] == 0:
                if point[axis] != 0:
                    on_line = False
                    break
            else:
                candidate = point[axis] / direction[axis]
                if alpha is None:
                    alpha = candidate
                elif candidate != alpha:
                    on_line = False
                    break
        if on_line and alpha not in (None, 0):
            points.append((label, alpha))
    return points


def phase_exponent(alpha, direction, translation):
    """`alpha * (v . t)`; the Bloch phase of the line irrep at that point."""
    return alpha * sum(a * b for a, b in zip(direction, translation))


def solve_complex(matrix, rhs, tolerance=1e-9):
    """Gaussian elimination over the complex numbers (small systems)."""
    rows = [
        [complex(value) for value in row] + [complex(rhs[index])]
        for index, row in enumerate(matrix)
    ]
    width = len(rows[0]) - 1
    pivot_row = 0
    pivots = []
    for column in range(width):
        chosen = None
        best = tolerance
        for candidate in range(pivot_row, len(rows)):
            if abs(rows[candidate][column]) > best:
                chosen = candidate
                best = abs(rows[candidate][column])
        if chosen is None:
            continue
        rows[pivot_row], rows[chosen] = rows[chosen], rows[pivot_row]
        pivot = rows[pivot_row][column]
        rows[pivot_row] = [value / pivot for value in rows[pivot_row]]
        for other in range(len(rows)):
            if other != pivot_row and abs(rows[other][column]) > tolerance:
                factor = rows[other][column]
                rows[other] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(rows[other], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    solution = [0j] * width
    for row, column in enumerate(pivots):
        solution[column] = rows[row][-1]
    consistent = all(
        all(abs(row[column]) <= 1e-7 for column in range(width))
        and abs(row[-1]) <= 1e-7
        for row in rows[pivot_row:]
    )
    return solution, len(pivots), consistent


def run_iso(sg, commands, timeout=300):
    """Run the official program for one parent space group."""
    binary = os.path.join(ISO_DIR, "iso")
    if not os.path.isfile(binary):
        raise FileNotFoundError(
            "the ISOTROPY binary is not extracted; run\n"
            f"  unzip -o {os.path.join(ISO_DIR, 'iso.zip')} -d {ISO_DIR}"
        )
    full = [
        "PAGE 1000",
        "SC 1000",
        # `VALUE ELEMENT` accepts the ITA spelling, which is the one
        # `SHOW ELEMENTS` prints under this label setting.
        "LABEL ELEMENT INTERNATIONAL",
        f"VALUE PARENT {sg}",
    ] + commands + ["QUIT"]
    result = subprocess.run(
        [binary],
        input="\n".join(full) + "\n",
        capture_output=True,
        text=True,
        cwd=ISO_DIR,
        env=dict(os.environ, ISODATA=ISO_DIR + os.sep),
        timeout=timeout,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"iso exited with {result.returncode} for SG {sg}: {result.stderr.strip()[:200]}"
        )
    return result.stdout


def section(text, marker):
    """The text between a table header containing `marker` and the next prompt."""
    if marker not in text:
        raise RuntimeError(f"missing {marker!r} in iso output: {text[-300:]}")
    body = text.split(marker, 1)[1]
    return body.split("*", 1)[0]


ELEMENT_TERM_RE = re.compile(r"([+-]?)(\d+(?:/\d+)?)?([xyz]?)")


def parse_elements(text):
    """`(x,y,z), (x,-y+3/4,-z+3/4)` -> [(rotation, translation), ...].

    The program prints a space-group element in the International Tables form,
    i.e. as the image of a general point, so a screw or glide shows up as a
    constant inside a component (`-y+3/4`).  The `(R|t)` layout is accepted too.
    """
    body = section(text, "Elements")
    elements = []
    for chunk in ELEMENT_RE.findall(body):
        parts = [part.strip() for part in chunk.split(",")]
        translation_parts = None
        if len(parts) >= 6 and all(re.fullmatch(r"[-\d/]+", part) for part in parts[3:6]):
            translation_parts = [Fraction(part) for part in parts[3:6]]
        rotation = []
        translation = []
        for part in parts[:3]:
            row = [Fraction(0), Fraction(0), Fraction(0)]
            constant = Fraction(0)
            for sign, number, letter in ELEMENT_TERM_RE.findall(part.replace(" ", "")):
                if not number and not letter:
                    continue
                value = Fraction(number) if number else Fraction(1)
                if sign == "-":
                    value = -value
                if letter:
                    row["xyz".index(letter)] += value
                else:
                    constant += value
            rotation.append(tuple(row))
            translation.append(constant)
        if translation_parts is not None:
            translation = translation_parts
        elements.append((tuple(rotation), tuple(translation), chunk))
    return elements


def parse_k_vector(text, k_label):
    """The representative `(0,2a,0)` of the selected k domain."""
    for line in text.splitlines():
        stripped = line.strip()
        match = re.match(rf"^{re.escape(k_label)}\s+(\(.*?\))\s*$", stripped)
        if match:
            return match.group(1)
    raise RuntimeError(f"no k vector row for {k_label}: {text[-300:]}")


TERM_RE = re.compile(r"([+-]?)(\d+(?:/\d+)?)?([abg]?)")


def parse_direction(vector):
    """"(0,2a,0)" -> (Fraction(0), Fraction(2), Fraction(0)) at parameters = 1."""
    direction = []
    for part in vector.strip("()").split(","):
        total = Fraction(0)
        for sign, number, letter in TERM_RE.findall(part.replace(" ", "")):
            if not number and not letter:
                continue
            value = Fraction(number) if number else Fraction(1)
            total += -value if sign == "-" else value
        direction.append(total)
    return tuple(direction)


def parse_irrep_rows(text):
    """`DISPLAY IRREP` rows -> [(label, k vector, dimension), ...]."""
    rows = []
    for line in text.splitlines():
        match = IRREP_ROW_RE.match(line.strip())
        if match:
            rows.append((match.group(1), match.group(2), int(match.group(3))))
    return rows


def parse_compat(text):
    """`GM4: DT1 DT2 DT2` -> ["DT1", "DT2", "DT2"].

    The element/character column shares the row with the compatibility list, so
    everything from the first parenthesis on is dropped (labels carry none).
    """
    body = section(text, "Compat")
    match = re.search(r":\s*(.*)$", body.strip(), re.S)
    if not match:
        return []
    labels = match.group(1).split("(")[0]
    return labels.split()


def parse_characters(text):
    """`(x,y,z) 3.000` rows -> [(element string, Fraction), ...]."""
    characters = []
    for line in text.splitlines():
        # The element/character pair shares its row with the compatibility list
        # when both tables are displayed in one session.
        match = CHARACTER_ROW_RE.search(line.strip().rstrip("*").strip())
        if match:
            characters.append((match.group(1), Fraction(match.group(2))))
    return characters


def apply_rotation(rotation, vector):
    return tuple(
        sum(rotation[i][j] * vector[j] for j in range(3)) for i in range(3)
    )


def solve(matrix, rhs):
    """Exact rational solve; returns (solution, rank, consistent)."""
    rows = [
        [Fraction(value) for value in row] + [Fraction(rhs[i])]
        for i, row in enumerate(matrix)
    ]
    width = len(rows[0]) - 1
    pivot_row = 0
    pivots = []
    for column in range(width):
        chosen = None
        for candidate in range(pivot_row, len(rows)):
            if rows[candidate][column] != 0:
                chosen = candidate
                break
        if chosen is None:
            continue
        rows[pivot_row], rows[chosen] = rows[chosen], rows[pivot_row]
        pivot = rows[pivot_row][column]
        rows[pivot_row] = [value / pivot for value in rows[pivot_row]]
        for other in range(len(rows)):
            if other != pivot_row and rows[other][column] != 0:
                factor = rows[other][column]
                rows[other] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(rows[other], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    solution = [Fraction(0)] * width
    for row, column in enumerate(pivots):
        solution[column] = rows[row][-1]
    consistent = all(
        all(row[column] == 0 for column in range(width)) and row[-1] == 0
        for row in rows[pivot_row:]
    )
    return solution, len(pivots), consistent


def extra_point_equations(sg, k_label, direction, little, phase_sign=1):
    """Equations from the other special points of the same k line.

    At a point `k' = alpha * v` the line irrep's character on `(R, t)` is
    `exp(phase_sign * 2i pi alpha v.t) * D(R)`, so a point irrep's compatibility
    row gives `sum_i m_i D_i(R) = chi(R, t) * exp(-phase_sign * 2i pi alpha v.t)`.
    """
    k_rows = parse_k_list(
        run_iso(sg, ["SHOW KPOINT", "SHOW STAR", "DISPLAY KPOINT"])
    )
    equations = []
    points = []
    for label, alpha in points_on_line(k_rows, direction):
        irreps = parse_irrep_rows(
            run_iso(
                sg,
                [
                    f"VALUE KPOINT {label}",
                    "SHOW IRREP",
                    "SHOW DIMENSION",
                    "SHOW KPOINT",
                    "DISPLAY IRREP",
                ],
            )
        )
        used = 0
        for irrep_label, _, _ in irreps:
            commands = [f"VALUE IRREP {irrep_label}", "SHOW COMPATIBILITY"]
            commands += [f"VALUE COMPATIBILITY {k_label}", "DISPLAY IRREP"]
            commands += ["SHOW CHARACTER"]
            for element in little:
                commands += [
                    f"VALUE ELEMENT {element[2].upper().replace(',', ' ')}",
                    "DISPLAY IRREP",
                ]
            text = run_iso(sg, commands)
            compat = parse_compat(text)
            if not compat:
                continue
            characters = parse_characters(text)
            if len(characters) != len(little):
                continue
            multiplicity = {}
            for entry in compat:
                multiplicity[entry] = multiplicity.get(entry, 0) + 1
            equations.append(
                {
                    "point": label,
                    "alpha": str(alpha),
                    "irrep": irrep_label,
                    "multiplicity": multiplicity,
                    "characters": [complex(value) for _, value in characters],
                    "phases": [
                        cmath.exp(
                            -phase_sign
                            * 2j
                            * math.pi
                            * float(phase_exponent(alpha, direction, element[1]))
                        )
                        for element in little
                    ],
                }
            )
            used += 1
        if used:
            points.append((label, str(alpha), used))
    return equations, points


def resolve_from_extra_points(
    unknown_labels, equations, extra, wanted, expected_identity=None, tolerance=1e-6
):
    """Solve the combined system for the requested sources (complex arithmetic)."""
    matrix = []
    rhs = []
    for equation in equations:
        row = [
            Fraction(equation["multiplicity"].get(label, 0))
            for label in unknown_labels
        ]
        if all(value == 0 for value in row):
            continue
        matrix.append([complex(value) for value in row])
        rhs.append([complex(value) for value in equation["characters"]])
    for equation in extra:
        row = [
            Fraction(equation["multiplicity"].get(label, 0))
            for label in unknown_labels
        ]
        if all(value == 0 for value in row):
            continue
        matrix.append([complex(value) for value in row])
        rhs.append(
            [
                character / phase
                for character, phase in zip(equation["characters"], equation["phases"])
            ]
        )
    if not matrix:
        return {}, ["no equations"]
    # A^T, so the dual system `A^T y = v` fixes the requested linear functional
    # even when the individual components are not determined.
    transposed = [
        [matrix[row][column] for row in range(len(matrix))]
        for column in range(len(unknown_labels))
    ]
    results = {}
    for label in wanted:
        # A compound source covers several components.
        components = COMPONENT_RE.findall(label)
        vector = [
            Fraction(len([c for c in components if c == unknown]))
            for unknown in unknown_labels
        ]
        solution, _, consistent = solve_complex(transposed, vector)
        if not consistent:
            continue
        check = [
            sum(
                transposed[row][column] * solution[column]
                for column in range(len(solution))
            )
            for row in range(len(vector))
        ]
        if any(abs(a - complex(b)) > tolerance for a, b in zip(check, vector)):
            continue
        values = []
        ok = True
        for op_index in range(len(rhs[0])):
            total = sum(
                complex(dual) * values_op
                for dual, values_op in zip(solution, [row[op_index] for row in rhs])
            )
            rounded = round(total.real)
            if abs(total - rounded) > 1e-6:
                ok = False
                break
            values.append(rounded)
        if ok and expected_identity is not None:
            # The identity character is the dimension of the source.  The extra
            # points of a line have a star larger than one, so the program's
            # `SHOW CHARACTER` prints the *full* irrep's character there and not
            # the little character the compatibility row refers to; a solution
            # that fails this gate is such a contaminated read, not a result.
            expected = expected_identity.get(label)
            if expected is not None and values and values[0] != expected:
                continue
        if ok:
            results[label] = values
    if not results:
        return {}, ["the extra points do not determine any requested source"]
    return results, []


def is_identity(rotation):
    return all(
        rotation[row][column] == (1 if row == column else 0)
        for row in range(3)
        for column in range(3)
    )


def multiply(left, right):
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(3))
            for column in range(3)
        ]
        for row in range(3)
    ]


def cogroup_pair_route(little, labels, unknown_labels, matrix, rhs, duals):
    """The last two line irreps, which no Gamma row separates.

    SG 202/203/209/210 `DT3`/`DT4` are one-dimensional irreps of the little
    cogroup of the `DT` line -- `C2v = {E, C2, m1, m2}` for the two F-centered
    cubic groups with a mirror class, `C4 = {E, C2, R, R^3}` for the two without
    -- so every Gamma compatibility row, which sees the pair only through the
    combinations that `SHOW COMPATIBILITY` lists, fixes `DT3 + DT4` and leaves
    the difference free.  Two pinned facts close the gap:

    * the summed character `DT3 + DT4` *is* a determined functional of the
      Gamma system (the dual system `A^T y = e3 + e4` is consistent), which
      fixes the `C2` entry, and
    * the pinned data orders the four sources of such a domain the way the
      standard character table does: the two sources the Gamma data *does*
      determine come out as the first two rows `(+,+,+)` and `(+,-,-)` over
      `(C2, o2, o3)`, and they are the first two sources of the domain.

    The continuation of that ordering is `DT3 = (s, m, -m)` and
    `DT4 = (s, -m, m)` with `2 s` the summed `C2` character and `m = 1` for
    `C2v` (the two remaining sign characters) and `m = i` for `C4` (the
    conjugate pair, whose characters on the order-four generator are `+-i`).
    The identity character is the (one-dimensional) dimension of the source.

    These eight tables are the only frozen data that uses the source ordering
    instead of a Gamma row.  The ordering is checked per space group against the
    determined sources of the same domain, and the result is cross-checked
    against the archived CIR characters of the compatible points in
    `tests/w_little_characters.rs`.
    """
    problems = []
    if len(little) != 4:
        return {}, [f"the little cogroup has {len(little)} operations, not 4"]
    rotations = [
        [[int(value) for value in row] for row in element[0]] for element in little
    ]
    identity = [[1 if row == column else 0 for column in range(3)] for row in range(3)]
    if rotations[0] != identity:
        return {}, ["the little-group operations do not start with the identity"]
    involutions = [rotation for rotation in rotations if multiply(rotation, rotation) == identity]
    if len(involutions) == 4:
        kind = "C2v"
        if multiply(rotations[2], rotations[3]) != rotations[1]:
            return {}, ["the little cogroup is not C2v"]
    elif len(involutions) == 2:
        # C4 with the order-four generator at index 2 (its inverse at index 3).
        kind = "C4"
        if multiply(rotations[2], rotations[2]) != rotations[1] or multiply(
            rotations[2], rotations[3]
        ) != identity:
            return {}, ["the little cogroup is not C4"]
    else:
        return {}, ["the little cogroup is neither C2v nor C4"]

    solved = [label for label in labels if label in duals]
    missing = [label for label in labels if label not in duals]
    pair_labels = labels[-2:]
    # Only the last two sources of the domain may be open, and the pair's
    # summed character is what the Gamma rows determine.  A partially solved
    # pair (one sibling fixed by an extra point) is completed here.
    if len(labels) < 4 or not [label for label in labels[-2:] if label not in duals] or any(
        label not in duals for label in labels[:-2]
    ):
        return {}, [
            "the open sources are not the last two of the domain: "
            f"solved {solved}, missing {missing}"
        ]
    expected = [(1, 1, 1), (1, -1, -1)]
    for position, label in enumerate(labels[:2]):
        if label not in duals:
            continue
        # `rhs` holds the Gamma characters per operation; a source's own
        # character is its dual functional applied to them.
        own = [
            sum(duals[label][row] * rhs[row][op] for row in range(len(rhs)))
            for op in range(len(little))
        ]
        pattern = tuple(int(value) for value in own[1:])
        if pattern != expected[position]:
            return {}, [
                f"the determined source {label} is {pattern}, not the expected "
                f"cogroup pattern {expected[position]}: the source ordering "
                "convention does not hold for this domain"
            ]

    vector = [Fraction(0)] * len(unknown_labels)
    for label in pair_labels:
        vector[unknown_labels.index(label)] += 1
    transposed = [
        [matrix[row][column] for row in range(len(matrix))]
        for column in range(len(unknown_labels))
    ]
    dual, _, consistent = solve(transposed, vector)
    if not consistent:
        return {}, ["the summed character of the pair is not determined"]
    check = [
        sum(transposed[i][j] * dual[j] for j in range(len(dual)))
        for i in range(len(vector))
    ]
    if check != [Fraction(value) for value in vector]:
        return {}, ["the summed character of the pair is not determined"]
    total = [
        sum(dual[row] * rhs[row][op] for row in range(len(rhs)))
        for op in range(len(little))
    ]
    if total[0] != 2 or total[2] != 0 or total[3] != 0 or total[1] not in (2, -2):
        return {}, [f"the summed character of the pair is {total}, not (2,+-2,0,0)"]
    if kind == "C2v":
        pair = [(Fraction(1), Fraction(0)), (Fraction(-1), Fraction(0))]
    else:
        # The conjugate pair of C4: the order-four generator carries +-i.
        pair = [(Fraction(0), Fraction(1)), (Fraction(0), Fraction(-1))]
    sign = total[1] / 2
    values = {}
    for offset, label in enumerate(labels[-2:]):
        first, second = pair if offset == 0 else (pair[1], pair[0])
        if label in duals:
            # A sibling the extra-point route already fixed: the ordering
            # prediction must agree with it, otherwise the convention is wrong
            # for this domain and nothing is returned.
            own = [
                sum(duals[label][row] * rhs[row][op] for row in range(len(rhs)))
                for op in range(len(little))
            ]
            predicted = [
                Fraction(1),
                sign,
                first[0],
                second[0],
            ]
            if [Fraction(value) for value in own] != predicted or any(
                part[1] != 0 for part in (first, second)
            ):
                return {}, [
                    f"the ordered pair disagrees with the determined source {label}"
                ]
            continue
        values[label] = [
            (Fraction(1), Fraction(0)),
            (sign, Fraction(0)),
            first,
            second,
        ]
    return values, problems


def derive_line(sg, k_label, labels, verbose=False):
    """Solve the little-group character tables of every irrep on one k domain.

    The sources of a line share a little group, and a Gamma irrep's
    compatibility list mixes the line's irreps (`GM4: DT1 DT2 DT2`), so the
    system is solved jointly and split per source afterwards.
    """
    elements = parse_elements(run_iso(sg, ["SHOW ELEMENTS", "DISPLAY PARENT"]))
    k_text = run_iso(sg, [f"VALUE KPOINT {k_label}", "SHOW KPOINT", "DISPLAY KPOINT"])
    k_vector = parse_k_vector(k_text, k_label)
    direction = parse_direction(k_vector)
    little = [
        element
        for element in elements
        if apply_rotation(element[0], direction) == direction
    ]

    gamma = parse_irrep_rows(
        run_iso(
            sg,
            [
                "VALUE KPOINT GM",
                "SHOW IRREP",
                "SHOW DIMENSION",
                "SHOW KPOINT",
                "DISPLAY IRREP",
            ],
        )
    )
    lines = parse_irrep_rows(
        run_iso(
            sg,
            [
                f"VALUE KPOINT {k_label}",
                "SHOW IRREP",
                "SHOW DIMENSION",
                "SHOW KPOINT",
                "DISPLAY IRREP",
            ],
        )
    )
    line_labels = [label for label, _, _ in lines]
    printed = set(line_labels)
    for label in labels:
        # A compound source may be printed either as its own row (`DT3DT4`) or
        # as separate rows (`DT3`, `DT4`).
        if label in printed:
            continue
        missing = [c for c in COMPONENT_RE.findall(label) if c not in printed]
        if missing:
            return {}, [f"SG {sg}: {missing} not printed for k {k_label}"]

    equations = []
    for gamma_label, _, gamma_dim in gamma:
        commands = [f"VALUE IRREP {gamma_label}", "SHOW COMPATIBILITY"]
        commands += [f"VALUE COMPATIBILITY {k_label}", "DISPLAY IRREP"]
        commands += ["SHOW CHARACTER"]
        for element in little:
            commands += [
                f"VALUE ELEMENT {element[2].upper().replace(',', ' ')}",
                "DISPLAY IRREP",
            ]
        text = run_iso(sg, commands)
        multiplicity = {}
        for label in parse_compat(text):
            multiplicity[label] = multiplicity.get(label, 0) + 1
        count = component_count(gamma_label)
        characters = parse_characters(text)
        if len(characters) != len(little):
            return {}, [
                f"SG {sg} {gamma_label}: {len(characters)} characters for "
                f"{len(little)} little-group operations"
            ]
        equations.append(
            {
                "gamma": gamma_label,
                "dim": gamma_dim,
                "components": count,
                "compat": parse_compat(text),
                "characters": [str(value) for _, value in characters],
                "multiplicity": {
                    label: value * count for label, value in multiplicity.items()
                },
            }
        )

    # Every little irrep the compatibility data mentions is an unknown, so the
    # system is closed: a Gamma row mixing `DT1 DT2 DT2` cannot be truncated.
    unknown_labels = []
    for equation in equations:
        for label in equation["compat"]:
            if label not in unknown_labels:
                unknown_labels.append(label)

    matrix = []
    rhs = []
    used = []
    for equation in equations:
        row = [
            Fraction(equation["multiplicity"].get(label, 0))
            for label in unknown_labels
        ]
        if all(value == 0 for value in row):
            continue
        matrix.append(row)
        rhs.append([Fraction(value) for value in equation["characters"]])
        used.append(equation)

    # A source is only usable if its *own* character (the sum over its
    # components) is determined by the equations, which is weaker than every
    # single component being determined: a compound printed as one row has one
    # unknown, and conjugate components may only be fixed jointly.
    def functional(labels_of_source):
        vector = [Fraction(0)] * len(unknown_labels)
        for component in COMPONENT_RE.findall(labels_of_source):
            if component not in unknown_labels:
                return None
            vector[unknown_labels.index(component)] += 1
        return vector

    transposed = [
        [matrix[row][column] for row in range(len(matrix))]
        for column in range(len(unknown_labels))
    ]
    duals = {}
    problems = []
    for label in labels:
        vector = functional(label)
        if vector is None:
            problems.append(
                f"SG {sg} k {k_label}: {label} is not in the compatibility data"
            )
            continue
        dual, _, consistent = solve(transposed, vector)
        if not consistent:
            problems.append(
                f"SG {sg} k {k_label}: character of {label} is not determined by "
                "the compatibility data"
            )
            continue
        # `dual` solves A^T y = v, so v.x = y.(A x) = y.b for every solution x:
        # the source character is determined even when single components are not.
        check = [
            sum(transposed[i][j] * dual[j] for j in range(len(dual)))
            for i in range(len(vector))
        ]
        if check != [Fraction(v) for v in vector]:
            problems.append(
                f"SG {sg} k {k_label}: character of {label} is not determined by "
                "the compatibility data"
            )
            continue
        duals[label] = dual

    extra_solved = {}
    extra_used = []
    routes = {}
    missing = [label for label in labels if label not in duals]
    reported = set(problems)
    if missing:
        full_dim = {label: dim for label, _, dim in lines}
        star = len(elements) // max(len(little), 1)
        expected_identity = {}
        for label in missing:
            if label in full_dim and star and full_dim[label] % star == 0:
                expected_identity[label] = full_dim[label] // star
        extra, extra_used = extra_point_equations(sg, k_label, direction, little)
        if extra:
            extra_solved, _ = resolve_from_extra_points(
                unknown_labels, equations, extra, missing, expected_identity
            )
        remaining = [label for label in missing if label not in extra_solved]
        if remaining:
            # A pair of one-dimensional line irreps that only ever appears in
            # the Gamma rows through its sum is closed by the little cogroup's
            # standard source ordering; see `cogroup_pair_route`.
            pair_solved, pair_problems = cogroup_pair_route(
                little, labels, unknown_labels, matrix, rhs, duals
            )
            pair_solved = {
                label: values
                for label, values in pair_solved.items()
                if label in remaining
            }
            extra_solved.update(pair_solved)
            if pair_solved:
                routes.update({label: "cogroup_pair" for label in pair_solved})
                remaining = [label for label in missing if label not in extra_solved]
            for message in pair_problems:
                if remaining and message not in reported:
                    problems.append(
                        f"SG {sg} k {k_label}: character of {remaining[0]} is not "
                        f"determined by the compatibility data ({message})"
                    )
                    reported.add(message)
        for label in remaining:
            message = (
                f"SG {sg} k {k_label}: character of {label} is not determined "
                "by the compatibility data"
            )
            if message not in reported:
                problems.append(message)
                reported.add(message)
        routes.update({label: "extra_point" for label in extra_solved if label in missing})

    operations = []
    for index, element in enumerate(little):
        values = [values[index] for values in rhs]
        solution, rank, consistent = solve(matrix, values)
        if not consistent:
            return {}, [f"SG {sg} k {k_label}: incompatible characters at {element[2]}"]
        characters = {}
        for label in duals:
            characters[label] = (
                sum(a * b for a, b in zip(duals[label], values)),
                Fraction(0),
            )
        for label, solved in extra_solved.items():
            value = solved[index]
            # The Gamma route yields rational characters, the little-cogroup
            # route yields exact (real, imaginary) pairs.
            characters[label] = value if isinstance(value, tuple) else (Fraction(value), Fraction(0))
        operations.append(
            {
                "element": element[2],
                "rotation": [[str(value) for value in row] for row in element[0]],
                "translation": [str(value) for value in element[1]],
                "characters": {
                    label: [str(real), str(imaginary)]
                    for label, (real, imaginary) in characters.items()
                },
            }
        )
    if verbose:
        print(
            f"  SG {sg} k {k_label} (direction {k_vector}): {len(little)} "
            f"little-group operations, {len(equations)} Gamma irreps, "
            f"irreps {unknown_labels}",
            flush=True,
        )

    tables = {}
    expected_dimension = {}
    full_dims = {label: dim for label, _, dim in lines}
    star_size = len(elements) // max(len(little), 1)
    for label in labels:
        if label in full_dims and star_size and full_dims[label] % star_size == 0:
            expected_dimension[label] = full_dims[label] // star_size
    for label in labels:
        if label not in duals and label not in extra_solved:
            continue
        components = COMPONENT_RE.findall(label)
        tables[label] = {
            "space_group": sg,
            "label": label,
            "dimension": expected_dimension.get(label),
            "k_label": k_label,
            "k_vector": k_vector,
            "direction": [str(value) for value in direction],
            "components": components,
            "route": routes.get(label, "gamma"),
            "equations": equations,
            "extra_points": [
                {"point": point, "alpha": alpha, "irreps": used}
                for point, alpha, used in extra_used
            ]
            if routes.get(label) == "extra_point"
            else [],
            "operations": [
                {
                    "element": operation["element"],
                    "rotation": operation["rotation"],
                    "translation": operation["translation"],
                    "character": operation["characters"][label],
                }
                for operation in operations
            ],
        }
    # Drop the "not determined" messages of sources a later route solved.
    unresolved = [
        label
        for label in labels
        if label not in duals and label not in extra_solved
    ]
    stale = [
        f"SG {sg} k {k_label}: character of {label} is not determined by the "
        "compatibility data"
        for label in labels
        if label not in unresolved
    ]
    problems = [
        message
        for message in problems
        if not any(
            message == note or message.startswith(note + " (") for note in stale
        )
    ]
    return tables, problems


def character_parts(value):
    """Exact (real, imaginary) parts of a stored character, as integers."""
    real, imaginary = Fraction(value[0]), Fraction(value[1])
    for part in (real, imaginary):
        if part.denominator != 1:
            raise ValueError(f"character {value} is not a Gaussian integer")
    return str(real.numerator), str(imaginary.numerator)


def emit_rust(tables, missing, path):
    """Write the frozen character tables as a generated Rust data module."""
    lines = [
        "//! Little-group character tables of the parametric-k irreps behind the",
        "//! `other_wave_vector_subduction` rows.",
        "//!",
        "//! `DO NOT EDIT`: regenerate with",
        "//! `python3 scripts/freeze_w_little_characters.py --json <full.json> --rust src/irrep/w_little_characters_data.rs`.",
        "//!",
        "//! Every table is solved from the pinned ISOTROPY archive",
        "//! `isotropy_subgroup/iso.zip`",
        "//! (SHA-256 `568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b`),",
        "//! official `iso` 9.6.1: the Gamma irreps' compatibility rows for the line",
        "//! (`SHOW COMPATIBILITY`) fix the little irreps' characters through",
        "//! `sum_i n_Gamma mult(Gamma -> i) D_i(R) = chi_Gamma(R)` at k = Gamma, where",
        "//! no Bloch phase enters.  65 of the 73 sources are fixed that way.",
        "//! The remaining eight (SG 202/203/209/210 `DT3`/`DT4`) appear in every Gamma",
        "//! row only through `DT3 + DT4`, which the same rows do determine; the split",
        "//! uses the pinned source ordering of the little cogroup `C2v` (`DT1`, `DT2`",
        "//! come out as `A1`, `A2`, and the ordering check for `B1`, `B2` is verified",
        "//! against those two sources per space group).  `cogroup_pair_route` in",
        "//! `scripts/freeze_w_little_characters.py` documents the gate, and",
        "//! `tests/w_little_characters.rs` cross-checks the resulting tables against",
        "//! the archived CIR characters of the compatible points.  See",
        "//! `docs/isotropy-data-semantics.md` section 4.",
        "",
        "/// One little-group operation of a parametric-k irrep.",
        "#[derive(Debug, Clone, Copy)]",
        "pub struct LittleOperation {",
        "    /// International-Tables spelling, e.g. `-x,y,-z`.",
        "    pub element: &'static str,",
        "    /// Rotation part, row major, in the parent's conventional basis.",
        "    pub rotation: [[i8; 3]; 3],",
        "    /// Translation part in the parent's conventional cell, as fractions.",
        "    pub translation: [&'static str; 3],",
        "    /// Character of the little irrep on this operation, as exact",
        "    /// (real, imaginary) integer parts; the cubic `DT3`/`DT4` pair of",
        "    /// SG 209/210 is the conjugate pair with +-i on the order-four",
        "    /// generator, every other source is real.",
        "    pub character: [i32; 2],",
        "}",
        "",
        "/// Character table of one little irrep at a parametric k domain.",
        "#[derive(Debug, Clone, Copy)]",
        "pub struct LittleCharacterTable {",
        "    /// Parent space-group number.",
        "    pub space_group: u8,",
        "    /// Compact little-irrep label, e.g. `DT1`.",
        "    pub label: &'static str,",
        "    /// Program label of the k domain, e.g. `DT`.",
        "    pub k_label: &'static str,",
        "    /// Domain direction in the parent's conventional reciprocal basis.",
        "    pub direction: [&'static str; 3],",
        "    /// Dimension of the little irrep (the identity character).",
        "    pub dimension: u8,",
        "    /// Little-group operations and their characters.",
        "    pub operations: &'static [LittleOperation],",
        "}",
        "",
        "/// The 73 little irreps behind the parametric-k subduction rows.",
        "pub static W_LITTLE_CHARACTERS: &[LittleCharacterTable] = &[",
    ]
    for table in tables:
        ops = table["operations"]
        lines.append("    LittleCharacterTable {")
        lines.append(f"        space_group: {table['space_group']},")
        lines.append(f"        label: \"{table['label']}\",")
        lines.append(f"        k_label: \"{table['k_label']}\",")
        lines.append(
            "        direction: ["
            + ", ".join(f"\"{value}\"" for value in table["direction"])
            + "],"
        )
        lines.append(f"        dimension: {table['dimension']},")
        lines.append("        operations: &[")
        for op in ops:
            rotation = ", ".join(
                "[" + ", ".join(str(int(value)) for value in row) + "]"
                for row in op["rotation"]
            )
            translation = ", ".join(f"\"{value}\"" for value in op["translation"])
            lines.append(
                "            LittleOperation { element: \""
                + op["element"]
                + "\", rotation: ["
                + rotation
                + "], translation: ["
                + translation
                + "], character: ["
                + ", ".join(character_parts(op["character"]))
                + "] },"
            )
        lines.append("        ],")
        lines.append("    },")
    lines.append("];")
    lines.append("")
    lines.append("/// Sources the frozen data does not cover (empty when every source is")
    lines.append("/// covered by the Gamma rows or by the little-cogroup ordering route).")
    lines.append("pub static W_LITTLE_CHARACTERS_UNRESOLVED: &[(u8, &str)] = &[")
    for label, sg in sorted(missing):
        lines.append(f"    ({sg}, \"{label}\"),")
    lines.append("];")
    lines.append("")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    return len(tables)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", dest="json_path", default=None)
    parser.add_argument("--rust", dest="rust_path", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    try:
        sources, dimensions, _, _ = pinned_sources()
    except FileNotFoundError as error:
        print(f"pinned input missing: {error}", file=sys.stderr)
        return 2
    if not os.path.isfile(os.path.join(ISO_DIR, "iso")):
        print(
            "the ISOTOPY binary is not extracted; run\n"
            f"  unzip -o {os.path.join(ISO_DIR, 'iso.zip')} -d {ISO_DIR}",
            file=sys.stderr,
        )
        return 2
    if args.limit is not None:
        sources = sources[: args.limit]

    tables = []
    failures = []
    groups = {}
    for source in sources:
        groups.setdefault((source.sg, source.k_label), []).append(source.label)
    for (sg, k_label), labels in sorted(groups.items()):
        line_tables, problems = derive_line(sg, k_label, labels, verbose=args.verbose)
        failures.extend(problems)
        for label in labels:
            if label in line_tables:
                tables.append(line_tables[label])
    print(f"sources solved: {len(tables)} | failures: {len(failures)}")
    if args.rust_path:
        missing = []
        solved_keys = {(table["space_group"], table["label"]) for table in tables}
        for source in sources:
            if (source.sg, source.label) not in solved_keys:
                missing.append((source.label, source.sg))
        written = emit_rust(tables, missing, args.rust_path)
        print(f"wrote {written} tables to {args.rust_path}")
    for message in failures[:20]:
        print(f"FAIL {message}")
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "sources": len(sources),
                    "solved": len(tables),
                    "failures": failures,
                    "tables": tables,
                },
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
