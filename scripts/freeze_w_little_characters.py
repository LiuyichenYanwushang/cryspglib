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
            extra_solved, extra_problems = resolve_from_extra_points(
                unknown_labels, equations, extra, missing, expected_identity
            )
            for label in missing:
                if label in extra_solved:
                    continue
                message = (
                    f"SG {sg} k {k_label}: character of {label} is not determined "
                    "by the compatibility data"
                )
                if message not in reported:
                    problems.append(message)
                    reported.add(message)
        else:
            for label in missing:
                message = (
                    f"SG {sg} k {k_label}: character of {label} is not determined "
                    "by the compatibility data"
                )
                if message not in reported:
                    problems.append(message)
                    reported.add(message)

    operations = []
    for index, element in enumerate(little):
        values = [values[index] for values in rhs]
        solution, rank, consistent = solve(matrix, values)
        if not consistent:
            return {}, [f"SG {sg} k {k_label}: incompatible characters at {element[2]}"]
        characters = {}
        for label in duals:
            characters[label] = sum(a * b for a, b in zip(duals[label], values))
        for label, solved in extra_solved.items():
            characters[label] = Fraction(solved[index])
        operations.append(
            {
                "element": element[2],
                "rotation": [[str(value) for value in row] for row in element[0]],
                "translation": [str(value) for value in element[1]],
                "characters": {
                    label: str(value) for label, value in characters.items()
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
    for label in labels:
        if label not in duals and label not in extra_solved:
            continue
        components = COMPONENT_RE.findall(label)
        tables[label] = {
            "space_group": sg,
            "label": label,
            "k_label": k_label,
            "k_vector": k_vector,
            "direction": [str(value) for value in direction],
            "components": components,
            "equations": equations,
            "extra_points": [
                {"point": point, "alpha": alpha, "irreps": used}
                for point, alpha, used in extra_used
            ]
            if label in extra_solved
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
    return tables, problems


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", dest="json_path", default=None)
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
