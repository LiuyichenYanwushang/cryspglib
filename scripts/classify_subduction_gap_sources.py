#!/usr/bin/env python3
"""Classify where the target characters of the remaining subduction gaps come from.

The input is the gap manifest produced by
``examples/census_subduction_gaps.rs`` (one row per folded child star that the
pinned discrete child table cannot answer).  For every distinct
``(child space group, canonical star, setting)`` group this tool reports, with
exact rational arithmetic and without touching the production solver:

* the matching archived **PIR** records (little-group irreps of parametric k
  domains) and the parameter values that place the domain on the star;
* the records that are *less* constrained than the best match, i.e. the
  general-position records that also pass through the same point and must not be
  mistaken for the target there;
* the exact little co-group at the star point, computed from the space group
  operations in the same primitive frame as the archive;
* the projective factor system ``s_i s_j = T_{L_ij} s_k`` ->
  ``omega_ij = exp(+2 pi i q.L_ij)``, with the nontrivial values listed, so
  nonsymmorphic (screw/glide) little groups are visible instead of silently
  reusing the symmorphic character table;
* a classification of the group as an analytic general-position target (pure
  translation little group, the child-#1 route), a parameterized archived source
  that only needs its parameter substituted, or a real gap with no source.

Usage::

    python3 scripts/classify_subduction_gap_sources.py gaps.tsv > groups.tsv

The TSV goes to stdout, the summary to stderr.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction
import sys
from typing import List, Optional, Tuple

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))

from iso_irrep_exact import (  # noqa: E402
    ExactKArm,
    ExactSourceRecord,
    ExactSpaceGroupUniverse,
    load_exact_iso_irrep_sources,
)

Vector = Tuple[Fraction, Fraction, Fraction]
Rotation = Tuple[Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]]


HALF = Fraction(1, 2)
CENTERING_VECTORS = {
    "P": (),
    "A": ((0, HALF, HALF),),
    "B": ((HALF, 0, HALF),),
    "C": ((HALF, HALF, 0),),
    "I": ((HALF, HALF, HALF),),
    "F": ((0, HALF, HALF), (HALF, 0, HALF), (HALF, HALF, 0)),
    "R": ((Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)),
          (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3))),
}


def in_primitive_lattice(vector: Vector, centering: str) -> bool:
    """Whether ``vector`` (conventional direct coordinates) is a lattice vector.

    The archive lists one representative per centring coset, so an operation
    product can differ from the chosen representative by a centring vector; the
    primitive lattice is ``Z^3`` plus those centring vectors.
    """
    vectors = CENTERING_VECTORS[centering]
    ranges = [
        range(int(max(1 / value for value in (1,))) + 1) for _ in ()
    ]
    # Enumerate the centring coefficients over one period (their denominators
    # are 2 or 3, so a small box is exhaustive).
    from itertools import product as _product

    for coefficients in _product(range(3), repeat=len(vectors)):
        residual = list(vector)
        for coefficient, centring in zip(coefficients, vectors):
            for axis in range(3):
                residual[axis] -= coefficient * centring[axis]
        if all(value.denominator == 1 for value in residual):
            return True
    return False


def in_primitive_reciprocal(vector: Vector, centering: str) -> bool:
    """Whether ``vector`` is a reciprocal lattice vector of the primitive cell."""
    if any(value.denominator != 1 for value in vector):
        return False
    for centring in CENTERING_VECTORS[centering]:
        if sum(vector[axis] * centring[axis] for axis in range(3)).denominator != 1:
            return False
    return True


def parse_star(text: str) -> list[Vector]:
    """One star as the manifest writes it: points separated by ``;``."""
    points = []
    for point in text.split(";"):
        components = []
        for part in point.split(","):
            numerator, _, denominator = part.partition("/")
            components.append(Fraction(int(numerator), int(denominator) if denominator else 1))
        if len(components) != 3:
            raise ValueError(f"star point {point!r} is not three-dimensional")
        points.append((components[0], components[1], components[2]))
    return points


def free_directions(arm: ExactKArm) -> list[int]:
    """Axes whose parameter vector is a real direction of the k domain.

    The archive always carries three parameter slots; a discrete point stores
    zero vectors there, so "free" is decided by the vectors, not by the slots.
    """
    return [
        axis
        for axis in range(3)
        if arm.parameters[axis] is not None
        and any(value != 0 for value in arm.parameters[axis])
    ]


def solve_exact(columns: list[Vector], rhs: Vector) -> tuple[Fraction, ...] | None:
    """Solve ``sum_j t_j * columns[j] = rhs`` exactly, or return ``None``."""
    width = len(columns)
    rows = [[columns[j][i] for j in range(width)] + [rhs[i]] for i in range(3)]
    pivots: list[int] = []
    row = 0
    for column in range(width):
        pivot = next((i for i in range(row, 3) if rows[i][column] != 0), None)
        if pivot is None:
            continue
        rows[row], rows[pivot] = rows[pivot], rows[row]
        divisor = rows[row][column]
        rows[row] = [value / divisor for value in rows[row]]
        for i in range(3):
            if i != row and rows[i][column] != 0:
                factor = rows[i][column]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[row])]
        pivots.append(column)
        row += 1
        if row == 3:
            break
    for i in range(row, 3):
        if all(rows[i][j] == 0 for j in range(width)) and rows[i][width] != 0:
            return None
    solution = [Fraction(0)] * width
    for i, column in enumerate(pivots):
        solution[column] = rows[i][width]
    for i in range(3):
        total = sum(solution[j] * columns[j][i] for j in range(width))
        if total != rhs[i]:
            return None
    return tuple(solution)


def arm_point(arm: ExactKArm, free: list[int], parameters: tuple[Fraction, ...]) -> Vector:
    """The arm's k at the solved parameter values."""
    out = []
    for axis in range(3):
        value = arm.constant[axis]
        if axis in free:
            direction = arm.parameters[axis]
            assert direction is not None
            value += direction[axis] * parameters[free.index(axis)]
        out.append(value)
    return (out[0], out[1], out[2])


@lru_cache(maxsize=None)
def reciprocal_offsets(centering: str, span: int = 2):
    """Primitive reciprocal lattice vectors near the origin.

    An arm may represent the same physical k as the query with a different
    reciprocal representative (including the centring extinctions), so matching
    happens modulo this lattice and never by rounding.
    """
    for candidate in (
        (Fraction(i), Fraction(j), Fraction(k))
        for i in range(-span, span + 1)
        for j in range(-span, span + 1)
        for k in range(-span, span + 1)
    ):
        if in_primitive_reciprocal(candidate, centering):
            yield candidate


@lru_cache(maxsize=None)
def class_period(direction: Vector, centering: str) -> Fraction:
    """Smallest positive ``t`` with ``t * direction`` a primitive reciprocal vector.

    The parameter of a k domain is periodic modulo the reciprocal lattice, and
    for a centred lattice that period can exceed one: along ``(0, 1, 0)`` in a
    C-centred group only ``t = 2`` returns to the same class.  Sampling a
    shorter interval would miss q classes and misreport them as sourceless.
    """
    scale = 1
    for value in direction:
        scale = scale * value.denominator // __import__("math").gcd(scale, value.denominator)
    step = scale
    while True:
        candidate = tuple(value * step for value in direction)
        if in_primitive_reciprocal(candidate, centering):
            return Fraction(step)
        step += scale


@lru_cache(maxsize=None)
def arm_box(arm: ExactKArm, centering: str) -> Tuple[Vector, Vector]:
    """Bounding box of the arm's k over one full period of its parameters.

    Each parameter is sampled over ``[0, D]`` with ``D`` the common denominator
    of its direction, which covers every class of the domain modulo the
    reciprocal lattice.  Only offsets that land inside this box can match, so
    the search never enumerates the whole neighbourhood of the point.
    """
    lower = list(arm.constant)
    upper = list(arm.constant)
    for axis in range(3):
        direction = arm.parameters[axis]
        if direction is None or all(value == 0 for value in direction):
            continue
        period = class_period(direction, centering)
        for component in range(3):
            step = direction[component] * period
            if step > 0:
                upper[component] += step
            elif step < 0:
                lower[component] += step
    return (lower[0], lower[1], lower[2]), (upper[0], upper[1], upper[2])


def _ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


@lru_cache(maxsize=None)
def match_arm(
    arm: ExactKArm, point: Vector, centering: str
) -> Tuple[Optional[Tuple[Fraction, ...]], List[int]]:
    """Return the solved parameters placing ``arm`` on ``point``, if it does.

    ``point`` is a class modulo the primitive reciprocal lattice.  The candidate
    offsets are exactly the reciprocal translations that bring the point into
    the arm's bounding box; every candidate is then verified by an exact solve,
    so the match never rounds and never reports a neighbour.
    """
    free = free_directions(arm)
    lower, upper = arm_box(arm, centering)
    ranges = []
    for axis in range(3):
        start = _ceil_fraction(lower[axis] - point[axis])
        stop = (upper[axis] - point[axis]).numerator // (upper[axis] - point[axis]).denominator
        ranges.append(range(start, stop + 1))
    candidates = [
        (Fraction(i), Fraction(j), Fraction(k))
        for i in ranges[0]
        for j in ranges[1]
        for k in ranges[2]
    ]
    candidates.sort(key=lambda offset: sum(abs(value) for value in offset))
    columns = [arm.parameters[axis] for axis in free] if free else []
    for offset in candidates:
        if not in_primitive_reciprocal(offset, centering):
            continue
        target = tuple(point[i] + offset[i] for i in range(3))
        if not free:
            if tuple(arm.constant) == target:
                return (), free
            continue
        assert all(column is not None for column in columns)
        parameters = solve_exact(
            columns, tuple(target[i] - arm.constant[i] for i in range(3))
        )
        if parameters is None:
            continue
        if arm_point(arm, free, parameters) == target:
            return parameters, free
    return None, free


def matmul(left: Rotation, right: Rotation):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def rotation_inverse(rotation: Rotation) -> Rotation:
    """Inverse of an integer GL(3,Z) rotation, by adjugate over the determinant."""
    a, b, c = rotation
    cofactors = (
        (b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2], b[0] * c[1] - b[1] * c[0]),
        (c[1] * a[2] - c[2] * a[1], c[2] * a[0] - c[0] * a[2], c[0] * a[1] - c[1] * a[0]),
        (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]),
    )
    determinant = (
        a[0] * cofactors[0][0] + a[1] * cofactors[0][1] + a[2] * cofactors[0][2]
    )
    return tuple(
        tuple(Fraction(cofactors[i][j], determinant) for j in range(3)) for i in range(3)
    )


def preserves_q(rotation: Rotation, q: Vector, centering: str) -> bool:
    """Whether the rotation fixes ``q`` modulo the primitive reciprocal lattice.

    Direct fractional coordinates transform as ``x -> R x``, so reciprocal
    fractional coordinates transform as ``q -> (R^-1)^T q``; the difference must
    be a reciprocal vector of the *primitive* cell, which is exactly the set of
    integral vectors obeying the centring extinctions.
    """
    inverse = rotation_inverse(rotation)
    difference = tuple(
        sum(inverse[other][axis] * q[other] for other in range(3)) - q[axis]
        for axis in range(3)
    )
    return in_primitive_reciprocal(difference, centering)


def apply_rotation(rotation: Rotation, vector: Vector) -> Vector:
    return tuple(
        sum(Fraction(rotation[i][j]) * vector[j] for j in range(3)) for i in range(3)
    )


def little_group(
    universe: ExactSpaceGroupUniverse, q: Vector
) -> tuple[list[tuple[Rotation, Vector]], int]:
    """The little group of ``q`` and its co-group order.

    The returned operations are the archive's own operations (primitive frame)
    whose rotation fixes ``q``; the co-group order counts the distinct rotations
    among them.
    """
    centering = universe.centering.value
    operations = [
        (operation.rotation, operation.translation)
        for operation in universe.operations
        if preserves_q(operation.rotation, q, centering)
    ]
    return operations, len({rotation for rotation, _ in operations})


def factor_system(
    operations: list[tuple[Rotation, Vector]], q: Vector, centering: str
) -> tuple[list[tuple[int, int, Fraction]], int]:
    """The projective cocycle of the little group at ``q``.

    One representative per rotation is kept (the first occurrence).  For every
    ordered pair, ``s_i s_j = T_L s_k`` is solved in the group and the phase
    ``omega_ij`` is reported as a rational number of turns, ``q.L``; the second
    return value counts the operations whose translation is not a lattice
    vector (screw/glide in the primitive frame).
    """
    representatives: list[tuple[Rotation, Vector]] = []
    seen: set[Rotation] = set()
    for rotation, translation in operations:
        if rotation in seen:
            continue
        seen.add(rotation)
        representatives.append((rotation, translation))
    nontrivial: list[tuple[int, int, Fraction]] = []
    for first_index, (first_rotation, first_translation) in enumerate(representatives):
        for second_index, (second_rotation, second_translation) in enumerate(representatives):
            product_rotation = matmul(first_rotation, second_rotation)
            product_translation = tuple(
                apply_rotation(first_rotation, second_translation)[i] + first_translation[i]
                for i in range(3)
            )
            target = next(
                (
                    (rotation, translation)
                    for rotation, translation in representatives
                    if rotation == product_rotation
                ),
                None,
            )
            if target is None:
                raise AssertionError("little group is not closed under multiplication")
            lattice = tuple(product_translation[i] - target[1][i] for i in range(3))
            if not in_primitive_lattice(lattice, centering):
                raise AssertionError("operation product left the primitive lattice")
            turns = sum(q[i] * lattice[i] for i in range(3))
            turns -= turns.numerator // turns.denominator  # into [0, 1)
            if turns != 0:
                nontrivial.append((first_index, second_index, turns))
    nonsymmorphic = sum(
        1 for _, translation in representatives if any(value.denominator != 1 for value in translation)
    )
    return nontrivial, nonsymmorphic


def classify(
    matches: List[Tuple[ExactSourceRecord, List[int], Tuple[Fraction, ...]]],
    co_group_order: int,
) -> Tuple[str, int]:
    """The source classification and the number of excluded generic records.

    ``co_group_order == 1`` means the little group at the star point is the
    translation group, so the target is a one-dimensional Bloch phase and no
    archived character data is needed at all (the child-#1 route).  That holds
    whether or not some archived domain passes through the point.

    With a larger little group the archived parametric records are the natural
    source.  A general-position record that passes through the point is *not* a
    source there -- its little group is smaller than the actual one -- so a
    point where only general-position records match is reported as a special
    value with no archived domain, never silently answered from the general
    record.
    """
    if not matches:
        return "no_source_needs_algorithm", 0
    best = min(len(free) for _, free, _ in matches)
    kept = [entry for entry in matches if len(entry[1]) == best]
    excluded = len(matches) - len(kept)
    if co_group_order == 1:
        return "analytic_general_position", excluded
    if best == 3:
        return "special_value_no_source", excluded
    return "parameterized_source", excluded


def analyse_group(
    universe: ExactSpaceGroupUniverse,
    records: list[ExactSourceRecord],
    points: list[Vector],
) -> dict:
    """Every reported quantity of one ``(child, star, setting)`` group."""
    operations, co_group_order = little_group(universe, points[0])
    centering = universe.centering.value
    matches: list[tuple[ExactSourceRecord, list[int], tuple[Fraction, ...]]] = []
    for record in records:
        best_for_record: Optional[tuple[list[int], tuple[Fraction, ...]]] = None
        for arm in record.k_arms:
            for point in points:
                parameters, free = match_arm(arm, point, centering)
                if parameters is None:
                    continue
                if best_for_record is None or len(free) < len(best_for_record[0]):
                    best_for_record = (free, parameters)
                break
        if best_for_record is not None:
            matches.append((record, best_for_record[0], best_for_record[1]))
    classification, excluded = classify(matches, co_group_order)
    best = min((len(free) for _, free, _ in matches), default=-1)
    kept = [entry for entry in matches if len(entry[1]) == best]
    nontrivial, nonsymmorphic = factor_system(operations, points[0], centering)
    return {
        "matched_irnumbers": sorted({record.irnumber for record, _, _ in kept}),
        "matched_labels": sorted({record.irrep_label for record, _, _ in kept}),
        "matched_dimensions": sorted({record.dimension for record, _, _ in kept}),
        "matched_parameters": sorted({tuple(str(value) for value in parameters) for _, _, parameters in kept}),
        "excluded_generic": excluded,
        "min_free_directions": best,
        "little_group_ops": len(operations),
        "little_co_group_order": co_group_order,
        "factor_system_nontrivial": len(nontrivial),
        "factor_system_values": sorted({str(value) for _, _, value in nontrivial}),
        "nonsymmorphic_ops": nonsymmorphic,
        "classification": classification,
    }


HEADER = (
    "child_sg\tcanonical_q\tsetting_numerator\tsetting_denominator\tchild_shift\twitness_ordinal"
    "\twitness_probe\tmatched_irnumbers\tmatched_labels\tmatched_dimensions\tmatched_parameters"
    "\texcluded_generic\tmin_free_directions\tlittle_group_ops\tlittle_co_group_order"
    "\tfactor_system_nontrivial\tfactor_system_values\tnonsymmorphic_ops\tclassification"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", help="gap manifest TSV from census_subduction_gaps")
    parser.add_argument("--limit", type=int, default=0, help="only the first N groups")
    args = parser.parse_args()

    import csv

    rows = [row for row in csv.DictReader(open(args.manifest), delimiter="\t") if row["status"] == "missing_discrete_scalar_data"]
    groups: dict[tuple, dict] = {}
    for row in rows:
        key = (
            int(row["child_sg"]),
            row["canonical_q"],
            row["setting_numerator"],
            row["setting_denominator"],
            row["child_shift"],
        )
        groups.setdefault(key, row)

    database = load_exact_iso_irrep_sources()
    records_by_sg: dict[int, list[ExactSourceRecord]] = defaultdict(list)
    for record in database.pir_records:
        records_by_sg[record.spacegroup].append(record)

    keys = sorted(groups)
    if args.limit:
        keys = keys[: args.limit]

    print(HEADER)
    summary: Counter = Counter()
    per_child: dict[int, Counter] = defaultdict(Counter)
    for key in keys:
        child_sg, q_text, setting_numerator, setting_denominator, child_shift = key
        row = groups[key]
        points = parse_star(q_text)
        universe = database.source_universe(child_sg)
        report = analyse_group(universe, records_by_sg[child_sg], points)
        print(
            "\t".join(
                [
                    str(child_sg),
                    q_text,
                    setting_numerator,
                    setting_denominator,
                    child_shift,
                    row["ordinal"],
                    row["probe_cdml"],
                    ",".join(str(number) for number in report["matched_irnumbers"]),
                    ",".join(report["matched_labels"]),
                    ",".join(str(value) for value in report["matched_dimensions"]),
                    ";".join(",".join(values) for values in report["matched_parameters"]),
                    str(report["excluded_generic"]),
                    str(report["min_free_directions"]),
                    str(report["little_group_ops"]),
                    str(report["little_co_group_order"]),
                    str(report["factor_system_nontrivial"]),
                    ",".join(report["factor_system_values"]),
                    str(report["nonsymmorphic_ops"]),
                    report["classification"],
                ]
            )
        )
        summary[report["classification"]] += 1
        summary["factor_system_nontrivial"] += 1 if report["factor_system_nontrivial"] else 0
        summary["nonsymmorphic_children"] += 1 if report["nonsymmorphic_ops"] else 0
        summary["excluded_generic_groups"] += 1 if report["excluded_generic"] else 0
        per_child[child_sg][report["classification"]] += 1

    total = len(keys)
    print(f"groups={total}", file=sys.stderr)
    for name, count in sorted(summary.items()):
        print(f"  {name}={count}", file=sys.stderr)
    print(f"  children={len(per_child)}", file=sys.stderr)
    for child_sg in sorted(per_child, key=lambda sg: -sum(per_child[sg].values()))[:10]:
        print(f"  child #{child_sg}: {dict(per_child[child_sg])}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
