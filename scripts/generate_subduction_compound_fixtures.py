#!/usr/bin/env python3
"""Raw CIR full-matrix witnesses for compound and conjugate-star restrictions.

Each entry is a complex source constituent, never a character assembled by
the production compound adapter. The shared extractor verifies the archive
checksum and reads the actual full-matrix and individual arm-block traces.
"""

import argparse
from collections import Counter, deque
import cmath
from fractions import Fraction
import hashlib
import math
from pathlib import Path
import re
import zipfile

try:
    from .generate_subduction_star_fixtures import render
    from . import generate_irrep_data as source
    from .iso_irrep_exact import CIR_ARCHIVE_SHA256
    from .verify_isotropy_oracle import PRIMITIVE_BASIS
except ImportError:
    from generate_subduction_star_fixtures import render
    import generate_irrep_data as source
    from iso_irrep_exact import CIR_ARCHIVE_SHA256
    from verify_isotropy_oracle import PRIMITIVE_BASIS


CASES = {
    (19, "R1"): 559,
    (23, "W1"): 716,
    (45, "W1"): 1805,
    (45, "W2"): 1806,
    (46, "W1"): 1839,
    (46, "W2"): 1840,
    (76, "Z1"): 3852,
    (76, "Z3"): 3854,
    (83, "GM3+"): 4077,
    (83, "GM4+"): 4078,
    (167, "T1"): 8293,
    (167, "T2"): 8294,
    (182, "H1"): 8988,
    (182, "H2"): 8989,
    (198, "R2"): 9843,
    (219, "L3"): 10640,
    (220, "P1"): 10677,
}
ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "tests/data/subduction_compound_cir.rs"


def check_self_conjugate_characters(header, raw_arms, index, lines):
    """Compare D and D* on every source coset and finite translation phase.

    Each phase vector lists the eigenvalues of one translation on all arms.
    The three primitive translations generate their finite image. Exhausting
    that image checks arbitrary lattice shifts, including centring, without
    treating the full-star trace as a single Bloch eigenvalue.
    """
    arms = []
    for offset in range(0, len(raw_arms), 16):
        raw = raw_arms[offset:offset + 16]
        if any(raw[i] for i in (4, 5, 6, 8, 9, 10, 12, 13, 14)):
            raise ValueError(f"parameterized type-2 source {header['irnumber']}")
        arms.append(tuple(Fraction(raw[i], raw[3]) for i in range(3)))
    basis = PRIMITIVE_BASIS[header["space_group_symbol"][0]]
    generators = [tuple(sum(a * b for a, b in zip(k, row)) % 1 for k in arms) for row in basis]
    zero = (Fraction(0),) * len(arms)
    phases, queue = {zero}, deque([zero])
    while queue:
        phase = queue.popleft()
        for generator in generators:
            new = tuple((a + b) % 1 for a, b in zip(phase, generator))
            if new not in phases:
                phases.add(new)
                queue.append(new)
    eigenvalues = [[cmath.exp(2j * math.pi * float(value)) for value in phase] for phase in phases]
    dim = header["dim"]
    little_dim = dim // len(arms)
    comparisons = 0
    for operation in range(header["opcount"]):
        context = f"type-2 source {header['irnumber']} operation {operation}"
        _, index = source._parse_cir_operation_row(lines, index, context)
        matrix, _, index = source._parse_cir_complex_block(lines, index, dim * dim, context)
        blocks = []
        for arm in range(len(arms)):
            trace = 0j
            for i in range(arm * little_dim, (arm + 1) * little_dim):
                re_part, im_part = matrix[i * dim + i]
                trace += complex(re_part.materialize(), im_part.materialize())
            blocks.append(trace)
        for phase in eigenvalues:
            trace = sum(value * eigenvalue for value, eigenvalue in zip(blocks, phase))
            if abs(trace.imag) > 1e-9:
                raise ValueError(f"non-self-conjugate type-2 character: {context}, {trace}")
            comparisons += 1
    return comparisons


def audit_compound_sources():
    """Check the same-k constituent assumption against all frozen CIR links.

    Metadata supplies exact source IDs; no ML-label splitting is performed.
    Compare the complete affine k expressions, including any parameters. This
    checks source binding/arm geometry, not every character of all 672 rows.
    """
    archive = ROOT / "isotropy_subgroup/CIR_data.zip"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != CIR_ARCHIVE_SHA256:
        raise ValueError("CIR archive checksum differs from the pinned source")
    with zipfile.ZipFile(archive) as zipped:
        lines = zipped.read("CIR_data.txt").decode("ascii").splitlines()
    generated = (ROOT / "src/irrep/generated_data.rs").read_text()
    blocks = re.findall(r"CompoundMetadata \{([^}]+)\}", generated)
    if len(blocks) != 672:
        raise ValueError(f"compound metadata census changed: {len(blocks)}")
    entries = []
    for block in blocks:
        sg = int(re.search(r"sg: (\d+)", block)[1])
        ids = tuple(map(int, re.search(r"cir_irnumbers: \[(\d+), (\d+)\]", block).groups()))
        labels = re.search(r'cir_labels: \["([^"]+)", "([^"]+)"\]', block).groups()
        dims = tuple(map(int, re.search(r"cir_dimensions: \[(\d+), (\d+)\]", block).groups()))
        kind = re.search(r"semantics: CompoundCharacterSemantics::(\w+)", block)[1]
        entries.append((sg, ids, labels, dims, kind))
    needed = {irnumber for _, ids, _, _, _ in entries for irnumber in ids}
    headers, matrix_starts = {}, {}
    for index, line in enumerate(lines):
        if '"' not in line:
            continue
        header = source._parse_cir_header(line, index + 1)
        if header["irnumber"] not in needed:
            continue
        arms, next_index = source._read_exact_cir_block(
            lines, index + 1, 16 * header["kcount"], "compound source arms",
            source._parse_cir_integer,
        )
        headers[header["irnumber"]] = (header, arms)
        matrix_starts[header["irnumber"]] = next_index
    census = Counter()
    equivalent_characters = 0
    for sg, ids, labels, dims, kind in entries:
        constituents = [headers[irnumber] for irnumber in ids]
        for (header, _), label, dim in zip(constituents, labels, dims):
            if (header["sg"], header["label"]) != (sg, label):
                raise ValueError(f"source identity mismatch: {sg} {ids}")
            if header["dim"] != dim * header["kcount"]:
                raise ValueError(f"source dimension mismatch: {sg} {ids}")
        if constituents[0][1] != constituents[1][1] or dims[0] != dims[1]:
            raise ValueError(f"different constituent arm expressions/dimensions: {sg} {ids}")
        first, second = (item[0]["irtype"] for item in constituents)
        if kind == "DistinctComponentSum":
            if ids[0] == ids[1] or (first, second) != (3, 3):
                raise ValueError(f"distinct source semantics mismatch: {sg} {ids}")
            census["distinct"] += 1
        elif kind == "ConjugateRealification":
            if ids[0] != ids[1] or first not in (2, 3):
                raise ValueError(f"realification source semantics mismatch: {sg} {ids}")
            census[f"realification_irtype_{first}"] += 1
            if first == 2:
                equivalent_characters += check_self_conjugate_characters(
                    constituents[0][0], constituents[0][1], matrix_starts[ids[0]], lines,
                )
        else:
            raise ValueError(f"unknown compound semantics: {kind}")
    expected = {"distinct": 519, "realification_irtype_2": 41, "realification_irtype_3": 112}
    if census != expected:
        raise ValueError(f"compound source census changed: {dict(census)}")
    if equivalent_characters != 1484:
        raise ValueError(f"type-2 character coverage changed: {equivalent_characters}")
    return dict(census), equivalent_characters


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    census, equivalent_characters = audit_compound_sources()
    expected = render(CASES, "scripts/generate_subduction_compound_fixtures.py --write")
    if args.write:
        OUTPUT.write_text(expected)
    elif OUTPUT.read_text() != expected:
        raise SystemExit("compound source fixtures differ; regenerate with --write")
    print(f"compound CIR fixtures: {len(CASES)} source records, archive checksum verified")
    print(f"all compound source bindings and arm expressions checked: {census}")
    print(f"type-2 full characters agree with their conjugates: {equivalent_characters} coset/translation comparisons")
