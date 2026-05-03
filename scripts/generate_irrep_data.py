#!/usr/bin/env python3
"""
Parse iso_data files and generate Rust source code for cryspglib.

INPUT:  iso_data/data_irreps.txt, data_isotropy.txt, data_little.txt,
        data_images.txt, data_space.txt
OUTPUT: src/irrep/generated_data.rs (flat arrays for programmatic use)
        Updates src/irrep/{triclinic,monoclinic,...,cubic}.rs with rustdoc tables

All arrays in data_irreps.txt and data_isotropy.txt are PARALLEL:
position N in one array corresponds to position N in all others.
"""

import re, sys, os, zipfile, io
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(__file__)
ISO_DIR = os.path.join(SCRIPT_DIR, "..", "isotropy_subgroup")
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "src", "irrep")

# Import the direction mapping module from this directory
sys.path.insert(0, os.path.dirname(__file__))
from direction_map import build_direction_map

# ── zip-based file reading ──────────────────────────────────────────────────

def _open_zip_path(zip_name, inner_path):
    """Open a file from inside a zip archive, or from extracted dir if exists."""
    extracted = os.path.join(ISO_DIR, zip_name.replace('.zip', ''), inner_path)
    if os.path.exists(extracted):
        return open(extracted)
    zip_path = os.path.join(ISO_DIR, zip_name)
    if os.path.exists(zip_path):
        zf = zipfile.ZipFile(zip_path)
        # Find matching file inside zip
        for name in zf.namelist():
            if name.endswith(inner_path) or name == inner_path:
                return io.TextIOWrapper(zf.open(name))
        raise FileNotFoundError(f"{inner_path} not found in {zip_name}")
    raise FileNotFoundError(f"Neither {extracted} nor {zip_path} found")

def read_file(inner_path, zip_name="iso.zip"):
    """Read lines from a file inside a zip archive."""
    with _open_zip_path(zip_name, inner_path) as f:
        return f.readlines()

def get_sections(lines):
    """Return {section_name: line_index} for all section headers."""
    sections = {}
    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        # section headers: only lowercase letters and underscores
        if s[0].isalpha() and s[0].islower():
            # exclude quoted strings
            if not s.startswith('"') and not s.startswith("'"):
                sections[s] = i
    return sections

def parse_labels(lines, sections, name):
    """Extract '...' or \"...\" quoted labels from a section."""
    start = sections[name] + 1
    keys = list(sections.keys())
    idx = keys.index(name)
    end = sections[keys[idx + 1]] if idx + 1 < len(keys) else len(lines)
    data = []
    for line in lines[start:end]:
        for m in re.finditer(r'"([^"]*)"', line):
            data.append(m.group(1).strip())
        for m in re.finditer(r"'([^']*)'", line):
            val = m.group(1).strip()
            if val and val != '***':
                data.append(val)
    return data

def parse_ints(lines, sections, name):
    """Extract space-separated integers from a section."""
    start = sections[name] + 1
    keys = list(sections.keys())
    idx = keys.index(name)
    end = sections[keys[idx + 1]] if idx + 1 < len(keys) else len(lines)
    data = []
    for line in lines[start:end]:
        for token in line.split():
            try:
                data.append(int(token))
            except ValueError:
                pass
    return data

def parse_floats(lines, sections, name):
    """Extract float values from a section."""
    start = sections[name] + 1
    keys = list(sections.keys())
    idx = keys.index(name)
    end = sections[keys[idx + 1]] if idx + 1 < len(keys) else len(lines)
    data = []
    for line in lines[start:end]:
        for token in line.split():
            try:
                data.append(float(token))
            except ValueError:
                pass
    return data

# ── LaTeX label conversion ───────────────────────────────────────────────────

# Mapping from k-point prefix to LaTeX Greek letter
KPOINT_LATEX = {
    "GM": r"\Gamma",   # Gamma point
    "G":  r"\Gamma",
    "X":  "X",
    "M":  "M",
    "R":  "R",
    "A":  "A",
    "H":  "H",
    "K":  "K",
    "L":  "L",
    "Y":  "Y",
    "Z":  "Z",
    "T":  "T",
    "S":  "S",
    "U":  "U",
    "V":  "V",
    "W":  "W",
    "DT": r"\Delta",   # Delta line (Γ→X)
    "LD": r"\Lambda",  # Lambda line (Γ→R)
    "SM": r"\Sigma",   # Sigma line (Γ→M)
    "F":  "F",
    "B":  "B",
    "C":  "C",
    "D":  "D",
    "E":  "E",
    "N":  "N",
    "P":  "P",
    "Q":  "Q",
    "GP": "GP",        # general point
}

def label_to_latex(label):
    """Convert a CDML irrep label to LaTeX math notation.

    Examples:
        GM4+  → \\Gamma_4^+
        X3-   → X_3^-
        DT1   → \\Delta_1
        LD2   → \\Lambda_2
        SM3   → \\Sigma_3
        R1+   → R_1^+
        T4    → T_4
        GP1   → GP_1
    """
    # Find the k-point prefix (letters) and the rest (digits + signs)
    m = re.match(r'^([A-Za-z]+)(\d*)([+-]?)(.*)$', label)
    if not m:
        return label

    prefix = m.group(1)
    number = m.group(2)
    sign   = m.group(3)
    rest   = m.group(4)

    latex_prefix = KPOINT_LATEX.get(prefix, prefix)

    if number:
        result = f"{latex_prefix}_{{{number}}}"
    else:
        result = latex_prefix

    if sign == '+':
        result += "^+"
    elif sign == '-':
        result += "^-"

    if rest:
        result += rest  # e.g. combined labels like "H2H3"

    return result

# Dimension: image label first letter encodes dimension
# A=1, B=2, C=3, D=4, E=6, F=8, G=12, H=16, J=24
IMAGE_DIM = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 6, 'F': 8, 'G': 12, 'H': 16, 'J': 24}

# ── PIR k-vector parsing ──────────────────────────────────────────────────────

def _read_pir_lines():
    """Read PIR_data.txt from the zip archive."""
    zip_path = os.path.join(ISO_DIR, "PIR_data.zip")
    zf = zipfile.ZipFile(zip_path)
    with zf.open("PIR_data.txt") as f:
        return io.TextIOWrapper(f).readlines()


def _parse_pir_kvectors():
    """Parse PIR_data.txt and return dict: (SG#, ML_label) -> (kx, ky, kz, denom)."""
    lines = _read_pir_lines()

    kvec_map = {}
    i = 3  # skip 3 header lines
    while i < len(lines):
        line = lines[i].strip()
        if not line or not line.split()[0].isdigit():
            i += 1
            continue

        tokens = line.split()
        sg_num = int(tokens[1])
        quoted = re.findall(r'"([^"]*)"', line)
        if len(quoted) < 2:
            i += 1
            continue
        label = quoted[1].strip()

        # Get pmkcount from after the second quote
        after = line[line.rindex('"') + 1:].strip().split()
        if len(after) < 5:
            i += 1
            continue
        pmkcount = int(after[3])

        # Next line: k-vectors
        i += 1
        if i >= len(lines):
            break
        kvals = [int(x) for x in lines[i].strip().split()]

        if len(kvals) >= 4:
            key = (sg_num, label)
            if key not in kvec_map:
                kvec_map[key] = (kvals[0], kvals[1], kvals[2], kvals[3])

        # Skip to next irrep (find next line with double quotes and digits)
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if nxt and nxt.split()[0].isdigit() and '"' in nxt:
                break
            i += 1

    return kvec_map


# ── PIR character parsing ─────────────────────────────────────────────────────

def _round_char(x, eps=1e-8):
    """Round a character value to clean floating point, removing numerical noise."""
    if abs(x) < eps:
        return 0.0
    r = round(x)
    if abs(x - r) < eps:
        return float(r)
    # Try common rational fractions
    for n in range(-12, 13):
        for d in (2, 3, 4, 6, 8):
            if abs(x - n / d) < eps:
                return n / d
    return round(x, 10)


def is_int_4(s):
    """Check if s is a small integer fitting Fortran I3 format (irtranslation)."""
    try:
        v = int(s)
        return -999 <= v <= 999
    except ValueError:
        return False


def _parse_pir_characters():
    """Parse PIR_data.txt and return two dicts:
	    chars_map:    (SG#, ML_label) -> [char1, char2, ..., charN]
	    matrices_map: (SG#, ML_label) -> [m11, m12, ..., mNN] flat values

    Parses the irrep matrix elements from PIR_data.txt, computes the character
    (trace) for each symmetry operation, and rounds to clean float values.
    """
    lines = _read_pir_lines()

    chars_map = {}
    matrices_map = {}
    i = 3  # skip 3 header lines

    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        tokens = line.split()
        if not tokens or not tokens[0].isdigit():
            i += 1
            continue

        # Parse header: ir# SG# "sym" "label" dim type kcount pmkcount opcount
        quoted = re.findall(r'"([^"]*)"', line)
        if len(quoted) < 2:
            i += 1
            continue

        sg = int(tokens[1])
        label = quoted[1].strip()
        after_line = line[line.rindex('"') + 1:].strip().split()
        if len(after_line) < 5:
            i += 1
            continue
        dim = int(after_line[0])
        pmkcount = int(after_line[3])
        opcount = int(after_line[4])

        # Skip k-vector data: pmkcount * 16 integers across lines.
        # When pmkcount > 1, the k-vector occupies multiple lines
        # and the old code only skipped one line, causing all
        # subsequent character data to be shifted.
        kvec_remaining = pmkcount * 16
        while kvec_remaining > 0 and i < len(lines) - 1:
            i += 1
            if i >= len(lines):
                break
            kvec_remaining -= len(lines[i].strip().split())

        # Read operator matrices + irrep matrices
        chars = []
        all_matrices = []  # flat: op0_row0, op0_row1, ..., op1_row0, ...
        for _op in range(opcount):
            # Advance to operator matrix line
            i += 1
            if i >= len(lines):
                break

            # Skip operator matrix (16 ints) and move to the irrep/irtranslation line
            i += 1
            if i >= len(lines):
                break

            # Check if this line is irtranslation (4 ints, non-special k only).
            # Irtranslation is written as '(4i3)' → exactly 4 small integers,
            # e.g. "  0  0  0  1".  1D irrep matrix values are single numbers
            # (e.g. "1" or "-1") and must NOT be mistaken for irtranslation.
            nxt = lines[i].strip().split()
            if (nxt and len(nxt) == 4
                    and not any('.' in x for x in nxt)
                    and all(is_int_4(x) for x in nxt)):
                # Verify it looks like irtranslation: last value (denominator) > 0
                try:
                    if int(nxt[3]) > 0:
                        i += 1  # skip irtranslation
                        if i >= len(lines):
                            break
                except ValueError:
                    pass

            # Read irrep matrix: dim^2 float values (raw floats in the file)
            matrix_vals = []
            needed = dim * dim
            while len(matrix_vals) < needed and i < len(lines):
                for token in lines[i].strip().split():
                    try:
                        matrix_vals.append(float(token))
                    except ValueError:
                        pass
                    if len(matrix_vals) >= needed:
                        break
                if len(matrix_vals) >= needed:
                    break
                i += 1

            # Store raw matrix values (keep full precision)
            all_matrices.extend(matrix_vals)

            # Compute trace (sum of diagonal elements)
            trace = 0.0
            for d in range(dim):
                idx = d * dim + d
                if idx < len(matrix_vals):
                    trace += matrix_vals[idx]
            chars.append(_round_char(trace))

        key = (sg, label)
        if key not in chars_map:
            chars_map[key] = chars
            matrices_map[key] = all_matrices

        # Advance to next irrep header
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if nxt and nxt.split()[0].isdigit() and '"' in nxt:
                break
            i += 1

    return chars_map, matrices_map


# ── CIR (Complex Irreducible Representations) parsing ────────────────────────

def _read_cir_lines():
    """Read CIR_data.txt from the zip archive."""
    zip_path = os.path.join(ISO_DIR, "CIR_data.zip")
    zf = zipfile.ZipFile(zip_path)
    with zf.open("CIR_data.txt") as f:
        return io.TextIOWrapper(f).readlines()


def _parse_complex(s):
    """Parse a complex number in (real,imag) format. Returns (real, imag)."""
    m = re.match(r'\(([^,]+),([^)]+)\)', s)
    if m:
        return float(m.group(1)), float(m.group(2))
    # Try plain float
    try:
        return float(s), 0.0
    except ValueError:
        return 0.0, 0.0


def _parse_cir_characters(needed_labels=None):
    """Parse CIR_data.txt and return character + matrix data.

    Args:
        needed_labels: optional set of (sg, label) to restrict matrix parsing.
                       If None, parses everything. Characters are always parsed.

    Returns:
        cir_chars: dict (sg, label) -> {'dim', 'opcount', 'chars': [(re,im,rounded_re)]}
        cir_matrices: dict (sg, label) -> flattened list of (real, imag) pairs
    """
    lines = _read_cir_lines()

    def _round_char(x, eps=1e-8):
        if abs(x) < eps:
            return 0.0
        r = round(x)
        if abs(x - r) < eps:
            return float(r)
        for n in range(-12, 13):
            for d in (2, 3, 4, 6, 8):
                if abs(x - n / d) < eps:
                    return n / d
        return round(x, 10)

    cir_chars = {}
    cir_matrices = {}  # (sg, label) -> [(re, im), ...] flattened
    i = 3  # skip 3 header lines

    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        tokens = s.split()
        if not tokens or not tokens[0].isdigit():
            i += 1
            continue

        quoted = re.findall(r'"([^"]*)"', s)
        if len(quoted) < 2:
            i += 1
            continue

        sg = int(tokens[1])
        label = quoted[1].strip()
        after_line = s[s.rindex('"') + 1:].strip().split()
        if len(after_line) < 5:
            i += 1
            continue

        dim = int(after_line[0])
        pmkcount = int(after_line[3])
        opcount = int(after_line[4])

        # Skip k-vector data
        kvec_remaining = pmkcount * 16
        while kvec_remaining > 0 and i < len(lines) - 1:
            i += 1
            if i >= len(lines):
                break
            kvec_remaining -= len(lines[i].strip().split())

        # Read operator matrices + complex irrep matrices
        chars = []
        all_matrices = []  # flattened complex matrix elements for all ops
        store_matrices = (needed_labels is None) or ((sg, label) in needed_labels)

        for _op in range(opcount):
            # Advance to operator matrix line (skip it)
            i += 1
            if i >= len(lines):
                break
            # Move past operator matrix to the complex irrep matrix line
            i += 1
            if i >= len(lines):
                break

            # CIR has no irtranslation
            # Read complex matrix: dim² complex numbers as (real,imag) pairs
            complex_vals = []  # list of (real, imag)
            needed = dim * dim
            while len(complex_vals) < needed and i < len(lines):
                for token in lines[i].strip().split():
                    re_val, im_val = _parse_complex(token)
                    complex_vals.append((re_val, im_val))
                    if len(complex_vals) >= needed:
                        break
                if len(complex_vals) >= needed:
                    break
                i += 1

            if store_matrices:
                all_matrices.extend(complex_vals)

            # Compute complex trace
            trace_re = 0.0
            trace_im = 0.0
            for d in range(dim):
                idx = d * dim + d
                if idx < len(complex_vals):
                    trace_re += complex_vals[idx][0]
                    trace_im += complex_vals[idx][1]

            chars.append((trace_re, trace_im, _round_char(trace_re)))

        key = (sg, label)
        if key not in cir_chars:
            cir_chars[key] = {
                'dim': dim,
                'opcount': opcount,
                'chars': chars,  # list of (re, im, rounded_re)
            }
            if store_matrices:
                cir_matrices[key] = all_matrices

        # Advance to next irrep header
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if nxt and nxt.split()[0].isdigit() and '"' in nxt:
                break
            i += 1

    return cir_chars, cir_matrices


def _build_real_matrix(cir_matrices, sg, parts):
    """Build the real matrix for a compound ISO irrep from CIR complex matrices.

    Args:
        cir_matrices: dict (sg, label) -> flattened list of (re, im) tuples
        sg: space group number
        parts: list of CIR label strings, e.g. ['H2', 'H3'] or ['R1', 'R1']

    Returns:
        flat list of f64 real matrix elements, or [] on failure
    """
    # Get complex matrices for each part
    cmats = []
    dims = []
    for part in parts:
        key = (sg, part)
        if key not in cir_matrices:
            return []
        cmats.append(cir_matrices[key])
        # Find dimension from cir_chars (passed separately or infer from matrix)
        # We infer dim from the stored matrix data: total elements / opcount
        # But we don't know opcount here. Let's figure dim from sqrt.
        # Actually we need cir_chars too. Let's compute dim from cir_matrices size.

    # Need cir_data for dim/opcount. We'll pass cir_data too.
    return []  # placeholder, will be reimplemented below


def _build_real_matrix_full(cir_data, cir_matrices, sg, ml):
    """Build real irrep matrix from CIR complex matrices.

    Handles two cases:
    1. Compound label D1D2 → D1 + D2 (different parts → conjugate pair)
       Real matrix = [[A, -B], [B, A]] where D1 = A+iB, D2 = A-iB
    2. Compound label R1R1 → 2×R1 (same part → self-pair)
       If R1 is essentially real (|imag| < eps): block diag [A, 0; 0, A]
       If R1 is complex: realify each diagonal element individually
    """
    key = (sg, ml)
    if key in cir_matrices and key in cir_data:
        # Exact match: single CIR entry
        entry = cir_data[key]
        cmat = cir_matrices[key]
        dim = entry['dim']
        opcount = entry['opcount']
        # Check if mostly real
        max_imag = max(abs(v[1]) for v in cmat) if cmat else 0.0
        if max_imag < 1e-10:
            # Already real, just extract real parts
            return [v[0] for v in cmat]
        else:
            # Complex single irrep: realify by building [[A,-B],[B,A]]
            return _realify_complex_matrix(cmat, dim, opcount)

    # Try decomposition
    parts = _decompose_compound_label(ml)
    if not parts or len(parts) < 2:
        return []

    # Collect matrices for each part
    part_mats = []
    part_dims = []
    for part in parts:
        pk = (sg, part)
        if pk not in cir_matrices or pk not in cir_data:
            return []
        part_mats.append(cir_matrices[pk])
        part_dims.append(cir_data[pk]['dim'])

    # Verify all parts have same opcount
    opcounts = [cir_data[(sg, p)]['opcount'] for p in parts]
    if len(set(opcounts)) != 1:
        return []
    opcount = opcounts[0]

    if parts[0] == parts[1]:
        # Same-label case: R1R1 → 2 copies of realified R1
        cmat = part_mats[0]
        dim = part_dims[0]
        # Check if essentially real
        max_imag = max(abs(v[1]) for v in cmat) if cmat else 0.0
        if max_imag < 1e-10:
            # Block diagonal: [[A,0],[0,A]]
            real_parts = [v[0] for v in cmat]
            return _block_diag_two(real_parts, dim, opcount)
        else:
            # Realify complex matrix then duplicate
            realified = _realify_complex_matrix(cmat, dim, opcount)
            # Now duplicate to 2×2 blocks
            new_dim = dim * 2  # realified is 2*dim, duplicate to 4*dim
            # Actually realified is already 2*dim. We need two copies.
            return _block_diag_two(realified, dim * 2, opcount)
    else:
        # Different-label case: H2H3 → H2 + H3 (conjugate pair)
        # For each operator, build [[A,-B],[B,A]] from H2 = A+iB
        cmat1 = part_mats[0]  # H2
        dim1 = part_dims[0]
        # Realify the first component (H2), assuming H3 is conjugate
        realified = _realify_complex_matrix(cmat1, dim1, opcount)
        # H2H3 is the direct sum of H2 (real form) and H3 (same real form, since H3=H2*)
        # So H2H3 = [[A,-B, 0, 0], [B, A, 0, 0], [0, 0, A, -B], [0, 0, B, A]]
        # Wait no. H2 is already realified as [[A,-B],[B,A]] (2*dim × 2*dim).
        # H3 is the complex conjugate, which has the SAME real form.
        # H2H3 = H2 ⊕ H3 = realified(H2) ⊕ realified(H2) = block_diag(realified, realified)
        return _block_diag_two(realified, dim1 * 2, opcount)


def _realify_complex_matrix(cmat, dim, opcount):
    """Convert a complex matrix to real block form [[A,-B],[B,A]].

    For each operator: the complex dim×dim matrix M = A + iB
    becomes a 2*dim × 2*dim real matrix:
        [ A, -B ]
        [ B,  A ]
    where A = Re(M), B = Im(M).
    """
    per_op = dim * dim
    result = []
    for op_idx in range(opcount):
        start = op_idx * per_op
        # Extract A and B from complex values
        a_vals = []  # dim×dim real
        b_vals = []  # dim×dim real
        for row in range(dim):
            for col in range(dim):
                idx = start + row * dim + col
                a_vals.append(cmat[idx][0])  # real part
                b_vals.append(cmat[idx][1])  # imag part

        new_dim = dim * 2
        # Build new_d × new_d matrix row by row
        for row in range(new_dim):
            for col in range(new_dim):
                if row < dim and col < dim:
                    # Top-left: A
                    result.append(a_vals[row * dim + col])
                elif row < dim and col >= dim:
                    # Top-right: -B
                    result.append(-b_vals[row * dim + (col - dim)])
                elif row >= dim and col < dim:
                    # Bottom-left: B
                    result.append(b_vals[(row - dim) * dim + col])
                else:
                    # Bottom-right: A
                    result.append(a_vals[(row - dim) * dim + (col - dim)])
    return result


def _block_diag_two(mat, dim, opcount):
    """Create block-diagonal matrix with two copies of mat on the diagonal.

    For each operator: result = [[mat, 0], [0, mat]]
    """
    per_op = dim * dim
    result = []
    for op_idx in range(opcount):
        start = op_idx * per_op
        new_dim = dim * 2
        for row in range(new_dim):
            for col in range(new_dim):
                if row < dim and col < dim:
                    # Top-left block
                    idx = row * dim + col
                    result.append(mat[start + idx])
                elif row >= dim and col >= dim:
                    # Bottom-right block
                    idx = (row - dim) * dim + (col - dim)
                    result.append(mat[start + idx])
                else:
                    # Zero off-diagonal blocks
                    result.append(0.0)
    return result


def _decompose_compound_label(ml):
    """Split a compound irrep label like 'D1D2' into ['D1', 'D2'].

    Handles patterns like:
      'D1D2'   → ['D1', 'D2']
      'H2H3'   → ['H2', 'H3']
      'A1A2'   → ['A1', 'A2']
      'R1R1'   → ['R1', 'R1']
      'M2+M3+' → ['M2+', 'M3+']
      'A2+A3+' → ['A2+', 'A3+']
      'M3+M4+' → ['M3+', 'M4+']
      'W2W2'   → ['W2', 'W2']
      'M1M2'   → ['M1', 'M2']
      'GM3+GM4+' → ['GM3+', 'GM4+']
    """
    # Try splitting on '+' or '-' boundaries
    # Pattern: letter(s) + digits + optional sign, repeated
    parts = re.findall(r'([A-Za-z]+\d+[+-]?)', ml)
    if len(parts) >= 2:
        # Verify the parts concatenate back to the original
        if ''.join(parts) == ml:
            return parts

    # Fallback: try without signs
    parts = re.findall(r'([A-Za-z]+\d+)', ml)
    if len(parts) >= 2 and ''.join(parts) == ml.rstrip('+-'):
        # Re-attach signs
        signs = re.findall(r'[+-]', ml)
        result = []
        for i, p in enumerate(parts):
            if i < len(signs):
                result.append(p + signs[i])
            else:
                result.append(p)
        if len(result) >= 2:
            return result

    return None  # can't decompose


def _lookup_cir_chars(cir_data, sg_num, ml_label):
    """Look up character data from CIR, handling compound labels.

    For compound labels like 'H2H3' = H2 ⊕ H3, the character is
    χ(H2H3) = 2 * Re(χ_CIR(H2)) assuming H3 is conjugate of H2.
    More generally, χ = sum of Re(χ) for each component.
    """
    # 1. Exact match in CIR
    key = (sg_num, ml_label)
    if key in cir_data:
        entry = cir_data[key]
        # Return real parts of complex characters
        return [c[2] for c in entry['chars']]

    # 2. Try decomposing compound label
    parts = _decompose_compound_label(ml_label)
    if parts and len(parts) >= 2:
        all_chars = []
        for part in parts:
            pk = (sg_num, part)
            if pk in cir_data:
                all_chars.append([c[2] for c in cir_data[pk]['chars']])
            else:
                # Part not found in CIR
                return []

        if len(all_chars) == len(parts):
            # Sum characters component-wise
            if all_chars:
                n_ops = len(all_chars[0])
                if all(len(ch) == n_ops for ch in all_chars):
                    summed = []
                    for op_idx in range(n_ops):
                        total = sum(ch[op_idx] for ch in all_chars)
                        summed.append(_round_char_cir(total))
                    return summed

    return []


def _round_char_cir(x, eps=1e-8):
    """Round character value (same as _round_char but accessible here)."""
    if abs(x) < eps:
        return 0.0
    r = round(x)
    if abs(x - r) < eps:
        return float(r)
    for n in range(-12, 13):
        for d in (2, 3, 4, 6, 8):
            if abs(x - n / d) < eps:
                return n / d
    return round(x, 10)


# ── main parsing ─────────────────────────────────────────────────────────────

def parse_all():
    print("Parsing data_irreps.txt...")
    irr_lines = read_file("data_irreps.txt")
    irr_sec = get_sections(irr_lines)

    ml_labels  = parse_labels(irr_lines, irr_sec, "irrep_label")
    bc_labels  = parse_labels(irr_lines, irr_sec, "irrep_label_bc")
    kov_labels = parse_labels(irr_lines, irr_sec, "irrep_label_kov")
    zak_labels = parse_labels(irr_lines, irr_sec, "irrep_label_zak")
    sg_numbers = parse_ints(irr_lines, irr_sec, "irrep_space_group")
    images     = parse_ints(irr_lines, irr_sec, "irrep_image")
    lifshitz   = parse_ints(irr_lines, irr_sec, "irrep_lifshitz")
    n_irreps = len(ml_labels)
    print(f"  {n_irreps} irreps, {len(bc_labels)} BC, {len(kov_labels)} Kov, {len(sg_numbers)} SG nums")

    # Ensure all arrays have same length
    assert len(bc_labels) == n_irreps, f"BC: {len(bc_labels)} != {n_irreps}"
    assert len(kov_labels) >= n_irreps, f"Kov: {len(kov_labels)} < {n_irreps}"
    assert len(sg_numbers) == n_irreps
    assert len(images) == n_irreps
    assert len(lifshitz) == n_irreps

    print("Parsing data_images.txt...")
    img_lines = read_file("data_images.txt")
    img_sec = get_sections(img_lines)
    img_labels = parse_labels(img_lines, img_sec, "image_label")
    # image code 1 → img_labels[0], code 2 → img_labels[1], etc.
    print(f"  {len(img_labels)} image labels")

    print("Parsing data_little.txt...")
    lit_lines = read_file("data_little.txt")
    lit_sec = get_sections(lit_lines)
    k_counts = parse_ints(lit_lines, lit_sec, "little_k_count")
    k_labels_all = parse_labels(lit_lines, lit_sec, "little_k_label")
    k_kov_all = parse_ints(lit_lines, lit_sec, "little_k_kov")
    print(f"  {len(k_counts)} SGs, {len(k_labels_all)} k-labels, {len(k_kov_all)} k-kov")

    print("Parsing data_isotropy.txt...")
    iso_lines = read_file("data_isotropy.txt")
    iso_sec = get_sections(iso_lines)

    iso_irrep       = parse_ints(iso_lines, iso_sec, "isotropy_irrep")
    iso_irrep_ptr   = parse_ints(iso_lines, iso_sec, "isotropy_irrep_pointer")
    iso_subgroups   = parse_ints(iso_lines, iso_sec, "isotropy_subgroup")
    iso_direction   = parse_ints(iso_lines, iso_sec, "isotropy_direction")
    iso_domains     = parse_ints(iso_lines, iso_sec, "isotropy_domain_count")
    iso_domain_type = parse_ints(iso_lines, iso_sec, "isotropy_domain_type_count")
    iso_arms        = parse_ints(iso_lines, iso_sec, "isotropy_arms")
    iso_order       = parse_ints(iso_lines, iso_sec, "isotropy_order")
    iso_ferroic     = parse_ints(iso_lines, iso_sec, "isotropy_ferroic")

    # basis and origin are floats
    iso_basis_raw   = parse_floats(iso_lines, iso_sec, "isotropy_basis")
    iso_origin_raw  = parse_floats(iso_lines, iso_sec, "isotropy_origin")

    # direction labels, dimension, and free parameter count for direction mapping
    iso_dir_labels  = parse_labels(iso_lines, iso_sec, "isotropy_orderparam_label")
    iso_dir_dim     = parse_ints(iso_lines, iso_sec, "isotropy_orderparam_dim")
    iso_dir_free    = parse_ints(iso_lines, iso_sec, "isotropy_orderparam_freeparam")

    print(f"  {len(iso_irrep)} iso entries, {len(iso_subgroups)} subgroups, {len(iso_irrep_ptr)} ptrs")
    print(f"  {len(iso_dir_labels)} direction labels")
    if iso_dir_labels:
        print(f"  direction label[1]={iso_dir_labels[1] if len(iso_dir_labels)>1 else 'N/A'}")
    print(f"  {len(iso_dir_dim)} dir dims, {len(iso_dir_free)} dir free params")

    # Build direction lookup using the comprehensive mapping
    dir_map = build_direction_map(iso_direction, iso_dir_dim, iso_dir_free, iso_dir_labels)
    print(f"  Built direction map with {len(dir_map)} entries")

    print("Parsing data_magnetic.txt (magnetic isotropy subgroups)...")
    mag_lines = read_file("data_magnetic.txt")
    mag_sec = get_sections(mag_lines)

    mag_iso_sg       = parse_ints(mag_lines, mag_sec, "mag_iso_subgroup")
    mag_iso_irrep    = parse_ints(mag_lines, mag_sec, "mag_iso_irrep")
    mag_iso_ptr      = parse_ints(mag_lines, mag_sec, "mag_iso_irrep_pointer")
    mag_nlabel       = parse_labels(mag_lines, mag_sec, "mag_nlabel")
    mag_bns_label    = parse_labels(mag_lines, mag_sec, "mag_bns_label")
    print(f"  {len(mag_iso_sg)} mag iso entries, {len(mag_iso_ptr)} ptrs, {len(mag_nlabel)} labels")

    # Direction labels for magnetic isotropy
    mag_iso_dir_labels = parse_labels(mag_lines, mag_sec, "mag_iso_orderparam_label")
    # Map direction codes to labels (similar to non-mag dir_map)
    mag_iso_dir_code  = parse_ints(mag_lines, mag_sec, "mag_iso_orderparam")
    mag_iso_dir_ptr   = parse_ints(mag_lines, mag_sec, "mag_iso_orderparam_pointer")
    # Build mag direction lookup: entry index → direction string
    mag_dir_by_entry = {}
    for entry_idx in range(len(mag_iso_sg)):
        if entry_idx < len(mag_iso_dir_ptr) and mag_iso_dir_ptr[entry_idx] > 0:
            ptr = mag_iso_dir_ptr[entry_idx] - 1  # 1-based → 0-based
            if ptr < len(mag_iso_dir_labels):
                mag_dir_by_entry[entry_idx] = mag_iso_dir_labels[ptr]
            else:
                mag_dir_by_entry[entry_idx] = f"dir{ptr}"
        else:
            mag_dir_by_entry[entry_idx] = "(a)"
    print(f"  {len(mag_dir_by_entry)} direction labels mapped")

    print("Parsing data_space.txt...")
    sp_lines = read_file("data_space.txt")
    sp_sec = get_sections(sp_lines)
    # We need SG number → symbol mapping
    # data_space has many sections; the key one for symbols
    sg_symbol_map = {}
    # Build from the preamble.rs SPACEGROUP_INDEX instead
    # (already has all 230 SG symbols)

    print("Parsing PIR_data.txt k-vectors...")
    kvec_map = _parse_pir_kvectors()
    print(f"  Parsed {len(kvec_map)} k-vector entries")

    print("Parsing PIR_data.txt characters and matrices...")
    chars_map, matrices_map = _parse_pir_characters()
    print(f"  Parsed {len(chars_map)} character table entries")
    print(f"  Parsed {len(matrices_map)} matrix data entries")

    # Determine which CIR labels need matrix data (missing from PIR)
    needed_cir = set()
    for i in range(len(ml_labels)):
        mm = _lookup_matrices(matrices_map, sg_numbers[i], ml_labels[i], kvec_map)
        if not mm:
            # This label needs CIR fallback
            needed_cir.add((sg_numbers[i], ml_labels[i]))
            # Also add decomposed parts
            parts = _decompose_compound_label(ml_labels[i])
            if parts:
                for p in parts:
                    needed_cir.add((sg_numbers[i], p))

    print(f"Parsing CIR_data.txt (fallback for {len(needed_cir)} needed labels)...")
    cir_data, cir_matrices = _parse_cir_characters(needed_labels=needed_cir)
    print(f"  Parsed {len(cir_data)} CIR character entries, {len(cir_matrices)} matrix entries")

    print("Parsing spinor (double-valued) irrep data from irrepTables...")
    try:
        from parse_spinor_data import parse_all_spinor
        spinor_irreps = parse_all_spinor()
        print(f"  Parsed {len(spinor_irreps)} spinor irreps")
    except (ImportError, FileNotFoundError) as e:
        print(f"  Skipping spinor data: {e}")
        spinor_irreps = []

    return {
        "n_irreps": n_irreps + len(spinor_irreps),
        "spinor_irreps": spinor_irreps,
        "ml_labels": ml_labels,
        "bc_labels": bc_labels,
        "kov_labels": kov_labels,
        "zak_labels": zak_labels,
        "sg_numbers": sg_numbers,
        "images": images,
        "lifshitz": lifshitz,
        "img_labels": img_labels,
        "k_counts": k_counts,
        "k_labels_all": k_labels_all,
        "k_kov_all": k_kov_all,
        "iso_irrep": iso_irrep,
        "iso_irrep_ptr": iso_irrep_ptr,
        "iso_subgroups": iso_subgroups,
        "iso_direction": iso_direction,
        "iso_domains": iso_domains,
        "iso_domain_type": iso_domain_type,
        "iso_arms": iso_arms,
        "iso_order": iso_order,
        "iso_basis_raw": iso_basis_raw,
        "iso_origin_raw": iso_origin_raw,
        "iso_ferroic": iso_ferroic,
        "dir_map": dir_map,
        "kvec_map": kvec_map,
        "chars_map": chars_map,
        "matrices_map": matrices_map,
        "cir_data": cir_data,
        "cir_matrices": cir_matrices,
        "mag_iso_sg": mag_iso_sg,
        "mag_iso_irrep": mag_iso_irrep,
        "mag_iso_ptr": mag_iso_ptr,
        "mag_nlabel": mag_nlabel,
        "mag_bns_label": mag_bns_label,
        "mag_dir_by_entry": mag_dir_by_entry,
    }

# ── data assembly ────────────────────────────────────────────────────────────

CRYSTAL_SYSTEMS = {
    "triclinic":    range(1, 3),
    "monoclinic":   range(3, 16),
    "orthorhombic": range(16, 75),
    "tetragonal":   range(75, 143),
    "trigonal":     range(143, 168),
    "hexagonal":    range(168, 195),
    "cubic":        range(195, 231),
}

# SG symbol lookup (from preamble.rs data)
SG_SYMBOLS = {
    1: ("P1", "C1^1"), 2: ("P-1", "Ci^1"),
    3: ("P2", "C2^1"), 4: ("P2_1", "C2^2"), 5: ("C2", "C2^3"),
    6: ("Pm", "Cs^1"), 7: ("Pc", "Cs^2"), 8: ("Cm", "Cs^3"), 9: ("Cc", "Cs^4"),
    10: ("P2/m", "C2h^1"), 11: ("P2_1/m", "C2h^2"), 12: ("C2/m", "C2h^3"),
    13: ("P2/c", "C2h^4"), 14: ("P2_1/c", "C2h^5"), 15: ("C2/c", "C2h^6"),
    16: ("P222", "D2^1"), 17: ("P222_1", "D2^2"), 18: ("P2_12_12", "D2^3"),
    19: ("P2_12_12_1", "D2^4"), 20: ("C222_1", "D2^5"), 21: ("C222", "D2^6"),
    22: ("F222", "D2^7"), 23: ("I222", "D2^8"), 24: ("I2_12_12_1", "D2^9"),
    25: ("Pmm2", "C2v^1"), 26: ("Pmc2_1", "C2v^2"), 27: ("Pcc2", "C2v^3"),
    28: ("Pma2", "C2v^4"), 29: ("Pca2_1", "C2v^5"), 30: ("Pnc2", "C2v^6"),
    31: ("Pmn2_1", "C2v^7"), 32: ("Pba2", "C2v^8"), 33: ("Pna2_1", "C2v^9"),
    34: ("Pnn2", "C2v^10"), 35: ("Cmm2", "C2v^11"), 36: ("Cmc2_1", "C2v^12"),
    37: ("Ccc2", "C2v^13"), 38: ("Amm2", "C2v^14"), 39: ("Abm2", "C2v^15"),
    40: ("Ama2", "C2v^16"), 41: ("Aba2", "C2v^17"), 42: ("Fmm2", "C2v^18"),
    43: ("Fdd2", "C2v^19"), 44: ("Imm2", "C2v^20"), 45: ("Iba2", "C2v^21"),
    46: ("Ima2", "C2v^22"), 47: ("Pmmm", "D2h^1"), 48: ("Pnnn", "D2h^2"),
    49: ("Pccm", "D2h^3"), 50: ("Pban", "D2h^4"), 51: ("Pmma", "D2h^5"),
    52: ("Pnna", "D2h^6"), 53: ("Pmna", "D2h^7"), 54: ("Pcca", "D2h^8"),
    55: ("Pbam", "D2h^9"), 56: ("Pccn", "D2h^10"), 57: ("Pbcm", "D2h^11"),
    58: ("Pnnm", "D2h^12"), 59: ("Pmmn", "D2h^13"), 60: ("Pbcn", "D2h^14"),
    61: ("Pbca", "D2h^15"), 62: ("Pnma", "D2h^16"), 63: ("Cmcm", "D2h^17"),
    64: ("Cmca", "D2h^18"), 65: ("Cmmm", "D2h^19"), 66: ("Cccm", "D2h^20"),
    67: ("Cmma", "D2h^21"), 68: ("Ccca", "D2h^22"), 69: ("Fmmm", "D2h^23"),
    70: ("Fddd", "D2h^24"), 71: ("Immm", "D2h^25"), 72: ("Ibam", "D2h^26"),
    73: ("Ibca", "D2h^27"), 74: ("Imma", "D2h^28"),
    75: ("P4", "C4^1"), 76: ("P4_1", "C4^2"), 77: ("P4_2", "C4^3"),
    78: ("P4_3", "C4^4"), 79: ("I4", "C4^5"), 80: ("I4_1", "C4^6"),
    81: ("P-4", "S4^1"), 82: ("I-4", "S4^2"), 83: ("P4/m", "C4h^1"),
    84: ("P4_2/m", "C4h^2"), 85: ("P4/n", "C4h^3"), 86: ("P4_2/n", "C4h^4"),
    87: ("I4/m", "C4h^5"), 88: ("I4_1/a", "C4h^6"), 89: ("P422", "D4^1"),
    90: ("P42_12", "D4^2"), 91: ("P4_122", "D4^3"), 92: ("P4_12_12", "D4^4"),
    93: ("P4_222", "D4^5"), 94: ("P4_22_12", "D4^6"), 95: ("P4_322", "D4^7"),
    96: ("P4_32_12", "D4^8"), 97: ("I422", "D4^9"), 98: ("I4_122", "D4^10"),
    99: ("P4mm", "C4v^1"), 100: ("P4bm", "C4v^2"), 101: ("P4_2cm", "C4v^3"),
    102: ("P4_2nm", "C4v^4"), 103: ("P4cc", "C4v^5"), 104: ("P4nc", "C4v^6"),
    105: ("P4_2mc", "C4v^7"), 106: ("P4_2bc", "C4v^8"), 107: ("I4mm", "C4v^9"),
    108: ("I4cm", "C4v^10"), 109: ("I4_1md", "C4v^11"), 110: ("I4_1cd", "C4v^12"),
    111: ("P-42m", "D2d^1"), 112: ("P-42c", "D2d^2"), 113: ("P-42_1m", "D2d^3"),
    114: ("P-42_1c", "D2d^4"), 115: ("P-4m2", "D2d^5"), 116: ("P-4c2", "D2d^6"),
    117: ("P-4b2", "D2d^7"), 118: ("P-4n2", "D2d^8"), 119: ("I-4m2", "D2d^9"),
    120: ("I-4c2", "D2d^10"), 121: ("I-42m", "D2d^11"), 122: ("I-42d", "D2d^12"),
    123: ("P4/mmm", "D4h^1"), 124: ("P4/mcc", "D4h^2"), 125: ("P4/nbm", "D4h^3"),
    126: ("P4/nnc", "D4h^4"), 127: ("P4/mbm", "D4h^5"), 128: ("P4/mnc", "D4h^6"),
    129: ("P4/nmm", "D4h^7"), 130: ("P4/ncc", "D4h^8"), 131: ("P4_2/mmc", "D4h^9"),
    132: ("P4_2/mcm", "D4h^10"), 133: ("P4_2/nbc", "D4h^11"), 134: ("P4_2/nnm", "D4h^12"),
    135: ("P4_2/mbc", "D4h^13"), 136: ("P4_2/mnm", "D4h^14"), 137: ("P4_2/nmc", "D4h^15"),
    138: ("P4_2/ncm", "D4h^16"), 139: ("I4/mmm", "D4h^17"), 140: ("I4/mcm", "D4h^18"),
    141: ("I4_1/amd", "D4h^19"), 142: ("I4_1/acd", "D4h^20"),
    143: ("P3", "C3^1"), 144: ("P3_1", "C3^2"), 145: ("P3_2", "C3^3"),
    146: ("R3", "C3^4"), 147: ("P-3", "C3i^1"), 148: ("R-3", "C3i^2"),
    149: ("P312", "D3^1"), 150: ("P321", "D3^2"), 151: ("P3_112", "D3^3"),
    152: ("P3_121", "D3^4"), 153: ("P3_212", "D3^5"), 154: ("P3_221", "D3^6"),
    155: ("R32", "D3^7"), 156: ("P3m1", "C3v^1"), 157: ("P31m", "C3v^2"),
    158: ("P3c1", "C3v^3"), 159: ("P31c", "C3v^4"), 160: ("R3m", "C3v^5"),
    161: ("R3c", "C3v^6"), 162: ("P-31m", "D3d^1"), 163: ("P-31c", "D3d^2"),
    164: ("P-3m1", "D3d^3"), 165: ("P-3c1", "D3d^4"), 166: ("R-3m", "D3d^5"),
    167: ("R-3c", "D3d^6"),
    168: ("P6", "C6^1"), 169: ("P6_1", "C6^2"), 170: ("P6_5", "C6^3"),
    171: ("P6_2", "C6^4"), 172: ("P6_4", "C6^5"), 173: ("P6_3", "C6^6"),
    174: ("P-6", "C3h^1"), 175: ("P6/m", "C6h^1"), 176: ("P6_3/m", "C6h^2"),
    177: ("P622", "D6^1"), 178: ("P6_122", "D6^2"), 179: ("P6_522", "D6^3"),
    180: ("P6_222", "D6^4"), 181: ("P6_422", "D6^5"), 182: ("P6_322", "D6^6"),
    183: ("P6mm", "C6v^1"), 184: ("P6cc", "C6v^2"), 185: ("P6_3cm", "C6v^3"),
    186: ("P6_3mc", "C6v^4"), 187: ("P-6m2", "D3h^1"), 188: ("P-6c2", "D3h^2"),
    189: ("P-62m", "D3h^3"), 190: ("P-62c", "D3h^4"), 191: ("P6/mmm", "D6h^1"),
    192: ("P6/mcc", "D6h^2"), 193: ("P6_3/mcm", "D6h^3"), 194: ("P6_3/mmc", "D6h^4"),
    195: ("P23", "T^1"), 196: ("F23", "T^2"), 197: ("I23", "T^3"),
    198: ("P2_13", "T^4"), 199: ("I2_13", "T^5"), 200: ("Pm-3", "Th^1"),
    201: ("Pn-3", "Th^2"), 202: ("Fm-3", "Th^3"), 203: ("Fd-3", "Th^4"),
    204: ("Im-3", "Th^5"), 205: ("Pa-3", "Th^6"), 206: ("Ia-3", "Th^7"),
    207: ("P432", "O^1"), 208: ("P4_232", "O^2"), 209: ("F432", "O^3"),
    210: ("F4_132", "O^4"), 211: ("I432", "O^5"), 212: ("P4_332", "O^6"),
    213: ("P4_132", "O^7"), 214: ("I4_132", "O^8"), 215: ("P-43m", "Td^1"),
    216: ("F-43m", "Td^2"), 217: ("I-43m", "Td^3"), 218: ("P-43n", "Td^4"),
    219: ("F-43c", "Td^5"), 220: ("I-43d", "Td^6"), 221: ("Pm-3m", "Oh^1"),
    222: ("Pn-3n", "Oh^2"), 223: ("Pm-3n", "Oh^3"), 224: ("Pn-3m", "Oh^4"),
    225: ("Fm-3m", "Oh^5"), 226: ("Fm-3c", "Oh^6"), 227: ("Fd-3m", "Oh^7"),
    228: ("Fd-3c", "Oh^8"), 229: ("Im-3m", "Oh^9"), 230: ("Ia-3d", "Oh^10"),
}

def get_sg_symbol(sg_num):
    """Get HM symbol and Schoenflies for a space group number."""
    if sg_num in SG_SYMBOLS:
        return SG_SYMBOLS[sg_num]
    return (f"SG{sg_num}", "")

def get_crystal_system(sg_num):
    for name, rng in CRYSTAL_SYSTEMS.items():
        if sg_num in rng:
            return name
    return "unknown"

# ── Rust code generation ─────────────────────────────────────────────────────

def escape_rust_str(s):
    """Escape a string for use in a Rust &'static str."""
    return s.replace("\\", "\\\\").replace("\"", "\\\"")

def latex_escape(s):
    """Escape a string for LaTeX math mode."""
    return s

def _lookup_kvec(kvec_map, sg_num, ml_label):
    """Look up k-vector for (SG#, ML_label) with fallback matching.

    First tries exact match. If that fails, tries progressively shorter labels
    for compound labels like "H2H3" or "P1P1" that may appear slightly differently
    in PIR_data.txt.
    """
    # 1. Exact match
    key = (sg_num, ml_label)
    if key in kvec_map:
        return kvec_map[key]

    # 2. Try without trailing characters for compound labels
    # e.g. "P1P1" might be "P1PA1" in PIR, try "P1" as prefix
    # Strip digits+signs from the end
    import re
    body = ml_label
    # Keep stripping trailing pieces until we find a match or the body is too short
    while len(body) > 2:
        # Try to strip trailing numeric-suffix groups
        # e.g. "H2H3" → "H2" then "H"
        new_body = re.sub(r'[0-9][+-]?$', '', body)
        if new_body == body:
            # No numeric suffix found, try stripping last char
            new_body = body[:-1]
        body = new_body
        if len(body) < 2:
            break
        # Try exact match with this shorter body
        key_short = (sg_num, body)
        if key_short in kvec_map:
            return kvec_map[key_short]
        # Also try as prefix match: find any key with this sg_num whose label
        # starts with body
        for k, v in kvec_map.items():
            if k[0] == sg_num and k[1].startswith(body):
                return v

    # 3. Try prefix matching on the original label
    # Find any kvec entry for this SG whose label starts the same way
    for k, v in kvec_map.items():
        if k[0] == sg_num and (k[1].startswith(ml_label) or ml_label.startswith(k[1])):
            return v

    return (0, 0, 0, 1)  # default: Gamma point


def _lookup_chars(chars_map, sg_num, ml_label, kvec_map=None):
    """Look up character table for (SG#, ML_label) with fallback matching.

    Uses progressively more aggressive matching strategies to handle
    label differences between data_irreps.txt and PIR_data.txt.
    """
    # 1. Exact match
    key = (sg_num, ml_label)
    if key in chars_map:
        return chars_map[key]

    # 2. Try without trailing characters for compound labels
    body = ml_label
    while len(body) > 2:
        new_body = re.sub(r'[0-9][+-]?$', '', body)
        if new_body == body:
            new_body = body[:-1]
        body = new_body
        if len(body) < 2:
            break
        key_short = (sg_num, body)
        if key_short in chars_map:
            return chars_map[key_short]
        for k, v in chars_map.items():
            if k[0] == sg_num and k[1].startswith(body):
                return v

    # 3. Try prefix matching on the original label
    for k, v in chars_map.items():
        if k[0] == sg_num and (k[1].startswith(ml_label) or ml_label.startswith(k[1])):
            return v

    # 4. K-vector based fallback: find PIR entries at the same k-point coords
    #    and match by numeric suffix (handles X/Y/Z ↔ A/B/C convention diffs).
    if kvec_map is not None:
        iso_kvec = _lookup_kvec(kvec_map, sg_num, ml_label)
        if iso_kvec != (0, 0, 0, 1):
            iso_num = _label_number(ml_label)
            iso_sign = '+' if '+' in ml_label else ('-' if '-' in ml_label else None)

            # Collect all PIR chars_map entries at same kvec
            same_kvec = [(k, v) for k, v in chars_map.items()
                         if k[0] == sg_num and kvec_map.get(k) == iso_kvec]

            # 4a. Same kvec + same numeric part
            for pir_key, pir_chars in same_kvec:
                pir_num = _label_number(pir_key[1])
                if iso_num is not None and pir_num is not None and iso_num == pir_num:
                    return pir_chars

            # 4b. Same kvec, only one PIR entry → unambiguous match
            if len(same_kvec) == 1:
                return same_kvec[0][1]

            # 4c. Same kvec + same sign pattern
            if iso_sign and len(same_kvec) > 0:
                sign_matches = [(k, v) for k, v in same_kvec if iso_sign in k[1]]
                if len(sign_matches) == 1:
                    return sign_matches[0][1]

            # 4d. Same kvec + same k-point letter + same sign (numbering offset).
            #     e.g. ISO "W2" vs PIR "W1" at same kvec → different numbering.
            iso_k_label = _kpoint_label_from_ml(ml_label)
            if iso_k_label:
                letter_sign_matches = []
                for pir_key, pir_chars in same_kvec:
                    pir_k_label = _kpoint_label_from_ml(pir_key[1])
                    if pir_k_label == iso_k_label:
                        pir_sign = '+' if '+' in pir_key[1] else ('-' if '-' in pir_key[1] else None)
                        if iso_sign == pir_sign:
                            letter_sign_matches.append((pir_key, pir_chars))
                if len(letter_sign_matches) == 1:
                    return letter_sign_matches[0][1]

    return []  # not found
def _label_number(label):
    """Extract the first numeric part from an irrep label. Returns int or None."""
    m = re.search(r'(\d+)', label)
    return int(m.group(1)) if m else None


def _kpoint_label_from_ml(ml):
    """Extract the k-point letter prefix from a Miller-Love label.
    E.g. 'GM4+' → 'GM', 'X3-' → 'X', 'W2W2' → 'W', 'P2P3' → 'P'."""
    body = ml.strip().rstrip('+-')
    # Find first digit
    m = re.search(r'\d', body)
    if m:
        return body[:m.start()]
    return body


def _lookup_matrices(matrices_map, sg_num, ml_label, kvec_map=None):
    """Look up matrix data for (SG#, ML_label) with the same matching as `_lookup_chars`."""
    # 1. Exact match
    key = (sg_num, ml_label)
    if key in matrices_map:
        return matrices_map[key]

    # 2. Try without trailing characters for compound labels
    body = ml_label
    while len(body) > 2:
        new_body = re.sub(r'[0-9][+-]?$', '', body)
        if new_body == body:
            new_body = body[:-1]
        body = new_body
        if len(body) < 2:
            break
        key_short = (sg_num, body)
        if key_short in matrices_map:
            return matrices_map[key_short]
        for k, v in matrices_map.items():
            if k[0] == sg_num and k[1].startswith(body):
                return v

    # 3. Try prefix matching on the original label
    for k, v in matrices_map.items():
        if k[0] == sg_num and (k[1].startswith(ml_label) or ml_label.startswith(k[1])):
            return v

    # 4. K-vector based fallback
    if kvec_map is not None:
        iso_kvec = _lookup_kvec(kvec_map, sg_num, ml_label)
        if iso_kvec != (0, 0, 0, 1):
            iso_num = _label_number(ml_label)
            iso_sign = '+' if '+' in ml_label else ('-' if '-' in ml_label else None)

            same_kvec = [(k, v) for k, v in matrices_map.items()
                         if k[0] == sg_num and kvec_map.get(k) == iso_kvec]

            # 4a. Same kvec + same numeric part
            for pir_key, pir_mat in same_kvec:
                pir_num = _label_number(pir_key[1])
                if iso_num is not None and pir_num is not None and iso_num == pir_num:
                    return pir_mat

            # 4b. Same kvec, only one PIR entry
            if len(same_kvec) == 1:
                return same_kvec[0][1]

            # 4c. Same kvec + same sign pattern
            if iso_sign and len(same_kvec) > 0:
                sign_matches = [(k, v) for k, v in same_kvec if iso_sign in k[1]]
                if len(sign_matches) == 1:
                    return sign_matches[0][1]

            # 4d. Same kvec + same k-point letter + same sign
            iso_k_label = _kpoint_label_from_ml(ml_label)
            if iso_k_label:
                letter_sign_matches = []
                for pir_key, pir_mat in same_kvec:
                    pir_k_label = _kpoint_label_from_ml(pir_key[1])
                    if pir_k_label == iso_k_label:
                        pir_sign = '+' if '+' in pir_key[1] else ('-' if '-' in pir_key[1] else None)
                        if iso_sign == pir_sign:
                            letter_sign_matches.append((pir_key, pir_mat))
                if len(letter_sign_matches) == 1:
                    return letter_sign_matches[0][1]

    return []  # not found


def generate_rust_data(data):
    """Generate the content of generated_data.rs."""
    ml  = data["ml_labels"]
    bc  = data["bc_labels"]
    kov = data["kov_labels"]
    sg  = data["sg_numbers"]
    img = data["images"]
    lif = data["lifshitz"]
    img_labels = data["img_labels"]

    # direction labels
    dir_map = data["dir_map"]
    # k-vector map: (SG#, ML_label) -> (kx, ky, kz, denom)
    kvec_map = data["kvec_map"]
    # character map: (SG#, ML_label) -> [char1, char2, ...]
    chars_map = data.get("chars_map", {})
    # matrix map: (SG#, ML_label) -> [flat_values...]
    matrices_map = data.get("matrices_map", {})

    cir_data = data.get("cir_data", {})

    # ── Build flat CHARACTERS array and per-irrep start/count ──
    chars_flat = []
    char_starts = []
    char_counts = []
    missing_chars = 0
    cir_filled = 0
    for i in range(len(ml)):
        ch = _lookup_chars(chars_map, sg[i], ml[i], kvec_map)
        if not ch and cir_data:
            # Fallback to CIR data
            ch = _lookup_cir_chars(cir_data, sg[i], ml[i])
            if ch:
                cir_filled += 1
        char_starts.append(len(chars_flat))
        char_counts.append(len(ch))
        chars_flat.extend(ch)
        if not ch:
            missing_chars += 1
    if missing_chars > 0:
        print(f"  Warning: {missing_chars}/{len(ml)} irreps have no character data")
    if cir_filled > 0:
        print(f"  (CIR fallback filled {cir_filled} character tables)")

    # ── Build flat MATRICES array and per-irrep start/count ──
    cir_mat = data.get("cir_matrices", {})
    cir_d = data.get("cir_data", {})
    matrices_flat = []
    mat_starts = []
    mat_counts = []
    missing_mat = 0
    cir_mat_filled = 0
    for i in range(len(ml)):
        mm = _lookup_matrices(matrices_map, sg[i], ml[i], kvec_map)
        if not mm and cir_mat:
            # Try to build real matrix from CIR data
            mm = _build_real_matrix_full(cir_d, cir_mat, sg[i], ml[i])
            if mm:
                cir_mat_filled += 1
        mat_starts.append(len(matrices_flat))
        mat_counts.append(len(mm))
        matrices_flat.extend(mm)
        if not mm:
            missing_mat += 1
    if missing_mat > 0:
        print(f"  Warning: {missing_mat}/{len(ml)} irreps have no matrix data")
    if cir_mat_filled > 0:
        print(f"  (CIR matrix conversion filled {cir_mat_filled} tables)")

    # ── Spinor (double-valued) irrep data ──
    spinor_irreps = data.get("spinor_irreps", [])
    spinor_starts = []
    spinor_counts = []
    for sir in spinor_irreps:
        spinor_starts.append(len(chars_flat))
        spinor_counts.append(len(sir["characters"]))
        chars_flat.extend(sir["characters"])
        # No matrices for spinor irreps
        mat_starts.append(len(matrices_flat))
        mat_counts.append(0)
    if spinor_irreps:
        print(f"  Added {len(spinor_irreps)} spinor irreps to database")

    lines = []
    lines.append("// Auto-generated from iso_data files by scripts/generate_irrep_data.py")
    lines.append("// DO NOT EDIT MANUALLY")
    lines.append("")
    lines.append("use crate::irrep::types::*;")
    lines.append("")
    lines.append("/// Flat array of all character values, indexed by IrrepRecord._char_start.")
    lines.append(f"pub static CHARACTERS: [f64; {len(chars_flat)}] = [")
    # Write in chunks of 10 values per line for readability.
    # Always produce valid f64 literals (integer values need ".0" suffix).
    def _fmt_char(v):
        if abs(v) < 1e-15:
            return "0.0"
        # If it's effectively an integer, use "N.0" format
        if abs(v - round(v)) < 1e-12:
            return f"{int(round(v))}.0"
        # Otherwise keep significant digits, ensuring "0." prefix for fractions < 1
        s = f"{v:.10}".rstrip("0")
        if s == "-0." or s == "0.":
            s = s + "0"
        return s

    for chunk_start in range(0, len(chars_flat), 10):
        chunk = chars_flat[chunk_start:chunk_start + 10]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── MATRICES flat array ──
    lines.append("/// Flat array of all irrep matrix elements, indexed by IrrepRecord._mat_start.")
    lines.append(f"pub static MATRICES: [f64; {len(matrices_flat)}] = [")
    for chunk_start in range(0, len(matrices_flat), 10):
        chunk = matrices_flat[chunk_start:chunk_start + 10]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── SG index ──
    # Build per-SG entry order: scalar entries first, then spinor, contiguous per SG
    sg_entries = defaultdict(list)
    for i in range(len(ml)):
        sg_entries[sg[i]].append(("scalar", i))
    for i, sir in enumerate(spinor_irreps):
        sg_entries[sir["sg"]].append(("spinor", i))

    lines.append("/// SG# → (start_index, count) into IRREPS")
    total_irreps = len(ml) + len(spinor_irreps)
    lines.append(f"pub static SG_IRREP_INDEX: [(u16, u16); 231] = [")
    lines.append("    (0, 0),  // dummy for index 0")
    irrep_idx = 0
    for s in range(1, 231):
        entries = sg_entries.get(s, [])
        if entries:
            start = irrep_idx
            count = len(entries)
            irrep_idx += count
        else:
            start = 0
            count = 0
        lines.append(f"    ({start}, {count}),  // SG {s}")
    lines.append("];")

    # Rebuild the IrrepRecord generation loop to emit entries in SG order
    # First, pre-compute all scalar IrrepRecord data (same as before)
    # ── Pre-compute isotropy subgroup ranges per irrep ──
    iso_starts = []
    iso_counts = []
    for i in range(len(ml)):
        if i < len(data["iso_irrep_ptr"]):
            s = data["iso_irrep_ptr"][i] - 1
        else:
            s = 0
        if i + 1 < len(data["iso_irrep_ptr"]):
            e = data["iso_irrep_ptr"][i + 1] - 1
        else:
            e = len(data["iso_subgroups"])
        iso_starts.append(s)
        iso_counts.append(max(0, e - s))

    mag_iso_ptr = data.get("mag_iso_ptr", [])
    mag_iso_sg_arr = data.get("mag_iso_sg", [])
    mag_iso_starts = []
    mag_iso_counts = []
    for i in range(len(ml)):
        if i < len(mag_iso_ptr):
            s = mag_iso_ptr[i] - 1
        else:
            s = 0
        if i + 1 < len(mag_iso_ptr):
            e = mag_iso_ptr[i + 1] - 1
        else:
            e = len(mag_iso_sg_arr)
        mag_iso_starts.append(s)
        mag_iso_counts.append(max(0, e - s))

    scalar_records = []  # list of dicts with all the generated fields
    for i in range(len(ml)):
        ml_label = ml[i]
        bc_label = bc[i]
        kov_label = kov[i]
        sg_num = sg[i]
        img_code = img[i]
        lif_val = lif[i]
        iso_s = iso_starts[i]
        iso_c = iso_counts[i]
        mag_iso_s = mag_iso_starts[i]
        mag_iso_c = mag_iso_counts[i]
        if 1 <= img_code <= len(img_labels):
            img_name = img_labels[img_code - 1]
        else:
            img_name = "?"
        if img_name and img_name[0] in IMAGE_DIM:
            dim = IMAGE_DIM[img_name[0]]
        else:
            dim = 1
        latex_bc = label_to_latex(bc_label)
        latex_kov = label_to_latex(kov_label)
        kx, ky, kz, kd = _lookup_kvec(kvec_map, sg_num, ml_label)
        char_s = char_starts[i]
        char_c = char_counts[i]
        mat_s = mat_starts[i]
        mat_c = mat_counts[i]
        scalar_records.append({
            "sg": sg_num, "ml": ml_label, "bc": latex_bc, "kov": latex_kov,
            "dim": dim, "img": img_name, "lifshitz": lif_val == 1,
            "spinor": False, "kx": kx, "ky": ky, "kz": kz, "kd": kd,
            "char_s": char_s, "char_c": char_c,
            "mat_s": mat_s, "mat_c": mat_c,
            "iso_s": iso_s, "iso_c": iso_c,
            "mag_iso_s": mag_iso_s, "mag_iso_c": mag_iso_c,
        })

    # Pre-compute spinor IrrepRecord data
    spinor_records = []
    for idx, sir in enumerate(spinor_irreps):
        latex_bc = label_to_latex(sir["ml_label"])
        spinor_records.append({
            "sg": sir["sg"], "ml": sir["ml_label"], "bc": latex_bc, "kov": "",
            "dim": sir["dim"], "img": "?", "lifshitz": False,
            "spinor": True,
            "kx": sir["kx"], "ky": sir["ky"], "kz": sir["kz"], "kd": sir["kd"],
            "char_s": spinor_starts[idx], "char_c": spinor_counts[idx],
            "mat_s": 0, "mat_c": 0,
            "iso_s": 0, "iso_c": 0,
            "mag_iso_s": 0, "mag_iso_c": 0,
        })

    # Now emit IrrepRecord entries in SG order
    lines.append("/// All irreducible representations (scalar + spinor), ordered by SG then k-point.")
    lines.append(f"pub static IRREPS: [IrrepRecord; {total_irreps}] = [")
    for s in range(1, 231):
        for entry_type, entry_idx in sg_entries.get(s, []):
            if entry_type == "scalar":
                r = scalar_records[entry_idx]
            else:
                r = spinor_records[entry_idx]
            lines.append(f"    IrrepRecord {{")
            lines.append(f"        sg: {r['sg']},")
            lines.append(f'        ml: "{escape_rust_str(r["ml"])}",')
            lines.append(f'        bc: "{escape_rust_str(r["bc"])}",')
            lines.append(f'        kov: "{escape_rust_str(r["kov"])}",')
            lines.append(f"        dim: {r['dim']},")
            lines.append(f'        image: "{r["img"]}",')
            lines.append(f"        lifshitz: {str(r['lifshitz']).lower()},")
            lines.append(f"        spinor: {str(r['spinor']).lower()},")
            lines.append(f"        kx: {r['kx']},")
            lines.append(f"        ky: {r['ky']},")
            lines.append(f"        kz: {r['kz']},")
            lines.append(f"        kd: {r['kd']},")
            lines.append(f"        _char_start: {r['char_s']},")
            lines.append(f"        _char_count: {r['char_c']},")
            lines.append(f"        _mat_start: {r['mat_s']},")
            lines.append(f"        _mat_count: {r['mat_c']},")
            lines.append(f"        _iso_start: {r['iso_s']},")
            lines.append(f"        _iso_count: {r['iso_c']},")
            lines.append(f"        _mag_iso_start: {r['mag_iso_s']},")
            lines.append(f"        _mag_iso_count: {r['mag_iso_c']},")
            lines.append(f"    }},")
    lines.append("];")
    lines.append("")

    # ── Isotropy subgroup records (flat, NOT deduplicated) ──
    iso_irrep = data["iso_irrep"]
    iso_irrep_ptr = data["iso_irrep_ptr"]
    iso_sg = data["iso_subgroups"]
    iso_dir = data["iso_direction"]
    iso_dom = data["iso_domains"]
    iso_arms = data["iso_arms"]

    n_irreps = len(ml)
    total_iso = len(iso_sg)

    lines.append(f"/// Isotropy subgroup records (flat, per-irrep ordering).")
    lines.append(f"pub static ISOTROPY_SUBGROUPS: [IsotropyRecord; {total_iso}] = [")
    for i in range(total_iso):
        sg_val = iso_sg[i] if i < len(iso_sg) else 0
        dir_val = iso_dir[i] if i < len(iso_dir) else 0
        dom_val = iso_dom[i] if i < len(iso_dom) else 1
        arms_val = iso_arms[i] if i < len(iso_arms) else 1

        dir_str = dir_map.get(dir_val, f"dir{dir_val}")
        symbol, sch = get_sg_symbol(sg_val)

        lines.append(f"    IsotropyRecord {{")
        lines.append(f"        sg: {sg_val},")
        lines.append(f'        symbol: "{symbol}",')
        lines.append(f'        schoenflies: "{sch}",')
        lines.append(f'        direction: "{dir_str}",')
        lines.append(f"        domains: {dom_val},")
        lines.append(f"        arms: {arms_val},")
        lines.append(f"    }},")
    lines.append("];")
    lines.append("")

    # ── Magnetic isotropy subgroup records ──
    mag_iso_sg_arr = data.get("mag_iso_sg", [])
    mag_nlabel = data.get("mag_nlabel", [])
    mag_bns_label = data.get("mag_bns_label", [])
    mag_dir_by_entry = data.get("mag_dir_by_entry", {})
    total_mag_iso = len(mag_iso_sg_arr)

    lines.append("/// Magnetic isotropy subgroup records (flat, per-irrep ordering).")
    lines.append(f"pub static MAGNETIC_ISOTROPY_SUBGROUPS: [MagneticIsotropyRecord; {total_mag_iso}] = [")
    for i in range(total_mag_iso):
        msg = mag_iso_sg_arr[i]
        # Magnetic SG label lookup
        iso_label = mag_nlabel[msg - 1] if 1 <= msg <= len(mag_nlabel) else f"{msg}"
        bns = mag_bns_label[msg - 1] if 1 <= msg <= len(mag_bns_label) else f"MSG{msg}"
        direction = mag_dir_by_entry.get(i, "(a)")

        lines.append(f"    MagneticIsotropyRecord {{")
        lines.append(f"        mag_sg: {msg},")
        lines.append(f'        bns_label: "{escape_rust_str(bns)}",')
        lines.append(f'        iso_label: "{escape_rust_str(iso_label)}",')
        lines.append(f'        direction: "{escape_rust_str(direction)}",')
        lines.append(f"    }},")
    lines.append("];")
    lines.append("")

    return "\n".join(lines)

# ── main ─────────────────────────────────────────────────────────────────────

def main():
    data = parse_all()

    # Print summary
    n_irreps = data["n_irreps"]
    sg_nums = sorted(set(data["sg_numbers"]))
    print(f"\nSummary: {n_irreps} irreps across {len(sg_nums)} space groups")
    print(f"SG range: {min(sg_nums)}-{max(sg_nums)}")

    # Count per crystal system
    for sys_name, sg_range in CRYSTAL_SYSTEMS.items():
        count = sum(1 for s in data["sg_numbers"] if s in sg_range)
        print(f"  {sys_name}: {count} irreps")

    # Generate Rust code
    print("\nGenerating Rust code...")
    rust_code = generate_rust_data(data)

    out_path = os.path.join(OUT_DIR, "generated_data.rs")
    with open(out_path, "w") as f:
        f.write(rust_code)
    print(f"  Written: {out_path} ({len(rust_code)} bytes)")

    print("\nDone!")

if __name__ == "__main__":
    main()
