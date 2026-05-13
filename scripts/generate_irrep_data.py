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

import re, sys, os, zipfile, io, math
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


# k-point prefixes for special (high-symmetry) points
_SPECIAL_KPTS = {
    "GM", "G", "X", "Y", "Z", "M", "R", "A", "H", "K", "L",
    "S", "T", "U", "V", "W", "N", "P", "F", "C", "D", "E",
    "Q", "B", "J",
}


def _is_special_kpoint(label):
    """Check if a k-point label is a special (high-symmetry) point.

    Special k-points have kspecial=true in the ISOTROPY data and do NOT
    have irtranslation lines after the operator matrix.
    Non-special k-points (lines DT, LD, SM, GP, etc.) have irtranslation.
    """
    # Extract the k-point prefix (letters before the first digit)
    m = re.match(r'^([A-Za-z]+)', label)
    if not m:
        return False
    prefix = m.group(1)
    return prefix in _SPECIAL_KPTS


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
    dim_map = {}     # (SG#, ML_label) -> dim (from PIR header)
    rots_map = {}     # (SG#, ML_label) -> [[r00..r22], ...] per operation
    trans_map = {}   # (SG#, ML_label) -> [[t0,t1,t2], ...] per operation
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
        rots = []         # rotation matrices: list of [r00..r22], 9 ints per op
        trans = []        # translations: list of [t0,t1,t2], 3 f64 per op
        all_matrices = []  # flat: op0_row0, op0_row1, ..., op1_row0, ...
        for _op in range(opcount):
            # Advance to operator matrix line
            i += 1
            if i >= len(lines):
                break

            # Parse operator matrix (16 ints): extract rotation matrix
            op_parts = lines[i].strip().split()
            if len(op_parts) >= 16:
                try:
                    op_nums = [int(x) for x in op_parts[:16]]
                    denom = op_nums[15] if op_nums[15] != 0 else 1
                    r00 = op_nums[0] // denom; r01 = op_nums[1] // denom; r02 = op_nums[2] // denom
                    r10 = op_nums[4] // denom; r11 = op_nums[5] // denom; r12 = op_nums[6] // denom
                    r20 = op_nums[8] // denom; r21 = op_nums[9] // denom; r22 = op_nums[10] // denom
                    rots.append([r00, r01, r02, r10, r11, r12, r20, r21, r22])
                    # operator matrix: 16 ints [r00,r01,r02,t0, r10,r11,r12,t1, r20,r21,r22,t2, 0,0,0,denom]
                    t0 = float(op_nums[3]) / float(denom)
                    t1 = float(op_nums[7]) / float(denom)
                    t2 = float(op_nums[11]) / float(denom)
                    trans.append([t0, t1, t2])
                except ValueError:
                    rots.append([0,0,0, 0,0,0, 0,0,0])
                    trans.append([0.0, 0.0, 0.0])
            else:
                rots.append([0,0,0, 0,0,0, 0,0,0])

            # Move to the irrep/irtranslation line
            i += 1
            if i >= len(lines):
                break

            # Check if this line is irtranslation.
            # Irtranslation is only present for NON-special k-points
            # (lines DT/LD/SM, planes, and general points GP).
            # Special k-points (GM, X, Z, M, R, A, H, K, L, etc.)
            # have kspecial=true → NO irtranslation.
            # We determine kspecial from the k-point label.
            nxt = lines[i].strip().split()
            if (not _is_special_kpoint(label) and nxt
                    and len(nxt) == 4
                    and not any('.' in x for x in nxt)
                    and all(is_int_4(x) for x in nxt)):
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
            rots_map[key] = rots
            trans_map[key] = trans
            dim_map[key] = dim

        # Advance to next irrep header
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if nxt and nxt.split()[0].isdigit() and '"' in nxt:
                break
            i += 1

    return chars_map, matrices_map, rots_map, dim_map, trans_map


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
        irtype = int(after_line[1]) if len(after_line) > 1 else 1  # 1=real, 2=complex
        pmkcount = int(after_line[3])
        opcount = int(after_line[4])
        has_conjugate = (irtype == 2)  # complex CIR irreps have conjugate data after ordinary

        # Skip k-vector data
        kvec_remaining = pmkcount * 16
        while kvec_remaining > 0 and i < len(lines) - 1:
            i += 1
            if i >= len(lines):
                break
            kvec_remaining -= len(lines[i].strip().split())

        # Read operator matrices + complex irrep matrices
        chars = []
        rots = []         # rotation matrices: list of [r00,r01,r02,r10,r11,r12,r20,r21,r22]
        all_matrices = []  # flattened complex matrix elements for all ops
        store_matrices = (needed_labels is None) or ((sg, label) in needed_labels)

        for _op in range(opcount):
            # Advance to operator augmented 4×4 matrix (16 ints).
            # For complex irreps, there may be conjugate complex values
            # before the next operator matrix — skip those.
            i += 1
            while i < len(lines):
                op_str = lines[i].strip()
                if not op_str:
                    i += 1
                    continue
                # Operator matrix line starts with a digit or minus sign
                if op_str[0].isdigit() or op_str[0] == '-':
                    op_parts = op_str.split()
                    if len(op_parts) >= 16:
                        try:
                            op_nums = [int(x) for x in op_parts[:16]]
                            denom = op_nums[15] if op_nums[15] != 0 else 1
                            r00 = op_nums[0] // denom; r01 = op_nums[1] // denom; r02 = op_nums[2] // denom
                            r10 = op_nums[4] // denom; r11 = op_nums[5] // denom; r12 = op_nums[6] // denom
                            r20 = op_nums[8] // denom; r21 = op_nums[9] // denom; r22 = op_nums[10] // denom
                            rots.append([r00, r01, r02, r10, r11, r12, r20, r21, r22])
                        except ValueError:
                            pass
                    break
                # Complex value line (conjugate from previous op): skip
                i += 1
            if i >= len(lines):
                break

            # Move to complex irrep matrix line (ordinary values)
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

            # Skip conjugate complex matrix if present (complex irreps only).
            # Conjugate has the same number of values (dim²) as ordinary.
            if has_conjugate:
                conj_needed = dim * dim
                conj_read = 0
                while conj_read < conj_needed and i < len(lines):
                    for token in lines[i].strip().split():
                        conj_read += 1
                        if conj_read >= conj_needed:
                            break
                    if conj_read >= conj_needed:
                        break
                    i += 1

        key = (sg, label)
        if key not in cir_chars:
            cir_chars[key] = {
                'dim': dim,
                'opcount': opcount,
                'chars': chars,  # list of (re, im, rounded_re)
                'rots': rots,   # list of [r00..r22], 9 ints per op
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
    img_dims = parse_ints(img_lines, img_sec, "image_dimension")
    # image code 1 → img_labels[0] / img_dims[0]
    print(f"  {len(img_labels)} image labels, {len(img_dims)} image dimensions")

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
    chars_map, matrices_map, rots_map, pir_dim_map, pir_trans_map = _parse_pir_characters()
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
        spinor_irreps, spinor_ops = parse_all_spinor()
        print(f"  Parsed {len(spinor_irreps)} spinor irreps, {len(spinor_ops)} SG spin op tables")
    except (ImportError, FileNotFoundError) as e:
        print(f"  Skipping spinor data: {e}")
        spinor_irreps = []
        spinor_ops = {}

    return {
        "n_irreps": n_irreps + len(spinor_irreps),
        "spinor_irreps": spinor_irreps,
        "spinor_ops": spinor_ops,
        "ml_labels": ml_labels,
        "bc_labels": bc_labels,
        "kov_labels": kov_labels,
        "zak_labels": zak_labels,
        "sg_numbers": sg_numbers,
        "images": images,
        "lifshitz": lifshitz,
        "img_labels": img_labels,
        "img_dims": img_dims,
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
        "rots_map": rots_map,
        "pir_trans_map": pir_trans_map,
        "pir_dim_map": pir_dim_map,
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


def _load_hall_operations():
    """Load spglib Hall symmetry operations exported by the Rust test."""
    import json
    hall_path = os.path.join(SCRIPT_DIR, "hall_operations.json")
    if not os.path.exists(hall_path):
        print("  Note: hall_operations.json not found — skipping spglib reorder")
        return {}
    with open(hall_path) as f:
        raw = json.load(f)
    # raw: {"517": {"sg": 221, "rots": [[...], ...], "trans": [[...], ...]}, ...}
    # Build sg → [(hall_number, rots_list, trans_list)] index
    sg_halls = defaultdict(list)
    for hall_str, hd in raw.items():
        sg_halls[hd["sg"]].append((int(hall_str), hd["rots"], hd.get("trans", [])))
    return sg_halls


def _reorder_to_spglib_order(
        sg, ml, chars_flat, char_starts, char_counts,
        matrices_flat, mat_starts, mat_counts,
        pir_rots_flat, pir_rot_starts, rots_map,
        spinor_irreps, spinor_starts, spinor_counts,
        cir_comp_flat=None, cir_comp_rots=None,
        cir_comp_starts=None, cir_comp_counts=None, cir_comp_ops=None,
        kvec_map=None):
    """Reorder per-irrep data from ISOTROPY order into spglib Hall order.

    For compound irreps without PIR rotations, falls back to CIR component
    rotations (CIR_ROTS) to build the mapping.

    Returns per-irrep list: None if unmapped, otherwise list[h_idx→pir_idx].
    Spinor entries appended after scalar entries.
    """
    sg_halls = _load_hall_operations()
    if not sg_halls:
        return [None] * (len(ml) + len(spinor_irreps))

    n_scalar = len(ml)
    reorder_results = []
    sg_hall_choice = {}  # sg → mapping list

    mapped_count = 0
    unmapped_count = 0

    for i in range(n_scalar):
        sg_num = sg[i]
        n_ops = char_counts[i]
        pir_rots = rots_map.get((sg_num, ml[i]), [])

        if n_ops == 0 or not pir_rots:
            reorder_results.append(None)
            unmapped_count += 1
            continue

        hall_candidates = sg_halls.get(sg_num, [])

        best_mapping = None
        for hall_num, hall_rots, hall_trans in hall_candidates:
            if len(hall_rots) == 0:
                continue
            mapping = []
            rot_cache = {}  # rotation_tuple → ISOTROPY index
            for h_rot in hall_rots:
                rot_key = tuple(h_rot)
                if rot_key in rot_cache:
                    mapping.append(rot_cache[rot_key])
                    continue
                found = None
                for p_idx, p_rot in enumerate(pir_rots):
                    if len(p_rot) != 9 or len(h_rot) != 9:
                        continue
                    if all(h_rot[d] == p_rot[d] for d in range(9)):
                        found = p_idx
                        break
                if found is not None:
                    rot_cache[rot_key] = found
                mapping.append(found)
            if sg_num == 9 and ml[i] == 'M1M2':
                print(f"  DEBUG SG9 M1M2 mapping: {mapping}")
                print(f"    hall_rots: {[r[:6] for r in hall_rots]}")
                print(f"    pir_rots (first 4): {[r[:6] for r in pir_rots[:4]]}")
            mapped_op_count = sum(1 for m in mapping if m is not None)
            if mapped_op_count == len(hall_rots):
                best_mapping = mapping
                break

        if best_mapping:
            if sg_num == 9 and ml[i] == 'M1M2':
                cs = char_starts[i]
                print(f"  DEBUG before reorder: n_ops={n_ops} char_start={cs}")
                for op in range(min(n_ops, 8)):
                    print(f"    ISO[{op}] = {chars_flat[cs + op]:.4f}")
            _apply_reorder(chars_flat, char_starts[i], n_ops, best_mapping, 1)
            if sg_num == 9 and ml[i] == 'M1M2':
                cs = char_starts[i]; cc = len(best_mapping)
                print(f"  DEBUG after reorder: char_start={cs} count={cc}")
                for op in range(cc):
                    print(f"    PIR[{op}] = {chars_flat[cs + op]:.4f} (from ISO[{best_mapping[op]}])")
            dim_sq = mat_counts[i] // n_ops if n_ops else 1
            if dim_sq > 0 and mat_counts[i] > 0:
                _apply_reorder(matrices_flat, mat_starts[i], n_ops, best_mapping, dim_sq)
            if n_ops > 0:
                _apply_reorder(pir_rots_flat, pir_rot_starts[i], n_ops, best_mapping, 9)
            char_counts[i] = len(best_mapping)
            sg_hall_choice[sg_num] = (hall_num, best_mapping, hall_trans)
            reorder_results.append(best_mapping)
            mapped_count += 1
        else:
            reorder_results.append(None)
            unmapped_count += 1

    # ── Reorder CIR component data for compound irreps ──
    cir_reordered = 0
    if cir_comp_flat is not None and cir_comp_rots is not None and cir_comp_starts is not None:
        for i in range(n_scalar):
            n_comp = cir_comp_counts[i] if i < len(cir_comp_counts) else 0
            if n_comp == 0:
                continue
            cir_ops = cir_comp_ops[i] if i < len(cir_comp_ops) else 0
            if cir_ops == 0:
                continue
            cir_start = cir_comp_starts[i] if i < len(cir_comp_starts) else 0
            if cir_start == 0:
                continue

            mapping = reorder_results[i]
            if mapping is None:
                continue

            # Reorder CIR data for ALL components using PIR mapping.
            # No Bloch phase needed for the first cir_ops positions:
            # both PIR and CIR share the same ISOTROPY translations,
            # so simple permutation preserves PIR=CIR_sum consistency.
            for comp in range(n_comp):
                comp_char_start = cir_start + comp * cir_ops * 2
                _apply_reorder(cir_comp_flat, comp_char_start, cir_ops, mapping, 2)
                comp_rot_start = (cir_start // 2) * 9 + comp * cir_ops * 9
                _apply_reorder(cir_comp_rots, comp_rot_start, cir_ops, mapping, 9)
            cir_reordered += 1

    if cir_reordered > 0:
        print(f"  CIR reorder: {cir_reordered} compound irreps")

    # Handle spinor irreps: inherit reordering from scalar at same k-point
    for spin_idx, sir in enumerate(spinor_irreps):
        sg_num = sir["sg"]
        n_ops = len(sir.get("op_indices", []))
        # Try to find scalar mapping for this SG
        choice = sg_hall_choice.get(sg_num)
        if choice and n_ops > 0 and spin_idx < len(spinor_starts):
            _, mapping, _ = choice
            # Spinor characters reorder same as scalar at this SG
            s_start = spinor_starts[spin_idx]
            s_count = spinor_counts[spin_idx]
            if s_count > 0:
                _apply_reorder(chars_flat, s_start, s_count, mapping, 1)
            reorder_results.append(mapping)
            mapped_count += 1
        else:
            reorder_results.append(None)
            unmapped_count += 1

    print(f"  Spglib reorder: {mapped_count} mapped, {unmapped_count} unmapped "
          f"(of {len(reorder_results)} irreps, {len(sg_halls)} SGs)")
    return reorder_results, sg_hall_choice


def _build_padding_plans(sg, ml, cir_comp_starts, cir_comp_counts, cir_comp_ops,
                          cir_comp_rots, reorder_results):
    """Build padding plans for unmapped compound irreps using CIR_ROTS matching.

    For unmapped entries (reorder_results[i] is None) that have CIR data,
    match the first component's CIR rotation matrices against Hall operations
    to build a cir_to_hall mapping. Each plan is (irrep_idx, hall_ops, cir_to_hall).

    Returns list of padding plans.
    """
    sg_halls = _load_hall_operations()
    if not sg_halls:
        return []

    padding_plans = []  # [(irrep_idx, hall_ops, cir_to_hall), ...]
    n_scalar = len(ml)

    for i in range(n_scalar):
        if reorder_results[i] is not None:
            continue  # Already mapped via PIR_ROTS
        if cir_comp_counts[i] == 0:
            continue  # No CIR data at all

        n_ops = cir_comp_ops[i]  # CIR ops count (little group size)
        if n_ops == 0:
            continue

        # Get first component's CIR rotation matrices
        rot_start = (cir_comp_starts[i] // 2) * 9
        cir_rots = []
        for op_idx in range(n_ops):
            r9 = cir_comp_rots[rot_start + op_idx * 9:rot_start + (op_idx + 1) * 9]
            cir_rots.append(r9)

        sg_num = sg[i]
        hall_candidates = sg_halls.get(sg_num, [])

        best_plan = None
        for hall_num, hall_rots, hall_trans in hall_candidates:
            hall_ops = len(hall_rots)
            if hall_ops == 0 or n_ops > hall_ops:
                continue

            # Match each CIR op → Hall op
            cir_to_hall = []  # cir_to_hall[ci] = hi
            for ci, c_rot in enumerate(cir_rots):
                found = None
                for hi, h_rot in enumerate(hall_rots):
                    if all(c_rot[d] == h_rot[d] for d in range(9)):
                        found = hi
                        break
                cir_to_hall.append(found)

            if None not in cir_to_hall:
                # Verify identity maps to identity
                if cir_to_hall[0] != 0:
                    print(f"  WARNING padding: SG{sg_num} {ml[i]}: "
                          f"CIR identity → Hall[{cir_to_hall[0]}], not 0")
                best_plan = (i, hall_ops, cir_to_hall)
                break  # Use first matching Hall setting

        if best_plan:
            padding_plans.append(best_plan)
        else:
            print(f"  WARNING padding: SG{sg_num} {ml[i]}: "
                  f"no Hall setting matches all CIR rots (n_ops={n_ops})")

    return padding_plans


def _apply_reorder(arr, start, count, mapping, stride):
    """Reorder `count` items in `arr` starting at `start`, each item of size `stride`.

    mapping[h_idx] = orig_idx.  Only reorders the first `len(mapping)` items.
    Extra items beyond len(mapping) are left untouched.
    """
    if stride == 0 or count == 0:
        return
    n_reorder = min(count, len(mapping))
    if n_reorder == 0:
        return
    old = arr[start:start + count * stride]
    for new_pos in range(n_reorder):
        old_pos = mapping[new_pos]
        if old_pos is not None and old_pos < count:
            src = start + old_pos * stride
            dst = start + new_pos * stride
            offset = src - start
            for d in range(stride):
                arr[dst + d] = old[offset + d]


def _apply_padding_plans(padding_plans, chars_flat, char_starts, char_counts,
                          matrices_flat, mat_starts, mat_counts,
                          pir_rots_flat, pir_rot_starts,
                          cir_comp_flat, cir_comp_rots, cir_comp_starts, cir_comp_counts, cir_comp_ops,
                          spinor_starts, spinor_counts):
    """Rebuild flat arrays with padded entries for compound irreps expanded to Hall size."""
    n_scalar = len(char_starts)
    plans_by_idx = {i: (hall_ops, cir_to_hall) for i, hall_ops, cir_to_hall in padding_plans}

    # Rebuild chars_flat
    new_chars = []
    new_char_starts = []
    new_char_counts = []
    for i in range(n_scalar):
        new_char_starts.append(len(new_chars))
        if i in plans_by_idx:
            hall_ops, cir_to_hall = plans_by_idx[i]
            old = chars_flat[char_starts[i]:char_starts[i] + char_counts[i]]
            for h in range(hall_ops):
                ci = None
                for cii, hi in enumerate(cir_to_hall):
                    if hi == h:
                        ci = cii
                        break
                if ci is not None and ci < len(old):
                    new_chars.append(old[ci])
                else:
                    new_chars.append(0.0)
            new_char_counts.append(hall_ops)
        else:
            n = char_counts[i]
            new_chars.extend(chars_flat[char_starts[i]:char_starts[i] + n])
            new_char_counts.append(n)

    # Rebuild matrices_flat
    new_mats = []
    new_mat_starts = []
    new_mat_counts = []
    for i in range(n_scalar):
        new_mat_starts.append(len(new_mats))
        if i in plans_by_idx:
            hall_ops, cir_to_hall = plans_by_idx[i]
            old_m = matrices_flat[mat_starts[i]:mat_starts[i] + mat_counts[i]]
            cir_ops = len(cir_to_hall)
            dim_sq = mat_counts[i] // cir_ops if cir_ops else 0
            for h in range(hall_ops):
                ci = None
                for cii, hi in enumerate(cir_to_hall):
                    if hi == h:
                        ci = cii
                        break
                if ci is not None and ci < cir_ops and dim_sq > 0:
                    new_mats.extend(old_m[ci * dim_sq:(ci + 1) * dim_sq])
                else:
                    new_mats.extend([0.0] * max(dim_sq, 0))
            new_mat_counts.append(hall_ops * max(dim_sq, 1))
        else:
            n = mat_counts[i]
            new_mats.extend(matrices_flat[mat_starts[i]:mat_starts[i] + n])
            new_mat_counts.append(n)

    # Rebuild pir_rots_flat
    new_rots = []
    new_rot_starts = []
    n_rot_entries = len(pir_rot_starts)
    for i in range(n_rot_entries):
        new_rot_starts.append(len(new_rots))
        orig_end = pir_rot_starts[i + 1] if i + 1 < n_rot_entries else len(pir_rots_flat)
        if i < n_scalar and i in plans_by_idx:
            hall_ops, cir_to_hall = plans_by_idx[i]
            cir_ops = len(cir_to_hall)
            old_r = pir_rots_flat[pir_rot_starts[i]:pir_rot_starts[i] + cir_ops * 9]
            for h in range(hall_ops):
                ci = None
                for cii, hi in enumerate(cir_to_hall):
                    if hi == h:
                        ci = cii
                        break
                if ci is not None and ci < cir_ops and len(old_r) > ci * 9:
                    new_rots.extend(old_r[ci * 9:(ci + 1) * 9])
                else:
                    new_rots.extend([0] * 9)
        else:
            new_rots.extend(pir_rots_flat[pir_rot_starts[i]:orig_end])

    # Rebuild cir_comp_flat and cir_comp_rots
    new_cir_flat = []
    new_cir_rots = []
    new_cir_starts = []
    for i in range(len(cir_comp_starts)):
        n_comp = cir_comp_counts[i]
        old_ops = cir_comp_ops[i]
        if n_comp == 0 or old_ops == 0:
            new_cir_starts.append(0)
            continue
        new_cir_starts.append(len(new_cir_flat))
        if i in plans_by_idx:
            hall_ops, cir_to_hall = plans_by_idx[i]
            for comp in range(n_comp):
                old_start = cir_comp_starts[i] + comp * old_ops * 2
                old_rot_start = (cir_comp_starts[i] // 2) * 9 + comp * old_ops * 9
                for h in range(hall_ops):
                    ci = None
                    for cii, hi in enumerate(cir_to_hall):
                        if hi == h:
                            ci = cii
                            break
                    if ci is not None and ci < old_ops:
                        new_cir_flat.append(cir_comp_flat[old_start + ci * 2])
                        new_cir_flat.append(cir_comp_flat[old_start + ci * 2 + 1])
                        new_cir_rots.extend(cir_comp_rots[old_rot_start + ci * 9:old_rot_start + (ci + 1) * 9])
                    else:
                        new_cir_flat.append(0.0)
                        new_cir_flat.append(0.0)
                        new_cir_rots.extend([0] * 9)
            cir_comp_ops[i] = hall_ops
        else:
            old_start = cir_comp_starts[i]
            total_chars = n_comp * old_ops * 2
            new_cir_flat.extend(cir_comp_flat[old_start:old_start + total_chars])
            old_rot_start = (cir_comp_starts[i] // 2) * 9
            total_rots = n_comp * old_ops * 9
            new_cir_rots.extend(cir_comp_rots[old_rot_start:old_rot_start + total_rots])

    # Copy back, preserving spinor data at the end.
    # Use spinor_starts[0] as the true scalar/spinor boundary, because
    # char_counts may have been truncated by _reorder_to_spglib_order
    if spinor_starts:
        old_scalar_chars_len = spinor_starts[0]
    else:
        old_scalar_chars_len = sum(char_counts)
    spinor_chars_tail = chars_flat[old_scalar_chars_len:] if old_scalar_chars_len < len(chars_flat) else []
    new_scalar_chars_len = len(new_chars)
    for j in range(len(spinor_starts)):
        old_spinor_start = spinor_starts[j]
        if old_spinor_start >= old_scalar_chars_len:
            spinor_starts[j] = new_scalar_chars_len + (old_spinor_start - old_scalar_chars_len)
    chars_flat[:] = new_chars + spinor_chars_tail
    char_starts[:] = new_char_starts
    char_counts[:] = new_char_counts

    old_scalar_mats_len = sum(mat_counts)
    spinor_mats_tail = matrices_flat[old_scalar_mats_len:] if old_scalar_mats_len < len(matrices_flat) else []
    new_scalar_mats_len = sum(new_mat_counts)
    for j in range(len(spinor_starts)):  # spinor mats at end
        pass  # spinor mat_starts are appended separately, not in mat_starts list
    matrices_flat[:] = new_mats + spinor_mats_tail
    mat_starts[:] = new_mat_starts
    mat_counts[:] = new_mat_counts

    pir_rots_flat[:] = new_rots
    pir_rot_starts[:] = new_rot_starts
    cir_comp_flat[:] = new_cir_flat
    cir_comp_rots[:] = new_cir_rots
    cir_comp_starts[:] = new_cir_starts


def generate_rust_data(data):
    """Generate the content of generated_data.rs."""
    ml  = data["ml_labels"]
    bc  = data["bc_labels"]
    kov = data["kov_labels"]
    sg  = data["sg_numbers"]
    img = data["images"]
    lif = data["lifshitz"]
    img_labels = data["img_labels"]
    img_dims = data.get("img_dims", [])
    pir_dim_map = data.get("pir_dim_map", {})

    # direction labels
    dir_map = data["dir_map"]
    # k-vector map: (SG#, ML_label) -> (kx, ky, kz, denom)
    kvec_map = data["kvec_map"]
    # character map: (SG#, ML_label) -> [char1, char2, ...]
    chars_map = data.get("chars_map", {})
    # matrix map: (SG#, ML_label) -> [flat_values...]
    matrices_map = data.get("matrices_map", {})

    cir_data = data.get("cir_data", {})

    # PIR rotation matrices for H_ops → PIR order mapping (Wigner test)
    rots_map = data.get("rots_map", {})

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

    # ── Build flat PIR_ROTS and PIR_TRANS arrays ──
    pir_trans_map = data.get("pir_trans_map", {})
    pir_rots_flat = []
    pir_trans_flat = []
    pir_rot_starts = []
    pir_trans_starts = []
    for i in range(len(ml)):
        rts = rots_map.get((sg[i], ml[i]), [])
        trs = pir_trans_map.get((sg[i], ml[i]), [])
        pir_rot_starts.append(len(pir_rots_flat))
        pir_trans_starts.append(len(pir_trans_flat))
        for r9 in rts:
            pir_rots_flat.extend(r9)
        for t3 in trs:
            pir_trans_flat.extend(t3)
        expected_ops = char_counts[i]
        needed_rots = expected_ops * 9
        needed_trans = expected_ops * 3
        current_rots = len(pir_rots_flat) - pir_rot_starts[-1]
        current_trans = len(pir_trans_flat) - pir_trans_starts[-1]
        if current_rots < needed_rots:
            pir_rots_flat.extend([0] * (needed_rots - current_rots))
        if current_trans < needed_trans:
            pir_trans_flat.extend([0.0] * (needed_trans - current_trans))
    # (Spinor irreps: PIR rot starts added after spinor_irreps variable is available)
    # For compound labels like Z1Z4 = Z1 ⊕ Z4, store the individual CIR
    # complex character tables as (re, im) pairs.  Used for Wigner test.
    cir_comp_flat = []   # (re, im) pairs, flattened
    cir_comp_rots = []   # rotation matrices (9 ints per op), same order as chars
    cir_comp_starts = []  # per-irrep start index (0 = not compound)
    cir_comp_counts = []  # number of CIR components (0 = not compound)
    cir_comp_ops = []     # operations per CIR component
    cir_comp_total = 0
    cir_comp_rejected = 0
    for i in range(len(ml)):
        parts = _decompose_compound_label(ml[i])
        if parts and len(parts) >= 2:
            comp_chars = []
            comp_rots = []
            n_ops = 0
            for p_idx, part in enumerate(parts):
                pk = (sg[i], part)
                if pk in cir_data:
                    entry = cir_data[pk]
                    n_ops = entry['opcount']
                    # For same-label compounds (e.g. P3P3 = P3 ⊕ P3*),
                    # conjugate the second occurrence so imaginary parts cancel.
                    # PIR irreps are always real; the CIR sum must be real.
                    conjugate_this = (p_idx > 0 and part == parts[0])
                    for (re_val, im_val, _) in entry['chars']:
                        comp_chars.append(re_val)
                        if conjugate_this:
                            comp_chars.append(-im_val)
                        else:
                            comp_chars.append(im_val)
                    # Rotation matrices (if available)
                    rots = entry.get('rots', [])
                    for r9 in rots:
                        comp_rots.extend(r9)
                else:
                    comp_chars = []
                    break

            # ── Full character table validation ──
            # For same-label compounds (e.g. P3P3), try both:
            #   A: straight sum  χ + χ
            #   B: conjugate sum χ + χ*  (imaginary parts cancel)
            # For different-label compounds, only straight sum.
            # Accept only if ALL operations match the PIR character table.
            valid = False
            if comp_chars:
                pir_ch = _lookup_chars(chars_map, sg[i], ml[i], kvec_map)
                if pir_ch and len(pir_ch) > 0 and len(pir_ch) == n_ops:
                    all_same = all(p == parts[0] for p in parts)
                    use_conjugate = False  # which candidate was selected

                    # Candidate 1: straight sum
                    straight_ok = True
                    for op in range(n_ops):
                        s_re = 0.0; s_im = 0.0
                        for p_idx in range(len(parts)):
                            s_re += comp_chars[2 * (p_idx * n_ops + op)]
                            s_im += comp_chars[2 * (p_idx * n_ops + op) + 1]
                        if abs(s_re - pir_ch[op]) > 0.01 or abs(s_im) > 0.01:
                            straight_ok = False
                            break
                    if straight_ok:
                        valid = True

                    # Candidate 2: conjugate second (only for same-label 2-part)
                    if not valid and all_same and len(parts) == 2:
                        conj_ok = True
                        for op in range(n_ops):
                            s_re = comp_chars[2 * op] + comp_chars[2 * (n_ops + op)]
                            s_im = comp_chars[2 * op + 1] - comp_chars[2 * (n_ops + op) + 1]
                            if abs(s_re - pir_ch[op]) > 0.01 or abs(s_im) > 0.01:
                                conj_ok = False
                                break
                        if conj_ok:
                            valid = True
                            use_conjugate = True

                    # Apply conjugation to comp_chars if that candidate won
                    if use_conjugate:
                        new_chars = []
                        for op in range(n_ops):
                            new_chars.append(comp_chars[2 * op])
                            new_chars.append(comp_chars[2 * op + 1])
                        for op in range(n_ops):
                            new_chars.append(comp_chars[2 * (n_ops + op)])
                            new_chars.append(-comp_chars[2 * (n_ops + op) + 1])
                        comp_chars = new_chars

            if comp_chars and valid:
                cir_comp_starts.append(len(cir_comp_flat))
                cir_comp_counts.append(len(parts))
                cir_comp_ops.append(n_ops)
                cir_comp_flat.extend(comp_chars)
                # Always extend rotations to keep arrays in sync.
                total_rots_needed = len(parts) * n_ops * 9
                if len(comp_rots) == total_rots_needed:
                    cir_comp_rots.extend(comp_rots)
                else:
                    cir_comp_rots.extend([0] * total_rots_needed)
                cir_comp_total += 1
            else:
                cir_comp_starts.append(0); cir_comp_counts.append(0); cir_comp_ops.append(0)
                cir_comp_rejected += 1
        else:
            cir_comp_starts.append(0); cir_comp_counts.append(0); cir_comp_ops.append(0)
    print(f"  CIR component chars: {cir_comp_total} compound irreps, {cir_comp_rejected} rejected (id mismatch), {len(cir_comp_flat)} values, {len(cir_comp_rots)} rotation ints")

    # ── Spinor (double-valued) irrep data ──
    spinor_irreps = data.get("spinor_irreps", [])
    # Build (SG, k_label) → rotation data lookup from scalar irreps.
    # All irreps at the same k-point share the same little group operations,
    # so spinor irreps can reuse scalar rotation data.
    kpoint_rots = {}  # (sg, k_label) -> flat [r00..r22, ...] per op
    for i in range(len(ml)):
        rts = rots_map.get((sg[i], ml[i]), [])
        if rts:
            k_label = _kpoint_label_from_ml(ml[i])
            key = (sg[i], k_label)
            if key not in kpoint_rots:
                flat = []
                for r9 in rts:
                    flat.extend(r9)
                kpoint_rots[key] = flat

    for sir in spinor_irreps:
        key = (sir['sg'], sir['k_label'])
        rts = kpoint_rots.get(key, [])
        pir_rot_starts.append(len(pir_rots_flat))
        if rts:
            pir_rots_flat.extend(rts)
        else:
            # Fallback: zero-pad based on character count
            n_ops = len(sir.get('characters', []))
            pir_rots_flat.extend([0] * (n_ops * 9))
    spinor_starts = []
    spinor_counts = []
    spin_extra_flat = []   # extra character values for spinor Wigner test
    spin_extra_starts = [] # per-irrep start index (0 = no extra)
    spin_extra_counts = [] # number of extra values (0 = no extra)
    for sir in spinor_irreps:
        spinor_starts.append(len(chars_flat))
        all_chars = sir["characters"]
        spinor_counts.append(len(all_chars))
        chars_flat.extend(all_chars)
        # Split extra characters: first n_lg are standard, rest are extra
        n_lg = len(sir.get("op_indices", []))
        extra = all_chars[n_lg:] if len(all_chars) > n_lg else []
        spin_extra_starts.append(len(spin_extra_flat))
        spin_extra_counts.append(len(extra))
        spin_extra_flat.extend(extra)
        # No matrices for spinor irreps
        mat_starts.append(len(matrices_flat))
        mat_counts.append(0)
    if spinor_irreps:
        n_with_extra = sum(1 for c in spin_extra_counts if c > 0)
        print(f"  Added {len(spinor_irreps)} spinor irreps ({n_with_extra} with extra chars)")

    # ── Spinor operation arrays ──
    spinor_ops_data = data.get("spinor_ops", {})
    spin_op_rots = []   # 9 i32 per op
    spin_op_trans = []  # 3 f64 per op
    spin_op_su2 = []    # 4 f64 per op
    spin_op_sg_start = [0] * 231  # 0-indexed, SG 0 unused
    spin_op_sg_count = [0] * 231
    for sg_num in range(1, 231):
        ops = spinor_ops_data.get(sg_num, [])
        spin_op_sg_start[sg_num] = len(spin_op_rots) // 9  # count in ops
        spin_op_sg_count[sg_num] = len(ops)
        for op in ops:
            spin_op_rots.extend(op['rot'])
            spin_op_trans.extend(op['trans'])
            spin_op_su2.extend(op['su2'])
    total_spin_ops = len(spin_op_rots) // 9

    # ── Spinor little-group character counts ──
    spin_lg_counts = []
    spin_lg_op_indices_flat = []  # op_indices flattened
    spin_lg_op_starts = []        # per-irrep start index
    spin_lg_op_counts = []        # per-irrep count (= len(op_indices))
    for sir in spinor_irreps:
        ops = sir.get("op_indices", [])
        spin_lg_counts.append(len(ops))
        spin_lg_op_starts.append(len(spin_lg_op_indices_flat))
        spin_lg_op_counts.append(len(ops))
        spin_lg_op_indices_flat.extend(ops)

    # ── Reorder characters/matrices/rots from ISOTROPY order to spglib order ──
    reorder_map_per_irrep, sg_hall_choice = _reorder_to_spglib_order(
        sg, ml, chars_flat, char_starts, char_counts,
        matrices_flat, mat_starts, mat_counts,
        pir_rots_flat, pir_rot_starts, rots_map,
        spinor_irreps, spinor_starts, spinor_counts,
        cir_comp_flat=cir_comp_flat, cir_comp_rots=cir_comp_rots,
        cir_comp_starts=cir_comp_starts, cir_comp_counts=cir_comp_counts,
        cir_comp_ops=cir_comp_ops,
        kvec_map=kvec_map)
    # CIR component data is also reordered in-place.
    # CIR components are used by the Wigner CIR path which matches rotations
    # at runtime via build_h_to_cir_map().  The characters stay in ISOTROPY order.
    # reorder_map_per_irrep[i] = None (unmapped) or list[h_idx→pir_idx] (mapped)
    # For spinor irreps: entries past len(ml) are the spinor reorder maps

    # ── Phase C: CIR padding for unmapped compound irreps ──
    padding_plans = _build_padding_plans(
        sg, ml, cir_comp_starts, cir_comp_counts, cir_comp_ops, cir_comp_rots,
        reorder_map_per_irrep)
    # CIR data is reordered in-place by _reorder_to_spglib_order above.
    # No expansion: CIR characters depend only on rotation, and CIR covers
    # only distinct rotation types (fewer ops than full little group).
    # PIR = Σ CIR_i * exp(i*k·t) accounts for translation differences;
    # CIR stays at the little co-group size.

    if padding_plans:
        print(f"  CIR padding: {len(padding_plans)} entries expanded to Hall size")
        _apply_padding_plans(padding_plans, chars_flat, char_starts, char_counts,
                             matrices_flat, mat_starts, mat_counts,
                             pir_rots_flat, pir_rot_starts,
                             cir_comp_flat, cir_comp_rots, cir_comp_starts,
                             cir_comp_counts, cir_comp_ops,
                             spinor_starts, spinor_counts)

    # ── Verify identity characters for ALL scalar entries ──
    bad_chars = []
    for i in range(len(ml)):
        cs = char_starts[i]
        cc = char_counts[i]
        if cc > 0 and cs < len(chars_flat):
            chi0 = chars_flat[cs]
            if chi0 <= 0.0:
                bad_chars.append(f"  SG{sg[i]} {ml[i]}: χ(E)={chi0} (n_ops={cc})")
    if bad_chars:
        print(f"  ⚠ {len(bad_chars)} ENTRIES WITH NON-POSITIVE χ(E) AFTER PADDING:")
        for line in bad_chars[:15]:
            print(line)
        if len(bad_chars) > 15:
            print(f"  ... and {len(bad_chars) - 15} more")
    else:
        print(f"  ✓ All {len(ml)} scalar entries have positive χ(E)")

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

    # ── PIR rotation matrices ──
    lines.append("/// Rotation matrices for PIR operations, 9 i32 per op, same order as CHARACTERS.")
    lines.append("/// Used to build H_ops → PIR index mapping for the Wigner test.")
    lines.append(f"pub static PIR_ROTS: [i32; {len(pir_rots_flat)}] = [")
    for chunk_start in range(0, len(pir_rots_flat), 9):
        chunk = pir_rots_flat[chunk_start:chunk_start + 9]
        vals = ", ".join(str(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── PIR translation vectors ──
    lines.append("/// Translation vectors for PIR operations, 3 f64 per op, same order as CHARACTERS and PIR_ROTS.")
    lines.append("/// Used with PIR_ROTS for full Seitz matching (rotation + translation).")
    lines.append(f"pub static PIR_TRANS: [f64; {len(pir_trans_flat)}] = [")
    for chunk_start in range(0, len(pir_trans_flat), 3):
        chunk = pir_trans_flat[chunk_start:chunk_start + 3]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── CIR component complex characters ──
    lines.append("/// Complex character tables for CIR components of compound irreps.")
    lines.append("/// Stored as (re, im) pairs.  Used by Wigner test for correct")
    lines.append("/// Type C classification of magnetic co-representations.")
    lines.append(f"pub static CIR_COMPONENT_CHARS: [f64; {len(cir_comp_flat)}] = [")
    for chunk_start in range(0, len(cir_comp_flat), 10):
        chunk = cir_comp_flat[chunk_start:chunk_start + 10]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── CIR rotation matrices (for runtime order matching) ──
    lines.append("/// Rotation matrices for CIR operations, 9 i32 per op, same order as CIR_COMPONENT_CHARS.")
    lines.append(f"pub static CIR_ROTS: [i32; {len(cir_comp_rots)}] = [")
    for chunk_start in range(0, len(cir_comp_rots), 9):
        chunk = cir_comp_rots[chunk_start:chunk_start + 9]
        vals = ", ".join(str(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── Spinor extra characters (Wigner test contributions) ──
    lines.append("/// Extra character values for spinor irreps at BZ boundary k-points.")
    lines.append("/// Sum over these values gives the Wigner indicator contribution.")
    lines.append("/// Indexed by IrrepRecord._spin_extra_start / _spin_extra_count.")
    lines.append(f"pub static SPIN_EXTRA_CHARS: [f64; {len(spin_extra_flat)}] = [")
    for chunk_start in range(0, len(spin_extra_flat), 10):
        chunk = spin_extra_flat[chunk_start:chunk_start + 10]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── Spinor little-group operation indices ──
    lines.append("/// Per-irrep little-group operation indices into SPIN_OP_* arrays.")
    lines.append("/// Maps global spin op index → local character table position.")
    lines.append("/// Indexed by IrrepRecord._spin_lg_op_start / _spin_lg_op_count.")
    lines.append(f"pub static SPIN_LG_OP_INDICES: [u16; {len(spin_lg_op_indices_flat)}] = [")
    for chunk_start in range(0, len(spin_lg_op_indices_flat), 12):
        chunk = spin_lg_op_indices_flat[chunk_start:chunk_start + 12]
        vals = ", ".join(str(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── Spinor (double-group) operation arrays ──
    lines.append("/// Spinor symmetry operations with SU(2) lifts, indexed by SG number.")
    lines.append("/// Use [`SPIN_OP_SG_INDEX`] to find start and count for each SG.")
    lines.append(f"pub static SPIN_OP_ROTS: [i32; {len(spin_op_rots)}] = [")
    for chunk_start in range(0, len(spin_op_rots), 9):
        chunk = spin_op_rots[chunk_start:chunk_start + 9]
        vals = ", ".join(str(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")
    lines.append(f"pub static SPIN_OP_TRANS: [f64; {len(spin_op_trans)}] = [")
    for chunk_start in range(0, len(spin_op_trans), 3):
        chunk = spin_op_trans[chunk_start:chunk_start + 3]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")
    lines.append(f"/// SU(2) Pauli coefficients [u₀,u₁,u₂,u₃] per spin operation.")
    lines.append(f"/// U = u₀·I + i(u₁·σx + u₂·σy + u₃·σz).")
    lines.append(f"/// Verified 229/229 SGs at 100% closure.")
    lines.append(f"pub static SPIN_OP_SU2: [f64; {len(spin_op_su2)}] = [")
    for chunk_start in range(0, len(spin_op_su2), 4):
        chunk = spin_op_su2[chunk_start:chunk_start + 4]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")
    lines.append(f"/// SG# → (operation_start, operation_count) into SPIN_OP_* arrays.")
    lines.append(f"pub static SPIN_OP_SG_INDEX: [(u16, u8); 231] = [")
    lines.append("    (0, 0),  // dummy for index 0")
    for s in range(1, 231):
        lines.append(f"    ({spin_op_sg_start[s]}, {spin_op_sg_count[s]}),  // SG {s}")
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
    dim_warnings = 0
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
        latex_bc = label_to_latex(bc_label)
        latex_kov = label_to_latex(kov_label)
        kx, ky, kz, kd = _lookup_kvec(kvec_map, sg_num, ml_label)
        char_s = char_starts[i]
        char_c = char_counts[i]

        # Determine dimension from the final character table (which already
        # went through _lookup_chars fallback matching), so dim and the
        # written CHARACTERS are guaranteed consistent.
        # Priority: χ(E) from CHARACTERS > PIR header > image_dimension > error
        pir_d = pir_dim_map.get((sg_num, ml_label))
        final_chars = chars_flat[char_s:char_s + char_c] if char_c > 0 else []
        if final_chars:
            dim = int(round(final_chars[0]))
        elif pir_d is not None and pir_d > 0:
            dim = pir_d
        elif 1 <= img_code <= len(img_dims):
            dim = img_dims[img_code - 1]
        else:
            raise ValueError(
                f"Cannot determine dim for SG{sg_num} {ml_label}, img={img_name}"
            )

        # Cross-check: warn if PIR header dim disagrees with final χ(E)
        if pir_d is not None and final_chars:
            chi_e = int(round(final_chars[0]))
            if pir_d != chi_e:
                print(f"  WARNING dim mismatch: SG{sg_num} {ml_label}: "
                      f"PIR header dim={pir_d}, χ(E)={final_chars[0]}, img={img_name}")
                dim_warnings += 1

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
            "cir_s": cir_comp_starts[i], "cir_c": cir_comp_counts[i], "cir_o": cir_comp_ops[i],
            "pir_rot_s": pir_rot_starts[i],
            "spin_lg_count": 0,
            "spin_lg_op_s": 0,
            "spin_lg_op_c": 0,
            "spin_extra_s": 0,
            "spin_extra_c": 0,
        })

    if dim_warnings > 0:
        print(f"  ⚠ {dim_warnings} irreps have PIR header dim ≠ χ(E) (fixed to χ(E))")

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
            "cir_s": 0, "cir_c": 0, "cir_o": 0,
            "pir_rot_s": pir_rot_starts[len(ml) + idx],
            "spin_lg_count": spin_lg_counts[idx],
            "spin_lg_op_s": spin_lg_op_starts[idx],
            "spin_lg_op_c": spin_lg_op_counts[idx],
            "spin_extra_s": spin_extra_starts[idx],
            "spin_extra_c": spin_extra_counts[idx],
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
            lines.append(f"        _pir_rot_start: {r['pir_rot_s']},")
            lines.append(f"        _char_start: {r['char_s']},")
            lines.append(f"        _char_count: {r['char_c']},")
            lines.append(f"        _mat_start: {r['mat_s']},")
            lines.append(f"        _mat_count: {r['mat_c']},")
            lines.append(f"        _iso_start: {r['iso_s']},")
            lines.append(f"        _iso_count: {r['iso_c']},")
            lines.append(f"        _mag_iso_start: {r['mag_iso_s']},")
            lines.append(f"        _mag_iso_count: {r['mag_iso_c']},")
            lines.append(f"        _cir_start: {r['cir_s']},")
            lines.append(f"        _cir_count: {r['cir_c']},")
            lines.append(f"        _cir_ops: {r['cir_o']},")
            lines.append(f"        _spin_lg_count: {r['spin_lg_count']},")
            lines.append(f"        _spin_lg_op_start: {r['spin_lg_op_s']},")
            lines.append(f"        _spin_lg_op_count: {r['spin_lg_op_c']},")
            lines.append(f"        _spin_extra_start: {r['spin_extra_s']},")
            lines.append(f"        _spin_extra_count: {r['spin_extra_c']},")
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
    lines.append("// SG setting data: basis matrices and origin shifts from ISOTROPY.")
    lines.append("include!(\"settings_data.rs\");")
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
