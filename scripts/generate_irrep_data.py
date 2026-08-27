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

import re, sys, os, zipfile, io, math, hashlib
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(__file__)
ISO_DIR = os.path.join(SCRIPT_DIR, "..", "isotropy_subgroup")
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "src", "irrep")

# These archives are the complete, pinned upstream inputs for generated data.
# Keep this list versioned with the generator: parsing an unverified archive
# would make generated metadata impossible to reproduce or audit.
PINNED_ARCHIVE_SHA256 = {
    "PIR_data.zip": "e909a4f0121688b0590ccaec10b0276171bc24619cf7eb562ba441268c01e121",
    "CIR_data.zip": "f4edcb2852b83a86d1b58f29fb862d9124a227cfc90f9e1ae17d2c97585264e6",
    "iso.zip": "568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b",
}


def _verify_pinned_archives():
    """Verify every required upstream archive before parsing any input."""
    for archive, expected in PINNED_ARCHIVE_SHA256.items():
        path = os.path.join(ISO_DIR, archive)
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"required pinned archive is missing: {path}"
            )
        digest = hashlib.sha256()
        with open(path, "rb") as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(chunk)
        actual = digest.hexdigest()
        if actual != expected:
            raise ValueError(
                f"pinned archive hash mismatch for {path}: "
                f"expected {expected}, got {actual}"
            )
    print("  Verified pinned PIR/CIR/ISO archive SHA-256 hashes")

# Import the direction mapping module from this directory
sys.path.insert(0, os.path.dirname(__file__))
from direction_map import build_direction_map

# ── zip-based file reading ──────────────────────────────────────────────────

def _open_zip_path(zip_name, inner_path):
    """Open a UTF-8 member from one pinned ZIP archive.

    Generated data has one authoritative input path: the verified ZIP archive.
    Extracted directories are deliberately ignored.  The member name must be
    unique when matched by archive-relative suffix.  The member is copied into
    memory before the ``ZipFile`` is closed, so the returned text stream stays
    valid for its caller's ``with`` block.
    """
    zip_path = os.path.join(ISO_DIR, zip_name)
    if not os.path.isfile(zip_path):
        raise FileNotFoundError(f"required ZIP archive not found: {zip_path}")
    with zipfile.ZipFile(zip_path) as zf:
        matches = [
            name for name in zf.namelist()
            if name == inner_path or name.endswith("/" + inner_path)
        ]
        if not matches:
            raise FileNotFoundError(f"{inner_path} not found in {zip_name}")
        if len(matches) > 1:
            raise ValueError(
                f"ambiguous archive member {inner_path!r} in {zip_name}: {matches}"
            )
        contents = zf.read(matches[0])
    return io.TextIOWrapper(io.BytesIO(contents), encoding="utf-8")

def read_file(inner_path, zip_name="iso.zip"):
    """Read lines from a member of a pinned ZIP archive."""
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


def _parse_pir_header(line, line_number):
    """Parse one PIR record header with source-line diagnostics."""
    tokens = line.strip().split()
    quoted = re.findall(r'"([^"]*)"', line)
    if not tokens or not tokens[0].isdigit() or len(quoted) < 2:
        raise ValueError(f"malformed PIR header at line {line_number}")
    after = line[line.rindex('"') + 1:].strip().split()
    if len(after) < 5:
        raise ValueError(f"truncated PIR header at line {line_number}")
    try:
        return (
            int(tokens[1]),
            quoted[1].strip(),
            int(after[0]),
            int(after[3]),
            int(after[4]),
        )
    except ValueError as error:
        raise ValueError(f"non-integer PIR header field at line {line_number}") from error


def _read_exact_pir_int_block(lines, start, count, context):
    """Read exactly ``count`` integer tokens, preserving line boundaries.

    PIR's Fortran reader consumes a fixed number of integer values for each
    record.  Rejecting extra tokens on the terminal line prevents a malformed
    k-vector block from consuming the following operation structure.
    """
    values = []
    line_index = start
    while len(values) < count:
        if line_index >= len(lines):
            raise ValueError(
                f"truncated PIR integer block for {context}: "
                f"expected {count}, got {len(values)}"
            )
        tokens = lines[line_index].strip().split()
        if not tokens:
            raise ValueError(
                f"empty PIR integer block line {line_index + 1} for {context}"
            )
        remaining = count - len(values)
        if len(tokens) > remaining:
            raise ValueError(
                f"extra PIR integer tokens at line {line_index + 1} for {context}: "
                f"expected {remaining}, got {len(tokens)}"
            )
        for token in tokens:
            try:
                values.append(int(token))
            except ValueError as error:
                raise ValueError(
                    f"non-integer PIR token {token!r} at line {line_index + 1} for {context}"
                ) from error
        line_index += 1
    return values, line_index


def _read_exact_pir_float_block(lines, start, count, context):
    """Read exactly ``count`` matrix scalar tokens with provenance checks."""
    values = []
    spellings = []
    line_index = start
    while len(values) < count:
        if line_index >= len(lines):
            raise ValueError(
                f"truncated PIR matrix block for {context}: "
                f"expected {count}, got {len(values)}"
            )
        tokens = lines[line_index].strip().split()
        if not tokens:
            raise ValueError(
                f"empty PIR matrix block line {line_index + 1} for {context}"
            )
        remaining = count - len(values)
        if len(tokens) > remaining:
            raise ValueError(
                f"extra PIR matrix tokens at line {line_index + 1} for {context}: "
                f"expected {remaining}, got {len(tokens)}"
            )
        for token in tokens:
            try:
                values.append(float(token))
            except ValueError as error:
                raise ValueError(
                    f"non-numeric PIR matrix token {token!r} at line "
                    f"{line_index + 1} for {context}"
                ) from error
            spellings.append(token)
        line_index += 1
    return values, spellings, line_index


def _read_pir_operation_row(lines, start, context):
    """Read one official 16-integer Seitz operation row."""
    if start >= len(lines):
        raise ValueError(f"truncated PIR operation row for {context}")
    tokens = lines[start].strip().split()
    if len(tokens) != 16:
        raise ValueError(
            f"PIR operation row at line {start + 1} for {context} has "
            f"{len(tokens)} integers, expected 16"
        )
    try:
        return [int(token) for token in tokens], start + 1
    except ValueError as error:
        raise ValueError(
            f"non-integer PIR operation row at line {start + 1} for {context}"
        ) from error


def _read_pir_irtranslation_row(lines, start, context):
    """Read one required four-integer nonspecial ``irtranslation`` row."""
    if start >= len(lines):
        raise ValueError(f"truncated PIR irtranslation row for {context}")
    tokens = lines[start].strip().split()
    if len(tokens) != 4:
        raise ValueError(
            f"PIR irtranslation row at line {start + 1} for {context} has "
            f"{len(tokens)} integers, expected 4"
        )
    try:
        values = [int(token) for token in tokens]
    except ValueError as error:
        raise ValueError(
            f"non-integer PIR irtranslation row at line {start + 1} for {context}"
        ) from error
    if values[3] <= 0:
        raise ValueError(
            f"invalid PIR irtranslation denominator {values[3]} at line "
            f"{start + 1} for {context}"
        )
    return values, start + 1


def _read_pir_operation_payload(lines, start, dim, special, context):
    """Read one PIR operation, its optional translation, and matrix.

    ``special`` is derived solely from the archived augmented k-vector by the
    official PIR reader rule.  Labels are intentionally not an input here.
    """
    op_nums, next_line = _read_pir_operation_row(lines, start, context)
    irtranslation = None
    if not special:
        irtranslation, next_line = _read_pir_irtranslation_row(
            lines, next_line, context
        )
    matrix_values, spellings, next_line = _read_exact_pir_float_block(
        lines, next_line, dim * dim, context
    )
    return op_nums, irtranslation, matrix_values, spellings, next_line


def _pir_kvector_is_special(kvector_values):
    """Apply PIR_data.f's ``kspecial`` test to the first augmented k-vector."""
    if len(kvector_values) < 16:
        raise ValueError("PIR k-vector record must contain at least 16 integers")
    # Fortran kvector(kp+4*k+m), converted to record-local Python offsets.
    component_offsets = (4, 5, 6, 8, 9, 10, 12, 13, 14)
    return all(kvector_values[offset] == 0 for offset in component_offsets)


def _parse_pir_kvectors():
    """Parse PIR_data.txt and return dict: (SG#, ML_label) -> (kx, ky, kz, denom)."""
    lines = _read_pir_lines()

    kvec_map = {}
    i = 3  # skip 3 header lines
    while i < len(lines):
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break
        sg_num, label, dim, pmkcount, opcount = _parse_pir_header(
            lines[i], i + 1
        )
        if dim <= 0 or pmkcount <= 0 or opcount <= 0:
            raise ValueError(
                f"invalid PIR dimensions/counts at line {i + 1}: "
                f"dim={dim}, pmkcount={pmkcount}, opcount={opcount}"
            )
        key = (sg_num, label)
        if key in kvec_map:
            raise ValueError(f"duplicate PIR source key {key!r}")
        kvals, i = _read_exact_pir_int_block(
            lines, i + 1, pmkcount * 16, f"SG{sg_num} {label!r}"
        )
        kvec_map[key] = (kvals[0], kvals[1], kvals[2], kvals[3])
        special = _pir_kvector_is_special(kvals)
        for op_index in range(opcount):
            _op_nums, _translation, _matrix, _spellings, i = (
                _read_pir_operation_payload(
                    lines,
                    i,
                    dim,
                    special,
                    f"SG{sg_num} {label!r} operation {op_index}",
                )
            )

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


def _format_rust_f64(value):
    """Format a finite float as a Rust ``f64`` literal without changing it.

    Python's default floating-point format may use scientific notation.  Do
    not trim trailing zeroes from the resulting string: in a value such as
    ``1.047197638e-10`` the final zero is part of the exponent, not redundant
    decimal padding.  Trimming it would amplify harmless 1e-10 phase noise to
    1e-1 in the generated database.
    """
    if not math.isfinite(value):
        raise ValueError(f"cannot emit non-finite Rust f64 literal: {value!r}")
    if abs(value) < 1e-15:
        return "0.0"
    if abs(value - round(value)) < 1e-12:
        return f"{int(round(value))}.0"

    literal = f"{value:.10}"
    parsed = float(literal)
    if not math.isclose(parsed, value, rel_tol=5e-10, abs_tol=1e-15):
        raise AssertionError(
            f"Rust f64 formatting changed {value!r} to {literal!r} ({parsed!r})"
        )
    return literal


def _parse_pir_characters():
    """Parse PIR_data.txt and return character, matrix, and operation maps:
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
    kvector_map = {} # (SG#, ML_label) -> all augmented k-vector integers
    irtranslation_rows = 0
    matrix_scalar_tokens = 0
    matrix_token_spellings = set()
    record_count = 0
    i = 3  # skip 3 header lines

    while i < len(lines):
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break

        sg, label, dim, pmkcount, opcount = _parse_pir_header(lines[i], i + 1)
        if dim <= 0 or pmkcount <= 0 or opcount <= 0:
            raise ValueError(
                f"invalid PIR dimensions/counts at line {i + 1}: "
                f"dim={dim}, pmkcount={pmkcount}, opcount={opcount}"
            )
        key = (sg, label)
        if key in chars_map:
            raise ValueError(f"duplicate PIR source key {key!r}")
        kvector_values, i = _read_exact_pir_int_block(
            lines, i + 1, pmkcount * 16, f"SG{sg} {label!r} k-vector"
        )
        kvector_map[key] = kvector_values
        special = _pir_kvector_is_special(kvector_values)

        # Read operator matrices + irrep matrices
        chars = []
        rots = []         # rotation matrices: list of [r00..r22], 9 ints per op
        trans = []        # translations: list of [t0,t1,t2], 3 f64 per op
        all_matrices = []  # flat: op0_row0, op0_row1, ..., op1_row0, ...
        for op_index in range(opcount):
            context = f"SG{sg} {label!r} operation {op_index}"
            op_nums, irtranslation, matrix_vals, spellings, i = (
                _read_pir_operation_payload(lines, i, dim, special, context)
            )
            if irtranslation is not None:
                irtranslation_rows += 1
            matrix_scalar_tokens += len(spellings)
            matrix_token_spellings.update(spellings)
            denom = op_nums[15]
            if denom <= 0:
                raise ValueError(f"invalid PIR operation denominator at {context}")
            r00 = op_nums[0] // denom; r01 = op_nums[1] // denom; r02 = op_nums[2] // denom
            r10 = op_nums[4] // denom; r11 = op_nums[5] // denom; r12 = op_nums[6] // denom
            r20 = op_nums[8] // denom; r21 = op_nums[9] // denom; r22 = op_nums[10] // denom
            rots.append([r00, r01, r02, r10, r11, r12, r20, r21, r22])
            trans.append([
                float(op_nums[3]) / float(denom),
                float(op_nums[7]) / float(denom),
                float(op_nums[11]) / float(denom),
            ])

            # Store raw matrix values (keep full precision)
            all_matrices.extend(matrix_vals)

            # Compute trace (sum of diagonal elements)
            trace = 0.0
            for d in range(dim):
                idx = d * dim + d
                if idx < len(matrix_vals):
                    trace += matrix_vals[idx]
            chars.append(_round_char(trace))

        chars_map[key] = chars
        matrices_map[key] = all_matrices
        rots_map[key] = rots
        trans_map[key] = trans
        dim_map[key] = dim
        record_count += 1

    census = {
        "records": record_count,
        "irtranslation_rows": irtranslation_rows,
        "matrix_scalar_tokens": matrix_scalar_tokens,
        "matrix_token_spellings": frozenset(matrix_token_spellings),
        "unmatched_structural_tokens": 0,
    }
    expected_spellings = frozenset({
        "-1", "-0.96593", "-0.86603", "-0.70711", "-0.68301",
        "-0.61237", "-0.50000", "-0.43301", "-0.35355", "-0.25882",
        "-0.25000", "-0.18301", "0", "0.18301", "0.25000",
        "0.25882", "0.35355", "0.43301", "0.50000", "0.61237",
        "0.68301", "0.70711", "0.86603", "0.96593", "1",
    })
    expected_census = {
        "records": 10294,
        "irtranslation_rows": 64588,
        "matrix_scalar_tokens": 8977752,
        "matrix_token_spellings": expected_spellings,
        "unmatched_structural_tokens": 0,
    }
    if census != expected_census:
        raise ValueError(
            "PIR archive structural census mismatch: "
            f"observed={census!r}, expected={expected_census!r}"
        )
    return chars_map, matrices_map, rots_map, dim_map, trans_map, kvector_map, census


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
    cir_irnumber_map = {}  # (SG#, ML_label) -> stable ISO-IR CIR irnumber
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
        irnumber = int(tokens[0])
        label = quoted[1].strip()
        after_line = s[s.rindex('"') + 1:].strip().split()
        if len(after_line) < 5:
            i += 1
            continue

        dim = int(after_line[0])
        star_count = int(after_line[2]) if len(after_line) > 2 else 1
        if star_count <= 0 or dim <= 0 or dim % star_count != 0:
            raise ValueError(
                f"invalid CIR dimension/star count for SG{sg} {label!r}: "
                f"dim={dim}, star_count={star_count}"
            )
        little_dim = dim // star_count
        opcount = int(after_line[4])

        # CIR stores one augmented 4x4 k-vector record per star arm.  The
        # fourth header field is the number of physically irreducible k
        # vectors, not an additional multiplier.  Using that field here leaves
        # a star-arm record in the stream whenever the two counts differ, where
        # it is then mistaken for the identity operation.
        kvec_remaining = star_count * 16
        while kvec_remaining > 0 and i < len(lines) - 1:
            i += 1
            if i >= len(lines):
                break
            kvec_remaining -= len(lines[i].strip().split())

        # Read operator matrices + complex irrep matrices
        chars = []
        little_chars = []  # trace of the first (stored-k) star-arm block
        rots = []         # rotation matrices: list of [r00,r01,r02,r10,r11,r12,r20,r21,r22]
        trans = []        # fractional translations: list of [t0,t1,t2]
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
                            trans.append([
                                float(op_nums[3]) / float(denom),
                                float(op_nums[7]) / float(denom),
                                float(op_nums[11]) / float(denom),
                            ])
                        except ValueError:
                            pass
                    break
                # Complex value line (conjugate from previous op): skip
                i += 1
            if i >= len(lines):
                break

            # Move to the optional irtranslation / complex irrep matrix line.
            i += 1
            if i >= len(lines):
                break

            # CIR records may contain a four-integer irtranslation even at
            # labels traditionally regarded as special k-points.  Complex
            # matrix values are always written as ``(real,imag)``, so this is
            # unambiguous and does not need a label heuristic.
            next_parts = lines[i].strip().split()
            if (len(next_parts) == 4
                    and all(re.fullmatch(r'-?\d+', value) for value in next_parts)):
                i += 1
                if i >= len(lines):
                    break

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

            little_trace_re = 0.0
            little_trace_im = 0.0
            for d in range(little_dim):
                idx = d * dim + d
                if idx < len(complex_vals):
                    little_trace_re += complex_vals[idx][0]
                    little_trace_im += complex_vals[idx][1]
            little_chars.append((little_trace_re, little_trace_im))

        if len(rots) != opcount or len(trans) != opcount or len(chars) != opcount:
            raise ValueError(
                f"incomplete CIR record SG{sg} {label!r}: "
                f"ops={len(rots)}/{opcount}, translations={len(trans)}/{opcount}, "
                f"characters={len(chars)}/{opcount}"
            )

        key = (sg, label)
        if key not in cir_chars:
            if irnumber in cir_irnumber_map.values():
                raise ValueError(f"duplicate CIR source irnumber {irnumber}")
            cir_chars[key] = {
                'dim': dim,
                'star_count': star_count,
                'little_dim': little_dim,
                'opcount': opcount,
                'chars': chars,  # list of (re, im, rounded_re)
                'little_chars': little_chars,
                'rots': rots,   # list of [r00..r22], 9 ints per op
                'trans': trans, # list of [t0,t1,t2], fractional per op
                'irnumber': irnumber,
            }
            if store_matrices:
                cir_matrices[key] = all_matrices
            cir_irnumber_map[key] = irnumber
        else:
            raise ValueError(f"duplicate CIR source key {key!r}")

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


COMPOUND_NAMING_GRAMMAR_VERSION = 1
COMPOUND_NAMING_PROVENANCE = "ISO-IR Miller-Love concatenation, resolver v1"


def _resolve_compound_constituents(sg, ml, cir_data):
    """Resolve one compound ML spelling to exactly two same-SG CIR records.

    ISO-IR does not carry an independent PIR-to-CIR relation.  The documented
    relation is the Miller--Love concatenation convention, so this resolver is
    the sole generation-time boundary where that convention is interpreted.
    Every component must resolve to one unique CIR key.  Operation alignment
    and accepted generated component data are checked by the caller before
    this identity is frozen into generated metadata.

    ``None`` denotes a label not recognized by the compound grammar.  A
    recognized compound with malformed, unresolved, ambiguous, or more than
    two constituents raises ``ValueError``; later generation checks likewise
    fail hard rather than silently omitting the record.
    """
    parts = _decompose_compound_label(ml)
    if parts is None:
        return None
    if len(parts) != 2:
        raise ValueError(
            f"compound SG{sg} {ml!r} has {len(parts)} constituents; expected exactly two"
        )
    entries = []
    for part in parts:
        entry = cir_data.get((sg, part))
        if entry is None:
            raise ValueError(
                f"compound SG{sg} {ml!r} has unresolved CIR constituent {part!r}"
            )
        entries.append(entry)
    if any(entry.get('irnumber') is None for entry in entries):
        raise ValueError(f"compound SG{sg} {ml!r} has CIR constituent without irnumber")
    if any(entry['opcount'] != entries[0]['opcount'] for entry in entries[1:]):
        raise ValueError(f"compound SG{sg} {ml!r} has mismatched CIR operation counts")
    semantics = "realification" if parts[0] == parts[1] else "distinct_sum"
    return {'parts': parts, 'entries': entries, 'semantics': semantics}


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
    _verify_pinned_archives()
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
    (chars_map, matrices_map, rots_map, pir_dim_map,
     pir_trans_map, pir_kvector_map, pir_census) = _parse_pir_characters()
    print(f"  Parsed {len(chars_map)} character table entries")
    print(f"  Parsed {len(matrices_map)} matrix data entries")
    if {
        key: tuple(value[:4]) for key, value in pir_kvector_map.items()
    } != kvec_map:
        raise ValueError("PIR k-vector parsers disagree on source records")
    print(
        "  PIR structural census: "
        f"{pir_census['records']} records, "
        f"{pir_census['irtranslation_rows']} irtranslation rows, "
        f"{pir_census['matrix_scalar_tokens']} matrix scalars"
    )

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
        "pir_kvector_map": pir_kvector_map,
        "pir_census": pir_census,
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


def _reorder_spin_ops_to_hall(spin_op_rots, spin_op_trans, spin_op_su2,
                                spin_op_sg_start, spin_op_sg_count,
                                spin_lg_op_indices_flat, spin_lg_op_starts,
                                spin_lg_op_counts, sg_hall_choice):
    """Reorder SPIN_OP data per SG from Bilbao order to spglib Hall order."""
    sg_halls = _load_hall_operations()
    if not sg_halls:
        return

    reordered_sgs = 0
    total_spin_ops = len(spin_op_rots) // 9
    # Build new arrays
    new_rots = []
    new_trans = []
    new_su2 = []
    new_sg_start = [0] * 231
    new_sg_count = [0] * 231
    # Build old→new mapping per SG. Both keys and values are SG-local
    # operation indices because spin_lg_op_indices are SG-local.
    sg_bilbao_to_new = {}  # sg_num → {old_local: new_local}

    for sg_num in range(1, 231):
        count = spin_op_sg_count[sg_num]
        if count == 0:
            continue
        old_start = spin_op_sg_start[sg_num]
        # Get canonical Hall for this SG
        hall_info = sg_hall_choice.get(sg_num)
        if hall_info is None:
            # No canonical Hall — keep original order
            new_pos = len(new_rots) // 9
            new_sg_start[sg_num] = new_pos
            new_sg_count[sg_num] = count
            mapping = {i: i for i in range(count)}
            for i in range(count):
                o = old_start + i
                new_rots.extend(spin_op_rots[o*9:(o+1)*9])
                new_trans.extend(spin_op_trans[o*3:(o+1)*3])
                new_su2.extend(spin_op_su2[o*4:(o+1)*4])
            sg_bilbao_to_new[sg_num] = mapping
            continue

        hall_num = hall_info[0]
        hall_rots = None
        for h_num, h_rots, h_trans in sg_halls.get(sg_num, []):
            if h_num == hall_num:
                hall_rots = h_rots
                break

        if hall_rots is None:
            # Fallback: keep original order
            new_pos = len(new_rots) // 9
            new_sg_start[sg_num] = new_pos
            new_sg_count[sg_num] = count
            mapping = {i: new_pos + i for i in range(count)}
            for i in range(count):
                o = old_start + i
                new_rots.extend(spin_op_rots[o*9:(o+1)*9])
                new_trans.extend(spin_op_trans[o*3:(o+1)*3])
                new_su2.extend(spin_op_su2[o*4:(o+1)*4])
            sg_bilbao_to_new[sg_num] = mapping
            continue

        # Build Bilbao→Hall position mapping
        n_hall = len(hall_rots)
        bilbao_to_hall = {}
        hall_to_bilbao = {}
        for bi in range(count):
            b_rot = spin_op_rots[(old_start + bi)*9:(old_start + bi)*9 + 9]
            for hi, h_rot in enumerate(hall_rots):
                if all(b_rot[d] == h_rot[d] for d in range(9)):
                    bilbao_to_hall[bi] = hi
                    hall_to_bilbao[hi] = bi
                    break

        # Build new arrays in Hall order
        new_pos = len(new_rots) // 9
        new_sg_start[sg_num] = new_pos
        holes = 0
        mapping = {}
        for hi in range(n_hall):
            bi = hall_to_bilbao.get(hi)
            if bi is not None:
                o = old_start + bi
                new_rots.extend(spin_op_rots[o*9:(o+1)*9])
                new_trans.extend(spin_op_trans[o*3:(o+1)*3])
                new_su2.extend(spin_op_su2[o*4:(o+1)*4])
                mapping[bi] = hi - holes
            else:
                holes += 1
        new_count = n_hall - holes
        new_sg_count[sg_num] = new_count
        sg_bilbao_to_new[sg_num] = mapping
        if holes > 0:
            reordered_sgs += 1  # count SGs with missing spin ops

    # Replace arrays
    spin_op_rots[:] = new_rots
    spin_op_trans[:] = new_trans
    spin_op_su2[:] = new_su2
    spin_op_sg_start[:] = new_sg_start
    spin_op_sg_count[:] = new_sg_count

    # Update spin_lg_op_indices using the mapping
    for i in range(len(spin_lg_op_indices_flat)):
        old_val = spin_lg_op_indices_flat[i]
        # Find which SG this index belongs to
        # spin_lg_op_indices are per-irrep; each irrep belongs to an SG
        # We need to find the irrep for index i and its SG
        # For now, use a simple approach: scan through spin_lg_op_starts to find the irrep
        # The SG is determined by the spinor_irreps order
        pass  # This is done separately (see below)

    print(f"  SPIN_OP reorder: {len(sg_bilbao_to_new)} SGs processed, "
          f"spin ops reordered to Hall order ({total_spin_ops} → {len(spin_op_rots)//9} total)")
    return sg_bilbao_to_new


def _reorder_to_spglib_order(
        sg, ml, chars_flat, char_starts, char_counts,
        matrices_flat, mat_starts, mat_counts,
        pir_rots_flat, pir_rot_starts, rots_map,
        little_chars_real=None, little_chars_imag=None, little_chars_valid=None,
        pir_trans_flat=None, pir_trans_starts=None,
        spinor_irreps=None, spinor_starts=None, spinor_counts=None,
        cir_comp_flat=None, cir_comp_rots=None, cir_comp_trans=None,
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
    hall_targets = [None] * n_scalar  # per scalar irrep: (Hall rotations, translations)
    sg_hall_choice = {}  # sg → mapping list
    orig_char_counts = list(char_counts)  # Save ISOTROPY sizes before reorder

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
        best_hall_target = None
        best_hall_num = None
        best_exact_count = -1
        pir_trans = []
        if pir_trans_flat is not None and pir_trans_starts is not None:
            trans_start = pir_trans_starts[i]
            raw_trans = pir_trans_flat[trans_start:trans_start + n_ops * 3]
            pir_trans = [raw_trans[j:j + 3] for j in range(0, len(raw_trans), 3)]
        for hall_num, hall_rots, hall_trans in hall_candidates:
            if len(hall_rots) == 0:
                continue
            mapping = []
            rot_cache = {}  # rotation_tuple → ISOTROPY index
            exact_count = 0
            for h_rot, h_trans in zip(hall_rots, hall_trans):
                rot_key = tuple(h_rot)
                rotation_matches = [
                    p_idx for p_idx, p_rot in enumerate(pir_rots)
                    if len(p_rot) == 9 and len(h_rot) == 9
                    and all(h_rot[d] == p_rot[d] for d in range(9))
                ]
                exact_matches = [
                    p_idx for p_idx in rotation_matches
                    if p_idx < len(pir_trans)
                    and all(
                        abs((h_trans[axis] - pir_trans[p_idx][axis])
                            - round(h_trans[axis] - pir_trans[p_idx][axis]))
                        < 1e-8
                        for axis in range(3)
                    )
                ]
                if exact_matches:
                    found = exact_matches[0]
                    exact_count += 1
                elif rot_key in rot_cache:
                    # Centered Hall tables can contain several translations
                    # for one rotation while a primitive source table stores
                    # only one representative.  Reuse its canonical source
                    # operation; Phase C expands the missing cosets.
                    found = rot_cache[rot_key]
                else:
                    found = rotation_matches[0] if rotation_matches else None
                if found is not None:
                    rot_cache[rot_key] = found
                mapping.append(found)
            mapped_op_count = sum(1 for m in mapping if m is not None)
            # Several Hall settings can have the same rotation multiset but
            # different origins.  Prefer the one whose full Seitz operations
            # match the PIR source table most closely.  This removes the
            # origin-setting ambiguity deterministically and preserves the
            # projective multiplication phases at Brillouin-zone boundaries.
            if mapped_op_count == len(hall_rots) and exact_count > best_exact_count:
                best_mapping = mapping
                best_hall_target = (hall_rots, hall_trans)
                best_hall_num = hall_num
                best_exact_count = exact_count

        if best_mapping:
            hall_targets[i] = best_hall_target
            needs_resize = len(best_mapping) != n_ops
            # A centered conventional Hall group can contain more operation
            # representatives than the primitive ISOTROPY table.  Defer all
            # mutations in that case: applying only the first n targets would
            # overwrite source entries that the later expansion still needs.
            if not needs_resize:
                _apply_reorder(chars_flat, char_starts[i], n_ops, best_mapping, 1)
                if little_chars_real is not None:
                    start = char_starts[i]
                    old_re = little_chars_real[start:start + n_ops]
                    old_im = little_chars_imag[start:start + n_ops]
                    old_valid = little_chars_valid[start:start + n_ops]
                    old_trans = pir_trans_flat[
                        pir_trans_starts[i]:pir_trans_starts[i] + n_ops * 3]
                    hall_trans = best_hall_target[1]
                    kx, ky, kz, kd = _lookup_kvec(kvec_map, sg_num, ml[i])
                    for h, source in enumerate(best_mapping):
                        if source is None or source >= n_ops:
                            little_chars_valid[start + h] = 0
                            continue
                        delta = [
                            hall_trans[h][axis] - old_trans[source * 3 + axis]
                            for axis in range(3)
                        ]
                        theta = 2.0 * math.pi * (
                            kx * delta[0] + ky * delta[1] + kz * delta[2]) / kd
                        phase_re = math.cos(theta)
                        phase_im = math.sin(theta)
                        little_chars_real[start + h] = (
                            old_re[source] * phase_re - old_im[source] * phase_im)
                        little_chars_imag[start + h] = (
                            old_re[source] * phase_im + old_im[source] * phase_re)
                        little_chars_valid[start + h] = old_valid[source]
            if not needs_resize:
                dim_sq = mat_counts[i] // n_ops if n_ops else 1
                if dim_sq > 0 and mat_counts[i] > 0:
                    _apply_reorder(matrices_flat, mat_starts[i], n_ops, best_mapping, dim_sq)
                if n_ops > 0:
                    _apply_reorder(pir_rots_flat, pir_rot_starts[i], n_ops, best_mapping, 9)
                    if pir_trans_flat is not None and pir_trans_starts is not None:
                        _apply_reorder(
                            pir_trans_flat, pir_trans_starts[i], n_ops, best_mapping, 3)
            char_counts[i] = len(best_mapping)
            sg_hall_choice[sg_num] = (best_hall_num, best_mapping, best_hall_target[1])
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

            mapping = reorder_results[i]
            if mapping is None:
                continue
            if len(mapping) != cir_ops:
                # Centered Hall expansion is handled in Phase C.  Mutating
                # only the first `cir_ops` destinations here would overwrite
                # source entries that the later full-size expansion still
                # needs (the same deferred-mutation rule as PIR data above).
                continue

            # Reorder selected-arm CIR data for every component.  Full Seitz
            # representatives can differ by a lattice vector between source
            # and Hall tables, so the corresponding Bloch phase is mandatory.
            hall_trans = hall_targets[i][1]
            kvec = _lookup_kvec(kvec_map, sg[i], ml[i])
            for comp in range(n_comp):
                comp_char_start = cir_start + comp * cir_ops * 2
                comp_trans_start = (cir_start // 2) * 3 + comp * cir_ops * 3
                old_chars = cir_comp_flat[
                    comp_char_start:comp_char_start + cir_ops * 2]
                old_trans = cir_comp_trans[
                    comp_trans_start:comp_trans_start + cir_ops * 3]
                for h, source in enumerate(mapping):
                    source_value = complex(
                        old_chars[source * 2], old_chars[source * 2 + 1])
                    source_translation = old_trans[source * 3:source * 3 + 3]
                    value = _phase_character(
                        source_value, hall_trans[h], source_translation, kvec)
                    cir_comp_flat[comp_char_start + h * 2] = value.real
                    cir_comp_flat[comp_char_start + h * 2 + 1] = value.imag
                comp_rot_start = (cir_start // 2) * 9 + comp * cir_ops * 9
                _apply_reorder(cir_comp_rots, comp_rot_start, cir_ops, mapping, 9)
                _apply_reorder(cir_comp_trans, comp_trans_start, cir_ops, mapping, 3)
            cir_reordered += 1

    if cir_reordered > 0:
        print(f"  CIR reorder: {cir_reordered} compound irreps")

    # Spinor irreps: after SPIN_OP reorder (Phase D), spin_lg_op_indices
    # point to Hall-ordered operations. The characters are indexed by local
    # position → spin_lg_op_indices[local] → operation, so no character
    # reordering is needed. Each spinor irrep gets an identity match record.
    for spin_idx, sir in enumerate(spinor_irreps):
        sg_num = sir["sg"]
        n_ops = len(sir.get("op_indices", []))
        if sg_hall_choice.get(sg_num) and n_ops > 0:
            reorder_results.append(list(range(n_ops)))  # identity mapping
            mapped_count += 1
        else:
            reorder_results.append(None)
            unmapped_count += 1

    print(f"  Spglib reorder: {mapped_count} mapped, {unmapped_count} unmapped "
          f"(of {len(reorder_results)} irreps, {len(sg_halls)} SGs)")
    return reorder_results, sg_hall_choice, orig_char_counts, hall_targets


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
            offset = src - start
            dst = start + new_pos * stride
            for d in range(stride):
                arr[dst + d] = old[offset + d]


def _translations_equal_mod_lattice(left, right, tolerance=1e-8):
    """Whether two fractional translations differ by an integer lattice vector."""
    return len(left) == len(right) == 3 and all(
        abs((left[axis] - right[axis]) - round(left[axis] - right[axis])) < tolerance
        for axis in range(3)
    )


def _phase_character(value, target_translation, source_translation, kvec):
    """Move a complex Bloch character between equivalent Seitz representatives."""
    kx, ky, kz, kd = kvec
    delta = [
        target_translation[axis] - source_translation[axis]
        for axis in range(3)
    ]
    theta = 2.0 * math.pi * (
        kx * delta[0] + ky * delta[1] + kz * delta[2]) / kd
    phase = complex(math.cos(theta), math.sin(theta))
    return value * phase


def _align_cir_characters(values, source_rots, source_trans,
                          target_rots, target_trans, kvec):
    """Align one CIR character row to a target full-Seitz operation order."""
    if not (len(values) == len(source_rots) == len(source_trans)):
        return None
    if len(target_rots) != len(target_trans) or len(values) != len(target_rots):
        return None

    aligned = []
    used = set()
    for target_rotation, target_translation in zip(target_rots, target_trans):
        matches = [
            index
            for index, (source_rotation, source_translation)
            in enumerate(zip(source_rots, source_trans))
            if index not in used
            and source_rotation == target_rotation
            and _translations_equal_mod_lattice(
                source_translation, target_translation)
        ]
        if len(matches) != 1:
            return None
        source = matches[0]
        used.add(source)
        aligned.append(_phase_character(
            values[source], target_translation, source_trans[source], kvec))
    return aligned


def _validate_compound_bindings(
        sg, ml, char_counts, pir_rots_flat, pir_rot_starts,
        pir_trans_flat, pir_trans_starts, cir_comp_starts,
        cir_comp_counts, cir_comp_ops, cir_comp_rots, cir_comp_trans):
    """Require every accepted compound row to share the final PIR Seitz order.

    CIR component rotations and generation-only translations are deliberately
    retained as auditable parallel arrays.  Once Hall reorder and padding are
    complete, both must equal the record's PIR arrays entry-for-entry; no
    prefix or rotation-only binding is acceptable.
    """
    accepted = 0
    for i, component_count in enumerate(cir_comp_counts):
        if component_count == 0:
            continue
        accepted += 1
        if pir_rot_starts[i] % 9 != 0:
            raise ValueError(
                f"compound SG{sg[i]} {ml[i]!r} has misaligned PIR rotation offset"
            )
        if pir_trans_starts[i] != (pir_rot_starts[i] // 9) * 3:
            raise ValueError(
                f"compound SG{sg[i]} {ml[i]!r} has PIR rotation/translation offset mismatch"
            )
        op_count = cir_comp_ops[i]
        if op_count <= 0 or op_count != char_counts[i]:
            raise ValueError(
                f"compound SG{sg[i]} {ml[i]!r} has CIR/PIR operation count "
                f"{op_count}/{char_counts[i]} after final ordering"
            )
        pir_rot_start = pir_rot_starts[i]
        pir_trans_start = pir_trans_starts[i]
        if (pir_rot_start + op_count * 9 > len(pir_rots_flat)
                or pir_trans_start + op_count * 3 > len(pir_trans_flat)):
            raise ValueError(
                f"compound SG{sg[i]} {ml[i]!r} has out-of-bounds final PIR Seitz row"
            )
        pir_rot_row = pir_rots_flat[pir_rot_start:pir_rot_start + op_count * 9]
        pir_trans_row = pir_trans_flat[pir_trans_start:pir_trans_start + op_count * 3]
        if len(pir_rot_row) != op_count * 9 or len(pir_trans_row) != op_count * 3:
            raise ValueError(
                f"compound SG{sg[i]} {ml[i]!r} has incomplete final PIR Seitz row"
            )
        cir_rot_base = (cir_comp_starts[i] // 2) * 9
        cir_trans_base = (cir_comp_starts[i] // 2) * 3
        for component in range(component_count):
            cir_rot_start = cir_rot_base + component * op_count * 9
            cir_trans_start = cir_trans_base + component * op_count * 3
            if (cir_rot_start + op_count * 9 > len(cir_comp_rots)
                    or cir_trans_start + op_count * 3 > len(cir_comp_trans)):
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} has out-of-bounds final CIR Seitz row"
                )
            cir_rot_row = cir_comp_rots[cir_rot_start:cir_rot_start + op_count * 9]
            cir_trans_row = cir_comp_trans[cir_trans_start:cir_trans_start + op_count * 3]
            if cir_rot_row != pir_rot_row:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} constituent {component} "
                    "CIR/PIR rotations differ after final ordering"
                )
            if cir_trans_row != pir_trans_row:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} constituent {component} "
                    "CIR/PIR translations differ after final ordering"
                )
    if accepted != 672:
        raise ValueError(
            f"final compound binding census expected 672 accepted records, got {accepted}"
        )


def _validate_pir_storage_alignment(
        sg, ml, char_starts, char_counts, chars_flat, pir_rots_flat,
        pir_rot_starts, pir_trans_flat, pir_trans_starts,
        little_chars_real, little_chars_imag, little_chars_valid):
    """Require every scalar PIR operation row to share one flat offset."""
    for i, op_count in enumerate(char_counts):
        if pir_rot_starts[i] % 9 != 0:
            raise ValueError(
                f"scalar SG{sg[i]} {ml[i]!r} has misaligned PIR rotation offset"
            )
        expected_translation_start = (pir_rot_starts[i] // 9) * 3
        if pir_trans_starts[i] != expected_translation_start:
            raise ValueError(
                f"scalar SG{sg[i]} {ml[i]!r} has PIR rotation/translation "
                "offset mismatch"
            )
        if char_starts[i] != pir_rot_starts[i] // 9:
            raise ValueError(
                f"scalar SG{sg[i]} {ml[i]!r} has character/PIR offset mismatch"
            )
        if (char_starts[i] + op_count > len(chars_flat)
                or char_starts[i] + op_count > len(little_chars_real)
                or char_starts[i] + op_count > len(little_chars_imag)
                or char_starts[i] + op_count > len(little_chars_valid)
                or pir_rot_starts[i] + op_count * 9 > len(pir_rots_flat)
                or expected_translation_start + op_count * 3 > len(pir_trans_flat)):
            raise ValueError(
                f"scalar SG{sg[i]} {ml[i]!r} has incomplete PIR operation row"
            )


def _apply_padding_plans(padding_plans, chars_flat, char_starts, char_counts,
                          little_chars_real, little_chars_imag, little_chars_valid,
                          sg, ml, kvec_map,
                          matrices_flat, mat_starts, mat_counts,
                          pir_rots_flat, pir_rot_starts,
                          pir_trans_flat, pir_trans_starts,
                          cir_comp_flat, cir_comp_rots, cir_comp_trans,
                          cir_comp_starts, cir_comp_counts, cir_comp_ops,
                          spinor_starts, spinor_counts,
                          reorder_map_per_irrep=None,
                          orig_char_counts=None,
                          hall_targets=None):
    """Rebuild flat arrays with padded entries for compound irreps expanded to Hall size."""
    n_scalar = len(char_starts)
    plans_by_idx = {i: (hall_ops, cir_to_hall) for i, hall_ops, cir_to_hall in padding_plans}

    def _needs_resize(i):
        """True if a mapped entry changes operation count in Hall order."""
        if reorder_map_per_irrep is None or orig_char_counts is None:
            return False
        m = reorder_map_per_irrep[i]
        if m is None:
            return False
        return len(m) != orig_char_counts[i]

    # Build resize plans for mapped entries whose Hall and source sizes differ.
    resize_plans = {}
    for i in range(n_scalar):
        if i not in plans_by_idx and _needs_resize(i):
            resize_plans[i] = reorder_map_per_irrep[i]

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
            if i in resize_plans:
                old_n = orig_char_counts[i]
                old = chars_flat[char_starts[i]:char_starts[i] + old_n]
                mapping = resize_plans[i]
                for h in range(n):
                    ci = mapping[h]
                    if ci is not None and ci < old_n:
                        new_chars.append(old[ci])
                    else:
                        new_chars.append(0.0)
            else:
                new_chars.extend(chars_flat[char_starts[i]:char_starts[i] + n])
            new_char_counts.append(n)

    def _rebuild_parallel_character_array(values, fill):
        rebuilt = []
        for i in range(n_scalar):
            if i in plans_by_idx:
                hall_ops, cir_to_hall = plans_by_idx[i]
                old = values[char_starts[i]:char_starts[i] + char_counts[i]]
                for h in range(hall_ops):
                    ci = next(
                        (cii for cii, hi in enumerate(cir_to_hall) if hi == h), None)
                    rebuilt.append(old[ci] if ci is not None and ci < len(old) else fill)
            elif i in resize_plans:
                old_n = orig_char_counts[i]
                old = values[char_starts[i]:char_starts[i] + old_n]
                for source in resize_plans[i]:
                    rebuilt.append(
                        old[source] if source is not None and source < old_n else fill)
            else:
                n = char_counts[i]
                rebuilt.extend(values[char_starts[i]:char_starts[i] + n])
        return rebuilt

    new_little_real = _rebuild_parallel_character_array(little_chars_real, 0.0)
    new_little_imag = _rebuild_parallel_character_array(little_chars_imag, 0.0)
    new_little_valid = _rebuild_parallel_character_array(little_chars_valid, 0)
    for i, mapping in resize_plans.items():
        target = hall_targets[i] if hall_targets is not None else None
        if target is None:
            continue
        _hall_rots, hall_trans = target
        old_n = orig_char_counts[i]
        old_re = little_chars_real[char_starts[i]:char_starts[i] + old_n]
        old_im = little_chars_imag[char_starts[i]:char_starts[i] + old_n]
        old_valid = little_chars_valid[char_starts[i]:char_starts[i] + old_n]
        old_trans = pir_trans_flat[
            pir_trans_starts[i]:pir_trans_starts[i] + old_n * 3]
        kx, ky, kz, kd = _lookup_kvec(kvec_map, sg[i], ml[i])
        new_start = new_char_starts[i]
        for h, source in enumerate(mapping):
            if source is None or source >= old_n:
                new_little_valid[new_start + h] = 0
                continue
            delta = [
                hall_trans[h][axis] - old_trans[source * 3 + axis]
                for axis in range(3)
            ]
            theta = 2.0 * math.pi * (
                kx * delta[0] + ky * delta[1] + kz * delta[2]) / kd
            phase_re = math.cos(theta)
            phase_im = math.sin(theta)
            new_little_real[new_start + h] = (
                old_re[source] * phase_re - old_im[source] * phase_im)
            new_little_imag[new_start + h] = (
                old_re[source] * phase_im + old_im[source] * phase_re)
            new_little_valid[new_start + h] = old_valid[source]

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
            if i in resize_plans:
                old_ops = orig_char_counts[i]
                mapping = resize_plans[i]
                old_m = matrices_flat[mat_starts[i]:mat_starts[i] + mat_counts[i]]
                dim_sq = mat_counts[i] // old_ops if old_ops else 0
                if dim_sq > 0:
                    for source in mapping:
                        if source is not None and source < old_ops:
                            new_mats.extend(old_m[source * dim_sq:(source + 1) * dim_sq])
                        else:
                            new_mats.extend([0.0] * dim_sq)
                    new_mat_counts.append(len(mapping) * dim_sq)
                else:
                    new_mat_counts.append(0)
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
        elif i < n_scalar and i in resize_plans:
            mapping = resize_plans[i]
            old_ops = orig_char_counts[i]
            old_r = pir_rots_flat[pir_rot_starts[i]:pir_rot_starts[i] + old_ops * 9]
            for source in mapping:
                if source is not None and source < old_ops:
                    new_rots.extend(old_r[source * 9:(source + 1) * 9])
                else:
                    new_rots.extend([0] * 9)
        else:
            new_rots.extend(pir_rots_flat[pir_rot_starts[i]:orig_end])

    # Rebuild scalar PIR translations in lockstep with PIR rotations.  This is
    # deliberately separate from the rotation offsets: Phase-C expansion can
    # change their lengths, and deriving one flat-array offset from the other
    # is valid only when every scalar record remains exactly aligned.
    new_trans = []
    new_trans_starts = []
    for i in range(n_scalar):
        new_trans_starts.append(len(new_trans))
        old_end = (pir_trans_starts[i + 1]
                   if i + 1 < n_scalar else len(pir_trans_flat))
        if i in plans_by_idx:
            hall_ops, cir_to_hall = plans_by_idx[i]
            old_ops = len(cir_to_hall)
            old_t = pir_trans_flat[pir_trans_starts[i]:pir_trans_starts[i] + old_ops * 3]
            for h in range(hall_ops):
                ci = next((cii for cii, hi in enumerate(cir_to_hall) if hi == h), None)
                if ci is not None and ci < old_ops:
                    new_trans.extend(old_t[ci * 3:(ci + 1) * 3])
                else:
                    new_trans.extend([0.0] * 3)
        elif i in resize_plans:
            target = hall_targets[i] if hall_targets is not None else None
            if target is None:
                raise ValueError(f"missing Hall target for resized irrep {i}")
            _hall_rots, hall_trans = target
            if len(hall_trans) != len(resize_plans[i]):
                raise ValueError(f"Hall translation count mismatch for irrep {i}")
            for translation in hall_trans:
                new_trans.extend(translation)
        else:
            new_trans.extend(pir_trans_flat[pir_trans_starts[i]:old_end])

    # Rebuild cir_comp_flat and cir_comp_rots
    new_cir_flat = []
    new_cir_rots = []
    new_cir_trans = []
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
                old_trans_start = (cir_comp_starts[i] // 2) * 3 + comp * old_ops * 3
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
                        new_cir_trans.extend(cir_comp_trans[old_trans_start + ci * 3:old_trans_start + (ci + 1) * 3])
                    else:
                        new_cir_flat.append(0.0)
                        new_cir_flat.append(0.0)
                        new_cir_rots.extend([0] * 9)
                        new_cir_trans.extend([0.0] * 3)
            cir_comp_ops[i] = hall_ops
        else:
            # Check if this mapped entry needs CIR expansion too
            if i in resize_plans:
                mapping = resize_plans[i]
                hall_ops = len(mapping)
                target = hall_targets[i] if hall_targets is not None else None
                if target is None:
                    raise ValueError(f"missing Hall target for resized CIR irrep {i}")
                hall_rots, hall_trans = target
                kvec = _lookup_kvec(kvec_map, sg[i], ml[i])
                for comp in range(n_comp):
                    old_start = cir_comp_starts[i] + comp * old_ops * 2
                    old_rot_start = (cir_comp_starts[i] // 2) * 9 + comp * old_ops * 9
                    old_trans_start = (cir_comp_starts[i] // 2) * 3 + comp * old_ops * 3
                    for h in range(hall_ops):
                        ci = mapping[h]
                        if ci is not None and ci < old_ops:
                            source_value = complex(
                                cir_comp_flat[old_start + ci * 2],
                                cir_comp_flat[old_start + ci * 2 + 1])
                            source_translation = cir_comp_trans[
                                old_trans_start + ci * 3:old_trans_start + (ci + 1) * 3]
                            value = _phase_character(
                                source_value, hall_trans[h], source_translation, kvec)
                            new_cir_flat.extend([value.real, value.imag])
                            new_cir_rots.extend(hall_rots[h])
                            new_cir_trans.extend(hall_trans[h])
                        else:
                            new_cir_flat.append(0.0)
                            new_cir_flat.append(0.0)
                            new_cir_rots.extend([0] * 9)
                            new_cir_trans.extend([0.0] * 3)
                cir_comp_ops[i] = hall_ops
            else:
                old_start = cir_comp_starts[i]
                total_chars = n_comp * old_ops * 2
                new_cir_flat.extend(cir_comp_flat[old_start:old_start + total_chars])
                old_rot_start = (cir_comp_starts[i] // 2) * 9
                total_rots = n_comp * old_ops * 9
                new_cir_rots.extend(cir_comp_rots[old_rot_start:old_rot_start + total_rots])
                old_trans_start = (cir_comp_starts[i] // 2) * 3
                total_trans = n_comp * old_ops * 3
                new_cir_trans.extend(
                    cir_comp_trans[old_trans_start:old_trans_start + total_trans])

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
    little_chars_real[:] = new_little_real
    little_chars_imag[:] = new_little_imag
    little_chars_valid[:] = new_little_valid
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
    pir_trans_flat[:] = new_trans
    pir_trans_starts[:] = new_trans_starts
    cir_comp_flat[:] = new_cir_flat
    cir_comp_rots[:] = new_cir_rots
    cir_comp_trans[:] = new_cir_trans
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

    # The selected-arm dimension for an ordinary scalar is an independent
    # CIR-header fact, not something inferred from a character value.  A
    # compound has no single CIR header for the physical row and is kept at
    # zero here; its semantic dimension is validated from both constituents
    # below.
    compound_resolutions = [
        _resolve_compound_constituents(sg[i], ml[i], cir_data)
        for i in range(len(ml))
    ]
    scalar_selected_dims = []
    scalar_full_dims = []
    for i, (sg_num, ml_label) in enumerate(zip(sg, ml)):
        resolved = compound_resolutions[i]
        if resolved is not None:
            scalar_selected_dims.append(0)
            entries = resolved["entries"]
            if resolved["semantics"] == "realification":
                full_dim = 2 * entries[0]["dim"]
                if entries[1]["dim"] != entries[0]["dim"]:
                    raise ValueError(
                        f"realification SG{sg_num} {ml_label!r} has unequal "
                        "constituent full dimensions"
                    )
            else:
                full_dim = sum(entry["dim"] for entry in entries)
            if full_dim <= 0:
                raise ValueError(
                    f"compound SG{sg_num} {ml_label!r} has invalid full dimension"
                )
            scalar_full_dims.append(full_dim)
            continue
        entry = cir_data.get((sg_num, ml_label))
        if entry is None:
            raise ValueError(
                f"ordinary scalar SG{sg_num} {ml_label!r} has no exact CIR header"
            )
        little_dim = entry["little_dim"]
        if little_dim <= 0 or little_dim > 255:
            raise ValueError(
                f"ordinary scalar SG{sg_num} {ml_label!r} has invalid "
                f"authoritative selected dimension {little_dim}"
            )
        scalar_selected_dims.append(little_dim)
        scalar_full_dims.append(entry["dim"])

    # PIR rotation matrices for H_ops → PIR order mapping (Wigner test)
    rots_map = data.get("rots_map", {})

    # ── Build flat CHARACTERS array and per-irrep start/count ──
    chars_flat = []
    little_chars_real = []
    little_chars_imag = []
    little_chars_valid = []
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
        cir_entry = cir_data.get((sg[i], ml[i]))
        little_cir = cir_entry.get('little_chars', []) if cir_entry else []
        cir_rots = cir_entry.get('rots', []) if cir_entry else []
        cir_trans = cir_entry.get('trans', []) if cir_entry else []
        pir_rots = rots_map.get((sg[i], ml[i]), [])
        pir_trans = data.get("pir_trans_map", {}).get((sg[i], ml[i]), [])
        little = []
        if (len(little_cir) == len(cir_rots) == len(cir_trans)
                == len(pir_rots) == len(pir_trans) == len(ch)
                and len(ch) > 0):
            for pir_rotation, pir_translation in zip(pir_rots, pir_trans):
                matches = [
                    index
                    for index, (cir_rotation, cir_translation)
                    in enumerate(zip(cir_rots, cir_trans))
                    if cir_rotation == pir_rotation
                    and all(
                        abs((cir_translation[axis] - pir_translation[axis])
                            - round(cir_translation[axis] - pir_translation[axis]))
                        < 1e-8
                        for axis in range(3)
                    )
                ]
                if len(matches) != 1:
                    little = []
                    break
                little.append(little_cir[matches[0]])
        valid = len(little) == len(ch) and len(ch) > 0
        if valid:
            little_chars_real.extend(value[0] for value in little)
            little_chars_imag.extend(value[1] for value in little)
            little_chars_valid.extend([1] * len(ch))
        else:
            little_chars_real.extend([0.0] * len(ch))
            little_chars_imag.extend([0.0] * len(ch))
            little_chars_valid.extend([0] * len(ch))
        if not ch:
            missing_chars += 1
    if missing_chars > 0:
        print(f"  Warning: {missing_chars}/{len(ml)} irreps have no character data")
    if cir_filled > 0:
        print(f"  (CIR fallback filled {cir_filled} character tables)")
    valid_little_tables = sum(
        1 for i in range(len(ml))
        if char_counts[i] > 0 and little_chars_valid[char_starts[i]] == 1
    )
    print(f"  CIR selected-arm characters: {valid_little_tables}/{len(ml)} tables")

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
    cir_comp_trans = []  # translations (3 floats per op), generation-time only
    cir_comp_starts = []  # per-irrep start index (0 = not compound)
    cir_comp_counts = []  # number of CIR components (0 = not compound)
    cir_comp_ops = []     # operations per CIR component
    cir_comp_total = 0
    compound_metadata = []
    compound_metadata_indices = [0] * len(ml)
    for i in range(len(ml)):
        resolved = compound_resolutions[i]
        if resolved is not None:
            parts = resolved['parts']
            comp_chars = []
            comp_full_chars = []
            comp_rots = []
            comp_trans = []
            n_ops = 0
            entries = resolved['entries']
            n_ops = entries[0]['opcount']

            # Prefer the compound PIR's exact Seitz order.  A small number of
            # compound records have no PIR operator list; their first CIR
            # component is the canonical source order and is copied into the
            # PIR parallel arrays below.
            target_rots = rots_map.get((sg[i], ml[i]), [])
            target_trans = pir_trans_map.get((sg[i], ml[i]), [])
            uses_pir_order = (
                len(target_rots) == len(target_trans) == n_ops and n_ops > 0)
            if entries and not uses_pir_order:
                target_rots = entries[0].get('rots', [])
                target_trans = entries[0].get('trans', [])
            if not (len(target_rots) == len(target_trans) == n_ops and n_ops > 0):
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} has no complete target "
                    f"Seitz order for {n_ops} operations"
                )

            kvec = _lookup_kvec(kvec_map, sg[i], ml[i])
            if entries:
                for p_idx, (part, entry) in enumerate(zip(parts, entries)):
                    little_values = [complex(re_val, im_val)
                                     for re_val, im_val in entry['little_chars']]
                    full_values = [complex(re_val, im_val)
                                   for re_val, im_val, _ in entry['chars']]
                    aligned_little = _align_cir_characters(
                        little_values, entry.get('rots', []), entry.get('trans', []),
                        target_rots, target_trans, kvec)
                    aligned_full = _align_cir_characters(
                        full_values, entry.get('rots', []), entry.get('trans', []),
                        target_rots, target_trans, kvec)
                    if aligned_little is None or aligned_full is None:
                        raise ValueError(
                            f"compound SG{sg[i]} {ml[i]!r} CIR constituent "
                            f"{part!r} failed Seitz operation alignment"
                        )

                    # Repeated CIR labels denote a representation together
                    # with its complex conjugate (for example P3P3).
                    if p_idx > 0 and resolved["semantics"] == "realification":
                        aligned_little = [value.conjugate() for value in aligned_little]
                        aligned_full = [value.conjugate() for value in aligned_full]

                    for value in aligned_little:
                        comp_chars.extend([value.real, value.imag])
                    comp_full_chars.append(aligned_full)
                    for rotation, translation in zip(target_rots, target_trans):
                        comp_rots.extend(rotation)
                        comp_trans.extend(translation)

            # Validate the stored selected-arm rows against their dimensions,
            # and the corresponding full CIR sums against every PIR character.
            if not comp_chars:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} produced no CIR character data"
                )
            pir_ch = _lookup_chars(chars_map, sg[i], ml[i], kvec_map)
            identity_rotation = [1, 0, 0, 0, 1, 0, 0, 0, 1]
            identities = [
                op for op, (rotation, translation)
                in enumerate(zip(target_rots, target_trans))
                if rotation == identity_rotation
                and all(abs(value - round(value)) < 1e-8 for value in translation)
            ]
            if len(identities) != 1:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} has {len(identities)} "
                    "identity operations after Seitz alignment"
                )
            identity = identities[0]
            for component, entry in enumerate(entries):
                value = complex(
                    comp_chars[2 * (component * n_ops + identity)],
                    comp_chars[2 * (component * n_ops + identity) + 1])
                if (abs(value.real - entry['little_dim']) > 0.01
                        or abs(value.imag) > 0.01):
                    raise ValueError(
                        f"compound SG{sg[i]} {ml[i]!r} constituent {component} "
                        "fails identity/dimension validation"
                    )

            if not pir_ch or len(pir_ch) != n_ops:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} has no PIR character row "
                    f"aligned to {n_ops} operations"
                )
            for op in range(n_ops):
                value = sum(component[op] for component in comp_full_chars)
                if (abs(value.real - pir_ch[op]) > 0.01
                        or abs(value.imag) > 0.01):
                    raise ValueError(
                        f"compound SG{sg[i]} {ml[i]!r} fails full character-sum "
                        f"validation at operation {op}"
                    )

            cir_comp_starts.append(len(cir_comp_flat))
            cir_comp_counts.append(len(parts))
            cir_comp_ops.append(n_ops)
            cir_comp_flat.extend(comp_chars)
            # Always extend rotations to keep arrays in sync.
            total_rots_needed = len(parts) * n_ops * 9
            if len(comp_rots) != total_rots_needed:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} has {len(comp_rots)} "
                    f"rotation values; expected {total_rots_needed}"
                )
            cir_comp_rots.extend(comp_rots)
            total_trans_needed = len(parts) * n_ops * 3
            if len(comp_trans) != total_trans_needed:
                raise ValueError(
                    f"compound SG{sg[i]} {ml[i]!r} has {len(comp_trans)} "
                    f"translation values; expected {total_trans_needed}"
                )
            cir_comp_trans.extend(comp_trans)

            # Make previously rotation-less compound PIR records part of
            # the same Hall-order pipeline as every other scalar record.
            if not uses_pir_order and char_counts[i] == n_ops:
                rots_map[(sg[i], ml[i])] = [list(rotation) for rotation in target_rots]
                rot_start = pir_rot_starts[i]
                trans_start = pir_trans_starts[i]
                for op, (rotation, translation) in enumerate(
                        zip(target_rots, target_trans)):
                    pir_rots_flat[rot_start + op * 9:rot_start + (op + 1) * 9] = rotation
                    pir_trans_flat[trans_start + op * 3:trans_start + (op + 1) * 3] = translation
            cir_comp_total += 1
            metadata_index = len(compound_metadata) + 1
            compound_metadata_indices[i] = metadata_index
            compound_metadata.append({
                "sg": sg[i],
                "record_label": ml[i],
                "cir_irnumbers": [entry['irnumber'] for entry in entries],
                "cir_labels": parts,
                "cir_dimensions": [entry['little_dim'] for entry in entries],
                "semantics": 1 if resolved["semantics"] == "realification" else 2,
            })
        else:
            cir_comp_starts.append(0); cir_comp_counts.append(0); cir_comp_ops.append(0)
    resolved_count = sum(resolved is not None for resolved in compound_resolutions)
    if resolved_count != 672 or cir_comp_total != resolved_count:
        raise ValueError(
            f"compound census expected 672 resolved/accepted records, got "
            f"{resolved_count} resolved and {cir_comp_total} accepted"
        )
    realification_count = sum(
        metadata["semantics"] == 1 for metadata in compound_metadata
    )
    distinct_count = sum(
        metadata["semantics"] == 2 for metadata in compound_metadata
    )
    if (realification_count, distinct_count) != (153, 519):
        raise ValueError(
            "compound semantics census expected 153 realifications and "
            f"519 distinct sums, got {realification_count} and {distinct_count}"
        )
    source_refs = 0
    for metadata in compound_metadata:
        sg_value = metadata["sg"]
        for cir_id, cir_label, cir_dim in zip(
                metadata["cir_irnumbers"], metadata["cir_labels"],
                metadata["cir_dimensions"]):
            entry = cir_data.get((sg_value, cir_label))
            if (entry is None or entry["irnumber"] != cir_id
                    or entry["little_dim"] != cir_dim or cir_dim <= 0):
                raise ValueError(
                    f"compound SG{sg_value} {metadata['record_label']!r} has "
                    f"invalid CIR source reference {cir_id}/{cir_label}/{cir_dim}"
                )
            source_refs += 1
    if source_refs != 1344:
        raise ValueError(
            f"compound CIR source-reference census expected 1344, got {source_refs}"
        )
    print(f"  CIR component chars: {cir_comp_total} compound irreps accepted, {len(cir_comp_flat)} values, {len(cir_comp_rots)} rotation ints")
    print(f"  Compound census: {resolved_count} resolved/accepted, {realification_count} realifications, {distinct_count} distinct sums, {source_refs} CIR refs")

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
    spin_extra_flat = []   # imaginary parts of spinor characters
    spin_extra_starts = [] # per-irrep start index
    spin_extra_counts = [] # number of imaginary character values
    for sir in spinor_irreps:
        spinor_starts.append(len(chars_flat))
        chars_real = sir["characters"]
        chars_imag = sir.get("characters_imag", [0.0] * len(chars_real))
        if len(chars_real) != len(chars_imag):
            raise ValueError(
                f"spinor character length mismatch for SG{sir['sg']} "
                f"{sir['ml_label']}: {len(chars_real)} real, {len(chars_imag)} imag"
            )
        spinor_counts.append(len(chars_real))
        chars_flat.extend(chars_real)
        spin_extra_starts.append(len(spin_extra_flat))
        spin_extra_counts.append(len(chars_imag))
        spin_extra_flat.extend(chars_imag)
        # No matrices for spinor irreps
        mat_starts.append(len(matrices_flat))
        mat_counts.append(0)
    if spinor_irreps:
        n_complex = sum(
            1 for sir in spinor_irreps
            if any(abs(v) > 1e-12 for v in sir.get("characters_imag", []))
        )
        print(f"  Added {len(spinor_irreps)} spinor irreps ({n_complex} with complex chars)")

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
    reorder_map_per_irrep, sg_hall_choice, orig_char_counts, hall_targets = _reorder_to_spglib_order(
        sg, ml, chars_flat, char_starts, char_counts,
        matrices_flat, mat_starts, mat_counts,
        pir_rots_flat, pir_rot_starts, rots_map,
        little_chars_real=little_chars_real,
        little_chars_imag=little_chars_imag,
        little_chars_valid=little_chars_valid,
        pir_trans_flat=pir_trans_flat, pir_trans_starts=pir_trans_starts,
        spinor_irreps=spinor_irreps, spinor_starts=spinor_starts,
        spinor_counts=spinor_counts,
        cir_comp_flat=cir_comp_flat, cir_comp_rots=cir_comp_rots,
        cir_comp_trans=cir_comp_trans,
        cir_comp_starts=cir_comp_starts, cir_comp_counts=cir_comp_counts,
        cir_comp_ops=cir_comp_ops,
        kvec_map=kvec_map)
    # CIR component data is also reordered in-place.
    # CIR components are selected-arm complex characters in data-Hall order.
    # Runtime consumes this order directly; CIR rotations remain as an older-
    # generated-data compatibility fallback.
    # reorder_map_per_irrep[i] = None (unmapped) or list[h_idx→pir_idx] (mapped)
    # For spinor irreps: entries past len(ml) are the spinor reorder maps

    # ── Phase C: CIR padding for unmapped compound irreps ──
    padding_plans = _build_padding_plans(
        sg, ml, cir_comp_starts, cir_comp_counts, cir_comp_ops, cir_comp_rots,
        reorder_map_per_irrep)
    # Data already mapped in Phase B stays in place.  Centered conventional
    # Hall groups are expanded here, with Bloch phases applied for every
    # changed Seitz representative.

    resize_count = sum(
        1 for i, mapping in enumerate(reorder_map_per_irrep[:len(ml)])
        if mapping is not None and len(mapping) != orig_char_counts[i]
    )
    if padding_plans or resize_count:
        print(f"  CIR padding: {len(padding_plans)} entries; "
              f"Hall resize: {resize_count} entries")
        _apply_padding_plans(padding_plans, chars_flat, char_starts, char_counts,
                             little_chars_real, little_chars_imag, little_chars_valid,
                             sg, ml, kvec_map,
                             matrices_flat, mat_starts, mat_counts,
                             pir_rots_flat, pir_rot_starts,
                             pir_trans_flat, pir_trans_starts,
                             cir_comp_flat, cir_comp_rots, cir_comp_trans,
                             cir_comp_starts,
                             cir_comp_counts, cir_comp_ops,
                             spinor_starts, spinor_counts,
                             reorder_map_per_irrep=reorder_map_per_irrep,
                             orig_char_counts=orig_char_counts,
                             hall_targets=hall_targets)

    # The final arrays are now in data-Hall order.  Freeze and verify the
    # operation binding only after every reorder/padding mutation is complete.
    _validate_pir_storage_alignment(
        sg, ml, char_starts, char_counts, chars_flat, pir_rots_flat,
        pir_rot_starts, pir_trans_flat, pir_trans_starts,
        little_chars_real, little_chars_imag, little_chars_valid)
    _validate_compound_bindings(
        sg, ml, char_counts, pir_rots_flat, pir_rot_starts,
        pir_trans_flat, pir_trans_starts, cir_comp_starts,
        cir_comp_counts, cir_comp_ops, cir_comp_rots, cir_comp_trans)

    # Ordinary scalar selected dimensions come from CIR headers above.  The
    # final Hall-ordered little character at the unique full-Seitz identity is
    # only an independent cross-check of that authoritative value.
    identity_rotation = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    for i, (sg_num, ml_label) in enumerate(zip(sg, ml)):
        if compound_resolutions[i] is not None:
            continue
        expected = scalar_selected_dims[i]
        op_count = char_counts[i]
        pir_rot_start = pir_rot_starts[i]
        pir_trans_start = pir_trans_starts[i]
        identities = []
        for op in range(op_count):
            rotation = pir_rots_flat[pir_rot_start + op * 9:pir_rot_start + (op + 1) * 9]
            translation = pir_trans_flat[pir_trans_start + op * 3:pir_trans_start + (op + 1) * 3]
            if rotation == identity_rotation and all(
                    abs(value - round(value)) <= 1e-8 for value in translation):
                identities.append(op)
        if len(identities) != 1:
            raise ValueError(
                f"ordinary scalar SG{sg_num} {ml_label!r} has "
                f"{len(identities)} final Seitz identities"
            )
        little_start = char_starts[i]
        little_end = little_start + op_count
        if little_end > len(little_chars_real) or little_end > len(little_chars_imag):
            raise ValueError(
                f"ordinary scalar SG{sg_num} {ml_label!r} has incomplete final little row"
            )
        identity = identities[0]
        little_re = little_chars_real[little_start + identity]
        little_im = little_chars_imag[little_start + identity]
        if abs(little_re - expected) > 1e-8 or abs(little_im) > 1e-8:
            raise ValueError(
                f"ordinary scalar SG{sg_num} {ml_label!r} CIR selected dimension "
                f"{expected} disagrees with final χ(E)=({little_re},{little_im})"
            )

    # ── Phase D: Reorder SPIN_OP data to spglib Hall order ──
    sg_bilbao_to_new = _reorder_spin_ops_to_hall(
        spin_op_rots, spin_op_trans, spin_op_su2,
        spin_op_sg_start, spin_op_sg_count,
        spin_lg_op_indices_flat, spin_lg_op_starts, spin_lg_op_counts,
        sg_hall_choice)

    # Update SG-local spin_lg_op_indices using the old→new local mapping.
    if sg_bilbao_to_new:
        updated_count = 0
        for i, sir in enumerate(spinor_irreps):
            sg_num = sir['sg']
            mapping = sg_bilbao_to_new.get(sg_num, {})
            if not mapping:
                continue
            start = spin_lg_op_starts[i]
            count = spin_lg_op_counts[i]
            for j in range(count):
                old_local = spin_lg_op_indices_flat[start + j]
                new_local = mapping.get(old_local, old_local)
                spin_lg_op_indices_flat[start + j] = new_local
                updated_count += 1
        print(f"  Updated {updated_count} spin_lg_op_indices after SPIN_OP reorder")

    # Spinor rows are a strict indexed view: all four counts are identical,
    # indices are a permutation of distinct SG-local operations, and every SG
    # table has complete rotation/translation/SU(2) parallel arrays.
    for sg_num in range(1, 231):
        ops = spinor_ops_data.get(sg_num, [])
        if (len(spin_op_rots) % 9 != 0 or len(spin_op_trans) % 3 != 0
                or len(spin_op_su2) % 4 != 0):
            raise ValueError("global spin operation arrays have invalid widths")
        sg_start = spin_op_sg_start[sg_num]
        sg_count = spin_op_sg_count[sg_num]
        if (len(ops) != sg_count or sg_start + sg_count > len(spin_op_rots) // 9
                or sg_start + sg_count > len(spin_op_trans) // 3
                or sg_start + sg_count > len(spin_op_su2) // 4):
            raise ValueError(f"spin operation table SG{sg_num} has inconsistent widths")

    for idx, sir in enumerate(spinor_irreps):
        count = spinor_counts[idx]
        if not (count > 0 and count == spin_lg_counts[idx]
                == spin_lg_op_counts[idx] == spin_extra_counts[idx]):
            raise ValueError(
                f"spinor SG{sir['sg']} {sir['ml_label']!r} has inconsistent "
                "character/index/imaginary counts"
            )
        start = spin_lg_op_starts[idx]
        indices = spin_lg_op_indices_flat[start:start + count]
        sg_count = spin_op_sg_count[sir['sg']]
        if len(indices) != count or len(set(indices)) != count or any(
                operation < 0 or operation >= sg_count for operation in indices):
            raise ValueError(
                f"spinor SG{sir['sg']} {sir['ml_label']!r} has invalid little-group indices"
            )
        sg_start = spin_op_sg_start[sir['sg']]
        identities = []
        for position, operation in enumerate(indices):
            operation = int(operation)
            rotation = spin_op_rots[(sg_start + operation) * 9:(sg_start + operation + 1) * 9]
            translation = spin_op_trans[(sg_start + operation) * 3:(sg_start + operation + 1) * 3]
            su2 = spin_op_su2[(sg_start + operation) * 4:(sg_start + operation + 1) * 4]
            if (rotation == identity_rotation
                    and all(abs(value - round(value)) <= 1e-8 for value in translation)
                    and all(abs(value - expected) <= 1e-8
                             for value, expected in zip(su2, [1.0, 0.0, 0.0, 0.0]))):
                identities.append(position)
        if len(identities) != 1:
            raise ValueError(
                f"spinor SG{sir['sg']} {sir['ml_label']!r} has "
                f"{len(identities)} canonical spin identities"
            )
        identity_value = chars_flat[
            spinor_starts[idx] + identities[0]]
        identity_imag = spin_extra_flat[spin_extra_starts[idx] + identities[0]]
        if abs(identity_value - sir['dim']) > 1e-8 or abs(identity_imag) > 1e-8:
            raise ValueError(
                f"spinor SG{sir['sg']} {sir['ml_label']!r} identity character "
                "disagrees with record dimension"
            )

    # ── Verify full scalar identity characters against authoritative dims ──
    for i, (sg_num, ml_label) in enumerate(zip(sg, ml)):
        op_count = char_counts[i]
        pir_rot_start = pir_rot_starts[i]
        pir_trans_start = pir_trans_starts[i]
        identities = []
        for op in range(op_count):
            rotation = pir_rots_flat[pir_rot_start + op * 9:pir_rot_start + (op + 1) * 9]
            translation = pir_trans_flat[pir_trans_start + op * 3:pir_trans_start + (op + 1) * 3]
            if rotation == identity_rotation and all(
                    abs(value - round(value)) <= 1e-8 for value in translation):
                identities.append(op)
        if len(identities) != 1:
            raise ValueError(
                f"scalar SG{sg_num} {ml_label!r} has "
                f"{len(identities)} final Seitz identities"
            )
        char_index = char_starts[i] + identities[0]
        if char_index >= len(chars_flat):
            raise ValueError(
                f"scalar SG{sg_num} {ml_label!r} has out-of-bounds identity character"
            )
        value = chars_flat[char_index]
        expected = scalar_full_dims[i]
        if (not math.isfinite(value) or value <= 0.0
                or abs(value - round(value)) > 1e-8
                or abs(value - expected) > 1e-8):
            raise ValueError(
                f"scalar SG{sg_num} {ml_label!r} identity character {value} "
                f"does not equal authoritative full dimension {expected}"
            )
    print(f"  ✓ All {len(ml)} scalar identity characters match authoritative dimensions")

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
        return _format_rust_f64(v)

    for chunk_start in range(0, len(chars_flat), 10):
        chunk = chars_flat[chunk_start:chunk_start + 10]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    lines.append("/// Complex characters of the first/stored k-star arm for scalar PIRs.")
    lines.append("/// Indices are PIR operation indices (`_pir_rot_start / 9`), not CHARACTERS indices.")
    lines.append(f"pub static SCALAR_LITTLE_CHARS_REAL: [f64; {len(little_chars_real)}] = [")
    for chunk_start in range(0, len(little_chars_real), 10):
        chunk = little_chars_real[chunk_start:chunk_start + 10]
        lines.append(f"    {', '.join(_fmt_char(v) for v in chunk)},")
    lines.append("];")
    lines.append(f"pub static SCALAR_LITTLE_CHARS_IMAG: [f64; {len(little_chars_imag)}] = [")
    for chunk_start in range(0, len(little_chars_imag), 10):
        chunk = little_chars_imag[chunk_start:chunk_start + 10]
        lines.append(f"    {', '.join(_fmt_char(v) for v in chunk)},")
    lines.append("];")
    lines.append(f"pub static SCALAR_LITTLE_CHARS_VALID: [u8; {len(little_chars_valid)}] = [")
    for chunk_start in range(0, len(little_chars_valid), 20):
        chunk = little_chars_valid[chunk_start:chunk_start + 20]
        lines.append(f"    {', '.join(str(v) for v in chunk)},")
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
    lines.append("/// Rotation matrices for final-Hall PIR operation metadata, 9 i32 per op.")
    lines.append("/// Used to build H_ops → PIR index mappings; not a phase-aligned pair with CHARACTERS.")
    lines.append(f"pub static PIR_ROTS: [i32; {len(pir_rots_flat)}] = [")
    for chunk_start in range(0, len(pir_rots_flat), 9):
        chunk = pir_rots_flat[chunk_start:chunk_start + 9]
        vals = ", ".join(str(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── PIR translation vectors ──
    lines.append("/// Translation vectors for final-Hall PIR operation metadata, 3 f64 per op.")
    lines.append("/// Used with PIR_ROTS for Seitz operation mapping; not phase-aligned with CHARACTERS.")
    lines.append(f"pub static PIR_TRANS: [f64; {len(pir_trans_flat)}] = [")
    for chunk_start in range(0, len(pir_trans_flat), 3):
        chunk = pir_trans_flat[chunk_start:chunk_start + 3]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── CIR component complex characters ──
    lines.append("/// Selected-arm complex character tables for CIR components of compound irreps.")
    lines.append("/// Stored as (re, im) pairs in data-Hall order and used by the Wigner test.")
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

    # ── Imaginary parts of spinor characters ──
    lines.append("/// Imaginary parts of spinor irrep characters.")
    lines.append("/// Indexed by IrrepRecord._spin_imag_start / _spin_imag_count.")
    lines.append(f"pub static SPIN_IMAG_CHARS: [f64; {len(spin_extra_flat)}] = [")
    for chunk_start in range(0, len(spin_extra_flat), 10):
        chunk = spin_extra_flat[chunk_start:chunk_start + 10]
        vals = ", ".join(_fmt_char(v) for v in chunk)
        lines.append(f"    {vals},")
    lines.append("];")
    lines.append("")

    # ── Spinor little-group operation indices ──
    lines.append("/// Per-irrep SG-local little-group operation indices into SPIN_OP_* arrays.")
    lines.append("/// Maps local character table position → SG-local spin operation index.")
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

    # ── SG_DATA_HALL: canonical Hall number per SG ──
    lines.append("/// Hall number identifying the final-Hall operation metadata for each SG (1-230).")
    lines.append("/// Use `SymmetryOps::from_database(SG_DATA_HALL[sg])` to get mapping H_ops.")
    lines.append("/// Legacy CHARACTERS may retain source Seitz representatives and are not phase-aligned by this table.")
    lines.append(f"pub static SG_DATA_HALL: [u16; 231] = [")
    lines.append("    0,  // dummy for index 0")
    for s in range(1, 231):
        hall = sg_hall_choice.get(s, (0, None, None))[0]
        lines.append(f"    {hall},  // SG {s}")
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
    ordinary_scalar_count = 0
    image_dimension_count = 0
    pir_dimension_count = 0
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

        # Dimensions are authoritative source metadata: exact CIR headers for
        # ordinary rows, or the resolved CIR constituent dimensions assembled
        # according to the frozen compound semantics. Character values are
        # checked against this value below, never used to define it.
        pir_d = pir_dim_map.get((sg_num, ml_label))
        dim = scalar_full_dims[i]
        if dim <= 0:
            raise ValueError(
                f"authoritative dimension is invalid for SG{sg_num} {ml_label!r}"
            )
        if not (1 <= img_code <= len(img_dims)):
            raise ValueError(
                f"missing authoritative image dimension for SG{sg_num} {ml_label!r}: "
                f"image code {img_code}"
            )
        image_dim = img_dims[img_code - 1]
        if image_dim != dim:
            raise ValueError(
                f"image/CIR authoritative dimension mismatch for SG{sg_num} "
                f"{ml_label!r}: image={image_dim}, CIR={dim}"
            )
        image_dimension_count += 1
        if pir_d is not None and pir_d != dim:
            raise ValueError(
                f"PIR/CIR authoritative dimension mismatch for SG{sg_num} "
                f"{ml_label!r}: PIR={pir_d}, CIR={dim}"
            )
        if pir_d is not None:
            pir_dimension_count += 1

        compound = compound_resolutions[i] is not None
        selected_dim = scalar_selected_dims[i]
        if compound:
            if selected_dim != 0:
                raise ValueError(
                    f"compound SG{sg_num} {ml_label!r} has nonzero scalar selected dimension"
                )
        else:
            ordinary_scalar_count += 1
            cir_entry = cir_data.get((sg_num, ml_label))
            if cir_entry is None or dim != cir_entry["dim"]:
                raise ValueError(
                    f"ordinary scalar SG{sg_num} {ml_label!r} PIR/CIR full dimension mismatch"
                )
            if dim != selected_dim * cir_entry["star_count"]:
                raise ValueError(
                    f"ordinary scalar SG{sg_num} {ml_label!r} dimension mismatch: "
                    f"record={dim}, CIR={cir_entry['dim']}, selected×star="
                    f"{selected_dim}×{cir_entry['star_count']}"
                )

        mat_s = mat_starts[i]
        mat_c = mat_counts[i]
        scalar_records.append({
            "sg": sg_num, "ml": ml_label, "bc": latex_bc, "kov": latex_kov,
            "dim": dim, "img": img_name, "lifshitz": lif_val == 1,
            "spinor": False, "kx": kx, "ky": ky, "kz": kz, "kd": kd,
            "scalar_selected_dim": selected_dim,
            "char_s": char_s, "char_c": char_c,
            "mat_s": mat_s, "mat_c": mat_c,
            "iso_s": iso_s, "iso_c": iso_c,
            "mag_iso_s": mag_iso_s, "mag_iso_c": mag_iso_c,
            "cir_s": cir_comp_starts[i], "cir_c": cir_comp_counts[i], "cir_o": cir_comp_ops[i],
            "compound_metadata_index": compound_metadata_indices[i],
            "pir_rot_s": pir_rot_starts[i],
            "spin_lg_count": 0,
            "spin_lg_op_s": 0,
            "spin_lg_op_c": 0,
            "spin_extra_s": 0,
            "spin_extra_c": 0,
        })

    if ordinary_scalar_count != 4105:
        raise ValueError(
            f"ordinary scalar CIR selected-dimension census expected 4105, got "
            f"{ordinary_scalar_count}"
        )
    if image_dimension_count != len(ml):
        raise ValueError(
            f"image dimension census expected {len(ml)}/{len(ml)}, got "
            f"{image_dimension_count}/{len(ml)}"
        )
    if pir_dimension_count != 4665:
        raise ValueError(
            f"PIR dimension census expected 4665/4777, got "
            f"{pir_dimension_count}/4777"
        )
    print(
        f"  Authoritative scalar dimensions: CIR {len(ml)}/4777, "
        f"PIR cross-check {pir_dimension_count}/4777, "
        f"image cross-check {image_dimension_count}/4777"
    )

    # Pre-compute spinor IrrepRecord data
    spinor_records = []
    for idx, sir in enumerate(spinor_irreps):
        latex_bc = label_to_latex(sir["ml_label"])
        spinor_records.append({
            "sg": sir["sg"], "ml": sir["ml_label"], "bc": latex_bc, "kov": "",
            "dim": sir["dim"], "img": "?", "lifshitz": False,
            "spinor": True,
            "scalar_selected_dim": 0,
            "kx": sir["kx"], "ky": sir["ky"], "kz": sir["kz"], "kd": sir["kd"],
            "char_s": spinor_starts[idx], "char_c": spinor_counts[idx],
            "mat_s": 0, "mat_c": 0,
            "iso_s": 0, "iso_c": 0,
            "mag_iso_s": 0, "mag_iso_c": 0,
            "cir_s": 0, "cir_c": 0, "cir_o": 0,
            "compound_metadata_index": 0,
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
            lines.append(f"        _scalar_selected_dim: {r['scalar_selected_dim']},")
            lines.append(f"        _mat_start: {r['mat_s']},")
            lines.append(f"        _mat_count: {r['mat_c']},")
            lines.append(f"        _iso_start: {r['iso_s']},")
            lines.append(f"        _iso_count: {r['iso_c']},")
            lines.append(f"        _mag_iso_start: {r['mag_iso_s']},")
            lines.append(f"        _mag_iso_count: {r['mag_iso_c']},")
            lines.append(f"        _cir_start: {r['cir_s']},")
            lines.append(f"        _cir_count: {r['cir_c']},")
            lines.append(f"        _cir_ops: {r['cir_o']},")
            lines.append(f"        _compound_metadata_index: {r['compound_metadata_index']},")
            lines.append(f"        _spin_lg_count: {r['spin_lg_count']},")
            lines.append(f"        _spin_lg_op_start: {r['spin_lg_op_s']},")
            lines.append(f"        _spin_lg_op_count: {r['spin_lg_op_c']},")
            lines.append(f"        _spin_imag_start: {r['spin_extra_s']},")
            lines.append(f"        _spin_imag_count: {r['spin_extra_c']},")
            lines.append(f"    }},")
    lines.append("];")
    lines.append("")

    lines.append("/// Frozen generation-time metadata for accepted compound records.")
    lines.append(f"pub static COMPOUND_METADATA: [CompoundMetadata; {len(compound_metadata)}] = [")
    for metadata in compound_metadata:
        cir_labels = ", ".join(
            f'"{escape_rust_str(label)}"' for label in metadata["cir_labels"]
        )
        cir_ids = ", ".join(str(value) for value in metadata["cir_irnumbers"])
        cir_dims = ", ".join(str(value) for value in metadata["cir_dimensions"])
        semantics = (
            "CompoundCharacterSemantics::ConjugateRealification"
            if metadata["semantics"] == 1
            else "CompoundCharacterSemantics::DistinctComponentSum"
        )
        lines.append("    CompoundMetadata {")
        lines.append(f"        sg: {metadata['sg']},")
        lines.append(f'        record_label: "{escape_rust_str(metadata["record_label"])}",')
        lines.append(f"        cir_irnumbers: [{cir_ids}],")
        lines.append(f"        cir_labels: [{cir_labels}],")
        lines.append(f"        cir_dimensions: [{cir_dims}],")
        lines.append(
            "        naming_grammar_version: "
            "COMPOUND_NAMING_GRAMMAR_VERSION,"
        )
        lines.append("        provenance: COMPOUND_NAMING_PROVENANCE,")
        lines.append(f"        semantics: {semantics},")
        lines.append("    },")
    lines.append("];\n")

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
