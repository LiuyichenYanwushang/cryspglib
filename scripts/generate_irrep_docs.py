#!/usr/bin/env python3
"""Generate rustdoc tables for tetragonal, trigonal, and hexagonal space groups.

Reads irrep data from src/irrep/generated_data.rs (produced by
scripts/generate_irrep_data.py) and writes three files:

  src/irrep/tetragonal.rs  – SG  75–142
  src/irrep/trigonal.rs    – SG 143–167
  src/irrep/hexagonal.rs   – SG 168–194

Format matches src/irrep/cubic.rs exactly.
"""

import os, re, sys

DATA = os.path.join(os.path.dirname(__file__), '..', 'src', 'irrep', 'generated_data.rs')
SRC_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'irrep')

with open(DATA, 'r') as f:
    content = f.read()

# ── Parse SG_IRREP_INDEX ──────────────────────────────────────────────────

sg_idx_start = content.index("pub static SG_IRREP_INDEX: [(u16, u16); 231] = [")
sg_idx_end = content.index("];\n\n", sg_idx_start)
sg_idx_text = content[sg_idx_start:sg_idx_end]
sg_irrep_index = {}
pairs = re.findall(r'\((\d+),\s*(\d+)\)', sg_idx_text)
for i, (start, count) in enumerate(pairs):
    sg_irrep_index[i] = (int(start), int(count))
print(f"SG_IRREP_INDEX entries: {len(sg_irrep_index)}")

# ── Parse IRREPS ──────────────────────────────────────────────────────────

def unescape(s):
    return s.replace('\\\\', '\\')

irreps_start = content.index("pub static IRREPS: [IrrepRecord; 4777] = [")
irreps_start += len("pub static IRREPS: [IrrepRecord; 4777] = [")
irreps_end = content.index("\n];", irreps_start)
irreps_text = content[irreps_start:irreps_end]

irreps = []
pattern = (
    r'IrrepRecord\s*\{\s*'
    r'sg:\s*(\d+),\s*'
    r'ml:\s*"([^"]*)",\s*'
    r'bc:\s*"([^"]*)",\s*'
    r'kov:\s*"([^"]*)",\s*'
    r'dim:\s*(\d+),\s*'
    r'image:\s*"([^"]*)",\s*'
    r'lifshitz:\s*(true|false),\s*'
    r'\},'
)
for m in re.finditer(pattern, irreps_text):
    irreps.append({
        'sg': int(m.group(1)),
        'ml': m.group(2),
        'bc': unescape(m.group(3)),
        'kov': unescape(m.group(4)),
        'dim': int(m.group(5)),
        'image': m.group(6),
        'lifshitz': m.group(7) == 'true',
    })
print(f"IRREPS parsed: {len(irreps)}")

# ── Parse ISOTROPY_SUBGROUPS ──────────────────────────────────────────────

iso_start = content.index("pub static ISOTROPY_SUBGROUPS: [IsotropyRecord; 15239] = [")
iso_start += len("pub static ISOTROPY_SUBGROUPS: [IsotropyRecord; 15239] = [")
iso_end = content.index("\n];", iso_start)
iso_text = content[iso_start:iso_end]

isos = []
iso_pattern = (
    r'IsotropyRecord\s*\{\s*'
    r'sg:\s*(\d+),\s*'
    r'symbol:\s*"([^"]*)",\s*'
    r'schoenflies:\s*"([^"]*)",\s*'
    r'direction:\s*"([^"]*)",\s*'
    r'domains:\s*(\d+),\s*'
    r'arms:\s*(\d+),\s*'
    r'\},'
)
for m in re.finditer(iso_pattern, iso_text):
    isos.append({
        'sg': int(m.group(1)),
        'symbol': unescape(m.group(2)),
        'schoenflies': unescape(m.group(3)),
        'direction': m.group(4),
        'domains': int(m.group(5)),
        'arms': int(m.group(6)),
    })
print(f"ISOTROPY_SUBGROUPS parsed: {len(isos)}")

# ── Parse IRREP_ISOTROPY_INDEX ────────────────────────────────────────────

iri_start = content.index("pub static IRREP_ISOTROPY_INDEX: [(u16, u16); 4777] = [")
iri_start += len("pub static IRREP_ISOTROPY_INDEX: [(u16, u16); 4777] = [")
iri_end = content.index("\n];", iri_start)
iri_text = content[iri_start:iri_end]
irrep_iso_index = []
for m in re.finditer(r'\((\d+),\s*(\d+)\)', iri_text):
    irrep_iso_index.append((int(m.group(1)), int(m.group(2))))
print(f"IRREP_ISOTROPY_INDEX entries: {len(irrep_iso_index)}")

# ── Common helpers ──────────────────────────────────────────────────────────

KPOINT_ORDER = ['GM', 'X', 'M', 'R', 'A', 'H', 'K', 'L', 'T', 'W', 'Z',
                'S', 'D', 'F', 'U', 'V', 'Y', 'N', 'P', 'Q', 'B', 'C',
                'E', 'G', 'I', 'O', 'J']

def kpoint_sort_key(prefix):
    upper = prefix.upper()
    try:
        return (0, KPOINT_ORDER.index(upper), '')
    except ValueError:
        return (1, 0, upper)

def extract_kpoint_prefix(ml_label):
    """Extract k-point prefix from ML label (first contiguous alphabetic run)."""
    prefix = []
    for c in ml_label:
        if c.isalpha():
            prefix.append(c)
        else:
            break
    return ''.join(prefix)

def sort_irreps(irrep_list):
    def key(r):
        return (kpoint_sort_key(r['kpoint_prefix']), r['ml'])
    return sorted(irrep_list, key=key)

def bc_latex(bc):
    return f"${bc}$"

def kov_latex(kov):
    return f"${kov}$"

# ── SG name tables (Hermann-Mauguin + Schoenflies) ────────────────────────

TETRAGONAL_NAMES = {
    75:  ("P4",       "C₄¹"),
    76:  ("P4₁",      "C₄²"),
    77:  ("P4₂",      "C₄³"),
    78:  ("P4₃",      "C₄⁴"),
    79:  ("I4",       "C₄⁵"),
    80:  ("I4₁",      "C₄⁶"),
    81:  ("P-4",      "S₄¹"),
    82:  ("I-4",      "S₄²"),
    83:  ("P4/m",     "C₄ₕ¹"),
    84:  ("P4₂/m",    "C₄ₕ²"),
    85:  ("P4/n",     "C₄ₕ³"),
    86:  ("P4₂/n",    "C₄ₕ⁴"),
    87:  ("I4/m",     "C₄ₕ⁵"),
    88:  ("I4₁/a",    "C₄ₕ⁶"),
    89:  ("P422",     "D₄¹"),
    90:  ("P42₁2",    "D₄²"),
    91:  ("P4₁22",    "D₄³"),
    92:  ("P4₁2₁2",   "D₄⁴"),
    93:  ("P4₂22",    "D₄⁵"),
    94:  ("P4₂2₁2",   "D₄⁶"),
    95:  ("P4₃22",    "D₄⁷"),
    96:  ("P4₃2₁2",   "D₄⁸"),
    97:  ("I422",     "D₄⁹"),
    98:  ("I4₁22",    "D₄¹⁰"),
    99:  ("P4mm",     "C₄ᵥ¹"),
    100: ("P4bm",     "C₄ᵥ²"),
    101: ("P4₂cm",    "C₄ᵥ³"),
    102: ("P4₂nm",    "C₄ᵥ⁴"),
    103: ("P4cc",     "C₄ᵥ⁵"),
    104: ("P4nc",     "C₄ᵥ⁶"),
    105: ("P4₂mc",    "C₄ᵥ⁷"),
    106: ("P4₂bc",    "C₄ᵥ⁸"),
    107: ("I4mm",     "C₄ᵥ⁹"),
    108: ("I4cm",     "C₄ᵥ¹⁰"),
    109: ("I4₁md",    "C₄ᵥ¹¹"),
    110: ("I4₁cd",    "C₄ᵥ¹²"),
    111: ("P-42m",    "D₂ₔ¹"),
    112: ("P-42c",    "D₂ₔ²"),
    113: ("P-42₁m",   "D₂ₔ³"),
    114: ("P-42₁c",   "D₂ₔ⁴"),
    115: ("P-4m2",    "D₂ₔ⁵"),
    116: ("P-4c2",    "D₂ₔ⁶"),
    117: ("P-4b2",    "D₂ₔ⁷"),
    118: ("P-4n2",    "D₂ₔ⁸"),
    119: ("I-4m2",    "D₂ₔ⁹"),
    120: ("I-4c2",    "D₂ₔ¹⁰"),
    121: ("I-42m",    "D₂ₔ¹¹"),
    122: ("I-42d",    "D₂ₔ¹²"),
    123: ("P4/mmm",   "D₄ₕ¹"),
    124: ("P4/mcc",   "D₄ₕ²"),
    125: ("P4/nbm",   "D₄ₕ³"),
    126: ("P4/nnc",   "D₄ₕ⁴"),
    127: ("P4/mbm",   "D₄ₕ⁵"),
    128: ("P4/mnc",   "D₄ₕ⁶"),
    129: ("P4/nmm",   "D₄ₕ⁷"),
    130: ("P4/ncc",   "D₄ₕ⁸"),
    131: ("P4₂/mmc",  "D₄ₕ⁹"),
    132: ("P4₂/mcm",  "D₄ₕ¹⁰"),
    133: ("P4₂/nbc",  "D₄ₕ¹¹"),
    134: ("P4₂/nnm",  "D₄ₕ¹²"),
    135: ("P4₂/mbc",  "D₄ₕ¹³"),
    136: ("P4₂/mnm",  "D₄ₕ¹⁴"),
    137: ("P4₂/nmc",  "D₄ₕ¹⁵"),
    138: ("P4₂/ncm",  "D₄ₕ¹⁶"),
    139: ("I4/mmm",   "D₄ₕ¹⁷"),
    140: ("I4/mcm",   "D₄ₕ¹⁸"),
    141: ("I4₁/amd",  "D₄ₕ¹⁹"),
    142: ("I4₁/acd",  "D₄ₕ²⁰"),
}

TRIGONAL_NAMES = {
    143: ("P3",       "C₃¹"),
    144: ("P3₁",      "C₃²"),
    145: ("P3₂",      "C₃³"),
    146: ("R3",       "C₃⁴"),
    147: ("P-3",      "C₃ᵢ¹"),
    148: ("R-3",      "C₃ᵢ²"),
    149: ("P312",     "D₃¹"),
    150: ("P321",     "D₃²"),
    151: ("P3₁12",    "D₃³"),
    152: ("P3₁21",    "D₃⁴"),
    153: ("P3₂12",    "D₃⁵"),
    154: ("P3₂21",    "D₃⁶"),
    155: ("R32",      "D₃⁷"),
    156: ("P3m1",     "C₃ᵥ¹"),
    157: ("P31m",     "C₃ᵥ²"),
    158: ("P3c1",     "C₃ᵥ³"),
    159: ("P31c",     "C₃ᵥ⁴"),
    160: ("R3m",      "C₃ᵥ⁵"),
    161: ("R3c",      "C₃ᵥ⁶"),
    162: ("P-31m",    "D₃ₔ¹"),
    163: ("P-31c",    "D₃ₔ²"),
    164: ("P-3m1",    "D₃ₔ³"),
    165: ("P-3c1",    "D₃ₔ⁴"),
    166: ("R-3m",     "D₃ₔ⁵"),
    167: ("R-3c",     "D₃ₔ⁶"),
}

HEXAGONAL_NAMES = {
    168: ("P6",       "C₆¹"),
    169: ("P6₁",      "C₆²"),
    170: ("P6₅",      "C₆³"),
    171: ("P6₂",      "C₆⁴"),
    172: ("P6₄",      "C₆⁵"),
    173: ("P6₃",      "C₆⁶"),
    174: ("P-6",      "C₃ₕ¹"),
    175: ("P6/m",     "C₆ₕ¹"),
    176: ("P6₃/m",    "C₆ₕ²"),
    177: ("P622",     "D₆¹"),
    178: ("P6₁22",    "D₆²"),
    179: ("P6₅22",    "D₆³"),
    180: ("P6₂22",    "D₆⁴"),
    181: ("P6₄22",    "D₆⁵"),
    182: ("P6₃22",    "D₆⁶"),
    183: ("P6mm",     "C₆ᵥ¹"),
    184: ("P6cc",     "C₆ᵥ²"),
    185: ("P6₃cm",    "C₆ᵥ³"),
    186: ("P6₃mc",    "C₆ᵥ⁴"),
    187: ("P-6m2",    "D₃ₕ¹"),
    188: ("P-6c2",    "D₃ₕ²"),
    189: ("P-62m",    "D₃ₕ³"),
    190: ("P-62c",    "D₃ₕ⁴"),
    191: ("P6/mmm",   "D₆ₕ¹"),
    192: ("P6/mcc",   "D₆ₕ²"),
    193: ("P6₃/mcm",  "D₆ₕ³"),
    194: ("P6₃/mmc",  "D₆ₕ⁴"),
}

# ── k-point descriptions per crystal system ───────────────────────────────

TETRAGONAL_KPOINTS = {
    'GM': r'\Gamma\ (0,0,0)',
    'X':  r'\mathrm{X}\ (1/2,0,0)',
    'M':  r'\mathrm{M}\ (1/2,1/2,0)',
    'Z':  r'\mathrm{Z}\ (0,0,1/2)',
    'R':  r'\mathrm{R}\ (0,1/2,1/2)',
    'A':  r'\mathrm{A}\ (1/2,1/2,1/2)',
}

TRIGONAL_KPOINTS = {
    'GM': r'\Gamma\ (0,0,0)',
    'A':  r'\mathrm{A}\ (0,0,1/2)',
    'K':  r'\mathrm{K}\ (1/3,1/3,0)',
    'H':  r'\mathrm{H}\ (1/3,1/3,1/2)',
    'M':  r'\mathrm{M}\ (1/2,0,0)',
    'L':  r'\mathrm{L}\ (1/2,0,1/2)',
}

HEXAGONAL_KPOINTS = {
    'GM': r'\Gamma\ (0,0,0)',
    'A':  r'\mathrm{A}\ (0,0,1/2)',
    'K':  r'\mathrm{K}\ (1/3,1/3,0)',
    'H':  r'\mathrm{H}\ (1/3,1/3,1/2)',
    'M':  r'\mathrm{M}\ (1/2,0,0)',
    'L':  r'\mathrm{L}\ (1/2,0,1/2)',
}

# ── Module-level doc comments (preserved from current placeholders) ───────

TETRAGONAL_HEADER = """//! # Tetragonal space groups (#75–#142)
//!
//! Tetragonal crystals have a = b ≠ c, α = β = γ = 90°.  The Brillouin
//! zone is a tetragonal prism.
//!
//! ## Common k-point labels (simple tetragonal)
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | 4/mmm (D₄ₕ) |
//! | X | (½, 0, 0) | 4/mmm (D₄ₕ) |
//! | M | (½, ½, 0) | 4/mmm (D₄ₕ) |
//! | Z | (0, 0, ½) | 4/mmm (D₄ₕ) |
//! | R | (0, ½, ½) | 4/mmm (D₄ₕ) |
//! | A | (½, ½, ½) | 4/mmm (D₄ₕ) |
//!
//! ---"""

TRIGONAL_HEADER = """//! # Trigonal space groups (#143–#167)
//!
//! Trigonal space groups use hexagonal axes (a = b ≠ c, γ = 120°).
//! The Brillouin zone is a hexagonal prism.
//!
//! ## Common k-point labels (hexagonal setting)
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | -3m (D₃ₔ) |
//! | A | (0, 0, ½) | -3m (D₃ₔ) |
//! | K | (⅓, ⅓, 0) | 3m (C₃ᵥ) |
//! | H | (⅓, ⅓, ½) | 3m (C₃ᵥ) |
//! | M | (½, 0, 0) | 2/m (C₂ₕ) |
//! | L | (½, 0, ½) | 2/m (C₂ₕ) |
//!
//! ---"""

HEXAGONAL_HEADER = """//! # Hexagonal space groups (#168–#194)
//!
//! Hexagonal crystals have a = b ≠ c, α = β = 90°, γ = 120°.
//! The Brillouin zone is a hexagonal prism.
//!
//! ## Common k-point labels
//!
//! | Label | Coords (fractional) | Little group |
//! |-------|---------------------|--------------|
//! | Γ | (0, 0, 0) | 6/mmm (D₆ₕ) |
//! | A | (0, 0, ½) | 6/mmm (D₆ₕ) |
//! | K | (⅓, ⅓, 0) | -6m2 (D₃ₕ) |
//! | H | (⅓, ⅓, ½) | -6m2 (D₃ₕ) |
//! | M | (½, 0, 0) | mmm (D₂ₕ) |
//! | L | (½, 0, ½) | mmm (D₂ₕ) |
//!
//! ---"""

# ── SG rustdoc generator ──────────────────────────────────────────────────

def generate_sg_rustdoc(sg_num, sg_names, kpoint_descriptions):
    """Generate rustdoc block for one space group."""
    hm, sch = sg_names[sg_num]
    start, count = sg_irrep_index[sg_num]

    lines = []
    lines.append(f'/// # {sg_num} {hm} ({sch})')
    lines.append('///')

    # Collect irreps for this SG
    sg_irreps = []
    for idx in range(start, start + count):
        r = dict(irreps[idx])
        r['_idx'] = idx
        sg_irreps.append(r)

    for r in sg_irreps:
        r['kpoint_prefix'] = extract_kpoint_prefix(r['ml'])

    # Group by k-point prefix
    kpoint_groups = {}
    for r in sg_irreps:
        prefix = r['kpoint_prefix']
        kpoint_groups.setdefault(prefix, []).append(r)

    sorted_prefixes = sorted(kpoint_groups.keys(), key=kpoint_sort_key)

    for prefix in sorted_prefixes:
        group = sort_irreps(kpoint_groups[prefix])
        kp_desc = kpoint_descriptions.get(prefix, prefix)
        lines.append(f'/// ## Irreps at ${kp_desc}$')
        lines.append('///')
        lines.append('/// | ML | B&C | Kovalev | Dim | Image | Lifshitz |')
        lines.append('/// |-----|-----|---------|-----|-------|----------|')
        for r in group:
            ml = r['ml']
            bc = bc_latex(r['bc'])
            kov = kov_latex(r['kov'])
            dim = r['dim']
            image = r['image']
            lif = 'yes' if r['lifshitz'] else 'no'
            lines.append(f'/// | {ml} | {bc} | {kov} | {dim} | {image} | {lif} |')
        lines.append('///')

        # Isotropy subgroup tables for each irrep
        for r in group:
            iso_s, iso_c = irrep_iso_index[r['_idx']]
            if iso_c == 0:
                continue
            lines.append(
                f'/// ### Isotropy subgroups for {r["ml"]} '
                f'({bc_latex(r["bc"])})')
            lines.append('///')
            lines.append('/// | Direction | SG | Symbol | Schoenflies | Domains | Arms |')
            lines.append('/// |-----------|-----|--------|-------------|---------|------|')
            for j in range(iso_s, iso_s + iso_c):
                iso = isos[j]
                direction = f"dir{{{iso['direction']}}}"
                sg = iso['sg']
                symbol = iso['symbol']
                schoenflies = iso['schoenflies']
                domains = iso['domains']
                arms = iso['arms']
                lines.append(
                    f'/// | {direction} | {sg} | {symbol} | {schoenflies} | {domains} | {arms} |')
            lines.append('///')

    lines.append(f'pub struct Sg{sg_num};')
    return '\n'.join(lines)


def write_system_file(filepath, header, sg_range, sg_names, kpoint_descriptions):
    """Write a complete crystal-system module file."""
    with open(filepath, 'w') as f:
        f.write(header)
        f.write('\n')
        for sg_num in sg_range:
            f.write('\n')
            f.write(generate_sg_rustdoc(sg_num, sg_names, kpoint_descriptions))
            f.write('\n')
    print(f"  Wrote {filepath}")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print("Generating rustdoc tables...\n")

    write_system_file(
        os.path.join(SRC_DIR, 'tetragonal.rs'),
        TETRAGONAL_HEADER,
        range(75, 143),
        TETRAGONAL_NAMES,
        TETRAGONAL_KPOINTS,
    )

    write_system_file(
        os.path.join(SRC_DIR, 'trigonal.rs'),
        TRIGONAL_HEADER,
        range(143, 168),
        TRIGONAL_NAMES,
        TRIGONAL_KPOINTS,
    )

    write_system_file(
        os.path.join(SRC_DIR, 'hexagonal.rs'),
        HEXAGONAL_HEADER,
        range(168, 195),
        HEXAGONAL_NAMES,
        HEXAGONAL_KPOINTS,
    )

    print("\nDone!")


if __name__ == '__main__':
    main()
