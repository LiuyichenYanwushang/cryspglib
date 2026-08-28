#!/usr/bin/env python3
"""Parse double-valued (spinor) irrep data from irrepTables package.

The irrepTables package (pip install irreptables) contains data files
sourced from the Bilbao Crystallographic Server. Each file covers one
space group:

  ~/.local/lib/python*/site-packages/irreptables/tables/irreps-SG=*-spin.dat

Output format:
- SG#, k-point label (GM, X, ...), coords, irrep label, dim,
  complex character table split into real and imaginary parts

## Bilbao SU(2) storage convention (verified by scripts/test_su2_closure.py)

Each symmetry operation line in spin.dat stores 20 values:

    rot[9] trans[3] amp[4] phase[4]

where:
  amp[4]   = |U_ij| — amplitudes of the 4 complex 2×2 matrix elements
  phase[4] = arg(U_ij)/π — phases in units of π

The complex SU(2) matrix is reconstructed as:

    U_ij = amp[ij] · exp(iπ · phase[ij])

Then decomposed into real Pauli coefficients (u₀, u₁, u₂, u₃):

    U = u₀·I + i(u₁·σx + u₂·σy + u₃·σz)

      = [[u₀ + iu₃,    u₂ + iu₁],
         [-u₂ + iu₁,    u₀ - iu₃]]

For crystallographic point-group rotations, the Pauli coefficients
take values only from the set {0, ±½, ±1/√2, ±√3/2, ±1}.
These are stored as f64 rounded to the nearest exact algebraic value.
"""
import base64
import hashlib
import math, os, re, sys, glob
from collections import defaultdict
from fractions import Fraction

IRREPTABLES_VERSION = "1.0.0"
IRREPTABLES_RECORD_SHA256 = (
    "726a945eb60ffeae968b8196e125494bdf870870f693c5a32c70b7761aa091bb"
)
IRREPTABLES_SG3_SOURCE_SHA256 = (
    "75020a002c7006503a15c13fe89040a229989f4c025dd0e90a1ea2f68818dbfc"
)

# ── Exact Pauli coefficient rounding ────────────────────────────────────
# For crystallographic double-group operations, the Pauli coefficients
# u₀,u₁,u₂,u₃ ∈ {0, ±½, ±1/√2, ±√3/2, ±1}.
# We round floating-point values to the nearest exact algebraic number
# to eliminate numerical noise from sin/cos computations.

_SQRT2 = math.sqrt(2.0)
_SQRT3 = math.sqrt(3.0)

_EXACT_PAULI_VALUES = [
    0.0,
    0.5, -0.5,
    1.0 / _SQRT2, -1.0 / _SQRT2,
    _SQRT3 / 2.0, -_SQRT3 / 2.0,
    1.0, -1.0,
]


def _round_to_exact_pauli(val, tol=1e-10):
    """Legacy-only rounding of a materialized Pauli coefficient."""
    for exact in _EXACT_PAULI_VALUES:
        if abs(val - exact) < tol:
            return exact
    return val


def _round_amplitude(val):
    """Legacy-only rounding of an already materialized amplitude.

    Uses a relaxed tolerance because spin.dat files store amplitudes with
    only ~5 significant digits (e.g. 0.70711 for 1/√2 ≈ 0.70710678).
    """
    for exact in [0.0, 1.0 / _SQRT2, 1.0]:
        if abs(val - exact) < 1e-4:
            return exact
    return val


def _amp_phase_to_pauli(amp, phase):
    """Convert Bilbao polar encoding (amp[4] + phase[4]/π) to Pauli coefficients.

    Returns [u₀, u₁, u₂, u₃] as exact f64 values.

    The amp values in spin.dat files are rounded to ~5 significant digits
    (e.g. 0.70711 for 1/√2).  We round them to the exact algebraic value
    {0, 1/√2, 1} before multiplying, so that uᵢ = amp × cos/sin(π·phase)
    comes out exact (e.g. 0.5 instead of 0.500002276).
    """
    # Round amplitudes to exact {0, 1/√2, 1} before multiplying.
    # This eliminates the ~5-digit precision loss in spin.dat amplitude values.
    exact_amp = [_round_amplitude(a) for a in amp]

    u0 = _round_to_exact_pauli(exact_amp[0] * math.cos(math.pi * phase[0]))
    u3 = _round_to_exact_pauli(exact_amp[0] * math.sin(math.pi * phase[0]))
    u2 = _round_to_exact_pauli(exact_amp[1] * math.cos(math.pi * phase[1]))
    u1 = _round_to_exact_pauli(exact_amp[1] * math.sin(math.pi * phase[1]))

    return [u0, u1, u2, u3]


def _rationalize_kvector(coords):
    """Convert decimal k coordinates to the smallest supported rational tuple.

    spin.dat stores thirds and sixths with six decimal places (for example
    0.333333), so their scaled roundoff can be slightly larger than 1e-6.
    """
    from math import gcd

    tolerance = 5e-6
    for kd in [1, 2, 3, 4, 6]:
        numerators = [int(round(v * kd)) for v in coords]
        if all(abs(v * kd - n) <= tolerance for v, n in zip(coords, numerators)):
            g = kd
            for n in numerators:
                g = gcd(g, abs(n))
            if g > 1:
                numerators = [n // g for n in numerators]
                kd //= g
            return numerators[0], numerators[1], numerators[2], kd

    raise ValueError(f"Unsupported spin.dat k-vector coordinates: {coords!r}")


def _package_roots():
    roots = list(sys.path)
    import site
    roots.extend(site.getusersitepackages().split(os.pathsep))
    roots.extend(site.getsitepackages())
    return [root for root in dict.fromkeys(roots) if root]


def _verified_irreptables_distribution():
    """Return the sole verified irreptables package root and RECORD manifest."""
    candidates = []
    for root in _package_roots():
        for dist_info in glob.glob(os.path.join(root, "irreptables-*.dist-info")):
            metadata_path = os.path.join(dist_info, "METADATA")
            record_path = os.path.join(dist_info, "RECORD")
            if not (os.path.isfile(metadata_path) and os.path.isfile(record_path)):
                continue
            with open(metadata_path, encoding="utf-8") as stream:
                metadata = stream.read()
            version = next(
                (line.split(":", 1)[1].strip()
                 for line in metadata.splitlines()
                 if line.startswith("Version:")),
                None,
            )
            if version == IRREPTABLES_VERSION:
                candidates.append((os.path.dirname(dist_info), record_path))
    if len(candidates) != 1:
        raise FileNotFoundError(
            "expected exactly one irreptables==1.0.0 distribution with RECORD, "
            f"found {len(candidates)}"
        )

    package_root, record_path = candidates[0]
    with open(record_path, "rb") as stream:
        record_bytes = stream.read()
    actual_record_hash = hashlib.sha256(record_bytes).hexdigest()
    if actual_record_hash != IRREPTABLES_RECORD_SHA256:
        raise ValueError(
            "irreptables RECORD hash mismatch: "
            f"expected {IRREPTABLES_RECORD_SHA256}, got {actual_record_hash}"
        )

    manifest = {}
    for line in record_bytes.decode("utf-8").splitlines():
        fields = line.split(",")
        if len(fields) < 3 or not fields[1].startswith("sha256="):
            continue
        path, digest = fields[0], fields[1][len("sha256="):]
        manifest[path] = bytes.hex(base64.urlsafe_b64decode(digest + "=="))
    return package_root, manifest


def _verify_spin_source_files(tables_dir, package_root, manifest, expected_count=230):
    files = sorted(glob.glob(os.path.join(tables_dir, "irreps-SG=*-spin.dat")))
    if len(files) != expected_count:
        raise ValueError(
            f"expected {expected_count} pinned spin source files, found {len(files)}"
        )
    source_hashes = {}
    for filepath in files:
        relpath = os.path.relpath(filepath, package_root).replace(os.sep, "/")
        expected = manifest.get(relpath)
        if expected is None:
            raise ValueError(f"spin source file {relpath} is absent from RECORD")
        with open(filepath, "rb") as stream:
            actual = hashlib.sha256(stream.read()).hexdigest()
        if actual != expected:
            raise ValueError(
                f"spin source file hash mismatch for {relpath}: "
                f"expected {expected}, got {actual}"
            )
        source_hashes[relpath] = actual
    return tables_dir, source_hashes


def _validate_spin_source_sgs(sg_numbers):
    expected = list(range(1, 231))
    if sorted(sg_numbers) != expected:
        raise ValueError(
            "spin source files must contain each SG exactly once; "
            f"found SGs {sorted(sg_numbers)}"
        )


def _verified_spin_source_manifest():
    package_root, manifest = _verified_irreptables_distribution()
    tables_dir = os.path.join(package_root, "irreptables", "tables")
    _tables_dir, source_hashes = _verify_spin_source_files(
        tables_dir, package_root, manifest
    )
    sg3_relpath = "irreptables/tables/irreps-SG=3-spin.dat"
    if source_hashes.get(sg3_relpath) != IRREPTABLES_SG3_SOURCE_SHA256:
        raise ValueError("SG3 spin source hash does not match the pinned provenance")
    return tables_dir, source_hashes


def find_tables_dir():
    """Locate and cryptographically verify the irreptables data directory."""
    tables_dir, _source_hashes = _verified_spin_source_manifest()
    return tables_dir


def verified_spin_source_paths():
    """Return the verified raw source paths consumed by legacy and exact parsers."""
    tables_dir, _source_hashes = _verified_spin_source_manifest()
    files = tuple(sorted(glob.glob(os.path.join(tables_dir, "irreps-SG=*-spin.dat"))))
    if len(files) != 230:
        raise ValueError(f"expected 230 pinned spin source files, found {len(files)}")
    return files


def _round_char(x, eps=1e-8):
    """Round character value to clean float."""
    if abs(x) < eps:
        return 0.0
    r = round(x)
    if abs(x - r) < eps:
        return float(r)
    # spin.dat phases are stored with five decimal places (e.g. 0.66667),
    # so decoded trigonometric values need a correspondingly relaxed snap
    # to the crystallographic character field Q(sqrt(2), sqrt(3)).
    exact_candidates = set()
    for n in range(-12, 13):
        exact_candidates.add(n / 2.0)
        exact_candidates.add(n / _SQRT2)
        exact_candidates.add(n * _SQRT3 / 2.0)
    for exact in exact_candidates:
        if abs(x - exact) < 5e-5:
            return exact
    for n in range(-12, 13):
        for d in (2, 3, 4, 6, 8):
            if abs(x - n / d) < eps:
                return n / d
    return round(x, 10)


def _decode_character_polar(values, n_ops):
    """Legacy-only decode of a materialized spin.dat character row.

    Rows contain either ``n`` real characters, or ``n`` amplitudes followed
    by ``n`` phases in units of pi.  The phase half was historically
    misinterpreted as auxiliary Wigner data.
    """
    if len(values) == n_ops:
        return ([_round_char(v) for v in values], [0.0] * n_ops)
    if len(values) != 2 * n_ops:
        raise ValueError(
            f"invalid spinor character row: {len(values)} values for {n_ops} operations"
        )

    amplitudes = values[:n_ops]
    phases = values[n_ops:]
    real = [
        _round_char(a * math.cos(math.pi * p))
        for a, p in zip(amplitudes, phases)
    ]
    imag = [
        _round_char(a * math.sin(math.pi * p))
        for a, p in zip(amplitudes, phases)
    ]
    return real, imag


def parse_spinor_file(filepath):
    """Legacy f64 parse of one spin.dat file.

    Returns:
        sg: int, space group number
        spin_ops: list of dicts with keys: rot[9], trans[3], su2[4]
                  su2[4] = Pauli coefficients (u₀, u₁, u₂, u₃), exact f64
        irreps: list of dicts with keys:
            k_label, kx/ky/kz/kd, ml_label, dim, characters, op_indices
    """
    with open(filepath) as f:
        lines = f.readlines()

    sg = None
    spin_ops = []  # global symmetry operations with SU(2) lifts
    irreps = []

    # Parse header
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("SG="):
            sg = int(line.split("=")[1])
        elif line.startswith("nsym="):
            pass
        elif line.startswith("spinor="):
            pass
        elif line.startswith("symmetries="):
            i += 1
            break
        i += 1

    # Parse symmetry operations (nsym lines)
    # Format: R(3x3) 9ints | t(3) 3floats | amp(4) 4floats | phase(4)/π 4floats
    # Converted to Pauli coefficients: SU(2) = u₀I + i(u₁σx + u₂σy + u₃σz)
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("kpoint"):
            break
        parts = line.split()
        if len(parts) >= 20:
            rot = [int(x) for x in parts[0:9]]
            trans = [float(x) for x in parts[9:12]]
            amp = [float(x) for x in parts[12:16]]
            phase = [float(x) for x in parts[16:20]]
            su2 = _amp_phase_to_pauli(amp, phase)
            spin_ops.append({
                "rot": rot, "trans": trans, "su2": su2,
                "raw_rotation": parts[0:9], "raw_translation": parts[9:12],
                "raw_amp": parts[12:16], "raw_phase": parts[16:20],
            })
        elif len(parts) >= 16:
            # Fallback for files without the extra 4 columns
            rot = [int(x) for x in parts[0:9]]
            trans = [float(x) for x in parts[9:12]]
            amp = [float(x) for x in parts[12:16]]
            phase = [0.0, 0.0, 0.0, 0.0]
            su2 = _amp_phase_to_pauli(amp, phase)
            spin_ops.append({
                "rot": rot, "trans": trans, "su2": su2,
                "raw_rotation": parts[0:9], "raw_translation": parts[9:12],
                "raw_amp": parts[12:16], "raw_phase": ["0.0"] * 4,
            })
        i += 1

    # Parse kpoints and irreps
    current_k = None
    current_kvec = None
    current_k_tokens = None
    current_op_indices = None
    source_row_ordinal = 0

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("kpoint"):
            # kpoint  GM : 0.0 0.0 0.0  : 1 2 3 4 ...
            parts = line.split(":")
            k_name = parts[0].replace("kpoint", "").strip()
            coords = [float(x) for x in parts[1].strip().split()]
            current_k_tokens = parts[1].strip().split()
            op_indices = [int(x) - 1 for x in parts[2].strip().split()]  # 0-indexed
            current_k = k_name
            current_kvec = coords
            current_op_indices = op_indices
        elif line.startswith("-"):
            # -GM6 2    2.0  0.0  ...
            parts = line.split()
            if len(parts) < 3:
                i += 1
                continue
            label = parts[0][1:]  # strip leading '-'
            dim = int(parts[1])
            chars_raw = [float(x) for x in parts[2:]]
            n_ops = len(current_op_indices or [])
            chars_real, chars_imag = _decode_character_polar(chars_raw, n_ops)

            # Compute k-vector denominator from coords.
            # Use the SMALLEST common denominator so that spinor and scalar
            # irreps at the same k-point get the same (kx,ky,kz,kd) tuple
            # and are correctly grouped by kpoints_of().
            if current_kvec:
                kx_i, ky_i, kz_i, kd = _rationalize_kvector(current_kvec)
            else:
                kx_i = ky_i = kz_i = 0
                kd = 1

            irreps.append({
                "sg": sg,
                "k_label": current_k,
                "kx": kx_i, "ky": ky_i, "kz": kz_i, "kd": kd,
                "ml_label": label,
                "dim": dim,
                "characters": chars_real,
                "characters_imag": chars_imag,
                "op_indices": current_op_indices,
                "source_row_ordinal": source_row_ordinal,
                "raw_k_label": current_k,
                "raw_k_tokens": current_k_tokens,
                "raw_operation_tokens": [str(x + 1) for x in current_op_indices or []],
                "raw_dimension": parts[1],
                "raw_character_tokens": parts[2:],
            })
            source_row_ordinal += 1

        i += 1

    return sg, spin_ops, irreps


def parse_all_spinor():
    """Parse all 230 spin.dat files.

    Returns:
        all_irreps: list of spinor irrep dicts
        all_spin_ops: dict SG# -> list of spin op dicts (rot, trans, su2)
                      su2[4] = Pauli coefficients (u₀, u₁, u₂, u₃)
    """
    import spinor_exact

    tables_dir, source_hashes = _verified_spin_source_manifest()
    files = list(verified_spin_source_paths())
    # The exact sidecar owns the verified source-path lookup and caches the
    # immutable default bundle, so repeated generator passes do not rerun the
    # full algebraic validator in the same process.
    exact_files = spinor_exact.parse_all_exact()
    spinor_exact.validate_exact_sources(exact_files)
    exact_hashes = spinor_exact.source_spelling_hashes(exact_files)
    # Keep the exact sidecar's golden assertions in the production gate.
    for name, expected in spinor_exact.DOMAIN_GOLDENS.items():
        if exact_hashes[name] != expected:
            raise ValueError(f"exact {name} spelling hash mismatch")
    for name, expected in (
        ("combined", spinor_exact.COMBINED_GOLDEN),
        ("su2_pair", spinor_exact.SU2_PAIR_GOLDEN),
        ("character_pair", spinor_exact.CHARACTER_PAIR_GOLDEN),
    ):
        if exact_hashes[name] != expected:
            raise ValueError(f"exact {name} spelling hash mismatch")

    all_irreps = []
    all_spin_ops = {}  # SG# -> list of spin ops
    exact_by_sg = {source.sg: source for source in exact_files}
    for f in files:
        sg, spin_ops, irreps = parse_spinor_file(f)
        relpath = os.path.relpath(f, os.path.dirname(os.path.dirname(tables_dir))).replace(
            os.sep, "/"
        )
        source_hash = source_hashes[relpath]
        for irrep in irreps:
            irrep["source_file"] = relpath
            irrep["source_file_sha256"] = source_hash
        all_spin_ops[sg] = spin_ops
        all_irreps.extend(irreps)

        exact_source = exact_by_sg.get(sg)
        if exact_source is None or len(exact_source.operations) != len(spin_ops):
            raise ValueError(f"legacy/exact SG {sg} operation linkage mismatch")
        for exact_op, legacy_op in zip(exact_source.operations, spin_ops):
            if (tuple(exact_op.rotation) != tuple(legacy_op["rot"])
                    or exact_op.raw_rotation != tuple(legacy_op["raw_rotation"])
                    or exact_op.raw_translation != tuple(legacy_op["raw_translation"])
                    or exact_op.raw_amp != tuple(legacy_op["raw_amp"])
                    or exact_op.raw_phase != tuple(legacy_op["raw_phase"])):
                raise ValueError(f"legacy/exact SG {sg} rotation linkage mismatch")
            if any(float(raw) != value for raw, value in
                   zip(exact_op.raw_translation, legacy_op["trans"])):
                raise ValueError(f"legacy/exact SG {sg} translation linkage mismatch")
        if len(exact_source.rows) != len(irreps):
            raise ValueError(f"legacy/exact SG {sg} row linkage mismatch")
        for exact_row, legacy_row in zip(exact_source.rows, irreps):
            legacy_k = tuple(Fraction(value, legacy_row["kd"]) for value in (
                legacy_row["kx"], legacy_row["ky"], legacy_row["kz"]
            ))
            if (exact_row.source_row_ordinal != legacy_row["source_row_ordinal"]
                    or exact_row.label != legacy_row["ml_label"]
                    or exact_row.dimension != legacy_row["dim"]
                    or exact_row.k != legacy_k
                    or exact_row.raw_k != tuple(legacy_row["raw_k_tokens"])
                    or exact_row.operation_indices != tuple(legacy_row["op_indices"])
                    or exact_row.raw_characters != tuple(legacy_row["raw_character_tokens"])):
                raise ValueError(f"legacy/exact SG {sg} row token linkage mismatch")

    _validate_spin_source_sgs(all_spin_ops)

    # Sort by SG then by k-label for contiguity
    all_irreps.sort(key=lambda x: (x["sg"], x["k_label"], x["ml_label"]))

    return all_irreps, all_spin_ops


if __name__ == "__main__":
    irreps, spin_ops = parse_all_spinor()
    print(f"Parsed {len(irreps)} spinor irreps from {len(spin_ops)} SGs")
    # Show sample ops
    for sg in sorted(spin_ops.keys())[:3]:
        ops = spin_ops[sg]
        print(f"  SG{sg}: {len(ops)} spin ops")
        if ops:
            print(f"    op[0]: rot={ops[0]['rot'][:3]}... trans={ops[0]['trans']} su2={ops[0]['su2']}")
    by_sg = defaultdict(int)
    for ir in irreps:
        by_sg[ir["sg"]] += 1
    print(f"  SG range: {min(by_sg)}-{max(by_sg)}")
    print(f"  Max per SG: {max(by_sg.values())} (SG {max(by_sg, key=by_sg.get)})")
    # Show sample
    for ir in irreps[:5]:
        print(f"  SG{ir['sg']} {ir['ml_label']} dim={ir['dim']} "
              f"k=({ir['kx']}/{ir['kd']},{ir['ky']}/{ir['kd']},{ir['kz']}/{ir['kd']}) "
              f"chars={ir['characters'][:5]}...")
