#!/usr/bin/env python3
"""Exact, generation-only validation of the pinned ``spin.dat`` sources.

This module is deliberately separate from :mod:`parse_spinor_data`.  The
legacy parser keeps its historical floating-point materialisation; this file
parses the same source tokens into a small exact algebra and never imports or
uses the legacy rounding helpers.
"""

from dataclasses import dataclass
from fractions import Fraction
from collections import Counter
import hashlib
import re


class ExactSpinSourceError(ValueError):
    """A source token or exact structural invariant is invalid."""


@dataclass(frozen=True)
class Radical24:
    """A+b*sqrt(2)+c*sqrt(3)+d*sqrt(6), with exact rational coefficients."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)
    c: Fraction = Fraction(0)
    d: Fraction = Fraction(0)

    def __post_init__(self):
        if not all(isinstance(x, Fraction) for x in (self.a, self.b, self.c, self.d)):
            raise TypeError("Radical24 coefficients must be Fraction")

    def __add__(self, other):
        return Radical24(
            self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d
        )

    def __sub__(self, other):
        return Radical24(
            self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d
        )

    def __neg__(self):
        return Radical24(-self.a, -self.b, -self.c, -self.d)

    def __mul__(self, other):
        a, b, c, d = self.a, self.b, self.c, self.d
        e, f, g, h = other.a, other.b, other.c, other.d
        return Radical24(
            a * e + 2 * b * f + 3 * c * g + 6 * d * h,
            a * f + b * e + 3 * c * h + 3 * d * g,
            a * g + c * e + 2 * b * h + 2 * d * f,
            a * h + d * e + b * g + c * f,
        )

    def __rmul__(self, other):
        if isinstance(other, int):
            other = Fraction(other)
        if isinstance(other, Fraction):
            return Radical24(
                other * self.a, other * self.b, other * self.c, other * self.d
            )
        return NotImplemented

    def is_zero(self):
        return self == ZERO


ZERO = Radical24()
ONE = Radical24(Fraction(1))
SQRT2 = Radical24(b=Fraction(1))
SQRT3 = Radical24(c=Fraction(1))
SQRT6 = Radical24(d=Fraction(1))


@dataclass(frozen=True)
class Complex24:
    re: Radical24 = ZERO
    im: Radical24 = ZERO

    def __add__(self, other):
        return Complex24(self.re + other.re, self.im + other.im)

    def __sub__(self, other):
        return Complex24(self.re - other.re, self.im - other.im)

    def __neg__(self):
        return Complex24(-self.re, -self.im)

    def __mul__(self, other):
        return Complex24(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    def conjugate(self):
        return Complex24(self.re, -self.im)

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Complex24 powers require a nonnegative integer")
        result = Complex24(ONE, ZERO)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result


ZETA24 = Complex24(
    Radical24(b=Fraction(1, 4), d=Fraction(1, 4)),
    Radical24(b=Fraction(-1, 4), d=Fraction(1, 4)),
)
ZETA24_POWERS = tuple(ZETA24**n for n in range(24))
if len(set(ZETA24_POWERS)) != 24 or ZETA24**24 != Complex24(ONE, ZERO):
    raise AssertionError("zeta24 does not have exact order 24")
if any(ZETA24**n == Complex24(ONE, ZERO) for n in range(1, 24)):
    raise AssertionError("zeta24 has a smaller positive period")


@dataclass(frozen=True)
class ExactSpinOperation:
    rotation: tuple
    translation: tuple
    su2: tuple
    file: str = ""
    line: int = 0
    raw_rotation: tuple = ()
    raw_translation: tuple = ()
    raw_amp: tuple = ()
    raw_phase: tuple = ()


@dataclass(frozen=True)
class ExactSpinRow:
    source_row_ordinal: int
    dimension: int
    k: tuple
    operation_indices: tuple
    characters: tuple
    file: str = ""
    line: int = 0
    raw_k: tuple = ()
    raw_operation_indices: tuple = ()
    raw_dimension: str = ""
    raw_characters: tuple = ()


@dataclass(frozen=True)
class ExactSpinFile:
    sg: int
    operations: tuple
    rows: tuple
    path: str


def _location(path, line, column, message):
    return ExactSpinSourceError(f"{path}:{line}:{column}: {message}")


def _tokens(text):
    return tuple((m.group(), m.start() + 1) for m in re.finditer(r"\S+", text))


def _lookup(token, table, path, line, column, label):
    try:
        return table[token]
    except KeyError as exc:
        raise _location(path, line, column, f"invalid {label} spelling {token!r}") from exc


def _fraction_table(values):
    return {key: Fraction(value) for key, value in values.items()}


K_TOKENS = _fraction_table({
    "-1.0": -1, "-0.5": Fraction(-1, 2), "0.0": 0,
    "0.333333": Fraction(1, 3), "0.5": Fraction(1, 2),
    "1.0": 1, "1.5": Fraction(3, 2),
})
TRANS_TOKENS = _fraction_table({
    "0.0": 0, "0.16667": Fraction(1, 6), "0.25": Fraction(1, 4),
    "0.33333": Fraction(1, 3), "0.5": Fraction(1, 2),
    "0.66667": Fraction(2, 3), "0.75": Fraction(3, 4),
    "0.83333": Fraction(5, 6),
})
AMP_TOKENS = {
    "0.0": ZERO, "0.70711": Radical24(b=Fraction(1, 2)), "1.0": ONE
}
DIM_TOKENS = {str(n): n for n in range(1, 5)}
ROT_TOKENS = {str(n): n for n in (-1, 0, 1)}
PHASE_EXPONENTS = {
    "-1": 12, "-5/6": 14, "-3/4": 15, "-2/3": 16,
    "-1/2": 18, "-1/3": 20, "-1/4": 21, "-1/6": 22,
    "0": 0, "1/6": 2, "1/4": 3, "1/3": 4, "1/2": 6,
    "2/3": 8, "3/4": 9, "5/6": 10, "1": 12,
}
PHASE_TOKENS = {
    "-1.0": 12, "-0.83333": 14, "-0.75": 15, "-0.66667": 16,
    "-0.5": 18, "-0.33333": 20, "-0.25": 21, "-0.16667": 22,
    "0.0": 0, "0.16667": 2, "0.25": 3, "0.33333": 4,
    "0.5": 6, "0.66667": 8, "0.75": 9, "0.83333": 10, "1.0": 12,
}

DIRECT_CHAR_TOKENS = {
    **{str(n) + ".0": Radical24(Fraction(n)) for n in range(-4, 5)},
    "-1.41421": -SQRT2, "1.41421": SQRT2,
    "-1.73205": -SQRT3, "1.73205": SQRT3,
}
DIRECT_CHAR_TOKENS.update({"1e-05": ZERO})

POLAR_AMP_TOKENS = {
    "0.0": ZERO, "1.0": ONE,
    "0.70711": Radical24(b=Fraction(1, 2)),
    "1.41421": SQRT2, "1.73205": SQRT3,
    "2.0": 2 * ONE, "3.0": 3 * ONE, "4.0": 4 * ONE,
}
SU2_POLAR_PAIRS = frozenset({
    ("0.0", "0.0"),
    *{("0.70711", phase) for phase in
      ("-0.25", "-0.5", "-0.75", "0.0", "0.25", "0.5", "0.75", "1.0")},
    *{("1.0", phase) for phase in
      ("-0.16667", "-0.25", "-0.33333", "-0.5", "-0.66667", "-0.75",
       "0.0", "0.16667", "0.25", "0.33333", "0.5", "0.66667", "0.75",
       "0.83333", "1.0")},
})
CHARACTER_POLAR_PAIRS = frozenset({
    ("0.0", phase) for phase in
    ("-0.16667", "-0.25", "0.0", "0.16667", "0.25", "0.5")
} | {
    ("1.0", phase) for phase in
    ("-0.16667", "-0.25", "-0.33333", "-0.5", "-0.66667", "-0.75",
     "-0.83333", "-1.0", "0.0", "0.16667", "0.25", "0.33333", "0.5",
     "0.66667", "0.75", "0.83333", "1.0")
} | {
    ("1.41421", phase) for phase in ("-0.5", "0.0", "0.5", "1.0")
} | {
    ("1.73205", phase) for phase in ("-0.5", "0.0", "0.5", "1.0")
} | {
    ("2.0", phase) for phase in ("-0.33333", "-0.5", "0.0", "0.33333", "0.5", "1.0")
} | {
    ("3.0", "0.0"), ("4.0", "0.0")
})


def _complex_from_polar(amp_token, phase_token, path, line, amp_col, phase_col):
    amp = _lookup(amp_token, POLAR_AMP_TOKENS, path, line, amp_col, "polar amplitude")
    exponent = _lookup(phase_token, PHASE_TOKENS, path, line, phase_col, "polar phase")
    if amp.is_zero():
        return Complex24()
    return Complex24(amp * ZETA24_POWERS[exponent].re, amp * ZETA24_POWERS[exponent].im)


def _parse_operation(tokens, path, line):
    if len(tokens) != 20:
        raise _location(path, line, 1, f"operation must have exactly 20 tokens, got {len(tokens)}")
    rot = tuple(_lookup(tok, ROT_TOKENS, path, line, col, "rotation") for tok, col in tokens[:9])
    trans = tuple(_lookup(tok, TRANS_TOKENS, path, line, col, "translation") for tok, col in tokens[9:12])
    entries = []
    for (amp, amp_col), (phase, phase_col) in zip(tokens[12:16], tokens[16:20]):
        if (amp, phase) not in SU2_POLAR_PAIRS:
            raise _location(path, line, amp_col, "invalid SU2 amplitude/phase pair")
        entries.append(_complex_from_polar(amp, phase, path, line, amp_col, phase_col))
    return ExactSpinOperation(
        rot, trans, tuple(entries), path, line,
        tuple(tok for tok, _ in tokens[:9]),
        tuple(tok for tok, _ in tokens[9:12]),
        tuple(tok for tok, _ in tokens[12:16]),
        tuple(tok for tok, _ in tokens[16:20]),
    )


def _parse_direct_character(token, path, line, column, sg, row, character_column):
    if token == "1e-05":
        allowed = {
            (193, 2, 13), (193, 2, 14), (194, 2, 13), (194, 2, 14)
        }
        if (sg, row, character_column) not in allowed:
            raise _location(path, line, column, "1e-05 is only allowed at pinned zero positions")
    return Complex24(_lookup(token, DIRECT_CHAR_TOKENS, path, line, column, "direct character"), ZERO)


def parse_spinor_file_exact(path):
    with open(path, encoding="utf-8") as stream:
        lines = stream.readlines()
    sg = None
    nsym = None
    sym_index = None
    for index, raw in enumerate(lines):
        text = raw.strip()
        if text.startswith("SG="):
            sg = int(text[3:])
        elif text.startswith("nsym="):
            nsym = int(text[5:].strip())
        elif text == "symmetries=":
            sym_index = index + 1
            break
    if sg is None or nsym is None or sym_index is None:
        raise ExactSpinSourceError(f"{path}: missing SG/nsym/symmetries header")
    operations = []
    index = sym_index
    while len(operations) < nsym:
        if index >= len(lines) or not lines[index].strip():
            raise ExactSpinSourceError(f"{path}:{index + 1}: missing symmetry operation")
        operations.append(_parse_operation(_tokens(lines[index]), path, index + 1))
        index += 1

    rows = []
    ordinal = 0
    current_k = None
    current_ops = None
    while index < len(lines):
        raw = lines[index]
        text = raw.strip()
        if not text:
            index += 1
            continue
        if text.startswith("kpoint"):
            pieces = text.split(":")
            if len(pieces) != 3:
                raise _location(path, index + 1, 1, "malformed kpoint line")
            k_tokens = _tokens(pieces[1])
            if len(k_tokens) != 3:
                raise _location(path, index + 1, 1, "kpoint must have three coordinates")
            current_k = tuple(_lookup(tok, K_TOKENS, path, index + 1, col, "k coordinate")
                              for tok, col in k_tokens)
            raw_k = tuple(tok for tok, _ in k_tokens)
            op_tokens = _tokens(pieces[2])
            if not op_tokens:
                raise _location(path, index + 1, 1, "kpoint has no operation indices")
            parsed_ops = []
            for tok, col in op_tokens:
                if not tok.isdigit() or str(int(tok)) != tok:
                    raise _location(path, index + 1, col, "invalid operation-index spelling")
                op = int(tok)
                if op < 1 or op > nsym:
                    raise _location(path, index + 1, col, "operation index out of range")
                parsed_ops.append(op - 1)
            if len(set(parsed_ops)) != len(parsed_ops):
                raise _location(path, index + 1, 1, "duplicate operation index")
            current_ops = tuple(parsed_ops)
            raw_operation_indices = tuple(tok for tok, _ in op_tokens)
        elif text.startswith("-"):
            if current_k is None or current_ops is None:
                raise _location(path, index + 1, 1, "row precedes a kpoint")
            tokens = _tokens(text)
            if len(tokens) < 3:
                raise _location(path, index + 1, 1, "truncated irrep row")
            label, _ = tokens[0]
            dim_token, dim_column = tokens[1]
            dim = _lookup(dim_token, DIM_TOKENS, path, index + 1, dim_column, "dimension")
            values = tokens[2:]
            if len(values) not in (len(current_ops), 2 * len(current_ops)):
                raise _location(path, index + 1, 1, "character count does not match operation columns")
            if len(values) == len(current_ops):
                chars = tuple(_parse_direct_character(tok, path, index + 1, col, sg, ordinal, colno)
                              for colno, (tok, col) in enumerate(values))
            else:
                amplitudes = values[:len(current_ops)]
                phases = values[len(current_ops):]
                chars_list = []
                for colno, ((amp, amp_col), (phase, phase_col)) in enumerate(zip(amplitudes, phases)):
                    if (amp, phase) not in CHARACTER_POLAR_PAIRS:
                        raise _location(path, index + 1, amp_col, "invalid character amplitude/phase pair")
                    chars_list.append(_complex_from_polar(
                        amp, phase, path, index + 1, amp_col, phase_col
                    ))
                chars = tuple(chars_list)
            rows.append(ExactSpinRow(
                ordinal, dim, current_k, current_ops, chars, path, index + 1,
                raw_k, raw_operation_indices, dim_token, tuple(tok for tok, _ in values),
            ))
            ordinal += 1
        else:
            raise _location(path, index + 1, 1, "unexpected source line")
        index += 1
    return ExactSpinFile(sg, tuple(operations), tuple(rows), path)


def source_paths():
    """Return the cryptographically verified pinned source paths."""
    import parse_spinor_data
    return parse_spinor_data.verified_spin_source_paths()


def parse_all_exact(paths=None):
    """Parse and validate every pinned source file without float conversion."""
    files = tuple(source_paths() if paths is None else paths)
    if len(files) != 230:
        raise ExactSpinSourceError(f"expected 230 source files, found {len(files)}")
    parsed = tuple(parse_spinor_file_exact(path) for path in files)
    if {item.sg for item in parsed} != set(range(1, 231)):
        raise ExactSpinSourceError("source files must contain each SG exactly once")
    return parsed


DOMAIN_ORDER = (
    "rot", "trans", "su2_amp", "su2_phase", "k_coord", "op_index", "dim",
    "char_real", "char_amp", "char_phase",
)
DOMAIN_GOLDENS = {
    "rot": "fee9f6b012fe5508ecbca11458cb4c575f568c1adeb6a2e2cf787c5e199d5335",
    "trans": "26f763565e8bdef02628b0cd5f277194167e516904c771900fcfdb9f0eee90ba",
    "su2_amp": "9147ebb64a2da5aeb69102bba4b689e1c1597d91359e03cd8dab38a20c6bee74",
    "su2_phase": "5fc4baa254d4161b08782b961d45526f746c13c6ea590ea834d2106598f847bb",
    "k_coord": "c8d0bda39edf2b57ad692b722372c4a3687f2fcf56b7c6ff9ea264ae3e1105e7",
    "op_index": "bc1a7bf5dda0caf6311f602b6b8519cc571a9afa82aae3f63a0f6614cc7d7850",
    "dim": "0f35794ba51b1a4a267ba3fd18f16a7368ea14f75ba768028b754487a918b949",
    "char_real": "517609d78b98ef123f99063bcd26d40df02b140173a00b7effafaf2599a86e7e",
    "char_amp": "eb6d238e4650ebfca5ecd8cfb2b9371b5b63e740ae2d60ac2349fc3a764dd312",
    "char_phase": "8460133dece11f009c2238a990744e75bf67d7bd51e44f7d5aa68bd0f5f290f8",
}
COMBINED_GOLDEN = "a6af59f9c591167ec7c83c98d354f4a6f12d79e91afefdb6f9f4c2e56654b906"
SU2_PAIR_GOLDEN = "4b6e853dbda39d54f67a9c6ef7ed6feceff8d3bad6ed5236a60084c057303d66"
CHARACTER_PAIR_GOLDEN = "8a2ec58566c578b8bb97d71fa4cf163a5c0baa64257cb2c4f1d8af9dca17005a"


def _digest_counter(counter):
    payload = "".join(
        f"{spelling}\t{counter[spelling]}\n" for spelling in sorted(counter)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def _digest_pair_counter(counter):
    payload = "".join(
        f"{left}\t{right}\t{counter[(left, right)]}\n"
        for left, right in sorted(counter)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def source_spelling_census(files):
    domains = {name: Counter() for name in DOMAIN_ORDER}
    su2_pairs = Counter()
    character_pairs = Counter()
    counts = Counter(files=len(files), operations=0, kblocks=0, rows=0,
                     real_rows=0, polar_rows=0, character_columns=0,
                     operation_indices=0, su2_pairs=0, character_pairs=0)
    dimensions = Counter()
    for source in files:
        counts["operations"] += len(source.operations)
        for operation in source.operations:
            domains["rot"].update(operation.raw_rotation)
            domains["trans"].update(operation.raw_translation)
            domains["su2_amp"].update(operation.raw_amp)
            domains["su2_phase"].update(operation.raw_phase)
            su2_pairs.update(zip(operation.raw_amp, operation.raw_phase))
            counts["su2_pairs"] += 4
        seen_blocks = set()
        for row in source.rows:
            counts["rows"] += 1
            dimensions[row.raw_dimension] += 1
            domains["dim"].update((row.raw_dimension,))
            block_key = (row.raw_k, row.raw_operation_indices)
            if block_key not in seen_blocks:
                seen_blocks.add(block_key)
                counts["kblocks"] += 1
                domains["k_coord"].update(row.raw_k)
                domains["op_index"].update(row.raw_operation_indices)
                counts["operation_indices"] += len(row.raw_operation_indices)
            counts["character_columns"] += len(row.operation_indices)
            if len(row.raw_characters) == len(row.operation_indices):
                counts["real_rows"] += 1
                domains["char_real"].update(row.raw_characters)
            else:
                counts["polar_rows"] += 1
                n = len(row.operation_indices)
                domains["char_amp"].update(row.raw_characters[:n])
                domains["char_phase"].update(row.raw_characters[n:])
                character_pairs.update(zip(row.raw_characters[:n], row.raw_characters[n:]))
                counts["character_pairs"] += n
    return counts, dimensions, domains, su2_pairs, character_pairs


def source_spelling_hashes(files):
    _counts, _dimensions, domains, su2_pairs, character_pairs = source_spelling_census(files)
    hashes = {name: _digest_counter(domains[name]) for name in DOMAIN_ORDER}
    hashes["combined"] = hashlib.sha256("".join(
        f"{domain}\t{spelling}\t{domains[domain][spelling]}\n"
        for domain in DOMAIN_ORDER for spelling in sorted(domains[domain])
    ).encode()).hexdigest()
    hashes["su2_pair"] = _digest_pair_counter(su2_pairs)
    hashes["character_pair"] = _digest_pair_counter(character_pairs)
    return hashes


def _matrix_multiply(left, right):
    return tuple(
        sum((left[2 * row + middle] * right[2 * middle + column]
             for middle in range(2)), Complex24())
        for row in range(2) for column in range(2)
    )


def _matrix_adjoint(matrix):
    return (matrix[0].conjugate(), matrix[2].conjugate(),
            matrix[1].conjugate(), matrix[3].conjugate())


def _matrix_identity():
    return (Complex24(ONE, ZERO), Complex24(), Complex24(), Complex24(ONE, ZERO))


def _matrix_is_negative(left, right):
    return all(a == -b for a, b in zip(left, right))


def _validate_source_row_ordinals(source):
    expected = tuple(range(len(source.rows)))
    actual = tuple(row.source_row_ordinal for row in source.rows)
    if actual != expected:
        raise ExactSpinSourceError(
            f"SG {source.sg}: source row ordinals are not pinned raw order"
        )


def _mod_one(value):
    return value % 1


def _translation_lattice_cosets(operations):
    """Generate the finite factor group Lambda/Z^3 from source Seitz data."""
    by_rotation = {}
    for operation in operations:
        by_rotation.setdefault(operation.rotation, []).append(operation.translation)
    generators = set()
    for translations in by_rotation.values():
        for left in translations:
            for right in translations:
                generators.add(tuple(_mod_one(a - b) for a, b in zip(left, right)))
    rotations = {operation.rotation: operation for operation in operations}
    for left in operations:
        for right in operations:
            rotation = tuple(sum(left.rotation[3 * row + m] * right.rotation[3 * m + col]
                                for m in range(3)) for row in range(3) for col in range(3))
            target = rotations.get(rotation)
            if target is None:
                raise ExactSpinSourceError("rotation product has no source representative")
            generators.add(tuple(_mod_one(left.translation[row] + sum(
                Fraction(left.rotation[3 * row + m]) * right.translation[m]
                for m in range(3)
            ) - target.translation[row]) for row in range(3)))
    cosets = {(Fraction(0), Fraction(0), Fraction(0))}
    frontier = list(cosets)
    while frontier:
        old = frontier.pop()
        for generator in generators:
            new = tuple(_mod_one(a + b) for a, b in zip(old, generator))
            if new not in cosets:
                cosets.add(new)
                frontier.append(new)
    return frozenset(cosets)


def validate_exact_sources(files):
    """Run exact source, SU(2), lattice, and character structural checks."""
    if len(files) != 230 or {source.sg for source in files} != set(range(1, 231)):
        raise ExactSpinSourceError("source file SG census is not exactly 1..230")
    products = 0
    same_lift = 0
    negative_lift = 0
    for source in files:
        _validate_source_row_ordinals(source)
        lattice_cosets = _translation_lattice_cosets(source.operations)
        rotations = {}
        for index, operation in enumerate(source.operations):
            if operation.rotation in rotations:
                raise ExactSpinSourceError(f"SG {source.sg}: duplicate source rotation")
            rotations[operation.rotation] = index
            if _matrix_multiply(_matrix_adjoint(operation.su2), operation.su2) != _matrix_identity():
                raise ExactSpinSourceError(f"SG {source.sg}: SU(2) unitarity failure")
            determinant = operation.su2[0] * operation.su2[3] - operation.su2[1] * operation.su2[2]
            if determinant != Complex24(ONE, ZERO):
                raise ExactSpinSourceError(f"SG {source.sg}: SU(2) determinant failure")
        identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        identity_ops = [
            operation for operation in source.operations
            if operation.rotation == identity and operation.translation == (Fraction(0),) * 3
            and operation.su2 == _matrix_identity()
        ]
        if len(identity_ops) != 1:
            raise ExactSpinSourceError(f"SG {source.sg}: missing or duplicate identity lift")
        for left in source.operations:
            for right in source.operations:
                products += 1
                rotation = tuple(sum(left.rotation[3 * row + m] * right.rotation[3 * m + col]
                                    for m in range(3)) for row in range(3) for col in range(3))
                target_index = rotations.get(rotation)
                if target_index is None:
                    raise ExactSpinSourceError(f"SG {source.sg}: rotation product missing")
                target = source.operations[target_index]
                shift = tuple(left.translation[row] + sum(
                    Fraction(left.rotation[3 * row + m]) * right.translation[m] for m in range(3)
                ) - target.translation[row] for row in range(3))
                if tuple(_mod_one(value) for value in shift) not in lattice_cosets:
                    raise ExactSpinSourceError(f"SG {source.sg}: non-lattice Seitz product")
                product = _matrix_multiply(left.su2, right.su2)
                if product == target.su2:
                    same_lift += 1
                elif _matrix_is_negative(product, target.su2):
                    negative_lift += 1
                else:
                    raise ExactSpinSourceError(f"SG {source.sg}: SU(2) lift product mismatch")
        for row in source.rows:
            if len(set(row.operation_indices)) != len(row.operation_indices):
                raise ExactSpinSourceError(f"{row.file}:{row.line}: duplicate row operation")
            if len(row.characters) != len(row.operation_indices):
                raise ExactSpinSourceError(f"{row.file}:{row.line}: character width mismatch")
            identity_columns = [
                index for index in row.operation_indices
                if source.operations[index].rotation == identity
                and source.operations[index].translation == (Fraction(0),) * 3
            ]
            if len(identity_columns) != 1 or row.characters[identity_columns[0]].im != ZERO:
                raise ExactSpinSourceError(f"{row.file}:{row.line}: character identity failure")
            if row.characters[identity_columns[0]].re != Radical24(Fraction(row.dimension)):
                raise ExactSpinSourceError(f"{row.file}:{row.line}: character dimension failure")
    return products, same_lift, negative_lift
