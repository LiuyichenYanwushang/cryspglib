# Rust finite-domain type modernization plan

## Goal

Replace primitive values only when the domain has stable semantics and the new
type can rule out invalid or accidentally mixed states. The migration is not a
mechanical conversion of every finite database column into an enum.

The design follows three representations:

1. Small closed sets that affect control flow use enums.
2. Large bounded numeric identifiers use validated newtypes.
3. Database symbols and labels remain table data, preferably borrowed static
   strings, rather than becoming enums with hundreds or thousands of variants.

## Classification

### Enums

- `OperationKind::{Unitary, Antiunitary}` replaces ambiguous time-reversal
  booleans at typed magnetic-operation boundaries.
- `TimeReversalPolicy::{Exclude, Include}` describes whether a search should
  generate the antiunitary coset.
- `TensorParity::{Polar, Axial}` replaces the positional `is_axial` boolean.
- Existing finite semantic enums such as `MagneticType`, `CorepType`,
  `AperiodicAxis`, `RotationClass`, and `MagneticCharacterTableColumns` remain
  enums.

Binary predicates that are genuinely yes/no metadata, such as whether a
database record satisfies the Lifshitz condition, may remain `bool`.

### Validated newtypes

- `SpaceGroupNumber`: 1 through 230.
- `HallNumber`: 1 through 530.
- `UniNumber`: 1 through 1651.

These domains are finite but too large for useful enums. Newtypes prevent
mixing identifiers and validate their ranges once at the input boundary.

### Database text

Hall symbols, international symbols, BNS numbers, OG numbers, k-point labels,
and irrep labels are data, not control-flow variants. A Hall-symbol enum would
have 530 setting variants and would duplicate the existing indexed database.
It would also make database regeneration and parsing harder without improving
lookup complexity.

Metadata returned directly from generated immutable tables should borrow
`&'static str` where this does not force owned allocation. Runtime-derived or
user-owned text remains `String`. Future parsing types such as `BnsNumber` may
store their numeric components, but they should not enumerate all 1651 labels.

## Stages

1. Add the semantic enums and bounded identifier newtypes with conversion,
   display, layout, and rejection tests.
2. Use `OperationKind` in the validated magnetic-operation layer while keeping
   explicit compatibility conversion at raw parallel-array/database edges.
3. Use `TimeReversalPolicy` and `TensorParity` inside magnetic tensor search so
   call sites cannot swap two booleans with different meanings.
4. Add typed lookup and identification entry points. Existing integer entry
   points validate and delegate during the compatibility window.
5. Borrow static symbol metadata from generated tables where the returned
   object is purely database-backed. Do not convert symbols into giant enums.
6. Run release tests, doctests, strict clippy, and downstream Rustb checks.

## Expected results

- Invalid Hall/UNI/space-group numbers cannot inhabit typed identifiers.
- Magnetic unitary/antiunitary intent is visible at call sites.
- Tensor parity and time-reversal search policy cannot be passed in the wrong
  order as indistinguishable booleans.
- Database metadata lookup avoids unnecessary string allocation where possible.
- Numeric kernels and generated lookup tables retain their compact indexed
  representation.

## Verification

- Unit tests cover every identifier boundary and enum conversion/composition.
- `size_of` assertions guard the intended one-byte/two-byte layouts.
- Full release library and integration tests cover all ordinary and magnetic
  groups, including the 4479 setting-aware magnetic round trips.
- Doctests and strict all-target clippy remain clean.
- Rustb is checked with its optional `cryspglib` integration enabled.
