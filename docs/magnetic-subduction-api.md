# High-level Γ magnetic subduction API (R30)

`cryspglib::irrep::magnetic_subduction_api` is a thin public request/answer
spelling of the already verified Γ magnetic pipeline.  It adds no mathematics:
the scalar path is R25 (`magnetic_input`) + R26 (`magnetic_subduction`), the
two-dimensional real path is R29 (`magnetic_matrix_subduction`), and both are
called unchanged.

The module is registered from `src/irrep/mod.rs`; the ordinary
`irrep::subduction_api` is not modified by this round.

## Why this layer exists

The R25/R26/R29 entry points are deliberately low level: each layer resolves
its own source/condensing labels, but a caller must still find the magnetic
isotropy row's **local index**, build `RecordIdentity`, choose the scalar or
matrix engine by the source dimension, and collect the label/provenance surface
by hand.  There was no public magnetic counterpart to `irrep::subduction_api`:
no entry point that:

* resolves a condensing record from an existing `IsotropyDirection` (or a UNI
  number) against the magnetic isotropy table;
* dispatches on the source dimension to the right verified engine;
* reports the **source**, **condensing**, **magnetic record** and **target**
  identities together, with both genuine ordinary-irrep label spellings;
* keeps the result payload usable (intertwiner, target operators, matrix copies)
  without a self-referential struct.

This module supplies exactly that and nothing else.

## Example

```rust
use cryspglib::irrep::magnetic_input::{EmbeddingEvidence, RealSign, SourceIdentity};
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticRecordSelector, MagneticSubductionRequest, subduce,
};
use cryspglib::irrep::{LabelConvention, isotropy::IsotropyDirection};

let request = MagneticSubductionRequest::new(
    SourceIdentity {
        parent_sg: 5,
        label: "GM1",
        convention: LabelConvention::Cdml,
    },
    MagneticRecordSelector::new(5, "GM2", LabelConvention::Cdml)
        .direction(IsotropyDirection::Label("P1"))
        .uni(22),
    Some(RealSign::Plus),
)
.with_embedding(EmbeddingEvidence::identity());
let report = subduce(&request).expect("the pinned SG 5 GM1 / GM2 witness");
assert_eq!(report.record().uni, 22);
assert_eq!(report.source().labels.cdml, "GM1");
assert_eq!(report.source().labels.bc.as_deref(), Some("Γ1"));
assert_eq!(report.condensing().labels.cdml, "GM2");
assert_eq!(report.condensing().labels.bc.as_deref(), Some("Γ2"));
let scalar = report.scalar().expect("a one-dimensional source");
assert_eq!(scalar.target.magnetic_label, "GM1");
assert_eq!(scalar.target.h_source.cdml, "GM1");
assert!(scalar.target.magnetic_bc.is_none());
```

## Request

`MagneticSubductionRequest::new(source, record, eta)`:

| part | type | notes |
|---|---|---|
| source | `SourceIdentity` | parent SG, label, **its own** `LabelConvention`; must resolve to exactly one irrep |
| record | `MagneticRecordSelector` | parent SG, condensing label, **its own** convention |
| `eta` | `Option<RealSign>` | the parent action's `η`; required exactly when the magnetic group has an antiunitary class |
| embedding | `EmbeddingEvidence` | caller's `map`, exact `shift`, `LatticePolicy`; default identity/containment |

`MagneticRecordSelector` accepts:

* `.direction(IsotropyDirection::Label(..))` — the ISOTROPY direction label
  (`"P1"`, `"C1"`, …), the only selector that works for both ordinary and
  magnetic records;
* `.direction(IsotropyDirection::Index(..))` — position in the condensing
  irrep's magnetic isotropy table;
* `.uni(..)` — the magnetic group's UNI number, accepted alone only when that
  UNI appears in exactly one row of the condensing irrep's magnetic table;
* both direction and UNI — they must select the same row.

The supported magnetic selectors are `Label`, `Index` and the UNI number.  The
magnetic tables store no official order-parameter component string, and this
module does not invent one: `IsotropyDirection::OfficialDescriptor` is passed
through and surfaces the existing `IsotropyError::NoOfficialDescriptor`.
`IsotropyDirection::Descriptor` keeps the layer's actual behaviour: its matching
field for a magnetic row is the row's `direction` field, which holds the
ISOTROPY label, so a `Descriptor` carrying that same label follows the existing
normalization and can select the row.  A genuine component string
(`"(a,0,0)"`) has no magnetic source and surfaces the existing
`IsotropyError::DirectionNotFound`; neither case is newly restricted or
reinterpreted.

The source dimension decides the engine:

| source `dim` | engine | target family |
|---|---|---|
| 1 | R25 verify + R26 `subduce_scalar_gamma` | one-dimensional Wigner type-A Γ corep, real `±1` source |
| 2 | R29 verify + `subduce_matrix_gamma` | one- or two-dimensional type-A real target |
| other | named `SourceDimensionUnsupported` | outside both contracts |

## Identities and labels

Four identities are reported separately:

| report field | meaning | label source |
|---|---|---|
| `source()` | the subduced parent ordinary irrep | `IrrepRecord::labels()`: CDML and genuine BC |
| `condensing()` | the parent ordinary irrep whose magnetic isotropy row was selected | the row's own irrep: CDML and genuine BC |
| `record()` | the magnetic group/setting: UNI, BNS, ISO label, direction label, basis, origin, local index, ordinal | generated magnetic isotropy row |
| `scalar()/matrix().*.target` | the magnetic Γ corepresentation(s) of `M` selected by the engine | magnetic corep summary at Γ |

The target identity deliberately splits into two independent names:

* `magnetic_label` — the **display label the magnetic summary constructs** for
  the corep (`"GM1"`, `"GM2"`, …).  The summary has no separate BC field for a
  corep label, so `magnetic_bc` is always `None`.  It is never filled from the
  `H` irrep's BC spelling and must not be read as a Bradley–Cracknell name.
* `h_source` — the ordinary irrep of the unitary subgroup `H` the corep is
  induced from, with its own verified CDML/BC labels (`IrrepLabels`).

A single-source one-dimensional corep can carry the same display string as its
`H` irrep's CDML label (e.g. UNI 22: both are `"GM1"`); that is a property of the
generated data, not a reason to equate the two identities.  The SG 11 `GM1-`
source / UNI 59 witness is the complementary case: the source is `GM1-`, the
target display label is `GM2`, and the target's `H = SG 6` irrep is `GM2` with
BC `Γ2`.  The current magnetic summary has no BC source for that `GM2` corep,
so no magnetic BC name is reported or fabricated.

`magnetic_bc` being `None` is also **not** a claim that BC labels are missing
elsewhere: `source()`, `condensing()` and `h_source` do report genuine BC
spellings when the database carries them.

## Ownership and reconstruction

`MagneticSubductionReport` owns:

* `MagneticVerifiedModel` — exactly the `ScalarGammaModel` or `MatrixGammaModel`
  the engine verified (embedding, both identities, `η`, source actions, classes,
  translation lattice);
* an owned snapshot of the result payload: labels, intertwiner(s), target
  operator tables and matrix copies.

There is no struct that stores a model and simultaneously borrows a result from
it.  For a per-class complex-vector reconstruction a caller can use:

* `report.classes()` and `report.translation_lattice()` for the class list and
  the reduction lattice;
* `report.scalar_model()` / `report.matrix_model()` for the verified source
  action (`actions()` or `source_operators()`);
* `report.scalar()` → `intertwiner` + `operators` (R26 `TargetOperator` with
  `apply`/`compose`), or `report.matrix()` → pieces with `copies` (row-major
  `2×d` intertwiners) and R29 `TargetMatrixOperator` tables.

For the scalar path the statement is `D(g) = U⁻¹·A(g)·U` on unitary classes and
`D(a) = U⁻¹·A(a)·conj(U)` on antiunitary classes.  For the matrix path each copy
`U` satisfies `S(g)∘φ = φ∘T(g)`, `φ(z) = U z`, with the antiunitary branch
conjugating the image vector; `U†U = I` and the concatenated copies have full
ℂ-rank, as checked by R29.

## Errors

`MagneticSubductionApiError` keeps the layers apart:

* request resolution: `UnsupportedSpaceGroup`, `ParentMismatch`,
  `SourceLabelNotFound`, `AmbiguousSourceLabel`, `SourceDimensionUnsupported`,
  `RecordSelectorMissing`, `RecordUniNotFound`, `AmbiguousRecordUni`,
  `RecordUniMismatch`, `Isotropy(..)`;
* target provenance: `HSourceIdentityUnresolved` when the engine's ordinary `H`
  irrep identity does not resolve to exactly one database row, so its verified
  CDML/BC labels cannot be reported (this is distinct from a genuine missing-BC
  datum, which is preserved as `bc: None`);
* engine refusals: `ScalarInput(..)` (R25), `Scalar(..)` (R26),
  `Matrix(..)` (R29 verify or subduce).

Every wrapped engine error preserves its original value and
`std::error::Error::source`.  There is no fallback to an identity-only answer, no
first-hit selection among ambiguous labels, and no missing-η default.

## Verified witnesses (current working tree)

Measured through this API's building blocks (`verify_*` + `subduce_*`) with the
default identity embedding unless noted; both η signs were checked where the
group has an antiunitary class.

| source | condensing record | UNI | `H` | target | notes |
|---|---|---|---|---|---|
| SG 5 `GM1` (BC `Γ1`) | SG 5 `GM2` (`Γ2`), direction `P1` | 22 | SG 1 | display `GM1`, `H` source `GM1` (`Γ1`) | source ≠ condensing, H ≠ parent; `U = 1` / `i`; anti group |
| SG 11 `GM1+` (`Γ1+`) | SG 11 `GM2-` (`Γ2-`), `P1` | 59 | SG 6 | display `GM1`, `H` source `GM1` (`Γ1`) | source ≠ condensing, H ≠ parent; both η |
| SG 11 `GM1-` (`Γ1-`) | SG 11 `GM2-` (`Γ2-`), `P1` | 59 | SG 6 | display `GM2`, `H` source `GM2` (`Γ2`) | target display label ≠ source label; both η |
| SG 89 `GM5` (`Γ5`) | SG 89 `M1` (`M1`) | 753 | SG 89 | display `GM5`, dim 2 | anti-translation, 8 anti classes, both η; map `[[-1,-1,0],[-1,1,0],[0,0,-1]]` |
| SG 89 `GM5` | SG 89 `GM5`, row 0 | 101 | SG 3 | display `GM1` + `GM2`, dim 1 + 1 | two inequivalent 1D pieces; H ≠ parent |
| SG 89 `GM5` | SG 89 `GM5`, row 2 | 10 | SG 1 | display `GM1`, dim 1, mult 2 | two independent copies |
| SG 156 `GM3` | SG 156 `GM3`, row 2 | 1 | SG 1 | display `GM1`, dim 1, mult 2 | crosses to `P1` |
| SG 191 `GM5+` (`Γ6+`) | SG 191 `GM1-` (`Γ1-`), row 0 | 1471 | SG 177 | display `GM5`, `H` source `GM5`, dim 2 | source ≠ condensing, H ≠ parent, 12 anti classes, both η |

These are selected witnesses, not a coverage rate.  The R26/R27/R29 suites
retain their own wider witness sets.

## Scope fences

Supported exactly:

* **Γ only** (`k = 0`), non-spinor, non-compound parent irreps;
* source real `±1` one-dimensional (R25/R26) or real two-dimensional (R29);
* targets: Wigner **type A**, single-source, one- or two-dimensional Γ
  magnetic corepresentations of `M`;
* caller-chosen `(map, shift, LatticePolicy)` and `η`, verified by the layers
  below.

Not supported and refused by name or out of shape (no fallback):

* non-Γ, full-star, spinor/double-group, type-B/C, higher-dimensional sources;
* complex scalar characters, compound sources, setting search, automatic
  embedding search;
* official order-parameter descriptors for magnetic records (no table source);
* ordinary (non-magnetic) subduction through this module — use
  `irrep::subduction_api`;
* coverage claims: the 16,721 magnetic rows are candidate records, not a
  denominator of successful Γ subductions.

## Verification commands

```text
CARGO_TARGET_DIR=target CARGO_BUILD_JOBS=8 RAYON_NUM_THREADS=8 \
  cargo test --offline --test magnetic_subduction_api
CARGO_TARGET_DIR=target CARGO_BUILD_JOBS=8 RAYON_NUM_THREADS=8 \
  cargo clippy --offline --lib --test magnetic_subduction_api -- -D warnings
```

The ordinary-API regression `tests/subduction_api.rs` and the R26/R27/R29
magnetic suites are run unchanged to show that registering this module does not
alter the layers it wraps.  No dependency, commit or push is part of this round.

## R32 increment: restricted paired type-C entry point

R32 adds one incremental request/answer pair without changing the R30 surface:

```text
pub fn subduce_type_c(
    request: &MagneticSubductionRequest<'_>,
) -> Result<MagneticTypeCReport, MagneticTypeCError>
```

It resolves the same `MagneticSubductionRequest` through the same private
resolver (parent check, source label/dimension, condensing row, four identities)
and the same `verify_matrix_gamma_input`, then calls the restricted
`magnetic_matrix_subduction::subduce_matrix_gamma_type_c`.  Only a two-
dimensional source is accepted (`MagneticTypeCError::SourceDimensionNotTwo`
otherwise).

### Ownership and identity shape

`MagneticTypeCReport` owns the actual verified `MatrixGammaModel` directly; it
has no scalar branch.  `model()` and `matrix_model()` both return that model, and
`source()`, `condensing()`, `record()`, `embedding()`, `classes()` and
`translation_lattice()` are the same identity surfaces as
`MagneticSubductionReport`.

The legacy `MagneticTargetIdentity` (single `h_source` field),
`MagneticOutcome` (two variants) and `subduce` are unchanged: no field or variant
was added, so existing struct literals and exhaustive matches stay valid.  The
paired identity is a new owned type:

```text
MagneticTypeCTargetIdentity {
    magnetic_label: String,
    magnetic_bc: Option<String>,   // always None: the summary has no corep BC field
    unitary_sg: u8,
    h_sources: [IrrepLabels; 2],   // both genuine CDML/BC ordinary-irrep identities
    dimension: u32,                // 2
    multiplicity: u32,             // 1
}
```

`MagneticTypeCOutcome` carries that identity plus the normalized `2×2` complex
copies, the reconstructed `TargetMatrixOperator` table and the checked
`covered_dimension` (2).  Neither `h_sources` entry is a fabricated magnetic BC
name; `magnetic_bc` stays `None`.

### Witnesses

Both requests use `source = SG 89 GM5`, `map = I`, `shift = 0`,
`LatticePolicy::RequireExact`, `direction = Index(0)`, and both `η` signs:

| condensing | UNI | H | `h_sources` (CDML / BC) | result |
|---|---|---|---|---|
| `GM2` | 749 | SG 16 P222 | `GM3`/`Γ4`, `GM4`/`Γ2` | display `GM3 + GM4`, dim 2, mult 1, covered 2 |
| `GM4` | 751 | SG 21 C222 | `GM3`/`Γ4`, `GM4`/`Γ2` | same pair; `D(a0²) = -I` branch |

The engine independently checks the target law, the non-zero complex
intertwiner, `U†U = I₂`, full complex rank, covered dimension 2 and the
complex-vector reconstruction for every class.  The two H identities are read
from their own database rows; the source (`GM5`) and condensing (`GM2`/`GM4`)
identities remain independent, and `H` is a different space group than the
source's parent.

The legacy `subduce` on the same UNI 749 request still refuses explicitly with
`Matrix(NoMatchingTarget { examined: ["GM1", "GM2"] })`; it does not silently
drop the second H identity or claim a one-source type-A success.  A Type-A
record (for example UNI 747) is refused by the paired entry with
`NoMatchingTarget` whose `examined` list is empty.

### Errors and scope

`MagneticTypeCError` is a separate enum (`Request`, `Engine`,
`SourceDimensionNotTwo`), wrapping `MagneticSubductionApiError` or
`MatrixTypeCError`; the R30 error enum keeps its original variants and sources.
The paired engine's refusals include finite-domain exclusions (spinor,
dimension, compound, non-Γ, finite non-real, off-grid, wrong source count) and
unresolved data (missing/ambiguous H row, `NonFinite`, equivalent H pair,
unresolved corep, provenance/identity mismatch).  The constants `dimension`,
`multiplicity` and `covered_dimension` are engine-checked, not copied from a
summary label.

Out of scope: type-A/B targets through this entry, complex H characters,
non-Γ/spinor sources, full-star, setting search and external corep oracles.
The type-C classification itself is the summary's Wigner result; R32 does not
recompute Wigner `W`.

### Verification

```text
CARGO_TARGET_DIR=target CARGO_BUILD_JOBS=8 cargo test --offline \
  --test magnetic_type_c_subduction
CARGO_TARGET_DIR=target CARGO_BUILD_JOBS=8 cargo test --offline --lib \
  magnetic_matrix_subduction::tests
CARGO_TARGET_DIR=target CARGO_BUILD_JOBS=8 cargo clippy --offline --lib \
  --test magnetic_type_c_subduction -- -D warnings
```

The legacy `tests/magnetic_subduction_api.rs`, `tests/magnetic_matrix_subduction.rs`
and the R25–R31 suites are run unchanged.  No dependency, commit or push is part
of R32.
