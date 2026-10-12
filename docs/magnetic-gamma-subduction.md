# The restricted scalar Γ magnetic subduction (R26 card 1)

`irrep::magnetic_subduction` is the first **subduction** step after the R25 input
contract (`docs/magnetic-gamma-input.md`).  It consumes one verified
[`ScalarGammaModel`](../src/irrep/magnetic_input.rs) — a parent irrep's real `±1`
Γ character restricted to an embedded magnetic group `M` — and answers which
**Γ magnetic co-representation of `M`** that restriction is, with the target's
operators reconstructed on the one-dimensional space and the complex intertwiner
between the two representations exhibited.

It is one entry point, `subduce_scalar_gamma(&ScalarGammaModel)`, with private
result fields.  It is **not** a general magnetic subduction: no full-star
decomposition, no spinor or higher-dimensional source, no complex scalar
character, no setting search, and no matrix framework.

## What the answer contains

| accessor | meaning |
|---|---|
| `target_label()` | the target corep's label as `magnetic_summary` names it |
| `target_source_label()` | the unitary subgroup `H`'s ordinary irrep the target is induced from |
| `target_unitary_sg()` | `H`'s space-group number |
| `target_dimension()` / `source_dimension()` / `covered_dimension()` | the dimension identity `dim(source) = multiplicity · dim(target)` |
| `multiplicity()` | the multiplicity proved by the intertwiner, `1` in this family |
| `intertwiner()` | the complex `U` with `D = U⁻¹·A·U` (unitary classes) and `D = U⁻¹·A·conj(U)` (antiunitary classes) |
| `operators()` / `operator_at(i)` | one reconstructed `1×1` operator per class of `M`, in the model's class order |
| `classes()` / `translation_lattice()` | the classes and lattice the operators refer to |
| `scope()` | the scope fence as one string |

The result **borrows** the model instead of copying it.  The embedding, both
identities and the `±1` action table remain the model's own data and cannot be
mutated behind the report; the report's fields are private and the only
construction path is the public function.

## Independence of the target identity

The summary lookup is `magnetic_summary::magnetic_irrep_summary_by_uni_partial`,
whose Γ coreps are computed from the magnetic group's unitary subgroup `H` and the
ordinary PIR rows of `H`, a different data path from the parent character row the
model used.  The target is chosen among the Γ coreps of dimension one with exactly
one non-spinor, one-dimensional unitary source, and the module reports *that*
identity — for the pinned SG 3 `GM2` → UNI 10 witness the source label is `GM2`
(SG 3) while the target identity is `GM1` (the scalar corep of `H = P1`), so the
answer is not a restatement of the input.

The summary's antiunitary columns of a type-A corep carry **no value** at all
(`CharacterCompleteness::TypeAAntiunitaryPending`).  This module never reads
them, never forms a character inner product from an antiunitary trace, never uses
one as operator evidence, and never treats a `None` there as "the target does not
exist".  The antiunitary operators are constructed instead.

## The construction, and what is checked on it

For a candidate with unitary characters `K` (read from the summary's unitary
columns, snapped exactly to `±1`) and canonical antiunitary representative `a₀`
(the first antiunitary class in the model's class order):

```text
D(g) = K(g)            for a unitary class g
D(a) = K(a₀⁻¹·a)       for an antiunitary class a      (so D(a₀) = 1·K)
```

and `D = U₀·K^{t}` with the canonical `U₀ = 1`.  `a₀⁻¹·a` is computed by the
*corep's own* class composition, not by an inverse formula: every unitary class is
multiplied by `a₀` with `magnetic_embedding::compose` and reduced with
`canonical_translation` modulo the model's translation lattice.

Before a candidate is accepted as the target, the module checks:

1. **the M-frame correspondence** — every Γ operation of the summary maps to a
   class of the model by *rotation, time-reversal flag and translation modulo the
   model's translation lattice*, every class is covered, and every unitary column
   of one class agrees with the class's character.  A mismatch is
   `TargetOperationUnmatched` / `TargetClassUncovered` / `TargetCharacterNotAClassFunction`,
   **never** a guessed setting alignment;
2. **the square** — `D(a₀)² = D(a₀²)` with `a₀²` unitary and `K(a₀²) = +1`, the
   one-dimensional case of Wigner's type A; a candidate that fails it is refused
   (`TargetSquareInconsistent`), not silently replaced;
3. **conjugation compatibility** — `D(h)·D(a₀) = D(h·a₀)` and
   `D(a₀)·conj(D(h)) = D(a₀·h)` for every unitary class `h`;
4. **the full corepresentation law** — `D(i)·[conj(D(j)) if t(i) else D(j)] = D(i·j)`
   on every ordered pair of classes, evaluated by *composing the reconstructed
   operators* (the operator composition conjugates the right factor exactly when
   the left one is antiunitary, so no branch is special-cased).

## The intertwiner and the multiplicity

`A` and `D` are equivalent iff a complex `U ≠ 0` satisfies

```text
D(g) = U⁻¹ · A(g) · U            (unitary g)
D(a) = U⁻¹ · A(a) · conj(U)      (antiunitary a)
```

Both sides are one-dimensional with `±1` values, so this is a two-unknown real-linear
system for `U = u₁ + i·u₂`: each condition is a real row, the accumulated solution
set is the span of the first row's null direction, and a later row keeps that span
exactly when it is **orthogonal** to the direction (testing the row for parallelism
with the direction instead would misread a *repeated, consistent* constraint as a
contradiction — the SG 10 witness below has two antiunitary classes repeating the
same equation).  The module solves the system, takes `U` from the solved direction
with a fixed canonical sign (exactly `1` or `i` on the module's real `±1` tables,
the gauge's own phase on a legal complex gauge), and **re-applies the contract to
every class of both branches before accepting the candidate**.  There is no
"substitute a phase" path: a candidate without a non-zero solution is skipped, and
the result is never reported for one.  The multiplicity follows from
`dim(source) = multiplicity · dim(target)` with `dim(source)` re-read from the
source irrep's own record; in this dimension-one family it is `1`.

Because `η` is the **parent action's** datum, `η = +1` and `η = −1` are a gauge of
the antiunitary coset: they give the *same* target and **different** intertwiners
(`U = 1` versus `U = i` when the source's antiunitary value is `+1` versus `−1`).
Both are legal inputs, and the tests pin both directions on the antiunitary
witnesses.

## Pinned witnesses

| kind | witness | target | evidence |
|---|---|---|---|
| Type I (ordinary) | SG 2 `GM1+`, record `GM1+` (UNI 4, `H = P-1`), identity embedding | `GM1+`, dim 1, multiplicity 1 | the ordinary `subduce` on the same subgroup names `GM1+` with dim 1 and multiplicity 1 |
| Type III | SG 2 `GM1-` → UNI 6 (`H = P1`), both `η` | `GM1`, dim 1, multiplicity 1 | same target for both `η`; intertwiner `i` for `η = +1` and `1` for `η = −1`; the two source tables differ on the antiunitary class |
| Type III | SG 3 `GM2` → UNI 10 (`H = P1`) | `GM1` of `P1` (not the source label) | target identity independent of the source; spinor candidates rejected |
| Type IV (anti-translation) | SG 1 `X1` → UNI 3, map `2I`, both `η` | `GM1`, dim 1, multiplicity 1 | the reconstructed second class is genuinely antiunitary (half translation, `D(a₀) = 1·K`); intertwiner `1` for `η = +1` and `i` for `η = −1` |
| frame not established | SG 5 `GM1` with the `GM2` record (UNI 22), SG 11 `GM1+` with the `GM2-` record (UNI 59) | — | the summary's Γ tables are refused by name (`TargetOperationUnmatched`) |
| repeated antiunitary constraint | SG 10 `GM1+` with the `GM2+` record (UNI 53; `M` has four classes, `H` has two unitary classes) | `GM1+`, dim 1, multiplicity 1 | `η = +1` gives `U = 1`, `η = −1` gives `U = i`; the two antiunitary classes repeat the same equation |
| four antiunitary classes, non-trivial `H` character | SG 12 `GM2+` with the `Y1+` record (UNI 56) | `GM2+`, dim 1, multiplicity 1 | `η = −1` gives `U = 1`, `η = +1` gives `U = i`; the unitary character takes both `+1` and `−1`, so the coset is not one repeated equation |

The tests re-derive the law, the square, conjugation compatibility and the
intertwiner contract from the public API — on scalars *and* on actual complex
vectors (`A(g)(U z or conj(U z)) = U D(g) z`) — and cross-check the reconstructed
unitary columns against the summary's own characters, so the module's internal
checks are not taken on faith.  The solver's refusal paths and the gauge
statements are additionally covered by unit tests inside the module
(`cargo test --lib magnetic_subduction`): repeated-consistent and contradictory
systems, a wrong unitarity branch, a legal complex gauge of the coset, an
alternative antiunitary representative `a₀`, and single-class corruptions that are
*not* a uniform coset phase.  What those tests do **not** claim: a mutation that is
semantically inert on real `±1` data (for example a sign in a coefficient that is
always zero there, or removing the final re-verification) is not discriminated by
the public witnesses; the mutation log records which changes were caught and which
were not.

## Boundaries (measured, not a coverage claim)

* **The M-frame correspondence is not established in general.**  A bounded witness
  search over the magnetic rows of the one-dimensional Γ sources of SG 1–12 found
  29 referenced UNIs; 24 pass the exact correspondence above and 5 (UNI 22, 41, 46,
  59, 69) do not, because the unitary subgroup's data-Hall frame used by the
  summary and the magnetic database's own settings differ — either in the rotations
  (UNI 22) or only in the translations (UNI 59).  Those requests are refused by
  name; **no alignment transform was searched for and none was refuted**.  This is
  witness selection, not a coverage or success rate of the interface, and it is the
  same kind of table-alignment boundary R25 measured for the parent character row.
* **Grey (Type II) records** are not part of R25's stated domain and no such record
  is reachable through the magnetic isotropy rows of the witnesses here; the module
  keys on the *target's* shape (dimension one, Wigner type A, non-spinor source),
  not on the magnetic-group type, so it makes no claim about them.
* **The unresolved path** (`UnresolvedTargetData`) is implemented so a Γ source
  whose corep could not be built is never silently absent, but no pinned witness
  exercised it: the bounded search found no unresolved one-dimensional Γ source.
* **No omission is claimed beyond the examined set.**  `NoMatchingTarget` carries
  the labels examined; a candidate whose unitary character is outside the real `±1`
  grid cannot be equivalent to this source and is skipped (examined, not omitted).
* The ordinary-SG, spinor, higher-dimensional and full-star paths are untouched:
  this module adds one file, one `mod.rs` registration, one test file and this
  document.
* **Correction history.**  A first draft of the solver tested a repeated constraint
  for *parallelism* with the accumulated null direction instead of orthogonality,
  which made SG 10 `GM1+`/`GM2+` (UNI 53) look unsolvable, and the call site then
  fell back to `U = 1` instead of skipping the candidate.  Both are fixed in the
  current code (orthogonality test, no phase substitution, per-class re-verification
  before acceptance) and SG 10 is a permanent regression.  The log
  `target/logs/r26-gamma/mutations-rework.txt` records that reverting the test to
  the determinant form fails that regression; the older
  `target/logs/r26-gamma/mutations.txt` was produced while the fallback was still
  present and is kept only as history.
