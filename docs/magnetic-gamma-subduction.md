# The restricted scalar Γ magnetic subduction (R26 card 1, R27 provenance card)

`irrep::magnetic_subduction` is the first **subduction** step after the R25 input
contract (`docs/magnetic-gamma-input.md`).  It consumes one verified
[`ScalarGammaModel`](../src/irrep/magnetic_input.rs) — a parent irrep's real `±1`
Γ character restricted to an embedded magnetic group `M` — and answers which
**Γ magnetic co-representation of `M`** that restriction is, with the target's
operators reconstructed on the one-dimensional space and the complex intertwiner
between the two representations exhibited.

R27 changes exactly one thing: the correspondence between the target summary's Γ
table and the classes of `M` is now resolved through the summary's raw operation
index instead of by comparing the table's displayed coordinates (see
[the provenance section](#the-raw-operation-provenance-r27)).

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

## The raw operation provenance (R27)

The Γ summary's operations are stored in the **H data-Hall frame**: their
`rotation`/`translation` fields are not coordinates of `M`'s own first-Hall
setting, so comparing them with the model's classes is a frame confusion.  Each
column carries `magnetic_operation_index`, the index of the raw first-Hall
operation it came from, and that index is the provenance link.

Before any index is resolved, the module rebuilds the exact first-Hall list with
`magnetic_embedding::magnetic_operations(uni)` — the same entry point the model's
classes come from — and requires it to be **element-wise identical** to the
`SymmetryOps::from_magnetic_database(uni)` list the summary is computed from:
same length, same order, same rotation, same exact translation, same
time-reversal flag, for every entry.  A divergence is
`TargetRawOperationsDiffer` / `TargetRawOperationsUnavailable`.

The index is then resolved into that list, and the resulting raw operation is
mapped to a class of `M` by rotation, time-reversal flag and translation modulo
the model's own translation lattice (the same exact
`canonical_translation`/`compose` arithmetic the model was verified with).  The
structure around it is checked as well:

* the summary lists as many Γ columns as the raw list has operations (the Γ little
  group is all of `M`), each index is in range and each raw operation is indexed
  exactly once, and each character column appears exactly once and in range —
  `TargetOperationCountMismatch`, `TargetOperationIndexOutOfRange`,
  `TargetOperationIndexRepeated`, `TargetColumnOutOfRange`, `TargetColumnRepeated`.
  An incomplete column set can only appear as a repeat or an out-of-range column and
  is named by one of those two;
* the **position association** is exact, not merely a permutation.  The summary
  builds its Γ columns by enumerating the little-group indices in increasing order,
  and Γ's little group is the whole group, so at position `i` the character column
  must be `i` and the raw operation index must be `i`.  This is checked at the
  filter level for all 1,651 UNIs (the Γ little group is `0..len` in raw order for
  every one) and on the pinned summaries; a swap of two legal, unique columns or
  indices — **including two operations of the same time-reversal flag** — is
  `TargetColumnMisassociated` / `TargetOperationIndexMisassociated`, never an
  accepted alternative correspondence;
* each column's displayed time-reversal flag must equal the indexed raw
  operation's, the one frame-independent field the two lists share —
  `TargetTimeReversalMismatch`;
* an off-grid displayed translation is still a data-integrity refusal
  (`TargetTranslationOffGrid`), even though the value is not used for the
  correspondence;
* every class must be covered (`TargetClassUncovered`) and every unitary column of
  one class must agree on the character (`TargetCharacterNotAClassFunction`).

There is no fallback: a failure in any of these checks refuses the request by name,
and the displayed frame is never used to patch a missing correspondence.

**What is not verifiable per call.**  A permutation *inside the raw operation list
itself* cannot be seen from a summary column, because the column's index is well
defined relative to whatever list the API returns.  That is not a residual hole in
the association: the list is re-read fresh from two public surfaces the module
proves element-wise identical (`magnetic_operations(uni)` and
`SymmetryOps::from_magnetic_database(uni)`), and the model's classes come from that
same list, so a reordered database reorders the summary's index, the raw list and
the model's classes consistently.  The association the module can and does verify
is the position link between summary columns and that list.

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

1. **the operation provenance and the M-frame correspondence** — every Γ column is
   resolved through its raw operation index (after the element-wise provenance check
   above), the index list must cover the raw list exactly once at the right
   positions, each character column must be in range, unique and the one at that
   position, each column's displayed time-reversal flag must be the indexed raw
   operation's, and the *raw* operation must map to a class of the model by
   rotation, time-reversal flag and translation modulo the model's translation
   lattice; every class must be covered and every unitary column of one class must
   agree with the class's character.  The displayed
   `rotation`/`translation` of a column are the H data-Hall frame and are **never**
   compared with `M`'s frame.  A mismatch is `TargetRawOperationsDiffer` /
   `TargetOperationCountMismatch` / `TargetOperationIndexOutOfRange` /
   `TargetOperationIndexRepeated` / `TargetOperationIndexMisassociated` /
   `TargetTimeReversalMismatch` / `TargetColumnOutOfRange` /
   `TargetColumnRepeated` / `TargetColumnMisassociated` /
   `TargetOperationUnmatched` / `TargetClassUncovered` /
   `TargetCharacterNotAClassFunction`, **never** a guessed setting alignment and
   never a fallback to the displayed frame;
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
| data-Hall frame, index-resolved | SG 5 `GM1` with the `GM2` record (UNI 22, `H = P1`); SG 11 `GM1+` with the `GM2-` record (UNI 59, `H = SG 6`); both `η` | `GM1` of `H` in both cases, dim 1, multiplicity 1 | the summary's Γ table is displayed in a frame that does not match `M`'s classes, so the correspondence goes through the raw operation index; the target's unitary characters are `+1` (its antiunitary columns are `None` and are never read), the square/conjugation/full law and the non-zero `U` are re-derived per class, `U = 1` for `η = +1` and `U = i` for `η = −1`, and the two `η` give the same operator table |
| raw operation surface | all 1,651 UNIs (`the_two_raw_operation_surfaces_agree_for_every_uni`) | — | `magnetic_operations(uni)` is element-wise identical to `SymmetryOps::from_magnetic_database(uni)`: length, order, rotation, exact translation and time-reversal flag |
| repeated antiunitary constraint | SG 10 `GM1+` with the `GM2+` record (UNI 53; `M` has four classes, `H` has two unitary classes) | `GM1+`, dim 1, multiplicity 1 | `η = +1` gives `U = 1`, `η = −1` gives `U = i`; the two antiunitary classes repeat the same equation |
| four antiunitary classes, non-trivial `H` character | SG 12 `GM2+` with the `Y1+` record (UNI 56) | `GM2+`, dim 1, multiplicity 1 | `η = −1` gives `U = 1`, `η = +1` gives `U = i`; the unitary character takes both `+1` and `−1`, so the coset is not one repeated equation |

The tests re-derive the law, the square, conjugation compatibility and the
intertwiner contract from the public API — on scalars *and* on actual complex
vectors (`A(g)(U z or conj(U z)) = U D(g) z`) — and cross-check the reconstructed
unitary columns against the summary's own characters *through the raw operation
index*, so the module's internal checks are not taken on faith.  The solver's
refusal paths, the provenance guards and the gauge statements are additionally
covered by unit tests inside the module
(`cargo test --lib magnetic_subduction`): repeated-consistent and contradictory
systems, a wrong unitarity branch, a legal complex gauge of the coset, an
alternative antiunitary representative `a₀`, single-class corruptions that are
*not* a uniform coset phase, and the R27 guards — a shortened summary list, a raw
list of the wrong length, an out-of-range or repeated operation index, a
time-reversal flag that disagrees with the index, a reordered or relabelled raw
list, a rewritten raw rotation and translation, a repeated or out-of-range
character column, and a same-flag index swap and a column swap (the position
mis-associations) — each asserted against its exact named variant.  What those tests do **not** claim: a mutation that is
semantically inert on real `±1` data (for example a sign in a coefficient that is
always zero there, or removing the final re-verification) is not discriminated by
the public witnesses; the mutation log records which changes were caught and which
were not.

## Boundaries (measured, not a coverage claim)

* **The R26 frame boundary is replaced by the measured index route.**  A bounded
  witness search over the magnetic rows of the one-dimensional Γ sources of SG 1–12
  (R26) found 29 referenced UNIs; 24 passed the exact displayed-frame correspondence
  and 5 (UNI 22, 41, 46, 59, 69) did not, because the unitary subgroup's data-Hall
  frame used by the summary and the magnetic database's own settings differ — either
  in the rotations (UNI 22) or only in the translations (UNI 59).  R27 resolves the
  correspondence through the raw operation index instead of the displayed
  coordinates, on the run-time proof that the two raw surfaces are element-wise the
  same list (checked for all 1,651 UNIs).  **Measured in R27**: UNI 22 and UNI 59 now
  succeed with an independently read target identity, checked characters and a
  verified non-zero intertwiner for both `η` signs.  **Not re-measured in R27**: UNI
  41, 46 and 69 — the index route is the same code path, but no witness for them was
  constructed, so this document does not claim they succeed.  The bounded walk
  remains witness selection, not a coverage or success rate of the interface.
* **The provenance link rests on the database's order, not on a setting search.**
  Both raw surfaces are read fresh for the requested UNI and compared element by
  element; the summary's Γ list enumerates the little group in increasing raw order
  (checked at the filter level for all 1,651 UNIs) and the module requires the
  position association `column == magnetic_operation_index == position`, so a swap
  of two legal, unique columns or indices is refused by name even when both
  operations carry the same time-reversal flag (see the provenance section above).
  No alignment transform was searched for, and the displayed frame is never used as
  evidence.
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
  R26 added one file, one `mod.rs` registration, one test file and this document,
  and R27 changes only the mapping inside this module, its tests and this document.
  The R25 input contract, the corep/summary/wigner computation, the summary data
  structures, `api.rs` and `mod.rs` are unchanged.
* **Correction history.**  A first draft of the solver tested a repeated constraint
  for *parallelism* with the accumulated null direction instead of orthogonality,
  which made SG 10 `GM1+`/`GM2+` (UNI 53) look unsolvable, and the call site then
  fell back to `U = 1` instead of skipping the candidate.  Both are fixed in the
  current code (orthogonality test, no phase substitution, per-class re-verification
  before acceptance) and SG 10 is a permanent regression.  The log
  `target/logs/r26-gamma/mutations-rework.txt` records that reverting the test to
  the determinant form fails that regression; the older
  `target/logs/r26-gamma/mutations.txt` was produced while the fallback was still
  present and is kept only as history.  R27 then replaced the R26 displayed-frame
  matching: the R26 entry point refused UNI 22/59 by name
  (`TargetOperationUnmatched`), which was correct given the frame mismatch but left
  the correspondence unestablished — the index route establishes it, and the two
  UNIs are permanent regressions of the new behavior.  A first R27 revision accepted
  any legal, unique permutation of columns and indices and called a same-flag index
  swap an unavoidable residual; the position association is in fact checkable (the
  summary enumerates Γ in increasing raw order), so the current revision refuses
  those swaps by name and the "residual" wording was removed.  No mutation battery
  was run for the R27 guards.  The **mapping-level** rejections listed above — count,
  index range/repeat/position, column range/repeat/position, the time-reversal
  cross-check and the rewritten raw rotation/translation — each have a focused unit
  test asserting that exact variant; a semantically inert mutation is not
  discriminated.  The raw-surface guard
  (`TargetRawOperationsDiffer` / `TargetRawOperationsUnavailable`) has **no** dynamic
  error witness: it is supported by the source-level fence and by the positive
  all-1,651-UNI surface equality check, not by a test that provokes it.  The other
  rejection branches (`TargetClassUncovered`, `TargetTranslationOffGrid`, the
  character-availability branches) are likewise not dynamically witnessed here, and
  nothing in this document claims otherwise.
