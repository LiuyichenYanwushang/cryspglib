# The scalar Γ magnetic input contract (R25 card 1)

`irrep::magnetic_input` turns one verified magnetic embedding into the exact
corepresentation a scalar Γ magnetic subduction starts from.  It is the input
layer of R10 (`docs/subduction-next-milestones.md`, "先实现 Γ 磁共表示分导"), not
the subduction.

## The model

Write `P` for the ordinary parent space group and `P × {1, Θ}` for its **grey**
group, with `Θ` the pure time reversal.  A magnetic subgroup `M` (Type I, III or
IV) is embedded into that grey parent by the verified map `φ`; `M` is **not**
required to contain `Θ`.  On a fixed **real one-dimensional** basis of the
parent's Γ space the parent's action is `ρ(p, Θ^t) = χ(p)·η^t·K^t`, so `Θ` acts
as

```text
Θ = η·K,      η = ±1,      K = complex conjugation
Θ² = η²·K² = +1
```

Let `χ` be the parent irrep's character — this contract covers the **real `±1`**
case explicitly, not scalar characters in general — and let `π(g)` be the parent
operation a magnetic class of `M` maps to under `φ`.  The model is the
**restriction of `ρ` to `φ(M)`**:

```text
A(g) = χ(π(g))          for unitary g
A(g) = η · χ(π(g))      for antiunitary g
```

Both halves are **computed**: `χ` from the parent's tabulated row, the images from
the embedding.  The caller states the embedding and `η`.

`η` is a datum of the **parent** action, so no carrier class is chosen: `A(g)` is
that action evaluated at `g`, and which listed class of `M` happens to be a pure
`Θ` (if any) is irrelevant.  The grey parent always contains `Θ`; a magnetic
subgroup of Type III/IV need not, and its antiunitary classes inherit `η` through
`φ`.  Both signs of `η` are legal and flipping `η` flips exactly the antiunitary
classes of the restriction.

## What is checked

`A` must be the restriction of a corepresentation of the grey parent:

```text
A(g) · [ conj(A(h))  if t(g) else  A(h) ] = A(g·h)
```

for every ordered pair of classes, with `t(g)` the time-reversal flag of the
**left** factor; the unitary branch keeps `A(h)` unchanged.  Every value is `±1`
and `conj(±1) = ±1`, so the conjugation is the identity on the values while the
branch is structural.

As the special case `g = h = Θ`, `A(Θ)² = A(Θ²) = +1`, because `Θ² = +1` is built
into the parent action.  **Time parity is not the square**: a time-odd scalar
(`η = -1`) is reported as `TimeParity::Odd` while `time_reversal_square()` stays
`+1`, and the model never reports a time-odd scalar as `T² = -1`.

`η` is a factor of the coset phases, not a per-class table: the antiunitary class
`h·Θ` carries `η·χ(π(h))`, so the coset values vary with the spatial part while
`η` itself is one number.

## Verified rather than trusted

| claim | how it is checked |
|---|---|
| the operation set is a group | `magnetic_embedding::verify_group` |
| it is the database's set | `magnetic_embedding::verify_against_database` |
| the record is the one named | the record is the row at the given index of the resolved condensing irrep's own magnetic isotropy table, and the geometry and operation set are derived from that record (UNI, basis, origin); the caller hands in no record |
| the caller's `(map, shift)` embeds the group | `magnetic_embedding::verify_embedding` — the public verifier, not a candidate family |
| the record-lattice requirement | `RequireExact` recomputed by `realises_record_lattice_checked`, which returns `Result`, so an exact-arithmetic limit is an error and never a silent `false` |
| the character is the parent's | read from the source row at the mapped rotation and translation, required real `±1`, identity `+1` |
| the table is a corepresentation | the product law on every ordered pair of classes |

Two identities are kept apart on purpose:

* the **source identity** (`SourceIdentity`: parent, label, convention) names the
  subduced parent irrep, which must be scalar, non-spinor and at Γ;
* the **record identity** (`RecordIdentity`: parent, condensing label, convention,
  row index) names the magnetic isotropy row.  The condensing irrep may sit at a
  line or plane; the source may not.

Each identity resolves under **its own** convention, and the condensing record is
then obtained by entering the isotropy query with the resolved irrep's canonical
CDML label instead of looking the caller's spelling up a second time, which is why
a Bradley-Cracknell request works.  The labels the model reports are the canonical
CDML ones, so `source_label()`/`source_convention()` and
`record_label()`/`record_convention()` each resolve back to the selected irrep.
`input_convention()` reports the convention the caller used to name the **source**
label only; no input string and no record input convention are retained.

The record-lattice requirement is a policy, not a guessed boolean.
`ContainmentOnly` is the proved property that the request always has: embedding in
the grey parent, checked by `verify_embedding`.  `RequireExact` adds `A·L_M =
L_record`, recomputed with `realises_record_lattice_checked`.  For `P-1'`, the map
`A = diag(2,1,1)` embeds every class but does not realise the record lattice, so it
is accepted under the first policy and refused under the second, and the accepted
model carries the recomputed verdict either way.

The returned `ScalarGammaModel` has private fields, no constructor and no `Clone`,
so the only path to one is `verify_scalar_gamma_input`; it carries the verified
`map`/`shift`, both identities, the policy and the action table.

A Type I (purely unitary) magnetic subgroup needs no `η` — its restriction never
uses one — but a stated `η` is preserved as the parent action's datum, and
`time_reversal_square()` is `None` because there is no antiunitary class to square.

## Scope fences

* **Real, non-spinor, one-dimensional Γ source.**  A complex character, a spinor
  record, a source `k ≠ Γ` and a source dimension other than 1 each have their own
  `MagneticInputError`; none is reported as another.  The isotropy tables contain
  no complex scalar Γ row, so the complex refusal is a contract on future data — a
  test scans and says so.
* **No uniqueness gate.**  Several verified embeddings of one record are
  legitimate; the caller selects one, it is checked, and the model records it.  The
  action table is a function of `χ` and the mapped classes; the tests pin that two
  embeddings of the same witness agree, and do not claim it as a general theorem.
* **No subduction, no matrix machinery, no dependency.**  The ordinary entry
  points are untouched.

## The one measured boundary: table alignment not established

The parent's operation table (`query::symmetry_operations_of`) and the source
character row (the PIR operation storage behind
`ordinary_scalar_selected_arm_block_trace`) are generated in their own settings.
A historical measurement (round 25, before this module existed) compared them
operation by operation modulo the parent's lattice and found them equal for 828 of
the 980 scalar Γ magnetic isotropy rows and different for the remaining rows — for
example SG 48 `GM1+`, whose parent table has `(1̄|(1/2,1/2,1/2))` where the row has
`(1̄|(1/2,1/2,0))`.  The rotation sets are equal in all 980 rows.  That number is a
measurement of the two tables' agreement, not a coverage or success rate of this
interface.

This round has **not established** how to align the two tables, so a row whose
mapped operations are not in the source row is refused by name
(`MagneticInputError::CharacterRowSettingMismatch`) rather than paired by
guesswork.  No alignment transform was searched for and none was refuted: the
refusal says what this contract does, not what is impossible.  SG 48 `GM1+` is
pinned as the regression for that refusal.
