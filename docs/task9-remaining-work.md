# Task 9 remaining work (handoff, 2026-09-22, updated round 41)

Task 9's ordinary identity-subduction table is closed by the engine, and the
pinned data behind the 5,756 `other_wave_vector_subduction` rows is now complete:
all 73 parametric-k sources have frozen little-group character tables.  The rows
themselves are still verified against the official program rather than computed
by the engine.  This note is the shortest complete description of what is left,
so the work can resume without re-deriving anything.

## Where things stand

| Item | Evidence |
|---|---|
| Ordinary table (15,239 records / 94,271 positives / 366,260 probes / Γ Frobenius 1,895) | `audit_irrep_subduction --require-complete` exits 0, `VERDICT complete scope=global`; re-run 2026-09-22 in 446.7 s |
| w rows (1,006 records / 5,756 rows) | `scripts/verify_w_subduction_oracle.py`: 28 groups, 1,150 records, 5,756 oracle rows = 5,756 pinned rows, 0 mismatches |
| w row characters | `cryspglib::irrep::w_little_characters_data::W_LITTLE_CHARACTERS` (73 tables) with `W_LITTLE_CHARACTERS_UNRESOLVED` empty |
| Frozen little-group characters | `src/irrep/w_little_characters_data.rs` (**73 of 73** sources; 65 from the Γ compatibility rows, 8 from the little-cogroup pair route), pinned by `tests/w_little_characters.rs` and `scripts/test_frozen_w_little_characters.py`; `--rust` regenerates the file byte-identically |
| Coverage split | `scripts/check_other_wave_vector_rows.py` prints `rows_with_frozen_table=5756 rows_blocked=0 unresolved_sources=[]` |
| Engine gate | `--require-w-complete` exits 2 (`w_scope: rows=5756 computed=0`) |

## Step 1 — compute the 5,408 rows that have frozen characters

The engine already implements everything except the arm character source
(`src/irrep/subduction_star_decompose.rs`, `trivial_content_with_embedding`):
arm enumeration, the fold-onto-child-Γ test, the child cell/reciprocal lattice
and the block assembly are validated. Add a sibling entry point that keeps that
machinery and evaluates the arm characters from
`w_little_characters_data.rs` instead of the CIR matrices:

```text
χ_{W'}(R, t) = exp(-2π i α (v · t)) · D(R),        α generic (e.g. 3/8)
frequency    = Σ_{H-orbits of arms folding to child Γ} mult(trivial_{H ∩ sG_Ls⁻¹}, W'^s)
```

* `v` is the table's `direction` (parent **conventional** reciprocal basis);
  `t` is the operation's translation in the parent conventional cell;
  `D(R)` is the table's character for that rotation (identity rotation first).
* The averaging over the child's lattice cosets is *equivalent* to the existing
  fold-onto-Γ test — do not re-derive it (the prototype that did produced 12
  folding arms where the pinned frequency is 4, see
  `docs/subduction-audit.md`).
* Sanity anchors: SG 196 `6D1` (P1 child) → `DT1 = 6`, `DT2 = 6`, `SM1 = 12`;
  `w_scope` must move from `computed=0` to `computed=5408`.

Exact integration point (checked 2026-09-22): `build_block` in
`src/irrep/subduction_star_decompose.rs` builds the q-block characters with

```rust
parent_characters.push(star.q_block_character(point.arm_indices(), operation)?);
```

so the only coupling to the CIR data is `ScalarStar::q_block_character`. A
variant that answers from `w_little_characters_data` (rotation lookup in the
table plus the Bloch phase) plugged in there gives the w frequency through the
existing, validated code path; `little_group_operations`, `identity_position`,
`select_representative` and the child-cell handling stay untouched. Sanity gates
to keep: `chi_q(E) = q_block_dimension` (already asserted in `build_block`) and
the SG 196 `6D1` anchor above.

If instead the character sum is written directly (rather than through
`build_block`), these public accessors supply every input — checked 2026-09-22 in
`src/irrep/subduction.rs`:

| input | accessor |
|---|---|
| parent operations | `strict_sg_hall_ops(sg)?.operations()` |
| child operations modulo `L_H` (point-group order many) | `embedding.representatives()` |
| child operations in the parent frame | `embedding.operations()` |
| parent / subgroup translation lattices | `embedding.parent_lattice()` / `embedding.subgroup_lattice()` |

The little group is `{(R,t) : R v = v}` (rotation-only test); the only extra
averaging is over `L_H/L_parent`; the phase is `exp(-2 i pi alpha (v . t))` with a
generic rational `alpha`.  Beware: the naive sum *without* that averaging
overshoots supercell children (measured `6D1` = 12 ✓ but `4D1` = 12 where the
pinned value is 4), which is why the folding path of `build_block` is the safer
first target.

## Step 2 — the eight sources the Γ rows cannot separate (closed)

SG 202/203/209/210 `DT3`/`DT4` (348 rows) are frozen now.  The Γ compatibility
system sees the pair only through `DT3 + DT4`, which it *does* determine
(`cogroup_pair_route` solves that functional with the dual system); the split
comes from the pinned source ordering, and the route refuses to emit anything
unless that ordering is confirmed by the two sources the Γ rows fix by
themselves.

What the domain is.  The `DT` line's little cogroup is abelian of order four:

| parents | cogroup | `DT1`, `DT2` (= `A1`, `A2`) | `DT3`, `DT4` (the open pair) |
|---|---|---|---|
| SG 202, 203 | `C2v = {E, C2, m1, m2}` | `(1,1,1,1)`, `(1,1,-1,-1)` | `(1,-1,1,-1)`, `(1,-1,-1,1)` |
| SG 209, 210 | `C4 = {E, C2, R, R^3}` | `(1,1,1,1)`, `(1,1,-1,-1)` | `(1,-1,i,-i)`, `(1,-1,-i,i)` |

(Characters over the four little-group operations in the order `SHOW ELEMENTS`
prints them, identity first.  For the two `C4` parents the open pair is the
conjugate pair, so the frozen table stores `(real, imaginary)` integer parts --
`LittleOperation::character: [i32; 2]` -- and the sum `DT3 + DT4` is `(2,-2,0,0)`
in every one of the four parents.)

Evidence, all reproducible:

* `SG {sg} k DT` pair sum, from the pinned Γ system: `(2,-2,0,0)`;
* per-parent ordering check: the two determined sources come out `(1,1,1,1)` and
  `(1,1,-1,-1)`, i.e. `A1` and `A2`, in all four parents;
* an independent oracle route agrees wherever it reaches: the extra special point
  of the line (`X`) fixes `DT3` and `DT4` separately for SG 203/209/210, and its
  values are exactly the ones above (SG 202's `X` rows are inconsistent for
  `DT4`, which is why the ordering route exists);
* `tests/w_little_characters.rs::archived_cir_characters_confirm_the_sg202_dt_pairing`
  derives the SG 202 tables from the archived CIR characters of the compatible
  X-point irreps (`X3+ -> DT3`, `X2+ -> DT4`, `X2- -> DT3`, `X3- -> DT4`,
  `X1+/X4- -> DT1`, `X4+/X1- -> DT2`), dividing out the Bloch phase
  `exp(-2 pi i k_X . t)`: every one of the 4 rotations of every one of the eight
  compatibilities matches the frozen character.

The route is falsifiable and already exercised: SG 202 has 6 pinned records whose
`DT3` and `DT4` frequencies differ and SG 209 has 4, so the engine step below
will confirm or refute the ordering on real rows (SG 203/210 have none).

Rejected/closed routes, for the record: `SHOW CHARACTER` at another special point
prints the *full* irrep (star included), not the little irrep; `SHOW KERNEL`,
`SHOW MATRIX`, `VALUE KVALUE` on a parametric k domain print nothing usable (the
program refuses to select a parametric k irrep: "parameters not selected for k
vector"); the pinned `little_subduce` blocks are identical row for row for the
pair; `little_irr_full_matrices` is still undecoded (rounds 7-8).

## Step 3 — close the loop

* extend `examples/audit_irrep_subduction.rs` so the w rows are computed and
  compared like the ordinary rows (both directions: listed rows and zeros);
* make `--require-w-complete` pass and update the coverage statement in
  `docs/subduction-audit.md` from two tracks to one;
* keep the gates green: `cargo test --release -p cryspglib --tests` / `--doc`,
  `cargo clippy -p cryspglib --all-targets --release -- -D warnings`, and the
  five offline Python suites plus the three live oracles
  (`verify_isotropy_oracle.py`, `verify_w_subduction_oracle.py`,
  `check_other_wave_vector_rows.py`).
