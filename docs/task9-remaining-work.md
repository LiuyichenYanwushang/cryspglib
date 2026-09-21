# Task 9 remaining work (handoff, 2026-09-22)

Task 9's ordinary identity-subduction table is closed by the engine; the 5,756
`other_wave_vector_subduction` rows are verified against the official program but
not yet *computed* by the engine. This note is the shortest complete description
of what is left, so the work can resume without re-deriving anything.

## Where things stand

| Item | Evidence |
|---|---|
| Ordinary table (15,239 records / 94,271 positives / 366,260 probes / Γ Frobenius 1,895) | `audit_irrep_subduction --require-complete` exits 0, `VERDICT complete scope=global`; re-run 2026-09-22 in 446.7 s |
| w rows (1,006 records / 5,756 rows) | `scripts/verify_w_subduction_oracle.py`: 28 groups, 1,150 records, 5,756 oracle rows = 5,756 pinned rows, 0 mismatches |
| Frozen little-group characters | `src/irrep/w_little_characters_data.rs` (65 of 73 sources), pinned by `tests/w_little_characters.rs` and `scripts/test_frozen_w_little_characters.py` |
| Coverage split | `scripts/check_other_wave_vector_rows.py` prints `rows_with_frozen_table=5408 rows_blocked=348` |
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

## Step 2 — the eight unresolved sources (348 rows)

SG 202/203/209/210 `DT3`/`DT4`. The Γ compatibility data fixes only their sum;
every other official channel is now excluded:

* `SHOW CHARACTER` at another special point of the line prints the **full** irrep
  (star included), not the little irrep, so the compat rows cannot be used there;
* the pinned `little_subduce` blocks of `DT3` and `DT4` are identical row for row;
* `SHOW KERNEL` + `DISPLAY IRREP` prints nothing for parametric k either.

Remaining routes: `data_images.txt` (each image record carries the point-group
representation, hence the point characters, of a little irrep) or decoding the
pinned `little_irr_full_matrices` block of `data_little.txt` (rounds 7-8 recorded
what is known about the pointer chains; the value encoding is still open).

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
