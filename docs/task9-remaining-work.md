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

Probed 2026-09-22: every scalar field of the pinned little table is identical for
`DT1..DT4` (`little_irr_type`, `little_irr_cc`, `little_irr_lif`,
`little_irr_invar3`, `little_irr_full_dim`), but `little_irr_kov` is **not** the
identity: SG 202's DT slot reads `[1, 2, 4, 3]`, i.e. ISOTROPY's `DT3`/`DT4` are
Kovalev's irreps 4 and 3 respectively (SG 225 reads `[1, 3, 4, 2, 5]`).  So a
third route exists: look the two characters up in Kovalev's tables (the suite's
own ISO-KOV mapping page maps Kovalev onto CDML), which would settle the swap
without touching the undecoded matrix block.

**Best lead (round 37).** The engine already exposes the *little* representation of
a table irrep: `IrrepRecord::ordinary_scalar_selected_arm_block_trace()` (used by
`trivial_child_record`), whose values are exactly the selected arm's block trace,
i.e. the little rep on the little group's operations.  Combine that with the
official compatibility labels already used for the Γ system:

* SG 202's X point gives `X3+ -> DT3`, `X2+ -> DT4` (multiplicity 1 each);
* the little rep of `X3+` is one-dimensional there (`full_dim / star = 3 / 3`), so
  its block trace **is** the character of the line irrep `DT3` on the DT little
  group, evaluated at `k = X`;
* dividing out the Bloch phase `exp(-2 pi i k_X . t)` (the X point is at
  `alpha = 1/2` on the DT line) gives the Γ-point characters, i.e. exactly the two
  tables the frozen file is missing.

So the eight unresolved sources can be closed without `data_images.txt` and
without decoding `little_irr_full_matrices`: dump the selected-arm block traces of
the X irreps (`X3+`, `X2-` for `DT3`; `X2+`, `X1-` for `DT4`), read the phase, and
extend `freeze_w_little_characters.py` with that second route.  The remaining work
is a small read-only Rust dump of the block traces plus the op order they follow.

Confirmed 2026-09-22 with `examples/dump_little_trace.rs` (added, read-only):
`SG 202 X3+: little dim 1 over 96 operations`, i.e. the block trace really is the
one-dimensional little rep on the X little group, including its lattice-translation
cosets.  Note the loop matches every record whose label is `X3+`, so look at the
header line of the record at the X point (dim 1) rather than at the tail of a
later record (which can be 3-dimensional and shows zero traces).

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
