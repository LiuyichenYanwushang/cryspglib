# Task 9 setting pipeline (closed for the ordinary identity table)

The full-subduction engine maps every isotropy record's subgroup into the parent
frame of `SG_DATA_HALL`.  Two conventions are not in the stored tables and have
to be *derived* from the official ISOTROPY program before they can be frozen:

* `U` — the change of basis `W . P_parent . (P_sub . B_oracle)^-1` implied by the
  basis the program prints for the record,
* `child_shift` — the child-frame origin correction that reconciles the recorded
  ITA setting with the game frame the engine maps into.

Every derived convention is handed to the engine's own validation
(`examples/probe_subduction_settings.rs` -> `SubgroupEmbedding::probe_embedding`);
nothing is frozen on the strength of a re-implementation of that check.

## Pipeline

```bash
cd cryspglib
mkdir -p target/task9 && cd target/task9   # scripts write their evidence here

# 1. census: one official query per (parent, irrep), SET I ALL OR 1
python3 ../../../scripts/audit_subduction_settings.py \
    --output census_or1.jsonl --workers 8

# 2. which of those conventions does the engine already accept?
python3 ../../../scripts/task9/probe_census.py census_or1.jsonl --out rejected.jsonl

# 3. derive the child shift the official program implies (OR 1 vs OR 2), and keep
#    only the shifts the engine accepts
python3 ../../../scripts/task9/derive_shift_oracle.py \
    --census census_or1.jsonl --rejected rejected.jsonl --out derived_shift2.json

# 4. fall back to an engine-validated search where the derivation is rejected
python3 ../../../scripts/task9/solve_record_delta.py rejected.jsonl --out record_delta.json

# 5. the 195 records whose legacy compound spelling prints an empty table are
#    addressed through the compact little_irr_full_label spelling
python3 ../../../scripts/task9/fix_empty_labels.py \
    --census census_or1.jsonl --out empty_fixed.jsonl

# 6. emit `src/irrep/subduction_settings_data.rs`
python3 ../../../scripts/task9/build_table.py \
    --census census_or1.jsonl --shifts record_delta.json \
    --derived derived_shift2.json --legacy legacy_shifts.json \
    --empty empty_fixed.jsonl --out ../../../src/irrep/subduction_settings_data.rs
```

## Second pass: every record embedded (2026-09-21, second session)

Two changes closed the embedding table:

* the engine now takes the subgroup lattice from the candidate's own map
  (`U^-1 W P_parent`) instead of `W . P_parent`; the two differ exactly when `U`
  is not unimodular, which is why the fractional monoclinic conventions were
  rejected -- see commit "take the subgroup lattice from the embedding";
* the child shift is derived by **two** frame routes and each is probed
  separately (`scripts/task9/derive_shift_both.py`):
  `-(B^T)^-1 (o_OR1 - o_OR2)` for the ITA origin-choice difference (2,133
  records) and `(B^T)^-1 (o_printed - o_stored)` for the recorded-frame
  difference (172 more), the second being what the last monoclinic and
  compact-label records needed.

Result: **15,239 / 15,239 records embed** (0 rejected), stored identity
positives 93,720 passed / 391 mismatch, Gamma Frobenius 1,895/1,895,
absent_positive 380, hard_failures 771.

### The remaining 771 rows are *not* a setting problem

`examples/trace_embedding.rs` prints one record's subgroup lattice, transform and
every mapped coset representative.  For ordinal 12471 (SG 221 `X3+` `P2` -> #125,
where the engine reports 0 for the stored frequency 1), all 16 mapped
representatives are **exactly** the official `SHOW ELEMENTS` operations carried
into the parent frame modulo the parent lattice, and all 16 rotations are
distinct.  The embedding is therefore right and the disagreement is in the
trivial-multiplicity evaluation of the subduced block -- `parent_character_of` /
`solve_character_block` in `src/irrep/subduction.rs` -- not in
`subduction_settings_data.rs`.

Every failing ordinal has size > 1 (105 of size 2, 55 of size 4, 20 of size 8,
16 of size 32, 2 of size 6); the size-1 Gamma records are all covered by the
Frobenius gate, which passes 1,895/1,895.  That points at the supercell/folded-k
path: cosets whose parent-frame representative carries a lattice translation are
where the Bloch phase `chi(t + L) = chi(t) exp(2 pi i k . L)` has to be applied.

### Round 3: what the diagnostics established

* The engine's **folding convention matches the stored table**.  Ordinal 58
  (SG 5 `V1` -> #1) folds to two child k-points, `GM1` and `X1`, which differ by a
  vector of the subgroup's own reciprocal lattice; the stored table also records
  the trivial content as one (`GM1`), not the folded sum, so the "unfolded child
  k-point" description is the shared convention and the audit's reading is right.
* For ordinal 12471 the engine is **self-consistent**: its subduced character
  reproduces its own decomposition (`reconstruction` max error 2.4e-16), the
  totals agree with the stored table through the Frobenius index (both sides sum
  to 12), and the disagreement is purely *which* irreps of the star carry the
  trivial content -- engine `{X1-, X2-, M1+}` vs stored `{X3+, X4+, M4+}`.
  The stored rows for that record carry `domain` numbers 1 and 4, so part of the
  disagreement is a domain (conjugate-subgroup) mixture in the comparison, but
  the domain-1 subset still disagrees.
* `scan_character_norm` reports 3,030 `<chi,chi> != sum mult^2` rows; its doc
  records why that is a *granularity* artifact of the unfolded child-k
  convention rather than a wrong character, and is therefore not a bug gate as
  written.

**178 of the 198 failing ordinals fail on their own condensing irrep**: the stored
table says the condensing irrep contains the trivial representation of its
isotropy subgroup once, the engine says zero.  That needs no oracle to be a bug --
a subgroup is the isotropy subgroup of a direction precisely because the
condensing irrep contains its trivial representation -- so the full-star adapter
(`subduce_full_star_with_embedding`), not the frozen table, is where the next
round has to look.  A permanent gate should assert this condition directly
instead of only comparing against the stored frequencies.

Still open, reported explicitly and not counted as covered: 14,713 probes whose
folded child k is not in the shipped discrete table, 5,756 other-wave-vector
entries without resolved k parameters.

## Fourth pass: the engine's placement was not the official one (2026-09-21)

Three fixes closed 1,698 of the 1,758 consistency rows:

1. **The subgroup lattice has to follow the accepted map.**  `SubgroupEmbedding`
   kept `W . P_parent` in its lattice field while validating against the lattice
   the candidate's own affine map produces; with a fractional `U` those differ.
2. **Every record has to sit at the origin the official program prints.**
   `delta = (B^T)^-1 (o_printed - o_stored)` was only derived for records the
   engine *rejected*, but two placements can both be valid subgroups, so the
   engine's validation cannot tell the official one from the other.
   `--recorded` applies it to all 1,057 records whose stored origin differs from
   the printed one (523 engine-accepted), and it is oracle-free.
3. **The child's Hall origin choice has to be derived even where `delta = 0`
   validates.**  Re-deriving both frame routes for the 107 records the audit
   still flagged gave 103 accepted shifts (97 origin-choice, 6 recorded-frame).

| audit item | before the pass | after |
|---|---|---|
| embeddings | 15,239 / 15,239 | 15,239 / 15,239 |
| identity positives passed | 93,720 | **94,081** |
| identity mismatch | 391 | **30** |
| absent computed positive | 380 | **30** |
| hard failures | 771 | **60** |
| Gamma Frobenius | 1,895 / 1,895 | 1,895 / 1,895 |

Minimal witnesses, both fixed: ordinal 27 (SG 3 `A2` -> #3) where `delta = 0`
placed the subgroup a half-cell away and moved the trivial content onto `GM2`;
ordinal 1197 (SG 48 `R1+` -> #70) where the child's Hall row lives in the other
ITA origin choice and the engine reported `GM1-`.

The 60 remaining rows are 16 records of SG 5/12/15/67/68 (plus 13090/13106) whose
stored origin differs from the printed one by a parent **cell choice**;
`SET I <sg> CELL n` changes the direction set of those parents, so the records no
longer line up one-to-one and the parent setting needs another route.

## Sixth pass: the ordinary identity table is closed (2026-09-22)

The audit's 14,713 `uncomputed_missing_data` probes (160 of them stored
positives) were *not* a setting problem either.  `examples/trace_subduction` now
dumps the folded child stars when the full decomposition fails, and ordinal
10027 (SG 196 `W1` `P2` -> #24) shows why they are answerable:

```text
folded child stars: 3
  child star 0: q = (-2,-1/2,0), (2,1/2,0)          -> not Gamma, no #24 data
  child star 1: q = (0,-1,1/2), (0,1,-1/2)          -> not Gamma, no #24 data
  child star 2: q = (1,0,-1)                        -> Gamma, arms [0, 2]
```

Only the third star can carry the subgroup's trivial representation (a child
lattice translation `t` acts on a representation at `q` as `exp(-2 pi i q.t)`,
the trivial representation acts as `1`), so
`trivial_content_with_embedding` evaluates the Gamma star with the same
`build_block` stage and skips the other two exactly.  It returns 1 for `W1`
(stored frequency 1), 0 for the unlisted `W2`, and 0 for `L1`, whose stars carry
no Gamma point at all.

| audit item | fifth pass | sixth pass |
|---|---|---|
| embeddings | 15,239 / 15,239 | 15,239 / 15,239 |
| probes with an exact result | 351,547 of 366,260 | **366,260** (351,547 full + 14,713 identity-only) |
| stored identity positives | 94,111 passed | **94,271 passed**, 0 mismatch, 0 false positive |
| absent entries computed positive | 0 | 0 (271,989 zeros) |
| Gamma Frobenius | 1,895 / 1,895 | 1,895 / 1,895 |
| hard failures | 0 | 0 |
| verdict | `clean ... incomplete_categories=20629` | `--require-complete` exit 0, `VERDICT complete scope=global` |

The remaining `incomplete_categories` were the 5,756 `isotropy_w_subduce_*` rows.
Those name 73 parent irreps **at other wave vectors**; the pinned *irrep* table
stores them as label, space group, dimension and type only -- no k vector and no
character row -- so it cannot answer them.
`scripts/check_other_wave_vector_rows.py` makes that executable (11 offline
tests) and checks everything that table does pin; the audit reports the rows in
its own `w_scope` line and gates them with `--require-w-complete`, which exits 2
while they are uncomputed.  Their **little-group** tables are not missing, though:
`data_little.txt` carries all 73 sources (`little_irr_full_label` /
`little_irr_space_group` / `little_irr_full_dim`, dimensions equal to
`irrep_w_dimension`, e.g. SG 225 `DT1-DT4` = 6, `DT5` = `SM1-SM4` = 12), and the
little table's wave-vector decode shows every one of them to be a **parameterized
line**: `little_k` is 14 Bravais lattices x 27 k slots x 16 ints = a base point
plus up to three free directions per slot (the aP block reproduces
`data_space.txt`'s k points exactly: `Z=(0,0,1)/2 ... T=(0,1,1)/2`, and the
general position carries three directions).  All 73 sources are `k = Gamma + t*v`
with exactly one free parameter (cF: `DT=(1,0,1)`, `SM=(1,1,2)`; cI:
`DT=(1,-1,1)`, `SM=(0,0,1)`), so there is no numeric k to fold.  The remaining
work is the little-group character tables along those lines, then comparing the
5,756 stored frequencies line by line.  `scripts/check_other_wave_vector_rows.py`
asserts all of this from the pinned archive (15 offline tests).

## State after the first pass (2026-09-21)

The frozen table now addresses **all 15,239** pinned isotropy records.  The
production audit (`examples/audit_irrep_subduction`) reports:

| item | result |
|---|---|
| embeddings | 15,144 ok / 95 rejected |
| scalar probes | 348,932 complete / 14,713 missing child data / 760 engine errors |
| stored identity positives | 93,350 passed / 434 mismatch / 245 unavailable |
| Γ Frobenius | 1,895 / 1,895 |
| absent entries computed positive | 422 |
| other-wave-vector entries | 5,756 still without resolved k parameters |

`VERDICT inconsistent ... hard_failures=1698`, so this is progress, not
acceptance.

## Measured negative result: the monoclinic parents are *not* a setting problem

The 95 rejected embeddings concentrate in the monoclinic parents SG 3–15, whose
stored basis/origin do not match the default query frame (`SET I ALL OR 1`) for
about a third of their records.  Sweeping the parent setting
(`scripts/task9/sweep_mono.py`: `OR1`/`OR2` × `AXIS b|c` × `CELL 1|2|3`, 18
variants) and asking the engine about each variant's derived `U` gives:

| | |
|---|---|
| best variant | `OR1` for all 13 parents (ties broken by order) |
| records accepted | 172 / 262 |
| accepted set | **identical for every variant** (`variant_sets.py`, SG 14: all 18 variants accept the same 14 ordinals) |

So re-querying in the recorded ITA setting does not move the engine's answer:
the 90 rejected monoclinic records need a *child shift*, not a different parent
frame, and the `OR1`/`OR2` origin difference does not produce one for them
(the two settings print different cells, so `-(B^T)^-1 (o_OR1 - o_OR2)` is not
the frame difference).  Next entry point: derive the monoclinic child shift from
`SHOW ELEMENTS` (the generator's `origin_choice_shift` route) instead of from the
printed origins, and cross-check it against the stored identity frequencies.
