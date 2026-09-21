# Task 9 setting pipeline (work in progress)

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
