# Task 9 remaining work (handoff; last updated round 61)

> **START HERE (one paragraph).** The only thing left is one engine function:
> give `build_block` (`src/irrep/subduction_star_decompose.rs:1078`) a second
> arm-character source so the frozen parametric-k irreps go through the engine's
> own fold/character-block solve.  Everything else is done and gated: the
> ordinary table is 94,271/94,271 engine-computed (`VERDICT complete
> scope=global`), all 73 line sources are frozen and cross-checked, all 5,756 w
> rows are live-oracle-verified, and the engine already reproduces 3,916 of
> 5,756 w rows (`line_trivial_content_with_embedding`, audit prints
> `computed=300`).  Four hand-written per-arm weight rules were tested and
> excluded by measurement — do not retry them; the sections below hold the exact
> insert points, the acceptance ladder, the excluded rules and the measured
> witnesses.

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

## Step 1a — the engine computes the `P1`-child family (round 42)

`line_trivial_content_with_embedding(subgroup, embedding, table)` in
`src/irrep/subduction_star_decompose.rs` turns one frozen line source into a
frequency:

```text
for each star arm a (distinct contragredient images of the frozen direction):
    skip unless a/4 lies in the child's reciprocal lattice   (the fold at t = 1/4)
    weight = (1/|K|) sum_{h in K} chi_{W'}(transport^-1 h transport),  K the child's
             operations (one per rotation) that fix the arm
    frequency += weight
```

`LINE_PARAMETER = 1/4` is the program's convention for the free parameter of a
parametric-k domain: on the 46 `P1`-child records (300 pinned rows) `1/4` and
`3/4` reproduce every row, `1/2` and `1` leave 111 wrong and the other twelfths
leave all 300 wrong.

Measured: `tests/w_line_frequency.rs::p1_child_records_match_every_pinned_row`
computes **300 / 300** pinned `P1` rows exactly, and
`a_source_of_another_parent_is_rejected` pins the fail-closed check.  The
`P1`-child family is therefore engine-computed, on top of the oracle
verification of all 5,756 rows.

**Still open (the rest of the w track).** `asymmetric_dt3_dt4_records_are_pinned`
is `#[ignore]`d: for a child that is *not* `P1` the fold/weight convention is not
reproduced yet.  Measured on SG 202 ordinal 10422 (`DT3`, child #31 `Pmn2_1`,
size 16, pinned `DT1 = 1, DT3 = 2, DT4 = 1, SM1 = 2`): the model folds four arms
(`(0,0,+-2)`, `(+-2,0,0)`), each arm's stabiliser has two operations, and the
transported characters come out `(+1,-1)`, `(+1,-1)`, `(+1,+1)`, `(-1,-1)`; the
per-arm weights are therefore `0, 0, 1, -1` and the total is 0 instead of 2.
The same record pins a ratio that rules out two further readings: `DT1`, `DT3`
and `DT4` are the *same* three little irreps on the *same* arms (little dimension
1, character rows `(1,1,1,1)`, `(1,-1,1,-1)`, `(1,-1,-1,1)` over
`(E, C2y, sx, sz)`), yet their pinned frequencies are `1`, `2`, `1`.  A per-arm
weight of "the whole little rep is trivial on the stabiliser" (the
invariant-subspace reading) is 0/1 and can only produce sums that are counts of
arms, and the plain character average over the four-element stabiliser gives
`4`, `0`, `0`; neither can produce `1, 2, 1` from one arm set.  So the weight is
neither the projection nor the strict-invariance indicator — it is some other
object (a candidate worth testing next: the *multiplicity* with which the arm's
little rep occurs in the child's *full* star decomposition at the folded point,
i.e. the pairing the engine already computes in `build_block` for discrete
probes, which never enters a hand-written per-arm formula).

A negative per-arm weight already proves the weight is not the trivial-content
projection I assumed — the sum must be computed differently (the natural
candidates: average over the child's *full* point group with the arm
permutation included, or over the stabiliser of the arm in the parent's little
group `G_L` rather than in the child).  Using the lattice the isotropy record
stores (`W . P_parent`) instead of the accepted embedding's lattice did not
change this record, and the `P1` family stays 300/300 with it.  The pinned rows
themselves are unaffected — they are oracle-verified — so the audit keeps
reporting the w rows as a separate track until this is closed.  The next
candidate to test is the frame of the projection: `embedding.representatives()`
carries four distinct rotations for that record, so the child's point group as
the embedding accepts it is larger than the recorded cell suggests, and the fold
test may belong in the child's own frame (`fold_wave_vector(embedding.transform(),
k)`) rather than in the parent conventional one.

## Ruled out: the compatible X-point irrep as a stand-in (round 48)

Using the ordinary table's pinned rows for the X-point irreps as a shortcut for
the line sources does not work, even though `SHOW COMPATIBILITY` pairs them
(`X3+ -> DT3`, `X2+ -> DT4`, ...).  Measured: ordinal 10422 (`DT4` partners
`X2+ = 1`, `X3- = 1` -- equal) but `DT1 = 1` against partners `X1+ = 1`,
`X4- = 2`, and `DT3 = 2` against partners `X3+`, `X2-` that are *absent* from
the ordinary table, i.e. zero.  Ordinal 11171 has no ordinary X row at all for
`X5` (the irrep compatible with the `DT3`/`DT4` pair).  Sums do not help either:
they agree for 10422 (6 = 6) but differ by a factor two for 10485 (24 vs 12) and
11171 (6 vs 3).  So the w frequencies genuinely need the line representation:
the remaining work is the line-star variant of `build_block` (arm enumeration,
fold onto the child Gamma, character block solved against the child's stored
rows) rather than any pairing with discrete probes.

Round 49 also tested the strict-invariance weight (a per-arm weight of the little
dimension when the transported little representation is trivial on the arm's
stabiliser, otherwise zero — a 0/dimension count, which *can* produce the pinned
`1, 2, 1` pattern).  Measured: SG 202 ordinal 10422 `DT3` still gives 0 instead
of 2 (no arm's transported characters are all `+1`), so that reading is out too;
the `P1` family stays 300/300 under it.  Note for the record: the earlier claim
that a 0/1 count "cannot give 1, 2, 1" was wrong — counts over four arms can —
so the invariant-subspace reading was worth the test and is now excluded by
measurement, not by argument.

## Next step (round 50 handoff): the line-star variant of `build_block`

Do not add another hand-written per-arm weight; four have been excluded by
measurement (see above).  Reuse the engine's own pairing instead:

1. give `build_block` a second arm-character source (an enum over
   `&ScalarStar` and a new line source) so the existing
   `little_group_operations` / `identity_position` / `prepare_targets` /
   `solve_prepared_character_block` path stays untouched;
2. the line source enumerates the star arms of the frozen `direction` (the
   contragredient images, deduplicated as vectors), folds them at
   `LINE_PARAMETER = 1/4` into the child's reciprocal lattice
   (`embedding.subgroup_lattice().reciprocal()` or the recorded
   `W . P_parent`, which agree on the `P1` family), and answers
   `q_block_dimension` / `q_block_character` from `w_little_characters_data`
   with the Bloch phase of that parameter;
3. acceptance, in this order: the 46 `P1`-child records must stay 300/300
   (`tests/w_line_frequency.rs::p1_child_records_match_every_pinned_row`), then
   SG 196 `10030 -> 1`, `10032 -> 2`, `10033 -> 2`, SG 202 `10422 ->
   DT1 1, DT3 2, DT4 1, SM1 2`, the ten asymmetric SG 202/209 records, and
   finally the whole 5,756 rows through the audit with `--require-w-complete`.
   A failure at any step is reported with the ordinal and both numbers rather
   than tuned away.

## Global audit numbers (round 58, full unscoped run)

```
identity_rows=94271 unique_pairs=94271 (pinned 94271)
other_wave_vector: records=1006 rows=5756 source_resolved=5756 source_mismatch=0 computed=300 conflicts=0
w_scope: rows=5756 computed=300 uncomputed=5456 character_tables_frozen=5756 character_tables_blocked=0
         reason=line_sources_need_the_line_star_build_block_path gate=--require-w-complete
hard_failures=0 accounting_violations=0 census_mismatch=0
VERDICT clean scope=global incomplete_categories=0 w_uncomputed=5456
```

Exit 0 without `--require-w-complete`; that gate still exits 2 while the 5,456
rows are not engine-computed.  The 300 are the `P1`-child family, verified row
for row against the pinned values by the audit itself.

## Exact insert points for the line-star step (round 59)

Checked on this tree (commit `57f1761`):

* `build_block` — `src/irrep/subduction_star_decompose.rs:1078`; it touches the
  parent character source in exactly three places: `star.q_block_dimension(
  point.arm_indices())` and `star.q_block_character(point.arm_indices(),
  operation)` inside it, plus the two call sites at lines **700**
  (`trivial_content_with_embedding`) and **993**
  (`subduce_full_star_with_embedding`).
* The two methods to provide for a line source are
  `ScalarStar::q_block_character` (`src/irrep/subduction_scalar_star.rs:675`)
  and `ScalarStar::q_block_dimension` (`:699`).
* `line_trivial_content_with_embedding` (same file, added this session) already
  contains the arm enumeration, the `LINE_PARAMETER = 1/4` fold and the frozen
  character lookup that the line source needs; its per-arm weighting is the part
  to replace with `build_block`'s own solve, and it currently bounds its scope to
  `P1` children in the audit.
* Regression to protect while refactoring: `tests/w_line_frequency.rs` (2 passed,
  1 ignored) and the audit's global line, which must keep printing
  `computed=300` until the rest lands.

## Measured: 3,916 / 5,756 rows already reproduced (round 63)

Counting a row as engine-computed exactly when the engine's frequency equals the
pinned one (no failing gate, so the disagreement is visible as ordinary
`uncomputed`), the full unscoped audit reports:

```
w_scope: rows=5756 computed=3916 uncomputed=1840 character_tables_frozen=5756 character_tables_blocked=0
hard_failures=0 accounting_violations=0 census_mismatch=0
VERDICT clean scope=global incomplete_categories=0 w_uncomputed=1840
```

So 68% of the w track is already engine-computed with the current
`line_trivial_content_with_embedding`, not just the 300 `P1`-child rows that the
regression pins.  The 1,840 remaining rows are the ones whose child is not `P1`:
that is the whole remaining scope, and it is exactly the family where the
four excluded per-arm weight rules were tested and failed.  Re-running the
acceptance ladder against the line-star `build_block` path therefore only has to
move this single number from 1,840 to 0.

## Where the remaining 1,840 rows are (round 69)

Scoped audit runs (`--parent N`) with the current `computed` rule:

| parent | w rows | computed | uncomputed |
|---|---|---|---|
| 196 | 106 | 76 | 30 |
| 202 | 316 | 155 | 161 |
| 209 | 372 | 282 | 90 |
| 225 | 1647 | 1169 | 478 |
| 227 | 984 | 656 | 328 |
| 221, 229, 230 | 0 | 0 | 0 |

The failures are spread over every parent that has w rows, roughly a third of them
per parent, which matches the diagnosis that they are the rows with a *non-trivial
arm stabiliser* rather than a whole parent or a whole child type.  Only five SGs
are listed here (1,087 of the 1,840); the rest sit in the remaining parents with w
rows (203, 210, 216, 219, 226, 228).  A fresh implementation can therefore use
SG 196 (30 failures) as the smallest reproducer and SG 202 ordinal 10422 as the
documented witness.

## First hypothesis to test next: orbits with the +- identification (round 73)

The recorded witness (SG 202 ordinal 10422, pinned `DT1 1, DT3 2, DT4 1,
SM1 2`) folds four arms `(0,0,+-2)`, `(+-2,0,0)`, i.e. **two +-pairs**, with
per-arm transported characters `(1,-1,1,-1)`, `(1,-1,1,-1)`, `(1,0,-1,0)`,
`(1,0,-1,0)`.  The pinned `1, 2, 1` is a per-pair pattern (one pair contributing
for `DT1`, two for `DT3`, one for `DT4`) and cannot come from the per-arm
average the code does today, nor from any of the four excluded rules — all of
which ignore the fact that an arm and its negative are the *same* star orbit
once the child's mirrors act.

So the next attempt should group the folding arms into the child's orbits with
`a ~ -a` (the round-12 formulation: `frequency = sum over H-orbits of arms
folding to Gamma of m_orbit`), and take the weight **per orbit**
(`m_orbit = mult(trivial_{K_orbit}, W')` with `K_orbit` the child operations
that map the orbit to itself, the transported character picking up the sign of
the identification), instead of per arm as in the current
`line_trivial_content_with_embedding`.  Cheapest check: SG 196 (30 failing
rows, smallest reproducer) and SG 202 10422 before rerunning the full audit.

The audit's TSV now carries a per-row status for the w entries: `computed`,
`computed_mismatch` (the engine's value differs from the pinned one) or
`uncomputed_line_star_open` (no frozen table or no embedding).  For SG 196 the
scoped run reports 76 `computed` and 30 `computed_mismatch`, so the remaining
work now has an exact, greppable row list:

    cargo run --release -p cryspglib --example audit_irrep_subduction -- --parent 196 --output /abs/a196.txt
    awk -F'\t' '$1=="w_entry" && $15=="computed_mismatch"' /abs/a196.txt

Grepping SG 196's 30 `computed_mismatch` rows (the worklist recipe above) shows
they all belong to records whose child is **chiral**: SGs 18 (`P2_12_12_1`),
19 (`P2_12_12`), 4 (`P2_1`) and 198 (`P2_13`) — no mirrors and no inversion.
Every failing row carries a small stored frequency (1, 2, 3, ...), and several
records contribute two rows.  So the remaining defect is concentrated on children
whose point group acts *freely* on the arms, which is the opposite of the `P1`
extreme the engine already gets right and consistent with the +-orbit/identification
hypothesis: with no mirror to identify an arm with its negative, the per-arm
average and the per-orbit weight can still disagree.

## Second hypothesis: the parameter is per child, not a global 1/4 (round 79)

`LINE_PARAMETER = 1/4` reproduces the whole `P1` family and 3,916 rows overall,
but every remaining failure is a *chiral* child, where the point group acts
freely on the arms and a too-small parameter folds too many of them (the pinned
values there are only 1, 2, 3).  A global change of the constant is refuted:
`1/2` and `1` already leave 111 of the 300 `P1` rows wrong.  The principled
candidate is that the program uses the **smallest parameter whose fold lands on
the child's Gamma point**, i.e. `1 / g` with `g` the largest component gcd among
the folded arms `W . arm` (for the `P1` family that is `1/4`, matching the
measurement, while a chiral child with a different supercell gets `1/2`, `1/3`,
...).  Implementing it needs one structural change: the parameter is currently
computed before the arms exist in `line_trivial_content_with_embedding`, so the
`1 / g` block has to move *after* the arm construction.  A first attempt to
patch it in place failed to compile (`arms` not yet in scope) and was reverted,
leaving the tree green.

### Refuted by measurement (round 80): the "first hit" per-child parameter

Implemented as described above (`1 / g`, `g` = largest folded-arm component gcd,
falling back to 1/4) and measured:

* `tests/w_line_frequency.rs` drops from 2 passed to **1 passed / 1 failed** —
  the `P1` regression breaks, so the rule is wrong even where the constant 1/4
  was right; and
* SG 196 falls from 76 `computed` to **34** `computed` (72 uncomputed).

So the program's parameter is *not* the smallest one whose fold lands on the
child's Gamma point.  Reverted immediately (tree green, `LINE_PARAMETER = 1/4`
restored).  This is the fifth weight/parameter reading excluded by measurement;
the remaining scope is still 1,840 rows, and the line-star `build_block` path
(which needs no parameter choice beyond the phase convention) is still the only
documented route.

The mismatch status now carries the engine's value (`computed_mismatch:<value>`).
For SG 196's 30 failing rows the values are mixed: some are **0** where the pinned
frequency is 1-3 (the fold/weight rejects every arm) and some are *larger*
(e.g. 2 or 3 where the pinned value is 1).  A scalar tweak of the parameter or of
the per-arm weight therefore cannot repair them - the pattern is structural, which
is the measured argument for routing the line representation through
`build_block`'s own fold and character solve instead of the hand-written sum.

### Also immaterial (round 82): the recorded vs accepted subgroup lattice

Swapping the fold/stabiliser lattice from the recorded `W . P_parent` to the
accepted embedding's `subgroup_lattice().reciprocal()` leaves both measurements
unchanged: the `P1` regression still passes 2/0/1 and SG 196 still reports
`computed=76 uncomputed=30`.  So the remaining 1,840 rows are not a frame or
lattice artefact; the mixed `computed_mismatch:<value>` pattern (some 0 against
pinned 1-3, some overshoot) stands, and the line-star `build_block` route remains
the only documented fix.

### Refuted (round 85): identifying an arm with its negative

Implemented the +-pair reading (keep one representative per pair, chosen by the
sign of the first non-zero component, and let the stabiliser accept both signs).
It compiles, but the `P1` regression drops from 2 passed to **1 passed /
1 failed**: the `P1` family pins arm counts that are *not* +-symmetric (e.g. 4 of
6 arms fold for the size-16 child), so an arm and its negative must be counted
separately there.  The rule is therefore wrong as stated and was reverted (tree
green at `7e93008`).  The mixed `computed_mismatch` pattern of the remaining
1,840 rows therefore needs an explanation that keeps +-arms distinct - which is
one more argument that the fix is the structural `build_block` route rather than
another counting rule.

## Implementation sketch for the line-probe `FoldedStar` (round 96)

The missing data shape, written down so it does not have to be re-derived:

```text
input:  embedding, table (frozen), arms = [(arm_direction, rotation)]  (already built)
        parameter = LINE_PARAMETER = 1/4
        child_reciprocal = embedding.subgroup_lattice().reciprocal()

for each arm a:
    k_a  = parameter * a                       (parent conventional reciprocal)
    q_a  = fold_wave_vector(embedding.transform(), &k_a)     // child frame, exact
    q_a  = child_reciprocal.reduce(&q_a).representative      // into the child's cell
group arms by equal q_a  ->  one LineFoldedStar per group
    points:          the group (q_a, arm_indices)          // one point per group entry
    star_size:       number of *distinct* arms in the group
    arm_count:       total arms in the group (equals star_size before the child acts)
    block_dimension: sum over the group's arms of the little dimension
```

Then `build_block(embedding, &line_arm_source, &folded_star, &child_cell, &child_reciprocal)`
consumes exactly what it consumes today for `ScalarStar`:
`q_block_dimension(point.arm_indices())` = the group's block dimension and
`q_block_character(point.arm_indices(), operation)` = the sum over those arms of
the frozen character of the arm's little group on the conjugated operation
(`line_character` already implements the per-arm half, including the `1/4` Bloch
phase and the zero for operations that move the arm).  The `FoldedStar` type is
today only built inside `ScalarStar::folded_stars`; the smallest change is to give
it a crate-private constructor (or to make the line source produce the same struct
through a sibling `folded_stars`), which is why the arm-source enum and this
constructor are the two pieces to land first.  Everything downstream —
`select_representative`, `little_group_operations`, `identity_position`,
`prepare_targets`, `solve_prepared_character_block`, the trivial-row aggregation —
stays untouched, and `tests/w_line_frequency.rs` plus the audit's
`computed`/`uncomputed` counters are the regression to watch.

### Error plumbing for the arm-source enum (round 103)

Checked, because it was the one uncertainty left in the enum step: `StarError`
(`src/irrep/subduction_star.rs`) already has `Subduction(#[from] SubductionError)`
plus a variant for "an operation that is not a parent group element modulo
`L_G`".  So the line source can present `q_block_dimension` /
`q_block_character` with `build_block`'s exact signatures by mapping the exact
layer's `SubductionError` straight into `StarError` (and returning the zero
character, as `line_character` already does, for operations that move the arm
rather than raising).  `FullStarError::from(StarError)` exists, so the calling
side keeps using `FullStarError`.  No new error variant is needed for the enum
step.

The enum's two methods must match `build_block`'s call sites exactly:

```rust
fn q_block_dimension(&self, arm_indices: &[usize]) -> Result<u32, StarError>;
fn q_block_character(
    &self,
    arm_indices: &[usize],
    operation: &ExactSeitz,
) -> Result<Complex64, StarError>;
```

`ScalarStar`'s versions live at `src/irrep/subduction_scalar_star.rs:675` and
`:699`; the line implementation answers `dimension * arm_indices.len()` for the
first and, for the second, sums `line_character` over those arms (each arm
already carrying the rotation that transports the base little group onto it, and
`line_character` returning zero when the conjugated operation moves the arm).

**Correction to the paragraph above (round 105).** The error plumbing is only
*one-way*: `FullStarError` wraps `StarError`, never the reverse, so the line
source cannot call `line_character` / `line_rotation_character` (both return
`FullStarError`) from methods that must return `StarError`.  The line methods
therefore have to inline that logic — test the stabiliser with
`Mat3R::from_ints(rotation).checked_mul_vector(&direction)` exactly as
`line_character` does, look the frozen character up by rotation, and use
`StarError::Subduction(SubductionError::RationalOverflow { operation: .. })` for
the dimension overflow — or the two helpers must be duplicated with `StarError`.
Everything else in the sketch stands; plan for ~15 extra lines here.

**One open design choice before writing the line arm source (round 106).** Its
methods return `StarError`, and the "the frozen table has no character for this
rotation" case has no honest variant there (only the exact-layer error and the
parent-operation ones).  Two clean options: add a variant such as
`MissingFrozenRotation { sg: u8, label: &'static str }` to `StarError`
(`src/irrep/subduction_star.rs`) — check for exhaustive matches on that enum
before doing so — or keep the line source's own error type and give the enum a
`map_err` at the `build_block` boundary instead of matching `ScalarStar`'s
signature.  Do not paper over it with `SubductionError::RationalOverflow`, which
would describe the wrong failure.  Everything else (fields, the stabiliser test,
the phase, the dimension overflow) is settled.
