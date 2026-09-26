# Handoff: parameterized k-line subduction

Date: 2026-09-26
Repository: `cryspglib`
Code checkpoint: `6214a61` (`test: extend line coverage through denominator twelve`)
Previous implementation checkpoint: `64bb535` (`feat: add D4 projective line targets`)

This handoff describes the completed R6.4 work and the remaining route toward broader rational-parameter coverage. It does **not** claim that arbitrary rational parameters, magnetic corepresentations, or spinors are complete.

## Project objective

The library's larger goal is to support Landau and phase-transition workflows internally: select an isotropy subgroup from a parent space group, wave vector, irrep and order-parameter direction, then restrict other parent irreps to that subgroup and return all child irreps with multiplicities.

The active R6 branch of work concerns ordinary scalar irreps on parameterized k-line families. The immediate engineering goal is to close cases where the exact folding and little group are known but the projective child representation catalogue cannot yet supply a complete target table. Any unsupported case must continue to fail explicitly; never return a partial decomposition as complete.

## What is complete at this checkpoint

### D4 projective targets (`64bb535`)

- The projective catalogue supports order-eight D4/cubic `C4v` little co-groups when the cocycle is verified to be a coboundary.
- It constructs four gauge-twisted one-dimensional targets and the gauge-twisted ordinary two-dimensional standard target.
- Recognition checks the order profile, the D4 relations, generation, and the unique central involution. The coboundary check is fail-closed.
- The general one-dimensional solver enumerates finite roots from generator orders rather than searching a grid whose cost grows with cocycle denominator. `MAX_ORDER` is 8; higher-order groups and non-coboundary D4 cases still reject explicitly.

### Rational samples through denominator 12 (`6214a61`)

For every reduced `t = n/d` with `0 < t < 1`, `3 ≤ d ≤ 12`, and `t ∉ (1/4)Z`, the explicit CLI gate decomposes all 5,756 parameterized line rows. This is 42 exact parameters. On the final run every point reported:

```text
full=5756/5756, missing=0, other_errors=0, content_mismatches=0
```

Three D4 witnesses (ordinals 13543, 14106 and 13691) have a permanent library regression over all 42 points. It requires a constructed two-dimensional child target, zero trivial content, exact dimension coverage and per-operation reconstruction. The example test independently proves that its static CLI parameter list is exactly the reduced off-quarter set and checks the displayed fraction labels. The library regression generates the same domain from the mathematical rule, avoiding a second hard-coded list.

The full sweep is an explicit slow gate rather than part of the default example tests. Run it when changing parameterized subduction, the projective catalogue or folding logic, and before updating the 42-point coverage claim:

```bash
cargo run --release -p cryspglib --example line_family_coverage -- \
  --projective-sample-sweep --gate
```

The result is finite-sample evidence only. It does not establish support for every rational `t`.

## Verification completed for `6214a61`

- `cargo test --release -p cryspglib --examples`: **40 passed** across the example test binaries.
- `cargo test --release -p cryspglib --lib sampled_d4_gap_witnesses_decompose_with_the_constructed_two_dimensional_target`: **1 passed**.
- Exact example sample-list regression: **1 passed** after the final set/label assertion change.
- `cargo clippy -p cryspglib --all-targets --release -- -D warnings`: clean.
- The 42-point full-table CLI sweep with `--gate`: exit 0; all 42 rows of parameter checks complete with zero missing cases, other errors, or content mismatches; the existing line-family gate also passed.
- DSH read-only review: **no blockers**. It independently checked the finite fraction set, the three witness path, disjoint failure counters, fail-closed gate and documentation scope. Non-blocking wording nits found in that review were corrected before the code commit.

The full crate suite and the 366,260-probe general subduction audit were not rerun for this test/documentation-only extension; their last green results are recorded in `AGENTS.md`. Do not describe them as newly run on `6214a61` unless rerun.

At the previous implementation checkpoint `64bb535`, the global ordinary-scalar isotropy audit had passed **366,260/366,260 full decompositions**, with `identity_only=0` and no hard failures. The official `t=1/4` line-frequency check also covered all **5,756/5,756** rows. Those are prior-checkpoint results, not runs from this handoff turn.

## Current limits that remain

- Arbitrary rational `t` is not proved or fully covered. The 42-point scan covers only denominators 3–12 outside the quarter grid.
- The supported projective representation families remain bounded: generic one-dimensional targets up to group order 8, and the already verified higher-dimensional `C2×C2`, coboundary `D3`, and coboundary `D4` cases. Non-coboundary D4 and other/higher-order little co-groups still return a typed missing-target error.
- General-parameter nontrivial child irreps have no official complete-decomposition oracle in the shipped data. Current checks are internal exact/invariant checks (dimension coverage and operation-by-operation reconstruction); the official pinned `t=1/4` rows provide the external frequency anchor.
- The current parameterized source corpus is the 73 frozen ordinary scalar line sources (all corresponding to the shipped w table). This does not cover every line in all 230 space groups.
- Magnetic corepresentations, the 3,611 spinor records, and user-facing API promotion are outside the active R6 scope; see `AGENTS.md` §2.
- Keep the distinction between actual line irreps and formal induced values. `LineSubduction::parameter_kind()` labels `Formal` points such as `t=0` and some `t=1/2` cases; a successful formal decomposition is not evidence for an enhanced little-group line irrep.

## Suggested next work

Do not start by adding more arbitrary sample points. First replace sampling guesses with an exact parameter-domain census.

1. **Derive exact little-group transition parameters.** For each frozen source and each operation, solve the exact reciprocal-lattice condition for the operation to fix `k(t)=t·v`. Use the current conventional reciprocal coordinates and exact `Rat`/lattice helpers. Enumerate all exceptional rational `t` in a fundamental interval; separate the generic stabilizer from isolated enhanced-symmetry points.
2. **Classify the exact projective data on each domain.** For every resulting little co-group, compute its multiplication table and factor system as a function of `t`. Determine where the cocycle is gauge-trivial and where it represents a nontrivial projective class. Preserve nonsymmorphic Bloch phases; do not identify a continuous phase with its value at a nearby sample.
3. **Use the census to choose the next representation family.** Count affected `(source, subgroup, parameter-domain, co-group, cocycle-class)` cases. The next smallest family may be a higher-order dihedral group, a non-coboundary D4 class, or a parameterized target whose characters depend on `t`; decide from the exact census rather than from denominator samples.
4. **Implement one bounded family at a time.** Reuse the existing projective character solver and catalogue structures where possible. Keep construction exact, check character orthogonality and squared-dimension sum, then require the existing full-star dimension and per-operation reconstruction gates. Add a negative case proving unsupported classes still fail closed.
5. **Add one representative exact end-to-end gate and recount.** Report which parameter domains are covered and which remain unsupported. Only claim arbitrary-rational coverage if every domain from step 1 has a complete target source and reconstruction evidence.
   **Superseded (2026-09-26, external review):** "one representative parameter per domain" is withdrawn — the partition must cover the whole star (every folded arm's own little co-group and every pairwise arm identification), not the reference direction alone, and two parameters *inside* one interval must be compared rather than one representative trusted. See `AGENTS.md` section 3b, cards 3-6, and `docs/subduction-r6-coverage.md` section 4b; witnesses ordinal 10030 and 10038 both change at `t = 1/8`.

Useful source/document entry points:

- `src/irrep/subduction_star_decompose.rs`: exact full-star folding and decomposition.
- `src/irrep/subduction_catalogue.rs`: constructed projective little-group targets and finite character solver.
- `src/irrep/line_monodromy.rs`: exact transport under reciprocal shifts.
- `examples/line_family_coverage.rs`: 5,756-row summary, existing gate and 42-point extended gate.
- `docs/subduction-r6-coverage.md`: evidence levels, current bounds and performance caveats.
- `docs/subduction-r6-plan.md`: R6 task history and coordinate/frame conventions.
- `AGENTS.md`: active repository workflow, validation gates and known limits. Read its opening section first.

## Handoff workflow notes

- Do not run broad `cargo fmt` over the existing files; the repository contains pre-existing formatting differences. Use the current Clippy gate and keep diffs focused.
- Keep the 42-point full scan out of the default fast test path; invoke the explicit command above when changing relevant math or coverage claims.
- If using DeepSeek, follow the `dsh-worker` skill: headless DSH only, a scoped assignment, explicit file ownership for implementation, review its diff in Codex, and wait for a read-only audit to finish. Do not use the retired `v4_flash_worker` role.
- Local commits are expected at bounded task boundaries. Do not push unless the user asks.
