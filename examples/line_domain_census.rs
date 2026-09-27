//! Exact parameter-domain census of the frozen parametric-k line sources.
//!
//! R6.5 replaces denominator sampling with an exact partition of the parameter
//! line.  For every frozen source the module
//! [`cryspglib::irrep::subduction::star::line_domain`] computes, in exact
//! rational arithmetic, the parameters `t in [0, 1)` at which some parent
//! operation joins the little group of `k(t) = t . v`; the same computation is
//! repeated after folding `v` through each isotropy record's embedding, which is
//! where the *child* little co-group of the folded star grows.
//!
//! This example joins the two halves with the engine itself: at every parameter
//! of the partition (plus one generic sample per source) it runs the production
//! decomposition and records
//!
//! * the class of the answer -- `stored` (pinned child rows answered it),
//!   `constructed` (the projective/bloch catalogue built the targets),
//!   `unsupported` (`MissingChildStarData`), or `error`;
//! * the child little co-group order at that parameter from the census;
//! * the production `ParameterKind`, which must equal the census's notion of an
//!   enhanced parameter;
//! * the trivial content, so the anchor rows stay comparable with the pinned
//!   table.
//!
//! R6.7 card 4 binds every statistic to the individual **child-star block**
//! instead of to the whole probe.  A probe-level `stored`/`constructed` label is
//! a statement about *some* block: it says nothing about the block whose little
//! co-group or cocycle one is asking about, and the accepted external finding is
//! exactly that -- ordinal 13688 (SG 225 `SM1` -> child #134, `t = 1/8`) carries
//! an order-16 non-trivial block that is entirely **stored** while a *different*
//! block of the same probe holds constructed targets.  The block-level reading is
//! therefore the primary one here:
//!
//! * every `result.blocks()` entry yields a [`BlockStat`] (record, source,
//!   parameter, block index) with its own star size, arm set, dimension, little
//!   co-group, own cocycle class and own target source;
//! * the `--gate` boundary table pairs a non-trivial class at the reference
//!   folded point with the **reference-carrying block**'s own source, and keeps
//!   the old probe-level table only as a clearly labelled withdrawn convention;
//! * a full-star recount pass runs the production decomposition at every
//!   parameter the card-3 partition [`full_star_partition`] adds beyond the
//!   reference candidate set and classifies every block there, and now emits one
//!   `--output-recount` fingerprint row per probed parameter;
//! * `--domain-sweep` (R6.7 card 6) first runs the **production** decomposition
//!   at every parameter of `probe_parameters(parent)` -- the boundary points,
//!   which the interval comparison never touches -- and counts them
//!   (`boundary_points == boundary_successes + boundary_failures`, with every
//!   failure carrying the record, the label, the parameter and the engine's own
//!   error), and it recomputes the arm geometry **itself** from the folded arms
//!   and the child reciprocal lattice: every relation
//!   `t (R^-T v_i - v_j) in L*_H` is enumerated over the same triples the
//!   partition enumerates and decided with the census's own arithmetic (the
//!   difference vector, mapped into the child lattice basis, scaled by `t` and
//!   tested for integrality), the boundary set is re-derived from the relations'
//!   residue sets, the event counts at every boundary are compared with the
//!   partition's, no non-permanent relation may hold at an interior point of any
//!   interval, and the arm count is recomputed as an orbit-stabiliser count over
//!   the parent's rotations.  The partition's solver is never called on this
//!   path, so a partition that reports a wrong parameter set, a wrong count or a
//!   boundary no relation supports cannot agree with the audit by construction.
//! * It then compares two exact interior points of every interval of
//!   `probe_parameters(parent)` in two layers.  **Geometry**: blocks
//!   are matched one-to-one by their parent arm set -- the folded coordinates,
//!   the block order, the engine's representative and a constructed target's
//!   enumerated index all move with `t`, while which arms share a child star does
//!   not -- and the star size, arm counts, per-point arm grouping, dimensions and
//!   little co-group rotation sets must agree.  **Representation**: for each
//!   matched block the targets are compared by their little-group characters at a
//!   common seed arm, with the card-1 gauge `exp(2 pi i q(t) . tau)` divided out
//!   ([`FullStarBlock::target_little_character`]), matched one-to-one by
//!   character inner product and then compared multiplicity by multiplicity.
//!
//! The sweep compares **characters**, never a sorted `(dimension, multiplicity)`
//! table: two distinct one-dimensional targets whose multiplicities are swapped
//! leave that table unchanged.  It also never compares the induced full-star
//! character directly -- that trace sums over the star's arms, each with its own
//! parameter-dependent Bloch phase, so it is not related to its value at another
//! parameter by any single gauge (measured: the matched blocks of the corpus
//! score 0.73-1.0 that way, while the per-arm reading scores 1 - 2.2e-16).
//!
//! Usage:
//!
//! ```text
//! line_domain_census                 # summary + parameter partition
//! line_domain_census --gate          # exit 1 unless every invariant holds
//! line_domain_census --full-star-recount
//!                                    # rerun the gate's card-3 recount pass alone
//! line_domain_census --domain-sweep  # the card-5 interval sweep (own failures
//!                                    # are gate violations and exit 1)
//! line_domain_census --output out.tsv
//! line_domain_census --output-blocks blocks.tsv
//! line_domain_census --output-recount recount.tsv
//! line_domain_census --sequential    # one thread (default: rayon over records)
//! ```
//!
//! The gate deliberately accepts `unsupported > 0`: the census exists to *name*
//! the unsupported domains, not to pretend they are closed.  `--require-covered`
//! turns a non-empty unsupported set into a failure, for the round that closes
//! them all.
//!
//! Measured scope of the sweep on this corpus (8 threads, 1 m 39 s - 1 m 52 s
//! over repeated runs of `--gate --require-covered --domain-sweep`): 5,756
//! (record, label) pairs,
//! 46,048 intervals, 92,096 interior decompositions, 46,048 comparisons, 0
//! failures; the card-6 boundary pass answers all **46,048** probe parameters
//! (0 failures, 429,888 parent dimension, equal to the covered dimension and to
//! the closed form `parameters x little dimension x arms`); the independent
//! audit enumerates **2,477,298** relations (the closed form `n^2 r - n`),
//! **205,044** of them identically true (the partition's own `permanent_counts`
//! summed), re-derives the same **46,048** boundary parameters, evaluates
//! **18,178,032** predicate instances at the boundaries and **36,356,064** at
//! the interior points (twice as many, and exactly
//! `evaluated x pairs == boundaries x (relations - permanent)`), finds no
//! relation at any interior point and **0** disagreements, and reproduces the
//! **50,226** arms as an orbit-stabiliser count (5,756 checks, 0 disagreements);
//! 168,408 matched block pairs / 174,672 matched target pairs; every
//! interior target is **constructed** (0 stored, 0 mixed), which is why the
//! stored/constructed branch is pinned by a module test; 6,264 interior blocks
//! carry two targets and all of them have the shape `[(1, 1), (1, 1)]`, so a
//! multiplicity swap between them is a no-op on this corpus (measured, and
//! reported as `same_dimension_swaps = 0`).

use cryspglib::irrep::line_monodromy::{line_direction, line_table};
use cryspglib::irrep::subduction::star::decompose::{
    FullStarBlock, FullStarError, FullStarTarget, LineSubduction, ParameterKind,
    official_line_parameter, subduce_line_at_parameter,
};
use cryspglib::irrep::subduction::star::line_domain::{
    FoldedArm, FullStarPartition, GammaParameters, ParentDomain, StarBoundary, StarEventKind,
    child_cocycle_is_a_coboundary, child_exceptional_parameters, full_star_partition,
    little_co_group_order, minimal_parameter_step, minimal_parameter_step_via_coordinates,
    parent_domain, point_cocycle_is_a_coboundary, point_little_co_group_order,
    reciprocal_lattice, require_reciprocal_direction, rotation_set, verify_against_grid,
};
use cryspglib::irrep::generated_data::SG_DATA_HALL;
use cryspglib::irrep::subduction::{
    ExactSeitz, Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R, bloch_phase, fold_wave_vector,
};
use num_complex::Complex64;
use cryspglib::irrep::w_little_characters_data::{LittleCharacterTable, W_LITTLE_CHARACTERS};
use cryspglib::irrep::{LabelConvention, isotropy, query};
use cryspglib::mathfunc::Mat3I;
use cryspglib::{HallNumber, SymmetryOps};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{BufWriter, Write as _};
use std::process::ExitCode;
use std::sync::Mutex;

/// The `--output-blocks` header, shared by the writer and the read-back so the
/// column names are validated too (verification review F3: swapping two header
/// names left the gate green).
const BLOCK_TSV_HEADER: &str = "ordinal\tparent_sg\tchild_sg\tlabel\tparameter\tindex\t\
star_size\tarm_count\tarm_indices\tblock_dimension\tlittle_co_group\t\
cocycle_trivial\tblock_source\tterms\tcarries_reference\tcarries_gamma";

const USAGE: &str = "\
line_domain_census [--gate] [--require-covered] [--full-star-recount] [--domain-sweep]
                   [--sequential] [--output <path>] [--output-blocks <path>]
                   [--output-recount <path>]

  --gate               exit 1 unless the census invariants hold (includes the
                       full-star recount pass)
  --require-covered    additionally fail when any (record, parameter) is unsupported
  --full-star-recount  run the card-3 full-star recount pass on its own, with
                       the gate's failure semantics
  --domain-sweep       R6.7 cards 5-6: run the production decomposition at
                       every parameter of `probe_parameters(parent)` (the
                       boundary points, counted and conserved), recompute the
                       full-star arm geometry independently of the partition's
                       solver (relations, boundary set, event counts, interior
                       relation-freeness, orbit-stabiliser arm count), and
                       compare two exact interior points of every interval by
                       stable block identity (parent arm sets) and, per matched
                       block, by gauge-unified little-group characters; failures
                       are gate violations
  --sequential         one thread (default parallel over isotropy records)
  --output <path>      write the per-probe table as TSV
  --output-blocks <path>
                       write the per-block table as TSV
  --output-recount <path>
                       write the per-recount-probe geometry fingerprint as TSV
                       (implies the recount pass; its row count is asserted)
";

/// The `--output-recount` header, written and re-read independently so a changed
/// column cannot move both sides at once (the card-4 audit lesson for
/// `--output-blocks`).
const RECOUNT_TSV_HEADER: &str = "ordinal\tparent_sg\tchild_sg\tlabel\tparameter\t\
block_count\tblocks";

/// The exact fraction of an interval's span used for its two interior points:
/// the trisection points `left + span/3` and `left + 2 span/3`, both strictly
/// inside because `span > 0`.
const INTERIOR_FRACTION: (i128, i128) = (1, 3);

/// Tolerance of the gauge-unified character matching.  The corpus agreement is
/// at machine precision (measured worst score `1 - 2.3e-16` over the whole
/// corpus), so this only absorbs `f64` trigonometry, not a modelling gap.
const SWEEP_TOLERANCE: f64 = 1e-9;

/// Matched block pairs whose **reduced** representative `stored_k` differs
/// between the two interior points of one interval.
///
/// Measured: **every** pair (168,408 of 168,408) — the folded point moves
/// continuously, so its reduced representative moves with it, and `stored_k` is
/// therefore not a cross-parameter identity (which is why the comparison is on
/// gauge-unified content).  The card-5 math review found the stronger, separate
/// phenomenon of a representative that jumps **discontinuously** inside an
/// interval (a reduction wall: 96 of them on a denser grid, witness ordinal
/// 13719 SG 225 `SM1` at `t = 1/6`, where the gauge-unified content still agrees
/// to 3.3e-16); this counter does not measure those walls, it pins the volume
/// fact that `stored_k` is not usable as the cross-parameter key.
const REPRESENTATIVE_SHIFT_PIN: usize = 168_408;

/// The identity rotation: every aligned little co-group contains exactly one, so
/// the non-identity part of the fixity control's checks is predicted by the
/// little-co-group histogram (the operand binding of verification round 2).
const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

// ── R6.7 card 6 pinned totals ────────────────────────────────────────────────
//
// Every number below is printed by a plain run before the gate asserts it, so a
// reviewer can read the measured value and the pinned value side by side.  Each
// has an independent source next to its assertion in `check_invariants`.

/// `Σ partition.probe_parameters(parent).len()`: the boundary pass runs the
/// **production** decomposition once per probe parameter per pair, so this is
/// also `boundary_successes + boundary_failures` (asserted).  The sweep's probe
/// set is `probe_parameters`, which on this corpus is the eighth grid of every
/// pair (`1/7` is a main-probe sample, not a boundary): 8 x 5,756.
const BOUNDARY_POINT_PIN: usize = 46_048;

/// Relations the audit enumerates over the corpus: `Σ (n^2 r - n)` with `n` the
/// partition's arm count and `r` its child rotation count.  The gate compares
/// this against the closed form accumulated independently per pair.
const RELATION_PIN: usize = 2_477_298;

/// Relations whose difference vector is exactly zero (true at every parameter),
/// and the same number as the partition itself accounts for (`permanent_counts`
/// summed).
const PERMANENT_RELATION_PIN: usize = 205_044;

/// Predicate instances `t w in L*_H` the audit evaluates at the partition's
/// boundary parameters: every boundary evaluates the pair's whole non-permanent
/// relation list, so this is `8 x (relations - permanent)` on this corpus (the
/// gate also asserts that identity and the pinned value).
const BOUNDARY_RELATION_EVALUATION_PIN: usize = 18_178_032;

/// Predicate instances the audit evaluates at the two interior points of every
/// interval: `16 x (relations - permanent)` on this corpus, twice the boundary
/// total because each pair has two interior points per interval against one
/// boundary parameter each.
const INTERIOR_RELATION_EVALUATION_PIN: usize = 36_356_064;

/// `Σ result.parent_dimension()` over the answered boundary parameters: the
/// engine's own dimension identity summed over the card-6 boundary pass.  The
/// gate pins it (it is also asserted equal to the covered-dimension sum), so a
/// pass that does not really run the engine cannot keep its counters green.
const BOUNDARY_DIMENSION_PIN: u64 = 429_888;

/// Synthetic partitions the audit's gate-level self-check drives: the honest
/// fixture, the lowered-count fixture and the moved-boundary fixture.
const SELF_CHECK_CASE_PIN: usize = 3;

/// `Σ partition.arms.len()`: the arm count the orbit–stabiliser route predicts
/// from the parent's rotation group alone (73 sources, 5,756 pairs).
const ARM_ORBIT_PIN: usize = 50_226;

/// The generic parameter used to prove that a non-exceptional point is answered
/// from the frozen line little group.  It is off the quarter grid and not one of
/// the exceptional parameters of any frozen source (`{0, 1/2}`).
const GENERIC_SAMPLE: (i128, i128) = (1, 7);

/// How the engine answered one (record, label, parameter) probe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum TargetClass {
    /// Pinned child rows carried the folded star.
    Stored,
    /// The constructed catalogue (Bloch phase or projective targets) answered.
    Constructed,
    /// `MissingChildStarData`: the census must explain it as an enhanced child
    /// little co-group.
    Unsupported,
    /// Any other engine error: a hard failure of the census run.
    Error,
}

impl TargetClass {
    const fn label(self) -> &'static str {
        match self {
            Self::Stored => "stored",
            Self::Constructed => "constructed",
            Self::Unsupported => "unsupported",
            Self::Error => "error",
        }
    }
}

/// How the targets of one **child-star block** were sourced.
///
/// This is the R6.7 card 4 replacement for the probe-level [`TargetClass`]: the
/// old label answered "does *some* target of this probe come from a pinned child
/// row", which is not a statement about the block whose little co-group or
/// cocycle is under discussion.  The accepted witness is ordinal 13688 (SG 225
/// `SM1` -> child #134, `t = 1/8`), where the non-trivial order-16 block is
/// entirely stored while another block of the same probe holds constructed
/// targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum BlockSource {
    /// Every target of the block came from a pinned child row.
    Stored,
    /// No target of the block came from a pinned child row.
    Constructed,
    /// The block mixes both sources.
    Mixed,
}

impl BlockSource {
    /// Position in a counts array.
    const fn index(self) -> usize {
        match self {
            Self::Stored => 0,
            Self::Constructed => 1,
            Self::Mixed => 2,
        }
    }

    /// Short human-readable name, for the report.
    const fn label(self) -> &'static str {
        match self {
            Self::Stored => "stored",
            Self::Constructed => "constructed",
            Self::Mixed => "mixed",
        }
    }
}

/// Classify one block from **its own targets and nothing else**.
///
/// The signature is the control: it receives one slice, so adding, removing or
/// reclassifying an unrelated block of the same probe cannot move this block's
/// label.  An empty slice is `Stored` ("every target is stored" holds
/// vacuously); the probe-level rule below shows why that is the harmless
/// reading rather than a claim about data that does not exist.
fn classify_block(targets: &[FullStarTarget]) -> BlockSource {
    classify_counts(
        targets.iter().filter(|target| target.irnumber.is_some()).count(),
        targets.len(),
    )
}

/// The classification rule as a function of the two counts the block itself
/// determines.
///
/// Split out so the gate can recompute a recorded label from the counts it
/// stored while reading the block, instead of calling [`classify_block`] on the
/// same slice again: with a single function the gate check would be a
/// restatement, and a mutation that classifies every block by the probe-level
/// rule would pass it.
const fn classify_counts(stored: usize, total: usize) -> BlockSource {
    if stored == total {
        BlockSource::Stored
    } else if stored == 0 {
        BlockSource::Constructed
    } else {
        BlockSource::Mixed
    }
}

/// The legacy **probe-level** rule, kept only to measure what it would have
/// said.  It is the rule the card-4 finding withdraws: `stored` means "every
/// target of the whole probe is stored", and one constructed target in any
/// block makes every block of the probe `constructed`.
///
/// It has exactly two outcomes -- a probe is never `mixed` under it -- which is
/// the shape of the objection: a probe that carries both kinds is reported as if
/// all of it were constructed, and the stored block's own answer is lost.
fn classify_probe(targets: &[&FullStarTarget]) -> BlockSource {
    if targets.iter().all(|target| target.irnumber.is_some()) {
        BlockSource::Stored
    } else {
        BlockSource::Constructed
    }
}

/// The little co-group order and cocycle class of one exact child point, memoized.
///
/// Both quantities are functions of the child space group and the **exact**
/// point, and the same point is asked for over and over: at `t = 0` every probe
/// of a child space group asks about `q = (0, 0, 0)`, where the little co-group
/// is the whole child point group and the class decision is the expensive one
/// (the order-independent syzygy solve of M2.2).  The memo is keyed by the exact
/// rational point, so a hit means an identical input -- no reduction, no coset
/// identity, is assumed here.
#[derive(Default)]
struct PointClassCache {
    entries: Mutex<PointClassTable>,
}

/// The memo's key: the child space group and the **exact** rational point.  A
/// named alias because the type is otherwise unreadable at the use sites.
type PointClassTable = HashMap<(u8, [Rat; 3]), Result<PointClass, String>>;

/// What one child point contributes to a block's statistics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PointClass {
    /// Order of the little co-group fixing the point modulo the child
    /// reciprocal lattice.
    little_co_group: usize,
    /// Whether the point's own factor system is a coboundary.
    cocycle_trivial: bool,
}

impl PointClassCache {
    /// The class data of one point of one child group.
    fn classify(&self, child_sg: u8, point: &Vec3R) -> Result<PointClass, String> {
        let key = (child_sg, *point.as_array());
        let lock = || {
            // Defensive guard, not a swallowed corpus result: a poisoned memo
            // mutex means another worker panicked, and the memo holds only
            // recomputable point classes, so the inner table is taken as is.
            self.entries
                .lock()
                .unwrap_or_else(|poison| poison.into_inner())
        };
        if let Some(hit) = lock().get(&key) {
            return hit.clone();
        }
        let computed = point_little_co_group_order(child_sg, point)
            .map_err(|error| error.to_string())
            .and_then(|little_co_group| {
                point_cocycle_is_a_coboundary(child_sg, point)
                    .map(|cocycle_trivial| PointClass {
                        little_co_group,
                        cocycle_trivial,
                    })
                    .map_err(|error| error.to_string())
            });
        lock().insert(key, computed.clone());
        computed
    }
}

/// Every statistic of one child-star block, bound to its
/// (record, source, parameter, block) key.
///
/// The block is the unit the card-4 finding asks for: its **own** arm set, star
/// size, dimension, little co-group, cocycle class and target source, none of
/// which may be read off another block or off the probe as a whole.
struct BlockStat {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    label: &'static str,
    parameter: Rat,
    /// Position in `result.blocks()`.
    index: usize,
    /// `block.points().len()`.
    star_size: usize,
    /// `Σ point.arm_count()`: `block.arm_count()`, recomputed here.
    arm_count: usize,
    /// `block.arm_indices()`: the parent arms of this child star.
    arm_indices: Vec<usize>,
    block_dimension: u32,
    little_co_group: usize,
    cocycle_trivial: bool,
    /// `classify_block(block.targets())`.
    source: BlockSource,
    /// `(dimension, multiplicity)` of every reported target.
    terms: Vec<(u8, u32)>,
    /// Whether the block carries the reference folded point `t . q1`, decided by
    /// the exact child reciprocal-lattice equivalence `same_mod`.
    carries_reference: bool,
    /// Whether one of the block's folded points is the child Gamma point, i.e.
    /// one of its arms reaches `Gamma` at this parameter.
    carries_gamma: bool,
    /// Targets read from a pinned child row.
    stored_targets: usize,
    /// Targets of the block.
    target_count: usize,
    /// `block.star_size()`, for the cross-check against `star_size`.
    declared_star_size: usize,
    /// `block.arm_count()`, for the cross-check against `arm_count`.
    declared_arm_count: usize,
    /// How many of the block's points were compared for little-co-group order
    /// **and** class agreement; must equal `star_size`.
    point_checks: usize,
    /// How many of those comparisons disagreed.  Disagreement is impossible for
    /// conjugate points and is a gate violation, not a warning: the block's
    /// reported co-group and class would otherwise describe one point of the
    /// star only.
    point_disagreements: usize,
    /// The block's own representative folded coordinate `q`.
    ///
    /// Kept so the gate can call the two point-level entry points again on the
    /// **block's** point instead of trusting the value the statistics were built
    /// with: a wrong point, a wrong child group or a stale memo entry then shows
    /// up as a disagreement.  It is not written to `--output-blocks`.
    representative_point: Vec3R,
}

/// The engine's own dimension bookkeeping of one probe, absent when the engine
/// did not answer at all.
struct BlockDimensions {
    /// `result.blocks().len()`.
    reported_blocks: usize,
    /// `result.covered_dimension()`.
    covered: u32,
    /// `result.parent_dimension()`.
    parent: u32,
    /// `parent / table.dimension`: the number of parent arms the probe carried,
    /// derived from the engine's own dimension identity rather than from the
    /// block list.
    arm_total: usize,
}

/// A borrowed view of [`GammaParameters`]: the `(parameter, arms)` entries and
/// the arms that are at Gamma at every parameter.
type GammaArms<'a> = (&'a [(Rat, Vec<usize>)], &'a [usize]);

/// What one `(record, label)` contributed to the card-3 Gamma report.
#[derive(Default)]
struct GammaReport {
    /// `(parameter, arms)` entries of [`FullStarPartition::gamma_parameters`].
    entries: usize,
    /// Entries whose parameter is a full-star boundary.
    contained: usize,
    /// Entries whose parameter is **not** a boundary and whose reaching arm is
    /// fixed exactly by the whole child point group (verified pointwise).
    exceptional: usize,
    /// Arms whose folded direction is exactly zero: at Gamma at every parameter,
    /// reported separately by the partition and counted here.
    zero_arms: usize,
    /// The distinct Gamma parameters, as report keys.
    shapes: BTreeMap<Vec<(String, Vec<usize>)>, usize>,
    /// Up to [`GAMMA_WITNESSES`] exceptional witnesses, for the printout.
    witnesses: Vec<String>,
}

/// What the full-star recount pass measured over one record.
#[derive(Default)]
struct RecountReport {
    /// The `(record, label, parameter)` probes the card-3 partition adds beyond
    /// the reference candidate set.
    probes: usize,
    /// Blocks classified in the pass.
    blocks: usize,
    /// Blocks by [`BlockSource`].
    sources: [usize; 3],
    /// Blocks whose **own** cocycle is non-trivial, by [`BlockSource`].
    non_trivial: [usize; 3],
    /// Little co-group order histogram of the non-trivial blocks.
    non_trivial_orders: BTreeMap<usize, usize>,
    /// Every distinct block geometry seen (the canonical partition of the arm
    /// list into child stars, by parent arm index).
    geometries: BTreeSet<Vec<Vec<usize>>>,
    /// `(record, label)` pairs whose geometry at an added parameter differs from
    /// the geometry at the generic sample `t = 1/7`.
    changed_pairs: BTreeSet<(usize, &'static str)>,
    /// One fingerprint row per probed parameter: the `--output-recount`
    /// artifact the card-4 review asked for (its row count is asserted against
    /// [`Self::probes`], so "the recount's set equality and its 28 geometries"
    /// can be checked from a file instead of only from an aggregate).
    rows: Vec<RecountRow>,
    /// Failures of this pass; they are gate violations.
    failures: Vec<String>,
}

/// One row of the `--output-recount` artifact: the probe key plus its
/// block-level geometry fingerprint.
#[derive(Debug, Clone, PartialEq, Eq)]
struct RecountRow {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    /// Owned so a row read back from disk is the same type as the row written.
    label: String,
    parameter: Rat,
    /// `(arm indices, arm count, block dimension, source)` per block, in block
    /// order.
    blocks: Vec<(Vec<usize>, usize, u32, BlockSource)>,
}

impl RecountRow {
    /// The row as TSV, exactly as the writer emits it.
    fn to_tsv(&self) -> String {
        let blocks = self
            .blocks
            .iter()
            .map(|(arms, arm_count, dimension, source)| {
                format!(
                    "{}|{arm_count}|{dimension}|{}",
                    arms.iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                    source.label()
                )
            })
            .collect::<Vec<_>>()
            .join(";");
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.ordinal,
            self.parent_sg,
            self.child_sg,
            self.label,
            self.parameter,
            self.blocks.len(),
            blocks
        )
    }

    /// Parse one written row back, field by field (never by re-running the
    /// writer): the card-4 audit's lesson is that a row count alone does not
    /// notice a corrupted row.
    fn parse(line: &str, row: usize) -> Result<Self, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 7 {
            return Err(format!("{row}: {} field(s), expected 7", fields.len()));
        }
        let ordinal = fields[0]
            .parse::<usize>()
            .map_err(|error| format!("{row}: ordinal {}: {error}", fields[0]))?;
        let parent_sg = fields[1]
            .parse::<u8>()
            .map_err(|error| format!("{row}: parent space group {}: {error}", fields[1]))?;
        let child_sg = fields[2]
            .parse::<u8>()
            .map_err(|error| format!("{row}: child space group {}: {error}", fields[2]))?;
        let parameter = parse_parameter(fields[4])
            .map_err(|error| format!("{row}: parameter: {error}"))?;
        let block_count = fields[5]
            .parse::<usize>()
            .map_err(|error| format!("{row}: block count {}: {error}", fields[5]))?;
        let mut blocks = Vec::new();
        if !fields[6].is_empty() {
            for entry in fields[6].split(';') {
                let parts: Vec<&str> = entry.split('|').collect();
                if parts.len() != 4 {
                    return Err(format!("{row}: block {entry:?} has {} field(s)", parts.len()));
                }
                let arms = if parts[0].is_empty() {
                    Vec::new()
                } else {
                    parts[0]
                        .split(',')
                        .map(|value| {
                            value
                                .parse::<usize>()
                                .map_err(|error| format!("{row}: arm {value}: {error}"))
                        })
                        .collect::<Result<Vec<_>, _>>()?
                };
                let arm_count = parts[1]
                    .parse::<usize>()
                    .map_err(|error| format!("{row}: arm count {}: {error}", parts[1]))?;
                let dimension = parts[2]
                    .parse::<u32>()
                    .map_err(|error| format!("{row}: dimension {}: {error}", parts[2]))?;
                let source = match parts[3] {
                    "stored" => BlockSource::Stored,
                    "constructed" => BlockSource::Constructed,
                    "mixed" => BlockSource::Mixed,
                    other => return Err(format!("{row}: block source {other:?}")),
                };
                blocks.push((arms, arm_count, dimension, source));
            }
        }
        if blocks.len() != block_count {
            return Err(format!(
                "{row}: {block_count} block(s) announced but {} written",
                blocks.len()
            ));
        }
        Ok(Self {
            ordinal,
            parent_sg,
            child_sg,
            label: fields[3].to_string(),
            parameter,
            blocks,
        })
    }
}

/// What the R6.7 card-5 `--domain-sweep` measured over one record.
///
/// Every counter is a checkable object: `intervals == comparisons +
/// failed_intervals` is the count conservation the card asks for, so an interval
/// that silently disappeared cannot be hidden behind a green gate.
#[derive(Default)]
struct SweepReport {
    /// `(record, label)` pairs swept.
    pairs: usize,
    /// Open intervals of `probe_parameters(parent)` over all pairs.
    intervals: usize,
    /// Production decompositions run at interior points: two per interval that
    /// reached the engine.
    points: usize,
    /// Intervals whose two interior points were answered **and** whose geometry
    /// and representation comparison reported no failure.
    comparisons: usize,
    /// Intervals that produced at least one failure; every one of them has a
    /// message in [`Self::failures`].
    failed_intervals: usize,
    /// Matched block pairs / matched target pairs.
    blocks: usize,
    targets: usize,
    /// Interior blocks by provenance ([`BlockSource`] order).  Measured, not
    /// assumed: on this corpus no interior point of any interval has a stored
    /// target, which is what makes the stored/constructed matching a synthetic
    /// branch here (see the module tests).
    sources: [usize; 3],
    /// Matched target pairs whose multiplicities disagreed.
    multiplicity_mismatches: usize,
    /// Interior blocks that contain two targets of the same dimension with
    /// different multiplicities: the blocks on which the card-5
    /// "swap two distinct one-dimensional targets" mutation could be observed on
    /// this corpus at all.  Measured: **none** -- every one of the 6,264
    /// multi-target interior blocks has the shape `[(1, 1), (1, 1)]`, two
    /// distinct one-dimensional targets with **equal** multiplicity, so a swap
    /// between them is a no-op on this corpus.  That exact case is therefore
    /// pinned by the module test
    /// `the_matcher_uses_characters_and_not_a_sorted_term_table` (which is also
    /// what the sorted-list mutation of the comparison fails), and the
    /// `multi_shapes` histogram below is what makes the corpus's shape visible.
    same_dimension_swaps: usize,
    /// The target shapes (ascending `(dimension, multiplicity)` lists) of the
    /// interior blocks that carry more than one target.
    multi_shapes: BTreeMap<Vec<(u8, u32)>, usize>,
    /// Matchings that were not one-to-one: a target with no partner, two
    /// candidates at tolerance, or an unmatched dimension.
    ambiguities: usize,
    /// Targets per interior block, and little co-group order at the seed points:
    /// the two shapes the sweep has to be able to exercise.
    target_counts: BTreeMap<usize, usize>,
    orders: BTreeMap<usize, usize>,
    /// Matched block pairs whose **reduced** representative point
    /// ([`FullStarBlock::stored_k`]) differs between the two interior points:
    /// measured **every** pair, because the folded point moves with the parameter
    /// and its representative moves with it.  This is the volume fact that makes
    /// `stored_k` unusable as a cross-parameter key; it does **not** measure the
    /// discontinuous reduction walls the card-5 math review found (P1-2: 96
    /// strictly inside intervals on a denser grid, witness ordinal 13719 SG 225
    /// `SM1` at `t = 1/6`).  The card-5 comparison is on the gauge-unified
    /// content, which is invariant to both.
    representative_shifts: usize,
    /// Exact-fixity checks of the card-1 gauge's hypothesis: for every matched
    /// block, every aligned little-group operation and both parameters, the
    /// rotation must fix the seed arm's point **exactly**, not only modulo the
    /// child reciprocal lattice (that weaker condition is what the little group
    /// is built from).  A non-exact fixity is a partition boundary, so a
    /// mismatch here is both a partition error and a broken gauge hypothesis.
    little_fixity_checks: usize,
    /// The subset of those checks whose rotation is **not** the identity.  The
    /// count pin alone does not bind the operand: replacing the checked rotation
    /// by the identity keeps the loop running 549,184 times and stays green
    /// (verification of `b2d39f3`, P0).  Every aligned little co-group contains
    /// exactly one identity, so this count is predicted independently as
    /// `2 x (sum(order x blocks) - block pairs)`, and a constant non-identity
    /// rotation cannot hide either: it fails the exact-fixity test at the first
    /// point it does not fix, which the mismatch counter reports.
    little_fixity_non_identity_checks: usize,
    little_fixity_mismatches: usize,
    /// Probe parameters that are **not** full-star boundaries, i.e. the entries
    /// the parent's own formal parameter set contributes to the union.  Measured
    /// (corpus: **0**), because the parent's set is `{0, 1/2}` and both are
    /// already eighth-grid boundaries -- so the union this sweep is built on is
    /// implemented but not exercised by this corpus, and that is a scope
    /// statement rather than a reason to drop it.
    parent_only_parameters: usize,
    /// The worst gauge-unified character score seen (`1.0` when nothing was
    /// compared); the corpus agreement is at machine precision.
    worst_score: f64,
    /// The R6.7 card-6 boundary pass: the number of `probe_parameters(parent)`
    /// entries at which the **production** `subduce_line_at_parameter` was run,
    /// and how many of them were answered / failed.  Every parameter is counted
    /// before the engine is asked, so `boundary_points == boundary_successes +
    /// boundary_failures` holds by construction and a lost item is impossible;
    /// the gate asserts it and pins all three.
    boundary_points: usize,
    boundary_successes: usize,
    boundary_failures: usize,
    /// `Σ result.parent_dimension()` and `Σ result.covered_dimension()` over the
    /// answered boundary parameters: two numbers taken **from the returned
    /// decomposition**, so a boundary pass that stops calling the engine -- or
    /// calls it and drops the answer -- moves them while every count stays where
    /// it is.  They must be equal (the engine's own dimension identity) and are
    /// pinned.
    boundary_parent_dimension: u64,
    boundary_covered_dimension: u64,
    /// The R6.7 card-6 independent recomputation of the same pairs' arm geometry
    /// from the folded arms and the child reciprocal lattice.
    geometry: GeometryReport,
    /// Failures of this pass; they are gate violations.
    failures: Vec<String>,
}

/// Keep the Gamma exception witnesses bounded: the condition is measured, and
/// the report shows the first few rather than growing with the corpus.
const GAMMA_WITNESSES: usize = 8;

/// What the engine's own [`FullStarBlock`] says about one block, read directly
/// from the block object.
///
/// The card-4 audit found that the gate's only provenance check was a relation
/// among values written by one `block_stat` call, so a mutation that swapped
/// `source`, `stored_targets`, `target_count` and `terms` between two blocks of
/// the same probe -- for every ordinal except the two pinned witnesses -- left
/// the gate green (P1).  Keeping the engine's own reading and comparing the
/// recorded statistics against it binds every block statistic to the block it
/// claims to describe; the same reading covers the per-target (dimension,
/// multiplicity) pairs, which the audit could also bump by one unnoticed.
#[derive(Debug, Clone, PartialEq, Eq)]
struct EngineBlock {
    star_size: usize,
    arm_count: usize,
    arm_indices: Vec<usize>,
    block_dimension: u32,
    /// `Σ dimension × multiplicity` through the engine's **own** bookkeeping
    /// (`FullStarBlock::little_dimension`), which is computed from the solved
    /// targets rather than re-summed here: the independent route to the
    /// per-target terms the verification review asked for.
    little_dimension: u32,
    source: BlockSource,
    stored_targets: usize,
    target_count: usize,
    terms: Vec<(u8, u32)>,
    /// The representative the engine reports, and the folded coordinate it was
    /// read at.  Without these a post-construction move of the recorded
    /// `representative_point` was invisible (verification review F4).
    representative: usize,
    q: Vec3R,
    /// Whether the block carries the reference folded point / reaches Gamma,
    /// recomputed here from the engine's own points.  Both were unbound before
    /// (verification review F1/F5), and the unbound `carries_reference` was what
    /// let a swap move the headline corrected table from 4138/192/0 to
    /// 4146/184/0 with a green gate.
    carries_reference: bool,
    carries_gamma: bool,
}

/// Read one engine block's own statistics.
///
/// Every lattice query propagates its error: an overflow or an inconsistent
/// lattice used to be swallowed into `false` on **both** sides of the compare,
/// which would have turned a real failure into a matching pair of `false`s
/// (card-5 audit F7).
fn engine_block(
    block: &FullStarBlock,
    reference_point: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<EngineBlock, String> {
    let mut carries_reference = false;
    let mut carries_gamma = false;
    for point in block.points() {
        carries_gamma |= child_reciprocal
            .contains(point.q())
            .map_err(|error| format!("the Gamma containment of {}: {error}", point.q()))?;
        carries_reference |= child_reciprocal
            .same_mod(point.q(), reference_point)
            .map_err(|error| format!("the reference comparison of {}: {error}", point.q()))?;
    }
    Ok(EngineBlock {
        star_size: block.points().len(),
        arm_count: block.points().iter().map(|point| point.arm_count()).sum(),
        arm_indices: block.arm_indices(),
        block_dimension: block.block_dimension(),
        little_dimension: block.little_dimension(),
        source: classify_block(block.targets()),
        stored_targets: block
            .targets()
            .iter()
            .filter(|target| target.irnumber.is_some())
            .count(),
        target_count: block.targets().len(),
        terms: block
            .targets()
            .iter()
            .map(|target| (target.dimension, target.multiplicity))
            .collect(),
        representative: block.representative(),
        q: *block.q(),
        carries_reference,
        carries_gamma,
    })
}

/// One probe of the census.
struct Probe {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    label: &'static str,
    parameter: Rat,
    parent_order: usize,
    parent_added: usize,
    child_order: usize,
    class: TargetClass,
    parameter_kind: Option<ParameterKind>,
    content: Option<u32>,
    detail: String,
    /// One [`BlockStat`] per `result.blocks()` entry, in block order; empty when
    /// the engine did not answer.
    blocks: Vec<BlockStat>,
    /// The engine's own reading of the same blocks, in the same order.  The gate
    /// compares every recorded statistic against it, so a corruption of the
    /// census-side bookkeeping cannot pass on its own.
    engine_blocks: Vec<EngineBlock>,
    /// The engine's own dimension bookkeeping, absent for an unsupported or
    /// failed probe.
    dimensions: Option<BlockDimensions>,
    /// `FullStarPartition::arms.len()` for this `(record, label)`, when the
    /// card-3 partition was built; the gate compares it with the arm total the
    /// engine's dimension identity gives.
    partition_arms: Option<usize>,
}

/// One isotropy record that carries at least one parametric-k row.
struct Record {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    subgroup: isotropy::IsotropySubgroup,
    /// Frozen source label and the pinned frequency of its row, so the anchor
    /// probe can be compared with the pinned table instead of only printed.
    labels: Vec<(&'static str, u16)>,
}

/// What one record contributed to the census.
struct RecordReport {
    probes: Vec<Probe>,
    /// `(checks, mismatches)` of the child-side grid cross-check.
    child_grid: (usize, usize),
    /// `(checks, mismatches)` of the child-side two-algorithm step comparison.
    child_algorithms: (usize, usize),
    /// The record's candidate parameters (parent union child).
    child_parameters: Vec<Rat>,
    /// The card-3 Gamma-reaching report of this record's `(record, label)` pairs.
    gamma: GammaReport,
    /// The card-3 full-star recount pass of this record.
    recount: RecountReport,
    /// The R6.7 card-5 interval sweep of this record.
    sweep: SweepReport,
    failures: Vec<String>,
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => code,
        Err(message) => {
            eprintln!("line_domain_census: {message}");
            ExitCode::from(3)
        }
    }
}

fn run() -> Result<ExitCode, String> {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let gate = arguments.iter().any(|argument| argument == "--gate");
    let require_covered = arguments.iter().any(|argument| argument == "--require-covered");
    // The full-star recount pass is a gate pass: `--gate` runs it, and
    // `--full-star-recount` runs it (and its assertions) on its own, with the
    // same failure semantics -- a flag that runs a check and then exits 0 on its
    // failures would be a failure path that does not fail.
    let full_star_recount = arguments.iter().any(|argument| argument == "--full-star-recount");
    // R6.7 card 5: the sweep is its own pass.  `--gate` alone does not run it
    // (the acceptance command names both flags), but when it runs its failures
    // are gate violations and exit 1 under `--gate` -- and under the flag alone,
    // like `--full-star-recount`, because a check that reports and then exits 0
    // would be a failure path that does not fail.
    let domain_sweep = arguments.iter().any(|argument| argument == "--domain-sweep");
    let sequential = arguments.iter().any(|argument| argument == "--sequential");
    if arguments.iter().any(|argument| argument == "--help" || argument == "-h") {
        print!("{USAGE}");
        return Ok(ExitCode::SUCCESS);
    }
    let output = arguments
        .iter()
        .position(|argument| argument == "--output")
        .map(|index| {
            arguments
                .get(index + 1)
                .cloned()
                .ok_or_else(|| "--output needs a path".to_string())
        })
        .transpose()?;
    let output_blocks = arguments
        .iter()
        .position(|argument| argument == "--output-blocks")
        .map(|index| {
            arguments
                .get(index + 1)
                .cloned()
                .ok_or_else(|| "--output-blocks needs a path".to_string())
        })
        .transpose()?;
    let output_recount = arguments
        .iter()
        .position(|argument| argument == "--output-recount")
        .map(|index| {
            arguments
                .get(index + 1)
                .cloned()
                .ok_or_else(|| "--output-recount needs a path".to_string())
        })
        .transpose()?;
    // Asking for the artifact asks for the pass that produces it, with the same
    // failure semantics as `--full-star-recount`.
    let recount = gate || full_star_recount || output_recount.is_some();

    let domains = source_domains()?;
    // Card 6: a subgroup record the census cannot read is a counted failure, not
    // a silent skip, so it starts the violation list.
    let (records, record_failures) = records()?;
    let official = official_line_parameter().map_err(|error| error.to_string())?;
    let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1).map_err(|e| e.to_string())?;
    let cache = PointClassCache::default();

    let per_record: Vec<RecordReport> = if sequential {
        records
            .iter()
            .map(|record| {
                probe_record(
                    record,
                    &domains,
                    official,
                    generic,
                    recount,
                    domain_sweep,
                    &cache,
                )
            })
            .collect()
    } else {
        records
            .par_iter()
            .map(|record| {
                probe_record(
                    record,
                    &domains,
                    official,
                    generic,
                    recount,
                    domain_sweep,
                    &cache,
                )
            })
            .collect()
    };
    let mut probes: Vec<Probe> = Vec::new();
    let mut probe_errors: Vec<String> = Vec::new();
    let mut child_grid = (0usize, 0usize);
    let mut child_algorithms = (0usize, 0usize);
    let mut child_union: Vec<Rat> = Vec::new();
    let mut gamma = GammaReport::default();
    let mut recount_report = RecountReport::default();
    let mut sweep_report = SweepReport {
        worst_score: 1.0,
        ..SweepReport::default()
    };
    for (record, report) in records.iter().zip(per_record) {
        if report.probes.is_empty() {
            probe_errors.push(format!("ordinal {} produced no probe", record.ordinal));
        }
        probe_errors.extend(report.failures);
        child_grid.0 += report.child_grid.0;
        child_grid.1 += report.child_grid.1;
        child_algorithms.0 += report.child_algorithms.0;
        child_algorithms.1 += report.child_algorithms.1;
        for parameter in report.child_parameters {
            if !child_union.contains(&parameter) {
                child_union.push(parameter);
            }
        }
        merge_gamma(&mut gamma, report.gamma);
        merge_recount(&mut recount_report, report.recount);
        merge_sweep(&mut sweep_report, report.sweep);
        probes.extend(report.probes);
    }
    child_union = sorted_parameters(child_union);

    if let Some(path) = &output {
        let mut writer = BufWriter::new(
            std::fs::File::create(path).map_err(|error| format!("cannot create {path}: {error}"))?,
        );
        writeln!(
            writer,
            "ordinal\tparent_sg\tchild_sg\tlabel\tparameter\tparent_order\tparent_added\t\
             child_order\tclass\tparameter_kind\tcontent\tdetail"
        )
        .map_err(|error| error.to_string())?;
        for probe in &probes {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                probe.ordinal,
                probe.parent_sg,
                probe.child_sg,
                probe.label,
                probe.parameter,
                probe.parent_order,
                probe.parent_added,
                probe.child_order,
                probe.class.label(),
                match probe.parameter_kind {
                    Some(ParameterKind::LineIrrep) => "line",
                    Some(ParameterKind::Formal) => "formal",
                    None => "-",
                },
                probe
                    .content
                    .map_or_else(|| "-".to_string(), |value| value.to_string()),
                probe.detail,
            )
            .map_err(|error| error.to_string())?;
        }
    }

    // The per-block table.  `--output` above is deliberately unchanged; this is
    // a new artifact with one row per child-star block, so a block's own source,
    // co-group and cocycle are readable without re-running the engine.  The
    // emission loop counts its rows on every run and the written file is counted
    // back from disk, so a writer that drops rows cannot pass either check.
    let mut block_rows = 0usize;
    {
        let mut writer = match &output_blocks {
            None => None,
            Some(path) => Some(BufWriter::new(
                std::fs::File::create(path)
                    .map_err(|error| format!("cannot create {path}: {error}"))?,
            )),
        };
        if let Some(writer) = writer.as_mut() {
            writeln!(writer, "{BLOCK_TSV_HEADER}").map_err(|error| error.to_string())?;
        }
        for probe in &probes {
            for block in &probe.blocks {
                if let Some(writer) = writer.as_mut() {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        block.ordinal,
                        block.parent_sg,
                        block.child_sg,
                        block.label,
                        block.parameter,
                        block.index,
                        block.star_size,
                        block.arm_count,
                        block
                            .arm_indices
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                        block.block_dimension,
                        block.little_co_group,
                        if block.cocycle_trivial { "trivial" } else { "non-trivial" },
                        block.source.label(),
                        block
                            .terms
                            .iter()
                            .map(|(dimension, multiplicity)| format!("{dimension}x{multiplicity}"))
                            .collect::<Vec<_>>()
                            .join(","),
                        block.carries_reference,
                        block.carries_gamma,
                    )
                    .map_err(|error| error.to_string())?;
                }
                block_rows += 1;
            }
        }
        if let Some(writer) = writer.as_mut() {
            writer.flush().map_err(|error| error.to_string())?;
        }
    }
    let mut block_file_failures: Vec<String> = Vec::new();
    let block_file_rows = match &output_blocks {
        None => None,
        Some(path) => {
            let text = std::fs::read_to_string(path)
                .map_err(|error| format!("cannot read back {path}: {error}"))?;
            // One header line plus one line per block; `str::lines` does not
            // invent a last empty line, so the count is the file's row count.
            // Every data row is also **parsed** and compared with the statistic
            // it came from: a row count alone did not notice the card-4 audit's
            // `star_size + 1` mutation, which left the file the right length.
            // A separate literal, deliberately **not** the writer's own const: two
            // independent definition sources are what makes a changed header
            // visible (with the shared const, swapping two names moved both sides
            // together and stayed green).
            const EXPECTED_HEADER: &str = "ordinal\tparent_sg\tchild_sg\tlabel\tparameter\t\
index\tstar_size\tarm_count\tarm_indices\tblock_dimension\tlittle_co_group\t\
cocycle_trivial\tblock_source\tterms\tcarries_reference\tcarries_gamma";
            if text.lines().next() != Some(EXPECTED_HEADER) {
                block_file_failures.push(format!(
                    // `unwrap_or("")` is a message default, not a swallowed
                    // result: the branch itself is the failure, and an empty
                    // header cannot equal the literal above.
                    "the --output-blocks header is {:?}, expected {EXPECTED_HEADER:?}",
                    text.lines().next().unwrap_or("")
                ));
            }
            let mut written = text.lines().skip(1);
            let mut rows = 0usize;
            for (probe, block) in probes
                .iter()
                .flat_map(|probe| probe.blocks.iter().map(move |block| (probe, block)))
            {
                let Some(line) = written.next() else {
                    block_file_failures.push(format!(
                        "the --output-blocks file ends after {rows} row(s), before block \
                         ({} {} t={} block {})",
                        probe.ordinal, probe.label, probe.parameter, block.index
                    ));
                    break;
                };
                rows += 1;
                if let Err(error) = check_block_row(line, block, rows) {
                    block_file_failures.push(format!(
                        "the --output-blocks row for ({} {} t={} block {}): {error}",
                        probe.ordinal, probe.label, probe.parameter, block.index
                    ));
                }
            }
            if let Some(extra) = written.next() {
                block_file_failures.push(format!(
                    "the --output-blocks file has a data row beyond the {rows} collected block(s):                      {:?}",
                    &extra[..extra.len().min(60)]
                ));
            }
            Some(rows)
        }
    };

    // The per-recount-probe artifact, the residual the card-4 review left open
    // (there was no per-probe product, so the recount's set equality and its 28
    // geometries could only be re-checked by re-running the pass).  The row
    // count is asserted against the pass's own probe count and every row is
    // parsed back and compared field by field.
    let mut recount_rows_written = 0usize;
    let mut recount_file_failures: Vec<String> = Vec::new();
    let recount_file_rows = match &output_recount {
        None => None,
        Some(path) => {
            let mut writer = BufWriter::new(
                std::fs::File::create(path)
                    .map_err(|error| format!("cannot create {path}: {error}"))?,
            );
            writeln!(writer, "{RECOUNT_TSV_HEADER}").map_err(|error| error.to_string())?;
            for row in &recount_report.rows {
                writeln!(writer, "{}", row.to_tsv()).map_err(|error| error.to_string())?;
                recount_rows_written += 1;
            }
            writer.flush().map_err(|error| error.to_string())?;
            let text = std::fs::read_to_string(path)
                .map_err(|error| format!("cannot read back {path}: {error}"))?;
            // A separate literal, deliberately **not** the writer's own const:
            // two independent definition sources are what makes a changed header
            // visible (card-4 audit, verification round 2).
            const EXPECTED_HEADER: &str = "ordinal\tparent_sg\tchild_sg\tlabel\tparameter\t\
block_count\tblocks";
            if text.lines().next() != Some(EXPECTED_HEADER) {
                recount_file_failures.push(format!(
                    // A message default, not a swallowed result: the branch is
                    // the failure (see the `--output-blocks` header above).
                    "the --output-recount header is {:?}, expected {EXPECTED_HEADER:?}",
                    text.lines().next().unwrap_or("")
                ));
            }
            let mut rows = 0usize;
            for line in text.lines().skip(1) {
                rows += 1;
                match RecountRow::parse(line, rows) {
                    Ok(parsed) => match recount_report.rows.get(rows - 1) {
                        Some(expected) if *expected == parsed => {}
                        Some(expected) => recount_file_failures.push(format!(
                            "the --output-recount row {rows} is {parsed:?}, expected {expected:?}"
                        )),
                        None => recount_file_failures.push(format!(
                            "the --output-recount file has a row {rows} beyond the {} collected",
                            recount_report.rows.len()
                        )),
                    },
                    Err(error) => {
                        recount_file_failures.push(format!("the --output-recount row {rows}: {error}"))
                    }
                }
            }
            Some(rows)
        }
    };

    let (algorithm_checks, algorithm_mismatches, algorithm_failures) =
        step_algorithm_cross_check(&domains);
    let evidence = CensusEvidence {
        child_grid,
        child_algorithms,
        child_union: &child_union,
        algorithm_checks,
        algorithm_mismatches,
        block_rows,
        block_file_rows,
        gamma: &gamma,
        recount: &recount_report,
        recount_ran: recount,
        recount_rows_written,
        recount_file_rows,
        sweep: &sweep_report,
        sweep_ran: domain_sweep,
    };
    report(
        &domains,
        &records,
        &probes,
        &probe_errors,
        &evidence,
        &gamma,
        &recount_report,
    );

    let mut violations: Vec<String> = probe_errors;
    // Card 6: the records the walk could not census join the violations, so they
    // are visible in the gate rather than only in the corpus total.
    violations.extend(record_failures.iter().cloned());
    violations.extend(algorithm_failures);
    violations.extend(block_file_failures);
    violations.extend(recount_file_failures.iter().cloned());
    violations.extend(recount_report.failures.iter().cloned());
    violations.extend(sweep_report.failures.iter().cloned());
    check_invariants(&domains, &records, &probes, &evidence, &mut violations);

    let unsupported = probes
        .iter()
        .filter(|probe| probe.class == TargetClass::Unsupported)
        .count();
    let errors = probes.iter().filter(|probe| probe.class == TargetClass::Error).count();
    println!(
        "probes={} stored={} constructed={} unsupported={} errors={}",
        probes.len(),
        probes.iter().filter(|p| p.class == TargetClass::Stored).count(),
        probes.iter().filter(|p| p.class == TargetClass::Constructed).count(),
        unsupported,
        errors
    );
    for violation in &violations {
        eprintln!("census violation: {violation}");
    }
    if gate {
        // Boundary census: the projective class at the child's **own** exceptional
        // parameters, where the little co-group is strictly larger than the exact
        // stabiliser of the direction and the domain theorem therefore does not
        // apply.  The point is to measure the boundary rather than assume it: the
        // engine needs a non-coboundary family exactly at the parameters counted
        // as NON-TRIVIAL here.
        //
        // Two tables are printed.  The **corrected** one (R6.7 card 4) pairs the
        // non-trivial class of the reference folded point with the *own* source
        // and little co-group of the block that carries that point.  The
        // **withdrawn** one is the old probe-level table, whose `Constructed`
        // column only ever meant "some target of some block of this probe was
        // constructed" -- the mixing that made "292/44/4 = 340 non-trivial
        // constructed classes" wrong (witness ordinal 13688).  It stays in the
        // printout, labelled, so its numbers can be compared with the corrected
        // ones on the same run; nothing is pinned to it.
        let mut legacy: BTreeMap<(usize, bool, TargetClass), usize> = BTreeMap::new();
        let mut corrected: BTreeMap<(usize, BlockSource), usize> = BTreeMap::new();
        // Non-trivial-class probes by source, in the two conventions.
        let mut legacy_totals = [0usize; 3];
        let mut corrected_totals = [0usize; 3];
        let mut boundary_probes = 0usize;
        let mut boundary_errors = 0usize;
        // A non-trivial probe whose reference point no block carries: the
        // corrected table would be silently short, so it is a counted failure.
        let mut missing_reference_blocks = 0usize;
        let answered: BTreeMap<String, &Probe> = probes
            .iter()
            .map(|probe| {
                (
                    format!("{}|{}|{}", probe.ordinal, probe.label, probe.parameter),
                    probe,
                )
            })
            .collect();
        for record in &records {
            let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup) else {
                boundary_errors += 1;
                violations.push(format!(
                    "ordinal {}: the boundary census cannot build the embedding",
                    record.ordinal
                ));
                continue;
            };
            let Ok(reciprocal) = reciprocal_lattice(record.child_sg) else {
                boundary_errors += 1;
                violations.push(format!(
                    "ordinal {}: the boundary census found no reciprocal lattice for child #{}",
                    record.ordinal, record.child_sg
                ));
                continue;
            };
            let Ok(rotations) = rotation_set(record.child_sg) else {
                boundary_errors += 1;
                violations.push(format!(
                    "ordinal {}: the boundary census found no rotation set for child #{}",
                    record.ordinal, record.child_sg
                ));
                continue;
            };
            for (label, _) in &record.labels {
                let Some(table) = line_table(record.parent_sg, label) else {
                    boundary_errors += 1;
                    violations.push(format!(
                        "ordinal {}: the boundary census found no frozen table for SG {} {label}",
                        record.ordinal, record.parent_sg
                    ));
                    continue;
                };
                let Ok(candidates) = child_exceptional_parameters(&embedding, table) else {
                    boundary_errors += 1;
                    violations.push(format!(
                        "ordinal {} {label}: the boundary census cannot enumerate the child's \
                         exceptional parameters",
                        record.ordinal
                    ));
                    continue;
                };
                let Some(direction) = line_direction(table) else {
                    boundary_errors += 1;
                    violations.push(format!(
                        "ordinal {} {label}: the boundary census cannot parse the frozen direction",
                        record.ordinal
                    ));
                    continue;
                };
                let Ok(folded) = fold_wave_vector(embedding.transform(), &direction) else {
                    boundary_errors += 1;
                    violations.push(format!(
                        "ordinal {} {label}: the boundary census cannot fold the direction",
                        record.ordinal
                    ));
                    continue;
                };
                for (parameter, _) in &candidates {
                    let Ok(order) =
                        little_co_group_order(&reciprocal, &folded, &rotations, *parameter)
                    else {
                        boundary_errors += 1;
                        continue;
                    };
                    boundary_probes += 1;
                    let key = format!("{}|{}|{}", record.ordinal, label, parameter);
                    let probe = answered.get(&key).copied();
                    if probe.is_none() {
                        boundary_errors += 1;
                        violations.push(format!(
                            "ordinal {} {label} t={parameter}: the boundary census has no probe for \
                             a child-exceptional parameter",
                            record.ordinal
                        ));
                    }
                    let class = probe.map_or(TargetClass::Error, |probe| probe.class);
                    match child_cocycle_is_a_coboundary(record.child_sg, &folded, *parameter) {
                        Ok(trivial) => {
                            *legacy.entry((order, trivial, class)).or_insert(0usize) += 1;
                            if trivial {
                                continue;
                            }
                        }
                        Err(error) => {
                            boundary_errors += 1;
                            violations.push(format!(
                                "ordinal {} {label} t={parameter}: the boundary census cannot \
                                 decide the class: {error}",
                                record.ordinal
                            ));
                            continue;
                        }
                    }
                    legacy_totals[match class {
                        TargetClass::Stored => 0,
                        TargetClass::Constructed => 1,
                        TargetClass::Unsupported | TargetClass::Error => 2,
                    }] += 1;
                    match probe.and_then(|probe| {
                        probe.blocks.iter().find(|block| block.carries_reference)
                    }) {
                        Some(block) => {
                            corrected_totals[block.source.index()] += 1;
                            *corrected
                                .entry((block.little_co_group, block.source))
                                .or_insert(0usize) += 1;
                        }
                        None => missing_reference_blocks += 1,
                    }
                }
            }
        }
        println!(
            "withdrawn convention -- projective class at the child's exceptional parameters, \
             probe-level: {boundary_probes} probes, {boundary_errors} error(s)"
        );
        for ((order, trivial, class), count) in &legacy {
            println!(
                "    child order {order} {} engine {class:?}: {count}",
                if *trivial { "trivial     " } else { "NON-TRIVIAL " }
            );
        }
        let non_trivial = legacy_totals.iter().sum::<usize>();
        println!(
            "    non-trivial-class probes: stored={} constructed={} other={} \
             (the withdrawn `292/44/4 = 340` family)",
            legacy_totals[0], legacy_totals[1], legacy_totals[2]
        );
        println!(
            "corrected convention (R6.7 card 4) -- the reference-carrying block's own source: \
             {non_trivial} non-trivial-class probe(s), {missing_reference_blocks} without a \
             reference block"
        );
        for ((order, source), count) in &corrected {
            println!(
                "    child order {order} NON-TRIVIAL block {}: {count}",
                source.label()
            );
        }
        println!(
            "    non-trivial-class blocks: stored={} constructed={} mixed={}",
            corrected_totals[0], corrected_totals[1], corrected_totals[2]
        );
        if missing_reference_blocks > 0 {
            violations.push(format!(
                "the corrected boundary census lost {missing_reference_blocks} non-trivial-class \
                 probe(s) whose reference point no block carries"
            ));
        }
        // The corrected table is the card-4 deliverable, so it is **pinned**, not
        // only printed: the verification review's `carries_reference` swap rewrote
        // it to 4146/184/0 (per-child-order 2/4/8/16 rows moved) with a green gate.
        if non_trivial != 4_330 || corrected_totals != [4_138, 192, 0] {
            violations.push(format!(
                "the corrected boundary table is {}/{}/{} over {non_trivial} non-trivial-class \
                 probe(s), expected 4138/192/0 over 4330",
                corrected_totals[0], corrected_totals[1], corrected_totals[2]
            ));
        }
        for (order, expected) in [(4usize, [3_068usize, 192, 0]), (8, [936, 0, 0]), (16, [134, 0, 0])] {
            let mut found = [0usize; 3];
            for ((found_order, source), count) in &corrected {
                if *found_order == order {
                    found[source.index()] += count;
                }
            }
            if found != expected {
                violations.push(format!(
                    "the corrected boundary table at child order {order} is {}/{}/{}, expected \
                     {}/{}/{}",
                    found[0], found[1], found[2], expected[0], expected[1], expected[2]
                ));
            }
        }
    }

    // `--domain-sweep` joins the flags that fail on their own failures: a pass
    // that runs a check, reports a violation and then exits 0 would be a failure
    // path that does not fail (the same rule `--full-star-recount` follows).
    if gate || require_covered || full_star_recount || domain_sweep {
        if !violations.is_empty() {
            println!("gate: FAILED ({} violation(s))", violations.len());
            return Ok(ExitCode::from(1));
        }
        if require_covered && unsupported > 0 {
            println!("gate: FAILED ({unsupported} unsupported probe(s), --require-covered)");
            return Ok(ExitCode::from(1));
        }
        println!("gate: ok");
    }
    Ok(ExitCode::SUCCESS)
}

/// The exact domain of every frozen source, checked against its own table.
fn source_domains() -> Result<BTreeMap<(u8, &'static str), ParentDomain>, String> {
    let mut domains = BTreeMap::new();
    for table in W_LITTLE_CHARACTERS {
        let domain = parent_domain(table).map_err(|error| {
            format!("SG {} {}: parameter domain: {error}", table.space_group, table.label)
        })?;
        if domain.generic_order != table.operations.len() {
            return Err(format!(
                "SG {} {}: the generic stabiliser has order {} but the frozen table lists {} \
                 operation(s)",
                table.space_group,
                table.label,
                domain.generic_order,
                table.operations.len()
            ));
        }
        domains.insert((table.space_group, table.label), domain);
    }
    Ok(domains)
}

/// Every isotropy record that carries parametric-k rows, with its labels, and the
/// subgroup records that could not be censused at all.
///
/// The second half is the card-6 replacement for the two silent skips of
/// [`records_of`]: a record the census cannot read is a **counted** failure, so
/// the corpus total and the gate cannot disagree about what was censused.
fn records() -> Result<(Vec<Record>, Vec<String>), String> {
    let mut out = Vec::new();
    let mut failures = Vec::new();
    for sg in 1..=230u8 {
        out.extend(records_of(sg, &mut failures)?);
    }
    Ok((out, failures))
}

/// The same walk, restricted to one parent space group.
///
/// Split out so the card-4 regression tests can build one record without
/// enumerating the whole corpus first; both callers go through the same code, so
/// a test cannot exercise a record the census would not have built.
///
/// `failures` collects the subgroup records the census **cannot** enter.  Both of
/// them used to be a silent `continue` (card-6 finding): a record whose pinned
/// rows cannot be read, and a record whose child space group number does not fit
/// the census's `u8` key.  Dropping either would shrink the corpus without
/// moving any counter except the total, so they are counted failures now.
fn records_of(sg: u8, failures: &mut Vec<String>) -> Result<Vec<Record>, String> {
    let mut out = Vec::new();
    for record in query::irreps_of(sg) {
        if record.spinor || record.subgroups().is_empty() {
            continue;
        }
        let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
            .map_err(|error| format!("SG {sg} {}: {error}", record.ml))?;
        for subgroup in subgroups {
            let rows = match subgroup.other_wave_vector_subduction() {
                Ok(rows) => rows,
                Err(error) => {
                    failures.push(format!(
                        "ordinal {} (SG {sg}): the pinned wave-vector rows cannot be read, so the \
                         record is not censused: {error}",
                        subgroup.ordinal
                    ));
                    continue;
                }
            };
            let child_sg = match u8::try_from(subgroup.record.sg) {
                Ok(child_sg) => child_sg,
                Err(error) => {
                    failures.push(format!(
                        "ordinal {} (SG {sg}): the child space group {} does not fit the census's \
                         space group key: {error}",
                        subgroup.ordinal, subgroup.record.sg
                    ));
                    continue;
                }
            };
            let mut labels: Vec<(&'static str, u16)> = Vec::new();
            for row in rows {
                match labels.iter().find(|(label, _)| *label == row.parent_ml) {
                    Some((_, frequency)) => {
                        if *frequency != row.frequency {
                            return Err(format!(
                                "ordinal {}: label {} carries two pinned frequencies \
                                 ({frequency} and {})",
                                subgroup.ordinal, row.parent_ml, row.frequency
                            ));
                        }
                    }
                    None => labels.push((row.parent_ml, row.frequency)),
                }
            }
            if labels.is_empty() {
                continue;
            }
            out.push(Record {
                ordinal: subgroup.ordinal,
                parent_sg: sg,
                child_sg,
                subgroup,
                labels,
            });
        }
    }
    Ok(out)
}

fn gcd(left: i128, right: i128) -> i128 {
    let (mut a, mut b) = (left.abs(), right.abs());
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    a
}

/// Exact ordering key for the census grid: `Rat` deliberately has no `Ord`, so
/// the parameters are compared through a common denominator.  The census grid is
/// bounded by the source directions (denominators at most 12 in the frozen
/// corpus), so the multiplication cannot overflow; a violation would be a
/// finding, not a silent wrap.
fn sorted_parameters(mut values: Vec<Rat>) -> Vec<Rat> {
    let mut common: i128 = 1;
    for value in &values {
        let divisor = gcd(common, value.denominator());
        common = (common / divisor)
            .checked_mul(value.denominator())
            .expect("the census parameter grid stays small");
    }
    values.sort_by_key(|value| value.numerator() * (common / value.denominator()));
    values
}

/// The parameters of one record: `(probed parameters with the census's child
/// little co-group order there, the child's own candidate parameters)`.
///
/// The order is zero for a parameter the child census does not list, which the
/// probe path reads as "no census-side expectation at this parameter".
type RecordParameters = (Vec<(Rat, usize)>, Vec<Rat>);

/// The partition of one record: the parent's exceptional parameters, the folded
/// child's exceptional parameters, the official anchor and one generic sample.
fn record_parameters(
    record: &Record,
    domain: &ParentDomain,
    table: &'static LittleCharacterTable,
    official: Rat,
    generic: Rat,
) -> Result<RecordParameters, String> {
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup)
        .map_err(|error| format!("ordinal {}: {error}", record.ordinal))?;
    let mut parameters: Vec<(Rat, usize)> = Vec::new();
    let push = |parameter: Rat, order: usize, parameters: &mut Vec<(Rat, usize)>| {
        match parameters.iter_mut().find(|(value, _)| *value == parameter) {
            Some((_, slot)) => *slot = (*slot).max(order),
            None => parameters.push((parameter, order)),
        }
    };
    for entry in &domain.exceptional {
        push(entry.parameter, 0, &mut parameters);
    }
    let folded = child_exceptional_parameters(&embedding, table)
        .map_err(|error| format!("ordinal {} child domain: {error}", record.ordinal))?;
    let child_candidates: Vec<Rat> = folded.iter().map(|(parameter, _)| *parameter).collect();
    // The child order the **census itself** computes at each of its candidate
    // parameters.  It is carried into the probe so the gate can compare it with
    // the order recomputed from `reciprocal_lattice(child_sg)`: without that
    // comparison the child-frame assertions below can only compare that function
    // with itself, and a `reciprocal_lattice` corrupted for one space group is
    // invisible (fifth review round).
    for (parameter, order) in &folded {
        push(*parameter, *order, &mut parameters);
    }
    for parameter in [official, generic] {
        push(parameter, 0, &mut parameters);
    }
    let mut ordered: Vec<(Rat, usize)> = Vec::with_capacity(parameters.len());
    let sorted = sorted_parameters(parameters.iter().map(|(value, _)| *value).collect());
    for value in sorted {
        let order = parameters
            .iter()
            .find(|(parameter, _)| *parameter == value)
            .map_or(0, |(_, order)| *order);
        ordered.push((value, order));
    }
    Ok((ordered, child_candidates))
}

/// Whether the step solver reports "no constraint, or the integrality step one",
/// which is what a rotation that fixes the whole direction must produce.
fn step_is_vacuous_or_one(step: &Option<Rat>) -> bool {
    match step {
        None => true,
        Some(value) => *value == Rat::new(1, 1).expect("one"),
    }
}

/// The three child-frame values one probe needs, or `None` with the reason
/// recorded in `failures`.
///
/// All three are required before anything is probed, so the failure is reported
/// once per record instead of silently turning into an empty probe list.
fn probe_frame(
    record: &Record,
    failures: &mut Vec<String>,
) -> Option<(SubgroupEmbedding, Lattice, Vec<Mat3I>)> {
    let embedding = match SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup) {
        Ok(embedding) => embedding,
        Err(error) => {
            failures.push(format!("ordinal {}: embedding rejected: {error}", record.ordinal));
            return None;
        }
    };
    let reciprocal = match reciprocal_lattice(record.child_sg) {
        Ok(lattice) => lattice,
        Err(error) => {
            failures.push(format!(
                "ordinal {}: child #{} has no reciprocal lattice: {error}",
                record.ordinal, record.child_sg
            ));
            return None;
        }
    };
    let rotations = match rotation_set(record.child_sg) {
        Ok(rotations) => rotations,
        Err(error) => {
            failures.push(format!(
                "ordinal {}: child #{} has no rotation set: {error}",
                record.ordinal, record.child_sg
            ));
            return None;
        }
    };
    Some((embedding, reciprocal, rotations))
}

/// `factor . vector`, componentwise.  `line_domain` keeps its own scaling helper
/// private, so the census scales the reference folded direction itself; the
/// result is compared against the engine's folded points with `same_mod`, not
/// with equality, so a componentwise rational product is all that is needed.
fn scale_point(vector: &Vec3R, factor: &Rat) -> Result<Vec3R, String> {
    let mut values = [Rat::ZERO; 3];
    for (axis, value) in vector.as_array().iter().enumerate() {
        values[axis] = value
            .checked_mul(*factor)
            .map_err(|error| format!("scaling {vector} by {factor}: {error}"))?;
    }
    Ok(Vec3R::new(values))
}

/// Everything one probe's [`BlockStat`]s need beyond the engine result itself.
struct BlockContext<'a> {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    label: &'static str,
    child_reciprocal: &'a Lattice,
    /// `q1 = T^T v`: the reference folded **direction**; the reference folded
    /// point of a probe is `t . q1`.
    reference_direction: Vec3R,
    /// The partition's reference arm, when the card-3 partition was built.  The
    /// `same_mod` reading of `carries_reference` is compared with it: the arm
    /// identity and the point identity must agree.
    reference_arm: Option<usize>,
    /// The partition's Gamma enumeration, when it was built.  Used to
    /// cross-check the blocks' own `carries_gamma` against the partition's
    /// arithmetic.
    gamma: Option<GammaArms<'a>>,
    cache: &'a PointClassCache,
}

/// One [`BlockStat`] per block of one answered probe, plus any consistency
/// failure found while reading them.
///
/// The block is described by its **own** data only: its points give the star
/// size, the arm count and the arm set; its folded coordinate gives the little
/// co-group and the cocycle class through the point-based entry points, and the
/// same two quantities are recomputed at every other point of the star -- star
/// points are conjugate, so a disagreement means the block's single reported
/// value would describe only one of its points.
fn block_stats(
    context: &BlockContext,
    parameter: Rat,
    result: &LineSubduction,
    dimensions: &BlockDimensions,
) -> (Vec<BlockStat>, Vec<String>) {
    let mut stats = Vec::with_capacity(result.blocks().len());
    let mut failures = Vec::new();
    let reference_point = match scale_point(&context.reference_direction, &parameter) {
        Ok(point) => point,
        Err(error) => {
            failures.push(format!(
                "ordinal {} {} t={parameter}: the reference folded point: {error}",
                context.ordinal, context.label
            ));
            context.reference_direction
        }
    };
    for (index, block) in result.blocks().iter().enumerate() {
        let stat = block_stat(context, &reference_point, parameter, index, block, &mut failures);
        stats.push(stat);
    }
    // The reference point must be carried by exactly one block, and the arms of
    // the blocks must partition the arm list the engine's dimension identity
    // gives.  Both are checked here as well as in `check_invariants`, where the
    // stored statistics are asserted instead of the slice they came from.
    let reference_blocks = stats.iter().filter(|stat| stat.carries_reference).count();
    if reference_blocks > 1 {
        failures.push(format!(
            "ordinal {} {} t={parameter}: {reference_blocks} blocks carry the reference folded \
             point",
            context.ordinal,
            context.label
        ));
    }
    let mut arm_seen = vec![0usize; dimensions.arm_total];
    for stat in &stats {
        for arm in &stat.arm_indices {
            match arm_seen.get_mut(*arm) {
                Some(slot) => *slot += 1,
                None => failures.push(format!(
                    "ordinal {} {} t={parameter}: block {} carries arm {arm} outside the \
                     {}-arm list",
                    context.ordinal,
                    context.label,
                    stat.index,
                    dimensions.arm_total
                )),
            }
        }
    }
    for (arm, count) in arm_seen.iter().enumerate() {
        if *count != 1 {
            failures.push(format!(
                "ordinal {} {} t={parameter}: arm {arm} appears in {count} blocks, not exactly \
                 one",
                context.ordinal,
                context.label
            ));
        }
    }
    (stats, failures)
}

/// One block's statistics; `failures` collects the pointwise disagreements.
fn block_stat(
    context: &BlockContext,
    reference_point: &Vec3R,
    parameter: Rat,
    index: usize,
    block: &FullStarBlock,
    failures: &mut Vec<String>,
) -> BlockStat {
    let arms = block.arm_indices();
    let star_size = block.points().len();
    let arm_count: usize = block.points().iter().map(|point| point.arm_count()).sum();
    let mut little_co_group = 0usize;
    let mut cocycle_trivial = true;
    let mut point_checks = 0usize;
    let mut point_disagreements = 0usize;
    let mut carries_gamma = false;
    let mut carries_reference = false;
    for point in block.points() {
        match context.child_reciprocal.contains(point.q()) {
            Ok(true) => carries_gamma = true,
            Ok(false) => {}
            Err(error) => failures.push(format!(
                "ordinal {} {} t={parameter} block {index}: the Gamma containment test failed: \
                 {error}",
                context.ordinal, context.label
            )),
        }
        match context
            .child_reciprocal
            .same_mod(point.q(), reference_point)
        {
            Ok(true) => carries_reference = true,
            Ok(false) => {}
            Err(error) => failures.push(format!(
                "ordinal {} {} t={parameter} block {index}: the reference-point test failed: \
                 {error}",
                context.ordinal, context.label
            )),
        }
        match context.cache.classify(context.child_sg, point.q()) {
            Ok(class) => {
                if point_checks == 0 {
                    little_co_group = class.little_co_group;
                    cocycle_trivial = class.cocycle_trivial;
                } else if class.little_co_group != little_co_group
                    || class.cocycle_trivial != cocycle_trivial
                {
                    point_disagreements += 1;
                    failures.push(format!(
                        "ordinal {} {} t={parameter} block {index}: the conjugate points of one \
                         child star disagree (order {} / {:?} against {} / {:?})",
                        context.ordinal,
                        context.label,
                        class.little_co_group,
                        class.cocycle_trivial,
                        little_co_group,
                        cocycle_trivial
                    ));
                }
                point_checks += 1;
            }
            Err(error) => failures.push(format!(
                "ordinal {} {} t={parameter} block {index}: point class: {error}",
                context.ordinal, context.label
            )),
        }
    }
    if let (Some(reference_arm), true) = (context.reference_arm, carries_reference)
        && !arms.contains(&reference_arm)
    {
        failures.push(format!(
            "ordinal {} {} t={parameter} block {index}: the block carries the reference folded \
             point but not the reference arm {reference_arm}",
            context.ordinal, context.label
        ));
    }
    if let Some((entries, zero_arms)) = context.gamma {
        let gamma = gamma_arms(entries, zero_arms, &parameter);
        let expected = arms.iter().any(|arm| gamma.contains(arm));
        if expected != carries_gamma {
            failures.push(format!(
                "ordinal {} {} t={parameter} block {index}: the block's folded points say \
                 carries_gamma={carries_gamma} but the partition's Gamma arithmetic says \
                 {expected} for arms {arms:?}",
                context.ordinal, context.label
            ));
        }
    }
    let stored_targets = block
        .targets()
        .iter()
        .filter(|target| target.irnumber.is_some())
        .count();
    BlockStat {
        ordinal: context.ordinal,
        parent_sg: context.parent_sg,
        child_sg: context.child_sg,
        label: context.label,
        parameter,
        index,
        star_size,
        arm_count,
        arm_indices: arms,
        block_dimension: block.block_dimension(),
        little_co_group,
        cocycle_trivial,
        source: classify_block(block.targets()),
        terms: block
            .targets()
            .iter()
            .map(|target| (target.dimension, target.multiplicity))
            .collect(),
        carries_reference,
        carries_gamma,
        stored_targets,
        target_count: block.targets().len(),
        declared_star_size: block.star_size(),
        declared_arm_count: block.arm_count(),
        point_checks,
        point_disagreements,
        representative_point: *block.q(),
    }
}

/// The arms at Gamma at one parameter: the partition's entries for it plus the
/// arms whose folded direction is exactly zero (at Gamma at every parameter).
fn gamma_arms(
    entries: &[(Rat, Vec<usize>)],
    zero_arms: &[usize],
    parameter: &Rat,
) -> Vec<usize> {
    let mut arms = zero_arms.to_vec();
    if let Some((_, listed)) = entries.iter().find(|(value, _)| value == parameter) {
        arms.extend(listed.iter().copied());
    }
    arms.sort_unstable();
    arms.dedup();
    arms
}

/// The engine's own dimension bookkeeping of one answered probe.
///
/// `covered` is the engine's own sum of block dimensions and `parent` its
/// `dimension x arms` identity; dividing the latter by the frozen little
/// dimension gives the number of parent arms the probe carried, which the block
/// arm sets must partition.  An inexact division is reported, never rounded.
fn block_dimensions(
    table: &LittleCharacterTable,
    result: &LineSubduction,
    failures: &mut Vec<String>,
) -> BlockDimensions {
    let parent = result.parent_dimension();
    let dimension = u32::from(table.dimension);
    let arm_total = if dimension == 0 || !parent.is_multiple_of(dimension) {
        failures.push(format!(
            "ordinal {} {}: the parent dimension {parent} is not a multiple of the frozen \
             little dimension {dimension}",
            result.ordinal(),
            result.label()
        ));
        0
    } else {
        // The conversion cannot fail on a target with 32-bit or wider `usize`;
        // it is a counted failure rather than a `unwrap_or(0)` so a narrow target
        // (or a future wider `u32`) shows up as the exact division it is instead
        // of as "no arm count" (card-6 sweep of the error-swallowing sites).
        match usize::try_from(parent / dimension) {
            Ok(arm_total) => arm_total,
            Err(error) => {
                failures.push(format!(
                    "ordinal {} {}: the arm count {} does not fit a usize: {error}",
                    result.ordinal(),
                    result.label(),
                    parent / dimension
                ));
                0
            }
        }
    };
    BlockDimensions {
        reported_blocks: result.blocks().len(),
        covered: result.covered_dimension(),
        parent,
        arm_total,
    }
}

/// Whether the whole child point group fixes one folded direction **exactly**
/// (`w_R = 0` for every child rotation).
///
/// This is the premise of the M2 local theorem and the only documented reason a
/// Gamma-reaching parameter can be interior to a partition interval: where every
/// child rotation fixes the arm's direction, the arm's little co-group never
/// grows, so no boundary of the partition is created by the arm arriving at
/// Gamma -- the condition changes *where the answer is read from*, not the
/// geometry.
fn exactly_fixed_arm(arm: &FoldedArm, rotations: &[Mat3I]) -> Result<bool, String> {
    for rotation in rotations {
        let image = Mat3R::from_ints(*rotation)
            .inverse()
            .and_then(|matrix| matrix.transpose().checked_mul_vector(&arm.direction))
            .map_err(|error| error.to_string())?;
        if image != arm.direction {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Fill the Gamma-reaching report of one `(record, label)` from the card-3
/// partition, returning the `(parameter, arms)` entries and the always-Gamma
/// arms for the block-side cross-check.
///
/// The condition reported here is a **provenance** one (`t v_i in L*_H`), not a
/// geometric one: where an arm reaches Gamma the engine may read stored child
/// Gamma data instead of building a target.  It is therefore *not* part of
/// [`FullStarPartition::boundaries`], and the report asserts the containment
/// that does hold -- every Gamma parameter is a boundary unless the reaching arm
/// is fixed exactly by the whole child point group, which is verified pointwise
/// before the exception is granted.
fn report_gamma(
    partition: &FullStarPartition,
    ordinal: usize,
    report: &mut GammaReport,
    failures: &mut Vec<String>,
) -> GammaParameters {
    let (entries, zero_arms) = match partition.gamma_parameters() {
        Ok(entries) => entries,
        Err(error) => {
            failures.push(format!(
                "ordinal {} {}: the Gamma enumeration failed: {error}",
                ordinal, partition.label
            ));
            return (Vec::new(), Vec::new());
        }
    };
    report.entries += entries.len();
    report.zero_arms += zero_arms.len();
    let shape: Vec<(String, Vec<usize>)> = entries
        .iter()
        .map(|(parameter, arms)| (parameter.to_string(), arms.clone()))
        .collect();
    *report.shapes.entry(shape).or_insert(0) += 1;
    for (parameter, arms) in &entries {
        if partition.is_boundary(parameter) {
            report.contained += 1;
            continue;
        }
        // Not a boundary: the documented exception must hold, and it is verified
        // rather than asserted -- for an arm listed at this parameter, every
        // child rotation must fix the arm's direction exactly.
        let mut explained = None;
        for arm in arms {
            let Some(folded) = partition.arms.get(*arm) else {
                failures.push(format!(
                    "ordinal {} {} t={parameter}: the Gamma arm {arm} is not in the arm list",
                    ordinal, partition.label
                ));
                continue;
            };
            match exactly_fixed_arm(folded, &partition.child_rotations) {
                Ok(true) => {
                    explained = Some((*arm, folded.parent_rotation));
                    break;
                }
                Ok(false) => {}
                Err(error) => failures.push(format!(
                    "ordinal {} {} t={parameter}: the exact-fixity test of arm {arm} failed: \
                     {error}",
                    ordinal, partition.label
                )),
            }
        }
        match explained {
            Some((arm, _)) => {
                report.exceptional += 1;
                if report.witnesses.len() < GAMMA_WITNESSES {
                    report.witnesses.push(format!(
                        "ordinal {} {} t={parameter}: arm {arm} (direction ({}, {}, {})) reaches \
                         Gamma inside an interval and is fixed exactly by all {} child rotations",
                        ordinal,
                        partition.label,
                        partition.arms[arm].direction.get(0),
                        partition.arms[arm].direction.get(1),
                        partition.arms[arm].direction.get(2),
                        partition.child_rotations.len()
                    ));
                }
            }
            None => failures.push(format!(
                "ordinal {} {} t={parameter}: arms {arms:?} reach Gamma but the parameter is not \
                 a full-star boundary and no reaching arm is fixed exactly by the whole child \
                 point group",
                ordinal, partition.label
            )),
        }
    }
    (entries, zero_arms)
}

/// The canonical geometry of a block list: the parent arm indices grouped by
/// child star, ascending.
///
/// Block **order**, folded coordinates and provenance all move with the
/// parameter; which arms share a child star does not (inside one partition
/// interval).  Card 5 matches blocks across parameters by exactly this key.
fn block_geometry(blocks: &[BlockStat]) -> Vec<Vec<usize>> {
    let mut geometry: Vec<Vec<usize>> = blocks
        .iter()
        .map(|block| block.arm_indices.clone())
        .collect();
    geometry.sort();
    geometry
}

// ── R6.7 card 5: the two interior points of one interval ─────────────────────

/// Parse one written rational back (`n` or `n/d`).
fn parse_parameter(text: &str) -> Result<Rat, String> {
    let (numerator, denominator) = match text.split_once('/') {
        Some((numerator, denominator)) => (numerator, denominator),
        None => (text, "1"),
    };
    let numerator = numerator
        .trim()
        .parse::<i128>()
        .map_err(|error| format!("{text:?}: numerator: {error}"))?;
    let denominator = denominator
        .trim()
        .parse::<i128>()
        .map_err(|error| format!("{text:?}: denominator: {error}"))?;
    Rat::new(numerator, denominator).map_err(|error| format!("{text:?}: {error}"))
}

/// Stable sort key of one exact rotation, so the two parameters' little groups
/// can be aligned by rotation without relying on either list's order.
fn rotation_key(rotation: &Mat3I) -> [i32; 9] {
    let mut key = [0i32; 9];
    for (slot, value) in key.iter_mut().enumerate() {
        *value = rotation[slot / 3][slot % 3];
    }
    key
}

/// The open intervals of one `(record, label)`'s probe set, wrapping around
/// `t = 1` like [`FullStarPartition::intervals`] but built from
/// `probe_parameters(parent)`: the parent's formal boundaries are a separate set
/// from the full-star boundaries, and on this corpus they happen to be nested
/// inside the eighth grid, so `intervals()` alone cannot exercise their union.
fn probe_intervals(parameters: &[Rat]) -> Result<Vec<(Rat, Rat)>, String> {
    let Some(first) = parameters.first().copied() else {
        return Ok(vec![(Rat::ZERO, Rat::ONE)]);
    };
    let mut out = Vec::with_capacity(parameters.len());
    for (index, start) in parameters.iter().enumerate() {
        let end = match parameters.get(index + 1) {
            Some(next) => *next,
            None => first.checked_add(Rat::ONE).map_err(|error| error.to_string())?,
        };
        out.push((*start, end));
    }
    Ok(out)
}

/// The two exact interior points of one interval: its trisection points.
///
/// A non-positive interval is a typed error rather than two copies of the same
/// boundary point: `interior_points` is what makes two *different* parameters of
/// the interval comparable, and a zero-length interval would compare a
/// decomposition with itself while still counting as a checked interval.
fn interior_points(left: &Rat, right: &Rat) -> Result<(Rat, Rat, Rat), String> {
    let span = right.checked_sub(*left).map_err(|error| error.to_string())?;
    // `Rat` is normalized with a positive denominator, so the sign is the
    // numerator's.
    if span.numerator() <= 0 {
        return Err(format!("the interval ({left}, {right}) is not positive"));
    }
    let fraction = Rat::new(INTERIOR_FRACTION.0, INTERIOR_FRACTION.1)
        .map_err(|error| error.to_string())?;
    let step = span.checked_mul(fraction).map_err(|error| error.to_string())?;
    let first = left.checked_add(step).map_err(|error| error.to_string())?;
    let second = first.checked_add(step).map_err(|error| error.to_string())?;
    Ok((step, first, second))
}

// ── R6.7 card 6: the independent arm-geometry audit ─────────────────────────

/// One full-star relation, recomputed by the census from the folded arms and the
/// child reciprocal lattice **itself**.
///
/// The predicate is the card-3 one, `t (R^-T v_arm - v_image) in L*_H`, and the
/// triples are exactly the partition's: every `(arm, image, rotation)` with the
/// reflexive `(i, i, identity)` excluded.  The **decision** is not the
/// partition's: the difference vector is formed explicitly, mapped once into the
/// child reciprocal-lattice basis, and evaluated at a parameter by scaling that
/// coordinate vector and asking for integrality.  The partition's own solver
/// (`minimal_parameter_step`) and its residue enumeration are never called here,
/// so a partition that reports a wrong parameter set, a wrong event count or a
/// boundary that no relation supports cannot agree with this audit by
/// construction.
struct CensusRelation {
    /// Which part of the folded geometry the relation can change.
    kind: StarEventKind,
    /// `C_H (R^-T v_arm - v_image)`: `t w in L*_H` is exactly `t * coordinates in
    /// Z^3`, because a vector is `rows^T . coordinates` in the lattice basis.
    coordinates: Vec3R,
    /// Whether the difference vector is exactly zero: then the relation holds at
    /// **every** parameter, which is why the partition records it as permanent
    /// instead of cutting the domain at it.
    permanent: bool,
}

impl CensusRelation {
    /// Whether `t w in L*_H` at this parameter.
    fn holds(&self, parameter: &Rat) -> Result<bool, String> {
        if self.permanent {
            return Ok(true);
        }
        for axis in 0..3 {
            let component = self.coordinates.get(axis);
            if component.is_zero() {
                continue;
            }
            let scaled = parameter.checked_mul(component).map_err(|error| {
                format!("the relation coordinate at axis {axis} times t={parameter}: {error}")
            })?;
            if scaled.denominator() != 1 {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// What the card-6 independent arm-geometry audit measured over the corpus.
#[derive(Default)]
struct GeometryReport {
    /// `(record, label)` pairs whose relation list was rebuilt.
    pairs: usize,
    /// Relations enumerated (the partition's triple set, minus the reflexive
    /// `(i, i, identity)` relation).
    relations: usize,
    /// The closed form the enumeration has to match, `n^2 r - n`, accumulated
    /// from the partition's own arm and rotation counts.
    predicted_relations: usize,
    /// Relations whose difference vector is exactly zero.
    permanent: usize,
    /// `Σ partition.permanent_counts`: the partition's own account of the same
    /// relation set.
    permanent_partition: usize,
    /// Boundary parameters of the partition at which the whole relation list was
    /// evaluated and compared with the recorded event counts.
    boundaries: usize,
    /// Boundary parameters the independent residue computation derives.
    recomputed_boundaries: usize,
    /// Predicate instances evaluated at boundary parameters.
    evaluated: usize,
    /// Interior points of the sweep's intervals at which the relation list was
    /// evaluated (two per interval).
    interior_points: usize,
    /// Predicate instances evaluated at those interior points.
    interior_evaluated: usize,
    /// Interior points at which a non-permanent relation held: a cut the
    /// partition does not report, so a violation rather than a note.
    interior_holds: usize,
    /// Disagreements between the audit and the partition, of any kind.
    disagreements: usize,
    /// Orbit–stabiliser arm-count checks, and how many of them disagreed with
    /// `partition.arms.len()`.
    arm_orbit_checks: usize,
    arm_orbit_disagreements: usize,
    /// `Σ partition.arms.len()`, the arm count the orbit–stabiliser route
    /// predicts.
    arms: usize,
    /// `Σ probe parameter count x frozen little dimension x arm count`: the
    /// closed form of the boundary pass' reported parent dimension, accumulated
    /// from the partition's own counts.  It is the independent prediction of
    /// [`SweepReport::boundary_parent_dimension`], which comes from the returned
    /// decompositions themselves.
    predicted_boundary_dimension: u64,
}

/// The census's own enumeration of one pair's full-star relations.
///
/// The rotation action is applied explicitly (`R^-T v = (R^T)^-1 v` on the arm's
/// folded direction) and the difference is expressed in the child reciprocal
/// lattice basis; nothing here reads the partition's events, counts or steps.
fn census_relations(
    partition: &FullStarPartition,
    key: &str,
) -> Result<Vec<CensusRelation>, String> {
    let arms = &partition.arms;
    let rotations = &partition.child_rotations;
    let mut out = Vec::with_capacity(arms.len() * arms.len() * rotations.len());
    for (arm_index, arm) in arms.iter().enumerate() {
        for rotation in rotations {
            let action = Mat3R::from_ints(*rotation)
                .inverse()
                .map_err(|error| format!("{key}: the rotation {rotation:?} is singular: {error}"))?
                .transpose();
            let image = action
                .checked_mul_vector(&arm.direction)
                .map_err(|error| format!("{key}: the image of arm {arm_index}: {error}"))?;
            for (image_index, image_arm) in arms.iter().enumerate() {
                if arm_index == image_index && *rotation == IDENTITY_ROTATION {
                    // The reflexive relation: an arm is itself at every
                    // parameter.  The partition excludes it for the same reason.
                    continue;
                }
                let difference = image.checked_sub(&image_arm.direction).map_err(|error| {
                    format!("{key}: the difference of arms {arm_index} and {image_index}: {error}")
                })?;
                let coordinates = partition
                    .child_reciprocal
                    .coordinates(&difference)
                    .map_err(|error| {
                        format!("{key}: the child lattice coordinates of the difference: {error}")
                    })?;
                let kind = if arm_index == image_index {
                    StarEventKind::LittleCoGroupGrowth
                } else if *rotation == IDENTITY_ROTATION {
                    StarEventKind::ArmMerge
                } else {
                    StarEventKind::OrbitIdentification
                };
                out.push(CensusRelation {
                    kind,
                    coordinates,
                    permanent: difference.is_zero(),
                });
            }
        }
    }
    Ok(out)
}

/// Evaluate the whole relation list at one parameter, counting by kind.
///
/// **Only the non-permanent relations are counted**: an identically-true
/// relation holds at every parameter by construction, so counting it here would
/// make every interior point look like a missed cut.  Those relations are
/// accounted for separately, by comparing their number with the partition's own
/// `permanent_counts`.
///
/// `evaluated` is incremented once per relation **before** the predicate is
/// decided, so a relation that fails to be decided is still counted as
/// evaluated.
fn relation_counts(
    relations: &[CensusRelation],
    parameter: &Rat,
    evaluated: &mut usize,
) -> Result<[usize; 3], String> {
    let mut counts = [0usize; 3];
    for relation in relations {
        if relation.permanent {
            continue;
        }
        *evaluated += 1;
        if relation.holds(parameter)? {
            counts[relation.kind.index()] += 1;
        }
    }
    Ok(counts)
}

/// `t mod 1` of a normalized rational, in `[0, 1)`.
fn reduce_modulo_one(value: Rat) -> Rat {
    let denominator = value.denominator();
    Rat::new(value.numerator().rem_euclid(denominator), denominator)
        .expect("a residue of a normalized rational is a rational")
}

/// Sorted, deduplicated exact parameters: the census-side set representation,
/// because `Rat` deliberately has no `Ord` and a `BTreeSet` is therefore not
/// available.
fn parameter_set(values: Vec<Rat>) -> Vec<Rat> {
    let mut sorted = sorted_parameters(values);
    sorted.dedup();
    sorted
}

/// The values of the sorted `left` that are not in the sorted `right`.
fn parameter_difference(left: &[Rat], right: &[Rat]) -> Vec<Rat> {
    left.iter()
        .copied()
        .filter(|value| !right.contains(value))
        .collect()
}

/// A parameter list for a failure message, through `Rat`'s `Display`: its
/// `Debug` prints the private fields, which would hide the parameter's value.
fn parameters_text(values: &[Rat]) -> String {
    values
        .iter()
        .map(Rat::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}

/// Least positive rational that is an integer multiple of both arguments.
///
/// For normalized positive `p1/q1` and `p2/q2` the intersection of the subgroups
/// they generate is generated by `lcm(p1, p2) / gcd(q1, q2)`.
fn rational_lcm(left: &Rat, right: &Rat) -> Result<Rat, String> {
    let left_numerator = left
        .numerator()
        .checked_abs()
        .ok_or_else(|| "the relation step overflows".to_string())?;
    let right_numerator = right
        .numerator()
        .checked_abs()
        .ok_or_else(|| "the relation step overflows".to_string())?;
    let numerator = if left_numerator == 0 || right_numerator == 0 {
        0
    } else {
        (left_numerator / gcd(left_numerator, right_numerator))
            .checked_mul(right_numerator)
            .ok_or_else(|| "the relation step overflows".to_string())?
    };
    Rat::new(numerator, gcd(left.denominator(), right.denominator()))
        .map_err(|error| format!("the relation step: {error}"))
}

/// The parameters of `[0, 1)` at which one relation holds, from its coordinate
/// vector alone.
///
/// `t w in L*` holds iff `t c_i in Z` for every component `c_i = p_i / q_i` of
/// the coordinate vector, i.e. iff `t` is a multiple of `q_i / |p_i|` for every
/// nonzero component, so the solution set is the subgroup generated by the
/// rational least common multiple of those generators and its residues in
/// `[0, 1)` are the multiples of the reduced generator.
fn relation_residues(relation: &CensusRelation) -> Result<Vec<Rat>, String> {
    if relation.permanent {
        return Ok(Vec::new());
    }
    let mut step: Option<Rat> = None;
    for axis in 0..3 {
        let component = relation.coordinates.get(axis);
        if component.is_zero() {
            continue;
        }
        let numerator = component
            .numerator()
            .checked_abs()
            .ok_or_else(|| "a relation coordinate overflows".to_string())?;
        let generator = Rat::new(component.denominator(), numerator)
            .map_err(|error| format!("the relation generator: {error}"))?;
        step = Some(match step {
            None => generator,
            Some(current) => rational_lcm(&current, &generator)?,
        });
    }
    let Some(step) = step else {
        // Every component is zero, which is the `permanent` case above; kept as a
        // defensive guard rather than a silent empty answer.
        return Err("the relation has no nonzero coordinate but is not permanent".to_string());
    };
    let capacity = usize::try_from(step.denominator())
        .map_err(|error| format!("the relation residue count: {error}"))?;
    let mut out = Vec::with_capacity(capacity);
    for multiple in 0..step.denominator() {
        let factor = Rat::new(multiple, 1)
            .map_err(|error| format!("the relation residue {multiple}: {error}"))?;
        let value = step
            .checked_mul(factor)
            .map_err(|error| format!("the relation residue {multiple}: {error}"))?;
        out.push(reduce_modulo_one(value));
    }
    Ok(out)
}

/// The arm count as an **orbit–stabiliser** count.
///
/// The arms are the distinct images `R^-T v` of the frozen direction under the
/// parent's rotation group, so `|orbit| = |G| / |Stab_G(v)|`.  The group is read
/// here through [`rotation_set`] -- a different route from the deduplicated Hall
/// operations `arm_images` enumerates -- and the stabiliser is counted by exact
/// vector equality of the images.
fn arm_orbit_count(parent_sg: u8, direction: &Vec3R) -> Result<usize, String> {
    let rotations =
        rotation_set(parent_sg).map_err(|error| format!("SG {parent_sg}: {error}"))?;
    let mut stabiliser = 0usize;
    for rotation in &rotations {
        let image = Mat3R::from_ints(*rotation)
            .inverse()
            .map_err(|error| format!("SG {parent_sg}: the rotation {rotation:?} is singular: {error}"))?
            .transpose()
            .checked_mul_vector(direction)
            .map_err(|error| format!("SG {parent_sg}: the rotation image: {error}"))?;
        if image == *direction {
            stabiliser += 1;
        }
    }
    if stabiliser == 0 || !rotations.len().is_multiple_of(stabiliser) {
        return Err(format!(
            "the orbit–stabiliser count of SG {parent_sg} is not exact ({} rotation(s), \
             {stabiliser} fixing the direction)",
            rotations.len()
        ));
    }
    Ok(rotations.len() / stabiliser)
}

/// Recompute one pair's arm geometry from the folded arms and the child
/// reciprocal lattice, and check the partition against it.
///
/// Returns the relation list for the interior audit of the sweep's intervals; it
/// is empty when the enumeration itself failed, which is already a counted
/// failure.
fn audit_arm_geometry(
    table: &'static LittleCharacterTable,
    partition: &FullStarPartition,
    probe_parameters: usize,
    key: &str,
    report: &mut SweepReport,
) -> Vec<CensusRelation> {
    let relations = match census_relations(partition, key) {
        Ok(relations) => relations,
        Err(error) => {
            report.geometry.disagreements += 1;
            report
                .failures
                .push(format!("{key}: the independent relation enumeration failed: {error}"));
            return Vec::new();
        }
    };
    report.geometry.pairs += 1;
    // The closed form `n^2 r - n` (the reflexive `(i, i, identity)` triples are
    // excluded) is an independent count of the same enumeration, read from the
    // partition's own arm and rotation counts: a truncated or duplicated
    // enumeration cannot match it without a matching change in the partition.
    let predicted = partition
        .arms
        .len()
        .saturating_mul(partition.arms.len())
        .saturating_mul(partition.child_rotations.len())
        .saturating_sub(partition.arms.len());
    report.geometry.relations += relations.len();
    report.geometry.predicted_relations += predicted;
    if relations.len() != predicted {
        report.geometry.disagreements += 1;
        report.failures.push(format!(
            "{key}: the independent relation enumeration produced {} relation(s) but the closed \
             form n^2 r - n of the partition's own counts predicts {predicted}",
            relations.len()
        ));
    }
    let permanent = relations
        .iter()
        .filter(|relation| relation.permanent)
        .count();
    let partition_permanent: usize = partition.permanent_counts.iter().sum();
    report.geometry.permanent += permanent;
    report.geometry.permanent_partition += partition_permanent;
    if permanent != partition_permanent {
        report.geometry.disagreements += 1;
        report.failures.push(format!(
            "{key}: the audit finds {permanent} identically-true relation(s) but the partition \
             records {partition_permanent}"
        ));
    }
    // The boundary set, recomputed as the union of the relations' residue sets.
    // `Rat` deliberately has no `Ord`, so the set is a sorted, deduplicated
    // vector rather than a `BTreeSet`.
    let mut recomputed: Vec<Rat> = Vec::new();
    for relation in &relations {
        match relation_residues(relation) {
            Ok(residues) => recomputed.extend(residues),
            Err(error) => {
                report.geometry.disagreements += 1;
                report
                    .failures
                    .push(format!("{key}: the independent residue computation: {error}"));
            }
        }
    }
    let recomputed = parameter_set(recomputed);
    let claimed_raw = partition.boundary_parameters();
    let claimed = parameter_set(claimed_raw.clone());
    report.geometry.recomputed_boundaries += recomputed.len();
    if claimed.len() != claimed_raw.len() {
        report.geometry.disagreements += 1;
        report.failures.push(format!(
            "{key}: the partition lists {} boundary parameter(s) but only {} distinct one(s)",
            claimed_raw.len(),
            claimed.len()
        ));
    }
    if recomputed != claimed {
        report.geometry.disagreements += 1;
        report.failures.push(format!(
            "{key}: the independently recomputed full-star boundary set is not the partition's: \
             the partition cuts [{}] that no relation supports, and does not cut [{}] where a \
             relation holds",
            parameters_text(&parameter_difference(&claimed, &recomputed)),
            parameters_text(&parameter_difference(&recomputed, &claimed))
        ));
    }
    // At every boundary the partition claims, the relations that hold must be
    // exactly the ones it records: the counts are exact on both sides.
    for parameter in &claimed {
        report.geometry.boundaries += 1;
        let mut evaluated = report.geometry.evaluated;
        let counts = match relation_counts(&relations, parameter, &mut evaluated) {
            Ok(counts) => counts,
            Err(error) => {
                report.geometry.evaluated = evaluated;
                report.geometry.disagreements += 1;
                report.failures.push(format!(
                    "{key} t={parameter}: the boundary relation audit: {error}"
                ));
                continue;
            }
        };
        report.geometry.evaluated = evaluated;
        let recorded = partition
            .boundary(parameter)
            .map_or([0usize; 3], |boundary| boundary.counts);
        if counts != recorded {
            report.geometry.disagreements += 1;
            report.failures.push(format!(
                "{key} t={parameter}: the audit finds {}/{}/{} relation(s) holding (little-co-group \
                 growth / arm merge / orbit identification) but the partition records {}/{}/{}",
                counts[0], counts[1], counts[2], recorded[0], recorded[1], recorded[2]
            ));
        }
    }
    // The arm count by orbit and stabiliser, against the partition's own list.
    // The closed form of the boundary pass' reported parent dimension: every
    // probe parameter contributes `frozen little dimension x arms`, and the
    // engine's own `parent_dimension()` is what the pass accumulated.
    report.geometry.predicted_boundary_dimension += (probe_parameters as u64)
        * u64::from(table.dimension)
        * (partition.arms.len() as u64);
    report.geometry.arm_orbit_checks += 1;
    report.geometry.arms += partition.arms.len();
    match line_direction(table) {
        Some(direction) => match arm_orbit_count(table.space_group, &direction) {
            Ok(orbit) if orbit == partition.arms.len() => {}
            Ok(orbit) => {
                report.geometry.arm_orbit_disagreements += 1;
                report.failures.push(format!(
                    "{key}: the orbit–stabiliser count over the parent's rotations gives {orbit} \
                     arm(s) but the partition lists {}",
                    partition.arms.len()
                ));
            }
            Err(error) => {
                report.geometry.arm_orbit_disagreements += 1;
                report.failures.push(format!("{key}: the orbit–stabiliser arm count: {error}"));
            }
        },
        None => {
            report.geometry.arm_orbit_disagreements += 1;
            report.failures.push(format!(
                "{key}: the frozen direction does not parse, so the arm count cannot be recomputed"
            ));
        }
    }
    relations
}

/// The child little group at one exact child point, in the child **Hall** frame
/// [`FullStarBlock::target_little_character`] evaluates in.
///
/// The whole `embedding.representatives()` list is unfolded with the engine's
/// own `unmap_operation` and the embedding's `child_shift`, then filtered by the
/// exact fixity test `R^-T q = q (mod L*_child)` — the same filter, frame and
/// shift `decompose::little_group_operations` applies internally.  The result is
/// one operation per rotation [measured: the group's size equals
/// `point_little_co_group_order` at every seed point of the sweep], and it is
/// the **same** list at both parameters because the filter only depends on the
/// rotation set, which the interval keeps fixed.
fn child_little_group(
    embedding: &SubgroupEmbedding,
    child_reciprocal: &Lattice,
    point: &Vec3R,
) -> Result<Vec<ExactSeitz>, String> {
    let shift = embedding
        .child_shift()
        .checked_neg()
        .map_err(|error| error.to_string())?;
    let mut out = Vec::new();
    for operation in embedding.representatives() {
        let child = embedding
            .transform()
            .unmap_operation(operation)
            .map_err(|error| error.to_string())?;
        match child_reciprocal.preserves(child.rotation(), point) {
            Ok(true) => {}
            Ok(false) => continue,
            Err(error) => return Err(error.to_string()),
        }
        let rotated = Mat3R::from_ints(child.rotation())
            .checked_mul_vector(&shift)
            .map_err(|error| error.to_string())?;
        let translation = child
            .translation()
            .checked_add(&shift)
            .and_then(|value| value.checked_sub(&rotated))
            .map_err(|error| error.to_string())?;
        out.push(ExactSeitz::new(child.rotation(), translation));
    }
    if out.is_empty() {
        return Err("the child little group of the point is empty".to_string());
    }
    Ok(out)
}

/// The child rotations fixing one exact child point modulo the child reciprocal
/// lattice: the little co-group as a **set**, not only its order.
fn little_co_group_rotations(
    child_reciprocal: &Lattice,
    rotations: &[Mat3I],
    point: &Vec3R,
) -> Result<Vec<Mat3I>, String> {
    let mut out = Vec::new();
    for rotation in rotations {
        match child_reciprocal.preserves(*rotation, point) {
            Ok(true) => out.push(*rotation),
            Ok(false) => {}
            Err(error) => return Err(error.to_string()),
        }
    }
    out.sort_by_key(rotation_key);
    out.dedup();
    Ok(out)
}

/// Layer 1 of the card-5 comparison: the parameter-stable geometry of one block.
#[derive(Debug, Clone)]
struct SweepBlock {
    /// `FullStarBlock::arm_indices()`: the block's identity across parameters.
    arm_indices: Vec<usize>,
    /// Per folded point, the parent arms folding onto it (each ascending); the
    /// multiset is compared as a whole because the star's point **order** moves
    /// with the parameter.
    points: Vec<Vec<usize>>,
    star_size: usize,
    arm_count: usize,
    block_dimension: u32,
    /// Little co-group of the block's own representative point.
    rotations: Vec<Mat3I>,
}

impl SweepBlock {
    /// Read one block's geometry from the engine's own block object, **bound to
    /// the parameter whose decomposition produced it**.
    ///
    /// Every folded point is checked against `t * direction` of each parent arm
    /// it carries: the point of a block is a function of the parameter, so a
    /// block read from one parameter's decomposition and compared as if it were
    /// the other parameter's is a lost check that every other guard of the sweep
    /// accepts (measured: taking both sides of the comparison from the first
    /// parameter's decomposition — `M5` in the card-5 report — left the gate at
    /// exit 0 with 46,048 "comparisons" and no failure).  With this binding the
    /// same mutation fails on the first interval, before any matching.
    fn read(
        context: &SweepContext,
        parameter: Rat,
        block: &FullStarBlock,
    ) -> Result<Self, String> {
        for point in block.points() {
            for arm in point.arm_indices() {
                let folded = context.arms.get(*arm).ok_or_else(|| {
                    format!(
                        "block {:?} carries the arm index {arm}, which the partition's arm table \
                         does not have",
                        block.arm_indices()
                    )
                })?;
                let expected = scale_point(&folded.direction, &parameter)
                    .map_err(|error| format!("the folded direction of arm {arm}: {error}"))?;
                // Modulo the child reciprocal lattice, the equivalence the star
                // constructor merges arms by: at a **boundary** two merged arms
                // fold onto one point and the engine reports a class
                // representative, so exact equality would be wrong there.
                let same = context
                    .child_reciprocal
                    .same_mod(&expected, point.q())
                    .map_err(|error| format!("the point of arm {arm}: {error}"))?;
                if !same {
                    return Err(format!(
                        "at t={parameter} the block with arms {:?} reports the point {} for arm \
                         {arm}, but that arm's folded direction at this parameter is {expected}, \
                         which is not the same point modulo the child reciprocal lattice",
                        block.arm_indices(),
                        point.q()
                    ));
                }
            }
        }
        let mut points: Vec<Vec<usize>> = block
            .points()
            .iter()
            .map(|point| {
                let mut arms = point.arm_indices().to_vec();
                arms.sort_unstable();
                arms
            })
            .collect();
        points.sort();
        Ok(Self {
            arm_indices: block.arm_indices(),
            arm_count: block.points().iter().map(|point| point.arm_count()).sum(),
            star_size: block.points().len(),
            block_dimension: block.block_dimension(),
            points,
            rotations: little_co_group_rotations(
                context.child_reciprocal,
                context.child_rotations,
                block.q(),
            )?,
        })
    }
}

/// Match the blocks of two parameters one-to-one by arm set and compare every
/// parameter-stable quantity, naming the offending block on any mismatch.
fn match_block_geometry(
    left: &[SweepBlock],
    right: &[SweepBlock],
    key: &str,
) -> Result<Vec<(usize, usize)>, String> {
    if left.len() != right.len() {
        return Err(format!(
            "{key}: {} block(s) at the first parameter against {} at the second",
            left.len(),
            right.len()
        ));
    }
    let mut matching = Vec::with_capacity(left.len());
    let mut taken = vec![false; right.len()];
    for (index, block) in left.iter().enumerate() {
        let candidates: Vec<usize> = right
            .iter()
            .enumerate()
            .filter(|(slot, other)| !taken[*slot] && other.arm_indices == block.arm_indices)
            .map(|(slot, _)| slot)
            .collect();
        let partner = match candidates.as_slice() {
            [only] => *only,
            [] => {
                return Err(format!(
                    "{key}: block {index} with arms {:?} has no partner at the second parameter",
                    block.arm_indices
                ));
            }
            many => {
                return Err(format!(
                    "{key}: block {index} with arms {:?} matches {many:?} blocks at the second \
                     parameter, so the identity is not one-to-one",
                    block.arm_indices
                ));
            }
        };
        taken[partner] = true;
        let other = &right[partner];
        for (name, left_value, right_value) in [
            ("star size", block.star_size, other.star_size),
            ("arm count", block.arm_count, other.arm_count),
            (
                "block dimension",
                block.block_dimension as usize,
                other.block_dimension as usize,
            ),
        ] {
            if left_value != right_value {
                return Err(format!(
                    "{key}: block {index} (arms {:?}) has {name} {left_value} at the first \
                     parameter and {right_value} at the second",
                    block.arm_indices
                ));
            }
        }
        if block.points != other.points {
            return Err(format!(
                "{key}: block {index} (arms {:?}) folds its parent arms onto different points \
                 ({:?} against {:?})",
                block.arm_indices, block.points, other.points
            ));
        }
        if block.rotations != other.rotations {
            return Err(format!(
                "{key}: block {index} (arms {:?}) has little co-group order {} at the first \
                 parameter and {} at the second",
                block.arm_indices,
                block.rotations.len(),
                other.rotations.len()
            ));
        }
        matching.push((index, partner));
    }
    Ok(matching)
}

/// Layer 2 of the card-5 comparison: one reported target with its
/// **gauge-unified** little-group character vector.
#[derive(Debug, Clone, PartialEq)]
struct SweepTarget {
    dimension: u8,
    multiplicity: u32,
    /// Provenance, reported but deliberately **not** part of the identity: a
    /// stored and a constructed presentation of the same target must match.
    stored: bool,
    /// `chi_target(h) * conj(exp(2 pi i q(tau_seed) . tau_h))` over the common
    /// operation list.
    vector: Vec<Complex64>,
}

/// Normalized **signed** inner product of two character vectors: `1` for equal
/// characters, `-1` for one equal to the other's negative, `0` for orthogonal
/// ones.
///
/// The absolute value this used to take only tested collinearity, which accepted
/// a global `-1` as "the same target" (external review of the card-5 line, P2:
/// multiplying every second-side vector by `-1` left the whole gate green).  Once
/// the card-1 gauge is divided out there is no residual phase freedom left — the
/// remaining object is an honest little-group character, and its value at the
/// identity is the **dimension** (a positive real), which
/// [`gauge_target_slice`] checks separately.  Comparing the real part is
/// therefore the right test, and it is what makes `[1, 1]` and `[-1, -1]`
/// different targets.
fn character_score(left: &[Complex64], right: &[Complex64]) -> f64 {
    let mut inner = Complex64::new(0.0, 0.0);
    let mut left_norm = 0.0;
    let mut right_norm = 0.0;
    for (first, second) in left.iter().zip(right) {
        inner += first * second.conj();
        left_norm += first.norm_sqr();
        right_norm += second.norm_sqr();
    }
    if left_norm == 0.0 || right_norm == 0.0 {
        return 0.0;
    }
    inner.re / (left_norm.sqrt() * right_norm.sqrt())
}

/// Match two blocks' target lists one-to-one by gauge-unified character, then
/// compare the matched multiplicities.
///
/// This is the card-5 core.  Matching is by the character vector rather than by
/// any sorted `(dimension, multiplicity)` list, so two distinct one-dimensional
/// targets whose multiplicities are swapped are caught: the swap leaves the
/// sorted list unchanged but attaches the wrong multiplicity to each character.
/// A target with no partner, a tie at the tolerance, a dimension disagreement
/// and a matched pair with different multiplicities are all explicit errors.
///
/// Provenance is **not** part of the identity: [`SweepTarget::stored`] is
/// carried for the report only, so a stored and a constructed presentation of
/// the same character match (the corpus has no interior stored target, so this
/// branch is covered by the module tests instead).
fn match_targets(
    left: &[SweepTarget],
    right: &[SweepTarget],
    key: &str,
) -> Result<(Vec<(usize, usize)>, f64), TargetMismatch> {
    if left.len() != right.len() {
        return Err(TargetMismatch::Ambiguous(format!(
            "{key}: {} target(s) at the first parameter against {} at the second",
            left.len(),
            right.len()
        )));
    }
    let mut matching = Vec::with_capacity(left.len());
    let mut taken = vec![false; right.len()];
    let mut worst = 1.0f64;
    for (index, target) in left.iter().enumerate() {
        // Every candidate at tolerance, without a preference order: a tie is an
        // error rather than a silent pick.
        let hits: Vec<(usize, f64)> = right
            .iter()
            .enumerate()
            .filter(|(slot, _)| !taken[*slot])
            .map(|(slot, other)| (slot, character_score(&target.vector, &other.vector)))
            .filter(|(_, score)| *score >= 1.0 - SWEEP_TOLERANCE)
            .collect();
        let partner = match hits.as_slice() {
            [(only, score)] => {
                worst = worst.min(*score);
                *only
            }
            [] => {
                return Err(TargetMismatch::Ambiguous(format!(
                    "{key}: target {index} (dimension {}, multiplicity {}, {}) has no character \
                     partner at the second parameter",
                    target.dimension,
                    target.multiplicity,
                    if target.stored { "stored" } else { "constructed" }
                )));
            }
            many => {
                return Err(TargetMismatch::Ambiguous(format!(
                    "{key}: target {index} (dimension {}, multiplicity {}) has {} character \
                     partners at the second parameter, so the identity is not one-to-one",
                    target.dimension,
                    target.multiplicity,
                    many.len()
                )));
            }
        };
        taken[partner] = true;
        let other = &right[partner];
        if target.dimension != other.dimension {
            return Err(TargetMismatch::Ambiguous(format!(
                "{key}: matched target {index} has dimension {} at the first parameter and {} at \
                 the second",
                target.dimension, other.dimension
            )));
        }
        if target.multiplicity != other.multiplicity {
            return Err(TargetMismatch::Multiplicity(format!(
                "{key}: the matched dimension-{} target {index} has multiplicity {} at the first \
                 parameter and {} at the second (the character match is the identity, so this is \
                 a real multiplicity change, not a relabelling)",
                target.dimension, target.multiplicity, other.multiplicity
            )));
        }
        matching.push((index, partner));
    }
    Ok((matching, worst))
}

/// The gauge-unified character vectors of one block's targets, at one of the
/// block's own folded points (`point`, the common seed arm's point), over the
/// child little group of that point.
///
/// The gauge is the card-1 explicit cochain: the reported little-group character
/// at the operation `h` carries the Bloch phase `exp(2 pi i q(t) . tau_h)`,
/// which moves with the parameter, so dividing by it is what makes two
/// parameters of one interval comparable.  Derivation and the measured
/// confirmation are in the module documentation of `subduction_star_decompose`
/// (`FullStarBlock::target_little_character`) and in the card-5 report.
///
/// The derivation is for a **constructed** target, whose Bloch phase is the
/// block's folded point.  A stored target is evaluated at its pinned `k`, where
/// the same division is not proved (R6.7 card 5 math review, P1-1), so this
/// fails closed instead of comparing an unproved quantity.  On this corpus the
/// branch is unreachable — every interior point of every interval reports
/// constructed targets only (measured: `stored=0` of 174,672) — which is why the
/// module test `the_stored_branch_of_the_gauge_fails_closed` is the only
/// witness that the guard fires.
fn gauge_unified_targets(
    block: &FullStarBlock,
    point: &Vec3R,
    operations: &[ExactSeitz],
) -> Result<Vec<SweepTarget>, String> {
    gauge_target_slice(block.targets(), point, operations, |term, operation| {
        block
            .target_little_character(term, point, operation)
            .map_err(|error| error.to_string())
    })
}

/// [`gauge_unified_targets`] with the character reading injected.
///
/// The fail-closed guard against a stored target lives **inside** this function,
/// on the same path the sweep uses, so removing it breaks the guard test — a
/// separately callable `refuse_stored_targets(..)` left the call site unbound
/// (audit of `d66db48`, P0-3: deleting the call kept the gate and every test
/// green).  The injected reader also lets the test prove the guard fires
/// **before** any character is read.
fn gauge_target_slice(
    targets: &[FullStarTarget],
    point: &Vec3R,
    operations: &[ExactSeitz],
    mut read: impl FnMut(usize, &ExactSeitz) -> Result<Complex64, String>,
) -> Result<Vec<SweepTarget>, String> {
    // The identity of the little group, for the character-normalization check
    // below: an honest little-group character takes the value **dimension** there,
    // and the card-1 gauge divides by `exp(2 pi i q . 0) = 1`.
    let identity = operations
        .iter()
        .position(|operation| {
            operation.rotation() == IDENTITY_ROTATION && operation.translation().is_zero()
        })
        .ok_or_else(|| {
            "the aligned little group carries no zero-translation identity, so the character \
             normalization cannot be checked"
                .to_string()
        })?;
    let mut out = Vec::with_capacity(targets.len());
    for (term, target) in targets.iter().enumerate() {
        // The card-1 gauge is derived for constructed targets read at the block's
        // folded point only; a stored target is evaluated at its pinned `k`, where
        // the same division is unproved (R6.7 card 5 math review, P1-1).  The
        // corpus never reaches this branch (interior stored targets 0 of 174,672),
        // so this guard and its test are the only witnesses that it fires.
        if let Some(irnumber) = target.irnumber {
            return Err(format!(
                "the card-1 gauge is derived for constructed targets read at the block's folded \
                 point only, and this block's target {term} is the stored irrep {irnumber}: \
                 divide out its own pinned k phase (or restrict the sweep to the interior of the \
                 stored-k coincidence set) before comparing it across parameters"
            ));
        }
        let mut vector = Vec::with_capacity(operations.len());
        for operation in operations {
            let value = read(term, operation)?;
            let gauge =
                bloch_phase(point, operation.translation()).map_err(|error| error.to_string())?;
            vector.push(value * gauge.conj());
        }
        // Character validity, not just shape: every entry finite, and `chi(E)`
        // equal to the reported dimension.  Without this the matcher's normalized
        // inner product would accept a vector that is not a character at all.
        for (slot, value) in vector.iter().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(format!(
                    "target {term} has a non-finite character entry {value} at operation {slot}"
                ));
            }
        }
        let expected = Complex64::new(f64::from(target.dimension), 0.0);
        if (vector[identity] - expected).norm() > SWEEP_TOLERANCE {
            return Err(format!(
                "target {term} has chi(E) = {} but its reported dimension is {}",
                vector[identity], target.dimension
            ));
        }
        out.push(SweepTarget {
            dimension: target.dimension,
            multiplicity: target.multiplicity,
            stored: false,
            vector,
        });
    }
    Ok(out)
}

/// [`exactly_fixes`] with the control's counters incremented from the rotation
/// that is **actually passed in**.
///
/// Counting the operations the loop *intends* to check does not bind the operand:
/// replacing the argument by the identity left the count and the gate green
/// (verification of `b2d39f3`, P0).  Here the counter sees the argument itself,
/// so the same substitution moves `little_fixity_non_identity_checks` to zero and
/// trips the gate's prediction.
fn check_exact_fixity(
    report: &mut SweepReport,
    rotation: &Mat3I,
    point: &Vec3R,
) -> Result<bool, String> {
    report.little_fixity_checks += 1;
    if *rotation != IDENTITY_ROTATION {
        report.little_fixity_non_identity_checks += 1;
    }
    exactly_fixes(rotation, point)
}

/// The little group built at `point` must fix it **exactly**, not only modulo
/// the reciprocal lattice.
///
/// This is the hypothesis that makes the card-1 gauge parameter-independent at
/// the seed point even when the seed arm is not the block's representative arm
/// (R6.7 card 5 math review, P1-3: the derivation in
/// `FullStarBlock::target_little_character` needs `R_h^T q = q`, and the little
/// group itself is built from the weaker `R_h^{-T} q = q (mod L*)`).  Inside a
/// partition interval the two coincide — a non-exact fixity is exactly a
/// boundary of the partition — so this is a control on the partition as much as
/// on the gauge, and it is counted rather than assumed.
fn exactly_fixes(rotation: &Mat3I, point: &Vec3R) -> Result<bool, String> {
    let action = Mat3R::from_ints(*rotation)
        .inverse()
        .map_err(|error| error.to_string())?
        .transpose();
    let image = action
        .checked_mul_vector(point)
        .map_err(|error| error.to_string())?;
    Ok(image == *point)
}

/// The common seed arm of two matched blocks: the smallest arm index they share
/// (the arm sets are equal by construction, so this is the whole set's minimum).
///
/// An arm index is the stable identity here — folded coordinates, block order,
/// the engine's representative point and the constructed target's enumerated
/// index all move with the parameter.
fn seed_point(block: &FullStarBlock, arm: usize, key: &str) -> Result<Vec3R, String> {
    let mut found: Option<Vec3R> = None;
    for point in block.points() {
        if point.arm_indices().contains(&arm) {
            if found.is_some() {
                return Err(format!(
                    "{key}: arm {arm} folds onto more than one point of one block"
                ));
            }
            found = Some(*point.q());
        }
    }
    found.ok_or_else(|| format!("{key}: the seed arm {arm} is in no point of its own block"))
}

/// Everything the sweep of one `(record, label)` needs beyond the partition.
struct SweepContext<'a> {
    ordinal: usize,
    label: &'static str,
    child_sg: u8,
    child_reciprocal: &'a Lattice,
    child_rotations: &'a [Mat3I],
    /// The production folded directions, in arm order: `arms[a].direction`
    /// scaled by `t` is arm `a`'s folded point at that parameter.
    arms: &'a [FoldedArm],
}

/// Why one matched target pair was rejected; the two kinds are counted
/// separately so the report can say which invariant moved.
#[derive(Debug)]
enum TargetMismatch {
    /// A matched pair whose multiplicities disagree (the character identified
    /// the two targets as the same representation).
    Multiplicity(String),
    /// No partner, a tie at the tolerance, or a dimension disagreement.
    Ambiguous(String),
}

/// One `(record, label)`'s card-5 sweep, appended to the record's report.
///
/// Failures are reported, never skipped: a partition, probe-parameter, rational,
/// engine, geometry or character failure all end up in `report.failures` and
/// make the interval a counted failure, so the gate cannot claim to have checked
/// an interval it lost.  `intervals == comparisons + failed_intervals` is then
/// checkable from the report alone.
fn sweep_label(
    subgroup: &isotropy::IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    partition: &FullStarPartition,
    domain: &ParentDomain,
    context: &SweepContext,
    report: &mut SweepReport,
) {
    let parameters = match partition.probe_parameters(domain) {
        Ok(parameters) => parameters,
        Err(error) => {
            report.failures.push(format!(
                "ordinal {} {}: the sweep cannot list the probe parameters: {error}",
                context.ordinal, context.label
            ));
            return;
        }
    };
    let key = format!("ordinal {} {}", context.ordinal, context.label);
    // A probe set with no boundary at all is the single interval `(0, 1)`, which
    // `probe_intervals` handles; an **empty** one would silently sweep nothing.
    if parameters.is_empty() {
        report.failures.push(format!(
            "ordinal {} {}: the sweep probe set is empty, so no interval would be checked",
            context.ordinal, context.label
        ));
        return;
    }
    let intervals = match probe_intervals(&parameters) {
        Ok(intervals) => intervals,
        Err(error) => {
            report.failures.push(format!(
                "ordinal {} {}: the sweep cannot split the probe set: {error}",
                context.ordinal, context.label
            ));
            return;
        }
    };
    // The probe set has to **be** the partition's, not merely come from the same
    // function: a boundary that is *moved* (not dropped) leaves the interval
    // count, every volume counter and every interval comparison identical, so
    // the card-3 witness would silently stop being cut (card-5 audit F1:
    // replacing 10038 `DT1`'s `1/8` by `1/7` kept the gate at exit 0 with only
    // the printed `parent_only_parameters` moving 0 -> 1).
    let mut expected: Vec<Rat> = partition.boundary_parameters();
    for entry in &domain.exceptional {
        if !expected.contains(&entry.parameter) {
            expected.push(entry.parameter);
        }
    }
    let expected = sorted_parameters(expected);
    if parameters != expected {
        report.failures.push(format!(
            "ordinal {} {}: the sweep's probe parameter set is not the full-star boundaries plus \
             the parent's exceptional parameters ({} entries against {})",
            context.ordinal,
            context.label,
            parameters.len(),
            expected.len()
        ));
        return;
    }
    report.pairs += 1;
    for parameter in &parameters {
        if !partition.is_boundary(parameter) {
            report.parent_only_parameters += 1;
        }
    }
    // R6.7 card 6, part 1: the production decomposition at **every** probe
    // parameter, not only at the interval interiors.  This is what makes the
    // boundary parameters part of the swept evidence instead of only of the
    // separate main-probe pass.  The counters are incremented around the call
    // itself, with no early exit between them, so
    // `boundary_points == boundary_successes + boundary_failures` cannot lose an
    // item (the gate asserts it and pins the total).
    sweep_boundary_pass(subgroup, embedding, table, &parameters, context, report);
    // R6.7 card 6, part 2: the same pair's arm geometry, recomputed from the
    // folded arms and the child reciprocal lattice without the partition's
    // solver; the returned relation list is evaluated again at the interior
    // points of every interval below.
    let relations = audit_arm_geometry(table, partition, parameters.len(), &key, report);
    for (left, right) in intervals {
        report.intervals += 1;
        let mut failures: Vec<String> = Vec::new();
        if let Err(error) = sweep_interval(
            subgroup,
            embedding,
            table,
            context,
            &relations,
            left,
            right,
            report,
            &mut failures,
        ) {
            failures.push(error);
        }
        if failures.is_empty() {
            report.comparisons += 1;
        } else {
            report.failed_intervals += 1;
            report.failures.extend(failures);
        }
    }
}

/// R6.7 card 6: run the **production** decomposition at every probe parameter of
/// one `(record, label)`.
///
/// `boundary_points` is incremented before the engine is asked and exactly one of
/// the two outcome counters after it, so every parameter is accounted for
/// whether it is answered or fails.  A failure carries the record, the label,
/// the parameter and the engine's own error text.
fn sweep_boundary_pass(
    subgroup: &isotropy::IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    parameters: &[Rat],
    context: &SweepContext,
    report: &mut SweepReport,
) {
    for parameter in parameters {
        report.boundary_points += 1;
        match subduce_line_at_parameter(subgroup, embedding, table, *parameter) {
            Ok(result) => {
                report.boundary_successes += 1;
                report.boundary_parent_dimension += u64::from(result.parent_dimension());
                report.boundary_covered_dimension += u64::from(result.covered_dimension());
                let mut failures = Vec::new();
                let dimensions = block_dimensions(table, &result, &mut failures);
                if dimensions.covered != dimensions.parent {
                    failures.push(format!(
                        "ordinal {} {} t={parameter}: the boundary decomposition covers {} of the \
                         parent dimension {}",
                        context.ordinal, context.label, dimensions.covered, dimensions.parent
                    ));
                }
                report.failures.extend(failures);
            }
            Err(error) => {
                report.boundary_failures += 1;
                report.failures.push(format!(
                    "ordinal {} {} t={parameter}: the boundary decomposition failed: {error}",
                    context.ordinal, context.label
                ));
            }
        }
    }
}

/// One interval: two interior points, both decompositions, geometry and
/// representation comparison.
#[allow(clippy::too_many_arguments)]
fn sweep_interval(
    subgroup: &isotropy::IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    context: &SweepContext,
    relations: &[CensusRelation],
    left: Rat,
    right: Rat,
    report: &mut SweepReport,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    let key = format!(
        "ordinal {} {} interval ({left}, {right})",
        context.ordinal, context.label
    );
    let (_, first, second) =
        interior_points(&left, &right).map_err(|error| format!("{key}: interior points: {error}"))?;
    // R6.7 card 6, part 3: the interior half of the arm-geometry audit.  The
    // partition's exactness says that no relation with a nonzero difference
    // vector holds strictly inside an interval; here that is checked by
    // evaluating the whole independently enumerated relation list at both
    // interior points.  Both points are counted before either is evaluated, so
    // the two-per-interval conservation holds even when a relation fails to be
    // decided.
    report.geometry.interior_points += 2;
    for parameter in [first, second] {
        let mut evaluated = report.geometry.interior_evaluated;
        let counts = match relation_counts(relations, &parameter, &mut evaluated) {
            Ok(counts) => counts,
            Err(error) => {
                report.geometry.interior_evaluated = evaluated;
                return Err(format!(
                    "{key}: the interior relation audit at t={parameter}: {error}"
                ));
            }
        };
        report.geometry.interior_evaluated = evaluated;
        let holding: usize = counts.iter().sum();
        if holding > 0 {
            report.geometry.interior_holds += 1;
            failures.push(format!(
                "{key}: at the interior parameter t={parameter} {holding} full-star relation(s) \
                 hold (little-co-group growth {}, arm merge {}, orbit identification {}), so the \
                 parameter is a cut the partition does not report",
                counts[0], counts[1], counts[2]
            ));
        }
    }
    let mut results = Vec::with_capacity(2);
    for parameter in [first, second] {
        report.points += 1;
        match subduce_line_at_parameter(subgroup, embedding, table, parameter) {
            Ok(result) => results.push((parameter, result)),
            Err(error) => {
                return Err(format!(
                    "{key}: the production decomposition failed at t={parameter}: {error}"
                ));
            }
        }
    }
    let (t1, r1) = &results[0];
    let (t2, r2) = &results[1];
    let key = format!("{key} t={t1} against t={t2}");
    // Layer 1: the production geometry at both parameters, read from the
    // engine's own blocks.  `line_star_geometry` is deliberately not consulted
    // here -- the comparison has to see the blocks the decomposition reported.
    let mut geometry = Vec::with_capacity(2);
    for (parameter, result) in [(t1, r1), (t2, r2)] {
        let mut blocks = Vec::with_capacity(result.blocks().len());
        for block in result.blocks() {
            blocks.push(SweepBlock::read(context, *parameter, block)?);
        }
        geometry.push(blocks);
    }
    let matching = match_block_geometry(&geometry[0], &geometry[1], &key)?;
    for (index, partner) in &matching {
        report.blocks += 1;
        let block = &r1.blocks()[*index];
        let other = &r2.blocks()[*partner];
        for target in block.targets() {
            report.sources[usize::from(target.irnumber.is_none())] += 1;
        }
        *report.target_counts.entry(block.targets().len()).or_insert(0) += 1;
        // The seed arm is the smallest arm of the (identical) arm sets; each
        // block's own folded point for it is that side's gauge reference point.
        let Some(arm) = block.arm_indices().first().copied() else {
            failures.push(format!("{key}: block {index} carries no parent arm"));
            continue;
        };
        let q1 = seed_point(block, arm, &key)?;
        let q2 = seed_point(other, arm, &key)?;
        // The partition's arm table and the engine's folded points have to be
        // the same object: `arms[arm].direction * t` must be exactly the point
        // the block reports for that arm -- on **both** sides, each against its
        // own parameter.  Side 2 is the control the earlier version lacked: with
        // both sides read from one decomposition (the `M5` vacuity mutation) the
        // second parameter's point is not `t_2 * direction`, and this fires
        // instead of the comparison quietly agreeing with itself.
        let mut bound = true;
        if q1 == q2 {
            failures.push(format!(
                "{key}: block {index} arm {arm} is at the same point {q1} at both interior \
                 parameters, so comparing this block against itself would be vacuous"
            ));
            bound = false;
        }
        for (side, parameter, point) in [("the first", t1, q1), ("the second", t2, q2)] {
            let expected = scale_point(&context.arms[arm].direction, parameter)
                .map_err(|error| format!("{key}: the seed arm direction: {error}"))?;
            let same = context
                .child_reciprocal
                .same_mod(&expected, &point)
                .map_err(|error| format!("{key}: the seed arm point: {error}"))?;
            if !same {
                failures.push(format!(
                    "{key}: at {side} parameter t={parameter} block {index} arm {arm} is at \
                     {point} but the partition's folded direction says {expected}"
                ));
                bound = false;
            }
        }
        if !bound {
            continue;
        }
        // The same finite group at both ends, aligned by rotation: the
        // operations are parameter-independent child Hall operations, and only
        // the little co-group filter could have moved with the parameter.
        let aligned = match align_little_groups(
            embedding,
            context,
            &q1,
            &q2,
            block.arm_indices().as_slice(),
            &key,
            failures,
        )? {
            Some(aligned) => aligned,
            None => continue,
        };
        *report.orders.entry(aligned.len()).or_insert(0) += 1;
        // The card-1 gauge hypothesis, counted rather than assumed (R6.7 card 5
        // math review, P1-3): each aligned little-group rotation has to fix the
        // seed point **exactly**.  The little group itself is built from the
        // weaker mod-`L*` fixity, so this is the step that makes the gauge
        // parameter-independent even though the seed arm need not be the block's
        // representative arm.
        for (first_operation, second_operation) in &aligned {
            for (side, parameter, point, operation) in [
                ("the first", t1, q1, first_operation),
                ("the second", t2, q2, second_operation),
            ] {
                match check_exact_fixity(report, &operation.rotation(), &point) {
                    Ok(true) => {}
                    Ok(false) => {
                        report.little_fixity_mismatches += 1;
                        failures.push(format!(
                            "{key}: at {side} parameter t={parameter} the aligned little-group \
                             rotation fixes the seed point {point} only modulo the child \
                             reciprocal lattice, so the card-1 gauge is not exact there"
                        ));
                    }
                    Err(error) => failures
                        .push(format!("{key}: the exact-fixity control at t={parameter}: {error}")),
                }
            }
        }
        // The engine's canonical representative is allowed to jump inside an
        // interval (a reduction wall); the comparison is on the gauge-unified
        // content, so this is counted and reported rather than required to be
        // constant.
        if block.stored_k() != other.stored_k() {
            report.representative_shifts += 1;
        }
        // The child little-group order the census computes through
        // `line_domain` must agree with the operation list built here.
        match point_little_co_group_order(context.child_sg, &q1) {
            Ok(order) if order == aligned.len() => {}
            Ok(order) => failures.push(format!(
                "{key}: the census child order at {q1} is {order} but the little group built for \
                 the sweep has {} operation(s)",
                aligned.len()
            )),
            Err(error) => failures.push(format!("{key}: the census child order at {q1}: {error}")),
        }
        // Layer 2: the targets' gauge-unified little-group characters.
        let first_operations: Vec<ExactSeitz> =
            aligned.iter().map(|(operation, _)| *operation).collect();
        let second_operations: Vec<ExactSeitz> =
            aligned.iter().map(|(_, counterpart)| *counterpart).collect();
        // The block's arm set joins the key: a target mismatch has to name the
        // offending block, not only the interval.
        let block_key = format!("{key}: block {index} (arms {:?})", block.arm_indices());
        let targets = gauge_unified_targets(block, &q1, &first_operations)
            .map_err(|error| format!("{block_key}: the first parameter's targets: {error}"))?;
        let other_targets = gauge_unified_targets(other, &q2, &second_operations)
            .map_err(|error| format!("{block_key}: the second parameter's targets: {error}"))?;
        // The per-arm identity that makes the little-group reading faithful: the
        // sum over the block's own points of the little-group characters is the
        // induced character the engine's solver and reconstruction use.  A point
        // that resolves to no arm, or an arm the block's points miss, shows up
        // here instead of silently scoring zero.
        for (parameter, result, operations, side_index) in [
            (t1, r1, &first_operations, *index),
            (t2, r2, &second_operations, *partner),
        ] {
            // The **matched partner's** block on its own side: the pairing is by
            // parent arm set, so the two indices need not be equal (card-5 audit
            // F2 measured `index != partner` on the corpus).  Using the left
            // index on the right decomposition checks a block the pair does not
            // name.
            let side_block = &result.blocks()[side_index];
            for (term, _) in side_block.targets().iter().enumerate() {
                for operation in operations {
                    let mut total = Complex64::new(0.0, 0.0);
                    for star_point in side_block.points() {
                        total += side_block
                            .target_little_character(term, star_point.q(), operation)
                            .map_err(|error| {
                                format!("{key}: the little-group character: {error}")
                            })?;
                    }
                    let induced = side_block
                        .target_character(term, operation)
                        .map_err(|error| format!("{key}: the induced character: {error}"))?;
                    if (total - induced).norm() > SWEEP_TOLERANCE {
                        failures.push(format!(
                            "{key}: at t={parameter} the block's points sum to {total} but the \
                             induced character of target {term} is {induced}"
                        ));
                    }
                }
            }
        }
        // What the multi-target blocks of this corpus actually look like, and
        // whether a multiplicity swap between two distinct targets is observable
        // on them at all: measured rather than assumed, because a mutation that
        // the corpus cannot see must not be reported as caught.
        let mut shapes: Vec<(u8, u32)> = targets
            .iter()
            .map(|target| (target.dimension, target.multiplicity))
            .collect();
        shapes.sort_unstable();
        if shapes.len() > 1 {
            *report.multi_shapes.entry(shapes.clone()).or_insert(0) += 1;
        }
        for (slot, target) in targets.iter().enumerate() {
            if targets[slot + 1..]
                .iter()
                .any(|other| {
                    other.dimension == target.dimension
                        && other.multiplicity != target.multiplicity
                })
            {
                report.same_dimension_swaps += 1;
                break;
            }
        }
        match match_targets(&targets, &other_targets, &block_key) {
            Ok((pairs, matched_worst)) => {
                report.worst_score = report.worst_score.min(matched_worst);
                report.targets += pairs.len();
            }
            Err(TargetMismatch::Multiplicity(error)) => {
                report.multiplicity_mismatches += 1;
                failures.push(error);
            }
            Err(TargetMismatch::Ambiguous(error)) => {
                report.ambiguities += 1;
                failures.push(error);
            }
        }
    }
    Ok(())
}

/// The child little groups at the two seed points, aligned by rotation.
///
/// Returns `None` when the two sides expose different rotation sets (a counted
/// failure is pushed) and an error when a group cannot be built at all.
fn align_little_groups(
    embedding: &SubgroupEmbedding,
    context: &SweepContext,
    q1: &Vec3R,
    q2: &Vec3R,
    arms: &[usize],
    key: &str,
    failures: &mut Vec<String>,
) -> Result<Option<Vec<(ExactSeitz, ExactSeitz)>>, String> {
    let operations = child_little_group(embedding, context.child_reciprocal, q1)
        .map_err(|error| format!("{key}: the child little group at {q1}: {error}"))?;
    let other_operations = child_little_group(embedding, context.child_reciprocal, q2)
        .map_err(|error| format!("{key}: the child little group at {q2}: {error}"))?;
    let mut other_by_rotation: BTreeMap<[i32; 9], ExactSeitz> = other_operations
        .iter()
        .map(|operation| (rotation_key(&operation.rotation()), *operation))
        .collect();
    let mut aligned = Vec::with_capacity(operations.len());
    for operation in &operations {
        match other_by_rotation.remove(&rotation_key(&operation.rotation())) {
            Some(counterpart) => aligned.push((*operation, counterpart)),
            None => {
                failures.push(format!(
                    "{key}: block with arms {arms:?} has a little-group rotation at the first \
                     parameter that the second parameter's little group at {q2} does not have"
                ));
                return Ok(None);
            }
        }
    }
    if !other_by_rotation.is_empty() {
        failures.push(format!(
            "{key}: block with arms {arms:?} has {} little-group rotation(s) at the second \
             parameter that the first parameter's little group at {q1} does not have",
            other_by_rotation.len()
        ));
        return Ok(None);
    }
    Ok(Some(aligned))
}

/// The card-4 full-star recount pass of one `(record, label)`.
///
/// Runs the **production** decomposition at every parameter the card-3 partition
/// adds beyond the reference candidate set and classifies every block there, so
/// the parameter set that the legacy reference-only partition never cut is
/// measured instead of assumed to be equal.  Every failure is pushed into the
/// report, which the caller turns into gate violations.
#[allow(clippy::too_many_arguments)]
fn recount_label(
    subgroup: &isotropy::IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    partition: &FullStarPartition,
    child_candidates: &[Rat],
    generic: Rat,
    context: &BlockContext,
    probes: &[Probe],
    report: &mut RecountReport,
) {
    let mut added: Vec<Rat> = Vec::new();
    for parameter in partition.boundary_parameters() {
        if !child_candidates.contains(&parameter) {
            added.push(parameter);
        }
    }
    // The geometry at the generic sample, taken from the production blocks of the
    // probe the census already ran there: the comparison is engine output against
    // engine output, so a recount parameter that does not change the block
    // structure is visible as such.
    let generic_geometry = probes
        .iter()
        .find(|probe| probe.label == table.label && probe.parameter == generic)
        .map(|probe| block_geometry(&probe.blocks));
    for parameter in added {
        report.probes += 1;
        let result = match subduce_line_at_parameter(subgroup, embedding, table, parameter) {
            Ok(result) => result,
            Err(error) => {
                report.failures.push(format!(
                    "ordinal {} {} t={parameter}: the full-star recount decomposition failed: \
                     {error}",
                    context.ordinal, table.label
                ));
                continue;
            }
        };
        let mut failures = Vec::new();
        let dimensions = block_dimensions(table, &result, &mut failures);
        let (blocks, block_failures) = block_stats(context, parameter, &result, &dimensions);
        failures.extend(block_failures);
        if dimensions.covered != dimensions.parent {
            failures.push(format!(
                "ordinal {} {} t={parameter}: the recount blocks cover {} of the parent \
                 dimension {}",
                context.ordinal, table.label, dimensions.covered, dimensions.parent
            ));
        }
        for failure in failures {
            report.failures.push(failure);
        }
        report.blocks += blocks.len();
        report.rows.push(RecountRow {
            ordinal: context.ordinal,
            parent_sg: context.parent_sg,
            child_sg: context.child_sg,
            label: table.label.to_string(),
            parameter,
            blocks: blocks
                .iter()
                .map(|block| {
                    (
                        block.arm_indices.clone(),
                        block.arm_count,
                        block.block_dimension,
                        block.source,
                    )
                })
                .collect(),
        });
        for block in &blocks {
            report.sources[block.source.index()] += 1;
            if !block.cocycle_trivial {
                report.non_trivial[block.source.index()] += 1;
                *report.non_trivial_orders.entry(block.little_co_group).or_insert(0) += 1;
            }
        }
        let geometry = block_geometry(&blocks);
        if generic_geometry.as_ref().is_some_and(|generic| *generic != geometry) {
            report.changed_pairs.insert((context.ordinal, table.label));
        }
        report.geometries.insert(geometry);
    }
}

/// How many blocks one recount witness decomposition may have; the pinned
/// witnesses are small and a changed structure must not blow this up silently.
const WITNESS_BLOCKS: usize = 64;

/// The two card-3 witnesses, asserted through the card-4 per-block statistics.
///
/// Both are properties of the **engine's own block structure**, measured on the
/// production answer, so a partition that no longer describes the engine cannot
/// satisfy them:
///
/// * ordinal 10038 (SG 196 `DT1` -> P1): six one-armed blocks at `t = 1/9` and
///   four blocks with arm counts `2, 1, 1, 2` at `t = 1/8`.  `1/9` is *not* a
///   partition boundary, so the witness runs the engine there explicitly: the
///   merge at `1/8` is only a change if the neighbouring generic parameter is
///   really unmerged.
/// * ordinal 10030 (SG 196 `DT1` -> #18): at `t = 1/8` a block whose **own**
///   cocycle is non-trivial and whose little co-group has order four -- the
///   change the reference partition cannot see because the reference arm itself
///   stays trivial there.
fn recount_witnesses(
    record: &Record,
    embedding: &SubgroupEmbedding,
    cache: &PointClassCache,
) -> Vec<String> {
    let pinned = match record.ordinal {
        10_038 => (196, "DT1", 1),
        10_030 => (196, "DT1", 18),
        _ => return Vec::new(),
    };
    let mut failures = Vec::new();
    if record.parent_sg != pinned.0 || !record.labels.iter().any(|(label, _)| *label == pinned.1) {
        failures.push(format!(
            "ordinal {}: the pinned witness source SG {} {} is gone (found SG {} with labels {:?})",
            record.ordinal,
            pinned.0,
            pinned.1,
            record.parent_sg,
            record.labels.iter().map(|(label, _)| *label).collect::<Vec<_>>()
        ));
        return failures;
    }
    let label = pinned.1;
    let Some(table) = line_table(record.parent_sg, label) else {
        failures.push(format!("ordinal {}: no frozen table for {label}", record.ordinal));
        return failures;
    };
    if record.child_sg != pinned.2 {
        failures.push(format!(
            "ordinal {}: the pinned witness child is #{}, not #{}",
            record.ordinal, record.child_sg, pinned.2
        ));
    }
    let Ok(child_reciprocal) = reciprocal_lattice(record.child_sg) else {
        failures.push(format!("ordinal {}: no child reciprocal lattice", record.ordinal));
        return failures;
    };
    let Some(direction) = line_direction(table) else {
        failures.push(format!("ordinal {}: no parsable direction for {label}", record.ordinal));
        return failures;
    };
    let Ok(reference_direction) = fold_wave_vector(embedding.transform(), &direction) else {
        failures.push(format!("ordinal {}: folding the direction failed", record.ordinal));
        return failures;
    };
    let context = BlockContext {
        ordinal: record.ordinal,
        parent_sg: record.parent_sg,
        child_sg: record.child_sg,
        label,
        child_reciprocal: &child_reciprocal,
        reference_direction,
        // The witnesses are about the block structure itself, so the reference
        // arm and the Gamma enumeration are not needed; the partition-side
        // versions of both are checked by the recount pass above.
        reference_arm: None,
        gamma: None,
        cache,
    };
    let eighth = Rat::new(1, 8).expect("1/8");
    let ninth = Rat::new(1, 9).expect("1/9");
    for parameter in [ninth, eighth] {
        let result =
            match subduce_line_at_parameter(&record.subgroup, embedding, table, parameter) {
                Ok(result) => result,
                Err(error) => {
                    failures.push(format!(
                        "ordinal {} {label} t={parameter}: the witness decomposition failed: \
                         {error}",
                        record.ordinal
                    ));
                    continue;
                }
            };
        if result.blocks().len() > WITNESS_BLOCKS {
            failures.push(format!(
                "ordinal {} {label} t={parameter}: {} blocks, more than the witness bound {}",
                record.ordinal,
                result.blocks().len(),
                WITNESS_BLOCKS
            ));
            continue;
        }
        let mut dimension_failures = Vec::new();
        let dimensions = block_dimensions(table, &result, &mut dimension_failures);
        failures.extend(dimension_failures);
        let (blocks, block_failures) = block_stats(&context, parameter, &result, &dimensions);
        failures.extend(block_failures);
        if record.ordinal == 10_038 {
            // Six one-armed blocks at the generic `1/9`, four blocks with arm
            // counts `2, 1, 1, 2` at the merge `1/8`.  The expectation is padded
            // with zeros, so a lost block cannot pass as a shorter list.
            let expected: [usize; 6] = if parameter == ninth {
                [1, 1, 1, 1, 1, 1]
            } else {
                [1, 1, 2, 2, 0, 0]
            };
            let mut arms: Vec<usize> = blocks.iter().map(|block| block.arm_count).collect();
            arms.sort_unstable();
            let mut actual = arms.clone();
            actual.resize(expected.len(), 0);
            if actual != expected.to_vec() {
                failures.push(format!(
                    "ordinal {} {label} t={parameter}: the witness expects block arm counts {:?} \
                     but the engine reports {:?} ({} block(s))",
                    record.ordinal,
                    expected,
                    arms,
                    blocks.len()
                ));
            }
        }
        if record.ordinal == 10_030 {
            // The class change is the witness here: order four and non-trivial
            // at `1/8`, and nothing non-trivial at the generic `1/9`.
            let non_trivial: Vec<(usize, bool, BlockSource)> = blocks
                .iter()
                .filter(|block| !block.cocycle_trivial)
                .map(|block| (block.little_co_group, block.cocycle_trivial, block.source))
                .collect();
            if parameter == ninth {
                if !non_trivial.is_empty() {
                    failures.push(format!(
                        "ordinal {} {label} t={parameter}: the generic control already carries a \
                         non-trivial block class: {non_trivial:?}",
                        record.ordinal
                    ));
                }
            } else {
                match non_trivial.iter().find(|(order, _, _)| *order == 4) {
                    Some((_, _, source)) if *source != BlockSource::Stored => failures.push(format!(
                        "ordinal {} {label} t={parameter}: the order-four non-trivial block is {}, \
                         not stored",
                        record.ordinal,
                        source.label()
                    )),
                    Some(_) => {}
                    None => failures.push(format!(
                        "ordinal {} {label} t={parameter}: no block has a non-trivial cocycle of \
                         order four (non-trivial blocks: {non_trivial:?})",
                        record.ordinal
                    )),
                }
            }
        }
    }
    failures
}
/// Probe one record at every parameter of its partition.
///
/// The child little co-group order is recomputed at **every** probed parameter
/// with [`little_co_group_order`], not only at the parameters the child census
/// lists, so no probe ever carries an unknown order.  The child-side enumeration
/// is also cross-checked against the group on a uniform grid right here.
///
/// Every answered probe also carries one [`BlockStat`] per `result.blocks()`
/// entry, and -- when `recount` is set -- the card-3 full-star partition drives
/// an extra pass over the parameters it adds to the reference candidate set.
fn probe_record(
    record: &Record,
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    official: Rat,
    generic: Rat,
    recount: bool,
    domain_sweep: bool,
    cache: &PointClassCache,
) -> RecordReport {
    let mut out = Vec::new();
    let mut failures = Vec::new();
    let mut child_grid = (0usize, 0usize);
    let mut child_algorithms = (0usize, 0usize);
    let mut child_parameters: Vec<Rat> = Vec::new();
    let mut gamma = GammaReport::default();
    let mut recount_report = RecountReport::default();
    let mut sweep_report = SweepReport {
        worst_score: 1.0,
        ..SweepReport::default()
    };
    let frame = probe_frame(record, &mut failures);
    let Some((embedding, child_reciprocal, child_rotations)) = frame else {
        return RecordReport {
            probes: out,
            child_grid,
            child_algorithms: (0, 0),
            child_parameters,
            gamma,
            recount: recount_report,
            sweep: sweep_report,
            failures,
        };
    };
    // The child frame must be the child's own: the lattice and the rotation set
    // are rebuilt from the embedding's subgroup, the lattice must be invariant
    // under every rotation used, and the record's subgroup number must agree.
    // Neither the grid check nor the two-algorithm comparison can see a wrong
    // frame (both are fed the same lattice), so this is the control that does;
    // the mutations of the third review round are what it exists for.
    //
    // Two limits, both measured.  (1) The census consumes the child frame through
    // the reciprocal lattice, the rotation set and the embedding's transform (the
    // transform is not a function of the child number, so a sibling swap keeps
    // it); a *sibling* space group that shares the lattice and the rotation set —
    // 109 of the 116 child space groups have one — is therefore
    // indistinguishable here.  The engine does consume the whole embedding, but a
    // sibling frame cannot reach it: `SubgroupEmbedding::build` refuses a record
    // whose numbers disagree (`StaleIsotropyRecord`), which also makes the
    // child-number comparison just below a restatement rather than a control.
    // The **parent** frame has the same blind spot and no independent source
    // inside this gate: dropping a rotation from `rotation_set(parent_sg)` outside
    // the generic stabiliser leaves the gate at exit 0 with an unchanged summary
    // line (measured: parent comparisons 2,580 -> 2,571, and the published parent
    // orders in the TSV change), and only the module tests
    // `every_space_group_rotation_set_is_its_stored_hall_settings` and
    // `every_frozen_source_is_exceptional_only_at_zero_and_one_half` catch it.
    // (2) Both frame values below are compared with *themselves*: after the two
    // space group numbers are known to agree, `reciprocal_lattice(child_sg)` and
    // `reciprocal_lattice(embedding.subgroup_sg())` are the same function of the
    // same argument, and so are the two `rotation_set` calls.  Neither comparison
    // can see its function being wrong.  That is covered for the lattice by the
    // probe-side order comparison against the census's own child order (the only
    // place the two lattice constructions meet), and for both values by the
    // per-space-group unit tests `every_space_group_reciprocal_lattice_is_...` and
    // `every_space_group_rotation_set_is_its_stored_hall_settings`.  Measured in
    // the fifth review round: an A/C swap of the SG 38 lattice and SG 38 -> Z^3
    // both pass this gate unchanged (exit 0, identical counts) and are caught by
    // the unit tests and by the order comparison; dropping one rotation from
    // SG 38's rotation set passes this gate unchanged and is caught by the
    // rotation-set unit test alone.
    if record.child_sg != embedding.subgroup_sg() {
        failures.push(format!(
            "ordinal {}: record child #{} != embedding subgroup #{}",
            record.ordinal,
            record.child_sg,
            embedding.subgroup_sg()
        ));
    }
    match reciprocal_lattice(embedding.subgroup_sg()) {
        Ok(expected) => {
            for (name, left, right) in [
                ("rebuilt", &child_reciprocal, &expected),
                ("expected", &expected, &child_reciprocal),
            ] {
                for row in 0..3 {
                    let vector = Vec3R::new(*right.rows().row(row));
                    // Card 6: this lattice error used to be swallowed into
                    // `false` by `unwrap_or(false)`, which would have reported
                    // "not the embedding's lattice" instead of the real error.
                    // It is propagated now.
                    match left.contains(&vector) {
                        Ok(true) => {}
                        Ok(false) => failures.push(format!(
                            "ordinal {}: the {name} child lattice is not the embedding's \
                             reciprocal lattice (row {row})",
                            record.ordinal
                        )),
                        Err(error) => failures.push(format!(
                            "ordinal {}: the {name} child lattice test (row {row}): {error}",
                            record.ordinal
                        )),
                    }
                }
            }
        }
        Err(error) => failures.push(format!(
            "ordinal {}: no reciprocal lattice for child #{}: {error}",
            record.ordinal,
            embedding.subgroup_sg()
        )),
    }
    match rotation_set(embedding.subgroup_sg()) {
        Ok(expected) => {
            if child_rotations != expected {
                failures.push(format!(
                    "ordinal {}: the child rotation set is not the embedding subgroup's",
                    record.ordinal
                ));
            }
        }
        Err(error) => failures.push(format!(
            "ordinal {}: no rotation set for child #{}: {error}",
            record.ordinal,
            embedding.subgroup_sg()
        )),
    }
    for rotation in &child_rotations {
        let Ok(action) = cryspglib::irrep::subduction::Mat3R::from_ints(*rotation).inverse()
        else {
            failures.push(format!("ordinal {}: child rotation is singular", record.ordinal));
            continue;
        };
        let Ok(action) = action.transpose().inverse() else {
            failures.push(format!("ordinal {}: child rotation is singular", record.ordinal));
            continue;
        };
        for row in 0..3 {
            let basis = Vec3R::new(*child_reciprocal.rows().row(row));
            match action.checked_mul_vector(&basis) {
                Ok(image) => match child_reciprocal.contains(&image) {
                    Ok(true) => {}
                    Ok(false) => failures.push(format!(
                        "ordinal {}: child rotation {rotation:?} does not preserve the child \
                         lattice (row {row})",
                        record.ordinal
                    )),
                    Err(error) => failures.push(format!(
                        "ordinal {}: child lattice test failed: {error}",
                        record.ordinal
                    )),
                },
                Err(error) => failures.push(format!(
                    "ordinal {}: child rotation action failed: {error}",
                    record.ordinal
                )),
            }
        }
    }
    for (label, pinned) in &record.labels {
        let Some(table) = line_table(record.parent_sg, label) else {
            failures.push(format!(
                "ordinal {}: no frozen table for SG {} {label}",
                record.ordinal, record.parent_sg
            ));
            continue;
        };
        let Some(domain) = domains.get(&(record.parent_sg, label)) else {
            failures.push(format!(
                "ordinal {}: no census domain for SG {} {label}",
                record.ordinal, record.parent_sg
            ));
            continue;
        };
        let Ok(direction) = line_direction(table).ok_or(()) else {
            failures.push(format!(
                "ordinal {}: SG {} {label} has no parsable direction",
                record.ordinal, record.parent_sg
            ));
            continue;
        };
        // The folded direction is `q(t) = t . q1`; both the grid cross-check and
        // the per-parameter order use `q1`.
        let Ok(folded) = fold_wave_vector(embedding.transform(), &direction) else {
            failures.push(format!("ordinal {}: folding failed", record.ordinal));
            continue;
        };
        if !folded.is_zero()
            && let Err(error) =
                require_reciprocal_direction(&child_reciprocal, &folded, record.child_sg, "")
        {
            failures.push(format!(
                "ordinal {}: the folded direction is out of scope: {error}",
                record.ordinal
            ));
            continue;
        }
        // The child frame gets the same arithmetic-independent control as the
        // parent: the centring scan and the lattice coordinate map must produce
        // the same step for every child rotation.  This is what makes the
        // child-side claim independent rather than a second reading of the same
        // algebra (a wrong *frame* is caught by the frame checks, a wrong step
        // by this comparison).
        //
        // The same loop also checks the premise of the domain theorem (see the
        // module documentation of `line_domain`): a rotation fixes the *whole*
        // direction, `w_R = R^{-T} q1 - q1 in L*_child`, exactly when its step is
        // one.  Those are the rotations whose factor system stays a cocycle for
        // every real parameter, which is what makes the projective class trivial
        // on the whole domain.  The count is then compared with the child order at
        // the generic sample below, which is the same number by an independent
        // route (group membership rather than the step solver).
        let mut generic_rotations = 0usize;
        if !folded.is_zero() {
            for rotation in &child_rotations {
                let image = match cryspglib::irrep::subduction::Mat3R::from_ints(*rotation)
                    .inverse()
                    .and_then(|matrix| matrix.transpose().checked_mul_vector(&folded))
                {
                    Ok(image) => image,
                    Err(error) => {
                        failures.push(format!(
                            "ordinal {}: child rotation image: {error}",
                            record.ordinal
                        ));
                        continue;
                    }
                };
                let Ok(w) = image.checked_sub(&folded) else {
                    failures.push(format!(
                        "ordinal {}: child rotation difference failed",
                        record.ordinal
                    ));
                    continue;
                };
                // The theorem's premise is **exact fixity** of the whole direction,
                // `w_R = 0`: only then is `chi_{t q1} . lambda` a factor system for
                // every real parameter, which is what makes the projective class
                // trivial on the whole domain.  Membership of `w_R` in `L*_child`
                // was the first, wrong formulation of this check; on this corpus it
                // is **vacuous** rather than weaker (the scope fence plus every
                // rotation preserving the lattice makes it true for every rotation:
                // 28,713/28,713), and it does not make the factor system a cocycle
                // for every real parameter -- on that rotation set the cocycle
                // identity fails at 4,877 of 51,804 (pair, parameter) probes, while
                // on the exact stabiliser it holds 51,804/51,804.  The figure
                // "24,430 violations" quoted here before is not reproducible as
                // described (the described count is 5,117 pairs) and is withdrawn.
                // The step solver is not used as
                // the premise either: `Some(1)` is also returned for `w` in `Z^3`
                // outside `L*`, where the constraint binds exactly at the integer
                // parameters and the residues `{0}` are correct.
                if w.is_zero() {
                    generic_rotations += 1;
                }
                let step = match minimal_parameter_step(&child_reciprocal, &w) {
                    Ok(step) => step,
                    Err(error) => {
                        failures.push(format!(
                            "ordinal {} {label}: child rotation step: {error}",
                            record.ordinal
                        ));
                        continue;
                    }
                };
                if w.is_zero() && !step_is_vacuous_or_one(&step) {
                    failures.push(format!(
                        "ordinal {} {label}: the rotation {rotation:?} fixes the whole direction \
                         but the step solver returns {step:?}",
                        record.ordinal
                    ));
                }
                match (
                    minimal_parameter_step(&child_reciprocal, &w),
                    minimal_parameter_step_via_coordinates(&child_reciprocal, &w),
                ) {
                    (Ok(scan), Ok(coordinates)) => {
                        child_algorithms.0 += 1;
                        if scan != coordinates {
                            child_algorithms.1 += 1;
                            failures.push(format!(
                                "ordinal {} child rotation {rotation:?}: centring scan {:?} != \
                                 coordinate map {:?}",
                                record.ordinal,
                                scan.map(|value| value.to_string()),
                                coordinates.map(|value| value.to_string())
                            ));
                        }
                    }
                    (scan, coordinates) => failures.push(format!(
                        "ordinal {} child rotation {rotation:?}: step failed ({:?} / {:?})",
                        record.ordinal,
                        scan.err().map(|error| error.to_string()),
                        coordinates.err().map(|error| error.to_string())
                    )),
                }
            }
        }
        // The premise count against the census's own generic order: the rotations
        // that fix the whole direction are exactly the generic little co-group, so
        // the child order at the generic sample must equal `generic_rotations`.
        match little_co_group_order(&child_reciprocal, &folded, &child_rotations, generic) {
            Ok(order) => {
                if order != generic_rotations {
                    failures.push(format!(
                        "ordinal {} {label}: the generic child order {order} differs from the \
                         {generic_rotations} rotations that fix the whole direction",
                        record.ordinal
                    ));
                }
            }
            Err(error) => failures.push(format!(
                "ordinal {} {label}: generic child order: {error}",
                record.ordinal
            )),
        }
        for denominator in [24i128, 120] {
            match verify_against_grid(&child_reciprocal, &folded, &child_rotations, denominator) {
                Ok(check) => {
                    child_grid.0 += check.checks;
                    child_grid.1 += check.mismatches;
                }
                Err(error) => failures.push(format!(
                    "ordinal {} child grid 1/{denominator}: {error}",
                    record.ordinal
                )),
            }
        }
        let (parameters, child_candidates) =
            match record_parameters(record, domain, table, official, generic) {
                Ok(partition) => partition,
                Err(error) => {
                    // Counted, with the reason: this used to read "parameter
                    // partition failed" and drop the cause, which is the
                    // difference between a diagnosable corpus gap and a silent
                    // one (card 6).
                    failures.push(format!(
                        "ordinal {} {label}: the parameter partition failed: {error}",
                        record.ordinal
                    ));
                    continue;
                }
            };
        for parameter in &child_candidates {
            if !child_parameters.contains(parameter) {
                child_parameters.push(*parameter);
            }
        }
        // The card-3 full-star partition of this (record, label): the boundary
        // set the recount pass runs on and the Gamma-reaching arms the block
        // statistics are cross-checked against.  A failure here is a counted
        // failure, never a skipped pass.
        let mut partition: Option<FullStarPartition> = None;
        let mut gamma_entries: Option<GammaParameters> = None;
        if recount || domain_sweep {
            match full_star_partition(&record.subgroup, &embedding, table) {
                Ok(built) => {
                    let (entries, zero_arms) =
                        report_gamma(&built, record.ordinal, &mut gamma, &mut failures);
                    gamma_entries = Some((entries, zero_arms));
                    partition = Some(built);
                }
                Err(error) => failures.push(format!(
                    "ordinal {} {label}: the full-star partition failed: {error}",
                    record.ordinal
                )),
            }
        }
        let block_context = BlockContext {
            ordinal: record.ordinal,
            parent_sg: record.parent_sg,
            child_sg: record.child_sg,
            label,
            child_reciprocal: &child_reciprocal,
            reference_direction: folded,
            reference_arm: partition.as_ref().map(|built| built.reference_arm),
            gamma: gamma_entries
                .as_ref()
                .map(|(entries, zero_arms)| (entries.as_slice(), zero_arms.as_slice())),
            cache,
        };
        for (parameter, census_order) in parameters {
            // The reference folded point the engine binding needs, computed with
            // the same helper the statistics use (`scale_point`).
            let reference_point_for_engine = match scale_point(&folded, &parameter) {
                Ok(point) => point,
                Err(error) => {
                    failures.push(format!(
                        "ordinal {} {label} t={parameter}: the reference folded point: {error}",
                        record.ordinal
                    ));
                    continue;
                }
            };
            let child_order = if folded.is_zero() {
                child_rotations.len()
            } else {
                match little_co_group_order(&child_reciprocal, &folded, &child_rotations, parameter)
                {
                    Ok(order) => order,
                    Err(error) => {
                        failures.push(format!(
                            "ordinal {} {label} t={parameter}: child order: {error}",
                            record.ordinal
                        ));
                        continue;
                    }
                }
            };
            // The probe-side order comes from `reciprocal_lattice(child_sg)`, the
            // census-side one from `child_candidates`, which builds its lattice
            // inline from `exact_primitive_basis`.  This is the only place in the
            // gate where the two constructions meet, so it is what makes a wrong
            // child lattice observable here (a parameter the child census does not
            // list carries order zero and is skipped; the per-space-group
            // correctness of the function itself is pinned by
            // `line_domain::tests::every_space_group_reciprocal_lattice_is_integral_with_small_exponent`).
            if census_order > 0 && census_order != child_order {
                failures.push(format!(
                    "ordinal {} {label} t={parameter}: probe child order {child_order} != census \
                     child order {census_order}",
                    record.ordinal
                ));
            }
            // The projective class of the folded direction at this parameter.  At a
            // parameter that is **not** one of the child's own exceptional
            // parameters the little co-group is the exact stabiliser of the folded
            // direction, and the domain theorem (`line_domain` module
            // documentation) then makes the class trivial — the same conclusion for
            // every real parameter, not just this one.  That is a direct test of
            // the theorem on every generic probe; the exceptional parameters carry
            // a larger little co-group and get the separate measurement below.
            match child_cocycle_is_a_coboundary(record.child_sg, &folded, parameter) {
                Ok(trivial) => {
                    if census_order == 0 && !trivial {
                        failures.push(format!(
                            "ordinal {} {label} t={parameter}: the projective class is non-trivial \
                             at a generic parameter, contradicting the domain theorem",
                            record.ordinal
                        ));
                    }
                }
                Err(error) => failures.push(format!(
                    "ordinal {} {label} t={parameter}: projective class: {error}",
                    record.ordinal
                )),
            }
            let parent = domain
                .exceptional
                .iter()
                .find(|entry| entry.parameter == parameter);
            let result = subduce_line_at_parameter(
                &record.subgroup,
                &embedding,
                table,
                parameter,
            );
            let (class, parameter_kind, content, detail, blocks, dimensions, engine_blocks) =
                match result {
                Ok(result) => {
                    let dimensions = block_dimensions(table, &result, &mut failures);
                    let (blocks, block_failures) =
                        block_stats(&block_context, parameter, &result, &dimensions);
                    failures.extend(block_failures);
                    // The engine's own reading of the same blocks, kept so the
                    // gate can bind every recorded statistic to it.
                    let mut engine_blocks: Vec<EngineBlock> = Vec::with_capacity(result.blocks().len());
                    for block in result.blocks() {
                        match engine_block(
                            block,
                            &reference_point_for_engine,
                            block_context.child_reciprocal,
                        ) {
                            Ok(engine) => engine_blocks.push(engine),
                            Err(error) => failures.push(format!(
                                "ordinal {} {label} t={parameter}: the engine's own block reading: \
                                 {error}",
                                record.ordinal
                            )),
                        }
                    }
                    let targets: Vec<_> = result
                        .blocks()
                        .iter()
                        .flat_map(|block| block.targets())
                        .collect();
                    // The **withdrawn** probe-level rule, measured alongside the
                    // per-block sources so the two readings can be compared on
                    // the same run (`--output-blocks` carries the block one).
                    let stored = classify_probe(&targets) == BlockSource::Stored;
                    let content = result.trivial_content().map_err(|error| error.to_string());
                    // `content.as_ref().ok()` is not a swallowed error: the
                    // `Err` arm is written into `detail` just below, and
                    // `check_invariants` turns an answered probe without a
                    // trivial content into a violation (guard 5), so a content
                    // failure cannot pass as a missing number.
                    (
                        if stored { TargetClass::Stored } else { TargetClass::Constructed },
                        Some(result.parameter_kind()),
                        content.as_ref().ok().copied(),
                        format!(
                            "blocks={} targets={} {}",
                            result.blocks().len(),
                            targets.len(),
                            match &content {
                                Ok(value) => format!("content={value}"),
                                Err(error) => format!("content error: {error}"),
                            }
                        ),
                        blocks,
                        Some(dimensions),
                        engine_blocks,
                    )
                }
                Err(FullStarError::MissingChildStarData { sg, points, .. }) => (
                    TargetClass::Unsupported,
                    None,
                    None,
                    format!("missing child data: sg={sg} points={points}"),
                    Vec::new(),
                    None,
                    Vec::new(),
                ),
                Err(error) => (
                    TargetClass::Error,
                    None,
                    None,
                    error.to_string(),
                    Vec::new(),
                    None,
                    Vec::new(),
                ),
            };
            let anchor_mismatch = if parameter == official {
                match (content, *pinned) {
                    (Some(value), expected) if value != u32::from(expected) => Some(format!(
                        "anchor content {value} != pinned {expected}"
                    )),
                    (None, _) if class != TargetClass::Error => {
                        Some("anchor has no trivial content".to_string())
                    }
                    _ => None,
                }
            } else {
                None
            };
            out.push(Probe {
                ordinal: record.ordinal,
                parent_sg: record.parent_sg,
                child_sg: record.child_sg,
                label,
                parameter,
                parent_order: parent.map_or(domain.generic_order, |entry| entry.order),
                parent_added: parent.map_or(0, |entry| entry.added),
                child_order,
                class,
                parameter_kind,
                content,
                detail: match anchor_mismatch {
                    Some(note) => format!("{note}; {}", detail),
                    None => detail,
                },
                blocks,
                engine_blocks,
                dimensions,
                partition_arms: partition.as_ref().map(|built| built.arms.len()),
            });
        }
        // The recount pass: every parameter the card-3 partition adds beyond the
        // reference candidate set, every block of the production answer there.
        if let Some(partition) = &partition {
            recount_label(
                &record.subgroup,
                &embedding,
                table,
                partition,
                &child_candidates,
                generic,
                &block_context,
                &out,
                &mut recount_report,
            );
            // R6.7 card 5: two interior points of every interval of
            // `probe_parameters(parent)`, compared by block identity and by
            // gauge-unified little-group characters.
            if domain_sweep {
                let sweep_context = SweepContext {
                    ordinal: record.ordinal,
                    label,
                    child_sg: record.child_sg,
                    child_reciprocal: &child_reciprocal,
                    child_rotations: &child_rotations,
                    arms: &partition.arms,
                };
                sweep_label(
                    &record.subgroup,
                    &embedding,
                    table,
                    partition,
                    domain,
                    &sweep_context,
                    &mut sweep_report,
                );
            }
        }
    }
    if recount {
        failures.extend(recount_witnesses(record, &embedding, cache));
    }
    child_parameters.sort_by_key(|value| value.numerator() * 10_080 / value.denominator());
    child_parameters.dedup();
    RecordReport {
        probes: out,
        child_grid,
        child_algorithms,
        child_parameters,
        gamma,
        recount: recount_report,
        sweep: sweep_report,
        failures,
    }
}

/// Human-readable summary: the partition, the classes and the open boundary.
/// The census evidence that does not live in a single probe.
struct CensusEvidence<'a> {
    /// `(checks, mismatches)` of the child-side grid cross-check over all records.
    child_grid: (usize, usize),
    /// `(checks, mismatches)` of the child-side two-algorithm comparison.
    child_algorithms: (usize, usize),
    /// The union of the records' child candidate parameters.
    child_union: &'a [Rat],
    /// How often the centring scan and the coordinate map route were compared,
    /// and how often they disagreed.
    algorithm_checks: usize,
    algorithm_mismatches: usize,
    /// Data rows the `--output-blocks` emission loop produced.  Checked against
    /// the collected block count on **every** run, whether or not a file was
    /// asked for, so dropping a row in that loop cannot pass unnoticed.
    block_rows: usize,
    /// Data rows read back from the written `--output-blocks` file, when one was
    /// written; the file is re-read instead of trusting the writer's own count,
    /// and each row is parsed and compared with its statistic.
    block_file_rows: Option<usize>,
    /// The Gamma-reaching report, asserted rather than only printed: the card-4
    /// audit showed that dropping the `t = 0` entries left a self-contradictory
    /// printout and a green gate.
    gamma: &'a GammaReport,
    /// The full-star recount aggregate, and whether the pass was requested at all.
    /// Without the "ran" flag a mutation that made the pass vacuous left the gate
    /// green while it printed `full-star recount: not run` (verification review F6).
    recount: &'a RecountReport,
    recount_ran: bool,
    /// Data rows the `--output-recount` emission loop produced, and the rows read
    /// back from the written file (parsed and compared, not only counted).
    recount_rows_written: usize,
    recount_file_rows: Option<usize>,
    /// The R6.7 card-5 sweep aggregate, and whether `--domain-sweep` asked for it.
    sweep: &'a SweepReport,
    sweep_ran: bool,
}

/// Compare one written `--output-blocks` row against the statistic it came from.
///
/// Field by field, by **parsing** the text: re-formatting the expected row with
/// the writer's own formatter would only re-run the writer.
fn check_block_row(line: &str, block: &BlockStat, row: usize) -> Result<(), String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 16 {
        return Err(format!("{row}: {} field(s), expected 16", fields.len()));
    }
    let arms = block
        .arm_indices
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let terms = block
        .terms
        .iter()
        .map(|(dimension, multiplicity)| format!("{dimension}x{multiplicity}"))
        .collect::<Vec<_>>()
        .join(",");
    let expected: [String; 16] = [
        block.ordinal.to_string(),
        block.parent_sg.to_string(),
        block.child_sg.to_string(),
        block.label.to_string(),
        block.parameter.to_string(),
        block.index.to_string(),
        block.star_size.to_string(),
        block.arm_count.to_string(),
        arms,
        block.block_dimension.to_string(),
        block.little_co_group.to_string(),
        if block.cocycle_trivial {
            "trivial".to_string()
        } else {
            "non-trivial".to_string()
        },
        block.source.label().to_string(),
        terms,
        block.carries_reference.to_string(),
        block.carries_gamma.to_string(),
    ];
    for (slot, (found, want)) in fields.iter().zip(expected.iter()).enumerate() {
        if found != want {
            return Err(format!("{row} field {slot}: {found:?} != {want:?}"));
        }
    }
    Ok(())
}

/// Add one record's Gamma report to the corpus-wide one.
fn merge_gamma(total: &mut GammaReport, part: GammaReport) {
    total.entries += part.entries;
    total.contained += part.contained;
    total.exceptional += part.exceptional;
    total.zero_arms += part.zero_arms;
    for (shape, count) in part.shapes {
        *total.shapes.entry(shape).or_insert(0) += count;
    }
    for witness in part.witnesses {
        if total.witnesses.len() < GAMMA_WITNESSES {
            total.witnesses.push(witness);
        }
    }
}

/// Add one record's card-5 sweep report to the corpus-wide one.
fn merge_sweep(total: &mut SweepReport, part: SweepReport) {
    total.pairs += part.pairs;
    total.intervals += part.intervals;
    total.points += part.points;
    total.comparisons += part.comparisons;
    total.failed_intervals += part.failed_intervals;
    total.parent_only_parameters += part.parent_only_parameters;
    total.blocks += part.blocks;
    total.targets += part.targets;
    for (slot, count) in total.sources.iter_mut().zip(part.sources) {
        *slot += count;
    }
    total.multiplicity_mismatches += part.multiplicity_mismatches;
    total.same_dimension_swaps += part.same_dimension_swaps;
    for (shape, count) in part.multi_shapes {
        *total.multi_shapes.entry(shape).or_insert(0) += count;
    }
    total.ambiguities += part.ambiguities;
    for (count, seen) in part.target_counts {
        *total.target_counts.entry(count).or_insert(0) += seen;
    }
    for (order, seen) in part.orders {
        *total.orders.entry(order).or_insert(0) += seen;
    }
    total.representative_shifts += part.representative_shifts;
    total.little_fixity_checks += part.little_fixity_checks;
    total.little_fixity_non_identity_checks += part.little_fixity_non_identity_checks;
    total.little_fixity_mismatches += part.little_fixity_mismatches;
    total.boundary_points += part.boundary_points;
    total.boundary_successes += part.boundary_successes;
    total.boundary_failures += part.boundary_failures;
    total.boundary_parent_dimension += part.boundary_parent_dimension;
    total.boundary_covered_dimension += part.boundary_covered_dimension;
    merge_geometry(&mut total.geometry, part.geometry);
    total.worst_score = total.worst_score.min(part.worst_score);
    total.failures.extend(part.failures);
}

/// Add one record's card-6 arm-geometry audit to the corpus-wide one.
fn merge_geometry(total: &mut GeometryReport, part: GeometryReport) {
    total.pairs += part.pairs;
    total.relations += part.relations;
    total.predicted_relations += part.predicted_relations;
    total.permanent += part.permanent;
    total.permanent_partition += part.permanent_partition;
    total.boundaries += part.boundaries;
    total.recomputed_boundaries += part.recomputed_boundaries;
    total.evaluated += part.evaluated;
    total.interior_points += part.interior_points;
    total.interior_evaluated += part.interior_evaluated;
    total.interior_holds += part.interior_holds;
    total.disagreements += part.disagreements;
    total.arm_orbit_checks += part.arm_orbit_checks;
    total.arm_orbit_disagreements += part.arm_orbit_disagreements;
    total.arms += part.arms;
    total.predicted_boundary_dimension += part.predicted_boundary_dimension;
}

/// Add one record's recount report to the corpus-wide one.
fn merge_recount(total: &mut RecountReport, part: RecountReport) {
    total.probes += part.probes;
    total.blocks += part.blocks;
    for (slot, count) in total.sources.iter_mut().zip(part.sources) {
        *slot += count;
    }
    for (slot, count) in total.non_trivial.iter_mut().zip(part.non_trivial) {
        *slot += count;
    }
    for (order, count) in part.non_trivial_orders {
        *total.non_trivial_orders.entry(order).or_insert(0) += count;
    }
    total.geometries.extend(part.geometries);
    total.changed_pairs.extend(part.changed_pairs);
    total.rows.extend(part.rows);
    total.failures.extend(part.failures);
}

fn report(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    records: &[Record],
    probes: &[Probe],
    errors: &[String],
    evidence: &CensusEvidence,
    gamma: &GammaReport,
    recount: &RecountReport,
) {
    println!(
        "sources={} records={} probes={} probe_errors={}",
        domains.len(),
        records.len(),
        probes.len(),
        errors.len()
    );
    println!(
        "child-side grid cross-check: {} predicates, {} mismatch(es); candidate parameters {{{}}}",
        evidence.child_grid.0,
        evidence.child_grid.1,
        evidence
            .child_union
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
    println!(
        "probes without a child order: {}",
        probes.iter().filter(|probe| probe.child_order == 0).count()
    );
    println!(
        "step algorithms: parent {} rotation comparisons / {} disagreement(s); \
         child {} / {}",
        evidence.algorithm_checks,
        evidence.algorithm_mismatches,
        evidence.child_algorithms.0,
        evidence.child_algorithms.1
    );
    let mut parent_shapes: BTreeMap<Vec<(String, usize, usize)>, usize> = BTreeMap::new();
    for domain in domains.values() {
        let shape: Vec<(String, usize, usize)> = domain
            .exceptional
            .iter()
            .map(|entry| (entry.parameter.to_string(), entry.order, entry.added))
            .collect();
        *parent_shapes.entry(shape).or_insert(0) += 1;
    }
    println!("parent-exceptional parameter shapes (source count):");
    for (shape, count) in &parent_shapes {
        let text: Vec<String> = shape
            .iter()
            .map(|(parameter, order, added)| format!("t={parameter} order={order} +{added}"))
            .collect();
        println!("  {count:>3} source(s): {}", text.join(", "));
    }
    let mut per_source: BTreeMap<(&'static str, u8), Vec<Rat>> = BTreeMap::new();
    for probe in probes {
        let entry = per_source.entry((probe.label, probe.parent_sg)).or_default();
        if !entry.contains(&probe.parameter) {
            entry.push(probe.parameter);
        }
    }
    for parameters in per_source.values_mut() {
        *parameters = sorted_parameters(std::mem::take(parameters));
    }
    let mut partition_shapes: BTreeMap<Vec<String>, usize> = BTreeMap::new();
    for parameters in per_source.values() {
        let shape: Vec<String> = parameters.iter().map(|value| value.to_string()).collect();
        *partition_shapes.entry(shape).or_insert(0) += 1;
    }
    println!("probed parameter sets (source count):");
    for (shape, count) in partition_shapes.iter().take(12) {
        println!("  {count:>3} source(s): {{{}}}", shape.join(", "));
    }
    if partition_shapes.len() > 12 {
        println!("  ... {} more distinct sets", partition_shapes.len() - 12);
    }
    let mut classes: BTreeMap<(&'static str, String), usize> = BTreeMap::new();
    for probe in probes {
        *classes
            .entry((probe.class.label(), probe.parameter.to_string()))
            .or_insert(0) += 1;
    }
    println!("class by parameter:");
    for ((class, parameter), count) in &classes {
        println!("  t={parameter:<6} {class:<12} {count}");
    }
    let unsupported: Vec<&Probe> = probes
        .iter()
        .filter(|probe| probe.class == TargetClass::Unsupported)
        .collect();
    if !unsupported.is_empty() {
        let mut by_child: BTreeMap<(u8, usize), usize> = BTreeMap::new();
        for probe in &unsupported {
            *by_child.entry((probe.child_sg, probe.child_order)).or_insert(0) += 1;
        }
        println!("unsupported probes by (child space group, child co-group order):");
        for ((child, order), count) in &by_child {
            println!("  child #{child} order={order}: {count}");
        }
        for probe in unsupported.iter().take(5) {
            println!(
                "  witness ordinal {} {} t={} child #{}: {}",
                probe.ordinal, probe.label, probe.parameter, probe.child_sg, probe.detail
            );
        }
    }
    // The card-4 block statistics: totals by source, the legacy probe-level
    // reading on the same probes, and the two card-3 passes.
    let mut block_sources = [0usize; 3];
    let mut block_classes = [0usize; 3];
    let mut reference_blocks = 0usize;
    let mut blocks_with_gamma = 0usize;
    let mut empty_blocks = 0usize;
    for probe in probes {
        for block in &probe.blocks {
            block_sources[block.source.index()] += 1;
            if !block.cocycle_trivial {
                block_classes[block.source.index()] += 1;
            }
            if block.carries_reference {
                reference_blocks += 1;
            }
            if block.carries_gamma {
                blocks_with_gamma += 1;
            }
            if block.target_count == 0 {
                empty_blocks += 1;
            }
        }
    }
    let mut probe_classes = [0usize; 3];
    for probe in probes {
        probe_classes[match probe.class {
            TargetClass::Stored => 0,
            TargetClass::Constructed => 1,
            TargetClass::Unsupported | TargetClass::Error => 2,
        }] += 1;
    }
    println!(
        "child-star blocks: {} block(s) over {} probe(s), {} row(s) emitted for --output-blocks \
         ({} empty, {} carry the reference point, {} reach Gamma)",
        block_sources.iter().sum::<usize>(),
        probes.len(),
        evidence.block_rows,
        empty_blocks,
        reference_blocks,
        blocks_with_gamma
    );
    println!(
        "  block source: stored={} constructed={} mixed={}",
        block_sources[0], block_sources[1], block_sources[2]
    );
    println!(
        "  non-trivial own cocycle: stored={} constructed={} mixed={}",
        block_classes[0], block_classes[1], block_classes[2]
    );
    println!(
        "  withdrawn probe-level source on the same probes: stored={} constructed={} \
         other={}",
        probe_classes[0], probe_classes[1], probe_classes[2]
    );
    if recount.probes == 0 {
        println!("full-star recount: not run (pass --gate or --full-star-recount)");
    } else {
        println!(
            "full-star recount: {} added parameter probe(s), {} block(s), {} distinct block \
             geometr{}",
            recount.probes,
            recount.blocks,
            recount.geometries.len(),
            if recount.geometries.len() == 1 { "y" } else { "ies" }
        );
        println!(
            "  block source: stored={} constructed={} mixed={}",
            recount.sources[0], recount.sources[1], recount.sources[2]
        );
        println!(
            "  non-trivial own cocycle: stored={} constructed={} mixed={} by little co-group \
             order {:?}",
            recount.non_trivial[0],
            recount.non_trivial[1],
            recount.non_trivial[2],
            recount.non_trivial_orders
        );
        println!(
            "  (record, label) pairs whose recount geometry differs from the generic sample: {}",
            recount.changed_pairs.len()
        );
        println!(
            "Gamma-reaching arms: {} (parameter, arms) entr(ies), {} contained in the full-star \
             boundaries, {} exceptional (arm fixed exactly by the whole child point group), \
             {} always-Gamma arm(s)",
            gamma.entries, gamma.contained, gamma.exceptional, gamma.zero_arms
        );
        println!("  Gamma parameter-set shapes (pair count):");
        for (shape, count) in gamma.shapes.iter().take(12) {
            println!(
                "    {count:>4} pair(s): {{{}}}",
                shape
                    .iter()
                    .map(|(parameter, arms)| format!(
                        "t={parameter} arms[{}]",
                        arms.iter().map(usize::to_string).collect::<Vec<_>>().join(",")
                    ))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        if gamma.shapes.len() > 12 {
            println!("    ... {} more distinct shape(s)", gamma.shapes.len() - 12);
        }
        for witness in &gamma.witnesses {
            println!("    exception: {witness}");
        }
    }
    // R6.7 card 5: the interval sweep, printed whether or not it ran so a
    // vacuous pass is visible in the report as well as in the counters.
    if !evidence.sweep_ran {
        println!("domain sweep: not run (pass --domain-sweep)");
    } else {
        let sweep = evidence.sweep;
        println!(
            "domain sweep: {} (record, label) pair(s), {} interval(s), {} interior point(s), \
             {} comparison(s), {} failed interval(s)",
            sweep.pairs, sweep.intervals, sweep.points, sweep.comparisons, sweep.failed_intervals
        );
        println!(
            "  matched {} block pair(s) / {} target pair(s); interior block source \
             stored={} constructed={} mixed={}",
            sweep.blocks, sweep.targets, sweep.sources[0], sweep.sources[1], sweep.sources[2]
        );
        println!(
            "  interior targets per block {:?}, child little co-group order at the seed \
             points {:?}",
            sweep.target_counts, sweep.orders
        );
        println!(
            "  worst matched gauge-unified character score {:.17}, multiplicity mismatch(es) \
             {}, ambiguous matching(s) {}",
            sweep.worst_score, sweep.multiplicity_mismatches, sweep.ambiguities
        );
        println!(
            "  multi-target block shapes {:?}; blocks with a swappable same-dimension pair {}; \
             probe parameter(s) that are not full-star boundaries {}",
            sweep.multi_shapes, sweep.same_dimension_swaps, sweep.parent_only_parameters
        );
        println!(
            "  matched block pair(s) whose reduced representative differs between the two \
             interior points {} (of {}); exact-fixity check(s) of the card-1 gauge {} ({} \
             non-identity) with {} mismatch(es)",
            sweep.representative_shifts,
            sweep.blocks,
            sweep.little_fixity_checks,
            sweep.little_fixity_non_identity_checks,
            sweep.little_fixity_mismatches
        );
        // R6.7 card 6: the boundary pass and the independent arm-geometry audit.
        println!(
            "  boundary pass: {} boundary point(s) at the production engine, {} answered, {} \
             failed; answered decompositions {} parent / {} covered dimension",
            sweep.boundary_points,
            sweep.boundary_successes,
            sweep.boundary_failures,
            sweep.boundary_parent_dimension,
            sweep.boundary_covered_dimension
        );
        let geometry = &sweep.geometry;
        println!(
            "  independent arm-geometry audit: {} pair(s), {} relation(s) enumerated (closed form \
             {}), {} permanent (partition records {}), boundary set recomputed {} (partition {}), \
             {} disagreement(s)",
            geometry.pairs,
            geometry.relations,
            geometry.predicted_relations,
            geometry.permanent,
            geometry.permanent_partition,
            geometry.recomputed_boundaries,
            geometry.boundaries,
            geometry.disagreements
        );
        println!(
            "  relation predicate(s) evaluated: {} at the boundary parameters, {} at the {} \
             interior point(s), {} interior point(s) where a relation holds; orbit–stabiliser arm \
             count {} check(s) / {} disagreement(s) over {} arm(s)",
            geometry.evaluated,
            geometry.interior_evaluated,
            geometry.interior_points,
            geometry.interior_holds,
            geometry.arm_orbit_checks,
            geometry.arm_orbit_disagreements,
            geometry.arms
        );
    }
    if let Some(rows) = evidence.recount_file_rows {
        println!(
            "--output-recount: {} row(s) emitted, {rows} row(s) read back and parsed \
             against the collected fingerprint(s)",
            evidence.recount_rows_written
        );
    }
}

/// Compare the centring scan with the coordinate-map route for every source
/// rotation: an arithmetic-independent control of the step, since the two
/// algorithms share no code beyond the rational helpers.
fn step_algorithm_cross_check(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
) -> (usize, usize, Vec<String>) {
    let mut checks = 0usize;
    let mut mismatches = 0usize;
    let mut failures = Vec::new();
    for ((sg, label), domain) in domains {
        let Ok(lattice) = reciprocal_lattice(*sg) else {
            failures.push(format!("SG {sg} {label}: no reciprocal lattice"));
            continue;
        };
        let Ok(rotations) = rotation_set(*sg) else {
            failures.push(format!("SG {sg} {label}: no rotation set"));
            continue;
        };
        for rotation in &rotations {
            let image = match cryspglib::irrep::subduction::Mat3R::from_ints(*rotation)
                .inverse()
                .and_then(|matrix| matrix.transpose().checked_mul_vector(&domain.direction))
            {
                Ok(image) => image,
                Err(error) => {
                    failures.push(format!("SG {sg} {label}: rotation image: {error}"));
                    continue;
                }
            };
            let Ok(w) = image.checked_sub(&domain.direction) else {
                failures.push(format!("SG {sg} {label}: rotation difference failed"));
                continue;
            };
            match (
                minimal_parameter_step(&lattice, &w),
                minimal_parameter_step_via_coordinates(&lattice, &w),
            ) {
                (Ok(scan), Ok(coordinates)) => {
                    checks += 1;
                    if scan != coordinates {
                        mismatches += 1;
                        failures.push(format!(
                            "SG {sg} {label} rotation {rotation:?}: centring scan {:?} != \
                             coordinate map {:?}",
                            scan.map(|value| value.to_string()),
                            coordinates.map(|value| value.to_string())
                        ));
                    }
                }
                (scan, coordinates) => failures.push(format!(
                    "SG {sg} {label} rotation {rotation:?}: step failed ({:?} / {:?})",
                    scan.err().map(|error| error.to_string()),
                    coordinates.err().map(|error| error.to_string())
                )),
            }
        }
    }
    (checks, mismatches, failures)
}

/// The count conservations of the R6.7 card-6 sweep, as a gate-level check.
///
/// Split out of [`check_invariants`] so a fault-injection test can assert the
/// gate's own reaction to a report instead of re-implementing it: the two
/// conservations say that a parameter or an interval that was counted is either
/// answered or an explicitly counted failure, and the two zero-checks turn a
/// failure into a gate violation.
fn sweep_conservation(sweep: &SweepReport, violations: &mut Vec<String>) {
    if sweep.intervals != sweep.comparisons + sweep.failed_intervals {
        violations.push(format!(
            "the domain sweep checked {} interval(s) as {} comparison(s) + {} failure(s)",
            sweep.intervals, sweep.comparisons, sweep.failed_intervals
        ));
    }
    if sweep.boundary_points != sweep.boundary_successes + sweep.boundary_failures {
        violations.push(format!(
            "the domain sweep's boundary pass ran {} parameter(s) as {} success(es) + {} failure(s)",
            sweep.boundary_points, sweep.boundary_successes, sweep.boundary_failures
        ));
    }
    if sweep.failed_intervals > 0 {
        violations.push(format!(
            "the domain sweep failed on {} interval(s)",
            sweep.failed_intervals
        ));
    }
    if sweep.boundary_failures > 0 {
        violations.push(format!(
            "the domain sweep's boundary pass failed on {} parameter(s)",
            sweep.boundary_failures
        ));
    }
}

/// Gate-level synthetic self-check of the audit's own comparisons.
///
/// The corpus is **clean**: every recomputed boundary set equals the partition's,
/// every recorded event count equals the recomputed one and no interior point
/// carries a relation.  A comparison branch that is deleted or that always
/// answers "disagreement" is therefore invisible on the corpus, which is exactly
/// the failure mode the previous review rounds found (a control that is only
/// ever exercised in its passing direction).  Two synthetic partitions -- the
/// fields are public -- drive the same [`audit_arm_geometry`] code:
///
/// * the honest one, whose boundary set and event counts are derived from the
///   same relations: the audit must report neither a set nor a count
///   disagreement;
/// * the same partition with one recorded count lowered to `[0, 1, 0]`: the
///   count comparison must fire;
/// * the same partition with its second boundary moved from `1/2` to `1/4`: the
///   set comparison must fire.
///
/// The synthetic arms are `(1, 0, 0)` and `(-1, 0, 0)` with the identity
/// rotation, so the two arm-merge relations `t . (+-2, 0, 0) in Z^3` hold exactly
/// at `t = 0` and `t = 1/2`.  Only the two messages under test are asserted: the
/// arm count of this fixture is deliberately not the orbit of any real source, so
/// the orbit-stabiliser check reports its own (expected) disagreement and is not
/// part of this self-check.
fn audit_self_check() -> (usize, Vec<String>) {
    let mut violations = Vec::new();
    // How many synthetic partitions were actually driven through the audit.  The
    // gate pins it, so a self-check reduced to its honest fixture (the one that
    // cannot fail) is visible as a moved count instead of as a green gate.
    let mut cases = 0usize;
    let arms = vec![
        FoldedArm {
            direction: Vec3R::from_ints([1, 0, 0]),
            parent_rotation: IDENTITY_ROTATION,
        },
        FoldedArm {
            direction: Vec3R::from_ints([-1, 0, 0]),
            parent_rotation: IDENTITY_ROTATION,
        },
    ];
    let identity: Mat3I = IDENTITY_ROTATION;
    let base = FullStarPartition {
        source_sg: 196,
        child_sg: 1,
        label: "SELFTEST",
        arms,
        reference_arm: 0,
        child_reciprocal: Lattice::integer(),
        child_rotations: vec![identity],
        boundaries: Vec::new(),
        permanent_counts: [0; 3],
        permanent_witnesses: Vec::new(),
    };
    let zero = Rat::ZERO;
    let half = Rat::new(1, 2).expect("1/2");
    let boundary = |parameter: Rat, counts: [usize; 3]| StarBoundary {
        parameter,
        counts,
        witnesses: Vec::new(),
    };
    let mut audit = |partition: &FullStarPartition| -> Vec<String> {
        cases += 1;
        let table = line_table(196, "DT1").expect("the frozen DT1 table");
        let mut report = SweepReport {
            worst_score: 1.0,
            ..SweepReport::default()
        };
        audit_arm_geometry(table, partition, 0, "the audit self-check", &mut report);
        report.failures
    };
    let reports = |failures: &[String], needle: &str| -> bool {
        failures.iter().any(|failure| failure.contains(needle))
    };
    // The honest fixture: two arm-merge relations holding at 0 and 1/2.
    let mut honest = base.clone();
    honest.boundaries = vec![
        boundary(zero, [0, 2, 0]),
        boundary(half, [0, 2, 0]),
    ];
    let failures = audit(&honest);
    for needle in ["boundary set is not the partition's", "the audit finds"] {
        if reports(&failures, needle) {
            violations.push(format!(
                "the arm-geometry audit's self-check reported {needle:?} on an honest synthetic \
                 partition: {failures:?}"
            ));
        }
    }
    // The recorded count lowered: the comparison has to fire.
    let mut wrong_count = honest.clone();
    wrong_count.boundaries[1].counts = [0, 1, 0];
    let failures = audit(&wrong_count);
    if !reports(&failures, "but the partition records 0/1/0") {
        violations.push(format!(
            "the arm-geometry audit's self-check did not notice a recorded event count of 0/1/0 \
             where the relations hold 0/2/0: {failures:?}"
        ));
    }
    // The boundary moved: the set comparison has to fire.
    let mut moved = honest.clone();
    moved.boundaries[1] = boundary(Rat::new(1, 4).expect("1/4"), [0, 0, 0]);
    let failures = audit(&moved);
    if !reports(&failures, "does not cut [1/2]") {
        violations.push(format!(
            "the arm-geometry audit's self-check did not notice a boundary moved from 1/2 to 1/4: \
             {failures:?}"
        ));
    }
    (cases, violations)
}

/// Every invariant the census claims, checked against the probes.
fn check_invariants(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    records: &[Record],
    probes: &[Probe],
    evidence: &CensusEvidence,
    violations: &mut Vec<String>,
) {
    // The child-side step comparison runs once per (record, label) per child
    // rotation, so its total is `sum over (record, label) of |rotation_set(child)|`
    // — with the rotation set counted through the crate's **independent** database
    // route (the stored Hall setting, as in
    // `every_space_group_rotation_set_is_its_stored_hall_settings`).  Without this
    // the total is only printed: dropping one rotation from `rotation_set` leaves
    // the gate's summary line intact and merely lowers the printed count
    // (measured: 28,713 -> 28,666 when SG 38 loses one).
    let mut expected_child_comparisons = 0usize;
    for record in records {
        let stored = SG_DATA_HALL[usize::from(record.child_sg)];
        let independent = match SymmetryOps::from_hall_number(
            HallNumber::try_from(usize::from(stored)).expect("stored Hall number"),
        ) {
            Ok(operations) => operations,
            Err(error) => {
                violations.push(format!(
                    "child #{}: independent rotation lookup failed: {error}",
                    record.child_sg
                ));
                continue;
            }
        };
        let mut rotations: Vec<[[i32; 3]; 3]> = Vec::new();
        for operation in &independent.operations {
            if !rotations.contains(&operation.rotation) {
                rotations.push(operation.rotation);
            }
        }
        expected_child_comparisons = expected_child_comparisons
            .saturating_add(rotations.len().saturating_mul(record.labels.len()));
    }
    if evidence.child_algorithms.0 != expected_child_comparisons {
        violations.push(format!(
            "the child step comparison ran {} time(s) but the stored Hall settings of the \
             records' child space groups predict {expected_child_comparisons}",
            evidence.child_algorithms.0
        ));
    }
    // 1. The frozen tables are the generic stabiliser (already enforced while
    //    building the domains; repeated here so a future refactor cannot drop it).
    for ((sg, label), domain) in domains {
        let Some(table) = line_table(*sg, label) else {
            violations.push(format!("SG {sg} {label}: no frozen table"));
            continue;
        };
        if domain.generic_order != table.operations.len() {
            violations.push(format!(
                "SG {sg} {label}: generic order {} != frozen operations {}",
                domain.generic_order,
                table.operations.len()
            ));
        }
        if domain.exceptional.is_empty() {
            violations.push(format!("SG {sg} {label}: no exceptional parameter at all"));
        }
        if !domain
            .exceptional
            .iter()
            .any(|entry| entry.parameter.is_zero())
        {
            violations.push(format!("SG {sg} {label}: t = 0 is not in the domain"));
        }
        for entry in &domain.exceptional {
            if entry.order < domain.generic_order {
                violations.push(format!(
                    "SG {sg} {label}: parameter {} has order {} below the generic {}",
                    entry.parameter, entry.order, domain.generic_order
                ));
            }
        }
    }
    // 2. The production `ParameterKind` agrees with the census's notion of an
    //    enhanced parameter.
    for probe in probes {
        let Some(kind) = probe.parameter_kind else {
            continue;
        };
        let enhanced = domains
            .get(&(probe.parent_sg, probe.label))
            .is_some_and(|domain| domain.is_exceptional(&probe.parameter));
        let formal = kind == ParameterKind::Formal;
        if formal != enhanced {
            violations.push(format!(
                "ordinal {} {} t={}: ParameterKind {:?} but census enhanced={enhanced}",
                probe.ordinal, probe.label, probe.parameter, kind
            ));
        }
    }
    // 3. Every probe carries a real child little co-group order (never the
    //    sentinel zero), a trivial co-group is always answerable, and an
    //    unsupported probe must be explained by an enhanced child co-group.
    for probe in probes {
        if probe.child_order == 0 {
            violations.push(format!(
                "ordinal {} {} t={}: no child little co-group order was computed",
                probe.ordinal, probe.label, probe.parameter
            ));
        }
        if probe.class == TargetClass::Unsupported && probe.child_order <= 1 {
            violations.push(format!(
                "ordinal {} {} t={}: unsupported although the child co-group has order {}",
                probe.ordinal, probe.label, probe.parameter, probe.child_order
            ));
        }
    }
    // 3b. An engine error is a hard failure of the census run, under every flag:
    //     `errors` is not allowed to be a pure diagnostic counter.
    for probe in probes.iter().filter(|probe| probe.class == TargetClass::Error) {
        violations.push(format!(
            "ordinal {} {} t={}: engine error: {}",
            probe.ordinal, probe.label, probe.parameter, probe.detail
        ));
    }
    // 3c. The child-side census is cross-checked against the child group on a
    //     uniform grid, and its corpus-wide candidate set is exactly the eighth
    //     grid measured for this corpus.  Both are assertions, not report lines.
    if evidence.child_grid.1 > 0 {
        violations.push(format!(
            "the child enumeration disagrees with the child group in {} of {} grid predicates",
            evidence.child_grid.1, evidence.child_grid.0
        ));
    }
    if evidence.child_grid.0 == 0 {
        violations.push("the child grid cross-check never ran".to_string());
    }
    let expected_child: Vec<Rat> = [
        (0i128, 1i128),
        (1, 8),
        (1, 4),
        (3, 8),
        (1, 2),
        (5, 8),
        (3, 4),
        (7, 8),
    ]
    .iter()
    .map(|(numerator, denominator)| Rat::new(*numerator, *denominator).expect("rational"))
    .collect();
    if evidence.child_union != expected_child.as_slice() {
        violations.push(format!(
            "the child candidate parameters are {:?}, expected the eighth grid {:?}",
            evidence
                .child_union
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>(),
            expected_child
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
        ));
    }
    // The Gamma-reaching condition, asserted rather than only printed (card-4
    // audit F5): the mutation that dropped the `t = 0` entries printed a
    // self-contradictory line ("17,268 entries, 23,024 contained") and still left
    // the gate green.  The four counts and the shape total are corpus
    // measurements; `contained == entries` is the measured claim that every
    // Gamma-reaching parameter is a full-star boundary (the documented exception
    // has no corpus witness).
    if evidence.gamma.entries != 23_024 {
        violations.push(format!(
            "the Gamma enumeration reports {} (parameter, arms) entr(ies), expected 23024",
            evidence.gamma.entries
        ));
    }
    if evidence.gamma.contained != evidence.gamma.entries {
        violations.push(format!(
            "the Gamma enumeration reports {} contained of {} entr(ies); every entry must be a \
             full-star boundary",
            evidence.gamma.contained, evidence.gamma.entries
        ));
    }
    if evidence.gamma.exceptional != 0 {
        violations.push(format!(
            "{} Gamma entr(ies) took the exactly-fixed-arm exception, which has no corpus witness",
            evidence.gamma.exceptional
        ));
    }
    if evidence.gamma.zero_arms != 0 {
        violations.push(format!(
            "{} arm(s) have a zero folded direction and are at Gamma everywhere",
            evidence.gamma.zero_arms
        ));
    }
    let gamma_pairs: usize = evidence.gamma.shapes.values().sum();
    if gamma_pairs != 5_756 {
        violations.push(format!(
            "the Gamma parameter-set shapes cover {gamma_pairs} (record, label) pairs, expected 5756"
        ));
    }
    if evidence.gamma.shapes.len() != 8 {
        violations.push(format!(
            "the Gamma parameter-set shapes number {}, expected 8",
            evidence.gamma.shapes.len()
        ));
    }
    // The full-star recount pass, when it was requested, must have run and must
    // reproduce its corpus totals -- otherwise a mutation can turn it into a no-op
    // and the gate keeps passing (verification review F6).
    if evidence.recount_ran {
        if evidence.recount.probes != 11_009 {
            violations.push(format!(
                "the full-star recount probed {} parameter(s), expected 11009",
                evidence.recount.probes
            ));
        }
        if evidence.recount.blocks != 27_507 {
            violations.push(format!(
                "the full-star recount classified {} block(s), expected 27507",
                evidence.recount.blocks
            ));
        }
        if evidence.recount.geometries.len() != 28 {
            violations.push(format!(
                "the full-star recount saw {} distinct block geometr(ies), expected 28",
                evidence.recount.geometries.len()
            ));
        }
        if evidence.recount.changed_pairs.len() != 1_692 {
            violations.push(format!(
                "the full-star recount changed {} (record, label) pair(s) against t = 1/7, \
                 expected 1692",
                evidence.recount.changed_pairs.len()
            ));
        }
        // The card-4 residual: the recount now has a per-probe artifact, and its
        // row count is the pass's own probe count -- not a separately written
        // number that could drift.  The rows are also parsed back and compared
        // field by field (in `run`), so a corrupted row cannot pass on the count.
        if evidence.recount.rows.len() != evidence.recount.probes {
            violations.push(format!(
                "the recount fingerprint carries {} row(s) for {} probe(s)",
                evidence.recount.rows.len(),
                evidence.recount.probes
            ));
        }
        if let Some(rows) = evidence.recount_file_rows
            && rows != evidence.recount.probes
        {
            violations.push(format!(
                "the --output-recount file has {rows} row(s) for {} probe(s)",
                evidence.recount.probes
            ));
        }
        // The emission-loop counter is only meaningful when the file was asked
        // for; without the flag it is zero by construction (the check below runs
        // on the file's own row count instead).
        if evidence.recount_file_rows.is_some()
            && evidence.recount.rows.len() != evidence.recount_rows_written
        {
            violations.push(format!(
                "the --output-recount emission loop wrote {} row(s) for {} collected fingerprint(s)",
                evidence.recount_rows_written,
                evidence.recount.rows.len()
            ));
        }
    }
    // R6.7 card 5/6: the interval sweep and the card-6 boundary pass.  Count
    // conservation comes first: every interval is either compared or an
    // explicitly reported failure, and every probe parameter either reached the
    // engine or is a counted failure, so a lost item cannot hide behind the
    // failure list or behind a green gate.
    if evidence.sweep_ran {
        sweep_conservation(evidence.sweep, violations);
        // The audit's own comparison branches, driven synthetically: the corpus
        // only ever exercises their passing direction.  The number of fixtures
        // actually audited is pinned, so a self-check that stops driving the two
        // failing fixtures is visible.
        let (self_check_cases, self_check_violations) = audit_self_check();
        violations.extend(self_check_violations);
        if self_check_cases != SELF_CHECK_CASE_PIN {
            violations.push(format!(
                "the arm-geometry audit's self-check drove {self_check_cases} synthetic \
                 partition(s), expected {SELF_CHECK_CASE_PIN}"
            ));
        }
        if evidence.sweep.multiplicity_mismatches > 0 {
            violations.push(format!(
                "the domain sweep matched {} target pair(s) whose multiplicities disagree",
                evidence.sweep.multiplicity_mismatches
            ));
        }
        if evidence.sweep.ambiguities > 0 {
            violations.push(format!(
                "the domain sweep left {} matching(s) ambiguous or incomplete",
                evidence.sweep.ambiguities
            ));
        }
        if evidence.sweep.worst_score < 1.0 - 1e-9 {
            violations.push(format!(
                "the domain sweep's worst matched gauge-unified character score is {}",
                evidence.sweep.worst_score
            ));
        }
        // The corpus totals, pinned: a sweep that silently stops sweep
        // (a wrong interval set, a shortened parameter list, a dropped block)
        // moves one of them.  All of them are measured on the same run that
        // establishes the pass.  The card-6 boundary pass and geometry audit add
        // their own pinned totals in the block right below.
        for (name, found, expected) in [
            ("(record, label) pair(s)", evidence.sweep.pairs, 5_756usize),
            ("interval(s)", evidence.sweep.intervals, 46_048),
            ("interior point(s)", evidence.sweep.points, 92_096),
            ("comparison(s)", evidence.sweep.comparisons, 46_048),
            ("matched block pair(s)", evidence.sweep.blocks, 168_408),
            ("matched target pair(s)", evidence.sweep.targets, 174_672),
            ("boundary point(s)", evidence.sweep.boundary_points, BOUNDARY_POINT_PIN),
            (
                "answered boundary point(s)",
                evidence.sweep.boundary_successes,
                BOUNDARY_POINT_PIN,
            ),
            ("failed boundary point(s)", evidence.sweep.boundary_failures, 0),
        ] {
            if found != expected {
                violations.push(format!(
                    "the domain sweep covered {found} {name}, expected {expected}"
                ));
            }
        }
        // The two dimension sums are `u64` (the engine's own `u32` dimensions
        // accumulated over 46,048 parameters), so they are compared and pinned
        // in their own loop rather than through a lossy conversion.
        for (name, found, expected) in [
            (
                "boundary parent dimension",
                evidence.sweep.boundary_parent_dimension,
                BOUNDARY_DIMENSION_PIN,
            ),
            (
                "boundary covered dimension",
                evidence.sweep.boundary_covered_dimension,
                BOUNDARY_DIMENSION_PIN,
            ),
        ] {
            if found != expected {
                violations.push(format!(
                    "the domain sweep covered {found} {name}, expected {expected}"
                ));
            }
        }
        // The two dimension sums of the boundary pass are the engine's own
        // bookkeeping on every answered parameter; they must agree, which is the
        // only census-side trace of the returned decompositions themselves.  The
        // parent sum is then compared with the closed form the audit accumulates
        // from the partition's own counts, so a pass that stops asking the engine
        // cannot keep the pin.
        if evidence.sweep.boundary_parent_dimension
            != evidence.sweep.geometry.predicted_boundary_dimension
        {
            violations.push(format!(
                "the boundary pass reports {} parent dimension over its probe parameters, the closed \
                 form (parameters x little dimension x arms) predicts {}",
                evidence.sweep.boundary_parent_dimension,
                evidence.sweep.geometry.predicted_boundary_dimension
            ));
        }
        if evidence.sweep.boundary_parent_dimension != evidence.sweep.boundary_covered_dimension {
            violations.push(format!(
                "the boundary pass' answered decompositions cover {} of the {} parent dimension \
                 they report",
                evidence.sweep.boundary_covered_dimension, evidence.sweep.boundary_parent_dimension
            ));
        }
        // R6.7 card 6: the independent arm-geometry audit.  Its own conservation
        // relations come first -- the enumerated relation count against the closed
        // form read from the partition's arm and rotation counts, and the
        // identically-true relations against the partition's own accounting -- and
        // then the pinned totals, each of which has an independent predicted value
        // next to it.
        let geometry = &evidence.sweep.geometry;
        if geometry.relations != geometry.predicted_relations {
            violations.push(format!(
                "the arm-geometry audit enumerated {} relation(s), the closed form n^2 r - n of the \
                 partition's own counts predicts {}",
                geometry.relations, geometry.predicted_relations
            ));
        }
        if geometry.permanent != geometry.permanent_partition {
            violations.push(format!(
                "the arm-geometry audit finds {} identically-true relation(s), the partition \
                 records {}",
                geometry.permanent, geometry.permanent_partition
            ));
        }
        if geometry.recomputed_boundaries != geometry.boundaries {
            violations.push(format!(
                "the arm-geometry audit recomputed {} boundary parameter(s), the partition claims \
                 {}",
                geometry.recomputed_boundaries, geometry.boundaries
            ));
        }
        match evidence
            .sweep
            .boundary_points
            .checked_sub(evidence.sweep.parent_only_parameters)
        {
            Some(expected) if expected == geometry.boundaries => {}
            Some(expected) => violations.push(format!(
                "the arm-geometry audit checked {} full-star boundary parameter(s), but the sweep \
                 probed {expected} parameter(s) that are not parent-only",
                geometry.boundaries
            )),
            None => violations.push(
                "the domain sweep counted more parent-only parameters than probe parameters"
                    .to_string(),
            ),
        }
        if geometry.interior_points != 2 * evidence.sweep.intervals {
            violations.push(format!(
                "the arm-geometry audit evaluated {} interior point(s) for {} interval(s), expected \
                 two per interval",
                geometry.interior_points, evidence.sweep.intervals
            ));
        }
        if geometry.interior_holds > 0 {
            violations.push(format!(
                "the arm-geometry audit found a full-star relation holding at {} interior point(s) \
                 of the sweep",
                geometry.interior_holds
            ));
        }
        if geometry.disagreements > 0 {
            violations.push(format!(
                "the arm-geometry audit disagreed with the partition {} time(s)",
                geometry.disagreements
            ));
        }
        if geometry.evaluated == 0 || geometry.interior_evaluated == 0 {
            violations.push(format!(
                "the arm-geometry audit evaluated {} boundary and {} interior relation predicate(s)",
                geometry.evaluated, geometry.interior_evaluated
            ));
        }
        // The evaluation totals are bound to the relation list rather than only
        // pinned: every boundary parameter evaluates its pair's whole
        // non-permanent relation list, and every interval evaluates the same list
        // at both of its interior points, so
        // `evaluated x pairs == boundaries x (relations - permanent)` and
        // `interior_evaluated x pairs == 2 x intervals x (relations - permanent)`.
        // A loop that stops early, skips a relation or evaluates a stale list
        // moves one side and not the other.  (On this corpus both reduce to
        // 8x and 16x the non-permanent relation count.)
        let non_permanent = match geometry.relations.checked_sub(geometry.permanent) {
            Some(value) => value,
            None => {
                violations.push(format!(
                    "the arm-geometry audit counts {} identically-true relation(s) more than its {} \
                     enumerated relation(s)",
                    geometry.permanent, geometry.relations
                ));
                0
            }
        };
        for (name, left, right) in [
            (
                "boundary",
                (geometry.evaluated as u128) * (geometry.pairs as u128),
                (geometry.boundaries as u128) * (non_permanent as u128),
            ),
            (
                "interior",
                (geometry.interior_evaluated as u128) * (geometry.pairs as u128),
                (2 * evidence.sweep.intervals as u128) * (non_permanent as u128),
            ),
        ] {
            if left != right {
                violations.push(format!(
                    "the arm-geometry audit's {name} evaluation count is {left} relation-pair \
                     evaluation(s), but every parameter evaluates the whole non-permanent relation \
                     list, which gives {right}"
                ));
            }
        }
        if geometry.arm_orbit_checks != evidence.sweep.pairs {
            violations.push(format!(
                "the orbit–stabiliser arm count ran {} time(s) for {} pair(s)",
                geometry.arm_orbit_checks, evidence.sweep.pairs
            ));
        }
        if geometry.arm_orbit_disagreements > 0 {
            violations.push(format!(
                "the orbit–stabiliser arm count disagreed with the partition {} time(s)",
                geometry.arm_orbit_disagreements
            ));
        }
        for (name, found, expected) in [
            (
                "(record, label) pair(s) with a recomputed relation list",
                geometry.pairs,
                5_756usize,
            ),
            ("relation(s) enumerated", geometry.relations, RELATION_PIN),
            (
                "identically-true relation(s)",
                geometry.permanent,
                PERMANENT_RELATION_PIN,
            ),
            (
                "boundary parameter(s) audited",
                geometry.boundaries,
                BOUNDARY_POINT_PIN,
            ),
            (
                "boundary relation predicate(s) evaluated",
                geometry.evaluated,
                BOUNDARY_RELATION_EVALUATION_PIN,
            ),
            (
                "interior relation predicate(s) evaluated",
                geometry.interior_evaluated,
                INTERIOR_RELATION_EVALUATION_PIN,
            ),
            (
                "arm(s) counted by the orbit–stabiliser route",
                geometry.arms,
                ARM_ORBIT_PIN,
            ),
            (
                "orbit–stabiliser arm count check(s)",
                geometry.arm_orbit_checks,
                5_756,
            ),
        ] {
            if found != expected {
                violations.push(format!(
                    "the arm-geometry audit measured {found} {name}, expected {expected}"
                ));
            }
        }
        // Measured corpus facts, not assumptions: no interior point of any
        // interval has a stored or a mixed target (so the stored/constructed
        // match is a synthetic branch on this corpus), the little co-group at
        // the seed points is order 1, 2, 4 or 8, and every multi-target block
        // carries one target per dimension.
        if evidence.sweep.sources != [0, 174_672, 0] {
            violations.push(format!(
                "the domain sweep's interior blocks are stored/constructed/mixed {}/{}/{}, \
                 expected 0/174672/0",
                evidence.sweep.sources[0], evidence.sweep.sources[1], evidence.sweep.sources[2]
            ));
        }
        let expected_target_counts: BTreeMap<usize, usize> =
            [(1usize, 162_144usize), (2, 6_264)].into_iter().collect();
        if evidence.sweep.target_counts != expected_target_counts {
            violations.push(format!(
                "the domain sweep's interior targets per block are {:?}, expected {:?}",
                evidence.sweep.target_counts, expected_target_counts
            ));
        }
        let expected_orders: BTreeMap<usize, usize> = [
            (1usize, 91_984usize),
            (2, 62_408),
            (4, 13_584),
            (8, 432),
        ]
        .into_iter()
        .collect();
        if evidence.sweep.orders != expected_orders {
            violations.push(format!(
                "the domain sweep's seed-point little co-group orders are {:?}, expected {:?}",
                evidence.sweep.orders, expected_orders
            ));
        }
        // Pinned measured facts (card-5 audit F3): every multi-target interior
        // block has the shape [(1,1),(1,1)], and no probe parameter comes from
        // the parent side alone (the parent's formal boundaries are nested in
        // the eighth grid on this corpus; the set equality above is what keeps
        // this from hiding a moved child boundary).
        let expected_shapes: BTreeMap<Vec<(u8, u32)>, usize> =
            [(vec![(1u8, 1u32), (1, 1)], 6_264usize)].into_iter().collect();
        if evidence.sweep.multi_shapes != expected_shapes {
            violations.push(format!(
                "the domain sweep's multi-target block shapes are {:?}, expected {:?}",
                evidence.sweep.multi_shapes, expected_shapes
            ));
        }
        if evidence.sweep.parent_only_parameters != 0 {
            violations.push(format!(
                "the domain sweep has {} probe parameter(s) that are not full-star boundaries,                  expected 0",
                evidence.sweep.parent_only_parameters
            ));
        }
        if evidence.sweep.same_dimension_swaps != 0 {
            violations.push(format!(
                "the domain sweep found {} interior block(s) with two same-dimension targets of \
                 different multiplicities, which the pinned corpus shape says do not exist",
                evidence.sweep.same_dimension_swaps
            ));
        }
        // The card-1 gauge's hypothesis has to be exercised and to hold: every
        // aligned little-group rotation fixes its seed point exactly (a
        // non-exact fixity would be a partition boundary, so this is also a
        // partition control).
        if evidence.sweep.little_fixity_mismatches > 0 {
            violations.push(format!(
                "the exact-fixity control of the card-1 gauge failed {} time(s)",
                evidence.sweep.little_fixity_mismatches
            ));
        }
        // Non-vacuousness has to hold **at gate level**, not only in the unit
        // test: the count is predicted from the independently measured little
        // co-group histogram as `2 x sum(order x blocks)` (one check per aligned
        // operation per side), so a control that stops checking operations is
        // caught here, and a synthetic pair separates exact from mod-`L*` fixity,
        // so one that always answers `true` is caught too (audit of `d66db48`,
        // P0-2: reducing the loop to the identity operation kept the gate green,
        // and `exactly_fixes` constant was invisible to it).
        let predicted_fixity: usize = 2 * evidence
            .sweep
            .orders
            .iter()
            .map(|(order, count)| order * count)
            .sum::<usize>();
        if evidence.sweep.little_fixity_checks != predicted_fixity {
            violations.push(format!(
                "the exact-fixity control of the card-1 gauge ran {} time(s), predicted \
                 2 x sum(little-co-group order x block pairs) = {}",
                evidence.sweep.little_fixity_checks, predicted_fixity
            ));
        }
        // The **operand** has to be pinned too: each aligned little co-group has
        // exactly one identity, so the number of checks whose rotation is not the
        // identity is `2 x (sum(order x blocks) - block pairs)`.  Replacing the
        // checked rotation by the identity keeps the loop and the count green
        // (verification of `b2d39f3`, P0) and moves this counter.
        // `2 * sum - 2 * blocks`: written out because `2 * (sum).saturating_sub(..)`
        // parses the other way round and silently predicted 0.
        let predicted_non_identity = 2 * evidence
            .sweep
            .orders
            .iter()
            .map(|(order, count)| order * count)
            .sum::<usize>()
            - 2 * evidence.sweep.blocks;
        if evidence.sweep.little_fixity_non_identity_checks != predicted_non_identity {
            violations.push(format!(
                "the exact-fixity control checked {} non-identity rotation(s), predicted \
                 2 x (sum(order x blocks) - block pairs) = {}",
                evidence.sweep.little_fixity_non_identity_checks, predicted_non_identity
            ));
        }
        if predicted_fixity == 0 {
            violations.push(
                "the exact-fixity control of the card-1 gauge never ran (no aligned little-group \
                 operation was checked)"
                    .to_string(),
            );
        }
        let selftest_rotation: Mat3I = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let half = Rat::new(1, 2).expect("the constant 1/2");
        let third = Rat::new(1, 3).expect("the constant 1/3");
        let selftest_lattice = reciprocal_lattice(1).expect("SG 1 reciprocal lattice");
        for (label, point, expected) in [
            (
                "fixed only modulo L*",
                Vec3R::new([half, half, Rat::ZERO]),
                false,
            ),
            ("a point on the axis", Vec3R::new([Rat::ZERO, Rat::ZERO, third]), true),
        ] {
            // The input is bound independently of the answer: it must be fixed
            // modulo `L*` (the little group's own predicate), otherwise this pair
            // of cases could be replaced by two copies of the other point.
            match selftest_lattice.preserves(selftest_rotation, &point) {
                Ok(true) => {}
                Ok(false) => violations.push(format!(
                    "the exact-fixity control's synthetic witness ({label}) is not fixed modulo \
                     the reciprocal lattice, so it does not separate exact from mod-L* fixity"
                )),
                Err(error) => violations.push(format!(
                    "the exact-fixity control's synthetic witness ({label}) lattice test: {error}"
                )),
            }
            match exactly_fixes(&selftest_rotation, &point) {
                Ok(found) if found == expected => {}
                Ok(found) => violations.push(format!(
                    "the exact-fixity control's own synthetic witness ({label}) answers {found}, \
                     expected {expected}"
                )),
                Err(error) => violations.push(format!(
                    "the exact-fixity control's own synthetic witness ({label}): {error}"
                )),
            }
        }
        // Measured, not assumed: the engine's reduced representative does move
        // between the two interior points of some intervals (a reduction wall,
        // R6.7 card 5 math review P1-2), while the gauge-unified content the
        // sweep compares stays constant.
        if evidence.sweep.representative_shifts != REPRESENTATIVE_SHIFT_PIN {
            violations.push(format!(
                "the domain sweep saw {} matched block pair(s) whose reduced representative \
                 differs between the two interior points, expected {} (every matched pair)",
                evidence.sweep.representative_shifts, REPRESENTATIVE_SHIFT_PIN
            ));
        }
    }
    if evidence.algorithm_mismatches > 0 {
        violations.push(format!(
            "the centring scan and the coordinate map disagree in {} of {} rotations",
            evidence.algorithm_mismatches, evidence.algorithm_checks
        ));
    }
    if evidence.algorithm_checks < 2_000 {
        violations.push(format!(
            "the two step algorithms were compared only {} times",
            evidence.algorithm_checks
        ));
    }
    if evidence.child_algorithms.1 > 0 {
        violations.push(format!(
            "the child step algorithms disagree in {} of {} comparisons",
            evidence.child_algorithms.1, evidence.child_algorithms.0
        ));
    }
    if evidence.child_algorithms.0 == 0 {
        violations.push("the child two-algorithm comparison never ran".to_string());
    }
    // 3e. The frozen directions are in scope for the residue representation.
    for ((sg, label), domain) in domains {
        match reciprocal_lattice(*sg).and_then(|lattice| {
            require_reciprocal_direction(&lattice, &domain.direction, *sg, label)
        }) {
            Ok(()) => {}
            Err(error) => violations.push(format!("SG {sg} {label}: {error}")),
        }
    }
    // 4. The official anchor is answered for every pinned row, and its trivial
    //    content equals the pinned frequency (compared while probing; the label
    //    detail carries the mismatch).
    //
    // `official_line_parameter().ok()` deliberately maps the error to `None`:
    // with it, no probe matches the anchor below and `anchor == 0` is itself the
    // violation two checks further down, so the failure is counted either way.
    let official = official_line_parameter().ok();
    let mut anchor = 0usize;
    for probe in probes {
        if Some(probe.parameter) != official {
            continue;
        }
        anchor += 1;
        if probe.class == TargetClass::Unsupported {
            violations.push(format!(
                "ordinal {} {}: the official anchor is unsupported",
                probe.ordinal, probe.label
            ));
        }
        if probe.detail.contains("!= pinned") {
            violations.push(format!(
                "ordinal {} {}: {detail}",
                probe.ordinal,
                probe.label,
                detail = probe.detail
            ));
        }
    }
    if anchor != 5_756 {
        violations.push(format!("the anchor probes {anchor} rows, expected 5756"));
    }
    // 5. Every probe that the engine answered carries a trivial content.
    for probe in probes {
        if probe.class != TargetClass::Error
            && probe.class != TargetClass::Unsupported
            && probe.content.is_none()
        {
            violations.push(format!(
                "ordinal {} {} t={}: answered without a trivial content",
                probe.ordinal, probe.label, probe.parameter
            ));
        }
    }
    // 6. No probe was dropped: every record contributed at least one.
    if records.is_empty() {
        violations.push("no isotropy record carries parametric-k rows".to_string());
    }
    // 7b. The frozen table is the generic stabiliser as a *set*, not only in
    //     order: compare its rotations with those that fix the direction exactly.
    for ((sg, label), domain) in domains {
        let Ok(rotations) = rotation_set(*sg) else {
            violations.push(format!("SG {sg} {label}: no rotation set"));
            continue;
        };
        let mut generic_flat: Vec<Vec<i32>> = Vec::new();
        for rotation in &rotations {
            let Ok(image) = cryspglib::irrep::subduction::Mat3R::from_ints(*rotation)
                .inverse()
                .and_then(|matrix| matrix.transpose().checked_mul_vector(&domain.direction))
            else {
                violations.push(format!("SG {sg} {label}: rotation image failed"));
                continue;
            };
            if image == domain.direction {
                generic_flat.push(rotation.iter().flatten().copied().collect());
            }
        }
        let Some(table) = line_table(*sg, label) else {
            violations.push(format!("SG {sg} {label}: no frozen table"));
            continue;
        };
        let mut frozen: Vec<Vec<i32>> = table
            .operations
            .iter()
            .map(|operation| {
                operation
                    .rotation
                    .iter()
                    .flatten()
                    .map(|value| i32::from(*value))
                    .collect()
            })
            .collect();
        generic_flat.sort();
        frozen.sort();
        if generic_flat != frozen {
            violations.push(format!(
                "SG {sg} {label}: the frozen rotation set is not the generic stabiliser \
                 (generic {} vs frozen {} rotations)",
                generic_flat.len(),
                frozen.len()
            ));
        }
    }
    // 8. The exact enumeration is cross-checked against the group itself on a
    //    uniform grid.  The two methods share only the lattice membership test.
    let mut grid_checks = 0usize;
    let mut grid_mismatches = 0usize;
    for ((sg, label), domain) in domains {
        let Ok(lattice) = reciprocal_lattice(*sg) else {
            violations.push(format!("SG {sg} {label}: no parent reciprocal lattice"));
            continue;
        };
        let Ok(rotations) = rotation_set(*sg) else {
            violations.push(format!("SG {sg} {label}: no rotation set"));
            continue;
        };
        for denominator in [24i128, 120] {
            match verify_against_grid(&lattice, &domain.direction, &rotations, denominator) {
                Ok(check) => {
                    grid_checks += check.checks;
                    grid_mismatches += check.mismatches;
                }
                Err(error) => violations.push(format!(
                    "SG {sg} {label}: grid check 1/{denominator}: {error}"
                )),
            }
        }
    }
    if grid_mismatches > 0 {
        violations.push(format!(
            "the enumeration disagrees with the group in {grid_mismatches} of \
             {grid_checks} grid predicates"
        ));
    }
    if grid_checks == 0 {
        violations.push("the parent grid cross-check never ran".to_string());
    }
    println!("grid cross-check: {grid_checks} predicates, {grid_mismatches} mismatch(es)");
    // 7. The generic sample is never a formal parameter and never unsupported.
    //    A generic sample that cannot even be built is a violation, not a
    //    silently skipped check (card 6: the `if let Ok` here used to be the
    //    only thing standing between a broken constant and a green gate).
    let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1).map_err(|error| {
        violations.push(format!(
            "the generic sample {}/{} is not a rational: {error}",
            GENERIC_SAMPLE.0, GENERIC_SAMPLE.1
        ));
        error
    });
    if let Ok(generic) = generic {
        for probe in probes.iter().filter(|probe| probe.parameter == generic) {
            if probe.parameter_kind == Some(ParameterKind::Formal) {
                violations.push(format!(
                    "ordinal {} {}: the generic sample is labelled Formal",
                    probe.ordinal, probe.label
                ));
            }
            if probe.class == TargetClass::Unsupported {
                violations.push(format!(
                    "ordinal {} {}: the generic sample is unsupported",
                    probe.ordinal, probe.label
                ));
            }
        }
    }
    check_blocks(probes, evidence, violations);
}

/// R6.7 card 4: every per-block statistic, asserted on the **stored** data.
///
/// The checks here read the [`BlockStat`]s the probe carries, not the
/// `result.blocks()` slice they were read from, so a block list that is dropped,
/// reordered, mislabelled or given another block's data between `probe_record`
/// and this point is visible.  The two point-level entry points are called a
/// second time on each block's own representative point: that catches a wrong
/// point, a wrong child group or a stale memo entry (the plumbing), but *not* an
/// error inside `point_little_co_group_order` /
/// `point_cocycle_is_a_coboundary` itself -- those are pinned by the
/// `line_domain` unit tests, and no check here can see into them.
fn check_blocks(probes: &[Probe], evidence: &CensusEvidence, violations: &mut Vec<String>) {
    // A **fresh** memo for the second reading: the recorded value came from the
    // run's shared memo, so reusing it here would answer with the same entry
    // whatever point is passed in, and the comparison could not fail on a wrong
    // point.
    let recompute = PointClassCache::default();
    let mut block_total = 0usize;
    let mut answered_probes = 0usize;
    for probe in probes {
        let Some(dimensions) = &probe.dimensions else {
            if !probe.blocks.is_empty() {
                violations.push(format!(
                    "ordinal {} {} t={}: the engine did not answer but the probe carries {} block \
                     statistic(s)",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    probe.blocks.len()
                ));
            }
            continue;
        };
        answered_probes += 1;
        if dimensions.reported_blocks != probe.blocks.len() {
            violations.push(format!(
                "ordinal {} {} t={}: the engine reports {} block(s) but {} statistic(s) were \
                 collected",
                probe.ordinal,
                probe.label,
                probe.parameter,
                dimensions.reported_blocks,
                probe.blocks.len()
            ));
        }
        if dimensions.covered != dimensions.parent {
            violations.push(format!(
                "ordinal {} {} t={}: the blocks cover {} of the parent dimension {}",
                probe.ordinal, probe.label, probe.parameter, dimensions.covered, dimensions.parent
            ));
        }
        if dimensions.arm_total == 0 {
            violations.push(format!(
                "ordinal {} {} t={}: the parent dimension {} gives no arm count",
                probe.ordinal, probe.label, probe.parameter, dimensions.parent
            ));
        }
        if let Some(arms) = probe.partition_arms
            && arms != dimensions.arm_total
        {
            violations.push(format!(
                "ordinal {} {} t={}: the full-star partition lists {arms} arm(s) but the engine's \
                 dimension identity gives {}",
                probe.ordinal, probe.label, probe.parameter, dimensions.arm_total
            ));
        }
        let mut dimension_sum = 0u32;
        let mut arm_seen = vec![0usize; dimensions.arm_total];
        let mut reference_blocks = 0usize;
        for (index, block) in probe.blocks.iter().enumerate() {
            if block.index != index {
                violations.push(format!(
                    "ordinal {} {} t={}: block statistic {index} carries block index {}",
                    probe.ordinal, probe.label, probe.parameter, block.index
                ));
            }
            if block.ordinal != probe.ordinal
                || block.parent_sg != probe.parent_sg
                || block.child_sg != probe.child_sg
                || block.label != probe.label
                || block.parameter != probe.parameter
            {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} is bound to another key \
                     ({}, {}, #{}, {}, t={})",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.ordinal,
                    block.parent_sg,
                    block.child_sg,
                    block.label,
                    block.parameter
                ));
            }
            dimension_sum += block.block_dimension;
            // `star_size` / `arm_count` / `arm_indices` are recomputed from the
            // block's points; the engine's own accessors must agree.
            if block.star_size != block.declared_star_size {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} has {} point(s) but reports star size {}",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.star_size,
                    block.declared_star_size
                ));
            }
            if block.arm_count != block.declared_arm_count {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} sums {} arm(s) over its points but reports \
                     {}",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.arm_count,
                    block.declared_arm_count
                ));
            }
            if block.arm_indices.len() != block.arm_count {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} lists {} distinct arm(s) for {} arm(s)",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.arm_indices.len(),
                    block.arm_count
                ));
            }
            if !block.arm_indices.windows(2).all(|pair| pair[0] < pair[1]) {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} arm list {:?} is not strictly ascending",
                    probe.ordinal, probe.label, probe.parameter, block.arm_indices
                ));
            }
            if block.point_checks != block.star_size {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} compared {} of its {} point(s)",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.point_checks,
                    block.star_size
                ));
            }
            if block.point_disagreements != 0 {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} has {} conjugate-point disagreement(s): its \
                     reported co-group and class describe one point of the star only",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.point_disagreements
                ));
            }
            if block.terms.len() != block.target_count {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} lists {} term(s) for {} target(s)",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.terms.len(),
                    block.target_count
                ));
            }
            // The classifier recomputed from the block's own target counts.  A
            // single small mutation -- classifying every block by the whole
            // probe, the withdrawn convention -- makes this disagree on every
            // probe that mixes sources, and the pinned witness below covers the
            // case where it matters most.
            let recomputed = classify_counts(block.stored_targets, block.target_count);
            if recomputed != block.source {
                violations.push(format!(
                    "ordinal {} {} t={}: block {index} is recorded as {} but its own \
                     {}/{} stored target(s) classify it as {}",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.source.label(),
                    block.stored_targets,
                    block.target_count,
                    recomputed.label()
                ));
            }
            // The engine's own reading of the same block (card-4 audit P1/P2).
            // Without this, a swap of `source`, `stored_targets`, `target_count`
            // and `terms` between two blocks of one probe passed the whole gate
            // for every ordinal except the pinned witnesses, and a
            // `multiplicity + 1` mutation left no trace at all.
            match probe.engine_blocks.get(index) {
                None => violations.push(format!(
                    "ordinal {} {} t={}: block {index} has no engine reading to bind to",
                    probe.ordinal, probe.label, probe.parameter
                )),
                Some(engine) => {
                    let recorded_little: u32 = block
                        .terms
                        .iter()
                        .map(|(dimension, multiplicity)| {
                            u32::from(*dimension).saturating_mul(*multiplicity)
                        })
                        .sum();
                    if block.star_size != engine.star_size
                        || block.arm_count != engine.arm_count
                        || block.arm_indices != engine.arm_indices
                        || block.block_dimension != engine.block_dimension
                        || block.source != engine.source
                        || block.stored_targets != engine.stored_targets
                        || block.target_count != engine.target_count
                        || block.terms != engine.terms
                        || block.representative_point != engine.q
                        || block.carries_reference != engine.carries_reference
                        || block.carries_gamma != engine.carries_gamma
                    {
                        violations.push(format!(
                            "ordinal {} {} t={}: block {index} records star/arms/dimension \
                             {}/{}/{} and source {}/{} of {} with terms {:?}, but the engine \
                             block is {}/{}/{} and {}/{} of {} with terms {:?}",
                            probe.ordinal,
                            probe.label,
                            probe.parameter,
                            block.star_size,
                            block.arm_count,
                            block.block_dimension,
                            block.source.label(),
                            block.stored_targets,
                            block.target_count,
                            block.terms,
                            engine.star_size,
                            engine.arm_count,
                            engine.block_dimension,
                            engine.source.label(),
                            engine.stored_targets,
                            engine.target_count,
                            engine.terms,
                        ));
                        violations.push(format!(
                            "ordinal {} {} t={}: block {index} differs from the engine in \
                             carries_reference {}/{} and carries_gamma {}/{} and representative \
                             point {:?}/{:?}",
                            probe.ordinal,
                            probe.label,
                            probe.parameter,
                            block.carries_reference,
                            engine.carries_reference,
                            block.carries_gamma,
                            engine.carries_gamma,
                            block.representative_point,
                            engine.q,
                        ));
                    }
                    // The engine's own little-group bookkeeping is a second route
                    // to the recorded per-target terms (`Σ dim × mult`), so a
                    // consistent `multiplicity + 1` in both readings still shows
                    // up here (verification review F2).
                    if recorded_little != engine.little_dimension {
                        violations.push(format!(
                            "ordinal {} {} t={}: block {index} records terms {:?} summing to \
                             {recorded_little}, but the engine's little-group dimension is {}",
                            probe.ordinal,
                            probe.label,
                            probe.parameter,
                            block.terms,
                            engine.little_dimension
                        ));
                    }
                    // A representative outside the star it indexes would make the
                    // point-class re-read below meaningless (verification review
                    // F4: moving `representative_point` to the last point of the
                    // block was invisible).
                    if block.star_size == 0 || engine.representative >= block.star_size {
                        violations.push(format!(
                            "ordinal {} {} t={}: block {index} has star size {} but the engine's \
                             representative is {}",
                            probe.ordinal, probe.label, probe.parameter, block.star_size,
                            engine.representative
                        ));
                    }
                }
            }
            match recompute.classify(probe.child_sg, &block.representative_point) {
                Ok(class)
                    if class.little_co_group == block.little_co_group
                        && class.cocycle_trivial == block.cocycle_trivial => {}
                Ok(class) => violations.push(format!(
                    "ordinal {} {} t={}: block {index} records little co-group {} / cocycle \
                     trivial={} but its own point has {} / {}",
                    probe.ordinal,
                    probe.label,
                    probe.parameter,
                    block.little_co_group,
                    block.cocycle_trivial,
                    class.little_co_group,
                    class.cocycle_trivial
                )),
                Err(error) => violations.push(format!(
                    "ordinal {} {} t={}: block {index} point-class recomputation failed: {error}",
                    probe.ordinal, probe.label, probe.parameter
                )),
            }
            for arm in &block.arm_indices {
                match arm_seen.get_mut(*arm) {
                    Some(slot) => *slot += 1,
                    None => violations.push(format!(
                        "ordinal {} {} t={}: block {index} carries arm {arm} outside the \
                         {}-arm list",
                        probe.ordinal, probe.label, probe.parameter, dimensions.arm_total
                    )),
                }
            }
            if block.carries_reference {
                reference_blocks += 1;
            }
        }
        if dimension_sum != dimensions.covered {
            violations.push(format!(
                "ordinal {} {} t={}: the block dimensions sum to {dimension_sum} but the engine \
                 covers {}",
                probe.ordinal, probe.label, probe.parameter, dimensions.covered
            ));
        }
        for (arm, count) in arm_seen.iter().enumerate() {
            if *count != 1 {
                violations.push(format!(
                    "ordinal {} {} t={}: arm {arm} appears in {count} block(s), not exactly one",
                    probe.ordinal, probe.label, probe.parameter
                ));
            }
        }
        // "At most one block carries the reference point, and it must exist":
        // the reference arm belongs to exactly one child star, so a count of two
        // would mean two blocks share a point and a count of zero that the
        // engine's folding lost the reference arm.
        if reference_blocks != 1 {
            violations.push(format!(
                "ordinal {} {} t={}: {reference_blocks} block(s) carry the reference folded point, \
                 expected exactly one",
                probe.ordinal, probe.label, probe.parameter
            ));
        }
        block_total += probe.blocks.len();
    }
    if answered_probes == 0 {
        violations.push("no probe was answered, so no block was classified".to_string());
    }
    if evidence.block_rows != block_total {
        violations.push(format!(
            "the --output-blocks emission loop produced {} row(s) for {block_total} collected \
             block(s)",
            evidence.block_rows
        ));
    }
    if let Some(rows) = evidence.block_file_rows
        && rows != block_total
    {
        violations.push(format!(
            "the --output-blocks file holds {rows} data row(s) for {block_total} collected block(s)"
        ));
    }
    // 10. The pinned card-4 witness (ordinal 13688, SG 225 `SM1` -> child #134,
    //     `t = 1/8`).  The legacy rule labels the whole probe by whether *any*
    //     block carries a constructed target; here the order-16 non-trivial block
    //     is entirely stored while a different block of the same probe is not.
    //     Without this the block-level correction would be untested on the case
    //     that motivated it.
    let eighth = Rat::new(1, 8).expect("1/8");
    let witness: Vec<&Probe> = probes
        .iter()
        .filter(|probe| {
            probe.ordinal == 13_688 && probe.label == "SM1" && probe.parameter == eighth
        })
        .collect();
    if witness.is_empty() {
        violations.push(
            "ordinal 13688 SM1 t=1/8 is gone: the pinned card-4 witness cannot be checked"
                .to_string(),
        );
    }
    for probe in witness {
        if probe.child_sg != 134 {
            violations.push(format!(
                "ordinal 13688 SM1: the child is #{} but the pinned witness is #134",
                probe.child_sg
            ));
        }
        if probe.class != TargetClass::Constructed {
            violations.push(format!(
                "ordinal 13688 SM1 t=1/8: the legacy probe-level rule says {:?}, so the witness no \
                 longer shows the two conventions disagreeing",
                probe.class
            ));
        }
        let clean = probe
            .blocks
            .iter()
            .find(|block| !block.cocycle_trivial && block.little_co_group == 16);
        match clean {
            Some(block) if block.source == BlockSource::Stored => {}
            Some(block) => violations.push(format!(
                "ordinal 13688 SM1 t=1/8: the order-16 non-trivial block is {}",
                block.source.label()
            )),
            None => violations.push(
                "ordinal 13688 SM1 t=1/8: no block has a non-trivial order-16 cocycle".to_string(),
            ),
        }
        if !probe
            .blocks
            .iter()
            .any(|block| block.source != BlockSource::Stored)
        {
            violations.push(
                "ordinal 13688 SM1 t=1/8: every block is stored, so the probe-level rule cannot \
                 have been mislabelling one of them"
                    .to_string(),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cryspglib::irrep::subduction::SubductionComponent;
    use cryspglib::irrep::subduction::star::decompose::line_star_geometry;
    use cryspglib::irrep::subduction::star::line_domain::StarBoundary;

    fn rational(numerator: i128, denominator: i128) -> Rat {
        Rat::new(numerator, denominator).expect("rational")
    }

    /// The census record of one ordinal, through the same walk `records` uses.
    ///
    /// The walk's counted failures are asserted empty here: a test that silently
    /// built its record from a walk that could not read half the corpus would be
    /// testing a different census (card 6).
    fn record_of(parent_sg: u8, ordinal: usize) -> Record {
        let mut failures = Vec::new();
        let records = records_of(parent_sg, &mut failures).expect("the parent's records");
        assert!(failures.is_empty(), "the witness walk must be clean: {failures:?}");
        records
            .into_iter()
            .find(|record| record.ordinal == ordinal)
            .unwrap_or_else(|| panic!("SG {parent_sg} has no parametric-k record {ordinal}"))
    }

    /// The census report of one record, exactly as `run` builds it.
    fn census_report(
        parent_sg: u8,
        ordinal: usize,
        recount: bool,
        domain_sweep: bool,
    ) -> RecordReport {
        let domains = source_domains().expect("the frozen domains");
        let official = official_line_parameter().expect("the official parameter");
        let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1).expect("the generic sample");
        let cache = PointClassCache::default();
        let record = record_of(parent_sg, ordinal);
        probe_record(
            &record,
            &domains,
            official,
            generic,
            recount,
            domain_sweep,
            &cache,
        )
    }

    /// One probe of a report, by source label and parameter.
    fn probe_of<'a>(report: &'a RecordReport, label: &str, parameter: Rat) -> &'a Probe {
        report
            .probes
            .iter()
            .find(|probe| probe.label == label && probe.parameter == parameter)
            .unwrap_or_else(|| panic!("no probe {label} t={parameter}"))
    }

    /// A synthetic target: `irnumber` is the whole stored/constructed
    /// distinction [`classify_block`] reads.
    fn target(irnumber: Option<u32>) -> FullStarTarget {
        FullStarTarget {
            sg: 1,
            ml: None,
            bc: None,
            row_ml: None,
            component: if irnumber.is_some() {
                SubductionComponent::Ordinary
            } else {
                SubductionComponent::Constructed { q: [Rat::ZERO; 3], index: 0 }
            },
            dimension: 1,
            multiplicity: 1,
            irnumber,
        }
    }

    /// **R6.7 card 4, regression (b).**  The classifier is a function of the
    /// slice it is handed and of nothing else, which is the API-level control
    /// for "adding or removing an unrelated block cannot change this block's
    /// label".  The withdrawal the card turns on is shown at the same time: the
    /// legacy probe-level rule applied to the **union** of the two blocks says
    /// `constructed` for a block that is entirely stored.
    #[test]
    fn the_block_classifier_takes_only_its_own_slice() {
        let stored = [target(Some(1)), target(Some(2))];
        let constructed = [target(None)];
        let mixed = [target(Some(1)), target(None)];
        assert_eq!(classify_block(&stored), BlockSource::Stored);
        assert_eq!(classify_block(&constructed), BlockSource::Constructed);
        assert_eq!(classify_block(&mixed), BlockSource::Mixed);
        // The empty slice is `Stored` ("every target is stored" holds
        // vacuously); the report counts empty blocks so the reading is visible.
        assert_eq!(classify_block(&[]), BlockSource::Stored);

        let before = classify_block(&stored);
        // An unrelated block appears, is classified, and disappears again.
        let unrelated = [target(None), target(None)];
        assert_eq!(classify_block(&unrelated), BlockSource::Constructed);
        assert_eq!(
            classify_block(&stored),
            before,
            "classifying another slice must not move this block's label"
        );
        assert_eq!(classify_block(&stored), BlockSource::Stored);

        // The legacy rule on the whole probe: one constructed target anywhere
        // makes the whole probe `constructed`, and with it every block of it.
        let union: Vec<&FullStarTarget> = stored.iter().chain(constructed.iter()).collect();
        assert_eq!(classify_probe(&union), BlockSource::Constructed);
        assert_ne!(
            classify_probe(&union),
            classify_block(&stored),
            "the withdrawn probe-level rule must disagree with the block's own slice"
        );
    }

    /// **R6.7 card 4, regression (a): the pinned witness.**  Ordinal 13688
    /// (SG 225 `SM1` -> child #134) at `t = 1/8` carries an order-16 block whose
    /// own cocycle is non-trivial and whose targets are **all stored**, while a
    /// different block of the same probe carries constructed targets.  The
    /// legacy probe-level rule therefore labels this probe `constructed` and
    /// would have said the same about the non-trivial block -- the accepted
    /// finding this card exists to remove.
    #[test]
    fn the_13688_witness_is_a_block_level_correction() {
        let report = census_report(225, 13_688, false, false);
        assert!(
            report.failures.is_empty(),
            "the record must probe cleanly: {:?}",
            report.failures
        );
        let probe = probe_of(&report, "SM1", rational(1, 8));
        assert_eq!(probe.child_sg, 134, "the pinned witness child");
        assert_eq!(
            probe.class,
            TargetClass::Constructed,
            "the legacy probe-level rule is what the finding withdraws"
        );
        let own: Vec<BlockSource> = probe.blocks.iter().map(|block| block.source).collect();
        assert!(
            own.contains(&BlockSource::Stored) && own.iter().any(|s| *s != BlockSource::Stored),
            "the probe must mix block sources, got {own:?}"
        );
        let clean = probe
            .blocks
            .iter()
            .find(|block| !block.cocycle_trivial && block.little_co_group == 16)
            .expect("the order-16 non-trivial block");
        assert_eq!(
            clean.source,
            BlockSource::Stored,
            "the non-trivial block is entirely stored and must not inherit the probe's label"
        );
        assert!(clean.carries_reference, "the reference block is the order-16 one");
        assert_eq!(
            clean.little_co_group, probe.child_order,
            "the reference-carrying block's own co-group is the census's order at the reference point"
        );
        // The block's own class is the class of the point it carries, read
        // again from that point: the plane the finding is about is exactly this
        // one -- the reference point is non-trivial, and the other blocks are not
        // (one of them is constructed), so a single probe-level label cannot
        // describe both.
        assert!(
            !point_cocycle_is_a_coboundary(134, &clean.representative_point)
                .expect("the block's own class"),
            "the order-16 block's own point is non-trivial"
        );
        assert!(
            probe.blocks.len() >= 2,
            "the witness needs at least two blocks, got {}",
            probe.blocks.len()
        );
        assert!(
            probe
                .blocks
                .iter()
                .any(|block| block.index != clean.index && block.source != clean.source),
            "a different block of the same probe carries the other source"
        );
    }

    /// **Card 4, the Gamma exception branch, exercised synthetically.**  Every
    /// Gamma-reaching parameter of the corpus is a full-star boundary (measured:
    /// 23,024/23,024 contained, 0 exceptional), so the documented exception --
    /// an arm fixed *exactly* by the whole child point group, whose Gamma
    /// parameters can then be interior to an interval -- has no corpus witness,
    /// and neither has the failure next to it.  A synthetic partition (the
    /// fields are public) exercises both, so the branch is not dead code that
    /// only looks like a check.
    #[test]
    fn the_gamma_exception_requires_an_exactly_fixed_arm() {
        let child_reciprocal = reciprocal_lattice(1).expect("P1's reciprocal lattice");
        let quarter_turn: Mat3I = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let partition = |rotations: Vec<Mat3I>| FullStarPartition {
            source_sg: 196,
            child_sg: 1,
            label: "DT1",
            // `(2, 0, 0)` reaches Gamma exactly at `t = 0` and `t = 1/2` in Z^3.
            arms: vec![FoldedArm {
                direction: Vec3R::from_ints([2, 0, 0]),
                parent_rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            }],
            reference_arm: 0,
            child_reciprocal,
            child_rotations: rotations,
            // `t = 0` is cut, `t = 1/2` is deliberately interior: the child is
            // P1, so no rotation ever enters the little co-group and the
            // partition has nothing to cut at `1/2`.
            boundaries: vec![StarBoundary {
                parameter: Rat::ZERO,
                counts: [0; 3],
                witnesses: Vec::new(),
            }],
            permanent_counts: [0; 3],
            permanent_witnesses: Vec::new(),
        };

        // P1's point group is trivial, so every arm is fixed by all of it: the
        // interior parameter is explained and counted as the exception.
        let mut report = GammaReport::default();
        let mut failures = Vec::new();
        let (entries, zero_arms) = report_gamma(
            &partition(vec![[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]),
            196,
            &mut report,
            &mut failures,
        );
        assert!(failures.is_empty(), "the explained exception is not a failure: {failures:?}");
        assert!(zero_arms.is_empty(), "the arm does not fold to zero");
        assert_eq!(entries.len(), 2, "t = 0 and t = 1/2 both reach Gamma");
        assert_eq!(report.entries, 2);
        assert_eq!(report.contained, 1, "t = 0 is a boundary");
        assert_eq!(report.exceptional, 1, "t = 1/2 is interior and explained");
        assert_eq!(report.witnesses.len(), 1, "the exception keeps its witness");

        // The same interior parameter without exact fixity: a rotation of the
        // child moves the arm, so the exception does not apply and the report
        // must fail instead of quietly counting it.
        let mut report = GammaReport::default();
        let mut failures = Vec::new();
        let (entries, _) = report_gamma(&partition(vec![quarter_turn]), 196, &mut report, &mut failures);
        assert_eq!(entries.len(), 2);
        assert_eq!(report.exceptional, 0, "a moved arm is not an exception");
        assert_eq!(failures.len(), 1, "the unexplained interior parameter is a failure");
        assert!(
            failures[0].contains("no reaching arm is fixed exactly"),
            "the failure names the missing premise: {failures:?}"
        );
    }

    /// **R6.7 card 4, regression (c): the arm merge of ordinal 10038.**
    /// SG 196 `DT1` -> P1 has a trivial child little co-group everywhere, so no
    /// cocycle changes; the geometry still does: six one-armed child stars at
    /// the generic `t = 1/9` and four stars with arm counts `2, 1, 1, 2` at
    /// `t = 1/8`.  The production decomposition is checked on the same two
    /// parameters, so the geometry cannot drift away from the reported blocks.
    #[test]
    fn the_10038_witness_merges_its_arms_at_one_eighth() {
        let record = record_of(196, 10_038);
        let embedding =
            SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup).expect("embedding");
        assert_eq!(embedding.subgroup_sg(), 1, "the P1 child of the witness");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen DT1 table");
        for (parameter, expected) in [
            (rational(1, 9), vec![1usize, 1, 1, 1, 1, 1]),
            (rational(1, 7), vec![1, 1, 1, 1, 1, 1]),
            (rational(1, 8), vec![1, 1, 2, 2]),
        ] {
            let geometry = line_star_geometry(&record.subgroup, &embedding, table, parameter)
                .expect("the production geometry");
            assert_eq!(geometry.len(), expected.len(), "child stars at t={parameter}");
            let mut arms: Vec<usize> = geometry.iter().map(|star| star.arm_count()).collect();
            arms.sort_unstable();
            assert_eq!(arms, expected, "arm counts at t={parameter}");
            assert!(
                geometry.iter().all(|star| star.star_size() == 1),
                "each child star of this line is a single folded point"
            );
        }
        // The reported decomposition of the same two parameters.  `1/8` is not a
        // parameter the reference-partition census probes for this record -- that
        // is the card-3 finding -- so it is asked for through the production
        // entry point, and the census side is checked through the recount report
        // just below.
        for (parameter, expected) in [
            (rational(1, 9), vec![1usize, 1, 1, 1, 1, 1]),
            (rational(1, 7), vec![1, 1, 1, 1, 1, 1]),
            (rational(1, 8), vec![1, 1, 2, 2]),
        ] {
            let result =
                subduce_line_at_parameter(&record.subgroup, &embedding, table, parameter)
                    .expect("the production decomposition");
            let mut counts: Vec<usize> =
                result.blocks().iter().map(|block| block.arm_count()).collect();
            counts.sort_unstable();
            assert_eq!(counts, expected, "the engine's blocks at t={parameter}");
            // The reported terms are pinned too, so the census test does not defer
            // the whole structure to the library test (verification review F10).
            let expected_terms: Vec<(u8, u32)> = if expected.len() == 6 {
                vec![(1, 1); 6]
            } else {
                vec![(1, 1), (1, 1), (1, 2), (1, 2)]
            };
            let mut terms: Vec<(u8, u32)> = result
                .blocks()
                .iter()
                .flat_map(|block| {
                    block
                        .targets()
                        .iter()
                        .map(|target| (target.dimension, target.multiplicity))
                })
                .collect();
            terms.sort_unstable();
            assert_eq!(terms, expected_terms, "the engine's terms at t={parameter}");
            // No per-block class or order assertion here: the child is P1, so
            // `point_cocycle_is_a_coboundary(1, _)` and
            // `point_little_co_group_order(1, _)` are constant at *every* point by
            // construction and would assert nothing (card-4 audit F8, the same
            // tautology the card-3 math review removed from the library test).
        }
        let report = census_report(196, 10_038, true, false);
        assert!(
            report.failures.is_empty() && report.recount.failures.is_empty(),
            "the recount of the witness must be clean: {:?} {:?}",
            report.failures,
            report.recount.failures
        );
        assert!(
            report.recount.changed_pairs.contains(&(10_038, "DT1")),
            "the recount must see the reference partition's blind spot on this pair"
        );
        assert!(
            report.recount.probes >= 2,
            "the recount must add the merge parameter and its neighbours, got {}",
            report.recount.probes
        );
    }
    // ── R6.7 card 5: the interval sweep ──────────────────────────────────────

    /// A synthetic target for the matcher: only the fields the comparison reads.
    fn sweep_target(
        dimension: u8,
        multiplicity: u32,
        stored: bool,
        vector: [f64; 2],
    ) -> SweepTarget {
        SweepTarget {
            dimension,
            multiplicity,
            stored,
            vector: vector
                .iter()
                .map(|value| Complex64::new(*value, 0.0))
                .collect(),
        }
    }

    /// **Card 5, regression (c): two distinct one-dimensional targets with
    /// swapped multiplicities must fail.**  The two characters (`[1, 1]` and
    /// `[1, -1]`) are orthogonal, so the match is by character and one-to-one;
    /// the swap leaves the multiset `{1x2, 1x1}` unchanged, so a comparison that
    /// sorted the `(dimension, multiplicity)` table would see nothing.  The
    /// positive control reverses the order **and** flips the provenance: the
    /// same targets must then match, which is the "a stored and a constructed
    /// representation of the same target must match" regression.
    #[test]
    fn the_matcher_uses_characters_and_not_a_sorted_term_table() {
        let left = [
            sweep_target(1, 2, true, [1.0, 1.0]),
            sweep_target(1, 1, false, [1.0, -1.0]),
        ];
        let swapped = [
            sweep_target(1, 1, false, [1.0, 1.0]),
            sweep_target(1, 2, true, [1.0, -1.0]),
        ];
        match match_targets(&left, &swapped, "the swap witness") {
            Err(TargetMismatch::Multiplicity(message)) => {
                assert!(
                    message.contains("multiplicity 2 at the first parameter and 1"),
                    "the failure must name the swapped multiplicity: {message}"
                );
            }
            other => panic!("the swapped multiplicities must fail, got {other:?}"),
        }

        // The same two targets in the other order, with the provenance moved
        // with them: matched by character, not by position or by label.
        let reordered = [
            sweep_target(1, 1, true, [1.0, -1.0]),
            sweep_target(1, 2, false, [1.0, 1.0]),
        ];
        let (pairs, worst) =
            match_targets(&left, &reordered, "the reordered control").expect("the control matches");
        assert_eq!(pairs, vec![(0, 1), (1, 0)], "the pairing is by character");
        assert!(worst >= 1.0 - SWEEP_TOLERANCE, "the matched score is 1, got {worst}");

        // Two targets with the same character are not distinguishable: an
        // explicit error, never a silent pick.
        let duplicated = [
            sweep_target(1, 1, false, [1.0, 1.0]),
            sweep_target(1, 1, true, [1.0, 1.0]),
        ];
        match match_targets(&left, &duplicated, "the ambiguity witness") {
            Err(TargetMismatch::Ambiguous(message)) => {
                assert!(message.contains("partners"), "the tie must be named: {message}");
            }
            other => panic!("a character tie must be an explicit error, got {other:?}"),
        }

        // External review of the card-5 line, P2: a global `-1` is **not** the
        // same target.  The matcher used to take the absolute value of the inner
        // product, which accepted `[1, 1]` against `[-1, -1]` with score ~1; the
        // signed score is `-1` there, so there is no partner at all.
        let flipped = [
            sweep_target(1, 2, true, [-1.0, -1.0]),
            sweep_target(1, 1, false, [-1.0, 1.0]),
        ];
        match match_targets(&left, &flipped, "the sign-flip witness") {
            Err(TargetMismatch::Ambiguous(message)) => {
                assert!(message.contains("no character partner"), "{message}");
            }
            other => panic!("a global sign flip must not match, got {other:?}"),
        }
        // The reviewer's exact corpus mutation is the same thing at vector level:
        // one honest vector against its own negation scores `-1`, not `1`.
        let honest = sweep_target(1, 1, false, [1.0, 1.0]);
        let negated = sweep_target(1, 1, false, [-1.0, -1.0]);
        assert!(
            character_score(&honest.vector, &negated.vector) < -0.9,
            "the signed score must separate a character from its negative"
        );
        assert!(
            character_score(&honest.vector, &honest.vector) > 0.9,
            "and it must still accept the character itself"
        );

        // A target with no character partner is an error too, not a dropped row.
        let missing = [
            sweep_target(1, 2, true, [1.0, 1.0]),
            sweep_target(1, 1, false, [0.0, 1.0]),
        ];
        match match_targets(&left, &missing, "the missing witness") {
            Err(TargetMismatch::Ambiguous(message)) => {
                assert!(message.contains("no character partner"), "{message}");
            }
            other => panic!("an unmatched target must be an explicit error, got {other:?}"),
        }

        // A dimension disagreement between character-equal targets is caught.
        let wrong_dimension = [
            sweep_target(2, 2, true, [1.0, 1.0]),
            sweep_target(1, 1, false, [1.0, -1.0]),
        ];
        match match_targets(&left, &wrong_dimension, "the dimension witness") {
            Err(TargetMismatch::Ambiguous(message)) => {
                assert!(message.contains("dimension 1 at the first parameter and 2"), "{message}");
            }
            other => panic!("a dimension change must be an explicit error, got {other:?}"),
        }
    }

    /// One synthetic block geometry for the layer-1 comparison.
    fn sweep_block(
        arm_indices: Vec<usize>,
        points: Vec<Vec<usize>>,
        rotations: Vec<Mat3I>,
    ) -> SweepBlock {
        let star_size = points.len();
        let arm_count = points.iter().map(Vec::len).sum();
        SweepBlock {
            arm_indices,
            points,
            star_size,
            arm_count,
            block_dimension: u32::try_from(arm_count).expect("small"),
            rotations,
        }
    }

    /// **Card 5, layer 1: blocks are matched by parent arm set, never by block
    /// index.**  The two parameters' block orders really do move (the folded
    /// coordinates, and with them the canonical order, are functions of `t`), so
    /// an index-based comparison would pair `[1, 2]` with `[3]` and report a
    /// mismatch on a pair that is actually the same geometry -- or, worse, pair
    /// two different blocks and call it a match.  The pairing itself is asserted,
    /// which is what an index-based matcher cannot satisfy.
    #[test]
    fn the_geometry_layer_matches_blocks_by_arm_set_and_not_by_index() {
        let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let left = [
            sweep_block(vec![1, 2], vec![vec![1], vec![2]], vec![identity]),
            sweep_block(vec![3], vec![vec![3]], vec![identity]),
        ];
        let reordered = [
            sweep_block(vec![3], vec![vec![3]], vec![identity]),
            sweep_block(vec![1, 2], vec![vec![1], vec![2]], vec![identity]),
        ];
        let matching = match_block_geometry(&left, &reordered, "the reordered control")
            .expect("the reordered blocks match by arm set");
        assert_eq!(
            matching,
            vec![(0, 1), (1, 0)],
            "the identity has to be the arm set, not the position"
        );

        // The star's internal grouping moves too: the same arms and the same
        // star size, folded onto different points.
        let wide_left = [
            sweep_block(vec![1, 2, 3, 4], vec![vec![1, 2], vec![3, 4]], vec![identity]),
            sweep_block(vec![5], vec![vec![5]], vec![identity]),
        ];
        let regrouped = [
            sweep_block(vec![1, 2, 3, 4], vec![vec![1, 3], vec![2, 4]], vec![identity]),
            sweep_block(vec![5], vec![vec![5]], vec![identity]),
        ];
        let error = match_block_geometry(&wide_left, &regrouped, "the regrouped witness")
            .expect_err("a regrouped star is not the same block");
        assert!(
            error.contains("folds its parent arms onto different points"),
            "the failure must name the offending block: {error}"
        );

        // A little co-group that grew between the two parameters.
        let grown = [
            sweep_block(
                vec![1, 2],
                vec![vec![1], vec![2]],
                vec![[[0, -1, 0], [1, 0, 0], [0, 0, 1]]],
            ),
            sweep_block(vec![3], vec![vec![3]], vec![identity]),
        ];
        let error = match_block_geometry(&left, &grown, "the co-group witness")
            .expect_err("a grown little co-group is not the same block");
        assert!(error.contains("little co-group order"), "{error}");

        // A different block count is named before any pairing.
        let short = [sweep_block(vec![1, 2], vec![vec![1], vec![2]], vec![identity])];
        let error = match_block_geometry(&left, &short, "the count witness")
            .expect_err("a different block count is a mismatch");
        assert!(error.contains("2 block(s)"), "{error}");

        // Card-5 audit F4: the volume comparisons need a falsifier of their own
        // (deleting them used to leave every test and the whole corpus green).
        // Each is flipped in isolation, with the same arms, points and
        // co-group, so only the named quantity differs.
        let base: SweepBlock = sweep_block(vec![1, 2], vec![vec![1], vec![2]], vec![identity]);
        let mut grown_star = base.clone();
        grown_star.star_size += 1;
        let error = match_block_geometry(
            std::slice::from_ref(&base),
            std::slice::from_ref(&grown_star),
            "the star-size witness",
        )
        .expect_err("a different star size is a mismatch");
        assert!(error.contains("star size 2 at the first parameter and 3"), "{error}");
        let mut heavier = base.clone();
        heavier.arm_count += 1;
        let error = match_block_geometry(
            std::slice::from_ref(&base),
            std::slice::from_ref(&heavier),
            "the arm-count witness",
        )
        .expect_err("a different arm count is a mismatch");
        assert!(error.contains("arm count 2 at the first parameter and 3"), "{error}");
        let mut fatter = base.clone();
        fatter.block_dimension += 1;
        let error = match_block_geometry(
            std::slice::from_ref(&base),
            std::slice::from_ref(&fatter),
            "the block-dimension witness",
        )
        .expect_err("a different block dimension is a mismatch");
        assert!(
            error.contains("block dimension 2 at the first parameter and 3"),
            "{error}"
        );
    }

    /// The frozen record's own `(record, label)` pair list, for the sweep tests.
    fn sweep_context_of(
        record: &Record,
        label: &'static str,
    ) -> (SubgroupEmbedding, FullStarPartition) {
        let embedding =
            SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup).expect("embedding");
        let table = line_table(record.parent_sg, label).expect("the frozen table");
        let partition = full_star_partition(&record.subgroup, &embedding, table)
            .expect("the full-star partition");
        (embedding, partition)
    }

    /// **Card 5, regressions (a) and (b): the gauge is load-bearing, and a pair
    /// that differs only by the moving constructed point must not fail.**
    ///
    /// Ordinal 10030 `DT1`'s interval `(1/8, 1/4)` is the measured witness: the
    /// block with arms `[4, 5]` reports a *different* constructed point at
    /// `t = 1/6` (`q = 1/3`) and `t = 5/24` (`q = 1/6`), and its order-two little
    /// co-group moves the induced full-star character from `-1` to `-sqrt(3)`.
    /// The gauge-unified per-arm characters agree exactly, which is the identity
    /// the sweep's layer 2 rests on; without the gauge the same comparison fails
    /// (measured in the card-5 mutation run).
    #[test]
    fn the_gauge_unifies_a_moving_constructed_point() {
        let record = record_of(196, 10_030);
        let (embedding, partition) = sweep_context_of(&record, "DT1");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let t1 = rational(1, 6);
        let t2 = rational(5, 24);
        let first = subduce_line_at_parameter(&record.subgroup, &embedding, table, t1)
            .expect("the decomposition at 1/6");
        let second = subduce_line_at_parameter(&record.subgroup, &embedding, table, t2)
            .expect("the decomposition at 5/24");
        let block = first
            .blocks()
            .iter()
            .find(|block| block.arm_indices() == vec![4, 5])
            .expect("the witness block at 1/6");
        let other = second
            .blocks()
            .iter()
            .find(|block| block.arm_indices() == vec![4, 5])
            .expect("the witness block at 5/24");
        let q1 = seed_point(block, 4, "the witness").expect("the seed point at 1/6");
        let q2 = seed_point(other, 4, "the witness").expect("the seed point at 5/24");
        assert_ne!(q1, q2, "the constructed point must move between the two parameters");
        assert_eq!(q1, Vec3R::new([rational(-2, 3), Rat::ZERO, Rat::ZERO]));
        assert_eq!(q2, Vec3R::new([rational(-5, 6), Rat::ZERO, Rat::ZERO]));

        let operations = child_little_group(&embedding, &partition.child_reciprocal, &q1)
            .expect("the little group at 1/6");
        assert_eq!(operations.len(), 2, "the witness little co-group has order two");
        let other_operations = child_little_group(&embedding, &partition.child_reciprocal, &q2)
            .expect("the little group at 5/24");
        assert_eq!(
            operations
                .iter()
                .map(|operation| rotation_key(&operation.rotation()))
                .collect::<Vec<_>>(),
            other_operations
                .iter()
                .map(|operation| rotation_key(&operation.rotation()))
                .collect::<Vec<_>>(),
            "the same finite group at both parameters, aligned by rotation"
        );

        // The raw presentation moves: the induced character at the non-identity
        // operation is `-1` at `1/6` and `-sqrt(3)` at `5/24`.
        let nontrivial = operations[1];
        let raw1 = block
            .target_character(0, &nontrivial)
            .expect("the induced character at 1/6");
        let raw2 = other
            .target_character(0, &nontrivial)
            .expect("the induced character at 5/24");
        assert!(
            (raw1 - raw2).norm() > 0.1,
            "the raw presentation must move: {raw1} against {raw2}"
        );
        // The little-group reading at the common seed arm is what the sweep
        // compares, and the gauge removes the motion.
        let targets = gauge_unified_targets(block, &q1, &operations).expect("the targets at 1/6");
        let other_targets =
            gauge_unified_targets(other, &q2, &other_operations).expect("the targets at 5/24");
        assert_eq!(targets.len(), 1);
        let score = character_score(&targets[0].vector, &other_targets[0].vector);
        assert!(
            score >= 1.0 - 1e-12,
            "the gauge-unified characters must agree, score {score}"
        );
        let (pairs, worst) = match_targets(&targets, &other_targets, "the witness")
            .expect("the pair must match");
        assert_eq!(pairs, vec![(0, 0)]);
        assert!(worst >= 1.0 - 1e-12);
        // The move is not a no-op: the ungauged per-arm values differ.
        let ungauged1 = block
            .target_little_character(0, &q1, &nontrivial)
            .expect("the little-group character at 1/6");
        let ungauged2 = other
            .target_little_character(0, &q2, &nontrivial)
            .expect("the little-group character at 5/24");
        assert!(
            (ungauged1 - ungauged2).norm() > 1e-9,
            "the gauge has to be doing work: {ungauged1} against {ungauged2}"
        );
    }

    /// **Card 5, regression (d): the interval logic is load-bearing, and the
    /// 10038 structures are pinned through the sweep.**
    ///
    /// Ordinal 10038 `DT1` merges its arms exactly at `t = 1/8`: six one-armed
    /// blocks at `1/9` and `1/7`, four blocks with arm counts `2, 1, 1, 2` at
    /// `1/8`.  The sweep only compares **inside** an interval, so comparing the
    /// two parameters across the boundary must fail -- which is asserted here on
    /// the very same production decompositions, so a sweep that used the wrong
    /// interval set (or none) cannot pass this test.
    #[test]
    fn the_sweep_intervals_are_what_keeps_the_10038_merge_out() {
        let record = record_of(196, 10_038);
        let (_, partition) = sweep_context_of(&record, "DT1");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let domains = source_domains().expect("the frozen domains");
        let domain = domains
            .get(&(record.parent_sg, "DT1"))
            .expect("the DT1 domain");
        let parameters = partition.probe_parameters(domain).expect("the probe set");
        let intervals = probe_intervals(&parameters).expect("the intervals");
        // `Rat` deliberately has no `Ord`; compare by cross-multiplication with
        // the (normalized, positive) denominators.
        let greater = |first: &Rat, second: &Rat| {
            first.numerator() * second.denominator() > second.numerator() * first.denominator()
        };
        let inside = |parameter: &Rat| -> Option<(Rat, Rat)> {
            intervals
                .iter()
                .find(|(left, right)| greater(parameter, left) && greater(right, parameter))
                .copied()
        };
        let ninth = rational(1, 9);
        let seventh = rational(1, 7);
        let eighth = rational(1, 8);
        // Measured: `probe_parameters(parent)` is the eighth grid, so 1/9 and
        // 1/7 are interior to *different* intervals -- `(0, 1/8)` and
        // `(1/8, 1/4)`.  The merge parameter 1/8 is the endpoint between them and
        // is therefore never compared by the sweep.
        assert_eq!(inside(&ninth), Some((Rat::ZERO, eighth)));
        assert_eq!(inside(&seventh), Some((eighth, rational(1, 4))));
        assert_ne!(inside(&ninth), inside(&seventh));
        assert!(
            inside(&eighth).is_none(),
            "1/8 is the merge boundary and must not be interior to any interval"
        );
        assert!(
            intervals
                .iter()
                .any(|(left, right)| *left == eighth || *right == eighth),
            "1/8 must be an endpoint of an interval"
        );

        // Comparisons across the boundary really do have different geometry, so
        // the interval logic is what keeps them out of the sweep.
        let keep = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup).expect("embedding");
        let decompose = |parameter: Rat| {
            subduce_line_at_parameter(&record.subgroup, &keep, table, parameter)
                .expect("the production decomposition")
        };
        let at_boundary = decompose(eighth);
        let at_inside = decompose(ninth);
        let context = SweepContext {
            ordinal: record.ordinal,
            label: "DT1",
            child_sg: record.child_sg,
            child_reciprocal: &partition.child_reciprocal,
            child_rotations: &partition.child_rotations,
            arms: &partition.arms,
        };
        let read = |parameter: Rat, result: &LineSubduction| -> Vec<SweepBlock> {
            result
                .blocks()
                .iter()
                .map(|block| {
                    SweepBlock::read(&context, parameter, block).expect("the block geometry")
                })
                .collect()
        };
        let error = match_block_geometry(
            &read(ninth, &at_inside),
            &read(eighth, &at_boundary),
            "the cross-boundary witness",
        )
        .expect_err("the merge must be visible across the boundary");
        assert!(
            error.contains("6 block(s) at the first parameter against 4"),
            "the merge is a block-count change: {error}"
        );

        // And the sweep's own report for the witness: eight intervals, all of
        // them compared, no failures.
        let report = census_report(196, 10_038, true, true);
        assert!(
            report.failures.is_empty() && report.sweep.failures.is_empty(),
            "the sweep of the witness must be clean: {:?} {:?}",
            report.failures,
            report.sweep.failures
        );
        assert_eq!(report.sweep.pairs, record.labels.len());
        assert_eq!(report.sweep.intervals, 8 * record.labels.len());
        assert_eq!(report.sweep.points, 2 * report.sweep.intervals);
        assert_eq!(report.sweep.comparisons, report.sweep.intervals);
        assert_eq!(report.sweep.failed_intervals, 0);
        assert_eq!(report.sweep.sources[0], 0, "this child has no stored interior target");
    }

    /// **Card 5: every block is bound to the parameter whose decomposition
    /// produced it, and the two interior points are two different parameters.**
    ///
    /// A comparison whose two sides are read from one decomposition satisfies
    /// every other guard of the sweep: measured on the corpus (`M5`, author's
    /// mutation, card-5 report) taking both sides from the first parameter left
    /// the gate at exit 0 with 46,048 "comparisons" and 0 failures.  The binding
    /// that closes it is asserted here at unit level, together with the
    /// degenerate-interval guard that keeps a zero-length interval from comparing
    /// one boundary point with itself.
    #[test]
    fn the_sweep_binds_every_block_to_its_own_parameter() {
        let record = record_of(196, 10_030);
        let (embedding, partition) = sweep_context_of(&record, "DT1");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let context = SweepContext {
            ordinal: record.ordinal,
            label: "DT1",
            child_sg: record.child_sg,
            child_reciprocal: &partition.child_reciprocal,
            child_rotations: &partition.child_rotations,
            arms: &partition.arms,
        };
        let first = subduce_line_at_parameter(&record.subgroup, &embedding, table, rational(1, 6))
            .expect("the decomposition at 1/6");
        let block = first
            .blocks()
            .iter()
            .find(|block| block.arm_indices() == vec![4, 5])
            .expect("the witness block at 1/6");
        SweepBlock::read(&context, rational(1, 6), block)
            .expect("the block is accepted at its own parameter");
        let error = SweepBlock::read(&context, rational(5, 24), block)
            .expect_err("a block of one parameter must not pass as the other's");
        assert!(error.contains("folded direction"), "{error}");

        // Two different interior points, strictly inside the interval.
        let (step, low, high) = interior_points(&Rat::ZERO, &rational(1, 4)).expect("trisection");
        assert_eq!(step, rational(1, 12));
        assert_eq!(low, rational(1, 12));
        assert_eq!(high, rational(1, 6));
        assert_ne!(low, high, "the two interior points must be different parameters");

        // A zero-length or reversed interval is a typed error, never a
        // self-comparison counted as a checked interval.
        for (left, right) in [
            (rational(1, 4), rational(1, 4)),
            (rational(1, 4), rational(1, 8)),
        ] {
            let error = interior_points(&left, &right).expect_err("a non-positive interval");
            assert!(error.contains("not positive"), "{error}");
        }
    }

    /// **Card 5 math review P1-1: the card-1 gauge is derived for constructed
    /// targets only, so a stored target fails closed instead of being compared
    /// against an unproved quantity.**  The corpus never reaches the branch (no
    /// interior stored target exists: `stored=0` of 174,672), so this test is the
    /// guard's only witness.
    #[test]
    fn the_stored_branch_of_the_gauge_fails_closed() {
        let point = Vec3R::new([Rat::ZERO, Rat::ZERO, Rat::ZERO]);
        let operations = [ExactSeitz::identity()];
        let constructed = [target(None)];
        gauge_target_slice(&constructed, &point, &operations, |_, _| {
            Ok(Complex64::new(1.0, 0.0))
        })
        .expect("a constructed target is gaugeable");
        // The guard sits on the path `gauge_unified_targets` uses, and it fires
        // before any character is read: the reader panics if it is reached
        // (audit of `d66db48`, P0-3 -- a deletable `refuse_stored_targets` call
        // left this branch unbound).
        // The stored target comes first so the reader below cannot run for a
        // constructed neighbour before the guard is reached.
        let stored = [target(Some(7)), target(None)];
        let error = gauge_target_slice(&stored, &point, &operations, |_, _| {
            panic!("the stored guard must fail before any character is read")
        })
        .expect_err("a stored target must fail closed");
        assert!(error.contains("stored irrep 7"), "{error}");
        assert!(error.contains("card-1 gauge"), "{error}");
    }

    /// **Card 5 math review P1-3: the control that pins the gauge's exact-fixity
    /// hypothesis is not vacuous.**  `Lattice::preserves` builds the little group
    /// from mod-`L*` fixity; the gauge derivation needs exact fixity.  The same
    /// point separates the two: a quarter turn about `z` fixes `(1/2, 1/2, 0)`
    /// modulo `Z^3` but not exactly, and fixes `(0, 0, 1/3)` both ways.
    #[test]
    fn the_exact_fixity_control_distinguishes_mod_lattice_fixity() {
        let lattice = reciprocal_lattice(1).expect("SG 1 reciprocal lattice");
        let rotation: Mat3I = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let class_point = Vec3R::new([rational(1, 2), rational(1, 2), Rat::ZERO]);
        assert!(
            lattice
                .preserves(rotation, &class_point)
                .expect("the mod-lattice fixity test"),
            "the point is fixed modulo L*, which is why the little group contains the rotation"
        );
        assert!(
            !exactly_fixes(&rotation, &class_point).expect("the exact fixity test"),
            "the rotation must not count as an exact fixity of this point"
        );
        let axis_point = Vec3R::new([Rat::ZERO, Rat::ZERO, rational(1, 3)]);
        assert!(
            lattice
                .preserves(rotation, &axis_point)
                .expect("the mod-lattice fixity test")
        );
        assert!(
            exactly_fixes(&rotation, &axis_point).expect("the exact fixity test"),
            "a point on the rotation axis is fixed exactly"
        );
        let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        assert!(exactly_fixes(&identity, &class_point).expect("the identity is exact"));
    }

    /// **External review of the card-5 line, P2: the gauge-unified vector has to
    /// be a character, not merely a vector of the right length.**  `chi(E)` is the
    /// reported dimension (a positive real) and every entry is finite; a vector
    /// that fails either check is an explicit error instead of something the
    /// normalized inner product could still match against an honest character.
    #[test]
    fn the_gauge_validates_the_character_at_the_identity() {
        let point = Vec3R::new([Rat::ZERO, Rat::ZERO, Rat::ZERO]);
        // The aligned little group of a point always contains the identity.
        let operations = [ExactSeitz::identity()];
        let constructed = [target(None)];
        let honest = gauge_target_slice(&constructed, &point, &operations, |_, _| {
            Ok(Complex64::new(1.0, 0.0))
        })
        .expect("chi(E) = dimension is accepted");
        assert_eq!(honest.len(), 1);
        let error = gauge_target_slice(&constructed, &point, &operations, |_, _| {
            Ok(Complex64::new(-1.0, 0.0))
        })
        .expect_err("chi(E) = -1 is not the dimension of a one-dimensional target");
        assert!(error.contains("chi(E)"), "{error}");
        let error = gauge_target_slice(&constructed, &point, &operations, |_, _| {
            Ok(Complex64::new(f64::NAN, 0.0))
        })
        .expect_err("a non-finite character entry is rejected");
        assert!(error.contains("non-finite"), "{error}");
        // And a little group without a zero-translation identity is an error
        // rather than a silent skip of the normalization check.
        let shifted = [ExactSeitz::new(
            IDENTITY_ROTATION,
            Vec3R::new([Rat::ZERO, Rat::ZERO, Rat::new(1, 3).unwrap()]),
        )];
        let error = gauge_target_slice(&constructed, &point, &shifted, |_, _| {
            Ok(Complex64::new(1.0, 0.0))
        })
        .expect_err("a missing zero-translation identity is an error");
        assert!(error.contains("zero-translation identity"), "{error}");
    }

    /// **Card 4 residual: the `--output-recount` fingerprint round-trips.**  The
    /// gate writes the file and re-reads it; this pins the two halves at the unit
    /// level, including that a corrupted field changes the parsed row (a row
    /// count alone would not notice).
    #[test]
    fn the_recount_fingerprint_rows_round_trip() {
        let report = census_report(196, 10_038, true, false);
        assert_eq!(
            report.recount.rows.len(),
            report.recount.probes,
            "one fingerprint row per recount probe"
        );
        assert!(report.recount.probes > 0, "the witness has recount probes");
        for (index, row) in report.recount.rows.iter().enumerate() {
            let text = row.to_tsv();
            let parsed = RecountRow::parse(&text, index + 1).expect("the row parses back");
            assert_eq!(&parsed, row, "row {index} must round-trip");
        }
        let mut corrupted: Vec<String> = report.recount.rows[0]
            .to_tsv()
            .split('\t')
            .map(String::from)
            .collect();
        corrupted[5] = (report.recount.rows[0].blocks.len() + 1).to_string();
        let error = RecountRow::parse(&corrupted.join("\t"), 1)
            .expect_err("an inconsistent block count must be rejected");
        assert!(error.contains("announced"), "{error}");
    }

    // ── R6.7 card 6: the boundary pass and the arm-geometry audit ────────────

    /// One `(record, label)`'s sweep under an arbitrary partition: the frozen
    /// one, or a deliberately mutated copy of it.
    fn sweep_of(
        record: &Record,
        table: &'static LittleCharacterTable,
        embedding: &SubgroupEmbedding,
        partition: &FullStarPartition,
    ) -> SweepReport {
        let domains = source_domains().expect("the frozen domains");
        let domain = domains
            .get(&(record.parent_sg, table.label))
            .unwrap_or_else(|| panic!("no frozen domain for SG {} {}", record.parent_sg, table.label));
        let context = SweepContext {
            ordinal: record.ordinal,
            label: table.label,
            child_sg: record.child_sg,
            child_reciprocal: &partition.child_reciprocal,
            child_rotations: &partition.child_rotations,
            arms: &partition.arms,
        };
        let mut report = SweepReport {
            worst_score: 1.0,
            ..SweepReport::default()
        };
        sweep_label(
            &record.subgroup,
            embedding,
            table,
            partition,
            domain,
            &context,
            &mut report,
        );
        report
    }

    /// The three child-frame values `sweep_boundary_pass` needs, for a witness
    /// record.
    fn witness_frame(record: &Record) -> (SubgroupEmbedding, Lattice) {
        let embedding =
            SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup).expect("embedding");
        let reciprocal = reciprocal_lattice(record.child_sg).expect("the child lattice");
        (embedding, reciprocal)
    }

    /// **Card 6, positive control: the audit is not vacuous, and it reproduces
    /// the partition on the frozen witness.**
    ///
    /// The numbers here are the single-pair form of the corpus totals: the
    /// relation count is the closed form `n^2 r - n`, the identically-true
    /// relations are the partition's own `permanent_counts`, the recomputed
    /// boundary set is the partition's, the two interior points of each of the
    /// eight intervals hold no relation, and the arm count agrees with the
    /// orbit–stabiliser count over the parent's rotations.
    #[test]
    fn the_arm_geometry_audit_reproduces_the_frozen_partition() {
        let record = record_of(196, 10_030);
        let (embedding, partition) = sweep_context_of(&record, "DT1");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let report = sweep_of(&record, table, &embedding, &partition);
        assert!(
            report.failures.is_empty(),
            "the unmutated witness must audit cleanly: {:?}",
            report.failures
        );
        let geometry = &report.geometry;
        let arms = partition.arms.len();
        let rotations = partition.child_rotations.len();
        assert_eq!(geometry.pairs, 1);
        assert!(
            arms > 1 && rotations > 1,
            "the witness must have a real star and rotations: {arms} arm(s), {rotations} rotation(s)"
        );
        assert_eq!(geometry.relations, arms * arms * rotations - arms);
        assert_eq!(geometry.predicted_relations, geometry.relations);
        assert_eq!(
            geometry.permanent,
            partition.permanent_counts.iter().sum::<usize>()
        );
        assert!(geometry.permanent > 0, "the witness has identically-true relations");
        assert_eq!(geometry.boundaries, partition.boundary_parameters().len());
        assert!(geometry.boundaries > 0, "the witness has boundaries");
        assert_eq!(geometry.recomputed_boundaries, geometry.boundaries);
        assert_eq!(geometry.disagreements, 0);
        assert_eq!(
            geometry.evaluated,
            geometry.boundaries * (geometry.relations - geometry.permanent),
            "every boundary evaluates the whole non-permanent relation list"
        );
        assert!(geometry.evaluated > 0);
        assert_eq!(
            geometry.interior_evaluated,
            2 * report.intervals * (geometry.relations - geometry.permanent),
            "every interval evaluates the whole list at both interior points"
        );
        assert!(geometry.interior_evaluated > 0);
        assert_eq!(geometry.interior_points, 2 * report.intervals);
        assert_eq!(geometry.interior_holds, 0);
        assert_eq!(geometry.arm_orbit_checks, 1);
        assert_eq!(geometry.arm_orbit_disagreements, 0);
        assert_eq!(geometry.arms, arms);
        // The boundary pass ran at every probe parameter and answered all of
        // them, with the conservation the gate asserts.
        assert_eq!(report.boundary_points, report.pairs * geometry.boundaries);
        assert_eq!(report.boundary_successes, report.boundary_points);
        assert_eq!(report.boundary_failures, 0);
        assert_eq!(
            report.boundary_points,
            report.boundary_successes + report.boundary_failures
        );
        // The dimension sums come from the returned decompositions and match the
        // closed form the audit accumulates from the partition's own counts.
        assert_eq!(
            report.boundary_parent_dimension, report.boundary_covered_dimension,
            "every answered boundary decomposition must cover its parent dimension"
        );
        assert_eq!(
            report.boundary_parent_dimension,
            (report.boundary_points as u64) * u64::from(table.dimension) * (arms as u64)
        );
        assert_eq!(
            report.geometry.predicted_boundary_dimension, report.boundary_parent_dimension
        );
    }

    /// **Card 6: the audit's own comparisons, at unit level.**  The corpus only
    /// ever exercises them in their passing direction (every count and every set
    /// agrees), so the gate drives them with synthetic partitions; this is the
    /// same call, kept in the fast test battery so a comparison that stops
    /// working is visible without the minute-long corpus sweep.
    #[test]
    fn the_arm_geometry_audit_checks_its_own_comparisons() {
        let (cases, violations) = audit_self_check();
        assert_eq!(cases, SELF_CHECK_CASE_PIN);
        assert!(violations.is_empty(), "{violations:?}");
    }

    /// **Card 6, fault injection (a): dropping a boundary is caught by the
    /// independent recomputation.**
    ///
    /// The mutation is a synthetic data mutation on the census's own copy: the
    /// frozen partition of ordinal 10038 `DT1` is cloned and the arm-merge cut at
    /// `t = 1/8` is removed from `boundaries` (`FullStarPartition`'s fields are
    /// public; the library builder is untouched).  `probe_parameters` follows the
    /// mutated list, so the sweep's own set-equality check cannot see the drop:
    /// what sees it is the audit, which derives the boundary set from the
    /// relations themselves.
    ///
    /// The **binding** half of the card's "boundary/parameter binding *or*
    /// geometry integrity" is not reachable by mutating a `FullStarPartition`:
    /// `probe_parameters` and `boundary_parameters` are both functions of the
    /// same `boundaries` list, so any data mutation moves them together.  Its
    /// witness is therefore the library-level one the card-5 audit recorded
    /// (moving `1/8` inside `boundary_parameters`), and the duplicate-entry case
    /// is asserted here instead: a partition that lists a parameter twice is
    /// reported by the audit and turns into a zero-length interval the interval
    /// guard rejects.
    #[test]
    fn dropping_a_boundary_is_caught_by_the_geometry_audit() {
        let record = record_of(196, 10_038);
        let (embedding, partition) = sweep_context_of(&record, "DT1");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let eighth = rational(1, 8);
        assert!(
            partition.is_boundary(&eighth),
            "the witness cut must exist before the mutation"
        );
        let clean = sweep_of(&record, table, &embedding, &partition);
        assert!(
            clean.failures.is_empty(),
            "the unmutated witness must sweep cleanly: {:?}",
            clean.failures
        );
        assert_eq!(clean.geometry.disagreements, 0);

        let mut mutated = partition.clone();
        mutated
            .boundaries
            .retain(|boundary| boundary.parameter != eighth);
        let report = sweep_of(&record, table, &embedding, &mutated);
        assert!(
            report.geometry.recomputed_boundaries > report.geometry.boundaries,
            "the audit must still derive the dropped cut: {} against {}",
            report.geometry.recomputed_boundaries,
            report.geometry.boundaries
        );
        let audit = report
            .failures
            .iter()
            .find(|failure| failure.contains("independently recomputed full-star boundary set"))
            .unwrap_or_else(|| panic!("the drop must be reported: {:?}", report.failures));
        assert!(
            audit.contains("does not cut [1/8]"),
            "the failure must name the missing parameter: {audit}"
        );
        // Count conservation is untouched by the mutation: the parameters that
        // are still there were all probed and answered, so the gate's arithmetic
        // still balances and the *finding* is the audit's, not a lost count.
        let mut violations = Vec::new();
        sweep_conservation(&report, &mut violations);
        assert!(
            report.boundary_points == report.boundary_successes + report.boundary_failures,
            "a dropped boundary must not unbalance the boundary counters"
        );

        // A duplicated boundary entry is a different kind of corrupted partition
        // list: the audit reports the duplicate instead of silently accepting the
        // shortened set.
        let mut duplicated = partition.clone();
        let first = duplicated.boundaries[0].clone();
        duplicated.boundaries.push(first);
        let report = sweep_of(&record, table, &embedding, &duplicated);
        assert!(
            report
                .failures
                .iter()
                .any(|failure| failure.contains("boundary parameter(s) but only")),
            "a duplicated boundary must be reported: {:?}",
            report.failures
        );
        assert!(
            report
                .failures
                .iter()
                .any(|failure| failure.contains("is not positive")),
            "the zero-length interval it creates must be rejected: {:?}",
            report.failures
        );
    }

    /// **Card 6, fault injection (b): the dropped cut is a 6-against-4
    /// arm-grouping change.**
    ///
    /// Same mutation as (a).  Two things are measured here, and the second is
    /// what the card asks for:
    ///
    /// 1. The sweep's own interval `(0, 1/4)` of the mutated partition is
    ///    compared at its trisection points `1/12` and `1/6`.  Both are
    ///    **generic**: the merge at `1/8` is an isolated point event, so no
    ///    interval comparison can see a dropped cut.  (This is asserted, not
    ///    assumed: those are the measured block counts.)
    /// 2. The cut really is an arm-grouping change of the kind the interval check
    ///    is about: the production decompositions on either side of it -- the
    ///    generic `1/12` and the dropped `1/8` -- are compared with the sweep's
    ///    own [`match_block_geometry`], and the mismatch is the pinned `6` against
    ///    `4` block count.
    #[test]
    fn the_dropped_merge_is_a_six_against_four_arm_grouping_change() {
        let record = record_of(196, 10_038);
        let (embedding, partition) = sweep_context_of(&record, "DT1");
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let eighth = rational(1, 8);
        let twelfth = rational(1, 12);
        let mut mutated = partition.clone();
        mutated
            .boundaries
            .retain(|boundary| boundary.parameter != eighth);
        let report = sweep_of(&record, table, &embedding, &mutated);
        assert!(
            !report
                .failures
                .iter()
                .any(|failure| failure.contains("at the first parameter against")),
            "no interval comparison of the mutated sweep may see the dropped cut: {:?}",
            report.failures
        );
        assert!(report.failures.len() >= 1, "the drop is still a failure");
        assert_eq!(
            report.intervals, 7,
            "dropping one of the eight grid parameters leaves seven intervals"
        );

        let context = SweepContext {
            ordinal: record.ordinal,
            label: "DT1",
            child_sg: record.child_sg,
            child_reciprocal: &mutated.child_reciprocal,
            child_rotations: &mutated.child_rotations,
            arms: &mutated.arms,
        };
        let read = |parameter: Rat| -> Vec<SweepBlock> {
            let result =
                subduce_line_at_parameter(&record.subgroup, &embedding, table, parameter)
                    .expect("the production decomposition");
            result
                .blocks()
                .iter()
                .map(|block| {
                    SweepBlock::read(&context, parameter, block).expect("the block geometry")
                })
                .collect()
        };
        let generic = read(twelfth);
        let merged = read(eighth);
        assert_eq!(
            generic.len(),
            6,
            "the generic geometry of this line is six one-armed stars"
        );
        assert_eq!(merged.len(), 4, "the merge at 1/8 gives four stars");
        // The two trisection points of the mutated interval are both generic,
        // which is why the interval machinery alone cannot see the drop.
        assert_eq!(
            read(rational(1, 6)).len(),
            generic.len(),
            "1/6 is generic as well: the two interior points of (0, 1/4) agree"
        );
        let error = match_block_geometry(&generic, &merged, "the dropped 10038 merge")
            .expect_err("the dropped cut is an arm-grouping change");
        assert!(
            error.contains("6 block(s) at the first parameter against 4"),
            "the merge is a block-count change: {error}"
        );
    }

    /// **Card 6, fault injection (c): a block's source taken from an unrelated
    /// block of the same probe is caught by the per-block binding.**
    ///
    /// The mutation is a synthetic data mutation on the census's own report: the
    /// pinned ordinal 13688 `SM1` `t = 1/8` probe mixes stored and constructed
    /// blocks, so moving one block's `source` to the other's value is visible
    /// both in the block's own counts (`classify_counts`) and against the
    /// engine's own reading of that block.
    #[test]
    fn a_borrowed_block_source_is_caught_by_the_engine_binding() {
        let mut report = census_report(225, 13_688, false, false);
        assert!(report.failures.is_empty(), "the witness must probe cleanly");
        let eighth = rational(1, 8);
        let probe = report
            .probes
            .iter_mut()
            .find(|probe| probe.label == "SM1" && probe.parameter == eighth)
            .expect("the pinned card-4 witness probe");
        let stored = probe
            .blocks
            .iter()
            .position(|block| block.source == BlockSource::Stored)
            .expect("a stored block");
        let constructed = probe
            .blocks
            .iter()
            .position(|block| block.source != BlockSource::Stored)
            .expect("a constructed block");
        assert_ne!(stored, constructed);
        probe.blocks[stored].source = probe.blocks[constructed].source;

        let block_rows = report
            .probes
            .iter()
            .map(|probe| probe.blocks.len())
            .sum::<usize>();
        let evidence = CensusEvidence {
            child_grid: report.child_grid,
            child_algorithms: report.child_algorithms,
            child_union: &report.child_parameters,
            algorithm_checks: 0,
            algorithm_mismatches: 0,
            block_rows,
            block_file_rows: None,
            gamma: &report.gamma,
            recount: &report.recount,
            recount_ran: false,
            recount_rows_written: 0,
            recount_file_rows: None,
            sweep: &report.sweep,
            sweep_ran: false,
        };
        let mut violations = Vec::new();
        check_blocks(&report.probes, &evidence, &mut violations);
        assert!(
            violations
                .iter()
                .any(|violation| violation.contains("is recorded as")
                    && violation.contains("classify it as")),
            "the block's own counts must reject the borrowed label: {violations:?}"
        );
        assert!(
            violations
                .iter()
                .any(|violation| violation.contains("but the engine block is")),
            "the engine's own reading must reject it too: {violations:?}"
        );
        assert!(
            violations.iter().any(|violation| violation.contains("block")
                && violation.contains("1/8")),
            "the failure must name the probe: {violations:?}"
        );
    }

    /// **Card 6, fault injection (d): a decomposition that returns an error is
    /// counted and fails the gate.**
    ///
    /// Every production entry point rejects a record whose embedding belongs to
    /// another isotropy record, so pairing ordinal 10038's `DT1` subgroup with
    /// ordinal 10030's embedding makes every boundary decomposition fail.  The
    /// mutation is expressed directly on the census's own inputs (the two
    /// records are both SG 196 `DT1` records of the frozen corpus); the boundary
    /// pass is then run on its own, and the test asserts both halves of the card:
    /// the item is still **counted** (`boundary_points == successes + failures`,
    /// failures == parameters) and the **gate fails** (`sweep_conservation`
    /// reports the failed parameters, and the runtime turns the failure list into
    /// violations).
    #[test]
    fn a_failing_boundary_decomposition_is_counted_and_fails() {
        let record = record_of(196, 10_038);
        let other = record_of(196, 10_030);
        let (embedding, partition) = sweep_context_of(&record, "DT1");
        let (wrong_embedding, _) = witness_frame(&other);
        let table = line_table(record.parent_sg, "DT1").expect("the frozen table");
        let domains = source_domains().expect("the frozen domains");
        let domain = domains.get(&(196, "DT1")).expect("the DT1 domain");
        let parameters = partition.probe_parameters(domain).expect("the probe set");
        assert_eq!(parameters.len(), 8, "the eighth grid of the witness");
        // The engine error this mutation produces, read from the engine itself
        // rather than described: the failure message has to carry it.
        let engine_error = subduce_line_at_parameter(
            &record.subgroup,
            &wrong_embedding,
            table,
            parameters[0],
        )
        .expect_err("the mismatched embedding must be rejected")
        .to_string();
        let context = SweepContext {
            ordinal: record.ordinal,
            label: "DT1",
            child_sg: record.child_sg,
            child_reciprocal: &partition.child_reciprocal,
            child_rotations: &partition.child_rotations,
            arms: &partition.arms,
        };
        let mut report = SweepReport {
            worst_score: 1.0,
            ..SweepReport::default()
        };
        // The positive control on the very same call: with the record's own
        // embedding every parameter is answered.
        sweep_boundary_pass(
            &record.subgroup,
            &embedding,
            table,
            &parameters,
            &context,
            &mut report,
        );
        assert_eq!(report.boundary_successes, parameters.len());
        assert_eq!(report.boundary_failures, 0);
        assert!(report.failures.is_empty(), "{:?}", report.failures);
        // The mutation: the other record's embedding.
        let mut report = SweepReport {
            worst_score: 1.0,
            ..SweepReport::default()
        };
        sweep_boundary_pass(
            &record.subgroup,
            &wrong_embedding,
            table,
            &parameters,
            &context,
            &mut report,
        );
        assert_eq!(report.boundary_points, parameters.len());
        assert_eq!(
            report.boundary_points,
            report.boundary_successes + report.boundary_failures,
            "no parameter may be lost when the engine fails"
        );
        assert_eq!(report.boundary_successes, 0);
        assert_eq!(report.boundary_failures, parameters.len());
        assert_eq!(
            report.boundary_parent_dimension, 0,
            "a failed parameter contributes no engine dimension"
        );
        assert_eq!(report.failures.len(), parameters.len());
        for (index, failure) in report.failures.iter().enumerate() {
            assert!(failure.contains("ordinal 10038"), "{failure}");
            assert!(failure.contains("DT1"), "{failure}");
            assert!(
                failure.contains(&parameters[index].to_string()),
                "the failure must name its parameter: {failure}"
            );
            assert!(
                failure.contains(&engine_error),
                "the failure must carry the engine's error {engine_error:?}: {failure}"
            );
        }
        // The gate half: the conservation still holds (nothing was lost) and the
        // failed parameters are a violation, so a run under `--gate` or
        // `--domain-sweep` exits 1.
        let mut violations = Vec::new();
        sweep_conservation(&report, &mut violations);
        assert!(
            violations
                .iter()
                .any(|violation| violation.contains("boundary pass failed on 8 parameter(s)")),
            "the failed boundary parameters must fail the gate: {violations:?}"
        );
        assert!(
            !violations
                .iter()
                .any(|violation| violation.contains("as 0 success(es) + 8 failure(s)")),
            "the conservation itself must hold: {violations:?}"
        );
    }
}
