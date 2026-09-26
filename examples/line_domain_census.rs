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
//!   reference candidate set and classifies every block there.
//!
//! Usage:
//!
//! ```text
//! line_domain_census                 # summary + parameter partition
//! line_domain_census --gate          # exit 1 unless every invariant holds
//! line_domain_census --full-star-recount
//!                                    # rerun the gate's card-3 recount pass alone
//! line_domain_census --output out.tsv
//! line_domain_census --output-blocks blocks.tsv
//! line_domain_census --sequential    # one thread (default: rayon over records)
//! ```
//!
//! The gate deliberately accepts `unsupported > 0`: the census exists to *name*
//! the unsupported domains, not to pretend they are closed.  `--require-covered`
//! turns a non-empty unsupported set into a failure, for the round that closes
//! them all.

use cryspglib::irrep::line_monodromy::{line_direction, line_table};
use cryspglib::irrep::subduction::star::decompose::{
    FullStarBlock, FullStarError, FullStarTarget, LineSubduction, ParameterKind,
    official_line_parameter, subduce_line_at_parameter,
};
use cryspglib::irrep::subduction::star::line_domain::{
    FoldedArm, FullStarPartition, GammaParameters, ParentDomain, child_cocycle_is_a_coboundary, child_exceptional_parameters,
    full_star_partition, little_co_group_order, minimal_parameter_step,
    minimal_parameter_step_via_coordinates, parent_domain, point_cocycle_is_a_coboundary,
    point_little_co_group_order, reciprocal_lattice, require_reciprocal_direction, rotation_set,
    verify_against_grid,
};
use cryspglib::irrep::generated_data::SG_DATA_HALL;
use cryspglib::irrep::subduction::{
    Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R, fold_wave_vector,
};
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
line_domain_census [--gate] [--require-covered] [--full-star-recount] [--sequential]
                   [--output <path>] [--output-blocks <path>]

  --gate               exit 1 unless the census invariants hold (includes the
                       full-star recount pass)
  --require-covered    additionally fail when any (record, parameter) is unsupported
  --full-star-recount  run the card-3 full-star recount pass on its own, with
                       the gate's failure semantics
  --sequential         one thread (default parallel over isotropy records)
  --output <path>      write the per-probe table as TSV
  --output-blocks <path>
                       write the per-block table as TSV
";

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
fn engine_block(block: &FullStarBlock, reference_point: &Vec3R, child_reciprocal: &Lattice) -> EngineBlock {
    let mut carries_reference = false;
    let mut carries_gamma = false;
    for point in block.points() {
        carries_gamma |= child_reciprocal.contains(point.q()).unwrap_or(false);
        carries_reference |= child_reciprocal
            .same_mod(point.q(), reference_point)
            .unwrap_or(false);
    }
    EngineBlock {
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
    }
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
    let recount = gate || full_star_recount;
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

    let domains = source_domains()?;
    let records = records()?;
    let official = official_line_parameter().map_err(|error| error.to_string())?;
    let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1).map_err(|e| e.to_string())?;
    let cache = PointClassCache::default();

    let per_record: Vec<RecordReport> = if sequential {
        records
            .iter()
            .map(|record| probe_record(record, &domains, official, generic, recount, &cache))
            .collect()
    } else {
        records
            .par_iter()
            .map(|record| probe_record(record, &domains, official, generic, recount, &cache))
            .collect()
    };
    let mut probes: Vec<Probe> = Vec::new();
    let mut probe_errors: Vec<String> = Vec::new();
    let mut child_grid = (0usize, 0usize);
    let mut child_algorithms = (0usize, 0usize);
    let mut child_union: Vec<Rat> = Vec::new();
    let mut gamma = GammaReport::default();
    let mut recount_report = RecountReport::default();
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
    violations.extend(algorithm_failures);
    violations.extend(block_file_failures);
    violations.extend(recount_report.failures.iter().cloned());
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

    if gate || require_covered || full_star_recount {
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

/// Every isotropy record that carries parametric-k rows, with its labels.
fn records() -> Result<Vec<Record>, String> {
    let mut out = Vec::new();
    for sg in 1..=230u8 {
        out.extend(records_of(sg)?);
    }
    Ok(out)
}

/// The same walk, restricted to one parent space group.
///
/// Split out so the card-4 regression tests can build one record without
/// enumerating the whole corpus first; both callers go through the same code, so
/// a test cannot exercise a record the census would not have built.
fn records_of(sg: u8) -> Result<Vec<Record>, String> {
    let mut out = Vec::new();
    for record in query::irreps_of(sg) {
        if record.spinor || record.subgroups().is_empty() {
            continue;
        }
        let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
            .map_err(|error| format!("SG {sg} {}: {error}", record.ml))?;
        for subgroup in subgroups {
            let Ok(rows) = subgroup.other_wave_vector_subduction() else {
                continue;
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
                child_sg: u8::try_from(subgroup.record.sg).unwrap_or(0),
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
        if context.child_reciprocal.contains(point.q()).unwrap_or(false) {
            carries_gamma = true;
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
        usize::try_from(parent / dimension).unwrap_or(0)
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
    cache: &PointClassCache,
) -> RecordReport {
    let mut out = Vec::new();
    let mut failures = Vec::new();
    let mut child_grid = (0usize, 0usize);
    let mut child_algorithms = (0usize, 0usize);
    let mut child_parameters: Vec<Rat> = Vec::new();
    let mut gamma = GammaReport::default();
    let mut recount_report = RecountReport::default();
    let frame = probe_frame(record, &mut failures);
    let Some((embedding, child_reciprocal, child_rotations)) = frame else {
        return RecordReport {
            probes: out,
            child_grid,
            child_algorithms: (0, 0),
            child_parameters,
            gamma,
            recount: recount_report,
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
                    if !left.contains(&vector).unwrap_or(false) {
                        failures.push(format!(
                            "ordinal {}: the {name} child lattice is not the embedding's \
                             reciprocal lattice (row {row})",
                            record.ordinal
                        ));
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
        let Ok((parameters, child_candidates)) =
            record_parameters(record, domain, table, official, generic)
        else {
            failures.push(format!(
                "ordinal {}: parameter partition failed for {label}",
                record.ordinal
            ));
            continue;
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
        if recount {
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
                    let engine_blocks: Vec<EngineBlock> = result
                        .blocks()
                        .iter()
                        .map(|block| {
                            engine_block(
                                block,
                                &reference_point_for_engine,
                                block_context.child_reciprocal,
                            )
                        })
                        .collect();
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
    let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1);
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
    fn record_of(parent_sg: u8, ordinal: usize) -> Record {
        records_of(parent_sg)
            .expect("the parent's records")
            .into_iter()
            .find(|record| record.ordinal == ordinal)
            .unwrap_or_else(|| panic!("SG {parent_sg} has no parametric-k record {ordinal}"))
    }

    /// The census report of one record, exactly as `run` builds it.
    fn census_report(parent_sg: u8, ordinal: usize, recount: bool) -> RecordReport {
        let domains = source_domains().expect("the frozen domains");
        let official = official_line_parameter().expect("the official parameter");
        let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1).expect("the generic sample");
        let cache = PointClassCache::default();
        let record = record_of(parent_sg, ordinal);
        probe_record(&record, &domains, official, generic, recount, &cache)
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
        let report = census_report(225, 13_688, false);
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
        let report = census_report(196, 10_038, true);
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
}
