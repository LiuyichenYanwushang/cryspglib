//! Task 9: exhaustive ordinary isotropy / identity-subduction audit.
//!
//! This example is a reproducible census of the pinned ordinary isotropy table.
//! It reuses the production embedding (`SubgroupEmbedding::from_isotropy_subgroup`)
//! and the production scalar full-star decomposition
//! (`subduce_full_star_with_embedding`) and never re-implements a solver.
//!
//! ```text
//! cargo run --release -p cryspglib --example audit_irrep_subduction -- --output target/task9/full.tsv
//! cargo run --release -p cryspglib --example audit_irrep_subduction -- --parent 221
//! cargo run --release -p cryspglib --example audit_irrep_subduction -- --ordinal 13345
//! ```
//!
//! Contract of the report:
//!
//! * Every non-spinor (`ordinary scalar`) record of an in-scope parent space
//!   group is a probe of every condensate record.  Whenever the isotropy record
//!   embeds, the probe is decomposed by the production full-star engine -- even
//!   when the exact folded-star geometry already proves that the trivial child
//!   irrep cannot appear.  Geometry is an independent check on the computed
//!   result, never a reason to skip the decomposition.
//! * When that full decomposition has no child data for a folded star, the probe
//!   is *not* written off: the engine's `trivial_content_with_embedding` answers
//!   the identity content alone, exactly, by skipping the stars whose wave
//!   vector is not the child's Gamma point (they provably cannot carry the
//!   trivial representation; see that function's contract).  Such rows are
//!   `identity_only` and are compared with the stored table exactly like a full
//!   result; the row's detail keeps the full-decomposition error verbatim.
//! * One TSV row is emitted per scalar probe, including probes whose embedding
//!   is unavailable.  Such rows are `uncomputed` and can never count as passes.
//! * The probe partition is
//!   `full_success + identity_only + missing + error + uncomputed`, with
//!   `uncomputed = uncomputed_embedding + uncomputed_no_trivial`; it is checked
//!   against `probes_total` both per record and globally.
//! * The positive stored entries (the pinned 94271 identity rows) form a
//!   separate partition: `passed + mismatch + embedding_unavailable +
//!   no_trivial + missing + error == unique stored rows`, while duplicates and
//!   unresolvable labels are reported separately.  Rows are never summed: a
//!   repeated pair with the same frequency is deduplicated (`duplicate_same`)
//!   and a repeated pair with different frequencies is a conflict
//!   (`duplicate_conflict`), never an addition.
//! * `full_success` requires the engine's own dimension, integrality and
//!   reconstruction guarantees plus this example's independent check that every
//!   reported target has a unique frozen CIR source; a result that fails any of
//!   them is an `error` (exit 1), never a pass.
//! * Missing child data is `missing` (incomplete) and genuine engine
//!   inconsistency is `error` (exit 1); neither is ignored or downgraded.
//! * `other_wave_vector_subduction` rows (5756 over 1006 records) are a
//!   separate single-valued table naming parent irreps **at other wave
//!   vectors**.  Each row is resolved here against the frozen source list
//!   (`IRREP_W_LABELS`/`IRREP_W_SPACE_GROUP`, 73 entries, reached independently
//!   of the accessor that produced the row) and its space group must be the
//!   record's parent; a failure of either half is a hard failure.  Their
//!   frequencies cannot be computed from the pinned data as it stands: the
//!   irrep table stores them as label, space group, dimension and type only, and
//!   the archived little table shows all 73 to be **parameterized line** wave
//!   vectors `k = Gamma + t*v` (`scripts/check_other_wave_vector_rows.py`
//!   decodes and asserts that), so there is no single numeric k to fold.  They
//!   are reported in their own `w_scope` line and their own completeness counter
//!   (`w_incomplete`), which
//!   `--require-complete` deliberately does not gate on and
//!   `--require-w-complete` does.  Identity closure never substitutes for their
//!   missing parameters.
//! * spinors are a separate unsupported space: spinor records are counted and
//!   any spinor reference in these tables is reported, never computed.
//! * the Frobenius check runs only on condensates at `Gamma` whose stored rows
//!   are all at `Gamma` and whose complex component dimensions come from the
//!   row's own source metadata.  `ConjugateRealification` rows and the pinned
//!   `DistinctComponentSum` rows (whose two component dimensions are equal)
//!   contribute `frequency * sum(component_dims) / component_count`; a row with
//!   unequal component dimensions is explicitly unevaluated.  The independent
//!   target is `|G| / |H|` from the parent point-group order and the embedding's
//!   coset representatives.
//!
//! Production checks (the API's own dimension, integrality, tiling and
//! reconstruction checks plus the source-availability cross-check of every
//! reported CIR target) are reported separately from the independent evidence
//! (stored-frequency comparison, folded-star geometry zeros and Frobenius
//! reciprocity).  The identity comparison checks multiplicities only; it does
//! not validate non-trivial label correctness, and this example does not claim
//! otherwise.
//!
//! Exit codes: `0` clean, `1` a genuine mismatch/inconsistency, `2`
//! `--require-complete` with incomplete ordinary-table coverage (a probe with
//! neither a full result nor an exact identity content, an uncompared stored
//! row, or an unevaluated geometry/Frobenius gate), or `--require-w-complete`
//! with other-wave-vector rows left uncomputed, `3` CLI/IO failure.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::process::ExitCode;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;

use cryspglib::irrep::isotropy::{self, IsotropySubgroup, parent_primitive_basis, subgroup_size};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, FullStarSubduction, subduce_full_star_with_embedding,
    trivial_content_with_embedding,
};
use cryspglib::irrep::subduction::star::scalar_star::ScalarStar;
use cryspglib::irrep::subduction::{
    Lattice, Mat3R, Rat, SubductionComponent, SubductionError, SubgroupEmbedding,
};
use cryspglib::irrep::types::{
    CompoundCharacterSemantics, IrrepRecord, IrrepSourceIdentity, KVector,
};
use cryspglib::irrep::{LabelConvention, generated_data};
use rayon::prelude::*;

const EXPECTED_CONDENSATES: usize = 15_239;
const EXPECTED_IDENTITY_ROWS: usize = 94_271;
const EXPECTED_OTHER_WAVE_ROWS: usize = 5_756;
const EXPECTED_OTHER_WAVE_RECORDS: usize = 1_006;
const EXPECTED_SPINOR_RECORDS: usize = 3_611;
const EXPECTED_DISTINCT_PROBES: usize = 4_777;

const HEADER: &str = "kind\tordinal\tparent_sg\tsubgroup_sg\tdirection\tdomain\tarms\tsize\t\
probe_ml\tprobe_sg\tprobe_k\tprobe_src\tstored\tcomputed\tstatus\tdetail";

const RECONSTRUCTION_TOLERANCE: f64 = 1e-7;
const MISMATCH_PRINT_LIMIT: usize = 40;

// ── Options ──────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Options {
    parent: Option<u8>,
    ordinal: Option<usize>,
    output: Option<String>,
    require_complete: bool,
    require_w_complete: bool,
    /// Require a *full decomposition* for every ordinary probe in scope, not
    /// just an exact identity content.
    require_full_decomposition: bool,
    progress: usize,
    help: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            parent: None,
            ordinal: None,
            output: None,
            require_complete: false,
            require_w_complete: false,
            require_full_decomposition: false,
            progress: 500,
            help: false,
        }
    }
}

impl Options {
    fn parse(args: impl Iterator<Item = String>) -> Result<Self, String> {
        let mut options = Options::default();
        let mut args = args;
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--parent" => {
                    let value = args.next().ok_or("--parent needs a space group number")?;
                    let parent: u8 = value
                        .parse()
                        .map_err(|_| format!("invalid --parent value {value}"))?;
                    if parent == 0 || parent > 230 {
                        return Err(format!("--parent {parent} is outside 1-230"));
                    }
                    options.parent = Some(parent);
                }
                "--ordinal" => {
                    let value = args.next().ok_or("--ordinal needs a record index")?;
                    let ordinal: usize = value
                        .parse()
                        .map_err(|_| format!("invalid --ordinal value {value}"))?;
                    if ordinal >= EXPECTED_CONDENSATES {
                        return Err(format!("--ordinal {ordinal} is outside 0..15239"));
                    }
                    options.ordinal = Some(ordinal);
                }
                "--output" => {
                    options.output = Some(args.next().ok_or("--output needs a path")?.to_string());
                }
                "--require-complete" => options.require_complete = true,
                "--require-w-complete" => options.require_w_complete = true,
                "--require-full-decomposition" => options.require_full_decomposition = true,
                "--progress" => {
                    let value = args.next().ok_or("--progress needs a record count")?;
                    options.progress = value
                        .parse()
                        .map_err(|_| format!("invalid --progress value {value}"))?;
                }
                "--help" | "-h" => options.help = true,
                other => return Err(format!("unknown argument {other}")),
            }
        }
        Ok(options)
    }

    fn scoped(&self) -> bool {
        self.parent.is_some() || self.ordinal.is_some()
    }
}

const USAGE: &str = "\
usage: audit_irrep_subduction [--parent N] [--ordinal N] [--output PATH]
                              [--require-complete] [--require-full-decomposition]
                              [--require-w-complete] [--progress N]

  --parent N          audit only the condensate records of space group N (1-230)
  --ordinal N         audit only the isotropy record with this global ordinal
  --output PATH       write the machine-readable TSV to PATH (default stdout)
  --require-complete  exit 2 unless the ordinary identity-subduction table is
                      complete in scope: every embedding is frozen, every stored
                      positive row is reproduced, every other probe has an exact
                      result (full decomposition or exact identity content),
                      and no geometry/Frobenius check is left unevaluated
  --require-full-decomposition
                      exit 2 unless every ordinary probe in scope has a *full*
                      decomposition.  An exact identity-only content is not
                      enough: that category stays incomplete here even though
                      --require-complete accepts it.  A scoped run reports a
                      scoped result only and never claims global coverage
  --require-w-complete
                      additionally require the 5756 other-wave-vector rows to be
                      computed.  This flag is the explicit gate for that separate
                      question; see the `w_scope` summary line for its state
  --progress N        print a progress line to stderr every N records (0 off)

The TSV goes to stdout (or --output); the terse summary goes to stderr.

Exit codes: 1 = a demonstrated inconsistency (computation error, frequency
conflict, accounting violation); 2 = a requested coverage gate is not met;
0 = every requested gate passed within the reported scope.";

// ── Accounting partitions ────────────────────────────────────────────────────

/// Catalog of one decomposition attempt.  The two tables are disjoint: a probe
/// that could not be attempted is `uncomputed`, never a success.
#[derive(Debug, Default, Clone)]
struct ProbeCounts {
    total: usize,
    full_success: usize,
    /// Probes the engine answered exactly for the identity content alone: the
    /// full decomposition has no data for a folded star that provably cannot
    /// contribute to that content.
    identity_only: usize,
    missing: usize,
    error: usize,
    uncomputed_embedding: usize,
    uncomputed_no_trivial: usize,
}

macro_rules! add_count_fields {
    ($target:ident; $($field:ident),+ $(,)?) => {
        $($target.$field += *$field;)+
    };
}

impl ProbeCounts {
    /// Probes handed to the production engine.
    fn attempted(&self) -> usize {
        self.full_success + self.identity_only + self.missing + self.error
    }

    /// Probes with a valid embedding that the engine could not answer at all:
    /// neither a full decomposition nor an exact identity content.
    fn incomplete(&self) -> usize {
        self.missing
    }

    /// Probes no engine call was possible for.
    fn uncomputed(&self) -> usize {
        self.uncomputed_embedding + self.uncomputed_no_trivial
    }

    fn partition_sum(&self) -> usize {
        self.attempted() + self.uncomputed()
    }

    fn partition_error(&self) -> Option<String> {
        if self.partition_sum() != self.total {
            return Some(format!(
                "probe partition {}+{}+{}+{}+{}+{} = {} != probes_total {}",
                self.full_success,
                self.identity_only,
                self.missing,
                self.error,
                self.uncomputed_embedding,
                self.uncomputed_no_trivial,
                self.partition_sum(),
                self.total
            ));
        }
        None
    }

    fn add(&mut self, other: &ProbeCounts) {
        let ProbeCounts {
            total,
            full_success,
            identity_only,
            missing,
            error,
            uncomputed_embedding,
            uncomputed_no_trivial,
        } = other;
        add_count_fields!(
            self;
            total,
            full_success,
            identity_only,
            missing,
            error,
            uncomputed_embedding,
            uncomputed_no_trivial,
        );
    }
}

/// Catalog of the positive stored identity entries (the pinned 94271 rows).
#[derive(Debug, Default, Clone)]
struct EntryCounts {
    entries: usize,
    unique: usize,
    unresolved: usize,
    duplicate_same: usize,
    duplicate_conflict: usize,
    passed: usize,
    mismatch: usize,
    embedding_unavailable: usize,
    no_trivial: usize,
    missing: usize,
    error: usize,
}

impl EntryCounts {
    fn compared(&self) -> usize {
        self.passed + self.mismatch
    }

    fn unique_partition_sum(&self) -> usize {
        self.passed
            + self.mismatch
            + self.embedding_unavailable
            + self.no_trivial
            + self.missing
            + self.error
    }

    fn partition_error(&self) -> Option<String> {
        if self.unique_partition_sum() != self.unique {
            return Some(format!(
                "positive-entry partition {}+{}+{}+{}+{}+{} = {} != unique stored rows {}",
                self.passed,
                self.mismatch,
                self.embedding_unavailable,
                self.no_trivial,
                self.missing,
                self.error,
                self.unique_partition_sum(),
                self.unique
            ));
        }
        let total_sum =
            self.unique + self.unresolved + self.duplicate_same + self.duplicate_conflict;
        if total_sum != self.entries {
            return Some(format!(
                "stored entry partition {}+{}+{}+{} = {} != stored entries {}",
                self.unique,
                self.unresolved,
                self.duplicate_same,
                self.duplicate_conflict,
                total_sum,
                self.entries
            ));
        }
        None
    }

    fn incomplete(&self) -> usize {
        self.embedding_unavailable
            + self.no_trivial
            + self.missing
            + self.unresolved
            + self.duplicate_same
    }

    fn hard(&self) -> usize {
        self.mismatch + self.error + self.duplicate_conflict
    }

    fn add(&mut self, other: &EntryCounts) {
        let EntryCounts {
            entries,
            unique,
            unresolved,
            duplicate_same,
            duplicate_conflict,
            passed,
            mismatch,
            embedding_unavailable,
            no_trivial,
            missing,
            error,
        } = other;
        add_count_fields!(
            self;
            entries,
            unique,
            unresolved,
            duplicate_same,
            duplicate_conflict,
            passed,
            mismatch,
            embedding_unavailable,
            no_trivial,
            missing,
            error,
        );
    }
}

/// Per-record tally; absorbed into [`Counts`] after the record is complete so
/// every emission path updates one object and the partitions tile exactly.
#[derive(Debug, Default, Clone)]
struct RecordTally {
    records: usize,
    records_embedding_ok: usize,
    records_embedding_failed: usize,
    records_no_trivial: usize,
    probes: ProbeCounts,
    entries: EntryCounts,
    absent_zero_engine: usize,
    absent_zero_geometry: usize,
    absent_positive: usize,
    absent_attempted: usize,
    geometry_reject_absent: usize,
    geometry_eligible_absent: usize,
    geometry_filter_errors: usize,
    geometry_contradictions: usize,
    basis_errors: usize,
}

/// Global accounting.  `probes` and `entries` are two independent partitions:
/// every scalar probe and every positive stored entry lands in exactly one
/// category of its own partition.
#[derive(Debug, Default, Clone)]
struct Counts {
    records: usize,
    records_embedding_ok: usize,
    records_embedding_failed: usize,
    records_no_trivial: usize,
    probes: ProbeCounts,
    entries: EntryCounts,
    absent_zero_engine: usize,
    absent_zero_geometry: usize,
    absent_positive: usize,
    absent_attempted: usize,
    geometry_reject_absent: usize,
    geometry_eligible_absent: usize,
    geometry_filter_errors: usize,
    geometry_contradictions: usize,
    basis_errors: usize,
    probe_rows_emitted: usize,
    stored_k_mismatch: usize,
    stored_domain_out_of_range: usize,
    target_source_unmatched: usize,
    production_dim_mismatch: usize,
    production_integrality_mismatch: usize,
    production_recon_mismatch: usize,
    label_source_disagreement: usize,
    frobenius_records: usize,
    frobenius_pass: usize,
    frobenius_mismatch: usize,
    frobenius_unevaluated: usize,
    frobenius_strict: usize,
    frobenius_realification: usize,
    frobenius_compound: usize,
    gamma_record_non_gamma_probe: usize,
    w_records: usize,
    w_entries: usize,
    w_rows_emitted: usize,
    /// Rows whose `(label, space group)` pair is one of the 73 frozen
    /// other-wave-vector sources and whose space group is the record's parent.
    w_source_resolved: usize,
    /// Rows that fail either half of that resolution: a real inconsistency in
    /// the pinned pairing, not a missing-parameter question.
    w_source_mismatch: usize,
    /// Resolved rows whose source has a frozen little-group character table
    /// (`w_little_characters_data`), i.e. rows the engine can compute once the
    /// arm character source is wired up.
    w_character_frozen: usize,
    /// Resolved rows whose source character table is still unresolved, so the
    /// row stays reported instead of computed.
    w_character_blocked: usize,
    w_frequency_mismatch: usize,
    /// Transport failures: a row whose **complete** decomposition at `t = 5/4`
    /// differs from the decomposition its monodromy image `M_v(label)` has at
    /// `t = 1/4`.
    ///
    /// Every frozen direction is a parent reciprocal lattice vector, so the two
    /// parameters describe the same band at the same point of the parent's zone
    /// and the step is the relabelling `(k + K, M_K(alpha)) ~ (k, alpha)`, not a
    /// gauge (`docs/subduction-conventions.md` §16).  R6.1 compared the *same*
    /// label at both parameters -- the over-strong `M == 1` reading -- and
    /// canonicalized the wave vector to silence the 40 rows (SG 210/227/228) that
    /// disagreed.
    w_parameter_shift_mismatch: usize,
    w_parameter_shift_checked: usize,
    /// Rows whose monodromy image the frozen fingerprints do not determine, so
    /// the transport comparison is skipped instead of guessed.
    w_parameter_shift_skipped: usize,
    /// Rows where the revoked same-label reading answers differently from the
    /// transported one: the measured non-vacuity of the gate.
    w_parameter_shift_witnesses: usize,
    /// Block-route calls that returned `Err`: an arithmetic/context failure in
    /// the computation itself, not a missing pinned input.  A demonstrated
    /// error must fail the run under every flag combination.
    w_engine_error: usize,
    w_conflict: usize,
    w_computed: usize,
    spinor_records: usize,
    spinor_probe_rows: usize,
    accounting_violations: usize,
    census_mismatch: usize,
}

/// Run classification; the exit code is derived from it, never from identity
/// closure alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    Clean,
    Incomplete,
    Inconsistent,
}

/// The three independent completeness gates.  Each answers a different
/// question and they are never merged into one boolean:
///
/// * `complete` -- the ordinary identity-subduction table is closed in scope
///   (a full decomposition *or* an exact identity content answers every probe);
/// * `full_decomposition` -- every ordinary probe in scope has a full
///   decomposition, so an identity-only result is a gap here;
/// * `w_complete` -- the parameterized other-wave-vector rows are computed.
///
/// A hard failure (computation error, frequency conflict, accounting
/// violation) outranks every gate: the run is inconsistent, exit 1.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Gates {
    complete: bool,
    full_decomposition: bool,
    w_complete: bool,
}

impl Gates {
    fn from_options(options: &Options) -> Self {
        Self {
            complete: options.require_complete,
            full_decomposition: options.require_full_decomposition,
            w_complete: options.require_w_complete,
        }
    }

    fn any(self) -> bool {
        self.complete || self.full_decomposition || self.w_complete
    }

    fn label(self) -> String {
        let mut names = Vec::new();
        if self.complete {
            names.push("--require-complete");
        }
        if self.full_decomposition {
            names.push("--require-full-decomposition");
        }
        if self.w_complete {
            names.push("--require-w-complete");
        }
        if names.is_empty() {
            "none".to_string()
        } else {
            names.join(",")
        }
    }
}

fn line_image_skip_reason(
    image: Option<&cryspglib::irrep::line_monodromy::LabelImage>,
) -> (&'static str, &'static str) {
    use cryspglib::irrep::line_monodromy::LabelImage;

    match image {
        Some(LabelImage::UnsupportedShift) => (
            "computed_transport_skipped:unsupported_shift",
            "exact little-group rotation criterion for the phase twist is not met",
        ),
        Some(LabelImage::Ambiguous(_)) => (
            "computed_transport_skipped:ambiguous_image",
            "frozen fingerprints give multiple image labels",
        ),
        Some(LabelImage::Missing) | None => (
            "computed_transport_skipped:missing_image",
            "frozen data has no image label",
        ),
        Some(LabelImage::Unique(_)) => {
            unreachable!("line_image_skip_reason called for a unique image")
        }
    }
}

impl Counts {
    fn add(&mut self, other: &Counts) {
        let Counts {
            records,
            records_embedding_ok,
            records_embedding_failed,
            records_no_trivial,
            probes,
            entries,
            absent_zero_engine,
            absent_zero_geometry,
            absent_positive,
            absent_attempted,
            geometry_reject_absent,
            geometry_eligible_absent,
            geometry_filter_errors,
            geometry_contradictions,
            basis_errors,
            probe_rows_emitted,
            stored_k_mismatch,
            stored_domain_out_of_range,
            target_source_unmatched,
            production_dim_mismatch,
            production_integrality_mismatch,
            production_recon_mismatch,
            label_source_disagreement,
            frobenius_records,
            frobenius_pass,
            frobenius_mismatch,
            frobenius_unevaluated,
            frobenius_strict,
            frobenius_realification,
            frobenius_compound,
            gamma_record_non_gamma_probe,
            w_records,
            w_entries,
            w_rows_emitted,
            w_source_resolved,
            w_source_mismatch,
            w_character_frozen,
            w_character_blocked,
            w_frequency_mismatch,
            w_parameter_shift_mismatch,
            w_parameter_shift_checked,
            w_parameter_shift_skipped,
            w_parameter_shift_witnesses,
            w_engine_error,
            w_conflict,
            w_computed,
            spinor_records,
            spinor_probe_rows,
            accounting_violations,
            census_mismatch,
        } = other;
        self.probes.add(probes);
        self.entries.add(entries);
        add_count_fields!(
            self;
            records,
            records_embedding_ok,
            records_embedding_failed,
            records_no_trivial,
            absent_zero_engine,
            absent_zero_geometry,
            absent_positive,
            absent_attempted,
            geometry_reject_absent,
            geometry_eligible_absent,
            geometry_filter_errors,
            geometry_contradictions,
            basis_errors,
            probe_rows_emitted,
            stored_k_mismatch,
            stored_domain_out_of_range,
            target_source_unmatched,
            production_dim_mismatch,
            production_integrality_mismatch,
            production_recon_mismatch,
            label_source_disagreement,
            frobenius_records,
            frobenius_pass,
            frobenius_mismatch,
            frobenius_unevaluated,
            frobenius_strict,
            frobenius_realification,
            frobenius_compound,
            gamma_record_non_gamma_probe,
            w_records,
            w_entries,
            w_rows_emitted,
            w_source_resolved,
            w_source_mismatch,
            w_character_frozen,
            w_character_blocked,
            w_frequency_mismatch,
            w_parameter_shift_mismatch,
            w_parameter_shift_checked,
            w_parameter_shift_skipped,
            w_parameter_shift_witnesses,
            w_engine_error,
            w_conflict,
            w_computed,
            spinor_records,
            spinor_probe_rows,
            accounting_violations,
            census_mismatch,
        );
    }

    fn absorb(&mut self, tally: &RecordTally) {
        self.records += tally.records;
        self.records_embedding_ok += tally.records_embedding_ok;
        self.records_embedding_failed += tally.records_embedding_failed;
        self.records_no_trivial += tally.records_no_trivial;
        self.probes.add(&tally.probes);
        self.entries.add(&tally.entries);
        self.absent_zero_engine += tally.absent_zero_engine;
        self.absent_zero_geometry += tally.absent_zero_geometry;
        self.absent_positive += tally.absent_positive;
        self.absent_attempted += tally.absent_attempted;
        self.geometry_reject_absent += tally.geometry_reject_absent;
        self.geometry_eligible_absent += tally.geometry_eligible_absent;
        self.geometry_filter_errors += tally.geometry_filter_errors;
        self.geometry_contradictions += tally.geometry_contradictions;
        self.basis_errors += tally.basis_errors;
    }

    /// The row and candidate partitions must tile exactly.  Any mismatch means
    /// a row was dropped, double counted or summed twice.
    fn accounting_error(&self) -> Option<String> {
        if self.probe_rows_emitted != self.probes.total {
            return Some(format!(
                "emitted probe rows {} != probes_total {}",
                self.probe_rows_emitted, self.probes.total
            ));
        }
        if self.w_rows_emitted != self.w_entries {
            return Some(format!(
                "emitted other-wave-vector rows {} != w_entries {}",
                self.w_rows_emitted, self.w_entries
            ));
        }
        self.probes
            .partition_error()
            .or_else(|| self.entries.partition_error())
    }

    fn absent_zero(&self) -> usize {
        self.absent_zero_engine + self.absent_zero_geometry
    }

    fn positive_comparisons(&self) -> usize {
        self.entries.compared()
    }

    /// Categories that make the ordinary identity-subduction table incomplete
    /// in scope: a probe with neither a full result nor an exact identity
    /// content, a stored row left uncompared, an unchecked absent-zero probe,
    /// or an unevaluated geometry/Frobenius gate.
    ///
    /// The other-wave-vector rows are *not* part of this scope: their k vectors
    /// and character rows are absent from the pinned archive, so no engine can
    /// compute them from it.  They are counted by [`Counts::w_incomplete`] and
    /// gated separately by `--require-w-complete`.
    fn incomplete(&self) -> usize {
        self.probes.incomplete()
            + self.probes.uncomputed()
            + self.entries.incomplete()
            + self.geometry_filter_errors
            + self.basis_errors
            + self.frobenius_unevaluated
    }

    /// Other-wave-vector rows left uncomputed, with the reason reported beside
    /// them: the pinned archive carries no k vector and no character row for
    /// the 73 irreps these rows name.
    fn w_incomplete(&self) -> usize {
        self.w_entries.saturating_sub(self.w_computed)
    }

    /// Categories that leave the *full decomposition* of the ordinary table
    /// incomplete in scope.  Everything [`Counts::incomplete`] already counts,
    /// plus the identity-only results: those answer the trivial multiplicity
    /// exactly but are not a full decomposition, so the stronger gate must not
    /// accept them.
    fn full_decomposition_incomplete(&self) -> usize {
        self.incomplete() + self.probes.identity_only
    }

    /// Categories that are genuine inconsistencies and fail the default run.
    fn hard_failures(&self) -> usize {
        self.probes.error
            + self.entries.hard()
            + self.absent_positive
            + self.geometry_contradictions
            + self.basis_errors
            + self.stored_k_mismatch
            + self.stored_domain_out_of_range
            + self.gamma_record_non_gamma_probe
            + self.target_source_unmatched
            + self.production_dim_mismatch
            + self.production_integrality_mismatch
            + self.production_recon_mismatch
            + self.label_source_disagreement
            + self.frobenius_mismatch
            + self.w_conflict
            + self.w_source_mismatch
            + self.w_frequency_mismatch
            + self.w_parameter_shift_mismatch
            + self.w_parameter_shift_skipped
            + self.w_engine_error
            + self.accounting_violations
            + self.census_mismatch
    }

    fn verdict(&self, gates: Gates) -> Verdict {
        if self.hard_failures() > 0 {
            Verdict::Inconsistent
        } else if (gates.complete && self.incomplete() > 0)
            || (gates.full_decomposition && self.full_decomposition_incomplete() > 0)
            || (gates.w_complete && self.w_incomplete() > 0)
        {
            Verdict::Incomplete
        } else {
            Verdict::Clean
        }
    }

    fn exit_code(&self, gates: Gates) -> u8 {
        match self.verdict(gates) {
            Verdict::Clean => 0,
            Verdict::Inconsistent => 1,
            Verdict::Incomplete => 2,
        }
    }
}

// ── Small data types ─────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy)]
struct StoredRow {
    frequency: u16,
    domain: u16,
    ml: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Geometry {
    /// The folded star has a point in the child reciprocal lattice modulo the
    /// child lattice: the trivial child irrep may appear.
    Eligible,
    /// The folded star has no such point: the trivial child irrep cannot
    /// appear by exact lattice arithmetic.
    Reject,
    /// The independent geometry check could not run.
    Unknown,
    /// The probe has a stored positive entry, so the zero-candidate geometry
    /// filter does not apply; the engine result is compared to the stored
    /// frequency instead.
    NotChecked,
}

impl Geometry {
    fn label(self) -> &'static str {
        match self {
            Geometry::Eligible => "eligible",
            Geometry::Reject => "reject",
            Geometry::Unknown => "unknown",
            Geometry::NotChecked => "not_checked",
        }
    }
}

#[derive(Debug, Default, Clone)]
struct ChildSources {
    ordinary: HashMap<u32, u8>,
    constituent: HashMap<u32, u8>,
}

impl ChildSources {
    fn build(sg: u8) -> Self {
        let mut sources = ChildSources::default();
        for record in query::irreps_of(sg) {
            match record.source_identity() {
                IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => {
                    // The target dimension is the selected-arm (little-group)
                    // complex dimension, not the full-star row dimension.
                    if let Ok(row) = record.ordinary_scalar_selected_arm_block_trace()
                        && let Ok(dimension) = u8::try_from(row.dimension())
                    {
                        sources.ordinary.insert(cir_irnumber, dimension);
                    }
                }
                IrrepSourceIdentity::Compound { .. } => {
                    if let Some(metadata) = record.compound_metadata() {
                        for (index, irnumber) in metadata.cir_irnumbers.iter().enumerate() {
                            sources
                                .constituent
                                .insert(*irnumber, metadata.cir_dimensions[index]);
                        }
                    }
                }
                IrrepSourceIdentity::Spin { .. } => {}
            }
        }
        sources
    }

    fn has(&self, irnumber: u32, dimension: u8) -> bool {
        self.ordinary.get(&irnumber) == Some(&dimension)
            || self.constituent.get(&irnumber) == Some(&dimension)
    }
}

const SPACE_GROUP_CACHE_SIZE: usize = u8::MAX as usize + 1;

fn cache_cells<T>() -> Vec<OnceLock<T>> {
    (0..SPACE_GROUP_CACHE_SIZE)
        .map(|_| OnceLock::new())
        .collect()
}

struct SharedChildCaches {
    trivial: Vec<OnceLock<Option<&'static IrrepRecord>>>,
    sources: Vec<OnceLock<ChildSources>>,
    reciprocal: Vec<OnceLock<Option<Lattice>>>,
}

impl SharedChildCaches {
    fn new() -> Self {
        Self {
            trivial: cache_cells(),
            sources: cache_cells(),
            reciprocal: cache_cells(),
        }
    }

    fn child_sources(&self, sg: u8) -> &ChildSources {
        self.sources[usize::from(sg)].get_or_init(|| ChildSources::build(sg))
    }
}

static SHARED_CHILD_CACHES: OnceLock<SharedChildCaches> = OnceLock::new();

fn shared_child_caches() -> &'static SharedChildCaches {
    SHARED_CHILD_CACHES.get_or_init(SharedChildCaches::new)
}

#[derive(Debug, Default, Clone)]
struct CallOutcome {
    total: u32,
    by_label: u32,
    targets: usize,
    targets_without_source: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CallErrorKind {
    MissingChildData,
    MissingIrrepData,
    Hard,
}

impl CallErrorKind {
    fn is_missing(self) -> bool {
        matches!(
            self,
            CallErrorKind::MissingChildData | CallErrorKind::MissingIrrepData
        )
    }

    fn label(self) -> &'static str {
        match self {
            CallErrorKind::MissingChildData => "missing_child_data",
            CallErrorKind::MissingIrrepData => "missing_irrep_data",
            CallErrorKind::Hard => "hard_error",
        }
    }
}

// ── Auditor ──────────────────────────────────────────────────────────────────

#[derive(Clone, Default)]
struct ReportBuffer(Arc<Mutex<Vec<u8>>>);

impl ReportBuffer {
    fn contents(&self) -> Vec<u8> {
        self.0.lock().expect("report buffer lock").clone()
    }
}

impl Write for ReportBuffer {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        self.0
            .lock()
            .map_err(|_| io::Error::other("report buffer lock poisoned"))?
            .extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct AuditPart {
    counts: Counts,
    error_detail: BTreeMap<String, usize>,
    mismatches: Vec<String>,
    anomalies: Vec<String>,
    distinct_probes: HashSet<(u8, u16)>,
    output: Vec<u8>,
}

struct ProgressReporter {
    completed: std::sync::atomic::AtomicUsize,
    printed_through: Mutex<usize>,
    interval: usize,
    started: Instant,
}

impl ProgressReporter {
    fn new(interval: usize, started: Instant) -> Self {
        Self {
            completed: std::sync::atomic::AtomicUsize::new(0),
            printed_through: Mutex::new(0),
            interval,
            started,
        }
    }

    fn record_complete(&self) {
        use std::sync::atomic::Ordering;

        let completed = self.completed.fetch_add(1, Ordering::Relaxed) + 1;
        if self.interval == 0 || !completed.is_multiple_of(self.interval) {
            return;
        }
        let Ok(mut printed_through) = self.printed_through.lock() else {
            return;
        };
        let current = self.completed.load(Ordering::Relaxed);
        let mut next = *printed_through + self.interval;
        while next <= current {
            eprintln!(
                "progress records={next} elapsed={:.1}s",
                self.started.elapsed().as_secs_f64()
            );
            *printed_through = next;
            next += self.interval;
        }
    }
}

struct Auditor {
    options: Options,
    counts: Counts,
    error_detail: BTreeMap<String, usize>,
    mismatches: Vec<String>,
    anomalies: Vec<String>,
    distinct_probes: HashSet<(u8, u16)>,
    order_cache: HashMap<u8, usize>,
    writer: Box<dyn Write>,
    started: Instant,
    progress: Option<Arc<ProgressReporter>>,
    check_census: bool,
    io_error: Option<String>,
}

impl Auditor {
    fn new(options: Options, writer: Box<dyn Write>) -> Self {
        Self {
            options,
            counts: Counts::default(),
            error_detail: BTreeMap::new(),
            mismatches: Vec::new(),
            anomalies: Vec::new(),
            distinct_probes: HashSet::new(),
            order_cache: HashMap::new(),
            writer,
            started: Instant::now(),
            progress: None,
            check_census: true,
            io_error: None,
        }
    }

    fn bump_error(&mut self, key: &str) {
        *self.error_detail.entry(key.to_string()).or_insert(0) += 1;
    }

    fn note_w_transport_skipped(&mut self, ordinal: usize, label: &str, reason: &str) {
        self.counts.w_parameter_shift_skipped += 1;
        self.mismatch(format!(
            "ordinal {ordinal}: other-wave-vector row {label} transport comparison skipped: {reason}"
        ));
    }

    fn anomaly(&mut self, message: String) {
        if self.anomalies.len() < MISMATCH_PRINT_LIMIT {
            self.anomalies.push(message.clone());
        }
        eprintln!("ANOMALY {message}");
    }

    fn mismatch(&mut self, message: String) {
        if self.mismatches.len() < MISMATCH_PRINT_LIMIT {
            self.mismatches.push(message.clone());
        }
        eprintln!("MISMATCH {message}");
    }

    fn absorb_part(&mut self, part: AuditPart) {
        self.counts.add(&part.counts);
        for (key, value) in part.error_detail {
            *self.error_detail.entry(key).or_default() += value;
        }
        self.anomalies.extend(
            part.anomalies
                .into_iter()
                .take(MISMATCH_PRINT_LIMIT.saturating_sub(self.anomalies.len())),
        );
        self.mismatches.extend(
            part.mismatches
                .into_iter()
                .take(MISMATCH_PRINT_LIMIT.saturating_sub(self.mismatches.len())),
        );
        self.distinct_probes.extend(part.distinct_probes);

        if self.io_error.is_none()
            && let Err(error) = self.writer.write_all(&part.output)
        {
            self.io_error = Some(error.to_string());
        }
    }

    fn finish_part(mut self, output: Vec<u8>) -> Result<AuditPart, String> {
        if let Some(error) = self.io_error.take() {
            return Err(format!("writing the report failed: {error}"));
        }
        self.writer
            .flush()
            .map_err(|error| format!("flushing the report failed: {error}"))?;
        Ok(AuditPart {
            counts: self.counts,
            error_detail: self.error_detail,
            mismatches: self.mismatches,
            anomalies: self.anomalies,
            distinct_probes: self.distinct_probes,
            output,
        })
    }

    fn emit(&mut self, line: String) {
        if self.io_error.is_some() {
            return;
        }
        if let Err(error) = writeln!(self.writer, "{line}") {
            self.io_error = Some(error.to_string());
        }
    }

    /// One row per scalar probe, counted so `accounting_error` can prove that
    /// no probe was silently dropped from the report.
    #[allow(clippy::too_many_arguments)]
    fn emit_probe(
        &mut self,
        context: &RecordContext<'_>,
        kind: &str,
        probe: &IrrepRecord,
        stored: Option<StoredRow>,
        computed: Option<&CallOutcome>,
        status: &str,
        detail: &str,
    ) {
        self.counts.probe_rows_emitted += 1;
        self.emit(identity_line(
            context, kind, probe, stored, computed, status, detail,
        ));
    }

    /// Compare one computed trivial content with the stored table, update the
    /// per-record tallies, and return the status of the emitted row.
    ///
    /// `success_status` separates a full decomposition from the identity-only
    /// result of a probe whose non-Gamma folding star has no stored data.
    #[allow(clippy::too_many_arguments)]
    fn compare_trivial_content(
        &mut self,
        tally: &mut RecordTally,
        context: &RecordContext<'_>,
        probe: &IrrepRecord,
        stored_row: Option<StoredRow>,
        total: u32,
        geometry: Geometry,
        success_status: &'static str,
    ) -> &'static str {
        let ordinal = context.ordinal;
        let sg = context.sg;
        let child_sg = context.child_sg;
        let direction = context.subgroup.record.direction_label;
        let status = if let Some(row) = stored_row {
            if total == u32::from(row.frequency) {
                tally.entries.passed += 1;
                success_status
            } else {
                tally.entries.mismatch += 1;
                self.mismatch(format!(
                    "ordinal {ordinal} parent {sg} subgroup {child_sg} direction {direction} \
                     probe {}: stored frequency {} but the decomposition gives {total}",
                    probe.ml, row.frequency
                ));
                "mismatch"
            }
        } else if total == 0 {
            if geometry == Geometry::Reject {
                tally.absent_zero_geometry += 1;
            } else {
                tally.absent_zero_engine += 1;
            }
            "absent_zero"
        } else {
            tally.absent_positive += 1;
            self.mismatch(format!(
                "ordinal {ordinal} parent {sg} subgroup {child_sg} direction {direction} \
                 probe {}: not listed but the decomposition gives {total}",
                probe.ml
            ));
            "absent_positive"
        };
        if geometry == Geometry::Reject && total > 0 {
            tally.geometry_contradictions += 1;
            self.mismatch(format!(
                "ordinal {ordinal} probe {}: geometry proves zero but the decomposition gives \
                 {total} trivial terms",
                probe.ml
            ));
        }
        status
    }

    fn trivial_child(&mut self, sg: u8) -> Option<&'static IrrepRecord> {
        shared_child_caches()
            .trivial
            .get(usize::from(sg))?
            .get_or_init(|| find_trivial_child(sg))
            .as_ref()
            .copied()
    }

    fn child_reciprocal(&mut self, sg: u8) -> Option<Lattice> {
        shared_child_caches()
            .reciprocal
            .get(usize::from(sg))?
            .get_or_init(|| build_child_reciprocal(sg).ok())
            .as_ref()
            .copied()
    }

    fn parent_point_group_order(&mut self, sg: u8) -> Option<usize> {
        if let Some(cached) = self.order_cache.get(&sg) {
            return Some(*cached);
        }
        let operations = query::symmetry_operations_of(sg).ok()?;
        let rotations: HashSet<[[i32; 3]; 3]> = operations
            .iter()
            .map(|operation| operation.rotation)
            .collect();
        let order = rotations.len();
        self.order_cache.insert(sg, order);
        Some(order)
    }

    fn audit_sg(&mut self, sg: u8) -> Result<(), String> {
        if let Some(parent) = self.options.parent
            && parent != sg
        {
            return Ok(());
        }
        let records = query::irreps_of(sg);
        if records.is_empty() {
            return Err(format!("space group {sg} has no irrep records"));
        }
        self.counts.spinor_records += records.iter().filter(|record| record.spinor).count();

        // Wrap the stored isotropy records once per irrep, then keep only the
        // records selected by --ordinal; out-of-scope space groups never build
        // stars or touch the embedding engine.
        let mut wrapped: Vec<(usize, Vec<IsotropySubgroup>)> = Vec::new();
        for (record_index, record) in records.iter().enumerate() {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
                .map_err(|error| format!("space group {sg} irrep {}: {error}", record.ml))?;
            if subgroups.len() != record.subgroups().len() {
                self.anomaly(format!(
                    "space group {sg} irrep {}: {} wrapped isotropy records for {} stored",
                    record.ml,
                    subgroups.len(),
                    record.subgroups().len()
                ));
                continue;
            }
            wrapped.push((record_index, subgroups));
        }
        if let Some(ordinal) = self.options.ordinal {
            wrapped.retain(|(_, subgroups)| {
                subgroups.iter().any(|subgroup| subgroup.ordinal == ordinal)
            });
        }
        if wrapped.is_empty() {
            return Ok(());
        }
        let probes: Vec<&'static IrrepRecord> =
            records.iter().filter(|record| !record.spinor).collect();
        let mut probe_by_ml: HashMap<&'static str, usize> = HashMap::new();
        for (index, probe) in probes.iter().enumerate() {
            if probe_by_ml.insert(probe.ml, index).is_some() {
                self.anomaly(format!(
                    "space group {sg} has more than one record labelled {}; labels are ambiguous",
                    probe.ml
                ));
            }
        }
        let stars: Vec<Result<ScalarStar, String>> = probes
            .iter()
            .map(|probe| ScalarStar::new(probe).map_err(|error| error.to_string()))
            .collect();

        for (record_index, subgroups) in &wrapped {
            for subgroup in subgroups {
                if let Some(ordinal) = self.options.ordinal
                    && subgroup.ordinal != ordinal
                {
                    continue;
                }
                self.audit_record(sg, *record_index, subgroup, &probes, &probe_by_ml, &stars)?;
                if let Some(progress) = &self.progress {
                    progress.record_complete();
                }
            }
        }
        Ok(())
    }

    fn audit_record(
        &mut self,
        sg: u8,
        record_index: usize,
        subgroup: &IsotropySubgroup,
        probes: &[&'static IrrepRecord],
        probe_by_ml: &HashMap<&'static str, usize>,
        stars: &[Result<ScalarStar, String>],
    ) -> Result<(), String> {
        let ordinal = subgroup.ordinal;
        let child_sg = u8::try_from(subgroup.record.sg).map_err(|_| {
            format!(
                "ordinal {ordinal}: subgroup number {} does not fit u8",
                subgroup.record.sg
            )
        })?;
        let record = &query::irreps_of(sg)[record_index];
        let stored = isotropy::identity_subduction(ordinal)
            .map_err(|error| format!("ordinal {ordinal}: identity subduction failed: {error}"))?;
        let wave_entries = isotropy::other_wave_vector_subduction(ordinal).map_err(|error| {
            format!("ordinal {ordinal}: other-wave-vector subduction failed: {error}")
        })?;

        let mut tally = RecordTally {
            records: 1,
            ..RecordTally::default()
        };
        tally.probes.total = probes.len();
        // The stored basis determinant is already computed in i128 by the
        // engine; keep the full error (including its i128 payload) instead of
        // silently narrowing or dropping it.
        let size = match subgroup_size(subgroup.record.basis) {
            Ok(size) => Some(size),
            Err(error) => {
                tally.basis_errors += 1;
                self.anomaly(format!("ordinal {ordinal}: subgroup size: {error}"));
                None
            }
        };
        let ctx = RecordContext {
            ordinal,
            sg,
            child_sg,
            subgroup,
            size,
        };

        // ── Repeat / domain audit before any row is treated as a term ────────
        let mut stored_map: BTreeMap<usize, StoredRow> = BTreeMap::new();
        for entry in &stored {
            tally.entries.entries += 1;
            if entry.parent_sg != sg {
                tally.entries.unresolved += 1;
                self.anomaly(format!(
                    "ordinal {ordinal}: stored row parent space group {} != {sg}",
                    entry.parent_sg
                ));
                continue;
            }
            let Some(&probe_index) = probe_by_ml.get(entry.parent_ml) else {
                tally.entries.unresolved += 1;
                let spinor = query::irreps_of(sg)
                    .iter()
                    .any(|candidate| candidate.ml == entry.parent_ml && candidate.spinor);
                if spinor {
                    self.counts.spinor_probe_rows += 1;
                    self.bump_error("stored_probe_spinor_unsupported");
                } else {
                    self.bump_error("stored_probe_label_unresolved");
                    self.anomaly(format!(
                        "ordinal {ordinal}: stored probe {} is not a record of space group {sg}",
                        entry.parent_ml
                    ));
                }
                continue;
            };
            let probe = probes[probe_index];
            if usize::from(entry.domain) > subgroup.record.domains {
                self.counts.stored_domain_out_of_range += 1;
                self.mismatch(format!(
                    "ordinal {ordinal}: stored probe {} has domain {} but the record declares {} domains",
                    probe.ml, entry.domain, subgroup.record.domains
                ));
            }
            if !same_k(probe.k_vector(), entry.parent_k) {
                self.counts.stored_k_mismatch += 1;
                self.mismatch(format!(
                    "ordinal {ordinal}: probe {} stored k {} but the record has k {}",
                    probe.ml,
                    format_k(entry.parent_k),
                    format_k(probe.k_vector())
                ));
            }
            let row = StoredRow {
                frequency: entry.frequency,
                domain: entry.domain,
                ml: entry.parent_ml,
            };
            if let Err(message) =
                insert_stored_row(&mut stored_map, probe_index, row, &mut tally.entries)
            {
                self.mismatch(format!("ordinal {ordinal}: {message}"));
            }
        }
        tally.entries.unique = stored_map.len();
        for probe_index in stored_map.keys() {
            self.distinct_probes.insert((sg, *probe_index as u16));
        }
        for entry in &stored {
            if entry.parent_k.numerators != [0, 0, 0] {
                self.counts.gamma_record_non_gamma_probe += usize::from(is_gamma(record));
            }
        }

        // ── Embedding ────────────────────────────────────────────────────────
        let embedding = match SubgroupEmbedding::from_isotropy_subgroup(subgroup) {
            Ok(embedding) => {
                tally.records_embedding_ok = 1;
                Some(embedding)
            }
            Err(error) => {
                tally.records_embedding_failed = 1;
                self.bump_error(&format!("embedding:{}", embedding_label(&error)));
                None
            }
        };
        let trivial = if embedding.is_some() {
            match self.trivial_child(child_sg) {
                Some(record) => Some(record),
                None => {
                    tally.records_no_trivial = 1;
                    self.bump_error("child_trivial_record_missing");
                    None
                }
            }
        } else {
            None
        };
        let trivial_cir = trivial.and_then(|record| match record.source_identity() {
            IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => Some(cir_irnumber),
            _ => None,
        });
        let subgroup_order = embedding
            .as_ref()
            .map(|embedding| embedding.representatives().len());

        // ── Every scalar probe: geometry evidence plus the production call ───
        let mut computed: Vec<Option<u32>> = vec![None; probes.len()];
        for (probe_index, probe) in probes.iter().enumerate() {
            let stored_row = stored_map.get(&probe_index).copied();
            let kind = if stored_row.is_some() {
                "identity"
            } else {
                "absent"
            };
            let (Some(embedding), Some(trivial), Some(trivial_cir)) =
                (embedding.as_ref(), trivial, trivial_cir)
            else {
                // No engine call is possible.  This probe is `uncomputed` and
                // still gets a row: identity closure elsewhere is no evidence
                // for it.
                let (status, detail) = if embedding.is_none() {
                    tally.probes.uncomputed_embedding += 1;
                    if stored_row.is_some() {
                        tally.entries.embedding_unavailable += 1;
                    }
                    (
                        "uncomputed_embedding_unavailable",
                        "uncomputed: no valid embedding for this isotropy record".to_string(),
                    )
                } else {
                    tally.probes.uncomputed_no_trivial += 1;
                    if stored_row.is_some() {
                        tally.entries.no_trivial += 1;
                    }
                    (
                        "uncomputed_no_trivial_child",
                        "uncomputed: the child table has no unique trivial Gamma record"
                            .to_string(),
                    )
                };
                self.emit_probe(&ctx, kind, probe, stored_row, None, status, &detail);
                continue;
            };

            // Exact folded-star geometry for absent candidates: independent
            // evidence, but never a gate on the engine call.
            let geometry = if stored_row.is_none() {
                tally.absent_attempted += 1;
                match (&stars[probe_index], self.child_reciprocal(child_sg)) {
                    (Ok(star), Some(reciprocal)) => {
                        match probe_is_gamma_eligible(star, embedding, &reciprocal) {
                            Ok(true) => {
                                tally.geometry_eligible_absent += 1;
                                Geometry::Eligible
                            }
                            Ok(false) => {
                                tally.geometry_reject_absent += 1;
                                Geometry::Reject
                            }
                            Err(error) => {
                                tally.geometry_filter_errors += 1;
                                self.bump_error("geometry_filter_error");
                                self.anomaly(format!(
                                    "ordinal {ordinal}, probe {}: geometry filter failed: {error}",
                                    probe.ml
                                ));
                                Geometry::Unknown
                            }
                        }
                    }
                    _ => {
                        tally.geometry_filter_errors += 1;
                        Geometry::Unknown
                    }
                }
            } else {
                Geometry::NotChecked
            };

            match subduce_full_star_with_embedding(subgroup, embedding, probe) {
                Ok(result) => match self.inspect_result(&result, trivial, trivial_cir, child_sg) {
                    Ok(outcome) => {
                        tally.probes.full_success += 1;
                        computed[probe_index] = Some(outcome.total);
                        let detail = production_detail(&outcome, geometry);
                        let status = self.compare_trivial_content(
                            &mut tally,
                            &ctx,
                            probe,
                            stored_row,
                            outcome.total,
                            geometry,
                            "passed",
                        );
                        self.emit_probe(
                            &ctx,
                            kind,
                            probe,
                            stored_row,
                            Some(&outcome),
                            status,
                            &detail,
                        );
                    }
                    Err(detail) => {
                        tally.probes.error += 1;
                        if stored_row.is_some() {
                            tally.entries.error += 1;
                        }
                        self.mismatch(format!(
                            "ordinal {ordinal} probe {}: engine result failed its own checks: {detail}",
                            probe.ml
                        ));
                        self.emit_probe(
                            &ctx,
                            kind,
                            probe,
                            stored_row,
                            None,
                            "engine_self_check_failed",
                            &format!("geometry={} {detail}", geometry.label()),
                        );
                    }
                },
                Err(error) => {
                    let call_kind = classify_call_error(&error);
                    // A folded child star without stored data cannot carry the
                    // subgroup's trivial representation (its wave vector is not
                    // the child's Gamma point), so the identity content the
                    // pinned table stores stays computable without it.  The
                    // full decomposition is impossible, and the row says so.
                    if call_kind.is_missing()
                        && let Ok(content) =
                            trivial_content_with_embedding(subgroup, embedding, probe)
                    {
                        tally.probes.identity_only += 1;
                        computed[probe_index] = Some(content.total);
                        let detail = format!(
                            "geometry={} identity_only gamma_stars={} skipped_stars={} \
                             trivial_by_label={} full_decomposition={error}",
                            geometry.label(),
                            content.gamma_stars,
                            content.skipped_stars,
                            content.by_label,
                        );
                        let status = self.compare_trivial_content(
                            &mut tally,
                            &ctx,
                            probe,
                            stored_row,
                            content.total,
                            geometry,
                            "identity_only",
                        );
                        let outcome = CallOutcome {
                            total: content.total,
                            by_label: content.by_label,
                            targets: 0,
                            targets_without_source: 0,
                        };
                        self.emit_probe(
                            &ctx,
                            kind,
                            probe,
                            stored_row,
                            Some(&outcome),
                            status,
                            &detail,
                        );
                        continue;
                    }
                    self.bump_error(&format!("engine:{}", call_kind.label()));
                    let detail = format!("geometry={} {error}", geometry.label());
                    let status = if call_kind.is_missing() {
                        tally.probes.missing += 1;
                        if stored_row.is_some() {
                            tally.entries.missing += 1;
                        }
                        "uncomputed_missing_data"
                    } else {
                        tally.probes.error += 1;
                        if stored_row.is_some() {
                            tally.entries.error += 1;
                        }
                        self.mismatch(format!(
                            "ordinal {ordinal} probe {}: engine error: {error}",
                            probe.ml
                        ));
                        "engine_error"
                    };
                    self.emit_probe(&ctx, kind, probe, stored_row, None, status, &detail);
                }
            }
        }

        // ── Record line, partitions, Frobenius, w rows ───────────────────────
        self.counts.absorb(&tally);
        self.check_record_partition(ordinal, &tally);
        let embedding_status = if tally.records_embedding_ok == 1 {
            "embedding_ok"
        } else {
            "embedding_failed"
        };
        let direction = subgroup.record.direction_label;
        let origin = subgroup.record.origin;
        self.emit(format!(
            "record\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{}\t{}\t\t\t\t\t{}\t{}\t{}\t\
probes={} full_success={} missing={} error={} uncomputed={} (embedding={} no_trivial={}) \
stored={} passed={} mismatch={} absent_zero={} geometry_reject_absent={} domains={} arms={} \
origin={},{},{},{}",
            subgroup.record.arms,
            size.map_or_else(|| "unset".to_string(), |value| value.to_string()),
            tally.probes.full_success,
            tally.entries.passed,
            embedding_status,
            tally.probes.total,
            tally.probes.full_success,
            tally.probes.missing,
            tally.probes.error,
            tally.probes.uncomputed(),
            tally.probes.uncomputed_embedding,
            tally.probes.uncomputed_no_trivial,
            tally.entries.unique,
            tally.entries.passed,
            tally.entries.mismatch,
            tally.absent_zero_engine + tally.absent_zero_geometry,
            tally.geometry_reject_absent,
            subgroup.record.domains,
            subgroup.record.arms,
            origin[0],
            origin[1],
            origin[2],
            origin[3],
        ));

        // ── Frobenius reciprocity for complete Gamma character rows ──────────
        if is_gamma(record) {
            self.frobenius(
                ordinal,
                sg,
                child_sg,
                subgroup,
                probes,
                &stored_map,
                &computed,
                size,
                subgroup_order,
            );
        }

        // ── Other-wave-vector rows: frozen source resolution only ────────────
        //
        // These rows name parent irreps **at other wave vectors** (the `DT`/`SM`
        // star labels of the cubic groups).  The pinned *irrep* table stores those
        // 73 irreps as label, space group, dimension and type only -- no k vector
        // and no character row -- and the pinned little table shows every one of
        // them to sit on a line `k = Gamma + t*v` rather than a point
        // (`scripts/check_other_wave_vector_rows.py` decodes and asserts it), so
        // there is no numeric k to fold.  The
        // audit therefore resolves each row against the frozen source list
        // itself (independent of the accessor that produced the entry) and says
        // so per row; `w_scope` in the summary reports the outcome separately
        // from the ordinary table's completeness.
        if !wave_entries.is_empty() {
            self.counts.w_records += 1;
            let mut seen: BTreeMap<&'static str, u16> = BTreeMap::new();
            for entry in &wave_entries {
                self.counts.w_entries += 1;
                let mut w_status = "uncomputed_line_star_open".to_string();
                let parent_sg_match = usize::from(entry.parent_sg) == usize::from(sg);
                let frozen = generated_data::IRREP_W_LABELS
                    .iter()
                    .zip(generated_data::IRREP_W_SPACE_GROUP.iter())
                    .any(|(label, source_sg)| {
                        *label == entry.parent_ml && usize::from(*source_sg) == usize::from(sg)
                    });
                if frozen && parent_sg_match {
                    self.counts.w_source_resolved += 1;
                    let table = cryspglib::irrep::w_little_characters_data::W_LITTLE_CHARACTERS
                        .iter()
                        .find(|table| {
                            usize::from(table.space_group) == usize::from(entry.parent_sg)
                                && table.label == entry.parent_ml
                        });
                    match table {
                        Some(table) => {
                            self.counts.w_character_frozen += 1;
                            // R6.1: the **complete** decomposition at the
                            // official parameter is *the* computation; the pinned
                            // value is only compared against its trivial content.
                            // Selecting an algorithm by its agreement with the
                            // expected answer would hide exactly the errors this
                            // audit exists to find.
                            if let Some(embedding) = embedding.as_ref() {
                                let parameter = cryspglib::irrep::subduction::star::decompose::
                                    official_line_parameter()
                                    .map_err(|error| error.to_string())?;
                                let decomposition =
                                    cryspglib::irrep::subduction::star::decompose::
                                        subduce_line_at_parameter(
                                            subgroup, embedding, table, parameter,
                                        );
                                match decomposition {
                                    Ok(result) => {
                                        let computed = result.trivial_content();
                                        match computed {
                                            Ok(value) if value == u32::from(entry.frequency) => {
                                                self.counts.w_computed += 1;
                                                w_status = "computed".to_string();
                                                // R6.2 monodromy transport, at audit
                                                // scale: `t = 5/4` is one parent
                                                // reciprocal lattice vector along the
                                                // frozen direction, so the band the row
                                                // names is the one whose frozen table is
                                                // the *monodromy image* `M_v(label)`.
                                                // The decomposition at `t = 5/4` must
                                                // therefore equal the decomposition of the
                                                // image at `t = 1/4`.  Comparing the same
                                                // label at both parameters (the revoked
                                                // R6.1 reading) is what produced the 40
                                                // SG 210/227/228 rows and the canonical wave
                                                // vector that silenced them.
                                                let shift =
                                                    match cryspglib::irrep::line_monodromy::
                                                        line_direction(table)
                                                    {
                                                        Some(direction) => Some(direction),
                                                        None => {
                                                            w_status =
                                                                "computed_transport_skipped:unparsable_direction"
                                                                    .to_string();
                                                            self.note_w_transport_skipped(
                                                                ordinal,
                                                                entry.parent_ml,
                                                                "unparsable frozen direction",
                                                            );
                                                            None
                                                        }
                                                    };
                                                if let Some(shift) = shift {
                                                    let map =
                                                        cryspglib::irrep::line_monodromy::monodromy(
                                                            entry.parent_sg,
                                                            &shift,
                                                        );
                                                    match map {
                                                        Err(error) => {
                                                            w_status =
                                                                "computed_transport_skipped:invalid_reciprocal_shift"
                                                                    .to_string();
                                                            self.note_w_transport_skipped(
                                                                ordinal,
                                                                entry.parent_ml,
                                                                &format!("invalid reciprocal shift: {error}"),
                                                            );
                                                        }
                                                        Ok(map) => match map
                                                            .unique_image(entry.parent_ml)
                                                        {
                                                            None => {
                                                                let (status, reason) =
                                                                    line_image_skip_reason(
                                                                        map.image(entry.parent_ml),
                                                                    );
                                                                w_status = status.to_string();
                                                                self.note_w_transport_skipped(
                                                                    ordinal,
                                                                    entry.parent_ml,
                                                                    reason,
                                                                );
                                                            }
                                                            Some(image) => {
                                                                let image_table =
                                                            cryspglib::irrep::line_monodromy::
                                                                line_table(
                                                                    entry.parent_sg,
                                                                    image,
                                                                );
                                                                let shifted =
                                                            cryspglib::irrep::subduction::Rat::new(
                                                                5, 4,
                                                            )
                                                            .map_err(|error| error.to_string())?;
                                                                // The errors are narrowed to `String`
                                                                // right here: `FullStarError` is large
                                                                // and must not travel through a closure.
                                                                let transported =
                                                            cryspglib::irrep::subduction::star::decompose::
                                                                subduce_line_at_parameter(
                                                                    subgroup,
                                                                    embedding,
                                                                    table,
                                                                    shifted,
                                                                )
                                                                .map_err(|error| error.to_string());
                                                                let reference = image_table.map(|image_table| {
                                                            cryspglib::irrep::subduction::star::decompose::
                                                                subduce_line_at_parameter(
                                                                    subgroup,
                                                                    embedding,
                                                                    image_table,
                                                                    parameter,
                                                                )
                                                                .map_err(|error| error.to_string())
                                                        });
                                                                match (transported, reference) {
                                                                    (
                                                                        Ok(shifted_result),
                                                                        Some(Ok(reference)),
                                                                    ) => {
                                                                        self.counts
                                                                        .w_parameter_shift_checked += 1;
                                                                        if line_decomposition_key(
                                                                        &shifted_result,
                                                                    ) != line_decomposition_key(
                                                                        &reference,
                                                                    ) {
                                                                        self.counts
                                                                        .w_parameter_shift_mismatch += 1;
                                                                        self.mismatch(format!(
                                                                        "ordinal {ordinal}: \
                                                                         other-wave-vector row {} \
                                                                         at t = 5/4 does not \
                                                                         transport to its monodromy \
                                                                         image {image} at t = 1/4",
                                                                        entry.parent_ml
                                                                    ));
                                                                    }
                                                                        // The revoked reading, measured
                                                                        // rather than assumed: where it
                                                                        // differs, the gate above is not
                                                                        // vacuous.
                                                                        if line_decomposition_key(
                                                                        &result,
                                                                    ) != line_decomposition_key(
                                                                        &shifted_result,
                                                                    ) {
                                                                        self.counts
                                                                        .w_parameter_shift_witnesses += 1;
                                                                    }
                                                                    }
                                                                    (Ok(_), None) => {
                                                                        w_status =
                                                                        "computed_transport_skipped:missing_image_table"
                                                                            .to_string();
                                                                        self.note_w_transport_skipped(
                                                                        ordinal,
                                                                        entry.parent_ml,
                                                                        "monodromy image has no frozen character table",
                                                                    );
                                                                    }
                                                                    (Ok(_), Some(Err(error))) => {
                                                                        self.counts
                                                                    .w_parameter_shift_mismatch += 1;
                                                                        self.bump_error(&format!(
                                                                            "w-line-image:{error}"
                                                                        ));
                                                                        self.mismatch(format!(
                                                                            "ordinal {ordinal}: \
                                                                     monodromy image {image} of \
                                                                     {} has no decomposition at \
                                                                     t = 1/4: {error}",
                                                                            entry.parent_ml
                                                                        ));
                                                                    }
                                                                    (Err(error), _) => {
                                                                        self.counts
                                                                    .w_parameter_shift_mismatch += 1;
                                                                        self.bump_error(&format!(
                                                                            "w-line-shift:{error}"
                                                                        ));
                                                                        self.mismatch(format!(
                                                                            "ordinal {ordinal}: \
                                                                     other-wave-vector row {} \
                                                                     engine error at t = 5/4: \
                                                                     {error}",
                                                                            entry.parent_ml
                                                                        ));
                                                                    }
                                                                }
                                                            }
                                                        },
                                                    }
                                                }
                                            }
                                            Ok(value) => {
                                                self.counts.w_frequency_mismatch += 1;
                                                w_status = format!("computed_mismatch:{value}");
                                                self.mismatch(format!(
                                                    "ordinal {ordinal}: other-wave-vector row {} \
                                                     computed {value}, pinned {}",
                                                    entry.parent_ml, entry.frequency
                                                ));
                                            }
                                            Err(error) => {
                                                self.counts.w_engine_error += 1;
                                                w_status = format!("blocks_err:{error}");
                                                self.bump_error(&format!("w-line:{error}"));
                                                self.mismatch(format!(
                                                    "ordinal {ordinal}: other-wave-vector row {} \
                                                     engine error: {error}",
                                                    entry.parent_ml
                                                ));
                                            }
                                        }
                                    }
                                    Err(error) => {
                                        // A `Err` from the block route is a
                                        // computation failure, not a missing
                                        // pinned input: it is a hard failure and
                                        // is reported per row as well as in the
                                        // engine-error detail map.
                                        self.counts.w_engine_error += 1;
                                        w_status = format!("blocks_err:{error}");
                                        self.bump_error(&format!("w-line:{error}"));
                                        self.mismatch(format!(
                                            "ordinal {ordinal}: other-wave-vector row {} engine \
                                             error: {error}",
                                            entry.parent_ml
                                        ));
                                    }
                                }
                            }
                        }
                        None => self.counts.w_character_blocked += 1,
                    }
                } else {
                    self.counts.w_source_mismatch += 1;
                    self.mismatch(format!(
                        "ordinal {ordinal}: other-wave-vector row {} of space group {} is not a \
                         frozen source of parent space group {sg}",
                        entry.parent_ml, entry.parent_sg
                    ));
                }
                if entry.frequency == 0 {
                    self.counts.w_source_mismatch += 1;
                    self.mismatch(format!(
                        "ordinal {ordinal}: other-wave-vector row {} carries frequency 0",
                        entry.parent_ml
                    ));
                }
                match seen.get(entry.parent_ml) {
                    Some(frequency) if *frequency != entry.frequency => {
                        self.counts.w_conflict += 1;
                        self.mismatch(format!(
                            "ordinal {ordinal}: other-wave-vector label {} appears with both {} and {}",
                            entry.parent_ml, frequency, entry.frequency
                        ));
                    }
                    Some(_) => {}
                    None => {
                        seen.insert(entry.parent_ml, entry.frequency);
                    }
                }
                self.counts.w_rows_emitted += 1;
                self.emit(format!(
                    "w_entry\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{}\t{}\t{}\t{sg}\t\tother_wave_vector\t{}\t\t{w_status}\tstored_parent_sg={} parent_sg_match={parent_sg_match} frozen_source={frozen}",
                    subgroup.record.arms,
                    size.map_or_else(|| "unset".to_string(), |value| value.to_string()),
                    entry.parent_ml,
                    entry.frequency,
                    entry.parent_sg,
                ));
            }
        }
        Ok(())
    }

    /// Cross-check one production result.  Success means the engine's dimension,
    /// integrality/tiling and reconstruction guarantees hold and every reported
    /// target resolves to its frozen CIR source.
    fn inspect_result(
        &mut self,
        result: &FullStarSubduction,
        trivial: &'static IrrepRecord,
        trivial_cir: u32,
        child_sg: u8,
    ) -> Result<CallOutcome, String> {
        let mut outcome = CallOutcome::default();
        {
            let sources = shared_child_caches().child_sources(child_sg);
            for block in result.blocks() {
                let term_dimension: u32 = block
                    .targets()
                    .iter()
                    .map(|target| u32::from(target.dimension) * target.multiplicity)
                    .sum();
                if term_dimension != block.little_dimension() {
                    self.counts.production_integrality_mismatch += 1;
                    return Err(format!(
                        "block terms carry {term_dimension} but its little dimension is {}",
                        block.little_dimension()
                    ));
                }
                let carried = u64::from(block.little_dimension()) * block.star_size() as u64;
                if carried != u64::from(block.block_dimension()) {
                    self.counts.production_integrality_mismatch += 1;
                    return Err(format!(
                        "little dimension {} x star size {} = {carried} != block dimension {}",
                        block.little_dimension(),
                        block.star_size(),
                        block.block_dimension()
                    ));
                }
                outcome.by_label += block.multiplicity(trivial.ml);
                for target in block.targets() {
                    outcome.targets += 1;
                    // A constructed target has no frozen CIR source: its
                    // identity is the exact point and the constructed little
                    // group, checked below, not a stored row.
                    match target.irnumber {
                        Some(irnumber) => {
                            if !sources.has(irnumber, target.dimension) {
                                outcome.targets_without_source += 1;
                            }
                        }
                        None => {
                            if !matches!(target.component, SubductionComponent::Constructed { .. })
                            {
                                outcome.targets_without_source += 1;
                            }
                        }
                    }
                    if target.dimension == trivial.dim && target.irnumber == Some(trivial_cir) {
                        outcome.total += target.multiplicity;
                    }
                }
            }
        }
        if result.covered_dimension() != result.parent_dimension() {
            self.counts.production_dim_mismatch += 1;
            return Err(format!(
                "blocks carry {} of the parent full-star dimension {}",
                result.covered_dimension(),
                result.parent_dimension()
            ));
        }
        let (parent, reconstructed) = result.reconstruction();
        if parent.len() != reconstructed.len()
            || parent
                .iter()
                .zip(reconstructed.iter())
                .any(|(left, right)| (left - right).norm() > RECONSTRUCTION_TOLERANCE)
        {
            self.counts.production_recon_mismatch += 1;
            return Err("the reported child irreps do not reconstruct the parent star".to_string());
        }
        if outcome.targets_without_source > 0 {
            self.counts.target_source_unmatched += 1;
            return Err(format!(
                "{} reported targets have no frozen CIR source",
                outcome.targets_without_source
            ));
        }
        if outcome.by_label != outcome.total {
            self.counts.label_source_disagreement += 1;
            return Err(format!(
                "trivial multiplicity by label {} disagrees with its CIR source {}",
                outcome.by_label, outcome.total
            ));
        }
        Ok(outcome)
    }

    fn check_record_partition(&mut self, ordinal: usize, tally: &RecordTally) {
        if let Some(error) = tally.probes.partition_error() {
            self.counts.accounting_violations += 1;
            self.anomaly(format!("ordinal {ordinal}: {error}"));
        }
        if let Some(error) = tally.entries.partition_error() {
            self.counts.accounting_violations += 1;
            self.anomaly(format!("ordinal {ordinal}: {error}"));
        }
        let geometry_partition = tally.geometry_reject_absent
            + tally.geometry_eligible_absent
            + tally.geometry_filter_errors;
        if geometry_partition != tally.absent_attempted {
            self.counts.accounting_violations += 1;
            self.anomaly(format!(
                "ordinal {ordinal}: geometry partition {geometry_partition} != attempted absent probes {}",
                tally.absent_attempted
            ));
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn frobenius(
        &mut self,
        ordinal: usize,
        sg: u8,
        child_sg: u8,
        subgroup: &IsotropySubgroup,
        probes: &[&'static IrrepRecord],
        stored_map: &BTreeMap<usize, StoredRow>,
        computed: &[Option<u32>],
        size: Option<u32>,
        subgroup_order: Option<usize>,
    ) {
        let gamma: Vec<usize> = probes
            .iter()
            .enumerate()
            .filter(|(_, probe)| is_gamma(probe))
            .map(|(index, _)| index)
            .collect();
        let direction = subgroup.record.direction_label;
        let arms = subgroup.record.arms;
        let size_text = size.map_or_else(|| "unset".to_string(), |value| value.to_string());
        self.counts.frobenius_records += 1;

        // The complex Frobenius convention only applies to Gamma condensates
        // with an unchanged translation lattice.
        if size != Some(1) {
            self.counts.frobenius_unevaluated += 1;
            self.emit(format!(
                "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_non_unit_size\tsize={size_text}"
            ));
            return;
        }
        if gamma.is_empty() {
            self.counts.frobenius_unevaluated += 1;
            self.emit(format!(
                "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_no_gamma_rows\t"
            ));
            return;
        }

        let mut category = "strict_complex";
        let mut numerator2: i128 = 0;
        let mut stored2: i128 = 0;
        for index in &gamma {
            let probe = probes[*index];
            let Some(value) = computed[*index] else {
                self.counts.frobenius_unevaluated += 1;
                self.emit(format!(
                    "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_probe_uncomputed\tprobe={}",
                    probe.ml
                ));
                return;
            };
            let Ok((component_sum, component_count)) = complex_component_dimensions(probe) else {
                self.counts.frobenius_unevaluated += 1;
                self.emit(format!(
                    "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_component_dimensions\tprobe={}",
                    probe.ml
                ));
                return;
            };
            // A Gamma probe's star is a single arm, so the row dimension is the
            // complex component sum of its source metadata.
            if component_sum != probe.dim {
                self.counts.frobenius_unevaluated += 1;
                self.emit(format!(
                    "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_gamma_dimension\tprobe={} row_dim={} component_sum={component_sum}",
                    probe.ml, probe.dim
                ));
                return;
            }
            if let IrrepSourceIdentity::Compound { .. } = probe.source_identity() {
                let distinct = matches!(
                    probe.compound_character_semantics(),
                    Some(CompoundCharacterSemantics::DistinctComponentSum)
                );
                category = if distinct {
                    "distinct_compound"
                } else if category == "strict_complex" {
                    "realification"
                } else {
                    category
                };
            }
            numerator2 +=
                2 * i128::from(value) * i128::from(component_sum) / i128::from(component_count);
            let stored_value = stored_map
                .get(index)
                .map_or(0, |row| i128::from(row.frequency));
            stored2 += 2 * stored_value * i128::from(component_sum) / i128::from(component_count);
        }
        let Some(parent_order) = self.parent_point_group_order(sg) else {
            self.counts.frobenius_unevaluated += 1;
            self.emit(format!(
                "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_parent_order\t"
            ));
            return;
        };
        let subgroup_order = subgroup_order.unwrap_or_default();
        if subgroup_order == 0 || parent_order % subgroup_order != 0 {
            self.counts.frobenius_unevaluated += 1;
            self.emit(format!(
                "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t\t\tunevaluated_index\tparent_order={parent_order} subgroup_order={subgroup_order}"
            ));
            return;
        }
        let index = parent_order / subgroup_order;
        let target2 = 2 * index as i128;
        let status = if numerator2 == target2 {
            self.counts.frobenius_pass += 1;
            "passed"
        } else {
            self.counts.frobenius_mismatch += 1;
            self.mismatch(format!(
                "ordinal {ordinal} parent {sg} subgroup {child_sg} direction {direction}: Frobenius \
                 sum {numerator2}/2 != index {index}"
            ));
            "mismatch"
        };
        match category {
            "strict_complex" => self.counts.frobenius_strict += 1,
            "realification" => self.counts.frobenius_realification += 1,
            _ => self.counts.frobenius_compound += 1,
        }
        self.emit(format!(
            "frobenius\t{ordinal}\t{sg}\t{child_sg}\t{direction}\t0\t{arms}\t{size_text}\t\t\t\t\t{index}\t{}\t{status}\tcategory={category} computed2={numerator2} stored2={stored2} parent_order={parent_order} subgroup_order={subgroup_order} gamma_probes={}",
            numerator2 / 2,
            gamma.len(),
        ));
    }

    fn finish(mut self) -> Result<(u8, Counts), String> {
        if let Some(error) = self.io_error.take() {
            return Err(format!("writing the report failed: {error}"));
        }
        self.writer
            .flush()
            .map_err(|error| format!("flushing the report failed: {error}"))?;
        let distinct_probes = self.distinct_probes.len();
        if self.check_census {
            self.census_check(distinct_probes);
        }
        if let Some(error) = self.counts.accounting_error() {
            self.counts.accounting_violations += 1;
            eprintln!("ACCOUNTING {error}");
        }
        let counts = self.counts.clone();
        let elapsed = self.started.elapsed().as_secs_f64();
        let scope = match (self.options.parent, self.options.ordinal) {
            (None, None) => "global".to_string(),
            (Some(parent), None) => format!("parent={parent}"),
            (None, Some(ordinal)) => format!("ordinal={ordinal}"),
            (Some(parent), Some(ordinal)) => format!("parent={parent} ordinal={ordinal}"),
        };
        eprintln!("=== task9 ordinary isotropy audit ===");
        eprintln!("scope={scope}");
        eprintln!(
            "condensate_records={} (pinned {EXPECTED_CONDENSATES}) embedding_ok={} embedding_failed={} no_trivial={}",
            counts.records,
            counts.records_embedding_ok,
            counts.records_embedding_failed,
            counts.records_no_trivial
        );
        eprintln!(
            "identity_rows={} unique_pairs={} (pinned {EXPECTED_IDENTITY_ROWS})",
            counts.entries.entries, counts.entries.unique
        );
        eprintln!(
            "other_wave_vector_records={} rows={} (pinned {EXPECTED_OTHER_WAVE_RECORDS}/{EXPECTED_OTHER_WAVE_ROWS})",
            counts.w_records, counts.w_entries
        );
        eprintln!(
            "spinor_records={} spinor_rows_in_identity={}",
            counts.spinor_records, counts.spinor_probe_rows
        );
        eprintln!("probes_total={}", counts.probes.total);
        eprintln!(
            "probe_partition: full_success={} identity_only={} missing={} error={} uncomputed_embedding={} uncomputed_no_trivial={} uncomputed={} partition_sum={} rows_emitted={}",
            counts.probes.full_success,
            counts.probes.identity_only,
            counts.probes.missing,
            counts.probes.error,
            counts.probes.uncomputed_embedding,
            counts.probes.uncomputed_no_trivial,
            counts.probes.uncomputed(),
            counts.probes.partition_sum(),
            counts.probe_rows_emitted
        );
        eprintln!(
            "engine_calls={} (full_success+missing+error)",
            counts.probes.attempted()
        );
        eprintln!(
            "positive_entries: entries={} unique={} passed={} mismatch={} embedding_unavailable={} no_trivial={} missing={} error={} unresolved={} duplicates_same={} duplicates_conflict={}",
            counts.entries.entries,
            counts.entries.unique,
            counts.entries.passed,
            counts.entries.mismatch,
            counts.entries.embedding_unavailable,
            counts.entries.no_trivial,
            counts.entries.missing,
            counts.entries.error,
            counts.entries.unresolved,
            counts.entries.duplicate_same,
            counts.entries.duplicate_conflict
        );
        eprintln!(
            "absent_zero={} (engine={} geometry={}) absent_positive={} geometry_reject_absent={} geometry_eligible_absent={} geometry_contradictions={} geometry_filter_errors={}",
            counts.absent_zero(),
            counts.absent_zero_engine,
            counts.absent_zero_geometry,
            counts.absent_positive,
            counts.geometry_reject_absent,
            counts.geometry_eligible_absent,
            counts.geometry_contradictions,
            counts.geometry_filter_errors
        );
        eprintln!(
            "frobenius: records={} passed={} mismatch={} unevaluated={} strict={} realification={} distinct_compound={}",
            counts.frobenius_records,
            counts.frobenius_pass,
            counts.frobenius_mismatch,
            counts.frobenius_unevaluated,
            counts.frobenius_strict,
            counts.frobenius_realification,
            counts.frobenius_compound
        );
        eprintln!(
            "production_checks: dimension_mismatch={} integrality_mismatch={} reconstruction_mismatch={} target_source_unmatched={} label_source_disagreement={}",
            counts.production_dim_mismatch,
            counts.production_integrality_mismatch,
            counts.production_recon_mismatch,
            counts.target_source_unmatched,
            counts.label_source_disagreement
        );
        eprintln!(
            "other_wave_vector: records={} rows={} emitted={} source_resolved={} source_mismatch={} computed={} engine_errors={} conflicts={}",
            counts.w_records,
            counts.w_entries,
            counts.w_rows_emitted,
            counts.w_source_resolved,
            counts.w_source_mismatch,
            counts.w_computed,
            counts.w_engine_error,
            counts.w_conflict
        );
        eprintln!(
            "w_scope: rows={} emitted={} computed={} uncomputed={} mismatched={} engine_errors={} character_tables_frozen={} character_tables_blocked={} reason=none gate=--require-w-complete",
            counts.w_entries,
            counts.w_rows_emitted,
            counts.w_computed,
            counts.w_incomplete(),
            counts.w_frequency_mismatch,
            counts.w_engine_error,
            counts.w_character_frozen,
            counts.w_character_blocked
        );
        eprintln!(
            "w_parameter_shift: checked={} mismatched={} comparison_skipped={} same_label_witnesses={} (t = 5/4 against the monodromy image at t = 1/4, complete decomposition per block and target)",
            counts.w_parameter_shift_checked,
            counts.w_parameter_shift_mismatch,
            counts.w_parameter_shift_skipped,
            counts.w_parameter_shift_witnesses
        );
        eprintln!(
            "completeness: missing_probes={} uncomputed_probes={} uncomputed_entries={} unresolved_entries={} duplicate_same={} geometry_filter_errors={} frobenius_unevaluated={} basis_errors={} w_uncomputed={}",
            counts.probes.missing,
            counts.probes.uncomputed(),
            counts.entries.incomplete(),
            counts.entries.unresolved,
            counts.entries.duplicate_same,
            counts.geometry_filter_errors,
            counts.frobenius_unevaluated,
            counts.basis_errors,
            counts.w_entries.saturating_sub(counts.w_computed)
        );
        eprintln!(
            "coverage: embedded_records={}/{} positive_stored_compared={}/{} probe_full_success={}/{} probe_identity_only={}/{} probe_answered={}/{} absent_zero={} geometry_reject_absent={} frobenius_evaluated={}/{} w_computed={}/{}",
            counts.records_embedding_ok,
            counts.records,
            counts.positive_comparisons(),
            counts.entries.unique,
            counts.probes.full_success,
            counts.probes.total,
            counts.probes.identity_only,
            counts.probes.total,
            counts.probes.full_success + counts.probes.identity_only,
            counts.probes.total,
            counts.absent_zero(),
            counts.geometry_reject_absent,
            counts.frobenius_pass,
            counts.frobenius_records,
            counts.w_computed,
            counts.w_entries
        );
        eprintln!(
            "distinct_referenced_probes={distinct_probes} stored_k_mismatch={} stored_domain_out_of_range={} gamma_record_non_gamma_probe={} hard_failures={} accounting_violations={} census_mismatch={}",
            counts.stored_k_mismatch,
            counts.stored_domain_out_of_range,
            counts.gamma_record_non_gamma_probe,
            counts.hard_failures(),
            counts.accounting_violations,
            counts.census_mismatch
        );
        let mut detail: Vec<String> = self
            .error_detail
            .iter()
            .map(|(key, value)| format!("{key}={value}"))
            .collect();
        detail.sort();
        if !detail.is_empty() {
            eprintln!("error_detail: {}", detail.join(" "));
        }
        for message in &self.anomalies {
            eprintln!("ANOMALY {message}");
        }
        for message in &self.mismatches {
            eprintln!("MISMATCH {message}");
        }
        eprintln!("elapsed={elapsed:.1}s");
        // The full-decomposition state is reported on its own line, for the
        // scope this run actually covers.  A scoped run can be locally complete
        // and still say nothing about the global denominator, so that
        // distinction is printed rather than implied.
        let scoped = self.options.scoped();
        eprintln!(
            "full_decomposition: scope={scope} probes={} full_success={} identity_only={} \
             missing={} error={} uncomputed={} incomplete={} global={} gate=--require-full-decomposition",
            counts.probes.total,
            counts.probes.full_success,
            counts.probes.identity_only,
            counts.probes.missing,
            counts.probes.error,
            counts.probes.uncomputed(),
            counts.full_decomposition_incomplete(),
            if scoped { "not_established" } else { "covered" }
        );
        let gates = Gates::from_options(&self.options);
        let verdict = counts.verdict(gates);
        let incomplete = counts.incomplete();
        let hard = counts.hard_failures();
        match verdict {
            Verdict::Inconsistent => {
                eprintln!("VERDICT inconsistent scope={scope} hard_failures={hard}");
            }
            Verdict::Incomplete => {
                eprintln!(
                    "VERDICT incomplete scope={scope} gates={} incomplete_categories={incomplete} \
                     full_decomposition_incomplete={} w_uncomputed={}",
                    gates.label(),
                    counts.full_decomposition_incomplete(),
                    counts.w_incomplete()
                );
            }
            Verdict::Clean if gates.any() => {
                if scoped {
                    eprintln!(
                        "VERDICT complete scope={scope} gates={} global_coverage=not_established \
                         (a scoped run cannot certify the global table)",
                        gates.label()
                    );
                } else {
                    eprintln!(
                        "VERDICT complete scope={scope} gates={} full_decomposition={}",
                        gates.label(),
                        if gates.full_decomposition {
                            "complete"
                        } else {
                            "not_gated"
                        }
                    );
                }
            }
            Verdict::Clean => {
                eprintln!(
                    "VERDICT clean scope={scope} gates=none incomplete_categories={incomplete} \
                     full_decomposition_incomplete={} w_uncomputed={}",
                    counts.full_decomposition_incomplete(),
                    counts.w_incomplete()
                );
            }
        }
        Ok((counts.exit_code(gates), counts))
    }

    fn census_check(&mut self, distinct_probes: usize) {
        let counts = &mut self.counts;
        let mut mismatch = false;
        let report = |name: &str, found: usize, pinned: usize, mismatch: &mut bool| {
            if found != pinned {
                eprintln!("CENSUS {name} {found} != pinned {pinned}");
                *mismatch = true;
            }
        };
        report(
            "condensate records",
            counts.records,
            EXPECTED_CONDENSATES,
            &mut mismatch,
        );
        report(
            "identity rows",
            counts.entries.entries,
            EXPECTED_IDENTITY_ROWS,
            &mut mismatch,
        );
        report(
            "unique identity rows",
            counts.entries.unique,
            EXPECTED_IDENTITY_ROWS,
            &mut mismatch,
        );
        report(
            "other-wave rows",
            counts.w_entries,
            EXPECTED_OTHER_WAVE_ROWS,
            &mut mismatch,
        );
        report(
            "other-wave records",
            counts.w_records,
            EXPECTED_OTHER_WAVE_RECORDS,
            &mut mismatch,
        );
        report(
            "spinor records",
            counts.spinor_records,
            EXPECTED_SPINOR_RECORDS,
            &mut mismatch,
        );
        report(
            "distinct referenced probes",
            distinct_probes,
            EXPECTED_DISTINCT_PROBES,
            &mut mismatch,
        );
        if mismatch {
            counts.census_mismatch = 1;
        }
    }
}

// ── Pure helpers ─────────────────────────────────────────────────────────────

fn insert_stored_row(
    map: &mut BTreeMap<usize, StoredRow>,
    index: usize,
    row: StoredRow,
    entries: &mut EntryCounts,
) -> Result<(), String> {
    match map.get(&index) {
        None => {
            map.insert(index, row);
            Ok(())
        }
        Some(existing) if existing.frequency == row.frequency => {
            entries.duplicate_same += 1;
            Ok(())
        }
        Some(existing) => {
            entries.duplicate_conflict += 1;
            Err(format!(
                "probe {} has conflicting stored frequencies {} and {} (domains {} and {}); \
                 repeated rows must not be summed",
                row.ml, existing.frequency, row.frequency, existing.domain, row.domain
            ))
        }
    }
}

fn find_trivial_child(sg: u8) -> Option<&'static IrrepRecord> {
    let mut found: Option<&'static IrrepRecord> = None;
    for record in query::irreps_of(sg) {
        if record.spinor || !is_gamma(record) || record.dim != 1 {
            continue;
        }
        let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() else {
            continue;
        };
        if row.dimension() == 1 && row.values().iter().all(|value| (value - 1.0).norm() < 1e-9) {
            if found.is_some() {
                return None;
            }
            found = Some(record);
        }
    }
    found
}

fn is_gamma(record: &IrrepRecord) -> bool {
    record.kx == 0 && record.ky == 0 && record.kz == 0
}

fn build_child_reciprocal(sg: u8) -> Result<Lattice, String> {
    let basis = parent_primitive_basis(sg).map_err(|error| error.to_string())?;
    let mut rows = [[Rat::ZERO; 3]; 3];
    for (row, values) in basis.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            rows[row][column] = Rat::from_grid(*value, 12).map_err(|error| error.to_string())?;
        }
    }
    let cell = Lattice::new(Mat3R::new(rows)).map_err(|error| error.to_string())?;
    cell.reciprocal().map_err(|error| error.to_string())
}

fn probe_is_gamma_eligible(
    star: &ScalarStar,
    embedding: &SubgroupEmbedding,
    child_reciprocal: &Lattice,
) -> Result<bool, String> {
    let folded = star
        .folded_stars(embedding)
        .map_err(|error| error.to_string())?;
    for child_star in &folded {
        for point in child_star.points() {
            if child_reciprocal
                .contains(point.q())
                .map_err(|error| error.to_string())?
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn classify_call_error(error: &FullStarError) -> CallErrorKind {
    match error {
        FullStarError::MissingChildStarData { .. } => CallErrorKind::MissingChildData,
        FullStarError::Subduction(SubductionError::MissingIrrepData { .. }) => {
            CallErrorKind::MissingIrrepData
        }
        _ => CallErrorKind::Hard,
    }
}

fn embedding_label(error: &SubductionError) -> &'static str {
    match error {
        SubductionError::NoValidEmbedding { .. } => "no_valid_embedding",
        SubductionError::AmbiguousEmbedding { .. } => "ambiguous_embedding",
        SubductionError::StrictHallUnavailable { .. } => "strict_hall_unavailable",
        SubductionError::HallLoadFailed { .. } => "hall_load_failed",
        SubductionError::StaleIsotropyRecord { .. } => "stale_isotropy_record",
        SubductionError::IsotropyOrdinalOutOfRange { .. } => "ordinal_out_of_range",
        SubductionError::InvalidSubgroupNumber { .. } => "invalid_subgroup_number",
        SubductionError::FrozenEmbeddingRejected { .. } => "frozen_embedding_rejected",
        _ => "other",
    }
}

/// True complex component dimensions and their count for one parent probe.
///
/// Ordinary rows have one component; compound rows carry exactly two CIR
/// constituents in their frozen metadata.
fn complex_component_dimensions(probe: &IrrepRecord) -> Result<(u8, u8), String> {
    match probe.source_identity() {
        IrrepSourceIdentity::OrdinaryScalar { .. } => Ok((probe.dim, 1)),
        IrrepSourceIdentity::Compound { .. } => {
            let metadata = probe
                .compound_metadata()
                .ok_or_else(|| format!("compound record {} has no metadata", probe.ml))?;
            let [first, second] = metadata.cir_dimensions;
            if first != second {
                return Err(format!(
                    "compound record {} has unequal component dimensions {first} and {second}",
                    probe.ml
                ));
            }
            Ok((first + second, 2))
        }
        IrrepSourceIdentity::Spin { .. } => Err(format!("spinor record {}", probe.ml)),
    }
}

fn normalize_k(k: KVector) -> Option<([i32; 3], i32)> {
    if k.denominator == 0 {
        return None;
    }
    let mut denominator = i32::from(k.denominator);
    let mut numerators = [
        i32::from(k.numerators[0]),
        i32::from(k.numerators[1]),
        i32::from(k.numerators[2]),
    ];
    if denominator < 0 {
        denominator = -denominator;
        numerators = [-numerators[0], -numerators[1], -numerators[2]];
    }
    let mut divisor = denominator;
    for value in numerators {
        divisor = gcd_i32(divisor, value.abs());
    }
    if divisor == 0 {
        divisor = 1;
    }
    Some((
        [
            numerators[0] / divisor,
            numerators[1] / divisor,
            numerators[2] / divisor,
        ],
        denominator / divisor,
    ))
}

fn same_k(left: KVector, right: KVector) -> bool {
    match (normalize_k(left), normalize_k(right)) {
        (Some(left), Some(right)) => left == right,
        _ => false,
    }
}

fn gcd_i32(mut left: i32, mut right: i32) -> i32 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.abs()
}

fn format_k(k: KVector) -> String {
    match normalize_k(k) {
        Some((numerators, denominator)) => format!(
            "{}/{},{}/{},{}/{}",
            numerators[0], denominator, numerators[1], denominator, numerators[2], denominator
        ),
        None => "invalid/0".to_string(),
    }
}

fn probe_source_label(probe: &IrrepRecord) -> &'static str {
    match probe.source_identity() {
        IrrepSourceIdentity::OrdinaryScalar { .. } => "ordinary",
        IrrepSourceIdentity::Compound { .. } => "compound",
        IrrepSourceIdentity::Spin { .. } => "spin",
    }
}

/// A parameter-independent summary of one line decomposition: every block with
/// its stored child `k`, star size, arm count, dimensions and every target with
/// its multiplicity, label, CIR number and component identity.
///
/// Blocks are compared as a **multiset**: a parameter shift may reorder the same
/// blocks (measured: ordinal 10030 `SM1` swaps its two blocks between `t = 1/4`
/// and `t = 3/4` with identical contents), and the order is not part of the
/// decomposition.  Targeting `q` is deliberately excluded: it is the unreduced
/// folded coordinate, so it may legitimately differ by a child reciprocal
/// vector between equivalent parameters.
fn line_decomposition_key(
    result: &cryspglib::irrep::subduction::star::decompose::LineSubduction,
) -> Vec<String> {
    let mut keys: Vec<String> = result
        .blocks()
        .iter()
        .map(|block| {
            let targets: Vec<String> = block
                .targets()
                .iter()
                .map(|target| {
                    format!(
                        "{}x{} {:?} {:?} {:?}",
                        target.multiplicity,
                        target.dimension,
                        target.ml,
                        target.irnumber,
                        target.component
                    )
                })
                .collect();
            format!(
                "k=({},{},{}) star={} arms={} dim={} little={} [{}]",
                block.stored_k().get(0),
                block.stored_k().get(1),
                block.stored_k().get(2),
                block.star_size(),
                block.arm_count(),
                block.block_dimension(),
                block.little_dimension(),
                targets.join("; ")
            )
        })
        .collect();
    keys.sort();
    keys
}

fn production_detail(outcome: &CallOutcome, geometry: Geometry) -> String {
    format!(
        "geometry={} blocks_targets={} targets_without_source={} trivial_by_label={}",
        geometry.label(),
        outcome.targets,
        outcome.targets_without_source,
        outcome.by_label
    )
}

/// Common coordinates of one isotropy record's probe rows.
struct RecordContext<'a> {
    ordinal: usize,
    sg: u8,
    child_sg: u8,
    subgroup: &'a IsotropySubgroup,
    size: Option<u32>,
}

fn identity_line(
    context: &RecordContext<'_>,
    kind: &str,
    probe: &IrrepRecord,
    stored: Option<StoredRow>,
    computed: Option<&CallOutcome>,
    status: &str,
    detail: &str,
) -> String {
    let RecordContext {
        ordinal,
        sg,
        child_sg,
        subgroup,
        size,
    } = *context;
    format!(
        "{kind}\t{ordinal}\t{sg}\t{child_sg}\t{}\t{}\t{}\t{}\t{}\t{sg}\t{}\t{}\t{}\t{}\t{status}\t{detail}",
        subgroup.record.direction_label,
        stored.map_or(0, |row| row.domain),
        subgroup.record.arms,
        size.map_or_else(|| "unset".to_string(), |value| value.to_string()),
        probe.ml,
        format_k(probe.k_vector()),
        probe_source_label(probe),
        stored.map_or_else(String::new, |row| row.frequency.to_string()),
        computed.map_or_else(String::new, |outcome| outcome.total.to_string()),
    )
}

// ── Entry point ──────────────────────────────────────────────────────────────

fn audit_space_group(
    options: &Options,
    sg: u8,
    progress: Option<Arc<ProgressReporter>>,
) -> Result<AuditPart, String> {
    let report = ReportBuffer::default();
    let mut auditor = Auditor::new(options.clone(), Box::new(report.clone()));
    auditor.progress = progress;
    auditor.audit_sg(sg)?;
    auditor.finish_part(report.contents())
}

fn collect_audit_parts(
    options: &Options,
    space_groups: &[u8],
    parallel: bool,
    progress: Option<Arc<ProgressReporter>>,
) -> Vec<Result<AuditPart, String>> {
    if parallel {
        space_groups
            .par_iter()
            .map(|&sg| audit_space_group(options, sg, progress.clone()))
            .collect()
    } else {
        space_groups
            .iter()
            .map(|&sg| audit_space_group(options, sg, progress.clone()))
            .collect()
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => code,
        Err(message) => {
            eprintln!("audit_irrep_subduction: {message}");
            ExitCode::from(3)
        }
    }
}

fn run() -> Result<ExitCode, String> {
    let options = Options::parse(env::args().skip(1))?;
    if options.help {
        println!("{USAGE}");
        return Ok(ExitCode::SUCCESS);
    }
    let writer: Box<dyn Write> = match &options.output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).map_err(|error| format!("cannot create {path}: {error}"))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let (exit_code, _counts) = run_audit(options, writer)?;
    Ok(ExitCode::from(exit_code))
}

fn run_audit(options: Options, writer: Box<dyn Write>) -> Result<(u8, Counts), String> {
    let space_groups = match options.parent {
        Some(parent) => vec![parent],
        None => (1..=230u8).collect(),
    };
    let check_census = !options.scoped();
    run_audit_with_space_groups(options, writer, &space_groups, true, check_census)
}

fn run_audit_with_space_groups(
    options: Options,
    writer: Box<dyn Write>,
    space_groups: &[u8],
    parallel: bool,
    check_census: bool,
) -> Result<(u8, Counts), String> {
    let started = Instant::now();
    let progress =
        (options.progress > 0).then(|| Arc::new(ProgressReporter::new(options.progress, started)));
    let mut auditor = Auditor::new(options, writer);
    auditor.started = started;
    auditor.check_census = check_census;
    auditor.emit(HEADER.to_string());
    // A slice is an indexed Rayon iterator, so collect preserves the input
    // order while the independent space-group jobs run concurrently.
    let parts = collect_audit_parts(&auditor.options, space_groups, parallel, progress);
    for part in parts {
        auditor.absorb_part(part?);
    }
    // An empty scope must never read as a clean run: `--parent 2 --ordinal 0`
    // selects nothing (ordinal 0 belongs to SG 1), and a CI loop over
    // (parent, ordinal) pairs would otherwise get `VERDICT clean` with zero
    // probes for a typo.  The pinned ordinal space is checked by ownership
    // here, not only by bounds.
    if auditor.options.scoped() && auditor.counts.records == 0 {
        let scope = match (auditor.options.parent, auditor.options.ordinal) {
            (Some(parent), Some(ordinal)) => format!("--parent {parent} --ordinal {ordinal}"),
            (Some(parent), None) => format!("--parent {parent}"),
            (None, Some(ordinal)) => format!("--ordinal {ordinal}"),
            (None, None) => unreachable!("an unscoped run is never empty"),
        };
        return Err(format!(
            "empty scope: {scope} selects no non-spinor isotropy record; the ordinal does \
             not belong to that parent (or the parent has no condensates), so a report \
             over zero probes would prove nothing"
        ));
    }
    auditor.finish()
}

// ── Permanent tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    #[derive(Clone, Default)]
    struct SharedSink(Arc<Mutex<Vec<u8>>>);

    impl SharedSink {
        fn text(&self) -> String {
            String::from_utf8(self.0.lock().unwrap().clone()).expect("report is UTF-8")
        }
    }

    impl Write for SharedSink {
        fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buffer);
            Ok(buffer.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn audit_ordinal_gates(ordinal: usize, gates: Gates) -> (u8, Counts, String) {
        let sink = SharedSink::default();
        let options = Options {
            ordinal: Some(ordinal),
            require_complete: gates.complete,
            require_full_decomposition: gates.full_decomposition,
            require_w_complete: gates.w_complete,
            progress: 0,
            ..Options::default()
        };
        let (exit_code, counts) =
            run_audit(options, Box::new(sink.clone())).expect("scoped audit runs");
        (exit_code, counts, sink.text())
    }

    fn audit_ordinal(ordinal: usize, require_complete: bool) -> (u8, Counts, String) {
        audit_ordinal_gates(
            ordinal,
            Gates {
                complete: require_complete,
                ..Gates::default()
            },
        )
    }

    #[test]
    fn parallel_space_groups_match_serial_run_output_and_merged_counts() {
        let options = Options {
            progress: 0,
            ..Options::default()
        };
        // Deliberately use a non-sorted order: collection must preserve the
        // caller's order even when groups finish at different times.
        let space_groups = [24, 1, 5];
        let parallel_sink = SharedSink::default();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .expect("parallel test pool");
        let (parallel_exit, parallel_counts) = pool
            .install(|| {
                run_audit_with_space_groups(
                    options.clone(),
                    Box::new(parallel_sink.clone()),
                    &space_groups,
                    true,
                    false,
                )
            })
            .expect("parallel space-group audit succeeds");
        assert!(
            parallel_sink
                .text()
                .lines()
                .any(|line| line.starts_with("record\t")),
            "parallel report must contain emitted rows"
        );

        let serial_sink = SharedSink::default();
        let (serial_exit, serial_counts) = run_audit_with_space_groups(
            options,
            Box::new(serial_sink.clone()),
            &space_groups,
            false,
            false,
        )
        .expect("serial space-group audit succeeds");

        assert_eq!(parallel_sink.text(), serial_sink.text());
        assert_eq!(format!("{parallel_counts:?}"), format!("{serial_counts:?}"));
        assert_eq!(parallel_exit, serial_exit);
    }

    #[test]
    fn child_source_cache_initializes_safely_under_concurrent_first_access() {
        let caches = SharedChildCaches::new();
        let barrier = Arc::new(std::sync::Barrier::new(3));
        let caches = &caches;
        let addresses = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..3)
                .map(|_| {
                    let barrier = Arc::clone(&barrier);
                    let caches = caches;
                    scope.spawn(move || {
                        barrier.wait();
                        caches.child_sources(1) as *const ChildSources as usize
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|handle| handle.join().expect("cache worker panicked"))
                .collect::<Vec<_>>()
        });

        assert!(addresses.iter().all(|address| *address == addresses[0]));
    }

    /// A scope that selects nothing is an error, not a clean run of zero probes.
    ///
    /// Ordinal 0 belongs to SG 1, so `--parent 2 --ordinal 0` is a contradiction;
    /// before this guard the audit reported `VERDICT clean` with
    /// `probe_full_success=0/0` and exit 0 (reproduced before the fix).
    #[test]
    fn an_empty_scope_is_rejected_instead_of_reporting_clean() {
        let sink = SharedSink::default();
        let options = Options {
            parent: Some(2),
            ordinal: Some(0),
            ..Options::default()
        };
        let error =
            run_audit(options, Box::new(sink.clone())).expect_err("an empty scope must fail");
        assert!(error.contains("empty scope"), "{error}");
        // The same ordinal under its real parent still runs, and so does a
        // parent-only scope, so the guard rejects only the empty case.
        let options = Options {
            parent: Some(1),
            ordinal: Some(0),
            ..Options::default()
        };
        assert!(run_audit(options, Box::new(SharedSink::default())).is_ok());
        let options = Options {
            parent: Some(1),
            ..Options::default()
        };
        assert!(run_audit(options, Box::new(SharedSink::default())).is_ok());
    }

    /// Every combination of the three independent gates.
    fn all_gate_combinations() -> Vec<Gates> {
        let mut combinations = Vec::new();
        for complete in [false, true] {
            for full_decomposition in [false, true] {
                for w_complete in [false, true] {
                    combinations.push(Gates {
                        complete,
                        full_decomposition,
                        w_complete,
                    });
                }
            }
        }
        combinations
    }

    fn gates_with_full_decomposition() -> Gates {
        Gates {
            full_decomposition: true,
            ..Gates::default()
        }
    }

    fn all_gates() -> Gates {
        Gates {
            complete: true,
            full_decomposition: true,
            w_complete: true,
        }
    }

    fn sg16_r2_p1() -> IsotropySubgroup {
        isotropy::isotropy_subgroup_for_direction(
            16,
            "R2",
            LabelConvention::Cdml,
            isotropy::IsotropyDirection::Label("P1"),
        )
        .expect("SG 16 R2 P1 resolves")
    }

    #[derive(Debug)]
    struct EmittedProbeRow {
        probe: String,
        computed: String,
        status: String,
        detail: String,
    }

    /// One row per scalar probe, whatever the outcome; `record`, `frobenius` and
    /// `w_entry` rows are not probe rows.
    fn emitted_probe_rows(tsv: &str) -> Vec<EmittedProbeRow> {
        tsv.lines()
            .filter_map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() != 16 || !matches!(fields[0], "identity" | "absent") {
                    return None;
                }
                Some(EmittedProbeRow {
                    probe: fields[8].to_string(),
                    computed: fields[13].to_string(),
                    status: fields[14].to_string(),
                    detail: fields[15].to_string(),
                })
            })
            .collect()
    }

    // ── Verdict plumbing ─────────────────────────────────────────────────────

    #[test]
    fn non_gamma_stored_probe_for_gamma_condensate_fails_the_verdict() {
        let counts = Counts {
            gamma_record_non_gamma_probe: 1,
            ..Counts::default()
        };
        assert_eq!(counts.exit_code(Gates::default()), 1);
        assert_eq!(
            counts.exit_code(Gates {
                complete: true,
                full_decomposition: true,
                w_complete: true,
            }),
            1
        );
    }

    /// The three gates are separate switches: the new flag must parse, default
    /// off, and stay visible in the usage text.
    #[test]
    fn the_full_decomposition_flag_parses_and_defaults_off() {
        assert!(!Options::default().require_full_decomposition);
        let options = Options::parse(["--require-full-decomposition".to_string()].into_iter())
            .expect("the flag parses");
        assert!(options.require_full_decomposition);
        assert!(!options.require_complete && !options.require_w_complete);
        assert!(USAGE.contains("--require-full-decomposition"));
        assert_eq!(
            Gates::from_options(&options).label(),
            "--require-full-decomposition"
        );
    }

    /// A computed frequency that disagrees with the pinned one is a *known*
    /// error, so it must fail the run under every combination of the strict
    /// switches -- not only under the one that also counts uncomputed rows.
    #[test]
    fn a_frequency_mismatch_fails_under_every_flag_combination() {
        let counts = Counts {
            w_frequency_mismatch: 1,
            ..Counts::default()
        };
        assert!(counts.hard_failures() >= 1, "a mismatch is a hard failure");
        for gates in all_gate_combinations() {
            assert_eq!(counts.exit_code(gates), 1, "gates={}", gates.label());
        }
    }

    /// An `Err` from the block route is an arithmetic/context failure of the
    /// computation itself.  It must be a hard failure under every flag
    /// combination, exactly like a computed value that disagrees with the pinned
    /// one -- a demonstrated error is neither a clean run nor missing data.
    #[test]
    fn a_w_engine_error_fails_under_every_flag_combination() {
        let counts = Counts {
            w_engine_error: 1,
            w_entries: 1,
            ..Counts::default()
        };
        assert!(
            counts.hard_failures() >= 1,
            "an engine error is a hard failure"
        );
        // The failed row is also still "uncomputed", so the w gate sees it; the
        // hard failure must dominate the incompleteness verdict.
        assert_eq!(counts.w_incomplete(), 1);
        for gates in all_gate_combinations() {
            assert_eq!(counts.exit_code(gates), 1, "gates={}", gates.label());
        }
    }

    /// A transport mismatch is a hard failure under every gate combination.
    /// This is the counter that makes the R6.2 monodromy class (a parameter step
    /// that does not land on the image label, or a wave vector canonicalized back
    /// onto the original parameter) impossible to ship again: the identity
    /// frequency alone does not see it, and neither does a gate that only asks
    /// whether *some* decomposition came out.
    #[test]
    fn a_line_parameter_shift_mismatch_fails_under_every_flag_combination() {
        let counts = Counts {
            w_parameter_shift_mismatch: 1,
            w_entries: 1,
            w_computed: 1,
            w_parameter_shift_checked: 1,
            ..Counts::default()
        };
        assert!(
            counts.hard_failures() >= 1,
            "a parameter step that does not transport is a hard failure"
        );
        // Every row is computed, so no completeness gate sees it: only the hard
        // failure can stop this run.
        assert_eq!(counts.w_incomplete(), 0);
        for gates in all_gate_combinations() {
            assert_eq!(counts.exit_code(gates), 1, "gates={}", gates.label());
        }
        // An undetermined image is a hard failure: the audit cannot certify
        // transport without knowing the monodromy image, and must not silently
        // accept a skipped comparison.
        let skipped = Counts {
            w_parameter_shift_skipped: 3,
            w_entries: 1,
            w_rows_emitted: 1,
            w_computed: 1,
            ..Counts::default()
        };
        assert_eq!(skipped.hard_failures(), 3);
        assert_eq!(skipped.w_incomplete(), 0);
        for gates in all_gate_combinations() {
            assert_eq!(skipped.exit_code(gates), 1, "gates={}", gates.label());
        }
    }

    #[test]
    fn line_image_diagnostics_distinguish_unsupported_ambiguous_and_missing() {
        use cryspglib::irrep::line_monodromy::LabelImage;

        assert_eq!(
            line_image_skip_reason(Some(&LabelImage::UnsupportedShift)),
            (
                "computed_transport_skipped:unsupported_shift",
                "exact little-group rotation criterion for the phase twist is not met"
            )
        );
        assert_eq!(
            line_image_skip_reason(Some(&LabelImage::Ambiguous(vec!["DT1", "DT2"]))),
            (
                "computed_transport_skipped:ambiguous_image",
                "frozen fingerprints give multiple image labels"
            )
        );
        assert_eq!(
            line_image_skip_reason(Some(&LabelImage::Missing)),
            (
                "computed_transport_skipped:missing_image",
                "frozen data has no image label"
            )
        );
    }

    #[test]
    fn a_skipped_line_transport_names_the_affected_row() {
        let mut auditor = Auditor::new(Options::default(), Box::new(SharedSink::default()));
        auditor.note_w_transport_skipped(10030, "DT1", "no unique monodromy image");

        assert_eq!(auditor.counts.w_parameter_shift_skipped, 1);
        assert_eq!(auditor.counts.hard_failures(), 1);
        assert!(auditor.mismatches[0].contains("ordinal 10030"));
        assert!(auditor.mismatches[0].contains("DT1"));
        assert!(auditor.mismatches[0].contains("no unique monodromy image"));
    }

    #[test]
    fn missing_w_report_rows_are_an_accounting_error() {
        let counts = Counts {
            w_entries: 1,
            ..Counts::default()
        };
        assert_eq!(
            counts.accounting_error().as_deref(),
            Some("emitted other-wave-vector rows 0 != w_entries 1")
        );
    }

    /// An identity-only result closes the *identity* gate and is a gap for the
    /// *full-decomposition* gate.  The two questions must stay separate: this is
    /// the distinction the old two gates could not express.
    #[test]
    fn a_full_decomposition_gap_fails_only_the_full_gate() {
        let counts = Counts {
            probes: ProbeCounts {
                total: 3,
                full_success: 2,
                identity_only: 1,
                ..ProbeCounts::default()
            },
            ..Counts::default()
        };
        assert_eq!(counts.hard_failures(), 0);
        assert_eq!(
            counts.incomplete(),
            0,
            "identity-only closes the identity gate"
        );
        assert_eq!(counts.full_decomposition_incomplete(), 1);
        assert_eq!(
            counts.verdict(Gates {
                complete: true,
                ..Gates::default()
            }),
            Verdict::Clean
        );
        assert_eq!(counts.exit_code(Gates::default()), 0);
        // The new gate reports the gap as incompleteness, never as an error...
        assert_eq!(
            counts.verdict(gates_with_full_decomposition()),
            Verdict::Incomplete
        );
        assert_eq!(counts.exit_code(gates_with_full_decomposition()), 2);
        // ...and it is the only switch that reacts to it.
        for gates in all_gate_combinations() {
            let expected = if gates.full_decomposition { 2 } else { 0 };
            assert_eq!(counts.exit_code(gates), expected, "gates={}", gates.label());
        }
    }

    /// A demonstrated error outranks every coverage gate, including the new one.
    #[test]
    fn a_hard_failure_outranks_the_full_decomposition_gate() {
        let counts = Counts {
            probes: ProbeCounts {
                total: 2,
                identity_only: 1,
                error: 1,
                ..ProbeCounts::default()
            },
            ..Counts::default()
        };
        assert!(counts.hard_failures() >= 1);
        assert_eq!(counts.full_decomposition_incomplete(), 1);
        for gates in all_gate_combinations() {
            assert_eq!(counts.exit_code(gates), 1, "gates={}", gates.label());
        }
    }

    /// The wiring, not the production increments: every counted production
    /// violation is summed into `hard_failures` and turns the exit code into 1
    /// under **every** gate combination, so a run cannot trade one of those
    /// invariants for coverage.
    ///
    /// Honest scope (reviewer B): this test sets the counters directly and
    /// therefore does not prove that the production paths increment them -- it
    /// pins `Counts::hard_failures` and `Counts::exit_code` only.  Two of the
    /// five are unreachable on this corpus anyway and the other three compare
    /// objects the engine has already validated itself; the report's
    /// independence table says which evidence is real.
    #[test]
    fn every_counted_production_violation_is_a_hard_failure() {
        for (name, counts) in [
            (
                "dimension",
                Counts {
                    production_dim_mismatch: 1,
                    ..Counts::default()
                },
            ),
            (
                "integrality",
                Counts {
                    production_integrality_mismatch: 1,
                    ..Counts::default()
                },
            ),
            (
                "reconstruction",
                Counts {
                    production_recon_mismatch: 1,
                    ..Counts::default()
                },
            ),
            (
                "source identity",
                Counts {
                    target_source_unmatched: 1,
                    ..Counts::default()
                },
            ),
            (
                "label agreement",
                Counts {
                    label_source_disagreement: 1,
                    ..Counts::default()
                },
            ),
        ] {
            assert!(counts.hard_failures() >= 1, "{name} must be a hard failure");
            assert_eq!(counts.exit_code(Gates::default()), 1, "{name}: exit code");
            for gates in all_gate_combinations() {
                assert_eq!(
                    counts.exit_code(gates),
                    1,
                    "{name}: gates={}",
                    gates.label()
                );
            }
        }
    }

    #[test]
    fn a_run_without_mismatches_still_exits_zero() {
        let counts = Counts::default();
        assert_eq!(counts.hard_failures(), 0);
        assert_eq!(counts.w_engine_error, 0);
        assert_eq!(counts.full_decomposition_incomplete(), 0);
        for gates in all_gate_combinations() {
            assert_eq!(counts.exit_code(gates), 0, "gates={}", gates.label());
        }
    }

    // ── Pinned census and the per-record acceptance witnesses ────────────────

    #[test]
    fn pinned_census_and_identity_rows_are_unique_per_probe() {
        let mut condensates = 0usize;
        let mut identity_rows = 0usize;
        let mut wave_rows = 0usize;
        let mut wave_records = 0usize;
        let mut spinor_records = 0usize;
        let mut spinor_rows = 0usize;
        let mut distinct_probes: HashSet<(u8, u16)> = HashSet::new();
        let mut wave_labels: HashSet<(u8, &'static str)> = HashSet::new();
        for sg in 1..=230u8 {
            let records = query::irreps_of(sg);
            spinor_records += records.iter().filter(|record| record.spinor).count();
            for (index, record) in records.iter().enumerate() {
                if record.spinor {
                    continue;
                }
                for subgroup in
                    isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml).unwrap()
                {
                    condensates += 1;
                    let stored = subgroup.identity_subduction().unwrap();
                    identity_rows += stored.len();
                    let mut seen: HashSet<usize> = HashSet::new();
                    for entry in &stored {
                        let probe = records
                            .iter()
                            .position(|candidate| candidate.ml == entry.parent_ml)
                            .expect("stored probe resolves");
                        if records[probe].spinor {
                            spinor_rows += 1;
                        }
                        assert!(
                            seen.insert(probe),
                            "ordinal {} repeats probe {}",
                            subgroup.ordinal,
                            entry.parent_ml
                        );
                    }
                    let wave = subgroup.other_wave_vector_subduction().unwrap();
                    if !wave.is_empty() {
                        wave_records += 1;
                    }
                    for entry in &wave {
                        wave_rows += 1;
                        wave_labels.insert((sg, entry.parent_ml));
                    }
                }
                distinct_probes.insert((sg, index as u16));
            }
        }
        assert_eq!(condensates, EXPECTED_CONDENSATES);
        assert_eq!(identity_rows, EXPECTED_IDENTITY_ROWS);
        assert_eq!(wave_rows, EXPECTED_OTHER_WAVE_ROWS);
        assert_eq!(wave_records, EXPECTED_OTHER_WAVE_RECORDS);
        assert_eq!(spinor_records, EXPECTED_SPINOR_RECORDS);
        assert_eq!(spinor_rows, 0);
        assert_eq!(distinct_probes.len(), EXPECTED_DISTINCT_PROBES);
        // The other-wave labels are a separate grammar: none of them resolves
        // to a record of its own space group.
        for (sg, label) in wave_labels {
            assert!(
                !query::irreps_of(sg).iter().any(|record| record.ml == label),
                "other-wave label {label} unexpectedly resolves in space group {sg}"
            );
        }
    }

    #[test]
    fn every_space_group_has_a_unique_trivial_child_record() {
        for sg in 1..=230u8 {
            let trivial = find_trivial_child(sg)
                .unwrap_or_else(|| panic!("space group {sg} has no trivial child record"));
            assert!(is_gamma(trivial), "trivial child of {sg} is not at Gamma");
            assert_eq!(trivial.dim, 1);
            assert!(!trivial.spinor);
        }
    }

    #[test]
    fn ordinal_13345_answers_every_probe_with_constructed_targets() {
        let (exit_code, counts, tsv) = audit_ordinal(13345, true);
        assert_eq!(counts.probes.total, 31);
        assert_eq!(
            counts.probes.attempted(),
            31,
            "every probe of a valid embedding must be answered, geometry-zero or not"
        );
        assert_eq!(counts.probes.uncomputed(), 0);
        assert_eq!(counts.probes.error, 0, "no engine inconsistency expected");
        assert!(
            counts.probes.partition_error().is_none(),
            "the probe partition must tile"
        );
        assert_eq!(
            counts.probes.full_success, 31,
            "R4 batch 1 answers every folded star of this context: the ones with \
             pinned child rows from the table and the trivial-little-co-group ones \
             from the constructed Bloch phase"
        );
        assert_eq!(
            counts.probes.identity_only, 0,
            "the five W probes used to be identity-only and are now complete"
        );
        assert_eq!(counts.probes.missing, 0);
        assert_eq!(
            exit_code, 0,
            "every probe has an exact answer, so --require-complete cannot fail"
        );
        assert_eq!(
            counts.verdict(Gates {
                complete: true,
                w_complete: true,
                ..Gates::default()
            }),
            Verdict::Clean,
            "the identity gate accepts the exact identity-only answers"
        );
        // This context has no gap left, so even the full-decomposition gate and
        // the w gate pass at this scope.
        assert_eq!(counts.verdict(all_gates()), Verdict::Clean);
        assert_eq!(counts.exit_code(all_gates()), 0);

        let rows = emitted_probe_rows(&tsv);
        assert_eq!(rows.len(), counts.probes.total, "one row per scalar probe");
        assert!(
            rows.iter().all(|row| !row.detail.contains("identity_only")),
            "no probe of this context is identity-only any more"
        );
        // The five W probes carry no Gamma folded star, so their absent zeroes
        // come from the engine and not from the geometry filter.
        for row in rows.iter().filter(|row| row.probe.starts_with('W')) {
            assert_eq!(
                row.computed, "0",
                "probe {} has no Gamma folded star, so the trivial content is zero",
                row.probe
            );
        }
        assert_eq!(counts.entries.passed, 8);
        assert_eq!(counts.entries.mismatch, 0);
        assert_eq!(counts.entries.unique, 8);
        assert!(
            rows.iter()
                .any(|row| row.detail.contains("geometry=reject")),
            "the independent geometry evidence must be emitted"
        );
    }

    /// The pinned baseline in miniature: ordinal 13346 carries five
    /// identity-only probes (the parametric `B`-line stars of the R3
    /// `parameterized_source` batch), so the full-decomposition gate fails there
    /// while the identity gate keeps passing.  Ordinal 13345 is the batch-1
    /// contrast: identical child and embedding, and no gap left.
    ///
    /// This is the guard against reporting "identity content is closed" as
    /// "full decomposition is closed".
    #[test]
    fn the_full_decomposition_gate_accepts_the_closed_batches() {
        // The witnesses of the three R4 batches: 13345 (trivial co-groups,
        // batch 1), 13346 (two-fold co-groups, batch 2a) and 3988 (four-element
        // co-group with one omega-regular class, batch 2b).  All three are
        // complete now, so the full-decomposition gate passes at their scope.
        // The gate's own separation logic is pinned by the synthetic tallies in
        // the unit tests below (`identity_only: 1` fixtures).
        for (ordinal, probes) in [(13345, 31), (13346, 31), (3988, 12)] {
            let (exit_code, counts, _) =
                audit_ordinal_gates(ordinal, gates_with_full_decomposition());
            assert_eq!(exit_code, 0, "ordinal {ordinal} has no gap left");
            assert_eq!(counts.probes.full_success, probes, "ordinal {ordinal}");
            assert_eq!(counts.probes.identity_only, 0, "ordinal {ordinal}");
        }

        // SG 16 R2 P1 -> #22 decomposes all 32 probes, so the same gate passes
        // there -- for that record's scope only.
        let (covered_exit, covered, _) = audit_ordinal_gates(314, gates_with_full_decomposition());
        assert_eq!(covered_exit, 0);
        assert_eq!(covered.probes.full_success, 32);
        assert_eq!(covered.probes.identity_only, 0);
        assert_eq!(covered.full_decomposition_incomplete(), 0);
    }

    #[test]
    fn ordinal_12400_decomposes_all_forty_sources_not_only_the_eligible_ones() {
        let (exit_code, counts, tsv) = audit_ordinal(12400, true);
        assert_eq!(counts.probes.total, 40);
        assert_eq!(
            counts.probes.attempted(),
            40,
            "the geometry filter must not skip decomposition"
        );
        assert_eq!(counts.probes.full_success, 40);
        assert_eq!(counts.probes.identity_only, 0);
        assert_eq!(counts.probes.missing, 0);
        assert_eq!(counts.probes.error, 0);
        assert_eq!(counts.probes.uncomputed(), 0);
        assert_eq!(counts.geometry_reject_absent, 30);
        assert_eq!(counts.geometry_eligible_absent, 7);
        assert_eq!(counts.absent_zero(), 37);
        assert_eq!(counts.geometry_contradictions, 0);
        assert_eq!(counts.entries.passed, 3);
        assert_eq!(exit_code, 0);
        assert_eq!(
            counts.verdict(Gates {
                complete: true,
                ..Gates::default()
            }),
            Verdict::Clean
        );

        let rows = emitted_probe_rows(&tsv);
        assert_eq!(rows.len(), 40);
        for row in &rows {
            assert!(
                !row.computed.is_empty(),
                "probe {} has no full result: {}",
                row.probe,
                row.status
            );
        }
        let rejected: Vec<&EmittedProbeRow> = rows
            .iter()
            .filter(|row| row.detail.contains("geometry=reject"))
            .collect();
        assert_eq!(rejected.len(), 30);
        for row in rejected {
            assert_eq!(row.status, "absent_zero");
            assert_eq!(row.computed, "0");
        }
    }

    /// The emitter must give every probe a row even when no engine call is
    /// possible, and those rows must be counted `uncomputed`, never silently
    /// skipped.
    ///
    /// The removed witness for this invariant was ordinal 0, whose embedding
    /// used to fail.  After the task-9 embedding fixes all 15,239 pinned records
    /// embed, so the synthetic tally below drives the same emission path with
    /// the shape [`Auditor::audit_record`] builds for it, and the last block
    /// pins the repaired witness instead of leaving the invariant untested.
    #[test]
    fn probes_without_an_embedding_are_reported_uncomputed_never_skipped() {
        let subgroup = sg16_r2_p1();
        let size = subgroup_size(subgroup.record.basis).ok();
        let ctx = RecordContext {
            ordinal: subgroup.ordinal,
            sg: 16,
            child_sg: u8::try_from(subgroup.record.sg).expect("child sg fits u8"),
            subgroup: &subgroup,
            size,
        };
        let sink = SharedSink::default();
        let mut auditor = Auditor::new(Options::default(), Box::new(sink.clone()));
        let probes: Vec<&'static IrrepRecord> = query::irreps_of(16)
            .iter()
            .filter(|record| !record.spinor)
            .take(8)
            .collect();
        assert_eq!(probes.len(), 8);
        let mut tally = RecordTally {
            records: 1,
            ..RecordTally::default()
        };
        tally.probes.total = probes.len();
        for probe in &probes {
            tally.probes.uncomputed_embedding += 1;
            auditor.emit_probe(
                &ctx,
                "absent",
                probe,
                None,
                None,
                "uncomputed_embedding_unavailable",
                "uncomputed: no valid embedding for this isotropy record",
            );
        }
        auditor.counts.absorb(&tally);
        assert_eq!(auditor.counts.probes.total, 8);
        assert_eq!(auditor.counts.probes.attempted(), 0);
        assert_eq!(auditor.counts.probes.uncomputed(), 8);
        assert!(auditor.counts.probes.partition_error().is_none());
        assert_eq!(
            auditor.counts.accounting_error(),
            None,
            "one emitted row per probe"
        );
        assert_eq!(
            auditor.counts.exit_code(Gates {
                complete: true,
                ..Gates::default()
            }),
            2,
            "--require-complete must report the uncomputed probes"
        );
        assert_eq!(auditor.counts.verdict(Gates::default()), Verdict::Clean);
        let rows = emitted_probe_rows(&sink.text());
        assert_eq!(
            rows.len(),
            8,
            "every probe must have a row even without an embedding"
        );
        for row in &rows {
            assert_eq!(row.status, "uncomputed_embedding_unavailable");
            assert!(row.computed.is_empty());
        }

        // The old witness: ordinal 0 now embeds and every probe is decomposed.
        let (exit_code, counts, _) = audit_ordinal(0, true);
        assert_eq!(counts.records_embedding_ok, 1);
        assert_eq!(counts.records_embedding_failed, 0);
        assert_eq!(counts.probes.total, 8);
        assert_eq!(counts.probes.full_success, 8);
        assert_eq!(counts.probes.uncomputed_embedding, 0);
        assert_eq!(exit_code, 0);
    }

    // ── Accounting ───────────────────────────────────────────────────────────

    #[test]
    fn injected_accounting_inconsistencies_are_rejected() {
        // Probe partition: one dropped probe.
        let mut counts = Counts::default();
        counts.probes.total = 5;
        counts.probes.full_success = 3;
        counts.probes.missing = 1;
        counts.probes.error = 1;
        assert!(counts.probes.partition_error().is_none());
        counts.probes.full_success = 2;
        let error = counts
            .probes
            .partition_error()
            .expect("a dropped probe must be detected");
        assert!(error.contains("probe partition"));
        assert!(counts.accounting_error().is_some());
        // A probe whose categories tile but whose row was never emitted is
        // still a dropped probe.
        let mut counts = Counts::default();
        counts.probes.total = 1;
        counts.probes.full_success = 1;
        let error = counts
            .accounting_error()
            .expect("a missing probe row must be detected");
        assert!(error.contains("emitted probe rows"));
        counts.probe_rows_emitted = 1;
        assert!(counts.accounting_error().is_none());
        // A genuine engine error is exit 1 even with --require-complete.
        let mut counts = Counts::default();
        counts.probes.total = 1;
        counts.probes.error = 1;
        assert_eq!(counts.verdict(Gates::default()), Verdict::Inconsistent);
        assert_eq!(
            counts.exit_code(Gates {
                complete: true,
                ..Gates::default()
            }),
            1
        );
        // Missing data alone is incomplete, and only under --require-complete.
        let mut counts = Counts::default();
        counts.probes.total = 2;
        counts.probes.full_success = 1;
        counts.probes.missing = 1;
        assert_eq!(counts.verdict(Gates::default()), Verdict::Clean);
        assert_eq!(
            counts.verdict(Gates {
                complete: true,
                ..Gates::default()
            }),
            Verdict::Incomplete
        );
        assert_eq!(
            counts.exit_code(Gates {
                complete: true,
                ..Gates::default()
            }),
            2
        );
        // Uncomputed probes are incompleteness, never silent coverage.
        let mut counts = Counts::default();
        counts.probes.total = 3;
        counts.probes.full_success = 1;
        counts.probes.uncomputed_embedding = 2;
        assert_eq!(
            counts.verdict(Gates {
                complete: true,
                ..Gates::default()
            }),
            Verdict::Incomplete
        );
        // The other-wave rows are a scope of their own: uncomputed w rows are
        // incomplete under the w gate, not under the ordinary one.
        let mut counts = Counts::default();
        counts.probes.total = 1;
        counts.probes.full_success = 1;
        counts.w_entries = 4;
        assert_eq!(counts.verdict(Gates::default()), Verdict::Clean);
        assert_eq!(
            counts.verdict(Gates {
                complete: true,
                ..Gates::default()
            }),
            Verdict::Clean
        );
        assert_eq!(
            counts.verdict(Gates {
                w_complete: true,
                ..Gates::default()
            }),
            Verdict::Incomplete
        );
        assert_eq!(
            counts.exit_code(Gates {
                w_complete: true,
                ..Gates::default()
            }),
            2
        );
        counts.w_computed = 4;
        assert_eq!(
            counts.verdict(Gates {
                w_complete: true,
                ..Gates::default()
            }),
            Verdict::Clean
        );
        // A w computation error is a hard failure, whatever the gates say.
        let mut counts = Counts::default();
        counts.w_entries = 4;
        counts.w_computed = 4;
        counts.w_engine_error = 1;
        assert_eq!(counts.verdict(Gates::default()), Verdict::Inconsistent);
        assert_eq!(
            counts.exit_code(Gates {
                w_complete: true,
                ..Gates::default()
            }),
            1
        );
        // Positive-entry partition: a dropped stored row is caught.
        let mut entries = EntryCounts::default();
        entries.entries = 2;
        entries.unique = 2;
        entries.passed = 1;
        entries.embedding_unavailable = 1;
        assert!(entries.partition_error().is_none());
        entries.embedding_unavailable = 0;
        assert!(entries.partition_error().is_some());
        // A double-counted entry total is caught.
        let mut entries = EntryCounts::default();
        entries.entries = 3;
        entries.unique = 2;
        entries.passed = 2;
        entries.duplicate_same = 2;
        assert!(entries.partition_error().is_some());
        // Geometry contradictions are hard failures, not incompleteness.
        let mut counts = Counts::default();
        counts.probes.total = 1;
        counts.probes.full_success = 1;
        counts.geometry_contradictions = 1;
        assert_eq!(counts.verdict(Gates::default()), Verdict::Inconsistent);
        assert_eq!(
            counts.exit_code(Gates {
                complete: true,
                ..Gates::default()
            }),
            1
        );
    }

    #[test]
    fn repeated_stored_rows_are_deduplicated_not_summed() {
        let mut map = BTreeMap::new();
        let mut entries = EntryCounts::default();
        let row = StoredRow {
            frequency: 2,
            domain: 1,
            ml: "GM1",
        };
        insert_stored_row(&mut map, 7, row, &mut entries).unwrap();
        insert_stored_row(&mut map, 7, row, &mut entries).unwrap();
        assert_eq!(entries.duplicate_same, 1);
        assert_eq!(map.len(), 1);
        assert_eq!(map.values().filter(|entry| entry.frequency == 2).count(), 1);
        let conflicting = StoredRow {
            frequency: 4,
            ..row
        };
        assert!(insert_stored_row(&mut map, 7, conflicting, &mut entries).is_err());
        assert_eq!(entries.duplicate_conflict, 1);
        assert_eq!(
            map[&7].frequency, 2,
            "the first copy wins, nothing is summed"
        );
    }

    #[test]
    fn frobenius_weight_uses_true_complex_component_dimensions() {
        // Ordinary, realification and distinct-component rows from the pinned
        // tables; each must report its complex component dimension sum.
        let mut saw_realification = false;
        let mut saw_distinct = false;
        for sg in 1..=230u8 {
            for record in query::irreps_of(sg) {
                if record.spinor {
                    continue;
                }
                let (sum, count) = complex_component_dimensions(record).unwrap();
                match record.source_identity() {
                    IrrepSourceIdentity::OrdinaryScalar { .. } => {
                        assert_eq!((sum, count), (record.dim, 1));
                    }
                    IrrepSourceIdentity::Compound { .. } => {
                        assert_eq!(count, 2);
                        if is_gamma(record) {
                            // A Gamma star is one arm, so the row dimension is
                            // the complex component sum.
                            assert_eq!(u16::from(sum), u16::from(record.dim));
                        } else {
                            // Other stars carry arm_count copies of the sum.
                            assert_eq!(u16::from(record.dim) % u16::from(sum), 0);
                        }
                        match record.compound_character_semantics() {
                            Some(CompoundCharacterSemantics::DistinctComponentSum) => {
                                saw_distinct = true;
                            }
                            Some(CompoundCharacterSemantics::ConjugateRealification) => {
                                saw_realification = true;
                            }
                            None => panic!("compound row {} has no semantics", record.ml),
                        }
                    }
                    IrrepSourceIdentity::Spin { .. } => {
                        panic!("spinor row in the scalar sweep")
                    }
                }
            }
        }
        assert!(saw_realification && saw_distinct);
    }

    #[test]
    fn geometry_zero_is_an_independent_check_never_a_skip() {
        // SG 16 R2 P1 -> #22: every probe is decomposed, and a geometry-rejected
        // probe must come back with exactly zero trivial terms.
        let (exit_code, counts, tsv) = audit_ordinal(314, true);
        assert_eq!(counts.probes.total, 32);
        assert_eq!(
            counts.probes.attempted(),
            32,
            "the geometry filter must not skip decomposition"
        );
        assert_eq!(counts.probes.full_success, 32);
        assert_eq!(counts.probes.uncomputed(), 0);
        assert_eq!(counts.probes.missing, 0);
        assert_eq!(counts.probes.error, 0);
        assert_eq!(counts.geometry_reject_absent, 24);
        assert_eq!(counts.geometry_eligible_absent, 6);
        assert_eq!(counts.absent_zero(), 30);
        assert_eq!(counts.geometry_contradictions, 0);
        assert_eq!(counts.entries.passed, 2);
        assert_eq!(exit_code, 0);
        assert_eq!(counts.verdict(all_gates()), Verdict::Clean);

        // The independent geometry filter must agree with the computed content
        // on all 32 probes: a rejected probe has exactly zero trivial terms and
        // every positive probe is eligible.  The eight eligible probes are the
        // six absent ones plus the two stored positives.
        let subgroup = sg16_r2_p1();
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let reciprocal = build_child_reciprocal(embedding.subgroup_sg()).unwrap();
        let trivial = find_trivial_child(embedding.subgroup_sg()).unwrap();
        let trivial_cir = match trivial.source_identity() {
            IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => cir_irnumber,
            _ => panic!("trivial child is not ordinary"),
        };
        let (mut probes, mut eligible, mut positive, mut rejected) =
            (0usize, 0usize, 0usize, 0usize);
        for probe in query::irreps_of(16).iter().filter(|record| !record.spinor) {
            probes += 1;
            let star = ScalarStar::new(probe).unwrap();
            let filter = probe_is_gamma_eligible(&star, &embedding, &reciprocal).unwrap();
            eligible += usize::from(filter);
            let result = match subduce_full_star_with_embedding(&subgroup, &embedding, probe) {
                Ok(result) => result,
                Err(error) => panic!("probe {} must be decomposed: {error}", probe.ml),
            };
            let computed: u32 = result
                .blocks()
                .iter()
                .flat_map(|block| block.targets())
                .filter(|target| {
                    target.irnumber == Some(trivial_cir) && target.dimension == trivial.dim
                })
                .map(|target| target.multiplicity)
                .sum();
            positive += usize::from(computed > 0);
            if !filter {
                rejected += 1;
                assert_eq!(
                    computed, 0,
                    "probe {} is geometry-rejected but has trivial terms",
                    probe.ml
                );
            }
            assert!(
                computed == 0 || filter,
                "probe {} is positive but the geometry filter rejects it",
                probe.ml
            );
        }
        assert_eq!(probes, 32);
        assert_eq!(eligible, 8);
        assert_eq!(positive, 2);
        assert_eq!(rejected, 24);

        // The row-level evidence: the 24 rejects are emitted with computed=0 and
        // every probe row carries its full result.
        let rows = emitted_probe_rows(&tsv);
        assert_eq!(rows.len(), 32);
        for row in &rows {
            assert!(
                !row.computed.is_empty(),
                "probe {} has no full result: {}",
                row.probe,
                row.status
            );
        }
        let rejected_rows: Vec<&EmittedProbeRow> = rows
            .iter()
            .filter(|row| row.detail.contains("geometry=reject"))
            .collect();
        assert_eq!(rejected_rows.len(), 24);
        for row in rejected_rows {
            assert_eq!(row.status, "absent_zero");
            assert_eq!(row.computed, "0");
        }
    }
}
