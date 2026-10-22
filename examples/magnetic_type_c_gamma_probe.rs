//! R31 finite Type-C Γ feasibility gate (probe only, no production API).
//!
//! This probe answers one bounded question: for three **fixed** two-dimensional
//! Γ sources of their parent space groups, do the existing magnetic isotropy
//! records expose a magnetic Γ corepresentation that is
//!
//! * Wigner **type C** with dimension **2**, and
//! * built from **two inequivalent, non-spinor, one-dimensional H irreps with
//!   real `±1` characters**?
//!
//! It is deliberately **not** a Type-C implementation and it does not touch the
//! production subduction entry points, which remain Type-A-only.  The fixed
//! sources are `SG 89 GM5`, `SG 156 GM3` and `SG 191 GM5+` (all CDML).
//!
//! # Population and provenance
//!
//! * A **row** is one row of the generated magnetic isotropy tables.  Every
//!   scalar irrep record of a parent contributes its own `magnetic_subgroups()`
//!   rows; two rows that share a UNI are still two records and are never
//!   merged.  The source label and the condensing label are independent
//!   identities: a row is never dropped because it "is" the source row.
//! * A **UNI** is the magnetic group of a row.  Its Γ corepresentations are
//!   read once, through the public partial summary
//!   [`magnetic_summary::magnetic_irrep_summary_by_uni_partial`].  The partial
//!   form matters: a failed source at a **non-Γ** k-point is counted separately
//!   and cannot turn a resolved Γ into "zero candidates"; a Γ failure is an
//!   explicit unresolved count, never a negative result.
//! * The corep's `source_irreps` are the **H** (unitary subgroup) source
//!   identities reported by the summary.  Their characters are re-read from the
//!   typed selected-arm row of the matching H record and required to be real `±1`
//!   on every operation.  The corep **label** is display text and is never used
//!   as an identity.  When a source is an individually named constituent of a
//!   compound (`CIR`) H record and has no standalone `IrrepRecord`, its typed
//!   constituent row is recovered from the compound view; it is not dropped and
//!   not relabelled.  A constituent whose row is complex is an explicit
//!   real-`±1`-domain exclusion.
//! * Spinor **parent** records have no generated magnetic isotropy table (the
//!   accessor itself refuses them).  They are counted as a named skip, so the
//!   row total is compared with the public row query rather than silently
//!   shrinking the inventory.
//!
//! # Layers
//!
//! 1. inventory: rows per parent, distinct UNIs, summary/Γ/errors;
//! 2. `type_c_dim2` coreps;
//! 3. the real-H-pair screen (non-spinor, dimension 1, real `±1`, unequal rows);
//! 4. (only when layer 3 yields candidates) the first exact-realising R24
//!    embedding, the certified source model under both `η`, the unitary
//!    restriction `trace A(h) = χ1(h)+χ2(h)` class by class through the R27
//!    raw-index → class binding, the a-twisted pairing
//!    `χ1(a0⁻¹ h a0) = χ2(h)` and `χ2(a0⁻¹ h a0) = χ1(h)`, and the matrix-level
//!    exchange of the two H eigen-spaces under every antiunitary class.  If
//!    layer 3 yields zero candidates the probe stops there and says so: no
//!    later layer is claimed to have run.
//!
//! The layer-4 source model is `verify_matrix_gamma_input` (the R25/R29 input
//! certifier), **not** `subduce_matrix_gamma`: no Type-C target is constructed
//! and the accepted production subduction remains type-A-only.
//!
//! # Exit semantics
//!
//! * `0` — every layer that applies completed, all conservation identities
//!   hold, and every candidate is either a witness or a definite
//!   incompatibility under the fixed embedding.  Zero candidates is a valid
//!   finite result and is printed as such.
//! * `1` — an internal invariant or the built-in self-check failed.  This is a
//!   bug or a corrupted reading, not a physical negative result.
//! * `2` — the fixed domain is **incomplete**: a UNI summary failed, Γ is
//!   missing, a Γ source is unresolved, a type-C candidate could not be
//!   screened because its H source/character data is unavailable, its R24
//!   embedding is absent/undecided, or the source model could not be verified.
//!   Such a run never reports "no candidate".
//!
//! `--self-check` runs only the synthetic screen checks: a positive type-C
//! dimension-2 real pair, and the named classifications of type-A, complex,
//! spinor, dimension-2, equivalent, missing-character and single-source cases.
//! `--tsv PATH` writes a stable, tab-separated candidate table.
//!
//! # Scope fences
//!
//! Fixed three parents, fixed three sources, existing tables only.  No setting
//! search, no new dependency, no fallback H label or character, no claim that a
//! zero is a statement about all magnetic groups or all settings.  The witness
//! is source/record/embedding compatibility for the fixed R24 family, not a
//! complete band assignment and not a proof about maps outside that family.

use std::collections::{BTreeMap, BTreeSet};
use std::process::ExitCode;

use num_complex::Complex64;

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::corep::CorepType;
use cryspglib::irrep::isotropy;
use cryspglib::irrep::magnetic_embedding::{
    self, MagneticOperation, MagneticOperationSet, ParentEmbedding, canonical_translation, compose,
    geometry_of, magnetic_operations, parent_embedding_verdict,
};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::{
    MatrixGammaModel, MatrixGammaRequest, verify_matrix_gamma_input,
};
use cryspglib::irrep::magnetic_summary::{self, MagneticCorepSummary};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::types::{
    CharacterRow, CharacterViewError, CompoundSelectedArmCharacter, SeitzOperation,
};

/// The published R30 baseline.
const BASELINE: &str = "cf0a4baab95ce999cb96084c9571b604d0596807";

/// Γ as the magnetic summary encodes k-point coordinates.
const GAMMA: (i8, i8, i8, i8) = (0, 0, 0, 1);

/// Tolerance of the finite floating-point `±1` classification.  This is a
/// tolerance decision on the generated `f64` values, not an exact proof of
/// their binary representation.
const GRID_TOLERANCE: f64 = 1e-9;

/// Quantization of Seitz translations when two typed H character rows are
/// compared operation by operation.  The generated rows carry about five
/// decimal places, so `1e-6` is finer than the data and coarse enough to be a
/// stable key.
const OPERATION_TRANSLATION_QUANTIZATION: f64 = 1.0e6;

/// The fixed sources.
const FIXED_SOURCES: [(u8, &str); 3] = [(89, "GM5"), (156, "GM3"), (191, "GM5+")];

fn main() -> ExitCode {
    let mut tsv_path: Option<String> = None;
    let mut self_check_only = false;
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--self-check" => self_check_only = true,
            "--tsv" => {
                index += 1;
                let Some(path) = args.get(index) else {
                    eprintln!("--tsv needs a path");
                    return ExitCode::from(1);
                };
                tsv_path = Some(path.clone());
            }
            other => {
                eprintln!("unknown argument {other:?}");
                eprintln!("usage: magnetic_type_c_gamma_probe [--self-check] [--tsv PATH]");
                return ExitCode::from(1);
            }
        }
        index += 1;
    }

    let self_check_violations = self_check();
    for violation in &self_check_violations {
        eprintln!("SELF-CHECK FAILED: {violation}");
    }
    if self_check_only {
        println!("self_check violations={}", self_check_violations.len());
        return if self_check_violations.is_empty() {
            ExitCode::from(0)
        } else {
            ExitCode::from(1)
        };
    }

    let report = scan();
    print_report(&report);
    if let Some(path) = tsv_path {
        match std::fs::write(&path, report.candidate_tsv()) {
            Ok(()) => println!("tsv {path}"),
            Err(error) => {
                eprintln!("cannot write {path}: {error}");
                return ExitCode::from(1);
            }
        }
    }

    if !self_check_violations.is_empty() || !report.invariant_violations.is_empty() {
        return ExitCode::from(1);
    }
    if report.incomplete() {
        ExitCode::from(2)
    } else {
        ExitCode::from(0)
    }
}

// ---------------------------------------------------------------------------
// Inventory
// ---------------------------------------------------------------------------

/// One magnetic isotropy row of one parent.
#[derive(Clone, Copy)]
struct Row {
    parent_sg: u8,
    /// CDML label of the condensing (owning) parent irrep.
    condensing: &'static str,
    /// Position in that irrep's magnetic isotropy list.
    direction_index: usize,
    /// Magnetic space group of the row.
    uni: usize,
}

/// The fixed source record, validated before the scan.
struct FixedSource {
    sg: u8,
    label: &'static str,
    spinor: bool,
    selected_dimension: usize,
}

fn fixed_sources() -> Result<Vec<FixedSource>, String> {
    let mut out = Vec::new();
    for (sg, label) in FIXED_SOURCES {
        let matches = query::find_irreps(sg, label, LabelConvention::Cdml);
        let [ir] = matches.as_slice() else {
            return Err(format!(
                "fixed source SG {sg} {label} resolved to {} records",
                matches.len()
            ));
        };
        if ir.spinor {
            return Err(format!("fixed source SG {sg} {label} is spinor"));
        }
        if ir.dim != 2 {
            return Err(format!(
                "fixed source SG {sg} {label} has generated dim {}, expected 2",
                ir.dim
            ));
        }
        let dimension = match ir.ordinary_scalar_selected_arm_block_trace() {
            Ok(row) => row.dimension(),
            Err(error) => {
                return Err(format!(
                    "fixed source SG {sg} {label} has no typed selected-arm row: {error}"
                ));
            }
        };
        if dimension != 2 {
            return Err(format!(
                "fixed source SG {sg} {label} has selected-arm dimension {dimension}, expected 2"
            ));
        }
        out.push(FixedSource {
            sg,
            label,
            spinor: false,
            selected_dimension: dimension,
        });
    }
    Ok(out)
}

/// Enumerate every magnetic isotropy row of a parent.
///
/// Spinor records are inspected too: they have no magnetic isotropy data in the
/// current tables, but the count is checked against the public row query so a
/// future non-empty spinor table cannot make the inventory silently smaller.
fn rows_of_parent(sg: u8) -> (Vec<Row>, usize, usize) {
    let mut rows = Vec::new();
    let mut spinor_records_skipped = 0usize;
    for irrep in query::irreps_of(sg) {
        if irrep.spinor {
            // The generated accessor itself refuses spinor isotropy data, and
            // the public row query skips these records too.  They are counted
            // as a named boundary instead of being read through a panic.
            spinor_records_skipped += 1;
            continue;
        }
        for (direction_index, record) in irrep.magnetic_subgroups().iter().enumerate() {
            rows.push(Row {
                parent_sg: sg,
                condensing: irrep.ml,
                direction_index,
                uni: record.mag_sg,
            });
        }
    }
    let public_count = query::magnetic_isotropy_subgroups_of(sg).len();
    (rows, public_count, spinor_records_skipped)
}

// ---------------------------------------------------------------------------
// Γ corepresentations of one UNI
// ---------------------------------------------------------------------------

/// Why one H source row could not be used as a real `±1` one-dimensional row.
#[derive(Clone, Debug, PartialEq, Eq)]
enum CharacterProblem {
    MissingHSource,
    AmbiguousHSource,
    DimensionMismatch { summary: u8, typed: usize },
    Unavailable(&'static str),
    NonFinite,
    NonReal,
    OffGrid,
    IdentityNotOne,
    ZeroDimension,
    MalformedRow,
}

impl CharacterProblem {
    fn name(&self) -> String {
        match self {
            CharacterProblem::MissingHSource => "MissingHSource".to_string(),
            CharacterProblem::AmbiguousHSource => "AmbiguousHSource".to_string(),
            CharacterProblem::DimensionMismatch { summary, typed } => {
                format!("DimensionMismatch({summary}!={typed})")
            }
            CharacterProblem::Unavailable(variant) => format!("Unavailable({variant})"),
            CharacterProblem::NonFinite => "NonFinite".to_string(),
            CharacterProblem::NonReal => "NonReal".to_string(),
            CharacterProblem::OffGrid => "OffGrid".to_string(),
            CharacterProblem::IdentityNotOne => "IdentityNotOne".to_string(),
            CharacterProblem::ZeroDimension => "ZeroDimension".to_string(),
            CharacterProblem::MalformedRow => "MalformedRow".to_string(),
        }
    }
}

/// A one-dimensional real `±1` character row, keyed by a quantized operation identity.
#[derive(Clone, Debug, PartialEq, Eq)]
struct RealOneDimRow {
    values: BTreeMap<([i32; 9], [i64; 3]), i8>,
}

impl RealOneDimRow {
    /// Operation-by-operation equality, including the operation domains.
    fn equivalent(&self, other: &Self) -> bool {
        self.values == other.values
    }
}

/// One H source as reported by the summary, plus the re-read typed row.
#[derive(Clone)]
struct SourceFacts {
    sg: u8,
    label: &'static str,
    /// The verified BC label of the same source record, when it has one.
    bc: Option<String>,
    spinor: bool,
    dimension: u8,
    character: Result<RealOneDimRow, CharacterProblem>,
}

fn source_facts(source: &cryspglib::irrep::magnetic_summary::SourceIrrepSummary) -> SourceFacts {
    let mut facts = SourceFacts {
        sg: source.sg,
        label: source.ml,
        bc: source.labels().bc,
        spinor: source.spinor,
        dimension: source.dim,
        character: Err(CharacterProblem::MissingHSource),
    };
    let matches = query::find_irreps(source.sg, source.ml, LabelConvention::Cdml);
    let row = match matches.as_slice() {
        [ir] => {
            if ir.spinor != source.spinor {
                facts.character = Err(CharacterProblem::Unavailable("spinor_flag_mismatch"));
                return facts;
            }
            match ir.ordinary_scalar_selected_arm_block_trace() {
                Ok(row) => row,
                Err(CharacterViewError::NotApplicable) if ir.compound_metadata().is_some() => {
                    match compound_constituent_row(source.sg, source.ml, source.dim) {
                        Ok(Some(row)) => row,
                        Ok(None) => return facts,
                        Err(problem) => {
                            facts.character = Err(problem);
                            return facts;
                        }
                    }
                }
                Err(error) => {
                    facts.character = Err(CharacterProblem::Unavailable(character_variant(&error)));
                    return facts;
                }
            }
        }
        [] => match compound_constituent_row(source.sg, source.ml, source.dim) {
            Ok(Some(row)) => row,
            Ok(None) => return facts,
            Err(problem) => {
                facts.character = Err(problem);
                return facts;
            }
        },
        _ => {
            facts.character = Err(CharacterProblem::AmbiguousHSource);
            return facts;
        }
    };
    if row.dimension() != source.dim as usize {
        facts.character = Err(CharacterProblem::DimensionMismatch {
            summary: source.dim,
            typed: row.dimension(),
        });
        return facts;
    }
    facts.character = real_one_dim_row(&row);
    facts
}

fn character_variant(error: &CharacterViewError) -> &'static str {
    match error {
        CharacterViewError::NotApplicable => "NotApplicable",
        CharacterViewError::MissingData => "MissingData",
        CharacterViewError::InvalidDimension => "InvalidDimension",
        CharacterViewError::NonFiniteData => "NonFiniteData",
        CharacterViewError::MisalignedStorage => "MisalignedStorage",
        CharacterViewError::LengthMismatch => "LengthMismatch",
        CharacterViewError::InvalidOperationIndex => "InvalidOperationIndex",
        CharacterViewError::InvalidComponent => "InvalidComponent",
        CharacterViewError::DuplicateOperationIndex => "DuplicateOperationIndex",
        CharacterViewError::NonUniqueIdentity => "NonUniqueIdentity",
        CharacterViewError::InvalidIdentityDimension => "InvalidIdentityDimension",
        CharacterViewError::DimensionMismatch => "DimensionMismatch",
        CharacterViewError::OperationOrderMismatch => "OperationOrderMismatch",
    }
}

/// A compound record is a real physical source whose irreducible CIR
/// constituents may not have their own standalone `IrrepRecord` row.  The
/// summary reports such a constituent by its CIR label; recover the constituent
/// row from the owning compound record's typed view.  A label that occurs in
/// two compound records (or twice in one view) is ambiguous, not picked.
fn compound_constituent_row(
    sg: u8,
    label: &str,
    dimension: u8,
) -> Result<Option<CharacterRow>, CharacterProblem> {
    let mut found: Option<CharacterRow> = None;
    for irrep in query::irreps_of(sg) {
        if irrep.spinor || irrep.compound_metadata().is_none() {
            continue;
        }
        let Ok(view) = irrep.compound_selected_arm_view() else {
            continue;
        };
        let mut candidates = Vec::new();
        match &view {
            CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => {
                candidates.push(seed);
            }
            CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } => {
                candidates.push(first);
                candidates.push(second);
            }
        }
        for candidate in candidates {
            if candidate.label == label && candidate.dimension as u8 == dimension {
                if found.is_some() {
                    return Err(CharacterProblem::AmbiguousHSource);
                }
                found = Some(candidate.row.clone());
            }
        }
    }
    Ok(found)
}

fn real_one_dim_row(row: &CharacterRow) -> Result<RealOneDimRow, CharacterProblem> {
    if row.dimension() == 0 {
        return Err(CharacterProblem::ZeroDimension);
    }
    if row.dimension() != 1 {
        return Err(CharacterProblem::DimensionMismatch {
            summary: 1,
            typed: row.dimension(),
        });
    }
    if row.len() != row.operations().len() || row.is_empty() {
        return Err(CharacterProblem::MalformedRow);
    }
    let pairs: Vec<(Complex64, SeitzOperation)> = row
        .values()
        .iter()
        .copied()
        .zip(row.operations().iter().copied())
        .collect();
    real_one_dim_row_from_pairs(&pairs)
}

/// The numeric half of [`real_one_dim_row`], separated so the self-check can
/// exercise the real `f64 -> ±1 / NonFinite / NonReal / OffGrid` path without
/// constructing a private `CharacterRow`.
fn real_one_dim_row_from_pairs(
    pairs: &[(Complex64, SeitzOperation)],
) -> Result<RealOneDimRow, CharacterProblem> {
    if pairs.is_empty() {
        return Err(CharacterProblem::MalformedRow);
    }
    let mut values = BTreeMap::new();
    let mut saw_identity = false;
    for (value, operation) in pairs {
        let real = grid_sign(value)?;
        let key = operation_key(operation)?;
        if values.insert(key, real).is_some() {
            return Err(CharacterProblem::MalformedRow);
        }
        if is_identity(operation) {
            saw_identity = true;
            if real != 1 {
                return Err(CharacterProblem::IdentityNotOne);
            }
        }
    }
    if !saw_identity {
        return Err(CharacterProblem::MalformedRow);
    }
    Ok(RealOneDimRow { values })
}

/// Classify one generated character value on the finite `±1` grid.
///
/// * non-finite real or imaginary part → [`CharacterProblem::NonFinite`]
///   (unresolved, never a domain exclusion);
/// * finite but with a non-zero imaginary part → [`CharacterProblem::NonReal`];
/// * finite real value that is not within tolerance of `+1`/`-1` →
///   [`CharacterProblem::OffGrid`].
fn grid_sign(value: &Complex64) -> Result<i8, CharacterProblem> {
    if !value.re.is_finite() || !value.im.is_finite() {
        return Err(CharacterProblem::NonFinite);
    }
    if value.im.abs() > GRID_TOLERANCE {
        return Err(CharacterProblem::NonReal);
    }
    if (value.re - 1.0).abs() <= GRID_TOLERANCE {
        Ok(1)
    } else if (value.re + 1.0).abs() <= GRID_TOLERANCE {
        Ok(-1)
    } else {
        Err(CharacterProblem::OffGrid)
    }
}

fn operation_key(operation: &SeitzOperation) -> Result<([i32; 9], [i64; 3]), CharacterProblem> {
    let mut translation = [0i64; 3];
    for (index, &value) in operation.translation.iter().enumerate() {
        let quantized = value * OPERATION_TRANSLATION_QUANTIZATION;
        if !quantized.is_finite() || quantized.abs() > i64::MAX as f64 {
            return Err(CharacterProblem::MalformedRow);
        }
        translation[index] = quantized.round() as i64;
    }
    Ok((operation.rotation, translation))
}

fn is_identity(operation: &SeitzOperation) -> bool {
    operation.rotation == [1, 0, 0, 0, 1, 0, 0, 0, 1]
        && operation
            .translation
            .iter()
            .all(|value| (value - value.round()).abs() <= 1e-8)
}

/// One character column of a Γ corep, with the summary's raw first-Hall index
/// and the data-Hall rotation/translation in which its character was computed.
#[derive(Clone)]
struct GammaColumn {
    /// Index into the raw `magnetic_operations(uni).operations` list (R27).
    raw_index: usize,
    rotation: [i32; 9],
    translation: [f64; 3],
    time_reversal: bool,
}

/// One Γ corep reduced to the facts the screen and the pairing checks need.
#[derive(Clone)]
struct GammaCorep {
    label: String,
    corep_type: CorepType,
    dim: usize,
    sources: Vec<SourceFacts>,
    columns: Vec<GammaColumn>,
    /// Character value of each column; `None` only where the summary has no
    /// constructed value (type-A antiunitary columns).
    characters: Vec<Option<Complex64>>,
}

/// The Γ reading of one UNI.
enum UniGamma {
    Resolved {
        unitary_sg: u8,
        coreps: Vec<GammaCorep>,
        /// Unresolved coreps whose summary k-label is Γ.
        unresolved_gamma: usize,
        /// Unresolved coreps at other k-points: reported, never mixed into Γ.
        unresolved_other: usize,
        gamma_columns: usize,
    },
    SummaryFailed(String),
    GammaMissing,
}

fn gamma_of(uni: usize) -> UniGamma {
    let partial = match magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni) {
        Ok(partial) => partial,
        Err(error) => return UniGamma::SummaryFailed(format!("{error:?}")),
    };
    let Some(kpoint) = partial
        .summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.coords == GAMMA)
    else {
        return UniGamma::GammaMissing;
    };
    let unresolved_gamma = partial
        .unresolved_coreps
        .iter()
        .filter(|entry| entry.k_label == "GM")
        .count();
    let unresolved_other = partial.unresolved_coreps.len() - unresolved_gamma;
    let coreps = kpoint
        .coreps
        .iter()
        .map(|corep| gamma_corep(corep, &kpoint.operations))
        .collect();
    UniGamma::Resolved {
        unitary_sg: partial.summary.unitary_sg,
        coreps,
        unresolved_gamma,
        unresolved_other,
        gamma_columns: kpoint.operations.len(),
    }
}

fn gamma_corep(
    corep: &MagneticCorepSummary,
    operations: &[cryspglib::irrep::magnetic_summary::MagneticLittleGroupOperation],
) -> GammaCorep {
    let columns = operations
        .iter()
        .map(|operation| GammaColumn {
            raw_index: operation.magnetic_operation_index,
            rotation: [
                operation.rotation[0][0],
                operation.rotation[0][1],
                operation.rotation[0][2],
                operation.rotation[1][0],
                operation.rotation[1][1],
                operation.rotation[1][2],
                operation.rotation[2][0],
                operation.rotation[2][1],
                operation.rotation[2][2],
            ],
            translation: operation.translation,
            time_reversal: operation.time_reversal,
        })
        .collect();
    GammaCorep {
        label: corep.label.clone(),
        corep_type: corep.corep_type,
        dim: corep.dim,
        sources: corep.source_irreps.iter().map(source_facts).collect(),
        columns,
        characters: corep.characters.clone(),
    }
}

// ---------------------------------------------------------------------------
// Layer 3: the real-H-pair screen
// ---------------------------------------------------------------------------

/// Why a type-C dimension-2 corep is **not** in the fixed real-H-pair domain.
#[derive(Clone, Debug, PartialEq, Eq)]
enum PairExclusion {
    Spinor(usize),
    Dimension { index: usize, dimension: u8 },
    NonReal(usize),
}

/// Why a type-C dimension-2 corep could **not be screened**.  This is not a
/// negative result and is never counted as one.
#[derive(Clone, Debug, PartialEq, Eq)]
enum PairUnresolved {
    NotTypeCDim2 {
        corep_type: CorepType,
        dim: usize,
    },
    SourceCount(usize),
    Character {
        index: usize,
        problem: CharacterProblem,
    },
    EquivalentPair,
    SameIdentity,
}

impl PairUnresolved {
    fn name(&self) -> String {
        match self {
            PairUnresolved::NotTypeCDim2 { corep_type, dim } => {
                format!("NotTypeCDim2({corep_type:?},{dim})")
            }
            PairUnresolved::SourceCount(count) => format!("SourceCount({count})"),
            PairUnresolved::Character { index, problem } => {
                format!("Character[{index}]({})", problem.name())
            }
            PairUnresolved::EquivalentPair => "EquivalentPair".to_string(),
            PairUnresolved::SameIdentity => "SameIdentity".to_string(),
        }
    }
}

/// The result of the single real-H-pair screen.
enum PairOutcome {
    Candidate { first: usize, second: usize },
    Excluded(PairExclusion),
    Unresolved(PairUnresolved),
}

/// The single screen used both by the real scan and by the self-check.
fn screen_pair(corep_type: CorepType, dim: usize, sources: &[SourceFacts]) -> PairOutcome {
    if corep_type != CorepType::C || dim != 2 {
        return PairOutcome::Unresolved(PairUnresolved::NotTypeCDim2 { corep_type, dim });
    }
    if sources.len() != 2 {
        return PairOutcome::Unresolved(PairUnresolved::SourceCount(sources.len()));
    }
    for (index, source) in sources.iter().enumerate() {
        if source.spinor {
            return PairOutcome::Excluded(PairExclusion::Spinor(index));
        }
        if source.dimension != 1 {
            return PairOutcome::Excluded(PairExclusion::Dimension {
                index,
                dimension: source.dimension,
            });
        }
    }
    if (sources[0].sg, sources[0].label) == (sources[1].sg, sources[1].label) {
        return PairOutcome::Unresolved(PairUnresolved::SameIdentity);
    }
    let mut rows = Vec::with_capacity(2);
    for (index, source) in sources.iter().enumerate() {
        match &source.character {
            Ok(row) => rows.push(row),
            Err(CharacterProblem::NonReal) | Err(CharacterProblem::OffGrid) => {
                return PairOutcome::Excluded(PairExclusion::NonReal(index));
            }
            Err(CharacterProblem::NonFinite) => {
                return PairOutcome::Unresolved(PairUnresolved::Character {
                    index,
                    problem: CharacterProblem::NonFinite,
                });
            }
            Err(problem) => {
                return PairOutcome::Unresolved(PairUnresolved::Character {
                    index,
                    problem: problem.clone(),
                });
            }
        }
    }
    if rows[0].equivalent(rows[1]) {
        return PairOutcome::Unresolved(PairUnresolved::EquivalentPair);
    }
    PairOutcome::Candidate {
        first: 0,
        second: 1,
    }
}

// ---------------------------------------------------------------------------
// The scan
// ---------------------------------------------------------------------------

/// Layer counters, kept separate so every report line is traceable.
#[derive(Default)]
struct Counters {
    rows: usize,
    public_rows: usize,
    spinor_records_skipped: usize,
    distinct_unis: usize,
    summary_failed: usize,
    gamma_missing: usize,
    unresolved_gamma: usize,
    unresolved_other: usize,
    gamma_columns: usize,
    first_summary_error: Option<String>,
    gamma_coreps: usize,
    type_a: usize,
    type_b: usize,
    type_c: usize,
    type_c_dim2_unique: usize,
    type_c_dim2_rows: usize,
    pair_candidates_unique: usize,
    pair_candidates_rows: usize,
    pair_excluded: BTreeMap<String, usize>,
    pair_unresolved: BTreeMap<String, usize>,
    candidate_rows: usize,
    per_parent: BTreeMap<u8, ParentCounters>,
}

#[derive(Default)]
struct ParentCounters {
    rows: usize,
    distinct_unis: usize,
    type_c_dim2_rows: usize,
    candidate_rows: usize,
}

/// The layer-4 evidence of one metadata candidate.
#[derive(Clone, Default)]
struct CandidateEvidence {
    /// `witness`, `not_established` or `unexecuted`.
    stage: String,
    /// Named reason when the candidate is not a witness.
    reason: String,
    map: Option<Mat3R>,
    shift: Option<[Rat; 3]>,
    embedding_kind: String,
    eta_plus: String,
    eta_minus: String,
    restriction: String,
    pairing: String,
    exchange_plus: String,
    exchange_minus: String,
    /// The layer-4 result is a definite incompatibility, not an unverified
    /// interface/embedding gap.
    definite: bool,
    /// A layer-4 prerequisite could not be evaluated.
    incomplete: bool,
}

/// One metadata candidate, after the row expansion.
#[derive(Clone)]
struct MetadataCandidate {
    row: Row,
    source_sg: u8,
    source_label: &'static str,
    /// The summary's unitary subgroup of the record's magnetic group.
    unitary_sg: u8,
    corep: GammaCorep,
    h_sources: [(u8, &'static str); 2],
    evidence: CandidateEvidence,
}

/// The complete metadata gate.
struct Report {
    fixed_sources: Vec<FixedSource>,
    counters: Counters,
    candidates: Vec<MetadataCandidate>,
    invariant_violations: Vec<String>,
    /// `true` when the metadata layer has no candidates and therefore no later
    /// layer was executed.
    stopped_at_metadata: bool,
    /// One line per screen that could not be decided, with its provenance.
    pair_unresolved_details: Vec<String>,
    layer4_executed: bool,
    witnesses: usize,
    candidate_incomplete: bool,
}

impl Report {
    fn incomplete(&self) -> bool {
        self.counters.summary_failed > 0
            || self.counters.gamma_missing > 0
            || self.counters.unresolved_gamma > 0
            || !self.counters.pair_unresolved.is_empty()
            || self.candidate_incomplete
    }

    fn candidate_tsv(&self) -> String {
        let mut out = String::from(
            "parent_sg\tsource_sg\tsource\tcondensing\tdirection_index\tuni\tunitary_sg\tcorep_label\
             \th1_sg\th1\th1_bc\th2_sg\th2\th2_bc\th1_signs\th2_signs\tcorep_type\tcorep_dim\
             \tstage\treason\tembedding_kind\tmap\tshift\
             \teta_plus\teta_minus\trestriction\tpairing\texchange_plus\texchange_minus\n",
        );
        for candidate in &self.candidates {
            let evidence = &candidate.evidence;
            let map = evidence
                .map
                .map(format_map)
                .unwrap_or_else(|| "-".to_string());
            let shift = evidence
                .shift
                .map(format_shift)
                .unwrap_or_else(|| "-".to_string());
            let fields = [
                candidate.row.parent_sg.to_string(),
                candidate.source_sg.to_string(),
                candidate.source_label.to_string(),
                candidate.row.condensing.to_string(),
                candidate.row.direction_index.to_string(),
                candidate.row.uni.to_string(),
                candidate.unitary_sg.to_string(),
                candidate.corep.label.clone(),
                candidate.h_sources[0].0.to_string(),
                candidate.h_sources[0].1.to_string(),
                candidate.corep.sources[0].bc.clone().unwrap_or_default(),
                candidate.h_sources[1].0.to_string(),
                candidate.h_sources[1].1.to_string(),
                candidate.corep.sources[1].bc.clone().unwrap_or_default(),
                sign_summary(&candidate.corep.sources[0]),
                sign_summary(&candidate.corep.sources[1]),
                format!("{:?}", candidate.corep.corep_type),
                candidate.corep.dim.to_string(),
                evidence.stage.clone(),
                evidence.reason.clone(),
                evidence.embedding_kind.clone(),
                map,
                shift,
                evidence.eta_plus.clone(),
                evidence.eta_minus.clone(),
                evidence.restriction.clone(),
                evidence.pairing.clone(),
                evidence.exchange_plus.clone(),
                evidence.exchange_minus.clone(),
            ];
            out.push_str(
                &fields
                    .iter()
                    .map(|field| sanitize(field))
                    .collect::<Vec<_>>()
                    .join("\t"),
            );
            out.push('\n');
        }
        out
    }
}

fn scan() -> Report {
    let fixed_sources = match fixed_sources() {
        Ok(sources) => sources,
        Err(error) => {
            return Report {
                fixed_sources: Vec::new(),
                counters: Counters::default(),
                candidates: Vec::new(),
                invariant_violations: vec![format!("fixed source validation: {error}")],
                stopped_at_metadata: true,
                pair_unresolved_details: Vec::new(),
                layer4_executed: false,
                witnesses: 0,
                candidate_incomplete: false,
            };
        }
    };

    let mut counters = Counters::default();
    let mut rows: Vec<Row> = Vec::new();
    for &(sg, _) in &FIXED_SOURCES {
        let (parent_rows, public_count, spinor_records) = rows_of_parent(sg);
        counters.public_rows += public_count;
        counters.spinor_records_skipped += spinor_records;
        counters.per_parent.entry(sg).or_default().rows = parent_rows.len();
        counters.rows += parent_rows.len();
        rows.extend(parent_rows);
    }

    let unis: BTreeSet<usize> = rows.iter().map(|row| row.uni).collect();
    counters.distinct_unis = unis.len();
    for (&sg, parent) in counters.per_parent.iter_mut() {
        parent.distinct_unis = rows
            .iter()
            .filter(|row| row.parent_sg == sg)
            .map(|row| row.uni)
            .collect::<BTreeSet<_>>()
            .len();
    }

    let mut gamma: BTreeMap<usize, UniGamma> = BTreeMap::new();
    for &uni in &unis {
        gamma.insert(uni, gamma_of(uni));
    }

    // Per-UNI layer-2 and layer-3 classification.  These are the **unique**
    // readings: one UNI contributes one set of Γ coreps no matter how many rows
    // point at it.
    let mut uni_screens: BTreeMap<usize, TypeCScreen> = BTreeMap::new();
    let mut pair_unresolved_details: Vec<String> = Vec::new();
    for (&uni, reading) in &gamma {
        match reading {
            UniGamma::SummaryFailed(detail) => {
                counters.summary_failed += 1;
                if counters.first_summary_error.is_none() {
                    counters.first_summary_error = Some(detail.clone());
                }
            }
            UniGamma::GammaMissing => counters.gamma_missing += 1,
            UniGamma::Resolved {
                unitary_sg: _,
                coreps,
                unresolved_gamma,
                unresolved_other,
                gamma_columns,
            } => {
                counters.unresolved_gamma += unresolved_gamma;
                counters.unresolved_other += unresolved_other;
                counters.gamma_columns += gamma_columns;
                let mut screen = TypeCScreen::default();
                for corep in coreps {
                    counters.gamma_coreps += 1;
                    match corep.corep_type {
                        CorepType::A => counters.type_a += 1,
                        CorepType::B => counters.type_b += 1,
                        CorepType::C => counters.type_c += 1,
                    }
                    if corep.corep_type == CorepType::C && corep.dim == 2 {
                        counters.type_c_dim2_unique += 1;
                        match screen_pair(corep.corep_type, corep.dim, &corep.sources) {
                            PairOutcome::Candidate { first, second } => {
                                counters.pair_candidates_unique += 1;
                                let _ = (first, second);
                                screen.candidates.push(corep.clone());
                            }
                            PairOutcome::Excluded(reason) => {
                                *counters
                                    .pair_excluded
                                    .entry(format!("{reason:?}"))
                                    .or_insert(0) += 1;
                            }
                            PairOutcome::Unresolved(reason) => {
                                *counters.pair_unresolved.entry(reason.name()).or_insert(0) += 1;
                                pair_unresolved_details.push(format!(
                                    "uni={uni} corep={} sources=[{}] problem={}",
                                    corep.label,
                                    corep
                                        .sources
                                        .iter()
                                        .map(|source| format!(
                                            "{}:{} dim={} spinor={} gamma_labels=[{}]",
                                            source.sg,
                                            source.label,
                                            source.dimension,
                                            source.spinor,
                                            query::irreps_of(source.sg)
                                                .iter()
                                                .filter(|ir| ir.k_label() == "GM")
                                                .map(|ir| ir.ml)
                                                .collect::<Vec<_>>()
                                                .join(",")
                                        ))
                                        .collect::<Vec<_>>()
                                        .join(" | "),
                                    reason.name()
                                ));
                            }
                        }
                    }
                }
                uni_screens.insert(uni, screen);
            }
        }
    }

    // Row expansion: every row sees the UNI's candidates.  A UNI shared by
    // several rows is several records; the row identity is not collapsed.
    let mut candidates = Vec::new();
    let mut non_candidate_rows = 0usize;
    for row in &rows {
        let screen = uni_screens.get(&row.uni).cloned().unwrap_or_default();
        let row_type_c = match gamma.get(&row.uni) {
            Some(UniGamma::Resolved { coreps, .. }) => coreps
                .iter()
                .filter(|corep| corep.corep_type == CorepType::C && corep.dim == 2)
                .count(),
            _ => 0,
        };
        counters.type_c_dim2_rows += row_type_c;
        counters.pair_candidates_rows += screen.candidates.len();
        counters.candidate_rows += screen.candidates.len();
        non_candidate_rows += row_type_c - screen.candidates.len();
        if let Some(parent) = counters.per_parent.get_mut(&row.parent_sg) {
            parent.type_c_dim2_rows += row_type_c;
            parent.candidate_rows += screen.candidates.len();
        }
        let source = fixed_sources
            .iter()
            .find(|source| source.sg == row.parent_sg)
            .expect("the row parent has a fixed source");
        let unitary_sg = match gamma.get(&row.uni) {
            Some(UniGamma::Resolved { unitary_sg, .. }) => *unitary_sg,
            _ => 0,
        };
        for corep in screen.candidates {
            let pair = [
                (corep.sources[0].sg, corep.sources[0].label),
                (corep.sources[1].sg, corep.sources[1].label),
            ];
            candidates.push(MetadataCandidate {
                row: *row,
                source_sg: source.sg,
                source_label: source.label,
                unitary_sg,
                corep,
                h_sources: pair,
                evidence: CandidateEvidence::default(),
            });
        }
    }

    let mut invariant_violations = Vec::new();
    if counters.rows != counters.public_rows {
        invariant_violations.push(format!(
            "rows {} != public row query {} (spinor records are a named skip, not rows)",
            counters.rows, counters.public_rows
        ));
    }
    let parent_row_sum: usize = counters.per_parent.values().map(|parent| parent.rows).sum();
    if parent_row_sum != counters.rows {
        invariant_violations.push(format!(
            "per-parent rows {parent_row_sum} != total {}",
            counters.rows
        ));
    }
    let resolved_unis = gamma
        .values()
        .filter(|reading| matches!(reading, UniGamma::Resolved { .. }))
        .count();
    if counters.distinct_unis != counters.summary_failed + counters.gamma_missing + resolved_unis {
        invariant_violations.push("distinct UNIs do not partition into readings".to_string());
    }
    let (excluded, unresolved): (usize, usize) = (
        counters.pair_excluded.values().sum(),
        counters.pair_unresolved.values().sum(),
    );
    if counters.type_c_dim2_unique != counters.pair_candidates_unique + excluded + unresolved {
        invariant_violations.push(format!(
            "unique type-C dim2 coreps {} != candidates {} + excluded {} + unresolved {}",
            counters.type_c_dim2_unique, counters.pair_candidates_unique, excluded, unresolved
        ));
    }
    if counters.type_c_dim2_rows != counters.pair_candidates_rows + non_candidate_rows {
        invariant_violations.push(format!(
            "row type-C dim2 exposures {} != candidate rows {} + excluded/unresolved rows {}",
            counters.type_c_dim2_rows, counters.pair_candidates_rows, non_candidate_rows
        ));
    }
    if candidates.len() != counters.candidate_rows {
        invariant_violations.push(format!(
            "candidate list {} != candidate rows {}",
            candidates.len(),
            counters.candidate_rows
        ));
    }

    let layer4_executed = !candidates.is_empty();
    let mut witnesses = 0usize;
    let mut candidate_incomplete = false;
    if layer4_executed {
        for candidate in &mut candidates {
            let evidence = verify_candidate(candidate);
            if evidence.stage == "witness" {
                witnesses += 1;
            }
            if evidence.incomplete {
                candidate_incomplete = true;
            }
            candidate.evidence = evidence;
        }
    }
    let stopped_at_metadata = counters.candidate_rows == 0;
    Report {
        fixed_sources,
        counters,
        candidates,
        invariant_violations,
        stopped_at_metadata,
        pair_unresolved_details,
        layer4_executed,
        witnesses,
        candidate_incomplete,
    }
}

#[derive(Clone, Default)]
struct TypeCScreen {
    candidates: Vec<GammaCorep>,
}

// ---------------------------------------------------------------------------
// Layer 4: explicit embedding, source model and a-twisted pairing
// ---------------------------------------------------------------------------

const RECONSTRUCTION_TOLERANCE: f64 = 1e-7;

/// A checked statement about one layer of one candidate.
enum CheckFailure {
    /// The candidate is incompatible under the fixed embedding.
    Definite(String),
    /// A prerequisite (operation binding, character column, separating class)
    /// is missing, so the layer is not established either way.
    Incomplete(String),
}

/// The fixed first embedding of the R24 candidate family, following the R28
/// policy: `RequireExact`, first realising entry under a source-independent
/// total key.
enum FixedEmbedding {
    Fixed(Box<ParentEmbedding>),
    NoRealising,
    CertifiedAbsence,
    UndecidedResidual,
}

type ExactMapKey = [(i128, i128); 9];
type ExactShiftKey = [(i128, i128); 3];

fn embedding_key(entry: &ParentEmbedding) -> (bool, ExactMapKey, ExactShiftKey) {
    let mut map = [(0i128, 1i128); 9];
    for row in 0..3 {
        for column in 0..3 {
            let value = entry.map.row(row)[column];
            map[row * 3 + column] = (value.numerator(), value.denominator());
        }
    }
    let mut shift = [(0i128, 1i128); 3];
    for (axis, value) in entry.shift.iter().enumerate() {
        shift[axis] = (value.numerator(), value.denominator());
    }
    (entry.map != Mat3R::identity(), map, shift)
}

fn fixed_embedding(sg: u8, row: &Row) -> Result<FixedEmbedding, String> {
    let rows = isotropy::magnetic_isotropy_subgroups(sg, row.condensing, LabelConvention::Cdml)
        .map_err(|error| format!("{error:?}"))?;
    let Some(subgroup) = rows.get(row.direction_index) else {
        return Err(format!(
            "direction {} is outside the condensing row list",
            row.direction_index
        ));
    };
    let geometry = geometry_of(sg, &subgroup.record).map_err(|error| format!("{error:?}"))?;
    let set = magnetic_operations(row.uni).map_err(|error| format!("{error:?}"))?;
    let verdict =
        parent_embedding_verdict(&geometry, &set).map_err(|error| format!("{error:?}"))?;
    let mut realising: Vec<&ParentEmbedding> = verdict
        .embeddings
        .iter()
        .filter(|entry| entry.realises_record_lattice)
        .collect();
    if realising.is_empty() {
        return Ok(if !verdict.embeddings.is_empty() {
            FixedEmbedding::NoRealising
        } else if verdict.absence_is_certified {
            FixedEmbedding::CertifiedAbsence
        } else {
            FixedEmbedding::UndecidedResidual
        });
    }
    realising.sort_by_key(|entry| embedding_key(entry));
    Ok(FixedEmbedding::Fixed(Box::new(realising[0].clone())))
}

fn verify_model(
    candidate: &MetadataCandidate,
    map: Mat3R,
    shift: [Rat; 3],
    eta: RealSign,
) -> Result<MatrixGammaModel, String> {
    let request = MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: candidate.row.parent_sg,
            label: candidate.source_label,
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: candidate.row.parent_sg,
            condensing_label: candidate.row.condensing,
            convention: LabelConvention::Cdml,
            direction_index: candidate.row.direction_index,
        },
        Some(eta),
    )
    .with_embedding(EmbeddingEvidence {
        map,
        shift,
        lattice_policy: LatticePolicy::RequireExact,
    });
    verify_matrix_gamma_input(&request).map_err(|error| format!("{error:?}"))
}

fn verify_candidate(candidate: &MetadataCandidate) -> CandidateEvidence {
    let mut evidence = CandidateEvidence {
        stage: "unexecuted".to_string(),
        eta_plus: "not_run".to_string(),
        eta_minus: "not_run".to_string(),
        restriction: "not_run".to_string(),
        pairing: "not_run".to_string(),
        exchange_plus: "not_run".to_string(),
        exchange_minus: "not_run".to_string(),
        ..CandidateEvidence::default()
    };
    let embedding = match fixed_embedding(candidate.row.parent_sg, &candidate.row) {
        Ok(embedding) => embedding,
        Err(detail) => {
            evidence.reason = format!("embedding_error {detail}");
            evidence.incomplete = true;
            return evidence;
        }
    };
    let ParentEmbedding {
        map, shift, kind, ..
    } = match embedding {
        FixedEmbedding::Fixed(entry) => *entry,
        FixedEmbedding::NoRealising => {
            evidence.reason = "no_realising_candidate".to_string();
            evidence.incomplete = true;
            return evidence;
        }
        FixedEmbedding::CertifiedAbsence => {
            evidence.reason = "certified_absence_in_the_R24_family".to_string();
            evidence.definite = true;
            return evidence;
        }
        FixedEmbedding::UndecidedResidual => {
            evidence.reason = "undecided_residual".to_string();
            evidence.incomplete = true;
            return evidence;
        }
    };
    evidence.stage = "not_established".to_string();
    evidence.map = Some(map);
    evidence.shift = Some(shift);
    evidence.embedding_kind = format!("{kind:?}");

    let plus = verify_model(candidate, map, shift, RealSign::Plus);
    let minus = verify_model(candidate, map, shift, RealSign::Minus);
    evidence.eta_plus = match &plus {
        Ok(_) => "ok".to_string(),
        Err(error) => format!("error {error}"),
    };
    evidence.eta_minus = match &minus {
        Ok(_) => "ok".to_string(),
        Err(error) => format!("error {error}"),
    };
    let (plus, minus) = match (plus, minus) {
        (Ok(plus), Ok(minus)) => (plus, minus),
        _ => {
            evidence.reason = "source_model_not_verified".to_string();
            evidence.incomplete = true;
            return evidence;
        }
    };

    let set = match magnetic_operations(candidate.row.uni) {
        Ok(set) => set,
        Err(error) => {
            evidence.reason = format!("operations_error {error:?}");
            evidence.incomplete = true;
            return evidence;
        }
    };

    match check_source_restriction(&plus, &candidate.corep, &set) {
        Ok(classes) => evidence.restriction = format!("ok classes={classes}"),
        Err(CheckFailure::Definite(detail)) => {
            evidence.restriction = detail;
            evidence.reason = "source_restriction_mismatch".to_string();
            evidence.definite = true;
            return evidence;
        }
        Err(CheckFailure::Incomplete(detail)) => {
            evidence.restriction = detail;
            evidence.reason = "source_restriction_unresolved".to_string();
            evidence.incomplete = true;
            return evidence;
        }
    }

    match check_a_twisted_pairing(&candidate.corep, &set) {
        Ok((classes, discriminating)) => {
            evidence.pairing = format!("ok classes={classes} discriminating={discriminating}");
        }
        Err(CheckFailure::Definite(detail)) => {
            evidence.pairing = detail;
            evidence.reason = "a_twisted_pairing_mismatch".to_string();
            evidence.definite = true;
            return evidence;
        }
        Err(CheckFailure::Incomplete(detail)) => {
            evidence.pairing = detail;
            evidence.reason = "a_twisted_pairing_unresolved".to_string();
            evidence.incomplete = true;
            return evidence;
        }
    }

    match check_antiunitary_exchange(&plus, &candidate.corep) {
        Ok((antiunitary, unitary)) => {
            evidence.exchange_plus = format!("ok anti={antiunitary} unitary={unitary}");
        }
        Err(CheckFailure::Definite(detail)) => {
            evidence.exchange_plus = detail;
            evidence.reason = "antiunitary_exchange_mismatch_eta_plus".to_string();
            evidence.definite = true;
            return evidence;
        }
        Err(CheckFailure::Incomplete(detail)) => {
            evidence.exchange_plus = detail;
            evidence.reason = "antiunitary_exchange_unresolved_eta_plus".to_string();
            evidence.incomplete = true;
            return evidence;
        }
    }
    match check_antiunitary_exchange(&minus, &candidate.corep) {
        Ok((antiunitary, unitary)) => {
            evidence.exchange_minus = format!("ok anti={antiunitary} unitary={unitary}");
        }
        Err(CheckFailure::Definite(detail)) => {
            evidence.exchange_minus = detail;
            evidence.reason = "antiunitary_exchange_mismatch_eta_minus".to_string();
            evidence.definite = true;
            return evidence;
        }
        Err(CheckFailure::Incomplete(detail)) => {
            evidence.exchange_minus = detail;
            evidence.reason = "antiunitary_exchange_unresolved_eta_minus".to_string();
            evidence.incomplete = true;
            return evidence;
        }
    }

    evidence.stage = "witness".to_string();
    evidence.reason = "ok".to_string();
    evidence.definite = true;
    evidence
}

type ClassKey = ([i32; 9], [(i128, i128); 3], bool);

fn class_key(operation: &MagneticOperation, lattice: &[[Rat; 3]]) -> ClassKey {
    let translation = canonical_translation(operation.translation, lattice);
    (
        [
            operation.rotation[0][0],
            operation.rotation[0][1],
            operation.rotation[0][2],
            operation.rotation[1][0],
            operation.rotation[1][1],
            operation.rotation[1][2],
            operation.rotation[2][0],
            operation.rotation[2][1],
            operation.rotation[2][2],
        ],
        [
            (translation[0].numerator(), translation[0].denominator()),
            (translation[1].numerator(), translation[1].denominator()),
            (translation[2].numerator(), translation[2].denominator()),
        ],
        operation.time_reversal,
    )
}

fn check_source_restriction(
    model: &MatrixGammaModel,
    corep: &GammaCorep,
    set: &MagneticOperationSet,
) -> Result<usize, CheckFailure> {
    if model.source_operators().len() != model.classes().len() {
        return Err(CheckFailure::Incomplete(
            "the model action table does not have one entry per class".to_string(),
        ));
    }
    if corep.sources.len() != 2 {
        return Err(CheckFailure::Incomplete(
            "the paired corep does not expose exactly two H sources".to_string(),
        ));
    }
    let first = match &corep.sources[0].character {
        Ok(row) => row,
        Err(problem) => {
            return Err(CheckFailure::Incomplete(format!(
                "H source {}:{} has no real row: {}",
                corep.sources[0].sg,
                corep.sources[0].label,
                problem.name()
            )));
        }
    };
    let second = match &corep.sources[1].character {
        Ok(row) => row,
        Err(problem) => {
            return Err(CheckFailure::Incomplete(format!(
                "H source {}:{} has no real row: {}",
                corep.sources[1].sg,
                corep.sources[1].label,
                problem.name()
            )));
        }
    };
    let lattice = magnetic_embedding::translation_lattice(set);
    let mut class_of_raw: BTreeMap<ClassKey, usize> = BTreeMap::new();
    for (index, class) in model.classes().iter().enumerate() {
        class_of_raw.insert(class_key(class, &lattice), index);
    }
    let mut expected: BTreeMap<usize, f64> = BTreeMap::new();
    let mut expected_from_h: BTreeMap<usize, f64> = BTreeMap::new();
    for (column_index, column) in corep.columns.iter().enumerate() {
        if column.time_reversal {
            continue;
        }
        let Some(Some(character)) = corep.characters.get(column_index) else {
            return Err(CheckFailure::Incomplete(format!(
                "unitary column {column_index} has no corep character"
            )));
        };
        if character.im.abs() > RECONSTRUCTION_TOLERANCE {
            return Err(CheckFailure::Definite(format!(
                "unitary corep character {column_index} is not real: {character:?}"
            )));
        }
        let Some(raw) = set.operations.get(column.raw_index) else {
            return Err(CheckFailure::Incomplete(format!(
                "raw magnetic operation {} is out of range",
                column.raw_index
            )));
        };
        let Some(&class_index) = class_of_raw.get(&class_key(raw, &lattice)) else {
            return Err(CheckFailure::Incomplete(format!(
                "raw operation {} is not a class of the model",
                column.raw_index
            )));
        };
        let h_key = column_key(column)?;
        let (Some(&first_h), Some(&second_h)) =
            (first.values.get(&h_key), second.values.get(&h_key))
        else {
            return Err(CheckFailure::Incomplete(format!(
                "unitary column {column_index} has no H operation key in the typed row"
            )));
        };
        let from_h = first_h as f64 + second_h as f64;
        if (from_h - character.re).abs() > RECONSTRUCTION_TOLERANCE {
            return Err(CheckFailure::Incomplete(format!(
                "class {class_index}: summary chi1+chi2 {} != typed H row sum {from_h}",
                character.re
            )));
        }
        if let Some(previous) = expected.insert(class_index, character.re)
            && (previous - character.re).abs() > RECONSTRUCTION_TOLERANCE
        {
            return Err(CheckFailure::Definite(format!(
                "class {class_index} has two different unitary characters"
            )));
        }
        if let Some(previous) = expected_from_h.insert(class_index, from_h)
            && (previous - from_h).abs() > RECONSTRUCTION_TOLERANCE
        {
            return Err(CheckFailure::Definite(format!(
                "class {class_index} has two different typed H row sums"
            )));
        }
    }
    let mut checked = 0usize;
    for (index, class) in model.classes().iter().enumerate() {
        if class.time_reversal {
            continue;
        }
        let Some(&value) = expected.get(&index) else {
            return Err(CheckFailure::Incomplete(format!(
                "unitary class {index} is not covered by the summary columns"
            )));
        };
        let Some(&from_h) = expected_from_h.get(&index) else {
            return Err(CheckFailure::Incomplete(format!(
                "unitary class {index} has no typed H row sum"
            )));
        };
        let entries = model.source_operators()[index].entries();
        let trace = entries[0] + entries[3];
        if (trace - value).abs() > RECONSTRUCTION_TOLERANCE
            || (trace - from_h).abs() > RECONSTRUCTION_TOLERANCE
        {
            return Err(CheckFailure::Definite(format!(
                "class {index}: trace(source action) {trace} != chi1+chi2 {value}                  (typed H rows {from_h})"
            )));
        }
        checked += 1;
    }
    Ok(checked)
}

fn column_key(column: &GammaColumn) -> Result<([i32; 9], [i64; 3]), CheckFailure> {
    operation_key(&SeitzOperation {
        rotation: column.rotation,
        translation: column.translation,
    })
    .map_err(|problem| CheckFailure::Incomplete(format!("H operation key: {}", problem.name())))
}

fn check_a_twisted_pairing(
    corep: &GammaCorep,
    set: &MagneticOperationSet,
) -> Result<(usize, usize), CheckFailure> {
    if corep.sources.len() != 2 {
        return Err(CheckFailure::Incomplete(
            "the paired corep does not expose exactly two H sources".to_string(),
        ));
    }
    let first = match &corep.sources[0].character {
        Ok(row) => row,
        Err(problem) => {
            return Err(CheckFailure::Incomplete(format!(
                "H source {}:{} has no real row: {}",
                corep.sources[0].sg,
                corep.sources[0].label,
                problem.name()
            )));
        }
    };
    let second = match &corep.sources[1].character {
        Ok(row) => row,
        Err(problem) => {
            return Err(CheckFailure::Incomplete(format!(
                "H source {}:{} has no real row: {}",
                corep.sources[1].sg,
                corep.sources[1].label,
                problem.name()
            )));
        }
    };
    let lattice = magnetic_embedding::translation_lattice(set);
    let mut column_of_class: BTreeMap<ClassKey, usize> = BTreeMap::new();
    for (index, column) in corep.columns.iter().enumerate() {
        let Some(raw) = set.operations.get(column.raw_index) else {
            return Err(CheckFailure::Incomplete(format!(
                "raw magnetic operation {} is out of range",
                column.raw_index
            )));
        };
        column_of_class
            .entry(class_key(raw, &lattice))
            .or_insert(index);
    }
    let Some(a0_column) = corep.columns.iter().find(|column| column.time_reversal) else {
        return Err(CheckFailure::Incomplete(
            "the type-C corep has no antiunitary column".to_string(),
        ));
    };
    let Some(a0) = set.operations.get(a0_column.raw_index) else {
        return Err(CheckFailure::Incomplete(
            "the antiunitary representative is out of range".to_string(),
        ));
    };
    let Some(a0_inverse) = inverse_operation(a0, set) else {
        return Err(CheckFailure::Incomplete(
            "the antiunitary representative has no inverse in the operation list".to_string(),
        ));
    };

    let mut checked = 0usize;
    let mut discriminating = 0usize;
    for column in corep.columns.iter().filter(|column| !column.time_reversal) {
        let Some(h_raw) = set.operations.get(column.raw_index) else {
            return Err(CheckFailure::Incomplete(format!(
                "raw magnetic operation {} is out of range",
                column.raw_index
            )));
        };
        let h_conjugated = compose(h_raw, a0)
            .and_then(|right| compose(&a0_inverse, &right))
            .map_err(|error| {
                CheckFailure::Incomplete(format!("a0 conjugation failed: {error:?}"))
            })?;
        let Some(&partner_column) = column_of_class.get(&class_key(&h_conjugated, &lattice)) else {
            return Err(CheckFailure::Incomplete(
                "the a0-conjugated unitary class has no summary column".to_string(),
            ));
        };
        let h_key = column_key(column)?;
        let partner_key = column_key(&corep.columns[partner_column])?;
        let (Some(&first_h), Some(&second_h)) =
            (first.values.get(&h_key), second.values.get(&h_key))
        else {
            return Err(CheckFailure::Incomplete(
                "a unitary column has no H operation key in the typed row".to_string(),
            ));
        };
        let (Some(&first_partner), Some(&second_partner)) = (
            first.values.get(&partner_key),
            second.values.get(&partner_key),
        ) else {
            return Err(CheckFailure::Incomplete(
                "the a0-conjugated column has no H operation key in the typed row".to_string(),
            ));
        };
        checked += 1;
        if first_h != second_h {
            discriminating += 1;
        }
        if first_partner != second_h || second_partner != first_h {
            return Err(CheckFailure::Definite(format!(
                "column {}: chi1(a0^-1 h a0)={first_partner}, chi2(h)={second_h}; \
                 chi2(a0^-1 h a0)={second_partner}, chi1(h)={first_h}",
                column.raw_index
            )));
        }
    }
    Ok((checked, discriminating))
}

fn inverse_operation(
    operation: &MagneticOperation,
    set: &MagneticOperationSet,
) -> Option<MagneticOperation> {
    let lattice = magnetic_embedding::translation_lattice(set);
    for candidate in &set.operations {
        let Ok(left) = compose(operation, candidate) else {
            continue;
        };
        let Ok(right) = compose(candidate, operation) else {
            continue;
        };
        if identity_class(&left, &lattice) && identity_class(&right, &lattice) {
            return Some(*candidate);
        }
    }
    None
}

fn identity_class(operation: &MagneticOperation, lattice: &[[Rat; 3]]) -> bool {
    !operation.time_reversal
        && operation.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        && canonical_translation(operation.translation, lattice)
            .iter()
            .all(|value| value.is_zero())
}

/// The common real eigenbasis of the unitary subgroup action: one vector for
/// `+1` and one for `-1` on a class where the two H characters differ.
type ExchangeBasis = (f64, f64, f64, f64);

fn separating_basis(
    model: &MatrixGammaModel,
    corep: &GammaCorep,
) -> Result<ExchangeBasis, CheckFailure> {
    let mut saw_unitary = false;
    for (index, class) in model.classes().iter().enumerate() {
        if class.time_reversal {
            continue;
        }
        saw_unitary = true;
        let entries = model.source_operators()[index].entries();
        let (a, b, c, d) = (entries[0], entries[1], entries[2], entries[3]);
        if (a + d).abs() > RECONSTRUCTION_TOLERANCE
            || (a * d - b * c + 1.0).abs() > RECONSTRUCTION_TOLERANCE
        {
            continue;
        }
        let (mut x, mut y) = if b.abs() > 1e-12 {
            (b, 1.0 - a)
        } else if a > 0.0 {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        };
        let norm = (x * x + y * y).sqrt();
        if !norm.is_finite() || norm < 1e-12 {
            continue;
        }
        x /= norm;
        y /= norm;
        let (mx, my) = (-y, x);
        return Ok((x, y, mx, my));
    }
    let _ = corep;
    Err(CheckFailure::Incomplete(format!(
        "no separating unitary class among {saw_unitary} unitary classes"
    )))
}

fn conjugate_matrix(entries: [f64; 4], basis: ExchangeBasis) -> [f64; 4] {
    let (px, py, mx, my) = basis;
    let (a, b, c, d) = (entries[0], entries[1], entries[2], entries[3]);
    [
        px * px * a + px * py * b + py * px * c + py * py * d,
        px * mx * a + px * my * b + py * mx * c + py * my * d,
        mx * px * a + mx * py * b + my * px * c + my * py * d,
        mx * mx * a + mx * my * b + my * mx * c + my * my * d,
    ]
}

fn check_antiunitary_exchange(
    model: &MatrixGammaModel,
    corep: &GammaCorep,
) -> Result<(usize, usize), CheckFailure> {
    let basis = separating_basis(model, corep)?;
    let mut antiunitary = 0usize;
    let mut unitary = 0usize;
    for (index, class) in model.classes().iter().enumerate() {
        let matrix = conjugate_matrix(model.source_operators()[index].entries(), basis);
        if class.time_reversal {
            if matrix[0].abs() > RECONSTRUCTION_TOLERANCE
                || matrix[3].abs() > RECONSTRUCTION_TOLERANCE
            {
                return Err(CheckFailure::Definite(format!(
                    "antiunitary class {index} keeps a diagonal entry in the H eigenbasis: {matrix:?}"
                )));
            }
            antiunitary += 1;
            if matrix[1].abs() < RECONSTRUCTION_TOLERANCE
                || matrix[2].abs() < RECONSTRUCTION_TOLERANCE
            {
                return Err(CheckFailure::Definite(format!(
                    "antiunitary class {index} does not exchange the two H constituents: {matrix:?}"
                )));
            }
        } else {
            if matrix[1].abs() > RECONSTRUCTION_TOLERANCE
                || matrix[2].abs() > RECONSTRUCTION_TOLERANCE
            {
                return Err(CheckFailure::Definite(format!(
                    "unitary class {index} is not diagonal in the H eigenbasis: {matrix:?}"
                )));
            }
            unitary += 1;
        }
    }
    Ok((antiunitary, unitary))
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

/// The checked real `±1` content of one H source row, for the witness line.
fn sign_summary(source: &SourceFacts) -> String {
    match &source.character {
        Ok(row) => {
            let plus = row.values.values().filter(|value| **value == 1).count();
            let minus = row.values.values().filter(|value| **value == -1).count();
            format!("dim={} +{plus}/-{minus}", source.dimension)
        }
        Err(problem) => format!("dim={} unavailable:{}", source.dimension, problem.name()),
    }
}

fn print_report(report: &Report) {
    println!("R31 finite TypeC Gamma feasibility gate");
    println!("baseline {BASELINE}");
    println!("scope fixed_parents=3 fixed_sources=3 production_api=none");
    for source in &report.fixed_sources {
        println!(
            "fixed_source sg={} label={} dim={} spinor={}",
            source.sg, source.label, source.selected_dimension, source.spinor
        );
    }
    let counters = &report.counters;
    for (sg, parent) in &counters.per_parent {
        println!(
            "rows sg={} rows={} distinct_unis={} type_c_dim2_row_exposure={} \
             metadata_candidates={}",
            sg, parent.rows, parent.distinct_unis, parent.type_c_dim2_rows, parent.candidate_rows
        );
    }
    println!(
        "layer1 rows={} public_rows={} spinor_records_skipped={} distinct_unis={}",
        counters.rows,
        counters.public_rows,
        counters.spinor_records_skipped,
        counters.distinct_unis
    );
    println!(
        "layer1 uni_summary failed={} gamma_missing={} unresolved_gamma={} \
         unresolved_other_k={} gamma_columns={}",
        counters.summary_failed,
        counters.gamma_missing,
        counters.unresolved_gamma,
        counters.unresolved_other,
        counters.gamma_columns
    );
    if let Some(detail) = &counters.first_summary_error {
        println!("layer1 first_summary_error {detail}");
    }
    println!(
        "layer2 gamma_coreps={} type_a={} type_b={} type_c={} \
         type_c_dim2_unique={} type_c_dim2_row_exposure={}",
        counters.gamma_coreps,
        counters.type_a,
        counters.type_b,
        counters.type_c,
        counters.type_c_dim2_unique,
        counters.type_c_dim2_rows
    );
    println!(
        "layer3 metadata_candidates_unique={} metadata_candidates_rows={} \
         excluded_unique={} unresolved_unique={}",
        counters.pair_candidates_unique,
        counters.pair_candidates_rows,
        counters.pair_excluded.values().sum::<usize>(),
        counters.pair_unresolved.values().sum::<usize>()
    );
    for (reason, count) in &counters.pair_excluded {
        println!("exclude {reason:?} count={count}");
    }
    for (reason, count) in &counters.pair_unresolved {
        println!("unresolved {reason} count={count}");
    }
    for detail in &report.pair_unresolved_details {
        println!("unresolved_detail {detail}");
    }
    for violation in &report.invariant_violations {
        println!("INVARIANT VIOLATION {violation}");
    }
    if report.layer4_executed {
        println!(
            "layer4 executed=true candidates={} witnesses={} incomplete={}",
            report.candidates.len(),
            report.witnesses,
            report.candidate_incomplete
        );
    } else {
        println!("layer4 not_executed reason=zero_metadata_candidates");
    }
    println!(
        "status {}",
        if report.incomplete() {
            "incomplete"
        } else if report.stopped_at_metadata {
            "complete_metadata_zero_candidates"
        } else if report.witnesses == report.candidates.len() {
            "complete_all_candidates_are_witnesses"
        } else if report.witnesses == 0 {
            "complete_candidates_no_witness"
        } else {
            "complete_mixed_candidates"
        }
    );
    for candidate in &report.candidates {
        let evidence = &candidate.evidence;
        let corep_type = format!("{:?}", candidate.corep.corep_type);
        let map = evidence
            .map
            .map(format_map)
            .unwrap_or_else(|| "-".to_string());
        let shift = evidence
            .shift
            .map(format_shift)
            .unwrap_or_else(|| "-".to_string());
        println!(
            "candidate sg={} source_sg={} source={} condensing={} dir={} uni={} corep={} \
             unitary_sg={} corep_type={} corep_dim={} h1={}:{} ({}) h2={}:{} ({}) h1_signs={} h2_signs={} \
             stage={} reason={} kind={} map={} shift={} eta_plus={} eta_minus={} \
             restriction={} pairing={} exchange_plus={} exchange_minus={}",
            candidate.row.parent_sg,
            candidate.source_sg,
            candidate.source_label,
            candidate.row.condensing,
            candidate.row.direction_index,
            candidate.row.uni,
            candidate.corep.label,
            candidate.unitary_sg,
            corep_type,
            candidate.corep.dim,
            candidate.h_sources[0].0,
            candidate.h_sources[0].1,
            candidate.corep.sources[0].bc.as_deref().unwrap_or("-"),
            candidate.h_sources[1].0,
            candidate.h_sources[1].1,
            candidate.corep.sources[1].bc.as_deref().unwrap_or("-"),
            sign_summary(&candidate.corep.sources[0]),
            sign_summary(&candidate.corep.sources[1]),
            evidence.stage,
            evidence.reason,
            evidence.embedding_kind,
            map,
            shift,
            evidence.eta_plus,
            evidence.eta_minus,
            evidence.restriction,
            evidence.pairing,
            evidence.exchange_plus,
            evidence.exchange_minus,
        );
    }
}

fn format_map(map: Mat3R) -> String {
    let mut out = String::new();
    for row in 0..3 {
        if row > 0 {
            out.push(';');
        }
        for column in 0..3 {
            if column > 0 {
                out.push(',');
            }
            let value = map.row(row)[column];
            out.push_str(&format!("{}/{}", value.numerator(), value.denominator()));
        }
    }
    out
}

fn format_shift(shift: [Rat; 3]) -> String {
    shift
        .iter()
        .map(|value| format!("{}/{}", value.numerator(), value.denominator()))
        .collect::<Vec<_>>()
        .join(",")
}

fn sanitize(text: &str) -> String {
    text.replace(['\t', '\n', '\r'], " ")
}

// ---------------------------------------------------------------------------
// Self-check
// ---------------------------------------------------------------------------

fn synthetic_row(pairs: &[([i32; 9], [i64; 3], i8)]) -> RealOneDimRow {
    RealOneDimRow {
        values: pairs
            .iter()
            .map(|(rotation, translation, value)| ((*rotation, *translation), *value))
            .collect(),
    }
}

fn synthetic_source(
    label: &'static str,
    spinor: bool,
    dimension: u8,
    character: Result<RealOneDimRow, CharacterProblem>,
) -> SourceFacts {
    SourceFacts {
        sg: 1,
        label,
        bc: None,
        spinor,
        dimension,
        character,
    }
}

fn synthetic_source_from_values(
    label: &'static str,
    spinor: bool,
    dimension: u8,
    values: &[(Complex64, SeitzOperation)],
) -> SourceFacts {
    SourceFacts {
        sg: 1,
        label,
        bc: None,
        spinor,
        dimension,
        character: real_one_dim_row_from_pairs(values),
    }
}

/// The built-in checks.  The first one is the positive control: a type-C
/// dimension-2 pair of inequivalent real `±1` rows must pass the same screen
/// used by the scan.  The remaining cases are classified by name, so a
/// type-A corep, a complex/spinor/wrong-dimensional source, a missing
/// character row or a non-finite value can never silently become "zero
/// candidates".  The non-finite cases go through the real numeric
/// `grid_sign` -> row -> screen path, not through a pre-fabricated error.
fn self_check() -> Vec<String> {
    let identity = ([1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0]);
    let swap = ([-1, 0, 0, 0, -1, 0, 0, 0, 1], [0, 0, 0]);
    let plus = synthetic_row(&[(identity.0, identity.1, 1), (swap.0, swap.1, 1)]);
    let minus = synthetic_row(&[(identity.0, identity.1, 1), (swap.0, swap.1, -1)]);
    let real_pair = [
        synthetic_source("GM1", false, 1, Ok(plus.clone())),
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    let mut violations = Vec::new();

    match screen_pair(CorepType::C, 2, &real_pair) {
        PairOutcome::Candidate { first, second } if (first, second) == (0, 1) => {}
        other => violations.push(format!(
            "positive control did not become a candidate: {}",
            outcome_name(&other)
        )),
    }
    match screen_pair(CorepType::A, 1, &real_pair) {
        PairOutcome::Unresolved(PairUnresolved::NotTypeCDim2 { .. }) => {}
        other => violations.push(format!(
            "type-A negative control was misclassified: {}",
            outcome_name(&other)
        )),
    }
    let complex_pair = [
        synthetic_source("GM1", false, 1, Err(CharacterProblem::NonReal)),
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &complex_pair) {
        PairOutcome::Excluded(PairExclusion::NonReal(0)) => {}
        other => violations.push(format!(
            "complex source was not excluded as NonReal: {}",
            outcome_name(&other)
        )),
    }
    let identity_operation = SeitzOperation {
        rotation: identity.0,
        translation: [0.0, 0.0, 0.0],
    };
    let swap_operation = SeitzOperation {
        rotation: swap.0,
        translation: [0.0, 0.0, 0.0],
    };
    let non_finite_cases = [
        ("NaN", Complex64::new(f64::NAN, 0.0)),
        ("+Inf", Complex64::new(f64::INFINITY, 0.0)),
        ("-Inf", Complex64::new(f64::NEG_INFINITY, 0.0)),
        ("imagNaN", Complex64::new(1.0, f64::NAN)),
    ];
    for (name, bad) in non_finite_cases {
        let source = synthetic_source_from_values(
            "GM1",
            false,
            1,
            &[
                (Complex64::new(1.0, 0.0), identity_operation),
                (bad, swap_operation),
            ],
        );
        if !matches!(source.character, Err(CharacterProblem::NonFinite)) {
            violations.push(format!(
                "{name} did not classify as NonFinite: {:?}",
                source.character
            ));
            continue;
        }
        let pair = [source, synthetic_source("GM2", false, 1, Ok(minus.clone()))];
        match screen_pair(CorepType::C, 2, &pair) {
            PairOutcome::Unresolved(PairUnresolved::Character {
                index: 0,
                problem: CharacterProblem::NonFinite,
            }) => {}
            other => violations.push(format!(
                "{name} was not unresolved through the screen: {}",
                outcome_name(&other)
            )),
        }
    }
    let finite_complex = synthetic_source_from_values(
        "GM1",
        false,
        1,
        &[
            (Complex64::new(1.0, 0.0), identity_operation),
            (Complex64::new(1.0, 1.0), swap_operation),
        ],
    );
    let finite_complex_pair = [
        finite_complex,
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &finite_complex_pair) {
        PairOutcome::Excluded(PairExclusion::NonReal(0)) => {}
        other => violations.push(format!(
            "finite non-real value was not excluded as NonReal: {}",
            outcome_name(&other)
        )),
    }
    let finite_off_grid = synthetic_source_from_values(
        "GM1",
        false,
        1,
        &[
            (Complex64::new(1.0, 0.0), identity_operation),
            (Complex64::new(0.5, 0.0), swap_operation),
        ],
    );
    let finite_off_grid_pair = [
        finite_off_grid,
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &finite_off_grid_pair) {
        PairOutcome::Excluded(PairExclusion::NonReal(0)) => {}
        other => violations.push(format!(
            "finite 0.5 value was not excluded as OffGrid/NonReal: {}",
            outcome_name(&other)
        )),
    }
    let spinor_pair = [
        synthetic_source("GM1", true, 1, Ok(plus.clone())),
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &spinor_pair) {
        PairOutcome::Excluded(PairExclusion::Spinor(0)) => {}
        other => violations.push(format!(
            "spinor source was not excluded as Spinor: {}",
            outcome_name(&other)
        )),
    }
    let dimension_pair = [
        synthetic_source("GM1", false, 2, Ok(plus.clone())),
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &dimension_pair) {
        PairOutcome::Excluded(PairExclusion::Dimension {
            index: 0,
            dimension: 2,
        }) => {}
        other => violations.push(format!(
            "dimension-2 source was not excluded as Dimension: {}",
            outcome_name(&other)
        )),
    }
    let equivalent_pair = [
        synthetic_source("GM1", false, 1, Ok(plus.clone())),
        synthetic_source("GM2", false, 1, Ok(plus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &equivalent_pair) {
        PairOutcome::Unresolved(PairUnresolved::EquivalentPair) => {}
        other => violations.push(format!(
            "equivalent pair was not unresolved: {}",
            outcome_name(&other)
        )),
    }
    let missing_pair = [
        synthetic_source(
            "GM1",
            false,
            1,
            Err(CharacterProblem::Unavailable("MissingData")),
        ),
        synthetic_source("GM2", false, 1, Ok(minus.clone())),
    ];
    match screen_pair(CorepType::C, 2, &missing_pair) {
        PairOutcome::Unresolved(PairUnresolved::Character { index: 0, .. }) => {}
        other => violations.push(format!(
            "missing character row was not unresolved: {}",
            outcome_name(&other)
        )),
    }
    let single_source = [synthetic_source("GM1", false, 1, Ok(plus.clone()))];
    match screen_pair(CorepType::C, 2, &single_source) {
        PairOutcome::Unresolved(PairUnresolved::SourceCount(1)) => {}
        other => violations.push(format!(
            "single-source type C was not unresolved: {}",
            outcome_name(&other)
        )),
    }
    violations
}

fn outcome_name(outcome: &PairOutcome) -> String {
    match outcome {
        PairOutcome::Candidate { first, second } => {
            format!("Candidate({first},{second})")
        }
        PairOutcome::Excluded(reason) => format!("Excluded({reason:?})"),
        PairOutcome::Unresolved(reason) => format!("Unresolved({})", reason.name()),
    }
}
