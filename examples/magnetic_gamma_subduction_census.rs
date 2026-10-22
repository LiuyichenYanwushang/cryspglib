//! R28 restricted scalar-Γ magnetic-subduction **coverage diagnostic**.
//!
//! This example is a census, not a new algorithm: it points the accepted R25
//! input contract ([`verify_scalar_gamma_input`]) and the accepted R26/R27
//! restricted subduction ([`subduce_scalar_gamma`]) at **every** magnetic
//! isotropy record of every one of the 230 parent space groups and at every
//! Γ source of that parent, and reports what the real entry points answer.
//!
//! # The measured population
//!
//! * A **record** is one row of the generated magnetic isotropy tables (the same
//!   16,721 rows `magnetic_embedding_census` walks).  Its identity in this
//!   census is `(sg, condensing CDML label, per-irrep direction index, UNI)`;
//!   two rows that share a UNI are two records and are never merged.
//! * A **source** is a Γ irrep of the parent.  The domain screen is stated
//!   explicitly: non-spinor, dimension 1, and real `±1` character.  Every Γ
//!   irrep excluded by the screen is counted under its own reason; a spinor, a
//!   higher-dimensional or a complex source is never silently dropped and never
//!   probed as if it were in the family.
//! * A **(record, source) pair** is probed only with **one fixed embedding**
//!   chosen from the R24 search (`parent_embedding_verdict`) with the exact
//!   lattice policy `RequireExact`: the entries whose `realises_record_lattice`
//!   flag is set are sorted by a total, source-independent key and the first is
//!   taken.  This **reuses the R24 candidate family as it stands** and neither
//!   adds nor widens any setting search.  The choice is a function of the record
//!   alone, so it can never be picked to make a subduction succeed.  A record
//!   with no realising entry is
//!   reported as such and not probed; whether that absence is *certified*, a
//!   plain "no realising candidate among the contained ones", or an
//!   *undecided arithmetic residual* is kept apart, because an undecided answer
//!   is not a proof of absence.  **Scope fence**: the R24 verdict is a statement
//!   about its candidate family `{S·B}` (48 signed permutations times the five
//!   tabulated bases), so `certified_absence` means "no member of that family
//!   embeds the record", never "no setting of the record embeds it".
//! * Every probed pair is called **twice**, with `η = +1` and `η = -1`, through
//!   the real contract and the real subduction, and every probe lands in exactly
//!   one stage: `verify_input`, `subduce`, `success_check` or `success`.
//!   Successful results are re-checked class by class: the dimension identity
//!   `dim(source) = multiplicity·dim(target)`, a finite non-zero unit `U`, and
//!   the reconstruction `A(g)·(U z) = U·(D(g) z)` on complex vectors (the source
//!   action of the model on the left, the reconstructed target operator on the
//!   right, conjugated when `g` is antiunitary), so a "success" printed here is a
//!   reconstructed one, not an unexamined `Ok`.
//!
//! The census therefore measures a **conditional** denominator: the probes are
//! "over the fixed first realising embedding of the R24 family, under
//! `RequireExact`", and its success count is not a success rate over all
//! embeddings or settings of a record, which this round does not enumerate.
//!
//! # Usage
//!
//! ```text
//! magnetic_gamma_subduction_census --inventory        # denominator only
//! magnetic_gamma_subduction_census --smoke            # SGs 1-12 plus 221
//! magnetic_gamma_subduction_census --sgs 1-12,221     # explicit scope
//! magnetic_gamma_subduction_census --tsv target/logs/r28-census-worker/census.tsv
//! magnetic_gamma_subduction_census --gate             # exit 1 unless every invariant holds
//! ```
//!
//! `RAYON_NUM_THREADS` limits the parallelism (one space group per task).  The
//! TSV is written in space-group order, so its bytes do not depend on the thread
//! count.
//!
//! # What the gate checks
//!
//! * the traversal is complete: the record count of this census equals the count
//!   reached independently through `magnetic_isotropy_subgroups_of`;
//! * **the actual TSV rows**, parsed per space group, are exactly one row per
//!   expected `(sg, condensing irrep, direction index, UNI, source, η)` key with
//!   no duplicate, unknown or missing key, and the parsed rows -- not the
//!   counters -- are what the printed stage totals are compared against (see
//!   [`validate_tsv`]); every record's `(map, shift)` is one fixed pair across
//!   all its sources and both `η`;
//! * the population partitions: records = probed + unprobed, pairs = probed +
//!   record-without-fixed-embedding + source-domain-excluded;
//! * the pinned R25/R26/R27 witnesses are re-run through the same entry points
//!   with their own pinned embeddings, and the `RequireExact` refusal of the
//!   SG 1 `X1`/UNI 3 witness keeps its named variant;
//! * the two η probes of a pair agree on their status;
//! * a two-dimensional Γ source exists in the tables and the contract refuses it
//!   with `NotScalarSource` (the boundary evidence of the next card); the first
//!   such source is reported together with whether it is a compound/realified row
//!   and with the matrix data of the first ordinary `dim = 2` source.
//!
//! The gate fails on a violated identity, not on an unfavourable measurement: a
//! real support failure is reported as a number and a named variant.

use cryspglib::irrep::magnetic_embedding::{
    ParentEmbedding, geometry_of, magnetic_operations, parent_embedding_verdict,
};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, MagneticInputError, RealSign, RecordIdentity,
    ScalarGammaModel, ScalarGammaRequest, SourceIdentity, verify_scalar_gamma_input,
};
use cryspglib::irrep::magnetic_subduction::{
    MagneticSubduction, MagneticSubductionError, subduce_scalar_gamma,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::types::{CharacterViewError, IrrepRecord};
use cryspglib::irrep::{LabelConvention, isotropy, magnetic_summary};
use num_complex::Complex64;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::process::ExitCode;

/// Denominators a real `±1` character is snapped to, mirroring the R25
/// contract (`magnetic_input::MAX_EXACT_DENOMINATOR`).
const MAX_EXACT_DENOMINATOR: i128 = 12;

/// Tolerance of that snap, mirroring the R25 contract.
const GRID_TOLERANCE: f64 = 1e-9;

/// Γ in the generated k-vector encoding: the numerators vanish (the contract's
/// own test), whatever the common denominator is.
const GAMMA_NUMERATORS: [i8; 3] = [0, 0, 0];

/// Γ as the magnetic summary encodes a k-point's coordinates.
const GAMMA: (i8, i8, i8, i8) = (0, 0, 0, 1);

/// A second ordinary `dim = 2` Γ source, named explicitly as a cross-check of
/// the scan: a trigonal parent whose Γ `GM3` is a genuine (non-compound)
/// two-dimensional row.  The scan itself is not limited to it.
const NAMED_TWO_DIMENSIONAL_GAMMA: (u8, &str) = (156, "GM3");

// ---------------------------------------------------------------------------
// Inventory: records and sources of one parent space group
// ---------------------------------------------------------------------------

/// One Γ source of a parent, with the domain screen's verdict.
#[derive(Clone)]
struct SourceRow {
    /// Canonical CDML label, as the contract reports it back.
    label: &'static str,
    /// Generated dimension.
    dim: u8,
    /// Whether the k-vector numerators vanish.
    is_gamma: bool,
    /// `None` when the source is inside the probe domain; otherwise the named
    /// screen that excludes it.
    exclusion: Option<String>,
}

/// One magnetic isotropy record of a parent.
#[derive(Clone, Copy)]
struct RecordRow {
    /// CDML label of the **condensing** irrep that owns the row.
    condensing: &'static str,
    /// Position of the row in that irrep's magnetic isotropy list.
    direction_index: usize,
    /// UNI of the magnetic subgroup the row describes.
    uni: usize,
}

/// Everything the inventory reaches for one parent space group.
#[derive(Clone, Default)]
struct SgInventory {
    sg: u8,
    records: Vec<RecordRow>,
    sources: Vec<SourceRow>,
    /// Irreps that are not at Γ (**spinor ones included**): out of scope by
    /// definition, counted only so the source inventory is traceable.
    irreps_not_gamma: usize,
    /// Γ irreps whose `k_label()` is `"GM"`: the label-based reading, printed
    /// next to the data-based one, never used to decide the domain.
    gamma_by_label: usize,
    /// Γ irreps by the contract's own test (`numerators == [0, 0, 0]`).
    gamma_by_data: usize,
    /// Records reached independently through `magnetic_isotropy_subgroups_of`.
    records_by_query: usize,
}

/// The independent record reading, for the traversal-completeness gate.
fn records_by_query(sg: u8) -> usize {
    query::magnetic_isotropy_subgroups_of(sg).len()
}

/// Build the inventory of one parent space group.  No probe runs here.
fn inventory_of(sg: u8) -> SgInventory {
    let mut out = SgInventory {
        sg,
        ..SgInventory::default()
    };
    for irrep in query::irreps_of(sg) {
        if !irrep.spinor {
            for (direction_index, record) in irrep.magnetic_subgroups().iter().enumerate() {
                out.records.push(RecordRow {
                    condensing: irrep.ml,
                    direction_index,
                    uni: record.mag_sg,
                });
            }
        }
        let is_gamma = irrep.k_vector().numerators == GAMMA_NUMERATORS;
        if is_gamma {
            out.gamma_by_data += 1;
        }
        if irrep.k_label() == "GM" {
            out.gamma_by_label += 1;
        }
        if irrep.spinor || !is_gamma {
            // Spinor irreps are not Γ sources of this scalar family, and a
            // non-Γ irrep is a different k-point: both are inventory entries,
            // not probe-domain candidates.
            if !is_gamma {
                out.irreps_not_gamma += 1;
            }
            if is_gamma {
                out.sources.push(SourceRow {
                    label: irrep.ml,
                    dim: irrep.dim,
                    is_gamma,
                    exclusion: Some("spinor_source".to_string()),
                });
            }
            continue;
        }
        out.sources.push(classify_scalar_gamma_source(irrep));
    }
    out.records_by_query = records_by_query(sg);
    out
}

/// The named domain screen of one scalar Γ irrep: non-spinor, one-dimensional,
/// real `±1` character.
fn classify_scalar_gamma_source(irrep: &'static IrrepRecord) -> SourceRow {
    let row = SourceRow {
        label: irrep.ml,
        dim: irrep.dim,
        is_gamma: true,
        exclusion: None,
    };
    if irrep.dim != 1 {
        return SourceRow {
            exclusion: Some("dimension_gt_one".to_string()),
            ..row
        };
    }
    match irrep.ordinary_scalar_selected_arm_block_trace() {
        Err(error) => SourceRow {
            exclusion: Some(format!(
                "character_unavailable_{}",
                character_view_variant(&error)
            )),
            ..row
        },
        Ok(characters) => {
            for value in characters.values() {
                let Some(real) = exact_fraction(value.re) else {
                    return SourceRow {
                        exclusion: Some("character_off_grid".to_string()),
                        ..row
                    };
                };
                let Some(imaginary) = exact_fraction(value.im) else {
                    return SourceRow {
                        exclusion: Some("character_off_grid".to_string()),
                        ..row
                    };
                };
                let negative_one = Rat::ONE.checked_neg().expect("one is representable");
                if !imaginary.is_zero() || (real != Rat::ONE && real != negative_one) {
                    return SourceRow {
                        exclusion: Some("non_real_character".to_string()),
                        ..row
                    };
                }
            }
            row
        }
    }
}

/// The contract's own exact small-denominator snap, mirrored here so the screen
/// and the contract read the same grid.
fn exact_fraction(value: f64) -> Option<Rat> {
    if !value.is_finite() {
        return None;
    }
    for denominator in 1..=MAX_EXACT_DENOMINATOR {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < GRID_TOLERANCE {
            return Rat::new(scaled.round() as i128, denominator).ok();
        }
    }
    None
}

/// Name of a `CharacterViewError` variant, for the exclusion histogram.
fn character_view_variant(error: &CharacterViewError) -> &'static str {
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

// ---------------------------------------------------------------------------
// The fixed embedding of a record
// ---------------------------------------------------------------------------

/// A chosen embedding: the map and the exact origin shift.
struct FixedCandidate {
    map: Mat3R,
    shift: [Rat; 3],
}

/// The record's fixed embedding, or why there is none.
enum FixedEmbedding {
    /// The first realising candidate under the total key, with `RequireExact`.
    Fixed(Box<FixedCandidate>),
    /// The search returned contained entries, none of which realises the
    /// record's lattice exactly.
    NoRealising,
    /// The search returned nothing and every candidate was refuted exactly.
    CertifiedAbsence,
    /// The search returned nothing and some candidate is **undecided**
    /// (arithmetic limit or component cap): not a proof of absence.
    UndecidedResidual,
    /// The record could not even be prepared; the reason is named.
    Error {
        variant: &'static str,
        detail: String,
    },
}

/// The exact entries of a candidate map, as `(numerator, denominator)` pairs.
type ExactMapKey = [(i128, i128); 9];

/// The exact entries of a candidate shift, as `(numerator, denominator)` pairs.
type ExactShiftKey = [(i128, i128); 3];

/// The total, source-independent ordering key of one returned entry.
///
/// The first component prefers the **tabulated frame** (`map = I`, with whatever
/// shift the search solved): when the record is tabulated directly in the
/// parent's conventional frame, that is the embedding every accepted R25/R26
/// witness uses, so the census's frame choice stays aligned with them.  The
/// remaining components order the other candidates by their exact entries, so
/// the choice is a function of the record alone and never of the source or of
/// the subduction outcome.
fn embedding_key(entry: &ParentEmbedding) -> (bool, ExactMapKey, ExactShiftKey) {
    let mut map = [(0i128, 1i128); 9];
    for row in 0..3 {
        for column in 0..3 {
            map[row * 3 + column] = value_of(entry.map, row, column);
        }
    }
    let mut shift = [(0i128, 1i128); 3];
    for (axis, value) in entry.shift.iter().enumerate() {
        shift[axis] = (value.numerator(), value.denominator());
    }
    (entry.map != Mat3R::identity(), map, shift)
}

/// One entry of a `Mat3R`, as `(numerator, denominator)`.
fn value_of(matrix: Mat3R, row: usize, column: usize) -> (i128, i128) {
    let value = matrix.row(row)[column];
    (value.numerator(), value.denominator())
}

/// Pick the record's fixed embedding.  The choice reads only the record.
fn fixed_embedding(sg: u8, record: &RecordRow) -> FixedEmbedding {
    let geometry = match geometry_of(sg, &record_of(sg, record)) {
        Ok(geometry) => geometry,
        Err(error) => {
            return FixedEmbedding::Error {
                variant: "geometry",
                detail: sanitize(&error.to_string()),
            };
        }
    };
    let set = match magnetic_operations(record.uni) {
        Ok(set) => set,
        Err(error) => {
            return FixedEmbedding::Error {
                variant: "operations",
                detail: sanitize(&error.to_string()),
            };
        }
    };
    let verdict = match parent_embedding_verdict(&geometry, &set) {
        Ok(verdict) => verdict,
        Err(error) => {
            return FixedEmbedding::Error {
                variant: "verdict",
                detail: sanitize(&error.to_string()),
            };
        }
    };
    let mut realising: Vec<&ParentEmbedding> = verdict
        .embeddings
        .iter()
        .filter(|entry| entry.realises_record_lattice)
        .collect();
    if realising.is_empty() {
        return if !verdict.embeddings.is_empty() {
            FixedEmbedding::NoRealising
        } else if verdict.absence_is_certified {
            FixedEmbedding::CertifiedAbsence
        } else {
            FixedEmbedding::UndecidedResidual
        };
    }
    realising.sort_by_key(|entry| embedding_key(entry));
    let chosen = realising[0];
    FixedEmbedding::Fixed(Box::new(FixedCandidate {
        map: chosen.map,
        shift: chosen.shift,
    }))
}

/// The record the inventory row points at, read back from the database.
fn record_of(sg: u8, record: &RecordRow) -> cryspglib::irrep::types::MagneticIsotropyRecord {
    isotropy::magnetic_isotropy_subgroups(sg, record.condensing, LabelConvention::Cdml)
        .expect("the inventory's condensing label resolves")[record.direction_index]
        .record
}

// ---------------------------------------------------------------------------
// Probes
// ---------------------------------------------------------------------------

/// One probe's outcome: exactly one stage.
#[derive(Clone)]
enum Probe {
    /// The R25 input contract refused it.
    VerifyInput {
        variant: &'static str,
        detail: String,
    },
    /// The R26 subduction refused a model the contract accepted.
    Subduce {
        variant: &'static str,
        detail: String,
    },
    /// Both answered, but the reconstruction did not hold.  Must never happen.
    SuccessCheck { variant: String, detail: String },
    /// Contract, subduction and class-by-class reconstruction all hold.
    Success {
        target_label: String,
        target_sg: u8,
        multiplicity: u32,
        operators: Vec<cryspglib::irrep::magnetic_subduction::TargetOperator>,
        intertwiner: Complex64,
    },
}

impl Probe {
    /// The stage name, one per probe.
    fn stage(&self) -> &'static str {
        match self {
            Probe::VerifyInput { .. } => "verify_input",
            Probe::Subduce { .. } => "subduce",
            Probe::SuccessCheck { .. } => "success_check",
            Probe::Success { .. } => "success",
        }
    }

    /// The named variant or `ok`.
    fn variant(&self) -> String {
        match self {
            Probe::VerifyInput { variant, .. } | Probe::Subduce { variant, .. } => {
                (*variant).to_string()
            }
            Probe::SuccessCheck { variant, .. } => variant.clone(),
            Probe::Success { .. } => "ok".to_string(),
        }
    }

    /// The single-line detail.
    fn detail(&self) -> String {
        match self {
            Probe::VerifyInput { detail, .. }
            | Probe::Subduce { detail, .. }
            | Probe::SuccessCheck { detail, .. } => detail.clone(),
            Probe::Success { .. } => "-".to_string(),
        }
    }
}

/// Run one real probe: R25 contract, then R26 subduction, then the class-level
/// reconstruction check.
fn probe(
    sg: u8,
    record: &RecordRow,
    source: &'static str,
    eta: RealSign,
    map: Mat3R,
    shift: [Rat; 3],
) -> Probe {
    let request = ScalarGammaRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: sg,
            condensing_label: record.condensing,
            convention: LabelConvention::Cdml,
            direction_index: record.direction_index,
        },
        Some(eta),
    )
    .with_embedding(EmbeddingEvidence {
        map,
        shift,
        lattice_policy: LatticePolicy::RequireExact,
    });
    let model = match verify_scalar_gamma_input(&request) {
        Ok(model) => model,
        Err(error) => {
            return Probe::VerifyInput {
                variant: input_variant(&error),
                detail: sanitize(&error.to_string()),
            };
        }
    };
    let report = match subduce_scalar_gamma(&model) {
        Ok(report) => report,
        Err(error) => {
            return Probe::Subduce {
                variant: subduction_variant(&error),
                detail: sanitize(&error.to_string()),
            };
        }
    };
    if let Some(failure) = reconstruction_failure(&model, &report) {
        return Probe::SuccessCheck {
            variant: "reconstruction".to_string(),
            detail: failure,
        };
    }
    Probe::Success {
        target_label: report.target_label().to_string(),
        target_sg: report.target_unitary_sg(),
        multiplicity: report.multiplicity(),
        operators: report.operators().to_vec(),
        intertwiner: report.intertwiner(),
    }
}

/// The independent reconstruction check on a successful probe.  Returns the
/// name of the first violated claim, or `None`.
fn reconstruction_failure(
    model: &ScalarGammaModel,
    report: &MagneticSubduction<'_>,
) -> Option<String> {
    if report.operators().len() != model.classes().len() {
        return Some("the operator table has one entry per class".to_string());
    }
    if report.multiplicity() == 0
        || report.multiplicity() * report.target_dimension() != report.source_dimension()
        || report.covered_dimension() != report.source_dimension()
    {
        return Some("dim(source) = multiplicity * dim(target)".to_string());
    }
    // A positive predicate: a NaN or non-finite `U` must fail here rather than
    // slip through the complement of an ordinary comparison.
    let intertwiner = report.intertwiner();
    if !is_finite_unit(intertwiner.norm_sqr()) {
        return Some("the intertwiner is a finite non-zero unit scalar".to_string());
    }
    for (index, action) in model.actions().iter().enumerate() {
        if report.operators()[index].is_antiunitary() != action.class.time_reversal {
            return Some(format!("class {index} keeps its time-reversal flag"));
        }
        for probe_value in [Complex64::new(1.0, 2.0), Complex64::new(-0.7, 0.4)] {
            let coordinate = intertwiner * probe_value;
            let source_side = Complex64::new(f64::from(action.value.value()), 0.0)
                * if action.class.time_reversal {
                    coordinate.conj()
                } else {
                    coordinate
                };
            let target_side = intertwiner * report.operators()[index].apply(probe_value);
            // Also a positive predicate, so a NaN operator value is a failure.
            if !within_tolerance((source_side - target_side).norm()) {
                return Some(format!(
                    "class {index}: A(g)(U z) = U D(g) z failed for z = {probe_value:?}"
                ));
            }
        }
    }
    None
}

/// Whether a squared modulus is a finite non-zero unit within tolerance.
///
/// Written as a positive predicate so `NaN` and infinities are failures; the
/// complement of `|x - 1| < eps` would accept a `NaN`.
fn is_finite_unit(value: f64) -> bool {
    value.is_finite() && value > 0.0 && (value - 1.0).abs() < 1e-12
}

/// Whether a residual is finite and within the reconstruction tolerance.
fn within_tolerance(value: f64) -> bool {
    value.is_finite() && value.abs() <= 1e-9
}

/// Name of a `MagneticInputError` variant.
fn input_variant(error: &MagneticInputError) -> &'static str {
    match error {
        MagneticInputError::UnknownParentSpaceGroup { .. } => "UnknownParentSpaceGroup",
        MagneticInputError::LabelNotFound { .. } => "LabelNotFound",
        MagneticInputError::AmbiguousLabel { .. } => "AmbiguousLabel",
        MagneticInputError::SpinorSource { .. } => "SpinorSource",
        MagneticInputError::NotScalarSource { .. } => "NotScalarSource",
        MagneticInputError::SourceNotGamma { .. } => "SourceNotGamma",
        MagneticInputError::ParentMismatch { .. } => "ParentMismatch",
        MagneticInputError::DirectionOutOfRange { .. } => "DirectionOutOfRange",
        MagneticInputError::DatabaseMismatch(_) => "DatabaseMismatch",
        MagneticInputError::GroupProof(_) => "GroupProof",
        MagneticInputError::ClassListInconsistent { .. } => "ClassListInconsistent",
        MagneticInputError::OffGrid { .. } => "OffGrid",
        MagneticInputError::ParentTranslationOffGrid { .. } => "ParentTranslationOffGrid",
        MagneticInputError::NonRealCharacter { .. } => "NonRealCharacter",
        MagneticInputError::IdentityCharacterNotOne { .. } => "IdentityCharacterNotOne",
        MagneticInputError::CharacterRowSettingMismatch { .. } => "CharacterRowSettingMismatch",
        MagneticInputError::EmbeddingRejected(_) => "EmbeddingRejected",
        MagneticInputError::LatticePolicyUnsatisfied { .. } => "LatticePolicyUnsatisfied",
        MagneticInputError::MissingTimeReversal { .. } => "MissingTimeReversal",
        MagneticInputError::NotACorepresentation { .. } => "NotACorepresentation",
        MagneticInputError::Arithmetic(_) => "Arithmetic",
    }
}

/// Name of a `MagneticSubductionError` variant.
fn subduction_variant(error: &MagneticSubductionError) -> &'static str {
    match error {
        MagneticSubductionError::SourceIdentityMismatch { .. } => "SourceIdentityMismatch",
        MagneticSubductionError::SummaryFailed { .. } => "SummaryFailed",
        MagneticSubductionError::NoGammaKPoint { .. } => "NoGammaKPoint",
        MagneticSubductionError::UnresolvedTargetData { .. } => "UnresolvedTargetData",
        MagneticSubductionError::TargetRawOperationsUnavailable { .. } => {
            "TargetRawOperationsUnavailable"
        }
        MagneticSubductionError::TargetRawOperationsDiffer { .. } => "TargetRawOperationsDiffer",
        MagneticSubductionError::TargetOperationCountMismatch { .. } => {
            "TargetOperationCountMismatch"
        }
        MagneticSubductionError::TargetOperationIndexOutOfRange { .. } => {
            "TargetOperationIndexOutOfRange"
        }
        MagneticSubductionError::TargetOperationIndexRepeated { .. } => {
            "TargetOperationIndexRepeated"
        }
        MagneticSubductionError::TargetTimeReversalMismatch { .. } => "TargetTimeReversalMismatch",
        MagneticSubductionError::TargetColumnOutOfRange { .. } => "TargetColumnOutOfRange",
        MagneticSubductionError::TargetColumnRepeated { .. } => "TargetColumnRepeated",
        MagneticSubductionError::TargetColumnMisassociated { .. } => "TargetColumnMisassociated",
        MagneticSubductionError::TargetOperationIndexMisassociated { .. } => {
            "TargetOperationIndexMisassociated"
        }
        MagneticSubductionError::TargetOperationUnmatched { .. } => "TargetOperationUnmatched",
        MagneticSubductionError::TargetClassUncovered { .. } => "TargetClassUncovered",
        MagneticSubductionError::TargetTranslationOffGrid { .. } => "TargetTranslationOffGrid",
        MagneticSubductionError::TargetCharacterUnavailable { .. } => "TargetCharacterUnavailable",
        MagneticSubductionError::TargetCharacterNotAClassFunction { .. } => {
            "TargetCharacterNotAClassFunction"
        }
        MagneticSubductionError::TargetSquareInconsistent { .. } => "TargetSquareInconsistent",
        MagneticSubductionError::TargetCosetInconsistent { .. } => "TargetCosetInconsistent",
        MagneticSubductionError::TargetConjugationInconsistent { .. } => {
            "TargetConjugationInconsistent"
        }
        MagneticSubductionError::TargetLawViolated { .. } => "TargetLawViolated",
        MagneticSubductionError::NoMatchingTarget { .. } => "NoMatchingTarget",
        MagneticSubductionError::AmbiguousTarget { .. } => "AmbiguousTarget",
        MagneticSubductionError::DimensionNotConserved { .. } => "DimensionNotConserved",
        MagneticSubductionError::Arithmetic { .. } => "Arithmetic",
    }
}

/// One line, tab-free, so the TSV stays parseable.
fn sanitize(text: &str) -> String {
    text.replace(['\t', '\n', '\r'], " ")
}

// ---------------------------------------------------------------------------
// Per-space-group census
// ---------------------------------------------------------------------------

/// Counters of the probe pass, merged across space groups.
#[derive(Clone, Default)]
struct Counters {
    records_total: usize,
    records_fixed: usize,
    /// Of the fixed embeddings, those whose map is exactly the identity (the
    /// tabulated frame); the shift may still be a solved non-zero representative.
    records_fixed_identity_map: usize,
    /// Of those, the ones whose solved shift representative is zero.
    records_fixed_identity_zero_shift: usize,
    records_no_fixed: usize,
    records_no_realising: usize,
    records_certified_absence: usize,
    records_undecided_residual: usize,
    records_error: usize,
    sources_gamma: usize,
    sources_eligible: usize,
    sources_spinor: usize,
    sources_dimension_gt_one: usize,
    sources_character_unavailable: usize,
    sources_non_real_character: usize,
    sources_character_off_grid: usize,
    pairs_total: usize,
    pairs_probed: usize,
    pairs_record_without_fixed: usize,
    pairs_source_domain_excluded: usize,
    probes: usize,
    success: usize,
    verify_input_failures: usize,
    subduce_failures: usize,
    success_check_failures: usize,
    tsv_rows: usize,
    tsv_rows_record_without_fixed: usize,
    eta_status_mismatches: usize,
    eta_target_label_differs: usize,
    eta_operator_table_differs: usize,
    eta_intertwiner_differs: usize,
    verify_variants: BTreeMap<String, usize>,
    subduce_variants: BTreeMap<String, usize>,
    no_fixed_variants: BTreeMap<String, usize>,
    exclusion_variants: BTreeMap<String, usize>,
}

impl Counters {
    /// Fold another space group's counters in.  BTreeMaps make the merge order
    /// irrelevant, so the printed histograms are thread-count independent.
    fn merge(&mut self, other: &Counters) {
        macro_rules! add {
            ($($field:ident),*) => { $( self.$field += other.$field; )* };
        }
        add!(
            records_total,
            records_fixed,
            records_fixed_identity_map,
            records_fixed_identity_zero_shift,
            records_no_fixed,
            records_no_realising,
            records_certified_absence,
            records_undecided_residual,
            records_error,
            sources_gamma,
            sources_eligible,
            sources_spinor,
            sources_dimension_gt_one,
            sources_character_unavailable,
            sources_non_real_character,
            sources_character_off_grid,
            pairs_total,
            pairs_probed,
            pairs_record_without_fixed,
            pairs_source_domain_excluded,
            probes,
            success,
            verify_input_failures,
            subduce_failures,
            success_check_failures,
            tsv_rows,
            tsv_rows_record_without_fixed,
            eta_status_mismatches,
            eta_target_label_differs,
            eta_operator_table_differs,
            eta_intertwiner_differs
        );
        for (target, source) in [
            (&mut self.verify_variants, &other.verify_variants),
            (&mut self.subduce_variants, &other.subduce_variants),
            (&mut self.no_fixed_variants, &other.no_fixed_variants),
            (&mut self.exclusion_variants, &other.exclusion_variants),
        ] {
            for (key, value) in source {
                *target.entry(key.clone()).or_insert(0) += value;
            }
        }
    }

    /// One line of the per-space-group table.
    fn line(&self, sg: u8) -> String {
        format!(
            "sg {sg} records {} fixed {} nofix {} sources {} probes {} success {} \
             verify_input {} subduce {}",
            self.records_total,
            self.records_fixed,
            self.records_no_fixed,
            self.sources_eligible,
            self.probes,
            self.success,
            self.verify_input_failures,
            self.subduce_failures
        )
    }
}

/// One space group's census: counters and its TSV chunk.
#[derive(Clone, Default)]
struct SgReport {
    sg: u8,
    counters: Counters,
    tsv: String,
}

/// Census one space group: fix one embedding per record, probe each eligible
/// source under both η, and emit the TSV rows in database order.
fn census_sg(inventory: &SgInventory) -> SgReport {
    let sg = inventory.sg;
    let mut counters = Counters::default();
    let mut tsv = String::new();
    let eligible: Vec<&SourceRow> = inventory
        .sources
        .iter()
        .filter(|source| source.exclusion.is_none())
        .collect();
    counters.records_total = inventory.records.len();
    counters.sources_gamma = inventory.sources.iter().filter(|s| s.is_gamma).count();
    counters.sources_eligible = eligible.len();
    for source in &inventory.sources {
        if let Some(exclusion) = &source.exclusion {
            *counters
                .exclusion_variants
                .entry(exclusion.clone())
                .or_insert(0) += 1;
            if exclusion == "spinor_source" {
                counters.sources_spinor += 1;
            } else if exclusion == "dimension_gt_one" {
                counters.sources_dimension_gt_one += 1;
            } else if exclusion == "non_real_character" {
                counters.sources_non_real_character += 1;
            } else if exclusion == "character_off_grid" {
                counters.sources_character_off_grid += 1;
            } else if exclusion.starts_with("character_unavailable_") {
                counters.sources_character_unavailable += 1;
            }
        }
    }
    counters.pairs_total = counters.records_total * counters.sources_gamma;
    counters.pairs_source_domain_excluded =
        counters.records_total * (counters.sources_gamma - counters.sources_eligible);

    // Records without a fixed embedding are collected first: their TSV rows
    // carry the same (record, source, η) shape as a probe, with a named reason.
    let mut no_fixed: Vec<(RecordRow, &'static str, String)> = Vec::new();
    for record in &inventory.records {
        match fixed_embedding(sg, record) {
            FixedEmbedding::Fixed(candidate) => {
                let map = candidate.map;
                let shift = candidate.shift;
                counters.records_fixed += 1;
                if map == Mat3R::identity() {
                    counters.records_fixed_identity_map += 1;
                    if shift.iter().all(|value| value.is_zero()) {
                        counters.records_fixed_identity_zero_shift += 1;
                    }
                }
                for source in &eligible {
                    counters.pairs_probed += 1;
                    let plus = probe(sg, record, source.label, RealSign::Plus, map, shift);
                    let minus = probe(sg, record, source.label, RealSign::Minus, map, shift);
                    count_probe(&mut counters, &plus);
                    count_probe(&mut counters, &minus);
                    compare_eta(&mut counters, &plus, &minus);
                    push_row(
                        &mut tsv,
                        sg,
                        record,
                        source.label,
                        RealSign::Plus,
                        &plus,
                        map,
                        shift,
                    );
                    push_row(
                        &mut tsv,
                        sg,
                        record,
                        source.label,
                        RealSign::Minus,
                        &minus,
                        map,
                        shift,
                    );
                }
            }
            other => {
                let (variant, detail) = classify_no_fixed(&other);
                *counters
                    .no_fixed_variants
                    .entry(variant.to_string())
                    .or_insert(0) += 1;
                counters.records_no_fixed += 1;
                match other {
                    FixedEmbedding::NoRealising => counters.records_no_realising += 1,
                    FixedEmbedding::CertifiedAbsence => counters.records_certified_absence += 1,
                    FixedEmbedding::UndecidedResidual => counters.records_undecided_residual += 1,
                    FixedEmbedding::Error { .. } => counters.records_error += 1,
                    FixedEmbedding::Fixed(_) => unreachable!("handled above"),
                }
                no_fixed.push((*record, variant, detail));
            }
        }
    }

    for (record, variant, detail) in &no_fixed {
        counters.pairs_record_without_fixed += eligible.len();
        for source in &eligible {
            for eta in [RealSign::Plus, RealSign::Minus] {
                push_no_fixed_row(&mut tsv, sg, record, source.label, eta, variant, detail);
            }
        }
    }
    counters.tsv_rows = counters.pairs_probed * 2 + counters.pairs_record_without_fixed * 2;
    counters.tsv_rows_record_without_fixed = counters.pairs_record_without_fixed * 2;

    SgReport { sg, counters, tsv }
}

/// Attribute one probe to its single stage.
fn count_probe(counters: &mut Counters, outcome: &Probe) {
    counters.probes += 1;
    match outcome {
        Probe::Success { .. } => counters.success += 1,
        Probe::VerifyInput { variant, .. } => {
            counters.verify_input_failures += 1;
            *counters
                .verify_variants
                .entry((*variant).to_string())
                .or_insert(0) += 1;
        }
        Probe::Subduce { variant, .. } => {
            counters.subduce_failures += 1;
            *counters
                .subduce_variants
                .entry((*variant).to_string())
                .or_insert(0) += 1;
        }
        Probe::SuccessCheck { variant, .. } => {
            counters.success_check_failures += 1;
            *counters.verify_variants.entry(variant.clone()).or_insert(0) += 1;
        }
    }
}

/// The two η probes must agree on their status; their target identity and
/// intertwiner are reported separately because the antiunitary sign may select
/// a different (still valid) representative.
fn compare_eta(counters: &mut Counters, plus: &Probe, minus: &Probe) {
    match (plus, minus) {
        (
            Probe::Success {
                target_label: plus_label,
                target_sg: plus_sg,
                multiplicity: plus_mult,
                operators: plus_ops,
                intertwiner: plus_u,
            },
            Probe::Success {
                target_label: minus_label,
                target_sg: minus_sg,
                multiplicity: minus_mult,
                operators: minus_ops,
                intertwiner: minus_u,
            },
        ) => {
            if plus_label != minus_label || plus_sg != minus_sg || plus_mult != minus_mult {
                counters.eta_target_label_differs += 1;
            }
            if plus_ops != minus_ops {
                counters.eta_operator_table_differs += 1;
            }
            if plus_u != minus_u {
                counters.eta_intertwiner_differs += 1;
            }
        }
        (Probe::VerifyInput { variant: a, .. }, Probe::VerifyInput { variant: b, .. })
            if a == b => {}
        (Probe::Subduce { variant: a, .. }, Probe::Subduce { variant: b, .. }) if a == b => {}
        (Probe::SuccessCheck { variant: a, .. }, Probe::SuccessCheck { variant: b, .. })
            if a == b => {}
        _ => counters.eta_status_mismatches += 1,
    }
}

/// The named reason a record has no fixed embedding.
fn classify_no_fixed(embedding: &FixedEmbedding) -> (&'static str, String) {
    match embedding {
        FixedEmbedding::NoRealising => (
            "no_realising_candidate",
            "the search returned contained entries, none realising the record lattice".to_string(),
        ),
        FixedEmbedding::CertifiedAbsence => (
            "certified_absence",
            "every candidate map of the R24 family was refuted by an exact condition of an \
             embedding; this is a statement about that family, not about all settings"
                .to_string(),
        ),
        FixedEmbedding::UndecidedResidual => (
            "undecided_residual",
            "some candidate map is undecided (arithmetic limit or component cap): the answer is \
             unknown, not a proof of absence"
                .to_string(),
        ),
        FixedEmbedding::Error { variant, detail } => (*variant, detail.clone()),
        FixedEmbedding::Fixed(_) => unreachable!("fixed embeddings are probed"),
    }
}

/// One TSV row of a probed or unprobed pair.
#[allow(clippy::too_many_arguments)]
fn push_row(
    tsv: &mut String,
    sg: u8,
    record: &RecordRow,
    source: &'static str,
    eta: RealSign,
    outcome: &Probe,
    map: Mat3R,
    shift: [Rat; 3],
) {
    let mut row = String::new();
    row.push_str(&format!(
        "{sg}\t{}\t{}\t{}\t{source}\t{:+}\t{}",
        record.condensing,
        record.direction_index,
        record.uni,
        eta.value(),
        outcome.stage()
    ));
    row.push('\t');
    row.push_str(&outcome.variant());
    row.push('\t');
    row.push_str(&format_map(map));
    row.push('\t');
    row.push_str(&format_shift(shift));
    row.push('\t');
    match outcome {
        Probe::Success {
            target_label,
            target_sg,
            multiplicity,
            operators,
            intertwiner,
        } => {
            row.push_str(&format!(
                "{target_label}\t{target_sg}\t{multiplicity}\t{:.17}\t{:.17}\t{op_count}",
                intertwiner.re,
                intertwiner.im,
                op_count = operators.len()
            ));
        }
        _ => row.push_str("-\t-\t-\t-\t-\t-"),
    }
    row.push('\t');
    row.push_str(&outcome.detail());
    row.push('\n');
    tsv.push_str(&row);
}

/// One TSV row of a pair whose record has no fixed embedding.
fn push_no_fixed_row(
    tsv: &mut String,
    sg: u8,
    record: &RecordRow,
    source: &'static str,
    eta: RealSign,
    variant: &str,
    detail: &str,
) {
    tsv.push_str(&format!(
        "{sg}\t{}\t{}\t{}\t{source}\t{:+}\tno_fixed_embedding\t{variant}\t-\t-\t-\t-\t-\t-\t-\t-\t{detail}\n",
        record.condensing,
        record.direction_index,
        record.uni,
        eta.value()
    ));
}

/// A matrix as nine `num/den` values, row-major.
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

/// A vector as three `num/den` values.
fn format_shift(shift: [Rat; 3]) -> String {
    shift
        .iter()
        .map(|value| format!("{}/{}", value.numerator(), value.denominator()))
        .collect::<Vec<_>>()
        .join(",")
}

// ---------------------------------------------------------------------------
// TSV validation: the gate reads the actual rows, not the counters
// ---------------------------------------------------------------------------

/// Number of tab-separated fields of one data row.
const TSV_FIELDS: usize = 17;

/// The stage names a data row may carry.
const STAGES: [&str; 5] = [
    "no_fixed_embedding",
    "verify_input",
    "subduce",
    "success_check",
    "success",
];

/// The named reasons a record may have no fixed embedding.
const NO_FIXED_REASONS: [&str; 6] = [
    "no_realising_candidate",
    "certified_absence",
    "undecided_residual",
    "geometry",
    "operations",
    "verdict",
];

/// The identity of one TSV row: the exact `(record, source, η)` key.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RowKey {
    sg: u8,
    condensing: String,
    direction_index: usize,
    uni: usize,
    source: String,
    eta: i32,
}

/// One parsed data row.
struct ParsedRow {
    key: RowKey,
    stage: String,
    /// The variant column, read from the row itself.
    variant: String,
    map: String,
    shift: String,
}

/// Parse TSV data rows, reporting every line that is not well formed.
fn parse_tsv(text: &str) -> (Vec<ParsedRow>, Vec<String>) {
    let mut rows = Vec::new();
    let mut errors = Vec::new();
    for (line_number, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with("sg\tcondensing\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != TSV_FIELDS {
            errors.push(format!(
                "line {}: {} fields, expected {TSV_FIELDS}",
                line_number + 1,
                fields.len()
            ));
            continue;
        }
        let (Ok(sg), Ok(direction_index), Ok(uni)) = (
            fields[0].parse::<u8>(),
            fields[2].parse::<usize>(),
            fields[3].parse::<usize>(),
        ) else {
            errors.push(format!(
                "line {}: unreadable numeric fields",
                line_number + 1
            ));
            continue;
        };
        let eta = match fields[5] {
            "+1" => 1i32,
            "-1" => -1i32,
            other => {
                errors.push(format!(
                    "line {}: eta {other:?} is neither +1 nor -1",
                    line_number + 1
                ));
                continue;
            }
        };
        rows.push(ParsedRow {
            key: RowKey {
                sg,
                condensing: fields[1].to_string(),
                direction_index,
                uni,
                source: fields[4].to_string(),
                eta,
            },
            stage: fields[6].to_string(),
            variant: fields[7].to_string(),
            map: fields[8].to_string(),
            shift: fields[9].to_string(),
        });
    }
    (rows, errors)
}

/// The rows one space group must produce: one per `(record, eligible source, η)`.
fn expected_keys(inventory: &SgInventory) -> Vec<RowKey> {
    let mut keys = Vec::new();
    for record in &inventory.records {
        for source in &inventory.sources {
            if source.exclusion.is_some() {
                continue;
            }
            for eta in [1i32, -1i32] {
                keys.push(RowKey {
                    sg: inventory.sg,
                    condensing: record.condensing.to_string(),
                    direction_index: record.direction_index,
                    uni: record.uni,
                    source: source.label.to_string(),
                    eta,
                });
            }
        }
    }
    keys
}

/// Validate one space group's actual TSV rows against the expected key set and
/// the counters.  Everything here reads the emitted text: a dropped row, a
/// duplicated η, an unknown or missing key, a record whose `(map, shift)` is not
/// one fixed pair, or a stage count that disagrees with the counters, is a
/// violation.
fn validate_tsv(inventory: &SgInventory, text: &str, counters: &Counters) -> Vec<String> {
    let sg = inventory.sg;
    let mut violations = Vec::new();
    let (rows, parse_errors) = parse_tsv(text);
    for error in parse_errors {
        violations.push(format!("SG {sg}: {error}"));
    }
    let mut seen: BTreeMap<RowKey, usize> = BTreeMap::new();
    for row in &rows {
        *seen.entry(row.key.clone()).or_insert(0) += 1;
        if !STAGES.contains(&row.stage.as_str()) {
            violations.push(format!("SG {sg}: unknown stage {:?}", row.stage));
        }
        if row.stage == "no_fixed_embedding" && !NO_FIXED_REASONS.contains(&row.variant.as_str()) {
            violations.push(format!(
                "SG {sg}: unknown no-fixed-embedding reason {:?}",
                row.variant
            ));
        }
    }
    let expected: std::collections::BTreeSet<RowKey> =
        expected_keys(inventory).into_iter().collect();
    for (key, count) in &seen {
        if *count > 1 {
            violations.push(format!("SG {sg}: duplicated row {key:?} x{count}"));
        }
        if !expected.contains(key) {
            violations.push(format!("SG {sg}: unknown row {key:?}"));
        }
    }
    for key in &expected {
        if !seen.contains_key(key) {
            violations.push(format!("SG {sg}: missing row {key:?}"));
        }
    }

    // The counters are compared against what the text actually says.
    let no_fixed_rows = rows
        .iter()
        .filter(|row| row.stage == "no_fixed_embedding")
        .count();
    let probed_rows = rows.len() - no_fixed_rows;
    let stage_count = |stage: &str| rows.iter().filter(|row| row.stage == stage).count();
    if rows.len() != counters.tsv_rows {
        violations.push(format!(
            "SG {sg}: {} parsed rows != tsv_rows {}",
            rows.len(),
            counters.tsv_rows
        ));
    }
    if no_fixed_rows != counters.tsv_rows_record_without_fixed {
        violations.push(format!(
            "SG {sg}: {no_fixed_rows} no-fixed rows != {}",
            counters.tsv_rows_record_without_fixed
        ));
    }
    if probed_rows != counters.pairs_probed * 2 {
        violations.push(format!(
            "SG {sg}: {probed_rows} probed rows != 2 x pairs_probed {}",
            counters.pairs_probed
        ));
    }
    for (stage, expected_count) in [
        ("success", counters.success),
        ("verify_input", counters.verify_input_failures),
        ("subduce", counters.subduce_failures),
        ("success_check", counters.success_check_failures),
    ] {
        let actual = stage_count(stage);
        if actual != expected_count {
            violations.push(format!(
                "SG {sg}: {actual} rows staged {stage} != counter {expected_count}"
            ));
        }
    }

    // One record has one fixed `(map, shift)`, shared by all its sources and
    // both η; a record without one carries `-` in both columns.
    type RecordIdentityKey = (String, usize);
    type RecordShape = (
        std::collections::BTreeSet<(String, String)>,
        std::collections::BTreeSet<bool>,
    );
    let mut per_record: BTreeMap<RecordIdentityKey, RecordShape> = BTreeMap::new();
    for row in &rows {
        let entry = per_record
            .entry((row.key.condensing.clone(), row.key.direction_index))
            .or_default();
        entry.0.insert((row.map.clone(), row.shift.clone()));
        entry.1.insert(row.stage == "no_fixed_embedding");
    }
    for ((condensing, direction_index), (pairs, kinds)) in &per_record {
        if kinds.len() != 1 {
            violations.push(format!(
                "SG {sg}: record {condensing}[{direction_index}] mixes probed and unprobed rows"
            ));
        }
        if kinds.contains(&true) {
            if pairs.len() != 1 || !pairs.contains(&("-".to_string(), "-".to_string())) {
                violations.push(format!(
                    "SG {sg}: record {condensing}[{direction_index}] without a fixed embedding \
                     must carry '-' once, found {pairs:?}"
                ));
            }
        } else {
            if pairs.len() != 1 {
                violations.push(format!(
                    "SG {sg}: record {condensing}[{direction_index}] has {} distinct (map, shift) \
                     values, expected one",
                    pairs.len()
                ));
            }
            if pairs.contains(&("-".to_string(), "-".to_string())) {
                violations.push(format!(
                    "SG {sg}: record {condensing}[{direction_index}] is probed but carries '-'"
                ));
            }
        }
    }
    violations
}

/// Validate the concatenated document as well: no duplicate key across space
/// groups, and the row count is the one the counters report.
fn validate_global_tsv(text: &str, counters: &Counters) -> Vec<String> {
    let mut violations = Vec::new();
    let (rows, parse_errors) = parse_tsv(text);
    for error in parse_errors {
        violations.push(format!("global: {error}"));
    }
    if rows.len() != counters.tsv_rows {
        violations.push(format!(
            "global: {} parsed rows != tsv_rows {}",
            rows.len(),
            counters.tsv_rows
        ));
    }
    let mut seen: std::collections::BTreeSet<RowKey> = std::collections::BTreeSet::new();
    for row in &rows {
        if !seen.insert(row.key.clone()) {
            violations.push(format!("global: duplicated row {:?}", row.key));
        }
    }
    violations
}

// ---------------------------------------------------------------------------
// Pinned witnesses (R25/R26/R27)
// ---------------------------------------------------------------------------

/// Expected outcome of one pinned witness.
enum AnchorExpect {
    /// The contract, the subduction and the reconstruction must all succeed.
    Success {
        /// Target unitary subgroup, when the witness pins it.
        target_sg: Option<u8>,
        /// Target source label of the unitary subgroup, when pinned.
        target_source_label: Option<&'static str>,
    },
    /// The contract must refuse the pinned embedding with this named variant.
    Refused(&'static str),
}

/// One pinned witness, with the embedding the accepted round used.
struct Anchor {
    name: &'static str,
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    uni: usize,
    map: Mat3R,
    policy: LatticePolicy,
    etas: &'static [Option<RealSign>],
    expect: AnchorExpect,
}

/// The witnesses the ledger pins, re-run through the same entry points.
fn anchors() -> Vec<Anchor> {
    let identity = Mat3R::identity();
    let two = Mat3R::from_ints([[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
    let both = &[Some(RealSign::Plus), Some(RealSign::Minus)][..];
    vec![
        Anchor {
            name: "sg2_typeI_uni4",
            sg: 2,
            source: "GM1+",
            condensing: "GM1+",
            uni: 4,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: &[None],
            expect: AnchorExpect::Success {
                target_sg: Some(2),
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg2_typeIII_uni6",
            sg: 2,
            source: "GM1-",
            condensing: "GM1-",
            uni: 6,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg1_typeIV_uni3_containment",
            sg: 1,
            source: "GM1",
            condensing: "X1",
            uni: 3,
            map: two,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg1_typeIV_uni3_require_exact",
            sg: 1,
            source: "GM1",
            condensing: "X1",
            uni: 3,
            map: two,
            policy: LatticePolicy::RequireExact,
            etas: both,
            expect: AnchorExpect::Refused("LatticePolicyUnsatisfied"),
        },
        Anchor {
            name: "sg3_uni10_independent_target",
            sg: 3,
            source: "GM2",
            condensing: "GM2",
            uni: 10,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: &[Some(RealSign::Minus)][..],
            expect: AnchorExpect::Success {
                target_sg: Some(1),
                target_source_label: Some("GM1"),
            },
        },
        Anchor {
            name: "sg10_uni53_repeated_constraints",
            sg: 10,
            source: "GM1+",
            condensing: "GM2+",
            uni: 53,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg12_uni56_four_antiunitary",
            sg: 12,
            source: "GM2+",
            condensing: "Y1+",
            uni: 56,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg5_uni22_raw_indices",
            sg: 5,
            source: "GM1",
            condensing: "GM2",
            uni: 22,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg11_uni59_raw_indices",
            sg: 11,
            source: "GM1+",
            condensing: "GM2-",
            uni: 59,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
        Anchor {
            name: "sg11_uni59_nontrivial_source",
            sg: 11,
            source: "GM1-",
            condensing: "GM2-",
            uni: 59,
            map: identity,
            policy: LatticePolicy::ContainmentOnly,
            etas: both,
            expect: AnchorExpect::Success {
                target_sg: None,
                target_source_label: None,
            },
        },
    ]
}

/// Re-run every pinned witness.  Returns the violations, empty when all hold.
fn run_anchors() -> Vec<String> {
    let mut violations = Vec::new();
    for anchor in anchors() {
        let rows = match isotropy::magnetic_isotropy_subgroups(
            anchor.sg,
            anchor.condensing,
            LabelConvention::Cdml,
        ) {
            Ok(rows) => rows,
            Err(error) => {
                violations.push(format!("{}: condensing row: {error}", anchor.name));
                continue;
            }
        };
        let Some(index) = rows.iter().position(|row| row.record.mag_sg == anchor.uni) else {
            violations.push(format!(
                "{}: UNI {} is not a magnetic isotropy row of SG {} {}",
                anchor.name, anchor.uni, anchor.sg, anchor.condensing
            ));
            continue;
        };
        for eta in anchor.etas {
            let request = ScalarGammaRequest::new(
                SourceIdentity {
                    parent_sg: anchor.sg,
                    label: anchor.source,
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: anchor.sg,
                    condensing_label: anchor.condensing,
                    convention: LabelConvention::Cdml,
                    direction_index: index,
                },
                *eta,
            )
            .with_embedding(EmbeddingEvidence {
                map: anchor.map,
                shift: [Rat::ZERO; 3],
                lattice_policy: anchor.policy,
            });
            let model = match verify_scalar_gamma_input(&request) {
                Ok(model) => model,
                Err(error) => {
                    let named = input_variant(&error);
                    match &anchor.expect {
                        AnchorExpect::Refused(expected) if named == *expected => continue,
                        _ => {
                            violations.push(format!(
                                "{} (eta {eta:?}): refused as {named}: {error}",
                                anchor.name
                            ));
                            continue;
                        }
                    }
                }
            };
            if let AnchorExpect::Refused(expected) = &anchor.expect {
                violations.push(format!(
                    "{} (eta {eta:?}): expected the refusal {expected}, but the contract accepted \
                     the embedding",
                    anchor.name
                ));
                continue;
            }
            let report = match subduce_scalar_gamma(&model) {
                Ok(report) => report,
                Err(error) => {
                    violations.push(format!(
                        "{} (eta {eta:?}): subduction refused as {}: {error}",
                        anchor.name,
                        subduction_variant(&error)
                    ));
                    continue;
                }
            };
            if let Some(failure) = reconstruction_failure(&model, &report) {
                violations.push(format!("{} (eta {eta:?}): {failure}", anchor.name));
            }
            if let AnchorExpect::Success {
                target_sg,
                target_source_label,
            } = &anchor.expect
            {
                if let Some(expected) = target_sg
                    && report.target_unitary_sg() != *expected
                {
                    violations.push(format!(
                        "{} (eta {eta:?}): target unitary SG {} != {expected}",
                        anchor.name,
                        report.target_unitary_sg()
                    ));
                }
                if let Some(expected) = target_source_label
                    && report.target_source_label() != *expected
                {
                    violations.push(format!(
                        "{} (eta {eta:?}): target source label {} != {expected}",
                        anchor.name,
                        report.target_source_label()
                    ));
                }
                if model.uni() != anchor.uni {
                    violations.push(format!(
                        "{} (eta {eta:?}): record UNI {} != {}",
                        anchor.name,
                        model.uni(),
                        anchor.uni
                    ));
                }
            }
        }
    }
    violations
}

/// One Γ source of dimension ≥ 2 found by the inventory scan, with the data
/// the next card would have to consume.
#[derive(Clone, Copy, Debug)]
struct BoundarySource {
    sg: u8,
    label: &'static str,
    dim: u8,
    /// Whether the row is a compound (CIR/realified) row rather than an ordinary
    /// single-valued irreducible row.
    compound: bool,
    /// Length of the tabulated matrix storage of the row.
    matrices_len: usize,
    /// Number of tabulated character values of the row.
    characters_len: usize,
    /// Naming prefix this source is reported under.
    kind: &'static str,
}

/// The explicitly named cross-check source, read from the tables like the scan
/// reads them; `None` when that parent has no such Γ row.
fn named_source() -> Option<BoundarySource> {
    let (sg, label) = NAMED_TWO_DIMENSIONAL_GAMMA;
    let irrep = query::irreps_of(sg).iter().find(|irrep| {
        irrep.ml == label
            && !irrep.spinor
            && irrep.dim >= 2
            && irrep.k_vector().numerators == GAMMA_NUMERATORS
    })?;
    Some(BoundarySource {
        sg,
        label: irrep.ml,
        dim: irrep.dim,
        compound: irrep.compound_metadata().is_some(),
        matrices_len: irrep.matrices().len(),
        characters_len: irrep.characters().len(),
        kind: "two_dimensional_gamma_named",
    })
}

/// Boundary evidence for the next card: the first Γ source of dimension ≥ 2 and
/// the first **ordinary** (non-compound) one, each really offered to the contract
/// so the refusal is the contract's own, plus the Γ corep dimensions of that
/// parent's magnetic summaries on the target side.
///
/// The two sources are reported separately because the compound row is a
/// realified/compound label and must not be read as an ordinary complex
/// irreducible of dimension 2.
fn boundary_two_dimensional_gamma() -> (Vec<BoundarySource>, Vec<String>) {
    let mut first: Option<BoundarySource> = None;
    let mut ordinary: Option<BoundarySource> = None;
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor || irrep.dim < 2 || irrep.k_vector().numerators != GAMMA_NUMERATORS {
                continue;
            }
            let candidate = BoundarySource {
                sg,
                label: irrep.ml,
                dim: irrep.dim,
                compound: irrep.compound_metadata().is_some(),
                matrices_len: irrep.matrices().len(),
                characters_len: irrep.characters().len(),
                kind: "two_dimensional_gamma",
            };
            if first.is_none_or(|current| (candidate.sg, candidate.dim) < (current.sg, current.dim))
            {
                first = Some(candidate);
            }
            if !candidate.compound
                && ordinary
                    .is_none_or(|current| (candidate.sg, candidate.dim) < (current.sg, current.dim))
            {
                ordinary = Some(BoundarySource {
                    kind: "two_dimensional_gamma_ordinary",
                    ..candidate
                });
            }
        }
    }
    let sources: Vec<BoundarySource> = [first, ordinary, named_source()]
        .into_iter()
        .flatten()
        .collect();
    let mut notes = Vec::new();
    if sources.is_empty() {
        notes.push("two_dimensional_gamma: none".to_string());
        return (sources, notes);
    }
    for source in &sources {
        let BoundarySource {
            sg,
            label,
            dim,
            compound,
            matrices_len,
            characters_len,
            kind,
        } = *source;
        let record = query::magnetic_isotropy_subgroups_of(sg)
            .into_iter()
            .next()
            .map(|entry| (entry.ml_label, entry.subgroup.mag_sg));
        let Some((condensing_label, _uni)) = record else {
            notes.push(format!(
                "{kind} sg {sg} {label} dim {dim} compound {compound}: SG {sg} has no magnetic \
                 isotropy row to name as a condensing record"
            ));
            continue;
        };
        // A real call, with a real isotropy row of the same parent as the record
        // identity, so the refusal is the contract's own and not an inference
        // from the dimension in the table.
        let request = ScalarGammaRequest::new(
            SourceIdentity {
                parent_sg: sg,
                label,
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: sg,
                condensing_label,
                convention: LabelConvention::Cdml,
                direction_index: 0,
            },
            None,
        );
        match verify_scalar_gamma_input(&request) {
            Ok(_) => notes.push(format!(
                "{kind} sg {sg} {label} dim {dim} compound {compound} matrices_len \
                 {matrices_len} characters_len {characters_len}: the contract accepted it, which \
                 contradicts the stated domain"
            )),
            Err(error) => notes.push(format!(
                "{kind} sg {sg} {label} dim {dim} compound {compound} matrices_len \
                 {matrices_len} characters_len {characters_len}: refused as {}",
                input_variant(&error)
            )),
        }
        notes.extend(boundary_gamma_target_dimensions(sg, kind));
    }
    (sources, notes)
}

/// Evidence about the **target** side of a 2D Γ case: every magnetic group of
/// that parent is read at Γ through the magnetic summary, and its corep
/// dimensions are counted, so the report can say whether higher-dimensional
/// targets exist in the shipped data or whether the card lacks both halves.
fn boundary_gamma_target_dimensions(sg: u8, kind: &'static str) -> Vec<String> {
    let mut notes = Vec::new();
    let rows = query::magnetic_isotropy_subgroups_of(sg);
    if rows.is_empty() {
        notes.push(format!(
            "{kind}_target sg {sg}: no magnetic isotropy row to read"
        ));
        return notes;
    }
    let mut unis: Vec<usize> = rows.iter().map(|entry| entry.subgroup.mag_sg).collect();
    unis.sort_unstable();
    unis.dedup();
    let mut dimensions: BTreeMap<usize, usize> = BTreeMap::new();
    let mut summary_failed = 0usize;
    let mut gamma_missing = 0usize;
    let mut with_higher_dimension = 0usize;
    let mut higher_dimension_single_non_spinor_source = 0usize;
    for uni in &unis {
        let Ok(partial) = magnetic_summary::magnetic_irrep_summary_by_uni_partial(*uni) else {
            summary_failed += 1;
            continue;
        };
        let Some(gamma) = partial
            .summary
            .kpoints
            .iter()
            .find(|kpoint| kpoint.coords == GAMMA)
        else {
            gamma_missing += 1;
            continue;
        };
        let mut higher = false;
        for corep in &gamma.coreps {
            *dimensions.entry(corep.dim).or_insert(0) += 1;
            if corep.dim >= 2 {
                higher = true;
                if corep.source_irreps.len() == 1 && !corep.source_irreps[0].spinor {
                    higher_dimension_single_non_spinor_source += 1;
                }
            }
        }
        if higher {
            with_higher_dimension += 1;
        }
    }
    notes.push(format!(
        "{kind}_target sg {sg} rows {} distinct_unis {} summary_failed {} gamma_missing {} \
         gammas_with_higher_dimension {with_higher_dimension} \
         higher_dim_single_non_spinor_source_coreps \
         {higher_dimension_single_non_spinor_source} gamma_corep_dimensions {dimensions:?}",
        rows.len(),
        unis.len(),
        summary_failed,
        gamma_missing
    ));
    notes
}

// ---------------------------------------------------------------------------
// Command line and summary
// ---------------------------------------------------------------------------

/// The requested scope: SG numbers in increasing order.
fn parse_scope(spec: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    for part in spec.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((start, end)) = part.split_once('-') {
            let start: u8 = start
                .trim()
                .parse()
                .map_err(|_| format!("bad range {part}"))?;
            let end: u8 = end
                .trim()
                .parse()
                .map_err(|_| format!("bad range {part}"))?;
            if start == 0 || end > 230 || start > end {
                return Err(format!("out-of-range {part}"));
            }
            out.extend(start..=end);
        } else {
            let sg: u8 = part.parse().map_err(|_| format!("bad SG {part}"))?;
            if sg == 0 || sg > 230 {
                return Err(format!("out-of-range SG {part}"));
            }
            out.push(sg);
        }
    }
    out.sort_unstable();
    out.dedup();
    if out.is_empty() {
        return Err("the scope is empty".to_string());
    }
    Ok(out)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut gate = false;
    let mut inventory_only = false;
    let mut scope_spec: Option<String> = None;
    let mut smoke = false;
    let mut tsv_path: Option<String> = None;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--gate" => gate = true,
            "--inventory" => inventory_only = true,
            "--smoke" => smoke = true,
            "--tsv" => {
                index += 1;
                match args.get(index) {
                    Some(path) => tsv_path = Some(path.clone()),
                    None => {
                        eprintln!("--tsv needs a path");
                        return ExitCode::from(2);
                    }
                }
            }
            "--sgs" => {
                index += 1;
                match args.get(index) {
                    Some(spec) => scope_spec = Some(spec.clone()),
                    None => {
                        eprintln!("--sgs needs a specification");
                        return ExitCode::from(2);
                    }
                }
            }
            other => {
                eprintln!("unknown argument {other}");
                return ExitCode::from(2);
            }
        }
        index += 1;
    }
    if smoke && scope_spec.is_none() {
        scope_spec = Some("1-12,221".to_string());
    }
    let scope: Vec<u8> = match scope_spec.as_deref() {
        Some(spec) => match parse_scope(spec) {
            Ok(scope) => scope,
            Err(error) => {
                eprintln!("--sgs: {error}");
                return ExitCode::from(2);
            }
        },
        None => (1..=230u8).collect(),
    };

    let started = std::time::Instant::now();
    let inventories: Vec<SgInventory> = (1..=230u8).map(inventory_of).collect();

    // Whole-corpus denominator, always printed: the probe pass may be a subset.
    let mut denominator_records = 0usize;
    let mut denominator_gamma_sources = 0usize;
    let mut denominator_eligible = 0usize;
    let mut denominator_not_gamma = 0usize;
    let mut denominator_gamma_by_label = 0usize;
    let mut denominator_records_by_query = 0usize;
    let mut denominator_pairs = 0usize;
    let mut denominator_eligible_pairs = 0usize;
    let mut exclusions: BTreeMap<String, usize> = BTreeMap::new();
    let mut dimension_histogram: BTreeMap<u8, usize> = BTreeMap::new();
    for inventory in &inventories {
        let gamma_sources = inventory.sources.iter().filter(|s| s.is_gamma).count();
        let eligible = inventory
            .sources
            .iter()
            .filter(|s| s.is_gamma && s.exclusion.is_none())
            .count();
        denominator_records += inventory.records.len();
        denominator_records_by_query += inventory.records_by_query;
        denominator_gamma_sources += gamma_sources;
        denominator_eligible += eligible;
        denominator_not_gamma += inventory.irreps_not_gamma;
        denominator_gamma_by_label += inventory.gamma_by_label;
        // The pair count is a per-parent product, not the product of the two
        // corpus totals: the eligible sources of one parent pair with that
        // parent's records only.
        denominator_pairs += inventory.records.len() * gamma_sources;
        denominator_eligible_pairs += inventory.records.len() * eligible;
        for source in &inventory.sources {
            if let Some(exclusion) = &source.exclusion {
                *exclusions.entry(exclusion.clone()).or_insert(0) += 1;
            }
            if source.exclusion.as_deref() != Some("spinor_source") {
                *dimension_histogram.entry(source.dim).or_insert(0) += 1;
            }
        }
    }
    println!("scope_sgs {}", scope.len());
    if scope.len() != 230 {
        println!(
            "scope_list {}",
            scope
                .iter()
                .map(|sg| sg.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
    }
    println!("denominator_records {denominator_records}");
    println!("denominator_records_by_query {denominator_records_by_query}");
    println!("denominator_gamma_sources {denominator_gamma_sources}");
    println!("denominator_eligible_sources {denominator_eligible}");
    println!("denominator_not_gamma_irreps {denominator_not_gamma}");
    println!("denominator_gamma_by_label {denominator_gamma_by_label}");
    println!("source_exclusions {exclusions:?}");
    println!("gamma_nonspinor_dimension_histogram {dimension_histogram:?}");
    println!("denominator_pairs {denominator_pairs}");
    println!("denominator_eligible_pairs {denominator_eligible_pairs}");

    if inventory_only {
        println!("elapsed {:.1}s", started.elapsed().as_secs_f64());
        return ExitCode::SUCCESS;
    }

    // Probe pass over the selected space groups, one task per SG.
    let reports: Vec<SgReport> = scope
        .par_iter()
        .map(|sg| census_sg(&inventories[usize::from(*sg) - 1]))
        .collect();

    let mut counters = Counters::default();
    for report in &reports {
        counters.merge(&report.counters);
    }

    // The TSV is assembled in SG order, so its bytes do not depend on threads.
    let mut tsv = String::new();
    if tsv_path.is_some() {
        tsv.push_str(
            "sg\tcondensing\tdirection_index\tuni\tsource\teta\tstage\tvariant\tmap\tshift\t\
             target\ttarget_sg\tmultiplicity\tintertwiner_re\tintertwiner_im\toperators\tdetail\n",
        );
    }
    for report in &reports {
        tsv.push_str(&report.tsv);
    }
    if let Some(path) = &tsv_path {
        if let Err(error) = std::fs::write(path, tsv.as_bytes()) {
            eprintln!("cannot write {path}: {error}");
            return ExitCode::from(2);
        }
        println!("tsv {path}");
        println!("tsv_bytes {}", tsv.len());
    }

    // Per-space-group lines, in SG order.
    for report in &reports {
        println!("{}", report.counters.line(report.sg));
    }

    // Boundary evidence for the next card: the first 2D Γ source and the first
    // ordinary (non-compound) one.
    let (boundary_sources, boundary_notes) = boundary_two_dimensional_gamma();
    for note in &boundary_notes {
        println!("boundary {note}");
    }

    // Pinned witnesses.
    let anchor_violations = run_anchors();

    println!("records_total {}", counters.records_total);
    println!("records_with_fixed_embedding {}", counters.records_fixed);
    println!(
        "records_fixed_with_identity_map {}",
        counters.records_fixed_identity_map
    );
    println!(
        "records_fixed_with_identity_zero_shift {}",
        counters.records_fixed_identity_zero_shift
    );
    println!(
        "records_without_fixed_embedding {}",
        counters.records_no_fixed
    );
    println!(
        "records_no_realising_candidate {}",
        counters.records_no_realising
    );
    println!(
        "records_certified_absence {}",
        counters.records_certified_absence
    );
    println!(
        "records_undecided_residual {}",
        counters.records_undecided_residual
    );
    println!("records_record_error {}", counters.records_error);
    println!("sources_gamma {}", counters.sources_gamma);
    println!("sources_eligible {}", counters.sources_eligible);
    println!(
        "sources_excluded_total {}",
        counters.sources_gamma - counters.sources_eligible
    );
    println!("pairs_total {}", counters.pairs_total);
    println!("pairs_probed {}", counters.pairs_probed);
    println!(
        "pairs_record_without_fixed {}",
        counters.pairs_record_without_fixed
    );
    println!(
        "pairs_source_domain_excluded {}",
        counters.pairs_source_domain_excluded
    );
    println!("probes {}", counters.probes);
    println!("success {}", counters.success);
    println!("verify_input_failures {}", counters.verify_input_failures);
    println!("subduce_failures {}", counters.subduce_failures);
    println!("success_check_failures {}", counters.success_check_failures);
    println!("tsv_rows {}", counters.tsv_rows);
    println!(
        "tsv_rows_record_without_fixed {}",
        counters.tsv_rows_record_without_fixed
    );
    println!("eta_status_mismatches {}", counters.eta_status_mismatches);
    println!(
        "eta_target_label_differs {}",
        counters.eta_target_label_differs
    );
    println!(
        "eta_operator_table_differs {}",
        counters.eta_operator_table_differs
    );
    println!(
        "eta_intertwiner_differs {}",
        counters.eta_intertwiner_differs
    );
    println!("verify_variants {:?}", counters.verify_variants);
    println!("subduce_variants {:?}", counters.subduce_variants);
    println!("no_fixed_variants {:?}", counters.no_fixed_variants);
    println!("exclusion_variants {:?}", counters.exclusion_variants);
    println!("anchor_violations {}", anchor_violations.len());
    println!("elapsed {:.1}s", started.elapsed().as_secs_f64());

    if gate {
        let mut violations = anchor_violations.clone();
        // The gate reads the rows that were actually emitted, per space group and
        // as one document: dropped rows, duplicated η, unknown or missing keys,
        // a record whose fixed (map, shift) is not one pair, and stage counts that
        // disagree with the counters all show up here.
        for report in &reports {
            violations.extend(validate_tsv(
                &inventories[usize::from(report.sg) - 1],
                &report.tsv,
                &report.counters,
            ));
        }
        violations.extend(validate_global_tsv(&tsv, &counters));
        // Coverage-independent traversal completeness.
        let scoped_records: usize = scope
            .iter()
            .map(|sg| inventories[usize::from(*sg) - 1].records_by_query)
            .sum();
        if counters.records_total != scoped_records {
            violations.push(format!(
                "traversal: {} records enumerated, {scoped_records} reached through \
                 magnetic_isotropy_subgroups_of",
                counters.records_total
            ));
        }
        // Population partition.  The pair identities are per parent, because
        // the pair count is a sum of per-parent products; checking the two
        // aggregate products against each other would be a different (and
        // false) identity.
        for report in &reports {
            let c = &report.counters;
            let sg = report.sg;
            if c.records_fixed + c.records_no_fixed != c.records_total {
                violations.push(format!(
                    "SG {sg}: records do not partition into fixed / without fixed"
                ));
            }
            if c.records_no_realising
                + c.records_certified_absence
                + c.records_undecided_residual
                + c.records_error
                != c.records_no_fixed
            {
                violations.push(format!(
                    "SG {sg}: the records without a fixed embedding do not partition by reason"
                ));
            }
            if c.sources_spinor
                + c.sources_dimension_gt_one
                + c.sources_character_unavailable
                + c.sources_non_real_character
                + c.sources_character_off_grid
                != c.sources_gamma - c.sources_eligible
            {
                violations.push(format!("SG {sg}: the source exclusions do not add up"));
            }
            if c.pairs_source_domain_excluded
                != c.records_total * (c.sources_gamma - c.sources_eligible)
            {
                violations.push(format!(
                    "SG {sg}: the domain-excluded pair count is not records x excluded"
                ));
            }
            if c.pairs_probed + c.pairs_record_without_fixed
                != c.pairs_total - c.pairs_source_domain_excluded
            {
                violations.push(format!("SG {sg}: the pair partition does not close"));
            }
            if c.records_total * c.sources_eligible != c.pairs_probed + c.pairs_record_without_fixed
            {
                violations.push(format!(
                    "SG {sg}: eligible pairs are not records x eligible"
                ));
            }
            if c.probes != c.pairs_probed * 2 {
                violations.push(format!(
                    "SG {sg}: each probed pair must carry both eta probes"
                ));
            }
            if c.probes
                != c.success
                    + c.verify_input_failures
                    + c.subduce_failures
                    + c.success_check_failures
            {
                violations.push(format!("SG {sg}: the probes do not partition into stages"));
            }
            if c.tsv_rows != (c.pairs_probed + c.pairs_record_without_fixed) * 2 {
                violations.push(format!(
                    "SG {sg}: the TSV row count does not match the pairs"
                ));
            }
        }
        // The aggregate probe identities are sums of the per-parent ones; they
        // are checked again so the printed totals cannot drift apart.
        if counters.probes != counters.pairs_probed * 2 {
            violations.push("each probed pair must carry both eta probes".to_string());
        }
        if counters.probes
            != counters.success
                + counters.verify_input_failures
                + counters.subduce_failures
                + counters.success_check_failures
        {
            violations.push("the probes do not partition into stages".to_string());
        }
        if counters.tsv_rows != (counters.pairs_probed + counters.pairs_record_without_fixed) * 2 {
            violations.push("the TSV row count does not match the pairs".to_string());
        }
        if counters.pairs_probed + counters.pairs_record_without_fixed
            != counters.pairs_total - counters.pairs_source_domain_excluded
        {
            violations.push("the aggregate pair partition does not close".to_string());
        }
        if counters.success_check_failures != 0 {
            violations.push(format!(
                "successful probes failed the reconstruction: {}",
                counters.success_check_failures
            ));
        }
        if counters.eta_status_mismatches != 0 {
            violations.push(format!(
                "the two eta probes disagree on status for {} pairs",
                counters.eta_status_mismatches
            ));
        }
        for name in counters.no_fixed_variants.keys() {
            if ![
                "no_realising_candidate",
                "certified_absence",
                "undecided_residual",
                "geometry",
                "operations",
                "verdict",
            ]
            .contains(&name.as_str())
            {
                violations.push(format!("unnamed no-fixed-embedding reason {name}"));
            }
        }
        if boundary_sources.is_empty() {
            violations.push("no two-dimensional Gamma source was found".to_string());
        }
        for source in &boundary_sources {
            let prefix = format!("{} sg {} {}", source.kind, source.sg, source.label);
            if !boundary_notes.iter().any(|note| {
                note.starts_with(&prefix) && note.contains("refused as NotScalarSource")
            }) {
                violations.push(format!(
                    "the {prefix} boundary source was not refused as NotScalarSource"
                ));
            }
        }
        if !violations.is_empty() {
            for violation in &violations {
                eprintln!("magnetic_gamma_subduction_census: {violation}");
            }
            return ExitCode::from(1);
        }
        println!("VERDICT complete scope=restricted-scalar-gamma-subduction-census gate=--gate");
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use cryspglib::irrep::magnetic_embedding::EmbeddingKind;

    fn rational(numerator: i128, denominator: i128) -> Rat {
        Rat::new(numerator, denominator).expect("fraction")
    }

    /// One record, one eligible Γ source, and the counters a single probed pair
    /// must produce.  No database is touched.
    fn fixture_inventory() -> SgInventory {
        SgInventory {
            sg: 1,
            records: vec![RecordRow {
                condensing: "GM1",
                direction_index: 0,
                uni: 7,
            }],
            sources: vec![SourceRow {
                label: "GM1",
                dim: 1,
                is_gamma: true,
                exclusion: None,
            }],
            irreps_not_gamma: 0,
            gamma_by_label: 1,
            gamma_by_data: 1,
            records_by_query: 1,
        }
    }

    /// The counters of one probed pair (both η succeed).
    fn fixture_probed_counters() -> Counters {
        Counters {
            records_total: 1,
            records_fixed: 1,
            pairs_total: 1,
            pairs_probed: 1,
            probes: 2,
            success: 2,
            tsv_rows: 2,
            ..Counters::default()
        }
    }

    /// The counters of one pair whose record has no fixed embedding.
    fn fixture_unprobed_counters() -> Counters {
        Counters {
            records_total: 1,
            records_no_fixed: 1,
            records_certified_absence: 1,
            pairs_total: 1,
            pairs_record_without_fixed: 1,
            tsv_rows: 2,
            tsv_rows_record_without_fixed: 2,
            ..Counters::default()
        }
    }

    /// A successful fixture outcome, formatted exactly like a real one.
    fn fixture_outcome() -> Probe {
        Probe::Success {
            target_label: "GM1".to_string(),
            target_sg: 1,
            multiplicity: 1,
            operators: Vec::new(),
            intertwiner: Complex64::new(1.0, 0.0),
        }
    }

    /// Two probed rows of the fixture record.
    fn fixture_probed_tsv() -> String {
        let record = RecordRow {
            condensing: "GM1",
            direction_index: 0,
            uni: 7,
        };
        let mut tsv = String::new();
        for eta in [RealSign::Plus, RealSign::Minus] {
            push_row(
                &mut tsv,
                1,
                &record,
                "GM1",
                eta,
                &fixture_outcome(),
                Mat3R::identity(),
                [Rat::ZERO; 3],
            );
        }
        tsv
    }

    /// Two rows of the fixture record without a fixed embedding.
    fn fixture_unprobed_tsv() -> String {
        let record = RecordRow {
            condensing: "GM1",
            direction_index: 0,
            uni: 7,
        };
        let mut tsv = String::new();
        for eta in [RealSign::Plus, RealSign::Minus] {
            push_no_fixed_row(
                &mut tsv,
                1,
                &record,
                "GM1",
                eta,
                "certified_absence",
                "every candidate map of the R24 family was refuted",
            );
        }
        tsv
    }

    /// The row validator must accept the normal two-η output and reject a
    /// dropped row, a same-length duplicated η, and unknown identities.
    #[test]
    fn the_tsv_validator_rejects_dropped_duplicated_and_unknown_rows() {
        let inventory = fixture_inventory();
        let counters = fixture_probed_counters();
        let tsv = fixture_probed_tsv();
        assert!(
            validate_tsv(&inventory, &tsv, &counters).is_empty(),
            "{:?}",
            validate_tsv(&inventory, &tsv, &counters)
        );

        // A dropped row: the same counters, one line fewer.
        let dropped: String = tsv
            .lines()
            .take(1)
            .map(|line| format!("{line}\n"))
            .collect();
        let violations = validate_tsv(&inventory, &dropped, &counters);
        assert!(
            violations.iter().any(|v| v.contains("missing row")),
            "a dropped row must be a missing key: {violations:?}"
        );

        // The same number of rows, but one η written twice.
        let duplicated = tsv.replace("\tGM1\t-1\t", "\tGM1\t+1\t");
        assert_eq!(
            duplicated.lines().count(),
            tsv.lines().count(),
            "the mutation must keep the row count"
        );
        let violations = validate_tsv(&inventory, &duplicated, &counters);
        assert!(
            violations.iter().any(|v| v.contains("duplicated row")),
            "a duplicated η must be caught: {violations:?}"
        );
        assert!(
            violations.iter().any(|v| v.contains("missing row")),
            "the lost η must be a missing key: {violations:?}"
        );

        // An unknown record identity, at the same row count.
        let unknown_record = tsv.replace("\t7\t", "\t999\t");
        let violations = validate_tsv(&inventory, &unknown_record, &counters);
        assert!(
            violations.iter().any(|v| v.contains("unknown row")),
            "an unknown record must be caught: {violations:?}"
        );

        // An unknown source identity.
        let unknown_source = tsv.replace("\tGM1\t+1\t", "\tX1\t+1\t");
        let violations = validate_tsv(&inventory, &unknown_source, &counters);
        assert!(
            violations.iter().any(|v| v.contains("unknown row")),
            "an unknown source must be caught: {violations:?}"
        );

        // One η carries a different (map, shift) than the other: the record's
        // fixed embedding is no longer one pair.
        let inconsistent = tsv.replacen(
            "1/1,0/1,0/1;0/1,1/1,0/1;0/1,0/1,1/1",
            "2/1,0/1,0/1;0/1,1/1,0/1;0/1,0/1,1/1",
            1,
        );
        let violations = validate_tsv(&inventory, &inconsistent, &counters);
        assert!(
            violations
                .iter()
                .any(|v| v.contains("distinct (map, shift)")),
            "an inconsistent fixed embedding must be caught: {violations:?}"
        );

        // A stage count that the text does not support.
        let mut wrong_counters = fixture_probed_counters();
        wrong_counters.success = 1;
        let violations = validate_tsv(&inventory, &tsv, &wrong_counters);
        assert!(
            violations.iter().any(|v| v.contains("rows staged success")),
            "a counter must not survive the text: {violations:?}"
        );
    }

    /// The same validator accepts an unprobed record's rows and rejects an
    /// unknown reason.
    #[test]
    fn the_tsv_validator_covers_the_unprobed_rows() {
        let inventory = fixture_inventory();
        let counters = fixture_unprobed_counters();
        let tsv = fixture_unprobed_tsv();
        assert!(
            validate_tsv(&inventory, &tsv, &counters).is_empty(),
            "{:?}",
            validate_tsv(&inventory, &tsv, &counters)
        );
        let mutated = tsv.replace("certified_absence", "made_up_reason");
        let violations = validate_tsv(&inventory, &mutated, &counters);
        assert!(
            violations
                .iter()
                .any(|v| v.contains("no-fixed-embedding reason")),
            "an unknown reason must be caught: {violations:?}"
        );
    }

    /// The two reconstruction predicates are positive but not the same shape:
    /// `is_finite_unit` rejects NaN, infinities **and zero** (the intertwiner
    /// must be a non-zero unit), while `within_tolerance` accepts a legitimate
    /// zero residual and rejects only non-finite or too-large values.
    #[test]
    fn the_finite_predicates_reject_nan_and_non_units() {
        assert!(is_finite_unit(1.0));
        assert!(is_finite_unit(1.0 + 1e-13));
        assert!(!is_finite_unit(f64::NAN));
        assert!(!is_finite_unit(f64::INFINITY));
        assert!(!is_finite_unit(f64::NEG_INFINITY));
        assert!(!is_finite_unit(0.0));
        assert!(!is_finite_unit(-1.0));
        assert!(!is_finite_unit(1.0 + 1e-6));
        assert!(within_tolerance(0.0));
        assert!(within_tolerance(1e-10));
        assert!(!within_tolerance(f64::NAN));
        assert!(!within_tolerance(f64::INFINITY));
        assert!(!within_tolerance(1e-6));
    }

    /// The fixed embedding must not depend on the order the family produced its
    /// entries in: the key is total and the first entry under it is the choice.
    #[test]
    fn the_fixed_candidate_key_is_total_and_order_independent() {
        let entry = |map: Mat3R, shift: [Rat; 3], kind: EmbeddingKind| ParentEmbedding {
            map,
            shift,
            kind,
            realises_record_lattice: true,
        };
        let identity = Mat3R::identity();
        let negated = Mat3R::from_ints([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]);
        let two = Mat3R::from_ints([[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
        let a = entry(identity, [Rat::ZERO; 3], EmbeddingKind::Plain);
        let b = entry(
            two,
            [rational(1, 2), Rat::ZERO, Rat::ZERO],
            EmbeddingKind::PrimitiveBasis,
        );
        let c = entry(
            identity,
            [rational(1, 3), Rat::ZERO, Rat::ZERO],
            EmbeddingKind::ConventionalBasis,
        );
        let mut forward = vec![&a, &b, &c];
        let mut backward = vec![&c, &b, &a];
        forward.sort_by_key(|entry| embedding_key(entry));
        backward.sort_by_key(|entry| embedding_key(entry));
        let forward: Vec<_> = forward.iter().map(|entry| embedding_key(entry)).collect();
        let backward: Vec<_> = backward.iter().map(|entry| embedding_key(entry)).collect();
        assert_eq!(forward, backward, "the sort is not total");
        // The tabulated frame wins when it is among the candidates, and among
        // the identity maps the smallest solved shift does.
        assert!(!forward[0].0, "the identity map must be preferred");
        assert_eq!(forward[0].2, [(0, 1), (0, 1), (0, 1)]);
        assert_eq!(forward[0].1[0], (1, 1));
        assert_eq!(forward[1].2, [(1, 3), (0, 1), (0, 1)]);
        // Without an identity candidate the exact entries order the rest, and a
        // negated map sorts before the doubled one.
        let d = entry(negated, [Rat::ZERO; 3], EmbeddingKind::Plain);
        let e = entry(two, [Rat::ZERO; 3], EmbeddingKind::Plain);
        let mut without_identity = vec![&e, &d];
        without_identity.sort_by_key(|entry| embedding_key(entry));
        let keys: Vec<_> = without_identity
            .iter()
            .map(|entry| embedding_key(entry))
            .collect();
        assert!(keys[0].0 && keys[1].0);
        assert_eq!(keys[0].1[0], (-1, 1));
        assert_eq!(keys[1].1[0], (2, 1));
    }

    /// The domain screen and the contract must read the same exact grid.
    #[test]
    fn the_domain_screen_reads_the_contract_grid() {
        assert_eq!(exact_fraction(0.0), Some(Rat::ZERO));
        assert_eq!(exact_fraction(1.0), Some(Rat::ONE));
        assert_eq!(exact_fraction(-1.0), Some(Rat::ONE.checked_neg().unwrap()));
        assert_eq!(
            exact_fraction(-1.0 + 1e-10),
            Some(Rat::ONE.checked_neg().unwrap())
        );
        assert_eq!(exact_fraction(0.5), Some(rational(1, 2)));
        assert_eq!(exact_fraction(1.0 / 3.0), Some(rational(1, 3)));
        assert_eq!(exact_fraction(0.05), None, "1/20 is off the exact grid");
        assert_eq!(exact_fraction(f64::NAN), None);
        let plus = Complex64::new(1.0, 0.0);
        let minus = Complex64::new(-1.0, 0.0);
        let complex = Complex64::new(0.0, 1.0);
        let fractional = Complex64::new(1.0 / 3.0, 0.0);
        assert!(is_real_plus_minus(plus));
        assert!(is_real_plus_minus(minus));
        assert!(!is_real_plus_minus(complex));
        assert!(!is_real_plus_minus(fractional));
    }

    /// The screen's `±1` predicate, shared with the test above.
    fn is_real_plus_minus(value: Complex64) -> bool {
        let (Some(real), Some(imaginary)) = (exact_fraction(value.re), exact_fraction(value.im))
        else {
            return false;
        };
        let negative_one = Rat::ONE.checked_neg().expect("one is representable");
        imaginary.is_zero() && (real == Rat::ONE || real == negative_one)
    }

    /// The record inventory must agree with the independent query and must add
    /// up to the corpus the R24 census pinned.
    #[test]
    fn the_record_inventory_matches_the_independent_query() {
        let mut total = 0usize;
        let mut by_query = 0usize;
        for sg in 1..=230u8 {
            let inventory = inventory_of(sg);
            assert_eq!(
                inventory.records.len(),
                inventory.records_by_query,
                "SG {sg}: {} vs {}",
                inventory.records.len(),
                inventory.records_by_query
            );
            total += inventory.records.len();
            by_query += inventory.records_by_query;
        }
        assert_eq!(total, 16_721, "the corpus record count moved");
        assert_eq!(by_query, 16_721);
    }

    /// Every Γ irrep found through the k-vector is a source candidate, and the
    /// per-source exclusion reasons cover the reported totals.
    #[test]
    fn the_source_inventory_partitions_by_reason() {
        let mut gamma = 0usize;
        let mut eligible = 0usize;
        let mut excluded = BTreeMap::new();
        for sg in 1..=230u8 {
            let inventory = inventory_of(sg);
            for source in &inventory.sources {
                assert!(source.is_gamma, "a non-Γ irrep reached the Γ source list");
                gamma += 1;
                match &source.exclusion {
                    None => eligible += 1,
                    Some(reason) => *excluded.entry(reason.clone()).or_insert(0) += 1,
                }
            }
        }
        assert_eq!(gamma, eligible + excluded.values().sum::<usize>());
        assert!(eligible > 0, "no Γ source is inside the family");
        println!("eligible {eligible} exclusions {excluded:?}");
    }

    /// The pinned witnesses must keep their outcomes through the same entry
    /// points, including the named `RequireExact` refusal.
    #[test]
    fn the_pinned_anchors_hold() {
        let violations = run_anchors();
        assert!(violations.is_empty(), "{violations:#?}");
    }

    /// The two-dimensional Γ boundary exists and is refused by name, once for
    /// the first source and once for the first ordinary (non-compound) one.
    #[test]
    fn the_two_dimensional_gamma_boundary_is_named() {
        let (sources, notes) = boundary_two_dimensional_gamma();
        assert!(!sources.is_empty(), "{notes:?}");
        assert!(
            sources.iter().any(|source| source.compound),
            "the first 2D source is expected to be a compound row: {sources:?}"
        );
        assert!(
            sources.iter().any(|source| !source.compound),
            "an ordinary 2D source is expected: {sources:?}"
        );
        let named = named_source().expect("the named cross-check source exists");
        assert_eq!(named.dim, 2);
        assert!(!named.compound, "{named:?}");
        assert!(
            sources
                .iter()
                .any(|source| source.kind == "two_dimensional_gamma_named"),
            "{sources:?}"
        );
        for source in &sources {
            let prefix = format!("{} sg {} {}", source.kind, source.sg, source.label);
            assert!(
                notes
                    .iter()
                    .any(|note| note.starts_with(&prefix)
                        && note.contains("refused as NotScalarSource")),
                "{prefix} not refused by name: {notes:?}"
            );
        }
    }
}
