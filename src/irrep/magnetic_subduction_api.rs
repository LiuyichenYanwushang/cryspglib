//! The supported high-level entry point for Γ magnetic subduction (R30).
//!
//! The R25–R29 modules are the verified layers of the Γ magnetic pipeline:
//! [`crate::irrep::magnetic_input`] turns one explicit embedding into a verified
//! parent action, [`crate::irrep::magnetic_subduction`] identifies a
//! one-dimensional magnetic co-representation, and
//! [`crate::irrep::magnetic_matrix_subduction`] does the same for a genuine
//! two-dimensional real source.  Each layer resolves its own two labels, but a
//! caller still has to find the magnetic isotropy row's **local index**, choose
//! the scalar or matrix engine, and collect the label/provenance surface by
//! hand.  There is no public magnetic counterpart to `irrep::subduction_api`.
//!
//! This module is only the public request/answer spelling of that pipeline.  It
//! adds no mathematics:
//!
//! * one request names the source irrep (with its own [`LabelConvention`]), the
//!   condensing record (with its own convention) and either an existing
//!   [`IsotropyDirection`] or the magnetic group's UNI number;
//! * the request is resolved against the generated magnetic isotropy table, the
//!   source dimension picks the scalar (1) or matrix (2) entry point, and the
//!   engine below is called **unchanged**;
//! * the answer carries the verified model (owned, not borrowed through a
//!   self-referential struct) plus an owned copy of the result payload, so the
//!   intertwiners, target operators and matrix copies stay usable for a per-class
//!   complex-vector reconstruction.
//!
//! # Identity and label semantics
//!
//! Four identities are deliberately kept apart, because three different tables
//! contribute them; only the first two are ordinary irrep identities:
//!
//! | field | what it is | source |
//! |---|---|---|
//! | [`MagneticSourceIdentity`] | the parent ordinary irrep that is subduced | `query::find_irreps`, CDML **and** BC both from [`IrrepRecord::labels`] |
//! | [`MagneticCondensingIdentity`] | the parent ordinary irrep whose magnetic isotropy row was selected | `isotropy::magnetic_isotropy_subgroups`, CDML **and** BC from the row's own irrep |
//! | [`MagneticRecordIdentity`] | the magnetic group and setting: UNI, BNS, ISO label, direction label, basis, origin | the generated magnetic isotropy row |
//! | [`MagneticTargetIdentity`] | the magnetic Γ corepresentation of `M` selected by the engine | magnetic corep summary at Γ, plus the target's `H` ordinary source irrep |
//!
//! The source and condensing identities are ordinary-irrep labels and therefore
//! have genuine CDML and Bradley–Cracknell spellings; the record identity is a
//! magnetic-table identity (UNI/BNS/ISO/basis/origin), and the target identity
//! is a magnetic corep identity together with the ordinary `H` irrep it comes
//! from.  The magnetic corepresentation label
//! ([`MagneticTargetIdentity::magnetic_label`]) is the **display label the
//! magnetic summary constructs** for the corep; the summary has no separate BC
//! field for a corep label, so [`MagneticTargetIdentity::magnetic_bc`] is `None`
//! and is never filled from the `H` irrep's BC spelling.  The target's ordinary
//! `H` irrep is reported independently in
//! [`MagneticTargetIdentity::h_source`] with its own verified CDML/BC labels.
//! A one-dimensional corep of a single-source `H` irrep may carry the same
//! string as that irrep's CDML label; that is data, not a licence to reuse the
//! label as a magnetic BC name, and the two fields stay separate.
//!
//! Selector boundary: the supported magnetic selectors are
//! [`IsotropyDirection::Label`], [`IsotropyDirection::Index`] and the UNI
//! number.  The magnetic tables store no official order-parameter component
//! string, so [`IsotropyDirection::OfficialDescriptor`] is *not* invented here:
//! it is passed through and surfaces the existing
//! [`IsotropyError::NoOfficialDescriptor`].  [`IsotropyDirection::Descriptor`]
//! keeps the layer's actual behaviour: its matching field for a magnetic row is
//! the row's `direction` field, which holds the ISOTROPY label, so a
//! `Descriptor` carrying that same label follows the existing normalization
//! and can select the row.  A genuine component string (`"(a,0,0)"`)
//! has no magnetic source and surfaces the existing
//! [`IsotropyError::DirectionNotFound`]; neither case is newly restricted or
//! reinterpreted here.
//!
//! # Ownership
//!
//! [`MagneticSubductionReport`] owns the verified model and an owned snapshot of
//! the result data.  It does **not** own a `model` field and a result that
//! borrows that field at the same time: the result payload (labels, intertwiners,
//! operators, matrix copies) is copied out of the engine's borrow, and every
//! remaining piece of provenance is available through the report's own getters.
//! [`MagneticSubductionReport::model`] exposes the exact scalar or matrix model
//! the engine verified, so a caller can re-run its own checks on the source
//! action, classes and translation lattice.
//!
//! # Scope
//!
//! Exactly the R25/R26/R29 input domain: non-spinor, Γ, real `±1`
//! one-dimensional or real non-compound two-dimensional parent irreps; targets
//! are Wigner type-A, single-source, one- or two-dimensional Γ
//! corepresentations.  There is no non-Γ, full-star, spinor, type-B/C,
//! higher-dimensional or setting-search path, and no ordinary-API change.
//! Errors from the layers below are returned with their original source, never
//! replaced by an identity-only answer or by a guessed row.
//!
//! R32 adds a **separate** restricted paired type-C entry point,
//! [`subduce_type_c`], with its own [`MagneticTypeCReport`],
//! [`MagneticTypeCTargetIdentity`] (two `h_sources`), [`MagneticTypeCOutcome`]
//! and [`MagneticTypeCError`].  The legacy [`MagneticTargetIdentity`],
//! [`MagneticOutcome`], [`MagneticSubductionReport`] and [`subduce`] keep their
//! exact fields, variants and type-A-only matrix domain: a type-C-only record is
//! still refused by [`subduce`] with its existing [`MagneticSubductionApiError`],
//! so no caller's exhaustive match is broken and no second H identity is
//! silently dropped.  [`subduce_type_c`] reuses the same private request resolver
//! and the same [`matrix::verify_matrix_gamma_input`] model; only its target
//! family is the two-inequivalent-real-`±1`-H type-C slice.
//!
//! ```
//! use cryspglib::irrep::magnetic_input::{EmbeddingEvidence, RealSign, SourceIdentity};
//! use cryspglib::irrep::magnetic_subduction_api::{
//!     MagneticRecordSelector, MagneticSubductionRequest, subduce,
//! };
//! use cryspglib::irrep::{LabelConvention, isotropy::IsotropyDirection};
//!
//! let request = MagneticSubductionRequest::new(
//!     SourceIdentity {
//!         parent_sg: 5,
//!         label: "GM1",
//!         convention: LabelConvention::Cdml,
//!     },
//!     MagneticRecordSelector::new(5, "GM2", LabelConvention::Cdml)
//!         .direction(IsotropyDirection::Label("P1"))
//!         .uni(22),
//!     Some(RealSign::Plus),
//! )
//! .with_embedding(EmbeddingEvidence::identity());
//! let report = subduce(&request).expect("the pinned SG 5 GM1 / GM2 witness");
//! assert_eq!(report.record().uni, 22);
//! assert_eq!(report.source().labels.cdml, "GM1");
//! assert_eq!(report.condensing().labels.cdml, "GM2");
//! let scalar = report.scalar().expect("a one-dimensional source");
//! assert_eq!(scalar.target.h_source.cdml, "GM1");
//! assert!(scalar.target.magnetic_bc.is_none());
//! assert!(!scalar.operators.is_empty());
//! ```

use num_complex::Complex64;

use crate::irrep::isotropy::{self, IsotropyDirection, IsotropyError, MagneticIsotropySubgroup};
use crate::irrep::labels::{IrrepLabels, LabelConvention};
use crate::irrep::magnetic_embedding::MagneticOperation;
use crate::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, MagneticInputError, RealSign, RecordIdentity,
    ScalarGammaModel, ScalarGammaRequest, SourceIdentity, TimeParity, verify_scalar_gamma_input,
};
use crate::irrep::magnetic_matrix_subduction::{
    self as matrix, MatrixGammaModel, MatrixGammaRequest, MatrixSubductionError, MatrixTypeCError,
    TargetMatrixOperator,
};
use crate::irrep::magnetic_subduction::{self as scalar, MagneticSubductionError, TargetOperator};
use crate::irrep::query;
use crate::irrep::subduction::{Mat3R, Rat};
use crate::irrep::types::IrrepRecord;

/// Why a [`MagneticSubductionRequest`] could not be answered.
///
/// The variants split into request resolution (labels, rows, direction), and
/// the layered errors of the verified engine.  Every engine error keeps its
/// original value and [`std::error::Error::source`], so a caller can tell a
/// missing label from a failed group law.
#[derive(Debug)]
pub enum MagneticSubductionApiError {
    /// A parent space group number is outside 1–230.
    UnsupportedSpaceGroup { sg: u8 },
    /// The source and the condensing record name different parents.
    ParentMismatch { source_sg: u8, record_sg: u8 },
    /// No irrep of the parent carries the source label under its convention.
    SourceLabelNotFound {
        sg: u8,
        label: String,
        convention: LabelConvention,
    },
    /// Several parent irreps carry the source label under that convention.
    AmbiguousSourceLabel {
        sg: u8,
        label: String,
        convention: LabelConvention,
        candidates: Vec<&'static str>,
    },
    /// The source irrep exists but is neither the one- nor the two-dimensional
    /// real Γ source of the two supported entry points.
    SourceDimensionUnsupported {
        sg: u8,
        label: &'static str,
        dimension: u8,
    },
    /// The engine named an ordinary `H` irrep for the target, but that identity
    /// does not resolve to exactly one database row, so its CDML/BC labels
    /// cannot be reported.
    HSourceIdentityUnresolved {
        sg: u8,
        label: &'static str,
        matches: usize,
    },
    /// Neither a direction nor a UNI number was given, so the record cannot be
    /// selected.
    RecordSelectorMissing { sg: u8 },
    /// No magnetic isotropy row of the condensing irrep has that UNI number.
    RecordUniNotFound {
        sg: u8,
        label: String,
        convention: LabelConvention,
        uni: usize,
    },
    /// Several magnetic isotropy rows of the condensing irrep have that UNI
    /// number; the direction must be given as well.
    AmbiguousRecordUni {
        sg: u8,
        uni: usize,
        candidates: Vec<(usize, &'static str)>,
    },
    /// An explicit direction and a UNI number were both given and they select
    /// different rows.
    RecordUniMismatch {
        sg: u8,
        direction_index: usize,
        requested: usize,
        found: usize,
    },
    /// The magnetic isotropy lookup failed (unknown/ambiguous irrep label,
    /// unknown direction, unsupported official descriptor, index out of range).
    Isotropy(IsotropyError),
    /// The R25 scalar input contract refused the request.
    ScalarInput(Box<MagneticInputError>),
    /// The R26 scalar subduction refused the verified model.
    Scalar(Box<MagneticSubductionError>),
    /// The R29 matrix input contract or matrix subduction refused the request.
    Matrix(Box<MatrixSubductionError>),
}

impl std::fmt::Display for MagneticSubductionApiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSpaceGroup { sg } => {
                write!(f, "parent space group {sg} is outside 1-230")
            }
            Self::ParentMismatch {
                source_sg,
                record_sg,
            } => write!(
                f,
                "the source's parent is SG {source_sg} but the record's is SG {record_sg}"
            ),
            Self::SourceLabelNotFound {
                sg,
                label,
                convention,
            } => write!(
                f,
                "SG {sg} has no irrep labelled {label:?} under {convention:?}"
            ),
            Self::AmbiguousSourceLabel {
                sg,
                label,
                convention,
                candidates,
            } => write!(
                f,
                "SG {sg}: {convention:?} label {label:?} matches {} irreps ({}), so the source \
                 dimension is not defined",
                candidates.len(),
                candidates.join(", ")
            ),
            Self::SourceDimensionUnsupported {
                sg,
                label,
                dimension,
            } => write!(
                f,
                "SG {sg} source {label:?} has dimension {dimension}; this entry point covers the \
                 one-dimensional scalar and two-dimensional real matrix Γ sources only"
            ),
            Self::HSourceIdentityUnresolved { sg, label, matches } => write!(
                f,
                "SG {sg}: the target's H source {label:?} resolves to {matches} database rows;                  its CDML/BC labels cannot be reported"
            ),
            Self::RecordSelectorMissing { sg } => write!(
                f,
                "SG {sg}: the condensing record needs a direction or a UNI number"
            ),
            Self::RecordUniNotFound {
                sg,
                label,
                convention,
                uni,
            } => write!(
                f,
                "SG {sg} {convention:?} {label:?} has no magnetic isotropy row with UNI {uni}"
            ),
            Self::AmbiguousRecordUni {
                sg,
                uni,
                candidates,
            } => write!(
                f,
                "SG {sg}: UNI {uni} appears in {} magnetic isotropy rows ({}); select one with a \
                 direction",
                candidates.len(),
                candidates
                    .iter()
                    .map(|(index, direction)| format!("{direction} (index {index})"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            Self::RecordUniMismatch {
                sg,
                direction_index,
                requested,
                found,
            } => write!(
                f,
                "SG {sg}: direction index {direction_index} selects UNI {found}, not the requested \
                 UNI {requested}"
            ),
            Self::Isotropy(error) => write!(f, "{error}"),
            Self::ScalarInput(error) => write!(f, "{error}"),
            Self::Scalar(error) => write!(f, "{error}"),
            Self::Matrix(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for MagneticSubductionApiError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Isotropy(error) => Some(error),
            Self::ScalarInput(error) => Some(error),
            Self::Scalar(error) => Some(error),
            Self::Matrix(error) => Some(error),
            Self::UnsupportedSpaceGroup { .. }
            | Self::ParentMismatch { .. }
            | Self::SourceLabelNotFound { .. }
            | Self::AmbiguousSourceLabel { .. }
            | Self::SourceDimensionUnsupported { .. }
            | Self::HSourceIdentityUnresolved { .. }
            | Self::RecordSelectorMissing { .. }
            | Self::RecordUniNotFound { .. }
            | Self::AmbiguousRecordUni { .. }
            | Self::RecordUniMismatch { .. } => None,
        }
    }
}

impl From<IsotropyError> for MagneticSubductionApiError {
    fn from(error: IsotropyError) -> Self {
        Self::Isotropy(error)
    }
}

impl From<MagneticInputError> for MagneticSubductionApiError {
    fn from(error: MagneticInputError) -> Self {
        Self::ScalarInput(Box::new(error))
    }
}

impl From<MagneticSubductionError> for MagneticSubductionApiError {
    fn from(error: MagneticSubductionError) -> Self {
        Self::Scalar(Box::new(error))
    }
}

impl From<MatrixSubductionError> for MagneticSubductionApiError {
    fn from(error: MatrixSubductionError) -> Self {
        Self::Matrix(Box::new(error))
    }
}

/// Which magnetic isotropy record to use, and how to find it.
///
/// The condensing label has its **own** [`LabelConvention`], independent of the
/// source label's convention.  [`Self::direction`] and [`Self::uni`] are
/// selectors for the same table:
///
/// * a direction resolves through
///   [`isotropy::magnetic_isotropy_subgroup_for_direction`] and accepts
///   [`IsotropyDirection::Label`] / [`IsotropyDirection::Index`]; the existing
///   official-descriptor boundary of the magnetic table is preserved;
/// * a UNI number alone is accepted only when the condensing irrep lists that
///   magnetic group exactly once;
/// * when both are given they must select the same row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MagneticRecordSelector<'a> {
    parent_sg: u8,
    condensing_label: &'a str,
    convention: LabelConvention,
    direction: Option<IsotropyDirection<'a>>,
    uni: Option<usize>,
}

impl<'a> MagneticRecordSelector<'a> {
    /// A selector for one condensing parent irrep; add a direction or a UNI
    /// before calling [`subduce`].
    pub fn new(parent_sg: u8, condensing_label: &'a str, convention: LabelConvention) -> Self {
        Self {
            parent_sg,
            condensing_label,
            convention,
            direction: None,
            uni: None,
        }
    }

    /// Select by an existing isotropy direction of the condensing irrep.
    pub fn direction(mut self, direction: IsotropyDirection<'a>) -> Self {
        self.direction = Some(direction);
        self
    }

    /// Pin the magnetic isotropy record by its UNI number.
    pub fn uni(mut self, uni: usize) -> Self {
        self.uni = Some(uni);
        self
    }

    /// The parent space group number.
    pub fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// The condensing irrep label as requested.
    pub fn condensing_label(&self) -> &'a str {
        self.condensing_label
    }

    /// The requested label convention.
    pub fn convention(&self) -> LabelConvention {
        self.convention
    }

    /// The requested direction, if any.
    pub fn requested_direction(&self) -> Option<IsotropyDirection<'a>> {
        self.direction
    }

    /// The requested UNI number, if any.
    pub fn requested_uni(&self) -> Option<usize> {
        self.uni
    }
}

/// One high-level Γ magnetic subduction request.
///
/// The source and the condensing record are independent identities with their
/// own conventions; `time_reversal` is the parent action's `η` and is required
/// exactly when the selected magnetic group has an antiunitary class.
#[derive(Debug, Clone, Copy)]
pub struct MagneticSubductionRequest<'a> {
    source: SourceIdentity<'a>,
    record: MagneticRecordSelector<'a>,
    embedding: EmbeddingEvidence,
    time_reversal: Option<RealSign>,
}

impl<'a> MagneticSubductionRequest<'a> {
    /// A request for one source irrep and one condensing magnetic record.
    pub fn new(
        source: SourceIdentity<'a>,
        record: MagneticRecordSelector<'a>,
        time_reversal: Option<RealSign>,
    ) -> Self {
        Self {
            source,
            record,
            embedding: EmbeddingEvidence::identity(),
            time_reversal,
        }
    }

    /// The caller's explicit embedding (`map`, `shift`, lattice policy).
    pub fn with_embedding(mut self, embedding: EmbeddingEvidence) -> Self {
        self.embedding = embedding;
        self
    }

    /// The source irrep identity.
    pub fn source(&self) -> SourceIdentity<'a> {
        self.source
    }

    /// The condensing-record selector.
    pub fn record(&self) -> MagneticRecordSelector<'a> {
        self.record
    }

    /// The caller's embedding.
    pub fn embedding(&self) -> EmbeddingEvidence {
        self.embedding
    }

    /// The parent action's `η`, if stated.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.time_reversal
    }
}

/// The verified ordinary-irrep identity of the subduced source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagneticSourceIdentity {
    /// Parent space group number.
    pub parent_sg: u8,
    /// Both genuine label spellings of the source irrep.
    pub labels: IrrepLabels,
    /// Dimension of the source irrep.
    pub dimension: u8,
    /// Whether the source row is double valued.
    pub spinor: bool,
    /// The convention the caller used to name the source.
    pub input_convention: LabelConvention,
}

/// The verified ordinary-irrep identity of the condensing parent irrep.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagneticCondensingIdentity {
    /// Parent space group number.
    pub parent_sg: u8,
    /// Both genuine label spellings of the condensing irrep.
    pub labels: IrrepLabels,
    /// Dimension of the condensing irrep.
    pub dimension: u8,
    /// The convention the caller used to name the condensing irrep.
    pub input_convention: LabelConvention,
}

/// The magnetic isotropy record the request selected.
///
/// `direction_label` is the record's `direction` field, which for a magnetic
/// row is the ISOTROPY direction label (`"P1"`, `"C1"`, …); magnetic rows carry
/// no component descriptor.  `direction_index` is the position in the
/// condensing irrep's magnetic table (the value [`RecordIdentity`] takes),
/// while `ordinal` is the global generated-table index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MagneticRecordIdentity {
    /// Parent space group number.
    pub parent_sg: u8,
    /// Magnetic space group UNI number.
    pub uni: usize,
    /// BNS symbol of the magnetic group.
    pub bns_label: &'static str,
    /// ISOTROPY label of the magnetic group.
    pub iso_label: &'static str,
    /// ISOTROPY direction label.
    pub direction_label: &'static str,
    /// Dimension of the order-parameter space.
    pub direction_dim: u8,
    /// Number of free parameters along the direction.
    pub direction_free: u8,
    /// Position in the condensing irrep's magnetic isotropy table.
    pub direction_index: usize,
    /// Position in the generated magnetic isotropy table.
    pub ordinal: usize,
    /// Primitive basis of the magnetic subgroup in the parent's primitive frame.
    pub basis: [[i32; 3]; 3],
    /// Origin of the magnetic subgroup setting, `[x, y, z, d]` = `(x/d, y/d, z/d)`.
    pub origin: [i32; 4],
}

/// The embedding the verified model actually accepted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MagneticEmbeddingIdentity {
    /// The verified map.
    pub map: Mat3R,
    /// The verified origin shift.
    pub shift: [Rat; 3],
    /// The lattice policy the model was checked against.
    pub lattice_policy: LatticePolicy,
    /// Whether the map realises the record lattice exactly, always recomputed.
    pub realises_record_lattice: bool,
}

/// One magnetic Γ corepresentation of `M` selected by the engine.
///
/// `magnetic_label` is the magnetic summary's display label for the corep.  It
/// is **not** a Bradley–Cracknell name: the summary has no BC field for corep
/// labels, so `magnetic_bc` stays `None`.  The ordinary irrep of the unitary
/// subgroup `H` that the corep is induced from is reported separately in
/// `h_source`, with its own verified CDML/BC labels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagneticTargetIdentity {
    /// Magnetic corep display label from the magnetic summary.
    pub magnetic_label: String,
    /// Magnetic corep BC label.  Always `None`: the summary provides no corep
    /// BC spelling, and this field is never filled from [`Self::h_source`].
    pub magnetic_bc: Option<String>,
    /// Space group number of the target's unitary subgroup `H`.
    pub unitary_sg: u8,
    /// Verified CDML/BC labels of the `H` ordinary irrep the target comes from.
    pub h_source: IrrepLabels,
    /// Dimension of the magnetic corepresentation.
    pub dimension: u32,
    /// Multiplicity of this target in the restriction.
    pub multiplicity: u32,
}

impl MagneticTargetIdentity {
    /// `multiplicity · dimension`, the source dimension this target covers.
    pub fn covered_dimension(&self) -> u32 {
        self.multiplicity * self.dimension
    }
}

/// The one-dimensional result: target identity, complex intertwiner and the
/// reconstructed operator of every magnetic class.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarGammaOutcome {
    /// The selected magnetic corepresentation.
    pub target: MagneticTargetIdentity,
    /// The complex `U` with `D = U⁻¹·A·U` (unitary) / `U⁻¹·A·conj(U)` (anti).
    pub intertwiner: Complex64,
    /// The reconstructed target operator of every class of `M`, in class order.
    pub operators: Vec<TargetOperator>,
}

/// One target piece of a two-dimensional source, with its independent
/// intertwiners and reconstructed target operators.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixGammaPieceOutcome {
    /// The selected magnetic corepresentation.
    pub target: MagneticTargetIdentity,
    /// The `2×d` complex intertwiners, `multiplicity` of them, row-major.
    pub copies: Vec<Vec<Complex64>>,
    /// The reconstructed target operator of every class of `M`, in class order.
    pub operators: Vec<TargetMatrixOperator>,
}

/// The two-dimensional result: one piece per matched target, plus the dimension
/// accounting checked by the engine.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixGammaOutcome {
    /// The matched target pieces.
    pub pieces: Vec<MatrixGammaPieceOutcome>,
    /// Dimension of the two-dimensional source.
    pub source_dimension: u32,
    /// `Σ multiplicityᵢ · dimensionᵢ`, checked equal to `source_dimension`.
    pub covered_dimension: u32,
}

/// The owned result payload of one successful subduction.
#[derive(Debug, Clone, PartialEq)]
pub enum MagneticOutcome {
    /// A one-dimensional source was subduced (R25/R26).
    Scalar(ScalarGammaOutcome),
    /// A two-dimensional real source was subduced (R29).
    Matrix(MatrixGammaOutcome),
}

/// The verified scalar or matrix model the engine built.
///
/// The model is owned by the report; it is the exact object the layer below
/// verified (embedding, both identities, `η`, source action, classes and
/// translation lattice).  Both variants expose the same common accessors here,
/// and the concrete model is available through [`Self::as_scalar`] /
/// [`Self::as_matrix`].
#[derive(Debug)]
pub enum MagneticVerifiedModel {
    /// The R25 scalar Γ model.
    Scalar(ScalarGammaModel),
    /// The R29 two-dimensional real Γ model.
    Matrix(MatrixGammaModel),
}

impl MagneticVerifiedModel {
    /// Parent space group of the source and the record.
    pub fn parent_sg(&self) -> u8 {
        match self {
            Self::Scalar(model) => model.parent_sg(),
            Self::Matrix(model) => model.parent_sg(),
        }
    }

    /// The source parent irrep label, in canonical CDML spelling.
    pub fn source_label(&self) -> &'static str {
        match self {
            Self::Scalar(model) => model.source_label(),
            Self::Matrix(model) => model.source_label(),
        }
    }

    /// The convention the canonical source label resolves under.
    pub fn source_convention(&self) -> LabelConvention {
        match self {
            Self::Scalar(model) => model.source_convention(),
            Self::Matrix(model) => model.source_convention(),
        }
    }

    /// The condensing irrep label, in canonical CDML spelling.
    pub fn record_label(&self) -> &'static str {
        match self {
            Self::Scalar(model) => model.record_label(),
            Self::Matrix(model) => model.record_label(),
        }
    }

    /// The convention the canonical condensing label resolves under.
    pub fn record_convention(&self) -> LabelConvention {
        match self {
            Self::Scalar(model) => model.record_convention(),
            Self::Matrix(model) => model.record_convention(),
        }
    }

    /// Position of the record in the condensing irrep's magnetic table.
    pub fn direction_index(&self) -> usize {
        match self {
            Self::Scalar(model) => model.direction_index(),
            Self::Matrix(model) => model.direction_index(),
        }
    }

    /// Position of the record in the generated magnetic table.
    pub fn ordinal(&self) -> usize {
        match self {
            Self::Scalar(model) => model.ordinal(),
            Self::Matrix(model) => model.ordinal(),
        }
    }

    /// The magnetic space group (UNI) number.
    pub fn uni(&self) -> usize {
        match self {
            Self::Scalar(model) => model.uni(),
            Self::Matrix(model) => model.uni(),
        }
    }

    /// The verified embedding map.
    pub fn map(&self) -> &Mat3R {
        match self {
            Self::Scalar(model) => model.map(),
            Self::Matrix(model) => model.map(),
        }
    }

    /// The verified origin shift.
    pub fn shift(&self) -> &[Rat; 3] {
        match self {
            Self::Scalar(model) => model.shift(),
            Self::Matrix(model) => model.shift(),
        }
    }

    /// The lattice policy the model was checked against.
    pub fn lattice_policy(&self) -> LatticePolicy {
        match self {
            Self::Scalar(model) => model.lattice_policy(),
            Self::Matrix(model) => model.lattice_policy(),
        }
    }

    /// Whether the map realises the record lattice exactly, always recomputed.
    pub fn realises_record_lattice(&self) -> bool {
        match self {
            Self::Scalar(model) => model.realises_record_lattice(),
            Self::Matrix(model) => model.realises_record_lattice(),
        }
    }

    /// The parent action's `η`, or `None` when none was stated.
    pub fn time_reversal(&self) -> Option<RealSign> {
        match self {
            Self::Scalar(model) => model.time_reversal(),
            Self::Matrix(model) => model.time_reversal(),
        }
    }

    /// The classes of the magnetic group, in the engine's class order.
    pub fn classes(&self) -> &[MagneticOperation] {
        match self {
            Self::Scalar(model) => model.classes(),
            Self::Matrix(model) => model.classes(),
        }
    }

    /// The group's translation lattice.
    pub fn translation_lattice(&self) -> &[[Rat; 3]] {
        match self {
            Self::Scalar(model) => model.translation_lattice(),
            Self::Matrix(model) => model.translation_lattice(),
        }
    }

    /// The concrete scalar model, when the source was one-dimensional.
    pub fn as_scalar(&self) -> Option<&ScalarGammaModel> {
        match self {
            Self::Scalar(model) => Some(model),
            Self::Matrix(_) => None,
        }
    }

    /// The concrete matrix model, when the source was two-dimensional.
    pub fn as_matrix(&self) -> Option<&MatrixGammaModel> {
        match self {
            Self::Scalar(_) => None,
            Self::Matrix(model) => Some(model),
        }
    }
}

/// The complete answer to one [`MagneticSubductionRequest`].
///
/// The report owns the verified model and an owned snapshot of the result.  It
/// exposes the identities resolved from the request, the selected magnetic
/// record, the accepted embedding and the engine's result payload; no field is
/// borrowed through a model stored inside the same struct.
#[derive(Debug)]
pub struct MagneticSubductionReport {
    model: MagneticVerifiedModel,
    source: MagneticSourceIdentity,
    condensing: MagneticCondensingIdentity,
    record: MagneticRecordIdentity,
    embedding: MagneticEmbeddingIdentity,
    outcome: MagneticOutcome,
}

impl MagneticSubductionReport {
    /// The verified model this report was computed from.
    pub fn model(&self) -> &MagneticVerifiedModel {
        &self.model
    }

    /// The concrete scalar model, when the source was one-dimensional.
    pub fn scalar_model(&self) -> Option<&ScalarGammaModel> {
        self.model.as_scalar()
    }

    /// The concrete matrix model, when the source was two-dimensional.
    pub fn matrix_model(&self) -> Option<&MatrixGammaModel> {
        self.model.as_matrix()
    }

    /// The verified subduced source identity.
    pub fn source(&self) -> &MagneticSourceIdentity {
        &self.source
    }

    /// The verified condensing irrep identity.
    pub fn condensing(&self) -> &MagneticCondensingIdentity {
        &self.condensing
    }

    /// The selected magnetic isotropy record.
    pub fn record(&self) -> &MagneticRecordIdentity {
        &self.record
    }

    /// The accepted embedding.
    pub fn embedding(&self) -> &MagneticEmbeddingIdentity {
        &self.embedding
    }

    /// The engine's result payload.
    pub fn outcome(&self) -> &MagneticOutcome {
        &self.outcome
    }

    /// The one-dimensional payload, when the source was scalar.
    pub fn scalar(&self) -> Option<&ScalarGammaOutcome> {
        match &self.outcome {
            MagneticOutcome::Scalar(outcome) => Some(outcome),
            MagneticOutcome::Matrix(_) => None,
        }
    }

    /// The two-dimensional payload, when the source was a real matrix.
    pub fn matrix(&self) -> Option<&MatrixGammaOutcome> {
        match &self.outcome {
            MagneticOutcome::Scalar(_) => None,
            MagneticOutcome::Matrix(outcome) => Some(outcome),
        }
    }

    /// Dimension of the source irrep.
    pub fn source_dimension(&self) -> u32 {
        u32::from(self.source.dimension)
    }

    /// Dimension covered by the result: the source dimension the target pieces
    /// account for, as checked by the engine.
    pub fn covered_dimension(&self) -> u32 {
        match &self.outcome {
            MagneticOutcome::Scalar(outcome) => outcome.target.covered_dimension(),
            MagneticOutcome::Matrix(outcome) => outcome.covered_dimension,
        }
    }

    /// The parent action's `η`, or `None` when none was stated.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.model.time_reversal()
    }

    /// The parity of the parent action's `η`, which is `None` exactly when no
    /// `η` was stated.  This is not the square of `T`; see the R25 contract.
    pub fn time_parity(&self) -> Option<TimeParity> {
        self.time_reversal().map(|sign| match sign {
            RealSign::Plus => TimeParity::Even,
            RealSign::Minus => TimeParity::Odd,
        })
    }

    /// The classes of the magnetic group, in the result's class order.
    pub fn classes(&self) -> &[MagneticOperation] {
        self.model.classes()
    }

    /// The translation lattice the class comparisons reduced modulo.
    pub fn translation_lattice(&self) -> &[[Rat; 3]] {
        self.model.translation_lattice()
    }
}

/// Resolve one source label to the unique irrep it names.
fn resolve_source(
    source: SourceIdentity<'_>,
) -> Result<&'static IrrepRecord, MagneticSubductionApiError> {
    if !(1..=230).contains(&source.parent_sg) {
        return Err(MagneticSubductionApiError::UnsupportedSpaceGroup {
            sg: source.parent_sg,
        });
    }
    let mut matches = query::find_irreps(source.parent_sg, source.label, source.convention);
    match matches.len() {
        0 => Err(MagneticSubductionApiError::SourceLabelNotFound {
            sg: source.parent_sg,
            label: source.label.to_string(),
            convention: source.convention,
        }),
        1 => Ok(matches.pop().expect("exactly one match")),
        _ => {
            let mut candidates: Vec<&'static str> =
                matches.iter().map(|record| record.ml).collect();
            candidates.sort_unstable();
            candidates.dedup();
            Err(MagneticSubductionApiError::AmbiguousSourceLabel {
                sg: source.parent_sg,
                label: source.label.to_string(),
                convention: source.convention,
                candidates,
            })
        }
    }
}

/// Resolve the condensing record under its own convention and the requested
/// direction/UNI selector.
fn resolve_record(
    selector: &MagneticRecordSelector<'_>,
) -> Result<(MagneticIsotropySubgroup, usize), MagneticSubductionApiError> {
    if !(1..=230).contains(&selector.parent_sg) {
        return Err(MagneticSubductionApiError::UnsupportedSpaceGroup {
            sg: selector.parent_sg,
        });
    }
    let rows = isotropy::magnetic_isotropy_subgroups(
        selector.parent_sg,
        selector.condensing_label,
        selector.convention,
    )?;
    let local_index = |subgroup: &MagneticIsotropySubgroup| {
        rows.iter()
            .position(|row| row.ordinal == subgroup.ordinal)
            .ok_or(MagneticSubductionApiError::Isotropy(
                IsotropyError::InvalidIsotropyRecord {
                    ordinal: subgroup.ordinal,
                    len: rows.len(),
                },
            ))
    };
    match (selector.direction, selector.uni) {
        (None, None) => Err(MagneticSubductionApiError::RecordSelectorMissing {
            sg: selector.parent_sg,
        }),
        (Some(direction), uni) => {
            let subgroup = isotropy::magnetic_isotropy_subgroup_for_direction(
                selector.parent_sg,
                selector.condensing_label,
                selector.convention,
                direction,
            )?;
            let direction_index = local_index(&subgroup)?;
            if let Some(uni) = uni
                && subgroup.record.mag_sg != uni
            {
                return Err(MagneticSubductionApiError::RecordUniMismatch {
                    sg: selector.parent_sg,
                    direction_index,
                    requested: uni,
                    found: subgroup.record.mag_sg,
                });
            }
            Ok((subgroup, direction_index))
        }
        (None, Some(uni)) => {
            let candidates: Vec<(usize, &'static str)> = rows
                .iter()
                .enumerate()
                .filter(|(_, row)| row.record.mag_sg == uni)
                .map(|(index, row)| (index, row.record.direction))
                .collect();
            match candidates.as_slice() {
                [] => Err(MagneticSubductionApiError::RecordUniNotFound {
                    sg: selector.parent_sg,
                    label: selector.condensing_label.to_string(),
                    convention: selector.convention,
                    uni,
                }),
                [(index, _)] => Ok((rows[*index], *index)),
                _ => Err(MagneticSubductionApiError::AmbiguousRecordUni {
                    sg: selector.parent_sg,
                    uni,
                    candidates,
                }),
            }
        }
    }
}

/// CDML/BC labels of the target's ordinary `H` irrep, read from the database.
///
/// The engine reports the `H` irrep by its canonical CDML label; that identity
/// must resolve to **exactly one** row before it may be reported.  A missing or
/// non-unique lookup is a named refusal: it means the reported `H` identity
/// cannot be anchored in the database, which is a different statement from "this
/// row has no BC mapping".  A legitimate `bc: None` from
/// [`IrrepRecord::labels`] (the `"***"` placeholder) is preserved as the normal
/// missing-BC datum, not turned into an error.
fn h_source_labels(sg: u8, ml: &'static str) -> Result<IrrepLabels, MagneticSubductionApiError> {
    let matches = query::find_irreps(sg, ml, LabelConvention::Cdml);
    match matches.as_slice() {
        [record] => Ok(record.labels()),
        _ => Err(MagneticSubductionApiError::HSourceIdentityUnresolved {
            sg,
            label: ml,
            matches: matches.len(),
        }),
    }
}

/// The request-resolved identities and engine inputs shared by the legacy
/// dispatch and the restricted type-C entry point.
struct ResolvedMagneticRequest<'a> {
    source: SourceIdentity<'a>,
    source_identity: MagneticSourceIdentity,
    condensing_identity: MagneticCondensingIdentity,
    record_identity: MagneticRecordIdentity,
    engine_record: RecordIdentity<'a>,
    embedding: EmbeddingEvidence,
    time_reversal: Option<RealSign>,
}

/// Resolve one request exactly as [`subduce`] always has: parent check, source
/// label, source dimension, condensing row and the four identities.
fn resolve_magnetic_request<'a>(
    request: &MagneticSubductionRequest<'a>,
) -> Result<ResolvedMagneticRequest<'a>, MagneticSubductionApiError> {
    let source = request.source;
    let record = request.record;
    if source.parent_sg != record.parent_sg {
        return Err(MagneticSubductionApiError::ParentMismatch {
            source_sg: source.parent_sg,
            record_sg: record.parent_sg,
        });
    }
    let source_record = resolve_source(source)?;
    // The source dimension decides which verified engine can answer at all, so
    // an out-of-family source is refused before record resolution: the answer
    // does not depend on a condensing row that neither entry point could use.
    if source_record.dim != 1 && source_record.dim != 2 {
        return Err(MagneticSubductionApiError::SourceDimensionUnsupported {
            sg: source.parent_sg,
            label: source_record.ml,
            dimension: source_record.dim,
        });
    }
    let (subgroup, direction_index) = resolve_record(&record)?;

    let source_identity = MagneticSourceIdentity {
        parent_sg: source.parent_sg,
        labels: source_record.labels(),
        dimension: source_record.dim,
        spinor: source_record.spinor,
        input_convention: source.convention,
    };
    let condensing_identity = MagneticCondensingIdentity {
        parent_sg: record.parent_sg,
        labels: subgroup.labels(),
        dimension: subgroup.irrep_dim,
        input_convention: record.convention,
    };
    let record_identity = MagneticRecordIdentity {
        parent_sg: record.parent_sg,
        uni: subgroup.record.mag_sg,
        bns_label: subgroup.record.bns_label,
        iso_label: subgroup.record.iso_label,
        direction_label: subgroup.record.direction,
        direction_dim: subgroup.record.direction_dim,
        direction_free: subgroup.record.direction_free,
        direction_index,
        ordinal: subgroup.ordinal,
        basis: subgroup.record.basis,
        origin: subgroup.record.origin,
    };

    // The record identity is passed to the engine in canonical CDML form with
    // the local row index already resolved; the caller's spelling is preserved
    // in `condensing_identity.input_convention`, so a BC request reaches the
    // same row without the engine having to re-do the selector.
    let engine_record = RecordIdentity {
        parent_sg: record.parent_sg,
        condensing_label: subgroup.irrep_ml,
        convention: LabelConvention::Cdml,
        direction_index,
    };
    Ok(ResolvedMagneticRequest {
        source,
        source_identity,
        condensing_identity,
        record_identity,
        engine_record,
        embedding: request.embedding,
        time_reversal: request.time_reversal,
    })
}

/// Answer one [`MagneticSubductionRequest`].
///
/// The source dimension decides which verified engine is called: dimension one
/// goes through [`magnetic_input::verify_scalar_gamma_input`] then
/// [`magnetic_subduction::subduce_scalar_gamma`]; dimension two goes through
/// [`magnetic_matrix_subduction::verify_matrix_gamma_input`] then
/// [`magnetic_matrix_subduction::subduce_matrix_gamma`].  Any other dimension is
/// refused by name, because neither engine's mathematical scope contains it.
pub fn subduce(
    request: &MagneticSubductionRequest<'_>,
) -> Result<MagneticSubductionReport, MagneticSubductionApiError> {
    let ResolvedMagneticRequest {
        source,
        source_identity,
        condensing_identity,
        record_identity,
        engine_record,
        embedding,
        time_reversal,
    } = resolve_magnetic_request(request)?;

    let (model, outcome) = match source_identity.dimension {
        1 => {
            let scalar_model = verify_scalar_gamma_input(
                &ScalarGammaRequest::new(source, engine_record, time_reversal)
                    .with_embedding(embedding),
            )?;
            let outcome = {
                let result = scalar::subduce_scalar_gamma(&scalar_model)?;
                let target = MagneticTargetIdentity {
                    magnetic_label: result.target_label().to_string(),
                    magnetic_bc: None,
                    unitary_sg: result.target_unitary_sg(),
                    h_source: h_source_labels(
                        result.target_unitary_sg(),
                        result.target_source_label(),
                    )?,
                    dimension: result.target_dimension(),
                    multiplicity: result.multiplicity(),
                };
                MagneticOutcome::Scalar(ScalarGammaOutcome {
                    target,
                    intertwiner: result.intertwiner(),
                    operators: result.operators().to_vec(),
                })
            };
            (MagneticVerifiedModel::Scalar(scalar_model), outcome)
        }
        2 => {
            let matrix_model = matrix::verify_matrix_gamma_input(
                &MatrixGammaRequest::new(source, engine_record, time_reversal)
                    .with_embedding(embedding),
            )?;
            let outcome = {
                let result = matrix::subduce_matrix_gamma(&matrix_model)?;
                let mut pieces = Vec::with_capacity(result.pieces().len());
                for piece in result.pieces() {
                    pieces.push(MatrixGammaPieceOutcome {
                        target: MagneticTargetIdentity {
                            magnetic_label: piece.label().to_string(),
                            magnetic_bc: None,
                            unitary_sg: piece.unitary_sg(),
                            h_source: h_source_labels(piece.unitary_sg(), piece.source_label())?,
                            dimension: piece.dimension(),
                            multiplicity: piece.multiplicity(),
                        },
                        copies: piece.copies().to_vec(),
                        operators: piece.operators().to_vec(),
                    });
                }
                MagneticOutcome::Matrix(MatrixGammaOutcome {
                    pieces,
                    source_dimension: result.source_dimension(),
                    covered_dimension: result.covered_dimension(),
                })
            };
            (MagneticVerifiedModel::Matrix(matrix_model), outcome)
        }
        dimension => {
            return Err(MagneticSubductionApiError::SourceDimensionUnsupported {
                sg: source.parent_sg,
                label: source_identity.labels.cdml,
                dimension,
            });
        }
    };

    let embedding_identity = MagneticEmbeddingIdentity {
        map: *model.map(),
        shift: *model.shift(),
        lattice_policy: model.lattice_policy(),
        realises_record_lattice: model.realises_record_lattice(),
    };
    Ok(MagneticSubductionReport {
        model,
        source: source_identity,
        condensing: condensing_identity,
        record: record_identity,
        embedding: embedding_identity,
        outcome,
    })
}

// ---------------------------------------------------------------------------
// R32: restricted paired type-C entry point
// ---------------------------------------------------------------------------

/// One magnetic Γ type-C corepresentation of `M`, reported with **both** of its
/// unitary-subgroup `H` source irreps.
///
/// The legacy [`MagneticTargetIdentity`] has a single `h_source` field and is
/// deliberately unchanged: adding a required partner field or a new
/// [`MagneticOutcome`] variant would break existing struct literals and
/// exhaustive matches.  This separate identity therefore carries the real pair.
/// The magnetic corep display label stays a display label; `magnetic_bc` is
/// always `None` because the summary has no BC field for corep labels.  The two
/// `h_sources` entries are genuine ordinary-irrep identities with their own
/// CDML and BC spellings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagneticTypeCTargetIdentity {
    /// Magnetic corep display label from the magnetic summary.
    pub magnetic_label: String,
    /// Magnetic corep BC label.  Always `None`; never filled from an H source.
    pub magnetic_bc: Option<String>,
    /// Space group number of the target's unitary subgroup `H`.
    pub unitary_sg: u8,
    /// The two verified H ordinary-irrep identities, in the summary's source
    /// order.
    pub h_sources: [IrrepLabels; 2],
    /// Dimension of the magnetic corepresentation (two).
    pub dimension: u32,
    /// Multiplicity of this target in the restriction (one).
    pub multiplicity: u32,
}

impl MagneticTypeCTargetIdentity {
    /// `multiplicity · dimension`, the source dimension this target covers.
    pub fn covered_dimension(&self) -> u32 {
        self.multiplicity * self.dimension
    }
}

/// The owned payload of one restricted type-C subduction.
#[derive(Debug, Clone, PartialEq)]
pub struct MagneticTypeCOutcome {
    /// The selected magnetic type-C corepresentation with both H identities.
    pub target: MagneticTypeCTargetIdentity,
    /// The normalized `2×2` complex intertwiners, `multiplicity` of them,
    /// row-major.
    pub copies: Vec<Vec<Complex64>>,
    /// The reconstructed target operator of every class of `M`, in class order.
    pub operators: Vec<TargetMatrixOperator>,
    /// `Σ multiplicityᵢ · dimensionᵢ`, checked equal to the source dimension.
    pub covered_dimension: u32,
}

/// Why a restricted type-C request could not be answered.
///
/// This is a new enum, so [`MagneticSubductionApiError`] and its exhaustive
/// matches are untouched.  The request-resolution layer and the paired-type-C
/// engine keep their own errors as sources.
#[derive(Debug)]
pub enum MagneticTypeCError {
    /// The request could not be resolved into a source, a record and the four
    /// identities.
    Request(Box<MagneticSubductionApiError>),
    /// The paired type-C engine refused the verified model.
    Engine(Box<MatrixTypeCError>),
    /// The resolved source is not the two-dimensional real matrix source this
    /// entry point requires.
    SourceDimensionNotTwo {
        sg: u8,
        label: &'static str,
        dimension: u8,
    },
}

impl std::fmt::Display for MagneticTypeCError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Request(error) => write!(f, "{error}"),
            Self::Engine(error) => write!(f, "{error}"),
            Self::SourceDimensionNotTwo {
                sg,
                label,
                dimension,
            } => write!(
                f,
                "SG {sg} source {label:?} has dimension {dimension}; the paired type-C entry point \
                 requires the two-dimensional real matrix source"
            ),
        }
    }
}

impl std::error::Error for MagneticTypeCError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Request(error) => Some(error.as_ref()),
            Self::Engine(error) => Some(error.as_ref()),
            Self::SourceDimensionNotTwo { .. } => None,
        }
    }
}

impl From<MagneticSubductionApiError> for MagneticTypeCError {
    fn from(error: MagneticSubductionApiError) -> Self {
        Self::Request(Box::new(error))
    }
}

impl From<MatrixTypeCError> for MagneticTypeCError {
    fn from(error: MatrixTypeCError) -> Self {
        Self::Engine(Box::new(error))
    }
}

impl From<MatrixSubductionError> for MagneticTypeCError {
    fn from(error: MatrixSubductionError) -> Self {
        Self::Request(Box::new(MagneticSubductionApiError::Matrix(Box::new(
            error,
        ))))
    }
}

/// The complete answer to one restricted paired type-C request.
///
/// The report owns the actual verified [`MatrixGammaModel`] (there is no scalar
/// branch here) and the same four identity surfaces as
/// [`MagneticSubductionReport`].  The legacy report and outcome types are
/// unchanged.
#[derive(Debug)]
pub struct MagneticTypeCReport {
    model: MatrixGammaModel,
    source: MagneticSourceIdentity,
    condensing: MagneticCondensingIdentity,
    record: MagneticRecordIdentity,
    embedding: MagneticEmbeddingIdentity,
    outcome: MagneticTypeCOutcome,
}

impl MagneticTypeCReport {
    /// The verified two-dimensional real Γ matrix model.
    pub fn model(&self) -> &MatrixGammaModel {
        &self.model
    }

    /// The verified two-dimensional real Γ matrix model.
    pub fn matrix_model(&self) -> &MatrixGammaModel {
        &self.model
    }

    /// The verified subduced source identity.
    pub fn source(&self) -> &MagneticSourceIdentity {
        &self.source
    }

    /// The verified condensing irrep identity.
    pub fn condensing(&self) -> &MagneticCondensingIdentity {
        &self.condensing
    }

    /// The selected magnetic isotropy record.
    pub fn record(&self) -> &MagneticRecordIdentity {
        &self.record
    }

    /// The accepted embedding.
    pub fn embedding(&self) -> &MagneticEmbeddingIdentity {
        &self.embedding
    }

    /// The paired type-C payload.
    pub fn outcome(&self) -> &MagneticTypeCOutcome {
        &self.outcome
    }

    /// The selected type-C target with both H identities.
    pub fn target(&self) -> &MagneticTypeCTargetIdentity {
        &self.outcome.target
    }

    /// Dimension of the source irrep (two).
    pub fn source_dimension(&self) -> u32 {
        u32::from(self.source.dimension)
    }

    /// `multiplicity · dimension`, checked equal to the source dimension.
    pub fn covered_dimension(&self) -> u32 {
        self.outcome.covered_dimension
    }

    /// The parent action's `η`, or `None` when none was stated.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.model.time_reversal()
    }

    /// The parity of the parent action's `η`; `None` exactly when no `η` was
    /// stated.  This is not the square of `T`.
    pub fn time_parity(&self) -> Option<TimeParity> {
        self.time_reversal().map(|sign| match sign {
            RealSign::Plus => TimeParity::Even,
            RealSign::Minus => TimeParity::Odd,
        })
    }

    /// The classes of the magnetic group, in the result's class order.
    pub fn classes(&self) -> &[MagneticOperation] {
        self.model.classes()
    }

    /// The translation lattice the class comparisons reduced modulo.
    pub fn translation_lattice(&self) -> &[[Rat; 3]] {
        self.model.translation_lattice()
    }
}

/// Answer one restricted paired type-C request.
///
/// This is the incremental companion of [`subduce`]: it resolves the same
/// request through the same private resolver and the same
/// [`matrix::verify_matrix_gamma_input`] model, then calls the restricted
/// [`matrix::subduce_matrix_gamma_type_c`] engine.  Only a two-dimensional
/// source is accepted; the legacy [`subduce`] keeps its type-A-only matrix
/// domain and refuses a type-C record with its existing named error.
pub fn subduce_type_c(
    request: &MagneticSubductionRequest<'_>,
) -> Result<MagneticTypeCReport, MagneticTypeCError> {
    let ResolvedMagneticRequest {
        source,
        source_identity,
        condensing_identity,
        record_identity,
        engine_record,
        embedding,
        time_reversal,
    } = resolve_magnetic_request(request)?;
    if source_identity.dimension != 2 {
        return Err(MagneticTypeCError::SourceDimensionNotTwo {
            sg: source.parent_sg,
            label: source_identity.labels.cdml,
            dimension: source_identity.dimension,
        });
    }
    let matrix_model = matrix::verify_matrix_gamma_input(
        &MatrixGammaRequest::new(source, engine_record, time_reversal).with_embedding(embedding),
    )?;
    let outcome = {
        let result = matrix::subduce_matrix_gamma_type_c(&matrix_model)?;
        let piece = result.piece();
        let h_sources = piece.h_sources();
        let target = MagneticTypeCTargetIdentity {
            magnetic_label: piece.label().to_string(),
            magnetic_bc: None,
            unitary_sg: piece.unitary_sg(),
            h_sources: [
                h_source_labels(piece.unitary_sg(), h_sources[0])?,
                h_source_labels(piece.unitary_sg(), h_sources[1])?,
            ],
            dimension: piece.dimension(),
            multiplicity: piece.multiplicity(),
        };
        MagneticTypeCOutcome {
            target,
            copies: piece.copies().to_vec(),
            operators: piece.operators().to_vec(),
            covered_dimension: result.covered_dimension(),
        }
    };
    let embedding_identity = MagneticEmbeddingIdentity {
        map: *matrix_model.map(),
        shift: *matrix_model.shift(),
        lattice_policy: matrix_model.lattice_policy(),
        realises_record_lattice: matrix_model.realises_record_lattice(),
    };
    Ok(MagneticTypeCReport {
        model: matrix_model,
        source: source_identity,
        condensing: condensing_identity,
        record: record_identity,
        embedding: embedding_identity,
        outcome,
    })
}
