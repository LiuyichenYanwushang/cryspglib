//! R29: genuine two-dimensional real Γ magnetic **matrix** subduction.
//!
//! The R25/R26 entry points are limited to real `±1`, one-dimensional Γ sources.
//! This module is the matrix analogue for a **true, non-compound two-dimensional
//! real Γ parent irrep**: it certifies operation-aware source matrices, builds the
//! parent grey action on the verified embedding, decomposes the restriction into
//! the magnetic group's Γ corepresentations, and exhibits matrix intertwiners
//! with a full-space reconstruction.
//!
//! # Why the raw matrix storage needs its own certification
//!
//! [`crate::irrep::types::IrrepRecord::matrices`] is the legacy raw/source
//! storage.  Its own documentation refuses to call it an operation-aware matrix
//! row: the source representatives and basis are not *proven* aligned with
//! [`crate::irrep::types::IrrepRecord::pir_rotations`] /
//! [`crate::irrep::types::IrrepRecord::pir_translations`].  A matching length is
//! not a proof.  This module therefore never zips the two arrays silently: for
//! every accepted source (and every accepted target H-irrep) it runs
//!
//! 1. shape/finiteness of the `m`×`d²` storage against the `m` operation rows;
//! 2. exact orthogonality `Mᵀ M = I` of every stored matrix — this is what
//!    justifies the Euclidean metric and the `U†U = I` normalization later, so
//!    it is *checked*, not assumed;
//! 3. `trace(Mᵢ) = χ(i)` against the legacy character row **and** against the
//!    typed selected-arm Γ row
//!    ([`crate::irrep::types::IrrepRecord::ordinary_scalar_selected_arm_block_trace`]),
//!    which is the character surface the production corep path reads;
//! 4. the **complete ordered-pair matrix group law** on the stored operation
//!    list, reduced modulo that list's own translation lattice (this is the
//!    discriminating check: a same-χ permutation of two operation rows keeps the
//!    traces and breaks the law);
//! 5. irreducibility over the reals (the real commutant is one-dimensional).
//!
//! What that establishes is a **Γ-equivalent certified representation basis** on
//! the stored operation list — not a proof of the original PIR provenance of
//! every stored byte.  The storage caveat covers the source representatives
//! **and** the basis/gauge, and `k = 0` only makes the representative half's
//! Bloch phase trivial; it proves nothing about the basis.  This module therefore
//! does not rest on that phase argument: the χ agreement on two storage surfaces,
//! the complete group law and the class binding below **independently** establish
//! an equivalent real representation on the stored operation list, and Γ is the
//! condition of that claim (non-Γ sources are refused), not a proof of the
//! storage semantics.  A mis-associated row fails these checks by name instead
//! of being read.
//!
//! # Operations, provenance and the grey action
//!
//! Every magnetic class is mapped into the parent's frame with the *same*
//! geometry as R25 ([`crate::irrep::magnetic_input`]'s parent frame and class
//! mapping, reused crate-visibly) and is then bound to the unique stored
//! operation row with the same rotation whose translation differs by a parent
//! lattice vector.  The binding is not a fallback: nothing matches, or two rows
//! match, is a named refusal.
//!
//! The parent grey action on the certified basis is
//!
//! ```text
//! A(p)     = D(p)                  for a unitary class
//! A(p, Θ)  = η · D(p) ∘ conj       for an antiunitary class of sign η = ±1
//! ```
//!
//! and whether that really is a corepresentation of `M` is not assumed: the
//! constructor checks `A(g)·[conj(A(h)) if t(g)] = A(gh)` on every ordered pair
//! of the model's classes, on the exact class composition of the magnetic group.
//! (For a real `D` the conjugation is an involution on the entries, so this is
//! exactly the statement that the chosen parent basis is real on the restricted
//! group.)  Both `η` signs are legal inputs and behave as in R26: they are a
//! gauge of the antiunitary coset and generally give different intertwiners.
//!
//! # Targets and the decomposition
//!
//! The target identity is read **independently** from the magnetic corep summary
//! at Γ (R26/R27's `MagneticCorepSummary`: label, dimension, Wigner type,
//! source H-irrep).  For each candidate the unitary matrices are rebuilt from
//! *that H-irrep's own* stored matrices through the H data-Hall frame used by the
//! summary ([`crate::irrep::corep`]'s unitary-subgroup transform, `find_seitz`,
//! `build_h_to_irrep_op_map`), and every unitary class is cross-checked against
//! the summary's own character value.  The antiunitary part is **constructed
//! from the actual multiplication of `M`**, never read from the pending
//! antiunitary columns:
//!
//! * solve the intertwining system `U_h X = X U_{σ(h)}` for the coset matrix `X`
//!   (`σ` is the conjugation of a unitary class by the canonical antiunitary
//!   representative `a₀`), then
//! * fix the scale by the square `X² = U_{a₀²}`,
//! * define `D(a₀ h) = X · U_h` and verify the square, the conjugation
//!   compatibility `D(h)D(a₀) = D(h a₀)` and the full ordered-pair law.
//!
//! The restriction is then decomposed: every candidate whose intertwiner space is
//! non-zero contributes `m` independent, individually normalized `2×d` copies
//! (`U†U = I_d`, taken greedily by ℂ-rank from the real-linear solution space,
//! which carries the conjugation of the antiunitary branch).  A successful result
//! must satisfy `Σ mᵢ dᵢ = dim(source) = 2` **and** the concatenated copies must
//! have full ℂ-rank, so a partial decomposition or a single identity-like example
//! is refused by name.  The reported operators are re-verified against the source
//! action class by class on complex vectors,
//! `S(g)∘φ = φ∘T(g)` with `φ(z) = U z` for every copy, where the antiunitary
//! branch conjugates the **image vector `U z`**, not merely `z`.
//!
//! # Scope fences
//!
//! Source: non-spinor, non-compound, `k = 0`, real two-dimensional Γ parent
//! irrep whose stored matrix table passes the certification above.  Target:
//! Wigner **type A**, single non-spinor source H-irrep of dimension 1 or 2 whose
//! stored matrices pass the same certification.  Spinor, compound, non-Γ,
//! higher-dimensional, type B/C and complex targets are refused by name or
//! skipped as out of shape; no setting search, no phase substitution, no
//! fallback `U`, and no new matrix framework is introduced.
//!
//! # R32: restricted paired type-C target
//!
//! [`subduce_matrix_gamma_type_c`] is the incremental companion for the type-C
//! slice of the same certified source model: the target is a two-dimensional
//! Γ corepresentation whose unitary restriction is **exactly two inequivalent
//! non-spinor real `±1` one-dimensional H irreps**, exchanged by the
//! antiunitary coset (both `D(a₀²) = +I` and `-I` are legal).  Its target table
//! is constructed from the two independent H character tables and the actual
//! `a₀`, `a₀²` and `a₀⁻¹ h a₀`; the twisted pairing, the square and every
//! ordered pair are checked before the complex intertwiner is matched.  The
//! real-linear Hom space of this complex (Schur-C) corep is two-dimensional
//! (`I` and `diag(i, -i)` both intertwine it), so its dimension is deliberately
//! **not** the multiplicity: the copy count is the number of normalized,
//! pairwise-orthogonal copies under the source-dimension cap, and the success
//! evidence is `U†U = I₂`, full complex rank, covered dimension 2 and the
//! complex-vector reconstruction.  The legacy [`MatrixSubductionPiece`],
//! [`MatrixSubductionError`] and [`subduce_matrix_gamma`] are unchanged; a
//! type-C-only record is still refused by the legacy entry point with
//! [`MatrixSubductionError::NoMatchingTarget`].

use num_complex::Complex64;

use crate::irrep::LabelConvention;
use crate::irrep::corep::{self, CorepType};
use crate::irrep::magnetic_embedding::{
    self, GeometryError, MagneticOperation, canonical_translation, compose, magnetic_operations,
};
use crate::irrep::magnetic_input::{
    self, EmbeddingEvidence, LatticePolicy, MappedClass, ParentFrame, RealSign, RecordIdentity,
    SourceIdentity,
};
use crate::irrep::magnetic_subduction;
use crate::irrep::magnetic_summary::{self, MagneticCorepSummary, MagneticKPointSummary};
use crate::irrep::query;
use crate::irrep::subduction::{Mat3I, Rat};
use crate::irrep::types::{CharacterRow, IrrepRecord};
use crate::irrep::wigner;

/// The Γ (k = 0) k-point coordinates the summary tables use.
const GAMMA: (i8, i8, i8, i8) = (0, 0, 0, 1);

/// The source dimension this entry point accepts.
const SOURCE_DIMENSION: usize = 2;

/// Tolerance of the exact-grid promotions and of the certified arithmetic.
const GRID_TOLERANCE: f64 = 1e-9;

/// Denominator grid a stored translation is snapped to.
const MAX_TRANSLATION_DENOMINATOR: i128 = 96;

/// The one-line scope fence of [`subduce_matrix_gamma`].
pub const SCOPE: &str = "Two-dimensional real non-spinor non-compound Gamma parent irrep (certified \
operation-aware matrix table) restricted to an embedded magnetic group; target: type-A real \
single-source 1- or 2-dimensional Gamma corepresentation; spinor, compound, non-Gamma, complex, \
higher-dimensional, type-B/C and setting-search cases are outside this entry point.";

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Why one stored irrep matrix table could not be certified as an
/// operation-aware real Γ representation.
///
/// Every variant is a specific refusal.  A same-χ permutation of two operation
/// rows keeps [`Self::TraceMismatch`] silent and fails
/// [`Self::LawViolated`] instead, so the two checks are not substitutes.
#[derive(Debug, Clone, PartialEq)]
pub enum MatrixCertificationProblem {
    /// The record's declared dimension is not the requested one.
    DeclaredDimensionMismatch { declared: usize, expected: usize },
    /// The requested dimension is outside the 1/2 matrix family.
    DimensionOutOfScope { dimension: usize },
    /// The operation rows, matrix entries or traces do not have the shape
    /// `operators × dimension²` / `operators`.
    ShapeMismatch {
        operators: usize,
        matrices: usize,
        traces: usize,
        dimension: usize,
    },
    /// A stored matrix or translation is not finite.
    NonFinite { operation: usize },
    /// A stored translation is off the exact small-denominator grid.
    TranslationOffGrid { operation: usize },
    /// `Mᵀ M` is not the identity for this operation, so the stored basis is not
    /// orthogonal and the Euclidean metric may not be used unqualified.
    NotOrthogonal { operation: usize },
    /// The operation list has no identity rotation with zero canonical
    /// translation.
    IdentityMissing,
    /// A unique identity operation exists but its matrix is not the identity.
    IdentityNotIdentity { operation: usize },
    /// `trace(Mᵢ)` disagrees with the stored character row.
    TraceMismatch {
        operation: usize,
        matrix_trace: f64,
        stored_trace: f64,
    },
    /// The ordered product of two operation rows lands on no stored row
    /// (relative to the list's own translation lattice).
    LawProductNotUnique {
        left: usize,
        right: usize,
        products: usize,
    },
    /// Several lattice-equivalent stored rows represent the product but do not
    /// carry the same matrix, so the storage is not a function of the group
    /// element even at Γ.
    LawProductAmbiguous {
        left: usize,
        right: usize,
        products: usize,
    },
    /// The stored matrices do not multiply like the stored operations.
    LawViolated {
        left: usize,
        right: usize,
        product: usize,
    },
    /// The real commutant is not one-dimensional, so the matrices are not one
    /// real irreducible representation.
    Reducible { commutant_dimension: usize },
    /// The typed selected-arm Γ row is unavailable for this record.
    TypedRowMissing,
    /// The typed row's shape disagrees with the stored operation list.
    TypedRowShapeMismatch {
        row: usize,
        operators: usize,
        dimension: usize,
    },
    /// A typed-row operation does not equal the stored operation at that index.
    TypedRowOperationMismatch { operation: usize },
    /// The typed row's value disagrees with the matrix trace.
    TypedRowMismatch {
        operation: usize,
        matrix_trace: f64,
        row_value: f64,
    },
    /// The typed row's value is not real at this Γ row.
    TypedRowNotReal { operation: usize, imaginary: f64 },
}

impl std::fmt::Display for MatrixCertificationProblem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::DeclaredDimensionMismatch { declared, expected } => write!(
                f,
                "declared dimension {declared} is not the requested {expected}"
            ),
            Self::DimensionOutOfScope { dimension } => {
                write!(f, "dimension {dimension} is outside the 1/2 matrix family")
            }
            Self::ShapeMismatch {
                operators,
                matrices,
                traces,
                dimension,
            } => write!(
                f,
                "shape mismatch: operators={operators} matrices={matrices} traces={traces} dimension={dimension}"
            ),
            Self::NonFinite { operation } => {
                write!(
                    f,
                    "operation {operation} has a non-finite matrix or translation"
                )
            }
            Self::TranslationOffGrid { operation } => write!(
                f,
                "operation {operation} has a translation off the exact grid"
            ),
            Self::NotOrthogonal { operation } => write!(
                f,
                "operation {operation} matrix is not orthogonal (M^T M != I)"
            ),
            Self::IdentityMissing => write!(f, "no identity operation in the stored list"),
            Self::IdentityNotIdentity { operation } => write!(
                f,
                "identity operation {operation} has a non-identity matrix"
            ),
            Self::TraceMismatch {
                operation,
                matrix_trace,
                stored_trace,
            } => write!(
                f,
                "operation {operation}: matrix trace {matrix_trace} != stored character {stored_trace}"
            ),
            Self::LawProductNotUnique {
                left,
                right,
                products,
            } => write!(
                f,
                "product of operations {left} and {right} has {products} stored representatives"
            ),
            Self::LawProductAmbiguous {
                left,
                right,
                products,
            } => write!(
                f,
                "product of operations {left} and {right} has {products} stored representatives with different matrices"
            ),
            Self::LawViolated {
                left,
                right,
                product,
            } => write!(
                f,
                "matrices do not multiply like operations: M({left})M({right}) != M({product})"
            ),
            Self::Reducible {
                commutant_dimension,
            } => write!(
                f,
                "the real commutant has dimension {commutant_dimension}, not 1 (reducible)"
            ),
            Self::TypedRowMissing => write!(f, "typed selected-arm Gamma row is unavailable"),
            Self::TypedRowShapeMismatch {
                row,
                operators,
                dimension,
            } => write!(
                f,
                "typed row has {row} entries, dimension {dimension}, but the stored list has {operators} operations"
            ),
            Self::TypedRowOperationMismatch { operation } => write!(
                f,
                "typed row operation {operation} differs from the stored operation"
            ),
            Self::TypedRowMismatch {
                operation,
                matrix_trace,
                row_value,
            } => write!(
                f,
                "operation {operation}: matrix trace {matrix_trace} != typed row value {row_value}"
            ),
            Self::TypedRowNotReal {
                operation,
                imaginary,
            } => write!(
                f,
                "operation {operation}: typed row value has imaginary part {imaginary}"
            ),
        }
    }
}

/// Why a two-dimensional real Γ magnetic matrix subduction was refused.
#[derive(Debug)]
pub enum MatrixSubductionError {
    /// A parent space group number is outside 1–230.
    UnknownParentSpaceGroup { sg: u8 },
    /// The source and the record name different parents.
    ParentMismatch { source_sg: u8, record_sg: u8 },
    /// No irrep of that parent has this label under this convention.
    LabelNotFound {
        sg: u8,
        label: String,
        convention: LabelConvention,
    },
    /// The label exists but does not name a unique irrep.
    AmbiguousLabel { sg: u8, label: String },
    /// The source irrep is double-valued.
    SourceIsSpinor { sg: u8, label: &'static str },
    /// The source irrep is a compound (realified/direct-sum) record.
    SourceIsCompound { sg: u8, label: &'static str },
    /// The source irrep is not two-dimensional.
    SourceDimensionNotTwo {
        sg: u8,
        label: &'static str,
        dimension: u8,
    },
    /// The source irrep's wave vector is not Γ.
    SourceNotGamma {
        sg: u8,
        label: &'static str,
        numerators: [i8; 3],
        denominator: i8,
    },
    /// The record's direction index is outside the condensing irrep's table.
    DirectionOutOfRange {
        sg: u8,
        label: String,
        index: usize,
        len: usize,
    },
    /// A fresh database read disagrees with the record or the operation set.
    DatabaseMismatch { reason: String },
    /// The operation set does not pass the group-axiom proof.
    GroupProof { reason: String },
    /// The class list and the set's own class count disagree.
    ClassListInconsistent { classes: usize, recorded: usize },
    /// The caller's `(map, shift)` does not embed the group.
    EmbeddingRejected(GeometryError),
    /// The selected map does not satisfy the requested lattice policy.
    LatticePolicyUnsatisfied { policy: LatticePolicy },
    /// Exact arithmetic failed inside the checker.
    Arithmetic { reason: String },
    /// The magnetic group has an antiunitary coset but no `η` was supplied.
    MissingTimeReversal { uni: usize },
    /// The parent operation frame could not be built.
    ParentFrameUnavailable { sg: u8, reason: String },
    /// A magnetic class could not be mapped into the parent frame.
    ClassMapping { reason: String },
    /// The source's stored matrix table is not a certified operation-aware real
    /// Γ representation.
    SourceMatrixUncertified {
        sg: u8,
        label: &'static str,
        problem: MatrixCertificationProblem,
    },
    /// A mapped magnetic class has no stored source operation with its rotation
    /// and translation class.
    SourceOperationUnmatched {
        sg: u8,
        label: &'static str,
        rotation: Mat3I,
        time_reversal: bool,
    },
    /// More than one stored source operation matches a mapped class, so the
    /// operation-aware binding is ambiguous.
    SourceOperationAmbiguous {
        sg: u8,
        label: &'static str,
        rotation: Mat3I,
    },
    /// The product of two classes is not a class of the model's group.
    SourceLawUnmatched { left: usize, right: usize },
    /// The certified source action violates the corepresentation law.
    SourceLawViolated {
        left: usize,
        right: usize,
        expected_class: usize,
        expected: SourceMatrixOperator,
        found: SourceMatrixOperator,
    },
    /// The magnetic corep summary of the model's UNI could not be computed.
    SummaryFailed { uni: usize, reason: String },
    /// The summary has no Γ k-point.
    NoGammaKPoint { uni: usize },
    /// A possibly one- or two-dimensional non-spinor Γ corep could not be
    /// constructed, so the candidate list is not known to be complete.
    UnresolvedTargetData {
        uni: usize,
        source: String,
        minimum_dimension: Option<usize>,
        reason: String,
    },
    /// The R27 raw-operation provenance of the summary's Γ columns failed.
    TargetProvenance { uni: usize, reason: String },
    /// The target's source H-irrep could not be resolved as this module
    /// requires.
    TargetSourceUnresolved {
        uni: usize,
        sg: u8,
        label: String,
        reason: String,
    },
    /// The target H-irrep's stored matrix table is not certified.
    TargetMatrixUncertified {
        uni: usize,
        sg: u8,
        label: &'static str,
        problem: MatrixCertificationProblem,
    },
    /// The target's characters, antiunitary construction or multiplication
    /// disagree with the requested corepresentation.
    TargetCorepInconsistent {
        uni: usize,
        label: String,
        reason: String,
    },
    /// An in-shape target H-irrep could not be bound to the magnetic classes.
    TargetOperationUnmatched {
        uni: usize,
        label: String,
        reason: String,
    },
    /// No examined candidate is equivalent to the restriction.
    NoMatchingTarget {
        parent_sg: u8,
        source: &'static str,
        uni: usize,
        examined: Vec<String>,
    },
    /// The matched pieces do not cover the whole source space.
    IncompleteDecomposition {
        source_dimension: u32,
        covered_dimension: u32,
    },
    /// The exhibited copies are not one orthonormal basis of the source space
    /// (`V†V = I` or `Σ U U† = I₂` failed), so a metric-explicit completeness
    /// statement is not available.
    SourceSpaceNotIsometric { covered_dimension: u32 },
}

impl std::fmt::Display for MatrixSubductionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownParentSpaceGroup { sg } => write!(f, "unknown parent {sg}"),
            Self::ParentMismatch {
                source_sg,
                record_sg,
            } => write!(f, "source parent {source_sg} != record parent {record_sg}"),
            Self::LabelNotFound {
                sg,
                label,
                convention,
            } => write!(f, "SG{sg} has no irrep {label:?} under {convention:?}"),
            Self::AmbiguousLabel { sg, label } => write!(f, "SG{sg} irrep {label:?} is ambiguous"),
            Self::SourceIsSpinor { sg, label } => {
                write!(f, "SG{sg} {label} is a spinor source")
            }
            Self::SourceIsCompound { sg, label } => {
                write!(f, "SG{sg} {label} is a compound source")
            }
            Self::SourceDimensionNotTwo {
                sg,
                label,
                dimension,
            } => write!(f, "SG{sg} {label} has dimension {dimension}, not 2"),
            Self::SourceNotGamma {
                sg,
                label,
                numerators,
                denominator,
            } => write!(
                f,
                "SG{sg} {label} is not the Gamma point ({numerators:?}/{denominator})"
            ),
            Self::DirectionOutOfRange {
                sg,
                label,
                index,
                len,
            } => write!(
                f,
                "SG{sg} {label} has no direction index {index} ({len} rows)"
            ),
            Self::DatabaseMismatch { reason } => write!(f, "database mismatch: {reason}"),
            Self::GroupProof { reason } => write!(f, "group proof failed: {reason}"),
            Self::ClassListInconsistent { classes, recorded } => {
                write!(f, "class list has {classes} entries but records {recorded}")
            }
            Self::EmbeddingRejected(error) => write!(f, "embedding rejected: {error}"),
            Self::LatticePolicyUnsatisfied { policy } => {
                write!(f, "lattice policy {policy:?} is not satisfied")
            }
            Self::Arithmetic { reason } => write!(f, "exact arithmetic failed: {reason}"),
            Self::MissingTimeReversal { uni } => {
                write!(f, "UNI {uni} has an antiunitary coset but no eta was given")
            }
            Self::ParentFrameUnavailable { sg, reason } => {
                write!(f, "parent frame of SG{sg} is unavailable: {reason}")
            }
            Self::ClassMapping { reason } => write!(f, "class mapping failed: {reason}"),
            Self::SourceMatrixUncertified { sg, label, problem } => {
                write!(f, "SG{sg} {label} source matrix table: {problem}")
            }
            Self::SourceOperationUnmatched {
                sg,
                label,
                rotation,
                time_reversal,
            } => write!(
                f,
                "SG{sg} {label}: no stored operation for rotation {rotation:?} (time reversal {time_reversal})"
            ),
            Self::SourceOperationAmbiguous {
                sg,
                label,
                rotation,
            } => write!(
                f,
                "SG{sg} {label}: several stored operations match rotation {rotation:?}"
            ),
            Self::SourceLawUnmatched { left, right } => write!(
                f,
                "the product of classes {left} and {right} is not a class"
            ),
            Self::SourceLawViolated {
                left,
                right,
                expected_class,
                expected,
                found,
            } => write!(
                f,
                "source corepresentation law failed at ({left}, {right}): class {expected_class} expects {expected:?}, composition gives {found:?}"
            ),
            Self::SummaryFailed { uni, reason } => {
                write!(f, "UNI {uni} magnetic corep summary failed: {reason}")
            }
            Self::NoGammaKPoint { uni } => write!(f, "UNI {uni} has no Gamma k-point"),
            Self::UnresolvedTargetData {
                uni,
                source,
                minimum_dimension,
                reason,
            } => write!(
                f,
                "UNI {uni} source {source} could not be classified (minimum dimension {minimum_dimension:?}): {reason}"
            ),
            Self::TargetProvenance { uni, reason } => {
                write!(f, "UNI {uni} target operation provenance failed: {reason}")
            }
            Self::TargetSourceUnresolved {
                uni,
                sg,
                label,
                reason,
            } => write!(
                f,
                "UNI {uni} target source SG{sg} {label} is unresolved: {reason}"
            ),
            Self::TargetMatrixUncertified {
                uni,
                sg,
                label,
                problem,
            } => write!(
                f,
                "UNI {uni} target source SG{sg} {label} matrix table: {problem}"
            ),
            Self::TargetCorepInconsistent { uni, label, reason } => {
                write!(f, "UNI {uni} corep {label} is inconsistent: {reason}")
            }
            Self::TargetOperationUnmatched { uni, label, reason } => {
                write!(f, "UNI {uni} corep {label}: {reason}")
            }
            Self::NoMatchingTarget {
                parent_sg,
                source,
                uni,
                examined,
            } => write!(
                f,
                "no Gamma corep of UNI {uni} matches SG{parent_sg} {source} (examined {examined:?})"
            ),
            Self::IncompleteDecomposition {
                source_dimension,
                covered_dimension,
            } => write!(
                f,
                "the matched pieces cover {covered_dimension} of {source_dimension} source dimensions"
            ),
            Self::SourceSpaceNotIsometric { covered_dimension } => write!(
                f,
                "the {covered_dimension} covered dimensions are not one orthonormal copy basis of the source space"
            ),
        }
    }
}

impl std::error::Error for MatrixSubductionError {}

impl From<GeometryError> for MatrixSubductionError {
    fn from(error: GeometryError) -> Self {
        Self::EmbeddingRejected(error)
    }
}

// ---------------------------------------------------------------------------
// Request
// ---------------------------------------------------------------------------

/// Everything a caller must state to obtain a verified two-dimensional real Γ
/// matrix model.
///
/// The identities, the embedding evidence and `η` are the *same* public types
/// the R25 scalar contract uses, so a caller's provenance vocabulary does not
/// fork; only the source dimension and the output differ.
pub struct MatrixGammaRequest<'a> {
    /// The source parent irrep identity.
    pub source: SourceIdentity<'a>,
    /// The condensing magnetic isotropy record identity.
    pub record: RecordIdentity<'a>,
    /// The caller's explicit embedding.
    pub embedding: EmbeddingEvidence,
    /// The parent action's `η`, required when the magnetic group has an
    /// antiunitary coset.
    pub time_reversal: Option<RealSign>,
}

impl<'a> MatrixGammaRequest<'a> {
    /// A request for one two-dimensional source and one condensing record.
    pub fn new(
        source: SourceIdentity<'a>,
        record: RecordIdentity<'a>,
        time_reversal: Option<RealSign>,
    ) -> Self {
        Self {
            source,
            record,
            embedding: EmbeddingEvidence::identity(),
            time_reversal,
        }
    }

    /// The same request with a different embedding.
    pub fn with_embedding(mut self, embedding: EmbeddingEvidence) -> Self {
        self.embedding = embedding;
        self
    }
}

// ---------------------------------------------------------------------------
// Operators
// ---------------------------------------------------------------------------

/// One certified `2×2` real source action, with its unitarity branch.
///
/// A unitary operator acts as `v ↦ M v`; an antiunitary one as
/// `v ↦ M conj(v)`, exactly the parent grey action `D(p) η^{t} K^{t}`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SourceMatrixOperator {
    entries: [f64; 4],
    antiunitary: bool,
}

impl SourceMatrixOperator {
    /// Row-major matrix entry, or `None` outside the `2×2` block.
    pub fn entry(&self, row: usize, column: usize) -> Option<f64> {
        (row < 2 && column < 2).then(|| self.entries[row * 2 + column])
    }

    /// All four row-major entries.
    pub fn entries(&self) -> [f64; 4] {
        self.entries
    }

    /// Whether the class is antiunitary.
    pub fn is_antiunitary(&self) -> bool {
        self.antiunitary
    }

    /// Apply the operator to a two-component complex vector.
    pub fn apply(&self, vector: [Complex64; 2]) -> [Complex64; 2] {
        let source = if self.antiunitary {
            [vector[0].conj(), vector[1].conj()]
        } else {
            vector
        };
        [
            Complex64::new(self.entries[0], 0.0) * source[0]
                + Complex64::new(self.entries[1], 0.0) * source[1],
            Complex64::new(self.entries[2], 0.0) * source[0]
                + Complex64::new(self.entries[3], 0.0) * source[1],
        ]
    }
}

/// One reconstructed `d×d` real target operator on the target space.
#[derive(Debug, Clone, PartialEq)]
pub struct TargetMatrixOperator {
    dimension: usize,
    entries: Vec<f64>,
    antiunitary: bool,
}

impl TargetMatrixOperator {
    /// Dimension of the target operator.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Row-major matrix entry, or `None` outside the block.
    pub fn entry(&self, row: usize, column: usize) -> Option<f64> {
        (row < self.dimension && column < self.dimension)
            .then(|| self.entries[row * self.dimension + column])
    }

    /// The row-major matrix entries.
    pub fn entries(&self) -> &[f64] {
        &self.entries
    }

    /// Whether the class is antiunitary.
    pub fn is_antiunitary(&self) -> bool {
        self.antiunitary
    }

    /// Apply the operator to a `d`-component complex vector.
    pub fn apply(&self, vector: &[Complex64]) -> Vec<Complex64> {
        let source: Vec<Complex64> = if self.antiunitary {
            vector.iter().map(|value| value.conj()).collect()
        } else {
            vector.to_vec()
        };
        (0..self.dimension)
            .map(|row| {
                (0..self.dimension)
                    .map(|column| {
                        Complex64::new(self.entries[row * self.dimension + column], 0.0)
                            * source[column]
                    })
                    .sum()
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// The model
// ---------------------------------------------------------------------------

/// A verified two-dimensional real Γ magnetic model: the parent source action on
/// the classes of one magnetic group, with the embedding and both identities.
///
/// Fields are private and the only construction path is
/// [`verify_matrix_gamma_input`], so a successful object cannot be forged from
/// outside; the report borrows it instead of copying it.
#[derive(Debug)]
pub struct MatrixGammaModel {
    parent_sg: u8,
    source_label: &'static str,
    source_convention: LabelConvention,
    input_convention: LabelConvention,
    record_label: &'static str,
    record_convention: LabelConvention,
    direction_index: usize,
    ordinal: usize,
    uni: usize,
    map: crate::irrep::subduction::Mat3R,
    shift: [Rat; 3],
    lattice_policy: LatticePolicy,
    realises_record_lattice: bool,
    time_reversal: Option<RealSign>,
    actions: Vec<SourceMatrixOperator>,
    classes: Vec<MagneticOperation>,
    lattice: Vec<[Rat; 3]>,
}

impl MatrixGammaModel {
    /// Parent space group of the source and the record.
    pub fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// The source parent irrep label, in its canonical CDML spelling.
    pub fn source_label(&self) -> &'static str {
        self.source_label
    }

    /// The convention [`Self::source_label`] is written in (`Cdml`).
    pub fn source_convention(&self) -> LabelConvention {
        self.source_convention
    }

    /// The convention the caller used to name the source label.
    pub fn input_convention(&self) -> LabelConvention {
        self.input_convention
    }

    /// The condensing irrep label of the record, in its canonical CDML spelling.
    pub fn record_label(&self) -> &'static str {
        self.record_label
    }

    /// The convention [`Self::record_label`] is written in (`Cdml`).
    pub fn record_convention(&self) -> LabelConvention {
        self.record_convention
    }

    /// The record's row index in the condensing irrep's magnetic isotropy table.
    pub fn direction_index(&self) -> usize {
        self.direction_index
    }

    /// The record's ordinal in the generated magnetic isotropy table.
    pub fn ordinal(&self) -> usize {
        self.ordinal
    }

    /// The record's magnetic space group (UNI) number.
    pub fn uni(&self) -> usize {
        self.uni
    }

    /// The verified map of the embedding.
    pub fn map(&self) -> &crate::irrep::subduction::Mat3R {
        &self.map
    }

    /// The verified origin shift of the embedding.
    pub fn shift(&self) -> &[Rat; 3] {
        &self.shift
    }

    /// The lattice policy the model was checked against.
    pub fn lattice_policy(&self) -> LatticePolicy {
        self.lattice_policy
    }

    /// Whether the selected map realises the record's lattice exactly.
    pub fn realises_record_lattice(&self) -> bool {
        self.realises_record_lattice
    }

    /// The parent action's `η`, or `None` when the caller stated none.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.time_reversal
    }

    /// The source dimension, always two for this entry point.
    pub fn source_dimension(&self) -> u32 {
        SOURCE_DIMENSION as u32
    }

    /// The certified source action of every class, in the set's class order.
    pub fn source_operators(&self) -> &[SourceMatrixOperator] {
        &self.actions
    }

    /// The certified source action of one class, by class index.
    pub fn source_operator_at(&self, index: usize) -> Option<&SourceMatrixOperator> {
        self.actions.get(index)
    }

    /// The classes of the magnetic group, in the set's class order.
    pub fn classes(&self) -> &[MagneticOperation] {
        &self.classes
    }

    /// The group's translation lattice, as the class comparisons reduced it.
    pub fn translation_lattice(&self) -> &[[Rat; 3]] {
        &self.lattice
    }

    /// One line stating what this entry point covers.
    pub fn scope(&self) -> &'static str {
        SCOPE
    }
}

// ---------------------------------------------------------------------------
// The result
// ---------------------------------------------------------------------------

/// One target Γ corepresentation of the magnetic group, with the independent
/// intertwiner copies that realise it inside the two-dimensional source.
#[derive(Debug, Clone)]
pub struct MatrixSubductionPiece {
    label: String,
    source_label: &'static str,
    unitary_sg: u8,
    dimension: u32,
    multiplicity: u32,
    copies: Vec<Vec<Complex64>>,
    operators: Vec<TargetMatrixOperator>,
}

impl MatrixSubductionPiece {
    /// The target corep label as the magnetic summary names it.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// The unitary subgroup's ordinary irrep label the target is induced from.
    pub fn source_label(&self) -> &'static str {
        self.source_label
    }

    /// Space group number of the target's unitary subgroup `H`.
    pub fn unitary_sg(&self) -> u8 {
        self.unitary_sg
    }

    /// Dimension of this target corepresentation.
    pub fn dimension(&self) -> u32 {
        self.dimension
    }

    /// Multiplicity of this target in the restriction.
    pub fn multiplicity(&self) -> u32 {
        self.multiplicity
    }

    /// `multiplicity · dimension`, the source dimensions this piece covers.
    pub fn covered_dimension(&self) -> u32 {
        self.multiplicity * self.dimension
    }

    /// The independently normalized `2×d` complex intertwiners, `multiplicity`
    /// of them, each row-major with `2·d` entries.
    pub fn copies(&self) -> &[Vec<Complex64>] {
        &self.copies
    }

    /// The reconstructed target operator of every class of `M`, in the model's
    /// class order.
    pub fn operators(&self) -> &[TargetMatrixOperator] {
        &self.operators
    }

    /// The reconstructed target operator of one class.
    pub fn operator_at(&self, index: usize) -> Option<&TargetMatrixOperator> {
        self.operators.get(index)
    }
}

/// The result of one two-dimensional real Γ matrix subduction.
///
/// The fields are private and the only construction path is
/// [`subduce_matrix_gamma`].  The report **borrows** the verified model, so the
/// embedding, both identities, `η` and the certified `2×2` action table stay the
/// model's own data.
#[derive(Debug)]
pub struct MatrixSubduction<'a> {
    model: &'a MatrixGammaModel,
    pieces: Vec<MatrixSubductionPiece>,
    source_dimension: u32,
    covered_dimension: u32,
}

impl MatrixSubduction<'_> {
    /// The verified input model this report is about.
    pub fn model(&self) -> &MatrixGammaModel {
        self.model
    }

    /// Dimension of the subduced source corepresentation (two).
    pub fn source_dimension(&self) -> u32 {
        self.source_dimension
    }

    /// `Σ multiplicityᵢ · dimensionᵢ`, checked equal to `source_dimension`.
    pub fn covered_dimension(&self) -> u32 {
        self.covered_dimension
    }

    /// The matched pieces, in the summary's Γ candidate order.
    pub fn pieces(&self) -> &[MatrixSubductionPiece] {
        &self.pieces
    }

    /// One matched piece.
    pub fn piece(&self, index: usize) -> Option<&MatrixSubductionPiece> {
        self.pieces.get(index)
    }

    /// One line stating what this entry point covers.
    pub fn scope(&self) -> &'static str {
        SCOPE
    }
}

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// A stored matrix table certified as an operation-aware real Γ representation.
#[derive(Debug, Clone)]
struct CertifiedTable {
    dimension: usize,
    operations: Vec<(Mat3I, [Rat; 3])>,
    matrices: Vec<Vec<f64>>,
}

impl CertifiedTable {
    fn dimension(&self) -> usize {
        self.dimension
    }
}

/// Exact small-denominator conversion of one stored `f64`.
fn exact_fraction(value: f64) -> Option<Rat> {
    if !value.is_finite() {
        return None;
    }
    for denominator in 1..=MAX_TRANSLATION_DENOMINATOR {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < GRID_TOLERANCE {
            return Rat::new(scaled.round() as i128, denominator).ok();
        }
    }
    None
}

/// The identity rotation.
fn identity_rotation() -> Mat3I {
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
}

/// Row-major `d×d` matrix product of two `d²` slices.
fn matrix_product(left: &[f64], right: &[f64], dimension: usize) -> Vec<f64> {
    let mut out = vec![0.0; dimension * dimension];
    for row in 0..dimension {
        for column in 0..dimension {
            out[row * dimension + column] = (0..dimension)
                .map(|k| left[row * dimension + k] * right[k * dimension + column])
                .sum();
        }
    }
    out
}

/// Trace of a row-major `d×d` matrix.
fn matrix_trace(entries: &[f64], dimension: usize) -> f64 {
    (0..dimension)
        .map(|index| entries[index * dimension + index])
        .sum()
}

/// Whether `Mᵀ M = I`.
/// Positive finite closeness: both sides finite and within tolerance.
///
/// A `NaN` or infinite value is never "close" (a negation-style `> tol` test
/// would silently accept it), and an exact zero residual passes.
fn close(left: f64, right: f64) -> bool {
    left.is_finite() && right.is_finite() && (left - right).abs() <= GRID_TOLERANCE
}

/// Positive finite closeness of two row-major matrices.
fn matrices_close(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len() && left.iter().zip(right).all(|(a, b)| close(*a, *b))
}

/// Positive finite "small": a finite norm at most `tolerance`.
fn small(value: f64, tolerance: f64) -> bool {
    value.is_finite() && value <= tolerance
}

fn is_orthogonal(entries: &[f64], dimension: usize) -> bool {
    for row in 0..dimension {
        for column in 0..dimension {
            let value: f64 = (0..dimension)
                .map(|k| entries[k * dimension + row] * entries[k * dimension + column])
                .sum();
            let expected = if row == column { 1.0 } else { 0.0 };
            if !close(value, expected) {
                return false;
            }
        }
    }
    true
}

/// A basis of the real nullspace of `rows` (width columns).
fn real_kernel(rows: &[Vec<f64>], width: usize) -> Vec<Vec<f64>> {
    let mut reduced: Vec<Vec<f64>> = rows.to_vec();
    let mut pivots: Vec<usize> = Vec::new();
    let mut pivot_row = 0usize;
    for column in 0..width {
        if pivot_row >= reduced.len() {
            break;
        }
        let Some(found) = (pivot_row..reduced.len()).find(|&row| reduced[row][column].abs() > 1e-9)
        else {
            continue;
        };
        reduced.swap(pivot_row, found);
        let scale = 1.0 / reduced[pivot_row][column];
        for value in reduced[pivot_row][column..].iter_mut() {
            *value *= scale;
        }
        for row in 0..reduced.len() {
            if row == pivot_row {
                continue;
            }
            let factor = reduced[row][column];
            if factor.abs() <= 1e-12 {
                continue;
            }
            let pivot_values: Vec<f64> = reduced[pivot_row][column..].to_vec();
            for (offset, value) in pivot_values.into_iter().enumerate() {
                reduced[row][column + offset] -= factor * value;
            }
        }
        pivots.push(column);
        pivot_row += 1;
    }
    let free: Vec<usize> = (0..width)
        .filter(|column| !pivots.contains(column))
        .collect();
    let mut basis = Vec::with_capacity(free.len());
    for free_column in free {
        let mut vector = vec![0.0; width];
        vector[free_column] = 1.0;
        for (index, pivot) in pivots.iter().enumerate() {
            vector[*pivot] = -reduced[index][free_column];
        }
        basis.push(vector);
    }
    basis
}

/// The dimension of the **real** commutant of a real matrix table: the space of
/// real `d×d` matrices `X` with `X D = D X` for every stored matrix.
fn real_commutant_dimension(matrices: &[Vec<f64>], dimension: usize) -> usize {
    let width = dimension * dimension;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for matrix in matrices {
        for row in 0..dimension {
            for column in 0..dimension {
                let mut equation = vec![0.0; width];
                for k in 0..dimension {
                    equation[row * dimension + k] += matrix[k * dimension + column];
                    equation[k * dimension + column] -= matrix[row * dimension + k];
                }
                rows.push(equation);
            }
        }
    }
    real_kernel(&rows, width).len()
}

/// The positive scale `c` with `(c·base)² = square`.
///
/// `base` is the kernel matrix and `square = D(a0²)`.  The kernel satisfies
/// `square = λ·base²` with `λ = square/base²` **positive** and finite (the
/// non-unit fixture has `λ = 1/4`), so `c = √λ` and the rescaled coset matrix is
/// `c·base = √λ·base` — dividing by `√λ` would invert the correction (with
/// `λ = 1/4` it gives `2·base`, whose square is `16·square`).  The ratio is
/// checked entry by entry on the whole `d×d` square, and the rescaled matrix's
/// square must reproduce `square` entry by entry.
fn coset_scale(base_square: &[f64], square: &[f64], dimension: usize) -> Option<f64> {
    if base_square.len() != dimension * dimension || square.len() != dimension * dimension {
        return None;
    }
    let mut scale: Option<f64> = None;
    for index in 0..base_square.len() {
        if !base_square[index].is_finite() || !square[index].is_finite() {
            return None;
        }
        if base_square[index].abs() <= GRID_TOLERANCE {
            // A vanishing entry of base² must have a vanishing square entry.
            if square[index].abs() > GRID_TOLERANCE {
                return None;
            }
            continue;
        }
        let ratio = square[index] / base_square[index];
        if !ratio.is_finite() || ratio <= 0.0 {
            return None;
        }
        match scale {
            None => scale = Some(ratio),
            Some(known) if close(known, ratio) => {}
            Some(_) => return None,
        }
    }
    let scale = scale?;
    let magnitude = scale.sqrt();
    if !magnitude.is_finite() || magnitude <= 0.0 {
        return None;
    }
    // The whole rescaled square must reproduce D(a0²); `scaled[k] = c²·base²[k]`.
    for index in 0..base_square.len() {
        if !close(magnitude * magnitude * base_square[index], square[index]) {
            return None;
        }
    }
    Some(magnitude)
}

/// The real corepresentation of one `d`-dimensional table, as the intertwiner
/// solver sees it.
#[derive(Debug, Clone)]
struct RealCorep {
    dimension: usize,
    matrices: Vec<Vec<f64>>,
    antiunitary: Vec<bool>,
}

impl RealCorep {
    fn from_source(operators: &[SourceMatrixOperator]) -> Self {
        Self {
            dimension: SOURCE_DIMENSION,
            matrices: operators
                .iter()
                .map(|operator| operator.entries.to_vec())
                .collect(),
            antiunitary: operators
                .iter()
                .map(|operator| operator.antiunitary)
                .collect(),
        }
    }

    fn from_target(operators: &[TargetMatrixOperator]) -> Self {
        Self {
            dimension: operators[0].dimension,
            matrices: operators
                .iter()
                .map(|operator| operator.entries.clone())
                .collect(),
            antiunitary: operators
                .iter()
                .map(|operator| operator.antiunitary)
                .collect(),
        }
    }
}

/// The real-linear system whose nullspace is `Hom(left, right)`.
///
/// The unknown is a complex `left_dim × right_dim` matrix `U`, laid out as
/// `(row·right_dim + column)·2 + {0 = re, 1 = im}`.  A unitary class gives
/// `L V = V R`; an antiunitary class gives `L conj(V) = V R`, which only flips
/// the sign of the imaginary coupling of the left factor.  That conjugation is
/// what makes the solution space real-linear rather than complex-linear.
fn hom_rows(left: &RealCorep, right: &RealCorep) -> (Vec<Vec<f64>>, usize) {
    let rows_dim = left.dimension;
    let columns_dim = right.dimension;
    let width = 2 * rows_dim * columns_dim;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for class in 0..left.matrices.len() {
        let l = &left.matrices[class];
        let r = &right.matrices[class];
        let sign = if left.antiunitary[class] { -1.0 } else { 1.0 };
        for row in 0..rows_dim {
            for column in 0..columns_dim {
                let mut re_row = vec![0.0; width];
                let mut im_row = vec![0.0; width];
                for k in 0..rows_dim {
                    re_row[(k * columns_dim + column) * 2] += l[row * rows_dim + k];
                    im_row[(k * columns_dim + column) * 2 + 1] += sign * l[row * rows_dim + k];
                }
                for j in 0..columns_dim {
                    re_row[(row * columns_dim + j) * 2] -= r[j * columns_dim + column];
                    im_row[(row * columns_dim + j) * 2 + 1] -= r[j * columns_dim + column];
                }
                rows.push(re_row);
                rows.push(im_row);
            }
        }
    }
    (rows, width)
}

/// Rebuild a complex `rows × columns` matrix from one nullspace vector.
fn intertwiner_from_vector(vector: &[f64], columns: usize) -> Vec<Complex64> {
    (0..(vector.len() / 2))
        .map(|entry| Complex64::new(vector[2 * entry], vector[2 * entry + 1]))
        .collect::<Vec<_>>()
        .chunks(columns)
        .flat_map(|chunk| chunk.to_vec())
        .collect()
}

/// The certified matrix table of one Γ irrep record.
fn certify_irrep_matrices(
    irrep: &IrrepRecord,
    expected_dimension: usize,
) -> Result<CertifiedTable, MatrixCertificationProblem> {
    let declared = irrep.dim as usize;
    if declared != expected_dimension {
        return Err(MatrixCertificationProblem::DeclaredDimensionMismatch {
            declared,
            expected: expected_dimension,
        });
    }
    if expected_dimension == 0 || expected_dimension > 2 {
        return Err(MatrixCertificationProblem::DimensionOutOfScope {
            dimension: expected_dimension,
        });
    }
    let rotations = irrep.pir_rotations();
    let translations = irrep.pir_translations();
    let matrices = irrep.matrices();
    let characters = irrep.characters();
    if !rotations.len().is_multiple_of(9) {
        return Err(MatrixCertificationProblem::ShapeMismatch {
            operators: rotations.len() / 9,
            matrices: matrices.len(),
            traces: characters.len(),
            dimension: expected_dimension,
        });
    }
    let operators = rotations.len() / 9;
    let entries = expected_dimension * expected_dimension;
    if operators == 0
        || translations.len() != operators * 3
        || characters.len() != operators
        || matrices.len() != operators * entries
    {
        return Err(MatrixCertificationProblem::ShapeMismatch {
            operators,
            matrices: matrices.len(),
            traces: characters.len(),
            dimension: expected_dimension,
        });
    }

    // Exact operation table and its own translation lattice.
    let mut operations: Vec<(Mat3I, [Rat; 3])> = Vec::with_capacity(operators);
    let mut lattice: Vec<[Rat; 3]> = Vec::new();
    for index in 0..operators {
        let cell = &rotations[index * 9..index * 9 + 9];
        let rotation = [
            [cell[0], cell[1], cell[2]],
            [cell[3], cell[4], cell[5]],
            [cell[6], cell[7], cell[8]],
        ];
        let mut translation = [Rat::ZERO; 3];
        for axis in 0..3 {
            let value = translations[index * 3 + axis];
            if !value.is_finite() {
                return Err(MatrixCertificationProblem::NonFinite { operation: index });
            }
            translation[axis] = exact_fraction(value)
                .ok_or(MatrixCertificationProblem::TranslationOffGrid { operation: index })?;
        }
        if rotation == identity_rotation()
            && translation.iter().any(|value| !value.is_zero())
            && !lattice.contains(&translation)
        {
            lattice.push(translation);
        }
        operations.push((rotation, translation));
    }

    let mut table = Vec::with_capacity(operators);
    let mut traces = Vec::with_capacity(operators);
    for index in 0..operators {
        let matrix: Vec<f64> = matrices[index * entries..(index + 1) * entries].to_vec();
        if matrix.iter().any(|value| !value.is_finite()) {
            return Err(MatrixCertificationProblem::NonFinite { operation: index });
        }
        if !is_orthogonal(&matrix, expected_dimension) {
            return Err(MatrixCertificationProblem::NotOrthogonal { operation: index });
        }
        let trace = matrix_trace(&matrix, expected_dimension);
        if !close(trace, characters[index]) {
            return Err(MatrixCertificationProblem::TraceMismatch {
                operation: index,
                matrix_trace: trace,
                stored_trace: characters[index],
            });
        }
        traces.push(trace);
        table.push(matrix);
    }

    // The identity must be present and act as the identity.  A centred group
    // lists redundant rows per coset; at Γ they all denote the same group
    // element, so the identity may appear more than once as long as every such
    // row carries the identity matrix.
    let identities: Vec<usize> = operations
        .iter()
        .enumerate()
        .filter(|(_, (rotation, translation))| {
            *rotation == identity_rotation() && translation.iter().all(|value| value.is_zero())
        })
        .map(|(index, _)| index)
        .collect();
    let Some(&identity) = identities.first() else {
        return Err(MatrixCertificationProblem::IdentityMissing);
    };
    let identity_matrix: Vec<f64> = (0..entries)
        .map(|index| {
            if index % (expected_dimension + 1) == 0 {
                1.0
            } else {
                0.0
            }
        })
        .collect();
    if identities
        .iter()
        .any(|&index| table[index] != identity_matrix)
    {
        return Err(MatrixCertificationProblem::IdentityNotIdentity {
            operation: identity,
        });
    }

    // Complete ordered-pair group law on the stored operation list, modulo its
    // own translation lattice.  Redundant rows of one coset are allowed, but
    // they must carry the same matrix at Γ — otherwise the storage is not a
    // well-defined function of the group element.
    check_stored_matrix_law(&operations, &table, expected_dimension, &lattice)?;

    // Irreducibility over the reals: the real commutant is one-dimensional.
    // (`hom_rows` treats the real and imaginary parts of a *complex* unknown as
    // independent, so it deliberately does not serve as the real commutant.)
    let commutant = real_commutant_dimension(&table, expected_dimension);
    if commutant != 1 {
        return Err(MatrixCertificationProblem::Reducible {
            commutant_dimension: commutant,
        });
    }

    // Cross-check against the typed selected-arm Γ row the production corep
    // path reads, so the certified traces do not rest on one storage array.
    let row = irrep
        .ordinary_scalar_selected_arm_block_trace()
        .map_err(|_| MatrixCertificationProblem::TypedRowMissing)?;
    if row.len() != operators || row.dimension() != expected_dimension {
        return Err(MatrixCertificationProblem::TypedRowShapeMismatch {
            row: row.len(),
            operators,
            dimension: row.dimension(),
        });
    }
    for index in 0..operators {
        let operation = row.operations()[index];
        let flat = [
            operation.rotation[0],
            operation.rotation[1],
            operation.rotation[2],
            operation.rotation[3],
            operation.rotation[4],
            operation.rotation[5],
            operation.rotation[6],
            operation.rotation[7],
            operation.rotation[8],
        ];
        if [
            [flat[0], flat[1], flat[2]],
            [flat[3], flat[4], flat[5]],
            [flat[6], flat[7], flat[8]],
        ] != operations[index].0
        {
            return Err(MatrixCertificationProblem::TypedRowOperationMismatch { operation: index });
        }
        let value = row.values()[index];
        if !close(value.im, 0.0) {
            return Err(MatrixCertificationProblem::TypedRowNotReal {
                operation: index,
                imaginary: value.im,
            });
        }
        if !close(value.re, traces[index]) {
            return Err(MatrixCertificationProblem::TypedRowMismatch {
                operation: index,
                matrix_trace: traces[index],
                row_value: value.re,
            });
        }
    }

    Ok(CertifiedTable {
        dimension: expected_dimension,
        operations,
        matrices: table,
    })
}

/// The complete ordered-pair matrix law of one stored operation list.
///
/// This is the discriminating check of the certification: a permutation of two
/// operation rows that keeps every trace (and the commutation structure) still
/// fails here, because the stored matrices then do not multiply like the stored
/// operations.
fn check_stored_matrix_law(
    operations: &[(Mat3I, [Rat; 3])],
    matrices: &[Vec<f64>],
    dimension: usize,
    lattice: &[[Rat; 3]],
) -> Result<(), MatrixCertificationProblem> {
    let entries = dimension * dimension;
    let operators = operations.len();
    for left in 0..operators {
        for right in 0..operators {
            let product = compose(
                &operation_of(&operations[left]),
                &operation_of(&operations[right]),
            )
            .map_err(|_| MatrixCertificationProblem::LawProductNotUnique {
                left,
                right,
                products: 0,
            })?;
            let matches: Vec<usize> = (0..operators)
                .filter(|&index| {
                    operations[index].0 == product.rotation
                        && differs_by_lattice(&product.translation, &operations[index].1, lattice)
                })
                .collect();
            if matches.is_empty() {
                return Err(MatrixCertificationProblem::LawProductNotUnique {
                    left,
                    right,
                    products: matches.len(),
                });
            }
            let expected = &matrices[matches[0]];
            if matches
                .iter()
                .any(|&index| !matrices_close(&matrices[index], expected))
            {
                return Err(MatrixCertificationProblem::LawProductAmbiguous {
                    left,
                    right,
                    products: matches.len(),
                });
            }
            let found = matrix_product(&matrices[left], &matrices[right], dimension);
            if (0..entries).any(|index| !close(found[index], expected[index])) {
                return Err(MatrixCertificationProblem::LawViolated {
                    left,
                    right,
                    product: matches[0],
                });
            }
        }
    }
    Ok(())
}

/// A listed operation as the magnetic composition sees it.
fn operation_of(operation: &(Mat3I, [Rat; 3])) -> MagneticOperation {
    MagneticOperation {
        rotation: operation.0,
        translation: operation.1,
        time_reversal: false,
    }
}

/// [`crate::irrep::magnetic_input`]'s lattice-aware translation equality.
fn differs_by_lattice(left: &[Rat; 3], right: &[Rat; 3], lattice: &[[Rat; 3]]) -> bool {
    magnetic_input::differs_by_lattice(left, right, lattice)
}

// ---------------------------------------------------------------------------
// The input contract
// ---------------------------------------------------------------------------

/// Verify one two-dimensional real Γ source and its embedding, certify the
/// source's operation-aware matrices, and build the parent grey action on the
/// magnetic group's classes.
pub fn verify_matrix_gamma_input(
    request: &MatrixGammaRequest<'_>,
) -> Result<MatrixGammaModel, MatrixSubductionError> {
    let source_identity = request.source;
    let record_identity = request.record;
    if source_identity.parent_sg == 0 || source_identity.parent_sg > 230 {
        return Err(MatrixSubductionError::UnknownParentSpaceGroup {
            sg: source_identity.parent_sg,
        });
    }
    if source_identity.parent_sg != record_identity.parent_sg {
        return Err(MatrixSubductionError::ParentMismatch {
            source_sg: source_identity.parent_sg,
            record_sg: record_identity.parent_sg,
        });
    }
    let parent_sg = source_identity.parent_sg;
    let source_irrep = resolve_matrix_source(source_identity)?;

    let condensing = resolve_convention_irrep(
        record_identity.parent_sg,
        record_identity.condensing_label,
        record_identity.convention,
    )?;
    let rows = crate::irrep::isotropy::magnetic_isotropy_subgroups(
        parent_sg,
        condensing.ml,
        LabelConvention::Cdml,
    )
    .map_err(|error| MatrixSubductionError::DatabaseMismatch {
        reason: error.to_string(),
    })?;
    let Some(subgroup) = rows.get(record_identity.direction_index).copied() else {
        return Err(MatrixSubductionError::DirectionOutOfRange {
            sg: parent_sg,
            label: condensing.ml.to_string(),
            index: record_identity.direction_index,
            len: rows.len(),
        });
    };
    let subgroup_record = subgroup.record;
    let geometry = magnetic_embedding::geometry_of(parent_sg, &subgroup_record)?;

    let set = magnetic_operations(subgroup_record.mag_sg).map_err(|error| {
        MatrixSubductionError::DatabaseMismatch {
            reason: error.to_string(),
        }
    })?;
    magnetic_embedding::verify_group(&set).map_err(|error| MatrixSubductionError::GroupProof {
        reason: error.to_string(),
    })?;
    magnetic_embedding::verify_against_database(&set).map_err(|error| {
        MatrixSubductionError::DatabaseMismatch {
            reason: error.to_string(),
        }
    })?;
    if set.class_count() != set.classes.len() {
        return Err(MatrixSubductionError::ClassListInconsistent {
            classes: set.classes.len(),
            recorded: set.class_count(),
        });
    }

    let evidence = request.embedding;
    if evidence.map.inverse().is_err() {
        return Err(MatrixSubductionError::Arithmetic {
            reason: "the selected map is singular".to_string(),
        });
    }
    magnetic_embedding::verify_embedding(&geometry, &set, evidence.map, evidence.shift)
        .map_err(MatrixSubductionError::EmbeddingRejected)?;
    let realises_record_lattice =
        magnetic_embedding::realises_record_lattice_checked(&geometry, &set, &evidence.map)?;
    if evidence.lattice_policy == LatticePolicy::RequireExact && !realises_record_lattice {
        return Err(MatrixSubductionError::LatticePolicyUnsatisfied {
            policy: evidence.lattice_policy,
        });
    }

    let time_reversal = match (set.antiunitary_classes, request.time_reversal) {
        (0, given) => given,
        (_, Some(sign)) => Some(sign),
        (_, None) => {
            return Err(MatrixSubductionError::MissingTimeReversal { uni: set.uni });
        }
    };

    // The source's own stored matrices, certified before any of them is used.
    let certified = certify_irrep_matrices(source_irrep, SOURCE_DIMENSION).map_err(|problem| {
        MatrixSubductionError::SourceMatrixUncertified {
            sg: parent_sg,
            label: source_irrep.ml,
            problem,
        }
    })?;
    if certified.dimension() != SOURCE_DIMENSION {
        return Err(MatrixSubductionError::SourceDimensionNotTwo {
            sg: parent_sg,
            label: source_irrep.ml,
            dimension: certified.dimension() as u8,
        });
    }

    let parent = ParentFrame::of(parent_sg).map_err(|error| {
        MatrixSubductionError::ParentFrameUnavailable {
            sg: parent_sg,
            reason: format!("{error:?}"),
        }
    })?;
    let mapped = magnetic_input::mapped_classes(&set, &evidence, &parent).map_err(|error| {
        MatrixSubductionError::ClassMapping {
            reason: format!("{error:?}"),
        }
    })?;

    let eta = time_reversal.map(|sign| sign.value() as f64).unwrap_or(1.0);
    let mut actions = Vec::with_capacity(mapped.len());
    for mapped_class in &mapped {
        let index = bind_source_operation(
            &certified,
            source_irrep.ml,
            parent_sg,
            mapped_class,
            &parent,
        )?;
        let mut entries = [
            certified.matrices[index][0],
            certified.matrices[index][1],
            certified.matrices[index][2],
            certified.matrices[index][3],
        ];
        if mapped_class.time_reversal {
            for value in entries.iter_mut() {
                *value *= eta;
            }
        }
        actions.push(SourceMatrixOperator {
            entries,
            antiunitary: mapped_class.time_reversal,
        });
    }

    let classes = set.classes.clone();
    let lattice = magnetic_embedding::translation_lattice(&set);
    check_source_law(&classes, &actions, &lattice)?;

    Ok(MatrixGammaModel {
        parent_sg,
        source_label: source_irrep.ml,
        source_convention: LabelConvention::Cdml,
        input_convention: source_identity.convention,
        record_label: condensing.ml,
        record_convention: LabelConvention::Cdml,
        direction_index: record_identity.direction_index,
        ordinal: subgroup.ordinal,
        uni: set.uni,
        map: evidence.map,
        shift: evidence.shift,
        lattice_policy: evidence.lattice_policy,
        realises_record_lattice,
        time_reversal,
        actions,
        classes,
        lattice,
    })
}

/// Resolve one label in one convention, refusing a missing or ambiguous one.
fn resolve_convention_irrep(
    sg: u8,
    label: &str,
    convention: LabelConvention,
) -> Result<&'static IrrepRecord, MatrixSubductionError> {
    let candidates: Vec<&'static IrrepRecord> = query::irreps_of(sg)
        .iter()
        .filter(|irrep| irrep.label(convention).as_deref() == Some(label))
        .collect();
    match candidates.len() {
        0 => Err(MatrixSubductionError::LabelNotFound {
            sg,
            label: label.to_string(),
            convention,
        }),
        1 => Ok(candidates[0]),
        _ => Err(MatrixSubductionError::AmbiguousLabel {
            sg,
            label: label.to_string(),
        }),
    }
}

/// Resolve the two-dimensional real non-compound Γ source, refusing the rest of
/// the parent's irreps by name.
fn resolve_matrix_source(
    source: SourceIdentity<'_>,
) -> Result<&'static IrrepRecord, MatrixSubductionError> {
    let irrep = resolve_convention_irrep(source.parent_sg, source.label, source.convention)?;
    if irrep.spinor {
        return Err(MatrixSubductionError::SourceIsSpinor {
            sg: source.parent_sg,
            label: irrep.ml,
        });
    }
    if irrep.compound_metadata().is_some() {
        return Err(MatrixSubductionError::SourceIsCompound {
            sg: source.parent_sg,
            label: irrep.ml,
        });
    }
    if irrep.dim != SOURCE_DIMENSION as u8 {
        return Err(MatrixSubductionError::SourceDimensionNotTwo {
            sg: source.parent_sg,
            label: irrep.ml,
            dimension: irrep.dim,
        });
    }
    let k = irrep.k_vector();
    if k.numerators != [0, 0, 0] || k.denominator != 1 {
        return Err(MatrixSubductionError::SourceNotGamma {
            sg: source.parent_sg,
            label: irrep.ml,
            numerators: k.numerators,
            denominator: k.denominator,
        });
    }
    Ok(irrep)
}

/// Bind one mapped magnetic class to its unique stored source operation row.
///
/// A centred source lists several rows per coset; they must all carry the same
/// matrix at Γ, and then any of them names the same group element.
fn bind_source_operation(
    certified: &CertifiedTable,
    label: &'static str,
    sg: u8,
    mapped: &MappedClass,
    parent: &ParentFrame,
) -> Result<usize, MatrixSubductionError> {
    let rotation = mapped.operation.0;
    let matches: Vec<usize> = certified
        .operations
        .iter()
        .enumerate()
        .filter(|(_, (candidate, translation))| {
            *candidate == rotation
                && differs_by_lattice(&mapped.operation.1, translation, &parent.lattice)
        })
        .map(|(index, _)| index)
        .collect();
    let Some(&index) = matches.first() else {
        return Err(MatrixSubductionError::SourceOperationUnmatched {
            sg,
            label,
            rotation,
            time_reversal: mapped.time_reversal,
        });
    };
    if matches
        .iter()
        .any(|&other| certified.matrices[other] != certified.matrices[index])
    {
        return Err(MatrixSubductionError::SourceOperationAmbiguous {
            sg,
            label,
            rotation,
        });
    }
    Ok(index)
}

/// The ordered composition of two certified source operators.
fn compose_source_operators(
    left: &SourceMatrixOperator,
    right: &SourceMatrixOperator,
) -> SourceMatrixOperator {
    let factor = if left.antiunitary {
        // conj(A(h)) on real entries is A(h) itself; the branch is structural.
        right.entries
    } else {
        right.entries
    };
    let mut entries = [0.0f64; 4];
    for row in 0..2 {
        for column in 0..2 {
            entries[row * 2 + column] = (0..2)
                .map(|k| left.entries[row * 2 + k] * factor[k * 2 + column])
                .sum();
        }
    }
    SourceMatrixOperator {
        entries,
        antiunitary: left.antiunitary ^ right.antiunitary,
    }
}

/// `A(g)·[conj(A(h)) if t(g)] = A(gh)` on every ordered pair of classes.
fn check_source_law(
    classes: &[MagneticOperation],
    actions: &[SourceMatrixOperator],
    lattice: &[[Rat; 3]],
) -> Result<(), MatrixSubductionError> {
    for left in 0..classes.len() {
        for right in 0..classes.len() {
            let product = compose(&classes[left], &classes[right]).map_err(|error| {
                MatrixSubductionError::Arithmetic {
                    reason: error.to_string(),
                }
            })?;
            let expected_class = find_model_class(&product, classes, lattice)
                .ok_or(MatrixSubductionError::SourceLawUnmatched { left, right })?;
            let found = compose_source_operators(&actions[left], &actions[right]);
            if found.antiunitary != actions[expected_class].antiunitary
                || !matrices_close(&found.entries, &actions[expected_class].entries)
            {
                return Err(MatrixSubductionError::SourceLawViolated {
                    left,
                    right,
                    expected_class,
                    expected: actions[expected_class],
                    found,
                });
            }
        }
    }
    Ok(())
}

/// The class an operation denotes, or `None` when it is not a class.
fn find_model_class(
    operation: &MagneticOperation,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
) -> Option<usize> {
    let canonical = canonical_translation(operation.translation, lattice);
    classes.iter().position(|candidate| {
        candidate.rotation == operation.rotation
            && candidate.time_reversal == operation.time_reversal
            && canonical_translation(candidate.translation, lattice) == canonical
    })
}

/// The class of `left · right`, or `None` when the product or its class is
/// missing.
fn product_class(
    left: &MagneticOperation,
    right: &MagneticOperation,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
) -> Option<usize> {
    compose(left, right)
        .ok()
        .and_then(|product| find_model_class(&product, classes, lattice))
}

// ---------------------------------------------------------------------------
// The subduction
// ---------------------------------------------------------------------------

/// Decompose the model's two-dimensional restriction into the magnetic group's Γ
/// corepresentations, exhibiting matrix intertwiners with a full-space
/// reconstruction.
pub fn subduce_matrix_gamma(
    model: &MatrixGammaModel,
) -> Result<MatrixSubduction<'_>, MatrixSubductionError> {
    let source_irrep = resolve_matrix_source(SourceIdentity {
        parent_sg: model.parent_sg,
        label: model.source_label,
        convention: model.source_convention,
    })?;
    let uni = model.uni;
    let partial =
        magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni).map_err(|error| {
            MatrixSubductionError::SummaryFailed {
                uni,
                reason: format!("{error:?}"),
            }
        })?;
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.coords == GAMMA)
        .ok_or(MatrixSubductionError::NoGammaKPoint { uni })?;
    for unresolved in &partial.unresolved_coreps {
        if unresolved.k_label == gamma.label
            && !unresolved.spinor
            && unresolved
                .minimum_dimension
                .is_none_or(|dimension| dimension <= SOURCE_DIMENSION)
        {
            return Err(MatrixSubductionError::UnresolvedTargetData {
                uni,
                source: unresolved.source_irrep.clone(),
                minimum_dimension: unresolved.minimum_dimension,
                reason: unresolved.reason.clone(),
            });
        }
    }

    // The R27 proof that the summary's Γ columns really are the raw first-Hall
    // operations, and that position, index and column associate exactly.
    let class_columns =
        magnetic_subduction::proven_gamma_class_columns(uni, gamma, &model.classes, &model.lattice)
            .map_err(|error| MatrixSubductionError::TargetProvenance {
                uni,
                reason: error.to_string(),
            })?;

    let context = TargetContext::build(uni, gamma, model, &class_columns)?;
    let source_corep = RealCorep::from_source(&model.actions);
    let mut examined: Vec<String> = Vec::new();
    let mut pieces: Vec<MatrixSubductionPiece> = Vec::new();
    for corep in &gamma.coreps {
        if !is_matrix_candidate(corep) {
            continue;
        }
        examined.push(corep.label.clone());
        let tables = build_target_tables(corep, &context, uni)?;
        for table in tables {
            if let Some(piece) = match_piece(corep, &table, &source_corep, &context, uni)? {
                if pieces.iter().any(|kept| kept.operators == piece.operators) {
                    continue;
                }
                pieces.push(piece);
                break;
            }
        }
    }

    if pieces.is_empty() {
        return Err(MatrixSubductionError::NoMatchingTarget {
            parent_sg: model.parent_sg,
            source: source_irrep.ml,
            uni,
            examined,
        });
    }
    let covered: u32 = pieces.iter().map(|piece| piece.covered_dimension()).sum();
    if covered != SOURCE_DIMENSION as u32 {
        return Err(MatrixSubductionError::IncompleteDecomposition {
            source_dimension: SOURCE_DIMENSION as u32,
            covered_dimension: covered,
        });
    }
    // The exhibited copies must cover the whole two-dimensional source space,
    // not merely add up dimensionally: `V†V = I` and `Σ U U† = I₂`.
    if !verify_full_isometry(&pieces) {
        return Err(MatrixSubductionError::SourceSpaceNotIsometric {
            covered_dimension: covered,
        });
    }
    for piece in &pieces {
        verify_reconstruction(model, piece)?;
    }

    Ok(MatrixSubduction {
        model,
        pieces,
        source_dimension: SOURCE_DIMENSION as u32,
        covered_dimension: covered,
    })
}

/// Whether a Γ corep is inside this module's target family.
fn is_matrix_candidate(corep: &MagneticCorepSummary) -> bool {
    corep.corep_type == CorepType::A
        && corep.source_irreps.len() == 1
        && !corep.source_irreps[0].spinor
        && corep.source_irreps[0].dim as usize == corep.dim
        && (corep.dim == 1 || corep.dim == 2)
}

/// Everything the target construction needs from the magnetic group and the
/// summary's operation columns.
struct TargetContext {
    uni: usize,
    h_sg: u8,
    h_seitz: Vec<wigner::SeitzOp>,
    setting: Option<wigner::SettingTransform>,
    classes: Vec<MagneticOperation>,
    lattice: Vec<[Rat; 3]>,
    class_columns: Vec<Vec<usize>>,
}

impl TargetContext {
    fn build(
        uni: usize,
        gamma: &MagneticKPointSummary,
        model: &MatrixGammaModel,
        class_columns: &[Vec<usize>],
    ) -> Result<Self, MatrixSubductionError> {
        let info = corep::cached_unitary_subgroup_with_hall(uni).ok_or(
            MatrixSubductionError::TargetSourceUnresolved {
                uni,
                sg: model.parent_sg,
                label: "<unitary subgroup>".to_string(),
                reason: "the unitary subgroup of the UNI could not be identified".to_string(),
            },
        )?;
        let h_seitz = wigner::ops_to_seitz(&info.ops_from_hall);
        let _ = gamma;
        Ok(Self {
            uni,
            h_sg: info.sg as u8,
            h_seitz,
            setting: info.msg_to_data.clone(),
            classes: model.classes.clone(),
            lattice: model.lattice.clone(),
            class_columns: class_columns.to_vec(),
        })
    }

    /// The H operation index of a unitary magnetic class, in the H data-Hall
    /// frame the summary itself uses.
    fn unitary_operation_index(&self, class: usize) -> Result<usize, MatrixSubductionError> {
        let operation = &self.classes[class];
        let translation = [
            operation.translation[0].to_f64(),
            operation.translation[1].to_f64(),
            operation.translation[2].to_f64(),
        ];
        let (rotation, translation) = match &self.setting {
            Some(setting) => setting
                .transform_seitz(&operation.rotation, &translation)
                .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
                    uni: self.uni,
                    label: String::new(),
                    reason: format!("class {class} does not transform to the H frame"),
                })?,
            None => (operation.rotation, translation),
        };
        let matched =
            wigner::find_seitz(&rotation, &translation, &self.h_seitz).ok_or_else(|| {
                MatrixSubductionError::TargetCorepInconsistent {
                    uni: self.uni,
                    label: String::new(),
                    reason: format!(
                        "class {class} has no operation in the unitary subgroup's Hall list"
                    ),
                }
            })?;
        Ok(matched.op_index)
    }
}

/// Resolve, certify and (where the group has an antiunitary coset) extend one
/// candidate's target table, returning every legal coset gauge.
fn build_target_tables(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
    uni: usize,
) -> Result<Vec<TargetTable>, MatrixSubductionError> {
    let source = &corep.source_irreps[0];
    let candidates: Vec<&'static IrrepRecord> = query::irreps_of(source.sg)
        .iter()
        .filter(|irrep| irrep.ml == source.ml)
        .collect();
    let irrep = match candidates.as_slice() {
        [irrep] => *irrep,
        _ => {
            return Err(MatrixSubductionError::TargetSourceUnresolved {
                uni,
                sg: source.sg,
                label: source.ml.to_string(),
                reason: format!("{} database rows carry this label", candidates.len()),
            });
        }
    };
    if irrep.spinor || irrep.compound_metadata().is_some() {
        return Err(MatrixSubductionError::TargetSourceUnresolved {
            uni,
            sg: source.sg,
            label: source.ml.to_string(),
            reason: "the target source is spinor or compound".to_string(),
        });
    }
    let k = irrep.k_vector();
    if k.numerators != [0, 0, 0] || k.denominator != 1 {
        return Err(MatrixSubductionError::TargetSourceUnresolved {
            uni,
            sg: source.sg,
            label: source.ml.to_string(),
            reason: "the target source is not a Gamma irrep".to_string(),
        });
    }
    let certified = certify_irrep_matrices(irrep, corep.dim).map_err(|problem| {
        MatrixSubductionError::TargetMatrixUncertified {
            uni,
            sg: source.sg,
            label: irrep.ml,
            problem,
        }
    })?;
    let dimension = certified.dimension();

    let h_to_pir = wigner::build_h_to_irrep_op_map(
        &context.h_seitz,
        irrep.pir_rotations(),
        irrep.pir_translations(),
    )
    .ok_or(MatrixSubductionError::TargetMatrixUncertified {
        uni,
        sg: source.sg,
        label: irrep.ml,
        problem: MatrixCertificationProblem::TypedRowShapeMismatch {
            row: context.h_seitz.len(),
            operators: certified.operations.len(),
            dimension,
        },
    })?;

    // Unitary part: H's own matrices, cross-checked against the summary's
    // character values for the same class.
    let mut unitary: Vec<Option<Vec<f64>>> = vec![None; context.classes.len()];
    for (class, slot) in unitary.iter_mut().enumerate() {
        if context.classes[class].time_reversal {
            continue;
        }
        let h_index = context
            .unitary_operation_index(class)
            .map_err(|error| label_error(error, uni, &corep.label))?;
        let pir_index = *h_to_pir.get(h_index).ok_or_else(|| {
            MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!("H operation {h_index} has no PIR row"),
            }
        })?;
        let matrix = certified.matrices[pir_index].clone();
        let trace = matrix_trace(&matrix, dimension);
        let mut seen = None;
        for &column in &context.class_columns[class] {
            let character = corep
                .characters
                .get(column)
                .copied()
                .flatten()
                .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
                    uni,
                    label: corep.label.clone(),
                    reason: format!("unitary class {class} has no summary character"),
                })?;
            if !close(character.im, 0.0) || !close(character.re, trace) {
                return Err(MatrixSubductionError::TargetCorepInconsistent {
                    uni,
                    label: corep.label.clone(),
                    reason: format!(
                        "class {class}: certified target trace {trace} != summary character {character:?}"
                    ),
                });
            }
            seen = Some(character);
        }
        if seen.is_none() {
            return Err(MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!("unitary class {class} is not covered by the summary columns"),
            });
        }
        *slot = Some(matrix);
    }

    let first_anti = context.classes.iter().position(|class| class.time_reversal);
    let Some(first_anti) = first_anti else {
        let operators = (0..context.classes.len())
            .map(|class| TargetMatrixOperator {
                dimension,
                entries: unitary[class]
                    .clone()
                    .expect("a Type I ordinary group has only unitary classes"),
                antiunitary: false,
            })
            .collect();
        return Ok(vec![TargetTable {
            dimension,
            operators,
        }]);
    };

    // The coset matrix X solves U_h X = X U_sigma(h) for every unitary class,
    // where `sigma` conjugates a unitary class by the canonical antiunitary
    // representative; the antiunitary table is then D(a) = X U_{a0^{-1} a}.
    let identities = find_identity_class(&context.classes).ok_or_else(|| {
        MatrixSubductionError::TargetCorepInconsistent {
            uni,
            label: corep.label.clone(),
            reason: "the class list has no identity class".to_string(),
        }
    })?;
    let inverse = inverse_class(&context.classes, &context.lattice, first_anti, identities)
        .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
            uni,
            label: corep.label.clone(),
            reason: "the antiunitary representative has no inverse class".to_string(),
        })?;
    let mut sigma = vec![0usize; context.classes.len()];
    for (class, slot) in sigma.iter_mut().enumerate() {
        if context.classes[class].time_reversal {
            continue;
        }
        let left =
            compose(&context.classes[inverse], &context.classes[class]).map_err(|error| {
                MatrixSubductionError::Arithmetic {
                    reason: error.to_string(),
                }
            })?;
        let conjugated = compose(&left, &context.classes[first_anti]).map_err(|error| {
            MatrixSubductionError::Arithmetic {
                reason: error.to_string(),
            }
        })?;
        let sigma_class = find_model_class(&conjugated, &context.classes, &context.lattice)
            .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!("the conjugate of unitary class {class} is not a class"),
            })?;
        if context.classes[sigma_class].time_reversal {
            return Err(MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!("the conjugate of unitary class {class} is antiunitary"),
            });
        }
        *slot = sigma_class;
    }
    let mut anti_source = vec![usize::MAX; context.classes.len()];
    let mut seen_unitary = vec![false; context.classes.len()];
    for (class, slot) in anti_source.iter_mut().enumerate() {
        if !context.classes[class].time_reversal {
            continue;
        }
        let product =
            compose(&context.classes[inverse], &context.classes[class]).map_err(|error| {
                MatrixSubductionError::Arithmetic {
                    reason: error.to_string(),
                }
            })?;
        let unitary_class = find_model_class(&product, &context.classes, &context.lattice)
            .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!("the coset of class {class} has no unitary representative"),
            })?;
        if context.classes[unitary_class].time_reversal {
            return Err(MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!("the coset of class {class} is not unitary"),
            });
        }
        if std::mem::replace(&mut seen_unitary[unitary_class], true) {
            return Err(MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: corep.label.clone(),
                reason: format!(
                    "unitary class {unitary_class} is the coset representative of several antiunitary classes"
                ),
            });
        }
        *slot = unitary_class;
    }

    // X's linear system.
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let width = dimension * dimension;
    for class in 0..context.classes.len() {
        if context.classes[class].time_reversal {
            continue;
        }
        let u = unitary[class]
            .as_ref()
            .expect("a unitary class has a matrix");
        let s = unitary[sigma[class]]
            .as_ref()
            .expect("sigma of a unitary class is unitary");
        for row in 0..dimension {
            for column in 0..dimension {
                let mut equation = vec![0.0; width];
                for k in 0..dimension {
                    equation[k * dimension + column] += u[row * dimension + k];
                }
                for k in 0..dimension {
                    equation[row * dimension + k] -= s[k * dimension + column];
                }
                rows.push(equation);
            }
        }
    }
    let kernel = real_kernel(&rows, width);
    if kernel.is_empty() {
        return Err(MatrixSubductionError::TargetCorepInconsistent {
            uni,
            label: corep.label.clone(),
            reason: "no real coset matrix satisfies U_h X = X U_sigma(h)".to_string(),
        });
    }
    if kernel.len() != 1 {
        return Err(MatrixSubductionError::TargetCorepInconsistent {
            uni,
            label: corep.label.clone(),
            reason: format!(
                "the coset matrix is not unique (kernel dimension {})",
                kernel.len()
            ),
        });
    }
    let base = &kernel[0];
    let square = product_class(
        &context.classes[first_anti],
        &context.classes[first_anti],
        &context.classes,
        &context.lattice,
    )
    .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
        uni,
        label: corep.label.clone(),
        reason: "the antiunitary representative does not square to a class".to_string(),
    })?;
    let square_matrix = unitary[square]
        .as_ref()
        .expect("the square of an antiunitary representative is unitary");
    let base_square = matrix_product(base, base, dimension);
    let magnitude = coset_scale(&base_square, square_matrix, dimension).ok_or_else(|| {
        let reason = "X^2 is not a positive finite multiple of D(a0^2), so no real coset matrix squares to the unitary representative"
            .to_string();
        MatrixSubductionError::TargetCorepInconsistent {
            uni,
            label: corep.label.clone(),
            reason,
        }
    })?;

    let mut tables = Vec::new();
    for sign in [1.0f64, -1.0] {
        let x: Vec<f64> = base.iter().map(|value| sign * magnitude * value).collect();
        let mut operators: Vec<TargetMatrixOperator> = Vec::with_capacity(context.classes.len());
        for class in 0..context.classes.len() {
            let entries = if context.classes[class].time_reversal {
                matrix_product(
                    &x,
                    &unitary[anti_source[class]].clone().expect("coset source"),
                    dimension,
                )
            } else {
                unitary[class].clone().expect("unitary class")
            };
            if !is_orthogonal(&entries, dimension) {
                return Err(MatrixSubductionError::TargetCorepInconsistent {
                    uni,
                    label: corep.label.clone(),
                    reason: format!("the antiunitary operator of class {class} is not orthogonal"),
                });
            }
            operators.push(TargetMatrixOperator {
                dimension,
                entries,
                antiunitary: context.classes[class].time_reversal,
            });
        }
        let table = TargetTable {
            dimension,
            operators,
        };
        if check_target_law(&table, context).is_ok() {
            tables.push(table);
        }
    }
    if tables.is_empty() {
        return Err(MatrixSubductionError::TargetCorepInconsistent {
            uni,
            label: corep.label.clone(),
            reason: "no coset gauge satisfies the full target corepresentation law".to_string(),
        });
    }
    Ok(tables)
}

/// Re-label a target error with the corep label it was raised for.
fn label_error(error: MatrixSubductionError, uni: usize, label: &str) -> MatrixSubductionError {
    match error {
        MatrixSubductionError::TargetCorepInconsistent { reason, .. } => {
            MatrixSubductionError::TargetCorepInconsistent {
                uni,
                label: label.to_string(),
                reason,
            }
        }
        other => other,
    }
}

/// The identity class index, or `None`.
fn find_identity_class(classes: &[MagneticOperation]) -> Option<usize> {
    classes.iter().position(|class| {
        !class.time_reversal
            && class.rotation == identity_rotation()
            && class.translation.iter().all(|value| value.is_zero())
    })
}

/// The inverse class of one class, or `None` when it is not in the list.
fn inverse_class(
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
    class: usize,
    identity: usize,
) -> Option<usize> {
    (0..classes.len()).find(|&candidate| {
        product_class(&classes[candidate], &classes[class], classes, lattice) == Some(identity)
    })
}

/// One fully constructed target table for the magnetic group's classes.
#[derive(Debug, Clone, PartialEq)]
struct TargetTable {
    dimension: usize,
    operators: Vec<TargetMatrixOperator>,
}

/// The full ordered-pair corepresentation law of one target table.
fn check_target_law(
    table: &TargetTable,
    context: &TargetContext,
) -> Result<(), MatrixSubductionError> {
    for left in 0..context.classes.len() {
        for right in 0..context.classes.len() {
            let expected = product_class(
                &context.classes[left],
                &context.classes[right],
                &context.classes,
                &context.lattice,
            )
            .ok_or_else(|| MatrixSubductionError::TargetCorepInconsistent {
                uni: context.uni,
                label: String::new(),
                reason: format!("the product of classes {left} and {right} is not a class"),
            })?;
            let left_operator = &table.operators[left];
            let right_operator = &table.operators[right];
            let factor: Vec<f64> = if left_operator.antiunitary {
                // conj(A(h)) is A(h) itself on the real certified entries; the
                // branch is structural.
                right_operator.entries.to_vec()
            } else {
                right_operator.entries.clone()
            };
            let entries = matrix_product(&left_operator.entries, &factor, table.dimension);
            let found = TargetMatrixOperator {
                dimension: table.dimension,
                entries,
                antiunitary: left_operator.antiunitary ^ right_operator.antiunitary,
            };
            if found.antiunitary != table.operators[expected].antiunitary
                || !matrices_close(&found.entries, &table.operators[expected].entries)
            {
                return Err(MatrixSubductionError::TargetCorepInconsistent {
                    uni: context.uni,
                    label: String::new(),
                    reason: format!("target law fails at classes ({left}, {right}) -> {expected}"),
                });
            }
        }
    }
    Ok(())
}

/// The best available `2×d` intertwiner copies for one target table.
/// The real-linear intertwiner kernel of `source → target`, and the normalized
/// independent copies taken from it.
///
/// This is deliberately shared by the Type-A and restricted Type-C paths: the
/// copy count is **not** the kernel dimension.  Copies are individually
/// normalized (`U†U = I`) and pairwise orthogonal, and the caller still has to
/// prove that they cover the source space.
fn matching_copies(
    source: &RealCorep,
    table: &TargetTable,
) -> (Vec<Vec<f64>>, Vec<Vec<Complex64>>) {
    let dimension = table.dimension;
    let target = RealCorep::from_target(&table.operators);
    let (rows, width) = hom_rows(source, &target);
    let kernel = real_kernel(&rows, width);
    let mut chosen: Vec<Vec<Complex64>> = Vec::new();
    for vector in &kernel {
        let matrix = intertwiner_from_vector(vector, dimension);
        let Some(normalized) = normalize_intertwiner(&matrix, dimension) else {
            continue;
        };
        if chosen
            .iter()
            .any(|kept| !copies_orthogonal(kept, &normalized, dimension))
        {
            continue;
        }
        chosen.push(normalized);
        if chosen.len() * dimension >= SOURCE_DIMENSION {
            break;
        }
    }
    (kernel, chosen)
}

fn match_piece(
    corep: &MagneticCorepSummary,
    table: &TargetTable,
    source: &RealCorep,
    context: &TargetContext,
    uni: usize,
) -> Result<Option<MatrixSubductionPiece>, MatrixSubductionError> {
    let dimension = table.dimension;
    let (_, chosen) = matching_copies(source, table);
    if chosen.is_empty() {
        return Ok(None);
    }
    let multiplicity = chosen.len() as u32;
    if (multiplicity as usize) * dimension > SOURCE_DIMENSION {
        return Ok(None);
    }
    let operators = table.operators.clone();
    let piece = MatrixSubductionPiece {
        label: corep.label.clone(),
        source_label: corep.source_irreps[0].ml,
        unitary_sg: context.h_sg,
        dimension: dimension as u32,
        multiplicity,
        copies: chosen,
        operators,
    };
    let _ = uni;
    Ok(Some(piece))
}

/// Whether the `d×d` cross block `left† right` vanishes.
fn copies_orthogonal(left: &[Complex64], right: &[Complex64], dimension: usize) -> bool {
    for row in 0..dimension {
        for column in 0..dimension {
            let block: Complex64 = (0..SOURCE_DIMENSION)
                .map(|k| left[k * dimension + row].conj() * right[k * dimension + column])
                .sum();
            if !small(block.norm(), 1e-8) {
                return false;
            }
        }
    }
    true
}

/// The explicit metric statement of a successful decomposition: the exhibited
/// copies form one orthonormal basis of the source space,
/// `V†V = I` **and** `Σ U U† = V V† = I₂`.
fn verify_full_isometry(pieces: &[MatrixSubductionPiece]) -> bool {
    let total: usize = pieces
        .iter()
        .map(|piece| piece.copies.len() * piece.dimension as usize)
        .sum();
    if total != SOURCE_DIMENSION {
        return false;
    }
    let mut columns: Vec<Vec<Complex64>> = Vec::with_capacity(total);
    for piece in pieces {
        let dimension = piece.dimension as usize;
        for copy in &piece.copies {
            for column in 0..dimension {
                columns.push(vec![copy[column], copy[dimension + column]]);
            }
        }
    }
    // V†V = I_total
    for left in 0..total {
        for right in 0..total {
            let entry: Complex64 = (0..SOURCE_DIMENSION)
                .map(|row| columns[left][row].conj() * columns[right][row])
                .sum();
            let expected = if left == right {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            if !small((entry - expected).norm(), 1e-8) {
                return false;
            }
        }
    }
    // Σ U U† = V V† = I₂
    for row in 0..SOURCE_DIMENSION {
        for column in 0..SOURCE_DIMENSION {
            let entry: Complex64 = (0..total)
                .map(|index| columns[index][row] * columns[index][column].conj())
                .sum();
            let expected = if row == column {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            if !small((entry - expected).norm(), 1e-8) {
                return false;
            }
        }
    }
    true
}

/// Normalize one intertwiner so that `U†U = I`, or `None` when the columns are
/// not orthogonal-by-Schur.
fn normalize_intertwiner(matrix: &[Complex64], dimension: usize) -> Option<Vec<Complex64>> {
    let mut gram = vec![Complex64::new(0.0, 0.0); dimension * dimension];
    for row in 0..dimension {
        for column in 0..dimension {
            gram[row * dimension + column] = (0..SOURCE_DIMENSION)
                .map(|k| matrix[k * dimension + row].conj() * matrix[k * dimension + column])
                .sum();
        }
    }
    let scalar = gram[0];
    // The Gram scalar must be finite, real and strictly positive; a `NaN` or a
    // non-positive value is a refusal, not a normalization.
    // `scalar.re` is finite here, so `<= 0.0` is the total-order complement of
    // "strictly positive" and does not hide a NaN.
    if !scalar.re.is_finite()
        || !scalar.im.is_finite()
        || scalar.re <= 0.0
        || !close(scalar.im, 0.0)
    {
        return None;
    }
    for row in 0..dimension {
        for column in 0..dimension {
            let expected = if row == column {
                scalar
            } else {
                Complex64::new(0.0, 0.0)
            };
            if !small(
                (gram[row * dimension + column] - expected).norm(),
                GRID_TOLERANCE,
            ) {
                return None;
            }
        }
    }
    let scale = 1.0 / scalar.re.sqrt();
    let normalized: Vec<Complex64> = matrix.iter().map(|value| value * scale).collect();
    if normalized
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return None;
    }
    Some(normalized)
}

/// `S(g)∘φ = φ∘T(g)` for every class and every copy, on complex vectors.
fn verify_reconstruction(
    model: &MatrixGammaModel,
    piece: &MatrixSubductionPiece,
) -> Result<(), MatrixSubductionError> {
    verify_reconstruction_data(
        model,
        &piece.label,
        piece.dimension as usize,
        &piece.copies,
        &piece.operators,
    )
}

/// The shared complex-vector reconstruction check: it only needs the class
/// actions, the target label for diagnostics, the copies and their operators.
fn verify_reconstruction_data(
    model: &MatrixGammaModel,
    label: &str,
    dimension: usize,
    copies: &[Vec<Complex64>],
    operators: &[TargetMatrixOperator],
) -> Result<(), MatrixSubductionError> {
    let mut probes: Vec<Vec<Complex64>> = (0..dimension)
        .map(|index| {
            (0..dimension)
                .map(|column| {
                    if column == index {
                        Complex64::new(1.0, 0.0)
                    } else {
                        Complex64::new(0.0, 0.0)
                    }
                })
                .collect()
        })
        .collect();
    probes.push(
        (0..dimension)
            .map(|column| Complex64::new(1.0 + column as f64, 1.0 - 2.0 * column as f64))
            .collect(),
    );

    for copy in copies {
        for probe in &probes {
            let image: Vec<Complex64> = (0..SOURCE_DIMENSION)
                .map(|row| {
                    (0..dimension)
                        .map(|column| copy[row * dimension + column] * probe[column])
                        .sum::<Complex64>()
                })
                .collect();
            for (class, operator) in operators.iter().enumerate().take(model.classes.len()) {
                let source = model.actions[class].apply([image[0], image[1]]);
                let target = operator.apply(probe);
                let expected: Vec<Complex64> = (0..SOURCE_DIMENSION)
                    .map(|row| {
                        (0..dimension)
                            .map(|column| copy[row * dimension + column] * target[column])
                            .sum::<Complex64>()
                    })
                    .collect();
                if (0..SOURCE_DIMENSION)
                    .any(|row| !small((source[row] - expected[row]).norm(), 1e-8))
                {
                    return Err(MatrixSubductionError::TargetCorepInconsistent {
                        uni: model.uni,
                        label: label.to_string(),
                        reason: format!(
                            "the complex reconstruction fails at class {class} for a copy"
                        ),
                    });
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// R32: restricted type-C paired-target subduction
// ---------------------------------------------------------------------------

/// The one-line scope fence of [`subduce_matrix_gamma_type_c`].
pub const TYPE_C_SCOPE: &str = "Two-dimensional real non-spinor non-compound Gamma parent irrep \
(certified operation-aware matrix table) restricted to an embedded magnetic group; target: one \
Wigner type-C two-dimensional Gamma corepresentation whose unitary subgroup restriction is \
exactly two inequivalent non-spinor real ±1 one-dimensional H irreps; the antiunitary coset \
exchanges the pair.  Spinor, compound, non-Gamma, complex H characters, type-A/B targets, \
other H-pair shapes and setting-search cases are outside this entry point.";

/// Why one generated `H` character value is not usable as the real `±1` datum
/// this restricted path requires.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TypeCCharacterProblem {
    /// NaN or infinity: unresolved data, not a domain exclusion.
    NonFinite,
    /// Finite but with a non-zero imaginary part.
    NonReal,
    /// Finite and real, but not within tolerance of `+1` or `-1`.
    OffGrid,
}

impl TypeCCharacterProblem {
    fn name(self) -> &'static str {
        match self {
            Self::NonFinite => "NonFinite",
            Self::NonReal => "NonReal",
            Self::OffGrid => "OffGrid",
        }
    }
}

/// Classify a stored real matrix entry on the finite `±1` grid.
fn type_c_real_sign(value: f64) -> Result<i8, TypeCCharacterProblem> {
    if !value.is_finite() {
        return Err(TypeCCharacterProblem::NonFinite);
    }
    if (value - 1.0).abs() <= GRID_TOLERANCE {
        Ok(1)
    } else if (value + 1.0).abs() <= GRID_TOLERANCE {
        Ok(-1)
    } else {
        Err(TypeCCharacterProblem::OffGrid)
    }
}

/// Classify one generated complex character value on the finite real `±1` grid.
///
/// Non-finite values are unresolved; a finite non-zero imaginary part and a
/// finite real value off the grid are domain exclusions.  This is a tolerance
/// classification of generated `f64` data, not a proof about its exact binary
/// value.
fn type_c_complex_sign(value: Complex64) -> Result<i8, TypeCCharacterProblem> {
    if !value.re.is_finite() || !value.im.is_finite() {
        return Err(TypeCCharacterProblem::NonFinite);
    }
    if value.im.abs() > GRID_TOLERANCE {
        return Err(TypeCCharacterProblem::NonReal);
    }
    type_c_real_sign(value.re)
}

/// Why the restricted paired type-C subduction was refused.
///
/// This is a separate enum from [`MatrixSubductionError`] so the existing
/// public error type and its exhaustive matches are untouched.  Exclusions
/// (spinor, compound, non-real finite, off-grid, dimension, non-Γ) are reported
/// in [`Self::NoMatchingTarget`]'s `rejected` list; unavailable or non-finite
/// data is an unresolved refusal, never a silent zero.
#[derive(Debug)]
pub enum MatrixTypeCError {
    /// The magnetic corep summary of the model's UNI could not be computed.
    SummaryFailed { uni: usize, reason: String },
    /// The summary has no Γ k-point.
    NoGammaKPoint { uni: usize },
    /// A possibly two-dimensional non-spinor Γ corep could not be constructed,
    /// so the target list is not known to be complete.
    UnresolvedTargetData {
        uni: usize,
        source: String,
        minimum_dimension: Option<usize>,
        reason: String,
    },
    /// The R27 raw-operation provenance of the summary's Γ columns failed.
    TargetProvenance { uni: usize, reason: String },
    /// One H source record is missing, ambiguous or unusable as data.
    HSource {
        sg: u8,
        label: String,
        reason: String,
    },
    /// A finite generated character is outside the real `±1` domain or is
    /// non-finite (see `problem`).
    Character {
        sg: u8,
        label: &'static str,
        problem: &'static str,
    },
    /// The two H sources are operation-by-operation equivalent, so the type-C
    /// exchange is not meaningful: unresolved, not "no candidate".
    EquivalentHSources {
        uni: usize,
        label: String,
        first: &'static str,
        second: &'static str,
    },
    /// The constructed target or its binding to `M` is inconsistent.
    Inconsistent {
        uni: usize,
        label: String,
        stage: &'static str,
        reason: String,
    },
    /// No examined type-C candidate is equivalent to the restriction.
    NoMatchingTarget {
        parent_sg: u8,
        source: &'static str,
        uni: usize,
        examined: Vec<String>,
        rejected: Vec<String>,
    },
    /// More than one type-C candidate matches the restriction.
    AmbiguousTarget { uni: usize, labels: Vec<String> },
    /// The matched piece does not cover the whole two-dimensional source.
    IncompleteDecomposition {
        source_dimension: u32,
        covered_dimension: u32,
    },
    /// The exhibited copy is not one orthonormal full-rank basis of the source.
    SourceSpaceNotIsometric { covered_dimension: u32 },
    /// Exact arithmetic failed inside the checker.
    Arithmetic { reason: String },
}

impl std::fmt::Display for MatrixTypeCError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SummaryFailed { uni, reason } => {
                write!(f, "UNI {uni}: the magnetic corep summary failed: {reason}")
            }
            Self::NoGammaKPoint { uni } => write!(f, "UNI {uni}: the summary has no Γ k-point"),
            Self::UnresolvedTargetData {
                uni,
                source,
                minimum_dimension,
                reason,
            } => write!(
                f,
                "UNI {uni}: the Γ corep of source {source} (minimum dimension {minimum_dimension:?}) \
                 is unresolved: {reason}"
            ),
            Self::TargetProvenance { uni, reason } => write!(
                f,
                "UNI {uni}: the R27 raw-operation provenance of the Γ columns failed: {reason}"
            ),
            Self::HSource { sg, label, reason } => write!(
                f,
                "SG {sg}: the paired type-C H source {label:?} is unusable: {reason}"
            ),
            Self::Character { sg, label, problem } => write!(
                f,
                "SG {sg}: the H character of {label:?} is {problem}; the paired type-C domain \
                 requires finite real ±1 values"
            ),
            Self::EquivalentHSources {
                uni,
                label,
                first,
                second,
            } => write!(
                f,
                "UNI {uni}: the type-C corep {label} pairs {first} and {second}, whose character \
                 rows are operation-by-operation equal; the exchange cannot be decided"
            ),
            Self::Inconsistent {
                uni,
                label,
                stage,
                reason,
            } => write!(
                f,
                "UNI {uni}: the type-C target {label} is inconsistent at {stage}: {reason}"
            ),
            Self::NoMatchingTarget {
                parent_sg,
                source,
                uni,
                examined,
                rejected,
            } => write!(
                f,
                "SG {parent_sg} source {source} / UNI {uni}: none of the examined type-C targets \
                 ({}) is realised; rejected: {}",
                examined.join(", "),
                rejected.join("; ")
            ),
            Self::AmbiguousTarget { uni, labels } => write!(
                f,
                "UNI {uni}: several type-C targets match the restriction: {}",
                labels.join(", ")
            ),
            Self::IncompleteDecomposition {
                source_dimension,
                covered_dimension,
            } => write!(
                f,
                "the type-C piece covers {covered_dimension} of {source_dimension} source dimensions"
            ),
            Self::SourceSpaceNotIsometric { covered_dimension } => write!(
                f,
                "the {covered_dimension} covered dimensions are not one orthonormal full-rank copy \
                 of the two-dimensional source space"
            ),
            Self::Arithmetic { reason } => write!(f, "exact arithmetic failed: {reason}"),
        }
    }
}

impl std::error::Error for MatrixTypeCError {}

/// One restricted type-C target of the magnetic group, with both independent H
/// source identities and the complex intertwiner that realises it.
#[derive(Debug, Clone)]
pub struct MatrixTypeCPiece {
    label: String,
    unitary_sg: u8,
    h_sources: [&'static str; 2],
    dimension: u32,
    multiplicity: u32,
    copies: Vec<Vec<Complex64>>,
    operators: Vec<TargetMatrixOperator>,
}

impl MatrixTypeCPiece {
    /// The target corep display label as the magnetic summary names it.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Space group number of the target's unitary subgroup `H`.
    pub fn unitary_sg(&self) -> u8 {
        self.unitary_sg
    }

    /// The two canonical CDML labels of the H irreps, in the summary's source
    /// order.  Both are real data identities; neither is a magnetic BC label.
    pub fn h_sources(&self) -> [&'static str; 2] {
        self.h_sources
    }

    /// Dimension of this target corepresentation (two).
    pub fn dimension(&self) -> u32 {
        self.dimension
    }

    /// Multiplicity of this target in the restriction (one independent complex
    /// copy for the restricted two-dimensional real target).
    pub fn multiplicity(&self) -> u32 {
        self.multiplicity
    }

    /// `multiplicity · dimension`, the source dimensions this piece covers.
    pub fn covered_dimension(&self) -> u32 {
        self.multiplicity * self.dimension
    }

    /// The independently normalized `2×2` complex intertwiners, `multiplicity`
    /// of them, row-major with four entries each.
    pub fn copies(&self) -> &[Vec<Complex64>] {
        &self.copies
    }

    /// The reconstructed target operator of every class of `M`, in the model's
    /// class order.
    pub fn operators(&self) -> &[TargetMatrixOperator] {
        &self.operators
    }

    /// The reconstructed target operator of one class.
    pub fn operator_at(&self, index: usize) -> Option<&TargetMatrixOperator> {
        self.operators.get(index)
    }
}

/// The result of one restricted type-C Γ subduction.
#[derive(Debug)]
pub struct MatrixTypeCSubduction<'a> {
    model: &'a MatrixGammaModel,
    piece: MatrixTypeCPiece,
    source_dimension: u32,
    covered_dimension: u32,
}

impl MatrixTypeCSubduction<'_> {
    /// The verified input model this report is about.
    pub fn model(&self) -> &MatrixGammaModel {
        self.model
    }

    /// The single paired type-C piece.
    pub fn piece(&self) -> &MatrixTypeCPiece {
        &self.piece
    }

    /// Dimension of the subduced source corepresentation (two).
    pub fn source_dimension(&self) -> u32 {
        self.source_dimension
    }

    /// `multiplicity · dimension`, checked equal to `source_dimension`.
    pub fn covered_dimension(&self) -> u32 {
        self.covered_dimension
    }

    /// One line stating what this entry point covers.
    pub fn scope(&self) -> &'static str {
        TYPE_C_SCOPE
    }
}

/// One H source of a screened type-C pair, with its certified one-dimensional
/// matrix table and the H-operation → PIR-row map.
#[derive(Clone, Debug)]
struct TypeCPairSource {
    sg: u8,
    label: &'static str,
    h_to_pir: Vec<usize>,
    certified: CertifiedTable,
    typed: CharacterRow,
}

impl TypeCPairSource {
    /// The real `±1` character of one H operation row, from the certified
    /// matrix table and the typed row.  A disagreement is a data inconsistency.
    fn sign_at_pir(
        &self,
        pir: usize,
        corep_label: &str,
        uni: usize,
    ) -> Result<Result<i8, TypeCCharacterProblem>, MatrixTypeCError> {
        let certified = type_c_real_sign(self.certified.matrices[pir][0]);
        let typed = type_c_complex_sign(self.typed.values()[pir]);
        match (certified, typed) {
            (Ok(certified), Ok(typed)) if certified == typed => Ok(Ok(certified)),
            (Ok(_), Ok(_)) => Err(MatrixTypeCError::HSource {
                sg: self.sg,
                label: self.label.to_string(),
                reason: format!(
                    "the certified and typed character signs differ at PIR row {pir} of {corep_label} (UNI {uni})"
                ),
            }),
            (Err(problem), _) | (_, Err(problem)) => Ok(Err(problem)),
        }
    }
}

/// The screened pair: two independent H sources and their class characters.
#[derive(Clone, Debug)]
struct TypeCPair {
    first: TypeCPairSource,
    second: TypeCPairSource,
    /// `characters[source][class]`, `0` on antiunitary classes.
    characters: [Vec<i8>; 2],
}

/// Screen one summary type-C candidate against the restricted real `±1` pair
/// domain.
///
/// `Ok(Ok(pair))` is a candidate, `Ok(Err(reason))` is a named domain exclusion
/// and `Err(_)` is an unresolved/unusable-data refusal.
fn screen_type_c_pair(
    corep: &MagneticCorepSummary,
    context: &TargetContext,
    uni: usize,
) -> Result<Result<TypeCPair, String>, MatrixTypeCError> {
    if corep.corep_type != CorepType::C || corep.dim != 2 {
        return Ok(Err(format!(
            "not-type-C-dim-2({}:{})",
            corep.dim, corep.corep_type as u8
        )));
    }
    let sources = &corep.source_irreps;
    if sources.len() != 2 {
        return Ok(Err(format!("source-count({})", sources.len())));
    }
    let mut facts: Vec<TypeCPairSource> = Vec::with_capacity(2);
    for (position, source) in sources.iter().enumerate() {
        if source.spinor {
            return Ok(Err(format!("spinor({position})")));
        }
        if source.dim != 1 {
            return Ok(Err(format!("dimension({position}:{})", source.dim)));
        }
        if source.sg != context.h_sg {
            return Err(MatrixTypeCError::HSource {
                sg: source.sg,
                label: source.ml.to_string(),
                reason: format!(
                    "the summary H source belongs to SG {} but the magnetic group's unitary subgroup is SG {}",
                    source.sg, context.h_sg
                ),
            });
        }
        let matches = query::find_irreps(source.sg, source.ml, LabelConvention::Cdml);
        let record = match matches.as_slice() {
            [record] => *record,
            [] => {
                return Err(MatrixTypeCError::HSource {
                    sg: source.sg,
                    label: source.ml.to_string(),
                    reason: "no database row carries this CDML label".to_string(),
                });
            }
            _ => {
                return Err(MatrixTypeCError::HSource {
                    sg: source.sg,
                    label: source.ml.to_string(),
                    reason: format!("{} database rows carry this CDML label", matches.len()),
                });
            }
        };
        if record.spinor != source.spinor {
            return Err(MatrixTypeCError::HSource {
                sg: source.sg,
                label: source.ml.to_string(),
                reason: "the summary and database spinor flags disagree".to_string(),
            });
        }
        if record.compound_metadata().is_some() {
            return Ok(Err(format!("compound({position})")));
        }
        if record.dim as usize != 1 {
            return Err(MatrixTypeCError::HSource {
                sg: source.sg,
                label: source.ml.to_string(),
                reason: format!(
                    "the summary selected-arm dimension is 1 but the record declares dimension {}",
                    record.dim
                ),
            });
        }
        let k = record.k_vector();
        if k.numerators != [0, 0, 0] || k.denominator != 1 {
            return Ok(Err(format!("non-Gamma({position})")));
        }
        let certified =
            certify_irrep_matrices(record, 1).map_err(|problem| MatrixTypeCError::HSource {
                sg: source.sg,
                label: source.ml.to_string(),
                reason: format!("the one-dimensional matrix table is not certified: {problem}"),
            })?;
        let typed = record
            .ordinary_scalar_selected_arm_block_trace()
            .map_err(|error| MatrixTypeCError::HSource {
                sg: source.sg,
                label: source.ml.to_string(),
                reason: format!("the typed selected-arm Γ row is unavailable: {error:?}"),
            })?;
        if typed.dimension() != 1 || typed.len() != certified.operations.len() {
            return Err(MatrixTypeCError::HSource {
                sg: source.sg,
                label: source.ml.to_string(),
                reason: format!(
                    "the typed row has dimension {} and {} entries, but the certified table has {} operations",
                    typed.dimension(),
                    typed.len(),
                    certified.operations.len()
                ),
            });
        }
        for index in 0..certified.matrices.len() {
            let certified_sign = match type_c_real_sign(certified.matrices[index][0]) {
                Ok(sign) => sign,
                Err(TypeCCharacterProblem::NonFinite) => {
                    return Err(MatrixTypeCError::Character {
                        sg: source.sg,
                        label: source.ml,
                        problem: "NonFinite",
                    });
                }
                Err(problem) => return Ok(Err(format!("{}({position})", problem.name()))),
            };
            match type_c_complex_sign(typed.values()[index]) {
                Ok(typed_sign) if typed_sign == certified_sign => {}
                Ok(_) => {
                    return Err(MatrixTypeCError::HSource {
                        sg: source.sg,
                        label: source.ml.to_string(),
                        reason: format!(
                            "the typed and certified character signs differ at H row {index}"
                        ),
                    });
                }
                Err(TypeCCharacterProblem::NonFinite) => {
                    return Err(MatrixTypeCError::Character {
                        sg: source.sg,
                        label: source.ml,
                        problem: "NonFinite",
                    });
                }
                Err(problem) => return Ok(Err(format!("{}({position})", problem.name()))),
            }
        }
        let h_to_pir = wigner::build_h_to_irrep_op_map(
            &context.h_seitz,
            record.pir_rotations(),
            record.pir_translations(),
        )
        .ok_or_else(|| MatrixTypeCError::HSource {
            sg: source.sg,
            label: source.ml.to_string(),
            reason: "the H operations do not map one-to-one to the irrep's PIR rows".to_string(),
        })?;
        facts.push(TypeCPairSource {
            sg: source.sg,
            label: source.ml,
            h_to_pir,
            certified,
            typed,
        });
    }

    let mut characters: [Vec<i8>; 2] = [
        vec![0; context.classes.len()],
        vec![0; context.classes.len()],
    ];
    for (class, operation) in context.classes.iter().enumerate() {
        if operation.time_reversal {
            continue;
        }
        let h_index = context.unitary_operation_index(class).map_err(|error| {
            MatrixTypeCError::Inconsistent {
                uni,
                label: corep.label.clone(),
                stage: "H-class-binding",
                reason: error.to_string(),
            }
        })?;
        for (position, source) in facts.iter().enumerate() {
            let pir =
                *source
                    .h_to_pir
                    .get(h_index)
                    .ok_or_else(|| MatrixTypeCError::Inconsistent {
                        uni,
                        label: corep.label.clone(),
                        stage: "H-class-binding",
                        reason: format!("H operation {h_index} has no PIR row"),
                    })?;
            let sign = match source.sign_at_pir(pir, &corep.label, uni)? {
                Ok(sign) => sign,
                Err(TypeCCharacterProblem::NonFinite) => {
                    return Err(MatrixTypeCError::Character {
                        sg: source.sg,
                        label: source.label,
                        problem: "NonFinite",
                    });
                }
                Err(problem) => return Ok(Err(format!("{}({position})", problem.name()))),
            };
            characters[position][class] = sign;
        }
    }
    let separating = (0..context.classes.len()).any(|class| {
        !context.classes[class].time_reversal && characters[0][class] != characters[1][class]
    });
    if !separating {
        return Err(MatrixTypeCError::EquivalentHSources {
            uni,
            label: corep.label.clone(),
            first: facts[0].label,
            second: facts[1].label,
        });
    }

    // The summary's unitary TypeC character column is independently compared
    // with χ1+χ2 on the same class.  Both routes must agree; neither replaces
    // the other.
    for (class, operation) in context.classes.iter().enumerate() {
        if operation.time_reversal {
            continue;
        }
        let columns = &context.class_columns[class];
        if columns.is_empty() {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: corep.label.clone(),
                stage: "summary-character",
                reason: format!("unitary class {class} has no summary column"),
            });
        }
        let expected = f64::from(characters[0][class]) + f64::from(characters[1][class]);
        for &column in columns {
            let value = corep
                .characters
                .get(column)
                .copied()
                .flatten()
                .ok_or_else(|| MatrixTypeCError::Inconsistent {
                    uni,
                    label: corep.label.clone(),
                    stage: "summary-character",
                    reason: format!("unitary class {class} has no summary character"),
                })?;
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(MatrixTypeCError::Inconsistent {
                    uni,
                    label: corep.label.clone(),
                    stage: "summary-character",
                    reason: format!(
                        "unitary class {class} has the non-finite summary character {value:?}"
                    ),
                });
            }
            if value.im.abs() > GRID_TOLERANCE {
                return Err(MatrixTypeCError::Inconsistent {
                    uni,
                    label: corep.label.clone(),
                    stage: "summary-character",
                    reason: format!("unitary class {class} has the non-real character {value:?}"),
                });
            }
            if (value.re - expected).abs() > GRID_TOLERANCE {
                return Err(MatrixTypeCError::Inconsistent {
                    uni,
                    label: corep.label.clone(),
                    stage: "summary-character",
                    reason: format!(
                        "unitary class {class}: summary character {} != χ1+χ2 {expected}",
                        value.re
                    ),
                });
            }
        }
    }

    let second = facts.pop().expect("two H facts");
    let first = facts.pop().expect("two H facts");
    if first.sg != second.sg {
        return Err(MatrixTypeCError::Inconsistent {
            uni,
            label: corep.label.clone(),
            stage: "H-space-group",
            reason: format!(
                "the paired H sources belong to different space groups {} and {}",
                first.sg, second.sg
            ),
        });
    }
    Ok(Ok(TypeCPair {
        first,
        second,
        characters,
    }))
}

/// Build the canonical type-C target table from the two H characters and the
/// actual antiunitary representative `a0`.
///
/// The unitary part is `diag(χ1, χ2)`; the coset matrix is the real
/// `X = [[0, 1], [s, 0]]` with `s = χ1(a0²) = χ2(a0²)`, so `X² = D(a0²)` (this
/// includes the `-I` case).  The twisted pairing `χ1(a0⁻¹ h a0) = χ2(h)` and the
/// full ordered-pair law are checked before the table is returned.
fn build_type_c_table(
    pair: &TypeCPair,
    context: &TargetContext,
    uni: usize,
    label: &str,
    a0: usize,
) -> Result<TargetTable, MatrixTypeCError> {
    if !context.classes[a0].time_reversal {
        return Err(MatrixTypeCError::Inconsistent {
            uni,
            label: label.to_string(),
            stage: "a0",
            reason: format!("class {a0} is not antiunitary"),
        });
    }
    let identity =
        find_identity_class(&context.classes).ok_or_else(|| MatrixTypeCError::Inconsistent {
            uni,
            label: label.to_string(),
            stage: "identity",
            reason: "the class list has no identity class".to_string(),
        })?;
    let inverse =
        inverse_class(&context.classes, &context.lattice, a0, identity).ok_or_else(|| {
            MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "a0-inverse",
                reason: "the antiunitary representative has no inverse class".to_string(),
            }
        })?;
    let square = product_class(
        &context.classes[a0],
        &context.classes[a0],
        &context.classes,
        &context.lattice,
    )
    .ok_or_else(|| MatrixTypeCError::Inconsistent {
        uni,
        label: label.to_string(),
        stage: "square",
        reason: "the antiunitary representative does not square to a class".to_string(),
    })?;
    if context.classes[square].time_reversal {
        return Err(MatrixTypeCError::Inconsistent {
            uni,
            label: label.to_string(),
            stage: "square",
            reason: "the square of the antiunitary representative is antiunitary".to_string(),
        });
    }
    let s = pair.characters[0][square];
    if pair.characters[1][square] != s {
        return Err(MatrixTypeCError::Inconsistent {
            uni,
            label: label.to_string(),
            stage: "square",
            reason: format!(
                "χ1(a0²)={} and χ2(a0²)={} differ, so no diagonal square is possible",
                pair.characters[0][square], pair.characters[1][square]
            ),
        });
    }
    // X² = s·I = diag(χ1(a0²), χ2(a0²)).
    let x = [0.0, 1.0, f64::from(s), 0.0];
    let square_matrix = matrix_product(&x, &x, 2);
    let expected_square = [f64::from(s), 0.0, 0.0, f64::from(s)];
    if !matrices_close(&square_matrix, &expected_square) {
        return Err(MatrixTypeCError::Inconsistent {
            uni,
            label: label.to_string(),
            stage: "square",
            reason: "the canonical coset matrix does not square to D(a0²)".to_string(),
        });
    }
    // The twisted pairing: σ(h) = a0⁻¹ h a0 is a unitary class and the two
    // characters exchange under it.
    for (class, operation) in context.classes.iter().enumerate() {
        if operation.time_reversal {
            continue;
        }
        let left =
            compose(&context.classes[inverse], &context.classes[class]).map_err(|error| {
                MatrixTypeCError::Arithmetic {
                    reason: error.to_string(),
                }
            })?;
        let conjugated =
            compose(&left, &context.classes[a0]).map_err(|error| MatrixTypeCError::Arithmetic {
                reason: error.to_string(),
            })?;
        let sigma =
            find_model_class(&conjugated, &context.classes, &context.lattice).ok_or_else(|| {
                MatrixTypeCError::Inconsistent {
                    uni,
                    label: label.to_string(),
                    stage: "twist",
                    reason: format!("the a0-conjugate of unitary class {class} is not a class"),
                }
            })?;
        if context.classes[sigma].time_reversal {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "twist",
                reason: format!("the a0-conjugate of unitary class {class} is antiunitary"),
            });
        }
        if pair.characters[0][sigma] != pair.characters[1][class]
            || pair.characters[1][sigma] != pair.characters[0][class]
        {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "twist",
                reason: format!(
                    "class {class}: χ1(a0⁻¹ h a0)={}, χ2(h)={}; χ2(a0⁻¹ h a0)={}, χ1(h)={}",
                    pair.characters[0][sigma],
                    pair.characters[1][class],
                    pair.characters[1][sigma],
                    pair.characters[0][class]
                ),
            });
        }
    }

    // Every antiunitary class is a0 times one unitary class; the coset
    // representative is unique.
    let mut anti_source = vec![usize::MAX; context.classes.len()];
    let mut seen_unitary = vec![false; context.classes.len()];
    for (class, operation) in context.classes.iter().enumerate() {
        if !operation.time_reversal {
            continue;
        }
        let product =
            compose(&context.classes[inverse], &context.classes[class]).map_err(|error| {
                MatrixTypeCError::Arithmetic {
                    reason: error.to_string(),
                }
            })?;
        let unitary_class = find_model_class(&product, &context.classes, &context.lattice)
            .ok_or_else(|| MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "coset",
                reason: format!("the coset of class {class} has no unitary representative"),
            })?;
        if context.classes[unitary_class].time_reversal {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "coset",
                reason: format!("the coset of class {class} is not unitary"),
            });
        }
        if std::mem::replace(&mut seen_unitary[unitary_class], true) {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "coset",
                reason: format!(
                    "unitary class {unitary_class} is the coset representative of several antiunitary classes"
                ),
            });
        }
        anti_source[class] = unitary_class;
    }

    let mut operators: Vec<TargetMatrixOperator> = Vec::with_capacity(context.classes.len());
    for (class, operation) in context.classes.iter().enumerate() {
        let entries = if operation.time_reversal {
            let h = anti_source[class];
            let d = [
                f64::from(pair.characters[0][h]),
                0.0,
                0.0,
                f64::from(pair.characters[1][h]),
            ];
            matrix_product(&x, &d, 2)
        } else {
            vec![
                f64::from(pair.characters[0][class]),
                0.0,
                0.0,
                f64::from(pair.characters[1][class]),
            ]
        };
        if !is_orthogonal(&entries, 2) {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: label.to_string(),
                stage: "operator",
                reason: format!("the target operator of class {class} is not orthogonal"),
            });
        }
        operators.push(TargetMatrixOperator {
            dimension: 2,
            entries,
            antiunitary: operation.time_reversal,
        });
    }
    let table = TargetTable {
        dimension: 2,
        operators,
    };
    check_target_law(&table, context).map_err(|error| MatrixTypeCError::Inconsistent {
        uni,
        label: label.to_string(),
        stage: "target-law",
        reason: error.to_string(),
    })?;
    Ok(table)
}

/// The explicit metric statement of one restricted type-C piece: `U†U = I₂`
/// and `|det U| = 1`, so the single copy is a full-rank isometry of the whole
/// two-dimensional source space.  This is the complex-dimension accounting; the
/// real-linear kernel length is never used as a multiplicity.
fn verify_type_c_isometry(piece: &MatrixTypeCPiece) -> bool {
    if piece.dimension != 2 || piece.multiplicity != 1 || piece.copies.len() != 1 {
        return false;
    }
    let copy = &piece.copies[0];
    if copy.len() != 4 {
        return false;
    }
    for row in 0..2 {
        for column in 0..2 {
            let entry: Complex64 = (0..2)
                .map(|k| copy[k * 2 + row].conj() * copy[k * 2 + column])
                .sum();
            let expected = if row == column {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            if !small((entry - expected).norm(), 1e-8) {
                return false;
            }
        }
    }
    let determinant = copy[0] * copy[3] - copy[1] * copy[2];
    if !determinant.re.is_finite()
        || !determinant.im.is_finite()
        || !small((determinant.norm() - 1.0).abs(), 1e-8)
    {
        return false;
    }
    true
}

/// Decompose the model's two-dimensional restriction into a magnetic group's
/// **type-C** Γ corepresentation whose unitary restriction is exactly two
/// inequivalent real `±1` one-dimensional H irreps, exhibiting the complex
/// intertwiner and the reconstructed target operators.
pub fn subduce_matrix_gamma_type_c(
    model: &MatrixGammaModel,
) -> Result<MatrixTypeCSubduction<'_>, MatrixTypeCError> {
    let uni = model.uni;
    let partial =
        magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni).map_err(|error| {
            MatrixTypeCError::SummaryFailed {
                uni,
                reason: format!("{error:?}"),
            }
        })?;
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.coords == GAMMA)
        .ok_or(MatrixTypeCError::NoGammaKPoint { uni })?;
    for unresolved in &partial.unresolved_coreps {
        if unresolved.k_label == gamma.label
            && !unresolved.spinor
            && unresolved
                .minimum_dimension
                .is_none_or(|dimension| dimension <= SOURCE_DIMENSION)
        {
            return Err(MatrixTypeCError::UnresolvedTargetData {
                uni,
                source: unresolved.source_irrep.clone(),
                minimum_dimension: unresolved.minimum_dimension,
                reason: unresolved.reason.clone(),
            });
        }
    }
    let class_columns =
        magnetic_subduction::proven_gamma_class_columns(uni, gamma, &model.classes, &model.lattice)
            .map_err(|error| MatrixTypeCError::TargetProvenance {
                uni,
                reason: error.to_string(),
            })?;
    let context = TargetContext::build(uni, gamma, model, &class_columns).map_err(|error| {
        MatrixTypeCError::Inconsistent {
            uni,
            label: String::new(),
            stage: "context",
            reason: error.to_string(),
        }
    })?;
    let first_anti = match context.classes.iter().position(|class| class.time_reversal) {
        Some(index) => index,
        None => {
            // A type-C corep needs an antiunitary coset that exchanges the two
            // H components.  A group with no antiunitary class is outside the
            // domain, not a data inconsistency: report the examined type-C
            // candidates (if any) and the boundary reason.
            let examined = gamma
                .coreps
                .iter()
                .filter(|corep| corep.corep_type == CorepType::C && corep.dim == 2)
                .map(|corep| corep.label.clone())
                .collect();
            return Err(MatrixTypeCError::NoMatchingTarget {
                parent_sg: model.parent_sg,
                source: model.source_label,
                uni,
                examined,
                rejected: vec!["the magnetic group has no antiunitary coset".to_string()],
            });
        }
    };
    let source_corep = RealCorep::from_source(&model.actions);
    let mut examined: Vec<String> = Vec::new();
    let mut rejected: Vec<String> = Vec::new();
    let mut pieces: Vec<MatrixTypeCPiece> = Vec::new();
    for corep in &gamma.coreps {
        if corep.corep_type != CorepType::C || corep.dim != 2 {
            continue;
        }
        examined.push(corep.label.clone());
        let pair = match screen_type_c_pair(corep, &context, uni)? {
            Ok(pair) => pair,
            Err(reason) => {
                rejected.push(format!("{}: {reason}", corep.label));
                continue;
            }
        };
        let table = build_type_c_table(&pair, &context, uni, &corep.label, first_anti)?;
        // The complex-entry real-linear Hom space of this type-C corep is
        // **two**-dimensional: `I` and `diag(i, -i)` both satisfy the unitary
        // and the antiunitary semi-linear self-intertwiner equations.  That
        // kernel dimension is not a multiplicity: `matching_copies` keeps the
        // independent copies selected by complex-column orthogonality under a
        // fixed source dimension, so the two self-intertwiners do not become
        // two copies of the same target.  The multiplicity below is therefore
        // `copies.len()`, never `kernel.len()`.
        let (_kernel, copies) = matching_copies(&source_corep, &table);
        if copies.is_empty() {
            rejected.push(format!("{}: no-intertwiner", corep.label));
            continue;
        }
        let multiplicity = copies.len() as u32;
        if multiplicity != 1 {
            return Err(MatrixTypeCError::Inconsistent {
                uni,
                label: corep.label.clone(),
                stage: "multiplicity",
                reason: format!(
                    "the type-C target has {multiplicity} independent copies in a two-dimensional source"
                ),
            });
        }
        let piece = MatrixTypeCPiece {
            label: corep.label.clone(),
            unitary_sg: context.h_sg,
            h_sources: [pair.first.label, pair.second.label],
            dimension: 2,
            multiplicity,
            copies,
            operators: table.operators.clone(),
        };
        if pieces
            .iter()
            .any(|kept| kept.operators == piece.operators && kept.h_sources == piece.h_sources)
        {
            continue;
        }
        if let Some(kept) = pieces.iter().find(|kept| kept.operators == piece.operators) {
            return Err(MatrixTypeCError::AmbiguousTarget {
                uni,
                labels: vec![kept.label.clone(), piece.label.clone()],
            });
        }
        pieces.push(piece);
    }
    if pieces.is_empty() {
        return Err(MatrixTypeCError::NoMatchingTarget {
            parent_sg: model.parent_sg,
            source: model.source_label,
            uni,
            examined,
            rejected,
        });
    }
    if pieces.len() > 1 {
        return Err(MatrixTypeCError::AmbiguousTarget {
            uni,
            labels: pieces.iter().map(|piece| piece.label.clone()).collect(),
        });
    }
    let piece = pieces.pop().expect("one matched type-C piece");
    let covered = piece.covered_dimension();
    if covered != SOURCE_DIMENSION as u32 {
        return Err(MatrixTypeCError::IncompleteDecomposition {
            source_dimension: SOURCE_DIMENSION as u32,
            covered_dimension: covered,
        });
    }
    if !verify_type_c_isometry(&piece) {
        return Err(MatrixTypeCError::SourceSpaceNotIsometric {
            covered_dimension: covered,
        });
    }
    verify_reconstruction_data(
        model,
        &piece.label,
        piece.dimension as usize,
        &piece.copies,
        &piece.operators,
    )
    .map_err(|error| MatrixTypeCError::Inconsistent {
        uni,
        label: piece.label.clone(),
        stage: "reconstruction",
        reason: error.to_string(),
    })?;
    Ok(MatrixTypeCSubduction {
        model,
        piece,
        source_dimension: SOURCE_DIMENSION as u32,
        covered_dimension: covered,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The stored operation table, matrix table and translation lattice of one
    /// record, built exactly like the certification reads them.
    struct StoredFixture {
        operations: Vec<(Mat3I, [Rat; 3])>,
        matrices: Vec<Vec<f64>>,
        lattice: Vec<[Rat; 3]>,
        dimension: usize,
    }

    fn stored_table(sg: u8, label: &str) -> StoredFixture {
        let irrep = crate::irrep::query::irreps_of(sg)
            .iter()
            .find(|irrep| irrep.ml == label)
            .expect("the witness record exists");
        let dimension = irrep.dim as usize;
        let entries = dimension * dimension;
        let operators = irrep.characters().len();
        let mut operations = Vec::with_capacity(operators);
        let mut lattice: Vec<[Rat; 3]> = Vec::new();
        for index in 0..operators {
            let cell = &irrep.pir_rotations()[index * 9..index * 9 + 9];
            let rotation = [
                [cell[0], cell[1], cell[2]],
                [cell[3], cell[4], cell[5]],
                [cell[6], cell[7], cell[8]],
            ];
            let mut translation = [Rat::ZERO; 3];
            for (axis, value) in translation.iter_mut().enumerate() {
                *value = exact_fraction(irrep.pir_translations()[index * 3 + axis]).expect("grid");
            }
            if rotation == identity_rotation()
                && translation.iter().any(|value| !value.is_zero())
                && !lattice.contains(&translation)
            {
                lattice.push(translation);
            }
            operations.push((rotation, translation));
        }
        let matrices = (0..operators)
            .map(|index| irrep.matrices()[index * entries..(index + 1) * entries].to_vec())
            .collect();
        StoredFixture {
            operations,
            matrices,
            lattice,
            dimension,
        }
    }

    #[test]
    fn the_certification_accepts_both_witness_tables() {
        for (sg, label) in [(89u8, "GM5"), (156, "GM3")] {
            let irrep = crate::irrep::query::irreps_of(sg)
                .iter()
                .find(|irrep| irrep.ml == label)
                .expect("record");
            let certified = certify_irrep_matrices(irrep, 2).expect("certified table");
            assert_eq!(certified.dimension(), 2);
            assert_eq!(certified.matrices.len(), irrep.characters().len());
        }
    }

    /// The law is the discriminating check: a permutation of two operation rows
    /// keeps every trace (the permutation the two witness space groups admit)
    /// and still has to fail here.
    #[test]
    fn a_same_character_permutation_breaks_the_law_but_not_the_traces() {
        for (sg, label, left, right) in [(89u8, "GM5", 1usize, 3usize), (156, "GM3", 1, 2)] {
            let StoredFixture {
                operations,
                matrices,
                lattice,
                dimension,
            } = stored_table(sg, label);
            let traces: Vec<f64> = matrices
                .iter()
                .map(|matrix| matrix_trace(matrix, dimension))
                .collect();
            assert!(
                check_stored_matrix_law(&operations, &matrices, dimension, &lattice).is_ok(),
                "SG{sg} {label} stores a group law"
            );

            let mut permuted = matrices.clone();
            permuted.swap(left, right);
            let permuted_traces: Vec<f64> = permuted
                .iter()
                .map(|matrix| matrix_trace(matrix, dimension))
                .collect();
            assert_eq!(traces, permuted_traces, "the permutation keeps every trace");
            assert!(
                matches!(
                    check_stored_matrix_law(&operations, &permuted, dimension, &lattice),
                    Err(MatrixCertificationProblem::LawViolated { .. })
                ),
                "SG{sg} {label}: the same-character permutation must break the law"
            );
        }
    }

    #[test]
    fn the_certification_rejects_a_non_orthogonal_table() {
        let StoredFixture {
            operations,
            mut matrices,
            lattice,
            dimension,
        } = stored_table(89, "GM5");
        for value in matrices[1].iter_mut() {
            *value *= 2.0;
        }
        assert!(matches!(
            check_stored_matrix_law(&operations, &matrices, dimension, &lattice),
            Err(MatrixCertificationProblem::LawViolated { .. })
        ));
        assert!(!is_orthogonal(&matrices[1], dimension));
    }

    #[test]
    fn the_real_commutant_is_the_irreducibility_invariant() {
        let StoredFixture {
            matrices,
            dimension,
            ..
        } = stored_table(89, "GM5");
        assert_eq!(real_commutant_dimension(&matrices, dimension), 1);
        let identity_only = vec![vec![1.0, 0.0, 0.0, 1.0]];
        assert_eq!(real_commutant_dimension(&identity_only, 2), 4);
    }

    #[test]
    fn the_real_kernel_exposes_the_expected_nullspace() {
        // x1 + x2 = 0 and x2 = 0 -> kernel {0}
        let rows = vec![vec![1.0, 1.0], vec![0.0, 1.0]];
        assert!(real_kernel(&rows, 2).is_empty());
        // x1 + x2 = 0 -> one-dimensional kernel spanned by (1, -1)
        let rows = vec![vec![1.0, 1.0]];
        let kernel = real_kernel(&rows, 2);
        assert_eq!(kernel.len(), 1);
        // The basis vector spans (1, -1); its normalization is not pinned.
        assert!(kernel[0][0].abs() > 1e-9);
        assert!((kernel[0][0] + kernel[0][1]).abs() < 1e-9);
    }

    #[test]
    fn a_normalized_intertwiner_rejects_a_non_scalar_gram() {
        // Columns (1, 0) and (0, 1) are orthogonal; the identity is normalized.
        let identity = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        assert!(normalize_intertwiner(&identity, 2).is_some());
        // Columns (1, 1) and (0, 1) are not orthogonal -> rejected.
        let skewed = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        assert!(normalize_intertwiner(&skewed, 2).is_none());
    }

    /// The square rescaling must reproduce `D(a0²)` for a **non-unit** ratio.
    ///
    /// `base = [[1, -√3], [√3, 1]]` has `base² = 4·R120`, so `ratio = 1/4` and
    /// the correct coset matrix is `base/2` (`X² = R120`).  Dividing by `√ratio`
    /// instead gives `2·base`, whose square is `16·R120` — the discriminating
    /// case this fixture pins.
    #[test]
    fn the_coset_scale_reproduces_the_square_for_a_non_unit_ratio() {
        let root3 = 3.0f64.sqrt();
        let base = [1.0, -root3, root3, 1.0];
        let base_square = matrix_product(&base, &base, 2);
        let r120 = [-0.5, -root3 / 2.0, root3 / 2.0, -0.5];
        for index in 0..4 {
            assert!(close(base_square[index], 4.0 * r120[index]));
        }
        let magnitude =
            coset_scale(&base_square, &r120, 2).expect("ratio 1/4 is a positive multiple");
        assert!(close(magnitude, 0.5));
        let x: Vec<f64> = base.iter().map(|value| magnitude * value).collect();
        let x_square = matrix_product(&x, &x, 2);
        for index in 0..4 {
            assert!(close(x_square[index], r120[index]), "X² must be D(a0²)");
        }
        assert!(is_orthogonal(&x, 2));
        // The inverted convention (`base/√ratio = 2·base`) does not square to it.
        let inverted: Vec<f64> = base.iter().map(|value| value / magnitude).collect();
        let inverted_square = matrix_product(&inverted, &inverted, 2);
        assert!(!matrices_close(&inverted_square, &r120));
    }

    #[test]
    fn the_coset_scale_rejects_zero_slot_mismatch_and_non_finite_inputs() {
        // A vanishing entry of base² must not hide a non-vanishing D(a0²).
        assert!(coset_scale(&[0.0, 1.0, 0.0, 1.0], &[1.0, 0.0, 0.0, 0.0], 2).is_none());
        // Non-finite entries and non-positive ratios are refusals.
        assert!(coset_scale(&[f64::NAN, 1.0, 0.0, 1.0], &[1.0, 1.0, 0.0, 1.0], 2).is_none());
        assert!(coset_scale(&[1.0, 1.0, 0.0, 1.0], &[f64::INFINITY, 1.0, 0.0, 1.0], 2).is_none());
        assert!(coset_scale(&[1.0, 1.0, 0.0, 1.0], &[-1.0, -1.0, 0.0, -1.0], 2).is_none());
        // A consistent unit ratio still passes.
        assert!(coset_scale(&[2.0, 0.0, 0.0, 2.0], &[0.5, 0.0, 0.0, 0.5], 2).is_some());
    }

    /// The success guards are positive-finite predicates: a `NaN`/infinite value
    /// is a refusal, while a legal zero residual passes.
    #[test]
    fn the_finite_guards_reject_nan_and_infinity_but_accept_zero_residual() {
        assert!(close(0.0, 0.0));
        assert!(close(1.0, 1.0 + 1e-16));
        assert!(!close(f64::NAN, 0.0));
        assert!(!close(0.0, f64::NAN));
        assert!(!close(f64::INFINITY, f64::INFINITY));
        assert!(small(0.0, 1e-8));
        assert!(!small(f64::NAN, 1e-8));
        assert!(!small(f64::INFINITY, 1e-8));
        assert!(!is_orthogonal(&[f64::NAN, 0.0, 0.0, 1.0], 2));
        assert!(!is_orthogonal(&[f64::INFINITY, 0.0, 0.0, 1.0], 2));
        // A zero intertwiner has a non-positive Gram scalar, a NaN one is not
        // finite: both are refusals, and the identity still normalizes.
        let zero = vec![Complex64::new(0.0, 0.0); 2];
        assert!(normalize_intertwiner(&zero, 1).is_none());
        let nan = vec![Complex64::new(f64::NAN, 0.0), Complex64::new(0.0, 0.0)];
        assert!(normalize_intertwiner(&nan, 1).is_none());
        let unit = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
        assert!(normalize_intertwiner(&unit, 1).is_some());
        let nan_copy = vec![Complex64::new(f64::NAN, 0.0), Complex64::new(0.0, 0.0)];
        assert!(!copies_orthogonal(&nan_copy, &unit, 1));
    }

    // -----------------------------------------------------------------
    // R32 restricted paired type-C tests
    // -----------------------------------------------------------------

    /// The two frozen paired witnesses, built through the public input path.
    fn type_c_witness_model(condensing: &'static str, eta: RealSign) -> MatrixGammaModel {
        use crate::irrep::subduction::Mat3R;
        let request = MatrixGammaRequest::new(
            SourceIdentity {
                parent_sg: 89,
                label: "GM5",
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: 89,
                condensing_label: condensing,
                convention: LabelConvention::Cdml,
                direction_index: 0,
            },
            Some(eta),
        )
        .with_embedding(EmbeddingEvidence {
            map: Mat3R::from_ints([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            shift: [Rat::ZERO; 3],
            lattice_policy: LatticePolicy::RequireExact,
        });
        verify_matrix_gamma_input(&request).expect("the frozen paired witness model verifies")
    }

    struct TypeCWitness {
        model: MatrixGammaModel,
        context: TargetContext,
        corep: MagneticCorepSummary,
        pair: TypeCPair,
        table: TargetTable,
    }

    fn type_c_witness(condensing: &'static str) -> TypeCWitness {
        let model = type_c_witness_model(condensing, RealSign::Plus);
        let uni = model.uni();
        let partial =
            magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni).expect("the UNI summary");
        let gamma = partial
            .summary
            .kpoints
            .iter()
            .find(|kpoint| kpoint.coords == GAMMA)
            .expect("the Gamma k-point");
        let class_columns = magnetic_subduction::proven_gamma_class_columns(
            uni,
            gamma,
            &model.classes,
            &model.lattice,
        )
        .expect("the raw provenance");
        let context =
            TargetContext::build(uni, gamma, &model, &class_columns).expect("the target context");
        let corep = gamma
            .coreps
            .iter()
            .find(|corep| corep.corep_type == CorepType::C && corep.dim == 2)
            .expect("the fixed witness has a type-C dim-2 corep")
            .clone();
        let pair = match screen_type_c_pair(&corep, &context, uni).expect("screenable data") {
            Ok(pair) => pair,
            Err(reason) => panic!("the fixed witness must screen: {reason}"),
        };
        let first_anti = context
            .classes
            .iter()
            .position(|class| class.time_reversal)
            .expect("an antiunitary class");
        let table = build_type_c_table(&pair, &context, uni, &corep.label, first_anti)
            .expect("the fixed witness target table");
        TypeCWitness {
            model,
            context,
            corep,
            pair,
            table,
        }
    }

    /// Type C is the complex (Schur-C) corep type: the real-linear Hom space of
    /// the target with itself is **two**-dimensional, generated by `I` and
    /// `diag(i, -i)`.  Both are full-rank isometries, but they are not two
    /// independent copies: the copy selection keeps one because their cross
    /// Gram block does not vanish.  The production multiplicity is the copy
    /// count, never the kernel dimension.
    #[test]
    fn the_type_c_hom_space_is_two_dimensional_but_the_multiplicity_is_one() {
        for condensing in ["GM2", "GM4"] {
            let witness = type_c_witness(condensing);
            let target = RealCorep::from_target(&witness.table.operators);
            let (rows, width) = hom_rows(&target, &target);
            let kernel = real_kernel(&rows, width);
            assert_eq!(
                kernel.len(),
                2,
                "the type-C self-Hom space must be two-dimensional, not the copy count"
            );
            // Layout: (row·2 + column)·2 + {0 = re, 1 = im}.
            let identity = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
            let f = vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0];
            for vector in [&identity, &f] {
                for row in &rows {
                    let residual: f64 = row
                        .iter()
                        .zip(vector.iter())
                        .map(|(left, right)| left * right)
                        .sum();
                    assert!(
                        residual.abs() < 1e-9,
                        "I and diag(i, -i) must satisfy every self-Hom equation: {residual}"
                    );
                }
            }
            let source = RealCorep::from_source(&witness.model.actions);
            let (kernel, copies) = matching_copies(&source, &witness.table);
            assert_eq!(kernel.len(), 2);
            assert_eq!(copies.len(), 1, "one independent complex copy, not two");
        }
    }

    /// Both fixed witnesses have a real `±1` square sign and the canonical
    /// off-diagonal coset matrix squares to it (the `-I` case is accepted).
    #[test]
    fn the_type_c_coset_matrix_square_matches_the_actual_a0_square() {
        // UNI 749 has D(a0²) = +I; UNI 751 is the real -I branch.
        for (condensing, expected_sign) in [("GM2", 1i8), ("GM4", -1i8)] {
            let witness = type_c_witness(condensing);
            let first_anti = witness
                .context
                .classes
                .iter()
                .position(|class| class.time_reversal)
                .expect("anti class");
            let square = product_class(
                &witness.context.classes[first_anti],
                &witness.context.classes[first_anti],
                &witness.context.classes,
                &witness.context.lattice,
            )
            .expect("a0 squared is a class");
            let s = witness.pair.characters[0][square];
            assert_eq!(witness.pair.characters[1][square], s);
            assert_eq!(s, expected_sign);
            let anti = witness
                .table
                .operators
                .iter()
                .position(|operator| operator.is_antiunitary())
                .expect("an antiunitary operator");
            let x = witness.table.operators[anti].entries();
            let product = matrix_product(x, x, 2);
            let expected = [f64::from(s), 0.0, 0.0, f64::from(s)];
            assert!(
                matrices_close(&product, &expected),
                "the canonical X must square to s·I (s = {s})"
            );
        }
    }

    /// The actual frozen H pairs satisfy `χ1(a0⁻¹ h a0) = χ2(h)` and
    /// `χ2(a0⁻¹ h a0) = χ1(h)` on every unitary class, with separating classes
    /// present (the exchange is not vacuous).
    #[test]
    fn the_actual_witnesses_satisfy_the_a_twisted_pairing() {
        for condensing in ["GM2", "GM4"] {
            let witness = type_c_witness(condensing);
            let first_anti = witness
                .context
                .classes
                .iter()
                .position(|class| class.time_reversal)
                .expect("anti class");
            let identity = find_identity_class(&witness.context.classes).expect("identity");
            let inverse = inverse_class(
                &witness.context.classes,
                &witness.context.lattice,
                first_anti,
                identity,
            )
            .expect("a0 inverse");
            let mut checked = 0usize;
            let mut discriminating = 0usize;
            for class in 0..witness.context.classes.len() {
                if witness.context.classes[class].time_reversal {
                    continue;
                }
                let left = compose(
                    &witness.context.classes[inverse],
                    &witness.context.classes[class],
                )
                .expect("a0^-1 h");
                let conjugated =
                    compose(&left, &witness.context.classes[first_anti]).expect("a0^-1 h a0");
                let sigma = find_model_class(
                    &conjugated,
                    &witness.context.classes,
                    &witness.context.lattice,
                )
                .expect("sigma is a class");
                checked += 1;
                if witness.pair.characters[0][class] != witness.pair.characters[1][class] {
                    discriminating += 1;
                }
                assert_eq!(
                    witness.pair.characters[0][sigma], witness.pair.characters[1][class],
                    "{condensing}: χ1(σ(h)) != χ2(h) at class {class}"
                );
                assert_eq!(
                    witness.pair.characters[1][sigma], witness.pair.characters[0][class],
                    "{condensing}: χ2(σ(h)) != χ1(h) at class {class}"
                );
            }
            assert!(checked >= 4);
            assert!(
                discriminating > 0,
                "{condensing}: the two H characters must be discriminated on at least one class"
            );
        }
    }

    #[test]
    fn type_c_character_classification_names_finite_and_non_finite_problems() {
        let tol = GRID_TOLERANCE;
        assert_eq!(type_c_real_sign(1.0), Ok(1));
        assert_eq!(type_c_real_sign(-1.0), Ok(-1));
        assert_eq!(type_c_real_sign(0.5), Err(TypeCCharacterProblem::OffGrid));
        assert_eq!(
            type_c_real_sign(f64::NAN),
            Err(TypeCCharacterProblem::NonFinite)
        );
        assert_eq!(
            type_c_real_sign(f64::INFINITY),
            Err(TypeCCharacterProblem::NonFinite)
        );
        assert_eq!(
            type_c_real_sign(f64::NEG_INFINITY),
            Err(TypeCCharacterProblem::NonFinite)
        );
        assert_eq!(
            type_c_complex_sign(Complex64::new(0.0, tol * 10.0)),
            Err(TypeCCharacterProblem::NonReal)
        );
        assert_eq!(type_c_complex_sign(Complex64::new(1.0, tol * 0.5)), Ok(1));
        assert_eq!(
            type_c_complex_sign(Complex64::new(f64::NAN, 0.0)),
            Err(TypeCCharacterProblem::NonFinite)
        );
        assert_eq!(
            type_c_complex_sign(Complex64::new(0.0, f64::NAN)),
            Err(TypeCCharacterProblem::NonFinite)
        );
    }

    #[test]
    fn the_same_h_pair_is_unresolved_not_no_candidate() {
        let witness = type_c_witness("GM2");
        let mut fake = witness.corep.clone();
        fake.source_irreps[1] = fake.source_irreps[0].clone();
        match screen_type_c_pair(&fake, &witness.context, 749) {
            Err(MatrixTypeCError::EquivalentHSources {
                uni, first, second, ..
            }) => {
                assert_eq!(uni, 749);
                assert_eq!(first, "GM3");
                assert_eq!(second, "GM3");
            }
            other => panic!("expected EquivalentHSources, got {other:?}"),
        }
    }

    #[test]
    fn type_c_metadata_boundaries_are_named_and_not_silently_skipped() {
        let witness = type_c_witness("GM2");
        // Spinor and dimension are finite-domain exclusions.
        let mut spinor = witness.corep.clone();
        spinor.source_irreps[0].spinor = true;
        match screen_type_c_pair(&spinor, &witness.context, 749) {
            Ok(Err(reason)) => assert_eq!(reason, "spinor(0)"),
            other => panic!("expected a spinor exclusion, got {other:?}"),
        }
        let mut dimension = witness.corep.clone();
        dimension.source_irreps[0].dim = 2;
        match screen_type_c_pair(&dimension, &witness.context, 749) {
            Ok(Err(reason)) => assert_eq!(reason, "dimension(0:2)"),
            other => panic!("expected a dimension exclusion, got {other:?}"),
        }
        // A summary H source from another space group is an identity
        // inconsistency, not a mathematical negative.
        let mut wrong_group = witness.corep.clone();
        wrong_group.source_irreps[0].sg = 1;
        match screen_type_c_pair(&wrong_group, &witness.context, 749) {
            Err(MatrixTypeCError::HSource { reason, .. }) => {
                assert!(reason.contains("unitary subgroup"), "{reason}")
            }
            other => panic!("expected an HSource identity refusal, got {other:?}"),
        }
    }

    #[test]
    fn a_twisted_character_tamper_is_refused_by_the_twist_guard() {
        let witness = type_c_witness("GM2");
        let first_anti = witness
            .context
            .classes
            .iter()
            .position(|class| class.time_reversal)
            .expect("anti class");
        let mut saw_twist = false;
        for class in 0..witness.context.classes.len() {
            if witness.context.classes[class].time_reversal
                || witness.pair.characters[0][class] == witness.pair.characters[1][class]
            {
                continue;
            }
            let mut tampered = witness.pair.clone();
            tampered.characters[0][class] = witness.pair.characters[1][class];
            tampered.characters[1][class] = witness.pair.characters[0][class];
            let built = build_type_c_table(
                &tampered,
                &witness.context,
                749,
                &witness.corep.label,
                first_anti,
            );
            if let Err(MatrixTypeCError::Inconsistent { stage: "twist", .. }) = built {
                saw_twist = true;
                break;
            }
        }
        assert!(
            saw_twist,
            "a swapped character at a separating class must break the twisted pairing"
        );
    }

    #[test]
    fn a_single_anti_entry_tamper_breaks_the_law_and_the_reconstruction() {
        let witness = type_c_witness("GM2");
        let mut tampered = witness.table.clone();
        let anti = tampered
            .operators
            .iter()
            .position(|operator| operator.is_antiunitary())
            .expect("an antiunitary operator");
        tampered.operators[anti].entries[1] = 0.0;
        assert!(
            check_target_law(&tampered, &witness.context).is_err(),
            "a single tampered off-diagonal entry must break the ordered-pair law"
        );
        let source = RealCorep::from_source(&witness.model.actions);
        let (_, copies) = matching_copies(&source, &witness.table);
        assert_eq!(copies.len(), 1);
        assert!(
            verify_reconstruction_data(
                &witness.model,
                &witness.corep.label,
                2,
                &copies,
                &tampered.operators,
            )
            .is_err(),
            "the tampered operator table must fail the complex reconstruction"
        );
    }

    #[test]
    fn another_antiunitary_representative_builds_an_equivalent_target() {
        let witness = type_c_witness("GM2");
        let anti_classes: Vec<usize> = witness
            .context
            .classes
            .iter()
            .enumerate()
            .filter(|(_, class)| class.time_reversal)
            .map(|(index, _)| index)
            .collect();
        assert!(anti_classes.len() > 1);
        let source = RealCorep::from_source(&witness.model.actions);
        for a0 in &anti_classes {
            let table = build_type_c_table(
                &witness.pair,
                &witness.context,
                749,
                &witness.corep.label,
                *a0,
            )
            .expect("every antiunitary representative builds a lawful target");
            check_target_law(&table, &witness.context).expect("the target law");
            let (_, copies) = matching_copies(&source, &table);
            assert_eq!(copies.len(), 1);
            assert!(normalize_intertwiner(&copies[0], 2).is_some());
        }
    }

    #[test]
    fn a_wrong_raw_class_does_not_fall_back() {
        let mut model = type_c_witness_model("GM2", RealSign::Plus);
        assert_eq!(model.uni(), 749);
        model.classes[0].rotation[0][0] = 2;
        match subduce_matrix_gamma_type_c(&model) {
            Err(MatrixTypeCError::TargetProvenance { uni, .. }) => assert_eq!(uni, 749),
            other => panic!("expected the raw provenance refusal, got {other:?}"),
        }
    }
}
