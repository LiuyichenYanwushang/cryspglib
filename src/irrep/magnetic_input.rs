//! Round 25 card 1: the **input contract** of a scalar Γ magnetic subduction.
//!
//! Milestone `docs/subduction-next-milestones.md` R10 ("先实现 Γ 磁共表示分导")
//! prescribes the order of work: first the unitary restriction, then the
//! antiunitary pairing "按母群表示真实的反幺正作用", and only then one scalar Γ
//! example with an explicit time parity.  R10 explicitly forbids replacing that
//! action with a summary classification: the classification is computed *from*
//! an action, it is not a substitute for one.
//!
//! This module delivers **only that input layer**, and deliberately not the
//! subduction itself.  It turns one verified embedding into the exact
//! corepresentation a scalar Γ magnetic subduction starts from, and adds no
//! matrix machinery or algorithm.
//!
//! # The model, in one paragraph
//!
//! Write `P` for the ordinary parent space group and let
//!
//! ```text
//! P × {1, Θ}
//! ```
//!
//! be its **grey** group, `Θ` the pure time reversal.  The magnetic subgroup `M`
//! (Type I, III or IV) is embedded into that grey parent by the verified map
//! `φ`, and `M` is **not** required to contain `Θ`.  The parent's action on the
//! fixed real one-dimensional basis of the parent's Γ space is
//!
//! ```text
//! ρ(p, Θ^t) = χ(p) · η^t · K^t,     χ(p) = ±1,   η = ±1,   K = complex conjugation
//! ```
//!
//! with `Θ² = +1` by this fixed-basis construction, and the model records the
//! **restriction of `ρ` to `φ(M)`**:
//!
//! ```text
//! A(g) = η^{t(g)} · χ(π(g))
//! ```
//!
//! where `π(g)` is the parent operation the class maps to and `t(g)` its
//! time-reversal flag.  Both `χ` and the images are **computed** — `χ` from the
//! parent's tabulated scalar Γ row, the images from the embedding.  The caller
//! states the embedding, the parent source and `η`.
//!
//! This explicit parent action is why no carrier class is chosen: `η` is a datum
//! of the *parent* action, `A(g)` is its restriction evaluated at `g`, and which
//! listed class of `M` happens to be a pure `Θ` (if any) is irrelevant.  A grey
//! parent always contains `Θ`; a magnetic subgroup of Type III/IV need not, and
//! its antiunitary classes inherit `η` through `φ`.
//!
//! The family is the **explicitly real `±1`** one.  That is a restriction of this
//! contract, not a property of scalar irreps in general: a one-dimensional
//! character can be complex (`χ = ±i`), and such a source is refused by name.
//!
//! # What is checked
//!
//! `A` must be the restriction of a corepresentation of the grey parent:
//! for every ordered pair of classes,
//!
//! ```text
//! A(g) · [ conj(A(h))  if t(g) else  A(h) ] = A(g·h)
//! ```
//!
//! The bracket is decided by the time-reversal flag of the **left** factor, and
//! the unitary branch keeps `A(h)` unchanged.  Every value is `±1` and
//! `conj(±1) = ±1`, so the conjugation is the identity on the values while the
//! branch is structural.
//!
//! As the special case `g = h = Θ`, `A(Θ)² = A(Θ²) = +1`, because `Θ² = +1` is
//! built into the parent action.  **Time parity is not the square**: a time-odd
//! scalar (`η = -1`) is reported as [`TimeParity::Odd`] while
//! [`ScalarGammaModel::time_reversal_square`] stays `+1`, and the model never
//! reports a time-odd scalar as `T² = -1`.
//!
//! Both signs of `η` are legal.  Flipping `η` is not a tampering the module
//! rejects: it flips exactly the antiunitary classes' values, and the tests check
//! that (plus that the unitary half is untouched).
//!
//! # What is verified rather than trusted
//!
//! * the operation set is a group and matches a fresh database read
//!   ([`crate::irrep::magnetic_embedding::verify_group`],
//!   [`crate::irrep::magnetic_embedding::verify_against_database`]);
//! * the **record identity** is read fresh from `(parent, condensing label,
//!   convention, row index)`, never handed in: the row is looked up in the
//!   resolved irrep's own magnetic isotropy table and the geometry is derived from
//!   that record (UNI, basis, origin), so there is nothing for a caller to
//!   cross-wire;
//! * the **source identity** is separate from the condensing one: the subduced
//!   parent irrep is named by its own label and convention and must be a real
//!   `±1` scalar Γ irrep, while the condensing record may sit anywhere;
//! * the caller's `(map, shift)` embeds the group into the grey parent
//!   ([`crate::irrep::magnetic_embedding::verify_embedding`]); under
//!   [`LatticePolicy::RequireExact`] the exact record-lattice equality is
//!   additionally recomputed with
//!   [`crate::irrep::magnetic_embedding::realises_record_lattice_checked`], which
//!   returns `Result` so an exact-arithmetic limit is an error and never a silent
//!   `false`.  The exact verdict is always computed and carried by the model,
//!   whatever the policy.
//!
//! [`ScalarGammaModel`]'s fields are private and it is not `Clone`, so the only
//! way to obtain one is this verification; it carries the verified `(map, shift)`,
//! both identities, the lattice verdict and the action table for the subduction
//! layer to use.
//!
//! # Scope fences
//!
//! * **Real, one-dimensional, non-spinor Γ source.**  A complex character, a
//!   spinor record, a source k-point other than Γ and a source dimension other
//!   than 1 each have their own refusal.  The real-`±1` limit is a restriction of
//!   this contract (scalar characters can be complex); no scalar Γ row of the
//!   shipped tables is complex, which a test scans and reports.
//! * **No uniqueness gate and no general equivalence claim.**  Several verified
//!   embeddings of one record are legitimate; the caller selects one, it is
//!   checked, and the model records it.  The action table is a function of `χ`
//!   and the mapped classes.
//! * **No subduction.**  The ordinary entry points are untouched.

use crate::irrep::isotropy;
use crate::irrep::labels::LabelConvention;
use crate::irrep::magnetic_embedding::{
    self, GeometryError, MagneticOperation, MagneticOperationSet, verify_embedding,
};
use crate::irrep::subduction::{Mat3I, Mat3R, Rat, Vec3R};
use crate::irrep::types::{CharacterRow, IrrepRecord, SeitzOperation};

/// Denominators a parent value is snapped to.
///
/// A one-dimensional unitary character takes values in the roots of unity; the
/// corpus's scalar little-group characters sit on denominators up to twelve, the
/// same grid the module's translation contract uses.  A value off the grid is an
/// error, never a rounded value.
const MAX_EXACT_DENOMINATOR: i128 = 12;

/// Tolerance of the grid snap.  The parent row is stored as `f64`; this decides
/// which exact rational it stands for, and the exact value is what is then used.
const GRID_TOLERANCE: f64 = 1e-9;

const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// The isotropy tables are queried by the irrep's own CDML label, whatever
/// spelling the caller used to name it.
const CDML: LabelConvention = LabelConvention::Cdml;

/// The scalar the pure time reversal acts with on the fixed real one-dimensional
/// basis: `Θ = η·K`, with `K` complex conjugation, hence `Θ² = +1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RealSign {
    /// `η = +1`.
    Plus,
    /// `η = -1`.
    Minus,
}

impl RealSign {
    /// `+1` or `-1` as an integer.
    pub fn value(self) -> i32 {
        match self {
            Self::Plus => 1,
            Self::Minus => -1,
        }
    }
}

/// The parity of the time reversal on the parent's Γ space: `η = +1` is even,
/// `η = -1` is odd.
///
/// This is the parity of the **action**.  In this model `Θ² = +1`, so the parity
/// is not the sign of `T²`: the model never reports a time-odd scalar as
/// `T² = -1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TimeParity {
    /// `η = +1`.
    Even,
    /// `η = -1`.
    Odd,
}

/// One class of the magnetic group with the `±1` the verified model assigns it.
#[derive(Clone, Debug)]
pub struct ClassAction {
    /// The class representative, as the database lists the set.
    pub class: MagneticOperation,
    /// `A(class)`: `χ(π(class))` for a unitary class, `η·χ(π(class))` for an
    /// antiunitary one.
    pub value: RealSign,
}

/// The identity of the parent representation being subduced.
///
/// This is **not** the condensing record: the record names an isotropy subgroup
/// and may sit at a line or plane, while this names the parent irrep whose
/// character the model uses.  Keeping them apart is what lets the source be a Γ
/// irrep while the condensing record is not.
#[derive(Clone, Copy, Debug)]
pub struct SourceIdentity<'a> {
    /// Parent space group number (1–230).
    pub parent_sg: u8,
    /// Parent irrep label, read under [`Self::convention`].
    pub label: &'a str,
    /// Which spelling of the labels the caller uses.
    pub convention: LabelConvention,
}

/// The identity of the magnetic isotropy record the embedding is stated against.
#[derive(Clone, Copy, Debug)]
pub struct RecordIdentity<'a> {
    /// Parent space group number (1–230); must equal the source's.
    pub parent_sg: u8,
    /// Label of the **condensing** parent irrep, read under [`Self::convention`].
    ///
    /// This convention is the record's own: the source and the record are
    /// independent requests and may be named in different spellings.  The row is
    /// obtained by entering the isotropy query with the resolved irrep's canonical
    /// CDML label, so either spelling reaches the same record.
    pub condensing_label: &'a str,
    /// Which spelling of this label the caller uses.
    pub convention: LabelConvention,
    /// Position of the row in that irrep's magnetic isotropy table.
    pub direction_index: usize,
}

/// What the caller requires of the record's lattice under the selected map.
///
/// The weaker policy is not a second geometric condition: embedding the group in
/// the parent is already proved by
/// [`crate::irrep::magnetic_embedding::verify_embedding`] for every request, so
/// `ContainmentOnly` means "that containment is enough" and adds nothing.  The
/// exact record-lattice equality is an **extra** requirement the caller may ask
/// for with `RequireExact`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LatticePolicy {
    /// Embedding in the parent is enough; the record's lattice may be a proper
    /// superlattice of the map's image.
    ContainmentOnly,
    /// Additionally require `A·L_M = L_record` exactly.
    RequireExact,
}

/// The caller's chosen embedding.
#[derive(Clone, Copy, Debug)]
pub struct EmbeddingEvidence {
    /// The map applied to the magnetic group's coordinates.
    pub map: Mat3R,
    /// The exact origin shift, in the parent's conventional coordinates.
    pub shift: [Rat; 3],
    /// What the caller requires of the record's lattice.
    pub lattice_policy: LatticePolicy,
}

impl EmbeddingEvidence {
    /// The identity map, a zero shift and the containment policy.
    pub fn identity() -> Self {
        Self {
            map: Mat3R::identity(),
            shift: [Rat::ZERO; 3],
            lattice_policy: LatticePolicy::ContainmentOnly,
        }
    }

    /// The same embedding with a different lattice policy.
    pub fn with_lattice_policy(mut self, policy: LatticePolicy) -> Self {
        self.lattice_policy = policy;
        self
    }
}

/// Everything a caller must state to obtain a verified scalar Γ model.
pub struct ScalarGammaRequest<'a> {
    source: SourceIdentity<'a>,
    record: RecordIdentity<'a>,
    embedding: EmbeddingEvidence,
    time_reversal: Option<RealSign>,
}

impl<'a> ScalarGammaRequest<'a> {
    /// A request for one source representation and one condensing record.
    ///
    /// `time_reversal` is the parent action's `η`.  It is required for a magnetic
    /// subgroup with an antiunitary class and may be omitted for a purely unitary
    /// (Type I) subgroup, whose restriction to the unitary part does not use it;
    /// a stated value is preserved either way.
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

    /// The embedding the caller selected.
    pub fn with_embedding(mut self, embedding: EmbeddingEvidence) -> Self {
        self.embedding = embedding;
        self
    }

    /// The source representation identity.
    pub fn source(&self) -> SourceIdentity<'a> {
        self.source
    }

    /// The condensing record identity.
    pub fn record(&self) -> RecordIdentity<'a> {
        self.record
    }

    /// The selected embedding.
    pub fn embedding(&self) -> EmbeddingEvidence {
        self.embedding
    }

    /// The requested time-reversal sign, if any.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.time_reversal
    }
}

/// Why a scalar Γ magnetic input was refused.
///
/// Every variant is a specific refusal: the cases outside this round's scope
/// each have their own name and none is reported as another.
#[derive(Debug)]
pub enum MagneticInputError {
    /// A parent space group number is outside 1–230.
    UnknownParentSpaceGroup { sg: u8 },
    /// No irrep of that parent has this label under this convention.
    LabelNotFound {
        sg: u8,
        label: String,
        convention: LabelConvention,
    },
    /// The label exists but does not name a unique irrep.
    AmbiguousLabel { sg: u8, label: String },
    /// The source irrep is double-valued.
    SpinorSource { sg: u8, label: &'static str },
    /// The source irrep is not one-dimensional.
    NotScalarSource {
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
    /// The source and the record name different parents.
    ParentMismatch { source_sg: u8, record_sg: u8 },
    /// The record's direction index is outside the condensing irrep's table.
    DirectionOutOfRange {
        sg: u8,
        label: String,
        index: usize,
        len: usize,
    },
    /// The operation set does not match a fresh database read.
    DatabaseMismatch(String),
    /// The set does not pass the group-axiom check.
    GroupProof(String),
    /// The class list and the set's own class count disagree.
    ClassListInconsistent { classes: usize, recorded: usize },
    /// A value is outside the exact grid the contract represents.
    OffGrid { sg: u8, rotation: Mat3I, value: f64 },
    /// A parent operation's translation is outside that grid.
    ParentTranslationOffGrid { sg: u8, rotation: Mat3I, value: f64 },
    /// A character value is not the exact real `±1` this family allows.
    NonRealCharacter {
        sg: u8,
        rotation: Mat3I,
        real: f64,
        imaginary: f64,
    },
    /// The identity class's character is not `+1`.
    IdentityCharacterNotOne { sg: u8, value: f64 },
    /// A mapped class has no entry in the source's character row, i.e. the row
    /// and the parent's operation table describe different settings.
    ///
    /// The parent table ([`crate::irrep::query::symmetry_operations_of`]) and the
    /// source row (the PIR operation storage behind
    /// [`crate::irrep::types::IrrepRecord::ordinary_scalar_selected_arm_block_trace`])
    /// are generated in their own settings.  They agree operation by operation
    /// modulo the parent's lattice for 828 of the 980 scalar Γ rows; for the
    /// other 152 (24 space groups, e.g. SG 48 and SG 227) they do not, while
    /// their **rotation sets are equal** in all 980 rows.  This round has not
    /// established how to align the two tables, so such a row is refused instead
    /// of being read through a guessed pairing; the message claims no more than
    /// that, and in particular does not claim that no alignment exists.
    CharacterRowSettingMismatch { sg: u8, rotation: Mat3I, uni: usize },
    /// The caller's `(map, shift)` does not embed the group, or the exact check
    /// hit an arithmetic limit.
    EmbeddingRejected(GeometryError),
    /// The selected map does not satisfy the requested lattice policy; the
    /// predicate was recomputed and this is a real failure.
    LatticePolicyUnsatisfied { policy: LatticePolicy },
    /// The set has an antiunitary coset but no `η` was supplied.
    MissingTimeReversal { uni: usize },
    /// The group law failed for a specific pair of classes.
    NotACorepresentation {
        uni: usize,
        left: usize,
        right: usize,
        expected: i32,
        found: i32,
    },
    /// Exact arithmetic failed inside the checker.
    Arithmetic(String),
}

impl std::fmt::Display for MagneticInputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownParentSpaceGroup { sg } => {
                write!(f, "parent space group {sg} is outside 1-230")
            }
            Self::LabelNotFound {
                sg,
                label,
                convention,
            } => write!(
                f,
                "SG {sg} has no irrep with label {label:?} under {convention:?}"
            ),
            Self::AmbiguousLabel { sg, label } => {
                write!(f, "SG {sg} has several irreps for label {label:?}")
            }
            Self::SpinorSource { sg, label } => write!(
                f,
                "SG {sg} source irrep {label:?} is double-valued (spinor); this contract covers \
                 the scalar case only"
            ),
            Self::NotScalarSource {
                sg,
                label,
                dimension,
            } => write!(
                f,
                "SG {sg} source irrep {label:?} has dimension {dimension}; this contract covers \
                 dimension 1 only"
            ),
            Self::SourceNotGamma {
                sg,
                label,
                numerators,
                denominator,
            } => write!(
                f,
                "SG {sg} source irrep {label:?} sits at k = {numerators:?}/{denominator}, not at \
                 Gamma"
            ),
            Self::ParentMismatch {
                source_sg,
                record_sg,
            } => write!(
                f,
                "the source's parent is SG {source_sg} but the record's is SG {record_sg}"
            ),
            Self::DirectionOutOfRange {
                sg,
                label,
                index,
                len,
            } => write!(
                f,
                "SG {sg} condensing irrep {label:?} has {len} magnetic isotropy rows; index \
                 {index} is out of range"
            ),
            Self::DatabaseMismatch(detail) => {
                write!(f, "the database read or comparison failed: {detail}")
            }
            Self::GroupProof(detail) => write!(f, "the operation set is not a group: {detail}"),
            Self::ClassListInconsistent { classes, recorded } => write!(
                f,
                "the set lists {classes} classes but records {recorded}"
            ),
            Self::OffGrid { sg, rotation, value } => write!(
                f,
                "SG {sg} value at rotation {rotation:?} is {value}, which is not on the exact grid \
                 the contract represents"
            ),
            Self::ParentTranslationOffGrid {
                sg,
                rotation,
                value,
            } => write!(
                f,
                "SG {sg}: translation {value} of the operation with rotation {rotation:?} is not \
                 on the exact grid the contract represents"
            ),
            Self::NonRealCharacter {
                sg,
                rotation,
                real,
                imaginary,
            } => write!(
                f,
                "SG {sg} character at rotation {rotation:?} is {real} + {imaginary}i; this \
                 contract covers the real +-1 characters only"
            ),
            Self::IdentityCharacterNotOne { sg, value } => write!(
                f,
                "SG {sg} character at the identity class is {value}, not one"
            ),
            Self::CharacterRowSettingMismatch {
                sg,
                rotation,
                uni,
            } => write!(
                f,
                "SG {sg} (UNI {uni}): the source character row and the parent operation table do \
                 not describe the same setting (rotation {rotation:?} has no row entry modulo the \
                 parent lattice); the setting transform is not part of this contract"
            ),
            Self::EmbeddingRejected(error) => {
                write!(f, "the selected (map, shift) is not an embedding: {error}")
            }
            Self::LatticePolicyUnsatisfied { policy } => write!(
                f,
                "the selected map does not satisfy the requested record-lattice policy {policy:?}"
            ),
            Self::MissingTimeReversal { uni } => write!(
                f,
                "UNI {uni} has an antiunitary coset, so eta must be supplied"
            ),
            Self::NotACorepresentation {
                uni,
                left,
                right,
                expected,
                found,
            } => write!(
                f,
                "UNI {uni}: the product of classes {left} and {right} requires {expected}, but \
                 the table gives {found}"
            ),
            Self::Arithmetic(detail) => write!(f, "exact arithmetic failed: {detail}"),
        }
    }
}

impl std::error::Error for MagneticInputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::EmbeddingRejected(error) => Some(error),
            _ => None,
        }
    }
}

impl From<GeometryError> for MagneticInputError {
    fn from(error: GeometryError) -> Self {
        Self::EmbeddingRejected(error)
    }
}

/// The verified scalar Γ magnetic model: the exact corepresentation a magnetic
/// subduction of this parent irrep starts from, with the evidence it was checked
/// against.
///
/// All fields are private and the type is not `Clone`, so the only construction
/// path is [`verify_scalar_gamma_input`]; every accessor is read-only.
#[derive(Debug)]
pub struct ScalarGammaModel {
    parent_sg: u8,
    source_label: &'static str,
    record_label: &'static str,
    input_convention: LabelConvention,
    direction_index: usize,
    ordinal: usize,
    uni: usize,
    map: Mat3R,
    shift: [Rat; 3],
    lattice_policy: LatticePolicy,
    realises_record_lattice: bool,
    time_reversal: Option<RealSign>,
    actions: Vec<ClassAction>,
    classes: Vec<MagneticOperation>,
    lattice: Vec<[Rat; 3]>,
}

impl ScalarGammaModel {
    /// Parent space group of the source and the record.
    pub fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// The source (subduced) parent irrep label, in its canonical CDML spelling.
    ///
    /// The label and [`Self::source_convention`] always resolve back to the
    /// selected irrep: a request may name a Bradley-Cracknell spelling, and the
    /// model reports the canonical one instead of a spelling that would not
    /// resolve under the convention it is paired with.  What the caller typed is
    /// available as [`Self::input_convention`].
    pub fn source_label(&self) -> &'static str {
        self.source_label
    }

    /// The convention [`Self::source_label`] is written in (`Cdml`).
    pub fn source_convention(&self) -> LabelConvention {
        LabelConvention::Cdml
    }

    /// The condensing irrep label of the record, in its canonical CDML spelling.
    pub fn record_label(&self) -> &'static str {
        self.record_label
    }

    /// The convention [`Self::record_label`] is written in (`Cdml`).
    pub fn record_convention(&self) -> LabelConvention {
        LabelConvention::Cdml
    }

    /// The convention the caller used to name the **source** label.
    ///
    /// Provenance only: there is no field for the record's input spelling, and
    /// the reported labels are canonical, so this accessor is not the convention
    /// that makes them resolvable.
    pub fn input_convention(&self) -> LabelConvention {
        self.input_convention
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

    /// The verified map of the embedding the model was built with.
    pub fn map(&self) -> &Mat3R {
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

    /// Whether the selected map realises the record's lattice exactly, always
    /// recomputed from the map and the record.
    ///
    /// Under [`LatticePolicy::ContainmentOnly`] this may be `false`: the model is
    /// still valid, because the chosen map already embeds the group in the
    /// parent, and the value is reported instead of being hidden.
    pub fn realises_record_lattice(&self) -> bool {
        self.realises_record_lattice
    }

    /// The parent action's `η`, or `None` when the caller stated none.
    ///
    /// For a Type I (purely unitary) magnetic subgroup the restriction never uses
    /// `η`, so a request may omit it; the value is still preserved when given.
    pub fn time_reversal(&self) -> Option<RealSign> {
        self.time_reversal
    }

    /// The parity of `η`, which is `None` exactly when no `η` was stated.
    ///
    /// A Type I subgroup may state one (it is the parent datum) and then reports a
    /// parity of the parent action, so this is not "no antiunitary coset" in
    /// disguise.
    pub fn time_parity(&self) -> Option<TimeParity> {
        self.time_reversal.map(|sign| match sign {
            RealSign::Plus => TimeParity::Even,
            RealSign::Minus => TimeParity::Odd,
        })
    }

    /// The square of the **parent** action's time reversal: `ρ(Θ)² = ρ(Θ²) = +1`.
    ///
    /// This is a property of the parent action on the fixed real basis, so it is
    /// `Some(Plus)` exactly when the caller stated an `η` (Type I included) and
    /// `None` when none was stated.  It is **not** the square of an arbitrary
    /// antiunitary class of the magnetic subgroup: such a class may have a spatial
    /// part of order four, and its square is then a non-identity unitary class,
    /// which the pairwise law checks without a separate claim here.  So a
    /// time-odd scalar (`η = -1`) reports `Odd` with a `+1` square, which is the
    /// case a naive "parity is `T²`" reading would get wrong.
    pub fn time_reversal_square(&self) -> Option<RealSign> {
        self.time_reversal.map(|_| RealSign::Plus)
    }

    /// Every class of the magnetic group with its verified `±1`, in the set's
    /// class order.
    pub fn actions(&self) -> &[ClassAction] {
        &self.actions
    }

    /// The action of one class, by index in the set's class order.
    pub fn action_at(&self, index: usize) -> Option<&ClassAction> {
        self.actions.get(index)
    }

    /// The classes the model was verified on, in the set's class order.
    pub fn classes(&self) -> &[MagneticOperation] {
        &self.classes
    }

    /// The group's translation lattice, as the class comparisons reduced it.
    pub fn translation_lattice(&self) -> &[[Rat; 3]] {
        &self.lattice
    }
}

/// Build and verify the scalar Γ magnetic model for one request.
///
/// This is the whole input contract of R25 card 1: it reads the source and the
/// record fresh from the database, proves the operation set is a group, derives
/// the geometry and compares it with the record, checks the caller's embedding,
/// recomputes the requested lattice policy, reads the parent's real `±1`
/// character, builds `A(g) = η^{t(g)}·χ(π(g))` and checks the corepresentation
/// law on every ordered pair of classes.
pub fn verify_scalar_gamma_input(
    request: &ScalarGammaRequest<'_>,
) -> Result<ScalarGammaModel, MagneticInputError> {
    let source_identity = request.source;
    let record_identity = request.record;
    if source_identity.parent_sg == 0 || source_identity.parent_sg > 230 {
        return Err(MagneticInputError::UnknownParentSpaceGroup {
            sg: source_identity.parent_sg,
        });
    }
    if source_identity.parent_sg != record_identity.parent_sg {
        return Err(MagneticInputError::ParentMismatch {
            source_sg: source_identity.parent_sg,
            record_sg: record_identity.parent_sg,
        });
    }
    let parent_sg = source_identity.parent_sg;
    let source_irrep = resolve_source(source_identity)?;

    // Each identity resolves under its **own** convention (they are independent
    // requests); the isotropy query is then entered with the resolved irrep's own
    // CDML label, so a Bradley-Cracknell spelling works, and the returned row
    // carries its own global ordinal.
    let condensing = resolve_condensing(record_identity)?;
    let rows = isotropy::magnetic_isotropy_subgroups(parent_sg, condensing.ml, CDML)
        .map_err(|error| MagneticInputError::DatabaseMismatch(error.to_string()))?;
    let Some(subgroup) = rows.get(record_identity.direction_index).copied() else {
        return Err(MagneticInputError::DirectionOutOfRange {
            sg: parent_sg,
            label: condensing.ml.to_string(),
            index: record_identity.direction_index,
            len: rows.len(),
        });
    };
    let subgroup_record = subgroup.record;
    let geometry = magnetic_embedding::geometry_of(parent_sg, &subgroup_record)
        .map_err(MagneticInputError::from)?;

    // The operation set comes from the record's own UNI, so a caller cannot hand
    // in operations belonging to another group.
    let set = magnetic_embedding::magnetic_operations(subgroup_record.mag_sg)
        .map_err(|error| MagneticInputError::DatabaseMismatch(error.to_string()))?;
    magnetic_embedding::verify_group(&set)
        .map_err(|error| MagneticInputError::GroupProof(error.to_string()))?;
    magnetic_embedding::verify_against_database(&set)
        .map_err(|error| MagneticInputError::DatabaseMismatch(error.to_string()))?;
    if set.class_count() != set.classes.len() {
        return Err(MagneticInputError::ClassListInconsistent {
            classes: set.classes.len(),
            recorded: set.class_count(),
        });
    }

    let evidence = request.embedding;
    if evidence.map.inverse().is_err() {
        return Err(MagneticInputError::Arithmetic(
            "the selected map is singular".to_string(),
        ));
    }
    verify_embedding(&geometry, &set, evidence.map, evidence.shift)
        .map_err(MagneticInputError::EmbeddingRejected)?;
    // The exact verdict is always recomputed and carried; the policy only decides
    // whether a `false` is a refusal.  A `true` under `ContainmentOnly` is the
    // same proved embedding, not a weaker geometry.
    let realises_record_lattice = magnetic_embedding::realises_record_lattice_checked(
        &geometry,
        &set,
        &evidence.map,
    )?;
    if evidence.lattice_policy == LatticePolicy::RequireExact && !realises_record_lattice {
        return Err(MagneticInputError::LatticePolicyUnsatisfied {
            policy: evidence.lattice_policy,
        });
    }

    // `η` is the parent action's datum, so it is preserved whenever the caller
    // states it, Type I included (the restriction to a unitary subgroup does not
    // use it, which is why a Type I request may omit it).  A magnetic subgroup
    // with an antiunitary class cannot be described without it.
    let time_reversal = match (set.antiunitary_classes, request.time_reversal) {
        (0, given) => given,
        (_, Some(sign)) => Some(sign),
        (_, None) => return Err(MagneticInputError::MissingTimeReversal { uni: set.uni }),
    };

    let row = source_irrep
        .ordinary_scalar_selected_arm_block_trace()
        .map_err(|error| MagneticInputError::DatabaseMismatch(error.to_string()))?;
    let parent = ParentFrame::of(parent_sg)?;
    let mapped = mapped_classes(&set, &evidence, &parent)?;
    let character = read_character(parent_sg, &row, &set, &mapped, &parent)?;

    let actions: Vec<ClassAction> = set
        .classes
        .iter()
        .cloned()
        .zip(character)
        .zip(mapped.iter())
        .map(|((class, chi), mapped_class)| ClassAction {
            class,
            value: if mapped_class.time_reversal {
                combine(chi, time_reversal.expect("an antiunitary class implies eta"))
            } else {
                chi
            },
        })
        .collect();

    let lattice = magnetic_embedding::translation_lattice(&set);
    check_corepresentation_law(&set, &actions, &lattice)?;

    Ok(ScalarGammaModel {
        parent_sg,
        source_label: source_irrep.ml,
        record_label: condensing.ml,
        input_convention: source_identity.convention,
        direction_index: record_identity.direction_index,
        ordinal: subgroup.ordinal,
        uni: set.uni,
        map: evidence.map,
        shift: evidence.shift,
        lattice_policy: evidence.lattice_policy,
        realises_record_lattice,
        time_reversal,
        actions,
        classes: set.classes.clone(),
        lattice,
    })
}

/// Multiply two `±1`.
fn combine(left: RealSign, right: RealSign) -> RealSign {
    if left.value() * right.value() > 0 {
        RealSign::Plus
    } else {
        RealSign::Minus
    }
}

/// The index of the class an operation denotes, or `None` if it is not a class.
pub(crate) fn find_class(
    set: &MagneticOperationSet,
    operation: &MagneticOperation,
    lattice: &[[Rat; 3]],
) -> Option<usize> {
    let canonical = magnetic_embedding::canonical_translation(operation.translation, lattice);
    set.classes.iter().position(|candidate| {
        candidate.rotation == operation.rotation
            && candidate.time_reversal == operation.time_reversal
            && magnetic_embedding::canonical_translation(candidate.translation, lattice)
                == canonical
    })
}

/// Resolve the source irrep and refuse everything outside the first family.
fn resolve_source(source: SourceIdentity<'_>) -> Result<&'static IrrepRecord, MagneticInputError> {
    let irrep = select_irrep(source.parent_sg, source.label, source.convention)?;
    if irrep.spinor {
        return Err(MagneticInputError::SpinorSource {
            sg: source.parent_sg,
            label: irrep.ml,
        });
    }
    if irrep.dim != 1 {
        return Err(MagneticInputError::NotScalarSource {
            sg: source.parent_sg,
            label: irrep.ml,
            dimension: irrep.dim,
        });
    }
    let k = irrep.k_vector();
    if k.numerators != [0, 0, 0] {
        return Err(MagneticInputError::SourceNotGamma {
            sg: source.parent_sg,
            label: irrep.ml,
            numerators: k.numerators,
            denominator: k.denominator,
        });
    }
    Ok(irrep)
}

/// Resolve the condensing irrep under the record's own convention.
fn resolve_condensing(
    record: RecordIdentity<'_>,
) -> Result<&'static IrrepRecord, MagneticInputError> {
    select_irrep(record.parent_sg, record.condensing_label, record.convention)
}

/// Resolve one label in one convention, refusing a missing or ambiguous one.
fn select_irrep(
    sg: u8,
    label: &str,
    convention: LabelConvention,
) -> Result<&'static IrrepRecord, MagneticInputError> {
    let candidates: Vec<&'static IrrepRecord> = crate::irrep::query::irreps_of(sg)
        .iter()
        .filter(|irrep| irrep.label(convention).as_deref() == Some(label))
        .collect();
    match candidates.len() {
        0 => Err(MagneticInputError::LabelNotFound {
            sg,
            label: label.to_string(),
            convention,
        }),
        1 => Ok(candidates[0]),
        _ => Err(MagneticInputError::AmbiguousLabel {
            sg,
            label: label.to_string(),
        }),
    }
}

/// Exact small-denominator conversion of one `f64` value.
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

// ---------------------------------------------------------------------------
// Parent frame
// ---------------------------------------------------------------------------

/// The parent's operations and translation-lattice generators, exactly.
///
/// Same conventions as the embedding layer: the operation table lists only the
/// non-zero pure translations, so `Z^3` stays implicit in the lattice.
///
/// Crate-visible so the R29 matrix entry point can map magnetic classes into the
/// same parent frame this scalar contract uses instead of re-deriving it.
pub(crate) struct ParentFrame {
    pub(crate) operations: Vec<(Mat3I, [Rat; 3])>,
    pub(crate) lattice: Vec<[Rat; 3]>,
}

impl ParentFrame {
    pub(crate) fn of(sg: u8) -> Result<Self, MagneticInputError> {
        let parent = crate::irrep::query::symmetry_operations_of(sg)
            .map_err(|error| MagneticInputError::DatabaseMismatch(error.to_string()))?;
        let mut operations = Vec::with_capacity(parent.operations.len());
        let mut lattice: Vec<[Rat; 3]> = Vec::new();
        for operation in &parent.operations {
            let mut translation = [Rat::ZERO; 3];
            for (axis, value) in operation.translation.iter().enumerate() {
                translation[axis] =
                    exact_fraction(*value).ok_or(MagneticInputError::ParentTranslationOffGrid {
                        sg,
                        rotation: operation.rotation,
                        value: *value,
                    })?;
            }
            if operation.rotation == IDENTITY_ROTATION
                && translation.iter().any(|value| !value.is_zero())
                && !lattice.contains(&translation)
            {
                lattice.push(translation);
            }
            operations.push((operation.rotation, translation));
        }
        Ok(Self {
            operations,
            lattice,
        })
    }

    /// Whether some parent operation has this rotation and lands in this
    /// translation's class modulo the parent's lattice.
    pub(crate) fn has_operation(&self, rotation: Mat3I, translation: [Rat; 3]) -> bool {
        self.operations.iter().any(|(candidate, parent)| {
            *candidate == rotation && differs_by_lattice(&translation, parent, &self.lattice)
        })
    }
}

/// Whether two exact points differ by a lattice vector.
pub(crate) fn differs_by_lattice(left: &[Rat; 3], right: &[Rat; 3], lattice: &[[Rat; 3]]) -> bool {
    let mut difference = [Rat::ZERO; 3];
    for axis in 0..3 {
        let Ok(value) = left[axis].checked_sub(right[axis]) else {
            return false;
        };
        difference[axis] = value;
    }
    magnetic_embedding::canonical_translation(difference, lattice)
        .iter()
        .all(|value| value.is_zero())
}

// ---------------------------------------------------------------------------
// Mapped classes
// ---------------------------------------------------------------------------

/// One class of the magnetic group, mapped into the parent's frame.
pub(crate) struct MappedClass {
    /// The `(rotation, translation)` the class lands on.
    pub(crate) operation: (Mat3I, [Rat; 3]),
    /// The time-reversal flag, which a spatial map does not change.
    pub(crate) time_reversal: bool,
}

/// Map every class into the parent with the caller's `(map, shift)`.
///
/// The translation rule is the one
/// [`crate::irrep::magnetic_embedding::verify_embedding`] checks:
/// `A τ + (I − R_p) δ` with `R_p = A R A^{-1}`.  A class that hits no parent
/// operation is reported: for a verified embedding it cannot happen, and if it
/// does the caller must not be told a plausible table.
pub(crate) fn mapped_classes(
    set: &MagneticOperationSet,
    evidence: &EmbeddingEvidence,
    parent: &ParentFrame,
) -> Result<Vec<MappedClass>, MagneticInputError> {
    let inverse = evidence
        .map
        .inverse()
        .map_err(|_| MagneticInputError::Arithmetic("the selected map is singular".to_string()))?;
    let mut mapped = Vec::with_capacity(set.classes.len());
    for class in &set.classes {
        let rotation = Mat3R::from_ints(class.rotation);
        let conjugated = evidence
            .map
            .checked_mul(&rotation)
            .and_then(|product| product.checked_mul(&inverse))
            .map_err(|_| {
                MagneticInputError::Arithmetic("conjugating a class rotation".to_string())
            })?;
        let mapped_rotation = conjugated.to_int_matrix().map_err(|_| {
            MagneticInputError::Arithmetic("the mapped rotation is not integral".to_string())
        })?;
        let moved = evidence
            .map
            .checked_mul_vector(&Vec3R::new(class.translation))
            .map_err(|_| {
                MagneticInputError::Arithmetic("mapping a class translation".to_string())
            })?;
        let mut translation = [Rat::ZERO; 3];
        for row in 0..3 {
            let mut sum = Rat::ZERO;
            for (column, shift) in evidence.shift.iter().enumerate() {
                let entry = i128::from(if row == column { 1 } else { 0 })
                    - i128::from(mapped_rotation[row][column]);
                let term = shift.checked_mul(Rat::from_integer(entry)).map_err(|_| {
                    MagneticInputError::Arithmetic("applying the origin shift".to_string())
                })?;
                sum = sum.checked_add(term).map_err(|_| {
                    MagneticInputError::Arithmetic("applying the origin shift".to_string())
                })?;
            }
            translation[row] = moved.get(row).checked_add(sum).map_err(|_| {
                MagneticInputError::Arithmetic("placing a class on a parent operation".to_string())
            })?;
        }
        if !parent.has_operation(mapped_rotation, translation) {
            return Err(MagneticInputError::Arithmetic(format!(
                "the mapped class with rotation {mapped_rotation:?} hits no parent operation"
            )));
        }
        mapped.push(MappedClass {
            operation: (mapped_rotation, translation),
            time_reversal: class.time_reversal,
        });
    }
    Ok(mapped)
}

// ---------------------------------------------------------------------------
// The parent's real +-1 character
// ---------------------------------------------------------------------------

/// The character the parent's tabulated row gives to every mapped class.
///
/// The row is the exact `f64` storage snapped to rationals.  A value that is not
/// real `±1` is refused with its own error; the identity class's value must be
/// `+1`, which is what makes the product checks meaningful.
fn read_character(
    sg: u8,
    row: &CharacterRow,
    set: &MagneticOperationSet,
    mapped: &[MappedClass],
    parent: &ParentFrame,
) -> Result<Vec<RealSign>, MagneticInputError> {
    let negative_one = Rat::ONE.checked_neg().expect("one is representable");
    let mut out = Vec::with_capacity(mapped.len());
    for (index, class) in mapped.iter().enumerate() {
        let (rotation, translation) = &class.operation;
        let mut entry = None;
        for (operation, value) in row.operations().iter().zip(row.values()) {
            if flatten_rotation(&operation.rotation) != *rotation {
                continue;
            }
            match row_translation_differs_by_lattice(operation, translation, &parent.lattice) {
                Some(true) => entry = Some(*value),
                Some(false) => continue,
                // A row translation off the exact grid makes the row unusable
                // here, and reporting "no entry" would be a different (and
                // false) statement, so it has its own error.
                None => {
                    let offending = operation
                        .translation
                        .iter()
                        .find(|value| exact_fraction(**value).is_none())
                        .copied()
                        .unwrap_or(0.0);
                    return Err(MagneticInputError::ParentTranslationOffGrid {
                        sg,
                        rotation: *rotation,
                        value: offending,
                    });
                }
            }
        }
        let value = entry.ok_or(MagneticInputError::CharacterRowSettingMismatch {
            sg,
            rotation: *rotation,
            uni: set.uni,
        })?;
        let real = exact_fraction(value.re).ok_or(MagneticInputError::OffGrid {
            sg,
            rotation: *rotation,
            value: value.re,
        })?;
        let imaginary = exact_fraction(value.im).ok_or(MagneticInputError::OffGrid {
            sg,
            rotation: *rotation,
            value: value.im,
        })?;
        if !imaginary.is_zero() || (real != Rat::ONE && real != negative_one) {
            return Err(MagneticInputError::NonRealCharacter {
                sg,
                rotation: *rotation,
                real: value.re,
                imaginary: value.im,
            });
        }
        let is_identity_class = !set.classes[index].time_reversal
            && *rotation == IDENTITY_ROTATION
            && translation.iter().all(|value| value.is_zero());
        if is_identity_class && real != Rat::ONE {
            return Err(MagneticInputError::IdentityCharacterNotOne {
                sg,
                value: value.re,
            });
        }
        out.push(if real == Rat::ONE {
            RealSign::Plus
        } else {
            RealSign::Minus
        });
    }
    Ok(out)
}

/// The character row stores rotations row-major in a flat array of nine.
fn flatten_rotation(rotation: &[i32; 9]) -> Mat3I {
    [
        [rotation[0], rotation[1], rotation[2]],
        [rotation[3], rotation[4], rotation[5]],
        [rotation[6], rotation[7], rotation[8]],
    ]
}

/// Whether a row operation's translation lands in the mapped translation's
/// class; `None` when the row translation is not on the exact grid.
fn row_translation_differs_by_lattice(
    operation: &SeitzOperation,
    translation: &[Rat; 3],
    lattice: &[[Rat; 3]],
) -> Option<bool> {
    let mut row_translation = [Rat::ZERO; 3];
    for (axis, value) in operation.translation.iter().enumerate() {
        row_translation[axis] = exact_fraction(*value)?;
    }
    Some(differs_by_lattice(translation, &row_translation, lattice))
}

// ---------------------------------------------------------------------------
// The corepresentation law
// ---------------------------------------------------------------------------

/// Check `A(g)·[conj(A(h)) if t(g) else A(h)] = A(g·h)` on every ordered pair of
/// classes.
///
/// The bracket is decided by the time-reversal flag of the **left** factor; the
/// unitary branch keeps `A(h)` unchanged.  Every value is `±1` and
/// `conj(±1) = ±1`, so the conjugation is the identity on the values while the
/// branch is written structurally.
///
/// What this proves is exactly one thing: a phase table that violates the law is
/// refused, with the failing pair.  It does not claim to detect every wrong
/// embedding or every wrong `η` — both signs of `η` satisfy the law by
/// construction, so a table that closes is a valid corepresentation of the grey
/// parent and not a statement that the caller's other choices were right.
fn check_corepresentation_law(
    set: &MagneticOperationSet,
    actions: &[ClassAction],
    lattice: &[[Rat; 3]],
) -> Result<(), MagneticInputError> {
    for left in 0..set.classes.len() {
        for right in 0..set.classes.len() {
            let product = magnetic_embedding::compose(&set.classes[left], &set.classes[right])
                .map_err(|error| MagneticInputError::Arithmetic(error.to_string()))?;
            let product_index = find_class(set, &product, lattice).ok_or_else(|| {
                MagneticInputError::Arithmetic(format!(
                    "the product of classes {left} and {right} is not a class"
                ))
            })?;
            // `A(g)·conj^{t(g)}(A(h))`: the conjugation slot is decided by the
            // left factor's flag, and for the `±1` phases used here
            // `conj(z) = z`, so a single arithmetic path serves both branches.
            let _conjugation_slot = set.classes[left].time_reversal;
            let expected = actions[left].value.value() * actions[right].value.value();
            let found = actions[product_index].value.value();
            if expected != found {
                return Err(MagneticInputError::NotACorepresentation {
                    uni: set.uni,
                    left,
                    right,
                    expected,
                    found,
                });
            }
        }
    }
    Ok(())
}
