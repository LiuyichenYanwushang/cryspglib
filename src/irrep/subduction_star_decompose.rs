//! Staged full-star decomposition of a scalar parent irrep (tasks 8b/8c of
//! `docs/full-irrep-subduction-plan.md`).
//!
//! [`subduce_full_star_with_embedding`] takes the full parent star built by
//! [`ScalarStar`], folds its arms into the embedding with
//! [`ScalarStar::folded_stars`], and decomposes every folded child star:
//!
//! * the child table stores one representative **arm** per irrep, so the
//!   representative `q` of a child star is the folded point whose class matches
//!   a stored child `k` exactly modulo the **child** reciprocal lattice
//!   (centring extinctions included), searched over the whole folded orbit —
//!   never `points()[0]` by assumption and never a nearest match.  A
//!   `ConjugateRealification` row is reachable through its conjugate arm `-k`
//!   too, even when no `IrrepRecord` stores that arm;
//! * the child target catalogue enumerates **complex components**, not physical
//!   rows: an ordinary row is one component, a `DistinctComponentSum` row keeps
//!   both stored CIR constituents, and a `ConjugateRealification` row keeps its
//!   seed at `k` and the conjugate at `-k`.  A component whose own arm differs
//!   from the representative `q` inside the same child star is transported
//!   exactly by a child group element before its row is compared;
//! * the seed and conjugate of one realification are compared once, narrowly
//!   inside that record: equal unit rows are **one** complex target
//!   (canonicalized to the seed, so the solver accumulates the multiplicity),
//!   orthogonal rows stay two targets, and anything in between is an explicit
//!   error.  Unrelated CIR identities are never merged;
//! * on `H_q`, the little group of `q`, the q-block character is the induced
//!   one restricted to the parent arms that fold onto `q`:
//!   `chi_q(h) = sum_i chi_component(g_i^-1 h g_i)`, with every Seitz
//!   translation kept for the row lookup and every arm contributing its own
//!   component dimension.  Fixing `q` can still permute parent arms that fold
//!   onto it (e.g. SG 139 `X1+` restricted to #126); moved arms contribute
//!   zero. The full-star trace is never divided by an arm count;
//! * the q-block is decomposed against the prepared complex child rows by the
//!   shared `solve_prepared_character_block` solver, so the Gram identity, the
//!   non-negative integral multiplicities, the little-dimension sum and the
//!   per-operation reconstruction are the same checks the single arm entry
//!   point uses;
//! * the whole subduced full-star character on **all** subgroup coset
//!   representatives is then rebuilt from the reported child irreps induced
//!   over their own **child** stars and compared with the parent full-star
//!   character.  That reconstruction never calls the q-block builder again, so
//!   a wrong frame transport between parent arms and child stars cannot cancel
//!   out.
//!
//! Boundaries: scalar parents (ordinary and canonical compound rows) only;
//! magnetic groups and spinor probes stay out of this stage.  Missing child
//! `k`/irrep data is a typed error: a partial set of blocks is never returned.

use num_complex::Complex64;

use crate::irrep::isotropy::IsotropySubgroup;
use crate::irrep::query;
use crate::irrep::types::{
    CharacterRow, CompoundCharacterSemantics, CompoundSelectedArmCharacter, IrrepRecord,
    IrrepSourceIdentity,
};
use crate::mathfunc::Mat3I;

use super::super::{
    ComplexTarget, ExactSeitz, Lattice, Rat, SUBDUCTION_TOLERANCE, SubductionComponent,
    SubductionError, SubductionTarget, SubgroupEmbedding, Vec3R, character_of,
    exact_primitive_basis, inline_k_vector, shift_operations, solve_prepared_character_block,
    strict_sg_hall_ops, validate_subduction_context,
};
use super::scalar_star::{ComponentStar, ScalarStar};
use super::{FoldedStar, OrdinaryStar, StarError, arm_wave_vector};

/// The identity rotation, as stored in every Hall operation table.
const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

// ── Errors ───────────────────────────────────────────────────────────────────

/// Errors from the staged full-star decomposition.
///
/// The exact subduction layer and the star adapter keep their own error types
/// ([`FullStarError::Subduction`], [`FullStarError::Star`]); the variants below
/// are specific to the block and star bookkeeping this stage adds.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FullStarError {
    /// An error from the exact rational/affine layer.
    #[error(transparent)]
    Subduction(#[from] SubductionError),
    /// An error from the full-star / folded-star adapter.
    #[error(transparent)]
    Star(#[from] StarError),
    /// No folded point of one child star matches stored child irrep data.
    ///
    /// The child table only lists a representative arm, so every point of the
    /// orbit is tried before this is reported; a partial result is never
    /// produced.
    #[error(
        "no stored irrep of space group {sg} matches any of the {points} folded points of a \
         child star (first q = ({}, {}, {}))",
        q[0],
        q[1],
        q[2]
    )]
    MissingChildStarData {
        /// Child space group number.
        sg: u8,
        /// First folded point of the unmatched star, in the child frame.
        q: [Rat; 3],
        /// Number of folded points in the star.
        points: usize,
    },
    /// No child operation transports a component's effective arm to the block
    /// representative inside the child star.
    ///
    /// Reachability is per component: an ordinary or `DistinctComponentSum`
    /// component sits at its stored `k`, a `ConjugateRealification` component
    /// at its seed `k` or at `-k`.  The variant exists so an incoherent folded
    /// geometry fails closed instead of comparing a row at the wrong arm.
    #[error(
        "no operation of space group {sg} transports the effective arm ({}, {}, {}) to the \
         block representative ({}, {}, {})",
        k.get(0),
        k.get(1),
        k.get(2),
        q.get(0),
        q.get(1),
        q.get(2)
    )]
    MissingComponentTransport {
        /// Child space group number.
        sg: u8,
        /// Component effective arm, in the child frame.
        k: Box<Vec3R>,
        /// Block representative, in the child frame.
        q: Box<Vec3R>,
    },
    /// A realification seed or conjugate row is not a unit complex irrep over
    /// `H_q`, so the equivalence test is meaningless.
    #[error("realification {ml} of space group {sg} has Gram norm {norm} on H_q instead of 1")]
    RealificationNormMismatch {
        /// Child space group number.
        sg: u8,
        /// Realification row label.
        ml: &'static str,
        /// Computed norm.
        norm: f64,
    },
    /// A realification's transported seed and conjugate rows are neither equal
    /// (one complex target) nor orthogonal (two complex targets).
    #[error(
        "realification {ml} of space group {sg}: transported seed/conjugate rows are neither \
         orthogonal nor pointwise equal (overlap {inner})"
    )]
    RealificationOverlap {
        /// Child space group number.
        sg: u8,
        /// Realification row label.
        ml: &'static str,
        /// Inner product over `H_q`.
        inner: Complex64,
    },
    /// The selected little group has no identity operation modulo the child
    /// lattice, so `chi_q(E)` cannot be read at the true identity.
    #[error("the selected little group has no identity operation modulo the child lattice")]
    MissingIdentityOperation,
    /// A folded point reached the block builder without a parent arm.
    ///
    /// [`FoldedStar`] never produces such a point; the variant exists so a
    /// mutated or incoherent geometry fails closed instead of asking for a
    /// zero-arm character.
    #[error("child star q = ({}, {}, {}) has no parent arms", q[0], q[1], q[2])]
    EmptyQBlock {
        /// Representative folded point.
        q: [Rat; 3],
    },
    /// `chi_q(E)` does not reproduce the q-block geometry.
    #[error(
        "child star q = ({}, {}, {}): chi_q(E) = {identity} does not equal the q-block \
         dimension {expected}",
        q[0],
        q[1],
        q[2]
    )]
    QBlockIdentityMismatch {
        /// Representative folded point.
        q: [Rat; 3],
        /// Character evaluated at the identity operation.
        identity: Complex64,
        /// `arms_at_q × selected-arm dimension`.
        expected: u32,
    },
    /// The shared solver's dimension disagrees with the q-block geometry.
    #[error(
        "child star q = ({}, {}, {}): solved dimension {solved} does not equal the q-block \
         dimension {expected}",
        q[0],
        q[1],
        q[2]
    )]
    QBlockDimensionMismatch {
        /// Representative folded point.
        q: [Rat; 3],
        /// Dimension returned by the block solver.
        solved: u32,
        /// `arms_at_q × selected-arm dimension`.
        expected: u32,
    },
    /// `Σ multiplicity × little dimension × star size` does not equal the
    /// parent subspace the child star carries.
    #[error(
        "child star q = ({}, {}, {}): {little_dimension} little dimension x star size \
         {star_size} = {carried} does not equal the folded block dimension {block_dimension}",
        q[0],
        q[1],
        q[2]
    )]
    StarDimensionMismatch {
        /// Representative folded point.
        q: [Rat; 3],
        /// `Σ multiplicity × child little dimension`.
        little_dimension: u32,
        /// Number of folded points in the star.
        star_size: usize,
        /// `little_dimension × star_size`.
        carried: u32,
        /// `FoldedStar::block_dimension()`.
        block_dimension: u32,
    },
    /// The child stars do not tile the parent full-star dimension.
    #[error("child stars carry {found} of the parent full-star dimension {expected}")]
    TotalDimensionMismatch {
        /// `ScalarStar::dimension()`.
        expected: u32,
        /// Sum of the reported block dimensions.
        found: u32,
    },
    /// No unique child record backed a reported target.
    #[error("child target {ml} of space group {sg} has no unique source record in its q block")]
    AmbiguousTargetSource {
        /// Child space group number.
        sg: u8,
        /// Target label.
        ml: &'static str,
    },
    /// A reported constituent does not agree with its stored CIR source.
    #[error("child target {ml} of space group {sg} disagrees with its stored CIR source")]
    TargetSourceMismatch {
        /// Child space group number.
        sg: u8,
        /// Target label.
        ml: &'static str,
    },
    /// Rebuilding the subduced character from the reported child irreps differs
    /// from the parent full-star character.
    #[error(
        "restricted full-star reconstruction failed on subgroup representative {index}: \
         {found} != {expected}"
    )]
    ReconstructionMismatch {
        /// Position in the embedding's representative list.
        index: usize,
        /// Character rebuilt from the reported child irreps.
        found: Complex64,
        /// `ScalarStar::character` of the same representative.
        expected: Complex64,
    },
}

// ── Result types ─────────────────────────────────────────────────────────────

/// One non-zero child irrep term of a folded child star.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FullStarTarget {
    /// Child space group number.
    pub sg: u8,
    /// Target label: the ordinary row's ML label, or the CIR constituent label.
    /// A realification seed and its conjugate share the one CIR label.
    pub ml: &'static str,
    /// Bradley-Cracknell label, for display only.
    pub bc: &'static str,
    /// Physical child row this term came from.
    pub row_ml: &'static str,
    /// Which complex constituent of that row this term is.  The variant
    /// distinguishes a realification seed from its conjugate even though both
    /// carry the same label and source number.
    pub component: SubductionComponent,
    /// Complex dimension of the child little-group irrep.
    pub dimension: u8,
    /// Multiplicity in the child little-group representation at `q`.
    pub multiplicity: u32,
    /// Frozen CIR source number of this term (`record.source_identity()` for an
    /// ordinary row, the constituent's own `irnumber` for a compound row, the
    /// shared seed `irnumber` for both realification components).
    pub irnumber: u32,
}

/// The decomposition of one folded child star.
#[derive(Debug, Clone)]
pub struct FullStarBlock {
    q: Vec3R,
    stored_k: Vec3R,
    star_size: usize,
    arm_count: usize,
    block_dimension: u32,
    little_dimension: u32,
    targets: Vec<FullStarTarget>,
    /// One child full-star character evaluator per target, used by the
    /// independent reconstruction; not part of the reported decomposition.
    evaluators: Vec<ChildStarEvaluator>,
}

impl FullStarBlock {
    /// The representative folded point whose class matched stored child data.
    pub const fn q(&self) -> &Vec3R {
        &self.q
    }

    /// The **effective** component wave vector the representative `q` matched
    /// modulo the child reciprocal lattice.
    ///
    /// For an ordinary or `DistinctComponentSum` block this is the stored child
    /// `k` of the child table's representative arm.  For a
    /// `ConjugateRealification` block reached through the conjugate arm it is
    /// the exact negation of the record's stored `k`, even though no
    /// `IrrepRecord` stores that arm; it is never reduced or wrapped to a
    /// positive representative.
    pub const fn stored_k(&self) -> &Vec3R {
        &self.stored_k
    }

    /// Number of folded points in this child star.
    pub const fn star_size(&self) -> usize {
        self.star_size
    }

    /// Number of parent component arms in this child star. Coincident arms of
    /// different components are counted separately.
    pub const fn arm_count(&self) -> usize {
        self.arm_count
    }

    /// Parent full-star subspace carried by this block:
    /// `arm_count × complex-component selected-arm dimension`.
    pub const fn block_dimension(&self) -> u32 {
        self.block_dimension
    }

    /// `Σ multiplicity × child little dimension` on `H_q`.
    pub const fn little_dimension(&self) -> u32 {
        self.little_dimension
    }

    /// Non-zero target terms, in the child table's order.
    pub fn targets(&self) -> &[FullStarTarget] {
        &self.targets
    }

    /// Multiplicity of one target label, or `0` when it does not appear.
    pub fn multiplicity(&self, ml: &str) -> u32 {
        self.targets
            .iter()
            .find(|target| target.ml == ml)
            .map_or(0, |target| target.multiplicity)
    }
}

/// The complete decomposition of one scalar parent full star on one embedded
/// subgroup.
#[derive(Debug, Clone)]
pub struct FullStarSubduction {
    parent_sg: u8,
    parent_ml: &'static str,
    parent_bc: &'static str,
    parent_source_identity: IrrepSourceIdentity,
    parent_irnumber: Option<u32>,
    parent_dimension: u32,
    subgroup_sg: u8,
    ordinal: usize,
    setting: Mat3I,
    seed_k: Vec3R,
    blocks: Vec<FullStarBlock>,
    representatives: Vec<ExactSeitz>,
    parent_characters: Vec<Complex64>,
    reconstructed: Vec<Complex64>,
    tolerance: f64,
}

impl FullStarSubduction {
    /// Parent space group number.
    pub const fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// Miller-Love label of the subduced parent irrep.
    pub const fn parent_ml(&self) -> &'static str {
        self.parent_ml
    }

    /// Bradley-Cracknell label of the subduced parent irrep.
    pub const fn parent_bc(&self) -> &'static str {
        self.parent_bc
    }

    /// Frozen source identity of the parent probe: `OrdinaryScalar` for an
    /// ordinary row, `Compound` for a compound row.
    pub const fn parent_source_identity(&self) -> IrrepSourceIdentity {
        self.parent_source_identity
    }

    /// Frozen CIR source number of the parent probe, when it has exactly one:
    /// `Some` for an ordinary row, `None` for a compound row (which is a sum of
    /// two CIR sources and therefore has no single source number).
    pub const fn parent_irnumber(&self) -> Option<u32> {
        self.parent_irnumber
    }

    /// Full-star dimension of the parent probe (`OrdinaryStar::dimension()`).
    pub const fn parent_dimension(&self) -> u32 {
        self.parent_dimension
    }

    /// Subgroup space group number.
    pub const fn subgroup_sg(&self) -> u8 {
        self.subgroup_sg
    }

    /// Isotropy record ordinal of the embedding.
    pub const fn ordinal(&self) -> usize {
        self.ordinal
    }

    /// The setting transform of the embedding this result belongs to.
    pub const fn setting(&self) -> Mat3I {
        self.setting
    }

    /// The parent probe's stored selected-arm wave vector.
    ///
    /// This is the seed component's `base_k`; a conjugated component's
    /// **effective** arm is its exact negation, and each block reports the
    /// effective arm it was decomposed at through [`FullStarBlock::stored_k`].
    pub const fn seed_k(&self) -> &Vec3R {
        &self.seed_k
    }

    /// Folded child-star blocks, in canonical `q` order.
    pub fn blocks(&self) -> &[FullStarBlock] {
        &self.blocks
    }

    /// Subgroup coset representatives the reconstruction was checked on.
    pub fn representatives(&self) -> &[ExactSeitz] {
        &self.representatives
    }

    /// The parent full-star character on the subgroup representatives, and its
    /// reconstruction by inducing the reported child irreps over their own
    /// child stars: equal entries are the witness that the decomposition is
    /// complete in the parent frame, not only per `q`.
    pub fn reconstruction(&self) -> (&[Complex64], &[Complex64]) {
        (&self.parent_characters, &self.reconstructed)
    }

    /// Tolerance used for the multiplicity, dimension and reconstruction
    /// checks.
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Total parent subspace accounted for by the reported blocks; equal to
    /// [`Self::parent_dimension`] by construction.
    pub fn covered_dimension(&self) -> u32 {
        self.blocks.iter().map(FullStarBlock::block_dimension).sum()
    }
}

// ── Child full-star character evaluators (independent reconstruction) ────────

/// The child full-star character of one reported target.
#[derive(Debug, Clone)]
enum ChildStarEvaluator {
    /// An ordinary child row: the record's own full star.
    Ordinary(OrdinaryStar),
    /// One complex component of a compound child row (a
    /// `DistinctComponentSum` constituent, or a realification seed/conjugate
    /// with its effective `k` and conjugation flag).
    Component(ComponentStar),
}

impl ChildStarEvaluator {
    fn character(&self, operation: &ExactSeitz) -> Result<Complex64, FullStarError> {
        match self {
            Self::Ordinary(star) => Ok(star.character(operation)?),
            Self::Component(star) => Ok(star.character(operation)?),
        }
    }
}

// ── Entry point ──────────────────────────────────────────────────────────────

/// Decompose the full star of a scalar parent irrep on the subgroup selected by
/// an isotropy record and a prebuilt embedding.
///
/// `subgroup`, `embedding` and `probe` are revalidated as one context before any
/// block is built (see `validate_subduction_context`).  Ordinary and canonical
/// compound scalar probes are expanded by [`ScalarStar`]; spinor probes are
/// rejected with [`SubductionError::UnsupportedCharacterSpace`], and a folded
/// child star without stored child data is a
/// [`FullStarError::MissingChildStarData`] error rather than a partial result.
pub fn subduce_full_star_with_embedding(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<FullStarSubduction, FullStarError> {
    validate_subduction_context(subgroup, embedding, probe)?;
    if probe.spinor {
        return Err(SubductionError::UnsupportedCharacterSpace {
            sg: probe.sg,
            ml: probe.ml.to_string(),
        }
        .into());
    }
    let star = ScalarStar::new(probe)?;
    decompose_scalar_star(embedding, &star)
}

/// The reusable core: decompose an already validated full parent star.
///
/// Split out so the tests can drive a differently ordered but equivalent
/// transversal through [`ScalarStar::from_transporters`] and a mutated folded
/// geometry through [`decompose_folded_stars`]; both are private, so no caller
/// can reach them with an unvalidated context.
fn decompose_scalar_star(
    embedding: &SubgroupEmbedding,
    star: &ScalarStar,
) -> Result<FullStarSubduction, FullStarError> {
    let folded = star.folded_stars(embedding)?;
    decompose_folded_stars(embedding, star, &folded)
}

/// Decompose one folded geometry of an already validated star.
fn decompose_folded_stars(
    embedding: &SubgroupEmbedding,
    star: &ScalarStar,
    folded: &[FoldedStar],
) -> Result<FullStarSubduction, FullStarError> {
    let child_sg = embedding.subgroup_sg();
    let child_cell = Lattice::new(exact_primitive_basis(child_sg)?)?;
    let child_reciprocal = child_cell.reciprocal()?;

    let mut blocks = Vec::with_capacity(folded.len());
    let mut covered = 0u32;
    for folded_star in folded {
        let block = build_block(embedding, star, folded_star, &child_cell, &child_reciprocal)?;
        covered = covered.checked_add(block.block_dimension()).ok_or(
            SubductionError::RationalOverflow {
                operation: "full-star dimension",
            },
        )?;
        blocks.push(block);
    }
    if covered != star.dimension() {
        return Err(FullStarError::TotalDimensionMismatch {
            expected: star.dimension(),
            found: covered,
        });
    }

    let (parent_characters, reconstructed) = reconstruct(embedding, star, &blocks)?;
    Ok(FullStarSubduction {
        parent_sg: star.parent_sg(),
        parent_ml: star.probe_ml(),
        parent_bc: star.probe().bc,
        parent_source_identity: star.source_identity(),
        parent_irnumber: parent_irnumber(star.probe()),
        parent_dimension: star.dimension(),
        subgroup_sg: child_sg,
        ordinal: embedding.ordinal(),
        setting: embedding.setting(),
        seed_k: *star.seed_k(),
        blocks,
        representatives: embedding.representatives().to_vec(),
        parent_characters,
        reconstructed,
        tolerance: SUBDUCTION_TOLERANCE,
    })
}

/// The single CIR source number of a parent probe, when it has one.
fn parent_irnumber(record: &'static IrrepRecord) -> Option<u32> {
    match record.source_identity() {
        IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => Some(cir_irnumber),
        IrrepSourceIdentity::Compound { .. } | IrrepSourceIdentity::Spin { .. } => None,
    }
}

// ── One folded child star ────────────────────────────────────────────────────

/// One complex component of a child record, evaluated at its own stored row.
#[derive(Debug, Clone)]
struct ChildComponent {
    record: &'static IrrepRecord,
    component: SubductionComponent,
    label: &'static str,
    bc: &'static str,
    irnumber: u32,
    dimension: u8,
    base_k: Vec3R,
    /// `base_k`, or its exact negation for a realification conjugate.
    effective_k: Vec3R,
    conjugate: bool,
    row: CharacterRow,
}

impl ChildComponent {
    /// Stable identity used to keep one occurrence of one component per orbit.
    fn key(&self) -> (u16, SubductionComponent) {
        (self.record.id().index(), self.component)
    }
}

/// The representative point of one folded child star and every complex child
/// component whose effective arm lies in that star.
#[derive(Debug)]
struct Representative {
    /// Position in the folded star's point list.
    point: usize,
    /// The representative folded wave vector (child frame).
    q: Vec3R,
    /// The effective component arm that matched this representative.
    effective_k: Vec3R,
    /// Every reachable complex component, in the child table's order.
    components: Vec<ChildComponent>,
}

/// Decompose one folded child star: pick its representative arm, build the
/// `H_q` q-block character and solve it against the prepared complex child
/// components.
fn build_block(
    embedding: &SubgroupEmbedding,
    star: &ScalarStar,
    folded_star: &FoldedStar,
    child_cell: &Lattice,
    child_reciprocal: &Lattice,
) -> Result<FullStarBlock, FullStarError> {
    let child_sg = embedding.subgroup_sg();
    let representative = select_representative(child_sg, folded_star, child_reciprocal)?;
    let point = &folded_star.points()[representative.point];
    let q = representative.q;
    let q_key = [q.get(0), q.get(1), q.get(2)];

    let (parent_operations, pulled_back) =
        little_group_operations(embedding, &q, child_reciprocal)?;
    let identity = identity_position(&pulled_back, child_cell)?;

    // The q-block character: exactly the parent arms folding onto this q, each
    // contributing its own component's row and selected dimension.  No trace is
    // divided by an arm count.
    if point.arm_indices().is_empty() {
        return Err(FullStarError::EmptyQBlock { q: q_key });
    }
    let q_block_dimension = star.q_block_dimension(point.arm_indices())?;
    let mut parent_characters = Vec::with_capacity(parent_operations.len());
    for operation in &parent_operations {
        parent_characters.push(star.q_block_character(point.arm_indices(), operation)?);
    }

    // `chi_q(E)` has to be exactly the q-block geometry, read at the actual
    // identity operation of the aligned list.
    let identity_value = parent_characters[identity];
    if (identity_value - Complex64::new(f64::from(q_block_dimension), 0.0)).norm()
        > SUBDUCTION_TOLERANCE
    {
        return Err(FullStarError::QBlockIdentityMismatch {
            q: q_key,
            identity: identity_value,
            expected: q_block_dimension,
        });
    }

    let prepared = prepare_targets(
        child_sg,
        &representative,
        child_cell,
        child_reciprocal,
        &pulled_back,
    )?;
    let solved =
        solve_prepared_character_block(child_sg, prepared, &parent_characters, &pulled_back)?;
    if u32::from(solved.dimension) != q_block_dimension {
        return Err(FullStarError::QBlockDimensionMismatch {
            q: q_key,
            solved: u32::from(solved.dimension),
            expected: q_block_dimension,
        });
    }

    let mut targets = Vec::with_capacity(solved.targets.len());
    let mut evaluators = Vec::with_capacity(solved.targets.len());
    let mut little_dimension = 0u32;
    for target in &solved.targets {
        let component = component_for_target(child_sg, &representative.components, target)?;
        let evaluator = build_evaluator(child_sg, component)?;
        little_dimension = little_dimension
            .checked_add(
                u32::from(target.dimension)
                    .checked_mul(target.multiplicity)
                    .ok_or(SubductionError::RationalOverflow {
                        operation: "child little-group dimension",
                    })?,
            )
            .ok_or(SubductionError::RationalOverflow {
                operation: "child little-group dimension",
            })?;
        targets.push(FullStarTarget {
            sg: child_sg,
            ml: target.ml,
            bc: target.bc,
            row_ml: target.row_ml,
            component: target.component,
            dimension: target.dimension,
            multiplicity: target.multiplicity,
            irnumber: component.irnumber,
        });
        evaluators.push(evaluator);
    }

    let star_size = folded_star.star_size();
    let carried = u32::try_from(star_size)
        .ok()
        .and_then(|size| little_dimension.checked_mul(size))
        .ok_or(SubductionError::RationalOverflow {
            operation: "child star dimension",
        })?;
    if carried != folded_star.block_dimension() {
        return Err(FullStarError::StarDimensionMismatch {
            q: q_key,
            little_dimension,
            star_size,
            carried,
            block_dimension: folded_star.block_dimension(),
        });
    }
    Ok(FullStarBlock {
        q,
        stored_k: representative.effective_k,
        star_size,
        arm_count: folded_star.arm_count(),
        block_dimension: folded_star.block_dimension(),
        little_dimension,
        targets,
        evaluators,
    })
}

// ── The complex child target catalogue ───────────────────────────────────────

/// `H_q`: the subgroup coset representatives whose child-frame rotation fixes
/// `q` modulo the child reciprocal lattice, paired with their parent-frame
/// images.
///
/// Lattice translations stay in the operations; only the `child_shift` of the
/// embedding is undone, to land in the child Hall frame the shipped rows live
/// in.  The returned lists are aligned: `parent_operations[i]` is the parent
/// representative of `pulled_back[i]`.
fn little_group_operations(
    embedding: &SubgroupEmbedding,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<(Vec<ExactSeitz>, Vec<ExactSeitz>), FullStarError> {
    let mut parent_operations: Vec<ExactSeitz> = Vec::new();
    let mut child_operations: Vec<ExactSeitz> = Vec::new();
    for operation in embedding.representatives() {
        let child = embedding.transform().unmap_operation(operation)?;
        if !child_reciprocal.preserves(child.rotation(), q)? {
            continue;
        }
        parent_operations.push(*operation);
        child_operations.push(child);
    }
    if parent_operations.is_empty() {
        return Err(FullStarError::MissingIdentityOperation);
    }
    let pulled_back = shift_operations(&child_operations, &embedding.child_shift().checked_neg()?)?;
    Ok((parent_operations, pulled_back))
}

/// Reachability: which child records have a complex component whose **effective**
/// arm folds onto `q`.
///
/// The child table stores one representative arm per irrep, so an ordinary or
/// `DistinctComponentSum` row is reachable at its stored `k`, and a
/// `ConjugateRealification` row is reachable at its stored `k` **and** at `-k`,
/// even when no `IrrepRecord` stores the negated arm.  Spinor rows are a
/// different representation space and are not targets of a scalar parent.
fn child_components_at(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<Vec<ChildComponent>, FullStarError> {
    let mut out = Vec::new();
    for record in query::irreps_of(child_sg) {
        for component in child_components(record)? {
            if child_reciprocal.same_mod(q, &component.effective_k)? {
                out.push(component);
            }
        }
    }
    Ok(out)
}

/// Expand one non-spinor child record into its complex components.
fn child_components(record: &'static IrrepRecord) -> Result<Vec<ChildComponent>, FullStarError> {
    if record.spinor {
        return Ok(Vec::new());
    }
    let base_k = inline_k_vector(record)?;
    let unsupported = || SubductionError::UnsupportedCharacterSpace {
        sg: record.sg,
        ml: record.ml.to_string(),
    };
    match record.compound_metadata() {
        None => {
            let row = record
                .ordinary_scalar_selected_arm_block_trace()
                .map_err(|_| unsupported())?;
            let dimension = component_dimension(row.dimension(), record.sg, record.ml)?;
            let irnumber = ordinary_irnumber(record)?;
            Ok(vec![ChildComponent {
                record,
                component: SubductionComponent::Ordinary,
                label: record.ml,
                bc: record.bc,
                irnumber,
                dimension,
                base_k,
                effective_k: base_k,
                conjugate: false,
                row,
            }])
        }
        Some(_) => {
            let view = record
                .compound_selected_arm_view()
                .map_err(|_| unsupported())?;
            match view {
                CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } => {
                    Ok(vec![
                        constituent_component(record, 0, first, base_k, false)?,
                        constituent_component(record, 1, second, base_k, false)?,
                    ])
                }
                CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => {
                    // The conjugate stores the *same* row and source number as
                    // the seed and only flips the conjugation flag, so both
                    // components come from the one stored constituent.
                    let conjugate = constituent_component(record, 0, seed.clone(), base_k, true)?;
                    let seed_component = constituent_component(record, 0, seed, base_k, false)?;
                    Ok(vec![seed_component, conjugate])
                }
            }
        }
    }
}

/// One stored CIR constituent as a complex child component.
fn constituent_component(
    record: &'static IrrepRecord,
    index: u8,
    constituent: crate::irrep::types::CompoundConstituentCharacter,
    base_k: Vec3R,
    conjugate: bool,
) -> Result<ChildComponent, FullStarError> {
    let dimension = component_dimension(constituent.dimension, record.sg, constituent.label)?;
    let component = if conjugate {
        SubductionComponent::RealificationConjugate {
            irnumber: constituent.irnumber,
        }
    } else if record.compound_metadata().is_some_and(|metadata| {
        metadata.semantics == CompoundCharacterSemantics::DistinctComponentSum
    }) {
        SubductionComponent::Constituent {
            index,
            irnumber: constituent.irnumber,
        }
    } else {
        SubductionComponent::RealificationSeed {
            irnumber: constituent.irnumber,
        }
    };
    let effective_k = if conjugate {
        base_k.checked_neg()?
    } else {
        base_k
    };
    Ok(ChildComponent {
        record,
        component,
        label: constituent.label,
        bc: record.bc,
        irnumber: constituent.irnumber,
        dimension,
        base_k,
        effective_k,
        conjugate,
        row: constituent.row,
    })
}

/// A component dimension that must fit the `u8` target field.
fn component_dimension(dimension: usize, sg: u8, ml: &'static str) -> Result<u8, FullStarError> {
    u8::try_from(dimension)
        .ok()
        .filter(|value| *value > 0)
        .ok_or_else(|| {
            SubductionError::UnsupportedCharacterSpace {
                sg,
                ml: ml.to_string(),
            }
            .into()
        })
}

/// The representative point of a child star and every component reachable from
/// it.
///
/// The whole orbit is searched in canonical order; the first point with at
/// least one reachable component fixes the frame the block is decomposed in.
/// Components are collected over the **whole** orbit, because a component may
/// sit on another arm of the same child star and then has to be transported to
/// the representative before it is comparable.  Both `q` and the component arms
/// stay unreduced and are compared modulo the child reciprocal lattice.
fn select_representative(
    child_sg: u8,
    folded_star: &FoldedStar,
    child_reciprocal: &Lattice,
) -> Result<Representative, FullStarError> {
    let mut representative: Option<(usize, Vec3R)> = None;
    let mut components: Vec<ChildComponent> = Vec::new();
    for (index, point) in folded_star.points().iter().enumerate() {
        let at_point = child_components_at(child_sg, point.q(), child_reciprocal)?;
        if let Some(first) = at_point.first().filter(|_| representative.is_none()) {
            representative = Some((index, first.effective_k));
        }
        for component in at_point {
            if !components
                .iter()
                .any(|existing| existing.key() == component.key())
            {
                components.push(component);
            }
        }
    }
    let (point, effective_k) =
        representative.ok_or_else(|| FullStarError::MissingChildStarData {
            sg: child_sg,
            q: folded_star
                .points()
                .first()
                .map_or([Rat::ZERO; 3], |point| {
                    [point.q().get(0), point.q().get(1), point.q().get(2)]
                }),
            points: folded_star.points().len(),
        })?;
    Ok(Representative {
        point,
        q: *folded_star.points()[point].q(),
        effective_k,
        components,
    })
}

/// Evaluate every reachable component's little-group row over `H_q`,
/// transported to the representative when its own arm is another point of the
/// child star.
fn prepare_targets(
    child_sg: u8,
    representative: &Representative,
    child_cell: &Lattice,
    child_reciprocal: &Lattice,
    pulled_back: &[ExactSeitz],
) -> Result<Vec<ComplexTarget>, FullStarError> {
    let mut values: Vec<Vec<Complex64>> = Vec::with_capacity(representative.components.len());
    for component in &representative.components {
        values.push(transported_values(
            child_sg,
            component,
            &representative.q,
            child_reciprocal,
            pulled_back,
            child_cell,
        )?);
    }

    // The narrow realification merge: only the seed and conjugate of the *same*
    // record are compared, and only over the full H_q list.  Distinct CIR
    // identities are never candidates.
    let mut merged = vec![false; representative.components.len()];
    for (index, component) in representative.components.iter().enumerate() {
        if merged[index] {
            continue;
        }
        if !matches!(
            component.component,
            SubductionComponent::RealificationSeed { .. }
        ) {
            continue;
        }
        let partner = representative.components.iter().position(|other| {
            other.record.id() == component.record.id()
                && matches!(
                    other.component,
                    SubductionComponent::RealificationConjugate { .. }
                )
        });
        let Some(partner) = partner else {
            continue;
        };
        match realification_relation(child_sg, component.label, &values[index], &values[partner])? {
            RealificationRelation::Equivalent => merged[partner] = true,
            RealificationRelation::Orthogonal => {}
        }
    }

    let mut targets = Vec::with_capacity(representative.components.len());
    for (index, component) in representative.components.iter().enumerate() {
        if merged[index] {
            continue;
        }
        targets.push(ComplexTarget {
            ml: component.label,
            bc: component.bc,
            row_ml: component.record.ml,
            irnumber: component.irnumber,
            dimension: component.dimension,
            component: component.component,
            values: values[index].clone(),
        });
    }
    Ok(targets)
}

/// How the transported seed and conjugate rows of one realification relate.
enum RealificationRelation {
    /// One complex target: the rows agree, so the solver accumulates the
    /// multiplicity of both parent contributions on the seed component.
    Equivalent,
    /// Two complex targets: the rows are orthogonal little-group irreps.
    Orthogonal,
}

/// Compare a realification's transported seed and conjugate rows over `H_q`.
///
/// Both rows must be unit complex irreps; equal unit rows are one target,
/// orthogonal rows are two, and any other overlap is an explicit error rather
/// than a silent merge or a Gram failure.
fn realification_relation(
    sg: u8,
    ml: &'static str,
    seed: &[Complex64],
    conjugate: &[Complex64],
) -> Result<RealificationRelation, FullStarError> {
    let count = seed.len();
    let scale = 1.0 / count as f64;
    for values in [seed, conjugate] {
        let norm = values.iter().map(|value| value.norm_sqr()).sum::<f64>() * scale;
        if (norm - 1.0).abs() > SUBDUCTION_TOLERANCE {
            return Err(FullStarError::RealificationNormMismatch { sg, ml, norm });
        }
    }
    let inner: Complex64 = seed
        .iter()
        .zip(conjugate)
        .map(|(left, right)| left * right.conj())
        .sum::<Complex64>()
        * scale;
    if (inner - Complex64::new(1.0, 0.0)).norm() <= SUBDUCTION_TOLERANCE
        && seed
            .iter()
            .zip(conjugate)
            .all(|(left, right)| (left - right).norm() <= SUBDUCTION_TOLERANCE)
    {
        return Ok(RealificationRelation::Equivalent);
    }
    if inner.norm() <= SUBDUCTION_TOLERANCE {
        return Ok(RealificationRelation::Orthogonal);
    }
    Err(FullStarError::RealificationOverlap { sg, ml, inner })
}

/// A component's row over `H_q`, transported to `q` when its own effective arm
/// is another point of the same child star.
fn transported_values(
    child_sg: u8,
    component: &ChildComponent,
    q: &Vec3R,
    child_reciprocal: &Lattice,
    pulled_back: &[ExactSeitz],
    child_cell: &Lattice,
) -> Result<Vec<Complex64>, FullStarError> {
    let transporter = transport_operation(child_sg, &component.effective_k, q, child_reciprocal)?;
    let mut values = Vec::with_capacity(pulled_back.len());
    for (index, operation) in pulled_back.iter().enumerate() {
        let conjugated = match &transporter {
            Some(transporter) => transporter
                .inverse()?
                .compose(operation)?
                .compose(transporter)?,
            None => *operation,
        };
        let value = character_of(
            &component.row,
            component.label,
            &conjugated,
            child_cell,
            &component.base_k,
            index,
        )?;
        values.push(if component.conjugate {
            value.conj()
        } else {
            value
        });
    }
    Ok(values)
}

/// A child operation transporting an effective arm to the representative,
/// exactly; `None` means the arm already *is* the representative class.
fn transport_operation(
    child_sg: u8,
    effective_k: &Vec3R,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<Option<ExactSeitz>, FullStarError> {
    if child_reciprocal.same_mod(effective_k, q)? {
        return Ok(None);
    }
    for operation in &strict_sg_hall_ops(child_sg)?.operations {
        let image = arm_wave_vector(effective_k, operation)?;
        if child_reciprocal.same_mod(&image, q)? {
            return Ok(Some(*operation));
        }
    }
    Err(FullStarError::MissingComponentTransport {
        sg: child_sg,
        k: Box::new(*effective_k),
        q: Box::new(*q),
    })
}

/// The identity operation of an aligned little-group list, modulo the child
/// lattice.
fn identity_position(
    operations: &[ExactSeitz],
    child_cell: &Lattice,
) -> Result<usize, FullStarError> {
    let position = operations
        .iter()
        .position(|operation| operation.rotation() == IDENTITY_ROTATION)
        .ok_or(FullStarError::MissingIdentityOperation)?;
    if !child_cell.contains(operations[position].translation())? {
        return Err(FullStarError::MissingIdentityOperation);
    }
    Ok(position)
}

/// The prepared complex component behind one reported target.
///
/// The target identity is taken from the stored CIR source, never from a name
/// synthesized out of the row label: the component variant (which carries the
/// CIR source number for every compound variant), the stable label, the physical
/// row and the frozen dimension must all match one reachable component.
fn component_for_target<'a>(
    child_sg: u8,
    components: &'a [ChildComponent],
    target: &SubductionTarget,
) -> Result<&'a ChildComponent, FullStarError> {
    let found = components.iter().find(|component| {
        component.component == target.component
            && component.label == target.ml
            && component.record.ml == target.row_ml
            && component.dimension == target.dimension
    });
    found.ok_or(FullStarError::TargetSourceMismatch {
        sg: child_sg,
        ml: target.ml,
    })
}

/// The frozen CIR source number of an ordinary scalar record.
fn ordinary_irnumber(record: &'static IrrepRecord) -> Result<u32, FullStarError> {
    match record.source_identity() {
        IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => Ok(cir_irnumber),
        _ => Err(SubductionError::UnsupportedCharacterSpace {
            sg: record.sg,
            ml: record.ml.to_string(),
        }
        .into()),
    }
}

/// The child full-star evaluator of one reported target.
///
/// Ordinary rows keep the child's own [`OrdinaryStar`]; every compound
/// component (a `DistinctComponentSum` constituent, or a realification
/// seed/conjugate) is induced over the effective wave vector's star with the
/// component's own conjugation flag.
fn build_evaluator(
    child_sg: u8,
    component: &ChildComponent,
) -> Result<ChildStarEvaluator, FullStarError> {
    match component.component {
        SubductionComponent::Ordinary => Ok(ChildStarEvaluator::Ordinary(OrdinaryStar::new(
            component.record,
        )?)),
        SubductionComponent::Constituent { .. }
        | SubductionComponent::RealificationSeed { .. }
        | SubductionComponent::RealificationConjugate { .. } => {
            Ok(ChildStarEvaluator::Component(ComponentStar::new(
                child_sg,
                component.label,
                component.component,
                component.irnumber,
                component.base_k,
                component.row.clone(),
                component.conjugate,
                None,
            )?))
        }
    }
}

// ── Independent full reconstruction ──────────────────────────────────────────

/// Rebuild the subduced character on every subgroup representative by inducing
/// the reported child irreps over their **child** stars, and compare it with
/// the parent full-star character.
///
/// This deliberately never reruns the q-block builder: the child stars come
/// from the child group's own Hall operations and the child rows' own stored
/// `k`, so a wrong fold or a wrong frame transport cannot cancel out between
/// the two sides.
fn reconstruct(
    embedding: &SubgroupEmbedding,
    star: &ScalarStar,
    blocks: &[FullStarBlock],
) -> Result<(Vec<Complex64>, Vec<Complex64>), FullStarError> {
    let representatives = embedding.representatives();
    let mut parent_characters = Vec::with_capacity(representatives.len());
    let mut pulled_back = Vec::with_capacity(representatives.len());
    for operation in representatives {
        parent_characters.push(star.character(operation)?);
        let child = embedding.transform().unmap_operation(operation)?;
        pulled_back.push(
            shift_operations(
                std::slice::from_ref(&child),
                &embedding.child_shift().checked_neg()?,
            )?[0],
        );
    }
    let mut reconstructed = vec![Complex64::new(0.0, 0.0); representatives.len()];
    for block in blocks {
        for (target, evaluator) in block.targets.iter().zip(&block.evaluators) {
            let weight = f64::from(target.multiplicity);
            for (index, operation) in pulled_back.iter().enumerate() {
                reconstructed[index] += evaluator.character(operation)? * weight;
            }
        }
    }
    for (index, (found, expected)) in reconstructed.iter().zip(&parent_characters).enumerate() {
        if (found - expected).norm() > SUBDUCTION_TOLERANCE {
            return Err(FullStarError::ReconstructionMismatch {
                index,
                found: *found,
                expected: *expected,
            });
        }
    }
    Ok((parent_characters, reconstructed))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
    use crate::irrep::subduction::subduce_irrep_with_embedding;

    fn probe(sg: u8, ml: &str) -> &'static IrrepRecord {
        query::irreps_of(sg)
            .iter()
            .find(|record| record.ml == ml && !record.spinor)
            .unwrap_or_else(|| panic!("SG {sg} has no scalar irrep {ml}"))
    }

    fn embedding(parent: u8, ml: &str, direction: &str) -> SubgroupEmbedding {
        let subgroup =
            isotropy_subgroup_for_direction(parent, ml, IsotropyDirection::Label(direction))
                .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction}: {error}"));
        SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
            .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction} embedding: {error}"))
    }

    fn rat(num: i128, den: i128) -> Rat {
        Rat::new(num, den).expect("a rational with a non-zero denominator")
    }

    fn subgroup_of(parent: u8, ml: &str, direction: &str) -> IsotropySubgroup {
        isotropy_subgroup_for_direction(parent, ml, IsotropyDirection::Label(direction))
            .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction}: {error}"))
    }

    fn result_of(
        parent: u8,
        condensing: &str,
        direction: &str,
        probe_ml: &str,
    ) -> (SubgroupEmbedding, FullStarSubduction) {
        let subgroup = subgroup_of(parent, condensing, direction);
        let built = embedding(parent, condensing, direction);
        let record = probe(parent, probe_ml);
        let result =
            subduce_full_star_with_embedding(&subgroup, &built, record).unwrap_or_else(|error| {
                panic!("SG {parent} {condensing} {direction} {probe_ml}: {error}")
            });
        (built, result)
    }

    /// `(label, dimension, multiplicity, component, irnumber)` of one block.
    fn terms(block: &FullStarBlock) -> Vec<(&'static str, u8, u32, SubductionComponent, u32)> {
        block
            .targets()
            .iter()
            .map(|target| {
                (
                    target.ml,
                    target.dimension,
                    target.multiplicity,
                    target.component,
                    target.irnumber,
                )
            })
            .collect()
    }

    /// `(q, stored k, star size, arm count, block dimension, little dimension)`.
    fn shape(block: &FullStarBlock) -> (Vec3R, Vec3R, usize, usize, u32, u32) {
        (
            *block.q(),
            *block.stored_k(),
            block.star_size(),
            block.arm_count(),
            block.block_dimension(),
            block.little_dimension(),
        )
    }

    fn constituent(index: u8, irnumber: u32) -> SubductionComponent {
        SubductionComponent::Constituent { index, irnumber }
    }

    /// The stored child identity behind one reported target, looked up in the
    /// generated table directly (never through this module).
    fn assert_target_identity(child_sg: u8, target: &FullStarTarget) {
        match target.component {
            SubductionComponent::Ordinary => {
                let record = query::irreps_of(child_sg)
                    .iter()
                    .find(|record| record.ml == target.ml && record.compound_metadata().is_none())
                    .unwrap_or_else(|| {
                        panic!("SG {child_sg} has no ordinary {ml}", ml = target.ml)
                    });
                assert!(matches!(
                    record.source_identity(),
                    IrrepSourceIdentity::OrdinaryScalar { cir_irnumber }
                        if cir_irnumber == target.irnumber
                ));
                let row = record
                    .ordinary_scalar_selected_arm_block_trace()
                    .expect("ordinary row");
                assert_eq!(u32::from(target.dimension), row.dimension() as u32);
                assert_eq!(target.ml, target.row_ml);
            }
            SubductionComponent::Constituent { index, irnumber } => {
                assert_eq!(target.irnumber, irnumber);
                let record = query::irreps_of(child_sg)
                    .iter()
                    .find(|record| record.ml == target.row_ml)
                    .unwrap_or_else(|| panic!("SG {child_sg} has no {ml}", ml = target.row_ml));
                let metadata = record.compound_metadata().expect("compound metadata");
                assert_eq!(metadata.cir_irnumbers[usize::from(index)], irnumber);
                assert_eq!(metadata.cir_labels[usize::from(index)], target.ml);
                assert_eq!(
                    u32::from(metadata.cir_dimensions[usize::from(index)]),
                    u32::from(target.dimension)
                );
                let view = record.compound_selected_arm_view().expect("compound view");
                let CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } = view
                else {
                    panic!("{ml} is not a distinct component sum", ml = target.row_ml);
                };
                let stored = if index == 0 { first } else { second };
                assert_eq!(stored.irnumber, irnumber);
                assert_eq!(stored.label, target.ml);
                assert_eq!(stored.dimension as u32, u32::from(target.dimension));
            }
            SubductionComponent::RealificationSeed { irnumber }
            | SubductionComponent::RealificationConjugate { irnumber } => {
                // Both components share the single stored CIR seed of the row.
                let record = query::irreps_of(child_sg)
                    .iter()
                    .find(|record| record.ml == target.row_ml)
                    .unwrap_or_else(|| panic!("SG {child_sg} has no {ml}", ml = target.row_ml));
                let metadata = record.compound_metadata().expect("compound metadata");
                assert_eq!(
                    metadata.semantics,
                    CompoundCharacterSemantics::ConjugateRealification
                );
                assert_eq!(metadata.cir_irnumbers[0], irnumber);
                assert_eq!(metadata.cir_labels[0], target.ml);
                assert_eq!(
                    u32::from(metadata.cir_dimensions[0]),
                    u32::from(target.dimension)
                );
                let view = record.compound_selected_arm_view().expect("compound view");
                let CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } = view else {
                    panic!("{ml} is not a realification", ml = target.row_ml);
                };
                assert_eq!(seed.irnumber, irnumber);
                assert_eq!(seed.label, target.ml);
                assert_eq!(seed.dimension as u32, u32::from(target.dimension));
            }
        }
    }

    /// Every structural, dimension, source-identity and reconstruction check of
    /// one full-star result.
    fn assert_invariants(result: &FullStarSubduction, parent_sg: u8, child_sg: u8) {
        assert_eq!(result.parent_sg(), parent_sg);
        assert_eq!(result.subgroup_sg(), child_sg);
        assert!(result.parent_dimension() > 0);
        assert!(!result.blocks().is_empty());
        assert_eq!(result.covered_dimension(), result.parent_dimension());
        match result.parent_source_identity() {
            IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => {
                assert_eq!(result.parent_irnumber(), Some(cir_irnumber));
            }
            IrrepSourceIdentity::Compound { .. } => {
                assert_eq!(
                    result.parent_irnumber(),
                    None,
                    "a compound parent has no single CIR source number"
                );
            }
            IrrepSourceIdentity::Spin { .. } => panic!("a spinor reached the scalar decompose"),
        }

        let child_cell =
            Lattice::new(exact_primitive_basis(child_sg).expect("child basis")).expect("lattice");
        let child_reciprocal = child_cell.reciprocal().expect("reciprocal");
        let mut carried = 0u32;
        for block in result.blocks() {
            assert!(block.star_size() > 0);
            assert!(block.arm_count() > 0);
            assert!(block.block_dimension() > 0);
            assert!(!block.targets().is_empty());
            assert!(
                child_reciprocal
                    .same_mod(block.q(), block.stored_k())
                    .expect("same_mod"),
                "q {q:?} and effective k {k:?} differ modulo the child reciprocal lattice",
                q = block.q(),
                k = block.stored_k()
            );
            // The effective arm is the stored k of a child record, or the exact
            // negation of one for a `ConjugateRealification` conjugate
            // component; nothing else may back a block.
            let mut backed = false;
            for record in query::irreps_of(child_sg).iter().filter(|r| !r.spinor) {
                let stored = inline_k_vector(record).expect("k");
                let negated = stored.checked_neg().expect("negation");
                let realification = record.compound_metadata().is_some_and(|metadata| {
                    metadata.semantics == CompoundCharacterSemantics::ConjugateRealification
                });
                if child_reciprocal
                    .same_mod(&stored, block.stored_k())
                    .expect("same_mod")
                    || (realification
                        && child_reciprocal
                            .same_mod(&negated, block.stored_k())
                            .expect("same_mod"))
                {
                    backed = true;
                    break;
                }
            }
            assert!(
                backed,
                "effective k {k:?} is neither a stored child arm nor a realification's -k",
                k = block.stored_k()
            );
            let mut little = 0u32;
            for target in block.targets() {
                assert_eq!(target.sg, child_sg);
                assert!(target.multiplicity > 0);
                assert!(target.dimension > 0);
                little += u32::from(target.dimension) * target.multiplicity;
                assert_target_identity(child_sg, target);
            }
            assert_eq!(little, block.little_dimension());
            assert_eq!(
                little * u32::try_from(block.star_size()).expect("star size"),
                block.block_dimension()
            );
            carried += block.block_dimension();
        }
        assert_eq!(carried, result.parent_dimension());

        let (parent, rebuilt) = result.reconstruction();
        assert!(!parent.is_empty());
        assert_eq!(parent.len(), result.representatives().len());
        assert_eq!(parent.len(), rebuilt.len());
        for (index, (expected, found)) in parent.iter().zip(rebuilt).enumerate() {
            assert!(
                (expected - found).norm() <= result.tolerance(),
                "representative {index}: {found} != {expected}"
            );
        }
    }

    /// 221 `GM4+` P1 -> #83: the **Gamma** probes, including the compound
    /// constituent pair, agree with the stored single-arm decomposition.
    #[test]
    fn golden_gamma_case_reports_the_stored_child_terms() {
        let subgroup = subgroup_of(221, "GM4+", "P1");
        let built = embedding(221, "GM4+", "P1");
        assert_eq!((subgroup.ordinal, subgroup.record.sg), (12400, 83));

        let gm3 = subduce_full_star_with_embedding(&subgroup, &built, probe(221, "GM3+"))
            .expect("GM3+ full star");
        assert_eq!(gm3.blocks().len(), 1);
        assert_eq!(shape(&gm3.blocks()[0]).2, 1);
        assert_eq!(gm3.blocks()[0].q(), &Vec3R::zero());
        assert_eq!(
            terms(&gm3.blocks()[0]),
            [
                ("GM1+", 1, 1, SubductionComponent::Ordinary, 4075),
                ("GM2+", 1, 1, SubductionComponent::Ordinary, 4076),
            ]
        );
        assert_invariants(&gm3, 221, 83);

        let gm4 = subduce_full_star_with_embedding(&subgroup, &built, probe(221, "GM4+"))
            .expect("GM4+ full star");
        assert_eq!(gm4.parent_dimension(), 3);
        assert_eq!(
            terms(&gm4.blocks()[0]),
            [
                ("GM1+", 1, 1, SubductionComponent::Ordinary, 4075),
                ("GM3+", 1, 1, constituent(0, 4077), 4077),
                ("GM4+", 1, 1, constituent(1, 4078), 4078),
            ]
        );
        assert_eq!(gm4.blocks()[0].targets()[1].row_ml, "GM3+GM4+");
        assert_eq!(gm4.blocks()[0].targets()[2].row_ml, "GM3+GM4+");
        assert_invariants(&gm4, 221, 83);
    }

    /// 221 `GM4+` P1 -> #83, `X1+` (source 10719, dimension 3): the three
    /// parent arms split into a two-point child star with a two-arm block and a
    /// one-point star, and the reported terms are the stored child rows.
    #[test]
    fn sg221_x1_star_splits_into_a_two_point_and_a_one_point_star() {
        let (_, result) = result_of(221, "GM4+", "P1", "X1+");
        assert_eq!(result.parent_irnumber(), Some(10719));
        assert_eq!(result.parent_dimension(), 3);
        assert_eq!(result.covered_dimension(), 3);
        assert_eq!(
            *result.seed_k(),
            Vec3R::new([Rat::ZERO, rat(1, 2), Rat::ZERO])
        );
        assert_eq!(
            result.blocks().iter().map(shape).collect::<Vec<_>>(),
            [
                (
                    Vec3R::new([Rat::ZERO, rat(-1, 2), Rat::ZERO]),
                    Vec3R::new([Rat::ZERO, rat(1, 2), Rat::ZERO]),
                    2,
                    2,
                    2,
                    1,
                ),
                (
                    Vec3R::new([Rat::ZERO, Rat::ZERO, rat(-1, 2)]),
                    Vec3R::new([Rat::ZERO, Rat::ZERO, rat(1, 2)]),
                    1,
                    1,
                    1,
                    1,
                ),
            ]
        );
        assert_eq!(
            terms(&result.blocks()[0]),
            [("X1+", 1, 1, SubductionComponent::Ordinary, 4099)]
        );
        assert_eq!(
            terms(&result.blocks()[1]),
            [("Z1+", 1, 1, SubductionComponent::Ordinary, 4115)]
        );
        assert_invariants(&result, 221, 83);
    }

    /// 221 `GM4+` P1 -> #83, `X5+` (source 10723, dimension 6): a child star
    /// with little dimension 2, one ordinary pair and one `DistinctComponentSum`
    /// pair whose constituents keep their own CIR source numbers.
    #[test]
    fn sg221_x5_star_exercises_little_dimension_two_and_constituents() {
        let (_, result) = result_of(221, "GM4+", "P1", "X5+");
        assert_eq!(result.parent_irnumber(), Some(10723));
        assert_eq!(result.parent_dimension(), 6);
        assert_eq!(
            result.blocks().iter().map(shape).collect::<Vec<_>>(),
            [
                (
                    Vec3R::new([Rat::ZERO, rat(-1, 2), Rat::ZERO]),
                    Vec3R::new([Rat::ZERO, rat(1, 2), Rat::ZERO]),
                    2,
                    2,
                    4,
                    2,
                ),
                (
                    Vec3R::new([Rat::ZERO, Rat::ZERO, rat(-1, 2)]),
                    Vec3R::new([Rat::ZERO, Rat::ZERO, rat(1, 2)]),
                    1,
                    1,
                    2,
                    2,
                ),
            ]
        );
        assert_eq!(
            terms(&result.blocks()[0]),
            [
                ("X1+", 1, 1, SubductionComponent::Ordinary, 4099),
                ("X2+", 1, 1, SubductionComponent::Ordinary, 4100),
            ]
        );
        assert_eq!(
            terms(&result.blocks()[1]),
            [
                ("Z3+", 1, 1, constituent(0, 4117), 4117),
                ("Z4+", 1, 1, constituent(1, 4118), 4118),
            ]
        );
        assert_eq!(result.blocks()[1].targets()[0].row_ml, "Z3+Z4+");
        assert_invariants(&result, 221, 83);
    }

    /// 221 `GM4+` P2 -> #12, `X5+`: a child star of size 2 carries `V1+` with
    /// multiplicity 2, so `2 x 1 x 2 = 4` of the six parent dimensions come from
    /// one term of one block.
    #[test]
    fn sg221_p2_x5_star_reports_a_twofold_child_multiplicity() {
        let (_, result) = result_of(221, "GM4+", "P2", "X5+");
        assert_eq!(result.parent_irnumber(), Some(10723));
        assert_eq!(result.parent_dimension(), 6);
        assert_eq!(result.subgroup_sg(), 12);
        assert_eq!(result.ordinal(), 12401);
        let mut reported: Vec<(&str, u32, u8, usize, u32)> = result
            .blocks()
            .iter()
            .flat_map(|block| {
                block.targets().iter().map(|target| {
                    (
                        target.ml,
                        target.irnumber,
                        target.dimension,
                        block.star_size(),
                        target.multiplicity,
                    )
                })
            })
            .collect();
        reported.sort_unstable();
        assert_eq!(
            reported,
            [
                ("A1+", 265, 1, 1, 1),
                ("A2+", 266, 1, 1, 1),
                ("V1+", 269, 1, 2, 2),
            ]
        );
        assert_invariants(&result, 221, 12);
    }

    /// 139 `M1-` P1 -> #126 (ordinal 6213, frozen child shift `(1/4)^3`):
    /// `X1+` merges two parent arms onto one child `q` as one `M1` (dimension
    /// 2), and `N1+` spreads four arms over two `q` points of one star with two
    /// arms each.
    #[test]
    fn sg139_merged_arms_decompose_on_the_merged_little_group() {
        let (built, x1) = result_of(139, "M1-", "P1", "X1+");
        assert_eq!(x1.parent_irnumber(), Some(7301));
        assert_eq!(x1.parent_dimension(), 2);
        assert_eq!(
            shape(&x1.blocks()[0]),
            (
                Vec3R::new([rat(1, 2), rat(1, 2), Rat::ZERO]),
                Vec3R::new([rat(1, 2), rat(1, 2), Rat::ZERO]),
                1,
                2,
                2,
                2,
            )
        );
        assert_eq!(
            terms(&x1.blocks()[0]),
            [("M1", 2, 1, SubductionComponent::Ordinary, 6344)]
        );
        assert_invariants(&x1, 139, 126);
        // Every H operation fixes the sole child q, but some permute the two
        // parent arms at that q. Such operations have zero q-block trace.
        let parent = OrdinaryStar::new(probe(139, "X1+")).unwrap();
        let moved = built
            .representatives()
            .iter()
            .find(|h| {
                !parent
                    .parent_reciprocal
                    .preserves(h.rotation(), parent.seed_k())
                    .unwrap()
            })
            .expect("H_q permutes the parent arms");
        assert!(parent.character(moved).unwrap().norm() < SUBDUCTION_TOLERANCE);

        let (_, n1) = result_of(139, "M1-", "P1", "N1+");
        assert_eq!(n1.parent_irnumber(), Some(7319));
        assert_eq!(n1.parent_dimension(), 4);
        assert_eq!(
            shape(&n1.blocks()[0]),
            (
                Vec3R::new([Rat::ZERO, rat(-1, 2), rat(1, 2)]),
                Vec3R::new([Rat::ZERO, rat(1, 2), rat(1, 2)]),
                2,
                4,
                4,
                2,
            )
        );
        assert_eq!(
            terms(&n1.blocks()[0]),
            [("R1", 2, 1, SubductionComponent::Ordinary, 6354)]
        );
        assert_invariants(&n1, 139, 126);
    }

    /// 167 `GM3+` P1 -> #15 (ordinal 7651): `F1+` splits into a one-point and a
    /// two-point child star.  The folded `q = (1,0,0)` of the first star is
    /// **not** the stored representative arm `(0,1,0)`; only matching modulo the
    /// child reciprocal lattice finds the stored row.  The two points of the
    /// second star differ by `(0,1,0)`, which is *not* a #15 reciprocal lattice
    /// vector (the C-centring extinction), so a per-coordinate `mod 1` merge
    /// would collapse them and break the dimension.
    #[test]
    fn sg167_representative_arm_is_found_modulo_the_child_reciprocal_lattice() {
        let (built, result) = result_of(167, "GM3+", "P1", "F1+");
        assert_eq!(result.parent_irnumber(), Some(8288));
        assert_eq!(result.parent_dimension(), 3);
        assert_eq!(result.subgroup_sg(), 15);
        assert_eq!(result.ordinal(), 7651);
        assert_eq!(
            result.blocks().iter().map(shape).collect::<Vec<_>>(),
            [
                (
                    Vec3R::from_ints([1, 0, 0]),
                    Vec3R::from_ints([0, 1, 0]),
                    1,
                    1,
                    1,
                    1,
                ),
                (
                    Vec3R::new([rat(1, 2), rat(1, 2), rat(1, 2)]),
                    Vec3R::new([rat(1, 2), rat(1, 2), rat(1, 2)]),
                    2,
                    2,
                    2,
                    1,
                ),
            ]
        );
        let child_cell = Lattice::new(exact_primitive_basis(15).expect("basis")).expect("lattice");
        let child_reciprocal = child_cell.reciprocal().expect("reciprocal");
        // The first point of the two-point star has no stored row; selecting
        // only points()[0] would falsely report missing data for this case.
        let folded = ScalarStar::new(probe(167, "F1+"))
            .unwrap()
            .folded_stars(&built)
            .unwrap();
        let two_point = folded.iter().find(|star| star.star_size() == 2).unwrap();
        assert!(
            child_components_at(15, two_point.points()[0].q(), &child_reciprocal)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            select_representative(15, two_point, &child_reciprocal)
                .unwrap()
                .point,
            1
        );
        // The representative-arm witness: the folded q differs from the stored
        // arm by a child reciprocal lattice vector, not by chance.
        let difference = result.blocks()[0]
            .q()
            .checked_sub(result.blocks()[0].stored_k())
            .expect("difference");
        assert_eq!(difference, Vec3R::from_ints([1, -1, 0]));
        assert!(child_reciprocal.contains(&difference).expect("contains"));
        // The centred-class witness: the two points of the second star stay
        // distinct because (0,1,0) is extinct for the C-centred reciprocal
        // lattice, while an integer-cell test would merge them.
        let second = result.blocks()[1].q();
        let first_point = Vec3R::new([rat(1, 2), rat(-1, 2), rat(1, 2)]);
        let separation = second.checked_sub(&first_point).expect("separation");
        assert_eq!(separation, Vec3R::from_ints([0, 1, 0]));
        assert!(!child_reciprocal.contains(&separation).expect("contains"));
        assert!(
            Lattice::integer()
                .contains(&separation)
                .expect("integer lattice contains")
        );
        assert_eq!(
            terms(&result.blocks()[0]),
            [("Y1-", 1, 1, SubductionComponent::Ordinary, 356)]
        );
        assert_eq!(
            terms(&result.blocks()[1]),
            [("L1-", 1, 1, SubductionComponent::Ordinary, 362)]
        );
        assert_invariants(&result, 167, 15);
    }

    /// The phase-sensitive non-Gamma regression of task 7 stays `T3` (not
    /// `T2`) through the full-star entry point: 16 `R2` P1 -> #22, probe `X1`.
    #[test]
    fn phase_sensitive_case_16_r2_to_22_keeps_t3() {
        let (_, result) = result_of(16, "R2", "P1", "X1");
        assert_eq!(result.parent_irnumber(), Some(379));
        assert_eq!(result.subgroup_sg(), 22);
        assert_eq!(
            terms(&result.blocks()[0]),
            [("T3", 1, 1, SubductionComponent::Ordinary, 673)]
        );
        assert_eq!(result.blocks()[0].multiplicity("T2"), 0);
        assert_invariants(&result, 16, 22);
    }

    /// The parent-side characters the reconstruction is checked against are the
    /// archived CIR traces, pinned here from
    /// `tests/data/subduction_star_cir.rs` (which reads the raw matrices).  The
    /// decomposition target is therefore source-backed, not another induced
    /// evaluation.
    #[test]
    fn pinned_parent_star_characters_match_archived_cir_traces() {
        // `(space group, probe label, pinned rotation/trace pairs)`.
        type PinnedTrace = (u8, &'static str, &'static [(Mat3I, f64)]);
        let cases: [PinnedTrace; 4] = [
            (
                221,
                "X1+",
                &[
                    (IDENTITY_ROTATION, 3.0),
                    ([[1, 0, 0], [0, -1, 0], [0, 0, -1]], 3.0),
                    ([[0, 0, 1], [1, 0, 0], [0, 1, 0]], 0.0),
                    ([[0, -1, 0], [1, 0, 0], [0, 0, 1]], 1.0),
                ],
            ),
            (
                221,
                "X5+",
                &[
                    (IDENTITY_ROTATION, 6.0),
                    ([[1, 0, 0], [0, -1, 0], [0, 0, -1]], -2.0),
                    ([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], -2.0),
                    ([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], -2.0),
                ],
            ),
            (
                139,
                "X1+",
                &[
                    (IDENTITY_ROTATION, 2.0),
                    ([[1, 0, 0], [0, -1, 0], [0, 0, -1]], 0.0),
                    ([[0, 1, 0], [1, 0, 0], [0, 0, -1]], 2.0),
                    ([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], 2.0),
                ],
            ),
            (
                139,
                "N1+",
                &[
                    (IDENTITY_ROTATION, 4.0),
                    ([[1, 0, 0], [0, -1, 0], [0, 0, -1]], 2.0),
                    ([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], 2.0),
                    ([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], 0.0),
                ],
            ),
        ];
        for (sg, ml, pinned) in cases {
            let star = OrdinaryStar::new(probe(sg, ml)).expect("star");
            for (rotation, expected) in pinned {
                let operation = strict_sg_hall_ops(sg)
                    .expect("Hall operations")
                    .operations
                    .iter()
                    .copied()
                    .find(|operation| operation.rotation() == *rotation)
                    .unwrap_or_else(|| panic!("SG {sg} has no operation {rotation:?}"));
                assert_eq!(operation.translation(), &Vec3R::zero());
                let found = star.character(&operation).expect("character");
                assert!(
                    (found - Complex64::new(*expected, 0.0)).norm() <= SUBDUCTION_TOLERANCE,
                    "SG {sg} {ml} {rotation:?}: {found} != {expected}"
                );
            }
        }
    }

    /// One pinned full-star result per frozen context goes through every
    /// invariant, so the checks are exercised on ordinary, compound-constituent,
    /// merged-arm and split-star cases.
    #[test]
    fn every_pinned_case_passes_dimension_identity_and_reconstruction_checks() {
        let cases = [
            (221u8, "GM4+", "P1", "X1+", 83u8),
            (221, "GM4+", "P1", "X5+", 83),
            (221, "GM4+", "P1", "GM4+", 83),
            (221, "GM4+", "P2", "GM4+", 12),
            (221, "GM4+", "P2", "X1+", 12),
            (139, "M1-", "P1", "X1+", 126),
            (139, "M1-", "P1", "N1+", 126),
            (167, "GM3+", "P1", "F1+", 15),
            (16, "R1", "P1", "X1", 22),
            (16, "R2", "P1", "X1", 22),
        ];
        for (parent, condensing, direction, ml, child) in cases {
            let (_, result) = result_of(parent, condensing, direction, ml);
            assert_invariants(&result, parent, child);
        }
    }

    /// Reordering the validated transversal changes neither the canonical arm
    /// set nor any block of the decomposition, including for compound parents
    /// whose two components carry their own transversals.
    #[test]
    fn arm_order_invariance() {
        for (parent, condensing, direction, ml) in [
            (221u8, "GM4+", "P1", "X1+"),
            (221, "GM4+", "P1", "X5+"),
            (139, "M1-", "P1", "N1+"),
            (167, "GM3+", "P1", "F1+"),
            (83, "GM1+", "P1", "GM3+GM4+"),
            (23, "GM1", "P1", "W1W1"),
            (45, "GM1", "P1", "W1W1"),
        ] {
            let built = embedding(parent, condensing, direction);
            let record = probe(parent, ml);
            let canonical = ScalarStar::new(record).expect("star");
            let mut transporters: Vec<Vec<ExactSeitz>> = canonical
                .components()
                .iter()
                .map(|component| {
                    component
                        .arms()
                        .iter()
                        .map(|arm| *arm.transporter())
                        .collect::<Vec<_>>()
                })
                .collect();
            for list in &mut transporters {
                list.reverse();
            }
            let permuted =
                ScalarStar::from_transporters(record, &transporters).expect("permuted star");
            assert_eq!(canonical.arms(), permuted.arms());
            assert_eq!(canonical.dimension(), permuted.dimension());
            let expected = decompose_scalar_star(&built, &canonical).expect("canonical");
            let found = decompose_scalar_star(&built, &permuted).expect("permuted");
            assert_eq!(expected.parent_dimension(), found.parent_dimension());
            assert_eq!(expected.representatives(), found.representatives());
            let expected_shape: Vec<_> = expected.blocks().iter().map(shape).collect();
            let found_shape: Vec<_> = found.blocks().iter().map(shape).collect();
            assert_eq!(expected_shape, found_shape, "SG {parent} {ml}");
            let expected_terms: Vec<_> = expected.blocks().iter().map(terms).collect();
            let found_terms: Vec<_> = found.blocks().iter().map(terms).collect();
            assert_eq!(expected_terms, found_terms, "SG {parent} {ml}");
            assert_eq!(
                expected.reconstruction(),
                found.reconstruction(),
                "SG {parent} {ml}"
            );
        }
    }

    /// Listing the two current components of a compound parent in the opposite
    /// order changes neither the parent character nor the reported blocks: the
    /// decomposition is ordered by the child table, not by the parent's
    /// constituent order.
    #[test]
    fn component_order_invariance() {
        for (parent, condensing, direction, ml) in [
            (83u8, "GM1+", "P1", "GM3+GM4+"),
            (23, "GM1", "P1", "W1W1"),
            (45, "GM1", "P1", "W1W1"),
            (167, "GM3+", "P1", "T1T2"),
        ] {
            let built = embedding(parent, condensing, direction);
            let record = probe(parent, ml);
            let canonical = ScalarStar::new(record).expect("star");
            let swapped = canonical.with_components_reversed();
            assert_eq!(canonical.arms().len(), swapped.arms().len());
            assert_eq!(canonical.dimension(), swapped.dimension());
            for operation in &strict_sg_hall_ops(parent).expect("hall").operations {
                let left = canonical.character(operation).expect("character");
                let right = swapped.character(operation).expect("character");
                assert!(
                    (left - right).norm() <= SUBDUCTION_TOLERANCE,
                    "SG {parent} {ml}: {left} != {right}"
                );
            }
            let expected = decompose_scalar_star(&built, &canonical).expect("canonical");
            let found = decompose_scalar_star(&built, &swapped).expect("swapped");
            let expected_terms: Vec<_> = expected.blocks().iter().map(terms).collect();
            let found_terms: Vec<_> = found.blocks().iter().map(terms).collect();
            assert_eq!(expected_terms, found_terms, "SG {parent} {ml}");
        }
    }

    /// Every supported single-arm probe reproduces the established entry point,
    /// and every multi-arm probe still gets `UnsupportedMultiArmStar` from that
    /// entry point: the new stage extends the old one, it does not change it.
    #[test]
    fn single_arm_results_agree_with_the_existing_entry_point() {
        let mut compared = 0usize;
        let mut multi_arm = 0usize;
        for (parent, condensing, direction) in [
            (221u8, "GM4+", "P1"),
            (221, "GM4+", "P2"),
            (16, "R1", "P1"),
            (167, "GM3+", "P1"),
            (139, "M1-", "P1"),
        ] {
            let subgroup = subgroup_of(parent, condensing, direction);
            let built = embedding(parent, condensing, direction);
            for record in query::irreps_of(parent) {
                if record.spinor || record.compound_metadata().is_some() {
                    continue;
                }
                let star = OrdinaryStar::new(record).expect("star");
                if star.arm_count() == 1 {
                    let old = subduce_irrep_with_embedding(&subgroup, &built, record)
                        .unwrap_or_else(|error| panic!("SG {parent} {}: {error}", record.ml));
                    let new = subduce_full_star_with_embedding(&subgroup, &built, record)
                        .unwrap_or_else(|error| panic!("SG {parent} {}: {error}", record.ml));
                    assert_eq!(new.blocks().len(), 1, "SG {parent} {}", record.ml);
                    let block = &new.blocks()[0];
                    assert_eq!(block.star_size(), 1);
                    assert_eq!(u32::from(old.parent_dimension()), block.block_dimension());
                    let old_terms: Vec<(&str, u8, u32, SubductionComponent)> = old
                        .targets()
                        .iter()
                        .map(|target| {
                            (
                                target.ml,
                                target.dimension,
                                target.multiplicity,
                                target.component,
                            )
                        })
                        .collect();
                    let new_terms: Vec<(&str, u8, u32, SubductionComponent)> = block
                        .targets()
                        .iter()
                        .map(|target| {
                            (
                                target.ml,
                                target.dimension,
                                target.multiplicity,
                                target.component,
                            )
                        })
                        .collect();
                    assert_eq!(old_terms, new_terms, "SG {parent} {}", record.ml);
                    assert_eq!(
                        old.folded_k(),
                        [block.q().get(0), block.q().get(1), block.q().get(2)]
                    );
                    compared += 1;
                } else {
                    assert!(
                        matches!(
                            subduce_irrep_with_embedding(&subgroup, &built, record),
                            Err(SubductionError::UnsupportedMultiArmStar { .. })
                        ),
                        "SG {parent} {} must stay unsupported in the single-arm entry point",
                        record.ml
                    );
                    multi_arm += 1;
                }
            }
        }
        assert_eq!(
            (compared, multi_arm),
            (99, 62),
            "single-arm agreement and multi-arm rejection counts"
        );
    }

    /// A changed context, a foreign probe and the unsupported representation
    /// spaces are rejected before any block is built.
    #[test]
    fn context_mismatch_and_unsupported_families_are_rejected() {
        let p1 = subgroup_of(221, "GM4+", "P1");
        let p2 = subgroup_of(221, "GM4+", "P2");
        let p1_embedding = embedding(221, "GM4+", "P1");
        let error = subduce_full_star_with_embedding(&p2, &p1_embedding, probe(221, "GM3+"))
            .expect_err("P2 in the P1 embedding");
        assert!(matches!(
            error,
            FullStarError::Subduction(SubductionError::EmbeddingContextMismatch {
                parent_sg: 221,
                ordinal: 12401,
                subgroup_sg: 12,
                ..
            })
        ));

        let mut tampered = p1;
        tampered.record.origin = [1, 0, 0, 1];
        assert!(matches!(
            subduce_full_star_with_embedding(&tampered, &p1_embedding, probe(221, "GM3+")),
            Err(FullStarError::Subduction(
                SubductionError::StaleIsotropyRecord { ordinal: 12400 }
            ))
        ));

        let foreign = probe(83, "GM1+");
        assert!(matches!(
            subduce_full_star_with_embedding(&p1, &p1_embedding, foreign),
            Err(FullStarError::Subduction(SubductionError::ForeignProbe {
                sg: 221,
                ..
            }))
        ));

        // A compound row of another space group is a foreign probe first: the
        // context check runs before the representation-space adapter.
        let foreign_compound = query::irreps_of(83)
            .iter()
            .find(|record| record.ml == "GM3+GM4+")
            .expect("SG 83 compound row");
        assert!(matches!(
            subduce_full_star_with_embedding(&p1, &p1_embedding, foreign_compound),
            Err(FullStarError::Subduction(SubductionError::ForeignProbe {
                sg: 221,
                ..
            }))
        ));

        let spinor = query::irreps_of(221)
            .iter()
            .find(|record| record.spinor)
            .expect("spinor row");
        assert!(matches!(
            subduce_full_star_with_embedding(&p1, &p1_embedding, spinor),
            Err(FullStarError::Subduction(
                SubductionError::UnsupportedCharacterSpace { sg: 221, .. }
            ))
        ));
    }

    /// The four pinned scalar self-restrictions of task 8c report the frozen
    /// complex components: SG 83 and SG 167 keep both distinct CIR sources,
    /// SG 23 separates the conjugate into its own child star at `-k`, SG 19
    /// merges the equivalent conjugate into the seed with the multiplicity
    /// accumulated, and SG 45 transports the same-star conjugate to the stored
    /// seed arm before comparing.
    #[test]
    fn pinned_scalar_self_restrictions_report_the_complex_components() {
        // ── SG 83 `GM3+GM4+` (DistinctComponentSum, dim 2) ──
        let (_, result) = result_of(83, "GM1+", "P1", "GM3+GM4+");
        assert_eq!(result.subgroup_sg(), 83);
        assert_eq!(
            result.parent_source_identity(),
            IrrepSourceIdentity::Compound {
                metadata_index: 178
            }
        );
        assert_eq!(result.parent_irnumber(), None);
        assert_eq!(result.parent_dimension(), 2);
        assert_eq!(result.blocks().len(), 1);
        assert_eq!(
            terms(&result.blocks()[0]),
            [
                ("GM3+", 1, 1, constituent(0, 4077), 4077),
                ("GM4+", 1, 1, constituent(1, 4078), 4078),
            ]
        );
        assert_invariants(&result, 83, 83);

        // ── SG 23 `W1W1` (ConjugateRealification of a CIR at k with -k in a
        // different star): two child stars, one seed and one conjugate. ──
        let (built, w1) = result_of(23, "GM1", "P1", "W1W1");
        assert_eq!(w1.subgroup_sg(), 23);
        assert_eq!(w1.parent_irnumber(), None);
        assert_eq!(w1.parent_dimension(), 2);
        assert_eq!(built.ordinal(), w1.ordinal());
        assert_eq!(w1.blocks().len(), 2, "k and -k are different child stars");
        let seed_k = Vec3R::new([rat(1, 2), rat(1, 2), rat(1, 2)]);
        let negated = seed_k.checked_neg().expect("negation");
        assert_eq!(w1.blocks()[0].q(), &negated);
        assert_eq!(w1.blocks()[1].q(), &seed_k);
        assert_eq!(w1.blocks()[0].stored_k(), &negated);
        assert_eq!(w1.blocks()[1].stored_k(), &seed_k);
        assert_eq!(
            terms(&w1.blocks()[0]),
            [(
                "W1",
                1,
                1,
                SubductionComponent::RealificationConjugate { irnumber: 716 },
                716
            )]
        );
        assert_eq!(
            terms(&w1.blocks()[1]),
            [(
                "W1",
                1,
                1,
                SubductionComponent::RealificationSeed { irnumber: 716 },
                716
            )]
        );
        for block in w1.blocks() {
            assert_eq!(block.targets()[0].row_ml, "W1W1");
            assert_eq!(block.star_size(), 1);
            assert_eq!(block.arm_count(), 1);
            assert_eq!(block.block_dimension(), 1);
        }
        assert_invariants(&w1, 23, 23);

        // ── SG 19 `R1R1` (self-conjugate star, one arm): the conjugate is an
        // equivalent complex irrep, so the block reports ONE seed target whose
        // multiplicity accumulates both parent contributions. ──
        let (_, r1) = result_of(19, "GM1", "P1", "R1R1");
        assert_eq!(r1.parent_dimension(), 4);
        assert_eq!(r1.blocks().len(), 1);
        assert_eq!(
            terms(&r1.blocks()[0]),
            [(
                "R1",
                2,
                2,
                SubductionComponent::RealificationSeed { irnumber: 559 },
                559
            )]
        );
        assert_eq!(r1.blocks()[0].targets()[0].row_ml, "R1R1");
        assert_eq!(r1.blocks()[0].little_dimension(), 4);
        assert_invariants(&r1, 19, 19);

        // ── SG 45 `W1W1`/`W2W2` (self-conjugate two-arm star): the conjugate
        // lives on the other arm of the SAME child star and has to be
        // transported to the stored seed arm before the rows are compared. ──
        for (probe_ml, cir_label, irnumber) in [("W1W1", "W1", 1805u32), ("W2W2", "W2", 1806)] {
            let (_, w) = result_of(45, "GM1", "P1", probe_ml);
            assert_eq!(w.parent_dimension(), 4);
            assert_eq!(w.blocks().len(), 1);
            assert_eq!(w.blocks()[0].star_size(), 2);
            assert_eq!(w.blocks()[0].arm_count(), 4);
            assert_eq!(
                terms(&w.blocks()[0]),
                [(
                    cir_label,
                    1,
                    2,
                    SubductionComponent::RealificationSeed { irnumber },
                    irnumber
                )]
            );
            assert_eq!(w.blocks()[0].targets()[0].row_ml, probe_ml);
            assert_eq!(w.blocks()[0].little_dimension(), 2);
            assert_invariants(&w, 45, 45);
        }

        // ── SG 167 `T1T2` (DistinctComponentSum, dim 4) through the existing
        // `GM3+` P1 -> #15 embedding. ──
        // ── SG 167 `T1T2` (DistinctComponentSum, dim 4) through the existing
        // `GM3+` P1 -> #15 embedding: both parent constituents restrict onto the
        // same child little-group irrep, so the block reports `M1` twice.  The
        // physical dimension is 2 x 2 = 4. ──
        let (_, t) = result_of(167, "GM3+", "P1", "T1T2");
        assert_eq!(t.subgroup_sg(), 15);
        assert_eq!(t.parent_irnumber(), None);
        assert_eq!(t.parent_dimension(), 4);
        assert_eq!(
            t.blocks().iter().map(terms).collect::<Vec<_>>(),
            [[("M1", 2, 2, SubductionComponent::Ordinary, 363)]]
        );
        assert_invariants(&t, 167, 15);
    }

    /// A mutated folded geometry fails on the dimension checks or on the
    /// missing-data error; it never returns a smaller partial decomposition.
    #[test]
    fn lost_blocks_and_missing_data_fail_instead_of_partial_results() {
        let built = embedding(221, "GM4+", "P1");
        let star = ScalarStar::new(probe(221, "X1+")).expect("star");
        let folded = star.folded_stars(&built).expect("folded");
        assert_eq!(folded.len(), 2);
        assert_eq!(folded[0].star_size(), 2);

        // A lost whole star is a dimension failure: the covered subspace no
        // longer reaches `chi(E)`.
        let dropped = [folded[1].clone()];
        assert!(matches!(
            decompose_folded_stars(&built, &star, &dropped),
            Err(FullStarError::TotalDimensionMismatch {
                expected: 3,
                found: 1
            })
        ));
        assert!(matches!(
            decompose_folded_stars(&built, &star, &[]),
            Err(FullStarError::TotalDimensionMismatch {
                expected: 3,
                found: 0
            })
        ));

        // A lost point inside a star breaks the geometry the kept block claims.
        let mut truncated = folded[0].clone();
        truncated.points.truncate(1);
        assert!(matches!(
            decompose_folded_stars(&built, &star, &[truncated, folded[1].clone()]),
            Err(FullStarError::StarDimensionMismatch { .. })
        ));

        // A lost arm index also breaks the block: no parent arm is left to
        // carry the q-block, and that fails closed instead of asking for a
        // zero-arm character.
        let mut armless = folded[0].clone();
        armless.points[0].arm_indices.clear();
        assert!(matches!(
            decompose_folded_stars(&built, &star, &[armless, folded[1].clone()]),
            Err(FullStarError::EmptyQBlock { .. })
        ));

        // A folded point that matches no stored child record is a typed error,
        // not a zero block.  The shift is not a child reciprocal lattice
        // vector, unlike an integer one.
        let mut foreign = folded[0].clone();
        foreign.points[0].q = foreign.points[0]
            .q
            .checked_add(&Vec3R::new([rat(1, 3), Rat::ZERO, Rat::ZERO]))
            .expect("shift");
        assert!(matches!(
            decompose_folded_stars(&built, &star, &[foreign, folded[1].clone()]),
            Err(FullStarError::MissingChildStarData {
                sg: 83,
                points: 2,
                ..
            })
        ));
    }

    /// The nine fixed contexts scanned probe by probe: the five embedding
    /// contexts of task 8b plus the four scalar self-restrictions of task 8c.
    /// Every scalar probe — ordinary **and** compound — either returns a
    /// complete decomposition that passes every invariant, or is one of the
    /// explicitly counted missing cases; skipped spinor rows are never counted
    /// as successes.  The counts are pinned, so a newly unsupported probe fails
    /// this gate instead of silently shrinking the tested coverage.
    #[test]
    fn frozen_contexts_report_success_missing_and_unsupported_separately() {
        let contexts = [
            (221u8, "GM4+", "P1", 83u8),
            (221, "GM4+", "P2", 12),
            (16, "R1", "P1", 22),
            (167, "GM3+", "P1", 15),
            (139, "M1-", "P1", 126),
            (19, "GM1", "P1", 19),
            (23, "GM1", "P1", 23),
            (45, "GM1", "P1", 45),
            (83, "GM1+", "P1", 83),
        ];
        let mut total_ordinary = 0usize;
        let mut total_compound = 0usize;
        let mut total_spinor = 0usize;
        let mut total_missing = 0usize;
        let mut total_multi_arm = 0usize;
        let mut total_single_arm = 0usize;
        // The first five contexts are the task-8b set, whose ordinary split is
        // also pinned by `single_arm_results_agree_with_the_existing_entry_point`.
        let mut task8b_multi_arm = 0usize;
        let mut task8b_single_arm = 0usize;
        for (context_index, (parent, condensing, direction, child)) in
            contexts.into_iter().enumerate()
        {
            let subgroup = subgroup_of(parent, condensing, direction);
            let built = embedding(parent, condensing, direction);
            assert_eq!(built.subgroup_sg(), child);
            let mut ordinary = 0usize;
            let mut compound = 0usize;
            let mut spinor = 0usize;
            let mut missing = 0usize;
            let mut other = Vec::new();
            for record in query::irreps_of(parent) {
                if record.spinor {
                    spinor += 1;
                    continue;
                }
                let is_compound = record.compound_metadata().is_some();
                match subduce_full_star_with_embedding(&subgroup, &built, record) {
                    Ok(result) => {
                        assert_invariants(&result, parent, child);
                        assert_eq!(result.parent_irnumber().is_none(), is_compound);
                        match ScalarStar::new(record).expect("star").arm_count() {
                            1 => {
                                total_single_arm += 1;
                                if context_index < 5 {
                                    task8b_single_arm += 1;
                                }
                            }
                            _ => {
                                total_multi_arm += 1;
                                if context_index < 5 {
                                    task8b_multi_arm += 1;
                                }
                            }
                        }
                        if is_compound {
                            compound += 1;
                        } else {
                            ordinary += 1;
                        }
                    }
                    Err(FullStarError::MissingChildStarData { .. }) => missing += 1,
                    Err(error) => other.push((record.ml, error)),
                }
            }
            assert!(
                other.is_empty(),
                "SG {parent} {condensing} {direction}: unexpected errors {other:?}"
            );
            assert_eq!(
                (ordinary, compound, spinor, missing),
                match (parent, direction) {
                    (221, "P1") => (40, 0, 20, 0),
                    (221, "P2") => (40, 0, 20, 0),
                    (16, "P1") => (32, 0, 8, 0),
                    (167, "P1") => (12, 1, 14, 0),
                    (139, "P1") => (37, 0, 16, 0),
                    (19, "P1") => (7, 7, 20, 0),
                    (23, "P1") => (14, 4, 9, 0),
                    (45, "P1") => (10, 4, 10, 0),
                    (83, "P1") => (24, 8, 40, 0),
                    _ => unreachable!(),
                },
                "context SG {parent} {condensing} {direction}"
            );
            total_ordinary += ordinary;
            total_compound += compound;
            total_spinor += spinor;
            total_missing += missing;
        }
        assert_eq!(
            (total_ordinary, total_compound, total_spinor, total_missing),
            (216, 24, 157, 0),
            "nine-context scalar decomposition census"
        );
        // Parent-arm coverage: the task-8b five contexts keep their pinned
        // 99 single-arm / 62 multi-arm ordinary split plus the one compound
        // probe (two component arms); the four stronger self-restrictions add
        // the rest.
        assert_eq!((task8b_single_arm, task8b_multi_arm), (99, 63));
        assert_eq!((total_single_arm, total_multi_arm), (138, 102));
    }

    /// The conjugate component's Bloch phase is the `-k` phase, witnessed per
    /// component: SG 23 `W1W1` at the I-centring representative
    /// `(1/2,1/2,1/2)` gives `-i` for the seed and `+i` for the conjugate, and
    /// the physical sum cancels to zero.  Evaluating the unconjugated row at
    /// `-k` instead would give `-i` for both and a spurious `-2i` total.
    #[test]
    fn realification_conjugate_uses_the_negative_k_phase() {
        let star = ScalarStar::new(probe(23, "W1W1")).expect("star");
        let operation = strict_sg_hall_ops(23)
            .expect("Hall operations")
            .operations
            .iter()
            .copied()
            .find(|operation| {
                operation.rotation() == IDENTITY_ROTATION
                    && *operation.translation() == Vec3R::new([rat(1, 2), rat(1, 2), rat(1, 2)])
            })
            .expect("the I-centring representative");
        let seed = star.components()[0].trace(&operation).expect("seed trace");
        let conjugate = star.components()[1]
            .trace(&operation)
            .expect("conjugate trace");
        assert!(
            (seed - Complex64::new(0.0, -1.0)).norm() <= SUBDUCTION_TOLERANCE,
            "seed trace {seed} is not -i"
        );
        assert!(
            (conjugate - Complex64::new(0.0, 1.0)).norm() <= SUBDUCTION_TOLERANCE,
            "conjugate trace {conjugate} is not +i"
        );
        assert!((seed + conjugate).norm() <= SUBDUCTION_TOLERANCE);
        // The physical full-star character reduces the representative modulo
        // the parent lattice, where both components contribute their dimension.
        let identity = star
            .character(&ExactSeitz::identity())
            .expect("identity character");
        assert!((identity - Complex64::new(2.0, 0.0)).norm() <= SUBDUCTION_TOLERANCE);
    }

    #[test]
    fn near_unit_overlap_does_not_replace_pointwise_character_equality() {
        let seed = [Complex64::new(1.0, 0.0); 4];
        let phase = Complex64::from_polar(1.0, 1e-4);
        let corrupted = [seed[0], phase, phase.conj(), seed[3]];
        let inner: Complex64 = seed
            .iter()
            .zip(corrupted)
            .map(|(left, right)| left * right.conj())
            .sum::<Complex64>()
            / 4.0;
        // Unit Gram norms and an overlap within tolerance do not bound each
        // entry to that tolerance: here the residual is of order sqrt(epsilon).
        assert!((inner - 1.0).norm() < SUBDUCTION_TOLERANCE);
        assert!((phase - 1.0).norm() > SUBDUCTION_TOLERANCE);
        assert!(matches!(
            realification_relation(19, "R1", &seed, &corrupted),
            Err(FullStarError::RealificationOverlap { .. })
        ));

        // Match the shared solver's Gram diagonal tolerance, without taking a
        // square root that would allow almost twice its error on a discarded row.
        let inflated = [Complex64::new(1.0 + 0.75 * SUBDUCTION_TOLERANCE, 0.0); 4];
        assert!(matches!(
            realification_relation(19, "R1", &seed, &inflated),
            Err(FullStarError::RealificationNormMismatch { .. })
        ));
    }

    /// SG 19 `R1R1`: the transported conjugate row equals the seed row, so the
    /// adapter reports one target.  Feeding both rows to the shared solver
    /// instead must fail the Gram identity rather than double-counting.
    #[test]
    fn equivalent_realification_rows_merge_and_duplicates_fail_the_gram_check() {
        let built = embedding(19, "GM1", "P1");
        let star = ScalarStar::new(probe(19, "R1R1")).expect("star");
        let folded = star.folded_stars(&built).expect("folded");
        assert_eq!(folded.len(), 1);
        let child_cell = Lattice::new(exact_primitive_basis(19).expect("basis")).expect("cell");
        let child_reciprocal = child_cell.reciprocal().expect("reciprocal");
        let representative =
            select_representative(19, &folded[0], &child_reciprocal).expect("representative");
        assert_eq!(representative.components.len(), 2);
        let (parent_operations, pulled_back) =
            little_group_operations(&built, &representative.q, &child_reciprocal).expect("H_q");
        let parent_characters: Vec<Complex64> = parent_operations
            .iter()
            .map(|operation| {
                star.q_block_character(
                    folded[0].points()[representative.point].arm_indices(),
                    operation,
                )
                .expect("q-block character")
            })
            .collect();
        let merged = prepare_targets(
            19,
            &representative,
            &child_cell,
            &child_reciprocal,
            &pulled_back,
        )
        .expect("prepared targets");
        assert_eq!(merged.len(), 1, "the equivalent conjugate is one target");
        assert!(matches!(
            merged[0].component,
            SubductionComponent::RealificationSeed { irnumber: 559 }
        ));

        // The same one target twice is not an orthogonal set.
        let mut duplicated = Vec::new();
        for _ in 0..2 {
            duplicated.push(ComplexTarget {
                ml: merged[0].ml,
                bc: merged[0].bc,
                row_ml: merged[0].row_ml,
                irnumber: merged[0].irnumber,
                dimension: merged[0].dimension,
                component: merged[0].component,
                values: merged[0].values.clone(),
            });
        }
        assert!(matches!(
            solve_prepared_character_block(19, duplicated, &parent_characters, &pulled_back),
            Err(SubductionError::TargetRowsNotOrthogonal { .. })
        ));
        let solved = solve_prepared_character_block(19, merged, &parent_characters, &pulled_back)
            .expect("merged solution");
        assert_eq!(solved.targets.len(), 1);
        assert_eq!(solved.targets[0].multiplicity, 2);
    }

    /// Distinct CIR identities are never merged, even inside one
    /// `DistinctComponentSum` row: SG 83 reports both constituents, and the
    /// realification merge is bound to the seed/conjugate pair of one record.
    #[test]
    fn distinct_component_sources_are_never_merged() {
        let built = embedding(83, "GM1+", "P1");
        let star = ScalarStar::new(probe(83, "GM3+GM4+")).expect("star");
        let folded = star.folded_stars(&built).expect("folded");
        assert_eq!(folded.len(), 1);
        let child_cell = Lattice::new(exact_primitive_basis(83).expect("basis")).expect("cell");
        let child_reciprocal = child_cell.reciprocal().expect("reciprocal");
        let representative =
            select_representative(83, &folded[0], &child_reciprocal).expect("representative");
        // Every Gamma record of SG 83 is reachable from the Gamma star; the
        // block's catalogue therefore holds all of them.
        assert_eq!(representative.components.len(), 8);
        let (_, pulled_back) =
            little_group_operations(&built, &representative.q, &child_reciprocal).expect("H_q");
        let targets = prepare_targets(
            83,
            &representative,
            &child_cell,
            &child_reciprocal,
            &pulled_back,
        )
        .expect("prepared targets");
        let pair: Vec<&ComplexTarget> = targets
            .iter()
            .filter(|target| target.row_ml == "GM3+GM4+")
            .collect();
        assert_eq!(pair.len(), 2, "the two distinct constituents stay separate");
        assert_eq!(pair[0].irnumber, 4077);
        assert_eq!(pair[1].irnumber, 4078);
        assert_eq!(
            pair.iter()
                .map(|target| target.component)
                .collect::<Vec<_>>(),
            [
                SubductionComponent::Constituent {
                    index: 0,
                    irnumber: 4077
                },
                SubductionComponent::Constituent {
                    index: 1,
                    irnumber: 4078
                }
            ]
        );
    }

    /// A parent star with one component or one arm removed fails the block
    /// tiling check instead of returning a smaller decomposition.
    #[test]
    fn lost_components_and_arms_fail_the_block_tiling() {
        let built = embedding(83, "GM1+", "P1");
        let star = ScalarStar::new(probe(83, "GM3+GM4+")).expect("star");
        assert_eq!(star.dimension(), 2);
        for mutated in [star.without_last_component(), star.without_last_arm()] {
            assert!(matches!(
                decompose_scalar_star(&built, &mutated),
                Err(FullStarError::TotalDimensionMismatch {
                    expected: 2,
                    found: 1
                })
            ));
        }
    }
}
