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
use crate::irrep::w_little_characters_data::LittleCharacterTable;
use crate::irrep::types::{
    CharacterRow, CompoundCharacterSemantics, CompoundSelectedArmCharacter, IrrepRecord,
    IrrepSourceIdentity,
};
use crate::mathfunc::Mat3I;

use super::super::{
    ComplexTarget, ExactSeitz, Lattice, Mat3R, Rat, SUBDUCTION_TOLERANCE, SubductionComponent,
    SubductionError, SubductionTarget, SubgroupEmbedding, Vec3R, bloch_phase, character_of,
    exact_primitive_basis, inline_k_vector, shift_operations,
    solve_prepared_character_block, validate_record_and_embedding,
    strict_sg_hall_ops, validate_subduction_context,
};
use super::scalar_star::{ComponentStar, ConstructedStar, ScalarStar};
use super::{
    ConstructedLittleRep, FoldArm, FoldedStar, OrdinaryStar, StarError, arm_wave_vector,
    catalogue, fold_arms,
};

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
    /// A parametric-k source table does not belong to the embedding it was
    /// asked to answer for.
    #[error(
        "ordinal {ordinal}: the frozen line source {label} of space group {source_sg} does not \
         belong to parent space group {parent}"
    )]
    LineSourceMismatch {
        /// Isotropy record ordinal of the embedding.
        ordinal: usize,
        /// Parent space group of the embedding.
        parent: u8,
        /// Space group the frozen table belongs to.
        source_sg: u8,
        /// Frozen source label.
        label: &'static str,
    },
    /// A frozen line source carries a direction or a rotation the little group
    /// of the parent's operation list does not reproduce.
    #[error("frozen line source {label} of space group {sg} is inconsistent with the parent's little group")]
    MissingLineRotation {
        /// Parent space group number.
        sg: u8,
        /// Frozen source label.
        label: &'static str,
    },
    /// The frozen direction of a line source is not a parsable rational vector.
    #[error("frozen line source {label} of space group {sg} has an unreadable direction")]
    InvalidLineDirection {
        /// Parent space group number.
        sg: u8,
        /// Frozen source label.
        label: &'static str,
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
    ///
    /// Table-corruption guard, **unreachable through the public API on the pinned
    /// data**: both readings compared here are derived from the same target list
    /// (`ml` and `irnumber` come from one `ChildComponent`), so the branch only
    /// fires if the frozen tables are inconsistent with each other.  No negative
    /// test can be written against it without fabricating such a table.
    #[error("child target {ml} of space group {sg} disagrees with its stored CIR source")]
    TargetSourceMismatch {
        /// Child space group number.
        sg: u8,
        /// Target label.
        ml: &'static str,
    },
    /// The child table has no unique trivial row, so the trivial content of a
    /// subduction onto that child is not defined by the table.
    ///
    /// Table-corruption guard, unreachable through the public API on the pinned
    /// data: an all-ones one-dimensional Gamma row is the trivial representation
    /// of the child group in its own stored setting, and exactly one exists for
    /// every space group (pinned for all 230 by
    /// `every_space_group_has_one_trivial_gamma_row`), so this variant reports a
    /// broken child table rather than a physical possibility.
    #[error("space group {sg} has no unique trivial Gamma irrep")]
    MissingChildTrivialIrrep {
        /// Child space group number.
        sg: u8,
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
    /// A realification seed and its conjugate share the one CIR label.  `None`
    /// for a constructed target: no pinned row names it, so no label is
    /// borrowed or invented.
    pub ml: Option<&'static str>,
    /// Bradley-Cracknell label, for display only; `None` when the target has no
    /// sourced name.
    pub bc: Option<&'static str>,
    /// Physical child row this term came from; `None` for a constructed target.
    pub row_ml: Option<&'static str>,
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
    /// shared seed `irnumber` for both realification components); `None` for a
    /// constructed target, whose identity is its component and exact point.
    pub irnumber: Option<u32>,
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
    ///
    /// This is the **unreduced** folded coordinate `T^T k`: two parameters that
    /// differ by a parent reciprocal vector give child cosets that differ by a
    /// child reciprocal vector, so `q()` itself is not comparable across
    /// parameters (`t = 1/4` and `t = 5/4` differ here on 5,716 of the 5,756
    /// pinned rows while every other field agrees).  Compare
    /// [`LineSubduction::wave_vector`]
    /// or the reported `(stored_k, star, dim, arm, targets)` tuple instead, or
    /// reduce `q()` into the child's own reciprocal cell first.
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
    /// positive representative.  For a **constructed** block -- the ordinary
    /// case at a generic parameter, where no child table row exists at all --
    /// there is no stored `k` to report and this is the block's folded
    /// coordinate reduced into the child's reciprocal cell, i.e. the canonical
    /// child coset the constructed Bloch target lives on.
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

    /// Multiplicity of one sourced target label, or `0` when it does not appear.
    /// A constructed target has no label and is never matched here.
    pub fn multiplicity(&self, ml: &str) -> u32 {
        self.targets
            .iter()
            .find(|target| target.ml == Some(ml))
            .map_or(0, |target| target.multiplicity)
    }

    /// Multiplicity of the constructed target with this **exact identity**, or
    /// `0` when the block reports no such target or the identity is not constructed.
    ///
    /// A constructed target's identity is the pair
    /// [`SubductionComponent::Constructed`] carries: its point **reduced modulo
    /// the child reciprocal lattice** and its index in the constructed list at
    /// that point, exactly as reported in [`FullStarTarget::component`].  The
    /// point is not the raw folded representative: `(-1/4, -1/4, -1)` and
    /// `(3/4, 3/4, 0)` are the same target, and only the canonical form is the
    /// identity, so looking one up by a raw folded coordinate would silently
    /// answer zero.
    pub fn constructed_multiplicity(&self, identity: SubductionComponent) -> u32 {
        if !matches!(identity, SubductionComponent::Constructed { .. }) {
            return 0;
        }
        self.targets
            .iter()
            .find(|target| target.component == identity)
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
    /// Shared denominator of `setting`: `U = setting / setting_denominator`.
    setting_denominator: i32,
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

    /// Numerator of the setting transform of the embedding this result belongs
    /// to; the exact matrix is `setting() / setting_denominator()`.
    pub const fn setting(&self) -> Mat3I {
        self.setting
    }

    /// Denominator of the setting transform, always positive.
    pub const fn setting_denominator(&self) -> i32 {
        self.setting_denominator
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
    Ordinary(Box<OrdinaryStar>),
    /// One complex component of a compound child row (a
    /// `DistinctComponentSum` constituent, or a realification seed/conjugate
    /// with its effective `k` and conjugation flag).
    Component(Box<ComponentStar>),
    /// A little-group representation constructed at an exact folded point,
    /// induced over the child's own star.  For a child with a trivial point
    /// group the star is the point itself, so the evaluator is the
    /// little-group character.
    Constructed(Box<ConstructedStar>),
}

impl ChildStarEvaluator {
    fn character(&self, operation: &ExactSeitz) -> Result<Complex64, FullStarError> {
        match self {
            Self::Ordinary(star) => Ok(star.character(operation)?),
            Self::Component(star) => Ok(star.character(operation)?),
            Self::Constructed(star) => Ok(star.character(operation)?),
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

/// How often the trivial representation of the subgroup appears in the
/// subduction of one scalar parent probe.
///
/// This is the quantity the pinned `isotropy_subduce_*` table stores as the
/// subduction frequency of a parent irrep, and it is the only quantity that
/// table constrains.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrivialContent {
    /// Multiplicity read through the frozen CIR source number of the child's
    /// trivial row: the sum over the evaluated Gamma blocks of the multiplicity
    /// of every target that carries that source.  This is the number the stored
    /// table lists.
    pub total: u32,
    /// The same sum read through the child row's Miller-Love label.  The engine
    /// fails closed with [`FullStarError::TargetSourceMismatch`] unless the two
    /// readings agree, so a caller may use either.
    pub by_label: u32,
    /// Folded child stars evaluated because the subgroup's Gamma point is one of
    /// their points.
    pub gamma_stars: usize,
    /// Folded child stars skipped because none of their points is the
    /// subgroup's Gamma point.  They cannot contribute; see the theorem in
    /// [`trivial_content_with_embedding`].
    pub skipped_stars: usize,
}

/// The star arms of a parametric-k line as a `build_block` character source.
///
/// Mirrors the two methods `build_block` calls on `ScalarStar`: the dimension of
/// a q-block is the little dimension times its arms, and its character sums the
/// transported little character over those arms (`line_character`'s logic,
/// inlined here because `FullStarError` does not convert into `StarError`).
struct LineArmSource<'a> {
    table: &'static LittleCharacterTable,
    direction: Vec3R,
    wave_vector: Vec3R,
    arms: &'a [(Vec3R, Mat3I)],
}

impl LineArmSource<'_> {
    /// Dimension of the parent q-block carried by `arm_indices`.
    fn q_block_dimension(&self, arm_indices: &[usize]) -> Result<u32, StarError> {
        let arms = u32::try_from(arm_indices.len()).map_err(|_| {
            StarError::Subduction(SubductionError::RationalOverflow {
                operation: "line q-block arms",
            })
        })?;
        u32::from(self.table.dimension)
            .checked_mul(arms)
            .ok_or(StarError::Subduction(SubductionError::RationalOverflow {
                operation: "line q-block dimension",
            }))
    }

    /// Character of the parent q-block on one parent operation.
    fn q_block_character(
        &self,
        arm_indices: &[usize],
        operation: &ExactSeitz,
    ) -> Result<Complex64, StarError> {
        let mut value = Complex64::new(0.0, 0.0);
        for index in arm_indices {
            // The index guards an internal caller and cannot fail in practice.
            let Some((_, rotation)) = self.arms.get(*index) else {
                return Err(StarError::OperationNotInParentGroup {
                    sg: self.table.space_group,
                });
            };
            let transport = ExactSeitz::new(*rotation, Vec3R::new([Rat::ZERO; 3]));
            let conjugate = transport
                .inverse()?
                .compose(operation)?
                .compose(&transport)?;
            let image =
                Mat3R::from_ints(conjugate.rotation()).checked_mul_vector(&self.direction)?;
            if image != self.direction {
                // The operation moves this arm; it contributes no diagonal term.
                continue;
            }
            let character = self
                .table
                .operations
                .iter()
                .find(|frozen| {
                    let frozen_rotation: Mat3I = frozen.rotation.map(|row| row.map(i32::from));
                    frozen_rotation == conjugate.rotation()
                })
                .map(|frozen| frozen.character)
                .ok_or(StarError::MissingFrozenRotation {
                    sg: self.table.space_group,
                    label: self.table.label,
                })?;
            value += bloch_phase(&self.wave_vector, conjugate.translation())?
                * Complex64::new(f64::from(character[0]), f64::from(character[1]));
        }
        Ok(value)
    }
}

/// Multiplicity of the subgroup's **trivial** representation in the subduction
/// of one scalar parent probe.
///
/// Unlike [`subduce_full_star_with_embedding`] this stays exact when a folded
/// child star has no stored child data, because such a star provably cannot
/// contribute:
///
/// * every representation carried at the folded wave vector `q` acts on a child
///   lattice translation `t` as the scalar `exp(-2 pi i q.t)`, so the parent
///   star's character on that coset is a phase times its value at `q = 0`;
/// * the trivial representation acts on every translation as `1`;
/// * a common irreducible constituent of the two therefore needs
///   `exp(-2 pi i q.t) = 1` for every translation `t` of the child lattice,
///   i.e. `q = 0` modulo the **child** reciprocal lattice (centring extinctions
///   included, exactly as [`Lattice::contains`] decides it);
///
/// so every folded star whose points are not the subgroup's Gamma point
/// contributes zero to the trivial content -- no child `k`/character data is
/// needed for it, and skipping it is not an approximation.  The Gamma blocks
/// themselves are built by the same `build_block` stage the full decomposition
/// uses, with the same character pairing and the same frozen child origin
/// shift.
///
/// `subgroup`, `embedding` and `probe` are revalidated as one context first.
/// Spinor probes are rejected exactly as in the full entry point, and a missing
/// or ambiguous trivial child row is a
/// [`FullStarError::MissingChildTrivialIrrep`] error rather than a silent zero.
pub fn trivial_content_with_embedding(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<TrivialContent, FullStarError> {
    validate_subduction_context(subgroup, embedding, probe)?;
    if probe.spinor {
        return Err(SubductionError::UnsupportedCharacterSpace {
            sg: probe.sg,
            ml: probe.ml.to_string(),
        }
        .into());
    }
    let child_sg = embedding.subgroup_sg();
    let trivial = trivial_child_record(child_sg)
        .ok_or(FullStarError::MissingChildTrivialIrrep { sg: child_sg })?;
    let trivial_cir = match trivial.source_identity() {
        IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => cir_irnumber,
        IrrepSourceIdentity::Compound { .. } | IrrepSourceIdentity::Spin { .. } => {
            return Err(FullStarError::MissingChildTrivialIrrep { sg: child_sg });
        }
    };

    let star = ScalarStar::new(probe)?;
    let folded = star.folded_stars(embedding)?;
    let child_cell = Lattice::new(exact_primitive_basis(child_sg)?)?;
    let child_reciprocal = child_cell.reciprocal()?;

    let mut content = TrivialContent {
        total: 0,
        by_label: 0,
        gamma_stars: 0,
        skipped_stars: 0,
    };
    for folded_star in &folded {
        let mut is_gamma = false;
        for point in folded_star.points() {
            if child_reciprocal.contains(point.q())? {
                is_gamma = true;
                break;
            }
        }
        if !is_gamma {
            content.skipped_stars += 1;
            continue;
        }
        content.gamma_stars += 1;
        let block = build_block(
            embedding,
            &ArmCharacterSource::Cir(&star),
            folded_star,
            &child_cell,
            &child_reciprocal,
        )?;
        content.by_label += block.multiplicity(trivial.ml);
        for target in block.targets() {
            if target.dimension == trivial.dim && target.irnumber == Some(trivial_cir) {
                content.total += target.multiplicity;
            }
        }
    }
    if content.by_label != content.total {
        return Err(FullStarError::TargetSourceMismatch {
            sg: child_sg,
            ml: trivial.ml,
        });
    }
    Ok(content)
}

// ── Parametric-k sources (`other_wave_vector_subduction`) ────────────────────

/// Free-parameter value the official program uses for a parametric-k domain.
///
/// A line irrep's character on a little-group element is the Gamma-point value
/// times a Bloch phase, `chi(R, T) = D(R) * exp(2 pi i t (v . T))`, so the
/// frequencies the pinned `isotropy_w_subduce_*` rows carry belong to one
/// particular `t`.  Measured on the 46 `P1`-child w records (300 pinned rows):
/// `t = 1/4` and `t = 3/4` reproduce every row, `t = 1/2` and `t = 1` leave
/// 111 rows wrong, and `t = 1/8, 3/8, 1/6, 1/12` leave all 300 wrong.  A
/// different parameter is a different representation of the *parent* group, so
/// this is the program's convention rather than an engine choice; it is pinned
/// here and the whole pinned table validates it in the audit.
pub const OFFICIAL_LINE_PARAMETER: (i128, i128) = (1, 4);

/// The official parameter as an exact rational.
///
/// The constant cannot fail `Rat::new`; the `Result` exists so callers never
/// unwrap.  The audit and the R6 tests all build the value from here, so moving
/// the constant is felt by every pinned comparison instead of being pinned twice.
pub fn official_line_parameter() -> Result<Rat, SubductionError> {
    Rat::new(OFFICIAL_LINE_PARAMETER.0, OFFICIAL_LINE_PARAMETER.1)
}

/// The frozen table's direction as an exact rational vector.
fn line_direction(table: &LittleCharacterTable) -> Result<Vec3R, FullStarError> {
    let mut values = [Rat::ZERO; 3];
    for (axis, text) in table.direction.iter().enumerate() {
        let text = text.trim();
        let (numerator, denominator) = match text.split_once('/') {
            Some((numerator, denominator)) => (numerator, denominator),
            None => (text, "1"),
        };
        let numerator = numerator
            .parse::<i128>()
            .map_err(|_| FullStarError::InvalidLineDirection {
                sg: table.space_group,
                label: table.label,
            })?;
        let denominator = denominator
            .parse::<i128>()
            .map_err(|_| FullStarError::InvalidLineDirection {
                sg: table.space_group,
                label: table.label,
            })?;
        values[axis] = Rat::new(numerator, denominator)
            .map_err(|_| FullStarError::InvalidLineDirection {
                sg: table.space_group,
                label: table.label,
            })?;
    }
    Ok(Vec3R::new(values))
}



/// Multiplicity of the subgroup's **trivial** representation in the subduction
/// of one parametric-k parent source (`other_wave_vector_subduction`).
///
/// Delegates to [`line_trivial_content_via_blocks`], which folds the line's arms
/// through the same `build_block` stage the discrete probes use.  The earlier
/// hand-written per-arm sum is gone, and with it the possibility of a caller
/// picking whichever of two answers it liked: that implementation disagreed with
/// the pinned rows on 1,599 values and returned errors on 241 more (measured
/// during the R5 review round and recorded in `docs/task9-remaining-work.md`,
/// "Review fixes (post-113)"; the code itself has since been deleted, so the two
/// numbers are a historical record rather than a recomputable check).
pub fn line_trivial_content_with_embedding(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
) -> Result<u32, FullStarError> {
    line_trivial_content_via_blocks(subgroup, embedding, table)
}

/// The arms of a parametric-k line: the images of its direction under the
/// parent's own operations, deduplicated by **exact vector equality** -- not
/// modulo the parent reciprocal lattice, so `a` and `a + G` are two arms.
///
/// That is the frozen convention the pinned rows were produced with: for the 73
/// frozen sources the two readings do differ (2,574 arm pairs differ by a parent
/// reciprocal vector), and the exact one is the one whose
/// `little dimension x arms` reproduces the pinned full-star dimensions 73/73.
/// Every arm still lies in the centre's coset modulo the parent reciprocal
/// lattice (50,226/50,226 arms over all 5,756 rows), which is why the folded
/// child stars are unaffected by the choice.
///
/// The direction comes from the frozen little-character table (the frame the
/// official program prints it in); the same vector is what the parameter
/// multiplies in [`line_wave_vector`], so arms, characters and folding all speak
/// one frame.
fn line_arms(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    direction: &Vec3R,
) -> Result<Vec<(Vec3R, Mat3I)>, FullStarError> {
    let parent_lattice = embedding.parent_lattice();
    let parent_ops =
        parent_lattice.deduplicate(&strict_sg_hall_ops(subgroup.parent_sg)?.operations)?;
    let mut arms: Vec<(Vec3R, Mat3I)> = Vec::new();
    for operation in &parent_ops {
        let rotation = operation.rotation();
        let action = Mat3R::from_ints(rotation).inverse()?.transpose();
        let image = action.checked_mul_vector(direction)?;
        if arms.iter().any(|(arm, _)| *arm == image) {
            continue;
        }
        arms.push((image, rotation));
    }
    if arms.is_empty() {
        return Err(FullStarError::MissingLineRotation {
            sg: table.space_group,
            label: table.label,
        });
    }
    Ok(arms)
}

/// `k = t . direction` as an exact vector; `t` is the caller's rational
/// parameter and `direction` the frozen table's own direction.
///
/// **This is the wave vector the frozen table is evaluated at, unreduced.**  The
/// frozen little co-group matrices `D` are solved at `k = Gamma` (no Bloch
/// factor), so the character of the band on an operation `(R, T)` is
/// `D(R) * exp(2 pi i k.T)`: the parameter *is* part of the representation, and
/// reducing `k` into the parent's fundamental cell while keeping `D` would
/// silently evaluate a different band.
///
/// What a parameter step `t -> t + n` does is therefore a **relabelling**, not a
/// gauge: `chi^t+n = chi^t * Phi_{n v}` with `Phi_K(R) = exp(2 pi i K.T_R)` a
/// one-dimensional little-group character, so the shifted table is the table of
/// the image under the monodromy map `M_{n v}` of
/// [`crate::irrep::line_monodromy`] and
///
/// ```text
/// decompose(alpha, t + n)  ==  decompose(M_{n v}(alpha), t).
/// ```
///
/// R6.1 shipped the opposite reading (reduce `k`, keep the label) and the R6.2
/// audit revoked it.  Witness (ordinal 11328, SG 210, measured on the current
/// build): `M_v(DT3) = DT4`; the reduced engine answered `Z1` for `DT3` at both
/// parameters, while the correct transport gives `Z2` -- the decomposition the
/// frozen table of `DT4` has at `t = 1/4`, and exactly what `DT3` at `t = 5/4`
/// now returns (`tests/line_monodromy.rs` checks this on all 5,756 rows).
fn line_wave_vector(direction: &Vec3R, parameter: &Rat) -> Result<Vec3R, FullStarError> {
    let mut scaled = [Rat::ZERO; 3];
    for (axis, value) in scaled.iter_mut().enumerate() {
        *value = parameter.checked_mul(direction.get(axis))?;
    }
    Ok(Vec3R::new(scaled))
}

/// Whether a parameter keeps the frozen line little co-group or makes distinct
/// line arms coincide at the same parent wave vector.
fn line_parameter_kind(
    parent_lattice: &Lattice,
    table: &'static LittleCharacterTable,
    direction: &Vec3R,
    wave_vector: &Vec3R,
    parameter: &Rat,
    arms: &[(Vec3R, Mat3I)],
) -> Result<ParameterKind, FullStarError> {
    if !arms.iter().any(|(arm, _)| arm == direction) {
        return Err(FullStarError::MissingLineRotation {
            sg: table.space_group,
            label: table.label,
        });
    }
    for (arm, _) in arms {
        if arm == direction {
            continue;
        }
        let arm_wave_vector = line_wave_vector(arm, parameter)?;
        let difference = arm_wave_vector.checked_sub(wave_vector)?;
        // For a direct lattice with row basis B, q belongs to its reciprocal
        // lattice exactly when B q has integer coordinates. This preserves
        // centring extinctions without constructing (B^-1)^T for every query.
        if parent_lattice
            .rows()
            .checked_mul_vector(&difference)?
            .is_integral()
        {
            return Ok(ParameterKind::Formal);
        }
    }
    Ok(ParameterKind::LineIrrep)
}

/// Validate the context shared by every line entry point.
fn validate_line_context(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
) -> Result<(), FullStarError> {
    // Reuse the production context validation: the record must still be the
    // stored one (a mutated basis, origin, subgroup number or irrep context is a
    // `StaleIsotropyRecord`) and the embedding must have been built from it.
    validate_record_and_embedding(subgroup, embedding)?;
    // The line-specific half: the frozen table has to belong to this parent.
    if usize::from(table.space_group) != usize::from(subgroup.parent_sg) {
        return Err(FullStarError::LineSourceMismatch {
            ordinal: subgroup.ordinal,
            parent: subgroup.parent_sg,
            source_sg: table.space_group,
            label: table.label,
        });
    }
    Ok(())
}

/// Fold a parametric-k line's arms into the child's zone and partition them into
/// **child stars**.
///
/// The line's arm list is `(image of the direction, parent rotation)` at one
/// parameter value; this turns it into the same `FoldArm` geometry the discrete
/// star adapters use, so the folded blocks are child *orbits* (one block per
/// orbit, each carrying every arm that folds onto one of its points) and the
/// reconstruction invariant is the discrete one.  Treating each folded q as its
/// own star instead would count a multi-point child orbit once per point: the
/// R5 Gamma-only path never reconstructed, so it could not see that; the first
/// full decomposition at a generic `t` did (identity character 12 or 8 instead of
/// 6 in the probe that motivated this helper).
fn line_folded_arms(
    arms: &[(Vec3R, Mat3I)],
    wave_vector: &Vec3R,
    embedding: &SubgroupEmbedding,
    little_dimension: u8,
) -> Result<Vec<FoldedStar>, FullStarError> {
    // Each arm's wave vector is the image of the centre wave vector under the
    // arm's own parent rotation: `R^-T k(t)`, i.e. exactly `t . arm`.  It is fed
    // to the fold **unreduced**: the parent's Bloch phases and the child's are
    // read at wave vectors that differ by (K, T^T K) together, and reducing only
    // one of the two halves would pair the frozen table with a wave vector it
    // does not describe (see [`line_wave_vector`]).
    let mut folded = Vec::with_capacity(arms.len());
    for (_, rotation) in arms {
        let action = Mat3R::from_ints(*rotation).inverse()?.transpose();
        let arm_k = action.checked_mul_vector(wave_vector)?;
        folded.push(FoldArm {
            wave_vector: arm_k,
            dimension: usize::from(little_dimension),
        });
    }
    Ok(fold_arms(embedding, embedding.parent_sg(), &folded)?)
}

/// The same frequency computed through `build_block` instead of the hand-written
/// per-arm sum.
///
/// Folds the line's arms into the child's zone with the shared [`fold_arms`], reads
/// the child Gamma block through the same `build_block` stage the discrete probes
/// use, and extracts the child's trivial row exactly as
/// [`trivial_content_with_embedding`] does.  It is kept as the R5 route so the
/// 5,756-row audit history stays reproducible, and as a second reading of
/// [`subduce_line_at_parameter`]`(...).trivial_content()`: it skips the non-Gamma
/// blocks and the reconstruction, so agreement between the two is a real (if
/// narrow) check rather than a delegation.
pub fn line_trivial_content_via_blocks(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
) -> Result<u32, FullStarError> {
    validate_line_context(subgroup, embedding, table)?;

    let direction = line_direction(table)?;
    let arms = line_arms(subgroup, embedding, table, &direction)?;
    let parameter = Rat::new(OFFICIAL_LINE_PARAMETER.0, OFFICIAL_LINE_PARAMETER.1)?;
    let wave_vector = line_wave_vector(&direction, &parameter)?;
    let arms: &[(Vec3R, Mat3I)] = &arms;
    let wave_vector = &wave_vector;
    let direction = &direction;
    let child_sg = embedding.subgroup_sg();
    let trivial = trivial_child_record(child_sg)
        .ok_or(FullStarError::MissingChildTrivialIrrep { sg: child_sg })?;
    let trivial_cir = match trivial.source_identity() {
        IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => cir_irnumber,
        IrrepSourceIdentity::Compound { .. } | IrrepSourceIdentity::Spin { .. } => {
            return Err(FullStarError::MissingChildTrivialIrrep { sg: child_sg });
        }
    };
    let child_cell = Lattice::new(exact_primitive_basis(child_sg)?)?;
    let child_reciprocal = child_cell.reciprocal()?;
    let source = LineArmSource {
        table,
        direction: *direction,
        wave_vector: *wave_vector,
        arms,
    };
    let stars = line_folded_arms(arms, wave_vector, embedding, table.dimension)?;
    let mut total = 0u32;
    for star in &stars {
        let mut is_gamma = false;
        for point in star.points() {
            if child_reciprocal.contains(point.q())? {
                is_gamma = true;
                break;
            }
        }
        if !is_gamma {
            continue;
        }
        let block = build_block(
            embedding,
            &ArmCharacterSource::Line(&source),
            star,
            &child_cell,
            &child_reciprocal,
        )?;
        for target in block.targets() {
            if target.dimension == trivial.dim && target.irnumber == Some(trivial_cir) {
                total += target.multiplicity;
            }
        }
    }
    Ok(total)
}

/// The complete subduction of one parametric-k line at **one explicit
/// parameter value** (`k = t . direction`).
///
/// This is R6's *ability A*: the whole folded-star decomposition of the parent
/// line irrep at a rational `t`, with the same invariants the discrete entry
/// point enforces (dimension conservation, integral multiplicities, per-arm
/// reconstruction of the parent character).  It is **not** a statement about a
/// parameter range: a different `t` is a different parent representation and has
/// to be asked for separately (see `docs/subduction-r6-plan.md` §1).
/// A successful result reports whether this is an actual frozen-line irrep
/// (`ParameterKind::LineIrrep`) or a formal induction at an enhanced-symmetry
/// parameter (`ParameterKind::Formal`).
///
/// Failure semantics are the discrete ones: a folded child star whose little
/// co-group is outside the constructed families and whose `q` matches no stored
/// child row is a [`FullStarError::MissingChildStarData`] error, never a zero
/// result.  A table from another parent is [`FullStarError::LineSourceMismatch`],
/// and a stale record/embedding pair fails exactly as in the discrete path.
pub fn subduce_line_at_parameter(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    parameter: Rat,
) -> Result<LineSubduction, FullStarError> {
    validate_line_context(subgroup, embedding, table)?;
    let direction = line_direction(table)?;
    let arms = line_arms(subgroup, embedding, table, &direction)?;
    let wave_vector = line_wave_vector(&direction, &parameter)?;
    let parameter_kind = line_parameter_kind(
        embedding.parent_lattice(),
        table,
        &direction,
        &wave_vector,
        &parameter,
        &arms,
    )?;
    let child_sg = embedding.subgroup_sg();
    let child_cell = Lattice::new(exact_primitive_basis(child_sg)?)?;
    let child_reciprocal = child_cell.reciprocal()?;
    let source = LineArmSource {
        table,
        direction,
        wave_vector,
        arms: &arms,
    };
    let stars = line_folded_arms(&arms, &wave_vector, embedding, table.dimension)?;
    // The parent full-star dimension of a line irrep is `little dim x arms`: the
    // frozen table's `dimension` is the little dimension and the arms are the
    // star of the line.
    let arm_count = u32::try_from(arms.len()).map_err(|_| {
        SubductionError::RationalOverflow {
            operation: "line arm count",
        }
    })?;
    let parent_dimension =
        u32::from(table.dimension)
            .checked_mul(arm_count)
            .ok_or(SubductionError::RationalOverflow {
                operation: "line full-star dimension",
            })?;

    let source = ArmCharacterSource::Line(&source);
    let mut blocks = Vec::with_capacity(stars.len());
    let mut covered = 0u32;
    for star in &stars {
        let block = build_block(embedding, &source, star, &child_cell, &child_reciprocal)?;
        covered = covered.checked_add(block.block_dimension()).ok_or(
            SubductionError::RationalOverflow {
                operation: "line full-star dimension",
            },
        )?;
        blocks.push(block);
    }
    if covered != parent_dimension {
        return Err(FullStarError::TotalDimensionMismatch {
            expected: parent_dimension,
            found: covered,
        });
    }
    let (parent_characters, reconstructed) = reconstruct(embedding, &source, &blocks)?;
    Ok(LineSubduction {
        parent_sg: subgroup.parent_sg,
        label: table.label,
        k_label: table.k_label,
        parameter,
        wave_vector,
        subgroup_sg: child_sg,
        ordinal: embedding.ordinal(),
        setting: embedding.setting(),
        setting_denominator: embedding.setting_denominator(),
        parent_dimension,
        parameter_kind,
        blocks,
        representatives: embedding.representatives().to_vec(),
        parent_characters,
        reconstructed,
        tolerance: SUBDUCTION_TOLERANCE,
    })
}

/// Classification of a parameter on a frozen parametric-k line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParameterKind {
    /// No distinct frozen line arm is reciprocal-equivalent to this wave vector;
    /// the frozen line little co-group is the actual little co-group.
    LineIrrep,
    /// Distinct frozen line arms coincide modulo the parent reciprocal lattice.
    /// The result is a formal induction from the frozen line little co-group and
    /// must not be interpreted as an irrep at the enhanced-symmetry wave vector.
    Formal,
}

/// The full decomposition of one parametric-k line irrep at one parameter
/// value; see [`subduce_line_at_parameter`].
#[derive(Debug, Clone)]
pub struct LineSubduction {
    parent_sg: u8,
    label: &'static str,
    k_label: &'static str,
    parameter: Rat,
    wave_vector: Vec3R,
    subgroup_sg: u8,
    ordinal: usize,
    setting: Mat3I,
    setting_denominator: i32,
    parent_dimension: u32,
    parameter_kind: ParameterKind,
    blocks: Vec<FullStarBlock>,
    representatives: Vec<ExactSeitz>,
    parent_characters: Vec<Complex64>,
    reconstructed: Vec<Complex64>,
    tolerance: f64,
}

impl LineSubduction {
    /// Parent space group number.
    pub const fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// Compact little-irrep label of the source, e.g. `DT1`.
    pub const fn label(&self) -> &'static str {
        self.label
    }

    /// Program label of the k domain, e.g. `DT`.
    pub const fn k_label(&self) -> &'static str {
        self.k_label
    }

    /// The exact parameter `t` this result belongs to.
    pub const fn parameter(&self) -> &Rat {
        &self.parameter
    }

    /// The raw `t . direction` in the frozen table's own frame, **unreduced**
    /// (see [`line_wave_vector`]).
    ///
    /// Two parameters one step apart report wave vectors differing by exactly the
    /// frozen direction, which is a parent reciprocal lattice vector: the same
    /// point of the parent's zone, and the reason the transported query is the
    /// monodromy image of the label rather than the label again.
    pub const fn wave_vector(&self) -> &Vec3R {
        &self.wave_vector
    }

    /// Subgroup space group number.
    pub const fn subgroup_sg(&self) -> u8 {
        self.subgroup_sg
    }

    /// Isotropy record ordinal of the embedding.
    pub const fn ordinal(&self) -> usize {
        self.ordinal
    }

    /// Numerator of the embedding's setting transform; the exact matrix is
    /// `setting() / setting_denominator()`.
    pub const fn setting(&self) -> Mat3I {
        self.setting
    }

    /// Denominator of the setting transform, always positive.
    pub const fn setting_denominator(&self) -> i32 {
        self.setting_denominator
    }

    /// Full-star dimension of the line irrep: little dimension x arms.
    pub const fn parent_dimension(&self) -> u32 {
        self.parent_dimension
    }

    /// Whether this is an actual line irrep or a formal induction at an
    /// enhanced-symmetry parameter value. `Formal` results preserve the
    /// decomposition from the frozen line little group; they are not claims
    /// about irreps of the larger little group at that wave vector.
    pub const fn parameter_kind(&self) -> ParameterKind {
        self.parameter_kind
    }

    /// Folded child-star blocks, one per child **orbit** (star), not one per
    /// folded `q`: a block carries every folded point of its orbit together with
    /// every parent arm that lands on one of them.  At a generic parameter most
    /// orbits have more than one point, so `blocks().len()` is strictly smaller
    /// than the number of folded points.
    pub fn blocks(&self) -> &[FullStarBlock] {
        &self.blocks
    }

    /// Subgroup coset representatives the reconstruction was checked on.
    pub fn representatives(&self) -> &[ExactSeitz] {
        &self.representatives
    }

    /// The parent line character on the subgroup representatives and its
    /// reconstruction from the reported child irreps (equal entries are the
    /// witness that the decomposition is complete in the parent frame).
    pub fn reconstruction(&self) -> (&[Complex64], &[Complex64]) {
        (&self.parent_characters, &self.reconstructed)
    }

    /// Tolerance used for the multiplicity, dimension and reconstruction checks.
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Total parent subspace accounted for by the reported blocks; equal to
    /// [`Self::parent_dimension`] by construction.
    pub fn covered_dimension(&self) -> u32 {
        self.blocks.iter().map(FullStarBlock::block_dimension).sum()
    }

    /// Multiplicity of the subgroup's trivial representation in this
    /// decomposition, read through the child's frozen CIR source number and
    /// cross-checked against the child row's label -- the same two readings the
    /// R5 audit compares against the pinned identity table.
    ///
    /// Errors instead of guessing: a missing or ambiguous trivial child row is
    /// [`FullStarError::MissingChildTrivialIrrep`], and two readings that
    /// disagree are a [`FullStarError::TargetSourceMismatch`].
    pub fn trivial_content(&self) -> Result<u32, FullStarError> {
        let trivial = trivial_child_record(self.subgroup_sg)
            .ok_or(FullStarError::MissingChildTrivialIrrep { sg: self.subgroup_sg })?;
        let cir = match trivial.source_identity() {
            IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => cir_irnumber,
            IrrepSourceIdentity::Compound { .. } | IrrepSourceIdentity::Spin { .. } => {
                return Err(FullStarError::MissingChildTrivialIrrep { sg: self.subgroup_sg });
            }
        };
        let mut total = 0u32;
        for block in &self.blocks {
            for target in block.targets() {
                if target.dimension == trivial.dim && target.irnumber == Some(cir) {
                    total += target.multiplicity;
                }
            }
        }
        let by_label: u32 = self
            .blocks
            .iter()
            .map(|block| block.multiplicity(trivial.ml))
            .sum();
        if total != by_label {
            return Err(FullStarError::TargetSourceMismatch {
                sg: self.subgroup_sg,
                ml: trivial.ml,
            });
        }
        Ok(total)
    }
}

/// The child table's trivial row: one-dimensional, at Gamma, and `+1` on every
/// operation of its own stored setting.
///
/// `None` when the child has no such row or more than one of them, so the
/// trivial content of a subduction onto it is not defined by the table.
fn trivial_child_record(sg: u8) -> Option<&'static IrrepRecord> {
    let mut found: Option<&'static IrrepRecord> = None;
    for record in query::irreps_of(sg) {
        if record.spinor || record.dim != 1 || record.k_vector().numerators != [0, 0, 0] {
            continue;
        }
        let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() else {
            continue;
        };
        if row.dimension() == 1
            && row
                .values()
                .iter()
                .all(|value| (value - 1.0).norm() < SUBDUCTION_TOLERANCE)
        {
            if found.is_some() {
                return None;
            }
            found = Some(record);
        }
    }
    found
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

    let source = ArmCharacterSource::Cir(star);
    let mut blocks = Vec::with_capacity(folded.len());
    let mut covered = 0u32;
    for folded_star in folded {
        let block = build_block(
            embedding,
            &source,
            folded_star,
            &child_cell,
            &child_reciprocal,
        )?;
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

    let (parent_characters, reconstructed) = reconstruct(embedding, &source, &blocks)?;
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
        setting_denominator: embedding.setting_denominator(),
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

/// One complex component of a child target.
///
/// A stored component is expanded from a pinned child record; a constructed one
/// is computed from the subgroup operations at an exact folded point and has no
/// record, no label and no CIR number at all.
#[derive(Debug, Clone)]
struct ChildComponent {
    /// Pinned child row this component was expanded from.
    record: Option<&'static IrrepRecord>,
    component: SubductionComponent,
    label: Option<&'static str>,
    bc: Option<&'static str>,
    irnumber: Option<u32>,
    dimension: u8,
    base_k: Vec3R,
    /// `base_k`, or its exact negation for a realification conjugate.
    effective_k: Vec3R,
    conjugate: bool,
    characters: ComponentCharacters,
}

/// Where one child component's little-group characters come from.
///
/// [`ConstructedLittleRep`] lives with the arm bookkeeping in
/// [`super`](crate::irrep::subduction::star), next to the induced-star character
/// source it shares with stored rows.
#[derive(Debug, Clone)]
enum ComponentCharacters {
    /// A pinned row, indexed by its own Seitz representatives.
    Stored(CharacterRow),
    /// A representation constructed from the subgroup operations and an exact
    /// folded point.
    Constructed(ConstructedLittleRep),
}

impl ChildComponent {
    /// Stable identity used to keep one occurrence of one component per orbit.
    fn key(&self) -> (Option<u16>, SubductionComponent) {
        (self.record.map(|record| record.id().index()), self.component)
    }

    /// The physical row behind this component, when it has one.
    fn row_ml(&self) -> Option<&'static str> {
        self.record.map(|record| record.ml)
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
/// The arm-character source a parent q-block is read from.
///
/// A discrete probe's star and a parametric-k line's arms answer the same two
/// questions, so `build_block` is generic over them without changing what it
/// does downstream.
enum ArmCharacterSource<'a> {
    /// A discrete probe's star.
    Cir(&'a ScalarStar),
    /// A parametric-k line's arms.
    Line(&'a LineArmSource<'a>),
}

impl ArmCharacterSource<'_> {
    fn q_block_dimension(&self, arm_indices: &[usize]) -> Result<u32, StarError> {
        match self {
            Self::Cir(star) => star.q_block_dimension(arm_indices),
            Self::Line(source) => source.q_block_dimension(arm_indices),
        }
    }

    fn q_block_character(
        &self,
        arm_indices: &[usize],
        operation: &ExactSeitz,
    ) -> Result<Complex64, StarError> {
        match self {
            Self::Cir(star) => star.q_block_character(arm_indices, operation),
            Self::Line(source) => source.q_block_character(arm_indices, operation),
        }
    }

    /// Character of the **whole** parent star on one parent operation.
    ///
    /// `build_block` reads a single q-block (the arms folding onto one child
    /// point); the independent reconstruction instead needs the parent star's
    /// value on each subgroup representative, which is the sum over every arm.
    /// Both sources answer that from the same per-arm data, so the line path
    /// gets the same reconstruction gate as the discrete one.
    fn full_character(&self, operation: &ExactSeitz) -> Result<Complex64, FullStarError> {
        match self {
            Self::Cir(star) => Ok(star.character(operation)?),
            Self::Line(source) => {
                let mut value = Complex64::new(0.0, 0.0);
                for index in 0..source.arms.len() {
                    value += source.q_block_character(&[index], operation)?;
                }
                Ok(value)
            }
        }
    }
}

fn build_block(
    embedding: &SubgroupEmbedding,
    source: &ArmCharacterSource,
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
    let q_block_dimension = source.q_block_dimension(point.arm_indices())?;
    let mut parent_characters = Vec::with_capacity(parent_operations.len());
    for operation in &parent_operations {
        parent_characters.push(source.q_block_character(point.arm_indices(), operation)?);
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

/// Does **any** point of this folded star reach a pinned child row?
///
/// A star is the reachability unit for the constructed fallback: one physical
/// child irrep lives on the whole star, so a star that already has a stored
/// component must never *also* get a constructed one for another of its arms --
/// that would present one irrep twice, with two identities the solver cannot
/// match.  Only a star with no stored component anywhere is constructed.
fn star_has_stored_components(
    child_sg: u8,
    star: &FoldedStar,
    child_reciprocal: &Lattice,
) -> Result<bool, FullStarError> {
    for point in star.points() {
        if !stored_child_components_at(child_sg, point.q(), child_reciprocal)?.is_empty() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Pinned child rows whose effective arm folds onto `q`.
/// The pinned components that live at `q`, modulo the child reciprocal lattice.
///
/// A child record the engine cannot expand (an unsupported character space, e.g.
/// a compound row with no readable selected-arm trace) is **not** skipped here:
/// it fails the whole lookup loudly.  Turning it into `continue` would report the
/// star as `MissingChildStarData`, i.e. as an archive gap, when the truth is that
/// the shipped record could not be read -- a generation bug wearing a coverage
/// gap's clothes.  The distinction is unreachable on the pinned corpus (reviewer
/// B checked all 4,105 ordinary, 672 compound and 3,611 spinor records: zero
/// unexpandable), and the loud branch is what keeps it visible if that changes.
fn stored_child_components_at(
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

/// Targets constructed from the subgroup operations and the exact folded point,
/// for the child space groups where that is enough.
///
/// The criterion is the geometry, not the child number: when every child
/// rotation that fixes `q` modulo the child reciprocal lattice is the identity,
/// the little co-group at `q` is trivial, the little group is the translation
/// group, and the one-dimensional Bloch phase `D(T_L) = exp(+2 pi i q.L)` *is*
/// the little-group irrep.  No character table and no archived source is needed,
/// so every such child is answered.  A point with a non-trivial little co-group
/// returns an empty list and keeps the existing `MissingChildStarData`
/// behaviour for the batch that supplies those sources.
///
/// The identity carries the point reduced modulo the child reciprocal lattice:
/// two folded points differing by a reciprocal lattice vector are the same
/// representation, and using both would present the same row twice.
fn constructed_child_components_at(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<Vec<ChildComponent>, FullStarError> {
    if !has_trivial_little_co_group(child_sg, q, child_reciprocal)? {
        return constructed_projective_components(child_sg, q, child_reciprocal);
    }
    let exact = child_reciprocal.reduce(q)?.representative;
    let q_key = [exact.get(0), exact.get(1), exact.get(2)];
    Ok(vec![ChildComponent {
        record: None,
        component: SubductionComponent::Constructed { q: q_key, index: 0 },
        label: None,
        bc: None,
        irnumber: None,
        dimension: 1,
        base_k: exact,
        effective_k: exact,
        conjugate: false,
        characters: ComponentCharacters::Constructed(ConstructedLittleRep::BlochPhase { q: exact }),
    }])
}

/// The one-dimensional projective catalogue of a non-trivial little co-group.
///
/// The **same** reduced point is used for the cocycle, for the constants and for
/// the evaluation: a reconstructed little-group operation is not a lattice
/// translate of its representative, so the two halves of
/// `D(R, T) = exp(2 pi i (psi_R - q.t_R + q.T))` have to share one `q`.
/// The catalogue is only returned when the solver finds exactly `|P_q|`
/// characters, i.e. when every irreducible projective representation of the
/// co-group is one-dimensional; otherwise the caller keeps
/// `MissingChildStarData`.
fn constructed_projective_components(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<Vec<ChildComponent>, FullStarError> {
    let exact = child_reciprocal.reduce(q)?.representative;
    let co_group = catalogue::little_co_group(child_sg, &exact, child_reciprocal)?;
    let characters = catalogue::one_dimensional_characters(&co_group)?;
    let q_key = [exact.get(0), exact.get(1), exact.get(2)];
    let mut components = Vec::new();
    if characters.len() == co_group.order() {
        for (index, psi) in characters.iter().enumerate() {
            let constants = catalogue::constants(&co_group, psi)?;
            components.push(constructed_component(
                q_key,
                exact,
                index,
                1,
                ConstructedLittleRep::Projective { q: exact, constants },
            )?);
        }
        return Ok(components);
    }
    // R4 batch 2b: the small higher-dimensional families (one two-dimensional
    // irrep for a non-degenerate four-element abelian co-group, the gauged
    // ordinary `D3` irreps otherwise).  The catalogue is gated by its own
    // structural and orthogonality checks and is empty for anything else.
    let targets = catalogue::projective_targets(child_sg, &exact, child_reciprocal)?;
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    for (index, target) in targets.into_iter().enumerate() {
        components.push(constructed_component(
            q_key,
            exact,
            index,
            target.dimension,
            ConstructedLittleRep::ProjectiveTable {
                q: exact,
                constants: target.constants,
            },
        )?);
    }
    Ok(components)
}

/// One constructed child component with its canonical identity.
fn constructed_component(
    q_key: [Rat; 3],
    exact: Vec3R,
    index: usize,
    dimension: u8,
    characters: ConstructedLittleRep,
) -> Result<ChildComponent, FullStarError> {
    Ok(ChildComponent {
        record: None,
        component: SubductionComponent::Constructed {
            q: q_key,
            index: u16::try_from(index).map_err(|_| SubductionError::RationalOverflow {
                operation: "constructed catalogue index",
            })?,
        },
        label: None,
        bc: None,
        irnumber: None,
        dimension,
        base_k: exact,
        effective_k: exact,
        conjugate: false,
        characters: ComponentCharacters::Constructed(characters),
    })
}

/// Whether the little co-group of `q` is trivial in `child_sg`: no non-identity
/// child rotation of the child's own data-Hall operations fixes `q` modulo the
/// child reciprocal lattice (centring extinctions included).
fn has_trivial_little_co_group(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<bool, FullStarError> {
    for operation in &strict_sg_hall_ops(child_sg)?.operations {
        let rotation = operation.rotation();
        if rotation != IDENTITY_ROTATION && child_reciprocal.preserves(rotation, q)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether a folded child star with no pinned child row is answerable by a
/// **constructed** target at its canonical first point.
///
/// This is the single production predicate behind the constructed fallback:
/// R4 batch 1 answers a trivial little co-group with the one-dimensional Bloch
/// phase, R4 batch 2a answers a non-trivial one whose projective irreps are all
/// one-dimensional (complete exactly when the solver finds `|P_q|` characters).
/// Everything else needs the higher-dimensional batch and is not answerable yet.
///
/// The offline census uses it so its per-star reachability model cannot drift
/// from what the entry point actually does.
pub fn constructed_targets_available(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<bool, FullStarError> {
    Ok(!constructed_child_components_at(child_sg, q, child_reciprocal)?.is_empty())
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
                record: Some(record),
                component: SubductionComponent::Ordinary,
                label: Some(record.ml),
                bc: Some(record.bc),
                irnumber: Some(irnumber),
                dimension,
                base_k,
                effective_k: base_k,
                conjugate: false,
                characters: ComponentCharacters::Stored(row),
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
        record: Some(record),
        component,
        label: Some(constituent.label),
        bc: Some(record.bc),
        irnumber: Some(constituent.irnumber),
        dimension,
        base_k,
        effective_k,
        conjugate,
        characters: ComponentCharacters::Stored(constituent.row),
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
///
/// The constructed fallback is chosen **per star**, not per point: a star that
/// reaches any pinned row is answered from pinned rows alone (a constructed
/// component would be the same physical irrep under a second identity), and only
/// a star with no pinned row anywhere is constructed -- once, at the star's
/// canonical first point, because a trivial little co-group carries exactly one
/// little-group irrep on the whole star.
fn select_representative(
    child_sg: u8,
    folded_star: &FoldedStar,
    child_reciprocal: &Lattice,
) -> Result<Representative, FullStarError> {
    let constructed = !star_has_stored_components(child_sg, folded_star, child_reciprocal)?;
    let mut representative: Option<(usize, Vec3R)> = None;
    let mut components: Vec<ChildComponent> = Vec::new();
    if constructed {
        if let Some(point) = folded_star.points().first() {
            components = constructed_child_components_at(child_sg, point.q(), child_reciprocal)?;
            if let Some(first) = components.first() {
                representative = Some((0, first.effective_k));
            }
        }
    } else {
        for (index, point) in folded_star.points().iter().enumerate() {
            let at_point = stored_child_components_at(child_sg, point.q(), child_reciprocal)?;
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
        let (SubductionComponent::RealificationSeed { .. }, Some(record)) =
            (component.component, component.record)
        else {
            continue;
        };
        let partner = representative.components.iter().position(|other| {
            other.record.is_some_and(|other| other.id() == record.id())
                && matches!(
                    other.component,
                    SubductionComponent::RealificationConjugate { .. }
                )
        });
        let Some(partner) = partner else {
            continue;
        };
        match realification_relation(child_sg, record.ml, &values[index], &values[partner])? {
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
            row_ml: component.row_ml(),
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
        let value = match &component.characters {
            ComponentCharacters::Stored(row) => character_of(
                row,
                component.label.unwrap_or("stored component"),
                &conjugated,
                child_cell,
                &component.base_k,
                index,
            )?,
            ComponentCharacters::Constructed(rep) => rep.character(&conjugated)?,
        };
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
            && component.row_ml() == target.row_ml
            && component.dimension == target.dimension
    });
    found.ok_or(FullStarError::TargetSourceMismatch {
        sg: child_sg,
        ml: target
            .ml
            .unwrap_or("<constructed target with no matching component>"),
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
    match (&component.component, &component.characters) {
        (SubductionComponent::Ordinary, ComponentCharacters::Stored(_)) => {
            let record = component.record.ok_or(FullStarError::TargetSourceMismatch {
                sg: child_sg,
                ml: "stored ordinary component without a row",
            })?;
            Ok(ChildStarEvaluator::Ordinary(Box::new(OrdinaryStar::new(record)?)))
        }
        (
            SubductionComponent::Constituent { .. }
            | SubductionComponent::RealificationSeed { .. }
            | SubductionComponent::RealificationConjugate { .. },
            ComponentCharacters::Stored(row),
        ) => Ok(ChildStarEvaluator::Component(Box::new(ComponentStar::new(
            child_sg,
            component.label.unwrap_or("stored component"),
            component.component,
            component.irnumber.unwrap_or(0),
            component.base_k,
            row.clone(),
            component.conjugate,
            None,
        )?))),
        (SubductionComponent::Constructed { .. }, ComponentCharacters::Constructed(rep)) => {
            Ok(ChildStarEvaluator::Constructed(Box::new(
                ConstructedStar::new(child_sg, component.base_k, rep.clone(), component.dimension)?,
            )))
        }
        _ => Err(FullStarError::TargetSourceMismatch {
            sg: child_sg,
            ml: "component and character source disagree",
        }),
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
    source: &ArmCharacterSource,
    blocks: &[FullStarBlock],
) -> Result<(Vec<Complex64>, Vec<Complex64>), FullStarError> {
    let representatives = embedding.representatives();
    let mut parent_characters = Vec::with_capacity(representatives.len());
    let mut pulled_back = Vec::with_capacity(representatives.len());
    for operation in representatives {
        parent_characters.push(source.full_character(operation)?);
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
    use crate::irrep::LabelConvention;
    use crate::irrep::line_monodromy;
    use crate::irrep::subduction::star::FoldedPoint;
    use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
    use crate::irrep::subduce_irrep;
    use crate::irrep::subduction::subduce_irrep_with_embedding;

    fn probe(sg: u8, ml: &str) -> &'static IrrepRecord {
        query::irreps_of(sg)
            .iter()
            .find(|record| record.ml == ml && !record.spinor)
            .unwrap_or_else(|| panic!("SG {sg} has no scalar irrep {ml}"))
    }

    fn embedding(parent: u8, ml: &str, direction: &str) -> SubgroupEmbedding {
        let subgroup = isotropy_subgroup_for_direction(
            parent,
            ml,
            LabelConvention::Cdml,
            IsotropyDirection::Label(direction),
        )
        .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction}: {error}"));
        SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
            .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction} embedding: {error}"))
    }

    /// Every space group of the pinned table has exactly one trivial Gamma row,
    /// so the `MissingChildTrivialIrrep` guard is unreachable data-wise and can
    /// never turn a real subduction into a silent zero.
    #[test]
    fn every_space_group_has_one_trivial_gamma_row() {
        for sg in 1u8..=230 {
            let found = trivial_child_record(sg)
                .unwrap_or_else(|| panic!("space group {sg} has no trivial Gamma row"));
            assert_eq!(found.dim, 1, "space group {sg}");
            assert_eq!(found.k_vector().numerators, [0, 0, 0], "space group {sg}");
            assert!(!found.spinor, "space group {sg}");
        }
    }

    /// The identity-only entry point reproduces the full decomposition's trivial
    /// content whenever the full decomposition has data, and answers exactly
    /// (with every non-Gamma star skipped) when it does not.
    #[test]
    fn trivial_content_skips_only_stars_that_cannot_carry_the_trivial_rep() {
        // SG 221 GM4+ P1 -> #83 has a full decomposition for every probe.
        let subgroup = subgroup_of(221, "GM4+", "P1");
        let built = embedding(221, "GM4+", "P1");
        let trivial = trivial_child_record(built.subgroup_sg()).unwrap();
        for record in query::irreps_of(221)
            .iter()
            .filter(|record| !record.spinor)
        {
            let content = trivial_content_with_embedding(&subgroup, &built, record)
                .unwrap_or_else(|error| panic!("SG 221 GM4+ P1 {}: {error}", record.ml));
            let full = subduce_full_star_with_embedding(&subgroup, &built, record).unwrap();
            let full_total: u32 = full
                .blocks()
                .iter()
                .map(|block| block.multiplicity(trivial.ml))
                .sum();
            assert_eq!(content.total, full_total, "probe {}", record.ml);
            assert_eq!(content.total, content.by_label, "probe {}", record.ml);
            let folded = ScalarStar::new(record).unwrap().folded_stars(&built).unwrap();
            assert_eq!(
                content.gamma_stars + content.skipped_stars,
                folded.len(),
                "probe {}",
                record.ml
            );
        }

        // SG 196 W1 -> #24 folds onto child stars the #24 table has no k for.
        // R4 batch 2a constructs them from the exact cocycle, so the full entry
        // now agrees with the identity-only entry on the trivial content.
        let subgroup = subgroup_of(196, "W1", "P2");
        let built = embedding(196, "W1", "P2");
        let full = subduce_full_star_with_embedding(&subgroup, &built, probe(196, "W1"))
            .expect("the one-dimensional catalogue answers these stars");
        let content = trivial_content_with_embedding(&subgroup, &built, probe(196, "W1")).unwrap();
        assert_eq!(
            (content.total, content.by_label, content.gamma_stars, content.skipped_stars),
            (1, 1, 1, 2)
        );
        let full_total: u32 = full
            .blocks()
            .iter()
            .map(|block| block.multiplicity(trivial_child_record(24).unwrap().ml))
            .sum();
        assert_eq!(full_total, content.total);
    }

    fn rat(num: i128, den: i128) -> Rat {
        Rat::new(num, den).expect("a rational with a non-zero denominator")
    }

    fn subgroup_of(parent: u8, ml: &str, direction: &str) -> IsotropySubgroup {
        isotropy_subgroup_for_direction(
            parent,
            ml,
            LabelConvention::Cdml,
            IsotropyDirection::Label(direction),
        )
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
                    target.ml.expect("stored target ML label"),
                    target.dimension,
                    target.multiplicity,
                    target.component,
                    target.irnumber.expect("stored target CIR number"),
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
                    .find(|record| record.ml == target.ml.expect("stored target ML label") && record.compound_metadata().is_none())
                    .unwrap_or_else(|| {
                        panic!("SG {child_sg} has no ordinary {ml}", ml = target.ml.expect("stored target ML label"))
                    });
                assert!(matches!(
                    record.source_identity(),
                    IrrepSourceIdentity::OrdinaryScalar { cir_irnumber }
                        if cir_irnumber == target.irnumber.expect("stored target CIR number")
                ));
                let row = record
                    .ordinary_scalar_selected_arm_block_trace()
                    .expect("ordinary row");
                assert_eq!(u32::from(target.dimension), row.dimension() as u32);
                assert_eq!(target.ml.expect("stored target ML label"), target.row_ml.expect("stored target row"));
            }
            SubductionComponent::Constituent { index, irnumber } => {
                assert_eq!(target.irnumber.expect("stored target CIR number"), irnumber);
                let record = query::irreps_of(child_sg)
                    .iter()
                    .find(|record| record.ml == target.row_ml.expect("stored target row"))
                    .unwrap_or_else(|| panic!("SG {child_sg} has no {ml}", ml = target.row_ml.expect("stored target row")));
                let metadata = record.compound_metadata().expect("compound metadata");
                assert_eq!(metadata.cir_irnumbers[usize::from(index)], irnumber);
                assert_eq!(metadata.cir_labels[usize::from(index)], target.ml.expect("stored target ML label"));
                assert_eq!(
                    u32::from(metadata.cir_dimensions[usize::from(index)]),
                    u32::from(target.dimension)
                );
                let view = record.compound_selected_arm_view().expect("compound view");
                let CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } = view
                else {
                    panic!("{ml} is not a distinct component sum", ml = target.row_ml.expect("stored target row"));
                };
                let stored = if index == 0 { first } else { second };
                assert_eq!(stored.irnumber, irnumber);
                assert_eq!(stored.label, target.ml.expect("stored target ML label"));
                assert_eq!(stored.dimension as u32, u32::from(target.dimension));
            }
            SubductionComponent::Constructed { .. } => {
                panic!("a constructed target cannot appear in a pinned stored fixture")
            }
            SubductionComponent::RealificationSeed { irnumber }
            | SubductionComponent::RealificationConjugate { irnumber } => {
                // Both components share the single stored CIR seed of the row.
                let record = query::irreps_of(child_sg)
                    .iter()
                    .find(|record| record.ml == target.row_ml.expect("stored target row"))
                    .unwrap_or_else(|| panic!("SG {child_sg} has no {ml}", ml = target.row_ml.expect("stored target row")));
                let metadata = record.compound_metadata().expect("compound metadata");
                assert_eq!(
                    metadata.semantics,
                    CompoundCharacterSemantics::ConjugateRealification
                );
                assert_eq!(metadata.cir_irnumbers[0], irnumber);
                assert_eq!(metadata.cir_labels[0], target.ml.expect("stored target ML label"));
                assert_eq!(
                    u32::from(metadata.cir_dimensions[0]),
                    u32::from(target.dimension)
                );
                let view = record.compound_selected_arm_view().expect("compound view");
                let CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } = view else {
                    panic!("{ml} is not a realification", ml = target.row_ml.expect("stored target row"));
                };
                assert_eq!(seed.irnumber, irnumber);
                assert_eq!(seed.label, target.ml.expect("stored target ML label"));
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
            gm3.blocks()[0].constructed_multiplicity(SubductionComponent::Ordinary),
            0,
            "a stored component must never match a constructed-target query"
        );
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
        assert_eq!(gm4.blocks()[0].targets()[1].row_ml, Some("GM3+GM4+"));
        assert_eq!(gm4.blocks()[0].targets()[2].row_ml, Some("GM3+GM4+"));
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
        assert_eq!(result.blocks()[1].targets()[0].row_ml, Some("Z3+Z4+"));
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
                        target.ml.expect("stored target ML label"),
                        target.irnumber.expect("stored target CIR number"),
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
        // The star as a whole does reach a pinned row, so the constructed
        // fallback must stay out of it.
        let folded = ScalarStar::new(probe(167, "F1+"))
            .unwrap()
            .folded_stars(&built)
            .unwrap();
        let two_point = folded.iter().find(|star| star.star_size() == 2).unwrap();
        assert!(
            stored_child_components_at(15, two_point.points()[0].q(), &child_reciprocal)
                .unwrap()
                .is_empty()
        );
        assert!(star_has_stored_components(15, two_point, &child_reciprocal).unwrap());
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
                                target.ml.expect("stored target ML label"),
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
                                target.ml.expect("stored target ML label"),
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
            assert_eq!(block.targets()[0].row_ml, Some("W1W1"));
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
        assert_eq!(r1.blocks()[0].targets()[0].row_ml, Some("R1R1"));
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
            assert_eq!(w.blocks()[0].targets()[0].row_ml, Some(probe_ml));
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

        // A folded point that matches no stored child record is answered by the
        // constructed fallback whenever the little co-group is in scope: ordinal
        // 3988 (SG 109 `GM3` `P1` -> #43) folds onto a four-element co-group
        // whose projective table R4 batch 2b supplies, so the whole probe
        // decomposes.  The out-of-scope case (a cubic Gamma point, order 48) is
        // pinned by `constructed_targets_have_their_own_identity_and_no_borrowed_labels`.
        let two_dimensional = subgroup_of(109, "GM3", "P1");
        let built_43 = embedding(109, "GM3", "P1");
        assert!(subduce_full_star_with_embedding(&two_dimensional, &built_43, probe(109, "P1"))
            .is_ok());
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
            .filter(|target| target.row_ml.expect("stored target row") == "GM3+GM4+")
            .collect();
        assert_eq!(pair.len(), 2, "the two distinct constituents stay separate");
        assert_eq!(pair[0].irnumber.expect("stored target CIR number"), 4077);
        assert_eq!(pair[1].irnumber.expect("stored target CIR number"), 4078);
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

    // ── Constructed targets (R1/R2) ─────────────────────────────────────────

    /// The constructed P1 target is the Bloch phase with the sign convention the
    /// stored rows already use: `D(T_L) = exp(+2 pi i q.L)`.  Hand-computed
    /// values, including a non-zero translation and rational wave-vector
    /// components, pin the sign and the arithmetic without any table.
    #[test]
    fn constructed_bloch_phase_pins_the_positive_sign_convention() {
        let q = Vec3R::new([rat(1, 4), rat(-1, 2), rat(2, 3)]);
        let rep = ConstructedLittleRep::BlochPhase { q };
        let cases = [
            ([0, 0, 0], 0.0),
            ([1, 0, 0], 1.0 / 4.0),
            ([0, 1, 0], -1.0 / 2.0),
            ([0, 0, 1], 2.0 / 3.0),
            ([1, -2, 3], 1.0 / 4.0 + 1.0 + 2.0),
        ];
        for (translation, expected_turns) in cases {
            let operation = ExactSeitz::new(
                IDENTITY_ROTATION,
                Vec3R::new(translation.map(|value| rat(i128::from(value), 1))),
            );
            let found = rep.character(&operation).expect("translation character");
            let expected =
                Complex64::from_polar(1.0, std::f64::consts::TAU * expected_turns);
            assert!(
                (found - expected).norm() < SUBDUCTION_TOLERANCE,
                "translation {translation:?}: {found} != {expected}"
            );
        }
        // The conjugate point gives the conjugate phase, never the same one.
        let conjugated = ConstructedLittleRep::BlochPhase {
            q: q.checked_neg().expect("negation"),
        };
        let operation = ExactSeitz::new(
            IDENTITY_ROTATION,
            Vec3R::new([rat(1, 1), rat(0, 1), rat(0, 1)]),
        );
        assert_eq!(
            conjugated.character(&operation).expect("phase"),
            rep.character(&operation).expect("phase").conj()
        );
        // A rotation outside the class the rep was built from is an error, never
        // a silently wrong phase.
        let rotation = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let operation = ExactSeitz::new(rotation, Vec3R::new([rat(0, 1); 3]));
        assert!(matches!(
            rep.character(&operation),
            Err(StarError::ConstructedRotationNotCovered { .. })
        ));
    }

    /// A fractional frozen setting has to reach the result object **with its
    /// denominator**: ordinal 26 is `U = [[1,2,1],[-1,2,-1],[-1,0,1]] / 2`, and
    /// a caller that only received the numerator would transform coordinates
    /// with the wrong matrix.  Both result types carry `(setting, denominator)`
    /// now; this pins the real denominator-2 record rather than a synthetic one.
    #[test]
    fn a_fractional_setting_reaches_the_results_with_its_denominator() {
        let built = |parent: u8, ml: &str, direction: &str| embedding(parent, ml, direction);
        let mut found = None;
        for irrep in query::irreps_of(3).iter().filter(|irrep| !irrep.spinor) {
            let Ok(subgroups) =
                crate::irrep::isotropy::isotropy_subgroups(3, irrep.ml, LabelConvention::Cdml)
            else {
                continue;
            };
            for subgroup in subgroups {
                if subgroup.ordinal == 26 {
                    found = Some((subgroup, irrep));
                }
            }
        }
        let (subgroup, probe) = found.expect("ordinal 26 is owned by an SG 3 irrep");
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).expect("embedding");
        assert_eq!(embedding.setting_denominator(), 2);
        assert_eq!(embedding.setting(), [[1, 2, 1], [-1, 2, -1], [-1, 0, 1]]);
        let full = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
            .expect("the fractional record decomposes");
        assert_eq!(full.setting(), embedding.setting());
        assert_eq!(full.setting_denominator(), embedding.setting_denominator());
        // No Gamma-owned record carries a fractional setting on this corpus
        // (scanned over all 230 parents when this test was written), so the
        // compact result type is pinned on a Gamma context with denominator one
        // plus the same accessor wiring.
        let gamma = subgroup_of(221, "GM4+", "P1");
        let gamma_embedding = built(221, "GM4+", "P1");
        assert_eq!(gamma_embedding.setting_denominator(), 1);
        let compact = subduce_irrep(&gamma, "GM3+").expect("golden Gamma case");
        assert_eq!(compact.setting(), gamma_embedding.setting());
        assert_eq!(
            compact.setting_denominator(),
            gamma_embedding.setting_denominator()
        );
    }

    /// Every ordinary isotropy record of the pinned table, keyed by ordinal.
    fn subgroups() -> std::collections::BTreeMap<usize, IsotropySubgroup> {
        let mut table = std::collections::BTreeMap::new();
        for sg in 1..=230u8 {
            for record in query::irreps_of(sg) {
                if record.spinor || record.subgroups().is_empty() {
                    continue;
                }
                let Ok(records) =
                    crate::irrep::isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
                else {
                    continue;
                };
                for subgroup in records {
                    table.insert(subgroup.ordinal, subgroup);
                }
            }
        }
        table
    }

    // ── R6.1: the parametric-k line at one explicit parameter value ──────────

    /// The frozen little-character table of one parent and source label.
    fn line_table(parent: u8, label: &str) -> &'static LittleCharacterTable {
        crate::irrep::w_little_characters_data::W_LITTLE_CHARACTERS
            .iter()
            .find(|table| {
                usize::from(table.space_group) == usize::from(parent) && table.label == label
            })
            .unwrap_or_else(|| panic!("SG {parent} has no frozen line source {label}"))
    }

    /// The number of parent arms whose `t . arm` folds onto the child's Gamma
    /// point, times the little dimension: the closed form the trivial content of
    /// a **trivial-little-co-group child** must equal.
    ///
    /// This deliberately does not touch the character/multiplicity machinery: it
    /// folds each arm and counts, which is the geometric statement.  (The child
    /// must have a trivial little co-group at every folded point for the closed
    /// form to apply, which is why it is only used with child #1.)
    ///
    /// Scope, so the comparison is not over-read: only the **multiplicity
    /// solver** (`build_block` and its character inner product) is replaced.  The
    /// arm set, the frame and the Gamma test come from the same helpers
    /// (`line_direction`, `line_arms`, `fold_wave_vector`, `Lattice::contains`)
    /// the entry point uses, so a frame or arm-set error is absorbed on both
    /// sides.  It is a second reading of the multiplicity stage, not an
    /// independent oracle; only the pinned `t = 1/4` frequencies (and the audit
    /// over all 5,756 rows) are external.
    fn arms_folding_to_child_gamma(
        subgroup: &IsotropySubgroup,
        embedding: &SubgroupEmbedding,
        table: &'static LittleCharacterTable,
        parameter: &Rat,
    ) -> u32 {
        let direction = line_direction(table).expect("direction");
        let arms = line_arms(subgroup, embedding, table, &direction).expect("arms");
        let child_cell = Lattice::new(exact_primitive_basis(embedding.subgroup_sg()).unwrap())
            .expect("child cell");
        let child_reciprocal = child_cell.reciprocal().expect("child reciprocal");
        let mut count = 0u32;
        for (arm, _) in &arms {
            let mut scaled = [Rat::ZERO; 3];
            for (axis, value) in scaled.iter_mut().enumerate() {
                *value = parameter.checked_mul(arm.get(axis)).expect("scaled arm");
            }
            let q = crate::irrep::subduction::fold_wave_vector(
                embedding.transform(),
                &Vec3R::new(scaled),
            )
            .expect("folded arm");
            if child_reciprocal.contains(&q).expect("gamma test") {
                count += 1;
            }
        }
        count * u32::from(table.dimension)
    }

    /// One target per block as `(multiplicity, dimension, label, source)`, in the
    /// block order, so two decompositions can be compared structurally.
    fn target_rows(result: &LineSubduction) -> Vec<(u32, u8, Option<&'static str>, Option<u32>)> {
        result
            .blocks()
            .iter()
            .flat_map(|block| {
                block
                    .targets()
                    .iter()
                    .map(|target| {
                        (target.multiplicity, target.dimension, target.ml, target.irnumber)
                    })
            })
            .collect()
    }

    /// The official parameter, read through the public accessor so that moving
    /// the constant is felt by these tests instead of being pinned twice.
    fn official() -> Rat {
        official_line_parameter().expect("the official parameter is a valid rational")
    }

    /// Re-assert the two invariants the constructor already enforces before it
    /// returns `Ok` (dimension conservation and per-arm reconstruction).  This is
    /// a **tripwire for future refactors that would drop those checks**, not
    /// independent evidence: it re-runs the same predicates the same tolerance
    /// gates, so it cannot fail for a result that was built successfully.
    fn assert_line_invariants(result: &LineSubduction) {
        assert_eq!(
            result.covered_dimension(),
            result.parent_dimension(),
            "the blocks must cover the line's full-star dimension"
        );
        let (parent, reconstructed) = result.reconstruction();
        assert_eq!(parent.len(), reconstructed.len());
        for (found, expected) in reconstructed.iter().zip(parent) {
            assert!(
                (found - expected).norm() <= result.tolerance(),
                "reconstruction {found} != {expected}"
            );
        }
    }

    /// R6.1 known-point regression: at the official parameter `t = 1/4` the
    /// **complete** decomposition of every pinned SG 196 line row carries the
    /// pinned identity frequency (the pinned rows come from the official
    /// program, so this is an external comparison, not internal consistency).
    #[test]
    fn the_official_line_parameter_reproduces_the_pinned_frequencies_of_sg_196() {
        let contexts = subgroups();
        let mut rows = 0usize;
        for (ordinal, subgroup) in &contexts {
            if subgroup.parent_sg != 196 {
                continue;
            }
            let Ok(pinned) = subgroup.other_wave_vector_subduction() else {
                continue;
            };
            if pinned.is_empty() {
                continue;
            }
            let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
            for row in &pinned {
                let table = line_table(196, row.parent_ml);
                let parameter = official();
                let result = subduce_line_at_parameter(subgroup, &embedding, table, parameter)
                    .unwrap_or_else(|error| {
                        panic!("ordinal {ordinal} {} t=1/4: {error}", row.parent_ml)
                    });
                assert_line_invariants(&result);
                let trivial = result.trivial_content().expect("trivial content");
                assert_eq!(
                    trivial,
                    u32::from(row.frequency),
                    "ordinal {ordinal} {} t=1/4",
                    row.parent_ml
                );
                rows += 1;
            }
        }
        // SG 196's pinned line rows; pinned so a shrinking sweep fails.
        assert_eq!(rows, 106, "SG 196 pinned line rows");
    }

    /// A generic parameter is answered by the constructed Bloch phases: no stored
    /// child row is used, every folded q is its own one-dimensional child irrep,
    /// and the trivial content is the **geometric** count of arms folding onto the
    /// child's Gamma point (zero here, so this particular assertion is `0 == 0`;
    /// the nonzero anchor for the same closed form is the pinned `t = 1/4`
    /// comparison in the tests around it and in the full audit).
    #[test]
    fn a_generic_line_parameter_decomposes_into_constructed_targets() {
        let subgroup = subgroup_of(196, "W1", "4D1");
        let embedding = embedding(196, "W1", "4D1");
        assert_eq!(embedding.subgroup_sg(), 1);
        let table = line_table(196, "DT1");
        let parameter = Rat::new(1, 7).unwrap();
        let result = subduce_line_at_parameter(&subgroup, &embedding, table, parameter)
            .expect("a generic parameter decomposes");
        assert_line_invariants(&result);
        assert_eq!(result.parent_dimension(), 6);
        // Child #1 has a trivial point group, so every orbit is a singleton and
        // the block count still equals the number of folded points here.
        assert_eq!(result.blocks().len(), 6, "one block per child orbit");
        assert!(
            result
                .blocks()
                .iter()
                .all(|block| block.targets().iter().all(|target| matches!(
                    target.component,
                    SubductionComponent::Constructed { .. }
                ))),
            "a generic point has no stored child row"
        );
        let hand = arms_folding_to_child_gamma(&subgroup, &embedding, table, &parameter);
        assert_eq!(hand, 0, "no arm folds onto the child Gamma at t = 1/7");
        assert_eq!(result.trivial_content().unwrap(), hand);
    }

    /// At enhanced-symmetry values the frozen line representation is only a
    /// formal induction. The type tag comes from exact reciprocal equivalence
    /// of distinct line arms, not from special-casing parameter text.
    #[test]
    fn enhanced_line_parameters_are_marked_formal() {
        let subgroup = subgroup_of(196, "W1", "4D1");
        let embedding = embedding(196, "W1", "4D1");
        let table = line_table(196, "DT1");

        for parameter in [Rat::ZERO, Rat::new(1, 2).unwrap()] {
            let result = subduce_line_at_parameter(&subgroup, &embedding, table, parameter)
                .unwrap_or_else(|error| panic!("t = {parameter} should decompose: {error}"));
            assert_eq!(
                result.parameter_kind(),
                ParameterKind::Formal,
                "t = {parameter} has an enhanced parent little group"
            );
            assert_line_invariants(&result);
        }

        for parameter in [official(), Rat::new(1, 6).unwrap(), Rat::new(1, 3).unwrap()] {
            let result = subduce_line_at_parameter(&subgroup, &embedding, table, parameter)
                .unwrap_or_else(|error| panic!("t = {parameter} should decompose: {error}"));
            assert_eq!(
                result.parameter_kind(),
                ParameterKind::LineIrrep,
                "t = {parameter} keeps the frozen line little group"
            );
        }
    }

    /// The special value `t = 1/4` and its two neighbours: the decomposition is
    /// complete on both sides, but only the special value carries the pinned
    /// frequency, and the trivial content equals the geometric arm count in
    /// all three cases (the `t = 1/4` one is nonzero and externally pinned; the
    /// two neighbours are `0 == 0`).
    #[test]
    fn a_special_line_parameter_differs_from_its_neighbours() {
        let subgroup = subgroup_of(196, "W1", "4D1");
        let embedding = embedding(196, "W1", "4D1");
        let table = line_table(196, "DT1");
        let pinned = 4u32;

        let special = official();
        let result = subduce_line_at_parameter(&subgroup, &embedding, table, special)
            .expect("the official parameter decomposes");
        assert_line_invariants(&result);
        assert_eq!(result.trivial_content().unwrap(), pinned);
        assert_eq!(
            arms_folding_to_child_gamma(&subgroup, &embedding, table, &special),
            pinned,
            "the hand count must agree with the pinned frequency"
        );
        // Stored child rows answer the special point: `Z1` x2 and `GM1` x4.
        let rows = target_rows(&result);
        assert_eq!(rows.len(), 2, "{rows:?}");
        assert!(rows.iter().any(|row| row.2 == Some("Z1") && row.0 == 2));
        assert!(rows.iter().any(|row| row.2 == Some("GM1") && row.0 == 4));

        for (numerator, denominator) in [(1i128, 6i128), (1, 3)] {
            let parameter = Rat::new(numerator, denominator).unwrap();
            let result =
                subduce_line_at_parameter(&subgroup, &embedding, table, parameter)
                    .expect("both sides of the special value decompose");
            assert_line_invariants(&result);
            let hand = arms_folding_to_child_gamma(&subgroup, &embedding, table, &parameter);
            assert_eq!(hand, 0, "t = {parameter} is generic");
            assert_eq!(
                result.trivial_content().unwrap(),
                hand,
                "t = {parameter} must not reproduce the pinned frequency"
            );
        }
    }

    /// A `1/2` parameter shift that keeps the star's arm set is a case-level
    /// equivalence, not a general rule: `1/4`, `3/4` and `5/4` give the same
    /// decomposition for this source.  Characterising the general shift belongs
    /// to R6.2, so this test pins the observed case and nothing more.
    #[test]
    fn a_half_shift_of_the_line_parameter_can_leave_the_decomposition_unchanged() {
        let subgroup = subgroup_of(196, "W1", "4D1");
        let embedding = embedding(196, "W1", "4D1");
        let table = line_table(196, "DT1");
        let reference = subduce_line_at_parameter(
            &subgroup,
            &embedding,
            table,
            official(),
        )
        .expect("t = 1/4");
        for (numerator, denominator) in [(3i128, 4i128), (5, 4)] {
            let result = subduce_line_at_parameter(
                &subgroup,
                &embedding,
                table,
                Rat::new(numerator, denominator).unwrap(),
            )
            .expect("shifted parameter");
            assert_line_invariants(&result);
            assert_eq!(target_rows(&result), target_rows(&reference));
            assert_eq!(
                result.trivial_content().unwrap(),
                reference.trivial_content().unwrap()
            );
        }
    }

    /// Failure semantics: a frozen table of another parent is refused, and a
    /// folded point outside the constructed families with no stored row is a
    /// missing-data error, never a zero multiplicity.
    #[test]
    fn a_line_refuses_foreign_sources_and_out_of_scope_points() {
        let subgroup = subgroup_of(196, "W1", "4D1");
        let embedding = embedding(196, "W1", "4D1");
        let foreign = line_table(202, "DT1");
        assert!(matches!(
            subduce_line_at_parameter(
                &subgroup,
                &embedding,
                foreign,
                official()
            ),
            Err(FullStarError::LineSourceMismatch { .. })
        ));

        // Ordinal 13543 (SG 225 `W5` -> child #136): at a generic parameter two
        // folded points of one child star have no stored row and their little
        // co-group is outside the constructed families.
        let contexts = subgroups();
        let subgroup = contexts
            .values()
            .find(|subgroup| subgroup.ordinal == 13543)
            .expect("ordinal 13543");
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
        assert_eq!(embedding.subgroup_sg(), 136);
        let table = line_table(subgroup.parent_sg, "DT5");
        assert!(matches!(
            subduce_line_at_parameter(
                subgroup,
                &embedding,
                table,
                Rat::new(1, 7).unwrap()
            ),
            Err(FullStarError::MissingChildStarData { sg: 136, .. })
        ));
        // The same source at the official parameter does have the data, so the
        // error above is about the parameter, not about the record.
        let official = subduce_line_at_parameter(
            subgroup,
            &embedding,
            table,
            official(),
        )
        .expect("the official parameter is covered");
        assert_line_invariants(&official);
    }

    /// R6.1 regression for the orbit-shaped folding: at a generic parameter the
    /// folded child stars really are **orbits** with more than one point, and the
    /// block tiling stays exact.
    ///
    /// Before the orbit fix each reduced `q` was its own star, which the
    /// Gamma-only R5 path never revealed; these contexts were the witnesses
    /// (reviewer D: replacing the fix with the old grouping left every other test
    /// and the whole audit green).
    #[test]
    fn a_generic_parameter_folds_into_multi_point_child_stars() {
        // SG 196 `W1` -> child #18 `DT1` (ordinal 10030): three folded stars of
        // two points each.
        let contexts = subgroups();
        let subgroup = contexts
            .values()
            .find(|subgroup| subgroup.ordinal == 10030)
            .expect("ordinal 10030");
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
        assert_eq!(embedding.subgroup_sg(), 18);
        let result = subduce_line_at_parameter(
            subgroup,
            &embedding,
            line_table(196, "DT1"),
            Rat::new(1, 7).unwrap(),
        )
        .expect("the generic point decomposes");
        assert_line_invariants(&result);
        assert_eq!(result.blocks().len(), 3, "three folded orbits");
        assert!(
            result
                .blocks()
                .iter()
                .all(|block| block.star_size() == 2 && block.arm_count() == 2),
            "{:?}",
            result
                .blocks()
                .iter()
                .map(|block| (block.star_size(), block.arm_count(), block.block_dimension()))
                .collect::<Vec<_>>()
        );
        assert_eq!(result.parent_dimension(), 6);
        assert_eq!(result.trivial_content().unwrap(), 0);

        // SG 225 `W5` -> child #136 `SM1`: orbits of four and eight points.
        let subgroup = contexts
            .values()
            .find(|subgroup| subgroup.ordinal == 13543)
            .expect("ordinal 13543");
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
        let result = subduce_line_at_parameter(
            subgroup,
            &embedding,
            line_table(subgroup.parent_sg, "SM1"),
            Rat::new(1, 7).unwrap(),
        )
        .expect("the generic point decomposes");
        assert_line_invariants(&result);
        let mut stars: Vec<usize> = result.blocks().iter().map(FullStarBlock::star_size).collect();
        stars.sort_unstable();
        assert_eq!(stars, vec![4, 8], "multi-point orbits");
    }

    /// A parameter step by a parent reciprocal lattice vector is the
    /// **monodromy image** of the label, not a gauge of the label itself:
    ///
    /// ```text
    /// decompose(alpha, t + n)  ==  decompose(M_{n v}(alpha), t)
    /// ```
    ///
    /// The 40 contexts below are the rows of reviewer D's P1 (DT rows of SG
    /// 210/227/228) where the two readings differ.  R6.1 read them as a gauge
    /// slip and canonicalized the wave vector, which answers `Z1` for SG 210
    /// `DT3` at *both* parameters; the correct transport answers `Z2`, the
    /// decomposition of the pinned `DT4` at `t = 1/4`, because `M_v(DT3) = DT4`
    /// (see `tests/line_monodromy.rs` for the corpus-wide version).
    #[test]
    fn a_reciprocal_vector_shift_of_the_parameter_is_the_monodromy_image() {
        let contexts = subgroups();
        let cases: [(usize, &str); 40] = [
            (11328, "DT3"),
            (11329, "DT3"),
            (11330, "DT3"),
            (11331, "DT3"),
            (11332, "DT3"),
            (11333, "DT3"),
            (14430, "DT1"),
            (14432, "DT1"),
            (14504, "DT1"),
            (14506, "DT1"),
            (14723, "DT1"),
            (14726, "DT1"),
            (11328, "DT1"),
            (11329, "DT1"),
            (11330, "DT1"),
            (11331, "DT1"),
            (11332, "DT1"),
            (11333, "DT1"),
            (14430, "DT2"),
            (14432, "DT2"),
            (14504, "DT2"),
            (14506, "DT2"),
            (14723, "DT2"),
            (14726, "DT2"),
            (11328, "DT4"),
            (11329, "DT4"),
            (11330, "DT4"),
            (11331, "DT4"),
            (11332, "DT4"),
            (11333, "DT4"),
            (14430, "DT4"),
            (14432, "DT4"),
            (14504, "DT4"),
            (14506, "DT4"),
            (14723, "DT4"),
            (14726, "DT4"),
            (14430, "DT3"),
            (14432, "DT3"),
            (14723, "DT4"),
            (14726, "DT3"),
        ];
        let mut checked = 0usize;
        for (ordinal, label) in cases {
            let Some(subgroup) = contexts.get(&ordinal) else {
                continue;
            };
            let embedding =
                SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
            let table = line_table(subgroup.parent_sg, label);
            let direction = line_direction(table).expect("the frozen direction parses");
            let map = line_monodromy::monodromy(subgroup.parent_sg, &direction)
                .expect("frozen line direction is reciprocal");
            for (steps, numerator) in [(1usize, 5i128), (2, 9), (3, 13)] {
                let orbit = map
                    .orbit(label, steps)
                    .unwrap_or_else(|| panic!("ordinal {ordinal} {label}: orbit {steps}"));
                let image = orbit[steps];
                let image_table = line_table(subgroup.parent_sg, image);
                let reference = subduce_line_at_parameter(
                    subgroup,
                    &embedding,
                    image_table,
                    official(),
                )
                .unwrap_or_else(|error| {
                    panic!("ordinal {ordinal} {label} image {image} t=1/4: {error}")
                });
                let shifted = subduce_line_at_parameter(
                    subgroup,
                    &embedding,
                    table,
                    Rat::new(numerator, 4).unwrap(),
                )
                .unwrap_or_else(|error| panic!("ordinal {ordinal} {label} t={numerator}/4: {error}"));
                // The premise of the contract: the two parameters differ by an
                // integer multiple of the frozen direction, i.e. by a parent
                // reciprocal lattice vector.
                let delta = shifted
                    .wave_vector()
                    .checked_sub(reference.wave_vector())
                    .expect("the wave vectors are both exact");
                let factor = Rat::new(numerator, 4)
                    .unwrap()
                    .checked_sub(official())
                    .expect("the parameter step");
                let expected = Vec3R::new([
                    factor.checked_mul(direction.get(0)).unwrap(),
                    factor.checked_mul(direction.get(1)).unwrap(),
                    factor.checked_mul(direction.get(2)).unwrap(),
                ]);
                assert_eq!(
                    delta, expected,
                    "ordinal {ordinal} {label} t={numerator}/4 must be {factor} direction steps \
                     away from the official parameter"
                );
                assert_eq!(
                    target_rows(&shifted),
                    target_rows(&reference),
                    "ordinal {ordinal} {label} t={numerator}/4 must transport to {image}"
                );
                assert_eq!(
                    shifted
                        .blocks()
                        .iter()
                        .map(|block| (block.star_size(), block.arm_count(), block.block_dimension()))
                        .collect::<Vec<_>>(),
                    reference
                        .blocks()
                        .iter()
                        .map(|block| (block.star_size(), block.arm_count(), block.block_dimension()))
                        .collect::<Vec<_>>(),
                    "ordinal {ordinal} {label} t={numerator}/4 block geometry"
                );
            }
            checked += 1;
        }
        assert_eq!(checked, 40, "every gauge-regression context must be present");
    }

    /// A constructed star reports the little dimension of its own
    /// representation.  The two-dimensional family used to be reported as one
    /// (the method was hard-coded), which no production path read but which
    /// would mislead any caller that asked.
    #[test]
    fn a_two_dimensional_constructed_star_reports_its_dimension() {
        let cell = Lattice::new(exact_primitive_basis(43).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        let q = Vec3R::new([rat(0, 1), rat(1, 1), rat(1, 2)]);
        let components = constructed_child_components_at(43, &q, &reciprocal).unwrap();
        assert_eq!(components.len(), 1, "one two-dimensional irrep");
        let component = &components[0];
        assert_eq!(component.dimension, 2);
        let ComponentCharacters::Constructed(rep) = &component.characters else {
            panic!("a constructed component carries a constructed source");
        };
        let star = ConstructedStar::new(43, component.base_k, rep.clone(), component.dimension)
            .expect("child #43 constructed star");
        assert_eq!(star.dimension(), 2);
        assert!(star.arm_count() > 0);
    }

    /// The fail-closed boundary itself, end to end in the block builder: a
    /// folded child star whose little co-group is outside the constructed
    /// families must report `MissingChildStarData` -- not a partial block, and
    /// not a silent zero.
    ///
    /// Closed coverage means no pinned probe reaches this point any more, so the
    /// star is built by hand at the point `constructed_targets_available`
    /// rejects and fed to the same `build_block` stage the entry point uses.
    /// This is the Err path `select_representative` owns; the audits that used
    /// to carry it were replaced by positive tests once the batches closed the
    /// gaps, which left the boundary untested until this witness.
    #[test]
    fn an_out_of_scope_co_group_still_reports_missing_child_star_data() {
        let child = 221;
        let cell = Lattice::new(exact_primitive_basis(child).expect("child #221 basis"))
            .expect("child lattice");
        let reciprocal = cell.reciprocal().expect("child reciprocal");
        // Child #221 at (0, 0, 1/2): a sixteen-element little co-group, beyond
        // the one-dimensional batch (`MAX_ORDER = 6`) and outside both gated
        // projective families.
        let q = Vec3R::new([rat(0, 1), rat(0, 1), rat(1, 2)]);
        assert!(
            stored_child_components_at(child, &q, &reciprocal)
                .expect("stored lookup")
                .is_empty(),
            "the witness must not be answerable from pinned rows"
        );
        assert!(
            !constructed_targets_available(child, &q, &reciprocal).expect("catalogue lookup"),
            "the witness must stay outside the constructed families"
        );
        let star = FoldedStar::from_parts(vec![FoldedPoint::from_parts(q, vec![0])], 1, 1, 1);
        let built = embedding(225, "X1+", "P3");
        assert_eq!(built.subgroup_sg(), child);
        let parent_star = ScalarStar::new(probe(225, "X1+")).expect("parent star");
        let error = build_block(
            &built,
            &ArmCharacterSource::Cir(&parent_star),
            &star,
            &cell,
            &reciprocal,
        )
        .expect_err("an out-of-scope co-group must fail closed");
        assert!(
            matches!(error, FullStarError::MissingChildStarData { sg, .. } if sg == child),
            "{error:?}"
        );
    }

    /// A constructed target is induced over the **child's own star**, not
    /// answered by its little-group character.
    ///
    /// Child #2 (P-1) at a general `q = (1/4, 1/3, 0)` has a two-arm star
    /// `{q, -q}` and a trivial little co-group, so the induced character is
    /// hand-computable: the identity is the star size, a lattice translation is
    /// `2 cos(2 pi q.L)`, and the inversion moves both arms and therefore
    /// contributes zero.
    #[test]
    fn a_constructed_star_induces_over_the_child_star() {
        let q = Vec3R::new([rat(1, 4), rat(1, 3), rat(0, 1)]);
        let star = ConstructedStar::new(2, q, ConstructedLittleRep::BlochPhase { q }, 1)
            .expect("child #2 constructed star");
        assert_eq!(star.arm_count(), 2, "P-1 sends q to -q, two arms");
        assert_eq!(star.dimension(), 1);
        let translation = |values: [i32; 3]| {
            ExactSeitz::new(
                IDENTITY_ROTATION,
                Vec3R::new(values.map(|value| rat(i128::from(value), 1))),
            )
        };
        let identity = star.character(&translation([0, 0, 0])).expect("identity");
        assert!((identity - Complex64::new(2.0, 0.0)).norm() < SUBDUCTION_TOLERANCE);
        // exp(2 pi i (1/4)) + exp(-2 pi i (1/4)) = 0 exactly.
        let quarter = star.character(&translation([1, 0, 0])).expect("phase");
        assert!(quarter.norm() < SUBDUCTION_TOLERANCE, "{quarter}");
        // exp(2 pi i (1/3)) + exp(-2 pi i (1/3)) = -1 exactly.
        let third = star.character(&translation([0, 1, 0])).expect("phase");
        assert!((third - Complex64::new(-1.0, 0.0)).norm() < SUBDUCTION_TOLERANCE);
        // The inversion exchanges the two arms, so neither is fixed.
        let inversion = star
            .character(&ExactSeitz::new(
                [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
                Vec3R::new([rat(0, 1); 3]),
            ))
            .expect("moved arms contribute zero");
        assert!(inversion.norm() < SUBDUCTION_TOLERANCE, "{inversion}");
    }

    /// The constructed source is only a fallback: it exists for child #1, it
    /// carries no borrowed row, label or CIR number, and points that differ by a
    /// child reciprocal lattice vector are one identity, not two.
    #[test]
    fn constructed_targets_have_their_own_identity_and_no_borrowed_labels() {
        let reciprocal = Lattice::new(exact_primitive_basis(1).expect("child #1 basis"))
            .expect("lattice")
            .reciprocal()
            .expect("reciprocal");
        let q = Vec3R::new([rat(-1, 4), rat(-1, 4), rat(-1, 1)]);
        let components = constructed_child_components_at(1, &q, &reciprocal).expect("components");
        assert_eq!(components.len(), 1);
        let component = &components[0];
        assert!(component.record.is_none());
        assert!(component.label.is_none());
        assert!(component.bc.is_none());
        assert!(component.irnumber.is_none());
        assert_eq!(component.dimension, 1);
        assert!(matches!(
            component.component,
            SubductionComponent::Constructed { .. }
        ));
        // No pinned row is reachable there, so the stored source stays empty.
        assert!(
            stored_child_components_at(1, &q, &reciprocal)
                .expect("stored lookup")
                .is_empty()
        );
        // q and q + G are the same representation: one canonical identity.
        let shifted = q
            .checked_add(&Vec3R::new([rat(1, 1), rat(0, 1), rat(0, 1)]))
            .expect("shift");
        let other = constructed_child_components_at(1, &shifted, &reciprocal).expect("components");
        assert_eq!(other.len(), 1);
        assert_eq!(other[0].component, component.component);
        let operation = ExactSeitz::new(
            IDENTITY_ROTATION,
            Vec3R::new([rat(1, 1), rat(2, 1), rat(-3, 1)]),
        );
        let (ComponentCharacters::Constructed(first), ComponentCharacters::Constructed(second)) =
            (&component.characters, &other[0].characters)
        else {
            panic!("constructed components carry constructed characters");
        };
        let left = first.character(&operation).expect("phase");
        let right = second.character(&operation).expect("phase");
        assert!((left - right).norm() < SUBDUCTION_TOLERANCE);
        // The construction criterion is the little co-group, not the child
        // number.  Child #45 (P-4m2) at a general point has a trivial one, so
        // the same exact Bloch phase answers it.
        let reciprocal_45 = Lattice::new(exact_primitive_basis(45).expect("child #45 basis"))
            .expect("lattice")
            .reciprocal()
            .expect("reciprocal");
        let general_45 = constructed_child_components_at(45, &q, &reciprocal_45)
            .expect("general point construction");
        assert_eq!(general_45.len(), 1);
        assert_eq!(general_45[0].dimension, 1);
        assert!(matches!(
            general_45[0].component,
            SubductionComponent::Constructed { .. }
        ));
        // A non-trivial little co-group is constructed too when all of its
        // projective irreps are one-dimensional: child #45's Gamma point has a
        // four-element abelian co-group with a coboundary cocycle, so the
        // catalogue is complete there.
        let gamma = constructed_child_components_at(45, &Vec3R::new([rat(0, 1); 3]), &reciprocal_45)
            .expect("gamma lookup");
        assert_eq!(gamma.len(), 4, "four one-dimensional characters");
        assert!(gamma.iter().all(|component| component.dimension == 1
            && matches!(component.component, SubductionComponent::Constructed { .. })));
        // Beyond the batch's order cap nothing is constructed: a cubic child's
        // Gamma point is fixed by all 48 rotations.
        let reciprocal_221 = Lattice::new(exact_primitive_basis(221).expect("child #221 basis"))
            .expect("lattice")
            .reciprocal()
            .expect("reciprocal");
        assert!(
            constructed_child_components_at(221, &Vec3R::new([rat(0, 1); 3]), &reciprocal_221)
                .expect("gamma lookup")
                .is_empty(),
            "a co-group beyond the one-dimensional batch must not be constructed"
        );
    }

    /// Ordinal 1045 (SG 45 `S1S2` / C1, probe `W1W1` -> child #1) used to stop
    /// at `MissingChildStarData`; the constructed source now answers every folded
    /// star, and the entry point's own dimension and reconstruction checks pass.
    #[test]
    fn a_child_p1_gap_is_answered_by_constructed_targets() {
        let subgroup = subgroup_of(45, "S1S2", "C1");
        let built = embedding(45, "S1S2", "C1");
        assert_eq!(built.subgroup_sg(), 1);
        let probe = probe(45, "W1W1");
        let result = subduce_full_star_with_embedding(&subgroup, &built, probe)
            .expect("child #1 folded stars are constructed");
        assert!(!result.blocks().is_empty());
        let mut constructed = 0usize;
        for block in result.blocks() {
            for target in block.targets() {
                assert!(
                    target.ml.is_none() && target.row_ml.is_none() && target.irnumber.is_none(),
                    "a constructed target must not borrow a stored identity"
                );
                assert!(matches!(
                    target.component,
                    SubductionComponent::Constructed { .. }
                ));
                assert_eq!(target.dimension, 1);
                constructed += 1;
            }
        }
        assert!(constructed > 0);
        // The entry point already compared the reconstruction with the parent
        // character; the dimension conservation is asserted again here.
        assert_eq!(result.covered_dimension(), result.parent_dimension());

        // A constructed target's identity is its point **reduced modulo the
        // child reciprocal lattice**, not the raw folded representative: the
        // lookup only accepts the canonical identity, so an equivalent
        // coordinate can never answer a false zero.
        let reciprocal = Lattice::new(exact_primitive_basis(1).expect("child #1 basis"))
            .expect("lattice")
            .reciprocal()
            .expect("reciprocal");
        let mut pinned = Vec::new();
        for block in result.blocks() {
            for target in block.targets() {
                let SubductionComponent::Constructed { q, index } = target.component else {
                    continue;
                };
                let identity = Vec3R::new(q);
                let reduced = reciprocal
                    .reduce(block.q())
                    .expect("reduce the folded point")
                    .representative;
                assert_eq!(
                    reduced, identity,
                    "the identity must be the canonical point of the block"
                );
                assert_eq!(
                    block.constructed_multiplicity(target.component),
                    target.multiplicity
                );
                // Any other identity answers zero, never a stale hit.
                assert_eq!(
                    block.constructed_multiplicity(SubductionComponent::Constructed {
                        q: [rat(1, 2), rat(1, 2), rat(1, 2)],
                        index,
                    }),
                    0
                );
                pinned.push(([q[0], q[1], q[2]], target.multiplicity));
            }
        }
        // Two folded stars, each one constructed target of multiplicity 2; the
        // raw representatives reduce to these canonical identities.
        assert_eq!(
            pinned,
            [
                ([rat(3, 4), rat(3, 4), rat(0, 1)], 2),
                ([rat(1, 4), rat(1, 4), rat(0, 1)], 2)
            ],
            "ordinal 1045 constructed targets"
        );
    }
}
