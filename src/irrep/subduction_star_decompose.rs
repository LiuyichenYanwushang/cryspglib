//! Staged full-star decomposition of an ordinary scalar parent irrep (task 8b
//! of `docs/full-irrep-subduction-plan.md`).
//!
//! [`subduce_full_star_with_embedding`] takes the full parent star built by
//! [`OrdinaryStar`], folds its arms into the embedding with
//! [`OrdinaryStar::folded_stars`], and decomposes every folded child star:
//!
//! * the child table stores one representative **arm** per irrep, so the
//!   representative `q` of a child star is the folded point whose class matches
//!   a stored child `k` exactly modulo the **child** reciprocal lattice
//!   (centring extinctions included), searched over the whole folded orbit —
//!   never `points()[0]` by assumption and never a nearest match;
//! * on `H_q`, the little group of `q`, the q-block character is the induced
//!   one restricted to the parent arms that fold onto `q`:
//!   `chi_q(h) = sum_i chi_seed(g_i^-1 h g_i)`, with every Seitz translation
//!   kept for the row lookup. Fixing `q` can still permute parent arms that
//!   fold onto it (e.g. SG 139 `X1+` restricted to #126); moved arms contribute
//!   zero. The full-star trace is never divided by an arm count;
//! * the q-block is decomposed against every complex child row stored at `q` by
//!   the shared `solve_character_block` solver, so the Gram
//!   identity, the non-negative integral multiplicities, the little-dimension
//!   sum and the per-operation reconstruction are the same checks the single
//!   arm entry point uses;
//! * the whole subduced full-star character on **all** subgroup coset
//!   representatives is then rebuilt from the reported child irreps induced
//!   over their own **child** stars and compared with the parent full-star
//!   character.  That reconstruction never calls the q-block builder again, so
//!   a wrong frame transport between parent arms and child stars cannot cancel
//!   out.
//!
//! Boundaries: ordinary scalar **parents** only; child ordinary rows and
//! `DistinctComponentSum` constituents are supported, while
//! `ConjugateRealification` child rows are rejected until the k/-k source-star
//! grouping is implemented.  Magnetic groups and spinor/compound parents stay
//! out of this stage.  Missing child `k`/irrep data is a typed error: a partial
//! set of blocks is never returned.

use num_complex::Complex64;

use crate::irrep::isotropy::IsotropySubgroup;
use crate::irrep::query;
use crate::irrep::types::{
    CharacterRow, CompoundCharacterSemantics, CompoundSelectedArmCharacter, IrrepRecord,
    IrrepSourceIdentity,
};
use crate::mathfunc::Mat3I;

use super::super::{
    ExactSeitz, Lattice, Rat, SUBDUCTION_TOLERANCE, SubductionComponent, SubductionError,
    SubductionTarget, SubgroupEmbedding, Vec3R, exact_primitive_basis, inline_k_vector,
    reduce_operations, shift_operations, solve_character_block, strict_sg_hall_ops,
    validate_subduction_context,
};
use super::{FoldedStar, OrdinaryStar, StarArm, StarError, collect_arms, induced_character};

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
    /// A needed child row is a `ConjugateRealification`; its k/-k source-star
    /// grouping is a follow-up stage.
    #[error(
        "child irrep {ml} of space group {sg} is a ConjugateRealification row; k/-k source \
         stars are not handled by this stage"
    )]
    ConjugateRealificationTarget {
        /// Child space group number.
        sg: u8,
        /// Child row label.
        ml: &'static str,
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
        /// `OrdinaryStar::dimension()`.
        expected: u32,
        /// Sum of the reported block dimensions.
        found: u32,
    },
    /// No unique child record backed a reported ordinary target.
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
        /// `OrdinaryStar::character` of the same representative.
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
    pub ml: &'static str,
    /// Bradley-Cracknell label, for display only.
    pub bc: &'static str,
    /// Physical child row this term came from.
    pub row_ml: &'static str,
    /// Which complex constituent of that row this term is.
    pub component: SubductionComponent,
    /// Complex dimension of the child little-group irrep.
    pub dimension: u8,
    /// Multiplicity in the child little-group representation at `q`.
    pub multiplicity: u32,
    /// Frozen CIR source number of this term (`record.source_identity()` for an
    /// ordinary row, the constituent's own `irnumber` for a compound row).
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

    /// The stored child `k` that `q` matched modulo the child reciprocal
    /// lattice (the child table's representative arm).
    pub const fn stored_k(&self) -> &Vec3R {
        &self.stored_k
    }

    /// Number of folded points in this child star.
    pub const fn star_size(&self) -> usize {
        self.star_size
    }

    /// Number of parent arms in this child star.
    pub const fn arm_count(&self) -> usize {
        self.arm_count
    }

    /// Parent full-star subspace carried by this block:
    /// `arm_count × selected-arm dimension`.
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

/// The complete decomposition of one ordinary scalar parent full star on one
/// embedded subgroup.
#[derive(Debug, Clone)]
pub struct FullStarSubduction {
    parent_sg: u8,
    parent_ml: &'static str,
    parent_bc: &'static str,
    parent_irnumber: u32,
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

    /// Frozen CIR source number of the parent probe.
    pub const fn parent_irnumber(&self) -> u32 {
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

/// A child full-star character rebuilt from one **constituent** row.
///
/// The child's own `OrdinaryStar` handles ordinary rows; CIR constituents of a
/// `DistinctComponentSum` row are not records of the child table, so their star
/// is assembled from the same child Hall transporters and evaluated with the
/// same induced-character formula.
#[derive(Debug, Clone)]
struct ConstituentStar {
    sg: u8,
    ml: &'static str,
    seed_k: Vec3R,
    row: CharacterRow,
    arms: Vec<StarArm>,
    lattice: Lattice,
    reciprocal: Lattice,
    operations: Vec<ExactSeitz>,
}

impl ConstituentStar {
    fn new(
        sg: u8,
        seed_k: &Vec3R,
        row: &CharacterRow,
        ml: &'static str,
    ) -> Result<Self, FullStarError> {
        let lattice = Lattice::new(exact_primitive_basis(sg)?)?;
        let reciprocal = lattice.reciprocal()?;
        let hall = strict_sg_hall_ops(sg)?;
        let operations = reduce_operations(&hall.operations, &lattice)?;
        let arms = collect_arms(sg, seed_k, &reciprocal, &hall.operations, false)?;
        Ok(Self {
            sg,
            ml,
            seed_k: *seed_k,
            row: row.clone(),
            arms,
            lattice,
            reciprocal,
            operations,
        })
    }

    fn character(&self, operation: &ExactSeitz) -> Result<Complex64, FullStarError> {
        let reduced = operation.reduce(&self.lattice)?;
        if !self.operations.contains(&reduced) {
            return Err(StarError::OperationNotInParentGroup { sg: self.sg }.into());
        }
        Ok(induced_character(
            &self.arms,
            &self.row,
            self.ml,
            &self.lattice,
            &self.reciprocal,
            &self.seed_k,
            operation,
        )?)
    }
}

/// The child full-star character of one reported target.
#[derive(Debug, Clone)]
enum ChildStarEvaluator {
    Ordinary(OrdinaryStar),
    Constituent(ConstituentStar),
}

impl ChildStarEvaluator {
    fn character(&self, operation: &ExactSeitz) -> Result<Complex64, FullStarError> {
        match self {
            Self::Ordinary(star) => Ok(star.character(operation)?),
            Self::Constituent(star) => star.character(operation),
        }
    }
}

// ── Entry point ──────────────────────────────────────────────────────────────

/// Decompose the full star of an ordinary scalar parent irrep on the subgroup
/// selected by an isotropy record and a prebuilt embedding.
///
/// `subgroup`, `embedding` and `probe` are revalidated as one context before any
/// block is built (see `validate_subduction_context`).  Compound and spinor
/// probes are rejected with [`SubductionError::UnsupportedCharacterSpace`], and
/// a folded child star without stored child data is a
/// [`FullStarError::MissingChildStarData`] error rather than a partial result.
pub fn subduce_full_star_with_embedding(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<FullStarSubduction, FullStarError> {
    validate_subduction_context(subgroup, embedding, probe)?;
    if probe.spinor || probe.compound_metadata().is_some() {
        return Err(SubductionError::UnsupportedCharacterSpace {
            sg: probe.sg,
            ml: probe.ml.to_string(),
        }
        .into());
    }
    let star = OrdinaryStar::new(probe)?;
    decompose_ordinary_star(embedding, &star)
}

/// The reusable core: decompose an already validated full parent star.
///
/// Split out so the tests can drive a differently ordered but equivalent
/// transversal through [`OrdinaryStar::from_transporters`] and a mutated folded
/// geometry through [`decompose_folded_stars`]; both are private, so no caller
/// can reach them with an unvalidated context.
fn decompose_ordinary_star(
    embedding: &SubgroupEmbedding,
    star: &OrdinaryStar,
) -> Result<FullStarSubduction, FullStarError> {
    let folded = star.folded_stars(embedding)?;
    decompose_folded_stars(embedding, star, &folded)
}

/// Decompose one folded geometry of an already validated star.
fn decompose_folded_stars(
    embedding: &SubgroupEmbedding,
    star: &OrdinaryStar,
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
        parent_irnumber: ordinary_irnumber(star.probe())?,
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

// ── One folded child star ────────────────────────────────────────────────────

/// Decompose one folded child star: pick its stored representative arm, build
/// the `H_q` q-block character and solve it against the stored child rows.
fn build_block(
    embedding: &SubgroupEmbedding,
    star: &OrdinaryStar,
    folded_star: &FoldedStar,
    child_cell: &Lattice,
    child_reciprocal: &Lattice,
) -> Result<FullStarBlock, FullStarError> {
    let child_sg = embedding.subgroup_sg();
    let (point_index, records) = select_representative(child_sg, folded_star, child_reciprocal)?;
    let point = &folded_star.points()[point_index];
    let q = *point.q();
    let q_key = [q.get(0), q.get(1), q.get(2)];

    // A `ConjugateRealification` row is a sum over a k/-k pair; expanding it as
    // seed + conjugate here would silently answer a question this stage does
    // not implement.
    ensure_no_realification(child_sg, &records)?;

    // `H_q`: the subgroup coset representatives whose child-frame rotation
    // fixes `q` modulo the child reciprocal lattice.  Lattice translations stay
    // in the operations; only the `child_shift` of the embedding is undone, to
    // land
    // in the child Hall frame the shipped rows live in.
    let mut parent_operations: Vec<ExactSeitz> = Vec::new();
    let mut child_operations: Vec<ExactSeitz> = Vec::new();
    for operation in embedding.representatives() {
        let child = embedding.transform().unmap_operation(operation)?;
        if !child_reciprocal.preserves(child.rotation(), &q)? {
            continue;
        }
        parent_operations.push(*operation);
        child_operations.push(child);
    }
    if parent_operations.is_empty() {
        return Err(FullStarError::MissingIdentityOperation);
    }
    let pulled_back = shift_operations(&child_operations, &embedding.child_shift().checked_neg()?)?;
    let identity = identity_position(&pulled_back, child_cell)?;

    // The q-block character: only the parent arms folding onto this q, each
    // fixed or moved modulo the parent reciprocal lattice by the same
    // `induced_character` formula the full star uses.
    let arms_at_q: Vec<StarArm> = point
        .arm_indices()
        .iter()
        .map(|index| star.arms()[*index])
        .collect();
    if arms_at_q.is_empty() {
        return Err(FullStarError::EmptyQBlock { q: q_key });
    }
    let seed_dimension = star.selected_row().dimension();
    let q_block_dimension = dimension_product(point.arm_count(), seed_dimension)?;
    let mut parent_characters = Vec::with_capacity(parent_operations.len());
    for operation in &parent_operations {
        parent_characters.push(induced_character(
            &arms_at_q,
            star.selected_row(),
            star.probe_ml(),
            &star.parent_lattice,
            &star.parent_reciprocal,
            star.seed_k(),
            operation,
        )?);
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

    let solved = solve_character_block(
        child_sg,
        &records,
        &parent_characters,
        &pulled_back,
        child_cell,
        &q,
    )?;
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
        let irnumber = target_irnumber(child_sg, &records, target)?;
        let evaluator = build_evaluator(child_sg, &records, target)?;
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
            irnumber,
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
    let stored_k = inline_k_vector(records[0])?;
    Ok(FullStarBlock {
        q,
        stored_k,
        star_size,
        arm_count: folded_star.arm_count(),
        block_dimension: folded_star.block_dimension(),
        little_dimension,
        targets,
        evaluators,
    })
}

/// The folded point of a child star whose class matches stored child data.
///
/// The child table stores a representative arm, so the whole orbit is searched
/// in canonical order; the first point with stored data fixes the frame the
/// block is decomposed in. Both `q` and the stored `k` use child coordinates;
/// they remain unreduced and are compared modulo the child reciprocal lattice.
fn select_representative(
    child_sg: u8,
    folded_star: &FoldedStar,
    child_reciprocal: &Lattice,
) -> Result<(usize, Vec<&'static IrrepRecord>), FullStarError> {
    for (index, point) in folded_star.points().iter().enumerate() {
        let records = child_records_at(child_sg, point.q(), child_reciprocal)?;
        if !records.is_empty() {
            return Ok((index, records));
        }
    }
    let q = folded_star
        .points()
        .first()
        .map_or([Rat::ZERO; 3], |point| {
            [point.q().get(0), point.q().get(1), point.q().get(2)]
        });
    Err(FullStarError::MissingChildStarData {
        sg: child_sg,
        q,
        points: folded_star.points().len(),
    })
}

/// Every non-spinor child record stored at `q` modulo the child reciprocal
/// lattice.
///
/// Spinor rows are a different representation space; the scalar parent's
/// restriction is single-valued, so they are not targets here.  No other record
/// is skipped.
fn child_records_at(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<Vec<&'static IrrepRecord>, FullStarError> {
    let mut records = Vec::new();
    for record in query::irreps_of(child_sg) {
        if record.spinor {
            continue;
        }
        let stored = inline_k_vector(record)?;
        if child_reciprocal.contains(&q.checked_sub(&stored)?)? {
            records.push(record);
        }
    }
    Ok(records)
}

/// Reject a `ConjugateRealification` row before any character is evaluated.
fn ensure_no_realification(
    child_sg: u8,
    records: &[&'static IrrepRecord],
) -> Result<(), FullStarError> {
    for record in records {
        let is_realification = record.compound_metadata().is_some_and(|metadata| {
            metadata.semantics == CompoundCharacterSemantics::ConjugateRealification
        });
        if is_realification {
            return Err(FullStarError::ConjugateRealificationTarget {
                sg: child_sg,
                ml: record.ml,
            });
        }
    }
    Ok(())
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

/// `arms × dimension` as a `u32`, with checked arithmetic.
fn dimension_product(arms: usize, dimension: usize) -> Result<u32, FullStarError> {
    u32::try_from(arms)
        .ok()
        .and_then(|arms| {
            u32::try_from(dimension)
                .ok()
                .and_then(|dimension| arms.checked_mul(dimension))
        })
        .ok_or_else(|| {
            SubductionError::RationalOverflow {
                operation: "q-block dimension",
            }
            .into()
        })
}

/// The unique ordinary child record behind a reported ordinary target.
fn ordinary_record(
    child_sg: u8,
    records: &[&'static IrrepRecord],
    ml: &'static str,
) -> Result<&'static IrrepRecord, FullStarError> {
    let mut found: Option<&'static IrrepRecord> = None;
    for record in records {
        if record.ml == ml && record.compound_metadata().is_none() {
            if found.is_some() {
                return Err(FullStarError::AmbiguousTargetSource { sg: child_sg, ml });
            }
            found = Some(record);
        }
    }
    found.ok_or(FullStarError::AmbiguousTargetSource { sg: child_sg, ml })
}

/// The frozen CIR source number of a reported target.
fn target_irnumber(
    child_sg: u8,
    records: &[&'static IrrepRecord],
    target: &SubductionTarget,
) -> Result<u32, FullStarError> {
    match target.component {
        SubductionComponent::Ordinary => {
            ordinary_irnumber(ordinary_record(child_sg, records, target.ml)?)
        }
        SubductionComponent::Constituent { irnumber, .. } => Ok(irnumber),
        SubductionComponent::RealificationSeed { .. }
        | SubductionComponent::RealificationConjugate { .. } => {
            Err(FullStarError::ConjugateRealificationTarget {
                sg: child_sg,
                ml: target.ml,
            })
        }
    }
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
fn build_evaluator(
    child_sg: u8,
    records: &[&'static IrrepRecord],
    target: &SubductionTarget,
) -> Result<ChildStarEvaluator, FullStarError> {
    match target.component {
        SubductionComponent::Ordinary => Ok(ChildStarEvaluator::Ordinary(OrdinaryStar::new(
            ordinary_record(child_sg, records, target.ml)?,
        )?)),
        SubductionComponent::Constituent { index, irnumber } => {
            let record = records
                .iter()
                .copied()
                .find(|record| record.ml == target.row_ml)
                .ok_or(FullStarError::AmbiguousTargetSource {
                    sg: child_sg,
                    ml: target.row_ml,
                })?;
            let view = record.compound_selected_arm_view().map_err(|_| {
                SubductionError::UnsupportedCharacterSpace {
                    sg: child_sg,
                    ml: record.ml.to_string(),
                }
            })?;
            let CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } = view
            else {
                return Err(FullStarError::ConjugateRealificationTarget {
                    sg: child_sg,
                    ml: record.ml,
                });
            };
            let constituent = if index == 0 { first } else { second };
            // The reported identity is taken from the stored CIR source, not
            // from the row label; a mismatch here means the block was built
            // from a different row than the one it is reported against.
            if constituent.irnumber != irnumber
                || constituent.dimension != usize::from(target.dimension)
            {
                return Err(FullStarError::TargetSourceMismatch {
                    sg: child_sg,
                    ml: constituent.label,
                });
            }
            let seed_k = inline_k_vector(record)?;
            Ok(ChildStarEvaluator::Constituent(ConstituentStar::new(
                child_sg,
                &seed_k,
                &constituent.row,
                constituent.label,
            )?))
        }
        SubductionComponent::RealificationSeed { .. }
        | SubductionComponent::RealificationConjugate { .. } => {
            Err(FullStarError::ConjugateRealificationTarget {
                sg: child_sg,
                ml: target.ml,
            })
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
    star: &OrdinaryStar,
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
            SubductionComponent::RealificationSeed { .. }
            | SubductionComponent::RealificationConjugate { .. } => {
                panic!("a realification target leaked into a full-star block")
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
                "q {q:?} and stored k {k:?} differ modulo the child reciprocal lattice",
                q = block.q(),
                k = block.stored_k()
            );
            assert!(
                query::irreps_of(child_sg).iter().any(|record| {
                    !record.spinor
                        && child_reciprocal
                            .same_mod(&inline_k_vector(record).expect("k"), block.stored_k())
                            .expect("same_mod")
                }),
                "stored k {k:?} is not a child record",
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
        assert_eq!(result.parent_irnumber(), 10719);
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
        assert_eq!(result.parent_irnumber(), 10723);
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
        assert_eq!(result.parent_irnumber(), 10723);
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
        assert_eq!(x1.parent_irnumber(), 7301);
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
        assert_eq!(n1.parent_irnumber(), 7319);
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
        assert_eq!(result.parent_irnumber(), 8288);
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
        let folded = OrdinaryStar::new(probe(167, "F1+"))
            .unwrap()
            .folded_stars(&built)
            .unwrap();
        let two_point = folded.iter().find(|star| star.star_size() == 2).unwrap();
        assert!(
            child_records_at(15, two_point.points()[0].q(), &child_reciprocal)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            select_representative(15, two_point, &child_reciprocal)
                .unwrap()
                .0,
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
        assert_eq!(result.parent_irnumber(), 379);
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
    /// set nor any block of the decomposition.
    #[test]
    fn arm_order_invariance() {
        for (parent, condensing, direction, ml) in [
            (221u8, "GM4+", "P1", "X1+"),
            (221, "GM4+", "P1", "X5+"),
            (139, "M1-", "P1", "N1+"),
            (167, "GM3+", "P1", "F1+"),
        ] {
            let built = embedding(parent, condensing, direction);
            let record = probe(parent, ml);
            let canonical = OrdinaryStar::new(record).expect("star");
            let mut transporters: Vec<ExactSeitz> = canonical
                .arms()
                .iter()
                .map(|arm| *arm.transporter())
                .collect();
            transporters.reverse();
            let permuted =
                OrdinaryStar::from_transporters(record, &transporters).expect("permuted star");
            assert_eq!(canonical.arms(), permuted.arms());
            let expected = decompose_ordinary_star(&built, &canonical).expect("canonical");
            let found = decompose_ordinary_star(&built, &permuted).expect("permuted");
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

        let compound = query::irreps_of(167)
            .iter()
            .find(|record| record.ml == "T1T2")
            .expect("compound row");
        assert!(compound.compound_metadata().is_some());
        let context_167 = subgroup_of(167, "GM3+", "P1");
        let embedding_167 = embedding(167, "GM3+", "P1");
        assert!(matches!(
            subduce_full_star_with_embedding(&context_167, &embedding_167, compound),
            Err(FullStarError::Subduction(
                SubductionError::UnsupportedCharacterSpace { sg: 167, .. }
            ))
        ));

        // A compound row of another space group is a foreign probe first: the
        // context check runs before the representation-space check.
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

    /// `ConjugateRealification` child rows are rejected with a typed error, in
    /// the block pre-scan and in the source-identity layer; the frozen contexts
    /// contain no such row, so the real SG 23 record is driven directly.
    #[test]
    fn conjugate_realification_child_rows_are_rejected() {
        let realification = query::irreps_of(23)
            .iter()
            .find(|record| record.ml == "W1W1")
            .expect("SG 23 W1W1");
        assert!(matches!(
            ensure_no_realification(23, &[realification]),
            Err(FullStarError::ConjugateRealificationTarget { sg: 23, ml: "W1W1" })
        ));
        let distinct = query::irreps_of(83)
            .iter()
            .find(|record| record.ml == "GM3+GM4+")
            .expect("compound row");
        assert!(ensure_no_realification(83, &[distinct]).is_ok());

        let target = SubductionTarget {
            sg: 23,
            ml: "W1W1",
            bc: realification.bc,
            row_ml: "W1W1",
            dimension: 1,
            component: SubductionComponent::RealificationSeed { irnumber: 716 },
            multiplicity: 1,
        };
        assert!(matches!(
            build_evaluator(23, &[realification], &target),
            Err(FullStarError::ConjugateRealificationTarget { sg: 23, ml: "W1W1" })
        ));
        assert!(matches!(
            target_irnumber(23, &[realification], &target),
            Err(FullStarError::ConjugateRealificationTarget { sg: 23, ml: "W1W1" })
        ));
    }

    /// A mutated folded geometry fails on the dimension checks or on the
    /// missing-data error; it never returns a smaller partial decomposition.
    #[test]
    fn lost_blocks_and_missing_data_fail_instead_of_partial_results() {
        let built = embedding(221, "GM4+", "P1");
        let star = OrdinaryStar::new(probe(221, "X1+")).expect("star");
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

    /// The five frozen contexts scanned probe by probe.  Every ordinary scalar
    /// probe either returns a complete decomposition that passes every
    /// invariant, or is one of the explicitly counted missing/unsupported
    /// cases; skipped spinor and compound parent rows are never counted as
    /// successes.  The counts are pinned, so a newly unsupported probe fails
    /// this gate instead of silently shrinking the tested coverage.
    #[test]
    fn frozen_contexts_report_success_missing_and_unsupported_separately() {
        let contexts = [
            (221u8, "GM4+", "P1", 83u8),
            (221, "GM4+", "P2", 12),
            (16, "R1", "P1", 22),
            (167, "GM3+", "P1", 15),
            (139, "M1-", "P1", 126),
        ];
        let mut total_ok = 0usize;
        let mut total_skipped = 0usize;
        let mut total_multi_arm = 0usize;
        let mut total_single_arm = 0usize;
        for (parent, condensing, direction, child) in contexts {
            let subgroup = subgroup_of(parent, condensing, direction);
            let built = embedding(parent, condensing, direction);
            assert_eq!(built.subgroup_sg(), child);
            let mut ok = 0usize;
            let mut skipped = 0usize;
            let mut missing = 0usize;
            let mut unsupported = 0usize;
            let mut other = Vec::new();
            for record in query::irreps_of(parent) {
                if record.spinor || record.compound_metadata().is_some() {
                    skipped += 1;
                    continue;
                }
                match subduce_full_star_with_embedding(&subgroup, &built, record) {
                    Ok(result) => {
                        assert_invariants(&result, parent, child);
                        if OrdinaryStar::new(record).expect("star").arm_count() > 1 {
                            total_multi_arm += 1;
                        } else {
                            total_single_arm += 1;
                        }
                        ok += 1;
                    }
                    Err(FullStarError::MissingChildStarData { .. }) => missing += 1,
                    Err(FullStarError::ConjugateRealificationTarget { .. })
                    | Err(FullStarError::Subduction(
                        SubductionError::UnsupportedCharacterSpace { .. },
                    )) => unsupported += 1,
                    Err(error) => other.push((record.ml, error)),
                }
            }
            assert!(
                other.is_empty(),
                "SG {parent} {condensing} {direction}: unexpected errors {other:?}"
            );
            assert_eq!(
                (ok, missing, unsupported, skipped),
                match (parent, direction) {
                    (221, "P1") => (40, 0, 0, 20),
                    (221, "P2") => (40, 0, 0, 20),
                    (16, "P1") => (32, 0, 0, 8),
                    (167, "P1") => (12, 0, 0, 15),
                    (139, "P1") => (37, 0, 0, 16),
                    _ => unreachable!(),
                },
                "context SG {parent} {condensing} {direction}"
            );
            total_ok += ok;
            total_skipped += skipped;
        }
        assert_eq!(total_ok, 161);
        assert_eq!(total_skipped, 79);
        // Pin actual single-arm and multi-arm coverage separately.
        assert_eq!(total_multi_arm, 62);
        assert_eq!(total_single_arm, 99);
    }
}
