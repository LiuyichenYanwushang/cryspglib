//! Full-star induction and folded child-star geometry for ordinary scalar
//! probes (task 8a of `docs/full-irrep-subduction-plan.md`).
//!
//! This stage is deliberately narrow:
//!
//! * only **ordinary scalar** parent probes are accepted.  Compound rows and
//!   spinors fail with [`SubductionError::UnsupportedCharacterSpace`] before any
//!   typed row is built, because their full-star trace is a different object.
//! * the parent operations come from [`strict_sg_hall_ops`] — the recorded data
//!   Hall setting — never from `SymmetryOps::from_sg`, and the parent
//!   translation lattice is the exact primitive basis of the same frame.
//!   Every arm is the contragredient image `R^-T k` of the stored selected `k`
//!   by a full Seitz transporter of that authoritative group, one transporter
//!   per equivalence class modulo the parent reciprocal lattice.
//! * the full-star character is the induced one: an arm that the operation
//!   moves to another arm contributes **zero**, a fixed arm contributes
//!   `chi_seed(g_i^-1 h g_i)`.  The conjugated Seitz operation is evaluated
//!   with the **seed** wave vector through [`super::character_of`], so the
//!   Bloch phase of every stored representative survives.
//!
//! [`decompose`] restricts the ordinary full star to an embedded subgroup and
//! returns child irreps and multiplicities. [`folded_stars`] provides the
//! underlying geometry: it groups parent arms into child `k`-stars and checks
//! their dimensions. Neither path divides a full-star trace by a star size to
//! obtain a selected-arm character.
//!
//! [`folded_stars`]: OrdinaryStar::folded_stars

/// Exact little co-groups and their one-dimensional projective characters
/// (R4 batch 2a).  The catalogue is only used when the solver finds exactly
/// `|P_q|` characters, i.e. when every irreducible projective representation of
/// the co-group is one-dimensional; see the module documentation.
#[path = "subduction_catalogue.rs"]
pub mod catalogue;
#[path = "subduction_star_decompose.rs"]
pub mod decompose;
#[path = "subduction_scalar_star.rs"]
pub mod scalar_star;

use num_complex::Complex64;

use crate::irrep::query;
use crate::irrep::types::{CharacterRow, IrrepRecord};
use crate::mathfunc::Mat3I;

use super::{
    ExactSeitz, Lattice, Mat3R, Rat, SUBDUCTION_TOLERANCE, SubductionError, SubgroupEmbedding,
    Vec3R, bloch_phase, character_of, exact_primitive_basis, fold_wave_vector, inline_k_vector,
    reduce_operations, strict_sg_hall_ops,
};

/// The identity rotation, as stored in every Hall operation table.
const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

// ── Errors ───────────────────────────────────────────────────────────────────

/// Errors from the full-star and folded-star adapter.
///
/// Everything the exact subduction layer can report keeps its own variant
/// ([`StarError::Subduction`]); the variants below are specific to the arm
/// bookkeeping this stage adds.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum StarError {
    /// An error from the exact rational/affine layer.
    #[error(transparent)]
    Subduction(#[from] SubductionError),
    /// A frozen parametric-k source has no character for one of its own
    /// little-group rotations, so the line arm source cannot answer.
    #[allow(dead_code)]
    #[error("frozen line source {label} of space group {sg} has no character for its own rotation")]
    MissingFrozenRotation { sg: u8, label: &'static str },
    /// The transport list does not contain the identity operation itself.
    ///
    /// The list is a *transversal*, and the identity is the canonical
    /// representative of the little group coset, so a list that only carries
    /// `{I|L}` for a lattice vector `L` is rejected rather than silently
    /// re-based.
    #[error("the transport list for space group {sg} has no identity operation")]
    MissingIdentityTransporter { sg: u8 },
    /// Two transporters select the same arm class.
    #[error("transporters {first} and {second} select the same star arm")]
    DuplicateArm { first: usize, second: usize },
    /// A transporter is not an operation of the parent group modulo `L_G`.
    #[error(
        "transporter {index} is not an operation of space group {sg} modulo the parent lattice"
    )]
    TransporterNotInParentGroup { sg: u8, index: usize },
    /// An operation handed to [`OrdinaryStar::character`] is not a parent
    /// group element modulo `L_G`.
    #[error("operation is not in space group {sg} modulo the parent lattice")]
    OperationNotInParentGroup { sg: u8 },
    /// A constructed little-group representation was asked for the character of
    /// a rotation it was not built from.
    ///
    /// Only a trivial little co-group is constructed, and a contribution is only
    /// ever requested for an operation that fixes the arm it is conjugated into;
    /// a non-identity rotation therefore means the caller used the wrong little
    /// group, and that fails closed instead of answering a phase.
    #[error(
        "constructed little-group representation at q = ({}, {}, {}) cannot be evaluated on a \
         rotation it was not built from",
        q[0],
        q[1],
        q[2]
    )]
    ConstructedRotationNotCovered { q: [Rat; 3] },
    /// The child rotations fixing an exact point do not close under products,
    /// or the product of two of them leaves the child lattice.
    #[error(
        "the child operations fixing q = ({}, {}, {}) do not form a closed little co-group",
        q[0],
        q[1],
        q[2]
    )]
    LittleCoGroupNotClosed { q: [Rat; 3] },
    /// The transversals do not cover the whole parent star.
    #[error(
        "{represented} transporters of {ml} do not cover all {operations} operations of \
         space group {sg}; first uncovered image ({}, {}, {})",
        missing_k[0], missing_k[1], missing_k[2]
    )]
    IncompleteStar {
        sg: u8,
        ml: &'static str,
        represented: usize,
        operations: usize,
        missing_k: Box<[Rat; 3]>,
    },
    /// The induced dimension does not match the selected row and the record.
    #[error(
        "identity character {identity_trace} of the full star of {ml} does not match \
         {selected_dimension} x {arms} = {product} or the record dimension {probe_dim}"
    )]
    StarDimensionMismatch {
        ml: &'static str,
        probe_dim: u8,
        selected_dimension: usize,
        arms: usize,
        product: i64,
        identity_trace: Complex64,
    },
    /// The embedding belongs to another parent.
    #[error(
        "embedding parent {embedding_parent_sg} does not match the star parent {star_parent_sg}"
    )]
    EmbeddingParentMismatch {
        star_parent_sg: u8,
        embedding_parent_sg: u8,
    },
    /// A folded `q` leaves the folded point set under a subgroup rotation.
    #[error(
        "child star at q = ({}, {}, {}) is not closed under the subgroup rotations",
        q[0], q[1], q[2]
    )]
    ChildStarNotClosed { q: [Rat; 3] },
    /// One orbit carries different arm counts at different `q` points.
    #[error("child star has non-uniform arm blocks: {arms} and {other_arms} parent arms")]
    NonUniformArmBlock { arms: usize, other_arms: usize },
    /// The scalar components of one parent record do not share one selected-arm
    /// dimension, so the folded geometry's uniform block dimension is undefined.
    ///
    /// The frozen scalar table has equal constituent dimensions for every
    /// compound record (pinned by the full-table census); a record that breaks
    /// that invariant fails closed instead of being folded with a mixed block.
    #[error(
        "scalar components carry different selected dimensions {first} and {other}; the folded \
         geometry requires one dimension"
    )]
    NonUniformComponentDimension { first: usize, other: usize },
    /// The assembled scalar full-star dimension disagrees with the record.
    #[error(
        "scalar full-star dimension of {ml} is {found} (identity character {identity}) but the \
         record dimension is {record}"
    )]
    ScalarDimensionMismatch {
        ml: &'static str,
        record: u32,
        found: u32,
        identity: Complex64,
    },
    /// A component transversal handed to the scalar adapter is not a parent
    /// group element modulo the parent lattice.
    #[error(
        "component {component} transporter {index} is not an operation of space group {sg} \
         modulo the parent lattice"
    )]
    ComponentTransporterNotInParentGroup {
        sg: u8,
        component: usize,
        index: usize,
    },
    /// The scalar adapter was handed the wrong number of component transversals.
    #[error("scalar adapter needs {expected} component transversals but got {found}")]
    ComponentTransporterCount { expected: usize, found: usize },
    /// A folded arm index does not address a parent scalar arm.
    #[error("folded arm index {index} is out of range for {arms} parent scalar arms")]
    UnknownArmIndex { index: usize, arms: usize },
    /// The folded stars do not account for every parent arm.
    #[error("folded stars cover {folded} of {arms} parent arms")]
    FoldLostArms { folded: usize, arms: usize },
}

// ── Exact arm bookkeeping ────────────────────────────────────────────────────

/// One arm of a full parent star.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StarArm {
    wave_vector: Vec3R,
    transporter: ExactSeitz,
}

impl StarArm {
    /// The exact arm wave vector `R_i^-T k` (unreduced).
    pub const fn wave_vector(&self) -> &Vec3R {
        &self.wave_vector
    }

    /// The parent Hall operation that produced it.
    pub const fn transporter(&self) -> &ExactSeitz {
        &self.transporter
    }
}

/// The contragredient image `R^-T k` of one transporter.
fn arm_wave_vector(seed: &Vec3R, transporter: &ExactSeitz) -> Result<Vec3R, SubductionError> {
    Mat3R::from_ints(transporter.rotation())
        .inverse()?
        .transpose()
        .checked_mul_vector(seed)
}

/// Stable total order on the exact rational encoding (not numerical order).
/// Keep denominators: scaling each vector separately and retaining only its
/// numerators would give `(1/2,0,0)` and `(1,0,0)` the same key.
fn canonical_key(vector: &Vec3R) -> [(i128, i128); 3] {
    vector
        .as_array()
        .map(|value| (value.numerator(), value.denominator()))
}

/// One arm per equivalence class of `transporters` modulo the parent
/// reciprocal lattice, identity first and canonically ordered afterwards.
///
/// `duplicates_are_errors` separates the public transversal constructor from
/// [`OrdinaryStar::new`], which scans the complete Hall list and therefore sees
/// every element of a little-group coset.
fn collect_arms(
    parent_sg: u8,
    seed: &Vec3R,
    parent_reciprocal: &Lattice,
    transporters: &[ExactSeitz],
    duplicates_are_errors: bool,
) -> Result<Vec<StarArm>, StarError> {
    let identity_position = transporters
        .iter()
        .position(|transporter| *transporter == ExactSeitz::identity())
        .ok_or(StarError::MissingIdentityTransporter { sg: parent_sg })?;

    let mut ordered: Vec<(usize, ExactSeitz)> = Vec::with_capacity(transporters.len());
    ordered.push((identity_position, ExactSeitz::identity()));
    ordered.extend(
        transporters
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != identity_position)
            .map(|(index, transporter)| (index, *transporter)),
    );

    let mut collected: Vec<(usize, StarArm)> = Vec::new();
    for (input_index, transporter) in ordered {
        let wave_vector = arm_wave_vector(seed, &transporter)?;
        let mut existing: Option<usize> = None;
        for (first_index, arm) in &collected {
            if parent_reciprocal.same_mod(&arm.wave_vector, &wave_vector)? {
                existing = Some(*first_index);
                break;
            }
        }
        match existing {
            Some(first) if duplicates_are_errors => {
                return Err(StarError::DuplicateArm {
                    first,
                    second: input_index,
                });
            }
            Some(_) => continue,
            None => collected.push((
                input_index,
                StarArm {
                    wave_vector,
                    transporter,
                },
            )),
        }
    }
    // The identity transporter was first in `ordered`, so it is the first arm.
    let identity_arm = collected
        .first()
        .filter(|(index, _)| *index == identity_position)
        .map(|(_, arm)| *arm)
        .ok_or(StarError::MissingIdentityTransporter { sg: parent_sg })?;
    let mut keyed: Vec<_> = collected
        .into_iter()
        .skip(1)
        .map(|(_, arm)| (canonical_key(&arm.wave_vector), arm))
        .collect();
    keyed.sort_by_key(|entry| entry.0);
    let mut arms = Vec::with_capacity(keyed.len() + 1);
    arms.push(identity_arm);
    arms.extend(keyed.into_iter().map(|(_, arm)| arm));
    Ok(arms)
}

// ── The parent full star ─────────────────────────────────────────────────────

/// The full star of one ordinary scalar parent probe, with its induced
/// character.
#[derive(Debug, Clone)]
pub struct OrdinaryStar {
    parent_sg: u8,
    probe: &'static IrrepRecord,
    seed_k: Vec3R,
    selected_row: CharacterRow,
    arms: Vec<StarArm>,
    parent_lattice: Lattice,
    parent_reciprocal: Lattice,
    /// Parent Hall operations reduced modulo `L_G`: the membership test for
    /// arbitrary operations handed to [`Self::character`].
    parent_operations: Vec<ExactSeitz>,
    /// Complete parent star size times the selected-arm dimension.
    dimension: u32,
}

impl OrdinaryStar {
    /// Build the full star of a canonical ordinary scalar probe.
    ///
    /// The probe must be the generated table's own record (identity checked
    /// against [`query::irreps_of`]); compound rows and spinors are rejected
    /// with [`SubductionError::UnsupportedCharacterSpace`] **before** the typed
    /// row is requested, and the row itself must reproduce the record
    /// dimension: `chi(E)` has to equal `selected dimension × arms` and the
    /// record's `dim` at the same time, otherwise the constructor fails closed.
    pub fn new(probe: &'static IrrepRecord) -> Result<Self, StarError> {
        let (seed_k, parent_lattice, parent_reciprocal, selected_row) =
            Self::validate_probe(probe)?;
        let hall = strict_sg_hall_ops(probe.sg)?;
        let arms = collect_arms(
            probe.sg,
            &seed_k,
            &parent_reciprocal,
            &hall.operations,
            false,
        )?;
        Self::assemble(
            probe,
            seed_k,
            parent_lattice,
            parent_reciprocal,
            selected_row,
            arms,
            &hall.operations,
        )
    }

    /// Build the star from an explicit transversal (one transporter per arm).
    ///
    /// This is the reusable constructor behind [`Self::new`].  It revalidates
    /// everything a caller could inject: every transporter must be a parent
    /// group element modulo `L_G`, no two may select the same arm, and the
    /// resulting arms must cover the images of **all** parent operations, so a
    /// missing arm is rejected as [`StarError::IncompleteStar`] instead of
    /// producing a smaller star.
    pub fn from_transporters(
        probe: &'static IrrepRecord,
        transporters: &[ExactSeitz],
    ) -> Result<Self, StarError> {
        let (seed_k, parent_lattice, parent_reciprocal, selected_row) =
            Self::validate_probe(probe)?;
        let hall = strict_sg_hall_ops(probe.sg)?;
        let reduced = reduce_operations(&hall.operations, &parent_lattice)?;
        for (index, transporter) in transporters.iter().enumerate() {
            let image = transporter.reduce(&parent_lattice)?;
            if !reduced.contains(&image) {
                return Err(StarError::TransporterNotInParentGroup {
                    sg: probe.sg,
                    index,
                });
            }
        }
        let arms = collect_arms(probe.sg, &seed_k, &parent_reciprocal, transporters, true)?;
        Self::assemble(
            probe,
            seed_k,
            parent_lattice,
            parent_reciprocal,
            selected_row,
            arms,
            &hall.operations,
        )
    }

    /// Probe membership, character space and exact geometry.
    fn validate_probe(
        probe: &'static IrrepRecord,
    ) -> Result<(Vec3R, Lattice, Lattice, CharacterRow), StarError> {
        if !query::irreps_of(probe.sg)
            .iter()
            .any(|record| std::ptr::eq(record, probe))
        {
            return Err(SubductionError::ForeignProbe {
                sg: probe.sg,
                ml: probe.ml,
            }
            .into());
        }
        if probe.spinor || probe.compound_metadata().is_some() {
            return Err(SubductionError::UnsupportedCharacterSpace {
                sg: probe.sg,
                ml: probe.ml.to_string(),
            }
            .into());
        }
        let selected_row = probe
            .ordinary_scalar_selected_arm_block_trace()
            .map_err(|_| SubductionError::UnsupportedCharacterSpace {
                sg: probe.sg,
                ml: probe.ml.to_string(),
            })?;
        let seed_k = inline_k_vector(probe)?;
        let parent_lattice = Lattice::new(exact_primitive_basis(probe.sg)?)?;
        let parent_reciprocal = parent_lattice.reciprocal()?;
        Ok((seed_k, parent_lattice, parent_reciprocal, selected_row))
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble(
        probe: &'static IrrepRecord,
        seed_k: Vec3R,
        parent_lattice: Lattice,
        parent_reciprocal: Lattice,
        selected_row: CharacterRow,
        arms: Vec<StarArm>,
        hall_operations: &[ExactSeitz],
    ) -> Result<Self, StarError> {
        // Completeness: the arms are exactly the star of the seed, so every
        // parent operation has to land on one of them.
        for operation in hall_operations {
            let image = arm_wave_vector(&seed_k, operation)?;
            let mut covered = false;
            for arm in &arms {
                if parent_reciprocal.same_mod(&arm.wave_vector, &image)? {
                    covered = true;
                    break;
                }
            }
            if !covered {
                return Err(StarError::IncompleteStar {
                    sg: probe.sg,
                    ml: probe.ml,
                    represented: arms.len(),
                    operations: hall_operations.len(),
                    missing_k: Box::new([image.get(0), image.get(1), image.get(2)]),
                });
            }
        }
        let parent_operations = reduce_operations(hall_operations, &parent_lattice)?;
        let product = i64::try_from(selected_row.dimension())
            .ok()
            .and_then(|dimension| dimension.checked_mul(arms.len() as i64))
            .ok_or(SubductionError::RationalOverflow {
                operation: "full-star dimension",
            })?;
        let star = Self {
            parent_sg: probe.sg,
            probe,
            seed_k,
            selected_row,
            arms,
            parent_lattice,
            parent_reciprocal,
            parent_operations,
            dimension: u32::try_from(product).map_err(|_| SubductionError::RationalOverflow {
                operation: "full-star dimension",
            })?,
        };
        let identity_trace = star.character(&ExactSeitz::identity())?;
        let expected = Complex64::new(product as f64, 0.0);
        if (identity_trace - expected).norm() > SUBDUCTION_TOLERANCE
            || i64::from(probe.dim) != product
        {
            return Err(StarError::StarDimensionMismatch {
                ml: probe.ml,
                probe_dim: probe.dim,
                selected_dimension: star.selected_row.dimension(),
                arms: star.arms.len(),
                product,
                identity_trace,
            });
        }
        Ok(star)
    }

    /// Parent space group number.
    pub const fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// The canonical probe record this star belongs to.
    pub const fn probe(&self) -> &'static IrrepRecord {
        self.probe
    }

    /// Miller-Love label of the probe.
    pub const fn probe_ml(&self) -> &'static str {
        self.probe.ml
    }

    /// The selected arm's exact stored wave vector.
    pub const fn seed_k(&self) -> &Vec3R {
        &self.seed_k
    }

    /// The selected-arm character row the induction evaluates.
    pub const fn selected_row(&self) -> &CharacterRow {
        &self.selected_row
    }

    /// One arm per little-group coset, identity arm first.
    pub fn arms(&self) -> &[StarArm] {
        &self.arms
    }

    /// Number of arms: the index of the little group in the parent group.
    pub fn arm_count(&self) -> usize {
        self.arms.len()
    }

    /// Full-star dimension `chi(E)`, verified against the record dimension.
    pub const fn dimension(&self) -> u32 {
        self.dimension
    }

    /// The induced full-star character of one parent operation.
    ///
    /// Only arms fixed by `operation` modulo the parent reciprocal lattice
    /// contribute; the others contribute exactly zero.  A fixed arm
    /// contributes `chi_seed(g_i^-1 h g_i)`, evaluated on the exact conjugated
    /// Seitz operation with the **seed** wave vector so the stored Bloch
    /// phases are preserved.
    pub fn character(&self, operation: &ExactSeitz) -> Result<Complex64, StarError> {
        let reduced = operation.reduce(&self.parent_lattice)?;
        if !self.parent_operations.contains(&reduced) {
            return Err(StarError::OperationNotInParentGroup { sg: self.parent_sg });
        }
        induced_character(
            &self.arms,
            &LittleCharacter::Stored {
                row: &self.selected_row,
                ml: self.probe.ml,
                row_k: &self.seed_k,
            },
            &self.parent_lattice,
            &self.parent_reciprocal,
            operation,
        )
    }

    /// Fold every arm into the child frame and partition the folded points into
    /// child stars.
    ///
    /// Grouping is exact: equal folded `q` is equality modulo the **child**
    /// reciprocal lattice built from the child's own primitive basis, and each
    /// group keeps every parent arm index that folded onto it.  The `q` points
    /// are then partitioned into orbits of the child point group pulled back
    /// through `T` (rotations only, so the frozen child origin shift plays no
    /// role); every orbit is checked for closure and for a uniform arm count.
    ///
    /// The result reports geometry only.  It deliberately performs no character
    /// aggregation: the child multiplicities need the decomposition stage,
    /// which this adapter does not implement.
    pub fn folded_stars(
        &self,
        embedding: &SubgroupEmbedding,
    ) -> Result<Vec<FoldedStar>, StarError> {
        let dimension = self.selected_row.dimension();
        let arms: Vec<FoldArm> = self
            .arms
            .iter()
            .map(|arm| FoldArm {
                wave_vector: *arm.wave_vector(),
                dimension,
            })
            .collect();
        fold_arms(embedding, self.parent_sg, &arms)
    }
}

/// One parent arm as the folded geometry sees it: its exact wave vector and the
/// selected-arm dimension of the component it belongs to.
///
/// The folded geometry is shared by [`OrdinaryStar`] (one component, one
/// dimension) and the scalar adapter (several components whose frozen selected
/// dimensions are all equal); a mixed-dimension arm list fails closed with
/// [`StarError::NonUniformComponentDimension`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct FoldArm {
    wave_vector: Vec3R,
    dimension: usize,
}

/// Fold one arm list into the child frame and partition it into child stars.
///
/// This is the geometry both star adapters share; see
/// [`OrdinaryStar::folded_stars`] for the exactness and grouping contract.
pub(super) fn fold_arms(
    embedding: &SubgroupEmbedding,
    parent_sg: u8,
    arms: &[FoldArm],
) -> Result<Vec<FoldedStar>, StarError> {
    if embedding.parent_sg() != parent_sg {
        return Err(StarError::EmbeddingParentMismatch {
            star_parent_sg: parent_sg,
            embedding_parent_sg: embedding.parent_sg(),
        });
    }
    let selected_dimension = arms
        .first()
        .map(|arm| arm.dimension)
        .ok_or(StarError::FoldLostArms { folded: 0, arms: 0 })?;
    if let Some(other) = arms.iter().find(|arm| arm.dimension != selected_dimension) {
        return Err(StarError::NonUniformComponentDimension {
            first: selected_dimension,
            other: other.dimension,
        });
    }
    let child_cell = Lattice::new(exact_primitive_basis(embedding.subgroup_sg())?)?;
    let child_reciprocal = child_cell.reciprocal()?;

    // Fold every arm and merge equal q, preserving all arm indices.
    let mut points: Vec<FoldedPoint> = Vec::new();
    for (arm_index, arm) in arms.iter().enumerate() {
        let q = fold_wave_vector(embedding.transform(), &arm.wave_vector)?;
        let mut found: Option<usize> = None;
        for (slot, point) in points.iter().enumerate() {
            if child_reciprocal.same_mod(&point.q, &q)? {
                found = Some(slot);
                break;
            }
        }
        match found {
            Some(slot) => points[slot].arm_indices.push(arm_index),
            None => points.push(FoldedPoint {
                q,
                arm_indices: vec![arm_index],
            }),
        }
    }

    // Child rotations pulled back through T: rotations are unaffected by
    // the frozen child shift, which is a translation.
    let mut rotations: Vec<Mat3I> = Vec::new();
    for operation in embedding.representatives() {
        let child = embedding.transform().unmap_operation(operation)?;
        if !rotations.contains(&child.rotation()) {
            rotations.push(child.rotation());
        }
    }

    // Orbit partition with closure checked at every image.
    let mut assignment: Vec<Option<usize>> = vec![None; points.len()];
    let mut orbits: Vec<Vec<usize>> = Vec::new();
    for start in 0..points.len() {
        if assignment[start].is_some() {
            continue;
        }
        let star_index = orbits.len();
        let mut orbit = vec![start];
        assignment[start] = Some(star_index);
        let mut cursor = 0;
        while cursor < orbit.len() {
            let current = orbit[cursor];
            cursor += 1;
            let q = points[current].q;
            for rotation in &rotations {
                let image = Mat3R::from_ints(*rotation)
                    .inverse()?
                    .transpose()
                    .checked_mul_vector(&q)?;
                let mut target: Option<usize> = None;
                for (slot, point) in points.iter().enumerate() {
                    if child_reciprocal.same_mod(&point.q, &image)? {
                        target = Some(slot);
                        break;
                    }
                }
                let target = target.ok_or_else(|| StarError::ChildStarNotClosed {
                    q: [q.get(0), q.get(1), q.get(2)],
                })?;
                match assignment[target] {
                    Some(existing) if existing != star_index => {
                        return Err(StarError::ChildStarNotClosed {
                            q: [q.get(0), q.get(1), q.get(2)],
                        });
                    }
                    Some(_) => {}
                    None => {
                        assignment[target] = Some(star_index);
                        orbit.push(target);
                    }
                }
            }
        }
        orbits.push(orbit);
    }

    if assignment.iter().any(Option::is_none) {
        return Err(StarError::ChildStarNotClosed {
            q: [Rat::ZERO, Rat::ZERO, Rat::ZERO],
        });
    }

    // Canonically order the points inside every orbit, check the uniform
    // arm block and derive the block dimension.
    let mut stars: Vec<(Vec<usize>, usize, u32)> = Vec::new();
    let mut folded_arms = 0usize;
    for orbit in orbits {
        let mut order: Vec<_> = orbit
            .into_iter()
            .map(|slot| (canonical_key(&points[slot].q), slot))
            .collect();
        order.sort_by_key(|entry| entry.0);
        let ordered: Vec<usize> = order.into_iter().map(|(_, slot)| slot).collect();
        let arms_per_point = points[ordered[0]].arm_indices.len();
        for slot in &ordered {
            if points[*slot].arm_indices.len() != arms_per_point {
                return Err(StarError::NonUniformArmBlock {
                    arms: arms_per_point,
                    other_arms: points[*slot].arm_indices.len(),
                });
            }
        }
        let arm_count = arms_per_point * ordered.len();
        folded_arms += arm_count;
        let block_dimension = u32::try_from(arm_count * selected_dimension).map_err(|_| {
            SubductionError::RationalOverflow {
                operation: "folded star dimension",
            }
        })?;
        stars.push((ordered, arm_count, block_dimension));
    }
    if folded_arms != arms.len() {
        return Err(StarError::FoldLostArms {
            folded: folded_arms,
            arms: arms.len(),
        });
    }

    // Canonical star order uses the same exact encoding as point order.
    let mut keyed: Vec<_> = stars
        .into_iter()
        .map(|(ordered, arm_count, block_dimension)| {
            (
                canonical_key(&points[ordered[0]].q),
                ordered,
                arm_count,
                block_dimension,
            )
        })
        .collect();
    keyed.sort_by_key(|entry| entry.0);
    Ok(keyed
        .into_iter()
        .map(|(_, ordered, arm_count, block_dimension)| FoldedStar {
            points: ordered
                .into_iter()
                .map(|slot| points[slot].clone())
                .collect(),
            seed_dimension: selected_dimension,
            arm_count,
            block_dimension,
        })
        .collect())
}

// ── Little-group character sources ───────────────────────────────────────────

/// A little-group representation computed from the subgroup operations and an
/// exact folded point, with no pinned table row behind it.
///
/// The character is a function of the little-group operation itself, so it can
/// be transported to another arm of the same star by conjugating the operation,
/// exactly like a stored row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum ConstructedLittleRep {
    /// Trivial little co-group: every little-group operation is a translation,
    /// so the representation at the exact point `q` is the one-dimensional
    /// Bloch phase `D(T_L) = exp(+2 pi i q.L)` in the child's own frame (the
    /// sign convention [`bloch_phase`] already fixes for stored rows).
    BlochPhase { q: Vec3R },
    /// One one-dimensional allowed irrep of a **non-trivial** little co-group:
    /// `D(R, T) = exp(2 pi i (constant(R) + q.T))`.
    ///
    /// `q` is the same exact point the constants were built with, which is what
    /// makes the two halves of the phase agree: a reconstructed little-group
    /// operation is **not** a lattice translate of its representative (its
    /// translation can carry quarters), so mixing a raw and a reduced `q`
    /// between the constants and the evaluation silently changes the character.
    Projective {
        /// The exact folded point the catalogue belongs to.
        q: Vec3R,
        /// `(rotation, constant)` per little co-group rotation.
        constants: Vec<(Mat3I, Rat)>,
    },
}

impl ConstructedLittleRep {
    /// Character of one little-group operation.
    ///
    /// A contribution is only ever asked for an operation that fixes the arm it
    /// is conjugated into, so a non-identity rotation means the caller used the
    /// wrong little group: that fails closed instead of answering a phase.
    pub(super) fn character(&self, operation: &ExactSeitz) -> Result<Complex64, StarError> {
        match self {
            Self::BlochPhase { q } => {
                if operation.rotation() != IDENTITY_ROTATION {
                    return Err(StarError::ConstructedRotationNotCovered {
                        q: [q.get(0), q.get(1), q.get(2)],
                    });
                }
                Ok(bloch_phase(q, operation.translation())?)
            }
            Self::Projective { q, constants } => {
                Ok(catalogue::character_value(constants, q, operation)?)
            }
        }
    }
}

/// The little-group character of one arm's conjugated operation.
///
/// A stored row is looked up by its own Seitz representatives with the **stored**
/// wave vector, so every pinned Bloch phase survives; a constructed rep answers
/// from the little-group operation itself. Both keep the arm's own fixity test
/// and transporter conjugation in [`arm_character`].
pub(super) enum LittleCharacter<'a> {
    /// A pinned row, evaluated at the wave vector the row belongs to.
    Stored {
        row: &'a CharacterRow,
        ml: &'static str,
        row_k: &'a Vec3R,
    },
    /// A representation built at an exact folded point.
    Constructed(&'a ConstructedLittleRep),
}

impl LittleCharacter<'_> {
    fn at(
        &self,
        operation: &ExactSeitz,
        lattice: &Lattice,
        index: usize,
    ) -> Result<Complex64, StarError> {
        match self {
            Self::Stored { row, ml, row_k } => {
                Ok(character_of(row, ml, operation, lattice, row_k, index)?)
            }
            Self::Constructed(rep) => rep.character(operation),
        }
    }
}

/// The induced trace: fixed arms contribute their conjugated seed character,
/// moved arms contribute zero.
fn induced_character(
    arms: &[StarArm],
    character: &LittleCharacter<'_>,
    parent_lattice: &Lattice,
    parent_reciprocal: &Lattice,
    operation: &ExactSeitz,
) -> Result<Complex64, StarError> {
    induced_component_character(
        arms,
        character,
        parent_lattice,
        parent_reciprocal,
        false,
        operation,
    )
}

/// The induced trace of one complex component.
///
/// The arm fixity test uses the arms' own (effective) wave vectors, while the
/// row is always evaluated at `row_k`, the wave vector the stored row belongs
/// to.  `conjugate` then conjugates the whole sum, which is the character of
/// the complex-conjugate component: its Bloch phase is the `-row_k` phase, not
/// the `+row_k` phase, and conjugating the evaluated value (phase included)
/// produces exactly that.
#[allow(clippy::too_many_arguments)]
fn induced_component_character(
    arms: &[StarArm],
    character: &LittleCharacter<'_>,
    lattice: &Lattice,
    reciprocal: &Lattice,
    conjugate: bool,
    operation: &ExactSeitz,
) -> Result<Complex64, StarError> {
    let mut total = Complex64::new(0.0, 0.0);
    for (index, arm) in arms.iter().enumerate() {
        total += arm_character(arm, character, lattice, reciprocal, conjugate, index, operation)?;
    }
    Ok(total)
}

/// One arm's contribution to an induced trace: zero unless `operation` fixes
/// the arm modulo the reciprocal lattice.
fn arm_character(
    arm: &StarArm,
    character: &LittleCharacter<'_>,
    lattice: &Lattice,
    reciprocal: &Lattice,
    conjugate: bool,
    index: usize,
    operation: &ExactSeitz,
) -> Result<Complex64, StarError> {
    if !reciprocal.preserves(operation.rotation(), &arm.wave_vector)? {
        return Ok(Complex64::new(0.0, 0.0));
    }
    let conjugated = arm
        .transporter
        .inverse()?
        .compose(operation)?
        .compose(&arm.transporter)?;
    let value = character.at(&conjugated, lattice, index)?;
    Ok(if conjugate { value.conj() } else { value })
}

// ── Folded child stars ───────────────────────────────────────────────────────

/// One folded child wave vector with every parent arm that folded onto it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FoldedPoint {
    q: Vec3R,
    arm_indices: Vec<usize>,
}

impl FoldedPoint {
    /// The exact folded wave vector `T^T k_i` (unreduced).
    pub const fn q(&self) -> &Vec3R {
        &self.q
    }

    /// Indices into the source adapter's arm list, ascending. For a scalar
    /// compound, coincident arms retain separate component indices.
    pub fn arm_indices(&self) -> &[usize] {
        &self.arm_indices
    }

    /// Number of parent component arms folded onto this `q`.
    pub fn arm_count(&self) -> usize {
        self.arm_indices.len()
    }

    /// Build one folded point from a folded wave vector and its parent arms.
    ///
    /// The parametric-k line source folds the arms of a line rather than the
    /// star of a discrete `k`, so it has to build these directly; see
    /// `docs/task9-remaining-work.md` for the data shape.
    #[allow(dead_code)]
    pub(crate) fn from_parts(q: Vec3R, arm_indices: Vec<usize>) -> Self {
        Self { q, arm_indices }
    }
}

/// One star of the child group: an orbit of folded `q` points under the child
/// point group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FoldedStar {
    points: Vec<FoldedPoint>,
    seed_dimension: usize,
    arm_count: usize,
    block_dimension: u32,
}

impl FoldedStar {
    /// The orbit's points, canonically ordered by `q`.
    pub fn points(&self) -> &[FoldedPoint] {
        &self.points
    }

    /// Number of distinct folded `q` points in the orbit.
    pub fn star_size(&self) -> usize {
        self.points.len()
    }

    /// Total number of parent component arms in the orbit, including coincident
    /// arms of different components.
    pub fn arm_count(&self) -> usize {
        self.arm_count
    }

    /// Common selected-arm dimension of one complex parent component.
    pub const fn seed_dimension(&self) -> usize {
        self.seed_dimension
    }

    /// Dimension of the parent full-star block carried by this child star:
    /// `arm_count × seed_dimension`.
    pub const fn block_dimension(&self) -> u32 {
        self.block_dimension
    }

    /// Build one folded child star from its points and dimensions.
    ///
    /// Same reason as [`FoldedPoint::from_parts`]: the parametric-k line source
    /// produces folded stars that no `ScalarStar` ever builds.
    #[allow(dead_code)]
    pub(crate) fn from_parts(
        points: Vec<FoldedPoint>,
        seed_dimension: usize,
        arm_count: usize,
        block_dimension: u32,
    ) -> Self {
        Self {
            points,
            seed_dimension,
            arm_count,
            block_dimension,
        }
    }

    /// Every parent arm index of the orbit, ascending.
    pub fn parent_arm_indices(&self) -> Vec<usize> {
        let mut out: Vec<usize> = self
            .points
            .iter()
            .flat_map(|point| point.arm_indices.iter().copied())
            .collect();
        out.sort_unstable();
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::LabelConvention;
    use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
    use crate::irrep::subduction::IDENTITY_SETTING;

    fn pair(num: i128, den: i128) -> Rat {
        Rat::new(num, den).expect("non-zero denominator")
    }

    fn probe(sg: u8, ml: &str) -> &'static IrrepRecord {
        query::irreps_of(sg)
            .iter()
            .find(|record| record.ml == ml && !record.spinor)
            .unwrap_or_else(|| panic!("SG {sg} has no scalar irrep {ml}"))
    }

    fn embedding(parent: u8, ml: &str, direction: &str) -> SubgroupEmbedding {
        let subgroup =
            isotropy_subgroup_for_direction(
                parent,
                ml,
                LabelConvention::Cdml,
                IsotropyDirection::Label(direction),
            )
                .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction}: {error}"));
        SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
            .unwrap_or_else(|error| panic!("SG {parent} {ml} {direction} embedding: {error}"))
    }

    fn identity_trace(star: &OrdinaryStar) -> Complex64 {
        star.character(&ExactSeitz::identity())
            .expect("identity is a parent operation")
    }

    #[test]
    fn every_ordinary_scalar_probe_of_the_frozen_parents_folds_completely() {
        // The four frozen (parent, condensing irrep, direction) contexts of
        // tasks 4-7, plus the second 221 direction.  Every ordinary scalar
        // probe of the parent must build its complete star and fold all of its
        // arms into child stars whose block dimensions tile the star.
        let contexts = [
            (221, "GM4+", "P1", 40),
            (221, "GM4+", "P2", 40),
            (16, "R1", "P1", 32),
            (167, "GM3+", "P1", 12),
            (139, "M1-", "P1", 37),
        ];
        for (parent, condensing, direction, expected) in contexts {
            let embedding = embedding(parent, condensing, direction);
            let mut folded = 0usize;
            for record in query::irreps_of(parent) {
                if record.spinor || record.compound_metadata().is_some() {
                    continue;
                }
                let star = OrdinaryStar::new(record).unwrap_or_else(|error| {
                    panic!(
                        "SG {parent} {}: star construction failed: {error}",
                        record.ml
                    )
                });
                assert_eq!(u32::from(record.dim), star.dimension());
                let stars = star.folded_stars(&embedding).unwrap_or_else(|error| {
                    panic!("SG {parent} {}: folding failed: {error}", record.ml)
                });
                let arm_total: usize = stars.iter().map(FoldedStar::arm_count).sum();
                let block_total: u32 = stars.iter().map(FoldedStar::block_dimension).sum();
                assert_eq!(arm_total, star.arm_count(), "SG {parent} {}", record.ml);
                assert_eq!(block_total, star.dimension(), "SG {parent} {}", record.ml);
                folded += 1;
            }
            assert_eq!(
                folded, expected,
                "context {parent} {condensing} {direction} changed probe count"
            );
        }
    }

    #[test]
    fn canonical_order_keeps_rational_denominators() {
        assert_ne!(
            canonical_key(&Vec3R::new([pair(1, 2), Rat::ZERO, Rat::ZERO])),
            canonical_key(&Vec3R::from_ints([1, 0, 0])),
        );
    }

    /// Two parent arms fold onto **one** child `q`.
    ///
    /// Pinned from the frozen (139 `M1-` `P1`, ordinal 6213) embedding: the two
    /// arms of `X1+` (`k = (1/2,1/2,0)`) both fold to `q = (1/2,1/2,0)` modulo
    /// the child reciprocal lattice, so the whole parent star is one child star
    /// with a two-arm block.
    #[test]
    fn several_parent_arms_fold_onto_the_same_child_q() {
        let star = OrdinaryStar::new(probe(139, "X1+")).expect("X1+ star");
        assert_eq!(star.arm_count(), 2);
        assert_eq!(star.dimension(), 2);
        assert_eq!(star.selected_row().dimension(), 1);
        let embedding = embedding(139, "M1-", "P1");
        assert_eq!(embedding.ordinal(), 6213);
        assert_eq!(embedding.subgroup_sg(), 126);

        let stars = star.folded_stars(&embedding).expect("folded");
        assert_eq!(stars.len(), 1);
        let folded = &stars[0];
        assert_eq!(folded.star_size(), 1);
        assert_eq!(folded.arm_count(), 2);
        assert_eq!(folded.seed_dimension(), 1);
        assert_eq!(folded.block_dimension(), 2);
        assert_eq!(folded.parent_arm_indices(), [0, 1]);
        assert_eq!(folded.points().len(), 1);
        assert_eq!(
            *folded.points()[0].q(),
            Vec3R::new([pair(1, 2), pair(1, 2), Rat::ZERO])
        );
        assert_eq!(folded.points()[0].arm_indices(), [0, 1]);
    }

    /// A richer merge: 139 `N1+` has four arms that fold to **two** child `q`
    /// points, two arms each, inside a single child star (ordinal 6213 again).
    #[test]
    fn merged_arms_group_by_folded_q_and_stay_in_one_star() {
        let star = OrdinaryStar::new(probe(139, "N1+")).expect("N1+ star");
        assert_eq!(star.arm_count(), 4);
        assert_eq!(star.dimension(), 4);
        let stars = star
            .folded_stars(&embedding(139, "M1-", "P1"))
            .expect("folded");
        assert_eq!(stars.len(), 1);
        assert_eq!(stars[0].star_size(), 2);
        assert_eq!(stars[0].arm_count(), 4);
        assert_eq!(stars[0].block_dimension(), 4);
        let points: Vec<([Rat; 3], Vec<usize>)> = stars[0]
            .points()
            .iter()
            .map(|point| {
                (
                    [point.q().get(0), point.q().get(1), point.q().get(2)],
                    point.arm_indices().to_vec(),
                )
            })
            .collect();
        assert_eq!(
            points,
            [
                ([Rat::ZERO, pair(-1, 2), pair(1, 2)], vec![2, 3]),
                ([pair(1, 2), Rat::ZERO, pair(1, 2)], vec![0, 1]),
            ],
            "canonical q order, every parent arm preserved"
        );
    }

    /// Three parent arms split into **two** child stars.
    ///
    /// Pinned from the frozen (221 `GM4+` `P1`, ordinal 12400) embedding: the
    /// three `X1+` arms (one selected-arm dimension each) fold to two distinct
    /// child `q` points in one star plus a third `q` in a second star, with
    /// block dimensions 2 and 1.
    #[test]
    fn parent_arms_split_into_several_child_stars() {
        let star = OrdinaryStar::new(probe(221, "X1+")).expect("X1+ star");
        assert_eq!(star.arm_count(), 3);
        assert_eq!(star.dimension(), 3);
        assert_eq!(star.selected_row().dimension(), 1);
        let embedding = embedding(221, "GM4+", "P1");
        assert_eq!(embedding.ordinal(), 12400);
        assert_eq!(embedding.subgroup_sg(), 83);

        let stars = star.folded_stars(&embedding).expect("folded");
        let shape: Vec<(usize, usize, u32)> = stars
            .iter()
            .map(|star| (star.star_size(), star.arm_count(), star.block_dimension()))
            .collect();
        assert_eq!(shape, [(2, 2, 2), (1, 1, 1)]);
        let points: Vec<([Rat; 3], Vec<usize>)> = stars
            .iter()
            .flat_map(|star| star.points())
            .map(|point| {
                (
                    [point.q().get(0), point.q().get(1), point.q().get(2)],
                    point.arm_indices().to_vec(),
                )
            })
            .collect();
        assert_eq!(
            points,
            [
                ([Rat::ZERO, pair(-1, 2), Rat::ZERO], vec![0]),
                ([pair(1, 2), Rat::ZERO, Rat::ZERO], vec![2]),
                ([Rat::ZERO, Rat::ZERO, pair(-1, 2)], vec![1]),
            ]
        );
        // Block dimensions tile the full-star dimension; nothing is lost.
        assert_eq!(
            stars.iter().map(FoldedStar::block_dimension).sum::<u32>(),
            star.dimension()
        );
        assert_eq!(
            stars.iter().map(FoldedStar::arm_count).sum::<usize>(),
            star.arm_count()
        );
    }

    /// The same split in the second frozen parent, 167 `GM3+` `P1` ->
    /// #15 (ordinal 7651), whose folded `q` values include `(1,0,0)`.
    #[test]
    fn second_parent_arms_split_into_several_child_stars() {
        let star = OrdinaryStar::new(probe(167, "F1+")).expect("F1+ star");
        assert_eq!(star.arm_count(), 3);
        assert_eq!(star.dimension(), 3);
        let embedding = embedding(167, "GM3+", "P1");
        assert_eq!(embedding.ordinal(), 7651);
        assert_eq!(embedding.subgroup_sg(), 15);
        let stars = star.folded_stars(&embedding).expect("folded");
        let points: Vec<([Rat; 3], Vec<usize>)> = stars
            .iter()
            .flat_map(|star| star.points())
            .map(|point| {
                (
                    [point.q().get(0), point.q().get(1), point.q().get(2)],
                    point.arm_indices().to_vec(),
                )
            })
            .collect();
        assert_eq!(
            points,
            [
                ([Rat::from_integer(1), Rat::ZERO, Rat::ZERO], vec![0]),
                ([pair(1, 2), pair(-1, 2), pair(1, 2)], vec![2]),
                ([pair(1, 2), pair(1, 2), pair(1, 2)], vec![1]),
            ]
        );
    }

    #[test]
    fn folded_stars_reject_a_foreign_parent_embedding() {
        let star = OrdinaryStar::new(probe(221, "X1+")).expect("star");
        let foreign = embedding(167, "GM3+", "P1");
        assert!(matches!(
            star.folded_stars(&foreign),
            Err(StarError::EmbeddingParentMismatch {
                star_parent_sg: 221,
                embedding_parent_sg: 167,
            })
        ));
    }

    /// The operation with this rotation, exactly as the data Hall table stores
    /// it (translations included).
    fn hall_operation(parent: u8, rotation: Mat3I) -> ExactSeitz {
        strict_sg_hall_ops(parent)
            .expect("Hall operations")
            .operations
            .iter()
            .copied()
            .find(|operation| operation.rotation() == rotation)
            .unwrap_or_else(|| panic!("SG {parent} has no operation with rotation {rotation:?}"))
    }

    fn assert_close(found: Complex64, expected: f64) {
        assert!(
            (found - Complex64::new(expected, 0.0)).norm() <= SUBDUCTION_TOLERANCE,
            "{found} != {expected}"
        );
    }

    /// Only arms fixed modulo `L_G*` contribute; every moved arm contributes
    /// exactly zero.
    ///
    /// In SG 221, `R4z` fixes only the `(0,0,1/2)` arm of `X1+` (the other two
    /// evaluate to zero), so `chi = 1`.  In SG 139 the same rotation swaps the
    /// two `X1+` arms, so no arm is fixed and `chi = 0`.
    #[test]
    fn only_fixed_arms_contribute_to_the_induced_trace() {
        let r4z = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];

        let star = OrdinaryStar::new(probe(221, "X1+")).expect("star");
        let operation = hall_operation(221, r4z);
        let mut fixed = 0usize;
        for (index, arm) in star.arms().iter().enumerate() {
            let conjugated = arm
                .transporter()
                .inverse()
                .expect("inverse")
                .compose(&operation)
                .expect("compose")
                .compose(arm.transporter())
                .expect("compose");
            let value = character_of(
                star.selected_row(),
                star.probe.ml,
                &conjugated,
                &star.parent_lattice,
                &star.seed_k,
                index,
            )
            .expect("conjugated character");
            if star
                .parent_reciprocal
                .preserves(r4z, arm.wave_vector())
                .expect("preserves")
            {
                fixed += 1;
                assert_close(value, 1.0);
            } else {
                assert_close(value, 0.0);
            }
        }
        assert_eq!(fixed, 1);
        assert_close(star.character(&operation).expect("character"), 1.0);

        let star = OrdinaryStar::new(probe(139, "X1+")).expect("star");
        let operation = hall_operation(139, r4z);
        let fixed = star
            .arms()
            .iter()
            .filter(|arm| {
                star.parent_reciprocal
                    .preserves(r4z, arm.wave_vector())
                    .expect("preserves")
            })
            .count();
        assert_eq!(fixed, 0);
        assert_close(star.character(&operation).expect("character"), 0.0);
    }

    /// The identity transporter is first, the remaining arms are canonically
    /// ordered by their exact wave vector, and the images are stored
    /// unreduced: `R4z` sends `(0,1/2,0)` to `(-1/2,0,0)`, not to `(1/2,0,0)`.
    #[test]
    fn arms_are_identity_first_canonical_and_unreduced() {
        let star = OrdinaryStar::new(probe(221, "X1+")).expect("star");
        assert_eq!(star.arm_count(), 3);
        assert_eq!(star.dimension(), 3);
        assert_eq!(u32::from(probe(221, "X1+").dim), 3);
        let arms: Vec<(Vec3R, Mat3I)> = star
            .arms()
            .iter()
            .map(|arm| (*arm.wave_vector(), arm.transporter().rotation()))
            .collect();
        assert_eq!(
            arms,
            [
                (
                    Vec3R::new([Rat::ZERO, pair(1, 2), Rat::ZERO]),
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                ),
                (
                    Vec3R::new([pair(-1, 2), Rat::ZERO, Rat::ZERO]),
                    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                ),
                (
                    Vec3R::new([Rat::ZERO, Rat::ZERO, pair(1, 2)]),
                    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
                ),
            ]
        );
        assert_eq!(star.arms()[0].transporter(), &ExactSeitz::identity());
        assert_eq!(*star.arms()[0].wave_vector(), *star.seed_k());
    }

    /// Reordering the transversal changes neither the canonical arm set, nor
    /// any induced character, nor the folded child-star partition.
    #[test]
    fn arm_permutation_preserves_trace_and_folded_groups() {
        let embedding = embedding(221, "GM4+", "P1");
        let star = OrdinaryStar::new(probe(221, "X1+")).expect("star");
        let mut transporters: Vec<ExactSeitz> =
            star.arms().iter().map(|arm| *arm.transporter()).collect();
        transporters.reverse();
        let permuted =
            OrdinaryStar::from_transporters(probe(221, "X1+"), &transporters).expect("permuted");
        assert_eq!(permuted.arms(), star.arms(), "canonical arm order");
        assert_eq!(permuted.dimension(), star.dimension());
        for operation in strict_sg_hall_ops(221).expect("hall").operations {
            let expected = star.character(&operation).expect("character");
            let found = permuted.character(&operation).expect("character");
            assert!(
                (found - expected).norm() <= SUBDUCTION_TOLERANCE,
                "{operation:?}: {found} != {expected}"
            );
        }
        // Both sides are emitted in canonical star/point order, so the
        // structural comparison is the canonical comparison.
        assert_eq!(
            permuted.folded_stars(&embedding).expect("folded"),
            star.folded_stars(&embedding).expect("folded")
        );
    }

    /// Removing one arm from the transversal is rejected by the completeness
    /// check, not silently folded into a smaller star.  The removed
    /// transporter is itself a valid parent operation, so the failure is the
    /// missing coverage of the parent star.
    #[test]
    fn a_missing_arm_fails_completeness() {
        for (sg, ml) in [(221u8, "X1+"), (139, "X1+")] {
            let star = OrdinaryStar::new(probe(sg, ml)).expect("star");
            assert!(star.arm_count() > 1);
            let mut transporters: Vec<ExactSeitz> =
                star.arms().iter().map(|arm| *arm.transporter()).collect();
            let removed = transporters.remove(1);
            assert!(
                star.parent_operations
                    .contains(&removed.reduce(&star.parent_lattice).expect("reduce")),
                "the removed transporter must be a parent operation"
            );
            let expected_operations = strict_sg_hall_ops(sg).expect("hall").operations.len();
            match OrdinaryStar::from_transporters(probe(sg, ml), &transporters) {
                Err(StarError::IncompleteStar {
                    sg: found_sg,
                    ml: found_ml,
                    represented,
                    operations,
                    ..
                }) => {
                    assert_eq!(found_sg, sg);
                    assert_eq!(found_ml, ml);
                    assert_eq!(represented, star.arm_count() - 1);
                    assert_eq!(operations, expected_operations);
                }
                other => panic!("expected IncompleteStar, got {other:?}"),
            }
        }
    }

    /// A transport list entry that is not a parent operation modulo `L_G` is
    /// rejected before any arm is built.
    #[test]
    fn a_non_parent_transporter_is_rejected() {
        let forged = ExactSeitz::new(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            Vec3R::new([Rat::ZERO, Rat::ZERO, pair(1, 2)]),
        );
        let transporters = [ExactSeitz::identity(), forged];
        match OrdinaryStar::from_transporters(probe(221, "X1+"), &transporters) {
            Err(StarError::TransporterNotInParentGroup { sg, index }) => {
                assert_eq!(sg, 221);
                assert_eq!(index, 1);
            }
            other => panic!("expected TransporterNotInParentGroup, got {other:?}"),
        }
    }

    #[test]
    fn gamma_star_has_one_arm_and_the_record_dimension() {
        let star = OrdinaryStar::new(probe(221, "GM4+")).expect("gamma star");
        assert_eq!(star.arm_count(), 1);
        assert_eq!(star.dimension(), 3);
        assert_eq!(u32::from(probe(221, "GM4+").dim), 3);
        assert_close(identity_trace(&star), 3.0);
        assert_eq!(*star.seed_k(), Vec3R::zero());
        assert_eq!(
            *star.arms()[0].wave_vector(),
            Vec3R::zero(),
            "the identity arm keeps the seed k exactly"
        );
        assert_eq!(*star.arms()[0].transporter(), ExactSeitz::identity());
    }

    #[test]
    fn compound_and_spinor_probes_fail_unsupported() {
        // SG 83 GM3+GM4+ is a `DistinctComponentSum` compound row.
        let compound = query::irreps_of(83)
            .iter()
            .find(|record| record.ml == "GM3+GM4+")
            .expect("pinned compound record");
        assert!(compound.compound_metadata().is_some());
        assert!(matches!(
            OrdinaryStar::new(compound),
            Err(StarError::Subduction(
                SubductionError::UnsupportedCharacterSpace { sg: 83, .. }
            ))
        ));

        let spinor = query::irreps_of(221)
            .iter()
            .find(|record| record.spinor)
            .expect("a spinor record in SG 221");
        assert!(matches!(
            OrdinaryStar::new(spinor),
            Err(StarError::Subduction(
                SubductionError::UnsupportedCharacterSpace { sg: 221, .. }
            ))
        ));
    }

    #[test]
    fn foreign_probe_is_rejected_by_identity() {
        // A hand-built copy is a different object even with identical fields.
        let original = probe(221, "GM4+");
        let copy = IrrepRecord { ..*original };
        assert!(matches!(
            OrdinaryStar::new(Box::leak(Box::new(copy))),
            Err(StarError::Subduction(SubductionError::ForeignProbe {
                sg: 221,
                ..
            }))
        ));
    }

    #[test]
    fn character_rejects_non_parent_operations() {
        let star = OrdinaryStar::new(probe(221, "GM4+")).expect("gamma star");
        let forged = ExactSeitz::new(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            Vec3R::new([pair(0, 1), pair(0, 1), pair(1, 2)]),
        );
        assert!(matches!(
            star.character(&forged),
            Err(StarError::OperationNotInParentGroup { sg: 221 })
        ));
    }

    #[test]
    fn probe_membership_survives_the_table_lookup() {
        // `query::irreps_of` hands out the very reference the star validates,
        // and the star keeps it instead of copying it.
        let record = probe(167, "GM3+");
        let star = OrdinaryStar::new(record).expect("gamma star");
        assert_eq!(star.probe_ml(), "GM3+");
        assert!(std::ptr::eq(star.probe(), record));
    }

    /// The four scalar self-restriction embeddings of task 8c are frozen with
    /// the identity setting and no shift, and their mapped subgroup operations
    /// reproduce the official `SET I ALL OR 1` `SHOW ELEMENTS` rows exactly.
    #[test]
    fn frozen_scalar_self_embeddings_pin_the_official_operations() {
        /// `(rotation, translation)` with the translation as reduced
        /// `(numerator, denominator)` pairs, exactly as the fixture stores it.
        type Expected = (Mat3I, [(i128, i128); 3]);
        let zero = [(0, 1), (0, 1), (0, 1)];
        let half_z = [(0, 1), (0, 1), (1, 2)];
        let half_xy = [(1, 2), (1, 2), (0, 1)];
        let half_yz = [(0, 1), (1, 2), (1, 2)];
        let half_xz = [(1, 2), (0, 1), (1, 2)];
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let c2x = [[1, 0, 0], [0, -1, 0], [0, 0, -1]];
        let c2y = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]];
        let c2z = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];
        let sgx = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let sgy = [[1, 0, 0], [0, -1, 0], [0, 0, 1]];
        let c4z_plus = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let c4z_minus = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];
        let sgz = [[1, 0, 0], [0, 1, 0], [0, 0, -1]];
        let s4z_minus = [[0, 1, 0], [-1, 0, 0], [0, 0, -1]];
        let s4z_plus = [[0, -1, 0], [1, 0, 0], [0, 0, -1]];
        let cases: [(u8, &str, u8, &[Expected]); 4] = [
            (
                19,
                "GM1",
                19,
                &[
                    (identity, zero),
                    (c2x, half_xy),
                    (c2y, half_yz),
                    (c2z, half_xz),
                ],
            ),
            (
                23,
                "GM1",
                23,
                &[(identity, zero), (c2x, zero), (c2y, zero), (c2z, zero)],
            ),
            (
                45,
                "GM1",
                45,
                &[(identity, zero), (c2z, zero), (sgx, half_z), (sgy, half_z)],
            ),
            (
                83,
                "GM1+",
                83,
                &[
                    (identity, zero),
                    (c2z, zero),
                    (c4z_plus, zero),
                    (c4z_minus, zero),
                    (inversion, zero),
                    (sgz, zero),
                    (s4z_minus, zero),
                    (s4z_plus, zero),
                ],
            ),
        ];
        for (parent, condensing, subgroup_sg, expected) in cases {
            let built = embedding(parent, condensing, "P1");
            assert_eq!(built.subgroup_sg(), subgroup_sg);
            assert_eq!(built.setting(), IDENTITY_SETTING);
            assert_eq!(built.child_shift(), &Vec3R::zero());
            assert_eq!(built.candidate_count(), 1);
            assert_eq!(built.representatives().len(), expected.len());
            let mut found: Vec<(Mat3I, [(i128, i128); 3])> = built
                .representatives()
                .iter()
                .map(|operation| {
                    (
                        operation.rotation(),
                        operation
                            .translation()
                            .as_array()
                            .map(|value| (value.numerator(), value.denominator())),
                    )
                })
                .collect();
            let mut expected: Vec<(Mat3I, [(i128, i128); 3])> = expected.to_vec();
            found.sort_by_key(|entry| (entry.0, entry.1));
            expected.sort_by_key(|entry| (entry.0, entry.1));
            assert_eq!(found, expected, "SG {parent} {condensing}");
        }
    }
}
