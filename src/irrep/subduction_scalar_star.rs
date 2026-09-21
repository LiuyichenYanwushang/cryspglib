//! Scalar parent full stars: the compound-aware source adapter of task 8c.
//!
//! [`OrdinaryStar`](super::OrdinaryStar) deliberately accepts only ordinary
//! scalar probes.  [`ScalarStar`] is the adapter the decomposition stage uses:
//! it accepts every canonical non-spinor record of the generated table and
//! expands the frozen compound metadata into **complex components**:
//!
//! * an ordinary row is one component, evaluated over the star of its stored
//!   `k`;
//! * a `DistinctComponentSum` row is two components, both at the record's
//!   stored `k`, each keeping its own CIR source number and selected-arm row;
//! * a `ConjugateRealification` row is its CIR seed at `k` and the conjugate of
//!   that seed at `-k`.  The two components are kept even when they turn out to
//!   be equivalent complex irreps; only the child target catalogue may merge
//!   them, and only after transporting both to one arm.
//!
//! Every component evaluates exactly the same induced trace as
//! [`OrdinaryStar`](super::OrdinaryStar): an arm that the operation moves
//! contributes zero, a fixed arm contributes `chi_seed(g_i^-1 h g_i)` with the
//! **stored** row and wave vector.  A conjugated component conjugates that
//! whole evaluated value — phase included — which is what makes its Bloch
//! phase the `-k` phase instead of the `+k` phase.  The full parent character
//! is the sum over the components, never a trace divided by an arm count, and
//! the assembled dimension is checked against `IrrepRecord::dim`.
//!
//! The folded geometry is shared with `OrdinaryStar` through
//! [`fold_arms`](super::fold_arms): it takes one annotated arm per (component,
//! arm) pair, so two components that fold onto one `q` keep their identities,
//! and it requires the frozen equality of the component dimensions (pinned for
//! every compound record by the full-table census).

use num_complex::Complex64;

use crate::irrep::query;
use crate::irrep::types::{
    CharacterRow, CompoundSelectedArmCharacter, IrrepRecord, IrrepSourceIdentity,
};

use super::super::{
    ExactSeitz, Lattice, SUBDUCTION_TOLERANCE, SubductionComponent, SubductionError,
    SubgroupEmbedding, Vec3R, exact_primitive_basis, inline_k_vector, reduce_operations,
    strict_sg_hall_ops,
};
use super::{
    FoldArm, FoldedStar, StarArm, StarError, arm_character, arm_wave_vector, collect_arms,
    fold_arms, induced_component_character,
};

// ── One complex component ────────────────────────────────────────────────────

/// One complex irreducible component of a scalar record, induced over its own
/// full star.
///
/// The row is always the **stored** row of its wave vector `base_k`; a
/// conjugated component stores `base_k` too and carries `conjugate = true`, so
/// the Bloch phase of the evaluation stays the stored `+base_k` phase and the
/// conjugation of the evaluated value supplies the `-base_k` phase exactly.
#[derive(Debug, Clone)]
pub struct ComponentStar {
    sg: u8,
    label: &'static str,
    component: SubductionComponent,
    irnumber: u32,
    base_k: Vec3R,
    effective_k: Vec3R,
    conjugate: bool,
    row: CharacterRow,
    arms: Vec<StarArm>,
    lattice: Lattice,
    reciprocal: Lattice,
    operations: Vec<ExactSeitz>,
}

impl ComponentStar {
    /// Build one component at `base_k`.
    ///
    /// `transporters` selects the transversal used to collect the arms of the
    /// component's **effective** wave vector: `None` scans the complete Hall
    /// operation list (the canonical constructor), `Some` validates an explicit
    /// one-arm-per-arm transversal with duplicates rejected.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        sg: u8,
        label: &'static str,
        component: SubductionComponent,
        irnumber: u32,
        base_k: Vec3R,
        row: CharacterRow,
        conjugate: bool,
        transporters: Option<&[ExactSeitz]>,
    ) -> Result<Self, StarError> {
        let lattice = Lattice::new(exact_primitive_basis(sg)?)?;
        let reciprocal = lattice.reciprocal()?;
        let hall = strict_sg_hall_ops(sg)?;
        let effective_k = if conjugate {
            base_k.checked_neg()?
        } else {
            base_k
        };
        let (transporters, duplicates_are_errors) = match transporters {
            Some(list) => (list, true),
            None => (hall.operations.as_slice(), false),
        };
        let arms = collect_arms(
            sg,
            &effective_k,
            &reciprocal,
            transporters,
            duplicates_are_errors,
        )?;
        let operations = reduce_operations(&hall.operations, &lattice)?;
        Ok(Self {
            sg,
            label,
            component,
            irnumber,
            base_k,
            effective_k,
            conjugate,
            row,
            arms,
            lattice,
            reciprocal,
            operations,
        })
    }

    /// Space group this component lives in (the parent for a parent record, the
    /// child for a child evaluator).
    pub const fn sg(&self) -> u8 {
        self.sg
    }

    /// Stable source label: the CIR label of a constituent, or the row label of
    /// an ordinary component.
    pub const fn label(&self) -> &'static str {
        self.label
    }

    /// Which complex constituent of the physical row this component is.
    pub const fn component(&self) -> SubductionComponent {
        self.component
    }

    /// Frozen CIR source number of this component.
    pub const fn irnumber(&self) -> u32 {
        self.irnumber
    }

    /// The stored wave vector the row belongs to.
    pub const fn base_k(&self) -> &Vec3R {
        &self.base_k
    }

    /// The effective wave vector whose star this component is induced over:
    /// `base_k`, or its exact negation for a conjugated component.
    pub const fn effective_k(&self) -> &Vec3R {
        &self.effective_k
    }

    /// Whether the evaluated trace is conjugated (the `-base_k` component).
    pub const fn conjugate(&self) -> bool {
        self.conjugate
    }

    /// The stored selected-arm character row, never conjugated in place.
    pub const fn row(&self) -> &CharacterRow {
        &self.row
    }

    /// One arm per class of the effective wave vector's star.
    pub fn arms(&self) -> &[StarArm] {
        &self.arms
    }

    /// Number of arms of this component's star.
    pub fn arm_count(&self) -> usize {
        self.arms.len()
    }

    /// Selected-arm dimension of this component.
    pub fn dimension(&self) -> usize {
        self.row.dimension()
    }

    /// Full-star dimension of this component: `arm_count × dimension`.
    pub fn full_dimension(&self) -> u64 {
        (self.arms.len() as u64) * (self.row.dimension() as u64)
    }

    /// The exact parent (child) lattice of this component's frame.
    pub const fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    /// The reciprocal lattice of [`Self::lattice`].
    pub const fn reciprocal(&self) -> &Lattice {
        &self.reciprocal
    }

    /// The complete operation universe reduced modulo the lattice.
    pub fn operations(&self) -> &[ExactSeitz] {
        &self.operations
    }

    /// The induced full-star character, after checking membership.
    pub fn character(&self, operation: &ExactSeitz) -> Result<Complex64, StarError> {
        let reduced = operation.reduce(&self.lattice)?;
        if !self.operations.contains(&reduced) {
            return Err(StarError::OperationNotInParentGroup { sg: self.sg });
        }
        self.trace(operation)
    }

    /// The induced full-star character without the membership check.
    pub(super) fn trace(&self, operation: &ExactSeitz) -> Result<Complex64, StarError> {
        induced_component_character(
            &self.arms,
            &self.row,
            self.label,
            &self.lattice,
            &self.reciprocal,
            &self.base_k,
            self.conjugate,
            operation,
        )
    }

    /// One arm's induced contribution, with the component's conjugation rule.
    pub(super) fn arm_trace(
        &self,
        arm: usize,
        operation: &ExactSeitz,
    ) -> Result<Complex64, StarError> {
        let star_arm = self.arms.get(arm).ok_or(StarError::UnknownArmIndex {
            index: arm,
            arms: self.arms.len(),
        })?;
        arm_character(
            star_arm,
            &self.row,
            self.label,
            &self.lattice,
            &self.reciprocal,
            &self.base_k,
            self.conjugate,
            arm,
            operation,
        )
    }
}

// ── The scalar full star ─────────────────────────────────────────────────────

/// One arm of the assembled scalar star, tagged with its component identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScalarArm {
    wave_vector: Vec3R,
    component: usize,
    arm: usize,
}

impl ScalarArm {
    /// The exact arm wave vector (of the component's effective wave vector).
    pub const fn wave_vector(&self) -> &Vec3R {
        &self.wave_vector
    }

    /// Index into [`ScalarStar::components`].
    pub const fn component(&self) -> usize {
        self.component
    }

    /// Index into that component's own arm list.
    pub const fn arm(&self) -> usize {
        self.arm
    }
}

/// The full star of one canonical non-spinor scalar record: ordinary,
/// `DistinctComponentSum` or `ConjugateRealification`.
#[derive(Debug, Clone)]
pub struct ScalarStar {
    parent_sg: u8,
    probe: &'static IrrepRecord,
    source_identity: IrrepSourceIdentity,
    seed_k: Vec3R,
    components: Vec<ComponentStar>,
    arms: Vec<ScalarArm>,
    parent_lattice: Lattice,
    parent_reciprocal: Lattice,
    parent_operations: Vec<ExactSeitz>,
    dimension: u32,
}

impl ScalarStar {
    /// Build the full star of a canonical non-spinor probe.
    ///
    /// The probe must be the generated table's own record (identity checked
    /// against [`query::irreps_of`]); spinors and rows without a scalar typed
    /// character space fail with [`SubductionError::UnsupportedCharacterSpace`]
    /// before any arm is built.  The assembled dimension must equal
    /// `IrrepRecord::dim` *and* the identity character, otherwise construction
    /// fails closed.
    pub fn new(probe: &'static IrrepRecord) -> Result<Self, StarError> {
        Self::validate_probe(probe)?;
        let seed_k = inline_k_vector(probe)?;
        let components = Self::expand_components(probe, &seed_k, None)?;
        Self::assemble(probe, seed_k, components)
    }

    /// Build the star from one explicit transversal per complex component.
    ///
    /// This is the reusable constructor behind [`Self::new`] and the arm-order
    /// invariance witness.  Every transporter must be a parent group element
    /// modulo `L_G`, duplicates within a component are rejected, and the
    /// assembled arms must still cover every parent operation's image, so a
    /// missing arm fails as [`StarError::IncompleteStar`] rather than folding a
    /// smaller star.
    pub fn from_transporters(
        probe: &'static IrrepRecord,
        transporters: &[Vec<ExactSeitz>],
    ) -> Result<Self, StarError> {
        Self::validate_probe(probe)?;
        let seed_k = inline_k_vector(probe)?;
        let expected = Self::component_count(probe)?;
        if transporters.len() != expected {
            return Err(StarError::ComponentTransporterCount {
                expected,
                found: transporters.len(),
            });
        }
        let parent_lattice = Lattice::new(exact_primitive_basis(probe.sg)?)?;
        let parent_operations =
            reduce_operations(&strict_sg_hall_ops(probe.sg)?.operations, &parent_lattice)?;
        for (component, list) in transporters.iter().enumerate() {
            for (index, transporter) in list.iter().enumerate() {
                let reduced = transporter.reduce(&parent_lattice)?;
                if !parent_operations.contains(&reduced) {
                    return Err(StarError::ComponentTransporterNotInParentGroup {
                        sg: probe.sg,
                        component,
                        index,
                    });
                }
            }
        }
        let components = Self::expand_components(probe, &seed_k, Some(transporters))?;
        Self::assemble(probe, seed_k, components)
    }

    /// Probe membership and representation space.
    fn validate_probe(probe: &'static IrrepRecord) -> Result<(), StarError> {
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
        if probe.spinor {
            return Err(SubductionError::UnsupportedCharacterSpace {
                sg: probe.sg,
                ml: probe.ml.to_string(),
            }
            .into());
        }
        Ok(())
    }

    /// Number of complex components of a scalar record.
    fn component_count(probe: &'static IrrepRecord) -> Result<usize, StarError> {
        if probe.compound_metadata().is_some() {
            Ok(2)
        } else if probe.ordinary_scalar_selected_arm_block_trace().is_ok() {
            Ok(1)
        } else {
            Err(SubductionError::UnsupportedCharacterSpace {
                sg: probe.sg,
                ml: probe.ml.to_string(),
            }
            .into())
        }
    }

    /// Expand the frozen metadata into complex components at `seed_k`.
    ///
    /// A `ConjugateRealification` keeps **both** components: the seed at `k` and
    /// its conjugate at `-k`, each with the same CIR source number and row, and
    /// the second one flagged `conjugate`.
    fn expand_components(
        probe: &'static IrrepRecord,
        seed_k: &Vec3R,
        transporters: Option<&[Vec<ExactSeitz>]>,
    ) -> Result<Vec<ComponentStar>, StarError> {
        let explicit = |index: usize| -> Option<&[ExactSeitz]> {
            transporters.map(|lists| lists[index].as_slice())
        };
        let unsupported = || SubductionError::UnsupportedCharacterSpace {
            sg: probe.sg,
            ml: probe.ml.to_string(),
        };
        match probe.compound_metadata() {
            None => {
                let row = probe
                    .ordinary_scalar_selected_arm_block_trace()
                    .map_err(|_| unsupported())?;
                let irnumber = match probe.source_identity() {
                    IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => cir_irnumber,
                    _ => return Err(unsupported().into()),
                };
                Ok(vec![ComponentStar::new(
                    probe.sg,
                    probe.ml,
                    SubductionComponent::Ordinary,
                    irnumber,
                    *seed_k,
                    row,
                    false,
                    explicit(0),
                )?])
            }
            Some(_) => {
                let view = probe
                    .compound_selected_arm_view()
                    .map_err(|_| unsupported())?;
                match view {
                    CompoundSelectedArmCharacter::DistinctComponentSum {
                        first, second, ..
                    } => Ok(vec![
                        ComponentStar::new(
                            probe.sg,
                            first.label,
                            SubductionComponent::Constituent {
                                index: 0,
                                irnumber: first.irnumber,
                            },
                            first.irnumber,
                            *seed_k,
                            first.row,
                            false,
                            explicit(0),
                        )?,
                        ComponentStar::new(
                            probe.sg,
                            second.label,
                            SubductionComponent::Constituent {
                                index: 1,
                                irnumber: second.irnumber,
                            },
                            second.irnumber,
                            *seed_k,
                            second.row,
                            false,
                            explicit(1),
                        )?,
                    ]),
                    CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => Ok(vec![
                        ComponentStar::new(
                            probe.sg,
                            seed.label,
                            SubductionComponent::RealificationSeed {
                                irnumber: seed.irnumber,
                            },
                            seed.irnumber,
                            *seed_k,
                            seed.row.clone(),
                            false,
                            explicit(0),
                        )?,
                        ComponentStar::new(
                            probe.sg,
                            seed.label,
                            SubductionComponent::RealificationConjugate {
                                irnumber: seed.irnumber,
                            },
                            seed.irnumber,
                            *seed_k,
                            seed.row,
                            true,
                            explicit(1),
                        )?,
                    ]),
                }
            }
        }
    }

    /// Check the assembles invariants and build the star.
    fn assemble(
        probe: &'static IrrepRecord,
        seed_k: Vec3R,
        components: Vec<ComponentStar>,
    ) -> Result<Self, StarError> {
        let parent_lattice = Lattice::new(exact_primitive_basis(probe.sg)?)?;
        let parent_reciprocal = parent_lattice.reciprocal()?;
        let hall = strict_sg_hall_ops(probe.sg)?;
        let parent_operations = reduce_operations(&hall.operations, &parent_lattice)?;

        // The folded geometry needs one arm dimension; the frozen table has
        // equal constituent dimensions for every compound record.
        let selected = components
            .first()
            .map(ComponentStar::dimension)
            .ok_or(StarError::FoldLostArms { folded: 0, arms: 0 })?;
        if let Some(other) = components
            .iter()
            .find(|component| component.dimension() != selected)
        {
            return Err(StarError::NonUniformComponentDimension {
                first: selected,
                other: other.dimension(),
            });
        }

        // Completeness per component: its arms are exactly the star of its
        // effective wave vector, so every parent operation lands on one arm.
        for component in &components {
            for operation in &hall.operations {
                let image = arm_wave_vector(component.effective_k(), operation)?;
                let mut covered = false;
                for arm in component.arms() {
                    if parent_reciprocal.same_mod(arm.wave_vector(), &image)? {
                        covered = true;
                        break;
                    }
                }
                if !covered {
                    return Err(StarError::IncompleteStar {
                        sg: probe.sg,
                        ml: probe.ml,
                        represented: component.arm_count(),
                        operations: hall.operations.len(),
                        missing_k: Box::new([image.get(0), image.get(1), image.get(2)]),
                    });
                }
            }
        }

        let mut arms = Vec::new();
        let mut dimension = 0u32;
        for (index, component) in components.iter().enumerate() {
            let full = u32::try_from(component.full_dimension()).map_err(|_| {
                SubductionError::RationalOverflow {
                    operation: "scalar full-star dimension",
                }
            })?;
            dimension = dimension
                .checked_add(full)
                .ok_or(SubductionError::RationalOverflow {
                    operation: "scalar full-star dimension",
                })?;
            for (arm, star_arm) in component.arms().iter().enumerate() {
                arms.push(ScalarArm {
                    wave_vector: *star_arm.wave_vector(),
                    component: index,
                    arm,
                });
            }
        }

        let star = Self {
            parent_sg: probe.sg,
            probe,
            source_identity: probe.source_identity(),
            seed_k,
            components,
            arms,
            parent_lattice,
            parent_reciprocal,
            parent_operations,
            dimension,
        };
        let identity = star.character(&ExactSeitz::identity())?;
        let expected = Complex64::new(f64::from(dimension), 0.0);
        if u32::from(probe.dim) != dimension || (identity - expected).norm() > SUBDUCTION_TOLERANCE
        {
            return Err(StarError::ScalarDimensionMismatch {
                ml: probe.ml,
                record: u32::from(probe.dim),
                found: dimension,
                identity,
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

    /// Frozen source identity of the probe.
    pub const fn source_identity(&self) -> IrrepSourceIdentity {
        self.source_identity
    }

    /// The probe's stored selected-arm wave vector.
    pub const fn seed_k(&self) -> &Vec3R {
        &self.seed_k
    }

    /// Every complex component, in the record's constituent order.
    pub fn components(&self) -> &[ComponentStar] {
        &self.components
    }

    /// One component by index.
    pub fn component(&self, index: usize) -> Option<&ComponentStar> {
        self.components.get(index)
    }

    /// Every arm of every component, identity arm first per component.
    pub fn arms(&self) -> &[ScalarArm] {
        &self.arms
    }

    /// One arm by index.
    pub fn arm(&self, index: usize) -> Option<&ScalarArm> {
        self.arms.get(index)
    }

    /// Total number of arms over all components.
    pub fn arm_count(&self) -> usize {
        self.arms.len()
    }

    /// The parent translation lattice `L_G`.
    pub const fn parent_lattice(&self) -> &Lattice {
        &self.parent_lattice
    }

    /// The parent reciprocal lattice.
    pub const fn parent_reciprocal(&self) -> &Lattice {
        &self.parent_reciprocal
    }

    /// Full-star dimension: `Σ_component arm_count × selected_dimension`.
    pub const fn dimension(&self) -> u32 {
        self.dimension
    }

    /// The induced full-star character of the whole physical row, checked to be
    /// a parent operation modulo `L_G`.
    pub fn character(&self, operation: &ExactSeitz) -> Result<Complex64, StarError> {
        let reduced = operation.reduce(&self.parent_lattice)?;
        if !self.parent_operations.contains(&reduced) {
            return Err(StarError::OperationNotInParentGroup { sg: self.parent_sg });
        }
        let mut total = Complex64::new(0.0, 0.0);
        for component in &self.components {
            total += component.trace(operation)?;
        }
        Ok(total)
    }

    /// The q-block character on one operation: the sum of the induced
    /// contributions of exactly the folded arms at that `q`.
    ///
    /// Duplicate arms of several components are preserved; the trace is never
    /// divided by an arm count.
    pub fn q_block_character(
        &self,
        arm_indices: &[usize],
        operation: &ExactSeitz,
    ) -> Result<Complex64, StarError> {
        let mut total = Complex64::new(0.0, 0.0);
        for index in arm_indices {
            let arm = self.arm(*index).ok_or(StarError::UnknownArmIndex {
                index: *index,
                arms: self.arms.len(),
            })?;
            let component = self
                .component(arm.component())
                .ok_or(StarError::UnknownArmIndex {
                    index: *index,
                    arms: self.arms.len(),
                })?;
            total += component.arm_trace(arm.arm(), operation)?;
        }
        Ok(total)
    }

    /// The q-block dimension: the sum of the selected-arm dimensions of exactly
    /// the folded arms at that `q`.
    pub fn q_block_dimension(&self, arm_indices: &[usize]) -> Result<u32, StarError> {
        let mut total = 0u32;
        for index in arm_indices {
            let arm = self.arm(*index).ok_or(StarError::UnknownArmIndex {
                index: *index,
                arms: self.arms.len(),
            })?;
            let component = self
                .component(arm.component())
                .ok_or(StarError::UnknownArmIndex {
                    index: *index,
                    arms: self.arms.len(),
                })?;
            total = total
                .checked_add(u32::try_from(component.dimension()).map_err(|_| {
                    SubductionError::RationalOverflow {
                        operation: "q-block dimension",
                    }
                })?)
                .ok_or(SubductionError::RationalOverflow {
                    operation: "q-block dimension",
                })?;
        }
        Ok(total)
    }

    /// Fold every component arm into the child frame and partition the folded
    /// points into child stars.
    ///
    /// The geometry itself is [`fold_arms`](super::fold_arms); the annotations
    /// keep the component identity of every arm, so duplicate arms that merge
    /// on one `q` survive the merge.
    pub fn folded_stars(
        &self,
        embedding: &SubgroupEmbedding,
    ) -> Result<Vec<FoldedStar>, StarError> {
        let arms: Vec<FoldArm> = self
            .arms
            .iter()
            .map(|arm| FoldArm {
                wave_vector: *arm.wave_vector(),
                dimension: self.components[arm.component()].dimension(),
            })
            .collect();
        fold_arms(embedding, self.parent_sg, &arms)
    }

    /// The same star with its complex components listed in the opposite order.
    /// Test-only witness that the decomposition is ordered by the child table
    /// and not by the parent constituent order; the arm annotations are rebuilt
    /// so the geometry stays coherent.
    #[cfg(test)]
    pub(super) fn with_components_reversed(&self) -> Self {
        let mut swapped = self.clone();
        swapped.components.reverse();
        swapped.arms = swapped
            .components
            .iter()
            .enumerate()
            .flat_map(|(index, component)| {
                component
                    .arms()
                    .iter()
                    .enumerate()
                    .map(move |(arm, star_arm)| ScalarArm {
                        wave_vector: *star_arm.wave_vector(),
                        component: index,
                        arm,
                    })
            })
            .collect();
        swapped
    }

    /// The same star with its last complex component removed.
    ///
    /// Test-only witness that a lost component fails the block tiling instead
    /// of silently returning a smaller decomposition; the declared dimension is
    /// deliberately left stale.
    #[cfg(test)]
    pub(super) fn without_last_component(&self) -> Self {
        let mut reduced = self.clone();
        reduced.components.pop();
        let kept = reduced.components.len();
        reduced.arms.retain(|arm| arm.component < kept);
        reduced
    }

    /// The same star with its last arm removed.
    ///
    /// Test-only witness that a lost arm fails the block tiling instead of
    /// silently returning a smaller decomposition.
    #[cfg(test)]
    pub(super) fn without_last_arm(&self) -> Self {
        let mut reduced = self.clone();
        reduced.arms.pop();
        reduced
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::subduction::Rat;
    use crate::irrep::types::CompoundCharacterSemantics;

    /// Construct the scalar star of every canonical non-spinor record of the
    /// frozen table and pin the full census.
    ///
    /// The counts are the structural ones of the frozen metadata: 4105 ordinary
    /// rows and 672 compound rows, of which 519 are `DistinctComponentSum` and
    /// 153 `ConjugateRealification`; the realifications split into 112 whose
    /// `k` and `-k` stars are disjoint and 41 whose `-k` lies in the same star.
    /// Every record's assembled dimension equals `IrrepRecord::dim` and its
    /// identity character, every compound pair has equal frozen constituent
    /// dimensions with `(d0 + d1) x arms == dim`, and every distinct pair keeps
    /// two components on exactly the same star.
    #[test]
    fn the_scalar_adapter_constructs_every_frozen_scalar_record() {
        let mut ordinary = 0usize;
        let mut distinct = 0usize;
        let mut realification = 0usize;
        let mut disjoint = 0usize;
        let mut same = 0usize;
        for sg in 1..=230u8 {
            for record in query::irreps_of(sg) {
                if record.spinor {
                    continue;
                }
                let star = ScalarStar::new(record).unwrap_or_else(|error| {
                    panic!(
                        "SG {sg} {}: scalar star construction failed: {error}",
                        record.ml
                    )
                });
                assert_eq!(
                    star.dimension(),
                    u32::from(record.dim),
                    "SG {sg} {}",
                    record.ml
                );
                let identity = star
                    .character(&ExactSeitz::identity())
                    .expect("identity is a parent operation");
                assert!(
                    (identity - Complex64::new(f64::from(star.dimension()), 0.0)).norm()
                        <= SUBDUCTION_TOLERANCE,
                    "SG {sg} {}: {identity}",
                    record.ml
                );
                match record.compound_metadata() {
                    None => {
                        ordinary += 1;
                        assert_eq!(star.components().len(), 1);
                        assert_eq!(
                            star.components()[0].component(),
                            SubductionComponent::Ordinary
                        );
                        assert!(matches!(
                            star.source_identity(),
                            IrrepSourceIdentity::OrdinaryScalar { .. }
                        ));
                        assert_eq!(star.components()[0].label(), record.ml);
                    }
                    Some(metadata) => {
                        assert_eq!(star.components().len(), 2, "SG {sg} {}", record.ml);
                        assert_eq!(
                            metadata.cir_dimensions[0], metadata.cir_dimensions[1],
                            "SG {sg} {}: unequal frozen constituent dimensions",
                            record.ml
                        );
                        let selected = usize::from(metadata.cir_dimensions[0]);
                        assert_eq!(
                            (2 * selected * star.components()[0].arm_count()) as u32,
                            u32::from(record.dim),
                            "SG {sg} {}",
                            record.ml
                        );
                        assert!(matches!(
                            star.source_identity(),
                            IrrepSourceIdentity::Compound { .. }
                        ));
                        match metadata.semantics {
                            CompoundCharacterSemantics::DistinctComponentSum => {
                                distinct += 1;
                                let first = &star.components()[0];
                                let second = &star.components()[1];
                                assert_eq!(
                                    first.component(),
                                    SubductionComponent::Constituent {
                                        index: 0,
                                        irnumber: metadata.cir_irnumbers[0],
                                    },
                                    "SG {sg} {}",
                                    record.ml
                                );
                                assert_eq!(
                                    second.component(),
                                    SubductionComponent::Constituent {
                                        index: 1,
                                        irnumber: metadata.cir_irnumbers[1],
                                    },
                                    "SG {sg} {}",
                                    record.ml
                                );
                                assert_eq!(first.label(), metadata.cir_labels[0]);
                                assert_eq!(second.label(), metadata.cir_labels[1]);
                                // Both constituents are induced over exactly
                                // the same star.
                                let first_arms: Vec<Vec3R> =
                                    first.arms().iter().map(|arm| *arm.wave_vector()).collect();
                                let second_arms: Vec<Vec3R> =
                                    second.arms().iter().map(|arm| *arm.wave_vector()).collect();
                                assert_eq!(first_arms, second_arms, "SG {sg} {}", record.ml);
                            }
                            CompoundCharacterSemantics::ConjugateRealification => {
                                realification += 1;
                                let seed = &star.components()[0];
                                let conjugate = &star.components()[1];
                                assert_eq!(
                                    seed.component(),
                                    SubductionComponent::RealificationSeed {
                                        irnumber: metadata.cir_irnumbers[0],
                                    }
                                );
                                assert_eq!(
                                    conjugate.component(),
                                    SubductionComponent::RealificationConjugate {
                                        irnumber: metadata.cir_irnumbers[0],
                                    }
                                );
                                assert!(!seed.conjugate());
                                assert!(conjugate.conjugate());
                                assert_eq!(seed.label(), conjugate.label());
                                assert_eq!(seed.irnumber(), conjugate.irnumber());
                                assert_eq!(seed.arm_count(), conjugate.arm_count());
                                assert_eq!(
                                    seed.dimension(),
                                    usize::from(metadata.cir_dimensions[0])
                                );
                                let negated = seed.base_k().checked_neg().expect("negation");
                                assert_eq!(conjugate.effective_k(), &negated);
                                let negated_in_star = seed.arms().iter().any(|arm| {
                                    star.parent_reciprocal()
                                        .same_mod(arm.wave_vector(), &negated)
                                        .expect("same_mod")
                                });
                                if negated_in_star {
                                    same += 1;
                                } else {
                                    disjoint += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(ordinary, 4105, "ordinary scalar records");
        assert_eq!(distinct, 519, "DistinctComponentSum records");
        assert_eq!(realification, 153, "ConjugateRealification records");
        assert_eq!(disjoint, 112, "realifications with disjoint k/-k stars");
        assert_eq!(same, 41, "realifications with k/-k in one star");
    }

    /// A conjugated component evaluates the conjugate of the whole induced
    /// trace, phase included, at every parent operation.
    #[test]
    fn the_conjugate_component_is_the_exact_complex_conjugate() {
        for (sg, ml) in [(19u8, "R1R1"), (23, "W1W1"), (45, "W1W1"), (45, "W2W2")] {
            let star = ScalarStar::new(
                query::irreps_of(sg)
                    .iter()
                    .find(|record| record.ml == ml)
                    .expect("pinned record"),
            )
            .expect("star");
            assert_eq!(star.components().len(), 2);
            for operation in &strict_sg_hall_ops(sg).expect("hall").operations {
                let seed = star.components()[0].trace(operation).expect("seed");
                let conjugate = star.components()[1].trace(operation).expect("conjugate");
                assert!(
                    (conjugate - seed.conj()).norm() <= SUBDUCTION_TOLERANCE,
                    "SG {sg} {ml}: {conjugate} != conj({seed})"
                );
            }
        }
    }

    /// Explicit transversals are validated: one per component, identity first,
    /// parent operations only, and no missing arm.
    #[test]
    fn explicit_component_transversals_are_validated() {
        let record = query::irreps_of(83)
            .iter()
            .find(|record| record.ml == "GM3+GM4+")
            .expect("pinned record");
        let canonical = ScalarStar::new(record).expect("star");
        let transporters: Vec<Vec<ExactSeitz>> = canonical
            .components()
            .iter()
            .map(|component| {
                component
                    .arms()
                    .iter()
                    .map(|arm| *arm.transporter())
                    .collect()
            })
            .collect();
        let rebuilt = ScalarStar::from_transporters(record, &transporters).expect("rebuilt");
        assert_eq!(rebuilt.arms(), canonical.arms());
        assert_eq!(rebuilt.dimension(), canonical.dimension());

        // A wrong component count is a typed error.
        assert!(matches!(
            ScalarStar::from_transporters(record, &transporters[..1]),
            Err(StarError::ComponentTransporterCount {
                expected: 2,
                found: 1
            })
        ));

        // A forged second component arm is not a parent operation.  A lattice
        // translation would be, so the forged one is a genuine non-lattice
        // fractional translation this group does not contain.
        let forged = ExactSeitz::new(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            Vec3R::new([Rat::new(1, 3).expect("exact third"), Rat::ZERO, Rat::ZERO]),
        );
        let mut forged_transporters = transporters.clone();
        forged_transporters[1] = vec![ExactSeitz::identity(), forged];
        assert!(matches!(
            ScalarStar::from_transporters(record, &forged_transporters),
            Err(StarError::ComponentTransporterNotInParentGroup {
                sg: 83,
                component: 1,
                index: 1
            })
        ));
    }

    /// `OrdinaryStar` keeps rejecting compound rows while the scalar adapter
    /// accepts them: the two adapters have different, explicit spaces.
    #[test]
    fn ordinary_star_still_rejects_compound_and_scalar_accepts_it() {
        let compound = query::irreps_of(83)
            .iter()
            .find(|record| record.ml == "GM3+GM4+")
            .expect("pinned record");
        assert!(matches!(
            super::super::OrdinaryStar::new(compound),
            Err(StarError::Subduction(
                SubductionError::UnsupportedCharacterSpace { sg: 83, .. }
            ))
        ));
        let star = ScalarStar::new(compound).expect("scalar star");
        assert_eq!(star.dimension(), 2);
    }
}
