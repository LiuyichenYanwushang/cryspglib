//! No-SOC collinear spin symmetry, independently of magnetic space groups.
//!
//! The base detection/validation API returns a **unitary finite quotient**, not
//! a complete spin space group or an SU(2) representation. `Swap` is a fixed spin
//! pi rotation perpendicular to the magnetic axis, not time reversal. Its
//! momentum action is `k -> W^{-T} k`, without complex conjugation.
//! [`ValidatedCollinearSpinGroup::with_spin_only_time_reversal`] includes the
//! ordinary collinear symmetry `A=[C2_perp T||E]`: A is antiunitary but PRESERVES
//! channels, conjugates orbital coefficients and sends `k -> -k`. Its coset gA
//! has the same channel action as g and momentum `-W^{-T}k`. Spatial rotations
//! never rotate the independent spin axis.
//!
//! The continuous rotations about that axis need not be sampled for an already
//! block-diagonal `diag(H_up, H_down)`. SOC, external magnetic fields and orbital
//! flux are outside this geometry-only certificate; callers must exclude them
//! or separately prove Hamiltonian covariance. Wannier gauges, orbital sewing
//! matrices and Hamiltonian projection belong to the consuming code.
//!
//! ```
//! use cryspglib::{Crystal, CollinearMoments, CollinearSpinTolerances,
//!     SpinChannelAction, detect_collinear_spin_group};
//! let crystal = Crystal::new(
//!     [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]],
//!     vec![[0.; 3], [0.5, 0., 0.]], vec![26, 26],
//! )?;
//! let group = detect_collinear_spin_group(
//!     &crystal, CollinearMoments::Signed(&[2., -2.]),
//!     CollinearSpinTolerances::default(),
//! )?;
//! assert!(group.operations().iter().any(|g|
//!     g.operation.channel_action == SpinChannelAction::Swap));
//! // Ordinary no-SOC collinear model, with orbital time-reversal covariance:
//! let full = group.with_spin_only_time_reversal();
//! let a = full.operation(full.spin_only_time_reversal_index()).unwrap();
//! assert!(a.kind.is_antiunitary());
//! assert_eq!(a.operation.channel_action, SpinChannelAction::Preserve);
//! assert_eq!(a.momentum_sign(), -1);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::cell::{Cell, TensorRank};
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_get_determinant_d3, mat_get_metric, mat_inverse_matrix_d3,
};
use crate::symmetry::get_lattice_symmetry;
use crate::{Crystal, OperationKind};

const IDENTITY: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
// Above this bound f64 cannot reliably distinguish neighbouring cell shifts.
const MAX_CELL_INTEGER: f64 = (1_u64 << 52) as f64;
const MAX_IMAGES: u64 = 4096;

/// Permutation of the two spin blocks; unrelated to antiunitarity/time reversal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpinChannelAction {
    Preserve,
    Swap,
}

impl SpinChannelAction {
    /// The sign multiplying the signed moment, or the `H_z` component.
    pub const fn sign(self) -> i8 {
        match self {
            Self::Preserve => 1,
            Self::Swap => -1,
        }
    }

    /// Channel action of `self` after `right` (multiplication of the two signs).
    pub const fn compose(self, right: Self) -> Self {
        if self.sign() == right.sign() {
            Self::Preserve
        } else {
            Self::Swap
        }
    }
}

/// A unitary representative `[spin || {W|t}]` of the collinear finite quotient.
///
/// Input translations are retained, including their integer parts. The spin
/// representative is identity for `Preserve`, a fixed perpendicular pi rotation
/// for `Swap`. Its spinor central phase is deliberately not represented here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CollinearSpinOperation {
    pub rotation: Mat3I,
    pub translation: Vec3,
    pub channel_action: SpinChannelAction,
}

impl CollinearSpinOperation {
    pub const fn identity() -> Self {
        Self {
            rotation: IDENTITY,
            translation: [0.; 3],
            channel_action: SpinChannelAction::Preserve,
        }
    }
}

/// Magnetic input, in the original atom order. No MSG axial-vector rule is used.
#[derive(Debug, Clone, Copy)]
pub enum CollinearMoments<'a> {
    /// Signed magnitudes relative to a caller-chosen common spin axis.
    Signed(&'a [f64]),
    /// Cartesian vectors. A common axis is inferred from a largest moment;
    /// every perpendicular residual is checked. Signed magnitudes retain the
    /// actual vector norms; only their signs come from the axis projection.
    Vectors(&'a [Vec3]),
}

/// Independent absolute tolerances. No sign-only magnetic comparison is made.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CollinearSpinTolerances {
    /// Cartesian distance, also used for lattice metric and Seitz comparisons.
    /// Must be finite and strictly positive.
    pub position: f64,
    /// Absolute signed-magnitude residual, in the input moment units (>= 0).
    pub moment: f64,
    /// Cartesian perpendicular moment norm, in the input moment units (>= 0).
    pub collinearity: f64,
}

impl Default for CollinearSpinTolerances {
    fn default() -> Self {
        Self {
            position: 1e-5,
            moment: 1e-5,
            collinearity: 1e-5,
        }
    }
}

/// Mapping of one source atom (its index in the enclosing mapping array).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CollinearAtomMapping {
    pub target: usize,
    /// Integer shift in the original, not standardized, lattice basis.
    pub cell_shift: [i64; 3],
    /// Fractional residual: `W r_source + t = r_target + cell_shift + residual`.
    /// Its Cartesian norm is at most `tolerances.position`.
    pub position_residual: Vec3,
    /// `m_target - channel_action.sign() * m_source`, not an absolute value.
    pub moment_residual: f64,
}

/// An operation and its bijective atom mapping in the original input frame.
#[derive(Debug, Clone, PartialEq)]
pub struct MappedCollinearSpinOperation {
    pub operation: CollinearSpinOperation,
    pub atom_mapping: Vec<CollinearAtomMapping>,
}

/// Witness for `g h = {I|cell_shift + translation_residual} p` (first h, then g).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CollinearSpinProduct {
    pub operation_index: usize,
    pub cell_shift: [i64; 3],
    /// Fractional Seitz mismatch, explicitly retained rather than rounded away.
    pub translation_residual: Vec3,
}

/// A numerical group/geometry certificate, not an orbital/Hamiltonian certificate.
///
/// Construction checks the entire supplied set without removing or adding
/// entries. Rotation, permutation, channel multiplication and integer-shift
/// composition are exact; geometric and moment residuals use the stated bounds.
/// Only immutable access is provided so the certificate cannot be invalidated.
#[derive(Debug, Clone)]
pub struct ValidatedCollinearSpinGroup {
    operations: Vec<MappedCollinearSpinOperation>,
    multiplication_table: Vec<Vec<CollinearSpinProduct>>,
    identity_index: usize,
    inverse_indices: Vec<usize>,
    tolerances: CollinearSpinTolerances,
    signed_moments: Vec<f64>,
    spin_axis: Option<Vec3>,
    collinearity_residuals: Vec<f64>,
}

impl ValidatedCollinearSpinGroup {
    pub fn operations(&self) -> &[MappedCollinearSpinOperation] {
        &self.operations
    }
    /// Row g, column h is the product g after h, in the original operation order.
    pub fn multiplication_table(&self) -> &[Vec<CollinearSpinProduct>] {
        &self.multiplication_table
    }
    pub fn identity_index(&self) -> usize {
        self.identity_index
    }
    pub fn inverse_indices(&self) -> &[usize] {
        &self.inverse_indices
    }
    pub fn tolerances(&self) -> CollinearSpinTolerances {
        self.tolerances
    }
    pub fn signed_moments(&self) -> &[f64] {
        &self.signed_moments
    }
    /// Cartesian unit axis for vector input; `None` for caller-signed input.
    pub fn spin_axis(&self) -> Option<Vec3> {
        self.spin_axis
    }
    /// Per-atom perpendicular norms (zero for signed input).
    pub fn collinearity_residuals(&self) -> &[f64] {
        &self.collinearity_residuals
    }

    /// Include `A = [C2_perp T || E]`, the channel-preserving antiunitary spin
    /// symmetry of an ordinary no-SOC collinear Hamiltonian. This explicit opt-in
    /// asserts the absence of orbital time-reversal breaking; geometry cannot
    /// certify that assumption or the consumer's Wannier sewing matrices.
    ///
    /// The unitary subgroup is retained unchanged. Its central extension by A
    /// reuses all mappings and product witnesses without copying a doubled table.
    pub fn with_spin_only_time_reversal(self) -> CollinearSpinTimeReversalGroup {
        CollinearSpinTimeReversalGroup { unitary: self }
    }
}

/// A borrowed semilinear action `g A^kappa`. The underlying unitary factor g
/// supplies the spatial action and the TOTAL channel action: A preserves channels.
/// `kind` is independent of Preserve/Swap and is not an MSG axial-vector rule.
#[derive(Debug, Clone, Copy)]
pub struct CollinearSpinOperationRef<'a> {
    pub operation: &'a CollinearSpinOperation,
    pub atom_mapping: &'a [CollinearAtomMapping],
    pub kind: OperationKind,
}

impl CollinearSpinOperationRef<'_> {
    /// Momentum is `momentum_sign() * W^{-T} k`. Only antiunitarity adds a minus.
    pub const fn momentum_sign(self) -> i8 {
        if self.kind.is_antiunitary() { -1 } else { 1 }
    }
}

/// Conditional central extension `G x {1,A}` with `A=[C2_perp T||E]`.
///
/// A is antiunitary, preserves up/dn, and squares to +I for spin 1/2. It implies
/// same-channel `E_sigma(k)=E_sigma(-k)`, NOT up/down or Kramers degeneracy.
/// Other antiunitary elements can have nontrivial spinor/Bloch squares: this
/// quotient, like G, does not encode SU(2) central phases or orbital sewing.
///
/// Indices `0..n` retain G's original order; `n..2*n` represent `g A` in that same
/// order. Full-group products are derived exactly from G and antiunitary XOR;
/// no geometric or numerical checks are weakened or repeated.
#[derive(Debug, Clone)]
pub struct CollinearSpinTimeReversalGroup {
    unitary: ValidatedCollinearSpinGroup,
}

impl CollinearSpinTimeReversalGroup {
    pub fn unitary_subgroup(&self) -> &ValidatedCollinearSpinGroup {
        &self.unitary
    }
    pub fn len(&self) -> usize {
        2 * self.unitary.operations.len()
    }
    pub fn is_empty(&self) -> bool {
        self.unitary.operations.is_empty()
    }
    pub fn identity_index(&self) -> usize {
        self.unitary.identity_index
    }
    /// A's quotient representative. Original identity Seitz/mapping representatives
    /// are retained, so this need not have literal zero translation for explicit G.
    pub fn spin_only_time_reversal_index(&self) -> usize {
        self.unitary.operations.len() + self.unitary.identity_index
    }
    pub fn operation(&self, index: usize) -> Option<CollinearSpinOperationRef<'_>> {
        if index >= self.len() {
            return None;
        }
        let n = self.unitary.operations.len();
        let mapped = &self.unitary.operations[index % n];
        Some(CollinearSpinOperationRef {
            operation: &mapped.operation,
            atom_mapping: &mapped.atom_mapping,
            kind: if index < n {
                OperationKind::Unitary
            } else {
                OperationKind::Antiunitary
            },
        })
    }
    pub fn operations(&self) -> impl ExactSizeIterator<Item = CollinearSpinOperationRef<'_>> {
        (0..self.len()).map(|i| self.operation(i).expect("index in the validated extension"))
    }
    /// First right, then left. Antiunitarity conjugates orbital coefficients,
    /// not the REAL spatial W,t, integer cell shifts, or channel signs.
    pub fn product(&self, left: usize, right: usize) -> Option<CollinearSpinProduct> {
        let kind = self
            .operation(left)?
            .kind
            .compose(self.operation(right)?.kind);
        let n = self.unitary.operations.len();
        let mut product = self.unitary.multiplication_table[left % n][right % n];
        if kind.is_antiunitary() {
            product.operation_index += n;
        }
        Some(product)
    }
    pub fn inverse_index(&self, index: usize) -> Option<usize> {
        let kind = self.operation(index)?.kind;
        let n = self.unitary.operations.len();
        Some(self.unitary.inverse_indices[index % n] + if kind.is_antiunitary() { n } else { 0 })
    }
}

/// Witness-bearing refusal; no failure is converted to a silently smaller group.
#[derive(Debug, Clone, thiserror::Error)]
pub enum CollinearSpinError {
    #[error("invalid {name} tolerance: {value}")]
    InvalidTolerance { name: &'static str, value: f64 },
    #[error("invalid collinear spin input: {reason}")]
    InvalidInput { reason: &'static str },
    #[error(
        "all moments are zero within the moment tolerance; no magnetic channel labelling is determined"
    )]
    AllZeroMoments,
    #[error(
        "atom {atom_index} is noncollinear: perpendicular moment residual {perpendicular_residual}"
    )]
    NonCollinear {
        atom_index: usize,
        perpendicular_residual: f64,
    },
    #[error("numerical range exceeded while {context}")]
    NumericalRange { context: &'static str },
    #[error(
        "periodic image search exceeds {MAX_IMAGES} candidates; reduce the lattice basis or position tolerance"
    )]
    ImageSearchLimit,
    #[error("more than one lattice image lies within the position tolerance")]
    AmbiguousPeriodicImage,
    #[error("structural lattice symmetry search failed")]
    StructuralSearchFailed,
    #[error("collinear spin operation set is empty")]
    Empty,
    #[error("operation {operation_index} is not a unimodular lattice isometry within tolerance")]
    InvalidRotation { operation_index: usize },
    #[error("operation {operation_index} contains a nonfinite or unresolvable translation")]
    InvalidTranslation { operation_index: usize },
    #[error(
        "operations {first_index} and {duplicate_index} are duplicates modulo lattice translations"
    )]
    Duplicate {
        first_index: usize,
        duplicate_index: usize,
    },
    #[error("operation set is missing the Preserve identity")]
    MissingIdentity,
    #[error("operation {operation_index} cannot map atom {atom_index} to a same-species atom")]
    AtomMappingMissing {
        operation_index: usize,
        atom_index: usize,
    },
    #[error("operation {operation_index} maps atom {atom_index} to multiple same-species atoms")]
    AtomMappingAmbiguous {
        operation_index: usize,
        atom_index: usize,
    },
    #[error("operation {operation_index} maps multiple atoms to target {target_index}")]
    NotPermutation {
        operation_index: usize,
        target_index: usize,
    },
    #[error(
        "operation {operation_index}, atom {atom_index}: signed moment residual {residual} exceeds tolerance"
    )]
    MomentMismatch {
        operation_index: usize,
        atom_index: usize,
        residual: f64,
    },
    #[error("product of operations {left_index} and {right_index} is absent: {product:?}")]
    NotClosed {
        left_index: usize,
        right_index: usize,
        product: CollinearSpinOperation,
    },
    #[error(
        "product of operations {left_index} and {right_index} has multiple matching representatives"
    )]
    AmbiguousProduct {
        left_index: usize,
        right_index: usize,
    },
    #[error("operation {operation_index} has no two-sided inverse")]
    MissingInverse { operation_index: usize },
    #[error("identity does not act trivially on atom {atom_index} modulo its integer translation")]
    IdentityMapping { atom_index: usize },
    #[error(
        "atom {atom_index}: mapping/cell-shift composition fails for operations {left_index} and {right_index}"
    )]
    MappingComposition {
        left_index: usize,
        right_index: usize,
        atom_index: usize,
    },
    #[error(
        "multiplication/translation cocycle is not associative at operations {left_index}, {middle_index}, {right_index}"
    )]
    NotAssociative {
        left_index: usize,
        middle_index: usize,
        right_index: usize,
    },
}

/// Discover the finite unitary collinear quotient in the original input cell.
///
/// Structural lattice rotations are obtained from the existing nonmagnetic
/// lattice search. Same-species translations are tested with a unique geometric
/// bijection, then BOTH signed-magnitude conditions `m_target = +/- m_source`
/// are tested independently of W. The complete surviving set is validated;
/// nonclosure is an error, not a request to prune operations.
///
/// All-zero (or numerically zero) moments and 2D layer cells are explicitly
/// rejected. Magnetic moments already attached to `crystal` are not implicitly
/// used: pass them as `CollinearMoments::Vectors` if desired.
pub fn detect_collinear_spin_group(
    crystal: &Crystal,
    moments: CollinearMoments<'_>,
    tolerances: CollinearSpinTolerances,
) -> Result<ValidatedCollinearSpinGroup, CollinearSpinError> {
    let prepared = Prepared::new(crystal, moments, tolerances)?;
    // Reject geometrically ambiguous sites before using any of them as anchors.
    prepared.map_atoms(&CollinearSpinOperation::identity(), 0)?;
    let rotations = get_lattice_symmetry(&prepared.cell, tolerances.position, -1.0);
    if rotations.is_empty() {
        return Err(CollinearSpinError::StructuralSearchFailed);
    }
    let types = crystal.types();
    let anchor = (0..types.len())
        .min_by_key(|&i| types.iter().filter(|&&t| t == types[i]).count())
        .ok_or(CollinearSpinError::InvalidInput {
            reason: "empty structure",
        })?;
    let mut operations = Vec::new();
    for rotation in rotations.rot {
        for target in 0..types.len() {
            if types[target] != types[anchor] {
                continue;
            }
            let affine = bounded_affine(
                &rotation,
                prepared.cell.position[anchor],
                [0.; 3],
                prepared.cell.position[target],
            );
            if !prepared.within_position(&affine.map(|v| Bounded {
                value: 0.,
                error: v.error * 16.,
            }))? {
                return Err(CollinearSpinError::NumericalRange {
                    context: "constructing structural candidates",
                });
            }
            let translation = affine.map(|v| (-v.value).rem_euclid(1.0));
            let mut operation = CollinearSpinOperation {
                rotation,
                translation,
                channel_action: SpinChannelAction::Preserve,
            };
            let mapping = match prepared.map_atoms(&operation, operations.len()) {
                Ok(mapping) => mapping,
                Err(CollinearSpinError::AtomMappingMissing { .. }) => continue,
                Err(error) => return Err(error),
            };
            for action in [SpinChannelAction::Preserve, SpinChannelAction::Swap] {
                operation.channel_action = action;
                match prepared.check_moments(&operation, &mapping, operations.len()) {
                    Ok(()) => operations.push(operation),
                    Err(CollinearSpinError::MomentMismatch { .. }) => {}
                    Err(error) => return Err(error),
                }
            }
        }
    }
    prepared.validate(&operations)
}

/// Validate an explicit target group, retaining every operation and its order.
///
/// Checks identity, uniqueness, all ordered Seitz/channel products, inverses,
/// associativity, structural/moment compatibility, and the original-frame atom
/// permutation/integer-shift composition. Approximate Seitz discrepancies are
/// returned in the product witnesses; they are never absorbed into cell shifts.
/// Passing only generators is an error: supply the complete target operation set.
pub fn validate_collinear_spin_group(
    crystal: &Crystal,
    moments: CollinearMoments<'_>,
    operations: &[CollinearSpinOperation],
    tolerances: CollinearSpinTolerances,
) -> Result<ValidatedCollinearSpinGroup, CollinearSpinError> {
    Prepared::new(crystal, moments, tolerances)?.validate(operations)
}

struct Prepared<'a> {
    crystal: &'a Crystal,
    cell: Cell,
    inverse_row_bounds: Vec3,
    positions: Vec<SplitVector>,
    signed_moments: Vec<f64>,
    moment_errors: Vec<f64>,
    spin_axis: Option<Vec3>,
    collinearity_residuals: Vec<f64>,
    tolerances: CollinearSpinTolerances,
}

impl<'a> Prepared<'a> {
    fn new(
        crystal: &'a Crystal,
        moments: CollinearMoments<'_>,
        tolerances: CollinearSpinTolerances,
    ) -> Result<Self, CollinearSpinError> {
        for (name, value) in [
            ("position", tolerances.position),
            ("moment", tolerances.moment),
            ("collinearity", tolerances.collinearity),
        ] {
            if !value.is_finite() || value < 0.0 || (name == "position" && value == 0.0) {
                return Err(CollinearSpinError::InvalidTolerance { name, value });
            }
        }
        if crystal.natom() == 0 || crystal.aperiodic_axis().is_some() {
            return Err(CollinearSpinError::InvalidInput {
                reason: "need a nonempty 3D periodic cell",
            });
        }
        let lattice = crystal.lattice();
        if lattice.iter().flatten().any(|x| !x.is_finite())
            || crystal
                .positions()
                .iter()
                .flatten()
                .any(|x| !x.is_finite() || x.abs() > MAX_CELL_INTEGER)
        {
            return Err(CollinearSpinError::InvalidInput {
                reason: "nonfinite lattice or unresolvable positions",
            });
        }
        let determinant = mat_get_determinant_d3(&lattice);
        if !determinant.is_finite() || determinant == 0.0 {
            return Err(CollinearSpinError::InvalidInput {
                reason: "singular or unresolvable lattice",
            });
        }
        let inverse_lattice =
            mat_inverse_matrix_d3(&lattice, 0.0).map_err(|_| CollinearSpinError::InvalidInput {
                reason: "singular lattice",
            })?;
        if inverse_lattice.iter().flatten().any(|x| !x.is_finite())
            || mat_get_metric(&lattice)
                .iter()
                .flatten()
                .any(|x| !x.is_finite())
        {
            return Err(CollinearSpinError::InvalidInput {
                reason: "unresolvable lattice metric",
            });
        }
        let (signed_moments, moment_errors, spin_axis, collinearity_residuals) = match moments {
            CollinearMoments::Signed(values) => {
                if values.len() != crystal.natom() || values.iter().any(|v| !v.is_finite()) {
                    return Err(CollinearSpinError::InvalidInput {
                        reason: "invalid signed moment array",
                    });
                }
                (
                    values.to_vec(),
                    vec![0.; values.len()],
                    None,
                    vec![0.; values.len()],
                )
            }
            CollinearMoments::Vectors(values) => {
                if values.len() != crystal.natom()
                    || values.iter().flatten().any(|v| !v.is_finite())
                {
                    return Err(CollinearSpinError::InvalidInput {
                        reason: "invalid vector moment array",
                    });
                }
                let measured: Vec<_> = values
                    .iter()
                    .map(|v| bounded_norm(v.map(Bounded::exact)))
                    .collect();
                let lengths: Vec<_> = measured.iter().map(|v| v.value).collect();
                if measured
                    .iter()
                    .any(|v| !v.value.is_finite() || !v.error.is_finite())
                {
                    return Err(CollinearSpinError::NumericalRange {
                        context: "measuring vector moments",
                    });
                }
                let largest = (0..lengths.len())
                    .max_by(|&i, &j| lengths[i].total_cmp(&lengths[j]))
                    .ok_or(CollinearSpinError::AllZeroMoments)?;
                if lengths[largest] <= tolerances.moment {
                    return Err(CollinearSpinError::AllZeroMoments);
                }
                if lengths[largest].is_subnormal() {
                    return Err(CollinearSpinError::NumericalRange {
                        context: "normalizing subnormal vector moments",
                    });
                }
                let mut axis = values[largest].map(|v| v / lengths[largest]);
                // Fix the global up/down naming without using any spatial rotation.
                if axis.iter().find(|&&v| v != 0.0).is_some_and(|&v| v < 0.0) {
                    axis = axis.map(|v| -v);
                }
                // Enclose normalization of the largest original vector, not
                // merely its rounded axis. Pure coordinate axes are exact.
                let pure_axis = values[largest].iter().filter(|&&v| v != 0.).count() == 1;
                let bounded_axis = std::array::from_fn(|j| Bounded {
                    value: axis[j],
                    error: if pure_axis || values[largest][j] == 0. {
                        0.
                    } else {
                        (axis[j].abs() * 16. * f64::EPSILON).next_up()
                    },
                });
                let mut signed = Vec::with_capacity(values.len());
                let mut residuals = Vec::with_capacity(values.len());
                for (atom_index, value) in values.iter().enumerate() {
                    let projection = bounded_dot(value.map(Bounded::exact), bounded_axis);
                    let residual = bounded_norm(std::array::from_fn(|j| {
                        Bounded::exact(value[j]).add(projection.mul(bounded_axis[j]).neg())
                    }));
                    if !projection.value.is_finite()
                        || !projection.error.is_finite()
                        || !residual.value.is_finite()
                        || !residual.error.is_finite()
                    {
                        return Err(CollinearSpinError::NumericalRange {
                            context: "projecting vector moments",
                        });
                    }
                    if (residual.value - residual.error).next_down() > tolerances.collinearity {
                        return Err(CollinearSpinError::NonCollinear {
                            atom_index,
                            perpendicular_residual: residual.value,
                        });
                    }
                    if !residual.within(tolerances.collinearity)
                        || (lengths[atom_index] != 0. && projection.value.abs() <= projection.error)
                    {
                        return Err(CollinearSpinError::NumericalRange {
                            context: "resolving vector collinearity or channel sign",
                        });
                    }
                    // Collinearity tolerance must not hide a size mismatch.
                    signed.push(lengths[atom_index].copysign(projection.value));
                    residuals.push(residual.value);
                }
                (
                    signed,
                    measured.iter().map(|v| v.error).collect(),
                    Some(axis),
                    residuals,
                )
            }
        };
        if signed_moments.iter().all(|m| m.abs() <= tolerances.moment) {
            return Err(CollinearSpinError::AllZeroMoments);
        }
        // Never add a small residual to a large integer site representative.
        // Keep integral transports exact and do geometry on centred fractions.
        let positions: Vec<_> = crystal
            .positions()
            .iter()
            .map(|&r| SplitVector::new(r))
            .collect::<Result<_, _>>()?;
        let fractions: Vec<_> = positions.iter().map(|r| r.fractional).collect();
        let inverse_row_bounds = inverse_row_bounds(&lattice, &inverse_lattice)?;
        let mut cell = Cell::new(crystal.natom(), TensorRank::NoSpin);
        cell.set_cell(&lattice, &fractions, crystal.types())
            .map_err(|_| CollinearSpinError::InvalidInput {
                reason: "invalid structural arrays",
            })?;
        Ok(Self {
            crystal,
            cell,
            inverse_row_bounds,
            positions,
            signed_moments,
            moment_errors,
            spin_axis,
            collinearity_residuals,
            tolerances,
        })
    }

    /// All integer images within the Cartesian tolerance, not componentwise
    /// fractional nearest-integer rounding (which is wrong for skew cells).
    fn image(&self, difference: Vec3) -> Result<Option<([i64; 3], Vec3)>, CollinearSpinError> {
        self.image_bounded(difference.map(Bounded::exact))
    }

    fn image_bounded(
        &self,
        difference: [Bounded; 3],
    ) -> Result<Option<([i64; 3], Vec3)>, CollinearSpinError> {
        let mut bounds = [[0_i64; 2]; 3];
        let mut count = 1_u64;
        for j in 0..3 {
            let radius = up_add(
                (self.inverse_row_bounds[j] * self.tolerances.position).next_up(),
                difference[j].error,
            );
            let lo = integer((difference[j].value - radius).next_down().ceil())?;
            let hi = integer((difference[j].value + radius).next_up().floor())?;
            if hi < lo {
                return Ok(None);
            }
            count = count
                .checked_mul((hi - lo + 1) as u64)
                .ok_or(CollinearSpinError::ImageSearchLimit)?;
            if count > MAX_IMAGES {
                return Err(CollinearSpinError::ImageSearchLimit);
            }
            bounds[j] = [lo, hi];
        }
        // ponytail: bounded 3D enumeration; reduce pathological bases before
        // raising the ceiling or introducing a general closest-vector solver.
        let mut found = None;
        for x in bounds[0][0]..=bounds[0][1] {
            for y in bounds[1][0]..=bounds[1][1] {
                for z in bounds[2][0]..=bounds[2][1] {
                    let shift = [x, y, z];
                    let residual = std::array::from_fn(|j| {
                        difference[j].add(Bounded::exact(-(shift[j] as f64)))
                    });
                    if self.within_position(&residual)? {
                        if found.is_some() {
                            return Err(CollinearSpinError::AmbiguousPeriodicImage);
                        }
                        found = Some((shift, residual.map(|v| v.value)));
                    }
                }
            }
        }
        Ok(found)
    }

    fn within_position(&self, value: &[Bounded; 3]) -> Result<bool, CollinearSpinError> {
        let cart = self
            .cell
            .lattice
            .map(|row| bounded_dot(row.map(Bounded::exact), *value));
        let length = bounded_norm(cart);
        if !length.value.is_finite() || !length.error.is_finite() {
            return Err(CollinearSpinError::NumericalRange {
                context: "measuring a geometric residual",
            });
        }
        if length.within(self.tolerances.position) {
            return Ok(true);
        }
        if (length.value - length.error).next_down() > self.tolerances.position {
            return Ok(false);
        }
        Err(CollinearSpinError::NumericalRange {
            context: "resolving distance at the tolerance boundary",
        })
    }

    fn map_atoms(
        &self,
        operation: &CollinearSpinOperation,
        operation_index: usize,
    ) -> Result<Vec<CollinearAtomMapping>, CollinearSpinError> {
        let mut mapping = Vec::with_capacity(self.crystal.natom());
        let mut used = vec![false; self.crystal.natom()];
        let translation = SplitVector::new(operation.translation)?;
        for (atom_index, position) in self.positions.iter().enumerate() {
            let integral = rotate_shift(&operation.rotation, position.integral);
            let mut found = None;
            for (target, target_position) in self.positions.iter().enumerate() {
                if self.crystal.types()[atom_index] != self.crystal.types()[target] {
                    continue;
                }
                let difference = bounded_affine(
                    &operation.rotation,
                    position.fractional,
                    translation.fractional,
                    target_position.fractional,
                );
                if let Some((image_shift, position_residual)) = self.image_bounded(difference)? {
                    let cell_shift = checked_shift(std::array::from_fn(|j| {
                        integral[j] + i128::from(translation.integral[j])
                            - i128::from(target_position.integral[j])
                            + i128::from(image_shift[j])
                    }))?;
                    if found.is_some() {
                        return Err(CollinearSpinError::AtomMappingAmbiguous {
                            operation_index,
                            atom_index,
                        });
                    }
                    found = Some(CollinearAtomMapping {
                        target,
                        cell_shift,
                        position_residual,
                        moment_residual: 0.,
                    });
                }
            }
            let found = found.ok_or(CollinearSpinError::AtomMappingMissing {
                operation_index,
                atom_index,
            })?;
            if used[found.target] {
                return Err(CollinearSpinError::NotPermutation {
                    operation_index,
                    target_index: found.target,
                });
            }
            used[found.target] = true;
            mapping.push(found);
        }
        Ok(mapping)
    }

    fn check_moments(
        &self,
        operation: &CollinearSpinOperation,
        mapping: &[CollinearAtomMapping],
        operation_index: usize,
    ) -> Result<(), CollinearSpinError> {
        for (atom_index, atom) in mapping.iter().enumerate() {
            if atom.target == atom_index && operation.channel_action == SpinChannelAction::Preserve
            {
                continue; // The same input magnitude, not two independent norm errors.
            }
            let source = Bounded {
                value: self.signed_moments[atom_index],
                error: self.moment_errors[atom_index],
            };
            let target = Bounded {
                value: self.signed_moments[atom.target],
                error: self.moment_errors[atom.target],
            };
            let residual =
                target.add(source.mul(Bounded::exact(-f64::from(operation.channel_action.sign()))));
            if !residual.value.is_finite()
                || (residual.value.abs() - residual.error).next_down() > self.tolerances.moment
            {
                return Err(CollinearSpinError::MomentMismatch {
                    operation_index,
                    atom_index,
                    residual: residual.value,
                });
            }
            if !residual.error.is_finite() || !residual.within(self.tolerances.moment) {
                return Err(CollinearSpinError::NumericalRange {
                    context: "comparing vector moment magnitudes",
                });
            }
        }
        Ok(())
    }

    fn validate_rotation(
        &self,
        rotation: &Mat3I,
        operation_index: usize,
    ) -> Result<(), CollinearSpinError> {
        let det = crate::mathfunc::mat_get_determinant_i3_wide(rotation);
        if det != 1 && det != -1 {
            return Err(CollinearSpinError::InvalidRotation { operation_index });
        }
        let columns: [[Bounded; 3]; 3] =
            std::array::from_fn(|j| self.cell.lattice.map(|row| Bounded::exact(row[j])));
        let transformed: [[Bounded; 3]; 3] = std::array::from_fn(|j| {
            self.cell.lattice.map(|row| {
                bounded_dot(
                    row.map(Bounded::exact),
                    rotation.map(|r| Bounded::exact(f64::from(r[j]))),
                )
            })
        });
        let lengths = columns.map(|c| norm(&c.map(|v| v.value)));
        for i in 0..3 {
            for j in 0..3 {
                let bound = self.tolerances.position * (lengths[i] + lengths[j])
                    + self.tolerances.position.powi(2);
                let residual = bounded_dot(columns[i], columns[j])
                    .add(bounded_dot(transformed[i], transformed[j]).neg());
                if !bound.is_finite() || !residual.value.is_finite() || !residual.error.is_finite()
                {
                    return Err(CollinearSpinError::NumericalRange {
                        context: "certifying the lattice metric",
                    });
                }
                if up_add(residual.value.abs(), residual.error)
                    <= (bound * (1. - 16. * f64::EPSILON)).next_down()
                {
                    continue;
                }
                if (residual.value.abs() - residual.error).next_down()
                    > (bound * (1. + 16. * f64::EPSILON)).next_up()
                {
                    return Err(CollinearSpinError::InvalidRotation { operation_index });
                }
                return Err(CollinearSpinError::NumericalRange {
                    context: "certifying the lattice metric",
                });
            }
        }
        Ok(())
    }

    fn validate(
        self,
        operations: &[CollinearSpinOperation],
    ) -> Result<ValidatedCollinearSpinGroup, CollinearSpinError> {
        if operations.is_empty() {
            return Err(CollinearSpinError::Empty);
        }
        let mut mapped = Vec::with_capacity(operations.len());
        let mut identity = None;
        for (operation_index, operation) in operations.iter().enumerate() {
            if operation
                .translation
                .iter()
                .any(|v| !v.is_finite() || v.abs() > MAX_CELL_INTEGER)
            {
                return Err(CollinearSpinError::InvalidTranslation { operation_index });
            }
            self.validate_rotation(&operation.rotation, operation_index)?;
            let translation = SplitVector::new(operation.translation)?;
            for (first_index, earlier) in operations[..operation_index].iter().enumerate() {
                let earlier_translation = SplitVector::new(earlier.translation)?;
                if earlier.rotation == operation.rotation
                    && earlier.channel_action == operation.channel_action
                    && self
                        .image_bounded(bounded_affine(
                            &IDENTITY,
                            [0.; 3],
                            earlier_translation.fractional,
                            translation.fractional,
                        ))?
                        .is_some()
                {
                    return Err(CollinearSpinError::Duplicate {
                        first_index,
                        duplicate_index: operation_index,
                    });
                }
            }
            let mut atom_mapping = self.map_atoms(operation, operation_index)?;
            self.check_moments(operation, &atom_mapping, operation_index)?;
            for (i, atom) in atom_mapping.iter_mut().enumerate() {
                atom.moment_residual = self.signed_moments[atom.target]
                    - f64::from(operation.channel_action.sign()) * self.signed_moments[i];
            }
            if operation.rotation == IDENTITY
                && operation.channel_action == SpinChannelAction::Preserve
                && let Some((image_shift, _)) = self.image(translation.fractional)?
            {
                let shift = checked_shift(std::array::from_fn(|j| {
                    i128::from(translation.integral[j]) + i128::from(image_shift[j])
                }))?;
                for (atom_index, atom) in atom_mapping.iter().enumerate() {
                    if atom.target != atom_index || atom.cell_shift != shift {
                        return Err(CollinearSpinError::IdentityMapping { atom_index });
                    }
                }
                identity = Some(operation_index);
            }
            mapped.push(MappedCollinearSpinOperation {
                operation: *operation,
                atom_mapping,
            });
        }
        let identity_index = identity.ok_or(CollinearSpinError::MissingIdentity)?;
        let mut table = Vec::with_capacity(operations.len());
        // ponytail: exhaustive finite-set certification, O(|G|^3 + |G|^2 N);
        // index rotations if measured large-supercell workloads need it.
        for (left_index, left) in operations.iter().enumerate() {
            let mut row = Vec::with_capacity(operations.len());
            for (right_index, right) in operations.iter().enumerate() {
                let product = compose(left, right)?;
                let mut found = None;
                for (operation_index, candidate) in operations.iter().enumerate() {
                    if candidate.rotation != product.operation.rotation
                        || candidate.channel_action != product.operation.channel_action
                    {
                        continue;
                    }
                    let candidate_translation = SplitVector::new(candidate.translation)?;
                    if let Some((image_shift, translation_residual)) =
                        self.image_bounded(std::array::from_fn(|j| {
                            Bounded {
                                value: product.operation.translation[j],
                                error: product.translation_error[j],
                            }
                            .add(Bounded::exact(-candidate_translation.fractional[j]))
                        }))?
                    {
                        let cell_shift = checked_shift(std::array::from_fn(|j| {
                            product.integral[j] - i128::from(candidate_translation.integral[j])
                                + i128::from(image_shift[j])
                        }))?;
                        if found.is_some() {
                            return Err(CollinearSpinError::AmbiguousProduct {
                                left_index,
                                right_index,
                            });
                        }
                        found = Some(CollinearSpinProduct {
                            operation_index,
                            cell_shift,
                            translation_residual,
                        });
                    }
                }
                let product = found.ok_or(CollinearSpinError::NotClosed {
                    left_index,
                    right_index,
                    product: product.operation,
                })?;
                self.check_mapping_product(&mapped, left_index, right_index, product)?;
                row.push(product);
            }
            table.push(row);
        }
        let mut inverse_indices = Vec::with_capacity(operations.len());
        for (i, row) in table.iter().enumerate() {
            if table[identity_index][i].operation_index != i
                || row[identity_index].operation_index != i
            {
                return Err(CollinearSpinError::MissingIdentity);
            }
            let inverse = (0..operations.len())
                .find(|&j| {
                    row[j].operation_index == identity_index
                        && table[j][i].operation_index == identity_index
                })
                .ok_or(CollinearSpinError::MissingInverse { operation_index: i })?;
            inverse_indices.push(inverse);
        }
        for g in 0..operations.len() {
            for h in 0..operations.len() {
                for k in 0..operations.len() {
                    let gh = table[g][h];
                    let hk = table[h][k];
                    let lhs = table[gh.operation_index][k];
                    let rhs = table[g][hk.operation_index];
                    let rotated_shift = rotate_shift(&operations[g].rotation, hk.cell_shift);
                    if lhs.operation_index != rhs.operation_index
                        || (0..3).any(|j| {
                            i128::from(gh.cell_shift[j]) + i128::from(lhs.cell_shift[j])
                                != rotated_shift[j] + i128::from(rhs.cell_shift[j])
                        })
                    {
                        return Err(CollinearSpinError::NotAssociative {
                            left_index: g,
                            middle_index: h,
                            right_index: k,
                        });
                    }
                }
            }
        }
        Ok(ValidatedCollinearSpinGroup {
            operations: mapped,
            multiplication_table: table,
            identity_index,
            inverse_indices,
            tolerances: self.tolerances,
            signed_moments: self.signed_moments,
            spin_axis: self.spin_axis,
            collinearity_residuals: self.collinearity_residuals,
        })
    }

    fn check_mapping_product(
        &self,
        mapped: &[MappedCollinearSpinOperation],
        g: usize,
        h: usize,
        product: CollinearSpinProduct,
    ) -> Result<(), CollinearSpinError> {
        for (atom_index, hi) in mapped[h].atom_mapping.iter().enumerate() {
            let gi = mapped[g].atom_mapping[hi.target];
            let pi = mapped[product.operation_index].atom_mapping[atom_index];
            let shift = rotate_shift(&mapped[g].operation.rotation, hi.cell_shift);
            let residual = bounded_affine(
                &mapped[g].operation.rotation,
                hi.position_residual,
                gi.position_residual,
                pi.position_residual,
            );
            let residual = std::array::from_fn(|j| {
                residual[j].add(Bounded::exact(-product.translation_residual[j]))
            });
            if gi.target != pi.target
                || (0..3).any(|j| {
                    i128::from(gi.cell_shift[j]) + shift[j]
                        != i128::from(pi.cell_shift[j]) + i128::from(product.cell_shift[j])
                })
                || !self.within_position(&residual)?
            {
                return Err(CollinearSpinError::MappingComposition {
                    left_index: g,
                    right_index: h,
                    atom_index,
                });
            }
        }
        Ok(())
    }
}

fn norm(v: &Vec3) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

// Local roundoff witnesses, not a general interval/matrix framework. FMA
// encloses multiplication error and TwoSum encloses addition cancellation.
#[derive(Clone, Copy)]
struct Bounded {
    value: f64,
    error: f64,
}

fn up_add(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}

impl Bounded {
    fn exact(value: f64) -> Self {
        Self { value, error: 0. }
    }
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            ..self
        }
    }
    fn within(self, tolerance: f64) -> bool {
        if self.error == 0. {
            self.value.abs() <= tolerance
        } else {
            up_add(self.value.abs(), self.error) <= tolerance
        }
    }
    fn add(self, right: Self) -> Self {
        if self.value == 0. && self.error == 0. {
            return right;
        }
        if right.value == 0. && right.error == 0. {
            return self;
        }
        if self.value == -right.value && self.error == 0. && right.error == 0. {
            return Self::exact(0.);
        }
        let value = self.value + right.value;
        let virtual_right = value - self.value;
        let remainder = (self.value - (value - virtual_right)) + (right.value - virtual_right);
        Self {
            value,
            error: up_add(up_add(self.error, right.error), remainder.abs().next_up()),
        }
    }
    fn mul(self, right: Self) -> Self {
        if (self.value == 0. && self.error == 0.) || (right.value == 0. && right.error == 0.) {
            return Self::exact(0.);
        }
        if self.error == 0. && self.value.abs() == 1. {
            return if self.value > 0. { right } else { right.neg() };
        }
        if right.error == 0. && right.value.abs() == 1. {
            return if right.value > 0. { self } else { self.neg() };
        }
        let value = self.value * right.value;
        let remainder = self.value.mul_add(right.value, -value).abs().next_up();
        let propagated = up_add(
            (self.value.abs() * right.error).next_up(),
            (right.value.abs() * self.error).next_up(),
        );
        Self {
            value,
            error: up_add(
                remainder,
                up_add(propagated, (self.error * right.error).next_up()),
            ),
        }
    }
}

fn bounded_dot(left: [Bounded; 3], right: [Bounded; 3]) -> Bounded {
    left.into_iter()
        .zip(right)
        .fold(Bounded::exact(0.), |sum, (a, b)| sum.add(a.mul(b)))
}

fn bounded_norm(value: [Bounded; 3]) -> Bounded {
    if value.iter().all(|v| v.value.is_finite() && v.error == 0.)
        && value.iter().filter(|v| v.value != 0.).count() <= 1
    {
        return Bounded::exact(value.iter().map(|v| v.value.abs()).fold(0_f64, f64::max));
    }
    let length = norm(&value.map(|v| v.value));
    Bounded {
        value: length,
        error: up_add(
            norm(&value.map(|v| v.error)) * (1. + 8. * f64::EPSILON),
            (length * 8. * f64::EPSILON).next_up(),
        ),
    }
}

fn bounded_affine(
    rotation: &Mat3I,
    position: Vec3,
    translation: Vec3,
    target: Vec3,
) -> [Bounded; 3] {
    std::array::from_fn(|j| {
        bounded_dot(
            rotation[j].map(|v| Bounded::exact(f64::from(v))),
            position.map(Bounded::exact),
        )
        .add(Bounded::exact(translation[j]))
        .add(Bounded::exact(-target[j]))
    })
}

fn integer(value: f64) -> Result<i64, CollinearSpinError> {
    if !value.is_finite() || value.abs() > MAX_CELL_INTEGER {
        return Err(CollinearSpinError::NumericalRange {
            context: "resolving integer cell shifts",
        });
    }
    Ok(value as i64)
}

fn rotate_shift(rotation: &Mat3I, shift: [i64; 3]) -> [i128; 3] {
    rotation.map(|row| {
        (0..3)
            .map(|j| i128::from(row[j]) * i128::from(shift[j]))
            .sum()
    })
}

/// Exact integer transport plus a small fractional part. Nearest-integer
/// subtraction is exact for finite input f64s in range (Sterbenz's lemma).
#[derive(Clone, Copy)]
struct SplitVector {
    integral: [i64; 3],
    fractional: Vec3,
}

impl SplitVector {
    fn new(value: Vec3) -> Result<Self, CollinearSpinError> {
        let mut integral = [0; 3];
        for j in 0..3 {
            integral[j] = integer(value[j].round())?;
        }
        Ok(Self {
            integral,
            fractional: std::array::from_fn(|j| value[j] - integral[j] as f64),
        })
    }
}

fn checked_shift(value: [i128; 3]) -> Result<[i64; 3], CollinearSpinError> {
    let mut shift = [0; 3];
    for (dst, component) in shift.iter_mut().zip(value) {
        *dst = i64::try_from(component).map_err(|_| CollinearSpinError::NumericalRange {
            context: "combining exact integer cell shifts",
        })?;
    }
    Ok(shift)
}

// Certify B A = I-E with ||E||_infinity < 1/4, including a conservative
// floating dot-product error bound. A finite cofactor inverse alone does NOT
// establish invertibility: a singular large matrix can have a nonzero float det.
// Neumann's bound then encloses each true inverse row, so image intervals do
// not silently shrink when the computed inverse is slightly inaccurate.
fn inverse_row_bounds(lattice: &Mat3, inverse: &Mat3) -> Result<Vec3, CollinearSpinError> {
    let mut error = 0_f64;
    for (i, row) in inverse.iter().enumerate() {
        let mut row_error = 0.;
        for (j, _) in lattice[0].iter().enumerate() {
            let terms = std::array::from_fn::<_, 3, _>(|k| row[k] * lattice[k][j]);
            let sum: f64 = terms.iter().sum();
            let magnitude: f64 = terms.iter().map(|v| v.abs()).sum();
            row_error += (sum - f64::from(i == j)).abs()
                + 16. * f64::EPSILON * (magnitude + 1.)
                + 16. * f64::MIN_POSITIVE;
        }
        if !row_error.is_finite() || row_error >= 0.25 {
            return Err(CollinearSpinError::InvalidInput {
                reason: "lattice inverse cannot be certified",
            });
        }
        error = error.max(row_error);
    }
    let rows = inverse.map(|row| norm(&row) * (1. + 8. * f64::EPSILON));
    let largest = rows.into_iter().fold(0_f64, f64::max);
    let padding = error / (1. - error) * largest;
    let bounds = rows.map(|n| (n + padding) * (1. + 8. * f64::EPSILON));
    if bounds.iter().any(|n| !n.is_finite()) {
        return Err(CollinearSpinError::NumericalRange {
            context: "bounding the inverse lattice",
        });
    }
    Ok(bounds)
}

struct ComposedOperation {
    /// Translation is the centred fractional representative, not a large sum.
    operation: CollinearSpinOperation,
    integral: [i128; 3],
    translation_error: Vec3,
}

fn compose(
    left: &CollinearSpinOperation,
    right: &CollinearSpinOperation,
) -> Result<ComposedOperation, CollinearSpinError> {
    let mut rotation = [[0; 3]; 3];
    for (i, row) in rotation.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            let sum: i128 = (0..3)
                .map(|k| i128::from(left.rotation[i][k]) * i128::from(right.rotation[k][j]))
                .sum();
            *value = i32::try_from(sum).map_err(|_| CollinearSpinError::NumericalRange {
                context: "composing integer rotations",
            })?;
        }
    }
    let left_translation = SplitVector::new(left.translation)?;
    let right_translation = SplitVector::new(right.translation)?;
    let raw_translation = bounded_affine(
        &left.rotation,
        right_translation.fractional,
        left_translation.fractional,
        [0.; 3],
    );
    let translation = SplitVector::new(raw_translation.map(|v| v.value))?;
    let rotated_integral = rotate_shift(&left.rotation, right_translation.integral);
    Ok(ComposedOperation {
        operation: CollinearSpinOperation {
            rotation,
            translation: translation.fractional,
            channel_action: left.channel_action.compose(right.channel_action),
        },
        integral: std::array::from_fn(|j| {
            i128::from(left_translation.integral[j])
                + rotated_integral[j]
                + i128::from(translation.integral[j])
        }),
        translation_error: raw_translation.map(|v| v.error),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composition_overflow_is_typed_but_wide_intermediate_cancellation_is_supported() {
        let shear = CollinearSpinOperation {
            rotation: [[1, i32::MAX, 0], [0, 1, 0], [0, 0, 1]],
            ..CollinearSpinOperation::identity()
        };
        assert!(matches!(
            compose(&shear, &shear),
            Err(CollinearSpinError::NumericalRange {
                context: "composing integer rotations"
            })
        ));
        let n = 46_340;
        let involution = CollinearSpinOperation {
            rotation: [[n, 1 - n * n, 0], [1, -n, 0], [0, 0, 1]],
            ..CollinearSpinOperation::identity()
        };
        assert_eq!(
            compose(&involution, &involution).unwrap().operation,
            CollinearSpinOperation::identity()
        );
        let large = CollinearSpinOperation {
            translation: [MAX_CELL_INTEGER, 0., 0.],
            ..CollinearSpinOperation::identity()
        };
        assert_eq!(
            compose(&large, &large).unwrap().integral,
            [1_i128 << 53, 0, 0]
        );
        assert!(checked_shift([i128::from(i64::MAX) + 1, 0, 0]).is_err());
        assert!(checked_shift([i128::from(i64::MIN) - 1, 0, 0]).is_err());
        assert!(integer(f64::NAN).is_err());
        assert!(integer(2. * MAX_CELL_INTEGER).is_err());
        assert_eq!(integer(-MAX_CELL_INTEGER).unwrap(), -(1_i64 << 52));
    }

    #[test]
    fn altered_permutation_shift_and_residual_fail_mapping_composition() {
        let crystal = Crystal::new(
            [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]],
            vec![[0.; 3], [0.5, 0., 0.]],
            vec![26, 26],
        )
        .unwrap();
        let operations = [
            CollinearSpinOperation::identity(),
            CollinearSpinOperation {
                translation: [0.5, 0., 0.],
                channel_action: SpinChannelAction::Swap,
                ..CollinearSpinOperation::identity()
            },
        ];
        let group = validate_collinear_spin_group(
            &crystal,
            CollinearMoments::Signed(&[2., -2.]),
            &operations,
            Default::default(),
        )
        .unwrap();
        let prepared = Prepared::new(
            &crystal,
            CollinearMoments::Signed(&[2., -2.]),
            Default::default(),
        )
        .unwrap();
        let product = group.multiplication_table[1][1];
        for corruption in 0..3 {
            let mut mapped = group.operations.clone();
            let atom = &mut mapped[1].atom_mapping[0];
            match corruption {
                0 => atom.target = 0,
                1 => atom.cell_shift[0] += 1,
                _ => atom.position_residual[0] += 0.01,
            }
            assert!(matches!(
                prepared.check_mapping_product(&mapped, 1, 1, product),
                Err(CollinearSpinError::MappingComposition {
                    left_index: 1,
                    right_index: 1,
                    ..
                })
            ));
        }
    }
}
