//! Spin-group symmetry independent of magnetic space groups.
//!
//! A spin-group element is written
//!
//! ```text
//! g = [A || {W|t}] T^theta,   theta in {0,1},
//! ```
//!
//! where `{W|t}` is a lattice isometry in the input cell, `A` is the **net
//! action on the magnetic order parameter** (already including the sign of
//! time reversal), and `theta` is the unitary/anti-unitary flag.  The two
//! independent pieces are never collapsed into one Boolean:
//!
//! * a unitary `Swap` exchanges the two spin blocks without conjugating orbital
//!   matrices and without an extra momentum minus;
//! * an anti-unitary `Preserve` conjugates and reverses momentum while keeping
//!   the spin block.
//!
//! This module deliberately does not use the magnetic-space-group
//! `time_reversal` flag or the axial-vector rule `det(W) W m`.  It reuses only
//! the lattice-isometry search and atom mapping machinery, and returns
//! finite representative sets modulo the continuous spin kernel appropriate to
//! each magnetic-order category.
//!
//! Three categories are exposed:
//!
//! * [`collinear`] — signed scalar moments on a common axis; the finite action
//!   is a sign (`Preserve`/`Swap`);
//! * coplanar — moments in a common plane; the finite action is an `O(2)`
//!   matrix (interface reserved);
//! * non-coplanar — generic moments; the finite action is an `O(3)` matrix
//!   (interface reserved).
//!
//! The module is the long-term home of this capability.  The older
//! `collinear_spin` entry points are retired in a follow-up cut-over once the
//! three categories and the oriented-spin-space-group database identification
//! are in place; see `docs/spin-group.md`.

use crate::OperationKind;
use crate::mathfunc::{Mat3I, Vec3};

pub mod collinear;
pub mod database;
pub use collinear::{
    CollinearAction, CollinearMoments, CollinearOrder, CollinearSpinGroup,
    detect_collinear_spin_group, validate_collinear_spin_group,
};

/// Independent tolerances for spin-group structure/moment checks.
///
/// `position` is a Cartesian length; `moment` is an absolute residual in the
/// moment units; `direction` is an angular (radian) tolerance for collinearity
/// and, later, coplanarity/spin-rotation checks.  No tolerance is borrowed
/// from a quantity with different units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpinTolerances {
    pub position: f64,
    pub moment: f64,
    pub direction: f64,
}

impl Default for SpinTolerances {
    fn default() -> Self {
        Self {
            position: 1e-5,
            moment: 1e-5,
            direction: 1e-5,
        }
    }
}

impl SpinTolerances {
    fn validate(self) -> Result<(), SpinGroupError> {
        for (name, value, strictly_positive) in [
            ("position", self.position, true),
            ("moment", self.moment, false),
            ("direction", self.direction, false),
        ] {
            if !value.is_finite() || value < 0.0 || (strictly_positive && value == 0.0) {
                return Err(SpinGroupError::InvalidTolerance { name, value });
            }
        }
        Ok(())
    }
}

/// Magnetic-order category of a validated spin group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpinGroupKind {
    Collinear,
    Coplanar,
    NonCoplanar,
}

/// One spin-group operation.
///
/// `action` is the net action on the magnetic order parameter; `kind` carries
/// only unitarity.  The spatial part is always in the original input cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpinOperation<A> {
    pub rotation: Mat3I,
    pub translation: Vec3,
    pub action: A,
    pub kind: OperationKind,
}

/// Atom mapping of one operation in the original input cell.
///
/// The position residual is the **fractional** residual of
/// `W r_source + t = r_target + cell_shift + residual`; its Cartesian norm is
/// at most `SpinTolerances::position`.  `M` is the category-specific moment
/// residual (a signed scalar for collinear, a vector for the other modes).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpinAtomMapping<M> {
    pub target: usize,
    pub cell_shift: [i64; 3],
    pub position_residual: Vec3,
    pub moment_residual: M,
}

/// Witness-bearing product `g h` (first `h`, then `g`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpinGroupProduct {
    pub operation_index: usize,
    pub cell_shift: [i64; 3],
    pub translation_residual: Vec3,
    pub action_residual: f64,
}

/// Net action on the magnetic order parameter.
///
/// Implementations must compose like the target action group, not like the
/// full spin-space group: time reversal lives in [`OperationKind`] and is not
/// part of an action.
pub trait SpinAction: Copy + core::fmt::Debug + PartialEq {
    fn identity() -> Self;
    fn compose(self, right: Self) -> Self;
    fn inverse(self) -> Self;
    /// Zero for exact actions; a tolerance-scaled residual for floating actions.
    fn residual(self, other: Self) -> f64;
}

/// Immutable finite spin-group certificate.
///
/// `A` is the category action, `O` the validated order-parameter metadata and
/// `M` the moment-residual type of the atom mappings.  Construction is only
/// through the `detect_*`/`validate_*` entry points of the category modules.
#[derive(Debug, Clone)]
pub struct SpinGroup<A: SpinAction, O, M> {
    kind: SpinGroupKind,
    operations: Vec<SpinOperation<A>>,
    atom_mappings: Vec<Vec<SpinAtomMapping<M>>>,
    multiplication_table: Vec<Vec<SpinGroupProduct>>,
    identity_index: usize,
    inverse_indices: Vec<usize>,
    tolerances: SpinTolerances,
    order: O,
}

impl<A: SpinAction, O, M> SpinGroup<A, O, M> {
    pub fn kind(&self) -> SpinGroupKind {
        self.kind
    }
    pub fn len(&self) -> usize {
        self.operations.len()
    }
    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }
    pub fn operations(&self) -> &[SpinOperation<A>] {
        &self.operations
    }
    pub fn operation(&self, index: usize) -> Option<&SpinOperation<A>> {
        self.operations.get(index)
    }
    pub fn atom_mappings(&self, index: usize) -> Option<&[SpinAtomMapping<M>]> {
        self.atom_mappings.get(index).map(Vec::as_slice)
    }
    pub fn multiplication_table(&self) -> &[Vec<SpinGroupProduct>] {
        &self.multiplication_table
    }
    /// Product `g h` (first `h`, then `g`) in the validated operation order.
    pub fn product(&self, left: usize, right: usize) -> Option<SpinGroupProduct> {
        self.multiplication_table
            .get(left)
            .and_then(|row| row.get(right))
            .copied()
    }
    pub fn identity_index(&self) -> usize {
        self.identity_index
    }
    pub fn inverse_indices(&self) -> &[usize] {
        &self.inverse_indices
    }
    pub fn inverse_index(&self, index: usize) -> Option<usize> {
        self.inverse_indices.get(index).copied()
    }
    pub fn tolerances(&self) -> SpinTolerances {
        self.tolerances
    }
    pub fn order(&self) -> &O {
        &self.order
    }
}

/// Errors of the spin-group layer.  Every failure is witness-bearing; no
/// operation is silently dropped to obtain a smaller group.
#[derive(Debug, Clone, thiserror::Error)]
pub enum SpinGroupError {
    #[error("invalid {name} tolerance: {value}")]
    InvalidTolerance { name: &'static str, value: f64 },
    #[error("invalid spin-group input: {reason}")]
    InvalidInput { reason: &'static str },
    #[error("all moments are zero within the moment tolerance; no order parameter is defined")]
    AllZeroMoments,
    #[error("atom {atom_index} is not collinear: perpendicular residual {residual}")]
    NonCollinear { atom_index: usize, residual: f64 },
    #[error("atom {atom_index} is not coplanar: out-of-plane residual {residual}")]
    NonCoplanar { atom_index: usize, residual: f64 },
    #[error("numerical range exceeded while {context}")]
    NumericalRange { context: &'static str },
    #[error("periodic image search exceeds its bound")]
    ImageSearchLimit,
    #[error("more than one lattice image lies within the position tolerance")]
    AmbiguousPeriodicImage,
    #[error("structural lattice symmetry search failed")]
    StructuralSearchFailed,
    #[error("spin operation set is empty")]
    Empty,
    #[error("operation {operation_index} is not a lattice isometry within tolerance")]
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
    #[error("operation set is missing the identity")]
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
        "operation {operation_index}, atom {atom_index}: moment residual {residual} exceeds tolerance"
    )]
    MomentMismatch {
        operation_index: usize,
        atom_index: usize,
        residual: f64,
    },
    #[error("product of spin operations {left_index} and {right_index} is absent")]
    NotClosed {
        left_index: usize,
        right_index: usize,
    },
    #[error(
        "product of spin operations {left_index} and {right_index} has multiple matching representatives"
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
        "multiplication/translation cocycle is not associative at ({left_index},{middle_index},{right_index})"
    )]
    NotAssociative {
        left_index: usize,
        middle_index: usize,
        right_index: usize,
    },
    #[error(
        "collinear group requires equally many unitary and anti-unitary operations (found {unitary} and {antiunitary})"
    )]
    KindCount { unitary: usize, antiunitary: usize },
    #[error(
        "anti-unitary operation {operation_index} has no unitary spatial partner in the target set"
    )]
    AntiunitaryPartner { operation_index: usize },
    #[error(
        "operation {operation_index}: spin action is not orthogonal within tolerance (residual {residual})"
    )]
    NonOrthogonalSpinAction {
        operation_index: usize,
        residual: f64,
    },
    #[error("ambiguous magnetic-order classification: {reason}")]
    AmbiguousClassification { reason: &'static str },
}

impl From<crate::collinear_spin::CollinearSpinError> for SpinGroupError {
    fn from(error: crate::collinear_spin::CollinearSpinError) -> Self {
        use crate::collinear_spin::CollinearSpinError as E;
        match error {
            E::InvalidTolerance { name, value } => Self::InvalidTolerance { name, value },
            E::InvalidInput { reason } => Self::InvalidInput { reason },
            E::AllZeroMoments => Self::AllZeroMoments,
            E::NonCollinear {
                atom_index,
                perpendicular_residual,
            } => Self::NonCollinear {
                atom_index,
                residual: perpendicular_residual,
            },
            E::NumericalRange { context } => Self::NumericalRange { context },
            E::ImageSearchLimit => Self::ImageSearchLimit,
            E::AmbiguousPeriodicImage => Self::AmbiguousPeriodicImage,
            E::StructuralSearchFailed => Self::StructuralSearchFailed,
            E::Empty => Self::Empty,
            E::InvalidRotation { operation_index } => Self::InvalidRotation { operation_index },
            E::InvalidTranslation { operation_index } => {
                Self::InvalidTranslation { operation_index }
            }
            E::Duplicate {
                first_index,
                duplicate_index,
            } => Self::Duplicate {
                first_index,
                duplicate_index,
            },
            E::MissingIdentity => Self::MissingIdentity,
            E::AtomMappingMissing {
                operation_index,
                atom_index,
            } => Self::AtomMappingMissing {
                operation_index,
                atom_index,
            },
            E::AtomMappingAmbiguous {
                operation_index,
                atom_index,
            } => Self::AtomMappingAmbiguous {
                operation_index,
                atom_index,
            },
            E::NotPermutation {
                operation_index,
                target_index,
            } => Self::NotPermutation {
                operation_index,
                target_index,
            },
            E::MomentMismatch {
                operation_index,
                atom_index,
                residual,
            } => Self::MomentMismatch {
                operation_index,
                atom_index,
                residual,
            },
            E::NotClosed {
                left_index,
                right_index,
                ..
            } => Self::NotClosed {
                left_index,
                right_index,
            },
            E::AmbiguousProduct {
                left_index,
                right_index,
            } => Self::AmbiguousProduct {
                left_index,
                right_index,
            },
            E::MissingInverse { operation_index } => Self::MissingInverse { operation_index },
            E::IdentityMapping { atom_index } => Self::IdentityMapping { atom_index },
            E::MappingComposition {
                left_index,
                right_index,
                atom_index,
            } => Self::MappingComposition {
                left_index,
                right_index,
                atom_index,
            },
            E::NotAssociative {
                left_index,
                middle_index,
                right_index,
            } => Self::NotAssociative {
                left_index,
                middle_index,
                right_index,
            },
        }
    }
}
