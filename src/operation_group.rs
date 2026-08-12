//! Validated magnetic operation sets and spin-$1/2$ lifts.
//!
//! This module is the boundary between a numerical Hamiltonian symmetry test
//! and the crystallographic database machinery.  A caller first filters the
//! structural operations by the Hamiltonian covariance equation, then creates
//! a [`ValidatedMagneticOperationSet`].  Validation is deliberately strict:
//! database identification is attempted only after the surviving operations
//! have been proved to form a finite magnetic group modulo lattice
//! translations.

use crate::api::{SymmetryOp, SymmetryOps};
use crate::irrep::wigner::{SeitzOp, compose_seitz};
use crate::mathfunc::{Mat3, Mat3I, Vec3, mat_get_determinant_d3};
use crate::symmetry::MagneticSymmetry;
use crate::{MagneticType, SymError};

const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// A witness-bearing failure while validating or identifying a magnetic
/// operation set.
///
/// Unlike the low-level [`SymError`] codes, these variants preserve the
/// offending operation indices (and, for closure, the missing product).  This
/// makes a failed Hamiltonian covariance check diagnosable without parsing a
/// log or mapping the failure to a sentinel group.
#[derive(Debug, Clone, thiserror::Error)]
pub enum MagneticOperationSetError {
    /// The comparison tolerance must be finite and lie below half a lattice
    /// period, so modular equality remains separating.
    #[error("operation-set tolerance must be finite and lie in (0, 0.5), got {tolerance}")]
    InvalidTolerance { tolerance: f64 },
    /// A magnetic group cannot be represented by an empty operation set.
    #[error("magnetic operation set is empty")]
    Empty,
    /// An operation contains a non-finite fractional translation component.
    #[error(
        "operation {operation_index} has non-finite translation component {component}: {value}"
    )]
    NonFiniteTranslation {
        operation_index: usize,
        component: usize,
        value: f64,
    },
    /// A crystallographic rotation must be unimodular.
    #[error(
        "operation {operation_index} has rotation determinant {determinant}, expected +1 or -1"
    )]
    InvalidRotationDeterminant {
        operation_index: usize,
        determinant: i128,
    },
    /// Composing two otherwise unimodular integer rotations would leave the
    /// public `i32` matrix representation.
    #[error(
        "rotation product for operations {left_index} and {right_index} overflows i32 at ({row}, {column})"
    )]
    RotationProductOverflow {
        left_index: usize,
        right_index: usize,
        row: usize,
        column: usize,
    },
    /// Two entries represent the same magnetic coset operation.
    #[error(
        "operations {first_index} and {duplicate_index} are duplicates modulo a lattice translation"
    )]
    Duplicate {
        first_index: usize,
        duplicate_index: usize,
    },
    /// The unprimed identity $\{I|0\}$ is absent.
    #[error("operation set does not contain the unprimed identity")]
    MissingIdentity,
    /// The indexed operation has no two-sided inverse in the set.
    #[error("operation {operation_index} has no two-sided inverse in the set")]
    MissingInverse { operation_index: usize },
    /// A product is missing from the supplied set.
    #[error(
        "operation set is not closed: product of entries {left_index} and {right_index} is absent ({product:?})"
    )]
    NotClosed {
        left_index: usize,
        right_index: usize,
        product: SymmetryOp,
    },
    /// The symmetry precision supplied to database identification is invalid.
    #[error("symmetry precision must be finite and positive, got {symprec}")]
    InvalidSymmetryPrecision { symprec: f64 },
    /// The optional structural-supergroup Hall provenance is out of range.
    #[error("structural-supergroup Hall number {hall_number} is outside 1..=530")]
    InvalidStructuralSupergroupHall { hall_number: usize },
    /// The family-space-group Hall setting could not be obtained from the
    /// time-reversal-erased spatial projection.
    #[error("failed to derive the family-space-group Hall setting: {source}")]
    FamilyHallDerivationFailed {
        #[source]
        source: SymError,
    },
    /// Magnetic database identification failed even with the derived family
    /// Hall setting.  In particular, UNI-setting ambiguity is propagated.
    #[error(
        "failed to identify magnetic group with derived family Hall {family_hall_number}: {source}"
    )]
    MagneticIdentificationFailed {
        family_hall_number: usize,
        #[source]
        source: SymError,
    },
}

/// Complete, setting-aware magnetic-space-group identification.
///
/// `family_hall_number` is derived from the spatial projections of the
/// *surviving Hamiltonian operations*.  `structural_supergroup_hall` is only
/// provenance: it records the group from which candidates were tested, and is
/// never supplied to magnetic identification as a parent-family hint.
#[derive(Debug, Clone)]
pub struct MagneticGroupIdentification {
    /// UNI magnetic-space-group number.
    pub uni_number: usize,
    /// Litvin magnetic-space-group number.
    pub litvin_number: usize,
    /// Crystallographic family-space-group number in `1..=230`.
    pub spacegroup_number: usize,
    /// Magnetic group type (I, II, III, or IV).
    pub magnetic_type: MagneticType,
    /// BNS number, for example `"123.345"`.
    pub bns_number: String,
    /// OG number.
    pub og_number: String,
    /// Hall setting returned by magnetic database matching.
    pub hall_number: usize,
    /// Hall setting independently derived from the effective family group.
    pub family_hall_number: usize,
    /// Hall setting of the original structural supergroup, if supplied.
    pub structural_supergroup_hall: Option<usize>,
    /// Input-to-standard transformation returned by magnetic identification.
    pub transformation_matrix: Mat3,
    /// Origin shift accompanying [`Self::transformation_matrix`].
    pub origin_shift: Vec3,
    /// Cartesian rotation aligning the standardized lattice.
    pub std_rotation_matrix: Mat3,
}

/// An owned magnetic operation set whose group axioms have been checked.
///
/// Operations act on fractional coordinates by
///
/// $$
/// g\mathbf{x}=W_g\mathbf{x}+\mathbf{w}_g,
/// \qquad W_g\in GL(3,\mathbb{Z}),
/// $$
///
/// and carry a time-reversal bit $\epsilon_g\in\{0,1\}$.  Equality is
/// crystallographic equality,
///
/// $$
/// (W,\mathbf{w},\epsilon)\sim(W',\mathbf{w}',\epsilon')
/// \iff W=W',\;\epsilon=\epsilon',\;
/// \mathbf{w}-\mathbf{w}'\in\mathbb{Z}^3,
/// $$
///
/// evaluated componentwise with the construction tolerance.  Composition is
/// delegated to the same Seitz implementation used by the irrep code:
///
/// $$
/// \{W_1|\mathbf{w}_1\}\{W_2|\mathbf{w}_2\}
/// =\{W_1W_2|\mathbf{w}_1+W_1\mathbf{w}_2\},\qquad
/// \epsilon_{12}=\epsilon_1\mathbin{\mathrm{xor}}\epsilon_2.
/// $$
///
/// The constructor verifies non-emptiness, finite translations,
/// $\det W=\pm1$, uniqueness, the unprimed identity, two-sided inverses, and
/// multiplication closure.  Translations are normalized into $[0,1)$ only
/// after all scalar inputs have passed the finiteness check.
#[derive(Debug, Clone)]
pub struct ValidatedMagneticOperationSet {
    operations: SymmetryOps,
    tolerance: f64,
}

impl ValidatedMagneticOperationSet {
    /// Validate and take an owned copy of `operations`.
    ///
    /// `tolerance` is used only for fractional translations modulo integer
    /// lattice vectors; rotations and time-reversal flags are compared
    /// exactly.  Every failure reports a witness through
    /// [`MagneticOperationSetError`].
    pub fn try_from_symmetry_ops(
        operations: &SymmetryOps,
        tolerance: f64,
    ) -> Result<Self, MagneticOperationSetError> {
        if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 0.5 {
            return Err(MagneticOperationSetError::InvalidTolerance { tolerance });
        }
        if operations.is_empty() {
            return Err(MagneticOperationSetError::Empty);
        }

        let mut normalized = Vec::with_capacity(operations.len());
        for (operation_index, operation) in operations.iter().copied().enumerate() {
            for (component, value) in operation.translation.into_iter().enumerate() {
                if !value.is_finite() {
                    return Err(MagneticOperationSetError::NonFiniteTranslation {
                        operation_index,
                        component,
                        value,
                    });
                }
            }
            let determinant = determinant_i128(&operation.rotation);
            if determinant != 1 && determinant != -1 {
                return Err(MagneticOperationSetError::InvalidRotationDeterminant {
                    operation_index,
                    determinant,
                });
            }
            normalized.push(SymmetryOp {
                translation: operation.translation.map(normalize_fractional),
                ..operation
            });
        }

        for duplicate_index in 0..normalized.len() {
            for first_index in 0..duplicate_index {
                if operations_equivalent(
                    &normalized[first_index],
                    &normalized[duplicate_index],
                    tolerance,
                    false,
                ) {
                    return Err(MagneticOperationSetError::Duplicate {
                        first_index,
                        duplicate_index,
                    });
                }
            }
        }

        let identity = SymmetryOp {
            rotation: IDENTITY_ROTATION,
            translation: [0.0; 3],
            time_reversal: false,
        };
        if !normalized
            .iter()
            .any(|operation| operations_equivalent(operation, &identity, tolerance, false))
        {
            return Err(MagneticOperationSetError::MissingIdentity);
        }

        // `compose_seitz` intentionally uses the library-wide `i32` rotation
        // representation.  Prove its integer arithmetic representable before
        // calling it so adversarial, but unimodular, matrices cannot trigger a
        // debug-overflow panic.
        for (left_index, left) in normalized.iter().enumerate() {
            for (right_index, right) in normalized.iter().enumerate() {
                ensure_rotation_product_representable(
                    &left.rotation,
                    &right.rotation,
                    left_index,
                    right_index,
                )?;
            }
        }

        for (operation_index, operation) in normalized.iter().enumerate() {
            let has_inverse = normalized.iter().any(|candidate| {
                let left = compose(operation, candidate);
                let right = compose(candidate, operation);
                operations_equivalent(&left, &identity, tolerance, false)
                    && operations_equivalent(&right, &identity, tolerance, false)
            });
            if !has_inverse {
                return Err(MagneticOperationSetError::MissingInverse { operation_index });
            }
        }

        for (left_index, left) in normalized.iter().enumerate() {
            for (right_index, right) in normalized.iter().enumerate() {
                let product = compose(left, right);
                if !normalized
                    .iter()
                    .any(|candidate| operations_equivalent(candidate, &product, tolerance, false))
                {
                    return Err(MagneticOperationSetError::NotClosed {
                        left_index,
                        right_index,
                        product,
                    });
                }
            }
        }

        Ok(Self {
            operations: SymmetryOps {
                operations: normalized,
            },
            tolerance,
        })
    }

    /// Borrow the normalized, validated operations.
    pub fn operations(&self) -> &SymmetryOps {
        &self.operations
    }

    /// Number of magnetic operations.
    pub fn len(&self) -> usize {
        self.operations.len()
    }

    /// A validated group is never empty; this method exists for collection
    /// API symmetry and always returns `false`.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Translation-equivalence tolerance used during validation.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Identify the effective magnetic space group in a setting-safe way.
    ///
    /// First, time reversal is erased and equal spatial projections are
    /// deduplicated:
    ///
    /// $$
    /// F(M)=\left\{\{W_g|\mathbf{w}_g\}:g\in M\right\}.
    /// $$
    ///
    /// Because the input $M$ has already passed closure validation, $F(M)$ is
    /// its family space group.  Its Hall setting is derived using the supplied
    /// lattice and then passed to
    /// [`crate::magnetic_spacegroup::msg_identify_with_parent_hall`].  This is
    /// crucial for Type-IV settings: the Hall number of a *larger structural
    /// group* is not generally the Hall number of $F(M)$ and therefore cannot
    /// safely disambiguate a UNI number.
    ///
    /// `structural_supergroup_hall` is retained only in the returned
    /// [`MagneticGroupIdentification`] as provenance.  It has no effect on the
    /// identification calculation.
    pub fn identify(
        &self,
        lattice: &Mat3,
        structural_supergroup_hall: Option<usize>,
        symprec: f64,
    ) -> Result<MagneticGroupIdentification, MagneticOperationSetError> {
        if !symprec.is_finite() || symprec <= 0.0 {
            return Err(MagneticOperationSetError::InvalidSymmetryPrecision { symprec });
        }
        if let Some(hall_number) = structural_supergroup_hall
            && !(1..=530).contains(&hall_number)
        {
            return Err(MagneticOperationSetError::InvalidStructuralSupergroupHall { hall_number });
        }

        let family = self.family_spatial_projection();
        let rotations: Vec<_> = family.iter().map(|operation| operation.rotation).collect();
        let translations: Vec<_> = family
            .iter()
            .map(|operation| operation.translation)
            .collect();
        let family_hall_number =
            crate::identify_hall_number(&rotations, &translations, lattice, true, symprec)
                .map_err(
                    |source| MagneticOperationSetError::FamilyHallDerivationFailed { source },
                )?;

        let magnetic_symmetry = self.as_magnetic_symmetry();
        let dataset = crate::magnetic_spacegroup::msg_identify_with_parent_hall(
            lattice,
            &magnetic_symmetry,
            Some(family_hall_number),
            symprec,
        )
        .map_err(
            |source| MagneticOperationSetError::MagneticIdentificationFailed {
                family_hall_number,
                source,
            },
        )?;
        let metadata = crate::msg_database::msgdb_get_magnetic_spacegroup_type(dataset.uni_number);

        Ok(MagneticGroupIdentification {
            uni_number: dataset.uni_number,
            litvin_number: metadata.litvin_number,
            spacegroup_number: metadata.number,
            magnetic_type: dataset.msg_type,
            bns_number: metadata.bns_number.trim().to_string(),
            og_number: metadata.og_number.trim().to_string(),
            hall_number: dataset.hall_number,
            family_hall_number,
            structural_supergroup_hall,
            transformation_matrix: dataset.transformation_matrix,
            origin_shift: dataset.origin_shift,
            std_rotation_matrix: dataset.std_rotation_matrix,
        })
    }

    fn family_spatial_projection(&self) -> SymmetryOps {
        let mut operations: Vec<SymmetryOp> = Vec::with_capacity(self.len());
        for operation in self.operations.iter() {
            let projection = SymmetryOp {
                time_reversal: false,
                ..*operation
            };
            if !operations
                .iter()
                .any(|existing| operations_equivalent(existing, &projection, self.tolerance, true))
            {
                operations.push(projection);
            }
        }
        SymmetryOps { operations }
    }

    fn as_magnetic_symmetry(&self) -> MagneticSymmetry {
        let mut symmetry = MagneticSymmetry::new(self.len());
        for (index, operation) in self.operations.iter().enumerate() {
            symmetry.rot[index] = operation.rotation;
            symmetry.trans[index] = operation.translation;
            symmetry.timerev[index] = operation.time_reversal;
        }
        symmetry
    }
}

fn normalize_fractional(value: f64) -> f64 {
    let normalized = value.rem_euclid(1.0);
    if normalized == 1.0 || normalized == -0.0 {
        0.0
    } else {
        normalized
    }
}

fn determinant_i128(matrix: &Mat3I) -> i128 {
    let value = |row: usize, column: usize| i128::from(matrix[row][column]);
    value(0, 0) * (value(1, 1) * value(2, 2) - value(1, 2) * value(2, 1))
        + value(0, 1) * (value(1, 2) * value(2, 0) - value(1, 0) * value(2, 2))
        + value(0, 2) * (value(1, 0) * value(2, 1) - value(1, 1) * value(2, 0))
}

fn ensure_rotation_product_representable(
    left: &Mat3I,
    right: &Mat3I,
    left_index: usize,
    right_index: usize,
) -> Result<(), MagneticOperationSetError> {
    for (row, left_row) in left.iter().enumerate() {
        for (column, _) in right[0].iter().enumerate() {
            let products = (0..3).map(|inner| {
                i64::from(left_row[inner]).checked_mul(i64::from(right[inner][column]))
            });
            let mut sum = 0_i64;
            for product in products {
                let Some(product) = product else {
                    return Err(MagneticOperationSetError::RotationProductOverflow {
                        left_index,
                        right_index,
                        row,
                        column,
                    });
                };
                let Some(next_sum) = sum.checked_add(product) else {
                    return Err(MagneticOperationSetError::RotationProductOverflow {
                        left_index,
                        right_index,
                        row,
                        column,
                    });
                };
                sum = next_sum;
            }
            if i32::try_from(sum).is_err()
                || (0..3).any(|inner| left_row[inner].checked_mul(right[inner][column]).is_none())
            {
                return Err(MagneticOperationSetError::RotationProductOverflow {
                    left_index,
                    right_index,
                    row,
                    column,
                });
            }

            // Match the exact left-associative evaluation in
            // `mat_multiply_matrix_i3`, including intermediate additions.
            let first = left_row[0] * right[0][column];
            let second = left_row[1] * right[1][column];
            let third = left_row[2] * right[2][column];
            if first
                .checked_add(second)
                .and_then(|partial| partial.checked_add(third))
                .is_none()
            {
                return Err(MagneticOperationSetError::RotationProductOverflow {
                    left_index,
                    right_index,
                    row,
                    column,
                });
            }
        }
    }
    Ok(())
}

fn operations_equivalent(
    left: &SymmetryOp,
    right: &SymmetryOp,
    tolerance: f64,
    ignore_time_reversal: bool,
) -> bool {
    left.rotation == right.rotation
        && (ignore_time_reversal || left.time_reversal == right.time_reversal)
        && left
            .translation
            .into_iter()
            .zip(right.translation)
            .all(|(left, right)| {
                let difference = left - right;
                (difference - difference.round()).abs() <= tolerance
            })
}

fn compose(left: &SymmetryOp, right: &SymmetryOp) -> SymmetryOp {
    let left = SeitzOp::new(left.rotation, left.translation, left.time_reversal);
    let right = SeitzOp::new(right.rotation, right.translation, right.time_reversal);
    let (product, _) = compose_seitz(&left, &right);
    SymmetryOp {
        rotation: product.rot,
        translation: product.trans,
        time_reversal: product.timerev,
    }
}

/// Failure to lift an orthogonal spatial transformation to spin $1/2$.
#[derive(Debug, Clone, thiserror::Error)]
pub enum SpinLiftError {
    /// The matrix-validation tolerance must be finite and positive.
    #[error("spin-lift tolerance must be finite and positive, got {tolerance}")]
    InvalidTolerance { tolerance: f64 },
    /// One matrix entry is NaN or infinite.
    #[error("matrix entry ({row}, {column}) is non-finite: {value}")]
    NonFiniteMatrixEntry {
        row: usize,
        column: usize,
        value: f64,
    },
    /// $Q^TQ$ differs from the identity beyond tolerance.
    #[error("matrix is not orthogonal; maximum |Q^T Q - I| is {max_residual}")]
    NotOrthogonal { max_residual: f64 },
    /// An orthogonal matrix must have determinant $+1$ or $-1$.
    #[error("matrix determinant {determinant} is not within tolerance of +1 or -1")]
    InvalidDeterminant { determinant: f64 },
    /// The quaternion extraction became numerically singular.
    #[error("failed to extract a finite normalized spin quaternion")]
    NumericalFailure,
}

/// Lift a Cartesian polar transformation to its axial spin-$1/2$ action.
///
/// The input is an orthogonal Cartesian matrix $Q\in O(3)$ acting on polar
/// vectors.  Spin is axial, so its proper rotation is
///
/// $$
/// R_{\mathrm{spin}}=(\det Q)Q\in SO(3).
/// $$
///
/// The returned coefficients are the normalized quaternion
/// $(q_0,q_x,q_y,q_z)$ in the Pauli convention
///
/// $$
/// U(Q)=q_0 I-i(q_x\sigma_x+q_y\sigma_y+q_z\sigma_z),
/// \qquad q_0^2+q_x^2+q_y^2+q_z^2=1.
/// $$
///
/// Thus a proper rotation by angle $\theta$ about the unit vector
/// $\hat{\mathbf n}$ gives
/// $q=(\cos\frac{\theta}{2},\hat{\mathbf n}\sin\frac{\theta}{2})$.
/// The double-cover sign is canonicalized by choosing the first non-negligible
/// component of $(q_0,q_x,q_y,q_z)$ positive.  For inversion $Q=-I$, axial
/// spin is unchanged and the result is `[1, 0, 0, 0]`.
///
/// # Errors
///
/// Every entry must be finite, $Q^TQ$ must agree with $I$, and $\det Q$ must
/// agree with $\pm1$, all to `tolerance`.
pub fn axial_spin_half_lift(matrix: &Mat3, tolerance: f64) -> Result<[f64; 4], SpinLiftError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(SpinLiftError::InvalidTolerance { tolerance });
    }
    for (row, values) in matrix.iter().enumerate() {
        for (column, value) in values.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(SpinLiftError::NonFiniteMatrixEntry { row, column, value });
            }
        }
    }

    let mut max_residual = 0.0_f64;
    for row in 0..3 {
        for column in 0..3 {
            let dot = (0..3)
                .map(|index| matrix[index][row] * matrix[index][column])
                .sum::<f64>();
            let expected = f64::from(row == column);
            max_residual = max_residual.max((dot - expected).abs());
        }
    }
    if max_residual > tolerance {
        return Err(SpinLiftError::NotOrthogonal { max_residual });
    }

    let determinant = mat_get_determinant_d3(matrix);
    let determinant_sign = if (determinant - 1.0).abs() <= tolerance {
        1.0
    } else if (determinant + 1.0).abs() <= tolerance {
        -1.0
    } else {
        return Err(SpinLiftError::InvalidDeterminant { determinant });
    };
    let rotation = matrix.map(|row| row.map(|value| determinant_sign * value));
    let trace = rotation[0][0] + rotation[1][1] + rotation[2][2];

    let mut quaternion = if trace > 0.0 {
        let scale = 2.0 * (trace + 1.0).max(0.0).sqrt();
        if scale <= f64::EPSILON {
            return Err(SpinLiftError::NumericalFailure);
        }
        [
            0.25 * scale,
            (rotation[2][1] - rotation[1][2]) / scale,
            (rotation[0][2] - rotation[2][0]) / scale,
            (rotation[1][0] - rotation[0][1]) / scale,
        ]
    } else if rotation[0][0] > rotation[1][1] && rotation[0][0] > rotation[2][2] {
        let scale = 2.0
            * (1.0 + rotation[0][0] - rotation[1][1] - rotation[2][2])
                .max(0.0)
                .sqrt();
        if scale <= f64::EPSILON {
            return Err(SpinLiftError::NumericalFailure);
        }
        [
            (rotation[2][1] - rotation[1][2]) / scale,
            0.25 * scale,
            (rotation[0][1] + rotation[1][0]) / scale,
            (rotation[0][2] + rotation[2][0]) / scale,
        ]
    } else if rotation[1][1] > rotation[2][2] {
        let scale = 2.0
            * (1.0 + rotation[1][1] - rotation[0][0] - rotation[2][2])
                .max(0.0)
                .sqrt();
        if scale <= f64::EPSILON {
            return Err(SpinLiftError::NumericalFailure);
        }
        [
            (rotation[0][2] - rotation[2][0]) / scale,
            (rotation[0][1] + rotation[1][0]) / scale,
            0.25 * scale,
            (rotation[1][2] + rotation[2][1]) / scale,
        ]
    } else {
        let scale = 2.0
            * (1.0 + rotation[2][2] - rotation[0][0] - rotation[1][1])
                .max(0.0)
                .sqrt();
        if scale <= f64::EPSILON {
            return Err(SpinLiftError::NumericalFailure);
        }
        [
            (rotation[1][0] - rotation[0][1]) / scale,
            (rotation[0][2] + rotation[2][0]) / scale,
            (rotation[1][2] + rotation[2][1]) / scale,
            0.25 * scale,
        ]
    };

    let norm = quaternion
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    if !norm.is_finite() || norm <= f64::EPSILON {
        return Err(SpinLiftError::NumericalFailure);
    }
    for value in &mut quaternion {
        *value /= norm;
    }
    if quaternion
        .iter()
        .copied()
        .find(|component| component.abs() > tolerance)
        .is_some_and(|component| component < 0.0)
    {
        for value in &mut quaternion {
            *value = -*value;
        }
    }
    Ok(quaternion)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity(time_reversal: bool) -> SymmetryOp {
        SymmetryOp {
            rotation: IDENTITY_ROTATION,
            translation: [0.0; 3],
            time_reversal,
        }
    }

    fn ordinary_c2() -> SymmetryOps {
        SymmetryOps {
            operations: vec![
                identity(false),
                SymmetryOp {
                    rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    translation: [0.0; 3],
                    time_reversal: false,
                },
            ],
        }
    }

    fn assert_quaternion_close(actual: [f64; 4], expected: [f64; 4]) {
        for (actual, expected) in actual.into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
        }
    }

    #[test]
    fn validates_closed_ordinary_and_grey_groups() {
        let ordinary = ordinary_c2();
        let validated =
            ValidatedMagneticOperationSet::try_from_symmetry_ops(&ordinary, 1e-8).unwrap();
        assert_eq!(validated.len(), 2);

        let grey = ordinary.grey_extension().unwrap();
        let validated = ValidatedMagneticOperationSet::try_from_symmetry_ops(&grey, 1e-8).unwrap();
        assert_eq!(validated.len(), 4);
    }

    #[test]
    fn reports_a_nonclosure_witness() {
        let c4 = SymmetryOp {
            rotation: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            translation: [0.0; 3],
            time_reversal: false,
        };
        let c4_inverse = SymmetryOp {
            rotation: [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
            translation: [0.0; 3],
            time_reversal: false,
        };
        let operations = SymmetryOps {
            operations: vec![identity(false), c4, c4_inverse],
        };

        let error =
            ValidatedMagneticOperationSet::try_from_symmetry_ops(&operations, 1e-8).unwrap_err();
        assert!(matches!(
            error,
            MagneticOperationSetError::NotClosed {
                left_index: 1,
                right_index: 1,
                ..
            }
        ));
    }

    #[test]
    fn reports_large_rotation_products_without_panicking() {
        let shear = SymmetryOp {
            rotation: [[1, i32::MAX, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0.0; 3],
            time_reversal: false,
        };
        let inverse_shear = SymmetryOp {
            rotation: [[1, -i32::MAX, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0.0; 3],
            time_reversal: false,
        };
        let operations = SymmetryOps {
            operations: vec![identity(false), shear, inverse_shear],
        };

        assert!(matches!(
            ValidatedMagneticOperationSet::try_from_symmetry_ops(&operations, 1e-8),
            Err(MagneticOperationSetError::RotationProductOverflow {
                left_index: 1,
                right_index: 1,
                row: 0,
                column: 1,
            })
        ));
    }

    #[test]
    fn rejects_true_duplicates_but_distinguishes_time_reversal() {
        let mut translated_identity = identity(false);
        translated_identity.translation = [1.0, -2.0, 3.0];
        let duplicate = SymmetryOps {
            operations: vec![identity(false), translated_identity],
        };
        assert!(matches!(
            ValidatedMagneticOperationSet::try_from_symmetry_ops(&duplicate, 1e-8),
            Err(MagneticOperationSetError::Duplicate {
                first_index: 0,
                duplicate_index: 1
            })
        ));

        let grey_identity = SymmetryOps {
            operations: vec![identity(false), identity(true)],
        };
        ValidatedMagneticOperationSet::try_from_symmetry_ops(&grey_identity, 1e-8).unwrap();
    }

    #[test]
    fn identifies_a_known_group_lowered_from_cubic_structural_provenance() {
        let uni_number = 1005;
        let hall_number = crate::api::find_first_hall_for_uni(uni_number).unwrap();
        let operations =
            SymmetryOps::from_magnetic_database_with_hall(uni_number, hall_number).unwrap();
        let validated =
            ValidatedMagneticOperationSet::try_from_symmetry_ops(&operations, 1e-7).unwrap();
        let cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        let identification = validated.identify(&cubic_lattice, Some(517), 1e-5).unwrap();

        assert_eq!(identification.uni_number, uni_number);
        assert_eq!(identification.structural_supergroup_hall, Some(517));
        assert_eq!(identification.family_hall_number, hall_number);
        assert_eq!(identification.bns_number, "123.345");
    }

    #[test]
    fn axial_spin_lift_handles_identity_c2_and_inversion() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c2z = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]];
        let inversion = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];

        assert_quaternion_close(
            axial_spin_half_lift(&identity, 1e-12).unwrap(),
            [1.0, 0.0, 0.0, 0.0],
        );
        assert_quaternion_close(
            axial_spin_half_lift(&c2z, 1e-12).unwrap(),
            [0.0, 0.0, 0.0, 1.0],
        );
        assert_quaternion_close(
            axial_spin_half_lift(&inversion, 1e-12).unwrap(),
            [1.0, 0.0, 0.0, 0.0],
        );
    }

    #[test]
    fn axial_spin_lift_rejects_a_nonorthogonal_matrix() {
        let shear = [[1.0, 0.2, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(matches!(
            axial_spin_half_lift(&shear, 1e-8),
            Err(SpinLiftError::NotOrthogonal { .. })
        ));
    }
}
