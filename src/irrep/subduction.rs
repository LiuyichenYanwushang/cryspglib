//! Exact rational affine and lattice layer for full irrep subduction.
//!
//! Staging area for tasks 3-12 of `docs/full-irrep-subduction-plan.md`; the
//! conventions implemented here are pinned in `docs/subduction-conventions.md`.
//! The module is deliberately `#[doc(hidden)]` until task 12 gives it a
//! supported public surface.
//!
//! Why not reuse the existing algebra:
//!
//! * `wigner::ExactSeitzOp` fixes translations to a **twelfth** grid, while the
//!   isotropy origins reach denominator 16 and affine images reach 48, so a
//!   rounded grid cannot represent the embeddings this engine must verify.
//! * `SymmetryOps` stores `f64` translations.  Every conversion here goes
//!   through a checked grid snap ([`Rat::from_grid`]) instead of a cast, and
//!   any value that is not exactly on the grid is an error.
//! * lattice membership is defined **per lattice**: an operation that is equal
//!   modulo the parent lattice need not be equal modulo the subgroup lattice,
//!   and merging those is a real bug (see [`Lattice::deduplicate`]).
//!
//! The layer is exact and total: rational arithmetic is checked, matrix
//! inversion reports singular input, and rotations are only converted to
//! integers after the integrality check.

use crate::SymError;
use crate::api::SymmetryOps;
use crate::irrep::types::generated_data::SG_DATA_HALL;
use crate::mathfunc::Mat3I;

/// Denominator of the operation tables shipped in `data_space.txt`/Hall data.
pub const SOURCE_TRANSLATION_GRID: i128 = 12;

/// Tolerance for the grid snap, in units of the grid denominator.  The tables
/// are exact; this only absorbs `f64` round-off from the stored values (the
/// algebra itself never uses floats).
const GRID_SNAP_TOLERANCE: f64 = 1e-9;

// ── Errors ───────────────────────────────────────────────────────────────────

/// Errors from the exact subduction layer.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SubductionError {
    /// Space group number outside 1-230.
    #[error("space group {sg} is outside 1-230")]
    InvalidSpaceGroup { sg: u8 },
    /// `SG_DATA_HALL` has no Hall number for this space group.
    ///
    /// The strict source never falls back to the first Hall setting: the
    /// isotropy tables are indexed in the data Hall frame, so a silent fallback
    /// would produce operations in a different setting.
    #[error("space group {sg} has no data Hall setting; no fallback is attempted")]
    StrictHallUnavailable { sg: u8 },
    /// The recorded Hall setting could not be loaded.
    #[error("space group {sg}: loading Hall {hall} failed: {source}")]
    HallLoadFailed {
        sg: u8,
        hall: usize,
        source: SymError,
    },
    /// A rational had a zero denominator.
    #[error("rational with zero denominator")]
    ZeroDenominator,
    /// Division by a rational zero.
    #[error("division by zero rational")]
    DivisionByZero,
    /// A rational numerator or denominator overflowed `i128`.
    #[error("rational overflow during {operation}")]
    RationalOverflow { operation: &'static str },
    /// A non-finite value was offered to a conversion.
    #[error("non-finite value {value}")]
    NonFiniteValue { value: f64 },
    /// The grid denominator was not positive (or did not fit).
    #[error("invalid grid denominator {denominator}")]
    InvalidGrid { denominator: i128 },
    /// A stored translation is not exactly on the assumed grid.
    #[error("value {value} is not on the 1/{denominator} grid")]
    OffGridValue { value: f64, denominator: i128 },
    /// A rational was required to be an integer and was not.
    #[error("value {value} is not an integer")]
    NotAnInteger { value: Rat },
    /// An integer did not fit the target type.
    #[error("value {value} does not fit the target integer type")]
    IntegerOutOfRange { value: i128 },
    /// A matrix that had to be invertible was singular.
    #[error("singular matrix")]
    SingularMatrix,
    /// `T R T^-1` was not an integer matrix, so the affine map is not a
    /// symmetry embedding of the requested operation set.
    ///
    /// Boxed: an exact 3x3 rational matrix is 288 bytes, and an error that big
    /// would be paid on every `Result` in this module.
    #[error("image of the rotation is not an integer matrix")]
    NonIntegralRotationImage { rotation: Box<Mat3R> },
}

// ── Checked rationals ────────────────────────────────────────────────────────

/// Exact rational number, always normalized (`den > 0`, `gcd(num, den) = 1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Rat {
    num: i128,
    den: i128,
}

impl Rat {
    /// `0`.
    pub const ZERO: Self = Self { num: 0, den: 1 };
    /// `1`.
    pub const ONE: Self = Self { num: 1, den: 1 };

    /// Build `num / den`, normalizing the sign and the common divisor.
    pub fn new(num: i128, den: i128) -> Result<Self, SubductionError> {
        if den == 0 {
            return Err(SubductionError::ZeroDenominator);
        }
        let (num, den) = if den < 0 {
            (
                num.checked_neg()
                    .ok_or(SubductionError::RationalOverflow {
                        operation: "normalize",
                    })?,
                den.checked_neg()
                    .ok_or(SubductionError::RationalOverflow {
                        operation: "normalize",
                    })?,
            )
        } else {
            (num, den)
        };
        let divisor = gcd_positive(num.unsigned_abs(), den.unsigned_abs());
        let divisor = i128::try_from(divisor).map_err(|_| SubductionError::RationalOverflow {
            operation: "normalize",
        })?;
        Ok(Self {
            num: num / divisor,
            den: den / divisor,
        })
    }

    /// Exact integer as a rational.
    pub const fn from_integer(value: i128) -> Self {
        Self { num: value, den: 1 }
    }

    /// Snap a floating value to the `1 / denominator` grid.
    ///
    /// This is the only float entry point: the value must be finite and
    /// integral *on the grid* within [`GRID_SNAP_TOLERANCE`], otherwise the
    /// conversion reports [`SubductionError::OffGridValue`].  There is no
    /// truncation and no implicit denominator.
    pub fn from_grid(value: f64, denominator: i128) -> Result<Self, SubductionError> {
        if denominator <= 0 {
            return Err(SubductionError::InvalidGrid { denominator });
        }
        if !value.is_finite() {
            return Err(SubductionError::NonFiniteValue { value });
        }
        let scale = f64::from(
            i32::try_from(denominator)
                .map_err(|_| SubductionError::InvalidGrid { denominator })?,
        );
        let scaled = value * scale;
        let rounded = scaled.round();
        if (scaled - rounded).abs() > GRID_SNAP_TOLERANCE {
            return Err(SubductionError::OffGridValue { value, denominator });
        }
        let numerator =
            exact_i128(rounded).ok_or(SubductionError::OffGridValue { value, denominator })?;
        Self::new(numerator, denominator)
    }

    /// Numerator (sign-carrying, normalized).
    pub const fn numerator(self) -> i128 {
        self.num
    }

    /// Denominator (always positive).
    pub const fn denominator(self) -> i128 {
        self.den
    }

    /// Whether the value is exactly zero.
    pub const fn is_zero(self) -> bool {
        self.num == 0
    }

    /// Whether the value is an integer.
    pub const fn is_integer(self) -> bool {
        self.den == 1
    }

    /// Exact integer value, or [`SubductionError::NotAnInteger`].
    pub fn to_integer(self) -> Result<i128, SubductionError> {
        if self.den == 1 {
            Ok(self.num)
        } else {
            Err(SubductionError::NotAnInteger { value: self })
        }
    }

    /// Exact integer value that must also fit `i32`.
    pub fn to_i32(self) -> Result<i32, SubductionError> {
        let value = self.to_integer()?;
        i32::try_from(value).map_err(|_| SubductionError::IntegerOutOfRange { value })
    }

    /// Checked addition.
    pub fn checked_add(self, other: Self) -> Result<Self, SubductionError> {
        let shared = gcd_positive(self.den.unsigned_abs(), other.den.unsigned_abs());
        let shared = i128::try_from(shared).map_err(|_| SubductionError::RationalOverflow {
            operation: "add",
        })?;
        let num = self
            .num
            .checked_mul(other.den / shared)
            .and_then(|left| {
                other
                    .num
                    .checked_mul(self.den / shared)
                    .and_then(|right| left.checked_add(right))
            })
            .ok_or(SubductionError::RationalOverflow { operation: "add" })?;
        let den = (self.den / shared)
            .checked_mul(other.den)
            .ok_or(SubductionError::RationalOverflow { operation: "add" })?;
        Self::new(num, den)
    }

    /// Checked subtraction.
    pub fn checked_sub(self, other: Self) -> Result<Self, SubductionError> {
        self.checked_add(other.checked_neg()?)
    }

    /// Checked negation.
    pub fn checked_neg(self) -> Result<Self, SubductionError> {
        Ok(Self {
            num: self
                .num
                .checked_neg()
                .ok_or(SubductionError::RationalOverflow {
                    operation: "negate",
                })?,
            den: self.den,
        })
    }

    /// Checked multiplication, cross-reduced before multiplying so that the
    /// intermediate products stay as small as the result allows.
    pub fn checked_mul(self, other: Self) -> Result<Self, SubductionError> {
        let left = gcd_positive(self.num.unsigned_abs(), other.den.unsigned_abs());
        let right = gcd_positive(other.num.unsigned_abs(), self.den.unsigned_abs());
        let left = i128::try_from(left)
            .map_err(|_| SubductionError::RationalOverflow { operation: "multiply" })?;
        let right = i128::try_from(right)
            .map_err(|_| SubductionError::RationalOverflow { operation: "multiply" })?;
        let num = (self.num / left)
            .checked_mul(other.num / right)
            .ok_or(SubductionError::RationalOverflow {
                operation: "multiply",
            })?;
        let den = (self.den / right)
            .checked_mul(other.den / left)
            .ok_or(SubductionError::RationalOverflow {
                operation: "multiply",
            })?;
        Self::new(num, den)
    }

    /// Checked division.
    pub fn checked_div(self, other: Self) -> Result<Self, SubductionError> {
        if other.num == 0 {
            return Err(SubductionError::DivisionByZero);
        }
        let num = self
            .num
            .checked_mul(other.den)
            .ok_or(SubductionError::RationalOverflow { operation: "divide" })?;
        let den = self
            .den
            .checked_mul(other.num)
            .ok_or(SubductionError::RationalOverflow { operation: "divide" })?;
        Self::new(num, den)
    }
}

impl std::fmt::Display for Rat {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.den == 1 {
            write!(formatter, "{}", self.num)
        } else {
            write!(formatter, "{}/{}", self.num, self.den)
        }
    }
}

impl std::fmt::Display for Vec3R {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "({},{},{})", self.0[0], self.0[1], self.0[2])
    }
}

impl std::fmt::Display for Mat3R {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "[{},{},{};{},{},{};{},{},{}]",
            self.0[0][0],
            self.0[0][1],
            self.0[0][2],
            self.0[1][0],
            self.0[1][1],
            self.0[1][2],
            self.0[2][0],
            self.0[2][1],
            self.0[2][2]
        )
    }
}

fn gcd_positive(left: u128, right: u128) -> u128 {
    let (mut a, mut b) = (left, right);
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
}

/// Convert an integral `f64` within the exactly representable range.
///
/// Floating point is confined to the grid snap; the value is verified to be
/// integral and small enough that the conversion cannot round or truncate.
fn exact_i128(value: f64) -> Option<i128> {
    /// 2^53: the largest magnitude with unit spacing in `f64`.
    const EXACT_LIMIT: f64 = 9_007_199_254_740_992.0;
    if !value.is_finite() || value.fract() != 0.0 || value.abs() > EXACT_LIMIT {
        return None;
    }
    // Exact by the checks above: integral and inside the f64-integer range.
    Some(i128::from(value as i64))
}

// ── Rational vectors and matrices ────────────────────────────────────────────

/// Exact three-vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Vec3R([Rat; 3]);

impl Vec3R {
    /// Componentwise constructor.
    pub const fn new(values: [Rat; 3]) -> Self {
        Self(values)
    }

    /// The zero vector.
    pub const fn zero() -> Self {
        Self([Rat::ZERO; 3])
    }

    /// Exact integer vector.
    pub fn from_ints(values: [i32; 3]) -> Self {
        Self(values.map(|value| Rat::from_integer(i128::from(value))))
    }

    /// Components.
    pub const fn as_array(&self) -> &[Rat; 3] {
        &self.0
    }

    /// One component.
    pub fn get(&self, index: usize) -> Rat {
        self.0[index]
    }

    /// Whether every component is zero.
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|value| value.is_zero())
    }

    /// Whether every component is an integer.
    pub fn is_integral(&self) -> bool {
        self.0.iter().all(|value| value.is_integer())
    }

    /// Exact componentwise sum.
    pub fn checked_add(&self, other: &Self) -> Result<Self, SubductionError> {
        Ok(Self([
            self.0[0].checked_add(other.0[0])?,
            self.0[1].checked_add(other.0[1])?,
            self.0[2].checked_add(other.0[2])?,
        ]))
    }

    /// Exact componentwise difference.
    pub fn checked_sub(&self, other: &Self) -> Result<Self, SubductionError> {
        Ok(Self([
            self.0[0].checked_sub(other.0[0])?,
            self.0[1].checked_sub(other.0[1])?,
            self.0[2].checked_sub(other.0[2])?,
        ]))
    }

    /// Exact negation.
    pub fn checked_neg(&self) -> Result<Self, SubductionError> {
        Ok(Self([
            self.0[0].checked_neg()?,
            self.0[1].checked_neg()?,
            self.0[2].checked_neg()?,
        ]))
    }

    /// Exact integer components, if every component is an `i32`.
    pub fn to_ints(&self) -> Result<[i32; 3], SubductionError> {
        Ok([self.0[0].to_i32()?, self.0[1].to_i32()?, self.0[2].to_i32()?])
    }
}

/// Exact 3x3 matrix, rows first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mat3R([[Rat; 3]; 3]);

impl Mat3R {
    /// Componentwise constructor.
    pub const fn new(rows: [[Rat; 3]; 3]) -> Self {
        Self(rows)
    }

    /// The identity matrix.
    pub const fn identity() -> Self {
        Self([
            [Rat::ONE, Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, Rat::ONE, Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ])
    }

    /// Exact integer matrix.
    pub fn from_ints(rows: Mat3I) -> Self {
        Self(rows.map(|row| row.map(|value| Rat::from_integer(i128::from(value)))))
    }

    /// Diagonal matrix from exact integers.
    pub fn diagonal(values: [i32; 3]) -> Self {
        Self::from_ints([
            [values[0], 0, 0],
            [0, values[1], 0],
            [0, 0, values[2]],
        ])
    }

    /// Rows.
    pub const fn rows(&self) -> &[[Rat; 3]; 3] {
        &self.0
    }

    /// One row.
    pub const fn row(&self, index: usize) -> &[Rat; 3] {
        &self.0[index]
    }

    /// Whether every entry is an integer.
    pub fn is_integral(&self) -> bool {
        self.0
            .iter()
            .all(|row| row.iter().all(|value| value.is_integer()))
    }

    /// Integer entries, or [`SubductionError::NotAnInteger`].
    pub fn to_int_matrix(&self) -> Result<Mat3I, SubductionError> {
        let mut out = [[0i32; 3]; 3];
        for (row_index, row) in self.0.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                out[row_index][column_index] = value.to_i32()?;
            }
        }
        Ok(out)
    }

    /// Exact matrix product.
    pub fn checked_mul(&self, other: &Self) -> Result<Self, SubductionError> {
        let mut out = [[Rat::ZERO; 3]; 3];
        for (row_index, out_row) in out.iter_mut().enumerate() {
            for (column, cell) in out_row.iter_mut().enumerate() {
                let mut sum = Rat::ZERO;
                for inner in 0..3 {
                    sum = sum.checked_add(
                        self.0[row_index][inner].checked_mul(other.0[inner][column])?,
                    )?;
                }
                *cell = sum;
            }
        }
        Ok(Self(out))
    }

    /// Exact matrix-vector product.
    pub fn checked_mul_vector(&self, vector: &Vec3R) -> Result<Vec3R, SubductionError> {
        let mut out = [Rat::ZERO; 3];
        for (row, cell) in out.iter_mut().enumerate() {
            let mut sum = Rat::ZERO;
            for (column, value) in self.0[row].iter().enumerate() {
                sum = sum.checked_add(value.checked_mul(vector.0[column])?)?;
            }
            *cell = sum;
        }
        Ok(Vec3R(out))
    }

    /// Transpose.
    pub fn transpose(&self) -> Self {
        let mut out = [[Rat::ZERO; 3]; 3];
        for (row, values) in self.0.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                out[column][row] = *value;
            }
        }
        Self(out)
    }

    /// Exact determinant.
    pub fn determinant(&self) -> Result<Rat, SubductionError> {
        let minor = |row: usize, column: usize| -> Result<Rat, SubductionError> {
            let rows: Vec<usize> = (0..3).filter(|index| *index != row).collect();
            let columns: Vec<usize> = (0..3).filter(|index| *index != column).collect();
            let first = self.0[rows[0]][columns[0]].checked_mul(self.0[rows[1]][columns[1]])?;
            let second = self.0[rows[0]][columns[1]].checked_mul(self.0[rows[1]][columns[0]])?;
            first.checked_sub(second)
        };
        let first = self.0[0][0].checked_mul(minor(0, 0)?)?;
        let second = self.0[0][1].checked_mul(minor(0, 1)?)?;
        let third = self.0[0][2].checked_mul(minor(0, 2)?)?;
        first.checked_sub(second)?.checked_add(third)
    }

    /// Exact inverse, or [`SubductionError::SingularMatrix`].
    pub fn inverse(&self) -> Result<Self, SubductionError> {
        let determinant = self.determinant()?;
        if determinant.is_zero() {
            return Err(SubductionError::SingularMatrix);
        }
        let mut out = [[Rat::ZERO; 3]; 3];
        for (row, out_row) in out.iter_mut().enumerate() {
            for (column, cell) in out_row.iter_mut().enumerate() {
                // Cofactor (row, column) of the transpose, i.e. the adjugate.
                let rows: Vec<usize> = (0..3).filter(|index| *index != column).collect();
                let columns: Vec<usize> = (0..3).filter(|index| *index != row).collect();
                let first =
                    self.0[rows[0]][columns[0]].checked_mul(self.0[rows[1]][columns[1]])?;
                let second =
                    self.0[rows[0]][columns[1]].checked_mul(self.0[rows[1]][columns[0]])?;
                let cofactor = first.checked_sub(second)?;
                let sign = if (row + column) % 2 == 0 { 1 } else { -1 };
                let cofactor = if sign < 0 {
                    cofactor.checked_neg()?
                } else {
                    cofactor
                };
                *cell = cofactor.checked_div(determinant)?;
            }
        }
        Ok(Self(out))
    }
}

// ── Affine transform and exact Seitz operations ──────────────────────────────

/// Exact Seitz operation `{R|t}` with an integer rotation.
///
/// `t` is a *representative*; equality of operations only makes sense modulo a
/// lattice, which is why [`ExactSeitz::reduce`] takes the lattice explicitly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExactSeitz {
    rotation: Mat3I,
    translation: Vec3R,
}

impl ExactSeitz {
    /// Construct from an integer rotation and an exact translation.
    pub const fn new(rotation: Mat3I, translation: Vec3R) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    /// The identity operation.
    pub fn identity() -> Self {
        Self::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], Vec3R::zero())
    }

    /// Rotation part.
    pub const fn rotation(&self) -> Mat3I {
        self.rotation
    }

    /// Translation part (as stored, not reduced).
    pub const fn translation(&self) -> &Vec3R {
        &self.translation
    }

    /// Apply the operation to an exact point.
    pub fn apply(&self, point: &Vec3R) -> Result<Vec3R, SubductionError> {
        Mat3R::from_ints(self.rotation)
            .checked_mul_vector(point)?
            .checked_add(&self.translation)
    }

    /// Compose: apply `other` first, then `self`.
    pub fn compose(&self, other: &Self) -> Result<Self, SubductionError> {
        let rotation = Mat3R::from_ints(self.rotation)
            .checked_mul(&Mat3R::from_ints(other.rotation))?
            .to_int_matrix()?;
        let translation = Mat3R::from_ints(self.rotation)
            .checked_mul_vector(&other.translation)?
            .checked_add(&self.translation)?;
        Ok(Self::new(rotation, translation))
    }

    /// Inverse operation.
    pub fn inverse(&self) -> Result<Self, SubductionError> {
        let rotation_matrix = Mat3R::from_ints(self.rotation).inverse()?;
        let rotation = rotation_matrix.to_int_matrix()?;
        let translation = rotation_matrix
            .checked_mul_vector(&self.translation)?
            .checked_neg()?;
        Ok(Self::new(rotation, translation))
    }

    /// The same operation with its translation canonical modulo `lattice`.
    pub fn reduce(&self, lattice: &Lattice) -> Result<Self, SubductionError> {
        let reduction = lattice.reduce(&self.translation)?;
        Ok(Self::new(self.rotation, reduction.representative))
    }
}

/// Exact affine map of the subgroup conventional frame into the parent
/// conventional frame: `x_G = T x_H + o`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeitzTransform {
    matrix: Mat3R,
    origin: Vec3R,
}

impl SeitzTransform {
    /// Construct from a matrix and an origin shift.
    pub const fn new(matrix: Mat3R, origin: Vec3R) -> Self {
        Self { matrix, origin }
    }

    /// The identity map.
    pub const fn identity() -> Self {
        Self {
            matrix: Mat3R::identity(),
            origin: Vec3R::zero(),
        }
    }

    /// `T`.
    pub const fn matrix(&self) -> &Mat3R {
        &self.matrix
    }

    /// `o`.
    pub const fn origin(&self) -> &Vec3R {
        &self.origin
    }

    /// `x_G = T x_H + o`.
    pub fn map_point(&self, point: &Vec3R) -> Result<Vec3R, SubductionError> {
        self.matrix
            .checked_mul_vector(point)?
            .checked_add(&self.origin)
    }

    /// `R_G = T R_H T^-1`, rejected unless the image is an integer matrix.
    pub fn map_rotation(&self, rotation: Mat3I) -> Result<Mat3I, SubductionError> {
        let image = self
            .matrix
            .checked_mul(&Mat3R::from_ints(rotation))?
            .checked_mul(&self.matrix.inverse()?)?;
        if !image.is_integral() {
            return Err(SubductionError::NonIntegralRotationImage {
                rotation: Box::new(image),
            });
        }
        image.to_int_matrix()
    }

    /// `t_G = T t_H + o - R_G o`.
    pub fn map_translation(
        &self,
        rotation_image: Mat3I,
        translation: &Vec3R,
    ) -> Result<Vec3R, SubductionError> {
        let rotated_origin = Mat3R::from_ints(rotation_image).checked_mul_vector(&self.origin)?;
        self.matrix
            .checked_mul_vector(translation)?
            .checked_add(&self.origin)?
            .checked_sub(&rotated_origin)
    }

    /// Map a complete operation into the parent frame.
    pub fn map_operation(&self, operation: &ExactSeitz) -> Result<ExactSeitz, SubductionError> {
        let rotation = self.map_rotation(operation.rotation)?;
        let translation = self.map_translation(rotation, &operation.translation)?;
        Ok(ExactSeitz::new(rotation, translation))
    }

    /// Inverse map `x_H = T^-1 (x_G - o)`.
    pub fn inverse(&self) -> Result<Self, SubductionError> {
        let matrix = self.matrix.inverse()?;
        let origin = matrix
            .checked_mul_vector(&self.origin)?
            .checked_neg()?;
        Ok(Self { matrix, origin })
    }

    /// Composition: `self` after `other`, i.e. `T = T_self T_other`,
    /// `o = T_self o_other + o_self`.
    pub fn compose(&self, other: &Self) -> Result<Self, SubductionError> {
        let matrix = self.matrix.checked_mul(&other.matrix)?;
        let origin = self
            .matrix
            .checked_mul_vector(&other.origin)?
            .checked_add(&self.origin)?;
        Ok(Self { matrix, origin })
    }
}

// ── Lattices ─────────────────────────────────────────────────────────────────

/// Lattice spanned by the rows of an invertible exact matrix.
///
/// Membership and reduction are always taken with respect to *this* lattice:
/// the parent lattice `L_G` decides whether an operation belongs to the parent
/// group, while the subgroup lattice `L_H` decides when two subgroup
/// operations are the same element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Lattice {
    rows: Mat3R,
    /// `(rows^T)^-1`: maps a point to its coordinates in this row basis.
    ///
    /// The transpose matters.  With the rows of `rows` as basis vectors, a
    /// point is `rows^T . coordinates`, so the coordinate map is
    /// `(rows^T)^-1` and *not* `rows^-1`; the two agree only for symmetric
    /// bases, which is why a diagonal-only test suite cannot tell them apart.
    coordinates: Mat3R,
}

impl Lattice {
    /// Build from row vectors, rejecting a singular basis.
    pub fn new(rows: Mat3R) -> Result<Self, SubductionError> {
        let coordinates = rows.inverse()?.transpose();
        Ok(Self { rows, coordinates })
    }

    /// The integer lattice `Z^3`.
    pub fn integer() -> Self {
        Self {
            rows: Mat3R::identity(),
            coordinates: Mat3R::identity(),
        }
    }

    /// Matrix that maps a point to its coordinates in this lattice basis.
    pub const fn coordinate_matrix(&self) -> &Mat3R {
        &self.coordinates
    }

    /// Row vectors.
    pub const fn rows(&self) -> &Mat3R {
        &self.rows
    }

    /// Exact determinant (signed volume).
    pub fn determinant(&self) -> Result<Rat, SubductionError> {
        self.rows.determinant()
    }

    /// Express `vector` in lattice coordinates: `vector = rows^T . coordinates`.
    pub fn coordinates(&self, vector: &Vec3R) -> Result<Vec3R, SubductionError> {
        self.coordinates.checked_mul_vector(vector)
    }

    /// The point with the given lattice coordinates: `rows^T . coordinates`.
    pub fn point(&self, coordinates: &Vec3R) -> Result<Vec3R, SubductionError> {
        self.rows.transpose().checked_mul_vector(coordinates)
    }

    /// Reduce `vector` modulo the lattice.
    ///
    /// Returns the representative inside the half-open fundamental cell
    /// (`0 <= coordinate < 1` in lattice coordinates) together with the integer
    /// combination of lattice rows that was removed:
    /// `vector = shift . rows + representative`.
    pub fn reduce(&self, vector: &Vec3R) -> Result<LatticeReduction, SubductionError> {
        let coordinates = self.coordinates(vector)?;
        let mut shift = [0i128; 3];
        let mut representative = *vector;
        for (axis, shift_entry) in shift.iter_mut().enumerate() {
            let floor = floor_rational(coordinates.get(axis))?;
            *shift_entry = floor;
            if floor != 0 {
                let mut row_multiple = [Rat::ZERO; 3];
                for (column, entry) in self.rows.row(axis).iter().enumerate() {
                    row_multiple[column] = entry.checked_mul(Rat::from_integer(floor))?;
                }
                representative = representative.checked_sub(&Vec3R::new(row_multiple))?;
            }
        }
        Ok(LatticeReduction {
            representative,
            shift,
        })
    }

    /// Whether `vector` lies in the lattice.
    pub fn contains(&self, vector: &Vec3R) -> Result<bool, SubductionError> {
        Ok(self.reduce(vector)?.representative.is_zero())
    }

    /// Whether two vectors agree modulo the lattice.
    pub fn same_mod(&self, left: &Vec3R, right: &Vec3R) -> Result<bool, SubductionError> {
        self.contains(&left.checked_sub(right)?)
    }

    /// Canonical representatives of `operations` modulo this lattice, keeping
    /// the first occurrence of each class.
    ///
    /// Deduplicating modulo the *parent* lattice instead is a real error: for a
    /// supercell embedding two distinct subgroup elements can differ by a
    /// parent lattice vector.
    pub fn deduplicate(&self, operations: &[ExactSeitz]) -> Result<Vec<ExactSeitz>, SubductionError> {
        let mut out: Vec<ExactSeitz> = Vec::with_capacity(operations.len());
        for operation in operations {
            let reduced = operation.reduce(self)?;
            if !out.contains(&reduced) {
                out.push(reduced);
            }
        }
        Ok(out)
    }
}

/// Result of reducing a vector modulo a lattice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LatticeReduction {
    /// Representative inside the half-open fundamental cell.
    pub representative: Vec3R,
    /// Integer coefficients of the removed lattice combination.
    pub shift: [i128; 3],
}

/// Floor of an exact rational.
fn floor_rational(value: Rat) -> Result<i128, SubductionError> {
    let quotient = value
        .numerator()
        .checked_div(value.denominator())
        .ok_or(SubductionError::RationalOverflow {
            operation: "floor",
        })?;
    let remainder = value
        .numerator()
        .checked_rem(value.denominator())
        .ok_or(SubductionError::RationalOverflow {
            operation: "floor",
        })?;
    if remainder < 0 {
        quotient
            .checked_sub(1)
            .ok_or(SubductionError::RationalOverflow { operation: "floor" })
    } else {
        Ok(quotient)
    }
}

// ── Strict Hall operation source ─────────────────────────────────────────────

/// Operations of one space group in its recorded data Hall setting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SgHallOperations {
    /// Space group number (1-230).
    pub sg: u8,
    /// Hall number taken from `SG_DATA_HALL`.
    pub hall: usize,
    /// Operations, translations snapped exactly to the source grid.
    pub operations: Vec<ExactSeitz>,
}

/// Load `sg`'s operations from the recorded data Hall setting, without any
/// fallback.
///
/// Unlike [`crate::irrep::bridge::canonical_hall_ops`], a missing or unreadable
/// Hall setting is an error: the isotropy tables and the character data are
/// indexed in this frame, so silently switching to the first Hall setting would
/// mix frames.
pub fn strict_sg_hall_ops(sg: u8) -> Result<SgHallOperations, SubductionError> {
    if sg == 0 || sg > 230 {
        return Err(SubductionError::InvalidSpaceGroup { sg });
    }
    let hall = SG_DATA_HALL[usize::from(sg)];
    if hall == 0 {
        return Err(SubductionError::StrictHallUnavailable { sg });
    }
    let hall = usize::from(hall);
    load_hall_operations(sg, hall)
}

fn load_hall_operations(sg: u8, hall: usize) -> Result<SgHallOperations, SubductionError> {
    let ops = SymmetryOps::from_database(hall).map_err(|source| SubductionError::HallLoadFailed {
        sg,
        hall,
        source,
    })?;
    let mut operations = Vec::with_capacity(ops.len());
    for op in ops.iter() {
        let translation = Vec3R::new([
            Rat::from_grid(op.translation[0], SOURCE_TRANSLATION_GRID)?,
            Rat::from_grid(op.translation[1], SOURCE_TRANSLATION_GRID)?,
            Rat::from_grid(op.translation[2], SOURCE_TRANSLATION_GRID)?,
        ]);
        operations.push(ExactSeitz::new(op.rotation, translation));
    }
    Ok(SgHallOperations {
        sg,
        hall,
        operations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rat(num: i128, den: i128) -> Rat {
        Rat::new(num, den).expect("valid rational")
    }

    fn vec3(values: [i32; 3]) -> Vec3R {
        Vec3R::from_ints(values)
    }

    #[test]
    fn rational_arithmetic_is_exact() {
        assert_eq!(rat(1, 16).checked_add(rat(1, 16)).unwrap(), rat(1, 8));
        assert_eq!(rat(1, 3).checked_mul(rat(3, 1)).unwrap(), Rat::ONE);
        assert_eq!(rat(5, 2).checked_sub(rat(1, 2)).unwrap(), rat(2, 1));
        assert_eq!(rat(1, -2), rat(-1, 2));
        assert_eq!(rat(4, 6), rat(2, 3));
        assert_eq!(rat(-1, 16).numerator(), -1);
        assert_eq!(rat(-1, 16).denominator(), 16);
        // 1/16 is exactly representable: nothing is truncated to twelfths.
        assert_eq!(rat(1, 16).checked_mul(rat(16, 1)).unwrap(), Rat::ONE);
        assert!(rat(1, 16).to_integer().is_err());
        assert_eq!(rat(-7, 1).to_i32().unwrap(), -7);
        assert_eq!(rat(3, 2).to_i32(), Err(SubductionError::NotAnInteger { value: rat(3, 2) }));
        assert_eq!(
            Rat::from_integer(i128::from(i32::MAX) + 1).to_i32(),
            Err(SubductionError::IntegerOutOfRange {
                value: i128::from(i32::MAX) + 1
            })
        );
    }

    #[test]
    fn rational_errors_are_reported_not_wrapped() {
        assert_eq!(Rat::new(1, 0), Err(SubductionError::ZeroDenominator));
        assert_eq!(
            rat(1, 2).checked_div(Rat::ZERO),
            Err(SubductionError::DivisionByZero)
        );
        let big = Rat::from_integer(i128::MAX);
        assert!(matches!(
            big.checked_add(big),
            Err(SubductionError::RationalOverflow { .. })
        ));
        assert!(matches!(
            big.checked_mul(big),
            Err(SubductionError::RationalOverflow { .. })
        ));
        assert!(matches!(
            Rat::from_integer(i128::MIN).checked_neg(),
            Err(SubductionError::RationalOverflow { .. })
        ));
        assert!(matches!(
            rat(1, 2).checked_add(Rat::from_integer(i128::MAX)),
            Err(SubductionError::RationalOverflow { .. })
        ));
    }

    #[test]
    fn grid_snap_is_exact_and_rejects_off_grid_values() {
        assert_eq!(Rat::from_grid(1.0 / 3.0, 12).unwrap(), rat(1, 3));
        assert_eq!(Rat::from_grid(-1e-16, 12).unwrap(), Rat::ZERO);
        assert_eq!(Rat::from_grid(0.25 + 1e-12, 12).unwrap(), rat(1, 4));
        assert_eq!(Rat::from_grid(5.0 / 2.0, 12).unwrap(), rat(5, 2));
        assert_eq!(Rat::from_grid(1.0 / 3.0 + 1e-13, 12).unwrap(), rat(1, 3));
        // Values that are merely *close* to a grid point are still rejected.
        assert!(matches!(
            Rat::from_grid(1.0 / 3.0 + 1e-6, 12),
            Err(SubductionError::OffGridValue { .. })
        ));
        assert!(matches!(
            Rat::from_grid(0.1, 12),
            Err(SubductionError::OffGridValue { .. })
        ));
        assert!(matches!(
            Rat::from_grid(f64::NAN, 12),
            Err(SubductionError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            Rat::from_grid(f64::INFINITY, 12),
            Err(SubductionError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            Rat::from_grid(1.0, 0),
            Err(SubductionError::InvalidGrid { .. })
        ));
    }

    #[test]
    fn matrix_inverse_determinant_and_integrality() {
        let asymmetric = Mat3R::from_ints([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]);
        assert_eq!(asymmetric.determinant().unwrap(), rat(2, 1));
        let inverse = asymmetric.inverse().unwrap();
        assert_eq!(
            inverse,
            Mat3R::new([
                [rat(1, 2), rat(-1, 2), Rat::ZERO],
                [rat(1, 2), rat(1, 2), Rat::ZERO],
                [Rat::ZERO, Rat::ZERO, Rat::ONE],
            ])
        );
        assert_eq!(asymmetric.checked_mul(&inverse).unwrap(), Mat3R::identity());
        assert!(!inverse.is_integral());
        assert_eq!(
            Mat3R::identity().to_int_matrix().unwrap(),
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        );
        assert_eq!(Mat3R::identity().transpose(), Mat3R::identity());
        let singular = Mat3R::diagonal([1, 1, 0]);
        assert_eq!(singular.determinant().unwrap(), Rat::ZERO);
        assert_eq!(singular.inverse(), Err(SubductionError::SingularMatrix));
        assert!(matches!(
            inverse.to_int_matrix(),
            Err(SubductionError::NotAnInteger { .. })
        ));
    }

    #[test]
    fn transform_round_trips_include_asymmetric_and_negative_determinant() {
        let transforms = [
            // 221 GM4+ P2 -> #12 C2/m: asymmetric, det = 2.
            SeitzTransform::new(
                Mat3R::from_ints([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]),
                vec3([0, 0, 0]),
            ),
            // Axis swap: det = -1, with an origin shift.
            SeitzTransform::new(
                Mat3R::from_ints([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
                vec3([1, 0, 0]),
            ),
            // Thirds appear in the 167 GM3+ P1 embedding (det = -2/9).
            SeitzTransform::new(
                Mat3R::new([
                    [rat(-1, 3), rat(2, 3), Rat::ZERO],
                    [rat(1, 3), rat(1, 3), rat(-1, 3)],
                    [Rat::ONE, Rat::ZERO, Rat::ZERO],
                ]),
                Vec3R::new([rat(1, 16), rat(3, 16), rat(-1, 8)]),
            ),
        ];
        let points = [
            vec3([0, 0, 0]),
            Vec3R::new([rat(1, 2), rat(-1, 4), rat(3, 16)]),
            Vec3R::new([rat(7, 3), rat(1, 6), rat(-5, 12)]),
        ];
        for transform in transforms {
            let inverse = transform.inverse().unwrap();
            assert_eq!(
                transform.compose(&inverse).unwrap(),
                SeitzTransform::identity()
            );
            assert_eq!(
                inverse.compose(&transform).unwrap(),
                SeitzTransform::identity()
            );
            for point in points {
                let mapped = transform.map_point(&point).unwrap();
                assert_eq!(inverse.map_point(&mapped).unwrap(), point);
            }
        }
    }

    #[test]
    fn rotation_images_must_be_integral() {
        // C4z in a frame stretched by 2 along x is not a lattice symmetry.
        let stretched = SeitzTransform::new(Mat3R::diagonal([2, 1, 1]), Vec3R::zero());
        assert!(matches!(
            stretched.map_rotation([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
            Err(SubductionError::NonIntegralRotationImage { .. })
        ));
        // The 221 GM4+ P1 embedding maps the subgroup's own C4z to C4x.
        let transform = SeitzTransform::new(
            Mat3R::from_ints([[0, 0, 1], [0, -1, 0], [1, 0, 0]]),
            Vec3R::zero(),
        );
        assert_eq!(
            transform
                .map_rotation([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
                .unwrap(),
            [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
        );
    }

    #[test]
    fn sixteenth_denominators_survive_a_transform() {
        let transform = SeitzTransform::new(
            Mat3R::identity(),
            Vec3R::new([rat(1, 16), rat(3, 16), rat(-1, 8)]),
        );
        // Identity rotation: `o` cancels and the 1/16 is carried through.
        let mapped = transform
            .map_operation(&ExactSeitz::new(
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                Vec3R::new([rat(1, 16), Rat::ZERO, Rat::ZERO]),
            ))
            .unwrap();
        assert_eq!(mapped.translation().as_array(), &[rat(1, 16), Rat::ZERO, Rat::ZERO]);
        // A mirror combines it with the origin instead of cancelling it.
        let mapped = transform
            .map_operation(&ExactSeitz::new(
                [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                Vec3R::new([Rat::ZERO, Rat::ZERO, rat(1, 16)]),
            ))
            .unwrap();
        assert_eq!(
            mapped.translation().as_array(),
            &[Rat::ZERO, rat(3, 8), rat(1, 16)]
        );
        // Mixed denominators in one operation.  With the origin above and an
        // off-diagonal mirror:
        //   t_G = T t + o - R o = (0,0,1/12) + (1/16,3/16,-1/8) - (3/16,1/16,1/8)
        //       = (-1/8, 1/8, -1/6),
        // so eighths, sixteenths and twelfths all pass through unchanged
        // (a twelfth grid truncation, or rounding to 1/12, would corrupt it).
        let mapped = transform
            .map_operation(&ExactSeitz::new(
                [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
                Vec3R::new([Rat::ZERO, Rat::ZERO, rat(1, 12)]),
            ))
            .unwrap();
        assert_eq!(
            mapped.translation().as_array(),
            &[rat(-1, 8), rat(1, 8), rat(-1, 6)]
        );
        assert_eq!(mapped.translation().get(0).denominator(), 8);
        assert_eq!(mapped.translation().get(2).denominator(), 6);
        // The origin denominators themselves are preserved by the inverse map:
        // x_H = T^-1 (x_G - o) must not round the 1/16 away.
        let inverse = transform.inverse().unwrap();
        assert_eq!(
            inverse.map_point(&vec3([0, 0, 0])).unwrap().as_array(),
            &[rat(-1, 16), rat(-3, 16), rat(1, 8)]
        );
    }

    #[test]
    fn lattice_reduction_returns_representative_and_removed_shift() {
        let superlattice = Lattice::new(Mat3R::diagonal([2, 2, 2])).unwrap();
        assert_eq!(superlattice.determinant().unwrap(), rat(8, 1));
        let reduction = superlattice.reduce(&vec3([3, 0, 0])).unwrap();
        assert_eq!(reduction.representative, vec3([1, 0, 0]));
        assert_eq!(reduction.shift, [1, 0, 0]);
        assert!(!superlattice.contains(&vec3([1, 0, 0])).unwrap());
        assert!(superlattice.contains(&vec3([0, -2, 4])).unwrap());
        // v = shift . rows + representative must hold exactly ...
        let value = Vec3R::new([rat(7, 2), rat(-3, 2), rat(1, 16)]);
        let reduction = superlattice.reduce(&value).unwrap();
        let rows = *superlattice.rows();
        let mut reconstructed = reduction.representative;
        for (axis, coefficient) in reduction.shift.iter().enumerate() {
            let mut row_multiple = [Rat::ZERO; 3];
            for (column, entry) in rows.row(axis).iter().enumerate() {
                row_multiple[column] = entry.checked_mul(Rat::from_integer(*coefficient)).unwrap();
            }
            reconstructed = reconstructed
                .checked_add(&Vec3R::new(row_multiple))
                .unwrap();
        }
        assert_eq!(reconstructed, value);
        // ... and the representative must lie in the half-open cell.
        for coordinate in superlattice
            .coordinates(&reduction.representative)
            .unwrap()
            .as_array()
        {
            assert!(coordinate.numerator() >= 0);
            assert!(coordinate.numerator() < coordinate.denominator());
        }
        assert_eq!(
            Lattice::new(Mat3R::diagonal([1, 1, 0])),
            Err(SubductionError::SingularMatrix)
        );
    }

    #[test]
    fn lattice_coordinates_use_the_transposed_inverse() {
        // Non-symmetric, non-diagonal basis: the C-centring of the 167 GM3+ P1
        // embedding.  `rows^-1` and `(rows^T)^-1` differ here, so a transposed
        // coordinate map cannot pass this test.
        let lattice = Lattice::new(Mat3R::new([
            [rat(1, 2), rat(1, 2), Rat::ZERO],
            [rat(-1, 2), rat(1, 2), Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ]))
        .unwrap();
        // (C^-1)^T, not C^-1 = [[1,-1,0],[1,1,0],[0,0,1]].
        assert_eq!(
            lattice.coordinate_matrix(),
            &Mat3R::from_ints([[1, 1, 0], [-1, 1, 0], [0, 0, 1]])
        );
        assert_eq!(
            lattice.coordinates(&vec3([0, 1, 0])).unwrap(),
            vec3([1, 1, 0])
        );
        assert_eq!(
            lattice.coordinates(&vec3([1, 0, 0])).unwrap(),
            Vec3R::from_ints([1, -1, 0])
        );
        for point in [
            vec3([1, 0, 0]),
            vec3([0, 1, 0]),
            Vec3R::new([rat(1, 3), rat(-2, 5), rat(3, 4)]),
        ] {
            let coordinates = lattice.coordinates(&point).unwrap();
            assert_eq!(lattice.point(&coordinates).unwrap(), point);
        }
        // Integer coordinates are exactly the lattice.
        assert!(lattice.contains(&vec3([1, 0, 0])).unwrap());
        assert!(lattice.contains(&vec3([0, 1, 0])).unwrap());
        assert!(lattice.contains(&vec3([0, 0, 1])).unwrap());
        assert!(!lattice
            .contains(&Vec3R::new([rat(1, 3), Rat::ZERO, Rat::ZERO]))
            .unwrap());
        let half = Vec3R::new([rat(1, 2), rat(1, 2), Rat::ZERO]);
        assert!(lattice.contains(&half).unwrap());
        assert_eq!(
            lattice.coordinates(&half).unwrap(),
            Vec3R::from_ints([1, 0, 0])
        );
        // Reduction returns the integer combination it removed.
        let reduction = lattice.reduce(&vec3([0, 1, 0])).unwrap();
        assert!(reduction.representative.is_zero());
        assert_eq!(reduction.shift, [1, 1, 0]);
        let reduction = lattice
            .reduce(&Vec3R::new([Rat::ZERO, Rat::ZERO, rat(1, 2)]))
            .unwrap();
        assert_eq!(reduction.shift, [0, 0, 0]);
        assert_eq!(
            reduction.representative,
            Vec3R::new([Rat::ZERO, Rat::ZERO, rat(1, 2)])
        );
    }

    #[test]
    fn lattice_reduction_is_exact_for_a_centred_non_orthogonal_basis() {
        // Child lattice of 167 GM3+ P1 -> #15 C2/c in parent conventional
        // coordinates: rows = C_centring . B_printed (an R lattice).
        let child = Lattice::new(Mat3R::new([
            [rat(-1, 3), rat(1, 3), rat(1, 3)],
            [rat(-2, 3), rat(-1, 3), rat(-1, 3)],
            [rat(-1, 3), rat(-2, 3), rat(1, 3)],
        ]))
        .unwrap();
        assert_eq!(child.determinant().unwrap(), rat(1, 3));
        // The printed inversion translation is a primitive lattice vector of
        // this lattice, so it must reduce to zero.
        let centring = Vec3R::new([rat(-1, 3), rat(1, 3), rat(1, 3)]);
        assert!(child.contains(&centring).unwrap());
        for point in [
            Vec3R::new([rat(1, 6), rat(1, 3), rat(1, 3)]),
            Vec3R::new([rat(-5, 6), rat(2, 3), rat(-1, 2)]),
        ] {
            let reduction = child.reduce(&point).unwrap();
            let reconstructed = child
                .point(&Vec3R::new(reduction.shift.map(Rat::from_integer)))
                .unwrap()
                .checked_add(&reduction.representative)
                .unwrap();
            assert_eq!(reconstructed, point);
        }
    }

    #[test]
    fn lattice_membership_is_per_lattice() {
        let parent = Lattice::integer();
        let subgroup = Lattice::new(Mat3R::diagonal([2, 2, 2])).unwrap();
        let rotation = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let operation = ExactSeitz::new(rotation, Vec3R::zero());
        let shifted = ExactSeitz::new(rotation, Vec3R::from_ints([1, 0, 0]));
        // Equal modulo L_G ...
        assert!(parent
            .same_mod(operation.translation(), shifted.translation())
            .unwrap());
        // ... but different modulo L_H, so they must not be merged.
        assert!(!subgroup
            .same_mod(operation.translation(), shifted.translation())
            .unwrap());
        assert_eq!(parent.deduplicate(&[operation, shifted]).unwrap().len(), 1);
        assert_eq!(subgroup.deduplicate(&[operation, shifted]).unwrap().len(), 2);
    }

    #[test]
    fn operation_composition_and_inverse_are_exact() {
        let lattice = Lattice::integer();
        let operation = ExactSeitz::new(
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            Vec3R::from_ints([1, 0, 0]),
        );
        let inverse = operation.inverse().unwrap();
        let identity = operation.compose(&inverse).unwrap();
        assert_eq!(identity.rotation(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        assert!(lattice
            .same_mod(identity.translation(), &Vec3R::zero())
            .unwrap());
        assert_eq!(operation.apply(&Vec3R::zero()).unwrap(), vec3([1, 0, 0]));
        let point = Vec3R::new([rat(1, 3), rat(2, 5), rat(-1, 7)]);
        let mapped = operation.apply(&point).unwrap();
        assert_eq!(inverse.apply(&mapped).unwrap(), point);
        // Applying the operation twice equals composing it with itself.
        assert_eq!(
            operation.apply(&mapped).unwrap(),
            operation.compose(&operation).unwrap().apply(&point).unwrap()
        );
    }

    #[test]
    fn strict_hall_source_covers_all_space_groups_on_the_source_grid() {
        let mut operations_total = 0usize;
        for sg in 1..=230u8 {
            let loaded = strict_sg_hall_ops(sg).unwrap_or_else(|error| {
                panic!("SG {sg}: {error}");
            });
            assert_eq!(loaded.sg, sg);
            assert_eq!(loaded.hall, usize::from(SG_DATA_HALL[usize::from(sg)]));
            assert!(!loaded.operations.is_empty(), "SG {sg} has no operations");
            for operation in &loaded.operations {
                for value in operation.translation().as_array() {
                    assert!(value.denominator() > 0);
                    assert!(
                        SOURCE_TRANSLATION_GRID % value.denominator() == 0,
                        "SG {sg} has translation {value:?} off the \
                         1/{SOURCE_TRANSLATION_GRID} grid"
                    );
                }
            }
            operations_total += loaded.operations.len();
        }
        // Pinned against the shipped tables (230 conventional settings,
        // centring expansions included); a change here means the Hall source
        // moved, which is exactly what this layer must notice.
        assert_eq!(operations_total, 4425);
    }

    #[test]
    fn strict_hall_source_rejects_bad_input_without_falling_back() {
        assert_eq!(
            strict_sg_hall_ops(0),
            Err(SubductionError::InvalidSpaceGroup { sg: 0 })
        );
        assert_eq!(
            strict_sg_hall_ops(231),
            Err(SubductionError::InvalidSpaceGroup { sg: 231 })
        );
        // Hall 0 is a real table entry for the dummy index 0; loading it must
        // fail loudly rather than select the first Hall setting of some SG.
        assert!(matches!(
            load_hall_operations(1, 0),
            Err(SubductionError::HallLoadFailed { sg: 1, hall: 0, .. })
        ));
        assert_eq!(
            SubductionError::StrictHallUnavailable { sg: 42 }.to_string(),
            "space group 42 has no data Hall setting; no fallback is attempted"
        );
    }
}
