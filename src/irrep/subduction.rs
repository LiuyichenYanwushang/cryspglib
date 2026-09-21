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
use crate::irrep::isotropy::{IsotropySubgroup, parent_primitive_basis};
use crate::irrep::query;
use crate::irrep::types::generated_data::{ISOTROPY_SUBGROUPS, SG_DATA_HALL};
use crate::irrep::types::{
    CharacterRow, CompoundSelectedArmCharacter, IrrepRecord, IsotropyRecord, KVector,
    SeitzOperation,
};
use crate::mathfunc::Mat3I;
use num_complex::Complex64;

#[path = "subduction_star.rs"]
pub mod star;

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
    /// The isotropy ordinal is outside the generated table.
    #[error("isotropy record {ordinal} is outside the table of {len} records")]
    IsotropyOrdinalOutOfRange { ordinal: usize, len: usize },
    /// The record handed in is not the one stored at that ordinal, or the
    /// ordinal is not listed for the claimed (parent, irrep) context.
    ///
    /// `IsotropySubgroup` has public fields, so a caller can build one by hand;
    /// the embedding never trusts those fields without rechecking them against
    /// the generated table.
    #[error("isotropy record {ordinal} does not match the generated table")]
    StaleIsotropyRecord { ordinal: usize },
    /// A cached embedding and the isotropy record handed in with it belong to
    /// different contexts.
    ///
    /// Reusing one [`SubgroupEmbedding`] across probes of the *same* record is
    /// the point of [`subduce_irrep_with_embedding`]; reusing it for another
    /// record silently mixes two contexts (subgroup number, setting and child
    /// origin shift all come from the cached record), so the mismatch is
    /// rejected instead of folded.
    #[error(
        "embedding context mismatch: subgroup (parent {parent_sg}, ordinal {ordinal}, \
         subgroup {subgroup_sg}) but the cached embedding is (parent {embedding_parent_sg}, \
         ordinal {embedding_ordinal}, subgroup {embedding_subgroup_sg})"
    )]
    EmbeddingContextMismatch {
        parent_sg: u8,
        ordinal: usize,
        subgroup_sg: u8,
        embedding_parent_sg: u8,
        embedding_ordinal: usize,
        embedding_subgroup_sg: u8,
    },
    /// The probe is not a member of the parent's generated irrep table.
    ///
    /// This staged API requires a reference into the generated `IRREPS` table:
    /// a copied record can retain its source ID while its public fields change.
    /// Only a record obtained from
    /// [`query::irreps_of(parent_sg)`](crate::irrep::query::irreps_of) can be
    /// interpreted in the parent's frame.
    #[error("probe {ml} is not a record of the irrep table of space group {sg}")]
    ForeignProbe { sg: u8, ml: &'static str },
    /// The parent space group has no scalar irrep with the requested label.
    #[error("space group {sg} has no irrep with label {ml}")]
    UnknownParentIrrep { sg: u8, ml: &'static str },
    /// The subgroup number stored in the record is not a space group number.
    #[error("subgroup number {sg} is outside 1-230")]
    InvalidSubgroupNumber { sg: usize },
    /// No setting candidate reproduced the subgroup inside the parent.
    #[error("no setting of subgroup {subgroup_sg} embeds into the parent ({candidates} tried)")]
    NoValidEmbedding { subgroup_sg: u8, candidates: usize },
    /// The parent has no irrep with the requested label at Gamma.
    #[error("space group {sg} has no Gamma irrep {ml}")]
    ProbeIrrepNotFound { sg: u8, ml: String },
    /// The requested probe irrep is not at Gamma; that needs the k-folding and
    /// Bloch-phase stage (task 7).
    #[error("irrep {ml} of space group {sg} is not at Gamma (k = {k:?})")]
    ProbeNotAtGamma { sg: u8, ml: String, k: KVector },
    /// The representation space is not the scalar ordinary one (compound rows
    /// are task 6, spinors task 11).
    #[error("irrep {ml} of space group {sg} has an unsupported character space")]
    UnsupportedCharacterSpace { sg: u8, ml: String },
    /// An operation of the embedding has no counterpart in a character row.
    #[error("no character value for operation {index} ({rotation:?}) in {ml}")]
    OperationNotInCharacterRow {
        ml: &'static str,
        index: usize,
        rotation: [i32; 9],
    },
    /// Matching the same operation twice with different values.
    #[error("character row {ml} disagrees with itself on operation {index}")]
    InconsistentCharacterRow { ml: &'static str, index: usize },
    /// The computed multiplicity is not a non-negative integer.
    #[error("multiplicity of {ml} is {value} (not a non-negative integer)")]
    NonIntegralMultiplicity { ml: &'static str, value: Complex64 },
    /// The target rows are not a usable (orthogonal) set.
    #[error("target rows {first} and {second} are not orthogonal ({value})")]
    TargetRowsNotOrthogonal {
        first: &'static str,
        second: &'static str,
        value: Complex64,
    },
    /// The multiplicities do not reproduce the subduced dimension.
    #[error("decomposed dimension {found} does not equal the parent dimension {expected}")]
    DimensionSumMismatch { expected: u32, found: i64 },
    /// The multiplicities do not reproduce the parent characters.
    #[error("character reconstruction failed at operation {index}: {found} != {expected}")]
    CharacterMismatch {
        index: usize,
        found: Complex64,
        expected: Complex64,
    },
    /// Two rows produced the same canonical source identity.
    #[error("target source {ml} (irnumber {irnumber}) appears twice")]
    DuplicateTargetSource { ml: &'static str, irnumber: u32 },
    /// A complex target row is not irreducible over the coset set.
    #[error("target {ml} has norm {norm} instead of 1; it is not a complex irrep")]
    TargetNotIrreducible { ml: &'static str, norm: f64 },
    /// The folded wave vector is not among the subgroup's stored k-points.
    ///
    /// The reduction is done by exact rational arithmetic and matched modulo the
    /// subgroup's reciprocal lattice (centring extinctions included); a nearest
    /// match or an arbitrary rounding is never used.
    #[error("subgroup {sg} has no irrep data at the folded k = ({}, {}, {})", k[0], k[1], k[2])]
    MissingIrrepData { sg: u8, k: [Rat; 3] },
    /// The parent star has more than one arm, so the subduced representation
    /// splits across several subgroup stars (task 8).
    #[error("parent star of k = ({}, {}, {}) has several arms", k[0], k[1], k[2])]
    UnsupportedMultiArmStar { sg: u8, k: [Rat; 3] },
    /// A frozen setting did not reproduce the subgroup inside the parent.
    #[error(
        "recorded setting {setting:?}/{denominator} for subgroup {subgroup_sg} failed validation"
    )]
    FrozenEmbeddingRejected {
        subgroup_sg: u8,
        setting: Mat3I,
        denominator: i32,
    },
    /// Several setting candidates are consistent with the parent group, so the
    /// label correspondence cannot be decided from the stored data alone.
    #[error(
        "subgroup {subgroup_sg} has {candidates} consistent settings; the canonical one \
         needs a recorded setting transform"
    )]
    AmbiguousEmbedding { subgroup_sg: u8, candidates: usize },
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
                num.checked_neg().ok_or(SubductionError::RationalOverflow {
                    operation: "normalize",
                })?,
                den.checked_neg().ok_or(SubductionError::RationalOverflow {
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
            i32::try_from(denominator).map_err(|_| SubductionError::InvalidGrid { denominator })?,
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

    /// Lossy `f64` value, for trigonometric phase evaluation only.
    ///
    /// Every comparison in this module is exact; this exists solely to turn an
    /// exact phase angle into a unit complex number.
    pub fn to_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }

    /// Checked addition.
    pub fn checked_add(self, other: Self) -> Result<Self, SubductionError> {
        let shared = gcd_positive(self.den.unsigned_abs(), other.den.unsigned_abs());
        let shared = i128::try_from(shared)
            .map_err(|_| SubductionError::RationalOverflow { operation: "add" })?;
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
        let left = i128::try_from(left).map_err(|_| SubductionError::RationalOverflow {
            operation: "multiply",
        })?;
        let right = i128::try_from(right).map_err(|_| SubductionError::RationalOverflow {
            operation: "multiply",
        })?;
        let num = (self.num / left).checked_mul(other.num / right).ok_or(
            SubductionError::RationalOverflow {
                operation: "multiply",
            },
        )?;
        let den = (self.den / right).checked_mul(other.den / left).ok_or(
            SubductionError::RationalOverflow {
                operation: "multiply",
            },
        )?;
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
            .ok_or(SubductionError::RationalOverflow {
                operation: "divide",
            })?;
        let den = self
            .den
            .checked_mul(other.num)
            .ok_or(SubductionError::RationalOverflow {
                operation: "divide",
            })?;
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
        Ok([
            self.0[0].to_i32()?,
            self.0[1].to_i32()?,
            self.0[2].to_i32()?,
        ])
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
        Self::from_ints([[values[0], 0, 0], [0, values[1], 0], [0, 0, values[2]]])
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
                let first = self.0[rows[0]][columns[0]].checked_mul(self.0[rows[1]][columns[1]])?;
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

    /// Map an operation **back** into the subgroup frame.
    ///
    /// `R_H = T^-1 R_G T` (rejected unless integral) and
    /// `t_H = T^-1 (t_G + R_G o - o)`, the inverse of [`Self::map_operation`].
    pub fn unmap_operation(&self, operation: &ExactSeitz) -> Result<ExactSeitz, SubductionError> {
        let inverse = self.inverse()?;
        let rotation = inverse
            .matrix()
            .checked_mul(&Mat3R::from_ints(operation.rotation()))?
            .checked_mul(&self.matrix)?;
        if !rotation.is_integral() {
            return Err(SubductionError::NonIntegralRotationImage {
                rotation: Box::new(rotation),
            });
        }
        let rotation = rotation.to_int_matrix()?;
        let rotated_origin =
            Mat3R::from_ints(operation.rotation()).checked_mul_vector(&self.origin)?;
        let translation =
            inverse.map_point(&operation.translation().checked_add(&rotated_origin)?)?;
        Ok(ExactSeitz::new(rotation, translation))
    }

    /// Inverse map `x_H = T^-1 (x_G - o)`.
    pub fn inverse(&self) -> Result<Self, SubductionError> {
        let matrix = self.matrix.inverse()?;
        let origin = matrix.checked_mul_vector(&self.origin)?.checked_neg()?;
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

    /// The reciprocal lattice of this direct lattice.
    ///
    /// Row vectors are `(L^-1)^T`, so a vector `h` is a reciprocal lattice
    /// vector exactly when `h . t` is an integer for every direct lattice
    /// vector `t`.  For a centred cell this keeps the centring extinctions
    /// (`h1 + h2` even for C-centring, and so on), which is why wave vector
    /// equality must use this lattice rather than per-coordinate reduction.
    pub fn reciprocal(&self) -> Result<Self, SubductionError> {
        Self::new(self.rows.inverse()?.transpose())
    }

    /// Whether a direct-space rotation `R` fixes the wave vector `k` modulo
    /// **this** lattice.
    ///
    /// Convention, because the two matrices agree only for orthogonal `R`:
    ///
    /// * `self` is a **reciprocal** lattice `L*` (the caller passes
    ///   `Lattice::reciprocal()` of the direct lattice `L` whose symmetry
    ///   operations are being tested).  Membership is decided by `self` so a
    ///   centred cell keeps its extinctions: for C-centring the difference must
    ///   be an integer combination of `(1,1,0)`, `(-1,1,0)`, `(0,0,1)`, not of
    ///   the primitive `Z^3` basis.
    /// * `rotation` is a **direct-space column-action** rotation, `x' = R x`,
    ///   exactly as stored in the Hall operations and in [`ExactSeitz`].
    /// * `wave_vector` is in reciprocal-space coordinates (the same coordinates
    ///   `k_H = T^T k_G` folds, i.e. conjugate to the direct fractional
    ///   coordinates through `k . x`).
    ///
    /// The action on those coordinates is the contragredient one,
    /// `k' = R^-T k`, *not* `R k`: `R k` is the action on direct-space
    /// coordinates and is wrong for wave vectors.  A hexagonal `C3` such as
    /// SG 143's `[[0,-1,0],[1,-1,0],[0,0,1]]` is not orthogonal, and there the
    /// two differ: at `K = (1/3,1/3,0)` the contragredient image differs from
    /// `k` by `(-1,0,0) ∈ L*`, while `R k` differs by `(-2/3,-1/3,0) ∉ L*`.
    pub fn preserves(&self, rotation: Mat3I, wave_vector: &Vec3R) -> Result<bool, SubductionError> {
        let action = Mat3R::from_ints(rotation).inverse()?.transpose();
        let image = action.checked_mul_vector(wave_vector)?;
        self.contains(&image.checked_sub(wave_vector)?)
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
    pub fn deduplicate(
        &self,
        operations: &[ExactSeitz],
    ) -> Result<Vec<ExactSeitz>, SubductionError> {
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
        .ok_or(SubductionError::RationalOverflow { operation: "floor" })?;
    let remainder = value
        .numerator()
        .checked_rem(value.denominator())
        .ok_or(SubductionError::RationalOverflow { operation: "floor" })?;
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
    let ops = SymmetryOps::from_database(hall)
        .map_err(|source| SubductionError::HallLoadFailed { sg, hall, source })?;
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

// ── Subgroup embedding ───────────────────────────────────────────────────────

/// The generated per-ordinal setting metadata (task 9).
///
/// One entry per covered isotropy record ordinal, derived from the official
/// ISOTROPY program by `scripts/generate_subduction_settings.py`; the module
/// header records the pinned archive hashes and the derivation.  The table is
/// keyed by ordinal, so the lookup can fail closed on a parent/subgroup
/// mismatch instead of validating a setting that was recorded for another
/// embedding.
#[path = "subduction_settings_data.rs"]
mod settings_data;

use settings_data::FROZEN_EMBEDDING_SETTINGS;

/// The identity setting transform.
///
/// The generated table carries its own constant, so this test-only spelling
/// exists for the star-module tests that compare an embedding against the
/// identity setting.
#[cfg(test)]
const IDENTITY_SETTING: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// The recorded setting of the isotropy record `ordinal`, if the fixtures cover
/// it.
///
/// A covered ordinal whose pinned record is a different `(parent, subgroup)`
/// pair is a hard error, never a reason to fall back to the candidate search:
/// the same ordinal under another parent is a different embedding, and a
/// setting that validates there can still pair the irreps wrongly.
fn frozen_setting_for(
    ordinal: usize,
    parent_sg: u8,
    subgroup_sg: u8,
) -> Result<Option<&'static settings_data::FrozenEmbeddingSetting>, SubductionError> {
    let Some(entry) = FROZEN_EMBEDDING_SETTINGS
        .iter()
        .find(|entry| entry.0 == ordinal)
    else {
        return Ok(None);
    };
    if entry.1 != parent_sg || entry.2 != subgroup_sg {
        return Err(SubductionError::StaleIsotropyRecord { ordinal });
    }
    Ok(Some(entry))
}

/// The 48 signed permutation matrices, the search space for `U`.
pub fn signed_permutations() -> Vec<Mat3I> {
    let mut out = Vec::with_capacity(48);
    for permutation in [
        [0usize, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        for signs in [
            [1i32, 1, 1],
            [1, 1, -1],
            [1, -1, 1],
            [1, -1, -1],
            [-1, 1, 1],
            [-1, 1, -1],
            [-1, -1, 1],
            [-1, -1, -1],
        ] {
            let mut matrix = [[0i32; 3]; 3];
            for row in 0..3 {
                matrix[row][permutation[row]] = signs[row];
            }
            out.push(matrix);
        }
    }
    out
}

/// An isotropy subgroup embedded in its parent's conventional frame.
#[derive(Debug, Clone)]
pub struct SubgroupEmbedding {
    parent_sg: u8,
    subgroup_sg: u8,
    ordinal: usize,
    setting: Mat3I,
    setting_denominator: i32,
    child_shift: Vec3R,
    transform: SeitzTransform,
    parent_lattice: Lattice,
    subgroup_lattice: Lattice,
    operations: Vec<ExactSeitz>,
    representatives: Vec<ExactSeitz>,
    candidates: usize,
}

impl SubgroupEmbedding {
    /// Build and validate the embedding of one isotropy subgroup.
    ///
    /// The record's public fields are rechecked against the generated table,
    /// the stored basis `W` is converted through the parent's primitive basis,
    /// and every signed-permutation setting candidate is tested by mapping the
    /// subgroup's own Hall operations into the parent: the mapped operations
    /// must all be parent operations (modulo the parent lattice), their coset
    /// representatives (modulo the subgroup lattice) must close under
    /// multiplication and inversion, and their count must equal the subgroup's
    /// point group order.  A single surviving candidate is used; several
    /// survivors are resolved by the frozen embedding metadata, otherwise the
    /// construction reports ambiguity instead of guessing.
    pub fn from_isotropy_subgroup(subgroup: &IsotropySubgroup) -> Result<Self, SubductionError> {
        Self::build(subgroup, None)
    }

    /// Validate one explicit `(setting, child_shift)` convention for a record.
    ///
    /// This runs exactly the construction and validation
    /// [`Self::from_isotropy_subgroup`] performs, but with the convention handed
    /// in instead of read from the frozen table (and instead of the candidate
    /// search).  The task-9 provenance tools use it to test oracle-derived
    /// conventions against the engine's own parent frame *before* freezing them,
    /// so a setting is never frozen on the strength of a re-implementation of
    /// this check.  A rejected convention reports the same errors a frozen row
    /// would, and the returned embedding is never cached.
    pub fn probe_embedding(
        subgroup: &IsotropySubgroup,
        setting: Mat3I,
        setting_denominator: i32,
        child_shift: [i32; 4],
    ) -> Result<Self, SubductionError> {
        if setting_denominator <= 0 {
            return Err(SubductionError::InvalidGrid {
                denominator: i128::from(setting_denominator),
            });
        }
        Self::build(subgroup, Some((setting, setting_denominator, child_shift)))
    }

    fn build(
        subgroup: &IsotropySubgroup,
        explicit: Option<(Mat3I, i32, [i32; 4])>,
    ) -> Result<Self, SubductionError> {
        let parent_sg = subgroup.parent_sg;
        if parent_sg == 0 || parent_sg > 230 {
            return Err(SubductionError::InvalidSpaceGroup { sg: parent_sg });
        }
        let stored = ISOTROPY_SUBGROUPS.get(subgroup.ordinal).ok_or(
            SubductionError::IsotropyOrdinalOutOfRange {
                ordinal: subgroup.ordinal,
                len: ISOTROPY_SUBGROUPS.len(),
            },
        )?;
        if !same_record(stored, &subgroup.record) {
            return Err(SubductionError::StaleIsotropyRecord {
                ordinal: subgroup.ordinal,
            });
        }
        let irrep = query::irreps_of(parent_sg)
            .iter()
            .find(|record| record.ml == subgroup.irrep_ml)
            .ok_or(SubductionError::UnknownParentIrrep {
                sg: parent_sg,
                ml: subgroup.irrep_ml,
            })?;
        if irrep.spinor
            || irrep.k_vector() != subgroup.k
            || !irrep
                .subgroups()
                .iter()
                .any(|record| std::ptr::eq(record, stored))
        {
            return Err(SubductionError::StaleIsotropyRecord {
                ordinal: subgroup.ordinal,
            });
        }
        let subgroup_sg = u8::try_from(subgroup.record.sg).map_err(|_| {
            SubductionError::InvalidSubgroupNumber {
                sg: subgroup.record.sg,
            }
        })?;
        if subgroup_sg == 0 {
            return Err(SubductionError::InvalidSubgroupNumber {
                sg: subgroup.record.sg,
            });
        }

        let parent_primitive = exact_primitive_basis(parent_sg)?;
        let subgroup_primitive = exact_primitive_basis(subgroup_sg)?;
        let parent_lattice = Lattice::new(parent_primitive)?;
        let stored_basis = Mat3R::from_ints(subgroup.record.basis);
        let basis_conventional = stored_basis.checked_mul(&parent_primitive)?;
        let subgroup_lattice = Lattice::new(basis_conventional)?;
        let origin = exact_origin(&subgroup.record.origin, &parent_primitive)?;
        let parent_operations =
            reduce_operations(&strict_sg_hall_ops(parent_sg)?.operations, &parent_lattice)?;
        let frozen = match explicit {
            Some(pair) => Some(pair),
            None => frozen_setting_for(subgroup.ordinal, parent_sg, subgroup_sg)?
                .map(|entry| (entry.3, entry.4, entry.5)),
        };
        let shift = match frozen {
            Some((_, _, delta)) => exact_origin(&delta, &Mat3R::identity())?,
            None => Vec3R::zero(),
        };
        let subgroup_operations =
            shift_operations(&strict_sg_hall_ops(subgroup_sg)?.operations, &shift)?;
        let expected = distinct_rotations(&subgroup_operations);
        if expected == 0 {
            return Err(SubductionError::NoValidEmbedding {
                subgroup_sg,
                candidates: 0,
            });
        }
        let to_conventional = subgroup_primitive.inverse()?;
        let transform_for =
            |setting: Mat3I, denominator: i32| -> Result<SeitzTransform, SubductionError> {
            let setting_matrix = rational_setting(setting, denominator)?;
            let basis = to_conventional
                .checked_mul(&setting_matrix.inverse()?.checked_mul(&basis_conventional)?)?;
            Ok(SeitzTransform::new(basis.transpose(), origin))
        };

        let (setting, setting_denominator, transform, operations, representatives, candidate_count) =
            match frozen {
            // A recorded setting is a convention, not a search result: it is
            // validated like any other candidate and a failure is reported
            // instead of silently falling back to a different setting.
            Some((setting, denominator, _)) => {
                let transform = transform_for(setting, denominator)?;
                match validate_candidate(
                    &subgroup_operations,
                    &parent_operations,
                    &transform,
                    &parent_lattice,
                    &subgroup_lattice,
                    expected,
                )? {
                    Some((operations, representatives)) => {
                        (setting, denominator, transform, operations, representatives, 1)
                    }
                    None => {
                        return Err(SubductionError::FrozenEmbeddingRejected {
                            subgroup_sg,
                            setting,
                            denominator,
                        });
                    }
                }
            }
            None => {
                let candidates = signed_permutations();
                let mut accepted: Vec<(Mat3I, SeitzTransform, Vec<ExactSeitz>, Vec<ExactSeitz>)> =
                    Vec::new();
                for setting in &candidates {
                    let transform = transform_for(*setting, 1)?;
                    if let Some((operations, representatives)) = validate_candidate(
                        &subgroup_operations,
                        &parent_operations,
                        &transform,
                        &parent_lattice,
                        &subgroup_lattice,
                        expected,
                    )? {
                        accepted.push((*setting, transform, operations, representatives));
                    }
                }
                match accepted.len() {
                    0 => {
                        return Err(SubductionError::NoValidEmbedding {
                            subgroup_sg,
                            candidates: candidates.len(),
                        });
                    }
                    1 => {
                        let (setting, transform, operations, representatives) = accepted.remove(0);
                        (setting, 1, transform, operations, representatives, 1)
                    }
                    count => {
                        return Err(SubductionError::AmbiguousEmbedding {
                            subgroup_sg,
                            candidates: count,
                        });
                    }
                }
            }
        };
        Ok(Self {
            parent_sg,
            subgroup_sg,
            ordinal: subgroup.ordinal,
            setting,
            setting_denominator,
            child_shift: shift,
            transform,
            parent_lattice,
            subgroup_lattice,
            operations,
            representatives,
            candidates: candidate_count,
        })
    }

    /// Parent space group number.
    pub const fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// Subgroup space group number.
    pub const fn subgroup_sg(&self) -> u8 {
        self.subgroup_sg
    }

    /// Isotropy record ordinal this embedding was built from.
    pub const fn ordinal(&self) -> usize {
        self.ordinal
    }

    /// The setting transform `U` that produced the accepted candidate.
    ///
    /// The convention is the exact rational matrix `setting() / denominator()`:
    /// most records are integral (denominator one), but the stored isotropy
    /// basis of a few monoclinic records is consistent with the official cell
    /// only through a rational change of basis, and both parts are frozen so the
    /// engine reproduces the official cell exactly.
    pub const fn setting(&self) -> Mat3I {
        self.setting
    }

    /// Denominator of [`Self::setting`]; always positive.
    pub const fn setting_denominator(&self) -> i32 {
        self.setting_denominator
    }

    /// The subgroup-frame origin shift applied to the shipped operations.
    ///
    /// Non-zero when the isotropy record uses a different ITA origin choice than
    /// the subgroup's Hall setting (#126 is the recorded example).  Mapping an
    /// operation *back* into the subgroup frame therefore has to undo it before
    /// the shipped rows can be looked up.
    pub const fn child_shift(&self) -> &Vec3R {
        &self.child_shift
    }

    /// `x_G = T x_H + o`.
    pub const fn transform(&self) -> &SeitzTransform {
        &self.transform
    }

    /// Parent translation lattice `L_G`.
    pub const fn parent_lattice(&self) -> &Lattice {
        &self.parent_lattice
    }

    /// Subgroup translation lattice `L_H` (parent conventional frame).
    pub const fn subgroup_lattice(&self) -> &Lattice {
        &self.subgroup_lattice
    }

    /// Subgroup operations mapped into the parent frame (reduced mod `L_G`).
    pub fn operations(&self) -> &[ExactSeitz] {
        &self.operations
    }

    /// Coset representatives of the subgroup's translation subgroup, reduced
    /// modulo `L_H`; their number is the subgroup's point group order.
    pub fn representatives(&self) -> &[ExactSeitz] {
        &self.representatives
    }

    /// How many setting candidates were consistent with the parent group.
    pub const fn candidate_count(&self) -> usize {
        self.candidates
    }
}

/// Compare the fields the embedding relies on, so a hand-built record cannot
/// pass as the stored one.
fn same_record(left: &IsotropyRecord, right: &IsotropyRecord) -> bool {
    left.sg == right.sg
        && left.basis == right.basis
        && left.origin == right.origin
        && left.direction_label == right.direction_label
        && left.direction == right.direction
        && left.domains == right.domains
        && left.arms == right.arms
}

/// The stored primitive basis as exact rationals.
fn exact_primitive_basis(sg: u8) -> Result<Mat3R, SubductionError> {
    let basis =
        parent_primitive_basis(sg).map_err(|_| SubductionError::InvalidSpaceGroup { sg })?;
    let mut rows = [[Rat::ZERO; 3]; 3];
    for (row, values) in basis.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            rows[row][column] = Rat::from_grid(*value, SOURCE_TRANSLATION_GRID)?;
        }
    }
    Ok(Mat3R::new(rows))
}

/// The stored origin `(x, y, z, d)` as a point in the parent conventional frame.
/// `U = numerator / denominator` as an exact rational matrix.
///
/// The frozen conventions are stored as an integer numerator plus one shared
/// positive denominator because a static table cannot build [`Rat`] values with
/// [`Rat::new`].  A denominator of one is the common case; the rational form is
/// needed only where the stored isotropy basis reaches the official subgroup
/// cell through a fractional change of basis.
fn rational_setting(numerator: Mat3I, denominator: i32) -> Result<Mat3R, SubductionError> {
    if denominator <= 0 {
        return Err(SubductionError::InvalidGrid {
            denominator: i128::from(denominator),
        });
    }
    if denominator == 1 {
        return Ok(Mat3R::from_ints(numerator));
    }
    let denominator = i128::from(denominator);
    let mut rows = [[Rat::ZERO; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            rows[row][column] = Rat::new(i128::from(numerator[row][column]), denominator)?;
        }
    }
    Ok(Mat3R::new(rows))
}

fn exact_origin(origin: &[i32; 4], parent_primitive: &Mat3R) -> Result<Vec3R, SubductionError> {
    if origin[3] <= 0 {
        return Err(SubductionError::InvalidGrid {
            denominator: i128::from(origin[3]),
        });
    }
    let denominator = i128::from(origin[3]);
    let primitive = Vec3R::new([
        Rat::new(i128::from(origin[0]), denominator)?,
        Rat::new(i128::from(origin[1]), denominator)?,
        Rat::new(i128::from(origin[2]), denominator)?,
    ]);
    parent_primitive.transpose().checked_mul_vector(&primitive)
}

/// Re-express subgroup operations after an origin shift `delta` in the
/// subgroup's own frame: `t' = t + delta - R delta`.
///
/// The isotropy tables and the shipped Hall setting of a subgroup can use
/// different ITA origin choices; the shift is frozen per isotropy ordinal in
/// `settings_data::FROZEN_EMBEDDING_SETTINGS` and validated like everything
/// else.
fn shift_operations(
    operations: &[ExactSeitz],
    shift: &Vec3R,
) -> Result<Vec<ExactSeitz>, SubductionError> {
    let mut out = Vec::with_capacity(operations.len());
    for operation in operations {
        let rotated = Mat3R::from_ints(operation.rotation()).checked_mul_vector(shift)?;
        let translation = operation
            .translation()
            .checked_add(shift)?
            .checked_sub(&rotated)?;
        out.push(ExactSeitz::new(operation.rotation(), translation));
    }
    Ok(out)
}

/// Canonical representatives of a mapped operation set.
fn reduce_operations(
    operations: &[ExactSeitz],
    lattice: &Lattice,
) -> Result<Vec<ExactSeitz>, SubductionError> {
    let mut out = Vec::with_capacity(operations.len());
    for operation in operations {
        out.push(operation.reduce(lattice)?);
    }
    Ok(out)
}

/// Number of distinct rotations in an operation set: the point group order for
/// a conventional cell, centring translations included.
fn distinct_rotations(operations: &[ExactSeitz]) -> usize {
    let mut rotations: Vec<Mat3I> = Vec::new();
    for operation in operations {
        if !rotations.contains(&operation.rotation()) {
            rotations.push(operation.rotation());
        }
    }
    rotations.len()
}

/// A validated candidate: mapped operations and their coset representatives.
type CandidateOperations = (Vec<ExactSeitz>, Vec<ExactSeitz>);

/// Map the subgroup's own operations into the parent frame and check that they
/// form the subgroup there.  `None` means "this setting does not work".
fn validate_candidate(
    subgroup_operations: &[ExactSeitz],
    parent_operations: &[ExactSeitz],
    transform: &SeitzTransform,
    parent_lattice: &Lattice,
    subgroup_lattice: &Lattice,
    expected: usize,
) -> Result<Option<CandidateOperations>, SubductionError> {
    let mut operations = Vec::with_capacity(subgroup_operations.len());
    // The mapped operations are kept **before** the parent-lattice reduction as
    // well: reducing mod `L_G` can remove a vector in `L_G \ L_H`, and
    // deduplicating those already-reduced representatives mod `L_H` would then
    // report a class that the subgroup does not have (SG 139 `M1-` P1 -> #126:
    // the inversion `(-I | 5/2,5/2,5/2)` reduces to translation zero mod `L_G`
    // because `(1/2,1/2,1/2)` is an I-centring vector, but mod `L_H = Z^3` its
    // representative is `(1/2,1/2,1/2)`, not zero).
    let mut mapped_operations = Vec::with_capacity(subgroup_operations.len());
    for operation in subgroup_operations {
        let mapped = match transform.map_operation(operation) {
            Ok(mapped) => mapped,
            Err(SubductionError::NonIntegralRotationImage { .. }) => return Ok(None),
            Err(error) => return Err(error),
        };
        let reduced = mapped.reduce(parent_lattice)?;
        if !parent_operations.contains(&reduced) {
            return Ok(None);
        }
        operations.push(reduced);
        mapped_operations.push(mapped);
    }
    // Representatives are the subgroup's cosets modulo `L_H`, derived from the
    // original mapped operations so `L_G \ L_H` translations survive; the
    // parent-reduced `operations` above stay the membership view exposed by
    // [`SubgroupEmbedding::operations`].
    let representatives = subgroup_lattice.deduplicate(&mapped_operations)?;
    if representatives.len() != expected {
        return Ok(None);
    }
    for left in &representatives {
        let inverse = left.inverse()?.reduce(subgroup_lattice)?;
        if !representatives.contains(&inverse) {
            return Ok(None);
        }
        for right in &representatives {
            let product = left.compose(right)?.reduce(subgroup_lattice)?;
            if !representatives.contains(&product) {
                return Ok(None);
            }
        }
    }
    Ok(Some((operations, representatives)))
}

// ── First complete Gamma decomposition ───────────────────────────────────────

/// Tolerance for character inner products, multiplicities and reconstruction.
///
/// The subduced characters are algebraic numbers built from the shipped tables;
/// this only absorbs the floating-point summation error of the inner products.
const SUBDUCTION_TOLERANCE: f64 = 1e-7;

/// Which complex constituent of the subgroup's character data a target is.
///
/// A compound record is a *sum* of complex irreps, so it cannot be used as a
/// single target.  The component identifies the constituent by its stable CIR
/// source (`irnumber`), never by a name synthesized from the row label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubductionComponent {
    /// The record is an ordinary (single complex) irrep.
    Ordinary,
    /// Constituent of a `DistinctComponentSum` row.
    Constituent { index: u8, irnumber: u32 },
    /// The CIR seed of a `ConjugateRealification` row.
    RealificationSeed { irnumber: u32 },
    /// The conjugate of that seed: same source identity, conjugated characters.
    RealificationConjugate { irnumber: u32 },
}

/// One target irrep of a subduced representation.
#[derive(Debug, Clone, PartialEq)]
pub struct SubductionTarget {
    /// Subgroup space group number.
    pub sg: u8,
    /// Stable source label: the CIR label for a constituent, otherwise the
    /// Miller-Love label of the row.
    pub ml: &'static str,
    /// Bradley-Cracknell label, for display only.
    pub bc: &'static str,
    /// The physical row this constituent came from.
    pub row_ml: &'static str,
    /// Complex dimension of this irreducible constituent.
    pub dimension: u8,
    /// Which complex constituent of that row this target is.
    pub component: SubductionComponent,
    /// Multiplicity in the subduced representation.
    pub multiplicity: u32,
}

/// The decomposition of one parent irrep on one isotropy subgroup.
///
/// Produced by [`subduce_irrep`]; the embedded subgroup is part of the result
/// because the same subgroup number reached through a different direction is a
/// different decomposition.
#[derive(Debug, Clone)]
pub struct IrrepSubduction {
    parent_sg: u8,
    parent_ml: &'static str,
    parent_bc: &'static str,
    parent_dimension: u8,
    subgroup_sg: u8,
    ordinal: usize,
    setting: Mat3I,
    folded_k: [Rat; 3],
    targets: Vec<SubductionTarget>,
    parent_characters: Vec<Complex64>,
    reconstructed: Vec<Complex64>,
    tolerance: f64,
}

impl IrrepSubduction {
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

    /// Complex dimension of the restricted parent character space.
    pub const fn parent_dimension(&self) -> u8 {
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

    /// The probe's wave vector folded into the subgroup frame, `T^T k_G`.
    ///
    /// Exact: the subgroup k-block this result decomposes is the one equivalent
    /// to this vector modulo the subgroup's reciprocal lattice.
    pub const fn folded_k(&self) -> [Rat; 3] {
        self.folded_k
    }

    /// Non-zero terms, in the subgroup's irrep order.
    pub fn targets(&self) -> &[SubductionTarget] {
        &self.targets
    }

    /// Multiplicity of one target label, or `0` when it does not appear.
    pub fn multiplicity(&self, ml: &str) -> u32 {
        self.targets
            .iter()
            .find(|target| target.ml == ml)
            .map_or(0, |target| target.multiplicity)
    }

    /// Subduced characters on the coset representatives, and their
    /// reconstruction from the reported multiplicities: equal entries are the
    /// witness that the decomposition is complete.
    pub fn reconstruction(&self) -> (&[Complex64], &[Complex64]) {
        (&self.parent_characters, &self.reconstructed)
    }

    /// Tolerance used for the multiplicity and reconstruction checks.
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

/// Decompose one parent irrep on the subgroup selected by an isotropy record.
///
/// `subgroup` carries both the condensing context (parent, irrep, direction,
/// ordinal) and the embedding; `probe_ml` is the **subduced** parent irrep.
/// This convenience entry point accepts Gamma scalar probes, including compound
/// rows. Use [`subduce_irrep_with_embedding`] for supported non-Gamma probes.
pub fn subduce_irrep(
    subgroup: &IsotropySubgroup,
    probe_ml: &str,
) -> Result<IrrepSubduction, SubductionError> {
    let parent_sg = subgroup.parent_sg;
    if subgroup.k.numerators != [0, 0, 0] {
        return Err(SubductionError::ProbeNotAtGamma {
            sg: parent_sg,
            ml: subgroup.irrep_ml.to_string(),
            k: subgroup.k,
        });
    }
    let probe = match query::irreps_of(parent_sg)
        .iter()
        .find(|record| record.ml == probe_ml && is_gamma(record))
    {
        Some(record) if record.spinor => {
            return Err(SubductionError::UnsupportedCharacterSpace {
                sg: parent_sg,
                ml: record.ml.to_string(),
            });
        }
        Some(record) => record,
        None => {
            // Distinguish "this label exists at another k" from "no such label".
            return match query::irreps_of(parent_sg)
                .iter()
                .find(|record| record.ml == probe_ml)
            {
                Some(other) => Err(SubductionError::ProbeNotAtGamma {
                    sg: parent_sg,
                    ml: other.ml.to_string(),
                    k: KVector {
                        numerators: [other.kx, other.ky, other.kz],
                        denominator: other.kd,
                    },
                }),
                None => Err(SubductionError::ProbeIrrepNotFound {
                    sg: parent_sg,
                    ml: probe_ml.to_string(),
                }),
            };
        }
    };
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup)?;
    subduce_irrep_with_embedding(subgroup, &embedding, probe)
}

/// Check that `subgroup`, `embedding` and `probe` still describe one context.
///
/// [`SubgroupEmbedding`] is expensive to build (it searches setting candidates),
/// so it is cached and reused across the probes of one isotropy record.  The
/// cache is only sound while the three inputs agree; this revalidates them
/// without rebuilding the embedding:
///
/// * `IsotropySubgroup` has public fields, so its ordinal, its record and its
///   irrep context are rechecked against the generated table.  The record must
///   be the stored one **and** must be listed by the parent irrep that
///   `irrep_ml` names, which is what ties parent, irrep and record together; an
///   ordinal-only guard would accept a record whose fields were mutated.
/// * the cached embedding must have been built from that same record: same
///   parent, same ordinal, same subgroup number.  Its setting, child shift and
///   operation lists all come from that record, so a mismatch would fold one
///   record's k in another record's frame.
/// * the probe must be the parent table's *own* record, by identity rather than
///   by label: a label is not unique across space groups and a copy is a
///   different object.
fn validate_subduction_context(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<(), SubductionError> {
    let parent_sg = subgroup.parent_sg;
    if parent_sg == 0 || parent_sg > 230 {
        return Err(SubductionError::InvalidSpaceGroup { sg: parent_sg });
    }
    let ordinal = subgroup.ordinal;
    let stored =
        ISOTROPY_SUBGROUPS
            .get(ordinal)
            .ok_or(SubductionError::IsotropyOrdinalOutOfRange {
                ordinal,
                len: ISOTROPY_SUBGROUPS.len(),
            })?;
    if !same_record(stored, &subgroup.record) {
        return Err(SubductionError::StaleIsotropyRecord { ordinal });
    }
    let subgroup_sg = u8::try_from(subgroup.record.sg)
        .ok()
        .filter(|sg| *sg != 0)
        .ok_or(SubductionError::InvalidSubgroupNumber {
            sg: subgroup.record.sg,
        })?;

    // The record must belong to the irrep the context names.  `subgroups()`
    // panics on spinor irreps, so that case is rejected before the lookup.
    let irrep = query::irreps_of(parent_sg)
        .iter()
        .find(|record| record.ml == subgroup.irrep_ml)
        .ok_or(SubductionError::UnknownParentIrrep {
            sg: parent_sg,
            ml: subgroup.irrep_ml,
        })?;
    if irrep.spinor
        || irrep.k_vector() != subgroup.k
        || !irrep
            .subgroups()
            .iter()
            .any(|record| std::ptr::eq(record, stored))
    {
        return Err(SubductionError::StaleIsotropyRecord { ordinal });
    }

    if !query::irreps_of(parent_sg)
        .iter()
        .any(|record| std::ptr::eq(record, probe))
    {
        return Err(SubductionError::ForeignProbe {
            sg: parent_sg,
            ml: probe.ml,
        });
    }

    if embedding.parent_sg() != parent_sg
        || embedding.ordinal() != ordinal
        || embedding.subgroup_sg() != subgroup_sg
    {
        return Err(SubductionError::EmbeddingContextMismatch {
            parent_sg,
            ordinal,
            subgroup_sg,
            embedding_parent_sg: embedding.parent_sg(),
            embedding_ordinal: embedding.ordinal(),
            embedding_subgroup_sg: embedding.subgroup_sg(),
        });
    }
    Ok(())
}

/// [`subduce_irrep`] with an embedding that was already built and validated.
///
/// Scanning many probes of the same subgroup should not rebuild the embedding
/// (and re-run its candidate validation) for every probe.
///
/// Reuse is only sound while `subgroup`, `embedding` and `probe` still name one
/// context, so every call revalidates the cheap identities first (see
/// [`validate_subduction_context`]): the public `IsotropySubgroup` fields are
/// rechecked against the generated table, the cached embedding must have been
/// built from that exact record, and the probe must be the parent's own table
/// record.  A stale record, a foreign probe or a mismatched embedding is an
/// error, never a silently reinterpreted decomposition.
///
/// Wave vectors are folded exactly: `k_H = T^T k_G` with checked rational
/// arithmetic, then matched against the subgroup's stored k-points **modulo the
/// subgroup's reciprocal lattice**, which keeps the centring extinctions.  Only
/// single-arm parent stars are supported (the little group must be the whole
/// point group); everything else is an explicit error, and a folded k with no
/// stored data is reported as missing rather than truncated or rounded onto the
/// nearest point.
pub fn subduce_irrep_with_embedding(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<IrrepSubduction, SubductionError> {
    validate_subduction_context(subgroup, embedding, probe)?;
    let parent_sg = subgroup.parent_sg;
    let parent_lattice = *embedding.parent_lattice();

    // Fold the probe's wave vector into the subgroup frame.
    let wave_vector = inline_k_vector(probe)?;
    let folded = fold_wave_vector(embedding.transform(), &wave_vector)?;
    let parent_reciprocal = parent_lattice.reciprocal()?;
    // Single arm: every parent operation fixes k_G.  A larger star folds onto
    // several subgroup stars and needs the full-star adapter (task 8).
    for operation in &strict_sg_hall_ops(parent_sg)?.operations {
        if !parent_reciprocal.preserves(operation.rotation(), &wave_vector)? {
            return Err(SubductionError::UnsupportedMultiArmStar {
                sg: parent_sg,
                k: [folded.get(0), folded.get(1), folded.get(2)],
            });
        }
    }

    // The subgroup block at the folded k, matched modulo its reciprocal lattice.
    let child_cell = Lattice::new(exact_primitive_basis(embedding.subgroup_sg())?)?;
    let child_reciprocal = child_cell.reciprocal()?;
    let mut block: Vec<&'static IrrepRecord> = Vec::new();
    for record in query::irreps_of(embedding.subgroup_sg()) {
        if record.spinor {
            continue;
        }
        let stored = inline_k_vector(record)?;
        if child_reciprocal.contains(&folded.checked_sub(&stored)?)? {
            block.push(record);
        }
    }
    if block.is_empty() {
        return Err(SubductionError::MissingIrrepData {
            sg: embedding.subgroup_sg(),
            k: [folded.get(0), folded.get(1), folded.get(2)],
        });
    }

    // Active coset representatives: those in the subgroup's little group of the
    // folded k.  Characters of a selected-arm row are only meaningful there, so
    // they are looked up on this set (the row itself is indexed by the complete
    // operation universe).
    // Parent-side characters and the active representative list do not depend on
    // the subgroup-frame convention: rotations are unaffected by an origin shift
    // and the little-group test only uses them.
    let mut parent_characters = Vec::new();
    let mut active: Vec<ExactSeitz> = Vec::new();
    for operation in embedding.representatives() {
        // `unmap_operation` inverts internally: hand it the forward transform.
        // The image keeps its translation exactly; reducing here would discard
        // child lattice vectors and with them the Bloch phase the shipped row
        // stores for the representative (see `character_of`, which pairs modulo
        // the lattice and corrects the phase).
        let child = embedding.transform().unmap_operation(operation)?;
        if !child_reciprocal.preserves(child.rotation(), &folded)? {
            continue;
        }
        let index = active.len();
        parent_characters.push(parent_character_of(
            probe,
            &parent_lattice,
            operation,
            &wave_vector,
            index,
        )?);
        active.push(child);
    }
    if parent_characters.is_empty() {
        return Err(SubductionError::MissingIrrepData {
            sg: embedding.subgroup_sg(),
            k: [folded.get(0), folded.get(1), folded.get(2)],
        });
    }

    // The embedding's operation list was built from the subgroup's shipped Hall
    // operations by applying `child_shift`; the shipped rows live in the frame
    // *before* that shift, so undo it exactly and deterministically.  There is
    // no second origin to choose between: `child_shift` is frozen per embedding
    // and validated when the embedding is built.
    let pulled_back = shift_operations(&active, &embedding.child_shift().checked_neg()?)?;
    decompose_active_block(
        subgroup,
        embedding,
        probe,
        block.as_slice(),
        &folded,
        &parent_characters,
        &pulled_back,
        &child_cell,
    )
}

/// Multiplicities and reconstruction witnesses of one solved character block.
struct SolvedCharacterBlock {
    /// `chi(E)` of the block: the dimension of the subduced space.
    dimension: u8,
    /// Non-zero target terms, in the child table's order.
    targets: Vec<SubductionTarget>,
    /// The targets re-summed on the same operations the block was evaluated on.
    reconstructed: Vec<Complex64>,
}

/// Decompose the active block with one concrete subgroup-frame reading.
#[allow(clippy::too_many_arguments)]
fn decompose_active_block(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
    block: &[&'static IrrepRecord],
    folded: &Vec3R,
    parent_characters: &[Complex64],
    pulled_back: &[ExactSeitz],
    child_cell: &Lattice,
) -> Result<IrrepSubduction, SubductionError> {
    let solved = solve_character_block(
        embedding.subgroup_sg(),
        block,
        parent_characters,
        pulled_back,
        child_cell,
        folded,
    )?;
    Ok(IrrepSubduction {
        parent_sg: subgroup.parent_sg,
        parent_ml: probe.ml,
        parent_bc: probe.bc,
        parent_dimension: solved.dimension,
        subgroup_sg: embedding.subgroup_sg(),
        ordinal: embedding.ordinal(),
        setting: embedding.setting(),
        folded_k: [folded.get(0), folded.get(1), folded.get(2)],
        targets: solved.targets,
        parent_characters: parent_characters.to_vec(),
        reconstructed: solved.reconstructed,
        tolerance: SUBDUCTION_TOLERANCE,
    })
}

/// Decompose an aligned character block into the block's complex irreps.
///
/// `parent_characters[i]` corresponds to `pulled_back[i]` in the child frame.
/// Both the identity dimension and the child characters use this aligned list,
/// never the embedding's unfiltered representative list.
fn solve_character_block(
    subgroup_sg: u8,
    block: &[&'static IrrepRecord],
    parent_characters: &[Complex64],
    pulled_back: &[ExactSeitz],
    child_cell: &Lattice,
    folded: &Vec3R,
) -> Result<SolvedCharacterBlock, SubductionError> {
    let targets = complex_targets(subgroup_sg, block, pulled_back, child_cell, folded)?;
    solve_prepared_character_block(subgroup_sg, targets, parent_characters, pulled_back)
}

/// The shared decomposition core: Gram identity, integral multiplicities,
/// dimension sum and per-operation reconstruction of one prepared target list.
///
/// `targets` are already evaluated as complex-irreducible rows on the same
/// aligned `pulled_back` operation list.  [`solve_character_block`] prepares
/// them from one child `k` block; the full-star adapter prepares its own list,
/// because its child target rows can live on another arm of the child star and
/// then need exact transport before they are comparable.
fn solve_prepared_character_block(
    subgroup_sg: u8,
    targets: Vec<ComplexTarget>,
    parent_characters: &[Complex64],
    pulled_back: &[ExactSeitz],
) -> Result<SolvedCharacterBlock, SubductionError> {
    let parent_dimension = complex_dimension(parent_characters, pulled_back)?;
    let count = parent_characters.len();
    let scale = 1.0 / count as f64;

    // Complex-irreducibility of every target over the active set, and pairwise
    // orthogonality: the Gram matrix must be the identity on the k-block.  A
    // compound row fed in as one target fails here (it has norm 2), which is the
    // whole point of expanding constituents first.
    for (index, first) in targets.iter().enumerate() {
        let first_values = &first.values;
        let norm = (first_values
            .iter()
            .map(|value| value.norm_sqr())
            .sum::<f64>()
            * scale)
            .sqrt();
        if (norm - 1.0).abs() > SUBDUCTION_TOLERANCE {
            return Err(SubductionError::TargetNotIrreducible { ml: first.ml, norm });
        }
        for second in targets.iter().skip(index + 1) {
            let second_values = &second.values;
            let inner: Complex64 = first_values
                .iter()
                .zip(second_values)
                .map(|(left, right)| left * right.conj())
                .sum::<Complex64>()
                * scale;
            if inner.norm() > SUBDUCTION_TOLERANCE {
                return Err(SubductionError::TargetRowsNotOrthogonal {
                    first: first.ml,
                    second: second.ml,
                    value: inner,
                });
            }
        }
    }

    let mut dimension_sum = 0i64;
    let mut reconstructed = vec![Complex64::new(0.0, 0.0); count];
    let mut reported: Vec<SubductionTarget> = Vec::new();
    let mut seen: Vec<(SubductionComponent, &'static str, u32)> = Vec::new();
    for target in targets.iter() {
        let inner: Complex64 = parent_characters
            .iter()
            .zip(&target.values)
            .map(|(parent, value)| parent * value.conj())
            .sum::<Complex64>()
            * scale;
        let multiplicity = integral_multiplicity(target.ml, inner)?;
        if multiplicity == 0 {
            continue;
        }
        let key = (target.component, target.ml, target.irnumber);
        if seen.contains(&key) {
            return Err(SubductionError::DuplicateTargetSource {
                ml: target.ml,
                irnumber: target.irnumber,
            });
        }
        seen.push(key);
        dimension_sum += i64::from(multiplicity) * i64::from(target.dimension);
        for (slot, value) in reconstructed.iter_mut().zip(&target.values) {
            *slot += value * f64::from(multiplicity);
        }
        reported.push(SubductionTarget {
            sg: subgroup_sg,
            ml: target.ml,
            bc: target.bc,
            row_ml: target.row_ml,
            dimension: target.dimension,
            component: target.component,
            multiplicity,
        });
    }
    if dimension_sum != i64::from(parent_dimension) {
        return Err(SubductionError::DimensionSumMismatch {
            expected: u32::from(parent_dimension),
            found: dimension_sum,
        });
    }
    for (index, (found, expected)) in reconstructed.iter().zip(parent_characters).enumerate() {
        if (found - expected).norm() > SUBDUCTION_TOLERANCE {
            return Err(SubductionError::CharacterMismatch {
                index,
                found: *found,
                expected: *expected,
            });
        }
    }
    Ok(SolvedCharacterBlock {
        dimension: parent_dimension,
        targets: reported,
        reconstructed,
    })
}

/// An irrep record's wave vector as exact rationals.
fn inline_k_vector(record: &IrrepRecord) -> Result<Vec3R, SubductionError> {
    if record.kd <= 0 {
        return Err(SubductionError::InvalidGrid {
            denominator: i128::from(record.kd),
        });
    }
    exact_k_vector(KVector {
        numerators: [record.kx, record.ky, record.kz],
        denominator: record.kd,
    })
}

/// A wave vector as exact rationals; the narrow `i8` fields are widened, never
/// truncated.
fn exact_k_vector(k: KVector) -> Result<Vec3R, SubductionError> {
    let denominator = i128::from(k.denominator);
    if denominator <= 0 {
        return Err(SubductionError::InvalidGrid { denominator });
    }
    Ok(Vec3R::new([
        Rat::new(i128::from(k.numerators[0]), denominator)?,
        Rat::new(i128::from(k.numerators[1]), denominator)?,
        Rat::new(i128::from(k.numerators[2]), denominator)?,
    ]))
}

/// `k_H = T^T k_G`, exactly.
pub fn fold_wave_vector(
    transform: &SeitzTransform,
    wave_vector: &Vec3R,
) -> Result<Vec3R, SubductionError> {
    transform
        .matrix()
        .transpose()
        .checked_mul_vector(wave_vector)
}

/// The complex dimension of the parent representation: its character on the
/// identity operation.  This is *not* `IrrepRecord::dim` for a compound probe,
/// where the physical row is a realification of a complex irrep.
fn complex_dimension(
    characters: &[Complex64],
    representatives: &[ExactSeitz],
) -> Result<u8, SubductionError> {
    let identity_rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let index = representatives
        .iter()
        .position(|operation| operation.rotation() == identity_rotation)
        .ok_or(SubductionError::TargetNotIrreducible {
            ml: "identity operation",
            norm: 0.0,
        })?;
    let value = characters[index].re;
    let rounded = value.round();
    if (value - rounded).abs() > SUBDUCTION_TOLERANCE
        || !(0.0..=f64::from(u8::MAX)).contains(&rounded)
    {
        return Err(SubductionError::NonIntegralMultiplicity {
            ml: "identity character",
            value: Complex64::new(value, 0.0),
        });
    }
    // Exact by the checks above.
    Ok(rounded as u8)
}

/// The character row of an irrep, whichever typed space it lives in.
enum ProbeRow {
    Ordinary(CharacterRow),
    Compound(Box<CompoundSelectedArmCharacter>),
}

fn probe_row(record: &'static IrrepRecord) -> Result<ProbeRow, SubductionError> {
    match record.ordinary_scalar_selected_arm_block_trace() {
        Ok(row) => Ok(ProbeRow::Ordinary(row)),
        Err(crate::irrep::types::CharacterViewError::NotApplicable) => record
            .compound_selected_arm_view()
            .map(|view| ProbeRow::Compound(Box::new(view)))
            .map_err(|_| SubductionError::UnsupportedCharacterSpace {
                sg: record.sg,
                ml: record.ml.to_string(),
            }),
        Err(_) => Err(SubductionError::UnsupportedCharacterSpace {
            sg: record.sg,
            ml: record.ml.to_string(),
        }),
    }
}

/// The probe's character on one parent-frame operation.
///
/// Compound probes contribute the **sum** of their complex constituents
/// (`first + second`, or `seed + conjugate(seed)` for a realification), which is
/// exactly the character of the physical representation.
fn parent_character_of(
    record: &'static IrrepRecord,
    lattice: &Lattice,
    operation: &ExactSeitz,
    wave_vector: &Vec3R,
    index: usize,
) -> Result<Complex64, SubductionError> {
    match probe_row(record)? {
        ProbeRow::Ordinary(row) => {
            character_of(&row, record.ml, operation, lattice, wave_vector, index)
        }
        ProbeRow::Compound(view) => match view.as_ref() {
            CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } => {
                let left = character_of(
                    &first.row,
                    first.label,
                    operation,
                    lattice,
                    wave_vector,
                    index,
                )?;
                let right = character_of(
                    &second.row,
                    second.label,
                    operation,
                    lattice,
                    wave_vector,
                    index,
                )?;
                Ok(left + right)
            }
            CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => {
                let value = character_of(
                    &seed.row,
                    seed.label,
                    operation,
                    lattice,
                    wave_vector,
                    index,
                )?;
                Ok(Complex64::new(2.0 * value.re, 0.0))
            }
        },
    }
}

/// A complex irreducible constituent of the subgroup's character data.
struct ComplexTarget {
    ml: &'static str,
    bc: &'static str,
    row_ml: &'static str,
    irnumber: u32,
    dimension: u8,
    component: SubductionComponent,
    values: Vec<Complex64>,
}

/// Evaluate one row over the pulled-back operations.
fn evaluate(
    row: &CharacterRow,
    ml: &'static str,
    pulled_back: &[ExactSeitz],
    child_cell: &Lattice,
    wave_vector: &Vec3R,
) -> Result<Vec<Complex64>, SubductionError> {
    pulled_back
        .iter()
        .enumerate()
        .map(|(index, operation)| character_of(row, ml, operation, child_cell, wave_vector, index))
        .collect()
}

/// Expand the subgroup's Gamma rows into complex irreducible constituents.
///
/// There is deliberately no comparison against the stored block trace: that
/// trace is *assembled from these same constituents* by
/// `compound_selected_arm_view`, so gating on it was circular and proved
/// nothing.  Independent physical-character fixtures for compound rows live in
/// `tests/compound_subduction_regressions.rs`.
///
/// Ordinary rows are already single complex irreps.  A `DistinctComponentSum`
/// row contributes both stored CIR constituents; a `ConjugateRealification` row
/// contributes its CIR seed **and** the conjugate of that seed, because the
/// subduced representation is complex and the two multiplicities need not be
/// equal.  Feeding a compound row in as one target would violate the complex
/// orthogonality relation (its norm is 2), which is why it is expanded here.
fn complex_targets(
    subgroup_sg: u8,
    block: &[&'static IrrepRecord],
    pulled_back: &[ExactSeitz],
    child_cell: &Lattice,
    wave_vector: &Vec3R,
) -> Result<Vec<ComplexTarget>, SubductionError> {
    let mut out = Vec::new();
    for record in block {
        let record = *record;
        match record.ordinary_scalar_selected_arm_block_trace() {
            Ok(row) => {
                let dimension = u8::try_from(row.dimension()).map_err(|_| {
                    SubductionError::UnsupportedCharacterSpace {
                        sg: subgroup_sg,
                        ml: record.ml.to_string(),
                    }
                })?;
                out.push(ComplexTarget {
                    ml: record.ml,
                    bc: record.bc,
                    row_ml: record.ml,
                    irnumber: 0,
                    dimension,
                    component: SubductionComponent::Ordinary,
                    values: evaluate(&row, record.ml, pulled_back, child_cell, wave_vector)?,
                });
            }
            Err(crate::irrep::types::CharacterViewError::NotApplicable) => {
                let view = record.compound_selected_arm_view().map_err(|_| {
                    SubductionError::UnsupportedCharacterSpace {
                        sg: subgroup_sg,
                        ml: record.ml.to_string(),
                    }
                })?;
                match view {
                    CompoundSelectedArmCharacter::DistinctComponentSum {
                        first,
                        second,
                        block_trace: _,
                    } => {
                        let first_values = evaluate(
                            &first.row,
                            first.label,
                            pulled_back,
                            child_cell,
                            wave_vector,
                        )?;
                        let second_values = evaluate(
                            &second.row,
                            second.label,
                            pulled_back,
                            child_cell,
                            wave_vector,
                        )?;
                        for (index, (constituent, values)) in
                            [(first, first_values), (second, second_values)]
                                .into_iter()
                                .enumerate()
                        {
                            let dimension = u8::try_from(constituent.dimension).map_err(|_| {
                                SubductionError::UnsupportedCharacterSpace {
                                    sg: subgroup_sg,
                                    ml: constituent.label.to_string(),
                                }
                            })?;
                            out.push(ComplexTarget {
                                ml: constituent.label,
                                bc: record.bc,
                                row_ml: record.ml,
                                irnumber: constituent.irnumber,
                                dimension,
                                component: SubductionComponent::Constituent {
                                    index: index as u8,
                                    irnumber: constituent.irnumber,
                                },
                                values,
                            });
                        }
                    }
                    CompoundSelectedArmCharacter::ConjugateRealification {
                        seed,
                        block_trace: _,
                    } => {
                        let seed_values =
                            evaluate(&seed.row, seed.label, pulled_back, child_cell, wave_vector)?;
                        let conjugate: Vec<Complex64> =
                            seed_values.iter().map(|value| value.conj()).collect();
                        let dimension = u8::try_from(seed.dimension).map_err(|_| {
                            SubductionError::UnsupportedCharacterSpace {
                                sg: subgroup_sg,
                                ml: seed.label.to_string(),
                            }
                        })?;
                        out.push(ComplexTarget {
                            ml: seed.label,
                            bc: record.bc,
                            row_ml: record.ml,
                            irnumber: seed.irnumber,
                            dimension,
                            component: SubductionComponent::RealificationSeed {
                                irnumber: seed.irnumber,
                            },
                            values: seed_values,
                        });
                        out.push(ComplexTarget {
                            ml: seed.label,
                            bc: record.bc,
                            row_ml: record.ml,
                            irnumber: seed.irnumber,
                            dimension,
                            component: SubductionComponent::RealificationConjugate {
                                irnumber: seed.irnumber,
                            },
                            values: conjugate,
                        });
                    }
                }
            }
            Err(_) => {
                return Err(SubductionError::UnsupportedCharacterSpace {
                    sg: subgroup_sg,
                    ml: record.ml.to_string(),
                });
            }
        }
    }
    Ok(out)
}

/// Whether an irrep record sits at Gamma.
fn is_gamma(record: &IrrepRecord) -> bool {
    record.kx == 0 && record.ky == 0 && record.kz == 0
}

/// Bloch phase `exp(+2 pi i k . delta)`.
///
/// Sign convention, pinned by the shipped data (see
/// `shipped_rows_carry_the_bloch_phase_of_their_representatives`): a stored row
/// entry for the representative `t` and the same element represented by
/// `t + L` (with `L` a lattice vector) satisfy
/// `chi(t + L) = chi(t) * exp(+2 pi i k . L)`.  The phase is evaluated with
/// `f64` trigonometry from exact rational angles; the pairing itself stays
/// exact.
pub fn bloch_phase(wave_vector: &Vec3R, delta: &Vec3R) -> Result<Complex64, SubductionError> {
    let mut angle = 0.0f64;
    for axis in 0..3 {
        angle += wave_vector.get(axis).checked_mul(delta.get(axis))?.to_f64();
    }
    Ok(Complex64::from_polar(1.0, std::f64::consts::TAU * angle))
}

/// The character of `operation` in a typed row.
///
/// The row is indexed by its own Seitz representatives, which may differ from
/// the embedding's by a lattice vector; such entries describe the same element
/// and are related by the Bloch phase, so the lookup matches rotation plus
/// translation *modulo the lattice* and then corrects the phase.  Matching is
/// deliberately not by array order.
///
/// If several representatives are congruent they must all give the same
/// corrected value: that is also what validates the phase sign, since a wrong
/// sign makes them disagree instead of agreeing.
fn character_of(
    row: &CharacterRow,
    ml: &'static str,
    operation: &ExactSeitz,
    lattice: &Lattice,
    wave_vector: &Vec3R,
    index: usize,
) -> Result<Complex64, SubductionError> {
    let mut found: Option<Complex64> = None;
    for (value, candidate) in row.values().iter().zip(row.operations()) {
        if candidate.rotation != flatten_rotation(operation.rotation()) {
            continue;
        }
        let candidate_translation = *exact_operation(candidate)?.translation();
        let delta = operation
            .translation()
            .checked_sub(&candidate_translation)?;
        if !lattice.contains(&delta)? {
            continue;
        }
        let corrected = value * bloch_phase(wave_vector, &delta)?;
        match found {
            Some(existing) if (existing - corrected).norm() > SUBDUCTION_TOLERANCE => {
                return Err(SubductionError::InconsistentCharacterRow { ml, index });
            }
            _ => found = Some(corrected),
        }
    }
    found.ok_or(SubductionError::OperationNotInCharacterRow {
        ml,
        index,
        rotation: flatten_rotation(operation.rotation()),
    })
}

/// Row-major rotation of an exact operation.
fn flatten_rotation(rotation: Mat3I) -> [i32; 9] {
    let mut out = [0i32; 9];
    for (index, value) in out.iter_mut().enumerate() {
        *value = rotation[index / 3][index % 3];
    }
    out
}

/// A seeded character row operation as an exact operation.
fn exact_operation(operation: &SeitzOperation) -> Result<ExactSeitz, SubductionError> {
    let mut rotation = [[0i32; 3]; 3];
    for (row, values) in rotation.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = operation.rotation[3 * row + column];
        }
    }
    let translation = Vec3R::new([
        Rat::from_grid(operation.translation[0], SOURCE_TRANSLATION_GRID)?,
        Rat::from_grid(operation.translation[1], SOURCE_TRANSLATION_GRID)?,
        Rat::from_grid(operation.translation[2], SOURCE_TRANSLATION_GRID)?,
    ]);
    Ok(ExactSeitz::new(rotation, translation))
}

/// A multiplicity that must be a non-negative integer.
fn integral_multiplicity(ml: &'static str, value: Complex64) -> Result<u32, SubductionError> {
    if value.im.abs() > SUBDUCTION_TOLERANCE || value.re < -SUBDUCTION_TOLERANCE {
        return Err(SubductionError::NonIntegralMultiplicity { ml, value });
    }
    let rounded = value.re.round();
    if (value.re - rounded).abs() > SUBDUCTION_TOLERANCE || rounded > f64::from(u32::MAX) {
        return Err(SubductionError::NonIntegralMultiplicity { ml, value });
    }
    // Exact by the checks above: the value is integral and inside u32.
    Ok(rounded as u32)
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
        assert_eq!(
            rat(3, 2).to_i32(),
            Err(SubductionError::NotAnInteger { value: rat(3, 2) })
        );
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
        assert_eq!(
            mapped.translation().as_array(),
            &[rat(1, 16), Rat::ZERO, Rat::ZERO]
        );
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
        assert!(
            !lattice
                .contains(&Vec3R::new([rat(1, 3), Rat::ZERO, Rat::ZERO]))
                .unwrap()
        );
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
        assert!(
            parent
                .same_mod(operation.translation(), shifted.translation())
                .unwrap()
        );
        // ... but different modulo L_H, so they must not be merged.
        assert!(
            !subgroup
                .same_mod(operation.translation(), shifted.translation())
                .unwrap()
        );
        assert_eq!(parent.deduplicate(&[operation, shifted]).unwrap().len(), 1);
        assert_eq!(
            subgroup.deduplicate(&[operation, shifted]).unwrap().len(),
            2
        );
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
        assert!(
            lattice
                .same_mod(identity.translation(), &Vec3R::zero())
                .unwrap()
        );
        assert_eq!(operation.apply(&Vec3R::zero()).unwrap(), vec3([1, 0, 0]));
        let point = Vec3R::new([rat(1, 3), rat(2, 5), rat(-1, 7)]);
        let mapped = operation.apply(&point).unwrap();
        assert_eq!(inverse.apply(&mapped).unwrap(), point);
        // Applying the operation twice equals composing it with itself.
        assert_eq!(
            operation.apply(&mapped).unwrap(),
            operation
                .compose(&operation)
                .unwrap()
                .apply(&point)
                .unwrap()
        );
    }

    #[test]
    fn wave_vector_folding_is_exact_and_does_not_truncate() {
        // A scale far outside what the narrow `i8` k fields can hold: the fold
        // must stay exact rather than saturate or round.
        let transform = SeitzTransform::new(Mat3R::diagonal([1024, 1, 1]), Vec3R::zero());
        let wave_vector = Vec3R::new([rat(1, 3), rat(-1, 6), rat(7, 12)]);
        let folded = fold_wave_vector(&transform, &wave_vector).unwrap();
        assert_eq!(folded.get(0), rat(1024, 3));
        assert_eq!(folded.get(1), rat(-1, 6));
        assert_eq!(folded.get(2), rat(7, 12));
        // And the reciprocal lattice keeps the centring extinctions: for
        // C-centring, (1/2,1/2,0) is a direct lattice vector but (1,0,0) is not
        // a reciprocal one.
        let direct = Lattice::new(Mat3R::new([
            [rat(1, 2), rat(1, 2), Rat::ZERO],
            [rat(-1, 2), rat(1, 2), Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, Rat::ONE],
        ]))
        .unwrap();
        let reciprocal = direct.reciprocal().unwrap();
        assert!(reciprocal.contains(&Vec3R::from_ints([1, 1, 0])).unwrap());
        assert!(!reciprocal.contains(&Vec3R::from_ints([1, 0, 0])).unwrap());
        assert!(reciprocal.contains(&Vec3R::from_ints([2, 0, 0])).unwrap());
    }

    #[test]
    fn folding_sends_a_zone_boundary_point_to_the_subgroup_block() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        // 16 R1 -> #22 F222 is a 2x2x2 supercell, so the parent's X point folds
        // to a non-Gamma block of the subgroup and the parent's zone corners
        // fold onto subgroup points equivalent modulo the F reciprocal lattice.
        let subgroup = isotropy_subgroup_for_direction(16, "R1", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("R1 record");
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).expect("embedding");
        let probe = query::irreps_of(16)
            .iter()
            .find(|record| record.ml == "X1")
            .expect("X1");
        let result = subduce_irrep_with_embedding(&subgroup, &embedding, probe).expect("fold");
        assert_eq!(result.folded_k(), [rat(1, 1), Rat::ZERO, Rat::ZERO]);
        assert_eq!(result.parent_dimension(), 1);
        assert_eq!(result.multiplicity("T1"), 1);
        assert_eq!(result.targets().len(), 1);
    }

    #[test]
    fn unmap_operation_is_the_exact_inverse() {
        let transform = SeitzTransform::new(
            Mat3R::from_ints([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]),
            Vec3R::zero(),
        );
        let operation = ExactSeitz::new([[0, 1, 0], [1, 0, 0], [0, 0, -1]], Vec3R::zero());
        let back = transform.unmap_operation(&operation).unwrap();
        assert_eq!(back.rotation(), [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]);
        // Forward and back are inverse on both parts.
        let forward = transform.map_operation(&back).unwrap();
        assert_eq!(forward, operation);
        // And `unmap_operation` is the *only* inversion: feeding it the
        // already-inverted transform would apply the forward map again.
        let child = ExactSeitz::new([[1, 0, 0], [0, -1, 0], [0, 0, -1]], Vec3R::zero());
        let mapped = transform.map_operation(&child).unwrap();
        assert_eq!(transform.unmap_operation(&mapped).unwrap(), child);
        assert_ne!(
            transform
                .inverse()
                .unwrap()
                .unmap_operation(&mapped)
                .unwrap(),
            child,
            "the already-inverted transform must not be handed to unmap_operation"
        );
    }

    #[test]
    fn gamma_subduction_golden_case() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        let subgroup = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("golden record");
        let result = subduce_irrep(&subgroup, "GM3+").expect("golden decomposition");

        assert_eq!(result.parent_sg(), 221);
        assert_eq!(result.parent_ml(), "GM3+");
        assert_eq!(result.parent_dimension(), 2);
        assert_eq!(result.subgroup_sg(), 83);
        let terms: Vec<(&str, u8, u32)> = result
            .targets()
            .iter()
            .map(|target| (target.ml, target.dimension, target.multiplicity))
            .collect();
        assert_eq!(terms, [("GM1+", 1, 1), ("GM2+", 1, 1)]);
        assert_eq!(result.multiplicity("GM1+"), 1);
        assert_eq!(result.multiplicity("GM2+"), 1);
        // Everything else is zero, including the compound rows of #83.
        for label in ["GM1-", "GM2-", "GM3+GM4+", "GM3-GM4-"] {
            assert_eq!(result.multiplicity(label), 0, "{label} must not appear");
        }
        // 2 = 1 + 1, and the reconstruction is per operation, not just in sum.
        let dimension_sum: u32 = result
            .targets()
            .iter()
            .map(|target| u32::from(target.dimension) * target.multiplicity)
            .sum();
        assert_eq!(dimension_sum, 2);
        let (parent, rebuilt) = result.reconstruction();
        assert_eq!(parent.len(), 8, "one value per coset representative");
        for (index, (expected, found)) in parent.iter().zip(rebuilt).enumerate() {
            assert!(
                (expected - found).norm() <= result.tolerance(),
                "operation {index}: {found} != {expected}"
            );
        }
    }

    /// Multiplicities recomputed straight from the public character rows, by
    /// rotation only.  Valid here because both 221 and #83 are primitive, so a
    /// rotation identifies an operation modulo the lattice; this keeps the
    /// check independent of the embedding's pairing code.
    fn independent_multiplicities(
        subgroup: &crate::irrep::isotropy::IsotropySubgroup,
        probe_ml: &str,
    ) -> Vec<(&'static str, u32)> {
        use crate::irrep::query;
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
        let probe = query::irreps_of(subgroup.parent_sg)
            .iter()
            .find(|record| record.ml == probe_ml)
            .expect("probe");
        let parent_row = probe
            .ordinary_scalar_selected_arm_block_trace()
            .expect("parent row");
        let forward = embedding.transform();
        let child_cell = Lattice::new(exact_primitive_basis(embedding.subgroup_sg()).unwrap())
            .expect("child cell");
        let mut pairs = Vec::new();
        for operation in embedding.representatives() {
            let parent_index = parent_row
                .operations()
                .iter()
                .position(|candidate| candidate.rotation == flatten_rotation(operation.rotation()))
                .expect("parent character");
            let child = forward
                .unmap_operation(operation)
                .expect("pull back")
                .reduce(&child_cell)
                .expect("reduce");
            pairs.push((parent_row.values()[parent_index], child));
        }
        // Every target row expanded into its complex constituents, exactly as
        // the metadata describes them.
        let mut rows: Vec<(&'static str, CharacterRow, bool)> = Vec::new();
        for record in query::irreps_of(embedding.subgroup_sg()) {
            if record.spinor || !is_gamma(record) {
                continue;
            }
            if let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() {
                rows.push((record.ml, row, false));
                continue;
            }
            if let Ok(view) = record.compound_selected_arm_view() {
                match view {
                    crate::irrep::types::CompoundSelectedArmCharacter::DistinctComponentSum {
                        first,
                        second,
                        ..
                    } => {
                        rows.push((first.label, first.row, false));
                        rows.push((second.label, second.row, false));
                    }
                    crate::irrep::types::CompoundSelectedArmCharacter::ConjugateRealification {
                        seed,
                        ..
                    } => {
                        rows.push((seed.label, seed.row.clone(), false));
                        rows.push((seed.label, seed.row, true));
                    }
                }
            }
        }
        let mut out = Vec::new();
        for (label, row, conjugate) in rows {
            let mut sum = Complex64::new(0.0, 0.0);
            for (parent_value, child) in &pairs {
                let index = row
                    .operations()
                    .iter()
                    .position(|candidate| candidate.rotation == flatten_rotation(child.rotation()))
                    .expect("child character");
                let value = row.values()[index];
                let value = if conjugate { value.conj() } else { value };
                sum += parent_value * value.conj();
            }
            let value = sum / pairs.len() as f64;
            let rounded = value.re.round();
            if value.im.abs() < 1e-7 && (value.re - rounded).abs() < 1e-7 && rounded > 0.5 {
                out.push((label, rounded as u32));
            }
        }
        out
    }

    #[test]
    fn subduction_is_character_driven_not_hard_coded() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        let subgroup = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("golden record");
        let mut probes = 0;
        for record in crate::irrep::query::irreps_of(221) {
            if record.spinor || !is_gamma(record) {
                continue;
            }
            if record.ordinary_scalar_selected_arm_block_trace().is_err() {
                continue;
            }
            let result = subduce_irrep(&subgroup, record.ml)
                .unwrap_or_else(|error| panic!("{}: {error}", record.ml));
            let expected = independent_multiplicities(&subgroup, record.ml);
            let found: Vec<(&str, u32)> = result
                .targets()
                .iter()
                .map(|target| (target.ml, target.multiplicity))
                .collect();
            assert_eq!(found, expected, "{} multiplicities", record.ml);
            let dimension_sum: u32 = result
                .targets()
                .iter()
                .map(|target| u32::from(target.dimension) * target.multiplicity)
                .sum();
            assert_eq!(
                dimension_sum,
                u32::from(record.dim),
                "{} dimension sum",
                record.ml
            );
            // The golden case is one probe among many, not a special path.
            if record.ml == "GM3+" {
                assert_eq!(found, [("GM1+", 1), ("GM2+", 1)]);
            }
            probes += 1;
        }
        assert!(
            probes >= 9,
            "checked {probes} scalar Gamma probes of SG 221"
        );
    }

    #[test]
    fn compound_rows_are_expanded_into_complex_constituents() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        use crate::irrep::types::CompoundSelectedArmCharacter;
        let subgroup = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("golden record");
        // 221 GM4+ is a 3-dimensional ordinary irrep; #83 realises its A_u part
        // as an ordinary row and its E_u part only as a compound row, so this
        // case fails unless the compound row is expanded by its metadata.
        let result = subduce_irrep(&subgroup, "GM4+").expect("decomposition of GM4+");
        let terms: Vec<(&str, &str, u32, SubductionComponent)> = result
            .targets()
            .iter()
            .map(|target| {
                (
                    target.ml,
                    target.row_ml,
                    target.multiplicity,
                    target.component,
                )
            })
            .collect();
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].0, "GM1+");
        assert_eq!(terms[0].3, SubductionComponent::Ordinary);
        assert_eq!(terms[0].2, 1);
        for (index, expected_irnumber) in [(0usize, 4077u32), (1, 4078)] {
            let term = terms[index + 1];
            assert_eq!(
                term.1, "GM3+GM4+",
                "the compound row label is kept as provenance"
            );
            assert_eq!(term.2, 1);
            assert_eq!(
                term.3,
                SubductionComponent::Constituent {
                    index: index as u8,
                    irnumber: expected_irnumber,
                }
            );
        }
        // The two constituents are distinct complex irreps with distinct source
        // numbers, which is what the metadata says the row is.
        let view = crate::irrep::query::irreps_of(83)
            .iter()
            .find(|record| record.ml == "GM3+GM4+")
            .expect("compound row")
            .compound_selected_arm_view()
            .expect("compound view");
        let CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } = view else {
            panic!("GM3+GM4+ must be a distinct component sum");
        };
        assert_eq!(first.irnumber, 4077);
        assert_eq!(second.irnumber, 4078);
        assert_ne!(first.label, second.label);
        // 3 = 1 + 1 + 1 and every operation is rebuilt.
        let dimension_sum: u32 = result
            .targets()
            .iter()
            .map(|target| u32::from(target.dimension) * target.multiplicity)
            .sum();
        assert_eq!(dimension_sum, 3);
        let (parent, rebuilt) = result.reconstruction();
        for (expected, found) in parent.iter().zip(rebuilt) {
            assert!((expected - found).norm() <= result.tolerance());
        }
    }

    #[test]
    fn a_realification_reports_seed_and_conjugate_separately() {
        // Adapter-level witness for the expansion of a genuine
        // `ConjugateRealification` row on an explicit four-element translation
        // quotient: SG 23 `W1W1` at k = (1/2,1/2,1/2), probed on the pure
        // translations 0, L, 2L, 3L with L = (1/2,1/2,1/2).  The row is a
        // realification, so the expansion must yield the CIR seed *and* its
        // conjugate as separate complex constituents.
        //
        // Boundary, stated honestly: this exercises component expansion (and
        // the unequal multiplicities it makes representable), NOT a full
        // supported subgroup decomposition.  There is no full-star grouping of
        // k and -k here, and the nine frozen child space groups contain zero
        // realification rows (and zero Gamma realification rows), so the
        // end-to-end path cannot reach this case yet.
        use crate::irrep::query;

        let record = query::irreps_of(23)
            .iter()
            .find(|record| record.ml == "W1W1" && !record.spinor)
            .expect("SG 23 W1W1");
        assert_eq!((record.kx, record.ky, record.kz, record.kd), (1, 1, 1, 2));
        let CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } =
            record.compound_selected_arm_view().expect("compound view")
        else {
            panic!("SG 23 W1W1 must be a conjugate realification");
        };
        assert_eq!(seed.dimension, 1);

        let half = Rat::new(1, 2).expect("1/2");
        let mut translations = Vec::new();
        for step in 0..4i128 {
            let scale = Rat::from_integer(step);
            let component = half.checked_mul(scale).expect("n/2");
            translations.push(Vec3R::new([component; 3]));
        }
        let pulled_back: Vec<ExactSeitz> = translations
            .iter()
            .map(|translation| ExactSeitz::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], *translation))
            .collect();
        let block = [record];
        let wave_vector = Vec3R::new([half; 3]);
        let targets = complex_targets(
            23,
            &block,
            &pulled_back,
            &Lattice::new(exact_primitive_basis(23).unwrap()).unwrap(),
            &wave_vector,
        )
        .expect("expansion of the realification row");
        assert_eq!(targets.len(), 2, "seed and conjugate are separate targets");

        let mut seed_values = None;
        let mut conjugate_values = None;
        for target in &targets {
            assert_eq!(target.row_ml, "W1W1");
            assert_eq!(target.irnumber, seed.irnumber);
            assert_eq!(usize::from(target.dimension), seed.dimension);
            match target.component {
                SubductionComponent::RealificationSeed { irnumber } => {
                    assert_eq!(irnumber, seed.irnumber);
                    seed_values = Some(target.values.clone());
                }
                SubductionComponent::RealificationConjugate { irnumber } => {
                    assert_eq!(irnumber, seed.irnumber);
                    conjugate_values = Some(target.values.clone());
                }
                other => panic!("unexpected component {other:?}"),
            }
        }
        let seed_values = seed_values.expect("seed component");
        let conjugate_values = conjugate_values.expect("conjugate component");

        // Expected characters: chi(t = nL) = exp(2 pi i k . nL) with
        // k . L = 3/4 gives 1, -i, -1, +i; the conjugate flips the sign of i.
        let expected_seed = [
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, -1.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(0.0, 1.0),
        ];
        let expected_conjugate: Vec<Complex64> =
            expected_seed.iter().map(|value| value.conj()).collect();
        for (found, expected) in seed_values.iter().zip(&expected_seed) {
            assert!(
                (found - expected).norm() < SUBDUCTION_TOLERANCE,
                "{found} != {expected}"
            );
        }
        for (found, expected) in conjugate_values.iter().zip(&expected_conjugate) {
            assert!(
                (found - expected).norm() < SUBDUCTION_TOLERANCE,
                "{found} != {expected}"
            );
        }

        // Inner products over this cyclic translation quotient Z4.  A seed
        // probe has (1,0), its conjugate (0,1): the two components are
        // orthogonal and each is normalised.
        let overlap = |probe: &[Complex64], values: &[Complex64]| -> Complex64 {
            probe
                .iter()
                .zip(values)
                .map(|(left, right)| left * right.conj())
                .sum::<Complex64>()
                * (1.0 / probe.len() as f64)
        };
        let seed_probe = overlap(&seed_values, &seed_values);
        let seed_against_conjugate = overlap(&seed_values, &conjugate_values);
        assert!(
            (seed_probe - 1.0).norm() < SUBDUCTION_TOLERANCE,
            "{seed_probe}"
        );
        assert!(
            seed_against_conjugate.norm() < SUBDUCTION_TOLERANCE,
            "{seed_against_conjugate}"
        );
        let conjugate_probe = overlap(&conjugate_values, &conjugate_values);
        let conjugate_against_seed = overlap(&conjugate_values, &seed_values);
        assert!(
            (conjugate_probe - 1.0).norm() < SUBDUCTION_TOLERANCE,
            "{conjugate_probe}"
        );
        assert!(
            conjugate_against_seed.norm() < SUBDUCTION_TOLERANCE,
            "{conjugate_against_seed}"
        );
    }

    #[test]
    fn subduction_depends_on_the_embedding_context() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        // The same parent irrep on a different condensing direction and
        // subgroup: the decomposition must follow the embedding, not the label.
        let first = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("P1 record");
        let second = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P2"))
            .expect("P2 record");
        let one = subduce_irrep(&first, "GM3+").expect("decomposition on #83");
        let two = subduce_irrep(&second, "GM3+").expect("decomposition on #12");
        assert_eq!(one.subgroup_sg(), 83);
        assert_eq!(two.subgroup_sg(), 12);
        assert_eq!(one.ordinal(), 12400);
        assert_eq!(two.ordinal(), 12401);
        // The labels of the two subgroup point groups coincide (`GM1+`, `GM2+`
        // of #83 versus #12 are different irreps), so the decomposition has to
        // be compared through the embedding: the coset sets differ in size and
        // each result must reproduce its own parent characters.
        assert_eq!(one.reconstruction().0.len(), 8);
        assert_eq!(two.reconstruction().0.len(), 4);
        assert_ne!(one.setting(), two.setting());
        assert_ne!(one.ordinal(), two.ordinal());
        for result in [&one, &two] {
            let (parent, rebuilt) = result.reconstruction();
            for (expected, found) in parent.iter().zip(rebuilt) {
                assert!((expected - found).norm() <= result.tolerance());
            }
        }
        for result in [&one, &two] {
            let sum: u32 = result
                .targets()
                .iter()
                .map(|target| u32::from(target.dimension) * target.multiplicity)
                .sum();
            assert_eq!(sum, 2);
        }
    }

    #[test]
    fn subduction_rejects_unsupported_probes_and_contexts() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        let subgroup = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("golden record");
        assert!(matches!(
            subduce_irrep(&subgroup, "NOPE"),
            Err(SubductionError::ProbeIrrepNotFound { .. })
        ));
        // A label that exists only away from Gamma is reported as such.
        let away = crate::irrep::query::irreps_of(221)
            .iter()
            .find(|record| !is_gamma(record) && !record.spinor)
            .expect("a non-Gamma irrep");
        assert!(matches!(
            subduce_irrep(&subgroup, away.ml),
            Err(SubductionError::ProbeNotAtGamma { .. })
        ));
        // A spinor label at Gamma is refused instead of silently skipped.
        if let Some(spinor) = crate::irrep::query::irreps_of(221)
            .iter()
            .find(|record| record.spinor && is_gamma(record))
        {
            assert!(matches!(
                subduce_irrep(&subgroup, spinor.ml),
                Err(SubductionError::UnsupportedCharacterSpace { .. })
            ));
        }
        // A condensing irrep away from Gamma is out of scope here (task 7).
        let away_subgroup =
            crate::irrep::isotropy::isotropy_subgroups(221, away.ml, LabelConvention::Cdml)
            .expect("subgroups")
            .into_iter()
            .next();
        if let Some(away_subgroup) = away_subgroup {
            assert!(matches!(
                subduce_irrep(&away_subgroup, "GM3+"),
                Err(SubductionError::ProbeNotAtGamma { .. })
            ));
        }
    }

    #[test]
    fn subgroup_embeddings_build_for_the_fixture_cases() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        let cases = [
            (221u8, "GM4+", "P1"),
            (221, "GM4+", "P2"),
            (221, "GM4+", "P3"),
            (221, "GM3+", "P1"),
            (221, "GM3+", "C1"),
            (225, "GM4-", "C1"),
            (225, "GM4-", "C2"),
            (16, "R1", "P1"),
            (167, "GM3+", "P1"),
            (139, "M1-", "P1"),
        ];
        for (sg, ml, label) in cases {
            let subgroup = isotropy_subgroup_for_direction(sg, ml, LabelConvention::Cdml, IsotropyDirection::Label(label))
                .unwrap_or_else(|error| panic!("SG {sg} {ml} {label}: {error}"));
            let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
                .unwrap_or_else(|error| panic!("SG {sg} {ml} {label}: {error}"));
            println!(
                "SG {sg} {ml} {label} -> #{} candidates {} setting {:?} ops {} reps {}",
                embedding.subgroup_sg(),
                embedding.candidate_count(),
                embedding.setting(),
                embedding.operations().len(),
                embedding.representatives().len(),
            );
            assert_eq!(embedding.parent_sg(), sg);
            assert_eq!(embedding.subgroup_sg(), subgroup.record.sg as u8);
            assert!(!embedding.representatives().is_empty());
        }
    }

    #[test]
    fn frozen_settings_address_their_own_records() {
        use super::settings_data::{FROZEN_EMBEDDING_SETTINGS, FrozenEmbeddingSetting};
        let entries: &[FrozenEmbeddingSetting] = FROZEN_EMBEDDING_SETTINGS;
        assert!(!entries.is_empty());
        let mut ordinals = Vec::new();
        let mut pairs = Vec::new();
        let mut shifted = Vec::new();
        for entry in entries {
            let (ordinal, parent, child, setting, setting_denominator, child_shift) = *entry;
            assert!(setting_denominator > 0, "ordinal {ordinal}: denominator");
            assert!(!ordinals.contains(&ordinal), "ordinal {ordinal} repeats");
            ordinals.push(ordinal);
            if !pairs.contains(&(parent, child)) {
                pairs.push((parent, child));
            }
            // The ordinal addresses the generated isotropy record, and an irrep
            // of the claimed parent really owns that record.
            let stored = &ISOTROPY_SUBGROUPS[ordinal];
            assert_eq!(stored.sg, usize::from(child), "ordinal {ordinal}");
            assert!(
                query::irreps_of(parent).iter().any(|irrep| {
                    !irrep.spinor
                        && irrep
                            .subgroups()
                            .iter()
                            .any(|record| std::ptr::eq(record, stored))
                }),
                "ordinal {ordinal}: no irrep of SG {parent} owns the record"
            );
            // `U = setting / setting_denominator` is an exact unimodular
            // rational matrix; the generated table is not restricted to signed
            // permutations, and a few monoclinic records need the fraction.
            let matrix = super::rational_setting(setting, setting_denominator)
                .expect("frozen setting is a valid rational matrix");
            let determinant = matrix.determinant().expect("integer matrix");
            assert!(
                determinant == Rat::ONE || determinant == Rat::from_integer(-1),
                "ordinal {ordinal}: det U = {determinant:?}"
            );
            let [x, y, z, denominator] = child_shift;
            assert!(denominator > 0, "ordinal {ordinal}");
            let divisor = gcd_positive(
                gcd_positive(
                    i128::from(x).unsigned_abs(),
                    i128::from(y).unsigned_abs(),
                ),
                gcd_positive(
                    i128::from(z).unsigned_abs(),
                    i128::from(denominator).unsigned_abs(),
                ),
            );
            assert_eq!(divisor, 1, "ordinal {ordinal}: unreduced child shift");
            if child_shift != [0, 0, 0, 1] {
                shifted.push((parent, child, child_shift));
            }
        }
        // #126 is the documented ITA origin-choice correction; it must stay in
        // the table whatever else the generator freezes.
        assert!(
            shifted.contains(&(139, 126, [1, 1, 1, 4])),
            "the #126 origin-choice shift disappeared: {shifted:?}"
        );
        assert!(!pairs.is_empty());
    }

    #[test]
    fn frozen_setting_lookup_fails_closed_on_a_mismatch() {
        let covered = FROZEN_EMBEDDING_SETTINGS[0];
        let (ordinal, parent, child) = (covered.0, covered.1, covered.2);
        assert_eq!(
            frozen_setting_for(ordinal, parent, child).expect("lookup"),
            Some(&covered)
        );
        // A covered ordinal under another parent or subgroup is an error, not a
        // reason to search: the recorded setting belongs to this embedding only.
        let other_parent = if parent == 1 { 2 } else { parent - 1 };
        let other_child = if child == 1 { 2 } else { child - 1 };
        assert_eq!(
            frozen_setting_for(ordinal, other_parent, child),
            Err(SubductionError::StaleIsotropyRecord { ordinal })
        );
        assert_eq!(
            frozen_setting_for(ordinal, parent, other_child),
            Err(SubductionError::StaleIsotropyRecord { ordinal })
        );
        // The committed table covers every pinned record, so the candidate
        // search is a safety net rather than the normal path.  An ordinal no
        // pinned record owns is not an error either.
        let out_of_range = ISOTROPY_SUBGROUPS.len();
        assert_eq!(
            frozen_setting_for(out_of_range, 221, 83).expect("lookup"),
            None
        );
        assert_eq!(
            FROZEN_EMBEDDING_SETTINGS.len(),
            ISOTROPY_SUBGROUPS.len(),
            "the frozen table must address every pinned isotropy record"
        );
    }

    #[test]
    fn forged_isotropy_records_are_rejected() {
        use crate::irrep::LabelConvention;
        use crate::irrep::isotropy::{IsotropyDirection, isotropy_subgroup_for_direction};
        let subgroup = isotropy_subgroup_for_direction(221, "GM4+", LabelConvention::Cdml, IsotropyDirection::Label("P1"))
            .expect("golden record");
        assert!(SubgroupEmbedding::from_isotropy_subgroup(&subgroup).is_ok());

        let mut moved = subgroup;
        moved.ordinal = 12401;
        assert!(matches!(
            SubgroupEmbedding::from_isotropy_subgroup(&moved),
            Err(SubductionError::StaleIsotropyRecord { ordinal: 12401 })
        ));

        let mut tampered = subgroup;
        tampered.record.basis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        assert!(matches!(
            SubgroupEmbedding::from_isotropy_subgroup(&tampered),
            Err(SubductionError::StaleIsotropyRecord { ordinal: 12400 })
        ));

        let mut wrong_irrep = subgroup;
        wrong_irrep.irrep_ml = "GM3+";
        assert!(matches!(
            SubgroupEmbedding::from_isotropy_subgroup(&wrong_irrep),
            Err(SubductionError::StaleIsotropyRecord { ordinal: 12400 })
        ));

        let mut out_of_range = subgroup;
        out_of_range.ordinal = ISOTROPY_SUBGROUPS.len();
        assert!(matches!(
            SubgroupEmbedding::from_isotropy_subgroup(&out_of_range),
            Err(SubductionError::IsotropyOrdinalOutOfRange { .. })
        ));
    }

    #[test]
    fn signed_permutations_are_the_48_unimodular_axis_permutations() {
        let permutations = signed_permutations();
        assert_eq!(permutations.len(), 48);
        let mut unique: Vec<Mat3I> = Vec::new();
        for matrix in &permutations {
            assert!(
                Mat3R::from_ints(*matrix)
                    .determinant()
                    .unwrap()
                    .numerator()
                    .abs()
                    == 1
            );
            assert!(!unique.contains(matrix), "duplicate signed permutation");
            unique.push(*matrix);
        }
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
