//! Common types for irreducible representation data.
//!
//! Based on Stokes & Hatch (1988), *Isotropy Subgroups of the 230
//! Crystallographic Space Groups*.  Three irrep labeling conventions are
//! supported.
//!
//! # Labeling conventions
//!
//! | Convention | Reference | Example |
//! |-----------|-----------|---------|
//! | Miller & Love | Miller & Love (1967) | `GM1+`, `X2-` |
//! | Kovalev | Kovalev (1986) | `τ1`, `k6τ2` |
//! | Bradley & Cracknell | B&C (1972) | `Γ1+`, `X1` |

use num_complex::Complex64;

// ── Compact record types (flat-array storage) ───────────────────────────────

/// How a scalar physical compound record is assembled from CIR rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompoundCharacterSemantics {
    /// A CIR row and its complex conjugate are realified: χ = 2 Re(χ_CIR).
    ConjugateRealification,
    /// Distinct CIR constituents are directly summed: χ = Σ χ_CIR.
    DistinctComponentSum,
}

/// Version of the generation-time Miller--Love compound naming grammar.
pub const COMPOUND_NAMING_GRAMMAR_VERSION: u8 = 1;

/// Provenance of the generation-time PIR/CIR compound association.
pub const COMPOUND_NAMING_PROVENANCE: &str = "ISO-IR Miller-Love concatenation, resolver v1";

/// Generation-time provenance for a compound physical irrep.
///
/// ISO-IR does not provide an independent PIR-to-CIR link. The two constituent
/// identities are resolved from the documented Miller--Love concatenation
/// grammar (version 1) against unique same-SG CIR source records, then frozen
/// here as generated read-only metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompoundMetadata {
    /// Space group number shared by the data_irreps record and both CIR rows.
    pub sg: u8,
    /// Original compound Miller--Love label from data_irreps.txt and the
    /// runtime [`IrrepRecord`]. This is not a PIR_data.txt source label.
    pub record_label: &'static str,
    /// Stable source `irnumber`s from CIR_data.txt, in constituent order.
    pub cir_irnumbers: [u32; 2],
    /// Constituent Miller--Love labels resolved from CIR_data.txt.
    pub cir_labels: [&'static str; 2],
    /// Dimensions of the two selected-arm CIR constituents.
    pub cir_dimensions: [u8; 2],
    /// Version of the naming grammar used for the resolution.
    pub naming_grammar_version: u8,
    /// Global provenance string; generated rows reference
    /// [`COMPOUND_NAMING_PROVENANCE`] rather than duplicating its text.
    pub provenance: &'static str,
    /// Algebra used to assemble the physical character row.
    pub semantics: CompoundCharacterSemantics,
}

/// Error returned when stored irrep matrices cannot be aligned with a supplied
/// symmetry-operation order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum MatrixReorderError {
    /// No unambiguous rotation mapping exists between the two operation lists.
    #[error("could not map supplied symmetry operations to stored irrep rotations")]
    OperationMappingFailed,
}

/// Representation space occupied by a typed character row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RepresentationSpaceKind {
    /// The physical full-star/PIR character row in PIR operation order.
    FullStarPir,
    /// A selected star arm's block trace, still indexed by the complete PIR
    /// operation universe; callers must apply the arm little-group Seitz set.
    SelectedArmBlockTrace,
    /// One raw complex CIR constituent of a compound record.
    ConstituentCir,
}

/// A symmetry operation paired with one typed character value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeitzOperation {
    /// Rotation matrix in row-major order.
    pub rotation: [i32; 9],
    /// Fractional translation in the same setting as the character row.
    pub translation: [f64; 3],
}

/// Errors reported when a typed character row cannot be constructed safely.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum CharacterViewError {
    /// This representation space does not apply to the record.
    #[error("typed character space is not applicable to this irrep")]
    NotApplicable,
    /// Required generated character or operation data is absent.
    #[error("typed character data is missing")]
    MissingData,
    /// A dimension supplied by generated metadata is zero or invalid.
    #[error("typed character dimension is invalid")]
    InvalidDimension,
    /// A character, translation, or spin lift contains a non-finite value.
    #[error("typed character data contains a non-finite value")]
    NonFiniteData,
    /// A generated flat-array offset is not aligned with its element width.
    #[error("typed character storage offset is misaligned")]
    MisalignedStorage,
    /// Values and complete Seitz operations have different lengths.
    #[error("typed character values and operations have different lengths")]
    LengthMismatch,
    /// A generated operation index points outside the SG operation table.
    #[error("typed character operation index is out of bounds")]
    InvalidOperationIndex,
    /// A compound component selector is outside the two metadata components.
    #[error("typed compound component selector is invalid")]
    InvalidComponent,
    /// A spinor little-group operation index occurs more than once.
    #[error("typed spinor operation index occurs more than once")]
    DuplicateOperationIndex,
    /// The complete Seitz operation list does not contain exactly one identity.
    #[error("typed character operation list does not contain a unique identity")]
    NonUniqueIdentity,
    /// The identity character is not a positive integral dimension.
    #[error("typed character identity value is not a dimension")]
    InvalidIdentityDimension,
    /// The identity character disagrees with the expected dimension.
    #[error("typed character identity disagrees with expected dimension")]
    DimensionMismatch,
    /// A compound's constituent rows cannot share one operation order.
    #[error("compound constituent operation orders are inconsistent")]
    OperationOrderMismatch,
}

const SEITZ_IDENTITY_TOLERANCE: f64 = 1e-8;

/// Owned, operation-aware character values in one explicitly named space.
#[derive(Debug, Clone)]
pub struct CharacterRow {
    space: RepresentationSpaceKind,
    dimension: usize,
    values: Vec<Complex64>,
    operations: Vec<SeitzOperation>,
}

impl CharacterRow {
    fn from_parts(
        space: RepresentationSpaceKind,
        values: Vec<Complex64>,
        operations: Vec<SeitzOperation>,
        expected_dimension: Option<usize>,
    ) -> Result<Self, CharacterViewError> {
        if expected_dimension == Some(0) {
            return Err(CharacterViewError::InvalidDimension);
        }
        if values.is_empty() || operations.is_empty() {
            return Err(CharacterViewError::MissingData);
        }
        if values.len() != operations.len() {
            return Err(CharacterViewError::LengthMismatch);
        }
        if values
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
            || operations
                .iter()
                .any(|operation| operation.translation.iter().any(|value| !value.is_finite()))
        {
            return Err(CharacterViewError::NonFiniteData);
        }
        let identity_indices = operations
            .iter()
            .enumerate()
            .filter(|(_, operation)| is_identity_seitz(operation))
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        if identity_indices.len() != 1 {
            return Err(CharacterViewError::NonUniqueIdentity);
        }
        let identity = values[identity_indices[0]];
        if identity.im.abs() > SEITZ_IDENTITY_TOLERANCE
            || identity.re <= 0.0
            || (identity.re - identity.re.round()).abs() > SEITZ_IDENTITY_TOLERANCE
        {
            return Err(CharacterViewError::InvalidIdentityDimension);
        }
        let dimension = identity.re.round() as usize;
        if expected_dimension.is_some_and(|expected| expected != dimension) {
            return Err(CharacterViewError::DimensionMismatch);
        }
        Ok(Self {
            space,
            dimension,
            values,
            operations,
        })
    }

    /// Representation space of this row.
    pub fn representation_space(&self) -> RepresentationSpaceKind {
        self.space
    }

    /// Dimension of the representation in this row's space.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Number of values and paired complete Seitz operations.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether the row has no character values.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Borrow all complex character values in operation order.
    pub fn values(&self) -> &[Complex64] {
        &self.values
    }

    /// Get one complex character value by operation position.
    pub fn get(&self, index: usize) -> Option<Complex64> {
        self.values.get(index).copied()
    }

    /// Get one character together with its operation as an inseparable entry.
    pub fn entry(&self, index: usize) -> Option<(Complex64, SeitzOperation)> {
        Some((
            self.values.get(index).copied()?,
            self.operations.get(index).copied()?,
        ))
    }

    /// Borrow the complete Seitz operations paired with [`Self::values`].
    pub fn operations(&self) -> &[SeitzOperation] {
        &self.operations
    }

    /// Get one complete Seitz operation by character position.
    pub fn operation(&self, index: usize) -> Option<SeitzOperation> {
        self.operations.get(index).copied()
    }
}

fn is_identity_seitz(operation: &SeitzOperation) -> bool {
    operation.rotation == [1, 0, 0, 0, 1, 0, 0, 0, 1]
        && operation
            .translation
            .iter()
            .all(|value| (value - value.round()).abs() <= SEITZ_IDENTITY_TOLERANCE)
}

fn validate_spinor_operation_indices(
    indices: &[u16],
    operation_count: usize,
) -> Result<(), CharacterViewError> {
    let mut seen = vec![false; operation_count];
    for &index in indices {
        let index = index as usize;
        if index >= operation_count {
            return Err(CharacterViewError::InvalidOperationIndex);
        }
        if seen[index] {
            return Err(CharacterViewError::DuplicateOperationIndex);
        }
        seen[index] = true;
    }
    Ok(())
}

/// One compound CIR constituent with its frozen C1 provenance.
#[derive(Debug, Clone)]
pub struct CompoundConstituentCharacter {
    /// Zero-based constituent position in the compound record.
    pub component: usize,
    /// Stable CIR source identifier.
    pub irnumber: u32,
    /// CIR Miller--Love label.
    pub label: &'static str,
    /// Selected-arm CIR dimension.
    pub dimension: usize,
    /// Raw constituent row in the shared aligned Seitz order.
    pub row: CharacterRow,
}

/// A compound selected-arm view whose public shape follows its C1 semantics.
#[derive(Debug, Clone)]
pub enum CompoundSelectedArmCharacter {
    /// One authoritative CIR seed, its generated physical realification, and
    /// no second stored row exposed as an independent constituent.
    ConjugateRealification {
        seed: CompoundConstituentCharacter,
        block_trace: CharacterRow,
    },
    /// Two distinct CIR constituents and their direct physical sum.
    DistinctComponentSum {
        first: CompoundConstituentCharacter,
        second: CompoundConstituentCharacter,
        block_trace: CharacterRow,
    },
}

impl CompoundSelectedArmCharacter {
    /// Complete selected-arm block-trace row. Callers must match it to the
    /// desired little-group Seitz operation set.
    pub fn block_trace(&self) -> &CharacterRow {
        match self {
            Self::ConjugateRealification { block_trace, .. }
            | Self::DistinctComponentSum { block_trace, .. } => block_trace,
        }
    }
}

/// A spin-space operation: a Seitz operation paired with its SU(2) lift.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpinSeitzOperation {
    /// The spatial operation in the same setting as the character.
    pub seitz: SeitzOperation,
    /// Pauli coefficients `[u₀, u₁, u₂, u₃]` of the spin lift.
    pub su2: [f64; 4],
}

/// Operation-aware complex characters in the spinor selected-arm space.
#[derive(Debug, Clone)]
pub struct SpinCharacterRow {
    dimension: usize,
    values: Vec<Complex64>,
    operations: Vec<SpinSeitzOperation>,
}

impl SpinCharacterRow {
    fn from_parts(
        values: Vec<Complex64>,
        operations: Vec<SpinSeitzOperation>,
        expected_dimension: usize,
    ) -> Result<Self, CharacterViewError> {
        if expected_dimension == 0 {
            return Err(CharacterViewError::InvalidDimension);
        }
        if values.is_empty() || operations.is_empty() {
            return Err(CharacterViewError::MissingData);
        }
        if values.len() != operations.len() {
            return Err(CharacterViewError::LengthMismatch);
        }
        if values
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
            || operations.iter().any(|operation| {
                operation
                    .seitz
                    .translation
                    .iter()
                    .any(|value| !value.is_finite())
                    || operation.su2.iter().any(|value| !value.is_finite())
            })
        {
            return Err(CharacterViewError::NonFiniteData);
        }
        let identity_indices =
            operations
                .iter()
                .enumerate()
                .filter(|(_, operation)| {
                    is_identity_seitz(&operation.seitz)
                        && operation.su2.iter().zip([1.0, 0.0, 0.0, 0.0]).all(
                            |(value, expected)| {
                                (value - expected).abs() <= SEITZ_IDENTITY_TOLERANCE
                            },
                        )
                })
                .map(|(index, _)| index)
                .collect::<Vec<_>>();
        if identity_indices.len() != 1 {
            return Err(CharacterViewError::NonUniqueIdentity);
        }
        let identity = values[identity_indices[0]];
        if identity.im.abs() > SEITZ_IDENTITY_TOLERANCE
            || identity.re <= 0.0
            || (identity.re - identity.re.round()).abs() > SEITZ_IDENTITY_TOLERANCE
        {
            return Err(CharacterViewError::InvalidIdentityDimension);
        }
        let dimension = identity.re.round() as usize;
        if dimension != expected_dimension {
            return Err(CharacterViewError::DimensionMismatch);
        }
        Ok(Self {
            dimension,
            values,
            operations,
        })
    }

    /// Dimension of the spinor representation.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Number of paired spin characters and operations.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether this row is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Borrow complex spinor characters in little-group operation order.
    pub fn values(&self) -> &[Complex64] {
        &self.values
    }

    /// Borrow spin operations paired with [`Self::values`].
    pub fn operations(&self) -> &[SpinSeitzOperation] {
        &self.operations
    }

    /// Get one complex character and its full spin operation together.
    pub fn entry(&self, index: usize) -> Option<(Complex64, SpinSeitzOperation)> {
        Some((
            self.values.get(index).copied()?,
            self.operations.get(index).copied()?,
        ))
    }
}

/// Rational high-symmetry wave vector in reciprocal-lattice coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KVector {
    /// Cartesian reciprocal-coordinate numerators.
    pub numerators: [i8; 3],
    /// Common denominator for all three components.
    pub denominator: i8,
}

impl KVector {
    /// Construct a rational reciprocal vector.
    pub const fn new(numerators: [i8; 3], denominator: i8) -> Self {
        Self {
            numerators,
            denominator,
        }
    }
}

/// Compact irrep record for the generated flat array.
///
/// Field names are abbreviated to keep the generated code size manageable.
///
/// # Navigation
///
/// Use [`IrrepRecord::subgroups`] to get isotropy subgroups directly,
/// without needing to know global indices.
#[derive(Debug, Clone, Copy)]
pub struct IrrepRecord {
    /// Space group number (1–230)
    pub sg: u8,
    /// CDML / Miller-Love label: `"GM4+"`, `"X1-"`
    pub ml: &'static str,
    /// Bradley-Cracknell label (LaTeX): `"\\Gamma_4^+"`
    pub bc: &'static str,
    /// Kovalev label (LaTeX): `"k_{12}\\tau_{9}"`
    pub kov: &'static str,
    /// Dimension: 1, 2, 3, 4, 6, 8, 12, 16, 24
    pub dim: u8,
    /// Image symbol: `"A1a"`, `"C24c"`, `"B6a"`
    pub image: &'static str,
    /// Lifshitz condition satisfied (scalar irreps only)
    pub lifshitz: bool,
    /// Whether this is a double-valued (spinor) irrep
    pub spinor: bool,

    /// k-vector numerator x (fractional reciprocal coordinate)
    pub kx: i8,
    /// k-vector numerator y (fractional reciprocal coordinate)
    pub ky: i8,
    /// k-vector numerator z (fractional reciprocal coordinate)
    pub kz: i8,
    /// k-vector common denominator (actual coordinate = numerator / denominator)
    pub kd: i8,

    // ── internal: character table + matrix pointers ──
    /// Start index into [`CHARACTERS`]
    pub(crate) _char_start: u32,
    /// Number of operators (= number of character values)
    pub(crate) _char_count: u16,
    /// Authoritative selected-arm dimension from the exact CIR header for
    /// ordinary scalar irreps; zero for compounds and spinors.
    pub(crate) _scalar_selected_dim: u8,
    /// For spinor irreps: number of little-group operations, exactly equal to
    /// `_char_count` in the generated selected-arm data.
    pub(crate) _spin_lg_count: u8,
    /// Start index into [`MATRICES`] (u32: ~1M entries total)
    pub(crate) _mat_start: u32,
    /// Number of matrix elements = opcount × dim². Centered conventional
    /// Hall expansion can exceed `u16` for the largest induced irreps.
    pub(crate) _mat_count: u32,
    /// Start index into [`ISOTROPY_SUBGROUPS`]
    pub(crate) _iso_start: u16,
    /// Number of isotropy subgroups for this irrep
    pub(crate) _iso_count: u16,
    /// Start index into [`MAGNETIC_ISOTROPY_SUBGROUPS`]
    pub(crate) _mag_iso_start: u16,
    /// Number of magnetic isotropy subgroups for this irrep
    pub(crate) _mag_iso_count: u16,
    /// Start index into [`PIR_ROTS`] (9 i32 per op), for H_ops→PIR order mapping
    pub(crate) _pir_rot_start: u32,
    /// Start index into [`SPIN_LG_OP_INDICES`] (0 if no data)
    pub(crate) _spin_lg_op_start: u32,
    /// Number of little-group operation indices
    pub(crate) _spin_lg_op_count: u8,
    /// Start index into [`SPIN_IMAG_CHARS`] (imaginary spinor characters).
    pub(crate) _spin_imag_start: u32,
    /// Number of imaginary spinor character values.
    pub(crate) _spin_imag_count: u16,
    /// Start index into [`CIR_COMPONENT_CHARS`] (0 if not compound)
    pub(crate) _cir_start: u32,
    /// Number of CIR components (0 for non-compound irreps, 2 for Z1Z4 type)
    pub(crate) _cir_count: u8,
    /// Number of operations per CIR component
    pub(crate) _cir_ops: u8,
    /// One-based index into [`generated_data::COMPOUND_METADATA`], or zero.
    pub(crate) _compound_metadata_index: u16,
}

impl IrrepRecord {
    /// Rational wave vector associated with this irrep.
    pub const fn k_vector(&self) -> KVector {
        KVector::new([self.kx, self.ky, self.kz], self.kd)
    }

    /// For spinor irreps: number of characters corresponding exactly to the
    /// indexed little-group operations. Returns 0 for scalar irreps.
    pub fn spin_lg_char_count(&self) -> usize {
        self._spin_lg_count as usize
    }

    /// Little-group operation indices for spinor irreps.
    /// Maps local character position → SG-local SPIN_OP index.
    pub fn spin_lg_op_indices(&self) -> &'static [u16] {
        if self._spin_lg_op_count == 0 {
            return &[];
        }
        let start = self._spin_lg_op_start as usize;
        let len = self._spin_lg_op_count as usize;
        &super::generated_data::SPIN_LG_OP_INDICES[start..start + len]
    }

    /// Imaginary parts of the spinor character table.
    pub fn spin_character_imag(&self) -> &'static [f64] {
        if self._spin_imag_count == 0 {
            return &[];
        }
        let start = self._spin_imag_start as usize;
        let len = self._spin_imag_count as usize;
        &super::generated_data::SPIN_IMAG_CHARS[start..start + len]
    }

    /// Spin symmetry operations with SU(2) lifts for any space group.
    ///
    /// This is a standalone version — does not require an `IrrepRecord`.
    /// Get the ISOTROPY setting (basis matrix + origin shift) for a space group.
    ///
    /// Returns `(basis_3x3_row_major, origin_3vec)` as f64 slices.
    /// Basis is always identity (same axes as ITA), origin has 205/230 non-trivial.
    pub fn sg_setting(sg: u8) -> (&'static [f64], &'static [f64]) {
        let idx = sg.saturating_sub(1) as usize;
        if idx >= 230 {
            return (&[], &[]);
        }
        let b_start = idx * 9;
        let o_start = idx * 3;
        (
            &super::generated_data::SG_SETTING_BASIS[b_start..b_start + 9],
            &super::generated_data::SG_SETTING_ORIGIN[o_start..o_start + 3],
        )
    }

    pub fn spin_ops_for_sg(sg: u8) -> (&'static [i32], &'static [f64], &'static [f64]) {
        let sg_idx = sg as usize;
        if sg_idx == 0 || sg_idx > 230 {
            return (&[], &[], &[]);
        }
        let (start, count) = super::generated_data::SPIN_OP_SG_INDEX[sg_idx];
        let start = start as usize;
        let count = count as usize;
        let rots = &super::generated_data::SPIN_OP_ROTS[start * 9..(start + count) * 9];
        let trans = &super::generated_data::SPIN_OP_TRANS[start * 3..(start + count) * 3];
        let su2 = &super::generated_data::SPIN_OP_SU2[start * 4..(start + count) * 4];
        (rots, trans, su2)
    }

    /// Spin symmetry operations with SU(2) lifts for this irrep's space group.
    ///
    /// Returns `(rotations, translations, pauli_su2)` slices where:
    /// - `rotations`: 9 i32 per op (3×3 rotation matrix, row-major)
    /// - `translations`: 3 f64 per op
    /// - `pauli_su2`: 4 f64 per op — **Pauli coefficients** `[u₀, u₁, u₂, u₃]`
    ///
    /// ## SU(2) Pauli coefficient convention
    ///
    /// The spin-½ representation matrix is reconstructed as:
    ///
    /// ```text
    /// U = u₀·I + i(u₁·σx + u₂·σy + u₃·σz)
    ///   = [[u₀ + iu₃,    u₂ + iu₁],
    ///      [-u₂ + iu₁,    u₀ - iu₃]]
    /// ```
    ///
    /// For crystallographic point groups, each uᵢ ∈ {0, ±½, ±1/√2, ±√3/2, ±1}
    /// and is stored as an exact f64 value (no floating-point noise).
    ///
    /// Composition follows quaternion multiplication:
    /// `su2_compose()` in `wigner.rs`.
    ///
    /// Verified by `scripts/test_su2_closure.py`: 229/229 SGs at 100% closure.
    pub fn spin_ops(&self) -> (&'static [i32], &'static [f64], &'static [f64]) {
        let sg_idx = self.sg as usize;
        if sg_idx == 0 || sg_idx > 230 {
            return (&[], &[], &[]);
        }
        let (start, count) = super::generated_data::SPIN_OP_SG_INDEX[sg_idx];
        let start = start as usize;
        let count = count as usize;
        let rots = &super::generated_data::SPIN_OP_ROTS[start * 9..(start + count) * 9];
        let trans = &super::generated_data::SPIN_OP_TRANS[start * 3..(start + count) * 3];
        let su2 = &super::generated_data::SPIN_OP_SU2[start * 4..(start + count) * 4];
        (rots, trans, su2)
    }

    /// Translation vectors for PIR operations, 3 f64 per op, same order as [`Self::characters`].
    ///
    /// Together with [`Self::pir_rotations`], enables full Seitz matching.
    pub fn pir_translations(&self) -> &'static [f64] {
        let char_count = self._char_count as usize;
        if char_count == 0 {
            return &[];
        }
        let start = (self._pir_rot_start as usize) / 9 * 3;
        let len = char_count * 3;
        let total = super::generated_data::PIR_TRANS.len();
        if start >= total || start + len > total {
            return &[];
        }
        &super::generated_data::PIR_TRANS[start..start + len]
    }

    /// Rotation matrices for PIR operations, 9 i32 per op, same order as [`Self::characters`].
    ///
    /// Used to build H_ops→PIR index mapping for the Wigner test.
    pub fn pir_rotations(&self) -> &'static [i32] {
        let char_count = self._char_count as usize;
        if char_count == 0 {
            return &[];
        }
        let start = self._pir_rot_start as usize;
        let len = char_count * 9;
        &super::generated_data::PIR_ROTS[start..start + len]
    }

    fn typed_pir_operations(&self) -> Result<Vec<SeitzOperation>, CharacterViewError> {
        let count = self._char_count as usize;
        if count == 0 {
            return Err(CharacterViewError::MissingData);
        }
        let rotation_start = self._pir_rot_start as usize;
        if rotation_start % 9 != 0 {
            return Err(CharacterViewError::MisalignedStorage);
        }
        let rotation_len = count
            .checked_mul(9)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let translation_start = rotation_start
            .checked_div(9)
            .and_then(|start| start.checked_mul(3))
            .ok_or(CharacterViewError::LengthMismatch)?;
        let translation_len = count
            .checked_mul(3)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let rotation_end = rotation_start
            .checked_add(rotation_len)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let translation_end = translation_start
            .checked_add(translation_len)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if rotation_end > super::generated_data::PIR_ROTS.len()
            || translation_end > super::generated_data::PIR_TRANS.len()
        {
            return Err(CharacterViewError::MissingData);
        }
        Ok((0..count)
            .map(|index| SeitzOperation {
                rotation: super::generated_data::PIR_ROTS
                    [rotation_start + index * 9..rotation_start + (index + 1) * 9]
                    .try_into()
                    .expect("checked PIR rotation width"),
                translation: super::generated_data::PIR_TRANS
                    [translation_start + index * 3..translation_start + (index + 1) * 3]
                    .try_into()
                    .expect("checked PIR translation width"),
            })
            .collect())
    }

    /// The physical scalar full-star/PIR character row in paired Seitz order.
    ///
    /// This is the existing real PIR row represented as complex values with
    /// zero imaginary parts. It is unavailable for spinors, whose character
    /// data include a distinct indexed little-group representation.
    pub fn scalar_full_star_row(&self) -> Result<CharacterRow, CharacterViewError> {
        if self.spinor {
            return Err(CharacterViewError::NotApplicable);
        }
        let start = self._char_start as usize;
        let end = start
            .checked_add(self._char_count as usize)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if end > super::generated_data::CHARACTERS.len() {
            return Err(CharacterViewError::MissingData);
        }
        let values = self
            .checked_characters()
            .iter()
            .copied()
            .map(|value| Complex64::new(value, 0.0))
            .collect();
        CharacterRow::from_parts(
            RepresentationSpaceKind::FullStarPir,
            values,
            self.typed_pir_operations()?,
            Some(self.dim as usize),
        )
    }

    /// Raw complex block traces of the first (stored-k) star-arm block.
    ///
    /// The slices span the complete PIR operation universe, not just the
    /// operations stabilizing the selected arm. Use the typed ordinary scalar
    /// block-trace view when operation/value pairing is required.
    pub(crate) fn raw_scalar_selected_arm_block_traces(&self) -> (&'static [f64], &'static [f64]) {
        if self.spinor || self._char_count == 0 {
            return (&[], &[]);
        }
        let start = self._pir_rot_start as usize / 9;
        let len = self._char_count as usize;
        let end = start + len;
        if end > super::generated_data::SCALAR_LITTLE_CHARS_VALID.len()
            || super::generated_data::SCALAR_LITTLE_CHARS_VALID[start..end].contains(&0)
        {
            return (&[], &[]);
        }
        (
            &super::generated_data::SCALAR_LITTLE_CHARS_REAL[start..end],
            &super::generated_data::SCALAR_LITTLE_CHARS_IMAG[start..end],
        )
    }

    /// Number of CIR (complex) components this PIR irrep decomposes into.
    /// 0 = non-compound, 2 = compound like Z1Z4 = Z1 ⊕ Z4.
    pub(crate) fn raw_cir_component_count(&self) -> usize {
        self._cir_count as usize
    }

    // Kept crate-visible for existing non-irrep summary code; external users
    // must use the semantic compound view rather than this raw count.
    pub(crate) fn cir_component_count(&self) -> usize {
        self.raw_cir_component_count()
    }

    /// Generation-time metadata for this compound irrep.
    ///
    /// The checked index and SG/label identity guard make malformed generated
    /// records return `None`; runtime code never parses labels or infers
    /// semantics from character values.
    pub fn compound_metadata(&self) -> Option<&'static CompoundMetadata> {
        if self._compound_metadata_index == 0 || self._cir_count == 0 {
            return None;
        }
        let metadata = super::generated_data::COMPOUND_METADATA
            .get(self._compound_metadata_index as usize - 1)?;
        (metadata.sg == self.sg && metadata.record_label == self.ml).then_some(metadata)
    }

    /// Compound character semantics from frozen generation metadata.
    pub fn compound_character_semantics(&self) -> Option<CompoundCharacterSemantics> {
        self.compound_metadata().map(|metadata| metadata.semantics)
    }

    /// Construct an ordinary scalar representation on one selected k arm.
    ///
    /// The selected dimension is the authoritative exact CIR header value
    /// frozen in this record. The aligned scalar little-group characters are
    /// checked against it at their unique full-Seitz identity operation.
    pub fn ordinary_scalar_selected_arm_block_trace(
        &self,
    ) -> Result<CharacterRow, CharacterViewError> {
        if self.spinor || self._cir_count > 0 {
            return Err(CharacterViewError::NotApplicable);
        }
        let expected_dimension = self._scalar_selected_dim as usize;
        if expected_dimension == 0 {
            return Err(CharacterViewError::MissingData);
        }
        let count = self._char_count as usize;
        if count == 0 {
            return Err(CharacterViewError::MissingData);
        }
        let rotation_start = self._pir_rot_start as usize;
        if rotation_start % 9 != 0 {
            return Err(CharacterViewError::MisalignedStorage);
        }
        let start = rotation_start / 9;
        let end = start
            .checked_add(count)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if end > super::generated_data::SCALAR_LITTLE_CHARS_REAL.len()
            || end > super::generated_data::SCALAR_LITTLE_CHARS_IMAG.len()
            || end > super::generated_data::SCALAR_LITTLE_CHARS_VALID.len()
        {
            return Err(CharacterViewError::MissingData);
        }
        if super::generated_data::SCALAR_LITTLE_CHARS_VALID[start..end].contains(&0) {
            return Err(CharacterViewError::MissingData);
        }
        let values = super::generated_data::SCALAR_LITTLE_CHARS_REAL[start..end]
            .iter()
            .zip(&super::generated_data::SCALAR_LITTLE_CHARS_IMAG[start..end])
            .map(|(&real, &imag)| Complex64::new(real, imag))
            .collect::<Vec<_>>();
        let row = CharacterRow::from_parts(
            RepresentationSpaceKind::SelectedArmBlockTrace,
            values,
            self.typed_pir_operations()?,
            Some(expected_dimension),
        )?;
        Ok(row)
    }

    fn checked_characters(&self) -> &'static [f64] {
        let start = self._char_start as usize;
        let end = start.saturating_add(self._char_count as usize);
        if end > super::generated_data::CHARACTERS.len() {
            &[]
        } else {
            &super::generated_data::CHARACTERS[start..end]
        }
    }

    /// Construct a compound selected-arm view whose public shape preserves
    /// realification versus distinct-component provenance.
    pub fn compound_selected_arm_view(
        &self,
    ) -> Result<CompoundSelectedArmCharacter, CharacterViewError> {
        if self.spinor || self._cir_count == 0 {
            return Err(CharacterViewError::NotApplicable);
        }
        let metadata = self
            .compound_metadata()
            .ok_or(CharacterViewError::MissingData)?;
        match metadata.semantics {
            CompoundCharacterSemantics::ConjugateRealification => {
                let seed = self.compound_constituent_character(0)?;
                let values = seed
                    .row
                    .values()
                    .iter()
                    .map(|value| *value + value.conj())
                    .collect();
                let block_trace = CharacterRow::from_parts(
                    RepresentationSpaceKind::SelectedArmBlockTrace,
                    values,
                    seed.row.operations().to_vec(),
                    Some(2 * metadata.cir_dimensions[0] as usize),
                )?;
                Ok(CompoundSelectedArmCharacter::ConjugateRealification { seed, block_trace })
            }
            CompoundCharacterSemantics::DistinctComponentSum => {
                let first = self.compound_constituent_character(0)?;
                let second = self.compound_constituent_character(1)?;
                let values = first
                    .row
                    .values()
                    .iter()
                    .zip(second.row.values())
                    .map(|(first, second)| *first + *second)
                    .collect();
                let block_trace = CharacterRow::from_parts(
                    RepresentationSpaceKind::SelectedArmBlockTrace,
                    values,
                    first.row.operations().to_vec(),
                    Some(
                        metadata
                            .cir_dimensions
                            .iter()
                            .map(|&dimension| dimension as usize)
                            .sum(),
                    ),
                )?;
                Ok(CompoundSelectedArmCharacter::DistinctComponentSum {
                    first,
                    second,
                    block_trace,
                })
            }
        }
    }

    /// Internal constructor for the constituent rows embedded by
    /// [`Self::compound_selected_arm_view`]. CIR rotations are checked
    /// against the final PIR rotations here; the corresponding translations
    /// are guaranteed by the generator's post-padding assertion.
    fn compound_constituent_character(
        &self,
        component: usize,
    ) -> Result<CompoundConstituentCharacter, CharacterViewError> {
        let metadata = self
            .compound_metadata()
            .ok_or(CharacterViewError::NotApplicable)?;
        if component >= 2 || self._cir_count != 2 {
            return Err(CharacterViewError::InvalidComponent);
        }
        if self._cir_start % 2 != 0 {
            return Err(CharacterViewError::MisalignedStorage);
        }
        let operations = self.typed_pir_operations()?;
        if operations.len() != self._cir_ops as usize {
            return Err(CharacterViewError::OperationOrderMismatch);
        }
        let expected_values = operations
            .len()
            .checked_mul(2)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let raw_start = (self._cir_start as usize)
            .checked_add(
                component
                    .checked_mul(self._cir_ops as usize)
                    .and_then(|value| value.checked_mul(2))
                    .ok_or(CharacterViewError::LengthMismatch)?,
            )
            .ok_or(CharacterViewError::LengthMismatch)?;
        let raw_end = raw_start
            .checked_add(expected_values)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if raw_end > super::generated_data::CIR_COMPONENT_CHARS.len() {
            return Err(CharacterViewError::LengthMismatch);
        }
        let raw = &super::generated_data::CIR_COMPONENT_CHARS[raw_start..raw_end];
        let cir_rotation_start = (self._cir_start as usize / 2)
            .checked_mul(9)
            .and_then(|value| {
                component
                    .checked_mul(operations.len())
                    .and_then(|offset| offset.checked_mul(9))
                    .and_then(|offset| value.checked_add(offset))
            })
            .ok_or(CharacterViewError::LengthMismatch)?;
        let cir_rotation_end = cir_rotation_start
            .checked_add(operations.len() * 9)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if cir_rotation_end > super::generated_data::CIR_ROTS.len()
            || operations.iter().enumerate().any(|(index, operation)| {
                super::generated_data::CIR_ROTS
                    [cir_rotation_start + index * 9..cir_rotation_start + (index + 1) * 9]
                    != operation.rotation
            })
        {
            return Err(CharacterViewError::OperationOrderMismatch);
        }
        let values = raw
            .chunks_exact(2)
            .map(|pair| Complex64::new(pair[0], pair[1]))
            .collect();
        let row = CharacterRow::from_parts(
            RepresentationSpaceKind::ConstituentCir,
            values,
            operations,
            Some(metadata.cir_dimensions[component] as usize),
        )?;
        Ok(CompoundConstituentCharacter {
            component,
            irnumber: metadata.cir_irnumbers[component],
            label: metadata.cir_labels[component],
            dimension: metadata.cir_dimensions[component] as usize,
            row,
        })
    }

    /// Construct the spinor selected-arm row with indexed SU(2) operations.
    pub fn spinor_selected_arm_view(&self) -> Result<SpinCharacterRow, CharacterViewError> {
        if !self.spinor {
            return Err(CharacterViewError::NotApplicable);
        }
        let count = self._spin_lg_count as usize;
        if count == 0
            || self._char_count as usize != count
            || self._spin_lg_op_count as usize != count
            || self._spin_imag_count as usize != count
        {
            return Err(CharacterViewError::MissingData);
        }
        let char_start = self._char_start as usize;
        let char_end = char_start
            .checked_add(count)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let imag_start = self._spin_imag_start as usize;
        let imag_end = imag_start
            .checked_add(count)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if char_end > super::generated_data::CHARACTERS.len()
            || imag_end > super::generated_data::SPIN_IMAG_CHARS.len()
        {
            return Err(CharacterViewError::MissingData);
        }
        let index_start = self._spin_lg_op_start as usize;
        let index_end = index_start
            .checked_add(count)
            .ok_or(CharacterViewError::LengthMismatch)?;
        if index_end > super::generated_data::SPIN_LG_OP_INDICES.len() {
            return Err(CharacterViewError::InvalidOperationIndex);
        }
        let sg_index = self.sg as usize;
        if sg_index == 0 || sg_index > 230 {
            return Err(CharacterViewError::MissingData);
        }
        let (spin_start, spin_count) = super::generated_data::SPIN_OP_SG_INDEX[sg_index];
        let spin_start = spin_start as usize;
        let spin_count = spin_count as usize;
        let rotation_start = spin_start
            .checked_mul(9)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let rotation_end = spin_start
            .checked_add(spin_count)
            .and_then(|end| end.checked_mul(9))
            .ok_or(CharacterViewError::LengthMismatch)?;
        let translation_start = spin_start
            .checked_mul(3)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let translation_end = spin_start
            .checked_add(spin_count)
            .and_then(|end| end.checked_mul(3))
            .ok_or(CharacterViewError::LengthMismatch)?;
        if rotation_end > super::generated_data::SPIN_OP_ROTS.len()
            || translation_end > super::generated_data::SPIN_OP_TRANS.len()
        {
            return Err(CharacterViewError::MissingData);
        }
        let spin_rotations = &super::generated_data::SPIN_OP_ROTS[rotation_start..rotation_end];
        let spin_translations =
            &super::generated_data::SPIN_OP_TRANS[translation_start..translation_end];
        let su2_start = spin_start
            .checked_mul(4)
            .ok_or(CharacterViewError::LengthMismatch)?;
        let su2_end = spin_start
            .checked_add(spin_count)
            .and_then(|end| end.checked_mul(4))
            .ok_or(CharacterViewError::LengthMismatch)?;
        if su2_end > super::generated_data::SPIN_OP_SU2.len() {
            return Err(CharacterViewError::MissingData);
        }
        let spin_su2 = &super::generated_data::SPIN_OP_SU2[su2_start..su2_end];
        if spin_rotations.len() != spin_count * 9
            || spin_translations.len() != spin_count * 3
            || spin_su2.len() != spin_count * 4
        {
            return Err(CharacterViewError::LengthMismatch);
        }
        validate_spinor_operation_indices(
            &super::generated_data::SPIN_LG_OP_INDICES[index_start..index_end],
            spin_count,
        )?;
        let mut values = Vec::with_capacity(count);
        let mut operations = Vec::with_capacity(count);
        for index in index_start..index_end {
            let spin_index = super::generated_data::SPIN_LG_OP_INDICES[index] as usize;
            if spin_index >= spin_count {
                return Err(CharacterViewError::InvalidOperationIndex);
            }
            let value_index = char_start + (index - index_start);
            let imag_index = imag_start + (index - index_start);
            values.push(Complex64::new(
                super::generated_data::CHARACTERS[value_index],
                super::generated_data::SPIN_IMAG_CHARS[imag_index],
            ));
            operations.push(SpinSeitzOperation {
                seitz: SeitzOperation {
                    rotation: spin_rotations[spin_index * 9..(spin_index + 1) * 9]
                        .try_into()
                        .expect("checked spin rotation width"),
                    translation: spin_translations[spin_index * 3..(spin_index + 1) * 3]
                        .try_into()
                        .expect("checked spin translation width"),
                },
                su2: spin_su2[spin_index * 4..(spin_index + 1) * 4]
                    .try_into()
                    .expect("checked spin SU(2) width"),
            });
        }
        SpinCharacterRow::from_parts(values, operations, self.dim as usize)
    }

    /// Complex character table for a specific CIR component.
    ///
    /// Returns `(re, im)` pairs in the generated data-Hall operation order.
    /// [`Self::raw_cir_rotations`] contains the corresponding rotations.
    pub(crate) fn raw_cir_component_chars(&self, comp: usize) -> &'static [f64] {
        if comp >= self._cir_count as usize {
            return &[];
        }
        let start = self._cir_start as usize + comp * self._cir_ops as usize * 2;
        let len = self._cir_ops as usize * 2;
        &super::generated_data::CIR_COMPONENT_CHARS[start..start + len]
    }

    /// Rotation matrices for CIR operations of a specific component.
    ///
    /// Returns 9×n_ops i32 values (r00,r01,r02, r10,r11,r12, r20,r21,r22 per op).
    pub(crate) fn raw_cir_rotations(&self, comp: usize) -> &'static [i32] {
        if comp >= self._cir_count as usize {
            return &[];
        }
        // _cir_start indexes into CIR_COMPONENT_CHARS (2 f64 per op).
        // CIR_ROTS has 9 i32 per op in the SAME flat order.
        // Convert from f64-pair index to rotation index:
        let ops_before = self._cir_start as usize / 2;
        let start = (ops_before + comp * self._cir_ops as usize) * 9;
        let len = self._cir_ops as usize * 9;
        &super::generated_data::CIR_ROTS[start..start + len]
    }
}

impl IrrepRecord {
    /// Full-star/PIR character table: χ(g) = Tr(D(g)) for each stored
    /// space-group operator.
    ///
    /// The character χ(g) of a representation D is the trace of the
    /// representation matrix for each symmetry operation g.  The return
    /// slice is the stored PIR operation universe and each entry is a
    /// floating-point value (possibly negative, fractional, or zero); its
    /// length is not a selected-arm little-group contract. For operation-aware
    /// typed views, use the family-specific constructors.
    pub fn characters(&self) -> &'static [f64] {
        if self._char_count == 0 {
            return &[];
        }
        &self::generated_data::CHARACTERS
            [self._char_start as usize..(self._char_start as usize + self._char_count as usize)]
    }

    /// Full irrep matrices for each operator, flattened: op0(row0,row1,...), op1(...), ...
    ///
    /// **Order**: ISOTROPY (PIR_data.txt) order — NOT spglib H_ops order.
    /// Use [`Self::matrices_reordered`] to reorder to spglib H_ops order.
    ///
    /// The total number of elements is `opcount × dim²`.
    /// For operator `g`, the matrix D(g) is at offset `g × dim²` with
    /// row-major layout: D[0][0], D[0][1], ..., D[dim-1][dim-1].
    pub fn matrices(&self) -> &'static [f64] {
        if self._mat_count == 0 {
            return &[];
        }
        &self::generated_data::MATRICES
            [self._mat_start as usize..(self._mat_start + self._mat_count) as usize]
    }

    /// Full irrep matrices reordered to match spglib H_ops order.
    ///
    /// Only those H_ops that match a PIR operation (via rotation matrix)
    /// get matrix data.  Unmatched ops get zero-filled blocks.
    /// Returns an empty `Vec` if no matrix data are stored. Returns an error
    /// instead of silently treating PIR order as H_ops order when mapping fails.
    pub fn matrices_reordered(
        &self,
        h_seitz: &[crate::irrep::wigner::SeitzOp],
    ) -> Result<Vec<f64>, MatrixReorderError> {
        let mats = self.matrices();
        let rots = self.pir_rotations();
        if mats.is_empty() {
            return Ok(Vec::new());
        }
        if rots.is_empty() {
            return Err(MatrixReorderError::OperationMappingFailed);
        }
        let dim = self.dim as usize;
        let n_pir_ops = self._char_count as usize;
        let block_size = dim * dim;

        // Build partial H_ops → PIR map (only for ops in the little group)
        let h_to_pir = match crate::irrep::wigner::build_h_to_cir_map(h_seitz, rots) {
            Some(m) => m,
            None => {
                // Full mapping failed — try matching only ops present in PIR
                let n_cir = rots.len() / 9;
                let h_count = h_seitz.len().min(n_cir);
                if h_count == 0 {
                    return Err(MatrixReorderError::OperationMappingFailed);
                }
                crate::irrep::wigner::build_h_to_cir_map(&h_seitz[..h_count], rots)
                    .ok_or(MatrixReorderError::OperationMappingFailed)?
            }
        };

        let mut reordered = vec![0.0f64; mats.len()];
        for (h_idx, &pir_idx) in h_to_pir.iter().take(n_pir_ops).enumerate() {
            if pir_idx >= n_pir_ops {
                continue;
            }
            let src_start = pir_idx * block_size;
            let dst_start = h_idx * block_size;
            if src_start + block_size <= mats.len() && dst_start + block_size <= reordered.len() {
                reordered[dst_start..dst_start + block_size]
                    .copy_from_slice(&mats[src_start..src_start + block_size]);
            }
        }
        Ok(reordered)
    }

    /// Isotropy subgroups for this irrep — no index arithmetic needed.
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::irrep::query::irreps_of;
    ///
    /// for ir in irreps_of(221) {
    ///     if ir.ml == "GM4-" {
    ///         for sub in ir.subgroups() {
    ///             println!("#{} {}", sub.sg, sub.symbol);
    ///         }
    ///     }
    /// }
    /// ```
    pub fn subgroups(&self) -> &'static [IsotropyRecord] {
        &self::generated_data::ISOTROPY_SUBGROUPS
            [self._iso_start as usize..(self._iso_start + self._iso_count) as usize]
    }

    /// Magnetic isotropy subgroups for this irrep.
    ///
    /// When the order parameter of this irrep condenses, the system
    /// can lower its symmetry to one of these magnetic space groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::irrep::query::irreps_of;
    ///
    /// for ir in irreps_of(221) {
    ///     if ir.ml == "GM4-" {
    ///         for sub in ir.magnetic_subgroups() {
    ///             println!("{} {}", sub.bns_label, sub.direction);
    ///         }
    ///     }
    /// }
    /// ```
    pub fn magnetic_subgroups(&self) -> &'static [MagneticIsotropyRecord] {
        if self._mag_iso_count == 0 {
            return &[];
        }
        &self::generated_data::MAGNETIC_ISOTROPY_SUBGROUPS
            [self._mag_iso_start as usize..(self._mag_iso_start + self._mag_iso_count) as usize]
    }

    /// k-point label prefix extracted from the ML label.
    ///
    /// - `"GM4+"` → `"GM"` (Γ point)
    /// - `"X3-"` → `"X"`
    /// - `"DT1"` → `"DT"` (Δ line)
    pub fn k_label(&self) -> &'static str {
        let body = self.ml.trim_end_matches(['+', '-']);
        let end = body
            .find(|c: char| c.is_ascii_digit())
            .unwrap_or(body.len());
        &body[..end]
    }

    /// Whether this is a special k-point (not a line or plane).
    pub fn is_point(&self) -> bool {
        k_label_is_point(self.k_label())
    }
}

fn k_label_is_point(label: &str) -> bool {
    // Generated Miller–Love point labels are single-letter labels plus GM
    // for Gamma. Two-letter labels such as DT, LD, and SM denote lines, while
    // GP denotes a general-position manifold.
    label == "GM" || label.len() == 1
}

#[cfg(test)]
mod point_label_tests {
    use super::k_label_is_point;

    #[test]
    fn two_letter_line_labels_are_not_points() {
        assert!(k_label_is_point("GM"));
        assert!(k_label_is_point("X"));
        for label in ["DT", "LD", "SM", "GP"] {
            assert!(!k_label_is_point(label), "{label}");
        }
    }
}

impl IsotropyRecord {
    /// Human-readable one-line description.
    pub fn describe(&self) -> String {
        format!(
            "#{} {} ({}), domains={}, arms={}",
            self.sg, self.symbol, self.schoenflies, self.domains, self.arms
        )
    }
}

impl std::fmt::Display for IsotropyRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "#{} {} dir={} domains={} arms={}",
            self.sg, self.symbol, self.direction, self.domains, self.arms
        )
    }
}

impl std::fmt::Display for MagneticIsotropyRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "UNI {} ({}) dir={}",
            self.mag_sg, self.bns_label, self.direction
        )
    }
}

/// Compact isotropy subgroup record for the generated flat array.
#[derive(Debug, Clone, Copy)]
pub struct IsotropyRecord {
    /// Subgroup space group number (1–230)
    pub sg: usize,
    /// Hermann-Mauguin symbol
    pub symbol: &'static str,
    /// Schoenflies symbol
    pub schoenflies: &'static str,
    /// Order-parameter direction label
    pub direction: &'static str,
    /// Number of domains
    pub domains: usize,
    /// Number of arms in the star
    pub arms: usize,
}

/// A magnetic isotropy subgroup: the lower-symmetry magnetic space group
/// obtained when the order parameter condenses along a specific direction
/// for a given non-magnetic irrep.  Magnetic isotropy subgroups describe
/// the possible magnetic structures that can form.
#[derive(Debug, Clone, Copy)]
pub struct MagneticIsotropyRecord {
    /// Magnetic space group UNI number (1–1651)
    pub mag_sg: usize,
    /// BNS (Belov-Neronova-Smirnova) symbol, e.g. `"Pm'mm"`
    pub bns_label: &'static str,
    /// ISOTROPY label, e.g. `"47.251"`
    pub iso_label: &'static str,
    /// Order-parameter direction label
    pub direction: &'static str,
}

/// Auto-generated data from iso_data files.
///
/// This module is regenerated by `scripts/generate_irrep_data.py`.
/// Do not edit manually.
#[allow(missing_docs)]
pub mod generated_data {
    #![allow(clippy::all)]
    include!("generated_data.rs");
}

#[cfg(test)]
mod compound_metadata_tests {
    use num_complex::Complex64;

    use super::CompoundCharacterSemantics::{ConjugateRealification, DistinctComponentSum};

    #[test]
    fn generated_compound_metadata_is_complete_and_indexed() {
        assert_eq!(super::generated_data::COMPOUND_METADATA.len(), 672);
        let mut seen = 0;
        let mut realified = 0;
        let mut distinct = 0;
        for sg in 1..=230 {
            for irrep in crate::irrep::query::irreps_of(sg) {
                let metadata = irrep.compound_metadata();
                if irrep._cir_count == 0 {
                    assert!(metadata.is_none(), "noncompound SG{sg} {}", irrep.ml);
                    continue;
                }
                let metadata = metadata.expect("compound metadata");
                seen += 1;
                assert_eq!(metadata.sg, irrep.sg);
                assert_eq!(metadata.record_label, irrep.ml);
                assert_eq!(
                    metadata.naming_grammar_version,
                    super::COMPOUND_NAMING_GRAMMAR_VERSION
                );
                assert_eq!(metadata.provenance, super::COMPOUND_NAMING_PROVENANCE);
                assert!(
                    metadata
                        .cir_dimensions
                        .iter()
                        .all(|dimension| *dimension > 0)
                );
                match metadata.semantics {
                    ConjugateRealification => {
                        realified += 1;
                        assert_eq!(metadata.cir_irnumbers[0], metadata.cir_irnumbers[1]);
                        assert_eq!(metadata.cir_labels[0], metadata.cir_labels[1]);
                    }
                    DistinctComponentSum => {
                        distinct += 1;
                        assert_ne!(metadata.cir_irnumbers[0], metadata.cir_irnumbers[1]);
                        assert_ne!(metadata.cir_labels[0], metadata.cir_labels[1]);
                    }
                }
            }
        }
        assert_eq!(seen, 672);
        assert_eq!(realified, 153);
        assert_eq!(distinct, 519);
    }

    #[test]
    fn representative_compound_metadata_semantics() {
        for (sg, label, expected) in [
            (199, "P2P2", ConjugateRealification),
            (220, "P1P1", ConjugateRealification),
            (220, "P2P2", ConjugateRealification),
            (220, "P3P3", ConjugateRealification),
            (182, "H1H2", DistinctComponentSum),
            (46, "W1W2", DistinctComponentSum),
        ] {
            let irrep = crate::irrep::query::irreps_of(sg)
                .iter()
                .find(|irrep| irrep.ml == label)
                .expect("compound record");
            assert_eq!(irrep.compound_character_semantics(), Some(expected));
        }
    }

    #[test]
    fn noncompound_and_spinor_records_have_no_compound_metadata() {
        let scalar = crate::irrep::query::irreps_of(221)
            .iter()
            .find(|irrep| !irrep.spinor && irrep.ml == "GM1+")
            .expect("ordinary scalar record");
        assert!(scalar.compound_metadata().is_none());
        let spinor = crate::irrep::query::irreps_of(221)
            .iter()
            .find(|irrep| irrep.spinor)
            .expect("spinor record");
        assert!(spinor.compound_metadata().is_none());
    }

    #[test]
    fn selected_arm_character_census() {
        let mut scalar = 0;
        let mut ordinary = 0;
        let mut compound = 0;
        let mut spinor = 0;
        for sg in 1..=230 {
            for irrep in crate::irrep::query::irreps_of(sg) {
                if irrep.spinor {
                    let row = irrep.spinor_selected_arm_view().expect("spinor row");
                    assert_eq!(row.len(), row.operations().len());
                    assert!(row.dimension() > 0);
                    let identities = row
                        .operations()
                        .iter()
                        .enumerate()
                        .filter(|(_, operation)| {
                            super::is_identity_seitz(&operation.seitz)
                                && operation.su2 == [1.0, 0.0, 0.0, 0.0]
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(identities.len(), 1);
                    let identity = row.entry(identities[0].0).unwrap().0;
                    assert!((identity.re - row.dimension() as f64).abs() < 1e-10);
                    assert!(identity.im.abs() < 1e-10);
                    spinor += 1;
                } else if irrep.cir_component_count() > 0 {
                    assert_eq!(irrep._scalar_selected_dim, 0);
                    compound += 1;
                    let view = irrep.compound_selected_arm_view().expect("compound row");
                    let row = view.block_trace();
                    assert_eq!(row.len(), row.operations().len());
                    assert!(row.dimension() > 0);
                    let identities = row
                        .operations()
                        .iter()
                        .enumerate()
                        .filter(|(_, operation)| super::is_identity_seitz(operation))
                        .collect::<Vec<_>>();
                    assert_eq!(identities.len(), 1);
                    let identity = row.get(identities[0].0).unwrap();
                    assert!((identity.re - row.dimension() as f64).abs() < 1e-10);
                    assert!(identity.im.abs() < 1e-10);
                    scalar += 1;
                } else {
                    assert!(irrep._scalar_selected_dim > 0);
                    ordinary += 1;
                    let row = irrep
                        .ordinary_scalar_selected_arm_block_trace()
                        .expect("ordinary scalar row");
                    assert_eq!(row.len(), row.operations().len());
                    assert_eq!(row.dimension(), irrep._scalar_selected_dim as usize);
                    let identities = row
                        .operations()
                        .iter()
                        .enumerate()
                        .filter(|(_, operation)| super::is_identity_seitz(operation))
                        .collect::<Vec<_>>();
                    assert_eq!(identities.len(), 1);
                    let identity = row.get(identities[0].0).unwrap();
                    assert!((identity.re - row.dimension() as f64).abs() < 1e-10);
                    assert!(identity.im.abs() < 1e-10);
                    scalar += 1;
                }
            }
        }
        assert_eq!(scalar, 4777);
        assert_eq!(ordinary, 4105);
        assert_eq!(compound, 672);
        assert_eq!(spinor, 3611);
    }

    #[test]
    fn scalar_full_star_character_census() {
        let mut checked = 0;
        for sg in 1..=230 {
            for irrep in crate::irrep::query::irreps_of(sg)
                .iter()
                .filter(|irrep| !irrep.spinor)
            {
                let row = irrep.scalar_full_star_row().expect("scalar full-star row");
                assert_eq!(row.len(), row.operations().len());
                assert_eq!(row.dimension(), irrep.dim as usize);
                let identities = row
                    .operations()
                    .iter()
                    .enumerate()
                    .filter(|(_, operation)| super::is_identity_seitz(operation))
                    .collect::<Vec<_>>();
                assert_eq!(identities.len(), 1);
                let identity = row.entry(identities[0].0).unwrap().0;
                assert_eq!(identity.im, 0.0);
                assert_eq!(identity.re, irrep.dim as f64);
                checked += 1;
            }
        }
        assert_eq!(checked, 4777);
    }

    #[test]
    fn typed_character_family_applicability_is_explicit() {
        let ordinary = crate::irrep::query::irreps_of(1)
            .iter()
            .find(|irrep| !irrep.spinor)
            .expect("ordinary scalar");
        let compound = crate::irrep::query::irreps_of(23)
            .iter()
            .find(|irrep| irrep.raw_cir_component_count() > 0)
            .expect("compound scalar");
        let spinor = crate::irrep::query::irreps_of(1)
            .iter()
            .find(|irrep| irrep.spinor)
            .expect("spinor");
        assert!(matches!(
            ordinary.compound_selected_arm_view(),
            Err(super::CharacterViewError::NotApplicable)
        ));
        assert!(matches!(
            ordinary.spinor_selected_arm_view(),
            Err(super::CharacterViewError::NotApplicable)
        ));
        assert!(matches!(
            compound.ordinary_scalar_selected_arm_block_trace(),
            Err(super::CharacterViewError::NotApplicable)
        ));
        assert!(matches!(
            compound.spinor_selected_arm_view(),
            Err(super::CharacterViewError::NotApplicable)
        ));
        assert!(matches!(
            spinor.scalar_full_star_row(),
            Err(super::CharacterViewError::NotApplicable)
        ));
        assert!(matches!(
            spinor.ordinary_scalar_selected_arm_block_trace(),
            Err(super::CharacterViewError::NotApplicable)
        ));
        assert!(matches!(
            spinor.compound_selected_arm_view(),
            Err(super::CharacterViewError::NotApplicable)
        ));
    }

    #[test]
    fn typed_compound_rows_follow_metadata_semantics() {
        let mut compounds = 0;
        let mut realifications = 0;
        let mut distinct = 0;
        let mut references = 0;
        for sg in 1..=230 {
            for irrep in crate::irrep::query::irreps_of(sg) {
                if irrep.cir_component_count() == 0 {
                    continue;
                }
                compounds += 1;
                let metadata = irrep.compound_metadata().expect("compound metadata");
                let view = irrep.compound_selected_arm_view().expect("compound row");
                let (first, second, selected) = match &view {
                    super::CompoundSelectedArmCharacter::ConjugateRealification {
                        seed,
                        block_trace,
                    } => (seed, None, block_trace),
                    super::CompoundSelectedArmCharacter::DistinctComponentSum {
                        first,
                        second,
                        block_trace,
                    } => (first, Some(second), block_trace),
                };
                assert_eq!(first.row.len(), selected.len());
                assert_eq!(first.dimension, metadata.cir_dimensions[0] as usize);
                if let Some(second) = second {
                    assert_eq!(first.row.operations(), second.row.operations());
                    assert_eq!(second.dimension, metadata.cir_dimensions[1] as usize);
                }
                references += 2;
                assert_eq!(
                    selected.representation_space(),
                    super::RepresentationSpaceKind::SelectedArmBlockTrace
                );
                assert_eq!(selected.len(), first.row.len());
                for index in 0..selected.len() {
                    let expected = match metadata.semantics {
                        ConjugateRealification => {
                            first.row.get(index).unwrap() + first.row.get(index).unwrap().conj()
                        }
                        DistinctComponentSum => {
                            first.row.get(index).unwrap()
                                + second.expect("distinct second").row.get(index).unwrap()
                        }
                    };
                    assert!((selected.get(index).unwrap() - expected).norm() < 1e-10);
                }
                match metadata.semantics {
                    ConjugateRealification => realifications += 1,
                    DistinctComponentSum => distinct += 1,
                }
            }
        }
        assert_eq!(
            (compounds, realifications, distinct, references),
            (672, 153, 519, 1344)
        );
    }

    #[test]
    fn representative_typed_dimensions_and_realification_regression() {
        let sg76 = crate::irrep::query::irreps_of(76)
            .iter()
            .find(|irrep| irrep.ml == "R1R2")
            .expect("SG76 R1R2");
        assert_eq!(sg76.scalar_full_star_row().unwrap().dimension(), 4);
        assert_eq!(
            sg76.compound_selected_arm_view()
                .unwrap()
                .block_trace()
                .dimension(),
            2
        );

        let mut literal_sum_differs = false;
        let mut sg23_stored_differs = false;
        for (sg, label) in [
            (23, "W1W1"),
            (199, "P2P2"),
            (220, "P1P1"),
            (220, "P2P2"),
            (220, "P3P3"),
        ] {
            let irrep = crate::irrep::query::irreps_of(sg)
                .iter()
                .find(|irrep| irrep.ml == label)
                .expect("realification record");
            let view = irrep.compound_selected_arm_view().unwrap();
            let (first, selected) = match &view {
                super::CompoundSelectedArmCharacter::ConjugateRealification {
                    seed,
                    block_trace,
                } => (seed, block_trace),
                super::CompoundSelectedArmCharacter::DistinctComponentSum {
                    first,
                    block_trace,
                    ..
                } => (first, block_trace),
            };
            let stored_second = irrep.compound_constituent_character(1).unwrap();
            for index in 0..selected.len() {
                let expected = first.row.get(index).unwrap() + first.row.get(index).unwrap().conj();
                assert!((selected.get(index).unwrap() - expected).norm() < 1e-10);
                let literal = first.row.get(index).unwrap() + stored_second.row.get(index).unwrap();
                literal_sum_differs |= (literal - expected).norm() > 1e-10;
                if sg == 23 {
                    sg23_stored_differs |= (stored_second.row.get(index).unwrap()
                        - first.row.get(index).unwrap().conj())
                    .norm()
                        > 1e-10;
                }
            }
        }
        assert!(literal_sum_differs);
        assert!(sg23_stored_differs);
    }

    #[test]
    fn representative_distinct_sum_preserves_complex_phase() {
        let mut saw_complex = false;
        for (sg, label) in [(182, "H1H2"), (46, "W1W2")] {
            let irrep = crate::irrep::query::irreps_of(sg)
                .iter()
                .find(|irrep| irrep.ml == label)
                .expect("distinct-sum record");
            let view = irrep.compound_selected_arm_view().unwrap();
            let (first, second, selected) = match &view {
                super::CompoundSelectedArmCharacter::DistinctComponentSum {
                    first,
                    second,
                    block_trace,
                } => (first, second, block_trace),
                _ => panic!("expected distinct component sum"),
            };
            for index in 0..selected.len() {
                let expected = first.row.get(index).unwrap() + second.row.get(index).unwrap();
                assert!((selected.get(index).unwrap() - expected).norm() < 1e-10);
                saw_complex |= expected.im.abs() > 1e-10;
            }
        }
        assert!(saw_complex);
    }

    #[test]
    fn typed_spinor_rows_use_indexed_little_group_operations() {
        let mut found_nontrivial = false;
        for sg in 1..=230 {
            for irrep in crate::irrep::query::irreps_of(sg)
                .iter()
                .filter(|irrep| irrep.spinor)
            {
                let indices = irrep.spin_lg_op_indices();
                if !indices
                    .iter()
                    .enumerate()
                    .any(|(position, &index)| index as usize != position)
                {
                    continue;
                }
                let row = irrep.spinor_selected_arm_view().unwrap();
                let (rotations, translations, _) = irrep.spin_ops();
                for (position, &spin_index) in indices.iter().enumerate() {
                    let spin_index = spin_index as usize;
                    assert_eq!(
                        row.entry(position).unwrap().1.seitz.rotation,
                        rotations[spin_index * 9..(spin_index + 1) * 9]
                    );
                    assert_eq!(
                        row.entry(position).unwrap().1.seitz.translation,
                        translations[spin_index * 3..(spin_index + 1) * 3]
                    );
                }
                found_nontrivial = true;
                break;
            }
            if found_nontrivial {
                break;
            }
        }
        assert!(found_nontrivial);
    }

    #[test]
    fn typed_character_helpers_find_identity_without_first_column_assumption() {
        let identity = super::SeitzOperation {
            rotation: [1, 0, 0, 0, 1, 0, 0, 0, 1],
            translation: [1.0, 0.0, -2.0],
        };
        let nonidentity = super::SeitzOperation {
            rotation: [0, -1, 0, 1, 0, 0, 0, 0, 1],
            translation: [0.0, 0.0, 0.0],
        };
        let row = super::CharacterRow::from_parts(
            super::RepresentationSpaceKind::SelectedArmBlockTrace,
            vec![Complex64::new(0.0, 0.0), Complex64::new(2.0, 0.0)],
            vec![nonidentity, identity],
            Some(2),
        )
        .unwrap();
        assert_eq!(row.dimension(), 2);
        assert_eq!(row.get(1), Some(Complex64::new(2.0, 0.0)));
        assert!(matches!(
            super::CharacterRow::from_parts(
                super::RepresentationSpaceKind::SelectedArmBlockTrace,
                vec![Complex64::new(1.0, 0.0)],
                vec![identity],
                Some(2),
            ),
            Err(super::CharacterViewError::DimensionMismatch)
        ));
        assert!(matches!(
            super::CharacterRow::from_parts(
                super::RepresentationSpaceKind::SelectedArmBlockTrace,
                vec![Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0)],
                vec![identity, identity],
                None,
            ),
            Err(super::CharacterViewError::NonUniqueIdentity)
        ));
        assert!(matches!(
            super::CharacterRow::from_parts(
                super::RepresentationSpaceKind::SelectedArmBlockTrace,
                vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
                vec![identity],
                None,
            ),
            Err(super::CharacterViewError::LengthMismatch)
        ));
        let spin_identity = super::SpinSeitzOperation {
            seitz: identity,
            su2: [1.0, 0.0, 0.0, 0.0],
        };
        let spin_row = super::SpinCharacterRow::from_parts(
            vec![Complex64::new(2.0, 0.0)],
            vec![spin_identity],
            1,
        );
        assert!(matches!(
            spin_row,
            Err(super::CharacterViewError::DimensionMismatch)
        ));
        let spin_minus_identity = super::SpinSeitzOperation {
            seitz: identity,
            su2: [-1.0, 0.0, 0.0, 0.0],
        };
        assert!(matches!(
            super::SpinCharacterRow::from_parts(
                vec![Complex64::new(2.0, 0.0)],
                vec![spin_minus_identity],
                2,
            ),
            Err(super::CharacterViewError::NonUniqueIdentity)
        ));
        assert!(matches!(
            super::SpinCharacterRow::from_parts(
                vec![Complex64::new(2.0, 0.0), Complex64::new(0.0, 0.0)],
                vec![spin_identity],
                2,
            ),
            Err(super::CharacterViewError::LengthMismatch)
        ));
        assert!(matches!(
            super::validate_spinor_operation_indices(&[0, 0], 1),
            Err(super::CharacterViewError::DuplicateOperationIndex)
        ));
        assert!(matches!(
            super::validate_spinor_operation_indices(&[1], 1),
            Err(super::CharacterViewError::InvalidOperationIndex)
        ));
    }

    #[test]
    fn typed_character_rows_reject_nonfinite_data_and_bad_offsets() {
        let identity = super::SeitzOperation {
            rotation: [1, 0, 0, 0, 1, 0, 0, 0, 1],
            translation: [0.0, 0.0, 0.0],
        };
        let nonidentity = super::SeitzOperation {
            rotation: [0, -1, 0, 1, 0, 0, 0, 0, 1],
            translation: [0.0, 0.0, 0.0],
        };
        assert!(matches!(
            super::CharacterRow::from_parts(
                super::RepresentationSpaceKind::FullStarPir,
                vec![Complex64::new(f64::NAN, 0.0), Complex64::new(1.0, 0.0)],
                vec![identity, nonidentity],
                None,
            ),
            Err(super::CharacterViewError::NonFiniteData)
        ));
        assert!(matches!(
            super::CharacterRow::from_parts(
                super::RepresentationSpaceKind::FullStarPir,
                vec![Complex64::new(1.0, 0.0)],
                vec![super::SeitzOperation {
                    translation: [f64::INFINITY, 0.0, 0.0],
                    ..identity
                }],
                None,
            ),
            Err(super::CharacterViewError::NonFiniteData)
        ));
        assert!(matches!(
            super::CharacterRow::from_parts(
                super::RepresentationSpaceKind::FullStarPir,
                vec![Complex64::new(1.0, 0.0)],
                vec![identity],
                Some(0),
            ),
            Err(super::CharacterViewError::InvalidDimension)
        ));
        let spin_identity = super::SpinSeitzOperation {
            seitz: identity,
            su2: [1.0, 0.0, 0.0, 0.0],
        };
        assert!(matches!(
            super::SpinCharacterRow::from_parts(
                vec![Complex64::new(1.0, f64::INFINITY)],
                vec![spin_identity],
                1,
            ),
            Err(super::CharacterViewError::NonFiniteData)
        ));
        assert!(matches!(
            super::SpinCharacterRow::from_parts(
                vec![Complex64::new(1.0, 0.0)],
                vec![super::SpinSeitzOperation {
                    seitz: super::SeitzOperation {
                        translation: [f64::NAN, 0.0, 0.0],
                        ..identity
                    },
                    ..spin_identity
                }],
                1,
            ),
            Err(super::CharacterViewError::NonFiniteData)
        ));
        assert!(matches!(
            super::SpinCharacterRow::from_parts(
                vec![Complex64::new(1.0, 0.0)],
                vec![super::SpinSeitzOperation {
                    su2: [f64::NAN, 0.0, 0.0, 0.0],
                    ..spin_identity
                }],
                1,
            ),
            Err(super::CharacterViewError::NonFiniteData)
        ));

        let ordinary = crate::irrep::query::irreps_of(1)
            .iter()
            .find(|irrep| !irrep.spinor)
            .expect("ordinary scalar");
        let mut malformed = *ordinary;
        malformed._pir_rot_start += 1;
        assert!(matches!(
            malformed.ordinary_scalar_selected_arm_block_trace(),
            Err(super::CharacterViewError::MisalignedStorage)
        ));
        let compound = crate::irrep::query::irreps_of(23)
            .iter()
            .find(|irrep| irrep.cir_component_count() > 0)
            .expect("compound scalar");
        let mut malformed = *compound;
        malformed._cir_start += 1;
        assert!(matches!(
            malformed.compound_selected_arm_view(),
            Err(super::CharacterViewError::MisalignedStorage)
        ));
    }
}
