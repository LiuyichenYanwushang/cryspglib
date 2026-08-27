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

// ── Compact record types (flat-array storage) ───────────────────────────────

/// The representation space in which a character table is expressed.
///
/// A physical (PIR) table can be induced from one or more arms of a star of
/// **k**.  Its [`IrrepRecord::dim`] and [`IrrepRecord::characters`] describe
/// that full-star representation.  CIR and spinor tables, on the other hand,
/// describe a single selected arm and its little group.  Keeping this fact in
/// the type of a read view prevents accidentally using the full-star dimension
/// with a selected-arm character row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrrepRepresentationSpace {
    /// The full physical representation, including all star arms.
    FullStar,
    /// One selected **k** arm and its little-group representation.
    SelectedArm,
}

/// Operation ordering used by an [`IrrepCharacterView`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrrepOperationOrder {
    /// PIR/Hall operation order, as used by [`IrrepRecord::characters`].
    Pir,
    /// CIR operation order for one selected-arm component.
    Cir,
    /// Local little-group order for a spinor table.  `get(i)` corresponds to
    /// `spin_lg_op_indices()[i]`.
    SpinLittleGroup,
}

/// A checked, read-only character table view.
///
/// The view exposes the dimension and operation order alongside the values.
/// Construction is fallible: malformed or incomplete generated data (in
/// particular an identity character different from the dimension) yields
/// `None` rather than a view that can be misused.  Character values are
/// returned as complex numbers so scalar, CIR, and spinor tables share one
/// access pattern.
#[derive(Debug, Clone, Copy)]
pub struct IrrepCharacterView {
    space: IrrepRepresentationSpace,
    order: IrrepOperationOrder,
    dimension: usize,
    values: CharacterViewValues,
}

#[derive(Debug, Clone, Copy)]
enum CharacterViewValues {
    Real(&'static [f64]),
    /// Interleaved `(re, im)` values, as used by CIR component storage.
    Complex(&'static [f64]),
    /// Separate real and imaginary arrays, as used by scalar selected-arm
    /// and spinor storage.
    Split(&'static [f64], &'static [f64]),
    /// Concatenated CIR component rows.  `components` rows, each of
    /// `operations` interleaved values, are summed lazily by `get`.
    Compound {
        values: &'static [f64],
        components: usize,
        operations: usize,
    },
}

impl IrrepCharacterView {
    fn checked_real(
        space: IrrepRepresentationSpace,
        order: IrrepOperationOrder,
        dimension: usize,
        values: &'static [f64],
    ) -> Option<Self> {
        if values.is_empty()
            || !values.iter().all(|value| value.is_finite())
            || (values[0] - dimension as f64).abs() > 1e-6
        {
            return None;
        }
        Some(Self {
            space,
            order,
            dimension,
            values: CharacterViewValues::Real(values),
        })
    }

    fn checked_complex(
        space: IrrepRepresentationSpace,
        order: IrrepOperationOrder,
        dimension: usize,
        real: &'static [f64],
        imag: &'static [f64],
    ) -> Option<Self> {
        if real.is_empty()
            || real.len() != imag.len()
            || !real.iter().chain(imag).all(|value| value.is_finite())
            || (real[0] - dimension as f64).abs() > 1e-6
            || imag[0].abs() > 1e-6
        {
            return None;
        }
        Some(Self {
            space,
            order,
            dimension,
            values: CharacterViewValues::Split(real, imag),
        })
    }

    fn checked_interleaved_complex(
        space: IrrepRepresentationSpace,
        order: IrrepOperationOrder,
        dimension: usize,
        values: &'static [f64],
    ) -> Option<Self> {
        if values.is_empty()
            || !values.len().is_multiple_of(2)
            || !values.iter().all(|value| value.is_finite())
            || (values[0] - dimension as f64).abs() > 1e-6
            || values[1].abs() > 1e-6
        {
            return None;
        }
        Some(Self {
            space,
            order,
            dimension,
            values: CharacterViewValues::Complex(values),
        })
    }

    fn checked_compound(
        dimension: usize,
        values: &'static [f64],
        components: usize,
        operations: usize,
    ) -> Option<Self> {
        if components == 0
            || operations == 0
            || values.len() != components.checked_mul(operations)?.checked_mul(2)?
            || !values.iter().all(|value| value.is_finite())
        {
            return None;
        }
        let identity_real = (0..components)
            .map(|component| values[component * operations * 2])
            .sum::<f64>();
        let identity_imag = (0..components)
            .map(|component| values[component * operations * 2 + 1])
            .sum::<f64>();
        if (identity_real - dimension as f64).abs() > 1e-6 || identity_imag.abs() > 1e-6 {
            return None;
        }
        Some(Self {
            space: IrrepRepresentationSpace::SelectedArm,
            order: IrrepOperationOrder::Cir,
            dimension,
            values: CharacterViewValues::Compound {
                values,
                components,
                operations,
            },
        })
    }

    /// Whether this view describes the full star or a selected arm.
    pub const fn representation_space(&self) -> IrrepRepresentationSpace {
        self.space
    }

    /// Operation order associated with [`Self::get`].
    pub const fn operation_order(&self) -> IrrepOperationOrder {
        self.order
    }

    /// Dimension in the same representation space as this view.
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Number of character columns.
    pub const fn len(&self) -> usize {
        match self.values {
            CharacterViewValues::Real(values) => values.len(),
            CharacterViewValues::Complex(values) => values.len() / 2,
            CharacterViewValues::Split(real, _) => real.len(),
            CharacterViewValues::Compound { operations, .. } => operations,
        }
    }

    /// Whether this view has no character columns.
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get a character value without risking an out-of-bounds panic.
    pub fn get(&self, index: usize) -> Option<num_complex::Complex64> {
        match self.values {
            CharacterViewValues::Real(values) => values
                .get(index)
                .copied()
                .map(|value| num_complex::Complex64::new(value, 0.0)),
            CharacterViewValues::Complex(values) => {
                let start = index.checked_mul(2)?;
                Some(num_complex::Complex64::new(
                    *values.get(start)?,
                    *values.get(start + 1)?,
                ))
            }
            CharacterViewValues::Split(real, imag) => Some(num_complex::Complex64::new(
                *real.get(index)?,
                *imag.get(index)?,
            )),
            CharacterViewValues::Compound {
                values,
                components,
                operations,
            } => {
                if index >= operations {
                    return None;
                }
                let mut value = num_complex::Complex64::new(0.0, 0.0);
                for component in 0..components {
                    let start = (component * operations + index).checked_mul(2)?;
                    value.re += *values.get(start)?;
                    value.im += *values.get(start + 1)?;
                }
                Some(value)
            }
        }
    }
}

/// How a scalar physical (PIR) compound record is assembled from CIR rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompoundCharacterSemantics {
    /// A CIR row and its complex conjugate are realified: χ = 2 Re(χ_CIR).
    ConjugateRealification,
    /// Distinct CIR constituents are directly summed: χ = Σ χ_CIR.
    DistinctComponentSum,
}

/// Error returned when stored irrep matrices cannot be aligned with a supplied
/// symmetry-operation order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum MatrixReorderError {
    /// No unambiguous rotation mapping exists between the two operation lists.
    #[error("could not map supplied symmetry operations to stored irrep rotations")]
    OperationMappingFailed,
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
    /// For spinor irreps: number of little-group ops (≤ _char_count).
    /// Extra character values beyond this count are antiunitary/auxiliary.
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
}

impl IrrepRecord {
    /// Rational wave vector associated with this irrep.
    pub const fn k_vector(&self) -> KVector {
        KVector::new([self.kx, self.ky, self.kz], self.kd)
    }

    /// For spinor irreps: number of characters corresponding to the little-group
    /// operations (the first `n` values in [`Self::characters`]).
    /// Returns 0 for scalar irreps.
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

    /// Complex characters of the first (stored-k) star-arm block.
    ///
    /// These are generated from the complex ISO-IR matrices and are needed
    /// when a physically irreducible real PIR combines conjugate k arms. The
    /// returned slices are empty when no aligned CIR matrix was available.
    pub fn scalar_little_characters(&self) -> (&'static [f64], &'static [f64]) {
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

    /// Return the checked character view of the full physical (PIR) star.
    ///
    /// For scalar records this is [`IrrepRepresentationSpace::FullStar`] in
    /// PIR order.  Spinor records do not have a full-star PIR character row
    /// in this database and therefore return `None`; use
    /// [`Self::selected_arm_character_view`] for those records.
    pub fn full_star_character_view(&self) -> Option<IrrepCharacterView> {
        (!self.spinor).then(|| {
            IrrepCharacterView::checked_real(
                IrrepRepresentationSpace::FullStar,
                IrrepOperationOrder::Pir,
                self.dim as usize,
                self.characters(),
            )
        })?
    }

    /// Return a checked character view for one selected arm.
    ///
    /// For an ordinary scalar record `component` must be zero and the view is
    /// backed by the aligned `scalar_little_characters` arrays.  A compound
    /// scalar record exposes one CIR component per index, with interleaved
    /// `(real, imag)` values.  A spinor record has one component (index zero)
    /// and uses the first `spin_lg_char_count()` entries plus the stored
    /// imaginary parts.  Every view validates its identity character against
    /// its selected-arm dimension before being returned.
    pub fn selected_arm_character_view(&self, component: usize) -> Option<IrrepCharacterView> {
        if self.spinor {
            if component != 0 {
                return None;
            }
            let count = self.spin_lg_char_count();
            let chars = self.characters().get(..count)?;
            let imag = self.spin_character_imag().get(..count)?;
            let dimension = chars.first().copied()?;
            if !dimension.is_finite() || dimension.fract() != 0.0 || dimension < 0.0 {
                return None;
            }
            return IrrepCharacterView::checked_complex(
                IrrepRepresentationSpace::SelectedArm,
                IrrepOperationOrder::SpinLittleGroup,
                dimension as usize,
                chars,
                imag,
            );
        }

        if self.cir_component_count() > 0 {
            if component >= self.cir_component_count() {
                return None;
            }
            let operations = self._cir_ops as usize;
            let start =
                self._cir_start as usize + component.checked_mul(operations)?.checked_mul(2)?;
            let end = start.checked_add(operations.checked_mul(2)?)?;
            let values = super::generated_data::CIR_COMPONENT_CHARS.get(start..end)?;
            let dimension = values.first().copied()?;
            if !dimension.is_finite() || dimension.fract() != 0.0 || dimension < 0.0 {
                return None;
            }
            return IrrepCharacterView::checked_interleaved_complex(
                IrrepRepresentationSpace::SelectedArm,
                IrrepOperationOrder::Cir,
                dimension as usize,
                values,
            );
        }

        if component != 0 {
            return None;
        }
        let (real, imag) = self.scalar_little_characters();
        let dimension = real.first().copied()?;
        if !dimension.is_finite() || dimension.fract() != 0.0 || dimension < 0.0 {
            return None;
        }
        IrrepCharacterView::checked_complex(
            IrrepRepresentationSpace::SelectedArm,
            IrrepOperationOrder::Pir,
            dimension as usize,
            real,
            imag,
        )
    }

    /// Return the selected-arm character view of the complete scalar little
    /// representation.  For a compound record this is the direct sum of all
    /// CIR component rows, so its dimension is the sum of component
    /// dimensions (for example SG76 `R1R2` has full-star dimension 4 and
    /// selected-arm dimension 2).  For an ordinary scalar or spinor record it
    /// is equivalent to component zero from [`Self::selected_arm_character_view`].
    pub fn selected_arm_little_group_view(&self) -> Option<IrrepCharacterView> {
        if self.cir_component_count() == 0 {
            return self.selected_arm_character_view(0);
        }
        let operations = self._cir_ops as usize;
        let components = self._cir_count as usize;
        let total = components.checked_mul(operations)?.checked_mul(2)?;
        let start = self._cir_start as usize;
        let end = start.checked_add(total)?;
        let values = super::generated_data::CIR_COMPONENT_CHARS.get(start..end)?;
        let dimension = (0..components)
            .map(|component| values.get(component * operations * 2).copied())
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .sum::<f64>();
        if !dimension.is_finite() || dimension.fract() != 0.0 || dimension < 0.0 {
            return None;
        }
        IrrepCharacterView::checked_compound(
            dimension as usize,
            // CIR components are contiguous in the generated flat array.
            values,
            components,
            operations,
        )
    }

    /// Return either the full-star or selected-arm view requested by the
    /// caller.  `component` is ignored for the full-star case and must be
    /// zero for non-compound selected-arm data.
    pub fn character_view(
        &self,
        space: IrrepRepresentationSpace,
        component: usize,
    ) -> Option<IrrepCharacterView> {
        match space {
            IrrepRepresentationSpace::FullStar if component == 0 => self.full_star_character_view(),
            IrrepRepresentationSpace::SelectedArm => self.selected_arm_character_view(component),
            _ => None,
        }
    }

    /// Number of CIR (complex) components this PIR irrep decomposes into.
    /// 0 = non-compound, 2 = compound like Z1Z4 = Z1 ⊕ Z4.
    pub fn cir_component_count(&self) -> usize {
        self._cir_count as usize
    }

    /// Complex character table for a specific CIR component.
    ///
    /// Returns `(re, im)` pairs in the generated data-Hall operation order.
    /// [`Self::cir_rotations`] contains the corresponding rotations.
    pub fn cir_component_chars(&self, comp: usize) -> &'static [f64] {
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
    pub fn cir_rotations(&self, comp: usize) -> &'static [i32] {
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

    /// Semantics used to assemble this physical compound character table from
    /// its authoritative CIR constituents.
    ///
    /// Returns `None` for non-compound records and for a malformed compound
    /// label for which constituent identity cannot be recovered.  In
    /// particular, this accessor never guesses from character values or Gram
    /// products.  The generated CIR component count is used only to establish
    /// that this record is a compound record; constituent identity comes from
    /// the parsed constituent labels encoded in the ML record.
    pub fn compound_character_semantics(&self) -> Option<CompoundCharacterSemantics> {
        if self._cir_count == 0 {
            return None;
        }
        let (first, second) = split_compound_constituents(self.ml)?;
        Some(if first == second {
            CompoundCharacterSemantics::ConjugateRealification
        } else {
            CompoundCharacterSemantics::DistinctComponentSum
        })
    }
}

/// Parse the two constituent labels in the compact ML spelling used for a
/// compound record.  This accepts e.g. `P2P2`, `H1H2`, and `D1+D2+`, while
/// rejecting a single ordinary label.  It deliberately does not inspect any
/// numerical character values.
fn split_compound_constituents(label: &str) -> Option<(&str, &str)> {
    fn constituent_end(value: &str) -> Option<usize> {
        let bytes = value.as_bytes();
        let mut index = 0;
        while index < bytes.len() && bytes[index].is_ascii_alphabetic() {
            index += 1;
        }
        let digit_start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == 0 || index == digit_start {
            return None;
        }
        if bytes
            .get(index)
            .is_some_and(|byte| *byte == b'+' || *byte == b'-')
        {
            index += 1;
        }
        Some(index)
    }

    for split in 1..label.len() {
        let (first, second) = label.split_at(split);
        if constituent_end(first) == Some(first.len())
            && constituent_end(second) == Some(second.len())
        {
            return Some((first, second));
        }
    }
    None
}

impl IrrepRecord {
    /// Character table: χ(g) = Tr(D(g)) for each space-group operator.
    ///
    /// The character χ(g) of a representation D is the trace of the
    /// representation matrix for each symmetry operation g.  The return
    /// slice has length equal to the number of operators in the little
    /// co-group of the wave-vector, and each entry is a floating-point
    /// value (possibly negative, fractional, or zero).
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
mod character_view_tests {
    use super::{CompoundCharacterSemantics, IrrepRepresentationSpace};

    #[test]
    fn compound_semantics_are_available_for_every_generated_compound() {
        let mut compounds = 0;
        let mut realified = 0;
        let mut distinct = 0;
        for sg in 1..=230 {
            for irrep in crate::irrep::query::irreps_of(sg) {
                if irrep.cir_component_count() == 0 {
                    assert!(irrep.compound_character_semantics().is_none());
                    continue;
                }
                compounds += 1;
                match irrep.compound_character_semantics() {
                    Some(CompoundCharacterSemantics::ConjugateRealification) => realified += 1,
                    Some(CompoundCharacterSemantics::DistinctComponentSum) => distinct += 1,
                    None => panic!("missing semantics for SG{} {}", sg, irrep.ml),
                }
            }
        }
        assert_eq!(compounds, 672);
        assert_eq!(realified + distinct, compounds);
        assert!(realified > 0 && distinct > 0);
    }

    #[test]
    fn selected_arm_view_does_not_use_full_star_dimension() {
        let irrep = crate::irrep::query::irreps_of(76)
            .iter()
            .find(|irrep| !irrep.spinor && irrep.ml == "R1R2")
            .expect("SG76 R1R2 compound irrep");
        let full = irrep.full_star_character_view().expect("full PIR view");
        let selected = irrep
            .selected_arm_little_group_view()
            .expect("selected CIR view");
        assert_eq!(
            full.representation_space(),
            IrrepRepresentationSpace::FullStar
        );
        assert_eq!(full.dimension(), 4);
        assert_eq!(
            selected.representation_space(),
            IrrepRepresentationSpace::SelectedArm
        );
        assert_eq!(selected.dimension(), 2);
        assert_eq!(full.len(), selected.len());
        assert_eq!(full.get(usize::MAX), None);
    }

    #[test]
    fn representative_compound_records_use_their_constituent_identity() {
        use CompoundCharacterSemantics::{
            ConjugateRealification as Real, DistinctComponentSum as Sum,
        };
        for (sg, label, expected) in [
            (199, "P2P2", Real),
            (220, "P1P1", Real),
            (220, "P2P2", Real),
            (220, "P3P3", Real),
            (182, "H1H2", Sum),
            (46, "W1W2", Sum),
        ] {
            let irrep = crate::irrep::query::irreps_of(sg)
                .iter()
                .find(|irrep| irrep.ml == label)
                .expect("compound regression record");
            assert_eq!(irrep.compound_character_semantics(), Some(expected));
        }
    }
}
