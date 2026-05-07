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

/// A single irreducible representation at a k-point, with labels in three conventions.
///
/// The three label systems are cross-referenced from Stokes & Hatch Table 7.
#[derive(Debug, Clone)]
pub struct IrrepData {
    /// Miller & Love label (e.g. `"GM1+"`, `"X3-"`, `"R4+"`)
    pub ml_label: &'static str,
    /// Kovalev label (e.g. `"τ1"`, `"k6τ3"`)
    pub kovalev_label: &'static str,
    /// Bradley & Cracknell / CDML label (e.g. `"Γ1+"`, `"Γ4-"`)
    pub bc_label: &'static str,
    /// Dimension of the irrep: 1, 2, 3, 4, or 6
    pub dimension: u8,
    /// Stokes-Hatch image symbol (e.g. `"A1a"`, `"C24c"`, `"B4a"`)
    pub image: &'static str,
    /// Basis functions (e.g. `"1"`, `"x,y,z"`, `"Sx,Sy,Sz"`)
    pub basis_functions: &'static str,
}

/// A high-symmetry k-point in the Brillouin zone with its little co-group and irreps.
#[derive(Debug, Clone)]
pub struct KPointData {
    /// k-point label: `"Γ"`, `"X"`, `"M"`, `"R"`, `"A"`, `"H"`, `"K"`, `"L"`, etc.
    pub label: &'static str,
    /// Fractional reciprocal coordinates `[kx, ky, kz]`
    pub coords: [f64; 3],
    /// Little co-group (point group of the wave-vector): `"m-3m"`, `"4/mmm"`, etc.
    pub little_group: &'static str,
    /// Irreducible representations at this k-point
    pub irreps: &'static [IrrepData],
}

/// An isotropy subgroup: the lower-symmetry space group obtained when the order
/// parameter condenses along a specific direction for a given irrep.
#[derive(Debug, Clone)]
pub struct IsotropySubgroup {
    /// Space group number (1–230)
    pub sg_number: u16,
    /// Hermann-Mauguin symbol (e.g. `"Pm-3m"`, `"R-3m"`)
    pub symbol: &'static str,
    /// Schoenflies symbol (e.g. `"Oh^1"`, `"D3d^5"`)
    pub schoenflies: &'static str,
    /// Order-parameter direction (e.g. `"(a,0,0)"`, `"(a,a,a)"`)
    pub direction: &'static str,
    /// Number of domains
    pub domains: u8,
    /// Basis vectors of the subgroup cell relative to the parent cell
    pub basis: &'static str,
    /// Origin shift relative to the parent cell
    pub origin: &'static str,
}

// ── Machine-generated record types (flat-array storage) ─────────────────────

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
    /// Number of matrix elements = opcount × dim² (fits in u16: max ~27648)
    pub(crate) _mat_count: u16,
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
    /// Start index into [`SPIN_EXTRA_CHARS`] (0 if no extra)
    pub(crate) _spin_extra_start: u32,
    /// Number of extra character values (0 if no extra)
    pub(crate) _spin_extra_count: u16,
    /// Start index into [`CIR_COMPONENT_CHARS`] (0 if not compound)
    pub(crate) _cir_start: u32,
    /// Number of CIR components (0 for non-compound irreps, 2 for Z1Z4 type)
    pub(crate) _cir_count: u8,
    /// Number of operations per CIR component
    pub(crate) _cir_ops: u8,
}

impl IrrepRecord {
    /// For spinor irreps: number of characters corresponding to the little-group
    /// operations (the first `n` values in [`Self::characters`]).
    /// Returns 0 for scalar irreps.
    pub fn spin_lg_char_count(&self) -> usize {
        self._spin_lg_count as usize
    }

    /// Little-group operation indices for spinor irreps.
    /// Maps local character position → global SPIN_OP index.
    pub fn spin_lg_op_indices(&self) -> &'static [u16] {
        if self._spin_lg_op_count == 0 {
            return &[];
        }
        let start = self._spin_lg_op_start as usize;
        let len = self._spin_lg_op_count as usize;
        &super::generated_data::SPIN_LG_OP_INDICES[start..start + len]
    }

    /// Extra character values for spinor Wigner test (Bilbao pre-computed).
    /// Sum of these values gives the Wigner indicator.
    pub fn spin_extra_chars(&self) -> &'static [f64] {
        if self._spin_extra_count == 0 {
            return &[];
        }
        let start = self._spin_extra_start as usize;
        let len = self._spin_extra_count as usize;
        &super::generated_data::SPIN_EXTRA_CHARS[start..start + len]
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
        if idx >= 230 { return (&[], &[]); }
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
        if char_count == 0 { return &[]; }
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

    /// Number of CIR (complex) components this PIR irrep decomposes into.
    /// 0 = non-compound, 2 = compound like Z1Z4 = Z1 ⊕ Z4.
    pub fn cir_component_count(&self) -> usize {
        self._cir_count as usize
    }

    /// Complex character table for a specific CIR component.
    ///
    /// Returns `(re, im)` pairs in CIR/ISOTROPY operation order.
    /// Use `cir_rotation_at()` for the corresponding operation rotations
    /// and `build_cir_index_map()` to map to H_ops order.
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
            [self._mat_start as usize..(self._mat_start + self._mat_count as u32) as usize]
    }

    /// Full irrep matrices reordered to match spglib H_ops order.
    ///
    /// Only those H_ops that match a PIR operation (via rotation matrix)
    /// get matrix data.  Unmatched ops get zero-filled blocks.
    /// Returns an empty `Vec` if matrices or rotation data are unavailable.
    pub fn matrices_reordered(&self, h_seitz: &[crate::irrep::wigner::SeitzOp]) -> Vec<f64> {
        let mats = self.matrices();
        let rots = self.pir_rotations();
        if mats.is_empty() || rots.is_empty() {
            return mats.to_vec();
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
                    return mats.to_vec();
                }
                crate::irrep::wigner::build_h_to_cir_map(&h_seitz[..h_count], rots)
                    .unwrap_or_else(|| {
                        // Last resort: identity mapping for first n_pir_ops
                        (0..n_pir_ops).collect()
                    })
            }
        };

        let mut reordered = vec![0.0f64; mats.len()];
        for h_idx in 0..h_to_pir.len().min(n_pir_ops) {
            let pir_idx = h_to_pir[h_idx];
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
        reordered
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
        let body = self.ml.trim_end_matches(|c: char| c == '+' || c == '-');
        let end = body.find(|c: char| c.is_ascii_digit()).unwrap_or(body.len());
        &body[..end]
    }

    /// Whether this is a special k-point (not a line or plane).
    pub fn is_point(&self) -> bool {
        let k = self.k_label();
        // Lines and planes have longer prefixes (DT, LD, SM, etc.)
        // Points have short prefixes (GM, X, M, R, A, H, K, L, etc.)
        k.len() <= 2 && !matches!(k, "GP")
    }

    /// Look up the pre-computed h_to_pir mapping for this irrep.
    ///
    /// Returns a slice `h_to_pir` where `h_to_pir[h_idx] = pir_idx` maps each
    /// spglib H operation position to the corresponding ISOTROPY (PIR) character index.
    fn reorder_map(&self) -> &'static [u16] {
        let sg_idx = self.sg as usize;
        if sg_idx == 0 || sg_idx > 230 { return &[]; }
        let (reo_start, reo_count) = super::generated_data::SG_CHAR_REORDER_INDEX[sg_idx];
        let reo_start = reo_start as usize;
        let reo_count = reo_count as usize;
        // Find this irrep's position within its SG's irrep list
        let (ir_start, ir_count) = super::generated_data::SG_IRREP_INDEX[sg_idx];
        let sg_irreps = &super::generated_data::IRREPS[ir_start as usize..(ir_start + ir_count) as usize];
        let pos = sg_irreps.iter().position(|r| std::ptr::eq(r, self));
        match pos {
            Some(p) if p < reo_count => {
                let (perm_start, perm_count) =
                    super::generated_data::IRREP_CHAR_REORDER_INDEX[reo_start + p];
                let start = perm_start as usize;
                let len = perm_count as usize;
                if len > 0 {
                    &super::generated_data::PER_CHAR_REORDER[start..start + len]
                } else {
                    &[]
                }
            }
            _ => &[],
        }
    }

    /// Character table reordered to spglib H_ops order.
    ///
    /// Unlike [`Self::characters`] which returns characters in ISOTROPY (PIR) order,
    /// this method returns one character value per spglib H operation, using the
    /// pre-computed h_to_pir mapping.
    ///
    /// The returned `Vec` has length equal to the number of spglib H operations.
    pub fn characters_spglib(&self) -> Vec<f64> {
        let chars = self.characters();
        let h_to_pir = self.reorder_map();
        if h_to_pir.is_empty() || chars.is_empty() {
            return chars.to_vec();
        }
        h_to_pir.iter()
            .map(|&pir_idx| {
                let idx = pir_idx as usize;
                if idx < chars.len() { chars[idx] } else { 0.0 }
            })
            .collect()
    }
}

impl IsotropyRecord {
    /// Human-readable description of this isotropy subgroup.
    pub fn describe(&self) -> String {
        format!("#{} {} ({}), domains={}", self.sg, self.symbol, self.schoenflies, self.domains)
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
