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
    /// Start index into [`CIR_COMPONENT_CHARS`] (0 if not compound)
    pub(crate) _cir_start: u32,
    /// Number of CIR components (0 for non-compound irreps, 2 for Z1Z4 type)
    pub(crate) _cir_count: u8,
    /// Number of operations per CIR component
    pub(crate) _cir_ops: u8,
}

impl IrrepRecord {
    /// Number of CIR (complex) components this PIR irrep decomposes into.
    /// 0 = non-compound, 2 = compound like Z1Z4 = Z1 ⊕ Z4.
    pub fn cir_component_count(&self) -> usize {
        self._cir_count as usize
    }

    /// Complex character table for a specific CIR component.
    ///
    /// Returns `(re, im)` pairs.  Use only when `cir_component_count() > 0`.
    pub fn cir_component_chars(&self, comp: usize) -> &'static [f64] {
        if comp >= self._cir_count as usize || self._cir_start == 0 {
            return &[];
        }
        let start = self._cir_start as usize + comp * self._cir_ops as usize * 2;
        let len = self._cir_ops as usize * 2;
        // SAFETY: generated data guarantees valid bounds
        &super::generated_data::CIR_COMPONENT_CHARS[start..start + len]
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
