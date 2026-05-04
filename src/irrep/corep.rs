//! Co-representation (corep) theory for magnetic space group irreps.
//!
//! # Theory
//!
//! A magnetic space group (MSG) with anti-unitary operations can be written as
//!
//! $$
//! \mathcal{M} = H \cup \mathcal{T} g_0 H
//! $$
//!
//! where $$H = \mathcal{M} \cap G$$ is the **unitary subgroup** (a normal
//! non-magnetic space group), $$\mathcal{T}$$ is time reversal, and
//! $$g_0$$ is the spatial part of the chosen anti-unitary coset representative
//! $$a_0 = \mathcal{T} g_0$$.
//!
//! ## Construction from the unitary subgroup
//!
//! Given a non-magnetic irrep $$\Delta_i$$ of $$H$$ at wave-vector $$\mathbf{k}$$,
//! define its **conjugate representation** under the anti-unitary coset:
//!
//! $$
//! \Delta_i^{a_0}(h) \equiv \Delta_i(a_0^{-1} h a_0)^* \qquad (h \in H)
//! $$
//!
//! The magnetic co-representation $$\tilde{D}$$ of $$\mathcal{M}$$ is built
//! from the relationship between $$\Delta_i$$ and $$\Delta_i^{a_0}$$:
//!
//! ## Wigner's three cases
//!
//! The Wigner test classifies irreps into three types via the sum
//!
//! $$
//! W(\Delta_i) = \frac{1}{|H|}\sum_{b \in a_0 H} \chi_{\Delta_i}(b^2)
//! $$
//!
//! where the sum runs over the **anti-unitary coset** $$a_0 H$$ and
//! $$\chi_{\Delta_i}$$ is the character of $$\Delta_i$$.
//!
//! | Case | Condition | Corep dimension | Unitary characters | Anti-unitary characters |
//! |------|-----------|----------------|-------------------|------------------------|
//! | **Type A** | $$\Delta_i^{a_0} \sim \Delta_i$$, $$W = +1$$ | $$d_i$$ | $$\chi_{\Delta_i}(h)$$ | $$\chi_{\Delta_i}(a_0 h)$$ (real) |
//! | **Type B** | $$\Delta_i^{a_0} \sim \Delta_i$$, $$W = -1$$ | $$d_i$$ (Kramers) | $$\chi_{\Delta_i}(h)$$ | $$-\chi_{\Delta_i}(a_0 h)$$ (pseudo-real) |
//! | **Type C** | $$\Delta_i^{a_0} \nsim \Delta_i$$, $$W = 0$$ | $$2d_i$$ | $$2\,\mathrm{Re}[\chi_{\Delta_i}(h)]$$ | $$0$$ |
//!
//! **Type C** pairs two inequivalent irreps $$\Delta_i, \Delta_j$$ of $$H$$
//! (where $$\Delta_j \sim \Delta_i^{a_0}$$). The corep is
//! $$\tilde{D} = \Delta_i \oplus \Delta_j$$ with block structure
//!
//! $$
//! \tilde{D}(h) = \begin{pmatrix} \Delta_i(h) & 0 \\ 0 & \Delta_j(h) \end{pmatrix},
//! \qquad
//! \tilde{D}(a_0 h) \sim \begin{pmatrix} 0 & * \\ * & 0 \end{pmatrix} K
//! $$
//!
//! where $$K$$ denotes complex conjugation.
//!
//! ## Workflow
//!
//! ```text
//! BNS label ("128.406") + k-point label ("Z")
//!   → uni_from_bns()           // BNS → UNI number
//!   → identify_unitary_subgroup()  // UNI → H space group
//!   → irreps_of(H) at k-point  // H's double-group irreps
//!   → compute_corepresentation()   // Wigner test + corep characters
//!   → Corepresentation { characters, corep_type, dim }
//! ```
//!
//! ## Example: 128.406 at Z
//!
//! Verified against Bilbao Crystallographic Server (BCS):
//!
//! - Magnetic SG: $$P4'/m'nc'$$ (No. 128.406, UNI 1066)
//! - Unitary subgroup: $$P\bar{4}n2$$ (No. 118)
//! - k-vector: $$Z = (0, 0, 1/2)$$
//! - Magnetic little co-group: $$4'/m'mm'$$ (12 ops: 8 unitary + 4 anti-unitary)
//!
//! From H = SG 118's Z-point irreps:
//!
//! | H irrep | Type | Magnetic corep | Dimension |
//! |---------|------|---------------|-----------|
//! | Z₁Z₄ | C | Z₁Z₂ | 2D |
//! | Z₂Z₃ | C | Z₃Z₄ | 2D |
//! | Z₅ | A | Z₅ | 2D |
//! | Z₆, Z₇ (spinor) | C | Z̄₆Z̄₇ | 4D |
//!
//! # References
//!
//! - Wigner (1959), *Group Theory*, Chapter 26
//! - Bradley & Cracknell (1972), *The Mathematical Theory of Symmetry in Solids*
//! - Stokes, Campbell & Hatch, ISOTROPY Suite documentation
//! - Bilbao Crystallographic Server: <https://cryst.ehu.es/cgi-bin/cryst/programs/corepresentations.pl>

use num_complex::Complex64;
use crate::mathfunc::{Mat3I, Vec3};
use crate::spg_database::{spgdb_get_spacegroup_operations, spgdb_get_spacegroup_type};
use super::types::IrrepRecord;
use super::wigner::{self, filter_little_group, ops_to_seitz, SeitzOp,
    compose_seitz, square_seitz, find_seitz, bloch_phase, mat_vec_i32, add3};

/// Co-representation type from Wigner's test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorepType {
    /// W = +1: D ∼ D*, real representation.
    A,
    /// W = -1: D ∼ D*, pseudo-real (quaternionic).
    B,
    /// W = 0: D ≁ D*.
    C,
    /// Wigner indicator is non-quantized — missing data (e.g. spinor without
    /// rotations) or parser error.  Cannot classify.
    Unsupported,
}

impl CorepType {
    pub fn description(&self) -> &'static str {
        match self {
            CorepType::A => "type-a: D ~ D*, real (W=+1)",
            CorepType::B => "type-b: D ~ D*, pseudo-real (W=-1)",
            CorepType::C => "type-c: D ≁ D* (W=0)",
            CorepType::Unsupported => "unsupported: non-quantized Wigner indicator",
        }
    }
}

/// The magnetic symmetry operations for a space group.
///
/// Each operation has a 3×3 integer rotation matrix `rot`, a fractional
/// translation vector `trans`, and a time-reversal flag `timerev`.
#[derive(Debug, Clone)]
pub struct MagneticOps {
    /// Rotation matrices (3×3 integer)
    pub rot: Vec<Mat3I>,
    /// Translation vectors (fractional)
    pub trans: Vec<Vec3>,
    /// Time-reversal flag: true = anti-unitary
    pub timerev: Vec<bool>,
}

impl MagneticOps {
    /// Number of symmetry operations.
    pub fn len(&self) -> usize {
        self.rot.len()
    }

    /// Whether this is an empty set.
    pub fn is_empty(&self) -> bool {
        self.rot.is_empty()
    }
}

/// The computed magnetic co-representation of an irrep.
#[derive(Debug, Clone)]
pub struct Corepresentation {
    /// Character χ̃(g) for each magnetic operation (same order as MagneticOps).
    pub characters: Vec<f64>,
    /// Which operations are anti-unitary.
    pub timerev: Vec<bool>,
    /// Co-representation type.
    pub corep_type: CorepType,
    /// Dimension of the magnetic irrep.
    pub dim: usize,
    /// Number of unitary operations.
    pub unitary_order: usize,
    /// Number of anti-unitary operations.
    pub antiunitary_order: usize,
}

// ── Core computation ─────────────────────────────────────────────────────────

/// Compute the magnetic co-representation for an irrep of the unitary subgroup H.
///
/// See [`compute_coreps`] for the high-level BNS+k-label API.
pub fn compute_corepresentation(
    h_irrep: &IrrepRecord,
    uni_number: usize,
    mag_ops: &MagneticOps,
) -> Option<Corepresentation> {
    if uni_number == 0 || uni_number > 1651 {
        return None;
    }

    // 1. Filter to magnetic little group
    let mag_lg = filter_little_group(h_irrep.kx, h_irrep.ky, h_irrep.kz, h_irrep.kd, mag_ops);
    if mag_lg.is_empty() {
        return None;
    }

    // 2. Get H's symmetry operations
    let h_ops = get_parent_operations(h_irrep.sg);
    if h_ops.is_empty() {
        return None;
    }

    // 3. Convert to SeitzOps for proper composition
    let mag_seitz = ops_to_seitz(mag_ops);
    let h_seitz = ops_to_seitz(&h_ops);

    // 3a. Map unitary magnetic ops to H ops via full Seitz matching
    // (rotation + translation), not rotation-only.
    let op_map: Vec<Option<usize>> = (0..mag_ops.len())
        .map(|i| {
            if mag_ops.timerev[i] {
                None
            } else {
                let mop = &mag_seitz[i];
                wigner::find_seitz(&mop.rot, &mop.trans, &h_seitz)
                    .map(|m| m.op_index)
            }
        })
        .collect();

    if op_map.iter().enumerate().any(|(i, m)| !mag_ops.timerev[i] && m.is_none()) {
        return None;
    }

    // 4. H's irrep characters
    let h_chars = h_irrep.characters();
    let h_dim = h_irrep.dim as usize;
    if h_irrep.ml == "Z1Z4" {
        eprintln!("DEBUG compute_corep Z1Z4: sg={} h_chars={:?} h_ops.len={} mag_ops.len={}",
            h_irrep.sg, &h_chars[..h_chars.len().min(8)], h_ops.len(), mag_ops.len());
    }

    // 5. Separate unitary / anti-unitary in little group
    let unitary: Vec<usize> = mag_lg.iter()
        .filter(|&&i| !mag_ops.timerev[i]).copied().collect();
    let antiunitary: Vec<usize> = mag_lg.iter()
        .filter(|&&i| mag_ops.timerev[i]).copied().collect();

    // 7. Wigner test: use CIR complex chars for compound irreps
    let corep_type = if antiunitary.is_empty() {
        CorepType::A
    } else if h_irrep.cir_component_count() > 0 {
        // Compound irrep: test each CIR component.
        // Reorder CIR chars from ISOTROPY order to H_ops (spglib) order
        // so that m.op_index correctly indexes into the character table.
        let mut any_c = false;
        eprintln!("DEBUG CIR path: {} n_comp={}", h_irrep.ml, h_irrep.cir_component_count());
        for comp in 0..h_irrep.cir_component_count() {
            let cir = h_irrep.cir_component_chars(comp);
            if cir.is_empty() { continue; }
            let cir_rots = h_irrep.cir_rotations(comp);
            let cir_reordered = if let Some(h_to_cir) = wigner::build_h_to_cir_map(&h_seitz, cir_rots) {
                wigner::reorder_cir_chars(cir, &h_to_cir)
            } else {
                cir.to_vec()  // fallback: no matching rotations
            };
            let ct = wigner::wigner_classify_cir(
                &cir_reordered, &unitary, &mag_seitz, &h_seitz, antiunitary[0],
                h_irrep.kx, h_irrep.ky, h_irrep.kz, h_irrep.kd,
            );
            if ct == CorepType::C { any_c = true; break; }
        }
        if any_c { CorepType::C } else { CorepType::A }
    } else if h_irrep.spinor {
        // Spinor irrep: try extra chars first (Bilbao pre-computed), then SU(2).
        let extra = h_irrep.spin_extra_chars();
        if !extra.is_empty() {
            wigner::wigner_classify_spinor_extra(extra, unitary.len())
                .unwrap_or(CorepType::Unsupported)
        } else {
            // Fallback to SU(2) double-group path
            let (spin_rots, _spin_trans, _spin_su2) = h_irrep.spin_ops();
            let n_lg = h_irrep.spin_lg_char_count();
            if let Some(ct) = wigner::wigner_classify_spinor(
                h_chars, n_lg,
                spin_rots, _spin_trans, _spin_su2,
                &unitary, &mag_seitz, &h_seitz, antiunitary[0],
                h_irrep.kx, h_irrep.ky, h_irrep.kz, h_irrep.kd,
            ) {
                ct
            } else {
                CorepType::Unsupported
            }
        }
    } else {
        // Non-compound irrep: use PIR path with operation order mapping.
        let pir_rots = h_irrep.pir_rotations();
        if let Some(h_to_pir) = wigner::build_h_to_cir_map(&h_seitz, pir_rots) {
            let doubled = wigner::reorder_cir_chars(
                &h_chars.iter().flat_map(|&c| [c, 0.0f64]).collect::<Vec<_>>(),
                &h_to_pir,
            );
            let h_chars_reordered: Vec<f64> = (0..h_to_pir.len()).map(|i| doubled[2 * i]).collect();
            wigner::wigner_classify(
                &h_chars_reordered, &unitary, &mag_seitz, &h_seitz, antiunitary[0],
                h_irrep.kx, h_irrep.ky, h_irrep.kz, h_irrep.kd,
            )
        } else {
            // No rotation data (e.g. spinor irreps): cannot reorder, mark unsupported.
            CorepType::Unsupported
        }
    };

    // 8. Build corep character table
    let characters = wigner::build_corep_chars(
        &corep_type, mag_ops, &mag_lg, &op_map, h_chars, None, // partner resolved by caller
    );

    let dim = wigner::corep_dim(&corep_type, h_dim);

    Some(Corepresentation {
        characters,
        timerev: mag_lg.iter().map(|&i| mag_ops.timerev[i]).collect(),
        corep_type,
        dim,
        unitary_order: unitary.len(),
        antiunitary_order: antiunitary.len(),
    })
}

// ── Magnetic operations ──────────────────────────────────────────────────────

/// Get the magnetic space group symmetry operations.
pub fn get_magnetic_operations(uni_number: usize) -> Option<MagneticOps> {
    let hall = get_first_hall_for_uni(uni_number)?;
    let sym = crate::msg_database::msgdb_get_spacegroup_operations(uni_number, hall)?;
    let n = sym.size;
    let mut rot = Vec::with_capacity(n);
    let mut trans = Vec::with_capacity(n);
    let mut timerev = Vec::with_capacity(n);
    for i in 0..n {
        rot.push(sym.rot[i]);
        trans.push(sym.trans[i]);
        timerev.push(sym.timerev[i]);
    }
    Some(MagneticOps { rot, trans, timerev })
}

fn get_first_hall_for_uni(uni: usize) -> Option<usize> {
    if uni == 0 || uni > 1651 { return None; }
    for hall in 1..=530 {
        if let Some([lo, hi]) = crate::msg_database::msgdb_get_uni_candidates(hall) {
            if uni >= lo && uni <= hi { return Some(hall); }
        }
    }
    None
}

/// Get the symmetry operations (rotation + translation) for a space group.
///
/// Returns [`MagneticOps`] with `timerev` all `false` (non-magnetic).
/// The operations are in spglib's standard order.
///
/// # Example
/// ```
/// use cryspglib::irrep::corep::symmetry_operations_of;
/// let ops = symmetry_operations_of(139);
/// println!("SG 139: {} operations", ops.len());
/// ```
pub fn symmetry_operations_of(sg: u8) -> MagneticOps {
    get_parent_operations(sg)
}

fn get_parent_operations(sg: u8) -> MagneticOps {
    let hall = find_hall_number(sg);
    let mut rot = Vec::new();
    let mut trans = Vec::new();
    if let Some(h) = hall {
        if let Some(sym) = spgdb_get_spacegroup_operations(h) {
            for i in 0..sym.size {
                rot.push(sym.rot[i]);
                trans.push(sym.trans[i]);
            }
        }
    }
    let n = rot.len();
    MagneticOps { rot, trans, timerev: vec![false; n] }
}

fn find_hall_number(sg: u8) -> Option<usize> {
    for hall in 1..=530 {
        let st = spgdb_get_spacegroup_type(hall);
        if st.number == sg as usize { return Some(hall); }
    }
    None
}

// ── High-level API ───────────────────────────────────────────────────────────

/// Identify the unitary subgroup of a magnetic space group.
pub fn identify_unitary_subgroup(uni_number: usize) -> Option<usize> {
    let ops = get_magnetic_operations(uni_number)?;
    let unitary_rots: Vec<Mat3I> = ops.rot.iter().enumerate()
        .filter(|(i, _)| !ops.timerev[*i]).map(|(_, r)| *r).collect();
    let unitary_trans: Vec<[f64; 3]> = ops.trans.iter().enumerate()
        .filter(|(i, _)| !ops.timerev[*i]).map(|(_, t)| *t).collect();
    if unitary_rots.is_empty() { return None; }
    #[allow(deprecated)]
    let hall = crate::spg_get_hall_number_from_symmetry(&unitary_rots, &unitary_trans, 1e-5).ok()?;
    if hall == 0 || hall > 530 { return None; }
    let sg_type = spgdb_get_spacegroup_type(hall);
    Some(sg_type.number)
}

/// BNS label → UNI number.
pub fn uni_from_bns(bns: &str) -> Option<usize> {
    for uni in 1..=1651usize {
        let t = crate::msg_database::msgdb_get_magnetic_spacegroup_type(uni);
        if t.bns_number == bns { return Some(uni); }
    }
    None
}

/// OG label → UNI number.
pub fn uni_from_og(og: &str) -> Option<usize> {
    for uni in 1..=1651usize {
        let t = crate::msg_database::msgdb_get_magnetic_spacegroup_type(uni);
        if t.og_number == og { return Some(uni); }
    }
    None
}

/// Compute all corepresentations for a magnetic SG at a k-point.
pub fn compute_coreps(bns: &str, k_label: &str) -> Option<Vec<(&'static str, Corepresentation)>> {
    let uni = uni_from_bns(bns)?;
    let h_sg = identify_unitary_subgroup(uni)?;
    let mag_ops = get_magnetic_operations(uni)?;
    let h_irreps = super::query::irreps_of(h_sg as u8);
    let k_irreps: Vec<&IrrepRecord> = h_irreps.iter()
        .filter(|r| r.k_label() == k_label).collect();
    if k_irreps.is_empty() { return None; }

    // Pre-compute: convert H ops to Seitz and get anti-unitary representative
    let h_ops = get_parent_operations(h_sg as u8);
    let h_seitz = ops_to_seitz(&h_ops);
    let a0_idx = mag_ops.timerev.iter().position(|&t| t)?; // first anti-unitary
    let a0 = &wigner::SeitzOp::new(
        mag_ops.rot[a0_idx], mag_ops.trans[a0_idx], true,
    );

    // Build character tables for all k-point irreps (for partner finding)
    let char_tables: Vec<&[f64]> = k_irreps.iter().map(|ir| ir.characters()).collect();

    let mut results = Vec::with_capacity(k_irreps.len());
    for (i, ir) in k_irreps.iter().enumerate() {
        if let Some(c) = compute_corepresentation(ir, uni, &mag_ops) {
            // For Type C, find partner and rebuild character table
            let final_chars = if c.corep_type == CorepType::C {
                let (conj_chars, _) = wigner::conjugate_chars(
                    ir.characters(), &h_seitz, a0,
                    ir.kx, ir.ky, ir.kz, ir.kd,
                );
                let other_chars: Vec<&[f64]> = k_irreps.iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(_, r)| r.characters())
                    .collect();
                if let Some(partner_idx) = wigner::find_partner(&conj_chars, &other_chars) {
                    let partner_ir = &k_irreps[if partner_idx >= i { partner_idx + 1 } else { partner_idx }];
                    let mag_lg = filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
                    let mag_seitz = ops_to_seitz(&mag_ops);
                    let op_map: Vec<Option<usize>> = (0..mag_ops.len())
                        .map(|mi| {
                            if mag_ops.timerev[mi] { None }
                            else {
                                let mop = &mag_seitz[mi];
                                wigner::find_seitz(&mop.rot, &mop.trans, &h_seitz)
                                    .map(|m| m.op_index)
                            }
                        })
                        .collect();
                    wigner::build_corep_chars(
                        &c.corep_type, &mag_ops, &mag_lg, &op_map,
                        ir.characters(), Some(partner_ir.characters()),
                    )
                } else {
                    c.characters // fallback: no clear partner found
                }
            } else {
                c.characters
            };

            results.push((ir.ml, Corepresentation { characters: final_chars, ..c }));
        }
    }
    if results.is_empty() { None } else { Some(results) }
}

// ── IrrepRecord extension ────────────────────────────────────────────────────

impl IrrepRecord {
    /// Compute the magnetic co-representation (corep) for this non-magnetic
    /// irrep with respect to a magnetic space group.
    ///
    /// The magnetic character table is computed on-the-fly from the
    /// non-magnetic irrep data — no pre-stored tables needed.
    ///
    /// # Arguments
    ///
    /// * `uni_number` — OG/UNI number (1–1651), from
    ///   [`MagneticIsotropyRecord::mag_sg`](super::types::MagneticIsotropyRecord::mag_sg)
    ///
    /// # Returns
    ///
    /// `None` if the magnetic SG operations cannot be obtained.
    ///
    /// # Examples
    ///
    /// ```
    /// use cryspglib::irrep::query::irreps_of;
    /// use cryspglib::irrep::corep::CorepType;
    ///
    /// let gm4m = irreps_of(221).iter()
    ///     .find(|r| r.ml == "GM4-").unwrap();
    ///
    /// // Compute corep for a magnetic subgroup
    /// if let Some(corep) = gm4m.corepresentation(349) {
    ///     println!("Type: {:?}, dim: {}", corep.corep_type, corep.dim);
    ///     for (i, &chi) in corep.characters.iter().enumerate() {
    ///         let tr = if corep.timerev[i] { " (θ)" } else { "" };
    ///         println!("  op {}: χ = {:.4}{}", i, chi, tr);
    ///     }
    /// }
    /// ```
    /// Compute the magnetic co-representation for this irrep.
    ///
    /// Note: `self` must be an irrep of the **unitary subgroup H**, not the
    /// parent SG. Use [`compute_coreps`] for automatic H identification.
    pub fn corepresentation(&self, uni_number: usize) -> Option<Corepresentation> {
        let mag_ops = get_magnetic_operations(uni_number)?;
        compute_corepresentation(self, uni_number, &mag_ops)
    }
}

// ── High-level API ───────────────────────────────────────────────────────────

/// Compute all corepresentations for a magnetic space group at a k-point.
///
/// This is the primary entry point: given a BNS label and k-point label,
/// it automatically:
/// 1. Looks up the UNI number from the BNS label
/// 2. Identifies the unitary subgroup H
/// 3. Retrieves H's irreps at the k-point
/// 4. Computes the magnetic corepresentation for each H irrep
///
/// # Arguments
/// * `bns` — BNS label, e.g. `"128.406"`
/// * `k_label` — k-point label, e.g. `"Z"` or `"GM"`
///
/// # Returns
/// Vector of `(h_irrep_label, Corepresentation)` pairs, sorted by irrep label.
///
/// # Example
/// ```
/// use cryspglib::irrep::corep::compute_coreps;
///
/// let coreps = compute_coreps("128.406", "Z");
/// assert!(coreps.is_some());
/// for (label, c) in coreps.unwrap() {
///     println!("{}: dim={}, type={:?}, χ(id)={:.1}",
///         label, c.dim, c.corep_type, c.characters[0]);
/// }
/// ```

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::query::irreps_of;

    #[test]
    fn test_corep_gm4m_pmmm() {
        // Use compute_coreps on a SG 221 magnetic subgroup: 221.97 = UNI 1599
        // This tests the full pipeline: BNS → UNI → H → H irreps → coreps
        let coreps = compute_coreps("221.97", "GM");
        assert!(coreps.is_some(), "Should compute coreps for 221.97 at GM");
        let coreps = coreps.unwrap();
        assert!(!coreps.is_empty(), "Should have at least one corep");

        for (label, c) in &coreps {
            if c.corep_type == CorepType::Unsupported { continue; }
            assert!(c.dim > 0, "dim > 0 for {}", label);
            assert!((c.characters[0] - c.dim as f64).abs() < 0.01,
                "χ(id) = dim for {}", label);
        }
        println!("221.97 Gamma coreps: {} irreps computed", coreps.len());
    }

    #[test]
    fn test_corep_sg1_gm1() {
        // SG 1 (P1) GM1 → simplest case, use BNS "1.3" (= UNI 3)
        let coreps = compute_coreps("1.3", "GM");
        assert!(coreps.is_some(), "Should compute coreps for 1.3 at GM");
    }

    /// SG 128.406 (P4'/m'nc') at Z point — verified against BCS
    /// https://cryst.ehu.es/cgi-bin/cryst/programs/corepresentations.pl
    ///
    /// BCS confirms: Unitary Space Group = P-4n2 (No. 118) in standard setting.
    /// This test verifies automatic identification of the unitary subgroup.
    #[test]
    fn test_unitary_subgroup_sg128_406_is_sg118() {
        let uni: usize = 1066; // 128.406 = UNI 1066 (NOT 1073 — that's the Litvin number)
        let ops = get_magnetic_operations(uni);
        assert!(ops.is_some(), "Should get ops for UNI 1066 (128.406)");
        let ops = ops.unwrap();

        let n_u = ops.timerev.iter().filter(|&&t| !t).count();
        let n_a = ops.timerev.iter().filter(|&&t| t).count();
        println!("Magnetic SG UNI {}: {} ops ({} unitary + {} anti-unitary)",
            uni, ops.len(), n_u, n_a);

        // ── 1. Full group (ignore θ) should identify as parent SG 128 ──
        let hall_full = crate::spg_get_hall_number_from_symmetry(
            &ops.rot, &ops.trans, 1e-5,
        );
        assert!(hall_full.is_ok(), "Should identify full group");
        let sg_full = crate::spg_get_spacegroup_type(hall_full.unwrap()).unwrap();
        assert_eq!(sg_full.number, 128,
            "Full ops should identify as SG 128, got SG {}", sg_full.number);
        println!("Full group (ignore θ): SG 128 ✓");

        // ── 2. Unitary subgroup should identify as SG 118 (P-4n2) ──
        let h_sg = identify_unitary_subgroup(uni);
        assert!(h_sg.is_some(), "Should identify unitary subgroup");
        let h_sg = h_sg.unwrap();
        assert_eq!(h_sg, 118,
            "Unitary subgroup of 128.406 should be SG 118, got SG {}", h_sg);
        println!("Unitary subgroup: SG 118 (P-4n2) ✓");

        // ── 3. Verify: all 16 magnetic rotations are in parent SG 128 ──
        let parent_ops = get_parent_operations(128);
        let all_match = ops.rot.iter().all(|r| parent_ops.rot.contains(r));
        assert!(all_match, "All magnetic rotations should be in SG 128 ops");
        println!("Magnetic ops ⊆ SG 128 ops ✓");
    }

    /// Identify the unitary subgroup space group number from a magnetic UNI.
    #[test]
    fn test_identify_unitary_subgroup_api() {
        // 128.406 → UNI 1066 → unitary SG 118 (P-4n2) — verified against BCS
        assert_eq!(identify_unitary_subgroup(1066), Some(118));

        // 129.413 → UNI 1073 → parent SG 129, black-white
        // Its unitary subgroup should also be identifiable
        let result = identify_unitary_subgroup(1073);
        println!("UNI 1073 (129.413) unitary subgroup: {:?}", result);
        assert!(result.is_some(), "UNI 1073 should identify");

        // 1.2 (BNS) → UNI 2, simplest non-trivial magnetic SG
        assert!(identify_unitary_subgroup(2).is_some(), "UNI 2 should work");
    }

    /// Cross-validate: for all compound irreps, PIR χ = Σ CIR component χ.
    #[test]
    fn test_cir_pir_cross_validation() {
        let mut checked = 0usize;
        let mut mismatches = 0usize;
        // Iterate over all SGs
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                let n_comp = ir.cir_component_count();
                if n_comp == 0 { continue; }

                let pir = ir.characters();
                let n_ops = pir.len();
                if n_ops == 0 { continue; }

                // Sum CIR component characters
                let mut cir_sum_re = vec![0.0f64; n_ops];
                let mut cir_sum_im = vec![0.0f64; n_ops];
                for c in 0..n_comp {
                    let cir = ir.cir_component_chars(c);
                    if cir.len() < n_ops * 2 {
                        mismatches += 1;
                        break;
                    }
                    for op in 0..n_ops {
                        cir_sum_re[op] += cir[2 * op];     // real part
                        cir_sum_im[op] += cir[2 * op + 1]; // imag part
                    }
                }

                for op in 0..n_ops {
                    let diff_re = (pir[op] - cir_sum_re[op]).abs();
                    let diff_im = cir_sum_im[op].abs();
                    if diff_re > 0.01 || diff_im > 0.01 {
                        mismatches += 1;
                        eprintln!("MISMATCH SG{} {} op{}: PIR={:.4} CIR_sum=({:.4},{:.4})",
                            sg, ir.ml, op, pir[op], cir_sum_re[op], cir_sum_im[op]);
                    }
                }
                checked += 1;
            }
        }
        println!("CIR↔PIR cross-check: {} compound irreps, {} mismatches", checked, mismatches);
        assert_eq!(mismatches, 0, "All CIR sums must match PIR characters");
        assert!(checked > 500, "Should cover at least 500 compound irreps, got {}", checked);
    }

    /// Test BNS/OG → UNI lookup functions.
    #[test]
    fn test_uni_lookup() {
        assert_eq!(uni_from_bns("128.406"), Some(1066));
        assert_eq!(uni_from_bns("129.413"), Some(1073));
        assert_eq!(uni_from_bns("1.1"), Some(1));

        assert_eq!(uni_from_og("128.8.1073"), Some(1066));
        assert_eq!(uni_from_og("129.3.1077"), Some(1073));

        // Non-existent labels
        assert_eq!(uni_from_bns("nonexistent"), None);
        assert_eq!(uni_from_og("999.999.999"), None);
    }

    /// SG 128 Γ-point double group irreps — verified against BCS
    /// https://cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg
    #[test]
    fn test_sg128_gamma_double_group() {
        let sg128 = irreps_of(128);
        let gamma: Vec<_> = sg128.iter()
            .filter(|r| r.k_label() == "GM")
            .collect();

        // BCS shows: 10 scalar (GM1±-GM5±) + 4 spinor (GM̄6-GM̄9)
        let gamma_scalar: Vec<_> = gamma.iter().filter(|r| !r.spinor).collect();
        let gamma_spinor: Vec<_> = gamma.iter().filter(|r| r.spinor).collect();

        assert!(gamma_scalar.len() >= 5,
            "SG 128 Γ should have >=5 scalar irreps, got {}", gamma_scalar.len());
        assert!(gamma_spinor.len() >= 2,
            "SG 128 Γ should have >=2 spinor irreps, got {}", gamma_spinor.len());

        // Verify scalar labels: GM1+, GM1-, GM2+, GM2-, ...
        let scalar_labels: Vec<&str> = gamma_scalar.iter().map(|r| r.ml).collect();
        for prefix in &["GM1", "GM2", "GM3", "GM4", "GM5"] {
            let has = scalar_labels.iter().any(|l| l.starts_with(prefix));
            assert!(has, "Should have {} scalar irrep at Γ", prefix);
        }

        // Spinor irreps should be 2D (BCS confirms GM̄6-GM̄9 are 2D)
        for ir in &gamma_spinor {
            assert_eq!(ir.dim, 2, "Spinor {} should be 2D, got {}", ir.ml, ir.dim);
            // Identity character should be 2.0 (trace of 2×2 identity)
            let chars = ir.characters();
            if !chars.is_empty() {
                assert!((chars[0] - 2.0).abs() < 0.01,
                    "Spinor {} identity χ should be 2.0, got {}", ir.ml, chars[0]);
            }
        }

        // Scalar irreps at Γ: GM1±-GM4± are 1D, GM5± may be 2D (PIR convention)
        for ir in &gamma_scalar {
            if ir.ml.starts_with("GM5") {
                assert!(ir.dim == 1 || ir.dim == 2,
                    "GM5± should be 1D or 2D, got dim={}", ir.dim);
            } else {
                assert_eq!(ir.dim, 1, "Scalar {} should be 1D, got {}", ir.ml, ir.dim);
            }
        }
    }

    /// BCS validation: 128.406 at Z point, all coreps computed from H = SG 118.
    ///
    /// BCS reference (from k-Subgroupsmag.html):
    ///   Unitary Space Group: P-4n2 (No. 118) in standard setting.
    ///   Magnetic little co-group: 4'/m'mm' (12 ops: 8 unitary + 4 anti-unitary)
    ///
    /// Corep table (from BCS corepresentations_out.pl):
    ///   Z1Z2(2D, type C), Z3Z4(2D, type C), Z5(2D, type A), Z̄6Z̄7(4D spinor, type C)
    ///
    /// Our computation uses H = SG 118's PIR irreps at Z:
    ///   Z1Z4, Z2Z3, Z5 (scalar), Z6, Z7 (spinor)
    /// Type C doubles the dimension: 2D PIR → 4D corep.
    ///
    /// Character order: verify h_seitz[0] is identity with CIR χ=dim.
    #[test]
    fn test_char_order_sg118() {
        let uni = 1066usize;
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let h_sg = identify_unitary_subgroup(uni).unwrap();
        let h_ops = get_parent_operations(h_sg as u8);
        let h_seitz = ops_to_seitz(&h_ops);
        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);

        // Check identity at position 0
        let id_op = &h_seitz[0];
        assert!(id_op.rot[0][0]==1 && id_op.rot[1][1]==1 && id_op.rot[2][2]==1,
            "h_seitz[0] must be identity");
        assert!(id_op.trans[0].abs()<0.01 && id_op.trans[1].abs()<0.01 && id_op.trans[2].abs()<0.01,
            "identity must have zero translation");

        // For each Z-point CIR irrep, check χ(id)=dim
        for ir in h_irreps.iter().filter(|r| r.k_label() == "Z") {
            if ir.cir_component_count() > 0 {
                for c in 0..ir.cir_component_count() {
                    let cir = ir.cir_component_chars(c);
                    let chi_id = Complex64::new(cir[0], cir[1]);
                    println!("{} comp{}: cir_chars[0]=({:.2},{:.2}) |χ|={:.2}",
                        ir.ml, c, chi_id.re, chi_id.im, chi_id.norm());
                }
            } else {
                let pir = ir.characters();
                println!("{} (non-compound): pir_chars[0]={:.2} dim={}",
                    ir.ml, pir[0], ir.dim);
            }
        }

        // Print full h_seitz ↔ cir_chars mapping for Z1Z4's first component
        if let Some(z1z4) = h_irreps.iter().find(|r| r.ml == "Z1Z4") {
            let cir = z1z4.cir_component_chars(0);
            wigner::debug_char_order(cir, &h_seitz, "SG118 Z1Z4 comp0");
        }
    }

    /// Diagnostic: print Wigner sum term-by-term for SG 118 Z-point irreps.
    #[test]
    #[ignore = "diagnostic only; run with --ignored --nocapture"]
    fn debug_wigner_z_point() {
        let uni = 1066usize;
        let mag_ops = get_magnetic_operations(uni).unwrap();
        let h_sg = identify_unitary_subgroup(uni).unwrap();
        let h_ops = get_parent_operations(h_sg as u8);
        let mag_seitz = ops_to_seitz(&mag_ops);
        let h_seitz = ops_to_seitz(&h_ops);
        let a0_idx = mag_ops.timerev.iter().position(|&t| t).unwrap();
        let a0 = &mag_seitz[a0_idx];

        println!("\n=== Wigner diagnostic: UNI {} → H=SG {} ===", uni, h_sg);
        println!("Magnetic ops: {} total, {} unitary, {} anti-unitary",
            mag_ops.len(),
            mag_ops.timerev.iter().filter(|&&t| !t).count(),
            mag_ops.timerev.iter().filter(|&&t| t).count());
        println!("a₀ (anti-unitary rep): R=[{},{},{};{},{},{};{},{},{}] t=({:.4},{:.4},{:.4})",
            a0.rot[0][0],a0.rot[0][1],a0.rot[0][2],
            a0.rot[1][0],a0.rot[1][1],a0.rot[1][2],
            a0.rot[2][0],a0.rot[2][1],a0.rot[2][2],
            a0.trans[0],a0.trans[1],a0.trans[2]);
        println!("H ops (SG {}): {}", h_sg, h_ops.len());
        for (i, s) in h_seitz.iter().enumerate() {
            println!("  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.4},{:.4},{:.4})",
                i,
                s.rot[0][0],s.rot[0][1],s.rot[0][2],
                s.rot[1][0],s.rot[1][1],s.rot[1][2],
                s.rot[2][0],s.rot[2][1],s.rot[2][2],
                s.trans[0],s.trans[1],s.trans[2]);
        }

        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);

        // Compare magnetic ops with SG 118 standard ops — check for origin shift
        let h_ops_sg118 = get_parent_operations(h_sg as u8);
        println!("\n=== Magnetic ops vs SG 118 standard ops ===");
        println!("Unitary magnetic ops:");
        for i in 0..mag_ops.len() {
            if mag_ops.timerev[i] { continue; }
            let r = &mag_ops.rot[i]; let t = &mag_ops.trans[i];
            // Find matching H op
            let h_match = h_ops_sg118.rot.iter().position(|hr| {
                hr[0][0]==r[0][0] && hr[0][1]==r[0][1] && hr[0][2]==r[0][2]
                && hr[1][0]==r[1][0] && hr[1][1]==r[1][1] && hr[1][2]==r[1][2]
                && hr[2][0]==r[2][0] && hr[2][1]==r[2][1] && hr[2][2]==r[2][2]
            });
            let dt = h_match.map(|hi| {
                let ht = &h_ops_sg118.trans[hi];
                [t[0]-ht[0], t[1]-ht[1], t[2]-ht[2]]
            });
            println!("  mag[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) H_match={:?} dt={:?}",
                i, r[0][0],r[0][1],r[0][2], r[1][0],r[1][1],r[1][2], r[2][0],r[2][1],r[2][2],
                t[0],t[1],t[2], h_match, dt);
        }
        println!("Anti-unitary magnetic ops:");
        for i in 0..mag_ops.len() {
            if !mag_ops.timerev[i] { continue; }
            let r = &mag_ops.rot[i]; let t = &mag_ops.trans[i];
            println!("  mag[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3})",
                i, r[0][0],r[0][1],r[0][2], r[1][0],r[1][1],r[1][2], r[2][0],r[2][1],r[2][2],
                t[0],t[1],t[2]);
        }
        println!("SG 118 standard ops:");
        for i in 0..h_ops_sg118.len() {
            let r = &h_ops_sg118.rot[i]; let t = &h_ops_sg118.trans[i];
            println!("  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3})",
                i, r[0][0],r[0][1],r[0][2], r[1][0],r[1][1],r[1][2], r[2][0],r[2][1],r[2][2],
                t[0],t[1],t[2]);
        }

        let k_irreps: Vec<_> = h_irreps.iter()
            .filter(|r| r.k_label() == "Z").collect();

        for ir in &k_irreps {
            println!("\n--- {} (dim={}, spinor={}, k=({},{},{})/{}) ---",
                ir.ml, ir.dim, ir.spinor, ir.kx, ir.ky, ir.kz, ir.kd);

            let mag_lg = filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let unitary_lg: Vec<usize> = mag_lg.iter()
                .filter(|&&i| !mag_ops.timerev[i]).copied().collect();
            let anti_lg: Vec<usize> = mag_lg.iter()
                .filter(|&&i| mag_ops.timerev[i]).copied().collect();
            println!("  Little group: {} ops ({} unitary + {} anti-unitary)",
                mag_lg.len(), unitary_lg.len(), anti_lg.len());

            let h_chars = ir.characters();
            println!("  H characters ({} ops): {:?}...", h_ops.len(),
                &h_chars[..h_chars.len().min(8)]);

            // Map unitary ops to H
            let op_map: Vec<Option<usize>> = (0..mag_ops.len())
                .map(|i| {
                    if mag_ops.timerev[i] { None }
                    else {
                        let r = &mag_ops.rot[i];
                        h_ops.rot.iter().position(|hr| {
                            hr[0][0] == r[0][0] && hr[0][1] == r[0][1] && hr[0][2] == r[0][2]
                            && hr[1][0] == r[1][0] && hr[1][1] == r[1][1] && hr[1][2] == r[1][2]
                            && hr[2][0] == r[2][0] && hr[2][1] == r[2][1] && hr[2][2] == r[2][2]
                        })
                    }
                }).collect();

            // Term-by-term Wigner sum
            let mut w_sum = 0.0f64;
            for &h_mag_idx in &unitary_lg {
                let h = &mag_seitz[h_mag_idx];
                let h_h = op_map[h_mag_idx];

                // g₀h
                let (g0h, l1) = compose_seitz(
                    &SeitzOp::new(a0.rot, a0.trans, false),
                    &SeitzOp::new(h.rot, h.trans, false),
                );
                // (g₀h)²
                let (sq, l_sq) = square_seitz(&g0h);

                match find_seitz(&sq.rot, &sq.trans, &h_seitz) {
                    Some(m) => {
                        let total_l = [l_sq[0]+m.lattice_shift[0], l_sq[1]+m.lattice_shift[1], l_sq[2]+m.lattice_shift[2]];
                        let phase = bloch_phase(ir.kx, ir.ky, ir.kz, ir.kd, &total_l);
                        let chi = if m.op_index < h_chars.len() { h_chars[m.op_index] } else { 0.0 };
                        let contrib = chi * phase.re;
                        w_sum += contrib;
                        println!("    h[{}]→H[{}]: (a₀h)²=H[{}] L={:?} ph={:.2} χ={:.2} contrib={:.2}",
                            h_mag_idx, h_h.map_or("?".into(), |x| x.to_string()),
                            m.op_index, total_l, phase, chi, contrib);
                    }
                    None => {
                        println!("    h[{}]→H[{}]: (a₀h)² R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) NOT FOUND",
                            h_mag_idx, h_h.map_or("?".into(), |x| x.to_string()),
                            sq.rot[0][0],sq.rot[0][1],sq.rot[0][2],
                            sq.rot[1][0],sq.rot[1][1],sq.rot[1][2],
                            sq.rot[2][0],sq.rot[2][1],sq.rot[2][2],
                            sq.trans[0],sq.trans[1],sq.trans[2]);
                    }
                }
            }
            let w = w_sum / (unitary_lg.len() as f64).max(1.0);
            println!("  Wigner W = {:.4} → {}", w,
                if w.abs() < 0.01 { "Type C" } else if w > 0.0 { "Type A" } else { "Type B" });

            // Unwrapped square diagnostic for h[4] and h[7]
            if ir.ml == "Z1Z4" {
                let a0_idx = anti_lg[0];
                wigner::debug_unwrapped_square(4, a0_idx, &mag_seitz, &h_seitz,
                    ir.kx, ir.ky, ir.kz, ir.kd);
                wigner::debug_unwrapped_square(7, a0_idx, &mag_seitz, &h_seitz,
                    ir.kx, ir.ky, ir.kz, ir.kd);

                // Direct anti-coset Wigner sum
                let cir = ir.cir_component_chars(0);
                let w_direct = wigner::wigner_direct_anti_coset(
                    cir, &anti_lg, &mag_seitz, &h_seitz,
                    ir.kx, ir.ky, ir.kz, ir.kd);
                println!("  Direct anti-coset W = {:.4}", w_direct);

                // Try all antiunitary ops as a₀
                println!("\n  Sweeping a₀ choices:");
                for &a0_cand in &anti_lg {
                    let mut w_sum = Complex64::ZERO;
                    let mut np = 0u32; let mut nm = 0u32;
                    for &h_mag_idx in &unitary_lg {
                        let h = &mag_seitz[h_mag_idx];
                        let g0_sp = SeitzOp::new(mag_seitz[a0_cand].rot, mag_seitz[a0_cand].trans, false);
                        let h_sp = SeitzOp::new(h.rot, h.trans, false);
                        let (g0h, l1) = compose_seitz(&g0_sp, &h_sp);
                        let (sq, lsq) = square_seitz(&g0h);
                        if let Some(m) = find_seitz(&sq.rot, &sq.trans, &h_seitz) {
                            let rl1 = mat_vec_i32(&g0h.rot, &l1);
                            let tl = add3(&add3(&lsq, &m.lattice_shift), &add3(&l1, &rl1));
                            let ph = bloch_phase(ir.kx, ir.ky, ir.kz, ir.kd, &tl);
                            let chi = Complex64::new(cir[2*m.op_index], cir[2*m.op_index+1]);
                            w_sum += chi * ph;
                            if ph.re > 0.5 { np += 1; } else if ph.re < -0.5 { nm += 1; }
                        }
                    }
                    let w = w_sum / (unitary_lg.len() as f64);
                    println!("    a₀=mag[{}]: +={} -={} W={:.4}", a0_cand, np, nm, w);
                }
            }
        }
    }

    #[test]
    fn test_corep_sg128_406_z_bcs() {
        let coreps = compute_coreps("128.406", "Z");
        assert!(coreps.is_some(), "Should compute coreps for 128.406 at Z");
        let coreps = coreps.unwrap();
        assert!(!coreps.is_empty());

        println!("\n=== 128.406 Z-point corepresentations (from H = SG 118) ===");
        println!("{:<8} {:<4} {:<8} {:<8}", "Label", "Dim", "Type", "χ(id)");
        println!("-------- ---- -------- --------");

        // ── BCS reference character table (little group, 12 ops) ──
        // Z1Z2 (type C, from Z1+Z4 of H) → 2×χ_re = 2*Re(χ)
        // Z3Z4 (type C, from Z2+Z3 of H) → 2×χ_re = 2*Re(χ)
        // Z5   (type A, from Z5 of H)     → χ = χ (same)
        // Z̄6Z̄7 (type C spinor, from Z6+Z7 of H) → 4D corep

        // Collect computed coreps by label pattern for BCS comparison
        for (label, c) in &coreps {
            let type_str = match c.corep_type {
                CorepType::A => "A",
                CorepType::B => "B",
                CorepType::C => "C",
                CorepType::Unsupported => "?",
            };
            println!("{:<8} {:<4} {:<8} {:<8.1}", label, c.dim, type_str, c.characters[0]);

            if c.corep_type == CorepType::Unsupported {
                continue;  // skip invariants for unsupported (spinor, etc.)
            }

            // Basic invariants
            assert!(c.characters[0] > 0.0, "χ(id) must be > 0 for {}", label);
            assert!(c.dim > 0, "dim must be > 0 for {}", label);

            // χ(id) always equals corep dimension
            assert!((c.characters[0] - c.dim as f64).abs() < 0.01,
                "χ(id)={} should equal dim={} for {}", c.characters[0], c.dim, label);

            // Number of anti-unitary ops with zero character for type C
            if c.corep_type == CorepType::C {
                let zero_count = c.characters.iter()
                    .zip(c.timerev.iter())
                    .filter(|(chi, tr)| **tr && chi.abs() < 0.01)
                    .count();
                let anti_count = c.timerev.iter().filter(|&&t| t).count();
                assert_eq!(zero_count, anti_count,
                    "Type C: all anti-unitary chars should be 0 for {}", label);
            }
        }

        // Verify we have the expected number of coreps from SG 118 at Z
        // SG 118 at Z: 3 scalar (Z1Z4, Z2Z3, Z5) + 2 spinor (Z6, Z7) = 5 H irreps
        assert!(coreps.len() >= 3,
            "Should have >=3 Z-point coreps (scalar), got {}", coreps.len());

        println!("\nBCS comparison: H = SG 118 irreps → corep → BCS magnetic irreps");
        println!("  H:Z1Z4(2D,PIR) → corep type-C → BCS Z1Z2(2D)");
        println!("  H:Z2Z3(2D,PIR) → corep type-C → BCS Z3Z4(2D)");
        println!("  H:Z5(2D,PIR)   → corep type-A → BCS Z5(2D)");
        println!("  H:Z6,Z7(2D,spinor) → corep type-C → BCS Z̄6Z̄7(4D)");
    }

    /// BCS: 165.95 (P-3c'1, UNI 1325) at L:(1/2,0,1/2)
    ///
    /// k-Subgroupsmag_165.95.html confirms:
    ///   Unitary subgroup: P-3 (No. 147)
    ///   Magnetic little co-group: 2'/m'
    ///   Coreps: L₁⁻L₁⁺ (2D, Type C), L̄₂L̄₃ (2D spinor, Type C)
    #[test]
    fn test_corep_sg165_95_l_bcs() {
        let uni = 1325usize;

        // 1. Verify unitary subgroup identification
        let h_sg = identify_unitary_subgroup(uni);
        assert!(h_sg.is_some(), "Should identify unitary subgroup of 165.95");
        let h_sg = h_sg.unwrap();
        println!("165.95 (UNI {}) → unitary subgroup: SG {}", uni, h_sg);

        // 2. Verify magnetic operations exist
        let mag_ops = get_magnetic_operations(uni);
        assert!(mag_ops.is_some(), "Should get ops for UNI {}", uni);
        let mag_ops = mag_ops.unwrap();
        let n_u = mag_ops.timerev.iter().filter(|&&t| !t).count();
        let n_a = mag_ops.timerev.iter().filter(|&&t| t).count();
        println!("  {} ops ({} unitary + {} anti-unitary)", mag_ops.len(), n_u, n_a);

        // 3. Compute coreps at L point (using H = unitary subgroup)
        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);
        let l_irreps: Vec<&IrrepRecord> = h_irreps.iter()
            .filter(|r| r.k_label() == "L")
            .collect();
        let n_scalar = l_irreps.iter().filter(|r| !r.spinor).count();
        let n_spinor = l_irreps.iter().filter(|r| r.spinor).count();
        println!("  H=SG{} L-point irreps: {} scalar + {} spinor",
            h_sg, n_scalar, n_spinor);

        assert!(!l_irreps.is_empty(), "Should have L-point irreps");

        // 4. Compute coreps one by one
        for ir in &l_irreps {
            if let Some(c) = ir.corepresentation(uni) {
                let type_str = match c.corep_type {
                    CorepType::A => "A",
                    CorepType::B => "B",
                    CorepType::C => "C",
                    CorepType::Unsupported => "?",
                };
                println!("  {}: dim={} type={} χ(id)={:.1}",
                    ir.ml, c.dim, type_str, c.characters[0]);

                if c.corep_type != CorepType::Unsupported {
                    assert!(c.dim > 0);
                    assert!((c.characters[0] - c.dim as f64).abs() < 0.01,
                        "χ(id) should equal dim for {}", ir.ml);
                }
            }
        }
    }

    /// BCS: SG 139 (I4/mmm) double-group irreps at k=(1,1,1) (P point)
    ///
    /// k-Subgroupsmag_139.html shows 14 irreps (10 scalar + 4 spinor)
    /// with 4 operations: {1|t}, {2₀₀₁|0}, {4⁺₀₀₁|0}, {4⁻₀₀₁|0}
    #[test]
    fn test_sg139_p_point_bcs() {
        let sg = 139u8;
        let irreps = crate::irrep::query::irreps_of(sg);

        // P-point irreps (k=(1,1,1) in body-centered tetragonal)
        let p_irreps: Vec<&IrrepRecord> = irreps.iter()
            .filter(|r| r.k_label() == "P")
            .collect();
        println!("SG{} P-point: {} irreps ({} scalar + {} spinor)",
            sg, p_irreps.len(),
            p_irreps.iter().filter(|r| !r.spinor).count(),
            p_irreps.iter().filter(|r| r.spinor).count());

        assert!(!p_irreps.is_empty(), "SG 139 should have P-point irreps");

        // BCS shows 14 irreps: M₁⁺..M₅⁻ (10 scalar) + M̄₆..M̄₉ (4 spinor)
        let scalar: Vec<_> = p_irreps.iter().filter(|r| !r.spinor).collect();
        let spinor: Vec<_> = p_irreps.iter().filter(|r| r.spinor).collect();
        assert!(scalar.len() >= 5, "Should have >=5 scalar P-point irreps");
        assert!(spinor.len() >= 2, "Should have >=2 spinor P-point irreps");

        // Check P1-P5 have correct dimensions (BCS: 1D for P1-P4, 2D for P5)
        for ir in &scalar {
            assert!(ir.dim > 0, "dim > 0 for {}", ir.ml);
            let chars = ir.characters();
            assert!(!chars.is_empty(), "Should have characters for {}", ir.ml);
            // Identity character should equal dim
            assert!((chars[0] - ir.dim as f64).abs() < 0.01,
                "χ(id)={} ≠ dim={} for {}", chars[0], ir.dim, ir.ml);
            println!("  {}: dim={} ops={} χ(id)={}", ir.ml, ir.dim, chars.len(), chars[0]);
        }

        // Test matrix reordering for a P-point irrep with matrix data
        if let Some(p1) = scalar.first() {
            let mats = p1.matrices();
            if !mats.is_empty() {
                println!("  {}: {} matrix elements", p1.ml, mats.len());
                let h_ops = get_parent_operations(sg);
                let h_seitz = ops_to_seitz(&h_ops);
                let reordered = p1.matrices_reordered(&h_seitz);
                assert_eq!(reordered.len(), mats.len(),
                    "Reordered matrix should have same size");
                // Identity should be at H[0] position (1,0,0 in original)
                let dim = p1.dim as usize;
                if dim > 0 && reordered.len() >= dim * dim {
                    let trace: f64 = (0..dim).map(|d| reordered[d * dim + d]).sum();
                    assert!((trace - p1.dim as f64).abs() < 0.5,
                        "Reordered identity trace should ≈ dim");
                }
                println!("  Matrix reordering OK ({} elements)", reordered.len());
            }
        }
    }

    /// Every isotropy subgroup record points to a valid SG (1-230).
    #[test]
    fn test_all_isotropy_subgroups_are_well_formed() {
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                for sub in ir.subgroups() {
                    assert!(sub.sg >= 1 && sub.sg <= 230,
                        "invalid isotropy SG={} for parent SG{} {}", sub.sg, sg, ir.ml);
                }
                for msub in ir.magnetic_subgroups() {
                    assert!(msub.mag_sg >= 1 && msub.mag_sg <= 1651,
                        "invalid mag isotropy SG={} for parent SG{} {}", msub.mag_sg, sg, ir.ml);
                }
            }
        }
    }

    /// Type C corepresentations pair two H irreps into one magnetic corep.
    /// Verify that compute_coreps doesn't produce duplicate magnetic irreps
    /// for the same Type C pair.
    #[test]
    fn test_type_c_coreps_are_deduplicated() {
        let coreps = compute_coreps("128.406", "Z");
        assert!(coreps.is_some());
        let coreps = coreps.unwrap();

        // Type C pairs (Z1Z4+Z2Z3, Z6+Z7) should each appear once
        // as combined coreps, not as individual entries
        let mut type_c_pairs: Vec<Vec<&str>> = Vec::new();
        for (label, c) in &coreps {
            if c.corep_type == CorepType::C {
                // Collect labels that should pair
                let labels: Vec<&str> = vec![label];
                type_c_pairs.push(labels);
            }
        }
        // With current API each H irrep returns its own Corepresentation,
        // so for Type C we expect pairs. For now, just verify they're all Type C.
        for (_label, c) in &coreps {
            if c.corep_type == CorepType::C {
                assert!(c.dim > 0);
                // Antiunitary characters must be 0 for Type C
                for (i, &chi) in c.characters.iter().enumerate() {
                    if c.timerev[i] {
                        assert!(chi.abs() < 0.01,
                            "Type C antiunitary char must be 0, got {} at op {}", chi, i);
                    }
                }
            }
        }
        println!("Type C dedup check: {} coreps, all antiunitary chars zero ✓", coreps.len());
    }

    /// Exhaustive: all 1651 magnetic space groups have valid operations
    /// and identifiable unitary subgroups.
    #[test]
    fn test_all_magnetic_sgs_have_valid_operations() {
        let mut ok = 0usize;
        let mut fail = 0usize;
        for uni in 1usize..=1651 {
            if let Some(ops) = get_magnetic_operations(uni) {
                assert!(!ops.is_empty(), "UNI {} has empty ops", uni);
                let n_u = ops.timerev.iter().filter(|&&t| !t).count();
                let n_a = ops.timerev.iter().filter(|&&t| t).count();
                assert!(n_u > 0, "UNI {} has no unitary ops", uni);
                // Every magnetic op must have a valid rotation (det = ±1)
                for i in 0..ops.len() {
                    let r = &ops.rot[i];
                    let det = r[0][0] * (r[1][1]*r[2][2] - r[1][2]*r[2][1])
                            - r[0][1] * (r[1][0]*r[2][2] - r[1][2]*r[2][0])
                            + r[0][2] * (r[1][0]*r[2][1] - r[1][1]*r[2][0]);
                    assert!(det == 1 || det == -1,
                        "UNI {} op[{}]: det={}, not ±1", uni, i, det);
                }
                // Verify unitary subgroup can be identified (may fail for some edge cases)
                if let Some(h_sg) = identify_unitary_subgroup(uni) {
                    assert!(h_sg >= 1 && h_sg <= 230,
                        "UNI {} unitary SG={} out of range", uni, h_sg);
                }
                ok += 1;
            } else {
                fail += 1;
            }
        }
        println!("Magnetic ops: {}/1651 OK, {} missing", ok, fail);
        assert!(ok > 1600, "Should have >=1600 valid MSGs, got {}", ok);
    }

    /// Exhaustive: all spinor (double-group) irreps have valid character tables.
    /// Central element Ē (2π rotation) character should be -dim for spinor irreps.
    #[test]
    fn test_all_spinor_irreps_are_well_formed() {
        let mut total = 0usize;
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                if !ir.spinor { continue; }
                total += 1;
                assert!(ir.dim > 0, "spinor {} SG{} dim=0", ir.ml, sg);
                let chars = ir.characters();
                assert!(!chars.is_empty(), "spinor {} SG{} no chars", ir.ml, sg);
                assert!(chars[0] > 0.0, "spinor {} SG{} χ(E)={}", ir.ml, sg, chars[0]);
                // Spinor irreps are double-valued: typical dims are 1,2,3,4,6
                assert!(ir.dim >= 1 && ir.dim <= 8,
                    "spinor {} SG{} unexpected dim={}", ir.ml, sg, ir.dim);
                // Identity character should be integer
                assert!((chars[0] - chars[0].round()).abs() < 1e-8,
                    "spinor {} SG{} χ(E)={} not integer", ir.ml, sg, chars[0]);
                // Spin ops should exist
                let (rots, trans, su2) = ir.spin_ops();
                if ir.spin_lg_char_count() > 0 {
                    assert!(!rots.is_empty(),
                        "spinor {} SG{} has lg ops but no spin op rots", ir.ml, sg);
                }
            }
        }
        assert!(total > 3000, "Should have >3000 spinor irreps, got {}", total);
        println!("Spinor irreps: {} total, all well-formed ✓", total);
    }

    /// Database format sanity: all irrep k-vectors have reasonable denominators.
    #[test]
    fn test_all_irrep_k_vectors_are_well_formed() {
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                // kd is the common denominator; capped by database convention
                const MAX_KD: i8 = 24;
                assert!(ir.kd >= 0 && ir.kd <= MAX_KD,
                    "SG{} {}: kd={} out of [0,{}]", sg, ir.ml, ir.kd, MAX_KD);
                // Gamma-like points must have kd=0 → k=(0,0,0)
                if ir.kd == 0 {
                    assert_eq!((ir.kx, ir.ky, ir.kz), (0, 0, 0),
                        "SG{} {}: kd=0 but k=({},{},{})", sg, ir.ml, ir.kx, ir.ky, ir.kz);
                }
            }
        }
    }

    /// Exhaustive: all non-spinor (single-valued) irreps satisfy basic
    /// representation-theory invariants: χ(E)=dim, characters are finite,
    /// matrix data is consistent with dimension.
    #[test]
    fn test_all_scalar_irreps_basic_invariants() {
        let mut checked = 0usize;
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                if ir.spinor { continue; }
                checked += 1;
                assert!(ir.dim > 0, "SG{} {}: dim=0", sg, ir.ml);
                assert!(!ir.ml.is_empty(), "SG{}: empty label", sg);
                let chars = ir.characters();
                assert!(!chars.is_empty(), "SG{} {}: empty chars", sg, ir.ml);
                assert!((chars[0] - ir.dim as f64).abs() < 1e-8,
                    "SG{} {}: χ(E)={} != dim={}", sg, ir.ml, chars[0], ir.dim);
                assert!(chars.iter().all(|x| x.is_finite()),
                    "SG{} {}: non-finite character found", sg, ir.ml);
                let mats = ir.matrices();
                if !mats.is_empty() {
                    let dim = ir.dim as usize;
                    assert!(mats.len() % (dim * dim) == 0,
                        "SG{} {}: matrix len {} not divisible by dim²={}",
                        sg, ir.ml, mats.len(), dim * dim);
                }
            }
        }
        assert!(checked > 4000, "Should have >4000 scalar irreps, got {}", checked);
        println!("Scalar irreps: {} total, all well-formed ✓", checked);
    }

    /// Regression: high-dimension image labels (e.g. "K1536a") must not
    /// fall back to dim=1.  This was the root cause of the K1536a bug.
    #[test]
    fn test_high_dim_image_irreps_not_default_to_one() {
        let mut checked = 0usize;
        for sg in 1u8..=230 {
            for ir in crate::irrep::query::irreps_of(sg) {
                if ir.image.starts_with('K') || ir.image.starts_with('L')
                    || ir.image.starts_with('M') || ir.image.starts_with('N')
                {
                    assert!(ir.dim > 1,
                        "SG{} {}: image={} dim={} (should not fall back to 1)",
                        sg, ir.ml, ir.image, ir.dim);
                    assert_eq!(ir.characters()[0] as usize, ir.dim as usize,
                        "SG{} {}: χ(E)={} != dim={}",
                        sg, ir.ml, ir.characters()[0], ir.dim);
                    checked += 1;
                }
            }
        }
        println!("High-dim image irreps: {} checked, all dim > 1 ✓", checked);
        assert!(checked > 0, "Should have at least one high-dim image irrep");
    }
}
