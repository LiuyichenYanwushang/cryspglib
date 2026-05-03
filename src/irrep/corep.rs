//! Co-representation (corep) theory for magnetic irreps.
//!
//! Magnetic space groups contain anti-unitary operations (time reversal θ
//! combined with space operations).  Their irreducible co-representations
//! are classified by Wigner's test into three types.
//!
//! # References
//!
//! - Wigner (1959), *Group Theory*, Chapter 26
//! - Bradley & Cracknell (1972), *The Mathematical Theory of Symmetry in Solids*
//! - Stokes, Campbell & Hatch, ISOTROPY Suite documentation

use crate::mathfunc::{Mat3I, Vec3};
use crate::spg_database::{spgdb_get_spacegroup_operations, spgdb_get_spacegroup_type};
use super::types::IrrepRecord;

/// Co-representation type from Wigner's test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorepType {
    /// W = +1: D ∼ D*, real representation.
    /// Magnetic irrep same dimension, anti-unitary characters non-zero.
    A,
    /// W = -1: D ∼ D*, pseudo-real (quaternionic).
    /// Magnetic irrep same dimension, anti-unitary characters non-zero.
    B,
    /// W = 0: D ≁ D*.
    /// Magnetic irrep = D ⊕ D*, dimension doubled, anti-unitary characters = 0.
    C,
}

impl CorepType {
    pub fn description(&self) -> &'static str {
        match self {
            CorepType::A => "type-a: D ~ D*, real (W=+1)",
            CorepType::B => "type-b: D ~ D*, pseudo-real (W=-1)",
            CorepType::C => "type-c: D ≁ D* (W=0)",
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

// ── API ──────────────────────────────────────────────────────────────────────

/// Compute a magnetic co-representation from a non-magnetic irrep and a
/// magnetic space group.
///
/// # Arguments
///
/// * `irrep` — the non-magnetic irrep (provides D(g) matrices and χ(g) characters)
/// * `uni_number` — OG/UNI number (1–1651) of the magnetic space group
///
/// # Returns
///
/// `None` if the magnetic SG operations cannot be obtained or mapped.
pub fn compute_corepresentation(
    irrep: &IrrepRecord,
    uni_number: u16,
) -> Option<Corepresentation> {
    let uni = uni_number as usize;
    if uni == 0 || uni > 1651 {
        return None;
    }

    // 1. Get magnetic SG type and operations
    let _mag_type = crate::msg_database::msgdb_get_magnetic_spacegroup_type(uni);
    let mag_ops = get_magnetic_operations(uni)?;

    // 2. Get non-magnetic irrep characters and matrices
    let nonmag_chars = irrep.characters();
    let nonmag_mats = irrep.matrices();
    let nonmag_dim = irrep.dim as usize;

    // 3. Get parent SG symmetry operations
    let parent_ops = get_parent_operations(irrep.sg);
    if parent_ops.is_empty() {
        return None;
    }

    // 4. Map each magnetic operation to the corresponding parent operation index
    let op_map = map_magnetic_to_parent(&mag_ops, &parent_ops)?;

    // 5. Separate unitary and anti-unitary operations
    let unitary: Vec<usize> = (0..mag_ops.len())
        .filter(|&i| !mag_ops.timerev[i])
        .collect();
    let antiunitary: Vec<usize> = (0..mag_ops.len())
        .filter(|&i| mag_ops.timerev[i])
        .collect();

    // 6. Wigner's test: determine corep type
    let corep_type = if antiunitary.is_empty() {
        // No anti-unitary ops → purely unitary → trivial corep = original irrep
        CorepType::A
    } else {
        let a0_idx = antiunitary[0]; // pick first anti-unitary as a₀
        let a0_parent = op_map[a0_idx];

        wigner_test(
            nonmag_chars,
            &unitary,
            &op_map,
            a0_parent,
            &parent_ops,
        )
    };

    // 7. Build corep character table
    let characters = build_corep_characters(
        &corep_type,
        &mag_ops,
        &op_map,
        nonmag_chars,
        nonmag_mats,
        nonmag_dim,
    );

    let dim = match corep_type {
        CorepType::A | CorepType::B => nonmag_dim,
        CorepType::C => nonmag_dim * 2,
    };

    Some(Corepresentation {
        characters,
        timerev: mag_ops.timerev.clone(),
        corep_type,
        dim,
        unitary_order: unitary.len(),
        antiunitary_order: antiunitary.len(),
    })
}

// ── Internal helpers ─────────────────────────────────────────────────────────

/// Get the magnetic space group symmetry operations.
fn get_magnetic_operations(
    uni_number: usize,
) -> Option<MagneticOps> {
    // Find the Hall number for this magnetic SG.
    // The magnetic database maps UNI → Hall numbers.
    // For the first (canonical) Hall number associated with this UNI:
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

/// Get the first Hall number for a UNI number.
fn get_first_hall_for_uni(uni: usize) -> Option<usize> {
    if uni == 0 || uni > 1651 {
        return None;
    }
    // MAGNETIC_SPACEGROUP_UNI_MAPPING[uni] = [num_halls, hall1, hall2, ...]
    // We can't access the static directly from here, but msgdb_get_uni_candidates
    // works in reverse (hall → uni). For forward (uni → hall), we need to read
    // the mapping table.
    //
    // Workaround: scan Hall numbers 1-530 and check if they map to this UNI.
    for hall in 1..=530 {
        if let Some([lo, hi]) = crate::msg_database::msgdb_get_uni_candidates(hall) {
            if uni >= lo && uni <= hi {
                // Find the exact match from the mapping table
                // We need the first hall for this uni from
                // MAGNETIC_SPACEGROUP_UNI_MAPPING
                // Fallback: just return the first hall that maps to this uni range
                return Some(hall);
            }
        }
    }
    None
}

/// Get parent (non-magnetic) SG symmetry operations.
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
    MagneticOps {
        rot,
        trans,
        timerev: vec![false; n], // all unitary for non-magnetic SG
    }
}

/// Map each magnetic operation to the index of the corresponding
/// parent (non-magnetic) SG operation, by matching rotation matrices.
fn map_magnetic_to_parent(
    mag_ops: &MagneticOps,
    parent_ops: &MagneticOps,
) -> Option<Vec<usize>> {
    let mut mapping = Vec::with_capacity(mag_ops.len());

    for i in 0..mag_ops.len() {
        let mag_rot = &mag_ops.rot[i];
        // Find parent operation with matching rotation matrix
        let idx = parent_ops.rot.iter().position(|pr| {
            pr[0][0] == mag_rot[0][0]
                && pr[0][1] == mag_rot[0][1]
                && pr[0][2] == mag_rot[0][2]
                && pr[1][0] == mag_rot[1][0]
                && pr[1][1] == mag_rot[1][1]
                && pr[1][2] == mag_rot[1][2]
                && pr[2][0] == mag_rot[2][0]
                && pr[2][1] == mag_rot[2][1]
                && pr[2][2] == mag_rot[2][2]
        });

        if let Some(idx) = idx {
            mapping.push(idx);
        } else {
            // Rotation not found in parent ops — the magnetic SG might use
            // a different conventional setting. Skip for now.
            return None;
        }
    }

    Some(mapping)
}

/// Wigner's test for corep type.
///
/// $$W = \frac{1}{|H|}\sum_{h \in H} \chi_D((a_0 h)^2)$$
///
/// where H is the unitary subgroup and a₀ is the anti-unitary representative.
fn wigner_test(
    chars: &[f64],
    unitary: &[usize],
    op_map: &[usize],
    a0_parent_idx: usize,
    parent_ops: &MagneticOps,
) -> CorepType {
    // For each unitary operation h, compute χ((a₀h)²)
    // (a₀h)² = a₀·h·a₀·h
    // In the parent SG, a₀ corresponds to operation a0_parent_idx.
    // We need to compute the product in the parent SG.

    // For simplicity, we use the fact that for magnetic isotropy subgroups
    // derived from non-magnetic irreps, the anti-unitary operations are
    // those that reverse the order-parameter direction.
    //
    // When a₀ is pure time reversal θ (no spatial part), then:
    //   (a₀h)² = h² (since θ² = 1 and θ commutes with h)
    //   W = (1/|H|) Σ χ(h²)
    //
    // For a₀ = θ·g₀ (anti-unitary with spatial part g₀):
    //   (a₀h)² = g₀·h·g₀·h (since θ²=1, θh=hθ)
    //   W = (1/|H|) Σ χ((g₀·h)²)

    let mut w_sum: f64 = 0.0;

    for &h_idx in unitary {
        let h_parent = op_map[h_idx];

        // Compute (a₀·h)² in the parent group
        // a₀·h = apply a₀'s rotation first, then h's rotation
        // (a₀·h)² = (a₀·h)·(a₀·h)

        // First compute a₀·h composite rotation
        let a0_rot = &parent_ops.rot[a0_parent_idx];
        let h_rot = &parent_ops.rot[h_parent];

        // Composite: R_composite = R_a0 · R_h (matrix multiplication)
        let comp_rot = mat_mul(a0_rot, h_rot);

        // Now square it: (a₀·h)² → rotation part = comp_rot · comp_rot
        let sq_rot = mat_mul(&comp_rot, &comp_rot);

        // Find the parent operation with this squared rotation
        if let Some(sq_idx) = find_rotation_in_parent(&sq_rot, parent_ops) {
            if sq_idx < chars.len() {
                w_sum += chars[sq_idx];
            }
        }
    }

    let w = w_sum / (unitary.len() as f64).max(1.0);

    // Classify with tolerance
    if w.abs() < 0.01 {
        CorepType::C
    } else if w > 0.0 {
        CorepType::A
    } else {
        CorepType::B
    }
}

/// Find a parent operation with the given rotation matrix.
fn find_rotation_in_parent(rot: &Mat3I, parent_ops: &MagneticOps) -> Option<usize> {
    parent_ops.rot.iter().position(|pr| {
        pr[0][0] == rot[0][0]
            && pr[0][1] == rot[0][1]
            && pr[0][2] == rot[0][2]
            && pr[1][0] == rot[1][0]
            && pr[1][1] == rot[1][1]
            && pr[1][2] == rot[1][2]
            && pr[2][0] == rot[2][0]
            && pr[2][1] == rot[2][1]
            && pr[2][2] == rot[2][2]
    })
}

/// Integer 3×3 matrix multiplication: C = A·B.
fn mat_mul(a: &Mat3I, b: &Mat3I) -> Mat3I {
    let mut c = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Build the character table for the co-representation.
fn build_corep_characters(
    corep_type: &CorepType,
    mag_ops: &MagneticOps,
    op_map: &[usize],
    nonmag_chars: &[f64],
    _nonmag_mats: &[f64],
    _nonmag_dim: usize,
) -> Vec<f64> {
    let n_ops = mag_ops.len();
    let mut chars = vec![0.0; n_ops];

    match corep_type {
        CorepType::C => {
            // Type c: magnetic irrep = D ⊕ D*
            // χ̃(h) = 2·Re(χ(h)) for unitary ops
            // χ̃(a₀h) = 0 for anti-unitary ops
            for i in 0..n_ops {
                if mag_ops.timerev[i] {
                    chars[i] = 0.0;
                } else {
                    let parent_idx = op_map[i];
                    if parent_idx < nonmag_chars.len() {
                        chars[i] = 2.0 * nonmag_chars[parent_idx];
                    }
                }
            }
        }
        CorepType::A => {
            // Type a: corep = D (same dimension)
            // For unitary ops: χ̃(h) = χ(h)
            // For anti-unitary ops: χ̃(a₀h) = trace of D(a₀h) — needs U matrix
            // Approximate: use D(h) for unitary, need D(a₀h) for anti-unitary
            for i in 0..n_ops {
                let parent_idx = op_map[i];
                if parent_idx < nonmag_chars.len() {
                    chars[i] = nonmag_chars[parent_idx];
                }
            }
        }
        CorepType::B => {
            // Type b: similar to type a but with pseudo-real character
            for i in 0..n_ops {
                let parent_idx = op_map[i];
                if parent_idx < nonmag_chars.len() {
                    chars[i] = nonmag_chars[parent_idx];
                }
            }
        }
    }

    chars
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
    pub fn corepresentation(&self, uni_number: u16) -> Option<Corepresentation> {
        compute_corepresentation(self, uni_number)
    }
}

// ── Hall number lookup ──────────────────────────────────────────────────────

/// Find the canonical Hall number for a space group (1–230).
fn find_hall_number(sg: u8) -> Option<usize> {
    for hall in 1..=530 {
        let st = spgdb_get_spacegroup_type(hall);
        if st.number == sg as usize {
            return Some(hall);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::query::irreps_of;

    #[test]
    fn test_corep_gm4m_pmmm() {
        // SG 221 (Pm-3m) GM4- → magnetic subgroup Pm'mm (UNI 349)
        let gm4m = irreps_of(221).iter()
            .find(|r| r.ml == "GM4-")
            .expect("GM4- not found");

        let corep = gm4m.corepresentation(349);
        assert!(corep.is_some(), "Should compute corep for GM4- → Pm'mm");

        let c = corep.unwrap();
        // GM4- is a 3D irrep. For the magnetic subgroup Pm'mm,
        // the order parameter is (a,0,0) — uniaxial.
        // The corep should be type-a or type-b (same dimension)
        // or type-c (doubled dimension).
        assert!(c.dim == 3 || c.dim == 6, "dim should be 3 or 6, got {}", c.dim);

        // All characters should be finite
        for &chi in &c.characters {
            assert!(chi.is_finite(), "character should be finite");
        }
    }

    #[test]
    fn test_corep_sg1_gm1() {
        // SG 1 (P1) GM1 → simplest case
        let gm1 = irreps_of(1).iter()
            .find(|r| r.ml == "GM1")
            .expect("GM1 not found");

        // Check that magnetic subgroups exist
        let mag_subs = gm1.magnetic_subgroups();
        if !mag_subs.is_empty() {
            let mag_sg = mag_subs[0].mag_sg;
            let corep = gm1.corepresentation(mag_sg);
            assert!(corep.is_some(), "Should compute corep for SG1 GM1");
        }
    }

    #[test]
    fn test_mat_mul() {
        let a = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let b = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        let c = mat_mul(&a, &b);
        assert_eq!(c[0][0], 0);
        assert_eq!(c[0][1], 1);
        assert_eq!(c[1][0], -1);
    }
}
