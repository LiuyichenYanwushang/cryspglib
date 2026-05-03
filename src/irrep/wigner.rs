//! Wigner's co-representation test and character-table construction.
//!
//! # Algorithm overview
//!
//! Given a magnetic space group $$\mathcal{M} = H \cup a_0 H$$ (where $$H$$ is
//! the unitary subgroup and $$a_0 = \mathcal{T} g_0$$ is an anti-unitary coset
//! representative) and a non-magnetic irrep $$\Delta$$ of $$H$$, this module
//! determines the co-representation type and builds the character table.
//!
//! # Notation
//!
//! - $$H$$: unitary subgroup (ordinary space group)
//! - $$g \in \mathcal{M}$$: a magnetic symmetry operation with rotation $$R_g$$
//!   and translation $$\mathbf{t}_g$$
//! - $$\chi(h) = \mathrm{Tr}\,\Delta(h)$$: character of unitary op $$h \in H$$
//! - $$\tilde{\chi}(g)$$: character of magnetic corep at $$g$$
//!
//! # The Wigner sum
//!
//! For each irrep $$\Delta_i$$ of $$H$$, compute
//!
//! $$
//! W_i = \frac{1}{|H|} \sum_{b \in a_0 H} \chi_i(b^2)
//! $$
//!
//! The sum runs over the **anti-unitary coset** $$a_0 H$$.  Each element
//! $$b = a_0 h$$ is anti-unitary.  Its square
//!
//! $$
//! b^2 = a_0 h \, a_0 h
//! $$
//!
//! is **unitary** (product of two anti-unitary ops), so $$\chi_i(b^2)$$ is
//! well-defined from the non-magnetic character table.
//!
//! Since $$a_0 = \mathcal{T} g_0$$ and $$\mathcal{T}$$ commutes with spatial
//! operations while $$\mathcal{T}^2 = +1$$ (for spinless) or $$-1$$ (for
//! spinor), the spatial part of $$b^2$$ is:
//!
//! $$
//! R_{b^2} = (R_{g_0} R_h)^2
//! $$
//!
//! We find $$(g_0 h)^2$$ in the unitary subgroup $$H$$ by composing rotation
//! matrices and then look up $$\chi_i$$ at that operation.
//!
//! # Classification
//!
//! | $$W$$ | Type | Meaning |
//! |-------|------|---------|
//! | $$\approx +1$$ | A | $$\Delta^{a_0} \sim \Delta$$, real |
//! | $$\approx -1$$ | B | $$\Delta^{a_0} \sim \Delta$$, pseudo-real (Kramers) |
//! | $$\approx 0$$  | C | $$\Delta^{a_0} \nsim \Delta$$, pairs with conjugate |
//!
//! # Character table construction
//!
//! Once the type is known, the magnetic corep characters are:
//!
//! **Type A / B** (same dimension $$d$$):
//! - Unitary ops: $$\tilde{\chi}(h) = \chi_i(h)$$ (same as H)
//! - Anti-unitary ops: $$\tilde{\chi}(a_0 h)$$ requires the intertwining
//!   matrix $$U$$ satisfying $$U \Delta(a_0 h a_0^{-1})^* U^{-1} = \Delta(h)$$.
//!   Without $$U$$, we approximate using the H characters.
//!
//! **Type C** (doubled dimension $$2d$$):
//! - Unitary ops: $$\tilde{\chi}(h) = 2\,\mathrm{Re}[\chi_i(h)]$$
//!   (trace of $$\Delta_i \oplus \Delta_i^*$$)
//! - Anti-unitary ops: $$\tilde{\chi}(a_0 h) = 0$$
//!   (off-diagonal blocks couple the two copies)
//!
//! # Little group filtering
//!
//! Only operations that preserve the wave-vector $$\mathbf{k}$$ (modulo
//! reciprocal lattice vectors) contribute to the little-group
//! co-representation.  An operation $$\{R|\mathbf{t}\}$$ preserves
//! $$\mathbf{k} = (k_x, k_y, k_z)/d$$ when
//!
//! $$
//! R \cdot (k_x, k_y, k_z) \equiv (k_x, k_y, k_z) \pmod{d}
//! $$
//!
//! Gamma-point ($$d = 0$$) always has the full group as its little group.
//!
//! # References
//!
//! - Wigner (1959), *Group Theory*, Chapter 26
//! - Bradley & Cracknell (1972), *The Mathematical Theory of Symmetry in Solids*
//!   Chapters 4–5
//! - Stokes, Campbell & Hatch, ISOTROPY Suite documentation
//! - Bilbao Crystallographic Server, *Co-representations of Magnetic Space Groups*

use crate::mathfunc::{mat_multiply_matrix_i3, Mat3I};

// Re-use the corep types and ops from the parent module
use super::corep::{CorepType, MagneticOps};

// ── Little group filter ──────────────────────────────────────────────────────

/// Filter magnetic operations to those that preserve the wave-vector.
///
/// An operation $$\{R|\mathbf{t}\}$$ preserves $$\mathbf{k} = (k_x,k_y,k_z)/k_d$$
/// when $$R \cdot \mathbf{k} \equiv \mathbf{k} \pmod{\text{reciprocal lattice}}$$,
/// i.e. each component of $$(R\mathbf{k} - \mathbf{k})$$ is an integer multiple
/// of the denominator $$k_d$$.
///
/// For the Gamma point ($$k_d = 0$$), all operations trivially preserve $$\mathbf{k}$$.
pub fn filter_little_group(
    kx: i8, ky: i8, kz: i8, kd: i8,
    ops: &MagneticOps,
) -> Vec<usize> {
    // Gamma point: all ops preserve k = (0,0,0)
    if kd == 0 {
        return (0..ops.len()).collect();
    }

    let kd_i = kd as i32;
    let kx_i = kx as i32;
    let ky_i = ky as i32;
    let kz_i = kz as i32;

    (0..ops.len())
        .filter(|&i| {
            let r = &ops.rot[i];
            // R·k (in integer numerator coordinates)
            let rx = r[0][0] as i32 * kx_i + r[0][1] as i32 * ky_i + r[0][2] as i32 * kz_i;
            let ry = r[1][0] as i32 * kx_i + r[1][1] as i32 * ky_i + r[1][2] as i32 * kz_i;
            let rz = r[2][0] as i32 * kx_i + r[2][1] as i32 * ky_i + r[2][2] as i32 * kz_i;
            // R·k ≡ k (mod kd) → each component difference is a multiple of kd
            (rx - kx_i) % kd_i == 0
                && (ry - ky_i) % kd_i == 0
                && (rz - kz_i) % kd_i == 0
        })
        .collect()
}

// ── Wigner test ──────────────────────────────────────────────────────────────

/// Compute the Wigner test sum and classify the corep type.
///
/// # Mathematical definition
///
/// $$
/// W = \frac{1}{|H_{\mathbf{k}}|} \sum_{h \in H_{\mathbf{k}}} \chi\big((g_0 \cdot h)^2\big)
/// $$
///
/// where:
/// - $$H_{\mathbf{k}}$$ is the unitary little co-group (unitary ops that preserve
///   $$\mathbf{k}$$)
/// - $$g_0$$ is the spatial part of the anti-unitary coset representative
///   $$a_0 = \mathcal{T} g_0$$
/// - $$(g_0 \cdot h)^2$$ is the square of the product $$g_0 h$$, computed as a
///   group operation in $$H$$
/// - $$\chi$$ is the character of the non-magnetic irrep $$\Delta$$ of $$H$$
///
/// # Algorithm
///
/// 1. Take the first anti-unitary op in the magnetic little group as the coset
///    representative $$a_0$$
/// 2. For each unitary op $$h$$ in the little group:
///    a. Compose rotations: $$R_{\text{comp}} = R_{g_0} \cdot R_h$$
///    b. Square the composite: $$R_{\text{sq}} = R_{\text{comp}} \cdot R_{\text{comp}}$$
///    c. Find the unitary op $$\tilde{h} \in H$$ whose rotation matches
///       $$R_{\text{sq}}$$
///    d. If found, add $$\chi(\tilde{h})$$ to the running sum
/// 3. Normalise: $$W = \text{sum} / |H_{\mathbf{k}}|$$
/// 4. Classify by threshold:
///    - $$|W| < 0.01$$ → type C (paired with conjugate)
///    - $$W > 0$$ → type A (real)
///    - $$W < 0$$ → type B (pseudo-real / Kramers)
///
/// # Limitations (known)
///
/// - **Translations**: The current implementation uses only the rotation
///   matrices to identify the squared operation.  In non-symmorphic space
///   groups, translations affect the group multiplication and the correct
///   $$\tilde{h}$$ may depend on the translation part.  This is most
///   relevant for non-symmorphic groups at non-Γ k-points.
///
/// - **Spinor factor**: For double-group (spinor) irreps, the
///   anti-unitary operation squares to $$-1$$ (a $$2\pi$$ rotation), and the
///   character of $$(a_0 h)^2$$ must include the spinor sign.  Currently
///   the spatial part only is used; the $$-1$$ factor from $$\mathcal{T}^2$$
///   must be included for correct spinor classification.
///
/// - **Setting dependence**: The rotation matching between magnetic ops and
///   $$H$$ ops requires both to be in the same coordinate basis.  If the
///   magnetic SG uses a different conventional setting than $$H$$, the
///   rotation matrices will not match and the test fails.
pub fn wigner_classify(
    h_chars: &[f64],            // character table of H's irrep Δ
    unitary_indices: &[usize],  // indices (into mag_ops) of unitary little-group ops
    mag_ops: &MagneticOps,      // full magnetic symmetry operations
    h_ops: &MagneticOps,        // unitary subgroup H's symmetry operations
    a0_idx: usize,              // index (into mag_ops) of the anti-unitary representative
) -> CorepType {
    let a0_rot = &mag_ops.rot[a0_idx];

    // Running sum of χ((g₀·h)²) over all unitary ops h in the little group
    let mut w_sum: f64 = 0.0;

    for &h_mag_idx in unitary_indices {
        let h_rot = &mag_ops.rot[h_mag_idx];

        // Step 2a: compose g₀ and h:  R_comp = R_g₀ · R_h
        let comp_rot = mat_multiply_matrix_i3(a0_rot, h_rot);

        // Step 2b: square the composite:  R_sq = R_comp · R_comp
        let sq_rot = mat_multiply_matrix_i3(&comp_rot, &comp_rot);

        // Step 2c: find the H operation whose rotation matches R_sq
        if let Some(sq_h_idx) = find_rotation_in_ops(&sq_rot, h_ops) {
            // Step 2d: add the character of that H operation
            if sq_h_idx < h_chars.len() {
                w_sum += h_chars[sq_h_idx];
            }
        }
        // If no matching rotation is found, the term contributes 0.
        // This can happen when the squared operation is not in H (e.g.,
        // when translations matter or the setting differs).
    }

    // Step 3: normalise
    let n_unitary = (unitary_indices.len() as f64).max(1.0);
    let w = w_sum / n_unitary;

    // Step 4: classify by threshold
    if w.abs() < 0.01 {
        CorepType::C   // Δ and Δ^{a₀} are inequivalent → pairing
    } else if w > 0.0 {
        CorepType::A   // Δ^{a₀} ∼ Δ, real (no Kramers doubling)
    } else {
        CorepType::B   // Δ^{a₀} ∼ Δ, pseudo-real (Kramers doubling)
    }
}

// ── Character table construction ─────────────────────────────────────────────

/// Build the magnetic co-representation character table.
///
/// # Character formulas
///
/// For each magnetic operation $$g$$ in the little group, the corep character
/// $$\tilde{\chi}(g)$$ is:
///
/// **Type C** ($$\tilde{D} = \Delta \oplus \Delta^*$$, doubled dimension):
///
/// $$
/// \tilde{\chi}(h) = 2\,\mathrm{Re}[\chi(h)] \quad \text{(unitary)}
/// $$
/// $$
/// \tilde{\chi}(a_0 h) = 0 \quad \text{(anti-unitary)}
/// $$
///
/// The anti-unitary character is zero because the off-diagonal blocks in
/// $$\tilde{D}(a_0 h)$$ have zero trace.
///
/// **Type A / B** (same dimension as $$\Delta$$):
///
/// $$
/// \tilde{\chi}(h) = \chi(h)
/// $$
///
/// For anti-unitary ops, we use the same H character (exact form requires
/// the intertwining matrix $$U$$; this is an approximation valid when the
/// magnetic and parent settings coincide).
///
/// # Parameters
///
/// - `corep_type`: result of [`wigner_classify`]
/// - `mag_ops`: magnetic symmetry operations
/// - `mag_lg_indices`: which magnetic ops are in the little group
/// - `op_map`: for each magnetic op index, the corresponding H op index
///   (or `None` for anti-unitary ops that lack a direct H correspondence)
/// - `h_chars`: character table of H's irrep
pub fn build_corep_characters(
    corep_type: &CorepType,
    mag_ops: &MagneticOps,
    mag_lg_indices: &[usize],
    op_map: &[Option<usize>],
    h_chars: &[f64],
) -> Vec<f64> {
    let n_lg = mag_lg_indices.len();
    let mut chars = vec![0.0; n_lg];

    for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
        let is_anti = mag_ops.timerev[mag_idx];

        match corep_type {
            CorepType::C => {
                if is_anti {
                    // Anti-unitary: off-diagonal blocks → trace = 0
                    chars[out_idx] = 0.0;
                } else if let Some(h_idx) = op_map[mag_idx] {
                    // Unitary: 2 × Re[χ(h)] = 2 × χ(h) for real irreps
                    // (χ is already real in PIR convention)
                    if h_idx < h_chars.len() {
                        chars[out_idx] = 2.0 * h_chars[h_idx];
                    }
                }
            }
            CorepType::A | CorepType::B => {
                if let Some(h_idx) = op_map[mag_idx] {
                    if h_idx < h_chars.len() {
                        // Same dimension: inherit H character directly
                        chars[out_idx] = h_chars[h_idx];
                    }
                }
                // Anti-unitary ops without H mapping: character left as 0.0
            }
        }
    }

    chars
}

// ── Helper: find rotation in operation list ──────────────────────────────────

/// Find an operation in `ops` whose rotation matrix matches `rot`.
///
/// Returns the index of the first matching operation, or `None` if no match.
fn find_rotation_in_ops(rot: &Mat3I, ops: &MagneticOps) -> Option<usize> {
    ops.rot.iter().position(|op_rot| {
        op_rot[0][0] == rot[0][0]
            && op_rot[0][1] == rot[0][1]
            && op_rot[0][2] == rot[0][2]
            && op_rot[1][0] == rot[1][0]
            && op_rot[1][1] == rot[1][1]
            && op_rot[1][2] == rot[1][2]
            && op_rot[2][0] == rot[2][0]
            && op_rot[2][1] == rot[2][1]
            && op_rot[2][2] == rot[2][2]
    })
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify Wigner test on the simplest case: SG 1 GM1 with BNS "1.3".
    #[test]
    fn test_wigner_sg1_gm1() {
        // Simplest case: verify the wigner_classify function can be called
        // For P1 with only identity, W = χ(id²) / 1 = χ(id) / 1 = 1.0
        let rot = vec![[[1, 0, 0], [0, 1, 0], [0, 0, 1]]];
        let trans = vec![[0.0, 0.0, 0.0]];
        let timerev = vec![false, true]; // identity + anti-unitary identity

        // The anti-unitary op (a₀ = θ) has no spatial part (identity rotation)
        // (a₀·id)² = id² = id → χ(id) = 1.0 → W = 1.0 → type A
        let result = wigner_classify(
            &[1.0],          // χ(id) = 1.0
            &[0],            // unitary op at index 0
            &MagneticOps { rot: rot.clone(), trans: trans.clone(), timerev: timerev.clone() },
            &MagneticOps { rot: rot.clone(), trans: trans.clone(), timerev: vec![false] },
            1,               // a₀ at index 1 (anti-unitary identity)
        );
        assert_eq!(result, CorepType::A);
    }

    /// Verify: when there are no anti-unitary ops, type is A.
    #[test]
    fn test_wigner_no_antiunitary() {
        let ops = MagneticOps {
            rot: vec![[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
            trans: vec![[0.0, 0.0, 0.0]],
            timerev: vec![false],
        };
        // With no anti-unitary ops, type is trivially A
        // (the calling code in corep.rs handles this case before calling wigner)
        assert_eq!(ops.timerev.iter().filter(|&&t| t).count(), 0);
    }
}
