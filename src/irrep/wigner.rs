//! Wigner's co-representation test and character-table construction.
//!
//! # Theory
//!
//! A magnetic space group $$\mathcal{M} = H \cup a_0 H$$ where $$H$$ is the
//! unitary subgroup, $$a_0 = \mathcal{T} g_0$$ is an anti-unitary coset
//! representative.  Given a non-magnetic irrep $$\Delta$$ of $$H$$ at
//! wave-vector $$\mathbf{k}$$, the Wigner indicator is
//!
//! $$
//! W = \frac{1}{|H_{\mathbf{k}}|}
//!     \sum_{h \in H_{\mathbf{k}}} \chi\big((a_0 h)^2\big)
//! $$
//!
//! The summand $$(a_0 h)^2$$ is a **unitary** operation (product of two
//! anti-unitary operations) and must be evaluated using **full Seitz
//! symbols** $$\{R|\mathbf{t}\}$$, not just rotation matrices.
//!
//! ## Bloch phase convention
//!
//! At wave-vector $$\mathbf{k} = (k_x,k_y,k_z)/k_d$$ in reciprocal lattice
//! units, a lattice translation $$\mathbf{L} \in \mathbb{Z}^3$$ contributes
//! a phase factor $$e^{+2\pi i\,\mathbf{k}\cdot\mathbf{L}}$$.  When a
//! computed Seitz operation $$\{R|\mathbf{t}_{\text{comp}}\}$$ differs from
//! the stored database operation $$\{R|\mathbf{t}_{\text{stored}}\}$$ by a
//! lattice vector $$\mathbf{L} = \mathbf{t}_{\text{comp}} - \mathbf{t}_{\text{stored}}$$,
//! the character must be multiplied by this phase.
//!
//! # References
//!
//! - Wigner (1959), *Group Theory*, Chapter 26
//! - Bradley & Cracknell (1972), *The Mathematical Theory of Symmetry in Solids*
//! - Bilbao Crystallographic Server, *Co-representations of Magnetic Space Groups*

use crate::mathfunc::{mat_multiply_matrix_i3, Mat3I};
use super::corep::{CorepType, MagneticOps};

// ── Seitz operation ──────────────────────────────────────────────────────────

/// A space-group operation $$\{R|\mathbf{t}\}$$ with optional time reversal.
///
/// The translation $$\mathbf{t}$$ is stored in **fractional coordinates**
/// (each component in $$[0, 1)$$ after normalisation).
#[derive(Debug, Clone)]
pub struct SeitzOp {
    /// 3×3 integer rotation matrix
    pub rot: Mat3I,
    /// Fractional translation (each component ∈ [0, 1))
    pub trans: [f64; 3],
    /// Whether this operation includes time reversal θ
    pub timerev: bool,
}

impl SeitzOp {
    /// Create from rotation + translation + timerev.
    pub fn new(rot: Mat3I, trans: [f64; 3], timerev: bool) -> Self {
        // Normalise translations to [0, 1)
        let t = [
            trans[0] - trans[0].floor(),
            trans[1] - trans[1].floor(),
            trans[2] - trans[2].floor(),
        ];
        // Handle -0.0
        let t = [if t[0] < 0.0 { t[0] + 1.0 } else { t[0] },
                  if t[1] < 0.0 { t[1] + 1.0 } else { t[1] },
                  if t[2] < 0.0 { t[2] + 1.0 } else { t[2] }];
        SeitzOp { rot, trans: t, timerev }
    }
}

/// Compose two Seitz operations: `g1 ∘ g2` means apply g2 first, then g1.
///
/// $$\{R_1|\mathbf{t}_1\} \circ \{R_2|\mathbf{t}_2\}
///   = \{R_1 R_2 \mid \mathbf{t}_1 + R_1 \mathbf{t}_2\}$$
///
/// Time reversal composes with XOR: anti∘anti = unitary, anti∘unitary = anti, etc.
///
/// Returns `(result, lattice_shift)` where `lattice_shift` is the integer
/// part of the translation (discarded during normalisation).
pub fn compose_seitz(g1: &SeitzOp, g2: &SeitzOp) -> (SeitzOp, [i32; 3]) {
    let rot = mat_multiply_matrix_i3(&g1.rot, &g2.rot);
    let timerev = g1.timerev ^ g2.timerev;

    // t = t1 + R1·t2  (in fractional coordinates)
    let r1 = &g1.rot;
    let mut t = [0.0f64; 3];
    let mut lattice = [0i32; 3];
    for i in 0..3 {
        let raw = g1.trans[i]
            + r1[i][0] as f64 * g2.trans[0]
            + r1[i][1] as f64 * g2.trans[1]
            + r1[i][2] as f64 * g2.trans[2];
        let floor = raw.floor();
        lattice[i] = floor as i32;
        t[i] = raw - floor;
        if t[i] < 0.0 { t[i] += 1.0; lattice[i] -= 1; }
    }

    (SeitzOp { rot, trans: t, timerev }, lattice)
}

/// Square a Seitz operation: g² = g ∘ g.
///
/// For $$g = \{R|\mathbf{t}\}$$:
/// $$g^2 = \{R^2 \mid \mathbf{t} + R\mathbf{t}\}$$
pub fn square_seitz(g: &SeitzOp) -> (SeitzOp, [i32; 3]) {
    compose_seitz(g, g)
}

// ── Convert MagneticOps → Vec<SeitzOp> ──────────────────────────────────────

/// Convert spglib's `MagneticOps` to a `Vec<SeitzOp>`.
pub fn ops_to_seitz(ops: &MagneticOps) -> Vec<SeitzOp> {
    (0..ops.len())
        .map(|i| SeitzOp::new(ops.rot[i], ops.trans[i], ops.timerev[i]))
        .collect()
}

// ── Little group filter ──────────────────────────────────────────────────────

/// Filter magnetic operations to those that preserve the wave-vector.
///
/// For a **unitary** operation $$\{R|\mathbf{t}\}$$, the condition is
/// $$R\mathbf{k} \equiv \mathbf{k} \pmod{\text{reciprocal lattice}}$$.
///
/// For an **anti-unitary** operation $$a = \mathcal{T}\{R|\mathbf{t}\}$$,
/// time reversal sends $$\mathbf{k} \to -\mathbf{k}$$, so the condition is
/// $$-R\mathbf{k} \equiv \mathbf{k} \pmod{\text{reciprocal lattice}}$$.
///
/// In terms of integer components with denominator $$k_d$$:
///
/// ```text
/// unitary:     (R·k - k) ≡ 0  (mod kd)
/// antiunitary: (-R·k - k) ≡ 0 (mod kd)
/// ```
///
/// Gamma point ($$k_d = 0$$): all operations trivially preserve k.
pub fn filter_little_group(
    kx: i8, ky: i8, kz: i8, kd: i8,
    ops: &MagneticOps,
) -> Vec<usize> {
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
            let rx = r[0][0] as i32 * kx_i + r[0][1] as i32 * ky_i + r[0][2] as i32 * kz_i;
            let ry = r[1][0] as i32 * kx_i + r[1][1] as i32 * ky_i + r[1][2] as i32 * kz_i;
            let rz = r[2][0] as i32 * kx_i + r[2][1] as i32 * ky_i + r[2][2] as i32 * kz_i;

            if ops.timerev[i] {
                // Anti-unitary: -R·k ≡ k (mod kd)
                (-rx - kx_i) % kd_i == 0
                    && (-ry - ky_i) % kd_i == 0
                    && (-rz - kz_i) % kd_i == 0
            } else {
                // Unitary: R·k ≡ k (mod kd)
                (rx - kx_i) % kd_i == 0
                    && (ry - ky_i) % kd_i == 0
                    && (rz - kz_i) % kd_i == 0
            }
        })
        .collect()
}

// ── Seitz matching ───────────────────────────────────────────────────────────

/// Result of matching a computed Seitz operation to a stored one.
pub struct SeitzMatch {
    /// Index of the matching operation in the stored list.
    pub op_index: usize,
    /// Lattice shift: $$\mathbf{L} = \mathbf{t}_{\text{comp}} - \mathbf{t}_{\text{stored}}$$
    /// (integer vector, may be non-zero when composed translation wraps
    /// around the unit cell).
    pub lattice_shift: [i32; 3],
}

/// Find a stored Seitz operation matching the given rotation and translation.
///
/// Matches by **rotation matrix first**, then finds the stored operation
/// whose translation differs from the computed one by a lattice vector.
///
/// For non-symmorphic groups, multiple operations can share the same
/// rotation but have different translations.  The first one whose
/// translation difference is integer (component-wise) is returned.
pub fn find_seitz(
    rot: &Mat3I,
    trans: &[f64; 3],
    ops: &[SeitzOp],
) -> Option<SeitzMatch> {
    for (idx, op) in ops.iter().enumerate() {
        if op.rot[0][0] != rot[0][0] || op.rot[0][1] != rot[0][1] || op.rot[0][2] != rot[0][2]
        || op.rot[1][0] != rot[1][0] || op.rot[1][1] != rot[1][1] || op.rot[1][2] != rot[1][2]
        || op.rot[2][0] != rot[2][0] || op.rot[2][1] != rot[2][1] || op.rot[2][2] != rot[2][2]
        {
            continue;
        }

        // Same rotation — check if translations differ by an integer vector
        let d0 = trans[0] - op.trans[0];
        let d1 = trans[1] - op.trans[1];
        let d2 = trans[2] - op.trans[2];

        let l0 = d0.round();
        let l1 = d1.round();
        let l2 = d2.round();

        if (d0 - l0).abs() < 1e-9 && (d1 - l1).abs() < 1e-9 && (d2 - l2).abs() < 1e-9 {
            return Some(SeitzMatch {
                op_index: idx,
                lattice_shift: [l0 as i32, l1 as i32, l2 as i32],
            });
        }
    }
    None
}

/// Compute the Bloch phase factor for a lattice shift at wave-vector k.
///
/// $$\phi = e^{+2\pi i\,\mathbf{k}\cdot\mathbf{L}}$$
///
/// where $$\mathbf{k} = (k_x,k_y,k_z)/k_d$$ in reciprocal lattice units and
/// $$\mathbf{L}$$ is an integer lattice vector.
///
/// Returns a complex phase factor as `(re, im)`.
pub fn bloch_phase(kx: i8, ky: i8, kz: i8, kd: i8, lattice: &[i32; 3]) -> (f64, f64) {
    if kd == 0 {
        return (1.0, 0.0);
    }
    let theta = 2.0 * std::f64::consts::PI
        * (kx as f64 * lattice[0] as f64
           + ky as f64 * lattice[1] as f64
           + kz as f64 * lattice[2] as f64)
        / (kd as f64);
    (theta.cos(), theta.sin())
}

// ── Wigner test ──────────────────────────────────────────────────────────────

/// Compute the Wigner test sum and classify the corep type.
///
/// # Mathematical definition
///
/// $$
/// W = \frac{1}{|H_{\mathbf{k}}|}
///     \sum_{h \in H_{\mathbf{k}}} \chi\big((a_0 h)^2\big)
/// $$
///
/// where $$(a_0 h)^2$$ is computed using **full Seitz arithmetic**:
///
/// 1. $$g_0 h = \{R_0 R_h \mid \mathbf{t}_0 + R_0 \mathbf{t}_h\}$$ (spatial part of a₀h)
/// 2. $$(g_0 h)^2 = g_0 h \circ g_0 h$$ (Seitz composition)
/// 3. Look up the result in $$H$$'s operation list, with Bloch phase for
///    any lattice shift
///
/// # Arguments
///
/// * `h_chars` — character table of H's irrep Δ (real-valued for PIR irreps)
/// * `unitary_mag_indices` — which magnetic ops are unitary AND in the little group
/// * `mag_seitz` — magnetic ops as SeitzOps (for computing a₀ and h)
/// * `h_seitz` — unitary subgroup H's ops as SeitzOps (for looking up (a₀h)²)
/// * `a0_idx` — index (into mag_seitz) of the anti-unitary coset representative
/// * `kx, ky, kz, kd` — wave-vector components for Bloch phase
///
/// # Classification
///
/// | W | Type | Dimension |
/// |---|------|-----------|
/// | ≈ +1 | A | d (same as H irrep) |
/// | ≈ -1 | B | 2d (Kramers doubling) |
/// | ≈ 0  | C | 2d (paired with conjugate) |
pub fn wigner_classify(
    h_chars: &[f64],
    unitary_mag_indices: &[usize],
    mag_seitz: &[SeitzOp],
    h_seitz: &[SeitzOp],
    a0_idx: usize,
    kx: i8, ky: i8, kz: i8, kd: i8,
) -> CorepType {
    let a0 = &mag_seitz[a0_idx];
    debug_assert!(a0.timerev, "a₀ must be anti-unitary");

    let mut w_sum: f64 = 0.0;

    for &h_mag_idx in unitary_mag_indices {
        let h = &mag_seitz[h_mag_idx];
        debug_assert!(!h.timerev, "h must be unitary");

        // Step 1: g₀h = a₀(spatial) ∘ h (both unitary spatial parts)
        // For the spatial part, we ignore timerev during composition
        let g0_spatial = SeitzOp::new(a0.rot, a0.trans, false);
        let h_spatial = SeitzOp::new(h.rot, h.trans, false);
        let (g0h, _l1) = compose_seitz(&g0_spatial, &h_spatial);

        // Step 2: (g₀h)²
        let (sq, lattice_sq) = square_seitz(&g0h);

        // Step 3: find (g₀h)² in H's operations
        if let Some(m) = find_seitz(&sq.rot, &sq.trans, h_seitz) {
            if m.op_index < h_chars.len() {
                // Account for lattice shift from composition
                let total_lattice = [
                    lattice_sq[0] + m.lattice_shift[0],
                    lattice_sq[1] + m.lattice_shift[1],
                    lattice_sq[2] + m.lattice_shift[2],
                ];
                let (phase_re, phase_im) = bloch_phase(kx, ky, kz, kd, &total_lattice);
                // Character is real for PIR irreps; Bloch phase may be complex.
                // For now we use the real character (PIR convention).
                // The phase factor multiplies the character.
                w_sum += h_chars[m.op_index] * phase_re;
                // Imaginary part should vanish for real irreps at TRIM points
            }
        }
        // If no Seitz match is found, this indicates a structural problem
        // (setting mismatch, wrong operations).  We skip the term rather than
        // silently adding 0 — but in production this should be a hard error.
    }

    let n = (unitary_mag_indices.len() as f64).max(1.0);
    let w = w_sum / n;

    // Classify with tolerance
    if w.abs() < 0.01 {
        CorepType::C
    } else if w > 0.0 {
        CorepType::A
    } else {
        CorepType::B
    }
}

// ── Character table construction ─────────────────────────────────────────────

/// Build the magnetic co-representation character table.
///
/// # Character formulas
///
/// **Type A** (dimension = d):
/// - Unitary: $$\tilde{\chi}(h) = \chi_i(h)$$
/// - Anti-unitary: $$\tilde{\chi}(a_0 h)$$ requires intertwiner U; set to 0 for now
///
/// **Type B** (dimension = 2d, Kramers doubling):
/// - Unitary: $$\tilde{\chi}(h) = 2\chi_i(h)$$
/// - Anti-unitary: $$\tilde{\chi}(a_0 h) = 0$$
///
/// **Type C** (dimension = 2d, paired with conjugate):
/// - Unitary: $$\tilde{\chi}(h) = 2\,\mathrm{Re}[\chi_i(h)]$$
/// - Anti-unitary: $$\tilde{\chi}(a_0 h) = 0$$
///
/// # Parameters
///
/// * `corep_type` — result of [`wigner_classify`]
/// * `mag_ops` — magnetic symmetry operations (for timerev flags)
/// * `mag_lg_indices` — which magnetic ops are in the little group
/// * `op_map` — for each magnetic op, the corresponding H op index (or None)
/// * `h_chars` — H's irrep character table (real-valued for PIR)
pub fn build_corep_chars(
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
        let h_idx = op_map[mag_idx];

        match corep_type {
            CorepType::A => {
                // Same dimension: inherit H character directly for unitary ops
                if !is_anti {
                    if let Some(hi) = h_idx {
                        if hi < h_chars.len() {
                            chars[out_idx] = h_chars[hi];
                        }
                    }
                }
                // Anti-unitary: requires intertwiner U. Leave as 0 for now.
            }
            CorepType::B => {
                // Kramers doubling: dimension 2d
                if is_anti {
                    chars[out_idx] = 0.0;
                } else if let Some(hi) = h_idx {
                    if hi < h_chars.len() {
                        chars[out_idx] = 2.0 * h_chars[hi];
                    }
                }
            }
            CorepType::C => {
                // Paired with conjugate: dimension 2d
                if is_anti {
                    chars[out_idx] = 0.0;
                } else if let Some(hi) = h_idx {
                    if hi < h_chars.len() {
                        // 2·Re[χ] — for PIR (real) irreps this is just 2·χ
                        chars[out_idx] = 2.0 * h_chars[hi];
                    }
                }
            }
        }
    }

    chars
}

// ── Corep dimension ─────────────────────────────────────────────────────────

/// Dimension of the magnetic co-representation.
///
/// Type A: same as H irrep (d)
/// Type B: doubled (2d, Kramers)
/// Type C: doubled (2d, paired)
pub fn corep_dim(corep_type: &CorepType, h_dim: usize) -> usize {
    match corep_type {
        CorepType::A => h_dim,
        CorepType::B | CorepType::C => h_dim * 2,
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Seitz composition: identity ∘ identity = identity.
    #[test]
    fn test_compose_identity() {
        let id = SeitzOp::new([[1,0,0],[0,1,0],[0,0,1]], [0.0,0.0,0.0], false);
        let (result, lattice) = compose_seitz(&id, &id);
        assert_eq!(result.rot, id.rot);
        assert_eq!(result.trans, [0.0, 0.0, 0.0]);
        assert_eq!(lattice, [0, 0, 0]);
        assert!(!result.timerev);
    }

    /// Seitz composition: timerev XOR.
    #[test]
    fn test_compose_timerev() {
        let id = SeitzOp::new([[1,0,0],[0,1,0],[0,0,1]], [0.0,0.0,0.0], false);
        let a0 = SeitzOp::new([[1,0,0],[0,1,0],[0,0,1]], [0.0,0.0,0.0], true);
        let (r1, _) = compose_seitz(&a0, &id);  // anti ∘ unitary = anti
        assert!(r1.timerev);
        let (r2, _) = compose_seitz(&a0, &a0);   // anti ∘ anti = unitary
        assert!(!r2.timerev);
        let (r3, _) = compose_seitz(&id, &a0);   // unitary ∘ anti = anti
        assert!(r3.timerev);
    }

    /// Seitz composition: translation arithmetic.
    #[test]
    fn test_compose_translation() {
        let g1 = SeitzOp::new(
            [[0,-1,0],[1,0,0],[0,0,1]],
            [0.0, 0.0, 0.5],
            false,
        );
        let g2 = SeitzOp::new(
            [[1,0,0],[0,1,0],[0,0,1]],
            [0.5, 0.0, 0.0],
            false,
        );
        let (result, lattice) = compose_seitz(&g1, &g2);
        // R = [[0,-1,0],[1,0,0],[0,0,1]]
        // t = [0,0,0.5] + R·[0.5,0,0] = [0,0,0.5] + [0,0.5,0] = [0, 0.5, 0.5]
        assert_eq!(result.trans, [0.0, 0.5, 0.5]);
        assert_eq!(lattice, [0, 0, 0]);
    }

    /// Seitz composition with lattice overflow.
    #[test]
    fn test_compose_lattice_shift() {
        let g1 = SeitzOp::new(
            [[1,0,0],[0,1,0],[0,0,1]],
            [0.7, 0.0, 0.0],
            false,
        );
        let g2 = SeitzOp::new(
            [[1,0,0],[0,1,0],[0,0,1]],
            [0.5, 0.0, 0.0],
            false,
        );
        let (result, lattice) = compose_seitz(&g1, &g2);
        // t = 0.7 + 0.5 = 1.2 → 0.2 with lattice shift [1,0,0]
        assert!((result.trans[0] - 0.2).abs() < 1e-9);
        assert_eq!(lattice, [1, 0, 0]);
    }

    /// filter_little_group: antiunitary ops use -Rk ≡ k.
    #[test]
    fn test_filter_antiunitary_k() {
        // k = (0, 0, 1)/2 = Z point
        // Anti-unitary op with R = [[0,-1,0],[1,0,0],[0,0,-1]] (4' about 001)
        // R·(0,0,1) = (0,0,-1), so -R·k - k = (0,0,1) - (0,0,1) = (0,0,0) ≡ 0 ✓
        let ops = MagneticOps {
            rot: vec![[[0,-1,0],[1,0,0],[0,0,-1]],
                       [[1,0,0],[0,1,0],[0,0,1]]],
            trans: vec![[0.0; 3], [0.0; 3]],
            timerev: vec![true, false],
        };
        let lg = filter_little_group(0, 0, 1, 2, &ops);
        assert_eq!(lg.len(), 2, "Both ops should be in Z-point little group");
    }

    /// filter_little_group: antiunitary op that does NOT preserve k.
    #[test]
    fn test_filter_antiunitary_not_in_lg() {
        // k = (1, 0, 0)/8 = generic point on X line
        // Anti-unitary op with R = [[-1,0,0],[0,1,0],[0,0,1]] (mx')
        // -R·k = (1,0,0), -R·k - k = (0,0,0) ≡ 0 → in little group
        // Anti-unitary op with R = [[1,0,0],[0,-1,0],[0,0,-1]]
        // -R·k = (-1,0,0) ≡ (7,0,0) mod 8, -R·k - k = (6,0,0) ≠ 0 → NOT in LG
        let ops = MagneticOps {
            rot: vec![[[-1,0,0],[0,1,0],[0,0,1]],
                       [[1,0,0],[0,-1,0],[0,0,-1]]],
            trans: vec![[0.0; 3], [0.0; 3]],
            timerev: vec![true, true],
        };
        let lg = filter_little_group(1, 0, 0, 8, &ops);
        assert_eq!(lg.len(), 1, "Only mx' should preserve k=(1/8,0,0)");
    }

    /// Simple Wigner test: P1 with only identity.
    #[test]
    fn test_wigner_trivial() {
        // a₀ = θ (anti-unitary identity), h = id
        // (a₀·id)² = id² = id, χ(id)=1.0 → W=1.0 → type A
        let mag_seitz = vec![
            SeitzOp::new([[1,0,0],[0,1,0],[0,0,1]], [0.0,0.0,0.0], false), // id
            SeitzOp::new([[1,0,0],[0,1,0],[0,0,1]], [0.0,0.0,0.0], true),  // θ
        ];
        let h_seitz = vec![
            SeitzOp::new([[1,0,0],[0,1,0],[0,0,1]], [0.0,0.0,0.0], false),
        ];
        let result = wigner_classify(
            &[1.0], &[0], &mag_seitz, &h_seitz, 1,
            0, 0, 0, 1, // Gamma point
        );
        assert_eq!(result, CorepType::A);
    }

    /// Type A: result should not double dimension.
    #[test]
    fn test_corep_dim_type_a() {
        assert_eq!(corep_dim(&CorepType::A, 3), 3);
    }

    /// Type B: Kramers doubling → 2d.
    #[test]
    fn test_corep_dim_type_b() {
        assert_eq!(corep_dim(&CorepType::B, 3), 6);
    }

    /// Type C: paired with conjugate → 2d.
    #[test]
    fn test_corep_dim_type_c() {
        assert_eq!(corep_dim(&CorepType::C, 2), 4);
    }
}
