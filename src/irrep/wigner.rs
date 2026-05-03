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

use num_complex::Complex64;
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

// ── Lattice arithmetic helpers ──────────────────────────────────────────────

/// Multiply a 3×3 integer matrix by a 3D integer vector.
#[inline]
pub fn mat_vec_i32(r: &Mat3I, v: &[i32; 3]) -> [i32; 3] {
    [
        r[0][0] as i32 * v[0] + r[0][1] as i32 * v[1] + r[0][2] as i32 * v[2],
        r[1][0] as i32 * v[0] + r[1][1] as i32 * v[1] + r[1][2] as i32 * v[2],
        r[2][0] as i32 * v[0] + r[2][1] as i32 * v[1] + r[2][2] as i32 * v[2],
    ]
}

/// Add two [i32; 3] vectors.
#[inline]
pub fn add3(a: &[i32; 3], b: &[i32; 3]) -> [i32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
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
pub fn bloch_phase(kx: i8, ky: i8, kz: i8, kd: i8, lattice: &[i32; 3]) -> Complex64 {
    if kd == 0 {
        return Complex64::new(1.0, 0.0);
    }
    let theta = 2.0 * std::f64::consts::PI
        * (kx as f64 * lattice[0] as f64
           + ky as f64 * lattice[1] as f64
           + kz as f64 * lattice[2] as f64)
        / (kd as f64);
    Complex64::new(theta.cos(), theta.sin())
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

        let g0_spatial = SeitzOp::new(a0.rot, a0.trans, false);
        let h_spatial = SeitzOp::new(h.rot, h.trans, false);
        let (g0h, l1) = compose_seitz(&g0_spatial, &h_spatial);
        let (sq, lattice_sq) = square_seitz(&g0h);

        if let Some(m) = find_seitz(&sq.rot, &sq.trans, h_seitz) {
            if m.op_index < h_chars.len() {
                // Total lattice shift = L_sq + L_match + L1 + R_{g0h}·L1
                let r_l1 = mat_vec_i32(&g0h.rot, &l1);
                let total_lattice = add3(
                    &add3(&lattice_sq, &m.lattice_shift),
                    &add3(&l1, &r_l1),
                );
                let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
                let contrib = h_chars[m.op_index] * phase.re;
                w_sum += contrib;
                eprintln!("    wigner: h[{}]→H[{}] sq=H[{}] L={:?} ph={:.2} χ={:.2} → {:.2}",
                    h_mag_idx, "?", m.op_index, total_lattice, phase.re,
                    h_chars[m.op_index], contrib);
            }
        } else {
            eprintln!("    wigner: h[{}] sq R=[{},{},{};...] t=({:.3},{:.3},{:.3}) NOT FOUND",
                h_mag_idx, sq.rot[0][0],sq.rot[0][1],sq.rot[0][2],
                sq.trans[0],sq.trans[1],sq.trans[2]);
        }
    }

    let n = (unitary_mag_indices.len() as f64).max(1.0);
    let w = w_sum / n;

    // Classify with tolerance
    eprintln!("DEBUG wigner_classify: w_sum={:.4} n_unitary={} W={:.4} k=({},{},{})/{}",
        w_sum, unitary_mag_indices.len(), w, kx, ky, kz, kd);
    if w.abs() < 0.01 {
        CorepType::C
    } else if w > 0.0 {
        CorepType::A
    } else {
        CorepType::B
    }
}

// ── CIR-based Wigner test ───────────────────────────────────────────────────

/// Wigner test using CIR (complex) character tables.
///
/// For compound PIR irreps like Z1Z4 = Z1 ⊕ Z4, the underlying CIR
/// irreps Z1, Z4 are complex and may individually give Type C under
/// the antiunitary operation, even though the combined PIR gives Type A.
///
/// This function is called once per CIR component.  The caller (in
/// `compute_corepresentation`) should:
/// 1. Check `irrep.cir_component_count() > 0`
/// 2. Loop over components, call this function for each
/// 3. If any component gives Type C, the overall corep is Type C
///
/// $$ W = \frac{1}{|H|} \sum_{h \in H} \chi_{\text{CIR}}((a_0 h)^2) $$
///
/// where $$\chi_{\text{CIR}}$$ is complex-valued.  W is complex; we
/// classify by $$|W| < 0.01$$ → Type C, Re(W) > 0 → Type A, else Type B.
pub fn wigner_classify_cir(
    cir_chars: &[f64],  // (re, im) pairs for one CIR component
    unitary_mag_indices: &[usize],
    mag_seitz: &[SeitzOp],
    h_seitz: &[SeitzOp],
    a0_idx: usize,
    kx: i8, ky: i8, kz: i8, kd: i8,
) -> CorepType {
    let a0 = &mag_seitz[a0_idx];
    let mut w_sum = Complex64::new(0.0, 0.0);
    let mut n_plus = 0u32;
    let mut n_minus = 0u32;

    for &h_mag_idx in unitary_mag_indices {
        let h = &mag_seitz[h_mag_idx];
        let g0_spatial = SeitzOp::new(a0.rot, a0.trans, false);
        let h_spatial = SeitzOp::new(h.rot, h.trans, false);
        let (g0h, l1) = compose_seitz(&g0_spatial, &h_spatial);
        let (sq, lattice_sq) = square_seitz(&g0h);

        if let Some(m) = find_seitz(&sq.rot, &sq.trans, h_seitz) {
            let r_l1 = mat_vec_i32(&g0h.rot, &l1);
            let total_lattice = add3(
                &add3(&lattice_sq, &m.lattice_shift),
                &add3(&l1, &r_l1),
            );
            let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
            let chi = cir_char_at(cir_chars, m.op_index);
            w_sum += chi * phase;
            // Phase parity stats
            if phase.re > 0.5 { n_plus += 1; }
            else if phase.re < -0.5 { n_minus += 1; }
            eprintln!("    cir: h[{}]→H[{}] Lz_par={} ph={:.2} χ={:.2} → {:.2}",
                h_mag_idx, m.op_index,
                ((total_lattice[2] % 2) + 2) % 2,
                phase, chi, chi * phase);
        }
    }

    eprintln!("    phase stats: +={} -={}", n_plus, n_minus);
    let n = (unitary_mag_indices.len() as f64).max(1.0);
    let w = w_sum / n;
    eprintln!("DEBUG wigner_classify_cir: W=({:.8},{:.8}) |W|={:.4} k=({},{},{})/{}",
        w.re, w.im, w.norm(), kx, ky, kz, kd);

    let tol = 1e-6;
    if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        CorepType::A
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        CorepType::B
    } else if w.norm() < tol {
        CorepType::C
    } else {
        panic!(
            "Non-quantized Wigner indicator W=({:.8},{:.8}); expected 0, +1, or -1. \
             Check Seitz matching, character table ordering, and a₀ coset coverage.",
            w.re, w.im
        );
    }
}

/// Verify that `h_seitz` operation ordering matches the CIR character table.
/// Prints all operations with their characters for manual inspection.
#[cfg(test)]
pub fn debug_char_order(cir_chars: &[f64], h_seitz: &[SeitzOp], label: &str) {
    eprintln!("=== Character order check: {} ===", label);
    for (i, op) in h_seitz.iter().enumerate() {
        let re = cir_chars.get(2 * i).copied().unwrap_or(999.0);
        let im = cir_chars.get(2 * i + 1).copied().unwrap_or(999.0);
        let is_id = op.rot[0][0] == 1 && op.rot[0][1] == 0 && op.rot[0][2] == 0
                 && op.rot[1][0] == 0 && op.rot[1][1] == 1 && op.rot[1][2] == 0
                 && op.rot[2][0] == 0 && op.rot[2][1] == 0 && op.rot[2][2] == 1
                 && op.trans[0].abs() < 0.01 && op.trans[1].abs() < 0.01 && op.trans[2].abs() < 0.01;
        eprintln!("  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) chi=({:.3},{:.3}){}",
            i,
            op.rot[0][0],op.rot[0][1],op.rot[0][2],
            op.rot[1][0],op.rot[1][1],op.rot[1][2],
            op.rot[2][0],op.rot[2][1],op.rot[2][2],
            op.trans[0], op.trans[1], op.trans[2],
            re, im,
            if is_id { " ← ID" } else { "" },
        );
    }
}

/// Diagnostic: unwrapped Seitz square without intermediate normalization.
/// Computes (g₀h)² directly from raw translations, then compares
/// with the normalized+matched result.  Used to debug phase parity.
pub fn debug_unwrapped_square(
    h_mag_idx: usize,
    a0_idx: usize,
    mag_seitz: &[SeitzOp],
    h_seitz: &[SeitzOp],
    kx: i8, ky: i8, kz: i8, kd: i8,
) {
    let a0 = &mag_seitz[a0_idx];
    let h = &mag_seitz[h_mag_idx];

    // Step 1: g₀h raw (no normalization)
    let rc = mat_multiply_matrix_i3(&a0.rot, &h.rot);
    let r0_th = mat_vec_f64(&a0.rot, &h.trans);
    let tc_raw = [a0.trans[0] + r0_th[0], a0.trans[1] + r0_th[1], a0.trans[2] + r0_th[2]];

    // Step 2: (g₀h)² raw
    let rsq = mat_multiply_matrix_i3(&rc, &rc);
    let rc_tc = mat_vec_f64(&rc, &tc_raw);
    let tsq_raw = [tc_raw[0] + rc_tc[0], tc_raw[1] + rc_tc[1], tc_raw[2] + rc_tc[2]];

    eprintln!("=== unwrapped square: h[{}] ===", h_mag_idx);
    eprintln!("  a0: R={:?}, t={:?}", a0.rot, a0.trans);
    eprintln!("  h : R={:?}, t={:?}", h.rot, h.trans);
    eprintln!("  g0h raw: R={:?}, t={:?}", rc, tc_raw);
    eprintln!("  sq raw : R={:?}, t={:?}", rsq, tsq_raw);

    // Normalize for matching
    let (tsq_mod, l_reduce) = reduce01_with_lattice(&tsq_raw);
    eprintln!("  sq mod : t={:?}, L_reduce={:?}", tsq_mod, l_reduce);

    if let Some(m) = find_seitz(&rsq, &tsq_mod, h_seitz) {
        let stored_t = &h_seitz[m.op_index].trans;
        // Direct lattice difference: tsq_raw - stored_t
        let l_direct = [
            (tsq_raw[0] - stored_t[0]).round() as i32,
            (tsq_raw[1] - stored_t[1]).round() as i32,
            (tsq_raw[2] - stored_t[2]).round() as i32,
        ];
        let lz_par = ((l_direct[2] % 2) + 2) % 2;
        let phase = bloch_phase(kx, ky, kz, kd, &l_direct);

        eprintln!("  matched H[{}]: t_stored={:?}", m.op_index, stored_t);
        eprintln!("  L_direct={:?} Lz_par={} phase={:.2}", l_direct, lz_par, phase);
        eprintln!("  m.lattice_shift={:?} (from normalized match)", m.lattice_shift);
    } else {
        eprintln!("  NOT FOUND in h_seitz");
    }
}

/// Diagnostic: direct anti-coset Wigner sum.
/// Uses ALL antiunitary little-group ops b directly (not a₀h construction).
/// If this gives different phase parity than wigner_classify_cir,
/// the a₀h construction is wrong.
pub fn wigner_direct_anti_coset(
    cir_chars: &[f64],
    anti_lg_indices: &[usize],
    mag_seitz: &[SeitzOp],
    h_seitz: &[SeitzOp],
    kx: i8, ky: i8, kz: i8, kd: i8,
) -> Complex64 {
    let mut sum = Complex64::ZERO;
    let mut n_plus = 0u32;
    let mut n_minus = 0u32;

    for &b_idx in anti_lg_indices {
        let b = &mag_seitz[b_idx];
        let (sq, lattice_sq) = square_seitz(b);
        let m = find_seitz(&sq.rot, &sq.trans, h_seitz)
            .unwrap_or_else(|| panic!("direct anti: b[{}]^2 not found in H", b_idx));

        let total_lattice = add3(&lattice_sq, &m.lattice_shift);
        let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
        let chi = cir_char_at(cir_chars, m.op_index);
        let contrib = chi * phase;
        sum += contrib;

        if phase.re > 0.5 { n_plus += 1; }
        else if phase.re < -0.5 { n_minus += 1; }

        eprintln!("  direct: b[{}]^2→H[{}] L={:?} ph={:.2} χ={:.2} → {:.2}",
            b_idx, m.op_index, total_lattice, phase, chi, contrib);
    }
    let w = sum / (anti_lg_indices.len() as f64);
    eprintln!("  direct anti stats: +={} -={} W={:.4}", n_plus, n_minus, w);
    w
}

// ── Internal helpers ────────────────────────────────────────────────────────

/// f64 vector: R · v
fn mat_vec_f64(r: &Mat3I, v: &[f64; 3]) -> [f64; 3] {
    [
        r[0][0] as f64 * v[0] + r[0][1] as f64 * v[1] + r[0][2] as f64 * v[2],
        r[1][0] as f64 * v[0] + r[1][1] as f64 * v[1] + r[1][2] as f64 * v[2],
        r[2][0] as f64 * v[0] + r[2][1] as f64 * v[1] + r[2][2] as f64 * v[2],
    ]
}

/// Normalize translation to [0,1) and return discarded integer shift.
fn reduce01_with_lattice(t: &[f64; 3]) -> ([f64; 3], [i32; 3]) {
    let mut tr = [0.0f64; 3];
    let mut l = [0i32; 3];
    for i in 0..3 {
        let fl = t[i].floor();
        l[i] = fl as i32;
        tr[i] = t[i] - fl;
        if tr[i] < 0.0 { tr[i] += 1.0; l[i] -= 1; }
    }
    (tr, l)
}

/// Build index map from `h_seitz` (spglib H_ops order) to CIR operation order.
///
/// Returns `None` if any H operation cannot be matched (e.g., missing rotation data).
pub fn build_h_to_cir_map(h_seitz: &[SeitzOp], cir_rots: &[i32]) -> Option<Vec<usize>> {
    let n_ops = h_seitz.len();
    let n_cir_ops = cir_rots.len() / 9;
    if n_cir_ops == 0 {
        return None;
    }
    let mut map = vec![0usize; n_ops];

    eprintln!("build_h_to_cir_map: n_ops={} n_cir_ops={}", n_ops, n_cir_ops);
    for h_idx in 0..n_ops {
        let h_op = &h_seitz[h_idx];
        let r = &h_op.rot;
        let found = (0..n_cir_ops).find(|&c| {
            let off = c * 9;
            off + 8 < cir_rots.len()
            && cir_rots[off] == r[0][0] && cir_rots[off+1] == r[0][1] && cir_rots[off+2] == r[0][2]
            && cir_rots[off+3] == r[1][0] && cir_rots[off+4] == r[1][1] && cir_rots[off+5] == r[1][2]
            && cir_rots[off+6] == r[2][0] && cir_rots[off+7] == r[2][1] && cir_rots[off+8] == r[2][2]
        });
        match found {
            Some(c) => map[h_idx] = c,
            None => {
                eprintln!("build_h_to_cir_map: H[{}] R=[{},{},{};{},{},{};{},{},{}] not found in rots ({} ops)",
                    h_idx, r[0][0],r[0][1],r[0][2], r[1][0],r[1][1],r[1][2], r[2][0],r[2][1],r[2][2],
                    n_cir_ops);
                return None;
            }
        }
    }
    Some(map)
}

/// Reorder CIR complex characters from ISOTROPY order to H_ops (spglib) order.
///
/// After reordering, `out[2*h_idx]` and `out[2*h_idx+1]` give the (re, im)
/// character for the operation at `h_seitz[h_idx]`.
pub fn reorder_cir_chars(cir_chars: &[f64], h_to_cir: &[usize]) -> Vec<f64> {
    let n_ops = h_to_cir.len();
    let mut reordered = vec![0.0f64; n_ops * 2];
    for h_idx in 0..n_ops {
        let c_idx = h_to_cir[h_idx];
        if 2 * c_idx + 1 < cir_chars.len() {
            reordered[2 * h_idx] = cir_chars[2 * c_idx];
            reordered[2 * h_idx + 1] = cir_chars[2 * c_idx + 1];
        }
    }
    reordered
}

/// Helper: read a complex character from (re, im) pair array.
#[inline]
fn cir_char_at(cir_chars: &[f64], op_idx: usize) -> Complex64 {
    let i = 2 * op_idx;
    if i + 1 < cir_chars.len() {
        Complex64::new(cir_chars[i], cir_chars[i + 1])
    } else {
        Complex64::ZERO
    }
}

// ── Conjugate representation & partner finding ──────────────────────────────

/// Compute the conjugate character of Δ under anti-unitary operation a₀.
///
/// The conjugate representation is defined as:
///
/// $$ \chi^{a_0}(h) = \chi\big(a_0^{-1} h a_0\big)^* $$
///
/// where $$a_0^{-1} h a_0$$ is computed via Seitz composition: first conjugate
/// h by a₀'s spatial inverse, then apply to a₀'s spatial part.
///
/// For each unitary op h in H, we compute $$h' = g_0^{-1} \circ h \circ g_0$$
/// (where g₀ is the spatial part of a₀), find h' in H's operation list,
/// and read off $$\chi(h')^*$$.
///
/// # Arguments
/// * `h_chars` — character table of irrep Δ_i
/// * `h_seitz` — H's operations as SeitzOps
/// * `a0` — the anti-unitary coset representative a₀ = θ·g₀
/// * `kx, ky, kz, kd` — wave-vector for Bloch phases
///
/// # Returns
/// `(conj_chars, h_to_conj_map)` where `conj_chars[h_idx] = χ(a₀⁻¹ h a₀)*`
/// and `h_to_conj_map[h_idx]` is the H-index of the conjugated operation.
pub fn conjugate_chars(
    h_chars: &[f64],
    h_seitz: &[SeitzOp],
    a0: &SeitzOp,
    kx: i8, ky: i8, kz: i8, kd: i8,
) -> (Vec<f64>, Vec<Option<usize>>) {
    let n = h_seitz.len();
    let mut conj = vec![0.0f64; n];
    let mut h_to_conj = vec![None; n];

    // g₀⁻¹: inverse of the spatial part of a₀
    // For orthogonal matrices, R⁻¹ = Rᵀ
    let g0_inv_rot = {
        let r = &a0.rot;
        [
            [r[0][0], r[1][0], r[2][0]],
            [r[0][1], r[1][1], r[2][1]],
            [r[0][2], r[1][2], r[2][2]],
        ]
    };
    // g₀⁻¹ translation: t_inv = -R⁻¹·t
    let mut g0_inv_trans = [0.0f64; 3];
    for i in 0..3 {
        let s = g0_inv_rot[i][0] as f64 * a0.trans[0]
              + g0_inv_rot[i][1] as f64 * a0.trans[1]
              + g0_inv_rot[i][2] as f64 * a0.trans[2];
        g0_inv_trans[i] = -s;
    }
    let g0_inv = SeitzOp::new(g0_inv_rot, g0_inv_trans, false);

    for h_idx in 0..n {
        let h = &h_seitz[h_idx];

        // Compute h' = g₀⁻¹ ∘ h ∘ g₀
        let (hg0, _) = compose_seitz(h, a0);  // h ∘ g₀ (treat a0 as spatial only)
        let (h_prime, lattice) = compose_seitz(&g0_inv, &hg0); // g₀⁻¹ ∘ (h ∘ g₀)

        // Find h' in H's operation list
        if let Some(m) = find_seitz(&h_prime.rot, &h_prime.trans, h_seitz) {
            let total_lattice = [
                lattice[0] + m.lattice_shift[0],
                lattice[1] + m.lattice_shift[1],
                lattice[2] + m.lattice_shift[2],
            ];
            let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
            h_to_conj[h_idx] = Some(m.op_index);
            if m.op_index < h_chars.len() {
                // χ(a₀⁻¹ h a₀)* : complex conjugate of (χ · phase)
                // For real PIR characters: (χ * phase).conj() = χ * phase.conj()
                let val = Complex64::new(h_chars[m.op_index], 0.0) * phase.conj();
                conj[h_idx] = val.re; // PIR chars are real
            }
        }
    }

    (conj, h_to_conj)
}

/// Find the partner irrep for Type C by comparing conjugate characters.
///
/// For each irrep Δ_j of H at the same k-point, compute the overlap
///
/// $$ \text{overlap}_{ij} = \frac{1}{|H|}
///     \sum_{h \in H} \chi_i^{a_0}(h) \cdot \chi_j(h)^* $$
///
/// If Δ_j is equivalent to Δ_i^{a₀}, the overlap ≈ 1 (or the dimension d).
/// For Type C, the partner is the irrep with the highest overlap.
///
/// Returns the index of the partner irrep in `candidates`, or `None` if
/// no clear partner is found (Type A/B case).
pub fn find_partner(
    conj_chars: &[f64],
    candidate_chars: &[&[f64]],  // character tables of candidate irreps
) -> Option<usize> {
    let n_ops = conj_chars.len();
    if n_ops == 0 || candidate_chars.is_empty() {
        return None;
    }

    let norm = 1.0 / (n_ops as f64);

    // For each candidate irrep, compute overlap with conjugate
    let mut best_idx = None;
    let mut best_overlap = 0.0f64;

    for (j, chars_j) in candidate_chars.iter().enumerate() {
        if chars_j.len() < n_ops { continue; }
        let overlap: f64 = (0..n_ops)
            .map(|h| conj_chars[h] * chars_j[h])
            .sum();
        let overlap_norm = overlap.abs() * norm;

        if overlap_norm > best_overlap {
            best_overlap = overlap_norm;
            best_idx = Some(j);
        }
    }

    // Return partner if overlap is significantly above noise
    // For n-dimensional irrep, overlap ≈ d²/|H| per operation
    if best_overlap > 0.1 {
        best_idx
    } else {
        None
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
    partner_chars: Option<&[f64]>,  // for Type C: character table of paired irrep
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
                // Paired with conjugate irrep: dimension 2d
                if is_anti {
                    chars[out_idx] = 0.0;
                } else if let Some(hi) = h_idx {
                    let chi_i = if hi < h_chars.len() { h_chars[hi] } else { 0.0 };
                    let chi_partner = if let Some(pc) = partner_chars {
                        if hi < pc.len() { pc[hi] } else { 0.0 }
                    } else {
                        chi_i
                    };
                    chars[out_idx] = chi_i + chi_partner;
                }
            }
            CorepType::Unsupported => {
                chars[out_idx] = f64::NAN;
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
        CorepType::Unsupported => 0,
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
