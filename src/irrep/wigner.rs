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

use super::corep::CorepType;
use super::types::KVector;
use crate::SymmetryOps;
use crate::mathfunc::{
    Mat3, Mat3I, mat_get_determinant_i3, mat_inverse_matrix_d3, mat_multiply_matrix_d3,
    mat_multiply_matrix_i3,
};
use num_complex::Complex64;

/// Error returned when Wigner's test cannot classify a co-representation
/// as A, B, or C.
#[derive(Debug, Clone, PartialEq)]
pub struct WignerClassificationError {
    pub reason: String,
    pub wigner_value: Option<f64>,
}

impl WignerClassificationError {
    pub fn new(reason: impl Into<String>) -> Self {
        Self {
            reason: reason.into(),
            wigner_value: None,
        }
    }

    pub fn with_value(reason: impl Into<String>, w: f64) -> Self {
        Self {
            reason: reason.into(),
            wigner_value: Some(w),
        }
    }
}

impl std::fmt::Display for WignerClassificationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(w) = self.wigner_value {
            write!(f, "{} (W = {:.6})", self.reason, w)
        } else {
            write!(f, "{}", self.reason)
        }
    }
}

// ── Diagnostic counters for SU(2) central-element relation ──────────────────

use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

// ── Diagnostic counters for setting-transform origin solving ─────────────────

/// find_setting_transform was called
static XF_CALLED: AtomicUsize = AtomicUsize::new(0);
/// At least one valid (T,s) was found
static XF_FOUND: AtomicUsize = AtomicUsize::new(0);
/// Identity basis (T=I) worked
static XF_IDENTITY: AtomicUsize = AtomicUsize::new(0);
/// Non-identity basis (T≠I) worked
static XF_NON_IDENTITY: AtomicUsize = AtomicUsize::new(0);
/// Origin s is non-zero (|s| > 1e-8)
static XF_NONZERO_ORIGIN: AtomicUsize = AtomicUsize::new(0);
/// Multiple valid transforms found (ambiguous)
static XF_AMBIGUOUS: AtomicUsize = AtomicUsize::new(0);

pub fn read_xf_counters() -> (usize, usize, usize, usize, usize, usize) {
    (
        XF_CALLED.load(Ordering::Relaxed),
        XF_FOUND.load(Ordering::Relaxed),
        XF_IDENTITY.load(Ordering::Relaxed),
        XF_NON_IDENTITY.load(Ordering::Relaxed),
        XF_NONZERO_ORIGIN.load(Ordering::Relaxed),
        XF_AMBIGUOUS.load(Ordering::Relaxed),
    )
}

/// u_sq ≈ u_k  (same lift, no central element)
static SU2_REL_SAME: AtomicUsize = AtomicUsize::new(0);
/// u_sq ≈ -u_k (differs by Ebar = [-1,0,0,0])
static SU2_REL_EBAR: AtomicUsize = AtomicUsize::new(0);
/// u_sq not related to ±u_k (should never happen)
static SU2_REL_NONE: AtomicUsize = AtomicUsize::new(0);
static NONE_MATCH_OTHER_LG: AtomicUsize = AtomicUsize::new(0);
static NONE_MATCH_OTHER_GLOBAL: AtomicUsize = AtomicUsize::new(0);
static NONE_NO_MATCH_HAS_CAND: AtomicUsize = AtomicUsize::new(0);
static NONE_NO_CANDIDATE: AtomicUsize = AtomicUsize::new(0);
// Alternative antiunitary square formulas tested on NONE
static NONE_ALT_RAW: AtomicUsize = AtomicUsize::new(0);
static NONE_ALT_NEG_RAW: AtomicUsize = AtomicUsize::new(0);
static NONE_ALT_UUSTAR: AtomicUsize = AtomicUsize::new(0);
static NONE_ALT_NEG_UUSTAR: AtomicUsize = AtomicUsize::new(0);
static NONE_ALT_STARU: AtomicUsize = AtomicUsize::new(0);
static NONE_ALT_NEG_STARU: AtomicUsize = AtomicUsize::new(0);
static NONE_ALT_NONE: AtomicUsize = AtomicUsize::new(0);
// J-insertion antiunitary square: J = i*sigma_y = [0,0,1,0]
static NONE_JU_JU_STAR: AtomicUsize = AtomicUsize::new(0);
static NONE_NEG_JU_JU_STAR: AtomicUsize = AtomicUsize::new(0);
static NONE_UJ_UJ_STAR: AtomicUsize = AtomicUsize::new(0);
static NONE_NEG_UJ_UJ_STAR: AtomicUsize = AtomicUsize::new(0);
static NONE_J_NONE: AtomicUsize = AtomicUsize::new(0);

// MSG-gauge path triage (added 2026-06-18)
pub static MSG_GAUGE_OK: AtomicUsize = AtomicUsize::new(0);
pub static MSG_GAUGE_MAP_FAIL: AtomicUsize = AtomicUsize::new(0);
pub static MSG_GAUGE_W_FAIL: AtomicUsize = AtomicUsize::new(0);
pub static OLD_PATH_OK: AtomicUsize = AtomicUsize::new(0);
pub static OLD_PATH_FAIL: AtomicUsize = AtomicUsize::new(0);

// build_h_to_spin_map failure classification (added 2026-06-18)
pub static H2S_OK: AtomicUsize = AtomicUsize::new(0);
pub static H2S_AMBIGUOUS: AtomicUsize = AtomicUsize::new(0);
pub static H2S_MISSING: AtomicUsize = AtomicUsize::new(0);

// det distribution for NONE
static NONE_DET_A0_P1: AtomicUsize = AtomicUsize::new(0);
static NONE_DET_A0_M1: AtomicUsize = AtomicUsize::new(0);
static NONE_DET_G0H_P1: AtomicUsize = AtomicUsize::new(0);
static NONE_DET_G0H_M1: AtomicUsize = AtomicUsize::new(0);
// G-gauge oracle: compute central relation entirely in G spin database
static GGAUGE_SAME: AtomicUsize = AtomicUsize::new(0);
static GGAUGE_EBAR: AtomicUsize = AtomicUsize::new(0);
static GGAUGE_NONE: AtomicUsize = AtomicUsize::new(0);
static GGAUGE_H_LOOKUP_FAIL: AtomicUsize = AtomicUsize::new(0);
static GGAUGE_SQ_LOOKUP_FAIL: AtomicUsize = AtomicUsize::new(0);

/// Reset the SU(2) relation counters.
pub fn reset_su2_rel_counters() {
    SU2_REL_SAME.store(0, Ordering::Relaxed);
    SU2_REL_EBAR.store(0, Ordering::Relaxed);
    SU2_REL_NONE.store(0, Ordering::Relaxed);
    NONE_MATCH_OTHER_LG.store(0, Ordering::Relaxed);
    NONE_MATCH_OTHER_GLOBAL.store(0, Ordering::Relaxed);
    NONE_NO_MATCH_HAS_CAND.store(0, Ordering::Relaxed);
    NONE_NO_CANDIDATE.store(0, Ordering::Relaxed);
    NONE_ALT_RAW.store(0, Ordering::Relaxed);
    NONE_ALT_NEG_RAW.store(0, Ordering::Relaxed);
    NONE_ALT_UUSTAR.store(0, Ordering::Relaxed);
    NONE_ALT_NEG_UUSTAR.store(0, Ordering::Relaxed);
    NONE_ALT_STARU.store(0, Ordering::Relaxed);
    NONE_ALT_NEG_STARU.store(0, Ordering::Relaxed);
    NONE_ALT_NONE.store(0, Ordering::Relaxed);
    NONE_DET_A0_P1.store(0, Ordering::Relaxed);
    NONE_DET_A0_M1.store(0, Ordering::Relaxed);
    NONE_DET_G0H_P1.store(0, Ordering::Relaxed);
    NONE_DET_G0H_M1.store(0, Ordering::Relaxed);
    GGAUGE_SAME.store(0, Ordering::Relaxed);
    GGAUGE_EBAR.store(0, Ordering::Relaxed);
    GGAUGE_NONE.store(0, Ordering::Relaxed);
    GGAUGE_H_LOOKUP_FAIL.store(0, Ordering::Relaxed);
    GGAUGE_SQ_LOOKUP_FAIL.store(0, Ordering::Relaxed);
    NONE_JU_JU_STAR.store(0, Ordering::Relaxed);
    NONE_NEG_JU_JU_STAR.store(0, Ordering::Relaxed);
    NONE_UJ_UJ_STAR.store(0, Ordering::Relaxed);
    NONE_NEG_UJ_UJ_STAR.store(0, Ordering::Relaxed);
    NONE_J_NONE.store(0, Ordering::Relaxed);
}

/// Read the SU(2) relation counters: (same, ebar, none).
pub fn read_su2_rel_counters() -> (usize, usize, usize) {
    (
        SU2_REL_SAME.load(Ordering::Relaxed),
        SU2_REL_EBAR.load(Ordering::Relaxed),
        SU2_REL_NONE.load(Ordering::Relaxed),
    )
}

/// Read the NONE sub-category counters.
pub fn read_none_counters() -> (usize, usize, usize, usize) {
    (
        NONE_MATCH_OTHER_LG.load(Ordering::Relaxed),
        NONE_MATCH_OTHER_GLOBAL.load(Ordering::Relaxed),
        NONE_NO_MATCH_HAS_CAND.load(Ordering::Relaxed),
        NONE_NO_CANDIDATE.load(Ordering::Relaxed),
    )
}

/// Read NONE alternative formula counters.
pub fn read_none_alt_counters() -> (usize, usize, usize, usize, usize, usize, usize) {
    (
        NONE_ALT_RAW.load(Ordering::Relaxed),
        NONE_ALT_NEG_RAW.load(Ordering::Relaxed),
        NONE_ALT_UUSTAR.load(Ordering::Relaxed),
        NONE_ALT_NEG_UUSTAR.load(Ordering::Relaxed),
        NONE_ALT_STARU.load(Ordering::Relaxed),
        NONE_ALT_NEG_STARU.load(Ordering::Relaxed),
        NONE_ALT_NONE.load(Ordering::Relaxed),
    )
}

/// Read NONE det distribution.
pub fn read_none_det_counters() -> (usize, usize, usize, usize) {
    (
        NONE_DET_A0_P1.load(Ordering::Relaxed),
        NONE_DET_A0_M1.load(Ordering::Relaxed),
        NONE_DET_G0H_P1.load(Ordering::Relaxed),
        NONE_DET_G0H_M1.load(Ordering::Relaxed),
    )
}

/// Read G-gauge oracle counters.
pub fn read_ggauge_counters() -> (usize, usize, usize, usize, usize) {
    (
        GGAUGE_SAME.load(Ordering::Relaxed),
        GGAUGE_EBAR.load(Ordering::Relaxed),
        GGAUGE_NONE.load(Ordering::Relaxed),
        GGAUGE_H_LOOKUP_FAIL.load(Ordering::Relaxed),
        GGAUGE_SQ_LOOKUP_FAIL.load(Ordering::Relaxed),
    )
}

/// Read J-insertion oracle counters.
pub fn read_j_oracle_counters() -> (usize, usize, usize, usize, usize) {
    (
        NONE_JU_JU_STAR.load(Ordering::Relaxed),
        NONE_NEG_JU_JU_STAR.load(Ordering::Relaxed),
        NONE_UJ_UJ_STAR.load(Ordering::Relaxed),
        NONE_NEG_UJ_UJ_STAR.load(Ordering::Relaxed),
        NONE_J_NONE.load(Ordering::Relaxed),
    )
}

/// G→H spin frame parity for signed-permutation setting transforms.
///
/// When the spatial coordinate transform has det(P) = -1, spin (axial vector)
/// transforms under Q = det(P)·P ≠ P.  For a rotation R in G frame:
///
/// ```text
/// R' = P·R·P⁻¹          (spatial transform to H frame)
/// U_G→H = U_Q · U_G(R) · U_Q*   (spin lift in H frame)
/// ε(R) = sign(U_H(R') · U_G→H)  (G→H parity, ±1)
/// ```
///
/// This function computes ε(R) for a single rotation.  The Q matrix must
/// be a signed-permutation with det=1 (so it represents a proper rotation).
///
/// Returns `None` if either spin table lacks the rotation.
/// SU(2) lift of a signed-permutation matrix Q (det=1, entries ∈ {0,±1}).
/// Returns the unit quaternion [u0, u1, u2, u3] representing this rotation.
fn signed_perm_to_quat(q: &[[i32; 3]; 3]) -> Option<[f64; 4]> {
    let tr = q[0][0] + q[1][1] + q[2][2];
    let cos_half = ((tr + 1) as f64 / 4.0).sqrt();
    let sin_half = (1.0 - cos_half * cos_half).sqrt();
    let axis: [f64; 3] = if tr == -1 {
        // 180°: use first non-zero column of Q+I
        let mut ax = [0.0f64; 3];
        for (j, _) in q[0].iter().enumerate() {
            let v = if j == 0 {
                [q[0][j] as f64 + 1.0, q[1][j] as f64, q[2][j] as f64]
            } else if j == 1 {
                [q[0][j] as f64, q[1][j] as f64 + 1.0, q[2][j] as f64]
            } else {
                [q[0][j] as f64, q[1][j] as f64, q[2][j] as f64 + 1.0]
            };
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if n > 0.5 {
                ax = [v[0] / n, v[1] / n, v[2] / n];
                break;
            }
        }
        let n = (ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]).sqrt();
        if n < 0.5 {
            return None;
        }
        [ax[0] / n, ax[1] / n, ax[2] / n]
    } else if tr == 3 {
        [0.0, 0.0, 1.0]
    } else {
        let ax_x = (q[2][1] - q[1][2]) as f64;
        let ax_y = (q[0][2] - q[2][0]) as f64;
        let ax_z = (q[1][0] - q[0][1]) as f64;
        let n = (ax_x * ax_x + ax_y * ax_y + ax_z * ax_z).sqrt();
        if n < 0.5 {
            return None;
        }
        [ax_x / n, ax_y / n, ax_z / n]
    };
    Some([
        cos_half,
        sin_half * axis[0],
        sin_half * axis[1],
        sin_half * axis[2],
    ])
}

fn compute_signed_perm_spin_parity(
    q: &[[i32; 3]; 3],      // Q = det(P)·P (proper rotation)
    sq_rot: &[[i32; 3]; 3], // b² rotation in G/MSG frame
    g_spin_rots: &[i32],    // G spin table (9 per op)
    g_spin_su2: &[f64],     // G spin SU(2) (4 per op)
    h_spin_rots: &[i32],    // H spin table
    h_spin_su2: &[f64],     // H spin SU(2)
) -> Option<f64> {
    let n_g = g_spin_rots.len() / 9;
    let n_h = h_spin_rots.len() / 9;

    // 1. Find U_G(sq_rot) in G spin table
    let g_idx = (0..n_g).find(|&i| {
        let off = i * 9;
        g_spin_rots[off..off + 9]
            == [
                sq_rot[0][0],
                sq_rot[0][1],
                sq_rot[0][2],
                sq_rot[1][0],
                sq_rot[1][1],
                sq_rot[1][2],
                sq_rot[2][0],
                sq_rot[2][1],
                sq_rot[2][2],
            ]
    })?;
    let u_g = &[
        g_spin_su2[g_idx * 4],
        g_spin_su2[g_idx * 4 + 1],
        g_spin_su2[g_idx * 4 + 2],
        g_spin_su2[g_idx * 4 + 3],
    ];
    // P = det(P)·Q → since det(P) = -1, P = -Q
    // P·R·P⁻¹ = (-Q)·R·(-Q)⁻¹ = Q·R·Q⁻¹ (the -1 factors cancel)
    // For signed-permutation Q with det=1: Q⁻¹ = Q^T
    let q_mat: [[i32; 3]; 3] = [
        [q[0][0], q[0][1], q[0][2]],
        [q[1][0], q[1][1], q[1][2]],
        [q[2][0], q[2][1], q[2][2]],
    ];
    let r_mat: [[i32; 3]; 3] = *sq_rot;
    let q_inv: [[i32; 3]; 3] = [
        [q_mat[0][0], q_mat[1][0], q_mat[2][0]],
        [q_mat[0][1], q_mat[1][1], q_mat[2][1]],
        [q_mat[0][2], q_mat[1][2], q_mat[2][2]],
    ];
    let qr = crate::mathfunc::mat_multiply_matrix_i3(&q_mat, &r_mat);
    let r_h = crate::mathfunc::mat_multiply_matrix_i3(&qr, &q_inv);

    // 2. Find U_H(R') in H spin table
    let h_idx = (0..n_h).find(|&i| {
        let off = i * 9;
        h_spin_rots[off..off + 9]
            == [
                r_h[0][0], r_h[0][1], r_h[0][2], r_h[1][0], r_h[1][1], r_h[1][2], r_h[2][0],
                r_h[2][1], r_h[2][2],
            ]
    })?;
    let u_h = &[
        h_spin_su2[h_idx * 4],
        h_spin_su2[h_idx * 4 + 1],
        h_spin_su2[h_idx * 4 + 2],
        h_spin_su2[h_idx * 4 + 3],
    ];

    let u_q = signed_perm_to_quat(q)?;

    // Transform G lift: U_G→H = U_Q · U_G · U_Q⁻¹
    let u_q_inv = quat_conj(&u_q);
    let t1 = su2_compose(&u_q, u_g);
    let u_g_to_h = su2_compose(&t1, &u_q_inv);

    // Compare with H lift
    su2_lift_relation(&u_g_to_h, u_h).map(|rel| match rel {
        LiftRelation::Same => 1.0,
        LiftRelation::EBar => -1.0,
    })
}

/// Negate a Pauli coefficient vector (multiply by Ebar).
#[inline]
fn neg_pauli(v: &[f64; 4]) -> [f64; 4] {
    [-v[0], -v[1], -v[2], -v[3]]
}

/// Complex conjugate in Pauli convention: U = u0*I + i(u1*sx + u2*sy + u3*sz).
/// U* = u0*I + i((-u1)*sx + u2*sy + (-u3)*sz).
#[inline]
pub(crate) fn conj_pauli(v: &[f64; 4]) -> [f64; 4] {
    [v[0], -v[1], v[2], -v[3]]
}

/// Quaternion conjugate (inverse for unit quaternion): [u0, -u1, -u2, -u3].
#[inline]
fn quat_conj(v: &[f64; 4]) -> [f64; 4] {
    [v[0], -v[1], -v[2], -v[3]]
}

/// Antiunitary square for spin-1/2: A = Theta U = J K U, A^2 = (J U)(J U)*.
/// J = i*sigma_y = [0,0,1,0] in Pauli convention.
#[inline]
fn antiunitary_square_pauli(u: &[f64; 4]) -> [f64; 4] {
    let j = [0.0, 0.0, 1.0, 0.0];
    let ju = su2_compose(&j, u);
    su2_compose(&ju, &conj_pauli(&ju))
}

/// Kernel for antiunitary square in spinor Wigner test.
#[derive(Debug, Clone, Copy)]
pub enum SquareKernel {
    /// Old: U^2 (treats antiunitary as ordinary unitary square)
    OldU2,
    /// J-left: (J U)(J U)* with J = i*sigma_y
    JLeft,
}

impl SquareKernel {
    pub fn apply(&self, u: &[f64; 4]) -> [f64; 4] {
        match self {
            SquareKernel::OldU2 => su2_compose(u, u),
            SquareKernel::JLeft => antiunitary_square_pauli(u),
        }
    }
}

/// Find a Seitz operation in a spin database, preferring full Seitz match.
/// Returns `(index, is_minus_r_fallback)`.
pub(crate) fn find_spin_in_db(op: &SeitzOp, spin_seitz: &[SeitzOp]) -> Option<(usize, bool)> {
    // 1. Full Seitz match
    if let Some(idx) = spin_seitz
        .iter()
        .position(|s| same_seitz_mod_lattice(op, s))
    {
        return Some((idx, false));
    }
    // 2. Unique rotation match
    let rot_matches: Vec<usize> = spin_seitz
        .iter()
        .enumerate()
        .filter_map(|(i, s)| if s.rot == op.rot { Some(i) } else { None })
        .collect();
    if rot_matches.len() == 1 {
        return Some((rot_matches[0], false));
    }
    // 3. -R fallback
    let r_minus: Mat3I = [
        [-op.rot[0][0], -op.rot[0][1], -op.rot[0][2]],
        [-op.rot[1][0], -op.rot[1][1], -op.rot[1][2]],
        [-op.rot[2][0], -op.rot[2][1], -op.rot[2][2]],
    ];
    let minus_matches: Vec<usize> = spin_seitz
        .iter()
        .enumerate()
        .filter_map(|(i, s)| if s.rot == r_minus { Some(i) } else { None })
        .collect();
    if minus_matches.len() == 1 {
        return Some((minus_matches[0], true));
    }
    None
}

/// Infer the central parity eta_ebar = chi(Ebar) / chi(E) for a spinor irrep.
///
/// Looks for both identity lifts (E = [1,0,0,0] and Ebar = [-1,0,0,0]) in
/// `spin_lg_op_indices`.  If both are present, returns eta = chi(Ebar)/chi(E),
/// which is +1.0 for single-valued and -1.0 for genuine spinor irreps.
/// Returns `None` if one or both lifts are missing from the LG character table.
pub fn infer_eta_ebar(
    spin_chars: &[f64],
    spin_lg_op_indices: &[u16],
    h_spin_seitz: &[SeitzOp],
    h_spin_su2: &[f64],
) -> Option<f64> {
    let id_rot: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let u_e = [1.0, 0.0, 0.0, 0.0];
    let u_ebar = [-1.0, 0.0, 0.0, 0.0];

    let mut chi_e: Option<f64> = None;
    let mut chi_ebar: Option<f64> = None;

    for (local, &global) in spin_lg_op_indices.iter().enumerate() {
        let si = global as usize;
        let sop = h_spin_seitz.get(si)?;
        if sop.rot != id_rot {
            continue;
        }
        let u = spin_su2_at(h_spin_su2, si)?;

        if su2_same_up_to_sign(&u, &u_e) == Some(false) {
            chi_e = spin_chars.get(local).copied();
        }
        if su2_same_up_to_sign(&u, &u_ebar) == Some(false) {
            chi_ebar = spin_chars.get(local).copied();
        }
    }

    match (chi_e, chi_ebar) {
        (Some(e), Some(eb)) if e.abs() > 1e-9 => {
            let eta = eb / e;
            if (eta - 1.0).abs() < 1e-6 {
                Some(1.0)
            } else if (eta + 1.0).abs() < 1e-6 {
                Some(-1.0)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Spin-lift context for the Wigner test on spinor irreps.
///
/// For black-white (Type III) magnetic space groups $$\mathcal{M} = H \cup a_0 H$$,
/// the anti-unitary representative $$a_0 = \mathcal{T} g_0$$ has spatial part
/// $$g_0 \in G \setminus H$$ where $$G$$ is the parent space group.  Its SU(2)
/// lift therefore lives in $$G$$'s double group, not $$H$$'s.
///
/// This struct bundles both sets of spin operations so [`wigner_classify_spinor`]
/// can use $$H$$'s lifts for canonical little-group mapping and $$G$$'s lifts
/// for $$a_0$$ lookup.
#[derive(Debug, Clone)]
pub struct SpinLiftContext {
    /// H's spin ops (unitary subgroup): (rotations 9/op, translations 3/op, su2 4/op)
    pub h: (&'static [i32], &'static [f64], &'static [f64]),
    /// G's spin ops (parent spatial group): (rotations 9/op, translations 3/op, su2 4/op)
    /// Same as `h` for grey (Type II) and ordinary (Type I) groups.
    pub g: (&'static [i32], &'static [f64], &'static [f64]),
    /// SG number of H (1-230), for looking up ISOTROPY setting data.
    pub sg: u8,
}

/// Symmetry-operation inputs shared by Wigner classification paths.
#[derive(Debug, Clone, Copy)]
pub struct WignerGroupContext<'a> {
    /// Unitary magnetic-operation indices belonging to the little group.
    pub unitary_indices: &'a [usize],
    /// Magnetic operations in the magnetic database setting.
    pub magnetic_ops: &'a [SeitzOp],
    /// Operations of the unitary subgroup in its Hall setting.
    pub unitary_ops: &'a [SeitzOp],
    /// Index of the selected antiunitary coset representative.
    pub antiunitary_representative: usize,
}

/// Spinor character data and its rational wave vector.
#[derive(Debug, Clone, Copy)]
pub struct SpinorWignerInput<'a> {
    /// Real character components in little-group order.
    pub characters_real: &'a [f64],
    /// Imaginary character components in little-group order.
    pub characters_imag: &'a [f64],
    /// Local character position to global spin-operation index.
    pub operation_indices: &'a [u16],
    /// Rational reciprocal wave vector.
    pub k_vector: KVector,
}

macro_rules! debug_log {
    ($($arg:tt)*) => {
        #[cfg(feature = "debug-corep")]
        eprintln!($($arg)*);
    };
}

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
        let t = [
            if t[0] < 0.0 { t[0] + 1.0 } else { t[0] },
            if t[1] < 0.0 { t[1] + 1.0 } else { t[1] },
            if t[2] < 0.0 { t[2] + 1.0 } else { t[2] },
        ];
        SeitzOp {
            rot,
            trans: t,
            timerev,
        }
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
        if t[i] < 0.0 {
            t[i] += 1.0;
            lattice[i] -= 1;
        }
    }

    (
        SeitzOp {
            rot,
            trans: t,
            timerev,
        },
        lattice,
    )
}

/// Square a Seitz operation: g² = g ∘ g.
///
/// For $$g = \{R|\mathbf{t}\}$$:
/// $$g^2 = \{R^2 \mid \mathbf{t} + R\mathbf{t}\}$$
pub fn square_seitz(g: &SeitzOp) -> (SeitzOp, [i32; 3]) {
    compose_seitz(g, g)
}

// ── Convert MagneticOps → Vec<SeitzOp> ──────────────────────────────────────

/// Convert `SymmetryOps` to a `Vec<SeitzOp>`.
pub fn ops_to_seitz(ops: &SymmetryOps) -> Vec<SeitzOp> {
    (0..ops.len())
        .map(|i| SeitzOp::new(ops[i].rotation, ops[i].translation, ops[i].time_reversal))
        .collect()
}

// ── Little group filter ──────────────────────────────────────────────────────

/// Filter magnetic operations to those that preserve the wave-vector.
///
/// For a **unitary** operation $$\{R|\mathbf{t}\}$$, the condition is
/// $$R^{-T}\mathbf{k} \equiv \mathbf{k} \pmod{\text{reciprocal lattice}}$$.
///
/// For an **anti-unitary** operation $$a = \mathcal{T}\{R|\mathbf{t}\}$$,
/// time reversal sends $$\mathbf{k} \to -\mathbf{k}$$, so the condition is
/// $$-R^{-T}\mathbf{k} \equiv \mathbf{k} \pmod{\text{reciprocal lattice}}$$.
///
/// In terms of integer components with denominator $$k_d$$:
///
/// ```text
/// unitary:     (R⁻ᵀ·k - k) ∈ reciprocal lattice
/// antiunitary: (-R⁻ᵀ·k - k) ∈ reciprocal lattice
/// ```
///
/// For centered conventional cells, integer components alone are not
/// sufficient: the reciprocal shift must also have integer phase against
/// every pure translation in the unitary translation subgroup.
pub fn filter_little_group(kx: i8, ky: i8, kz: i8, kd: i8, ops: &SymmetryOps) -> Vec<usize> {
    filter_little_group_with_transform(kx, ky, kz, kd, ops, None, None)
}

/// Filter operations that preserve the k-vector.
///
/// When `canonical_pure_translations` is provided (from the canonical H Hall
/// setting), the k-preservation phase check uses the FULL translation subgroup.
/// This is essential for centered groups (F/I/C/A) where the MSG-derived pure
/// translations may be only a subset, causing antiunitary ops to incorrectly
/// pass the filter.
pub fn filter_little_group_with_transform(
    kx: i8,
    ky: i8,
    kz: i8,
    kd: i8,
    ops: &SymmetryOps,
    setting_xf: Option<&SettingTransform>,
    canonical_pure_translations: Option<&[[f64; 3]]>,
) -> Vec<usize> {
    if kd == 0 {
        return (0..ops.len()).collect();
    }

    // ── Atomic transform validation ──────────────────────────────────────
    // Before building the transformed op list, verify that the setting
    // transform is valid for ALL operations.  If any single operation fails
    // to transform (non-integer rotation), reject the entire transform and
    // fall back to the raw MSG-frame coordinates.  This prevents mixed-frame
    // little groups where some ops are in Hall frame and others in MSG frame.
    if let Some(xf) = setting_xf {
        let all_ok = ops.operations.iter().all(|op| {
            xf.transform_rotation(&op.rotation).is_some()
                && xf
                    .transform_translation(&op.rotation, &op.translation)
                    .is_some()
        });
        if !all_ok {
            // Transform is not universally valid — use MSG frame.
            debug_log!("  LG filter: setting transform invalid, using MSG frame");
            return filter_little_group_with_transform(
                kx,
                ky,
                kz,
                kd,
                ops,
                None,
                canonical_pure_translations,
            );
        }
    }

    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let transformed: Vec<(Mat3I, [f64; 3], bool)> = ops
        .operations
        .iter()
        .map(|op| {
            if let Some(xf) = setting_xf {
                // SAFETY: validated above that all ops transform successfully.
                let rot = xf.transform_rotation(&op.rotation).unwrap();
                let trans = xf
                    .transform_translation(&op.rotation, &op.translation)
                    .unwrap();
                (rot, trans, op.time_reversal)
            } else {
                (op.rotation, op.translation, op.time_reversal)
            }
        })
        .collect();

    // Use canonical H translation subgroup when available (centered groups).
    // The MSG-derived pure translations may be a subset, causing weaker
    // k-preservation filtering that incorrectly admits antiunitary ops.
    let pure_translations: Vec<[f64; 3]> = match canonical_pure_translations {
        Some(t) => t.to_vec(),
        None => transformed
            .iter()
            .filter(|(rotation, _, time_reversal)| !time_reversal && *rotation == identity)
            .map(|(_, translation, _)| *translation)
            .collect(),
    };

    (0..ops.len())
        .filter(|&i| {
            seitz_preserves_k(
                &transformed[i].0,
                transformed[i].2,
                &pure_translations,
                kx,
                ky,
                kz,
                kd,
            )
        })
        .collect()
}

fn inverse_transpose_unimodular(r: &Mat3I) -> Option<Mat3I> {
    let det = mat_get_determinant_i3(r);
    if det != 1 && det != -1 {
        return None;
    }
    let cofactors = [
        [
            r[1][1] * r[2][2] - r[1][2] * r[2][1],
            -(r[1][0] * r[2][2] - r[1][2] * r[2][0]),
            r[1][0] * r[2][1] - r[1][1] * r[2][0],
        ],
        [
            -(r[0][1] * r[2][2] - r[0][2] * r[2][1]),
            r[0][0] * r[2][2] - r[0][2] * r[2][0],
            -(r[0][0] * r[2][1] - r[0][1] * r[2][0]),
        ],
        [
            r[0][1] * r[1][2] - r[0][2] * r[1][1],
            -(r[0][0] * r[1][2] - r[0][2] * r[1][0]),
            r[0][0] * r[1][1] - r[0][1] * r[1][0],
        ],
    ];
    Some(cofactors.map(|row| row.map(|value| value / det)))
}

fn seitz_preserves_k(
    r: &Mat3I,
    time_reversal: bool,
    pure_translations: &[[f64; 3]],
    kx: i8,
    ky: i8,
    kz: i8,
    kd: i8,
) -> bool {
    if kd == 0 {
        return true;
    }
    let Some(reciprocal_rotation) = inverse_transpose_unimodular(r) else {
        return false;
    };
    let kd_i = kd as i32;
    let k = [kx as i32, ky as i32, kz as i32];
    let sign = if time_reversal { -1 } else { 1 };
    let mut reciprocal_shift = [0i32; 3];
    for i in 0..3 {
        let transformed = sign
            * (reciprocal_rotation[i][0] * k[0]
                + reciprocal_rotation[i][1] * k[1]
                + reciprocal_rotation[i][2] * k[2]);
        let delta = transformed - k[i];
        if delta % kd_i != 0 {
            return false;
        }
        reciprocal_shift[i] = delta / kd_i;
    }

    pure_translations.iter().all(|translation| {
        let phase = reciprocal_shift[0] as f64 * translation[0]
            + reciprocal_shift[1] as f64 * translation[1]
            + reciprocal_shift[2] as f64 * translation[2];
        (phase - phase.round()).abs() < 1e-8
    })
}

// ── Lattice arithmetic helpers ──────────────────────────────────────────────

/// Multiply a 3×3 integer matrix by a 3D integer vector.
#[inline]
pub fn mat_vec_i32(r: &Mat3I, v: &[i32; 3]) -> [i32; 3] {
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

/// Add two [i32; 3] vectors.
#[inline]
pub fn add3(a: &[i32; 3], b: &[i32; 3]) -> [i32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ── Seitz matching ───────────────────────────────────────────────────────────

/// Result of matching a computed Seitz operation to a stored one.
#[derive(Clone, Copy)]
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
pub fn find_seitz(rot: &Mat3I, trans: &[f64; 3], ops: &[SeitzOp]) -> Option<SeitzMatch> {
    for (idx, op) in ops.iter().enumerate() {
        if op.rot[0][0] != rot[0][0]
            || op.rot[0][1] != rot[0][1]
            || op.rot[0][2] != rot[0][2]
            || op.rot[1][0] != rot[1][0]
            || op.rot[1][1] != rot[1][1]
            || op.rot[1][2] != rot[1][2]
            || op.rot[2][0] != rot[2][0]
            || op.rot[2][1] != rot[2][1]
            || op.rot[2][2] != rot[2][2]
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

        if (d0 - l0).abs() < SEITZ_TRANS_TOL
            && (d1 - l1).abs() < SEITZ_TRANS_TOL
            && (d2 - l2).abs() < SEITZ_TRANS_TOL
        {
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
    bloch_phase_f64(
        kx,
        ky,
        kz,
        kd,
        &[lattice[0] as f64, lattice[1] as f64, lattice[2] as f64],
    )
}

/// Compute the Bloch phase for a general translation, including fractional
/// centering and nonsymmorphic shifts between equivalent operation
/// representatives.
pub fn bloch_phase_f64(kx: i8, ky: i8, kz: i8, kd: i8, translation: &[f64; 3]) -> Complex64 {
    if kd == 0 {
        return Complex64::new(1.0, 0.0);
    }
    let theta = 2.0
        * std::f64::consts::PI
        * (kx as f64 * translation[0] + ky as f64 * translation[1] + kz as f64 * translation[2])
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
    group: WignerGroupContext<'_>,
    k_vector: KVector,
) -> Result<CorepType, WignerClassificationError> {
    let WignerGroupContext {
        unitary_indices: unitary_mag_indices,
        magnetic_ops: mag_seitz,
        unitary_ops: h_seitz,
        antiunitary_representative: a0_idx,
    } = group;
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    if unitary_mag_indices.is_empty() {
        return Err(WignerClassificationError::new(
            "unitary little group is empty; Wigner indicator is undefined",
        ));
    }
    if h_chars.len() != h_seitz.len() {
        return Err(WignerClassificationError::new(format!(
            "character table length {} does not match unitary little-group order {}",
            h_chars.len(),
            h_seitz.len()
        )));
    }
    let a0 = mag_seitz.get(a0_idx).ok_or_else(|| {
        WignerClassificationError::new(format!(
            "antiunitary representative index {} is out of range",
            a0_idx
        ))
    })?;
    if !a0.timerev {
        return Err(WignerClassificationError::new(
            "antiunitary representative must have time reversal",
        ));
    }

    let mut w_sum: f64 = 0.0;

    for &h_mag_idx in unitary_mag_indices {
        let h = mag_seitz.get(h_mag_idx).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "unitary operation index {} is out of range",
                h_mag_idx
            ))
        })?;
        if h.timerev {
            return Err(WignerClassificationError::new(format!(
                "operation index {} must be unitary for the Wigner sum",
                h_mag_idx
            )));
        }

        let g0_spatial = SeitzOp::new(a0.rot, a0.trans, false);
        let h_spatial = SeitzOp::new(h.rot, h.trans, false);
        let (g0h, l1) = compose_seitz(&g0_spatial, &h_spatial);
        let (sq, lattice_sq) = square_seitz(&g0h);

        let matched = find_seitz(&sq.rot, &sq.trans, h_seitz).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "square of (a0·h) for unitary operation index {} is absent from the unitary little group",
                h_mag_idx
            ))
        })?;
        debug_assert!(matched.op_index < h_chars.len());

        // Total lattice shift = L_sq + L_match + L1 + R_{g0h}·L1
        let r_l1 = mat_vec_i32(&g0h.rot, &l1);
        let total_lattice = add3(
            &add3(&lattice_sq, &matched.lattice_shift),
            &add3(&l1, &r_l1),
        );
        let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
        let contrib = h_chars[matched.op_index] * phase.re;
        w_sum += contrib;
        debug_log!(
            "    wigner: h[{}]→H[{}] L={:?} ph={:.2} χ={:.2} → {:.2}",
            h_mag_idx,
            matched.op_index,
            total_lattice,
            phase.re,
            h_chars[matched.op_index],
            contrib
        );
    }

    let n = unitary_mag_indices.len() as f64;
    let w = w_sum / n;

    // Strict classification: W must be quantized to 0, +1, or -1.
    debug_log!(
        "DEBUG wigner_classify: w_sum={:.4} n_unitary={} W={:.4} k=({},{},{})/{}",
        w_sum,
        unitary_mag_indices.len(),
        w,
        kx,
        ky,
        kz,
        kd
    );
    let tol = 1e-6;
    if (w - 1.0).abs() < tol {
        Ok(CorepType::A)
    } else if (w + 1.0).abs() < tol {
        Ok(CorepType::B)
    } else if w.abs() < tol {
        Ok(CorepType::C)
    } else {
        debug_log!(
            "  Non-quantized Wigner indicator W={:.8}; expected 0, +1, or -1.",
            w
        );
        Err(WignerClassificationError::with_value(
            "Wigner indicator is non-quantized; expected 0, +1, or -1",
            w,
        ))
    }
}

// ── CIR-based Wigner test ───────────────────────────────────────────────────

/// Wigner test using CIR (complex) character tables.
///
/// For compound PIR irreps like Z1Z4 = Z1 ⊕ Z4, the underlying CIR
/// irreps Z1, Z4 are complex and may individually give Type C under
/// the antiunitary operation, even though the combined PIR gives Type A.
///
/// This function classifies one CIR component using a selected antiunitary
/// representative. Callers that already have the full antiunitary little-group
/// coset should prefer [`wigner_classify_cir_direct`].
///
/// A component-based caller should:
/// 1. Check `irrep.cir_component_count() > 0`
/// 2. Loop over components, call this function for each
/// 3. If any component gives Type C, the overall corep is Type C
///
/// $$ W = \frac{1}{|H|} \sum_{h \in H} \chi_{\text{CIR}}((a_0 h)^2) $$
///
/// where $$\chi_{\text{CIR}}$$ is complex-valued.  W is complex; we
/// classify by $$|W| < 0.01$$ → Type C, Re(W) > 0 → Type A, else Type B.
pub fn wigner_classify_cir(
    cir_chars: &[f64], // (re, im) pairs for one CIR component
    group: WignerGroupContext<'_>,
    k_vector: KVector,
) -> Result<CorepType, WignerClassificationError> {
    let WignerGroupContext {
        unitary_indices: unitary_mag_indices,
        magnetic_ops: mag_seitz,
        unitary_ops: h_seitz,
        antiunitary_representative: a0_idx,
    } = group;
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    if unitary_mag_indices.is_empty() {
        return Err(WignerClassificationError::new(
            "unitary little group is empty; Wigner indicator is undefined",
        ));
    }
    if cir_chars.len() != 2 * h_seitz.len() {
        return Err(WignerClassificationError::new(format!(
            "CIR character count {} does not match unitary little-group order {}",
            cir_chars.len() / 2,
            h_seitz.len()
        )));
    }
    let a0 = mag_seitz.get(a0_idx).ok_or_else(|| {
        WignerClassificationError::new(format!(
            "antiunitary representative index {} is out of range",
            a0_idx
        ))
    })?;
    if !a0.timerev {
        return Err(WignerClassificationError::new(
            "antiunitary representative must have time reversal",
        ));
    }
    let mut w_sum = Complex64::new(0.0, 0.0);
    #[cfg(feature = "debug-corep")]
    let mut n_plus = 0u32;
    #[cfg(feature = "debug-corep")]
    let mut n_minus = 0u32;

    for &h_mag_idx in unitary_mag_indices {
        let h = mag_seitz.get(h_mag_idx).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "unitary operation index {} is out of range",
                h_mag_idx
            ))
        })?;
        if h.timerev {
            return Err(WignerClassificationError::new(format!(
                "operation index {} must be unitary for the Wigner sum",
                h_mag_idx
            )));
        }
        let g0_spatial = SeitzOp::new(a0.rot, a0.trans, false);
        let h_spatial = SeitzOp::new(h.rot, h.trans, false);
        let (g0h, l1) = compose_seitz(&g0_spatial, &h_spatial);
        let (sq, lattice_sq) = square_seitz(&g0h);

        let matched = find_seitz(&sq.rot, &sq.trans, h_seitz).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "square of (a0·h) for unitary operation index {} is absent from the unitary little group",
                h_mag_idx
            ))
        })?;
        let r_l1 = mat_vec_i32(&g0h.rot, &l1);
        let total_lattice = add3(
            &add3(&lattice_sq, &matched.lattice_shift),
            &add3(&l1, &r_l1),
        );
        let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
        let chi = cir_char_at(cir_chars, matched.op_index);
        w_sum += chi * phase;
        // Phase parity stats
        #[cfg(feature = "debug-corep")]
        if phase.re > 0.5 {
            n_plus += 1;
        } else if phase.re < -0.5 {
            n_minus += 1;
        }
        debug_log!(
            "    cir: h[{}]→H[{}] Lz_par={} ph={:.2} χ={:.2} → {:.2}",
            h_mag_idx,
            matched.op_index,
            ((total_lattice[2] % 2) + 2) % 2,
            phase,
            chi,
            chi * phase
        );
    }

    debug_log!("    phase stats: +={} -={}", n_plus, n_minus);
    let n = unitary_mag_indices.len() as f64;
    let w = w_sum / n;
    debug_log!(
        "DEBUG wigner_classify_cir: W=({:.8},{:.8}) |W|={:.4} k=({},{},{})/{}",
        w.re,
        w.im,
        w.norm(),
        kx,
        ky,
        kz,
        kd
    );

    let tol = 1e-6;
    if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        Ok(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        Ok(CorepType::B)
    } else if w.norm() < tol {
        Ok(CorepType::C)
    } else {
        Err(WignerClassificationError::with_value(
            "Non-quantized Wigner indicator; expected 0, +1, or -1",
            w.norm(),
        ))
    }
}

/// Wigner test evaluated directly over the antiunitary little-group coset.
///
/// This is algebraically equivalent to the `a0 * H` form used by
/// [`wigner_classify_cir`], but it avoids choosing a coset representative and
/// repeatedly reducing the intermediate product.  That distinction matters
/// for nonsymmorphic boundary points, where losing an intermediate lattice
/// translation changes the Bloch phase.
pub fn wigner_classify_cir_direct(
    cir_chars: &[f64],
    antiunitary_mag_indices: &[usize],
    mag_seitz: &[SeitzOp],
    h_seitz: &[SeitzOp],
    k_vector: KVector,
) -> Result<CorepType, WignerClassificationError> {
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    if antiunitary_mag_indices.is_empty() {
        return Err(WignerClassificationError::new(
            "antiunitary little-group coset is empty",
        ));
    }
    if cir_chars.len() != 2 * h_seitz.len() {
        return Err(WignerClassificationError::new(format!(
            "CIR character count {} does not match unitary little-group order {}",
            cir_chars.len() / 2,
            h_seitz.len()
        )));
    }

    let mut sum = Complex64::ZERO;
    for &b_idx in antiunitary_mag_indices {
        let b = mag_seitz.get(b_idx).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "antiunitary operation index {} is out of range",
                b_idx
            ))
        })?;
        let (square, square_lattice) = square_seitz(b);
        let matched = find_seitz(&square.rot, &square.trans, h_seitz).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "square of antiunitary operation {} is absent from the unitary little group",
                b_idx
            ))
        })?;
        let total_lattice = add3(&square_lattice, &matched.lattice_shift);
        let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
        sum += cir_char_at(cir_chars, matched.op_index) * phase;
    }

    let w = sum / antiunitary_mag_indices.len() as f64;
    debug_log!(
        "DEBUG wigner_classify_cir_direct: W=({:.8},{:.8}) |W|={:.4} k=({},{},{})/{}",
        w.re,
        w.im,
        w.norm(),
        kx,
        ky,
        kz,
        kd
    );

    let tol = 1e-6;
    if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        Ok(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        Ok(CorepType::B)
    } else if w.norm() < tol {
        Ok(CorepType::C)
    } else {
        Err(WignerClassificationError::with_value(
            "Non-quantized direct antiunitary-coset Wigner indicator; expected 0, +1, or -1",
            w.norm(),
        ))
    }
}

// ── Phase 1: Setting-transform oracle (2026-06-19) ─────────────────────────

/// A basis-origin transformation between coordinate frames.
///
/// Convention (matching `magnetic_spacegroup.rs`):
///
/// ```text
/// x_hall = T · x_msg + s
/// R_hall = T · R_msg · T⁻¹
/// t_hall = s - R_hall·s + T·t_msg   (mod Z³)
/// ```
#[derive(Debug, Clone)]
pub struct SettingTransform {
    pub basis: Mat3,
    pub origin: [f64; 3],
}

impl SettingTransform {
    pub fn identity() -> Self {
        SettingTransform {
            basis: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            origin: [0.0; 3],
        }
    }

    /// Compose two transforms: `self` (A→B) followed by `other` (B→C),
    /// giving A→C.  The formula follows from `x_B = P_1·x_A + s_1` and
    /// `x_C = P_2·x_B + s_2`:
    ///
    /// ```text
    /// P_total = P_2 · P_1
    /// s_total = P_2 · s_1 + s_2   (mod Z³)
    /// ```
    ///
    /// Translation components are taken modulo 1 since they represent
    /// origin shifts in fractional coordinates.
    pub fn then(&self, other: &SettingTransform) -> SettingTransform {
        let p1 = self.basis;
        let s1 = self.origin;
        let p2 = other.basis;
        let s2 = other.origin;
        // P_total = P_2 · P_1
        let basis = mat_multiply_matrix_d3(&p2, &p1);
        // s_total = P_2 · s_1 + s_2 (mod Z³)
        let mut origin = [0.0f64; 3];
        for i in 0..3 {
            let ps: f64 = (0..3).map(|j| p2[i][j] * s1[j]).sum();
            origin[i] = (ps + s2[i]) % 1.0;
            if origin[i] < 0.0 {
                origin[i] += 1.0;
            }
        }
        SettingTransform { basis, origin }
    }

    /// Apply the forward transform to a rotation matrix.
    /// Returns `None` when the setting-transform basis is not invertible, or
    /// when the transform does not produce an integer rotation (i.e. the
    /// setting transform is invalid for this MSG frame).
    pub fn transform_rotation(&self, r_msg: &Mat3I) -> Option<Mat3I> {
        let t = self.basis;
        let t_inv = mat_inverse_matrix_d3(&t, 1e-10).ok()?;
        let mut transformed = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                transformed[i][j] = (0..3)
                    .flat_map(|a| (0..3).map(move |b| t[i][a] * r_msg[a][b] as f64 * t_inv[b][j]))
                    .sum();
            }
        }
        let mut result = [[0i32; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let rounded = transformed[i][j].round();
                if (transformed[i][j] - rounded).abs() > 1e-8 {
                    return None;
                }
                result[i][j] = rounded as i32;
            }
        }
        Some(result)
    }

    /// Apply the forward transform to a translation vector.
    /// Returns `None` when the rotation transform fails.
    pub fn transform_translation(&self, r_msg: &Mat3I, t_msg: &[f64; 3]) -> Option<[f64; 3]> {
        let r_hall = self.transform_rotation(r_msg)?;
        let s = self.origin;
        // t_hall = s - R_hall·s + T·t_msg
        let mut t_hall = [0.0f64; 3];
        for i in 0..3 {
            let rs: f64 = (0..3).map(|j| r_hall[i][j] as f64 * s[j]).sum();
            let tt: f64 = (0..3).map(|j| self.basis[i][j] * t_msg[j]).sum();
            t_hall[i] = (s[i] - rs + tt) % 1.0;
            if t_hall[i] < 0.0 {
                t_hall[i] += 1.0;
            }
        }
        Some(t_hall)
    }

    /// Transform a full Seitz operation (R|t). Returns `None` when the basis
    /// is not invertible or the rotation does not transform to an integer matrix.
    pub fn transform_seitz(&self, rot: &Mat3I, trans: &[f64; 3]) -> Option<(Mat3I, [f64; 3])> {
        let r_hall = self.transform_rotation(rot)?;
        let t_hall = self.transform_translation(rot, trans)?;
        Some((r_hall, t_hall))
    }
}

/// Multiply two 3×3 i32 matrices: C = A·B.
fn mat_multiply_3i_3i(a: &[[i32; 3]; 3], b: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    let mut c = [[0i32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

/// Inverse of a signed-permutation 3×3 matrix (det=±1, integer entries).
fn mat_inverse_3i(m: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    // For signed permutations, inverse = (1/det) * transpose of cofactor matrix.
    // Since det = ±1 and entries are integers, the inverse is also integer.
    let inv_det = det; // det = ±1
    let mut inv = [[0i32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let a = m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3]
                - m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3];
            inv[i][j] = inv_det * a;
        }
    }
    inv
}

/// Enumerate all 48 signed-permutation 3×3 matrices.
///
/// Each matrix has exactly one non-zero entry per row and column, with value ±1.
/// These are the orthogonal unimodular matrices — the only 3×3 integer matrices
/// with determinant ±1 whose inverse is also integer and whose rows are orthogonal.
pub fn enumerate_signed_permutations() -> Vec<[[i32; 3]; 3]> {
    let mut results = Vec::with_capacity(48);
    let signs = [
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [1, -1, -1],
        [-1, 1, 1],
        [-1, 1, -1],
        [-1, -1, 1],
        [-1, -1, -1],
    ];
    let perms: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for sign in &signs {
        for perm in &perms {
            let mut m = [[0i32; 3]; 3];
            for (row, &col) in perm.iter().enumerate() {
                m[row][col] = sign[row];
            }
            results.push(m);
        }
    }
    results
}

/// All unimodular 3×3 integer matrices with entries in {-1, 0, 1}.
///
/// Includes the 48 signed-permutations plus crystallographic shear
/// transformations (GL(3,Z) shears).  These are needed for Hall-setting
/// pairs that differ by more than axis permutation, e.g. monoclinic
/// unique-axis conversions (Hall 73↔72 for P2/c).
/// Cached unimodular bases — computed once, reused across all
/// `find_setting_transform` calls.  Avoids generating 6,960 matrices
/// per call (the main performance regression noted by codex review).
static UNIMODULAR_BASES: OnceLock<Vec<[[i32; 3]; 3]>> = OnceLock::new();

pub fn enumerate_unimodular_bases() -> &'static Vec<[[i32; 3]; 3]> {
    UNIMODULAR_BASES.get_or_init(compute_unimodular_bases)
}

fn compute_unimodular_bases() -> Vec<[[i32; 3]; 3]> {
    // Start with the 48 signed-permutations.
    let mut results = enumerate_signed_permutations();
    let seen: std::collections::HashSet<[[i32; 3]; 3]> = results.iter().copied().collect();
    // Generate all matrices with entries in {-1, 0, 1} and determinant ±1.
    let vals = [-1i32, 0, 1];
    for a00 in &vals {
        for a01 in &vals {
            for a02 in &vals {
                for a10 in &vals {
                    for a11 in &vals {
                        for a12 in &vals {
                            for a20 in &vals {
                                for a21 in &vals {
                                    for a22 in &vals {
                                        let m = [
                                            [*a00, *a01, *a02],
                                            [*a10, *a11, *a12],
                                            [*a20, *a21, *a22],
                                        ];
                                        if seen.contains(&m) {
                                            continue;
                                        }
                                        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                                            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                                            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
                                        if det == 1 || det == -1 {
                                            results.push(m);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    results
}

/// Find all setting transforms `(T, s)` that satisfy:
///
/// ```text
/// x_hall = T · x_msg + s
/// R_hall = T · R_msg · T⁻¹
/// t_hall = s - R_hall·s + T·t_msg   (mod Z³)
/// ```
///
/// Returns all valid `(T, s)` pairs.  If the rotation multisets already match
/// without any basis change, the identity `T=I` is tried with origin solving.
pub fn find_setting_transform(
    msg_rots: &[[[i32; 3]; 3]],
    msg_trans: &[[f64; 3]],
    hall_rots: &[[[i32; 3]; 3]],
    hall_trans: &[[f64; 3]],
) -> Vec<SettingTransform> {
    XF_CALLED.fetch_add(1, Ordering::Relaxed);
    let mut results = Vec::new();
    let validate_xf = |xf: &SettingTransform| -> bool {
        msg_rots.iter().zip(msg_trans.iter()).all(|(r, tm)| {
            let Some(xf_r) = xf.transform_rotation(r) else {
                return false;
            };
            let Some(xf_t) = xf.transform_translation(r, tm) else {
                return false;
            };
            hall_rots.iter().zip(hall_trans.iter()).any(|(hr, ht)| {
                *hr == xf_r && {
                    let d0 = xf_t[0] - ht[0];
                    let d1 = xf_t[1] - ht[1];
                    let d2 = xf_t[2] - ht[2];
                    (d0 - d0.round()).abs() < SEITZ_TRANS_TOL
                        && (d1 - d1.round()).abs() < SEITZ_TRANS_TOL
                        && (d2 - d2.round()).abs() < SEITZ_TRANS_TOL
                }
            })
        })
    };

    // The alternate-origin Hall pairs selected by the generated ISO-IR data
    // use eighth-cell fractions (quarters and 3/8 are the common cases). The
    // linear solver below assumes a greedy one-to-one pairing of equal
    // rotations; that pairing is underdetermined for centered settings and
    // can reject a perfectly valid origin. Keep this small exact grid as a
    // deterministic fallback for the identity-basis case.
    let identity_grid_origin = || -> Option<[f64; 3]> {
        let denominator = 8.0;
        for ix in 0..8 {
            for iy in 0..8 {
                for iz in 0..8 {
                    let origin = [
                        ix as f64 / denominator,
                        iy as f64 / denominator,
                        iz as f64 / denominator,
                    ];
                    let xf = SettingTransform {
                        basis: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        origin,
                    };
                    if validate_xf(&xf) {
                        return Some(origin);
                    }
                }
            }
        }
        None
    };

    // Try identity basis first.
    if rotation_multiset_eq(msg_rots, hall_rots) {
        if let Some(s) = solve_origin_for_t(
            &[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            msg_rots,
            msg_trans,
            hall_rots,
            hall_trans,
        ) {
            let xf = SettingTransform {
                basis: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                origin: s,
            };
            if validate_xf(&xf) {
                XF_FOUND.fetch_add(1, Ordering::Relaxed);
                XF_IDENTITY.fetch_add(1, Ordering::Relaxed);
                let origin_nz = s[0].abs() > 1e-8 || s[1].abs() > 1e-8 || s[2].abs() > 1e-8;
                if origin_nz {
                    XF_NONZERO_ORIGIN.fetch_add(1, Ordering::Relaxed);
                }
                results.push(xf);
                return results;
            }
        }
        if let Some(s) = identity_grid_origin() {
            XF_FOUND.fetch_add(1, Ordering::Relaxed);
            XF_IDENTITY.fetch_add(1, Ordering::Relaxed);
            if s.iter().any(|value| value.abs() > 1e-8) {
                XF_NONZERO_ORIGIN.fetch_add(1, Ordering::Relaxed);
            }
            results.push(SettingTransform {
                basis: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                origin: s,
            });
            return results;
        }
        // Identity basis matches rotations but either origin solving or the
        // full Seitz-set validation failed. Let the loop try other bases.
    }

    // Use enumerate_unimodular_bases() instead of signed-permutations only.
    // GL(3,Z) shear matrices (e.g. [[-1,0,-1],[0,-1,0],[0,0,-1]] for
    // monoclinic Hall 73→72) have entries in {-1,0,1} and det=±1, and
    // are included in this expanded candidate pool.
    for t in enumerate_unimodular_bases() {
        let t_inv = mat_inverse_3i(t);
        let xf_rots: Vec<[[i32; 3]; 3]> = msg_rots
            .iter()
            .map(|r| {
                let tmp = mat_multiply_3i_3i(t, r);
                mat_multiply_3i_3i(&tmp, &t_inv)
            })
            .collect();
        if !rotation_multiset_eq(&xf_rots, hall_rots) {
            continue;
        }
        if let Some(s) = solve_origin_for_t(t, msg_rots, msg_trans, hall_rots, hall_trans) {
            let xf = SettingTransform {
                basis: t.map(|row| row.map(|value| value as f64)),
                origin: s,
            };
            if !validate_xf(&xf) {
                continue;
            }
            XF_NON_IDENTITY.fetch_add(1, Ordering::Relaxed);
            let origin_nz = s[0].abs() > 1e-8 || s[1].abs() > 1e-8 || s[2].abs() > 1e-8;
            if origin_nz {
                XF_NONZERO_ORIGIN.fetch_add(1, Ordering::Relaxed);
            }
            results.push(xf);
        }
    }
    // If no candidate passed the origin solver, try zero-origin with full
    // Seitz-set validation.  This catches Hall pairs where the rotation
    // multisets match (same SG, same axis convention) but the origin solver
    // fails due to identity-row consistency not being checked (see
    // solve_origin_for_t).  For these cases, a zero origin may be correct
    // if the only difference is trivial axis permutation.
    if results.is_empty() {
        // Check if rotations already match (no permutation needed).
        if rotation_multiset_eq(msg_rots, hall_rots) {
            let xf = SettingTransform::identity();
            if validate_xf(&xf) {
                results.push(xf);
            }
        }
        // Use enumerate_unimodular_bases() for the same reason as above:
        // GL(3,Z) shears are needed for monoclinic unique-axis conversions.
        for t in enumerate_unimodular_bases() {
            let t_inv = mat_inverse_3i(t);
            let xf_rots: Vec<[[i32; 3]; 3]> = msg_rots
                .iter()
                .map(|r| {
                    let tmp = mat_multiply_3i_3i(t, r);
                    mat_multiply_3i_3i(&tmp, &t_inv)
                })
                .collect();
            if !rotation_multiset_eq(&xf_rots, hall_rots) {
                continue;
            }
            let xf = SettingTransform {
                basis: t.map(|row| row.map(|value| value as f64)),
                origin: [0.0; 3],
            };
            if validate_xf(&xf) {
                results.push(xf);
            }
        }
    }

    if !results.is_empty() {
        XF_FOUND.fetch_add(1, Ordering::Relaxed);
    }
    if results.len() > 1 {
        XF_AMBIGUOUS.fetch_add(1, Ordering::Relaxed);
    }
    results
}

/// For a given basis T, solve the origin s from the translation equations:
///
/// ```text
/// (I - R_hall) · s  ≡  t_hall - T·t_msg   (mod Z³)
/// ```
///
/// Builds a rotation correspondence (msg_idx → hall_idx) via greedy matching,
/// then solves the over-determined modulo-1 linear system by Gaussian
/// elimination over ℝ³ followed by a modulo-consistency check.
fn solve_origin_for_t(
    t: &[[i32; 3]; 3],
    msg_rots: &[[[i32; 3]; 3]],
    msg_trans: &[[f64; 3]],
    hall_rots: &[[[i32; 3]; 3]],
    hall_trans: &[[f64; 3]],
) -> Option<[f64; 3]> {
    // Build greedy rotation correspondence.
    let t_inv = mat_inverse_3i(t);
    let mut hall_used = vec![false; hall_rots.len()];
    let mut pairs: Vec<(usize, usize)> = Vec::new(); // (msg_idx, hall_idx)

    for (mi, mr) in msg_rots.iter().enumerate() {
        let tmp = mat_multiply_3i_3i(t, mr);
        let xf_r = mat_multiply_3i_3i(&tmp, &t_inv);
        let pos = hall_rots
            .iter()
            .enumerate()
            .position(|(j, hr)| *hr == xf_r && !hall_used[j]);
        if let Some(hi) = pos {
            hall_used[hi] = true;
            pairs.push((mi, hi));
        }
    }
    if pairs.is_empty() {
        return None;
    }

    // Collect scalar equations  A·s ≡ b (mod 1).
    // Each operation pair gives up to 3 independent scalar rows.
    let mut eqs: Vec<([f64; 3], f64)> = Vec::with_capacity(3 * pairs.len());

    for &(msg_idx, hall_idx) in &pairs {
        let rh = hall_rots[hall_idx];
        let a_row0 = [1.0 - rh[0][0] as f64, -rh[0][1] as f64, -rh[0][2] as f64];
        let a_row1 = [-rh[1][0] as f64, 1.0 - rh[1][1] as f64, -rh[1][2] as f64];
        let a_row2 = [-rh[2][0] as f64, -rh[2][1] as f64, 1.0 - rh[2][2] as f64];

        let th = hall_trans[hall_idx];
        let tm = msg_trans[msg_idx];
        let tt: [f64; 3] = [
            t[0][0] as f64 * tm[0] + t[0][1] as f64 * tm[1] + t[0][2] as f64 * tm[2],
            t[1][0] as f64 * tm[0] + t[1][1] as f64 * tm[1] + t[1][2] as f64 * tm[2],
            t[2][0] as f64 * tm[0] + t[2][1] as f64 * tm[1] + t[2][2] as f64 * tm[2],
        ];
        let b = [th[0] - tt[0], th[1] - tt[1], th[2] - tt[2]];

        // Only keep rows where A is not trivially zero (R ≠ I for that row).
        // For the identity rotation, A = 0 and we get 0 ≡ b, which is a
        // consistency condition rather than a constraint on s.
        let keep_row =
            |a: &[f64; 3]| a[0].abs() > 1e-12 || a[1].abs() > 1e-12 || a[2].abs() > 1e-12;
        if keep_row(&a_row0) {
            eqs.push((a_row0, b[0]));
        }
        if keep_row(&a_row1) {
            eqs.push((a_row1, b[1]));
        }
        if keep_row(&a_row2) {
            eqs.push((a_row2, b[2]));
        }
    }

    if eqs.is_empty() {
        // No non-trivial constraints — any origin works.
        return Some([0.0; 3]);
    }

    // Gaussian elimination over ℝ (augmented matrix: [A | b]).
    let mut aug: Vec<[f64; 4]> = eqs.iter().map(|(a, b)| [a[0], a[1], a[2], *b]).collect();
    let mut pivot_row = 0usize;

    for col in 0..3 {
        // Find pivot.
        let best = aug
            .iter()
            .enumerate()
            .skip(pivot_row)
            .find_map(|(row, values)| (values[col].abs() > 1e-12).then_some(row));
        let pr = match best {
            Some(r) => r,
            None => continue,
        };
        aug.swap(pivot_row, pr);

        // Normalize pivot row.
        let pv = aug[pivot_row][col];
        for value in aug[pivot_row].iter_mut().skip(col) {
            *value /= pv;
        }

        // Eliminate other rows.
        let pivot = aug[pivot_row];
        for (r, row) in aug.iter_mut().enumerate() {
            if r == pivot_row {
                continue;
            }
            let factor = row[col];
            if factor.abs() < 1e-12 {
                continue;
            }
            for (value, &pivot_value) in row.iter_mut().zip(&pivot).skip(col) {
                *value -= factor * pivot_value;
            }
        }
        pivot_row += 1;
    }

    // Read solution: for each column 0..2, if there's a pivot row with a 1,
    // read s[col] from the RHS; otherwise set to 0.
    let mut s = [0.0f64; 3];
    for col in 0..3 {
        let pr = (0..pivot_row).find(|&r| (aug[r][col] - 1.0).abs() < 1e-12);
        if let Some(r) = pr {
            s[col] = aug[r][3];
        }
        // else: free variable, keep 0
    }

    // Verify: for every original equation, (A·s - b) must be integer.
    for (a, b) in &eqs {
        let residual = a[0] * s[0] + a[1] * s[1] + a[2] * s[2] - b;
        let frac = residual - residual.round();
        if frac.abs() > 1e-8 {
            return None;
        }
    }

    // Normalize s to [0, 1).
    for value in &mut s {
        *value = (*value % 1.0 + 1.0) % 1.0;
        if value.abs() < 1e-10 {
            *value = 0.0;
        }
        if (*value - 1.0).abs() < 1e-10 {
            *value = 0.0;
        }
    }

    Some(s)
}

/// Compare two rotation multisets for equality (order-independent).
pub fn rotation_multiset_eq(a: &[[[i32; 3]; 3]], b: &[[[i32; 3]; 3]]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut b_used = vec![false; b.len()];
    for ra in a {
        let mut found = false;
        for (j, rb) in b.iter().enumerate() {
            if !b_used[j] && ra == rb {
                b_used[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

/// Verify that `h_seitz` operation ordering matches the CIR character table.
/// Prints all operations with their characters for manual inspection.
#[cfg(test)]
pub fn debug_char_order(cir_chars: &[f64], h_seitz: &[SeitzOp], _label: &str) {
    debug_log!("=== Character order check: {} ===", _label);
    for (i, op) in h_seitz.iter().enumerate() {
        let _re = cir_chars.get(2 * i).copied().unwrap_or(999.0);
        let _im = cir_chars.get(2 * i + 1).copied().unwrap_or(999.0);
        let _is_id = op.rot[0][0] == 1
            && op.rot[0][1] == 0
            && op.rot[0][2] == 0
            && op.rot[1][0] == 0
            && op.rot[1][1] == 1
            && op.rot[1][2] == 0
            && op.rot[2][0] == 0
            && op.rot[2][1] == 0
            && op.rot[2][2] == 1
            && op.trans[0].abs() < 0.01
            && op.trans[1].abs() < 0.01
            && op.trans[2].abs() < 0.01;
        debug_log!(
            "  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) chi=({:.3},{:.3}){}",
            i,
            op.rot[0][0],
            op.rot[0][1],
            op.rot[0][2],
            op.rot[1][0],
            op.rot[1][1],
            op.rot[1][2],
            op.rot[2][0],
            op.rot[2][1],
            op.rot[2][2],
            op.trans[0],
            op.trans[1],
            op.trans[2],
            _re,
            _im,
            if _is_id { " ← ID" } else { "" },
        );
    }
}

/// Diagnostic: unwrapped Seitz square without intermediate normalization.
/// Computes (g₀h)² directly from raw translations, then compares
/// with the normalized+matched result.  Used to debug phase parity.
pub fn debug_unwrapped_square(
    h_mag_idx: usize,
    group: WignerGroupContext<'_>,
    k_vector: KVector,
) -> Result<(), WignerClassificationError> {
    let a0_idx = group.antiunitary_representative;
    let mag_seitz = group.magnetic_ops;
    let h_seitz = group.unitary_ops;
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    let a0 = mag_seitz.get(a0_idx).ok_or_else(|| {
        WignerClassificationError::new("diagnostic a0 operation index is out of range")
    })?;
    let h = mag_seitz.get(h_mag_idx).ok_or_else(|| {
        WignerClassificationError::new("diagnostic H operation index is out of range")
    })?;
    if !a0.timerev || h.timerev {
        return Err(WignerClassificationError::new(
            "diagnostic operations have incorrect unitary/antiunitary roles",
        ));
    }

    // Step 1: g₀h raw (no normalization)
    let rc = mat_multiply_matrix_i3(&a0.rot, &h.rot);
    let r0_th = mat_vec_f64(&a0.rot, &h.trans);
    let tc_raw = [
        a0.trans[0] + r0_th[0],
        a0.trans[1] + r0_th[1],
        a0.trans[2] + r0_th[2],
    ];

    // Step 2: (g₀h)² raw
    let rsq = mat_multiply_matrix_i3(&rc, &rc);
    let rc_tc = mat_vec_f64(&rc, &tc_raw);
    let tsq_raw = [
        tc_raw[0] + rc_tc[0],
        tc_raw[1] + rc_tc[1],
        tc_raw[2] + rc_tc[2],
    ];

    debug_log!("=== unwrapped square: h[{}] ===", h_mag_idx);
    debug_log!("  a0: R={:?}, t={:?}", a0.rot, a0.trans);
    debug_log!("  h : R={:?}, t={:?}", h.rot, h.trans);
    debug_log!("  g0h raw: R={:?}, t={:?}", rc, tc_raw);
    debug_log!("  sq raw : R={:?}, t={:?}", rsq, tsq_raw);

    // Normalize for matching
    let (tsq_mod, _l_reduce) = reduce01_with_lattice(&tsq_raw);
    debug_log!("  sq mod : t={:?}, L_reduce={:?}", tsq_mod, _l_reduce);

    if let Some(m) = find_seitz(&rsq, &tsq_mod, h_seitz) {
        let stored_t = &h_seitz[m.op_index].trans;
        // Direct lattice difference: tsq_raw - stored_t
        let l_direct = [
            (tsq_raw[0] - stored_t[0]).round() as i32,
            (tsq_raw[1] - stored_t[1]).round() as i32,
            (tsq_raw[2] - stored_t[2]).round() as i32,
        ];
        let _lz_par = ((l_direct[2] % 2) + 2) % 2;
        let _phase = bloch_phase(kx, ky, kz, kd, &l_direct);

        debug_log!("  matched H[{}]: t_stored={:?}", m.op_index, stored_t);
        debug_log!(
            "  L_direct={:?} Lz_par={} phase={:.2}",
            l_direct,
            _lz_par,
            _phase
        );
        debug_log!(
            "  m.lattice_shift={:?} (from normalized match)",
            m.lattice_shift
        );
    } else {
        debug_log!("  NOT FOUND in h_seitz");
    }
    Ok(())
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
    k_vector: KVector,
) -> Result<Complex64, WignerClassificationError> {
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    let expected_chars = h_seitz
        .len()
        .checked_mul(2)
        .ok_or_else(|| WignerClassificationError::new("CIR character length overflow"))?;
    if anti_lg_indices.is_empty() {
        return Err(WignerClassificationError::new(
            "antiunitary little-group coset is empty",
        ));
    }
    if cir_chars.len() != expected_chars || cir_chars.iter().any(|value| !value.is_finite()) {
        return Err(WignerClassificationError::new(
            "CIR character table does not match H operations",
        ));
    }
    if anti_lg_indices
        .iter()
        .any(|&index| index >= mag_seitz.len())
    {
        return Err(WignerClassificationError::new(
            "antiunitary operation index is out of range",
        ));
    }

    let mut sum = Complex64::ZERO;
    #[cfg(feature = "debug-corep")]
    let mut n_plus = 0u32;
    #[cfg(feature = "debug-corep")]
    let mut n_minus = 0u32;

    for &b_idx in anti_lg_indices {
        let b = &mag_seitz[b_idx];
        let (sq, lattice_sq) = square_seitz(b);
        let m = find_seitz(&sq.rot, &sq.trans, h_seitz).ok_or_else(|| {
            WignerClassificationError::new(format!(
                "direct antiunitary square b[{b_idx}]² was not found in H"
            ))
        })?;

        let total_lattice = add3(&lattice_sq, &m.lattice_shift);
        let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
        let chi = cir_char_at(cir_chars, m.op_index);
        let contrib = chi * phase;
        sum += contrib;

        #[cfg(feature = "debug-corep")]
        if phase.re > 0.5 {
            n_plus += 1;
        } else if phase.re < -0.5 {
            n_minus += 1;
        }

        debug_log!(
            "  direct: b[{}]^2→H[{}] L={:?} ph={:.2} χ={:.2} → {:.2}",
            b_idx,
            m.op_index,
            total_lattice,
            phase,
            chi,
            contrib
        );
    }
    let w = sum / (anti_lg_indices.len() as f64);
    debug_log!("  direct anti stats: +={} -={} W={:.4}", n_plus, n_minus, w);
    Ok(w)
}

/// Spinor version of [`wigner_direct_anti_coset`]: directly iterates over
/// antiunitary little-group ops b ∈ M_k \ H_k instead of the a₀h construction.
///
/// For each antiunitary b, computes b² (guaranteed unitary in H_k by group
/// theory), looks up the spinor character via SU(2) composition, and sums to
/// get the Wigner indicator.
///
/// This avoids the a₀-selection and a₀h-composition issues that can cause
/// the main [`wigner_classify_spinor`] path to fail on Type III black-white
/// groups when (g₀R_h)² maps outside the little co-group lookup table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DirectAntiFailure {
    MissingSpinData,
    IndexOutOfRange,
    CharacterTableMismatch,
    SpinTableMismatch,
    SettingTransformFailed,
    SquareNotInSpinTable,
    SquareOutsideLittleGroup,
    SquareTranslationMismatch,
    AntiunitarySpinLookup,
    AntiunitarySu2Missing,
    SquareSu2Missing,
    Su2LiftMismatch,
    NonQuantized,
}

impl DirectAntiFailure {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::MissingSpinData => "missing_spin_data",
            Self::IndexOutOfRange => "index_out_of_range",
            Self::CharacterTableMismatch => "character_table_mismatch",
            Self::SpinTableMismatch => "spin_table_mismatch",
            Self::SettingTransformFailed => "setting_transform_failed",
            Self::SquareNotInSpinTable => "square_not_in_spin_table",
            Self::SquareOutsideLittleGroup => "square_outside_little_group",
            Self::SquareTranslationMismatch => "square_translation_mismatch",
            Self::AntiunitarySpinLookup => "antiunitary_spin_lookup",
            Self::AntiunitarySu2Missing => "antiunitary_su2_missing",
            Self::SquareSu2Missing => "square_su2_missing",
            Self::Su2LiftMismatch => "su2_lift_mismatch",
            Self::NonQuantized => "non_quantized",
        }
    }
}

/// How b² was matched to the H spin table.
/// 0=Exact (Seitz mod Z³), 1=Centering, 2=RotationOnly
pub type SpinMatchKind = u8;
pub const MATCH_EXACT: SpinMatchKind = 0;
pub const MATCH_CENTERING: SpinMatchKind = 1;
pub const MATCH_ROTATION_ONLY: SpinMatchKind = 2;

/// Per-term trace record for diagnosing non-quantized Wigner sums.
///
/// Each antiunitary coset element b ∈ a₀H produces one term in the Wigner sum:
/// χ((a₀b)²) with SU(2) central parity and Bloch phase.
#[derive(Debug, Clone)]
pub struct PerTermTrace {
    /// Index of this b in `mag_seitz`
    pub b_idx: usize,
    /// Rotation of b (in MSG frame)
    pub b_rot: Mat3I,
    /// Translation of b (in MSG frame)
    pub b_trans: [f64; 3],
    /// Rotation of b² = (a₀b)² (in Hall frame after setting transform)
    pub sq_rot: Mat3I,
    /// Translation of b² (reduced to [0,1))
    pub sq_trans: [f64; 3],
    /// Raw b² translation before reduction to [0,1)
    pub sq_trans_raw: [f64; 3],
    /// Translation of the matched H spin table entry
    pub sq_spin_trans: [f64; 3],
    /// Lattice difference: sq_raw - sq_spin_trans (rounded)
    pub trans_delta: [f64; 3],
    /// How the match was obtained
    pub match_kind: SpinMatchKind,
    /// Index of b² rotation in H spin table
    pub sq_spin_idx: usize,
    /// Local index in spin_lg_op_indices (for character lookup)
    pub sq_local_idx: usize,
    /// Spin character χ₀ at sq_local_idx (real part)
    pub chi0_re: f64,
    /// Spin character χ₀ at sq_local_idx (imag part, 0 if real-only)
    pub chi0_im: f64,
    /// true = central element (U_b² = -U_{b²}) → character sign flip
    pub central: bool,
    /// Bloch phase (real part)
    pub phase_re: f64,
    /// Bloch phase (imag part)
    pub phase_im: f64,
    /// Final term contribution = (central ? -χ₀ : χ₀) * exp(iφ) (real part)
    pub contrib_re: f64,
    /// Final term contribution (imag part)
    pub contrib_im: f64,
    /// Total translation used for Bloch phase (lattice + spin match shift)
    pub bloch_total_trans: [f64; 3],
    /// SU(2) lift of b: U_b = [a, b, c, d]
    pub u_b: [f64; 4],
    /// -U_b² (after neg_pauli, the actual value used for central parity)
    pub u_b_sq_actual: [f64; 4],
    /// SU(2) lift of b² in G spin table
    pub u_sq_g: [f64; 4],
}

static HALL_TO_SPIN_ORIGINS: OnceLock<[[f64; 3]; 231]> = OnceLock::new();
static HALL_TRANSLATION_LATTICES: OnceLock<Vec<Vec<[f64; 3]>>> = OnceLock::new();
const SEITZ_TRANS_TOL: f64 = 1e-5;

fn hall_to_spin_origin_for_sg(sg: u8) -> [f64; 3] {
    HALL_TO_SPIN_ORIGINS.get_or_init(build_hall_to_spin_origins)[sg as usize]
}

fn build_hall_to_spin_origins() -> [[f64; 3]; 231] {
    let mut origins = [[0.0; 3]; 231];
    for sg in 1u8..=230 {
        origins[sg as usize] = solve_hall_to_spin_origin_for_sg(sg).unwrap_or([0.0; 3]);
    }
    origins
}

fn solve_hall_to_spin_origin_for_sg(sg: u8) -> Option<[f64; 3]> {
    let hall = *super::generated_data::SG_DATA_HALL.get(sg as usize)? as usize;
    if hall == 0 {
        return None;
    }

    let hall_ops = SymmetryOps::from_database(hall).ok()?;
    let (spin_rots, spin_trans, _) = super::types::IrrepRecord::spin_ops_for_sg(sg);
    let spin_seitz = build_spin_seitz(spin_rots, spin_trans);
    if hall_ops.is_empty() || spin_seitz.is_empty() {
        return None;
    }

    let mut pairs: Vec<(SeitzOp, SeitzOp)> = Vec::new();
    if spin_seitz.len() <= hall_ops.len()
        && spin_seitz.iter().enumerate().all(|(i, spin)| {
            hall_ops
                .operations
                .get(i)
                .is_some_and(|hall| hall.rotation == spin.rot)
        })
    {
        for (i, spin) in spin_seitz.iter().enumerate() {
            let hall = &hall_ops.operations[i];
            pairs.push((
                SeitzOp::new(hall.rotation, hall.translation, false),
                spin.clone(),
            ));
        }
    } else {
        let mut used = vec![false; hall_ops.len()];
        for spin in &spin_seitz {
            if let Some((idx, hall)) = hall_ops
                .operations
                .iter()
                .enumerate()
                .find(|(idx, hall)| !used[*idx] && hall.rotation == spin.rot)
            {
                used[idx] = true;
                pairs.push((
                    SeitzOp::new(hall.rotation, hall.translation, false),
                    spin.clone(),
                ));
            }
        }
    }

    if pairs.is_empty() {
        return None;
    }

    let centering_shifts = centering_shifts_for_sg(sg);
    let denom = 24.0;
    let mut best = None;
    let mut best_score = 0usize;

    for ix in 0..24 {
        for iy in 0..24 {
            for iz in 0..24 {
                let origin = [ix as f64 / denom, iy as f64 / denom, iz as f64 / denom];
                let score = pairs
                    .iter()
                    .filter(|(hall, spin)| {
                        let transformed = apply_spin_origin_shift(hall.rot, hall.trans, origin);
                        translation_delta_in_lattice(
                            &[
                                transformed[0] - spin.trans[0],
                                transformed[1] - spin.trans[1],
                                transformed[2] - spin.trans[2],
                            ],
                            centering_shifts,
                        )
                    })
                    .count();
                if score > best_score {
                    best_score = score;
                    best = Some(origin);
                    if score == pairs.len() {
                        return best;
                    }
                }
            }
        }
    }

    if best_score == pairs.len() {
        best
    } else {
        None
    }
}

fn build_hall_translation_lattices() -> Vec<Vec<[f64; 3]>> {
    let mut lattices = vec![Vec::new(); 231];
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    for (sg, lattice) in lattices.iter_mut().enumerate().skip(1) {
        let Some(&hall) = super::generated_data::SG_DATA_HALL.get(sg) else {
            continue;
        };
        if hall == 0 {
            continue;
        }
        let Ok(hall_ops) = SymmetryOps::from_database(hall as usize) else {
            continue;
        };

        let mut shifts = Vec::new();
        for op in hall_ops
            .operations
            .iter()
            .filter(|op| op.rotation == identity)
        {
            let shift = normalize_translation(op.translation);
            if !shifts
                .iter()
                .any(|existing| translations_equal_mod_one(existing, &shift))
            {
                shifts.push(shift);
            }
        }
        *lattice = shifts;
    }

    lattices
}

fn normalize_frac(value: f64) -> f64 {
    let mut x = value % 1.0;
    if x < 0.0 {
        x += 1.0;
    }
    if x.abs() < 1e-10 || (x - 1.0).abs() < 1e-10 {
        0.0
    } else {
        x
    }
}

fn normalize_translation(trans: [f64; 3]) -> [f64; 3] {
    [
        normalize_frac(trans[0]),
        normalize_frac(trans[1]),
        normalize_frac(trans[2]),
    ]
}

fn translations_equal_mod_one(a: &[f64; 3], b: &[f64; 3]) -> bool {
    (0..3).all(|k| {
        let d = normalize_frac(a[k] - b[k]);
        d < SEITZ_TRANS_TOL || (d - 1.0).abs() < SEITZ_TRANS_TOL
    })
}

fn apply_spin_origin_shift(rot: Mat3I, trans: [f64; 3], origin: [f64; 3]) -> [f64; 3] {
    let mut t = trans;
    for i in 0..3 {
        let d: f64 = (0..3)
            .map(|j| {
                let delta = if i == j { 1.0 } else { 0.0 };
                (delta - rot[i][j] as f64) * origin[j]
            })
            .sum();
        t[i] = normalize_frac(t[i] - d);
    }
    t
}

fn translation_delta_in_lattice(delta: &[f64; 3], centering_shifts: &[[f64; 3]]) -> bool {
    let mut shifts: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]];
    for &cs in centering_shifts {
        if cs.iter().all(|&x| x.abs() < 1e-12) {
            continue;
        }
        let n = shifts.len();
        let cs_norm = [
            ((cs[0] % 1.0) + 1.0) % 1.0,
            ((cs[1] % 1.0) + 1.0) % 1.0,
            ((cs[2] % 1.0) + 1.0) % 1.0,
        ];
        for i in 0..n {
            shifts.push([
                (shifts[i][0] + cs_norm[0]) % 1.0,
                (shifts[i][1] + cs_norm[1]) % 1.0,
                (shifts[i][2] + cs_norm[2]) % 1.0,
            ]);
        }
    }

    let frac = [
        ((delta[0] % 1.0) + 1.0) % 1.0,
        ((delta[1] % 1.0) + 1.0) % 1.0,
        ((delta[2] % 1.0) + 1.0) % 1.0,
    ];
    shifts.iter().any(|s| {
        (0..3).all(|k| {
            let d = (frac[k] - s[k]).abs();
            d < SEITZ_TRANS_TOL || (d - 1.0).abs() < SEITZ_TRANS_TOL
        })
    })
}

pub fn wigner_classify_spinor_direct_anti(
    ctx: &SpinLiftContext,
    input: SpinorWignerInput<'_>,
    anti_lg_indices: &[usize],
    mag_seitz: &[SeitzOp],
    setting_xf: Option<&SettingTransform>,
) -> Result<CorepType, WignerClassificationError> {
    wigner_classify_spinor_direct_anti_diagnostic(
        ctx,
        input,
        anti_lg_indices,
        mag_seitz,
        setting_xf,
        None,
    )
    .map_err(|e| {
        WignerClassificationError::with_value(format!("direct anti path failed: {:?}", e), 0.0)
    })
}

pub fn wigner_classify_spinor_direct_anti_diagnostic(
    ctx: &SpinLiftContext,
    input: SpinorWignerInput<'_>,
    anti_lg_indices: &[usize],
    mag_seitz: &[SeitzOp],
    setting_xf: Option<&SettingTransform>,
    // Optional per-term trace collector (only filled when diagnosing failures)
    mut trace: Option<&mut Vec<PerTermTrace>>,
) -> Result<CorepType, DirectAntiFailure> {
    let SpinorWignerInput {
        characters_real: spin_chars_real,
        characters_imag: spin_chars_imag,
        operation_indices: spin_lg_op_indices,
        k_vector,
    } = input;
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    let (h_spin_rots, h_spin_trans, h_spin_su2) = ctx.h;
    let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;

    if h_spin_rots.is_empty() || h_spin_su2.is_empty() || anti_lg_indices.is_empty() {
        return Err(DirectAntiFailure::MissingSpinData);
    }
    if h_spin_rots.len() % 9 != 0
        || h_spin_trans.len() != h_spin_rots.len() / 3
        || h_spin_su2.len() != (h_spin_rots.len() / 9) * 4
        || g_spin_rots.len() % 9 != 0
        || g_spin_trans.len() != g_spin_rots.len() / 3
        || g_spin_su2.len() != (g_spin_rots.len() / 9) * 4
    {
        return Err(DirectAntiFailure::SpinTableMismatch);
    }
    if spin_chars_real.len() < spin_lg_op_indices.len() {
        return Err(DirectAntiFailure::CharacterTableMismatch);
    }
    if anti_lg_indices
        .iter()
        .any(|&index| index >= mag_seitz.len())
    {
        return Err(DirectAntiFailure::IndexOutOfRange);
    }

    let h_spin_seitz = build_spin_seitz(h_spin_rots, h_spin_trans);
    let g_spin_seitz = build_spin_seitz(g_spin_rots, g_spin_trans);
    if h_spin_seitz.is_empty() {
        return Err(DirectAntiFailure::MissingSpinData);
    }
    if spin_lg_op_indices
        .iter()
        .any(|&index| index as usize >= h_spin_seitz.len())
    {
        return Err(DirectAntiFailure::IndexOutOfRange);
    }

    let global_to_local: std::collections::HashMap<usize, usize> = spin_lg_op_indices
        .iter()
        .enumerate()
        .map(|(l, &g)| (g as usize, l))
        .collect();

    // Data-Hall → spin-table origin shift. Spin ops are rotation-reordered
    // into `SG_DATA_HALL[sg]` order but retain Bilbao/spin translations, so
    // the origin has to be solved from Hall ops vs spin ops, not taken from
    // `isotropy_origin` subgroup records.
    let origin = hall_to_spin_origin_for_sg(ctx.sg);
    let to_bilbao =
        |rot: Mat3I, trans: [f64; 3]| -> [f64; 3] { apply_spin_origin_shift(rot, trans, origin) };

    let n_anti = anti_lg_indices.len();
    let mut w_sum = Complex64::ZERO;

    // Centering shifts: fractional translations that together with Z³
    // span the full translation lattice (Body/I, Face/F, C, A, R, or empty for P).
    let centering_shifts = centering_shifts_for_sg(ctx.sg);

    for &b_idx in anti_lg_indices {
        let b = &mag_seitz[b_idx];

        // Apply setting transform if present: x_hall = T·x_msg + s.
        // This corrects for MSG-embedding basis ≠ canonical Hall basis.
        let (b_rot, b_trans) = if let Some(xf) = setting_xf {
            // SAFETY: the caller (via filter_little_group_with_transform)
            // validated atomically that the setting transform succeeds for
            // ALL MSG ops, so it must succeed for this antiunitary op too.
            let rot = xf
                .transform_rotation(&b.rot)
                .ok_or(DirectAntiFailure::SettingTransformFailed)?;
            let trans = xf
                .transform_translation(&b.rot, &b.trans)
                .ok_or(DirectAntiFailure::SettingTransformFailed)?;
            (rot, trans)
        } else {
            (b.rot, b.trans)
        };

        // Convert b to Bilbao setting, then square.
        let b_bilbao = SeitzOp::new(b_rot, to_bilbao(b_rot, b_trans), false);
        let (sq, lattice_sq) = square_seitz(&b_bilbao);

        // b² ∈ H₀ by group theory: b ∈ M_k ⇒ b² ∈ H_k ⇒ R_{b²} ∈ H₀.
        // Use LG-first matching to avoid picking a non-LG candidate.
        let (sq_spin_idx, sq_in_lg, match_kind) =
            match find_sq_spin_lg_first(&sq, &h_spin_seitz, spin_lg_op_indices, centering_shifts) {
                Some(v) => v,
                None => {
                    // Detailed diagnostic for primitive-group translation failures
                    {
                        static CNT: AtomicUsize = AtomicUsize::new(0);
                        let n = CNT.fetch_add(1, Ordering::Relaxed);
                        if n < 5 {
                            eprintln!("=== SQ_FAIL_DETAIL #{n}: SG{} ===", ctx.sg);
                            eprintln!(
                                "  b(before to_bilbao): rot={:?} trans=[{:.6},{:.6},{:.6}]",
                                b_rot, b_trans[0], b_trans[1], b_trans[2]
                            );
                            eprintln!(
                                "  b_bilbao: rot={:?} trans=[{:.6},{:.6},{:.6}]",
                                b_bilbao.rot,
                                b_bilbao.trans[0],
                                b_bilbao.trans[1],
                                b_bilbao.trans[2]
                            );
                            eprintln!(
                                "  sq: rot={:?} trans=[{:.6},{:.6},{:.6}]",
                                sq.rot, sq.trans[0], sq.trans[1], sq.trans[2]
                            );
                            // Show matching spin ops
                            eprintln!("  Spin ops with matching rot:");
                            let lg_ix: Vec<usize> =
                                spin_lg_op_indices.iter().map(|&x| x as usize).collect();
                            for (si, sop) in h_spin_seitz.iter().enumerate() {
                                if sop.rot == sq.rot {
                                    let in_lg = lg_ix.contains(&si);
                                    eprintln!(
                                        "    [{}] trans=[{:.6},{:.6},{:.6}] in_lg={}",
                                        si, sop.trans[0], sop.trans[1], sop.trans[2], in_lg
                                    );
                                }
                            }
                            eprintln!("  centering_shifts={:.6?}", centering_shifts);
                            eprintln!("  sg_setting_origin={:.6?}", origin);
                            // Show how to_bilbao transforms b
                            let tb = to_bilbao(b_rot, b_trans);
                            eprintln!("  to_bilbao(b): [{:.6},{:.6},{:.6}]", tb[0], tb[1], tb[2]);
                        }
                    }
                    debug_log!(
                        "  SPINOR_DIRECT_ANTI fail: b[{}]² rot={:?} not in H spin ops",
                        b_idx,
                        sq.rot
                    );
                    return Err(DirectAntiFailure::SquareNotInSpinTable);
                }
            };

        let sq_local_idx = if sq_in_lg {
            *global_to_local
                .get(&sq_spin_idx)
                .ok_or(DirectAntiFailure::SquareOutsideLittleGroup)?
        } else {
            debug_log!(
                "  SPINOR_DIRECT_ANTI fail: b[{}]² spin[{}] not in LG idxs",
                b_idx,
                sq_spin_idx
            );
            return Err(DirectAntiFailure::SquareOutsideLittleGroup);
        };
        let sq_spin = h_spin_seitz
            .get(sq_spin_idx)
            .ok_or(DirectAntiFailure::SquareNotInSpinTable)?;
        let spin_match_shift = [
            sq.trans[0] - sq_spin.trans[0],
            sq.trans[1] - sq_spin.trans[1],
            sq.trans[2] - sq_spin.trans[2],
        ];

        // SU(2) lift of b (rotation-only lookup in G spin ops — b may have
        // improper rotation from G \ H for Type III).
        //
        // IMPORTANT: G spin table is in the PARENT G's coordinate frame.
        // The MSG embedding is also in G's frame (MSG ops come from
        // get_spacegroup_operations with parent's Hall number).
        // Therefore G lookup must use the UNTRANSFORMED b.rot, NOT the
        // H-Hall-transformed b_rot.  The setting transform maps MSG→H,
        // but G spin is not in H's frame.
        let b_spin_idx = g_spin_seitz
            .iter()
            .position(|s| s.rot == b.rot)
            .or_else(|| {
                let neg: Mat3I = [
                    [-b.rot[0][0], -b.rot[0][1], -b.rot[0][2]],
                    [-b.rot[1][0], -b.rot[1][1], -b.rot[1][2]],
                    [-b.rot[2][0], -b.rot[2][1], -b.rot[2][2]],
                ];
                g_spin_seitz.iter().position(|s| s.rot == neg)
            });
        let Some(b_spin_idx) = b_spin_idx else {
            debug_log!(
                "  SPINOR_DIRECT_ANTI AntiunitarySpinLookup: H SG{} b_idx={} b.rot={:?} b_rot(H)={:?} G spin ops={} rotations={:?}",
                ctx.sg,
                b_idx,
                b.rot,
                b_rot,
                g_spin_seitz.len(),
                g_spin_seitz.iter().map(|op| op.rot).collect::<Vec<_>>()
            );
            return Err(DirectAntiFailure::AntiunitarySpinLookup);
        };
        let u_b =
            spin_su2_at(g_spin_su2, b_spin_idx).ok_or(DirectAntiFailure::AntiunitarySu2Missing)?;

        // SU(2) central detection: U_b² vs canonical U_{b²}.
        // Compute entirely in G frame to avoid cross-gauge comparison.
        // b.rot is in G/MSG frame (untransformed), and G spin table is
        // in G frame.  H spin table may use a different axis convention.
        // Antiunitary square for spin-1/2:
        //   b = Θ·g  →  D(b)² = -D(g)²  (since Θ² = -1)
        // Use neg_pauli(compose(U, U)) = -U² explicitly rather than
        // relying on the spin-table gauge to detect the Θ² sign.
        let u_b_sq = neg_pauli(&su2_compose(&u_b, &u_b));
        let b_sq_rot_ms_g = crate::mathfunc::mat_multiply_matrix_i3(&b.rot, &b.rot);
        let g_sq_idx = g_spin_seitz
            .iter()
            .position(|s| s.rot == b_sq_rot_ms_g)
            .or_else(|| {
                let neg: crate::mathfunc::Mat3I = [
                    [
                        -b_sq_rot_ms_g[0][0],
                        -b_sq_rot_ms_g[0][1],
                        -b_sq_rot_ms_g[0][2],
                    ],
                    [
                        -b_sq_rot_ms_g[1][0],
                        -b_sq_rot_ms_g[1][1],
                        -b_sq_rot_ms_g[1][2],
                    ],
                    [
                        -b_sq_rot_ms_g[2][0],
                        -b_sq_rot_ms_g[2][1],
                        -b_sq_rot_ms_g[2][2],
                    ],
                ];
                g_spin_seitz.iter().position(|s| s.rot == neg)
            });
        let u_sq_g = match g_sq_idx {
            Some(idx) => spin_su2_at(g_spin_su2, idx).ok_or(DirectAntiFailure::SquareSu2Missing)?,
            None => return Err(DirectAntiFailure::SquareNotInSpinTable),
        };
        // U_b² vs canonical U_{b²} in the G spin table.
        // Since u_b_sq = -U_b² (neg_pauli at :1842 accounts for Θ²=-1):
        //   LiftRelation::Same: -U_b² = +U_{b²} → U_b² = -U_{b²}
        //     → spatial square IS Ē → Θ² already compensated → central=false
        //   LiftRelation::EBar: -U_b² = -U_{b²} → U_b² = +U_{b²}
        //     → spatial square is NOT Ē → Θ² contributes alone → central=true
        let spatial_central =
            su2_lift_relation(&u_b_sq, &u_sq_g).ok_or(DirectAntiFailure::Su2LiftMismatch)?;
        let mut central = spatial_central == LiftRelation::EBar;

        // ── G→H spin frame parity ────────────────────────────────────────
        // For G≠H with signed-permutation setting transform, apply the
        // axial vector transform Q = det(P)·P via spin table comparison.
        if let Some(xf) = setting_xf {
            let rows_ok = xf.basis.iter().all(|row| {
                let nonzeros = row.iter().filter(|&&v| v.abs() > 0.1).count();
                nonzeros == 1
                    && row.iter().all(|&v| {
                        let r = v.round();
                        r == -1.0 || r == 0.0 || r == 1.0
                    })
            });
            let cols_ok =
                (0..3).all(|j| (0..3).filter(|&i| xf.basis[i][j].abs() > 0.1).count() == 1);
            let mut ppt = [[0i32; 3]; 3];
            let p_i32: [[i32; 3]; 3] = xf.basis.map(|row| row.map(|v| v.round() as i32));
            for i in 0..3 {
                for j in 0..3 {
                    ppt[i][j] = (0..3).map(|k| p_i32[i][k] * p_i32[j][k]).sum();
                }
            }
            if rows_ok && cols_ok && ppt == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                let det_p = p_i32[0][0] * (p_i32[1][1] * p_i32[2][2] - p_i32[1][2] * p_i32[2][1])
                    - p_i32[0][1] * (p_i32[1][0] * p_i32[2][2] - p_i32[1][2] * p_i32[2][0])
                    + p_i32[0][2] * (p_i32[1][0] * p_i32[2][1] - p_i32[1][1] * p_i32[2][0]);
                let q: [[i32; 3]; 3] = if det_p == -1 {
                    [
                        [-p_i32[0][0], -p_i32[0][1], -p_i32[0][2]],
                        [-p_i32[1][0], -p_i32[1][1], -p_i32[1][2]],
                        [-p_i32[2][0], -p_i32[2][1], -p_i32[2][2]],
                    ]
                } else {
                    p_i32
                };
                if q != [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                    // Use b_sq_rot_ms_g (G frame) NOT sq.rot (H frame).
                    if let Some(parity) = compute_signed_perm_spin_parity(
                        &q,
                        &b_sq_rot_ms_g,
                        g_spin_rots,
                        g_spin_su2,
                        h_spin_rots,
                        h_spin_su2,
                    ) && (parity - (-1.0)).abs() < 0.1
                    {
                        central = !central;
                    }
                }
            }
        }

        // Bloch phase includes both square reduction and the shift required
        // to match the canonical spin-table representative.
        let total_translation = [
            lattice_sq[0] as f64 + spin_match_shift[0],
            lattice_sq[1] as f64 + spin_match_shift[1],
            lattice_sq[2] as f64 + spin_match_shift[2],
        ];
        let phase = bloch_phase_f64(kx, ky, kz, kd, &total_translation);

        let chi0 = Complex64::new(
            *spin_chars_real
                .get(sq_local_idx)
                .ok_or(DirectAntiFailure::CharacterTableMismatch)?,
            spin_chars_imag.get(sq_local_idx).copied().unwrap_or(0.0),
        );
        let chi = if central { -chi0 } else { chi0 };
        let contrib = chi * phase;

        // Record per-term trace if the caller requested diagnostics
        if let Some(t) = trace.as_mut() {
            t.push(PerTermTrace {
                b_idx,
                b_rot: b.rot,
                b_trans: b.trans,
                sq_rot: sq.rot,
                sq_trans: sq.trans,
                sq_trans_raw: [
                    sq.trans[0] + lattice_sq[0] as f64,
                    sq.trans[1] + lattice_sq[1] as f64,
                    sq.trans[2] + lattice_sq[2] as f64,
                ],
                sq_spin_trans: sq_spin.trans,
                trans_delta: [
                    sq.trans[0] + lattice_sq[0] as f64 - sq_spin.trans[0],
                    sq.trans[1] + lattice_sq[1] as f64 - sq_spin.trans[1],
                    sq.trans[2] + lattice_sq[2] as f64 - sq_spin.trans[2],
                ],
                match_kind,
                sq_spin_idx,
                sq_local_idx,
                chi0_re: chi0.re,
                chi0_im: chi0.im,
                central,
                phase_re: phase.re,
                phase_im: phase.im,
                contrib_re: contrib.re,
                contrib_im: contrib.im,
                bloch_total_trans: total_translation,
                u_b,
                u_b_sq_actual: u_b_sq,
                u_sq_g,
            });
        }

        w_sum += contrib;
    }

    let w = w_sum / (n_anti as f64);

    // Dimension from identity canonical lift.
    #[cfg(feature = "debug-corep")]
    let id_rot: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    #[cfg(feature = "debug-corep")]
    let u_id = [1.0, 0.0, 0.0, 0.0];
    #[cfg(feature = "debug-corep")]
    let _h_dim = spin_lg_op_indices
        .iter()
        .map(|&x| x as usize)
        .find_map(|si| {
            let sop = h_spin_seitz.get(si)?;
            if sop.rot != id_rot {
                return None;
            }
            let u = spin_su2_at(h_spin_su2, si)?;
            if su2_same_up_to_sign(&u, &u_id) != Some(false) {
                return None;
            }
            let local = *global_to_local.get(&si)?;
            spin_chars_real
                .get(local)
                .map(|&c| c.abs().round().max(1.0))
        })
        .unwrap_or_else(|| {
            spin_chars_real
                .first()
                .map(|&c| c.abs().round().max(1.0))
                .unwrap_or(1.0)
        });

    // Tolerance for Wigner sum quantization checks.
    // 1e-5 is tight enough to reject true non-quantized results (W≥0.167)
    // while accommodating floating-point accumulation from irrational
    // character values (e.g. √3/2 ≈ 0.866) summed over ≤24 terms.
    // Wigner criterion: the normalised sum W = (1/|H|) Σ χ((a₀h)²)
    // must be 0 (type C), +1 (type A), or -1 (type B) for any irrep.
    // The ±dim branches are INCORRECT — W is independent of dimension.
    // See e.g. arXiv:2211.10740.
    let tol = 1e-5;
    if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        Ok(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        Ok(CorepType::B)
    } else if w.norm() < tol {
        Ok(CorepType::C)
    } else {
        Err(DirectAntiFailure::NonQuantized)
    }
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
        if tr[i] < 0.0 {
            tr[i] += 1.0;
            l[i] -= 1;
        }
    }
    (tr, l)
}

/// Build the strict index map from final-Hall `h_seitz` order to PIR operation order.
///
/// This is the only scalar character-operation pairing helper. Both operation
/// universes must be complete and equal-sized; rotations must match exactly
/// and fractional translations must match directly within
/// [`SEITZ_TRANS_TOL`]. Integer-lattice shifts and rotation-only matches are
/// intentionally rejected because they can carry a different Bloch phase.
pub fn build_h_to_irrep_op_map(
    h_seitz: &[SeitzOp],
    irrep_rots: &[i32],
    irrep_trans: &[f64],
) -> Option<Vec<usize>> {
    let n_ops = h_seitz.len();
    if n_ops == 0 || irrep_rots.len() % 9 != 0 {
        return None;
    }
    let n_ir_ops = irrep_rots.len() / 9;
    let expected_trans_len = n_ir_ops.checked_mul(3)?;
    if n_ir_ops == 0 || n_ops != n_ir_ops || irrep_trans.len() != expected_trans_len {
        return None;
    }
    let mut used = vec![false; n_ir_ops];
    let mut map = Vec::with_capacity(n_ops);
    for h in h_seitz {
        let mut found = None;
        for ir_idx in 0..n_ir_ops {
            let roff = ir_idx * 9;
            let toff = ir_idx * 3;
            let rotation_matches = irrep_rots[roff..roff + 9]
                == [
                    h.rot[0][0],
                    h.rot[0][1],
                    h.rot[0][2],
                    h.rot[1][0],
                    h.rot[1][1],
                    h.rot[1][2],
                    h.rot[2][0],
                    h.rot[2][1],
                    h.rot[2][2],
                ];
            let translation_matches = (0..3)
                .all(|axis| (h.trans[axis] - irrep_trans[toff + axis]).abs() < SEITZ_TRANS_TOL);
            if rotation_matches && translation_matches {
                if found.is_some() || used[ir_idx] {
                    return None;
                }
                found = Some(ir_idx);
            }
        }
        let ir_idx = found?;
        used[ir_idx] = true;
        map.push(ir_idx);
    }
    if used.iter().all(|matched| *matched) {
        Some(map)
    } else {
        None
    }
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

/// Pure translations from the Hall setting used by the generated irrep data.
///
/// The centering type is not enough here: the selected Hall setting may permute
/// conventional axes, so e.g. an A-centered ITA number can appear as a
/// `(1/2, 1/2, 0)` pure translation in the data-Hall coordinates. Derive the
/// finite translation lattice from the actual Hall operations instead.
fn centering_shifts_for_sg(sg: u8) -> &'static [[f64; 3]] {
    HALL_TRANSLATION_LATTICES
        .get_or_init(build_hall_translation_lattices)
        .get(sg as usize)
        .map(Vec::as_slice)
        .unwrap_or(&[])
}

/// Find the spin op index for a computed Seitz square `sq`, preferring
/// candidates **inside** `spin_lg_op_indices` over the full database.
///
/// Priority:
/// 1. Full Seitz match (rotation + translation mod lattice) inside LG candidates
/// 2. Centering-equivalent: translation differs by a lattice vector
///    (including centered-cell vectors like 1/2,1/2,1/2).
///
/// Rotation-only fallback is INTENTIONALLY REMOVED per codex review:
/// same rotation with different translations are different group elements.
pub(crate) fn find_sq_spin_lg_first(
    sq: &SeitzOp,
    h_spin_seitz: &[SeitzOp],
    spin_lg_op_indices: &[u16],
    centering_shifts: &[[f64; 3]],
) -> Option<(usize, bool, SpinMatchKind)> {
    let lg_cands: Vec<usize> = spin_lg_op_indices.iter().map(|&x| x as usize).collect();

    // 1. Full Seitz match (translation mod Z³) inside LG
    for &si in &lg_cands {
        if let Some(sop) = h_spin_seitz.get(si)
            && sop.rot == sq.rot
        {
            let delta = [
                sq.trans[0] - sop.trans[0],
                sq.trans[1] - sop.trans[1],
                sq.trans[2] - sop.trans[2],
            ];
            if translation_delta_in_lattice(&delta, centering_shifts) {
                return Some((si, true, MATCH_EXACT));
            }
        }
    }

    // 2. No match.  Log first few failures for diagnosis.
    {
        static CNT: AtomicUsize = AtomicUsize::new(0);
        let n = CNT.fetch_add(1, Ordering::Relaxed);
        if n < 20 {
            let best = lg_cands
                .iter()
                .filter_map(|&si| h_spin_seitz.get(si))
                .filter(|s| s.rot == sq.rot)
                .map(|s| {
                    let d = [
                        sq.trans[0] - s.trans[0],
                        sq.trans[1] - s.trans[1],
                        sq.trans[2] - s.trans[2],
                    ];
                    (
                        s.trans,
                        d,
                        (d[0] - d[0].round()).abs()
                            + (d[1] - d[1].round()).abs()
                            + (d[2] - d[2].round()).abs(),
                    )
                })
                .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
            if let Some((sp_trans, delta, err)) = best {
                eprintln!(
                    "  SQ_FAIL #{n}: sq_rot={:?} sq_t=[{:.12}, {:.12}, {:.12}] sp_t=[{:.12}, {:.12}, {:.12}] d=[{:.12}, {:.12}, {:.12}] err={:.2e}",
                    sq.rot,
                    sq.trans[0],
                    sq.trans[1],
                    sq.trans[2],
                    sp_trans[0],
                    sp_trans[1],
                    sp_trans[2],
                    delta[0],
                    delta[1],
                    delta[2],
                    err
                );
            } else {
                eprintln!(
                    "  SQ_FAIL #{n}: sq_rot={:?} sq_t=[{:.12}, {:.12}, {:.12}] NO rotation match in LG (n_lg={})",
                    sq.rot,
                    sq.trans[0],
                    sq.trans[1],
                    sq.trans[2],
                    lg_cands.len()
                );
            }
        }
    }
    None
}

// ── Spinor (double-group) operations ───────────────────────────────────────
//
// Bilbao spin.dat SU(2) convention (verified by scripts/test_su2_closure.py,
// 229/229 SGs pass at 100% closure):
//
//   rot[9] trans[3] amp[4] phase[4]/π
//   U_ij = amp[ij] · exp(iπ · phase[ij])
//
// Converted at generation time to real Pauli coefficients [u₀,u₁,u₂,u₃]:
//
//   U = u₀·I + i(u₁·σx + u₂·σy + u₃·σz)
//     = [[u₀ + iu₃,    u₂ + iu₁],
//        [-u₂ + iu₁,    u₀ - iu₃]]
//
// For crystallographic point groups the coefficients take values only
// from {0, ±½, ±1/√2, ±√3/2, ±1} and are stored as exact f64.
//
// Composition follows quaternion multiplication (isomorphic to SU(2)):
//   (u₀,u)·(v₀,v) = (u₀v₀ − u·v,  u₀v + v₀u + u×v)

/// Compose two SU(2) matrices in Pauli coefficient representation.
pub fn su2_compose(a: &[f64; 4], b: &[f64; 4]) -> [f64; 4] {
    let [u0, u1, u2, u3] = *a;
    let [v0, v1, v2, v3] = *b;
    // Quaternion multiply: (u₀, u₁, u₂, u₃) · (v₀, v₁, v₂, v₃)
    [
        u0 * v0 - u1 * v1 - u2 * v2 - u3 * v3,
        u0 * v1 + u1 * v0 + u2 * v3 - u3 * v2,
        u0 * v2 - u1 * v3 + u2 * v0 + u3 * v1,
        u0 * v3 + u1 * v2 - u2 * v1 + u3 * v0,
    ]
}

/// Check if two SU(2) Pauli coefficient vectors match up to sign (central element Ē).
///
/// Both `a` and `b` are `[u₀, u₁, u₂, u₃]` Pauli coefficients. The central
/// element of SU(2) is Ē = -I = [-1, 0, 0, 0], so `a = ±b` iff the unit-vector
/// dot product is ±1.
///
/// Returns `Some(false)` if a ≈ b, `Some(true)` if a ≈ -b (central differs),
/// `None` if they are unrelated.
/// Relationship between two SU(2) lifts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LiftRelation {
    /// U_a = +U_b (same lift, no central element Ē)
    Same,
    /// U_a = -U_b (opposite lift, Ē = nontrivial central element)
    EBar,
}

/// Compare two SU(2) lifts.  Returns `None` if the lifts are unrelated
/// (cos ≠ ±1), or `Some(LiftRelation)` indicating whether they are the
/// same lift or differ by the central element Ē = -I.
pub fn su2_lift_relation(a: &[f64; 4], b: &[f64; 4]) -> Option<LiftRelation> {
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    let na = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3]).sqrt();
    let nb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + b[3] * b[3]).sqrt();
    if na < 1e-10 || nb < 1e-10 {
        return None;
    }
    let cos = dot / (na * nb);
    if (cos - 1.0).abs() < 1e-6 {
        Some(LiftRelation::Same)
    } else if (cos + 1.0).abs() < 1e-6 {
        Some(LiftRelation::EBar)
    } else {
        None
    }
}

/// Legacy wrapper.  Prefer [`su2_lift_relation`].
/// Returns `Some(false)` for [`LiftRelation::Same`] and `Some(true)` for [`LiftRelation::EBar`].
pub fn su2_same_up_to_sign(a: &[f64; 4], b: &[f64; 4]) -> Option<bool> {
    su2_lift_relation(a, b).map(|r| match r {
        LiftRelation::Same => false,
        LiftRelation::EBar => true,
    })
}

// ── Spinor (double-group) Wigner test ──────────────────────────────────────
//
// For spinor irreps, each spatial operation {R|t} has two lifts in the
// double group: g and Ēg where Ē = 2π rotation = −1.  The character of
// the double-group element (a₀h)² is:
//
//   χ((a₀h)²) = {  χ(g_k)  if U_sq ≈ U_k
//                { −χ(g_k)  if U_sq ≈ −U_k  (central element Ē appears)
//
// where U_sq = (U_{a₀} U_h)² computed via Pauli-coefficient quaternion
// multiplication, and U_k is the canonical SU(2) lift of the spatial
// operation (a₀h)².

/// Extract a Seitz operation from the spin-op flat arrays.
fn spin_seitz_at(idx: usize, spin_op_rots: &[i32], spin_op_trans: &[f64]) -> Option<SeitzOp> {
    if 9 * idx + 8 >= spin_op_rots.len() || 3 * idx + 2 >= spin_op_trans.len() {
        return None;
    }
    let r = [
        [
            spin_op_rots[9 * idx],
            spin_op_rots[9 * idx + 1],
            spin_op_rots[9 * idx + 2],
        ],
        [
            spin_op_rots[9 * idx + 3],
            spin_op_rots[9 * idx + 4],
            spin_op_rots[9 * idx + 5],
        ],
        [
            spin_op_rots[9 * idx + 6],
            spin_op_rots[9 * idx + 7],
            spin_op_rots[9 * idx + 8],
        ],
    ];
    let t = [
        spin_op_trans[3 * idx],
        spin_op_trans[3 * idx + 1],
        spin_op_trans[3 * idx + 2],
    ];
    Some(SeitzOp::new(r, t, false))
}

/// Extract Pauli coefficients [u₀,u₁,u₂,u₃] from the spin-op flat array.
pub fn spin_su2_at(spin_op_su2: &[f64], idx: usize) -> Option<[f64; 4]> {
    if 4 * idx + 3 >= spin_op_su2.len() {
        return None;
    }
    Some([
        spin_op_su2[4 * idx],
        spin_op_su2[4 * idx + 1],
        spin_op_su2[4 * idx + 2],
        spin_op_su2[4 * idx + 3],
    ])
}

/// Check if two Seitz ops match modulo lattice translations.
fn same_seitz_mod_lattice(a: &SeitzOp, b: &SeitzOp) -> bool {
    if a.rot != b.rot {
        return false;
    }
    for i in 0..3 {
        let d = a.trans[i] - b.trans[i];
        if (d - d.round()).abs() > SEITZ_TRANS_TOL {
            return false;
        }
    }
    true
}

/// Build the legacy mapping from H_ops (spglib order) to spin-op index.
///
/// Strategy:
/// 1. Full Seitz matching (rotation + translation via `same_seitz_mod_lattice`)
///    — disambiguates g vs Ēg lifts at BZ-boundary k-points.
/// 2. Fallback: rotation-only, but only when exactly one candidate exists in
///    `spin_lg_op_indices`.  Multiple candidates with the same rotation but
///    different translations are UNRESOLVED → return None (caller skips).
///
/// Uses a `Vec` (not HashSet) for deterministic iteration order.  This helper
/// remains for Wigner classification and diagnostics only; it must not be
/// used to pair final spinor character values with Seitz operations.  Final
/// output uses [`build_h_to_spin_map_exact`].
pub(crate) fn build_h_to_spin_map(
    h_seitz: &[SeitzOp],
    spin_seitz: &[SeitzOp],
    spin_lg_op_indices: &[u16],
) -> Vec<Option<usize>> {
    let allowed: Vec<usize> = spin_lg_op_indices.iter().map(|&x| x as usize).collect();

    let mut map = Vec::with_capacity(h_seitz.len());
    for h in h_seitz {
        // 1. Full Seitz matching
        let full = allowed.iter().copied().find(|&si| {
            spin_seitz
                .get(si)
                .is_some_and(|sop| same_seitz_mod_lattice(h, sop))
        });
        if let Some(si) = full {
            H2S_OK.fetch_add(1, Ordering::Relaxed);
            map.push(Some(si));
            continue;
        }

        // 2. Rotation-only fallback — only when uniquely resolved
        let rots: Vec<usize> = allowed
            .iter()
            .copied()
            .filter(|&si| spin_seitz.get(si).is_some_and(|sop| sop.rot == h.rot))
            .collect();

        if rots.len() == 1 {
            H2S_OK.fetch_add(1, Ordering::Relaxed);
            map.push(Some(rots[0]));
        } else if rots.len() > 1 {
            H2S_AMBIGUOUS.fetch_add(1, Ordering::Relaxed);
            map.push(None);
        } else {
            H2S_MISSING.fetch_add(1, Ordering::Relaxed);
            map.push(None);
        }
    }
    map
}

/// Build the operation map used by final spinor character output.
///
/// Unlike [`build_h_to_spin_map`], this deliberately accepts only the same
/// representative: rotations must be equal and every translation component
/// must agree within [`SEITZ_TRANS_TOL`].  Integer lattice shifts are not
/// equivalent here because they carry a Bloch phase.  Candidates are limited
/// to `spin_lg_op_indices`; a missing or ambiguous match is represented by
/// `None` and must be rejected by the caller before constructing characters.
pub(crate) fn build_h_to_spin_map_exact(
    h_seitz: &[SeitzOp],
    spin_seitz: &[SeitzOp],
    spin_lg_op_indices: &[u16],
) -> Vec<Option<usize>> {
    let allowed: Vec<usize> = spin_lg_op_indices
        .iter()
        .map(|&index| index as usize)
        .collect();

    h_seitz
        .iter()
        .map(|h| {
            let candidates: Vec<usize> = allowed
                .iter()
                .copied()
                .filter(|&spin_index| {
                    spin_seitz.get(spin_index).is_some_and(|spin| {
                        spin.rot == h.rot
                            && spin
                                .trans
                                .iter()
                                .zip(h.trans)
                                .all(|(spin_t, h_t)| (spin_t - h_t).abs() < SEITZ_TRANS_TOL)
                    })
                })
                .collect();
            match candidates.as_slice() {
                [candidate] => Some(*candidate),
                _ => None,
            }
        })
        .collect()
}

/// Build a Vec<SeitzOp> from the spin-op flat arrays (public for testing).
pub fn build_spin_seitz(spin_op_rots: &[i32], spin_op_trans: &[f64]) -> Vec<SeitzOp> {
    let n = (spin_op_rots.len() / 9).min(spin_op_trans.len() / 3);
    (0..n)
        .filter_map(|i| spin_seitz_at(i, spin_op_rots, spin_op_trans))
        .collect()
}

fn spin_table_is_well_formed(rotations: &[i32], translations: &[f64], su2: &[f64]) -> bool {
    rotations.len().is_multiple_of(9)
        && translations.len() == rotations.len() / 3
        && su2.len() == (rotations.len() / 9) * 4
}

/// **DIAGNOSTIC ONLY — not an authoritative Wigner test.**
///
/// Bilbao spin.dat may contain extra character-like values at some k-points.
/// These values are NOT guaranteed to be term-by-term Wigner summands
/// χ((a₀h)²).  Counterexample: for a spin-½ grey group with a₀ = Θ,
/// the h = E term must be χ(Θ²) = χ(Ē) = -1, yet the stored extra value
/// can be 0 (see SG3 A3 at k=(½,0,½)).
///
/// This function only checks whether the raw sum accidentally gives a
/// quantized value (0, ±1, or ±|H|).  It must NOT be used as the primary
/// spinor Wigner test — use [`wigner_classify_spinor`] instead.
pub fn diagnostic_imag_sum(imag: &[f64]) -> f64 {
    imag.iter().sum()
}

/// Spinor Wigner test evaluated in the magnetic-group coordinate setting.
///
/// The magnetic operation `a0` and every `h` must be composed in one common
/// setting.  The previous implementation mixed `a0` from the MSG database
/// with canonical `h` operations from the standalone H spin table.  That
/// happens to work when the subgroup embedding is standard, but fails for
/// non-trivial embeddings (most visibly for cubic C2/C3 combinations).
fn wigner_classify_spinor_msg_gauge(
    ctx: &SpinLiftContext,
    input: SpinorWignerInput<'_>,
    group: WignerGroupContext<'_>,
) -> Result<CorepType, WignerClassificationError> {
    let SpinorWignerInput {
        characters_real: spin_chars_real,
        characters_imag: spin_chars_imag,
        operation_indices: spin_lg_op_indices,
        k_vector,
    } = input;
    let n_lg_ops = spin_lg_op_indices.len();
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    let WignerGroupContext {
        unitary_indices: unitary_mag_indices,
        magnetic_ops: mag_seitz,
        unitary_ops: h_seitz,
        antiunitary_representative: a0_idx,
    } = group;
    let (h_spin_rots, h_spin_trans, h_spin_su2) = ctx.h;
    let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;
    let h_spin_seitz = build_spin_seitz(h_spin_rots, h_spin_trans);
    let g_spin_seitz = build_spin_seitz(g_spin_rots, g_spin_trans);
    if h_spin_seitz.is_empty() || g_spin_seitz.is_empty() {
        return Err(WignerClassificationError::new(
            "msg_gauge: missing spin data",
        ));
    }

    let h_to_spin = build_h_to_spin_map(h_seitz, &h_spin_seitz, spin_lg_op_indices);
    let global_to_local: std::collections::HashMap<usize, usize> = spin_lg_op_indices
        .iter()
        .enumerate()
        .map(|(local, &global)| (global as usize, local))
        .collect();

    let mut spin_to_mag = std::collections::HashMap::<usize, usize>::new();
    for &mag_idx in unitary_mag_indices {
        let h_match = find_seitz(&mag_seitz[mag_idx].rot, &mag_seitz[mag_idx].trans, h_seitz)
            .ok_or_else(|| {
                WignerClassificationError::new("msg_gauge: unitary mag op not found in H seitz")
            })?;
        if let Some(Some(spin_idx)) = h_to_spin.get(h_match.op_index) {
            spin_to_mag.entry(*spin_idx).or_insert(mag_idx);
        }
    }
    let has_unmapped = spin_lg_op_indices
        .iter()
        .any(|&idx| !spin_to_mag.contains_key(&(idx as usize)));
    if has_unmapped {
        MSG_GAUGE_MAP_FAIL.fetch_add(1, Ordering::Relaxed);
    }

    let a0_spatial = SeitzOp::new(mag_seitz[a0_idx].rot, mag_seitz[a0_idx].trans, false);
    let (a0_spin_idx, _) = find_spin_in_db(&a0_spatial, &g_spin_seitz)
        .ok_or_else(|| WignerClassificationError::new("msg_gauge: a0 not found in G spin"))?;
    let u_a0 = spin_su2_at(g_spin_su2, a0_spin_idx)
        .ok_or_else(|| WignerClassificationError::new("msg_gauge: a0 SU(2) lift missing"))?;
    let eta_ebar = -1.0;
    let mut w_sum = Complex64::ZERO;
    let mut n_mapped: usize = 0;

    for &global_index in spin_lg_op_indices.iter().take(n_lg_ops) {
        let global_spin_idx = global_index as usize;
        let mag_idx = match spin_to_mag.get(&global_spin_idx) {
            Some(&m) => m,
            None => continue,
        };
        let h_msg = SeitzOp::new(mag_seitz[mag_idx].rot, mag_seitz[mag_idx].trans, false);

        let (g0h, l1) = compose_seitz(&a0_spatial, &h_msg);
        let (sq, lattice_sq) = square_seitz(&g0h);
        let sq_h_match = find_seitz(&sq.rot, &sq.trans, h_seitz).ok_or_else(|| {
            WignerClassificationError::new("msg_gauge: square not found in H seitz")
        })?;
        let sq_spin_idx = h_to_spin
            .get(sq_h_match.op_index)
            .copied()
            .flatten()
            .ok_or_else(|| WignerClassificationError::new("msg_gauge: square not in H→spin map"))?;
        let sq_local_idx = *global_to_local.get(&sq_spin_idx).ok_or_else(|| {
            WignerClassificationError::new("msg_gauge: spin→local lookup missing")
        })?;

        let (h_g_idx, _) = find_spin_in_db(&h_msg, &g_spin_seitz)
            .ok_or_else(|| WignerClassificationError::new("msg_gauge: h not found in G spin"))?;
        let u_h_g = spin_su2_at(g_spin_su2, h_g_idx)
            .ok_or_else(|| WignerClassificationError::new("msg_gauge: h SU(2) lift missing"))?;
        let u_g0h = su2_compose(&u_a0, &u_h_g);
        let u_sq_spatial = su2_compose(&u_g0h, &u_g0h);
        let u_sq_h = spin_su2_at(h_spin_su2, sq_spin_idx)
            .ok_or_else(|| WignerClassificationError::new("msg_gauge: sq SU(2) lift missing"))?;
        let spatial_central = su2_same_up_to_sign(&u_sq_spatial, &u_sq_h).ok_or_else(|| {
            WignerClassificationError::new("msg_gauge: SU(2) sign comparison failed")
        })?;

        let central = !spatial_central;

        let r_l1 = mat_vec_i32(&g0h.rot, &l1);
        let total_lattice = add3(
            &add3(&lattice_sq, &sq_h_match.lattice_shift),
            &add3(&l1, &r_l1),
        );
        let phase = bloch_phase(kx, ky, kz, kd, &total_lattice);
        let chi0 = Complex64::new(
            spin_chars_real[sq_local_idx],
            spin_chars_imag.get(sq_local_idx).copied().unwrap_or(0.0),
        );
        let chi = if central { eta_ebar * chi0 } else { chi0 };
        w_sum += chi * phase;
        n_mapped += 1;
    }

    if n_mapped == 0 {
        return Err(WignerClassificationError::new(
            "msg_gauge: no spin ops mapped to MSG",
        ));
    }
    let w = w_sum / (n_mapped as f64);
    let h_dim = spin_lg_op_indices
        .iter()
        .enumerate()
        .find_map(|(local, &global)| {
            let op = h_spin_seitz.get(global as usize)?;
            if op.rot == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
                spin_chars_real.get(local).map(|c| c.abs().round().max(1.0))
            } else {
                None
            }
        })
        .unwrap_or(1.0);

    let tol = 1e-6;
    if (w.re - h_dim).abs() < tol && w.im.abs() < tol {
        MSG_GAUGE_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::A)
    } else if (w.re + h_dim).abs() < tol && w.im.abs() < tol {
        MSG_GAUGE_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::B)
    } else if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        MSG_GAUGE_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        MSG_GAUGE_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::B)
    } else if w.norm() < tol {
        MSG_GAUGE_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::C)
    } else {
        let count = MSG_GAUGE_W_FAIL.fetch_add(1, Ordering::Relaxed);
        if count < 3 {
            eprintln!(
                "  W_FAIL#{}: sg={} k=({}/{},{}/{},{}/{}) dim={:.0} n_mapped={} w=({:.6},{:.6}) |w|={:.6}",
                count,
                ctx.sg,
                kx,
                kd,
                ky,
                kd,
                kz,
                kd,
                h_dim,
                n_mapped,
                w.re,
                w.im,
                w.norm()
            );
        }
        Err(WignerClassificationError::with_value(
            "msg_gauge: non-quantized Wigner indicator for spinor",
            w.norm(),
        ))
    }
}

/// Wigner test for spinor (double-valued) irreps using SU(2) composition.
///
/// Unlike scalar irreps, spinor irreps live in the double group where each
/// spatial operation {R|t} has two lifts: g and Ēg (Ē = 2π rotation = -1).
/// The spinor character table from spin.dat assigns characters to specific
/// double-group elements, indexed by SU(2) lift.
///
/// # Arguments
/// * `ctx` — [`SpinLiftContext`] with H's and G's spin operations.
///   For black-white MSGs, $$a_0$$'s SU(2) lift is looked up in G's spin ops
///   because $$g_0 \in G \setminus H$$.
/// * `spin_chars` — first `n_lg_ops` values are little-group characters
/// * `spin_lg_op_indices` — local char position → global spin op index
///
/// # Returns
/// `None` if spin ops are unavailable or result is non-quantized.
pub fn wigner_classify_spinor(
    ctx: &SpinLiftContext,
    input: SpinorWignerInput<'_>,
    group: WignerGroupContext<'_>,
    setting_xf: Option<&SettingTransform>,
    anti_lg_indices: Option<&[usize]>,
) -> Result<CorepType, WignerClassificationError> {
    let spin_chars_real = input.characters_real;
    let spin_chars_imag = input.characters_imag;
    let spin_lg_op_indices = input.operation_indices;
    let n_lg_ops = spin_lg_op_indices.len();
    let [kx, ky, kz] = input.k_vector.numerators;
    let kd = input.k_vector.denominator;
    let unitary_mag_indices = group.unitary_indices;
    let mag_seitz = group.magnetic_ops;
    let a0_idx = group.antiunitary_representative;
    let (h_spin_rots, h_spin_trans, h_spin_su2) = ctx.h;
    let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;
    if !spin_table_is_well_formed(h_spin_rots, h_spin_trans, h_spin_su2)
        || !spin_table_is_well_formed(g_spin_rots, g_spin_trans, g_spin_su2)
    {
        return Err(WignerClassificationError::new(
            "spinor operation tables have inconsistent lengths",
        ));
    }
    if n_lg_ops == 0
        || n_lg_ops != spin_lg_op_indices.len()
        || spin_chars_real.len() < n_lg_ops
        || (!spin_chars_imag.is_empty() && spin_chars_imag.len() < n_lg_ops)
    {
        return Err(WignerClassificationError::new(
            "spinor little-group indices and characters have inconsistent lengths",
        ));
    }
    let h_spin_count = h_spin_rots.len() / 9;
    if spin_lg_op_indices
        .iter()
        .any(|&index| index as usize >= h_spin_count)
    {
        return Err(WignerClassificationError::new(
            "spinor little-group operation index is out of range",
        ));
    }
    if a0_idx >= mag_seitz.len()
        || !mag_seitz[a0_idx].timerev
        || unitary_mag_indices
            .iter()
            .any(|&index| index >= mag_seitz.len() || mag_seitz[index].timerev)
    {
        return Err(WignerClassificationError::new(
            "magnetic operation index is out of range or has the wrong time-reversal role",
        ));
    }
    if let Some(indices) = anti_lg_indices
        && indices
            .iter()
            .any(|&index| index >= mag_seitz.len() || !mag_seitz[index].timerev)
    {
        return Err(WignerClassificationError::new(
            "antiunitary little-group operation index is out of range or unitary",
        ));
    }

    // ── Direct anti-coset path (frame-aware, primary) ────────────────────
    // This path uses setting_xf to transform all operations into the
    // ISOTROPY data-Hall frame, matching the spin table and character
    // conventions.  Tried first because the MSG-gauge primary path does
    // not use the setting transform and can give wrong results when the
    // MSG frame differs from the data Hall frame.
    //
    // Use pre-computed setting-aware indices when available (caller
    // already filtered with filter_little_group_with_transform).
    let indices: std::borrow::Cow<[usize]> = if let Some(ix) = anti_lg_indices {
        std::borrow::Cow::Borrowed(ix)
    } else {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let pure_translations: Vec<[f64; 3]> = mag_seitz
            .iter()
            .filter(|op| !op.timerev && op.rot == identity)
            .map(|op| op.trans)
            .collect();
        let v: Vec<usize> = mag_seitz
            .iter()
            .enumerate()
            .filter(|(_, op)| {
                op.timerev && seitz_preserves_k(&op.rot, true, &pure_translations, kx, ky, kz, kd)
            })
            .map(|(idx, _)| idx)
            .collect();
        std::borrow::Cow::Owned(v)
    };
    // ── Direct anti-coset path (primary, frame-aware) ────────────────────
    // Only MissingSpinData is allowed to fall back to the legacy path.
    // All other errors must propagate as classification errors.
    match wigner_classify_spinor_direct_anti_diagnostic(
        ctx, input, &indices, mag_seitz, setting_xf, None,
    ) {
        Ok(result) => return Ok(result),
        Err(DirectAntiFailure::MissingSpinData) => { /* fall through to legacy */ }
        Err(e) => {
            return Err(WignerClassificationError::new(format!(
                "spinor direct anti path failed: {:?}",
                e
            )));
        }
    }

    // ── Legacy MSG-gauge primary path (fallback) ─────────────────────────
    // Used only when the frame-aware direct path lacks spin data.
    wigner_classify_spinor_primary(ctx, input, group)
}

fn wigner_classify_spinor_primary(
    ctx: &SpinLiftContext,
    input: SpinorWignerInput<'_>,
    group: WignerGroupContext<'_>,
) -> Result<CorepType, WignerClassificationError> {
    let spin_chars_real = input.characters_real;
    let spin_chars_imag = input.characters_imag;
    let spin_lg_op_indices = input.operation_indices;
    let n_lg_ops = spin_lg_op_indices.len();
    let [kx, ky, kz] = input.k_vector.numerators;
    let kd = input.k_vector.denominator;
    let mag_seitz = group.magnetic_ops;
    let a0_idx = group.antiunitary_representative;
    // Try MSG-gauge path first (correct coordinate frame).
    if let Ok(result) = wigner_classify_spinor_msg_gauge(ctx, input, group) {
        return Ok(result);
    }

    let (h_spin_rots, h_spin_trans, h_spin_su2) = ctx.h;
    let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;

    if h_spin_rots.is_empty()
        || h_spin_trans.is_empty()
        || h_spin_su2.is_empty()
        || g_spin_rots.is_empty()
        || g_spin_trans.is_empty()
        || g_spin_su2.is_empty()
        || n_lg_ops == 0
        || spin_lg_op_indices.is_empty()
    {
        return Err(WignerClassificationError::new(
            "spinor primary: missing spin data",
        ));
    }

    // Spin Seitz ops in Bilbao setting — canonical little co-group representatives.
    let h_spin_seitz = build_spin_seitz(h_spin_rots, h_spin_trans);
    if h_spin_seitz.is_empty() {
        return Err(WignerClassificationError::new(
            "spinor primary: empty spin Seitz ops",
        ));
    }

    // H_op → spin global index mapping (for matching (a₀h)² back to spin ops)
    #[cfg(any(test, feature = "debug-corep"))]
    let _h_to_spin = build_h_to_spin_map(group.unitary_ops, &h_spin_seitz, spin_lg_op_indices);

    // global spin op index → local character table position
    let global_to_local: std::collections::HashMap<usize, usize> = spin_lg_op_indices
        .iter()
        .enumerate()
        .map(|(local, &global)| (global as usize, local))
        .collect();

    // Infer central parity eta_ebar = chi(Ebar)/chi(E) from the character table.
    // For genuine spinor irreps: -1.0.  For single-valued: +1.0.
    let eta_ebar = infer_eta_ebar(
        spin_chars_real,
        spin_lg_op_indices,
        &h_spin_seitz,
        h_spin_su2,
    )
    .unwrap_or(-1.0);

    let a0 = &mag_seitz[a0_idx];

    // a₀ SU(2) lift: rotation-only lookup in G's spin ops.
    let g_spin_seitz = build_spin_seitz(g_spin_rots, g_spin_trans);
    let a0_match = g_spin_seitz
        .iter()
        .position(|s| s.rot == a0.rot)
        .or_else(|| {
            // For improper rotations (det=-1, e.g. mirrors) that aren't
            // in the spin database, try the proper rotation counterpart
            // R_proper = -R (det=+1).  The SU(2) lift of an improper
            // rotation is ±U(proper_rotation_part), and the ± sign is
            // handled later by central element detection.
            let r: Mat3I = [
                [-a0.rot[0][0], -a0.rot[0][1], -a0.rot[0][2]],
                [-a0.rot[1][0], -a0.rot[1][1], -a0.rot[1][2]],
                [-a0.rot[2][0], -a0.rot[2][1], -a0.rot[2][2]],
            ];
            g_spin_seitz.iter().position(|s| s.rot == r)
        })
        .ok_or_else(|| WignerClassificationError::new("spinor primary: a0 not found in G spin"))?;
    let u_a0 = spin_su2_at(g_spin_su2, a0_match)
        .ok_or_else(|| WignerClassificationError::new("spinor primary: a0 SU(2) lift missing"))?;

    // Data-Hall → spin-table origin shift. a₀ comes from MSG/data-Hall and
    // needs conversion; h is from spin_seitz and is already in spin convention.
    let origin = hall_to_spin_origin_for_sg(ctx.sg);
    let to_bilbao =
        |rot: Mat3I, trans: [f64; 3]| -> [f64; 3] { apply_spin_origin_shift(rot, trans, origin) };

    let a0_bilbao = SeitzOp::new(a0.rot, to_bilbao(a0.rot, a0.trans), false);

    // ── Wigner sum over the little co-group ──
    // W = (1/|H₀|) Σ_{R∈H₀} χ̃(a₀·h_R)
    //
    // Co-character:
    //   χ̃(a₀ h) = ± χ_DG( (a₀h)² ) · exp(2πi k·L)
    //
    // Central sign from SU(2):
    //   (U_a₀·U_h)² ≈ ± U_{(a₀h)²}
    //   (+) → canonical lift;  (-) → Ē-lift, character flips sign.
    let mut w_sum = Complex64::ZERO;

    for &global_index in spin_lg_op_indices.iter().take(n_lg_ops) {
        let global_spin_idx = global_index as usize;

        // Canonical h in Bilbao (already in the correct setting).
        let h_spin = &h_spin_seitz[global_spin_idx];
        let u_h = spin_su2_at(h_spin_su2, global_spin_idx).ok_or_else(|| {
            WignerClassificationError::new("spinor primary: h SU(2) lift missing")
        })?;

        // Spatial: (a₀ h)² in Bilbao.
        let (g0h, l1) = compose_seitz(&a0_bilbao, h_spin);
        let (sq, lattice_sq) = square_seitz(&g0h);

        // Match square's rotation back to spin ops, preferring LG candidates first.
        // Priority: full Seitz in LG → unique rotation in LG → global rotation.
        // This avoids position() picking a non-LG candidate when an LG candidate exists.
        let (sq_spin_idx, sq_in_lg, _match_kind) = match find_sq_spin_lg_first(
            &sq,
            &h_spin_seitz,
            spin_lg_op_indices,
            centering_shifts_for_sg(ctx.sg),
        ) {
            Some(v) => v,
            None => {
                eprintln!("  WIGNER_SPINOR: sq_rot not in spin ops, aborting case");
                return Err(WignerClassificationError::new(
                    "spinor primary: sq_rot not in spin ops",
                ));
            }
        };

        // SU(2): (U_a₀·U_h)² vs canonical U_{(a₀h)²}.
        let u_g0h = su2_compose(&u_a0, &u_h);
        let u_sq = SquareKernel::OldU2.apply(&u_g0h);
        let u_k = spin_su2_at(h_spin_su2, sq_spin_idx).ok_or_else(|| {
            WignerClassificationError::new("spinor primary: sq SU(2) lift missing")
        })?;

        // Central element detection.
        // central=true: u_sq ≈ -u_k (differs by Ebar)
        // central=false: u_sq ≈ u_k (same lift)
        let central = match su2_same_up_to_sign(&u_sq, &u_k) {
            Some(false) => {
                SU2_REL_SAME.fetch_add(1, Ordering::Relaxed);
                false
            }
            Some(true) => {
                SU2_REL_EBAR.fetch_add(1, Ordering::Relaxed);
                true
            }
            None => {
                SU2_REL_NONE.fetch_add(1, Ordering::Relaxed);
                // Scan same-rotation candidates
                let lg_set: std::collections::HashSet<usize> =
                    spin_lg_op_indices.iter().map(|&x| x as usize).collect();
                let mut matched_other_lg = false;
                let mut matched_other_global = false;
                let mut has_cand = false;
                for (ci, cs) in h_spin_seitz.iter().enumerate() {
                    if cs.rot != sq.rot || ci == sq_spin_idx {
                        continue;
                    }
                    has_cand = true;
                    if let Some(uc) = spin_su2_at(h_spin_su2, ci)
                        && su2_same_up_to_sign(&u_sq, &uc).is_some()
                    {
                        if lg_set.contains(&ci) {
                            matched_other_lg = true;
                        } else {
                            matched_other_global = true;
                        }
                    }
                }
                if matched_other_lg {
                    NONE_MATCH_OTHER_LG.fetch_add(1, Ordering::Relaxed);
                } else if matched_other_global {
                    NONE_MATCH_OTHER_GLOBAL.fetch_add(1, Ordering::Relaxed);
                } else if has_cand {
                    NONE_NO_MATCH_HAS_CAND.fetch_add(1, Ordering::Relaxed);
                } else {
                    NONE_NO_CANDIDATE.fetch_add(1, Ordering::Relaxed);
                }
                // Det distribution
                let det_a0 = a0.rot[0][0]
                    * (a0.rot[1][1] * a0.rot[2][2] - a0.rot[1][2] * a0.rot[2][1])
                    - a0.rot[0][1] * (a0.rot[1][0] * a0.rot[2][2] - a0.rot[1][2] * a0.rot[2][0])
                    + a0.rot[0][2] * (a0.rot[1][0] * a0.rot[2][1] - a0.rot[1][1] * a0.rot[2][0]);
                if det_a0 > 0 {
                    NONE_DET_A0_P1.fetch_add(1, Ordering::Relaxed);
                } else {
                    NONE_DET_A0_M1.fetch_add(1, Ordering::Relaxed);
                }
                let det_g0h = g0h.rot[0][0]
                    * (g0h.rot[1][1] * g0h.rot[2][2] - g0h.rot[1][2] * g0h.rot[2][1])
                    - g0h.rot[0][1]
                        * (g0h.rot[1][0] * g0h.rot[2][2] - g0h.rot[1][2] * g0h.rot[2][0])
                    + g0h.rot[0][2]
                        * (g0h.rot[1][0] * g0h.rot[2][1] - g0h.rot[1][1] * g0h.rot[2][0]);
                if det_g0h > 0 {
                    NONE_DET_G0H_P1.fetch_add(1, Ordering::Relaxed);
                } else {
                    NONE_DET_G0H_M1.fetch_add(1, Ordering::Relaxed);
                }
                // Test 6 alternative antiunitary square formulas
                let u_cj = conj_pauli(&u_g0h);
                let alts = [
                    su2_compose(&u_g0h, &u_g0h),             // U^2 (raw)
                    neg_pauli(&su2_compose(&u_g0h, &u_g0h)), // -U^2
                    su2_compose(&u_g0h, &u_cj),              // U U*
                    neg_pauli(&su2_compose(&u_g0h, &u_cj)),  // -U U*
                    su2_compose(&u_cj, &u_g0h),              // U* U
                    neg_pauli(&su2_compose(&u_cj, &u_g0h)),  // -U* U
                ];
                let matches: Vec<bool> = alts
                    .iter()
                    .map(|a| su2_same_up_to_sign(a, &u_k).is_some())
                    .collect();
                if matches[0] {
                    NONE_ALT_RAW.fetch_add(1, Ordering::Relaxed);
                }
                if matches[1] {
                    NONE_ALT_NEG_RAW.fetch_add(1, Ordering::Relaxed);
                }
                if matches[2] {
                    NONE_ALT_UUSTAR.fetch_add(1, Ordering::Relaxed);
                }
                if matches[3] {
                    NONE_ALT_NEG_UUSTAR.fetch_add(1, Ordering::Relaxed);
                }
                if matches[4] {
                    NONE_ALT_STARU.fetch_add(1, Ordering::Relaxed);
                }
                if matches[5] {
                    NONE_ALT_NEG_STARU.fetch_add(1, Ordering::Relaxed);
                }
                if !matches.iter().any(|&m| m) {
                    NONE_ALT_NONE.fetch_add(1, Ordering::Relaxed);
                }
                // J-insertion antiunitary square: J = i*sigma_y = [0,0,1,0]
                let j = [0.0, 0.0, 1.0, 0.0];
                let ju = su2_compose(&j, &u_g0h);
                let uj = su2_compose(&u_g0h, &j);
                let sq_ju = su2_compose(&ju, &conj_pauli(&ju));
                let sq_uj = su2_compose(&uj, &conj_pauli(&uj));
                let j_matches = [
                    su2_same_up_to_sign(&sq_ju, &u_k).is_some(),
                    su2_same_up_to_sign(&neg_pauli(&sq_ju), &u_k).is_some(),
                    su2_same_up_to_sign(&sq_uj, &u_k).is_some(),
                    su2_same_up_to_sign(&neg_pauli(&sq_uj), &u_k).is_some(),
                ];
                if j_matches[0] {
                    NONE_JU_JU_STAR.fetch_add(1, Ordering::Relaxed);
                }
                if j_matches[1] {
                    NONE_NEG_JU_JU_STAR.fetch_add(1, Ordering::Relaxed);
                }
                if j_matches[2] {
                    NONE_UJ_UJ_STAR.fetch_add(1, Ordering::Relaxed);
                }
                if j_matches[3] {
                    NONE_NEG_UJ_UJ_STAR.fetch_add(1, Ordering::Relaxed);
                }
                if !j_matches.iter().any(|&m| m) {
                    NONE_J_NONE.fetch_add(1, Ordering::Relaxed);
                }
                // G-gauge oracle: all SU(2) in G spin database
                if let Some((h_g_idx, _)) = find_spin_in_db(h_spin, &g_spin_seitz) {
                    if let Some((sq_g_idx, _)) = find_spin_in_db(&sq, &g_spin_seitz) {
                        if let Some(u_h_g_val) = spin_su2_at(g_spin_su2, h_g_idx) {
                            if let Some(u_sq_g_table_val) = spin_su2_at(g_spin_su2, sq_g_idx) {
                                let u_g0h_g = su2_compose(&u_a0, &u_h_g_val);
                                let u_sq_g = su2_compose(&u_g0h_g, &u_g0h_g);
                                match su2_same_up_to_sign(&u_sq_g, &u_sq_g_table_val) {
                                    Some(false) => GGAUGE_SAME.fetch_add(1, Ordering::Relaxed),
                                    Some(true) => GGAUGE_EBAR.fetch_add(1, Ordering::Relaxed),
                                    None => GGAUGE_NONE.fetch_add(1, Ordering::Relaxed),
                                };
                            } else {
                                GGAUGE_NONE.fetch_add(1, Ordering::Relaxed);
                            }
                        } else {
                            GGAUGE_NONE.fetch_add(1, Ordering::Relaxed);
                        }
                    } else {
                        GGAUGE_SQ_LOOKUP_FAIL.fetch_add(1, Ordering::Relaxed);
                    }
                } else {
                    GGAUGE_H_LOOKUP_FAIL.fetch_add(1, Ordering::Relaxed);
                }
                return Err(WignerClassificationError::new(
                    "spinor primary: SU(2) square not comparable to canonical lift",
                ));
            }
        };

        // Character from LG table (if sq ∈ LG) or from extended table.
        let sq_local_idx = if sq_in_lg {
            *global_to_local.get(&sq_spin_idx).ok_or_else(|| {
                WignerClassificationError::new("spinor primary: spin→local lookup missing")
            })?
        } else {
            // sq outside LG: need extended character, abort case.
            eprintln!(
                "  WIGNER_SPINOR: sq[{}] not in LG idxs, aborting case",
                sq_spin_idx
            );
            return Err(WignerClassificationError::new(
                "spinor primary: square not in LG",
            ));
        };

        // Bloch phase from total lattice shift.
        let r_l1 = mat_vec_i32(&g0h.rot, &l1);
        let total_l = add3(&lattice_sq, &add3(&l1, &r_l1));
        let phase = bloch_phase(kx, ky, kz, kd, &total_l);

        let chi0 = Complex64::new(
            spin_chars_real[sq_local_idx],
            spin_chars_imag.get(sq_local_idx).copied().unwrap_or(0.0),
        );
        let chi = if central { eta_ebar * chi0 } else { chi0 };

        w_sum += chi * phase;
    }

    // W = w_sum / |H₀|  (little co-group order = n_lg_ops)
    let w = w_sum / (n_lg_ops as f64);

    // Robust dimension: find the identity canonical lift in spin_lg_op_indices.
    // spin_chars[0] may NOT be χ(E) — canonical lifts may be reordered.
    let id_rot: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let u_id = [1.0, 0.0, 0.0, 0.0];
    let h_dim = spin_lg_op_indices
        .iter()
        .map(|&x| x as usize)
        .find_map(|si| {
            let sop = h_spin_seitz.get(si)?;
            if sop.rot != id_rot {
                return None;
            }
            let u = spin_su2_at(h_spin_su2, si)?;
            if su2_same_up_to_sign(&u, &u_id) != Some(false) {
                return None;
            }
            let local = *global_to_local.get(&si)?;
            spin_chars_real
                .get(local)
                .map(|&c| c.abs().round().max(1.0))
        })
        .unwrap_or_else(|| {
            spin_chars_real
                .first()
                .map(|&c| c.abs().round().max(1.0))
                .unwrap_or(1.0)
        });

    // Classification.  For spinor irreps: W = +dim → Type A, W = -dim → Type B,
    // W = 0 → Type C.  Also accept W = ±1 for 1D irreps.
    let tol = 1e-6;
    if (w.re - h_dim).abs() < tol && w.im.abs() < tol {
        OLD_PATH_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::A)
    } else if (w.re + h_dim).abs() < tol && w.im.abs() < tol {
        OLD_PATH_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::B)
    } else if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        OLD_PATH_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        OLD_PATH_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::B)
    } else if w.norm() < tol {
        OLD_PATH_OK.fetch_add(1, Ordering::Relaxed);
        Ok(CorepType::C)
    } else {
        OLD_PATH_FAIL.fetch_add(1, Ordering::Relaxed);
        Err(WignerClassificationError::with_value(
            "spinor primary: non-quantized Wigner indicator",
            w.norm(),
        ))
    }
}

// ── Type A intertwiner + matrix utilities ────────────────────────────────────

include!("wigner_extra.rs");

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
/// - Unitary: $$\tilde{\chi}(h) = \chi_i(h) + \chi_{partner}(h)$$
///   (this becomes $2\,\mathrm{Re}[\chi_i(h)]$ only when the caller has
///   proved that the partner is the `a₀` conjugate under direct pure time
///   reversal)
/// - Anti-unitary: $$\tilde{\chi}(a_0 h) = 0$$
///
/// # Parameters
///
/// * `corep_type` — result of [`wigner_classify`]
/// * `mag_ops` — magnetic symmetry operations (for timerev flags)
/// * `mag_lg_indices` — which magnetic ops are in the little group
/// * `op_map` — for each magnetic op, the corresponding H op index (or None)
/// * `h_chars` — H's irrep character table (real-valued for PIR)
/// * `partner_chars` — Type C's operation-aware character table for the
///   paired irrep, supplied by a caller that established the pairing
pub(crate) fn build_corep_chars(
    corep_type: &CorepType,
    mag_ops: &SymmetryOps,
    mag_lg_indices: &[usize],
    op_map: &[Option<usize>],
    h_chars: &[f64],
    partner_chars: Option<&[f64]>, // Type C: character table of paired irrep
    au_chars: Option<&[f64]>,      // for Type A: pre-computed antiunitary chars
) -> Result<Vec<f64>, &'static str> {
    let n_lg = mag_lg_indices.len();
    if mag_lg_indices
        .iter()
        .any(|&index| index >= mag_ops.len() || index >= op_map.len())
    {
        return Err("magnetic little-group operation index is out of range");
    }
    if corep_type == &CorepType::A && au_chars.is_some_and(|chars| chars.len() < n_lg) {
        return Err("Type-A antiunitary character table is shorter than the magnetic little group");
    }
    if corep_type == &CorepType::C && partner_chars.is_none() {
        return Err("Type-C corep requires explicit partner characters");
    }
    let mut chars = vec![0.0; n_lg];

    for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
        let is_anti = mag_ops[mag_idx].time_reversal;
        let h_idx = op_map[mag_idx];

        match corep_type {
            CorepType::A => {
                if is_anti {
                    if let Some(ac) = au_chars
                        && out_idx < ac.len()
                    {
                        chars[out_idx] = ac[out_idx];
                    }
                } else {
                    let hi = h_idx.ok_or("unitary magnetic operation is missing its H mapping")?;
                    chars[out_idx] = *h_chars
                        .get(hi)
                        .ok_or("H-operation mapping exceeds the character table")?;
                }
            }
            CorepType::B => {
                // Kramers doubling: dimension 2d
                if is_anti {
                    chars[out_idx] = 0.0;
                } else {
                    let hi = h_idx.ok_or("unitary magnetic operation is missing its H mapping")?;
                    chars[out_idx] = 2.0
                        * h_chars
                            .get(hi)
                            .ok_or("H-operation mapping exceeds the character table")?;
                }
            }
            CorepType::C => {
                // Paired with conjugate irrep: dimension 2d
                if is_anti {
                    chars[out_idx] = 0.0;
                } else {
                    let hi = h_idx.ok_or("unitary magnetic operation is missing its H mapping")?;
                    let chi_i = *h_chars
                        .get(hi)
                        .ok_or("H-operation mapping exceeds the character table")?;
                    let pc =
                        partner_chars.ok_or("Type-C corep requires explicit partner characters")?;
                    let chi_partner = *pc
                        .get(hi)
                        .ok_or("partner character table is shorter than the H mapping")?;
                    chars[out_idx] = chi_i + chi_partner;
                }
            }
        }
    }

    Ok(chars)
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

    fn gamma() -> KVector {
        KVector::new([0, 0, 0], 1)
    }

    fn group<'a>(
        unitary_indices: &'a [usize],
        magnetic_ops: &'a [SeitzOp],
        unitary_ops: &'a [SeitzOp],
        antiunitary_representative: usize,
    ) -> WignerGroupContext<'a> {
        WignerGroupContext {
            unitary_indices,
            magnetic_ops,
            unitary_ops,
            antiunitary_representative,
        }
    }

    fn spinor_input<'a>(
        characters_real: &'a [f64],
        characters_imag: &'a [f64],
        operation_indices: &'a [u16],
    ) -> SpinorWignerInput<'a> {
        SpinorWignerInput {
            characters_real,
            characters_imag,
            operation_indices,
            k_vector: gamma(),
        }
    }

    /// Seitz composition: identity ∘ identity = identity.
    #[test]
    fn test_compose_identity() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0], false);
        let (result, lattice) = compose_seitz(&id, &id);
        assert_eq!(result.rot, id.rot);
        assert_eq!(result.trans, [0.0, 0.0, 0.0]);
        assert_eq!(lattice, [0, 0, 0]);
        assert!(!result.timerev);
    }

    /// Seitz composition: timerev XOR.
    #[test]
    fn test_compose_timerev() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0], false);
        let a0 = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0], true);
        let (r1, _) = compose_seitz(&a0, &id); // anti ∘ unitary = anti
        assert!(r1.timerev);
        let (r2, _) = compose_seitz(&a0, &a0); // anti ∘ anti = unitary
        assert!(!r2.timerev);
        let (r3, _) = compose_seitz(&id, &a0); // unitary ∘ anti = anti
        assert!(r3.timerev);
    }

    /// Seitz composition: translation arithmetic.
    #[test]
    fn test_compose_translation() {
        let g1 = SeitzOp::new([[0, -1, 0], [1, 0, 0], [0, 0, 1]], [0.0, 0.0, 0.5], false);
        let g2 = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.5, 0.0, 0.0], false);
        let (result, lattice) = compose_seitz(&g1, &g2);
        // R = [[0,-1,0],[1,0,0],[0,0,1]]
        // t = [0,0,0.5] + R·[0.5,0,0] = [0,0,0.5] + [0,0.5,0] = [0, 0.5, 0.5]
        assert_eq!(result.trans, [0.0, 0.5, 0.5]);
        assert_eq!(lattice, [0, 0, 0]);
    }

    /// Seitz composition with lattice overflow.
    #[test]
    fn test_compose_lattice_shift() {
        let g1 = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.7, 0.0, 0.0], false);
        let g2 = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.5, 0.0, 0.0], false);
        let (result, lattice) = compose_seitz(&g1, &g2);
        // t = 0.7 + 0.5 = 1.2 → 0.2 with lattice shift [1,0,0]
        assert!((result.trans[0] - 0.2).abs() < 1e-9);
        assert_eq!(lattice, [1, 0, 0]);
    }

    #[test]
    fn setting_transform_singular_basis_returns_none() {
        let transform = SettingTransform {
            basis: [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            origin: [0.25, 0.0, 0.0],
        };
        let rotation = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let translation = [0.5, 0.0, 0.0];

        assert!(transform.transform_rotation(&rotation).is_none());
        assert!(
            transform
                .transform_translation(&rotation, &translation)
                .is_none()
        );
        assert!(transform.transform_seitz(&rotation, &translation).is_none());
    }

    #[test]
    fn setting_transform_near_singular_basis_returns_none() {
        let transform = SettingTransform {
            basis: [[1e-6, 0.0, 0.0], [0.0, 1e-6, 0.0], [0.0, 0.0, 1e-6]],
            origin: [0.0; 3],
        };
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        assert!(transform.transform_rotation(&identity).is_none());
    }

    #[test]
    fn setting_transform_nontrivial_valid_control() {
        let transform = SettingTransform {
            basis: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            origin: [0.25, 0.0, 0.0],
        };
        let rotation = [[1, 0, 0], [0, -1, 0], [0, 0, -1]];
        let translation = [0.5, 0.0, 0.0];
        let expected_rotation = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]];
        let expected_translation = [0.5, 0.5, 0.0];

        assert_eq!(
            transform.transform_rotation(&rotation),
            Some(expected_rotation)
        );
        assert_eq!(
            transform.transform_translation(&rotation, &translation),
            Some(expected_translation)
        );
        assert_eq!(
            transform.transform_seitz(&rotation, &translation),
            Some((expected_rotation, expected_translation))
        );
    }

    #[test]
    fn setting_transform_rejects_invalid_identity_origin_candidate() {
        // SG7 Hall 23 and Hall 21 have the same rotation multiset but
        // different glide translations.  A greedy origin solve can pair the
        // rotations and suggest T=I even though the full Seitz sets differ.
        let source = SymmetryOps::from_database(23).unwrap();
        let target = SymmetryOps::from_database(21).unwrap();
        let source_rots: Vec<Mat3I> = source.operations.iter().map(|op| op.rotation).collect();
        let source_trans: Vec<[f64; 3]> =
            source.operations.iter().map(|op| op.translation).collect();
        let target_rots: Vec<Mat3I> = target.operations.iter().map(|op| op.rotation).collect();
        let target_trans: Vec<[f64; 3]> =
            target.operations.iter().map(|op| op.translation).collect();
        let target_seitz = ops_to_seitz(&target);

        let candidates =
            find_setting_transform(&source_rots, &source_trans, &target_rots, &target_trans);
        assert!(!candidates.is_empty());
        assert!(
            candidates
                .iter()
                .all(|candidate| candidate.basis != SettingTransform::identity().basis)
        );
        for candidate in candidates {
            for op in &source.operations {
                let (rotation, translation) = candidate
                    .transform_seitz(&op.rotation, &op.translation)
                    .unwrap();
                assert!(find_seitz(&rotation, &translation, &target_seitz).is_some());
            }
        }
    }

    /// filter_little_group: antiunitary ops use -Rk ≡ k.
    #[test]
    fn test_filter_antiunitary_k() {
        // k = (0, 0, 1)/2 = Z point
        // Anti-unitary op with R = [[0,-1,0],[1,0,0],[0,0,-1]] (4' about 001)
        // R·(0,0,1) = (0,0,-1), so -R·k - k = (0,0,1) - (0,0,1) = (0,0,0) ≡ 0 ✓
        let ops = SymmetryOps::from_parallel_owned(
            vec![
                [[0, -1, 0], [1, 0, 0], [0, 0, -1]],
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            ],
            vec![[0.0; 3], [0.0; 3]],
            vec![true, false],
        )
        .unwrap();
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
        let ops = SymmetryOps::from_parallel_owned(
            vec![
                [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
            ],
            vec![[0.0; 3], [0.0; 3]],
            vec![true, true],
        )
        .unwrap();
        let lg = filter_little_group(1, 0, 0, 8, &ops);
        assert_eq!(lg.len(), 1, "Only mx' should preserve k=(1/8,0,0)");
    }

    /// Simple Wigner test: P1 with only identity.
    #[test]
    fn test_wigner_trivial() {
        // a₀ = θ (anti-unitary identity), h = id
        // (a₀·id)² = id² = id, χ(id)=1.0 → W=1.0 → type A
        let mag_seitz = vec![
            SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0], false), // id
            SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0], true),  // θ
        ];
        let h_seitz = vec![SeitzOp::new(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [0.0, 0.0, 0.0],
            false,
        )];
        let result = wigner_classify(&[1.0], group(&[0], &mag_seitz, &h_seitz, 1), gamma());
        assert_eq!(result, Ok(CorepType::A));
    }

    #[test]
    fn test_scalar_wigner_rejects_empty_and_incomplete_inputs() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let theta = SeitzOp::new(id.rot, [0.0; 3], true);
        let mag = vec![id.clone(), theta];
        let h = vec![id];

        let real_empty = wigner_classify(&[1.0], group(&[], &mag, &h, 1), gamma()).unwrap_err();
        assert!(real_empty.reason.contains("unitary little group is empty"));
        let cir_empty =
            wigner_classify_cir(&[1.0, 0.0], group(&[], &mag, &h, 1), gamma()).unwrap_err();
        assert!(cir_empty.reason.contains("unitary little group is empty"));

        let h_two = vec![
            mag[0].clone(),
            SeitzOp::new([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], [0.0; 3], false),
        ];
        let real_short =
            wigner_classify(&[1.0], group(&[0], &mag, &h_two, 1), gamma()).unwrap_err();
        assert!(real_short.reason.contains("character table length"));
        let cir_short =
            wigner_classify_cir(&[1.0, 0.0], group(&[0], &mag, &h_two, 1), gamma()).unwrap_err();
        assert!(cir_short.reason.contains("CIR character count"));
    }

    #[test]
    fn test_scalar_wigner_rejects_invalid_indices_and_roles() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let theta = SeitzOp::new(id.rot, [0.0; 3], true);
        let mag = vec![id.clone(), theta];
        let h = vec![id];

        for result in [
            wigner_classify(&[1.0], group(&[0], &mag, &h, 2), gamma()),
            wigner_classify(&[1.0], group(&[2], &mag, &h, 1), gamma()),
            wigner_classify(&[1.0], group(&[1], &mag, &h, 1), gamma()),
        ] {
            assert!(result.is_err());
        }
        for result in [
            wigner_classify_cir(&[1.0, 0.0], group(&[0], &mag, &h, 2), gamma()),
            wigner_classify_cir(&[1.0, 0.0], group(&[2], &mag, &h, 1), gamma()),
            wigner_classify_cir(&[1.0, 0.0], group(&[1], &mag, &h, 1), gamma()),
        ] {
            assert!(result.is_err());
        }
    }

    #[test]
    fn test_scalar_wigner_rejects_missing_square_match() {
        let id_rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let c4z_rot = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let id = SeitzOp::new(id_rot, [0.0; 3], false);
        let mag = vec![id.clone(), SeitzOp::new(c4z_rot, [0.0; 3], true)];
        let h = vec![id];

        let real = wigner_classify(&[1.0], group(&[0], &mag, &h, 1), gamma()).unwrap_err();
        assert!(real.reason.contains("absent from the unitary little group"));
        let cir = wigner_classify_cir(&[1.0, 0.0], group(&[0], &mag, &h, 1), gamma()).unwrap_err();
        assert!(cir.reason.contains("absent from the unitary little group"));
    }

    #[test]
    fn direct_anti_coset_rejects_bad_indices_and_characters() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let theta = SeitzOp::new(id.rot, [0.0; 3], true);
        let mag = vec![theta];
        let h = vec![id];

        assert!(wigner_direct_anti_coset(&[1.0, 0.0], &[1], &mag, &h, gamma()).is_err());
        assert!(wigner_direct_anti_coset(&[1.0], &[0], &mag, &h, gamma()).is_err());
        assert!(wigner_direct_anti_coset(&[1.0, 0.0], &[], &mag, &h, gamma()).is_err());
    }

    #[test]
    fn direct_anti_coset_valid_type_a_control() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let theta = SeitzOp::new(id.rot, [0.0; 3], true);
        let value = wigner_direct_anti_coset(&[1.0, 0.0], &[0], &[theta], &[id], gamma()).unwrap();

        assert!((value.re - 1.0).abs() < 1e-12);
        assert!(value.im.abs() < 1e-12);
    }

    #[test]
    fn spinor_direct_anti_rejects_index_and_character_mismatches() {
        static ROTATIONS: [i32; 9] = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        static TRANSLATIONS: [f64; 3] = [0.0; 3];
        static SU2: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
        let ctx = SpinLiftContext {
            h: (&ROTATIONS, &TRANSLATIONS, &SU2),
            g: (&ROTATIONS, &TRANSLATIONS, &SU2),
            sg: 1,
        };
        let anti = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], true);

        assert_eq!(
            wigner_classify_spinor_direct_anti_diagnostic(
                &ctx,
                spinor_input(&[1.0], &[], &[0]),
                &[1],
                std::slice::from_ref(&anti),
                None,
                None,
            ),
            Err(DirectAntiFailure::IndexOutOfRange)
        );
        assert_eq!(
            wigner_classify_spinor_direct_anti_diagnostic(
                &ctx,
                spinor_input(&[], &[], &[0]),
                &[0],
                std::slice::from_ref(&anti),
                None,
                None,
            ),
            Err(DirectAntiFailure::CharacterTableMismatch)
        );

        let singular = SettingTransform {
            basis: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            origin: [0.0; 3],
        };
        assert_eq!(
            wigner_classify_spinor_direct_anti_diagnostic(
                &ctx,
                spinor_input(&[1.0], &[], &[0]),
                &[0],
                &[anti],
                Some(&singular),
                None,
            ),
            Err(DirectAntiFailure::SettingTransformFailed)
        );
    }

    #[test]
    fn spinor_classifier_rejects_legacy_fallback_index_contract_violations() {
        static ROTATIONS: [i32; 9] = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        static TRANSLATIONS: [f64; 3] = [0.0; 3];
        static SU2: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
        let ctx = SpinLiftContext {
            h: (&ROTATIONS, &TRANSLATIONS, &SU2),
            g: (&ROTATIONS, &TRANSLATIONS, &SU2),
            sg: 1,
        };
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let theta = SeitzOp::new(id.rot, [0.0; 3], true);
        let mag = vec![id.clone(), theta];
        let h = vec![id];

        let classify = |spin_indices: &[u16], unitary: &[usize], a0_idx| {
            wigner_classify_spinor(
                &ctx,
                spinor_input(&[1.0], &[], spin_indices),
                group(unitary, &mag, &h, a0_idx),
                None,
                Some(&[1]),
            )
        };

        assert!(classify(&[], &[0], 1).is_err());
        assert!(classify(&[1], &[0], 1).is_err());
        assert!(classify(&[0], &[99], 1).is_err());
        assert!(classify(&[0], &[0], usize::MAX).is_err());
        assert!(classify(&[0], &[1], 1).is_err());
    }

    #[test]
    fn exact_spin_map_rejects_lattice_and_ambiguous_matches() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let h = SeitzOp {
            rot: identity,
            trans: [0.5, 0.0, 0.0],
            timerev: false,
        };
        let exact = SeitzOp {
            rot: identity,
            trans: [0.5 + 1.0e-10, 0.0, 0.0],
            timerev: false,
        };
        assert_eq!(
            build_h_to_spin_map_exact(std::slice::from_ref(&h), std::slice::from_ref(&exact), &[0],),
            vec![Some(0)]
        );

        let integer_shift = SeitzOp {
            rot: identity,
            trans: [1.5, 0.0, 0.0],
            timerev: false,
        };
        assert_eq!(
            build_h_to_spin_map_exact(
                std::slice::from_ref(&h),
                std::slice::from_ref(&integer_shift),
                &[0],
            ),
            vec![None],
            "integer translation differences must not be matched modulo the lattice"
        );

        let different_representative = SeitzOp {
            rot: identity,
            trans: [0.0, 0.0, 0.0],
            timerev: false,
        };
        assert_eq!(
            build_h_to_spin_map_exact(
                std::slice::from_ref(&h),
                std::slice::from_ref(&different_representative),
                &[0],
            ),
            vec![None]
        );

        assert_eq!(
            build_h_to_spin_map_exact(&[h], &[exact.clone(), exact], &[0, 1]),
            vec![None],
            "duplicate exact candidates must be rejected as ambiguous"
        );
    }

    #[test]
    fn scalar_operation_map_requires_direct_full_seitz_bijection() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];
        let h = [
            SeitzOp {
                rot: identity,
                trans: [0.25, 0.0, 0.0],
                timerev: false,
            },
            SeitzOp {
                rot: inversion,
                trans: [0.0, 0.0, 0.0],
                timerev: false,
            },
        ];
        let flat_rotations = |ops: &[SeitzOp]| {
            ops.iter()
                .flat_map(|op| op.rot.iter().flat_map(|row| row.iter().copied()))
                .collect::<Vec<_>>()
        };
        let flat_translations =
            |ops: &[SeitzOp]| ops.iter().flat_map(|op| op.trans).collect::<Vec<_>>();
        let stored = [h[1].clone(), h[0].clone()];
        let mut stored_translations = flat_translations(&stored);
        stored_translations[3] += 1.0e-10;
        assert_eq!(
            build_h_to_irrep_op_map(&h, &flat_rotations(&stored), &stored_translations,),
            Some(vec![1, 0])
        );

        let mut integer_shift = flat_translations(&stored);
        integer_shift[3] += 1.0;
        assert_eq!(
            build_h_to_irrep_op_map(&h, &flat_rotations(&stored), &integer_shift),
            None
        );
        let mut different_fraction = flat_translations(&stored);
        different_fraction[3] = 0.5;
        assert_eq!(
            build_h_to_irrep_op_map(&h, &flat_rotations(&stored), &different_fraction),
            None
        );

        let duplicate_h = [h[0].clone(), h[0].clone()];
        let duplicate_stored = [h[0].clone(), h[0].clone()];
        assert_eq!(
            build_h_to_irrep_op_map(
                &duplicate_h,
                &flat_rotations(&duplicate_stored),
                &flat_translations(&duplicate_stored),
            ),
            None
        );

        assert_eq!(
            build_h_to_irrep_op_map(&h, &flat_rotations(&stored)[..8], &stored_translations),
            None
        );
        assert_eq!(
            build_h_to_irrep_op_map(&h, &flat_rotations(&stored), &stored_translations[..3]),
            None
        );
        assert_eq!(
            build_h_to_irrep_op_map(&h[..1], &flat_rotations(&stored), &stored_translations),
            None
        );
    }

    #[test]
    fn corep_character_builder_rejects_invalid_parallel_indices_and_maps() {
        let id = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let mag_ops = SymmetryOps::from_parallel_owned(
            vec![id, id],
            vec![[0.0; 3], [0.0; 3]],
            vec![false, true],
        )
        .unwrap();

        assert!(
            build_corep_chars(
                &CorepType::A,
                &mag_ops,
                &[2],
                &[Some(0), None],
                &[1.0],
                None,
                Some(&[0.0]),
            )
            .is_err()
        );
        assert!(
            build_corep_chars(&CorepType::A, &mag_ops, &[0], &[], &[1.0], None, None,).is_err()
        );
        assert!(
            build_corep_chars(
                &CorepType::A,
                &mag_ops,
                &[0, 1],
                &[Some(0), None],
                &[1.0],
                None,
                Some(&[]),
            )
            .is_err()
        );

        let valid = build_corep_chars(
            &CorepType::B,
            &mag_ops,
            &[0, 1],
            &[Some(0), None],
            &[1.0],
            None,
            None,
        )
        .unwrap();
        assert_eq!(valid, vec![2.0, 0.0]);

        assert!(
            build_corep_chars(
                &CorepType::C,
                &mag_ops,
                &[0, 1],
                &[Some(0), None],
                &[1.0],
                None,
                None,
            )
            .is_err()
        );
        let valid_type_c = build_corep_chars(
            &CorepType::C,
            &mag_ops,
            &[0, 1],
            &[Some(0), None],
            &[1.0],
            Some(&[2.0]),
            None,
        )
        .unwrap();
        assert_eq!(valid_type_c, vec![3.0, 0.0]);
    }

    #[test]
    fn type_a_character_helper_rejects_invalid_magnetic_indices() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let mag = vec![id.clone()];

        assert!(
            type_a_antiunitary_chars(&mag, &[0], &[1.0], std::slice::from_ref(&id), 0, gamma(),)
                .is_none()
        );
    }

    #[test]
    fn diagnostic_square_rejects_invalid_indices() {
        let id = SeitzOp::new([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0; 3], false);
        let theta = SeitzOp::new(id.rot, [0.0; 3], true);
        let magnetic_ops = [id.clone(), theta];

        assert!(
            debug_unwrapped_square(
                99,
                group(&[], &magnetic_ops, std::slice::from_ref(&id), 1),
                gamma(),
            )
            .is_err()
        );
    }

    #[test]
    fn test_scalar_wigner_known_b_and_c_remain_quantized() {
        let id_rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let c4_rot = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        let c2_rot = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];
        let c4_cubed_rot = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];

        let h_b = vec![
            SeitzOp::new(id_rot, [0.0; 3], false),
            SeitzOp::new(c2_rot, [0.0; 3], false),
        ];
        let mag_b = vec![
            h_b[0].clone(),
            h_b[1].clone(),
            SeitzOp::new(c4_rot, [0.0; 3], true),
        ];
        assert_eq!(
            wigner_classify(&[1.0, -1.0], group(&[0, 1], &mag_b, &h_b, 2), gamma()),
            Ok(CorepType::B)
        );
        assert_eq!(
            wigner_classify_cir(
                &[1.0, 0.0, -1.0, 0.0],
                group(&[0, 1], &mag_b, &h_b, 2),
                gamma(),
            ),
            Ok(CorepType::B)
        );

        let h_c = vec![
            SeitzOp::new(id_rot, [0.0; 3], false),
            SeitzOp::new(c4_rot, [0.0; 3], false),
            SeitzOp::new(c2_rot, [0.0; 3], false),
            SeitzOp::new(c4_cubed_rot, [0.0; 3], false),
        ];
        let mut mag_c = h_c.clone();
        mag_c.push(SeitzOp::new(id_rot, [0.0; 3], true));
        assert_eq!(
            wigner_classify(
                &[2.0, 0.0, -2.0, 0.0],
                group(&[0, 1, 2, 3], &mag_c, &h_c, 4),
                gamma(),
            ),
            Ok(CorepType::C)
        );
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

    /// Wigner type must be independent of which antiunitary op is chosen as a₀.
    #[test]
    fn test_wigner_classification_independent_of_a0() {
        use crate::irrep::corep::identify_unitary_subgroup;

        let uni = 1066usize;
        let mag_ops = crate::SymmetryOps::from_magnetic_database(uni).unwrap();
        let h_sg = identify_unitary_subgroup(uni).unwrap();
        let mag_seitz = ops_to_seitz(&mag_ops);

        let h_ops_raw = crate::SymmetryOps::from_sg(h_sg as u8).unwrap();
        let h_seitz = ops_to_seitz(&h_ops_raw);

        let h_irreps = crate::irrep::query::irreps_of(h_sg as u8);
        for ir in h_irreps.iter().filter(|r| r.k_label() == "Z" && !r.spinor) {
            let mag_lg = filter_little_group(ir.kx, ir.ky, ir.kz, ir.kd, &mag_ops);
            let unitary: Vec<usize> = mag_lg
                .iter()
                .copied()
                .filter(|&i| !mag_ops[i].time_reversal)
                .collect();
            let anti: Vec<usize> = mag_lg
                .iter()
                .copied()
                .filter(|&i| mag_ops[i].time_reversal)
                .collect();

            if anti.len() <= 1 {
                continue;
            }

            let mut types = Vec::new();
            for &a0 in &anti {
                let ty = wigner_classify(
                    ir.characters(),
                    group(&unitary, &mag_seitz, &h_seitz, a0),
                    ir.k_vector(),
                );
                types.push(ty);
            }
            // Only compare successful classifications
            let oks: Vec<_> = types.iter().filter_map(|r| r.as_ref().ok()).collect();
            assert!(
                !oks.is_empty(),
                "no Wigner classification succeeded for {}: {:?}",
                ir.ml,
                types
            );
            assert!(
                oks.iter().all(|&&x| x == *oks[0]),
                "Wigner type depends on a₀ for {}: {:?}",
                ir.ml,
                types
            );
        }
    }

    /// Invariant: ε(I) = +1 for any Q (identity always has parity +1).
    #[test]
    fn test_spin_parity_identity_is_plus_one() {
        // Test with Q = diag(-1,1,-1) (C2y, the SG92 case)
        let q = [[-1i32, 0, 0], [0, 1, 0], [0, 0, -1]];
        let i_rot = [[1i32, 0, 0], [0, 1, 0], [0, 0, 1]];
        // G table: identity lift is (1,0,0,0) = +I
        let g_rots: [i32; 9] = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        let g_su2: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
        // H table: identity lift also (1,0,0,0)
        let h_rots: [i32; 9] = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        let h_su2: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
        let parity = compute_signed_perm_spin_parity(&q, &i_rot, &g_rots, &g_su2, &h_rots, &h_su2)
            .expect("identity parity must be computable");
        assert!(
            (parity - 1.0).abs() < 0.1,
            "ε(I) must be +1, got {}",
            parity
        );
    }

    /// Invariant: U_Q · U_I · U_Q⁻¹ = U_I for any Q.
    #[test]
    fn test_spin_parity_conjugating_identity_gives_identity() {
        // Q = C2y = diag(-1,1,-1): 180° around y → quaternion (0,0,1,0)
        // U_Q · (1,0,0,0) · U_Q⁻¹ = (1,0,0,0) always
        let q = [[-1i32, 0, 0], [0, 1, 0], [0, 0, -1]];
        let i_rot = [[1i32, 0, 0], [0, 1, 0], [0, 0, 1]];
        let g_rots: [i32; 9] = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        let g_su2: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
        let h_rots: [i32; 9] = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        let h_su2: [f64; 4] = [1.0, 0.0, 0.0, 0.0];
        let parity = compute_signed_perm_spin_parity(&q, &i_rot, &g_rots, &g_su2, &h_rots, &h_su2)
            .expect("conjugating identity must succeed");
        assert!(
            (parity - 1.0).abs() < 0.1,
            "U_Q·U_I·U_Q⁻¹ = U_I, so ε(I) must be +1"
        );
    }

    /// Invariant: for Q = diag(-1,1,-1) (C2y), C2z in G frame transforms
    /// to C2z in H frame under Q·C2z·Q⁻¹.  With both tables using -k for
    /// C2z, the parity captures whether Q flips the SU(2) lift sign.
    ///
    /// Q=C2y lifts as j=(0,0,1,0).  j·(-k)·conj(j) = +k.  Since H table
    /// also uses -k for C2z, the two differ by sign → ε = -1.
    #[test]
    fn test_spin_parity_c2z_under_c2y() {
        let q = [[-1i32, 0, 0], [0, 1, 0], [0, 0, -1]];
        let c2z = [[-1i32, 0, 0], [0, -1, 0], [0, 0, 1]];
        let g_rots: [i32; 18] = [1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 1];
        let g_su2: [f64; 8] = [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, // -k
        ];
        let h_rots: [i32; 18] = [1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 1];
        let h_su2: [f64; 8] = [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, // -k
        ];
        let parity = compute_signed_perm_spin_parity(&q, &c2z, &g_rots, &g_su2, &h_rots, &h_su2)
            .expect("C2z parity must be computable");
        // j·(-k)·conj(j) = +k.  H table = -k.  Parity = -1.
        assert!(
            (parity - (-1.0)).abs() < 0.1,
            "C2z under C2y: j*(-k)*conj(j)=+k, H_table=-k, parity=-1. Got {}",
            parity
        );
    }

    /// Test: G table C2z = -k, H table C2z = +k.
    /// Q-transformed G lift: j*(-k)*conj(j) = +k.  H lift = +k.
    /// Same lift → ε = +1.
    #[test]
    fn test_spin_parity_c2z_opposite_lift() {
        let q = [[-1i32, 0, 0], [0, 1, 0], [0, 0, -1]];
        let c2z = [[-1i32, 0, 0], [0, -1, 0], [0, 0, 1]];
        let g_rots: [i32; 18] = [1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 1];
        let g_su2: [f64; 8] = [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, // -k
        ];
        // H table: C2z = +k
        let h_rots: [i32; 18] = [1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 1];
        let h_su2: [f64; 8] = [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, // +k
        ];
        let parity = compute_signed_perm_spin_parity(&q, &c2z, &g_rots, &g_su2, &h_rots, &h_su2)
            .expect("C2z parity must be computable");
        // j*(-k)*conj(j) = +k.  H table = +k.  Same → ε = +1.
        assert!(
            (parity - 1.0).abs() < 0.1,
            "C2z opposite lift: Q-transformed=+k, H_table=+k, parity=+1. Got {}",
            parity
        );
    }
}
