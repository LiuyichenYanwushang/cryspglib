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
use crate::SymmetryOps;
use super::corep::CorepType;

// ── Diagnostic counters for SU(2) central-element relation ──────────────────

use std::sync::atomic::{AtomicUsize, Ordering};

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
fn find_spin_in_db(op: &SeitzOp, spin_seitz: &[SeitzOp]) -> Option<(usize, bool)> {
    // 1. Full Seitz match
    if let Some(idx) = spin_seitz.iter().position(|s| same_seitz_mod_lattice(op, s)) {
        return Some((idx, false));
    }
    // 2. Unique rotation match
    let rot_matches: Vec<usize> = spin_seitz.iter().enumerate()
        .filter_map(|(i, s)| if s.rot == op.rot { Some(i) } else { None }).collect();
    if rot_matches.len() == 1 {
        return Some((rot_matches[0], false));
    }
    // 3. -R fallback
    let r_minus: Mat3I = [
        [-op.rot[0][0], -op.rot[0][1], -op.rot[0][2]],
        [-op.rot[1][0], -op.rot[1][1], -op.rot[1][2]],
        [-op.rot[2][0], -op.rot[2][1], -op.rot[2][2]],
    ];
    let minus_matches: Vec<usize> = spin_seitz.iter().enumerate()
        .filter_map(|(i, s)| if s.rot == r_minus { Some(i) } else { None }).collect();
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
        if sop.rot != id_rot { continue; }
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
            if (eta - 1.0).abs() < 1e-6 { Some(1.0) }
            else if (eta + 1.0).abs() < 1e-6 { Some(-1.0) }
            else { None }
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
    ops: &SymmetryOps,
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
            let r = &ops[i].rotation;
            let rx = r[0][0] as i32 * kx_i + r[0][1] as i32 * ky_i + r[0][2] as i32 * kz_i;
            let ry = r[1][0] as i32 * kx_i + r[1][1] as i32 * ky_i + r[1][2] as i32 * kz_i;
            let rz = r[2][0] as i32 * kx_i + r[2][1] as i32 * ky_i + r[2][2] as i32 * kz_i;

            if ops[i].time_reversal {
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
                debug_log!("    wigner: h[{}]→H[{}] sq=H[{}] L={:?} ph={:.2} χ={:.2} → {:.2}",
                    h_mag_idx, "?", m.op_index, total_lattice, phase.re,
                    h_chars[m.op_index], contrib);
            }
        } else {
            debug_log!("    wigner: h[{}] sq R=[{},{},{};...] t=({:.3},{:.3},{:.3}) NOT FOUND",
                h_mag_idx, sq.rot[0][0],sq.rot[0][1],sq.rot[0][2],
                sq.trans[0],sq.trans[1],sq.trans[2]);
        }
    }

    let n = (unitary_mag_indices.len() as f64).max(1.0);
    let w = w_sum / n;

    // Strict classification: W must be quantized to 0, +1, or -1.
    debug_log!("DEBUG wigner_classify: w_sum={:.4} n_unitary={} W={:.4} k=({},{},{})/{}",
        w_sum, unitary_mag_indices.len(), w, kx, ky, kz, kd);
    let tol = 1e-6;
    if (w - 1.0).abs() < tol {
        CorepType::A
    } else if (w + 1.0).abs() < tol {
        CorepType::B
    } else if w.abs() < tol {
        CorepType::C
    } else {
        debug_log!("  Non-quantized Wigner indicator W={:.8}; expected 0, +1, or -1.", w);
        CorepType::Unsupported
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
            debug_log!("    cir: h[{}]→H[{}] Lz_par={} ph={:.2} χ={:.2} → {:.2}",
                h_mag_idx, m.op_index,
                ((total_lattice[2] % 2) + 2) % 2,
                phase, chi, chi * phase);
        }
    }

    debug_log!("    phase stats: +={} -={}", n_plus, n_minus);
    let n = (unitary_mag_indices.len() as f64).max(1.0);
    let w = w_sum / n;
    debug_log!("DEBUG wigner_classify_cir: W=({:.8},{:.8}) |W|={:.4} k=({},{},{})/{}",
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
    debug_log!("=== Character order check: {} ===", label);
    for (i, op) in h_seitz.iter().enumerate() {
        let re = cir_chars.get(2 * i).copied().unwrap_or(999.0);
        let im = cir_chars.get(2 * i + 1).copied().unwrap_or(999.0);
        let is_id = op.rot[0][0] == 1 && op.rot[0][1] == 0 && op.rot[0][2] == 0
                 && op.rot[1][0] == 0 && op.rot[1][1] == 1 && op.rot[1][2] == 0
                 && op.rot[2][0] == 0 && op.rot[2][1] == 0 && op.rot[2][2] == 1
                 && op.trans[0].abs() < 0.01 && op.trans[1].abs() < 0.01 && op.trans[2].abs() < 0.01;
        debug_log!("  H[{}]: R=[{},{},{};{},{},{};{},{},{}] t=({:.3},{:.3},{:.3}) chi=({:.3},{:.3}){}",
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

    debug_log!("=== unwrapped square: h[{}] ===", h_mag_idx);
    debug_log!("  a0: R={:?}, t={:?}", a0.rot, a0.trans);
    debug_log!("  h : R={:?}, t={:?}", h.rot, h.trans);
    debug_log!("  g0h raw: R={:?}, t={:?}", rc, tc_raw);
    debug_log!("  sq raw : R={:?}, t={:?}", rsq, tsq_raw);

    // Normalize for matching
    let (tsq_mod, l_reduce) = reduce01_with_lattice(&tsq_raw);
    debug_log!("  sq mod : t={:?}, L_reduce={:?}", tsq_mod, l_reduce);

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

        debug_log!("  matched H[{}]: t_stored={:?}", m.op_index, stored_t);
        debug_log!("  L_direct={:?} Lz_par={} phase={:.2}", l_direct, lz_par, phase);
        debug_log!("  m.lattice_shift={:?} (from normalized match)", m.lattice_shift);
    } else {
        debug_log!("  NOT FOUND in h_seitz");
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

        debug_log!("  direct: b[{}]^2→H[{}] L={:?} ph={:.2} χ={:.2} → {:.2}",
            b_idx, m.op_index, total_lattice, phase, chi, contrib);
    }
    let w = sum / (anti_lg_indices.len() as f64);
    debug_log!("  direct anti stats: +={} -={} W={:.4}", n_plus, n_minus, w);
    w
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
pub fn wigner_classify_spinor_direct_anti(
    ctx: &SpinLiftContext,
    spin_chars: &[f64],
    spin_lg_op_indices: &[u16],
    anti_lg_indices: &[usize],
    mag_seitz: &[SeitzOp],
    kx: i8, ky: i8, kz: i8, kd: i8,
) -> Option<CorepType> {
    let (h_spin_rots, h_spin_trans, h_spin_su2) = ctx.h;
    let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;

    if h_spin_rots.is_empty() || h_spin_su2.is_empty() || anti_lg_indices.is_empty() {
        return None;
    }

    let h_spin_seitz = build_spin_seitz(h_spin_rots, h_spin_trans);
    let g_spin_seitz = build_spin_seitz(g_spin_rots, g_spin_trans);
    if h_spin_seitz.is_empty() { return None; }

    let global_to_local: std::collections::HashMap<usize, usize> = spin_lg_op_indices
        .iter().enumerate()
        .map(|(l, &g)| (g as usize, l))
        .collect();

    // ISOTROPY origin shift: spglib → Bilbao (same as in wigner_classify_spinor).
    let (_, origin) = super::types::IrrepRecord::sg_setting(ctx.sg);
    let to_bilbao = |rot: Mat3I, trans: [f64; 3]| -> [f64; 3] {
        if origin.len() < 3 { return trans; }
        let mut t = trans;
        for i in 0..3 {
            let d: f64 = (0..3).map(|j| {
                let delta = if i == j { 1.0 } else { 0.0 };
                (delta - rot[i][j] as f64) * origin[j]
            }).sum();
            t[i] = (t[i] - d) % 1.0;
            if t[i] < 0.0 { t[i] += 1.0; }
        }
        t
    };

    let n_anti = anti_lg_indices.len();
    let mut w_sum = Complex64::ZERO;

    for &b_idx in anti_lg_indices {
        let b = &mag_seitz[b_idx];

        // Convert b to Bilbao setting, then square.
        let b_bilbao = SeitzOp::new(b.rot, to_bilbao(b.rot, b.trans), false);
        let (sq, lattice_sq) = square_seitz(&b_bilbao);

        // b² ∈ H₀ by group theory: b ∈ M_k ⇒ b² ∈ H_k ⇒ R_{b²} ∈ H₀.
        // Use LG-first matching to avoid picking a non-LG candidate.
        let (sq_spin_idx, sq_in_lg) = match find_sq_spin_lg_first(
            &sq, &h_spin_seitz, spin_lg_op_indices) {
            Some(v) => v,
            None => {
                debug_log!("  SPINOR_DIRECT_ANTI fail: b[{}]² rot={:?} not in H spin ops",
                    b_idx, sq.rot);
                return None;
            }
        };

        let sq_local_idx = if sq_in_lg {
            *global_to_local.get(&sq_spin_idx)?
        } else {
            debug_log!("  SPINOR_DIRECT_ANTI fail: b[{}]² spin[{}] not in LG idxs",
                b_idx, sq_spin_idx);
            return None;
        };

        // SU(2) lift of b (rotation-only lookup in G spin ops — b may have
        // improper rotation from G \ H for Type III).
        let b_spin_idx = g_spin_seitz.iter().position(|s| s.rot == b.rot)
            .or_else(|| {
                let r: Mat3I = [
                    [-b.rot[0][0], -b.rot[0][1], -b.rot[0][2]],
                    [-b.rot[1][0], -b.rot[1][1], -b.rot[1][2]],
                    [-b.rot[2][0], -b.rot[2][1], -b.rot[2][2]],
                ];
                g_spin_seitz.iter().position(|s| s.rot == r)
            })?;
        let u_b = spin_su2_at(g_spin_su2, b_spin_idx)?;

        // SU(2) central detection: U_b² vs canonical U_{b²}.
        let u_b_sq = su2_compose(&u_b, &u_b);
        let u_sq = spin_su2_at(h_spin_su2, sq_spin_idx)?;
        let central = su2_same_up_to_sign(&u_b_sq, &u_sq)?;

        // Bloch phase from lattice shift (translation origin shift already
        // accounted for by to_bilbao).
        let phase = bloch_phase(kx, ky, kz, kd, &lattice_sq);

        let chi0 = spin_chars[sq_local_idx];
        let chi = if central { -chi0 } else { chi0 };
        w_sum += Complex64::new(chi, 0.0) * phase;
    }

    let w = w_sum / (n_anti as f64);

    // Dimension from identity canonical lift.
    let id_rot: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let u_id = [1.0, 0.0, 0.0, 0.0];
    let h_dim = spin_lg_op_indices.iter()
        .map(|&x| x as usize)
        .find_map(|si| {
            let sop = h_spin_seitz.get(si)?;
            if sop.rot != id_rot { return None; }
            let u = spin_su2_at(h_spin_su2, si)?;
            if su2_same_up_to_sign(&u, &u_id) != Some(false) { return None; }
            let local = *global_to_local.get(&si)?;
            spin_chars.get(local).map(|&c| c.abs().round().max(1.0))
        })
        .unwrap_or_else(|| {
            spin_chars.first().map(|&c| c.abs().round().max(1.0)).unwrap_or(1.0)
        });

    let tol = 1e-6;
    if (w.re - h_dim).abs() < tol && w.im.abs() < tol {
        Some(CorepType::A)
    } else if (w.re + h_dim).abs() < tol && w.im.abs() < tol {
        Some(CorepType::B)
    } else if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        Some(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        Some(CorepType::B)
    } else if w.norm() < tol {
        Some(CorepType::C)
    } else {
        None
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
        if tr[i] < 0.0 { tr[i] += 1.0; l[i] -= 1; }
    }
    (tr, l)
}

/// Build index map from `h_seitz` (spglib H_ops order) to CIR operation order.
///
/// Returns `None` if any H operation cannot be matched (e.g., missing rotation data).
/// Build mapping from H_ops to PIR/CIR operations using full Seitz matching.
///
/// Unlike [`build_h_to_cir_map`] (rotation-only), this matches by both
/// rotation AND translation (modulo lattice), eliminating ambiguity when
/// the same rotation appears with different fractional translations
/// (e.g. in nonsymmorphic groups).
pub fn build_h_to_irrep_op_map(
    h_seitz: &[SeitzOp],
    irrep_rots: &[i32],
    irrep_trans: &[f64],
) -> Option<Vec<usize>> {
    let n_ops = h_seitz.len();
    let n_ir_ops = irrep_rots.len() / 9;
    if n_ir_ops == 0 || irrep_trans.len() < n_ir_ops * 3 {
        return None;
    }
    let mut map = Vec::with_capacity(n_ops);
    for h_idx in 0..n_ops {
        let h = &h_seitz[h_idx];
        let mut found = None;
        for ir_idx in 0..n_ir_ops {
            let r_match = {
                let off = ir_idx * 9;
                off + 8 < irrep_rots.len()
                    && irrep_rots[off] == h.rot[0][0]
                    && irrep_rots[off+1] == h.rot[0][1]
                    && irrep_rots[off+2] == h.rot[0][2]
                    && irrep_rots[off+3] == h.rot[1][0]
                    && irrep_rots[off+4] == h.rot[1][1]
                    && irrep_rots[off+5] == h.rot[1][2]
                    && irrep_rots[off+6] == h.rot[2][0]
                    && irrep_rots[off+7] == h.rot[2][1]
                    && irrep_rots[off+8] == h.rot[2][2]
            };
            if !r_match { continue; }
            // Check translation modulo lattice
            let toff = ir_idx * 3;
            let t_ok = (0..3).all(|k| {
                let d = h.trans[k] - irrep_trans[toff + k];
                (d - d.round()).abs() < 1e-9
            });
            if t_ok {
                if found.is_some() {
                    return None; // ambiguous: shouldn't happen with full Seitz
                }
                found = Some(ir_idx);
            }
        }
        map.push(found?);
    }
    Some(map)
}

/// Build mapping from H_ops to PIR/CIR operations using rotation-only matching.
///
/// DEPRECATED: prefer [`build_h_to_irrep_op_map`] when translation data is available.
/// Rotation-only matching can be ambiguous for nonsymmorphic groups.
pub fn build_h_to_cir_map(h_seitz: &[SeitzOp], cir_rots: &[i32]) -> Option<Vec<usize>> {
    let n_ops = h_seitz.len();
    let n_cir_ops = cir_rots.len() / 9;
    if n_cir_ops == 0 {
        return None;
    }
    let mut map = vec![0usize; n_ops];

    debug_log!("build_h_to_cir_map: n_ops={} n_cir_ops={}", n_ops, n_cir_ops);
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
                debug_log!("build_h_to_cir_map: H[{}] R=[{},{},{};{},{},{};{},{},{}] not found in rots ({} ops)",
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

/// Find the spin op index for a computed Seitz square `sq`, preferring
/// candidates **inside** `spin_lg_op_indices` over the full database.
///
/// Priority:
/// 1. Full Seitz match (rotation + translation mod lattice) inside LG candidates
/// 2. Unique rotation-only match inside LG candidates
/// 3. Rotation-only match in full h_spin_seitz (fallback, may be outside LG)
pub(crate) fn find_sq_spin_lg_first(
    sq: &SeitzOp,
    h_spin_seitz: &[SeitzOp],
    spin_lg_op_indices: &[u16],
) -> Option<(usize, bool)> {
    let lg_cands: Vec<usize> = spin_lg_op_indices.iter().map(|&x| x as usize).collect();

    // 1. Full Seitz match inside LG
    for &si in &lg_cands {
        if let Some(sop) = h_spin_seitz.get(si) {
            if same_seitz_mod_lattice(sq, sop) {
                return Some((si, true));
            }
        }
    }

    // 2. Unique rotation-only match inside LG
    let rot_lg: Vec<usize> = lg_cands.iter().copied()
        .filter(|&si| h_spin_seitz.get(si).map_or(false, |s| s.rot == sq.rot))
        .collect();
    if rot_lg.len() == 1 {
        return Some((rot_lg[0], true));
    }

    // 3. Fallback: rotation-only in full database
    let global_idx = h_spin_seitz.iter().position(|s| s.rot == sq.rot)?;
    let in_lg = lg_cands.contains(&global_idx);
    Some((global_idx, in_lg))
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
pub fn su2_same_up_to_sign(a: &[f64; 4], b: &[f64; 4]) -> Option<bool> {
    let dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
    let na = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]).sqrt();
    let nb = (b[0]*b[0] + b[1]*b[1] + b[2]*b[2] + b[3]*b[3]).sqrt();
    if na < 1e-10 || nb < 1e-10 {
        return None;
    }
    let cos = dot / (na * nb);
    if (cos - 1.0).abs() < 1e-6 {
        Some(false)  // same lift
    } else if (cos + 1.0).abs() < 1e-6 {
        Some(true)   // opposite lift (central Ē)
    } else {
        None  // unrelated
    }
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
        [spin_op_rots[9 * idx + 0], spin_op_rots[9 * idx + 1], spin_op_rots[9 * idx + 2]],
        [spin_op_rots[9 * idx + 3], spin_op_rots[9 * idx + 4], spin_op_rots[9 * idx + 5]],
        [spin_op_rots[9 * idx + 6], spin_op_rots[9 * idx + 7], spin_op_rots[9 * idx + 8]],
    ];
    let t = [spin_op_trans[3 * idx], spin_op_trans[3 * idx + 1], spin_op_trans[3 * idx + 2]];
    Some(SeitzOp::new(r, t, false))
}

/// Extract Pauli coefficients [u₀,u₁,u₂,u₃] from the spin-op flat array.
pub fn spin_su2_at(spin_op_su2: &[f64], idx: usize) -> Option<[f64; 4]> {
    if 4 * idx + 3 >= spin_op_su2.len() {
        return None;
    }
    Some([
        spin_op_su2[4 * idx + 0],
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
        if (d - d.round()).abs() > 1e-9 {
            return false;
        }
    }
    true
}

/// Build mapping from H_ops (spglib order) to spin-op index.
///
/// Strategy:
/// 1. Full Seitz matching (rotation + translation via `same_seitz_mod_lattice`)
///    — disambiguates g vs Ēg lifts at BZ-boundary k-points.
/// 2. Fallback: rotation-only, but only when exactly one candidate exists in
///    `spin_lg_op_indices`.  Multiple candidates with the same rotation but
///    different translations are UNRESOLVED → return None (caller skips).
///
/// Uses a `Vec` (not HashSet) for deterministic iteration order.
pub fn build_h_to_spin_map(
    h_seitz: &[SeitzOp],
    spin_seitz: &[SeitzOp],
    spin_lg_op_indices: &[u16],
) -> Vec<Option<usize>> {
    let allowed: Vec<usize> = spin_lg_op_indices.iter().map(|&x| x as usize).collect();

    let mut map = Vec::with_capacity(h_seitz.len());
    for h in h_seitz {
        // 1. Full Seitz matching
        let full = allowed.iter().copied().find(|&si| {
            spin_seitz.get(si).map_or(false, |sop| same_seitz_mod_lattice(h, sop))
        });
        if let Some(si) = full { map.push(Some(si)); continue; }

        // 2. Rotation-only fallback — only when uniquely resolved
        let rots: Vec<usize> = allowed.iter().copied().filter(|&si| {
            spin_seitz.get(si).map_or(false, |sop| sop.rot == h.rot)
        }).collect();

        if rots.len() == 1 {
            map.push(Some(rots[0]));
        } else {
            map.push(None);
        }
    }
    map
}

/// Build a Vec<SeitzOp> from the spin-op flat arrays (public for testing).
pub fn build_spin_seitz(spin_op_rots: &[i32], spin_op_trans: &[f64]) -> Vec<SeitzOp> {
    let n = (spin_op_rots.len() / 9).min(spin_op_trans.len() / 3);
    (0..n)
        .filter_map(|i| spin_seitz_at(i, spin_op_rots, spin_op_trans))
        .collect()
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
pub fn diagnostic_extra_sum(extra: &[f64]) -> f64 {
    extra.iter().sum()
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
    spin_chars: &[f64],
    n_lg_ops: usize,
    spin_lg_op_indices: &[u16],
    _unitary_mag_indices: &[usize],
    mag_seitz: &[SeitzOp],
    h_seitz: &[SeitzOp],
    a0_idx: usize,
    kx: i8, ky: i8, kz: i8, kd: i8,
) -> Option<CorepType> {
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
        return None;
    }

    // Spin Seitz ops in Bilbao setting — canonical little co-group representatives.
    let h_spin_seitz = build_spin_seitz(h_spin_rots, h_spin_trans);
    if h_spin_seitz.is_empty() {
        return None;
    }

    // H_op → spin global index mapping (for matching (a₀h)² back to spin ops)
    let h_to_spin = build_h_to_spin_map(h_seitz, &h_spin_seitz, spin_lg_op_indices);

    // global spin op index → local character table position
    let global_to_local: std::collections::HashMap<usize, usize> = spin_lg_op_indices
        .iter()
        .enumerate()
        .map(|(local, &global)| (global as usize, local))
        .collect();

    // Infer central parity eta_ebar = chi(Ebar)/chi(E) from the character table.
    // For genuine spinor irreps: -1.0.  For single-valued: +1.0.
    let eta_ebar = infer_eta_ebar(
        spin_chars, spin_lg_op_indices, &h_spin_seitz, h_spin_su2,
    ).unwrap_or(-1.0);

    let a0 = &mag_seitz[a0_idx];

    // a₀ SU(2) lift: rotation-only lookup in G's spin ops.
    let g_spin_seitz = build_spin_seitz(g_spin_rots, g_spin_trans);
    let a0_match = g_spin_seitz.iter()
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
        })?;
    let u_a0 = spin_su2_at(g_spin_su2, a0_match)?;

    // ISOTROPY origin shift: spglib → Bilbao.  t_B = t_S - (I - R)·origin.
    // a₀ comes from spglib/MSG → needs conversion.
    // h is from spin_seitz (already Bilbao) → no conversion needed.
    let (_, origin) = super::types::IrrepRecord::sg_setting(ctx.sg);
    let to_bilbao = |rot: Mat3I, trans: [f64; 3]| -> [f64; 3] {
        if origin.len() < 3 { return trans; }
        let mut t = trans;
        for i in 0..3 {
            let d: f64 = (0..3).map(|j| {
                let delta = if i == j { 1.0 } else { 0.0 };
                (delta - rot[i][j] as f64) * origin[j]
            }).sum();
            t[i] = (t[i] - d) % 1.0;
            if t[i] < 0.0 { t[i] += 1.0; }
        }
        t
    };

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

    for local in 0..n_lg_ops {
        let global_spin_idx = spin_lg_op_indices[local] as usize;

        // Canonical h in Bilbao (already in the correct setting).
        let h_spin = &h_spin_seitz[global_spin_idx];
        let u_h = spin_su2_at(h_spin_su2, global_spin_idx)?;

        // Spatial: (a₀ h)² in Bilbao.
        let (g0h, l1) = compose_seitz(&a0_bilbao, h_spin);
        let (sq, lattice_sq) = square_seitz(&g0h);

        // Match square's rotation back to spin ops, preferring LG candidates first.
        // Priority: full Seitz in LG → unique rotation in LG → global rotation.
        // This avoids position() picking a non-LG candidate when an LG candidate exists.
        let (sq_spin_idx, sq_in_lg) = match find_sq_spin_lg_first(
            &sq, &h_spin_seitz, spin_lg_op_indices) {
            Some(v) => v,
            None => {
                eprintln!("  WIGNER_SPINOR: sq_rot not in spin ops, aborting case");
                return None;
            }
        };

        // SU(2): (U_a₀·U_h)² vs canonical U_{(a₀h)²}.
        let u_g0h = su2_compose(&u_a0, &u_h);
        let u_sq = SquareKernel::OldU2.apply(&u_g0h);
        let u_k = spin_su2_at(h_spin_su2, sq_spin_idx)?;

        // Central element detection.
        // central=true: u_sq ≈ -u_k (differs by Ebar)
        // central=false: u_sq ≈ u_k (same lift)
        let central = match su2_same_up_to_sign(&u_sq, &u_k) {
            Some(false) => { SU2_REL_SAME.fetch_add(1, Ordering::Relaxed); false }
            Some(true) => { SU2_REL_EBAR.fetch_add(1, Ordering::Relaxed); true }
            None => {
                SU2_REL_NONE.fetch_add(1, Ordering::Relaxed);
                // Scan same-rotation candidates
                let lg_set: std::collections::HashSet<usize> =
                    spin_lg_op_indices.iter().map(|&x| x as usize).collect();
                let mut matched_other_lg = false;
                let mut matched_other_global = false;
                let mut has_cand = false;
                for (ci, cs) in h_spin_seitz.iter().enumerate() {
                    if cs.rot != sq.rot || ci == sq_spin_idx { continue; }
                    has_cand = true;
                    if let Some(uc) = spin_su2_at(h_spin_su2, ci) {
                        if su2_same_up_to_sign(&u_sq, &uc).is_some() {
                            if lg_set.contains(&ci) { matched_other_lg = true; }
                            else { matched_other_global = true; }
                        }
                    }
                }
                if matched_other_lg { NONE_MATCH_OTHER_LG.fetch_add(1, Ordering::Relaxed); }
                else if matched_other_global { NONE_MATCH_OTHER_GLOBAL.fetch_add(1, Ordering::Relaxed); }
                else if has_cand { NONE_NO_MATCH_HAS_CAND.fetch_add(1, Ordering::Relaxed); }
                else { NONE_NO_CANDIDATE.fetch_add(1, Ordering::Relaxed); }
                // Det distribution
                let det_a0 = a0.rot[0][0]*(a0.rot[1][1]*a0.rot[2][2]-a0.rot[1][2]*a0.rot[2][1])
                    - a0.rot[0][1]*(a0.rot[1][0]*a0.rot[2][2]-a0.rot[1][2]*a0.rot[2][0])
                    + a0.rot[0][2]*(a0.rot[1][0]*a0.rot[2][1]-a0.rot[1][1]*a0.rot[2][0]);
                if det_a0 > 0 { NONE_DET_A0_P1.fetch_add(1, Ordering::Relaxed); }
                else { NONE_DET_A0_M1.fetch_add(1, Ordering::Relaxed); }
                let det_g0h = g0h.rot[0][0]*(g0h.rot[1][1]*g0h.rot[2][2]-g0h.rot[1][2]*g0h.rot[2][1])
                    - g0h.rot[0][1]*(g0h.rot[1][0]*g0h.rot[2][2]-g0h.rot[1][2]*g0h.rot[2][0])
                    + g0h.rot[0][2]*(g0h.rot[1][0]*g0h.rot[2][1]-g0h.rot[1][1]*g0h.rot[2][0]);
                if det_g0h > 0 { NONE_DET_G0H_P1.fetch_add(1, Ordering::Relaxed); }
                else { NONE_DET_G0H_M1.fetch_add(1, Ordering::Relaxed); }
                // Test 6 alternative antiunitary square formulas
                let u_cj = conj_pauli(&u_g0h);
                let alts = [
                    su2_compose(&u_g0h, &u_g0h),           // U^2 (raw)
                    neg_pauli(&su2_compose(&u_g0h, &u_g0h)), // -U^2
                    su2_compose(&u_g0h, &u_cj),             // U U*
                    neg_pauli(&su2_compose(&u_g0h, &u_cj)), // -U U*
                    su2_compose(&u_cj, &u_g0h),             // U* U
                    neg_pauli(&su2_compose(&u_cj, &u_g0h)), // -U* U
                ];
                let matches: Vec<bool> = alts.iter()
                    .map(|a| su2_same_up_to_sign(a, &u_k).is_some()).collect();
                if matches[0] { NONE_ALT_RAW.fetch_add(1, Ordering::Relaxed); }
                if matches[1] { NONE_ALT_NEG_RAW.fetch_add(1, Ordering::Relaxed); }
                if matches[2] { NONE_ALT_UUSTAR.fetch_add(1, Ordering::Relaxed); }
                if matches[3] { NONE_ALT_NEG_UUSTAR.fetch_add(1, Ordering::Relaxed); }
                if matches[4] { NONE_ALT_STARU.fetch_add(1, Ordering::Relaxed); }
                if matches[5] { NONE_ALT_NEG_STARU.fetch_add(1, Ordering::Relaxed); }
                if !matches.iter().any(|&m| m) { NONE_ALT_NONE.fetch_add(1, Ordering::Relaxed); }
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
                if j_matches[0] { NONE_JU_JU_STAR.fetch_add(1, Ordering::Relaxed); }
                if j_matches[1] { NONE_NEG_JU_JU_STAR.fetch_add(1, Ordering::Relaxed); }
                if j_matches[2] { NONE_UJ_UJ_STAR.fetch_add(1, Ordering::Relaxed); }
                if j_matches[3] { NONE_NEG_UJ_UJ_STAR.fetch_add(1, Ordering::Relaxed); }
                if !j_matches.iter().any(|&m| m) { NONE_J_NONE.fetch_add(1, Ordering::Relaxed); }
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
                            } else { GGAUGE_NONE.fetch_add(1, Ordering::Relaxed); }
                        } else { GGAUGE_NONE.fetch_add(1, Ordering::Relaxed); }
                    } else { GGAUGE_SQ_LOOKUP_FAIL.fetch_add(1, Ordering::Relaxed); }
                } else { GGAUGE_H_LOOKUP_FAIL.fetch_add(1, Ordering::Relaxed); }
                return None;
            }
        };

        // Character from LG table (if sq ∈ LG) or from extended table.
        let sq_local_idx = if sq_in_lg {
            *global_to_local.get(&sq_spin_idx)?
        } else {
            // sq outside LG: need extended character, abort case.
            eprintln!("  WIGNER_SPINOR: sq[{}] not in LG idxs, aborting case", sq_spin_idx);
            return None;
        };

        // Bloch phase from total lattice shift.
        let r_l1 = mat_vec_i32(&g0h.rot, &l1);
        let total_l = add3(&lattice_sq, &add3(&l1, &r_l1));
        let phase = bloch_phase(kx, ky, kz, kd, &total_l);

        let chi0 = spin_chars[sq_local_idx];
        let chi = if central { eta_ebar * chi0 } else { chi0 };

        w_sum += Complex64::new(chi, 0.0) * phase;
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
            spin_chars.get(local).map(|&c| c.abs().round().max(1.0))
        })
        .unwrap_or_else(|| {
            spin_chars
                .first()
                .map(|&c| c.abs().round().max(1.0))
                .unwrap_or(1.0)
        });

    // Classification.  For spinor irreps: W = +dim → Type A, W = -dim → Type B,
    // W = 0 → Type C.  Also accept W = ±1 for 1D irreps.
    let tol = 1e-6;
    if (w.re - h_dim).abs() < tol && w.im.abs() < tol {
        Some(CorepType::A)
    } else if (w.re + h_dim).abs() < tol && w.im.abs() < tol {
        Some(CorepType::B)
    } else if (w.re - 1.0).abs() < tol && w.im.abs() < tol {
        Some(CorepType::A)
    } else if (w.re + 1.0).abs() < tol && w.im.abs() < tol {
        Some(CorepType::B)
    } else if w.norm() < tol {
        Some(CorepType::C)
    } else {
        None // non-quantized → Unsupported
    }
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

    // Dimension of the conjugate irrep: χ(E)
    let conj_dim = conj_chars[0].abs();
    let norm = 1.0 / (n_ops as f64);

    // For each candidate irrep, compute overlap with conjugate
    let mut best_idx = None;
    let mut best_overlap = 0.0f64;

    for (j, chars_j) in candidate_chars.iter().enumerate() {
        if chars_j.len() < n_ops { continue; }
        // Dimension must match: Δⱼ ∼ Δᵢ^{a₀} ⇒ same dimension
        if (chars_j[0].abs() - conj_dim).abs() > 0.01 { continue; }
        let overlap: f64 = (0..n_ops)
            .map(|h| conj_chars[h] * chars_j[h])
            .sum();
        let overlap_norm = overlap.abs() * norm;

        if overlap_norm > best_overlap {
            best_overlap = overlap_norm;
            best_idx = Some(j);
        }
    }

    // Return partner if overlap is significantly above noise.
    // For n-dimensional irrep, overlap ≈ d²/|H| per operation.
    if best_overlap > 0.1 {
        best_idx
    } else {
        None
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
    mag_ops: &SymmetryOps,
    mag_lg_indices: &[usize],
    op_map: &[Option<usize>],
    h_chars: &[f64],
    partner_chars: Option<&[f64]>,  // for Type C: character table of paired irrep
    au_chars: Option<&[f64]>,       // for Type A: pre-computed antiunitary chars
) -> Vec<f64> {
    let n_lg = mag_lg_indices.len();
    let mut chars = vec![0.0; n_lg];

    for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
        let is_anti = mag_ops[mag_idx].time_reversal;
        let h_idx = op_map[mag_idx];

        match corep_type {
            CorepType::A => {
                if is_anti {
                    if let Some(ac) = au_chars {
                        if out_idx < ac.len() { chars[out_idx] = ac[out_idx]; }
                    }
                } else if let Some(hi) = h_idx {
                    if hi < h_chars.len() { chars[out_idx] = h_chars[hi]; }
                }
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
        let ops = SymmetryOps::from_parallel_owned(
            vec![[[0,-1,0],[1,0,0],[0,0,-1]],
                 [[1,0,0],[0,1,0],[0,0,1]]],
            vec![[0.0; 3], [0.0; 3]],
            vec![true, false],
        );
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
            vec![[[-1,0,0],[0,1,0],[0,0,1]],
                 [[1,0,0],[0,-1,0],[0,0,-1]]],
            vec![[0.0; 3], [0.0; 3]],
            vec![true, true],
        );
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
            let unitary: Vec<usize> = mag_lg.iter().copied()
                .filter(|&i| !mag_ops[i].time_reversal).collect();
            let anti: Vec<usize> = mag_lg.iter().copied()
                .filter(|&i| mag_ops[i].time_reversal).collect();

            if anti.len() <= 1 { continue; }

            let mut types = Vec::new();
            for &a0 in &anti {
                let ty = wigner_classify(
                    ir.characters(), &unitary, &mag_seitz, &h_seitz, a0,
                    ir.kx, ir.ky, ir.kz, ir.kd,
                );
                types.push(ty);
            }
            assert!(
                types.iter().all(|&x| x == types[0]),
                "Wigner type depends on a₀ for {}: {:?}", ir.ml, types
            );
        }
    }
}
