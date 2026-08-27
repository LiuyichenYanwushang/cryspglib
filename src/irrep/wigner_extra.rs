// Included by wigner.rs via include! — has access to all its imports and items.

/// Compute Type A anti-unitary characters for 1D irreps.
pub fn type_a_antiunitary_chars(
    mag_seitz: &[SeitzOp],
    mag_lg_indices: &[usize],
    h_chars: &[f64],
    h_seitz: &[SeitzOp],
    a0_idx: usize,
    k_vector: KVector,
) -> Option<(Vec<f64>, Complex64)> {
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    if a0_idx >= mag_seitz.len()
        || !mag_seitz[a0_idx].timerev
        || mag_lg_indices.iter().any(|&index| index >= mag_seitz.len())
    {
        return None;
    }
    let h_dim = h_chars.first().map(|&c| c.round() as usize).unwrap_or(1);
    if h_dim != 1 {
        return None;
    }
    let a0 = &mag_seitz[a0_idx];
    let g0 = SeitzOp::new(a0.rot, a0.trans, false);
    let (g0_sq, lattice_sq) = square_seitz(&g0);
    let m = find_seitz(&g0_sq.rot, &g0_sq.trans, h_seitz)?;
    let phase = bloch_phase(kx, ky, kz, kd, &lattice_sq);
    let chi_a0_sq = Complex64::new(*h_chars.get(m.op_index)?, 0.0) * phase;
    let u_val = if (chi_a0_sq.re - 1.0).abs() < 1e-6 {
        Complex64::new(1.0, 0.0)
    } else if (chi_a0_sq.re + 1.0).abs() < 1e-6 {
        Complex64::new(0.0, 1.0)
    } else {
        return None;
    };
    let mut au_chars = vec![0.0f64; mag_lg_indices.len()];
    for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
        if !mag_seitz[mag_idx].timerev {
            continue;
        }
        let mop = &mag_seitz[mag_idx];
        let mop_spatial = SeitzOp::new(mop.rot, mop.trans, false);
        if let Some(m) = find_seitz(&mop_spatial.rot, &mop_spatial.trans, h_seitz)
            && m.op_index < h_chars.len() {
                au_chars[out_idx] = (u_val * h_chars[m.op_index]).re;
            }
    }
    Some((au_chars, u_val))
}
