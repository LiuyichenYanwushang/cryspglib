// Included by wigner.rs via include! — has access to all its imports and items.

/// Compute Type A anti-unitary characters for 1D irreps.
pub(crate) fn type_a_antiunitary_chars(
    mag_seitz: &[SeitzOp],
    mag_lg_indices: &[usize],
    h_chars: &[f64],
    h_seitz: &[SeitzOp],
    h_dim: usize,
    a0_idx: usize,
    k_vector: KVector,
) -> Option<(Vec<f64>, Complex64)> {
    let [kx, ky, kz] = k_vector.numerators;
    let kd = k_vector.denominator;
    if h_dim != 1
        || h_chars.is_empty()
        || h_chars.len() != h_seitz.len()
        || a0_idx >= mag_seitz.len()
        || !mag_lg_indices.contains(&a0_idx)
        || mag_lg_indices.iter().any(|&index| index >= mag_seitz.len())
    {
        return None;
    }
    let a0 = &mag_seitz[a0_idx];
    if !a0.timerev
        || a0.rot != [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        || a0
            .trans
            .iter()
            .any(|translation| translation.abs() >= 1.0e-8)
    {
        return None;
    }
    let g0 = SeitzOp::new(a0.rot, a0.trans, false);
    let (g0_sq, lattice_sq) = square_seitz(&g0);
    let m = find_seitz(&g0_sq.rot, &g0_sq.trans, h_seitz)?;
    if m.lattice_shift != [0, 0, 0] {
        return None;
    }
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
        let mop = mag_seitz.get(mag_idx)?;
        if !mop.timerev {
            continue;
        }
        let mop_spatial = SeitzOp::new(mop.rot, mop.trans, false);
        let m = find_seitz(&mop_spatial.rot, &mop_spatial.trans, h_seitz)?;
        if m.lattice_shift != [0, 0, 0] {
            return None;
        }
        au_chars[out_idx] = (u_val * *h_chars.get(m.op_index)?).re;
    }
    Some((au_chars, u_val))
}
