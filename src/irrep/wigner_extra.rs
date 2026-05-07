// Included by wigner.rs via include! — has access to all its imports and items.

/// Invert a small complex matrix (d <= 8) using Gaussian elimination.
fn mat_inv(m: &[Complex64], dim: usize) -> Option<Vec<Complex64>> {
    let n = dim;
    if n == 0 { return Some(vec![]); }
    let mut a = m.to_vec();
    let mut inv = vec![Complex64::new(0.0, 0.0); n * n];
    for i in 0..n { inv[i * n + i] = Complex64::new(1.0, 0.0); }
    for col in 0..n {
        let mut pivot = col;
        let mut max_abs = a[col * n + col].norm();
        for row in (col+1)..n {
            let v = a[row * n + col].norm();
            if v > max_abs { max_abs = v; pivot = row; }
        }
        if max_abs < 1e-14 { return None; }
        if pivot != col {
            for j in 0..n { a.swap(col * n + j, pivot * n + j); inv.swap(col * n + j, pivot * n + j); }
        }
        let piv_val = a[col * n + col];
        for j in 0..n { a[col * n + j] = a[col * n + j] / piv_val; inv[col * n + j] = inv[col * n + j] / piv_val; }
        for row in 0..n {
            if row == col { continue; }
            let factor = a[row * n + col];
            if factor.norm() < 1e-14 { continue; }
            for j in 0..n {
                a[row * n + j] = a[row * n + j] - factor * a[col * n + j];
                inv[row * n + j] = inv[row * n + j] - factor * inv[col * n + j];
            }
        }
    }
    Some(inv)
}

fn mat_mul(a: &[Complex64], b: &[Complex64], dim: usize) -> Vec<Complex64> {
    let n = dim;
    let mut c = vec![Complex64::new(0.0, 0.0); n * n];
    for i in 0..n { for k in 0..n { let aik = a[i * n + k]; if aik.norm() < 1e-14 { continue; } for j in 0..n { c[i * n + j] = c[i * n + j] + aik * b[k * n + j]; } } }
    c
}

fn mat_conj(m: &[Complex64]) -> Vec<Complex64> { m.iter().map(|c| c.conj()).collect() }

/// Compute consquare root U such that U*conj(U) = M (Type A intertwiner).
fn consquare_root(m: &[Complex64], dim: usize) -> Option<Vec<Complex64>> {
    let n = dim; let n2 = n * n;
    let mut u = vec![Complex64::new(0.0, 0.0); n2];
    for i in 0..n { u[i * n + i] = Complex64::new(1.0, 0.0); }
    for _iter in 0..80 {
        let u_conj = mat_conj(&u);
        let u_conj_inv = mat_inv(&u_conj, n)?;
        let m_ucinv = mat_mul(m, &u_conj_inv, n);
        for i in 0..n2 { u[i] = Complex64::new((u[i].re + m_ucinv[i].re) * 0.5, (u[i].im + m_ucinv[i].im) * 0.5); }
        let u_conj = mat_conj(&u);
        let uuc = mat_mul(&u, &u_conj, n);
        let mut err = 0.0f64;
        for i in 0..n2 { err += (uuc[i] - m[i]).norm_sqr(); }
        if err < 1e-12 { return Some(u); }
    }
    None
}

/// Compute Type A anti-unitary characters for 1D irreps.
pub fn type_a_antiunitary_chars(
    mag_seitz: &[SeitzOp], mag_lg_indices: &[usize],
    _op_map: &[Option<usize>], h_chars: &[f64], h_seitz: &[SeitzOp],
    a0_idx: usize, kx: i8, ky: i8, kz: i8, kd: i8,
) -> Option<(Vec<f64>, Complex64)> {
    let h_dim = h_chars.first().map(|&c| c.round() as usize).unwrap_or(1);
    if h_dim != 1 { return None; }
    let a0 = &mag_seitz[a0_idx];
    let g0 = SeitzOp::new(a0.rot, a0.trans, false);
    let (g0_sq, lattice_sq) = square_seitz(&g0);
    let m = find_seitz(&g0_sq.rot, &g0_sq.trans, h_seitz)?;
    let phase = bloch_phase(kx, ky, kz, kd, &lattice_sq);
    let chi_a0_sq = Complex64::new(h_chars[m.op_index], 0.0) * phase;
    let u_val = if (chi_a0_sq.re - 1.0).abs() < 1e-6 { Complex64::new(1.0, 0.0) }
    else if (chi_a0_sq.re + 1.0).abs() < 1e-6 { Complex64::new(0.0, 1.0) }
    else { return None; };
    let mut au_chars = vec![0.0f64; mag_lg_indices.len()];
    for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
        if !mag_seitz[mag_idx].timerev { continue; }
        let mop = &mag_seitz[mag_idx];
        let mop_spatial = SeitzOp::new(mop.rot, mop.trans, false);
        if let Some(m) = find_seitz(&mop_spatial.rot, &mop_spatial.trans, h_seitz) {
            if m.op_index < h_chars.len() { au_chars[out_idx] = (u_val * h_chars[m.op_index]).re; }
        }
    }
    Some((au_chars, u_val))
}

/// Compute Type A anti-unitary characters for high-dimensional (>1D) irreps.
pub fn type_a_antiunitary_chars_high_dim(
    mag_seitz: &[SeitzOp], mag_lg_indices: &[usize],
    h_chars: &[f64], h_seitz: &[SeitzOp],
    a0_idx: usize, _kx: i8, _ky: i8, _kz: i8, _kd: i8,
    pir_mats: &[f64], pir_rots: &[i32],
) -> Option<Vec<f64>> {
    let h_dim = h_chars.first().map(|&c| c.round() as usize).unwrap_or(1);
    if h_dim <= 1 { return None; }
    let block = h_dim * h_dim;
    let n_pir_ops = pir_rots.len() / 9;
    if n_pir_ops == 0 || pir_mats.len() < n_pir_ops * block { return None; }
    let h_to_pir = build_h_to_cir_map(h_seitz, pir_rots)?;
    let a0 = &mag_seitz[a0_idx];
    let g0 = SeitzOp::new(a0.rot, a0.trans, false);
    let (g0_sq, _lattice_sq) = square_seitz(&g0);
    let a0sq_op = find_seitz(&g0_sq.rot, &g0_sq.trans, h_seitz)?;
    let a0sq_hpos = h_to_pir[a0sq_op.op_index];
    if a0sq_hpos >= n_pir_ops { return None; }
    let a0sq_mat = &pir_mats[a0sq_hpos * block..(a0sq_hpos + 1) * block];
    let a0sq_complex: Vec<Complex64> = a0sq_mat.iter().map(|&v| Complex64::new(v, 0.0)).collect();

    // D(a0^2) = -I: antiunitary chars = 0 (analogous to 1D U=i case)
    let mut is_minus_id = true;
    for i in 0..h_dim { for j in 0..h_dim {
        let v = a0sq_mat[i * h_dim + j];
        let expected = if i == j { -1.0 } else { 0.0 };
        if (v - expected).abs() > 1e-10 { is_minus_id = false; break; }
    }}
    if is_minus_id {
        let mut au_chars = vec![0.0f64; mag_lg_indices.len()];
        for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
            if !mag_seitz[mag_idx].timerev { continue; }
            au_chars[out_idx] = 0.0;
        }
        return Some(au_chars);
    }

    let u = consquare_root(&a0sq_complex, h_dim)?;
    let mut au_chars = vec![0.0f64; mag_lg_indices.len()];
    for (out_idx, &mag_idx) in mag_lg_indices.iter().enumerate() {
        if !mag_seitz[mag_idx].timerev { continue; }
        let mop = &mag_seitz[mag_idx];
        let mop_spatial = SeitzOp::new(mop.rot, mop.trans, false);
        if let Some(m) = find_seitz(&mop_spatial.rot, &mop_spatial.trans, h_seitz) {
            let hpos = h_to_pir[m.op_index];
            if hpos < n_pir_ops {
                let h_mat = &pir_mats[hpos * block..(hpos + 1) * block];
                let mut tr = Complex64::new(0.0, 0.0);
                for i in 0..h_dim { for j in 0..h_dim {
                    tr = tr + u[i * h_dim + j] * Complex64::new(h_mat[j * h_dim + i], 0.0);
                }}
                au_chars[out_idx] = tr.re;
            }
        }
    }
    Some(au_chars)
}
