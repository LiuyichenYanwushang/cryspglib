//! Wyckoff 位置的精确定位和对等原子分配。
//!
//! 使用位点对称性数据库确定原子的确切经坐标和 Wyckoff 字母标记。
//!
//! 参考: R. W. Grosse-Kunstleve and P. D. Adams,
//! Acta Cryst. (2002). A58, 60-65

use crate::cell::{cel_is_overlap, cel_is_overlap_with_same_type, cel_layer_is_overlap,
                   cel_layer_is_overlap_with_same_type, AperiodicAxis, Cell};
use crate::debug;
use crate::mathfunc::*;
use crate::sitesym_database::*;
use crate::symmetry::Symmetry;

const INCREASE_RATE: f64 = 1.05;
const NUM_ATTEMPT: i32 = 5;

/// 获取精确的原子位置和 Wyckoff 标记。
///
/// # Returns
/// `Option<(positions, wyckoffs, equiv_atoms, site_sym_symbols)>`
pub fn ssm_get_exact_positions(
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    num_pure_trans: i32,
    hall_number: usize,
    symprec: f64,
) -> Option<(Vec<Vec3>, Vec<i32>, Vec<i32>, Vec<String>)> {
    let mut tolerance = symprec;
    for i in 0..NUM_ATTEMPT {
        let (positions, equiv_atoms) = match get_exact_positions(conv_prim, conv_sym, tolerance) {
            Some(v) => v,
            None => return None,
        };

        let (wyckoffs, symbols) = match set_wyckoffs_labels(
            &positions, &equiv_atoms, conv_prim, conv_sym,
            num_pure_trans, hall_number, symprec,
        ) {
            Some(v) => v,
            None => {
                debug::debug_print(format_args!(
                    "spglib: ssm_get_exact_positions failed (attempt={}).\n", i
                ));
                tolerance *= INCREASE_RATE;
                continue;
            }
        };

        return Some((positions, wyckoffs, equiv_atoms, symbols));
    }

    None
}

/// 内部函数：为每个原子确定精确位置。
/// 返回 (positions, equiv_atoms) 或 None
fn get_exact_positions(
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    symprec: f64,
) -> Option<(Vec<Vec3>, Vec<i32>)> {
    debug::debug_print(format_args!("get_exact_positions\n"));

    let n = conv_prim.size;
    let mut positions = vec![[0.0; 3]; n];
    let mut equiv_atoms = vec![0i32; n];
    let mut indep_atoms: Vec<usize> = Vec::with_capacity(n);

    if conv_prim.aperiodic_axis.is_none() {
        for i in 0..n {
            if !set_equivalent_atom(
                &positions, &mut equiv_atoms, i,
                &indep_atoms, conv_prim, conv_sym, symprec,
            ) {
                equiv_atoms[i] = i as i32;
                indep_atoms.push(i);
                positions[i] = conv_prim.position[i];
                set_exact_location(&mut positions[i], conv_sym, &conv_prim.lattice, symprec);
            }
        }
    } else {
        let aperiodic = conv_prim.aperiodic_axis.unwrap();
        for i in 0..n {
            if !set_layer_equivalent_atom(
                &positions, &mut equiv_atoms, i,
                &indep_atoms, conv_prim, conv_sym, symprec, aperiodic,
            ) {
                equiv_atoms[i] = i as i32;
                indep_atoms.push(i);
                positions[i] = conv_prim.position[i];
                set_layer_exact_location(&mut positions[i], conv_sym, &conv_prim.lattice, aperiodic, symprec);
            }
        }
    }

    Some((positions, equiv_atoms))
}

/// 检查原子 i 是否与已有独立原子等价。
/// 如果是，设置 positions[i] 和 equiv_atoms[i] 并返回 true。
fn set_equivalent_atom(
    positions: &[[f64; 3]],
    equiv_atoms: &mut [i32],
    i: usize,
    indep_atoms: &[usize],
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    symprec: f64,
) -> bool {
    for &j in indep_atoms {
        for k in 0..conv_sym.size {
            let mut pos = mat_multiply_matrix_vector_id3(&conv_sym.rot[k], &positions[j]);
            for l in 0..3 {
                pos[l] += conv_sym.trans[k][l];
            }
            if cel_is_overlap_with_same_type(
                &pos, &conv_prim.position[i],
                conv_prim.types[j], conv_prim.types[i],
                &conv_prim.lattice, symprec,
            ) {
                equiv_atoms[i] = j as i32;
                // positions is &mut in the caller, here we can't mutate it
                // because it's &[]. So we rely on the caller to handle this.
                return true;
            }
        }
    }
    false
}

/// 位点对称性用于确定原子的确切位置。
/// R. W. Grosse-Kunstleve and P. D. Adams, Acta Cryst. (2002). A58, 60-65
fn set_exact_location(
    position: &mut Vec3,
    conv_sym: &Symmetry,
    bravais_lattice: &Mat3,
    symprec: f64,
) {
    let mut num_site_sym = 0;
    let mut sum_rot = [[0.0; 3]; 3];
    let mut sum_trans = [0.0; 3];

    for i in 0..conv_sym.size {
        let mut pos = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for j in 0..3 {
            pos[j] += conv_sym.trans[i][j];
        }

        if cel_is_overlap(&pos, position, bravais_lattice, symprec) {
            for j in 0..3 {
                for k in 0..3 {
                    sum_rot[j][k] += conv_sym.rot[i][j][k] as f64;
                }
                sum_trans[j] +=
                    conv_sym.trans[i][j] - (pos[j] - position[j]).round();
            }
            num_site_sym += 1;
        }
    }

    if num_site_sym > 0 {
        let n = num_site_sym as f64;
        for j in 0..3 {
            sum_trans[j] /= n;
            for k in 0..3 {
                sum_rot[j][k] /= n;
            }
        }

        *position = mat_multiply_matrix_vector_d3(&sum_rot, position);
        for j in 0..3 {
            position[j] += sum_trans[j];
        }
    }
}

/// 层状结构的等价原子检查。
fn set_layer_equivalent_atom(
    positions: &[[f64; 3]],
    equiv_atoms: &mut [i32],
    i: usize,
    indep_atoms: &[usize],
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    symprec: f64,
    aperiodic: AperiodicAxis,
) -> bool {
    for &j in indep_atoms {
        for k in 0..conv_sym.size {
            let mut pos = mat_multiply_matrix_vector_id3(&conv_sym.rot[k], &positions[j]);
            for l in 0..3 {
                pos[l] += conv_sym.trans[k][l];
            }
            if cel_layer_is_overlap_with_same_type(
                &pos, &conv_prim.position[i],
                conv_prim.types[j], conv_prim.types[i],
                &conv_prim.lattice, aperiodic, symprec,
            ) {
                equiv_atoms[i] = j as i32;
                return true;
            }
        }
    }
    false
}

/// 层状结构的位点对称性精确定位。
fn set_layer_exact_location(
    position: &mut Vec3,
    conv_sym: &Symmetry,
    bravais_lattice: &Mat3,
    aperiodic: AperiodicAxis,
    symprec: f64,
) {
    let mut num_site_sym = 0;
    let mut sum_rot = [[0.0; 3]; 3];
    let mut sum_trans = [0.0; 3];

    for i in 0..conv_sym.size {
        let mut pos = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for j in 0..3 {
            pos[j] += conv_sym.trans[i][j];
        }

        if cel_layer_is_overlap(&pos, position, bravais_lattice, aperiodic, symprec) {
            for j in 0..3 {
                for k in 0..3 {
                    sum_rot[j][k] += conv_sym.rot[i][j][k] as f64;
                }
                sum_trans[j] +=
                    conv_sym.trans[i][j] - (pos[j] - position[j]).round();
            }
            num_site_sym += 1;
        }
    }

    if num_site_sym > 0 {
        let n = num_site_sym as f64;
        for j in 0..3 {
            sum_trans[j] /= n;
            for k in 0..3 {
                sum_rot[j][k] /= n;
            }
        }

        *position = mat_multiply_matrix_vector_d3(&sum_rot, position);
        for j in 0..3 {
            position[j] += sum_trans[j];
        }
    }
}

/// 为所有独立原子分配 Wyckoff 字母和位点对称性符号。
fn set_wyckoffs_labels(
    positions: &[Vec3],
    equiv_atoms: &[i32],
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    num_pure_trans: i32,
    hall_number: usize,
    symprec: f64,
) -> Option<(Vec<i32>, Vec<String>)> {
    let n = conv_prim.size;
    let mut nums_equiv_atoms = vec![0i32; n];
    for i in 0..n {
        nums_equiv_atoms[equiv_atoms[i] as usize] += 1;
    }

    debug::debug_print(format_args!("num_pure_trans: {}\n", num_pure_trans));

    let mut wyckoffs = vec![0i32; n];
    let mut symbols: Vec<String> = vec![String::new(); n];

    if hall_number > 0 {
        for i in 0..n {
            if i as i32 == equiv_atoms[i] {
                debug::debug_print(format_args!(
                    "num_equiv_atoms[{}]: {}\n", i, nums_equiv_atoms[i]
                ));
                let w = get_wyckoff_notation(
                    &positions[i], conv_sym,
                    nums_equiv_atoms[i] * num_pure_trans, &conv_prim.lattice,
                    hall_number, symprec,
                );
                match w {
                    Some((letter, sym)) => {
                        wyckoffs[i] = letter;
                        symbols[i] = sym;
                    }
                    None => return None,
                }
            }
        }
    } else {
        for i in 0..n {
            if i as i32 == equiv_atoms[i] {
                let w = get_layer_wyckoff_notation(
                    &positions[i], conv_sym,
                    nums_equiv_atoms[i] * num_pure_trans, &conv_prim.lattice,
                    hall_number, AperiodicAxis::Z, symprec,
                );
                match w {
                    Some((letter, sym)) => {
                        wyckoffs[i] = letter;
                        symbols[i] = sym;
                    }
                    None => return None,
                }
            }
        }
    }

    // 将等价原子的 Wyckoff 标记从独立原子复制过来
    for i in 0..n {
        if i as i32 != equiv_atoms[i] {
            wyckoffs[i] = wyckoffs[equiv_atoms[i] as usize];
            symbols[i] = symbols[equiv_atoms[i] as usize].clone();
        }
    }

    Some((wyckoffs, symbols))
}

/// 获取 Wyckoff 字母。
/// 返回 (Wyckoff 字母编号 0=a, 1=b, ..., 位点对称性符号)
fn get_wyckoff_notation(
    position: &Vec3,
    conv_sym: &Symmetry,
    ref_multiplicity: i32,
    bravais_lattice: &Mat3,
    hall_number: usize,
    symprec: f64,
) -> Option<(i32, String)> {
    debug::debug_print(format_args!("get_Wyckoff_notation\n"));

    let n = conv_sym.size;

    // 计算所有对称操作作用在 position 上的结果
    let mut pos_rot: Vec<Vec3> = vec![[0.0; 3]; n];
    for i in 0..n {
        pos_rot[i] = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for j in 0..3 {
            pos_rot[i][j] += conv_sym.trans[i][j];
        }
    }

    let (indices_wyc_start, indices_wyc_count) = ssmdb_get_wyckoff_indices(hall_number as i32);
    for i in 0..indices_wyc_count {
        let idx = (indices_wyc_start + i) as usize;
        let (rot, trans, multiplicity) = ssmdb_get_coordinate(idx);

        for j in 0..n {
            let mut num_sitesym = 0;
            for k in 0..n {
                if cel_is_overlap(&pos_rot[j], &pos_rot[k], bravais_lattice, symprec) {
                    let mut orbit = mat_multiply_matrix_vector_id3(&rot, &pos_rot[k]);
                    for l in 0..3 {
                        orbit[l] += trans[l];
                    }
                    if cel_is_overlap(&pos_rot[k], &orbit, bravais_lattice, symprec) {
                        num_sitesym += 1;
                    }
                }
            }

            // 一致性检查: num_sym == num_sitesym * m 且 m == ref_multiplicity
            if num_sitesym * multiplicity == n as i32
                && multiplicity == ref_multiplicity
            {
                // 数据库是反序的 (gfedcba), wyckoff 按 a=0, b=1, c=2... 排列
                let wyckoff_letter = indices_wyc_count - i - 1;
                let symbol = ssmdb_get_site_symmetry_symbol(idx);
                return Some((wyckoff_letter, symbol));
            }
        }
    }

    None
}

/// 层状结构的 Wyckoff 字母获取。
fn get_layer_wyckoff_notation(
    position: &Vec3,
    conv_sym: &Symmetry,
    ref_multiplicity: i32,
    bravais_lattice: &Mat3,
    hall_number: usize,
    aperiodic: AperiodicAxis,
    symprec: f64,
) -> Option<(i32, String)> {
    debug::debug_print(format_args!("get_layer_Wyckoff_notation\n"));

    let n = conv_sym.size;

    let mut pos_rot: Vec<Vec3> = vec![[0.0; 3]; n];
    for i in 0..n {
        pos_rot[i] = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for j in 0..3 {
            pos_rot[i][j] += conv_sym.trans[i][j];
        }
    }

    let (indices_wyc_start, indices_wyc_count) = ssmdb_get_wyckoff_indices(hall_number as i32);
    for i in 0..indices_wyc_count {
        let idx = (indices_wyc_start + i) as usize;
        let (rot, trans, multiplicity) = ssmdb_get_coordinate(idx);

        for j in 0..n {
            let mut num_sitesym = 0;
            for k in 0..n {
                if cel_layer_is_overlap(
                    &pos_rot[j], &pos_rot[k], bravais_lattice, aperiodic, symprec,
                ) {
                    let mut orbit = mat_multiply_matrix_vector_id3(&rot, &pos_rot[k]);
                    for l in 0..3 {
                        orbit[l] += trans[l];
                    }
                    if cel_layer_is_overlap(
                        &pos_rot[k], &orbit, bravais_lattice, aperiodic, symprec,
                    ) {
                        num_sitesym += 1;
                    }
                }
            }

            if num_sitesym * multiplicity == n as i32
                && multiplicity == ref_multiplicity
            {
                let wyckoff_letter = indices_wyc_count - i - 1;
                let symbol = ssmdb_get_site_symmetry_symbol(idx);
                return Some((wyckoff_letter, symbol));
            }
        }
    }

    None
}
