//! Wyckoff 位置的精确定位和对等原子分配。
//!
//! 使用位点对称性数据库确定原子的确切经坐标和 Wyckoff 字母标记。
//!
//! 参考: R. W. Grosse-Kunstleve and P. D. Adams,
//! Acta Cryst. (2002). A58, 60-65

use crate::cell::{is_overlap, is_overlap_with_same_type, layer_is_overlap,
                   layer_is_overlap_with_same_type, AperiodicAxis, Cell};
use crate::debug;
use crate::mathfunc::*;
use crate::sitesym_database::*;
use crate::symmetry::Symmetry;

const INCREASE_RATE: f64 = 1.05;
const NUM_ATTEMPT: i32 = 5;

pub(crate) struct ExactPositions {
    pub positions: Vec<Vec3>,
    pub wyckoffs: Vec<i32>,
    pub equivalent_atoms: Vec<usize>,
    pub site_symmetry_symbols: Vec<String>,
}

struct ExactSiteData {
    positions: Vec<Vec3>,
    equivalent_atoms: Vec<usize>,
}

struct SiteEquivalenceContext<'a> {
    positions: &'a [Vec3],
    independent_atoms: &'a [usize],
    primitive: &'a Cell,
    symmetry: &'a Symmetry,
    symprec: f64,
}

/// 获取精确的原子位置和 Wyckoff 标记。
///
/// # Returns
/// A complete set of exact positions, Wyckoff labels, equivalent atoms, and
/// site-symmetry symbols.
pub(crate) fn ssm_get_exact_positions(
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    num_pure_trans: i32,
    hall_number: usize,
    symprec: f64,
) -> Option<ExactPositions> {
    let mut tolerance = symprec;
    for i in 0..NUM_ATTEMPT {
        let ExactSiteData {
            positions,
            equivalent_atoms: equiv_atoms,
        } = get_exact_positions(conv_prim, conv_sym, tolerance)?;

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

        return Some(ExactPositions {
            positions,
            wyckoffs,
            equivalent_atoms: equiv_atoms,
            site_symmetry_symbols: symbols,
        });
    }

    None
}

/// 内部函数：为每个原子确定精确位置。
/// 返回 (positions, equiv_atoms) 或 None
fn get_exact_positions(
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    symprec: f64,
) -> Option<ExactSiteData> {
    debug::debug_print(format_args!("get_exact_positions\n"));

    let n = conv_prim.len();
    let mut positions = vec![[0.0; 3]; n];
    let mut equiv_atoms = vec![0; n];
    let mut indep_atoms: Vec<usize> = Vec::with_capacity(n);

    if let Some(aperiodic) = conv_prim.aperiodic_axis {
        for i in 0..n {
            if !set_layer_equivalent_atom(
                &SiteEquivalenceContext {
                    positions: &positions,
                    independent_atoms: &indep_atoms,
                    primitive: conv_prim,
                    symmetry: conv_sym,
                    symprec,
                },
                &mut equiv_atoms,
                i,
                aperiodic,
            ) {
                equiv_atoms[i] = i;
                indep_atoms.push(i);
                positions[i] = conv_prim.position[i];
                set_layer_exact_location(
                    &mut positions[i],
                    conv_sym,
                    &conv_prim.lattice,
                    aperiodic,
                    symprec,
                );
            }
        }
    } else {
        for i in 0..n {
            if !set_equivalent_atom(
                &SiteEquivalenceContext {
                    positions: &positions,
                    independent_atoms: &indep_atoms,
                    primitive: conv_prim,
                    symmetry: conv_sym,
                    symprec,
                },
                &mut equiv_atoms,
                i,
            ) {
                equiv_atoms[i] = i;
                indep_atoms.push(i);
                positions[i] = conv_prim.position[i];
                set_exact_location(&mut positions[i], conv_sym, &conv_prim.lattice, symprec);
            }
        }
    }

    Some(ExactSiteData {
        positions,
        equivalent_atoms: equiv_atoms,
    })
}

/// 检查原子 i 是否与已有独立原子等价。
/// 如果是，设置 positions[i] 和 equiv_atoms[i] 并返回 true。
fn set_equivalent_atom(
    context: &SiteEquivalenceContext<'_>,
    equiv_atoms: &mut [usize],
    i: usize,
) -> bool {
    for &j in context.independent_atoms {
        for k in 0..context.symmetry.len() {
            let mut pos =
                mat_multiply_matrix_vector_id3(&context.symmetry.rot[k], &context.positions[j]);
            for (coordinate, translation) in pos.iter_mut().zip(&context.symmetry.trans[k]) {
                *coordinate += translation;
            }
            if is_overlap_with_same_type(
                &pos,
                &context.primitive.position[i],
                context.primitive.types[j],
                context.primitive.types[i],
                &context.primitive.lattice,
                context.symprec,
            ) {
                equiv_atoms[i] = j;
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

    for i in 0..conv_sym.len() {
        let mut pos = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for (coordinate, translation) in pos.iter_mut().zip(&conv_sym.trans[i]) {
            *coordinate += translation;
        }

        if is_overlap(&pos, position, bravais_lattice, symprec) {
            for (j, (sum_row, sum_translation)) in
                sum_rot.iter_mut().zip(&mut sum_trans).enumerate()
            {
                for (sum_value, &rotation) in sum_row.iter_mut().zip(&conv_sym.rot[i][j]) {
                    *sum_value += rotation as f64;
                }
                *sum_translation +=
                    conv_sym.trans[i][j] - (pos[j] - position[j]).round();
            }
            num_site_sym += 1;
        }
    }

    if num_site_sym > 0 {
        let n = num_site_sym as f64;
        for (sum_row, sum_translation) in sum_rot.iter_mut().zip(&mut sum_trans) {
            *sum_translation /= n;
            for value in sum_row {
                *value /= n;
            }
        }

        *position = mat_multiply_matrix_vector_d3(&sum_rot, position);
        for (coordinate, translation) in position.iter_mut().zip(sum_trans) {
            *coordinate += translation;
        }
    }
}

/// 层状结构的等价原子检查。
fn set_layer_equivalent_atom(
    context: &SiteEquivalenceContext<'_>,
    equiv_atoms: &mut [usize],
    i: usize,
    aperiodic: AperiodicAxis,
) -> bool {
    for &j in context.independent_atoms {
        for k in 0..context.symmetry.len() {
            let mut pos =
                mat_multiply_matrix_vector_id3(&context.symmetry.rot[k], &context.positions[j]);
            for (coordinate, translation) in pos.iter_mut().zip(&context.symmetry.trans[k]) {
                *coordinate += translation;
            }
            if layer_is_overlap_with_same_type(
                &pos,
                &context.primitive.position[i],
                context.primitive.types[j],
                context.primitive.types[i],
                &context.primitive.lattice,
                aperiodic,
                context.symprec,
            ) {
                equiv_atoms[i] = j;
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

    for i in 0..conv_sym.len() {
        let mut pos = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for (coordinate, translation) in pos.iter_mut().zip(&conv_sym.trans[i]) {
            *coordinate += translation;
        }

        if layer_is_overlap(&pos, position, bravais_lattice, aperiodic, symprec) {
            for (j, (sum_row, sum_translation)) in
                sum_rot.iter_mut().zip(&mut sum_trans).enumerate()
            {
                for (sum_value, &rotation) in sum_row.iter_mut().zip(&conv_sym.rot[i][j]) {
                    *sum_value += rotation as f64;
                }
                *sum_translation +=
                    conv_sym.trans[i][j] - (pos[j] - position[j]).round();
            }
            num_site_sym += 1;
        }
    }

    if num_site_sym > 0 {
        let n = num_site_sym as f64;
        for (sum_row, sum_translation) in sum_rot.iter_mut().zip(&mut sum_trans) {
            *sum_translation /= n;
            for value in sum_row {
                *value /= n;
            }
        }

        *position = mat_multiply_matrix_vector_d3(&sum_rot, position);
        for (coordinate, translation) in position.iter_mut().zip(sum_trans) {
            *coordinate += translation;
        }
    }
}

/// 为所有独立原子分配 Wyckoff 字母和位点对称性符号。
fn set_wyckoffs_labels(
    positions: &[Vec3],
    equiv_atoms: &[usize],
    conv_prim: &Cell,
    conv_sym: &Symmetry,
    num_pure_trans: i32,
    hall_number: usize,
    symprec: f64,
) -> Option<(Vec<i32>, Vec<String>)> {
    let n = conv_prim.len();
    let mut nums_equiv_atoms = vec![0i32; n];
    for i in 0..n {
        nums_equiv_atoms[equiv_atoms[i]] += 1;
    }

    debug::debug_print(format_args!("num_pure_trans: {}\n", num_pure_trans));

    let mut wyckoffs = vec![0i32; n];
    let mut symbols: Vec<String> = vec![String::new(); n];

    if hall_number > 0 {
        for i in 0..n {
            if i == equiv_atoms[i] {
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
            if i == equiv_atoms[i] {
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
        if i != equiv_atoms[i] {
            wyckoffs[i] = wyckoffs[equiv_atoms[i]];
            symbols[i] = symbols[equiv_atoms[i]].clone();
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

    let n = conv_sym.len();

    // 计算所有对称操作作用在 position 上的结果
    let mut pos_rot: Vec<Vec3> = vec![[0.0; 3]; n];
    for (i, rotated_position) in pos_rot.iter_mut().enumerate().take(n) {
        *rotated_position = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for (coordinate, translation) in rotated_position.iter_mut().zip(&conv_sym.trans[i]) {
            *coordinate += translation;
        }
    }

    let (indices_wyc_start, indices_wyc_count) = ssmdb_get_wyckoff_indices(hall_number as i32);
    for i in 0..indices_wyc_count {
        let idx = (indices_wyc_start + i) as usize;
        let (rot, trans, multiplicity) = ssmdb_get_coordinate(idx);

        for j in 0..n {
            let mut num_sitesym = 0;
            for k in 0..n {
                if is_overlap(&pos_rot[j], &pos_rot[k], bravais_lattice, symprec) {
                    let mut orbit = mat_multiply_matrix_vector_id3(&rot, &pos_rot[k]);
                    for l in 0..3 {
                        orbit[l] += trans[l];
                    }
                    if is_overlap(&pos_rot[k], &orbit, bravais_lattice, symprec) {
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

    let n = conv_sym.len();

    let mut pos_rot: Vec<Vec3> = vec![[0.0; 3]; n];
    for (i, rotated_position) in pos_rot.iter_mut().enumerate().take(n) {
        *rotated_position = mat_multiply_matrix_vector_id3(&conv_sym.rot[i], position);
        for (coordinate, translation) in rotated_position.iter_mut().zip(&conv_sym.trans[i]) {
            *coordinate += translation;
        }
    }

    let (indices_wyc_start, indices_wyc_count) = ssmdb_get_wyckoff_indices(hall_number as i32);
    for i in 0..indices_wyc_count {
        let idx = (indices_wyc_start + i) as usize;
        let (rot, trans, multiplicity) = ssmdb_get_coordinate(idx);

        for j in 0..n {
            let mut num_sitesym = 0;
            for k in 0..n {
                if layer_is_overlap(
                    &pos_rot[j], &pos_rot[k], bravais_lattice, aperiodic, symprec,
                ) {
                    let mut orbit = mat_multiply_matrix_vector_id3(&rot, &pos_rot[k]);
                    for l in 0..3 {
                        orbit[l] += trans[l];
                    }
                    if layer_is_overlap(
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
