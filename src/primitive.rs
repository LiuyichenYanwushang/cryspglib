//! 原胞查找。
//!
//! 利用检测到的纯平移对称性，将常规晶胞约化为原胞 (primitive cell)。
//! 如果存在多种可能的原胞，选择 Delaunay 约化后晶格最小的那个。
//!
//! 核心函数 [`prm_get_primitive`] 返回包含原胞晶格、位置、类型和
//! 纯平移操作的 [`Primitive`] 结构体，以及其对称操作。

use crate::cell::{Cell, cel_trim_cell};
use crate::debug;
use crate::delaunay::{del_delaunay_reduce, del_layer_delaunay_reduce};
use crate::mathfunc::{
    Mat3, Vec3, mat_cast_matrix_3d_to_3i, mat_cast_matrix_3i_to_3d, mat_check_identity_matrix_i3,
    mat_dabs, mat_dmod1, mat_get_determinant_d3,
    mat_get_determinant_i3, mat_inverse_matrix_d3, mat_multiply_matrix_d3, mat_multiply_matrix_di3,
    mat_multiply_matrix_vector_d3, mat_nint,
};
use crate::symmetry::{Symmetry, sym_get_pure_translation, sym_reduce_pure_translation};

const REDUCE_RATE: f64 = 0.95;
const NUM_ATTEMPT: i32 = 100;

pub struct Primitive {
    pub cell: Option<Cell>,
    pub mapping_table: Vec<i32>,
    pub size: usize,
    pub tolerance: f64,
    pub angle_tolerance: f64,
    pub orig_lattice: Mat3,
}

impl Primitive {
    pub fn new(size: usize) -> Self {
        Primitive {
            cell: None,
            mapping_table: vec![-1; size],
            size,
            tolerance: 0.0,
            angle_tolerance: -1.0,
            orig_lattice: [[0.0; 3]; 3],
        }
    }
}

pub fn prm_get_primitive(cell: &Cell, symprec: f64, angle_tolerance: f64) -> Option<Primitive> {
    get_primitive(cell, symprec, angle_tolerance)
}

pub fn prm_get_primitive_with_pure_trans(
    cell: &Cell,
    pure_trans: &[Vec3],
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Primitive> {
    let mut primitive = Primitive::new(cell.size);

    if pure_trans.len() == 1 {
        primitive.cell = get_cell_with_smallest_lattice(cell, symprec);
        if primitive.cell.is_none() {
            return None;
        }
        for i in 0..cell.size {
            primitive.mapping_table[i] = i as i32;
        }
    } else {
        let mut mapping_table_usize = vec![0; cell.size];
        primitive.cell = get_primitive_cell(
            &mut mapping_table_usize,
            cell,
            pure_trans,
            symprec,
            angle_tolerance,
        );
        if primitive.cell.is_none() {
            return None;
        }
        for i in 0..cell.size {
            primitive.mapping_table[i] = mapping_table_usize[i] as i32;
        }
    }

    primitive.tolerance = symprec;
    primitive.angle_tolerance = angle_tolerance;
    primitive.orig_lattice = cell.lattice;

    Some(primitive)
}

/// 从含中心化平移的对称操作提取原胞对称性和变换矩阵。
///
/// 给定可能包含中心化平移的对称操作集合，返回原胞表示：
/// `(t_mat, prim_symmetry)`，其中 `t_mat` 将原胞坐标变换到原始坐标。
///
/// # Arguments
/// * `symmetry` - 原始对称操作（可能含中心化平移，如 R/F/I/C 中心化）。
/// * `symprec` - 对称性检查容差。
///
/// # Returns
/// `Some((t_mat, prim_symmetry))` — `t_mat` 是从原胞到原始晶胞的变换矩阵，
/// `prim_symmetry` 包含原胞中的对称操作（旋转矩阵和分数平移）。
pub fn prm_get_primitive_symmetry(
    symmetry: &Symmetry,
    symprec: f64,
) -> Option<(Mat3, Symmetry)> {
    let pure_trans = collect_pure_translations(symmetry)?;
    let primsym_size = symmetry.size / pure_trans.len();

    let t_mat_inv = get_primitive_in_translation_space(&pure_trans, symmetry.size, symprec)?;
    let t_mat = mat_inverse_matrix_d3(&t_mat_inv, symprec).ok()?;

    let mut prim_symmetry = collect_primitive_symmetry(symmetry, primsym_size)?;

    for i in 0..prim_symmetry.size {
        let mut tmp_mat = mat_multiply_matrix_di3(&t_mat, &prim_symmetry.rot[i]);
        tmp_mat = mat_multiply_matrix_d3(&tmp_mat, &t_mat_inv);
        prim_symmetry.rot[i] = mat_cast_matrix_3d_to_3i(&tmp_mat);
        prim_symmetry.trans[i] = mat_multiply_matrix_vector_d3(&t_mat, &prim_symmetry.trans[i]);
    }

    for i in 0..prim_symmetry.size {
        debug::debug_print(format_args!("--- {} ---\n", i + 1));
        for j in 0..3 {
            debug::debug_print(format_args!(
                "{} {} {}\n",
                prim_symmetry.rot[i][j][0], prim_symmetry.rot[i][j][1], prim_symmetry.rot[i][j][2]
            ));
        }
        debug::debug_print(format_args!(
            "{} {} {}\n",
            prim_symmetry.trans[i][0], prim_symmetry.trans[i][1], prim_symmetry.trans[i][2]
        ));
    }

    Some((t_mat, prim_symmetry))
}

pub fn prm_get_primitive_lattice_vectors(
    cell: &Cell,
    pure_trans: &[Vec3],
    symprec: f64,
    angle_tolerance: f64,
) -> Option<(Mat3, usize)> {
    get_primitive_lattice_vectors(cell, pure_trans, symprec, angle_tolerance)
}

// --- Internal Functions ---

fn get_primitive(cell: &Cell, symprec: f64, angle_tolerance: f64) -> Option<Primitive> {
    debug::debug_print(format_args!("get_primitive (tolerance = {}):\n", symprec));

    let mut tolerance = symprec;

    for attempt in 0..NUM_ATTEMPT {
        debug::debug_print(format_args!("get_primitive (attempt = {}):\n", attempt));
        if let Some(pure_trans) = sym_get_pure_translation(cell, tolerance) {
            if let Some(primitive) = prm_get_primitive_with_pure_trans(
                cell,
                &pure_trans,
                tolerance,
                angle_tolerance,
            ) {
                return Some(primitive);
            }
        }

        tolerance *= REDUCE_RATE;
        debug::debug_print(format_args!("spglib: Reduce tolerance to {} ", tolerance));
    }

    None
}

fn get_cell_with_smallest_lattice(cell: &Cell, symprec: f64) -> Option<Cell> {
    debug::debug_print(format_args!("get_cell_with_smallest_lattice:\n"));

    let aperiodic_axis = cell.aperiodic_axis;

    let min_lat = if aperiodic_axis.is_none() {
        del_delaunay_reduce(&cell.lattice, symprec)?
    } else {
        del_layer_delaunay_reduce(&cell.lattice, aperiodic_axis, symprec)?
    };

    let inv_lat = mat_inverse_matrix_d3(&min_lat, 0.0).ok()?;
    let trans_mat = mat_multiply_matrix_d3(&inv_lat, &cell.lattice);

    let mut smallest_cell = Cell::new(cell.size, cell.tensor_rank);
    smallest_cell.lattice = min_lat;

    for i in 0..cell.size {
        smallest_cell.types[i] = cell.types[i];
        smallest_cell.position[i] = mat_multiply_matrix_vector_d3(&trans_mat, &cell.position[i]);
        for j in 0..3 {
            if aperiodic_axis.map_or(false, |ap| j == ap.axis_index()) {
                smallest_cell.aperiodic_axis = aperiodic_axis;
            } else {
                smallest_cell.position[i][j] = mat_dmod1(smallest_cell.position[i][j]);
            }
        }
    }

    Some(smallest_cell)
}

fn get_primitive_cell(
    mapping_table: &mut [usize],
    cell: &Cell,
    pure_trans: &[Vec3],
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Cell> {
    debug::debug_print(format_args!("get_primitive_cell:\n"));

    let Some((prim_lattice, _)) = get_primitive_lattice_vectors(
        cell,
        pure_trans,
        symprec,
        angle_tolerance,
    ) else {
        debug::debug_print(format_args!("spglib: Primitive cell could not be found\n"));
        return None;
    };

    cel_trim_cell(mapping_table, &prim_lattice, cell, symprec)
}

fn get_primitive_lattice_vectors(
    cell: &Cell,
    pure_trans: &[Vec3],
    symprec: f64,
    angle_tolerance: f64,
) -> Option<(Mat3, usize)> {
    let mut tolerance = symprec;
    let mut pure_trans_reduced = pure_trans.to_vec();

    for attempt in 0..NUM_ATTEMPT {
        let multi = pure_trans_reduced.len();
        let vectors = get_translation_candidates(&pure_trans_reduced);

        if let Some(found) = find_primitive_lattice_vectors(&vectors, cell, tolerance) {
            let lattice_copy = found;
            if let Some(reduced) = if cell.aperiodic_axis.is_none() {
                println!("reduced={:?}",del_delaunay_reduce(&lattice_copy, symprec));
                del_delaunay_reduce(&lattice_copy, symprec)
            } else {
                del_layer_delaunay_reduce(&lattice_copy, cell.aperiodic_axis, symprec)
            } {
                return Some((reduced, multi));
            } else {
                return None;
            }
        } else {
            if let Some(reduced) =
                sym_reduce_pure_translation(cell, &pure_trans_reduced, tolerance, angle_tolerance)
            {
                pure_trans_reduced = reduced;
                debug::debug_print(format_args!(
                    "spglib: Tolerance is reduced to {} ({}), num_pure_trans = {}\n",
                    tolerance,
                    attempt,
                    pure_trans_reduced.len()
                ));
                tolerance *= REDUCE_RATE;
            } else {
                return None;
            }
        }
    }
    None
}

fn find_primitive_lattice_vectors(
    vectors: &[Vec3],
    cell: &Cell,
    symprec: f64,
) -> Option<Mat3> {
    debug::debug_print(format_args!("find_primitive_lattice_vectors:\n"));

    let size = vectors.len();
    let initial_volume = mat_dabs(mat_get_determinant_d3(&cell.lattice));
    let aperiodic_axis = cell.aperiodic_axis;
    let mut min_vectors = [[0.0; 3]; 3];
    let mut found = false;

    if aperiodic_axis.is_none() {
        'outer: for i in 0..size {
            for j in (i + 1)..size {
                for k in (j + 1)..size {
                    let mut tmp_lattice = [[0.0; 3]; 3];
                    tmp_lattice[0] = mat_multiply_matrix_vector_d3(&cell.lattice, &vectors[i]);
                    tmp_lattice[1] = mat_multiply_matrix_vector_d3(&cell.lattice, &vectors[j]);
                    tmp_lattice[2] = mat_multiply_matrix_vector_d3(&cell.lattice, &vectors[k]);

                    let volume = mat_dabs(mat_get_determinant_d3(&tmp_lattice));
                    if volume > symprec {
                        if mat_nint(initial_volume / volume) == (size - 2) as i32 {
                            min_vectors[0] = vectors[i];
                            min_vectors[1] = vectors[j];
                            min_vectors[2] = vectors[k];
                            found = true;
                            break 'outer;
                        }
                    }
                }
            }
        }
    } else {
        let ap_idx = aperiodic_axis.unwrap().axis_index();
        let k_idx = size + ap_idx - 3;
        'outer_layer: for i in 0..size {
            for j in (i + 1)..size {
                if i != k_idx && j != k_idx {
                    let mut tmp_lattice = [[0.0; 3]; 3];
                    tmp_lattice[0] = mat_multiply_matrix_vector_d3(&cell.lattice, &vectors[i]);
                    tmp_lattice[1] = mat_multiply_matrix_vector_d3(&cell.lattice, &vectors[j]);
                    tmp_lattice[2] = mat_multiply_matrix_vector_d3(&cell.lattice, &vectors[k_idx]);

                    let volume = mat_dabs(mat_get_determinant_d3(&tmp_lattice));
                    if volume > symprec {
                        if mat_nint(initial_volume / volume) == (size - 2) as i32 {
                            min_vectors[0] = vectors[i];
                            min_vectors[1] = vectors[j];
                            if ap_idx == 2 {
                                min_vectors[2] = vectors[k_idx];
                            } else {
                                min_vectors[2] = min_vectors[ap_idx];
                                min_vectors[ap_idx] = vectors[k_idx];
                            }
                            found = true;
                            break 'outer_layer;
                        }
                    }
                }
            }
        }
    }

    if !found {
        debug::debug_print(format_args!(
            "spglib: Primitive lattice vectors could not be found\n"
        ));
        return None;
    }

    let mut relative_lattice = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            relative_lattice[j][i] = min_vectors[i][j];
        }
    }

    if let Ok(inv_mat_dbl) = mat_inverse_matrix_d3(&relative_lattice, 0.0) {
        let inv_mat_int = mat_cast_matrix_3d_to_3i(&inv_mat_dbl);
        if mat_get_determinant_i3(&inv_mat_int).abs() == (size - 2) as i32 {
            let inv_mat_dbl_clean = mat_cast_matrix_3i_to_3d(&inv_mat_int);
            if let Ok(rel) = mat_inverse_matrix_d3(&inv_mat_dbl_clean, 0.0) {
                relative_lattice = rel;
            }
        } else {
            debug::warning_print(format_args!(
                "spglib: Primitive lattice cleaning is incomplete\n"
            ));
        }
    }

    Some(mat_multiply_matrix_d3(&cell.lattice, &relative_lattice))
}

fn get_translation_candidates(pure_trans: &[Vec3]) -> Vec<Vec3> {
    let multi = pure_trans.len();
    let mut vectors = Vec::with_capacity(multi + 2);

    // Store pure translations (skipping the first one which is usually 0,0,0 or identity)
    // Note: C code loops i from 0 to multi-2, accessing pure_trans[i+1]
    for i in 0..(multi - 1) {
        vectors.push(pure_trans[i + 1]);
    }

    // Store lattice translations (1,0,0), (0,1,0), (0,0,1)
    for i in 0..3 {
        let mut v = [0.0; 3];
        v[i] = 1.0;
        vectors.push(v);
    }

    vectors
}

fn collect_pure_translations(symmetry: &Symmetry) -> Option<Vec<Vec3>> {
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let mut pure_trans = Vec::new();

    for i in 0..symmetry.size {
        if mat_check_identity_matrix_i3(&symmetry.rot[i], &identity) {
            pure_trans.push(symmetry.trans[i]);
        }
    }

    if pure_trans.is_empty() {
        None
    } else {
        Some(pure_trans)
    }
}

/// 在平移空间中找到原胞变换矩阵的逆。
///
/// 利用纯平移构建一个人工晶胞（单位晶格），搜索其原胞变换，
/// 返回原胞晶格矩阵作为 `t_mat_inv`。
///
/// # Returns
/// `Some(t_mat_inv)` — 从原胞到原始坐标的变换矩阵的逆。
fn get_primitive_in_translation_space(
    pure_trans: &[Vec3],
    symmetry_size: usize,
    symprec: f64,
) -> Option<Mat3> {
    let primsym_size = symmetry_size / pure_trans.len();
    if symmetry_size != primsym_size * pure_trans.len() {
        return None;
    }

    let mut cell = Cell::new(pure_trans.len(), crate::cell::TensorRank::NoSpin);
    for i in 0..pure_trans.len() {
        cell.types[i] = 1;
        cell.position[i] = pure_trans[i];
    }
    cell.lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let primitive = get_primitive(&cell, symprec, -1.0)?;
    let prim_cell = primitive.cell?;
    if prim_cell.size != 1 {
        return None;
    }
    Some(prim_cell.lattice)
}

fn collect_primitive_symmetry(symmetry: &Symmetry, primsym_size: usize) -> Option<Symmetry> {
    let mut prim_symmetry = Symmetry::new(primsym_size);
    let mut num_psym = 0;

    // First one
    prim_symmetry.rot[0] = symmetry.rot[0];
    prim_symmetry.trans[0] = symmetry.trans[0];
    num_psym += 1;

    for i in 1..symmetry.size {
        let mut is_found = true;
        for j in 0..num_psym {
            if mat_check_identity_matrix_i3(&prim_symmetry.rot[j], &symmetry.rot[i]) {
                is_found = false;
                break;
            }
        }
        if is_found {
            if num_psym == primsym_size {
                return None;
            }
            prim_symmetry.rot[num_psym] = symmetry.rot[i];
            prim_symmetry.trans[num_psym] = symmetry.trans[i];
            num_psym += 1;
        }
    }

    if num_psym != primsym_size {
        return None;
    }

    Some(prim_symmetry)
}
