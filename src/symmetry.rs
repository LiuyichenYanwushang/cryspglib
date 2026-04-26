//! 对称操作检测。
//!
//! 在给定精度下寻找晶胞的所有对称操作（旋转 + 平移）。
//! 核心函数 [`sym_get_operation`] 返回包含旋转矩阵（i32）和平移向量（f64）的 [`Symmetry`] 结构体。

use crate::cell::{AperiodicAxis, Cell, cel_is_overlap_with_same_type, cel_layer_is_overlap_with_same_type};
use crate::debug;
use crate::delaunay::{del_delaunay_reduce, del_layer_delaunay_reduce};
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_cast_matrix_3d_to_3i, mat_cast_matrix_3i_to_3d,
    mat_check_identity_matrix_i3,
    mat_dabs, mat_dmod1, mat_get_determinant_d3, mat_get_determinant_i3, mat_get_metric,
    mat_get_similar_matrix_d3, mat_inverse_matrix_d3, mat_is_int_matrix, mat_multiply_matrix_d3,
    mat_multiply_matrix_di3, mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_id3,
};
use crate::overlap::OverlapChecker;
use std::f64::consts::PI;

// 常量定义
const NUM_ATOMS_CRITERION_FOR_OPENMP: usize = 1000; // 虽然 Rust 使用 Rayon，但保留此常量以对应 C 逻辑
const ANGLE_REDUCE_RATE: f64 = 0.95;
const SIN_DTHETA2_CUTOFF: f64 = 1e-12;
const NUM_ATTEMPT: i32 = 100;

// 相对轴向量，用于生成所有可能的晶格基矢量变换矩阵 (3x3x3 - 1 = 26 个方向)
static RELATIVE_AXES: [[i32; 3]; 26] = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1],
    [0, 1, 1], [1, 0, 1], [1, 1, 0], [0, -1, -1], [-1, 0, -1], [-1, -1, 0],
    [0, 1, -1], [-1, 0, 1], [1, -1, 0], [0, -1, 1], [1, 0, -1], [-1, 1, 0],
    [1, 1, 1], [-1, -1, -1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1],
    [-1, 1, -1], [-1, -1, 1],
];

static IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// 对称操作集合。
///
/// 包含 `size` 个对称操作 `(W, w)`，其中 W 是 3x3 整数旋转矩阵，
/// w 是分数平移向量。对原子坐标 x，操作后为 `W*x + w`。
#[derive(Clone, Debug)]
pub struct Symmetry {
    /// 对称操作数量
    pub size: usize,
    /// 旋转矩阵列表（每个为 3x3 i32 矩阵）
    pub rot: Vec<Mat3I>,
    /// 平移向量列表（分数坐标）
    pub trans: Vec<Vec3>,
}

impl Symmetry {
    pub fn new(size: usize) -> Self {
        Symmetry {
            size,
            rot: vec![[[0; 3]; 3]; size],
            trans: vec![[0.0; 3]; size],
        }
    }
}

/// 点群对称性结构体
#[derive(Clone, Debug)]
pub struct PointSymmetry {
    pub size: usize,
    pub rot: Vec<Mat3I>,
}

impl PointSymmetry {
    pub fn new(size: usize) -> Self {
        PointSymmetry {
            size,
            rot: vec![[[0; 3]; 3]; size],
        }
    }
}

/// 磁性对称操作结构体
#[derive(Clone, Debug)]
pub struct MagneticSymmetry {
    pub size: usize,
    pub rot: Vec<Mat3I>,
    pub trans: Vec<Vec3>,
    pub timerev: Vec<bool>,
}

impl MagneticSymmetry {
    pub fn new(size: usize) -> Self {
        MagneticSymmetry {
            size,
            rot: vec![[[0; 3]; 3]; size],
            trans: vec![[0.0; 3]; size],
            timerev: vec![false; size],
        }
    }
}

// --- Public API ---

/// 获取晶胞的对称操作
pub fn sym_get_operation(primitive: &Cell, symprec: f64, angle_tolerance: f64) -> Option<Symmetry> {
    debug::debug_print(format_args!("sym_get_operations:\n"));
    get_operations(primitive, symprec, angle_tolerance)
}

/// 约化对称操作
pub fn sym_reduce_operation(
    primitive: &Cell,
    symmetry: &Symmetry,
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Symmetry> {
    reduce_operation(primitive, symmetry, symprec, angle_tolerance, false)
}

/// 获取纯平移操作
pub fn sym_get_pure_translation(cell: &Cell, symprec: f64) -> Option<Vec<Vec3>> {
    debug::debug_print(format_args!(
        "sym_get_pure_translation (tolerance = {}):\n",
        symprec
    ));

    let pure_trans = if cell.aperiodic_axis.is_none() {
        get_translation(&IDENTITY, cell, symprec, true)
    } else {
        get_layer_translation(&IDENTITY, cell, symprec, true)
    };

    if let Some(ref pt) = pure_trans {
        let multi = pt.len();
        // 检查原子数是否是平移重数的整数倍
        if (cell.size / multi) * multi == cell.size {
            debug::debug_print(format_args!(
                "spglib: sym_get_pure_translation: pure_trans->size = {}\n",
                multi
            ));
        } else {
            debug::warning_print(format_args!(
                "spglib: Finding pure translation failed.\n        cell->size {}, multi {}\n",
                cell.size, multi
            ));
        }
    } else {
        debug::debug_print(format_args!("spglib: get_translation failed.\n"));
    }

    pure_trans
}

/// 约化纯平移操作
pub fn sym_reduce_pure_translation(
    cell: &Cell,
    pure_trans: &[Vec3],
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Vec<Vec3>> {
    let multi = pure_trans.len();
    let mut symmetry = Symmetry::new(multi);

    for i in 0..multi {
        symmetry.rot[i] = IDENTITY;
        symmetry.trans[i] = pure_trans[i];
    }

    let symmetry_reduced = reduce_operation(cell, &symmetry, symprec, angle_tolerance, true)?;

    Some(symmetry_reduced.trans)
}

// --- Internal Functions ---

fn get_operations(primitive: &Cell, symprec: f64, angle_symprec: f64) -> Option<Symmetry> {
    debug::debug_print(format_args!("get_operations:\n"));

    let lattice_sym = get_lattice_symmetry(primitive, symprec, angle_symprec);
    if lattice_sym.size == 0 {
        return None;
    }

    get_space_group_operations(&lattice_sym, primitive, symprec)
}

fn reduce_operation(
    primitive: &Cell,
    symmetry: &Symmetry,
    symprec: f64,
    angle_symprec: f64,
    is_pure_trans: bool,
) -> Option<Symmetry> {
    debug::debug_print(format_args!("reduce_operation:\n"));

    let point_symmetry = if is_pure_trans {
        let mut ps = PointSymmetry::new(1);
        ps.rot[0] = IDENTITY;
        ps
    } else {
        let ps = get_lattice_symmetry(primitive, symprec, angle_symprec);
        if ps.size == 0 {
            return None;
        }
        ps
    };

    let mut rot_list = Vec::new();
    let mut trans_list = Vec::new();

    for i in 0..point_symmetry.size {
        for j in 0..symmetry.size {
            if mat_check_identity_matrix_i3(&point_symmetry.rot[i], &symmetry.rot[j]) {
                if is_overlap_all_atoms(
                    &symmetry.trans[j],
                    &symmetry.rot[j],
                    primitive,
                    symprec,
                    false,
                ) == 1
                {
                    rot_list.push(symmetry.rot[j]);
                    trans_list.push(symmetry.trans[j]);
                }
            }
        }
    }

    let mut sym_reduced = Symmetry::new(rot_list.len());
    sym_reduced.rot = rot_list;
    sym_reduced.trans = trans_list;

    Some(sym_reduced)
}

fn get_translation(rot: &Mat3I, cell: &Cell, symprec: f64, is_identity: bool) -> Option<Vec<Vec3>> {
    debug::debug_print(format_args!("get_translation (tolerance = {}):\n", symprec));

    let min_atom_index = get_index_with_least_atoms(cell);
    if min_atom_index == -1 {
        debug::debug_print(format_args!("spglib: get_index_with_least_atoms failed.\n"));
        return None;
    }
    let min_atom_index = min_atom_index as usize;

    let origin = mat_multiply_matrix_vector_id3(rot, &cell.position[min_atom_index]);

    let Some((is_found, num_trans)) = search_translation_part(
        cell,
        rot,
        min_atom_index,
        &origin,
        symprec,
        is_identity,
    ) else {
        return None;
    };

    if num_trans == 0 {
        return None;
    }

    let mut trans = Vec::with_capacity(num_trans);
    for i in 0..cell.size {
        if is_found[i] {
            let mut t = [0.0; 3];
            for j in 0..3 {
                t[j] = cell.position[i][j] - origin[j];
                t[j] = mat_dmod1(t[j]);
            }
            trans.push(t);
        }
    }

    Some(trans)
}

fn search_translation_part(
    cell: &Cell,
    rot: &Mat3I,
    min_atom_index: usize,
    origin: &Vec3,
    symprec: f64,
    is_identity: bool,
) -> Option<(Vec<bool>, usize)> {
    let mut checker = OverlapChecker::new(cell)?;
    let mut atoms_found = vec![false; cell.size];
    let mut num_trans = 0;

    for i in 0..cell.size {
        if atoms_found[i] {
            continue;
        }
        if cell.types[i] != cell.types[min_atom_index] {
            continue;
        }

        let mut trans = [0.0; 3];
        for j in 0..3 {
            trans[j] = cell.position[i][j] - origin[j];
        }

        let is_overlap = checker.check_total_overlap(&trans, rot, symprec, is_identity);
        if is_overlap == -1 {
            return None;
        } else if is_overlap == 1 {
            atoms_found[i] = true;
            num_trans += 1;
            if is_identity {
                num_trans += search_pure_translations(&mut atoms_found, cell, &trans, symprec);
            }
        }
    }
    Some((atoms_found, num_trans))
}

fn search_pure_translations(
    atoms_found: &mut [bool],
    cell: &Cell,
    trans: &Vec3,
    symprec: f64,
) -> usize {
    let mut num_trans = 0;
    let copy_atoms_found = atoms_found.to_vec();

    for initial_atom in 0..cell.size {
        if !copy_atoms_found[initial_atom] {
            continue;
        }

        let mut i_atom = initial_atom;
        for _ in 0..cell.size {
            let mut vec = [0.0; 3];
            for j in 0..3 {
                vec[j] = cell.position[i_atom][j] + trans[j];
            }

            for j in 0..cell.size {
                if cel_is_overlap_with_same_type(
                    &vec,
                    &cell.position[j],
                    cell.types[i_atom],
                    cell.types[j],
                    &cell.lattice,
                    symprec,
                ) {
                    if !atoms_found[j] {
                        atoms_found[j] = true;
                        num_trans += 1;
                    }
                    i_atom = j;
                    break;
                }
            }
            if i_atom == initial_atom {
                break;
            }
        }
    }
    num_trans
}

fn is_overlap_all_atoms(
    trans: &Vec3,
    rot: &Mat3I,
    cell: &Cell,
    symprec: f64,
    is_identity: bool,
) -> i32 {
    let mut checker = match OverlapChecker::new(cell) {
        Some(c) => c,
        None => return -1,
    };

    if cell.aperiodic_axis.is_none() {
        checker.check_total_overlap(trans, rot, symprec, is_identity)
    } else {
        checker.check_layer_total_overlap(trans, rot, symprec, is_identity)
    }
}

fn get_index_with_least_atoms(cell: &Cell) -> i32 {
    let mut mapping = vec![0; cell.size];
    for i in 0..cell.size {
        for j in 0..cell.size {
            if cell.types[i] == cell.types[j] {
                mapping[j] += 1;
                break;
            }
        }
    }

    let mut min = mapping[0];
    let mut min_index = 0;
    for i in 0..cell.size {
        if min > mapping[i] && mapping[i] > 0 {
            min = mapping[i];
            min_index = i;
        }
    }
    min_index as i32
}

fn get_layer_translation(
    rot: &Mat3I,
    cell: &Cell,
    symprec: f64,
    is_identity: bool,
) -> Option<Vec<Vec3>> {
    debug::debug_print(format_args!("get_translation (tolerance = {}):\n", symprec));

    let min_atom_index = get_index_with_least_atoms(cell);
    if min_atom_index == -1 {
        debug::debug_print(format_args!("spglib: get_index_with_least_atoms failed.\n"));
        return None;
    }
    let min_atom_index = min_atom_index as usize;

    let origin = mat_multiply_matrix_vector_id3(rot, &cell.position[min_atom_index]);

    let Some((is_found, num_trans)) = search_layer_translation_part(
        cell,
        rot,
        min_atom_index,
        &origin,
        symprec,
        is_identity,
    ) else {
        return None;
    };

    if num_trans == 0 {
        return None;
    }

    let mut trans = Vec::with_capacity(num_trans);
    for i in 0..cell.size {
        if is_found[i] {
            let mut t = [0.0; 3];
            for j in 0..3 {
                t[j] = cell.position[i][j] - origin[j];
                if cell.aperiodic_axis.map_or(true, |ap| j != ap.axis_index()) {
                    t[j] = mat_dmod1(t[j]);
                }
            }
            trans.push(t);
        }
    }
    Some(trans)
}

fn search_layer_translation_part(
    cell: &Cell,
    rot: &Mat3I,
    min_atom_index: usize,
    origin: &Vec3,
    symprec: f64,
    is_identity: bool,
) -> Option<(Vec<bool>, usize)> {
    let mut checker = OverlapChecker::new(cell)?;
    let mut atoms_found = vec![false; cell.size];
    let mut num_trans = 0;

    for i in 0..cell.size {
        if atoms_found[i] {
            continue;
        }
        if cell.types[i] != cell.types[min_atom_index] {
            continue;
        }

        let mut trans = [0.0; 3];
        for j in 0..3 {
            trans[j] = cell.position[i][j] - origin[j];
        }

        let is_overlap = checker.check_layer_total_overlap(&trans, rot, symprec, is_identity);
        if is_overlap == -1 {
            return None;
        } else if is_overlap == 1 {
            atoms_found[i] = true;
            num_trans += 1;
            if is_identity {
                let aperiodic = cell.aperiodic_axis.unwrap();
                num_trans += search_layer_pure_translations(
                    &mut atoms_found,
                    cell,
                    &trans,
                    aperiodic,
                    symprec,
                );
            }
        }
    }
    Some((atoms_found, num_trans))
}

fn search_layer_pure_translations(
    atoms_found: &mut [bool],
    cell: &Cell,
    trans: &Vec3,
    aperiodic: AperiodicAxis,
    symprec: f64,
) -> usize {
    let mut num_trans = 0;
    let copy_atoms_found = atoms_found.to_vec();

    for initial_atom in 0..cell.size {
        if !copy_atoms_found[initial_atom] {
            continue;
        }
        let mut i_atom = initial_atom;
        for _ in 0..cell.size {
            let mut vec = [0.0; 3];
            for j in 0..3 {
                vec[j] = cell.position[i_atom][j] + trans[j];
            }
            for j in 0..cell.size {
                if cel_layer_is_overlap_with_same_type(
                    &vec,
                    &cell.position[j],
                    cell.types[i_atom],
                    cell.types[j],
                    &cell.lattice,
                    aperiodic,
                    symprec,
                ) {
                    if !atoms_found[j] {
                        atoms_found[j] = true;
                        num_trans += 1;
                    }
                    i_atom = j;
                    break;
                }
            }
            if i_atom == initial_atom {
                break;
            }
        }
    }
    num_trans
}

fn get_space_group_operations(
    lattice_sym: &PointSymmetry,
    primitive: &Cell,
    symprec: f64,
) -> Option<Symmetry> {
    debug::debug_print(format_args!(
        "get_space_group_operations (tolerance = {}):\n",
        symprec
    ));

    let mut trans_vecs: Vec<Option<Vec<Vec3>>> = Vec::with_capacity(lattice_sym.size);
    let mut total_num_sym = 0;

    for i in 0..lattice_sym.size {
        let t = if primitive.aperiodic_axis.is_none() {
            get_translation(&lattice_sym.rot[i], primitive, symprec, false)
        } else {
            get_layer_translation(&lattice_sym.rot[i], primitive, symprec, false)
        };

        if let Some(ref v) = t {
            debug::debug_print(format_args!(
                "  match translation {}/{}; tolerance = {}\n",
                i + 1,
                lattice_sym.size,
                symprec
            ));
            total_num_sym += v.len();
        }
        trans_vecs.push(t);
    }

    let mut symmetry = Symmetry::new(total_num_sym);
    let mut num_sym = 0;

    for i in 0..lattice_sym.size {
        if let Some(ref vecs) = trans_vecs[i] {
            for v in vecs {
                symmetry.trans[num_sym] = *v;
                symmetry.rot[num_sym] = lattice_sym.rot[i];
                num_sym += 1;
            }
        }
    }

    Some(symmetry)
}

pub fn get_lattice_symmetry(cell: &Cell, symprec: f64, angle_symprec: f64) -> PointSymmetry {
    debug::debug_print(format_args!("get_lattice_symmetry:\n"));

    let mut lattice_sym = PointSymmetry::new(0); // Initially empty
    let aperiodic_axis = cell.aperiodic_axis;

    let Some(min_lattice) = (if aperiodic_axis.is_none() {
        del_delaunay_reduce(&cell.lattice, symprec)
    } else {
        del_layer_delaunay_reduce(&cell.lattice, aperiodic_axis, symprec)
    }) else {
        debug::debug_print(format_args!("get_lattice_symmetry failed.\n"));
        return lattice_sym;
    };

    let metric_orig = mat_get_metric(&min_lattice);
    let mut angle_tol = angle_symprec;

    // 使用循环标签 'attempt_loop 来模拟 C 代码中的 goto next_attempt 逻辑
    'attempt_loop: for _attempt in 0..NUM_ATTEMPT {
        let mut rot_list = Vec::new();
        let mut axes = [[0; 3]; 3];

        for i in 0..26 {
            for j in 0..26 {
                for k in 0..26 {
                    set_axes(&mut axes, i, j, k);

                    // Layer group checks
                    match aperiodic_axis {
                        Some(AperiodicAxis::Z) => {
                            if axes[0][2] != 0
                                || axes[1][2] != 0
                                || axes[2][0] != 0
                                || axes[2][1] != 0
                            {
                                continue;
                            }
                        }
                        Some(AperiodicAxis::X) => {
                            if axes[0][1] != 0
                                || axes[0][2] != 0
                                || axes[1][0] != 0
                                || axes[2][0] != 0
                            {
                                continue;
                            }
                        }
                        Some(AperiodicAxis::Y) => {
                            if axes[0][1] != 0
                                || axes[1][0] != 0
                                || axes[1][2] != 0
                                || axes[2][1] != 0
                            {
                                continue;
                            }
                        }
                        _ => {}
                    }

                    let det = mat_get_determinant_i3(&axes);
                    if det != 1 && det != -1 {
                        continue;
                    }

                    let lattice_trans = mat_multiply_matrix_di3(&min_lattice, &axes);
                    let metric = mat_get_metric(&lattice_trans);

                    if is_identity_metric(&metric, &metric_orig, symprec, angle_tol) {
                        if (aperiodic_axis.is_none() && rot_list.len() >= 48)
                            || (aperiodic_axis.is_some() && rot_list.len() >= 24)
                        {
                            debug::debug_print(format_args!(
                                "spglib: Too many lattice symmetries were found.\n"
                            ));
                            if angle_tol > 0.0 {
                                angle_tol *= ANGLE_REDUCE_RATE;
                                debug::debug_print(format_args!(
                                    "        Reducing angle tolerance to {}\n",
                                    angle_tol
                                ));
                                // 重新开始外层循环
                                continue 'attempt_loop;
                            }
                            // angle_tol <= 0, continue collecting symmetries
                        }
                        rot_list.push(axes);
                    }
                }
            }
        }

        if !rot_list.is_empty() {
            if (aperiodic_axis.is_none() && rot_list.len() <= 48)
                || (aperiodic_axis.is_some() && rot_list.len() <= 24)
                || angle_tol < 0.0
            {
                lattice_sym.size = rot_list.len();
                lattice_sym.rot = rot_list;
                return transform_pointsymmetry(&lattice_sym, &cell.lattice, &min_lattice);
            }
        }
    }

    debug::debug_print(format_args!("get_lattice_symmetry failed.\n"));
    lattice_sym
}

fn is_identity_metric(
    metric_rotated: &Mat3,
    metric_orig: &Mat3,
    symprec: f64,
    angle_symprec: f64,
) -> bool {
    let elem_sets = [[0, 1], [0, 2], [1, 2]];
    let mut length_orig = [0.0; 3];
    let mut length_rot = [0.0; 3];

    for i in 0..3 {
        length_orig[i] = metric_orig[i][i].sqrt();
        length_rot[i] = metric_rotated[i][i].sqrt();
        if mat_dabs(length_orig[i] - length_rot[i]) > symprec {
            return false;
        }
    }

    for i in 0..3 {
        let j = elem_sets[i][0];
        let k = elem_sets[i][1];
        if angle_symprec > 0.0 {
            if mat_dabs(get_angle(metric_orig, j, k) - get_angle(metric_rotated, j, k))
                > angle_symprec
            {
                return false;
            }
        } else {
            let cos1 = metric_orig[j][k] / length_orig[j] / length_orig[k];
            let cos2 = metric_rotated[j][k] / length_rot[j] / length_rot[k];
            let x = cos1 * cos2 + (1.0 - cos1 * cos1).sqrt() * (1.0 - cos2 * cos2).sqrt();
            let sin_dtheta2 = 1.0 - x * x;
            let length_ave2 =
                ((length_orig[j] + length_rot[j]) * (length_orig[k] + length_rot[k])) / 4.0;
            if sin_dtheta2 > SIN_DTHETA2_CUTOFF {
                if sin_dtheta2 * length_ave2 > symprec * symprec {
                    return false;
                }
            }
        }
    }
    true
}

fn get_angle(metric: &Mat3, i: usize, j: usize) -> f64 {
    let length_i = metric[i][i].sqrt();
    let length_j = metric[j][j].sqrt();
    (metric[i][j] / length_i / length_j).acos() / PI * 180.0
}

fn transform_pointsymmetry(
    lat_sym_orig: &PointSymmetry,
    new_lattice: &Mat3,
    original_lattice: &Mat3,
) -> PointSymmetry {
    let mut lat_sym_new = PointSymmetry::new(0);
    let mut rot_list = Vec::new();

    let inv_mat = mat_inverse_matrix_d3(original_lattice, 0.0).ok().unwrap_or([[0.0; 3]; 3]);
    let trans_mat = mat_multiply_matrix_d3(&inv_mat, new_lattice);

    for i in 0..lat_sym_orig.size {
        let mut drot = mat_cast_matrix_3i_to_3d(&lat_sym_orig.rot[i]);
        // 尝试获取相似矩阵，如果失败则跳过
        if let Ok(sim) = mat_get_similar_matrix_d3(&drot, &trans_mat, 0.0) {
            drot = sim;
            if mat_is_int_matrix(&drot, mat_dabs(mat_get_determinant_d3(&trans_mat)) / 10.0) {
                let rot_i = mat_cast_matrix_3d_to_3i(&drot);
                if mat_get_determinant_i3(&rot_i).abs() != 1 {
                    debug::warning_print(format_args!(
                        "spglib: A point symmetry operation is not unimodular.\n"
                    ));
                    return lat_sym_new; // Return empty on error
                }
                rot_list.push(rot_i);
            }
        }
    }

    if lat_sym_orig.size != rot_list.len() {
        debug::warning_print(format_args!(
            "spglib: Some of point symmetry operations were dropped.\n"
        ));
    }

    lat_sym_new.size = rot_list.len();
    lat_sym_new.rot = rot_list;
    lat_sym_new
}

fn set_axes(axes: &mut Mat3I, a1: usize, a2: usize, a3: usize) {
    for i in 0..3 {
        axes[i][0] = RELATIVE_AXES[a1][i];
        axes[i][1] = RELATIVE_AXES[a2][i];
        axes[i][2] = RELATIVE_AXES[a3][i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::TensorRank;

    #[test]
    fn test_sym_get_pure_translation_identity() {
        // 构造一个简单的立方晶胞
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0]];
        let types = [1];
        let mut cell = Cell::new(1, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);

        let trans = sym_get_pure_translation(&cell, 1e-5);
        assert!(trans.is_some());
        let t = trans.unwrap();
        assert_eq!(t.len(), 1);
        // 纯平移应包含 (0,0,0)
        assert!(t[0][0].abs() < 1e-5);
    }

    #[test]
    fn test_sym_get_pure_translation_supercell() {
        // 构造一个 2x1x1 超胞
        let lattice = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let types = [1, 1];
        let mut cell = Cell::new(2, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);

        let trans = sym_get_pure_translation(&cell, 1e-5);
        assert!(trans.is_some());
        let t = trans.unwrap();
        assert_eq!(t.len(), 2);
        // 应包含 (0,0,0) 和 (0.5,0,0)
        let has_zero = t.iter().any(|v| v[0].abs() < 1e-5);
        let has_half = t.iter().any(|v| (v[0] - 0.5).abs() < 1e-5);
        assert!(has_zero);
        assert!(has_half);
    }
}
