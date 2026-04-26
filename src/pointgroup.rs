// pointgroup.rs

use crate::cell::AperiodicAxis;
use crate::debug;
use crate::mathfunc::{
    Mat3I, mat_add_matrix_i3, mat_check_identity_matrix_i3,
    mat_get_determinant_i3, mat_get_trace_i3, mat_multiply_matrix_i3,
    mat_multiply_matrix_vector_i3, mat_norm_squared_i3,
};
use crate::symmetry::PointSymmetry;

const NUM_ROT_AXES: usize = 73;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Holohedry {
    None,
    Tricli,
    Monocli,
    Ortho,
    Trigo,
    Tetra,
    Hexa,
    Cubic,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Laue {
    None,
    Laue1,
    Laue2M,
    LaueMMM,
    Laue4M,
    Laue4MMM,
    Laue3,
    Laue3M,
    Laue6M,
    Laue6MMM,
    LaueM3,
    LaueM3M,
}

#[derive(Clone, Debug)]
pub struct Pointgroup {
    pub number: usize,
    pub symbol: String,
    pub schoenflies: String,
    pub holohedry: Holohedry,
    pub laue: Laue,
}

struct PointgroupType {
    table: [i32; 10],
    symbol: &'static str,
    schoenflies: &'static str,
    holohedry: Holohedry,
    laue: Laue,
}

// Point group data table (0-32)
static POINTGROUP_DATA: [PointgroupType; 33] = [
    PointgroupType {
        table: [0; 10],
        symbol: "     ",
        schoenflies: "   ",
        holohedry: Holohedry::None,
        laue: Laue::None,
    }, // 0
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        symbol: "1    ",
        schoenflies: "C1 ",
        holohedry: Holohedry::Tricli,
        laue: Laue::Laue1,
    }, // 1
    PointgroupType {
        table: [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
        symbol: "-1   ",
        schoenflies: "Ci ",
        holohedry: Holohedry::Tricli,
        laue: Laue::Laue1,
    }, // 2
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
        symbol: "2    ",
        schoenflies: "C2 ",
        holohedry: Holohedry::Monocli,
        laue: Laue::Laue2M,
    }, // 3
    PointgroupType {
        table: [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
        symbol: "m    ",
        schoenflies: "Cs ",
        holohedry: Holohedry::Monocli,
        laue: Laue::Laue2M,
    }, // 4
    PointgroupType {
        table: [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        symbol: "2/m  ",
        schoenflies: "C2h",
        holohedry: Holohedry::Monocli,
        laue: Laue::Laue2M,
    }, // 5
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 3, 0, 0, 0],
        symbol: "222  ",
        schoenflies: "D2 ",
        holohedry: Holohedry::Ortho,
        laue: Laue::LaueMMM,
    }, // 6
    PointgroupType {
        table: [0, 0, 0, 2, 0, 1, 1, 0, 0, 0],
        symbol: "mm2  ",
        schoenflies: "C2v",
        holohedry: Holohedry::Ortho,
        laue: Laue::LaueMMM,
    }, // 7
    PointgroupType {
        table: [0, 0, 0, 3, 1, 1, 3, 0, 0, 0],
        symbol: "mmm  ",
        schoenflies: "D2h",
        holohedry: Holohedry::Ortho,
        laue: Laue::LaueMMM,
    }, // 8
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 1, 0, 2, 0],
        symbol: "4    ",
        schoenflies: "C4 ",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4M,
    }, // 9
    PointgroupType {
        table: [0, 2, 0, 0, 0, 1, 1, 0, 0, 0],
        symbol: "-4   ",
        schoenflies: "S4 ",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4M,
    }, // 10
    PointgroupType {
        table: [0, 2, 0, 1, 1, 1, 1, 0, 2, 0],
        symbol: "4/m  ",
        schoenflies: "C4h",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4M,
    }, // 11
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 5, 0, 2, 0],
        symbol: "422  ",
        schoenflies: "D4 ",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4MMM,
    }, // 12
    PointgroupType {
        table: [0, 0, 0, 4, 0, 1, 1, 0, 2, 0],
        symbol: "4mm  ",
        schoenflies: "C4v",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4MMM,
    }, // 13
    PointgroupType {
        table: [0, 2, 0, 2, 0, 1, 3, 0, 0, 0],
        symbol: "-42m ",
        schoenflies: "D2d",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4MMM,
    }, // 14
    PointgroupType {
        table: [0, 2, 0, 5, 1, 1, 5, 0, 2, 0],
        symbol: "4/mmm",
        schoenflies: "D4h",
        holohedry: Holohedry::Tetra,
        laue: Laue::Laue4MMM,
    }, // 15
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 0, 2, 0, 0],
        symbol: "3    ",
        schoenflies: "C3 ",
        holohedry: Holohedry::Trigo,
        laue: Laue::Laue3,
    }, // 16
    PointgroupType {
        table: [0, 0, 2, 0, 1, 1, 0, 2, 0, 0],
        symbol: "-3   ",
        schoenflies: "C3i",
        holohedry: Holohedry::Trigo,
        laue: Laue::Laue3,
    }, // 17
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 3, 2, 0, 0],
        symbol: "32   ",
        schoenflies: "D3 ",
        holohedry: Holohedry::Trigo,
        laue: Laue::Laue3M,
    }, // 18
    PointgroupType {
        table: [0, 0, 0, 3, 0, 1, 0, 2, 0, 0],
        symbol: "3m   ",
        schoenflies: "C3v",
        holohedry: Holohedry::Trigo,
        laue: Laue::Laue3M,
    }, // 19
    PointgroupType {
        table: [0, 0, 2, 3, 1, 1, 3, 2, 0, 0],
        symbol: "-3m  ",
        schoenflies: "D3d",
        holohedry: Holohedry::Trigo,
        laue: Laue::Laue3M,
    }, // 20
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 1, 2, 0, 2],
        symbol: "6    ",
        schoenflies: "C6 ",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6M,
    }, // 21
    PointgroupType {
        table: [2, 0, 0, 1, 0, 1, 0, 2, 0, 0],
        symbol: "-6   ",
        schoenflies: "C3h",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6M,
    }, // 22
    PointgroupType {
        table: [2, 0, 2, 1, 1, 1, 1, 2, 0, 2],
        symbol: "6/m  ",
        schoenflies: "C6h",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6M,
    }, // 23
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 7, 2, 0, 2],
        symbol: "622  ",
        schoenflies: "D6 ",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6MMM,
    }, // 24
    PointgroupType {
        table: [0, 0, 0, 6, 0, 1, 1, 2, 0, 2],
        symbol: "6mm  ",
        schoenflies: "C6v",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6MMM,
    }, // 25
    PointgroupType {
        table: [2, 0, 0, 4, 0, 1, 3, 2, 0, 0],
        symbol: "-6m2 ",
        schoenflies: "D3h",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6MMM,
    }, // 26
    PointgroupType {
        table: [2, 0, 2, 7, 1, 1, 7, 2, 0, 2],
        symbol: "6/mmm",
        schoenflies: "D6h",
        holohedry: Holohedry::Hexa,
        laue: Laue::Laue6MMM,
    }, // 27
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 3, 8, 0, 0],
        symbol: "23   ",
        schoenflies: "T  ",
        holohedry: Holohedry::Cubic,
        laue: Laue::LaueM3,
    }, // 28
    PointgroupType {
        table: [0, 0, 8, 3, 1, 1, 3, 8, 0, 0],
        symbol: "m-3  ",
        schoenflies: "Th ",
        holohedry: Holohedry::Cubic,
        laue: Laue::LaueM3,
    }, // 29
    PointgroupType {
        table: [0, 0, 0, 0, 0, 1, 9, 8, 6, 0],
        symbol: "432  ",
        schoenflies: "O  ",
        holohedry: Holohedry::Cubic,
        laue: Laue::LaueM3M,
    }, // 30
    PointgroupType {
        table: [0, 6, 0, 6, 0, 1, 3, 8, 0, 0],
        symbol: "-43m ",
        schoenflies: "Td ",
        holohedry: Holohedry::Cubic,
        laue: Laue::LaueM3M,
    }, // 31
    PointgroupType {
        table: [0, 6, 8, 9, 1, 1, 9, 8, 6, 0],
        symbol: "m-3m ",
        schoenflies: "Oh ",
        holohedry: Holohedry::Cubic,
        laue: Laue::LaueM3M,
    }, // 32
];

static IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
static INVERSION: [[i32; 3]; 3] = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];

pub(crate) static ROT_AXES: [[i32; 3]; NUM_ROT_AXES] = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [0, 1, 1],
    [1, 0, 1],
    [1, 1, 0],
    [0, 1, -1],
    [-1, 0, 1],
    [1, -1, 0],
    [1, 1, 1], /* 10 */
    [-1, 1, 1],
    [1, -1, 1],
    [1, 1, -1],
    [0, 1, 2],
    [2, 0, 1],
    [1, 2, 0],
    [0, 2, 1],
    [1, 0, 2],
    [2, 1, 0],
    [0, -1, 2], /* 20 */
    [2, 0, -1],
    [-1, 2, 0],
    [0, -2, 1],
    [1, 0, -2],
    [-2, 1, 0],
    [2, 1, 1],
    [1, 2, 1],
    [1, 1, 2],
    [2, -1, -1],
    [-1, 2, -1], /* 30 */
    [-1, -1, 2],
    [2, 1, -1],
    [-1, 2, 1],
    [1, -1, 2],
    [2, -1, 1], /* 35 */
    [1, 2, -1],
    [-1, 1, 2],
    [3, 1, 2],
    [2, 3, 1],
    [1, 2, 3], /* 40 */
    [3, 2, 1],
    [1, 3, 2],
    [2, 1, 3],
    [3, -1, 2],
    [2, 3, -1], /* 45 */
    [-1, 2, 3],
    [3, -2, 1],
    [1, 3, -2],
    [-2, 1, 3],
    [3, -1, -2], /* 50 */
    [-2, 3, -1],
    [-1, -2, 3],
    [3, -2, -1],
    [-1, 3, -2],
    [-2, -1, 3], /* 55 */
    [3, 1, -2],
    [-2, 3, 1],
    [1, -2, 3],
    [3, 2, -1],
    [-1, 3, 2], /* 60 */
    [2, -1, 3],
    [1, 1, 3],
    [-1, 1, 3],
    [1, -1, 3],
    [-1, -1, 3], /* 65 */
    [1, 3, 1],
    [-1, 3, 1],
    [1, 3, -1],
    [-1, 3, -1],
    [3, 1, 1], /* 70 */
    [3, 1, -1],
    [3, -1, 1],
    [3, -1, -1],
];

/// 从旋转操作确定点群并计算变换矩阵。
///
/// # Arguments
/// * `rotations` - 旋转操作矩阵（整数矩阵）。
/// * `aperiodic_axis` - 层状结构的非周期轴（三维周期传 `None`）。
///
/// # Returns
/// `(transform_mat, pointgroup)` — `transform_mat` 将旋转操作变换到标准方向，
/// `pointgroup` 包含点群编号、符号等信息。若无法识别点群则 `pointgroup.number == 0`。
pub fn ptg_get_transformation_matrix(
    rotations: &[Mat3I],
    aperiodic_axis: Option<AperiodicAxis>,
) -> (Mat3I, Pointgroup) {
    debug::debug_print(format_args!("ptg_get_transformation_matrix:\n"));

    let mut transform_mat = [[0; 3]; 3];

    let pg_num = get_pointgroup_number_by_rotations(rotations);

    if pg_num > 0 && (aperiodic_axis.is_none() || pg_num < 28) {
        let pointgroup = ptg_get_pointgroup(pg_num as usize);
        let pointsym = ptg_get_pointsymmetry(rotations);
        let mut axes = [0; 3];
        get_axes(&mut axes, pointgroup.laue, &pointsym, aperiodic_axis);
        set_transformation_matrix(&mut transform_mat, &axes);
        (transform_mat, pointgroup)
    } else {
        debug::info_print(format_args!("spglib: No point group was found\n"));
        let mut pg = ptg_get_pointgroup(0);
        pg.number = 0;
        (transform_mat, pg)
    }
}

pub fn ptg_get_pointgroup(pointgroup_number: usize) -> Pointgroup {
    let idx = pointgroup_number as usize;
    let pg_type = &POINTGROUP_DATA[idx];

    Pointgroup {
        number: pointgroup_number,
        symbol: pg_type.symbol.trim().to_string(),
        schoenflies: pg_type.schoenflies.trim().to_string(),
        holohedry: pg_type.holohedry,
        laue: pg_type.laue,
    }
}

pub fn ptg_get_pointsymmetry(rotations: &[Mat3I]) -> PointSymmetry {
    let mut pointsym = PointSymmetry::new(0);
    // PointSymmetry in Rust uses Vec, so we can push
    // But the C struct has fixed size array or pointer.
    // In symmetry.rs, PointSymmetry has `rot: Vec<Mat3I>`.
    // We need to re-implement the logic to find unique rotations.

    let mut unique_rots = Vec::new();
    for rot in rotations {
        let mut found = false;
        for existing in &unique_rots {
            if mat_check_identity_matrix_i3(rot, existing) {
                found = true;
                break;
            }
        }
        if !found {
            unique_rots.push(*rot);
        }
    }

    pointsym.size = unique_rots.len();
    pointsym.rot = unique_rots;
    pointsym
}

// --- Internal Functions ---

fn get_pointgroup_number_by_rotations(rotations: &[Mat3I]) -> i32 {
    let pointsym = ptg_get_pointsymmetry(rotations);
    get_pointgroup_number(&pointsym)
}

fn get_pointgroup_number(pointsym: &PointSymmetry) -> i32 {
    debug::debug_print(format_args!("get_pointgroup_number:"));

    let mut table = [0; 10];
    if !get_pointgroup_class_table(&mut table, pointsym) {
        debug::debug_print(format_args!(" 0\n"));
        return 0;
    }

    for i in 1..33 {
        let pg_type = &POINTGROUP_DATA[i];
        let mut counter = 0;
        for j in 0..10 {
            if pg_type.table[j] == table[j] {
                counter += 1;
            }
        }
        if counter == 10 {
            debug::debug_print(format_args!(" {}\n", i));
            return i as i32;
        }
    }

    debug::debug_print(format_args!(" 0\n"));
    0
}

fn get_pointgroup_class_table(table: &mut [i32; 10], pointsym: &PointSymmetry) -> bool {
    for i in 0..10 {
        table[i] = 0;
    }
    for i in 0..pointsym.size {
        let rot_type = get_rotation_type(&pointsym.rot[i]);
        if rot_type == -1 {
            debug::warning_print(format_args!("spglib: No point group symbol found\n"));
            return false;
        } else {
            table[rot_type as usize] += 1;
        }
    }
    true
}

fn get_rotation_type(rot: &Mat3I) -> i32 {
    let det = mat_get_determinant_i3(rot);
    let trace = mat_get_trace_i3(rot);

    if det == -1 {
        match trace {
            -2 => 0, // -6
            -1 => 1, // -4
            0 => 2,  // -3
            1 => 3,  // -2
            -3 => 4, // -1
            _ => -1,
        }
    } else {
        match trace {
            3 => 5,  // 1
            -1 => 6, // 2
            0 => 7,  // 3
            1 => 8,  // 4
            2 => 9,  // 6
            _ => -1,
        }
    }
}

fn get_axes(axes: &mut [usize; 3], laue: Laue, pointsym: &PointSymmetry, aperiodic_axis: Option<AperiodicAxis>) {
    match laue {
        Laue::Laue1 => {
            axes[0] = 0;
            axes[1] = 1;
            axes[2] = 2;
        }
        Laue::Laue2M => {
            if aperiodic_axis.is_none() {
                laue2m(axes, pointsym);
            } else {
                layer_laue2m(axes, pointsym, aperiodic_axis.map(|a| a.axis_index() as i32).unwrap_or(-1));
            }
        }
        Laue::LaueMMM => {
            lauennn(axes, pointsym, 2, aperiodic_axis.map(|a| a.axis_index() as i32).unwrap_or(-1));
        }
        Laue::Laue4M | Laue::Laue4MMM => {
            laue_one_axis(axes, pointsym, 4);
        }
        Laue::Laue3 | Laue::Laue3M | Laue::Laue6M | Laue::Laue6MMM => {
            laue_one_axis(axes, pointsym, 3);
        }
        Laue::LaueM3 => {
            lauennn(axes, pointsym, 2, -1);
        }
        Laue::LaueM3M => {
            lauennn(axes, pointsym, 4, -1);
        }
        _ => {}
    }
}

fn laue2m(axes: &mut [usize; 3], pointsym: &PointSymmetry) -> bool {
    let mut prop_rot = [[0; 3]; 3];

    // Find unique axis b (first axis in axes array for now, will be swapped later if needed logic in C was commented out)
    // C code: axes[1] = get_rotation_axis
    for i in 0..pointsym.size {
        get_proper_rotation(&mut prop_rot, &pointsym.rot[i]);
        if mat_get_trace_i3(&prop_rot) == -1 {
            if let Some(axis) = get_rotation_axis(&prop_rot) {
                axes[1] = axis;
                break;
            }
        }
    }

    let mut ortho_axes = [0; NUM_ROT_AXES];
    let num_ortho_axis = get_orthogonal_axis(&mut ortho_axes, &prop_rot, 2);
    if num_ortho_axis == 0 {
        return false;
    }

    let mut min_norm = 8;
    let mut is_found = false;
    for i in 0..num_ortho_axis {
        let norm = mat_norm_squared_i3(&ROT_AXES[ortho_axes[i]]);
        if norm < min_norm {
            min_norm = norm;
            axes[0] = ortho_axes[i];
            is_found = true;
        }
    }
    if !is_found {
        return false;
    }

    min_norm = 8;
    is_found = false;
    for i in 0..num_ortho_axis {
        let norm = mat_norm_squared_i3(&ROT_AXES[ortho_axes[i]]);
        if norm < min_norm && ortho_axes[i] != axes[0] {
            min_norm = norm;
            axes[2] = ortho_axes[i];
            is_found = true;
        }
    }
    if !is_found {
        return false;
    }

    true
}

fn layer_laue2m(axes: &mut [usize; 3], pointsym: &PointSymmetry, aperiodic_axis: i32) -> bool {
    let mut prop_rot = [[0; 3]; 3];

    for i in 0..pointsym.size {
        get_proper_rotation(&mut prop_rot, &pointsym.rot[i]);
        if mat_get_trace_i3(&prop_rot) == -1 {
            if let Some(axis) = get_rotation_axis(&prop_rot) {
                axes[0] = axis;
                break;
            }
        }
    }

    let mut ortho_axes = [0; NUM_ROT_AXES];
    let num_ortho_axis = get_orthogonal_axis(&mut ortho_axes, &prop_rot, 2);
    if num_ortho_axis == 0 {
        return false;
    }

    let axis0_vec = ROT_AXES[axes[0]];
    let ap_idx = aperiodic_axis as usize;

    if axis0_vec[ap_idx] == 1 || axis0_vec[ap_idx] == -1 {
        // Monoclinic/Oblique
        let mut min_norm = 8;
        let mut is_found = false;
        for i in 0..num_ortho_axis {
            let norm = mat_norm_squared_i3(&ROT_AXES[ortho_axes[i]]);
            if norm < min_norm {
                min_norm = norm;
                axes[1] = ortho_axes[i];
                is_found = true;
            }
        }
        if !is_found {
            return false;
        }

        min_norm = 8;
        is_found = false;
        for i in 0..num_ortho_axis {
            let norm = mat_norm_squared_i3(&ROT_AXES[ortho_axes[i]]);
            if norm < min_norm && ortho_axes[i] != axes[1] {
                min_norm = norm;
                axes[2] = ortho_axes[i];
                is_found = true;
            }
        }
        if !is_found {
            return false;
        }
    } else if axis0_vec[ap_idx] == 0 {
        // Monoclinic/Rectangular
        let mut min_norm = 8;
        let mut is_found = false;
        for i in 0..num_ortho_axis {
            if ROT_AXES[ortho_axes[i]][ap_idx] == 0 {
                let norm = mat_norm_squared_i3(&ROT_AXES[ortho_axes[i]]);
                if norm < min_norm {
                    min_norm = norm;
                    axes[1] = ortho_axes[i];
                    is_found = true;
                }
            }
        }
        if !is_found {
            return false;
        }

        min_norm = 8;
        is_found = false;
        for i in 0..num_ortho_axis {
            let v = ROT_AXES[ortho_axes[i]][ap_idx];
            if v == 1 || v == -1 {
                let norm = mat_norm_squared_i3(&ROT_AXES[ortho_axes[i]]);
                if norm < min_norm {
                    min_norm = norm;
                    axes[2] = ortho_axes[i];
                    is_found = true;
                }
            }
        }
        if !is_found {
            return false;
        }
    } else {
        return false;
    }

    true
}

fn laue_one_axis(axes: &mut [usize; 3], pointsym: &PointSymmetry, rot_order: i32) -> bool {
    debug::debug_print(format_args!("laue_one_axis with rot_order {}\n", rot_order));

    let mut prop_rot = [[0; 3]; 3];

    for i in 0..pointsym.size {
        get_proper_rotation(&mut prop_rot, &pointsym.rot[i]);
        let trace = mat_get_trace_i3(&prop_rot);
        if (rot_order == 4 && trace == 1) || (rot_order == 3 && trace == 0) {
            if let Some(axis) = get_rotation_axis(&prop_rot) {
                axes[2] = axis;
                break;
            }
        }
    }

    let mut ortho_axes = [0; NUM_ROT_AXES];
    let num_ortho_axis = get_orthogonal_axis(&mut ortho_axes, &prop_rot, rot_order);
    if num_ortho_axis == 0 {
        return false;
    }

    let mut tmp_axes = [0; 3];
    tmp_axes[1] = usize::MAX; // Placeholder
    tmp_axes[2] = axes[2];

    for i in 0..num_ortho_axis {
        tmp_axes[0] = ortho_axes[i];
        let axis_vec = mat_multiply_matrix_vector_i3(&prop_rot, &ROT_AXES[tmp_axes[0]]);

        let mut is_found = 0;
        for j in 0..num_ortho_axis {
            is_found = is_exist_axis(&axis_vec, ortho_axes[j]);
            if is_found == 1 {
                tmp_axes[1] = ortho_axes[j];
                break;
            }
            if is_found == -1 {
                tmp_axes[1] = ortho_axes[j] + NUM_ROT_AXES;
                break;
            }
        }

        if is_found == 0 {
            continue;
        }

        let mut t_mat = [[0; 3]; 3];
        set_transformation_matrix(&mut t_mat, &tmp_axes);
        let det = mat_get_determinant_i3(&t_mat).abs();
        if det < 4 {
            axes[0] = tmp_axes[0];
            axes[1] = tmp_axes[1];


            return true;
        }
    }

    debug::warning_print(format_args!("spglib: Secondary axis is not found.\n"));
    false
}

fn lauennn(
    axes: &mut [usize; 3],
    pointsym: &PointSymmetry,
    rot_order: i32,
    aperiodic_axis: i32,
) -> bool {
    for i in 0..3 {
        axes[i] = usize::MAX;
    } // Using MAX as -1

    let mut count = 0;
    let mut prop_rot = [[0; 3]; 3];

    for i in 0..pointsym.size {
        get_proper_rotation(&mut prop_rot, &pointsym.rot[i]);
        let trace = mat_get_trace_i3(&prop_rot);

        if (trace == -1 && rot_order == 2) || (trace == 1 && rot_order == 4) {
            if let Some(axis) = get_rotation_axis(&prop_rot) {
                if axis != axes[0] && axis != axes[1] && axis != axes[2] {
                    axes[count] = axis;
                    count += 1;
                }
            }
        }
    }

    if aperiodic_axis == -1 {
        sort_axes(axes);
    } else {
        layer_check_and_sort_axes(axes, aperiodic_axis);
    }
    true
}

fn get_rotation_axis(proper_rot: &Mat3I) -> Option<usize> {
    if mat_check_identity_matrix_i3(proper_rot, &IDENTITY) {
        return None;
    }

    for i in 0..NUM_ROT_AXES {
        let vec = mat_multiply_matrix_vector_i3(proper_rot, &ROT_AXES[i]);
        if vec == ROT_AXES[i] {
            return Some(i);
        }
    }
    debug::debug_print(format_args!("rotation axis could not be found.\n"));
    None
}

fn get_orthogonal_axis(ortho_axes: &mut [usize], proper_rot: &Mat3I, rot_order: i32) -> usize {
    let mut num_ortho_axis = 0;
    let mut sum_rot = IDENTITY;
    let mut rot = IDENTITY;

    for _ in 0..(rot_order - 1) {
        rot = mat_multiply_matrix_i3(&rot, proper_rot);
        let temp = mat_add_matrix_i3(&sum_rot, &rot);
        sum_rot = temp;
    }

    for i in 0..NUM_ROT_AXES {
        let vec = mat_multiply_matrix_vector_i3(&sum_rot, &ROT_AXES[i]);
        if vec == [0, 0, 0] {
            ortho_axes[num_ortho_axis] = i;
            num_ortho_axis += 1;
        }
    }

    num_ortho_axis
}

fn get_proper_rotation(prop_rot: &mut Mat3I, rot: &Mat3I) {
    if mat_get_determinant_i3(rot) == -1 {
        *prop_rot = mat_multiply_matrix_i3(&INVERSION, rot);
    } else {
        *prop_rot = *rot;
    }
}

fn set_transformation_matrix(tmat: &mut Mat3I, axes: &[usize; 3]) {
    let mut s = [0; 3];
    for i in 0..3 {
        if axes[i] < NUM_ROT_AXES {
            s[i] = 1;
        } else {
            s[i] = -1;
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            tmat[i][j] = s[j] * ROT_AXES[axes[j] % NUM_ROT_AXES][i];
        }
    }
}

fn is_exist_axis(axis_vec: &[i32; 3], axis_index: usize) -> i32 {
    let target = ROT_AXES[axis_index];
    if *axis_vec == target {
        return 1;
    }
    if axis_vec[0] == -target[0] && axis_vec[1] == -target[1] && axis_vec[2] == -target[2] {
        return -1;
    }
    0
}

fn sort_axes(axes: &mut [usize; 3]) {
    if axes[1] > axes[2] {
        axes.swap(1, 2);
    }
    if axes[0] > axes[1] {
        axes.swap(0, 1);
    }
    if axes[1] > axes[2] {
        axes.swap(1, 2);
    }

    let mut t_mat = [[0; 3]; 3];
    set_transformation_matrix(&mut t_mat, axes);
    if mat_get_determinant_i3(&t_mat) < 0 {
        axes.swap(1, 2);
    }
}

fn layer_check_and_sort_axes(axes: &mut [usize; 3], aperiodic_axis: i32) {
    let ap_idx = aperiodic_axis as usize;
    let mut lattice_rank = 0;
    let mut arank = 0;
    let mut axis_i = 0;

    for i in 0..3 {
        let v = ROT_AXES[axes[i]][ap_idx];
        if v == 0 {
            lattice_rank += 1;
        } else if v == 1 || v == -1 {
            axis_i = i;
            arank += 1;
        }
    }

    if lattice_rank == 2 && arank == 1 {
        axes.swap(axis_i, 2);

        let mut t_mat = [[0; 3]; 3];
        set_transformation_matrix(&mut t_mat, axes);
        if mat_get_determinant_i3(&t_mat) < 0 {
            axes.swap(0, 1);
        }
    } else {
        debug::warning_print(format_args!("spglib: Invalid axes were found\n"));
        for i in 0..3 {
            axes[i] = 0;
        }
    }
}
