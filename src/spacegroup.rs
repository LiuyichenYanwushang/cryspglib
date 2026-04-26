//! 空间群搜索与识别。
//!
//! 本模块实现空间群识别的最后一步：对 Hall 符号匹配产生的候选编号，
//! 逐一尝试原点平移和常规晶胞变换，选择最佳匹配并返回完整的 [`Spacegroup`] 信息。
//!
//! # 入口函数
//!
//! - [`spa_search_spacegroup`] — 从原胞出发的完整空间群搜索
//! - [`spa_search_spacegroup_with_symmetry`] — 给定对称操作的空间群搜索
//! - [`spa_transform_to_primitive`] / [`spa_transform_from_primitive`] — 坐标变换辅助
//!
//! # 晶格矩阵约定
//!
//! 所有 [`Mat3`] 采用 `lattice[cart][vec]` 布局（行=笛卡尔分量, 列=晶格矢量）。
//! 详见 [`mathfunc`](crate::mathfunc) 模块文档。

use crate::cell::{AperiodicAxis, Cell, TensorRank, cel_trim_cell};
use crate::debug;
use crate::delaunay::del_layer_delaunay_reduce_2D;

use crate::hall_symbol::hal_match_hall_symbol_db;
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_cast_matrix_3d_to_3i, mat_cast_matrix_3i_to_3d,
    mat_check_identity_matrix_d3, mat_check_identity_matrix_id3,
    mat_dabs, mat_get_determinant_d3, mat_get_determinant_i3, mat_get_metric,
    mat_get_similar_matrix_d3, mat_inverse_matrix_d3, mat_is_int_matrix, mat_multiply_matrix_d3,
    mat_multiply_matrix_di3, mat_multiply_matrix_id3, mat_multiply_matrix_vector_d3,
};
use crate::niggli::niggli_reduce;
use crate::pointgroup::{
    Holohedry, Laue, ptg_get_pointsymmetry, ptg_get_transformation_matrix,
};
use crate::primitive::Primitive;
use crate::spg_database::{Centering, spgdb_get_spacegroup_type};
use crate::symmetry::{Symmetry, sym_get_operation, sym_reduce_operation};

const REDUCE_RATE: f64 = 0.95;
const NUM_ATTEMPT: i32 = 100;
const INT_PREC: f64 = 0.1;
const ZERO_PREC: f64 = 1e-10;

// --- Static Data (Matrices) ---

pub static M_BCC: Mat3I = [[0, 1, 1], [1, 0, 1], [1, 1, 0]];
pub static M_FCC: Mat3I = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]];
pub static M_AC: Mat3I = [[1, 0, 0], [0, 1, 1], [0, -1, 1]];
pub static M_BC: Mat3I = [[1, 0, 1], [0, 1, 0], [-1, 0, 1]];
pub static M_CC: Mat3I = [[1, -1, 0], [1, 1, 0], [0, 0, 1]];
pub static M_RC: Mat3I = [[1, 0, 1], [-1, 1, 1], [0, -1, 1]];

static M_BCC_INV: Mat3 = [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
static M_FCC_INV: Mat3 = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];
static M_AC_INV: Mat3 = [[1.0, 0.0, 0.0], [0.0, 0.5, -0.5], [0.0, 0.5, 0.5]];
static M_BC_INV: Mat3 = [[0.5, 0.0, -0.5], [0.0, 1.0, 0.0], [0.5, 0.0, 0.5]];
static M_CC_INV: Mat3 = [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]];
static M_RC_INV: Mat3 = [
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., -2. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];

static IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
static MONOCLI_I2C: Mat3 = [[1.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]];
static MONOCLI_A2C: Mat3 = [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]];
static RHOMBO_OBVERSE: Mat3 = [
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., -2. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];
static RHOMBO_REVERSE: Mat3 = [
    [1. / 3., -2. / 3., 1. / 3.],
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];
static A2C: Mat3 = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
static B2C: Mat3 = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]];

static A_MAT: Mat3 = [[1.0, 0.0, 0.0], [0.0, 0.5, -0.5], [0.0, 0.5, 0.5]];
static C_MAT: Mat3 = [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]];
static R_MAT: Mat3 = [
    [2. / 3., -1. / 3., -1. / 3.],
    [1. / 3., 1. / 3., -2. / 3.],
    [1. / 3., 1. / 3., 1. / 3.],
];
static I_MAT: Mat3 = [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
static F_MAT: Mat3 = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];

static CHANGE_OF_BASIS_MONOCLI: [Mat3; 48] = [
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    [[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]],
    [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
    [[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]],
    [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
    [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
    [[-1., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
    [[1., 0., -1.], [0., -1., 0.], [0., 0., -1.]],
    [[0., 1., -1.], [1., 0., 0.], [0., 0., -1.]],
    [[-1., -1., 0.], [0., 0., 1.], [-1., 0., 0.]],
    [[1., -1., 0.], [0., 0., 1.], [0., -1., 0.]],
    [[0., 1., 1.], [1., 0., 0.], [0., 1., 0.]],
    [[0., 0., -1.], [0., 1., 0.], [1., 0., -1.]],
    [[-1., 0., 0.], [0., -1., 0.], [-1., 0., 1.]],
    [[0., -1., 0.], [1., 0., 0.], [0., -1., 1.]],
    [[0., 1., 0.], [0., 0., 1.], [1., 1., 0.]],
    [[-1., 0., 0.], [0., 0., 1.], [-1., 1., 0.]],
    [[0., 0., -1.], [1., 0., 0.], [0., -1., -1.]],
    [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
    [[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]],
    [[0., 0., 1.], [-1., 0., 0.], [0., -1., 0.]],
    [[-1., 0., 0.], [0., 0., -1.], [0., -1., 0.]],
    [[0., 1., 0.], [0., 0., -1.], [-1., 0., 0.]],
    [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
    [[-1., 0., -1.], [0., -1., 0.], [-1., 0., 0.]],
    [[1., 0., 1.], [0., 1., 0.], [0., 0., 1.]],
    [[0., -1., -1.], [-1., 0., 0.], [0., 0., -1.]],
    [[1., -1., 0.], [0., 0., -1.], [1., 0., 0.]],
    [[-1., -1., 0.], [0., 0., -1.], [0., -1., 0.]],
    [[0., -1., 1.], [-1., 0., 0.], [0., -1., 0.]],
    [[0., 0., 1.], [0., -1., 0.], [1., 0., 1.]],
    [[-1., 0., 0.], [0., 1., 0.], [-1., 0., -1.]],
    [[0., 1., 0.], [-1., 0., 0.], [0., 1., 1.]],
    [[0., 1., 0.], [0., 0., -1.], [-1., 1., 0.]],
    [[1., 0., 0.], [0., 0., -1.], [1., 1., 0.]],
    [[0., 0., -1.], [-1., 0., 0.], [0., 1., -1.]],
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    [[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]],
    [[0., -1., 0.], [1., -1., 0.], [0., 0., 1.]],
    [[-1., 0., 0.], [-1., 1., 0.], [0., 0., -1.]],
    [[-1., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
    [[1., -1., 0.], [0., -1., 0.], [0., 0., -1.]],
    [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
    [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
    [[0., 1., 0.], [1., 1., 0.], [0., 0., -1.]],
    [[-1., 0., 0.], [-1., -1., 0.], [0., 0., 1.]],
    [[-1., -1., 0.], [-1., 0., 0.], [0., 0., -1.]],
    [[1., 1., 0.], [0., 1., 0.], [0., 0., 1.]],
];

// change_of_basis matrices for rhombohedral systems.
// C: spacegroup.c:241-252
static CHANGE_OF_BASIS_RHOMBO: [Mat3; 6] = [
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
    [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
    [[0., 0., -1.], [0., -1., 0.], [-1., 0., 0.]],
    [[0., -1., 0.], [-1., 0., 0.], [0., 0., -1.]],
    [[-1., 0., 0.], [0., 0., -1.], [0., -1., 0.]],
];

static CHANGE_OF_BASIS_RHOMBO_HEX: [Mat3; 6] = [
    [[1., 0., 1.], [-1., 1., 1.], [0., -1., 1.]],
    [[0., -1., 1.], [1., 0., 1.], [-1., 1., 1.]],
    [[-1., 1., 1.], [0., -1., 1.], [1., 0., 1.]],
    [[0., 1., -1.], [1., -1., -1.], [-1., 0., -1.]],
    [[1., -1., -1.], [-1., 0., -1.], [0., 1., -1.]],
    [[-1., 0., -1.], [0., 1., -1.], [1., -1., -1.]],
];

/// Rhombohedral Hall numbers with hexagonal (hP) setting.
/// C: match_hall_symbol_db_rhombo, spacegroup.c:1562
const RHOMBO_HEX_SETTING: [i32; 7] = [433, 436, 444, 450, 452, 458, 460];

pub(crate) static SPACEGROUP_TO_HALL_NUMBER: [i32; 230] = [
    1, 2, 3, 6, 9, 18, 21, 30, 39, 57, 60, 63, 72, 81, 90, 108, 109, 112, 115, 116, 119, 122, 123,
    124, 125, 128, 134, 137, 143, 149, 155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209,
    212, 215, 218, 221, 227, 228, 230, 233, 239, 245, 251, 257, 263, 266, 269, 275, 278, 284, 290,
    292, 298, 304, 310, 313, 316, 322, 334, 335, 337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
    355, 356, 357, 358, 359, 361, 363, 364, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376,
    377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395,
    396, 397, 398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 412, 413, 414, 416, 418, 419, 420,
    422, 424, 425, 426, 428, 430, 431, 432, 433, 435, 436, 438, 439, 440, 441, 442, 443, 444, 446,
    447, 448, 449, 450, 452, 454, 455, 456, 457, 458, 460, 462, 463, 464, 465, 466, 467, 468, 469,
    470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488,
    489, 490, 491, 492, 493, 494, 495, 497, 498, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509,
    510, 511, 512, 513, 514, 515, 516, 517, 518, 520, 521, 523, 524, 525, 527, 529, 530,
];

static LAYER_GROUP_TO_HALL_NUMBER: [i32; 80] = [
    -1, -2, -3, -4, -5, -8, -9, -12, -14, -16, -18, -20, -22, -24, -26, -28, -30, -32, -34, -35,
    -37, -38, -39, -40, -42, -43, -44, -46, -48, -50, -52, -54, -56, -58, -60, -62, -64, -65, -67,
    -68, -70, -72, -74, -76, -77, -79, -80, -81, -82, -83, -84, -85, -87, -88, -89, -90, -91, -92,
    -93, -94, -95, -96, -98, -99, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111,
    -112, -113, -114, -115, -116,
];

// --- Structs ---

/// 空间群识别结果。
///
/// 包含标准晶体学信息：空间群编号、Hall 编号、国际符号（短/长/完整）、
/// Schoenflies 符号、Bravais 晶格矩阵和原点平移向量。
///
/// `bravais_lattice` 采用 `[cart][vec]` 布局（行=笛卡尔, 列=Bravais 基矢量）。
#[derive(Clone, Debug)]
pub struct Spacegroup {
    pub number: usize,
    pub hall_number: usize,
    pub pointgroup_number: usize,
    pub schoenflies: String,
    pub hall_symbol: String,
    pub international: String,
    pub international_long: String,
    pub international_short: String,
    pub choice: String,
    pub bravais_lattice: Mat3,
    pub origin_shift: Vec3,
}

impl Spacegroup {
    pub fn new() -> Self {
        Spacegroup {
            number: 0,
            hall_number: 0,
            pointgroup_number: 0,
            schoenflies: String::new(),
            hall_symbol: String::new(),
            international: String::new(),
            international_long: String::new(),
            international_short: String::new(),
            choice: String::new(),
            bravais_lattice: [[0.0; 3]; 3],
            origin_shift: [0.0; 3],
        }
    }
}

// --- Public Functions ---

pub fn spa_search_spacegroup(
    primitive: &Primitive,
    hall_number: i32,
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Spacegroup> {
    debug::debug_print(format_args!(
        "search_spacegroup (tolerance = {}):\n",
        symprec
    ));

    let cell = primitive.cell.as_ref()?;
    let symmetry = sym_get_operation(cell, symprec, angle_tolerance)?;

    let candidates = if hall_number != 0 {
        vec![hall_number]
    } else if cell.aperiodic_axis.is_none() {
        SPACEGROUP_TO_HALL_NUMBER.to_vec()
    } else {
        LAYER_GROUP_TO_HALL_NUMBER.to_vec()
    };

    search_spacegroup_with_symmetry(primitive, &candidates, &symmetry, symprec, angle_tolerance)
}

pub fn spa_search_spacegroup_with_symmetry(
    symmetry: &Symmetry,
    prim_lat: &Mat3,
    symprec: f64,
) -> Option<Spacegroup> {
    let mut primitive = Primitive::new(1);
    let mut cell = Cell::new(1, TensorRank::NoSpin);
    cell.lattice = *prim_lat;
    cell.position[0] = [0.0; 3];
    primitive.cell = Some(cell);

    // Use all spacegroups (0-229)
    let candidates = SPACEGROUP_TO_HALL_NUMBER.to_vec();

    search_spacegroup_with_symmetry(&primitive, &candidates, symmetry, symprec, -1.0)
}

pub fn spa_transform_to_primitive(
    mapping_table: &mut [usize],
    cell: &Cell,
    trans_mat: &Mat3,
    centering: Centering,
    symprec: f64,
) -> Option<Cell> {
    let tmat_inv = match mat_inverse_matrix_d3(trans_mat, symprec) {
        Some(m) => m,
        None => return None,
    };

    let mut tmat = [[0.0; 3]; 3];
    match centering {
        Centering::Primitive => tmat = tmat_inv,
        Centering::AFace => tmat = mat_multiply_matrix_d3(&tmat_inv, &A_MAT),
        Centering::CFace => tmat = mat_multiply_matrix_d3(&tmat_inv, &C_MAT),
        Centering::Face => tmat = mat_multiply_matrix_d3(&tmat_inv, &F_MAT),
        Centering::Body => tmat = mat_multiply_matrix_d3(&tmat_inv, &I_MAT),
        Centering::RCenter => tmat = mat_multiply_matrix_d3(&tmat_inv, &R_MAT),
        _ => return None,
    }

    let prim_lat = mat_multiply_matrix_d3(&cell.lattice, &tmat);

    let primitive = cel_trim_cell(mapping_table, &prim_lat, cell, symprec);
    if primitive.is_none() {
        debug::warning_print(format_args!("spglib: cel_trim_cell failed.\n"));
    }
    primitive
}

pub fn spa_transform_from_primitive(
    primitive: &Cell,
    centering: Centering,
    symprec: f64,
) -> Option<Cell> {
    let mut tmat = [[0.0; 3]; 3];
    let mut inv_tmat = [[0.0; 3]; 3];

    match centering {
        Centering::Primitive => {}
        Centering::AFace => {
            tmat = A_MAT;
            inv_tmat = mat_inverse_matrix_d3(&A_MAT, 0.0).unwrap();
        }
        Centering::CFace => {
            tmat = C_MAT;
            inv_tmat = mat_inverse_matrix_d3(&C_MAT, 0.0).unwrap();
        }
        Centering::Face => {
            tmat = F_MAT;
            inv_tmat = mat_inverse_matrix_d3(&F_MAT, 0.0).unwrap();
        }
        Centering::Body => {
            tmat = I_MAT;
            inv_tmat = mat_inverse_matrix_d3(&I_MAT, 0.0).unwrap();
        }
        Centering::RCenter => {
            tmat = R_MAT;
            inv_tmat = mat_inverse_matrix_d3(&R_MAT, 0.0).unwrap();
        }
        _ => return None,
    }

    let mut shift = [[0.0; 3]; 3];
    let multi = get_centering_shifts(&mut shift, centering);

    let mut mapping_table = vec![0; primitive.size * multi];
    let mut std_cell = Cell::new(primitive.size * multi, primitive.tensor_rank);

    std_cell.lattice = mat_multiply_matrix_d3(&primitive.lattice, &inv_tmat);

    let mut num_atom = 0;
    for i in 0..primitive.size {
        std_cell.position[num_atom] = mat_multiply_matrix_vector_d3(&tmat, &primitive.position[i]);
        std_cell.types[num_atom] = primitive.types[i];
        num_atom += 1;
    }

    for i in 0..(multi - 1) {
        for j in 0..primitive.size {
            let src_pos = std_cell.position[j];
            std_cell.position[num_atom] = src_pos;
            for k in 0..3 {
                std_cell.position[num_atom][k] += shift[i][k];
            }
            std_cell.types[num_atom] = std_cell.types[j];
            num_atom += 1;
        }
    }

    let trimmed_cell = cel_trim_cell(&mut mapping_table, &std_cell.lattice, &std_cell, symprec);
    trimmed_cell
}

// --- Internal Functions ---

fn search_spacegroup_with_symmetry(
    primitive: &Primitive,
    candidates: &[i32],
    symmetry: &Symmetry,
    symprec: f64,
    angle_tolerance: f64,
) -> Option<Spacegroup> {
    debug::debug_print(format_args!(
        "search_spacegroup (tolerance = {}):\n",
        symprec
    ));

    let mut origin_shift = [0.0; 3];
    let mut conv_lattice = [[0.0; 3]; 3];

    let pointsym = ptg_get_pointsymmetry(&symmetry.rot);
    if pointsym.size < symmetry.size {
        debug::info_print(format_args!(
            "spglib: Point symmetry of primitive cell is broken.\n"
        ));
        return None;
    }

    let hall_number = iterative_search_hall_number(
        &mut origin_shift,
        &mut conv_lattice,
        candidates,
        primitive,
        symmetry,
        symprec,
        angle_tolerance,
    );

    if hall_number == 0 {
        return None;
    }

    get_spacegroup(hall_number, &origin_shift, &conv_lattice)
}

fn get_spacegroup(
    hall_number: i32,
    origin_shift: &Vec3,
    conv_lattice: &Mat3,
) -> Option<Spacegroup> {
    let spg_type = spgdb_get_spacegroup_type(hall_number as usize);

    let mut spacegroup = Spacegroup::new();
    spacegroup.bravais_lattice = *conv_lattice;
    spacegroup.origin_shift = *origin_shift;
    spacegroup.number = spg_type.number;
    spacegroup.hall_number = hall_number as usize;
    spacegroup.pointgroup_number = spg_type.pointgroup_number;
    spacegroup.schoenflies = spg_type.schoenflies;
    spacegroup.hall_symbol = spg_type.hall_symbol;
    spacegroup.international = spg_type.international;
    spacegroup.international_long = spg_type.international_full;
    spacegroup.international_short = spg_type.international_short;
    spacegroup.choice = spg_type.choice;

    Some(spacegroup)
}

fn iterative_search_hall_number(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    candidates: &[i32],
    primitive: &Primitive,
    symmetry: &Symmetry,
    symprec: f64,
    angle_tolerance: f64,
) -> i32 {
    debug::debug_print(format_args!("iterative_search_hall_number:\n"));

    let mut hall_number = search_hall_number(
        origin_shift,
        conv_lattice,
        candidates,
        primitive,
        symmetry,
        symprec,
    );

    if hall_number != 0 {
        return hall_number;
    }

    let mut tolerance = symprec;
    let mut current_symmetry = symmetry.clone();

    for attempt in 0..NUM_ATTEMPT {
        debug::debug_print(format_args!(
            "spglib: Attempt {} tolerance = {:e} failed",
            attempt, tolerance
        ));

        tolerance *= REDUCE_RATE;
        if let Some(sym_reduced) = sym_reduce_operation(
            primitive.cell.as_ref().unwrap(),
            &current_symmetry,
            tolerance,
            angle_tolerance,
        ) {
            hall_number = search_hall_number(
                origin_shift,
                conv_lattice,
                candidates,
                primitive,
                &sym_reduced,
                symprec,
            );
            current_symmetry = sym_reduced;
        }

        if hall_number != 0 {
            break;
        }
    }

    hall_number
}

fn search_hall_number(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    candidates: &[i32],
    primitive: &Primitive,
    symmetry: &Symmetry,
    symprec: f64,
) -> i32 {
    debug::debug_print(format_args!("search_hall_number:\n"));

    let aperiodic_axis = primitive.cell.as_ref().unwrap().aperiodic_axis;
    let mut tmat_int = [[0; 3]; 3];

    let pointgroup = ptg_get_transformation_matrix(&mut tmat_int, &symmetry.rot, aperiodic_axis);

    if pointgroup.number == 0 {
        return 0;
    }

    let mut conv_lattice_tmp = [[0.0; 3]; 3];

    if pointgroup.laue == Laue::Laue1 || pointgroup.laue == Laue::Laue2M {
        conv_lattice_tmp =
            mat_multiply_matrix_di3(&primitive.cell.as_ref().unwrap().lattice, &tmat_int);

        if pointgroup.laue == Laue::Laue1 {
            if !change_basis_tricli(
                &mut tmat_int,
                &conv_lattice_tmp,
                &primitive.cell.as_ref().unwrap().lattice,
                symprec,
                aperiodic_axis,
            ) {
                return 0;
            }
        }

        if pointgroup.laue == Laue::Laue2M {
            if !change_basis_monocli(
                &mut tmat_int,
                &conv_lattice_tmp,
                &primitive.cell.as_ref().unwrap().lattice,
                symprec,
                aperiodic_axis,
            ) {
                return 0;
            }
        }
    }

    let mut correction_mat = [[0.0; 3]; 3];
    let centering = get_centering(&mut correction_mat, &tmat_int, pointgroup.laue);
    if centering == Centering::Error {
        return 0;
    }

    let tmat = mat_multiply_matrix_id3(&tmat_int, &correction_mat);
    *conv_lattice = mat_multiply_matrix_d3(&primitive.cell.as_ref().unwrap().lattice, &tmat);

    let conv_symmetry = get_initial_conventional_symmetry(centering, &tmat, symmetry);
    if conv_symmetry.is_none() {
        return 0;
    }
    let conv_symmetry = conv_symmetry.unwrap();

    let holohedry = pointgroup.holohedry;
    let pg_number = pointgroup.number;

    for &cand in candidates {
        if match_hall_symbol_db(
            origin_shift,
            conv_lattice,
            cand,
            pg_number,
            holohedry,
            centering,
            &conv_symmetry,
            symprec,
        ) {
            debug::debug_print(format_args!("origin shift\n"));
            return cand;
        }
    }

    0
}

fn change_basis_tricli(
    tmat_int: &mut Mat3I,
    conv_lattice: &Mat3,
    primitive_lattice: &Mat3,
    symprec: f64,
    aperiodic_axis: Option<AperiodicAxis>,
) -> bool {
    let mut niggli_cell = *conv_lattice;

    if !niggli_reduce(&mut niggli_cell, symprec * symprec, aperiodic_axis) {
        return false;
    }

    let mut smallest_lattice = niggli_cell;
    if mat_get_determinant_d3(&smallest_lattice) < 0.0 {
        for i in 0..3 {
            for j in 0..3 {
                smallest_lattice[i][j] = -smallest_lattice[i][j];
            }
        }
    }

    let inv_lattice = mat_inverse_matrix_d3(primitive_lattice, 0.0).unwrap();
    let tmat = mat_multiply_matrix_d3(&inv_lattice, &smallest_lattice);
    *tmat_int = mat_cast_matrix_3d_to_3i(&tmat);

    true
}

fn change_basis_monocli(
    tmat_int: &mut Mat3I,
    conv_lattice: &Mat3,
    primitive_lattice: &Mat3,
    symprec: f64,
    aperiodic_axis_prim: Option<AperiodicAxis>,
) -> bool {
    let mut smallest_lattice = [[0.0; 3]; 3];
    let mut aperiodic_axis_conv: i32 = -1;
    let unique_axis;

    if aperiodic_axis_prim.is_none() {
        unique_axis = 1;
    } else {
        let ap_idx = aperiodic_axis_prim.unwrap().axis_index();
        for i in 0..3 {
            if tmat_int[ap_idx][i] != 0 {
                aperiodic_axis_conv = i as i32;
            }
        }
        unique_axis = 0;
    }

    if !del_layer_delaunay_reduce_2D(
        &mut smallest_lattice,
        conv_lattice,
        unique_axis,
        aperiodic_axis_conv,
        symprec,
    ) {
        return false;
    }

    if aperiodic_axis_conv == 0 {
        // CHANGE_OF_BASIS_MONOCLI[2] corresponds to index 2 in C array
        smallest_lattice = mat_multiply_matrix_d3(&smallest_lattice, &CHANGE_OF_BASIS_MONOCLI[2]);
    }

    let inv_lattice = mat_inverse_matrix_d3(primitive_lattice, 0.0).unwrap();
    let tmat = mat_multiply_matrix_d3(&inv_lattice, &smallest_lattice);
    *tmat_int = mat_cast_matrix_3d_to_3i(&tmat);

    true
}

pub(crate) fn get_initial_conventional_symmetry(
    centering: Centering,
    tmat: &Mat3,
    symmetry: &Symmetry,
) -> Option<Symmetry> {
    debug::debug_print(format_args!("get_initial_conventional_symmetry\n"));

    if centering == Centering::RCenter {
        get_conventional_symmetry(tmat, Centering::Primitive, symmetry)
    } else {
        get_conventional_symmetry(tmat, centering, symmetry)
    }
}

fn get_conventional_symmetry(
    tmat: &Mat3,
    centering: Centering,
    primitive_sym: &Symmetry,
) -> Option<Symmetry> {
    let size = primitive_sym.size;
    let mut symmetry;

    match centering {
        Centering::Face => symmetry = Symmetry::new(size * 4),
        Centering::RCenter => symmetry = Symmetry::new(size * 3),
        Centering::Body | Centering::AFace | Centering::BFace | Centering::CFace => {
            symmetry = Symmetry::new(size * 2)
        }
        _ => symmetry = Symmetry::new(size),
    }

    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0).unwrap_or([[0.0; 3]; 3]);

    for i in 0..size {
        let primitive_sym_rot_d3 = mat_cast_matrix_3i_to_3d(&primitive_sym.rot[i]);

        // C*S*C^-1
        let mut symmetry_rot_d3 = [[0.0; 3]; 3];
        if let Some(res) = mat_get_similar_matrix_d3(&primitive_sym_rot_d3, tmat, 0.0) {
            symmetry_rot_d3 = res;
        }
        symmetry.rot[i] = mat_cast_matrix_3d_to_3i(&symmetry_rot_d3);

        // translation: C = B^-1*P
        symmetry.trans[i] = mat_multiply_matrix_vector_d3(&inv_tmat, &primitive_sym.trans[i]);
    }

    let mut shift = [[0.0; 3]; 3];
    let multi = get_centering_shifts(&mut shift, centering);

    if centering != Centering::Primitive {
        for i in 0..(multi - 1) {
            for j in 0..size {
                //mat_copy_matrix_i3(&mut symmetry.rot[(i + 1) * size + j], &symmetry.rot[j]);
                let src_rot = symmetry.rot[j];
                symmetry.rot[(i + 1) * size + j] = src_rot;
                for k in 0..3 {
                    symmetry.trans[(i + 1) * size + j][k] = symmetry.trans[j][k] + shift[i][k];
                }
            }
        }
    }

    Some(symmetry)
}

pub(crate) fn get_centering(correction_mat: &mut Mat3, tmat: &Mat3I, laue: Laue) -> Centering {
    *correction_mat = IDENTITY;
    let det = mat_get_determinant_i3(tmat).abs();
    debug::debug_print(format_args!("laue class: {:?}\n", laue));
    debug::debug_print(format_args!("multiplicity: {}\n", det));

    match det {
        1 => Centering::Primitive,
        2 => {
            let mut centering = get_base_center(tmat);
            if centering == Centering::AFace {
                if laue == Laue::Laue2M {
                    debug::debug_print(format_args!("Monocli A to C\n"));
                    *correction_mat = MONOCLI_A2C;
                } else {
                    *correction_mat = A2C;
                }
                centering = Centering::CFace;
            }
            if centering == Centering::BFace {
                *correction_mat = B2C;
                centering = Centering::CFace;
            }
            if laue == Laue::Laue2M && centering == Centering::Body {
                debug::debug_print(format_args!("Monocli I to C\n"));
                *correction_mat = MONOCLI_I2C;
                centering = Centering::CFace;
            }
            centering
        }
        3 => {
            let mut trans_corr_mat = [[0.0; 3]; 3];
            trans_corr_mat = mat_multiply_matrix_id3(tmat, &RHOMBO_OBVERSE);
            if mat_is_int_matrix(&trans_corr_mat, INT_PREC) {
                *correction_mat = RHOMBO_OBVERSE;
                debug::debug_print(format_args!("R-center observe setting\n"));
                return Centering::RCenter;
            }
            trans_corr_mat = mat_multiply_matrix_id3(tmat, &RHOMBO_REVERSE);
            if mat_is_int_matrix(&trans_corr_mat, INT_PREC) {
                *correction_mat = RHOMBO_REVERSE;
                debug::debug_print(format_args!("R-center reverse setting\n"));
                return Centering::RCenter;
            }
            Centering::RCenter
        }
        4 => Centering::Face,
        _ => Centering::Error,
    }
}

fn get_base_center(tmat: &Mat3I) -> Centering {
    debug::debug_print(format_args!("lat_get_base_center\n"));

    // C center
    for i in 0..3 {
        if tmat[i][0] == 0 && tmat[i][1] == 0 && tmat[i][2].abs() == 1 {
            return Centering::CFace;
        }
    }

    // A center
    for i in 0..3 {
        if tmat[i][0].abs() == 1 && tmat[i][1] == 0 && tmat[i][2] == 0 {
            return Centering::AFace;
        }
    }

    // B center
    for i in 0..3 {
        if tmat[i][0] == 0 && tmat[i][1].abs() == 1 && tmat[i][2] == 0 {
            return Centering::BFace;
        }
    }

    // Body center
    if (tmat[0][0].abs() + tmat[0][1].abs() + tmat[0][2].abs() == 2)
        && (tmat[1][0].abs() + tmat[1][1].abs() + tmat[1][2].abs() == 2)
        && (tmat[2][0].abs() + tmat[2][1].abs() + tmat[2][2].abs() == 2)
    {
        return Centering::Body;
    }

    debug::warning_print(format_args!("spglib: No centring was found.\n"));
    Centering::Primitive
}

fn get_centering_shifts(shift: &mut [[f64; 3]; 3], centering: Centering) -> usize {
    let mut multi = 1;
    for i in 0..3 {
        shift[i] = [0.0; 3];
    }

    if centering != Centering::Primitive {
        if centering != Centering::Face && centering != Centering::RCenter {
            for i in 0..3 {
                shift[0][i] = 0.5;
            } // BODY
            if centering == Centering::AFace {
                shift[0][0] = 0.0;
            }
            if centering == Centering::BFace {
                shift[0][1] = 0.0;
            }
            if centering == Centering::CFace {
                shift[0][2] = 0.0;
            }
            multi = 2;
        }

        if centering == Centering::RCenter {
            shift[0][0] = 2.0 / 3.0;
            shift[0][1] = 1.0 / 3.0;
            shift[0][2] = 1.0 / 3.0;
            shift[1][0] = 1.0 / 3.0;
            shift[1][1] = 2.0 / 3.0;
            shift[1][2] = 2.0 / 3.0;
            multi = 3;
        }

        if centering == Centering::Face {
            shift[0][0] = 0.0;
            shift[0][1] = 0.5;
            shift[0][2] = 0.5;
            shift[1][0] = 0.5;
            shift[1][1] = 0.0;
            shift[1][2] = 0.5;
            shift[2][0] = 0.5;
            shift[2][1] = 0.5;
            shift[2][2] = 0.0;
            multi = 4;
        }
    }
    multi
}

fn is_equivalent_lattice(
    tmat: &mut Mat3,
    mode: i32,
    lattice: &Mat3,
    orig_lattice: &Mat3,
    symprec: f64,
) -> bool {
    if (mat_get_determinant_d3(lattice) - mat_get_determinant_d3(orig_lattice)).abs() > symprec {
        return false;
    }

    let mut inv_lat = [[0.0; 3]; 3];
    if mat_inverse_matrix_d3(lattice, symprec).is_none() {
        return false;
    }

    // orig_lattice == lattice @ tmat
    *tmat = mat_multiply_matrix_d3(&inv_lat, orig_lattice);

    match mode {
        0 => {
            if mat_check_identity_matrix_d3(&IDENTITY, tmat, symprec) {
                return true;
            }
        }
        1 => {
            let mut tmat_abs = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    tmat_abs[i][j] = mat_dabs(tmat[i][j]);
                }
            }
            if mat_check_identity_matrix_d3(&IDENTITY, &tmat_abs, symprec) {
                return true;
            }
        }
        2 => {
            let tmat_int = mat_cast_matrix_3d_to_3i(tmat);
            if mat_check_identity_matrix_id3(&tmat_int, tmat, symprec) {
                if mat_get_determinant_i3(&tmat_int) == 1 {
                    let orig_metric = mat_get_metric(orig_lattice);
                    let metric = mat_get_metric(lattice);
                    if mat_check_identity_matrix_d3(&orig_metric, &metric, symprec) {
                        return true;
                    }
                }
            }
        }
        _ => {}
    }

    false
}

// ============================================================================
// match_hall_symbol_db family — ported from C spacegroup.c:991-1679
// ============================================================================

/// Generic change-of-basis loop over candidate orientation matrices.
/// Expands symmetry with centering translations when `centering == RCenter`.
/// C: `match_hall_symbol_db_change_of_basis_loop`, spacegroup.c:1610
fn match_hall_symbol_db_change_of_basis_loop(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    change_of_basis: &[Mat3],
    hall_number: i32,
    centering: Centering,
    conv_symmetry: &Symmetry,
    symprec: f64,
) -> bool {
    let centering_for_symmetry = if centering == Centering::RCenter {
        Centering::RCenter
    } else {
        Centering::Primitive
    };

    for cob in change_of_basis {
        // Expand symmetry with centering translations.
        // RCenter: 12 ops → 36 ops via get_conventional_symmetry.
        let Some(changed_sym) =
            get_conventional_symmetry(cob, centering_for_symmetry, conv_symmetry)
        else {
            continue;
        };
        let mut changed_lat = mat_multiply_matrix_d3(&*conv_lattice, cob);
        // Pass original `centering` (not `centering_for_symmetry`) to matching:
        // for RCenter this means the DB expects RCenter operations.
        if hal_match_hall_symbol_db(
            origin_shift,
            &mut changed_lat,
            hall_number,
            centering,
            &changed_sym,
            symprec,
        ) {
            *conv_lattice = changed_lat;
            return true;
        }
    }
    false
}

/// Rhombohedral Hall symbol matching.
/// Handles both hexagonal (hP) and rhombohedral (a=b=c) settings.
/// C: `match_hall_symbol_db_rhombo`, spacegroup.c:1554
fn match_hall_symbol_db_rhombo(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    hall_number: i32,
    conv_symmetry: &Symmetry,
    symprec: f64,
) -> bool {
    if RHOMBO_HEX_SETTING.contains(&hall_number) {
        // Hexagonal setting: use rhombo_hex change-of-basis + RCenter expansion.
        // C: spacegroup.c:1565-1568
        match_hall_symbol_db_change_of_basis_loop(
            origin_shift,
            conv_lattice,
            &CHANGE_OF_BASIS_RHOMBO_HEX,
            hall_number,
            Centering::RCenter,
            conv_symmetry,
            symprec,
        )
    } else {
        // Rhombohedral (a=b=c) setting: primitive centering.
        // C: spacegroup.c:1576-1578
        match_hall_symbol_db_change_of_basis_loop(
            origin_shift,
            conv_lattice,
            &CHANGE_OF_BASIS_RHOMBO,
            hall_number,
            Centering::Primitive,
            conv_symmetry,
            symprec,
        )
    }
}

/// Top-level Hall symbol dispatcher with point-group filter and holohedry routing.
/// C: `match_hall_symbol_db`, spacegroup.c:991
fn match_hall_symbol_db(
    origin_shift: &mut Vec3,
    conv_lattice: &mut Mat3,
    hall_number: i32,
    pointgroup_number: usize,
    holohedry: Holohedry,
    centering: Centering,
    conv_symmetry: &Symmetry,
    symprec: f64,
) -> bool {
    // Point-group filter: skip Hall numbers whose point group doesn't match.
    let spg_type = spgdb_get_spacegroup_type(hall_number as usize);
    if pointgroup_number != spg_type.pointgroup_number {
        return false;
    }

    // Dispatch: rhombohedral RCenter cases get centering expansion;
    // all others delegate to hal_match_hall_symbol_db directly (as before).
    if holohedry == Holohedry::Trigo
        && centering == Centering::RCenter
        && RHOMBO_HEX_SETTING.contains(&hall_number)
    {
        match_hall_symbol_db_rhombo(origin_shift, conv_lattice, hall_number, conv_symmetry, symprec)
    } else {
        hal_match_hall_symbol_db(
            origin_shift,
            conv_lattice,
            hall_number,
            centering,
            conv_symmetry,
            symprec,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;

    const SYMPREC: f64 = 1e-5;
    const ANGLE_TOL: f64 = -1.0;

    /// 辅助函数：根据晶格和原子位置搜索空间群
    /// 流程：Cell → prm_get_primitive → spa_search_spacegroup
    fn search_spacegroup(
        lattice: &Mat3,
        positions: &[Vec3],
        types: &[i32],
    ) -> Option<Spacegroup> {
        let mut cell = Cell::new(types.len(), crate::cell::TensorRank::NoSpin);
        cell.set_cell(lattice, positions, types);
        cell.aperiodic_axis = None;

        // 第一步：获取原胞 (Primitive)
        let primitive = crate::primitive::prm_get_primitive(&cell, SYMPREC, ANGLE_TOL)?;

        // 第二步：搜索空间群
        spa_search_spacegroup(&primitive, 0, SYMPREC, ANGLE_TOL)
    }

    /// 验证空间群识别结果
    fn assert_spacegroup(
        spg: &Spacegroup,
        expected_number: usize,
        expected_short: &str,
        label: &str,
    ) {
        assert_eq!(
            spg.number, expected_number,
            "[{}] Expected space group #{}, got #{}",
            label, expected_number, spg.number
        );
        let actual_short = spg.international_short.trim();
        assert_eq!(
            actual_short, expected_short,
            "[{}] Expected short symbol '{}', got '{}'",
            label, expected_short, actual_short
        );
        // 验证 bravais_lattice 非零
        let det = mat_get_determinant_d3(&spg.bravais_lattice);
        assert!(det.abs() > 1e-10, "[{}] bravais_lattice determinant is zero", label);
    }

    // ==================== 立方晶系测试 ====================

    #[test]
    fn test_cubic_simple_pm3m() {
        // 简单立方晶格，1个原子在原点 → Pm-3m (#221)
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0]];
        let types = [1];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("Simple cubic: spacegroup search returned None");
        assert_spacegroup(&spg, 221, "Pm-3m", "simple cubic");
    }

    #[test]
    fn test_cubic_bcc_im3m() {
        // 体心立方：立方晶格 + 2个原子 → Im-3m (#229)
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let types = [1, 1];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("BCC: spacegroup search returned None");
        assert_spacegroup(&spg, 229, "Im-3m", "BCC");
    }

    #[test]
    fn test_cubic_fcc_fm3m() {
        // 面心立方：立方晶格 + 4个原子 → Fm-3m (#225)
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ];
        let types = [1, 1, 1, 1];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("FCC: spacegroup search returned None");
        assert_spacegroup(&spg, 225, "Fm-3m", "FCC");
    }

    #[test]
    fn test_cubic_diamond_fd3m() {
        // 金刚石结构：FCC + 2原子基元 → Fd-3m (#227)
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25],
            [0.0, 0.5, 0.5],
            [0.25, 0.75, 0.75],
            [0.5, 0.0, 0.5],
            [0.75, 0.25, 0.75],
            [0.5, 0.5, 0.0],
            [0.75, 0.75, 0.25],
        ];
        let types = [1, 1, 1, 1, 1, 1, 1, 1];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("Diamond: spacegroup search returned None");
        assert_spacegroup(&spg, 227, "Fd-3m", "diamond");
    }

    // ==================== 六方晶系测试 ====================

    /// 构建六方晶格矩阵（a=b, c 独立，γ=120°）
    /// 六方晶格在 cart×vec 约定下：每行是一个 Cartesian 分量，每列是一个晶格向量
    /// a = (a, 0, 0), b = (-a/2, a*sqrt(3)/2, 0), c = (0, 0, c)
    fn hexagonal_lattice(a: f64, c: f64) -> Mat3 {
        let sqrt3 = 3.0_f64.sqrt();
        [
            [a, -a / 2.0, 0.0],
            [0.0, a * sqrt3 / 2.0, 0.0],
            [0.0, 0.0, c],
        ]
    }

    #[test]
    fn test_graphene_p6mmm() {
        // 石墨烯：六方晶格，2个C原子，平面结构 → P6/mmm (#191)
        let a = 2.46;
        let c = 10.0;
        let lattice = hexagonal_lattice(a, c);
        // 蜂巢结构的两个原子：A 和 B 位点
        let positions = [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.0]];
        let types = [6, 6]; // 碳原子

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("Graphene: spacegroup search returned None");
        assert_spacegroup(&spg, 191, "P6/mmm", "graphene");
    }

    #[test]
    fn test_silicene_p3m1() {
        // 硅烯：六方晶格，2个Si原子，low-buckled结构 → P-3m1 (#164)
        // AB子晶格在z方向作相反方向buckling，保留了反演对称性但打破了水平镜面
        let a = 3.86;
        let c = 20.0;
        let lattice = hexagonal_lattice(a, c);
        let delta = 0.44 / c; // buckling 在分数坐标中
        let positions = [[1.0 / 3.0, 2.0 / 3.0, -delta], [2.0 / 3.0, 1.0 / 3.0, delta]];
        let types = [14, 14]; // 硅原子

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("Silicene: spacegroup search returned None");
        // D6h → D3d（保留3重轴和反演，失去6重轴和水平镜面）
        assert_spacegroup(&spg, 164, "P-3m1", "silicene");
    }

    #[test]
    fn test_hcp_p63mmc() {
        // HCP 结构：六方晶格，2个原子 → P6_3/mmc (#194)
        let a = 2.5;
        let c = a * (8.0_f64 / 3.0_f64).sqrt(); // 理想 c/a 比
        let lattice = hexagonal_lattice(a, c);
        let positions = [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.5]];
        let types = [1, 1];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("HCP: spacegroup search returned None");
        assert_spacegroup(&spg, 194, "P6_3/mmc", "HCP");
    }

    // ==================== 超胞测试 ====================

    #[test]
    fn test_supercell_2x2x2_simple_cubic() {
        // 2x2x2 超胞的简单立方 → 应仍返回 Pm-3m (#221)
        let a_super = 2.0;
        let lattice = [
            [a_super, 0.0, 0.0],
            [0.0, a_super, 0.0],
            [0.0, 0.0, a_super],
        ];
        // 8个原子在超胞的各个角上
        let mut positions = Vec::new();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    positions.push([
                        i as f64 * 0.5,
                        j as f64 * 0.5,
                        k as f64 * 0.5,
                    ]);
                }
            }
        }
        let types = vec![1; 8];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("Supercell 2x2x2: spacegroup search returned None");
        assert_spacegroup(&spg, 221, "Pm-3m", "supercell 2x2x2 simple cubic");
    }

    #[test]
    fn test_supercell_2x2x2_cscl() {
        // 2x2x2 超胞的 CsCl 结构 → 应仍返回 Pm-3m (#221)
        let a_super = 2.0;
        let lattice = [
            [a_super, 0.0, 0.0],
            [0.0, a_super, 0.0],
            [0.0, 0.0, a_super],
        ];
        // CsCl: 2种原子，Cs在角上，Cl在体心
        // 2x2x2 超胞：8个Cs + 8个Cl = 16个原子
        let mut positions = Vec::new();
        let mut types = Vec::new();
        // Cs 原子在角上 (type=1)
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    positions.push([
                        i as f64 * 0.5,
                        j as f64 * 0.5,
                        k as f64 * 0.5,
                    ]);
                    types.push(1);
                }
            }
        }
        // Cl 原子在体心位置 (type=2)
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    positions.push([
                        i as f64 * 0.5 + 0.25,
                        j as f64 * 0.5 + 0.25,
                        k as f64 * 0.5 + 0.25,
                    ]);
                    types.push(2);
                }
            }
        }

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("Supercell 2x2x2 CsCl: spacegroup search returned None");
        assert_spacegroup(&spg, 221, "Pm-3m", "supercell 2x2x2 CsCl");
    }

    #[test]
    fn test_supercell_2x2x1_graphene() {
        // 2x2x1 超胞的石墨烯 → 应仍返回 P6/mmm (#191)
        let a = 2.46;
        let c = 10.0;
        // 超胞: 2a x 2b x c
        let super_lattice = [
            [2.0 * a, -a, 0.0],
            [0.0, a * 3.0_f64.sqrt(), 0.0],
            [0.0, 0.0, c],
        ];
        // 原始石墨烯的2个原子 × (2×2×1) = 8个原子
        let base_positions = [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.0]];
        let mut positions = Vec::new();
        let mut types = Vec::new();
        for i in 0..2 {
            for j in 0..2 {
                for (base_pos, &base_type) in base_positions.iter().zip([6, 6].iter()) {
                    positions.push([
                        (base_pos[0] + i as f64) / 2.0,
                        (base_pos[1] + j as f64) / 2.0,
                        base_pos[2],
                    ]);
                    types.push(base_type);
                }
            }
        }

        let spg = search_spacegroup(&super_lattice, &positions, &types)
            .expect("Supercell 2x2x1 graphene: spacegroup search returned None");
        assert_spacegroup(&spg, 191, "P6/mmm", "supercell 2x2x1 graphene");
    }

    // ==================== 其他晶系测试 ====================

    #[test]
    fn test_rocksalt_fm3m() {
        // NaCl 岩盐结构 → Fm-3m (#225)
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [
            [0.0, 0.0, 0.0],   // Na
            [0.5, 0.5, 0.5],   // Cl
            [0.0, 0.5, 0.5],   // Na
            [0.5, 0.0, 0.5],   // Na
            [0.5, 0.5, 0.0],   // Na
            [0.0, 0.0, 0.5],   // Cl
            [0.0, 0.5, 0.0],   // Cl
            [0.5, 0.0, 0.0],   // Cl
        ];
        let types = [1, 2, 1, 1, 1, 2, 2, 2];

        let spg = search_spacegroup(&lattice, &positions, &types)
            .expect("NaCl: spacegroup search returned None");
        assert_spacegroup(&spg, 225, "Fm-3m", "NaCl");
    }

}
