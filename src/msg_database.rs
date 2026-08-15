//! 磁性空间群数据库。
//!
//! 包含磁性空间群类型、对称操作、Hall 映射等数据以及查询函数。
//! 数据部分由 `msg_database_gen.rs` 提供（自动转换自 C）。
//! 逻辑函数手动移植自 `msg_database.txt`。

include!("msg_database_gen.rs");

use crate::debug;
use crate::mathfunc::Mat3I;
use crate::symmetry::{MagneticSymmetry, Symmetry};
use crate::MagneticType;

/// 获取磁性空间群类型。
pub fn get_magnetic_spacegroup_type(uni_number: usize) -> &'static MagneticSpacegroupType {
    if uni_number > 0 && uni_number < MAGNETIC_SPACEGROUP_TYPES.len() {
        &MAGNETIC_SPACEGROUP_TYPES[uni_number]
    } else {
        &MAGNETIC_SPACEGROUP_TYPES[0]
    }
}

/// 返回给定 Hall 编号的最小和最大 UNI 编号。
/// `uni_number_range[0]` = 最小, `[1]` = 最大。
pub fn get_uni_candidates(hall_number: usize) -> Option<[usize; 2]> {
    if hall_number >= MAGNETIC_SPACEGROUP_HALL_MAPPING.len() {
        return None;
    }
    let raw = MAGNETIC_SPACEGROUP_HALL_MAPPING[hall_number];
    Some([raw[0] as usize, raw[1] as usize])
}

/// 获取磁性对称操作。
pub fn get_spacegroup_operations(
    uni_number: usize,
    hall_number: usize,
) -> Option<MagneticSymmetry> {
    let hall_number_offset = get_hall_number_offset(uni_number, hall_number)?;
    let order = MAGNETIC_SPACEGROUP_OPERATION_INDEX[uni_number][hall_number_offset][0];
    let start = MAGNETIC_SPACEGROUP_OPERATION_INDEX[uni_number][hall_number_offset][1];

    debug::debug_print(format_args!(
        "Load operations with UNI={}, hall={}: hall_number_offset={} order={} starts={}\n",
        uni_number, hall_number, hall_number_offset, order, start
    ));

    let mut sym = MagneticSymmetry::new(order as usize);
    for i in 0..order {
        let idx = (start + i) as usize;
        let (rot, trans, timerev) = decode_magnetic_operation(idx);
        sym.rot[i as usize] = rot;
        sym.trans[i as usize] = trans;
        sym.timerev[i as usize] = timerev;
    }

    Some(sym)
}

/// 获取标准变换。
pub fn get_std_transformations(
    uni_number: usize,
    hall_number: usize,
) -> Option<Symmetry> {
    let hall_number_offset = get_hall_number_offset(uni_number, hall_number)?;
    let mut sym = Symmetry::new(7);
    // Identity transformation as first element
    sym.rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    sym.trans[0] = [0.0; 3];

    for (i, &enc) in ALTERNATIVE_TRANSFORMATIONS[uni_number][hall_number_offset]
        .iter()
        .take(7)
        .enumerate()
    {
        if enc == 0 {
            sym.truncate(i + 1);
            break;
        }
        let (tmat, origin_shift) = crate::spg_database::decode_symmetry(enc);
        sym.rot[i + 1] = tmat;
        sym.trans[i + 1] = origin_shift;
    }

    Some(sym)
}

/// 解码磁性对称操作。
fn decode_magnetic_operation(op_number: usize) -> (Mat3I, [f64; 3], bool) {
    let enc = MAGNETIC_SYMMETRY_OPERATIONS[op_number];
    // timerev=true for anti operation, false for ordinary operation
    // 34012224 = 3^9 * 12^3
    let timerev = (enc / 34012224) != 0;
    let (rot, trans) = crate::spg_database::decode_symmetry(enc % 34012224);
    (rot, trans, timerev)
}

/// 获取 Hall 编号偏移量。
fn get_hall_number_offset(uni_number: usize, hall_number: usize) -> Option<usize> {
    if uni_number == 0 || uni_number > 1651 {
        return None;
    }

    let num_halls = MAGNETIC_SPACEGROUP_UNI_MAPPING[uni_number][0] as usize;
    let first_hall = MAGNETIC_SPACEGROUP_UNI_MAPPING[uni_number][1] as usize;

    let hall_number_offset = if hall_number > 0 && hall_number <= 530 {
        if hall_number < first_hall {
            return None;
        }
        hall_number - first_hall
    } else if hall_number == 0 {
        0
    } else {
        return None;
    };

    if hall_number_offset >= num_halls {
        return None;
    }

    Some(hall_number_offset)
}
