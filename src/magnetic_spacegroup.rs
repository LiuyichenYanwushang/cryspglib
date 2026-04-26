//! 磁性空间群识别。
//!
//! 使用磁性对称操作数据库识别磁性空间群类型。
//! 参考: Litvin, "Magnetic Group Tables" 2013

use crate::MagneticType;
use crate::mathfunc::{
    mat_cast_matrix_3i_to_3d, mat_check_identity_matrix_i3,
    mat_get_determinant_d3, mat_inverse_matrix_d3, mat_multiply_matrix_d3,
    mat_multiply_matrix_i3, mat_multiply_matrix_id3,
    mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_id3, mat_nint, Mat3, Mat3I, Vec3,
};
use crate::msg_database::{
    msgdb_get_magnetic_spacegroup_type, msgdb_get_spacegroup_operations,
    msgdb_get_std_transformations, msgdb_get_uni_candidates,
};
use crate::hall_symbol::hal_match_hall_symbol_db;
use crate::pointgroup::ptg_get_transformation_matrix;
use crate::primitive::prm_get_primitive_symmetry;
use crate::refinement::ref_find_similar_bravais_lattice;
use crate::spacegroup::{
    get_centering, get_initial_conventional_symmetry,
    spa_search_spacegroup_with_symmetry, Spacegroup,
};
use crate::spg_database::{spgdb_get_spacegroup_type, Centering};
use crate::symmetry::{MagneticSymmetry, Symmetry};

const MAX_DENOMINATOR: f64 = 100.0;

/// 磁性空间群识别结果的中间数据结构。
pub struct MagneticDataset {
    pub uni_number: usize,
    pub msg_type: MagneticType,
    pub hall_number: usize,
    pub transformation_matrix: Mat3,
    pub origin_shift: Vec3,
    pub std_rotation_matrix: Mat3,
}

/// 识别磁性空间群类型。
///
/// 给定晶格和磁性对称操作，返回识别出的磁性数据集。
/// 如果标准路径失败（磁对称性为非完整空间群），在 `parent_hall_number` 提供时
/// 使用非磁母空间群的 Hall 编号作为 fallback。
pub fn msg_identify_magnetic_space_group_type(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<MagneticDataset> {
    msg_identify_with_parent_hall(lattice, magnetic_symmetry, None, symprec)
}

/// 与 [`msg_identify_magnetic_space_group_type`] 相同，但可指定非磁母空间群的
/// Hall 编号作为 fallback，适用于磁对称性显著低于母空间群的场景。
pub fn msg_identify_with_parent_hall(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: Option<usize>,
    symprec: f64,
) -> Option<MagneticDataset> {
    // 先约化到原胞表示（去除纯平移导致的冗余操作）
    let prim_sym = reduce_to_primitive_magsym(magnetic_symmetry, symprec);

    // 标准路径: 从磁对称性中提取 FSG/XSG 并搜索空间群
    let (ref_sg, changed_symmetry, mut tmat, mut shift, msgtype_num) =
        match get_reference_space_group(&prim_sym, symprec) {
            Some(result) => {
                result
            }
            None => {
                // 标准路径失败 → 尝试 fallback（用母空间群的 Hall 编号）
                let hall_number = parent_hall_number?;
                let fb = build_fallback_reference(lattice, magnetic_symmetry, hall_number, symprec)?;
                fb
            }
        };

    // 对母空间群 Hall 编号和 FSG Hall 编号分别尝试 UNI 匹配。
    // 优先使用母空间群的 Hall（因为磁空间群是母空间群的子群）。
    let mut hall_numbers_try = vec![];
    if let Some(ph) = parent_hall_number {
        hall_numbers_try.push(ph);
    }
    if !hall_numbers_try.contains(&ref_sg.hall_number) {
        hall_numbers_try.push(ref_sg.hall_number);
    }

    let mut best_uni = 0;
    let mut best_msg_type = MagneticType::NonMagnetic;
    let mut best_hall_number = 0;

    for &hall_number in &hall_numbers_try {
        let range = match msgdb_get_uni_candidates(hall_number) {
            Some(r) => r,
            None => continue,
        };
        let min_uni = range[0];
        let max_uni = range[1];

        for uni_number in min_uni..=max_uni {
            let msgtype_db = msgdb_get_magnetic_spacegroup_type(uni_number);
            if msgtype_db.type_ != msgtype_num {
                continue;
            }
            let msg_uni = match msgdb_get_spacegroup_operations(uni_number, hall_number) {
                Some(u) => u,
                None => continue,
            };
            if changed_symmetry.size > msg_uni.size {
                continue;
            }

            let transformations = match msgdb_get_std_transformations(uni_number, hall_number) {
                Some(t) => t,
                None => continue,
            };

            let mut same = false;
            for trans_idx in 0..transformations.size {
                let tmat_cor = transformations.rot[trans_idx];
                let shift_cor = transformations.trans[trans_idx];

                let tmat_cor_d = mat_cast_matrix_3i_to_3d(&tmat_cor);

                let symmetry_cor = get_distinct_changed_magnetic_symmetry(
                    &tmat_cor_d, &shift_cor, &changed_symmetry,
                );

                let symmetry_cor = match symmetry_cor {
                    Some(s) => s,
                    None => continue,
                };

                let matched = is_subset(&symmetry_cor, &msg_uni, symprec);
                if matched {
                    same = true;
                    tmat = mat_multiply_matrix_d3(&tmat_cor_d, &tmat);
                    let shift_tmp = mat_multiply_matrix_vector_d3(&tmat_cor_d, &shift);
                    for s in 0..3 {
                        shift[s] = shift_tmp[s] + shift_cor[s];
                    }
                    break;
                }
            }

            if same {
                best_uni = uni_number;
                best_msg_type = msgtype_db.type_;
                best_hall_number = hall_number;
                break;
            }
        }

        if best_uni != 0 {
            break;
        }
    }

    if best_uni == 0 {
        return None;
    }

    let hall_number = best_hall_number;

    let _msgtype = msgdb_get_magnetic_spacegroup_type(best_uni);
    let mut ret = MagneticDataset {
        uni_number: best_uni,
        msg_type: best_msg_type,
        hall_number,
        transformation_matrix: [[0.0; 3]; 3],
        origin_shift: [0.0; 3],
        std_rotation_matrix: [[0.0; 3]; 3],
    };

    get_rigid_rotation(&mut ret.std_rotation_matrix, lattice, &tmat, &ref_sg);
    ret.transformation_matrix = tmat;
    ret.origin_shift = shift;

    Some(ret)
}

/// 获取参考空间群和变换后的磁性对称操作。
/// 从磁性对称中获取参考空间群、变换后的磁性对称操作、变换矩阵和类型。
///
/// 对应 C 原版的 `get_reference_space_group`。
fn get_reference_space_group(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<(Spacegroup, MagneticSymmetry, Mat3, Vec3, MagneticType)> {
    // 1. FSG = 所有操作忽略时间反演 → 空间群搜索
    let (mut fsg, sym_fsg) =
        match get_family_space_group_with_magnetic_symmetry(magnetic_symmetry, symprec) {
            Some(r) => r,
            None => {
                return None;
            }
        };

    // 2. XSG = 仅 timerev=0 操作 → 空间群搜索 (含完整的 Symmetry 用于 factor group)
    let (mut _xsg, sym_xsg) =
        match get_maximal_subspace_group_with_magnetic_symmetry(magnetic_symmetry, symprec) {
            Some(r) => r,
            None => {
                return None;
            }
        };

    // 3. 确定 MSG 类型 + 获取代表元
    let msgtype_num =
        get_magnetic_space_group_type(magnetic_symmetry, sym_fsg.size, sym_xsg.size, symprec)?;
    let representatives = build_representatives(msgtype_num, magnetic_symmetry)?;

    // 4. 选择参考设置: type-4 用 XSG, 其他用 FSG
    //    C 原版对 type-4 使用 xsg 作为 ref_sg
    let ref_sg = if msgtype_num == MagneticType::AntiTranslation {
        &mut _xsg
    } else {
        &mut fsg
    };

    // 5. 计算 tmat 和 shift
    let tmat = ref_sg.bravais_lattice;
    let shift = ref_sg.origin_shift;

    // 6. 合成变换后的磁性对称操作
    //    (C 原版: get_changed_magnetic_symmetry 分解 + 重合成)
    let changed_symmetry = get_changed_magnetic_symmetry(
        &tmat, &shift, &representatives, &sym_xsg,
        magnetic_symmetry, symprec,
    )?;

    // 7. 复制 ref_sg 用于返回
    let mut ref_sg_copy = Spacegroup::new();
    ref_sg_copy = ref_sg.clone();

    Some((ref_sg_copy, changed_symmetry, tmat, shift, msgtype_num))
}

/// Fallback when `get_reference_space_group` fails.
/// 使用提供的非磁 Hall 编号构建参考空间群，绕过空间群搜索。
fn build_fallback_reference(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: usize,
    symprec: f64,
) -> Option<(Spacegroup, MagneticSymmetry, Mat3, Vec3, MagneticType)> {
    // 1. 提取 FSG/XSG symmetry (只取操作，不搜索空间群)
    let sym_fsg = extract_symmetry(magnetic_symmetry, true, symprec)?;
    let sym_xsg = extract_symmetry(magnetic_symmetry, false, symprec)?;

    // 2. 确定磁性类型
    let msgtype_num =
        get_magnetic_space_group_type(magnetic_symmetry, sym_fsg.size, sym_xsg.size, symprec)?;

    // 3. 用非磁 Hall 编号构建参考 Spacegroup
    let spg_type = spgdb_get_spacegroup_type(parent_hall_number);
    let mut ref_sg = Spacegroup::new();
    ref_sg.hall_number = parent_hall_number;
    ref_sg.number = spg_type.number;
    ref_sg.pointgroup_number = spg_type.pointgroup_number;
    ref_sg.schoenflies = spg_type.schoenflies;
    ref_sg.hall_symbol = spg_type.hall_symbol;
    ref_sg.international = spg_type.international;
    ref_sg.international_long = spg_type.international_full;
    ref_sg.international_short = spg_type.international_short;
    ref_sg.choice = spg_type.choice;
    ref_sg.bravais_lattice = *lattice;
    // origin_shift 保持为 [0;3]（默认值）

    // 4. 计算 changed_symmetry: 使用完整的合成（representatives × pure_trans × factors）
    let tmat = ref_sg.bravais_lattice;
    let shift = ref_sg.origin_shift;
    let representatives = build_representatives(msgtype_num, magnetic_symmetry)?;
    let changed_symmetry = get_changed_magnetic_symmetry(
        &tmat, &shift, &representatives, &sym_xsg,
        magnetic_symmetry, symprec,
    )?;

    // 5. 复制 ref_sg 用于返回
    let ref_sg_copy = ref_sg.clone();

    Some((ref_sg_copy, changed_symmetry, tmat, shift, msgtype_num))
}

/// 从磁性对称操作中提取普通对称操作（不搜索空间群）。
/// `ignore_time_reversal=true` → FSG (所有操作)
/// `ignore_time_reversal=false` → XSG (仅 timerev=0)
pub(crate) fn extract_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    ignore_time_reversal: bool,
    symprec: f64,
) -> Option<Symmetry> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let num_sym_msg = magnetic_symmetry.size;

    // Check if MSG is type-II
    let is_type2 = magnetic_symmetry
        .rot
        .iter()
        .zip(magnetic_symmetry.trans.iter())
        .zip(magnetic_symmetry.timerev.iter())
        .any(|((rot, trans), &timerev)| {
            mat_check_identity_matrix_i3(&identity, rot)
                && trans[0].abs() < symprec
                && trans[1].abs() < symprec
                && trans[2].abs() < symprec
                && timerev
        });

    // Collect unique ops
    let mut sym = Symmetry::new(num_sym_msg);
    let mut num_sym = 0;
    for i in 0..num_sym_msg {
        if !ignore_time_reversal && magnetic_symmetry.timerev[i] {
            continue;
        }
        if is_type2 && magnetic_symmetry.timerev[i] {
            continue;
        }
        sym.rot[num_sym] = magnetic_symmetry.rot[i];
        sym.trans[num_sym] = magnetic_symmetry.trans[i];
        num_sym += 1;
    }
    sym.size = num_sym;

    // Deduplicate FSG ops: same (R, t) with different timerev produce duplicates
    if ignore_time_reversal || is_type2 {
        let mut dedup = Symmetry::new(num_sym);
        let mut n_dedup = 0;
        for i in 0..num_sym {
            let mut dup = false;
            for j in 0..n_dedup {
                if !mat_check_identity_matrix_i3(&dedup.rot[j], &sym.rot[i]) {
                    continue;
                }
                let mut diff = 0.0;
                for k in 0..3 {
                    let x = dedup.trans[j][k] - sym.trans[i][k];
                    diff += (x - x.round()).abs();
                }
                if diff < symprec {
                    dup = true;
                    break;
                }
            }
            if !dup {
                dedup.rot[n_dedup] = sym.rot[i];
                dedup.trans[n_dedup] = sym.trans[i];
                n_dedup += 1;
            }
        }
        dedup.size = n_dedup;
        sym = dedup;
    }

    if sym.size == 0 {
        None
    } else {
        Some(sym)
    }
}

/// Get family space group (FSG) and its symmetry.
fn get_family_space_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    get_space_group_with_magnetic_symmetry(magnetic_symmetry, true, symprec)
}

/// Get maximal subspace group (XSG) with space group search.
fn get_maximal_subspace_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    get_space_group_with_magnetic_symmetry(magnetic_symmetry, false, symprec)
}

/// Get maximal subspace group symmetry (XSG) — 仅提取对称操作, 不搜索空间群.
/// 保留作为 fallback 和其他不需要空间群搜索的场景。
fn get_maximal_subspace_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<Symmetry> {
    extract_symmetry(magnetic_symmetry, false, symprec)
}

/// Get space group from magnetic symmetry.
///
/// `ignore_time_reversal=true` → FSG (family space group, all ops regardless of timerev)
/// `ignore_time_reversal=false` → XSG (maximal subspace group, only ordinary ops)
///
/// Returns (Spacegroup, Symmetry) pair.
fn get_space_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    ignore_time_reversal: bool,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let unit_lat: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let num_sym_msg = magnetic_symmetry.size;

    // Check if MSG is type-II (has pure time-reversal operation (I, 0)1')
    let is_type2 = magnetic_symmetry
        .rot
        .iter()
        .zip(magnetic_symmetry.trans.iter())
        .zip(magnetic_symmetry.timerev.iter())
        .any(|((rot, trans), &timerev)| {
            mat_check_identity_matrix_i3(&identity, rot)
                && trans[0].abs() < symprec
                && trans[1].abs() < symprec
                && trans[2].abs() < symprec
                && timerev
        });

    // Extract operations
    let mut sym = Symmetry::new(num_sym_msg);
    let mut num_sym = 0;
    for i in 0..num_sym_msg {
        if !ignore_time_reversal && magnetic_symmetry.timerev[i] {
            continue;
        }
        // For type-II MSG, deduplicate: skip timerev=1 copies
        if is_type2 && magnetic_symmetry.timerev[i] {
            continue;
        }
        sym.rot[num_sym] = magnetic_symmetry.rot[i];
        sym.trans[num_sym] = magnetic_symmetry.trans[i];
        num_sym += 1;
    }
    sym.size = num_sym;

    if num_sym == 0 {
        return None;
    }

    // Deduplicate ops when ignoring time reversal: same (R, t) across timerev=0/1
    // duplicates would confuse the space group search (e.g. BCC AFM: 24 ops with
    // 12 ordinary + 12 anti → FSG should be 12 unique ops, not 24).
    if ignore_time_reversal || is_type2 {
        let mut dedup = Symmetry::new(num_sym);
        let mut n_dedup = 0;
        for i in 0..num_sym {
            let mut dup = false;
            for j in 0..n_dedup {
                if !mat_check_identity_matrix_i3(&dedup.rot[j], &sym.rot[i]) {
                    continue;
                }
                let mut diff = 0.0;
                for k in 0..3 {
                    let x = dedup.trans[j][k] - sym.trans[i][k];
                    diff += (x - x.round()).abs();
                }
                if diff < symprec {
                    dup = true;
                    break;
                }
            }
            if !dup {
                dedup.rot[n_dedup] = sym.rot[i];
                dedup.trans[n_dedup] = sym.trans[i];
                n_dedup += 1;
            }
        }
        dedup.size = n_dedup;
        sym = dedup;
    }

    // Get primitive symmetry: (a, b, c) = (a_prim, b_prim, c_prim) @ tmat
    let (tmat, prim_sym) = prm_get_primitive_symmetry(&sym, symprec)?;

    let mut spacegroup = match spa_search_spacegroup_with_symmetry(&prim_sym, &unit_lat, symprec) {
        Some(sg) => sg,
        None => {
            // 标准空间群搜索失败 → 使用 fallback
            return find_spacegroup_by_symmetry(&sym, &unit_lat, symprec)
                .map(|sg| (sg, sym));
        }
    };

    // Refine bravais lattice and origin_shift
    ref_find_similar_bravais_lattice(&mut spacegroup, symprec);

    // Change basis from primitive to original:
    // x = (tmat, 0)^-1 x_prim
    // => x_std = (P^-1, p) (tmat, 0) x = ( P^-1 @ tmat, p) x
    //    (a_std, b_std, c_std) = (a, b, c) @ tmat^-1 @ P
    let inv_tmat = mat_inverse_matrix_d3(&tmat, 0.0)?;
    spacegroup.bravais_lattice = mat_multiply_matrix_d3(&inv_tmat, &spacegroup.bravais_lattice);

    Some((spacegroup, sym))
}

/// Fallback: 直接从对称操作匹配 Hall 编号，绕过完整的空间群搜索。
/// 当 `spa_search_spacegroup_with_symmetry` 失败时使用。
fn find_spacegroup_by_symmetry(
    symmetry: &Symmetry,
    lattice: &Mat3,
    symprec: f64,
) -> Option<Spacegroup> {
    let mut origin_shift = [0.0; 3];
    let mut conv_lattice = [[0.0; 3]; 3];

    let (tmat_int, pointgroup) = ptg_get_transformation_matrix(&symmetry.rot, None);
    if pointgroup.number == 0 {
        return None;
    }

    let mut correction_mat = [[0.0; 3]; 3];
    let centering = get_centering(&mut correction_mat, &tmat_int, pointgroup.laue);
    if centering == Centering::Error {
        return None;
    }

    let tmat = mat_multiply_matrix_id3(&tmat_int, &correction_mat);
    conv_lattice = mat_multiply_matrix_d3(lattice, &tmat);

    let conv_symmetry = get_initial_conventional_symmetry(centering, &tmat, symmetry)?;

    // Try ALL 530 Hall numbers (not just the 230 representatives)
    for hall in 1..=530 {
        if hal_match_hall_symbol_db(
            &mut origin_shift,
            &conv_lattice,
            hall,
            centering,
            &conv_symmetry,
            symprec,
        ) {
            let spg_type = spgdb_get_spacegroup_type(hall as usize);
            let mut spacegroup = Spacegroup::new();
            spacegroup.bravais_lattice = conv_lattice;
            spacegroup.origin_shift = origin_shift;
            spacegroup.number = spg_type.number;
            spacegroup.hall_number = hall as usize;
            spacegroup.pointgroup_number = spg_type.pointgroup_number;
            spacegroup.schoenflies = spg_type.schoenflies;
            spacegroup.hall_symbol = spg_type.hall_symbol;
            spacegroup.international = spg_type.international;
            spacegroup.international_long = spg_type.international_full;
            spacegroup.international_short = spg_type.international_short;
            spacegroup.choice = spg_type.choice;
            return Some(spacegroup);
        }
    }

    None
}

/// Build coset representatives for the MSG type.
///
/// Type-1: identity only
/// Type-2: identity + identity with time reversal
/// Type-3/4: identity + one primed operation (via `get_representative`)
fn build_representatives(
    msgtype: MagneticType,
    magnetic_symmetry: &MagneticSymmetry,
) -> Option<MagneticSymmetry> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    match msgtype {
        MagneticType::Ordinary => {
            let mut rep = MagneticSymmetry::new(1);
            rep.rot[0] = identity;
            rep.trans[0] = [0.0; 3];
            rep.timerev[0] = false;
            rep.size = 1;
            Some(rep)
        }
        MagneticType::Grey => {
            let mut rep = MagneticSymmetry::new(2);
            rep.rot[0] = identity;
            rep.trans[0] = [0.0; 3];
            rep.timerev[0] = false;
            rep.rot[1] = identity;
            rep.trans[1] = [0.0; 3];
            rep.timerev[1] = true;
            rep.size = 2;
            Some(rep)
        }
        MagneticType::BlackWhite | MagneticType::AntiTranslation => {
            get_representative(magnetic_symmetry)
        }
        _ => None,
    }
}

/// Determine MSG type. Returns None if failed.
fn get_magnetic_space_group_type(
    magnetic_symmetry: &MagneticSymmetry,
    num_sym_fsg: usize,
    num_sym_xsg: usize,
    symprec: f64,
) -> Option<MagneticType> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    if num_sym_fsg == num_sym_xsg {
        let num_sym_msg = magnetic_symmetry.size;
        if num_sym_msg == num_sym_fsg {
            // Type-I: all operations are ordinary
            Some(MagneticType::Ordinary)
        } else if num_sym_msg == 2 * num_sym_fsg {
            // Type-II: has pure time-reversal operation
            Some(MagneticType::Grey)
        } else {
            None
        }
    } else if num_sym_fsg == 2 * num_sym_xsg {
        // Check if anti-operation is translation (type-IV) or rotation (type-III)
        // Note: in the reduced (primitive) representation, a genuine anti-translation
        // (I|t, timerev=true) with |t| ≈ 0 is a coset representative.  A centering+TR
        // combination (I|t, timerev=true) with |t| > 0 must be reduced to the primitive
        // cell to check — otherwise BCC AFM centering ops would be misclassified as
        // Type-IV.  The caller is responsible for passing already-reduced ops.
        let has_anti_translation = magnetic_symmetry
            .rot
            .iter()
            .zip(magnetic_symmetry.trans.iter())
            .zip(magnetic_symmetry.timerev.iter())
            .any(|((rot, trans), &timerev)| {
                mat_check_identity_matrix_i3(&identity, rot)
                    && trans[0].abs() < symprec
                    && trans[1].abs() < symprec
                    && trans[2].abs() < symprec
                    && timerev
            });

        if has_anti_translation {
            Some(MagneticType::AntiTranslation) // Type-IV: anti-translation
        } else {
            Some(MagneticType::BlackWhite) // Type-III: anti-rotation
        }
    } else {
        None
    }
}

/// Get coset representative of XSG in MSG.
/// For type-III and type-IV MSGs. Returns identity for type-I/II.
fn get_representative(magnetic_symmetry: &MagneticSymmetry) -> Option<MagneticSymmetry> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    let mut representative = MagneticSymmetry::new(2);

    // Set the first representative as identity
    representative.rot[0] = identity;
    representative.trans[0] = [0.0; 3];
    representative.timerev[0] = false;

    // Type-IV: anti-translation with identity rotation AND zero translation.
    // Centering translations with timerev=1 (e.g. in AFM structures) must NOT be
    // captured here — they belong to pure translations, not to the coset representative.
    for i in 0..magnetic_symmetry.size {
        if mat_check_identity_matrix_i3(&identity, &magnetic_symmetry.rot[i])
            && magnetic_symmetry.timerev[i]
            && magnetic_symmetry.trans[i][0].abs() < 1e-5
            && magnetic_symmetry.trans[i][1].abs() < 1e-5
            && magnetic_symmetry.trans[i][2].abs() < 1e-5
        {
            representative.rot[1] = magnetic_symmetry.rot[i];
            representative.trans[1] = magnetic_symmetry.trans[i];
            representative.timerev[1] = true;
            representative.size = 2;
            return Some(representative);
        }
    }

    // Type-III: anti-rotation.  Skip identity rotations — those with timerev=1
    // and non-zero translation are centering translations, not coset representatives.
    for i in 0..magnetic_symmetry.size {
        if !magnetic_symmetry.timerev[i] {
            continue;
        }
        if mat_check_identity_matrix_i3(&identity, &magnetic_symmetry.rot[i]) {
            continue; // skip identity rotations (centering translations with timerev)
        }
        representative.rot[1] = magnetic_symmetry.rot[i];
        representative.trans[1] = magnetic_symmetry.trans[i];
        representative.timerev[1] = true;
        representative.size = 2;
        return Some(representative);
    }

    // Type-I or II: no anti-operations
    representative.size = 1;
    Some(representative)
}

/// Apply transformation (tmat, shift) to magnetic symmetry, deduplicating.
///
/// R' = T^-1 * R * T
/// t' = T^-1 * (t + (R - I) * shift)
fn get_distinct_changed_magnetic_symmetry(
    tmat: &Mat3,
    shift: &Vec3,
    sym_msg: &MagneticSymmetry,
) -> Option<MagneticSymmetry> {
    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0)?;

    let mut changed = MagneticSymmetry::new(sym_msg.size);
    let mut size = 0usize;

    for i in 0..sym_msg.size {
        // R' = T^-1 * R * T
        let rot_f64 = mat_cast_matrix_3i_to_3d(&sym_msg.rot[i]);
        let tmp = mat_multiply_matrix_d3(&inv_tmat, &rot_f64);
        let r_new = mat_multiply_matrix_d3(&tmp, tmat);

        // Round to integer rotation matrix
        let rot_i = [
            [mat_nint(r_new[0][0]), mat_nint(r_new[0][1]), mat_nint(r_new[0][2])],
            [mat_nint(r_new[1][0]), mat_nint(r_new[1][1]), mat_nint(r_new[1][2])],
            [mat_nint(r_new[2][0]), mat_nint(r_new[2][1]), mat_nint(r_new[2][2])],
        ];

        // t' = T^-1 * (t + (R - I) * shift)
        let mut tmp_rot = sym_msg.rot[i];
        tmp_rot[0][0] -= 1;
        tmp_rot[1][1] -= 1;
        tmp_rot[2][2] -= 1;
        let shift_term = mat_multiply_matrix_vector_id3(&tmp_rot, shift);

        let mut t_new = [0.0; 3];
        for j in 0..3 {
            t_new[j] = sym_msg.trans[i][j] + shift_term[j];
        }
        t_new = mat_multiply_matrix_vector_d3(&inv_tmat, &t_new);

        // Check for uniqueness (same rotation, same translation, same timerev)
        let mut is_dup = false;
        for j in 0..size {
            if mat_check_identity_matrix_i3(&changed.rot[j], &rot_i) {
                let mut diff = [0.0; 3];
                for k in 0..3 {
                    diff[k] = changed.trans[j][k] - t_new[k];
                    diff[k] -= mat_nint(diff[k]) as f64;
                }
                if diff[0].abs() < 1e-5 && diff[1].abs() < 1e-5 && diff[2].abs() < 1e-5 {
                    if changed.timerev[j] == sym_msg.timerev[i] {
                        is_dup = true;
                        break;
                    }
                }
            }
        }

        if !is_dup {
            changed.rot[size] = rot_i;
            changed.trans[size] = t_new;
            changed.timerev[size] = sym_msg.timerev[i];
            size += 1;
        }
    }

    changed.size = size;
    Some(changed)
}

/// `a - round(a)`, 结果在 [0, 1)。
fn mat_dmod1(a: f64) -> f64 {
    let b = a - a.round();
    if b < 0.0 { b + 1.0 } else { b }
}

/// 检查旋转矩阵 `a` 是否已包含在 `sym_msg.rot[0..size]` 中。
fn is_contained_mat(a: &Mat3I, sym_msg: &MagneticSymmetry, size: usize) -> bool {
    for i in 0..size {
        if mat_check_identity_matrix_i3(a, &sym_msg.rot[i]) {
            return true;
        }
    }
    false
}

/// 检查向量 `v` 是否已包含在 `trans` 中 (tol = symprec)。
fn is_contained_vec(v: &Vec3, trans: &[Vec3], symprec: f64) -> bool {
    for t in trans {
        let mut eq = true;
        for s in 0..3 {
            if (v[s] - t[s]).abs() >= symprec {
                eq = false;
                break;
            }
        }
        if eq { return true; }
    }
    false
}

/// (I, w) = (tmat, shift)^-1 (I, w_std) (tmat, shift) — 纯平移变换。
/// 从输入晶胞的纯平移到参考设置的纯平移。
fn get_changed_pure_translations(
    tmat: &Mat3,
    pure_trans: &[Vec3],
    symprec: f64,
) -> Option<Vec<Vec3>> {
    let det = mat_get_determinant_d3(tmat);
    let size = mat_nint(pure_trans.len() as f64 / det.abs()) as usize;

    let mut changed: Vec<Vec3> = Vec::with_capacity(size);

    if (det - 1.0).abs() <= symprec {
        for pt in pure_trans {
            let trans = mat_multiply_matrix_vector_d3(tmat, pt);
            changed.push([mat_dmod1(trans[0]), mat_dmod1(trans[1]), mat_dmod1(trans[2])]);
        }
    } else {
        // 查找转动矩阵元素的最小公分母
        let mut denominator = 1;
        loop {
            let mut ok = true;
            for s in 0..3 {
                for t in 0..3 {
                    if (tmat[s][t] * denominator as f64
                        - mat_nint(tmat[s][t] * denominator as f64) as f64)
                        .abs()
                        > symprec
                    {
                        ok = false;
                        break;
                    }
                }
                if !ok { break; }
            }
            if ok { break; }
            denominator += 1;
            if denominator as f64 > MAX_DENOMINATOR {
                return None;
            }
        }


        // 为每个纯平移尝试额外的晶格矢量以恢复常规晶胞中的平移
        for n0 in 0..=denominator {
            for n1 in 0..=denominator {
                for n2 in 0..=denominator {
                    for pt in pure_trans {
                        let shifted = [
                            pt[0] + n0 as f64,
                            pt[1] + n1 as f64,
                            pt[2] + n2 as f64,
                        ];
                        let trans = mat_multiply_matrix_vector_d3(tmat, &shifted);
                        let t_mod = [mat_dmod1(trans[0]), mat_dmod1(trans[1]), mat_dmod1(trans[2])];

                        if !is_contained_vec(&t_mod, &changed, symprec) {
                            changed.push(t_mod);
                        }
                    }
                }
            }
        }
    }

    if changed.len() != size {
        // Fallback: just transform each pure_trans directly, mod1, and dedup
        changed.clear();
        for pt in pure_trans {
            let trans = mat_multiply_matrix_vector_d3(tmat, pt);
            let t_mod = [mat_dmod1(trans[0]), mat_dmod1(trans[1]), mat_dmod1(trans[2])];
            if !is_contained_vec(&t_mod, &changed, symprec) {
                changed.push(t_mod);
            }
        }
        // Accept whatever we got (may be more or less than expected due to basis change)
    }

    Some(changed)
}

/// 合成完整的变换磁性对称操作: representatives × pure_trans × factor_group。
///
/// 这是 C 原版 `get_changed_magnetic_symmetry` 的直接移植。
/// 与 `get_distinct_changed_magnetic_symmetry` 不同，本函数不只是简单地对每个操作
/// 做基变换，而是将其分解为代表元（representatives）、纯平移（pure_translations）
/// 和因子群（factor_group），在参考设置下重新合成，以匹配数据库的标准化表示。
fn get_changed_magnetic_symmetry(
    tmat: &Mat3,
    shift: &Vec3,
    representatives: &MagneticSymmetry,
    sym_xsg: &Symmetry,
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<MagneticSymmetry> {
    // 1. 代表元在参考设置下的形式
    let changed_representatives =
        match get_distinct_changed_magnetic_symmetry(tmat, shift, representatives) {
            Some(r) => r,
            None => return None,
        };

    // 2. 收集原始磁性对称中的纯平移（仅 timerev=0），变换到参考设置
    let pure_trans = match crate::spin::spn_collect_pure_translations_from_magnetic_symmetry(
        magnetic_symmetry,
    ) {
        Some(p) => {
            p
        }
        None => return None,
    };
    let changed_pure_trans = match get_changed_pure_translations(tmat, &pure_trans, symprec) {
        Some(p) => p,
        None => return None,
    };

    // 3. 从 XSG 对称性中收集因子群（仅去重旋转部分，timerev=0）
    let mut factors = MagneticSymmetry::new(sym_xsg.size);
    let mut num_factors = 0;
    for i in 0..sym_xsg.size {
        if is_contained_mat(&sym_xsg.rot[i], &factors, num_factors) {
            continue;
        }
        factors.rot[num_factors] = sym_xsg.rot[i];
        factors.trans[num_factors] = sym_xsg.trans[i];
        factors.timerev[num_factors] = false;
        num_factors += 1;
    }
    factors.size = num_factors;
    let changed_factors =
        match get_distinct_changed_magnetic_symmetry(tmat, shift, &factors) {
            Some(f) => f,
            None => return None,
        };

    // 4. 合成: (I, ti)(Pj, tj)(Pk, tk) = (Pj * Pk, Pj * tk + tj + ti)
    let size = changed_representatives.size * changed_pure_trans.len() * num_factors;
    let mut changed = MagneticSymmetry::new(size);
    let mut num_sym = 0;

    for i in 0..changed_pure_trans.len() {
        for j in 0..changed_representatives.size {
            for k in 0..num_factors {
                // R = Pj * Pk
                changed.rot[num_sym] = mat_multiply_matrix_i3(
                    &changed_representatives.rot[j],
                    &changed_factors.rot[k],
                );

                // t = Pj * tk + tj + ti
                let mut trans = mat_multiply_matrix_vector_id3(
                    &changed_representatives.rot[j],
                    &changed_factors.trans[k],
                );
                for s in 0..3 {
                    trans[s] += changed_representatives.trans[j][s]
                        + changed_pure_trans[i][s];
                    trans[s] = mat_dmod1(trans[s]);
                }
                changed.trans[num_sym] = trans;

                // timerev = changed_representatives.timerev XOR changed_factors.timerev
                // (factors 都是 ordinary，所以 XOR 就是 representatives 的 timerev)
                changed.timerev[num_sym] =
                    changed_representatives.timerev[j] != changed_factors.timerev[k];

                num_sym += 1;
            }
        }
    }

    changed.size = num_sym;
    Some(changed)
}

/// 检查两个磁性对称操作集合是否等价。
/// 用纯平移将磁对称操作约化到原胞表示。
/// 收集 identity 旋转的操作为晶格平移（非零 timerev=1 的平移也用于约化，
/// 但零平移 timerev=1 除外，避免破坏 Type-2 结构）。
pub(crate) fn reduce_to_primitive_magsym(sym: &MagneticSymmetry, symprec: f64) -> MagneticSymmetry {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let n_in = sym.size;

    // 1. 收集纯平移（identity 旋转）。
    //    timerev=1 的零平移 identity 是 Type-2 结构的一部分，不能用作约化平移。
    //    只有 timerev=1 且非零平移的（如 BCC AFM 中的 (0.5,0.5,0.5)）才用于约化。
    let mut lat_trans: Vec<(Vec3, bool)> = vec![([0.0; 3], false)];
    for i in 0..n_in {
        if mat_check_identity_matrix_i3(&identity, &sym.rot[i]) {
            let t = sym.trans[i];
            let tr = sym.timerev[i];
            if tr && t[0].abs() < symprec && t[1].abs() < symprec && t[2].abs() < symprec {
                continue;
            }
            let mut dup = false;
            for (lt, ltr) in &lat_trans {
                let mut d = 0.0;
                for k in 0..3 {
                    let mut x = t[k] - lt[k];
                    x -= x.round();
                    d += x.abs();
                }
                if d < symprec && *ltr == tr { dup = true; break; }
            }
            if !dup { lat_trans.push((t, tr)); }
        }
    }

    // 如果没有多余平移, 返回原样
    if lat_trans.len() <= 1 { return sym.clone(); }

    // 2. 约化和去重: 对每个操作, 尝试所有 lat_trans (±), 选最短 translation
    let mut out = MagneticSymmetry::new(n_in);
    let mut n = 0usize;

    'next: for i in 0..n_in {
        let rot_i = &sym.rot[i];
        let t_i = sym.trans[i];
        let tr_i = sym.timerev[i];

        let mut best_t = t_i;
        let mut best_tr = tr_i;
        let mut best_n2 = t_i[0] * t_i[0] + t_i[1] * t_i[1] + t_i[2] * t_i[2];
        for (lt, ltr) in &lat_trans {
            // 减去晶格平移: (R, t, tr) → (R, t - lt, tr XOR ltr)
            let cand = [t_i[0] - lt[0], t_i[1] - lt[1], t_i[2] - lt[2]];
            let n2 = cand[0] * cand[0] + cand[1] * cand[1] + cand[2] * cand[2];
            let new_tr = tr_i != *ltr;
            if n2 < best_n2 - 1e-10 || (n2 < best_n2 + 1e-10 && !new_tr && best_tr) {
                best_t = cand;
                best_n2 = n2;
                best_tr = new_tr;
            }
        }

        // 去重
        for j in 0..n {
            if out.timerev[j] != best_tr { continue; }
            if !mat_check_identity_matrix_i3(&out.rot[j], rot_i) { continue; }
            let mut d = 0.0;
            for k in 0..3 {
                let mut x = out.trans[j][k] - best_t[k];
                x -= x.round();
                d += x.abs();
            }
            if d < symprec { continue 'next; }
        }

        out.rot[n] = *rot_i;
        out.trans[n] = best_t;
        out.timerev[n] = best_tr;
        n += 1;
    }
    out.size = n;

    // 3. 收缩到实际大小
    let mut final_sym = MagneticSymmetry::new(n);
    for i in 0..n {
        final_sym.rot[i] = out.rot[i];
        final_sym.trans[i] = out.trans[i];
        final_sym.timerev[i] = out.timerev[i];
    }
    final_sym
}

/// 子集检查: sym1 的所有操作是否都能在 sym2 中找到。
/// sym1 可以是 sym2 的子集（size 不要求相等）。
fn is_subset(
    sym1: &MagneticSymmetry,
    sym2: &MagneticSymmetry,
    symprec: f64,
) -> bool {
    if sym1.size > sym2.size {
        return false; // sym1 不可能是 sym2 的子集
    }

    let mut found = vec![false; sym2.size];
    for i in 0..sym1.size {
        let mut matched = false;
        for j in 0..sym2.size {
            if found[j] {
                continue;
            }
            if !mat_check_identity_matrix_i3(&sym1.rot[i], &sym2.rot[j]) {
                continue;
            }
            if sym1.timerev[i] != sym2.timerev[j] {
                continue;
            }
            let mut diff = [0.0; 3];
            for k in 0..3 {
                diff[k] = sym1.trans[i][k] - sym2.trans[j][k];
                diff[k] -= mat_nint(diff[k]) as f64;
            }
            if diff[0].abs() < symprec && diff[1].abs() < symprec && diff[2].abs() < symprec {
                found[j] = true;
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }
    true
}

#[allow(dead_code)]
fn is_equal(
    sym1: &MagneticSymmetry,
    sym2: &MagneticSymmetry,
    symprec: f64,
) -> bool {
    if sym1.size != sym2.size {
        return false;
    }

    let mut found = vec![false; sym2.size];
    for i in 0..sym1.size {
        let mut matched = false;
        for j in 0..sym2.size {
            if found[j] {
                continue;
            }
            if !mat_check_identity_matrix_i3(&sym1.rot[i], &sym2.rot[j]) {
                continue;
            }
            if sym1.timerev[i] != sym2.timerev[j] {
                continue;
            }
            let mut diff = [0.0; 3];
            for k in 0..3 {
                diff[k] = sym1.trans[i][k] - sym2.trans[j][k];
                diff[k] -= mat_nint(diff[k]) as f64;
            }
            if diff[0].abs() < symprec && diff[1].abs() < symprec && diff[2].abs() < symprec {
                found[j] = true;
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }
    true
}

/// 计算刚性旋转矩阵。
fn get_rigid_rotation(
    rigid_rot: &mut Mat3,
    lattice: &Mat3,
    tmat: &Mat3,
    ref_sg: &Spacegroup,
) {
    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0);
    if let Some(inv) = inv_tmat {
        let tmp = mat_multiply_matrix_d3(&ref_sg.bravais_lattice, &inv);
        let inv_lat = mat_inverse_matrix_d3(lattice, 0.0);
        if let Some(inv_l) = inv_lat {
            let result = mat_multiply_matrix_d3(&inv_l, &tmp);
            *rigid_rot = result;
        }
    }
}

