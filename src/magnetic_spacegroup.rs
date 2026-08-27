//! 磁性空间群识别。
//!
//! 使用磁性对称操作数据库识别磁性空间群类型。
//! 参考: Litvin, "Magnetic Group Tables" 2013

use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

use crate::MagneticType;
use crate::SymError;
use crate::hall_symbol::hal_match_hall_symbol_db;
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_cast_matrix_3i_to_3d, mat_check_identity_matrix_i3, mat_dmod1,
    mat_get_determinant_d3, mat_inverse_matrix_d3, mat_multiply_matrix_d3, mat_multiply_matrix_i3,
    mat_multiply_matrix_id3, mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_id3,
    mat_nint,
};
use crate::msg_database::{
    get_magnetic_spacegroup_type, get_spacegroup_operations, get_std_transformations,
    get_uni_candidates,
};
use crate::pointgroup::get_transformation_matrix;
use crate::primitive::get_primitive_symmetry;
use crate::refinement::find_similar_bravais_lattice;
use crate::spacegroup::{
    Spacegroup, get_centering, get_initial_conventional_symmetry, search_spacegroup_with_symmetry,
};
use crate::spg_database::{Centering, get_spacegroup_type};
use crate::symmetry::{MagneticSymmetry, Symmetry};

const MAX_DENOMINATOR: f64 = 100.0;
const MAX_CHANGED_PURE_TRANSLATIONS: usize = 1_000_000;
const DATABASE_TRANSLATION_DENOMINATOR: i32 = 12;
const DATABASE_CANONICAL_SYMPREC: f64 = 1e-5;
const UNIT_LATTICE: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// 磁性空间群识别结果的中间数据结构。
pub struct MagneticDataset {
    pub uni_number: usize,
    pub msg_type: MagneticType,
    pub hall_number: usize,
    pub transformation_matrix: Mat3,
    pub origin_shift: Vec3,
    pub std_rotation_matrix: Mat3,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct DatabaseMagneticOperationKey {
    rotation: [i32; 9],
    translation_twelfths: [i32; 3],
    time_reversal: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct CanonicalMagneticKey {
    reference_hall_number: usize,
    operations: Vec<DatabaseMagneticOperationKey>,
}

#[derive(Clone)]
struct CanonicalizedMagneticSymmetry {
    key: CanonicalMagneticKey,
    reference_spacegroup: Spacegroup,
    transformation_matrix: Mat3,
    origin_shift: Vec3,
    msg_type: MagneticType,
}

#[derive(Clone)]
struct CanonicalDatabaseSetting {
    uni_number: usize,
    hall_number: usize,
    transformation_matrix: Mat3,
    origin_shift: Vec3,
}

#[derive(Clone)]
struct CanonicalDatabaseMatch {
    input: CanonicalizedMagneticSymmetry,
    candidates: Vec<CanonicalDatabaseSetting>,
}

static TYPE_IV_CANONICAL_INDEX: OnceLock<
    HashMap<CanonicalMagneticKey, Vec<CanonicalDatabaseSetting>>,
> = OnceLock::new();

/// 识别磁性空间群类型。
///
/// 给定晶格和磁性对称操作，返回识别出的磁性数据集。
/// 对只靠磁操作无法区分的 Type-IV BNS parent，返回
/// [`SymError::MagneticUniAmbiguous`]，而不是静默选择某个 UNI。
pub fn identify_magnetic_space_group_type(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Result<MagneticDataset, SymError> {
    identify_with_parent_hall(lattice, magnetic_symmetry, None, symprec)
}

/// 与 [`identify_magnetic_space_group_type`] 相同，但可指定非磁母空间群的
/// Hall 编号。规范数据库 setting 的严格匹配优先于自动标准化；对一般基变换和
/// 原点移动后的输入，则用母群空间群号筛选完整的 Type-IV 规范等价类。
pub fn identify_with_parent_hall(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: Option<usize>,
    symprec: f64,
) -> Result<MagneticDataset, SymError> {
    if let Some(hall_number) = parent_hall_number
        && let Some(dataset) = match_exact_parent_setting(magnetic_symmetry, hall_number, symprec)
    {
        return Ok(dataset);
    }

    // Type-IV standardization is not injective: distinct BNS parent groups
    // can have the same standardized XSG representation.  Build the complete
    // equivalence class from the database before the legacy single-Hall
    // search discards that information.  A supplied parent Hall filters the
    // class by its non-magnetic parent space-group number; without that hint,
    // a cross-UNI class is a genuine ambiguity and must not be guessed.
    let canonical_match = match_type_iv_canonical_class(magnetic_symmetry, symprec);
    let mut canonical_fallback = None;
    if let Some(canonical_match) = canonical_match {
        let mut candidates = canonical_match.candidates.clone();
        if let Some(parent_hall) = parent_hall_number {
            let parent_number = get_spacegroup_type(parent_hall).number;
            let parent_candidates: Vec<_> = candidates
                .iter()
                .filter(|candidate| {
                    get_magnetic_spacegroup_type(candidate.uni_number).number == parent_number
                })
                .cloned()
                .collect();
            if !parent_candidates.is_empty() {
                candidates = parent_candidates;
            }
        }

        let distinct_unis: HashSet<_> = candidates
            .iter()
            .map(|candidate| candidate.uni_number)
            .collect();
        if distinct_unis.len() > 1 {
            return Err(SymError::MagneticUniAmbiguous);
        }

        if let Some(candidate) =
            choose_canonical_candidate(&candidates, magnetic_symmetry, parent_hall_number, symprec)
        {
            canonical_fallback = dataset_from_canonical_candidate(
                lattice,
                magnetic_symmetry,
                &canonical_match.input,
                &candidate,
                symprec,
            );
            if parent_hall_number.is_some()
                && let Some(dataset) = canonical_fallback
            {
                return Ok(dataset);
            }
        }
    }

    let standardized = identify_in_single_reference_setting(
        lattice,
        magnetic_symmetry,
        parent_hall_number,
        symprec,
    );

    match (standardized, canonical_fallback) {
        (Ok(dataset), Some(canonical)) if dataset.uni_number != canonical.uni_number => {
            Ok(canonical)
        }
        (Ok(dataset), _) => Ok(dataset),
        (Err(_), Some(canonical)) => Ok(canonical),
        (Err(error), None) => Err(error),
    }
}

/// Original spglib-compatible single-reference-setting identification path.
fn identify_in_single_reference_setting(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: Option<usize>,
    symprec: f64,
) -> Result<MagneticDataset, SymError> {
    // 标准路径: 从磁对称性中提取 FSG/XSG 并搜索空间群
    let (ref_sg, changed_symmetry, mut tmat, mut shift, msgtype_num) =
        match get_reference_space_group(lattice, magnetic_symmetry, symprec) {
            Some(result) => result,
            None => {
                // 标准路径失败 → 尝试 fallback（用母空间群的 Hall 编号）
                let hall_number =
                    parent_hall_number.ok_or(SymError::MagneticReferenceGroupFailed)?;
                build_fallback_reference(lattice, magnetic_symmetry, hall_number, symprec)
                    .ok_or(SymError::MagneticFallbackReferenceFailed)?
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
        let range = match get_uni_candidates(hall_number) {
            Some(r) => r,
            None => continue,
        };
        let min_uni = range[0];
        let max_uni = range[1];

        for uni_number in min_uni..=max_uni {
            let msgtype_db = get_magnetic_spacegroup_type(uni_number);
            if msgtype_db.type_ != msgtype_num {
                continue;
            }
            let msg_uni = match get_spacegroup_operations(uni_number, hall_number) {
                Some(u) => u,
                None => continue,
            };
            if changed_symmetry.len() != msg_uni.len() {
                continue;
            }

            let transformations = match get_std_transformations(uni_number, hall_number) {
                Some(t) => t,
                None => continue,
            };

            let mut same = false;
            for trans_idx in 0..transformations.len() {
                let tmat_cor = transformations.rot[trans_idx];
                let shift_cor = transformations.trans[trans_idx];

                let tmat_cor_d = mat_cast_matrix_3i_to_3d(&tmat_cor);

                let symmetry_cor = get_distinct_changed_magnetic_symmetry(
                    &tmat_cor_d,
                    &shift_cor,
                    &changed_symmetry,
                );

                let symmetry_cor = match symmetry_cor {
                    Some(s) => s,
                    None => continue,
                };

                let matched = is_equal(&symmetry_cor, &msg_uni, symprec);
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
        return Err(SymError::MagneticUniMatchFailed);
    }

    let hall_number = best_hall_number;

    let _msgtype = get_magnetic_spacegroup_type(best_uni);
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

    Ok(ret)
}

/// Match an already-canonical magnetic operation set within an explicitly
/// supplied family-space-group setting.
///
/// This exact fast path preserves the caller's parent setting before Type-IV
/// XSG standardization. The general changed-basis path below performs the same
/// disambiguation through the database-derived canonical equivalence classes.
fn match_exact_parent_setting(
    magnetic_symmetry: &MagneticSymmetry,
    hall_number: usize,
    symprec: f64,
) -> Option<MagneticDataset> {
    let [min_uni, max_uni] = get_uni_candidates(hall_number)?;
    for uni_number in min_uni..=max_uni {
        let Some(database_symmetry) = get_spacegroup_operations(uni_number, hall_number) else {
            continue;
        };
        if !is_equal(magnetic_symmetry, &database_symmetry, symprec) {
            continue;
        }

        let msg_type = get_magnetic_spacegroup_type(uni_number).type_;
        return Some(MagneticDataset {
            uni_number,
            msg_type,
            hall_number,
            transformation_matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            origin_shift: [0.0; 3],
            std_rotation_matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        });
    }
    None
}

fn flatten_rotation(rotation: &Mat3I) -> [i32; 9] {
    [
        rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2],
    ]
}

fn quantize_database_translation(value: f64, symprec: f64) -> Option<i32> {
    if !value.is_finite() {
        return None;
    }

    let normalized = mat_dmod1(value);
    let scaled = normalized * DATABASE_TRANSLATION_DENOMINATOR as f64;
    let rounded = scaled.round();
    let nearest = rounded / DATABASE_TRANSLATION_DENOMINATOR as f64;
    let mut difference = normalized - nearest;
    difference -= difference.round();
    if difference.abs() >= symprec.max(1e-8) {
        return None;
    }

    Some((rounded as i32).rem_euclid(DATABASE_TRANSLATION_DENOMINATOR))
}

fn canonical_magnetic_key(
    reference_hall_number: usize,
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<CanonicalMagneticKey> {
    let mut operations = Vec::with_capacity(magnetic_symmetry.len());
    for operation in 0..magnetic_symmetry.len() {
        let mut translation_twelfths = [0; 3];
        for (axis, value) in magnetic_symmetry.trans[operation]
            .iter()
            .copied()
            .enumerate()
        {
            translation_twelfths[axis] = quantize_database_translation(value, symprec)?;
        }
        operations.push(DatabaseMagneticOperationKey {
            rotation: flatten_rotation(&magnetic_symmetry.rot[operation]),
            translation_twelfths,
            time_reversal: magnetic_symmetry.timerev[operation],
        });
    }
    operations.sort();

    Some(CanonicalMagneticKey {
        reference_hall_number,
        operations,
    })
}

fn canonicalize_magnetic_symmetry_for_database(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<CanonicalizedMagneticSymmetry> {
    let (reference_spacegroup, changed_symmetry, transformation_matrix, origin_shift, msg_type) =
        get_reference_space_group(&UNIT_LATTICE, magnetic_symmetry, symprec)?;
    let key = canonical_magnetic_key(reference_spacegroup.hall_number, &changed_symmetry, symprec)?;

    Some(CanonicalizedMagneticSymmetry {
        key,
        reference_spacegroup,
        transformation_matrix,
        origin_shift,
        msg_type,
    })
}

fn type_iv_canonical_index() -> &'static HashMap<CanonicalMagneticKey, Vec<CanonicalDatabaseSetting>>
{
    TYPE_IV_CANONICAL_INDEX.get_or_init(|| {
        let mut index: HashMap<CanonicalMagneticKey, Vec<CanonicalDatabaseSetting>> =
            HashMap::new();

        for uni_number in 1usize..=1651 {
            let metadata = get_magnetic_spacegroup_type(uni_number);
            if metadata.type_ != MagneticType::AntiTranslation {
                continue;
            }
            let [num_halls, first_hall] =
                crate::msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni_number];
            for hall_number in first_hall as usize..(first_hall + num_halls) as usize {
                let Some(database_symmetry) = get_spacegroup_operations(uni_number, hall_number)
                else {
                    continue;
                };
                let Some(canonical) = canonicalize_magnetic_symmetry_for_database(
                    &database_symmetry,
                    DATABASE_CANONICAL_SYMPREC,
                ) else {
                    continue;
                };
                index
                    .entry(canonical.key)
                    .or_default()
                    .push(CanonicalDatabaseSetting {
                        uni_number,
                        hall_number,
                        transformation_matrix: canonical.transformation_matrix,
                        origin_shift: canonical.origin_shift,
                    });
            }
        }

        for candidates in index.values_mut() {
            candidates.sort_by_key(|candidate| (candidate.uni_number, candidate.hall_number));
        }
        index
    })
}

fn match_type_iv_canonical_class(
    magnetic_symmetry: &MagneticSymmetry,
    symprec: f64,
) -> Option<CanonicalDatabaseMatch> {
    let input = canonicalize_magnetic_symmetry_for_database(magnetic_symmetry, symprec)?;
    if input.msg_type != MagneticType::AntiTranslation {
        return None;
    }
    let candidates = type_iv_canonical_index().get(&input.key)?.clone();
    Some(CanonicalDatabaseMatch { input, candidates })
}

fn choose_canonical_candidate(
    candidates: &[CanonicalDatabaseSetting],
    magnetic_symmetry: &MagneticSymmetry,
    parent_hall_number: Option<usize>,
    symprec: f64,
) -> Option<CanonicalDatabaseSetting> {
    if let Some(parent_hall) = parent_hall_number
        && let Some(candidate) = candidates
            .iter()
            .find(|candidate| candidate.hall_number == parent_hall)
    {
        return Some(candidate.clone());
    }

    for candidate in candidates {
        let Some(database_symmetry) =
            get_spacegroup_operations(candidate.uni_number, candidate.hall_number)
        else {
            continue;
        };
        if is_equal(magnetic_symmetry, &database_symmetry, symprec) {
            return Some(candidate.clone());
        }
    }

    candidates.first().cloned()
}

fn dataset_from_canonical_candidate(
    lattice: &Mat3,
    magnetic_symmetry: &MagneticSymmetry,
    input: &CanonicalizedMagneticSymmetry,
    candidate: &CanonicalDatabaseSetting,
    symprec: f64,
) -> Option<MagneticDataset> {
    let inverse_candidate = mat_inverse_matrix_d3(&candidate.transformation_matrix, 0.0).ok()?;
    let transformation_matrix =
        mat_multiply_matrix_d3(&inverse_candidate, &input.transformation_matrix);
    let shift_difference = [
        input.origin_shift[0] - candidate.origin_shift[0],
        input.origin_shift[1] - candidate.origin_shift[1],
        input.origin_shift[2] - candidate.origin_shift[2],
    ];
    let mut origin_shift = mat_multiply_matrix_vector_d3(&inverse_candidate, &shift_difference);
    for value in &mut origin_shift {
        *value = mat_dmod1(*value);
    }

    let transformed = get_distinct_changed_magnetic_symmetry(
        &transformation_matrix,
        &origin_shift,
        magnetic_symmetry,
    )?;
    let database_symmetry = get_spacegroup_operations(candidate.uni_number, candidate.hall_number)?;
    if !is_equal(&transformed, &database_symmetry, symprec) {
        return None;
    }

    let mut std_rotation_matrix = [[0.0; 3]; 3];
    get_rigid_rotation(
        &mut std_rotation_matrix,
        lattice,
        &input.transformation_matrix,
        &input.reference_spacegroup,
    );

    Some(MagneticDataset {
        uni_number: candidate.uni_number,
        msg_type: get_magnetic_spacegroup_type(candidate.uni_number).type_,
        hall_number: candidate.hall_number,
        transformation_matrix,
        origin_shift,
        std_rotation_matrix,
    })
}

/// 获取参考空间群和变换后的磁性对称操作。
/// 从磁性对称中获取参考空间群、变换后的磁性对称操作、变换矩阵和类型。
///
/// 对应 C 原版的 `get_reference_space_group`。
fn get_reference_space_group(
    lattice: &Mat3,
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
        get_magnetic_space_group_type(magnetic_symmetry, sym_fsg.len(), sym_xsg.len())?;
    let representatives = build_representatives(msgtype_num, magnetic_symmetry)?;

    // 4. 选择参考设置: type-4 用 XSG, 其他用 FSG
    //    C 原版对 type-4 使用 xsg 作为 ref_sg
    let ref_sg = if msgtype_num == MagneticType::AntiTranslation {
        &mut _xsg
    } else {
        &mut fsg
    };

    // 5. Refine the reference basis against the physical Cartesian metric,
    //    then form x_std = (tmat, shift) x. The Rust space-group search does
    //    not retain upstream's complete orig_lattice context, so this second
    //    refinement is required for non-cubic structure inputs.
    let lattice_inv = mat_inverse_matrix_d3(lattice, 0.0).ok()?;
    ref_sg.bravais_lattice = mat_multiply_matrix_d3(lattice, &ref_sg.bravais_lattice);
    find_similar_bravais_lattice(ref_sg, symprec);
    ref_sg.bravais_lattice = mat_multiply_matrix_d3(&lattice_inv, &ref_sg.bravais_lattice);
    let tmat = mat_inverse_matrix_d3(&ref_sg.bravais_lattice, 0.0).ok()?;
    let shift = ref_sg.origin_shift;

    // 6. 合成变换后的磁性对称操作
    //    (C 原版: get_changed_magnetic_symmetry 分解 + 重合成)
    let changed_symmetry = get_changed_magnetic_symmetry(
        &tmat,
        &shift,
        &representatives,
        &sym_xsg,
        magnetic_symmetry,
        symprec,
    )?;

    // 7. 复制 ref_sg 用于返回
    let ref_sg_copy = ref_sg.clone();

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
        get_magnetic_space_group_type(magnetic_symmetry, sym_fsg.len(), sym_xsg.len())?;

    // 3. 用非磁 Hall 编号构建参考 Spacegroup
    let spg_type = get_spacegroup_type(parent_hall_number);
    let ref_sg = Spacegroup::from_spg_type(parent_hall_number, [0.0; 3], *lattice, &spg_type);

    // 4. 计算 changed_symmetry: 使用完整的合成（representatives × pure_trans × factors）
    let tmat = ref_sg.bravais_lattice;
    let shift = ref_sg.origin_shift;
    let representatives = build_representatives(msgtype_num, magnetic_symmetry)?;
    let changed_symmetry = get_changed_magnetic_symmetry(
        &tmat,
        &shift,
        &representatives,
        &sym_xsg,
        magnetic_symmetry,
        symprec,
    )?;

    // 5. 复制 ref_sg 用于返回
    let ref_sg_copy = ref_sg.clone();

    Some((ref_sg_copy, changed_symmetry, tmat, shift, msgtype_num))
}

fn deduplicate_spatial_operations(sym: Symmetry, symprec: f64) -> Symmetry {
    let mut dedup = Symmetry::with_capacity(sym.len());
    for i in 0..sym.len() {
        let duplicate = (0..dedup.len()).any(|j| {
            if !mat_check_identity_matrix_i3(&dedup.rot[j], &sym.rot[i]) {
                return false;
            }
            let diff: f64 = (0..3)
                .map(|k| {
                    let x = dedup.trans[j][k] - sym.trans[i][k];
                    (x - x.round()).abs()
                })
                .sum();
            diff < symprec
        });
        if !duplicate {
            dedup.push(sym.rot[i], sym.trans[i]);
        }
    }
    dedup
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

    let mut sym = Symmetry::with_capacity(magnetic_symmetry.len());
    for i in 0..magnetic_symmetry.len() {
        if (!ignore_time_reversal || is_type2) && magnetic_symmetry.timerev[i] {
            continue;
        }
        sym.push(magnetic_symmetry.rot[i], magnetic_symmetry.trans[i]);
    }

    if ignore_time_reversal || is_type2 {
        sym = deduplicate_spatial_operations(sym, symprec);
    }

    if sym.is_empty() { None } else { Some(sym) }
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

/// Get space group from magnetic symmetry.
///
/// `ignore_time_reversal=true` → FSG (family space group, all ops regardless of timerev)
/// `ignore_time_reversal=false` → XSG (maximal subspace group, only ordinary ops)
///
/// Returns (Spacegroup, Symmetry) pair.
pub(crate) fn get_space_group_with_magnetic_symmetry(
    magnetic_symmetry: &MagneticSymmetry,
    ignore_time_reversal: bool,
    symprec: f64,
) -> Option<(Spacegroup, Symmetry)> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let unit_lat: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let num_sym_msg = magnetic_symmetry.len();

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

    // Extract operations. For type-II MSGs, primed copies are skipped
    // together with the time-reversal-erased duplicates.
    let mut sym = Symmetry::with_capacity(num_sym_msg);
    for i in 0..num_sym_msg {
        if (!ignore_time_reversal || is_type2) && magnetic_symmetry.timerev[i] {
            continue;
        }
        sym.push(magnetic_symmetry.rot[i], magnetic_symmetry.trans[i]);
    }

    if sym.is_empty() {
        return None;
    }

    // Deduplicate ops when ignoring time reversal: same (R, t) across timerev=0/1
    // duplicates would confuse the space group search (e.g. BCC AFM: 24 ops with
    // 12 ordinary + 12 anti → FSG should be 12 unique ops, not 24).
    if ignore_time_reversal || is_type2 {
        sym = deduplicate_spatial_operations(sym, symprec);
    }

    // Get primitive symmetry: (a, b, c) = (a_prim, b_prim, c_prim) @ tmat
    let (tmat, prim_sym) = get_primitive_symmetry(&sym, symprec)?;

    let mut spacegroup = match search_spacegroup_with_symmetry(&prim_sym, &unit_lat, symprec) {
        Ok(sg) => sg,
        Err(_) => {
            // 标准空间群搜索失败 → 使用 fallback
            return find_spacegroup_by_symmetry(&sym, &unit_lat, symprec).map(|sg| (sg, sym));
        }
    };

    // Refine bravais lattice and origin_shift
    find_similar_bravais_lattice(&mut spacegroup, symprec);

    // Change basis from primitive to original:
    // x = (tmat, 0)^-1 x_prim
    // => x_std = (P^-1, p) (tmat, 0) x = ( P^-1 @ tmat, p) x
    //    (a_std, b_std, c_std) = (a, b, c) @ tmat^-1 @ P
    let inv_tmat = mat_inverse_matrix_d3(&tmat, 0.0).ok()?;
    spacegroup.bravais_lattice = mat_multiply_matrix_d3(&inv_tmat, &spacegroup.bravais_lattice);

    Some((spacegroup, sym))
}

/// Fallback: 直接从对称操作匹配 Hall 编号，绕过完整的空间群搜索。
/// 当 `search_spacegroup_with_symmetry` 失败时使用。
fn find_spacegroup_by_symmetry(
    symmetry: &Symmetry,
    lattice: &Mat3,
    symprec: f64,
) -> Option<Spacegroup> {
    let mut origin_shift = [0.0; 3];

    let (tmat_int, pointgroup) = get_transformation_matrix(&symmetry.rot, None)?;

    let mut correction_mat = [[0.0; 3]; 3];
    let centering = get_centering(&mut correction_mat, &tmat_int, pointgroup.laue);
    if centering == Centering::Error {
        return None;
    }

    let tmat = mat_multiply_matrix_id3(&tmat_int, &correction_mat);
    let conv_lattice = mat_multiply_matrix_d3(lattice, &tmat);

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
            let spg_type = get_spacegroup_type(hall as usize);
            return Some(Spacegroup::from_spg_type(
                hall as usize,
                origin_shift,
                conv_lattice,
                &spg_type,
            ));
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
            let mut rep = MagneticSymmetry::with_capacity(1);
            rep.push(identity, [0.0; 3], false);
            Some(rep)
        }
        MagneticType::Grey => {
            let mut rep = MagneticSymmetry::with_capacity(2);
            rep.push(identity, [0.0; 3], false);
            rep.push(identity, [0.0; 3], true);
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
) -> Option<MagneticType> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    if num_sym_fsg == num_sym_xsg {
        let num_sym_msg = magnetic_symmetry.len();
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
        let representative = get_representative(magnetic_symmetry)?;
        if representative.len() != 2 {
            return None;
        }
        if mat_check_identity_matrix_i3(&identity, &representative.rot[1]) {
            Some(MagneticType::AntiTranslation)
        } else {
            Some(MagneticType::BlackWhite)
        }
    } else {
        None
    }
}

/// Get coset representative of XSG in MSG.
/// For type-III and type-IV MSGs. Returns identity for type-I/II.
fn get_representative(magnetic_symmetry: &MagneticSymmetry) -> Option<MagneticSymmetry> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    let mut representative = MagneticSymmetry::with_capacity(2);
    representative.push(identity, [0.0; 3], false);

    // A primed operation with identity linear part is the Type-IV
    // anti-translation representative. Its translation is generally nonzero.
    let antiunitary = magnetic_symmetry
        .timerev
        .iter()
        .position(|&timerev| timerev)?;
    let chosen = (0..magnetic_symmetry.len())
        .filter(|&i| magnetic_symmetry.timerev[i])
        .find(|&i| mat_check_identity_matrix_i3(&identity, &magnetic_symmetry.rot[i]))
        .unwrap_or(antiunitary);

    representative.push(
        magnetic_symmetry.rot[chosen],
        magnetic_symmetry.trans[chosen],
        true,
    );
    Some(representative)
}

/// Apply `x_std = (tmat, shift) x` to magnetic symmetry, deduplicating.
///
/// R_std = T * R * T^-1
/// t_std = shift - R_std * shift + T * t
fn get_distinct_changed_magnetic_symmetry(
    tmat: &Mat3,
    shift: &Vec3,
    sym_msg: &MagneticSymmetry,
) -> Option<MagneticSymmetry> {
    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0).ok()?;

    let mut changed = MagneticSymmetry::with_capacity(sym_msg.len());

    for i in 0..sym_msg.len() {
        // R_std = T * R * T^-1
        let rot_f64 = mat_cast_matrix_3i_to_3d(&sym_msg.rot[i]);
        let tmp = mat_multiply_matrix_d3(tmat, &rot_f64);
        let r_new = mat_multiply_matrix_d3(&tmp, &inv_tmat);

        // Round to integer rotation matrix
        let rot_i = [
            [
                mat_nint(r_new[0][0]),
                mat_nint(r_new[0][1]),
                mat_nint(r_new[0][2]),
            ],
            [
                mat_nint(r_new[1][0]),
                mat_nint(r_new[1][1]),
                mat_nint(r_new[1][2]),
            ],
            [
                mat_nint(r_new[2][0]),
                mat_nint(r_new[2][1]),
                mat_nint(r_new[2][2]),
            ],
        ];

        // t_std = shift - R_std * shift + T * t
        let rotated_shift = mat_multiply_matrix_vector_id3(&rot_i, shift);
        let transformed_trans = mat_multiply_matrix_vector_d3(tmat, &sym_msg.trans[i]);
        let mut t_new = [0.0; 3];
        for j in 0..3 {
            t_new[j] = mat_dmod1(shift[j] - rotated_shift[j] + transformed_trans[j]);
        }

        // Check for uniqueness (same rotation, same translation, same timerev)
        let is_dup = (0..changed.len()).any(|j| {
            if !mat_check_identity_matrix_i3(&changed.rot[j], &rot_i) {
                return false;
            }
            let mut diff = [0.0; 3];
            for k in 0..3 {
                diff[k] = changed.trans[j][k] - t_new[k];
                diff[k] -= mat_nint(diff[k]) as f64;
            }
            diff.iter().all(|value| value.abs() < 1e-5) && changed.timerev[j] == sym_msg.timerev[i]
        });

        if !is_dup {
            changed.push(rot_i, t_new, sym_msg.timerev[i]);
        }
    }

    Some(changed)
}

/// 检查旋转矩阵 `a` 是否已包含在 `sym_msg` 中。
fn is_contained_mat(a: &Mat3I, sym_msg: &MagneticSymmetry) -> bool {
    sym_msg
        .rot
        .iter()
        .any(|rot| mat_check_identity_matrix_i3(a, rot))
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
        if eq {
            return true;
        }
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
    if !det.is_finite() || det == 0.0 {
        return None;
    }

    let size_f = pure_trans.len() as f64 / det.abs();
    let rounded_size = size_f.round();
    let rounding_tolerance = symprec.max(16.0 * f64::EPSILON * size_f.abs());
    if !size_f.is_finite()
        || (size_f - rounded_size).abs() > rounding_tolerance
        || rounded_size < 0.0
        || rounded_size > MAX_CHANGED_PURE_TRANSLATIONS as f64
    {
        return None;
    }
    let size = rounded_size as usize;

    let mut changed: Vec<Vec3> = Vec::new();
    changed.try_reserve_exact(size).ok()?;

    let is_integer_unimodular = (det.abs() - 1.0).abs() <= symprec
        && tmat
            .iter()
            .flatten()
            .all(|value| value.is_finite() && (value - value.round()).abs() <= symprec);
    if is_integer_unimodular {
        for pt in pure_trans {
            let trans = mat_multiply_matrix_vector_d3(tmat, pt);
            changed.push([
                mat_dmod1(trans[0]),
                mat_dmod1(trans[1]),
                mat_dmod1(trans[2]),
            ]);
        }
    } else {
        // 查找转动矩阵元素的最小公分母
        let mut denominator = 1;
        loop {
            let mut ok = true;
            'matrix: for row in tmat {
                for &value in row {
                    if (value * denominator as f64 - mat_nint(value * denominator as f64) as f64)
                        .abs()
                        > symprec
                    {
                        ok = false;
                        break 'matrix;
                    }
                }
            }
            if ok {
                break;
            }
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
                        let shifted = [pt[0] + n0 as f64, pt[1] + n1 as f64, pt[2] + n2 as f64];
                        let trans = mat_multiply_matrix_vector_d3(tmat, &shifted);
                        let t_mod = [
                            mat_dmod1(trans[0]),
                            mat_dmod1(trans[1]),
                            mat_dmod1(trans[2]),
                        ];

                        if !is_contained_vec(&t_mod, &changed, symprec) {
                            changed.push(t_mod);
                        }
                    }
                }
            }
        }
    }

    if changed.len() != size {
        return None;
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
        get_distinct_changed_magnetic_symmetry(tmat, shift, representatives)?;

    // 2. 收集原始磁性对称中的纯平移（仅 timerev=0），变换到参考设置
    let pure_trans =
        crate::spin::collect_pure_translations_from_magnetic_symmetry(magnetic_symmetry);
    let changed_pure_trans = get_changed_pure_translations(tmat, &pure_trans, symprec)?;

    // 3. 从 XSG 对称性中收集因子群（仅去重旋转部分，timerev=0）
    let mut factors = MagneticSymmetry::with_capacity(sym_xsg.len());
    for i in 0..sym_xsg.len() {
        if !is_contained_mat(&sym_xsg.rot[i], &factors) {
            factors.push(sym_xsg.rot[i], sym_xsg.trans[i], false);
        }
    }
    let num_factors = factors.len();
    let changed_factors = get_distinct_changed_magnetic_symmetry(tmat, shift, &factors)?;

    // 4. 合成: (I, ti)(Pj, tj)(Pk, tk) = (Pj * Pk, Pj * tk + tj + ti)
    let size = changed_representatives.len() * changed_pure_trans.len() * num_factors;
    let mut changed = MagneticSymmetry::with_capacity(size);

    for pure_translation in &changed_pure_trans {
        for j in 0..changed_representatives.len() {
            for k in 0..num_factors {
                // R = Pj * Pk
                let rot = mat_multiply_matrix_i3(
                    &changed_representatives.rot[j],
                    &changed_factors.rot[k],
                );

                // t = Pj * tk + tj + ti
                let mut trans = mat_multiply_matrix_vector_id3(
                    &changed_representatives.rot[j],
                    &changed_factors.trans[k],
                );
                for (s, value) in trans.iter_mut().enumerate() {
                    *value += changed_representatives.trans[j][s] + pure_translation[s];
                    *value = mat_dmod1(*value);
                }

                // timerev = changed_representatives.timerev XOR changed_factors.timerev
                // (factors 都是 ordinary，所以 XOR 就是 representatives 的 timerev)
                let timerev = changed_representatives.timerev[j] != changed_factors.timerev[k];
                changed.push(rot, trans, timerev);
            }
        }
    }

    Some(changed)
}

/// 检查两个磁对称操作集合是否在周期平移意义下完全相等。
fn is_equal(sym1: &MagneticSymmetry, sym2: &MagneticSymmetry, symprec: f64) -> bool {
    if sym1.len() != sym2.len() {
        return false;
    }

    let mut found = vec![false; sym2.len()];
    for i in 0..sym1.len() {
        let mut matched = false;
        for (j, already_found) in found.iter_mut().enumerate() {
            if *already_found {
                continue;
            }
            if !mat_check_identity_matrix_i3(&sym1.rot[i], &sym2.rot[j]) {
                continue;
            }
            if sym1.timerev[i] != sym2.timerev[j] {
                continue;
            }
            let mut diff = [0.0; 3];
            for (k, value) in diff.iter_mut().enumerate() {
                *value = sym1.trans[i][k] - sym2.trans[j][k];
                *value -= mat_nint(*value) as f64;
            }
            if diff[0].abs() < symprec && diff[1].abs() < symprec && diff[2].abs() < symprec {
                *already_found = true;
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
fn get_rigid_rotation(rigid_rot: &mut Mat3, lattice: &Mat3, tmat: &Mat3, ref_sg: &Spacegroup) {
    let inv_tmat = mat_inverse_matrix_d3(tmat, 0.0).ok();
    if let Some(inv) = inv_tmat {
        let tmp = mat_multiply_matrix_d3(&ref_sg.bravais_lattice, &inv);
        let inv_lat = mat_inverse_matrix_d3(lattice, 0.0).ok();
        if let Some(inv_l) = inv_lat {
            let result = mat_multiply_matrix_d3(&inv_l, &tmp);
            *rigid_rot = result;
        }
    }
}

// ============================================================================
// 内部测试: DB 匹配算法、边界情况
// ============================================================================

#[cfg(test)]
mod tests {
    use crate::MagneticType;
    use crate::mathfunc::{Mat3, Mat3I, Vec3, is_proper};
    use crate::msg_database::get_magnetic_spacegroup_type;
    use crate::symmetry::MagneticSymmetry;

    const SYMPREC: f64 = 1e-5;

    fn cubic_lattice() -> Mat3 {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    fn pm3m_ops() -> Vec<(Mat3I, Vec3)> {
        let (count, start) = crate::spg_database::get_operation_index(517);
        (0..count)
            .filter_map(|i| crate::spg_database::get_operation_by_index(start + i))
            .collect()
    }

    fn make_mag_sym(timerev: &[bool], ops: &[(Mat3I, Vec3)]) -> MagneticSymmetry {
        assert_eq!(timerev.len(), ops.len());
        let mut sym = MagneticSymmetry::with_capacity(ops.len());
        for ((rot, trans), timerev) in ops.iter().zip(timerev) {
            sym.push(*rot, *trans, *timerev);
        }
        sym
    }

    #[test]
    fn changed_pure_translations_rejects_unsafe_determinants() {
        let translations = [[0.0, 0.0, 0.0]];
        let singular = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]];
        let nan = [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let infinite = [[f64::INFINITY, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let tiny = [[1e-12, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        for transform in [singular, nan, infinite, tiny] {
            assert!(
                super::get_changed_pure_translations(&transform, &translations, SYMPREC).is_none()
            );
        }
    }

    #[test]
    fn changed_pure_translations_handles_negative_determinants() {
        let translations = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let reflection = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let doubled = [[-2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        let reflected =
            super::get_changed_pure_translations(&reflection, &translations, SYMPREC).unwrap();
        assert_eq!(reflected.len(), 2);
        assert!(super::is_contained_vec(
            &[0.5, 0.0, 0.0],
            &reflected,
            SYMPREC
        ));

        let collapsed =
            super::get_changed_pure_translations(&doubled, &translations, SYMPREC).unwrap();
        assert_eq!(collapsed, vec![[0.0, 0.0, 0.0]]);
    }

    #[test]
    fn changed_pure_translations_expands_fractional_basis() {
        let translations = [[0.0, 0.0, 0.0]];
        let halved = [[0.5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let expanded =
            super::get_changed_pure_translations(&halved, &translations, SYMPREC).unwrap();

        assert_eq!(expanded.len(), 2);
        assert!(super::is_contained_vec(
            &[0.0, 0.0, 0.0],
            &expanded,
            SYMPREC
        ));
        assert!(super::is_contained_vec(
            &[0.5, 0.0, 0.0],
            &expanded,
            SYMPREC
        ));
    }

    #[test]
    fn changed_pure_translations_rejects_nonintegral_multiplicity() {
        let translations = [[0.0, 0.0, 0.0]];
        let transform = [[0.3, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        assert!(super::get_changed_pure_translations(&transform, &translations, SYMPREC).is_none());
    }

    #[test]
    fn changed_pure_translations_rejects_incompatible_fractional_unit_determinant() {
        let translations = [[0.0, 0.0, 0.0]];
        let transform = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]];

        assert!(super::get_changed_pure_translations(&transform, &translations, SYMPREC).is_none());
    }

    /// Type-1 (Ordinary): 所有 timerev=false
    #[test]
    fn test_db_type1() {
        let ops = pm3m_ops();
        let mag_sym = make_mag_sym(&vec![false; ops.len()], &ops);
        let ds = super::identify_magnetic_space_group_type(&cubic_lattice(), &mag_sym, SYMPREC)
            .expect("must match");
        assert_eq!(ds.msg_type, MagneticType::Ordinary);
        assert_eq!(ds.hall_number, 517);
        assert_eq!(
            get_magnetic_spacegroup_type(ds.uni_number).type_,
            MagneticType::Ordinary
        );
    }

    /// Type-2 (Grey): 每个操作加倍 (timerev=false + timerev=true)
    #[test]
    fn test_db_type2() {
        let ops = pm3m_ops();
        let mut mag_sym = MagneticSymmetry::with_capacity(ops.len() * 2);
        for &(rot, trans) in &ops {
            mag_sym.push(rot, trans, false);
            mag_sym.push(rot, trans, true);
        }
        let ds = super::identify_magnetic_space_group_type(&cubic_lattice(), &mag_sym, SYMPREC)
            .expect("must match");
        assert_eq!(ds.msg_type, MagneticType::Grey);
        assert_eq!(ds.hall_number, 517);
        assert_eq!(
            get_magnetic_spacegroup_type(ds.uni_number).type_,
            MagneticType::Grey
        );
    }

    /// Type-3 (BlackWhite): 非正当旋转带 timerev=true
    #[test]
    fn test_db_type3() {
        let ops = pm3m_ops();
        let timerev: Vec<bool> = ops.iter().map(|(r, _)| !is_proper(r)).collect();
        let mag_sym = make_mag_sym(&timerev, &ops);
        let ds = super::identify_magnetic_space_group_type(&cubic_lattice(), &mag_sym, SYMPREC)
            .expect("must match");
        assert_eq!(ds.msg_type, MagneticType::BlackWhite);
        assert_eq!(ds.hall_number, 517);
        assert_eq!(
            get_magnetic_spacegroup_type(ds.uni_number).type_,
            MagneticType::BlackWhite
        );
    }

    /// 空对称操作 → 返回 None
    #[test]
    fn test_empty_symmetry() {
        let mag_sym = MagneticSymmetry::new();
        assert!(
            super::identify_magnetic_space_group_type(&cubic_lattice(), &mag_sym, SYMPREC,)
                .is_err()
        );
    }

    /// 缺少单位操作 → 返回 Err
    #[test]
    fn test_no_identity() {
        let mut mag_sym = MagneticSymmetry::with_capacity(1);
        mag_sym.push([[0, -1, 0], [1, 0, 0], [0, 0, 1]], [0.0; 3], false);
        assert!(
            super::identify_magnetic_space_group_type(&cubic_lattice(), &mag_sym, SYMPREC,)
                .is_err()
        );
    }

    #[test]
    fn type_iv_parent_hall_disambiguates_a_changed_basis_and_origin() {
        let database = crate::msg_database::get_spacegroup_operations(282, 182).unwrap();
        let input_transform = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]];
        let input_shift = [0.137, 0.219, 0.311];
        let input = super::get_distinct_changed_magnetic_symmetry(
            &input_transform,
            &input_shift,
            &database,
        )
        .unwrap();

        assert!(matches!(
            super::identify_magnetic_space_group_type(&cubic_lattice(), &input, SYMPREC),
            Err(crate::SymError::MagneticUniAmbiguous)
        ));

        let bns_37 =
            super::identify_with_parent_hall(&cubic_lattice(), &input, Some(182), SYMPREC).unwrap();
        assert_eq!(bns_37.uni_number, 282);

        let bns_36 =
            super::identify_with_parent_hall(&cubic_lattice(), &input, Some(176), SYMPREC).unwrap();
        assert_eq!(bns_36.uni_number, 275);
    }

    #[test]
    fn type_iv_database_has_only_two_cross_uni_canonical_classes() {
        let mut cross_uni_classes = Vec::new();
        for candidates in super::type_iv_canonical_index().values() {
            let mut unis: Vec<_> = candidates
                .iter()
                .map(|candidate| candidate.uni_number)
                .collect();
            unis.sort_unstable();
            unis.dedup();
            if unis.len() > 1 {
                cross_uni_classes.push(unis);
            }
        }
        cross_uni_classes.sort();
        cross_uni_classes.dedup();

        assert_eq!(cross_uni_classes, vec![vec![275, 282], vec![277, 284]]);
    }

    #[test]
    fn type_iv_orthorhombic_metric_recovers_unique_283_and_reports_real_ambiguities() {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.3, 0.0], [0.0, 0.0, 1.7]];

        for uni in [282usize, 283, 284] {
            for hall in 182usize..=184 {
                let magnetic = crate::msg_database::get_spacegroup_operations(uni, hall).unwrap();
                let automatic =
                    super::identify_magnetic_space_group_type(&lattice, &magnetic, SYMPREC);
                if uni == 283 {
                    assert_eq!(automatic.unwrap().uni_number, 283, "input Hall {hall}");
                } else {
                    assert!(
                        matches!(automatic, Err(crate::SymError::MagneticUniAmbiguous)),
                        "UNI {uni} Hall {hall} must not be silently mapped to its analogue"
                    );
                }

                let with_parent =
                    super::identify_with_parent_hall(&lattice, &magnetic, Some(hall), SYMPREC)
                        .unwrap();
                assert_eq!(with_parent.uni_number, uni);
                assert_eq!(with_parent.hall_number, hall);
            }
        }
    }

    #[test]
    #[ignore = "focused upstream-comparison diagnostic"]
    fn diagnose_selected_database_reference_groups() {
        for uni in [132usize, 282, 667, 751, 890, 1338] {
            let hall = crate::msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni][1] as usize;
            let magnetic = crate::msg_database::get_spacegroup_operations(uni, hall).unwrap();
            let (fsg, sym_fsg) =
                super::get_family_space_group_with_magnetic_symmetry(&magnetic, SYMPREC).unwrap();
            let (xsg, sym_xsg) =
                super::get_maximal_subspace_group_with_magnetic_symmetry(&magnetic, SYMPREC)
                    .unwrap();
            let reference =
                super::get_reference_space_group(&cubic_lattice(), &magnetic, SYMPREC).unwrap();
            let result =
                super::identify_with_parent_hall(&cubic_lattice(), &magnetic, Some(hall), SYMPREC);

            eprintln!(
                "UNI {uni} input Hall {hall}: FSG Hall {} SG{} order {} P {:?} p {:?}; \
                 XSG Hall {} SG{} order {} P {:?} p {:?}; \
                 ref Hall {} type {:?} changed {} tmat {:?} shift {:?}",
                fsg.hall_number,
                fsg.number,
                sym_fsg.len(),
                fsg.bravais_lattice,
                fsg.origin_shift,
                xsg.hall_number,
                xsg.number,
                sym_xsg.len(),
                xsg.bravais_lattice,
                xsg.origin_shift,
                reference.0.hall_number,
                reference.4,
                reference.1.len(),
                reference.2,
                reference.3,
            );
            match result {
                Ok(dataset) => eprintln!(
                    "  result UNI {} type {:?} Hall {} tmat {:?} shift {:?}",
                    dataset.uni_number,
                    dataset.msg_type,
                    dataset.hall_number,
                    dataset.transformation_matrix,
                    dataset.origin_shift,
                ),
                Err(error) => eprintln!("  result error {error:?}"),
            }
        }
    }
}
