//! k 点地址生成。
//!
//! 在不可约布里渊区内生成 k 点坐标和权重，用于态密度和能带计算。

use crate::SymError;
use crate::kgrid;
use crate::mathfunc::{
    Mat3I, mat_check_identity_matrix_i3, mat_dabs, mat_get_determinant_d3, mat_multiply_matrix_i3,
    mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_i3, mat_multiply_matrix_vector_id3,
    mat_nint, mat_norm_squared_d3, mat_transpose_matrix_i3,
};

const KPT_NUM_BZ_SEARCH_SPACE: usize = 125;

// 静态搜索空间数组：用于在倒易空间寻找最近邻的格点
static BZ_SEARCH_SPACE: [[i32; 3]; KPT_NUM_BZ_SEARCH_SPACE] = [
    [0, 0, 0],
    [0, 0, 1],
    [0, 0, 2],
    [0, 0, -2],
    [0, 0, -1],
    [0, 1, 0],
    [0, 1, 1],
    [0, 1, 2],
    [0, 1, -2],
    [0, 1, -1],
    [0, 2, 0],
    [0, 2, 1],
    [0, 2, 2],
    [0, 2, -2],
    [0, 2, -1],
    [0, -2, 0],
    [0, -2, 1],
    [0, -2, 2],
    [0, -2, -2],
    [0, -2, -1],
    [0, -1, 0],
    [0, -1, 1],
    [0, -1, 2],
    [0, -1, -2],
    [0, -1, -1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 0, 2],
    [1, 0, -2],
    [1, 0, -1],
    [1, 1, 0],
    [1, 1, 1],
    [1, 1, 2],
    [1, 1, -2],
    [1, 1, -1],
    [1, 2, 0],
    [1, 2, 1],
    [1, 2, 2],
    [1, 2, -2],
    [1, 2, -1],
    [1, -2, 0],
    [1, -2, 1],
    [1, -2, 2],
    [1, -2, -2],
    [1, -2, -1],
    [1, -1, 0],
    [1, -1, 1],
    [1, -1, 2],
    [1, -1, -2],
    [1, -1, -1],
    [2, 0, 0],
    [2, 0, 1],
    [2, 0, 2],
    [2, 0, -2],
    [2, 0, -1],
    [2, 1, 0],
    [2, 1, 1],
    [2, 1, 2],
    [2, 1, -2],
    [2, 1, -1],
    [2, 2, 0],
    [2, 2, 1],
    [2, 2, 2],
    [2, 2, -2],
    [2, 2, -1],
    [2, -2, 0],
    [2, -2, 1],
    [2, -2, 2],
    [2, -2, -2],
    [2, -2, -1],
    [2, -1, 0],
    [2, -1, 1],
    [2, -1, 2],
    [2, -1, -2],
    [2, -1, -1],
    [-2, 0, 0],
    [-2, 0, 1],
    [-2, 0, 2],
    [-2, 0, -2],
    [-2, 0, -1],
    [-2, 1, 0],
    [-2, 1, 1],
    [-2, 1, 2],
    [-2, 1, -2],
    [-2, 1, -1],
    [-2, 2, 0],
    [-2, 2, 1],
    [-2, 2, 2],
    [-2, 2, -2],
    [-2, 2, -1],
    [-2, -2, 0],
    [-2, -2, 1],
    [-2, -2, 2],
    [-2, -2, -2],
    [-2, -2, -1],
    [-2, -1, 0],
    [-2, -1, 1],
    [-2, -1, 2],
    [-2, -1, -2],
    [-2, -1, -1],
    [-1, 0, 0],
    [-1, 0, 1],
    [-1, 0, 2],
    [-1, 0, -2],
    [-1, 0, -1],
    [-1, 1, 0],
    [-1, 1, 1],
    [-1, 1, 2],
    [-1, 1, -2],
    [-1, 1, -1],
    [-1, 2, 0],
    [-1, 2, 1],
    [-1, 2, 2],
    [-1, 2, -2],
    [-1, 2, -1],
    [-1, -2, 0],
    [-1, -2, 1],
    [-1, -2, 2],
    [-1, -2, -2],
    [-1, -2, -1],
    [-1, -1, 0],
    [-1, -1, 1],
    [-1, -1, 2],
    [-1, -1, -2],
    [-1, -1, -1],
];

/// 获取不可约倒易网格 (Irreducible Reciprocal Mesh)
///
/// # Arguments
/// * `grid_address` - 输出：网格点坐标
/// * `ir_mapping_table` - 输出：映射表，将每个网格点映射到其不可约代表点
/// * `mesh` - 网格尺寸 [Nx, Ny, Nz]
/// * `is_shift` - 网格位移 (Monkhorst-Pack shift)
/// * `rot_reciprocal` - 倒易空间中的点群旋转矩阵
pub(crate) fn get_irreducible_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &[Mat3I],
) -> Result<usize, SymError> {
    get_dense_irreducible_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        rot_reciprocal,
    )
}

pub(crate) fn get_dense_irreducible_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &[Mat3I],
) -> Result<usize, SymError> {
    get_dense_ir_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        rot_reciprocal,
    )
}

/// 获取考虑了时间反演和 q 点的稳定倒易网格
pub(crate) fn get_stabilized_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    is_time_reversal: i32,
    rotations: &[Mat3I],
    qpoints: &[[f64; 3]],
) -> Result<usize, SymError> {
    kgrid::validate_mesh(mesh)?;
    kgrid::validate_shift(is_shift)?;
    // 获取倒易空间点群（包含时间反演）
    let rot_reciprocal = get_point_group_reciprocal(rotations, is_time_reversal)
        .ok_or(SymError::ArraySizeShortage)?;

    // 计算容差
    let mesh_sum = i64::from(mesh[0]) + i64::from(mesh[1]) + i64::from(mesh[2]);
    let tolerance = 0.01 / mesh_sum as f64;

    // 获取稳定化 q 点后的点群
    let rot_reciprocal_q = get_point_group_reciprocal_with_q(&rot_reciprocal, tolerance, qpoints)
        .ok_or(SymError::ArraySizeShortage)?;

    get_dense_ir_reciprocal_mesh(
        grid_address,
        ir_mapping_table,
        mesh,
        is_shift,
        &rot_reciprocal_q,
    )
}

/// 将网格点重定位到第一布里渊区 (First Brillouin Zone)
pub(crate) fn relocate_bz_grid_address(
    bz_grid_address: &mut [[i32; 3]],
    bz_map: &mut [usize],
    grid_address: &[[i32; 3]],
    mesh: &[i32; 3],
    rec_lattice: &[[f64; 3]; 3],
    is_shift: &[i32; 3],
) -> Result<usize, SymError> {
    let total = kgrid::validate_mesh(mesh)?;
    kgrid::validate_shift(is_shift)?;
    let num_bz_map = total.checked_mul(8).ok_or(SymError::ArraySizeShortage)?;
    if bz_grid_address.len() < num_bz_map || bz_map.len() < num_bz_map || grid_address.len() < total
    {
        return Err(SymError::ArraySizeShortage);
    }
    // 使用 vec! 宏分配内存，对应 C 的 malloc
    let mut dense_bz_map = vec![0; num_bz_map];

    let num_bzgp = relocate_dense_bz_grid_address(
        bz_grid_address,
        &mut dense_bz_map,
        grid_address,
        mesh,
        rec_lattice,
        is_shift,
    )?;

    for i in 0..num_bz_map {
        if dense_bz_map[i] == num_bz_map {
            bz_map[i] = usize::MAX; // 对应 C 中的 -1 (size_t)
        } else {
            bz_map[i] = dense_bz_map[i];
        }
    }

    Ok(num_bzgp)
}

/// 对原始网格地址应用旋转，获取所有旋转后的双倍网格点索引。
pub(crate) fn get_dense_grid_points_by_rotations(
    rot_grid_points: &mut [usize],
    address_orig: &[i32; 3],
    rot_reciprocal: &[Mat3I],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
) -> Result<(), SymError> {
    kgrid::validate_mesh(mesh)?;
    kgrid::validate_shift(is_shift)?;
    if rot_grid_points.len() < rot_reciprocal.len() {
        return Err(SymError::ArraySizeShortage);
    }
    let mut address_double_orig = [0i32; 3];
    kgrid::get_grid_address_double_mesh(
        &mut address_double_orig,
        address_orig,
        mesh,
        is_shift,
    )?;
    for (output, rotation) in rot_grid_points
        .iter_mut()
        .zip(rot_reciprocal)
        .take(rot_reciprocal.len())
    {
        let address_double = mat_multiply_matrix_vector_i3(rotation, &address_double_orig);
        *output = kgrid::get_dense_grid_point_double_mesh(&address_double, mesh)?;
    }
    Ok(())
}

/// 对原始网格地址应用旋转，获取旋转后在 BZ 映射中的双倍网格点索引。
pub(crate) fn get_dense_bz_grid_points_by_rotations(
    rot_grid_points: &mut [usize],
    address_orig: &[i32; 3],
    rot_reciprocal: &[Mat3I],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    bz_map: &[usize],
) -> Result<(), SymError> {
    let total = kgrid::validate_mesh(mesh)?;
    kgrid::validate_shift(is_shift)?;
    let num_bz_map = total.checked_mul(8).ok_or(SymError::ArraySizeShortage)?;
    if rot_grid_points.len() < rot_reciprocal.len() || bz_map.len() < num_bz_map {
        return Err(SymError::ArraySizeShortage);
    }
    let mut address_double_orig = [0i32; 3];
    let bzmesh = [
        mesh[0].checked_mul(2).ok_or(SymError::ArraySizeShortage)?,
        mesh[1].checked_mul(2).ok_or(SymError::ArraySizeShortage)?,
        mesh[2].checked_mul(2).ok_or(SymError::ArraySizeShortage)?,
    ];
    kgrid::get_grid_address_double_mesh(
        &mut address_double_orig,
        address_orig,
        mesh,
        is_shift,
    )?;
    for (output, rotation) in rot_grid_points
        .iter_mut()
        .zip(rot_reciprocal)
        .take(rot_reciprocal.len())
    {
        let address_double = mat_multiply_matrix_vector_i3(rotation, &address_double_orig);
        let bz_index = kgrid::get_dense_grid_point_double_mesh(&address_double, &bzmesh)?;
        *output = bz_map[bz_index];
    }
    Ok(())
}

// --- Internal Logic ---

/// 获取倒易空间点群。
pub(crate) fn get_point_group_reciprocal(rotations: &[Mat3I], is_time_reversal: i32) -> Option<Vec<Mat3I>> {
    let inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];
    let size = if is_time_reversal != 0 {
        rotations.len() * 2
    } else {
        rotations.len()
    };

    let mut rot_reciprocal = vec![[[0; 3]; 3]; size];
    let mut unique_rot = vec![-1; size];

    for i in 0..rotations.len() {
        // 倒易空间的旋转矩阵是实空间旋转矩阵的转置
        let t = mat_transpose_matrix_i3(&rotations[i]);
        rot_reciprocal[i] = t;

        if is_time_reversal != 0 {
            let inv_rot = mat_multiply_matrix_i3(&inversion, &rot_reciprocal[i]);
            rot_reciprocal[rotations.len() + i] = inv_rot;
        }
    }

    // 筛选唯一旋转矩阵
    let mut num_rot = 0;
    for i in 0..rot_reciprocal.len() {
        let mut is_unique = true;
        if unique_rot[..num_rot].iter().any(|&index| {
            mat_check_identity_matrix_i3(
                &rot_reciprocal[index as usize],
                &rot_reciprocal[i],
            )
        }) {
            is_unique = false;
        }
        if is_unique {
            unique_rot[num_rot] = i as i32;
            num_rot += 1;
        }
    }

    let mut rot_return = vec![[[0; 3]; 3]; num_rot];
    for (output, &index) in rot_return.iter_mut().zip(&unique_rot).take(num_rot) {
        *output = rot_reciprocal[index as usize];
    }

    Some(rot_return)
}

/// 考虑 q 点的对称性
fn get_point_group_reciprocal_with_q(
    rot_reciprocal: &[Mat3I],
    symprec: f64,
    qpoints: &[[f64; 3]],
) -> Option<Vec<Mat3I>> {
    let mut ir_rot = vec![-1; rot_reciprocal.len()];
    let mut num_rot = 0;

    for (i, rot) in rot_reciprocal.iter().enumerate() {
        let mut is_all_ok = true;
        for qpoint in qpoints {
            let q_rot = mat_multiply_matrix_vector_id3(rot, qpoint);

            let mut found_diff = false;
            for candidate in qpoints {
                let mut diff = [0.0; 3];
                for l in 0..3 {
                    diff[l] = q_rot[l] - candidate[l];
                    diff[l] -= mat_nint(diff[l]) as f64;
                }
                if mat_dabs(diff[0]) < symprec
                    && mat_dabs(diff[1]) < symprec
                    && mat_dabs(diff[2]) < symprec
                {
                    found_diff = true;
                    break;
                }
            }

            if !found_diff {
                is_all_ok = false;
                break;
            }
        }

        if is_all_ok {
            ir_rot[num_rot] = i as i32;
            num_rot += 1;
        }
    }

    let mut rot_reciprocal_q = vec![[[0; 3]; 3]; num_rot];
    for (output, &index) in rot_reciprocal_q.iter_mut().zip(&ir_rot).take(num_rot) {
        *output = rot_reciprocal[index as usize];
    }

    Some(rot_reciprocal_q)
}

fn get_dense_ir_reciprocal_mesh(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &[Mat3I],
) -> Result<usize, SymError> {
    let total = kgrid::validate_mesh(mesh)?;
    kgrid::validate_shift(is_shift)?;
    if grid_address.len() < total || ir_mapping_table.len() < total {
        return Err(SymError::ArraySizeShortage);
    }
    if check_mesh_symmetry(mesh, is_shift, rot_reciprocal) {
        get_dense_ir_reciprocal_mesh_normal(
            grid_address,
            ir_mapping_table,
            mesh,
            is_shift,
            rot_reciprocal,
        )
    } else {
        get_dense_ir_reciprocal_mesh_distortion(
            grid_address,
            ir_mapping_table,
            mesh,
            is_shift,
            rot_reciprocal,
        )
    }
}

/// 普通网格约化（适用于正交或高对称性网格）
fn get_dense_ir_reciprocal_mesh_normal(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &[Mat3I],
) -> Result<usize, SymError> {
    kgrid::get_all_grid_addresses(grid_address, mesh)?;

    let total_pts = kgrid::validate_mesh(mesh)?;

    // Serial implementation matching the currently supported feature set.
    for i in 0..total_pts {
        let mut address_double = [0; 3];
        kgrid::get_grid_address_double_mesh(
            &mut address_double,
            &grid_address[i],
            mesh,
            is_shift,
        )?;

        ir_mapping_table[i] = i;

        for rot in rot_reciprocal {
            let address_double_rot = mat_multiply_matrix_vector_i3(rot, &address_double);
            let grid_point_rot =
                kgrid::get_dense_grid_point_double_mesh(&address_double_rot, mesh)?;

            if grid_point_rot < ir_mapping_table[i] {
                ir_mapping_table[i] = grid_point_rot;
                break;
            }
        }
    }

    get_dense_num_ir(ir_mapping_table, total_pts)
}

/// 畸变网格约化（适用于非正交网格或低对称性）
fn get_dense_ir_reciprocal_mesh_distortion(
    grid_address: &mut [[i32; 3]],
    ir_mapping_table: &mut [usize],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
    rot_reciprocal: &[Mat3I],
) -> Result<usize, SymError> {
    kgrid::get_all_grid_addresses(grid_address, mesh)?;

    let divisor = [
        i64::from(mesh[1]) * i64::from(mesh[2]),
        i64::from(mesh[2]) * i64::from(mesh[0]),
        i64::from(mesh[0]) * i64::from(mesh[1]),
    ];
    let total_pts = kgrid::validate_mesh(mesh)?;

    for i in 0..total_pts {
        let mut address_double = [0; 3];
        kgrid::get_grid_address_double_mesh(
            &mut address_double,
            &grid_address[i],
            mesh,
            is_shift,
        )?;

        let mut long_address_double = [0i64; 3];
        for j in 0..3 {
            long_address_double[j] = address_double[j] as i64 * divisor[j];
        }

        ir_mapping_table[i] = i;

        for rot in rot_reciprocal {
            let mut long_address_double_rot = [0i64; 3];
            for (k, rotated) in long_address_double_rot.iter_mut().enumerate() {
                *rotated = rot[k][0] as i64 * long_address_double[0]
                    + rot[k][1] as i64 * long_address_double[1]
                    + rot[k][2] as i64 * long_address_double[2];
            }

            let mut indivisible = false;
            let mut address_double_rot = [0; 3];

            for k in 0..3 {
                if long_address_double_rot[k] % divisor[k] != 0 {
                    indivisible = true;
                    break;
                }
                address_double_rot[k] = (long_address_double_rot[k] / divisor[k]) as i32;

                if (address_double_rot[k] % 2 != 0 && is_shift[k] == 0)
                    || (address_double_rot[k] % 2 == 0 && is_shift[k] == 1)
                {
                    indivisible = true;
                    break;
                }
            }

            if indivisible {
                continue;
            }

            let grid_point_rot =
                kgrid::get_dense_grid_point_double_mesh(&address_double_rot, mesh)?;

            if grid_point_rot < ir_mapping_table[i] {
                ir_mapping_table[i] = grid_point_rot;
                break;
            }
        }
    }

    get_dense_num_ir(ir_mapping_table, total_pts)
}

/// 统计不可约点数量并进行路径压缩
fn get_dense_num_ir(ir_mapping_table: &mut [usize], total_pts: usize) -> Result<usize, SymError> {
    if ir_mapping_table.len() < total_pts {
        return Err(SymError::ArraySizeShortage);
    }
    let mut num_ir = 0;

    for (i, &representative) in ir_mapping_table.iter().take(total_pts).enumerate() {
        if representative == i {
            num_ir += 1;
        }
    }

    // Full path compression: a single lookup is insufficient when symmetry
    // discovery creates chains such as 5 -> 3 -> 1 -> 0.
    let representatives = ir_mapping_table[..total_pts].to_vec();
    for representative in &mut ir_mapping_table[..total_pts] {
        let mut root = *representative;
        while representatives[root] != root {
            root = representatives[root];
        }
        *representative = root;
    }

    Ok(num_ir)
}

fn relocate_dense_bz_grid_address(
    bz_grid_address: &mut [[i32; 3]],
    bz_map: &mut [usize],
    grid_address: &[[i32; 3]],
    mesh: &[i32; 3],
    rec_lattice: &[[f64; 3]; 3],
    is_shift: &[i32; 3],
) -> Result<usize, SymError> {
    let total_num_gp = kgrid::validate_mesh(mesh)?;
    kgrid::validate_shift(is_shift)?;
    if rec_lattice.iter().flatten().any(|value| !value.is_finite())
        || mat_get_determinant_d3(rec_lattice).abs() <= f64::EPSILON
    {
        return Err(SymError::InvalidInput);
    }

    let num_bzmesh = total_num_gp
        .checked_mul(8)
        .ok_or(SymError::ArraySizeShortage)?;
    if bz_grid_address.len() < num_bzmesh
        || bz_map.len() < num_bzmesh
        || grid_address.len() < total_num_gp
    {
        return Err(SymError::ArraySizeShortage);
    }
    let tolerance = get_tolerance_for_bz_reduction(rec_lattice, mesh);
    let bzmesh = [
        mesh[0].checked_mul(2).ok_or(SymError::ArraySizeShortage)?,
        mesh[1].checked_mul(2).ok_or(SymError::ArraySizeShortage)?,
        mesh[2].checked_mul(2).ok_or(SymError::ArraySizeShortage)?,
    ];

    // 初始化 bz_map
    bz_map[..num_bzmesh].fill(num_bzmesh);

    let mut boundary_num_gp = 0usize;

    for (i, address) in grid_address.iter().take(total_num_gp).enumerate() {
        let mut distance = [0.0; KPT_NUM_BZ_SEARCH_SPACE];

        // 计算到所有邻近点的距离
        for (j, value) in distance.iter_mut().enumerate() {
            let mut q_vector = [0.0; 3];
            for k in 0..3 {
                let coordinate =
                    i64::from(address[k]) + i64::from(BZ_SEARCH_SPACE[j][k]) * i64::from(mesh[k]);
                q_vector[k] =
                    (coordinate * 2 + i64::from(is_shift[k])) as f64 / (mesh[k] as f64) / 2.0;
            }
            let q_vec_rec = mat_multiply_matrix_vector_d3(rec_lattice, &q_vector);
            *value = mat_norm_squared_d3(&q_vec_rec);
        }

        // 找到最小距离
        let mut min_distance = distance[0];
        let mut min_index = 0;
        for (j, &value) in distance.iter().enumerate().skip(1) {
            if value < min_distance {
                min_distance = value;
                min_index = j;
            }
        }

        // 标记所有在容差范围内的点（处理边界点）
        for (j, &value) in distance.iter().enumerate() {
            if value < min_distance + tolerance {
                let gp = if j == min_index {
                    i
                } else {
                    boundary_num_gp
                        .checked_add(total_num_gp)
                        .ok_or(SymError::ArraySizeShortage)?
                };
                if gp >= bz_grid_address.len() {
                    return Err(SymError::ArraySizeShortage);
                }

                let mut bz_address_double = [0; 3];
                for k in 0..3 {
                    let coordinate = i64::from(address[k])
                        + i64::from(BZ_SEARCH_SPACE[j][k]) * i64::from(mesh[k]);
                    bz_grid_address[gp][k] =
                        i32::try_from(coordinate).map_err(|_| SymError::InvalidInput)?;
                    let doubled = coordinate * 2 + i64::from(is_shift[k]);
                    bz_address_double[k] =
                        i32::try_from(doubled).map_err(|_| SymError::InvalidInput)?;
                }

                let bzgp =
                    kgrid::get_dense_grid_point_double_mesh(&bz_address_double, &bzmesh)?;
                bz_map[bzgp] = gp;

                if j != min_index {
                    boundary_num_gp += 1;
                }
            }
        }
    }

    boundary_num_gp
        .checked_add(total_num_gp)
        .ok_or(SymError::ArraySizeShortage)
}

fn get_tolerance_for_bz_reduction(rec_lattice: &[[f64; 3]; 3], mesh: &[i32; 3]) -> f64 {
    let mut length = [0.0; 3];
    for (i, value) in length.iter_mut().enumerate() {
        *value = rec_lattice.iter().map(|row| row[i] * row[i]).sum::<f64>()
            / ((mesh[i] as f64) * (mesh[i] as f64));
    }

    let tolerance = length.into_iter().fold(f64::NEG_INFINITY, f64::max);
    tolerance * 0.01
}

/// 检查网格对称性
///
/// 注意：提供的 C 代码中此函数存在复制粘贴错误（重复检查 column 0）。
/// 本 Rust 实现已修正此问题，正确检查 column 0, 1, 2 分别对应 a=b, b=c, c=a 的对称性。
fn check_mesh_symmetry(mesh: &[i32; 3], is_shift: &[i32; 3], rot_reciprocal: &[Mat3I]) -> bool {
    let mut eq = [false; 3]; // eq[0]: a=b, eq[1]: b=c, eq[2]: c=a

    for rot in rot_reciprocal {
        if rot.iter().flatten().map(|x| x.abs()).sum::<i32>() > 3 {
            return false;
        }
    }

    for rot in rot_reciprocal {
        // 检查 x <-> y 交换 (a=b)
        // 矩阵应为 [0 1 0; 1 0 0; 0 0 1] 或类似，关注列 0 变为 [0, 1, 0]
        if rot[0][0] == 0 && rot[1][0] == 1 && rot[2][0] == 0 {
            eq[0] = true;
        }
        // 检查 y <-> z 交换 (b=c)
        // 关注列 1 变为 [0, 0, 1]
        if rot[0][1] == 0 && rot[1][1] == 0 && rot[2][1] == 1 {
            eq[1] = true;
        }
        // 检查 z <-> x 交换 (c=a)
        // 关注列 2 变为 [1, 0, 0]
        if rot[0][2] == 1 && rot[1][2] == 0 && rot[2][2] == 0 {
            eq[2] = true;
        }
    }

    let cond1 = (eq[0] && mesh[0] == mesh[1] && is_shift[0] == is_shift[1]) || !eq[0];
    let cond2 = (eq[1] && mesh[1] == mesh[2] && is_shift[1] == is_shift[2]) || !eq[1];
    let cond3 = (eq[2] && mesh[2] == mesh[0] && is_shift[2] == is_shift[0]) || !eq[2];

    cond1 && cond2 && cond3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_mesh_symmetry() {
        let mesh = [4, 4, 4];
        let shift = [0, 0, 0];
        let mut rot = vec![[[0; 3]; 3]; 1];
        // Identity
        rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        assert!(check_mesh_symmetry(&mesh, &shift, &rot));
    }

    #[test]
    fn test_get_irreducible_reciprocal_mesh_simple() {
        // 简单的 2x2x2 网格，无位移，只有恒等操作
        let mesh = [2, 2, 2];
        let shift = [0, 0, 0];
        let mut rot = vec![[[0; 3]; 3]; 1];
        rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        let mut grid_address = vec![[0; 3]; 8];
        let mut map = vec![0; 8];

        let num_ir =
            get_irreducible_reciprocal_mesh(&mut grid_address, &mut map, &mesh, &shift, &rot)
                .unwrap();

        // 由于只有恒等操作，每个点都是不可约的
        assert_eq!(num_ir, 8);
        for (i, &representative) in map.iter().enumerate() {
            assert_eq!(representative, i);
        }
    }

    #[test]
    fn mesh_buffers_are_checked_before_indexing() {
        let mesh = [2, 2, 2];
        let shift = [0, 0, 0];
        let mut rot = vec![[[0; 3]; 3]; 1];
        rot[0] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        let mut short_grid = vec![[0; 3]; 7];
        let mut map = vec![0; 8];
        assert!(matches!(
            get_irreducible_reciprocal_mesh(&mut short_grid, &mut map, &mesh, &shift, &rot,),
            Err(SymError::ArraySizeShortage)
        ));

        let mut no_points = [];
        assert!(matches!(
            get_dense_grid_points_by_rotations(&mut no_points, &[0, 0, 0], &rot, &mesh, &shift,),
            Err(SymError::ArraySizeShortage)
        ));
    }

    #[test]
    fn irreducible_mapping_is_fully_path_compressed() {
        let mut mapping = vec![0, 0, 1, 2, 4, 4];
        let count = get_dense_num_ir(&mut mapping, 6).unwrap();
        assert_eq!(count, 2);
        assert_eq!(mapping, vec![0, 0, 0, 0, 4, 4]);
    }
}
