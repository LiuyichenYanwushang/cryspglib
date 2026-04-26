//! Delaunay 晶格约化。
//!
//! 将任意晶格约化为 Delaunay 标准形式（Niggli 约化的推广）。
//! 核心算法基于 Selling 变换（International Table A），通过迭代
//! 确保扩展基两两之间的点积 ≤ 0（Delaunay 条件）。
//!
//! 同时支持三维周期性晶格和二维层状晶格（有非周期性轴）。

use crate::cell::AperiodicAxis;
use crate::debug;
use crate::mathfunc::{
    Mat3, Mat3I, Vec3, mat_cast_matrix_3d_to_3i,
    mat_get_determinant_d3, mat_get_determinant_i3, mat_inverse_matrix_d3, mat_multiply_matrix_d3,
    mat_norm_squared_d3,
};
use std::env;

const ZERO_PREC: f64 = 1e-10;

/// 获取最大尝试次数，从环境变量 SPGLIB_NUM_ATTEMPTS 读取，默认为 1000
fn get_num_attempts() -> i32 {
    if let Ok(val_str) = env::var("SPGLIB_NUM_ATTEMPTS") {
        if let Ok(val) = val_str.parse::<i32>() {
            if val > 0 {
                return val;
            }
        }
        debug::warning_print(format_args!(
            "spglib: Could not parse SPGLIB_NUM_ATTEMPTS={}\n",
            val_str
        ));
    }
    1000
}

/// Delaunay 约化主入口
/// 返回 true 表示成功，false 表示失败
/// min_lattice: 输出的约化后晶格
/// lattice: 输入晶格
/// symprec: 对称性判定精度
pub fn del_delaunay_reduce(lattice: &Mat3, symprec: f64) -> Option<Mat3> {
    debug::debug_print(format_args!(
        "del_delaunay_reduce (tolerance = {}):\n",
        symprec
    ));
    delaunay_reduce(lattice, -1, symprec)
}

/// 层状结构的 Delaunay 约化
/// aperiodic_axis: 非周期轴
pub fn del_layer_delaunay_reduce(
    lattice: &Mat3,
    aperiodic_axis: Option<AperiodicAxis>,
    symprec: f64,
) -> Option<Mat3> {
    debug::debug_print(format_args!(
        "del_layer_delaunay_reduce (tolerance = {}):\n",
        symprec
    ));
    let ap_i32 = aperiodic_axis.map_or(-1, |ap| ap.axis_index() as i32);
    delaunay_reduce(lattice, ap_i32, symprec)
}

/// Delaunay 约化核心逻辑
/// Reference: International table A.
fn delaunay_reduce(
    lattice: &Mat3,
    aperiodic_axis: i32,
    symprec: f64,
) -> Option<Mat3> {
    let mut orig_lattice = [[0.0; 3]; 3];
    orig_lattice = *lattice;

    // 扩展基 basis[4][3]
    let mut basis: [[f64; 3]; 4] = [[0.0; 3]; 4];
    let lattice_rank = get_extended_basis(&mut basis, lattice, aperiodic_axis);

    let max_attempt = get_num_attempts();
    let mut succeeded = false;

    for attempt in 0..max_attempt {
        debug::debug_print(format_args!(
            "Trying delaunay_reduce_basis: attempt {}/{}\n",
            attempt, max_attempt
        ));
        // delaunay_reduce_basis 返回 true 表示已收敛（无需修改），false 表示进行了修改需继续迭代
        if delaunay_reduce_basis(&mut basis, lattice_rank, symprec) {
            succeeded = true;
            break;
        }
    }

    if !succeeded {
        return None;
    }

    get_delaunay_shortest_vectors(&mut basis, lattice_rank, symprec);

    // 将 basis 的前三个向量复制回 red_lattice
    // 注意：C 代码中 basis[j][i] -> red_lattice[i][j]，即转置关系
    // basis 是行向量列表 (4x3)，red_lattice 是列向量矩阵 (3x3)
    let mut red_lattice = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            red_lattice[i][j] = basis[j][i];
        }
    }

    // 如果是层状结构 (rank=2) 且非周期轴不是 Z 轴 (2)，需要把非周期轴移回原来的方向
    if lattice_rank == 2 && aperiodic_axis != 2 {
        let axis = aperiodic_axis as usize;
        for i in 0..3 {
            // 交换第 axis 列和第 2 列
            let temp = red_lattice[i][axis];
            red_lattice[i][axis] = red_lattice[i][2]; // basis[2] 是之前暂存非周期轴的地方
            red_lattice[i][2] = temp;
        }
    }

    let volume = mat_get_determinant_d3(&red_lattice);
    if volume.abs() < symprec {
        debug::info_print(format_args!("spglib: Minimum lattice has no volume.\n"));
        return None;
    }

    // 保持体积为正（右手系）
    if volume < 0.0 {
        for i in 0..3 {
            for j in 0..3 {
                red_lattice[i][j] = -red_lattice[i][j];
            }
        }
    }

    // 验证变换矩阵是否为幺模矩阵 (行列式为 +/- 1)
    // M_new = M_trans * M_orig => M_trans = M_new * M_orig^-1
    if let Ok(inv_red) = mat_inverse_matrix_d3(&red_lattice, symprec) {
        let tmp_mat = mat_multiply_matrix_d3(&inv_red, &orig_lattice);
        let tmp_mat_int = mat_cast_matrix_3d_to_3i(&tmp_mat);
        if mat_get_determinant_i3(&tmp_mat_int).abs() != 1 {
            debug::info_print(format_args!(
                "spglib: Determinant of Delaunay change of basis matrix has to be 1 or -1.\n"
            ));
            return None;
        }
    } else {
        return None;
    }

    Some(red_lattice)
}

/// 从 Delaunay 集合中搜索最短向量
/// 如果存在非周期轴，则在周期平面内找最短的两个，在平面外找最短的一个
fn get_delaunay_shortest_vectors(basis: &mut [[f64; 3]; 4], lattice_rank: usize, symprec: f64) {
    // 构造 7 个候选向量集合: {b1, b2, b1+b2, b3, b4, b2+b3, b3+b1}
    let mut b: [[f64; 3]; 7] = [[0.0; 3]; 7];

    // b[0] = basis[0], b[1] = basis[1]
    b[0] = basis[0];
    b[1] = basis[1];

    // b[2] = basis[0] + basis[1]
    for i in 0..3 {
        b[2][i] = basis[0][i] + basis[1][i];
    }

    // b[3] = basis[2], b[4] = basis[3]
    b[3] = basis[2];
    b[4] = basis[3];

    // b[5] = basis[1] + basis[2]
    for i in 0..3 {
        b[5][i] = basis[1][i] + basis[2][i];
    }

    // b[6] = basis[2] + basis[0]
    for i in 0..3 {
        b[6][i] = basis[2][i] + basis[0][i];
    }

    // 冒泡排序
    if lattice_rank == 3 {
        // 对所有 7 个向量排序 (索引 0..6)
        for _ in 0..6 {
            for j in 0..6 {
                if mat_norm_squared_d3(&b[j]) > mat_norm_squared_d3(&b[j + 1]) + ZERO_PREC {
                    let tmp = b[j];
                    b[j] = b[j + 1];
                    b[j + 1] = tmp;
                }
            }
        }
    } else {
        // 排序前三个 (b1, b2, b1+b2)
        for _ in 0..2 {
            for j in 0..2 {
                if mat_norm_squared_d3(&b[j]) > mat_norm_squared_d3(&b[j + 1]) + ZERO_PREC {
                    let tmp = b[j];
                    b[j] = b[j + 1];
                    b[j + 1] = tmp;
                }
            }
        }
        // 排序后四个 (b3, b4, b2+b3, b3+b1)
        for _ in 3..6 {  // 外层循环 3 次，对应 C 的 for (i = 3; i <= 5; i++)
            for j in 3..6 {
                if mat_norm_squared_d3(&b[j]) > mat_norm_squared_d3(&b[j + 1]) + ZERO_PREC {
                    let tmp = b[j];
                    b[j] = b[j + 1];
                    b[j + 1] = tmp;
                }
            }
        }
    }

    // 选取线性无关的三个向量
    let mut tmpmat = [[0.0; 3]; 3];
    for i in 2..7 {
        // 构造矩阵 [b[0], b[1], b[i]] (列向量)
        for j in 0..3 {
            tmpmat[j][0] = b[0][j];
            tmpmat[j][1] = b[1][j];
            tmpmat[j][2] = b[i][j];
        }

        if mat_get_determinant_d3(&tmpmat).abs() > symprec {
            // 找到线性无关组，更新 basis
            basis[0] = b[0];
            basis[1] = b[1];
            basis[2] = b[i];
            return;
        }
    }
}


/// 执行一次 Delaunay 约化步骤
/// 返回 true 表示已经满足 Delaunay 条件（无需修改），false 表示进行了修改
fn delaunay_reduce_basis(basis: &mut [[f64; 3]; 4], lattice_rank: usize, symprec: f64) -> bool {
    for i in 0..3 {
        for j in (i + 1)..4 {
            let mut dot_product = 0.0;
            for k in 0..3 {
                dot_product += basis[i][k] * basis[j][k];
            }

            if dot_product > symprec {
                if i < lattice_rank {
                    // 执行 Selling 变换
                    // b_k = b_k + b_i for k != i, j
                    for k in 0..4 {
                        if k != i && k != j {
                            for l in 0..3 {
                                basis[k][l] += basis[i][l];
                            }
                        }
                    }
                    // b_i = -b_i
                    for k in 0..3 {
                        basis[i][k] = -basis[i][k];
                    }
                    return false; // 发生了修改，需要重新迭代
                } else {
                    // lattice_rank == 2 且 dot(b3, b4) > 0 的情况
                    // 仅打印警告，不进行修改
                    debug::info_print(format_args!(
                        "spglib: Dot product between basis {}, {} larger than 0.\n",
                        i + 1,
                        j + 1
                    ));
                    debug::debug_print_vectors_d3(basis);
                }
            }
        }
    }
    true // 未发生修改，已收敛
}

/// 获取扩展基 (Extended Basis)
/// 如果 aperiodic_axis != -1，则将其暂时移动到 basis[2]
/// 返回 lattice_rank
fn get_extended_basis(basis: &mut [[f64; 3]; 4], lattice: &Mat3, aperiodic_axis: i32) -> usize {
    let lattice_rank;

    if aperiodic_axis == -1 {
        lattice_rank = 3;
        for i in 0..3 {
            for j in 0..3 {
                // 转置复制：lattice 是列向量矩阵，basis 存储行向量
                basis[i][j] = lattice[j][i];
            }
        }
    } else {
        lattice_rank = 2; // C 代码逻辑：初始设为 0 然后自增，结果是 2
        let mut current_rank = 0;
        for i in 0..3 {
            if i != aperiodic_axis {
                for j in 0..3 {
                    basis[current_rank][j] = lattice[j][i as usize];
                }
                current_rank += 1;
            }
        }
        // 将非周期轴放在 basis[2] (即 lattice_rank 位置)
        for j in 0..3 {
            basis[lattice_rank][j] = lattice[j][aperiodic_axis as usize];
        }
    }

    // basis[3] = -(basis[0] + basis[1] + basis[2])
    for i in 0..3 {
        basis[3][i] = -lattice[i][0] - lattice[i][1] - lattice[i][2];
    }

    lattice_rank
}

// --- 2D Delaunay Reduction Functions ---

/// 2D 层状结构的 Delaunay 约化
pub fn del_layer_delaunay_reduce_2D(
    red_lattice: &mut Mat3,
    lattice: &Mat3,
    unique_axis: i32,
    aperiodic_axis: i32,
    symprec: f64,
) -> bool {
    debug::debug_print(format_args!("del_layer_delaunay_reduce_2D:\n"));

    let mut j = -1;
    let mut k = -1;
    let lattice_rank;

    if aperiodic_axis == -1 || unique_axis == aperiodic_axis {
        // bulk or Monoclinic/oblique
        for i in 0..3 {
            if i != unique_axis {
                if j == -1 {
                    j = i;
                } else {
                    k = i;
                }
            }
        }
        lattice_rank = 2;
    } else {
        // Monoclinic/rectangular
        for i in 0..3 {
            if i != unique_axis && i != aperiodic_axis {
                j = i;
            }
        }
        k = aperiodic_axis;
        lattice_rank = 1;
    }

    if j < 0 || k < 0 {
        // 逻辑错误防护
        return false;
    }
    let j = j as usize;
    let k = k as usize;
    let unique_axis = unique_axis as usize;

    let mut unique_vec = [0.0; 3];
    let mut lattice_2d = [[0.0; 2]; 3];

    for i in 0..3 {
        unique_vec[i] = lattice[i][unique_axis];
        lattice_2d[i][0] = lattice[i][j];
        lattice_2d[i][1] = lattice[i][k];
    }

    let mut basis: [[f64; 3]; 3] = [[0.0; 3]; 3];
    get_extended_basis_2d(&mut basis, &lattice_2d);

    let max_attempt = get_num_attempts();
    let mut succeeded = false;

    for attempt in 0..max_attempt {
        debug::debug_print(format_args!(
            "Trying delaunay_reduce_basis_2D: attempt {}/{}\n",
            attempt, max_attempt
        ));
        if delaunay_reduce_basis_2d(&mut basis, lattice_rank, symprec) {
            succeeded = true;
            break;
        }
    }

    if !succeeded {
        return false;
    }

    get_delaunay_shortest_vectors_2d(&mut basis, &unique_vec, lattice_rank, symprec);

    for i in 0..3 {
        red_lattice[i][unique_axis] = lattice[i][unique_axis];
        red_lattice[i][j] = basis[0][i];
        red_lattice[i][k] = basis[1][i];
    }

    let volume = mat_get_determinant_d3(red_lattice);
    if volume.abs() < symprec {
        debug::info_print(format_args!("spglib: Minimum lattice has no volume.\n"));
        return false;
    }

    if volume < 0.0 {
        for i in 0..3 {
            red_lattice[i][unique_axis] = -red_lattice[i][unique_axis];
        }
    }

    true
}

fn delaunay_reduce_basis_2d(basis: &mut [[f64; 3]; 3], lattice_rank: usize, symprec: f64) -> bool {
    for i in 0..2 {
        for j in (i + 1)..3 {
            let mut dot_product = 0.0;
            for k in 0..3 {
                dot_product += basis[i][k] * basis[j][k];
            }

            if dot_product > symprec {
                if i < lattice_rank {
                    for k in 0..3 {
                        if k != i && k != j {
                            for l in 0..3 {
                                basis[k][l] += 2.0 * basis[i][l];
                            }
                            break; // 2D 只有 3 个向量，找到第三个即可
                        }
                    }
                    for k in 0..3 {
                        basis[i][k] = -basis[i][k];
                    }
                    return false;
                } else {
                    debug::info_print(format_args!(
                        "spglib: Dot product between basis {}, {} larger than 0.\n",
                        i + 1,
                        j + 1
                    ));
                    debug::debug_print_vectors_d3(basis);
                }
            }
        }
    }
    true
}

fn get_delaunay_shortest_vectors_2d(
    basis: &mut [[f64; 3]; 3],
    unique_vec: &[f64; 3],
    lattice_rank: usize,
    symprec: f64,
) {
    let mut b: [[f64; 3]; 4] = [[0.0; 3]; 4];

    // b[0..3] = basis[0..3]
    for i in 0..3 {
        b[i] = basis[i];
    }
    // b[3] = basis[0] + basis[1]
    for i in 0..3 {
        b[3][i] = basis[0][i] + basis[1][i];
    }

    // Bubble sort
    // 范围从 lattice_rank % 2 开始到 3
    let start_idx = lattice_rank % 2;
    for _ in start_idx..3 {
        for j in start_idx..3 {
            if mat_norm_squared_d3(&b[j]) > mat_norm_squared_d3(&b[j + 1]) + ZERO_PREC {
                let tmp = b[j];
                b[j] = b[j + 1];
                b[j + 1] = tmp;
            }
        }
    }

    let mut tmpmat = [[0.0; 3]; 3];
    // 设置 tmpmat 前两列
    for i in 0..3 {
        tmpmat[i][0] = b[0][i];
        tmpmat[i][1] = unique_vec[i];
    }

    // 寻找第三列
    for i in 1..4 {
        for j in 0..3 {
            tmpmat[j][2] = b[i][j];
        }
        if mat_get_determinant_d3(&tmpmat).abs() > symprec {
            basis[0] = b[0];
            basis[1] = b[i];
            break;
        }
    }
}

fn get_extended_basis_2d(basis: &mut [[f64; 3]; 3], lattice_2d: &[[f64; 2]; 3]) {
    for i in 0..2 {
        for j in 0..3 {
            basis[i][j] = lattice_2d[j][i];
        }
    }
    for i in 0..3 {
        basis[2][i] = -lattice_2d[i][0] - lattice_2d[i][1];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delaunay_reduce_simple() {
        // 构造一个简单的正交晶格，不需要约化
        let lattice: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let min_lattice = del_delaunay_reduce(&lattice, 1e-5);
        assert!(min_lattice.is_some());
        let min_lattice = min_lattice.unwrap();
        // 结果应该保持不变（或仅仅是基向量顺序变化，对于单位阵应该不变）
        // 注意：Delaunay 约化可能会改变基向量的符号或顺序，但对于单位阵，
        // 扩展基是 (1,0,0), (0,1,0), (0,0,1), (-1,-1,-1)。
        // 它们之间点积为 0 (前三个) 或负数。所以应该稳定。
        assert!((min_lattice[0][0] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_delaunay_reduce_skewed() {
        // 构造一个倾斜的晶格
        // a = (1, 0, 0)
        // b = (1, 1, 0) -> 应该被约化为 (0, 1, 0)
        // c = (0, 0, 1)
        let lattice: Mat3 = [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let min_lattice = del_delaunay_reduce(&lattice, 1e-5);
        assert!(min_lattice.is_some());
        let min_lattice = min_lattice.unwrap();

        // 检查结果是否更正交
        // 期望结果可能是 [[1, 0, 0], [0, 1, 0], [0, 0, 1]] 或者类似的等价形式
        let dot_prod = min_lattice[0][0] * min_lattice[0][1]
            + min_lattice[1][0] * min_lattice[1][1]
            + min_lattice[2][0] * min_lattice[2][1];
        assert!(dot_prod.abs() < 1e-5);
    }

    #[test]
    fn test_get_extended_basis() {
        let lattice: Mat3 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let mut basis: [[f64; 3]; 4] = [[0.0; 3]; 4];
        get_extended_basis(&mut basis, &lattice, -1);

        // basis[0] should be (1, 0, 0)
        assert_eq!(basis[0], [1.0, 0.0, 0.0]);
        // basis[3] should be (-1, -2, -3)
        assert_eq!(basis[3], [-1.0, -2.0, -3.0]);
    }
}
