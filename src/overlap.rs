//! 原子位置重叠检测。
//!
//! 提供高效的原子重叠判断，用于在对称性检测中确认
//! 两个原子位置是否在给定精度下等价。

use crate::cell::{AperiodicAxis, Cell};
use crate::mathfunc::{
    Mat3, Vec3, mat_multiply_matrix_vector_d3, mat_multiply_matrix_vector_id3,
    mat_nint, mat_norm_squared_d3,
};
use std::cmp::Ordering;

/// 用于排序的辅助结构体
#[derive(Clone, Copy)]
struct ValueWithIndex {
    value: f64,
    type_: i32,
    index: usize,
}

impl Default for ValueWithIndex {
    fn default() -> Self {
        Self {
            value: 0.0,
            type_: 0,
            index: 0,
        }
    }
}

/// 重叠检查器结构体
/// 对应 C 代码中的 OverlapChecker
pub struct OverlapChecker {
    pub size: usize,
    pub lattice: Mat3,
    pub types_sorted: Vec<i32>,
    pub pos_sorted: Vec<Vec3>,     // 已排序的原子位置
    pub periodic_axes: [usize; 2], // 用于层状结构
    
    // 缓存区，避免重复分配内存 (对应 C 中的 blob 和 argsort_work)
    perm_temp: Vec<usize>,
    distance_temp: Vec<f64>,
    pos_temp_1: Vec<Vec3>,
    pos_temp_2: Vec<Vec3>,
    argsort_work: Vec<ValueWithIndex>, // 优化：复用排序工作区
}

impl OverlapChecker {
    /// 初始化 OverlapChecker
    /// 对应 C: ovl_overlap_checker_init
    pub fn new(cell: &Cell) -> Option<Self> {
        let size = cell.size;
        // 预分配所有 Vec，避免后续 realloc
        let mut checker = OverlapChecker {
            size,
            lattice: [[0.0; 3]; 3],
            types_sorted: vec![0; size],
            pos_sorted: vec![[0.0; 3]; size],
            periodic_axes: [0, 0],
            perm_temp: vec![0; size],
            distance_temp: vec![0.0; size],
            pos_temp_1: vec![[0.0; 3]; size],
            pos_temp_2: vec![[0.0; 3]; size],
            argsort_work: vec![ValueWithIndex::default(); size],
        };

        checker.lattice = cell.lattice;

        // 获取排序排列
        // 注意：这里直接使用 checker 内部的 buffer
        if !argsort_by_lattice_point_distance(
            &mut checker.perm_temp,
            &cell.lattice,
            &cell.position,
            Some(&cell.types),
            &mut checker.distance_temp,
            &mut checker.argsort_work,
        ) {
            return None;
        }

        // 应用排序
        checker.pos_sorted = permute_vec3(&cell.position, &checker.perm_temp);
        checker.types_sorted = permute_i32(cell.types.as_slice(), &checker.perm_temp);

        // 设置周期性轴 (用于层状结构)
        let mut lattice_rank = 0;
        for i in 0..3 {
            let is_periodic = match cell.aperiodic_axis {
                None => true,
                Some(ap) => i != ap.axis_index(),
            };
            if is_periodic {
                if lattice_rank < 2 {
                    checker.periodic_axes[lattice_rank] = i;
                }
                lattice_rank += 1;
            }
        }

        Some(checker)
    }

    /// 检查完全重叠
    /// 返回: -1 (Error), 0 (False), 1 (True)
    /// 对应 C: ovl_check_total_overlap
    pub fn check_total_overlap(
        &mut self,
        test_trans: &Vec3,
        rot: &[[i32; 3]; 3],
        symprec: f64,
        is_identity: bool,
    ) -> i32 {
        // 快速检查
        if !self.check_possible_overlap(test_trans, rot, symprec) {
            return 0;
        }

        // 计算旋转和平移后的位置
        for i in 0..self.size {
            if is_identity {
                self.pos_temp_1[i] = self.pos_sorted[i];
            } else {
                self.pos_temp_1[i] = mat_multiply_matrix_vector_id3(rot, &self.pos_sorted[i]);
            }

            for k in 0..3 {
                self.pos_temp_1[i][k] += test_trans[k];
            }
        }

        // 对变换后的位置进行排序
        if !argsort_by_lattice_point_distance(
            &mut self.perm_temp,
            &self.lattice,
            &self.pos_temp_1,
            Some(&self.types_sorted),
            &mut self.distance_temp,
            &mut self.argsort_work,
        ) {
            return -1;
        }

        self.pos_temp_2 = permute_vec3(&self.pos_temp_1, &self.perm_temp);

        // 检查排序后的重叠
        check_total_overlap_for_sorted(
            &self.lattice,
            &self.pos_sorted,
            &self.pos_temp_2,
            &self.types_sorted,
            &self.types_sorted,
            self.size,
            symprec,
        )
    }

    /// 检查层状结构的完全重叠
    /// 对应 C: ovl_check_layer_total_overlap
    pub fn check_layer_total_overlap(
        &mut self,
        test_trans: &Vec3,
        rot: &[[i32; 3]; 3],
        symprec: f64,
        is_identity: bool,
    ) -> i32 {
        if !self.check_possible_overlap(test_trans, rot, symprec) {
            return 0;
        }

        for i in 0..self.size {
            if is_identity {
                self.pos_temp_1[i] = self.pos_sorted[i];
            } else {
                self.pos_temp_1[i] = mat_multiply_matrix_vector_id3(rot, &self.pos_sorted[i]);
            }
            for k in 0..3 {
                self.pos_temp_1[i][k] += test_trans[k];
            }
        }

        if !argsort_by_lattice_point_distance(
            &mut self.perm_temp,
            &self.lattice,
            &self.pos_temp_1,
            Some(&self.types_sorted),
            &mut self.distance_temp,
            &mut self.argsort_work,
        ) {
            return -1;
        }

        self.pos_temp_2 = permute_vec3(&self.pos_temp_1, &self.perm_temp);

        check_layer_total_overlap_for_sorted(
            &self.lattice,
            &self.pos_sorted,
            &self.pos_temp_2,
            &self.types_sorted,
            &self.types_sorted,
            self.size,
            &self.periodic_axes,
            symprec,
        )
    }

    /// 快速预检查
    /// 对应 C: check_possible_overlap
    fn check_possible_overlap(&self, test_trans: &Vec3, rot: &[[i32; 3]; 3], symprec: f64) -> bool {
        let max_search_num = 3;
        let search_num = if self.size <= max_search_num {
            self.size
        } else {
            max_search_num
        };

        for i_test in 0..search_num {
            let type_rot = self.types_sorted[i_test];
            let mut pos_rot = mat_multiply_matrix_vector_id3(rot, &self.pos_sorted[i_test]);
            for k in 0..3 {
                pos_rot[k] += test_trans[k];
            }

            let mut is_found = false;
            // 暴力搜索，因为只检查前几个原子，开销可控
            for i in 0..self.size {
                if has_overlap_with_same_type(
                    &pos_rot,
                    &self.pos_sorted[i],
                    type_rot,
                    self.types_sorted[i],
                    &self.lattice,
                    symprec,
                ) {
                    is_found = true;
                    break;
                }
            }

            if !is_found {
                return false;
            }
        }
        true
    }
}

// --- Helper Functions ---

#[inline]
fn cartesian_norm(lat: &Mat3, v: &Vec3) -> f64 {
    let mut temp = [0.0; 3];
    for i in 0..3 {
        temp[i] = lat[i][0] * v[0] + lat[i][1] * v[1] + lat[i][2] * v[2];
    }
    (temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]).sqrt()
}

#[inline]
fn has_overlap(a: &Vec3, b: &Vec3, lattice: &Mat3, symprec: f64) -> bool {
    let mut v_diff = [0.0; 3];
    for i in 0..3 {
        v_diff[i] = a[i] - b[i];
        v_diff[i] -= mat_nint(v_diff[i]) as f64;
    }
    cartesian_norm(lattice, &v_diff) <= symprec
}

#[inline]
fn has_overlap_with_same_type(
    a: &Vec3,
    b: &Vec3,
    type_a: i32,
    type_b: i32,
    lattice: &Mat3,
    symprec: f64,
) -> bool {
    if type_a == type_b {
        has_overlap(a, b, lattice, symprec)
    } else {
        false
    }
}

#[inline]
fn layer_has_overlap(
    a: &Vec3,
    b: &Vec3,
    lattice: &Mat3,
    periodic_axes: &[usize; 2],
    symprec: f64,
) -> bool {
    let mut v_diff = [0.0; 3];
    for i in 0..3 {
        v_diff[i] = a[i] - b[i];
    }
    // 仅在周期性方向上应用最小镜像约定
    v_diff[periodic_axes[0]] -= mat_nint(v_diff[periodic_axes[0]]) as f64;
    v_diff[periodic_axes[1]] -= mat_nint(v_diff[periodic_axes[1]]) as f64;

    cartesian_norm(lattice, &v_diff) <= symprec
}

#[inline]
fn layer_has_overlap_with_same_type(
    a: &Vec3,
    b: &Vec3,
    type_a: i32,
    type_b: i32,
    lattice: &Mat3,
    periodic_axes: &[usize; 2],
    symprec: f64,
) -> bool {
    if type_a == type_b {
        layer_has_overlap(a, b, lattice, periodic_axes, symprec)
    } else {
        false
    }
}

// --- Sorting Logic ---

/// 根据到最近格点的距离对原子进行排序
/// 对应 C: argsort_by_lattice_point_distance
fn argsort_by_lattice_point_distance(
    perm: &mut [usize],
    lattice: &Mat3,
    positions: &[Vec3],
    types: Option<&[i32]>,
    distance_temp: &mut [f64],
    work: &mut [ValueWithIndex], // 传入预分配的工作区
) -> bool {
    let size = positions.len();
    let mut diff = [0.0; 3];

    for i in 0..size {
        for k in 0..3 {
            let x = positions[i][k];
            diff[k] = x - mat_nint(x) as f64;
        }
        let diff_cart = mat_multiply_matrix_vector_d3(lattice, &diff);
        distance_temp[i] = mat_norm_squared_d3(&diff_cart);
    }

    // 填充工作区
    for i in 0..size {
        work[i].value = distance_temp[i];
        work[i].type_ = if let Some(t) = types { t[i] } else { 0 };
        work[i].index = i;
    }

    // 排序：先按类型降序，再按距离升序
    // 对应 C: ValueWithIndex_comparator
    work.sort_by(|a, b| {
        let type_cmp = b.type_.cmp(&a.type_); // Descending
        if type_cmp != Ordering::Equal {
            type_cmp
        } else {
            a.value.partial_cmp(&b.value).unwrap_or(Ordering::Equal) // Ascending
        }
    });

    for i in 0..size {
        perm[i] = work[i].index;
    }

    true
}

/// 按置换表重排 `Vec3` 切片，返回重排后的向量。
fn permute_vec3(data_in: &[Vec3], perm: &[usize]) -> Vec<Vec3> {
    perm.iter().map(|&idx| data_in[idx]).collect()
}

/// 按置换表重排 `i32` 切片，返回重排后的向量。
fn permute_i32(data_in: &[i32], perm: &[usize]) -> Vec<i32> {
    perm.iter().map(|&idx| data_in[idx]).collect()
}

/// 检查两个已排序的原子列表是否重叠
/// 对应 C: check_total_overlap_for_sorted
fn check_total_overlap_for_sorted(
    lattice: &Mat3,
    pos_original: &[Vec3],
    pos_rotated: &[Vec3],
    types_original: &[i32],
    types_rotated: &[i32],
    num_pos: usize,
    symprec: f64,
) -> i32 {
    let mut found = vec![false; num_pos];
    let mut search_start = 0;

    for i_orig in 0..num_pos {
        // 跳过开头已经匹配过的原子
        while search_start < num_pos && found[search_start] {
            search_start += 1;
        }

        let mut matched = false;
        for i_rot in search_start..num_pos {
            if found[i_rot] {
                continue;
            }

            if has_overlap_with_same_type(
                &pos_original[i_orig],
                &pos_rotated[i_rot],
                types_original[i_orig],
                types_rotated[i_rot],
                lattice,
                symprec,
            ) {
                found[i_rot] = true;
                matched = true;
                break;
            }
        }

        if !matched {
            return 0;
        }
    }
    1
}

/// 检查两个已排序的原子列表是否重叠 (层状结构)
/// 对应 C: check_layer_total_overlap_for_sorted
fn check_layer_total_overlap_for_sorted(
    lattice: &Mat3,
    pos_original: &[Vec3],
    pos_rotated: &[Vec3],
    types_original: &[i32],
    types_rotated: &[i32],
    num_pos: usize,
    periodic_axes: &[usize; 2],
    symprec: f64,
) -> i32 {
    let mut found = vec![false; num_pos];
    let mut search_start = 0;

    for i_orig in 0..num_pos {
        while search_start < num_pos && found[search_start] {
            search_start += 1;
        }

        let mut matched = false;
        for i_rot in search_start..num_pos {
            if found[i_rot] {
                continue;
            }

            if layer_has_overlap_with_same_type(
                &pos_original[i_orig],
                &pos_rotated[i_rot],
                types_original[i_orig],
                types_rotated[i_rot],
                lattice,
                periodic_axes,
                symprec,
            ) {
                found[i_rot] = true;
                matched = true;
                break;
            }
        }

        if !matched {
            return 0;
        }
    }
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{Cell, TensorRank};

    /// 创建一个简单的立方晶胞，包含一个原子
    fn simple_cell() -> Cell {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0]];
        let types = [1];
        let mut cell = Cell::new(1, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);
        cell
    }

    /// 创建包含两个原子的晶胞，用于测试纯平移
    fn dimer_cell() -> Cell {
        let lattice = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let types = [1, 1];
        let mut cell = Cell::new(2, TensorRank::NoSpin);
        cell.set_cell(&lattice, &positions, &types);
        cell
    }

    /// 创建层状结构 (c 轴非周期)
    fn layer_cell() -> Cell {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.1, 0.2, 0.3]];
        let types = [1];
        let mut cell = Cell::new(1, TensorRank::NoSpin);
        cell.set_layer_cell(&lattice, &positions, &types, Some(AperiodicAxis::Z));
        cell
    }

    #[test]
    fn test_has_overlap() {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let a = [0.1, 0.1, 0.1];
        let b = [0.100000001, 0.1, 0.1];
        let c = [0.5, 0.5, 0.5];

        assert!(has_overlap(&a, &b, &lattice, 1e-5));
        assert!(!has_overlap(&a, &c, &lattice, 1e-5));

        // 周期性边界
        let d = [1.1, 0.1, 0.1];
        assert!(has_overlap(&a, &d, &lattice, 1e-5));
    }

    #[test]
    fn test_has_overlap_with_same_type() {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let a = [0.1, 0.1, 0.1];
        let b = [0.100000001, 0.1, 0.1];
        let a_type = 1;
        let b_type = 1;
        let c_type = 2;

        assert!(has_overlap_with_same_type(&a, &b, a_type, b_type, &lattice, 1e-5));
        assert!(!has_overlap_with_same_type(&a, &b, a_type, c_type, &lattice, 1e-5));
    }

    #[test]
    fn test_layer_has_overlap() {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let periodic_axes = [0, 1]; // x, y 周期性
        let a = [0.1, 0.2, 0.3];
        let b = [1.1, 0.2, 0.3]; // 在 x 方向平移一个周期
        let c = [0.1, 1.2, 0.3]; // 在 y 方向平移一个周期
        let d = [0.1, 0.2, 1.3]; // 在 z 方向平移（非周期），不应重叠

        assert!(layer_has_overlap(&a, &b, &lattice, &periodic_axes, 1e-5));
        assert!(layer_has_overlap(&a, &c, &lattice, &periodic_axes, 1e-5));
        assert!(!layer_has_overlap(&a, &d, &lattice, &periodic_axes, 1e-5));
    }

    #[test]
    fn test_overlap_checker_simple() {
        let cell = simple_cell();
        let mut checker = OverlapChecker::new(&cell).expect("Failed to create checker");

        let identity_rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let zero_trans = [0.0, 0.0, 0.0];

        // 恒等操作应返回 1 (True)
        let result = checker.check_total_overlap(&zero_trans, &identity_rot, 1e-5, true);
        assert_eq!(result, 1);

        // 一个非对称平移应返回 0
        let bad_trans = [0.5, 0.5, 0.5];
        let result = checker.check_total_overlap(&bad_trans, &identity_rot, 1e-5, true);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_overlap_checker_dimer() {
        let cell = dimer_cell();
        let mut checker = OverlapChecker::new(&cell).expect("Failed to create checker");

        // 纯平移 (0.5, 0, 0) 应该交换两个原子，是一个对称操作
        let identity_rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let half_trans = [0.5, 0.0, 0.0];
        let result = checker.check_total_overlap(&half_trans, &identity_rot, 1e-5, true);
        assert_eq!(result, 1);

        // 恒等平移也应是对称的
        let zero_trans = [0.0, 0.0, 0.0];
        let result = checker.check_total_overlap(&zero_trans, &identity_rot, 1e-5, true);
        assert_eq!(result, 1);
    }

    #[test]
    fn test_overlap_checker_layer() {
        let cell = layer_cell();
        let mut checker = OverlapChecker::new(&cell).expect("Failed to create checker");

        // 层状结构中，沿 aperiodic_axis (z) 的平移不应是对称的
        let identity_rot = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let z_trans = [0.0, 0.0, 0.5];
        let result = checker.check_layer_total_overlap(&z_trans, &identity_rot, 1e-5, true);
        assert_eq!(result, 0); // 预期不是对称操作

        // 沿 x 方向平移整数倍应是周期性的
        let x_trans = [1.0, 0.0, 0.0];
        let result = checker.check_layer_total_overlap(&x_trans, &identity_rot, 1e-5, true);
        assert_eq!(result, 1);
    }
}
