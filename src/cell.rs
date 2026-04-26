//! 晶胞结构体和相关操作。
//!
//! [`Cell`] 是晶体结构的核心表示，包含晶格矢量、原子位置和类型。
//! 晶格矩阵采用 [`Mat3`] 的 `lattice[cart][vec]` 布局——
//! 行=笛卡尔分量 (x,y,z)，列=晶格矢量 (a,b,c)。
//! 详见 [`mathfunc`](crate::mathfunc) 模块文档。

use crate::debug;
use crate::mathfunc::{
    Mat3, Vec3, VecDBL, mat_cast_matrix_3d_to_3i, mat_dabs,
    mat_dmod1, mat_get_determinant_d3, mat_get_determinant_i3,
    mat_inverse_matrix_d3, mat_multiply_matrix_d3, mat_multiply_matrix_vector_d3,
    mat_multiply_matrix_vector_id3, mat_nint, mat_norm_squared_d3,
};

const INCREASE_RATE: f64 = 2.0;
const REDUCE_RATE: f64 = 0.95;
const NUM_ATTEMPT: i32 = 100;

/// 磁性张量类型
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum TensorRank {
    NoSpin = -1, // 对应 C 中的 NOSPIN
    Collinear = 0,
    NonCollinear = 1,
}

/// 层状结构的非周期轴。
///
/// `None` 表示三维完全周期性，`Some(AperiodicAxis)` 表示该方向无周期。
/// 例如 `Some(X)` 表示只有 YZ 面内周期，X 方向非周期。
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AperiodicAxis {
    /// X 轴非周期；YZ 面内周期
    X,
    /// Y 轴非周期；ZX 面内周期
    Y,
    /// Z 轴非周期；XY 面内周期
    Z,
}

impl AperiodicAxis {
    /// 返回两个周期轴的索引。
    pub fn periodic_indices(self) -> [usize; 2] {
        match self {
            AperiodicAxis::X => [1, 2],
            AperiodicAxis::Y => [0, 2],
            AperiodicAxis::Z => [0, 1],
        }
    }

    /// 格子秩（层状结构始终为 2）。
    pub fn lattice_rank(self) -> usize {
        2
    }

    /// 非周期轴的索引 (0, 1, 2)。
    pub fn axis_index(self) -> usize {
        match self {
            AperiodicAxis::X => 0,
            AperiodicAxis::Y => 1,
            AperiodicAxis::Z => 2,
        }
    }
}

/// 晶胞结构体。
///
/// `lattice` 采用 `[cart][vec]` 布局：`lattice[i][j]` = 第 i 个笛卡尔分量在第 j 个晶格矢量上的投影。
/// `position` 以分数坐标存储（相对于晶格矢量）。
///
/// 对应 C 中的 `Cell` 结构体。
#[derive(Clone, Debug)]
pub struct Cell {
    /// 原子数量
    pub size: usize,
    /// 3x3 晶格矩阵，布局 `[cart][vec]`（行=笛卡尔, 列=晶格矢量）
    pub lattice: Mat3,
    /// 原子类型（原子序数）
    pub types: Vec<i32>,
    /// 原子位置（分数坐标，相对于 lattice 矢量）
    pub position: Vec<Vec3>,
    /// 磁性张量（扁平化存储：Collinear 为 size，NonCollinear 为 size*3）
    pub tensors: Vec<f64>,
    /// 磁性张量秩
    pub tensor_rank: TensorRank,
    /// 非周期轴：None=三维周期, Some(X/Y/Z)=对应轴非周期
    pub aperiodic_axis: Option<AperiodicAxis>,
}

impl Cell {
    /// 创建一个新的 Cell
    /// 对应 C: cel_alloc_cell
    pub fn new(size: usize, tensor_rank: TensorRank) -> Self {
        let tensors_size = match tensor_rank {
            TensorRank::NoSpin => 0,
            TensorRank::Collinear => size,
            TensorRank::NonCollinear => size * 3,
        };

        Cell {
            size,
            lattice: [[0.0; 3]; 3],
            types: vec![0; size],
            position: vec![[0.0; 3]; size],
            tensors: vec![0.0; tensors_size],
            tensor_rank,
            aperiodic_axis: None,
        }
    }

    /// 设置 Cell 数据
    /// 对应 C: cel_set_cell
    pub fn set_cell(&mut self, lattice: &Mat3, position: &[Vec3], types: &[i32]) {
        self.lattice = *lattice;
        for i in 0..self.size {
            for j in 0..3 {
                // 确保位置在 [-0.5, 0.5) 区间内，或者 [0, 1) 取决于 mat_nint 实现
                // C 代码逻辑：position - Nint(position)
                self.position[i][j] = position[i][j] - mat_nint(position[i][j]) as f64;
            }
            self.types[i] = types[i];
        }
    }

    /// 设置层状 Cell 数据
    /// 对应 C: cel_set_layer_cell
    pub fn set_layer_cell(
        &mut self,
        lattice: &Mat3,
        position: &[Vec3],
        types: &[i32],
        aperiodic_axis: Option<AperiodicAxis>,
    ) {
        self.lattice = *lattice;
        for i in 0..self.size {
            for j in 0..3 {
                if aperiodic_axis.map_or(true, |ap| j != ap.axis_index()) {
                    self.position[i][j] = position[i][j] - mat_nint(position[i][j]) as f64;
                } else {
                    self.position[i][j] = position[i][j];
                }
            }
            self.types[i] = types[i];
        }
        self.aperiodic_axis = aperiodic_axis;
    }

    /// 设置带张量的 Cell 数据
    /// 对应 C: cel_set_cell_with_tensors
    pub fn set_cell_with_tensors(
        &mut self,
        lattice: &Mat3,
        position: &[Vec3],
        types: &[i32],
        tensors: &[f64],
    ) {
        self.set_cell(lattice, position, types);
        match self.tensor_rank {
            TensorRank::Collinear => {
                self.tensors[..self.size].copy_from_slice(&tensors[..self.size]);
            }
            TensorRank::NonCollinear => {
                self.tensors[..self.size * 3].copy_from_slice(&tensors[..self.size * 3]);
            }
            _ => {}
        }
    }
}

/// 检查两个原子位置是否重叠（考虑周期性边界）
/// 对应 C: cel_is_overlap
pub fn cel_is_overlap(a: &Vec3, b: &Vec3, lattice: &Mat3, symprec: f64) -> bool {
    let mut v_diff = [0.0; 3];
    for i in 0..3 {
        v_diff[i] = a[i] - b[i];
        v_diff[i] -= mat_nint(v_diff[i]) as f64;
    }

    let v_diff = mat_multiply_matrix_vector_d3(lattice, &v_diff);
    mat_norm_squared_d3(&v_diff).sqrt() < symprec
}

/// 检查两个相同类型的原子是否重叠
/// 对应 C: cel_is_overlap_with_same_type
pub fn cel_is_overlap_with_same_type(
    a: &Vec3,
    b: &Vec3,
    type_a: i32,
    type_b: i32,
    lattice: &Mat3,
    symprec: f64,
) -> bool {
    if type_a == type_b {
        cel_is_overlap(a, b, lattice, symprec)
    } else {
        false
    }
}

/// 检查 Cell 中是否存在任何重叠原子
/// 对应 C: cel_any_overlap
pub fn cel_any_overlap(cell: &Cell, symprec: f64) -> bool {
    // 使用 rayon 并行化可能在这里收益不大，因为通常 cell.size 较小，且有早期返回
    for i in 0..cell.size {
        for j in (i + 1)..cell.size {
            if cel_is_overlap(&cell.position[i], &cell.position[j], &cell.lattice, symprec) {
                return true;
            }
        }
    }
    false
}

/// 检查 Cell 中是否存在任何相同类型的重叠原子
/// 对应 C: cel_any_overlap_with_same_type
pub fn cel_any_overlap_with_same_type(cell: &Cell, symprec: f64) -> bool {
    for i in 0..cell.size {
        for j in (i + 1)..cell.size {
            if cel_is_overlap_with_same_type(
                &cell.position[i],
                &cell.position[j],
                cell.types[i],
                cell.types[j],
                &cell.lattice,
                symprec,
            ) {
                return true;
            }
        }
    }
    false
}

/// 层状结构的重叠检查（仅在两个周期性方向上应用周期性边界）
/// 对应 C: cel_layer_is_overlap
pub fn cel_layer_is_overlap(
    a: &Vec3,
    b: &Vec3,
    lattice: &Mat3,
    aperiodic: AperiodicAxis,
    symprec: f64,
) -> bool {
    let pa = aperiodic.periodic_indices();
    let mut v_diff = [0.0; 3];
    v_diff[0] = a[0] - b[0];
    v_diff[1] = a[1] - b[1];
    v_diff[2] = a[2] - b[2];

    v_diff[pa[0]] -= mat_nint(v_diff[pa[0]]) as f64;
    v_diff[pa[1]] -= mat_nint(v_diff[pa[1]]) as f64;

    let v_diff = mat_multiply_matrix_vector_d3(lattice, &v_diff);
    mat_norm_squared_d3(&v_diff).sqrt() < symprec
}

/// 对应 C: cel_layer_is_overlap_with_same_type
pub fn cel_layer_is_overlap_with_same_type(
    a: &Vec3,
    b: &Vec3,
    type_a: i32,
    type_b: i32,
    lattice: &Mat3,
    aperiodic: AperiodicAxis,
    symprec: f64,
) -> bool {
    if type_a == type_b {
        cel_layer_is_overlap(a, b, lattice, aperiodic, symprec)
    } else {
        false
    }
}

/// 对应 C: cel_layer_any_overlap_with_same_type
pub fn cel_layer_any_overlap_with_same_type(
    cell: &Cell,
    aperiodic: AperiodicAxis,
    symprec: f64,
) -> bool {
    for i in 0..cell.size {
        for j in (i + 1)..cell.size {
            if cel_layer_is_overlap_with_same_type(
                &cell.position[i],
                &cell.position[j],
                cell.types[i],
                cell.types[j],
                &cell.lattice,
                aperiodic,
                symprec,
            ) {
                return true;
            }
        }
    }
    false
}

/// 裁剪晶胞 (Trim Cell)
/// 将较大的晶胞 (cell) 投影到较小的晶格 (trimmed_lattice) 上。
/// mapping_table: 输出参数，映射原原子索引到新原子索引
/// 对应 C: cel_trim_cell
pub fn cel_trim_cell(
    mapping_table: &mut [usize],
    trimmed_lattice: &Mat3,
    cell: &Cell,
    symprec: f64,
) -> Option<Cell> {
    trim_cell(mapping_table, trimmed_lattice, cell, symprec)
}

/// 内部实现：trim_cell
/// 对应 C: static trim_cell
fn trim_cell(
    mapping_table: &mut [usize],
    trimmed_lattice: &Mat3,
    cell: &Cell,
    symprec: f64,
) -> Option<Cell> {
    debug::debug_print(format_args!("trim_cell\n"));

    let det_cell = mat_get_determinant_d3(&cell.lattice);
    let det_trimmed = mat_get_determinant_d3(trimmed_lattice);
    // 计算体积比，必须是整数
    let ratio = mat_nint(det_cell / det_trimmed).abs();

    // 验证变换矩阵的行列式是否等于体积比
    let inv_trimmed = mat_inverse_matrix_d3(trimmed_lattice, symprec).ok()?;
    let tmp_mat = mat_multiply_matrix_d3(&inv_trimmed, &cell.lattice);
    let tmp_mat_int = mat_cast_matrix_3d_to_3i(&tmp_mat);

    if mat_get_determinant_i3(&tmp_mat_int).abs() != ratio {
        debug::info_print(format_args!(
            "spglib: Determinant of change of basis matrix has to be same as volume ratio.\n"
        ));
        return None;
    }

    // 检查原子数是否能被比率整除
    if cell.size % ratio as usize != 0 {
        debug::info_print(format_args!("spglib: atom number ratio is inconsistent.\n"));
        return None;
    }

    let trimmed_size = cell.size / ratio as usize;
    let mut trimmed_cell = Cell::new(trimmed_size, cell.tensor_rank);

    // 将原子坐标转换到 trimmed_lattice 坐标系下
    // 对应 C: translate_atoms_in_trimmed_lattice
    let position_vec_dbl = translate_atoms_in_trimmed_lattice(cell, &tmp_mat_int)?;

    // 复制晶格和非周期轴信息
    trimmed_cell.lattice = *trimmed_lattice;
    trimmed_cell.aperiodic_axis = cell.aperiodic_axis;

    // 获取重叠表
    // 对应 C: get_overlap_table
    let overlap_table = get_overlap_table(
        &position_vec_dbl,
        cell.size,
        &cell.types,
        &trimmed_cell,
        symprec,
    )?;

    // 构建 mapping_table 并设置 trimmed_cell 的类型
    let mut index_atom = 0;
    for i in 0..cell.size {
        if overlap_table[i] == i {
            mapping_table[i] = index_atom;
            trimmed_cell.types[index_atom] = cell.types[i];
            index_atom += 1;
        } else {
            mapping_table[i] = mapping_table[overlap_table[i]];
        }
    }

    // 计算平均位置和张量
    // 对应 C: set_positions_and_tensors
    set_positions_and_tensors(
        &mut trimmed_cell,
        &position_vec_dbl,
        cell.tensor_rank,
        &cell.tensors,
        mapping_table,
        &overlap_table,
    );

    // Rust 中 VecDBL 会自动 Drop，不需要显式 free
    // mat_free_vecdbl(position_vec_dbl); 

    Some(trimmed_cell)
}

/// 平均重叠原子的位置和张量
/// 对应 C: set_positions_and_tensors
fn set_positions_and_tensors(
    trimmed_cell: &mut Cell,
    position: &VecDBL,
    tensor_rank: TensorRank,
    tensors: &[f64],
    mapping_table: &[usize],
    overlap_table: &[usize],
) {
    // 初始化位置和张量为 0
    for i in 0..trimmed_cell.size {
        trimmed_cell.position[i] = [0.0; 3];
        match tensor_rank {
            TensorRank::Collinear => trimmed_cell.tensors[i] = 0.0,
            TensorRank::NonCollinear => {
                for j in 0..3 {
                    trimmed_cell.tensors[3 * i + j] = 0.0;
                }
            }
            _ => {}
        }
    }

    // 累加位置
    for i in 0..position.size {
        let j = mapping_table[i];
        let k = overlap_table[i]; // k 是该组的代表原子索引
        for l in 0..3 {
            // 边界处理：如果两个重叠原子分别在 0.01 和 0.99，
            // 它们的差值 > 0.5，需要将其中一个移位 1.0 后再求和，以保证平均值正确
            if mat_dabs(position.vec[k][l] - position.vec[i][l]) > 0.5 {
                if position.vec[i][l] < position.vec[k][l] {
                    trimmed_cell.position[j][l] += position.vec[i][l] + 1.0;
                } else {
                    trimmed_cell.position[j][l] += position.vec[i][l] - 1.0;
                }
            } else {
                trimmed_cell.position[j][l] += position.vec[i][l];
            }
        }
    }

    // 累加张量
    for i in 0..position.size {
        let atom_idx = mapping_table[i];
        match tensor_rank {
            TensorRank::Collinear => {
                trimmed_cell.tensors[atom_idx] += tensors[i];
            }
            TensorRank::NonCollinear => {
                for j in 0..3 {
                    trimmed_cell.tensors[3 * atom_idx + j] += tensors[3 * i + j];
                }
            }
            _ => {}
        }
    }

    // 计算平均值
    let multi = (position.size / trimmed_cell.size) as f64;
    for i in 0..trimmed_cell.size {
        for j in 0..3 {
            trimmed_cell.position[i][j] /= multi;
            if trimmed_cell.aperiodic_axis.map_or(true, |ap| j != ap.axis_index()) {
                trimmed_cell.position[i][j] = mat_dmod1(trimmed_cell.position[i][j]);
            }
        }

        match tensor_rank {
            TensorRank::Collinear => {
                trimmed_cell.tensors[i] /= multi;
            }
            TensorRank::NonCollinear => {
                for j in 0..3 {
                    trimmed_cell.tensors[3 * i + j] /= multi;
                }
            }
            _ => {}
        }
    }
}

/// 将原子转换到 trimmed lattice 坐标系
/// 对应 C: translate_atoms_in_trimmed_lattice
fn translate_atoms_in_trimmed_lattice(cell: &Cell, tmat_p_i: &[[i32; 3]; 3]) -> Option<VecDBL> {
    let mut position = VecDBL::new(cell.size);

    for i in 0..cell.size {
        // 假设 mat_multiply_matrix_vector_id3 返回计算后的向量
        position.vec[i] = mat_multiply_matrix_vector_id3(tmat_p_i, &cell.position[i]);
        for j in 0..3 {
            if cell.aperiodic_axis.map_or(true, |ap| j != ap.axis_index()) {
                position.vec[i][j] = mat_dmod1(position.vec[i][j]);
            }
        }
    }
    Some(position)
}

/// 获取重叠表，包含自适应容差调整逻辑
/// 对应 C: get_overlap_table
fn get_overlap_table(
    position: &VecDBL,
    cell_size: usize,
    cell_types: &[i32],
    trimmed_cell: &Cell,
    symprec: f64,
) -> Option<Vec<usize>> {
    let lattice_rank = match trimmed_cell.aperiodic_axis {
        Some(ap) => ap.lattice_rank(), // 2
        None => 3,
    };

    let mut trim_tolerance = symprec;
    let ratio = cell_size / trimmed_cell.size;
    // 预分配内存，避免在循环中重复分配
    let mut overlap_table = vec![0; cell_size];

    // 尝试循环，对应 C 代码中的 goto cont 逻辑
    for attempt in 0..NUM_ATTEMPT {
        // 初始化 overlap_table
        for i in 0..cell_size {
            overlap_table[i] = i;
            for j in 0..cell_size {
                if cell_types[i] == cell_types[j] {
                    let is_overlap = if lattice_rank == 3 {
                        cel_is_overlap(
                            &position.vec[i],
                            &position.vec[j],
                            &trimmed_cell.lattice,
                            trim_tolerance,
                        )
                    } else {
                        cel_layer_is_overlap(
                            &position.vec[i],
                            &position.vec[j],
                            &trimmed_cell.lattice,
                            trimmed_cell.aperiodic_axis.unwrap(),
                            trim_tolerance,
                        )
                    };

                    if is_overlap {
                        if overlap_table[j] == j {
                            overlap_table[i] = j;
                            break;
                        }
                    }
                }
            }
        }

        // 检查重叠数量是否符合预期
        let mut success = true;
        let mut should_increase = false;
        let mut should_reduce = false;

        for i in 0..cell_size {
            if overlap_table[i] != i {
                continue;
            }

            let mut num_overlap = 0;
            for j in 0..cell_size {
                if i == overlap_table[j] {
                    num_overlap += 1;
                }
            }

            if num_overlap == ratio {
                continue;
            }

            // 如果任何一组的重叠数不符合比率，标记失败并决定调整方向
            if num_overlap < ratio {
                should_increase = true;
                success = false;
                break; // 对应 C 的 goto cont
            }
            if num_overlap > ratio {
                should_reduce = true;
                success = false;
                break; // 对应 C 的 goto cont
            }
        }

        if success {
            return Some(overlap_table);
        }

        // 调整容差并重试
        if should_increase {
            trim_tolerance *= INCREASE_RATE;
            debug::debug_print(format_args!(
                "spglib: Increase tolerance to {}\n",
                trim_tolerance
            ));
        } else if should_reduce {
            trim_tolerance *= REDUCE_RATE;
            debug::debug_print(format_args!(
                "spglib: Reduce tolerance to {} ({})\n",
                trim_tolerance, attempt
            ));
        }
    }

    debug::info_print(format_args!("spglib: Could not trim cell well\n"));
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cell_creation() {
        let cell = Cell::new(2, TensorRank::NoSpin);
        assert_eq!(cell.size, 2);
        assert_eq!(cell.tensors.len(), 0);

        let cell_col = Cell::new(2, TensorRank::Collinear);
        assert_eq!(cell_col.tensors.len(), 2);
    }

    #[test]
    fn test_overlap() {
        let lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let a = [0.1, 0.1, 0.1];
        let b = [0.100000001, 0.1, 0.1];
        let c = [0.5, 0.5, 0.5];

        assert!(cel_is_overlap(&a, &b, &lattice, 1e-5));
        assert!(!cel_is_overlap(&a, &c, &lattice, 1e-5));

        // Periodic boundary
        let d = [1.1, 0.1, 0.1];
        assert!(cel_is_overlap(&a, &d, &lattice, 1e-5));
    }

    #[test]
    fn test_trim_cell_simple() {
        // 构造一个 2x1x1 的超胞
        // Primitive: a=1, b=1, c=1. Atom at (0,0,0)
        // Supercell: a'=2, b'=1, c'=1. Atoms at (0,0,0) and (0.5, 0, 0)

        let trimmed_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let super_lattice = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        let mut cell = Cell::new(2, TensorRank::NoSpin);
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let types = [1, 1];
        cell.set_cell(&super_lattice, &positions, &types);

        let mut mapping_table = vec![0; 2];
        let trimmed = cel_trim_cell(&mut mapping_table, &trimmed_lattice, &cell, 1e-5);

        assert!(trimmed.is_some());
        let t = trimmed.unwrap();
        assert_eq!(t.size, 1);
        assert!((t.position[0][0] - 0.0).abs() < 1e-5);

        // Check mapping: both atoms map to 0
        assert_eq!(mapping_table[0], 0);
        assert_eq!(mapping_table[1], 0);
    }
}
