//! k 点网格生成。
//!
//! 生成 Monkhorst-Pack 风格的 k 点网格地址，用于能带计算和布里渊区采样。

use crate::SymError;

/// Upper bound used by allocating high-level grid APIs (20 MiB for one
/// address-plus-mapping pair, and at most 160 MiB for an eightfold BZ buffer).
pub(crate) const MAX_GRID_POINTS: usize = 1_000_000;

pub(crate) fn mesh_size(mesh: &[i32; 3]) -> Result<usize, SymError> {
    if mesh.iter().any(|&component| component <= 0) {
        return Err(SymError::InvalidInput);
    }

    mesh.iter()
        .try_fold(1usize, |product, &component| {
            product.checked_mul(component as usize)
        })
        .ok_or(SymError::ArraySizeShortage)
}

pub(crate) fn validate_mesh(mesh: &[i32; 3]) -> Result<usize, SymError> {
    let total = mesh_size(mesh)?;
    if total <= MAX_GRID_POINTS {
        Ok(total)
    } else {
        Err(SymError::ArraySizeShortage)
    }
}

pub(crate) fn validate_shift(is_shift: &[i32; 3]) -> Result<(), SymError> {
    if is_shift.iter().all(|&shift| shift == 0 || shift == 1) {
        Ok(())
    } else {
        Err(SymError::InvalidInput)
    }
}

/// 获取所有网格点的地址 (Monkhorst-Pack Grid Generation)
///
/// 该函数遍历由 `mesh` 定义的倒易空间网格，计算每个点的线性索引，
/// 并将对应的三维网格坐标（已约化到第一布里渊区附近）存储在 `grid_address` 中。
///
/// # Arguments
/// * `grid_address` - 输出参数，用于存储网格点坐标。大小必须 >= `mesh[0]*mesh[1]*mesh[2]`。
/// * `mesh` - 定义 k 点网格的尺寸 [Nx, Ny, Nz]。
pub(crate) fn get_all_grid_addresses(
    grid_address: &mut [[i32; 3]],
    mesh: &[i32; 3],
) -> Result<(), SymError> {
    let total_size = validate_mesh(mesh)?;
    if grid_address.len() < total_size {
        return Err(SymError::ArraySizeShortage);
    }

    let mut address = [0; 3];

    // 对应 C 代码中的三重循环
    // 遍历顺序：x -> y -> z (最内层)
    // 注意：get_grid_point_single_mesh 的默认实现是 z 变化最快 (Row-major like if mapped to 3D array)
    for i in 0..mesh[0] {
        address[0] = i;
        for j in 0..mesh[1] {
            address[1] = j;
            for k in 0..mesh[2] {
                address[2] = k;

                // 获取线性索引
                let grid_point = get_grid_point_single_mesh(&address, mesh);

                // 直接写入，Rust 会自动进行边界检查（由于前面已 assert，这里是安全的）
                // 复制当前坐标
                grid_address[grid_point] = address;

                // 执行边界约化，将坐标映射到以 Gamma 点为中心的区域
                // 例如 mesh=4: [0, 1, 2, 3] -> [0, 1, 2, -1]
                reduce_grid_address(&mut grid_address[grid_point], mesh);
            }
        }
    }
    Ok(())
}

/// 获取双倍网格下的线性索引
///
/// 对应 C: kgd_get_grid_point_double_mesh
pub(crate) fn get_grid_point_double_mesh(
    address_double: &[i32; 3],
    mesh: &[i32; 3],
) -> Result<usize, SymError> {
    mesh_size(mesh)?;
    Ok(get_grid_point_double_mesh_core(address_double, mesh))
}

/// 获取双倍网格下的线性索引 (Dense 版本)
///
/// 对应 C: kgd_get_dense_grid_point_double_mesh
/// 在当前实现中，逻辑与 get_grid_point_double_mesh 相同
pub(crate) fn get_dense_grid_point_double_mesh(
    address_double: &[i32; 3],
    mesh: &[i32; 3],
) -> Result<usize, SymError> {
    get_grid_point_double_mesh(address_double, mesh)
}

/// 计算双倍网格地址
///
/// 根据原始网格地址和位移（shift），计算在双倍密度网格上的坐标。
/// 公式: address_double = address * 2 + shift
///
/// # Arguments
/// * `address_double` - 输出：双倍网格坐标
/// * `address` - 输入：原始网格坐标
/// * `mesh` - 原始网格尺寸
/// * `is_shift` - 是否有半网格位移 (0 或 1)
pub(crate) fn get_grid_address_double_mesh(
    address_double: &mut [i32; 3],
    address: &[i32; 3],
    mesh: &[i32; 3],
    is_shift: &[i32; 3],
) -> Result<(), SymError> {
    mesh_size(mesh)?;
    validate_shift(is_shift)?;
    for i in 0..3 {
        let period = i64::from(mesh[i]) * 2;
        let mut reduced = (i64::from(address[i]) * 2 + i64::from(is_shift[i])).rem_euclid(period);
        if reduced > i64::from(mesh[i]) {
            reduced -= period;
        }
        address_double[i] = reduced as i32;
    }
    Ok(())
}

// --- Internal Helper Functions ---

/// 内部函数：从双倍网格坐标计算单倍网格的线性索引
#[inline]
fn get_grid_point_double_mesh_core(address_double: &[i32; 3], mesh: &[i32; 3]) -> usize {
    let mut address = [0; 3];
    for i in 0..3 {
        // 模拟整数除法的向下取整逻辑，还原原始网格坐标
        // 如果是偶数：val / 2
        // 如果是奇数：(val - 1) / 2
        address[i] = i64::from(address_double[i])
            .div_euclid(2)
            .rem_euclid(i64::from(mesh[i])) as i32;
    }
    get_grid_point_single_mesh(&address, mesh)
}

/// 内部函数：计算单倍网格的线性索引
///
/// 对应 C 中的 #ifndef GRID_ORDER_XYZ 分支
/// Index = z * (Nx * Ny) + y * (Nx) + x
#[inline]
fn get_grid_point_single_mesh(address: &[i32; 3], mesh: &[i32; 3]) -> usize {
    (address[2] as usize * mesh[0] as usize * mesh[1] as usize)
        + (address[1] as usize * mesh[0] as usize)
        + (address[0] as usize)
}

#[cfg(test)]
fn modulo_i3(values: &mut [i32; 3], modulus: &[i32; 3]) {
    for (value, &base) in values.iter_mut().zip(modulus) {
        *value %= base;
        if *value < 0 {
            *value += base;
        }
    }
}

/// 内部函数：单倍网格坐标约化
///
/// 将坐标映射到 (-mesh/2, mesh/2] 区间
/// 对应 C 中的 #ifndef GRID_BOUNDARY_AS_NEGATIVE 分支
#[inline]
fn reduce_grid_address(address: &mut [i32; 3], mesh: &[i32; 3]) {
    for i in 0..3 {
        // C: address[i] -= mesh[i] * (address[i] > mesh[i] / 2);
        if address[i] > mesh[i] / 2 {
            address[i] -= mesh[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_all_grid_addresses() {
        let mesh = [2, 2, 1];
        // Total size = 4
        let mut addresses = vec![[0; 3]; 4];
        get_all_grid_addresses(&mut addresses, &mesh).unwrap();

        // Expected:
        // (0,0,0) -> idx 0
        // (0,1,0) -> idx 2 (because idx = z*Nx*Ny + y*Nx + x -> 0 + 1*2 + 0 = 2)
        // (1,0,0) -> idx 1
        // (1,1,0) -> idx 3
        // Reduced:
        // 1 > 2/2 (false) -> remains 1.

        // Let's trace get_grid_point_single_mesh for [0,1,0], mesh=[2,2,1]
        // idx = 0*2*2 + 1*2 + 0 = 2. Correct.

        assert_eq!(addresses[0], [0, 0, 0]);
        assert_eq!(addresses[1], [1, 0, 0]);
        assert_eq!(addresses[2], [0, 1, 0]); // y=1, reduced? 1 > 2/2 (1>1 false) -> 1
        assert_eq!(addresses[3], [1, 1, 0]);
    }

    #[test]
    fn test_reduce_grid_address() {
        let mesh = [4, 4, 4];
        let mut addr = [0, 2, 3];
        reduce_grid_address(&mut addr, &mesh);
        // 0 -> 0
        // 2 > 4/2 (2>2 false) -> 2
        // 3 > 4/2 (3>2 true) -> 3 - 4 = -1
        assert_eq!(addr, [0, 2, -1]);
    }

    #[test]
    fn test_grid_point_single() {
        let mesh = [4, 4, 4];
        let addr = [1, 0, 0];
        assert_eq!(get_grid_point_single_mesh(&addr, &mesh), 1);
        let addr2 = [0, 1, 0];
        assert_eq!(get_grid_point_single_mesh(&addr2, &mesh), 4);
    }

    #[test]
    fn test_modulo() {
        let mesh = [4, 4, 4];
        let mut v = [-1, 5, 0];
        modulo_i3(&mut v, &mesh);
        assert_eq!(v, [3, 1, 0]);
    }

    #[test]
    fn test_double_mesh_logic() {
        let mesh = [4, 4, 4];
        // Case 1: Even address_double
        let addr_d = [0, 4, 8];
        // 0/2=0, 4/2=2, 8/2=4 -> mod 4 -> 0
        assert_eq!(
            get_grid_point_double_mesh_core(&addr_d, &mesh),
            get_grid_point_single_mesh(&[0, 2, 0], &mesh)
        );

        // Case 2: Odd address_double
        let addr_d_odd = [1, 5, 9];
        // (1-1)/2=0, (5-1)/2=2, (9-1)/2=4 -> mod 4 -> 0
        assert_eq!(
            get_grid_point_double_mesh_core(&addr_d_odd, &mesh),
            get_grid_point_single_mesh(&[0, 2, 0], &mesh)
        );
    }
}
