//! Niggli 晶格约化。
//!
//! 将三维晶格约化为 Niggli 标准形式。用于三斜和单斜晶系（对更高对称性的晶系，
//! 使用 Delaunay 约化）。Niggli 约化后的晶格满足一组严格的不等式条件，
//! 确保同一晶格的所有等价表示约化到相同的标准形式。

use crate::SymError;
use crate::cell::AperiodicAxis;
use crate::debug;
use crate::mathfunc::{Mat3, mat_get_metric, mat_multiply_matrix_d3};
use std::env;

// 版本常量
pub const NIGGLI_MAJOR_VERSION: i32 = 2;
pub const NIGGLI_MINOR_VERSION: i32 = 0;
pub const NIGGLI_MICRO_VERSION: i32 = 0;

/// Niggli 约化参数结构体
struct NiggliParams {
    a: f64,
    b: f64,
    c: f64,
    eta: f64,
    xi: f64,
    zeta: f64,
    eps: f64,
    l: i32,
    m: i32,
    n: i32,
    tmat: Mat3,    // 变换矩阵
    lattice: Mat3, // 当前晶格
}

impl NiggliParams {
    fn new(lattice: &Mat3, eps: f64) -> Self {
        NiggliParams {
            a: 0.0,
            b: 0.0,
            c: 0.0,
            eta: 0.0,
            xi: 0.0,
            zeta: 0.0,
            eps,
            l: 0,
            m: 0,
            n: 0,
            tmat: [[0.0; 3]; 3],
            lattice: *lattice,
        }
    }
}

/// 获取最大尝试次数
fn get_num_attempts() -> i32 {
    if let Ok(val_str) = env::var("SPGLIB_NUM_ATTEMPTS") {
        if let Ok(val) = val_str.parse::<i32>()
            && val > 0 {
                return val;
            }
        debug::warning_print(format_args!(
            "spglib: Could not parse SPGLIB_NUM_ATTEMPTS={}\n",
            val_str
        ));
    }
    1000
}

pub fn get_major_version() -> i32 {
    NIGGLI_MAJOR_VERSION
}
pub fn get_minor_version() -> i32 {
    NIGGLI_MINOR_VERSION
}
pub fn get_micro_version() -> i32 {
    NIGGLI_MICRO_VERSION
}

/// Niggli 约化主函数
/// lattice: 输入/输出晶格 (列向量矩阵)
/// eps: 容差
/// aperiodic_axis: None=体材料, Some(X/Y/Z)=非周期轴
/// 成功返回 `Ok(())`，失败返回 `Err(SymError::NiggliFailed)`。
pub fn reduce(lattice: &mut Mat3, eps: f64, aperiodic_axis: Option<AperiodicAxis>) -> Result<(), SymError> {
    let mut p = NiggliParams::new(lattice, eps);
    let mut succeeded = false;

    // Step 0: Move aperiodic axis to c for layer
    let configured = match aperiodic_axis {
        Some(axis @ (AperiodicAxis::X | AperiodicAxis::Y)) => layer_swap_axis(&mut p, axis),
        None | Some(AperiodicAxis::Z) => set_parameters(&mut p),
    };
    if !configured {
        return Err(SymError::NiggliFailed);
    }

    let max_attempts = get_num_attempts();
    for _ in 0..max_attempts {
        // C 代码逻辑：
        // 依次执行 Step 1 到 8。
        // 如果某个 Step 执行了变换 (返回 true)：
        //   1. 打印 debug 信息
        //   2. 调用 reset 更新参数
        //   3. 检查是否需要 break (重启外层循环)。
        //      - Step 1, 3, 4: 不需要 break，继续执行后续 Step。
        //      - Step 2, 5, 6, 7, 8: 需要 break，重新从 Step 1 开始。
        // 如果所有 Step 都没有触发 break (即顺利走到了 Step 8 之后)，则认为收敛成功。

        // Step 1
        if step1(&mut p) {
            debug_show(1, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            // No break for Step 1
        }

        // Step 2
        let do_step2 = if aperiodic_axis.is_some() {
            step2_for_layer(&mut p)
        } else {
            step2(&mut p)
        };
        if do_step2 {
            debug_show(2, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            continue; // Break inner loop, restart outer loop
        }

        // Step 3
        if step3(&mut p) {
            debug_show(3, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            // No break for Step 3
        }

        // Step 4
        if step4(&mut p) {
            debug_show(4, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            // No break for Step 4
        }

        // Step 5
        if step5(&mut p) {
            debug_show(5, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            continue; // Break inner loop, restart outer loop
        }

        // Step 6
        if step6(&mut p) {
            debug_show(6, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            continue; // Break inner loop, restart outer loop
        }

        // Step 7
        if step7(&mut p) {
            debug_show(7, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            continue; // Break inner loop, restart outer loop
        }

        // Step 8
        if step8(&mut p) {
            debug_show(8, &p);
            if !reset(&mut p) {
                return Err(SymError::NiggliFailed);
            }
            continue; // Break inner loop, restart outer loop
        }

        // 如果能运行到这里，说明没有触发任何 continue，即没有触发 Step 2, 5, 6, 7, 8 的变换。
        // 即使 Step 1, 3, 4 发生了变换，只要不触发重启条件的步骤，也视为一轮完成。
        // (C 代码中 j 循环结束意味着 j==8)
        succeeded = true;
        break;
    }

    debug_show(-1, &p);

    // Finalize: copy back to lattice
    *lattice = p.lattice;

    if succeeded { Ok(()) } else { Err(SymError::NiggliFailed) }
}

fn layer_swap_axis(p: &mut NiggliParams, aperiodic_axis: AperiodicAxis) -> bool {
    if aperiodic_axis == AperiodicAxis::X {
        p.tmat = [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]];
    } else { // Y
        p.tmat = [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]];
    }
    reset(p)
}

fn reset(p: &mut NiggliParams) -> bool {
    // lattice = lattice * tmat
    p.lattice = mat_multiply_matrix_d3(&p.lattice, &p.tmat);
    set_parameters(p)
}

fn set_parameters(p: &mut NiggliParams) -> bool {
    let g = mat_get_metric(&p.lattice);
    // G = [g00, g01, g02]
    //     [g10, g11, g12]
    //     [g20, g21, g22]
    // A=g00, B=g11, C=g22
    // xi=2*g12, eta=2*g02, zeta=2*g01

    p.a = g[0][0];
    p.b = g[1][1];
    p.c = g[2][2];
    p.xi = g[1][2] * 2.0;
    p.eta = g[0][2] * 2.0;
    p.zeta = g[0][1] * 2.0;

    set_angle_types(p);
    true
}

fn set_angle_types(p: &mut NiggliParams) {
    p.l = 0;
    p.m = 0;
    p.n = 0;
    if p.xi < -p.eps {
        p.l = -1;
    }
    if p.xi > p.eps {
        p.l = 1;
    }
    if p.eta < -p.eps {
        p.m = -1;
    }
    if p.eta > p.eps {
        p.m = 1;
    }
    if p.zeta < -p.eps {
        p.n = -1;
    }
    if p.zeta > p.eps {
        p.n = 1;
    }
}

fn debug_show(j: i32, p: &NiggliParams) {
    if j < 0 {
        debug::debug_print(format_args!("Finish: "));
    } else {
        debug::debug_print(format_args!("Step {}: ", j));
    }
    debug::debug_print(format_args!(
        "{:.6} {:.6} {:.6} {:.6} {:.6} {:.6}\n",
        p.a, p.b, p.c, p.xi, p.eta, p.zeta
    ));
}

// --- Steps ---

fn step1(p: &mut NiggliParams) -> bool {
    if p.a > p.b + p.eps || ((p.a - p.b).abs() <= p.eps && p.xi.abs() > p.eta.abs() + p.eps) {
        p.tmat = [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]];
        true
    } else {
        false
    }
}

fn step2(p: &mut NiggliParams) -> bool {
    if p.b > p.c + p.eps || ((p.b - p.c).abs() <= p.eps && p.eta.abs() > p.zeta.abs() + p.eps) {
        p.tmat = [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]];
        true
    } else {
        false
    }
}

fn step2_for_layer(p: &mut NiggliParams) -> bool {
    if p.b > p.c + p.eps || ((p.b - p.c).abs() <= p.eps && p.eta.abs() > p.zeta.abs() + p.eps) {
        debug::debug_print(format_args!(
            "niggli: B > C or B = C and |eta| > |zeta|. Please elongate the aperiodic axis.\n"
        ));
    }
    false
}

fn step3(p: &mut NiggliParams) -> bool {
    if p.l * p.m * p.n == 1 {
        let i = if p.l == -1 { -1.0 } else { 1.0 };
        let j = if p.m == -1 { -1.0 } else { 1.0 };
        let k = if p.n == -1 { -1.0 } else { 1.0 };
        p.tmat = [[i, 0.0, 0.0], [0.0, j, 0.0], [0.0, 0.0, k]];
        true
    } else {
        false
    }
}

fn step4(p: &mut NiggliParams) -> bool {
    if p.l == -1 && p.m == -1 && p.n == -1 {
        return false;
    }

    if p.l * p.m * p.n == 0 || p.l * p.m * p.n == -1 {
        let mut i = 1.0;
        let mut j = 1.0;
        let mut k = 1.0;
        let mut r = -1; // 0: i, 1: j, 2: k

        if p.l == 1 {
            i = -1.0;
        }
        if p.l == 0 {
            r = 0;
        }
        if p.m == 1 {
            j = -1.0;
        }
        if p.m == 0 {
            r = 1;
        }
        if p.n == 1 {
            k = -1.0;
        }
        if p.n == 0 {
            r = 2;
        }

        if i * j * k == -1.0 {
            if r == 0 {
                i = -1.0;
            }
            if r == 1 {
                j = -1.0;
            }
            if r == 2 {
                k = -1.0;
            }
        }

        p.tmat = [[i, 0.0, 0.0], [0.0, j, 0.0], [0.0, 0.0, k]];
        true
    } else {
        false
    }
}

fn step5(p: &mut NiggliParams) -> bool {
    if p.xi.abs() > p.b + p.eps
        || ((p.b - p.xi).abs() <= p.eps && 2.0 * p.eta < p.zeta - p.eps)
        || ((p.b + p.xi).abs() <= p.eps && p.zeta < -p.eps)
    {
        p.tmat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        if p.xi > 0.0 {
            p.tmat[1][2] = -1.0; // tmat[5] in C (row 1, col 2)
        }
        if p.xi < 0.0 {
            p.tmat[1][2] = 1.0;
        }
        true
    } else {
        false
    }
}

fn step6(p: &mut NiggliParams) -> bool {
    if p.eta.abs() > p.a + p.eps
        || ((p.a - p.eta).abs() <= p.eps && 2.0 * p.xi < p.zeta - p.eps)
        || ((p.a + p.eta).abs() <= p.eps && p.zeta < -p.eps)
    {
        p.tmat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        if p.eta > 0.0 {
            p.tmat[0][2] = -1.0; // tmat[2] in C (row 0, col 2)
        }
        if p.eta < 0.0 {
            p.tmat[0][2] = 1.0;
        }
        true
    } else {
        false
    }
}

fn step7(p: &mut NiggliParams) -> bool {
    if p.zeta.abs() > p.a + p.eps
        || ((p.a - p.zeta).abs() <= p.eps && 2.0 * p.xi < p.eta - p.eps)
        || ((p.a + p.zeta).abs() <= p.eps && p.eta < -p.eps)
    {
        p.tmat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        if p.zeta > 0.0 {
            p.tmat[0][1] = -1.0; // tmat[1] in C (row 0, col 1)
        }
        if p.zeta < 0.0 {
            p.tmat[0][1] = 1.0;
        }
        true
    } else {
        false
    }
}

fn step8(p: &mut NiggliParams) -> bool {
    if p.xi + p.eta + p.zeta + p.a + p.b < -p.eps
        || ((p.xi + p.eta + p.zeta + p.a + p.b).abs() <= p.eps
            && 2.0 * (p.a + p.eta) + p.zeta > p.eps)
    {
        p.tmat = [[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]];
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_identity() {
        let mut lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        reduce(&mut lattice, 1e-5, None).unwrap();
        // Identity should remain identity
        assert!((lattice[0][0] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_reduce_swap() {
        // a=2, b=1, c=1 -> should swap a and b (or b and c depending on steps)
        // Input: a=(2,0,0), b=(0,1,0), c=(0,0,1)
        // A=4, B=1, C=1.
        // Step 1: A > B? 4 > 1. Yes. Swap a, b.
        // New: a=(0,1,0), b=(2,0,0), c=(0,0,1). A=1, B=4, C=1.
        // Step 2: B > C? 4 > 1. Yes. Swap b, c.
        // New: a=(0,1,0), b=(0,0,1), c=(2,0,0). A=1, B=1, C=4.
        // Sorted.
        let mut lattice = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        reduce(&mut lattice, 1e-5, None).unwrap();

        let g = mat_get_metric(&lattice);
        assert!((g[0][0] - 1.0).abs() < 1e-5); // A=1
        assert!((g[1][1] - 1.0).abs() < 1e-5); // B=1
        assert!((g[2][2] - 4.0).abs() < 1e-5); // C=4
    }
}
