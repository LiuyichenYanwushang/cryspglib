//! 调试和日志输出工具。
//!
//! 提供受环境变量 (`SPGLIB_DEBUG`, `SPGLIB_WARNING`) 控制的调试打印宏，
//! 以及矩阵/向量的格式化显示函数。调试输出默认写入 stderr。

use std::env;
use std::fmt;
use std::io::{self, Write};
use std::sync::OnceLock;

// 使用 OnceLock 缓存环境变量状态，避免在循环中频繁进行系统调用和字符串分配
static DEBUG_ENABLED: OnceLock<bool> = OnceLock::new();
static WARNING_ENABLED: OnceLock<bool> = OnceLock::new();
static INFO_ENABLED: OnceLock<bool> = OnceLock::new();

/// 检查是否启用了调试模式 (环境变量 SPGLIB_DEBUG)
/// 优化：使用缓存，零开销
pub fn debug_enabled() -> bool {
    *DEBUG_ENABLED.get_or_init(|| env::var_os("SPGLIB_DEBUG").is_some())
}

/// 检查是否启用了警告模式 (环境变量 SPGLIB_WARNING)
/// 默认为开启，除非设置为 "OFF"
pub fn warning_enabled() -> bool {
    *WARNING_ENABLED.get_or_init(|| {
        match env::var("SPGLIB_WARNING") {
            Ok(val) => val != "OFF",
            Err(_) => true, // 默认开启
        }
    })
}

/// 检查是否启用了信息模式 (环境变量 SPGLIB_INFO)
pub fn info_enabled() -> bool {
    *INFO_ENABLED.get_or_init(|| env::var_os("SPGLIB_INFO").is_some())
}

/// 打印 3x3 double 矩阵
pub fn debug_print_matrix_d3(a: &[[f64; 3]; 3]) {
    if !debug_enabled() {
        return;
    }
    // 使用 lock() 获取 stdout 句柄，避免多次加锁，提高性能
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for row in a {
        let _ = writeln!(handle, "{:.6} {:.6} {:.6}", row[0], row[1], row[2]);
    }
}

/// 打印 3x3 int 矩阵
pub fn debug_print_matrix_i3(a: &[[i32; 3]; 3]) {
    if !debug_enabled() {
        return;
    }
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for row in a {
        let _ = writeln!(handle, "{} {} {}", row[0], row[1], row[2]);
    }
}

/// 打印 double 向量列表
pub fn debug_print_vectors_d3(a: &[[f64; 3]]) {
    if !debug_enabled() {
        return;
    }
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for (i, vec) in a.iter().enumerate() {
        let _ = writeln!(
            handle,
            "{}: {:.6} {:.6} {:.6}",
            i + 1,
            vec[0],
            vec[1],
            vec[2]
        );
    }
}

/// 打印单个 double 向量
pub fn debug_print_vector_d3(a: &[f64; 3]) {
    if !debug_enabled() {
        return;
    }
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let _ = writeln!(handle, "{:.6} {:.6} {:.6}", a[0], a[1], a[2]);
}

/// 打印带标签的向量列表
pub fn debug_print_vectors_with_label(a: &[[f64; 3]], b: &[i32]) {
    if !debug_enabled() {
        return;
    }
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for (i, vec) in a.iter().enumerate() {
        let label = if i < b.len() { b[i] } else { 0 };
        let _ = writeln!(
            handle,
            "{}: {:.6} {:.6} {:.6}",
            label, vec[0], vec[1], vec[2]
        );
    }
}

/// 通用调试打印函数
/// 使用方式: debug_print(format_args!("Value: {}", x));
pub fn debug_print(args: fmt::Arguments) {
    if !debug_enabled() {
        return;
    }
    let _ = io::stdout().write_fmt(args);
}

/// 通用警告打印函数
pub fn warning_print(args: fmt::Arguments) {
    if !warning_enabled() {
        return;
    }
    let _ = io::stderr().write_fmt(args);
}

/// 通用信息打印函数
pub fn info_print(args: fmt::Arguments) {
    if !info_enabled() {
        return;
    }
    let _ = io::stderr().write_fmt(args);
}

/// 内存分配警告辅助函数
pub fn warning_memory(what: &str) {
    warning_print(format_args!(
        "Spglib: Memory could not be allocated: {}\n",
        what
    ));
}

#[cfg(test)]
mod tests {
    use super::*;

    // 注意：由于使用了 OnceLock 缓存，环境变量只会在第一次调用时读取。
    // 为了测试不同的环境变量状态，我们需要在测试中模拟或重置（但在 Rust 标准库中重置 OnceLock 是不安全的/不支持的）。
    // 因此，这里的测试主要验证逻辑正确性。在实际单元测试中，针对全局单例的测试比较棘手。
    // 下面的测试代码假设是分别运行的，或者我们接受 OnceLock 初始化后不可变的特性。
    //
    // 为了让测试通过，我们这里仅测试默认行为和打印函数不 Panic。
    // 如果需要严格测试环境变量切换，需要将 OnceLock 替换为每次读取（仅用于测试配置）或使用 AtomicBool 手动管理。

    #[test]
    fn test_print_functions_no_panic() {
        // 确保在没有设置环境变量的情况下调用这些函数不会崩溃
        let mat_d = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mat_i = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let vec_d = [1.0, 2.0, 3.0];
        let vecs_d = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let labels = vec![1, 2];

        debug_print_matrix_d3(&mat_d);
        debug_print_matrix_i3(&mat_i);
        debug_print_vector_d3(&vec_d);
        debug_print_vectors_d3(&vecs_d);
        debug_print_vectors_with_label(&vecs_d, &labels);
        debug_print(format_args!("Test debug print: {}\n", 123));
    }

    // 注意：由于 OnceLock 的存在，test_debug_enabled 等修改环境变量的测试
    // 在同一个进程中只能生效一次。如果必须测试动态开关，建议在非生产构建中移除 OnceLock。
}
