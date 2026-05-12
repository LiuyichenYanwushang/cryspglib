//! # irrep 模块使用示例
//!
//! 演示三种查询不可约表示特征标的方式:
//! 1. 通过 Crystal → SpaceGroup 桥接 API
//! 2. 直接按 SG 编号查询
//! 3. 格式化特征标表

use cryspglib::irrep::query::*;
use cryspglib::SymmetryOps;
use cryspglib::Crystal;

fn main() {
    // ━━━ 方式一: 通过 Crystal 桥接 API ━━━
    // 适合: 已有晶体结构, 先识别空间群再查 irrep

    let al = Crystal::new(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
        vec![13, 13, 13, 13],
    );
    let ds = al.analyze().symprec(1e-5).dataset().unwrap();

    // 桥接 API: 从 SpaceGroup 直接查 irreps
    println!("SG {} ({}):", ds.spacegroup_number, ds.international_symbol);
    let gm_irreps = ds.irreps_at_k("GM");
    println!("  Gamma point: {} irreps", gm_irreps.len());

    let gm4m = gm_irreps.iter().find(|r| r.ml == "GM4-").unwrap();
    println!("  {}: dim={}, image={}", gm4m.ml, gm4m.dim, gm4m.image);

    // 对称操作 — spglib Hall 顺序
    let ops = ds.symmetry_ops();
    println!("  {} ops in spglib order", ops.len());

    // 特征标 — 直接在 spglib Hall 顺序中, chars[i] ↔ ops[i]
    let chars = gm4m.characters();
    println!();
    println!("  idx  χ(GM4-)");
    println!("  ---  -------");
    for i in 0..ops.len().min(6) {
        println!("  {:3}  {:6.1}", i, chars.get(i).copied().unwrap_or(f64::NAN));
    }

    // ━━━ 方式二: 直接按 SG 编号查询 (不需要 Crystal) ━━━
    // 适合: 已知空间群编号, 直接查 irrep 数据

    let irreps = irreps_of(221);
    let gm4m_direct = irreps.iter()
        .filter(|r| r.k_label() == "GM")
        .find(|r| r.ml == "GM4-").unwrap();

    // characters() 返回 spglib Hall 顺序的特征标 (100% 覆盖)
    println!("\n  Direct lookup characters: {:?}",
        &gm4m_direct.characters()[..gm4m_direct.characters().len().min(6)]);

    // ━━━ 方式三: 格式化特征标表 (人类可读) ━━━
    println!("\n{}", format_character_table(221, 0, 0, 0, 1));

    // ━━━ 查询 isotropy subgroup ━━━
    println!("Isotropy subgroups of GM4-:");
    for sub in gm4m_direct.subgroups() {
        println!("  SG #{} {}  direction={} domains={}",
            sub.sg, sub.symbol, sub.direction, sub.domains);
    }

    // ━━━ 磁共表示 (corepresentation) ━━━
    let corep = gm4m_direct.corepresentation(1599);
    if let Some(c) = corep {
        println!("\nMagnetic corep for 221.97 (UNI 1599) at GM:");
        println!("  type={:?}, dim={}, χ̃(id)={:.2}", c.corep_type, c.dim, c.characters[0]);
    }

    // ━━━ 遍历所有 k 点 ━━━
    println!("\nAll k-points in SG 221:");
    for kp in kpoints_of(221) {
        println!("  {} ({} irreps)", kp.label, kp.irreps.len());
    }

    // ━━━ 获取对称操作 (不通过 Crystal) ━━━
    let ops221 = SymmetryOps::from_sg(221).unwrap();
    println!("\nSG 221: {} symmetry operations", ops221.len());
    assert_eq!(ops221.len(), 48); // Pm-3m 有 48 个操作
}
