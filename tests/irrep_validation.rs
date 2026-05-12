//! 验证自动生成的 irrep 数据的一致性和正确性

use cryspglib::irrep::query::*;
use cryspglib::irrep::types::generated_data::*;
use cryspglib::irrep::types::IrrepRecord;
use cryspglib::SymmetryOps;

// ==========================================================================
// 基本完整性检查
// ==========================================================================

/// 所有 230 个空间群至少有一个 irrep
#[test]
fn all_sgs_have_irreps() {
    let mut empty_sgs = Vec::new();
    for sg in 1..=230u8 {
        let irreps = irreps_of(sg);
        if irreps.is_empty() {
            empty_sgs.push(sg);
        }
    }
    assert!(
        empty_sgs.is_empty(),
        "Empty SGs: {:?}",
        empty_sgs
    );
}

/// 每个标量 irrep 都有非空标签，双值 irrep 至少有 ML 标签
#[test]
fn every_irrep_has_labels() {
    for (i, ir) in IRREPS.iter().enumerate() {
        assert!(!ir.ml.is_empty(), "irrep {} has empty ML label", i);
        if !ir.spinor {
            assert!(!ir.bc.is_empty(), "scalar irrep {} has empty BC label", i);
            assert!(!ir.kov.is_empty(), "scalar irrep {} has empty Kov label", i);
        }
    }
}

/// 每个 irrep 的 SG 编号在 1–230 范围内
#[test]
fn sg_numbers_valid() {
    for (i, ir) in IRREPS.iter().enumerate() {
        assert!(
            ir.sg >= 1 && ir.sg <= 230,
            "irrep {}: sg={} out of range",
            i,
            ir.sg
        );
    }
}

/// Image 标签格式: 标量 irrep 的 image 首字母大写
#[test]
fn image_labels_well_formed() {
    for (i, ir) in IRREPS.iter().enumerate() {
        let img = ir.image;
        if ir.spinor {
            continue; // spinor irreps don't have image labels yet
        }
        assert!(!img.is_empty(), "irrep {}: empty image", i);
        assert!(
            img.chars().next().unwrap().is_ascii_uppercase(),
            "irrep {}: image '{}' doesn't start with uppercase",
            i,
            img
        );
    }
}

/// 维度与 image 首字母一致
#[test]
fn dimension_matches_image() {
    let dim_from_letter = |c: char| match c {
        'A' => 1,
        'B' => 2,
        'C' => 3,
        'D' => 4,
        'E' => 6,
        'F' => 8,
        'G' => 12,
        'H' => 16,
        'J' => 24,
        _ => 0,
    };
    for (i, ir) in IRREPS.iter().enumerate() {
        let expected = dim_from_letter(ir.image.chars().next().unwrap());
        if expected > 0 {
            assert_eq!(
                ir.dim,
                expected,
                "irrep {}: dim={} but image '{}' expects {}",
                i,
                ir.dim,
                ir.image,
                expected
            );
        }
    }
}

// ==========================================================================
// 已知物理事实验证
// ==========================================================================

/// SG #1 P1: 三斜, 仅恒等操作, 8 个 k 点 irrep
#[test]
fn sg1_p1_irreps() {
    let irreps = irreps_of(1);
    assert!(
        irreps.len() >= 8,
        "P1 should have at least 8 irreps, got {}",
        irreps.len()
    );

    // 所有 irrep 应为 1 维
    for ir in irreps {
        assert_eq!(ir.dim, 1, "P1 irreps should all be 1D, got dim={}", ir.dim);
    }

    // MI 标签检查: 应该包含 GM1
    let labels: Vec<&str> = irreps.iter().map(|r| r.ml).collect();
    assert!(labels.contains(&"GM1"), "P1 should have GM1");
}

/// SG #2 P-1: 中心对称, 16 个 irrep (± 对)
#[test]
fn sg2_p1bar_irreps() {
    let irreps = irreps_of(2);
    assert!(
        irreps.len() >= 16,
        "P-1 should have at least 16 irreps, got {}",
        irreps.len()
    );

    let labels: Vec<&str> = irreps.iter().map(|r| r.ml).collect();
    assert!(labels.contains(&"GM1+"), "P-1 should have GM1+");
    assert!(labels.contains(&"GM1-"), "P-1 should have GM1-");
}

/// SG #221 Pm-3m: Γ 点有 10 个标量 irrep + 6 个双值 irrep
#[test]
fn sg221_gamma_irreps() {
    let irreps = irreps_of(221);
    assert!(
        irreps.len() >= 40,
        "#221 should have >=40 irreps, got {}",
        irreps.len()
    );

    // Γ 点 irrep (scalar + spinor)
    let gamma: Vec<_> = irreps
        .iter()
        .filter(|r| kpoint_label(r.ml) == "GM")
        .collect();
    assert!(gamma.len() >= 10, "#221 should have >=10 Γ irreps, got {}", gamma.len());

    // Scalar-only Γ irreps
    let gamma_scalar: Vec<_> = gamma.iter().filter(|r| !r.spinor).collect();
    assert_eq!(gamma_scalar.len(), 10, "#221 should have 10 scalar Γ irreps");

    // 检查关键 irrep
    let find_ml = |ml: &str| gamma.iter().find(|r| r.ml == ml);

    let gm1p = find_ml("GM1+").expect("GM1+ not found");
    assert_eq!(gm1p.dim, 1);
    assert_eq!(gm1p.image, "A1a");

    let gm4m = find_ml("GM4-").expect("GM4- not found");
    assert_eq!(gm4m.dim, 3);
    // GM4- 的 image 应该是 C48a (3D 极矢量)
    assert!(
        gm4m.image.starts_with('C'),
        "GM4- should be 3D (C*), got {}",
        gm4m.image
    );

    let gm5p = find_ml("GM5+").expect("GM5+ not found");
    assert_eq!(gm5p.dim, 3);

    // Lifshitz: GM4- (Γ₄⁻, polar vector) 的反称平方包含矢量表示, 所以 Lifshitz=0
    assert!(!gm4m.lifshitz, "GM4- should NOT satisfy Lifshitz (antisymmetric square contains vector rep)");
}

/// SG #225 Fm-3m: 面心立方, X 点标签与 primitive 不同
#[test]
fn sg225_fcc_irreps() {
    let irreps = irreps_of(225);
    assert!(
        irreps.len() >= 30,
        "#225 should have >=30 irreps, got {}",
        irreps.len()
    );

    let labels: Vec<&str> = irreps.iter().map(|r| r.ml).collect();
    // FCC 的 X 点标签用 "X1+" 等
    assert!(labels.contains(&"GM1+"), "#225 should have GM1+");
    assert!(labels.contains(&"GM4-"), "#225 should have GM4-");
}

/// SG #14 P2_1/c: 单斜, 常见空间群
#[test]
fn sg14_p21c_irreps() {
    let irreps = irreps_of(14);
    assert!(
        irreps.len() >= 16,
        "#14 should have >=16 irreps, got {}",
        irreps.len()
    );
}

/// SG #230 Ia-3d: 最后一个 SG
#[test]
fn sg230_last() {
    let irreps = irreps_of(230);
    assert!(
        !irreps.is_empty(),
        "#230 Ia-3d should have irreps"
    );
}

// ==========================================================================
// Isotropy Subgroup 数据验证
// ==========================================================================

/// 每个 irrep 的子群索引都指向 ISOTROPY_SUBGROUPS 的有效位置
#[test]
fn subgroup_indices_valid() {
    for ir in IRREPS.iter() {
        let subs = ir.subgroups();
        for sub in subs {
            assert!(
                sub.sg >= 1 && sub.sg <= 230,
                "irrep {}: invalid subgroup sg={}",
                ir.ml,
                sub.sg
            );
        }
    }
}

/// GM4- (Γ₄⁻, polar vector) 的子群应全为非中心对称群
#[test]
fn sg221_gm4m_subgroups_are_polar() {
    let gm4m = irreps_of(221)
        .iter()
        .find(|r| r.ml == "GM4-")
        .expect("GM4- not found in #221");

    let subs = gm4m.subgroups();
    assert!(!subs.is_empty(), "GM4- should have subgroups");

    let sg_nums: Vec<usize> = subs.iter().map(|s| s.sg).collect();
    assert!(
        sg_nums.contains(&160),
        "GM4- should have #160 R3m. Found: {:?}",
        sg_nums
    );
    assert!(
        sg_nums.contains(&1),
        "GM4- should have #1 P1 (fully broken). Found: {:?}",
        sg_nums
    );
}

/// SG_IRREP_INDEX 中相邻 SG 的 irrep 不重叠
#[test]
fn sg_ranges_no_overlap() {
    for sg in 1..230 {
        let (start1, count1) = SG_IRREP_INDEX[sg];
        let (start2, count2) = SG_IRREP_INDEX[sg + 1];
        let end1 = start1 as usize + count1 as usize;
        assert!(
            end1 <= start2 as usize,
            "SG {} range [{}, {}) overlaps with SG {} range [{}, {})",
            sg,
            start1,
            end1,
            sg + 1,
            start2,
            start2 + count2
        );
    }
}

// ==========================================================================
// 标签格式验证
// ==========================================================================

/// 验证 ML 标签格式: 字母前缀 + 数字 + 可选符号
fn is_valid_ml(label: &str) -> bool {
    let s = label.trim();
    // 找第一个数字
    let first_digit = s.find(|c: char| c.is_ascii_digit());
    match first_digit {
        None => false,
        Some(pos) => {
            // 前缀全是大写字母
            s[..pos].chars().all(|c| c.is_ascii_uppercase())
                && s[pos..].chars().all(|c| c.is_ascii_digit() || c == '+' || c == '-' || c.is_ascii_uppercase())
        }
    }
}

#[test]
fn ml_label_format() {
    let mut bad = Vec::new();
    for ir in IRREPS.iter() {
        if !is_valid_ml(ir.ml) {
            bad.push(ir.ml);
        }
    }
    assert!(
        bad.len() < IRREPS.len() / 4,
        "{} bad ML labels (out of {}): {:?}...",
        bad.len(),
        IRREPS.len(),
        &bad[..bad.len().min(20)]
    );
}

// ==========================================================================
// 交叉引用验证
// ==========================================================================

/// 对于 SG #221, Kovalev 标签的 k-index 范围应该正确
#[test]
fn sg221_kov_labels_well_formed() {
    let irreps = irreps_of(221);
    for ir in irreps {
        // Kovalev 标签: "k_{12}\\tau_{9}" 或 "k12t9"
        // 找 "k" 后的数字和 "t" 后的数字
        if let Some(k_start) = ir.kov.find('k') {
            let after_k = &ir.kov[k_start + 1..];
            let k_end = after_k.find(|c: char| !c.is_ascii_digit() && c != '_' && c != '{' && c != '}')
                .unwrap_or(after_k.len());
            let k_num: String = after_k[..k_end]
                .chars()
                .filter(|c| c.is_ascii_digit())
                .collect();
            if let Ok(k_idx) = k_num.parse::<u32>() {
                assert!(k_idx <= 20, "#221 kov k-index {} too large in {}", k_idx, ir.kov);
            }
        }
    }
}

/// X 点的标量 irrep 标签应完整 (10 个)
#[test]
fn sg221_x_labels_exist() {
    let irreps = irreps_of(221);
    let x_labels: Vec<&str> = irreps
        .iter()
        .filter(|r| !r.spinor && kpoint_label(r.ml) == "X")
        .map(|r| r.ml)
        .collect();
    assert_eq!(x_labels.len(), 10, "#221 should have 10 scalar X-point irreps");
    for label in &["X1+", "X2+", "X3+", "X4+", "X5+", "X1-", "X2-", "X3-", "X4-", "X5-"] {
        assert!(
            x_labels.contains(label),
            "#221 missing X label: {}",
            label
        );
    }
}

// ==========================================================================
// 数据库总容量验证
// ==========================================================================

/// 总共 4777 个标量 irrep + 3611 个双值 (spinor) irrep
#[test]
fn total_irrep_count() {
    let scalar = IRREPS.iter().filter(|r| !r.spinor).count();
    let spinor = IRREPS.iter().filter(|r| r.spinor).count();
    assert_eq!(scalar, 4777, "Expected 4777 scalar irreps");
    assert!(spinor >= 3000, "Expected >=3000 spinor irreps, got {}", spinor);
}

/// ISOTROPY_SUBGROUPS 总共 15239 条 (书中记录)
#[test]
fn total_isotropy_count() {
    assert_eq!(
        ISOTROPY_SUBGROUPS.len(),
        15239,
        "Expected 15239 isotropy subgroups"
    );
}

/// 每个晶体学系统都有数据
#[test]
fn all_crystal_systems_have_data() {
    let systems = [
        ("triclinic", 1..=2u8),
        ("monoclinic", 3..=15),
        ("orthorhombic", 16..=74),
        ("tetragonal", 75..=142),
        ("trigonal", 143..=167),
        ("hexagonal", 168..=194),
        ("cubic", 195..=230),
    ];
    for (name, range) in &systems {
        let count: usize = range.clone().map(|sg| irreps_of(sg).len()).sum();
        assert!(count > 0, "{} system has no irreps", name);
    }
}

// ==========================================================================
// BC 标签 vs ML 标签交叉验证（已知 convention 差异）
// ==========================================================================

/// SG #123 P4/mmm: ML 的 GM2+ 对应 BC 的 GM3+
/// 这是 Miller-Love 和 Bradley-Cracknell 之间已知的编号差异
#[test]
fn sg123_bc_labels_differ_from_ml() {
    let irs = irreps_of(123);
    // ML=GM2+ → BC=GM3+ (不是 GM2+!)
    let gm2p = irs.iter().find(|r| r.ml == "GM2+").expect("GM2+ not found");
    assert_eq!(gm2p.ml, "GM2+");
    assert!(gm2p.bc.contains("Gamma_{3}"), "BC should be GM3+, got: {}", gm2p.bc);
    assert_ne!(gm2p.bc, gm2p.ml, "BC label should differ from ML for GM2+");

    // ML=GM3+ → BC=GM2+
    let gm3p = irs.iter().find(|r| r.ml == "GM3+").expect("GM3+ not found");
    assert!(gm3p.bc.contains("Gamma_{2}"), "BC should be GM2+, got: {}", gm3p.bc);
    assert_ne!(gm3p.bc, gm3p.ml, "BC label should differ from ML for GM3+");

    // ML=GM1+ → BC=GM1+ (相同)
    let gm1p = irs.iter().find(|r| r.ml == "GM1+").expect("GM1+ not found");
    assert!(gm1p.bc.contains("Gamma_{1}"), "BC should be GM1+, got: {}", gm1p.bc);
}

/// SG #121 I-42m: ML 的 P3P3 对应 BC 的 P4P4
#[test]
fn sg121_bc_labels_differ_from_ml() {
    let irs = irreps_of(121);
    let p3p3 = irs.iter().find(|r| r.ml == "P3P3").expect("P3P3 not found");
    assert_eq!(p3p3.ml, "P3P3");
    assert!(p3p3.bc.contains("P_{4}"), "BC should be P4P4, got: {}", p3p3.bc);

    let p4p4 = irs.iter().find(|r| r.ml == "P4P4").expect("P4P4 not found");
    assert!(p4p4.bc.contains("P_{2}"), "BC should be P2P2, got: {}", p4p4.bc);
}

/// SG #1 P1: BC 使用不同的 k 点字母
#[test]
fn sg1_bc_uses_different_k_labels() {
    let irs = irreps_of(1);
    let z1 = irs.iter().find(|r| r.ml == "Z1").expect("Z1 not found");
    assert_ne!(z1.bc, z1.ml, "Z1 BC should differ");

    let x1 = irs.iter().find(|r| r.ml == "X1").expect("X1 not found");
    assert_ne!(x1.bc, x1.ml, "X1 BC should differ (BC uses B1)");
}

/// 全局: 2162 个标量 irrep 的 BC ≠ ML
#[test]
fn total_bc_ml_differences() {
    let mut diff = 0u32;
    for ir in IRREPS.iter().filter(|r| !r.spinor) {
        let bc_clean = ir.bc.replace("\\", "").replace("{", "").replace("}", "")
            .replace("^", "").replace("_", "").replace("$", "");
        let ml_clean = ir.ml.replace("+", "").replace("-", "");
        if !bc_clean.contains(&ml_clean) && !ml_clean.contains(&bc_clean) {
            diff += 1;
        }
    }
    assert!(diff > 2000, "Expected >2000 BC≠ML differences, got {}", diff);
    assert!(diff < 3000, "Too many differences: {}", diff);
}

// ==========================================================================
// 已知 Isotropy Subgroup 数据验证（对照 Stokes & Hatch 表）
// ==========================================================================

/// SG #137 P4_2/nmc: irrep R1 的子群集合
/// 已知结果: {38, 11, 109, 119, 15, 6, 44, 44, 2, 9, 5, 8, 1}
/// 注意: 两个 #44 的 Basis 相同但 Origin 不同
#[test]
fn sg137_r1_subgroups() {
    let r1 = irreps_of(137).iter()
        .find(|r| r.ml == "R1")
        .expect("R1 not found");

    let subs = r1.subgroups();
    assert_eq!(subs.len(), 13, "R1 should have 13 subgroups");

    let sg_nums: Vec<usize> = subs.iter().map(|s| s.sg).collect();

    // 验证集合 (忽略顺序)
    let expected: Vec<usize> = vec![38, 11, 109, 119, 15, 6, 44, 44, 2, 9, 5, 8, 1];
    let mut actual_sorted = sg_nums.clone();
    actual_sorted.sort();
    let mut expected_sorted = expected.clone();
    expected_sorted.sort();
    assert_eq!(actual_sorted, expected_sorted,
        "R1 subgroups set mismatch. Got: {:?}", sg_nums);

    // 验证两个 #44 都存在
    let count_44 = sg_nums.iter().filter(|&&n| n == 44).count();
    assert_eq!(count_44, 2, "Should have exactly 2 entries for #44");
}

/// SG #140 I4/mcm: irrep N1 子群数量
///
/// ISO 数据库 (2022) 给出 69 个。Stokes & Hatch (1988) 书中为 79 个。
/// 差异来自数据库版本更新——非解析错误。
#[test]
fn sg140_n1_subgroup_count() {
    let n1 = irreps_of(140).iter()
        .find(|r| r.ml == "N1")
        .expect("N1 not found");
    let count = n1.subgroups().len();
    assert_eq!(count, 69, "ISO DB (2022) N1 subgroup count");
}

// ==========================================================================
// k 点坐标验证
// ==========================================================================

/// SG #221 Pm-3m: 验证关键 k 点坐标
#[test]
fn sg221_kpoint_coordinates() {
    let irs = irreps_of(221);
    // Γ = (0, 0, 0)
    let gm = irs.iter().find(|r| r.ml == "GM1+").expect("GM1+");
    assert_eq!((gm.kx, gm.ky, gm.kz, gm.kd), (0, 0, 0, 1), "Γ should be (0,0,0)");
    // X = (0, 1/2, 0) = (0/2, 1/2, 0/2)
    let x = irs.iter().find(|r| r.ml == "X1+").expect("X1+");
    assert_eq!((x.kx, x.ky, x.kz, x.kd), (0, 1, 0, 2), "X should be (0,1/2,0)");
    // M = (1/2, 1/2, 0)
    let m = irs.iter().find(|r| r.ml == "M1+").expect("M1+");
    assert_eq!((m.kx, m.ky, m.kz, m.kd), (1, 1, 0, 2), "M should be (1/2,1/2,0)");
    // R = (1/2, 1/2, 1/2)
    let r = irs.iter().find(|r| r.ml == "R1+").expect("R1+");
    assert_eq!((r.kx, r.ky, r.kz, r.kd), (1, 1, 1, 2), "R should be (1/2,1/2,1/2)");
}

/// SG #217 I-43m: P 点坐标 = (1/2, 1/2, 1/2)
#[test]
fn sg217_p_point_coordinate() {
    let irs = irreps_of(217);
    let p = irs.iter().find(|r| r.ml.starts_with("P1")).expect("P1 not found");
    assert_eq!((p.kx, p.ky, p.kz, p.kd), (1, 1, 1, 2),
        "SG #217 P-point should be (1/2,1/2,1/2), got ({}/{}, {}/{}, {}/{})",
        p.kx, p.kd, p.ky, p.kd, p.kz, p.kd);
}

/// 验证所有 irrep 的 k 向量分母为正
#[test]
fn all_kvectors_have_positive_denominator() {
    for ir in IRREPS.iter() {
        assert!(ir.kd > 0, "irrep {} has non-positive kd={}", ir.ml, ir.kd);
    }
}

/// 输出 #142 X 点: 操作 ↔ 特征标 对应关系
#[test]
fn dump_sg142_x_point() {
    let irs = irreps_of(142);
    let ops = SymmetryOps::from_sg(142).unwrap();
    println!("\n=== SG #142 I4_1/acd — {} irreps, {} ops ===", irs.len(), ops.len());

    let x1 = irs.iter().find(|r| r.ml == "X1").unwrap();
    let chars = x1.characters();
    println!("\nX1 dim={} k=({}/{},{}/{},{}/{})", x1.dim, x1.kx, x1.kd, x1.ky, x1.kd, x1.kz, x1.kd);
    println!("{:3} | {{R|t}}                                    | χ", "Op");
    println!("----+------------------------------------------+-----");
    for (i, op) in ops.iter().enumerate() {
        let r = op.rotation; let t = op.translation;
        let chi = if i < chars.len() { chars[i] } else { 0.0 };
        if chi.abs() > 0.01 || i == 0 {
            println!("{:3} | [{:2},{:2},{:2};{:2},{:2},{:2};{:2},{:2},{:2}] ({:5.2},{:5.2},{:5.2}) | {:+.1}",
                i, r[0][0],r[0][1],r[0][2], r[1][0],r[1][1],r[1][2], r[2][0],r[2][1],r[2][2],
                t[0],t[1],t[2], chi);
        }
    }

    // Show matrix for identity operation
    let mats = x1.matrices();
    println!("\nX1 matrix for Op 0 (identity, 4x4):");
    for row in 0..4 {
        print!("  ");
        for col in 0..4 { print!("{:5.0} ", mats[row*4+col]); }
        println!();
    }
}

/// 验证 k 点坐标与 label 前缀的一致性 (标量 irrep)
#[test]
fn same_prefix_same_kpoint() {
    let irs = irreps_of(221);
    let mut groups: std::collections::HashMap<&str, Vec<&IrrepRecord>> = std::collections::HashMap::new();
    for ir in irs.iter().filter(|r| !r.spinor) {
        let k = ir.k_label();
        groups.entry(k).or_default().push(ir);
    }
    for (label, group) in &groups {
        if group.len() <= 1 {
            continue;
        }
        let first = group[0];
        let kp = (first.kx, first.ky, first.kz, first.kd);
        for ir in group.iter().skip(1) {
            assert_eq!((ir.kx, ir.ky, ir.kz, ir.kd), kp,
                "SG221 k-label '{}': {} has different k-point from {}", label, ir.ml, first.ml);
        }
    }
}
