//! 磁性空间群集成测试。
//!
//! 所有测试走公共 API `Crystal` + `SymmetryAnalysis`，覆盖 Type-1/2/3 真实物理系统。

use cryspglib::{
    Crystal, MagneticSpaceGroupType, MagneticType,
};

const SYMPREC: f64 = 1e-5;

fn cubic_lattice() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn run_dataset(
    label: &str,
    lattice: &[[f64; 3]; 3],
    positions: &[[f64; 3]],
    types: &[i32],
    moments: Option<&[[f64; 3]]>,
) -> cryspglib::MagneticSymmetry {
    let mut cry = Crystal::new(*lattice, positions.to_vec(), types.to_vec());
    if let Some(m) = moments {
        cry = cry.with_magnetic(m.to_vec());
    }
    let result = cry.analyze().symprec(SYMPREC).magnetic_dataset()
        .unwrap_or_else(|| panic!("{}: magnetic_dataset returned None", label));
    eprintln!("=== {} ===", label);
    eprintln!("{}", result.to_string());
    result
}

// ====================================================================
// 公共 API: spg_get_magnetic_spacegroup_type_from_symmetry
// ====================================================================

fn pm3m_ops() -> (Vec<[[i32; 3]; 3]>, Vec<[f64; 3]>) {
    let (count, start) = cryspglib::spg_database::spgdb_get_operation_index(517);
    let mut rots = Vec::with_capacity(count);
    let mut trans = Vec::with_capacity(count);
    for i in 0..count {
        if let Some((r, t)) = cryspglib::spg_database::spgdb_get_operation_by_index(start + i) {
            rots.push(r);
            trans.push(t);
        }
    }
    (rots, trans)
}

#[test]
fn test_api_type1() {
    let (rots, trans) = pm3m_ops();
    let result = MagneticSpaceGroupType::classify(
        &rots, &trans, None, &cubic_lattice(), SYMPREC,
    );
    assert_eq!(result.type_, MagneticType::Ordinary);
    assert!(result.uni_number > 0);
}

#[test]
fn test_api_type2() {
    let (rots, trans) = pm3m_ops();
    let n = rots.len();
    let all_rots: Vec<_> = rots.iter().chain(rots.iter()).cloned().collect();
    let all_trans: Vec<_> = trans.iter().chain(trans.iter()).cloned().collect();
    let timerev: Vec<bool> = (0..n).map(|_| false).chain((0..n).map(|_| true)).collect();
    let result = MagneticSpaceGroupType::classify(
        &all_rots, &all_trans, Some(&timerev), &cubic_lattice(), SYMPREC,
    );
    assert_eq!(result.type_, MagneticType::Grey);
    assert!(result.uni_number > 0);
}

#[test]
fn test_api_type3() {
    let (rots, trans) = pm3m_ops();
    let timerev: Vec<bool> = rots.iter().map(|r| !cryspglib::mathfunc::is_proper(r)).collect();
    let result = MagneticSpaceGroupType::classify(
        &rots, &trans, Some(&timerev), &cubic_lattice(), SYMPREC,
    );
    assert_eq!(result.type_, MagneticType::BlackWhite);
    assert!(result.uni_number > 0);
}

// ====================================================================
// 物理系统测试 — 全部使用 Crystal API
// ====================================================================

/// Fe SC 体心, 非磁 (moments=None)
#[test]
fn test_fe_sc_nonmagnetic() {
    let lattice = cubic_lattice();
    let positions = [[0.5, 0.5, 0.5]];
    let types = [26];
    let r = run_dataset("Fe SC, non-magnetic", &lattice, &positions, &types, None);
    assert_eq!(r.spacegroup_number, 221);
    assert_eq!(r.magnetic_type, MagneticType::NonMagnetic);
    assert_eq!(r.uni_number, 0);
    assert!(r.num_operations > 0);
}

/// Fe SC 体心, 磁矩 [001] → P4/mmm (#123) BNS=123.345, UNI=1005
#[test]
fn test_fe_sc_001() {
    let lattice = cubic_lattice();
    let positions = [[0.5, 0.5, 0.5]];
    let types = [26];
    let moments = [[0.0, 0.0, 1.0]];
    let r = run_dataset("Fe SC [001]", &lattice, &positions, &types, Some(&moments));
    assert_eq!(r.spacegroup_number, 221, "non-mag: Pm-3m");
    assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
    assert!(r.uni_number > 0, "must match DB entry, not fallback");
    assert!(!r.bns_number.is_empty(), "must have BNS");
    assert_eq!(r.bns_number.trim(), "123.345");
    assert_eq!(r.uni_number, 1005);
}

/// Fe SC 体心, 磁矩 [100] → 应与 [001] 一样匹配 UNI=1005
#[test]
fn test_fe_sc_100() {
    let lattice = cubic_lattice();
    let positions = [[0.5, 0.5, 0.5]];
    let types = [26];
    let moments = [[1.0, 0.0, 0.0]];
    let r = run_dataset("Fe SC [100]", &lattice, &positions, &types, Some(&moments));
    assert_eq!(r.spacegroup_number, 221, "non-mag: Pm-3m");
    assert_eq!(r.uni_number, 1005, "[100] and [001] must match the same UNI");
    assert_eq!(r.bns_number.trim(), "123.345");
    assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
}

/// Fe BCC AFM [111]: 2 个 Fe 在 [0,0,0] 和 [0.5,0.5,0.5], 磁矩相反沿 [111]
/// 预期: R-3m (#166) type-3, BNS=166.101, UNI=1331
#[test]
fn test_fe_bcc_afm_111() {
    let lattice = cubic_lattice();
    let n = (3.0_f64).sqrt();
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
    let types = [26, 26];
    let moments = [
        [1.0 / n, 1.0 / n, 1.0 / n],
        [-1.0 / n, -1.0 / n, -1.0 / n],
    ];
    let r = run_dataset("Fe BCC AFM [111]", &lattice, &positions, &types, Some(&moments));
    assert_eq!(r.spacegroup_number, 229, "non-mag: Im-3m");
    assert!(r.uni_number > 0, "AFM [111] must match a DB entry");
    assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
    assert_eq!(r.bns_number.trim(), "166.101");
    assert_eq!(r.uni_number, 1331);
}

/// FCC FM [001]: 4 个原子, 全部磁矩沿 [001]
#[test]
fn test_fcc_fm_001() {
    let lattice = cubic_lattice();
    let positions = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let types = [26, 26, 26, 26];
    let moments = [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
    ];
    let r = run_dataset("FCC FM [001]", &lattice, &positions, &types, Some(&moments));
    assert_eq!(r.spacegroup_number, 225, "non-mag: Fm-3m");
    assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
    assert_eq!(r.uni_number, 1005);
    assert_eq!(r.bns_number.trim(), "123.345");
}

/// FCC FM [111]: 4 个原子, 全部磁矩沿 [111]
/// 预期: R-3m type-3, BNS=166.101, UNI=1331
#[test]
fn test_fcc_fm_111() {
    let lattice = cubic_lattice();
    let n = (3.0_f64).sqrt();
    let positions = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let types = [26, 26, 26, 26];
    let m = [1.0 / n, 1.0 / n, 1.0 / n];
    let moments = [m, m, m, m];
    let r = run_dataset("FCC FM [111]", &lattice, &positions, &types, Some(&moments));
    assert_eq!(r.spacegroup_number, 225, "non-mag: Fm-3m");
    assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
    assert!(r.uni_number > 0, "FCC FM [111] must match DB entry");
    assert_eq!(r.bns_number.trim(), "166.101");
    assert_eq!(r.uni_number, 1331);
}
