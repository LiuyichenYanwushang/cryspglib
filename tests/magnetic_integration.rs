//! 磁性空间群集成测试。
//!
//! 所有测试走公共 API `Crystal` + `SymmetryAnalysis`，覆盖 Type-1/2/3/4 真实物理系统。

use cryspglib::{Crystal, MagneticSpaceGroupType, MagneticType, SymError, SymmetryOps};

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
    let mut cry = Crystal::new(*lattice, positions.to_vec(), types.to_vec()).unwrap();
    if let Some(m) = moments {
        cry = cry.with_magnetic(m.to_vec()).unwrap();
    }
    let result = cry
        .analyze()
        .symprec(SYMPREC)
        .magnetic_dataset()
        .unwrap_or_else(|e| panic!("{}: magnetic_dataset failed: {:?}", label, e));
    eprintln!("=== {} ===", label);
    eprintln!("{}", result);
    result
}

// ====================================================================
// Public Rust API: MagneticSpaceGroupType::classify
// ====================================================================

fn pm3m_ops() -> (Vec<[[i32; 3]; 3]>, Vec<[f64; 3]>) {
    let operations = SymmetryOps::from_database(517).unwrap();
    let rots = operations
        .operations
        .iter()
        .map(|operation| operation.rotation)
        .collect();
    let trans = operations
        .operations
        .iter()
        .map(|operation| operation.translation)
        .collect();
    (rots, trans)
}

#[test]
fn test_api_type1() {
    let (rots, trans) = pm3m_ops();
    let result =
        MagneticSpaceGroupType::classify(&rots, &trans, None, &cubic_lattice(), SYMPREC).unwrap();
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
        &all_rots,
        &all_trans,
        Some(&timerev),
        &cubic_lattice(),
        SYMPREC,
    )
    .unwrap();
    assert_eq!(result.type_, MagneticType::Grey);
    assert!(result.uni_number > 0);
}

#[test]
fn test_api_type3() {
    let (rots, trans) = pm3m_ops();
    let timerev: Vec<bool> = rots
        .iter()
        .map(|r| !cryspglib::mathfunc::is_proper(r))
        .collect();
    let result =
        MagneticSpaceGroupType::classify(&rots, &trans, Some(&timerev), &cubic_lattice(), SYMPREC)
            .unwrap();
    assert_eq!(result.type_, MagneticType::BlackWhite);
    assert!(result.uni_number > 0);
}

#[test]
fn test_api_reports_operation_only_ambiguity() {
    let ops = SymmetryOps::from_magnetic_database(282).unwrap();
    let rotations: Vec<_> = ops.operations.iter().map(|op| op.rotation).collect();
    let translations: Vec<_> = ops.operations.iter().map(|op| op.translation).collect();
    let time_reversals: Vec<_> = ops.operations.iter().map(|op| op.time_reversal).collect();
    let lattice = [[1.0, 0.0, 0.0], [0.0, 1.3, 0.0], [0.0, 0.0, 1.7]];

    let result = MagneticSpaceGroupType::classify(
        &rotations,
        &translations,
        Some(&time_reversals),
        &lattice,
        SYMPREC,
    );

    assert!(matches!(result, Err(SymError::MagneticUniAmbiguous)));
}

#[test]
fn test_api_rejects_mismatched_operation_lengths() {
    let result = MagneticSpaceGroupType::classify(
        &[[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
        &[],
        None,
        &cubic_lattice(),
        SYMPREC,
    );

    assert!(matches!(result, Err(SymError::InvalidInput)));
}

#[test]
fn test_from_uni_rejects_invalid_identifiers() {
    for uni in [0, 1652, usize::MAX] {
        assert!(matches!(
            MagneticSpaceGroupType::from_uni(uni),
            Err(SymError::InvalidInput)
        ));
    }
}

#[test]
fn test_from_uni_accepts_database_boundaries() {
    let first = MagneticSpaceGroupType::from_uni(1).unwrap();
    assert_eq!(first.uni_number, 1);
    assert_eq!(first.bns_number.trim(), "1.1");
    assert_eq!(first.number, 1);
    assert_eq!(first.type_, MagneticType::Ordinary);

    let last = MagneticSpaceGroupType::from_uni(1651).unwrap();
    assert_eq!(last.uni_number, 1651);
    assert_eq!(last.bns_number.trim(), "230.149");
    assert_eq!(last.number, 230);
    assert_eq!(last.type_, MagneticType::BlackWhite);
}

#[test]
fn test_magnetic_dataset_rejects_invalid_crystal_inputs() {
    let lattice = cubic_lattice();
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];

    assert!(matches!(
        Crystal::new(lattice, positions.to_vec(), vec![26]),
        Err(SymError::InvalidInput)
    ));
    assert!(matches!(
        Crystal::new(lattice, positions.to_vec(), vec![26, 26])
            .unwrap()
            .with_magnetic(vec![[1.0, 0.0, 0.0]]),
        Err(SymError::InvalidInput)
    ));

    let empty = Crystal::new(lattice, vec![], vec![]).unwrap();
    assert!(matches!(
        empty.analyze().symprec(SYMPREC).magnetic_dataset(),
        Err(SymError::InvalidInput)
    ));
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
    assert_eq!(
        r.uni_number, 1005,
        "[100] and [001] must match the same UNI"
    );
    assert_eq!(r.bns_number.trim(), "123.345");
    assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
}

/// Fe BCC AFM [111]: 2 个 Fe 在 [0,0,0] 和 [0.5,0.5,0.5], 磁矩相反沿 [111]
/// 反幺正体心平移 (I|1/2,1/2,1/2)' 使其成为 type-4:
/// R-3c (#167), BNS=167.108, UNI=1338。
#[test]
fn test_fe_bcc_afm_111() {
    let lattice = cubic_lattice();
    let n = (3.0_f64).sqrt();
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
    let types = [26, 26];
    let moments = [[1.0 / n, 1.0 / n, 1.0 / n], [-1.0 / n, -1.0 / n, -1.0 / n]];
    let r = run_dataset(
        "Fe BCC AFM [111]",
        &lattice,
        &positions,
        &types,
        Some(&moments),
    );
    assert_eq!(r.spacegroup_number, 229, "non-mag: Im-3m");
    assert!(r.uni_number > 0, "AFM [111] must match a DB entry");
    assert_eq!(r.magnetic_type, MagneticType::AntiTranslation);
    assert_eq!(r.bns_number.trim(), "167.108");
    assert_eq!(r.uni_number, 1338);
}

/// FCC FM [001]: 4 个原子, 全部磁矩沿 [001]
/// FCC 中心化在四方标准 setting 中变为 I-centered:
/// I4/mmm (#139) type-3, BNS=139.537, UNI=1197。
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
    assert_eq!(r.uni_number, 1197);
    assert_eq!(r.bns_number.trim(), "139.537");
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

/// Graphene honeycomb lattice, 2 atoms per cell, AFM (opposite z-spins).
///
/// Lattice: hexagonal a1=(1,0,0), a2=(0.5,√3/2,0), a3=(0,0,10)
/// Atom A at (0,0,0) with spin +z, atom B at (1/3,2/3,0) with spin -z
///
/// Non-mag space group: P6/mmm (#191).
/// Magnetic space group: UNI 1466, BNS 191.236.
#[test]
fn test_graphene_afm_z() {
    let s3 = (3.0_f64).sqrt();
    // lattice[cart][vec]: rows=x/y/z, cols=a/b/c
    // a=(1,0,0), b=(1/2,√3/2,0), c=(0,0,2) → hexagonal #191
    let lattice = [
        [1.0, 0.5, 0.0],      // row x: a_x, b_x, c_x
        [0.0, s3 / 2.0, 0.0], // row y: a_y, b_y, c_y
        [0.0, 0.0, 2.0],      // row z: a_z, b_z, c_z
    ];
    // Sublattice A at origin, sublattice B at nearest-neighbor position
    // τ_B = (1/3, 1/3, 0) in fractional coords → C-C bond along δ₁ direction
    let positions = [[0.0, 0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0, 0.0]];
    let types = [6, 6];
    // AFM: opposite z-spins on A and B sublattices
    let moments = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
    // The identification must be stable across reasonable tolerances.
    for &sp in &[1e-3, 1e-4, 1e-5, 1e-6] {
        let cry = Crystal::new(lattice, positions.to_vec(), types.to_vec())
            .unwrap()
            .with_magnetic(moments.to_vec())
            .unwrap();
        let r = cry
            .analyze()
            .symprec(sp)
            .magnetic_dataset()
            .unwrap_or_else(|e| panic!("symprec={sp}: magnetic_dataset failed: {e:?}"));

        assert_eq!(r.spacegroup_number, 191, "symprec={sp}");
        assert_eq!(r.hall_number, 485, "symprec={sp}");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite, "symprec={sp}");
        assert_eq!(r.uni_number, 1466, "symprec={sp}");
        assert_eq!(r.bns_number.trim(), "191.236", "symprec={sp}");
        assert_eq!(r.num_operations, 24, "symprec={sp}");
        assert_eq!(
            r.time_reversals.iter().filter(|&&tr| !tr).count(),
            12,
            "symprec={sp}: unitary operation count",
        );
        assert_eq!(
            r.time_reversals.iter().filter(|&&tr| tr).count(),
            12,
            "symprec={sp}: anti-unitary operation count",
        );
    }
}

/// Bilayer graphene: two sublattices at slightly different z-heights (0.51 vs 0.49).
///
/// This breaks the perfect mirror symmetry of planar graphene (#191 → lower symmetry).
/// Three magnetic configurations are tested side-by-side to verify consistency.
#[test]
fn test_graphene_bilayer_z() {
    let s3 = (3.0_f64).sqrt();
    let lattice = [[1.0, 0.5, 0.0], [0.0, s3 / 2.0, 0.0], [0.0, 0.0, 2.0]];
    // Atom A at z=0.51, atom B at z=0.49 (fractional coords)
    let positions = [[0.0, 0.0, 0.51], [1.0 / 3.0, 1.0 / 3.0, 0.49]];
    let types = [6, 6];

    // --- Non-magnetic ---
    // Broken z-mirror: P6/mmm (#191) → P-3m1 (#164), 24→12 ops
    {
        let cry = Crystal::new(lattice, positions.to_vec(), types.to_vec()).unwrap();
        let r = cry
            .analyze()
            .symprec(1e-5)
            .magnetic_dataset()
            .unwrap_or_else(|e| panic!("bilayer non-mag: {e:?}"));
        assert_eq!(r.spacegroup_number, 164, "non-mag: P-3m1");
        assert_eq!(r.hall_number, 456);
        assert_eq!(r.num_operations, 12);
        assert_eq!(r.magnetic_type, MagneticType::NonMagnetic);
    }

    // --- FM: both moments along +z ---
    // Type-3 BlackWhite, UNI=1319, BNS=164.89
    {
        let moments = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]];
        let cry = Crystal::new(lattice, positions.to_vec(), types.to_vec())
            .unwrap()
            .with_magnetic(moments.to_vec())
            .unwrap();
        let r = cry
            .analyze()
            .symprec(1e-5)
            .magnetic_dataset()
            .unwrap_or_else(|e| panic!("bilayer FM: {e:?}"));
        assert_eq!(r.spacegroup_number, 164);
        assert_eq!(r.hall_number, 456);
        assert_eq!(r.uni_number, 1319);
        assert_eq!(r.bns_number.trim(), "164.89");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        assert_eq!(r.num_operations, 12);
    }

    // --- AFM: one +z, one -z ---
    // Type-3 BlackWhite, UNI=1318, BNS=164.88
    {
        let moments = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
        let cry = Crystal::new(lattice, positions.to_vec(), types.to_vec())
            .unwrap()
            .with_magnetic(moments.to_vec())
            .unwrap();
        let r = cry
            .analyze()
            .symprec(1e-5)
            .magnetic_dataset()
            .unwrap_or_else(|e| panic!("bilayer AFM: {e:?}"));
        assert_eq!(r.spacegroup_number, 164);
        assert_eq!(r.hall_number, 456);
        assert_eq!(r.uni_number, 1318);
        assert_eq!(r.bns_number.trim(), "164.88");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        assert_eq!(r.num_operations, 12);
    }
}

// ====================================================================
// MnAl2O4 spinel: non-magnetic Fd-3m (#227); Mn AFM along z → BNS 141.556
// ====================================================================
//
// Spinel structure from VASP POSCAR (`MnAl2O4`). Lattice is cubic with
// a = 8.2097997665 Å.  Direct coordinates: 8 Mn (8a), 16 Al (16d), 32 O (32e).
//
// Mn moments are collinear along z with the pattern −, +, −, +, −, +, −, +
// in POSCAR order, i.e. the Mn atoms at z = 0.25 and z = 0.75 carry −z
// moments and all remaining Mn carry +z.
//
// Expected results:
//   non-magnetic: SG 227 (Fd-3m);
//   magnetic:     UNI 1216, BNS 141.556 (Type-3 BlackWhite, parent I41/amd #141).

#[test]
fn test_mnal2o4_afm_z() {
    let lattice = [
        [8.2097997665, 0.0, 0.0],
        [0.0, 8.2097997665, 0.0],
        [0.0, 0.0, 8.2097997665],
    ];
    // 8 Mn + 16 Al + 32 O, direct coordinates from the VASP POSCAR.
    let positions = [
        // Mn
        [0.000000000, 0.000000000, 0.000000000],
        [0.250000000, 0.250000000, 0.750000000],
        [0.750000000, 0.250000000, 0.250000000],
        [0.500000000, 0.500000000, 0.000000000],
        [0.500000000, 0.000000000, 0.500000000],
        [0.250000000, 0.750000000, 0.250000000],
        [0.000000000, 0.500000000, 0.500000000],
        [0.750000000, 0.750000000, 0.750000000],
        // Al
        [0.375000000, 0.375000000, 0.375000000],
        [0.625000000, 0.875000000, 0.125000000],
        [0.875000000, 0.625000000, 0.125000000],
        [0.125000000, 0.625000000, 0.875000000],
        [0.375000000, 0.875000000, 0.875000000],
        [0.125000000, 0.125000000, 0.375000000],
        [0.125000000, 0.375000000, 0.125000000],
        [0.875000000, 0.375000000, 0.875000000],
        [0.375000000, 0.125000000, 0.125000000],
        [0.625000000, 0.125000000, 0.875000000],
        [0.625000000, 0.375000000, 0.625000000],
        [0.375000000, 0.625000000, 0.625000000],
        [0.625000000, 0.625000000, 0.375000000],
        [0.875000000, 0.125000000, 0.625000000],
        [0.875000000, 0.875000000, 0.375000000],
        [0.125000000, 0.875000000, 0.625000000],
        // O
        [0.140900001, 0.140900001, 0.140900001],
        [0.390900016, 0.109099999, 0.890900016],
        [0.109099999, 0.390900016, 0.890900016],
        [0.890900016, 0.390900016, 0.109099999],
        [0.609099984, 0.109099999, 0.109099999],
        [0.359099984, 0.359099984, 0.140900001],
        [0.640900016, 0.859099984, 0.359099984],
        [0.359099984, 0.140900001, 0.359099984],
        [0.109099999, 0.609099984, 0.109099999],
        [0.859099984, 0.640900016, 0.359099984],
        [0.140900001, 0.359099984, 0.359099984],
        [0.390900016, 0.890900016, 0.109099999],
        [0.609099984, 0.890900016, 0.890900016],
        [0.890900016, 0.609099984, 0.890900016],
        [0.640900016, 0.140900001, 0.640900016],
        [0.390900016, 0.609099984, 0.390900016],
        [0.609099984, 0.390900016, 0.390900016],
        [0.390900016, 0.390900016, 0.609099984],
        [0.609099984, 0.609099984, 0.609099984],
        [0.359099984, 0.859099984, 0.640900016],
        [0.640900016, 0.359099984, 0.859099984],
        [0.859099984, 0.140900001, 0.859099984],
        [0.859099984, 0.359099984, 0.640900016],
        [0.140900001, 0.640900016, 0.640900016],
        [0.640900016, 0.640900016, 0.140900001],
        [0.859099984, 0.859099984, 0.140900001],
        [0.890900016, 0.109099999, 0.390900016],
        [0.109099999, 0.890900016, 0.390900016],
        [0.890900016, 0.890900016, 0.609099984],
        [0.109099999, 0.109099999, 0.609099984],
        [0.140900001, 0.859099984, 0.859099984],
        [0.359099984, 0.640900016, 0.859099984],
    ];
    let types = [
        25, 25, 25, 25, 25, 25, 25, 25, // Mn
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, // Al
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, // O
    ];

    // --- Non-magnetic: spinel Fd-3m (#227) ---
    {
        let cry = Crystal::new(lattice, positions.to_vec(), types.to_vec()).unwrap();
        let r = cry
            .analyze()
            .symprec(SYMPREC)
            .dataset()
            .unwrap_or_else(|e| panic!("MnAl2O4 non-mag dataset failed: {e:?}"));
        assert_eq!(r.spacegroup_number, 227, "non-magnetic space group");
    }

    // --- AFM along z: UNI 1216 / BNS 141.556 ---
    {
        // Mn at z = 0.25 or z = 0.75 → −z; all other Mn → +z.
        let mut moments = vec![[0.0, 0.0, 0.0]; positions.len()];
        let mn_signs = [1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0];
        for (i, sign) in mn_signs.into_iter().enumerate() {
            moments[i] = [0.0, 0.0, sign];
        }

        let cry = Crystal::new(lattice, positions.to_vec(), types.to_vec())
            .unwrap()
            .with_magnetic(moments)
            .unwrap();
        let r = cry
            .analyze()
            .symprec(SYMPREC)
            .magnetic_dataset()
            .unwrap_or_else(|e| panic!("MnAl2O4 magnetic_dataset failed: {e:?}"));

        // `spacegroup_number` is the *structural* parent group detected before
        // the moments are considered (Fd-3m #227).  The magnetic group itself
        // is identified by UNI/BNS below (parent I41/amd #141).
        assert_eq!(r.spacegroup_number, 227);
        assert_eq!(r.hall_number, 525, "Fd-3m structural setting");
        assert_eq!(r.magnetic_type, MagneticType::BlackWhite);
        assert_eq!(r.uni_number, 1216);
        assert_eq!(r.bns_number.trim(), "141.556");
        // Type-3 black-white: 32 spatial operations of the conventional cell
        // and their time-reversed partners.
        assert_eq!(r.num_operations, 64);
        assert_eq!(
            r.time_reversals.iter().filter(|&&tr| !tr).count(),
            32,
            "unitary operation count",
        );
    }
}
