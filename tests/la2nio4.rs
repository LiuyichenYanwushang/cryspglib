//! La₂NiO₄ 空间群集成测试。
//!
//! 四方 P4₂/ncm (#138)，D₄ₕ 点群。晶胞含 8 La + 4 Ni + 16 O = 28 原子。

use cryspglib::{
    Crystal,
};

const SYMPREC: f64 = 1e-5;

/// La₂NiO₄ 四方 P4₂/ncm (#138)，D₄ₕ 点群，16 个对称操作。
#[test]
fn test_la2nio4() {
    let lattice: [[f64; 3]; 3] = [
        [5.4994997978, 0.0, 0.0],
        [0.0, 5.4994997978, 0.0],
        [0.0, 0.0, 12.5052003860],
    ];

    let positions: [[f64; 3]; 28] = [
        // La (8 atoms)
        [0.992799997, 0.992799997, 0.363900006],
        [0.492799997, 0.007200000, 0.136099994],
        [0.007200000, 0.492799997, 0.136099994],
        [0.507200003, 0.507200003, 0.363900006],
        [0.007200000, 0.007200000, 0.636099994],
        [0.507200003, 0.992799997, 0.863900006],
        [0.992799997, 0.507200003, 0.863900006],
        [0.492799997, 0.492799997, 0.636099994],
        // Ni (4 atoms)
        [0.000000000, 0.000000000, 0.000000000],
        [0.500000000, 0.000000000, 0.500000000],
        [0.000000000, 0.500000000, 0.500000000],
        [0.500000000, 0.500000000, 0.000000000],
        // O (16 atoms)
        [0.750000000, 0.750000000, 0.015500000],
        [0.250000000, 0.250000000, 0.484499991],
        [0.250000000, 0.250000000, 0.984499991],
        [0.750000000, 0.750000000, 0.515500009],
        [0.750000000, 0.250000000, 0.000000000],
        [0.250000000, 0.750000000, 0.500000000],
        [0.250000000, 0.750000000, 0.000000000],
        [0.750000000, 0.250000000, 0.500000000],
        [0.031399999, 0.031399999, 0.177100003],
        [0.531400025, 0.968599975, 0.322899997],
        [0.968599975, 0.531400025, 0.322899997],
        [0.468600005, 0.468600005, 0.177100003],
        [0.968599975, 0.968599975, 0.822899997],
        [0.468600005, 0.031399999, 0.677100003],
        [0.031399999, 0.468600005, 0.677100003],
        [0.531400025, 0.531400025, 0.822899997],
    ];

    let types: [i32; 28] = [
        57, 57, 57, 57, 57, 57, 57, 57, // La
        28, 28, 28, 28, // Ni
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // O
    ];

    let dataset = Crystal::new(lattice, positions.to_vec(), types.to_vec())
        .analyze()
        .symprec(SYMPREC)
        .dataset()
        .expect("dataset failed for La2NiO4");

    // P4₂/ncm (#138), Hall 422, 点群 4/mmm (D₄ₕ)
    assert_eq!(dataset.spacegroup_number, 138,
        "La₂NiO₄ expected P4₂/ncm (#138), got #{}", dataset.spacegroup_number);
    assert_eq!(dataset.hall_number, 422,
        "expected Hall 422, got {}", dataset.hall_number);

    let sg_type = SpaceGroupType::from_hall(dataset.hall_number)
        .expect("SpaceGroupType::from_hall failed");
    assert_eq!(sg_type.international_short.trim(), "P4_2/ncm");
    assert_eq!(sg_type.schoenflies.trim(), "D4h^16");

    // 点群: 4/mmm (D₄ₕ) = point group number 15
    let (_symbol, _transform, pg_number) =
        spg_get_pointgroup(&dataset.rotations).expect("spg_get_pointgroup failed");
    assert_eq!(pg_number, 15, "point group should be 4/mmm (D4h)");

    // D₄ₕ 点群: 16 个对称操作
    assert_eq!(dataset.n_operations, 16,
        "expected 16 symmetry ops (D4h), got {}", dataset.n_operations);

    // 原胞: 28 原子（无纯平移约化）
    assert_eq!(dataset.n_atoms, 28,
        "primitive should have 28 atoms, got {}", dataset.n_atoms);
}
