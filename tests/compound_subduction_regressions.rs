//! Independent physical-character references for compound subduction.

use cryspglib::irrep::query;
use cryspglib::irrep::types::CompoundSelectedArmCharacter;
use num_complex::Complex64;

#[test]
fn c4h_compounds_match_polar_and_axial_plane_representations() {
    // In the canonical C4h cell, (x,y) has character Rxx + Ryy.
    // The axial plane has character det(R) * (Rxx + Ryy).
    // These physical models supply a reference independent of summing CIR rows.
    for (label, axial, source_ids) in [
        ("GM3+GM4+", true, [4077, 4078]),
        ("GM3-GM4-", false, [4081, 4082]),
    ] {
        let record = query::irreps_of(83)
            .iter()
            .find(|record| record.ml == label)
            .unwrap();
        let CompoundSelectedArmCharacter::DistinctComponentSum {
            first,
            second,
            block_trace,
        } = record.compound_selected_arm_view().unwrap()
        else {
            panic!("{label} must have two distinct CIR constituents");
        };
        assert_eq!([first.irnumber, second.irnumber], source_ids);
        assert_eq!([first.dimension, second.dimension], [1, 1]);
        assert_eq!(block_trace.len(), 8);
        for (value, operation) in block_trace.values().iter().zip(block_trace.operations()) {
            let r = operation.rotation;
            let determinant = r[0] * (r[4] * r[8] - r[5] * r[7])
                - r[1] * (r[3] * r[8] - r[5] * r[6])
                + r[2] * (r[3] * r[7] - r[4] * r[6]);
            let expected = (r[0] + r[4]) * if axial { determinant } else { 1 };
            assert!(
                (value - Complex64::new(f64::from(expected), 0.0)).norm() < 1e-9,
                "{label}, {operation:?}: {value} != physical trace {expected}"
            );
        }
    }
}

#[test]
fn sg23_realification_has_opposite_translation_eigenvalues() {
    let record = query::irreps_of(23)
        .iter()
        .find(|record| record.ml == "W1W1")
        .unwrap();
    let CompoundSelectedArmCharacter::ConjugateRealification { seed, block_trace } =
        record.compound_selected_arm_view().unwrap()
    else {
        panic!("SG23 W1W1 must exercise realification");
    };
    assert_eq!((record.kx, record.ky, record.kz, record.kd), (1, 1, 1, 2));
    assert_eq!(seed.irnumber, 716);
    assert_eq!(seed.dimension, 1);
    assert_eq!(block_trace.dimension(), 2);
    // k.L = 3/4 for the I-centring translation. The seed eigenvalue is -i,
    // its conjugate is +i, and the physical trace is zero rather than -2i.
    let index = seed
        .row
        .operations()
        .iter()
        .position(|operation| {
            operation.rotation == [1, 0, 0, 0, 1, 0, 0, 0, 1] && operation.translation == [0.5; 3]
        })
        .unwrap();
    let value = seed.row.get(index).unwrap();
    assert!((value - Complex64::new(0.0, -1.0)).norm() < 1e-9);
    assert!((value.conj() - Complex64::new(0.0, 1.0)).norm() < 1e-9);
    assert!(block_trace.get(index).unwrap().norm() < 1e-9);
}
