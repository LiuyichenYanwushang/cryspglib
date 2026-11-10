//! End-to-end tests for standard-cell exact collinear spin-space-group
//! identification.

use cryspglib::Crystal;
use cryspglib::spin_group::database::collinear_spin_space_groups;
use cryspglib::spin_group::{
    CollinearMoments, SpinSpaceGroupIdentificationError, SpinTolerances,
    detect_collinear_spin_group, identify_collinear_signature, identify_collinear_spin_space_group,
};

const CUBIC: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const ORTHO: [[f64; 3]; 3] = [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]];

fn cubic_ferromagnet(moment: [f64; 3]) -> Crystal {
    Crystal::new(CUBIC, vec![[0.; 3]], vec![26])
        .unwrap()
        .with_magnetic(vec![moment])
        .unwrap()
}

fn orthorhombic_afm() -> Crystal {
    Crystal::new(ORTHO, vec![[0.; 3], [0.5, 0., 0.]], vec![26, 26])
        .unwrap()
        .with_magnetic(vec![[0., 0., 2.], [0., 0., -2.]])
        .unwrap()
}

fn detect(
    crystal: &Crystal,
    moments: CollinearMoments<'_>,
) -> cryspglib::spin_group::CollinearSpinGroup {
    detect_collinear_spin_group(crystal, moments, SpinTolerances::default()).unwrap()
}

#[test]
fn cubic_ferromagnet_identifies_the_pm3m_type_i_record() {
    let crystal = cubic_ferromagnet([0., 0., 1.]);
    let group = detect(&crystal, CollinearMoments::Vectors(&[[0., 0., 1.]]));
    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "221.221.1.1.L");
    assert_eq!(identification.number.g0(), 221);
    assert_eq!(identification.record.group_type, "typeI");
    assert_eq!(identification.record.spin_point_group, "∞m");
    assert_eq!(identification.record.basis, "G0");
    assert_eq!(identification.international_symbol(), "P1m1-31m∞m1");
    assert_eq!(identification.record.identifier, "221.221.1.1.L");
}

#[test]
fn signed_scalar_ferromagnet_identifies_the_same_record() {
    let crystal = cubic_ferromagnet([0., 0., 1.]);
    let group = detect(&crystal, CollinearMoments::Signed(&[2.0]));
    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "221.221.1.1.L");
}

#[test]
fn global_spin_rotation_does_not_change_the_number() {
    let crystal = cubic_ferromagnet([1., 1., 1.]);
    let group = detect(&crystal, CollinearMoments::Vectors(&[[1., 1., 1.]]));
    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "221.221.1.1.L");
}

#[test]
fn shifted_origin_is_explicitly_refused() {
    let crystal = Crystal::new(CUBIC, vec![[0.1, 0.1, 0.1]], vec![26])
        .unwrap()
        .with_magnetic(vec![[0., 0., 1.]])
        .unwrap();
    let group = detect(&crystal, CollinearMoments::Vectors(&[[0., 0., 1.]]));
    let error = identify_collinear_spin_space_group(&group).unwrap_err();
    assert!(matches!(
        error,
        SpinSpaceGroupIdentificationError::NoExactMatch
    ));
}

#[test]
fn orthorhombic_afm_identifies_a_pinned_record() {
    let crystal = orthorhombic_afm();
    let group = detect(&crystal, CollinearMoments::Signed(&[2., -2.]));
    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "47.47.2.1.L");
    assert_eq!(identification.number.g0(), 47);
    assert_eq!(identification.number.l0(), 47);
    assert_eq!(identification.record.operations.len(), 16);
    assert_eq!(identification.record.basis, "L0");
    assert_eq!(identification.international_symbol(), "P_a1m1m1m∞m1");
}

#[test]
fn every_pinned_quotient_signature_resolves_to_itself() {
    let mut count = 0usize;
    for record in collinear_spin_space_groups() {
        let identification = identify_collinear_signature(record.operations)
            .unwrap_or_else(|error| panic!("{}: {error}", record.identifier));
        assert_eq!(identification.number, record.number);
        assert_eq!(identification.record.identifier, record.identifier);
        count += 1;
    }
    assert_eq!(count, 1421);
}
