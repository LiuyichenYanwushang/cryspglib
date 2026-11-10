use cryspglib::spin_group::{
    CollinearAction, CollinearMoments, CollinearSpinGroup, SpinAction, SpinAtomMapping,
    SpinGroupError, SpinGroupKind, SpinTolerances, detect_collinear_spin_group,
    validate_collinear_spin_group,
};
use cryspglib::{Crystal, OperationKind};

const I: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const C4: [[i32; 3]; 3] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
const ORTHO: [[f64; 3]; 3] = [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]];
const CUBIC: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

fn afm() -> Crystal {
    Crystal::new(ORTHO, vec![[0.; 3], [0.5, 0., 0.]], vec![26, 26]).unwrap()
}

fn cart_length(lattice: [[f64; 3]; 3], v: [f64; 3]) -> f64 {
    lattice
        .map(|row| (0..3).map(|j| row[j] * v[j]).sum::<f64>())
        .into_iter()
        .fold(0., f64::hypot)
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}

/// Independent certificate check: geometry, bijection, residuals, full group
/// law (including the unitary/anti-unitary XOR) and mapping composition.
fn check_certificate(crystal: &Crystal, group: &CollinearSpinGroup) {
    let operations = group.operations();
    let moments = group.order().signed_moments();
    let tolerance = group.tolerances();
    assert_eq!(moments.len(), crystal.natom());
    assert_eq!(group.len(), operations.len());
    assert_eq!(group.kind(), SpinGroupKind::Collinear);
    for (g, op) in operations.iter().enumerate() {
        let mapping: &[SpinAtomMapping<f64>] =
            group.atom_mappings(g).expect("mapping per operation");
        assert_eq!(mapping.len(), crystal.natom());
        let mut targets: Vec<_> = mapping.iter().map(|a| a.target).collect();
        targets.sort_unstable();
        assert_eq!(targets, (0..crystal.natom()).collect::<Vec<_>>());
        for (i, atom) in mapping.iter().enumerate() {
            assert_eq!(crystal.types()[i], crystal.types()[atom.target]);
            let r = crystal.positions()[i];
            let target = crystal.positions()[atom.target];
            for (axis, &target_coordinate) in target.iter().enumerate() {
                let lhs = op.translation[axis]
                    + (0..3)
                        .map(|b| f64::from(op.rotation[axis][b]) * r[b])
                        .sum::<f64>();
                close(
                    lhs,
                    target_coordinate + atom.cell_shift[axis] as f64 + atom.position_residual[axis],
                );
            }
            assert!(cart_length(crystal.lattice(), atom.position_residual) <= tolerance.position);
            close(
                atom.moment_residual,
                moments[atom.target] - f64::from(op.action.sign()) * moments[i],
            );
            assert!(atom.moment_residual.abs() <= tolerance.moment);
        }
        let inverse = group.inverse_index(g).expect("inverse");
        assert_eq!(
            group.product(g, inverse).expect("product").operation_index,
            group.identity_index()
        );
        assert_eq!(
            group.product(inverse, g).expect("product").operation_index,
            group.identity_index()
        );
        for (h, right) in operations.iter().enumerate() {
            let product = group.product(g, h).expect("closed product");
            let p = &operations[product.operation_index];
            assert_eq!(op.action.compose(right.action), p.action);
            assert_eq!(op.kind.compose(right.kind), p.kind);
            for a in 0..3 {
                for b in 0..3 {
                    assert_eq!(
                        (0..3)
                            .map(|c| op.rotation[a][c] * right.rotation[c][b])
                            .sum::<i32>(),
                        p.rotation[a][b]
                    );
                }
                let raw_t = op.translation[a]
                    + (0..3)
                        .map(|b| f64::from(op.rotation[a][b]) * right.translation[b])
                        .sum::<f64>();
                close(
                    raw_t,
                    p.translation[a]
                        + product.cell_shift[a] as f64
                        + product.translation_residual[a],
                );
            }
            assert!(
                cart_length(crystal.lattice(), product.translation_residual) <= tolerance.position
            );
            assert_eq!(product.action_residual, 0.0);
            let right_mapping = group.atom_mappings(h).expect("mapping");
            let product_mapping = group
                .atom_mappings(product.operation_index)
                .expect("mapping");
            for (i, hi) in right_mapping.iter().enumerate() {
                let gi = mapping[hi.target];
                let pi = product_mapping[i];
                assert_eq!(gi.target, pi.target);
                for axis in 0..3 {
                    assert_eq!(
                        gi.cell_shift[axis]
                            + (0..3)
                                .map(|b| i64::from(op.rotation[axis][b]) * hi.cell_shift[b])
                                .sum::<i64>(),
                        pi.cell_shift[axis] + product.cell_shift[axis]
                    );
                    close(
                        gi.position_residual[axis]
                            + (0..3)
                                .map(|b| f64::from(op.rotation[axis][b]) * hi.position_residual[b])
                                .sum::<f64>(),
                        pi.position_residual[axis] + product.translation_residual[axis],
                    );
                }
                close(
                    gi.moment_residual + f64::from(op.action.sign()) * hi.moment_residual,
                    pi.moment_residual,
                );
            }
            for k in 0..operations.len() {
                let hk = group.product(h, k).expect("product");
                let lhs = group.product(product.operation_index, k).expect("product");
                let rhs = group.product(g, hk.operation_index).expect("product");
                assert_eq!(lhs.operation_index, rhs.operation_index);
                for axis in 0..3 {
                    assert_eq!(
                        product.cell_shift[axis] + lhs.cell_shift[axis],
                        rhs.cell_shift[axis]
                            + (0..3)
                                .map(|b| i64::from(op.rotation[axis][b]) * hk.cell_shift[b])
                                .sum::<i64>()
                    );
                }
            }
        }
    }
}

fn find_operation(
    group: &CollinearSpinGroup,
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
    action: CollinearAction,
    kind: OperationKind,
) -> Option<usize> {
    group.operations().iter().position(|op| {
        op.rotation == rotation
            && op
                .translation
                .iter()
                .zip(translation.iter())
                .all(|(a, b)| (a - b).abs() < 1e-12)
            && op.action == action
            && op.kind == kind
    })
}

#[test]
fn ferromagnet_full_group_has_unitary_and_anti_preserve_halves() {
    let crystal = afm();
    let vectors = [[1., 2., 3.], [1., 2., 3.]];
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Vectors(&vectors),
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(group.len(), 32);
    let unitary = group
        .operations()
        .iter()
        .filter(|op| op.kind == OperationKind::Unitary)
        .count();
    let anti = group.len() - unitary;
    assert_eq!((unitary, anti), (16, 16));
    assert!(
        group
            .operations()
            .iter()
            .all(|op| op.action == CollinearAction::Preserve)
    );
    let axis = group.order().axis().unwrap();
    for (a, expected) in axis.into_iter().zip([1., 2., 3.]) {
        close(a, expected / 14_f64.sqrt());
    }
    assert_eq!(group.order().signed_moments().len(), 2);
    check_certificate(&crystal, &group);
}

#[test]
fn antiferromagnetic_half_translation_is_unitary_swap_with_anti_partner() {
    let crystal = afm();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2.]),
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(group.len(), 32);
    let half = find_operation(
        &group,
        I,
        [0.5, 0., 0.],
        CollinearAction::Swap,
        OperationKind::Unitary,
    )
    .expect("unitary half translation");
    let mapping = group.atom_mappings(half).unwrap();
    assert_eq!(
        mapping
            .iter()
            .map(|a| (a.target, a.cell_shift))
            .collect::<Vec<_>>(),
        vec![(1, [0, 0, 0]), (0, [1, 0, 0])]
    );
    let square = group.product(half, half).unwrap();
    assert_eq!(square.operation_index, group.identity_index());
    assert_eq!(square.cell_shift, [1, 0, 0]);
    assert_eq!(square.translation_residual, [0.; 3]);

    let anti_half = find_operation(
        &group,
        I,
        [0.5, 0., 0.],
        CollinearAction::Swap,
        OperationKind::Antiunitary,
    )
    .expect("anti-unitary half translation");
    let anti_square = group.product(anti_half, anti_half).unwrap();
    let anti_square_op = &group.operations()[anti_square.operation_index];
    assert_eq!(anti_square_op.kind, OperationKind::Unitary);
    assert_eq!(anti_square_op.action, CollinearAction::Preserve);
    assert_eq!(anti_square.operation_index, group.identity_index());
    assert_eq!(anti_square.cell_shift, [1, 0, 0]);
    check_certificate(&crystal, &group);
}

#[test]
fn altermagnet_c4_rotation_swaps_without_same_k_requirement() {
    let crystal = Crystal::new(
        CUBIC,
        vec![[0.5, 0., 0.], [0., 0.5, 0.], [0.; 3]],
        vec![26, 26, 8],
    )
    .unwrap();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2., 0.]),
        SpinTolerances::default(),
    )
    .unwrap();
    let unitary_c4 = find_operation(
        &group,
        C4,
        [0.; 3],
        CollinearAction::Swap,
        OperationKind::Unitary,
    );
    assert!(unitary_c4.is_some());
    assert!(
        find_operation(
            &group,
            C4,
            [0.; 3],
            CollinearAction::Swap,
            OperationKind::Antiunitary
        )
        .is_some()
    );
    assert!(
        find_operation(
            &group,
            I,
            [0.; 3],
            CollinearAction::Swap,
            OperationKind::Unitary
        )
        .is_none()
    );
    check_certificate(&crystal, &group);
}

#[test]
fn signed_and_vector_inputs_describe_the_same_group() {
    let crystal = afm();
    let signed = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2.]),
        SpinTolerances::default(),
    )
    .unwrap();
    let vectors = [[0., 0., 2.], [0., 0., -2.]];
    let vector_group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Vectors(&vectors),
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(signed.len(), vector_group.len());
    for (a, b) in signed.operations().iter().zip(vector_group.operations()) {
        assert_eq!(a.rotation, b.rotation);
        assert_eq!(a.translation, b.translation);
        assert_eq!(a.action, b.action);
        assert_eq!(a.kind, b.kind);
    }
    assert!(signed.order().axis().is_none());
    let axis = vector_group.order().axis().unwrap();
    close(axis[2].abs(), 1.0);
    for atom in 0..2 {
        close(
            signed.order().signed_moments()[atom].abs(),
            vector_group.order().signed_moments()[atom].abs(),
        );
    }
    check_certificate(&crystal, &vector_group);
}

#[test]
fn explicit_full_set_is_validated_and_wrong_anti_partner_is_refused() {
    let crystal = afm();
    let moments = [2., -2.];
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        SpinTolerances::default(),
    )
    .unwrap();
    let validated = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        group.operations(),
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(validated.len(), group.len());

    let mut wrong = group.operations().to_vec();
    let anti = wrong
        .iter()
        .position(|op| op.kind == OperationKind::Antiunitary)
        .unwrap();
    wrong[anti].action = match wrong[anti].action {
        CollinearAction::Preserve => CollinearAction::Swap,
        CollinearAction::Swap => CollinearAction::Preserve,
    };
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &wrong,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AntiunitaryPartner { .. }));

    let unitary_only: Vec<_> = group
        .operations()
        .iter()
        .copied()
        .filter(|op| op.kind == OperationKind::Unitary)
        .collect();
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &unitary_only,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::KindCount { .. }));
}

#[test]
fn all_zero_noncollinear_and_invalid_tolerances_are_refused() {
    let crystal = afm();
    let error = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[0., 0.]),
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AllZeroMoments));

    let vectors = [[1., 0., 0.], [0., 1., 0.]];
    let error = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Vectors(&vectors),
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::NonCollinear { .. }));

    let error = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[1., -1.]),
        SpinTolerances {
            position: 0.0,
            ..SpinTolerances::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SpinGroupError::InvalidTolerance {
            name: "position",
            ..
        }
    ));
}

#[test]
fn identity_mapping_and_accessors_are_consistent() {
    let crystal = afm();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2.]),
        SpinTolerances::default(),
    )
    .unwrap();
    let identity = &group.operations()[group.identity_index()];
    assert_eq!(identity.rotation, I);
    assert_eq!(identity.translation, [0.; 3]);
    assert_eq!(identity.action, CollinearAction::Preserve);
    assert_eq!(identity.kind, OperationKind::Unitary);
    assert_eq!(group.tolerances(), SpinTolerances::default());
    let identity_mapping = group.atom_mappings(group.identity_index()).unwrap();
    for (index, mapping) in identity_mapping.iter().enumerate() {
        assert_eq!(mapping.target, index);
        assert_eq!(mapping.cell_shift, [0, 0, 0]);
        assert_eq!(mapping.position_residual, [0.; 3]);
        assert_eq!(mapping.moment_residual, 0.0);
    }
    check_certificate(&crystal, &group);
}
