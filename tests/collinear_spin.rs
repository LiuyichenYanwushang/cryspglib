use cryspglib::{
    CollinearMoments, CollinearSpinError, CollinearSpinOperation, CollinearSpinTolerances, Crystal,
    SpinChannelAction, ValidatedCollinearSpinGroup, detect_collinear_spin_group,
    validate_collinear_spin_group,
};
use num_complex::Complex64 as C;
use std::f64::consts::TAU;

const I: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const C4: [[i32; 3]; 3] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
const ORTHO: [[f64; 3]; 3] = [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]];
const CUBIC: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

fn operation(
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
    channel_action: SpinChannelAction,
) -> CollinearSpinOperation {
    CollinearSpinOperation {
        rotation,
        translation,
        channel_action,
    }
}
fn half_swap() -> CollinearSpinOperation {
    operation(I, [0.5, 0., 0.], SpinChannelAction::Swap)
}
fn afm() -> Crystal {
    Crystal::new(ORTHO, vec![[0.; 3], [0.5, 0., 0.]], vec![26, 26]).unwrap()
}
fn validate(
    crystal: &Crystal,
    moments: &[f64],
    operations: &[CollinearSpinOperation],
) -> Result<ValidatedCollinearSpinGroup, CollinearSpinError> {
    validate_collinear_spin_group(
        crystal,
        CollinearMoments::Signed(moments),
        operations,
        CollinearSpinTolerances::default(),
    )
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}
fn cart_length(lattice: [[f64; 3]; 3], v: [f64; 3]) -> f64 {
    lattice
        .map(|row| (0..3).map(|j| row[j] * v[j]).sum::<f64>())
        .into_iter()
        .fold(0., f64::hypot)
}

/// Independent original-coordinate oracle; no production math/mapping helpers.
fn check_certificate(crystal: &Crystal, group: &ValidatedCollinearSpinGroup) {
    let operations = group.operations();
    let table = group.multiplication_table();
    let moments = group.signed_moments();
    let tolerance = group.tolerances();
    assert_eq!(moments.len(), crystal.natom());
    assert_eq!(table.len(), operations.len());
    for (g, mapped) in operations.iter().enumerate() {
        let op = mapped.operation;
        assert_eq!(mapped.atom_mapping.len(), crystal.natom());
        let mut targets: Vec<_> = mapped.atom_mapping.iter().map(|a| a.target).collect();
        targets.sort_unstable();
        assert_eq!(targets, (0..crystal.natom()).collect::<Vec<_>>());
        for (i, atom) in mapped.atom_mapping.iter().enumerate() {
            assert_eq!(crystal.types()[i], crystal.types()[atom.target]);
            let r = crystal.positions()[i];
            let target = crystal.positions()[atom.target];
            for (a, &target_coordinate) in target.iter().enumerate() {
                let lhs = op.translation[a]
                    + (0..3)
                        .map(|b| f64::from(op.rotation[a][b]) * r[b])
                        .sum::<f64>();
                close(
                    lhs,
                    target_coordinate + atom.cell_shift[a] as f64 + atom.position_residual[a],
                );
            }
            assert!(cart_length(crystal.lattice(), atom.position_residual) <= tolerance.position);
            close(
                atom.moment_residual,
                moments[atom.target] - f64::from(op.channel_action.sign()) * moments[i],
            );
            assert!(atom.moment_residual.abs() <= tolerance.moment);
        }
        let inverse = group.inverse_indices()[g];
        assert_eq!(table[g][inverse].operation_index, group.identity_index());
        assert_eq!(table[inverse][g].operation_index, group.identity_index());
        for (h, right) in operations.iter().enumerate() {
            let product = table[g][h];
            let p = &operations[product.operation_index];
            assert_eq!(
                op.channel_action.compose(right.operation.channel_action),
                p.operation.channel_action
            );
            for a in 0..3 {
                for b in 0..3 {
                    assert_eq!(
                        (0..3)
                            .map(|c| op.rotation[a][c] * right.operation.rotation[c][b])
                            .sum::<i32>(),
                        p.operation.rotation[a][b]
                    );
                }
                let raw_t = op.translation[a]
                    + (0..3)
                        .map(|b| f64::from(op.rotation[a][b]) * right.operation.translation[b])
                        .sum::<f64>();
                close(
                    raw_t,
                    p.operation.translation[a]
                        + product.cell_shift[a] as f64
                        + product.translation_residual[a],
                );
            }
            assert!(
                cart_length(crystal.lattice(), product.translation_residual) <= tolerance.position
            );
            for (i, hi) in right.atom_mapping.iter().enumerate() {
                let gi = mapped.atom_mapping[hi.target];
                let pi = p.atom_mapping[i];
                assert_eq!(gi.target, pi.target);
                for a in 0..3 {
                    assert_eq!(
                        gi.cell_shift[a]
                            + (0..3)
                                .map(|b| i64::from(op.rotation[a][b]) * hi.cell_shift[b])
                                .sum::<i64>(),
                        pi.cell_shift[a] + product.cell_shift[a]
                    );
                    close(
                        gi.position_residual[a]
                            + (0..3)
                                .map(|b| f64::from(op.rotation[a][b]) * hi.position_residual[b])
                                .sum::<f64>(),
                        pi.position_residual[a] + product.translation_residual[a],
                    );
                }
                close(
                    gi.moment_residual + f64::from(op.channel_action.sign()) * hi.moment_residual,
                    pi.moment_residual,
                );
            }
            for k in 0..operations.len() {
                let hk = table[h][k];
                let lhs = table[product.operation_index][k];
                let rhs = table[g][hk.operation_index];
                assert_eq!(lhs.operation_index, rhs.operation_index);
                for a in 0..3 {
                    assert_eq!(
                        product.cell_shift[a] + lhs.cell_shift[a],
                        rhs.cell_shift[a]
                            + (0..3)
                                .map(|b| i64::from(op.rotation[a][b]) * hk.cell_shift[b])
                                .sum::<i64>()
                    );
                }
            }
        }
    }
}

#[test]
fn ferromagnet_preserves_channels_even_when_space_rotation_moves_the_spin_axis() {
    let crystal = afm();
    let vectors = [[1., 2., 3.], [1., 2., 3.]];
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Vectors(&vectors),
        CollinearSpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(group.operations().len(), 16);
    assert!(
        group
            .operations()
            .iter()
            .all(|g| g.operation.channel_action == SpinChannelAction::Preserve)
    );
    assert!(
        group
            .operations()
            .iter()
            .any(|g| g.operation.rotation == [[1, 0, 0], [0, -1, 0], [0, 0, -1]])
    );
    let axis = group.spin_axis().unwrap();
    for (a, expected) in axis.into_iter().zip([1., 2., 3.]) {
        close(a, expected / 14_f64.sqrt());
    }
    check_certificate(&crystal, &group);
}

#[test]
fn antiferromagnetic_half_translation_swaps_unitarily_and_squares_to_one_cell() {
    let crystal = afm();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2.]),
        CollinearSpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(group.operations().len(), 16);
    let g = group
        .operations()
        .iter()
        .position(|g| g.operation == half_swap())
        .unwrap();
    assert_eq!(
        group.operations()[g]
            .atom_mapping
            .iter()
            .map(|a| (a.target, a.cell_shift))
            .collect::<Vec<_>>(),
        vec![(1, [0, 0, 0]), (0, [1, 0, 0])]
    );
    let square = group.multiplication_table()[g][g];
    assert_eq!(square.operation_index, group.identity_index());
    assert_eq!(square.cell_shift, [1, 0, 0]);
    assert_eq!(square.translation_residual, [0.; 3]);
    check_certificate(&crystal, &group);
}

#[test]
fn altermagnet_rotation_swaps_but_does_not_force_same_k_degeneracy() {
    let crystal = Crystal::new(
        CUBIC,
        vec![[0.5, 0., 0.], [0., 0.5, 0.], [0.; 3]],
        vec![26, 26, 8],
    )
    .unwrap();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2., 0.]),
        CollinearSpinTolerances::default(),
    )
    .unwrap();
    assert!(
        group
            .operations()
            .iter()
            .any(|g| g.operation == operation(C4, [0.; 3], SpinChannelAction::Swap))
    );
    assert!(group.operations().iter().any(|g| g.operation.rotation
        == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        && g.operation.channel_action == SpinChannelAction::Preserve));
    // The nonmagnetic marker excludes the otherwise allowed (1/2,1/2,0) Swap translation.
    assert!(!group.operations().iter().any(
        |g| g.operation.rotation == I && g.operation.channel_action == SpinChannelAction::Swap
    ));
    let up = |k: [f64; 2]| 1.3 * (TAU * k[0]).cos() + 0.7 * (TAU * k[1]).cos();
    let down = |k: [f64; 2]| 0.7 * (TAU * k[0]).cos() + 1.3 * (TAU * k[1]).cos();
    let k = [0.17, 0.31];
    let q = [-k[1], k[0]]; // W^{-T}k = Wk for this C4, WITHOUT the antiunitary minus sign.
    close(up(k), down(q));
    assert!((up(k) - down(k)).abs() > 0.1);
    check_certificate(&crystal, &group);
}

type M2 = [[C; 2]; 2];
fn multiply(a: M2, b: M2) -> M2 {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..2).map(|k| a[i][k] * b[k][j]).sum()))
}
fn adjoint(a: M2) -> M2 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i].conj()))
}
fn hopping_blocks(k: f64) -> [M2; 2] {
    let theta = TAU * k;
    let x = C::from_polar(1., -theta);
    let f = C::i() * (C::new(1., 0.) - x);
    let odd = 2. * theta.sin();
    [
        [[C::new(odd + 1., 0.), f], [f.conj(), C::new(odd - 1., 0.)]],
        [[C::new(odd - 1., 0.), f], [f.conj(), C::new(odd + 1., 0.)]],
    ]
}

#[test]
fn complex_hopping_consumer_example_distinguishes_unitary_from_antiunitary_swap() {
    let group = validate(
        &afm(),
        &[2., -2.],
        &[CollinearSpinOperation::identity(), half_swap()],
    )
    .unwrap();
    let swap = group.operations()[1].operation;
    assert_eq!(swap.channel_action, SpinChannelAction::Swap);
    for k in [0.13, 0.25, -0.31] {
        let blocks = hopping_blocks(k);
        let p = [
            [C::new(0., 0.), C::from_polar(1., -TAU * k)],
            [C::new(1., 0.), C::new(0., 0.)],
        ];
        // Example orbital sewing supplied by a consumer, NOT inferred by cryspglib.
        let sewn = multiply(multiply(p, blocks[0]), adjoint(p));
        let wrong = multiply(
            multiply(p, blocks[0].map(|row| row.map(|v| v.conj()))),
            adjoint(p),
        );
        for i in 0..2 {
            for j in 0..2 {
                assert!((blocks[1][i][j] - sewn[i][j]).norm() < 1e-10);
                assert!((blocks[0][i][j] - blocks[0][j][i].conj()).norm() < 1e-10);
            }
        }
        assert!((wrong[0][1] - blocks[1][0][1]).norm() > 0.1);
        let opposite_k = hopping_blocks(-k);
        let trace = |a: M2| a[0][0] + a[1][1];
        assert!((trace(blocks[0]).conj() - trace(opposite_k[1])).norm() > 0.1);
        // Equivalently average H0 and eta*Hz; swapping never conjugates f.
        for (up_row, down_row) in blocks[0].iter().zip(blocks[1]) {
            for (&up, down) in up_row.iter().zip(down_row) {
                let h0 = (up + down) / 2.;
                let hz = (up - down) / 2.;
                assert!((h0 + f64::from(swap.channel_action.sign()) * hz - down).norm() < 1e-10);
            }
        }
    }
}

#[test]
fn screw_and_changed_site_and_seitz_representatives_preserve_all_integer_cocycles() {
    let rotation = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];
    let base_ops = [
        CollinearSpinOperation::identity(),
        operation(rotation, [0., 0., 0.5], SpinChannelAction::Swap),
    ];
    let positions = [[0.17, 0.23, 0.11], [-0.17, -0.23, 0.61]];
    let base = Crystal::new(ORTHO, positions.to_vec(), vec![26, 26]).unwrap();
    let group = validate(&base, &[3., -3.], &base_ops).unwrap();
    assert_eq!(group.multiplication_table()[1][1].cell_shift, [0, 0, 1]);
    check_certificate(&base, &group);
    let site_shifts = [[2, -3, 1], [-1, 4, -2]];
    let rep_shifts = [[2, -1, 3], [1, -2, 4]];
    let changed = Crystal::new(
        ORTHO,
        positions
            .iter()
            .enumerate()
            .map(|(i, r)| std::array::from_fn(|j| r[j] + f64::from(site_shifts[i][j])))
            .collect(),
        vec![26, 26],
    )
    .unwrap();
    let changed_ops: Vec<_> = base_ops
        .iter()
        .enumerate()
        .map(|(i, g)| CollinearSpinOperation {
            translation: std::array::from_fn(|j| g.translation[j] + f64::from(rep_shifts[i][j])),
            ..*g
        })
        .collect();
    let changed_group = validate(&changed, &[3., -3.], &changed_ops).unwrap();
    assert_eq!(
        changed_group
            .operations()
            .iter()
            .map(|g| g.operation)
            .collect::<Vec<_>>(),
        changed_ops
    );
    for g in 0..2 {
        for i in 0..2 {
            let atom = group.operations()[g].atom_mapping[i];
            for a in 0..3 {
                assert_eq!(
                    changed_group.operations()[g].atom_mapping[i].cell_shift[a],
                    atom.cell_shift[a] + i64::from(rep_shifts[g][a])
                        - i64::from(site_shifts[atom.target][a])
                        + (0..3)
                            .map(|b| i64::from(base_ops[g].rotation[a][b] * site_shifts[i][b]))
                            .sum::<i64>()
                );
            }
        }
    }
    for g in 0..2 {
        for h in 0..2 {
            let original = group.multiplication_table()[g][h];
            let changed = changed_group.multiplication_table()[g][h];
            assert_eq!(original.operation_index, changed.operation_index);
            for a in 0..3 {
                assert_eq!(
                    changed.cell_shift[a],
                    original.cell_shift[a] + i64::from(rep_shifts[g][a])
                        - i64::from(rep_shifts[original.operation_index][a])
                        + (0..3)
                            .map(|b| i64::from(base_ops[g].rotation[a][b] * rep_shifts[h][b]))
                            .sum::<i64>()
                );
            }
        }
    }
    check_certificate(&changed, &changed_group);
}

#[test]
fn approximate_position_moment_and_seitz_residuals_are_not_rounded_into_shifts() {
    let crystal = Crystal::new(
        ORTHO,
        vec![[0.123456789, 0., 0.], [0.623456989, 0., 0.]],
        vec![26, 26],
    )
    .unwrap();
    let swap = operation(I, [0.5000001, 0., 0.], SpinChannelAction::Swap);
    let group = validate(
        &crystal,
        &[2., -1.999999],
        &[CollinearSpinOperation::identity(), swap],
    )
    .unwrap();
    let map = &group.operations()[1].atom_mapping;
    assert_eq!(map[0].cell_shift, [0, 0, 0]);
    assert_eq!(map[1].cell_shift, [1, 0, 0]);
    close(map[0].position_residual[0], -1e-7);
    close(map[1].position_residual[0], 3e-7);
    close(map[0].moment_residual, 1e-6);
    let square = group.multiplication_table()[1][1];
    assert_eq!(square.cell_shift, [1, 0, 0]);
    close(square.translation_residual[0], 2e-7);
    check_certificate(&crystal, &group);
}

#[test]
fn skew_cell_mapping_uses_cartesian_images_beyond_the_nearest_27() {
    let lattice = [[1., 10., 0.], [0., 0.1, 0.], [0., 0., 1.]];
    let crystal = Crystal::new(lattice, vec![[0.; 3]], vec![26]).unwrap();
    // A deliberately approximate identity: the unique image within tolerance
    // is (5,0,0), not componentwise round(d) nor any of its 27 neighbours.
    let op = operation(I, [0.1, 0.49, 0.], SpinChannelAction::Preserve);
    let tolerance = CollinearSpinTolerances {
        position: 0.0495,
        ..Default::default()
    };
    let group =
        validate_collinear_spin_group(&crystal, CollinearMoments::Signed(&[2.]), &[op], tolerance)
            .unwrap();
    assert_eq!(group.operations()[0].atom_mapping[0].cell_shift, [5, 0, 0]);
    let delta = group.operations()[0].atom_mapping[0].position_residual;
    close(delta[0], -4.9);
    close(cart_length(lattice, delta), 0.049);
    check_certificate(&crystal, &group);
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2.]),
        &[op],
        CollinearSpinTolerances {
            position: 0.0511,
            ..tolerance
        },
    )
    .unwrap_err();
    assert!(matches!(error, CollinearSpinError::AmbiguousPeriodicImage));
}

#[test]
fn magnetic_magnitudes_are_compared_and_wrong_explicit_swap_is_refused() {
    let crystal = afm();
    let error = validate(
        &crystal,
        &[2., -3.],
        &[CollinearSpinOperation::identity(), half_swap()],
    )
    .unwrap_err();
    assert!(matches!(
        error,
        CollinearSpinError::MomentMismatch {
            operation_index: 1,
            atom_index: 0,
            ..
        }
    ));
    // Unequal magnitudes are valid magnetic input, but they remove the Swap.
    let detected = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -3.]),
        Default::default(),
    )
    .unwrap();
    assert!(
        detected
            .operations()
            .iter()
            .all(|g| g.operation.channel_action == SpinChannelAction::Preserve)
    );
    assert!(matches!(
        validate(
            &crystal,
            &[2., -2.],
            &[
                CollinearSpinOperation::identity(),
                operation(I, [0.5, 0., 0.], SpinChannelAction::Preserve)
            ]
        ),
        Err(CollinearSpinError::MomentMismatch { .. })
    ));
}

#[test]
fn zero_moment_atoms_still_have_to_map_structurally() {
    let crystal = Crystal::new(
        ORTHO,
        vec![[0.; 3], [0.5, 0., 0.], [0.25, 0.17, 0.19]],
        vec![26, 26, 8],
    )
    .unwrap();
    assert!(matches!(
        validate(
            &crystal,
            &[2., -2., 0.],
            &[CollinearSpinOperation::identity(), half_swap()]
        ),
        Err(CollinearSpinError::AtomMappingMissing {
            operation_index: 1,
            atom_index: 2
        })
    ));
}

#[test]
fn vectors_are_checked_before_projection_and_all_zero_has_no_channel_labelling() {
    let crystal = afm();
    let vectors = [[0., 0., 2.], [0.005, 0., -1.]];
    assert!(matches!(
        validate_collinear_spin_group(
            &crystal,
            CollinearMoments::Vectors(&vectors),
            &[CollinearSpinOperation::identity()],
            CollinearSpinTolerances {
                collinearity: 1e-4,
                ..Default::default()
            }
        ),
        Err(CollinearSpinError::NonCollinear { atom_index: 1, .. })
    ));
    let vectors = [[0., 0., 2.], [1e-6, 0., -1.]];
    let group = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Vectors(&vectors),
        &[CollinearSpinOperation::identity()],
        Default::default(),
    )
    .unwrap();
    assert_eq!(group.spin_axis(), Some([0., 0., 1.]));
    close(group.collinearity_residuals()[1], 1e-6);
    assert_eq!(group.signed_moments(), &[2., -1_f64.hypot(1e-6)]);
    // A loose collinearity bound cannot hide unequal actual vector magnitudes
    // behind equal axis projections. The largest-moment orbit fixes the z axis.
    let four_sites = Crystal::new(
        ORTHO,
        vec![[0.; 3], [0.5, 0., 0.], [0.1, 0.3, 0.], [0.6, 0.3, 0.]],
        vec![26, 26, 8, 8],
    )
    .unwrap();
    let vectors = [[0., 0., 1.], [0.01, 0., -1.], [0., 0., 3.], [0., 0., -3.]];
    assert!(matches!(
        validate_collinear_spin_group(
            &four_sites,
            CollinearMoments::Vectors(&vectors),
            &[CollinearSpinOperation::identity(), half_swap()],
            CollinearSpinTolerances {
                collinearity: 0.02,
                moment: 1e-6,
                ..Default::default()
            }
        ),
        Err(CollinearSpinError::MomentMismatch {
            operation_index: 1,
            atom_index: 0,
            ..
        })
    ));
    for moments in [
        CollinearMoments::Signed(&[0., 0.]),
        CollinearMoments::Signed(&[1e-6, -1e-6]),
        CollinearMoments::Vectors(&[[0.; 3]; 2]),
    ] {
        assert!(matches!(
            detect_collinear_spin_group(&crystal, moments, Default::default()),
            Err(CollinearSpinError::AllZeroMoments)
        ));
    }
}

#[test]
fn empty_duplicate_missing_identity_and_nonclosed_sets_report_witnesses() {
    let crystal = Crystal::new(CUBIC, vec![[0.; 3]], vec![26]).unwrap();
    let identity = CollinearSpinOperation::identity();
    assert!(matches!(
        validate(&crystal, &[2.], &[]),
        Err(CollinearSpinError::Empty)
    ));
    assert!(matches!(
        validate(
            &crystal,
            &[2.],
            &[
                identity,
                operation(I, [1., -2., 3.], SpinChannelAction::Preserve)
            ]
        ),
        Err(CollinearSpinError::Duplicate {
            first_index: 0,
            duplicate_index: 1
        })
    ));
    assert!(matches!(
        validate(&afm(), &[2., -2.], &[half_swap()]),
        Err(CollinearSpinError::MissingIdentity)
    ));
    let inverse_c4 = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
    let error = validate(
        &crystal,
        &[2.],
        &[
            identity,
            operation(C4, [0.; 3], SpinChannelAction::Preserve),
            operation(inverse_c4, [0.; 3], SpinChannelAction::Preserve),
        ],
    )
    .unwrap_err();
    match error {
        CollinearSpinError::NotClosed {
            left_index,
            right_index,
            product,
        } => {
            assert_eq!((left_index, right_index), (1, 1));
            assert_eq!(product.rotation, [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]);
            assert_eq!(product.channel_action, SpinChannelAction::Preserve);
        }
        other => panic!("expected concrete missing product, got {other:?}"),
    }
}

#[test]
fn detection_does_not_prune_a_nonclosed_tolerance_selected_set() {
    let crystal = Crystal::new(
        ORTHO,
        vec![[0.; 3], [0.25, 0., 0.], [0.5, 0., 0.], [0.75, 0., 0.]],
        vec![26; 4],
    )
    .unwrap();
    let error = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[1., 1.1, 1.2, 1.1]),
        CollinearSpinTolerances {
            moment: 0.15,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        matches!(error, CollinearSpinError::NotClosed { .. }),
        "{error:?}"
    );
}

#[test]
fn geometry_ambiguity_is_not_resolved_by_choosing_a_magnetic_match() {
    let crystal = Crystal::new(CUBIC, vec![[0.; 3], [1e-7, 0., 0.]], vec![26, 26]).unwrap();
    assert!(matches!(
        detect_collinear_spin_group(
            &crystal,
            CollinearMoments::Signed(&[2., -2.]),
            Default::default()
        ),
        Err(CollinearSpinError::AtomMappingAmbiguous { .. })
    ));
}

#[test]
fn invalid_dimensions_tolerances_lattice_positions_and_moments_are_refused() {
    let crystal = afm();
    let identity = [CollinearSpinOperation::identity()];
    for tolerance in [
        CollinearSpinTolerances {
            position: 0.,
            ..Default::default()
        },
        CollinearSpinTolerances {
            moment: -1.,
            ..Default::default()
        },
        CollinearSpinTolerances {
            collinearity: f64::NAN,
            ..Default::default()
        },
        CollinearSpinTolerances {
            position: f64::INFINITY,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            validate_collinear_spin_group(
                &crystal,
                CollinearMoments::Signed(&[2., -2.]),
                &identity,
                tolerance
            ),
            Err(CollinearSpinError::InvalidTolerance { .. })
        ));
    }
    for moments in [
        CollinearMoments::Signed(&[2.]),
        CollinearMoments::Signed(&[2., f64::INFINITY]),
        CollinearMoments::Vectors(&[[0., 0., f64::NAN], [0.; 3]]),
    ] {
        assert!(matches!(
            validate_collinear_spin_group(&crystal, moments, &identity, Default::default()),
            Err(CollinearSpinError::InvalidInput { .. })
        ));
    }
    for bad in [
        Crystal::new([[0.; 3]; 3], vec![[0.; 3]], vec![26]).unwrap(),
        Crystal::new(CUBIC, vec![[f64::NAN, 0., 0.]], vec![26]).unwrap(),
        Crystal::new(CUBIC, vec![[1e100, 0., 0.]], vec![26]).unwrap(),
        Crystal::new(CUBIC, vec![], vec![]).unwrap(),
        Crystal::new(CUBIC, vec![[0.; 3]], vec![26])
            .unwrap()
            .with_layer(cryspglib::cell::AperiodicAxis::Z),
    ] {
        assert!(matches!(
            validate(&bad, &[2.], &identity),
            Err(CollinearSpinError::InvalidInput { .. })
        ));
    }
}

#[test]
fn nonisometries_nonfinite_operations_and_unresolvable_image_search_are_refused() {
    let crystal = afm();
    for rotation in [
        [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
        C4,
        [[1, i32::MAX, 0], [0, 1, 0], [0, 0, 1]],
    ] {
        assert!(matches!(
            validate(
                &crystal,
                &[2., -2.],
                &[operation(rotation, [0.; 3], SpinChannelAction::Preserve)]
            ),
            Err(CollinearSpinError::InvalidRotation { operation_index: 0 })
        ));
    }
    for bad in [f64::NAN, f64::INFINITY, 1e100] {
        assert!(matches!(
            validate(
                &crystal,
                &[2., -2.],
                &[operation(I, [bad, 0., 0.], SpinChannelAction::Preserve)]
            ),
            Err(CollinearSpinError::InvalidTranslation { operation_index: 0 })
        ));
    }
    assert!(matches!(
        validate_collinear_spin_group(
            &crystal,
            CollinearMoments::Signed(&[2., -2.]),
            &[CollinearSpinOperation::identity()],
            CollinearSpinTolerances {
                position: 1e5,
                ..Default::default()
            }
        ),
        Err(CollinearSpinError::ImageSearchLimit)
    ));
}

#[test]
fn finite_lattice_isometry_can_have_wide_integer_intermediates() {
    let n = 46_340;
    let crystal = Crystal::new(
        [[1., -f64::from(n), 0.], [0., 1., 0.], [0., 0., 1.]],
        vec![[0.; 3]],
        vec![26],
    )
    .unwrap();
    let reflection = operation(
        [[n, 1 - n * n, 0], [1, -n, 0], [0, 0, 1]],
        [0.; 3],
        SpinChannelAction::Preserve,
    );
    let group = validate(
        &crystal,
        &[2.],
        &[CollinearSpinOperation::identity(), reflection],
    )
    .unwrap();
    assert_eq!(
        group.multiplication_table()[1][1].operation_index,
        group.identity_index()
    );
    // Use wide arithmetic in the independent oracle as well.
    for (i, row) in I.iter().enumerate() {
        for (j, &expected) in row.iter().enumerate() {
            let square: i128 = (0..3)
                .map(|k| {
                    i128::from(reflection.rotation[i][k]) * i128::from(reflection.rotation[k][j])
                })
                .sum();
            assert_eq!(square, i128::from(expected));
        }
    }
    assert!(
        group
            .operations()
            .iter()
            .all(|g| g.atom_mapping[0].cell_shift == [0; 3])
    );
}

#[test]
fn unequal_magnitude_orbits_and_arbitrary_origin_are_not_forced_onto_a_database_grid() {
    let origin = [0.137, 0.219, 0.317];
    let crystal = Crystal::new(
        ORTHO,
        [0., 0.25, 0.5, 0.75]
            .iter()
            .map(|x| [x + origin[0], origin[1], origin[2]])
            .collect(),
        vec![26; 4],
    )
    .unwrap();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[1., 2., -1., -2.]),
        Default::default(),
    )
    .unwrap();
    assert!(
        group
            .operations()
            .iter()
            .any(|g| g.operation == half_swap())
    );
    assert!(group.operations().iter().any(|g| {
        g.operation
            .translation
            .iter()
            .any(|t| (12. * t - (12. * t).round()).abs() > 1e-4)
    }));
    check_certificate(&crystal, &group);
}

#[test]
fn large_site_representatives_do_not_swallow_small_geometric_residuals() {
    let n = (1_i64 << 40) as f64;
    let reflection = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for x in [0., n] {
        let crystal = Crystal::new(CUBIC, vec![[x, 0., 0.]], vec![26]).unwrap();
        assert!(matches!(
            validate(
                &crystal,
                &[2.],
                &[
                    CollinearSpinOperation::identity(),
                    operation(reflection, [1e-4, 0., 0.], SpinChannelAction::Preserve)
                ]
            ),
            Err(CollinearSpinError::AtomMappingMissing {
                operation_index: 1,
                atom_index: 0
            })
        ));
        let group = validate(
            &crystal,
            &[2.],
            &[
                CollinearSpinOperation::identity(),
                operation(reflection, [1e-6, 0., 0.], SpinChannelAction::Preserve),
            ],
        )
        .unwrap();
        close(
            group.operations()[1].atom_mapping[0].position_residual[0],
            1e-6,
        );
        assert_eq!(
            group.operations()[1].atom_mapping[0].cell_shift[0],
            if x == 0. { 0 } else { -(1_i64 << 41) }
        );
    }
    let crystal = Crystal::new(CUBIC, vec![[n + 0.25, 0., 0.]], vec![26]).unwrap();
    let group = validate(
        &crystal,
        &[2.],
        &[
            operation(I, [-n, 0., 0.], SpinChannelAction::Preserve),
            operation(reflection, [n + 0.5, 0., 0.], SpinChannelAction::Preserve),
        ],
    )
    .unwrap();
    assert_eq!(
        group.operations()[1].atom_mapping[0].cell_shift[0],
        -(1_i64 << 40)
    );
    assert_eq!(
        group.operations()[1].atom_mapping[0].position_residual,
        [0.; 3]
    );
    assert_eq!(
        group.multiplication_table()[1][1].cell_shift[0],
        1_i64 << 40
    );
}

#[test]
fn singular_lattice_is_rejected_even_when_naive_float_determinant_is_nonzero() {
    let lattice = [[1., 2., 3.], [1., 2., 3.], [2., 3., 5.]].map(|row| row.map(|v| 1e11 * v));
    let crystal = Crystal::new(lattice, vec![[0.; 3]], vec![26]).unwrap();
    assert!(matches!(
        validate(&crystal, &[2.], &[CollinearSpinOperation::identity()]),
        Err(CollinearSpinError::InvalidInput { .. })
    ));
}

#[test]
fn unresolvable_subnormal_vector_axis_is_not_published_as_a_unit_vector() {
    let tiny = f64::from_bits(1);
    let crystal = Crystal::new(CUBIC, vec![[0.; 3]], vec![26]).unwrap();
    for collinearity in [0., 1e-5] {
        assert!(matches!(
            validate_collinear_spin_group(
                &crystal,
                CollinearMoments::Vectors(&[[tiny, tiny, 0.]]),
                &[CollinearSpinOperation::identity()],
                CollinearSpinTolerances {
                    moment: 0.,
                    collinearity,
                    ..Default::default()
                }
            ),
            Err(CollinearSpinError::NumericalRange {
                context: "normalizing subnormal vector moments"
            })
        ));
    }
}

#[test]
fn unresolved_fractional_metric_and_cartesian_roundoff_never_produce_certificates() {
    let n = 46_340;
    let scaled = 2_f64.powi(28);
    let large_w = Crystal::new(
        [
            [scaled, -scaled * f64::from(n), 0.],
            [0., scaled, 0.],
            [0., 0., scaled],
        ],
        vec![[0.06980541020493547, -0.34772512881875217, 0.]],
        vec![26],
    )
    .unwrap();
    let reflection = operation(
        [[n, 1 - n * n, 0], [1, -n, 0], [0, 0, 1]],
        [0.; 3],
        SpinChannelAction::Preserve,
    );
    assert!(matches!(
        validate(
            &large_w,
            &[2.],
            &[CollinearSpinOperation::identity(), reflection]
        ),
        Err(CollinearSpinError::NumericalRange { .. })
    ));
    let scale = 2_f64.powi(40);
    let lost_metric = Crystal::new(
        [
            [scale, scale, 0.],
            [0., 2_f64.powi(13), 0.],
            [0., 0., scale],
        ],
        vec![[0.; 3]],
        vec![26],
    )
    .unwrap();
    let swap_xy = operation(
        [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        [0.; 3],
        SpinChannelAction::Preserve,
    );
    assert!(matches!(
        validate(
            &lost_metric,
            &[2.],
            &[CollinearSpinOperation::identity(), swap_xy]
        ),
        Err(CollinearSpinError::NumericalRange { .. })
    ));
    let lost_cartesian = Crystal::new(
        [
            [1e20, 1.0000000000000002e20, 0.],
            [0., 1., 0.],
            [0., 0., 1.],
        ],
        vec![[0.; 3]],
        vec![26],
    )
    .unwrap();
    let approximate_identity = operation(
        I,
        [3e-6, -2.9999999999999997e-6, 0.],
        SpinChannelAction::Preserve,
    );
    assert!(matches!(
        validate(&lost_cartesian, &[2.], &[approximate_identity]),
        Err(CollinearSpinError::NumericalRange { .. })
    ));
    let two = Crystal::new(CUBIC, vec![[0.; 3], [0.25, 0., 0.]], vec![26, 26]).unwrap();
    let normal = [
        [1.2306185710887074, 2.4608322631435886, 0.],
        [0.24306543930618124, 0.486050908991771, 0.],
    ];
    for (vectors, collinearity) in [
        (normal, 0.),
        (normal.map(|v| v.map(|x| x * scale)), 1e-5),
        ([[2., 1., 0.], [f64::from_bits(1), 0., 0.]], 0.),
    ] {
        assert!(matches!(
            validate_collinear_spin_group(
                &two,
                CollinearMoments::Vectors(&vectors),
                &[CollinearSpinOperation::identity()],
                CollinearSpinTolerances {
                    moment: 0.,
                    collinearity,
                    ..Default::default()
                }
            ),
            Err(CollinearSpinError::NumericalRange { .. })
        ));
    }
    // Equal rounded hypot values are not proof of equal actual magnitudes.
    let s = 2_f64.powi(40);
    let vectors = [
        [s - 1., 2_f64.powi(21), 0.],
        [-(s - 1. + 2_f64.powi(-13)), -2_f64.powi(21), 0.],
    ];
    let ops = [CollinearSpinOperation::identity(), half_swap()];
    let tolerance = CollinearSpinTolerances {
        collinearity: 1.,
        ..Default::default()
    };
    assert!(matches!(
        validate_collinear_spin_group(&afm(), CollinearMoments::Vectors(&vectors), &ops, tolerance),
        Err(CollinearSpinError::NumericalRange {
            context: "comparing vector moment magnitudes"
        })
    ));
    assert!(matches!(
        detect_collinear_spin_group(&afm(), CollinearMoments::Vectors(&vectors), tolerance),
        Err(CollinearSpinError::NumericalRange { .. })
    ));
    validate(&afm(), &[s + 1., -s - 1.], &ops).unwrap(); // Caller-signed scalars are exact inputs.
}

#[test]
fn structural_search_handles_wide_determinants_and_rejects_unrepresentable_rotations() {
    // Both bases reduce to cubic. Their exact unimodular determinant needs
    // cancellation through products that overflow i32 in the old search.
    let small = Crystal::new(
        [[144., 233., 0.], [233., 377., 0.], [0., 0., 1.]],
        vec![[0.; 3]],
        vec![26],
    )
    .unwrap();
    let group =
        detect_collinear_spin_group(&small, CollinearMoments::Signed(&[2.]), Default::default())
            .unwrap();
    assert_eq!(group.operations().len(), 48);
    assert!(
        group
            .operations()
            .iter()
            .any(|g| g.operation.rotation == [[-54289, -87840, 0], [33553, 54289, 0], [0, 0, 1]])
    );
    let big = Crystal::new(
        [[46368., 75025., 0.], [75025., 121393., 0.], [0., 0., 1.]],
        vec![[0.; 3]],
        vec![26],
    )
    .unwrap();
    assert!(matches!(
        detect_collinear_spin_group(&big, CollinearMoments::Signed(&[2.]), Default::default()),
        Err(CollinearSpinError::StructuralSearchFailed)
    ));
}

#[test]
fn validation_retains_input_order_instead_of_inserting_or_normalizing_operations() {
    let operations = [
        half_swap(),
        operation(I, [2., 0., 0.], SpinChannelAction::Preserve),
    ];
    let group = validate(&afm(), &[2., -2.], &operations).unwrap();
    assert_eq!(group.identity_index(), 1);
    assert_eq!(
        group
            .operations()
            .iter()
            .map(|g| g.operation)
            .collect::<Vec<_>>(),
        operations
    );
    assert_eq!(group.multiplication_table()[0][0].cell_shift, [-1, 0, 0]);
    check_certificate(&afm(), &group);
}
