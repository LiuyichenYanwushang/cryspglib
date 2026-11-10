//! Independent regression tests for the S3 cut-over that retired the public
//! `cryspglib::collinear_spin` API.
//!
//! These tests are written only against the public `cryspglib::spin_group`
//! surface.  They do not reference private modules, do not read files and do
//! not use the network: every fixture is constructed in memory.  The retired
//! behaviours that must survive the cut-over are exercised through the new
//! names:
//!
//! * the cubic ferromagnet returns the full finite group (unitary half plus
//!   the anti-unitary `A0` coset), so callers no longer need
//!   `with_spin_only_time_reversal()`;
//! * the half-translation antiferromagnet keeps its integer cell shift and its
//!   `Swap^2 = identity (one cell)` quotient law;
//! * the C4 altermagnet keeps its unitary `Swap` channel and gains the matching
//!   anti-unitary partner without a same-k degeneracy claim;
//! * `CollinearMoments::Signed` and `CollinearMoments::Vectors` describe the
//!   same group, and the oriented axis is only metadata;
//! * `validate_collinear_spin_group` is strict about a missing or corrupted
//!   operation;
//! * all-zero moments, non-collinear vectors, non-finite input and invalid
//!   tolerances are refused with typed errors;
//! * numeric safety for a skew cell, a large integer cell shift and non-finite
//!   input;
//! * one standard-cell exact SSG identification fixture plus negative controls.
//!
//! The certificate checker below is deliberately independent of the library
//! internals: it recomputes products, atom images and residual laws from the
//! raw `rotation`/`translation` fields and only then compares them with the
//! published group tables.

use std::collections::BTreeMap;

use cryspglib::cell::AperiodicAxis;
use cryspglib::spin_group::database::{
    SpinSpaceGroupOperation, collinear_spin_space_group_by_identifier,
};
use cryspglib::spin_group::{
    CollinearAction, CollinearMoments, CollinearSpinGroup, SpinGroupError, SpinGroupKind,
    SpinOperation, SpinSpaceGroupIdentificationError, SpinTolerances, detect_collinear_spin_group,
    identify_collinear_signature, identify_collinear_spin_space_group,
    identify_collinear_spin_space_group_with_tolerance, validate_collinear_spin_group,
};
use cryspglib::{Crystal, OperationKind};

type Mat3 = [[f64; 3]; 3];
type Mat3I = [[i32; 3]; 3];
type Vec3 = [f64; 3];

const I3: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const C4Z: Mat3I = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
const C2Z: Mat3I = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];
const REFLECT_X: Mat3I = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
const CUBIC: Mat3 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const ORTHO: Mat3 = [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]];

fn cubic_ferromagnet() -> Crystal {
    Crystal::new(CUBIC, vec![[0.; 3]], vec![26]).unwrap()
}

fn orthorhombic_afm() -> Crystal {
    Crystal::new(ORTHO, vec![[0.; 3], [0.5, 0., 0.]], vec![26, 26]).unwrap()
}

fn detect(crystal: &Crystal, moments: CollinearMoments<'_>) -> CollinearSpinGroup {
    detect_collinear_spin_group(crystal, moments, SpinTolerances::default())
        .unwrap_or_else(|error| panic!("detection must succeed: {error}"))
}

fn close_at(left: f64, right: f64, epsilon: f64, what: &str) {
    assert!(
        (left - right).abs() <= epsilon,
        "{what}: {left} != {right} (epsilon {epsilon})"
    );
}

fn close_vec(left: Vec3, right: Vec3) {
    for axis in 0..3 {
        close_at(left[axis], right[axis], 1e-10, "vector component");
    }
}

fn cart_norm(lattice: Mat3, fractional: Vec3) -> f64 {
    lattice
        .map(|row| {
            (0..3)
                .map(|column| row[column] * fractional[column])
                .sum::<f64>()
        })
        .into_iter()
        .fold(0., f64::hypot)
}

fn mat_mul(left: Mat3I, right: Mat3I) -> Mat3I {
    std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            (0..3)
                .map(|inner| left[row][inner] * right[inner][column])
                .sum()
        })
    })
}

fn affine_apply(operation: &SpinOperation<CollinearAction>, position: Vec3) -> Vec3 {
    std::array::from_fn(|axis| {
        operation.translation[axis]
            + (0..3)
                .map(|inner| f64::from(operation.rotation[axis][inner]) * position[inner])
                .sum::<f64>()
    })
}

/// `t_left + W_left t_right`, i.e. the raw Seitz translation of the product
/// before the integer representative shift.
fn composed_translation(
    left: &SpinOperation<CollinearAction>,
    right: &SpinOperation<CollinearAction>,
) -> Vec3 {
    std::array::from_fn(|axis| {
        left.translation[axis]
            + (0..3)
                .map(|inner| f64::from(left.rotation[axis][inner]) * right.translation[inner])
                .sum::<f64>()
    })
}

fn flat_rotation(rotation: Mat3I) -> [i32; 9] {
    let mut flat = [0_i32; 9];
    for row in 0..3 {
        for column in 0..3 {
            flat[3 * row + column] = rotation[row][column];
        }
    }
    flat
}

fn operation_key(operation: &SpinOperation<CollinearAction>) -> ([i32; 9], i8, bool) {
    (
        flat_rotation(operation.rotation),
        operation.action.sign(),
        operation.kind.is_antiunitary(),
    )
}

fn same_spatial_action(
    left: &SpinOperation<CollinearAction>,
    right: &SpinOperation<CollinearAction>,
) -> bool {
    left.rotation == right.rotation
        && left.action == right.action
        && (0..3).all(|axis| (left.translation[axis] - right.translation[axis]).abs() <= 1e-12)
}

fn find_operation(
    group: &CollinearSpinGroup,
    rotation: Mat3I,
    translation: Vec3,
    action: CollinearAction,
    kind: OperationKind,
) -> Option<usize> {
    group.operations().iter().position(|operation| {
        operation.rotation == rotation
            && operation.action == action
            && operation.kind == kind
            && (0..3).all(|axis| (operation.translation[axis] - translation[axis]).abs() <= 1e-12)
    })
}

/// True when `raw - candidate` is (within the tolerance) an integer vector of
/// the input cell.  Used only for the orthogonal families in this file.
fn matches_modulo_lattice(lattice: Mat3, raw: Vec3, candidate: Vec3, tolerance: f64) -> bool {
    let mut difference = [0.0_f64; 3];
    for axis in 0..3 {
        let delta = raw[axis] - candidate[axis];
        difference[axis] = delta - delta.round();
    }
    cart_norm(lattice, difference) <= tolerance + 1e-12
}

/// Independent certificate check: recomputes the product table from the raw
/// operations, verifies closure/uniqueness/associativity/cocycles, verifies
/// every atom mapping and verifies that the anti-unitary coset is present.
fn check_certificate(crystal: &Crystal, group: &CollinearSpinGroup) {
    let len = group.len();
    let operations = group.operations();
    let tolerance = group.tolerances();
    let moments = group.order().signed_moments();
    let lattice = crystal.lattice();
    let identity = group.identity_index();

    assert_eq!(group.kind(), SpinGroupKind::Collinear);
    assert_eq!(len, operations.len());
    assert!(!group.is_empty());
    assert_eq!(moments.len(), crystal.natom());
    assert_eq!(group.inverse_indices().len(), len);
    assert!(group.operation(len).is_none());
    assert!(group.atom_mappings(len).is_none());
    assert!(group.product(len, 0).is_none());
    assert!(group.product(0, len).is_none());
    assert!(group.inverse_index(len).is_none());

    let identity_operation = &operations[identity];
    assert_eq!(identity_operation.rotation, I3);
    assert_eq!(identity_operation.action, CollinearAction::Preserve);
    assert_eq!(identity_operation.kind, OperationKind::Unitary);
    // The identity representative is the original input operation: it may
    // carry a lattice-vector translation instead of being normalized to zero.
    for axis in 0..3 {
        close_at(
            identity_operation.translation[axis],
            identity_operation.translation[axis].round(),
            1e-12,
            "identity translation must be an integer lattice vector",
        );
    }
    assert_eq!(group.inverse_index(identity), Some(identity));

    // The retired `with_spin_only_time_reversal()` extension is now implicit:
    // exactly half the operations are anti-unitary, and each unitary operation
    // has one anti-unitary partner with the same spatial part and same action.
    assert_eq!(len % 2, 0, "full collinear group has a unitary+anti half");
    let unitary: Vec<usize> = operations
        .iter()
        .enumerate()
        .filter(|(_, operation)| operation.kind.is_unitary())
        .map(|(index, _)| index)
        .collect();
    let antiunitary: Vec<usize> = operations
        .iter()
        .enumerate()
        .filter(|(_, operation)| operation.kind.is_antiunitary())
        .map(|(index, _)| index)
        .collect();
    assert_eq!(unitary.len(), len / 2);
    assert_eq!(antiunitary.len(), len / 2);
    for &index in &unitary {
        let partners = antiunitary
            .iter()
            .filter(|&&anti| same_spatial_action(&operations[index], &operations[anti]))
            .count();
        assert_eq!(partners, 1, "operation {index} must have one anti partner");
    }

    let mut by_key: BTreeMap<([i32; 9], i8, bool), Vec<usize>> = BTreeMap::new();
    for (index, operation) in operations.iter().enumerate() {
        by_key
            .entry(operation_key(operation))
            .or_default()
            .push(index);
    }

    for left_index in 0..len {
        let left = &operations[left_index];
        assert_eq!(
            group.product(left_index, identity).unwrap().operation_index,
            left_index
        );
        assert_eq!(
            group.product(identity, left_index).unwrap().operation_index,
            left_index
        );
        let inverse = group.inverse_index(left_index).unwrap();
        assert_eq!(
            group.product(left_index, inverse).unwrap().operation_index,
            identity
        );
        assert_eq!(
            group.product(inverse, left_index).unwrap().operation_index,
            identity
        );

        let mapping = group.atom_mappings(left_index).unwrap();
        assert_eq!(mapping.len(), crystal.natom());
        let mut targets: Vec<usize> = mapping.iter().map(|atom| atom.target).collect();
        targets.sort_unstable();
        assert_eq!(targets, (0..crystal.natom()).collect::<Vec<_>>());
        for (atom, mapped) in mapping.iter().enumerate() {
            assert_eq!(crystal.types()[atom], crystal.types()[mapped.target]);
            let expected = affine_apply(left, crystal.positions()[atom]);
            for (axis, expected_axis) in expected.iter().enumerate() {
                let rebuilt = crystal.positions()[mapped.target][axis]
                    + mapped.cell_shift[axis] as f64
                    + mapped.position_residual[axis];
                close_at(*expected_axis, rebuilt, 1e-12, "atom image reconstruction");
            }
            assert!(
                cart_norm(lattice, mapped.position_residual) <= tolerance.position + 1e-12,
                "position residual exceeds tolerance"
            );
            close_at(
                mapped.moment_residual,
                moments[mapped.target] - f64::from(left.action.sign()) * moments[atom],
                1e-12,
                "moment residual",
            );
            assert!(mapped.moment_residual.abs() <= tolerance.moment + 1e-12);
        }

        for right_index in 0..len {
            let right = &operations[right_index];
            let product = group.product(left_index, right_index).unwrap();
            let product_operation = &operations[product.operation_index];
            assert_eq!(
                product_operation.rotation,
                mat_mul(left.rotation, right.rotation),
                "rotation of {left_index} * {right_index}"
            );
            assert_eq!(
                product_operation.action.sign(),
                left.action.sign() * right.action.sign()
            );
            assert_eq!(
                product_operation.kind.is_antiunitary(),
                left.kind.is_antiunitary() ^ right.kind.is_antiunitary()
            );

            let raw = composed_translation(left, right);
            let mut residual = [0.0_f64; 3];
            for axis in 0..3 {
                residual[axis] = raw[axis]
                    - product_operation.translation[axis]
                    - product.cell_shift[axis] as f64;
            }
            assert!(
                cart_norm(lattice, residual) <= tolerance.position + 1e-12,
                "product translation residual for {left_index} * {right_index}"
            );

            // The stored representative must be the unique operation with this
            // rotation/action/kind whose translation equals the raw product.
            let candidates = by_key
                .get(&operation_key(product_operation))
                .expect("product representative is in the group");
            let mut matching = 0_usize;
            for &candidate in candidates {
                if matches_modulo_lattice(
                    lattice,
                    raw,
                    operations[candidate].translation,
                    tolerance.position,
                ) {
                    matching += 1;
                    assert_eq!(
                        candidate, product.operation_index,
                        "non-representative candidate for {left_index} * {right_index}"
                    );
                }
            }
            assert_eq!(
                matching, 1,
                "unique representative for {left_index} * {right_index}"
            );

            let right_mapping = group.atom_mappings(right_index).unwrap();
            let product_mapping = group.atom_mappings(product.operation_index).unwrap();
            for (atom, right_atom) in right_mapping.iter().enumerate() {
                let first = mapping[right_atom.target];
                let combined = product_mapping[atom];
                assert_eq!(first.target, combined.target);
                for (axis, first_shift) in first.cell_shift.iter().enumerate() {
                    assert_eq!(
                        *first_shift
                            + (0..3)
                                .map(|inner| {
                                    i64::from(left.rotation[axis][inner])
                                        * right_atom.cell_shift[inner]
                                })
                                .sum::<i64>(),
                        combined.cell_shift[axis] + product.cell_shift[axis],
                        "cell shift composition for {left_index} * {right_index}"
                    );
                    let composed = first.position_residual[axis]
                        + (0..3)
                            .map(|inner| {
                                f64::from(left.rotation[axis][inner])
                                    * right_atom.position_residual[inner]
                            })
                            .sum::<f64>();
                    let expected_composed = combined.position_residual[axis] + residual[axis];
                    assert!(
                        (composed - expected_composed).abs() <= 1e-12,
                        "position residual composition at ({left_index},{right_index},{atom}, axis {axis})"
                    );
                }
                close_at(
                    first.moment_residual
                        + f64::from(left.action.sign()) * right_atom.moment_residual,
                    combined.moment_residual,
                    1e-12,
                    "moment residual composition",
                );
            }
        }
    }

    // Every row of the multiplication table is a permutation of the group.
    for left_index in 0..len {
        let mut seen = vec![false; len];
        for right_index in 0..len {
            let index = group
                .product(left_index, right_index)
                .unwrap()
                .operation_index;
            assert!(!seen[index], "product map is not a bijection");
            seen[index] = true;
        }
        assert!(seen.iter().all(|&value| value));
    }

    // Associativity and the integer cocycle q(g,h) + q(gh,k) = W_g q(h,k) + q(g,hk).
    for (first, first_operation) in operations.iter().enumerate() {
        for second in 0..len {
            let first_second = group.product(first, second).unwrap();
            for third in 0..len {
                let left_product = group.product(first_second.operation_index, third).unwrap();
                let second_third = group.product(second, third).unwrap();
                let right_product = group.product(first, second_third.operation_index).unwrap();
                assert_eq!(
                    left_product.operation_index, right_product.operation_index,
                    "associativity at ({first},{second},{third})"
                );
                for axis in 0..3 {
                    let lhs = first_second.cell_shift[axis] + left_product.cell_shift[axis];
                    let rhs = right_product.cell_shift[axis]
                        + (0..3)
                            .map(|inner| {
                                i64::from(first_operation.rotation[axis][inner])
                                    * second_third.cell_shift[inner]
                            })
                            .sum::<i64>();
                    assert_eq!(lhs, rhs, "cocycle at ({first},{second},{third})");
                }
            }
        }
    }
}

/// The 48 signed permutation matrices, generated independently of the library.
fn signed_permutations_cubic() -> Vec<Mat3I> {
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let signs = [
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [1, -1, -1],
        [-1, 1, 1],
        [-1, 1, -1],
        [-1, -1, 1],
        [-1, -1, -1],
    ];
    let mut matrices = Vec::new();
    for permutation in permutations {
        for sign in signs {
            let mut matrix = [[0_i32; 3]; 3];
            for column in 0..3 {
                matrix[permutation[column]][column] = sign[column];
            }
            matrices.push(matrix);
        }
    }
    matrices
}

fn assert_same_spatial_group(left: &CollinearSpinGroup, right: &CollinearSpinGroup) {
    assert_eq!(left.len(), right.len());
    for (index, (left_operation, right_operation)) in
        left.operations().iter().zip(right.operations()).enumerate()
    {
        assert_eq!(
            left_operation.rotation, right_operation.rotation,
            "op {index}"
        );
        for axis in 0..3 {
            close_at(
                left_operation.translation[axis],
                right_operation.translation[axis],
                1e-12,
                "operation translation",
            );
        }
        assert_eq!(left_operation.action, right_operation.action);
        assert_eq!(left_operation.kind, right_operation.kind);
    }
}

#[test]
fn cubic_ferromagnet_returns_the_full_unitary_and_antiunitary_group() {
    let crystal = cubic_ferromagnet();
    let group = detect(&crystal, CollinearMoments::Signed(&[2.]));
    assert_eq!(group.kind(), SpinGroupKind::Collinear);
    // 48 spatial operations (m-3m) times the unitary+anti-unitary halves.
    assert_eq!(group.len(), 96);
    let unitary = group
        .operations()
        .iter()
        .filter(|operation| operation.kind.is_unitary())
        .count();
    assert_eq!((unitary, group.len() - unitary), (48, 48));
    assert!(
        group
            .operations()
            .iter()
            .all(|operation| operation.action == CollinearAction::Preserve),
        "a ferromagnet has no channel-swapping operation"
    );

    // The unitary spatial parts are exactly the 48 signed permutations.
    let mut expected: Vec<[i32; 9]> = signed_permutations_cubic()
        .into_iter()
        .map(flat_rotation)
        .collect();
    expected.sort_unstable();
    let mut actual: Vec<[i32; 9]> = group
        .operations()
        .iter()
        .filter(|operation| operation.kind.is_unitary())
        .map(|operation| flat_rotation(operation.rotation))
        .collect();
    actual.sort_unstable();
    assert_eq!(actual, expected);

    // The anti-unitary half is the A0 coset: the channel-preserving
    // [C2_perp T || E] representative is the identity spatial operation.
    let a0 = find_operation(
        &group,
        I3,
        [0.; 3],
        CollinearAction::Preserve,
        OperationKind::Antiunitary,
    )
    .expect("A0 = [C2_perp T || E] is present");
    assert_eq!(group.operations()[a0].rotation, I3);
    assert_eq!(group.operations()[a0].kind, OperationKind::Antiunitary);
    let square = group.product(a0, a0).unwrap();
    assert_eq!(square.operation_index, group.identity_index());
    assert_eq!(square.cell_shift, [0, 0, 0]);
    assert_eq!(square.translation_residual, [0.; 3]);

    // Every operation maps the single atom onto itself at zero residual.
    for index in 0..group.len() {
        let mapping = group.atom_mappings(index).unwrap();
        assert_eq!(mapping.len(), 1);
        assert_eq!(mapping[0].target, 0);
        assert_eq!(mapping[0].cell_shift, [0, 0, 0]);
        close_vec(mapping[0].position_residual, [0.; 3]);
        close_at(
            mapping[0].moment_residual,
            0.0,
            0.0,
            "ferromagnet moment residual",
        );
    }

    check_certificate(&crystal, &group);
}

#[test]
fn half_translation_antiferromagnet_keeps_the_integer_cell_shift() {
    let crystal = orthorhombic_afm();
    let group = detect(&crystal, CollinearMoments::Signed(&[2., -2.]));
    assert_eq!(group.len(), 32);
    let unitary = group
        .operations()
        .iter()
        .filter(|operation| operation.kind.is_unitary())
        .count();
    assert_eq!((unitary, group.len() - unitary), (16, 16));

    let half = find_operation(
        &group,
        I3,
        [0.5, 0., 0.],
        CollinearAction::Swap,
        OperationKind::Unitary,
    )
    .expect("unitary half translation with Swap action");
    let mapping = group.atom_mappings(half).unwrap();
    assert_eq!(
        mapping
            .iter()
            .map(|atom| (atom.target, atom.cell_shift))
            .collect::<Vec<_>>(),
        vec![(1, [0, 0, 0]), (0, [1, 0, 0])]
    );
    for atom in mapping {
        close_vec(atom.position_residual, [0.; 3]);
        close_at(atom.moment_residual, 0.0, 1e-12, "AFM moment residual");
    }

    // Swap^2 is the unitary identity quotient class but moves by one full cell.
    let square = group.product(half, half).unwrap();
    assert_eq!(square.operation_index, group.identity_index());
    assert_eq!(square.cell_shift, [1, 0, 0]);
    close_vec(square.translation_residual, [0.; 3]);
    let raw = composed_translation(&group.operations()[half], &group.operations()[half]);
    for (axis, raw_component) in raw.iter().enumerate() {
        close_at(
            *raw_component,
            group.operations()[group.identity_index()].translation[axis]
                + square.cell_shift[axis] as f64
                + square.translation_residual[axis],
            1e-12,
            "half-translation square reconstruction",
        );
    }

    // The anti-unitary partner has the same spatial part; its square is the
    // same unitary cell-shift class, not a second copy of the translation.
    let anti = find_operation(
        &group,
        I3,
        [0.5, 0., 0.],
        CollinearAction::Swap,
        OperationKind::Antiunitary,
    )
    .expect("anti-unitary half translation");
    let anti_square = group.product(anti, anti).unwrap();
    let anti_square_operation = &group.operations()[anti_square.operation_index];
    assert_eq!(anti_square.operation_index, group.identity_index());
    assert_eq!(anti_square_operation.kind, OperationKind::Unitary);
    assert_eq!(anti_square_operation.action, CollinearAction::Preserve);
    assert_eq!(anti_square.cell_shift, [1, 0, 0]);

    // The 16 unitary operations split into 8 Preserve and 8 Swap channels.
    let swap = group
        .operations()
        .iter()
        .filter(|operation| {
            operation.kind.is_unitary() && operation.action == CollinearAction::Swap
        })
        .count();
    assert_eq!((swap, unitary - swap), (8, 8));

    check_certificate(&crystal, &group);
}

#[test]
fn c4_altermagnet_keeps_the_swap_channel_and_its_anti_partner() {
    let crystal = Crystal::new(
        CUBIC,
        vec![[0.5, 0., 0.], [0., 0.5, 0.], [0.; 3]],
        vec![26, 26, 8],
    )
    .unwrap();
    let group = detect(&crystal, CollinearMoments::Signed(&[2., -2., 0.]));
    assert_eq!(group.len(), 32);
    let unitary = group
        .operations()
        .iter()
        .filter(|operation| operation.kind.is_unitary())
        .count();
    assert_eq!((unitary, group.len() - unitary), (16, 16));

    let c4 = find_operation(
        &group,
        C4Z,
        [0.; 3],
        CollinearAction::Swap,
        OperationKind::Unitary,
    )
    .expect("unitary C4 rotation with Swap action");
    let mapping = group.atom_mappings(c4).unwrap();
    assert_eq!(
        mapping
            .iter()
            .map(|atom| (atom.target, atom.cell_shift))
            .collect::<Vec<_>>(),
        vec![(1, [0, 0, 0]), (0, [-1, 0, 0]), (2, [0, 0, 0])]
    );
    for atom in mapping {
        close_vec(atom.position_residual, [0.; 3]);
        close_at(
            atom.moment_residual,
            0.0,
            1e-12,
            "altermagnet moment residual",
        );
    }

    // C4^2 = C2_z with the Preserve channel; the nonmagnetic marker excludes
    // the otherwise allowed (1/2,1/2,0) Swap translation.
    let square = group.product(c4, c4).unwrap();
    let expected_c2 = find_operation(
        &group,
        C2Z,
        [0.; 3],
        CollinearAction::Preserve,
        OperationKind::Unitary,
    )
    .expect("preserving C2_z");
    assert_eq!(square.operation_index, expected_c2);
    assert_eq!(
        group.operations()[square.operation_index].kind,
        OperationKind::Unitary
    );
    assert_eq!(
        group.operations()[square.operation_index].action,
        CollinearAction::Preserve
    );
    assert_eq!(square.cell_shift, [0, 0, 0]);
    assert!(
        find_operation(
            &group,
            I3,
            [0.; 3],
            CollinearAction::Swap,
            OperationKind::Unitary
        )
        .is_none(),
        "the nonmagnetic marker must exclude the identity-rotation Swap"
    );
    assert!(
        find_operation(
            &group,
            I3,
            [0.; 3],
            CollinearAction::Swap,
            OperationKind::Antiunitary
        )
        .is_none()
    );

    // The C4 channel is a genuine Swap and has an anti-unitary partner.
    assert_eq!(group.operations()[c4].action, CollinearAction::Swap);
    let anti_c4 = find_operation(
        &group,
        C4Z,
        [0.; 3],
        CollinearAction::Swap,
        OperationKind::Antiunitary,
    )
    .expect("anti-unitary C4 partner");
    assert_eq!(
        group.atom_mappings(anti_c4).unwrap()[0].target,
        group.atom_mappings(c4).unwrap()[0].target
    );

    check_certificate(&crystal, &group);
}

#[test]
fn signed_and_vector_moments_describe_the_same_group_and_axis_is_metadata() {
    let crystal = orthorhombic_afm();
    let signed = detect(&crystal, CollinearMoments::Signed(&[2., -2.]));
    assert_eq!(signed.order().signed_moments(), &[2., -2.]);
    assert!(signed.order().axis().is_none());
    assert_eq!(signed.order().collinearity_residuals(), &[0., 0.]);

    for (vectors, expected_axis) in [
        ([[0., 0., 2.], [0., 0., -2.]], [0., 0., 1.]),
        ([[2., 0., 0.], [-2., 0., 0.]], [1., 0., 0.]),
        ([[0., 2., 0.], [0., -2., 0.]], [0., 1., 0.]),
        ([[1., 1., 1.], [-1., -1., -1.]], [1. / 3_f64.sqrt(); 3]),
    ] {
        let vector_group = detect(&crystal, CollinearMoments::Vectors(&vectors));
        assert_same_spatial_group(&signed, &vector_group);
        for index in 0..signed.len() {
            let signed_mapping = signed.atom_mappings(index).unwrap();
            let vector_mapping = vector_group.atom_mappings(index).unwrap();
            for (left, right) in signed_mapping.iter().zip(vector_mapping) {
                assert_eq!(left.target, right.target);
                assert_eq!(left.cell_shift, right.cell_shift);
                close_vec(left.position_residual, right.position_residual);
            }
        }
        let axis = vector_group
            .order()
            .axis()
            .expect("vector input has an axis");
        close_vec(axis, expected_axis);
        let residuals = vector_group.order().collinearity_residuals();
        assert!(residuals.iter().all(|residual| *residual <= 1e-12));
        let magnitudes: Vec<f64> = vector_group
            .order()
            .signed_moments()
            .iter()
            .map(|moment| moment.abs())
            .collect();
        close_at(magnitudes[0], magnitudes[1], 1e-12, "AFM vector magnitudes");
        close_at(
            magnitudes[0],
            (vectors[0][0] * vectors[0][0]
                + vectors[0][1] * vectors[0][1]
                + vectors[0][2] * vectors[0][2])
                .sqrt(),
            1e-12,
            "vector norm",
        );
        check_certificate(&crystal, &vector_group);
    }
}

#[test]
fn strict_validation_refuses_missing_or_corrupted_operations() {
    let crystal = orthorhombic_afm();
    let moments = [2., -2.];
    let detected = detect(&crystal, CollinearMoments::Signed(&moments));
    let validated = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        detected.operations(),
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(validated.len(), detected.len());
    assert_eq!(validated.operations(), detected.operations());
    assert_eq!(validated.tolerances(), SpinTolerances::default());

    // Dropping an anti-unitary operation breaks the equal-half contract.
    let mut missing_anti = detected.operations().to_vec();
    missing_anti.pop();
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &missing_anti,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SpinGroupError::KindCount {
            unitary: 16,
            antiunitary: 15
        }
    ));

    // Dropping a unitary operation is refused symmetrically.
    let mut missing_unitary = detected.operations().to_vec();
    missing_unitary.remove(0);
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &missing_unitary,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SpinGroupError::KindCount {
            unitary: 15,
            antiunitary: 16
        }
    ));

    // Corrupting the channel of an anti-unitary operation is refused.
    let mut wrong_action = detected.operations().to_vec();
    let anti = wrong_action
        .iter()
        .position(|operation| operation.kind.is_antiunitary())
        .unwrap();
    wrong_action[anti].action = match wrong_action[anti].action {
        CollinearAction::Preserve => CollinearAction::Swap,
        CollinearAction::Swap => CollinearAction::Preserve,
    };
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &wrong_action,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AntiunitaryPartner { .. }));

    // Corrupting the translation of an anti-unitary operation is refused even
    // when the spatial rotation and action are unchanged.
    let mut wrong_translation = detected.operations().to_vec();
    let anti = wrong_translation
        .iter()
        .position(|operation| operation.kind.is_antiunitary())
        .unwrap();
    wrong_translation[anti].translation[0] += 0.25;
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &wrong_translation,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AntiunitaryPartner { .. }));

    // A hand-built set without the identity is refused by the engine
    // validator, not silently completed.
    let half = detected
        .operations()
        .iter()
        .find(|operation| {
            operation.rotation == I3
                && operation.action == CollinearAction::Swap
                && (operation.translation[0] - 0.5).abs() < 1e-12
                && operation.kind.is_unitary()
        })
        .copied()
        .unwrap();
    let anti_half = detected
        .operations()
        .iter()
        .find(|operation| {
            operation.rotation == I3
                && operation.action == CollinearAction::Swap
                && (operation.translation[0] - 0.5).abs() < 1e-12
                && operation.kind.is_antiunitary()
        })
        .copied()
        .unwrap();
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &[half, anti_half],
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::MissingIdentity));

    // A duplicated unitary operation is refused, not deduplicated.
    let mut duplicated = detected.operations().to_vec();
    duplicated.push(detected.operations()[0]);
    duplicated.push(
        detected
            .operations()
            .iter()
            .find(|operation| {
                operation.kind.is_antiunitary()
                    && same_spatial_action(operation, &detected.operations()[0])
            })
            .copied()
            .unwrap(),
    );
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &duplicated,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        SpinGroupError::Duplicate {
            first_index: 0,
            duplicate_index: 16
        }
    ));

    // A unitary operation that maps an atom to no same-species site is
    // refused with a witness instead of being dropped.
    let mut wrong_atom = detected.operations().to_vec();
    let unitary_half = wrong_atom
        .iter()
        .position(|operation| {
            operation.kind.is_unitary()
                && operation.rotation == I3
                && operation.action == CollinearAction::Swap
        })
        .unwrap();
    wrong_atom[unitary_half].translation[0] = 0.25;
    let error = validate_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&moments),
        &wrong_atom,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AtomMappingMissing { .. }));
}

#[test]
fn all_zero_noncollinear_and_invalid_inputs_are_refused() {
    let crystal = orthorhombic_afm();
    let default = SpinTolerances::default();

    for moments in [
        CollinearMoments::Signed(&[0., 0.]),
        CollinearMoments::Signed(&[1e-6, -1e-6]),
        CollinearMoments::Vectors(&[[0.; 3]; 2]),
    ] {
        let error = detect_collinear_spin_group(&crystal, moments, default).unwrap_err();
        assert!(matches!(error, SpinGroupError::AllZeroMoments));
    }

    let noncollinear = [[1., 0., 0.], [0., 1., 0.]];
    let error =
        detect_collinear_spin_group(&crystal, CollinearMoments::Vectors(&noncollinear), default)
            .unwrap_err();
    assert!(
        matches!(error, SpinGroupError::NonCollinear { .. }),
        "expected NonCollinear, got {error}"
    );

    let nearly_collinear = [[0., 0., 2.], [0.005, 0., -1.]];
    let error = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Vectors(&nearly_collinear),
        SpinTolerances {
            direction: 1e-4,
            ..default
        },
    )
    .unwrap_err();
    assert!(
        matches!(error, SpinGroupError::NonCollinear { .. }),
        "expected NonCollinear, got {error}"
    );

    for (tolerances, expected_name) in [
        (
            SpinTolerances {
                position: 0.0,
                ..default
            },
            "position",
        ),
        (
            SpinTolerances {
                moment: -1.0,
                ..default
            },
            "moment",
        ),
        (
            SpinTolerances {
                direction: f64::NAN,
                ..default
            },
            "direction",
        ),
        (
            SpinTolerances {
                direction: f64::INFINITY,
                ..default
            },
            "direction",
        ),
    ] {
        let error =
            detect_collinear_spin_group(&crystal, CollinearMoments::Signed(&[2., -2.]), tolerances)
                .unwrap_err();
        match error {
            SpinGroupError::InvalidTolerance { name, .. } => assert_eq!(name, expected_name),
            other => panic!("expected InvalidTolerance, got {other}"),
        }
    }

    for moments in [
        CollinearMoments::Signed(&[2., f64::INFINITY]),
        CollinearMoments::Signed(&[f64::NAN, -2.]),
        CollinearMoments::Vectors(&[[0., 0., f64::NAN], [0.; 3]]),
    ] {
        let error = detect_collinear_spin_group(&crystal, moments, default).unwrap_err();
        assert!(matches!(error, SpinGroupError::InvalidInput { .. }));
    }

    for broken in [
        Crystal::new([[0.; 3]; 3], vec![[0.; 3]], vec![26]).unwrap(),
        Crystal::new(CUBIC, vec![[f64::NAN, 0., 0.]], vec![26]).unwrap(),
        Crystal::new(CUBIC, vec![], vec![]).unwrap(),
        Crystal::new(CUBIC, vec![[0.; 3]], vec![26])
            .unwrap()
            .with_layer(AperiodicAxis::Z),
    ] {
        let error = detect_collinear_spin_group(&broken, CollinearMoments::Signed(&[2.]), default)
            .unwrap_err();
        assert!(matches!(error, SpinGroupError::InvalidInput { .. }));
    }

    let error = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2.]),
        SpinTolerances {
            position: 1e5,
            ..default
        },
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::ImageSearchLimit));
}

#[test]
fn numeric_safety_skew_cell_large_integer_shift_and_non_finite_input() {
    // Skew cell: the only lattice image within tolerance needs the integer
    // shift (5,0,0); componentwise rounding of the fractional translation would
    // give (0,0,0) and would be wrong.
    let skew_lattice = [[1., 10., 0.], [0., 0.1, 0.], [0., 0., 1.]];
    let skew_crystal = Crystal::new(skew_lattice, vec![[0.; 3]], vec![26]).unwrap();
    let unitary_skew = SpinOperation {
        rotation: I3,
        translation: [0.2, 0.49, 0.],
        action: CollinearAction::Preserve,
        kind: OperationKind::Unitary,
    };
    let anti_skew = SpinOperation {
        kind: OperationKind::Antiunitary,
        ..unitary_skew
    };
    let skew_group = validate_collinear_spin_group(
        &skew_crystal,
        CollinearMoments::Signed(&[2.]),
        &[unitary_skew, anti_skew],
        SpinTolerances {
            position: 0.112,
            ..SpinTolerances::default()
        },
    )
    .unwrap();
    assert_eq!(skew_group.len(), 2);
    let skew_mapping = skew_group.atom_mappings(0).unwrap();
    assert_eq!(skew_mapping.len(), 1);
    assert_eq!(skew_mapping[0].cell_shift, [5, 0, 0]);
    assert_eq!(skew_mapping[0].target, 0);
    close_at(
        skew_mapping[0].position_residual[0],
        -4.8,
        1e-12,
        "skew residual x",
    );
    close_at(
        skew_mapping[0].position_residual[1],
        0.49,
        1e-12,
        "skew residual y",
    );
    let expected_cart = (0.1_f64 * 0.1 + 0.049 * 0.049).sqrt();
    close_at(
        cart_norm(skew_lattice, skew_mapping[0].position_residual),
        expected_cart,
        1e-12,
        "skew Cartesian residual",
    );
    assert_eq!(
        skew_group.atom_mappings(1).unwrap()[0].cell_shift,
        skew_mapping[0].cell_shift,
        "the anti partner shares the spatial mapping"
    );
    // A slightly larger tolerance admits a second lattice image and the
    // validator refuses with AmbiguousPeriodicImage instead of choosing one.
    let error = validate_collinear_spin_group(
        &skew_crystal,
        CollinearMoments::Signed(&[2.]),
        &[unitary_skew, anti_skew],
        SpinTolerances {
            position: 0.114,
            ..SpinTolerances::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AmbiguousPeriodicImage));

    // Large integer shift: an atom at 2^39 with a lattice translation of
    // -2^39 must produce i64 cell shifts and products without i32 narrowing.
    // 2^32 and 2^-20 are both exactly representable, so the sub-tolerance
    // residual is not lost to the large-representative rounding of 1e-6.
    let large = 2_f64.powi(32);
    let tiny = 2_f64.powi(-20);
    let shifted_crystal = Crystal::new(CUBIC, vec![[large, 0., 0.]], vec![26]).unwrap();
    let identity_operation = SpinOperation {
        rotation: I3,
        translation: [-large, 0., 0.],
        action: CollinearAction::Preserve,
        kind: OperationKind::Unitary,
    };
    let reflection_operation = SpinOperation {
        rotation: REFLECT_X,
        translation: [large + tiny, 0., 0.],
        action: CollinearAction::Preserve,
        kind: OperationKind::Unitary,
    };
    let target = [
        identity_operation,
        reflection_operation,
        SpinOperation {
            kind: OperationKind::Antiunitary,
            ..identity_operation
        },
        SpinOperation {
            kind: OperationKind::Antiunitary,
            ..reflection_operation
        },
    ];
    let large_group = validate_collinear_spin_group(
        &shifted_crystal,
        CollinearMoments::Signed(&[2.]),
        &target,
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(large_group.len(), 4);
    let shift = -(1_i64 << 32);
    assert_eq!(
        large_group.atom_mappings(0).unwrap()[0].cell_shift,
        [shift, 0, 0]
    );
    assert_eq!(
        large_group.atom_mappings(1).unwrap()[0].cell_shift,
        [shift, 0, 0]
    );
    close_at(
        large_group.atom_mappings(1).unwrap()[0].position_residual[0],
        tiny,
        1e-15,
        "large-shift reflection residual",
    );
    assert_eq!(
        large_group.product(1, 1).unwrap().cell_shift,
        [1_i64 << 32, 0, 0]
    );
    assert_eq!(large_group.product(0, 0).unwrap().cell_shift, [shift, 0, 0]);
    // A full naive-float certificate check would lose the 2^-20 residual in
    // `2^32 + 2^-20 + 2^32` (the sum rounds to 2^33), so verify the exact
    // group laws that matter for this extreme representative directly.
    let kinds: Vec<bool> = large_group
        .operations()
        .iter()
        .map(|operation| operation.kind.is_antiunitary())
        .collect();
    assert_eq!(kinds, vec![false, false, true, true]);
    assert!(
        large_group
            .operations()
            .iter()
            .all(|operation| operation.action == CollinearAction::Preserve)
    );
    assert_eq!(large_group.identity_index(), 0);
    for index in 0..4 {
        let inverse = large_group.inverse_index(index).unwrap();
        assert_eq!(
            large_group.product(index, inverse).unwrap().operation_index,
            0
        );
        assert_eq!(
            large_group.product(inverse, index).unwrap().operation_index,
            0
        );
        let mut seen = [false; 4];
        for right in 0..4 {
            let product = large_group.product(index, right).unwrap().operation_index;
            assert!(product < 4);
            assert!(!seen[product]);
            seen[product] = true;
        }
        assert!(seen.iter().all(|&value| value));
    }
    assert_eq!(large_group.product(1, 0).unwrap().operation_index, 1);
    assert_eq!(
        large_group.product(1, 0).unwrap().cell_shift,
        [1_i64 << 32, 0, 0]
    );
    assert_eq!(large_group.product(0, 1).unwrap().operation_index, 1);
    assert_eq!(large_group.product(0, 1).unwrap().cell_shift, [shift, 0, 0]);
    assert_eq!(large_group.product(1, 2).unwrap().operation_index, 3);
    assert_eq!(large_group.product(2, 3).unwrap().operation_index, 1);

    // A residual larger than the position tolerance is refused.
    let mut too_far = target;
    too_far[1].translation[0] = large + 1e-4;
    too_far[3].translation[0] = large + 1e-4;
    let error = validate_collinear_spin_group(
        &shifted_crystal,
        CollinearMoments::Signed(&[2.]),
        &too_far,
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::AtomMappingMissing { .. }));

    // Non-finite operations are refused before any arithmetic.
    for bad in [f64::NAN, f64::INFINITY] {
        let broken = SpinOperation {
            rotation: I3,
            translation: [bad, 0., 0.],
            action: CollinearAction::Preserve,
            kind: OperationKind::Unitary,
        };
        let error = validate_collinear_spin_group(
            &large_crystal(),
            CollinearMoments::Signed(&[2.]),
            &[
                broken,
                SpinOperation {
                    kind: OperationKind::Antiunitary,
                    ..broken
                },
            ],
            SpinTolerances::default(),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            SpinGroupError::InvalidTranslation { operation_index: 0 }
        ));
    }
    let error = detect_collinear_spin_group(
        &Crystal::new(CUBIC, vec![[f64::NAN, 0., 0.]], vec![26]).unwrap(),
        CollinearMoments::Signed(&[2.]),
        SpinTolerances::default(),
    )
    .unwrap_err();
    assert!(matches!(error, SpinGroupError::InvalidInput { .. }));
}

fn large_crystal() -> Crystal {
    Crystal::new(CUBIC, vec![[0.; 3]], vec![26]).unwrap()
}

#[test]
fn standard_cell_identification_fixture_and_negative_controls() {
    let crystal = Crystal::new(ORTHO, vec![[0.; 3], [0.5, 0., 0.]], vec![26, 26])
        .unwrap()
        .with_magnetic(vec![[0., 0., 2.], [0., 0., -2.]])
        .unwrap();
    let group = detect(&crystal, CollinearMoments::Signed(&[2., -2.]));
    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "47.47.2.1.L");
    assert_eq!(identification.international_symbol(), "P_a1m1m1m∞m1");
    assert_eq!(identification.number.g0(), 47);
    assert_eq!(identification.number.l0(), 47);
    assert_eq!(identification.number.k_index(), 2);
    assert_eq!(identification.number.map_index(), 1);
    assert_eq!(identification.number.letter(), 'L');
    assert_eq!(identification.record.identifier, "47.47.2.1.L");
    assert_eq!(identification.record.basis, "L0");
    assert_eq!(
        identification.record.configuration,
        SpinGroupKind::Collinear
    );
    // The full detected group is exactly the unitary+anti double of the
    // pinned quotient signature.
    assert_eq!(identification.record.operations.len(), 16);
    assert_eq!(group.len(), 2 * identification.record.operations.len());

    // Identifying the pinned signature directly reproduces the same number.
    let direct = identify_collinear_signature(identification.record.operations).unwrap();
    assert_eq!(direct.number, identification.number);
    assert_eq!(direct.record.identifier, identification.record.identifier);

    // A global spin rotation and a lattice-vector origin shift do not change
    // the collinear number (no oriented magnetic axis).
    for moments in [
        CollinearMoments::Vectors(&[[0., 0., 1.]]),
        CollinearMoments::Vectors(&[[1., 1., 1.]]),
        CollinearMoments::Signed(&[2.]),
    ] {
        let rotated = detect(&cubic_ferromagnet(), moments);
        assert_eq!(
            identify_collinear_spin_space_group(&rotated)
                .unwrap()
                .identifier(),
            "221.221.1.1.L"
        );
    }
    let shifted_origin = Crystal::new(CUBIC, vec![[1., 0., 0.]], vec![26]).unwrap();
    let shifted = detect(&shifted_origin, CollinearMoments::Signed(&[2.]));
    assert_eq!(
        identify_collinear_spin_space_group(&shifted)
            .unwrap()
            .identifier(),
        "221.221.1.1.L"
    );

    // A supercell of the same magnetic structure is not the database standard
    // cell: identification must refuse instead of searching for a cell change.
    let supercell = Crystal::new(
        [[4., 0., 0.], [0., 3., 0.], [0., 0., 4.]],
        vec![[0.; 3], [0.25, 0., 0.], [0.5, 0., 0.], [0.75, 0., 0.]],
        vec![26; 4],
    )
    .unwrap();
    let supercell_group = detect(&supercell, CollinearMoments::Signed(&[2., -2., 2., -2.]));
    assert!(matches!(
        identify_collinear_spin_space_group(&supercell_group),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));

    // A non-standard origin (the operation data genuinely changes) is refused.
    let shifted_atom = Crystal::new(CUBIC, vec![[0.1, 0.1, 0.1]], vec![26]).unwrap();
    let shifted_atom_group = detect(&shifted_atom, CollinearMoments::Signed(&[2.]));
    assert!(matches!(
        identify_collinear_spin_space_group(&shifted_atom_group),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));

    // A corrupted pinned signature is refused.
    let mut corrupted = identification.record.operations.to_vec();
    corrupted[0].translation[0] = 0.25;
    assert!(matches!(
        identify_collinear_signature(&corrupted),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));

    // Identification errors carry their typed variants.
    assert!(matches!(
        identify_collinear_signature(&[]),
        Err(SpinSpaceGroupIdentificationError::EmptyGroup)
    ));
    let non_finite = [SpinSpaceGroupOperation {
        rotation: I3,
        translation: [f64::NAN, 0., 0.],
        channel: CollinearAction::Preserve,
    }];
    assert!(matches!(
        identify_collinear_signature(&non_finite),
        Err(SpinSpaceGroupIdentificationError::NonFiniteTranslation { .. })
    ));
    let conflicting = [
        SpinSpaceGroupOperation {
            rotation: I3,
            translation: [0., 0., 0.],
            channel: CollinearAction::Preserve,
        },
        SpinSpaceGroupOperation {
            rotation: I3,
            translation: [1e-9, 0., 0.],
            channel: CollinearAction::Swap,
        },
    ];
    assert!(matches!(
        identify_collinear_signature(&conflicting),
        Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
    ));
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY, 0.5] {
        assert!(matches!(
            identify_collinear_spin_space_group_with_tolerance(&group, tolerance),
            Err(SpinSpaceGroupIdentificationError::InvalidTolerance { .. })
        ));
    }

    // The documented lower-level lookup still finds the same pinned record.
    assert_eq!(
        collinear_spin_space_group_by_identifier("47.47.2.1.L")
            .unwrap()
            .identifier,
        "47.47.2.1.L"
    );
}
