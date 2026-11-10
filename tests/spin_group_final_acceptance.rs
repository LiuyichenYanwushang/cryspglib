//! Final independent acceptance tests for the collinear spin-group work.
//!
//! These tests were written by the final acceptance reviewer against the
//! *public* `cryspglib::spin_group` surface only.  They read no files, use no
//! network and do not import any private module: every expected value below is
//! either a pinned database constant or an independent construction from the
//! raw operation data.
//!
//! Acceptance surface covered:
//!
//! * pinned database count, release provenance and identifier format;
//! * the pinned table's structural invariants and independently counted
//!   aggregates (1,421 records / 29,457 quotient operations);
//! * all 1,421 pinned quotient signatures identifying themselves through the
//!   public lower-level entry point;
//! * detection -> identification end-to-end for the cubic ferromagnet, the
//!   orthorhombic two-atom antiferromagnet and the C4 altermagnet fixtures
//!   already used by the stage tests;
//! * the full group is exactly two halves (unitary + anti-unitary) of equal
//!   size, and unitary filtering reproduces the pinned quotient signature;
//! * the identification tolerance contract and negative controls;
//! * the exact values asserted by the `docs/spin-group.md` identification
//!   example;
//! * cross-cutting acceptance: the full-group certificate (group law, inverse,
//!   atom mappings, action/kind composition), `validate_collinear_spin_group`
//!   round-trip plus its refusal controls, and the documented fact that the
//!   anti-unitary `kind` is not part of the identification signature.

use std::collections::BTreeSet;

use cryspglib::Crystal;
use cryspglib::OperationKind;
use cryspglib::spin_group::database::{
    SpinSpaceGroupNumber, SpinSpaceGroupOperation, SpinSpaceGroupRecord,
    collinear_spin_space_group, collinear_spin_space_group_by_identifier,
    collinear_spin_space_groups, collinear_ssg_identifiers_sha256, collinear_ssg_manifest_sha256,
    collinear_ssg_record_count, collinear_ssg_release_id, collinear_ssg_release_name,
    collinear_ssg_release_source_sha256,
};
use cryspglib::spin_group::{
    CollinearAction, CollinearMoments, CollinearSpinGroup,
    DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE, SpinGroupError, SpinGroupKind, SpinOperation,
    SpinSpaceGroupIdentificationError, SpinTolerances, detect_collinear_spin_group,
    identify_collinear_signature, identify_collinear_signature_with_tolerance,
    identify_collinear_spin_space_group, identify_collinear_spin_space_group_with_tolerance,
    validate_collinear_spin_group,
};

type Mat3 = [[f64; 3]; 3];
type Mat3I = [[i32; 3]; 3];
type Vec3 = [f64; 3];

const CUBIC: Mat3 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const ORTHO: Mat3 = [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]];

const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const C4_Z: Mat3I = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
const C2_Z: Mat3I = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];

const RECORD_COUNT: usize = 1_421;
const QUOTIENT_OPERATION_COUNT: usize = 29_457;
const PRESERVE_COUNT: usize = 16_941;
const SWAP_COUNT: usize = 12_516;
const RECORDS_WITH_SWAP: usize = 1_191;
const ZERO_TOLERANCE: f64 = 1e-10;

// ---------------------------------------------------------------------------
// Fixtures (identical to the fixtures used by the stage tests)
// ---------------------------------------------------------------------------

fn cubic_ferromagnet(moment: Vec3) -> Crystal {
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

fn c4_altermagnet() -> Crystal {
    Crystal::new(
        CUBIC,
        vec![[0.5, 0., 0.], [0., 0.5, 0.], [0.; 3]],
        vec![26, 26, 8],
    )
    .unwrap()
    .with_magnetic(vec![[0., 0., 2.], [0., 0., -2.], [0., 0., 0.]])
    .unwrap()
}

fn detect(crystal: &Crystal, moments: CollinearMoments<'_>) -> CollinearSpinGroup {
    detect_collinear_spin_group(crystal, moments, SpinTolerances::default())
        .unwrap_or_else(|error| panic!("detection must succeed: {error}"))
}

fn record(identifier: &str) -> &'static SpinSpaceGroupRecord {
    collinear_spin_space_group_by_identifier(identifier)
        .unwrap_or_else(|| panic!("record {identifier} must exist"))
}

// ---------------------------------------------------------------------------
// Independent helpers
// ---------------------------------------------------------------------------

fn matrix_product(left: Mat3I, right: Mat3I) -> Mat3I {
    std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            (0..3)
                .map(|inner| left[row][inner] * right[inner][column])
                .sum()
        })
    })
}

/// Determinant of a `Mat3I` with enough headroom for signed-permutation or
/// small-integer rotations.
fn rotation_determinant(rotation: &Mat3I) -> i64 {
    let a = i64::from(rotation[0][0]);
    let b = i64::from(rotation[0][1]);
    let c = i64::from(rotation[0][2]);
    let d = i64::from(rotation[1][0]);
    let e = i64::from(rotation[1][1]);
    let f = i64::from(rotation[1][2]);
    let g = i64::from(rotation[2][0]);
    let h = i64::from(rotation[2][1]);
    let i = i64::from(rotation[2][2]);
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
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

fn close_at(left: f64, right: f64, tolerance: f64, what: &str) {
    assert!(
        (left - right).abs() <= tolerance,
        "{what}: {left} != {right} (tolerance {tolerance})"
    );
}

/// Fractional-translation equality modulo the input-cell lattice.
fn translation_close(left: Vec3, right: Vec3, tolerance: f64) -> bool {
    (0..3).all(|axis| {
        let difference = (left[axis] - right[axis]).rem_euclid(1.0);
        difference < tolerance || (1.0 - difference) < tolerance
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
            && translation_close(operation.translation, translation, ZERO_TOLERANCE)
    })
}

fn count_kind(group: &CollinearSpinGroup, want_unitary: bool) -> usize {
    group
        .operations()
        .iter()
        .filter(|operation| operation.kind.is_unitary() == want_unitary)
        .count()
}

/// Quotient signature of a detected group: deduplicate the full unitary +
/// anti-unitary operation list by `(rotation, translation)` and keep the
/// channel action.  This is exactly the object the database stores.
fn quotient_signature(group: &CollinearSpinGroup) -> Vec<SpinSpaceGroupOperation> {
    let mut signature: Vec<SpinSpaceGroupOperation> = Vec::new();
    for operation in group.operations() {
        let candidate = SpinSpaceGroupOperation {
            rotation: operation.rotation,
            translation: operation.translation,
            channel: operation.action,
        };
        match signature.iter().position(|existing| {
            existing.rotation == candidate.rotation
                && translation_close(existing.translation, candidate.translation, ZERO_TOLERANCE)
        }) {
            Some(index) => assert_eq!(
                signature[index].channel, candidate.channel,
                "unitary and anti-unitary partners must share the channel action"
            ),
            None => signature.push(candidate),
        }
    }
    signature
}

fn signature_matches_record(group: &CollinearSpinGroup, record: &SpinSpaceGroupRecord) {
    let signature = quotient_signature(group);
    assert_eq!(
        signature.len(),
        record.operations.len(),
        "{}: quotient-signature length",
        record.identifier
    );
    for expected in record.operations {
        let matches: Vec<_> = signature
            .iter()
            .filter(|actual| {
                actual.rotation == expected.rotation
                    && actual.channel == expected.channel
                    && translation_close(actual.translation, expected.translation, ZERO_TOLERANCE)
            })
            .collect();
        assert_eq!(
            matches.len(),
            1,
            "{}: expected operation {expected:?} must occur exactly once",
            record.identifier
        );
    }
    for actual in &signature {
        let matches = record
            .operations
            .iter()
            .filter(|expected| {
                expected.rotation == actual.rotation
                    && expected.channel == actual.channel
                    && translation_close(expected.translation, actual.translation, ZERO_TOLERANCE)
            })
            .count();
        assert_eq!(
            matches, 1,
            "{}: detected operation {actual:?} must occur in the pinned record",
            record.identifier
        );
    }
}

/// Independent full-group certificate check: kind split, partner bijection,
/// identity/inverse, the complete ordered product law (rotation, action, kind,
/// integer shift and Seitz residual), and every atom mapping.
fn check_full_group_certificate(crystal: &Crystal, group: &CollinearSpinGroup) {
    let operations = group.operations();
    let len = group.len();
    let identity = group.identity_index();
    let tolerance = group.tolerances();
    let lattice = crystal.lattice();
    let moments = group.order().signed_moments();

    assert_eq!(group.kind(), SpinGroupKind::Collinear);
    assert_eq!(len, operations.len());
    assert!(len > 0 && !group.is_empty());
    assert_eq!(len % 2, 0, "full collinear group is unitary + anti-unitary");
    assert_eq!(moments.len(), crystal.natom());
    assert_eq!(group.inverse_indices().len(), len);
    assert!(group.operation(len).is_none());
    assert!(group.atom_mappings(len).is_none());
    assert!(group.product(len, 0).is_none());
    assert!(group.product(0, len).is_none());
    assert!(group.inverse_index(len).is_none());

    // Exactly half unitary, half anti-unitary, with a bijective partner map.
    let unitary = count_kind(group, true);
    let antiunitary = count_kind(group, false);
    assert_eq!(unitary + antiunitary, len);
    assert_eq!(unitary, len / 2);
    assert_eq!(antiunitary, len / 2);
    let unitary_indices: Vec<usize> = operations
        .iter()
        .enumerate()
        .filter(|(_, operation)| operation.kind.is_unitary())
        .map(|(index, _)| index)
        .collect();
    for &index in &unitary_indices {
        let operation = &operations[index];
        let partners: Vec<usize> = operations
            .iter()
            .enumerate()
            .filter(|(_, candidate)| {
                candidate.kind.is_antiunitary()
                    && candidate.rotation == operation.rotation
                    && candidate.action == operation.action
                    && translation_close(
                        candidate.translation,
                        operation.translation,
                        ZERO_TOLERANCE,
                    )
            })
            .map(|(partner, _)| partner)
            .collect();
        assert_eq!(
            partners.len(),
            1,
            "unitary operation {index} must have exactly one anti-unitary partner"
        );
    }

    // Identity and inverses.
    let identity_operation = &operations[identity];
    assert_eq!(identity_operation.rotation, IDENTITY_ROTATION);
    assert_eq!(identity_operation.action, CollinearAction::Preserve);
    assert_eq!(identity_operation.kind, OperationKind::Unitary);
    assert_eq!(group.inverse_index(identity), Some(identity));
    for &index in &unitary_indices {
        let inverse = group.inverse_index(index).expect("unitary inverse");
        assert!(operations[inverse].kind.is_unitary());
        assert_eq!(
            group.product(index, inverse).unwrap().operation_index,
            identity
        );
        assert_eq!(
            group.product(inverse, index).unwrap().operation_index,
            identity
        );
    }

    // Full ordered product law with an independent recomputation.
    for (left_index, left) in operations.iter().enumerate() {
        assert_eq!(
            group.product(left_index, identity).unwrap().operation_index,
            left_index
        );
        assert_eq!(
            group.product(identity, left_index).unwrap().operation_index,
            left_index
        );
        let mapping = group.atom_mappings(left_index).expect("atom mapping");
        assert_eq!(mapping.len(), crystal.natom());
        let mut targets: Vec<usize> = mapping.iter().map(|atom| atom.target).collect();
        targets.sort_unstable();
        assert_eq!(targets, (0..crystal.natom()).collect::<Vec<_>>());
        for (atom_index, atom) in mapping.iter().enumerate() {
            assert_eq!(
                crystal.types()[atom_index],
                crystal.types()[atom.target],
                "atom species must be preserved"
            );
            let image = affine_apply(left, crystal.positions()[atom_index]);
            for (axis, image_axis) in image.iter().enumerate() {
                let rebuilt = crystal.positions()[atom.target][axis]
                    + atom.cell_shift[axis] as f64
                    + atom.position_residual[axis];
                close_at(*image_axis, rebuilt, 1e-12, "atom image reconstruction");
            }
            assert!(
                cart_norm(lattice, atom.position_residual) <= tolerance.position + 1e-12,
                "position residual exceeds the certificate tolerance"
            );
            close_at(
                atom.moment_residual,
                moments[atom.target] - f64::from(left.action.sign()) * moments[atom_index],
                1e-12,
                "moment residual",
            );
            assert!(atom.moment_residual.abs() <= tolerance.moment + 1e-12);
        }

        for (right_index, right) in operations.iter().enumerate() {
            let product = group.product(left_index, right_index).expect("product");
            let product_operation = &operations[product.operation_index];
            assert_eq!(
                product_operation.rotation,
                matrix_product(left.rotation, right.rotation),
                "rotation of {left_index} * {right_index}"
            );
            assert_eq!(
                product_operation.action.sign(),
                left.action.sign() * right.action.sign(),
                "channel of {left_index} * {right_index}"
            );
            assert_eq!(
                product_operation.kind.is_antiunitary(),
                left.kind.is_antiunitary() ^ right.kind.is_antiunitary(),
                "kind of {left_index} * {right_index}"
            );
            let raw = composed_translation(left, right);
            let residual: Vec3 = std::array::from_fn(|axis| {
                raw[axis] - product_operation.translation[axis] - product.cell_shift[axis] as f64
            });
            assert!(
                cart_norm(lattice, residual) <= tolerance.position + 1e-12,
                "translation residual of {left_index} * {right_index}"
            );
            assert!(product.action_residual.abs() <= 1e-12);
        }
    }
}

/// Every pinned rotation must be a unimodular integer matrix, i.e. an element
/// of `GL(3, Z)`: `det(W) = +-1`.  Integer orthogonality `W^T W = I` is *not*
/// asserted here: the table stores the rotation in the record's database basis
/// frame (`G0` or `L0`, see [`SpinSpaceGroupRecord::basis`]), and a standard
/// crystallographic setting with non-orthogonal axes (e.g. hexagonal) produces
/// unimodular but non-orthogonal integer matrices.  Only the lattice action is
/// used by the identification, so unimodularity is the basis-independent
/// invariant.
fn assert_unimodular_lattice_rotation(rotation: &Mat3I) {
    let determinant = rotation_determinant(rotation);
    assert_eq!(
        determinant.abs(),
        1,
        "rotation must be unimodular (GL(3,Z)), determinant {determinant}"
    );
}

fn assert_identifier_round_trip(record: &SpinSpaceGroupRecord) {
    let parsed = SpinSpaceGroupNumber::parse(record.identifier)
        .unwrap_or_else(|| panic!("{} must parse", record.identifier));
    assert_eq!(parsed, record.number);
    assert_eq!(record.number.to_string(), record.identifier);
    assert_eq!(parsed.g0(), record.number.g0());
    assert_eq!(parsed.l0(), record.number.l0());
    assert_eq!(parsed.k_index(), record.number.k_index());
    assert_eq!(parsed.map_index(), record.number.map_index());
    assert_eq!(parsed.letter(), record.number.letter());
    assert_eq!(record.number.letter(), 'L');
}

// ---------------------------------------------------------------------------
// 1. Database count and provenance
// ---------------------------------------------------------------------------

#[test]
fn database_count_and_provenance_are_pinned() {
    assert_eq!(collinear_ssg_record_count(), RECORD_COUNT);
    assert_eq!(collinear_spin_space_groups().len(), RECORD_COUNT);
    assert_eq!(collinear_ssg_release_id(), 2);
    assert_eq!(collinear_ssg_release_name(), "local-mysql8-output");
    assert_eq!(
        collinear_ssg_release_source_sha256(),
        "fcb74e809a9153e73e245c8b7dd63016d32da8589ec8588fb3266bcf09824505"
    );
    assert_eq!(
        collinear_ssg_manifest_sha256(),
        "4216cdee23364a51adeb46413e4c85e0eab18f0528c890498a4302c6d9d8fed7"
    );
    assert_eq!(
        collinear_ssg_identifiers_sha256(),
        "29fa98580fe4853e01d4354df905f300ed7ccba7d8e403a654754c190c99bd53"
    );
    for hash in [
        collinear_ssg_release_source_sha256(),
        collinear_ssg_manifest_sha256(),
        collinear_ssg_identifiers_sha256(),
    ] {
        assert_eq!(hash.len(), 64);
        assert!(hash.bytes().all(|byte| byte.is_ascii_hexdigit()));
        assert_eq!(hash, hash.to_ascii_lowercase());
    }

    // One type-I, one type-II and one C4-altermagnet representative.
    let pm3m_fm = record("221.221.1.1.L");
    assert_eq!(pm3m_fm.group_type, "typeI");
    assert_eq!(pm3m_fm.basis, "G0");
    assert_eq!(pm3m_fm.spin_point_group, "∞m");
    assert_eq!(pm3m_fm.international_symbol, "P1m1-31m∞m1");
    assert_eq!(pm3m_fm.num_all_elements, 48);
    assert_eq!(pm3m_fm.operations.len(), 48);

    let ortho = record("47.47.2.1.L");
    assert_eq!(ortho.group_type, "typeII");
    assert_eq!(ortho.basis, "L0");
    assert_eq!(ortho.international_symbol, "P_a1m1m1m∞m1");
    assert_eq!(ortho.operations.len(), 16);

    let smallest_type_ii = record("1.1.2.1.L");
    assert_eq!(smallest_type_ii.operations.len(), 2);
    assert_eq!(smallest_type_ii.operations[1].translation, [0.5, 0.0, 0.0]);
    assert_eq!(
        smallest_type_ii.operations[1].channel,
        CollinearAction::Swap
    );

    // The typed lookup, the identifier lookup and the iterator agree.
    for number in [
        record("221.221.1.1.L").number,
        record("221.200.1.1.L").number,
        record("229.229.1.1.L").number,
    ] {
        let by_number = collinear_spin_space_group(number).unwrap();
        let by_identifier = collinear_spin_space_group_by_identifier(by_number.identifier).unwrap();
        assert!(std::ptr::eq(by_number, by_identifier));
        assert_eq!(by_number.number, number);
    }
    assert!(collinear_spin_space_group_by_identifier("230.230.8.48.Z").is_none());

    // Canonical parsing is strict (the S1 review's L3 finding).
    for text in [
        "+1.1.1.1.L",
        "01.1.1.1.L",
        "1.1.1.1.l",
        "1.1.1.1.LL",
        "1.1.1.1",
    ] {
        assert!(
            SpinSpaceGroupNumber::parse(text).is_none(),
            "accepted non-canonical {text:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Pinned table structure and independently counted aggregates
// ---------------------------------------------------------------------------

#[test]
fn pinned_table_is_structurally_complete() {
    let mut identifiers = BTreeSet::new();
    let mut numbers = BTreeSet::new();
    let mut signatures: BTreeSet<Vec<(Mat3I, [i64; 3], i8)>> = BTreeSet::new();
    let mut previous_number: Option<SpinSpaceGroupNumber> = None;
    let mut operations = 0usize;
    let mut preserve = 0usize;
    let mut swap = 0usize;
    let mut records_with_swap = 0usize;
    let mut type_i = 0usize;
    let mut type_ii = 0usize;
    let mut basis_g0 = 0usize;
    let mut basis_l0 = 0usize;

    for record in collinear_spin_space_groups() {
        assert!(identifiers.insert(record.identifier));
        assert!(numbers.insert(record.number));
        assert_identifier_round_trip(record);
        assert_eq!(record.configuration, SpinGroupKind::Collinear);
        assert!(matches!(record.group_type, "typeI" | "typeII"));
        assert!(matches!(record.basis, "G0" | "L0"));
        assert_eq!(
            record.num_all_elements as usize,
            record.operations.len(),
            "{}",
            record.identifier
        );
        if let Some(previous) = previous_number {
            assert!(previous < record.number, "records must be strictly sorted");
        }
        previous_number = Some(record.number);
        match record.group_type {
            "typeI" => type_i += 1,
            "typeII" => type_ii += 1,
            other => panic!("unexpected group type {other}"),
        }
        match record.basis {
            "G0" => basis_g0 += 1,
            "L0" => basis_l0 += 1,
            other => panic!("unexpected basis {other}"),
        }

        let mut record_signature = Vec::with_capacity(record.operations.len());
        let mut record_has_swap = false;
        let mut seen = BTreeSet::new();
        for operation in record.operations {
            operations += 1;
            assert_unimodular_lattice_rotation(&operation.rotation);
            for component in operation.translation {
                assert!(component.is_finite());
                assert!((0.0..1.0).contains(&component));
            }
            match operation.channel {
                CollinearAction::Preserve => preserve += 1,
                CollinearAction::Swap => {
                    swap += 1;
                    record_has_swap = true;
                }
            }
            let key = (
                operation.rotation,
                [
                    (operation.translation[0] * 1e9).round() as i64,
                    (operation.translation[1] * 1e9).round() as i64,
                    (operation.translation[2] * 1e9).round() as i64,
                ],
                i64::from(operation.channel.sign()) as i8,
            );
            assert!(
                seen.insert(key),
                "{}: duplicate operation",
                record.identifier
            );
            record_signature.push(key);
        }
        if record_has_swap {
            records_with_swap += 1;
        }
        // The identity must be present and preserve the channel.
        let identity = record
            .operations
            .iter()
            .find(|operation| {
                operation.rotation == IDENTITY_ROTATION && operation.translation == [0.0, 0.0, 0.0]
            })
            .expect("identity operation");
        assert_eq!(identity.channel, CollinearAction::Preserve);

        record_signature.sort_unstable();
        assert!(
            signatures.insert(record_signature),
            "{}: duplicate quotient signature in the pinned table",
            record.identifier
        );
    }

    assert_eq!(identifiers.len(), RECORD_COUNT);
    assert_eq!(numbers.len(), RECORD_COUNT);
    assert_eq!(signatures.len(), RECORD_COUNT);
    assert_eq!(operations, QUOTIENT_OPERATION_COUNT);
    assert_eq!(preserve, PRESERVE_COUNT);
    assert_eq!(swap, SWAP_COUNT);
    assert_eq!(records_with_swap, RECORDS_WITH_SWAP);
    assert_eq!((type_i, type_ii), (904, 517));
    assert_eq!((basis_g0, basis_l0), (904, 517));
}

// ---------------------------------------------------------------------------
// 3. Every pinned quotient signature identifies itself
// ---------------------------------------------------------------------------

#[test]
fn every_pinned_signature_identifies_itself() {
    let mut count = 0usize;
    for pinned in collinear_spin_space_groups() {
        let identification = identify_collinear_signature(pinned.operations)
            .unwrap_or_else(|error| panic!("{}: {error}", pinned.identifier));
        assert_eq!(identification.number, pinned.number);
        assert_eq!(identification.identifier(), pinned.identifier);
        assert_eq!(
            identification.international_symbol(),
            pinned.international_symbol
        );
        assert!(std::ptr::eq(identification.record, pinned));

        // The explicit-tolerance entry point must agree with the default one.
        let with_tolerance = identify_collinear_signature_with_tolerance(
            pinned.operations,
            DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE,
        )
        .unwrap_or_else(|error| panic!("{}: {error}", pinned.identifier));
        assert_eq!(with_tolerance.number, pinned.number);
        count += 1;
    }
    assert_eq!(count, RECORD_COUNT);
}

// ---------------------------------------------------------------------------
// 4-6. Detection -> identification for the three stage fixtures
// ---------------------------------------------------------------------------

#[test]
fn cubic_ferromagnet_detects_and_identifies_full_group() {
    let crystal = cubic_ferromagnet([0., 0., 1.]);
    let group = detect(&crystal, CollinearMoments::Vectors(&[[0., 0., 1.]]));
    check_full_group_certificate(&crystal, &group);

    let pinned = record("221.221.1.1.L");
    assert_eq!(group.len(), 2 * pinned.operations.len());
    assert_eq!(count_kind(&group, true), pinned.operations.len());
    assert_eq!(count_kind(&group, false), pinned.operations.len());
    signature_matches_record(&group, pinned);
    // The ferromagnet has only Preserve channels in both halves.
    assert!(
        group
            .operations()
            .iter()
            .all(|operation| operation.action == CollinearAction::Preserve)
    );
    // A0 = [C2_perp T || E] is the anti-unitary identity-rotation partner.
    let a0 = find_operation(
        &group,
        IDENTITY_ROTATION,
        [0.; 3],
        CollinearAction::Preserve,
        OperationKind::Antiunitary,
    )
    .expect("A0 present");
    let a0_square = group.product(a0, a0).unwrap();
    assert_eq!(a0_square.operation_index, group.identity_index());
    assert_eq!(
        group.operations()[a0_square.operation_index].kind,
        OperationKind::Unitary
    );

    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "221.221.1.1.L");
    assert_eq!(identification.number.g0(), 221);
    assert_eq!(identification.number.l0(), 221);
    assert_eq!(identification.number.k_index(), 1);
    assert_eq!(identification.number.map_index(), 1);
    assert!(std::ptr::eq(identification.record, pinned));

    // A global spin rotation and signed-scalar input are equivalent: only the
    // quotient signature matters, no oriented axis is assumed.
    for direction in [[0., 0., 1.], [1., 1., 1.], [1., 2., 3.]] {
        let rotated = cubic_ferromagnet(direction);
        let rotated_group = detect(&rotated, CollinearMoments::Vectors(&[direction]));
        let id = identify_collinear_spin_space_group(&rotated_group).unwrap();
        assert_eq!(id.identifier(), "221.221.1.1.L");
    }
    let signed = cubic_ferromagnet([0., 0., 1.]);
    let signed_group = detect(&signed, CollinearMoments::Signed(&[2.]));
    assert_eq!(signed_group.len(), group.len());
    assert_eq!(
        identify_collinear_spin_space_group(&signed_group)
            .unwrap()
            .identifier(),
        "221.221.1.1.L"
    );
}

#[test]
fn orthorhombic_afm_detects_and_identifies_full_group() {
    let crystal = orthorhombic_afm();
    let group = detect(&crystal, CollinearMoments::Signed(&[2., -2.]));
    check_full_group_certificate(&crystal, &group);

    let pinned = record("47.47.2.1.L");
    assert_eq!(group.len(), 2 * pinned.operations.len());
    assert_eq!(group.len(), 32);
    assert_eq!(count_kind(&group, true), 16);
    assert_eq!(count_kind(&group, false), 16);
    signature_matches_record(&group, pinned);

    // The half-lattice spin translation is a unitary Swap; its square is the
    // identity with one integer cell shift, and its anti-unitary partner has
    // the same spatial part.
    let half_translation = find_operation(
        &group,
        IDENTITY_ROTATION,
        [0.5, 0., 0.],
        CollinearAction::Swap,
        OperationKind::Unitary,
    )
    .expect("unitary half translation");
    let square = group.product(half_translation, half_translation).unwrap();
    assert_eq!(square.operation_index, group.identity_index());
    assert_eq!(square.cell_shift, [1, 0, 0]);
    let anti_half = find_operation(
        &group,
        IDENTITY_ROTATION,
        [0.5, 0., 0.],
        CollinearAction::Swap,
        OperationKind::Antiunitary,
    )
    .expect("anti-unitary half translation");
    let anti_square = group.product(anti_half, anti_half).unwrap();
    assert_eq!(anti_square.operation_index, group.identity_index());
    assert_eq!(anti_square.cell_shift, [1, 0, 0]);
    assert_eq!(
        group.operations()[anti_square.operation_index].kind,
        OperationKind::Unitary
    );

    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "47.47.2.1.L");
    assert!(std::ptr::eq(identification.record, pinned));
}

#[test]
fn c4_altermagnet_detects_and_identifies_full_group() {
    let crystal = c4_altermagnet();
    let group = detect(&crystal, CollinearMoments::Signed(&[2., -2., 0.]));
    check_full_group_certificate(&crystal, &group);

    let pinned = record("123.47.1.1.L");
    assert_eq!(group.len(), 32);
    assert_eq!(group.len(), 2 * pinned.operations.len());
    assert_eq!(count_kind(&group, true), 16);
    assert_eq!(count_kind(&group, false), 16);
    signature_matches_record(&group, pinned);

    // The C4 rotation is a unitary Swap and its square is the preserving C2.
    let c4 = find_operation(
        &group,
        C4_Z,
        [0.; 3],
        CollinearAction::Swap,
        OperationKind::Unitary,
    )
    .expect("unitary C4 Swap");
    let c4_square = group.product(c4, c4).unwrap();
    assert_eq!(group.operations()[c4_square.operation_index].rotation, C2_Z);
    assert_eq!(
        group.operations()[c4_square.operation_index].action,
        CollinearAction::Preserve
    );
    assert!(
        find_operation(
            &group,
            C4_Z,
            [0.; 3],
            CollinearAction::Swap,
            OperationKind::Antiunitary
        )
        .is_some()
    );
    // The nonmagnetic marker atom excludes the identity-rotation Swap.
    assert!(
        find_operation(
            &group,
            IDENTITY_ROTATION,
            [0.; 3],
            CollinearAction::Swap,
            OperationKind::Unitary
        )
        .is_none()
    );
    // Inversion acting on the magnetic pair: Preserve.
    assert!(
        find_operation(
            &group,
            [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            [0.; 3],
            CollinearAction::Preserve,
            OperationKind::Unitary
        )
        .is_some()
    );

    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "123.47.1.1.L");
    assert!(std::ptr::eq(identification.record, pinned));
}

// ---------------------------------------------------------------------------
// 7. Identification tolerance contract and negative controls
// ---------------------------------------------------------------------------

#[test]
fn identification_tolerance_and_negative_controls() {
    assert_eq!(DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE, 1e-8);
    let pinned = record("221.221.1.1.L");

    // 5e-8 fractional noise exceeds the default, so the exact entry point must
    // refuse; the explicit 1e-7 tolerance accepts it.
    let noise = 5e-8;
    let noisy: Vec<SpinSpaceGroupOperation> = pinned
        .operations
        .iter()
        .map(|operation| SpinSpaceGroupOperation {
            rotation: operation.rotation,
            translation: operation.translation.map(|value| value + noise),
            channel: operation.channel,
        })
        .collect();
    assert!(matches!(
        identify_collinear_signature(&noisy),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));
    let relaxed = identify_collinear_signature_with_tolerance(&noisy, 1e-7)
        .expect("relaxed tolerance must absorb the noise");
    assert_eq!(relaxed.identifier(), "221.221.1.1.L");

    // Any physically meaningful explicit tolerance preserves every pinned
    // signature.  Distinct operations of one record with the same rotation are
    // separated by at least 1/3 in fractional translation in the pinned release
    // (trigonal records), so tolerances well below 1/3 cannot merge operations;
    // this is the acceptance boundary of the explicit-tolerance entry point.
    for tolerance in [1e-6, 1e-3, 1e-2] {
        for pinned_record in collinear_spin_space_groups() {
            let identified =
                identify_collinear_signature_with_tolerance(pinned_record.operations, tolerance)
                    .unwrap_or_else(|error| {
                        panic!(
                            "{} at tolerance {tolerance}: {error}",
                            pinned_record.identifier
                        )
                    });
            assert_eq!(identified.number, pinned_record.number);
        }
    }

    // Invalid tolerances are typed errors on both entry points.
    for tolerance in [0.0, -1e-8, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.5] {
        assert!(matches!(
            identify_collinear_signature_with_tolerance(pinned.operations, tolerance),
            Err(SpinSpaceGroupIdentificationError::InvalidTolerance { .. })
        ));
    }
    let crystal = cubic_ferromagnet([0., 0., 1.]);
    let group = detect(&crystal, CollinearMoments::Signed(&[2.]));
    for tolerance in [0.0, f64::NAN, 0.5] {
        assert!(matches!(
            identify_collinear_spin_space_group_with_tolerance(&group, tolerance),
            Err(SpinSpaceGroupIdentificationError::InvalidTolerance { .. })
        ));
    }
    assert_eq!(
        identify_collinear_spin_space_group_with_tolerance(
            &group,
            DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE
        )
        .unwrap()
        .identifier(),
        identify_collinear_spin_space_group(&group)
            .unwrap()
            .identifier()
    );

    // Empty signature.
    assert!(matches!(
        identify_collinear_signature(&[]),
        Err(SpinSpaceGroupIdentificationError::EmptyGroup)
    ));

    // Non-finite translation.
    let mut non_finite = pinned.operations.to_vec();
    non_finite[0].translation[1] = f64::NAN;
    assert!(matches!(
        identify_collinear_signature(&non_finite),
        Err(SpinSpaceGroupIdentificationError::NonFiniteTranslation { .. })
    ));

    // Conflicting channels for the same spatial operation (exact duplicate and
    // near duplicate within the matching tolerance).
    let identity = SpinSpaceGroupOperation {
        rotation: IDENTITY_ROTATION,
        translation: [0.; 3],
        channel: CollinearAction::Preserve,
    };
    for offset in [0.0, 5e-9] {
        let mut with_swap = vec![identity];
        with_swap.push(SpinSpaceGroupOperation {
            rotation: IDENTITY_ROTATION,
            translation: [offset, 0., 0.],
            channel: CollinearAction::Swap,
        });
        assert!(matches!(
            identify_collinear_signature(&with_swap),
            Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
        ));
    }

    // Corrupted channel: one operation flipped is not in the pinned release.
    let mut flipped = pinned.operations.to_vec();
    flipped[1].channel = match flipped[1].channel {
        CollinearAction::Preserve => CollinearAction::Swap,
        CollinearAction::Swap => CollinearAction::Preserve,
    };
    assert!(matches!(
        identify_collinear_signature(&flipped),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));

    // Removed and added operations must fail closed as well.
    let mut removed = pinned.operations.to_vec();
    removed.pop();
    assert!(matches!(
        identify_collinear_signature(&removed),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));
    let mut added = pinned.operations.to_vec();
    added.push(SpinSpaceGroupOperation {
        rotation: IDENTITY_ROTATION,
        translation: [0.25, 0.25, 0.25],
        channel: CollinearAction::Preserve,
    });
    assert!(matches!(
        identify_collinear_signature(&added),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));

    // A genuinely non-standard standard-cell origin is refused, not guessed.
    let shifted_origin = Crystal::new(CUBIC, vec![[0.1, 0.1, 0.1]], vec![26])
        .unwrap()
        .with_magnetic(vec![[0., 0., 1.]])
        .unwrap();
    let shifted_group = detect(&shifted_origin, CollinearMoments::Vectors(&[[0., 0., 1.]]));
    assert!(matches!(
        identify_collinear_spin_space_group(&shifted_group),
        Err(SpinSpaceGroupIdentificationError::NoExactMatch)
    ));

    // The three fixtures are pairwise distinguished (the identification number
    // is not a constant): this pins the discriminating power of the table.
    let ortho_group = detect(&orthorhombic_afm(), CollinearMoments::Signed(&[2., -2.]));
    let c4_group = detect(&c4_altermagnet(), CollinearMoments::Signed(&[2., -2., 0.]));
    let ids = [
        identify_collinear_spin_space_group(&group)
            .unwrap()
            .identifier(),
        identify_collinear_spin_space_group(&ortho_group)
            .unwrap()
            .identifier(),
        identify_collinear_spin_space_group(&c4_group)
            .unwrap()
            .identifier(),
    ];
    assert_eq!(ids, ["221.221.1.1.L", "47.47.2.1.L", "123.47.1.1.L"]);
}

// ---------------------------------------------------------------------------
// 8. docs/spin-group.md identification example
// ---------------------------------------------------------------------------

#[test]
fn docs_spin_group_identification_example_values() {
    // This is the exact construction of the identification example in
    // `docs/spin-group.md`.
    let crystal = Crystal::new(
        [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]],
        vec![[0., 0., 0.], [0.5, 0., 0.]],
        vec![26, 26],
    )
    .unwrap()
    .with_magnetic(vec![[0., 0., 2.], [0., 0., -2.]])
    .unwrap();
    let group = detect_collinear_spin_group(
        &crystal,
        CollinearMoments::Signed(&[2., -2.]),
        SpinTolerances::default(),
    )
    .unwrap();
    assert_eq!(group.kind(), SpinGroupKind::Collinear);
    let identification = identify_collinear_spin_space_group(&group).unwrap();
    assert_eq!(identification.identifier(), "47.47.2.1.L");
    assert_eq!(identification.international_symbol(), "P_a1m1m1m∞m1");

    // The documented tolerance and field/semantics statements.
    let tolerances = SpinTolerances::default();
    assert_eq!(
        (tolerances.position, tolerances.moment, tolerances.direction),
        (1e-5, 1e-5, 1e-5)
    );
    for operation in group.operations() {
        assert!(matches!(
            operation.action,
            CollinearAction::Preserve | CollinearAction::Swap
        ));
        assert!(matches!(
            operation.kind,
            OperationKind::Unitary | OperationKind::Antiunitary
        ));
        assert_ne!(operation.rotation, [[0; 3]; 3]);
    }
    // The record exposes exactly the fields the documentation lists.
    let pinned = identification.record;
    assert_eq!(pinned.identifier, "47.47.2.1.L");
    assert_eq!(pinned.configuration, SpinGroupKind::Collinear);
    assert_eq!(pinned.group_type, "typeII");
    assert_eq!(pinned.basis, "L0");
    assert_eq!(pinned.spin_point_group, "∞/mm");
    assert_eq!((pinned.tspin, pinned.gs), (-1, -1));
    assert_eq!(pinned.gso, "∞m");
    assert_eq!(pinned.num_all_elements as usize, pinned.operations.len());
    assert!(!pinned.international_notation_latex.is_empty());
}

// ---------------------------------------------------------------------------
// 9. Cross-cutting: full-group structure, validation and kind independence
// ---------------------------------------------------------------------------

#[test]
fn cross_cutting_full_group_validation_and_kind_independence() {
    let fixtures: [(&str, Crystal, CollinearMoments<'_>); 3] = [
        (
            "221.221.1.1.L",
            cubic_ferromagnet([0., 0., 1.]),
            CollinearMoments::Signed(&[2.]),
        ),
        (
            "47.47.2.1.L",
            orthorhombic_afm(),
            CollinearMoments::Signed(&[2., -2.]),
        ),
        (
            "123.47.1.1.L",
            c4_altermagnet(),
            CollinearMoments::Signed(&[2., -2., 0.]),
        ),
    ];

    for (identifier, crystal, moments) in fixtures {
        let group = detect(&crystal, moments);
        check_full_group_certificate(&crystal, &group);
        let pinned = record(identifier);
        signature_matches_record(&group, pinned);

        // Complete full target: the public validator accepts the group that
        // detection built and returns the same two-half certificate.
        let target: Vec<SpinOperation<CollinearAction>> = group.operations().to_vec();
        let validated =
            validate_collinear_spin_group(&crystal, moments, &target, SpinTolerances::default())
                .unwrap_or_else(|error| panic!("{identifier}: full target must validate: {error}"));
        assert_eq!(validated.len(), group.len());
        assert_eq!(count_kind(&validated, true), pinned.operations.len());
        signature_matches_record(&validated, pinned);
        assert_eq!(
            identify_collinear_spin_space_group(&validated)
                .unwrap()
                .identifier(),
            identifier
        );

        // Anti-unitarity is not part of the identification signature: feeding
        // only the anti-unitary half (same spatial data and channel) must
        // identify the same record.
        let anti_only: Vec<SpinSpaceGroupOperation> = group
            .operations()
            .iter()
            .filter(|operation| operation.kind.is_antiunitary())
            .map(|operation| SpinSpaceGroupOperation {
                rotation: operation.rotation,
                translation: operation.translation,
                channel: operation.action,
            })
            .collect();
        assert_eq!(
            identify_collinear_signature(&anti_only)
                .unwrap_or_else(|error| panic!("{identifier}: anti half: {error}"))
                .identifier(),
            identifier
        );

        // A missing anti-unitary half is a typed refusal, never a silently
        // smaller group.
        let mut truncated = target.clone();
        truncated.pop();
        assert!(matches!(
            validate_collinear_spin_group(&crystal, moments, &truncated, SpinTolerances::default()),
            Err(SpinGroupError::KindCount { .. })
        ));

        // A corrupted anti-unitary companion is refused as a missing partner.
        let mut corrupted = target.clone();
        let last = corrupted.len() - 1;
        corrupted[last].translation[0] += 0.1;
        assert!(matches!(
            validate_collinear_spin_group(&crystal, moments, &corrupted, SpinTolerances::default()),
            Err(SpinGroupError::AntiunitaryPartner { .. })
        ));

        // A corrupted unitary operation is refused as well (any typed error).
        let mut corrupted_unitary = target;
        corrupted_unitary[0].translation[0] += 0.1;
        assert!(
            validate_collinear_spin_group(
                &crystal,
                moments,
                &corrupted_unitary,
                SpinTolerances::default()
            )
            .is_err()
        );
    }
}
