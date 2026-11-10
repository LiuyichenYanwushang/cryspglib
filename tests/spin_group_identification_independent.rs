//! Independent review tests for the S2 standard-cell collinear
//! spin-space-group identification layer.
//!
//! This file is self-contained: it uses only the public `cryspglib` API, it
//! performs no I/O and it does not depend on `target/ssg-collinear-records-v1`
//! at runtime.  Every expected value below was derived by re-reading the
//! pinned FindSpinGroup exports (`<identifier>.json`) and converting the
//! symbolic spin matrices into the quotient channel, independently of the
//! implementation under test.  The raw exports use the spin axis `z`; the
//! spin-only matrices `I`, `R_z(phi)` and `m_100` act as `+1` (Preserve) and
//! `-I` acts as `-1` (Swap); identical `(W, t)` rows are collapsed.
//!
//! The derived fixtures are:
//!
//! * `221.221.1.1.L` (typeI, G0, 48 operations): the 48 signed permutation
//!   matrices of `m-3m`, all with `t = 0` and `Preserve` (cubic ferromagnet).
//! * `47.47.2.1.L` (typeII, L0, 16 operations): the eight `mmm` diagonal
//!   matrices, each with `t = 0`/`Preserve` and `t = (1/2, 0, 0)`/`Swap`
//!   (orthorhombic two-atom antiferromagnet in the doubled L0 cell).
//! * `1.1.2.1.L` (typeII, L0, 2 operations): `{I|0}`/`Preserve` and
//!   `{I|1/2 0 0}`/`Swap`, the smallest spin-translation record.
//! * `221.225.2.1.L` (typeII, L0, 384 operations): the 48 cubic rotations
//!   combined with all eight `(0/1/2)^3` spin translations; the channel is
//!   `Preserve` for an even number of half components and `Swap` otherwise.
//! * `225.225.1.1.L` / `229.229.1.1.L` (typeI, G0, 192 / 96 operations): the
//!   F- and I-centered cubic ferromagnets, i.e. the 48 cubic rotations with
//!   the 4 F-centering / 2 I-centering translations, all `Preserve`.

use std::collections::BTreeSet;
use std::ptr;

use cryspglib::Crystal;
use cryspglib::spin_group::database::{
    SpinSpaceGroupNumber, SpinSpaceGroupOperation, SpinSpaceGroupRecord,
    collinear_spin_space_group_by_identifier, collinear_spin_space_groups,
    collinear_ssg_record_count,
};
use cryspglib::spin_group::{
    CollinearAction, CollinearMoments, CollinearSpinGroup, CollinearSsgIdentification,
    SpinGroupKind, SpinSpaceGroupIdentificationError, SpinTolerances, detect_collinear_spin_group,
    identify_collinear_signature, identify_collinear_spin_space_group,
};

const CUBIC: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const ORTHORHOMBIC: [[f64; 3]; 3] = [[2., 0., 0.], [0., 3., 0.], [0., 0., 4.]];
/// A `2 x 1 x 1` supercell description of the simple-cubic ferromagnet.
const TETRAGONAL_SUPERCELL: [[f64; 3]; 3] = [[2., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

const IDENTITY_ROTATION: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const ZERO_TRANSLATION: [f64; 3] = [0., 0., 0.];
const HALF_A: [f64; 3] = [0.5, 0., 0.];
const HALF_BODY: [f64; 3] = [0.5, 0.5, 0.5];

const TYPE_II_384_RECORDS: [&str; 4] = [
    "221.225.2.1.L",
    "221.226.2.1.L",
    "224.227.2.1.L",
    "224.228.2.1.L",
];

/// A signature operation reduced to an exactly comparable key.  The quantized
/// translation uses the same `1e-9` fractional grid as the implementation's
/// internal index only to make *constructed fixtures* and pinned rows testable
/// for exact equality at the `1/12` rational resolution of the exports.
type CanonicalOperation = ([i32; 9], [i64; 3], i8);

/// A `Z2` character of the cubic point group: one of the four homomorphisms
/// `m-3m -> {+1, -1}`.
type ChannelCharacter = fn([[i32; 3]; 3]) -> i32;

fn operation(
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
    channel: CollinearAction,
) -> SpinSpaceGroupOperation {
    SpinSpaceGroupOperation {
        rotation,
        translation,
        channel,
    }
}

fn canonical_operation(operation: SpinSpaceGroupOperation) -> CanonicalOperation {
    let mut rotation = [0i32; 9];
    for (row, values) in operation.rotation.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            rotation[3 * row + column] = *value;
        }
    }
    let translation = operation
        .translation
        .map(|value| (value.rem_euclid(1.0) * 1e9).round() as i64);
    (rotation, translation, operation.channel.sign())
}

fn canonical(operations: &[SpinSpaceGroupOperation]) -> BTreeSet<CanonicalOperation> {
    operations
        .iter()
        .copied()
        .map(canonical_operation)
        .collect()
}

fn record(identifier: &str) -> &'static SpinSpaceGroupRecord {
    collinear_spin_space_group_by_identifier(identifier)
        .unwrap_or_else(|| panic!("pinned record {identifier} is missing"))
}

fn detect(crystal: &Crystal, moments: CollinearMoments<'_>) -> CollinearSpinGroup {
    detect_collinear_spin_group(crystal, moments, SpinTolerances::default())
        .unwrap_or_else(|error| panic!("detection failed: {error}"))
}

fn identify(
    signature: &[SpinSpaceGroupOperation],
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    identify_collinear_signature(signature)
}

/// The quotient part of a detected full group: the unitary and anti-unitary
/// partners share `(W, t, channel)`, so collapsing them must reproduce the
/// pinned quotient signature.
fn quotient_signature(group: &CollinearSpinGroup) -> Vec<SpinSpaceGroupOperation> {
    group
        .operations()
        .iter()
        .map(|operation| SpinSpaceGroupOperation {
            rotation: operation.rotation,
            translation: operation.translation,
            channel: operation.action,
        })
        .collect()
}

fn assert_identifies_to(
    result: Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError>,
    identifier: &str,
) {
    let identification = result.unwrap_or_else(|error| panic!("{identifier}: {error}"));
    let expected = record(identifier);
    assert_eq!(identification.identifier(), identifier);
    assert_eq!(identification.number, expected.number);
    assert_eq!(identification.record.number, expected.number);
    assert_eq!(
        identification.international_symbol(),
        expected.international_symbol
    );
    assert!(
        ptr::eq(identification.record, expected),
        "{identifier} did not return the canonical pinned record pointer"
    );
}

fn assert_no_exact_match(
    result: Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError>,
) {
    match result {
        Err(SpinSpaceGroupIdentificationError::NoExactMatch) => {}
        other => panic!("expected NoExactMatch, got {other:?}"),
    }
}

fn cubic_rotations() -> Vec<[[i32; 3]; 3]> {
    let permutations = [
        [0usize, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut rotations = Vec::with_capacity(48);
    for permutation in permutations {
        for signs in 0u32..8 {
            let mut rotation = [[0i32; 3]; 3];
            for (row, &column) in permutation.iter().enumerate() {
                rotation[row][column] = if signs & (1 << row) == 0 { 1 } else { -1 };
            }
            rotations.push(rotation);
        }
    }
    rotations
}

fn diagonal_rotations() -> Vec<[[i32; 3]; 3]> {
    let mut rotations = Vec::with_capacity(8);
    for sx in [1i32, -1] {
        for sy in [1i32, -1] {
            for sz in [1i32, -1] {
                rotations.push([[sx, 0, 0], [0, sy, 0], [0, 0, sz]]);
            }
        }
    }
    rotations
}

fn cubic_ferromagnet_signature() -> Vec<SpinSpaceGroupOperation> {
    cubic_rotations()
        .into_iter()
        .map(|rotation| operation(rotation, ZERO_TRANSLATION, CollinearAction::Preserve))
        .collect()
}

fn orthorhombic_afm_signature() -> Vec<SpinSpaceGroupOperation> {
    let mut operations = Vec::with_capacity(16);
    for rotation in diagonal_rotations() {
        operations.push(operation(
            rotation,
            ZERO_TRANSLATION,
            CollinearAction::Preserve,
        ));
        operations.push(operation(rotation, HALF_A, CollinearAction::Swap));
    }
    operations
}

fn centered_cubic_ferromagnet_signature(centering: &[[f64; 3]]) -> Vec<SpinSpaceGroupOperation> {
    let mut operations = Vec::with_capacity(48 * centering.len());
    for rotation in cubic_rotations() {
        for translation in centering {
            operations.push(operation(rotation, *translation, CollinearAction::Preserve));
        }
    }
    operations
}

/// `221.225.2.1.L`: F-lattice spin translations on top of the cubic rotations.
fn type_ii_spin_translation_signature() -> Vec<SpinSpaceGroupOperation> {
    let mut operations = Vec::with_capacity(384);
    for rotation in cubic_rotations() {
        for b0 in 0..2 {
            for b1 in 0..2 {
                for b2 in 0..2 {
                    let translation = [
                        0.5 * f64::from(b0),
                        0.5 * f64::from(b1),
                        0.5 * f64::from(b2),
                    ];
                    let channel = if (b0 + b1 + b2) % 2 == 0 {
                        CollinearAction::Preserve
                    } else {
                        CollinearAction::Swap
                    };
                    operations.push(operation(rotation, translation, channel));
                }
            }
        }
    }
    operations
}

/// `1.1.2.1.L`: one spin translation `{I | 1/2 0 0}` with a `Swap` channel.
fn smallest_spin_translation_signature() -> Vec<SpinSpaceGroupOperation> {
    vec![
        operation(
            IDENTITY_ROTATION,
            ZERO_TRANSLATION,
            CollinearAction::Preserve,
        ),
        operation(IDENTITY_ROTATION, HALF_A, CollinearAction::Swap),
    ]
}

fn shift_all(
    operations: &[SpinSpaceGroupOperation],
    shift: [f64; 3],
) -> Vec<SpinSpaceGroupOperation> {
    operations
        .iter()
        .map(|input| {
            operation(
                input.rotation,
                [
                    input.translation[0] + shift[0],
                    input.translation[1] + shift[1],
                    input.translation[2] + shift[2],
                ],
                input.channel,
            )
        })
        .collect()
}

fn flip_channels_of_identity(
    operations: &[SpinSpaceGroupOperation],
) -> Vec<SpinSpaceGroupOperation> {
    let mut flipped = operations.to_vec();
    for operation in &mut flipped {
        if operation.rotation == IDENTITY_ROTATION && operation.translation == ZERO_TRANSLATION {
            operation.channel = CollinearAction::Swap;
        }
    }
    flipped
}

fn cubic_half_shifted_origin_signature() -> Vec<SpinSpaceGroupOperation> {
    // Operations of the single-atom cubic crystal with its atom at
    // `r = (1/2, 0, 0)`: `{W | r - W r}`.  This is a pure origin choice, so
    // the identified number must not change by a silent guess; the exact
    // standard-cell layer must refuse it.
    let origin = HALF_A;
    cubic_rotations()
        .into_iter()
        .map(|rotation| {
            let mut translation = [0.0f64; 3];
            for (row, values) in rotation.iter().enumerate() {
                for (column, value) in values.iter().enumerate() {
                    translation[row] += f64::from(*value) * origin[column];
                }
            }
            let translation = [
                origin[0] - translation[0],
                origin[1] - translation[1],
                origin[2] - translation[2],
            ];
            operation(rotation, translation, CollinearAction::Preserve)
        })
        .collect()
}

struct ExpectedRecord {
    identifier: &'static str,
    number: [u16; 4],
    group_type: &'static str,
    basis: &'static str,
    num_all_elements: u32,
    international_symbol: &'static str,
    spin_point_group: &'static str,
    tspin: i8,
    gs: i8,
}

fn signed_permutation_parts(rotation: [[i32; 3]; 3]) -> ([usize; 3], i32) {
    let mut permutation = [0usize; 3];
    let mut sign = 1i32;
    for (row, values) in rotation.iter().enumerate() {
        let column = values
            .iter()
            .position(|value| *value != 0)
            .expect("a signed permutation matrix has one non-zero entry per row");
        permutation[row] = column;
        sign *= values[column];
    }
    (permutation, sign)
}

fn signed_permutation_parity(rotation: [[i32; 3]; 3]) -> i32 {
    let (permutation, _) = signed_permutation_parts(rotation);
    let mut parity = 1i32;
    for (index, &left) in permutation.iter().enumerate() {
        for &right in &permutation[index + 1..] {
            if left > right {
                parity = -parity;
            }
        }
    }
    parity
}

fn signed_permutation_determinant(rotation: [[i32; 3]; 3]) -> i32 {
    let (_, diagonal_sign) = signed_permutation_parts(rotation);
    diagonal_sign * signed_permutation_parity(rotation)
}

fn signed_permutation_determinant_times_parity(rotation: [[i32; 3]; 3]) -> i32 {
    signed_permutation_determinant(rotation) * signed_permutation_parity(rotation)
}

/// Positions of one generic m-3m orbit, in the same order as `cubic_rotations`.
fn cubic_generic_orbit(seed: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[[i32; 3]; 3]>) {
    let rotations = cubic_rotations();
    let mut positions = Vec::with_capacity(rotations.len());
    for rotation in &rotations {
        let mut position = [0.0f64; 3];
        for (row, values) in rotation.iter().enumerate() {
            position[row] = values
                .iter()
                .enumerate()
                .map(|(column, value)| f64::from(*value) * seed[column])
                .sum::<f64>()
                .rem_euclid(1.0);
        }
        positions.push(position);
    }
    (positions, rotations)
}

fn assert_record_metadata(expected: &ExpectedRecord) {
    let identifier = expected.identifier;
    let pinned = record(identifier);
    assert_eq!(pinned.number.g0(), expected.number[0], "{identifier} g0");
    assert_eq!(pinned.number.l0(), expected.number[1], "{identifier} l0");
    assert_eq!(
        pinned.number.k_index(),
        expected.number[2],
        "{identifier} k_index (raw ik)"
    );
    assert_eq!(
        pinned.number.map_index(),
        expected.number[3],
        "{identifier} map_index (raw Index)"
    );
    assert_eq!(pinned.number.letter(), 'L', "{identifier} letter");
    assert_eq!(pinned.number.to_string(), identifier);
    assert_eq!(SpinSpaceGroupNumber::parse(identifier), Some(pinned.number));
    assert_eq!(
        pinned.group_type, expected.group_type,
        "{identifier} group_type"
    );
    assert_eq!(pinned.basis, expected.basis, "{identifier} basis");
    assert_eq!(
        pinned.num_all_elements, expected.num_all_elements,
        "{identifier} num_all_elements"
    );
    assert_eq!(
        pinned.operations.len(),
        expected.num_all_elements as usize,
        "{identifier} operation count"
    );
    assert_eq!(
        pinned.international_symbol, expected.international_symbol,
        "{identifier} symbol"
    );
    assert_eq!(
        pinned.spin_point_group, expected.spin_point_group,
        "{identifier} spin point group"
    );
    assert_eq!(pinned.tspin, expected.tspin, "{identifier} tspin");
    assert_eq!(pinned.gs, expected.gs, "{identifier} gs");
}

#[test]
fn pinned_database_metadata_matches_the_reviewed_exports() {
    assert_eq!(collinear_ssg_record_count(), 1421);
    let records: Vec<&'static SpinSpaceGroupRecord> = collinear_spin_space_groups().collect();
    assert_eq!(records.len(), 1421);

    let mut identifiers = BTreeSet::new();
    let mut operation_count = 0usize;
    let mut records_with_swap = 0usize;
    let mut type_ii_384 = Vec::new();
    for pinned in &records {
        assert_eq!(pinned.configuration, SpinGroupKind::Collinear);
        assert_eq!(pinned.number.to_string(), pinned.identifier);
        assert_eq!(
            SpinSpaceGroupNumber::parse(pinned.identifier),
            Some(pinned.number),
            "{}",
            pinned.identifier
        );
        assert_eq!(
            pinned.num_all_elements as usize,
            pinned.operations.len(),
            "{}",
            pinned.identifier
        );
        identifiers.insert(pinned.identifier);
        operation_count += pinned.operations.len();
        let identity_operations = pinned
            .operations
            .iter()
            .filter(|operation| {
                operation.rotation == IDENTITY_ROTATION && operation.translation == ZERO_TRANSLATION
            })
            .count();
        assert_eq!(identity_operations, 1, "{}", pinned.identifier);
        let preserving_identities = pinned
            .operations
            .iter()
            .filter(|operation| {
                operation.rotation == IDENTITY_ROTATION
                    && operation.translation == ZERO_TRANSLATION
                    && operation.channel == CollinearAction::Preserve
            })
            .count();
        assert_eq!(preserving_identities, 1, "{}", pinned.identifier);
        if pinned
            .operations
            .iter()
            .any(|operation| operation.channel == CollinearAction::Swap)
        {
            records_with_swap += 1;
        }
        if pinned.operations.len() == 384 {
            type_ii_384.push(pinned.identifier);
        }
    }
    assert_eq!(identifiers.len(), 1421);
    // Aggregates of the pinned release, repeated from the raw exports.
    assert_eq!(operation_count, 29_457);
    assert_eq!(records_with_swap, 1_191);
    assert_eq!(type_ii_384, TYPE_II_384_RECORDS.to_vec());

    // Metadata of the fixtures, derived from the reviewed raw exports.
    assert_record_metadata(&ExpectedRecord {
        identifier: "221.221.1.1.L",
        number: [221, 221, 1, 1],
        group_type: "typeI",
        basis: "G0",
        num_all_elements: 48,
        international_symbol: "P1m1-31m∞m1",
        spin_point_group: "∞m",
        tspin: 1,
        gs: 1,
    });
    assert_record_metadata(&ExpectedRecord {
        identifier: "47.47.2.1.L",
        number: [47, 47, 2, 1],
        group_type: "typeII",
        basis: "L0",
        num_all_elements: 16,
        international_symbol: "P_a1m1m1m∞m1",
        spin_point_group: "∞/mm",
        tspin: -1,
        gs: -1,
    });
    assert_record_metadata(&ExpectedRecord {
        identifier: "1.1.2.1.L",
        number: [1, 1, 2, 1],
        group_type: "typeII",
        basis: "L0",
        num_all_elements: 2,
        international_symbol: "P_S11∞m1",
        spin_point_group: "∞/mm",
        tspin: -1,
        gs: -1,
    });
    assert_record_metadata(&ExpectedRecord {
        identifier: "225.225.1.1.L",
        number: [225, 225, 1, 1],
        group_type: "typeI",
        basis: "G0",
        num_all_elements: 192,
        international_symbol: "F1m1-31m∞m1",
        spin_point_group: "∞m",
        tspin: 1,
        gs: 1,
    });
    assert_record_metadata(&ExpectedRecord {
        identifier: "229.229.1.1.L",
        number: [229, 229, 1, 1],
        group_type: "typeI",
        basis: "G0",
        num_all_elements: 96,
        international_symbol: "I1m1-31m∞m1",
        spin_point_group: "∞m",
        tspin: 1,
        gs: 1,
    });
    for identifier in TYPE_II_384_RECORDS {
        let pinned = record(identifier);
        assert_eq!(pinned.group_type, "typeII");
        assert_eq!(pinned.basis, "L0");
        assert_eq!(pinned.num_all_elements, 384);
        assert_eq!(pinned.operations.len(), 384);
        assert_eq!(pinned.spin_point_group, "∞/mm");
        assert_eq!(pinned.tspin, -1);
        assert_eq!(pinned.gs, -1);
        assert_eq!(
            pinned
                .operations
                .iter()
                .filter(|operation| operation.channel == CollinearAction::Preserve)
                .count(),
            192,
            "{identifier} Preserve count"
        );
        assert_eq!(
            pinned
                .operations
                .iter()
                .filter(|operation| operation.channel == CollinearAction::Swap)
                .count(),
            192,
            "{identifier} Swap count"
        );
    }
}

#[test]
fn every_pinned_quotient_signature_resolves_to_its_own_record() {
    let mut count = 0usize;
    for pinned in collinear_spin_space_groups() {
        let identification = identify(pinned.operations)
            .unwrap_or_else(|error| panic!("{}: {error}", pinned.identifier));
        assert_eq!(
            identification.number, pinned.number,
            "{}",
            pinned.identifier
        );
        assert_eq!(
            identification.record.identifier, pinned.identifier,
            "{}",
            pinned.identifier
        );
        assert!(
            ptr::eq(identification.record, pinned),
            "{} did not resolve to the canonical pinned record",
            pinned.identifier
        );
        assert_eq!(identification.identifier(), pinned.identifier);
        assert_eq!(
            identification.international_symbol(),
            pinned.international_symbol
        );
        count += 1;
    }
    assert_eq!(count, 1421);
}

#[test]
fn pinned_signatures_stay_stable_under_subtolerance_noise_and_identity_channel_flips() {
    let noise = [5e-9, 0., 0.];
    let mut noisy_matches = 0usize;
    let mut flipped_identities = 0usize;
    for pinned in collinear_spin_space_groups() {
        // Sub-tolerance fractional noise in every translation must still
        // resolve to the same record (this exercises the tolerant fallback).
        let noisy = shift_all(pinned.operations, noise);
        let identification = identify(&noisy)
            .unwrap_or_else(|error| panic!("{} with noise: {error}", pinned.identifier));
        assert!(
            ptr::eq(identification.record, pinned),
            "{} with sub-tolerance noise",
            pinned.identifier
        );
        noisy_matches += 1;

        // No pinned record contains a channel-flipped identity operation
        // (checked in `pinned_database_metadata_matches_the_reviewed_exports`),
        // so flipping the identity channel of any record must be refused.
        let flipped = flip_channels_of_identity(pinned.operations);
        assert!(
            matches!(
                identify(&flipped),
                Err(SpinSpaceGroupIdentificationError::NoExactMatch)
            ),
            "{} with a channel-flipped identity",
            pinned.identifier
        );
        flipped_identities += 1;
    }
    assert_eq!(noisy_matches, 1421);
    assert_eq!(flipped_identities, 1421);
}

#[test]
fn type_ii_spin_translation_records_are_covered_by_the_signature_api() {
    let signature = type_ii_spin_translation_signature();
    assert_eq!(signature.len(), 384);
    assert_identifies_to(identify(&signature), "221.225.2.1.L");
    assert_eq!(
        canonical(&signature),
        canonical(record("221.225.2.1.L").operations),
        "the reconstructed 384-operation spin-translation signature must match the pinned row"
    );

    // All four pinned 384-operation type-II L0 records resolve through the
    // public signature entry point to exactly their own number.
    for identifier in TYPE_II_384_RECORDS {
        assert_identifies_to(identify(record(identifier).operations), identifier);
    }

    // Channel flip on a 384-operation signature is not silently absorbed.
    let mut flipped = signature.clone();
    let first = &mut flipped[0];
    first.channel = match first.channel {
        CollinearAction::Preserve => CollinearAction::Swap,
        CollinearAction::Swap => CollinearAction::Preserve,
    };
    assert_no_exact_match(identify(&flipped));

    // One translation displaced by a quarter of a cell is a different group.
    let mut displaced = signature;
    displaced[1].translation = [0.25, 0.25, 0.25];
    assert_no_exact_match(identify(&displaced));
}

#[test]
fn cubic_ferromagnet_identifies_the_pm3m_type_i_record() {
    let crystal = Crystal::new(CUBIC, vec![ZERO_TRANSLATION], vec![26]).unwrap();
    let expected = record("221.221.1.1.L");
    let constructed = cubic_ferromagnet_signature();
    assert_eq!(constructed.len(), 48);

    // Moment along z, along [111] and along a non-crystallographic direction:
    // no oriented magnetic axis is used, so all three quotient groups are the
    // same.  This is a global spin rotation of the same ferromagnet.
    for direction in [[0., 0., 1.], [1., 1., 1.], [1., 2., 3.]] {
        let moments = [direction];
        let group = detect(&crystal, CollinearMoments::Vectors(&moments));
        assert_eq!(group.len(), 96, "full unitary+anti-unitary group");
        assert_eq!(
            canonical(&quotient_signature(&group)),
            canonical(&constructed),
            "detected quotient signature for moment {direction:?}"
        );
        assert_identifies_to(identify_collinear_spin_space_group(&group), "221.221.1.1.L");
        assert_eq!(
            group.len(),
            2 * expected.operations.len(),
            "the anti-unitary half must collapse onto the same quotient"
        );
    }

    // Caller-signed scalar input describes the same ferromagnet.
    let signed = detect(&crystal, CollinearMoments::Signed(&[2.0]));
    assert_identifies_to(
        identify_collinear_spin_space_group(&signed),
        "221.221.1.1.L",
    );
    assert_eq!(
        canonical(&quotient_signature(&signed)),
        canonical(&constructed)
    );

    // The quotient signature alone resolves to the same record.
    assert_identifies_to(identify(&constructed), "221.221.1.1.L");
}

#[test]
fn orthorhombic_two_atom_afm_identifies_the_pinned_type_ii_record() {
    // Reordered atoms relative to the implementation's own test fixture: atom 0
    // is the moment at (1/2, 0, 0), atom 1 the moment at the origin.
    let crystal = Crystal::new(ORTHORHOMBIC, vec![HALF_A, ZERO_TRANSLATION], vec![26, 26]).unwrap();
    let expected = record("47.47.2.1.L");
    let constructed = orthorhombic_afm_signature();
    assert_eq!(constructed.len(), 16);

    // The signed-scalar order parameter is defined only up to a global sign
    // (the choice of the common axis orientation); both inputs are the same
    // physical antiferromagnet.
    for moments in [[-2., 2.], [2., -2.]] {
        let group = detect(&crystal, CollinearMoments::Signed(&moments));
        assert_eq!(group.len(), 32);
        assert_eq!(
            canonical(&quotient_signature(&group)),
            canonical(&constructed),
            "detected quotient signature for signed moments {moments:?}"
        );
        assert_identifies_to(identify_collinear_spin_space_group(&group), "47.47.2.1.L");
        assert_eq!(
            group.len(),
            2 * expected.operations.len(),
            "the anti-unitary half must collapse onto the same quotient"
        );
    }

    // A global spin rotation of the vector moments (here z -> x, which is not a
    // crystal symmetry direction for this orthorhombic cell) must not change
    // the identified number either.
    let rotated_moments = [[-2., 0., 0.], [2., 0., 0.]];
    let rotated = detect(&crystal, CollinearMoments::Vectors(&rotated_moments));
    assert_eq!(
        canonical(&quotient_signature(&rotated)),
        canonical(&constructed)
    );
    assert_identifies_to(identify_collinear_spin_space_group(&rotated), "47.47.2.1.L");

    assert_identifies_to(identify(&constructed), "47.47.2.1.L");
}

#[test]
fn pm3m_coupling_records_are_distinguished_by_the_channel_character() {
    // One generic `m-3m` orbit of 48 atoms with a sign pattern
    // `sigma_{W r0} = nu(W) * |m|`.  `nu` must be a homomorphism `m-3m -> Z2`;
    // the four possibilities are trivial, det, permutation parity and their
    // product, and each reconstructs the channel pattern of one of the four
    // 48-operation `Pm-3m` records (derived from the raw exports).
    let (positions, rotations) = cubic_generic_orbit([0.13, 0.21, 0.29]);
    let crystal = Crystal::new(CUBIC, positions, vec![26; rotations.len()]).unwrap();
    let characters: [(&str, ChannelCharacter); 4] = [
        ("221.221.1.1.L", |_| 1),
        ("221.207.1.1.L", signed_permutation_determinant),
        ("221.200.1.1.L", signed_permutation_parity),
        ("221.215.1.1.L", signed_permutation_determinant_times_parity),
    ];
    for (identifier, character) in characters {
        let moments: Vec<[f64; 3]> = rotations
            .iter()
            .map(|rotation| [0., 0., 2.0 * f64::from(character(*rotation))])
            .collect();
        let expected_signature: Vec<SpinSpaceGroupOperation> = rotations
            .iter()
            .map(|rotation| {
                let channel = if character(*rotation) > 0 {
                    CollinearAction::Preserve
                } else {
                    CollinearAction::Swap
                };
                operation(*rotation, ZERO_TRANSLATION, channel)
            })
            .collect();
        let group = detect(&crystal, CollinearMoments::Vectors(&moments));
        assert_eq!(group.len(), 96, "{identifier}");
        assert_eq!(
            canonical(&quotient_signature(&group)),
            canonical(&expected_signature),
            "{identifier} quotient signature"
        );
        assert_identifies_to(identify_collinear_spin_space_group(&group), identifier);
    }
}

#[test]
fn centered_cubic_ferromagnets_keep_the_centering_translations() {
    let i_centering = [ZERO_TRANSLATION, HALF_BODY];
    let i_crystal = Crystal::new(CUBIC, vec![ZERO_TRANSLATION, HALF_BODY], vec![26, 26]).unwrap();
    let i_moments = [[0., 0., 1.], [0., 0., 1.]];
    let i_group = detect(&i_crystal, CollinearMoments::Vectors(&i_moments));
    assert_eq!(i_group.len(), 192);
    assert_eq!(
        canonical(&quotient_signature(&i_group)),
        canonical(&centered_cubic_ferromagnet_signature(&i_centering))
    );
    assert_identifies_to(
        identify_collinear_spin_space_group(&i_group),
        "229.229.1.1.L",
    );

    let f_centering = [
        ZERO_TRANSLATION,
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.],
    ];
    let f_crystal = Crystal::new(CUBIC, f_centering.to_vec(), vec![26, 26, 26, 26]).unwrap();
    let f_moments = [[0., 0., 1.]; 4];
    let f_group = detect(&f_crystal, CollinearMoments::Vectors(&f_moments));
    assert_eq!(f_group.len(), 384);
    assert_eq!(
        canonical(&quotient_signature(&f_group)),
        canonical(&centered_cubic_ferromagnet_signature(&f_centering))
    );
    assert_identifies_to(
        identify_collinear_spin_space_group(&f_group),
        "225.225.1.1.L",
    );
}

#[test]
fn nonstandard_cells_are_refused_not_guessed() {
    // A translated origin of the cubic ferromagnet changes the fractional
    // translations of the detected operations; the exact standard-cell layer
    // must refuse instead of normalizing the origin silently.
    let shifted = Crystal::new(CUBIC, vec![[0.25, 0.125, 0.0625]], vec![26]).unwrap();
    let shifted_group = detect(&shifted, CollinearMoments::Vectors(&[[0., 0., 1.]]));
    assert_eq!(shifted_group.len(), 96);
    assert_no_exact_match(identify_collinear_spin_space_group(&shifted_group));

    // A half-lattice origin shift is still a nonstandard origin choice.
    let half_shifted = Crystal::new(CUBIC, vec![HALF_A], vec![26]).unwrap();
    let half_group = detect(&half_shifted, CollinearMoments::Vectors(&[[0., 0., 1.]]));
    assert_eq!(half_group.len(), 96);
    assert_no_exact_match(identify_collinear_spin_space_group(&half_group));

    // The same pure origin choice on the orthorhombic antiferromagnet.
    let shift = [0.1, 0.2, 0.3];
    let shifted_positions = vec![[HALF_A[0] + shift[0], shift[1], shift[2]], shift];
    let shifted_afm = Crystal::new(ORTHORHOMBIC, shifted_positions, vec![26, 26]).unwrap();
    let shifted_afm_group = detect(&shifted_afm, CollinearMoments::Signed(&[-2., 2.]));
    assert_eq!(shifted_afm_group.len(), 32);
    assert_no_exact_match(identify_collinear_spin_space_group(&shifted_afm_group));

    // A genuine 2x1x1 supercell description of the simple-cubic ferromagnet
    // (two basis atoms) is not the database standard cell either.
    let supercell = Crystal::new(
        TETRAGONAL_SUPERCELL,
        vec![ZERO_TRANSLATION, HALF_A],
        vec![26, 26],
    )
    .unwrap();
    let supercell_group = detect(
        &supercell,
        CollinearMoments::Vectors(&[[0., 0., 1.], [0., 0., 1.]]),
    );
    assert_eq!(supercell_group.len(), 64);
    assert_no_exact_match(identify_collinear_spin_space_group(&supercell_group));
}

#[test]
fn signature_negative_controls_are_refused() {
    // Empty signature.
    assert!(matches!(
        identify(&[]),
        Err(SpinSpaceGroupIdentificationError::EmptyGroup)
    ));

    // Conflicting channel for the same spatial operation.
    let mut conflicting = smallest_spin_translation_signature();
    conflicting.push(operation(
        IDENTITY_ROTATION,
        ZERO_TRANSLATION,
        CollinearAction::Swap,
    ));
    assert!(matches!(
        identify(&conflicting),
        Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
    ));
    let mut conflicting_translation = smallest_spin_translation_signature();
    conflicting_translation.push(operation(
        IDENTITY_ROTATION,
        HALF_A,
        CollinearAction::Preserve,
    ));
    assert!(matches!(
        identify(&conflicting_translation),
        Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
    ));

    // A duplicated operation with the identical channel is set-equivalent and
    // must still resolve to its own record.
    let mut duplicated = smallest_spin_translation_signature();
    duplicated.push(operation(
        IDENTITY_ROTATION,
        ZERO_TRANSLATION,
        CollinearAction::Preserve,
    ));
    assert_identifies_to(identify(&duplicated), "1.1.2.1.L");

    // The detected full group repeats every quotient operation once as its
    // anti-unitary companion with the same (W, t, channel).  That kind split is
    // not part of the signature, so duplicating every operation with the same
    // channel must still identify the same record.
    for (signature, identifier) in [
        (cubic_ferromagnet_signature(), "221.221.1.1.L"),
        (orthorhombic_afm_signature(), "47.47.2.1.L"),
        (type_ii_spin_translation_signature(), "221.225.2.1.L"),
    ] {
        let mut doubled = signature.clone();
        doubled.extend(signature);
        assert_identifies_to(identify(&doubled), identifier);
    }

    // Signature order is not part of the quotient data: reversing the pinned
    // operations must identify the same record.
    let mut reversed = record("47.47.2.1.L").operations.to_vec();
    reversed.reverse();
    assert_identifies_to(identify(&reversed), "47.47.2.1.L");
    let mut rotated = cubic_ferromagnet_signature();
    rotated.rotate_left(17);
    assert_identifies_to(identify(&rotated), "221.221.1.1.L");

    // A single flipped channel (the identity one, then the spin-translation
    // one) cannot match any pinned record: every pinned record has exactly one
    // identity and it is Preserve.
    assert_no_exact_match(identify(&flip_channels_of_identity(
        &cubic_ferromagnet_signature(),
    )));
    assert_no_exact_match(identify(&flip_channels_of_identity(
        &orthorhombic_afm_signature(),
    )));
    let mut flipped_spin_translation = smallest_spin_translation_signature();
    flipped_spin_translation[1].channel = CollinearAction::Preserve;
    assert_no_exact_match(identify(&flipped_spin_translation));

    // Removing, adding or changing an operation changes the group.
    let mut removed = cubic_ferromagnet_signature();
    removed.pop();
    assert_no_exact_match(identify(&removed));

    let mut added = cubic_ferromagnet_signature();
    added.push(operation(
        IDENTITY_ROTATION,
        [0.25, 0., 0.],
        CollinearAction::Preserve,
    ));
    assert_no_exact_match(identify(&added));

    let mut changed = cubic_ferromagnet_signature();
    let target = changed
        .iter_mut()
        .find(|operation| operation.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        .expect("the inversion operation is part of the cubic group");
    target.rotation = [[2, 0, 0], [0, 1, 0], [0, 0, 1]];
    assert_no_exact_match(identify(&changed));

    // Large translation offsets, applied to all operations or to one.
    assert_no_exact_match(identify(&shift_all(
        &cubic_ferromagnet_signature(),
        [0.1, 0., 0.],
    )));
    let mut one_offset = cubic_ferromagnet_signature();
    one_offset[0].translation = [0.1, 0.1, 0.1];
    assert_no_exact_match(identify(&one_offset));

    // A pure half-lattice origin shift of the ferromagnet.
    assert_no_exact_match(identify(&cubic_half_shifted_origin_signature()));
}

#[test]
fn origin_choices_that_preserve_the_operation_data_still_identify() {
    // `(1/2, 1/2, 1/2)` is an origin shift with `s - W s` integral for every
    // cubic `W`; the detected quotient data is therefore unchanged and the
    // correct number is returned.  This is a positive control for the exact
    // signature comparison (it must not "guess" a different record).
    let body_centered_origin = Crystal::new(CUBIC, vec![[0.5, 0.5, 0.5]], vec![26]).unwrap();
    let group = detect(
        &body_centered_origin,
        CollinearMoments::Vectors(&[[0., 0., 1.]]),
    );
    assert_eq!(
        canonical(&quotient_signature(&group)),
        canonical(&cubic_ferromagnet_signature())
    );
    assert_identifies_to(identify_collinear_spin_space_group(&group), "221.221.1.1.L");

    // The corresponding half-cell shift of the L0 frame of the orthorhombic
    // antiferromagnet (`s = (1/2, 0, 0)`) also leaves `s - W s` integral.
    let crystal = Crystal::new(
        ORTHORHOMBIC,
        vec![[1.0, 0., 0.], [0.5, 0., 0.]],
        vec![26, 26],
    )
    .unwrap();
    let group = detect(&crystal, CollinearMoments::Signed(&[-2., 2.]));
    assert_identifies_to(identify_collinear_spin_space_group(&group), "47.47.2.1.L");
}

#[test]
fn translation_tolerance_and_lattice_normalization() {
    let base = smallest_spin_translation_signature();

    // Fractions differing by a lattice vector are the same quotient operation.
    assert_identifies_to(identify(&shift_all(&base, [1.0, 1.0, 1.0])), "1.1.2.1.L");
    assert_identifies_to(identify(&shift_all(&base, [-1.0, -1.0, -1.0])), "1.1.2.1.L");

    // Clearly sub-tolerance noise is absorbed by the fallback comparison.
    assert_identifies_to(identify(&shift_all(&base, [5e-9, 0., 0.])), "1.1.2.1.L");
    assert_identifies_to(identify(&shift_all(&base, [-5e-9, 0., 0.])), "1.1.2.1.L");

    // The non-zero translation exercises the tolerant fallback path directly.
    let mut half_noise = base.clone();
    half_noise[1].translation = [HALF_A[0] + 9e-9, 0., 0.];
    assert_identifies_to(identify(&half_noise), "1.1.2.1.L");

    // Values beyond the frozen `1e-8` absolute fractional tolerance are
    // refused, for a near-zero translation and for the non-zero one.
    let mut near_zero = base.clone();
    near_zero[0].translation = [1.1e-8, 0., 0.];
    assert_no_exact_match(identify(&near_zero));

    let mut half_too_far = base.clone();
    half_too_far[1].translation = [HALF_A[0] + 1.1e-8, 0., 0.];
    assert_no_exact_match(identify(&half_too_far));
    assert_no_exact_match(identify(&shift_all(&base, [2e-8, 0., 0.])));

    // Zero-translation record: the documented tolerance squash must accept
    // sub-tolerance noise and reject noise above the tolerance.
    assert_identifies_to(
        identify(&shift_all(&cubic_ferromagnet_signature(), [9e-9, 0., 0.])),
        "221.221.1.1.L",
    );
    assert_no_exact_match(identify(&shift_all(
        &cubic_ferromagnet_signature(),
        [1.1e-8, 0., 0.],
    )));

    // A large translation offset is not inside the tolerance.
    assert_no_exact_match(identify(&shift_all(&base, [0.1, 0., 0.])));
    let mut one_large = base;
    one_large[1].translation = [HALF_A[0] + 0.1, 0., 0.];
    assert_no_exact_match(identify(&one_large));
}
