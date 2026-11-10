//! Acceptance tests for the pinned collinear spin-space-group database.
//!
//! The generated table is checked structurally and against a small set of
//! records whose metadata is also present in the public FindSpinGroup export.
//! No test reads `target/`: the generated table must be self-contained.

use std::collections::{BTreeSet, HashMap};

use cryspglib::spin_group::SpinGroupKind;
use cryspglib::spin_group::collinear::CollinearAction;
use cryspglib::spin_group::database::{
    SpinSpaceGroupNumber, SpinSpaceGroupOperation, collinear_spin_space_group,
    collinear_spin_space_group_by_identifier, collinear_spin_space_groups,
    collinear_ssg_identifiers_sha256, collinear_ssg_manifest_sha256, collinear_ssg_record_count,
    collinear_ssg_release_id, collinear_ssg_release_name, collinear_ssg_release_source_sha256,
};

const EXPECTED_COUNT: usize = 1_421;
const EXPECTED_RELEASE_ID: u32 = 2;
const EXPECTED_RELEASE_NAME: &str = "local-mysql8-output";
const EXPECTED_SOURCE_SHA256: &str =
    "fcb74e809a9153e73e245c8b7dd63016d32da8589ec8588fb3266bcf09824505";
const EXPECTED_MANIFEST_SHA256: &str =
    "4216cdee23364a51adeb46413e4c85e0eab18f0528c890498a4302c6d9d8fed7";
const EXPECTED_IDENTIFIERS_SHA256: &str =
    "29fa98580fe4853e01d4354df905f300ed7ccba7d8e403a654754c190c99bd53";

fn determinant(rotation: &[[i32; 3]; 3]) -> i64 {
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

fn translation_key(value: f64) -> i64 {
    let reduced = value.rem_euclid(1.0);
    let reduced = if (reduced - 1.0).abs() < 1e-12 {
        0.0
    } else {
        reduced
    };
    (reduced * 1_000_000_000.0).round() as i64
}

fn channel_sign(channel: CollinearAction) -> i32 {
    match channel {
        CollinearAction::Preserve => 1,
        CollinearAction::Swap => -1,
    }
}

type OperationKey = ([[i32; 3]; 3], [i64; 3]);

fn operation_key(operation: &SpinSpaceGroupOperation) -> OperationKey {
    (
        operation.rotation,
        operation.translation.map(translation_key),
    )
}

fn signature_key(operations: &[SpinSpaceGroupOperation]) -> Vec<(OperationKey, i32)> {
    let mut key: Vec<_> = operations
        .iter()
        .map(|operation| (operation_key(operation), channel_sign(operation.channel)))
        .collect();
    key.sort_unstable();
    key
}

fn panic_identifier(identifier: &str) -> ! {
    panic!("failed on pinned record {identifier}")
}

#[test]
fn record_count_and_provenance_are_pinned() {
    assert_eq!(collinear_ssg_record_count(), EXPECTED_COUNT);
    assert_eq!(collinear_spin_space_groups().len(), EXPECTED_COUNT);
    assert_eq!(collinear_ssg_release_id(), EXPECTED_RELEASE_ID);
    assert_eq!(collinear_ssg_release_name(), EXPECTED_RELEASE_NAME);
    assert_eq!(
        collinear_ssg_release_source_sha256(),
        EXPECTED_SOURCE_SHA256
    );
    assert_eq!(collinear_ssg_manifest_sha256(), EXPECTED_MANIFEST_SHA256);
    assert_eq!(
        collinear_ssg_identifiers_sha256(),
        EXPECTED_IDENTIFIERS_SHA256
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
}

#[test]
fn identifiers_are_unique_and_round_trip() {
    let mut identifiers = BTreeSet::new();
    let mut numbers = BTreeSet::new();
    for record in collinear_spin_space_groups() {
        assert!(
            identifiers.insert(record.identifier),
            "{}",
            record.identifier
        );
        assert!(numbers.insert(record.number), "{}", record.identifier);
        let parsed = SpinSpaceGroupNumber::parse(record.identifier)
            .unwrap_or_else(|| panic_identifier(record.identifier));
        assert_eq!(parsed, record.number, "{}", record.identifier);
        assert_eq!(parsed.to_string(), record.identifier);
        assert_eq!(record.number.to_string(), record.identifier);
        assert_eq!(record.number.g0(), parsed.g0());
        assert_eq!(record.number.l0(), parsed.l0());
        assert_eq!(record.number.it(), parsed.it());
        assert_eq!(record.number.ik(), parsed.ik());
        assert_eq!(record.number.letter(), parsed.letter());
    }
    assert_eq!(identifiers.len(), EXPECTED_COUNT);
    assert_eq!(numbers.len(), EXPECTED_COUNT);
}

#[test]
fn number_parsing_rejects_malformed_input() {
    for text in [
        "",
        "1",
        "1.1.1.1",
        "1.1.1.1.L.5",
        "0.1.1.1.L",
        "231.1.1.1.L",
        "1.0.1.1.L",
        "1.231.1.1.L",
        "1.1.0.1.L",
        "1.1.49.1.L",
        "1.1.1.0.L",
        "1.1.1.9.L",
        "1.1.1.1.LL",
        "1.1.1.1.1",
        "-1.1.1.1.L",
        "a.b.c.d.L",
    ] {
        assert!(
            SpinSpaceGroupNumber::parse(text).is_none(),
            "accepted {text:?}"
        );
    }
    let number = SpinSpaceGroupNumber::new(221, 200, 2, 1, 'L').unwrap();
    assert_eq!(number.to_string(), "221.200.2.1.L");
    assert!(SpinSpaceGroupNumber::new(0, 200, 2, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(221, 200, 2, 1, '1').is_none());
}

#[test]
fn metadata_is_in_range_and_configuration_is_collinear() {
    let mut previous = None;
    for record in collinear_spin_space_groups() {
        let number = record.number;
        assert!((1..=230).contains(&number.g0()), "{}", record.identifier);
        assert!((1..=230).contains(&number.l0()), "{}", record.identifier);
        assert!((1..=48).contains(&number.it()), "{}", record.identifier);
        assert!((1..=8).contains(&number.ik()), "{}", record.identifier);
        assert!(number.letter().is_ascii_alphabetic());
        assert_eq!(record.configuration, SpinGroupKind::Collinear);
        assert!(
            record.group_type.starts_with("type"),
            "{}",
            record.identifier
        );
        assert!(
            ["G0", "L0"].contains(&record.basis),
            "{}",
            record.identifier
        );
        assert!(!record.spin_point_group.is_empty(), "{}", record.identifier);
        assert!(!record.gso.is_empty(), "{}", record.identifier);
        assert!(matches!(record.tspin, -1 | 1), "{}", record.identifier);
        assert!(matches!(record.gs, -1 | 1), "{}", record.identifier);
        assert!(
            record.num_all_elements as usize >= record.operations.len(),
            "{}",
            record.identifier
        );
        assert!(!record.operations.is_empty(), "{}", record.identifier);
        assert!(!record.international_symbol.is_empty());
        assert!(
            !record.international_symbol.contains(['\\', '$', '{', '}']),
            "{} has non-plain symbol {:?}",
            record.identifier,
            record.international_symbol
        );
        assert!(!record.international_notation_latex.is_empty());
        if let Some(previous) = previous {
            assert!(
                previous < number,
                "records are not sorted: {previous} then {number}"
            );
        }
        previous = Some(number);
    }
}

#[test]
fn operations_are_well_formed_and_unique() {
    for record in collinear_spin_space_groups() {
        let mut seen = BTreeSet::new();
        for operation in record.operations {
            let determinant = determinant(&operation.rotation);
            assert_eq!(
                determinant.abs(),
                1,
                "{}: rotation has determinant {determinant}",
                record.identifier
            );
            for value in operation.translation {
                assert!(value.is_finite(), "{}", record.identifier);
                assert!(
                    (-1e-12..1.0 + 1e-12).contains(&value),
                    "{}: translation {value} out of [0, 1)",
                    record.identifier
                );
            }
            assert!(
                seen.insert(operation_key(operation)),
                "{}",
                record.identifier
            );
        }
        let identity: OperationKey = ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0]);
        let identity_operation = record
            .operations
            .iter()
            .find(|operation| operation_key(operation) == identity)
            .unwrap_or_else(|| panic_identifier(record.identifier));
        assert_eq!(
            identity_operation.channel,
            CollinearAction::Preserve,
            "{}",
            record.identifier
        );
    }
}

#[test]
fn spatial_projection_is_closed_and_channel_is_a_homomorphism() {
    for record in collinear_spin_space_groups() {
        let mut table: HashMap<OperationKey, i32> = HashMap::new();
        for operation in record.operations {
            table.insert(operation_key(operation), channel_sign(operation.channel));
        }
        for left in record.operations {
            for right in record.operations {
                let mut rotation = [[0i32; 3]; 3];
                for (row, rotation_row) in rotation.iter_mut().enumerate() {
                    for (column, cell) in rotation_row.iter_mut().enumerate() {
                        *cell = (0..3)
                            .map(|inner| {
                                i64::from(left.rotation[row][inner])
                                    * i64::from(right.rotation[inner][column])
                            })
                            .sum::<i64>() as i32;
                    }
                }
                let mut translation = [0.0f64; 3];
                for (row, component) in translation.iter_mut().enumerate() {
                    *component = left.translation[row]
                        + (0..3)
                            .map(|inner| {
                                f64::from(left.rotation[row][inner]) * right.translation[inner]
                            })
                            .sum::<f64>();
                }
                let key = (rotation, translation.map(translation_key));
                let product_channel = table.get(&key).unwrap_or_else(|| {
                    panic!(
                        "{}: spatial projection is not closed (missing product)",
                        record.identifier
                    )
                });
                assert_eq!(
                    *product_channel,
                    channel_sign(left.channel) * channel_sign(right.channel),
                    "{}: channel is not a homomorphism",
                    record.identifier
                );
            }
        }
    }
}

#[test]
fn lookup_by_number_and_identifier_agree() {
    for record in collinear_spin_space_groups() {
        let by_number = collinear_spin_space_group(record.number)
            .unwrap_or_else(|| panic_identifier(record.identifier));
        assert_eq!(by_number.identifier, record.identifier);
        let by_identifier = collinear_spin_space_group_by_identifier(record.identifier)
            .unwrap_or_else(|| panic_identifier(record.identifier));
        assert_eq!(by_identifier.identifier, record.identifier);
        assert_eq!(by_identifier.number, record.number);
    }
    let absent = SpinSpaceGroupNumber::new(230, 230, 48, 8, 'Z').unwrap();
    assert!(collinear_spin_space_group(absent).is_none());
    assert!(collinear_spin_space_group_by_identifier("230.230.48.8.Z").is_none());
    assert!(collinear_spin_space_group_by_identifier("not-an-identifier").is_none());
}

#[test]
fn quotient_signatures_are_unique() {
    let mut signatures = BTreeSet::new();
    for record in collinear_spin_space_groups() {
        assert!(
            signatures.insert(signature_key(record.operations)),
            "duplicate quotient signature for {}",
            record.identifier
        );
    }
    assert_eq!(signatures.len(), EXPECTED_COUNT);
}

#[test]
fn spot_checks_against_pinned_database_records() {
    let p1 = collinear_spin_space_group_by_identifier("1.1.1.1.L").unwrap();
    assert_eq!(p1.number.to_string(), "1.1.1.1.L");
    assert_eq!(p1.spin_point_group, "∞m");
    assert_eq!(p1.gso, "∞m");
    assert_eq!(p1.tspin, 1);
    assert_eq!(p1.gs, 1);
    assert_eq!(p1.basis, "G0");
    assert_eq!(p1.group_type, "typeI");
    assert_eq!(p1.num_all_elements, 1);
    assert_eq!(p1.international_symbol, "P11∞m1");
    assert_eq!(p1.operations.len(), 1);
    assert_eq!(p1.operations[0].rotation, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    assert_eq!(p1.operations[0].translation, [0.0, 0.0, 0.0]);
    assert_eq!(p1.operations[0].channel, CollinearAction::Preserve);

    let mirror_110 = [[0, -1, 0], [-1, 0, 0], [0, 0, 1]];
    let record_200 = collinear_spin_space_group_by_identifier("221.200.1.1.L").unwrap();
    assert_eq!(record_200.spin_point_group, "∞/mm");
    assert_eq!(record_200.gs, -1);
    assert_eq!(record_200.tspin, 1);
    assert_eq!(record_200.basis, "G0");
    assert_eq!(record_200.operations.len(), 48);
    let operation_200 = record_200
        .operations
        .iter()
        .find(|operation| operation.rotation == mirror_110 && operation.translation == [0.0; 3])
        .expect("221.200 must contain the spatial mirror m_110");
    assert_eq!(operation_200.channel, CollinearAction::Swap);

    let record_221 = collinear_spin_space_group_by_identifier("221.221.1.1.L").unwrap();
    assert_eq!(record_221.spin_point_group, "∞m");
    assert_eq!(record_221.gs, 1);
    assert_eq!(record_221.operations.len(), 48);
    let operation_221 = record_221
        .operations
        .iter()
        .find(|operation| operation.rotation == mirror_110 && operation.translation == [0.0; 3])
        .expect("221.221 must contain the spatial mirror m_110");
    assert_eq!(operation_221.channel, CollinearAction::Preserve);

    let record_225 = collinear_spin_space_group_by_identifier("221.225.2.1.L").unwrap();
    assert_eq!(record_225.group_type, "typeII");
    assert_eq!(record_225.spin_point_group, "∞/mm");
    assert_eq!(record_225.tspin, -1);
    assert_eq!(record_225.gs, -1);
    assert_eq!(record_225.basis, "L0");
    assert_eq!(record_225.num_all_elements, 384);
    assert_eq!(record_225.operations.len(), 384);
    assert!(record_225.international_symbol.contains("F_S"));
}
