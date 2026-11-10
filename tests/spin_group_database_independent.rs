//! Independent acceptance checks for the pinned collinear spin-space-group
//! database added in commit `8fadcfe`.
//!
//! This file was written by a reviewer, not by the importer.  Every expected
//! value below was produced by a separate Python parser of the raw
//! FindSpinGroup exports in `target/ssg-collinear-records-v1/` (release
//! `id = 2`, `name = local-mysql8-output`) and of the verified download
//! manifest; the importer's own helper functions were not used.  The fixtures
//! therefore pin the *raw* database content, while the library is exercised
//! only through its public API.  The test does not read `target/` at runtime,
//! touch the network, or depend on any file outside this crate.
//!
//! The stored quotient signature is one spatial operation `{W|t}` plus the net
//! action on the collinear order parameter (`Preserve`/`Swap`); the continuous
//! spin rotations about the common magnetic axis about z and the universal
//! spin-only stabilizer are the canonical kernel and are not rows of the
//! signature.  All 1,421 source records therefore reduce to 29,457 finite
//! quotient operations.
//!
//! Canonical record encoding used by the reviewer-side digest oracle (and
//! reproduced here byte for byte):
//!
//! * `u32` big-endian byte length + UTF-8 bytes, in order, for
//!   `identifier`, `international_symbol`, `international_notation_latex`,
//!   `group_type`, `basis`, `spin_point_group`, `gso`;
//! * one byte `0` for `configuration == Collinear`;
//! * `tspin` and `gs` as their signed 8-bit two's-complement representation;
//! * `u32` big-endian `num_all_elements`, then `u32` big-endian number of
//!   operations;
//! * each operation as `9 x i32` big-endian row-major rotation, `3 x i64`
//!   big-endian translation in units of `1/12` of a cell, and one channel
//!   byte (`1` = `Preserve`, `2` = `Swap`), the operations sorted
//!   lexicographically by these 61 bytes;
//! * the table digest hashes `u32` big-endian length + record bytes for every
//!   record in ascending four-index order.

use std::collections::{BTreeMap, BTreeSet};

use cryspglib::spin_group::database::{
    SpinSpaceGroupNumber, SpinSpaceGroupOperation, SpinSpaceGroupRecord,
    collinear_spin_space_group, collinear_spin_space_group_by_identifier,
    collinear_spin_space_groups, collinear_ssg_identifiers_sha256, collinear_ssg_manifest_sha256,
    collinear_ssg_record_count, collinear_ssg_release_id, collinear_ssg_release_name,
    collinear_ssg_release_source_sha256,
};
use cryspglib::spin_group::{CollinearAction, SpinGroupKind};

const EXPECTED_RECORD_COUNT: usize = 1_421;
const EXPECTED_RELEASE_ID: u32 = 2;
const EXPECTED_RELEASE_NAME: &str = "local-mysql8-output";
const EXPECTED_SOURCE_SHA256: &str =
    "fcb74e809a9153e73e245c8b7dd63016d32da8589ec8588fb3266bcf09824505";
const EXPECTED_MANIFEST_SHA256: &str =
    "4216cdee23364a51adeb46413e4c85e0eab18f0528c890498a4302c6d9d8fed7";
const EXPECTED_IDENTIFIERS_SHA256: &str =
    "29fa98580fe4853e01d4354df905f300ed7ccba7d8e403a654754c190c99bd53";

const EXPECTED_TOTAL_OPERATIONS: usize = 29_457;
const EXPECTED_PRESERVE_OPERATIONS: usize = 16_941;
const EXPECTED_SWAP_OPERATIONS: usize = 12_516;
const EXPECTED_RECORDS_WITH_SWAP: usize = 1_191;

/// All raw translations are rational multiples of `1/12` in the [0, 1) cube.
const TRANSLATION_UNITS_PER_CELL: i64 = 12;

const EXPECTED_OPERATION_COUNT_HISTOGRAM: &[(usize, usize)] = &[
    (1, 1),
    (2, 11),
    (3, 3),
    (4, 93),
    (6, 39),
    (8, 391),
    (9, 1),
    (12, 100),
    (16, 459),
    (18, 9),
    (24, 74),
    (32, 127),
    (36, 12),
    (48, 42),
    (64, 6),
    (72, 2),
    (96, 25),
    (192, 22),
    (384, 4),
];

const EXPECTED_GROUP_TYPE_COUNTS: &[(&str, usize)] = &[("typeI", 904), ("typeII", 517)];
const EXPECTED_BASIS_COUNTS: &[(&str, usize)] = &[("G0", 904), ("L0", 517)];
const EXPECTED_SPIN_POINT_GROUP_COUNTS: &[(&str, usize)] = &[("∞m", 230), ("∞/mm", 1_191)];
const EXPECTED_TSPIN_COUNTS: &[(i8, usize)] = &[(1, 904), (-1, 517)];
const EXPECTED_GS_COUNTS: &[(i8, usize)] = &[(1, 230), (-1, 1_191)];

/// SHA-256 of the canonical encoding of all 1,421 records.  Derived from the
/// raw exports by the reviewer-side parser, not by hashing the generated table.
const EXPECTED_TABLE_SHA256: &str =
    "d61f5a967b543a9c700f1e22baa3e42b5a068dab8fa14fbf8da26d3ec3f1be18";

/// One raw-export fixture.  `record_sha256` is the same canonical record
/// encoding as the table digest.
struct Fixture {
    identifier: &'static str,
    number: (u16, u16, u16, u16, char),
    international_symbol: &'static str,
    international_notation_latex: &'static str,
    group_type: &'static str,
    basis: &'static str,
    spin_point_group: &'static str,
    tspin: i8,
    gs: i8,
    num_all_elements: u32,
    preserve_ops: usize,
    swap_ops: usize,
    record_sha256: &'static str,
}

const FIXTURES: &[Fixture] = &[
    // Smallest record: P1 real-space group, one quotient operation.
    Fixture {
        identifier: "1.1.1.1.L",
        number: (1, 1, 1, 1, 'L'),
        international_symbol: "P11∞m1",
        international_notation_latex: "$P\\ce{^{1}{1}}\\ce{^{\\infty m}{1}}$",
        group_type: "typeI",
        basis: "G0",
        spin_point_group: "∞m",
        tspin: 1,
        gs: 1,
        num_all_elements: 1,
        preserve_ops: 1,
        swap_ops: 0,
        record_sha256: "11d9cb0801f462cc8ee2f103b7a71785096912c9313ffe382ba114161a1ba4f3",
    },
    // Smallest type-II record: an anti-translation (1/2, 0, 0) that swaps the
    // two spin directions, in the record's L0 standard setting.
    Fixture {
        identifier: "1.1.2.1.L",
        number: (1, 1, 2, 1, 'L'),
        international_symbol: "P_S11∞m1",
        international_notation_latex: "$\\ce{}{P_{S}}\\ce{^{1}{1}}\\ce{^{\\infty m}{1}}$",
        group_type: "typeII",
        basis: "L0",
        spin_point_group: "∞/mm",
        tspin: -1,
        gs: -1,
        num_all_elements: 2,
        preserve_ops: 1,
        swap_ops: 1,
        record_sha256: "7d17a5cf2a0079ff76e856795ea39b8300c3c62068c55456f4afce0817043cc8",
    },
    // 48-operation type-I record whose m_110 spatial mirror carries Swap.
    Fixture {
        identifier: "221.200.1.1.L",
        number: (221, 200, 1, 1, 'L'),
        international_symbol: "P1m1-3-1m∞m1",
        international_notation_latex: "$P\\ce{^{1}{m}}\\ce{^{1}{\\text{-}3}}\\ce{^{\\text{-}1}{m}}\\ce{^{\\infty m}{1}}$",
        group_type: "typeI",
        basis: "G0",
        spin_point_group: "∞/mm",
        tspin: 1,
        gs: -1,
        num_all_elements: 48,
        preserve_ops: 24,
        swap_ops: 24,
        record_sha256: "fcb7149f16c6281db6b4d730e98505802844eeaa83aeb856c840a7452f0dd141",
    },
    // 48-operation type-I record with no Swap operation at all.
    Fixture {
        identifier: "221.221.1.1.L",
        number: (221, 221, 1, 1, 'L'),
        international_symbol: "P1m1-31m∞m1",
        international_notation_latex: "$P\\ce{^{1}{m}}\\ce{^{1}{\\text{-}3}}\\ce{^{1}{m}}\\ce{^{\\infty m}{1}}$",
        group_type: "typeI",
        basis: "G0",
        spin_point_group: "∞m",
        tspin: 1,
        gs: 1,
        num_all_elements: 48,
        preserve_ops: 48,
        swap_ops: 0,
        record_sha256: "c64353590ca27fd5a9fd71f9993633371f23ed9a68e5017444f66d371076fefc",
    },
    // Largest record size (384 operations), type-II, L0 basis.
    Fixture {
        identifier: "221.225.2.1.L",
        number: (221, 225, 2, 1, 'L'),
        international_symbol: "F_S1m1-31m∞m1",
        international_notation_latex: "$\\ce{}{F_{S}}\\ce{^{1}{m}}\\ce{^{1}{\\text{-}3}}\\ce{^{1}{m}}\\ce{^{\\infty m}{1}}$",
        group_type: "typeII",
        basis: "L0",
        spin_point_group: "∞/mm",
        tspin: -1,
        gs: -1,
        num_all_elements: 384,
        preserve_ops: 192,
        swap_ops: 192,
        record_sha256: "d5b193a6d33d3643941c0d09aa07aae83c5fff7a057be476f60bb3b6090d3aef",
    },
    // Second 384-operation record with a different spatial signature.
    Fixture {
        identifier: "221.226.2.1.L",
        number: (221, 226, 2, 1, 'L'),
        international_symbol: "F_S1m1-31c∞m1",
        international_notation_latex: "$\\ce{}{F_{S}}\\ce{^{1}{m}}\\ce{^{1}{\\text{-}3}}\\ce{^{1}{c}}\\ce{^{\\infty m}{1}}$",
        group_type: "typeII",
        basis: "L0",
        spin_point_group: "∞/mm",
        tspin: -1,
        gs: -1,
        num_all_elements: 384,
        preserve_ops: 192,
        swap_ops: 192,
        record_sha256: "b76bd97f35e276acb6279fbae999b19f8be3c18ac077e535de5d1a605e8b427a",
    },
];

type OperationKey = ([[i32; 3]; 3], [i64; 3]);

fn sign(channel: CollinearAction) -> i32 {
    match channel {
        CollinearAction::Preserve => 1,
        CollinearAction::Swap => -1,
    }
}

/// Exact `1/12`-cell units of a stored translation.  The raw database values
/// are rational with denominators dividing 12, so this conversion is exact up
/// to the documented `1e-9` presentation budget and is otherwise lossless.
fn translation_units(translation: [f64; 3]) -> [i64; 3] {
    translation.map(|component| {
        assert!(
            component.is_finite() && (0.0..1.0).contains(&component),
            "translation component {component} is not in [0, 1)"
        );
        let scaled = component * TRANSLATION_UNITS_PER_CELL as f64;
        let rounded = scaled.round();
        assert!(
            (scaled - rounded).abs() < 1e-9,
            "translation component {component} is not a multiple of 1/{TRANSLATION_UNITS_PER_CELL}"
        );
        let units = rounded as i64;
        assert!(
            (0..TRANSLATION_UNITS_PER_CELL).contains(&units),
            "translation units {units} outside the cell"
        );
        units
    })
}

fn operation_key(operation: &SpinSpaceGroupOperation) -> OperationKey {
    (operation.rotation, translation_units(operation.translation))
}

fn rotation_determinant(rotation: &[[i32; 3]; 3]) -> i64 {
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

fn rotation_product(left: &[[i32; 3]; 3], right: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    let mut product = [[0i32; 3]; 3];
    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, cell) in product_row.iter_mut().enumerate() {
            let sum: i64 = (0..3)
                .map(|inner| i64::from(left[row][inner]) * i64::from(right[inner][column]))
                .sum();
            *cell = sum as i32;
        }
    }
    product
}

fn translation_product(
    left_rotation: &[[i32; 3]; 3],
    left_translation: [i64; 3],
    right_translation: [i64; 3],
) -> [i64; 3] {
    let mut result = [0i64; 3];
    for (row, value) in result.iter_mut().enumerate() {
        let mut units = left_translation[row];
        for (column, right_units) in right_translation.iter().enumerate() {
            units += i64::from(left_rotation[row][column]) * right_units;
        }
        *value = units.rem_euclid(TRANSLATION_UNITS_PER_CELL);
    }
    result
}

fn operation_bytes(operation: &SpinSpaceGroupOperation) -> [u8; 61] {
    let mut bytes = [0u8; 61];
    let mut offset = 0;
    for row in operation.rotation {
        for cell in row {
            bytes[offset..offset + 4].copy_from_slice(&cell.to_be_bytes());
            offset += 4;
        }
    }
    for units in translation_units(operation.translation) {
        bytes[offset..offset + 8].copy_from_slice(&units.to_be_bytes());
        offset += 8;
    }
    bytes[60] = match operation.channel {
        CollinearAction::Preserve => 1,
        CollinearAction::Swap => 2,
    };
    bytes
}

fn canonical_signature(record: &SpinSpaceGroupRecord) -> Vec<[u8; 61]> {
    let mut operations: Vec<[u8; 61]> = record.operations.iter().map(operation_bytes).collect();
    operations.sort_unstable();
    operations
}

fn canonical_record_bytes(record: &SpinSpaceGroupRecord) -> Vec<u8> {
    let mut bytes = Vec::new();
    let strings = [
        record.identifier,
        record.international_symbol,
        record.international_notation_latex,
        record.group_type,
        record.basis,
        record.spin_point_group,
        record.gso,
    ];
    for text in strings {
        let raw = text.as_bytes();
        bytes.extend_from_slice(
            &u32::try_from(raw.len())
                .expect("string too long")
                .to_be_bytes(),
        );
        bytes.extend_from_slice(raw);
    }
    bytes.push(match record.configuration {
        SpinGroupKind::Collinear => 0,
        SpinGroupKind::Coplanar => 1,
        SpinGroupKind::NonCoplanar => 2,
    });
    bytes.push(record.tspin.to_be_bytes()[0]);
    bytes.push(record.gs.to_be_bytes()[0]);
    bytes.extend_from_slice(&record.num_all_elements.to_be_bytes());
    let operations = canonical_signature(record);
    bytes.extend_from_slice(
        &u32::try_from(operations.len())
            .expect("operation count too large")
            .to_be_bytes(),
    );
    for operation in &operations {
        bytes.extend_from_slice(operation);
    }
    bytes
}

fn canonical_table_bytes() -> Vec<u8> {
    let mut bytes = Vec::new();
    for record in collinear_spin_space_groups() {
        let record_bytes = canonical_record_bytes(record);
        bytes.extend_from_slice(
            &u32::try_from(record_bytes.len())
                .expect("record too large")
                .to_be_bytes(),
        );
        bytes.extend_from_slice(&record_bytes);
    }
    bytes
}

const SHA256_K: [u32; 64] = [
    0x428a_2f98,
    0x7137_4491,
    0xb5c0_fbcf,
    0xe9b5_dba5,
    0x3956_c25b,
    0x59f1_11f1,
    0x923f_82a4,
    0xab1c_5ed5,
    0xd807_aa98,
    0x1283_5b01,
    0x2431_85be,
    0x550c_7dc3,
    0x72be_5d74,
    0x80de_b1fe,
    0x9bdc_06a7,
    0xc19b_f174,
    0xe49b_69c1,
    0xefbe_4786,
    0x0fc1_9dc6,
    0x240c_a1cc,
    0x2de9_2c6f,
    0x4a74_84aa,
    0x5cb0_a9dc,
    0x76f9_88da,
    0x983e_5152,
    0xa831_c66d,
    0xb003_27c8,
    0xbf59_7fc7,
    0xc6e0_0bf3,
    0xd5a7_9147,
    0x06ca_6351,
    0x1429_2967,
    0x27b7_0a85,
    0x2e1b_2138,
    0x4d2c_6dfc,
    0x5338_0d13,
    0x650a_7354,
    0x766a_0abb,
    0x81c2_c92e,
    0x9272_2c85,
    0xa2bf_e8a1,
    0xa81a_664b,
    0xc24b_8b70,
    0xc76c_51a3,
    0xd192_e819,
    0xd699_0624,
    0xf40e_3585,
    0x106a_a070,
    0x19a4_c116,
    0x1e37_6c08,
    0x2748_774c,
    0x34b0_bcb5,
    0x391c_0cb3,
    0x4ed8_aa4a,
    0x5b9c_ca4f,
    0x682e_6ff3,
    0x748f_82ee,
    0x78a5_636f,
    0x84c8_7814,
    0x8cc7_0208,
    0x90be_fffa,
    0xa450_6ceb,
    0xbef9_a3f7,
    0xc671_78f2,
];

fn sha256(data: &[u8]) -> [u8; 32] {
    let mut state = [
        0x6a09_e667u32,
        0xbb67_ae85,
        0x3c6e_f372,
        0xa54f_f53a,
        0x510e_527f,
        0x9b05_688c,
        0x1f83_d9ab,
        0x5be0_cd19,
    ];
    let mut message = data.to_vec();
    let bit_length = (data.len() as u64).wrapping_mul(8);
    message.push(0x80);
    while message.len() % 64 != 56 {
        message.push(0);
    }
    message.extend_from_slice(&bit_length.to_be_bytes());
    for chunk in message.as_chunks::<64>().0 {
        let mut schedule = [0u32; 64];
        for (word, bytes) in schedule.iter_mut().zip(chunk.as_chunks::<4>().0) {
            *word = u32::from_be_bytes(*bytes);
        }
        for index in 16..64 {
            let s0 = schedule[index - 15].rotate_right(7)
                ^ schedule[index - 15].rotate_right(18)
                ^ (schedule[index - 15] >> 3);
            let s1 = schedule[index - 2].rotate_right(17)
                ^ schedule[index - 2].rotate_right(19)
                ^ (schedule[index - 2] >> 10);
            schedule[index] = schedule[index - 16]
                .wrapping_add(s0)
                .wrapping_add(schedule[index - 7])
                .wrapping_add(s1);
        }
        let mut a = state[0];
        let mut b = state[1];
        let mut c = state[2];
        let mut d = state[3];
        let mut e = state[4];
        let mut f = state[5];
        let mut g = state[6];
        let mut h = state[7];
        for (index, word) in schedule.iter().enumerate() {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let choice = (e & f) ^ ((!e) & g);
            let t1 = h
                .wrapping_add(s1)
                .wrapping_add(choice)
                .wrapping_add(SHA256_K[index])
                .wrapping_add(*word);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let majority = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(majority);
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        state[0] = state[0].wrapping_add(a);
        state[1] = state[1].wrapping_add(b);
        state[2] = state[2].wrapping_add(c);
        state[3] = state[3].wrapping_add(d);
        state[4] = state[4].wrapping_add(e);
        state[5] = state[5].wrapping_add(f);
        state[6] = state[6].wrapping_add(g);
        state[7] = state[7].wrapping_add(h);
    }
    let mut digest = [0u8; 32];
    for (index, word) in state.iter().enumerate() {
        digest[4 * index..4 * index + 4].copy_from_slice(&word.to_be_bytes());
    }
    digest
}

fn sha256_hex(data: &[u8]) -> String {
    sha256(data)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn assert_counts<K: Ord + Clone + std::fmt::Debug>(
    actual: &BTreeMap<K, usize>,
    expected: &[(K, usize)],
) {
    let expected: BTreeMap<K, usize> = expected.iter().cloned().collect();
    assert_eq!(actual, &expected);
}

fn fixture_record(fixture: &Fixture) -> &'static SpinSpaceGroupRecord {
    let record =
        collinear_spin_space_group_by_identifier(fixture.identifier).unwrap_or_else(|| {
            panic!(
                "fixture record {} missing from the table",
                fixture.identifier
            )
        });
    let (g0, l0, third, fourth, letter) = fixture.number;
    let number = SpinSpaceGroupNumber::new(g0, l0, third, fourth, letter)
        .expect("fixture number must be inside the supported ranges");
    assert_eq!(record.number, number);
    assert_eq!(record.identifier, fixture.identifier);
    let by_number = collinear_spin_space_group(number)
        .expect("fixture number must resolve through the typed lookup");
    assert_eq!(by_number.identifier, fixture.identifier);
    record
}

fn find_operation(
    record: &SpinSpaceGroupRecord,
    rotation: [[i32; 3]; 3],
    units: [i64; 3],
    channel: CollinearAction,
) -> Option<&SpinSpaceGroupOperation> {
    record.operations.iter().find(|operation| {
        operation.rotation == rotation
            && translation_units(operation.translation) == units
            && operation.channel == channel
    })
}

#[test]
fn sha256_implementation_matches_published_vectors() {
    // Without a correct hash implementation the digest checks below would be
    // meaningless, so the helper is validated against NIST test vectors.
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
}

#[test]
fn provenance_constants_are_pinned_by_the_raw_release() {
    assert_eq!(collinear_ssg_record_count(), EXPECTED_RECORD_COUNT);
    assert_eq!(collinear_spin_space_groups().len(), EXPECTED_RECORD_COUNT);
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
    for digest in [
        collinear_ssg_release_source_sha256(),
        collinear_ssg_manifest_sha256(),
        collinear_ssg_identifiers_sha256(),
    ] {
        assert_eq!(digest.len(), 64);
        assert!(
            digest
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
            "digest is not lowercase hexadecimal: {digest}"
        );
    }
}

#[test]
fn table_shape_order_and_identifier_round_trips() {
    let mut identifiers = BTreeSet::new();
    let mut numbers = BTreeSet::new();
    let mut group_types: BTreeMap<&str, usize> = BTreeMap::new();
    let mut bases: BTreeMap<&str, usize> = BTreeMap::new();
    let mut spin_point_groups: BTreeMap<&str, usize> = BTreeMap::new();
    let mut tspin_counts: BTreeMap<i8, usize> = BTreeMap::new();
    let mut gs_counts: BTreeMap<i8, usize> = BTreeMap::new();
    let mut operation_counts: BTreeMap<usize, usize> = BTreeMap::new();
    let mut total_operations = 0usize;
    let mut preserve_operations = 0usize;
    let mut swap_operations = 0usize;
    let mut records_with_swap = 0usize;
    let mut previous = None;

    for record in collinear_spin_space_groups() {
        if let Some(previous) = previous {
            assert!(
                previous < record.number,
                "records are not strictly sorted: {previous} then {}",
                record.number
            );
        }
        previous = Some(record.number);
        assert!(
            identifiers.insert(record.identifier),
            "duplicate identifier {}",
            record.identifier
        );
        assert!(
            numbers.insert(record.number),
            "duplicate number {}",
            record.number
        );
        assert_eq!(record.identifier, record.number.to_string());
        let parsed = SpinSpaceGroupNumber::parse(record.identifier)
            .unwrap_or_else(|| panic!("unparsable identifier {}", record.identifier));
        assert_eq!(parsed, record.number);
        assert_eq!(parsed.g0(), record.number.g0());
        assert_eq!(parsed.l0(), record.number.l0());
        assert_eq!(parsed.k_index(), record.number.k_index());
        assert_eq!(parsed.map_index(), record.number.map_index());
        assert_eq!(parsed.letter(), record.number.letter());
        let by_number = collinear_spin_space_group(record.number)
            .unwrap_or_else(|| panic!("record {} is not reachable by number", record.identifier));
        assert_eq!(by_number.identifier, record.identifier);

        assert_eq!(record.configuration, SpinGroupKind::Collinear);
        assert!(!record.operations.is_empty(), "{}", record.identifier);
        assert_eq!(
            record.num_all_elements as usize,
            record.operations.len(),
            "{}: num_all_elements must equal the quotient signature length",
            record.identifier
        );
        assert!(!record.international_symbol.is_empty());
        assert!(
            record
                .international_symbol
                .chars()
                .all(|character| character.is_ascii_alphanumeric()
                    || matches!(character, '-' | '_' | '/' | '∞')),
            "{}: plain symbol {:?} contains unexpected characters",
            record.identifier,
            record.international_symbol
        );
        assert!(!record.international_notation_latex.is_empty());
        assert_eq!(record.gso, "∞m", "{}", record.identifier);
        assert!(matches!(record.tspin, -1 | 1), "{}", record.identifier);
        assert!(matches!(record.gs, -1 | 1), "{}", record.identifier);

        *group_types.entry(record.group_type).or_default() += 1;
        *bases.entry(record.basis).or_default() += 1;
        *spin_point_groups
            .entry(record.spin_point_group)
            .or_default() += 1;
        *tspin_counts.entry(record.tspin).or_default() += 1;
        *gs_counts.entry(record.gs).or_default() += 1;
        *operation_counts.entry(record.operations.len()).or_default() += 1;
        total_operations += record.operations.len();
        let mut has_swap = false;
        for operation in record.operations {
            match operation.channel {
                CollinearAction::Preserve => preserve_operations += 1,
                CollinearAction::Swap => {
                    swap_operations += 1;
                    has_swap = true;
                }
            }
        }
        if has_swap {
            records_with_swap += 1;
        }
    }

    assert_eq!(identifiers.len(), EXPECTED_RECORD_COUNT);
    assert_eq!(numbers.len(), EXPECTED_RECORD_COUNT);
    assert_eq!(total_operations, EXPECTED_TOTAL_OPERATIONS);
    assert_eq!(preserve_operations, EXPECTED_PRESERVE_OPERATIONS);
    assert_eq!(swap_operations, EXPECTED_SWAP_OPERATIONS);
    assert_eq!(records_with_swap, EXPECTED_RECORDS_WITH_SWAP);
    assert_counts(&group_types, EXPECTED_GROUP_TYPE_COUNTS);
    assert_counts(&bases, EXPECTED_BASIS_COUNTS);
    assert_counts(&spin_point_groups, EXPECTED_SPIN_POINT_GROUP_COUNTS);
    assert_counts(&tspin_counts, EXPECTED_TSPIN_COUNTS);
    assert_counts(&gs_counts, EXPECTED_GS_COUNTS);
    let histogram: Vec<(usize, usize)> = operation_counts.into_iter().collect();
    assert_eq!(histogram, EXPECTED_OPERATION_COUNT_HISTOGRAM.to_vec());
    assert_eq!(
        EXPECTED_OPERATION_COUNT_HISTOGRAM
            .iter()
            .map(|(operations, records)| operations * records)
            .sum::<usize>(),
        EXPECTED_TOTAL_OPERATIONS
    );
    assert_eq!(
        EXPECTED_OPERATION_COUNT_HISTOGRAM
            .iter()
            .map(|(_, records)| records)
            .sum::<usize>(),
        EXPECTED_RECORD_COUNT
    );
}

#[test]
fn operations_are_well_formed_and_unique_per_record() {
    let identity_rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for record in collinear_spin_space_groups() {
        let mut keys = BTreeSet::new();
        for operation in record.operations {
            let determinant = rotation_determinant(&operation.rotation);
            assert_eq!(
                determinant.abs(),
                1,
                "{}: spatial rotation is not unimodular ({determinant})",
                record.identifier
            );
            assert!(
                keys.insert(operation_key(operation)),
                "{}: duplicate spatial operation {:?}",
                record.identifier,
                operation.rotation
            );
        }
        let identity = find_operation(
            record,
            identity_rotation,
            [0, 0, 0],
            CollinearAction::Preserve,
        )
        .unwrap_or_else(|| {
            panic!(
                "{}: missing identity operation with Preserve channel",
                record.identifier
            )
        });
        assert_eq!(identity.rotation, identity_rotation);
    }
}

#[test]
fn spatial_projection_is_closed_and_channel_is_an_exact_homomorphism() {
    for record in collinear_spin_space_groups() {
        let mut table: BTreeMap<OperationKey, i32> = BTreeMap::new();
        for operation in record.operations {
            table.insert(operation_key(operation), sign(operation.channel));
        }
        let keys: Vec<OperationKey> = table.keys().copied().collect();
        for left in &keys {
            let left_sign = table[left];
            for right in &keys {
                let product_rotation = rotation_product(&left.0, &right.0);
                let product_translation = translation_product(&left.0, left.1, right.1);
                let product_key = (product_rotation, product_translation);
                let product_sign = table.get(&product_key).unwrap_or_else(|| {
                    panic!(
                        "{}: quotient signature is not closed under {:?} * {:?}",
                        record.identifier, left.0, right.0
                    )
                });
                assert_eq!(
                    *product_sign,
                    left_sign * table[right],
                    "{}: channel is not a homomorphism",
                    record.identifier
                );
            }
        }
    }
}

#[test]
fn quotient_signatures_are_globally_unique() {
    let mut signatures = BTreeSet::new();
    for record in collinear_spin_space_groups() {
        assert!(
            signatures.insert(canonical_signature(record)),
            "duplicate quotient signature for {}",
            record.identifier
        );
    }
    assert_eq!(signatures.len(), EXPECTED_RECORD_COUNT);
}

#[test]
fn representative_records_match_raw_export_fixtures() {
    let mut types = BTreeSet::new();
    let mut bases = BTreeSet::new();
    let mut sizes = BTreeSet::new();
    let mut saw_swap = false;
    let mut saw_pure_preserve = false;
    for fixture in FIXTURES {
        let record = fixture_record(fixture);
        types.insert(record.group_type);
        bases.insert(record.basis);
        sizes.insert(record.operations.len());
        assert_eq!(record.identifier, fixture.identifier);
        assert_eq!(
            (
                record.number.g0(),
                record.number.l0(),
                record.number.k_index(),
                record.number.map_index(),
                record.number.letter()
            ),
            fixture.number
        );
        assert_eq!(record.international_symbol, fixture.international_symbol);
        assert_eq!(
            record.international_notation_latex,
            fixture.international_notation_latex
        );
        assert_eq!(record.group_type, fixture.group_type);
        assert_eq!(record.basis, fixture.basis);
        assert_eq!(record.spin_point_group, fixture.spin_point_group);
        assert_eq!(record.tspin, fixture.tspin);
        assert_eq!(record.gs, fixture.gs);
        assert_eq!(record.num_all_elements, fixture.num_all_elements);
        assert_eq!(record.operations.len(), fixture.num_all_elements as usize);
        let mut preserve = 0usize;
        let mut swap = 0usize;
        for operation in record.operations {
            match operation.channel {
                CollinearAction::Preserve => preserve += 1,
                CollinearAction::Swap => swap += 1,
            }
        }
        assert_eq!(preserve, fixture.preserve_ops, "{}", fixture.identifier);
        assert_eq!(swap, fixture.swap_ops, "{}", fixture.identifier);
        if swap > 0 {
            saw_swap = true;
        }
        if preserve == record.operations.len() {
            saw_pure_preserve = true;
        }
        assert_eq!(
            sha256_hex(&canonical_record_bytes(record)),
            fixture.record_sha256,
            "{}: canonical record digest differs from the raw-export oracle",
            fixture.identifier
        );
    }

    // Coverage required by the review brief.
    assert!(types.contains("typeI") && types.contains("typeII"));
    assert!(bases.contains("G0") && bases.contains("L0"));
    assert!(sizes.contains(&1) && sizes.contains(&48) && sizes.contains(&384));
    assert!(saw_swap && saw_pure_preserve);

    // Explicit operations taken verbatim from the raw export: a swap mirror
    // and its preserve counterpart, an anti-translation, and an inversion.
    let swap_mirror = [[0, -1, 0], [-1, 0, 0], [0, 0, 1]];
    let inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let record_200 = collinear_spin_space_group_by_identifier("221.200.1.1.L").unwrap();
    assert!(find_operation(record_200, swap_mirror, [0, 0, 0], CollinearAction::Swap).is_some());
    assert!(find_operation(record_200, inversion, [0, 0, 0], CollinearAction::Preserve).is_some());
    let record_221 = collinear_spin_space_group_by_identifier("221.221.1.1.L").unwrap();
    assert!(
        find_operation(
            record_221,
            swap_mirror,
            [0, 0, 0],
            CollinearAction::Preserve
        )
        .is_some()
    );
    let record_anti_translation = collinear_spin_space_group_by_identifier("1.1.2.1.L").unwrap();
    assert!(
        find_operation(
            record_anti_translation,
            identity,
            [6, 0, 0],
            CollinearAction::Swap
        )
        .is_some()
    );
    assert!(
        find_operation(
            record_anti_translation,
            identity,
            [0, 0, 0],
            CollinearAction::Preserve
        )
        .is_some()
    );
    for identifier in ["221.225.2.1.L", "221.226.2.1.L"] {
        let record = collinear_spin_space_group_by_identifier(identifier).unwrap();
        assert!(
            find_operation(record, inversion, [0, 0, 0], CollinearAction::Preserve).is_some(),
            "{identifier}: missing inversion Preserve operation"
        );
        assert!(
            find_operation(record, inversion, [0, 0, 6], CollinearAction::Swap).is_some(),
            "{identifier}: missing half-translation inversion Swap operation"
        );
    }
}

#[test]
fn canonical_table_digest_matches_raw_export_oracle() {
    let digest = sha256_hex(&canonical_table_bytes());
    assert_eq!(digest.len(), 64);
    assert_eq!(
        digest, EXPECTED_TABLE_SHA256,
        "the full-table canonical digest differs from the raw-export oracle"
    );
}

#[test]
fn number_parsing_and_lookup_negative_controls() {
    for text in [
        "",
        "1",
        "1.1.1.1",
        "1.1.1.1.L.0",
        ".1.1.1.L",
        "1..1.1.L",
        "1.1.1.1.",
        "0.1.1.1.L",
        "1.0.1.1.L",
        "1.1.0.1.L",
        "1.1.1.0.L",
        "231.1.1.1.L",
        "1.231.1.1.L",
        "1.1.49.1.L",
        "1.1.9.1.L",
        "1.1.1.49.L",
        "1.1.1.1.LL",
        "1.1.1.1.1",
        "-1.1.1.1.L",
        "a.b.c.d.L",
        "1.1.1.1.L ",
        " 1.1.1.1.L",
    ] {
        assert!(
            SpinSpaceGroupNumber::parse(text).is_none(),
            "parse accepted malformed identifier {text:?}"
        );
    }

    assert!(SpinSpaceGroupNumber::new(0, 1, 1, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(231, 1, 1, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 0, 1, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 231, 1, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 1, 0, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 1, 49, 1, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 1, 1, 0, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 1, 1, 49, 'L').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 1, 1, 1, '1').is_none());
    assert!(SpinSpaceGroupNumber::new(1, 1, 1, 1, '∞').is_none());

    let edge = SpinSpaceGroupNumber::new(230, 230, 8, 48, 'Z').expect("edge number is in range");
    assert_eq!(edge.to_string(), "230.230.8.48.Z");
    assert!(collinear_spin_space_group(edge).is_none());
    assert!(collinear_spin_space_group_by_identifier("230.230.48.8.Z").is_none());
    // The table only contains the `L` family letter: a valid-range number with
    // a different letter must not resolve, nor may the identifier lookup be
    // case-insensitive.
    let absent_letter = SpinSpaceGroupNumber::new(1, 1, 1, 1, 'M').unwrap();
    assert!(collinear_spin_space_group(absent_letter).is_none());
    assert!(collinear_spin_space_group_by_identifier("1.1.1.1.M").is_none());
    assert!(collinear_spin_space_group_by_identifier("1.1.1.1.l").is_none());
    assert!(collinear_spin_space_group_by_identifier("not-an-identifier").is_none());

    let smallest = SpinSpaceGroupNumber::parse("1.1.1.1.L").expect("canonical identifier");
    assert_eq!(
        smallest,
        SpinSpaceGroupNumber::new(1, 1, 1, 1, 'L').unwrap()
    );
    assert_eq!(smallest.to_string(), "1.1.1.1.L");
    assert_eq!(
        collinear_spin_space_group(smallest).map(|record| record.identifier),
        Some("1.1.1.1.L")
    );
}
