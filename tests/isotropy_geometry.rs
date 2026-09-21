//! Whole-table invariants and oracle-pinned rows for the isotropy geometry.
//!
//! The expected values come from the bundled ISOTROPY program (see
//! `scripts/verify_isotropy_oracle.py`) and from Stokes & Hatch (1988); this
//! test keeps them from drifting when the generated data is regenerated.

use cryspglib::irrep::generated_data::{
    IRREP_W_LABELS, IRREP_W_SPACE_GROUP, IRREPS, ISOTROPY_DIRECTION_LABELS,
    ISOTROPY_SUBDUCE_DIRECTION, ISOTROPY_SUBDUCE_DOMAIN, ISOTROPY_SUBDUCE_FREQUENCY,
    ISOTROPY_SUBDUCE_IRREP, ISOTROPY_SUBDUCE_RANGES, ISOTROPY_SUBGROUPS,
    ISOTROPY_W_SUBDUCE_FREQUENCY, ISOTROPY_W_SUBDUCE_IRREP, ISOTROPY_W_SUBDUCE_RANGES,
    MAGNETIC_ISOTROPY_SUBGROUPS,
};
use cryspglib::irrep::isotropy::{
    IsotropyDirection, IsotropyError, basis_in_parent_conventional, centering_multiplicity,
    format_identity_subduction, format_isotropy_subgroups, format_magnetic_isotropy_subgroups,
    identity_subduction, isotropy_subgroup_for_direction, isotropy_subgroups,
    isotropy_subgroups_at_k, k_vectors_agree, magnetic_isotropy_subgroup_for_direction,
    magnetic_isotropy_subgroups, origin_shift_in_parent_conventional, other_wave_vector_subduction,
    parent_primitive_basis, subgroup_size,
};
use cryspglib::irrep::query;
use cryspglib::irrep::types::KVector;

fn det3(m: [[i32; 3]; 3]) -> i32 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

const ALLOWED_DENOMINATORS: [i32; 9] = [1, 2, 3, 4, 6, 8, 12, 16, 24];

/// Determinants that occur in the pinned isotropy tables.
const ALLOWED_DETERMINANTS: [i32; 8] = [1, 2, 3, 4, 6, 8, 16, 32];

fn gcd4(a: i32, b: i32, c: i32, d: i32) -> i32 {
    fn gcd(mut a: i32, mut b: i32) -> i32 {
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a.abs()
    }
    gcd(gcd(a.abs(), b.abs()), gcd(c.abs(), d.abs()))
}

#[test]
fn every_isotropy_record_has_valid_geometry() {
    assert_eq!(ISOTROPY_SUBGROUPS.len(), 15239);
    for (index, record) in ISOTROPY_SUBGROUPS.iter().enumerate() {
        assert!(
            (1..=230).contains(&record.sg),
            "record {index}: subgroup {} out of range",
            record.sg
        );
        assert!(!record.symbol.is_empty(), "record {index}: empty symbol");
        assert!(
            !record.direction.is_empty() && !record.direction_label.is_empty(),
            "record {index}: empty direction"
        );
        assert!(
            (1..=24).contains(&record.direction_dim),
            "record {index}: direction dimension {}",
            record.direction_dim
        );
        assert!(record.domains >= 1 && record.arms >= 1, "record {index}");
        let det = det3(record.basis);
        assert!(
            ALLOWED_DETERMINANTS.contains(&det),
            "record {index}: unexpected basis determinant {det}"
        );
        assert!(
            ALLOWED_DENOMINATORS.contains(&record.origin[3]),
            "record {index}: origin denominator {}",
            record.origin[3]
        );
        let o = record.origin;
        assert_eq!(
            gcd4(o[0], o[1], o[2], o[3]),
            1,
            "record {index}: origin {o:?} is not in lowest terms"
        );
        assert_eq!(subgroup_size(record.basis).unwrap(), det as u32);
    }
}

#[test]
fn every_magnetic_isotropy_record_has_valid_geometry() {
    assert_eq!(MAGNETIC_ISOTROPY_SUBGROUPS.len(), 16721);
    for (index, record) in MAGNETIC_ISOTROPY_SUBGROUPS.iter().enumerate() {
        assert!(
            (1..=1651).contains(&record.mag_sg),
            "record {index}: UNI {} out of range",
            record.mag_sg
        );
        assert!(!record.bns_label.is_empty() && !record.iso_label.is_empty());
        // A mis-tokenised label section shows up as a label carrying quote
        // marks or whitespace: the reader used to fold the next labels into
        // one token (`'" "P_S1        " "P-1         " "P-11'`), which left
        // 15546 of these records with the wrong BNS symbol.
        for (field, value) in [
            ("bns_label", record.bns_label),
            ("iso_label", record.iso_label),
        ] {
            assert!(
                !value.contains('"') && !value.chars().any(char::is_whitespace),
                "record {index}: {field} {value:?} is not a single symbol"
            );
        }
        assert!(!record.direction.is_empty());
        assert!(
            ALLOWED_DETERMINANTS.contains(&det3(record.basis)),
            "record {index}: unexpected basis determinant"
        );
        assert!(
            ALLOWED_DENOMINATORS.contains(&record.origin[3]),
            "record {index}: origin denominator {}",
            record.origin[3]
        );
        let o = record.origin;
        assert_eq!(gcd4(o[0], o[1], o[2], o[3]), 1, "record {index}: {o:?}");
        // Regression for the synthetic `dir<code>` placeholders that a
        // misaligned pointer produced for 72% of the magnetic records.
        assert!(
            !record.direction.starts_with("dir"),
            "record {index}: synthetic direction {:?}",
            record.direction
        );
    }
}

/// The BNS symbols are indexed by UNI, so a label list that is even one token
/// too long shifts every record.  Pin the values that the pre-fix reader
/// corrupted, so a regression cannot pass on a `!is_empty()` check.
#[test]
fn magnetic_bns_symbols_are_indexed_by_uni() {
    // UNI 6 is BNS 2.6 `P-1'` (the apostrophe is what broke the old reader).
    let bns_of = |uni: usize| {
        MAGNETIC_ISOTROPY_SUBGROUPS
            .iter()
            .find(|record| record.mag_sg == uni)
            .map(|record| record.bns_label)
            .unwrap_or_else(|| panic!("no magnetic isotropy record for UNI {uni}"))
    };
    assert_eq!(bns_of(3), "P_S1");
    assert_eq!(bns_of(6), "P-1'");
    assert_eq!(bns_of(1001), "P4/m'mm");
    assert_eq!(bns_of(1651), "Ia'-3'd'");
    // Every label list entry must be a single symbol, for all UNIs present.
    for record in MAGNETIC_ISOTROPY_SUBGROUPS {
        assert!(
            !record.bns_label.contains('"') && !record.bns_label.contains("  "),
            "UNI {}: {:?}",
            record.mag_sg,
            record.bns_label
        );
    }
}

#[test]
fn magnetic_direction_labels_are_unique_within_each_irrep() {
    let mut checked = 0usize;
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor {
                continue;
            }
            let subgroups = magnetic_isotropy_subgroups(sg, irrep.ml)
                .unwrap_or_else(|error| panic!("SG {sg} {}: {error}", irrep.ml));
            let mut labels: Vec<&str> = subgroups.iter().map(|s| s.record.direction).collect();
            let total = labels.len();
            labels.sort_unstable();
            labels.dedup();
            assert_eq!(
                labels.len(),
                total,
                "SG {sg} {}: duplicate magnetic direction labels",
                irrep.ml
            );
            checked += total;
        }
    }
    assert_eq!(checked, MAGNETIC_ISOTROPY_SUBGROUPS.len());
}

/// Source-indexed values for SG 9 `L1` (magnetic isotropy table): `P1` is
/// UNI 48, `P3` and `C1` are UNI 3 with different bases.  Before the pointer
/// fix, `Label("P3")` returned the `P1` record.
#[test]
fn magnetic_direction_selection_uses_the_per_record_label() {
    let p1 = magnetic_isotropy_subgroup_for_direction(9, "L1", IsotropyDirection::Label("P1"))
        .expect("L1 has a P1 direction");
    assert_eq!(p1.record.mag_sg, 48);
    assert_eq!(p1.record.basis, [[-1, -1, 1], [1, 1, 1], [-1, 1, 1]]);
    assert_eq!(p1.record.origin, [1, 1, 1, 4]);

    let p3 = magnetic_isotropy_subgroup_for_direction(9, "L1", IsotropyDirection::Label("P3"))
        .expect("L1 has a P3 direction");
    assert_eq!(p3.record.mag_sg, 3);
    assert_eq!(p3.record.basis, [[0, 0, 1], [1, -1, 0], [2, 0, 0]]);
    assert_eq!(p3.record.origin, [0, 0, 0, 1]);

    let c1 = magnetic_isotropy_subgroup_for_direction(9, "L1", IsotropyDirection::Label("C1"))
        .expect("L1 has a C1 direction");
    assert_eq!(c1.record.mag_sg, 3);
    assert_eq!(c1.record.basis, [[1, 1, 1], [-1, -1, 1], [2, 0, 0]]);
    assert_ne!(p3.ordinal, c1.ordinal);
}

/// Ordinals are 0-based positions in the flat isotropy table.  Pin a handful
/// against `data_isotropy.txt` (`isotropy_irrep_pointer` is 1-based, so an
/// irrep's records start at `pointer[k-1]-1`); a loop that compares the API
/// record with the very array the API sliced it from cannot fail, so this test
/// carries the coverage instead.
#[test]
fn pinned_ordinals_address_the_expected_records() {
    let cases: [(u8, &str, IsotropyDirection<'_>, usize, usize); 6] = [
        (16, "R1", IsotropyDirection::Label("P1"), 314, 22),
        (139, "M1-", IsotropyDirection::Label("P1"), 6213, 126),
        (221, "GM3+", IsotropyDirection::Label("P1"), 12397, 123),
        (221, "GM3+", IsotropyDirection::Label("C1"), 12398, 47),
        (221, "GM4+", IsotropyDirection::Label("P1"), 12400, 83),
        (221, "GM4+", IsotropyDirection::Label("S1"), 12402, 2),
    ];
    for (sg, ml, direction, ordinal, subgroup_sg) in cases {
        let subgroup = isotropy_subgroup_for_direction(sg, ml, direction)
            .unwrap_or_else(|error| panic!("SG {sg} {ml} {direction:?}: {error}"));
        assert_eq!(subgroup.ordinal, ordinal, "SG {sg} {ml} {direction:?}");
        assert_eq!(subgroup.record.sg, subgroup_sg, "SG {sg} {ml}");
        // The ordinal must address the same record in the flat generated table.
        let stored = &ISOTROPY_SUBGROUPS[ordinal];
        assert_eq!(stored.sg, subgroup_sg, "ordinal {ordinal}");
        assert_eq!(stored.basis, subgroup.record.basis, "ordinal {ordinal}");
        assert_eq!(stored.origin, subgroup.record.origin, "ordinal {ordinal}");
        // ... and `identity_subduction` must resolve that same ordinal.
        assert!(identity_subduction(ordinal).is_ok(), "ordinal {ordinal}");
    }
    // Consecutive records of one irrep keep their source order and labels.
    assert_eq!(ISOTROPY_SUBGROUPS[12397].direction_label, "P1");
    assert_eq!(ISOTROPY_SUBGROUPS[12398].direction_label, "C1");
}

#[test]
fn formatters_render_geometry_and_subduction() {
    let table = format_isotropy_subgroups(221, "GM3+").expect("ordinary table");
    assert!(table.contains("#123 P4/mmm"));
    assert!(
        !table.contains('?'),
        "formatter hid a conversion error: {table}"
    );

    let magnetic = format_magnetic_isotropy_subgroups(9, "L1").expect("magnetic table");
    assert!(magnetic.contains("| P3 |"), "{magnetic}");
    assert!(!magnetic.contains('?'), "{magnetic}");

    let subgroup =
        isotropy_subgroup_for_direction(221, "GM4+", IsotropyDirection::Descriptor("(a,0,0)"))
            .expect("Γ4+ (a,0,0)");
    let subduction = format_identity_subduction(subgroup.ordinal).expect("subduction table");
    assert!(subduction.contains("GM1+"), "{subduction}");
    assert!(subduction.contains("GM3+"), "{subduction}");
    assert!(subduction.contains("| 1 |"), "{subduction}");
}

#[test]
fn errors_report_their_variant_and_are_displayable() {
    // A message alone would accept any wrong variant, so assert the variant
    // (and the payload where it carries one) for every reachable error path.
    assert!(matches!(
        isotropy_subgroups(0, "GM1"),
        Err(IsotropyError::InvalidSpaceGroup(0))
    ));
    assert!(matches!(
        isotropy_subgroups(221, "NOPE"),
        Err(IsotropyError::UnknownIrrep { .. })
    ));
    assert!(matches!(
        isotropy_subgroups_at_k(221, KVector::new([1, 0, 0], 2), "GM3+"),
        Err(IsotropyError::IrrepNotAtKPoint { .. })
    ));
    assert!(matches!(
        isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(zzz)")),
        Err(IsotropyError::DirectionNotFound { .. })
    ));
    assert!(matches!(
        subgroup_size([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
        Err(IsotropyError::SingularSubgroupBasis { .. })
    ));
    assert!(matches!(
        subgroup_size([[1 << 21, 0, 0], [0, 1 << 21, 0], [0, 0, 1 << 22]]),
        Err(IsotropyError::SubgroupSizeOverflow { determinant }) if determinant == 1i128 << 64
    ));
    assert!(matches!(
        identity_subduction(usize::MAX),
        Err(IsotropyError::InvalidIsotropyRecord { .. })
    ));
    assert!(matches!(
        other_wave_vector_subduction(usize::MAX),
        Err(IsotropyError::InvalidIsotropyRecord { .. })
    ));
    assert!(matches!(
        origin_shift_in_parent_conventional(221, [0, 0, 0, 0]),
        Err(IsotropyError::InvalidOrigin { origin }) if origin == [0, 0, 0, 0]
    ));
    assert!(matches!(
        isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Index(usize::MAX)),
        Err(IsotropyError::SubgroupIndexOutOfRange { .. })
    ));
    // `MissingCentering` is unreachable for sg in 1..=230 (the centering is
    // derived from the Hall database, see the coverage test below); the spinor
    // path is asserted in the library's own unit tests.

    let errors = [
        isotropy_subgroups(0, "GM1").unwrap_err(),
        isotropy_subgroups(221, "NOPE").unwrap_err(),
        isotropy_subgroups_at_k(221, KVector::new([1, 0, 0], 2), "GM3+").unwrap_err(),
        isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(zzz)"))
            .unwrap_err(),
        subgroup_size([[0, 0, 0], [0, 0, 0], [0, 0, 0]]).unwrap_err(),
        subgroup_size([[1 << 21, 0, 0], [0, 1 << 21, 0], [0, 0, 1 << 22]]).unwrap_err(),
        identity_subduction(usize::MAX).unwrap_err(),
        origin_shift_in_parent_conventional(221, [0, 0, 0, 0]).unwrap_err(),
        isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Index(usize::MAX))
            .unwrap_err(),
    ];
    for error in errors {
        let text = error.to_string();
        assert!(!text.is_empty(), "{error:?}");
        assert!(!text.contains('?'), "message hides a value: {text}");
        let _: &dyn std::error::Error = &error;
    }
}

#[test]
fn wave_vector_helpers_reject_degenerate_input() {
    // A zero denominator is not a wave vector and never agrees with anything.
    assert!(!k_vectors_agree(
        KVector::new([0, 0, 1], 0),
        KVector::new([0, 0, 1], 1)
    ));
    assert!(isotropy_subgroups_at_k(221, KVector::new([0, 0, 0], 0), "GM3+").is_err());
    // Reduction and sign normalisation.
    assert!(k_vectors_agree(
        KVector::new([0, 0, 2], 4),
        KVector::new([0, 0, 1], 2)
    ));
    assert!(k_vectors_agree(
        KVector::new([0, 0, -1], -2),
        KVector::new([0, 0, 1], 2)
    ));
    assert!(!k_vectors_agree(
        KVector::new([0, 0, 1], 2),
        KVector::new([0, 0, 1], 3)
    ));
}

#[test]
fn isotropy_ranges_partition_the_table_per_irrep() {
    // Records live in ISOTROPY source order, which need not follow the emitted
    // irrep order, so check the partition as a set of ordinals.
    let mut ordinals = Vec::new();
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor {
                // Spinor irreps have no isotropy data and must fail closed
                // rather than panic.
                assert!(isotropy_subgroups(sg, irrep.ml).is_err());
                continue;
            }
            let subgroups = isotropy_subgroups(sg, irrep.ml)
                .unwrap_or_else(|error| panic!("SG {sg} {}: {error}", irrep.ml));
            for subgroup in &subgroups {
                ordinals.push(subgroup.ordinal);
            }
        }
    }
    ordinals.sort_unstable();
    assert_eq!(ordinals, (0..ISOTROPY_SUBGROUPS.len()).collect::<Vec<_>>());
}

/// SG 225 `W5` along `S60`: the ISOTROPY program's `SHOW FREQ DIR` prints nine
/// entries at the selected wave vector (ending `3 W5 S60(1)`) plus three at
/// *other* wave vectors (`3 DT5`, `3 SM3`, `3 SM4`).  Those three are
/// single-valued irreps on the `DT` and `SM` lines of SG 225 — not spinors —
/// which is what these entries are named after.
#[test]
fn other_wave_vector_subduction_matches_the_program() {
    let subgroup = isotropy_subgroup_for_direction(225, "W5", IsotropyDirection::Label("S60"))
        .expect("W5 has an S60 direction");
    let scalar = subgroup.identity_subduction().expect("scalar subduction");
    let actual: Vec<(&str, u16, &str)> = scalar
        .iter()
        .map(|entry| (entry.parent_ml, entry.frequency, entry.direction_label))
        .collect();
    assert_eq!(
        actual,
        vec![
            ("GM1+", 1, "P1"),
            ("GM2+", 1, "P1"),
            ("GM3+", 2, "C1"),
            ("L1-", 1, "P11"),
            ("L2-", 1, "P11"),
            ("L3-", 2, "C42"),
            ("X1+", 3, "S1"),
            ("X2+", 3, "S1"),
            ("W5", 3, "S60"),
        ]
    );
    let other = subgroup
        .other_wave_vector_subduction()
        .expect("other-wave-vector subduction");
    let actual: Vec<(&str, u8, u16)> = other
        .iter()
        .map(|entry| (entry.parent_ml, entry.parent_sg, entry.frequency))
        .collect();
    assert_eq!(
        actual,
        vec![("DT5", 225, 3), ("SM3", 225, 3), ("SM4", 225, 3)]
    );
    let table = format_identity_subduction(subgroup.ordinal).expect("table");
    assert!(table.contains("| S60 |"), "{table}");
    assert!(
        table.ends_with(
            "| Other-wave-vector parent irrep | i(G) |\n\
             |--------------------------------|------|\n\
             | DT5 (SG 225) | 3 |\n\
             | SM3 (SG 225) | 3 |\n\
             | SM4 (SG 225) | 3 |"
        ),
        "{table}"
    );
}

/// SHOW DIRECTION VECTOR pins the labels; data_isotropy.txt pins the raw
/// primitive bases.  A subgroup number alone cannot distinguish the two Cm
/// embeddings, and normalizing separators must not permute components.
#[test]
fn descriptor_aliases_preserve_the_direction_and_embedding() {
    let cases = [
        (
            225,
            "GM4-",
            "(a,b,0)",
            "C1",
            8,
            [[1, -1, 0], [0, -1, 1], [-1, 0, 0]],
        ),
        (
            225,
            "GM4-",
            "(a,a,b)",
            "C2",
            8,
            [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
        ),
        (
            177,
            "L1",
            "(a;b;0)",
            "C1",
            5,
            [[0, -2, 0], [1, -1, 1], [-1, 1, 1]],
        ),
        (
            177,
            "L1",
            "(a;b;a)",
            "C2",
            21,
            [[0, -2, 0], [2, 2, 0], [0, 0, 2]],
        ),
    ];
    for (parent, ml, descriptor, label, sg, basis) in cases {
        let aliases = [
            descriptor.to_string(),
            descriptor.replace(';', ","),
            descriptor.replace(',', ";"),
            format!(
                " \t{}\n",
                descriptor.replace(',', " ,\t").replace(';', " ;\t")
            ),
        ];
        for alias in aliases {
            let subgroup =
                isotropy_subgroup_for_direction(parent, ml, IsotropyDirection::Descriptor(&alias))
                    .unwrap_or_else(|error| panic!("SG {parent} {ml} {alias:?}: {error}"));
            assert_eq!(subgroup.record.direction_label, label, "{alias:?}");
            assert_eq!(subgroup.record.sg, sg, "{alias:?}");
            assert_eq!(subgroup.record.basis, basis, "{alias:?}");
        }
    }
    for (parent, ml, descriptor) in [
        (221, "GM4+", "(0,a,0)"),
        (225, "GM4-", "(a,0,b)"),
        (225, "GM4-", "(a,b,a)"),
        (177, "L1", "(a,a,b)"),
    ] {
        assert!(
            matches!(
                isotropy_subgroup_for_direction(
                    parent,
                    ml,
                    IsotropyDirection::Descriptor(descriptor)
                ),
                Err(IsotropyError::DirectionNotFound { .. })
            ),
            "SG {parent} {ml}: {descriptor} must not select a different component pattern"
        );
    }
}

#[test]
fn other_wave_vector_subduction_tables_tile() {
    assert_eq!(
        ISOTROPY_W_SUBDUCE_RANGES.len(),
        ISOTROPY_SUBGROUPS.len() + 1
    );
    assert!(
        ISOTROPY_W_SUBDUCE_RANGES
            .windows(2)
            .all(|pair| pair[0] <= pair[1])
    );
    assert_eq!(
        ISOTROPY_W_SUBDUCE_RANGES[ISOTROPY_SUBGROUPS.len()] as usize,
        ISOTROPY_W_SUBDUCE_IRREP.len()
    );
    assert_eq!(
        ISOTROPY_W_SUBDUCE_IRREP.len(),
        ISOTROPY_W_SUBDUCE_FREQUENCY.len()
    );
    assert_eq!(IRREP_W_LABELS.len(), IRREP_W_SPACE_GROUP.len());
    for index in 0..ISOTROPY_W_SUBDUCE_IRREP.len() {
        let irrep = ISOTROPY_W_SUBDUCE_IRREP[index] as usize;
        assert!((1..=IRREP_W_LABELS.len()).contains(&irrep));
        assert!(ISOTROPY_W_SUBDUCE_FREQUENCY[index] >= 1);
    }
    for ordinal in 0..ISOTROPY_SUBGROUPS.len() {
        let entries = other_wave_vector_subduction(ordinal).expect("resolves");
        let expected =
            (ISOTROPY_W_SUBDUCE_RANGES[ordinal + 1] - ISOTROPY_W_SUBDUCE_RANGES[ordinal]) as usize;
        assert_eq!(entries.len(), expected, "ordinal {ordinal}");
        // Compare content, not just the length derived from the same range: an
        // off-by-one in the range/pointer arithmetic must fail here.  The W
        // ranges are 0-based start offsets (unlike the 1-based scalar ones).
        for (offset, entry) in entries.iter().enumerate() {
            let packed = ISOTROPY_W_SUBDUCE_RANGES[ordinal] as usize + offset;
            let label_index = ISOTROPY_W_SUBDUCE_IRREP[packed] as usize - 1;
            assert_eq!(
                entry.parent_ml, IRREP_W_LABELS[label_index],
                "packed {packed}"
            );
            assert_eq!(
                entry.parent_sg, IRREP_W_SPACE_GROUP[label_index],
                "packed {packed}"
            );
            assert_eq!(
                entry.frequency, ISOTROPY_W_SUBDUCE_FREQUENCY[packed] as u16,
                "packed {packed}"
            );
        }
    }
}

#[test]
fn subduction_tables_are_parallel_and_in_range() {
    assert_eq!(ISOTROPY_SUBDUCE_RANGES.len(), ISOTROPY_SUBGROUPS.len() + 1);
    assert!(
        ISOTROPY_SUBDUCE_RANGES
            .windows(2)
            .all(|pair| pair[0] <= pair[1])
    );
    assert_eq!(
        ISOTROPY_SUBDUCE_IRREP.len(),
        ISOTROPY_SUBDUCE_FREQUENCY.len()
    );
    assert_eq!(ISOTROPY_SUBDUCE_IRREP.len(), ISOTROPY_SUBDUCE_DOMAIN.len());
    assert_eq!(
        ISOTROPY_SUBDUCE_RANGES[ISOTROPY_SUBGROUPS.len()] as usize - 1,
        ISOTROPY_SUBDUCE_IRREP.len()
    );
    assert!(
        ISOTROPY_DIRECTION_LABELS
            .iter()
            .all(|label| !label.is_empty())
    );
    for (index, irrep_index) in ISOTROPY_SUBDUCE_IRREP.iter().enumerate() {
        assert!((ISOTROPY_SUBDUCE_DIRECTION[index] as usize) < ISOTROPY_DIRECTION_LABELS.len());
        assert!(
            (1..=IRREPS.len()).contains(&(*irrep_index as usize)),
            "entry {index}: irrep index {irrep_index} out of range"
        );
        assert!(ISOTROPY_SUBDUCE_FREQUENCY[index] >= 1);
        assert!(!IRREPS[*irrep_index as usize - 1].spinor);
        // The record anchor behind `ISOTROPY_SUBDUCE_DIRECTION` is a record of
        // the *same parent space group* but usually of a different irrep
        // (15238 of 94271 entries), so "same irrep" is not the invariant.  The
        // anchor itself is not emitted; the generator gates it.
    }
}

/// #16 `P222`, irrep `R1` along `(a)` → #22 `F222`, as printed by the ISOTROPY
/// program: Size 2, conventional basis `(2,0,0),(0,2,0),(0,0,2)`.
#[test]
fn p222_r1_reaches_f222() {
    let subgroup = isotropy_subgroup_for_direction(16, "R1", IsotropyDirection::Label("P1"))
        .expect("R1 has a P1 direction");
    assert_eq!(subgroup.record.sg, 22);
    assert_eq!(subgroup.record.basis, [[0, 1, 1], [1, 0, 1], [1, 1, 0]]);
    assert_eq!(subgroup.record.origin, [0, 0, 0, 1]);
    assert_eq!(subgroup_size(subgroup.record.basis).unwrap(), 2);
    // The program prints the F-centred conventional cell
    // `(2,0,0),(0,2,0),(0,0,2)` for this transition.  Its primitive lattice is
    // that cell times the subgroup's own centring matrix, i.e. `P_F . printed`
    // (row `j` of `P_F` holds the coefficients of the printed rows), which must
    // reproduce the stored basis in the parent's (here primitive) frame.  This
    // pins the printed basis and the frame convention, where a determinant
    // identity alone would not.
    let printed = [[2.0f64, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    let p_f = parent_primitive_basis(22).unwrap();
    let mut expected = [[0.0f64; 3]; 3];
    for j in 0..3 {
        for i in 0..3 {
            expected[j][i] = (0..3).map(|t| p_f[j][t] * printed[t][i]).sum();
        }
    }
    let stored = basis_in_parent_conventional(16, subgroup.record.basis).unwrap();
    assert_eq!(
        stored, expected,
        "printed cell must rebuild the stored basis"
    );
    // Volume relation printed by the program: |det Basis| = Z(sub)*Size/Z(parent).
    assert_eq!(
        centering_multiplicity(22).unwrap() * subgroup_size(subgroup.record.basis).unwrap()
            / centering_multiplicity(16).unwrap(),
        8
    );
}

/// #167 `R-3c`, Γ3+ along `(a,0)` → #15 `C2/c`.  The stored primitive-frame
/// origin `(0, 1/2, 0)` is `(-1/6, 1/6, 1/6)` in hexagonal conventional
/// coordinates, which is what the ISO program prints.
#[test]
fn r3c_gm3plus_origin_converts_to_hexagonal_axes() {
    let subgroup = isotropy_subgroup_for_direction(167, "GM3+", IsotropyDirection::Label("P1"))
        .expect("Γ3+ has a P1 direction");
    assert_eq!(subgroup.record.sg, 15);
    assert_eq!(subgroup.record.basis, [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]);
    assert_eq!(subgroup.record.origin, [0, 1, 0, 2]);
    assert_eq!(subgroup_size(subgroup.record.basis).unwrap(), 1);
    let converted = origin_shift_in_parent_conventional(167, subgroup.record.origin).unwrap();
    let expected = [-1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];
    for (got, want) in converted.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-12, "{converted:?} != {expected:?}");
    }
}

/// #225 `Fm-3m`, Γ4- along `(a,0,0)` → #107 `I4mm` with the parent's lattice.
#[test]
fn fm3m_gm4minus_along_a0_keeps_the_lattice() {
    let subgroup =
        isotropy_subgroup_for_direction(225, "GM4-", IsotropyDirection::Descriptor("(a,0,0)"))
            .expect("Γ4- has an (a,0,0) direction");
    assert_eq!(subgroup.record.sg, 107);
    assert_eq!(subgroup.record.direction_label, "P1");
    assert_eq!(subgroup_size(subgroup.record.basis).unwrap(), 1);
    // Same lattice as the parent, expressed in the parent's primitive frame.
    assert_eq!(subgroup.record.basis, [[0, 1, 0], [-1, 0, 1], [1, -1, 0]]);
    assert_eq!(subgroup.record.origin, [0, 0, 0, 1]);
    assert_eq!(
        parent_primitive_basis(225).unwrap(),
        [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    );
}

/// #221 `Pm-3m`, Γ3+ along `(a,0)` → #123 `P4/mmm` with the parent's cell.
#[test]
fn pm3m_gm3plus_along_a0_keeps_the_cell() {
    let subgroup =
        isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(a,0)"))
            .expect("Γ3+ has an (a,0) direction");
    assert_eq!(subgroup.record.sg, 123);
    assert_eq!(subgroup.record.basis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    assert_eq!(subgroup.record.origin, [0, 0, 0, 1]);
    assert_eq!(subgroup_size(subgroup.record.basis).unwrap(), 1);
    assert_eq!(
        basis_in_parent_conventional(221, subgroup.record.basis).unwrap(),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    );
}

/// Magnetic isotropy records are keyed by non-magnetic parent irreps and end in
/// UNI numbers, exactly like the pinned `data_magnetic.txt`.
#[test]
fn magnetic_isotropy_geometry_is_available() {
    let subgroups = magnetic_isotropy_subgroups(221, "GM4+").expect("mag subgroups");
    assert!(!subgroups.is_empty());
    assert!(subgroups.iter().all(|sub| sub.parent_sg == 221));
    assert!(
        subgroups
            .iter()
            .all(|sub| (1..=1651).contains(&sub.record.mag_sg))
    );
    assert!(subgroups.iter().all(|sub| det3(sub.record.basis) != 0));
}

#[test]
fn centering_multiplicities_cover_every_space_group() {
    for sg in 1..=230u8 {
        let multiplicity =
            centering_multiplicity(sg).unwrap_or_else(|error| panic!("SG {sg}: {error}"));
        assert!((1..=4).contains(&multiplicity), "SG {sg}: Z={multiplicity}");
    }
    assert_eq!(centering_multiplicity(225).unwrap(), 4); // Fm-3m
    assert_eq!(centering_multiplicity(167).unwrap(), 3); // R-3c
    assert_eq!(centering_multiplicity(229).unwrap(), 2); // Im-3m
}

/// `scripts/verify_isotropy_oracle.py` carries its own hand-typed 230-entry
/// centering table, which selects the parent lattice used by every origin
/// comparison.  A drift between the two would silently validate against the
/// wrong lattice, so parse the Python table and tie it to the Rust
/// implementation instead of trusting two copies to stay equal.
#[test]
fn oracle_script_centering_table_matches_the_rust_implementation() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/scripts/verify_isotropy_oracle.py"
    );
    let script = std::fs::read_to_string(path).expect("the oracle script ships with the repo");
    let block = script
        .split("CENTERING_LETTER = {")
        .nth(1)
        .expect("centering table present")
        .split('}')
        .next()
        .expect("table terminator");
    let mut checked = 0usize;
    for entry in block.split(',') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let (key, value) = entry.split_once(':').expect("key: value");
        let sg: u8 = key.trim().parse().expect("space group number");
        let letter = value.trim().trim_matches('"');
        let expected = match letter {
            "P" => 1,
            "A" | "B" | "C" | "I" => 2,
            "R" => 3,
            "F" => 4,
            other => panic!("SG {sg}: unknown centering letter {other:?}"),
        };
        assert_eq!(
            centering_multiplicity(sg).unwrap(),
            expected,
            "SG {sg} is {letter}-centred in the oracle script"
        );
        checked += 1;
    }
    assert_eq!(
        checked, 230,
        "the oracle table must cover every space group"
    );
}

#[test]
fn out_of_range_space_groups_fail_closed() {
    assert!(isotropy_subgroups(0, "GM1").is_err());
    assert!(isotropy_subgroups(231, "GM1").is_err());
    assert!(parent_primitive_basis(231).is_err());
}
