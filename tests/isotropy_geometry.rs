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
    IsotropyDirection, basis_in_parent_conventional, centering_multiplicity,
    double_valued_subduction, format_identity_subduction, format_isotropy_subgroups,
    format_magnetic_isotropy_subgroups, identity_subduction, isotropy_subgroup_for_direction,
    isotropy_subgroups, isotropy_subgroups_at_k, k_vectors_agree,
    magnetic_isotropy_subgroup_for_direction, magnetic_isotropy_subgroups,
    origin_shift_in_parent_conventional, parent_primitive_basis, subgroup_size,
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

#[test]
fn ordinals_address_their_own_records() {
    for sg in 1..=230u8 {
        for irrep in query::irreps_of(sg) {
            if irrep.spinor {
                continue;
            }
            for subgroup in isotropy_subgroups(sg, irrep.ml).expect("subgroups") {
                let stored = &ISOTROPY_SUBGROUPS[subgroup.ordinal];
                assert_eq!(stored.sg, subgroup.record.sg, "SG {sg} {}", irrep.ml);
                assert_eq!(stored.direction, subgroup.record.direction);
                assert_eq!(stored.basis, subgroup.record.basis);
                assert_eq!(stored.origin, subgroup.record.origin);
            }
        }
    }
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
fn errors_are_displayable() {
    let errors = [
        isotropy_subgroups(0, "GM1").unwrap_err(),
        isotropy_subgroups(221, "NOPE").unwrap_err(),
        isotropy_subgroups_at_k(221, KVector::new([1, 0, 0], 2), "GM3+").unwrap_err(),
        isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(zzz)"))
            .unwrap_err(),
        subgroup_size([[0, 0, 0], [0, 0, 0], [0, 0, 0]]).unwrap_err(),
        identity_subduction(usize::MAX).unwrap_err(),
        origin_shift_in_parent_conventional(221, [0, 0, 0, 0]).unwrap_err(),
    ];
    for error in errors {
        let text = error.to_string();
        assert!(!text.is_empty(), "{error:?}");
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
            assert_eq!(subgroups.len(), irrep.subgroups().len());
            for subgroup in &subgroups {
                ordinals.push(subgroup.ordinal);
            }
        }
    }
    ordinals.sort_unstable();
    assert_eq!(ordinals, (0..ISOTROPY_SUBGROUPS.len()).collect::<Vec<_>>());
}

/// SG 225 `W5` along `S60`: the ISOTROPY program's `SHOW FREQ DIR` prints nine
/// scalar entries (ending `3 W5 S60(1)`) plus three double-valued ones
/// (`3 DT5`, `3 SM3`, `3 SM4`).  The double-valued table was missing from the
/// API before this gate existed.
#[test]
fn double_valued_subduction_matches_the_program() {
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
    let double = subgroup
        .double_valued_subduction()
        .expect("double-valued subduction");
    let actual: Vec<(&str, u8, u16)> = double
        .iter()
        .map(|entry| (entry.parent_ml, entry.parent_sg, entry.frequency))
        .collect();
    assert_eq!(
        actual,
        vec![("DT5", 225, 3), ("SM3", 225, 3), ("SM4", 225, 3)]
    );

    let table = format_identity_subduction(subgroup.ordinal).expect("table");
    assert!(table.contains("| S60 |"), "{table}");
    assert!(table.contains("DT5"), "{table}");
}

#[test]
fn double_valued_subduction_tables_tile() {
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
        let entries = double_valued_subduction(ordinal).expect("resolves");
        let expected =
            (ISOTROPY_W_SUBDUCE_RANGES[ordinal + 1] - ISOTROPY_W_SUBDUCE_RANGES[ordinal]) as usize;
        assert_eq!(entries.len(), expected, "ordinal {ordinal}");
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
    assert_eq!(subgroup_size(subgroup.record.basis).unwrap(), 2);
    // The ISOTROPY program prints the F-centred conventional cell
    // (2,0,0),(0,2,0),(0,0,2) for this transition, whose volume is
    // Z(F222) * Size / Z(P222) = 4 * 2 / 1 = 8 parent conventional cells.
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

#[test]
fn out_of_range_space_groups_fail_closed() {
    assert!(isotropy_subgroups(0, "GM1").is_err());
    assert!(isotropy_subgroups(231, "GM1").is_err());
    assert!(parent_primitive_basis(231).is_err());
}
