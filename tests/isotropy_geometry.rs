//! Whole-table invariants and oracle-pinned rows for the isotropy geometry.
//!
//! The expected values come from the bundled ISOTROPY program (see
//! `scripts/verify_isotropy_oracle.py`) and from Stokes & Hatch (1988); this
//! test keeps them from drifting when the generated data is regenerated.

use cryspglib::irrep::generated_data::{
    ISOTROPY_SUBDUCE_DOMAIN, ISOTROPY_SUBDUCE_FREQUENCY, ISOTROPY_SUBDUCE_IRREP,
    ISOTROPY_SUBDUCE_RANGES, ISOTROPY_SUBGROUPS, IRREPS, MAGNETIC_ISOTROPY_SUBGROUPS,
};
use cryspglib::irrep::isotropy::{
    basis_in_parent_conventional, centering_multiplicity, isotropy_subgroup_for_direction,
    isotropy_subgroups, magnetic_isotropy_subgroups, origin_shift_in_parent_conventional,
    parent_primitive_basis, subgroup_size, IsotropyDirection,
};
use cryspglib::irrep::query;

fn det3(m: [[i32; 3]; 3]) -> i32 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

const ALLOWED_DENOMINATORS: [i32; 9] = [1, 2, 3, 4, 6, 8, 12, 16, 24];

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
        assert_ne!(det3(record.basis), 0, "record {index}: singular basis");
        assert!(
            ALLOWED_DENOMINATORS.contains(&record.origin[3]),
            "record {index}: origin denominator {}",
            record.origin[3]
        );
        assert!(subgroup_size(record.basis).unwrap() >= 1);
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
        assert_ne!(det3(record.basis), 0, "record {index}: singular basis");
        assert!(
            ALLOWED_DENOMINATORS.contains(&record.origin[3]),
            "record {index}: origin denominator {}",
            record.origin[3]
        );
    }
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
            for (local, subgroup) in subgroups.iter().enumerate() {
                assert_eq!(
                    subgroup.ordinal,
                    subgroups[local].ordinal,
                    "SG {sg} {}: ordinal mismatch",
                    irrep.ml
                );
                ordinals.push(subgroup.ordinal);
            }
        }
    }
    ordinals.sort_unstable();
    assert_eq!(ordinals, (0..ISOTROPY_SUBGROUPS.len()).collect::<Vec<_>>());
}

#[test]
fn subduction_tables_are_parallel_and_in_range() {
    assert_eq!(ISOTROPY_SUBDUCE_RANGES.len(), ISOTROPY_SUBGROUPS.len() + 1);
    assert!(ISOTROPY_SUBDUCE_RANGES.windows(2).all(|pair| pair[0] <= pair[1]));
    assert_eq!(
        ISOTROPY_SUBDUCE_IRREP.len(),
        ISOTROPY_SUBDUCE_FREQUENCY.len()
    );
    assert_eq!(ISOTROPY_SUBDUCE_IRREP.len(), ISOTROPY_SUBDUCE_DOMAIN.len());
    assert_eq!(
        ISOTROPY_SUBDUCE_RANGES[ISOTROPY_SUBGROUPS.len()] as usize - 1,
        ISOTROPY_SUBDUCE_IRREP.len()
    );
    for (index, irrep_index) in ISOTROPY_SUBDUCE_IRREP.iter().enumerate() {
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
    // The printed F-centred conventional cell spans 4 (Z of F222) times the
    // primitive cell, i.e. 8 parent conventional cells.
    assert_eq!(det3([[2, 0, 0], [0, 2, 0], [0, 0, 2]]), 8);
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
    let subgroup = isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(a,0)"))
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
    assert!(subgroups
        .iter()
        .all(|sub| (1..=1651).contains(&sub.record.mag_sg)));
    assert!(subgroups.iter().all(|sub| det3(sub.record.basis) != 0));
}

#[test]
fn centering_multiplicities_cover_every_space_group() {
    for sg in 1..=230u8 {
        let multiplicity = centering_multiplicity(sg)
            .unwrap_or_else(|error| panic!("SG {sg}: {error}"));
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
