use cryspglib::irrep::isotropy::{self, IsotropyDirection};
use cryspglib::irrep::magnetic_embedding::{
    self as m, GeometryError, MagneticGeometry, MagneticOperationSet,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, query};
use rayon::prelude::*;

fn context(label: &str) -> (MagneticGeometry, MagneticOperationSet) {
    let row = isotropy::magnetic_isotropy_subgroup_for_direction(
        3,
        label,
        LabelConvention::Cdml,
        IsotropyDirection::Label("P1"),
    )
    .unwrap();
    (
        m::geometry_of(3, &row.record).unwrap(),
        m::magnetic_operations(row.record.mag_sg).unwrap(),
    )
}

#[test]
fn primitive_and_centred_source_rotations_give_opposite_obstructions() {
    for (label, uni, source_rank, parent_rank) in [("A1", 11, 0, 1), ("D1", 24, 1, 0)] {
        let (geometry, set) = context(label);
        assert_eq!(set.uni, uni);
        let proof = m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .unwrap();
        assert_eq!(proof.source_rank_mod_two, source_rank);
        assert_eq!(proof.parent_ranks_mod_two, vec![parent_rank]);
        assert_eq!((0..3).map(|i| proof.source_rotation[i][i]).sum::<i32>(), -1);
        assert_eq!(
            Mat3R::from_ints(proof.source_rotation)
                .determinant()
                .unwrap(),
            Rat::ONE
        );
    }
}

#[test]
fn parent_rotations_must_match_spectrum_and_preserve_the_entire_record_lattice() {
    let (mut geometry, set) = context("A1");
    geometry.parent_sg = 1;
    assert!(
        m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .unwrap()
            .parent_ranks_mod_two
            .is_empty()
    );
    geometry.parent_sg = 3;
    geometry.basis = *Mat3R::from_ints([[1, 1, 0], [0, 3, 0], [0, 0, 1]]).rows();
    geometry.index = 3;
    assert!(
        m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .unwrap()
            .parent_ranks_mod_two
            .is_empty()
    );
}

fn binary_matrices() -> impl Iterator<Item = [[i32; 3]; 3]> {
    (0..512u16).map(|bits| {
        std::array::from_fn(|i| std::array::from_fn(|j| i32::from((bits >> (3 * i + j)) & 1)))
    })
}
fn determinant(m: [[i32; 3]; 3]) -> i32 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
fn product_mod_two(a: [[i32; 3]; 3], b: [[i32; 3]; 3]) -> [[i32; 3]; 3] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum::<i32>().rem_euclid(2))
    })
}

#[test]
fn all_168_modular_changes_of_basis_independently_exclude_both_witnesses() {
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let target_a = [[1, 0, 0], [1, 1, 0], [0, 0, 1]];
    let source_d = [[0, 1, 0], [1, 0, 0], [0, 0, 1]];
    let mut checked = 0;
    for u in binary_matrices().filter(|u| determinant(*u).rem_euclid(2) == 1) {
        // S U = U R is tested directly, with neither the rank algorithm nor
        // inversion. No GL(3,F_2) conjugacy means no GL(3,Z) conjugacy.
        assert_ne!(product_mod_two(target_a, u), product_mod_two(u, identity));
        assert_ne!(product_mod_two(identity, u), product_mod_two(u, source_d));
        checked += 1;
    }
    assert_eq!(checked, 168);
}

#[test]
fn the_certificate_is_invariant_under_unimodular_rebases_of_the_same_record_lattice() {
    for label in ["A1", "D1", "E1"] {
        let (geometry, set) = context(label);
        let expected = m::record_lattice_obstruction(&geometry, &set).unwrap();
        let mut checked = 0;
        for u in binary_matrices().filter(|u| determinant(*u).abs() == 1) {
            let mut rebased = geometry.clone();
            rebased.basis = *Mat3R::from_ints(u)
                .checked_mul(&Mat3R::new(geometry.basis))
                .unwrap()
                .rows();
            assert_eq!(
                m::record_lattice_obstruction(&rebased, &set).unwrap(),
                expected
            );
            checked += 1;
        }
        assert_eq!(checked, 168);
    }
}

#[test]
fn no_modular_obstruction_is_not_an_embedding_witness() {
    let (geometry, set) = context("B1");
    assert!(
        m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .is_none()
    );
    assert!(
        m::embed_via_primitive_lattices(&geometry, &set)
            .unwrap()
            .is_empty()
    );
    assert!(
        m::parent_embedding_verdict(&geometry, &set)
            .unwrap()
            .embeddings
            .iter()
            .all(|e| !e.realises_record_lattice)
    );
}

#[test]
fn the_obstruction_does_not_reject_containment_with_a_different_lattice() {
    use cryspglib::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
    };
    use cryspglib::irrep::magnetic_subduction_api::{
        MagneticRecordSelector, MagneticSubductionRequest, subduce,
    };
    let (geometry, set) = context("A1");
    assert!(
        m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .is_some()
    );
    let map = Mat3R::diagonal([2, 2, 2]);
    let shift = [Rat::ZERO; 3];
    m::verify_embedding(&geometry, &set, map, shift).unwrap();
    assert!(!m::realises_record_lattice_checked(&geometry, &set, &map).unwrap());
    for label in ["GM1", "GM2"] {
        for eta in [RealSign::Plus, RealSign::Minus] {
            let request = MagneticSubductionRequest::new(
                SourceIdentity {
                    parent_sg: 3,
                    label,
                    convention: LabelConvention::Cdml,
                },
                MagneticRecordSelector::new(3, "A1", LabelConvention::Cdml)
                    .direction(IsotropyDirection::Label("P1"))
                    .uni(11),
                Some(eta),
            )
            .with_embedding(EmbeddingEvidence {
                map,
                shift,
                lattice_policy: LatticePolicy::ContainmentOnly,
            });
            let report = subduce(&request).unwrap();
            assert_eq!(
                (report.source_dimension(), report.covered_dimension()),
                (1, 1)
            );
            let strict = request.with_embedding(EmbeddingEvidence {
                map,
                shift,
                lattice_policy: LatticePolicy::RequireExact,
            });
            assert!(subduce(&strict).is_err());
        }
    }
}

#[test]
fn malformed_sources_and_unknown_parent_frames_are_errors_not_certificates() {
    let (mut geometry, set) = context("A1");
    geometry.uni += 1;
    assert!(matches!(
        m::record_lattice_obstruction(&geometry, &set),
        Err(GeometryError::PrimitiveBridgeUnavailable { .. })
    ));
    geometry.uni = set.uni;
    let mut altered = set.clone();
    altered.operations.swap(0, 1);
    assert!(matches!(
        m::record_lattice_obstruction(&geometry, &altered),
        Err(GeometryError::PrimitiveBridgeUnavailable { .. })
    ));
    altered = set.clone();
    altered.unitary_count += 1;
    assert!(matches!(
        m::record_lattice_obstruction(&geometry, &altered),
        Err(GeometryError::PrimitiveBridgeUnavailable { .. })
    ));
    geometry.parent_sg = 0;
    assert!(m::record_lattice_obstruction(&geometry, &set).is_err());
}

#[test]
fn arithmetic_and_singular_lattices_cannot_turn_into_absence_certificates() {
    let (mut geometry, set) = context("A1");
    geometry.basis = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                Rat::from_integer(i128::MAX)
            } else {
                Rat::ZERO
            }
        })
    });
    assert_eq!(
        m::record_lattice_obstruction(&geometry, &set),
        Err(GeometryError::Arithmetic)
    );
    geometry.basis = [[Rat::ZERO; 3]; 3];
    assert_eq!(
        m::record_lattice_obstruction(&geometry, &set),
        Err(GeometryError::Arithmetic)
    );
}

#[test]
fn all_records_separate_56_literal_lattice_obstructions_from_unresolved_cases() {
    let (records, obstructed) = (1u8..=230)
        .into_par_iter()
        .map(|sg| {
            let (mut records, mut obstructed) = (0, 0);
            for irrep in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
                for row in
                    isotropy::magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml)
                        .unwrap()
                {
                    records += 1;
                    let geometry = m::geometry_of(sg, &row.record).unwrap();
                    let set = m::magnetic_operations(row.record.mag_sg).unwrap();
                    if let Some(proof) = m::record_lattice_obstruction(&geometry, &set).unwrap() {
                        assert!((3..=15).contains(&sg));
                        assert!(
                            !proof
                                .parent_ranks_mod_two
                                .contains(&proof.source_rank_mod_two)
                        );
                        obstructed += 1;
                    }
                }
            }
            (records, obstructed)
        })
        .reduce(|| (0, 0), |(a, b), (c, d)| (a + c, b + d));
    assert_eq!((records, obstructed), (16721, 56));
}
