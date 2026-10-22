//! Public arithmetic regressions: exact limits are not non-containment proofs.
use cryspglib::irrep::magnetic_embedding::{
    GeometryError, MagneticContractError, canonical_translation, canonical_translation_checked,
    geometry_of, magnetic_operations, measure_parent_containment, search_parent_setting,
    verify_embedding, verify_group,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Mat3R, Rat};

fn geometry(parent: u8, uni: usize) -> cryspglib::irrep::magnetic_embedding::MagneticGeometry {
    let record = query::magnetic_isotropy_subgroups_of(parent)
        .into_iter()
        .find(|row| row.subgroup.mag_sg == uni)
        .unwrap()
        .subgroup;
    geometry_of(parent, &record).unwrap()
}

#[test]
fn reduction_discards_large_integer_parts_before_scaling() {
    for numerator in [i128::MAX, i128::MIN] {
        let point = [
            Rat::new(numerator, 3).unwrap(),
            Rat::new(1, 5).unwrap(),
            Rat::ZERO,
        ];
        let expected = [
            Rat::new(numerator.rem_euclid(3), 3).unwrap(),
            point[1],
            Rat::ZERO,
        ];
        assert_eq!(canonical_translation(point, &[]), expected);
    }
}

#[test]
fn representable_lattice_rejection_does_not_overflow_on_large_numerators() {
    let map = Mat3R::new([
        [Rat::new(i128::MAX, 3).unwrap(), Rat::ZERO, Rat::ZERO],
        [Rat::new(1, 5).unwrap(), Rat::ONE, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, Rat::ONE],
    ]);
    assert!(matches!(
        verify_embedding(
            &geometry(1, 1),
            &magnetic_operations(1).unwrap(),
            map,
            [Rat::ZERO; 3]
        ),
        Err(GeometryError::SourceLatticeNotInParent { uni: 1 })
    ));
}

#[test]
fn independent_fractional_axes_are_not_rounded_into_the_parent_lattice() {
    let map = Mat3R::new([
        [Rat::new(1, i128::MAX).unwrap(), Rat::ZERO, Rat::ZERO],
        [Rat::new(1, i128::MAX - 2).unwrap(), Rat::ONE, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, Rat::ONE],
    ]);
    assert_eq!(
        verify_embedding(
            &geometry(1, 1),
            &magnetic_operations(1).unwrap(),
            map,
            [Rat::ZERO; 3]
        ),
        Err(GeometryError::SourceLatticeNotInParent { uni: 1 })
    );
}

#[test]
fn group_verification_does_not_round_an_invalid_antiunitary_square_to_identity() {
    let mut set = magnetic_operations(2).unwrap();
    set.operations[1].translation = [
        Rat::new(1, i128::MAX).unwrap(),
        Rat::new(1, i128::MAX - 2).unwrap(),
        Rat::ZERO,
    ];
    assert!(matches!(
        verify_group(&set),
        Err(MagneticContractError::NotAGroup { uni: 2, .. })
    ));
}

#[test]
fn overflow_in_class_placement_is_not_reported_as_not_in_parent() {
    assert_eq!(
        verify_embedding(
            &geometry(2, 4),
            &magnetic_operations(4).unwrap(),
            Mat3R::identity(),
            [Rat::from_integer(i128::MAX), Rat::ZERO, Rat::ZERO],
        ),
        Err(GeometryError::Arithmetic)
    );
}

#[test]
fn valid_large_inputs_are_not_blanket_rejected() {
    let mut set = magnetic_operations(4).unwrap();
    // An inversion with any rational translation still squares to identity.
    set.operations[1].translation = [
        Rat::new(1, i128::MAX).unwrap(),
        Rat::new(1, i128::MAX - 2).unwrap(),
        Rat::ZERO,
    ];
    verify_group(&set).unwrap();
    verify_embedding(
        &geometry(1, 1),
        &magnetic_operations(1).unwrap(),
        Mat3R::identity(),
        [
            Rat::from_integer(i128::MAX),
            Rat::from_integer(-i128::MAX),
            Rat::ZERO,
        ],
    )
    .unwrap();
}

#[test]
fn an_unrepresentable_intermediate_negation_stays_unknown() {
    // Rat::checked_sub cannot negate MIN, even in MIN - MIN. Do not certify
    // non-containment when that existing exact-arithmetic boundary is reached.
    assert_eq!(
        verify_embedding(
            &geometry(1, 1),
            &magnetic_operations(1).unwrap(),
            Mat3R::identity(),
            [Rat::from_integer(i128::MIN), Rat::ZERO, Rat::ZERO]
        ),
        Err(GeometryError::Arithmetic)
    );
}

#[test]
fn modular_addition_near_i128_max_stays_exact() {
    let denominator = 5 * (1i128 << 124);
    let point = [
        Rat::new(denominator - 1, denominator).unwrap(),
        Rat::ZERO,
        Rat::ZERO,
    ];
    let lattice = [[Rat::new(4, 5).unwrap(), Rat::ZERO, Rat::ZERO]];
    assert_eq!(
        canonical_translation_checked(point, &lattice).unwrap(),
        [
            Rat::new(denominator / 5 - 1, denominator).unwrap(),
            Rat::ZERO,
            Rat::ZERO
        ]
    );
}

#[test]
fn checked_reduction_matches_independent_complete_coefficient_enumeration() {
    // Grid 60: these two coupled generators have exact orders six and twenty.
    // Enumerate their coefficient product directly, without the production BFS.
    let lattice = [[30, 20, 0], [0, 15, 12]].map(|g| g.map(|v| Rat::new(v, 60).unwrap()));
    for x in [-61i128, 0, 1, 29, 59, 120] {
        for y in [-61i128, 0, 1, 29, 59, 120] {
            for z in [-61i128, 0, 1, 29, 59, 120] {
                let mut best = [60; 3];
                for a in 0..6 {
                    for b in 0..20 {
                        best = best.min([
                            (x + 30 * a).rem_euclid(60),
                            (y + 20 * a + 15 * b).rem_euclid(60),
                            (z + 12 * b).rem_euclid(60),
                        ]);
                    }
                }
                let point = [x, y, z].map(|v| Rat::new(v, 60).unwrap());
                assert_eq!(
                    canonical_translation_checked(point, &lattice).unwrap(),
                    best.map(|v| Rat::new(v, 60).unwrap())
                );
            }
        }
    }
}

#[test]
fn a_coordinate_lcm_overflow_is_an_explicit_unknown_result() {
    let point = [Rat::new(1, i128::MAX).unwrap(), Rat::ZERO, Rat::ZERO];
    let lattice = [[Rat::new(1, 2).unwrap(), Rat::ZERO, Rat::ZERO]];
    assert_eq!(
        canonical_translation_checked(point, &lattice),
        Err(MagneticContractError::Arithmetic)
    );
}

#[test]
fn the_reduction_budget_accepts_its_boundary_and_reports_excess_as_unknown() {
    let lattice = |denominator| [[Rat::new(1, denominator).unwrap(), Rat::ZERO, Rat::ZERO]];
    assert_eq!(
        canonical_translation_checked([Rat::ZERO; 3], &lattice(65_536)),
        Ok([Rat::ZERO; 3])
    );
    assert_eq!(
        canonical_translation_checked([Rat::ZERO; 3], &lattice(65_537)),
        Err(MagneticContractError::ReductionCapacity { limit: 65_536 })
    );
    let mut set = magnetic_operations(1).unwrap();
    let mut translation = set.operations[0];
    translation.translation = lattice(65_537)[0];
    set.operations.push(translation);
    set.unitary_count += 1;
    assert_eq!(
        verify_group(&set),
        Err(MagneticContractError::ReductionCapacity { limit: 65_536 })
    );
}

#[test]
fn diagnostic_entry_points_propagate_arithmetic_limits() {
    let mut geometry = geometry(2, 4);
    geometry.origin = [Rat::from_integer(i128::MAX), Rat::ZERO, Rat::ZERO];
    let set = magnetic_operations(4).unwrap();
    assert_eq!(
        measure_parent_containment(&geometry, &set),
        Err(GeometryError::Arithmetic)
    );
    assert_eq!(
        search_parent_setting(&geometry, &set),
        Err(GeometryError::Arithmetic)
    );
}

#[test]
fn scalar_input_preserves_an_embedding_arithmetic_refusal() {
    use cryspglib::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, MagneticInputError, RealSign, RecordIdentity,
        ScalarGammaRequest, SourceIdentity, verify_scalar_gamma_input,
    };
    use cryspglib::irrep::{LabelConvention, isotropy};
    let direction_index = isotropy::magnetic_isotropy_subgroups(2, "GM1+", LabelConvention::Cdml)
        .unwrap()
        .iter()
        .position(|row| row.record.mag_sg == 4)
        .unwrap();
    let request = ScalarGammaRequest::new(
        SourceIdentity {
            parent_sg: 2,
            label: "GM1+",
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: 2,
            condensing_label: "GM1+",
            convention: LabelConvention::Cdml,
            direction_index,
        },
        Some(RealSign::Plus),
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::identity(),
        shift: [Rat::from_integer(i128::MAX), Rat::ZERO, Rat::ZERO],
        lattice_policy: LatticePolicy::ContainmentOnly,
    });
    assert!(matches!(
        verify_scalar_gamma_input(&request),
        Err(MagneticInputError::EmbeddingRejected(
            GeometryError::Arithmetic
        ))
    ));
}
