use cryspglib::irrep::isotropy::{self, IsotropyDirection};
use cryspglib::irrep::magnetic_embedding::{
    self as m, EmbeddingKind, GeometryError, MagneticGeometry, MagneticOperationSet,
};
use cryspglib::irrep::subduction::{Mat3R, Rat, SeitzTransform, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::{LabelConvention, query};

fn context(
    sg: u8,
    label: &str,
    direction: &str,
) -> (MagneticGeometry, MagneticOperationSet, SubgroupEmbedding) {
    let magnetic = isotropy::magnetic_isotropy_subgroup_for_direction(
        sg,
        label,
        LabelConvention::Cdml,
        IsotropyDirection::Label(direction),
    )
    .unwrap();
    let ordinary = isotropy::isotropy_subgroup_for_direction(
        sg,
        label,
        LabelConvention::Cdml,
        IsotropyDirection::Label(direction),
    )
    .unwrap();
    (
        m::geometry_of(sg, &magnetic.record).unwrap(),
        m::magnetic_operations(magnetic.record.mag_sg).unwrap(),
        SubgroupEmbedding::from_isotropy_subgroup(&ordinary).unwrap(),
    )
}

#[test]
fn a_real_bridge_witness_lies_outside_the_certified_r24_family() {
    let (geometry, set, ordinary) = context(3, "E1", "P1");
    assert_eq!(set.uni, 24);
    assert_eq!(ordinary.subgroup_sg(), 5);
    let old = m::parent_embedding_verdict(&geometry, &set).unwrap();
    assert!(old.embeddings.is_empty() && old.absence_is_certified);
    let entries = m::embed_via_unitary_subgroup(&geometry, &set, &ordinary).unwrap();
    let expected = Mat3R::from_ints([[0, 0, 1], [0, 2, 0], [-2, 0, 1]]);
    assert!(
        entries
            .iter()
            .any(|e| e.map == expected && e.shift == [Rat::ZERO; 3])
    );
    for entry in entries {
        assert_eq!(entry.kind, EmbeddingKind::UnitarySubgroup);
        assert!(entry.realises_record_lattice);
        assert!(m::realises_record_lattice_checked(&geometry, &set, &entry.map).unwrap());
        m::verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap();
    }
}

#[test]
fn nonzero_origins_and_centred_lattices_are_solved_not_copied() {
    for (sg, label, direction, denominator) in [(48, "X1", "P1", 4), (227, "GM4+", "P3", 8)] {
        let (geometry, set, ordinary) = context(sg, label, direction);
        let entries = m::embed_via_unitary_subgroup(&geometry, &set, &ordinary).unwrap();
        assert!(
            entries
                .iter()
                .any(|e| e.shift == [Rat::new(1, denominator).unwrap(); 3])
        );
        for entry in entries {
            m::verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap();
            assert!(m::realises_record_lattice_checked(&geometry, &set, &entry.map).unwrap());
        }
    }
}

#[test]
fn identities_and_raw_magnetic_data_cannot_be_substituted() {
    let (geometry, set, ordinary) = context(3, "E1", "P1");
    let (_, _, other_parent) = context(1, "GM1", "P1");
    let (_, _, other_h) = context(3, "GM1", "P1");
    for wrong in [&other_parent, &other_h] {
        assert!(matches!(
            m::embed_via_unitary_subgroup(&geometry, &set, wrong),
            Err(GeometryError::UnitaryBridgeUnavailable { .. })
        ));
    }
    let mut wrong = geometry.clone();
    wrong.uni += 1;
    assert!(matches!(
        m::embed_via_unitary_subgroup(&wrong, &set, &ordinary),
        Err(GeometryError::UnitaryBridgeUnavailable { .. })
    ));
    let mut reordered = set.clone();
    reordered.operations.swap(0, 1);
    assert!(matches!(
        m::embed_via_unitary_subgroup(&geometry, &reordered, &ordinary),
        Err(GeometryError::UnitaryBridgeUnavailable { .. })
    ));
}

#[test]
fn an_empty_bridge_does_not_disprove_other_exact_embeddings() {
    let (mut geometry, set, ordinary) = context(3, "E1", "P1");
    let original = m::embed_via_unitary_subgroup(&geometry, &set, &ordinary)
        .unwrap()
        .into_iter()
        .find(|e| e.shift == [Rat::ZERO; 3])
        .unwrap();
    for row in &mut geometry.basis {
        for value in row {
            *value = value.checked_mul(Rat::from_integer(2)).unwrap();
        }
    }
    geometry.index *= 8;
    assert!(
        m::embed_via_unitary_subgroup(&geometry, &set, &ordinary)
            .unwrap()
            .is_empty()
    );
    let another = original
        .map
        .checked_mul(&Mat3R::from_ints([[2, 0, 0], [0, 2, 0], [0, 0, 2]]))
        .unwrap();
    assert!(m::realises_record_lattice_checked(&geometry, &set, &another).unwrap());
    m::verify_embedding(&geometry, &set, another, [Rat::ZERO; 3]).unwrap();
}

#[test]
fn arithmetic_limits_are_not_empty_solution_sets() {
    let (mut geometry, set, ordinary) = context(3, "E1", "P1");
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
        m::embed_via_unitary_subgroup(&geometry, &set, &ordinary),
        Err(GeometryError::Arithmetic)
    );
}

#[test]
fn recovered_embeddings_reconstruct_nontrivial_gamma_actions_for_both_parities() {
    use cryspglib::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
    };
    use cryspglib::irrep::magnetic_subduction_api::{
        MagneticRecordSelector, MagneticSubductionRequest, subduce,
    };
    use num_complex::Complex64;
    let (geometry, set, ordinary) = context(3, "E1", "P1");
    let entry = m::embed_via_unitary_subgroup(&geometry, &set, &ordinary)
        .unwrap()
        .into_iter()
        .find(|e| e.shift == [Rat::ZERO; 3])
        .unwrap();
    let transform = SeitzTransform::new(entry.map, Vec3R::new(entry.shift));
    let condensing = query::irreps_of(3).iter().find(|r| r.ml == "E1").unwrap();
    let mut cases = 0;
    let mut negative_unitary = 0;
    let mut anti = 0;
    for source in query::irreps_of(3)
        .iter()
        .filter(|r| r.ml == "GM1" || r.ml == "GM2")
    {
        let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
        for eta in [RealSign::Plus, RealSign::Minus] {
            for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                let source_label = source.label(convention).unwrap();
                let record_label = condensing.label(convention).unwrap();
                let request = MagneticSubductionRequest::new(
                    SourceIdentity {
                        parent_sg: 3,
                        label: &source_label,
                        convention,
                    },
                    MagneticRecordSelector::new(3, &record_label, convention)
                        .direction(IsotropyDirection::Label("P1"))
                        .uni(24),
                    Some(eta),
                )
                .with_embedding(EmbeddingEvidence {
                    map: entry.map,
                    shift: entry.shift,
                    lattice_policy: LatticePolicy::RequireExact,
                });
                let report = subduce(&request).unwrap();
                let model = report.scalar_model().unwrap();
                let payload = report.scalar().unwrap();
                assert_eq!(model.source_label(), source.ml);
                assert_eq!(model.record_label(), "E1");
                assert_eq!(
                    (report.source_dimension(), report.covered_dimension()),
                    (1, 1)
                );
                assert_eq!(*model.map(), entry.map);
                for (index, action) in model.actions().iter().enumerate() {
                    let rotation = transform.map_rotation(action.class.rotation).unwrap();
                    let position = row
                        .operations()
                        .iter()
                        .position(|op| rotation.iter().flatten().copied().eq(op.rotation))
                        .unwrap();
                    let expected = row.values()[position].re
                        * if action.class.time_reversal {
                            f64::from(eta.value())
                        } else {
                            1.0
                        };
                    assert_eq!(f64::from(action.value.value()), expected);
                    negative_unitary +=
                        usize::from(!action.class.time_reversal && expected == -1.0);
                    anti += usize::from(action.class.time_reversal);
                    for z in [Complex64::new(1.0, 2.0), Complex64::new(-0.7, 0.4)] {
                        let v = payload.intertwiner * z;
                        let lhs = expected
                            * if action.class.time_reversal {
                                v.conj()
                            } else {
                                v
                            };
                        let rhs = payload.intertwiner * payload.operators[index].apply(z);
                        assert!((lhs - rhs).norm() < 1e-12);
                    }
                }
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 8);
    assert!(negative_unitary > 0 && anti > 0);
}

#[test]
fn frame_composition_order_is_not_commutative() {
    let (geometry, set, ordinary) = context(17, "U1", "P3");
    let entries = m::embed_via_unitary_subgroup(&geometry, &set, &ordinary).unwrap();
    let expected = Mat3R::from_ints([[0, 2, 0], [0, 0, -1], [-2, 0, 0]]);
    assert!(!entries.is_empty());
    assert!(entries.iter().all(|e| e.map == expected));
}

#[test]
fn an_h_embedding_and_exact_lattice_do_not_certify_the_antiunitary_coset() {
    let (geometry, set, ordinary) = context(1, "X1", "P1");
    let h = cryspglib::irrep::corep::identify_unitary_subgroup_with_hall(set.uni).unwrap();
    assert_eq!(
        h.msg_to_data.unwrap().basis,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    );
    let map = *ordinary.transform().matrix();
    assert!(m::realises_record_lattice_checked(&geometry, &set, &map).unwrap());
    // In P1 every unprimed operation is a lattice translation. The antiunitary
    // class has R=I and a non-integral image translation; no origin can alter it.
    for operation in set.classes.iter().filter(|op| !op.time_reversal) {
        assert_eq!(operation.rotation, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        assert!(
            map.checked_mul_vector(&Vec3R::new(operation.translation))
                .unwrap()
                .as_array()
                .iter()
                .all(|v| v.is_integer())
        );
    }
    let anti = set.classes.iter().find(|op| op.time_reversal).unwrap();
    assert_eq!(anti.rotation, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    assert!(
        map.checked_mul_vector(&Vec3R::new(anti.translation))
            .unwrap()
            .as_array()
            .iter()
            .any(|v| !v.is_integer())
    );
    assert!(
        m::embed_via_unitary_subgroup(&geometry, &set, &ordinary)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn the_full_record_census_keeps_the_bridge_domain_and_verified_count_explicit() {
    use rayon::prelude::*;
    let totals = (1u8..=230)
        .into_par_iter()
        .map(|sg| {
            let mut count = [0usize; 5];
            for irrep in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
                for magnetic in
                    isotropy::magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml)
                        .unwrap()
                {
                    count[0] += 1;
                    let Ok(ordinary) = isotropy::isotropy_subgroup_for_direction(
                        sg,
                        irrep.ml,
                        LabelConvention::Cdml,
                        IsotropyDirection::Label(magnetic.record.direction),
                    ) else {
                        count[1] += 1;
                        continue;
                    };
                    let h =
                        cryspglib::irrep::corep::identify_unitary_subgroup(magnetic.record.mag_sg)
                            .unwrap();
                    if h != ordinary.record.sg {
                        count[2] += 1;
                        continue;
                    }
                    let ordinary = SubgroupEmbedding::from_isotropy_subgroup(&ordinary).unwrap();
                    let geometry = m::geometry_of(sg, &magnetic.record).unwrap();
                    let set = m::magnetic_operations(magnetic.record.mag_sg).unwrap();
                    let entries =
                        m::embed_via_unitary_subgroup(&geometry, &set, &ordinary).unwrap();
                    if let Some(entry) = entries.first() {
                        count[4] += 1;
                        assert!(entry.realises_record_lattice);
                        m::verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap();
                    } else {
                        count[3] += 1;
                    }
                }
            }
            count
        })
        .reduce(|| [0; 5], |a, b| std::array::from_fn(|i| a[i] + b[i]));
    println!("records/no ordinary/different H/empty/verified: {totals:?}");
    assert_eq!(totals, [16721, 1482, 16, 2209, 13014]);
}
