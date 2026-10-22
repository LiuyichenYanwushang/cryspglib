use cryspglib::irrep::isotropy::{self, IsotropyDirection};
use cryspglib::irrep::magnetic_embedding::{
    self as m, EmbeddingKind, GeometryError, MagneticGeometry, MagneticOperationSet,
};
use cryspglib::irrep::subduction::{Mat3R, Rat, SeitzTransform, Vec3R};
use cryspglib::irrep::{LabelConvention, query};
use rayon::prelude::*;

fn context(sg: u8, label: &str, direction: &str) -> (MagneticGeometry, MagneticOperationSet) {
    let r = isotropy::magnetic_isotropy_subgroup_for_direction(
        sg,
        label,
        LabelConvention::Cdml,
        IsotropyDirection::Label(direction),
    )
    .unwrap();
    (
        m::geometry_of(sg, &r.record).unwrap(),
        m::magnetic_operations(r.record.mag_sg).unwrap(),
    )
}

#[test]
fn every_reference_primitive_lattice_matches_all_raw_source_translations() {
    (1..=1651).into_par_iter().for_each(|uni| {
        let set = m::magnetic_operations(uni).unwrap();
        let reference =
            u8::try_from(cryspglib::msg_database::get_magnetic_spacegroup_type(uni).number)
                .unwrap();
        let geometry = MagneticGeometry {
            parent_sg: reference,
            uni,
            basis: [
                [Rat::ONE, Rat::ZERO, Rat::ZERO],
                [Rat::ZERO, Rat::ONE, Rat::ZERO],
                [Rat::ZERO, Rat::ZERO, Rat::ONE],
            ],
            origin: [Rat::ZERO; 3],
            index: 1,
        };
        assert!(
            m::realises_record_lattice_checked(&geometry, &set, &Mat3R::identity()).unwrap(),
            "UNI {uni}"
        );
        // Empty is legitimate (e.g. a type-IV half translation cannot fit the
        // undoubled reference parent). This is NOT a 1651-embedding assertion.
        m::embed_via_primitive_lattices(&geometry, &set).unwrap();
    });
}

#[test]
fn centred_source_coordinates_require_the_source_primitive_inverse() {
    let (geometry, set) = context(21, "T1", "P1");
    assert_eq!(set.uni, 150);
    let naive = Mat3R::new(m::basis_in_parent_conventional(&geometry).unwrap()).transpose();
    assert!(!m::realises_record_lattice_checked(&geometry, &set, &naive).unwrap());
    let entries = m::embed_via_primitive_lattices(&geometry, &set).unwrap();
    let expected = Mat3R::from_ints([[0, -1, 0], [1, 0, 0], [0, 0, 2]]);
    assert!(
        entries
            .iter()
            .any(|e| e.map == expected && e.shift == [Rat::ZERO; 3])
    );
    for entry in entries {
        assert_eq!(entry.kind, EmbeddingKind::PrimitiveLattice);
        assert!(entry.realises_record_lattice);
        assert!(m::realises_record_lattice_checked(&geometry, &set, &entry.map).unwrap());
        m::verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap();
    }
}

#[test]
fn lattice_normalization_still_requires_a_nonzero_common_origin() {
    let (geometry, set) = context(21, "T3", "P1");
    let map = Mat3R::from_ints([[0, -1, 0], [1, 0, 0], [0, 0, 2]]);
    assert!(m::realises_record_lattice_checked(&geometry, &set, &map).unwrap());
    assert!(m::verify_embedding(&geometry, &set, map, [Rat::ZERO; 3]).is_err());
    let entries = m::embed_via_primitive_lattices(&geometry, &set).unwrap();
    assert!(
        entries
            .iter()
            .any(|e| e.map == map && e.shift == [Rat::ZERO, Rat::new(1, 4).unwrap(), Rat::ZERO])
    );
    for entry in entries {
        m::verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap();
    }
}

#[test]
fn raw_data_metadata_and_uni_are_verified_before_using_the_reference_cell() {
    let (mut geometry, set) = context(21, "T1", "P1");
    geometry.uni += 1;
    assert!(matches!(
        m::embed_via_primitive_lattices(&geometry, &set),
        Err(GeometryError::PrimitiveBridgeUnavailable { .. })
    ));
    geometry.uni = set.uni;
    let mut corrupt = set.clone();
    corrupt.operations.swap(0, 1);
    assert!(matches!(
        m::embed_via_primitive_lattices(&geometry, &corrupt),
        Err(GeometryError::PrimitiveBridgeUnavailable { .. })
    ));
    corrupt = set;
    corrupt.unitary_count += 1;
    assert!(matches!(
        m::embed_via_primitive_lattices(&geometry, &corrupt),
        Err(GeometryError::PrimitiveBridgeUnavailable { .. })
    ));
}

#[test]
fn no_witness_for_this_basis_order_is_not_global_nonexistence() {
    let (geometry, set) = context(19, "Y1", "P3");
    assert!(
        m::embed_via_primitive_lattices(&geometry, &set)
            .unwrap()
            .is_empty()
    );
    let other = Mat3R::from_ints([[0, -1, 0], [2, 0, 0], [0, 0, 1]]);
    let shift = [Rat::ZERO, Rat::new(1, 4).unwrap(), Rat::ZERO];
    assert!(m::realises_record_lattice_checked(&geometry, &set, &other).unwrap());
    m::verify_embedding(&geometry, &set, other, shift).unwrap();
}

#[test]
fn the_candidate_tracks_the_supplied_lattice_not_a_record_specific_lookup() {
    let (mut geometry, set) = context(21, "T1", "P1");
    for row in &mut geometry.basis {
        for value in row {
            *value = value.checked_mul(Rat::from_integer(2)).unwrap();
        }
    }
    geometry.index *= 8;
    let entries = m::embed_via_primitive_lattices(&geometry, &set).unwrap();
    let doubled = Mat3R::from_ints([[0, -2, 0], [2, 0, 0], [0, 0, 4]]);
    assert!(
        entries
            .iter()
            .any(|e| e.map == doubled && e.shift == [Rat::ZERO; 3])
    );
    assert!(m::realises_record_lattice_checked(&geometry, &set, &doubled).unwrap());
}

#[test]
fn arithmetic_limits_are_errors_not_negative_certificates() {
    let (mut geometry, set) = context(21, "T1", "P1");
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
        m::embed_via_primitive_lattices(&geometry, &set),
        Err(GeometryError::Arithmetic)
    );
}

#[test]
fn nontrivial_parent_actions_reconstruct_for_both_labels_and_parities() {
    use cryspglib::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
    };
    use cryspglib::irrep::magnetic_subduction_api::{
        MagneticRecordSelector, MagneticSubductionRequest, subduce,
    };
    use num_complex::Complex64;
    let (geometry, set) = context(21, "T1", "P1");
    let entry = m::embed_via_primitive_lattices(&geometry, &set)
        .unwrap()
        .into_iter()
        .find(|e| e.shift == [Rat::ZERO; 3])
        .unwrap();
    let transform = SeitzTransform::new(entry.map, Vec3R::new(entry.shift));
    let condensing = query::irreps_of(21).iter().find(|r| r.ml == "T1").unwrap();
    let (mut cases, mut negative_unitary, mut anti) = (0, 0, 0);
    for source in query::irreps_of(21)
        .iter()
        .filter(|r| ["GM1", "GM2", "GM3", "GM4"].contains(&r.ml))
    {
        let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
        for eta in [RealSign::Plus, RealSign::Minus] {
            for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                let source_label = source.label(convention).unwrap();
                let record_label = condensing.label(convention).unwrap();
                let request = MagneticSubductionRequest::new(
                    SourceIdentity {
                        parent_sg: 21,
                        label: &source_label,
                        convention,
                    },
                    MagneticRecordSelector::new(21, &record_label, convention)
                        .direction(IsotropyDirection::Label("P1"))
                        .uni(150),
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
                assert_eq!(model.record_label(), "T1");
                assert_eq!(model.source_label(), source.ml);
                assert_eq!(
                    (report.source_dimension(), report.covered_dimension()),
                    (1, 1)
                );
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
    assert_eq!(cases, 16);
    assert!(negative_unitary > 0 && anti > 0);
}

#[test]
fn the_complete_record_census_has_no_silent_source_lattice_failures() {
    let (records, verified) = (1u8..=230)
        .into_par_iter()
        .map(|sg| {
            let (mut records, mut verified) = (0, 0);
            for irrep in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
                for row in
                    isotropy::magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml)
                        .unwrap()
                {
                    records += 1;
                    let geometry = m::geometry_of(sg, &row.record).unwrap();
                    let set = m::magnetic_operations(row.record.mag_sg).unwrap();
                    let entries = m::embed_via_primitive_lattices(&geometry, &set).unwrap();
                    if let Some(entry) = entries.first() {
                        verified += 1;
                        assert!(entry.realises_record_lattice);
                        m::verify_embedding(&geometry, &set, entry.map, entry.shift).unwrap();
                    }
                }
            }
            (records, verified)
        })
        .reduce(|| (0, 0), |(a, b), (c, d)| (a + c, b + d));
    println!("primitive record census: records={records}, verified={verified}");
    assert_eq!(records, 16721);
    assert_eq!(verified, 11766);
}
