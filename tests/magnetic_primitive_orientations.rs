//! Exact lattice-aware setting candidates outside the R24 finite map family.
use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy;
use cryspglib::irrep::magnetic_embedding::{self as m, EmbeddingKind};
use cryspglib::irrep::subduction::{Mat3R, Rat};

fn context(sg: u8, label: &str, index: usize) -> (m::MagneticGeometry, m::MagneticOperationSet) {
    let rows = isotropy::magnetic_isotropy_subgroups(sg, label, LabelConvention::Cdml).unwrap();
    let record = &rows[index].record;
    (
        m::geometry_of(sg, record).unwrap(),
        m::magnetic_operations(record.mag_sg).unwrap(),
    )
}

#[test]
fn an_oriented_primitive_cell_finds_a_nonzero_origin_without_replacing_older_routes() {
    let (geometry, set) = context(85, "GM3+GM4+", 0);
    assert_eq!(set.uni, 76);
    assert!(
        m::parent_embedding_verdict(&geometry, &set)
            .unwrap()
            .embeddings
            .iter()
            .all(|e| !e.realises_record_lattice)
    );
    assert!(
        m::embed_via_primitive_lattices(&geometry, &set)
            .unwrap()
            .is_empty()
    );
    assert!(
        m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .is_none()
    );
    let map = Mat3R::from_ints([[-1, 0, 1], [0, 0, 1], [0, 1, 0]]);
    let shift = [Rat::new(1, 4).unwrap(), Rat::new(1, 4).unwrap(), Rat::ZERO];
    assert!(m::realises_record_lattice_checked(&geometry, &set, &map).unwrap());
    assert!(m::verify_embedding(&geometry, &set, map, [Rat::ZERO; 3]).is_err());
    m::verify_embedding(&geometry, &set, map, shift).unwrap();
    let entries = m::embed_via_oriented_primitive_lattices(&geometry, &set).unwrap();
    assert!(entries.iter().any(|e| e.map == map && e.shift == shift));
    for e in entries {
        assert_eq!(e.kind, EmbeddingKind::PrimitiveOrientation);
        assert!(e.realises_record_lattice);
        assert!(m::realises_record_lattice_checked(&geometry, &set, &e.map).unwrap());
        m::verify_embedding(&geometry, &set, e.map, e.shift).unwrap();
    }
}

#[test]
fn centred_and_nonprimitive_source_settings_keep_exact_full_lattice() {
    for (sg, label, index, uni) in [(185, "GM5", 1, 274), (211, "N3", 2, 24)] {
        let (geometry, set) = context(sg, label, index);
        assert_eq!(set.uni, uni);
        assert!(
            m::embed_via_primitive_lattices(&geometry, &set)
                .unwrap()
                .is_empty()
        );
        let entries = m::embed_via_oriented_primitive_lattices(&geometry, &set).unwrap();
        assert!(!entries.is_empty());
        for e in entries {
            assert!(m::realises_record_lattice_checked(&geometry, &set, &e.map).unwrap());
            m::verify_embedding(&geometry, &set, e.map, e.shift).unwrap();
        }
    }
}

#[test]
fn recovered_gamma_actions_reconstruct_for_both_parities_and_label_conventions() {
    use cryspglib::irrep::isotropy::IsotropyDirection;
    use cryspglib::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
    };
    use cryspglib::irrep::magnetic_subduction_api::{
        MagneticRecordSelector, MagneticSubductionRequest, subduce,
    };
    use cryspglib::irrep::query;
    use cryspglib::irrep::subduction::{SeitzTransform, Vec3R};
    use num_complex::Complex64;
    let (geometry, set) = context(85, "GM3+GM4+", 0);
    let entry = m::embed_via_oriented_primitive_lattices(&geometry, &set)
        .unwrap()
        .into_iter()
        .find(|e| {
            e.map == Mat3R::from_ints([[-1, 0, 1], [0, 0, 1], [0, 1, 0]])
                && e.shift == [Rat::new(1, 4).unwrap(), Rat::new(1, 4).unwrap(), Rat::ZERO]
        })
        .unwrap();
    let transform = SeitzTransform::new(entry.map, Vec3R::new(entry.shift));
    let condensing = query::irreps_of(85)
        .iter()
        .find(|r| r.ml == "GM3+GM4+")
        .unwrap();
    let mut cases = 0;
    let mut negative_unitary = 0;
    let mut anti = 0;
    for source in query::irreps_of(85)
        .iter()
        .filter(|r| ["GM1+", "GM2+", "GM1-", "GM2-"].contains(&r.ml))
    {
        let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
        for eta in [RealSign::Plus, RealSign::Minus] {
            for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                let source_label = source.label(convention).unwrap();
                let record_label = condensing.label(convention).unwrap();
                let request = MagneticSubductionRequest::new(
                    SourceIdentity {
                        parent_sg: 85,
                        label: &source_label,
                        convention,
                    },
                    MagneticRecordSelector::new(85, &record_label, convention)
                        .direction(IsotropyDirection::Index(0))
                        .uni(76),
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
                assert_eq!(
                    (report.source_dimension(), report.covered_dimension()),
                    (1, 1)
                );
                for (index, action) in model.actions().iter().enumerate() {
                    let mapped = transform.map_rotation(action.class.rotation).unwrap();
                    let position = row
                        .operations()
                        .iter()
                        .position(|op| mapped.iter().flatten().copied().eq(op.rotation))
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
fn literal_lattice_obstruction_is_not_overridden_by_new_orientation_proposals() {
    let (geometry, set) = context(3, "A1", 0);
    assert!(
        m::record_lattice_obstruction(&geometry, &set)
            .unwrap()
            .is_some()
    );
    assert!(
        m::embed_via_oriented_primitive_lattices(&geometry, &set)
            .unwrap()
            .is_empty()
    );
    let mut corrupted = set.clone();
    corrupted.operations.swap(0, 1);
    assert!(m::embed_via_oriented_primitive_lattices(&geometry, &corrupted).is_err());
}
