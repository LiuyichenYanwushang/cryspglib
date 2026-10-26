//! The magnetic raw primitive axes are not the data-Hall conventional axes.
use cryspglib::irrep::magnetic_embedding as m;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, isotropy, query};

fn context(sg: u8, label: &str, index: usize) -> (m::MagneticGeometry, m::MagneticOperationSet) {
    let row = isotropy::magnetic_isotropy_subgroups(sg, label, LabelConvention::Cdml).unwrap()
        [index]
        .record;
    (
        m::geometry_of(sg, &row).unwrap(),
        m::magnetic_operations(row.mag_sg).unwrap(),
    )
}

#[test]
fn monoclinic_parent_frames_match_official_display_parent_including_centering() {
    let half = Rat::new(1, 2).unwrap();
    for sg in 3..=15 {
        let row = query::magnetic_isotropy_subgroups_of(sg)[0].subgroup;
        let mut g = m::geometry_of(sg, &row).unwrap();
        g.basis = *Mat3R::identity().rows();
        g.origin = [
            Rat::new(1, 7).unwrap(),
            Rat::new(2, 7).unwrap(),
            Rat::new(3, 7).unwrap(),
        ];
        // Independently pinned by DISPLAY PARENT / SHOW BASIS, not inferred from
        // magnetic identities or the production centering helper.
        let expected = if [5, 8, 9, 12, 15].contains(&sg) {
            [
                [Rat::ZERO, Rat::ZERO, Rat::ONE],
                [half, half, Rat::ZERO],
                [Rat::new(-1, 2).unwrap(), half, Rat::ZERO],
            ]
        } else {
            *Mat3R::from_ints([[0, 0, 1], [1, 0, 0], [0, 1, 0]]).rows()
        };
        assert_eq!(
            m::basis_in_parent_conventional(&g).unwrap(),
            expected,
            "SG{sg}"
        );
        let expected_origin = if [5, 8, 9, 12, 15].contains(&sg) {
            [
                Rat::new(-1, 14).unwrap(),
                Rat::new(5, 14).unwrap(),
                Rat::new(1, 7).unwrap(),
            ]
        } else {
            [
                Rat::new(2, 7).unwrap(),
                Rat::new(3, 7).unwrap(),
                Rat::new(1, 7).unwrap(),
            ]
        };
        assert_eq!(
            m::origin_in_parent_conventional(&g).unwrap(),
            expected_origin,
            "SG{sg}"
        );
    }
}

#[test]
fn real_raw_rows_keep_provenance_but_convert_to_the_official_axes() {
    let (g, set) = context(3, "Y1", 0);
    assert_eq!(
        g.basis,
        *Mat3R::from_ints([[1, 0, 0], [0, 2, 0], [0, 0, 1]]).rows()
    );
    assert_eq!(
        m::basis_in_parent_conventional(&g).unwrap(),
        *Mat3R::from_ints([[0, 0, 1], [2, 0, 0], [0, 1, 0]]).rows()
    );
    assert_eq!(g.index, 2);
    assert!(m::record_lattice_obstruction(&g, &set).unwrap().is_none());
    let (a2, _) = context(3, "A2", 0);
    assert_eq!(
        m::origin_in_parent_conventional(&a2).unwrap(),
        [Rat::new(1, 2).unwrap(), Rat::ZERO, Rat::ZERO]
    );
    let rendered =
        isotropy::format_magnetic_isotropy_subgroups(3, "Y1", LabelConvention::Cdml).unwrap();
    assert!(rendered.contains("(0,0,1),(2,0,0),(0,1,0)"), "{rendered}");
}

#[test]
fn formerly_obstructed_records_have_exact_witnesses_in_the_correct_frame() {
    for (sg, label, index) in [
        (3, "Y1", 0),
        (3, "A1", 0),
        (3, "D1", 0),
        (7, "C1", 0),
        (13, "C1+", 0),
        (14, "C1", 1),
    ] {
        let (g, set) = context(sg, label, index);
        let entries = m::embed_via_oriented_primitive_lattices(&g, &set).unwrap();
        assert!(!entries.is_empty(), "SG{sg} {label}");
        for e in entries {
            assert!(m::realises_record_lattice_checked(&g, &set, &e.map).unwrap());
            m::verify_embedding(&g, &set, e.map, e.shift).unwrap();
        }
    }
}

#[test]
fn corrected_maps_reconstruct_high_level_gamma_actions_in_both_label_conventions() {
    use cryspglib::irrep::isotropy::IsotropyDirection;
    use cryspglib::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
    };
    use cryspglib::irrep::magnetic_subduction_api::{
        MagneticRecordSelector, MagneticSubductionRequest, subduce,
    };
    use cryspglib::irrep::subduction::{SeitzTransform, Vec3R};
    use num_complex::Complex64;
    let (mut cases, mut negative, mut anti) = (0, 0, 0);
    for (sg, label, index) in [
        (3, "Y1", 0),
        (3, "A2", 0),
        (5, "Y1", 0),
        (7, "C1", 0),
        (13, "C1+", 0),
        (14, "C1", 1),
    ] {
        let (g, set) = context(sg, label, index);
        let entry = m::embed_via_oriented_primitive_lattices(&g, &set)
            .unwrap()
            .remove(0);
        let transform = SeitzTransform::new(entry.map, Vec3R::new(entry.shift));
        let condensing = query::irreps_of(sg).iter().find(|r| r.ml == label).unwrap();
        for source in query::irreps_of(sg)
            .iter()
            .filter(|r| !r.spinor && r.ml.starts_with("GM"))
        {
            let trace = source.ordinary_scalar_selected_arm_block_trace().unwrap();
            for eta in [RealSign::Plus, RealSign::Minus] {
                for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                    let source_label = source.label(convention).unwrap();
                    let record_label = condensing.label(convention).unwrap();
                    let request = MagneticSubductionRequest::new(
                        SourceIdentity {
                            parent_sg: sg,
                            label: &source_label,
                            convention,
                        },
                        MagneticRecordSelector::new(sg, &record_label, convention)
                            .direction(IsotropyDirection::Index(index))
                            .uni(g.uni),
                        Some(eta),
                    )
                    .with_embedding(EmbeddingEvidence {
                        map: entry.map,
                        shift: entry.shift,
                        lattice_policy: LatticePolicy::RequireExact,
                    });
                    let report = subduce(&request).unwrap();
                    assert_eq!(
                        (report.source_dimension(), report.covered_dimension()),
                        (1, 1)
                    );
                    let model = report.scalar_model().unwrap();
                    let result = report.scalar().unwrap();
                    for (i, action) in model.actions().iter().enumerate() {
                        let rotation = transform.map_rotation(action.class.rotation).unwrap();
                        let j = trace
                            .operations()
                            .iter()
                            .position(|op| rotation.iter().flatten().copied().eq(op.rotation))
                            .unwrap();
                        let expected = trace.values()[j].re
                            * if action.class.time_reversal {
                                f64::from(eta.value())
                            } else {
                                1.0
                            };
                        assert_eq!(f64::from(action.value.value()), expected);
                        negative += usize::from(!action.class.time_reversal && expected == -1.0);
                        anti += usize::from(action.class.time_reversal);
                        for z in [Complex64::new(1.0, 2.0), Complex64::new(-0.7, 0.4)] {
                            let v = result.intertwiner * z;
                            let lhs = expected
                                * if action.class.time_reversal {
                                    v.conj()
                                } else {
                                    v
                                };
                            assert!(
                                (lhs - result.intertwiner * result.operators[i].apply(z)).norm()
                                    < 1e-12
                            );
                        }
                    }
                    cases += 1;
                }
            }
        }
    }
    assert_eq!(cases, 64);
    assert!(negative > 0 && anti > 0);
}
