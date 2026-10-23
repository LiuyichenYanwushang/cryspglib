//! Public proof boundaries for geometry-derived antitranslation proposals.
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::{self as m, EmbeddingKind, GeometryError};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
};
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticRecordSelector, MagneticSubductionRequest, subduce,
};
use cryspglib::irrep::subduction::{Mat3R, Rat, SeitzTransform, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::{LabelConvention, corep, isotropy, query};
use num_complex::Complex64;
use rayon::prelude::*;

fn context(
    sg: u8,
    label: &str,
    index: usize,
) -> (
    m::MagneticGeometry,
    m::MagneticOperationSet,
    SubgroupEmbedding,
) {
    let record =
        isotropy::magnetic_isotropy_subgroups(sg, label, LabelConvention::Cdml).unwrap()[index];
    let geometry = m::geometry_of(sg, &record.record).unwrap();
    let set = m::magnetic_operations(record.record.mag_sg).unwrap();
    let h = corep::identify_unitary_subgroup(set.uni).unwrap();
    let rows = isotropy::isotropy_subgroups(sg, label, LabelConvention::Cdml).unwrap();
    let row = rows
        .iter()
        .find(|r| r.record.sg == h && r.record.direction_label == record.record.direction)
        .or_else(|| rows.iter().find(|r| r.record.sg == h))
        .unwrap();
    (
        geometry,
        set,
        SubgroupEmbedding::from_isotropy_subgroup(row).unwrap(),
    )
}
fn checked(
    geometry: &m::MagneticGeometry,
    set: &m::MagneticOperationSet,
    entries: &[m::ParentEmbedding],
) {
    for e in entries {
        assert_eq!(e.kind, EmbeddingKind::Antitranslation);
        assert!(e.realises_record_lattice);
        assert!(m::realises_record_lattice_checked(geometry, set, &e.map).unwrap());
        m::verify_embedding(geometry, set, e.map, e.shift).unwrap();
    }
}

#[test]
fn primitive_centred_and_nonsymmorphic_sources_gain_exact_witnesses() {
    for (sg, label, uni) in [(3, "B1", 11), (3, "C1", 24), (4, "B1", 17)] {
        let (g, s, h) = context(sg, label, 0);
        assert_eq!(s.uni, uni);
        assert!(
            m::embed_via_unitary_subgroup(&g, &s, &h)
                .unwrap()
                .is_empty()
        );
        assert!(m::embed_via_primitive_lattices(&g, &s).unwrap().is_empty());
        let report = m::embed_via_antitranslations(&g, &s, &h).unwrap();
        assert!(!report.embeddings.is_empty());
        assert!(report.unresolved.is_empty());
        checked(&g, &s, &report.embeddings);
        for e in &report.embeddings {
            let shifts = m::embedding_shift_components(&g, &s, e.map).unwrap();
            for shift in shifts.representatives {
                assert!(
                    report
                        .embeddings
                        .iter()
                        .any(|other| other.map == e.map && other.shift == shift)
                );
            }
        }
    }
}

#[test]
fn original_screw_translations_cannot_be_replaced_by_the_auxiliary_group() {
    let (g, s, h) = context(4, "B1", 0);
    let report = m::embed_via_antitranslations(&g, &s, &h).unwrap();
    let e = &report.embeddings[0];
    let alpha = s
        .classes
        .iter()
        .find(|o| o.time_reversal && o.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        .unwrap()
        .translation;
    let mut virtual_group = s.clone();
    for o in virtual_group
        .operations
        .iter_mut()
        .chain(&mut virtual_group.classes)
    {
        o.translation = if o.time_reversal {
            alpha
        } else {
            [Rat::ZERO; 3]
        };
    }
    m::verify_group(&virtual_group).unwrap();
    assert!(matches!(
        m::verify_embedding(&g, &virtual_group, e.map, e.shift),
        Err(GeometryError::NotInParent { .. })
    ));
    m::verify_embedding(&g, &s, e.map, e.shift).unwrap();
}

#[test]
fn the_seven_cosets_are_recomputed_under_all_binary_unimodular_rebases() {
    let (g, s, h) = context(3, "B1", 0);
    let mut count = 0;
    for bits in 0..512 {
        let u = Mat3R::from_ints(std::array::from_fn(|i| {
            std::array::from_fn(|j| (bits >> (3 * i + j)) & 1)
        }));
        if u.determinant().unwrap().numerator().unsigned_abs() != 1 {
            continue;
        }
        let mut rebased = g.clone();
        rebased.basis = *Mat3R::new(g.basis)
            .transpose()
            .checked_mul(&u)
            .unwrap()
            .transpose()
            .rows();
        let report = m::embed_via_antitranslations(&rebased, &s, &h).unwrap();
        assert!(
            !report.embeddings.is_empty(),
            "bits={bits}, unresolved={:?}",
            report.unresolved
        );
        checked(&rebased, &s, &report.embeddings);
        count += 1;
    }
    assert_eq!(count, 168);
}

#[test]
fn malformed_identity_provenance_and_h_are_errors() {
    let (g, s, h) = context(3, "B1", 0);
    let mut bad = g.clone();
    bad.uni += 1;
    assert!(m::embed_via_antitranslations(&bad, &s, &h).is_err());
    let mut bad = g.clone();
    bad.parent_sg = 1;
    assert!(m::embed_via_antitranslations(&bad, &s, &h).is_err());
    let mut bad = s.clone();
    bad.operations.swap(0, 1);
    assert!(m::embed_via_antitranslations(&g, &bad, &h).is_err());
    let mut bad = s.clone();
    bad.unitary_classes += 1;
    assert!(m::embed_via_antitranslations(&g, &bad, &h).is_err());
    let (_, _, wrong_h) = context(3, "GM2", 0);
    assert_ne!(wrong_h.subgroup_sg(), h.subgroup_sg());
    assert!(m::embed_via_antitranslations(&g, &s, &wrong_h).is_err());
}

#[test]
fn arithmetic_limits_remain_errors_and_fractional_target_periods_are_rejected() {
    let (g, s, h) = context(3, "B1", 0);
    let mut bad = g.clone();
    bad.basis[0] = bad.basis[1];
    assert!(matches!(
        m::embed_via_antitranslations(&bad, &s, &h),
        Err(GeometryError::Arithmetic)
    ));
    let mut bad = g.clone();
    let huge = Rat::new(i128::MAX, 1).unwrap();
    bad.basis = [
        [huge, Rat::ZERO, Rat::ZERO],
        [Rat::ZERO, huge, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, huge],
    ];
    assert!(matches!(
        m::embed_via_antitranslations(&bad, &s, &h),
        Err(GeometryError::Arithmetic)
    ));
    let mut bad = g;
    bad.basis[0][1] = Rat::new(1, 2).unwrap();
    assert!(
        m::embed_via_antitranslations(&bad, &s, &h)
            .unwrap()
            .embeddings
            .is_empty()
    );
}

#[test]
fn inapplicable_colours_and_literal_lattice_obstructions_do_not_become_witnesses() {
    let rows = isotropy::isotropy_subgroups(1, "GM1", LabelConvention::Cdml).unwrap();
    let ordinary = SubgroupEmbedding::from_isotropy_subgroup(&rows[0]).unwrap();
    for uni in 1..=3 {
        let g = m::MagneticGeometry {
            parent_sg: 1,
            uni,
            basis: *Mat3R::identity().rows(),
            origin: [Rat::ZERO; 3],
            index: 1,
        };
        let report =
            m::embed_via_antitranslations(&g, &m::magnetic_operations(uni).unwrap(), &ordinary)
                .unwrap();
        assert!(report.embeddings.is_empty());
        assert!(report.unresolved.is_empty());
    }
    let (g, s, h) = context(3, "GM2", 0);
    assert!(
        m::embed_via_antitranslations(&g, &s, &h)
            .unwrap()
            .embeddings
            .is_empty()
    );
    for label in ["A1", "D1"] {
        let (g, s, h) = context(3, label, 0);
        assert!(m::record_lattice_obstruction(&g, &s).unwrap().is_some());
        assert!(
            m::embed_via_antitranslations(&g, &s, &h)
                .unwrap()
                .embeddings
                .is_empty()
        );
    }
}

#[test]
fn identifier_failure_is_preserved_as_an_unresolved_proposal() {
    let (g, s, h) = context(20, "T1", 0);
    let report = m::embed_via_antitranslations(&g, &s, &h).unwrap();
    assert!(report.embeddings.is_empty());
    assert!(report.unresolved.iter().any(|f| f.translation
        == [Rat::new(1, 2).unwrap(), Rat::ZERO, Rat::ZERO]
        && f.error == cryspglib::SymError::MagneticReferenceGroupFailed));
}

#[test]
fn high_level_sources_eta_and_labels_reconstruct_the_original_actions() {
    let (mut cases, mut negative_unitary, mut anti) = (0, 0, 0);
    for sg in [3, 4] {
        let (g, s, h) = context(sg, "B1", 0);
        let report = m::embed_via_antitranslations(&g, &s, &h).unwrap();
        let entry = report
            .embeddings
            .iter()
            .find(|e| e.shift == [Rat::ZERO; 3])
            .unwrap();
        let transform = SeitzTransform::new(entry.map, Vec3R::new(entry.shift));
        let condensing = query::irreps_of(sg).iter().find(|r| r.ml == "B1").unwrap();
        for source in query::irreps_of(sg)
            .iter()
            .filter(|r| ["GM1", "GM2"].contains(&r.ml))
        {
            let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
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
                            .direction(IsotropyDirection::Label("P1"))
                            .uni(s.uni),
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
                    let payload = report.scalar().unwrap();
                    for (index, action) in model.actions().iter().enumerate() {
                        let r = transform.map_rotation(action.class.rotation).unwrap();
                        let position = row
                            .operations()
                            .iter()
                            .position(|op| r.iter().flatten().copied().eq(op.rotation))
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
    }
    assert_eq!(cases, 16);
    assert!(negative_unitary > 0 && anti > 0);
}

#[test]
fn whole_record_inventory_separates_witnesses_and_proposal_uncertainty() {
    let counts = (1u8..=230)
        .into_par_iter()
        .map(|sg| {
            let mut counts = [0usize; 6];
            for irrep in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
                let ordinary =
                    isotropy::isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml).unwrap();
                for row in
                    isotropy::magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml)
                        .unwrap()
                {
                    counts[0] += 1;
                    let g = m::geometry_of(sg, &row.record).unwrap();
                    let s = m::magnetic_operations(g.uni).unwrap();
                    let h = corep::identify_unitary_subgroup(s.uni).unwrap();
                    let Some(row) = ordinary.iter().find(|r| r.record.sg == h) else {
                        counts[1] += 1;
                        continue;
                    };
                    let h = SubgroupEmbedding::from_isotropy_subgroup(row).unwrap();
                    match m::embed_via_antitranslations(&g, &s, &h) {
                        Ok(report) => {
                            checked(&g, &s, &report.embeddings);
                            counts[5] += report.unresolved.len();
                            counts[if report.embeddings.is_empty() { 2 } else { 3 }] += 1;
                        }
                        Err(GeometryError::AntitranslationBridgeUnavailable { .. }) => {
                            counts[4] += 1
                        }
                        Err(e) => panic!("SG{sg} {} UNI{}: {e}", irrep.ml, s.uni),
                    }
                }
            }
            counts
        })
        .reduce(
            || [0; 6],
            |mut a, b| {
                for i in 0..6 {
                    a[i] += b[i];
                }
                a
            },
        );
    println!(
        "antitranslation full inventory [records,no H,empty,positive,unavailable,unresolved]={counts:?}"
    );
    assert_eq!(counts, [16721, 12, 5184, 11525, 0, 130]);
}
