//! Exact public boundaries for finite, graded point-extension proposals.
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
    let set = m::magnetic_operations(geometry.uni).unwrap();
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
fn checked(g: &m::MagneticGeometry, s: &m::MagneticOperationSet, entries: &[m::ParentEmbedding]) {
    for e in entries {
        assert_eq!(e.kind, EmbeddingKind::PointExtension);
        assert!(e.realises_record_lattice);
        assert!(m::realises_record_lattice_checked(g, s, &e.map).unwrap());
        m::verify_embedding(g, s, e.map, e.shift).unwrap();
    }
}

#[test]
fn centred_and_nonsymmorphic_sources_have_complete_exact_origin_components() {
    for (sg, label, uni) in [
        (79, "GM3GM4", 22),
        (87, "GM3+GM4+", 69),
        (88, "GM3+GM4+", 96),
    ] {
        let (g, s, h) = context(sg, label, 0);
        assert_eq!(s.uni, uni);
        assert!(
            m::embed_via_unitary_subgroup(&g, &s, &h)
                .unwrap()
                .is_empty()
        );
        assert!(m::embed_via_primitive_lattices(&g, &s).unwrap().is_empty());
        assert!(
            m::embed_via_antitranslations(&g, &s, &h)
                .unwrap()
                .embeddings
                .is_empty()
        );
        let report = m::embed_via_point_extensions(&g, &s, &h).unwrap();
        assert!(!report.embeddings.is_empty());
        assert!(report.unresolved.is_empty());
        checked(&g, &s, &report.embeddings);
        if sg == 88 {
            assert!(report.embeddings.iter().all(|e| e.shift != [Rat::ZERO; 3]));
        }
        for e in &report.embeddings {
            for shift in m::embedding_shift_components(&g, &s, e.map)
                .unwrap()
                .representatives
            {
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
fn the_original_antiunitary_screw_is_not_replaced_by_the_zero_translation_auxiliary() {
    let (g, s, h) = context(4, "GM2", 0);
    assert_eq!(s.uni, 16);
    let report = m::embed_via_point_extensions(&g, &s, &h).unwrap();
    assert!(!report.embeddings.is_empty());
    let e = &report.embeddings[0];
    checked(&g, &s, &report.embeddings);
    let mut fake = s.clone();
    for op in fake.operations.iter_mut().chain(&mut fake.classes) {
        op.translation = [Rat::ZERO; 3];
    }
    m::verify_group(&fake).unwrap();
    assert!(matches!(
        m::verify_embedding(&g, &fake, e.map, e.shift),
        Err(GeometryError::NotInParent { .. })
    ));
    assert!(m::embed_via_point_extensions(&g, &fake, &h).is_err());
}

#[test]
fn rebases_preserve_exact_witnesses_and_find_all_binary_rebases() {
    let (g, s, h) = context(79, "GM3GM4", 0);
    let base = m::embed_via_point_extensions(&g, &s, &h).unwrap();
    assert!(!base.embeddings.is_empty());
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
        let report = m::embed_via_point_extensions(&rebased, &s, &h).unwrap();
        // The retained map proves existence under EVERY rebase independently
        // of whether the floating identifier proposes that map again.
        checked(&rebased, &s, &base.embeddings);
        checked(&rebased, &s, &report.embeddings);
        assert!(!report.embeddings.is_empty(), "bits={bits}");
        assert!(
            report.unresolved.is_empty(),
            "bits={bits}: {:?}",
            report.unresolved
        );
        count += 1;
    }
    assert_eq!(count, 168);
}

#[test]
fn wrong_identity_h_and_source_provenance_remain_errors() {
    let (g, s, h) = context(79, "GM3GM4", 0);
    let mut bad = g.clone();
    bad.uni += 1;
    assert!(m::embed_via_point_extensions(&bad, &s, &h).is_err());
    let mut bad = g.clone();
    bad.parent_sg = 1;
    assert!(m::embed_via_point_extensions(&bad, &s, &h).is_err());
    let mut bad = s.clone();
    bad.operations.swap(0, 1);
    assert!(m::embed_via_point_extensions(&g, &bad, &h).is_err());
    let mut bad = s.clone();
    bad.unitary_classes += 1;
    assert!(m::embed_via_point_extensions(&g, &bad, &h).is_err());
    let rows = isotropy::isotropy_subgroups(79, "GM1", LabelConvention::Cdml).unwrap();
    let wrong_h = SubgroupEmbedding::from_isotropy_subgroup(&rows[0]).unwrap();
    assert_ne!(wrong_h.subgroup_sg(), h.subgroup_sg());
    assert!(m::embed_via_point_extensions(&g, &s, &wrong_h).is_err());
}

#[test]
fn arithmetic_limits_are_not_reported_as_empty_searches() {
    let (g, s, h) = context(79, "GM3GM4", 0);
    let mut bad = g.clone();
    bad.basis[0] = bad.basis[1];
    assert!(matches!(
        m::embed_via_point_extensions(&bad, &s, &h),
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
        m::embed_via_point_extensions(&bad, &s, &h),
        Err(GeometryError::Arithmetic)
    ));
    let mut bad = g;
    bad.basis[0][1] = Rat::new(1, 2).unwrap();
    assert!(
        m::embed_via_point_extensions(&bad, &s, &h)
            .unwrap()
            .embeddings
            .is_empty()
    );
}

#[test]
fn inapplicable_colours_missing_parent_rotations_and_noninvariant_lattices_are_empty() {
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
            m::embed_via_point_extensions(&g, &m::magnetic_operations(uni).unwrap(), &ordinary)
                .unwrap();
        assert!(report.embeddings.is_empty());
        assert!(report.unresolved.is_empty());
    }
    let (mut g, s, _) = context(3, "GM2", 0);
    g.parent_sg = 1;
    assert!(
        m::embed_via_point_extensions(&g, &s, &ordinary)
            .unwrap()
            .embeddings
            .is_empty()
    );
    let (mut g, s, h) = context(10, "GM1-", 0);
    assert_eq!(h.subgroup_sg(), 3);
    g.basis = *Mat3R::from_ints([[3, 0, 0], [1, 1, 0], [0, 0, 1]]).rows();
    assert!(
        m::embed_via_point_extensions(&g, &s, &h)
            .unwrap()
            .embeddings
            .is_empty()
    );
    // Empty Type-III output cannot invalidate an independently proven Type-IV map.
    let (g, s, h) = context(3, "B1", 0);
    assert!(
        m::embed_via_point_extensions(&g, &s, &h)
            .unwrap()
            .embeddings
            .is_empty()
    );
    assert!(
        !m::embed_via_antitranslations(&g, &s, &h)
            .unwrap()
            .embeddings
            .is_empty()
    );
}

#[test]
fn two_sources_both_eta_and_both_label_conventions_reconstruct_every_action() {
    let (mut cases, mut negative_unitary, mut anti) = (0, 0, 0);
    for sg in [87, 88] {
        let (g, s, h) = context(sg, "GM3+GM4+", 0);
        let direction =
            isotropy::magnetic_isotropy_subgroups(sg, "GM3+GM4+", LabelConvention::Cdml).unwrap()
                [0]
            .record
            .direction;
        let report = m::embed_via_point_extensions(&g, &s, &h).unwrap();
        let entry = &report.embeddings[0];
        let transform = SeitzTransform::new(entry.map, Vec3R::new(entry.shift));
        let condensing = query::irreps_of(sg)
            .iter()
            .find(|r| r.ml == "GM3+GM4+")
            .unwrap();
        for source in query::irreps_of(sg)
            .iter()
            .filter(|r| ["GM1+", "GM1-"].contains(&r.ml))
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
                            .direction(IsotropyDirection::Label(direction))
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
fn whole_record_inventory_verifies_every_returned_witness_and_retains_uncertainty() {
    let counts = (1u8..=230)
        .into_par_iter()
        .map(|sg| {
            let mut counts = [0usize; 6];
            for irrep in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
                let ordinary =
                    isotropy::isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml).unwrap();
                for (index, row) in
                    isotropy::magnetic_isotropy_subgroups(sg, irrep.ml, LabelConvention::Cdml)
                        .unwrap()
                        .iter()
                        .enumerate()
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
                    match m::embed_via_point_extensions(&g, &s, &h) {
                        Ok(report) => {
                            checked(&g, &s, &report.embeddings);
                            counts[5] += report.unresolved.len();
                            counts[if report.embeddings.is_empty() { 2 } else { 3 }] += 1;
                            if !report.unresolved.is_empty() {
                                println!(
                                    "UNRESOLVED SG{sg} {} #{index} UNI{} {:?}",
                                    irrep.ml, s.uni, report.unresolved
                                );
                            }
                        }
                        Err(GeometryError::PointExtensionUnavailable { reason, .. }) => {
                            println!(
                                "UNAVAILABLE SG{sg} {} #{index} UNI{} {reason}",
                                irrep.ml, s.uni
                            );
                            counts[4] += 1;
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
        "point extension full inventory [records,no H,empty,positive,unavailable,unresolved]={counts:?}"
    );
    assert_eq!(counts, [16721, 12, 13837, 2872, 0, 0]);
}
