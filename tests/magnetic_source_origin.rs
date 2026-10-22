//! Origin-only alignment of Γ character rows; not a change of basis or a non-Γ phase rule.
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaRequest, SourceIdentity,
    verify_scalar_gamma_input,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, isotropy, query};

fn request(label: &str) -> ScalarGammaRequest<'_> {
    ScalarGammaRequest::new(
        SourceIdentity {
            parent_sg: 48,
            label,
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: 48,
            condensing_label: "GM1+",
            convention: LabelConvention::Cdml,
            direction_index: 0,
        },
        Some(RealSign::Plus),
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::identity(),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::ContainmentOnly,
    })
}

#[test]
fn sg48_has_one_origin_shift_that_matches_the_entire_source_row() {
    let source = query::irreps_of(48)
        .iter()
        .find(|r| r.ml == "GM1+")
        .unwrap();
    let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
    let parent = query::symmetry_operations_of(48).unwrap();
    assert_eq!(row.operations().len(), parent.operations.len());
    let mut origins = Vec::new();
    // Independent finite witness search on the eighth-cell grid, not production.
    // SG48 is primitive: its translation lattice is exactly Z^3.
    let units = |value: f64| {
        let scaled = value * 24.0;
        assert!((scaled - scaled.round()).abs() < 1e-10);
        scaled.round() as i64
    };
    for x in 0..8 {
        for y in 0..8 {
            for z in 0..8 {
                let shift = [3 * x, 3 * y, 3 * z];
                if row.operations().iter().all(|operation| {
                    parent.operations.iter().any(|p| {
                        p.rotation.iter().flatten().copied().eq(operation.rotation)
                            && (0..3).all(|i| {
                                (units(operation.translation[i]) + shift[i]
                                    - (0..3)
                                        .map(|j| {
                                            i64::from(operation.rotation[3 * i + j]) * shift[j]
                                        })
                                        .sum::<i64>()
                                    - units(p.translation[i]))
                                .rem_euclid(24)
                                    == 0
                            })
                    })
                }) {
                    origins.push(shift);
                }
            }
        }
    }
    assert!(
        !origins.is_empty(),
        "No joint origin; a per-operation match would be invalid"
    );
    println!("SG48 source-to-parent origins in units of 1/24: {origins:?}");
}

#[test]
fn sg48_aligned_scalar_input_reaches_the_public_contract() {
    let rows = isotropy::magnetic_isotropy_subgroups(48, "GM1+", LabelConvention::Cdml).unwrap();
    println!("SG48 GM1+ direction 0: UNI {}", rows[0].record.mag_sg);
    verify_scalar_gamma_input(&request("GM1+")).expect("a jointly origin-aligned source is usable");
}

// Check the reported coordinate change independently in integer 1/24 units.
// Full conventional parent operation lists include every centring representative.
fn assert_origin(sg: u8, label: &str, shift: &[Rat; 3]) {
    let source = query::irreps_of(sg).iter().find(|r| r.ml == label).unwrap();
    let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
    let parent = query::symmetry_operations_of(sg).unwrap();
    let units = |v: f64| {
        assert!((24.0 * v - (24.0 * v).round()).abs() < 1e-10);
        (24.0 * v).round() as i128
    };
    let shift = shift.map(|v| {
        let scaled = v.checked_mul(Rat::from_integer(24)).unwrap();
        assert!(scaled.is_integer());
        scaled.numerator()
    });
    for op in row.operations() {
        assert!(
            parent
                .operations
                .iter()
                .any(|p| p.rotation.iter().flatten().copied().eq(op.rotation)
                    && (0..3).all(|i| (units(op.translation[i]) + shift[i]
                        - (0..3)
                            .map(|j| i128::from(op.rotation[3 * i + j]) * shift[j])
                            .sum::<i128>()
                        - units(p.translation[i]))
                    .rem_euclid(24)
                        == 0)),
            "SG {sg}, {label}: {op:?}"
        );
    }
}

#[test]
fn high_level_subduction_preserves_labels_parities_and_nontrivial_source_actions() {
    use cryspglib::irrep::magnetic_subduction_api::{
        MagneticRecordSelector, MagneticSubductionRequest, subduce,
    };
    use num_complex::Complex64;
    let mut cases = 0;
    let mut anti_actions = 0;
    let mut negative_actions = 0;
    for source in query::irreps_of(48)
        .iter()
        .filter(|r| !r.spinor && r.dim == 1 && r.k_vector().numerators == [0; 3])
    {
        let row = source.ordinary_scalar_selected_arm_block_trace().unwrap();
        for condensing in ["GM1+", "GM1-"] {
            let record = query::irreps_of(48)
                .iter()
                .find(|r| r.ml == condensing)
                .unwrap();
            let rows = isotropy::magnetic_isotropy_subgroups(48, condensing, LabelConvention::Cdml)
                .unwrap();
            for eta in [RealSign::Plus, RealSign::Minus] {
                for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                    let source_label = source.label(convention).unwrap();
                    let record_label = record.label(convention).unwrap();
                    let selection = MagneticRecordSelector::new(48, &record_label, convention)
                        .direction(isotropy::IsotropyDirection::Label("P1"))
                        .uni(rows[0].record.mag_sg);
                    let request = MagneticSubductionRequest::new(
                        SourceIdentity {
                            parent_sg: 48,
                            label: &source_label,
                            convention,
                        },
                        selection,
                        Some(eta),
                    )
                    .with_embedding(EmbeddingEvidence::identity());
                    let report = subduce(&request).unwrap_or_else(|error| {
                        panic!(
                            "{} / {condensing} {eta:?} {convention:?}: {error}",
                            source.ml
                        )
                    });
                    let model = report.scalar_model().unwrap();
                    assert_eq!(model.source_label(), source.ml);
                    assert_eq!(model.record_label(), condensing);
                    assert_eq!(model.uni(), rows[0].record.mag_sg);
                    assert_eq!(model.time_reversal(), Some(eta));
                    assert_eq!(*model.source_origin_shift(), [Rat::new(1, 4).unwrap(); 3]);
                    assert_eq!(*model.shift(), [Rat::ZERO; 3]);
                    assert_eq!(
                        (report.source_dimension(), report.covered_dimension()),
                        (1, 1)
                    );
                    assert_origin(48, source.ml, model.source_origin_shift());
                    let payload = report.scalar().unwrap();
                    for (index, action) in model.actions().iter().enumerate() {
                        let position = row
                            .operations()
                            .iter()
                            .position(|op| {
                                action
                                    .class
                                    .rotation
                                    .iter()
                                    .flatten()
                                    .copied()
                                    .eq(op.rotation)
                            })
                            .unwrap();
                        let expected = row.values()[position].re
                            * if action.class.time_reversal {
                                f64::from(eta.value())
                            } else {
                                1.0
                            };
                        assert_eq!(f64::from(action.value.value()), expected);
                        anti_actions += usize::from(action.class.time_reversal);
                        negative_actions += usize::from(action.value == RealSign::Minus);
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
    assert_eq!(cases, 64);
    assert!(anti_actions > 0 && negative_actions > 0);
}

#[test]
fn a_centred_cubic_source_keeps_its_verified_nonzero_origin() {
    for label in ["GM1+", "GM2+"] {
        let request = ScalarGammaRequest::new(
            SourceIdentity {
                parent_sg: 227,
                label,
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: 227,
                condensing_label: "GM1+",
                convention: LabelConvention::Cdml,
                direction_index: 0,
            },
            Some(RealSign::Plus),
        );
        let model = verify_scalar_gamma_input(&request).unwrap();
        assert_ne!(*model.source_origin_shift(), [Rat::ZERO; 3]);
        assert_origin(227, label, model.source_origin_shift());
        cryspglib::irrep::magnetic_subduction::subduce_scalar_gamma(&model).unwrap();
    }
}

#[test]
fn an_already_aligned_source_keeps_the_zero_origin() {
    let request = ScalarGammaRequest::new(
        SourceIdentity {
            parent_sg: 1,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: 1,
            condensing_label: "GM1",
            convention: LabelConvention::Cdml,
            direction_index: 0,
        },
        None,
    );
    assert_eq!(
        *verify_scalar_gamma_input(&request)
            .unwrap()
            .source_origin_shift(),
        [Rat::ZERO; 3]
    );
}
