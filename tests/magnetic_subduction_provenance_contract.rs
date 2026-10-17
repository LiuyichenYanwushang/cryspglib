//! Independent raw-index/frame and complex reconstruction regressions.
use cryspglib::SymmetryOps;
use cryspglib::irrep::magnetic_embedding::{canonical_translation, compose, magnetic_operations};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaRequest, SourceIdentity,
    verify_scalar_gamma_input,
};
use cryspglib::irrep::magnetic_subduction::subduce_scalar_gamma;
use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_uni_partial;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, isotropy};
use num_complex::Complex64;

fn check(sg: u8, source: &str, record: &str, uni: usize) {
    let raw = SymmetryOps::from_magnetic_database(uni).unwrap();
    let exact = magnetic_operations(uni).unwrap();
    assert_eq!(raw.operations.len(), exact.operations.len());
    for (a, b) in raw.operations.iter().zip(&exact.operations) {
        assert_eq!(a.rotation, b.rotation);
        assert_eq!(a.time_reversal, b.time_reversal);
        for d in 0..3 {
            assert_eq!(a.translation[d], b.translation[d].to_f64());
        }
    }
    let partial = magnetic_irrep_summary_by_uni_partial(uni).unwrap();
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|k| k.coords == (0, 0, 0, 1))
        .unwrap();
    assert_eq!(gamma.operations.len(), raw.operations.len());
    // These witnesses actually have different display and raw coordinates.
    assert!(gamma.operations.iter().any(|op| {
        let r = &raw.operations[op.magnetic_operation_index];
        op.rotation != r.rotation || op.translation != r.translation
    }));
    let rows = isotropy::magnetic_isotropy_subgroups(sg, record, LabelConvention::Cdml).unwrap();
    let index = rows.iter().position(|r| r.record.mag_sg == uni).unwrap();
    let mut identity = None;
    for eta in [RealSign::Plus, RealSign::Minus] {
        let model = verify_scalar_gamma_input(
            &ScalarGammaRequest::new(
                SourceIdentity {
                    parent_sg: sg,
                    label: source,
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: sg,
                    condensing_label: record,
                    convention: LabelConvention::Cdml,
                    direction_index: index,
                },
                Some(eta),
            )
            .with_embedding(EmbeddingEvidence {
                map: Mat3R::identity(),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::ContainmentOnly,
            }),
        )
        .unwrap();
        let result = subduce_scalar_gamma(&model).unwrap();
        assert_eq!(result.multiplicity(), 1);
        assert_eq!(result.covered_dimension(), 1);
        if source == "GM1-" {
            assert_ne!(result.target_source_label(), "GM1");
            assert!(
                result
                    .operators()
                    .iter()
                    .any(|op| !op.is_antiunitary() && op.value().re == -1.0)
            );
        }
        let current = (
            result.target_label().to_owned(),
            result.target_source_label(),
            result.target_unitary_sg(),
        );
        if let Some(previous) = identity {
            assert_eq!(current, previous);
        }
        identity = Some(current);
        let target = gamma
            .coreps
            .iter()
            .find(|c| c.label == result.target_label())
            .unwrap();
        assert_eq!(target.dim, 1);
        assert_eq!(target.source_irreps.len(), 1);
        assert_eq!(target.source_irreps[0].ml, result.target_source_label());
        assert_eq!(target.source_irreps[0].sg, result.target_unitary_sg());
        for op in &gamma.operations {
            let r = &exact.operations[op.magnetic_operation_index];
            let i = model
                .classes()
                .iter()
                .position(|c| {
                    c.rotation == r.rotation
                        && c.time_reversal == r.time_reversal
                        && canonical_translation(c.translation, model.translation_lattice())
                            == canonical_translation(r.translation, model.translation_lattice())
                })
                .unwrap();
            assert_eq!(result.operators()[i].is_antiunitary(), r.time_reversal);
            if !r.time_reversal {
                assert!(
                    (result.operators()[i].value() - target.characters[op.column].unwrap()).norm()
                        < 1e-12
                );
            }
        }
        let u = result.intertwiner();
        assert!((u.norm_sqr() - 1.0).abs() < 1e-12);
        for (i, g) in model.classes().iter().enumerate() {
            for z in [Complex64::new(1.0, 2.0), Complex64::new(-0.7, 0.4)] {
                let v = u * z;
                let a = f64::from(model.actions()[i].value.value())
                    * if g.time_reversal { v.conj() } else { v };
                assert!((a - u * result.operators()[i].apply(z)).norm() < 1e-12);
            }
            for (j, h) in model.classes().iter().enumerate() {
                let p = compose(g, h).unwrap();
                let k = model
                    .classes()
                    .iter()
                    .position(|c| {
                        c.rotation == p.rotation
                            && c.time_reversal == p.time_reversal
                            && canonical_translation(c.translation, model.translation_lattice())
                                == canonical_translation(p.translation, model.translation_lattice())
                    })
                    .unwrap();
                assert_eq!(
                    result.operators()[i].compose(&result.operators()[j]),
                    result.operators()[k]
                );
            }
        }
    }
}

#[test]
fn centered_uni22_uses_raw_indices_not_data_hall_coordinates() {
    check(5, "GM1", "GM2", 22);
}

#[test]
fn uni59_uses_raw_indices_and_reconstructs_both_eta_choices() {
    check(11, "GM1+", "GM2-", 59);
}

#[test]
fn uni59_nontrivial_source_keeps_an_independent_target_identity() {
    check(11, "GM1-", "GM2-", 59);
}
