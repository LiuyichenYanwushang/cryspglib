//! Independent regressions: repeated constraints, separate identities and complex gauge.
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::{canonical_translation, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaModel,
    ScalarGammaRequest, SourceIdentity, verify_scalar_gamma_input,
};
use cryspglib::irrep::magnetic_subduction::{MagneticSubduction, subduce_scalar_gamma};
use cryspglib::irrep::subduction::{Mat3R, Rat, subduce_irrep};
use cryspglib::irrep::{LabelConvention, isotropy};
use num_complex::Complex64;

fn model(
    sg: u8,
    source: &str,
    condensing: &str,
    uni: usize,
    eta: Option<RealSign>,
) -> ScalarGammaModel {
    let rows =
        isotropy::magnetic_isotropy_subgroups(sg, condensing, LabelConvention::Cdml).unwrap();
    let index = rows
        .iter()
        .position(|row| row.record.mag_sg == uni)
        .unwrap();
    verify_scalar_gamma_input(
        &ScalarGammaRequest::new(
            SourceIdentity {
                parent_sg: sg,
                label: source,
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: sg,
                condensing_label: condensing,
                convention: LabelConvention::Cdml,
                direction_index: index,
            },
            eta,
        )
        .with_embedding(EmbeddingEvidence {
            map: Mat3R::identity(),
            shift: [Rat::ZERO; 3],
            lattice_policy: LatticePolicy::ContainmentOnly,
        }),
    )
    .unwrap()
}

fn check_intertwiner(model: &ScalarGammaModel, report: &MagneticSubduction<'_>) {
    assert_eq!(report.operators().len(), model.classes().len());
    assert_eq!(report.covered_dimension(), 1);
    let u = report.intertwiner();
    assert!((u.norm_sqr() - 1.0).abs() < 1e-12);
    // Check actual anti/linear maps on complex vectors, not just real factors.
    for z in [Complex64::new(1.0, 2.0), Complex64::new(-3.0, 0.5)] {
        for (index, action) in model.actions().iter().enumerate() {
            let source_coordinate = u * z;
            let left = Complex64::new(f64::from(action.value.value()), 0.0)
                * if action.class.time_reversal {
                    source_coordinate.conj()
                } else {
                    source_coordinate
                };
            let right = u * report.operators()[index].apply(z);
            assert!(
                (left - right).norm() < 1e-12,
                "SG {} record {} UNI {} class {index}: U={u:?}, source={left:?}, target={right:?}",
                model.parent_sg(),
                model.record_label(),
                model.uni()
            );
        }
    }
}

#[test]
fn repeated_antiunitary_constraints_do_not_lose_the_intertwiner() {
    // H has two classes; both antiunitary classes impose the same equation.
    // The first draft mistook that repetition for contradiction and silently
    // returned U=1, failing the eta=-1 reconstruction.
    let plus = model(10, "GM1+", "GM2+", 53, Some(RealSign::Plus));
    let minus = model(10, "GM1+", "GM2+", 53, Some(RealSign::Minus));
    assert_eq!(
        minus.classes().iter().filter(|c| c.time_reversal).count(),
        2
    );
    let plus_result = subduce_scalar_gamma(&plus).unwrap();
    let minus_result = subduce_scalar_gamma(&minus).unwrap();
    assert_eq!(plus_result.target_label(), minus_result.target_label());
    assert_eq!(plus_result.operators(), minus_result.operators());
    assert_eq!(plus_result.intertwiner(), Complex64::new(1.0, 0.0));
    assert_eq!(minus_result.intertwiner(), Complex64::new(0.0, 1.0));
    check_intertwiner(&plus, &plus_result);
    check_intertwiner(&minus, &minus_result);
}

#[test]
fn type_one_nontrivial_source_agrees_with_the_ordinary_engine() {
    let model = model(2, "GM1-", "GM1+", 4, None);
    let report = subduce_scalar_gamma(&model).unwrap();
    assert_eq!(report.target_source_label(), "GM1-");
    assert_ne!(report.target_source_label(), model.record_label());
    assert!(report.operators().iter().any(|op| op.value().re == -1.0));
    let ordinary_subgroup = isotropy::isotropy_subgroup_for_direction(
        2,
        "GM1+",
        LabelConvention::Cdml,
        IsotropyDirection::Label("P1"),
    )
    .unwrap();
    let ordinary = subduce_irrep(&ordinary_subgroup, "GM1-").unwrap();
    let [term] = ordinary.targets() else {
        panic!("exactly one ordinary target");
    };
    assert_eq!(term.ml, Some(report.target_source_label()));
    assert_eq!(term.sg, report.target_unitary_sg());
    assert_eq!(u32::from(term.dimension), report.target_dimension());
    assert_eq!(term.multiplicity, report.multiplicity());
    check_intertwiner(&model, &report);
}

#[test]
fn lawful_complex_gauge_preserves_all_operator_products_and_reconstruction() {
    let model = model(10, "GM1-", "GM2+", 53, Some(RealSign::Minus));
    let report = subduce_scalar_gamma(&model).unwrap();
    check_intertwiner(&model, &report);
    let q = Complex64::new(0.6, 0.8);
    let u = report.intertwiner() * q;
    let coefficients: Vec<_> = report
        .operators()
        .iter()
        .map(|op| {
            op.value()
                * if op.is_antiunitary() {
                    q.conj() / q
                } else {
                    Complex64::new(1.0, 0.0)
                }
        })
        .collect();
    let z = Complex64::new(1.0, 2.0);
    for (i, class) in model.classes().iter().enumerate() {
        let source_coordinate = u * z;
        let source = f64::from(model.actions()[i].value.value())
            * if class.time_reversal {
                source_coordinate.conj()
            } else {
                source_coordinate
            };
        let target = u * coefficients[i] * if class.time_reversal { z.conj() } else { z };
        assert!((source - target).norm() < 1e-12);
        for (j, right) in model.classes().iter().enumerate() {
            let product = compose(class, right).unwrap();
            let k = model
                .classes()
                .iter()
                .position(|candidate| {
                    candidate.rotation == product.rotation
                        && candidate.time_reversal == product.time_reversal
                        && canonical_translation(candidate.translation, model.translation_lattice())
                            == canonical_translation(
                                product.translation,
                                model.translation_lattice(),
                            )
                })
                .unwrap();
            let found = coefficients[i]
                * if class.time_reversal {
                    coefficients[j].conj()
                } else {
                    coefficients[j]
                };
            assert!((found - coefficients[k]).norm() < 1e-12);
        }
    }
    // Keeping this complex gauge but the old U breaks the fixed intertwiner,
    // although the gauged operators alone still obey the group law.
    let a = model
        .classes()
        .iter()
        .position(|c| c.time_reversal)
        .unwrap();
    let wrong = report.intertwiner();
    let left = f64::from(model.actions()[a].value.value()) * (wrong * z).conj();
    let right = wrong * coefficients[a] * z.conj();
    assert!((left - right).norm() > 1e-3);
}
