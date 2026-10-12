//! R26 regressions: the restricted scalar Γ magnetic subduction.
//!
//! Every test here re-derives its claim from the public API: the operator table is
//! re-checked against the actual multiplication of the magnetic group, the
//! intertwiner contract is re-applied class by class, and the target identity is
//! re-read from the magnetic corep summary.  The module's own checks are not taken
//! on faith.
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::MagneticOperation;
use cryspglib::irrep::magnetic_embedding::{canonical_translation, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaModel, ScalarGammaRequest,
    SourceIdentity, verify_scalar_gamma_input,
};
use cryspglib::irrep::magnetic_subduction::{
    MagneticSubduction, MagneticSubductionError, subduce_scalar_gamma,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::subduction_api::{SubductionRequest, subduce};
use cryspglib::irrep::{LabelConvention, isotropy, magnetic_summary};
use num_complex::Complex64;

fn exact(value: f64) -> Rat {
    for denominator in 1..=96i128 {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < 1e-9 {
            return Rat::new(scaled.round() as i128, denominator).unwrap();
        }
    }
    panic!("off grid: {value}");
}

fn class_of(
    operation: &MagneticOperation,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
) -> usize {
    let canonical = canonical_translation(operation.translation, lattice);
    classes
        .iter()
        .position(|class| {
            class.rotation == operation.rotation
                && class.time_reversal == operation.time_reversal
                && canonical_translation(class.translation, lattice) == canonical
        })
        .expect("the product of two classes must be a class")
}

/// Re-derive the law, the square, conjugation compatibility and the intertwiner
/// contract from the result and the model, using only public API.
fn verify_reconstruction(model: &ScalarGammaModel, report: &MagneticSubduction<'_>) {
    let classes = report.classes();
    let lattice = report.translation_lattice();
    let operators = report.operators();
    assert_eq!(operators.len(), classes.len(), "one operator per class");
    for (left, left_class) in classes.iter().enumerate() {
        for (right, right_class) in classes.iter().enumerate() {
            let product = compose(left_class, right_class).unwrap();
            let index = class_of(&product, classes, lattice);
            assert_eq!(
                operators[left].compose(&operators[right]),
                operators[index],
                "corepresentation law at classes ({left}, {right})"
            );
        }
    }
    if let Some(first_anti) = classes.iter().position(|class| class.time_reversal) {
        let square = class_of(
            &compose(&classes[first_anti], &classes[first_anti]).unwrap(),
            classes,
            lattice,
        );
        assert!(!classes[square].time_reversal, "a0^2 must be unitary");
        assert_eq!(
            operators[first_anti].compose(&operators[first_anti]),
            operators[square],
            "D(a0)^2 = D(a0^2)"
        );
        for (index, class) in classes.iter().enumerate() {
            if class.time_reversal {
                continue;
            }
            let left = class_of(&compose(class, &classes[first_anti]).unwrap(), classes, lattice);
            let right = class_of(&compose(&classes[first_anti], class).unwrap(), classes, lattice);
            assert_eq!(
                operators[index].compose(&operators[first_anti]),
                operators[left],
                "D(h) D(a0) = D(h a0) at class {index}"
            );
            assert_eq!(
                operators[first_anti].compose(&operators[index]),
                operators[right],
                "D(a0) conj(D(h)) = D(a0 h) at class {index}"
            );
        }
    }
    // The exhibited intertwiner satisfies its contract on every class.
    let intertwiner = report.intertwiner();
    assert!(intertwiner.norm() > 1e-12);
    let inverse = intertwiner.conj() / intertwiner.norm_sqr();
    for (index, action) in model.actions().iter().enumerate() {
        let source = Complex64::new(f64::from(action.value.value()), 0.0);
        let factor = if action.class.time_reversal {
            intertwiner.conj()
        } else {
            intertwiner
        };
        let transformed = inverse * source * factor;
        assert!(
            (transformed - operators[index].value()).norm() < 1e-12,
            "intertwiner contract at class {index}"
        );
        assert_eq!(
            operators[index].is_antiunitary(),
            action.class.time_reversal,
            "the operator's antiunitary flag is the class's"
        );
        assert!((operators[index].value().norm() - 1.0).abs() < 1e-12);
        assert!(operators[index].value().im.abs() < 1e-12);
    }
    // The same contract on actual complex vectors, not only on scalars:
    // A(g)(U z or conj(U z)) = U D(g) z for both branches.
    for z in [Complex64::new(1.0, 2.0), Complex64::new(-3.0, 0.5)] {
        for (index, action) in model.actions().iter().enumerate() {
            let source = Complex64::new(f64::from(action.value.value()), 0.0);
            let left = source
                * if action.class.time_reversal {
                    (intertwiner * z).conj()
                } else {
                    intertwiner * z
                };
            let right = intertwiner * operators[index].apply(z);
            assert!(
                (left - right).norm() < 1e-12,
                "intertwiner vector contract at class {index}"
            );
        }
    }
    assert_eq!(report.covered_dimension(), report.source_dimension());
}

/// Cross-check the reconstructed unitary columns against the summary's own
/// characters of the reported target label.
fn verify_against_summary(report: &MagneticSubduction<'_>) {
    let uni = report.model().uni();
    let partial = magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni).unwrap();
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.coords == (0, 0, 0, 1))
        .expect("a Gamma k-point");
    let corep = gamma
        .coreps
        .iter()
        .find(|corep| corep.label == report.target_label())
        .expect("the reported target label must be a summary label");
    assert_eq!(corep.dim, report.target_dimension() as usize);
    let classes = report.classes();
    let lattice = report.translation_lattice();
    for (column, operation) in gamma.operations.iter().enumerate() {
        if operation.time_reversal {
            continue;
        }
        let mapped = MagneticOperation {
            rotation: operation.rotation,
            translation: [
                exact(operation.translation[0]),
                exact(operation.translation[1]),
                exact(operation.translation[2]),
            ],
            time_reversal: operation.time_reversal,
        };
        let index = class_of(&mapped, classes, lattice);
        let character = corep.characters[column].expect("a unitary column has a character");
        assert!(
            (character - report.operators()[index].value()).norm() < 1e-12,
            "unitary column {column} disagrees with the summary"
        );
    }
}

fn gamma_model(
    parent_sg: u8,
    source: &str,
    condensing: &str,
    direction_index: usize,
    eta: Option<RealSign>,
    map: Mat3R,
) -> ScalarGammaModel {
    let request = ScalarGammaRequest::new(
        SourceIdentity {
            parent_sg,
            label: source,
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg,
            condensing_label: condensing,
            convention: LabelConvention::Cdml,
            direction_index,
        },
        eta,
    )
    .with_embedding(EmbeddingEvidence {
        map,
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::ContainmentOnly,
    });
    verify_scalar_gamma_input(&request).expect("the pinned R25 model must verify")
}

fn magnetic_row(parent_sg: u8, condensing: &str, uni: usize) -> usize {
    isotropy::magnetic_isotropy_subgroups(parent_sg, condensing, LabelConvention::Cdml)
        .unwrap()
        .iter()
        .position(|row| row.record.mag_sg == uni)
        .expect("the pinned magnetic isotropy row")
}

#[test]
fn type_one_ordinary_agrees_with_the_ordinary_subduction() {
    // SG 2 GM1+ with the identity embedding is UNI 4 (P-1, the parent itself).
    let index = magnetic_row(2, "GM1+", 4);
    let model = gamma_model(2, "GM1+", "GM1+", index, None, Mat3R::identity());
    let report = subduce_scalar_gamma(&model).unwrap();
    assert_eq!(report.target_unitary_sg(), 2);
    assert_eq!(report.target_dimension(), 1);
    assert_eq!(report.multiplicity(), 1);
    assert_eq!(report.source_dimension(), 1);
    assert_eq!(report.covered_dimension(), 1);
    assert!(!report.classes().iter().any(|class| class.time_reversal));
    assert_eq!(report.intertwiner(), Complex64::new(1.0, 0.0));
    verify_reconstruction(&model, &report);
    verify_against_summary(&report);

    // The ordinary subduction entry point on the same subgroup must name the same
    // target with the same dimension and multiplicity.
    let direction = isotropy::magnetic_isotropy_subgroups(2, "GM1+", LabelConvention::Cdml).unwrap()
        [index]
    .record
    .direction;
    let ordinary = subduce(
        &SubductionRequest::new(2, "GM1+", LabelConvention::Cdml)
            .direction(IsotropyDirection::Label(direction)),
    )
    .unwrap();
    let targets: Vec<_> = ordinary
        .blocks()
        .iter()
        .flat_map(|block| block.targets.iter())
        .collect();
    assert_eq!(ordinary.direction().subgroup_sg, 2);
    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].cdml, Some(report.target_label()));
    assert_eq!(u32::from(targets[0].dimension), report.target_dimension());
    assert_eq!(targets[0].multiplicity, report.multiplicity());
}

#[test]
fn type_three_signs_share_the_target_and_differ_by_the_intertwiner() {
    // SG 2 GM1- is the black-white group UNI 6 with H = P1: its 1D type-A corep
    // GM1 is induced from P1's scalar GM1, not from the source irrep.
    let index = magnetic_row(2, "GM1-", 6);
    let plus = gamma_model(2, "GM1-", "GM1-", index, Some(RealSign::Plus), Mat3R::identity());
    let minus = gamma_model(2, "GM1-", "GM1-", index, Some(RealSign::Minus), Mat3R::identity());
    let plus_report = subduce_scalar_gamma(&plus).unwrap();
    let minus_report = subduce_scalar_gamma(&minus).unwrap();

    assert_eq!(plus_report.target_label(), "GM1");
    assert_eq!(plus_report.target_source_label(), "GM1");
    assert_eq!(plus_report.target_unitary_sg(), 1);
    assert_eq!(plus_report.target_label(), minus_report.target_label());
    assert_eq!(
        plus_report.target_source_label(),
        minus_report.target_source_label()
    );
    assert_eq!(
        plus_report.target_dimension(),
        minus_report.target_dimension()
    );
    assert_eq!(plus_report.multiplicity(), minus_report.multiplicity());
    // Same target table, different necessary intertwiner: the eta sign is a gauge
    // of the antiunitary coset, not a different subduction.
    assert_eq!(plus_report.operators(), minus_report.operators());
    // The source character at the parent's reflected operation is -1, so for
    // eta = +1 the antiunitary class carries -1 and the intertwiner is i; for
    // eta = -1 it carries +1 and the intertwiner is 1.
    assert_eq!(plus_report.intertwiner(), Complex64::new(0.0, 1.0));
    assert_eq!(minus_report.intertwiner(), Complex64::new(1.0, 0.0));
    // The source actions do differ on the antiunitary class, so the two models are
    // genuinely different inputs to the same answer.
    assert_ne!(
        plus.actions()[1].value.value(),
        minus.actions()[1].value.value()
    );
    for (model, report) in [(&plus, &plus_report), (&minus, &minus_report)] {
        assert!(report
            .operators()
            .iter()
            .any(|operator| operator.is_antiunitary()));
        verify_reconstruction(model, report);
        verify_against_summary(report);
    }
}

#[test]
fn type_four_antiunitary_operator_is_reconstructed() {
    // SG 1 X1 condenses into the anti-translation group UNI 3, whose classes are
    // (E|0) and (E|(0,0,1/2), T); the map 2I is the pinned embedding there.
    let index = magnetic_row(1, "X1", 3);
    let two = Mat3R::from_ints([[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
    let plus = gamma_model(1, "GM1", "X1", index, Some(RealSign::Plus), two);
    let minus = gamma_model(1, "GM1", "X1", index, Some(RealSign::Minus), two);
    let plus_report = subduce_scalar_gamma(&plus).unwrap();
    let minus_report = subduce_scalar_gamma(&minus).unwrap();

    let antiunitary = minus.classes()[1];
    assert!(antiunitary.time_reversal);
    assert_ne!(antiunitary.translation, [Rat::ZERO; 3]);
    assert_eq!(minus_report.operators().len(), 2);
    assert!(!minus_report.operators()[0].is_antiunitary());
    assert!(minus_report.operators()[1].is_antiunitary());
    assert_eq!(
        minus_report.operators()[1].value(),
        Complex64::new(1.0, 0.0),
        "the canonical target makes the antiunitary representative act as 1.K"
    );
    assert_eq!(plus_report.intertwiner(), Complex64::new(1.0, 0.0));
    assert_eq!(minus_report.intertwiner(), Complex64::new(0.0, 1.0));
    for (model, report) in [(&plus, &plus_report), (&minus, &minus_report)] {
        verify_reconstruction(model, report);
        verify_against_summary(report);
    }
}

#[test]
fn the_target_identity_is_read_independently_of_the_source_label() {
    // SG 3 GM2 (a real +-1 source with chi(reflection) = -1) is UNI 10 with
    // H = P1: the target is labelled GM1, the scalar non-spinor corep of P1, so the
    // answer is not an echo of the source label and no spinor row is selected.
    let index = magnetic_row(3, "GM2", 10);
    let model = gamma_model(3, "GM2", "GM2", index, Some(RealSign::Minus), Mat3R::identity());
    let report = subduce_scalar_gamma(&model).unwrap();
    assert_eq!(model.source_label(), "GM2");
    assert_eq!(report.target_label(), "GM1");
    assert_eq!(report.target_source_label(), "GM1");
    assert_eq!(report.target_unitary_sg(), 1);
    assert_eq!(report.target_dimension(), 1);
    assert_eq!(report.multiplicity(), 1);
    verify_reconstruction(&model, &report);
    verify_against_summary(&report);

    // UNI 4's non-scalar candidates (GM2, GM3 of P-1) are spinor rows and must not
    // be selected even though they are one-dimensional.
    let type_one_index = magnetic_row(2, "GM1+", 4);
    let type_one = gamma_model(2, "GM1+", "GM1+", type_one_index, None, Mat3R::identity());
    let type_one_report = subduce_scalar_gamma(&type_one).unwrap();
    assert_eq!(type_one_report.target_label(), "GM1+");
    assert_eq!(type_one_report.target_source_label(), "GM1+");
}

#[test]
fn a_frame_whose_correspondence_is_not_established_is_refused_by_name() {
    // SG 5 GM1 with the GM2 record and the identity embedding builds the pinned
    // R25 model for UNI 22, but the summary's Γ operations for that UNI are given
    // in a setting whose rotations do not match the magnetic database's classes
    // (its antiunitary rotation is [[0,-1,0],[-1,0,0],[0,0,-1]] where the class has
    // [[-1,0,0],[0,1,0],[0,0,-1]]).  The correspondence is refused, not guessed, and
    // the refusal is not "the target does not exist".
    let index = magnetic_row(5, "GM2", 22);
    let model = gamma_model(
        5,
        "GM1",
        "GM2",
        index,
        Some(RealSign::Minus),
        Mat3R::identity(),
    );
    assert_eq!(model.uni(), 22);
    match subduce_scalar_gamma(&model) {
        Err(MagneticSubductionError::TargetOperationUnmatched {
            uni,
            time_reversal,
            ..
        }) => {
            assert_eq!(uni, 22);
            assert!(time_reversal, "the antiunitary column is the unmatched one");
        }
        other => panic!("expected a named frame refusal, got {other:?}"),
    }

    // SG 11 GM1+ with the GM2- record is UNI 59: there the rotations do match and
    // the translation-level mismatch is refused instead.
    let index = magnetic_row(11, "GM2-", 59);
    let model = gamma_model(
        11,
        "GM1+",
        "GM2-",
        index,
        Some(RealSign::Minus),
        Mat3R::identity(),
    );
    assert_eq!(model.uni(), 59);
    match subduce_scalar_gamma(&model) {
        Err(MagneticSubductionError::TargetOperationUnmatched { uni, .. }) => assert_eq!(uni, 59),
        other => panic!("expected a named frame refusal, got {other:?}"),
    }
}

#[test]
fn the_operators_are_the_actual_corepresentation_of_the_magnetic_classes() {
    // A direct operator-level statement on the antiunitary witness: composing the
    // reconstructed operators reproduces the group product on every ordered pair,
    // and the two unitarity branches really differ (the antiunitary operator
    // conjugates its input).
    let index = magnetic_row(1, "X1", 3);
    let two = Mat3R::from_ints([[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
    let model = gamma_model(1, "GM1", "X1", index, Some(RealSign::Minus), two);
    let report = subduce_scalar_gamma(&model).unwrap();
    let operators = report.operators();
    let z = Complex64::new(0.75, -0.25);
    let anti = operators[1];
    assert!((anti.apply(z) - anti.value() * z.conj()).norm() < 1e-12);
    let unitary = operators[0];
    assert!((unitary.apply(z) - unitary.value() * z).norm() < 1e-12);
    // (op1 . op1) lands on the identity class, exactly as the model's class
    // composition does.
    let classes = report.classes();
    let square = class_of(
        &compose(&classes[1], &classes[1]).unwrap(),
        classes,
        report.translation_lattice(),
    );
    assert_eq!(anti.compose(&anti), operators[square]);
    assert_eq!(square, 0);
    for every in operators {
        assert!(
            (every.apply(z).norm() - z.norm()).abs() < 1e-12,
            "an operator preserves the norm"
        );
    }
}

#[test]
fn repeated_antiunitary_constraints_keep_the_intertwiner() {
    // SG 10 GM1+ with the GM2+ record is UNI 53 (four classes, two antiunitary).
    // Both antiunitary classes repeat the same equation for U, so a solver that
    // mistakes a repeated consistent constraint for a contradiction loses the
    // equivalence; the pinned answers are U = 1 for eta = +1 and U = i for
    // eta = -1, and the vector-level contract is re-checked in
    // verify_reconstruction.
    for (eta, expected) in [
        (RealSign::Plus, Complex64::new(1.0, 0.0)),
        (RealSign::Minus, Complex64::new(0.0, 1.0)),
    ] {
        let model = gamma_model(10, "GM1+", "GM2+", 0, Some(eta), Mat3R::identity());
        assert_eq!(model.uni(), 53);
        assert_eq!(
            model.classes().iter().filter(|class| class.time_reversal).count(),
            2
        );
        let report = subduce_scalar_gamma(&model).unwrap();
        assert_eq!(report.target_label(), "GM1+");
        assert_eq!(report.target_dimension(), 1);
        assert_eq!(report.multiplicity(), 1);
        assert_eq!(report.intertwiner(), expected);
        verify_reconstruction(&model, &report);
        verify_against_summary(&report);
    }
}

#[test]
fn four_antiunitary_classes_with_a_nontrivial_h_character_are_reconstructed() {
    // SG 12 GM2+ with the Y1+ record is UNI 56: four antiunitary classes and a
    // unitary subgroup whose character takes both +1 and -1, so the coset is not a
    // single repeated equation.
    for (eta, expected) in [
        (RealSign::Minus, Complex64::new(1.0, 0.0)),
        (RealSign::Plus, Complex64::new(0.0, 1.0)),
    ] {
        let model = gamma_model(12, "GM2+", "Y1+", 0, Some(eta), Mat3R::identity());
        assert_eq!(model.uni(), 56);
        assert_eq!(
            model.classes().iter().filter(|class| class.time_reversal).count(),
            4
        );
        assert!(model
            .actions()
            .iter()
            .any(|action| !action.class.time_reversal && action.value.value() == -1));
        let report = subduce_scalar_gamma(&model).unwrap();
        assert_eq!(report.target_label(), "GM2+");
        assert_eq!(report.target_dimension(), 1);
        assert_eq!(report.multiplicity(), 1);
        assert_eq!(report.intertwiner(), expected);
        verify_reconstruction(&model, &report);
        verify_against_summary(&report);
    }
}
