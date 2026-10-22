//! R26/R27 regressions: the restricted scalar Γ magnetic subduction.
//!
//! Every test here re-derives its claim from the public API: the raw first-Hall
//! provenance of the summary's operation indices, the operator table against the
//! actual multiplication of the magnetic group, the intertwiner contract class by
//! class, and the target identity including its unitary characters from the magnetic
//! corep summary.  The module's own checks are not taken on faith.
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::MagneticOperation;
use cryspglib::irrep::magnetic_embedding::{canonical_translation, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaModel,
    ScalarGammaRequest, SourceIdentity, verify_scalar_gamma_input,
};
use cryspglib::irrep::magnetic_subduction::{MagneticSubduction, subduce_scalar_gamma};
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
            let left = class_of(
                &compose(class, &classes[first_anti]).unwrap(),
                classes,
                lattice,
            );
            let right = class_of(
                &compose(&classes[first_anti], class).unwrap(),
                classes,
                lattice,
            );
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

/// The Γ k-point summary of one UNI, cloned so a test can also tamper with it.
fn gamma_kpoint(uni: usize) -> magnetic_summary::MagneticKPointSummary {
    magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni)
        .unwrap()
        .summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.coords == (0, 0, 0, 1))
        .cloned()
        .expect("a Gamma k-point")
}

/// Re-derive the raw first-Hall provenance of a summary Γ table.
///
/// Checks that the two public raw surfaces are the same list element by element,
/// that the summary lists every raw operation exactly once (the Γ little group is the
/// whole group, so the index list is a permutation), that each column's displayed
/// time-reversal flag is the indexed raw operation's, and returns the class every
/// column maps to.
fn verify_raw_provenance(
    uni: usize,
    gamma: &magnetic_summary::MagneticKPointSummary,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
) -> Vec<usize> {
    use cryspglib::irrep::magnetic_embedding::magnetic_operations;
    let raw = magnetic_operations(uni).unwrap();
    let api = cryspglib::SymmetryOps::from_magnetic_database(uni).unwrap();
    assert_eq!(api.operations.len(), raw.operations.len(), "UNI {uni}");
    for (index, (api_operation, exact_operation)) in
        api.operations.iter().zip(&raw.operations).enumerate()
    {
        assert_eq!(
            api_operation.rotation, exact_operation.rotation,
            "UNI {uni} index {index}"
        );
        assert_eq!(
            api_operation.time_reversal, exact_operation.time_reversal,
            "UNI {uni} index {index}"
        );
        assert_eq!(
            [
                exact(api_operation.translation[0]),
                exact(api_operation.translation[1]),
                exact(api_operation.translation[2]),
            ],
            exact_operation.translation,
            "UNI {uni} index {index}"
        );
    }
    assert_eq!(
        gamma.operations.len(),
        raw.operations.len(),
        "UNI {uni}: the Γ little group is the whole group"
    );
    let mut seen = vec![false; raw.operations.len()];
    let mut mapped = Vec::with_capacity(gamma.operations.len());
    for (position, operation) in gamma.operations.iter().enumerate() {
        assert_eq!(
            operation.column, position,
            "UNI {uni}: the columns are the positions"
        );
        assert_eq!(
            operation.magnetic_operation_index, position,
            "UNI {uni}: the summary enumerates the raw list in order"
        );
        let raw_operation = raw.operations[operation.magnetic_operation_index];
        assert_eq!(
            raw_operation.time_reversal, operation.time_reversal,
            "UNI {uni} column {position}: the index no longer names this operation"
        );
        assert!(
            !seen[operation.magnetic_operation_index],
            "UNI {uni}: raw operation {} is indexed twice",
            operation.magnetic_operation_index
        );
        seen[operation.magnetic_operation_index] = true;
        mapped.push(class_of(&raw_operation, classes, lattice));
    }
    assert!(
        seen.iter().all(|seen| *seen),
        "UNI {uni}: every raw operation must be a column"
    );
    mapped
}

/// Cross-check the reconstructed operators against the summary's own characters.
///
/// The target identity is read from the summary, not echoed from the source: the
/// corep's own source row belongs to the summary's unitary subgroup `H`, is
/// one-dimensional and non-spinor, and the label is H's own Γ scalar irrep re-read
/// from the irrep database.  Every unitary column's character must equal the
/// reconstructed operator of the class its **operation index** names.
fn verify_against_summary(report: &MagneticSubduction<'_>) {
    let uni = report.model().uni();
    let partial = magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni).unwrap();
    let gamma = gamma_kpoint(uni);
    let corep = gamma
        .coreps
        .iter()
        .find(|corep| corep.label == report.target_label())
        .expect("the reported target label must be a summary label");
    assert_eq!(corep.dim, report.target_dimension() as usize);
    assert_eq!(corep.source_irreps.len(), 1);
    assert_eq!(corep.source_irreps[0].sg, report.target_unitary_sg());
    assert_eq!(corep.source_irreps[0].ml, report.target_source_label());
    assert_eq!(corep.source_irreps[0].dim, 1);
    assert!(!corep.source_irreps[0].spinor);
    assert_eq!(partial.summary.unitary_sg, report.target_unitary_sg());
    // Independent identity anchor: H itself must have exactly one scalar, non-spinor,
    // one-dimensional Γ irrep under the reported label.
    let candidates = cryspglib::irrep::query::find_irreps(
        report.target_unitary_sg(),
        report.target_source_label(),
        LabelConvention::Cdml,
    );
    let matches: Vec<_> = candidates
        .iter()
        .filter(|irrep| {
            irrep.dim == 1
                && !irrep.spinor
                && irrep.k_vector().numerators == [0, 0, 0]
                && irrep.k_vector().denominator == 1
        })
        .collect();
    assert_eq!(
        matches.len(),
        1,
        "the target's label must be H's own scalar Γ irrep"
    );

    let classes = report.classes();
    let lattice = report.translation_lattice();
    let raw = cryspglib::irrep::magnetic_embedding::magnetic_operations(uni).unwrap();
    for (column, operation) in gamma.operations.iter().enumerate() {
        let raw_operation = raw.operations[operation.magnetic_operation_index];
        assert_eq!(raw_operation.time_reversal, operation.time_reversal);
        if operation.time_reversal {
            continue;
        }
        let index = class_of(&raw_operation, classes, lattice);
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
    let direction = isotropy::magnetic_isotropy_subgroups(2, "GM1+", LabelConvention::Cdml)
        .unwrap()[index]
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
    let plus = gamma_model(
        2,
        "GM1-",
        "GM1-",
        index,
        Some(RealSign::Plus),
        Mat3R::identity(),
    );
    let minus = gamma_model(
        2,
        "GM1-",
        "GM1-",
        index,
        Some(RealSign::Minus),
        Mat3R::identity(),
    );
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
        assert!(
            report
                .operators()
                .iter()
                .any(|operator| operator.is_antiunitary())
        );
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
    let model = gamma_model(
        3,
        "GM2",
        "GM2",
        index,
        Some(RealSign::Minus),
        Mat3R::identity(),
    );
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

/// The Γ position association the module relies on, re-derived from the public API.
///
/// The summary builds its Γ columns by enumerating the little-group indices of the
/// operation list in increasing order, and the Γ little group is the whole group:
/// `filter_little_group_with_transform` returns `0..len` in order for every UNI (k = 0
/// is preserved by every operation, independently of the frame and of the pure
/// translations, since the reciprocal shift is zero).  So at position `i` the
/// character column must be `i` and the raw operation index must be `i` — on a real
/// summary too, which these pinned witnesses check.
#[test]
fn the_gamma_little_group_is_the_raw_list_in_order() {
    use cryspglib::irrep::wigner::filter_little_group_with_transform;
    let mut checked = 0usize;
    for uni in 1..=1651usize {
        let ops = cryspglib::SymmetryOps::from_magnetic_database(uni).unwrap();
        let little = filter_little_group_with_transform(0, 0, 0, 1, &ops, None, None);
        assert_eq!(
            little,
            (0..ops.operations.len()).collect::<Vec<_>>(),
            "UNI {uni}: the Γ little group is not the raw list in order"
        );
        checked += 1;
    }
    assert_eq!(checked, 1651, "every UNI was checked");
    for uni in [3usize, 4, 6, 10, 22, 53, 56, 59] {
        let gamma = gamma_kpoint(uni);
        assert!(!gamma.operations.is_empty());
        for (position, operation) in gamma.operations.iter().enumerate() {
            assert_eq!(operation.column, position, "UNI {uni} position {position}");
            assert_eq!(
                operation.magnetic_operation_index, position,
                "UNI {uni} position {position}"
            );
        }
    }
}

/// The raw first-Hall operation list is the same list on both public surfaces.
///
/// The summary's `magnetic_operation_index` is only a provenance link if the exact
/// list the classes are built from
/// ([`magnetic_embedding::magnetic_operations`]) is the very list the summary was
/// computed from ([`cryspglib::SymmetryOps::from_magnetic_database`]).  This re-checks
/// that claim for every UNI, element by element, instead of assuming it.
#[test]
fn the_two_raw_operation_surfaces_agree_for_every_uni() {
    use cryspglib::irrep::magnetic_embedding::magnetic_operations;
    let mut checked = 0usize;
    for uni in 1..=1651usize {
        let set = magnetic_operations(uni).unwrap_or_else(|error| panic!("UNI {uni}: {error}"));
        let api = cryspglib::SymmetryOps::from_magnetic_database(uni)
            .unwrap_or_else(|error| panic!("UNI {uni}: {error:?}"));
        assert_eq!(
            api.operations.len(),
            set.operations.len(),
            "UNI {uni}: the two raw lists differ in length"
        );
        for (index, (api_operation, exact_operation)) in
            api.operations.iter().zip(&set.operations).enumerate()
        {
            assert_eq!(
                api_operation.rotation, exact_operation.rotation,
                "UNI {uni} index {index}: rotation"
            );
            assert_eq!(
                api_operation.time_reversal, exact_operation.time_reversal,
                "UNI {uni} index {index}: time reversal"
            );
            assert_eq!(
                [
                    exact(api_operation.translation[0]),
                    exact(api_operation.translation[1]),
                    exact(api_operation.translation[2]),
                ],
                exact_operation.translation,
                "UNI {uni} index {index}: translation"
            );
        }
        checked += 1;
    }
    assert_eq!(checked, 1651, "every UNI was compared");
}

/// The index-resolved witnesses of the first restricted frame pair.
///
/// SG 5 `GM1` with the `GM2` record row 0 is UNI 22 and SG 11 `GM1+` with the `GM2-`
/// record row 0 is UNI 59.  Their summary Γ tables are displayed in a data-Hall frame
/// whose rotations/translations are **not** the magnetic classes' own frame, so the
/// correspondence is established through the operation index, and every claim is
/// re-derived here from the public API: the raw surface identity, the index
/// permutation, the unitary characters of the independently read target, the full
/// corepresentation law and the non-zero complex intertwiner.  Both `η` signs give the
/// same target table and differ only by the intertwiner.
#[test]
fn the_indexed_frame_provenance_yields_independent_character_targets() {
    for (sg, source, condensing, uni, unitary_sg) in [
        (5u8, "GM1", "GM2", 22usize, 1u8),
        (11, "GM1+", "GM2-", 59, 6),
    ] {
        let index = magnetic_row(sg, condensing, uni);
        let mut answers = Vec::new();
        let model = gamma_model(
            sg,
            source,
            condensing,
            index,
            Some(RealSign::Plus),
            Mat3R::identity(),
        );
        let gamma = gamma_kpoint(uni);
        let provenance =
            verify_raw_provenance(uni, &gamma, model.classes(), model.translation_lattice());
        assert!(
            provenance
                .iter()
                .all(|class| *class < model.classes().len()),
            "UNI {uni}: every column maps to a class"
        );
        assert!(
            (0..model.classes().len()).all(|class| provenance.contains(&class)),
            "UNI {uni}: every class is covered by a column"
        );
        let antiunitary_columns = gamma
            .operations
            .iter()
            .filter(|operation| operation.time_reversal)
            .count();
        assert!(antiunitary_columns > 0);
        for (eta, intertwiner) in [
            (RealSign::Plus, Complex64::new(1.0, 0.0)),
            (RealSign::Minus, Complex64::new(0.0, 1.0)),
        ] {
            let model = gamma_model(sg, source, condensing, index, Some(eta), Mat3R::identity());
            assert_eq!(model.uni(), uni);
            let report = subduce_scalar_gamma(&model)
                .unwrap_or_else(|error| panic!("UNI {uni} eta {eta:?}: {error}"));
            assert_eq!(report.target_label(), "GM1");
            assert_eq!(report.target_source_label(), "GM1");
            assert_eq!(
                report.target_unitary_sg(),
                unitary_sg,
                "UNI {uni}: the target's space group is read from H, not echoed from SG {sg}"
            );
            assert_ne!(report.target_unitary_sg(), sg);
            assert_eq!(report.target_dimension(), 1);
            assert_eq!(report.multiplicity(), 1);
            assert_eq!(report.source_dimension(), 1);
            assert_eq!(report.covered_dimension(), 1);
            assert_eq!(report.intertwiner(), intertwiner);
            assert!(
                report
                    .operators()
                    .iter()
                    .any(|operator| operator.is_antiunitary())
            );
            verify_reconstruction(&model, &report);
            verify_against_summary(&report);
            answers.push(report.operators().to_vec());
        }
        assert_eq!(
            answers[0], answers[1],
            "UNI {uni}: η is a gauge of the antiunitary coset, not a different target"
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
            model
                .classes()
                .iter()
                .filter(|class| class.time_reversal)
                .count(),
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
            model
                .classes()
                .iter()
                .filter(|class| class.time_reversal)
                .count(),
            4
        );
        assert!(
            model
                .actions()
                .iter()
                .any(|action| !action.class.time_reversal && action.value.value() == -1)
        );
        let report = subduce_scalar_gamma(&model).unwrap();
        assert_eq!(report.target_label(), "GM2+");
        assert_eq!(report.target_dimension(), 1);
        assert_eq!(report.multiplicity(), 1);
        assert_eq!(report.intertwiner(), expected);
        verify_reconstruction(&model, &report);
        verify_against_summary(&report);
    }
}
