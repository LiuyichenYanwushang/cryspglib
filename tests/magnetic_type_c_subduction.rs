//! R32 public regression: the restricted paired type-C Γ magnetic subduction.
//!
//! The fixed witnesses are SG 89 `GM5` condensed along `GM2` direction 0
//! (UNI 749, H = SG 16) and along `GM4` direction 0 (UNI 751, H = SG 21).  In
//! both records the two-dimensional real Γ source restricts to
//! `GM3 ⊕ GM4` of H (two inequivalent real `±1` one-dimensional irreps), and the
//! antiunitary coset exchanges the pair.  The tests exercise the low-level
//! paired engine and the high-level owned report, both `η` signs, the full
//! target law, the complex-vector reconstruction, the two independent CDML/BC
//! H identities, and the unchanged legacy dispatch.

use num_complex::Complex64;

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::{MagneticOperation, canonical_translation, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::{
    MatrixGammaModel, MatrixGammaRequest, MatrixSubductionPiece, MatrixTypeCPiece,
    subduce_matrix_gamma, subduce_matrix_gamma_type_c, verify_matrix_gamma_input,
};
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticRecordSelector, MagneticSubductionApiError, MagneticSubductionRequest,
    MagneticTypeCReport, subduce, subduce_type_c,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};

const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

fn matrix_request(
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    row: usize,
    eta: Option<RealSign>,
) -> MatrixGammaRequest<'static> {
    MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: sg,
            condensing_label: condensing,
            convention: LabelConvention::Cdml,
            direction_index: row,
        },
        eta,
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::from_ints(IDENTITY),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    })
}

fn paired_model(condensing: &'static str, eta: RealSign) -> MatrixGammaModel {
    verify_matrix_gamma_input(&matrix_request(89, "GM5", condensing, 0, Some(eta)))
        .expect("the frozen paired witness verifies")
}

fn high_level_request(
    condensing: &'static str,
    uni: usize,
    eta: RealSign,
) -> MagneticSubductionRequest<'static> {
    MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: "GM5",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(89, condensing, LabelConvention::Cdml)
            .direction(IsotropyDirection::Index(0))
            .uni(uni),
        Some(eta),
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::from_ints(IDENTITY),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    })
}

/// Apply one real target operator to a complex vector, with its unitarity flag.
fn target_apply(
    operator: &cryspglib::irrep::magnetic_matrix_subduction::TargetMatrixOperator,
    z: &[Complex64],
) -> Vec<Complex64> {
    operator.apply(z)
}

/// Apply one real source operator to a complex vector.
fn source_apply(
    operator: &cryspglib::irrep::magnetic_matrix_subduction::SourceMatrixOperator,
    z: [Complex64; 2],
) -> [Complex64; 2] {
    operator.apply(z)
}

/// Multiply one `2×2` copy by a target vector.
fn lift(copy: &[Complex64], z: &[Complex64]) -> Vec<Complex64> {
    (0..2)
        .map(|row| {
            (0..2)
                .map(|column| copy[row * 2 + column] * z[column])
                .sum::<Complex64>()
        })
        .collect()
}

/// One copy's Gram matrix `U†U`.
fn gram(copy: &[Complex64]) -> [Complex64; 4] {
    let mut out = [Complex64::new(0.0, 0.0); 4];
    for row in 0..2 {
        for column in 0..2 {
            out[row * 2 + column] = (0..2)
                .map(|k| copy[k * 2 + row].conj() * copy[k * 2 + column])
                .sum();
        }
    }
    out
}

fn class_index(model: &MatrixGammaModel, operation: &MagneticOperation) -> Option<usize> {
    let lattice = model.translation_lattice();
    model.classes().iter().position(|class| {
        class.rotation == operation.rotation
            && class.time_reversal == operation.time_reversal
            && canonical_translation(
                [
                    class.translation[0]
                        .checked_sub(operation.translation[0])
                        .expect("exact translation difference"),
                    class.translation[1]
                        .checked_sub(operation.translation[1])
                        .expect("exact translation difference"),
                    class.translation[2]
                        .checked_sub(operation.translation[2])
                        .expect("exact translation difference"),
                ],
                lattice,
            )
            .iter()
            .all(|value| value.is_zero())
    })
}

/// The complete ordered-pair target corepresentation law, recomputed from the
/// model's own class composition.
fn assert_target_law(model: &MatrixGammaModel, piece: &MatrixTypeCPiece) {
    let classes = model.classes();
    for left in 0..classes.len() {
        for right in 0..classes.len() {
            let product = compose(&classes[left], &classes[right]).expect("class product");
            let expected = class_index(model, &product).expect("product is a class");
            let left_operator = &piece.operators()[left];
            let right_operator = &piece.operators()[right];
            // conj(A(h)) = A(h) on the real certified entries.
            let factor = right_operator.entries().to_vec();
            let mut found = vec![0.0; 4];
            for row in 0..2 {
                for column in 0..2 {
                    found[row * 2 + column] = (0..2)
                        .map(|k| left_operator.entries()[row * 2 + k] * factor[k * 2 + column])
                        .sum();
                }
            }
            assert_eq!(
                left_operator.is_antiunitary() ^ right_operator.is_antiunitary(),
                piece.operators()[expected].is_antiunitary(),
                "the antiunitary flag of the product class"
            );
            for entry in 0..4 {
                assert!(
                    (found[entry] - piece.operators()[expected].entries()[entry]).abs() < 1e-9,
                    "target law fails at classes ({left}, {right}) -> {expected}: {found:?} vs {:?}",
                    piece.operators()[expected].entries()
                );
            }
        }
    }
}

/// `S(g)∘φ = φ∘T(g)` on complex vectors, for every class and every copy.
fn assert_complex_reconstruction(model: &MatrixGammaModel, piece: &MatrixTypeCPiece) {
    let probes = [
        [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
        [Complex64::new(0.7, 1.3), Complex64::new(-1.1, 0.4)],
    ];
    for copy in piece.copies() {
        assert_eq!(copy.len(), 4);
        for probe in &probes {
            let image = lift(copy, probe);
            for class in 0..model.classes().len() {
                let source = source_apply(&model.source_operators()[class], [image[0], image[1]]);
                let target = target_apply(&piece.operators()[class], probe);
                let expected = lift(copy, &target);
                for row in 0..2 {
                    assert!(
                        (source[row] - expected[row]).norm() < 1e-8,
                        "reconstruction fails at class {class}"
                    );
                }
            }
        }
    }
}

fn assert_piece_metric(piece: &MatrixTypeCPiece) {
    assert_eq!(piece.dimension(), 2);
    assert_eq!(piece.multiplicity(), 1);
    assert_eq!(piece.covered_dimension(), 2);
    assert_eq!(piece.copies().len(), 1);
    let copy = &piece.copies()[0];
    let value = gram(copy);
    for row in 0..2 {
        for column in 0..2 {
            let expected = if row == column {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            assert!((value[row * 2 + column] - expected).norm() < 1e-8);
        }
    }
    let determinant = copy[0] * copy[3] - copy[1] * copy[2];
    assert!(
        (determinant.norm() - 1.0).abs() < 1e-8,
        "the single copy is not full rank: |det| = {}",
        determinant.norm()
    );
}

fn assert_high_level_payload(report: &MagneticTypeCReport) {
    let outcome = report.outcome();
    assert_eq!(outcome.target.dimension, 2);
    assert_eq!(outcome.target.multiplicity, 1);
    assert_eq!(outcome.target.covered_dimension(), 2);
    assert_eq!(outcome.covered_dimension, 2);
    assert_eq!(outcome.copies.len(), 1);
    let copy = &outcome.copies[0];
    let value = gram(copy);
    for row in 0..2 {
        for column in 0..2 {
            let expected = if row == column {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            assert!((value[row * 2 + column] - expected).norm() < 1e-8);
        }
    }
    let probes = [
        [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        [Complex64::new(0.3, -1.7), Complex64::new(1.2, 0.9)],
    ];
    let model = report.model();
    for probe in &probes {
        let image = lift(copy, probe);
        for class in 0..model.classes().len() {
            let source = source_apply(&model.source_operators()[class], [image[0], image[1]]);
            let target = outcome.operators[class].apply(probe);
            let expected = lift(copy, &target);
            for row in 0..2 {
                assert!((source[row] - expected[row]).norm() < 1e-8);
            }
        }
    }
}

#[test]
fn both_fixed_witnesses_are_paired_type_c_targets_for_both_eta() {
    for (condensing, expected_uni, expected_h_sg) in
        [("GM2", 749usize, 16u8), ("GM4", 751usize, 21u8)]
    {
        let mut pieces = Vec::new();
        let mut copies: Vec<Vec<Complex64>> = Vec::new();
        for eta in [RealSign::Plus, RealSign::Minus] {
            let model = paired_model(condensing, eta);
            assert_eq!(model.uni(), expected_uni);
            assert_eq!(model.source_label(), "GM5");
            assert_eq!(model.record_label(), condensing);
            assert_ne!(model.source_label(), model.record_label());
            assert_eq!(model.classes().len(), 8);
            assert_eq!(
                model
                    .classes()
                    .iter()
                    .filter(|class| class.time_reversal)
                    .count(),
                4
            );
            let result = subduce_matrix_gamma_type_c(&model)
                .expect("the paired type-C witness reconstructs");
            let piece = result.piece();
            assert_eq!(result.source_dimension(), 2);
            assert_eq!(result.covered_dimension(), 2);
            assert_eq!(piece.unitary_sg(), expected_h_sg);
            assert_eq!(piece.h_sources(), ["GM3", "GM4"]);
            assert_eq!(piece.label(), "GM3 + GM4");
            assert_piece_metric(piece);
            assert_target_law(&model, piece);
            assert_complex_reconstruction(&model, piece);
            // The canonical type-C table has a diagonal unitary part and an
            // off-diagonal antiunitary exchange on every class.
            for (class, operator) in piece.operators().iter().enumerate() {
                if model.classes()[class].time_reversal {
                    assert_eq!(operator.entries()[0], 0.0);
                    assert_eq!(operator.entries()[3], 0.0);
                    assert_ne!(operator.entries()[1], 0.0);
                    assert_ne!(operator.entries()[2], 0.0);
                } else {
                    assert_eq!(operator.entries()[1], 0.0);
                    assert_eq!(operator.entries()[2], 0.0);
                }
            }
            pieces.push(piece.operators().to_vec());
            copies.push(piece.copies()[0].clone());
        }
        // The canonical target table is η-independent; only the intertwiner
        // absorbs the parent gauge, so the two η copies must differ.
        assert_eq!(pieces[0], pieces[1]);
        let difference: f64 = copies[0]
            .iter()
            .zip(copies[1].iter())
            .map(|(left, right)| (*left - *right).norm())
            .fold(0.0f64, f64::max);
        assert!(
            difference > 0.5,
            "{condensing}: η must change the copy, difference {difference}"
        );
    }
}

#[test]
fn the_high_level_report_carries_the_two_real_h_identities() {
    for (condensing, expected_uni, expected_h_sg) in
        [("GM2", 749usize, 16u8), ("GM4", 751usize, 21u8)]
    {
        for eta in [RealSign::Plus, RealSign::Minus] {
            let report = subduce_type_c(&high_level_request(condensing, expected_uni, eta))
                .expect("the high-level paired request is answered");
            assert_eq!(report.record().uni, expected_uni);
            assert_eq!(report.source().labels.cdml, "GM5");
            assert_eq!(report.condensing().labels.cdml, condensing);
            assert_eq!(report.source_dimension(), 2);
            assert_eq!(report.covered_dimension(), 2);
            let target = report.target();
            assert_eq!(target.magnetic_label, "GM3 + GM4");
            assert_eq!(target.magnetic_bc, None);
            assert_eq!(target.unitary_sg, expected_h_sg);
            assert_eq!(target.h_sources[0].cdml, "GM3");
            assert_eq!(target.h_sources[1].cdml, "GM4");
            assert_eq!(target.h_sources[0].bc.as_deref(), Some("Γ4"));
            assert_eq!(target.h_sources[1].bc.as_deref(), Some("Γ2"));
            assert_eq!(target.dimension, 2);
            assert_eq!(target.multiplicity, 1);
            assert_high_level_payload(&report);
            // The record's classes and lattice are exposed by the report.
            assert_eq!(report.classes().len(), report.model().classes().len());
            assert_eq!(
                report.translation_lattice(),
                report.model().translation_lattice()
            );
        }
    }
}

#[test]
fn the_legacy_type_a_dispatch_still_refuses_the_type_c_record_by_name() {
    // The old high-level dispatch keeps its type-A-only matrix domain: the
    // same request that the paired entry answers is refused explicitly, not
    // silently truncated to one H identity.
    let error = subduce(&high_level_request("GM2", 749, RealSign::Plus))
        .expect_err("the legacy matrix domain has no type-A target");
    match error {
        MagneticSubductionApiError::Matrix(inner) => match *inner {
            cryspglib::irrep::magnetic_matrix_subduction::MatrixSubductionError::NoMatchingTarget {
                uni,
                examined,
                ..
            } => {
                assert_eq!(uni, 749);
                assert_eq!(examined, vec!["GM1".to_string(), "GM2".to_string()]);
            }
            other => panic!("expected NoMatchingTarget, got {other}"),
        },
        other => panic!("expected the matrix error, got {other}"),
    }
}

#[test]
fn the_legacy_type_a_matrix_witnesses_are_unchanged() {
    // SG 89 GM5 / GM1 is the parent group itself (UNI 747): the old engine
    // still returns its one type-A two-dimensional piece.
    let model =
        verify_matrix_gamma_input(&matrix_request(89, "GM5", "GM1", 0, Some(RealSign::Plus)))
            .expect("the type-A witness verifies");
    assert_eq!(model.uni(), 747);
    let report = subduce_matrix_gamma(&model).expect("the type-A target still reconstructs");
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(report.pieces().len(), 1);
    let piece: &MatrixSubductionPiece = report.piece(0).expect("one piece");
    assert_eq!(piece.dimension(), 2);
    assert_eq!(piece.multiplicity(), 1);
    // The paired entry is not silently reusing the type-A result.
    assert!(subduce_matrix_gamma_type_c(&model).is_err());
}

#[test]
fn scalar_sources_are_refused_by_the_paired_entry_point() {
    let request = MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(89, "GM2", LabelConvention::Cdml)
            .direction(IsotropyDirection::Index(0))
            .uni(749),
        Some(RealSign::Plus),
    );
    match subduce_type_c(&request) {
        Err(
            cryspglib::irrep::magnetic_subduction_api::MagneticTypeCError::SourceDimensionNotTwo {
                sg,
                dimension,
                ..
            },
        ) => {
            assert_eq!(sg, 89);
            assert_eq!(dimension, 1);
        }
        other => panic!("expected SourceDimensionNotTwo, got {other:?}"),
    }
}

#[test]
fn a_type_a_record_has_no_paired_candidate_and_is_refused() {
    let report = subduce_type_c(&high_level_request("GM1", 747, RealSign::Plus))
        .expect_err("UNI 747 has no type-C target");
    match report {
        cryspglib::irrep::magnetic_subduction_api::MagneticTypeCError::Engine(inner) => match *inner
        {
            cryspglib::irrep::magnetic_matrix_subduction::MatrixTypeCError::NoMatchingTarget {
                uni,
                examined,
                ..
            } => {
                assert_eq!(uni, 747);
                assert!(examined.is_empty());
            }
            other => panic!("expected NoMatchingTarget, got {other}"),
        },
        other => panic!("expected the engine refusal, got {other}"),
    }
}
