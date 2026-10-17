//! R29 public regression: the two-dimensional real Γ magnetic **matrix**
//! subduction.
//!
//! The witnesses are frozen real records: `source` and `condensing` are
//! *independent* identities (as in R25), and each one is reached with an
//! explicit, verified embedding.  The tests exercise
//!
//! * a genuine two-dimensional target with a non-trivial antiunitary extension
//!   (SG 89 `GM5` condensed along `M1`, UNI 753, an anti-translation group with
//!   eight antiunitary classes) and both `η` signs;
//! * a Type I ordinary two-dimensional target (SG 89 `GM5` condensed along
//!   `GM1`, UNI 747 — the record where the magnetic group is the parent itself);
//! * the `2 → 1 + 1` split into two inequivalent one-dimensional targets
//!   (UNI 101) with both `η` signs;
//! * the `2 × 1` multiplicity-two case (UNI 10) and the SG 156 `GM3` cross
//!   success (UNI 1, `H = P1`);
//! * the named refusals of the out-of-family inputs.

use num_complex::Complex64;

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::{
    MatrixGammaModel, MatrixGammaRequest, MatrixSubductionError, MatrixSubductionPiece,
    subduce_matrix_gamma, verify_matrix_gamma_input,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};

const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

fn request(
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    row: usize,
    map: [[i32; 3]; 3],
    eta: Option<RealSign>,
    policy: LatticePolicy,
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
        map: Mat3R::from_ints(map),
        shift: [Rat::ZERO; 3],
        lattice_policy: policy,
    })
}

fn witness_model(
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    row: usize,
    map: [[i32; 3]; 3],
    eta: Option<RealSign>,
) -> MatrixGammaModel {
    verify_matrix_gamma_input(&request(
        sg,
        source,
        condensing,
        row,
        map,
        eta,
        LatticePolicy::RequireExact,
    ))
    .expect("the frozen witness request verifies")
}

/// `U†U` of one `2×d` copy.
fn gram(copy: &[Complex64], dimension: usize) -> Vec<Complex64> {
    let mut out = vec![Complex64::new(0.0, 0.0); dimension * dimension];
    for row in 0..dimension {
        for column in 0..dimension {
            out[row * dimension + column] = (0..2)
                .map(|k| copy[k * dimension + row].conj() * copy[k * dimension + column])
                .sum();
        }
    }
    out
}

fn assert_isometric_copy(copy: &[Complex64], dimension: usize) {
    let value = gram(copy, dimension);
    for row in 0..dimension {
        for column in 0..dimension {
            let expected = if row == column {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            assert!(
                (value[row * dimension + column] - expected).norm() < 1e-8,
                "the copy is not isometric: {value:?}"
            );
        }
    }
}

/// The public reconstruction check: `S(g)∘φ = φ∘T(g)` for every class and every
/// copy, on complex probe vectors.
fn assert_reconstruction(model: &MatrixGammaModel, piece: &MatrixSubductionPiece) {
    let dimension = piece.dimension() as usize;
    let probes: Vec<Vec<Complex64>> = (0..dimension)
        .map(|index| {
            (0..dimension)
                .map(|column| {
                    if column == index {
                        Complex64::new(1.0, 0.0)
                    } else {
                        Complex64::new(0.0, 0.0)
                    }
                })
                .collect()
        })
        .chain(std::iter::once(
            (0..dimension)
                .map(|column| Complex64::new(1.0 + column as f64, 1.0 - 2.0 * column as f64))
                .collect(),
        ))
        .collect();
    for copy in piece.copies() {
        for probe in &probes {
            let image: Vec<Complex64> = (0..2)
                .map(|row| {
                    (0..dimension)
                        .map(|column| copy[row * dimension + column] * probe[column])
                        .sum::<Complex64>()
                })
                .collect();
            for class in 0..model.classes().len() {
                let source = model
                    .source_operator_at(class)
                    .expect("class")
                    .apply([image[0], image[1]]);
                let target = piece.operator_at(class).expect("class").apply(probe);
                let expected: Vec<Complex64> = (0..2)
                    .map(|row| {
                        (0..dimension)
                            .map(|column| copy[row * dimension + column] * target[column])
                            .sum::<Complex64>()
                    })
                    .collect();
                for row in 0..2 {
                    assert!(
                        (source[row] - expected[row]).norm() < 1e-8,
                        "the reconstruction fails at class {class}"
                    );
                }
            }
        }
    }
}

/// The `2 × 1` copies of one one-dimensional piece span the source space.
fn assert_copies_span(copies: &[Vec<Complex64>], dimension: usize) {
    let columns: Vec<Vec<Complex64>> = copies
        .iter()
        .flat_map(|copy| {
            (0..dimension).map(move |column| vec![copy[column], copy[dimension + column]])
        })
        .collect();
    assert_eq!(columns.len(), 2);
    let determinant = columns[0][0] * columns[1][1] - columns[0][1] * columns[1][0];
    assert!(
        determinant.norm() > 0.5,
        "the copies do not span the source space"
    );
}

#[test]
fn the_type_i_dim2_target_is_rebuilt_from_the_parent_irrep() {
    // SG 89 GM5 condensed along GM1 gives M = P422 itself (UNI 747): a Type I
    // ordinary two-dimensional Γ target.  The unitary part equals the source in
    // this record by construction, which is exactly why this test pins it as the
    // ordinary end of the family and not as the antiunitary witness.
    let model = witness_model(89, "GM5", "GM1", 0, IDENTITY, Some(RealSign::Plus));
    assert_eq!(model.uni(), 747);
    assert_eq!(model.source_label(), "GM5");
    assert_eq!(model.source_dimension(), 2);
    assert_eq!(model.classes().len(), 8);
    assert!(model.realises_record_lattice());
    assert!(
        model
            .source_operators()
            .iter()
            .all(|op| !op.is_antiunitary())
    );

    let report = subduce_matrix_gamma(&model).expect("the ordinary record reduces");
    assert_eq!(report.source_dimension(), 2);
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(report.pieces().len(), 1);
    let piece = report.piece(0).expect("one piece");
    assert_eq!(piece.label(), "GM5");
    assert_eq!(piece.unitary_sg(), 89);
    assert_eq!(piece.source_label(), "GM5");
    assert_eq!(piece.dimension(), 2);
    assert_eq!(piece.multiplicity(), 1);
    assert_eq!(piece.copies().len(), 1);
    assert_isometric_copy(&piece.copies()[0], 2);
    assert_reconstruction(&model, piece);
    // M = the parent, so every target operator must equal the certified source
    // operator, and the intertwiner is the identity up to its canonical phase.
    for class in 0..model.classes().len() {
        let target = piece.operator_at(class).expect("class");
        let source = model.source_operator_at(class).expect("class");
        assert_eq!(target.is_antiunitary(), source.is_antiunitary());
        for entry in 0..4 {
            assert!((target.entries()[entry] - source.entries()[entry]).abs() < 1e-9);
        }
    }
    let identity = piece.copies()[0].clone();
    assert!((identity[0] - Complex64::new(1.0, 0.0)).norm() < 1e-9);
    assert!((identity[3] - Complex64::new(1.0, 0.0)).norm() < 1e-9);

    // A purely unitary group does not use eta, so both stated signs agree.
    let negative = witness_model(89, "GM5", "GM1", 0, IDENTITY, Some(RealSign::Minus));
    let report_negative = subduce_matrix_gamma(&negative).expect("both eta signs are legal");
    assert_eq!(report_negative.pieces().len(), report.pieces().len());
    for (left, right) in report.pieces().iter().zip(report_negative.pieces()) {
        assert_eq!(left.label(), right.label());
        assert_eq!(left.dimension(), right.dimension());
        assert_eq!(left.multiplicity(), right.multiplicity());
        assert_eq!(left.operators(), right.operators());
        for (a, b) in left
            .copies()
            .iter()
            .flatten()
            .zip(right.copies().iter().flatten())
        {
            assert!((a - b).norm() < 1e-9);
        }
    }
}

#[test]
fn the_anti_translation_dim2_target_has_a_nontrivial_antiunitary_extension() {
    // SG 89 GM5 condensed along M1 gives UNI 753 (P_C422, anti-translation),
    // whose unitary subgroup is the parent but whose antiunitary coset is new.
    // The embedding is the record's own tabulated frame; it is not the identity
    // map, so the target is not read off the source.
    const MAP: [[i32; 3]; 3] = [[-1, -1, 0], [-1, 1, 0], [0, 0, -1]];
    let plus = witness_model(89, "GM5", "M1", 0, MAP, Some(RealSign::Plus));
    assert_eq!(plus.uni(), 753);
    assert_eq!(plus.classes().len(), 16);
    assert!(plus.realises_record_lattice());

    let report_plus = subduce_matrix_gamma(&plus).expect("the anti-translation record reduces");
    assert_eq!(report_plus.covered_dimension(), 2);
    assert_eq!(report_plus.pieces().len(), 1);
    let piece = report_plus.piece(0).expect("one piece");
    assert_eq!(piece.label(), "GM5");
    assert_eq!(piece.unitary_sg(), 89);
    assert_eq!(piece.dimension(), 2);
    assert_eq!(piece.multiplicity(), 1);
    assert_isometric_copy(&piece.copies()[0], 2);
    assert_reconstruction(&plus, piece);
    let anti = piece
        .operators()
        .iter()
        .filter(|operator| operator.is_antiunitary())
        .count();
    assert_eq!(anti, 8, "the target carries the anti-translation coset");
    // Every target operator stays in the certified orthogonal metric.
    for operator in piece.operators() {
        let dimension = operator.dimension();
        for row in 0..dimension {
            for column in 0..dimension {
                let value: f64 = (0..dimension)
                    .map(|k| {
                        operator.entry(k, row).expect("entry")
                            * operator.entry(k, column).expect("entry")
                    })
                    .sum();
                let expected = if row == column { 1.0 } else { 0.0 };
                assert!((value - expected).abs() < 1e-9);
            }
        }
    }

    // The two eta signs are a gauge of the antiunitary coset: both are legal,
    // and they give different intertwiners.
    let minus = witness_model(89, "GM5", "M1", 0, MAP, Some(RealSign::Minus));
    let report_minus = subduce_matrix_gamma(&minus).expect("both eta signs are legal");
    let piece_minus = report_minus.piece(0).expect("one piece");
    assert_eq!(piece_minus.label(), piece.label());
    assert_eq!(piece_minus.operators(), piece.operators());
    let difference: f64 = piece
        .copies()
        .iter()
        .flatten()
        .zip(piece_minus.copies().iter().flatten())
        .map(|(left, right)| (left - right).norm())
        .fold(0.0f64, f64::max);
    assert!(
        difference > 0.5,
        "eta must change the intertwiner, not the target"
    );
}

#[test]
fn the_two_dimensional_restriction_splits_into_two_one_dimensional_pieces() {
    // SG 89 GM5 condensed along GM5 gives UNI 101 (P2'2'2, H = SG 3 = P2): the
    // restriction is a genuine 2 -> 1 + 1 split.
    let model = witness_model(89, "GM5", "GM5", 0, IDENTITY, Some(RealSign::Plus));
    assert_eq!(model.uni(), 101);
    assert_eq!(model.classes().len(), 4);
    let report = subduce_matrix_gamma(&model).expect("the split reduces");
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(report.pieces().len(), 2);
    let mut labels: Vec<&str> = report.pieces().iter().map(|piece| piece.label()).collect();
    labels.sort_unstable();
    assert_eq!(labels, vec!["GM1", "GM2"]);
    for piece in report.pieces() {
        assert_eq!(piece.unitary_sg(), 3);
        assert_eq!(piece.dimension(), 1);
        assert_eq!(piece.multiplicity(), 1);
        assert_eq!(piece.copies().len(), 1);
        assert_isometric_copy(&piece.copies()[0], 1);
        assert_reconstruction(&model, piece);
    }
    // The two pieces are inequivalent one-dimensional coreps, and their copies
    // together are one orthonormal basis of the two-dimensional source space.
    let copies: Vec<Vec<Complex64>> = report
        .pieces()
        .iter()
        .flat_map(|piece| piece.copies().iter().cloned())
        .collect();
    assert_copies_span(&copies, 1);

    // Both eta signs work and give different intertwiners.
    let negative = witness_model(89, "GM5", "GM5", 0, IDENTITY, Some(RealSign::Minus));
    let report_negative = subduce_matrix_gamma(&negative).expect("both eta signs are legal");
    assert_eq!(report_negative.pieces().len(), 2);
    let mut different = 0usize;
    for (plus_piece, minus_piece) in report.pieces().iter().zip(report_negative.pieces()) {
        assert_eq!(plus_piece.label(), minus_piece.label());
        let size = plus_piece.copies()[0]
            .iter()
            .zip(minus_piece.copies()[0].iter())
            .filter(|(left, right)| (**left - **right).norm() > 0.5)
            .count();
        if size > 0 {
            different += 1;
        }
    }
    assert_eq!(different, 2, "eta changes both intertwiners");
}

#[test]
fn the_p2_prime_record_gives_multiplicity_two_with_orthonormal_copies() {
    // SG 89 GM5 condensed along GM5 direction C1 gives UNI 10 (P2'): the
    // one-dimensional target appears twice.
    let model = witness_model(89, "GM5", "GM5", 2, IDENTITY, Some(RealSign::Plus));
    assert_eq!(model.uni(), 10);
    assert_eq!(model.classes().len(), 2);
    let report = subduce_matrix_gamma(&model).expect("the multiplicity-two record reduces");
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(report.pieces().len(), 1);
    let piece = report.piece(0).expect("one piece");
    assert_eq!(piece.label(), "GM1");
    assert_eq!(piece.unitary_sg(), 1);
    assert_eq!(piece.dimension(), 1);
    assert_eq!(piece.multiplicity(), 2);
    assert_eq!(piece.copies().len(), 2);
    for copy in piece.copies() {
        assert_isometric_copy(copy, 1);
    }
    assert_copies_span(piece.copies(), 1);
    assert_reconstruction(&model, piece);

    let negative = witness_model(89, "GM5", "GM5", 2, IDENTITY, Some(RealSign::Minus));
    let report_negative = subduce_matrix_gamma(&negative).expect("both eta signs are legal");
    let piece_negative = report_negative.piece(0).expect("one piece");
    assert_eq!(piece_negative.multiplicity(), 2);
    assert_reconstruction(&negative, piece_negative);
    let difference: f64 = piece
        .copies()
        .iter()
        .flatten()
        .zip(piece_negative.copies().iter().flatten())
        .map(|(left, right)| (left - right).norm())
        .fold(0.0f64, f64::max);
    assert!(difference > 0.5, "eta must change the intertwiner");
}

#[test]
fn the_sg156_gm3_source_crosses_to_p1_with_multiplicity_two() {
    let model = witness_model(156, "GM3", "GM3", 2, IDENTITY, Some(RealSign::Plus));
    assert_eq!(model.uni(), 1);
    assert_eq!(model.classes().len(), 1);
    let report = subduce_matrix_gamma(&model).expect("the SG 156 cross record reduces");
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(report.pieces().len(), 1);
    let piece = report.piece(0).expect("one piece");
    assert_eq!(piece.label(), "GM1");
    assert_eq!(piece.unitary_sg(), 1);
    assert_eq!(piece.dimension(), 1);
    assert_eq!(piece.multiplicity(), 2);
    for copy in piece.copies() {
        assert_isometric_copy(copy, 1);
    }
    assert_copies_span(piece.copies(), 1);
    assert_reconstruction(&model, piece);
}

#[test]
fn out_of_family_sources_are_refused_by_name() {
    // A spinor source.
    let error = verify_matrix_gamma_input(&request(
        89,
        "GM6",
        "GM1",
        0,
        IDENTITY,
        Some(RealSign::Plus),
        LatticePolicy::RequireExact,
    ))
    .expect_err("spinor sources are out of family");
    assert!(matches!(
        error,
        MatrixSubductionError::SourceIsSpinor { sg: 89, .. }
    ));

    // A one-dimensional source.
    let error = verify_matrix_gamma_input(&request(
        89,
        "GM1",
        "GM1",
        0,
        IDENTITY,
        None,
        LatticePolicy::RequireExact,
    ))
    .expect_err("one-dimensional sources are out of family");
    assert!(matches!(
        error,
        MatrixSubductionError::SourceDimensionNotTwo { dimension: 1, .. }
    ));

    // A non-Γ source.
    let error = verify_matrix_gamma_input(&request(
        89,
        "X1",
        "GM1",
        0,
        IDENTITY,
        Some(RealSign::Plus),
        LatticePolicy::RequireExact,
    ))
    .expect_err("non-Gamma sources are out of family");
    assert!(matches!(
        error,
        MatrixSubductionError::SourceNotGamma { .. }
    ));

    // A compound source (the realified pair row of SG 75).
    let error = verify_matrix_gamma_input(&request(
        75,
        "GM3GM4",
        "GM1",
        0,
        IDENTITY,
        Some(RealSign::Plus),
        LatticePolicy::RequireExact,
    ))
    .expect_err("compound sources are out of family");
    assert!(matches!(
        error,
        MatrixSubductionError::SourceIsCompound { sg: 75, .. }
    ));

    // Source and record naming different parents.
    let cross_parent = MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: "GM5",
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: 156,
            condensing_label: "GM3",
            convention: LabelConvention::Cdml,
            direction_index: 2,
        },
        Some(RealSign::Plus),
    );
    let error = verify_matrix_gamma_input(&cross_parent)
        .expect_err("different parents cannot be requested together");
    assert!(matches!(
        error,
        MatrixSubductionError::ParentMismatch {
            source_sg: 89,
            record_sg: 156
        }
    ));
}

#[test]
fn record_and_embedding_refusals_are_named() {
    // Direction index outside the table.
    let error = verify_matrix_gamma_input(&request(
        89,
        "GM5",
        "GM1",
        99,
        IDENTITY,
        Some(RealSign::Plus),
        LatticePolicy::RequireExact,
    ))
    .expect_err("the direction index is out of range");
    assert!(matches!(
        error,
        MatrixSubductionError::DirectionOutOfRange { index: 99, .. }
    ));

    // A missing eta on a record with an antiunitary coset.
    let error = verify_matrix_gamma_input(&request(
        89,
        "GM5",
        "GM5",
        0,
        IDENTITY,
        None,
        LatticePolicy::RequireExact,
    ))
    .expect_err("an antiunitary coset needs eta");
    assert!(matches!(
        error,
        MatrixSubductionError::MissingTimeReversal { uni: 101 }
    ));

    // A map that does not embed the group.
    let error = verify_matrix_gamma_input(&request(
        89,
        "GM5",
        "GM1",
        0,
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        Some(RealSign::Plus),
        LatticePolicy::RequireExact,
    ))
    .expect_err("the doubled map is not an embedding");
    assert!(
        !matches!(error, MatrixSubductionError::SourceMatrixUncertified { .. }),
        "the embedding must be refused before any source matrix claim: {error}"
    );
    assert!(
        matches!(
            error,
            MatrixSubductionError::EmbeddingRejected(_)
                | MatrixSubductionError::LatticePolicyUnsatisfied { .. }
        ),
        "unexpected error variant: {error}"
    );

    // An unknown label under the stated convention.
    let error = verify_matrix_gamma_input(&request(
        89,
        "GM99",
        "GM1",
        0,
        IDENTITY,
        Some(RealSign::Plus),
        LatticePolicy::RequireExact,
    ))
    .expect_err("the label does not resolve");
    assert!(matches!(
        error,
        MatrixSubductionError::LabelNotFound { sg: 89, .. }
    ));
}

#[test]
fn a_record_without_an_in_family_target_is_refused_with_its_examined_list() {
    // UNI 749 (P4'22') has a Type C two-dimensional branch and a spinor-source
    // Type A branch; neither is inside this module's family, so the refusal names
    // the candidates it examined instead of returning a partial success.
    let model = witness_model(89, "GM5", "GM2", 0, IDENTITY, Some(RealSign::Plus));
    assert_eq!(model.uni(), 749);
    let error = subduce_matrix_gamma(&model).expect_err("no in-family target");
    match error {
        MatrixSubductionError::NoMatchingTarget { examined, uni, .. } => {
            assert_eq!(uni, 749);
            assert_eq!(examined, vec!["GM1".to_string(), "GM2".to_string()]);
        }
        other => panic!("expected NoMatchingTarget, got {other}"),
    }
}

#[test]
fn the_full_trigonal_source_passes_its_numerical_matrix_law() {
    // SG 156 GM3 condensed along GM1 gives UNI 1279 (P3m1 itself): a full
    // trigonal unitary group whose matrices multiply the ±√3/2 entries, so the
    // source law is only numerically — not bitwise — equal to the product.  This
    // is the regression that a bitwise `!=` comparison used to refuse.
    for eta in [RealSign::Plus, RealSign::Minus] {
        let model = witness_model(156, "GM3", "GM1", 0, IDENTITY, Some(eta));
        assert_eq!(model.uni(), 1279);
        assert_eq!(model.classes().len(), 6);
        assert!(model.realises_record_lattice());
        assert!(
            model
                .source_operators()
                .iter()
                .all(|operator| !operator.is_antiunitary())
        );
        let report = subduce_matrix_gamma(&model).expect("the full trigonal target reconstructs");
        assert_eq!(report.covered_dimension(), 2);
        assert_eq!(report.pieces().len(), 1);
        let piece = report.piece(0).expect("one piece");
        assert_eq!(piece.label(), "GM3");
        assert_eq!(piece.unitary_sg(), 156);
        assert_eq!(piece.dimension(), 2);
        assert_eq!(piece.multiplicity(), 1);
        assert_isometric_copy(&piece.copies()[0], 2);
        assert_reconstruction(&model, piece);
    }
}

#[test]
fn a_hexagonal_dim2_anti_target_comes_from_another_space_group() {
    // SG 191 GM5+ condensed along GM1- gives UNI 1471: the two-dimensional
    // target GM5 is induced from H = SG 177 (P622), a *different* space group
    // than the source's parent, with twelve antiunitary classes.  Its coset
    // matrix has a non-unit square ratio (measured 1/4), so the rescaling
    // `√ratio · base` — not the inverted `base/√ratio` — is what closes the
    // square; this is the non-echo two-dimensional target witness.
    for eta in [RealSign::Plus, RealSign::Minus] {
        let model = witness_model(191, "GM5+", "GM1-", 0, IDENTITY, Some(eta));
        assert_eq!(model.uni(), 1471);
        assert_eq!(model.classes().len(), 24);
        assert!(model.realises_record_lattice());
        let report = subduce_matrix_gamma(&model).expect("the hexagonal target reconstructs");
        assert_eq!(report.covered_dimension(), 2);
        assert_eq!(report.pieces().len(), 1);
        let piece = report.piece(0).expect("one piece");
        assert_eq!(piece.label(), "GM5");
        assert_eq!(
            piece.unitary_sg(),
            177,
            "the target is induced from SG 177, not from the source's parent"
        );
        assert_eq!(piece.source_label(), "GM5");
        assert_eq!(piece.dimension(), 2);
        assert_eq!(piece.multiplicity(), 1);
        assert_eq!(
            piece
                .operators()
                .iter()
                .filter(|operator| operator.is_antiunitary())
                .count(),
            12
        );
        assert_isometric_copy(&piece.copies()[0], 2);
        assert_reconstruction(&model, piece);
    }
    // The two eta gauges give different intertwiners.
    let plus = witness_model(191, "GM5+", "GM1-", 0, IDENTITY, Some(RealSign::Plus));
    let minus = witness_model(191, "GM5+", "GM1-", 0, IDENTITY, Some(RealSign::Minus));
    let plus_report = subduce_matrix_gamma(&plus).expect("plus");
    let minus_report = subduce_matrix_gamma(&minus).expect("minus");
    let difference: f64 = plus_report
        .piece(0)
        .expect("plus piece")
        .copies()
        .iter()
        .flatten()
        .zip(
            minus_report
                .piece(0)
                .expect("minus piece")
                .copies()
                .iter()
                .flatten(),
        )
        .map(|(left, right)| (*left - *right).norm())
        .fold(0.0f64, f64::max);
    assert!(difference > 0.5, "eta must change the intertwiner");
}
