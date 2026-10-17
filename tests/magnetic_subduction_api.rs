//! Public regressions for the R30 high-level Γ magnetic subduction API.
//!
//! These tests exercise the public request/answer surface only.  The engine
//! layers R25/R26/R29 keep their own suites; here the point is that the new
//! entry point resolves the documented records, dispatches by source dimension,
//! reports source/condensing/record/target identities with real label sources,
//! and hands back payloads that reconstruct the operators on complex vectors.

use cryspglib::irrep::magnetic_embedding::{canonical_translation, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, MagneticInputError, RealSign, RecordIdentity,
    ScalarGammaRequest, SourceIdentity, verify_scalar_gamma_input,
};
use cryspglib::irrep::magnetic_matrix_subduction::MatrixGammaModel;
use cryspglib::irrep::magnetic_subduction::subduce_scalar_gamma;
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticOutcome, MagneticRecordSelector, MagneticSubductionApiError, MagneticSubductionReport,
    MagneticSubductionRequest, subduce,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, query};
use num_complex::Complex64;

const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const ANTI_TRANSLATION: [[i32; 3]; 3] = [[-1, -1, 0], [-1, 1, 0], [0, 0, -1]];
const TWO_I: [[i32; 3]; 3] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]];

fn request(
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    direction: cryspglib::irrep::isotropy::IsotropyDirection<'static>,
    uni: Option<usize>,
    eta: Option<RealSign>,
    embedding: EmbeddingEvidence,
) -> MagneticSubductionRequest<'static> {
    let mut record =
        MagneticRecordSelector::new(sg, condensing, LabelConvention::Cdml).direction(direction);
    if let Some(uni) = uni {
        record = record.uni(uni);
    }
    MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention: LabelConvention::Cdml,
        },
        record,
        eta,
    )
    .with_embedding(embedding)
}

fn embedding(map: [[i32; 3]; 3], policy: LatticePolicy) -> EmbeddingEvidence {
    EmbeddingEvidence {
        map: Mat3R::from_ints(map),
        shift: [Rat::ZERO; 3],
        lattice_policy: policy,
    }
}

fn p1() -> cryspglib::irrep::isotropy::IsotropyDirection<'static> {
    cryspglib::irrep::isotropy::IsotropyDirection::Label("P1")
}

fn close(left: Complex64, right: Complex64) -> bool {
    (left - right).norm() < 1e-12
}

/// Scalar relation `D(g)∘U = U∘A(g)` for every class, plus the full operator
/// product table.  Mirrors the independent R26 contract check.
fn assert_scalar_reconstruction(report: &MagneticSubductionReport) {
    let model = report.scalar_model().expect("a scalar model");
    let outcome = report.scalar().expect("a scalar payload");
    assert_eq!(outcome.operators.len(), model.actions().len());
    assert_eq!(outcome.operators.len(), report.classes().len());
    let u = outcome.intertwiner;
    assert!((u.norm_sqr() - 1.0).abs() < 1e-12, "U must be unimodular");
    for (index, class) in report.classes().iter().enumerate() {
        let operator = &outcome.operators[index];
        assert_eq!(
            operator.is_antiunitary(),
            class.time_reversal,
            "the antiunitary flag is the class's"
        );
        for z in [Complex64::new(1.0, 2.0), Complex64::new(-0.7, 0.4)] {
            let v = u * z;
            let source = f64::from(model.actions()[index].value.value())
                * if class.time_reversal { v.conj() } else { v };
            let target = u * operator.apply(z);
            assert!(
                close(source, target),
                "class {index}: source action {source:?} != U·D(g)·z {target:?}"
            );
        }
        for (other_index, other) in report.classes().iter().enumerate() {
            let product = compose(class, other).expect("class product");
            let product_index = report
                .classes()
                .iter()
                .position(|candidate| {
                    candidate.rotation == product.rotation
                        && candidate.time_reversal == product.time_reversal
                        && canonical_translation(
                            candidate.translation,
                            report.translation_lattice(),
                        ) == canonical_translation(
                            product.translation,
                            report.translation_lattice(),
                        )
                })
                .expect("complete class list");
            assert_eq!(
                outcome.operators[index].compose(&outcome.operators[other_index]),
                outcome.operators[product_index],
                "operator product ({index}, {other_index})"
            );
        }
    }
}

fn act(entries: &[f64], dimension: usize, antiunitary: bool, z: &[Complex64]) -> Vec<Complex64> {
    (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| {
                    entries[row * dimension + column]
                        * if antiunitary {
                            z[column].conj()
                        } else {
                            z[column]
                        }
                })
                .sum()
        })
        .collect()
}

fn lift(copy: &[Complex64], dimension: usize, z: &[Complex64]) -> Vec<Complex64> {
    (0..2)
        .map(|row| {
            (0..dimension)
                .map(|column| copy[row * dimension + column] * z[column])
                .sum()
        })
        .collect()
}

/// `V†V = I` and each copy's `S(g)∘φ = φ∘T(g)` on two complex probes.
fn assert_matrix_reconstruction(model: &MatrixGammaModel, report: &MagneticSubductionReport) {
    let outcome = report.matrix().expect("a matrix payload");
    assert_eq!(outcome.source_dimension, 2);
    assert_eq!(outcome.covered_dimension, 2);
    let total: usize = outcome
        .pieces
        .iter()
        .map(|piece| piece.target.multiplicity as usize * piece.target.dimension as usize)
        .sum();
    assert_eq!(total, 2, "the pieces must cover the whole source space");
    for piece in &outcome.pieces {
        let dimension = piece.target.dimension as usize;
        assert_eq!(piece.copies.len(), piece.target.multiplicity as usize);
        assert_eq!(piece.operators.len(), model.source_operators().len());
        for copy in &piece.copies {
            assert_eq!(copy.len(), 2 * dimension);
            for row in 0..dimension {
                for column in 0..dimension {
                    let entry: Complex64 = (0..2)
                        .map(|k| copy[k * dimension + row].conj() * copy[k * dimension + column])
                        .sum();
                    let expected = if row == column {
                        Complex64::new(1.0, 0.0)
                    } else {
                        Complex64::new(0.0, 0.0)
                    };
                    assert!(close(entry, expected), "V†V entry ({row}, {column})");
                }
            }
            for (class, source) in model.source_operators().iter().enumerate() {
                let target = &piece.operators[class];
                assert_eq!(target.is_antiunitary(), source.is_antiunitary());
                for z in [
                    vec![Complex64::new(0.31, 0.73); dimension],
                    (0..dimension)
                        .map(|i| Complex64::new(-0.51 + i as f64, -0.29))
                        .collect::<Vec<_>>(),
                ] {
                    let left = act(
                        &source.entries(),
                        2,
                        source.is_antiunitary(),
                        &lift(copy, dimension, &z),
                    );
                    let right = lift(
                        copy,
                        dimension,
                        &act(target.entries(), dimension, target.is_antiunitary(), &z),
                    );
                    assert!(left.iter().zip(&right).all(|(&a, b)| close(a, *b)));
                }
            }
        }
        // Independent copies are orthogonal in the `2×d` cross block.
        for (left_index, left) in piece.copies.iter().enumerate() {
            for right in piece.copies.iter().skip(left_index + 1) {
                for row in 0..dimension {
                    for column in 0..dimension {
                        let block: Complex64 = (0..2)
                            .map(|k| {
                                left[k * dimension + row].conj() * right[k * dimension + column]
                            })
                            .sum();
                        assert!(close(block, Complex64::new(0.0, 0.0)));
                    }
                }
            }
        }
    }
}

/// Mirror the request from the high-level report into the R26 engine and assert
/// that the payload it returns is exactly the one the wrapper reports.
fn assert_scalar_payload_matches_engine(report: &MagneticSubductionReport, eta: RealSign) {
    let sg = report.source().parent_sg;
    let source = report.source().labels.cdml;
    let condensing = report.condensing().labels.cdml;
    let index = report.record().direction_index;
    let model = verify_scalar_gamma_input(
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
            Some(eta),
        )
        .with_embedding((*report.embedding()).into_evidence()),
    )
    .expect("the engine verifies the same row");
    let engine = subduce_scalar_gamma(&model).expect("the engine answers");
    let outcome = report.scalar().expect("scalar payload");
    assert_eq!(outcome.target.magnetic_label, engine.target_label());
    assert_eq!(outcome.target.h_source.cdml, engine.target_source_label());
    assert_eq!(outcome.target.unitary_sg, engine.target_unitary_sg());
    assert_eq!(outcome.target.dimension, engine.target_dimension());
    assert_eq!(outcome.target.multiplicity, engine.multiplicity());
    assert_eq!(outcome.intertwiner, engine.intertwiner());
    assert_eq!(outcome.operators, engine.operators());
}

/// Small convenience so the engine-mirroring helper can rebuild an
/// `EmbeddingEvidence` from the public report identity.
trait EmbeddingEvidenceFromIdentity {
    fn into_evidence(self) -> EmbeddingEvidence;
}

impl EmbeddingEvidenceFromIdentity
    for cryspglib::irrep::magnetic_subduction_api::MagneticEmbeddingIdentity
{
    fn into_evidence(self) -> EmbeddingEvidence {
        EmbeddingEvidence {
            map: self.map,
            shift: self.shift,
            lattice_policy: self.lattice_policy,
        }
    }
}

#[test]
fn scalar_uni22_keeps_source_record_and_target_identities_apart() {
    // SG 5 GM1 (CDML Γ1) subduced with the GM2 (Γ2) records condensing into
    // UNI 22, whose unitary subgroup is P1.  Source ≠ condensing, H ≠ parent.
    let plus = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        Some(22),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("the pinned SG 5 witness");
    assert!(matches!(plus.outcome(), MagneticOutcome::Scalar(_)));
    assert!(plus.matrix().is_none());
    assert_eq!(plus.source().labels.cdml, "GM1");
    assert_eq!(plus.source().labels.bc.as_deref(), Some("Γ1"));
    assert_eq!(plus.condensing().labels.cdml, "GM2");
    assert_eq!(plus.condensing().labels.bc.as_deref(), Some("Γ2"));
    assert_ne!(
        plus.source().labels.cdml,
        plus.condensing().labels.cdml,
        "the source and the condensing record are distinct requests"
    );
    assert_eq!(plus.record().uni, 22);
    assert_eq!(plus.record().direction_label, "P1");
    assert_eq!(plus.record().direction_index, 0);
    assert_eq!(plus.source_dimension(), 1);
    assert_eq!(plus.covered_dimension(), 1);
    let target = &plus.scalar().expect("scalar answer").target;
    assert_eq!(target.magnetic_label, "GM1");
    assert!(target.magnetic_bc.is_none(), "no magnetic BC source exists");
    assert_eq!(target.h_source.cdml, "GM1");
    assert_eq!(target.h_source.bc.as_deref(), Some("Γ1"));
    assert_eq!(target.unitary_sg, 1);
    assert_ne!(
        u32::from(plus.condensing().parent_sg),
        u32::from(target.unitary_sg)
    );
    assert_eq!(target.dimension, 1);
    assert_eq!(target.multiplicity, 1);
    assert_scalar_reconstruction(&plus);
    assert_scalar_payload_matches_engine(&plus, RealSign::Plus);

    // η is a gauge of the antiunitary coset: same target table, different U.
    let minus = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        Some(22),
        Some(RealSign::Minus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("the η = -1 witness");
    let plus_payload = plus.scalar().expect("scalar");
    let minus_payload = minus.scalar().expect("scalar");
    assert_eq!(plus_payload.target, minus_payload.target);
    assert_eq!(plus_payload.operators, minus_payload.operators);
    assert_ne!(plus_payload.intertwiner, minus_payload.intertwiner);
    assert!(close(plus_payload.intertwiner, Complex64::new(1.0, 0.0)));
    assert!(close(minus_payload.intertwiner, Complex64::new(0.0, 1.0)));
    assert_scalar_reconstruction(&minus);
}

#[test]
fn scalar_uni59_nontrivial_source_has_an_independent_target_label() {
    // SG 11 GM1- (Γ1-) condenses with the GM2- (Γ2-) record into UNI 59, whose
    // unitary subgroup is SG 6.  The target is GM2 of H = SG 6, not an echo of
    // the GM1- source label; its BC Γ2 comes from the H irrep, not the corep.
    let report = subduce(&request(
        11,
        "GM1-",
        "GM2-",
        p1(),
        Some(59),
        Some(RealSign::Minus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("the pinned SG 11 GM1- witness");
    assert_eq!(report.source().labels.cdml, "GM1-");
    assert_eq!(report.source().labels.bc.as_deref(), Some("Γ1-"));
    assert_eq!(report.condensing().labels.cdml, "GM2-");
    assert_eq!(report.condensing().labels.bc.as_deref(), Some("Γ2-"));
    assert_eq!(report.record().uni, 59);
    let target = &report.scalar().expect("scalar").target;
    assert_eq!(target.magnetic_label, "GM2");
    assert_ne!(target.magnetic_label, report.source().labels.cdml);
    assert!(target.magnetic_bc.is_none());
    assert_eq!(target.h_source.cdml, "GM2");
    assert_eq!(target.h_source.bc.as_deref(), Some("Γ2"));
    assert_eq!(target.unitary_sg, 6);
    assert_ne!(target.unitary_sg, report.condensing().parent_sg);
    assert_ne!(report.source().labels.cdml, report.condensing().labels.cdml);
    assert_eq!(report.covered_dimension(), 1);
    assert_scalar_reconstruction(&report);
    assert_scalar_payload_matches_engine(&report, RealSign::Minus);

    // The same record with the GM1+ source keeps GM1 as the target: the record
    // alone does not determine the target, the source does.
    let plus_source = subduce(&request(
        11,
        "GM1+",
        "GM2-",
        p1(),
        Some(59),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("GM1+ source");
    assert_eq!(
        plus_source.scalar().expect("scalar").target.magnetic_label,
        "GM1"
    );
}

#[test]
fn matrix_sg191_uses_the_real_bc_of_both_ordinary_identities() {
    // SG 191 GM5+ (BC Γ6+) condenses with GM1- (BC Γ1-) row 0 into UNI 1471,
    // whose unitary subgroup is SG 177 with irrep GM5 (BC Γ5).  This is the
    // discriminative BC witness: replacing any genuine BC spelling by its CDML
    // string would change the asserted values.
    let source = query::find_irreps(191, "GM5+", LabelConvention::Cdml);
    assert_eq!(source.len(), 1);
    assert_eq!(source[0].label(LabelConvention::Bc).as_deref(), Some("Γ6+"));
    let h_source = query::find_irreps(177, "GM5", LabelConvention::Cdml);
    assert_eq!(h_source.len(), 1);
    assert_eq!(
        h_source[0].label(LabelConvention::Bc).as_deref(),
        Some("Γ5")
    );

    for eta in [RealSign::Plus, RealSign::Minus] {
        let report = subduce(&request(
            191,
            "GM5+",
            "GM1-",
            p1(),
            Some(1471),
            Some(eta),
            embedding(IDENTITY, LatticePolicy::ContainmentOnly),
        ))
        .expect("the pinned SG 191 witness");
        assert!(matches!(report.outcome(), MagneticOutcome::Matrix(_)));
        assert!(report.scalar().is_none());
        assert_eq!(report.source().labels.cdml, "GM5+");
        assert_eq!(report.source().labels.bc.as_deref(), Some("Γ6+"));
        assert_eq!(report.condensing().labels.cdml, "GM1-");
        assert_eq!(report.condensing().labels.bc.as_deref(), Some("Γ1-"));
        assert_eq!(report.record().uni, 1471);
        assert_eq!(report.record().direction_index, 0);
        assert_eq!(report.source_dimension(), 2);
        let outcome = report.matrix().expect("matrix payload");
        assert_eq!(outcome.pieces.len(), 1);
        let piece = &outcome.pieces[0];
        assert_eq!(piece.target.magnetic_label, "GM5");
        assert!(piece.target.magnetic_bc.is_none());
        assert_eq!(piece.target.h_source.cdml, "GM5");
        assert_eq!(piece.target.h_source.bc.as_deref(), Some("Γ5"));
        assert_eq!(piece.target.unitary_sg, 177);
        assert_ne!(piece.target.unitary_sg, report.condensing().parent_sg);
        assert_eq!(piece.target.dimension, 2);
        assert_eq!(piece.target.multiplicity, 1);
        assert_eq!(report.covered_dimension(), 2);
        let anti = piece
            .operators
            .iter()
            .filter(|operator| operator.is_antiunitary())
            .count();
        assert_eq!(anti, 12, "the target carries the antiunitary coset");
        let model = report.matrix_model().expect("matrix model");
        assert_eq!(model.source_label(), "GM5+");
        assert_eq!(model.record_label(), "GM1-");
        assert_matrix_reconstruction(model, &report);
    }
}

#[test]
fn matrix_type_i_split_and_multiplicity_witnesses_are_covered() {
    // Type I ordinary: SG 89 GM5 / GM1 is the parent itself (UNI 747).
    let type_i = subduce(&request(
        89,
        "GM5",
        "GM1",
        p1(),
        Some(747),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("UNI 747");
    let piece = &type_i.matrix().expect("matrix").pieces[0];
    assert_eq!(piece.target.h_source.cdml, "GM5");
    assert_eq!(piece.target.unitary_sg, 89);
    assert_eq!(piece.target.multiplicity, 1);
    assert_matrix_reconstruction(type_i.matrix_model().unwrap(), &type_i);

    // 2 → 1 + 1 split: SG 89 GM5 / GM5 row 0 is UNI 101 with H = SG 3.
    let split = subduce(&request(
        89,
        "GM5",
        "GM5",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(0),
        Some(101),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("UNI 101");
    let outcome = split.matrix().expect("matrix");
    assert_eq!(outcome.pieces.len(), 2);
    let mut labels: Vec<&str> = outcome
        .pieces
        .iter()
        .map(|piece| piece.target.magnetic_label.as_str())
        .collect();
    labels.sort_unstable();
    assert_eq!(labels, vec!["GM1", "GM2"]);
    for piece in &outcome.pieces {
        assert_eq!(piece.target.unitary_sg, 3);
        assert_eq!(piece.target.dimension, 1);
        assert_eq!(piece.target.multiplicity, 1);
        assert!(piece.target.magnetic_bc.is_none());
    }
    assert_matrix_reconstruction(split.matrix_model().unwrap(), &split);

    // Multiplicity two: SG 89 GM5 / GM5 row 2 is UNI 10 with H = SG 1.
    let mult = subduce(&request(
        89,
        "GM5",
        "GM5",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(2),
        Some(10),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("UNI 10");
    let piece = &mult.matrix().expect("matrix").pieces[0];
    assert_eq!(piece.target.h_source.cdml, "GM1");
    assert_eq!(piece.target.unitary_sg, 1);
    assert_eq!(piece.target.dimension, 1);
    assert_eq!(piece.target.multiplicity, 2);
    assert_eq!(piece.copies.len(), 2);
    assert_matrix_reconstruction(mult.matrix_model().unwrap(), &mult);

    // The same 2 × 1 shape across a different parent: SG 156 GM3 / GM3 row 2.
    let cross = subduce(&request(
        156,
        "GM3",
        "GM3",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(2),
        Some(1),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("UNI 1");
    let piece = &cross.matrix().expect("matrix").pieces[0];
    assert_eq!(piece.target.h_source.cdml, "GM1");
    assert_eq!(piece.target.unitary_sg, 1);
    assert_eq!(piece.target.multiplicity, 2);
    assert_matrix_reconstruction(cross.matrix_model().unwrap(), &cross);

    // The anti-translation target needs the pinned non-identity map; both η
    // signs give the same target.
    let anti_plus = subduce(&request(
        89,
        "GM5",
        "M1",
        p1(),
        Some(753),
        Some(RealSign::Plus),
        embedding(ANTI_TRANSLATION, LatticePolicy::ContainmentOnly),
    ))
    .expect("UNI 753 with its realising map");
    let anti_minus = subduce(&request(
        89,
        "GM5",
        "M1",
        p1(),
        Some(753),
        Some(RealSign::Minus),
        embedding(ANTI_TRANSLATION, LatticePolicy::ContainmentOnly),
    ))
    .expect("UNI 753, η = -1");
    let plus_piece = &anti_plus.matrix().expect("matrix").pieces[0];
    let minus_piece = &anti_minus.matrix().expect("matrix").pieces[0];
    assert_eq!(plus_piece.target, minus_piece.target);
    assert_eq!(plus_piece.operators, minus_piece.operators);
    assert_ne!(plus_piece.copies, minus_piece.copies);
    assert_eq!(plus_piece.target.unitary_sg, 89);
    assert_matrix_reconstruction(anti_plus.matrix_model().unwrap(), &anti_plus);
    assert_matrix_reconstruction(anti_minus.matrix_model().unwrap(), &anti_minus);
}

#[test]
fn direction_and_uni_selectors_agree_and_conflicts_are_named() {
    // Label, index and UNI select the same row of the same condensing irrep.
    let by_label = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("by label");
    let by_index = subduce(&request(
        5,
        "GM1",
        "GM2",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(0),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("by index");
    let by_uni = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        Some(22),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("by UNI");
    assert_eq!(by_label.record(), by_index.record());
    assert_eq!(by_label.record(), by_uni.record());
    assert_eq!(
        by_label.scalar().expect("scalar").target,
        by_uni.scalar().expect("scalar").target
    );

    // A descriptor carrying the magnetic row's own direction field follows the
    // layer's existing matching and selects the same row; a real component
    // string has no magnetic source.
    let by_descriptor = subduce(&request(
        5,
        "GM1",
        "GM2",
        cryspglib::irrep::isotropy::IsotropyDirection::Descriptor("P1"),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("the magnetic direction field is P1");
    assert_eq!(by_descriptor.record(), by_label.record());
    let component = subduce(&request(
        5,
        "GM1",
        "GM2",
        cryspglib::irrep::isotropy::IsotropyDirection::Descriptor("(a,0,0)"),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("a component descriptor has no magnetic source");
    assert!(matches!(
        component,
        MagneticSubductionApiError::Isotropy(
            cryspglib::irrep::isotropy::IsotropyError::DirectionNotFound { .. }
        )
    ));

    // The official-descriptor variant is not invented for magnetic rows.
    let official = subduce(&request(
        5,
        "GM1",
        "GM2",
        cryspglib::irrep::isotropy::IsotropyDirection::OfficialDescriptor("P1"),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("magnetic rows have no official descriptor table");
    assert!(matches!(
        official,
        MagneticSubductionApiError::Isotropy(
            cryspglib::irrep::isotropy::IsotropyError::NoOfficialDescriptor { .. }
        )
    ));

    // A direction and a UNI that name different rows are a real conflict.
    let mismatch = subduce(&request(
        89,
        "GM5",
        "GM5",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(0),
        Some(753),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("row 0 is UNI 101, not 753");
    assert!(matches!(
        mismatch,
        MagneticSubductionApiError::RecordUniMismatch {
            requested: 753,
            found: 101,
            ..
        }
    ));

    // A missing selector, an unknown UNI, an out-of-range index and an unknown
    // direction all keep their own names.
    let missing = subduce(&MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: "GM5",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(89, "GM5", LabelConvention::Cdml),
        Some(RealSign::Plus),
    ));
    assert!(matches!(
        missing,
        Err(MagneticSubductionApiError::RecordSelectorMissing { sg: 89 })
    ));
    let unknown_uni = subduce(&MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 5,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(5, "GM2", LabelConvention::Cdml).uni(59),
        Some(RealSign::Plus),
    ))
    .expect_err("SG 5 GM2 has no UNI 59 row");
    assert!(matches!(
        unknown_uni,
        MagneticSubductionApiError::RecordUniNotFound { uni: 59, .. }
    ));
    let out_of_range = subduce(&request(
        89,
        "GM5",
        "GM5",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(99),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("the table has three rows");
    assert!(matches!(
        out_of_range,
        MagneticSubductionApiError::Isotropy(
            cryspglib::irrep::isotropy::IsotropyError::SubgroupIndexOutOfRange { .. }
        )
    ));
    let unknown_direction = subduce(&request(
        89,
        "GM5",
        "GM5",
        cryspglib::irrep::isotropy::IsotropyDirection::Label("NOPE"),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("no such direction");
    assert!(matches!(
        unknown_direction,
        MagneticSubductionApiError::Isotropy(
            cryspglib::irrep::isotropy::IsotropyError::DirectionNotFound { .. }
        )
    ));
}

#[test]
fn a_uni_number_listed_twice_is_refused_as_ambiguous() {
    // SG 5 V1 lists UNI 3 at two directions (P3 and C1), so the UNI number
    // alone does not name one record.
    let ambiguous = subduce(&MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 5,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(5, "V1", LabelConvention::Cdml).uni(3),
        Some(RealSign::Plus),
    ))
    .expect_err("UNI 3 alone is ambiguous");
    match ambiguous {
        MagneticSubductionApiError::AmbiguousRecordUni {
            sg,
            uni,
            candidates,
        } => {
            assert_eq!(sg, 5);
            assert_eq!(uni, 3);
            assert_eq!(candidates.len(), 2);
            assert_eq!(candidates[0].0, 1);
            assert_eq!(candidates[1].0, 2);
        }
        other => panic!("expected AmbiguousRecordUni, got {other:?}"),
    }
}

#[test]
fn mixed_label_conventions_round_trip_to_the_same_identity() {
    // Source and condensing are independent requests and may use different
    // conventions.  Read the genuine BC spellings from the database and check
    // that the report's canonical identities equal the CDML request's.
    let source_bc = query::find_irreps(5, "GM1", LabelConvention::Cdml)[0]
        .label(LabelConvention::Bc)
        .expect("SG 5 GM1 has a BC label");
    let condensing_bc = query::find_irreps(5, "GM2", LabelConvention::Cdml)[0]
        .label(LabelConvention::Bc)
        .expect("SG 5 GM2 has a BC label");
    let mixed = subduce(
        &MagneticSubductionRequest::new(
            SourceIdentity {
                parent_sg: 5,
                label: &source_bc,
                convention: LabelConvention::Bc,
            },
            MagneticRecordSelector::new(5, &condensing_bc, LabelConvention::Bc)
                .direction(p1())
                .uni(22),
            Some(RealSign::Plus),
        )
        .with_embedding(EmbeddingEvidence::identity()),
    )
    .expect("mixed BC request");
    assert_eq!(mixed.source().labels.cdml, "GM1");
    assert_eq!(
        mixed.source().labels.bc.as_deref(),
        Some(source_bc.as_str())
    );
    assert_eq!(mixed.source().input_convention, LabelConvention::Bc);
    assert_eq!(mixed.condensing().labels.cdml, "GM2");
    assert_eq!(
        mixed.condensing().labels.bc.as_deref(),
        Some(condensing_bc.as_str())
    );
    assert_eq!(mixed.condensing().input_convention, LabelConvention::Bc);
    assert_eq!(mixed.record().uni, 22);
    let cdml = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        Some(22),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect("CDML request");
    assert_eq!(mixed.source().labels, cdml.source().labels);
    assert_eq!(mixed.source().dimension, cdml.source().dimension);
    assert_eq!(mixed.source().spinor, cdml.source().spinor);
    assert_eq!(mixed.condensing().labels, cdml.condensing().labels);
    assert_eq!(mixed.condensing().dimension, cdml.condensing().dimension);
    assert_eq!(mixed.record(), cdml.record());
    assert_eq!(
        mixed.scalar().expect("scalar").target,
        cdml.scalar().expect("scalar").target
    );

    // The matrix witness keeps the same property across its two BC labels, and
    // its H source keeps its own genuine BC Γ5.
    let source_bc = query::find_irreps(191, "GM5+", LabelConvention::Cdml)[0]
        .label(LabelConvention::Bc)
        .expect("SG 191 GM5+ BC");
    let condensing_bc = query::find_irreps(191, "GM1-", LabelConvention::Cdml)[0]
        .label(LabelConvention::Bc)
        .expect("SG 191 GM1- BC");
    assert_eq!(source_bc, "Γ6+");
    assert_eq!(condensing_bc, "Γ1-");
    let mixed = subduce(
        &MagneticSubductionRequest::new(
            SourceIdentity {
                parent_sg: 191,
                label: &source_bc,
                convention: LabelConvention::Bc,
            },
            MagneticRecordSelector::new(191, &condensing_bc, LabelConvention::Bc)
                .direction(p1())
                .uni(1471),
            Some(RealSign::Plus),
        )
        .with_embedding(EmbeddingEvidence::identity()),
    )
    .expect("mixed BC matrix request");
    assert_eq!(mixed.source().labels.cdml, "GM5+");
    assert_eq!(mixed.condensing().labels.cdml, "GM1-");
    let piece = &mixed.matrix().expect("matrix").pieces[0];
    assert_eq!(piece.target.h_source.cdml, "GM5");
    assert_eq!(piece.target.h_source.bc.as_deref(), Some("Γ5"));
    assert!(piece.target.magnetic_bc.is_none());
}

#[test]
fn engine_errors_and_source_domain_are_preserved() {
    // A three-dimensional source is outside both supported entries.
    let dimension = subduce(&request(
        221,
        "GM4+",
        "GM1",
        p1(),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("SG 221 GM4+ is three-dimensional");
    assert!(matches!(
        dimension,
        MagneticSubductionApiError::SourceDimensionUnsupported {
            sg: 221,
            dimension: 3,
            ..
        }
    ));

    // A non-Γ scalar source reaches R25 and keeps its own named refusal.
    let not_gamma = subduce(&request(
        1,
        "X1",
        "GM1",
        p1(),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("SG 1 X1 is not Γ");
    match not_gamma {
        MagneticSubductionApiError::ScalarInput(error) => {
            assert!(matches!(*error, MagneticInputError::SourceNotGamma { .. }))
        }
        other => panic!("expected the R25 SourceNotGamma, got {other:?}"),
    }

    // A spinor source is refused by the scalar contract, not silently dropped.
    let spinor = subduce(&request(
        1,
        "GM2",
        "GM1",
        p1(),
        None,
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ));
    // SG 1 has a spinor GM2 but no scalar GM1 condensing record; resolve the
    // source first to make sure the refusal is about the source domain.
    assert!(
        spinor.is_err(),
        "a spinor source must not produce an answer"
    );
    let source = query::irreps_of(1)
        .iter()
        .find(|record| record.spinor && record.ml == "GM2")
        .expect("SG 1 spinor GM2");
    let spinor = subduce(&request(
        1,
        source.ml,
        "X1",
        p1(),
        Some(3),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("spinor source");
    match spinor {
        MagneticSubductionApiError::ScalarInput(error) => {
            assert!(matches!(*error, MagneticInputError::SpinorSource { .. }))
        }
        other => panic!("expected SpinorSource, got {other:?}"),
    }

    // An antiunitary group without η is refused by the input contract.
    let missing_eta = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        Some(22),
        None,
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("UNI 22 has an antiunitary class");
    match missing_eta {
        MagneticSubductionApiError::ScalarInput(error) => {
            assert!(matches!(
                *error,
                MagneticInputError::MissingTimeReversal { uni: 22 }
            ))
        }
        other => panic!("expected MissingTimeReversal, got {other:?}"),
    }

    // The wrong map is a layer error, not a silently ignored embedding.
    let wrong_map = subduce(&request(
        89,
        "GM5",
        "M1",
        p1(),
        Some(753),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ))
    .expect_err("identity does not embed UNI 753");
    assert!(matches!(wrong_map, MagneticSubductionApiError::Matrix(_)));

    // RequireExact recomputes the record lattice; SG 1 X1 with 2I is contained
    // but does not realise the lattice exactly.
    let exact = subduce(&request(
        1,
        "GM1",
        "X1",
        cryspglib::irrep::isotropy::IsotropyDirection::Index(0),
        Some(3),
        Some(RealSign::Plus),
        embedding(TWO_I, LatticePolicy::RequireExact),
    ))
    .expect_err("2I is a proper superlattice embedding here");
    match exact {
        MagneticSubductionApiError::ScalarInput(error) => assert!(matches!(
            *error,
            MagneticInputError::LatticePolicyUnsatisfied {
                policy: LatticePolicy::RequireExact
            }
        )),
        other => panic!("expected LatticePolicyUnsatisfied, got {other:?}"),
    }

    // Unknown labels and parent mismatches are request errors.
    let unknown = subduce(&request(
        5,
        "NOPE",
        "GM2",
        p1(),
        Some(22),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ));
    assert!(matches!(
        unknown,
        Err(MagneticSubductionApiError::SourceLabelNotFound { sg: 5, .. })
    ));
    let mismatch = subduce(&request(
        5,
        "GM1",
        "GM2",
        p1(),
        Some(22),
        Some(RealSign::Plus),
        embedding(IDENTITY, LatticePolicy::ContainmentOnly),
    ));
    assert!(mismatch.is_ok());
    let mismatch = subduce(&MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 5,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(89, "GM5", LabelConvention::Cdml)
            .direction(p1())
            .uni(101),
        Some(RealSign::Plus),
    ));
    assert!(matches!(
        mismatch,
        Err(MagneticSubductionApiError::ParentMismatch {
            source_sg: 5,
            record_sg: 89
        })
    ));
}
