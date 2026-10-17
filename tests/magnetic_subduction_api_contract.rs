//! Leader-owned independent public-boundary checks; no production apply helpers.
use cryspglib::irrep::LabelConvention::{Bc, Cdml};
use cryspglib::irrep::isotropy::{self, IsotropyDirection, IsotropyError};
use cryspglib::irrep::magnetic_input::{EmbeddingEvidence, RealSign, SourceIdentity};
use cryspglib::irrep::magnetic_matrix_subduction::subduce_matrix_gamma;
use cryspglib::irrep::magnetic_subduction::subduce_scalar_gamma;
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticRecordSelector, MagneticSubductionApiError as ApiError, MagneticSubductionRequest,
    subduce,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::Mat3R;
use num_complex::Complex64 as C;

fn irrep(sg: u8, label: &str) -> Option<&'static cryspglib::irrep::types::IrrepRecord> {
    match query::find_irreps(sg, label, Cdml).as_slice() {
        [record] => Some(*record),
        _ => None,
    }
}

fn close(actual: C, expected: C) {
    assert!(actual.re.is_finite() && actual.im.is_finite());
    assert!(
        (actual - expected).norm() < 1e-8,
        "{actual:?} != {expected:?}"
    );
}

fn dense(matrix: &[f64], d: usize, anti: bool, z: &[C]) -> Vec<C> {
    (0..d)
        .map(|r| {
            (0..d)
                .map(|c| matrix[r * d + c] * if anti { z[c].conj() } else { z[c] })
                .sum()
        })
        .collect()
}

fn lift(u: &[C], d: usize, z: &[C]) -> [C; 2] {
    std::array::from_fn(|r| (0..d).map(|c| u[r * d + c] * z[c]).sum())
}

#[test]
fn scalar_payload_survives_request_and_keeps_independent_identities() {
    for (sg, source, condensing, uni, h, target) in [
        (5, "GM1", "GM2", 22, 1, "GM1"),
        (11, "GM1-", "GM2-", 59, 6, "GM2"),
    ] {
        for eta in [RealSign::Plus, RealSign::Minus] {
            // Owned answer must outlive both request and its input strings.
            let report = {
                let source_name = source.to_string();
                let condensing_bc = irrep(sg, condensing).unwrap().labels().bc.unwrap();
                let request = MagneticSubductionRequest::new(
                    SourceIdentity {
                        parent_sg: sg,
                        label: &source_name,
                        convention: Cdml,
                    },
                    MagneticRecordSelector::new(sg, &condensing_bc, Bc).uni(uni),
                    Some(eta),
                );
                subduce(&request).unwrap()
            };
            assert_eq!(report.source().labels, irrep(sg, source).unwrap().labels());
            assert_eq!(
                report.condensing().labels,
                irrep(sg, condensing).unwrap().labels()
            );
            assert_eq!(report.source().input_convention, Cdml);
            assert_eq!(report.condensing().input_convention, Bc);
            assert_eq!(report.record().uni, uni);
            assert_eq!(report.time_reversal(), Some(eta));
            assert_eq!(report.source_dimension(), 1);
            assert_eq!(report.covered_dimension(), 1);
            let model = report.scalar_model().unwrap();
            let scalar = report.scalar().unwrap();
            let direct = subduce_scalar_gamma(model).unwrap();
            assert_eq!(scalar.intertwiner, direct.intertwiner());
            assert_eq!(scalar.operators, direct.operators());
            assert_eq!(scalar.target.magnetic_label, target);
            assert_eq!(scalar.target.unitary_sg, h);
            assert_eq!(scalar.target.h_source, irrep(h, target).unwrap().labels());
            assert!(scalar.target.magnetic_bc.is_none());
            assert_eq!(report.record().ordinal, model.ordinal());
            assert_eq!(report.embedding().map, *model.map());
            assert_eq!(report.embedding().shift, *model.shift());
            assert_eq!(scalar.operators.len(), model.actions().len());
            close(
                scalar.intertwiner.conj() * scalar.intertwiner,
                C::new(1., 0.),
            );
            for (action, op) in model.actions().iter().zip(&scalar.operators) {
                assert_eq!(action.class.time_reversal, op.is_antiunitary());
                for z in [C::new(0.25, -0.7), C::new(-1.2, 0.6)] {
                    let image = scalar.intertwiner * z;
                    let lhs = f64::from(action.value.value())
                        * if action.class.time_reversal {
                            image.conj()
                        } else {
                            image
                        };
                    let rhs = scalar.intertwiner
                        * op.value()
                        * if op.is_antiunitary() { z.conj() } else { z };
                    close(lhs, rhs);
                }
            }
        }
    }
}

#[test]
fn matrix_payload_has_complete_copies_and_reconstructs_complex_vectors() {
    for (sg, source, condensing, row, uni, h, count, multiplicity) in [
        (89, "GM5", "GM1", 0, 747, 89, 1, 1),
        (89, "GM5", "M1", 0, 753, 89, 1, 1),
        (89, "GM5", "GM5", 0, 101, 3, 2, 1),
        (89, "GM5", "GM5", 2, 10, 1, 1, 2),
        (156, "GM3", "GM3", 2, 1, 1, 1, 2),
        (156, "GM3", "GM1", 0, 1279, 156, 1, 1),
        (191, "GM5+", "GM1-", 0, 1471, 177, 1, 1),
    ] {
        for eta in [RealSign::Plus, RealSign::Minus] {
            let source_label = if sg == 191 { "Γ6+" } else { source };
            let convention = if sg == 191 { Bc } else { Cdml };
            let embedding = if uni == 753 {
                EmbeddingEvidence {
                    map: Mat3R::from_ints([[-1, -1, 0], [-1, 1, 0], [0, 0, -1]]),
                    ..EmbeddingEvidence::identity()
                }
            } else {
                EmbeddingEvidence::identity()
            };
            let request = MagneticSubductionRequest::new(
                SourceIdentity {
                    parent_sg: sg,
                    label: source_label,
                    convention,
                },
                MagneticRecordSelector::new(sg, condensing, Cdml)
                    .direction(IsotropyDirection::Index(row))
                    .uni(uni),
                Some(eta),
            )
            .with_embedding(embedding);
            let report = subduce(&request).unwrap();
            assert_eq!(report.source().labels.cdml, source);
            if sg == 191 {
                assert_eq!(report.source().labels.bc.as_deref(), Some("Γ6+"));
            }
            assert_eq!(report.source().input_convention, convention);
            assert_eq!(report.record().uni, uni);
            let rows = isotropy::magnetic_isotropy_subgroups(sg, condensing, Cdml).unwrap();
            assert_eq!(report.record().ordinal, rows[row].ordinal);
            assert_eq!(report.record().basis, rows[row].record.basis);
            assert_eq!(report.record().origin, rows[row].record.origin);
            let model = report.matrix_model().unwrap();
            let outcome = report.matrix().unwrap();
            let direct = subduce_matrix_gamma(model).unwrap();
            assert_eq!(report.source_dimension(), 2);
            assert_eq!(report.covered_dimension(), 2);
            assert_eq!(outcome.pieces.len(), count);
            assert_eq!(outcome.pieces.len(), direct.pieces().len());
            if uni == 1471 {
                assert_eq!(
                    report.classes().iter().filter(|g| g.time_reversal).count(),
                    12
                );
            }
            let mut columns = Vec::new();
            for (piece, old) in outcome.pieces.iter().zip(direct.pieces()) {
                assert_eq!(piece.target.unitary_sg, h);
                assert_eq!(piece.target.multiplicity, multiplicity);
                assert!(piece.target.magnetic_bc.is_none());
                assert_eq!(
                    piece.target.h_source,
                    irrep(h, old.source_label()).unwrap().labels()
                );
                assert_eq!(piece.target.magnetic_label, old.label());
                assert_eq!(piece.copies, old.copies());
                assert_eq!(piece.operators, old.operators());
                let d = piece.target.dimension as usize;
                assert_eq!(piece.copies.len(), multiplicity as usize);
                assert_eq!(piece.operators.len(), report.classes().len());
                for u in &piece.copies {
                    assert_eq!(u.len(), 2 * d);
                    for c in 0..d {
                        columns.push([u[c], u[d + c]]);
                    }
                    let z: Vec<_> = (0..d)
                        .map(|j| C::new(0.3 + j as f64, -0.6 + 0.2 * j as f64))
                        .collect();
                    let image = lift(u, d, &z);
                    for (s, t) in model.source_operators().iter().zip(&piece.operators) {
                        assert_eq!(s.is_antiunitary(), t.is_antiunitary());
                        let lhs = dense(&s.entries(), 2, s.is_antiunitary(), &image);
                        let rhs = lift(u, d, &dense(t.entries(), d, t.is_antiunitary(), &z));
                        for r in 0..2 {
                            close(lhs[r], rhs[r]);
                        }
                    }
                }
            }
            assert_eq!(columns.len(), 2);
            for r in 0..2 {
                for c in 0..2 {
                    close(
                        columns.iter().map(|v| v[r] * v[c].conj()).sum(),
                        C::new(if r == c { 1. } else { 0. }, 0.),
                    );
                }
            }
        }
    }
}

#[test]
fn selectors_refuse_missing_conflicting_and_unavailable_identities() {
    let source = SourceIdentity {
        parent_sg: 5,
        label: "GM1",
        convention: Cdml,
    };
    let selector = MagneticRecordSelector::new(5, "GM2", Cdml);
    let error = |record| {
        subduce(&MagneticSubductionRequest::new(
            source,
            record,
            Some(RealSign::Plus),
        ))
        .unwrap_err()
    };
    assert!(matches!(
        error(selector),
        ApiError::RecordSelectorMissing { .. }
    ));
    assert!(matches!(
        error(selector.uni(usize::MAX)),
        ApiError::RecordUniNotFound { .. }
    ));
    assert!(matches!(
        error(selector.direction(IsotropyDirection::Index(0)).uni(59)),
        ApiError::RecordUniMismatch { .. }
    ));
    assert!(matches!(
        error(selector.direction(IsotropyDirection::OfficialDescriptor("(a)"))),
        ApiError::Isotropy(IsotropyError::NoOfficialDescriptor { .. })
    ));
    assert!(matches!(
        error(selector.direction(IsotropyDirection::Descriptor("(a)"))),
        ApiError::Isotropy(IsotropyError::DirectionNotFound { .. })
    ));
    let bad_source = SourceIdentity {
        label: "Γ1",
        convention: Cdml,
        ..source
    };
    assert!(matches!(
        subduce(&MagneticSubductionRequest::new(
            bad_source,
            selector.uni(22),
            Some(RealSign::Plus)
        )),
        Err(ApiError::SourceLabelNotFound { .. })
    ));
    let ambiguous = SourceIdentity {
        parent_sg: 1,
        label: "Γ1",
        convention: Bc,
    };
    assert!(matches!(
        subduce(&MagneticSubductionRequest::new(
            ambiguous,
            MagneticRecordSelector::new(1, "GM1", Cdml).direction(IsotropyDirection::Index(0)),
            None
        )),
        Err(ApiError::AmbiguousSourceLabel { .. })
    ));
    let missing_eta = subduce(&MagneticSubductionRequest::new(
        source,
        selector.uni(22),
        None,
    ))
    .unwrap_err();
    assert!(matches!(missing_eta, ApiError::ScalarInput(_)));
    assert!(std::error::Error::source(&missing_eta).is_some());
}
