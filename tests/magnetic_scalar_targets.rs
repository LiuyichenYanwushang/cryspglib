//! Complex scalar targets keep their CIR identities and reconstruct the original
//! real source; character-only or metadata-only agreement is not sufficient.
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::{canonical_translation_checked, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::{
    HCharacterStorage, MatrixGammaModel, MatrixGammaRequest, MatrixSubductionError,
    ScalarTargetPiece, subduce_matrix_gamma, subduce_matrix_gamma_scalar_targets,
    verify_matrix_gamma_input,
};
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticRecordSelector, MagneticSubductionRequest, subduce_scalar_targets,
};
use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_uni_partial;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, query};
use num_complex::Complex64 as C;

fn model(
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    row: usize,
    eta: RealSign,
) -> MatrixGammaModel {
    verify_matrix_gamma_input(
        &MatrixGammaRequest::new(
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
            Some(eta),
        )
        .with_embedding(EmbeddingEvidence {
            map: Mat3R::identity(),
            shift: [Rat::ZERO; 3],
            lattice_policy: LatticePolicy::RequireExact,
        }),
    )
    .unwrap()
}

fn check(model: &MatrixGammaModel, pieces: &[ScalarTargetPiece]) {
    let classes = model.classes();
    let canonical = |t| canonical_translation_checked(t, model.translation_lattice()).unwrap();
    let partial = magnetic_irrep_summary_by_uni_partial(model.uni()).unwrap();
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|k| k.coords == (0, 0, 0, 1))
        .unwrap();
    // Summary display coordinates are H data-Hall, not the model's first Hall.
    // Bind only through the independently enumerated original operation index.
    let raw = cryspglib::irrep::magnetic_embedding::magnetic_operations(model.uni()).unwrap();
    assert_eq!(raw.operations.len(), gamma.operations.len());
    for (position, op) in gamma.operations.iter().enumerate() {
        assert_eq!(op.column, position);
        assert_eq!(op.magnetic_operation_index, position);
        assert_eq!(op.time_reversal, raw.operations[position].time_reversal);
    }
    let mut columns = Vec::new();
    for piece in pieces {
        assert_eq!(piece.dimension(), 1);
        assert_eq!(piece.operators().len(), classes.len());
        let summary = gamma
            .coreps
            .iter()
            .find(|c| c.label == piece.label())
            .unwrap();
        for (i, g) in classes.iter().enumerate() {
            let op = piece.operators()[i];
            assert_eq!(op.is_antiunitary(), g.time_reversal);
            assert!((op.value().norm_sqr() - 1.0).abs() < 1e-9);
            if !g.time_reversal {
                let matches: Vec<_> = gamma
                    .operations
                    .iter()
                    .filter(|h| {
                        let original = &raw.operations[h.magnetic_operation_index];
                        !original.time_reversal
                            && original.rotation == g.rotation
                            && canonical(original.translation) == canonical(g.translation)
                    })
                    .collect();
                assert!(!matches.is_empty());
                for h in matches {
                    assert!((summary.characters[h.column].unwrap() - op.value()).norm() < 1e-9);
                }
            }
            for (j, h) in classes.iter().enumerate() {
                let product = compose(g, h).unwrap();
                let k = classes
                    .iter()
                    .position(|k| {
                        k.rotation == product.rotation
                            && k.time_reversal == product.time_reversal
                            && canonical(k.translation) == canonical(product.translation)
                    })
                    .unwrap();
                for z in [C::new(1.0, 0.0), C::new(0.0, 1.0), C::new(0.7, -0.2)] {
                    assert!(
                        (op.apply(piece.operators()[j].apply(z)) - piece.operators()[k].apply(z))
                            .norm()
                            < 1e-9
                    );
                }
            }
        }
        for u in piece.copies() {
            assert_eq!(u.len(), 2);
            for (g, source) in model.source_operators().iter().enumerate() {
                for z in [C::new(1.0, 0.0), C::new(0.0, 1.0), C::new(-0.3, 0.8)] {
                    let left = source.apply([u[0] * z, u[1] * z]);
                    let right = piece.operators()[g].apply(z);
                    for r in 0..2 {
                        assert!((left[r] - u[r] * right).norm() < 1e-9);
                    }
                }
            }
            columns.push(u);
        }
    }
    assert_eq!(columns.len(), 2);
    for i in 0..2 {
        for j in 0..2 {
            let gram: C = (0..2).map(|r| columns[i][r].conj() * columns[j][r]).sum();
            let projector: C = columns.iter().map(|u| u[i] * u[j].conj()).sum();
            let expected = C::new(f64::from(i == j), 0.0);
            assert!((gram - expected).norm() < 1e-9);
            assert!((projector - expected).norm() < 1e-9);
        }
    }
}

#[test]
fn complex_targets_are_authoritative_named_constituents() {
    for eta in [RealSign::Plus, RealSign::Minus] {
        let m = model(89, "GM5", "GM3", 0, eta);
        assert_eq!(m.uni(), 750);
        assert!(matches!(
            subduce_matrix_gamma(&m),
            Err(MatrixSubductionError::TargetSourceUnresolved { .. })
        ));
        let result = subduce_matrix_gamma_scalar_targets(&m).unwrap();
        assert_eq!(
            (result.source_dimension(), result.covered_dimension()),
            (2, 2)
        );
        assert_eq!(
            result.scope(),
            cryspglib::irrep::magnetic_matrix_subduction::SCALAR_TARGET_SCOPE
        );
        check(&m, result.pieces());
        assert_eq!(result.pieces().len(), 2);
        for p in result.pieces() {
            assert!(matches!(p.source().ml, "GM3" | "GM4"));
            assert!(query::find_irreps(75, p.source().ml, LabelConvention::Cdml).is_empty());
            assert!(p.source().labels().bc.is_none());
            let HCharacterStorage::Constituent {
                record,
                component,
                irnumber,
            } = p.storage()
            else {
                panic!("not a standalone row")
            };
            let owner = query::find_irreps(75, record.cdml, LabelConvention::Cdml);
            assert_eq!(owner.len(), 1);
            let metadata = owner[0].compound_metadata().unwrap();
            assert_eq!(metadata.cir_labels[*component], p.source().ml);
            assert_eq!(metadata.cir_irnumbers[*component], *irnumber);
            assert!(p.operators().iter().any(|o| o.value().im.abs() > 0.5));
        }
    }
}

#[test]
fn high_level_complex_targets_preserve_screws_centering_origins_and_labels() {
    for (sg, source, condensing, uni) in [
        (89, "GM5", "GM3", 750),
        (91, "GM5", "GM3", 766),
        (97, "GM5", "GM3", 814),
        (141, "GM5+", "GM3+", 1217),
    ] {
        let s = query::find_irreps(sg, source, LabelConvention::Cdml)[0];
        let c = query::find_irreps(sg, condensing, LabelConvention::Cdml)[0];
        for eta in [RealSign::Plus, RealSign::Minus] {
            for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                let sl = s.label(convention).unwrap();
                let cl = c.label(convention).unwrap();
                let request = MagneticSubductionRequest::new(
                    SourceIdentity {
                        parent_sg: sg,
                        label: &sl,
                        convention,
                    },
                    MagneticRecordSelector::new(sg, &cl, convention)
                        .direction(IsotropyDirection::Index(0))
                        .uni(uni),
                    Some(eta),
                )
                .with_embedding(EmbeddingEvidence {
                    map: Mat3R::identity(),
                    shift: [Rat::ZERO; 3],
                    lattice_policy: LatticePolicy::RequireExact,
                });
                let report = subduce_scalar_targets(&request).unwrap();
                assert_eq!(
                    report.scope(),
                    cryspglib::irrep::magnetic_matrix_subduction::SCALAR_TARGET_SCOPE
                );
                let m = report.matrix_model();
                assert_eq!(
                    (m.uni(), m.source_label(), m.record_label()),
                    (uni, source, condensing)
                );
                assert_eq!(m.time_reversal(), Some(eta));
                assert_eq!(report.covered_dimension(), 2);
                assert_eq!(report.source().labels.cdml, source);
                assert_eq!(report.condensing().labels.cdml, condensing);
                assert_eq!(*m.map(), Mat3R::identity());
                assert_eq!(*m.shift(), [Rat::ZERO; 3]);
                if sg == 141 {
                    assert_ne!(*m.source_origin_shift(), [Rat::ZERO; 3]);
                }
                if sg == 97 || sg == 141 {
                    assert!(!m.translation_lattice().is_empty());
                }
                if sg == 91 {
                    assert!(m.classes().iter().any(|g| g.translation != [Rat::ZERO; 3]));
                }
                check(m, report.pieces());
            }
        }
    }
}

#[test]
fn real_scalar_splits_and_multiplicities_remain_complete() {
    for (sg, source, condensing, row) in [
        (89, "GM5", "GM5", 0),
        (89, "GM5", "GM5", 2),
        (156, "GM3", "GM3", 2),
        (207, "GM3", "GM4", 0),
    ] {
        for eta in [RealSign::Plus, RealSign::Minus] {
            let m = model(sg, source, condensing, row, eta);
            let old = subduce_matrix_gamma(&m).unwrap();
            let new = subduce_matrix_gamma_scalar_targets(&m).unwrap();
            let a: Vec<_> = old
                .pieces()
                .iter()
                .map(|p| (p.label(), p.multiplicity()))
                .collect();
            let b: Vec<_> = new
                .pieces()
                .iter()
                .map(|p| (p.label(), p.multiplicity()))
                .collect();
            assert_eq!(a, b);
            assert!(
                new.pieces()
                    .iter()
                    .all(|p| matches!(p.storage(), HCharacterStorage::Standalone { .. }))
            );
            check(&m, new.pieces());
        }
    }
}

#[test]
fn non_scalar_targets_and_non_matrix_sources_are_not_silently_coerced() {
    for condensing in ["GM1", "GM2"] {
        let m = model(89, "GM5", condensing, 0, RealSign::Plus);
        assert!(matches!(
            subduce_matrix_gamma_scalar_targets(&m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
    }
    let req = MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(89, "GM1", LabelConvention::Cdml)
            .direction(IsotropyDirection::Index(0)),
        Some(RealSign::Plus),
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::identity(),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    });
    assert!(subduce_scalar_targets(&req).is_err());
}
