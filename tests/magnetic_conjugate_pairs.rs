//! Both the magnetic C carrier and the unitary packet retain named H rows.
use cryspglib::irrep::isotropy::IsotropyDirection;
use cryspglib::irrep::magnetic_embedding::{
    canonical_translation_checked, compose, magnetic_operations,
};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::{
    ConjugatePairKind, ConjugatePairPiece, HCharacterStorage, MatrixGammaModel, MatrixGammaRequest,
    MatrixSubductionError, subduce_matrix_gamma, subduce_matrix_gamma_conjugate_pair,
    subduce_matrix_gamma_scalar_targets, subduce_matrix_gamma_type_c, verify_matrix_gamma_input,
};
use cryspglib::irrep::magnetic_subduction_api::{
    MagneticRecordSelector, MagneticSubductionRequest, subduce_conjugate_pair,
};
use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_uni_partial;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, query, types::CompoundSelectedArmCharacter};
use num_complex::Complex64 as C;

#[derive(Clone, Copy)]
struct Case {
    sg: u8,
    source: &'static str,
    cond: &'static str,
    row: usize,
    uni: usize,
    map: [[i128; 3]; 3],
    kind: ConjugatePairKind,
}
const CASES: [Case; 4] = [
    Case {
        sg: 100,
        source: "GM5",
        cond: "A5",
        row: 0,
        uni: 688,
        map: [[-1, -1, 0], [-1, 1, 0], [0, 0, -2]],
        kind: ConjugatePairKind::MagneticTypeC,
    },
    Case {
        sg: 149,
        source: "GM3",
        cond: "H1",
        row: 2,
        uni: 1233,
        map: [[-2, 1, 0], [-1, 2, 0], [0, 0, -2]],
        kind: ConjugatePairKind::MagneticTypeC,
    },
    Case {
        sg: 149,
        source: "GM3",
        cond: "K1",
        row: 2,
        uni: 1231,
        map: [[-2, 1, 0], [-1, 2, 0], [0, 0, -1]],
        kind: ConjugatePairKind::UnitaryPair,
    },
    Case {
        sg: 141,
        source: "GM5+",
        cond: "P1",
        row: 4,
        uni: 702,
        map: [[1, -1, 0], [-1, -1, 0], [0, 0, -2]],
        kind: ConjugatePairKind::MagneticTypeC,
    },
];
fn evidence(c: Case) -> EmbeddingEvidence {
    EmbeddingEvidence {
        map: Mat3R::new(c.map.map(|r| r.map(|x| Rat::new(x, 1).unwrap()))),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    }
}
fn model(c: Case, eta: RealSign) -> MatrixGammaModel {
    verify_matrix_gamma_input(
        &MatrixGammaRequest::new(
            SourceIdentity {
                parent_sg: c.sg,
                label: c.source,
                convention: LabelConvention::Cdml,
            },
            RecordIdentity {
                parent_sg: c.sg,
                condensing_label: c.cond,
                convention: LabelConvention::Cdml,
                direction_index: c.row,
            },
            Some(eta),
        )
        .with_embedding(evidence(c)),
    )
    .unwrap()
}
fn check(m: &MatrixGammaModel, p: &ConjugatePairPiece) {
    assert_eq!((p.dimension(), p.multiplicity()), (2, 1));
    let u = &p.copies()[0];
    assert_eq!(u.len(), 4);
    assert_eq!(p.operators().len(), m.classes().len());
    for i in 0..2 {
        for j in 0..2 {
            let gram: C = (0..2).map(|r| u[r * 2 + i].conj() * u[r * 2 + j]).sum();
            let projection: C = (0..2).map(|r| u[i * 2 + r] * u[j * 2 + r].conj()).sum();
            let want = C::new(f64::from(i == j), 0.);
            assert!((gram - want).norm() < 1e-8);
            assert!((projection - want).norm() < 1e-8)
        }
    }
    let canon = |t| canonical_translation_checked(t, m.translation_lattice()).unwrap();
    let index = |op: &cryspglib::irrep::magnetic_embedding::MagneticOperation| {
        m.classes()
            .iter()
            .position(|x| {
                x.rotation == op.rotation
                    && x.time_reversal == op.time_reversal
                    && canon(x.translation) == canon(op.translation)
            })
            .unwrap()
    };
    for (i, s) in m.source_operators().iter().enumerate() {
        let t = &p.operators()[i];
        assert_eq!(t.is_antiunitary(), s.is_antiunitary());
        assert_eq!(t.entries().len(), 4);
        assert!((t.entries()[0] - t.entries()[3]).abs() < 1e-8);
        assert!((t.entries()[1] + t.entries()[2]).abs() < 1e-8);
        for z in [
            [C::new(1., 0.), C::new(0., 0.)],
            [C::new(0., 1.), C::new(0., 0.)],
            [C::new(0., 0.), C::new(1., 0.)],
            [C::new(0., 0.), C::new(0., 1.)],
            [C::new(0.3, -0.4), C::new(-0.7, 0.6)],
        ] {
            let lhs = s.apply([u[0] * z[0] + u[1] * z[1], u[2] * z[0] + u[3] * z[1]]);
            let v = t.apply(&z);
            let rhs = [u[0] * v[0] + u[1] * v[1], u[2] * v[0] + u[3] * v[1]];
            for r in 0..2 {
                assert!((lhs[r] - rhs[r]).norm() < 1e-8)
            }
        }
        for j in 0..m.classes().len() {
            let k = index(&compose(&m.classes()[i], &m.classes()[j]).unwrap());
            for z in [
                [C::new(1., 0.), C::new(0., 0.)],
                [C::new(0., 0.), C::new(0., 1.)],
            ] {
                let lhs = t.apply(&p.operators()[j].apply(&z));
                let rhs = p.operators()[k].apply(&z);
                for r in 0..2 {
                    assert!((lhs[r] - rhs[r]).norm() < 1e-8)
                }
            }
        }
    }
    let partial = magnetic_irrep_summary_by_uni_partial(m.uni()).unwrap();
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|g| g.coords == (0, 0, 0, 1))
        .unwrap();
    let target = gamma.coreps.iter().find(|x| x.label == p.label()).unwrap();
    assert_eq!(target.source_irreps.len(), 2);
    for (a, b) in p.h_sources().iter().zip(&target.source_irreps) {
        assert_eq!(
            (a.sg, a.ml, a.bc, a.dim, a.spinor),
            (b.sg, b.ml, b.bc, b.dim, b.spinor)
        );
    }
    let raw = magnetic_operations(m.uni()).unwrap();
    let mut rows = Vec::new();
    for (slot, storage) in p.storage().iter().enumerate() {
        let HCharacterStorage::Constituent {
            record,
            component,
            irnumber,
        } = storage
        else {
            panic!("fixture needs an authoritative CIR")
        };
        let owner =
            query::find_irreps(p.h_sources()[slot].sg, record.cdml, LabelConvention::Cdml)[0];
        let metadata = owner.compound_metadata().unwrap();
        assert_eq!(metadata.cir_labels[*component], p.h_sources()[slot].ml);
        assert_eq!(metadata.cir_irnumbers[*component], *irnumber);
        assert_eq!(p.h_sources()[slot].labels().bc, None);
        let CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } =
            owner.compound_selected_arm_view().unwrap()
        else {
            panic!()
        };
        rows.push(if *component == first.component {
            first.row
        } else {
            assert_eq!(*component, second.component);
            second.row
        });
    }
    let mut h_values = vec![None; m.classes().len()];
    for col in gamma.operations.iter().filter(|c| !c.time_reversal) {
        let i = index(&raw.operations[col.magnetic_operation_index]);
        let mut values = Vec::new();
        for row in &rows {
            let j = row
                .operations()
                .iter()
                .position(|h| {
                    h.rotation.as_slice() == col.rotation.as_flattened()
                        && h.translation
                            .iter()
                            .zip(col.translation)
                            .all(|(a, b)| ((a - b) - (a - b).round()).abs() < 1e-8)
                })
                .unwrap();
            values.push(row.values()[j]);
        }
        h_values[i] = Some(values[0]);
        assert!((values[0].conj() - values[1]).norm() < 1e-8);
        let a = p.operators()[i].entries();
        assert!((C::new(a[0], a[2]) - values[0]).norm() < 1e-8);
        assert!((values[0] + values[1] - target.characters[col.column].unwrap()).norm() < 1e-8);
    }
    // Use the authoritative H row and original M products independently of
    // the target matrices to test the canonical antiunitary coefficient.
    if let Some(a) = m.classes().iter().find(|g| g.time_reversal) {
        let square = index(&compose(a, a).unwrap());
        let lambda = h_values[square].unwrap().sqrt();
        let mut wrong_conjugation_differs = false;
        for (h, operation) in m
            .classes()
            .iter()
            .enumerate()
            .filter(|(_, g)| !g.time_reversal)
        {
            let chi = h_values[h].unwrap();
            let g = index(&compose(a, operation).unwrap());
            let want = lambda * chi;
            let wrong = lambda * chi.conj();
            let actual = p.operators()[g].entries();
            for (x, y) in actual.iter().zip([want.re, -want.im, want.im, want.re]) {
                assert!(
                    (x - y).abs() < 1e-8,
                    "UNI{}: T(a0 h) must use lambda*chi(h)",
                    m.uni()
                );
            }
            wrong_conjugation_differs |= (want - wrong).norm() > 1e-8;
        }
        assert!(
            wrong_conjugation_differs,
            "the witness must distinguish the wrong conjugated formula"
        );
    }
}
#[test]
fn authoritative_complex_pairs_reconstruct_both_magnetic_and_unitary_carriers() {
    for c in CASES {
        for eta in [RealSign::Plus, RealSign::Minus] {
            let m = model(c, eta);
            assert_eq!(m.uni(), c.uni);
            assert!(subduce_matrix_gamma(&m).is_err());
            assert!(subduce_matrix_gamma_type_c(&m).is_err());
            assert!(subduce_matrix_gamma_scalar_targets(&m).is_err());
            let out = subduce_matrix_gamma_conjugate_pair(&m).unwrap();
            assert_eq!(out.piece().kind(), c.kind);
            assert_eq!((out.source_dimension(), out.covered_dimension()), (2, 2));
            assert_eq!(
                out.scope(),
                cryspglib::irrep::magnetic_matrix_subduction::CONJUGATE_PAIR_SCOPE
            );
            assert_eq!(
                c.kind == ConjugatePairKind::UnitaryPair,
                m.classes().iter().all(|g| !g.time_reversal)
            );
            check(&m, out.piece());
        }
    }
}
#[test]
fn owned_pair_reports_preserve_both_label_conventions_and_original_geometry() {
    for c in CASES {
        let source = query::find_irreps(c.sg, c.source, LabelConvention::Cdml)[0];
        let cond = query::find_irreps(c.sg, c.cond, LabelConvention::Cdml)[0];
        for eta in [RealSign::Plus, RealSign::Minus] {
            for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                let sl = source.label(convention).unwrap();
                let cl = cond.label(convention).unwrap();
                let req = MagneticSubductionRequest::new(
                    SourceIdentity {
                        parent_sg: c.sg,
                        label: &sl,
                        convention,
                    },
                    MagneticRecordSelector::new(c.sg, &cl, convention)
                        .direction(IsotropyDirection::Index(c.row))
                        .uni(c.uni),
                    Some(eta),
                )
                .with_embedding(evidence(c));
                let out = subduce_conjugate_pair(&req).unwrap();
                let m = out.matrix_model();
                assert_eq!(out.piece().kind(), c.kind);
                assert_eq!(out.source().labels.cdml, c.source);
                assert_eq!(out.condensing().labels.cdml, c.cond);
                assert_eq!((out.source_dimension(), out.covered_dimension()), (2, 2));
                assert_eq!(*m.map(), evidence(c).map);
                assert_eq!(*m.shift(), [Rat::ZERO; 3]);
                assert_eq!(m.time_reversal(), Some(eta));
                assert!(out.embedding().realises_record_lattice);
                if c.sg == 141 {
                    assert_ne!(*m.source_origin_shift(), [Rat::ZERO; 3]);
                    assert!(!m.translation_lattice().is_empty());
                }
                check(m, out.piece());
            }
        }
    }
}
#[test]
fn real_pairs_scalar_a_and_irreducible_real_targets_stay_outside_this_entry() {
    for (cond, row) in [("GM1", 0), ("GM2", 0), ("GM3", 0), ("GM5", 0)] {
        let c = Case {
            sg: 89,
            source: "GM5",
            cond,
            row,
            uni: 0,
            map: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            kind: ConjugatePairKind::MagneticTypeC,
        };
        let m = model(c, RealSign::Plus);
        if cond == "GM2" {
            subduce_matrix_gamma_type_c(&m).unwrap();
        }
        assert!(matches!(
            subduce_matrix_gamma_conjugate_pair(&m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
    }
}

#[test]
fn the_new_entry_does_not_widen_the_source_domain() {
    let c = CASES[2];
    let req = MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: c.sg,
            label: "GM1",
            convention: LabelConvention::Cdml,
        },
        MagneticRecordSelector::new(c.sg, c.cond, LabelConvention::Cdml)
            .direction(IsotropyDirection::Index(c.row))
            .uni(c.uni),
        Some(RealSign::Plus),
    )
    .with_embedding(evidence(c));
    assert!(
        matches!(subduce_conjugate_pair(&req),Err(cryspglib::irrep::magnetic_subduction_api::MagneticSubductionApiError::Matrix(error)) if matches!(*error,MatrixSubductionError::SourceDimensionNotTwo{..}))
    );
}
