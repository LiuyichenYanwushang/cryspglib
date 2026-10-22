//! Leader-owned R32 boundary checks; independent dense algebra, no engine apply helpers.
use cryspglib::irrep::LabelConvention::{self, Bc, Cdml};
use cryspglib::irrep::magnetic_embedding::{
    MagneticOperation, canonical_translation, compose, magnetic_operations,
};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::MatrixGammaModel;
use cryspglib::irrep::magnetic_subduction_api as api;
use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_uni_partial;
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use num_complex::Complex64 as C;

#[test]
fn legacy_public_struct_literals_and_exhaustive_matches_remain_valid() {
    // This intentionally uses the old exact public struct shapes.  Adding a
    // required partner field or outcome variant must not silently break callers.
    let target = api::MagneticTargetIdentity {
        magnetic_label: "GM1".to_string(),
        magnetic_bc: None,
        unitary_sg: 1,
        h_source: query::find_irreps(1, "GM1", Cdml)[0].labels(),
        dimension: 1,
        multiplicity: 1,
    };
    let scalar = api::ScalarGammaOutcome {
        target: target.clone(),
        intertwiner: C::new(1.0, 0.0),
        operators: Vec::new(),
    };
    let piece = api::MatrixGammaPieceOutcome {
        target: api::MagneticTargetIdentity {
            multiplicity: 2,
            ..target
        },
        copies: vec![
            vec![C::new(1.0, 0.0), C::new(0.0, 0.0)],
            vec![C::new(0.0, 0.0), C::new(1.0, 0.0)],
        ],
        operators: Vec::new(),
    };
    let matrix = api::MatrixGammaOutcome {
        pieces: vec![piece],
        source_dimension: 2,
        covered_dimension: 2,
    };
    // The fabricated payloads above test layout compatibility only, not a
    // representation law.  Actual engine-produced payloads are checked below.
    for (outcome, dimension) in [
        (api::MagneticOutcome::Scalar(scalar), 1),
        (api::MagneticOutcome::Matrix(matrix), 2),
    ] {
        let covered = match outcome {
            api::MagneticOutcome::Scalar(result) => result.target.covered_dimension(),
            api::MagneticOutcome::Matrix(result) => result.covered_dimension,
        };
        assert_eq!(covered, dimension);
    }
    let report = api::subduce(&api::MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 1,
            label: "GM1",
            convention: Cdml,
        },
        api::MagneticRecordSelector::new(1, "GM1", Cdml).uni(1),
        Some(RealSign::Plus),
    ))
    .unwrap();
    match report.model() {
        api::MagneticVerifiedModel::Scalar(model) => assert_eq!(model.source_label(), "GM1"),
        api::MagneticVerifiedModel::Matrix(_) => panic!("scalar source changed legacy variant"),
    }
}

type M = [C; 4];
fn identity() -> M {
    [
        C::new(1.0, 0.0),
        C::new(0.0, 0.0),
        C::new(0.0, 0.0),
        C::new(1.0, 0.0),
    ]
}
fn multiply(a: M, b: M) -> M {
    std::array::from_fn(|i| (0..2).map(|k| a[(i / 2) * 2 + k] * b[k * 2 + i % 2]).sum())
}
fn conjugate(a: M) -> M {
    a.map(|v| v.conj())
}
fn dagger(a: M) -> M {
    [a[0].conj(), a[2].conj(), a[1].conj(), a[3].conj()]
}
fn same(a: M, b: M) -> bool {
    a.iter()
        .zip(b)
        .all(|(x, y)| x.re.is_finite() && x.im.is_finite() && (x - y).norm() < 1e-8)
}
fn real(a: &[f64]) -> M {
    assert_eq!(a.len(), 4);
    std::array::from_fn(|i| C::new(a[i], 0.0))
}
fn complex(a: &[C]) -> M {
    assert_eq!(a.len(), 4);
    [a[0], a[1], a[2], a[3]]
}
fn mv(a: M, z: [C; 2]) -> [C; 2] {
    [a[0] * z[0] + a[1] * z[1], a[2] * z[0] + a[3] * z[1]]
}
fn close_vec(a: [C; 2], b: [C; 2]) {
    for i in 0..2 {
        assert!(
            a[i].re.is_finite() && a[i].im.is_finite() && (a[i] - b[i]).norm() < 1e-8,
            "{a:?} != {b:?}"
        );
    }
}
fn class_of(g: &MagneticOperation, model: &MatrixGammaModel) -> usize {
    let found: Vec<_> = model
        .classes()
        .iter()
        .enumerate()
        .filter(|(_, h)| {
            g.rotation == h.rotation
                && g.time_reversal == h.time_reversal
                && canonical_translation(g.translation, model.translation_lattice())
                    == canonical_translation(h.translation, model.translation_lattice())
        })
        .map(|(i, _)| i)
        .collect();
    assert_eq!(found.len(), 1);
    found[0]
}
fn products(model: &MatrixGammaModel) -> Vec<Vec<usize>> {
    model
        .classes()
        .iter()
        .map(|g| {
            model
                .classes()
                .iter()
                .map(|h| class_of(&compose(g, h).unwrap(), model))
                .collect()
        })
        .collect()
}
fn obeys_law(table: &[M], anti: &[bool], product: &[Vec<usize>]) -> bool {
    (0..table.len()).all(|i| {
        (0..table.len()).all(|j| {
            same(
                multiply(
                    table[i],
                    if anti[i] {
                        conjugate(table[j])
                    } else {
                        table[j]
                    },
                ),
                table[product[i][j]],
            )
        })
    })
}
fn reconstructs(source: &[M], target: &[M], anti: &[bool], u: M) -> bool {
    (0..source.len()).all(|i| {
        same(
            multiply(source[i], if anti[i] { conjugate(u) } else { u }),
            multiply(u, target[i]),
        )
    })
}
fn request<'a>(
    source: &'a str,
    sc: LabelConvention,
    condensing: &'a str,
    cc: LabelConvention,
    uni: usize,
    eta: RealSign,
) -> api::MagneticSubductionRequest<'a> {
    api::MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: source,
            convention: sc,
        },
        api::MagneticRecordSelector::new(89, condensing, cc).uni(uni),
        Some(eta),
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::from_ints([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    })
}

#[test]
fn paired_public_payload_has_independent_laws_full_reconstruction_and_nontrivial_complex_gauge() {
    for (condensing, uni, hsg) in [("GM2", 749, 16), ("GM4", 751, 21)] {
        for eta in [RealSign::Plus, RealSign::Minus] {
            let req = request("GM5", Cdml, condensing, Cdml, uni, eta);
            let report = api::subduce_type_c(&req).unwrap();
            let model = report.model();
            let payload = report.outcome();
            assert_eq!(model.uni(), uni);
            assert_eq!(
                (
                    payload.target.unitary_sg,
                    payload.target.dimension,
                    payload.target.multiplicity
                ),
                (hsg, 2, 1)
            );
            assert_eq!(payload.copies.len(), 1);
            assert_eq!(payload.target.magnetic_bc, None);
            assert_eq!(payload.operators.len(), model.classes().len());
            assert_eq!(model.classes().len(), 8);
            let u = complex(&payload.copies[0]);
            assert!(same(multiply(dagger(u), u), identity()));
            assert!(same(multiply(u, dagger(u)), identity()));
            assert!((u[0] * u[3] - u[1] * u[2]).norm() > 0.9);
            let source: Vec<M> = model
                .source_operators()
                .iter()
                .map(|g| real(&g.entries()))
                .collect();
            let target: Vec<M> = payload
                .operators
                .iter()
                .map(|g| real(g.entries()))
                .collect();
            let anti: Vec<_> = model.classes().iter().map(|g| g.time_reversal).collect();
            for (g, t) in model.classes().iter().zip(&payload.operators) {
                assert_eq!(g.time_reversal, t.is_antiunitary());
            }
            let product = products(model);
            assert!(obeys_law(&source, &anti, &product));
            assert!(obeys_law(&target, &anti, &product));
            assert!(reconstructs(&source, &target, &anti, u));
            // Type C has two real Hom directions, not two independent copies:
            // F=diag(i,-i) is a non-scalar semilinear self-intertwiner with F²=-I.
            let f = [
                C::new(0.0, 1.0),
                C::new(0.0, 0.0),
                C::new(0.0, 0.0),
                C::new(0.0, -1.0),
            ];
            let uf = multiply(u, f);
            assert!(reconstructs(&source, &target, &anti, uf));
            assert!(same(multiply(dagger(u), uf), f));
            assert!(
                !same(multiply(dagger(u), uf), [C::new(0.0, 0.0); 4]),
                "the two solutions do not have orthogonal images"
            );
            for (i, a) in source.iter().enumerate() {
                // Full source operator, not just an invariant subspace:
                // A=U D U† (unitary), A=U D Uᵀ (antiunitary).
                assert!(same(
                    multiply(
                        multiply(u, target[i]),
                        if anti[i] {
                            conjugate(dagger(u))
                        } else {
                            dagger(u)
                        }
                    ),
                    *a
                ));
                for z in [
                    [C::new(1.0, 1.0), C::new(2.0, -3.0)],
                    [C::new(0.0, 1.0), C::new(-1.0, 0.5)],
                    [C::new(2.0, 0.0), C::new(0.0, 1.0)],
                ] {
                    let lifted = mv(u, z);
                    let lhs = mv(
                        *a,
                        if anti[i] {
                            lifted.map(|v| v.conj())
                        } else {
                            lifted
                        },
                    );
                    let rhs = mv(
                        u,
                        mv(target[i], if anti[i] { z.map(|v| v.conj()) } else { z }),
                    );
                    close_vec(lhs, rhs);
                }
            }
            // Independent complex basis change (the public target coefficients
            // are real; this exercises their general complex-equivalent action).
            let gauge = [
                C::from_polar(1.0, 0.37),
                C::new(0.0, 0.0),
                C::new(0.0, 0.0),
                C::from_polar(1.0, -0.53),
            ];
            let ug = multiply(u, gauge);
            let gauged: Vec<M> = target
                .iter()
                .enumerate()
                .map(|(i, d)| {
                    multiply(
                        multiply(dagger(gauge), *d),
                        if anti[i] { conjugate(gauge) } else { gauge },
                    )
                })
                .collect();
            assert!(obeys_law(&gauged, &anti, &product));
            assert!(reconstructs(&source, &gauged, &anti, ug));
            assert!(
                !reconstructs(&source, &target, &anti, ug),
                "wrong U without the coupled target gauge must fail"
            );
            let a = anti.iter().position(|a| *a).unwrap();
            let mut tampered = target.clone();
            tampered[a] = tampered[a].map(|z| -z);
            assert!(
                same(
                    multiply(tampered[a], conjugate(tampered[a])),
                    target[product[a][a]]
                ),
                "single-anti sign change preserves its square"
            );
            assert!(
                !obeys_law(&tampered, &anti, &product),
                "full cross-class law must still catch single-anti corruption"
            );
            assert!(!reconstructs(&source, &tampered, &anti, u));
            // Old single-H entry must never present a paired corep as one H.
            assert!(api::subduce(&req).is_err());
        }
    }
}

#[test]
fn paired_target_diagonal_is_anchored_to_both_typed_h_rows_via_raw_m_indices() {
    for (condensing, uni, hsg) in [("GM2", 749, 16), ("GM4", 751, 21)] {
        let report =
            api::subduce_type_c(&request("GM5", Cdml, condensing, Cdml, uni, RealSign::Plus))
                .unwrap();
        let payload = report.outcome();
        let model = report.model();
        let mut typed = Vec::new();
        for h in &payload.target.h_sources {
            let found = query::find_irreps(hsg, h.cdml, Cdml);
            assert_eq!(found.len(), 1);
            assert_eq!(*h, found[0].labels());
            typed.push(found[0].ordinary_scalar_selected_arm_block_trace().unwrap());
        }
        assert_eq!(payload.target.h_sources[0].cdml, "GM3");
        assert_eq!(payload.target.h_sources[0].bc.as_deref(), Some("Γ4"));
        assert_eq!(payload.target.h_sources[1].cdml, "GM4");
        assert_eq!(payload.target.h_sources[1].bc.as_deref(), Some("Γ2"));
        let partial = magnetic_irrep_summary_by_uni_partial(uni).unwrap();
        let gamma = partial
            .summary
            .kpoints
            .iter()
            .find(|k| k.coords == (0, 0, 0, 1))
            .unwrap();
        let raw = magnetic_operations(uni).unwrap();
        let mut seen = 0;
        for col in gamma.operations.iter().filter(|g| !g.time_reversal) {
            let i = class_of(&raw.operations[col.magnetic_operation_index], model);
            let rotation = [
                col.rotation[0][0],
                col.rotation[0][1],
                col.rotation[0][2],
                col.rotation[1][0],
                col.rotation[1][1],
                col.rotation[1][2],
                col.rotation[2][0],
                col.rotation[2][1],
                col.rotation[2][2],
            ];
            let mut characters = [C::new(0.0, 0.0); 2];
            for h in 0..2 {
                let js: Vec<_> = typed[h]
                    .operations()
                    .iter()
                    .enumerate()
                    .filter(|(_, op)| {
                        op.rotation == rotation
                            && (0..3).all(|k| (op.translation[k] - col.translation[k]).abs() < 1e-8)
                    })
                    .map(|(j, _)| j)
                    .collect();
                assert_eq!(js.len(), 1);
                characters[h] = typed[h].values()[js[0]];
                assert!(characters[h].re.is_finite() && characters[h].im.is_finite());
            }
            assert!(same(
                real(payload.operators[i].entries()),
                [
                    characters[0],
                    C::new(0.0, 0.0),
                    C::new(0.0, 0.0),
                    characters[1]
                ]
            ));
            seen += 1;
        }
        assert!(seen >= 4);
    }
}

#[test]
fn paired_owned_result_keeps_independent_input_conventions_and_two_genuine_h_identities() {
    for (condensing, uni) in [("GM2", 749), ("GM4", 751)] {
        for sc in [Cdml, Bc] {
            for cc in [Cdml, Bc] {
                for eta in [RealSign::Plus, RealSign::Minus] {
                    let report = {
                        let sl = query::find_irreps(89, "GM5", Cdml)[0]
                            .label(sc)
                            .unwrap()
                            .to_string();
                        let cl = query::find_irreps(89, condensing, Cdml)[0]
                            .label(cc)
                            .unwrap()
                            .to_string();
                        let req = request(&sl, sc, &cl, cc, uni, eta);
                        api::subduce_type_c(&req).unwrap()
                    };
                    assert_eq!(
                        report.source().labels,
                        query::find_irreps(89, "GM5", Cdml)[0].labels()
                    );
                    assert_eq!(
                        report.condensing().labels,
                        query::find_irreps(89, condensing, Cdml)[0].labels()
                    );
                    assert_eq!(report.source().input_convention, sc);
                    assert_eq!(report.condensing().input_convention, cc);
                    assert_eq!(report.record().uni, uni);
                    assert_eq!(report.model().time_reversal(), Some(eta));
                    assert_eq!(
                        *report.model().map(),
                        Mat3R::from_ints([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                    );
                    assert_eq!(*report.model().shift(), [Rat::ZERO; 3]);
                    assert_eq!(report.model().lattice_policy(), LatticePolicy::RequireExact);
                    assert!(report.model().realises_record_lattice());
                    assert_eq!(report.outcome().target.h_sources.len(), 2);
                    assert!(report.outcome().target.magnetic_bc.is_none());
                    assert_eq!(
                        report.outcome().target.dimension * report.outcome().target.multiplicity,
                        2
                    );
                }
            }
        }
    }
}
