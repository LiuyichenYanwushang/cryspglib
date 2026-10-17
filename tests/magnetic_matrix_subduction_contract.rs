//! Leader-owned independent checks of the R29 public contract.
use cryspglib::irrep::magnetic_embedding::{canonical_translation, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{LabelConvention, magnetic_matrix_subduction as matrix};
use num_complex::Complex64 as C;

fn request(
    sg: u8,
    source: &'static str,
    condensing: &'static str,
    index: usize,
    map: [[i32; 3]; 3],
    eta: RealSign,
) -> matrix::MatrixGammaRequest<'static> {
    matrix::MatrixGammaRequest::new(
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
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::from_ints(map),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    })
}
const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const ANTI_TRANSLATION: [[i32; 3]; 3] = [[-1, -1, 0], [-1, 1, 0], [0, 0, -1]];

fn close(a: C, b: C) -> bool {
    let error = (a - b).norm();
    error.is_finite() && error <= 1e-8
}

// Deliberately do not call the production apply/solver/reconstruction helpers.
fn act(a: &[f64], d: usize, anti: bool, z: &[C]) -> Vec<C> {
    assert_eq!(a.len(), d * d);
    assert_eq!(z.len(), d);
    (0..d)
        .map(|r| {
            (0..d)
                .map(|c| a[r * d + c] * if anti { z[c].conj() } else { z[c] })
                .sum()
        })
        .collect()
}
fn lift(u: &[C], d: usize, z: &[C]) -> Vec<C> {
    assert_eq!(u.len(), 2 * d);
    (0..2)
        .map(|r| (0..d).map(|c| u[r * d + c] * z[c]).sum())
        .collect()
}
fn reconstructs(
    model: &matrix::MatrixGammaModel,
    piece: &matrix::MatrixSubductionPiece,
    u: &[C],
) -> bool {
    let d = piece.dimension() as usize;
    let probes = [
        vec![C::new(0.31, 0.73); d],
        (0..d).map(|i| C::new(-0.51 + i as f64, -0.29)).collect(),
    ];
    for (g, source) in model.source_operators().iter().enumerate() {
        let target = &piece.operators()[g];
        for z in &probes {
            let lhs = act(
                &source.entries(),
                2,
                source.is_antiunitary(),
                &lift(u, d, z),
            );
            let rhs = lift(u, d, &act(target.entries(), d, target.is_antiunitary(), z));
            if !lhs.iter().zip(rhs).all(|(&a, b)| close(a, b)) {
                return false;
            }
        }
    }
    true
}

fn check_all(model: &matrix::MatrixGammaModel, result: &matrix::MatrixSubduction<'_>) {
    assert!(std::ptr::eq(model, result.model()));
    assert_eq!(result.source_dimension(), 2);
    assert_eq!(result.covered_dimension(), 2);
    let classes = model.classes();
    let n = classes.len();
    let products: Vec<Vec<usize>> = classes
        .iter()
        .map(|g| {
            classes
                .iter()
                .map(|h| {
                    let gh = compose(g, h).expect("exact operation product");
                    classes
                        .iter()
                        .position(|c| {
                            c.rotation == gh.rotation
                                && c.time_reversal == gh.time_reversal
                                && canonical_translation(c.translation, model.translation_lattice())
                                    == canonical_translation(
                                        gh.translation,
                                        model.translation_lattice(),
                                    )
                        })
                        .expect("complete M class list")
                })
                .collect()
        })
        .collect();
    // Source law and public apply agree with independent dense semilinear action.
    let z = vec![C::new(0.37, 0.81), C::new(-0.23, 0.41)];
    for g in 0..n {
        let s = &model.source_operators()[g];
        assert_eq!(s.is_antiunitary(), classes[g].time_reversal);
        assert!(s.entries().iter().all(|x| x.is_finite()));
        let entries = s.entries();
        for r in 0..2 {
            for c in 0..2 {
                let gram: f64 = (0..2)
                    .map(|k| entries[2 * k + r] * entries[2 * k + c])
                    .sum();
                assert!(close(
                    C::new(gram, 0.0),
                    C::new(if r == c { 1.0 } else { 0.0 }, 0.0)
                ));
            }
        }
        let dense = act(&s.entries(), 2, s.is_antiunitary(), &z);
        assert!(
            s.apply([z[0], z[1]])
                .iter()
                .zip(&dense)
                .all(|(&a, &b)| close(a, b))
        );
        for (h, t) in model.source_operators().iter().enumerate() {
            let p = &model.source_operators()[products[g][h]];
            assert_eq!(s.is_antiunitary() ^ t.is_antiunitary(), p.is_antiunitary());
            let lhs = act(
                &s.entries(),
                2,
                s.is_antiunitary(),
                &act(&t.entries(), 2, t.is_antiunitary(), &z),
            );
            let rhs = act(&p.entries(), 2, p.is_antiunitary(), &z);
            assert!(
                lhs.iter().zip(rhs).all(|(&a, b)| close(a, b)),
                "source pair {g},{h}"
            );
        }
    }
    let mut columns = Vec::new();
    for piece in result.pieces() {
        let d = piece.dimension() as usize;
        assert_eq!(piece.copies().len(), piece.multiplicity() as usize);
        assert_eq!(piece.operators().len(), n);
        for g in 0..n {
            let a = &piece.operators()[g];
            assert_eq!(a.dimension(), d);
            assert_eq!(a.is_antiunitary(), classes[g].time_reversal);
            assert_eq!(a.entries().len(), d * d);
            assert!(a.entries().iter().all(|x| x.is_finite()));
            for (h, b) in piece.operators().iter().enumerate() {
                let p = &piece.operators()[products[g][h]];
                for c in 0..d {
                    let mut v = vec![C::new(0.0, 0.0); d];
                    v[c] = C::new(0.29, 0.67);
                    let lhs = act(
                        a.entries(),
                        d,
                        a.is_antiunitary(),
                        &act(b.entries(), d, b.is_antiunitary(), &v),
                    );
                    let rhs = act(p.entries(), d, p.is_antiunitary(), &v);
                    assert!(
                        lhs.iter().zip(rhs).all(|(&x, y)| close(x, y)),
                        "target pair {g},{h}"
                    );
                }
            }
        }
        for copy in piece.copies() {
            assert!(copy.iter().all(|x| x.re.is_finite() && x.im.is_finite()));
            assert!(copy.iter().any(|x| x.norm() > 1e-9));
            assert!(reconstructs(model, piece, copy));
            for c in 0..d {
                columns.push([copy[c], copy[d + c]]);
            }
        }
    }
    assert_eq!(
        columns.len(),
        2,
        "actual copies, not only metadata dimensions"
    );
    // Independent V†V and VV† tests, including cross-target/cross-copy blocks.
    for i in 0..2 {
        for j in 0..2 {
            let expected = C::new(if i == j { 1.0 } else { 0.0 }, 0.0);
            let gram: C = (0..2).map(|r| columns[i][r].conj() * columns[j][r]).sum();
            let projector: C = (0..2).map(|k| columns[k][i] * columns[k][j].conj()).sum();
            assert!(close(gram, expected));
            assert!(close(projector, expected));
        }
    }
}

#[test]
fn frozen_witnesses_have_independent_laws_copies_and_whole_space_reconstruction() {
    let witnesses = [
        (89, "GM5", "GM1", 0, IDENTITY, 747, vec![(89, "GM5", 2, 1)]),
        (
            89,
            "GM5",
            "M1",
            0,
            ANTI_TRANSLATION,
            753,
            vec![(89, "GM5", 2, 1)],
        ),
        (
            89,
            "GM5",
            "GM5",
            0,
            IDENTITY,
            101,
            vec![(3, "GM1", 1, 1), (3, "GM2", 1, 1)],
        ),
        (89, "GM5", "GM5", 2, IDENTITY, 10, vec![(1, "GM1", 1, 2)]),
        (156, "GM3", "GM3", 2, IDENTITY, 1, vec![(1, "GM1", 1, 2)]),
        (
            156,
            "GM3",
            "GM1",
            0,
            IDENTITY,
            1279,
            vec![(156, "GM3", 2, 1)],
        ),
        (
            191,
            "GM5+",
            "GM1-",
            0,
            IDENTITY,
            1471,
            vec![(177, "GM5", 2, 1)],
        ),
    ];
    let mut pairs = 0;
    for (sg, src, cond, row, map, uni, mut expected) in witnesses {
        expected.sort();
        let plus =
            matrix::verify_matrix_gamma_input(&request(sg, src, cond, row, map, RealSign::Plus))
                .unwrap();
        let minus =
            matrix::verify_matrix_gamma_input(&request(sg, src, cond, row, map, RealSign::Minus))
                .unwrap();
        let results = [
            matrix::subduce_matrix_gamma(&plus).unwrap(),
            matrix::subduce_matrix_gamma(&minus).unwrap(),
        ];
        for result in &results {
            let model = result.model();
            assert_eq!(model.uni(), uni);
            assert_eq!(model.source_label(), src);
            assert_eq!(model.record_label(), cond);
            assert_eq!(model.direction_index(), row);
            assert!(model.realises_record_lattice());
            check_all(model, result);
            if uni == 1471 {
                let target = &result.pieces()[0];
                assert_eq!(
                    target
                        .operators()
                        .iter()
                        .filter(|x| x.is_antiunitary())
                        .count(),
                    12
                );
                // A carrier's anti square is D(a^2), not always pure Theta^2=I.
                assert!(
                    target
                        .operators()
                        .iter()
                        .filter(|x| x.is_antiunitary())
                        .any(|op| {
                            let a = op.entries();
                            (0..2).any(|r| {
                                (0..2).any(|c| {
                                    let square: f64 =
                                        (0..2).map(|k| a[2 * r + k] * a[2 * k + c]).sum();
                                    !close(
                                        C::new(square, 0.0),
                                        C::new(if r == c { 1.0 } else { 0.0 }, 0.0),
                                    )
                                })
                            })
                        })
                );
            }
            let mut actual: Vec<_> = result
                .pieces()
                .iter()
                .map(|p| {
                    (
                        p.unitary_sg(),
                        p.source_label(),
                        p.dimension(),
                        p.multiplicity(),
                    )
                })
                .collect();
            actual.sort();
            assert_eq!(actual, expected);
            pairs += model.classes().len().pow(2);
        }
        for g in 0..plus.classes().len() {
            let sign = if plus.classes()[g].time_reversal {
                -1.0
            } else {
                1.0
            };
            for (a, b) in plus.source_operators()[g]
                .entries()
                .into_iter()
                .zip(minus.source_operators()[g].entries())
            {
                assert!(close(C::new(a, 0.0), C::new(sign * b, 0.0)));
            }
        }
        // Target operators are independent of parent eta (the U changes instead).
        for (a, b) in results[0].pieces().iter().zip(results[1].pieces()) {
            assert_eq!(a.label(), b.label());
            assert_eq!(a.operators(), b.operators());
        }
    }
    assert_eq!(pairs, 1906);
}

#[test]
fn complex_anti_intertwiner_is_not_reconstructible_without_conjugating_its_image() {
    let model = matrix::verify_matrix_gamma_input(&request(
        89,
        "GM5",
        "M1",
        0,
        ANTI_TRANSLATION,
        RealSign::Plus,
    ))
    .unwrap();
    let result = matrix::subduce_matrix_gamma(&model).unwrap();
    let piece = &result.pieces()[0];
    let u = &piece.copies()[0];
    assert!(u.iter().any(|x| x.im.abs() > 0.1));
    assert!(reconstructs(&model, piece, u));
    let z = vec![C::new(0.29, 0.71), C::new(-0.43, 0.67)];
    let g = model
        .classes()
        .iter()
        .position(|op| op.time_reversal)
        .unwrap();
    let s = &model.source_operators()[g];
    let t = &piece.operators()[g];
    let conjugated: Vec<_> = z.iter().map(|x| x.conj()).collect();
    let wrong = act(&s.entries(), 2, false, &lift(u, 2, &conjugated)); // conjugates z, not U z
    let right = lift(u, 2, &act(t.entries(), 2, true, &z));
    assert!(wrong.iter().zip(right).any(|(&a, b)| !close(a, b)));
    let zero = vec![C::new(0.0, 0.0); 4];
    // Zero formally solves linear intertwining equations, but cannot be a copy.
    assert!(reconstructs(&model, piece, &zero));
    assert!(zero.iter().all(|x| x.norm() == 0.0));
    let mut changed = u.clone();
    changed[0] += C::new(0.47, 0.19);
    assert!(!reconstructs(&model, piece, &changed));
    let mut nonfinite = u.clone();
    nonfinite[0] = C::new(f64::NAN, 0.0);
    assert!(!reconstructs(&model, piece, &nonfinite));
}

#[test]
fn bc_source_and_condensing_identities_normalize_without_losing_source_identity() {
    use cryspglib::irrep::query;
    let source = query::irreps_of(89).iter().find(|x| x.ml == "GM5").unwrap();
    let condensing = query::irreps_of(89).iter().find(|x| x.ml == "M1").unwrap();
    let s = source.label(LabelConvention::Bc).unwrap();
    let c = condensing.label(LabelConvention::Bc).unwrap();
    let req = matrix::MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: 89,
            label: &s,
            convention: LabelConvention::Bc,
        },
        RecordIdentity {
            parent_sg: 89,
            condensing_label: &c,
            convention: LabelConvention::Bc,
            direction_index: 0,
        },
        Some(RealSign::Minus),
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::from_ints(ANTI_TRANSLATION),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    });
    let model = matrix::verify_matrix_gamma_input(&req).unwrap();
    assert_eq!(model.source_label(), "GM5");
    assert_eq!(model.record_label(), "M1");
    assert_eq!(model.source_convention(), LabelConvention::Cdml);
    assert_eq!(model.input_convention(), LabelConvention::Bc);
    check_all(&model, &matrix::subduce_matrix_gamma(&model).unwrap());
}

#[test]
fn nontrivial_trig_source_is_not_rejected_by_bitwise_float_equality() {
    // Unlike restriction to P1, the full trigonal point group actually multiplies
    // the +/-sqrt(3)/2 entries. A valid real representation is only numerically
    // equal under these products, not bitwise equal.
    let model =
        matrix::verify_matrix_gamma_input(&request(156, "GM3", "GM1", 0, IDENTITY, RealSign::Plus))
            .expect("full trigonal Gamma source must pass its numerical matrix law");
    assert!(model.classes().len() >= 6);
    let result =
        matrix::subduce_matrix_gamma(&model).expect("full trigonal target must reconstruct");
    assert_eq!(result.covered_dimension(), 2);
    assert_eq!(result.pieces().len(), 1);
    assert_eq!(result.pieces()[0].dimension(), 2);
}
