use cryspglib::irrep::magnetic_embedding::{canonical_translation_checked, compose};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::*;
use cryspglib::irrep::magnetic_subduction_api::{
    self as api, MagneticRecordSelector, MagneticSubductionRequest,
};
use cryspglib::irrep::{
    labels::LabelConvention,
    query,
    subduction::{Mat3I, Mat3R, Rat},
};
use num_complex::Complex64 as C;

const I: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
fn rat(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}
fn evidence(map: Mat3I, shift: [Rat; 3]) -> EmbeddingEvidence {
    EmbeddingEvidence {
        map: Mat3R::from_ints(map),
        shift,
        lattice_policy: LatticePolicy::RequireExact,
    }
}
fn low<'a>(
    sg: u8,
    source: &'a str,
    cond: &'a str,
    row: usize,
    eta: Option<RealSign>,
    e: EmbeddingEvidence,
) -> MatrixGammaRequest<'a> {
    MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention: LabelConvention::Cdml,
        },
        RecordIdentity {
            parent_sg: sg,
            condensing_label: cond,
            convention: LabelConvention::Cdml,
            direction_index: row,
        },
        eta,
    )
    .with_embedding(e)
}
fn check(m: &MatrixGammaModel, columns: [[C; 2]; 2], apply: impl Fn(usize, [C; 2]) -> [C; 2]) {
    assert_eq!(m.scope(), COMPOUND_SCOPE);
    let owner = query::irreps_of(m.parent_sg())
        .iter()
        .find(|r| r.ml == m.source_label())
        .unwrap();
    let meta = owner.compound_metadata().unwrap();
    for (i, p) in m.source_components().unwrap().iter().enumerate() {
        assert_eq!(
            (p.component, p.irnumber, p.label, p.dimension),
            (i, meta.cir_irnumbers[i], meta.cir_labels[i], 1)
        );
        assert_eq!(p.row.operations().len(), owner.characters().len());
        for (n, op) in p.row.operations().iter().enumerate() {
            assert_eq!(&op.rotation, &owner.pir_rotations()[n * 9..n * 9 + 9]);
            assert_eq!(&op.translation, &owner.pir_translations()[n * 3..n * 3 + 3]);
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            let gram: C = (0..2).map(|r| columns[i][r].conj() * columns[j][r]).sum();
            let proj: C = columns.iter().map(|u| u[i] * u[j].conj()).sum();
            let want = C::new(f64::from(i == j), 0.0);
            assert!((gram - want).norm() < 1e-8 && (proj - want).norm() < 1e-8);
        }
    }
    let vectors = [
        [C::new(1., 0.), C::new(0., 0.)],
        [C::new(0., 1.), C::new(0., 0.)],
        [C::new(0., 0.), C::new(1., 0.)],
        [C::new(0., 0.), C::new(0., 1.)],
        [C::new(0.3, -0.7), C::new(-0.4, 0.2)],
    ];
    let lift =
        |z: [C; 2]| std::array::from_fn::<_, 2, _>(|r| columns[0][r] * z[0] + columns[1][r] * z[1]);
    let canonical = |t| canonical_translation_checked(t, m.translation_lattice()).unwrap();
    for (i, s) in m.source_operators().iter().enumerate() {
        for z in vectors {
            let mut v = lift(z);
            if s.is_antiunitary() {
                v = v.map(|x| x.conj());
            }
            let a = s.entries();
            let left = [a[0] * v[0] + a[1] * v[1], a[2] * v[0] + a[3] * v[1]];
            let right = lift(apply(i, z));
            assert!(left.iter().zip(right).all(|(a, b)| (*a - b).norm() < 1e-8));
        }
        for (j, h) in m.classes().iter().enumerate() {
            let g = compose(&m.classes()[i], h).unwrap();
            let k = m
                .classes()
                .iter()
                .position(|x| {
                    x.rotation == g.rotation
                        && x.time_reversal == g.time_reversal
                        && canonical(x.translation) == canonical(g.translation)
                })
                .unwrap();
            for z in vectors {
                assert!(
                    apply(i, apply(j, z))
                        .iter()
                        .zip(apply(k, z))
                        .all(|(a, b)| (*a - b).norm() < 1e-8)
                );
            }
        }
    }
}
fn pair(m: &MatrixGammaModel, p: &ConjugatePairPiece) {
    assert_eq!((p.dimension(), p.multiplicity()), (2, 1));
    let u = &p.copies()[0];
    check(m, [[u[0], u[2]], [u[1], u[3]]], |i, mut z| {
        let op = &p.operators()[i];
        if op.is_antiunitary() {
            z = z.map(|x| x.conj());
        }
        let a = op.entries();
        [a[0] * z[0] + a[1] * z[1], a[2] * z[0] + a[3] * z[1]]
    });
}
fn scalar(m: &MatrixGammaModel, p: &[ScalarTargetPiece]) {
    let cols: Vec<_> = p
        .iter()
        .flat_map(|p| p.copies().iter().map(move |u| (p, u)))
        .collect();
    assert_eq!(cols.len(), 2);
    check(
        m,
        [[cols[0].1[0], cols[0].1[1]], [cols[1].1[0], cols[1].1[1]]],
        |i, z| {
            std::array::from_fn(|j| {
                let op = &cols[j].0.operators()[i];
                op.value()
                    * if op.is_antiunitary() {
                        z[j].conj()
                    } else {
                        z[j]
                    }
            })
        },
    );
}

#[test]
fn ten_shifted_compound_sources_keep_their_frames_labels_and_complete_decompositions() {
    let cases = [
        (
            85,
            ["GM3+GM4+", "GM3-GM4-"],
            "A1",
            1,
            [[2, 0, 1], [0, 0, -1], [0, 2, 0]],
            [rat(1, 4), rat(1, 4), Rat::ZERO],
            [rat(1, 4), rat(3, 4), Rat::ZERO],
        ),
        (
            86,
            ["GM3+GM4+", "GM3-GM4-"],
            "A3+A4+",
            0,
            [[0, 0, -1], [-2, 0, 1], [0, 2, 0]],
            [rat(1, 4); 3],
            [rat(1, 4); 3],
        ),
        (
            88,
            ["GM3+GM4+", "GM3-GM4-"],
            "GM3+GM4+",
            0,
            [[-1, 0, 0], [1, 0, 1], [0, 1, 0]],
            [Rat::ZERO, rat(1, 4), rat(1, 8)],
            [Rat::ZERO, rat(1, 4), rat(1, 8)],
        ),
        (
            201,
            ["GM2+GM3+", "GM2-GM3-"],
            "GM2+GM3+",
            0,
            I,
            [Rat::ZERO; 3],
            [rat(1, 4); 3],
        ),
        (
            203,
            ["GM2+GM3+", "GM2-GM3-"],
            "GM2+GM3+",
            0,
            I,
            [Rat::ZERO; 3],
            [rat(1, 8); 3],
        ),
    ];
    let mut calls = 0;
    for (sg, sources, scalar_cond, row, map, shift, origin) in cases {
        for source in sources {
            for eta in [RealSign::Plus, RealSign::Minus] {
                for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                    let owner = query::irreps_of(sg)
                        .iter()
                        .find(|x| x.ml == source)
                        .unwrap();
                    let sl = owner.label(convention).unwrap();
                    for is_pair in [false, true] {
                        let (cond, row, e) = if is_pair {
                            ("GM1+", 0, evidence(I, [Rat::ZERO; 3]))
                        } else {
                            (scalar_cond, row, evidence(map, shift))
                        };
                        let c = query::irreps_of(sg).iter().find(|x| x.ml == cond).unwrap();
                        let cl = c.label(convention).unwrap();
                        let req = MagneticSubductionRequest::new(
                            SourceIdentity {
                                parent_sg: sg,
                                label: &sl,
                                convention,
                            },
                            MagneticRecordSelector::new(sg, &cl, convention).direction(
                                cryspglib::irrep::isotropy::IsotropyDirection::Index(row),
                            ),
                            Some(eta),
                        )
                        .with_embedding(e);
                        if is_pair {
                            let out = api::subduce_compound_conjugate_pair(&req).unwrap();
                            assert_eq!(out.scope(), COMPOUND_CONJUGATE_PAIR_SCOPE);
                            assert_eq!(out.piece().kind(), ConjugatePairKind::UnitaryPair);
                            assert_eq!(out.matrix_model().source_origin_shift(), &origin);
                            assert_eq!(out.source().labels, owner.labels());
                            for h in out.piece().h_sources() {
                                assert!(h.labels().bc.is_none());
                            }
                            pair(out.matrix_model(), out.piece());
                        } else {
                            let out = api::subduce_compound_scalar_targets(&req).unwrap();
                            assert_eq!(out.scope(), COMPOUND_SCALAR_TARGET_SCOPE);
                            assert_eq!(out.matrix_model().source_origin_shift(), &origin);
                            assert_eq!(out.source().labels, owner.labels());
                            scalar(out.matrix_model(), out.pieces());
                        }
                        calls += 1;
                    }
                }
            }
        }
    }
    assert_eq!(calls, 80);
}

#[test]
fn compound_magnetic_pair_uses_the_existing_target_law() {
    for eta in [RealSign::Plus, RealSign::Minus] {
        let m = verify_compound_matrix_gamma_input(&low(
            75,
            "GM3GM4",
            "A1",
            0,
            Some(eta),
            evidence([[-1, 1, 0], [1, 1, 0], [0, 0, -2]], [Rat::ZERO; 3]),
        ))
        .unwrap();
        assert_eq!(m.uni(), 688);
        let out = subduce_matrix_gamma_conjugate_pair(&m).unwrap();
        assert_eq!(out.scope(), COMPOUND_CONJUGATE_PAIR_SCOPE);
        assert_eq!(out.piece().kind(), ConjugatePairKind::MagneticTypeC);
        pair(&m, out.piece());
    }
}

#[test]
fn b_target_gap_and_all_legacy_compound_refusals_remain_explicit() {
    for eta in [RealSign::Plus, RealSign::Minus] {
        let req = low(
            75,
            "GM3GM4",
            "GM2",
            0,
            Some(eta),
            evidence(I, [Rat::ZERO; 3]),
        );
        assert!(matches!(
            verify_matrix_gamma_input(&req),
            Err(MatrixSubductionError::SourceIsCompound { .. })
        ));
        let m = verify_compound_matrix_gamma_input(&req).unwrap();
        assert_eq!(m.uni(), 663);
        assert!(matches!(
            subduce_matrix_gamma_scalar_targets(&m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
        assert!(matches!(
            subduce_matrix_gamma_conjugate_pair(&m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
        let high = MagneticSubductionRequest::new(
            SourceIdentity {
                parent_sg: 75,
                label: "GM3GM4",
                convention: LabelConvention::Cdml,
            },
            MagneticRecordSelector::new(75, "GM2", LabelConvention::Cdml).uni(663),
            Some(eta),
        );
        for text in [
            format!("{:?}", api::subduce(&high).unwrap_err()),
            format!("{:?}", api::subduce_scalar_targets(&high).unwrap_err()),
            format!("{:?}", api::subduce_conjugate_pair(&high).unwrap_err()),
            format!("{:?}", api::subduce_type_c(&high).unwrap_err()),
        ] {
            assert!(text.contains("SourceIsCompound"), "{text}");
        }
    }
}

#[test]
fn explicit_compound_input_keeps_ordinary_spinor_non_gamma_and_missing_eta_fences() {
    let e = evidence(I, [Rat::ZERO; 3]);
    assert!(matches!(
        verify_compound_matrix_gamma_input(&low(75, "GM1", "GM2", 0, Some(RealSign::Plus), e)),
        Err(MatrixSubductionError::SourceNotCompound { .. })
    ));
    let spin = query::irreps_of(1).iter().find(|r| r.spinor).unwrap();
    assert!(matches!(
        verify_compound_matrix_gamma_input(&low(1, spin.ml, "GM1", 0, Some(RealSign::Plus), e)),
        Err(MatrixSubductionError::SourceIsSpinor { .. })
    ));
    let non = (1..=230)
        .flat_map(query::irreps_of)
        .find(|r| {
            !r.spinor
                && r.dim == 2
                && r.compound_metadata().is_some()
                && r.k_vector().numerators != [0; 3]
        })
        .unwrap();
    assert!(matches!(
        verify_compound_matrix_gamma_input(&low(non.sg, non.ml, "GM1", 0, Some(RealSign::Plus), e)),
        Err(MatrixSubductionError::SourceNotGamma { .. })
    ));
    assert!(verify_compound_matrix_gamma_input(&low(75, "GM3GM4", "GM2", 0, None, e)).is_err());
    assert!(
        verify_compound_matrix_gamma_input(&low(
            75,
            "GM3GM4",
            "GM2",
            0,
            Some(RealSign::Plus),
            evidence(I, [rat(1, 7), Rat::ZERO, Rat::ZERO])
        ))
        .is_err()
    );
}
