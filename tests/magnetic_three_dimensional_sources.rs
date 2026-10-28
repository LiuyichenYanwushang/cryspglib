use cryspglib::irrep::magnetic_embedding::{self as m, MagneticOperation as Op};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::*;
use cryspglib::irrep::subduction::{Mat3R, Rat, SeitzTransform, Vec3R};
use cryspglib::irrep::{LabelConvention as LC, query};
use num_complex::Complex64 as C;
use std::error::Error;

fn rat(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}
fn snap(x: f64) -> Rat {
    let n = (x * 48.).round();
    assert!(x.is_finite() && (x * 48. - n).abs() < 1e-10);
    rat(n as i128, 48)
}
fn request<'a>(
    sg: u8,
    source: &'a str,
    cond: &'a str,
    eta: Option<RealSign>,
) -> MatrixGammaRequest<'a> {
    MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention: LC::Cdml,
        },
        RecordIdentity {
            parent_sg: sg,
            condensing_label: cond,
            convention: LC::Cdml,
            direction_index: 0,
        },
        eta,
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::identity(),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    })
}
fn mixed_request(source: &str, eta: Option<RealSign>) -> MatrixGammaRequest<'_> {
    let mut r = request(220, source, "N1", eta);
    r.record.direction_index = 17;
    r.embedding.map = Mat3R::from_ints([[0, 0, 2], [0, 2, 0], [-2, 0, 0]]);
    r.embedding.shift = [rat(1, 8), Rat::ZERO, rat(3, 4)];
    r
}
fn same(a: &Op, b: &Op, lattice: &[[Rat; 3]]) -> bool {
    a.rotation == b.rotation
        && a.time_reversal == b.time_reversal
        && m::canonical_translation_checked(a.translation, lattice).unwrap()
            == m::canonical_translation_checked(b.translation, lattice).unwrap()
}
fn transform(t: &SeitzTransform, g: &Op) -> Op {
    let r = t.map_rotation(g.rotation).unwrap();
    Op {
        rotation: r,
        translation: *t
            .map_translation(r, &Vec3R::new(g.translation))
            .unwrap()
            .as_array(),
        time_reversal: false,
    }
}
fn close(a: [C; 3], b: [C; 3]) {
    for (x, y) in a.into_iter().zip(b) {
        assert!((x - y).norm() < 1e-9, "{a:?} != {b:?}")
    }
}
fn vectors() -> [[C; 3]; 7] {
    let mut v = [[C::new(0., 0.); 3]; 7];
    for i in 0..3 {
        v[2 * i][i] = C::new(1., 0.);
        v[2 * i + 1][i] = C::new(0., 1.);
    }
    v[6] = [C::new(0.3, -0.7), C::new(-0.4, 0.2), C::new(0.8, 0.1)];
    v
}
// Independent complete-parent alignment, raw-M binding (not displayed H coordinates),
// all ordered products and complex 1/i/mixed-vector action checks.
fn check(model: &ThreeDimensionalGammaModel) {
    assert_eq!(model.source_dimension(), 3);
    assert_eq!(model.scope(), THREE_DIMENSIONAL_SOURCE_SCOPE);
    assert!(model.realises_record_lattice());
    assert_eq!(model.lattice_policy(), LatticePolicy::RequireExact);
    assert_eq!(model.source_convention(), LC::Cdml);
    assert_eq!(model.record_convention(), LC::Cdml);
    let ir = query::irreps_of(model.parent_sg())
        .iter()
        .find(|r| r.ml == model.source_label())
        .unwrap();
    let shift = SeitzTransform::new(Mat3R::identity(), Vec3R::new(*model.source_origin_shift()));
    let source: Vec<_> = ir
        .pir_rotations()
        .as_chunks::<9>()
        .0
        .iter()
        .zip(ir.pir_translations().as_chunks::<3>().0)
        .map(|(r, t)| {
            transform(
                &shift,
                &Op {
                    rotation: std::array::from_fn(|i| std::array::from_fn(|j| r[3 * i + j])),
                    translation: [snap(t[0]), snap(t[1]), snap(t[2])],
                    time_reversal: false,
                },
            )
        })
        .collect();
    let parent = query::symmetry_operations_of(model.parent_sg()).unwrap();
    let parent: Vec<_> = parent
        .operations
        .iter()
        .map(|g| Op {
            rotation: g.rotation,
            translation: g.translation.map(snap),
            time_reversal: false,
        })
        .collect();
    let pl = m::lattice_of_operations(&parent);
    let sl = m::lattice_of_operations(&source);
    for (a, b) in [(&pl, &sl), (&sl, &pl)] {
        for t in a {
            assert!(
                m::canonical_translation_checked(*t, b)
                    .unwrap()
                    .iter()
                    .all(|x| x.is_zero())
            );
        }
    }
    for g in &source {
        assert!(parent.iter().any(|h| same(g, h, &pl)))
    }
    for g in &parent {
        assert!(source.iter().any(|h| same(g, h, &pl)))
    }
    let raw = m::magnetic_operations(model.uni()).unwrap();
    assert_eq!(model.classes(), raw.classes);
    let map = SeitzTransform::new(*model.map(), Vec3R::new(*model.shift()));
    for g in &raw.operations {
        let class = model
            .classes()
            .iter()
            .position(|h| same(g, h, model.translation_lattice()))
            .unwrap();
        let image = transform(&map, g);
        let rows: Vec<_> = source
            .iter()
            .enumerate()
            .filter(|(_, h)| same(&image, h, &pl))
            .map(|(i, _)| i)
            .collect();
        assert!(!rows.is_empty());
        let op = model.source_operator_at(class).unwrap();
        assert_eq!(op.is_antiunitary(), g.time_reversal);
        let eta = if g.time_reversal {
            f64::from(model.time_reversal().unwrap().value())
        } else {
            1.
        };
        for row in rows {
            let want = &ir.matrices()[9 * row..9 * row + 9];
            for (x, y) in op.entries().iter().zip(want) {
                assert!((*x - eta * y).abs() < 1e-10)
            }
            for z in vectors() {
                let v = if g.time_reversal {
                    z.map(|x| x.conj())
                } else {
                    z
                };
                close(
                    op.apply(z),
                    std::array::from_fn(|i| (0..3).map(|j| eta * want[3 * i + j] * v[j]).sum()),
                );
            }
        }
        assert_eq!(op.entry(3, 0), None);
        assert_eq!(op.entry(0, usize::MAX), None);
    }
    for (i, a) in model.classes().iter().enumerate() {
        for (j, b) in model.classes().iter().enumerate() {
            let p = m::compose(a, b).unwrap();
            let k = model
                .classes()
                .iter()
                .position(|g| same(&p, g, model.translation_lattice()))
                .unwrap();
            for z in vectors() {
                close(
                    model.source_operators()[i].apply(model.source_operators()[j].apply(z)),
                    model.source_operators()[k].apply(z),
                );
            }
        }
    }
    assert!(model.source_operator_at(model.classes().len()).is_none());
}

#[test]
fn all_87_real_three_dimensional_sources_keep_complete_frames_and_the_old_fence() {
    let mut count = 0;
    let mut shifted = 0;
    let mut bc = 0;
    for sg in 1..=230 {
        let sources: Vec<_> = query::irreps_of(sg)
            .iter()
            .filter(|r| {
                !r.spinor
                    && r.compound_metadata().is_none()
                    && r.dim == 3
                    && r.k_vector().numerators == [0, 0, 0]
                    && r.k_vector().denominator == 1
            })
            .collect();
        if sources.is_empty() {
            continue;
        }
        // These 36 full-parent GM1/GM1+ record-0 geometries were independently
        // measured as map=I, shift=0 in the frozen R47 table; no embedding search.
        let cond = if query::irreps_of(sg).iter().any(|r| r.ml == "GM1+") {
            "GM1+"
        } else {
            "GM1"
        };
        for source in sources {
            let r = request(sg, source.ml, cond, None);
            let model = verify_three_dimensional_gamma_input(&r).unwrap();
            check(&model);
            assert_eq!(model.time_reversal(), None);
            assert!(model.classes().iter().all(|g| !g.time_reversal));
            assert!(matches!(
                verify_matrix_gamma_input(&r),
                Err(MatrixSubductionError::SourceDimensionNotTwo { dimension: 3, .. })
            ));
            count += 1;
            if model.source_origin_shift().iter().any(|x| !x.is_zero()) {
                shifted += 1;
                let sl = source.label(LC::Bc).unwrap();
                let cl = query::irreps_of(sg)
                    .iter()
                    .find(|r| r.ml == cond)
                    .unwrap()
                    .label(LC::Bc)
                    .unwrap();
                let mut b = request(sg, &sl, &cl, None);
                b.source.convention = LC::Bc;
                b.record.convention = LC::Bc;
                let other = verify_three_dimensional_gamma_input(&b).unwrap();
                assert_eq!(other.input_convention(), LC::Bc);
                assert_eq!(other.source_label(), source.ml);
                assert_eq!(other.record_label(), cond);
                assert_eq!(other.source_origin_shift(), model.source_origin_shift());
                assert_eq!(other.source_operators(), model.source_operators());
                bc += 1;
            }
        }
    }
    assert_eq!((count, shifted, bc), (87, 20, 20));
}

#[test]
fn mixed_a_plus_b_source_has_a_complete_manual_complex_reconstruction() {
    // This is an independent fixture oracle, NOT a production target decomposition API.
    let mut calls = 0;
    for label in ["GM4", "GM5"] {
        for eta in [RealSign::Plus, RealSign::Minus] {
            for lc in [LC::Cdml, LC::Bc] {
                let source = query::irreps_of(220)
                    .iter()
                    .find(|r| r.ml == label)
                    .unwrap()
                    .label(lc)
                    .unwrap();
                let cond = query::irreps_of(220)
                    .iter()
                    .find(|r| r.ml == "N1")
                    .unwrap()
                    .label(lc)
                    .unwrap();
                let mut r = mixed_request(&source, Some(eta));
                r.source.convention = lc;
                r.record.convention = lc;
                r.record.condensing_label = &cond;
                let model = verify_three_dimensional_gamma_input(&r).unwrap();
                check(&model);
                assert_eq!(model.uni(), 701);
                assert_eq!(model.direction_index(), 17);
                assert_eq!(*model.source_origin_shift(), [Rat::ZERO; 3]);
                assert_ne!(model.shift(), model.source_origin_shift());
                let a = model
                    .classes()
                    .iter()
                    .position(|g| g.time_reversal)
                    .unwrap();
                let h = model
                    .classes()
                    .iter()
                    .position(|g| {
                        !g.time_reversal && g.rotation != [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                    })
                    .unwrap();
                let q = if (label == "GM4") == (eta == RealSign::Plus) {
                    C::new(0., 1.)
                } else {
                    C::new(1., 0.)
                };
                // U_A=q*e_x and U_B=q*(e_y,e_z): separate projectors have ranks1/2,
                // while their sum, and the combined Gram matrix, are I3.
                let u: [[C; 3]; 3] = std::array::from_fn(|i| {
                    std::array::from_fn(|j| if i == j { q } else { C::new(0., 0.) })
                });
                for i in 0..3 {
                    for j in 0..3 {
                        let gram: C = (0..3).map(|k| u[k][i].conj() * u[k][j]).sum();
                        let proj: C = (0..3).map(|k| u[i][k] * u[j][k].conj()).sum();
                        assert!((gram - C::new(f64::from(i == j), 0.)).norm() < 1e-10);
                        assert_eq!(gram, proj);
                    }
                }
                let ah = m::compose(&model.classes()[a], &model.classes()[h]).unwrap();
                let ah = model
                    .classes()
                    .iter()
                    .position(|g| same(g, &ah, model.translation_lattice()))
                    .unwrap();
                for (i, op) in model.source_operators().iter().enumerate() {
                    for z in vectors() {
                        let t = if i == a {
                            [z[0].conj(), -z[2].conj(), z[1].conj()]
                        } else if i == ah {
                            [z[0].conj(), z[2].conj(), -z[1].conj()]
                        } else if i == h {
                            [z[0], -z[1], -z[2]]
                        } else {
                            z
                        };
                        close(op.apply(z.map(|x| q * x)), t.map(|x| q * x));
                    }
                }
                let y = [C::new(0., 0.), C::new(1., 0.), C::new(0., 0.)];
                close(
                    model.source_operators()[a].apply(model.source_operators()[a].apply(y)),
                    y.map(|x| -x),
                );
                // iI K has square +I; it cannot stand in for the B plane's J K.
                let wrong = |v: [C; 3]| v.map(|z| C::new(0., 1.) * z.conj());
                assert!((wrong(wrong(y))[1] + y[1]).norm() > 1.);
                calls += 1;
            }
        }
    }
    assert_eq!(calls, 8);
}

#[test]
fn source_domain_identity_eta_and_arithmetic_failures_remain_distinct() {
    for (sg, label, cond, dim) in [(195, "GM1", "GM1", 1), (221, "GM3+", "GM1+", 2)] {
        assert!(
            matches!(verify_three_dimensional_gamma_input(&request(sg,label,cond,None)),Err(ThreeDimensionalGammaError::SourceDimensionNotThree{dimension,..}) if dimension==dim)
        );
    }
    let compound =
        verify_three_dimensional_gamma_input(&request(195, "GM2GM3", "GM1", None)).unwrap_err();
    assert!(
        matches!(compound,ThreeDimensionalGammaError::Validation(ref e) if matches!(**e,MatrixSubductionError::SourceIsCompound{..}))
    );
    let spin = query::irreps_of(195).iter().find(|r| r.spinor).unwrap();
    assert!(
        matches!(verify_three_dimensional_gamma_input(&request(195,spin.ml,"GM1",None)),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::SourceIsSpinor{..}))
    );
    let non_gamma = (1..=230)
        .flat_map(query::irreps_of)
        .find(|r| {
            !r.spinor
                && r.compound_metadata().is_none()
                && r.dim == 3
                && r.k_vector().numerators != [0, 0, 0]
        })
        .unwrap();
    assert!(
        matches!(verify_three_dimensional_gamma_input(&request(non_gamma.sg,non_gamma.ml,"GM1",None)),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::SourceNotGamma{..}))
    );
    assert!(
        matches!(verify_three_dimensional_gamma_input(&mixed_request("GM4",None)),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::MissingTimeReversal{uni:701}))
    );
    let mut bad = mixed_request("GM4", Some(RealSign::Plus));
    bad.source.parent_sg = 0;
    assert!(
        matches!(verify_three_dimensional_gamma_input(&bad),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::UnknownParentSpaceGroup{sg:0}))
    );
    bad = mixed_request("GM4", Some(RealSign::Plus));
    bad.record.parent_sg = 221;
    assert!(
        matches!(verify_three_dimensional_gamma_input(&bad),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::ParentMismatch{..}))
    );
    bad = mixed_request("GM4", Some(RealSign::Plus));
    bad.record.direction_index = usize::MAX;
    assert!(
        matches!(verify_three_dimensional_gamma_input(&bad),Err(ThreeDimensionalGammaError::Validation(e)) if matches!(*e,MatrixSubductionError::DirectionOutOfRange{..}))
    );
    bad = mixed_request("GM4", Some(RealSign::Plus));
    let max = Rat::from_integer(i128::MAX);
    bad.embedding.map = Mat3R::new([
        [max, Rat::ZERO, Rat::ZERO],
        [Rat::ZERO, max, Rat::ZERO],
        [Rat::ZERO, Rat::ZERO, Rat::ONE],
    ]);
    let error = verify_three_dimensional_gamma_input(&bad).unwrap_err();
    assert!(
        matches!(error,ThreeDimensionalGammaError::Validation(ref e) if matches!(**e,MatrixSubductionError::Arithmetic{..}))
    );
    assert!(
        error
            .source()
            .unwrap()
            .downcast_ref::<MatrixSubductionError>()
            .is_some()
    );
    assert!(format!("{error:?}").to_lowercase().contains("overflow"));
}

#[test]
fn shifted_sources_keep_antiunitary_action_in_nontrivial_embeddings() {
    let mut calls = 0;
    for sg in [201, 203, 222, 224, 227, 228] {
        for label in ["GM4+", "GM4-"] {
            for eta in [RealSign::Plus, RealSign::Minus] {
                for lc in [LC::Cdml, LC::Bc] {
                    let sl = query::irreps_of(sg)
                        .iter()
                        .find(|r| r.ml == label)
                        .unwrap()
                        .label(lc)
                        .unwrap();
                    let cl = query::irreps_of(sg)
                        .iter()
                        .find(|r| r.ml == "GM4+")
                        .unwrap()
                        .label(lc)
                        .unwrap();
                    let mut r = request(sg, &sl, &cl, Some(eta));
                    r.source.convention = lc;
                    r.record.convention = lc;
                    if sg >= 227 {
                        r.embedding.map = Mat3R::new([
                            [Rat::ZERO, Rat::ZERO, rat(-1, 1)],
                            [rat(1, 2), rat(-1, 2), Rat::ZERO],
                            [rat(-1, 2), rat(-1, 2), Rat::ZERO],
                        ]);
                        r.embedding.shift = [
                            if sg == 227 { rat(1, 4) } else { Rat::ZERO },
                            rat(1, 4),
                            rat(1, 4),
                        ];
                    }
                    let model = verify_three_dimensional_gamma_input(&r).unwrap();
                    check(&model);
                    assert!(model.classes().iter().any(|g| g.time_reversal));
                    assert!(model.source_origin_shift().iter().any(|x| !x.is_zero()));
                    assert_eq!(model.source_label(), label);
                    assert_eq!(model.input_convention(), lc);
                    calls += 1;
                }
            }
        }
    }
    assert_eq!(calls, 48);
}
