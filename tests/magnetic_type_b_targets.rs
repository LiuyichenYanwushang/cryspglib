//! Independent raw-M/H binding, quaternionic carrier and public boundary checks.
use cryspglib::irrep::corep::CorepType;
use cryspglib::irrep::isotropy::{IsotropyDirection, magnetic_isotropy_subgroups};
use cryspglib::irrep::magnetic_embedding::{
    MagneticOperation, canonical_translation_checked, compose, magnetic_operations,
};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::*;
use cryspglib::irrep::magnetic_subduction_api::{
    self as api, MagneticRecordSelector, MagneticSubductionRequest, MagneticTypeBReport,
};
use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_uni_partial;
use cryspglib::irrep::{
    LabelConvention, query,
    subduction::{Mat3R, Rat},
};
use num_complex::Complex64 as C;
const I: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
fn rat(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}
fn evidence(map: [[i32; 3]; 3], shift: [Rat; 3]) -> EmbeddingEvidence {
    EmbeddingEvidence {
        map: Mat3R::from_ints(map),
        shift,
        lattice_policy: LatticePolicy::RequireExact,
    }
}
fn request<'a>(
    sg: u8,
    source: &'a str,
    cond: &'a str,
    row: usize,
    convention: LabelConvention,
    eta: Option<RealSign>,
    e: EmbeddingEvidence,
) -> MagneticSubductionRequest<'a> {
    MagneticSubductionRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention,
        },
        MagneticRecordSelector::new(sg, cond, convention).direction(IsotropyDirection::Index(row)),
        eta,
    )
    .with_embedding(e)
}
fn act(a: &[f64], anti: bool, mut z: [C; 2]) -> [C; 2] {
    if anti {
        z = z.map(|v| v.conj());
    }
    [a[0] * z[0] + a[1] * z[1], a[2] * z[0] + a[3] * z[1]]
}
fn same(a: [C; 2], b: [C; 2]) {
    assert!(a.iter().zip(b).all(|(a, b)| (*a - b).norm() < 1e-8));
}
fn check(out: &MagneticTypeBReport) {
    let m = out.matrix_model();
    let p = out.piece();
    let h = out.h_source();
    assert_eq!(
        (
            out.source_dimension(),
            out.covered_dimension(),
            p.dimension(),
            p.multiplicity(),
            p.copies().len()
        ),
        (2, 2, 2, 1, 1)
    );
    assert_eq!(out.scope(), COMPOUND_TYPE_B_SCOPE);
    assert_eq!(h.dim, 1);
    assert!(!h.spinor);
    let ho = query::find_irreps(h.sg, h.ml, LabelConvention::Cdml)[0];
    assert!(matches!(out.storage(),HCharacterStorage::Standalone{record} if *record==ho.labels()));
    assert_eq!(h.labels(), ho.labels());
    assert_eq!(p.source_label(), h.ml);
    assert_eq!(p.unitary_sg(), h.sg);
    let row = ho.ordinary_scalar_selected_arm_block_trace().unwrap();
    let partial = magnetic_irrep_summary_by_uni_partial(m.uni()).unwrap();
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|k| k.coords == (0, 0, 0, 1))
        .unwrap();
    let c = gamma.coreps.iter().find(|c| c.label == p.label()).unwrap();
    assert_eq!(c.corep_type, CorepType::B);
    assert_eq!(c.source_irreps.len(), 1);
    assert_eq!(c.source_irreps[0].labels(), h.labels());
    let raw = magnetic_operations(m.uni()).unwrap();
    let canonical = |t| canonical_translation_checked(t, m.translation_lattice()).unwrap();
    let index = |g: &MagneticOperation| {
        m.classes()
            .iter()
            .position(|x| {
                x.rotation == g.rotation
                    && x.time_reversal == g.time_reversal
                    && canonical(x.translation) == canonical(g.translation)
            })
            .unwrap()
    };
    let identity = m
        .classes()
        .iter()
        .position(|g| {
            g.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                && !g.time_reversal
                && canonical(g.translation) == [Rat::ZERO; 3]
        })
        .unwrap();
    let mut chi = vec![None; m.classes().len()];
    for col in gamma.operations.iter().filter(|c| !c.time_reversal) {
        let i = index(&raw.operations[col.magnetic_operation_index]);
        let j = row
            .operations()
            .iter()
            .position(|o| {
                o.rotation.as_slice() == col.rotation.as_flattened()
                    && o.translation
                        .iter()
                        .zip(col.translation)
                        .all(|(a, b)| ((a - b) - (a - b).round()).abs() < 1e-8)
            })
            .unwrap();
        let value = row.values()[j];
        assert!(value.im.abs() < 1e-8 && (value.re.abs() - 1.).abs() < 1e-8);
        assert!((c.characters[col.column].unwrap() - 2. * value).norm() < 1e-8);
        assert!(
            p.operators()[i]
                .entries()
                .iter()
                .zip([value.re, 0., 0., value.re])
                .all(|(a, b)| (a - b).abs() < 1e-8)
        );
        if let Some(old) = chi[i] {
            assert_eq!(old, value.re);
        }
        chi[i] = Some(value.re);
    }
    assert!(
        m.classes()
            .iter()
            .enumerate()
            .all(|(i, g)| g.time_reversal || chi[i].is_some())
    );
    let a = m.classes().iter().position(|g| g.time_reversal).unwrap();
    let inv = m
        .classes()
        .iter()
        .position(|g| index(&compose(g, &m.classes()[a]).unwrap()) == identity)
        .unwrap();
    for (i, g) in m
        .classes()
        .iter()
        .enumerate()
        .filter(|(_, g)| g.time_reversal)
    {
        let h = index(&compose(&m.classes()[inv], g).unwrap());
        let x = chi[h].unwrap();
        assert!(
            p.operators()[i]
                .entries()
                .iter()
                .zip([0., -x, x, 0.])
                .all(|(a, b)| (a - b).abs() < 1e-8)
        );
        let square = index(&compose(g, g).unwrap());
        assert_eq!(chi[square], Some(-1.));
    }
    let u = &p.copies()[0];
    assert_eq!(u.len(), 4);
    for i in 0..2 {
        for j in 0..2 {
            let gram: C = (0..2).map(|r| u[2 * r + i].conj() * u[2 * r + j]).sum();
            let projector: C = (0..2).map(|r| u[2 * i + r] * u[2 * j + r].conj()).sum();
            let expected = C::new(if i == j { 1. } else { 0. }, 0.);
            assert!((gram - expected).norm() < 1e-8);
            assert!((projector - expected).norm() < 1e-8);
        }
    }
    let lift = |z: [C; 2]| [u[0] * z[0] + u[1] * z[1], u[2] * z[0] + u[3] * z[1]];
    let action = |i: usize, z| {
        act(
            p.operators()[i].entries(),
            p.operators()[i].is_antiunitary(),
            z,
        )
    };
    let vectors = [
        [C::new(1., 0.), C::new(0., 0.)],
        [C::new(0., 1.), C::new(0., 0.)],
        [C::new(0., 0.), C::new(1., 0.)],
        [C::new(0., 0.), C::new(0., 1.)],
        [C::new(0.3, -0.7), C::new(-0.4, 0.2)],
    ];
    assert_eq!(p.operators().len(), m.classes().len());
    for (i, g) in m.classes().iter().enumerate() {
        assert_eq!(p.operators()[i].is_antiunitary(), g.time_reversal);
        let source = &m.source_operators()[i];
        for z in vectors {
            same(
                act(&source.entries(), source.is_antiunitary(), lift(z)),
                lift(action(i, z)),
            );
        }
        for (j, h) in m.classes().iter().enumerate() {
            let k = index(&compose(g, h).unwrap());
            for z in vectors {
                same(action(i, action(j, z)), action(k, z));
            }
        }
    }
}
#[test]
fn all_twenty_b_targets_preserve_labels_frames_and_complex_reconstruction() {
    let cases = [
        (75, 663, 3),
        (76, 669, 4),
        (77, 675, 3),
        (78, 681, 4),
        (79, 687, 5),
        (80, 691, 5),
        (81, 695, 3),
        (82, 701, 5),
        (83, 705, 10),
        (84, 713, 10),
        (85, 721, 13),
        (86, 729, 13),
        (87, 737, 12),
        (88, 743, 15),
    ];
    let mut calls = 0;
    for (sg, uni, hsg) in cases {
        let cond = if sg < 83 { "GM2" } else { "GM2+" };
        let co = query::find_irreps(sg, cond, LabelConvention::Cdml)[0];
        let record = magnetic_isotropy_subgroups(sg, cond, LabelConvention::Cdml).unwrap()[0];
        let labels: &[&str] = if sg < 83 {
            &["GM3GM4"]
        } else {
            &["GM3+GM4+", "GM3-GM4-"]
        };
        for source in labels {
            let owner = query::find_irreps(sg, source, LabelConvention::Cdml)[0];
            for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
                for eta in [RealSign::Plus, RealSign::Minus] {
                    let sl = owner.label(convention).unwrap();
                    let cl = co.label(convention).unwrap();
                    let req = request(
                        sg,
                        &sl,
                        &cl,
                        0,
                        convention,
                        Some(eta),
                        evidence(I, [Rat::ZERO; 3]),
                    );
                    let out = api::subduce_compound_type_b(&req).unwrap();
                    assert_eq!(out.source().labels, owner.labels());
                    assert_eq!(out.condensing().labels, co.labels());
                    assert_eq!(out.source().input_convention, convention);
                    assert_eq!(out.condensing().input_convention, convention);
                    assert_eq!(
                        (
                            out.record().uni,
                            out.record().ordinal,
                            out.record().direction_index,
                            out.record().direction_label
                        ),
                        (uni, record.ordinal, 0, "P1")
                    );
                    assert_eq!(out.record().basis, record.record.basis);
                    assert_eq!(out.record().origin, record.record.origin);
                    assert_eq!(out.h_source().sg, hsg);
                    assert_eq!(
                        out.piece().label(),
                        if sg < 83 {
                            "GM2"
                        } else if source.contains('+') {
                            "GM2+"
                        } else {
                            "GM2-"
                        }
                    );
                    let origin = match sg {
                        85 => [rat(1, 4), rat(3, 4), Rat::ZERO],
                        86 => [rat(1, 4); 3],
                        88 => [Rat::ZERO, rat(1, 4), rat(1, 8)],
                        _ => [Rat::ZERO; 3],
                    };
                    assert_eq!(out.matrix_model().source_origin_shift(), &origin);
                    assert_eq!(out.embedding().shift, [Rat::ZERO; 3]);
                    assert_eq!(out.embedding().map, Mat3R::from_ints(I));
                    let meta = owner.compound_metadata().unwrap();
                    for (i, p) in out
                        .matrix_model()
                        .source_components()
                        .unwrap()
                        .iter()
                        .enumerate()
                    {
                        assert_eq!(
                            (p.component, p.label, p.irnumber, p.dimension),
                            (i, meta.cir_labels[i], meta.cir_irnumbers[i], 1)
                        );
                        assert!(query::find_irreps(sg, p.label, LabelConvention::Cdml).is_empty());
                    }
                    if sg < 83 {
                        assert_eq!(out.source().labels.bc.as_deref(), Some("Γ2Γ4"));
                        assert_eq!(out.condensing().labels.bc.as_deref(), Some("Γ3"));
                        assert_eq!(out.h_source().labels().bc.as_deref(), Some("Γ2"));
                    }
                    check(&out);
                    calls += 1;
                }
            }
        }
    }
    assert_eq!(calls, 80);
}
#[test]
fn nonidentity_embedding_and_nonzero_shift_keep_native_m_indices() {
    let map = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let shift = [rat(1, 2), rat(1, 2), Rat::ZERO];
    let owner = query::find_irreps(75, "GM3GM4", LabelConvention::Cdml)[0];
    let cond = query::find_irreps(75, "GM2", LabelConvention::Cdml)[0];
    for convention in [LabelConvention::Cdml, LabelConvention::Bc] {
        let sl = owner.label(convention).unwrap();
        let cl = cond.label(convention).unwrap();
        for eta in [RealSign::Plus, RealSign::Minus] {
            let old = api::subduce_compound_type_b(&request(
                75,
                &sl,
                &cl,
                0,
                convention,
                Some(eta),
                evidence(I, [Rat::ZERO; 3]),
            ))
            .unwrap();
            let out = api::subduce_compound_type_b(&request(
                75,
                &sl,
                &cl,
                0,
                convention,
                Some(eta),
                evidence(map, shift),
            ))
            .unwrap();
            assert_eq!(out.embedding().map, Mat3R::from_ints(map));
            assert_eq!(out.embedding().shift, shift);
            assert_eq!(out.matrix_model().source_origin_shift(), &[Rat::ZERO; 3]);
            assert_eq!(out.matrix_model().classes(), old.matrix_model().classes());
            for (s, t) in out
                .matrix_model()
                .source_operators()
                .iter()
                .zip(old.matrix_model().source_operators())
            {
                let sign = if s.is_antiunitary() { -1. } else { 1. };
                assert!(
                    s.entries()
                        .iter()
                        .zip(t.entries())
                        .all(|(a, b)| (a - sign * b).abs() < 1e-8)
                );
            }
            check(&out);
        }
    }
}
#[test]
fn b_is_explicit_and_source_fences_are_unchanged() {
    for eta in [RealSign::Plus, RealSign::Minus] {
        let req = request(
            75,
            "GM3GM4",
            "GM2",
            0,
            LabelConvention::Cdml,
            Some(eta),
            evidence(I, [Rat::ZERO; 3]),
        );
        let out = api::subduce_compound_type_b(&req).unwrap();
        let m = out.matrix_model();
        assert!(api::subduce(&req).is_err());
        assert!(api::subduce_type_c(&req).is_err());
        assert!(api::subduce_scalar_targets(&req).is_err());
        assert!(api::subduce_conjugate_pair(&req).is_err());
        assert!(api::subduce_type_b(&req).is_err());
        assert!(matches!(
            subduce_matrix_gamma(m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
        assert!(subduce_matrix_gamma_type_c(m).is_err());
        assert!(matches!(
            subduce_matrix_gamma_scalar_targets(m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
        assert!(matches!(
            subduce_matrix_gamma_conjugate_pair(m),
            Err(MatrixSubductionError::NoMatchingTarget { .. })
        ));
    }
    let sources = [
        query::find_irreps(75, "GM1", LabelConvention::Cdml)[0],
        query::irreps_of(75).iter().find(|r| r.spinor).unwrap(),
        query::irreps_of(75)
            .iter()
            .find(|r| r.compound_metadata().is_some() && r.k_vector().numerators != [0; 3])
            .unwrap(),
    ];
    for source in sources {
        assert!(
            api::subduce_compound_type_b(&request(
                75,
                source.ml,
                "GM2",
                0,
                LabelConvention::Cdml,
                Some(RealSign::Plus),
                evidence(I, [Rat::ZERO; 3])
            ))
            .is_err()
        );
    }
    assert!(
        api::subduce_compound_type_b(&request(
            195,
            "GM4",
            "GM1",
            0,
            LabelConvention::Cdml,
            Some(RealSign::Plus),
            evidence(I, [Rat::ZERO; 3])
        ))
        .is_err()
    );
    assert!(
        api::subduce_compound_type_b(&request(
            75,
            "GM3GM4",
            "GM2",
            0,
            LabelConvention::Cdml,
            None,
            evidence(I, [Rat::ZERO; 3])
        ))
        .is_err()
    );
    assert!(
        api::subduce_compound_type_b(&request(
            75,
            "GM3GM4",
            "GM2",
            0,
            LabelConvention::Cdml,
            Some(RealSign::Plus),
            evidence(I, [rat(1, 7), Rat::ZERO, Rat::ZERO])
        ))
        .is_err()
    );
}
#[test]
fn a_c_and_unitary_packets_are_not_misclassified_as_b() {
    for eta in [RealSign::Plus, RealSign::Minus] {
        let a = request(
            85,
            "GM3+GM4+",
            "A1",
            1,
            LabelConvention::Cdml,
            Some(eta),
            evidence(
                [[2, 0, 1], [0, 0, -1], [0, 2, 0]],
                [rat(1, 4), rat(1, 4), Rat::ZERO],
            ),
        );
        assert!(api::subduce_compound_scalar_targets(&a).is_ok());
        assert!(api::subduce_compound_type_b(&a).is_err());
        let c = request(
            75,
            "GM3GM4",
            "A1",
            0,
            LabelConvention::Cdml,
            Some(eta),
            evidence([[-1, 1, 0], [1, 1, 0], [0, 0, -2]], [Rat::ZERO; 3]),
        );
        assert_eq!(
            api::subduce_compound_conjugate_pair(&c)
                .unwrap()
                .piece()
                .kind(),
            ConjugatePairKind::MagneticTypeC
        );
        assert!(api::subduce_compound_type_b(&c).is_err());
        let u = request(
            75,
            "GM3GM4",
            "GM1",
            0,
            LabelConvention::Cdml,
            Some(eta),
            evidence(I, [Rat::ZERO; 3]),
        );
        assert_eq!(
            api::subduce_compound_conjugate_pair(&u)
                .unwrap()
                .piece()
                .kind(),
            ConjugatePairKind::UnitaryPair
        );
        assert!(api::subduce_compound_type_b(&u).is_err());
        let ordinary = request(
            149,
            "GM3",
            "K1",
            2,
            LabelConvention::Cdml,
            Some(eta),
            evidence([[-2, 1, 0], [-1, 2, 0], [0, 0, -1]], [Rat::ZERO; 3]),
        );
        assert!(api::subduce_conjugate_pair(&ordinary).is_ok());
        assert!(matches!(
            api::subduce_type_b(&ordinary),
            Err(api::MagneticSubductionApiError::Matrix(e)) if matches!(*e, MatrixSubductionError::NoMatchingTarget { .. })
        ));
    }
}
