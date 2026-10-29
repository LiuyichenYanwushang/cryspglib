//! Independent mixed-block, stored-H, full-law and complex reconstruction checks.
use cryspglib::irrep::magnetic_embedding::{
    canonical_translation_checked, compose, magnetic_operations,
};
use cryspglib::irrep::magnetic_input::{
    EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, SourceIdentity,
};
use cryspglib::irrep::magnetic_matrix_subduction::*;
use cryspglib::irrep::magnetic_summary::magnetic_irrep_summary_by_uni_partial;
use cryspglib::irrep::subduction::{Mat3R, Rat};
use cryspglib::irrep::{
    LabelConvention as LC, query,
    types::{CharacterRow, CompoundSelectedArmCharacter, IrrepRecord},
};
use num_complex::Complex64 as C;

fn rat(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}
fn request<'a>(
    sg: u8,
    source: &'a str,
    cond: &'a str,
    row: usize,
    lc: LC,
    eta: Option<RealSign>,
) -> MatrixGammaRequest<'a> {
    MatrixGammaRequest::new(
        SourceIdentity {
            parent_sg: sg,
            label: source,
            convention: lc,
        },
        RecordIdentity {
            parent_sg: sg,
            condensing_label: cond,
            convention: lc,
            direction_index: row,
        },
        eta,
    )
    .with_embedding(EmbeddingEvidence {
        map: Mat3R::identity(),
        shift: [Rat::ZERO; 3],
        lattice_policy: LatticePolicy::RequireExact,
    })
}
fn owner(sg: u8, label: &str) -> &'static IrrepRecord {
    let rows = query::find_irreps(sg, label, LC::Cdml);
    assert_eq!(rows.len(), 1);
    rows[0]
}
fn row(sg: u8, label: &str, storage: &ThreeDimensionalHStorage) -> CharacterRow {
    let record = match storage {
        ThreeDimensionalHStorage::Matrix { record }
        | ThreeDimensionalHStorage::Scalar(HCharacterStorage::Standalone { record }) => record,
        ThreeDimensionalHStorage::Scalar(HCharacterStorage::Constituent {
            record,
            component,
            irnumber,
        }) => {
            let ir = owner(sg, record.cdml);
            let meta = ir.compound_metadata().unwrap();
            assert_eq!(meta.cir_labels[*component], label);
            assert_eq!(meta.cir_irnumbers[*component], *irnumber);
            let parts = match ir.compound_selected_arm_view().unwrap() {
                CompoundSelectedArmCharacter::DistinctComponentSum { first, second, .. } => {
                    vec![first, second]
                }
                CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => vec![seed],
            };
            let part = parts
                .into_iter()
                .find(|p| p.component == *component)
                .unwrap();
            assert_eq!(part.label, label);
            return part.row;
        }
    };
    assert_eq!(record.cdml, label);
    let ir = owner(sg, record.cdml);
    assert!(ir.compound_metadata().is_none());
    ir.ordinary_scalar_selected_arm_block_trace().unwrap()
}
fn shape(r: &ThreeDimensionalSubduction<'_>) -> String {
    let mut v: Vec<_> = r
        .pieces()
        .iter()
        .map(|p| {
            format!(
                "{}{}x{}",
                match p.kind() {
                    ThreeDimensionalTargetKind::TypeA => "A",
                    ThreeDimensionalTargetKind::TypeB => "B",
                    ThreeDimensionalTargetKind::TypeC => "C",
                    ThreeDimensionalTargetKind::UnitaryPair => "UnitaryPair",
                },
                p.dimension(),
                p.multiplicity()
            )
        })
        .collect();
    v.sort();
    v.join("+")
}
fn check(r: &ThreeDimensionalSubduction<'_>) {
    let m = r.model();
    assert_eq!((r.source_dimension(), r.covered_dimension()), (3, 3));
    assert_eq!(r.scope(), THREE_DIMENSIONAL_TARGET_SCOPE);
    let columns: Vec<Vec<C>> = r
        .pieces()
        .iter()
        .flat_map(|p| {
            p.copies().iter().flat_map(move |u| {
                (0..p.dimension() as usize)
                    .map(move |j| (0..3).map(|i| u[i * p.dimension() as usize + j]).collect())
            })
        })
        .collect();
    assert_eq!(columns.len(), 3);
    for i in 0..3 {
        for j in 0..3 {
            let gram: C = (0..3).map(|k| columns[i][k].conj() * columns[j][k]).sum();
            let projector: C = (0..3).map(|k| columns[k][i] * columns[k][j].conj()).sum();
            let want = C::new(f64::from(i == j), 0.0);
            assert!((gram - want).norm() < 1e-8);
            assert!((projector - want).norm() < 1e-8);
        }
    }
    let canon = |t| canonical_translation_checked(t, m.translation_lattice()).unwrap();
    let index = |op: &cryspglib::irrep::magnetic_embedding::MagneticOperation| {
        m.classes()
            .iter()
            .position(|g| {
                g.rotation == op.rotation
                    && g.time_reversal == op.time_reversal
                    && canon(g.translation) == canon(op.translation)
            })
            .unwrap()
    };
    let raw = magnetic_operations(m.uni()).unwrap();
    let summary = magnetic_irrep_summary_by_uni_partial(m.uni()).unwrap();
    let gamma = summary
        .summary
        .kpoints
        .iter()
        .find(|k| k.coords == (0, 0, 0, 1))
        .unwrap();
    for p in r.pieces() {
        let target = gamma.coreps.iter().find(|c| c.label == p.label()).unwrap();
        let d = p.dimension() as usize;
        assert_eq!(p.h_sources().len(), p.storage().len());
        assert_eq!(p.h_sources().len(), target.source_irreps.len());
        let rows: Vec<_> = p
            .h_sources()
            .iter()
            .zip(p.storage())
            .map(|(s, st)| row(s.sg, s.ml, st))
            .collect();
        for (a, b) in p.h_sources().iter().zip(&target.source_irreps) {
            assert_eq!((a.sg, a.ml, a.dim, a.spinor), (b.sg, b.ml, b.dim, b.spinor));
        }
        for col in gamma.operations.iter().filter(|c| !c.time_reversal) {
            let i = index(&raw.operations[col.magnetic_operation_index]);
            let mut values = Vec::new();
            let mut indices = Vec::new();
            for h in &rows {
                let j = h
                    .operations()
                    .iter()
                    .position(|op| {
                        op.rotation.as_slice() == col.rotation.as_flattened()
                            && op
                                .translation
                                .iter()
                                .zip(col.translation)
                                .all(|(a, b)| ((a - b) - (a - b).round()).abs() < 1e-8)
                    })
                    .unwrap();
                indices.push(j);
                values.push(h.values()[j]);
            }
            let expected = if p.kind() == ThreeDimensionalTargetKind::TypeB {
                values[0] * 2.0
            } else {
                values.iter().sum()
            };
            let trace: C = (0..d).map(|j| p.operators()[i].entry(j, j).unwrap()).sum();
            assert!((expected - target.characters[col.column].unwrap()).norm() < 1e-8);
            assert!((trace - expected).norm() < 1e-8);
            if let ThreeDimensionalHStorage::Matrix { record } = &p.storage()[0] {
                let ir = owner(p.h_sources()[0].sg, record.cdml);
                let coefficients = &ir.matrices()[indices[0] * d * d..(indices[0] + 1) * d * d];
                for a in 0..d {
                    for b in 0..d {
                        assert!(
                            (p.operators()[i].entry(a, b).unwrap() - coefficients[a * d + b])
                                .norm()
                                < 1e-8
                        );
                    }
                }
            } else if d == 1 {
                assert!((p.operators()[i].entry(0, 0).unwrap() - values[0]).norm() < 1e-8);
            } else if p.kind() == ThreeDimensionalTargetKind::TypeB {
                for a in 0..2 {
                    for b in 0..2 {
                        assert!(
                            (p.operators()[i].entry(a, b).unwrap()
                                - if a == b { values[0] } else { C::new(0.0, 0.0) })
                            .norm()
                                < 1e-8
                        );
                    }
                }
            }
        }
        assert_eq!(p.operators().len(), m.classes().len());
        for (i, t) in p.operators().iter().enumerate() {
            assert_eq!(t.dimension(), d);
            assert_eq!(t.is_antiunitary(), m.classes()[i].time_reversal);
            assert!(t.entry(d, 0).is_none());
            assert!(t.apply(&vec![C::new(0.0, 0.0); d + 1]).is_none());
            for j in 0..m.classes().len() {
                let k = index(&compose(&m.classes()[i], &m.classes()[j]).unwrap());
                for axis in 0..d {
                    for z in [C::new(1.0, 0.0), C::new(0.0, 1.0)] {
                        let mut v = vec![C::new(0.0, 0.0); d];
                        v[axis] = z;
                        let left = t.apply(&p.operators()[j].apply(&v).unwrap()).unwrap();
                        let right = p.operators()[k].apply(&v).unwrap();
                        assert!(left.iter().zip(right).all(|(a, b)| (a - b).norm() < 1e-8));
                    }
                }
            }
            for u in p.copies() {
                assert_eq!(u.len(), 3 * d);
                for axis in 0..d {
                    for z in [C::new(1.0, 0.0), C::new(0.0, 1.0)] {
                        let mut v = vec![C::new(0.0, 0.0); d];
                        v[axis] = z;
                        let source_vector =
                            std::array::from_fn(|a| (0..d).map(|b| u[a * d + b] * v[b]).sum());
                        let left = m.source_operators()[i].apply(source_vector);
                        let target_vector = t.apply(&v).unwrap();
                        for a in 0..3 {
                            let right: C = (0..d).map(|b| u[a * d + b] * target_vector[b]).sum();
                            assert!((left[a] - right).norm() < 1e-8);
                        }
                    }
                }
            }
        }
    }
    let v = [C::new(0.3, -0.4), C::new(-0.7, 0.6), C::new(0.11, 0.23)];
    for (i, s) in m.source_operators().iter().enumerate() {
        let mut actual = [C::new(0.0, 0.0); 3];
        for p in r.pieces() {
            let d = p.dimension() as usize;
            for u in p.copies() {
                let coordinates: Vec<C> = (0..d)
                    .map(|b| (0..3).map(|a| u[a * d + b].conj() * v[a]).sum())
                    .collect();
                let image = p.operators()[i].apply(&coordinates).unwrap();
                for a in 0..3 {
                    actual[a] += (0..d).map(|b| u[a * d + b] * image[b]).sum::<C>();
                }
            }
        }
        assert!(
            actual
                .iter()
                .zip(s.apply(v))
                .all(|(a, b)| (a - b).norm() < 1e-8)
        );
    }
}

#[test]
fn all_eight_shapes_and_real_c_complex_a_keep_authoritative_sources() {
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let cases = [
        (195, "GM1", 0, identity, "A3x1"),
        (195, "GM2GM3", 0, identity, "A1x1+A1x1+A1x1"),
        (195, "GM4", 0, identity, "A1x1+A1x2"),
        (
            195,
            "GM4",
            1,
            [[0, 1, -1], [1, -1, -1], [-1, 0, -1]],
            "A1x1+UnitaryPair2x1",
        ),
        (195, "GM4", 2, identity, "A1x3"),
        (
            195,
            "R4",
            1,
            [[1, -1, 2], [0, 1, 2], [-1, 0, 2]],
            "A1x1+C2x1",
        ),
        (207, "GM3", 0, identity, "A1x1+A2x1"),
        (
            220,
            "N1",
            17,
            [[0, 0, 2], [0, 2, 0], [-2, 0, 0]],
            "A1x1+B2x1",
        ),
        (207, "GM3", 1, identity, "A1x1+C2x1"),
        (
            207,
            "GM4",
            2,
            [[1, 0, 1], [-1, 1, 1], [0, -1, 1]],
            "A1x1+A1x1+A1x1",
        ),
    ];
    let mut calls = 0;
    for (sg, cond, index, map, expected) in cases {
        for eta in [RealSign::Plus, RealSign::Minus] {
            for lc in [LC::Cdml, LC::Bc] {
                let sl = owner(sg, "GM4").label(lc).unwrap();
                let cl = owner(sg, cond).label(lc).unwrap();
                let mut req = request(sg, &sl, &cl, index, lc, Some(eta));
                req.embedding.map = Mat3R::from_ints(map);
                if sg == 220 {
                    req.embedding.shift = [rat(1, 8), Rat::ZERO, rat(3, 4)];
                }
                assert!(matches!(
                    verify_matrix_gamma_input(&req),
                    Err(MatrixSubductionError::SourceDimensionNotTwo { dimension: 3, .. })
                ));
                let m = verify_three_dimensional_gamma_input(&req).unwrap();
                let r = subduce_three_dimensional_gamma(&m)
                    .unwrap_or_else(|e| panic!("SG{sg} {cond}/{index} {eta:?}: {e:?}"));
                assert_eq!(shape(&r), expected);
                assert!(std::ptr::eq(r.model(), &m));
                check(&r);
                calls += 1;
            }
        }
    }
    assert_eq!(calls, 40);
}

#[test]
fn all_87_full_parent_targets_include_all_twenty_shifted_source_frames() {
    let mut count = 0;
    let mut shifted = 0;
    for sg in 1..=230 {
        let cond = if query::irreps_of(sg).iter().any(|r| r.ml == "GM1+") {
            "GM1+"
        } else {
            "GM1"
        };
        for s in query::irreps_of(sg).iter().filter(|r| {
            !r.spinor
                && r.compound_metadata().is_none()
                && r.dim == 3
                && r.k_vector().numerators == [0; 3]
                && r.k_vector().denominator == 1
        }) {
            let m =
                verify_three_dimensional_gamma_input(&request(sg, s.ml, cond, 0, LC::Cdml, None))
                    .unwrap();
            let r = subduce_three_dimensional_gamma(&m)
                .unwrap_or_else(|e| panic!("SG{sg} {}: {e:?}", s.ml));
            assert_eq!(shape(&r), "A3x1");
            check(&r);
            count += 1;
            if m.source_origin_shift().iter().any(|x| !x.is_zero()) {
                shifted += 1;
            }
        }
    }
    assert_eq!((count, shifted), (87, 20));
}

#[test]
fn shifted_sources_and_anti_embeddings_keep_complete_target_reconstruction() {
    let mut count = 0;
    for sg in [201, 203, 222, 224, 227, 228] {
        for source in ["GM4+", "GM4-"] {
            for eta in [RealSign::Plus, RealSign::Minus] {
                for lc in [LC::Cdml, LC::Bc] {
                    let sl = owner(sg, source).label(lc).unwrap();
                    let cl = owner(sg, "GM4+").label(lc).unwrap();
                    let mut req = request(sg, &sl, &cl, 0, lc, Some(eta));
                    if sg >= 227 {
                        req.embedding.map = Mat3R::new([
                            [Rat::ZERO, Rat::ZERO, rat(-1, 1)],
                            [rat(1, 2), rat(-1, 2), Rat::ZERO],
                            [rat(-1, 2), rat(-1, 2), Rat::ZERO],
                        ]);
                        req.embedding.shift = [
                            if sg == 227 { rat(1, 4) } else { Rat::ZERO },
                            rat(1, 4),
                            rat(1, 4),
                        ];
                    }
                    let m = verify_three_dimensional_gamma_input(&req).unwrap();
                    assert!(m.classes().iter().any(|g| g.time_reversal));
                    assert!(m.source_origin_shift().iter().any(|x| !x.is_zero()));
                    let r = subduce_three_dimensional_gamma(&m)
                        .unwrap_or_else(|e| panic!("SG{sg} {source}: {e:?}"));
                    check(&r);
                    count += 1;
                }
            }
        }
    }
    assert_eq!(count, 48);
}
