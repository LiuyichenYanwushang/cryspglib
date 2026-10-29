use cryspglib::irrep::subduction::star::decompose::physical_line::{
    PhysicalLineKind, PhysicalLineSubduction, partition_physical_line_with_embedding,
    physical_line_source, subduce_physical_line_with_embedding,
};
use cryspglib::irrep::subduction::star::decompose::{self, FullStarError};
use cryspglib::irrep::subduction::{
    ExactSeitz, Rat, SubductionComponent, SubgroupEmbedding, Vec3R,
};
use cryspglib::irrep::w_little_characters_data::{
    LittleCharacterTable, LittleOperation, W_LITTLE_CHARACTERS,
};
use cryspglib::irrep::{LabelConvention as L, isotropy, line_monodromy, query};
use num_complex::Complex64 as C;
const I: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const C2: [[i32; 3]; 3] = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]];
fn r(n: i128, d: i128) -> Rat {
    Rat::new(n, d).unwrap()
}
fn context(label: &str) -> isotropy::IsotropySubgroup {
    isotropy::isotropy_subgroups(3, label, L::Cdml).unwrap()[0]
}
fn near(a: C, b: C) {
    assert!((a - b).norm() < 1e-9, "{a} != {b}");
}
fn phase(b: Rat, op: &ExactSeitz) -> C {
    assert!(op.rotation() == I || op.rotation() == C2);
    assert!(
        op.translation()
            .as_array()
            .iter()
            .all(|x| x.denominator() == 1)
    );
    let theta = std::f64::consts::TAU
        * b.checked_mul(op.translation().get(1))
            .unwrap()
            .principal_angle()
            .to_f64();
    C::new(theta.cos(), theta.sin())
}
// Undo the frozen child-origin shift explicitly, independently of the private
// reconstruction helper. The initial unmap is the public affine geometry API.
fn child_operation(e: &SubgroupEmbedding, p: &ExactSeitz) -> ExactSeitz {
    let h = e.transform().unmap_operation(p).unwrap();
    let mut t = *h.translation().as_array();
    for (i, x) in t.iter_mut().enumerate() {
        *x = x.checked_sub(e.child_shift().get(i)).unwrap();
        for j in 0..3 {
            *x = x
                .checked_add(
                    r(i128::from(h.rotation()[i][j]), 1)
                        .checked_mul(e.child_shift().get(j))
                        .unwrap(),
                )
                .unwrap();
        }
    }
    ExactSeitz::new(h.rotation(), Vec3R::new(t))
}
fn verify(report: &PhysicalLineSubduction<'_>) -> usize {
    assert_eq!(report.source().dimension(), 2);
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(
        report.source().components(),
        [
            SubductionComponent::RealificationSeed { irnumber: 29 },
            SubductionComponent::RealificationConjugate { irnumber: 29 }
        ]
    );
    let e = report.embedding();
    let b = *report.parameter();
    assert_eq!(
        report.wave_vectors(),
        &[
            Vec3R::new([r(0, 1), b, r(0, 1)]),
            Vec3R::new([r(0, 1), b.checked_neg().unwrap(), r(0, 1)])
        ]
    );
    assert_eq!(
        report.conjugate_stars_coincide(),
        b.checked_mul(r(2, 1)).unwrap().denominator() == 1
    );
    let (source, reconstructed) = report.reconstruction();
    assert_eq!(source.len(), 27 * e.representatives().len());
    for ((op, a), v) in report
        .verification_operations()
        .iter()
        .zip(source)
        .zip(reconstructed)
    {
        let z = phase(b, op);
        near(*a, z + z.conj());
        near(*v, *a);
    }
    // Independent translates include negative and mixed coefficients outside
    // the production 0..2 cube. Retain all arm identities at coincident k.
    let samples = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [2, 0, 0],
        [-2, 1, 0],
        [1, -2, 1],
        [0, 0, -3],
        [3, 0, 1],
        [1, 1, 1],
    ];
    let mut operations = Vec::new();
    for op in e.representatives() {
        for coefficients in samples {
            let mut t = *op.translation().as_array();
            for (row, c) in e.subgroup_lattice().rows().rows().iter().zip(coefficients) {
                for (x, v) in t.iter_mut().zip(row) {
                    *x = x.checked_add(v.checked_mul(r(c, 1)).unwrap()).unwrap();
                }
            }
            let p = ExactSeitz::new(op.rotation(), Vec3R::new(t));
            let h = child_operation(e, &p);
            operations.push((p, h));
        }
    }
    let mut terms = Vec::new();
    let mut all_arms = Vec::new();
    for (bi, block) in report.blocks().iter().enumerate() {
        all_arms.extend(block.arm_indices());
        for (ti, target) in block.targets().iter().enumerate() {
            assert_eq!(target.dimension, 1);
            assert!(target.multiplicity > 0);
            if target.ml.is_none() {
                assert!(target.irnumber.is_none());
                assert!(matches!(
                    target.component,
                    SubductionComponent::Constructed { .. }
                ));
            }
            terms.push((bi, ti, target.multiplicity));
        }
    }
    all_arms.sort_unstable();
    assert_eq!(all_arms, [0, 1]);
    // Each native complex eigenline must consume exactly one actual target copy.
    // This handles both distinct targets and the endpoint's repeated scalar.
    let mut assigned = Vec::new();
    for conjugate in [false, true] {
        let slot = terms
            .iter()
            .position(|(bi, ti, m)| {
                *m > 0
                    && operations.iter().all(|(p, h)| {
                        let z = phase(b, p);
                        let z = if conjugate { z.conj() } else { z };
                        (report.blocks()[*bi].target_character(*ti, h).unwrap() - z).norm() < 1e-9
                    })
            })
            .expect("an entire native eigenline, not merely an identity character");
        terms[slot].2 -= 1;
        assigned.push((terms[slot].0, terms[slot].1));
    }
    assert!(terms.iter().all(|x| x.2 == 0));
    let s = 1.0 / 2.0_f64.sqrt();
    let u = [
        [C::new(s, 0.0), C::new(0.0, s)],
        [C::new(s, 0.0), C::new(0.0, -s)],
    ];
    for (p, h) in &operations {
        let z = phase(b, p);
        let native = [[z.re, z.im], [-z.im, z.re]];
        let mut matrix = [[C::new(0.0, 0.0); 2]; 2];
        for (column, (bi, ti)) in u.iter().zip(&assigned) {
            let value = report.blocks()[*bi].target_character(*ti, h).unwrap();
            for row in 0..2 {
                for col in 0..2 {
                    matrix[row][col] += column[row] * value * column[col].conj();
                }
            }
        }
        for row in 0..2 {
            for col in 0..2 {
                near(matrix[row][col], C::new(native[row][col], 0.0));
            }
        }
        // Check the native group law on a nontrivial translated partner too.
        let q = &operations[(operations.len() / 2 + 1) % operations.len()].0;
        near(phase(b, p) * phase(b, q), phase(b, &p.compose(q).unwrap()));
    }
    operations.len()
}

#[test]
fn physical_identity_is_one_native_seed_and_does_not_change_the_w_registry() {
    let s = physical_line_source(3, "LD1LE1").unwrap();
    assert_eq!(
        (
            s.parent_sg(),
            s.label(),
            s.pir_irnumber(),
            s.seed_irnumber(),
            s.dimension()
        ),
        (3, "LD1LE1", 29, 29, 2)
    );
    let t = s.seed_table();
    assert_eq!(
        (t.space_group, t.label, t.k_label, t.dimension),
        (3, "LD1", "LD", 1)
    );
    assert_eq!(t.direction, ["0", "1", "0"]);
    assert_eq!(t.operations.len(), 2);
    assert_eq!(t.operations[0].rotation, I.map(|a| a.map(|v| v as i8)));
    assert_eq!(t.operations[1].rotation, C2.map(|a| a.map(|v| v as i8)));
    for op in t.operations {
        assert_eq!(op.translation, ["0"; 3]);
        assert_eq!(op.character, [1, 0]);
    }
    for (sg, label) in [(3, "LD1"), (3, "LE1"), (3, "U1UA1"), (4, "LD1LE1")] {
        assert!(physical_line_source(sg, label).is_none());
    }
    assert_eq!(W_LITTLE_CHARACTERS.len(), 73);
    assert_eq!(line_monodromy::line_parents().len(), 11);
    assert!(line_monodromy::line_sources(3).is_empty());
    assert!(query::irreps_of(3).iter().all(|i| i.ml != "LD1LE1"));
}

#[test]
fn every_fixed_sg3_context_preserves_both_components_and_native_matrix_action() {
    let source = physical_line_source(3, "LD1LE1").unwrap();
    let mut calls = 0;
    let mut matrices = 0;
    let mut records = 0;
    for irrep in query::irreps_of(3).iter().filter(|r| !r.spinor) {
        for subgroup in isotropy::isotropy_subgroups(3, irrep.ml, L::Cdml).unwrap() {
            records += 1;
            let e = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
            for (n, d) in [
                (0, 1),
                (1, 7),
                (-1, 7),
                (1, 4),
                (-1, 4),
                (1, 2),
                (3, 4),
                (1, 1),
                (5, 4),
                (4_000_000_001, 4),
            ] {
                let report =
                    subduce_physical_line_with_embedding(&subgroup, &e, source, r(n, d)).unwrap();
                matrices += verify(&report);
                calls += 1;
            }
        }
    }
    assert_eq!((records, calls), (16, 160));
    println!("PHYSICAL_LINE calls={calls} native_matrices={matrices}");
}

// This intentionally wrong supplied table is legal for the old low-level
// single-complex-line API, but is NOT the physical PIR29 source.
static DOUBLE: LittleCharacterTable = LittleCharacterTable {
    space_group: 3,
    label: "synthetic-doubled-LD1",
    k_label: "synthetic",
    direction: ["0", "1", "0"],
    dimension: 2,
    operations: &[
        LittleOperation {
            element: "x,y,z",
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: ["0"; 3],
            character: [2, 0],
        },
        LittleOperation {
            element: "-x,y,-z",
            rotation: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
            translation: ["0"; 3],
            character: [2, 0],
        },
    ],
};
#[test]
fn screw_phase_and_endpoint_multiplicity_distinguish_pairing_from_doubling() {
    let source = physical_line_source(3, "LD1LE1").unwrap();
    let z2 = context("Z2");
    assert_eq!(z2.ordinal, 39);
    let e = SubgroupEmbedding::from_isotropy_subgroup(&z2).unwrap();
    let report = subduce_physical_line_with_embedding(&z2, &e, source, r(1, 4)).unwrap();
    verify(&report);
    assert!(!report.conjugate_stars_coincide());
    assert_eq!(report.blocks().len(), 1);
    let mut ids: Vec<_> = report.blocks()[0]
        .targets()
        .iter()
        .map(|t| (t.irnumber, t.ml, t.row_ml, t.multiplicity))
        .collect();
    ids.sort();
    assert_eq!(
        ids,
        [
            (Some(70), Some("Z1"), Some("Z1Z2"), 1),
            (Some(71), Some("Z2"), Some("Z1Z2"), 1)
        ]
    );
    let wrong = decompose::subduce_line_at_parameter(&z2, &e, &DOUBLE, r(1, 4)).unwrap();
    assert!(wrong.reconstruction().0.iter().any(|z| z.im.abs() > 1.9));
    assert!(report.reconstruction().0.iter().all(|z| z.im.abs() < 1e-12));
    let gm1 = context("GM1");
    let e = SubgroupEmbedding::from_isotropy_subgroup(&gm1).unwrap();
    for (b, label) in [(r(0, 1), "GM1"), (r(1, 2), "Z1"), (r(1, 1), "GM1")] {
        let v = subduce_physical_line_with_embedding(&gm1, &e, source, b).unwrap();
        verify(&v);
        assert!(v.conjugate_stars_coincide());
        assert_eq!(v.blocks().len(), 1);
        assert_eq!(v.blocks()[0].targets().len(), 1);
        assert_eq!(
            (
                v.blocks()[0].targets()[0].ml,
                v.blocks()[0].targets()[0].multiplicity
            ),
            (Some(label), 2)
        );
    }
}

#[test]
fn record_parent_and_arithmetic_errors_are_not_converted_to_unsupported() {
    let source = physical_line_source(3, "LD1LE1").unwrap();
    let a = context("GM1");
    let b = context("GM2");
    let e = SubgroupEmbedding::from_isotropy_subgroup(&a).unwrap();
    assert!(matches!(
        subduce_physical_line_with_embedding(&b, &e, source, r(1, 7)),
        Err(FullStarError::Subduction(_))
    ));
    let foreign = isotropy::isotropy_subgroups(1, "GM1", L::Cdml).unwrap()[0];
    let other = SubgroupEmbedding::from_isotropy_subgroup(&foreign).unwrap();
    assert!(matches!(
        subduce_physical_line_with_embedding(&foreign, &other, source, r(1, 7)),
        Err(FullStarError::LineSourceMismatch {
            source_sg: 3,
            parent: 1,
            ..
        })
    ));
    let error = subduce_physical_line_with_embedding(&a, &e, source, r(i128::MAX, 1)).unwrap_err();
    assert!(
        format!("{error:?}").contains("RationalOverflow"),
        "{error:?}"
    );
    let mut stale = a;
    stale.ordinal = usize::MAX;
    assert!(subduce_physical_line_with_embedding(&stale, &e, source, r(1, 7)).is_err());
}

fn fixed_contexts() -> Vec<isotropy::IsotropySubgroup> {
    query::irreps_of(3)
        .iter()
        .filter(|i| !i.spinor)
        .flat_map(|i| isotropy::isotropy_subgroups(3, i.ml, L::Cdml).unwrap())
        .collect()
}

#[test]
fn exact_physical_partition_covers_all_sixteen_contexts_without_sampling_roots() {
    use PhysicalLineKind::*;
    let source = physical_line_source(3, "LD1LE1").unwrap();
    let contexts = fixed_contexts();
    assert_eq!(contexts.len(), 16);
    let mut point_cells = 0;
    let mut intervals = 0;
    for c in &contexts {
        let e = SubgroupEmbedding::from_isotropy_subgroup(c).unwrap();
        let p = partition_physical_line_with_embedding(c, &e, source).unwrap();
        assert!(std::ptr::eq(p.source(), source));
        assert!(std::ptr::eq(p.subgroup(), c));
        assert!(std::ptr::eq(p.embedding(), &e));
        let g = if c.ordinal >= 38 { 2 } else { 1 };
        let h = if c.ordinal == 38 { 2 } else { 1 };
        assert_eq!(p.parent_collision_step(), r(1, 2));
        assert_eq!(p.momentum_collision_step(), r(1, 2 * g));
        assert_eq!(p.character_equivalence_step(), r(1, 2 * h));
        let expected: Vec<_> = (0..2 * g).map(|n| r(n, 2 * g)).collect();
        assert_eq!(
            p.boundaries()
                .iter()
                .map(|b| b.parameter)
                .collect::<Vec<_>>(),
            expected
        );
        let mut closed = expected.clone();
        closed.push(r(1, 1));
        assert_eq!(
            p.intervals(),
            closed.windows(2).map(|x| (x[0], x[1])).collect::<Vec<_>>()
        );
        point_cells += closed.len();
        intervals += p.intervals().len();
        for b in closed {
            let boundary = p.boundary(&b).unwrap();
            assert_eq!(
                boundary.parent_conjugate_stars_coincide,
                (2 * b.numerator()) % b.denominator() == 0
            );
            let kind = if c.ordinal == 39 && b.denominator() == 4 {
                DistinctCharactersAtOneMomentum
            } else {
                RepeatedCharacter
            };
            assert_eq!(boundary.kind, kind);
            assert_eq!(p.kind_at(&b), kind);
        }
        // Independent integer divisibility checks, not the production coordinate
        // solver. These supplement, rather than establish, the exact cut proof.
        for n in -194..=194 {
            let b = r(n, 97);
            let collision = (2 * n * g) % 97 == 0;
            let equal = (2 * n * h) % 97 == 0;
            let expected = if equal {
                RepeatedCharacter
            } else if collision {
                DistinctCharactersAtOneMomentum
            } else {
                DistinctMomenta
            };
            assert_eq!(p.kind_at(&b), expected);
        }
        // The old single-CIR partition has no cross-component relation at all.
        let old = cryspglib::irrep::subduction::star::line_domain::full_star_partition(
            c,
            &e,
            source.seed_table(),
        )
        .unwrap();
        assert_eq!(old.arms.len(), 1);
        assert!(old.boundaries.is_empty());
    }
    assert_eq!((point_cells, intervals), (52, 36));
}

#[test]
fn physical_partition_is_exact_near_cuts_and_for_huge_period_translates() {
    use PhysicalLineKind::*;
    let source = physical_line_source(3, "LD1LE1").unwrap();
    for label in ["Z1", "Z2"] {
        let c = context(label);
        let e = SubgroupEmbedding::from_isotropy_subgroup(&c).unwrap();
        let p = partition_physical_line_with_embedding(&c, &e, source).unwrap();
        let epsilon = r(1, 1_i128 << 90);
        for boundary in p.boundaries() {
            for b in [
                boundary.parameter.checked_add(epsilon).unwrap(),
                boundary.parameter.checked_sub(epsilon).unwrap(),
            ] {
                assert!(p.boundary(&b).is_none());
                assert_eq!(p.kind_at(&b), DistinctMomenta);
            }
        }
        let quarter = if label == "Z1" {
            RepeatedCharacter
        } else {
            DistinctCharactersAtOneMomentum
        };
        for b in [
            r(-1, 4),
            r((1_i128 << 90) + 1, 4),
            r(-((1_i128 << 90) + 1), 4),
        ] {
            assert_eq!(p.kind_at(&b), quarter);
            assert!(!p.boundary(&b).unwrap().parent_conjugate_stars_coincide);
        }
        assert_eq!(p.kind_at(&r(i128::MAX, 1)), RepeatedCharacter);
        assert!(
            p.boundary(&r(i128::MAX, 1))
                .unwrap()
                .parent_conjugate_stars_coincide
        );
        for b in [r((1_i128 << 90) + 1, 8), r(-((1_i128 << 90) + 1), 8)] {
            assert_eq!(p.kind_at(&b), DistinctMomenta);
        }
    }
}

#[test]
fn every_physical_point_cell_and_open_cell_matches_native_matrices_and_target_identities() {
    let source = physical_line_source(3, "LD1LE1").unwrap();
    let mut calls = 0;
    let mut matrices = 0;
    for c in fixed_contexts() {
        let e = SubgroupEmbedding::from_isotropy_subgroup(&c).unwrap();
        let p = partition_physical_line_with_embedding(&c, &e, source).unwrap();
        for b in p
            .boundaries()
            .iter()
            .map(|p| p.parameter)
            .chain(std::iter::once(r(1, 1)))
        {
            let v = subduce_physical_line_with_embedding(&c, &e, source, b).unwrap();
            calls += 1;
            matrices += verify(&v);
            assert_eq!(v.blocks().len(), 1);
            assert_eq!(v.blocks()[0].arm_indices(), vec![0, 1]);
            let mut got: Vec<_> = v.blocks()[0]
                .targets()
                .iter()
                .map(|t| (t.irnumber.unwrap(), t.multiplicity))
                .collect();
            got.sort();
            let at_half = b == r(1, 2);
            let at_quarter = b.denominator() == 4;
            let want = match c.ordinal {
                25 => vec![(if at_half { 7 } else { 1 }, 2)],
                30..=35 => vec![(if at_half { 85 } else { 81 }, 2)],
                38 => vec![(if at_quarter { 43 } else { 27 }, 2)],
                39 if at_quarter => vec![(70, 1), (71, 1)],
                39 => vec![(if at_half { 55 } else { 54 }, 2)],
                _ => vec![(if at_half { 43 } else { 27 }, 2)],
            };
            assert_eq!(got, want, "ordinal={} b={b}", c.ordinal);
            assert_eq!(
                v.conjugate_stars_coincide(),
                p.boundary(&b).unwrap().parent_conjugate_stars_coincide
            );
        }
        for (left, right) in p.intervals() {
            for (a, b) in [(1, 1), (2, 1), (1, 2)] {
                let t = left
                    .checked_mul(r(a, 1))
                    .unwrap()
                    .checked_add(right.checked_mul(r(b, 1)).unwrap())
                    .unwrap()
                    .checked_div(r(a + b, 1))
                    .unwrap();
                assert_eq!(p.kind_at(&t), PhysicalLineKind::DistinctMomenta);
                let v = subduce_physical_line_with_embedding(&c, &e, source, t).unwrap();
                calls += 1;
                matrices += verify(&v);
                assert_eq!(v.blocks().len(), 2);
                for block in v.blocks() {
                    assert_eq!(block.arm_indices().len(), 1);
                    assert_eq!(block.block_dimension(), 1);
                    assert_eq!(block.targets().len(), 1);
                    let target = &block.targets()[0];
                    assert_eq!(target.multiplicity, 1);
                    assert_eq!(target.irnumber, None);
                }
            }
        }
    }
    assert_eq!(calls, 160);
    println!(
        "PHYSICAL_PARTITION point_cells=52 open_cells=36 calls={calls} native_matrices={matrices}"
    );
}

#[test]
fn physical_partition_preserves_context_and_parent_refusals() {
    let source = physical_line_source(3, "LD1LE1").unwrap();
    let a = context("GM1");
    let b = context("GM2");
    let e = SubgroupEmbedding::from_isotropy_subgroup(&a).unwrap();
    assert!(partition_physical_line_with_embedding(&b, &e, source).is_err());
    let foreign = isotropy::isotropy_subgroups(1, "GM1", L::Cdml).unwrap()[0];
    let other = SubgroupEmbedding::from_isotropy_subgroup(&foreign).unwrap();
    assert!(matches!(
        partition_physical_line_with_embedding(&foreign, &other, source),
        Err(FullStarError::LineSourceMismatch {
            source_sg: 3,
            parent: 1,
            ..
        })
    ));
    let mut stale = a;
    stale.ordinal = usize::MAX;
    assert!(partition_physical_line_with_embedding(&stale, &e, source).is_err());
}
