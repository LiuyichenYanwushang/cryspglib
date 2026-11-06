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
    verify_native(report, 29)
}
// Exact native metadata, independently checked by the pinned archive oracle.
// (PIR label, id, CIR label, twice k0.x, twice k0.z, negative C2 point value).
const NATIVE_SOURCES: [(&str, u32, &str, i128, i128, bool); 8] = [
    ("LD1LE1", 29, "LD1", 0, 0, false),
    ("LD2LE2", 30, "LD2", 0, 0, true),
    ("U1UA1", 45, "U1", 1, 1, false),
    ("U2UA2", 46, "U2", 1, 1, true),
    ("V1VA1", 47, "V1", 0, 1, false),
    ("V2VA2", 48, "V2", 0, 1, true),
    ("W1WA1", 49, "W1", 1, 0, false),
    ("W2WA2", 50, "W2", 1, 0, true),
];
// PIR and CIR are separate namespaces: six SG4 pairs have different numbers.
const NATIVE_SCREW_SOURCES: [(&str, u32, u32, &str, i128, i128, bool); 8] = [
    ("LD1LE1", 56, 56, "LD1", 0, 0, false),
    ("LD2LE2", 57, 57, "LD2", 0, 0, true),
    ("U1UA1", 68, 72, "U1", 1, 1, false),
    ("U2UA2", 69, 73, "U2", 1, 1, true),
    ("V1VA1", 70, 74, "V1", 0, 1, false),
    ("V2VA2", 71, 75, "V2", 0, 1, true),
    ("W1WA1", 72, 76, "W1", 1, 0, false),
    ("W2WA2", 73, 77, "W2", 1, 0, true),
];
fn native_bits(irnumber: u32) -> (i128, i128, bool) {
    NATIVE_SOURCES
        .iter()
        .find(|s| s.1 == irnumber)
        .map(|s| (s.3, s.4, s.5))
        .or_else(|| {
            NATIVE_SCREW_SOURCES
                .iter()
                .find(|s| s.1 == irnumber)
                .map(|s| (s.4, s.5, s.6))
        })
        .unwrap()
}
// Do not read the implementation's offset or coefficient table. Point sign and
// the two transverse parity factors are independent, including mixed translates.
fn native_phase(irnumber: u32, b: Rat, op: &ExactSeitz) -> C {
    let z = if NATIVE_SCREW_SOURCES.iter().any(|s| s.1 == irnumber) {
        assert!(op.rotation() == I || op.rotation() == C2);
        assert_eq!(op.translation().get(0).denominator(), 1);
        assert_eq!(op.translation().get(2).denominator(), 1);
        let twice_y = op.translation().get(1).checked_mul(r(2, 1)).unwrap();
        assert_eq!(twice_y.denominator(), 1);
        assert_eq!(
            twice_y.numerator().rem_euclid(2),
            i128::from(op.rotation() == C2)
        );
        let theta = std::f64::consts::TAU
            * b.checked_mul(op.translation().get(1))
                .unwrap()
                .principal_angle()
                .to_f64();
        C::new(theta.cos(), theta.sin())
    } else {
        phase(b, op)
    };
    let (xbit, zbit, negative) = native_bits(irnumber);
    let t = op.translation();
    let odd = (xbit == 1 && t.get(0).numerator().rem_euclid(2) == 1)
        ^ (zbit == 1 && t.get(2).numerator().rem_euclid(2) == 1)
        ^ (negative && op.rotation() == C2);
    if odd { -z } else { z }
}
fn verify_native(report: &PhysicalLineSubduction<'_>, irnumber: u32) -> usize {
    assert_eq!(report.source().pir_irnumber(), irnumber);
    let cir = NATIVE_SCREW_SOURCES
        .iter()
        .find(|s| s.1 == irnumber)
        .map_or(irnumber, |s| s.2);
    assert_eq!(report.source().seed_irnumber(), cir);
    let (xbit, zbit, _) = native_bits(irnumber);
    let (xoffset, zoffset) = (r(xbit, 2), r(zbit, 2));
    assert_eq!(report.source().dimension(), 2);
    assert_eq!(report.covered_dimension(), 2);
    assert_eq!(
        report.source().components(),
        [
            SubductionComponent::RealificationSeed { irnumber: cir },
            SubductionComponent::RealificationConjugate { irnumber: cir }
        ]
    );
    let e = report.embedding();
    let b = *report.parameter();
    assert_eq!(
        report.wave_vectors(),
        &[
            Vec3R::new([xoffset, b, zoffset]),
            Vec3R::new([
                xoffset.checked_neg().unwrap(),
                b.checked_neg().unwrap(),
                zoffset.checked_neg().unwrap()
            ])
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
        let z = native_phase(irnumber, b, op);
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
                        let z = native_phase(irnumber, b, p);
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
        let z = native_phase(irnumber, b, p);
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
        near(
            native_phase(irnumber, b, p) * native_phase(irnumber, b, q),
            native_phase(irnumber, b, &p.compose(q).unwrap()),
        );
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
    let t = s.seed_table().unwrap();
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
    for (sg, label) in [(3, "LD1"), (3, "LE1"), (3, "UA1"), (5, "LD1LE1")] {
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
    fixed_contexts_for(3)
}
fn fixed_contexts_for(parent: u8) -> Vec<isotropy::IsotropySubgroup> {
    query::irreps_of(parent)
        .iter()
        .filter(|i| !i.spinor)
        .flat_map(|i| isotropy::isotropy_subgroups(parent, i.ml, L::Cdml).unwrap())
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
            source.seed_table().unwrap(),
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
#[test]
fn affine_source_has_real_identity_and_no_incomplete_zero_offset_table() {
    let source = physical_line_source(3, "U1UA1").unwrap();
    assert_eq!(
        (
            source.parent_sg(),
            source.label(),
            source.pir_irnumber(),
            source.seed_irnumber(),
            source.seed_label(),
            source.dimension()
        ),
        (3, "U1UA1", 45, 45, "U1", 2)
    );
    assert!(source.seed_table().is_none());
    assert_eq!(source.offset(), Vec3R::new([r(1, 2), r(0, 1), r(1, 2)]));
    assert_eq!(source.direction(), Vec3R::new([r(0, 1), r(1, 1), r(0, 1)]));
    let zero = physical_line_source(3, "LD1LE1").unwrap();
    assert_eq!(zero.offset(), Vec3R::new([r(0, 1); 3]));
    assert_eq!(zero.seed_label(), "LD1");
    assert!(zero.seed_table().is_some());
    assert!(physical_line_source(3, "U1").is_none());
    assert!(physical_line_source(3, "UA1").is_none());
    assert!(physical_line_source(5, "U1UA1").is_none());
    assert!(query::irreps_of(3).iter().all(|i| i.ml != "U1UA1"));
    assert!(line_monodromy::line_sources(3).is_empty());
    assert_eq!(W_LITTLE_CHARACTERS.len(), 73);
}

#[test]
fn affine_pair_reconstructs_native_matrices_in_every_fixed_context() {
    let source = physical_line_source(3, "U1UA1").unwrap();
    let mut calls = 0;
    let mut matrices = 0;
    let mut offset_sign_witnesses = 0;
    for subgroup in fixed_contexts() {
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let partition =
            partition_physical_line_with_embedding(&subgroup, &embedding, source).unwrap();
        for b in [
            r(0, 1),
            r(1, 4),
            r(1, 2),
            r(3, 4),
            r(1, 1),
            r(1, 3),
            r(-1, 4),
            r(7, 6),
            r(4_000_000_001, 4),
            r(-2_999_999_999, 3),
        ] {
            let report =
                subduce_physical_line_with_embedding(&subgroup, &embedding, source, b).unwrap();
            matrices += verify_native(&report, 45);
            calls += 1;
            for (op, ch) in report
                .verification_operations()
                .iter()
                .zip(report.reconstruction().0)
            {
                let missing_offset = phase(b, op);
                if (*ch - missing_offset - missing_offset.conj()).norm() > 1e-5 {
                    offset_sign_witnesses += 1;
                }
            }
            match partition.kind_at(&b) {
                PhysicalLineKind::DistinctMomenta => assert_eq!(report.blocks().len(), 2),
                PhysicalLineKind::DistinctCharactersAtOneMomentum => {
                    assert_eq!(report.blocks().len(), 1);
                    assert_eq!(report.blocks()[0].targets().len(), 2);
                    assert!(
                        report.blocks()[0]
                            .targets()
                            .iter()
                            .all(|t| t.multiplicity == 1)
                    );
                }
                PhysicalLineKind::RepeatedCharacter => {
                    assert_eq!(report.blocks().len(), 1);
                    assert_eq!(report.blocks()[0].targets().len(), 1);
                    assert_eq!(report.blocks()[0].targets()[0].multiplicity, 2);
                }
            }
        }
    }
    assert_eq!((calls, matrices), (160, 3100));
    assert!(offset_sign_witnesses > 100);
}
#[test]
fn affine_partition_cells_have_native_identities_not_the_zero_offset_labels() {
    let source = physical_line_source(3, "U1UA1").unwrap();
    let zero = physical_line_source(3, "LD1LE1").unwrap();
    // Independently decoded native CIR data, not inferred from this solver.
    let endpoints = [
        [("A1", 31), ("E1", 39)],
        [("U1", 4), ("R1", 6)],
        [("GM1", 27), ("Z1", 43)],
        [("GM2", 28), ("Z2", 44)],
        [("Y1", 41), ("C1", 35)],
        [("Y2", 42), ("C2", 36)],
        [("M1", 91), ("A1", 87)],
        [("M2", 92), ("A2", 88)],
        [("M1", 91), ("A1", 87)],
        [("M2", 92), ("A2", 88)],
        [("Y1", 85), ("GM1", 81)],
        [("Y2", 86), ("GM2", 82)],
        [("B1", 33), ("D1", 37)],
        [("B2", 34), ("D2", 38)],
    ];
    let mut points = 0;
    let mut intervals = 0;
    let mut calls = 0;
    let mut matrices = 0;
    for sg in fixed_contexts() {
        let e = SubgroupEmbedding::from_isotropy_subgroup(&sg).unwrap();
        let part = partition_physical_line_with_embedding(&sg, &e, source).unwrap();
        let old = partition_physical_line_with_embedding(&sg, &e, zero).unwrap();
        assert_eq!(part.boundaries(), old.boundaries());
        assert_eq!(
            part.momentum_collision_step(),
            old.momentum_collision_step()
        );
        assert_eq!(
            part.character_equivalence_step(),
            old.character_equivalence_step()
        );
        for b in part
            .boundaries()
            .iter()
            .map(|p| p.parameter)
            .chain([r(1, 1)])
        {
            let report = subduce_physical_line_with_embedding(&sg, &e, source, b).unwrap();
            matrices += verify_native(&report, 45);
            calls += 1;
            points += 1;
            let mut got: Vec<_> = report
                .blocks()
                .iter()
                .flat_map(|block| block.targets())
                .map(|t| (t.ml.unwrap(), t.irnumber.unwrap(), t.multiplicity))
                .collect();
            got.sort_unstable();
            let expected = match sg.ordinal {
                24..=37 => {
                    let (label, id) = endpoints[sg.ordinal - 24][usize::from(b == r(1, 2))];
                    vec![(label, id, 2)]
                }
                38 => {
                    if b == r(1, 4) || b == r(3, 4) {
                        vec![("E1", 39, 2)]
                    } else {
                        vec![("A1", 31, 2)]
                    }
                }
                39 => {
                    if b == r(1, 4) || b == r(3, 4) {
                        vec![("E1", 66, 1), ("E2", 67, 1)]
                    } else if b == r(1, 2) {
                        vec![("A2", 59, 2)]
                    } else {
                        vec![("A1", 58, 2)]
                    }
                }
                _ => unreachable!(),
            };
            assert_eq!(got, expected, "record{} b={b:?}", sg.ordinal);
        }
        for (left, right) in part.intervals() {
            intervals += 1;
            for numerator in [1, 2, 3] {
                let b = left
                    .checked_add(
                        right
                            .checked_sub(left)
                            .unwrap()
                            .checked_mul(r(numerator, 4))
                            .unwrap(),
                    )
                    .unwrap();
                let report = subduce_physical_line_with_embedding(&sg, &e, source, b).unwrap();
                matrices += verify_native(&report, 45);
                calls += 1;
                assert_eq!(report.blocks().len(), 2);
                assert_eq!(part.kind_at(&b), PhysicalLineKind::DistinctMomenta);
            }
        }
        for n in -194..=194 {
            let b = r(n, 97);
            assert_eq!(part.kind_at(&b), old.kind_at(&b));
        }
        for b in [
            r(i128::MIN, 1),
            r(i128::MAX, 1),
            r(1, 4).checked_add(r(1, 1i128 << 110)).unwrap(),
            r(1, 4).checked_sub(r(1, 1i128 << 110)).unwrap(),
        ] {
            assert_eq!(part.kind_at(&b), old.kind_at(&b));
        }
    }
    assert_eq!((points, intervals, calls, matrices), (52, 36, 160, 3110));
}
#[test]
fn all_scalar_beta_sources_keep_native_identity_and_genuine_zero_offset_tables() {
    for (label, id, cir, x, z, negative) in NATIVE_SOURCES {
        let source = physical_line_source(3, label).unwrap();
        assert_eq!(
            (
                source.label(),
                source.pir_irnumber(),
                source.seed_irnumber(),
                source.seed_label()
            ),
            (label, id, id, cir)
        );
        assert_eq!(source.offset(), Vec3R::new([r(x, 2), r(0, 1), r(z, 2)]));
        assert_eq!(source.direction(), Vec3R::new([r(0, 1), r(1, 1), r(0, 1)]));
        assert_eq!(source.seed_table().is_some(), x == 0 && z == 0);
        if let Some(table) = source.seed_table() {
            assert_eq!(table.label, cir);
            assert_eq!(
                table.operations[1].character,
                [if negative { -1 } else { 1 }, 0]
            );
        }
        assert!(physical_line_source(3, cir).is_none());
        assert!(physical_line_source(3, &label[cir.len()..]).is_none());
        assert!(physical_line_source(5, label).is_none());
        assert!(query::irreps_of(3).iter().all(|i| i.ml != label));
    }
    for label in ["F1", "G1", "GP1GQ1"] {
        assert!(physical_line_source(3, label).is_none());
    }
    assert!(line_monodromy::line_sources(3).is_empty());
    assert_eq!(W_LITTLE_CHARACTERS.len(), 73);
}

#[test]
fn all_eight_native_carriers_reconstruct_in_every_fixed_context() {
    let mut calls = 0;
    let mut matrices = 0;
    let mut wrong_point_sign = 0;
    for (label, id, _, _, _, negative) in NATIVE_SOURCES {
        let source = physical_line_source(3, label).unwrap();
        for sg in fixed_contexts() {
            let e = SubgroupEmbedding::from_isotropy_subgroup(&sg).unwrap();
            let part = partition_physical_line_with_embedding(&sg, &e, source).unwrap();
            for b in [
                r(0, 1),
                r(1, 4),
                r(1, 2),
                r(3, 4),
                r(1, 1),
                r(1, 3),
                r(-1, 4),
                r(7, 6),
                r(4_000_000_001, 4),
                r(-2_999_999_999, 3),
            ] {
                let report = subduce_physical_line_with_embedding(&sg, &e, source, b).unwrap();
                matrices += verify_native(&report, id);
                calls += 1;
                assert_eq!(
                    report.blocks().len() == 2,
                    part.kind_at(&b) == PhysicalLineKind::DistinctMomenta
                );
                if negative {
                    for (op, ch) in report
                        .verification_operations()
                        .iter()
                        .zip(report.reconstruction().0)
                    {
                        let positive = native_phase(id - 1, b, op);
                        wrong_point_sign +=
                            usize::from((*ch - positive - positive.conj()).norm() > 1e-5);
                    }
                }
            }
            for n in -194..=194 {
                let b = r(n, 97);
                let reduced = r(n.rem_euclid(97), 97);
                let g = if sg.ordinal >= 38 { 2 } else { 1 };
                let h = if sg.ordinal == 38 { 2 } else { 1 };
                let kind = if reduced.checked_mul(r(2 * h, 1)).unwrap().denominator() == 1 {
                    PhysicalLineKind::RepeatedCharacter
                } else if reduced.checked_mul(r(2 * g, 1)).unwrap().denominator() == 1 {
                    PhysicalLineKind::DistinctCharactersAtOneMomentum
                } else {
                    PhysicalLineKind::DistinctMomenta
                };
                assert_eq!(part.kind_at(&b), kind);
            }
            for b in [r(i128::MIN, 1), r(i128::MAX, 1)] {
                assert_eq!(part.kind_at(&b), PhysicalLineKind::RepeatedCharacter);
            }
        }
    }
    assert_eq!((calls, matrices), (1280, 24800));
    assert!(wrong_point_sign > 1000);
}
fn verify_native_boundary_fixture(
    parent: u8,
    fixture: &str,
    sources: &[(&str, u32)],
) -> (usize, usize, usize, usize) {
    use std::collections::BTreeSet;
    // Generated by an independent native-CIR/Fraction forward-affine audit,
    // not this solver. Both copies retain real owner/component provenance.
    let mut seen = BTreeSet::new();
    let mut matrices = 0;
    let mut calls = 0;
    let contexts = fixed_contexts_for(parent);
    for line in fixture.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 13);
        let id: u32 = fields[0].parse().unwrap();
        let ordinal: usize = fields[2].parse().unwrap();
        let b = if let Some((n, d)) = fields[4].split_once('/') {
            r(n.parse().unwrap(), d.parse().unwrap())
        } else {
            r(fields[4].parse().unwrap(), 1)
        };
        assert!(seen.insert((id, ordinal, b.numerator(), b.denominator())));
        let sg = contexts.iter().find(|s| s.ordinal == ordinal).unwrap();
        let e = SubgroupEmbedding::from_isotropy_subgroup(sg).unwrap();
        assert_eq!(e.subgroup_sg(), fields[3].parse::<u8>().unwrap());
        let report = subduce_physical_line_with_embedding(
            sg,
            &e,
            physical_line_source(parent, fields[1]).unwrap(),
            b,
        )
        .unwrap();
        matrices += verify_native(&report, id);
        calls += 1;
        let expected = |i: usize| {
            let irnumber = fields[i + 1].parse().unwrap();
            let component = if fields[i + 2].is_empty() {
                SubductionComponent::Ordinary
            } else {
                SubductionComponent::Constituent {
                    index: fields[i + 2].parse().unwrap(),
                    irnumber,
                }
            };
            (fields[i], irnumber, fields[i + 3], component)
        };
        let mut expected = [expected(5), expected(9)];
        expected.sort_by_key(|t| t.1);
        let mut got: Vec<_> = report
            .blocks()
            .iter()
            .flat_map(|block| block.targets())
            .flat_map(|t| {
                std::iter::repeat_n(
                    (
                        t.ml.unwrap(),
                        t.irnumber.unwrap(),
                        t.row_ml.unwrap(),
                        t.component,
                    ),
                    t.multiplicity as usize,
                )
            })
            .collect();
        got.sort_by_key(|t| t.1);
        assert_eq!(got, expected, "source{id} record{ordinal} b={b:?}");
    }
    let mut intervals = 0;
    for &(label, id) in sources {
        for sg in &contexts {
            let e = SubgroupEmbedding::from_isotropy_subgroup(sg).unwrap();
            let source = physical_line_source(parent, label).unwrap();
            let part = partition_physical_line_with_embedding(sg, &e, source).unwrap();
            for b in part
                .boundaries()
                .iter()
                .map(|p| p.parameter)
                .chain([source.parent_character_period()])
            {
                assert!(seen.contains(&(id, sg.ordinal, b.numerator(), b.denominator())));
            }
            for (left, right) in part.intervals() {
                intervals += 1;
                for n in [1, 2, 3] {
                    let b = left
                        .checked_add(
                            right
                                .checked_sub(left)
                                .unwrap()
                                .checked_mul(r(n, 4))
                                .unwrap(),
                        )
                        .unwrap();
                    let report = subduce_physical_line_with_embedding(sg, &e, source, b).unwrap();
                    matrices += verify_native(&report, id);
                    calls += 1;
                    assert_eq!(report.blocks().len(), 2);
                    assert!(
                        report
                            .blocks()
                            .iter()
                            .flat_map(|b| b.targets())
                            .all(|t| t.ml.is_none() && t.irnumber.is_none() && t.multiplicity == 1)
                    );
                }
            }
        }
    }
    (seen.len(), intervals, calls, matrices)
}

#[test]
fn all_scalar_line_boundary_cirs_and_owners_match_the_native_fixture() {
    let sources: Vec<_> = NATIVE_SOURCES.iter().map(|s| (s.0, s.1)).collect();
    assert_eq!(
        verify_native_boundary_fixture(3, include_str!("physical_line_boundaries.tsv"), &sources),
        (416, 288, 1280, 24880)
    );
}

#[test]
fn all_screw_boundary_cirs_and_owners_match_the_native_fixture() {
    let sources: Vec<_> = NATIVE_SCREW_SOURCES.iter().map(|s| (s.0, s.1)).collect();
    assert_eq!(
        verify_native_boundary_fixture(
            4,
            include_str!("physical_screw_line_boundaries.tsv"),
            &sources
        ),
        (512, 416, 1760, 27120)
    );
}

#[test]
fn screw_sources_keep_native_namespaces_and_full_parent_period() {
    for (label, pir, cir, seed, x, z, _) in NATIVE_SCREW_SOURCES {
        let source = physical_line_source(4, label).unwrap();
        assert_eq!((source.pir_irnumber(), source.seed_irnumber()), (pir, cir));
        assert_eq!((source.parent_sg(), source.seed_label()), (4, seed));
        assert_eq!(source.offset(), Vec3R::new([r(x, 2), Rat::ZERO, r(z, 2)]));
        assert_eq!(
            source.direction(),
            Vec3R::new([Rat::ZERO, Rat::ONE, Rat::ZERO])
        );
        assert_eq!(source.parent_character_period(), r(2, 1));
        assert_eq!(source.seed_table().is_some(), x == 0 && z == 0);
        assert_eq!(
            physical_line_source(3, label)
                .unwrap()
                .parent_character_period(),
            Rat::ONE
        );
    }
    for label in ["F1", "G1G1", "GP1GQ1", "LD1", "LE1", "UA1"] {
        assert!(physical_line_source(4, label).is_none());
    }
}

#[test]
fn screw_partition_uses_full_period_and_independent_integer_congruences() {
    use PhysicalLineKind::*;
    let gcd = |mut a: i128, mut b: i128| {
        a = a.abs();
        b = b.abs();
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    };
    let contexts = fixed_contexts_for(4);
    assert_eq!(contexts.len(), 12);
    let mut queries = 0;
    for c in &contexts {
        let e = SubgroupEmbedding::from_isotropy_subgroup(c).unwrap();
        let g = e.subgroup_lattice().rows().rows().iter().fold(0, |g, row| {
            assert_eq!(row[1].denominator(), 1);
            gcd(g, row[1].numerator())
        });
        assert!(g > 0);
        let a = e.representatives().iter().fold(2 * g, |a, op| {
            let y = op.translation().get(1).checked_mul(r(2, 1)).unwrap();
            assert_eq!(y.denominator(), 1);
            gcd(a, y.numerator())
        });
        assert!(a > 0);
        for (label, _, _, _, _, _, _) in NATIVE_SCREW_SOURCES {
            let source = physical_line_source(4, label).unwrap();
            let p = partition_physical_line_with_embedding(c, &e, source).unwrap();
            assert_eq!(p.parent_collision_step(), r(1, 2));
            assert_eq!(p.momentum_collision_step(), r(1, 2 * g));
            assert_eq!(p.character_equivalence_step(), r(1, a));
            assert_eq!(
                p.boundaries()
                    .iter()
                    .map(|b| b.parameter)
                    .collect::<Vec<_>>(),
                (0..4 * g).map(|n| r(n, 2 * g)).collect::<Vec<_>>()
            );
            assert_eq!(p.intervals().first().unwrap().0, Rat::ZERO);
            assert_eq!(p.intervals().last().unwrap().1, r(2, 1));
            assert_eq!(p.boundary(&Rat::ONE).unwrap().parameter, Rat::ONE);
            assert_eq!(p.boundary(&r(2, 1)).unwrap().parameter, Rat::ZERO);
            let mut samples: Vec<_> = (1..=24)
                .flat_map(|d| (-3 * d..=3 * d).map(move |n| r(n, d)))
                .collect();
            samples.extend([
                r(i128::MIN, 1),
                r(i128::MAX, 1),
                r(i128::MAX, 2),
                r(i128::MIN, 3),
                r(-1, i128::MAX),
                r(1, i128::MAX),
                r((1_i128 << 120) - 1, 1_i128 << 120),
                r((1_i128 << 120) + 1, 1_i128 << 120),
            ]);
            for b in samples {
                // Rat is reduced: m*b is integral iff its denominator divides m.
                // This oracle never multiplies the possibly extreme numerator.
                let collision = (2 * g) % b.denominator() == 0;
                let equal = a % b.denominator() == 0;
                let expected = if equal {
                    RepeatedCharacter
                } else if collision {
                    DistinctCharactersAtOneMomentum
                } else {
                    DistinctMomenta
                };
                assert_eq!(p.kind_at(&b), expected, "record{} {label} b={b}", c.ordinal);
                if let Some(point) = p.boundary(&b) {
                    let q = b.numerator() / b.denominator();
                    let rem = b.numerator() % b.denominator();
                    let floor = q - i128::from(rem < 0);
                    assert_eq!(
                        point.parameter.numerator() / point.parameter.denominator(),
                        floor.rem_euclid(2)
                    );
                    assert_eq!(
                        point.parent_conjugate_stars_coincide,
                        2 % b.denominator() == 0
                    );
                }
                queries += 1;
            }
        }
    }
    eprintln!(
        "SG4 exact classification: contexts={} queries={queries}",
        contexts.len()
    );
}

#[test]
fn screw_cohort_reconstructs_native_matrices_in_every_fixed_context() {
    let contexts = fixed_contexts_for(4);
    assert_eq!(contexts.len(), 12);
    let mut calls = 0;
    let mut matrices = 0;
    for (label, pir, _, _, _, _, _) in NATIVE_SCREW_SOURCES {
        let source = physical_line_source(4, label).unwrap();
        for c in &contexts {
            let e = SubgroupEmbedding::from_isotropy_subgroup(c).unwrap();
            let partition = partition_physical_line_with_embedding(c, &e, source).unwrap();
            let samples = [
                r(-1, 1),
                r(-3, 4),
                r(-1, 2),
                r(-1, 7),
                r(0, 1),
                r(1, 7),
                r(1, 4),
                r(1, 2),
                r(3, 4),
                r(1, 1),
                r(5, 4),
                r(3, 2),
                r(7, 4),
                r(2, 1),
                r(15, 7),
            ];
            for b in samples {
                let report = subduce_physical_line_with_embedding(c, &e, source, b).unwrap();
                matrices += verify_native(&report, pir);
                calls += 1;
                match partition.kind_at(&b) {
                    PhysicalLineKind::DistinctMomenta => assert_eq!(report.blocks().len(), 2),
                    PhysicalLineKind::DistinctCharactersAtOneMomentum => {
                        assert_eq!(report.blocks().len(), 1);
                        assert_eq!(report.blocks()[0].targets().len(), 2);
                        assert!(
                            report.blocks()[0]
                                .targets()
                                .iter()
                                .all(|t| t.multiplicity == 1)
                        );
                    }
                    PhysicalLineKind::RepeatedCharacter => {
                        assert_eq!(report.blocks().len(), 1);
                        assert_eq!(report.blocks()[0].targets().len(), 1);
                        assert_eq!(report.blocks()[0].targets()[0].multiplicity, 2);
                    }
                }
            }
        }
    }
    assert_eq!(calls, contexts.len() * 8 * 15);
    eprintln!(
        "SG4 full native reconstruction: contexts={} calls={calls} matrices={matrices}",
        contexts.len()
    );
}

#[test]
fn screw_integer_and_half_integer_parent_endpoints_are_not_identified() {
    let c = isotropy::isotropy_subgroups(4, "GM1", L::Cdml).unwrap()[0];
    let e = SubgroupEmbedding::from_isotropy_subgroup(&c).unwrap();
    assert_eq!(e.subgroup_sg(), 4);
    // Raw CIR parent endpoint pairs: integer +/-, then the compound's two halves.
    let endpoints = [
        (54, 55, 70, 71, "Z1Z2"),
        (58, 59, 66, 67, "E1E2"),
        (60, 61, 64, 65, "D1D2"),
        (68, 69, 62, 63, "C1C2"),
    ];
    for (index, (label, pir, _, _, _, _, negative)) in NATIVE_SCREW_SOURCES.into_iter().enumerate()
    {
        let source = physical_line_source(4, label).unwrap();
        let (plus, minus, first, second, owner) = endpoints[index / 2];
        let p = partition_physical_line_with_embedding(&c, &e, source).unwrap();
        for n in 0..=4 {
            let b = r(n, 2);
            let report = subduce_physical_line_with_embedding(&c, &e, source, b).unwrap();
            verify_native(&report, pir);
            assert!(report.conjugate_stars_coincide());
            let targets = report.blocks()[0].targets();
            if n % 2 == 0 {
                let id = if negative ^ (n / 2 % 2 == 1) {
                    minus
                } else {
                    plus
                };
                assert_eq!(p.kind_at(&b), PhysicalLineKind::RepeatedCharacter);
                assert_eq!(targets.len(), 1);
                assert_eq!(
                    (targets[0].irnumber, targets[0].multiplicity),
                    (Some(id), 2)
                );
            } else {
                assert_eq!(
                    p.kind_at(&b),
                    PhysicalLineKind::DistinctCharactersAtOneMomentum
                );
                let mut ids: Vec<_> = targets
                    .iter()
                    .map(|t| {
                        assert_eq!(t.row_ml, Some(owner));
                        assert_eq!(t.multiplicity, 1);
                        t.irnumber.unwrap()
                    })
                    .collect();
                ids.sort();
                assert_eq!(ids, [first, second]);
            }
        }
    }
}
