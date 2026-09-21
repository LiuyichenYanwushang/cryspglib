//! Compound-star witnesses from full CIR matrices, including individual phases.

use cryspglib::irrep::isotropy::{
    IsotropyDirection, isotropy_subgroup_for_direction, parent_primitive_basis,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
use cryspglib::irrep::subduction::star::scalar_star::ScalarStar;
use cryspglib::irrep::subduction::{
    ExactSeitz, Lattice, Mat3R, Rat, SubductionComponent, SubgroupEmbedding, Vec3R,
    strict_sg_hall_ops,
};
use cryspglib::irrep::types::{CompoundCharacterSemantics, IrrepSourceIdentity};
use num_complex::Complex64;

type RationalVector = [(i128, i128); 3];
type SourceOperation = ([i32; 9], RationalVector, (f64, f64), &'static [(f64, f64)]);

struct SourceCase {
    sg: u8,
    ml: &'static str,
    irnumber: u32,
    dimension: usize,
    arms: &'static [RationalVector],
    operations: &'static [SourceOperation],
}

include!("data/subduction_compound_cir.rs");

// Physical row -> two complex sources. The flag conjugates the WHOLE source
// character, including each arm's translation eigenvalue.
type CompoundCase = (u8, &'static str, [(u32, bool); 2]);
const COMPOUNDS: &[CompoundCase] = &[
    (19, "R1R1", [(559, false), (559, true)]),
    (23, "W1W1", [(716, false), (716, true)]),
    (45, "W1W1", [(1805, false), (1805, true)]),
    (45, "W2W2", [(1806, false), (1806, true)]),
    (46, "W1W2", [(1839, false), (1840, false)]),
    (76, "Z1Z3", [(3852, false), (3854, false)]),
    (83, "GM3+GM4+", [(4077, false), (4078, false)]),
    (167, "T1T2", [(8293, false), (8294, false)]),
    (182, "H1H2", [(8988, false), (8989, false)]),
    (198, "R2R2", [(9843, false), (9843, true)]),
    (219, "L3L3", [(10640, false), (10640, true)]),
    (220, "P1P1", [(10677, false), (10677, true)]),
];

fn vector(values: &RationalVector) -> Vec3R {
    Vec3R::new(values.map(|(n, d)| Rat::new(n, d).unwrap()))
}

fn rotation(r: &[i32; 9]) -> [[i32; 3]; 3] {
    [[r[0], r[1], r[2]], [r[3], r[4], r[5]], [r[6], r[7], r[8]]]
}

fn lattice(sg: u8) -> Lattice {
    Lattice::new(Mat3R::new(
        parent_primitive_basis(sg)
            .unwrap()
            .map(|row| row.map(|value| Rat::from_grid(value, 12).unwrap())),
    ))
    .unwrap()
}

fn source(irnumber: u32) -> &'static SourceCase {
    SOURCE_CASES
        .iter()
        .find(|s| s.irnumber == irnumber)
        .unwrap()
}

/// Full trace directly from archived diagonal arm blocks. Translation shifts
/// weight each arm separately; there is no single phase for a full-star trace.
fn archived_character(source: &SourceCase, h: &ExactSeitz, cell: &Lattice) -> Complex64 {
    let mut values = Vec::new();
    for (r, t, trace, blocks) in source.operations {
        if rotation(r) != h.rotation() {
            continue;
        }
        let delta = h.translation().checked_sub(&vector(t)).unwrap();
        if !cell.contains(&delta).unwrap() {
            continue;
        }
        assert_eq!(blocks.len(), source.arms.len());
        let raw: Complex64 = blocks.iter().map(|&(re, im)| Complex64::new(re, im)).sum();
        assert!((raw - Complex64::new(trace.0, trace.1)).norm() < 1e-12);
        let value: Complex64 = blocks
            .iter()
            .zip(source.arms)
            .map(|(&(re, im), k)| {
                let k = vector(k);
                let angle: f64 = (0..3)
                    .map(|i| k.get(i).checked_mul(delta.get(i)).unwrap().to_f64())
                    .sum();
                Complex64::new(re, im) * Complex64::from_polar(1.0, std::f64::consts::TAU * angle)
            })
            .sum();
        values.push(value);
    }
    let first = *values
        .first()
        .expect("operation must occur in the CIR fixture");
    assert!(values.iter().all(|value| (value - first).norm() < 1e-7));
    first
}

fn physical_character(case: &CompoundCase, h: &ExactSeitz, cell: &Lattice) -> Complex64 {
    case.2
        .iter()
        .map(|&(irnumber, conjugate)| {
            let value = archived_character(source(irnumber), h, cell);
            if conjugate { value.conj() } else { value }
        })
        .sum()
}

#[test]
fn individual_complex_components_and_physical_sums_match_raw_cir_matrices() {
    let (mut components_checked, mut sums_checked) = (0, 0);
    for case in COMPOUNDS {
        let &(sg, ml, sources) = case;
        let probe = query::irreps_of(sg).iter().find(|r| r.ml == ml).unwrap();
        let star = ScalarStar::new(probe).unwrap();
        let identities: Vec<_> = star
            .components()
            .iter()
            .map(|c| (c.irnumber(), c.conjugate()))
            .collect();
        assert_eq!(identities, sources, "SG{sg} {ml}");
        let cell = lattice(sg);
        let reciprocal = cell.reciprocal().unwrap();
        let hall = strict_sg_hall_ops(sg).unwrap();
        for component in star.components() {
            let fixture = source(component.irnumber());
            assert_eq!(fixture.sg, sg);
            assert_eq!(component.label(), fixture.ml);
            assert_eq!(component.full_dimension() as usize, fixture.dimension);
            assert_eq!(component.arms().len(), fixture.arms.len());
            for k in fixture.arms {
                let k = vector(k);
                let k = if component.conjugate() {
                    k.checked_neg().unwrap()
                } else {
                    k
                };
                assert!(
                    component
                        .arms()
                        .iter()
                        .any(|arm| reciprocal.same_mod(arm.wave_vector(), &k).unwrap())
                );
            }
            // Check source/Hall frame agreement before comparing any traces.
            for (r, t, _, _) in fixture.operations {
                assert!(hall.operations.iter().any(|h| h.rotation() == rotation(r)
                    && cell.same_mod(h.translation(), &vector(t)).unwrap()));
            }
        }
        assert_eq!(
            star.dimension() as usize,
            sources
                .iter()
                .map(|&(id, _)| source(id).dimension)
                .sum::<usize>()
        );
        let translations = [
            Vec3R::zero(),
            Vec3R::from_ints([2, -1, 3]),
            Vec3R::new(*cell.rows().row(0)),
            Vec3R::new(*cell.rows().row(1)),
        ];
        for (r, t, _, _) in source(sources[0].0).operations {
            for shift in &translations {
                assert!(cell.contains(shift).unwrap());
                let h = ExactSeitz::new(rotation(r), vector(t).checked_add(shift).unwrap());
                for component in star.components() {
                    let value = archived_character(source(component.irnumber()), &h, &cell);
                    let expected = if component.conjugate() {
                        value.conj()
                    } else {
                        value
                    };
                    assert!(
                        (component.character(&h).unwrap() - expected).norm() < 1e-7,
                        "SG{sg} {ml} component {:?}, operation {h:?}",
                        component.component()
                    );
                    components_checked += 1;
                }
                assert!(
                    (star.character(&h).unwrap() - physical_character(case, &h, &cell)).norm()
                        < 1e-7,
                    "SG{sg} {ml}, operation {h:?}"
                );
                sums_checked += 1;
            }
        }
    }
    assert_eq!(SOURCE_CASES.len(), 17);
    assert_eq!(
        SOURCE_CASES
            .iter()
            .map(|s| s.operations.len())
            .sum::<usize>(),
        156
    );
    assert_eq!((components_checked, sums_checked), (928, 464));
}

#[test]
fn centring_translation_distinguishes_seed_from_conjugate() {
    let probe = query::irreps_of(23)
        .iter()
        .find(|r| r.ml == "W1W1")
        .unwrap();
    let star = ScalarStar::new(probe).unwrap();
    let h = ExactSeitz::new(
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        vector(&[(1, 2), (1, 2), (1, 2)]),
    );
    assert_eq!(
        star.components()[0].component(),
        SubductionComponent::RealificationSeed { irnumber: 716 }
    );
    assert_eq!(
        star.components()[1].component(),
        SubductionComponent::RealificationConjugate { irnumber: 716 }
    );
    assert!(
        (star.components()[0].character(&h).unwrap() - Complex64::new(0.0, -1.0)).norm() < 1e-9
    );
    assert!((star.components()[1].character(&h).unwrap() - Complex64::new(0.0, 1.0)).norm() < 1e-9);
    // Checking only this real sum cannot detect a reversed Bloch sign.
    assert!(star.character(&h).unwrap().norm() < 1e-9);
}

#[test]
fn compound_restrictions_have_independently_pinned_complex_terms() {
    // (source ID, conjugate, little dimension, child star size, multiplicity).
    type Term = (u32, bool, u8, usize, u32);
    let cases: &[(u8, &str, &str, u8, &[Term])] = &[
        (19, "GM1", "R1R1", 19, &[(559, false, 2, 1, 2)]),
        (
            23,
            "GM1",
            "W1W1",
            23,
            &[(716, false, 1, 1, 1), (716, true, 1, 1, 1)],
        ),
        (45, "GM1", "W1W1", 45, &[(1805, false, 1, 2, 2)]),
        (45, "GM1", "W2W2", 45, &[(1806, false, 1, 2, 2)]),
        (
            83,
            "GM1+",
            "GM3+GM4+",
            83,
            &[(4077, false, 1, 1, 1), (4078, false, 1, 1, 1)],
        ),
        // A separate raw-CIR calculation gives T1 -> M1 and T2 -> M1.
        (167, "GM3+", "T1T2", 15, &[(363, false, 2, 1, 2)]),
    ];
    let mut operations_checked = 0;
    for &(sg, condensate, ml, child, terms) in cases {
        let subgroup =
            isotropy_subgroup_for_direction(sg, condensate, IsotropyDirection::Label("P1"))
                .unwrap();
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let probe = query::irreps_of(sg).iter().find(|r| r.ml == ml).unwrap();
        let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("SG{sg} {ml}: {error}"));
        assert_eq!(result.subgroup_sg(), child);
        assert_eq!(result.parent_irnumber(), None);
        assert_eq!(result.parent_source_identity(), probe.source_identity());
        let mut actual: Vec<Term> = result
            .blocks()
            .iter()
            .flat_map(|block| {
                block.targets().iter().map(|t| {
                    (
                        t.irnumber,
                        matches!(
                            t.component,
                            SubductionComponent::RealificationConjugate { .. }
                        ),
                        t.dimension,
                        block.star_size(),
                        t.multiplicity,
                    )
                })
            })
            .collect();
        actual.sort();
        let mut expected = terms.to_vec();
        expected.sort();
        assert_eq!(actual, expected, "SG{sg} {ml}");
        assert_eq!(result.covered_dimension(), u32::from(probe.dim));
        if sg == 23 {
            assert_eq!(result.blocks().len(), 2);
            let reciprocal = lattice(sg).reciprocal().unwrap();
            let k = vector(&[(1, 2), (1, 2), (1, 2)]);
            for block in result.blocks() {
                let expected_k = if matches!(
                    block.targets()[0].component,
                    SubductionComponent::RealificationConjugate { .. }
                ) {
                    k.checked_neg().unwrap()
                } else {
                    k
                };
                assert!(reciprocal.same_mod(block.stored_k(), &expected_k).unwrap());
            }
        }
        let case = COMPOUNDS.iter().find(|c| c.0 == sg && c.1 == ml).unwrap();
        let cell = lattice(sg);
        let (parent, reconstructed) = result.reconstruction();
        assert_eq!(result.representatives(), embedding.representatives());
        assert_eq!(parent.len(), result.representatives().len());
        assert_eq!(reconstructed.len(), parent.len());
        for ((h, actual), reconstructed) in result
            .representatives()
            .iter()
            .zip(parent)
            .zip(reconstructed)
        {
            let expected = physical_character(case, h, &cell);
            assert!(
                (actual - expected).norm() < 1e-7,
                "parent SG{sg} {ml}, {h:?}"
            );
            assert!(
                (reconstructed - expected).norm() < 1e-7,
                "child SG{sg} {ml}, {h:?}"
            );
            operations_checked += 1;
        }
    }
    assert_eq!(operations_checked, 28);
}

#[test]
fn self_restriction_preserves_every_complex_source_in_the_four_frozen_groups() {
    let mut census = Vec::new();
    for (sg, condensate) in [(19, "GM1"), (23, "GM1"), (45, "GM1"), (83, "GM1+")] {
        let subgroup =
            isotropy_subgroup_for_direction(sg, condensate, IsotropyDirection::Label("P1"))
                .unwrap();
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let cell = lattice(sg);
        let reciprocal = cell.reciprocal().unwrap();
        let hall = strict_sg_hall_ops(sg).unwrap();
        let (mut ordinary, mut compound, mut spinor) = (0, 0, 0);
        for probe in query::irreps_of(sg) {
            if probe.spinor {
                spinor += 1;
                continue;
            }
            let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
                .unwrap_or_else(|error| panic!("self restriction SG{sg} {}: {error}", probe.ml));
            let mut actual: Vec<_> = result
                .blocks()
                .iter()
                .flat_map(|block| {
                    block.targets().iter().map(|target| {
                        (
                            target.irnumber,
                            matches!(
                                target.component,
                                SubductionComponent::RealificationConjugate { .. }
                            ),
                            target.multiplicity,
                        )
                    })
                })
                .collect();
            let mut expected = match probe.source_identity() {
                IrrepSourceIdentity::OrdinaryScalar { cir_irnumber } => {
                    ordinary += 1;
                    vec![(cir_irnumber, false, 1)]
                }
                IrrepSourceIdentity::Compound { .. } => {
                    compound += 1;
                    let metadata = probe.compound_metadata().unwrap();
                    let [first, second] = metadata.cir_irnumbers;
                    match metadata.semantics {
                        CompoundCharacterSemantics::DistinctComponentSum => {
                            vec![(first, false, 1), (second, false, 1)]
                        }
                        CompoundCharacterSemantics::ConjugateRealification => {
                            // For these four source-pinned groups: #19/#45 have
                            // equivalent quaternionic seeds, #23 disjoint stars.
                            // Determine the k-star relation directly from Hall
                            // rotations, without using the scalar-star adapter.
                            let k =
                                Vec3R::new([probe.kx, probe.ky, probe.kz].map(|n| {
                                    Rat::new(i128::from(n), i128::from(probe.kd)).unwrap()
                                }));
                            let minus_k = k.checked_neg().unwrap();
                            let same_star = hall.operations.iter().any(|op| {
                                let rotation = Mat3R::from_ints(op.rotation());
                                let image = rotation
                                    .inverse()
                                    .unwrap()
                                    .transpose()
                                    .checked_mul_vector(&k)
                                    .unwrap();
                                reciprocal.same_mod(&image, &minus_k).unwrap()
                            });
                            if same_star {
                                vec![(first, false, 2)]
                            } else {
                                vec![(first, false, 1), (first, true, 1)]
                            }
                        }
                    }
                }
                IrrepSourceIdentity::Spin { .. } => unreachable!(),
            };
            actual.sort();
            expected.sort();
            assert_eq!(actual, expected, "self restriction SG{sg} {}", probe.ml);
            assert_eq!(result.covered_dimension(), u32::from(probe.dim));
        }
        census.push((sg, ordinary, compound, spinor));
    }
    assert_eq!(
        census,
        [
            (19, 7, 7, 20),
            (23, 14, 4, 9),
            (45, 10, 4, 10),
            (83, 24, 8, 40)
        ],
        "(SG, ordinary success, compound success, unsupported spinors)"
    );
}
