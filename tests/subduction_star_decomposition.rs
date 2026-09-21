//! Full-star decompositions pinned independently from archived CIR matrices.
//!
//! The expected terms below were calculated from the archived diagonal arm
//! blocks on H_q, with a separate character inner product. In particular, the
//! parent-side calculation did not use OrdinaryStar or its induction formula.

use cryspglib::irrep::isotropy::{
    IsotropyDirection, isotropy_subgroup_for_direction, parent_primitive_basis,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
use cryspglib::irrep::subduction::{ExactSeitz, Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R};
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

include!("data/subduction_star_cir.rs");

// (ML, CIR source identity, little dimension, child star size, multiplicity).
type Term = (&'static str, u32, u8, usize, u32);

struct DecompositionCase {
    sg: u8,
    condensate: &'static str,
    direction: &'static str,
    ordinal: usize,
    child: u8,
    probe: &'static str,
    terms: &'static [Term],
}

const DECOMPOSITIONS: &[DecompositionCase] = &[
    DecompositionCase {
        sg: 221,
        condensate: "GM4+",
        direction: "P1",
        ordinal: 12400,
        child: 83,
        probe: "X1+",
        terms: &[("X1+", 4099, 1, 2, 1), ("Z1+", 4115, 1, 1, 1)],
    },
    DecompositionCase {
        sg: 221,
        condensate: "GM4+",
        direction: "P1",
        ordinal: 12400,
        child: 83,
        probe: "X5+",
        terms: &[
            ("X1+", 4099, 1, 2, 1),
            ("X2+", 4100, 1, 2, 1),
            ("Z3+", 4117, 1, 1, 1),
            ("Z4+", 4118, 1, 1, 1),
        ],
    },
    DecompositionCase {
        sg: 221,
        condensate: "GM4+",
        direction: "P1",
        ordinal: 12400,
        child: 83,
        probe: "M5-",
        terms: &[
            ("M3-", 4097, 1, 1, 1),
            ("M4-", 4098, 1, 1, 1),
            ("R1-", 4113, 1, 2, 1),
            ("R2-", 4114, 1, 2, 1),
        ],
    },
    DecompositionCase {
        sg: 221,
        condensate: "GM4+",
        direction: "P2",
        ordinal: 12401,
        child: 12,
        probe: "X1+",
        terms: &[("A1+", 265, 1, 1, 1), ("V1+", 269, 1, 2, 1)],
    },
    DecompositionCase {
        sg: 221,
        condensate: "GM4+",
        direction: "P2",
        ordinal: 12401,
        child: 12,
        probe: "X5+",
        terms: &[
            ("A1+", 265, 1, 1, 1),
            ("A2+", 266, 1, 1, 1),
            ("V1+", 269, 1, 2, 2),
        ],
    },
    DecompositionCase {
        sg: 221,
        condensate: "GM4+",
        direction: "P2",
        ordinal: 12401,
        child: 12,
        probe: "M5-",
        terms: &[
            ("Y1-", 263, 1, 1, 1),
            ("Y2-", 264, 1, 1, 1),
            ("L1-", 272, 1, 2, 2),
        ],
    },
    DecompositionCase {
        sg: 139,
        condensate: "M1-",
        direction: "P1",
        ordinal: 6213,
        child: 126,
        probe: "X1+",
        terms: &[("M1", 6344, 2, 1, 1)],
    },
    DecompositionCase {
        sg: 139,
        condensate: "M1-",
        direction: "P1",
        ordinal: 6213,
        child: 126,
        probe: "N1+",
        terms: &[("R1", 6354, 2, 2, 1)],
    },
    DecompositionCase {
        sg: 139,
        condensate: "M1-",
        direction: "P1",
        ordinal: 6213,
        child: 126,
        probe: "P1",
        terms: &[("A1", 6350, 2, 1, 1)],
    },
    DecompositionCase {
        sg: 225,
        condensate: "GM4-",
        direction: "C1",
        ordinal: 13346,
        child: 8,
        probe: "L1+",
        terms: &[("V1", 154, 1, 2, 1), ("L1", 155, 1, 2, 1)],
    },
    DecompositionCase {
        sg: 225,
        condensate: "GM4-",
        direction: "C1",
        ordinal: 13346,
        child: 8,
        probe: "X1+",
        terms: &[
            ("Y1", 150, 1, 1, 1),
            ("A1", 152, 1, 1, 1),
            ("M1", 156, 1, 1, 1),
        ],
    },
    DecompositionCase {
        sg: 225,
        condensate: "GM4-",
        direction: "C2",
        ordinal: 13345,
        child: 8,
        probe: "L1+",
        terms: &[
            ("A1", 152, 1, 1, 1),
            ("V1", 154, 1, 2, 1),
            ("M1", 156, 1, 1, 1),
        ],
    },
    DecompositionCase {
        sg: 225,
        condensate: "GM4-",
        direction: "C2",
        ordinal: 13345,
        child: 8,
        probe: "X1+",
        terms: &[("Y1", 150, 1, 1, 1), ("L1", 155, 1, 2, 1)],
    },
];

fn vector(values: &RationalVector) -> Vec3R {
    Vec3R::new(values.map(|(n, d)| Rat::new(n, d).unwrap()))
}

fn rotation(r: &[i32; 9]) -> [[i32; 3]; 3] {
    [[r[0], r[1], r[2]], [r[3], r[4], r[5]], [r[6], r[7], r[8]]]
}

/// Direct matrix trace at the requested operation, weighting each source arm
/// separately when the translation differs from the archived representative.
fn archived_character(case: &SourceCase, h: &ExactSeitz, lattice: &Lattice) -> Complex64 {
    let mut values = Vec::new();
    for (r, t, trace, blocks) in case.operations {
        if rotation(r) != h.rotation() {
            continue;
        }
        let delta = h.translation().checked_sub(&vector(t)).unwrap();
        if !lattice.contains(&delta).unwrap() {
            continue;
        }
        let raw: Complex64 = blocks.iter().map(|&(re, im)| Complex64::new(re, im)).sum();
        assert!((raw - Complex64::new(trace.0, trace.1)).norm() < 1e-12);
        let value: Complex64 = blocks
            .iter()
            .zip(case.arms)
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
        .expect("mapped operation must occur in the CIR fixture");
    assert!(values.iter().all(|value| (value - first).norm() < 1e-7));
    first
}

#[test]
fn full_star_terms_and_reconstruction_match_independent_cir_witnesses() {
    let mut checked_operations = 0;
    for case in DECOMPOSITIONS {
        let subgroup = isotropy_subgroup_for_direction(
            case.sg,
            case.condensate,
            IsotropyDirection::Label(case.direction),
        )
        .unwrap();
        assert_eq!(subgroup.ordinal, case.ordinal);
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup).unwrap();
        let probe = query::irreps_of(case.sg)
            .iter()
            .find(|row| row.ml == case.probe)
            .unwrap();
        let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
            .unwrap_or_else(|error| panic!("ordinal {}, {}: {error}", case.ordinal, case.probe));
        let source = SOURCE_CASES
            .iter()
            .find(|source| source.sg == case.sg && source.ml == case.probe)
            .unwrap();
        assert_eq!(result.parent_irnumber(), source.irnumber);
        assert_eq!(result.parent_dimension() as usize, source.dimension);
        assert_eq!(result.subgroup_sg(), case.child);
        assert_eq!(result.ordinal(), case.ordinal);
        let mut actual: Vec<Term> = result
            .blocks()
            .iter()
            .flat_map(|block| {
                block.targets().iter().map(|target| {
                    (
                        target.ml,
                        target.irnumber,
                        target.dimension,
                        block.star_size(),
                        target.multiplicity,
                    )
                })
            })
            .collect();
        let mut expected = case.terms.to_vec();
        actual.sort();
        expected.sort();
        assert_eq!(
            actual, expected,
            "ordinal {}, probe {}",
            case.ordinal, case.probe
        );
        let carried: usize = actual
            .iter()
            .map(|&(_, _, dim, size, mult)| usize::from(dim) * size * mult as usize)
            .sum();
        assert_eq!(carried, source.dimension);

        let lattice = Lattice::new(Mat3R::new(
            parent_primitive_basis(case.sg)
                .unwrap()
                .map(|row| row.map(|x| Rat::from_grid(x, 12).unwrap())),
        ))
        .unwrap();
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
            let expected = archived_character(source, h, &lattice);
            assert!(
                (actual - expected).norm() < 1e-7,
                "parent ordinal {}, {}: {h:?}",
                case.ordinal,
                case.probe
            );
            assert!(
                (reconstructed - expected).norm() < 1e-7,
                "child ordinal {}, {}: {h:?}",
                case.ordinal,
                case.probe
            );
            checked_operations += 1;
        }
    }
    assert_eq!(DECOMPOSITIONS.len(), 13);
    // P4/m: 3×8, C2/m: 3×4, P4/nnc: 3×16, Cm: 4×2.
    assert_eq!(checked_operations, 92);
}
