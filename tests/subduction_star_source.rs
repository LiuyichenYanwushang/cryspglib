//! Independent full-star witnesses: archived CIR matrices, not induced traces.

use cryspglib::irrep::isotropy::parent_primitive_basis;
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::OrdinaryStar;
use cryspglib::irrep::subduction::{ExactSeitz, Lattice, Mat3R, Rat, Vec3R, strict_sg_hall_ops};
use cryspglib::irrep::types::IrrepSourceIdentity;
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

fn vector(values: &RationalVector) -> Vec3R {
    Vec3R::new(values.map(|(n, d)| Rat::new(n, d).unwrap()))
}

fn lattice(sg: u8) -> Lattice {
    let basis = parent_primitive_basis(sg).unwrap();
    Lattice::new(Mat3R::new(
        basis.map(|row| row.map(|value| Rat::from_grid(value, 12).unwrap())),
    ))
    .unwrap()
}

#[test]
fn induced_characters_match_archived_full_matrices_and_translated_blocks() {
    let mut checked = 0;
    for fixture in SOURCE_CASES {
        let probe = query::irreps_of(fixture.sg)
            .iter()
            .find(|record| record.ml == fixture.ml)
            .unwrap();
        assert_eq!(
            probe.source_identity(),
            IrrepSourceIdentity::OrdinaryScalar {
                cir_irnumber: fixture.irnumber,
            }
        );
        let star = OrdinaryStar::new(probe).unwrap();
        assert_eq!(star.dimension() as usize, fixture.dimension);
        assert_eq!(star.arms().len(), fixture.arms.len());
        let cell = lattice(fixture.sg);
        let reciprocal = cell.reciprocal().unwrap();
        // These fixtures use the same axes/origin in the source and data Hall
        // frames. Check the actual operations and arm sets before comparing
        // characters: source array positions alone are not a frame transform.
        for source_arm in fixture.arms {
            assert!(star.arms().iter().any(|arm| {
                reciprocal
                    .same_mod(arm.wave_vector(), &vector(source_arm))
                    .unwrap()
            }));
        }
        let hall = strict_sg_hall_ops(fixture.sg).unwrap();
        let translations = [
            Vec3R::zero(),
            Vec3R::from_ints([2, -1, 3]),
            Vec3R::new(*cell.rows().row(0)),
            Vec3R::new(*cell.rows().row(1)),
        ];
        for (rotation, translation, trace, blocks) in fixture.operations {
            let r = [
                [rotation[0], rotation[1], rotation[2]],
                [rotation[3], rotation[4], rotation[5]],
                [rotation[6], rotation[7], rotation[8]],
            ];
            let source_op = ExactSeitz::new(r, vector(translation));
            assert!(hall.operations.iter().any(|operation| {
                operation.rotation() == r
                    && cell
                        .same_mod(operation.translation(), source_op.translation())
                        .unwrap()
            }));
            assert_eq!(blocks.len(), fixture.arms.len());
            let raw_trace: Complex64 = blocks.iter().map(|&(re, im)| Complex64::new(re, im)).sum();
            assert!((raw_trace - Complex64::new(trace.0, trace.1)).norm() < 1e-12);
            for shift in &translations {
                assert!(cell.contains(shift).unwrap());
                // Left multiplying by {E|L} weights EACH archived diagonal
                // arm block by its own translation eigenvalue. A full-star
                // trace has no single Bloch phase that can be factored out.
                let expected: Complex64 = fixture
                    .arms
                    .iter()
                    .zip(*blocks)
                    .map(|(k, &(re, im))| {
                        let angle: f64 = k
                            .iter()
                            .enumerate()
                            .map(|(i, &(n, d))| {
                                Rat::new(n, d)
                                    .unwrap()
                                    .checked_mul(shift.get(i))
                                    .unwrap()
                                    .to_f64()
                            })
                            .sum();
                        Complex64::new(re, im)
                            * Complex64::from_polar(1.0, std::f64::consts::TAU * angle)
                    })
                    .sum();
                let operation =
                    ExactSeitz::new(r, source_op.translation().checked_add(shift).unwrap());
                let actual = star.character(&operation).unwrap();
                assert!(
                    (actual - expected).norm() < 1e-7,
                    "SG{} {} source#{} op={operation:?}: {actual} != {expected}",
                    fixture.sg,
                    fixture.ml,
                    fixture.irnumber
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 1232);
}

#[test]
fn every_ordinary_record_has_its_source_dimension_and_complete_star() {
    let mut checked = 0;
    for sg in 1..=230 {
        for probe in query::irreps_of(sg) {
            if !matches!(
                probe.source_identity(),
                IrrepSourceIdentity::OrdinaryScalar { .. }
            ) {
                continue;
            }
            let star = OrdinaryStar::new(probe)
                .unwrap_or_else(|error| panic!("SG{sg} {}: {error}", probe.ml));
            assert_eq!(star.dimension(), u32::from(probe.dim));
            assert_eq!(
                star.dimension() as usize,
                probe
                    .ordinary_scalar_selected_arm_block_trace()
                    .unwrap()
                    .dimension()
                    * star.arms().len()
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 4105);
}

#[test]
fn changing_transporter_translations_preserves_the_induced_character() {
    for (sg, ml) in [(92, "X1"), (139, "P1"), (198, "X1")] {
        let probe = query::irreps_of(sg)
            .iter()
            .find(|record| record.ml == ml)
            .unwrap();
        let star = OrdinaryStar::new(probe).unwrap();
        let cell = lattice(sg);
        let mut transporters = Vec::new();
        for (index, arm) in star.arms().iter().enumerate() {
            let transporter = arm.transporter();
            transporters.push(if index == 0 {
                // The public transversal contract pins the seed to {E|0}.
                *transporter
            } else {
                ExactSeitz::new(
                    transporter.rotation(),
                    transporter
                        .translation()
                        .checked_add(&Vec3R::new(*cell.rows().row(index % 3)))
                        .unwrap(),
                )
            });
        }
        transporters.reverse();
        let rephased = OrdinaryStar::from_transporters(probe, &transporters).unwrap();
        assert_eq!(star.dimension(), rephased.dimension());
        for operation in strict_sg_hall_ops(sg).unwrap().operations {
            let expected = star.character(&operation).unwrap();
            let actual = rephased.character(&operation).unwrap();
            assert!(
                (actual - expected).norm() < 1e-7,
                "SG{sg} {ml}: {operation:?}"
            );
        }
    }
}
