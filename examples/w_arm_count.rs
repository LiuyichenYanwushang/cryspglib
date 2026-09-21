//! Count the star arms of an other-wave-vector *line* that fold onto the child's Gamma point.
//!
//! Task-9 diagnostic for the 5,756 `other_wave_vector_subduction` rows.  Those
//! rows name parent irreps on parameterized lines `k = Gamma + t*v` (decoded
//! from the pinned little table by `scripts/check_other_wave_vector_rows.py`),
//! so their frequency cannot be read off a numeric k.  This tool tests the
//! geometric part of the obvious hypothesis -- that the frequency counts the
//! arms of the line's star whose folded direction lands on the child's Gamma
//! point -- and prints it next to the stored frequencies, so the hypothesis can
//! be checked against the pinned table instead of assumed.
//!
//! ```text
//! cargo run --release -p cryspglib --example w_arm_count -- 13824 DT1 1 0 1 1
//! ```
//!
//! Frames: the little table stores the direction in the parent's **conventional**
//! reciprocal basis, while the embedding folds parent **primitive** fractional
//! coordinates (`k_H = T^T k_G`).  Both readings are printed; which one the
//! engine's own path uses is visible from the arm count that matches the stored
//! frequency.

use std::collections::BTreeSet;

use cryspglib::irrep::isotropy::{self, IsotropySubgroup, parent_primitive_basis};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::LabelConvention;

fn find_subgroup(ordinal: usize) -> Result<IsotropySubgroup, String> {
    for sg in 1u8..=230 {
        for record in query::irreps_of(sg) {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
                .map_err(|error| format!("SG {sg} {}: {error}", record.ml))?;
            if let Some(found) = subgroups.into_iter().find(|s| s.ordinal == ordinal) {
                return Ok(found);
            }
        }
    }
    Err(format!("no ordinary isotropy record has ordinal {ordinal}"))
}

/// Child reciprocal lattice in the child's own fractional coordinates.
fn child_reciprocal(sg: u8) -> Result<Lattice, String> {
    let basis = parent_primitive_basis(sg).map_err(|error| error.to_string())?;
    let mut rows = [[Rat::ZERO; 3]; 3];
    for (row, values) in basis.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            rows[row][column] = Rat::from_grid(*value, 12).map_err(|error| error.to_string())?;
        }
    }
    Lattice::new(Mat3R::new(rows))
        .map_err(|error| error.to_string())?
        .reciprocal()
        .map_err(|error| error.to_string())
}

fn parse_rational(value: &str) -> Result<Rat, String> {
    match value.split_once('/') {
        Some((num, den)) => Ok(Rat::new(
            num.parse().map_err(|_| format!("bad numerator {num}"))?,
            den.parse().map_err(|_| format!("bad denominator {den}"))?,
        )
        .map_err(|error| error.to_string())?),
        None => Ok(Rat::from_integer(
            value.parse().map_err(|_| format!("bad value {value}"))?,
        )),
    }
}

fn format_q(value: &Vec3R) -> String {
    (0..3)
        .map(|axis| {
            let entry = value.get(axis);
            format!("{}/{}", entry.numerator(), entry.denominator())
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let ordinal: usize = args
        .next()
        .ok_or("usage: w_arm_count <ordinal> <label> <vx> <vy> <vz> [<den>]")?
        .parse()
        .map_err(|_| "ordinal must be a number".to_string())?;
    let label = args.next().ok_or("missing w-source label")?;
    let vx = parse_rational(&args.next().ok_or("missing vx")?)?;
    let vy = parse_rational(&args.next().ok_or("missing vy")?)?;
    let vz = parse_rational(&args.next().ok_or("missing vz")?)?;
    let v_conv = Vec3R::new([vx, vy, vz]);

    let subgroup = find_subgroup(ordinal)?;
    let embedding =
        SubgroupEmbedding::from_isotropy_subgroup(&subgroup).map_err(|e| e.to_string())?;
    let parent = subgroup.parent_sg;
    let child = embedding.subgroup_sg();
    let reciprocal = child_reciprocal(child)?;

    // Both candidate parent frames for the decoded direction.
    let primitive = parent_primitive_basis(parent).map_err(|error| error.to_string())?;
    let mut primitive_rows = [[Rat::ZERO; 3]; 3];
    for (row, values) in primitive.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            primitive_rows[row][column] =
                Rat::from_grid(*value, 12).map_err(|error| error.to_string())?;
        }
    }
    let to_primitive = Mat3R::new(primitive_rows);
    let v_primitive = to_primitive
        .checked_mul_vector(&v_conv)
        .map_err(|error| error.to_string())?;

    // Parent point group (unique rotations, column action on direct coordinates).
    let operations = query::symmetry_operations_of(parent).map_err(|error| error.to_string())?;
    let rotations: BTreeSet<[[i32; 3]; 3]> = operations
        .operations
        .iter()
        .map(|operation| operation.rotation)
        .collect();
    let transform = *embedding.transform().matrix();
    let transform_t = transform.transpose();

    println!(
        "ordinal {ordinal}: SG {parent} {label} -> #{child} | v_conv = {} | v_prim = {}",
        format_q(&v_conv),
        format_q(&v_primitive)
    );
    for (name, direction) in [("conv", &v_conv), ("prim", &v_primitive)] {
        let mut arms = 0usize;
        let mut gamma = 0usize;
        let mut first = None;
        for rotation in &rotations {
            // Contragredient action on wave vectors: v' = R^-T v.
            let action = Mat3R::from_ints(*rotation)
                .inverse()
                .map_err(|error| error.to_string())?
                .transpose();
            let arm = action
                .checked_mul_vector(direction)
                .map_err(|error| error.to_string())?;
            let folded = transform_t
                .checked_mul_vector(&arm)
                .map_err(|error| error.to_string())?;
            arms += 1;
            if reciprocal
                .contains(&folded)
                .map_err(|error| error.to_string())?
            {
                gamma += 1;
                if first.is_none() {
                    first = Some(format_q(&arm));
                }
            }
        }
        println!("  frame {name}: point-group arms {arms} | folding to child Gamma {gamma} (first arm {})", first.unwrap_or_else(|| "-".to_string()));
    }
    let stored = isotropy::other_wave_vector_subduction(ordinal).map_err(|e| e.to_string())?;
    println!("  stored rows for this record:");
    for entry in &stored {
        println!(
            "    {} x {} (dim {})",
            entry.frequency, entry.parent_ml, entry.parent_sg
        );
    }
    Ok(())
}
