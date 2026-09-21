//! Print one record's full-star subduction: blocks, stars, targets, multiplicities.
//!
//! Task-9 diagnostic companion to `trace_embedding`.  Where the audit reports an
//! identity-frequency mismatch, this shows *which* child star the block landed
//! on, how many arms it has, and which targets (with multiplicities) the solver
//! produced -- the information needed to tell a wrong star/arm from a wrong
//! character pairing.
//!
//! When the production call fails for missing child data, the tool instead dumps
//! the folded child-star geometry: every folded `q`, its residue modulo the
//! **child** reciprocal lattice, whether it is the child's Gamma point (the only
//! place the subgroup's trivial representation can appear), and the child's own
//! stored `k` points.  That is what separates a genuinely uncovered folded wave
//! vector from a Gamma star that the identity-only entry point can still answer.
//!
//! ```text
//! cargo run --release -p cryspglib --example trace_subduction -- 12471 X3+
//! cargo run --release -p cryspglib --example trace_subduction -- 10027 W1
//! ```

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
use cryspglib::irrep::subduction::star::scalar_star::ScalarStar;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, SubgroupEmbedding};
use cryspglib::irrep::{LabelConvention, types::IrrepRecord};

/// Child reciprocal lattice, built exactly from the child's own primitive basis.
fn child_reciprocal(sg: u8) -> Result<Lattice, String> {
    let basis = cryspglib::irrep::isotropy::parent_primitive_basis(sg)
        .map_err(|error| error.to_string())?;
    let mut rows = [[Rat::ZERO; 3]; 3];
    for (row, values) in basis.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            rows[row][column] = Rat::from_grid(*value, 12).map_err(|error| error.to_string())?;
        }
    }
    let cell = Lattice::new(Mat3R::new(rows)).map_err(|error| error.to_string())?;
    cell.reciprocal().map_err(|error| error.to_string())
}

/// Print the folded child-star geometry of one probe.
///
/// This is the diagnostic path for `MissingIrrepData`: it shows every folded
/// child star, whether each of its points is the child's Gamma point (the only
/// place a trivial representation can appear), and which stored child `k`
/// points exist at all.
fn trace_folded(
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
    child_sg: usize,
) -> Result<(), String> {
    let child_sg = u8::try_from(child_sg).map_err(|_| format!("child {child_sg} is not a space group"))?;
    let star = ScalarStar::new(probe).map_err(|error| error.to_string())?;
    let folded = star
        .folded_stars(embedding)
        .map_err(|error| error.to_string())?;
    let reciprocal = child_reciprocal(child_sg)?;
    println!("folded child stars: {}", folded.len());
    for (index, child_star) in folded.iter().enumerate() {
        println!(
            "  child star {index}: points {} arms {} block_dim {}",
            child_star.star_size(),
            child_star.arm_count(),
            child_star.block_dimension(),
        );
        for point in child_star.points() {
            let reduced = reciprocal
                .reduce(point.q())
                .map_err(|error| error.to_string())?;
            let is_gamma = reciprocal
                .contains(point.q())
                .map_err(|error| error.to_string())?;
            println!(
                "      q = {} arms {:?} | mod child reciprocal = {} | gamma {is_gamma}",
                format_q(point.q()),
                point.arm_indices(),
                format_q(&reduced.representative),
            );
        }
    }
    let mut stored: Vec<String> = query::irreps_of(child_sg)
        .iter()
        .map(|record| format!("{} {}", record.ml, format_k(record.k_vector())))
        .collect();
    stored.sort();
    stored.dedup();
    println!("stored child k points ({}): {}", stored.len(), stored.join(", "));
    Ok(())
}

fn format_q(value: &cryspglib::irrep::subduction::Vec3R) -> String {
    (0..3)
        .map(|axis| {
            let entry = value.get(axis);
            format!("{}/{}", entry.numerator(), entry.denominator())
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn format_k(k: cryspglib::irrep::types::KVector) -> String {
    format!(
        "({},{},{})/{}",
        k.numerators[0], k.numerators[1], k.numerators[2], k.denominator
    )
}

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

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let ordinal: usize = args
        .next()
        .ok_or("usage: trace_subduction <ordinal> [probe-ml]")?
        .parse()
        .map_err(|_| "ordinal must be a number".to_string())?;
    let subgroup = find_subgroup(ordinal)?;
    let probe_ml = args.next().unwrap_or_else(|| subgroup.irrep_ml.to_string());
    let probe: &'static IrrepRecord = query::irreps_of(subgroup.parent_sg)
        .iter()
        .find(|record| record.ml == probe_ml)
        .ok_or_else(|| format!("SG {} has no irrep {probe_ml}", subgroup.parent_sg))?;

    println!(
        "ordinal {ordinal}: SG {} {} -> #{} | probe {} (dim {}, k {:?})",
        subgroup.parent_sg,
        subgroup.irrep_ml,
        subgroup.record.sg,
        probe.ml,
        probe.dim,
        probe.k_vector(),
    );
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
        .map_err(|error| format!("embedding: {error}"))?;
    let result = match subduce_full_star_with_embedding(&subgroup, &embedding, probe) {
        Ok(result) => result,
        Err(error) => {
            println!("subduction failed: {error}");
            trace_folded(&embedding, probe, subgroup.record.sg)?;
            return Err(format!("subduction: {error}"));
        }
    };
    println!(
        "parent dimension {} | blocks {} | representatives {}",
        result.parent_dimension(),
        result.blocks().len(),
        result.representatives().len(),
    );
    for (index, block) in result.blocks().iter().enumerate() {
        println!(
            "  block {index}: q = {:?} stored_k = {:?} star_size {} arm_count {} \
             block_dim {} little_dim {}",
            block.q(),
            block.stored_k(),
            block.star_size(),
            block.arm_count(),
            block.block_dimension(),
            block.little_dimension(),
        );
        for target in block.targets() {
            println!(
                "      target {} (row {}) dim {} multiplicity {} irnumber {:?}",
                target.ml,
                target.row_ml,
                target.dimension,
                target.multiplicity,
                target.irnumber,
            );
        }
    }
    let (parent_characters, reconstructed) = result.reconstruction();
    let total: f64 = parent_characters.iter().map(|value| value.re).sum();
    println!(
        "reconstruction max |difference| = {:.3e}",
        parent_characters
            .iter()
            .zip(reconstructed)
            .map(|(left, right)| (left - right).norm())
            .fold(0.0_f64, f64::max)
    );
    println!(
        "whole-star character sum = {total:.6} over {} representatives -> trivial content {:.6}",
        parent_characters.len(),
        total / parent_characters.len() as f64
    );
    let norm: f64 = parent_characters
        .iter()
        .map(|value| value.norm() * value.norm())
        .sum::<f64>()
        / parent_characters.len() as f64;
    let squared: usize = result
        .blocks()
        .iter()
        .flat_map(|block| block.targets())
        .map(|target| usize::from(target.dimension).pow(2) * target.multiplicity as usize)
        .sum();
    println!("character norm = {norm:.6} | sum of target dim^2 = {squared}");
    println!("characters: {:?}", parent_characters
        .iter()
        .map(|value| (value.re * 1000.0).round() / 1000.0)
        .collect::<Vec<_>>());
    Ok(())
}
