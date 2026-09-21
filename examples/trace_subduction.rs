//! Print one record's full-star subduction: blocks, stars, targets, multiplicities.
//!
//! Task-9 diagnostic companion to `trace_embedding`.  Where the audit reports an
//! identity-frequency mismatch, this shows *which* child star the block landed
//! on, how many arms it has, and which targets (with multiplicities) the solver
//! produced -- the information needed to tell a wrong star/arm from a wrong
//! character pairing.
//!
//! ```text
//! cargo run --release -p cryspglib --example trace_subduction -- 12471 X3+
//! ```

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
use cryspglib::irrep::subduction::SubgroupEmbedding;
use cryspglib::irrep::{LabelConvention, types::IrrepRecord};

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
    let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
        .map_err(|error| format!("subduction: {error}"))?;
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
