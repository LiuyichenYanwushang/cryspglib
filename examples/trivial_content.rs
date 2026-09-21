//! Recompute a parent irrep's trivial content on a subgroup by direct lookup.
//!
//! Task-9 diagnostic.  The production path evaluates the *induced full-star*
//! character arm by arm (`arm_character` returns zero for an operation that
//! moves an arm).  This example instead evaluates the parent's operation-aware
//! row at each mapped coset representative, exactly as `character_of` pairs a
//! row entry (rotation, translation modulo the parent lattice, Bloch phase),
//! and averages over the representatives:
//!
//! ```text
//! mult(trivial_H, D|H) = (1/|H|) sum_h chi_D(h)
//! ```
//!
//! Comparing the two numbers tells a wrong star/arm/fixity decision from a wrong
//! character pairing.
//!
//! ```text
//! cargo run --release -p cryspglib --example trivial_content -- 12471 X3+
//! ```

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Rat, Vec3R, bloch_phase};
use cryspglib::irrep::subduction::SubgroupEmbedding;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
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
        .ok_or("usage: trivial_content <ordinal> [probe-ml]")?
        .parse()
        .map_err(|_| "ordinal must be a number".to_string())?;
    let subgroup = find_subgroup(ordinal)?;
    let probe_ml = args.next().unwrap_or_else(|| subgroup.irrep_ml.to_string());
    let probe: &'static IrrepRecord = query::irreps_of(subgroup.parent_sg)
        .iter()
        .find(|record| record.ml == probe_ml)
        .ok_or_else(|| format!("SG {} has no irrep {probe_ml}", subgroup.parent_sg))?;

    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
        .map_err(|error| format!("embedding: {error}"))?;
    let lattice = *embedding.parent_lattice();
    let wave_vector = Vec3R::new([
        Rat::new(i128::from(probe.k_vector().numerators[0]), i128::from(probe.k_vector().denominator))
            .map_err(|e| e.to_string())?,
        Rat::new(i128::from(probe.k_vector().numerators[1]), i128::from(probe.k_vector().denominator))
            .map_err(|e| e.to_string())?,
        Rat::new(i128::from(probe.k_vector().numerators[2]), i128::from(probe.k_vector().denominator))
            .map_err(|e| e.to_string())?,
    ]);
    let row = probe
        .ordinary_scalar_selected_arm_block_trace()
        .map_err(|error| format!("row: {error}"))?;

    let mut total = num_complex::Complex64::new(0.0, 0.0);
    let mut matched = 0usize;
    for operation in embedding.representatives() {
        let rotation = operation.rotation();
        for (value, candidate) in row.values().iter().zip(row.operations()) {
            let mut candidate_rotation = [[0i32; 3]; 3];
            for (index, slot) in candidate_rotation.iter_mut().enumerate() {
                *slot = [
                    candidate.rotation[index * 3],
                    candidate.rotation[index * 3 + 1],
                    candidate.rotation[index * 3 + 2],
                ];
            }
            if candidate_rotation != rotation {
                continue;
            }
            // The row's translation is on the source grid; the engine divides by
            // it when it builds the row, so mirror that here.
            let grid = 12.0_f64;
            let candidate_translation = Vec3R::new([
                Rat::new((candidate.translation[0] * grid).round() as i128, 12)
                    .map_err(|e| e.to_string())?,
                Rat::new((candidate.translation[1] * grid).round() as i128, 12)
                    .map_err(|e| e.to_string())?,
                Rat::new((candidate.translation[2] * grid).round() as i128, 12)
                    .map_err(|e| e.to_string())?,
            ]);
            let delta = operation
                .translation()
                .checked_sub(&candidate_translation)
                .map_err(|e| e.to_string())?;
            if !lattice.contains(&delta).map_err(|e| e.to_string())? {
                continue;
            }
            total += value * bloch_phase(&wave_vector, &delta).map_err(|e| e.to_string())?;
            matched += 1;
            break;
        }
    }
    let count = embedding.representatives().len();
    println!(
        "ordinal {ordinal} probe {probe_ml}: matched {matched} of {count} representatives, \
         sum = {:.6}, sum/|H| = {:.6}",
        total.re,
        total.re / count as f64
    );
    let stored = subgroup
        .identity_subduction()
        .map_err(|error| format!("stored: {error}"))?;
    println!("stored identity subduction ({} entries):", stored.len());
    for entry in &stored {
        println!(
            "   {} k = {:?} frequency {} domain {} direction {}",
            entry.parent_ml,
            entry.parent_k,
            entry.frequency,
            entry.domain,
            entry.direction_label,
        );
    }
    let result = subduce_full_star_with_embedding(&subgroup, &embedding, probe)
        .map_err(|error| format!("subduction: {error}"))?;
    let trivial = result
        .blocks()
        .iter()
        .map(|block| block.multiplicity("GM1+"))
        .sum::<u32>();
    println!("production full-star trivial multiplicity: {trivial}");
    Ok(())
}
