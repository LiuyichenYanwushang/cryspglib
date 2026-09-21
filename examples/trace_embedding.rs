//! Print the engine's mapped operations for one isotropy record.
//!
//! Task-9 diagnostic: the production audit reports an identity-frequency
//! mismatch for a record, and the first question is whether the *embedding* puts
//! the subgroup where the official program puts it.  This prints the transform,
//! the subgroup lattice and every mapped coset representative, so the answer can
//! be compared with the official `SHOW ELEMENTS` column instead of guessed.
//!
//! ```text
//! cargo run --release -p cryspglib --example trace_embedding -- 12471
//! ```

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::{LabelConvention, subduction::SubgroupEmbedding};

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
    let ordinal: usize = std::env::args()
        .nth(1)
        .ok_or("usage: trace_embedding <ordinal>")?
        .parse()
        .map_err(|_| "ordinal must be a number".to_string())?;
    let subgroup = find_subgroup(ordinal)?;
    println!(
        "ordinal {ordinal}: SG {} {} {} -> #{} size {}",
        subgroup.parent_sg,
        subgroup.irrep_ml,
        subgroup.record.direction,
        subgroup.record.sg,
        isotropy::subgroup_size(subgroup.record.basis).map_err(|e| e.to_string())?,
    );
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)
        .map_err(|error| format!("embedding: {error}"))?;
    println!(
        "setting {:?} / {}  candidates {}  representatives {}",
        embedding.setting(),
        embedding.setting_denominator(),
        embedding.candidate_count(),
        embedding.representatives().len(),
    );
    println!("child shift {:?}", embedding.child_shift());
    let transform = embedding.transform();
    println!("transform matrix {:?}", transform.matrix());
    println!("transform origin {:?}", transform.origin());
    println!("subgroup lattice {:?}", embedding.subgroup_lattice().rows());
    println!("parent lattice   {:?}", embedding.parent_lattice().rows());
    println!("--- mapped representatives (rotation rows, translation) ---");
    for operation in embedding.representatives() {
        let rotation = operation.rotation();
        let translation = operation.translation();
        let text: Vec<String> = rotation
            .iter()
            .map(|row| format!("[{}, {}, {}]", row[0], row[1], row[2]))
            .collect();
        let t: Vec<String> = translation
            .as_array()
            .iter()
            .map(|value| format!("{value:?}"))
            .collect();
        println!("  {} | t = {}", text.join(" "), t.join(", "));
    }
    Ok(())
}
