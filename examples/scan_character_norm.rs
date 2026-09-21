//! Compare the subduced character norm with the reported multiplicities.
//!
//! Task-9 diagnostic.  For irreps of the *subgroup* orthonormality gives
//!
//! ```text
//! (1/|H|) sum_h |chi(h)|^2  =  sum_i mult_i^2
//! ```
//!
//! The engine reports targets as **child tabulated k-points**, which is the
//! convention the stored identity-subduction table also uses: when the subgroup
//! is realized with a larger cell, several child k-points (differing by a vector
//! of the subgroup's own reciprocal lattice) describe the same subgroup irrep
//! and are reported separately.  The identity therefore does not have to hold as
//! written; a mismatch marks records where the unfolded child-k description is
//! coarser than the subgroup's own irrep labelling, which is where a naive
//! "does this irrep contain the trivial representation" test goes wrong.
//! Fixing the mismatch means folding the child k-points modulo the subgroup
//! lattice, not changing the character evaluation.
//!
//! ```text
//! cargo run --release -p cryspglib --example scan_character_norm
//! ```

use cryspglib::irrep::isotropy;
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::SubgroupEmbedding;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
use cryspglib::irrep::{LabelConvention, types::IrrepRecord};

fn main() -> Result<(), String> {
    let tolerance = 1e-6;
    let mut checked = 0usize;
    let mut violations = Vec::new();
    let mut errors = 0usize;
    for sg in 1u8..=230 {
        for record in query::irreps_of(sg) {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
                .map_err(|error| format!("SG {sg} {}: {error}", record.ml))?;
            // Only the clean case: an ordinary parent row.  Compound and
            // realification rows carry physical sums whose "norm" is not a
            // representation norm.
            if record.compound_metadata().is_some() {
                continue;
            }
            let probe: &'static IrrepRecord = record;
            for subgroup in subgroups {
                let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&subgroup) else {
                    errors += 1;
                    continue;
                };
                let Ok(result) =
                    subduce_full_star_with_embedding(&subgroup, &embedding, probe)
                else {
                    errors += 1;
                    continue;
                };
                let (characters, _reconstructed) = result.reconstruction();
                if characters.is_empty() {
                    continue;
                }
                let norm: f64 = characters
                    .iter()
                    .map(|value| value.norm() * value.norm())
                    .sum::<f64>()
                    / characters.len() as f64;
                // A compound row is the sum of two complex constituents: the
                // reported `dimension` is the row's, so its squared contribution
                // has to be folded by the constituent count exactly like the
                // Gamma Frobenius sum.
                let squared: usize = result
                    .blocks()
                    .iter()
                    .flat_map(|block| block.targets())
                    .map(|target| {
                        let row = query::irreps_of(target.sg)
                            .iter()
                            .find(|record| record.ml == target.row_ml);
                        if row.map(|record| record.compound_metadata().is_some())
                            .unwrap_or(true)
                        {
                            return 0;
                        }
                        // Orthonormality of irreps: <chi, chi> = sum of squared
                        // multiplicities, not sum of dim^2.
                        (target.multiplicity as usize).pow(2)
                    })
                    .sum();
                checked += 1;
                if (norm - squared as f64).abs() > tolerance {
                    violations.push((subgroup.ordinal, sg, probe.ml, subgroup.record.sg, norm, squared));
                }
            }
        }
    }
    println!(
        "checked {checked} (record, condensing irrep) pairs; \
         <chi,chi> violations {}; call errors {errors}",
        violations.len()
    );
    for (ordinal, sg, ml, child, norm, squared) in violations.iter().take(20) {
        println!(
            "   ordinal {ordinal}: SG {sg} {ml} -> #{child}  <chi,chi> {norm:.3} != sum mult^2 {squared}"
        );
    }
    Ok(())
}
