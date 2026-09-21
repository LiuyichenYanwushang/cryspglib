//! Dump one parent irrep's operation-aware character row.
//!
//! Task-9 diagnostic: lets the identity-frequency question ("does this parent
//! irrep contain the subgroup's trivial representation?") be recomputed outside
//! the engine from the shipped row plus the mapped coset representatives
//! printed by `trace_embedding`.
//!
//! ```text
//! cargo run --release -p cryspglib --example dump_character_row -- 221 X3+
//! ```

use cryspglib::irrep::query;
use cryspglib::irrep::types::IrrepRecord;

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let sg: u8 = args
        .next()
        .ok_or("usage: dump_character_row <sg> <ml>")?
        .parse()
        .map_err(|_| "sg must be a number".to_string())?;
    let ml = args.next().ok_or("usage: dump_character_row <sg> <ml>")?;
    let record: &'static IrrepRecord = query::irreps_of(sg)
        .iter()
        .find(|record| record.ml == ml)
        .ok_or_else(|| format!("SG {sg} has no irrep {ml}"))?;
    println!("k = {:?} dim = {}", record.k_vector(), record.dim);
    let row = record
        .ordinary_scalar_selected_arm_block_trace()
        .map_err(|error| format!("row: {error}"))?;
    for (operation, value) in row.operations().iter().zip(row.values()) {
        println!(
            "  rot {:?} trans {:?} -> {:.12} {:+.12}i",
            operation.rotation, operation.translation, value.re, value.im
        );
    }
    Ok(())
}
