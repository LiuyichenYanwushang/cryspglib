//! Dump every non-spinor character row of one space group, for offline checks.
//!
//! ```text
//! cargo run --release -p cryspglib --example dump_sg_rows -- 221
//! ```

use cryspglib::irrep::query;

fn main() -> Result<(), String> {
    let sg: u8 = std::env::args()
        .nth(1)
        .ok_or("usage: dump_sg_rows <sg>")?
        .parse()
        .map_err(|_| "sg must be a number".to_string())?;
    for record in query::irreps_of(sg) {
        if record.spinor {
            continue;
        }
        let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() else {
            continue;
        };
        let k = record.k_vector();
        println!(
            "IRREP {} k {} {} {} dim {} entries {}",
            record.ml,
            k.numerators[0],
            k.numerators[1],
            k.numerators[2],
            record.dim,
            row.values().len()
        );
        for (operation, value) in row.operations().iter().zip(row.values()) {
            println!(
                "  {} {} {} {} {} {} {} {} {} | {} {} {} | {:.9} {:.9}",
                operation.rotation[0],
                operation.rotation[1],
                operation.rotation[2],
                operation.rotation[3],
                operation.rotation[4],
                operation.rotation[5],
                operation.rotation[6],
                operation.rotation[7],
                operation.rotation[8],
                operation.translation[0],
                operation.translation[1],
                operation.translation[2],
                value.re,
                value.im,
            );
        }
    }
    Ok(())
}
