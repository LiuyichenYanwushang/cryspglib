//! Read-only dump of the selected-arm little-group block trace of one irrep.
//!
//! The other-wave-vector rows name irreps at parametric k domains whose
//! characters the pinned irrep table does not carry.  For an irrep at a k point
//! that lies on such a domain, the compatibility data already says which line
//! irrep it restricts to (`X3+: DT3` for SG 202), and the block trace printed
//! here *is* that little representation on the little group's operations.
//! Dividing out the Bloch phase at the point gives the line irrep's Gamma
//! characters, which is the route `docs/task9-remaining-work.md` records for the
//! eight sources the Gamma compatibility system cannot separate.
//!
//! Usage: `cargo run --release -p cryspglib --example dump_little_trace -- 202 X3+`

use cryspglib::irrep::query;

fn main() {
    let mut args = std::env::args().skip(1);
    let sg: u8 = args
        .next()
        .expect("usage: dump_little_trace <space group> <irrep label>")
        .parse()
        .expect("space group number");
    let label = args
        .next()
        .expect("usage: dump_little_trace <space group> <irrep label>");

    let mut found = false;
    for record in query::irreps_of(sg) {
        if record.ml != label {
            continue;
        }
        found = true;
        match record.ordinary_scalar_selected_arm_block_trace() {
            Ok(row) => {
                println!(
                    "SG {sg} {label}: little dim {} over {} operations",
                    row.dimension(),
                    row.len()
                );
                for index in 0..row.len() {
                    if let Some((value, operation)) = row.entry(index) {
                        println!(
                            "  {operation:?} char {:.6} {:+.6}i",
                            value.re, value.im
                        );
                    }
                }
            }
            Err(error) => println!("SG {sg} {label}: no scalar block trace ({error:?})"),
        }
    }
    if !found {
        eprintln!("SG {sg}: no irrep labelled {label}");
        std::process::exit(2);
    }
}
