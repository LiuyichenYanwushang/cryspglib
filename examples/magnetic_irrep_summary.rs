//! # Magnetic irrep summary example
//!
//! Demonstrates querying the little-group corep summary for a magnetic
//! space group by BNS label or UNI number.

use cryspglib::irrep::magnetic_summary::*;

fn main() {
    // ── BNS 128.406 (currently exposes an unsupported scalar path) ──
    println!("=== BNS 128.406 ===");
    match magnetic_irrep_summary_by_bns("128.406") {
        Ok(summary) => println!("{}", format_magnetic_irrep_summary(&summary)),
        Err(err) => println!("corep summary failed: {:?}", err),
    }

    // ── UNI 2 (grey P1, Type II) ──
    println!();
    println!("=== UNI 2 (grey P1) ===");
    let summary = magnetic_irrep_summary_by_uni(2).unwrap();
    for kp in &summary.kpoints {
        println!(
            "{} {:?} |LG|={} U={} A={}  coreps={}",
            kp.label,
            kp.coords,
            kp.little_group_order,
            kp.unitary_order,
            kp.antiunitary_order,
            kp.coreps.len()
        );
    }
}
