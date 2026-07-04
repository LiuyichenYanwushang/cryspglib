//! # Magnetic irrep summary example
//!
//! Demonstrates querying the complete magnetic irrep summary for a magnetic
//! space group by BNS label or UNI number.

use cryspglib::irrep::magnetic_summary::*;

fn main() {
    // ── BNS 128.406 (black-white, Type III) ──
    println!("=== BNS 128.406 ===");
    let summary = magnetic_irrep_summary_by_bns("128.406").unwrap();
    println!("{}", format_magnetic_irrep_summary(&summary));

    // Show details for Z point
    if let Some(z) = summary.kpoints.iter().find(|k| k.label == "Z") {
        println!();
        println!("--- Z point details ---");
        for c in &z.coreps {
            println!(
                "  {}  type={:?}  dim={}  source={:?}",
                c.label, c.corep_type, c.dim, c.source
            );
            for ic in &c.isotropy_candidates {
                println!(
                    "    isotropy from {} ({:?}): {} ordinary + {} magnetic",
                    ic.source_ml,
                    ic.relation,
                    ic.ordinary.len(),
                    ic.magnetic.len()
                );
            }
        }
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
