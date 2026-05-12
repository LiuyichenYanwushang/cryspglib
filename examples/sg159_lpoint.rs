//! SG159.62 (P31c1', UNI1292) L-point spinor irreps example.
//!
//! Demonstrates: magnetic little group, double-group irreps, Wigner classification.

use cryspglib::irrep::query::*;
use cryspglib::irrep::corep::*;
use cryspglib::irrep::wigner::*;
use cryspglib::irrep::types::IrrepRecord;
use cryspglib::MagneticSpaceGroupType;

fn main() {
    let uni = 1293;
    let msg = MagneticSpaceGroupType::from_uni(uni);
    println!("=== UNI{}: {} ({}), type={:?} ===", uni, msg.bns_number.trim(), msg.og_number.trim(), msg.type_);

    // ═══ Magnetic space group operations ═══
    let mag_ops = cryspglib::SymmetryOps::from_magnetic_database(uni).unwrap();
    let n_u = mag_ops.operations.iter().filter(|op| !op.time_reversal).count();
    let n_a = mag_ops.operations.iter().filter(|op| op.time_reversal).count();
    println!("\nFull magnetic group: {} ops ({} unitary + {} antiunitary)", mag_ops.len(), n_u, n_a);

    // ═══ H = unitary subgroup ═══
    let h_info = identify_unitary_subgroup_with_hall(uni).unwrap();
    println!("Unitary subgroup H: SG{} (Hall {})", h_info.sg, h_info.hall);

    // ═══ L-point (kx=1, ky=0, kz=1, kd=2) ═══
    let kx: i8 = 1; let ky: i8 = 0; let kz: i8 = 1; let kd: i8 = 2;
    println!("\n=== L-point k=({},{},{})/{} ===", kx, ky, kz, kd);

    // Filter little group
    let mag_seitz = ops_to_seitz(&mag_ops);
    let lg_indices = filter_little_group(kx, ky, kz, kd, &mag_ops);
    let unitary_lg: Vec<usize> = lg_indices.iter().filter(|&&i| !mag_ops[i].time_reversal).copied().collect();
    let antiunitary_lg: Vec<usize> = lg_indices.iter().filter(|&&i| mag_ops[i].time_reversal).copied().collect();
    println!("Magnetic little group at L: {} ops ({} unitary + {} antiunitary)",
        lg_indices.len(), unitary_lg.len(), antiunitary_lg.len());

    // Show all little group operations
    println!("\n  idx | rot                          | trans              | θ");
    println!("  ----+------------------------------+---------------------+---");
    for &i in &lg_indices {
        let op = &mag_ops[i];
        let r = op.rotation;
        println!("  {:3} | [{:2},{:2},{:2};{:2},{:2},{:2};{:2},{:2},{:2}] | ({:4.1},{:4.1},{:4.1}) | {}",
            i, r[0][0],r[0][1],r[0][2], r[1][0],r[1][1],r[1][2], r[2][0],r[2][1],r[2][2],
            op.translation[0], op.translation[1], op.translation[2],
            if op.time_reversal { "θ" } else { " " });
    }

    // ═══ Irreps at L-point ═══
    let h_ops = h_info.ops_from_msg;
    let h_seitz = ops_to_seitz(&h_ops);
    let h_sg = h_info.sg as u8;
    let l_irreps: Vec<&IrrepRecord> = irreps_of(h_sg).iter()
        .filter(|r| r.kx == kx && r.ky == ky && r.kz == kz && r.kd == kd)
        .collect();

    println!("\n=== {} irreps at L-point ({} scalar + {} spinor) ===",
        l_irreps.len(),
        l_irreps.iter().filter(|r| !r.spinor).count(),
        l_irreps.iter().filter(|r| r.spinor).count());

    for ir in &l_irreps {
        if ir.spinor {
            // Spinor: show Wigner classification
            let h_spin = ir.spin_ops();
            let g_sg = parent_spatial_sg(uni).unwrap_or(h_sg as usize) as u8;
            let g_spin = if g_sg == h_sg { h_spin } else { IrrepRecord::spin_ops_for_sg(g_sg) };
            let ctx = SpinLiftContext { h: h_spin, g: g_spin, sg: h_sg };
            println!("  g_sg={} H_spin_count={} G_spin_count={}", g_sg,
                h_spin.0.len()/9, g_spin.0.len()/9);
            let n_lg = ir.spin_lg_char_count();
            let indices = ir.spin_lg_op_indices();
            let a0_idx = antiunitary_lg.iter().copied()
                .find(|&i| mag_ops[i].rotation == [[1,0,0],[0,1,0],[0,0,1]])
                .unwrap_or(antiunitary_lg[0]);

            let ct = wigner_classify_spinor(
                &ctx, ir.characters(), n_lg, indices,
                &unitary_lg, &mag_seitz, &h_seitz, a0_idx,
                kx, ky, kz, kd,
            );

            // Debug: trace why None
            if ct.is_none() {
                let (g_spin_rots, g_spin_trans, g_spin_su2) = ctx.g;
                let g_spin_seitz = build_spin_seitz(g_spin_rots, g_spin_trans);
                let a0 = &mag_seitz[a0_idx];
                let a0_match = g_spin_seitz.iter().position(|s| s.rot == a0.rot);
                println!("      a0_rot={:?} a0_trans=({:.3},{:.3},{:.3})",
                    a0.rot, a0.trans[0], a0.trans[1], a0.trans[2]);
                println!("      a0_in_G_spin={:?}", a0_match);
                println!("      G spin ops rots: {:?}",
                    g_spin_seitz.iter().map(|s| s.rot).collect::<Vec<_>>());
                if let Some(am) = a0_match {
                    let u_a0 = spin_su2_at(g_spin_su2, am);
                    println!("      u_a0={:?}", u_a0);
                    // Try per-term: does (a0*h)^2 match back?
                    let (_, origin) = IrrepRecord::sg_setting(ctx.sg);
                    let a0_bilbao = SeitzOp::new(a0.rot,
                        [a0.trans[0], a0.trans[1], a0.trans[2]], false);
                    let h_spin_seitz = build_spin_seitz(ctx.h.0, ctx.h.1);
                    for local in 0..n_lg {
                        let gsi = indices[local] as usize;
                        let h_spin = &h_spin_seitz[gsi];
                        let (g0h, _l1) = compose_seitz(&a0_bilbao, h_spin);
                        let (sq, _lsq) = square_seitz(&g0h);
                        let sq_match = h_spin_seitz.iter().position(|s| s.rot == sq.rot);
                        println!("      [{}] h_rot={:?} sq_rot={:?} sq_in_H_spin={:?}",
                            local, h_spin.rot, sq.rot, sq_match);
                    }
                }
            }

            println!("  {:6} spinor dim={} spin_lg={:?} n_lg={} spin_chars={:?} Wigner={:?}",
                ir.ml, ir.dim, indices, n_lg,
                &ir.characters()[..n_lg.min(6)],
                ct);
        } else {
            // Scalar: show character table
            let chars = ir.characters();
            println!("  {:6} scalar dim={} image={} chars={:?}",
                ir.ml, ir.dim, ir.image,
                &chars[..chars.len().min(8)]);
        }
    }

    // ═══ Scalar irreps via corepresentation() ═══
    println!("\n=== Scalar corepresentations ===");
    for ir in &l_irreps {
        if ir.spinor { continue; }
        if let Some(corep) = ir.corepresentation(uni) {
            println!("  {}: type={:?} dim={} χ(id)={:.1}",
                ir.ml, corep.corep_type, corep.dim, corep.characters[0]);
        }
    }

    // ═══ Summary ═══
    println!("\n=== Summary ===");
    println!("  Magnetic group: {}", msg.og_number.trim());
    println!("  Unitary subgroup: SG{}", h_sg);
    println!("  L-point: ({},{},{})/{}", kx, ky, kz, kd);
    println!("  Little group: {} ops", lg_indices.len());
    println!("  Irreps at L: {} (scalar) + {} (spinor)",
        l_irreps.iter().filter(|r| !r.spinor).count(),
        l_irreps.iter().filter(|r| r.spinor).count());
}
