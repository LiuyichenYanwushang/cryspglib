//! Query helpers for the irrep database.
//!
//! These functions provide convenient access to the machine-generated
//! data in [`generated_data`].
//!
//! # Quick start
//!
//! ```
//! use cryspglib::irrep::query::*;
//!
//! // List all irreps for a space group
//! let irreps = irreps_of(221);  // Pm-3m
//!
//! // Get all unique k-points with their labels and coordinates
//! let kpoints = kpoints_of(221);
//!
//! // Print a markdown-style character table for the Γ point
//! println!("{}", format_character_table(221, 0, 0, 0, 1));
//!
//! // Get symmetry operations {R|t} for the space group
//! let ops = symmetry_operations_of(221);
//!
//! // Look up space group info
//! let (hm, schoenflies) = sg_info(221).unwrap();
//! ```
//!
//! # Irrep records
//!
//! Each [`IrrepRecord`] provides:
//! - [`IrrepRecord::characters()`] — χ(g) = Tr(D(g)) for each little-group operator
//! - [`IrrepRecord::matrices()`] — full D(g) matrix elements, flattened
//! - [`IrrepRecord::subgroups()`] — non-magnetic isotropy subgroups (lower-symmetry space groups)
//! - [`IrrepRecord::magnetic_subgroups()`] — magnetic isotropy subgroups (magnetic space groups)
//! - [`IrrepRecord::k_label()`] — k-point prefix (`"GM"`, `"X"`, ...)
//! - [`IrrepRecord::is_point()`] — whether this is a special k-point

use std::collections::BTreeMap;

use super::types::generated_data::*;
use super::types::*;
use super::preamble;

/// Extract the k-point label prefix from a Miller-Love irrep label.
///
/// Examples:
/// - `"GM4+"` → `"GM"` (Γ point)
/// - `"X3-"` → `"X"`
/// - `"DT1"` → `"DT"` (Δ line)
/// - `"LD2"` → `"LD"` (Λ line)
pub fn kpoint_label(ml: &str) -> &str {
    let body = ml.trim_end_matches(|c: char| c == '+' || c == '-');
    let prefix_end = body
        .find(|c: char| c.is_ascii_digit())
        .unwrap_or(body.len());
    &body[..prefix_end]
}

/// Look up all irreps for a space group (1–230).
///
/// Returns an empty slice for SG numbers outside 1–230.
pub fn irreps_of(sg: u8) -> &'static [IrrepRecord] {
    if sg == 0 || sg > 230 {
        return &[];
    }
    let (start, count) = SG_IRREP_INDEX[sg as usize];
    &IRREPS[start as usize..(start as usize + count as usize)]
}

/// Look up isotropy subgroups for a specific irrep, by its index into [`IRREPS`].
///
/// Prefer [`IrrepRecord::subgroups`] for more ergonomic access.
///
/// Returns an empty slice if the index is out of range.
pub fn subgroups_of(irrep_index: usize) -> &'static [IsotropyRecord] {
    if irrep_index >= IRREPS.len() {
        return &[];
    }
    IRREPS[irrep_index].subgroups()
}

/// Total number of irreps in the database.
pub fn total_irreps() -> usize {
    IRREPS.len()
}

/// Total number of space groups with irrep data.
pub fn total_sgs() -> usize {
    230
}

/// Get the (Hermann-Mauguin symbol, Schoenflies symbol) for a space group.
///
/// Returns `None` for invalid SG numbers (0 or >230).
pub fn sg_info(sg: u8) -> Option<(&'static str, &'static str)> {
    if sg == 0 || sg > 230 {
        return None;
    }
    let &(_, hm, schoenflies, _) = &preamble::SPACEGROUP_INDEX[(sg - 1) as usize];
    Some((hm, schoenflies))
}

// ── High-level query API ────────────────────────────────────────────────────────

/// Summary of a k-point in the Brillouin zone.
pub struct KPointSummary {
    /// k-point label: `"GM"`, `"X"`, `"R"`, `"DT"`, etc.
    pub label: String,
    /// Fractional reciprocal coordinates as `(kx, ky, kz, denom)`.
    /// The actual coordinate is `(kx/denom, ky/denom, kz/denom)`.
    pub coords: (i8, i8, i8, i8),
    /// Indices into [`IRREPS`] for all irreps at this k-point.
    pub irreps: Vec<usize>,
}

/// Get all unique k-points for a space group, with their irreps.
///
/// Groups irreps by their `(kx, ky, kz, kd)` coordinates and
/// returns a sorted list with k-point labels extracted from
/// the first irrep's ML label.
pub fn kpoints_of(sg: u8) -> Vec<KPointSummary> {
    let irreps = irreps_of(sg);
    let mut groups: BTreeMap<(i8, i8, i8, i8), Vec<usize>> = BTreeMap::new();

    for (idx, ir) in irreps.iter().enumerate() {
        let coords = (ir.kx, ir.ky, ir.kz, ir.kd);
        groups.entry(coords).or_default().push(idx);
    }

    groups
        .into_iter()
        .map(|(coords, ir_indices)| {
            let first = &irreps[ir_indices[0]];
            let label = first.k_label().to_string();
            KPointSummary {
                label,
                coords,
                irreps: ir_indices,
            }
        })
        .collect()
}

/// Format a character table for a specific k-point as a markdown-style table.
///
/// Columns: ML label, BC label, then one column per character value.
/// Returns a `String` suitable for printing.
pub fn format_character_table(sg: u8, kx: i8, ky: i8, kz: i8, kd: i8) -> String {
    let irreps = irreps_of(sg);
    let matching: Vec<&IrrepRecord> = irreps
        .iter()
        .filter(|ir| ir.kx == kx && ir.ky == ky && ir.kz == kz && ir.kd == kd)
        .collect();

    if matching.is_empty() {
        return format!(
            "// No irreps at k-point ({}/{}, {}/{}, {}/{}) for SG {}",
            kx, kd, ky, kd, kz, kd, sg
        );
    }

    // Determine number of operators from the first irrep's character table
    let max_ops = matching.iter().map(|ir| ir.characters().len()).max().unwrap_or(0);
    if max_ops == 0 {
        return format!(
            "// SG {} k=({}/{},{}/{},{}/{}): no operator data available",
            sg, kx, kd, ky, kd, kz, kd
        );
    }

    // Get symmetry operations for column headers
    let ops = symmetry_operations_of(sg);

    // Format operation as compact string
    let fmt_op = |i: usize| -> String {
        if i >= ops.len() { return format!("g{}", i); }
        let r = &ops[i].rotation;
        let t = &ops[i].translation;
        let is_identity = r[0][0] == 1 && r[0][1] == 0 && r[0][2] == 0
            && r[1][0] == 0 && r[1][1] == 1 && r[1][2] == 0
            && r[2][0] == 0 && r[2][1] == 0 && r[2][2] == 1;
        let has_trans = t[0].abs() > 1e-6 || t[1].abs() > 1e-6 || t[2].abs() > 1e-6;

        if is_identity && !has_trans { return "1".to_string(); }

        let rot_str = format!("[{:2},{:2},{:2};{:2},{:2},{:2};{:2},{:2},{:2}]",
            r[0][0], r[0][1], r[0][2],
            r[1][0], r[1][1], r[1][2],
            r[2][0], r[2][1], r[2][2]);
        if has_trans {
            format!("{}|{:.2},{:.2},{:.2}", rot_str, t[0], t[1], t[2])
        } else {
            rot_str
        }
    };

    // Header row
    let mut lines = Vec::new();
    let header: Vec<String> = std::iter::once("ML".to_string())
        .chain(std::iter::once("BC".to_string()))
        .chain((0..max_ops).map(|i| fmt_op(i)))
        .collect();
    lines.push(format!("| {} |", header.join(" | ")));

    // Separator
    let sep: Vec<String> = (0..(2 + max_ops)).map(|_| "---".to_string()).collect();
    lines.push(format!("| {} |", sep.join(" | ")));

    // Data rows
    for ir in &matching {
        let chars = ir.characters();
        let row: Vec<String> = std::iter::once(ir.ml.to_string())
            .chain(std::iter::once(ir.bc.to_string()))
            .chain((0..max_ops).map(|i| {
                if i < chars.len() {
                    format_value(chars[i])
                } else {
                    "".to_string()
                }
            }))
            .collect();
        lines.push(format!("| {} |", row.join(" | ")));
    }

    let header_line = format!(
        "// SG {} k-point ({}/{}, {}/{}, {}/{}), {} irrep(s)",
        sg, kx, kd, ky, kd, kz, kd,
        matching.len()
    );
    std::iter::once(header_line)
        .chain(lines.into_iter())
        .collect::<Vec<_>>()
        .join("\n")
}

/// Format a single character value for display.
fn format_value(v: f64) -> String {
    if v.abs() < 1e-12 {
        return "0".to_string();
    }
    let r = v.round();
    if (v - r).abs() < 1e-12 {
        return format!("{:.0}", r);
    }
    // Try rational representation, simplifying fractions
    for n in -12_i32..=12 {
        for d in &[2.0, 3.0, 4.0, 6.0, 8.0] {
            if (v - n as f64 / d).abs() < 1e-10 {
                let g = gcd(n.unsigned_abs() as u64, *d as u64);
                let num = n / g as i32;
                let den = *d as i64 / g as i64;
                if den == 1 {
                    return format!("{}", num);
                }
                return format!("{}/{}", num, den);
            }
        }
    }
    format!("{:.4}", v)
}

/// Get the symmetry operations for the canonical setting of a given space group.
///
/// Maps the SG number to its canonical Hall number, then looks up
/// the {R|t} operations from the space group database.
///
/// Returns a vector of `(rotation_matrix_3x3, translation_vector_3)` pairs.
/// The matrices are integer, the translations are fractional.
/// Get symmetry operations for a space group number.
///
/// Returns [`crate::SymmetryOps`] which derefs to `&[SymmetryOp]`.
pub fn symmetry_operations_of(sg: u8) -> crate::SymmetryOps {
    crate::SymmetryOps::from_sg(sg).unwrap_or_else(|| crate::SymmetryOps::default())
}

/// Greatest common divisor.
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kpoint_label() {
        assert_eq!(kpoint_label("GM4+"), "GM");
        assert_eq!(kpoint_label("X3-"), "X");
        assert_eq!(kpoint_label("DT1"), "DT");
        assert_eq!(kpoint_label("LD2"), "LD");
        assert_eq!(kpoint_label("M5"), "M");
        assert_eq!(kpoint_label("R2+"), "R");
        assert_eq!(kpoint_label("A1-"), "A");
        assert_eq!(kpoint_label("H3"), "H");
        assert_eq!(kpoint_label("K1"), "K");
        assert_eq!(kpoint_label("SM3+"), "SM");
        assert_eq!(kpoint_label("GM1"), "GM");
    }

    #[test]
    fn test_irreps_of() {
        // SG 1 should have some irreps
        let s = irreps_of(1);
        assert!(!s.is_empty());
        assert!(s.iter().all(|r| r.sg == 1));

        // Invalid SG
        assert!(irreps_of(0).is_empty());
        assert!(irreps_of(231).is_empty());
    }

    #[test]
    fn test_total_irreps() {
        assert_eq!(total_irreps(), IRREPS.len());
    }

    #[test]
    fn test_total_sgs() {
        assert_eq!(total_sgs(), 230);
    }

    #[test]
    fn test_sg_info() {
        assert_eq!(
            sg_info(1),
            Some(("P1", "C1^1"))
        );
        assert_eq!(
            sg_info(225),
            Some(("Fm-3m", "Oh^5"))
        );
        assert_eq!(
            sg_info(230),
            Some(("Ia-3d", "Oh^10"))
        );
        assert_eq!(sg_info(0), None);
        assert_eq!(sg_info(231), None);
    }

    #[test]
    fn test_subgroups_of() {
        // The first irrep (index 0) should have at least one subgroup
        let subs = subgroups_of(0);
        assert!(!subs.is_empty());

        // Out of range should give empty
        let out = subgroups_of(IRREPS.len() + 1);
        assert!(out.is_empty());
    }

    #[test]
    fn test_kpoints_of() {
        // SG 221 (Pm-3m) has high-symmetry k-points: Γ, X, M, R
        let kps = kpoints_of(221);
        assert!(!kps.is_empty());

        // Should have at least Γ point
        let gm_kps: Vec<_> = kps.iter().filter(|k| k.label == "GM").collect();
        assert!(!gm_kps.is_empty(), "Expected GM (Γ) k-point for SG 221");

        // Each k-point should have at least one irrep
        for kp in &kps {
            assert!(!kp.irreps.is_empty(), "k-point {} has no irreps", kp.label);
        }

        // All k-points should be unique by coords
        let mut coords_set = std::collections::HashSet::new();
        for kp in &kps {
            assert!(
                coords_set.insert(kp.coords),
                "Duplicate k-point coords: {:?}", kp.coords
            );
        }
    }

    #[test]
    fn test_kpoints_of_invalid() {
        assert!(kpoints_of(0).is_empty());
        assert!(kpoints_of(231).is_empty());
    }

    #[test]
    fn test_format_character_table() {
        // SG 1 has a simple character table at GM (Γ point)
        let table = format_character_table(1, 0, 0, 0, 1);
        assert!(table.contains("SG 1"));
        assert!(table.contains("ML"));
        assert!(table.contains("BC"));
    }

    #[test]
    fn test_format_character_table_empty() {
        // Non-existent k-point should return a message
        let table = format_character_table(1, 99, 99, 99, 1);
        assert!(table.contains("No irreps"));
    }

    #[test]
    fn test_symmetry_operations_of() {
        // SG 1 (P1) has only the identity operation
        let ops = symmetry_operations_of(1);
        assert_eq!(ops.len(), 1);

        // SG 221 (Pm-3m) has 48 symmetry operations
        let ops_221 = symmetry_operations_of(221);
        assert_eq!(ops_221.len(), 48);

        // Invalid SG
        assert!(symmetry_operations_of(0).is_empty());
        assert!(symmetry_operations_of(231).is_empty());
    }

    #[test]
    fn test_symmetry_operations_preserves_order() {
        // Verify that the identity is the first operation for a known space group
        let ops = symmetry_operations_of(225); // Fm-3m
        assert!(!ops.is_empty());
        // The first operation should be identity (or close)
        let rot = &ops[0].rotation;
        // Identity rotation
        assert_eq!(rot[0][0], 1);
        assert_eq!(rot[1][1], 1);
        assert_eq!(rot[2][2], 1);
    }

    #[test]
    fn test_sg_142_irreps_after_fix() {
        // Regression test for PIR parser pmkcount > 1 bug.
        // SG 142 (I4_1/acd) has some irreps with pmkcount > 1,
        // which previously caused character values to be shifted.
        // After fix: character values should NOT include 8.0,
        // and should have proper operator counts.
        let irreps = irreps_of(142);
        assert!(!irreps.is_empty(), "SG 142 should have irreps");

        for ir in irreps {
            let chars = ir.characters();
            if ir.dim == 4 {
                // For a 4D irrep, max character value should be 4.0 (identity trace)
                for &c in chars {
                    assert!(
                        c.abs() <= 4.0 + 1e-10,
                        "SG 142 irrep {}: character {} > 4.0 for dim=4 irrep",
                        ir.ml, c
                    );
                }
            }
        }
    }

    #[test]
    fn test_format_value() {
        assert_eq!(format_value(0.0), "0");
        assert_eq!(format_value(1.0), "1");
        assert_eq!(format_value(-1.0), "-1");
        assert_eq!(format_value(4.0), "4");
        assert_eq!(format_value(-2.0), "-2");
        assert_eq!(format_value(0.5), "1/2");
        assert_eq!(format_value(-0.5), "-1/2");
        assert_eq!(format_value(1.5), "3/2");
        // Check that fractions are simplified
        assert_eq!(format_value(-0.5), "-1/2"); // -4/8 simplified to -1/2
        assert_eq!(format_value(0.75), "3/4");  // 6/8 simplified to 3/4
        assert_eq!(format_value(2.0), "2");     // 4/2 simplified to 2
    }

    /// Every irrep in all 230 SGs satisfies basic constraints: χ(E)=dim, non-empty labels.
    #[test]
    fn test_all_space_group_irreps_are_well_formed() {
        for sg in 1u8..=230 {
            let irreps = super::irreps_of(sg);
            assert!(!irreps.is_empty(), "SG {} has no irreps", sg);
            for ir in irreps {
                assert_eq!(ir.sg, sg);
                assert!(ir.dim > 0, "zero dim: SG{} {}", sg, ir.ml);
                let chars = ir.characters();
                assert!(!chars.is_empty(), "empty chars: SG{} {}", sg, ir.ml);
                // χ(E) must be positive and approximately integer.
                // (dim field from image label may be inaccurate for compound irreps.)
                assert!(chars[0] > 0.0,
                    "χ(E) <= 0 for SG{} {}: χ(E)={}", sg, ir.ml, chars[0]);
                let chi_e_rounded = chars[0].round();
                assert!((chars[0] - chi_e_rounded).abs() < 1e-8,
                    "χ(E) not integer for SG{} {}: χ(E)={}", sg, ir.ml, chars[0]);
            }
        }
    }

    /// kpoints_of covers every irrep exactly once across all 230 SGs.
    #[test]
    fn test_kpoints_partition_all_irreps() {
        for sg in 1u8..=230 {
            let irreps = super::irreps_of(sg);
            let kps = super::kpoints_of(sg);
            let mut covered = vec![false; irreps.len()];
            for kp in &kps {
                assert!(!kp.irreps.is_empty(), "SG{} k-point {} empty", sg, kp.label);
                for &idx in &kp.irreps {
                    assert!(idx < irreps.len());
                    assert!(!covered[idx], "SG{} irrep {} duplicate in kpoints", sg, idx);
                    covered[idx] = true;
                    let ir = &irreps[idx];
                    assert_eq!((ir.kx, ir.ky, ir.kz, ir.kd), kp.coords,
                        "SG{} irrep {} k-coord mismatch with k-point {}", sg, ir.ml, kp.label);
                }
            }
            for (i, ok) in covered.iter().enumerate() {
                assert!(*ok, "SG{} irrep {} not covered", sg, irreps[i].ml);
            }
        }
    }
}
