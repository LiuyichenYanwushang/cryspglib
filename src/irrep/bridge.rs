//! Bridge connecting `SpaceGroup` (main spglib API) to the irrep module.
//!
//! Provides convenience methods on `SpaceGroup` for querying irreducible
//! representations, k-points, and character tables without manually
//! extracting the space group number.

use crate::SpaceGroup;
use crate::irrep::query;
use crate::irrep::query::{IsotropyEntry, MagneticIsotropyEntry};
use crate::irrep::types::IrrepRecord;
use crate::irrep::types::generated_data::SG_DATA_HALL;
use crate::SymmetryOps;

/// Get H_ops in the same order as the stored irrep data (CHARACTERS, PIR_ROTS, etc.)
///
/// For most SGs this matches the spglib default Hall setting.
/// Use this to avoid runtime rotation-matching when the data and H_ops are aligned.
pub fn canonical_hall_ops(sg: u8) -> Result<SymmetryOps, crate::SymError> {
    if sg == 0 || sg > 230 {
        return Err(crate::SymError::SpacegroupSearchFailed);
    }
    let hall = SG_DATA_HALL[sg as usize] as usize;
    if hall == 0 {
        // No canonical Hall recorded — fall back to SG-based lookup
        return SymmetryOps::from_sg(sg);
    }
    SymmetryOps::from_database(hall).or_else(|_| SymmetryOps::from_sg(sg))
}

impl SpaceGroup {
    /// All irreducible representations for this space group.
    pub fn irreps(&self) -> &'static [IrrepRecord] {
        query::irreps_of(self.spacegroup_number as u8)
    }

    /// Unique k-points and their irrep indices for this space group.
    pub fn kpoints(&self) -> Vec<query::KPointSummary> {
        query::kpoints_of(self.spacegroup_number as u8)
    }

    /// Irreps at a specific k-point label (e.g. `"GM"`, `"X"`, `"R"`).
    pub fn irreps_at_k(&self, label: &str) -> Vec<&'static IrrepRecord> {
        self.irreps()
            .iter()
            .filter(|r| r.k_label() == label)
            .collect()
    }

    /// Irreps at specific k-point coordinates (fractional, common denominator).
    pub fn irreps_at_coords(&self, kx: i8, ky: i8, kz: i8, kd: i8) -> Vec<&'static IrrepRecord> {
        self.irreps()
            .iter()
            .filter(|r| r.kx == kx && r.ky == ky && r.kz == kz && r.kd == kd)
            .collect()
    }

    /// Formatted character table at the given k-point coordinates.
    pub fn character_table(&self, kx: i8, ky: i8, kz: i8, kd: i8) -> String {
        query::format_character_table(self.spacegroup_number as u8, kx, ky, kz, kd)
    }

    /// Space group info — (Hermann-Mauguin symbol, Schoenflies symbol).
    pub fn sg_info(&self) -> Option<(&'static str, &'static str)> {
        query::sg_info(self.spacegroup_number as u8)
    }

    /// Symmetry operations in spglib order (time_reversal all `false`).
    pub fn symmetry_ops(&self) -> SymmetryOps {
        let n_op = self.n_operations;
        let mut operations = Vec::with_capacity(n_op);
        for i in 0..n_op {
            operations.push(crate::SymmetryOp {
                rotation: self.rotations[i],
                translation: self.translations[i],
                time_reversal: false,
            });
        }
        SymmetryOps { operations }
    }

    /// All isotropy subgroups across all scalar irreps of this space group.
    ///
    /// Each entry includes the source irrep (k-point, label, dimension).
    pub fn isotropy_subgroups(&self) -> Vec<IsotropyEntry> {
        query::isotropy_subgroups_of(self.spacegroup_number as u8)
    }

    /// All magnetic isotropy subgroups across all scalar irreps of this space group.
    pub fn magnetic_isotropy_subgroups(&self) -> Vec<MagneticIsotropyEntry> {
        query::magnetic_isotropy_subgroups_of(self.spacegroup_number as u8)
    }

    /// Isotropy subgroups at a specific k-point for this space group.
    pub fn isotropy_at_k(&self, kx: i8, ky: i8, kz: i8, kd: i8) -> String {
        query::format_isotropy_table(self.spacegroup_number as u8, kx, ky, kz, kd)
    }

    /// Magnetic isotropy subgroups at a specific k-point for this space group.
    pub fn magnetic_isotropy_at_k(&self, kx: i8, ky: i8, kz: i8, kd: i8) -> String {
        query::format_magnetic_isotropy_table(self.spacegroup_number as u8, kx, ky, kz, kd)
    }

    /// Irreps at a specific k-point label with their isotropy subgroups.
    pub fn irreps_with_isotropy_at_k(&self, label: &str) -> Vec<(&'static IrrepRecord, String)> {
        self.irreps_at_k(label)
            .into_iter()
            .map(|ir| {
                let subs = ir.subgroups();
                let desc = if subs.is_empty() {
                    "(none)".to_string()
                } else {
                    subs.iter()
                        .map(|s| format!("#{} {} dir={}", s.sg, s.symbol, s.direction))
                        .collect::<Vec<_>>()
                        .join("; ")
                };
                (ir, desc)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_hall_ops_reports_invalid_space_group() {
        assert!(matches!(
            canonical_hall_ops(0),
            Err(crate::SymError::SpacegroupSearchFailed)
        ));
        assert!(matches!(
            canonical_hall_ops(231),
            Err(crate::SymError::SpacegroupSearchFailed)
        ));
        assert_eq!(canonical_hall_ops(221).unwrap().len(), 48);
    }
}
