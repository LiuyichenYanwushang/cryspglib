//! Bridge connecting `SpaceGroup` (main spglib API) to the irrep module.
//!
//! Provides convenience methods on `SpaceGroup` for querying irreducible
//! representations, k-points, and character tables without manually
//! extracting the space group number.

use crate::SpaceGroup;
use crate::irrep::query;
use crate::irrep::types::IrrepRecord;
use crate::SymmetryOps;

impl SpaceGroup {
    /// All irreducible representations for this space group.
    ///
    /// Returns a static slice — the data is embedded in the binary.
    pub fn irreps(&self) -> &'static [IrrepRecord] {
        query::irreps_of(self.spacegroup_number as u8)
    }

    /// Unique k-points and their irrep counts for this space group.
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
}
