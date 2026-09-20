//! Isotropy subgroups for a chosen order-parameter direction.
//!
//! Given a parent space group, one of its irreps at a **k**-point, and an
//! order-parameter direction, this module returns the subgroup that leaves that
//! order parameter invariant — together with the geometry that embeds the
//! subgroup in the parent's coordinate system.
//!
//! # Data semantics
//!
//! The records come verbatim from the ISOTROPY suite (`data_isotropy.txt` and
//! `data_magnetic.txt` in the pinned `iso.zip` archive), the machine-readable
//! form of Stokes & Hatch, *Isotropy Subgroups of the 230 Crystallographic
//! Space Groups* (World Scientific, 1988).
//!
//! For every record:
//!
//! - [`IsotropyRecord::direction`] is the order-parameter direction in
//!   components, e.g. `"(a,0,0)"`; [`IsotropyRecord::direction_label`] is the
//!   ISOTROPY label (`"P1"`, `"C1"`, `"4D1"`, …) whose prefix is the subduction
//!   frequency `i(G)` defined in the book.
//! - [`IsotropyRecord::basis`] holds the **primitive** basis vectors of the
//!   subgroup lattice expressed in the parent's conventional basis (rows).  For
//!   centred subgroups this is *not* the conventional basis printed in the
//!   book: e.g. for #16 `P222` → #22 `F222` the stored basis is
//!   `(0,1,1), (1,0,1), (1,1,0)` while the book prints the F-centred
//!   conventional cell `(2,0,0), (0,2,0), (0,0,2)`; both span the same lattice.
//! - [`IsotropyRecord::origin`] is the subgroup origin relative to the parent
//!   origin, given in units of the parent's conventional basis and encoded
//!   exactly as `[x, y, z, d]` = `(x/d, y/d, z/d)`.  This is the book's
//!   "Origin" column (the origin shift, "中心点平移").
//!
//! A parent irrep can have several possible subgroups, one per direction.  The
//! direction strings and ISOTROPY labels are unique within an irrep.
//!
//! # Example
//!
//! ```
//! use cryspglib::irrep::isotropy::{isotropy_subgroup_for_direction, IsotropyDirection};
//!
//! // Pm-3m (221), Γ3+ order parameter along (a,0) → P4/mmm (123)
//! let sub = isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(a,0)"))
//!     .unwrap();
//! assert_eq!(sub.record.sg, 123);
//! assert_eq!(sub.record.symbol, "P4/mmm");
//! ```
//!
//! # Scope
//!
//! Only the directions listed in the ISOTROPY tables can be selected; the
//! stabiliser is looked up, not recomputed from representation matrices.
//! Landau "allowed" filtering (invariant polynomials) and the subduction of
//! parent irreps into the subgroup are separate steps.

use crate::irrep::query;
use crate::irrep::types::{IrrepRecord, IsotropyRecord, KVector, MagneticIsotropyRecord};

// ── Error type ────────────────────────────────────────────────────────────────

/// Errors from isotropy-subgroup queries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IsotropyError {
    /// Space group number out of range (1–230).
    InvalidSpaceGroup(u8),
    /// The space group has no scalar irrep with this Miller–Love label.
    UnknownIrrep { sg: u8, ml: String },
    /// The irrep exists but belongs to a different **k**-point.
    IrrepNotAtKPoint { sg: u8, ml: String, k: KVector },
    /// Double-valued irreps carry no isotropy-subgroup data.
    SpinorUnsupported { sg: u8, ml: String },
    /// No isotropy subgroup matches the requested direction.
    DirectionNotFound {
        sg: u8,
        ml: String,
        direction: String,
    },
    /// The direction matches more than one isotropy subgroup.
    DirectionAmbiguous {
        sg: u8,
        ml: String,
        direction: String,
        matches: usize,
    },
    /// Index out of range for this irrep's isotropy subgroup list.
    SubgroupIndexOutOfRange {
        sg: u8,
        ml: String,
        index: usize,
        len: usize,
    },
    /// The parent space group has no centering information in the database.
    MissingCentering { sg: u8 },
    /// The subgroup basis is singular and cannot describe a lattice.
    SingularSubgroupBasis { basis: [[i32; 3]; 3] },
}

impl std::fmt::Display for IsotropyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSpaceGroup(sg) => {
                write!(f, "space group number {sg} is out of range (1-230)")
            }
            Self::UnknownIrrep { sg, ml } => {
                write!(f, "space group {sg} has no scalar irrep labelled {ml}")
            }
            Self::IrrepNotAtKPoint { sg, ml, k } => write!(
                f,
                "irrep {ml} of space group {sg} is not at k = ({}/{}, {}/{}, {}/{})",
                k.numerators[0],
                k.denominator,
                k.numerators[1],
                k.denominator,
                k.numerators[2],
                k.denominator
            ),
            Self::SpinorUnsupported { sg, ml } => write!(
                f,
                "irrep {ml} of space group {sg} is double-valued; \
                 isotropy subgroups are defined for single-valued irreps only"
            ),
            Self::DirectionNotFound { sg, ml, direction } => write!(
                f,
                "irrep {ml} of space group {sg} has no isotropy subgroup for direction {direction}"
            ),
            Self::DirectionAmbiguous {
                sg,
                ml,
                direction,
                matches,
            } => write!(
                f,
                "direction {direction} of irrep {ml} (space group {sg}) matches {matches} \
                 isotropy subgroups; select one with IsotropyDirection::Index"
            ),
            Self::SubgroupIndexOutOfRange { sg, ml, index, len } => write!(
                f,
                "index {index} is out of range: irrep {ml} of space group {sg} has {len} \
                 isotropy subgroup(s)"
            ),
            Self::MissingCentering { sg } => {
                write!(f, "space group {sg} has no centering information")
            }
            Self::SingularSubgroupBasis { basis } => write!(
                f,
                "subgroup basis {basis:?} is singular and describes no lattice"
            ),
        }
    }
}

impl std::error::Error for IsotropyError {}

// ── Direction selector ────────────────────────────────────────────────────────

/// How to pick one isotropy subgroup of an irrep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsotropyDirection<'a> {
    /// Component description exactly as stored, e.g. `"(a,0,0)"`.
    Descriptor(&'a str),
    /// ISOTROPY direction label, e.g. `"P1"` or `"4D1"`.
    Label(&'a str),
    /// Position in this irrep's isotropy subgroup list (0-based, table order).
    Index(usize),
}

// ── Result types ──────────────────────────────────────────────────────────────

/// An isotropy subgroup with the irrep context it was derived from.
#[derive(Debug, Clone, Copy)]
pub struct IsotropySubgroup {
    /// Parent space group number (1–230).
    pub parent_sg: u8,
    /// Parent Hermann–Mauguin symbol.
    pub parent_symbol: &'static str,
    /// Wave vector of the condensing irrep.
    pub k: KVector,
    /// Miller–Love label of the condensing irrep.
    pub irrep_ml: &'static str,
    /// Bradley–Cracknell label of the condensing irrep.
    pub irrep_bc: &'static str,
    /// Dimension of the condensing irrep.
    pub irrep_dim: u8,
    /// Position of this record in the generated isotropy table.
    pub ordinal: usize,
    /// The subgroup record, including basis, origin, domains and arms.
    pub record: IsotropyRecord,
}

/// A magnetic isotropy subgroup with the irrep context it was derived from.
#[derive(Debug, Clone, Copy)]
pub struct MagneticIsotropySubgroup {
    /// Parent space group number (1–230) of the non-magnetic parent.
    pub parent_sg: u8,
    /// Parent Hermann–Mauguin symbol.
    pub parent_symbol: &'static str,
    /// Wave vector of the condensing irrep.
    pub k: KVector,
    /// Miller–Love label of the condensing irrep.
    pub irrep_ml: &'static str,
    /// Bradley–Cracknell label of the condensing irrep.
    pub irrep_bc: &'static str,
    /// Dimension of the condensing irrep.
    pub irrep_dim: u8,
    /// Position of this record in the generated magnetic isotropy table.
    pub ordinal: usize,
    /// The magnetic subgroup record, including basis, origin and UNI number.
    pub record: MagneticIsotropyRecord,
}

// ── Lookup helpers ────────────────────────────────────────────────────────────

/// Find a single-valued irrep of a space group by Miller–Love label.
pub fn irrep_by_label(sg: u8, ml: &str) -> Result<&'static IrrepRecord, IsotropyError> {
    if sg == 0 || sg > 230 {
        return Err(IsotropyError::InvalidSpaceGroup(sg));
    }
    let found = query::irreps_of(sg).iter().find(|ir| ir.ml == ml);
    let Some(irrep) = found else {
        return Err(IsotropyError::UnknownIrrep {
            sg,
            ml: ml.to_string(),
        });
    };
    if irrep.spinor {
        return Err(IsotropyError::SpinorUnsupported {
            sg,
            ml: ml.to_string(),
        });
    }
    Ok(irrep)
}

/// Whether two rational wave vectors denote the same point (up to reduction).
pub fn k_vectors_agree(a: KVector, b: KVector) -> bool {
    let (an, ad) = reduce_k(a);
    let (bn, bd) = reduce_k(b);
    an == bn && ad == bd
}

fn reduce_k(k: KVector) -> ([i32; 3], i32) {
    let d = k.denominator as i32;
    if d == 0 {
        return (
            [
                k.numerators[0] as i32,
                k.numerators[1] as i32,
                k.numerators[2] as i32,
            ],
            1,
        );
    }
    let mut g = d.abs();
    for n in k.numerators {
        g = gcd_i32(g, (n as i32).abs());
    }
    if g == 0 {
        g = 1;
    }
    let sign = if d < 0 { -1 } else { 1 };
    (
        [
            sign * k.numerators[0] as i32 / g,
            sign * k.numerators[1] as i32 / g,
            sign * k.numerators[2] as i32 / g,
        ],
        sign * d / g,
    )
}

fn gcd_i32(mut a: i32, mut b: i32) -> i32 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.abs()
}

fn wrap(irrep: &'static IrrepRecord, local: usize) -> IsotropySubgroup {
    IsotropySubgroup {
        parent_sg: irrep.sg,
        parent_symbol: query::sg_info(irrep.sg).map(|info| info.0).unwrap_or("?"),
        k: irrep.k_vector(),
        irrep_ml: irrep.ml,
        irrep_bc: irrep.bc,
        irrep_dim: irrep.dim,
        ordinal: irrep._iso_start as usize + local,
        record: irrep.subgroups()[local],
    }
}

fn wrap_magnetic(
    irrep: &'static IrrepRecord,
    local: usize,
    record: MagneticIsotropyRecord,
) -> MagneticIsotropySubgroup {
    MagneticIsotropySubgroup {
        parent_sg: irrep.sg,
        parent_symbol: query::sg_info(irrep.sg).map(|info| info.0).unwrap_or("?"),
        k: irrep.k_vector(),
        irrep_ml: irrep.ml,
        irrep_bc: irrep.bc,
        irrep_dim: irrep.dim,
        ordinal: irrep._mag_iso_start as usize + local,
        record,
    }
}

// ── Non-magnetic queries ──────────────────────────────────────────────────────

/// All isotropy subgroups of one irrep, in table order.
pub fn isotropy_subgroups(sg: u8, ml: &str) -> Result<Vec<IsotropySubgroup>, IsotropyError> {
    let irrep = irrep_by_label(sg, ml)?;
    Ok((0..irrep.subgroups().len())
        .map(|local| wrap(irrep, local))
        .collect())
}

/// All isotropy subgroups of one irrep, after checking the **k**-point.
pub fn isotropy_subgroups_at_k(
    sg: u8,
    k: KVector,
    ml: &str,
) -> Result<Vec<IsotropySubgroup>, IsotropyError> {
    let irrep = irrep_by_label(sg, ml)?;
    if !k_vectors_agree(irrep.k_vector(), k) {
        return Err(IsotropyError::IrrepNotAtKPoint {
            sg,
            ml: ml.to_string(),
            k,
        });
    }
    Ok((0..irrep.subgroups().len())
        .map(|local| wrap(irrep, local))
        .collect())
}

/// The isotropy subgroup selected by direction, ISOTROPY label, or index.
pub fn isotropy_subgroup_for_direction(
    sg: u8,
    ml: &str,
    direction: IsotropyDirection<'_>,
) -> Result<IsotropySubgroup, IsotropyError> {
    let irrep = irrep_by_label(sg, ml)?;
    let local = select_local_index(
        irrep,
        ml,
        direction,
        |record| record.direction,
        |record| record.direction_label,
    )?;
    Ok(wrap(irrep, local))
}

// ── Magnetic queries ──────────────────────────────────────────────────────────

/// All magnetic isotropy subgroups of one irrep, in table order.
pub fn magnetic_isotropy_subgroups(
    sg: u8,
    ml: &str,
) -> Result<Vec<MagneticIsotropySubgroup>, IsotropyError> {
    let irrep = irrep_by_label(sg, ml)?;
    Ok(irrep
        .magnetic_subgroups()
        .iter()
        .enumerate()
        .map(|(local, record)| wrap_magnetic(irrep, local, *record))
        .collect())
}

/// The magnetic isotropy subgroup selected by direction label or index.
pub fn magnetic_isotropy_subgroup_for_direction(
    sg: u8,
    ml: &str,
    direction: IsotropyDirection<'_>,
) -> Result<MagneticIsotropySubgroup, IsotropyError> {
    let irrep = irrep_by_label(sg, ml)?;
    let records = irrep.magnetic_subgroups();
    let local = match direction {
        IsotropyDirection::Index(index) => {
            if index >= records.len() {
                return Err(IsotropyError::SubgroupIndexOutOfRange {
                    sg,
                    ml: ml.to_string(),
                    index,
                    len: records.len(),
                });
            }
            index
        }
        IsotropyDirection::Descriptor(descriptor) => {
            select_unique(records, descriptor, |record| record.direction, sg, ml)?
        }
        IsotropyDirection::Label(label) => {
            select_unique(records, label, |record| record.direction, sg, ml)?
        }
    };
    Ok(wrap_magnetic(irrep, local, records[local]))
}

// ── Selection internals ───────────────────────────────────────────────────────

fn select_local_index(
    irrep: &'static IrrepRecord,
    ml: &str,
    direction: IsotropyDirection<'_>,
    descriptor_of: impl Fn(&IsotropyRecord) -> &'static str,
    label_of: impl Fn(&IsotropyRecord) -> &'static str,
) -> Result<usize, IsotropyError> {
    let records = irrep.subgroups();
    match direction {
        IsotropyDirection::Index(index) => {
            if index >= records.len() {
                return Err(IsotropyError::SubgroupIndexOutOfRange {
                    sg: irrep.sg,
                    ml: ml.to_string(),
                    index,
                    len: records.len(),
                });
            }
            Ok(index)
        }
        IsotropyDirection::Descriptor(descriptor) => {
            select_unique(records, descriptor, descriptor_of, irrep.sg, ml)
        }
        IsotropyDirection::Label(label) => select_unique(records, label, label_of, irrep.sg, ml),
    }
}

fn select_unique<T: Copy>(
    records: &[T],
    query: &str,
    key_of: impl Fn(&T) -> &'static str,
    sg: u8,
    ml: &str,
) -> Result<usize, IsotropyError> {
    let matches: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|(_, record)| key_of(record) == query)
        .map(|(index, _)| index)
        .collect();
    match matches.len() {
        0 => Err(IsotropyError::DirectionNotFound {
            sg,
            ml: ml.to_string(),
            direction: query.to_string(),
        }),
        1 => Ok(matches[0]),
        matches => Err(IsotropyError::DirectionAmbiguous {
            sg,
            ml: ml.to_string(),
            direction: query.to_string(),
            matches,
        }),
    }
}

// ── Derived geometry ──────────────────────────────────────────────────────────

/// Relative volume of the subgroup's primitive cell to the parent's.
///
/// The stored basis spans the subgroup's primitive cell inside the parent's
/// conventional cell, so the volume ratio is `|det W| · Z(parent)` where
/// `Z(parent)` is the parent's centering multiplicity.  This is the "Size"
/// column of Stokes & Hatch (1988).
pub fn subgroup_index_in_parent(parent_sg: u8, basis: [[i32; 3]; 3]) -> Result<u32, IsotropyError> {
    let det = det3(basis);
    if det == 0 {
        return Err(IsotropyError::SingularSubgroupBasis { basis });
    }
    let centering = centering_multiplicity(parent_sg)?;
    Ok(det.unsigned_abs() * centering)
}

/// Centering multiplicity `Z` of a space group (1, 2, 3 or 4).
pub fn centering_multiplicity(sg: u8) -> Result<u32, IsotropyError> {
    use crate::spg_database::Centering;
    if sg == 0 || sg > 230 {
        return Err(IsotropyError::InvalidSpaceGroup(sg));
    }
    let hall = crate::irrep::generated_data::SG_DATA_HALL[sg as usize] as usize;
    let centering = crate::spg_database::get_spacegroup_type(hall).centering;
    match centering {
        Centering::Primitive => Ok(1),
        Centering::Body | Centering::AFace | Centering::BFace | Centering::CFace => Ok(2),
        Centering::RCenter => Ok(3),
        Centering::Face => Ok(4),
        Centering::Error => Err(IsotropyError::MissingCentering { sg }),
    }
}

fn det3(m: [[i32; 3]; 3]) -> i32 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

// ── Formatting ────────────────────────────────────────────────────────────────

/// Markdown table of the isotropy subgroups of one irrep, including the
/// subgroup geometry (basis, origin, index).
pub fn format_isotropy_subgroups(sg: u8, ml: &str) -> Result<String, IsotropyError> {
    let subgroups = isotropy_subgroups(sg, ml)?;
    let irrep = irrep_by_label(sg, ml)?;
    let mut lines = Vec::new();
    lines.push(format!(
        "// SG {} {} irrep {} at k=({}/{}, {}/{}, {}/{}), {} isotropy subgroup(s)",
        sg,
        irrep.bc,
        irrep.ml,
        irrep.kx,
        irrep.kd,
        irrep.ky,
        irrep.kd,
        irrep.kz,
        irrep.kd,
        subgroups.len()
    ));
    lines.push(
        "| # | Subgroup | Label | Direction | Dir | Size | Basis | Origin | Domains | Arms |"
            .to_string(),
    );
    lines.push(
        "|---|----------|-------|-----------|-----|------|-------|--------|---------|------|"
            .to_string(),
    );
    for (index, sub) in subgroups.iter().enumerate() {
        let record = &sub.record;
        let size = subgroup_index_in_parent(sg, record.basis)
            .map(|value| value.to_string())
            .unwrap_or_else(|_| "?".to_string());
        lines.push(format!(
            "| {} | #{} {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            index,
            record.sg,
            record.symbol,
            record.direction_label,
            record.direction,
            record.direction_dim,
            size,
            format_basis(record.basis),
            format_origin(record.origin),
            record.domains,
            record.arms,
        ));
    }
    Ok(lines.join("\n"))
}

/// Markdown table of the magnetic isotropy subgroups of one irrep.
pub fn format_magnetic_isotropy_subgroups(sg: u8, ml: &str) -> Result<String, IsotropyError> {
    let subgroups = magnetic_isotropy_subgroups(sg, ml)?;
    let irrep = irrep_by_label(sg, ml)?;
    let mut lines = Vec::new();
    lines.push(format!(
        "// SG {} {} irrep {} at k=({}/{}, {}/{}, {}/{}), {} magnetic isotropy subgroup(s)",
        sg,
        irrep.bc,
        irrep.ml,
        irrep.kx,
        irrep.kd,
        irrep.ky,
        irrep.kd,
        irrep.kz,
        irrep.kd,
        subgroups.len()
    ));
    lines.push("| # | UNI | BNS | Label | Direction | Basis | Origin |".to_string());
    lines.push("|---|-----|-----|-------|-----------|-------|--------|".to_string());
    for (index, sub) in subgroups.iter().enumerate() {
        let record = &sub.record;
        lines.push(format!(
            "| {} | {} | {} | {} | {} | {} | {} |",
            index,
            record.mag_sg,
            record.bns_label,
            record.iso_label,
            record.direction,
            format_basis(record.basis),
            format_origin(record.origin),
        ));
    }
    Ok(lines.join("\n"))
}

fn format_basis(basis: [[i32; 3]; 3]) -> String {
    format!(
        "({},{},{}),({},{},{}),({},{},{})",
        basis[0][0],
        basis[0][1],
        basis[0][2],
        basis[1][0],
        basis[1][1],
        basis[1][2],
        basis[2][0],
        basis[2][1],
        basis[2][2],
    )
}

fn format_origin(origin: [i32; 4]) -> String {
    let shift = crate::irrep::types::decode_origin(origin);
    if origin[3] == 1 {
        return format!("({},{},{})", origin[0], origin[1], origin[2]);
    }
    format!(
        "({},{},{}) = ({}, {}, {})",
        origin[0], origin[1], origin[2], shift[0], shift[1], shift[2]
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pm3m_gm3plus_along_a0_gives_p4mmm() {
        let sub =
            isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(a,0)"))
                .expect("GM3+ has an (a,0) direction");
        assert_eq!(sub.record.sg, 123);
        assert_eq!(sub.record.symbol, "P4/mmm");
        assert_eq!(sub.record.direction_dim, 2);
        assert_eq!(sub.record.direction_free, 1);
        assert_eq!(sub.record.direction_label, "P1");
        assert_eq!(sub.record.basis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        assert_eq!(sub.record.origin, [0, 0, 0, 1]);
        assert_eq!(subgroup_index_in_parent(221, sub.record.basis).unwrap(), 1);
    }

    #[test]
    fn p222_r1_reaches_centred_f222_with_primitive_basis() {
        let sub = isotropy_subgroup_for_direction(16, "R1", IsotropyDirection::Index(0))
            .expect("R1 has an isotropy subgroup");
        assert_eq!(sub.record.sg, 22);
        assert_eq!(sub.record.basis, [[0, 1, 1], [1, 0, 1], [1, 1, 0]]);
        assert_eq!(sub.record.origin, [0, 0, 0, 1]);
        // The book prints Size 2 for this transition.
        assert_eq!(subgroup_index_in_parent(16, sub.record.basis).unwrap(), 2);
    }

    #[test]
    fn direction_lookup_can_use_the_isotropy_label() {
        let by_label = isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Label("P1"))
            .expect("P1 selects the (a,0) direction");
        let by_descriptor =
            isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(a,0)"))
                .expect("(a,0) selects the same record");
        assert_eq!(by_label.ordinal, by_descriptor.ordinal);
    }

    #[test]
    fn unknown_direction_and_irrep_fail_closed() {
        let missing =
            isotropy_subgroup_for_direction(221, "GM3+", IsotropyDirection::Descriptor("(zzz)"));
        assert!(matches!(
            missing,
            Err(IsotropyError::DirectionNotFound { .. })
        ));

        let unknown = isotropy_subgroups(221, "NOPE1");
        assert!(matches!(unknown, Err(IsotropyError::UnknownIrrep { .. })));

        let bad_sg = isotropy_subgroups(231, "GM1");
        assert!(matches!(bad_sg, Err(IsotropyError::InvalidSpaceGroup(231))));
    }

    #[test]
    fn spinor_irreps_are_rejected_without_panicking() {
        let spinor = query::irreps_of(1)
            .iter()
            .find(|irrep| irrep.spinor)
            .expect("SG 1 has a spinor irrep");
        let result = isotropy_subgroups(1, spinor.ml);
        assert!(matches!(
            result,
            Err(IsotropyError::SpinorUnsupported { .. })
        ));
    }

    #[test]
    fn k_point_mismatch_is_reported() {
        let wrong_k = KVector::new([1, 0, 0], 2);
        let result = isotropy_subgroups_at_k(221, wrong_k, "GM3+");
        assert!(matches!(
            result,
            Err(IsotropyError::IrrepNotAtKPoint { .. })
        ));
        let right_k = KVector::new([0, 0, 0], 1);
        assert!(isotropy_subgroups_at_k(221, right_k, "GM3+").is_ok());
    }

    #[test]
    fn magnetic_table_is_keyed_by_non_magnetic_parent_irreps() {
        let spinor_free = magnetic_isotropy_subgroups(221, "GM3+").expect("mag subgroups exist");
        assert!(!spinor_free.is_empty());
        assert!(
            spinor_free
                .iter()
                .all(|sub| (1..=1651).contains(&sub.record.mag_sg))
        );
    }

    #[test]
    fn centering_multiplicity_matches_known_space_groups() {
        assert_eq!(centering_multiplicity(1).unwrap(), 1); // P1
        assert_eq!(centering_multiplicity(16).unwrap(), 1); // P222
        assert_eq!(centering_multiplicity(22).unwrap(), 4); // F222
        assert_eq!(centering_multiplicity(225).unwrap(), 4); // Fm-3m
        assert_eq!(centering_multiplicity(167).unwrap(), 3); // R-3c
        assert_eq!(centering_multiplicity(229).unwrap(), 2); // Im-3m
    }

    #[test]
    fn formatting_includes_geometry_columns() {
        let table = format_isotropy_subgroups(221, "GM3+").expect("table renders");
        assert!(table.contains("| Size |"));
        assert!(table.contains("#123 P4/mmm"));
        assert!(table.contains("(1,0,0),(0,1,0),(0,0,1)"));
    }
}
