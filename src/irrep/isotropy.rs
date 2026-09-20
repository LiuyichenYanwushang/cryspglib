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
//!   ISOTROPY direction label (`"P1"`, `"C1"`, `"S1"`, `"4D1"`, …).  Its
//!   leading letter is the direction *type* prefix used by Stokes & Hatch; it is
//!   not a shortcut for the subduction frequency `i(G)` of the condensing irrep
//!   (use [`identity_subduction`] for that).
//! - [`IsotropyRecord::basis`] holds a **primitive** basis of the subgroup
//!   lattice (rows), and [`IsotropyRecord::origin`] the origin shift of the
//!   subgroup setting, both expressed in the **parent's primitive-cell frame**:
//!   the frame the ISOTROPY program works in internally.  For a primitive
//!   parent this coincides with the conventional basis, but for `C`, `A`, `B`,
//!   `I`, `F` and `R` lattices it does not — e.g. for #167 `R-3c` the stored
//!   origin `(0, 1/2, 0)` is `(-1/6, 1/6, 1/6)` in the hexagonal conventional
//!   basis, exactly as the book prints it.
//! - [`IsotropyRecord::origin`] is encoded exactly as `[x, y, z, d]` meaning
//!   `(x/d, y/d, z/d)`; use [`origin_shift_in_parent_conventional`] (or
//!   [`IsotropyRecord::origin_shift`] for the raw frame) to work with it.
//! - [`subgroup_size`] is the relative primitive-cell volume, i.e. `|det W|`,
//!   which is the "Size" column of the ISOTROPY tables.
//!
//! Both frames are available: [`basis_in_parent_conventional`] and
//! [`origin_shift_in_parent_conventional`] convert the stored values with the
//! parent's primitive basis, and are checked against the bundled ISOTROPY `iso`
//! binary by `scripts/verify_isotropy_oracle.py` (which keeps its own copy of
//! the primitive-basis table so the two can disagree).
//!
//! A parent irrep can have several possible subgroups, one per direction; both
//! the direction strings and the ISOTROPY labels are unique within an irrep
//! (checked for every irrep of both tables by `tests/isotropy_geometry.rs`).
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
//! The basis is a primitive basis of the subgroup lattice, not the subgroup's
//! ITA conventional cell (which additionally depends on the subgroup centering).
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
    /// Global isotropy record index out of range.
    InvalidIsotropyRecord { ordinal: usize, len: usize },
    /// Origin shift encoding with a non-positive denominator.
    InvalidOrigin { origin: [i32; 4] },
    /// Subgroup basis determinant that does not fit a lattice size.
    SubgroupSizeOverflow { determinant: i64 },
    /// A generated subduction table pointed outside its irrep table.
    SubductionIndexOutOfRange { entry: usize, index: usize },
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
            Self::InvalidIsotropyRecord { ordinal, len } => write!(
                f,
                "isotropy record index {ordinal} is out of range (0-{})",
                len.saturating_sub(1)
            ),
            Self::InvalidOrigin { origin } => {
                write!(f, "origin shift {origin:?} has a non-positive denominator")
            }
            Self::SubgroupSizeOverflow { determinant } => write!(
                f,
                "subgroup basis determinant {determinant} does not fit a lattice size"
            ),
            Self::SubductionIndexOutOfRange { entry, index } => write!(
                f,
                "subduction entry {entry} references irrep index {index}, which is \
                 outside the generated table"
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
    match (reduce_k(a), reduce_k(b)) {
        (Some((an, ad)), Some((bn, bd))) => an == bn && ad == bd,
        _ => false,
    }
}

/// Reduce a wave vector to lowest terms with a positive denominator.
///
/// Returns `None` for a zero denominator, which is not a representable wave
/// vector; callers must not treat it as an integer vector.
fn reduce_k(k: KVector) -> Option<([i32; 3], i32)> {
    let d = k.denominator as i32;
    if d == 0 {
        return None;
    }
    let mut g = d.abs();
    for n in k.numerators {
        g = gcd_i32(g, (n as i32).abs());
    }
    if g == 0 {
        g = 1;
    }
    let sign = if d < 0 { -1 } else { 1 };
    Some((
        [
            sign * k.numerators[0] as i32 / g,
            sign * k.numerators[1] as i32 / g,
            sign * k.numerators[2] as i32 / g,
        ],
        sign * d / g,
    ))
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

// ── Subduction of parent irreps into the subgroup ─────────────────────────────

/// A parent irrep whose restriction contains the trivial irrep of the isotropy
/// subgroup, with the multiplicity `i(G)` used by Landau theory.
///
/// This is the ISOTROPY "subduction frequency" table: the ISO program prints
/// exactly these rows for `SHOW FREQ` (with `SHOW FREQ DIR` adding the domain
/// index).  It answers "which parent irreps become totally
/// symmetric in the distorted phase".
///
/// The full decomposition of a parent irrep into *all* irreps of the subgroup
/// is a different quantity and is not part of this dataset.
#[derive(Debug, Clone, Copy)]
pub struct IdentitySubduction {
    /// Space group of the parent irrep.
    pub parent_sg: u8,
    /// Miller–Love label of the parent irrep.
    pub parent_ml: &'static str,
    /// Bradley–Cracknell label of the parent irrep.
    pub parent_bc: &'static str,
    /// Wave vector of the parent irrep.
    pub parent_k: KVector,
    /// Dimension of the parent irrep.
    pub parent_dim: u8,
    /// Number of times the parent irrep subduces the trivial irrep of the
    /// subgroup (`i(G)` in Stokes & Hatch, 1988).
    pub frequency: u16,
    /// Domain index of the subduction entry.
    pub domain: u16,
    /// Direction label the entry is anchored at (the `Dir` column of the
    /// ISOTROPY program's `SHOW FREQ DIR` output).
    pub direction_label: &'static str,
}

/// A double-valued (spinor) parent irrep whose restriction contains the trivial
/// irrep of the isotropy subgroup.
///
/// The ISO program prints these on the same `SHOW FREQ [DIR]` line as the
/// scalar entries; the double-valued tables carry Miller–Love labels only and
/// no wave vector, so they are reported separately.
#[derive(Debug, Clone, Copy)]
pub struct DoubleValuedSubduction {
    /// Space group of the double-valued parent irrep.
    pub parent_sg: u8,
    /// Miller–Love label of the double-valued parent irrep.
    pub parent_ml: &'static str,
    /// Subduction frequency `i(G)`.
    pub frequency: u16,
}

impl IsotropySubgroup {
    /// Parent irreps that subduce the trivial irrep of this subgroup.
    pub fn identity_subduction(&self) -> Result<Vec<IdentitySubduction>, IsotropyError> {
        identity_subduction(self.ordinal)
    }

    /// Double-valued (spinor) parent irreps that subduce the trivial irrep of
    /// this subgroup.
    pub fn double_valued_subduction(&self) -> Result<Vec<DoubleValuedSubduction>, IsotropyError> {
        double_valued_subduction(self.ordinal)
    }
}

/// Parent irreps that subduce the trivial irrep of an isotropy record.
///
/// `ordinal` is the global isotropy record index carried by
/// [`IsotropySubgroup::ordinal`]; [`IsotropySubgroup::identity_subduction`] is
/// the ergonomic entry point.
pub fn identity_subduction(ordinal: usize) -> Result<Vec<IdentitySubduction>, IsotropyError> {
    let ranges = &crate::irrep::generated_data::ISOTROPY_SUBDUCE_RANGES;
    let records = ranges.len().saturating_sub(1);
    if ordinal >= records {
        return Err(IsotropyError::InvalidIsotropyRecord {
            ordinal,
            len: records,
        });
    }
    // Ranges are 1-based into the subduction tables.
    let start = ranges[ordinal] as usize - 1;
    let end = ranges[ordinal + 1] as usize - 1;
    let irreps = &crate::irrep::generated_data::IRREPS;
    let mut result = Vec::with_capacity(end.saturating_sub(start));
    for entry in start..end {
        let irrep_index = crate::irrep::generated_data::ISOTROPY_SUBDUCE_IRREP[entry] as usize;
        let irrep =
            irreps
                .get(irrep_index - 1)
                .ok_or(IsotropyError::SubductionIndexOutOfRange {
                    entry,
                    index: irrep_index,
                })?;
        let label_index = crate::irrep::generated_data::ISOTROPY_SUBDUCE_DIRECTION[entry] as usize;
        let direction_label = crate::irrep::generated_data::ISOTROPY_DIRECTION_LABELS
            .get(label_index)
            .copied()
            .unwrap_or("?");
        result.push(IdentitySubduction {
            parent_sg: irrep.sg,
            parent_ml: irrep.ml,
            parent_bc: irrep.bc,
            parent_k: irrep.k_vector(),
            parent_dim: irrep.dim,
            frequency: crate::irrep::generated_data::ISOTROPY_SUBDUCE_FREQUENCY[entry] as u16,
            domain: crate::irrep::generated_data::ISOTROPY_SUBDUCE_DOMAIN[entry],
            direction_label,
        });
    }
    Ok(result)
}

/// Double-valued (spinor) parent irreps that subduce the trivial irrep of an
/// isotropy record.
pub fn double_valued_subduction(
    ordinal: usize,
) -> Result<Vec<DoubleValuedSubduction>, IsotropyError> {
    let ranges = &crate::irrep::generated_data::ISOTROPY_W_SUBDUCE_RANGES;
    let records = ranges.len().saturating_sub(1);
    if ordinal >= records {
        return Err(IsotropyError::InvalidIsotropyRecord {
            ordinal,
            len: records,
        });
    }
    let start = ranges[ordinal] as usize;
    let end = ranges[ordinal + 1] as usize;
    let labels = &crate::irrep::generated_data::IRREP_W_LABELS;
    let groups = &crate::irrep::generated_data::IRREP_W_SPACE_GROUP;
    let mut result = Vec::with_capacity(end.saturating_sub(start));
    for entry in start..end {
        let irrep_index = crate::irrep::generated_data::ISOTROPY_W_SUBDUCE_IRREP[entry] as usize;
        let (Some(label), Some(sg)) = (labels.get(irrep_index - 1), groups.get(irrep_index - 1))
        else {
            return Err(IsotropyError::SubductionIndexOutOfRange {
                entry,
                index: irrep_index,
            });
        };
        result.push(DoubleValuedSubduction {
            parent_sg: *sg,
            parent_ml: label,
            frequency: crate::irrep::generated_data::ISOTROPY_W_SUBDUCE_FREQUENCY[entry] as u16,
        });
    }
    Ok(result)
}

/// One-line-per-entry rendering of the subduction of an isotropy subgroup.
pub fn format_identity_subduction(ordinal: usize) -> Result<String, IsotropyError> {
    let entries = identity_subduction(ordinal)?;
    let mut lines = Vec::new();
    lines.push(format!(
        "// {} parent irrep(s) subduce the trivial irrep of isotropy record {}",
        entries.len(),
        ordinal + 1
    ));
    lines.push("| Parent irrep | BC | k | dim | i(G) | Dir | domain |".to_string());
    lines.push("|--------------|----|---|-----|------|-----|--------|".to_string());
    for entry in entries {
        lines.push(format!(
            "| {} | {} | ({}/{}, {}/{}, {}/{}) | {} | {} | {} | {} |",
            entry.parent_ml,
            entry.parent_bc,
            entry.parent_k.numerators[0],
            entry.parent_k.denominator,
            entry.parent_k.numerators[1],
            entry.parent_k.denominator,
            entry.parent_k.numerators[2],
            entry.parent_k.denominator,
            entry.parent_dim,
            entry.frequency,
            entry.direction_label,
            entry.domain,
        ));
    }
    let double_valued = double_valued_subduction(ordinal)?;
    if !double_valued.is_empty() {
        lines.push(String::new());
        lines.push("| Double-valued parent irrep | i(G) |".to_string());
        lines.push("|----------------------------|------|".to_string());
        for entry in double_valued {
            lines.push(format!(
                "| {} (SG {}) | {} |",
                entry.parent_ml, entry.parent_sg, entry.frequency
            ));
        }
    }
    Ok(lines.join("\n"))
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
/// The stored basis is a primitive basis of the subgroup lattice in the
/// parent's primitive frame, so the volume ratio is simply `|det W|`.  This is
/// the "Size" column of Stokes & Hatch (1988) and of the ISOTROPY program.
pub fn subgroup_size(basis: [[i32; 3]; 3]) -> Result<u32, IsotropyError> {
    let det = det3(basis);
    if det == 0 {
        return Err(IsotropyError::SingularSubgroupBasis { basis });
    }
    u32::try_from(det.unsigned_abs())
        .map_err(|_| IsotropyError::SubgroupSizeOverflow { determinant: det })
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

/// Primitive basis vectors of a space group's lattice in its conventional basis.
///
/// Rows are the primitive vectors in units of the conventional basis vectors,
/// using the ISOTROPY/ITA conventions (`I`: `(-1/2,1/2,1/2) …`, `F`:
/// `(0,1/2,1/2) …`, `C`: `(1/2,1/2,0) …`, `R`: obverse rhombohedral vectors in
/// hexagonal axes).  Pinned against the `P1` rows printed by the bundled
/// ISOTROPY binary; see `scripts/verify_isotropy_oracle.py`.
pub fn parent_primitive_basis(sg: u8) -> Result<[[f64; 3]; 3], IsotropyError> {
    use crate::spg_database::Centering;
    if sg == 0 || sg > 230 {
        return Err(IsotropyError::InvalidSpaceGroup(sg));
    }
    let hall = crate::irrep::generated_data::SG_DATA_HALL[sg as usize] as usize;
    let centering = crate::spg_database::get_spacegroup_type(hall).centering;
    let basis = match centering {
        Centering::Primitive => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        Centering::Body => [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]],
        Centering::Face => [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
        Centering::AFace => [[0.0, 0.5, 0.5], [0.0, -0.5, 0.5], [1.0, 0.0, 0.0]],
        Centering::BFace => [[0.5, 0.0, 0.5], [-0.5, 0.0, 0.5], [0.0, 1.0, 0.0]],
        Centering::CFace => [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]],
        Centering::RCenter => [
            [2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [-1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [-1.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
        ],
        Centering::Error => return Err(IsotropyError::MissingCentering { sg }),
    };
    Ok(basis)
}

/// Express the stored subgroup basis in the parent's conventional basis.
///
/// The result is still a primitive basis of the subgroup lattice (not the
/// subgroup's ITA conventional cell), now in parent conventional coordinates.
pub fn basis_in_parent_conventional(
    parent_sg: u8,
    basis: [[i32; 3]; 3],
) -> Result<[[f64; 3]; 3], IsotropyError> {
    let primitive = parent_primitive_basis(parent_sg)?;
    let mut converted = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            converted[i][j] = (0..3).map(|t| basis[i][t] as f64 * primitive[t][j]).sum();
        }
    }
    Ok(converted)
}

/// Express the stored origin shift in the parent's conventional basis.
///
/// This is a pure **frame conversion** of the stored value: `w_conventional =
/// w_primitive · P`, with `P` the parent's primitive basis.  Values are only
/// defined modulo the parent's lattice.
///
/// # Not the printed "Origin" column in general
///
/// The bundled ISOTROPY program prints its own Origin column, and a table-wide
/// sweep (`scripts/verify_isotropy_oracle.py` runs a sample of it) shows that
/// this frame-converted value differs from the printed one by a
/// **record-dependent** offset for a sizeable minority of records — the offset
/// is neither a parent-lattice vector nor a per-space-group constant
/// (counterexample: SG 139 `M1-` direction `P1` stores `(2,2,2)` in the
/// primitive frame, i.e. the parent origin, while the program prints
/// `(1/4,1/4,1/4)`).  The convention that maps one to the other is **not
/// pinned**; do not present either value as the other.  Use
/// [`IsotropyRecord::origin_rational`] / [`IsotropyRecord::origin_shift`] when
/// the verbatim upstream data is what you need.
pub fn origin_shift_in_parent_conventional(
    parent_sg: u8,
    origin: [i32; 4],
) -> Result<[f64; 3], IsotropyError> {
    let primitive = parent_primitive_basis(parent_sg)?;
    let w = decode_origin_checked(origin)?;
    let mut converted = [0.0f64; 3];
    for j in 0..3 {
        converted[j] = (0..3).map(|t| w[t] * primitive[t][j]).sum();
    }
    Ok(converted)
}

/// Decode an origin shift, rejecting non-positive denominators.
///
/// Only a frame-free decode of the stored encoding; see
/// [`origin_shift_in_parent_conventional`] for the frame conversion and its
/// caveats.
///
/// The generated tables always satisfy this; the check keeps the public
/// conversion helpers total for caller-supplied values.
pub(crate) fn decode_origin_checked(origin: [i32; 4]) -> Result<[f64; 3], IsotropyError> {
    crate::irrep::types::decode_origin_checked(origin)
        .ok_or(IsotropyError::InvalidOrigin { origin })
}

/// Determinant of a row-major 3x3 integer matrix.
///
/// Computed in `i64` so that caller-supplied bases cannot overflow `i32`; the
/// generated tables only contain small entries.
fn det3(m: [[i32; 3]; 3]) -> i64 {
    let m = m.map(|row| row.map(i64::from));
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
        "| # | Subgroup | Label | Direction | Dir | Size | Basis | Origin* | Domains | Arms |"
            .to_string(),
    );
    lines.push(
        "|---|----------|-------|-----------|-----|------|-------|--------|---------|------|"
            .to_string(),
    );
    for (index, sub) in subgroups.iter().enumerate() {
        let record = &sub.record;
        let size = subgroup_size(record.basis)
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
            format_basis(basis_in_parent_conventional(sg, record.basis)?),
            format_origin(origin_shift_in_parent_conventional(sg, record.origin)?),
            record.domains,
            record.arms,
        ));
    }
    lines.push(
        "\n*Origin* is the stored primitive-frame shift converted with the parent's \
         primitive basis; it equals the ISOTROPY program's printed Origin only \
         when the stored setting already matches the program's default (see \
         `origin_shift_in_parent_conventional`)."
            .to_string(),
    );
    Ok(lines.join("\n"))
}

/// Markdown table of the magnetic isotropy subgroups of one irrep.
///
/// The basis and origin are converted to the parent's conventional basis; the
/// basis remains a primitive basis of the magnetic subgroup lattice.
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
    lines.push("| # | UNI | BNS | Label | Direction | Size | Basis | Origin |".to_string());
    lines.push("|---|-----|-----|-------|-----------|------|-------|--------|".to_string());
    for (index, sub) in subgroups.iter().enumerate() {
        let record = &sub.record;
        let size = subgroup_size(record.basis)
            .map(|value| value.to_string())
            .unwrap_or_else(|_| "?".to_string());
        lines.push(format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} |",
            index,
            record.mag_sg,
            record.bns_label,
            record.iso_label,
            record.direction,
            size,
            format_basis(basis_in_parent_conventional(sg, record.basis)?),
            format_origin(origin_shift_in_parent_conventional(sg, record.origin)?),
        ));
    }
    Ok(lines.join("\n"))
}

fn format_basis(basis: [[f64; 3]; 3]) -> String {
    format!(
        "({},{},{}),({},{},{}),({},{},{})",
        format_component(basis[0][0]),
        format_component(basis[0][1]),
        format_component(basis[0][2]),
        format_component(basis[1][0]),
        format_component(basis[1][1]),
        format_component(basis[1][2]),
        format_component(basis[2][0]),
        format_component(basis[2][1]),
        format_component(basis[2][2]),
    )
}

fn format_component(value: f64) -> String {
    if (value - value.round()).abs() < 1e-9 {
        return format!("{}", value.round() as i64);
    }
    // Crystalline basis components are multiples of 1/2, 1/3 or 1/6; print the
    // exact fraction rather than a rounded decimal.
    for denominator in [2u32, 3, 4, 6] {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < 1e-9 {
            return format!("{}/{}", scaled.round() as i64, denominator);
        }
    }
    format!("{value:.6}")
}

fn format_origin(origin: [f64; 3]) -> String {
    format!(
        "({},{},{})",
        format_component(origin[0]),
        format_component(origin[1]),
        format_component(origin[2])
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
        assert_eq!(subgroup_size(sub.record.basis).unwrap(), 1);
    }

    #[test]
    fn p222_r1_reaches_centred_f222_with_primitive_basis() {
        let sub = isotropy_subgroup_for_direction(16, "R1", IsotropyDirection::Index(0))
            .expect("R1 has an isotropy subgroup");
        assert_eq!(sub.record.sg, 22);
        assert_eq!(sub.record.basis, [[0, 1, 1], [1, 0, 1], [1, 1, 0]]);
        assert_eq!(sub.record.origin, [0, 0, 0, 1]);
        // The ISO program prints Size 2 for this transition.
        assert_eq!(subgroup_size(sub.record.basis).unwrap(), 2);
        // P222 is primitive, so both frames coincide.
        assert_eq!(
            basis_in_parent_conventional(16, sub.record.basis).unwrap(),
            [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
        );
    }

    #[test]
    fn rhombohedral_origin_matches_the_conventional_frame() {
        // #167 R-3c, Γ3+ along (a,0) → #2 P-1 with a doubled cell.  The stored
        // primitive-frame origin (0, 1/2, 0) is (-1/6, 1/6, 1/6) in hexagonal
        // conventional coordinates, which is what the ISOTROPY program prints.
        let sub = isotropy_subgroup_for_direction(167, "GM3+", IsotropyDirection::Label("P1"))
            .expect("Γ3+ has a P1 direction");
        assert_eq!(sub.record.origin, [0, 1, 0, 2]);
        assert_eq!(sub.record.origin_shift(), Some([0.0, 0.5, 0.0]));
        let converted =
            origin_shift_in_parent_conventional(167, sub.record.origin).expect("conversion");
        let expected = [-1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];
        for (got, want) in converted.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-12, "{converted:?} != {expected:?}");
        }
    }

    #[test]
    fn centred_parent_origin_conversion_uses_the_primitive_basis() {
        // #229 Im-3m, Γ4- along (a,0,0) → #107 I4mm: same lattice, so the size
        // is 1 and the origin is the parent origin in both frames.
        let sub =
            isotropy_subgroup_for_direction(229, "GM4-", IsotropyDirection::Descriptor("(a,0,0)"))
                .expect("Γ4- has an (a,0,0) direction");
        assert_eq!(sub.record.sg, 107);
        assert_eq!(subgroup_size(sub.record.basis).unwrap(), 1);
        assert_eq!(
            origin_shift_in_parent_conventional(229, sub.record.origin).unwrap(),
            [0.0, 0.0, 0.0]
        );
    }

    #[test]
    fn parent_primitive_basis_matches_iso_conventions() {
        assert_eq!(
            parent_primitive_basis(1).unwrap(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(
            parent_primitive_basis(225).unwrap(),
            [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
        );
        assert_eq!(
            parent_primitive_basis(229).unwrap(),
            [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]
        );
        let r = parent_primitive_basis(167).unwrap();
        assert!((r[0][0] - 2.0 / 3.0).abs() < 1e-12);
        assert!((r[2][1] + 2.0 / 3.0).abs() < 1e-12);
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

    /// Expected values are the `SHOW FREQ DIR` output of the bundled
    /// ISOTROPY binary, converted to `(Miller-Love label, frequency, domain)`.
    #[test]
    fn identity_subduction_matches_the_isotropy_program() {
        let cases = [
            // (SG, irrep, direction, [(label, i(G), domain, Dir)])
            (
                221,
                "GM4+",
                "(a,0,0)",
                vec![
                    ("GM1+", 1u16, 1u16, "P1"),
                    ("GM3+", 1, 3, "P1"),
                    ("GM4+", 1, 1, "P1"),
                ],
            ),
            (
                221,
                "GM3+",
                "(a,0)",
                vec![("GM1+", 1, 1, "P1"), ("GM3+", 1, 1, "P1")],
            ),
            (
                16,
                "R1",
                "(a)",
                vec![("GM1", 1, 1, "P1"), ("R1", 1, 1, "P1")],
            ),
            (
                225,
                "GM4-",
                "(a,0,0)",
                vec![
                    ("GM1+", 1, 1, "P1"),
                    ("GM3+", 1, 3, "P1"),
                    ("GM4-", 1, 1, "P1"),
                ],
            ),
        ];
        for (sg, ml, direction, expected) in cases {
            let subgroup =
                isotropy_subgroup_for_direction(sg, ml, IsotropyDirection::Descriptor(direction))
                    .unwrap_or_else(|error| panic!("SG {sg} {ml} {direction}: {error}"));
            let entries = subgroup
                .identity_subduction()
                .unwrap_or_else(|error| panic!("SG {sg} {ml} {direction}: {error}"));
            let actual: Vec<(&str, u16, u16, &str)> = entries
                .iter()
                .map(|entry| {
                    (
                        entry.parent_ml,
                        entry.frequency,
                        entry.domain,
                        entry.direction_label,
                    )
                })
                .collect();
            assert_eq!(actual, expected, "SG {sg} {ml} {direction}");
            assert!(entries.iter().all(|entry| entry.parent_sg == sg));
        }
    }

    #[test]
    fn subduction_ranges_are_well_formed() {
        let ranges = &crate::irrep::generated_data::ISOTROPY_SUBDUCE_RANGES;
        assert_eq!(
            ranges.len(),
            crate::irrep::generated_data::ISOTROPY_SUBGROUPS.len() + 1
        );
        assert!(ranges.windows(2).all(|pair| pair[0] <= pair[1]));
        assert_eq!(
            ranges[ranges.len() - 1] as usize,
            crate::irrep::generated_data::ISOTROPY_SUBDUCE_IRREP.len() + 1
        );
        assert!(matches!(
            identity_subduction(ranges.len()),
            Err(IsotropyError::InvalidIsotropyRecord { .. })
        ));
    }
}
