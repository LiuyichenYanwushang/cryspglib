//! The supported entry point for ordinary irrep subduction (milestone R8).
//!
//! One request type, one computation entry point, and a report that carries its
//! own provenance.  The engine below it (`subduction::star::decompose`) stays
//! as it is; this module only *resolves* a request into exactly one
//! `(parent irrep, isotropy subgroup, embedding)` context and then reports what
//! the engine returned, with the identity of every input attached.
//!
//! Design rules this module enforces (each one is a permanent test):
//!
//! * **One entry point.** [`subduce`] is the only place that turns names into
//!   an engine call.  [`subduction_for_direction`] is a thin spelling of it.
//! * **Explicit label convention.** Every name lookup takes a
//!   [`LabelConvention`]; a BC label that matches more than one record is
//!   refused as [`SubductionApiError::AmbiguousCondensing`] instead of taking
//!   the first hit.
//! * **No silent fallback to identity-only.** If the engine cannot build the
//!   full decomposition the call returns an error.  The identity multiplicity
//!   keeps its own separate entry points
//!   ([`trivial_content`](crate::irrep::subduction::star::decompose::trivial_content_with_embedding)
//!   and its line counterpart); it is never substituted for a decomposition.
//! * **Parameterized k is explicit.** A request may carry an exact
//!   [`Rat`] parameter along the direction's own line; it is only accepted for
//!   the frozen line sources, and the report says whether the result is a line
//!   irrep or a formal induction ([`ParameterKind`]).
//!
//! ```
//! use cryspglib::irrep::subduction_api::{SubductionRequest, subduce};
//! use cryspglib::irrep::{LabelConvention, isotropy::IsotropyDirection};
//!
//! let request = SubductionRequest::new(221, "GM4+", LabelConvention::Cdml)
//!     .direction(IsotropyDirection::Label("P1"));
//! let report = subduce(&request).expect("SG 221 GM4+ has a P1 isotropy subgroup");
//! assert_eq!(report.condensing().cdml, "GM4+");
//! assert_eq!(report.direction().label, "P1");
//! assert!(report.covered_dimension() > 0);
//! assert!(!report.blocks().is_empty());
//! ```

use crate::irrep::isotropy::{IsotropyDirection, IsotropyError, IsotropySubgroup};
use crate::irrep::labels::LabelConvention;
use crate::irrep::line_monodromy::line_table;
use crate::irrep::query;
use crate::irrep::subduction::star::decompose::{
    FullStarBlock, FullStarError, FullStarSubduction, FullStarTarget, LineSubduction, ParameterKind,
    subduce_full_star_with_embedding, subduce_line_at_parameter,
};
use crate::irrep::subduction::{Rat, SubductionComponent, SubductionError, SubgroupEmbedding};
use crate::irrep::types::IrrepRecord;
use crate::irrep::{isotropy, w_little_characters_data::LittleCharacterTable};

/// The parameter classification carried by a parameterized report.
pub use crate::irrep::subduction::star::decompose::ParameterKind as ReportedParameterKind;

/// Why a [`SubductionRequest`] could not be answered.
///
/// Every variant carries the input that failed, so a caller never has to guess
/// which name was ambiguous or which direction did not exist.
#[derive(Debug, Clone, PartialEq)]
pub enum SubductionApiError {
    /// The parent space group number is outside 1–230.
    UnsupportedSpaceGroup { sg: u8 },
    /// No record of that parent carries the requested label under the requested
    /// convention.
    CondensingNotFound {
        sg: u8,
        label: String,
        convention: LabelConvention,
    },
    /// The requested label matches more than one record.  BC spellings are not
    /// unique in the generated tables, so this is a real input error, not a
    /// defensive branch; the candidates are the CDML labels of the matches.
    AmbiguousCondensing {
        sg: u8,
        label: String,
        convention: LabelConvention,
        candidates: Vec<&'static str>,
    },
    /// A parameter was requested for a parent irrep that is not one of the
    /// frozen parametric-k line sources.
    NotAParameterizedLine {
        sg: u8,
        cdml: &'static str,
    },
    /// A parameter was requested through a condensing irrep and direction
    /// whose subgroup carries several frozen line sources, so the request does
    /// not name one line.  Name the source label instead.
    AmbiguousLineSource {
        sg: u8,
        cdml: &'static str,
        candidates: Vec<&'static str>,
    },
    /// The isotropy lookup failed (unknown label/descriptor/index).
    Isotropy(IsotropyError),
    /// The embedding or the decomposition itself failed.
    ///
    /// Boxed: the engine errors are large structs and this is a cold path; the
    /// wrapper keeps `SubductionApiError` small enough to return by value.
    Subduction(Box<SubductionError>),
    /// The full-star decomposition failed (boxed for the same reason).
    FullStar(Box<FullStarError>),
}

impl std::fmt::Display for SubductionApiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSpaceGroup { sg } => {
                write!(f, "space group {sg} is outside 1-230")
            }
            Self::CondensingNotFound {
                sg,
                label,
                convention,
            } => write!(
                f,
                "space group {sg} has no {convention:?} label {label:?}"
            ),
            Self::AmbiguousCondensing {
                sg,
                label,
                convention,
                candidates,
            } => write!(
                f,
                "space group {sg}: {convention:?} label {label:?} matches {} records ({}); \
                 select one by its CDML label",
                candidates.len(),
                candidates.join(", ")
            ),
            Self::NotAParameterizedLine { sg, cdml } => write!(
                f,
                "space group {sg} {cdml} is not one of the frozen parametric-k line sources, \
                 so it takes no parameter"
            ),
            Self::AmbiguousLineSource {
                sg,
                cdml,
                candidates,
            } => write!(
                f,
                "space group {sg} {cdml}: this direction carries {} frozen line sources ({}); \
                 name the source label to select one",
                candidates.len(),
                candidates.join(", ")
            ),
            Self::Isotropy(error) => write!(f, "{error}"),
            Self::Subduction(error) => write!(f, "{error}"),
            Self::FullStar(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for SubductionApiError {}

impl From<IsotropyError> for SubductionApiError {
    fn from(error: IsotropyError) -> Self {
        Self::Isotropy(error)
    }
}

impl From<SubductionError> for SubductionApiError {
    fn from(error: SubductionError) -> Self {
        Self::Subduction(Box::new(error))
    }
}

impl From<FullStarError> for SubductionApiError {
    fn from(error: FullStarError) -> Self {
        Self::FullStar(Box::new(error))
    }
}

/// One subduction request: which parent irrep, along which direction, at which
/// exact parameter.
///
/// `parameter` is `None` for the direction's own (pinned) k point.  `Some(t)`
/// asks for the frozen parametric-k line at `t` and is only valid for the
/// sources of [`crate::irrep::line_monodromy::line_table`].
#[derive(Debug, Clone, Copy)]
pub struct SubductionRequest<'a> {
    parent_sg: u8,
    condensing: &'a str,
    convention: LabelConvention,
    direction: IsotropyDirection<'a>,
    parameter: Option<Rat>,
}

impl<'a> SubductionRequest<'a> {
    /// A request for the irrep `condensing` of `parent_sg`, read under
    /// `convention`, along the direction `direction`.
    pub fn new(parent_sg: u8, condensing: &'a str, convention: LabelConvention) -> Self {
        Self {
            parent_sg,
            condensing,
            convention,
            direction: IsotropyDirection::Index(0),
            parameter: None,
        }
    }

    /// Select the isotropy subgroup by ISOTROPY label, component descriptor or
    /// table index.
    pub fn direction(mut self, direction: IsotropyDirection<'a>) -> Self {
        self.direction = direction;
        self
    }

    /// Ask for the exact parameter `t` on this direction's frozen line.
    pub fn parameter(mut self, parameter: Rat) -> Self {
        self.parameter = Some(parameter);
        self
    }

    /// The parent space group number.
    pub fn parent_sg(&self) -> u8 {
        self.parent_sg
    }

    /// The condensing label as requested.
    pub fn condensing_label(&self) -> &'a str {
        self.condensing
    }

    /// The requested label convention.
    pub fn convention(&self) -> LabelConvention {
        self.convention
    }

    /// The requested direction.
    pub fn requested_direction(&self) -> IsotropyDirection<'a> {
        self.direction
    }

    /// The requested exact parameter, if any.
    pub fn requested_parameter(&self) -> Option<Rat> {
        self.parameter
    }
}

/// The condensing irrep's own identity, both label spellings included.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CondensingIdentity {
    /// Parent space group number.
    pub parent_sg: u8,
    /// CDML / Miller–Love label.
    pub cdml: &'static str,
    /// Bradley–Cracknell label, if the generated table carries one.
    pub bc: Option<String>,
    /// Kovalev label as stored.
    pub kov: &'static str,
    /// Dimension of the condensing irrep.
    pub dimension: u8,
    /// Whether the condensing irrep is double valued (spinor).
    pub spinor: bool,
    /// Stored k vector numerator triple and common denominator of the
    /// condensing irrep, exactly as the generated table holds them.
    pub k: [[i8; 3]; 2],
}

/// The isotropy subgroup's table provenance: everything needed to re-find the
/// same record without going through the same lookup again.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DirectionIdentity {
    /// ISOTROPY direction label (`"P1"`, `"4D1"`, …).
    pub label: &'static str,
    /// Component description as stored in the table.
    pub descriptor: &'static str,
    /// Dimension of the order-parameter space.
    pub dimension: u8,
    /// Number of free parameters along the direction.
    pub free: u8,
    /// Position of this record in the parent irrep's isotropy table.
    pub index_in_irrep: usize,
    /// Subgroup space group number.
    pub subgroup_sg: u8,
    /// Subgroup Hermann–Mauguin symbol.
    pub subgroup_symbol: &'static str,
    /// Subgroup Schoenflies symbol.
    pub subgroup_schoenflies: &'static str,
    /// Primitive basis of the subgroup in the parent's primitive frame.
    pub basis: [[i32; 3]; 3],
    /// Origin of the subgroup setting, `[x, y, z, d]` meaning `(x/d, y/d, z/d)`.
    pub origin: [i32; 4],
    /// Number of domains.
    pub domains: usize,
    /// Number of arms in the star.
    pub arms: usize,
    /// Position of this record in the generated isotropy table.
    pub ordinal: usize,
}

/// One target of one block, with both label spellings.
#[derive(Debug, Clone, PartialEq)]
pub struct TargetView {
    /// Subgroup space group number of the target.
    pub sg: u8,
    /// CDML label of the target, when the frozen table names one.
    pub cdml: Option<&'static str>,
    /// Bradley–Cracknell label of the target, when available.
    pub bc: Option<&'static str>,
    /// The physical child row this term came from (display identity).
    pub row_ml: Option<&'static str>,
    /// Complex dimension of the child little-group irrep.
    pub dimension: u8,
    /// Multiplicity in the child little-group representation at the folded point.
    pub multiplicity: u32,
    /// Frozen CIR source number, when the target is a pinned row.
    pub irnumber: Option<u32>,
    /// Which complex constituent of the source row this term is.
    pub component: SubductionComponent,
}

/// One block of the full-star decomposition (one folded child star).
#[derive(Debug, Clone, PartialEq)]
pub struct BlockView {
    /// The child point this block was matched against, exactly.
    pub stored_k: [Rat; 3],
    /// Size of the child star.
    pub star_size: usize,
    /// Number of parent arms in the block.
    pub arm_count: usize,
    /// The parent arm indices carried by the block (its stable identity).
    pub arm_indices: Vec<usize>,
    /// Dimension of the block.
    pub dimension: u32,
    /// Dimension of the child little-group irrep.
    pub little_dimension: u32,
    /// The block's targets.
    pub targets: Vec<TargetView>,
}

/// The complete answer to one [`SubductionRequest`].
#[derive(Debug, Clone, PartialEq)]
pub struct SubductionReport {
    condensing: CondensingIdentity,
    direction: DirectionIdentity,
    parameter: Option<Rat>,
    parameter_kind: Option<ParameterKind>,
    blocks: Vec<BlockView>,
    covered_dimension: u32,
}

impl SubductionReport {
    /// The condensing irrep's identity (both label spellings).
    pub fn condensing(&self) -> &CondensingIdentity {
        &self.condensing
    }

    /// The isotropy subgroup's table provenance.
    pub fn direction(&self) -> &DirectionIdentity {
        &self.direction
    }

    /// The exact parameter this report is for, if the request had one.
    pub fn parameter(&self) -> Option<Rat> {
        self.parameter
    }

    /// Whether the result is a real line irrep or a formal induction from the
    /// frozen line little group.  `None` for a non-parameterized request.
    pub fn parameter_kind(&self) -> Option<ParameterKind> {
        self.parameter_kind
    }

    /// The full-star blocks, in the engine's order.
    pub fn blocks(&self) -> &[BlockView] {
        &self.blocks
    }

    /// Dimension covered by the decomposition (sum over blocks of the block
    /// dimension times its multiplicity).
    pub fn covered_dimension(&self) -> u32 {
        self.covered_dimension
    }
}

/// Resolve the condensing irrep under an explicit convention.
fn resolve_condensing(
    parent_sg: u8,
    label: &str,
    convention: LabelConvention,
) -> Result<&'static IrrepRecord, SubductionApiError> {
    if !(1..=230).contains(&parent_sg) {
        return Err(SubductionApiError::UnsupportedSpaceGroup { sg: parent_sg });
    }
    let mut matches = query::find_irreps(parent_sg, label, convention);
    match matches.len() {
        0 => Err(SubductionApiError::CondensingNotFound {
            sg: parent_sg,
            label: label.to_string(),
            convention,
        }),
        1 => Ok(matches.pop().expect("exactly one match")),
        _ => {
            let mut candidates: Vec<&'static str> =
                matches.iter().map(|record| record.ml).collect();
            candidates.sort_unstable();
            candidates.dedup();
            Err(SubductionApiError::AmbiguousCondensing {
                sg: parent_sg,
                label: label.to_string(),
                convention,
                candidates,
            })
        }
    }
}

/// The position of this subgroup in its irrep's own isotropy list, which is
/// what [`IsotropyDirection::Index`] takes (the generated table's `ordinal` is
/// global and much larger).
fn index_in_irrep(subgroup: &IsotropySubgroup) -> Result<usize, SubductionApiError> {
    let list = isotropy::isotropy_subgroups(
        subgroup.parent_sg,
        subgroup.irrep_ml,
        LabelConvention::Cdml,
    )?;
    list.iter()
        .position(|candidate| candidate.ordinal == subgroup.ordinal)
        .ok_or(SubductionApiError::Isotropy(IsotropyError::SubgroupIndexOutOfRange {
            sg: subgroup.parent_sg,
            ml: subgroup.irrep_ml.to_string(),
            index: subgroup.ordinal,
            len: list.len(),
        }))
}

/// Find the context of a frozen parametric-k line source.
///
/// The frozen sources are labelled by their own source label (`"DT1"`), which
/// is not necessarily the CDML label of the condensing irrep: the rows of a
/// subgroup carry the frozen label, while the subgroup itself is built from
/// the irrep record.  The scan walks the parent's irreps (never the other way
/// round) and keeps the first subgroup whose rows name this source.
fn resolve_line_source(
    parent_sg: u8,
    label: &str,
    direction: IsotropyDirection<'_>,
) -> Result<(&'static IrrepRecord, IsotropySubgroup), SubductionApiError> {
    for record in query::irreps_of(parent_sg) {
        if record.spinor || record.subgroups().is_empty() {
            continue;
        }
        let Ok(subgroup) =
            isotropy::isotropy_subgroup_for_direction(parent_sg, record.ml, LabelConvention::Cdml, direction)
        else {
            continue;
        };
        let Ok(rows) = subgroup.other_wave_vector_subduction() else {
            continue;
        };
        if rows.iter().any(|row| row.parent_ml == label) {
            return Ok((record, subgroup));
        }
    }
    Err(SubductionApiError::CondensingNotFound {
        sg: parent_sg,
        label: label.to_string(),
        convention: LabelConvention::Cdml,
    })
}

fn direction_identity(
    subgroup: &IsotropySubgroup,
) -> Result<DirectionIdentity, SubductionApiError> {
    let record = &subgroup.record;
    Ok(DirectionIdentity {
        label: record.direction_label,
        descriptor: record.direction,
        dimension: record.direction_dim,
        free: record.direction_free,
        index_in_irrep: index_in_irrep(subgroup)?,
        subgroup_sg: record.sg as u8,
        subgroup_symbol: record.symbol,
        subgroup_schoenflies: record.schoenflies,
        basis: record.basis,
        origin: record.origin,
        domains: record.domains,
        arms: record.arms,
        ordinal: subgroup.ordinal,
    })
}

fn target_view(target: &FullStarTarget) -> TargetView {
    TargetView {
        sg: target.sg,
        cdml: target.ml,
        bc: target.bc,
        row_ml: target.row_ml,
        dimension: target.dimension,
        multiplicity: target.multiplicity,
        irnumber: target.irnumber,
        component: target.component,
    }
}

fn block_view(block: &FullStarBlock) -> BlockView {
    BlockView {
        stored_k: [
            block.stored_k().get(0),
            block.stored_k().get(1),
            block.stored_k().get(2),
        ],
        star_size: block.star_size(),
        arm_count: block.arm_count(),
        arm_indices: block.arm_indices(),
        dimension: block.block_dimension(),
        little_dimension: block.little_dimension(),
        targets: block.targets().iter().map(target_view).collect(),
    }
}

/// Answer one [`SubductionRequest`].
///
/// The request is resolved into exactly one isotropic subgroup and one
/// embedding; the full-star decomposition is then computed by the production
/// engine.  With a parameter the frozen line source is used instead, and the
/// returned [`ParameterKind`] says whether that parameter is a genuine line
/// irrep or a formal induction from the frozen line little group.
///
/// Errors are returned, never replaced by an identity-only answer: see
/// [`SubductionApiError`].
pub fn subduce(request: &SubductionRequest<'_>) -> Result<SubductionReport, SubductionApiError> {
    // A name is first read as an irrep label; if that fails, as a frozen
    // parametric-k **line source** label, whose spelling (`"DT1"`) need not be
    // the CDML label of the condensing irrep it belongs to.
    let (record, subgroup, line_source) =
        match resolve_condensing(request.parent_sg, request.condensing, request.convention) {
            Ok(record) => {
                let subgroup = isotropy::isotropy_subgroup_for_direction(
                    request.parent_sg,
                    record.ml,
                    LabelConvention::Cdml,
                    request.direction,
                )?;
                (record, subgroup, line_table(request.parent_sg, record.ml))
            }
            Err(SubductionApiError::CondensingNotFound { .. }) => {
                let (record, subgroup) = resolve_line_source(
                    request.parent_sg,
                    request.condensing,
                    request.direction,
                )?;
                let table = line_table(request.parent_sg, request.condensing);
                (record, subgroup, table)
            }
            Err(other) => return Err(other),
        };
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&subgroup)?;
    let condensing = CondensingIdentity {
        parent_sg: record.sg,
        cdml: record.ml,
        bc: record.label(LabelConvention::Bc),
        kov: record.kov,
        dimension: record.dim,
        spinor: record.spinor,
        k: [
            [record.kx, record.ky, record.kz],
            [record.kd, record.kd, record.kd],
        ],
    };
    let direction = direction_identity(&subgroup)?;
    let mut report = SubductionReport {
        condensing,
        direction,
        parameter: request.parameter,
        parameter_kind: None,
        blocks: Vec::new(),
        covered_dimension: 0,
    };
    match request.parameter {
        None => {
            let result: FullStarSubduction =
                subduce_full_star_with_embedding(&subgroup, &embedding, record)?;
            report.covered_dimension = result.covered_dimension();
            report.blocks = result.blocks().iter().map(block_view).collect();
        }
        Some(parameter) => {
            // A direction can carry several frozen lines through the same
            // subgroup (SG 196 `W1` `4D1` carries `DT1`, `DT2` and `SM1`), so
            // naming the irrep is not enough unless exactly one of them is a
            // frozen source.
            let table: &'static LittleCharacterTable = match line_source {
                Some(table) => table,
                None => {
                    let mut candidates: Vec<&'static str> = subgroup
                        .other_wave_vector_subduction()?
                        .iter()
                        .map(|row| row.parent_ml)
                        .filter(|label| line_table(request.parent_sg, label).is_some())
                        .collect();
                    candidates.sort_unstable();
                    candidates.dedup();
                    match candidates.len() {
                        0 => {
                            return Err(SubductionApiError::NotAParameterizedLine {
                                sg: request.parent_sg,
                                cdml: record.ml,
                            })
                        }
                        1 => line_table(request.parent_sg, candidates[0])
                            .expect("the candidate came from the table"),
                        _ => {
                            return Err(SubductionApiError::AmbiguousLineSource {
                                sg: request.parent_sg,
                                cdml: record.ml,
                                candidates,
                            })
                        }
                    }
                }
            };
            let result: LineSubduction =
                subduce_line_at_parameter(&subgroup, &embedding, table, parameter)?;
            report.parameter_kind = Some(result.parameter_kind());
            report.covered_dimension = result.covered_dimension();
            report.blocks = result.blocks().iter().map(block_view).collect();
        }
    }
    Ok(report)
}

/// Shorthand for `subduce(&SubductionRequest::new(sg, condensing, convention)
/// .direction(direction))`.
pub fn subduction_for_direction(
    parent_sg: u8,
    condensing: &str,
    convention: LabelConvention,
    direction: IsotropyDirection<'_>,
) -> Result<SubductionReport, SubductionApiError> {
    subduce(&SubductionRequest::new(parent_sg, condensing, convention).direction(direction))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_gamma_point_lookup_is_literal_under_cdml() {
        let record = resolve_condensing(221, "GM4+", LabelConvention::Cdml).expect("GM4+ exists");
        assert_eq!(record.ml, "GM4+");
        assert_eq!(record.dim, 3);
        // The Γ point stays "GM" in CDML; the BC spelling is a different string.
        assert!(resolve_condensing(221, "Γ4+", LabelConvention::Cdml).is_err());
        let bc = resolve_condensing(221, "Γ4+", LabelConvention::Bc).expect("BC spelling");
        assert_eq!(bc.ml, record.ml);
    }

    #[test]
    fn out_of_range_space_groups_are_rejected() {
        assert_eq!(
            resolve_condensing(0, "GM1", LabelConvention::Cdml).unwrap_err(),
            SubductionApiError::UnsupportedSpaceGroup { sg: 0 }
        );
        assert_eq!(
            resolve_condensing(231, "GM1", LabelConvention::Cdml).unwrap_err(),
            SubductionApiError::UnsupportedSpaceGroup { sg: 231 }
        );
    }
}
