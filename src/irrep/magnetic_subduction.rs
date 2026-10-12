//! R26: the first restricted scalar Γ magnetic subduction.
//!
//! [`subduce_scalar_gamma`] consumes the verified R25
//! [`ScalarGammaModel`](crate::irrep::magnetic_input::ScalarGammaModel) — one
//! parent's real `±1` Γ character restricted to an embedded magnetic group `M` —
//! and identifies the **Γ magnetic co-representation of `M`** that the restriction
//! is equivalent to, reconstructs that target's operators on the one-dimensional
//! space, and exhibits the complex intertwiner between the two.
//!
//! The target identity is **independent** of the model: it is read from
//! [`magnetic_summary::magnetic_irrep_summary_by_uni`]'s Γ coreps, whose unitary
//! characters come from the magnetic group's unitary subgroup `H` and its ordinary
//! PIR rows, not from the parent character row the model used.  The summary's
//! antiunitary Type-A columns carry **no** value at all
//! ([`crate::irrep::corep::CharacterCompleteness::TypeAAntiunitaryPending`]); this module never reads
//! them, never uses an antiunitary trace as a character inner product or as
//! operator evidence, and never treats a `None` there as "the target does not
//! exist".  The antiunitary operators are **constructed** from the target's
//! unitary characters and the actual multiplication of `M`.
//!
//! # What is computed
//!
//! For the summary's Γ k-point (the little group of Γ is all of `M`), every listed
//! operation is mapped to a class of the model *by rotation, time-reversal flag and
//! translation modulo the model's own translation lattice*; a class that no listed
//! operation covers, or an operation that no class covers, is a named
//! [`MagneticSubductionError::TargetOperationUnmatched`] /
//! [`MagneticSubductionError::TargetClassUncovered`] refusal — the correspondence
//! is **not** guessed by a setting search.
//!
//! For each candidate — a Γ corep of dimension one with exactly one unitary,
//! non-spinor, one-dimensional source irrep `K` — the canonical one-dimensional
//! co-representation is built:
//!
//! ```text
//! D(g) = K(g)              for a unitary class g
//! D(a) = K(a₀⁻¹·a)         for an antiunitary class a, with a fixed canonical a₀
//! ```
//!
//! the *choice* of `a₀` (the first antiunitary class in the model's class order)
//! makes the antiunitary representative act as `1·K` instead of introducing a free
//! phase.  The construction is then checked against the actual multiplication of
//! `M`: the square `D(a₀)² = D(a₀²)` (equivalently `K(a₀²) = +1`, the one-
//! dimensional case of Wigner's type A), conjugation compatibility
//! `D(h)·D(a₀) = D(h·a₀)` for every unitary class `h`, and the full ordered-pair
//! corepresentation law.  A candidate whose unitary part is not equivalent to the
//! model's restriction is skipped; a candidate whose unitary part *does* match but
//! whose own square or law fails is refused by name instead of being replaced by
//! another target.
//!
//! # The intertwiner and the multiplicity
//!
//! The source restriction `A` (from the model) and the target `D` are equivalent
//! iff there is a complex `U ≠ 0` with
//!
//! ```text
//! D(g) = U⁻¹ · A(g) · U            for a unitary class g
//! D(a) = U⁻¹ · A(a) · conj(U)      for an antiunitary class a
//! ```
//!
//! Both sides are one-dimensional and `A`, `D` take the values `±1` there, so this
//! is an exact two-unknown real-linear system for `U = u₁ + i·u₂`.  The module
//! solves it, keeps the target iff a non-zero solution exists, verifies the
//! exhibited `U` class by class, and reports the multiplicity from
//! `dim(source) = multiplicity · dim(target)`, which for this dimension-one family
//! is exactly what the equivalence proves.  `η = +1` and `η = −1` are a gauge of
//! the antiunitary coset: they give the same target and different intertwiners
//! (`U = 1` versus `U = i`), and both are legal inputs.
//!
//! # Scope fences
//!
//! The source is the R25 scalar Γ model (real `±1`, non-spinor, one-dimensional,
//! `k = 0`).  The target must be a one-dimensional Wigner **type A** Γ
//! co-representation of the model's magnetic group with a single non-spinor
//! one-dimensional unitary source — this includes a purely unitary (Type I
//! ordinary) group as the trivial case with no antiunitary class, and the
//! one-dimensional antiunitary classes of Type III/IV groups.  Complex scalar,
//! non-Γ, spinor, higher-dimensional, full-star and setting-search cases are
//! outside this entry point; they either cannot reach it (R25 refuses them) or have
//! no candidate of this shape (the type-B/type-C targets of dimension two and the
//! spinor rows are not reported as this module's question).  The
//! [`MagneticSubduction::scope`] string is the same fence in one line.

use num_complex::Complex64;

use crate::irrep::corep::CorepType;
use crate::irrep::magnetic_embedding::{MagneticOperation, canonical_translation, compose};
use crate::irrep::magnetic_input::ScalarGammaModel;
use crate::irrep::magnetic_summary::{self, MagneticCorepSummary, MagneticKPointSummary};
use crate::irrep::query;
use crate::irrep::subduction::{Mat3I, Rat};
use crate::irrep::types::IrrepRecord;

/// The Γ (k = 0) k-point coordinates the summary tables use.
const GAMMA: (i8, i8, i8, i8) = (0, 0, 0, 1);

/// Denominator grid a summary translation is snapped to, mirroring the exact
/// translation contract of [`crate::irrep::magnetic_embedding`].
const MAX_TRANSLATION_DENOMINATOR: i128 = 96;

/// Tolerance of the exact-grid snap of a summary `f64`.  The value is snapped and
/// the exact rational is what is compared, so an off-grid value is an error rather
/// than a rounded translation.
const GRID_TOLERANCE: f64 = 1e-9;

/// Why a scalar Γ magnetic subduction was refused.
///
/// Every variant is a specific refusal: an out-of-family candidate and a genuine
/// data gap never share a name, and a failed correspondence is never reported as
/// "the target does not exist".
#[derive(Debug)]
pub enum MagneticSubductionError {
    /// The model's own source label no longer resolves to one scalar Γ irrep.
    SourceIdentityMismatch { parent_sg: u8, label: &'static str },
    /// The magnetic corep summary of the model's UNI could not be computed.
    SummaryFailed { uni: usize, reason: String },
    /// The summary has no Γ k-point.
    NoGammaKPoint { uni: usize },
    /// A Γ source corep that could be one-dimensional could not be constructed,
    /// so the candidate list is not known to be complete.
    UnresolvedTargetData {
        uni: usize,
        source: String,
        minimum_dimension: Option<usize>,
        reason: String,
    },
    /// A summary Γ operation has no class in the model's magnetic group under the
    /// exact correspondence (rotation, time-reversal flag, translation modulo the
    /// model's translation lattice).  The two tables are generated in their own
    /// settings; this round has not established an alignment between them, so the
    /// request is refused instead of being paired by a guessed transform.  The
    /// message claims no more than that: it does not claim no alignment exists.
    TargetOperationUnmatched { uni: usize, rotation: Mat3I, time_reversal: bool },
    /// A class of the model's magnetic group is not covered by the summary's Γ
    /// operation list.
    TargetClassUncovered { uni: usize, class: usize },
    /// A summary Γ translation is outside the exact grid this contract represents.
    TargetTranslationOffGrid { uni: usize, column: usize, value: [f64; 3] },
    /// A unitary Γ column has no character value, so the candidate cannot be
    /// examined and is not silently dropped.
    TargetCharacterUnavailable { uni: usize, label: String, column: usize },
    /// Two columns of one unitary class disagree on the character.
    TargetCharacterNotAClassFunction { uni: usize, class: usize },
    /// The candidate's own square is incompatible with a one-dimensional
    /// co-representation: `K(a₀²) ≠ +1` for the canonical antiunitary `a₀`.
    TargetSquareInconsistent { uni: usize, label: String },
    /// The candidate's antiunitary coset does not act bijectively on the unitary
    /// classes through the canonical `a₀`.
    TargetCosetInconsistent { uni: usize, label: String },
    /// Conjugation compatibility failed for one unitary class.
    TargetConjugationInconsistent { uni: usize, label: String, class: usize },
    /// The reconstructed operator table is not a corepresentation for one ordered
    /// pair of classes.
    TargetLawViolated {
        uni: usize,
        label: String,
        left: usize,
        right: usize,
        expected: TargetOperator,
        found: TargetOperator,
    },
    /// No one-dimensional Γ target is equivalent to the model's restriction.  The
    /// labels examined are carried so a missing target is not confused with a
    /// target that was never considered.
    NoMatchingTarget {
        parent_sg: u8,
        source: &'static str,
        uni: usize,
        examined: Vec<String>,
    },
    /// Several distinct Γ targets are equivalent to the restriction.
    AmbiguousTarget { uni: usize, labels: Vec<String> },
    /// The dimension identity `dim(source) = multiplicity · dim(target)` failed.
    DimensionNotConserved {
        source: u32,
        target: u32,
        multiplicity: u32,
    },
    /// Exact arithmetic inside the checker failed.
    Arithmetic { reason: String },
}

impl std::fmt::Display for MagneticSubductionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SourceIdentityMismatch { parent_sg, label } => write!(
                f,
                "SG {parent_sg} no longer has exactly one scalar Gamma irrep labelled {label:?}"
            ),
            Self::SummaryFailed { uni, reason } => {
                write!(f, "UNI {uni}: the magnetic corep summary failed: {reason}")
            }
            Self::NoGammaKPoint { uni } => {
                write!(f, "UNI {uni}: the magnetic corep summary has no Gamma k-point")
            }
            Self::UnresolvedTargetData {
                uni,
                source,
                minimum_dimension,
                reason,
            } => write!(
                f,
                "UNI {uni}: the Gamma corep of source {source:?} is unresolved \
                 (minimum dimension {minimum_dimension:?}): {reason}"
            ),
            Self::TargetOperationUnmatched {
                uni,
                rotation,
                time_reversal,
            } => write!(
                f,
                "UNI {uni}: the summary's Gamma operation with rotation {rotation:?} and \
                 time reversal {time_reversal} has no class in the magnetic group; the table \
                 correspondence is not established and was not searched for"
            ),
            Self::TargetClassUncovered { uni, class } => write!(
                f,
                "UNI {uni}: class {class} of the magnetic group is not covered by the \
                 summary's Gamma operation list"
            ),
            Self::TargetTranslationOffGrid { uni, column, value } => write!(
                f,
                "UNI {uni}: column {column} has translation {value:?}, off the exact grid \
                 (denominators up to {MAX_TRANSLATION_DENOMINATOR})"
            ),
            Self::TargetCharacterUnavailable { uni, label, column } => write!(
                f,
                "UNI {uni}: target {label:?} has no unitary character at column {column}"
            ),
            Self::TargetCharacterNotAClassFunction { uni, class } => write!(
                f,
                "UNI {uni}: the columns of unitary class {class} disagree on the character"
            ),
            Self::TargetSquareInconsistent { uni, label } => write!(
                f,
                "UNI {uni}: target {label:?} fails the one-dimensional square K(a0^2) = +1"
            ),
            Self::TargetCosetInconsistent { uni, label } => write!(
                f,
                "UNI {uni}: target {label:?} does not act bijectively on the unitary classes \
                 from its canonical antiunitary a0"
            ),
            Self::TargetConjugationInconsistent { uni, label, class } => write!(
                f,
                "UNI {uni}: target {label:?} fails conjugation compatibility at unitary class \
                 {class}"
            ),
            Self::TargetLawViolated {
                uni,
                label,
                left,
                right,
                expected,
                found,
            } => write!(
                f,
                "UNI {uni}: target {label:?} violates the corepresentation law at classes \
                 ({left}, {right}): expected {expected:?}, found {found:?}"
            ),
            Self::NoMatchingTarget {
                parent_sg,
                source,
                uni,
                examined,
            } => write!(
                f,
                "SG {parent_sg} source {source:?} restricted to UNI {uni} is equivalent to none \
                 of the examined one-dimensional Gamma targets {examined:?}"
            ),
            Self::AmbiguousTarget { uni, labels } => write!(
                f,
                "UNI {uni}: several distinct Gamma targets are equivalent to the restriction: \
                 {labels:?}"
            ),
            Self::DimensionNotConserved {
                source,
                target,
                multiplicity,
            } => write!(
                f,
                "dimension is not conserved: source {source} != multiplicity {multiplicity} * \
                 target {target}"
            ),
            Self::Arithmetic { reason } => write!(f, "exact arithmetic failed: {reason}"),
        }
    }
}

impl std::error::Error for MagneticSubductionError {}

/// One reconstructed operator on the one-dimensional target space.
///
/// A unitary operator acts as `z ↦ value · z`; an antiunitary one as
/// `z ↦ value · conj(z)`.  Composition is map composition, which is exactly the
/// corepresentation product: the right factor is conjugated precisely when the
/// left factor is antiunitary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TargetOperator {
    value: Complex64,
    antiunitary: bool,
}

impl TargetOperator {
    /// The scalar factor of the operator.
    pub fn value(&self) -> Complex64 {
        self.value
    }

    /// Whether the operator is antiunitary.
    pub fn is_antiunitary(&self) -> bool {
        self.antiunitary
    }

    /// Apply the operator to one coordinate.
    pub fn apply(&self, z: Complex64) -> Complex64 {
        if self.antiunitary {
            self.value * z.conj()
        } else {
            self.value * z
        }
    }

    /// The composition `self ∘ inner` as an operator.
    pub fn compose(&self, inner: &TargetOperator) -> TargetOperator {
        let factor = if self.antiunitary {
            inner.value.conj()
        } else {
            inner.value
        };
        TargetOperator {
            value: self.value * factor,
            antiunitary: self.antiunitary ^ inner.antiunitary,
        }
    }
}

/// The result of one scalar Γ magnetic subduction.
///
/// The fields are private and the only construction path is
/// [`subduce_scalar_gamma`], so a successful object cannot be forged from outside.
/// The report **borrows** the verified model instead of copying it: the embedding,
/// both identities and the `±1` action table stay the model's own data, and nothing
/// can be mutated behind the report.
#[derive(Debug)]
pub struct MagneticSubduction<'a> {
    model: &'a ScalarGammaModel,
    target_label: String,
    target_source_label: &'static str,
    target_unitary_sg: u8,
    target_dimension: u32,
    multiplicity: u32,
    source_dimension: u32,
    intertwiner: Complex64,
    operators: Vec<TargetOperator>,
}

impl<'a> MagneticSubduction<'a> {
    /// The verified input model this report is about.
    pub fn model(&self) -> &'a ScalarGammaModel {
        self.model
    }

    /// The target corep's label as the magnetic summary names it.
    pub fn target_label(&self) -> &str {
        &self.target_label
    }

    /// The unitary subgroup's ordinary irrep label the target is induced from.
    pub fn target_source_label(&self) -> &'static str {
        self.target_source_label
    }

    /// Space group number of the target's unitary subgroup `H`.
    pub fn target_unitary_sg(&self) -> u8 {
        self.target_unitary_sg
    }

    /// Dimension of the target corepresentation (one in this family).
    pub fn target_dimension(&self) -> u32 {
        self.target_dimension
    }

    /// Multiplicity of the target in the restriction.
    pub fn multiplicity(&self) -> u32 {
        self.multiplicity
    }

    /// Dimension of the subduced source corepresentation, recomputed from the
    /// source irrep's own record.
    pub fn source_dimension(&self) -> u32 {
        self.source_dimension
    }

    /// `multiplicity · target_dimension`, checked equal to `source_dimension`.
    pub fn covered_dimension(&self) -> u32 {
        self.multiplicity * self.target_dimension
    }

    /// The complex intertwiner `U` with `D = U⁻¹·A·U` on unitary classes and
    /// `D = U⁻¹·A·conj(U)` on antiunitary ones.
    pub fn intertwiner(&self) -> Complex64 {
        self.intertwiner
    }

    /// The reconstructed target operator of every class of `M`, in the model's
    /// class order.
    pub fn operators(&self) -> &[TargetOperator] {
        &self.operators
    }

    /// The reconstructed target operator of one class.
    pub fn operator_at(&self, index: usize) -> Option<&TargetOperator> {
        self.operators.get(index)
    }

    /// The classes of `M` the operators refer to, in the model's class order.
    pub fn classes(&self) -> &'a [MagneticOperation] {
        self.model.classes()
    }

    /// The model's translation lattice, the one the correspondence reduced modulo.
    pub fn translation_lattice(&self) -> &'a [[Rat; 3]] {
        self.model.translation_lattice()
    }

    /// One line stating what this entry point covers.
    pub fn scope(&self) -> &'static str {
        SCOPE
    }
}

/// The scope fence of [`subduce_scalar_gamma`] in one line.
pub const SCOPE: &str = "One-dimensional Wigner type-A Gamma corepresentation of an R25 scalar \
Gamma magnetic model (Type I ordinary included, single non-spinor one-dimensional unitary source); \
complex scalar, non-Gamma, spinor, higher-dimensional, full-star and setting-search cases are \
outside this entry point.";

/// Identify the Γ magnetic co-representation of the model's magnetic group that the
/// model's restriction is equivalent to, and reconstruct its operators and the
/// complex intertwiner.
pub fn subduce_scalar_gamma(
    model: &ScalarGammaModel,
) -> Result<MagneticSubduction<'_>, MagneticSubductionError> {
    let parent_sg = model.parent_sg();
    let source = resolve_source(model)?;
    let uni = model.uni();
    let partial = magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni)
        .map_err(|error| MagneticSubductionError::SummaryFailed {
            uni,
            reason: format!("{error:?}"),
        })?;
    let gamma = partial
        .summary
        .kpoints
        .iter()
        .find(|kpoint| kpoint.coords == GAMMA)
        .ok_or(MagneticSubductionError::NoGammaKPoint { uni })?;

    // A possibly one-dimensional Gamma source whose corep could not be built means
    // the candidate list is not known to be complete; that is a named data refusal,
    // never a silently missing target.
    for unresolved in &partial.unresolved_coreps {
        if unresolved.k_label == gamma.label
            && !unresolved.spinor
            && unresolved.minimum_dimension.is_none_or(|dim| dim <= 1)
        {
            return Err(MagneticSubductionError::UnresolvedTargetData {
                uni,
                source: unresolved.source_irrep.clone(),
                minimum_dimension: unresolved.minimum_dimension,
                reason: unresolved.reason.clone(),
            });
        }
    }

    let classes = model.classes();
    let lattice = model.translation_lattice();
    let class_columns = map_operations(gamma, classes, lattice, uni)?;
    let first_anti = classes.iter().position(|class| class.time_reversal);
    let source_actions: Vec<SourceAction> = model
        .actions()
        .iter()
        .map(|action| (action.value.value(), action.class.time_reversal))
        .collect();

    let mut examined: Vec<String> = Vec::new();
    let mut matches: Vec<Candidate> = Vec::new();
    for corep in &gamma.coreps {
        if !is_one_dimensional_single_source(corep) {
            continue;
        }
        examined.push(corep.label.clone());
        let Some(unitary) = unitary_signs(corep, &class_columns, classes, uni)? else {
            // A unitary character that is not the real ±1 the source family takes
            // cannot be equivalent to this source; it is examined, not omitted.
            continue;
        };
        if !unitary_part_matches(&unitary, model) {
            continue;
        }
        let operators = build_target(corep, classes, lattice, &unitary, first_anti, uni)?;
        check_square(corep, classes, lattice, &operators, &unitary, uni)?;
        check_conjugation(corep, classes, lattice, &operators, uni)?;
        check_law(corep, classes, lattice, &operators, uni)?;
        // A candidate without a non-zero intertwiner is simply not this source's
        // target; it is skipped.  No phase is ever substituted for a missing
        // equivalence.
        let Some(intertwiner) = solve_intertwiner(&operators, &source_actions) else {
            continue;
        };
        matches.push(Candidate {
            label: corep.label.clone(),
            source_label: corep.source_irreps[0].ml,
            dimension: corep.dim as u32,
            operators,
            intertwiner,
        });
    }

    // Two examined candidates with the same canonical operator table are the same
    // corepresentation under two spellings; keeping the first and not duplicating it
    // is the "complex-equivalent targets are not counted twice" rule.  Two distinct
    // canonical tables that are both equivalent to the restriction cannot happen for
    // a one-dimensional source; the check names it instead of picking one.
    let mut distinct: Vec<Candidate> = Vec::new();
    for candidate in matches {
        if distinct
            .iter()
            .any(|kept| kept.operators == candidate.operators)
        {
            continue;
        }
        distinct.push(candidate);
    }
    let [target] = distinct.as_slice() else {
        return match distinct.len() {
            0 => Err(MagneticSubductionError::NoMatchingTarget {
                parent_sg,
                source: model.source_label(),
                uni,
                examined,
            }),
            _ => Err(MagneticSubductionError::AmbiguousTarget {
                uni,
                labels: distinct.iter().map(|one| one.label.clone()).collect(),
            }),
        };
    };

    let target_dimension = target.dimension;
    let source_dimension = u32::from(source.dim);
    let multiplicity = source_dimension.checked_div(target_dimension).filter(|multiplicity| {
        *multiplicity > 0 && multiplicity * target_dimension == source_dimension
    });
    let Some(multiplicity) = multiplicity else {
        return Err(MagneticSubductionError::DimensionNotConserved {
            source: source_dimension,
            target: target_dimension,
            multiplicity: 1,
        });
    };
    let target_unitary_sg = partial.summary.unitary_sg;

    Ok(MagneticSubduction {
        model,
        target_label: target.label.clone(),
        target_source_label: target.source_label,
        target_unitary_sg,
        target_dimension,
        multiplicity,
        source_dimension,
        intertwiner: target.intertwiner,
        operators: target.operators.clone(),
    })
}

/// One examined candidate that turned out equivalent to the restriction.
struct Candidate {
    label: String,
    source_label: &'static str,
    dimension: u32,
    operators: Vec<TargetOperator>,
    intertwiner: Complex64,
}

/// The source identity, re-resolved from the database: exactly one scalar,
/// non-spinor, one-dimensional Γ irrep must carry the model's canonical label.
fn resolve_source(model: &ScalarGammaModel) -> Result<&'static IrrepRecord, MagneticSubductionError> {
    let label = model.source_label();
    let matches = query::find_irreps(model.parent_sg(), label, model.source_convention());
    match matches.as_slice() {
        [irrep]
            if irrep.dim == 1
                && !irrep.spinor
                && irrep.k_vector().numerators == [0, 0, 0]
                && irrep.k_vector().denominator == 1 =>
        {
            Ok(irrep)
        }
        _ => Err(MagneticSubductionError::SourceIdentityMismatch {
            parent_sg: model.parent_sg(),
            label,
        }),
    }
}

/// Whether a summary corep is a candidate of this module's shape.
fn is_one_dimensional_single_source(corep: &MagneticCorepSummary) -> bool {
    corep.dim == 1
        && corep.corep_type == CorepType::A
        && corep.source_irreps.len() == 1
        && corep.source_irreps[0].dim == 1
        && !corep.source_irreps[0].spinor
}

/// Map every summary Γ operation to a class of the model, exactly.
fn map_operations(
    gamma: &MagneticKPointSummary,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
    uni: usize,
) -> Result<Vec<Vec<usize>>, MagneticSubductionError> {
    let mut class_columns = vec![Vec::new(); classes.len()];
    for operation in &gamma.operations {
        let translation = [
            exact_fraction(operation.translation[0]).ok_or(
                MagneticSubductionError::TargetTranslationOffGrid {
                    uni,
                    column: operation.column,
                    value: operation.translation,
                },
            )?,
            exact_fraction(operation.translation[1]).ok_or(
                MagneticSubductionError::TargetTranslationOffGrid {
                    uni,
                    column: operation.column,
                    value: operation.translation,
                },
            )?,
            exact_fraction(operation.translation[2]).ok_or(
                MagneticSubductionError::TargetTranslationOffGrid {
                    uni,
                    column: operation.column,
                    value: operation.translation,
                },
            )?,
        ];
        let mapped = MagneticOperation {
            rotation: operation.rotation,
            translation,
            time_reversal: operation.time_reversal,
        };
        let index = class_of(&mapped, classes, lattice).ok_or(
            MagneticSubductionError::TargetOperationUnmatched {
                uni,
                rotation: operation.rotation,
                time_reversal: operation.time_reversal,
            },
        )?;
        class_columns[index].push(operation.column);
    }
    for (index, columns) in class_columns.iter().enumerate() {
        if columns.is_empty() {
            return Err(MagneticSubductionError::TargetClassUncovered { uni, class: index });
        }
    }
    Ok(class_columns)
}

/// The sign of every unitary class, or `None` when a character is not the real
/// `±1` this family takes (such a candidate cannot be equivalent to the source).
fn unitary_signs(
    corep: &MagneticCorepSummary,
    class_columns: &[Vec<usize>],
    classes: &[MagneticOperation],
    uni: usize,
) -> Result<Option<Vec<i32>>, MagneticSubductionError> {
    let mut signs = vec![0i32; classes.len()];
    for (index, class) in classes.iter().enumerate() {
        if class.time_reversal {
            continue;
        }
        let mut value = None;
        for &column in &class_columns[index] {
            let character = corep
                .characters
                .get(column)
                .copied()
                .flatten()
                .ok_or(MagneticSubductionError::TargetCharacterUnavailable {
                    uni,
                    label: corep.label.clone(),
                    column,
                })?;
            match value {
                None => value = Some(character),
                Some(previous) if previous == character => {}
                Some(_) => {
                    return Err(MagneticSubductionError::TargetCharacterNotAClassFunction {
                        uni,
                        class: index,
                    })
                }
            }
        }
        let Some(character) = value else {
            return Err(MagneticSubductionError::TargetCharacterUnavailable {
                uni,
                label: corep.label.clone(),
                column: class_columns[index][0],
            });
        };
        match real_sign(character) {
            Some(sign) => signs[index] = sign,
            None => return Ok(None),
        }
    }
    Ok(Some(signs))
}

/// `±1` when the complex value is on the exact real grid, `None` otherwise.
fn real_sign(value: Complex64) -> Option<i32> {
    if value.im.abs() > GRID_TOLERANCE {
        return None;
    }
    for candidate in [1.0f64, -1.0] {
        if (value.re - candidate).abs() < GRID_TOLERANCE {
            return Some(candidate as i32);
        }
    }
    None
}

/// Whether the target's unitary signs equal the model's restriction on every
/// unitary class.
fn unitary_part_matches(signs: &[i32], model: &ScalarGammaModel) -> bool {
    model.actions().iter().enumerate().all(|(index, action)| {
        action.class.time_reversal || action.value.value() == signs[index]
    })
}

/// Build the canonical one-dimensional target table from the unitary signs and the
/// actual multiplication of `M`.
///
/// `first_anti` is the antiunitary representative the canonical gauge is taken
/// from.  The caller passes the first antiunitary class of the model's order; a
/// *different* antiunitary representative of the same coset is a legal gauge and
/// produces the same corepresentation up to a phase on the coset, which the tests
/// exercise through this parameter.
fn build_target(
    corep: &MagneticCorepSummary,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
    unitary: &[i32],
    first_anti: Option<usize>,
    uni: usize,
) -> Result<Vec<TargetOperator>, MagneticSubductionError> {
    let mut operators = vec![
        TargetOperator {
            value: Complex64::new(0.0, 0.0),
            antiunitary: false,
        };
        classes.len()
    ];
    for (index, class) in classes.iter().enumerate() {
        if !class.time_reversal {
            operators[index] = TargetOperator {
                value: Complex64::new(f64::from(unitary[index]), 0.0),
                antiunitary: false,
            };
        }
    }
    let Some(first_anti) = first_anti else {
        return Ok(operators);
    };
    if !classes[first_anti].time_reversal {
        return Err(MagneticSubductionError::TargetCosetInconsistent {
            uni,
            label: corep.label.clone(),
        });
    }
    let square = class_of(
        &compose(&classes[first_anti], &classes[first_anti]).map_err(|error| {
            MagneticSubductionError::Arithmetic {
                reason: error.to_string(),
            }
        })?,
        classes,
        lattice,
    )
    .ok_or_else(|| MagneticSubductionError::TargetCosetInconsistent {
        uni,
        label: corep.label.clone(),
    })?;
    if classes[square].time_reversal {
        return Err(MagneticSubductionError::TargetCosetInconsistent {
            uni,
            label: corep.label.clone(),
        });
    }
    for (index, class) in classes.iter().enumerate() {
        if class.time_reversal {
            continue;
        }
        let image = product_class(&classes[first_anti], class, classes, lattice).ok_or_else(|| {
            MagneticSubductionError::TargetCosetInconsistent {
                uni,
                label: corep.label.clone(),
            }
        })?;
        if !classes[image].time_reversal {
            return Err(MagneticSubductionError::TargetCosetInconsistent {
                uni,
                label: corep.label.clone(),
            });
        }
        operators[image] = TargetOperator {
            value: Complex64::new(f64::from(unitary[index]), 0.0),
            antiunitary: true,
        };
    }
    if operators
        .iter()
        .enumerate()
        .any(|(index, operator)| classes[index].time_reversal && operator.value == Complex64::new(0.0, 0.0))
    {
        return Err(MagneticSubductionError::TargetCosetInconsistent {
            uni,
            label: corep.label.clone(),
        });
    }
    Ok(operators)
}

/// The one-dimensional square check: the canonical antiunitary representative must
/// square to a unitary class whose character is `+1`.
fn check_square(
    corep: &MagneticCorepSummary,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
    operators: &[TargetOperator],
    unitary: &[i32],
    uni: usize,
) -> Result<(), MagneticSubductionError> {
    let Some(first_anti) = classes.iter().position(|class| class.time_reversal) else {
        return Ok(());
    };
    let square = product_class(&classes[first_anti], &classes[first_anti], classes, lattice)
        .ok_or_else(|| MagneticSubductionError::TargetCosetInconsistent {
            uni,
            label: corep.label.clone(),
        })?;
    let composite = operators[first_anti].compose(&operators[first_anti]);
    if unitary[square] != 1 || composite != operators[square] {
        return Err(MagneticSubductionError::TargetSquareInconsistent {
            uni,
            label: corep.label.clone(),
        });
    }
    Ok(())
}

/// `D(h)·D(a₀) = D(h·a₀)` and `D(a₀)·conj(D(h)) = D(a₀·h)` for every unitary class.
fn check_conjugation(
    corep: &MagneticCorepSummary,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
    operators: &[TargetOperator],
    uni: usize,
) -> Result<(), MagneticSubductionError> {
    let Some(first_anti) = classes.iter().position(|class| class.time_reversal) else {
        return Ok(());
    };
    for (index, class) in classes.iter().enumerate() {
        if class.time_reversal {
            continue;
        }
        let left = product_class(class, &classes[first_anti], classes, lattice).ok_or_else(|| {
            MagneticSubductionError::TargetCosetInconsistent {
                uni,
                label: corep.label.clone(),
            }
        })?;
        let right = product_class(&classes[first_anti], class, classes, lattice).ok_or_else(|| {
            MagneticSubductionError::TargetCosetInconsistent {
                uni,
                label: corep.label.clone(),
            }
        })?;
        if operators[index].compose(&operators[first_anti]) != operators[left]
            || operators[first_anti].compose(&operators[index]) != operators[right]
        {
            return Err(MagneticSubductionError::TargetConjugationInconsistent {
                uni,
                label: corep.label.clone(),
                class: index,
            });
        }
    }
    Ok(())
}

/// The full ordered-pair corepresentation law on the class set.
fn check_law(
    corep: &MagneticCorepSummary,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
    operators: &[TargetOperator],
    uni: usize,
) -> Result<(), MagneticSubductionError> {
    for (left, left_class) in classes.iter().enumerate() {
        for (right, right_class) in classes.iter().enumerate() {
            let product = product_class(left_class, right_class, classes, lattice).ok_or_else(|| {
                MagneticSubductionError::TargetCosetInconsistent {
                    uni,
                    label: corep.label.clone(),
                }
            })?;
            let expected = operators[product];
            let found = operators[left].compose(&operators[right]);
            if expected != found {
                return Err(MagneticSubductionError::TargetLawViolated {
                    uni,
                    label: corep.label.clone(),
                    left,
                    right,
                    expected,
                    found,
                });
            }
        }
    }
    Ok(())
}

/// The source `±1` action of one class: its value and its time-reversal flag.
type SourceAction = (i32, bool);

/// Solve `D = U⁻¹·A·U` / `D = U⁻¹·A·conj(U)` for the one-dimensional
/// intertwiner, then **verify the exhibited `U` on both branches class by class**.
///
/// Both sides are one-dimensional and of unit modulus, so every condition is a
/// real-linear row in the two real unknowns `U = u₁ + i·u₂` (`A` and `D` are
/// complex in general and `±1` for this module's own canonical tables).  The
/// accumulated solution set is the span of the first row's null direction, and a
/// later row keeps that span exactly when it is **orthogonal** to the direction; a
/// row that is not makes the system rank two, so there is no non-zero `U` and the
/// candidate is simply not this source's target.  Testing the row for parallelism
/// with the direction instead would reject a repeated, perfectly consistent
/// constraint as if it were a contradiction.
///
/// The returned phase is the solved direction, normalized and sign-fixed, which on
/// the real `±1` tables is exactly `U = 1` or `U = i` and on a legal complex gauge
/// of the coset is that gauge's own phase.  `None` is the only failure answer and
/// the caller never turns it into a phase: the final loop re-applies the contract
/// to every class.
fn solve_intertwiner(target: &[TargetOperator], source: &[SourceAction]) -> Option<Complex64> {
    if target.len() != source.len() {
        return None;
    }
    let mut rows: Vec<[f64; 2]> = Vec::new();
    for (index, &(a, time_reversal)) in source.iter().enumerate() {
        if target[index].antiunitary != time_reversal {
            // The table's unitarity branch has to be the class's own.
            return None;
        }
        let source_value = Complex64::new(f64::from(a), 0.0);
        let target_value = target[index].value;
        if time_reversal {
            // a·conj(U) − U·d = (a − d)·u₁ − i·(a + d)·u₂ = 0
            let difference = source_value - target_value;
            let sum = source_value + target_value;
            rows.push([difference.re, sum.im]);
            rows.push([difference.im, -sum.re]);
        } else {
            // (a − d)·U = 0
            let difference = source_value - target_value;
            rows.push([difference.re, -difference.im]);
            rows.push([difference.im, difference.re]);
        }
    }
    let mut direction: Option<[f64; 2]> = None;
    for row in rows {
        let row_norm = (row[0] * row[0] + row[1] * row[1]).sqrt();
        if row_norm < GRID_TOLERANCE {
            continue;
        }
        match direction {
            None => direction = Some([-row[1], row[0]]),
            Some(known) => {
                let known_norm = (known[0] * known[0] + known[1] * known[1]).sqrt();
                let orthogonal = row[0] * known[0] + row[1] * known[1];
                if orthogonal.abs() > GRID_TOLERANCE * row_norm * known_norm {
                    return None;
                }
            }
        }
    }
    let intertwiner = match direction {
        None => Complex64::new(1.0, 0.0),
        Some([u1, u2]) => {
            let value = Complex64::new(u1, u2);
            let norm = value.norm();
            if norm < GRID_TOLERANCE {
                return None;
            }
            // Canonical sign: the positive real direction, or the positive
            // imaginary one when the direction is purely imaginary.
            let sign = if u1 > 0.0 || (u1 == 0.0 && u2 >= 0.0) {
                1.0
            } else {
                -1.0
            };
            value * sign / norm
        }
    };
    if intertwiner.norm() < 0.5 {
        return None;
    }
    // The exhibited U is re-verified on every class of both branches.  The answer
    // is never taken from the solver's own confidence.
    let inverse = intertwiner.conj() / intertwiner.norm_sqr();
    for (index, &(a, time_reversal)) in source.iter().enumerate() {
        let source_value = Complex64::new(f64::from(a), 0.0);
        let factor = if time_reversal {
            intertwiner.conj()
        } else {
            intertwiner
        };
        let transformed = inverse * source_value * factor;
        if (transformed - target[index].value).norm() > GRID_TOLERANCE {
            return None;
        }
    }
    Some(intertwiner)
}

/// The class an operation denotes, or `None` when it is not a class.
fn class_of(
    operation: &MagneticOperation,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
) -> Option<usize> {
    let canonical = canonical_translation(operation.translation, lattice);
    classes.iter().position(|class| {
        class.rotation == operation.rotation
            && class.time_reversal == operation.time_reversal
            && canonical_translation(class.translation, lattice) == canonical
    })
}

/// The class of `left · right`, or `None` when the product or its class is missing.
fn product_class(
    left: &MagneticOperation,
    right: &MagneticOperation,
    classes: &[MagneticOperation],
    lattice: &[[Rat; 3]],
) -> Option<usize> {
    compose(left, right)
        .ok()
        .and_then(|product| class_of(&product, classes, lattice))
}

/// Exact conversion of one summary fraction with a small denominator.
fn exact_fraction(value: f64) -> Option<Rat> {
    for denominator in 1..=MAX_TRANSLATION_DENOMINATOR {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < GRID_TOLERANCE {
            return Rat::new(scaled.round() as i128, denominator).ok();
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::LabelConvention;
    use crate::irrep::magnetic_input::{
        EmbeddingEvidence, LatticePolicy, RealSign, RecordIdentity, ScalarGammaRequest,
        SourceIdentity, verify_scalar_gamma_input,
    };
    use crate::irrep::subduction::Mat3R;

    fn operator(value: f64, antiunitary: bool) -> TargetOperator {
        TargetOperator {
            value: Complex64::new(value, 0.0),
            antiunitary,
        }
    }

    /// The pinned R25 model of a witness, built through the public constructor.
    fn witness_model(
        sg: u8,
        source: &'static str,
        condensing: &'static str,
        index: usize,
        eta: RealSign,
    ) -> crate::irrep::magnetic_input::ScalarGammaModel {
        verify_scalar_gamma_input(
            &ScalarGammaRequest::new(
                SourceIdentity {
                    parent_sg: sg,
                    label: source,
                    convention: LabelConvention::Cdml,
                },
                RecordIdentity {
                    parent_sg: sg,
                    condensing_label: condensing,
                    convention: LabelConvention::Cdml,
                    direction_index: index,
                },
                Some(eta),
            )
            .with_embedding(EmbeddingEvidence {
                map: Mat3R::identity(),
                shift: [Rat::ZERO; 3],
                lattice_policy: LatticePolicy::ContainmentOnly,
            }),
        )
        .expect("the pinned witness model must verify")
    }

    fn source_actions(model: &ScalarGammaModel) -> Vec<SourceAction> {
        model
            .actions()
            .iter()
            .map(|action| (action.value.value(), action.class.time_reversal))
            .collect()
    }

    /// The candidate whose unitary part matches the model, with its unitary signs.
    fn matching_candidate(model: &ScalarGammaModel) -> (Vec<i32>, MagneticCorepSummary) {
        let uni = model.uni();
        let partial = magnetic_summary::magnetic_irrep_summary_by_uni_partial(uni).unwrap();
        let gamma = partial
            .summary
            .kpoints
            .iter()
            .find(|kpoint| kpoint.coords == GAMMA)
            .unwrap();
        let classes = model.classes();
        let lattice = model.translation_lattice();
        let columns = map_operations(gamma, classes, lattice, uni).unwrap();
        for corep in &gamma.coreps {
            if !is_one_dimensional_single_source(corep) {
                continue;
            }
            let Some(unitary) = unitary_signs(corep, &columns, classes, uni).unwrap() else {
                continue;
            };
            if unitary_part_matches(&unitary, model) {
                return (unitary, corep.clone());
            }
        }
        panic!("the pinned witness must have a matching candidate");
    }

    /// The intertwiner contract on actual complex vectors, both branches.
    fn check_vector_contract(target: &[TargetOperator], source: &[SourceAction], u: Complex64) {
        assert!(u.norm() > 1e-12, "the intertwiner must not be zero");
        for z in [Complex64::new(1.0, 2.0), Complex64::new(-3.0, 0.5)] {
            for (index, &(a, antiunitary)) in source.iter().enumerate() {
                assert_eq!(target[index].antiunitary, antiunitary);
                let value = Complex64::new(f64::from(a), 0.0);
                let left = value * if antiunitary { (u * z).conj() } else { u * z };
                let right = u * target[index].apply(z);
                assert!((left - right).norm() < 1e-12, "class {index}");
            }
        }
    }

    #[test]
    fn repeated_identical_antiunitary_constraints_keep_the_direction() {
        // Two antiunitary classes repeat A = -1 against a target acting as +1: the
        // repeated, consistent constraint must keep the imaginary direction.
        let target = [operator(1.0, false), operator(1.0, true), operator(1.0, true)];
        let source = [(1, false), (-1, true), (-1, true)];
        let u = solve_intertwiner(&target, &source).expect("a repeated constraint is not a contradiction");
        assert_eq!(u, Complex64::new(0.0, 1.0));
        check_vector_contract(&target, &source, u);
    }

    #[test]
    fn contradictory_antiunitary_constraints_have_no_solution() {
        // U would have to be real on one antiunitary class and imaginary on the
        // other at once.
        let target = [operator(1.0, false), operator(1.0, true), operator(1.0, true)];
        let source = [(1, false), (1, true), (-1, true)];
        assert_eq!(solve_intertwiner(&target, &source), None);
    }

    #[test]
    fn a_wrong_branch_or_character_is_not_an_intertwiner() {
        // The table's unitarity branch must be the class's own.
        let target = [operator(1.0, false), operator(1.0, false)];
        assert_eq!(solve_intertwiner(&target, &[(1, false), (1, true)]), None);
        // A unitary class with a different value cannot be conjugated away.
        assert_eq!(solve_intertwiner(&target, &[(1, false), (-1, false)]), None);
    }

    #[test]
    fn a_legal_complex_gauge_of_the_coset_is_still_equivalent() {
        // SG 10 GM1+ / GM2+ row 0 is UNI 53.  Multiplying the whole antiunitary
        // coset by i is a gauge of the same corepresentation, not new data: the
        // table stays a corepresentation and the restriction stays equivalent.
        let model = witness_model(10, "GM1+", "GM2+", 0, RealSign::Minus);
        let (unitary, corep) = matching_candidate(&model);
        let classes = model.classes();
        let lattice = model.translation_lattice();
        let uni = model.uni();
        let first_anti = classes.iter().position(|class| class.time_reversal);
        let canonical =
            build_target(&corep, classes, lattice, &unitary, first_anti, uni).unwrap();
        let source = source_actions(&model);
        let canonical_u = solve_intertwiner(&canonical, &source).expect("the canonical target");
        check_vector_contract(&canonical, &source, canonical_u);

        let gauge = Complex64::new(0.0, 1.0);
        let mut gauged = canonical.clone();
        for (index, class) in classes.iter().enumerate() {
            if class.time_reversal {
                gauged[index] = TargetOperator {
                    value: gauge * canonical[index].value(),
                    antiunitary: true,
                };
            }
        }
        assert_ne!(gauged, canonical);
        check_square(&corep, classes, lattice, &gauged, &unitary, uni).unwrap();
        check_conjugation(&corep, classes, lattice, &gauged, uni).unwrap();
        check_law(&corep, classes, lattice, &gauged, uni).unwrap();
        let gauged_u = solve_intertwiner(&gauged, &source).expect("a legal gauge is still equivalent");
        check_vector_contract(&gauged, &source, gauged_u);
    }

    #[test]
    fn an_alternative_antiunitary_representative_is_the_same_target_up_to_gauge() {
        // SG 12 GM2+ / Y1+ row 0 is UNI 56: four antiunitary classes and a unitary
        // character that takes -1, so choosing a different antiunitary a0 flips the
        // coset phase.  Every legal representative must stay a corepresentation of
        // the same M with the same unitary part and stay equivalent to the source.
        let model = witness_model(12, "GM2+", "Y1+", 0, RealSign::Minus);
        let (unitary, corep) = matching_candidate(&model);
        let classes = model.classes();
        let lattice = model.translation_lattice();
        let uni = model.uni();
        let source = source_actions(&model);
        let antiunitary: Vec<usize> = classes
            .iter()
            .enumerate()
            .filter(|(_, class)| class.time_reversal)
            .map(|(index, _)| index)
            .collect();
        assert!(antiunitary.len() >= 2);
        let canonical =
            build_target(&corep, classes, lattice, &unitary, Some(antiunitary[0]), uni).unwrap();
        let canonical_u = solve_intertwiner(&canonical, &source).expect("the canonical target");
        check_vector_contract(&canonical, &source, canonical_u);

        let mut saw_a_flipped_coset = false;
        for &representative in &antiunitary[1..] {
            let alternative =
                build_target(&corep, classes, lattice, &unitary, Some(representative), uni)
                    .unwrap();
            for (index, class) in classes.iter().enumerate() {
                if !class.time_reversal {
                    assert_eq!(alternative[index], canonical[index]);
                }
            }
            let ratio = alternative[representative].value() / canonical[representative].value();
            assert!((ratio.norm() - 1.0).abs() < 1e-12);
            for &index in &antiunitary {
                assert!(
                    (alternative[index].value() - ratio * canonical[index].value()).norm() < 1e-12,
                    "a representative change may only rephase the coset"
                );
            }
            check_square(&corep, classes, lattice, &alternative, &unitary, uni).unwrap();
            check_conjugation(&corep, classes, lattice, &alternative, uni).unwrap();
            check_law(&corep, classes, lattice, &alternative, uni).unwrap();
            let alternative_u =
                solve_intertwiner(&alternative, &source).expect("the same target under a gauge");
            check_vector_contract(&alternative, &source, alternative_u);
            if alternative != canonical {
                saw_a_flipped_coset = true;
            }
        }
        assert!(
            saw_a_flipped_coset,
            "the witness must contain a representative pair that flips the coset phase"
        );
    }

    #[test]
    fn a_single_changed_class_value_is_not_a_valid_target() {
        // A non-self-consistent table (one class changed, not a uniform coset
        // phase) must be refused by the square, conjugation or law check, or have no
        // intertwiner at all.
        let model = witness_model(10, "GM1+", "GM2+", 0, RealSign::Minus);
        let (unitary, corep) = matching_candidate(&model);
        let classes = model.classes();
        let lattice = model.translation_lattice();
        let uni = model.uni();
        let source = source_actions(&model);
        let first_anti = classes.iter().position(|class| class.time_reversal);
        let canonical =
            build_target(&corep, classes, lattice, &unitary, first_anti, uni).unwrap();
        for index in 0..classes.len() {
            let mut broken = canonical.clone();
            broken[index].value = -broken[index].value;
            let refused = check_square(&corep, classes, lattice, &broken, &unitary, uni).is_err()
                || check_conjugation(&corep, classes, lattice, &broken, uni).is_err()
                || check_law(&corep, classes, lattice, &broken, uni).is_err()
                || solve_intertwiner(&broken, &source).is_none();
            assert!(refused, "changing class {index} must not stay a valid target");
        }
    }
}
