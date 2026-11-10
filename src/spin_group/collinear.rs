//! Collinear spin group: signed scalar moments and a full unitary+anti-unitary
//! finite quotient.
//!
//! The magnetic order parameter is `m_i = sigma_i n`; only the signed scalar
//! `sigma_i` and the common axis `n` (up to a U(1) rotation about the axis) are
//! needed.  The finite action on this order parameter is a sign:
//! `Preserve` (`+1`) or `Swap` (`-1`).  The continuous spin rotations about
//! `n` are not sampled; the returned group is the finite quotient modulo that
//! kernel.
//!
//! Every physical collinear magnetic structure also has the channel-preserving
//! anti-unitary operation `A0 = [C2_perp T || E]`, and therefore the detected
//! group always contains both kinds: `n` unitary operations and their `n`
//! anti-unitary `g A0` partners.  This describes the spin symmetry of the
//! **magnetic structure**; whether the Hamiltonian actually realizes the
//! anti-unitary partners additionally requires the absence of SOC/orbital
//! time-reversal breaking, which the consumer (Rustb) must verify.

use super::*;
use crate::Crystal;
use crate::collinear_spin::{
    CollinearSpinOperation as LegacyOperation, CollinearSpinTimeReversalGroup as LegacyFullGroup,
    CollinearSpinTolerances as LegacyTolerances, detect_collinear_spin_group as detect_legacy,
    validate_collinear_spin_group as validate_legacy,
};

/// Net action on a collinear order parameter: the sign multiplying the moment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CollinearAction {
    Preserve,
    Swap,
}

impl CollinearAction {
    /// Net sign multiplying a signed scalar moment.
    pub const fn sign(self) -> i8 {
        match self {
            Self::Preserve => 1,
            Self::Swap => -1,
        }
    }

    const fn from_sign(sign: i8) -> Self {
        if sign >= 0 {
            Self::Preserve
        } else {
            Self::Swap
        }
    }
}

impl SpinAction for CollinearAction {
    fn identity() -> Self {
        Self::Preserve
    }

    fn compose(self, right: Self) -> Self {
        Self::from_sign(self.sign() * right.sign())
    }

    fn inverse(self) -> Self {
        self
    }

    fn residual(self, other: Self) -> f64 {
        if self == other { 0.0 } else { f64::INFINITY }
    }
}

/// Validated collinear order-parameter metadata.
#[derive(Debug, Clone)]
pub struct CollinearOrder {
    signed_moments: Vec<f64>,
    axis: Option<Vec3>,
    collinearity_residuals: Vec<f64>,
}

impl CollinearOrder {
    /// Signed moment magnitudes in the input atom order.
    pub fn signed_moments(&self) -> &[f64] {
        &self.signed_moments
    }
    /// Common Cartesian axis for vector input; `None` for caller-signed input.
    pub fn axis(&self) -> Option<Vec3> {
        self.axis
    }
    /// Per-atom perpendicular norms (zero for caller-signed input).
    pub fn collinearity_residuals(&self) -> &[f64] {
        &self.collinearity_residuals
    }
}

/// Magnetic input for the collinear detector.
#[derive(Debug, Clone, Copy)]
pub enum CollinearMoments<'a> {
    /// Signed magnitudes relative to a caller-chosen common axis.  Only the
    /// sign and magnitude enter the group condition; the axis is not needed.
    Signed(&'a [f64]),
    /// Cartesian vectors: the common axis is inferred and every perpendicular
    /// residual is checked; the signed scalar retains the actual vector norm.
    Vectors(&'a [[f64; 3]]),
}

/// Validated full collinear spin group (unitary and anti-unitary).
pub type CollinearSpinGroup = SpinGroup<CollinearAction, CollinearOrder, f64>;

fn legacy_moments(moments: CollinearMoments<'_>) -> crate::collinear_spin::CollinearMoments<'_> {
    match moments {
        CollinearMoments::Signed(values) => crate::collinear_spin::CollinearMoments::Signed(values),
        CollinearMoments::Vectors(values) => {
            crate::collinear_spin::CollinearMoments::Vectors(values)
        }
    }
}

fn legacy_tolerances(tolerances: SpinTolerances) -> LegacyTolerances {
    LegacyTolerances {
        position: tolerances.position,
        moment: tolerances.moment,
        collinearity: tolerances.direction,
    }
}

fn action_from_legacy(action: crate::collinear_spin::SpinChannelAction) -> CollinearAction {
    match action {
        crate::collinear_spin::SpinChannelAction::Preserve => CollinearAction::Preserve,
        crate::collinear_spin::SpinChannelAction::Swap => CollinearAction::Swap,
    }
}

fn action_to_legacy(action: CollinearAction) -> crate::collinear_spin::SpinChannelAction {
    match action {
        CollinearAction::Preserve => crate::collinear_spin::SpinChannelAction::Preserve,
        CollinearAction::Swap => crate::collinear_spin::SpinChannelAction::Swap,
    }
}

fn order_from_legacy(group: &crate::collinear_spin::ValidatedCollinearSpinGroup) -> CollinearOrder {
    CollinearOrder {
        signed_moments: group.signed_moments().to_vec(),
        axis: group.spin_axis(),
        collinearity_residuals: group.collinearity_residuals().to_vec(),
    }
}

fn group_from_legacy(
    full: LegacyFullGroup,
    order: CollinearOrder,
    tolerances: SpinTolerances,
) -> CollinearSpinGroup {
    let len = full.len();
    let mut operations = Vec::with_capacity(len);
    let mut atom_mappings = Vec::with_capacity(len);
    for index in 0..len {
        let operation = full
            .operation(index)
            .expect("validated spin-group extension index");
        operations.push(SpinOperation {
            rotation: operation.operation.rotation,
            translation: operation.operation.translation,
            action: action_from_legacy(operation.operation.channel_action),
            kind: operation.kind,
        });
        atom_mappings.push(
            operation
                .atom_mapping
                .iter()
                .map(|mapping| SpinAtomMapping {
                    target: mapping.target,
                    cell_shift: mapping.cell_shift,
                    position_residual: mapping.position_residual,
                    moment_residual: mapping.moment_residual,
                })
                .collect(),
        );
    }
    let mut multiplication_table = Vec::with_capacity(len);
    for left in 0..len {
        let mut row = Vec::with_capacity(len);
        for right in 0..len {
            let product = full
                .product(left, right)
                .expect("validated spin-group product");
            row.push(SpinGroupProduct {
                operation_index: product.operation_index,
                cell_shift: product.cell_shift,
                translation_residual: product.translation_residual,
                action_residual: 0.0,
            });
        }
        multiplication_table.push(row);
    }
    let mut inverse_indices = Vec::with_capacity(len);
    for index in 0..len {
        inverse_indices.push(
            full.inverse_index(index)
                .expect("validated spin-group inverse"),
        );
    }
    SpinGroup {
        kind: SpinGroupKind::Collinear,
        operations,
        atom_mappings,
        multiplication_table,
        identity_index: full.identity_index(),
        inverse_indices,
        tolerances,
        order,
    }
}

/// Detect the full collinear spin group (unitary and anti-unitary) from a
/// crystal and its moments.
///
/// The unitary subgroup is the finite quotient of the ordinary collinear spin
/// group; the anti-unitary half is the `g A0` coset with
/// `A0 = [C2_perp T || E]`.  The latter is a symmetry of every collinear moment
/// configuration; whether a Hamiltonian with SOC/orbital terms realizes it is
/// outside this geometric certificate.
pub fn detect_collinear_spin_group(
    crystal: &Crystal,
    moments: CollinearMoments<'_>,
    tolerances: SpinTolerances,
) -> Result<CollinearSpinGroup, SpinGroupError> {
    tolerances.validate()?;
    let unitary = detect_legacy(
        crystal,
        legacy_moments(moments),
        legacy_tolerances(tolerances),
    )?;
    let order = order_from_legacy(&unitary);
    let full = unitary.with_spin_only_time_reversal();
    Ok(group_from_legacy(full, order, tolerances))
}

/// Validate an explicit full collinear operation set.
///
/// The set must contain equally many unitary and anti-unitary operations and
/// every anti-unitary operation must be the `A0` partner of one unitary
/// operation.  The unitary subset is validated by the existing witness-bearing
/// validator; the whole set is never pruned to a smaller group.  The returned
/// certificate uses the canonical order (unitary operations in the validated
/// target order, then their anti-unitary partners).
pub fn validate_collinear_spin_group(
    crystal: &Crystal,
    moments: CollinearMoments<'_>,
    target: &[SpinOperation<CollinearAction>],
    tolerances: SpinTolerances,
) -> Result<CollinearSpinGroup, SpinGroupError> {
    tolerances.validate()?;
    let unitary_target = target
        .iter()
        .enumerate()
        .filter(|(_, operation)| operation.kind.is_unitary())
        .map(|(_, operation)| LegacyOperation {
            rotation: operation.rotation,
            translation: operation.translation,
            channel_action: action_to_legacy(operation.action),
        })
        .collect::<Vec<_>>();
    let antiunitary = target.len().saturating_sub(unitary_target.len());
    if unitary_target.is_empty() || antiunitary != unitary_target.len() {
        return Err(SpinGroupError::KindCount {
            unitary: unitary_target.len(),
            antiunitary,
        });
    }
    let validated = validate_legacy(
        crystal,
        legacy_moments(moments),
        &unitary_target,
        legacy_tolerances(tolerances),
    )?;
    let order = order_from_legacy(&validated);
    let full = validated.with_spin_only_time_reversal();
    let group = group_from_legacy(full, order, tolerances);
    debug_assert_eq!(group.len(), 2 * unitary_target.len());
    // Every anti-unitary target must be a unique companion of a unitary target.
    let mut matched = vec![false; group.len()];
    for (target_index, operation) in target.iter().enumerate() {
        if operation.kind.is_unitary() {
            continue;
        }
        let companion = group
            .operations()
            .iter()
            .enumerate()
            .skip(group.len() / 2)
            .find(|(index, candidate)| {
                !matched[*index]
                    && candidate.rotation == operation.rotation
                    && candidate.translation == operation.translation
                    && candidate.action == operation.action
            });
        let Some((index, _)) = companion else {
            return Err(SpinGroupError::AntiunitaryPartner {
                operation_index: target_index,
            });
        };
        matched[index] = true;
    }
    if matched.iter().skip(group.len() / 2).any(|matched| !matched) {
        return Err(SpinGroupError::KindCount {
            unitary: unitary_target.len(),
            antiunitary,
        });
    }
    Ok(group)
}
