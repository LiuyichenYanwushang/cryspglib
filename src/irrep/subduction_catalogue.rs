//! Exact little co-groups and their one-dimensional projective characters
//! (R4 batch 2a).
//!
//! A gap target whose folded point `q` has a non-trivial little co-group still
//! has an exactly computable target catalogue.  The child's little-group
//! irreps at `q` are the projective irreps of the little co-group `P_q` with
//! the factor system
//!
//! ```text
//! omega_ij = exp(2 pi i q . L_ij),   L_ij = t_i + R_i t_j - t_k
//! ```
//!
//! where `L_ij` is the **unreduced** defect of the product of the chosen
//! representatives `s_i = (R_i, t_i)` against the representative `s_k` of the
//! same rotation -- exactly the convention `docs/subduction-gap-sources.md`
//! pins for the R3 classifier, and the reason the reduction must not be applied
//! early: for `q` that is not a reciprocal lattice vector, reducing `L_ij`
//! changes the phase.
//!
//! This module computes, with exact rational arithmetic and without any
//! character table:
//!
//! * [`little_co_group`]: one child operation per distinct rotation that fixes
//!   `q` modulo the child primitive reciprocal lattice (centring extinctions
//!   included), with its factor system as rational turns;
//! * [`one_dimensional_characters`]: every one-dimensional projective character
//!   `psi` of `P_q`, i.e. every exact solution of
//!   `psi_i + psi_j - psi_k == omega_ij (mod 1)`.
//!
//! The caller only uses a catalogue when the solver returns exactly `|P_q|`
//! solutions.  The solution set is either empty or a coset of
//! `Hom(P_q, U(1))`, whose size is `|P_q / [P_q, P_q]|`; requiring that size to
//! be `|P_q|` therefore means `P_q` is abelian *and* the cocycle is a
//! coboundary, which is exactly the case where **every** irreducible projective
//! representation of `P_q` is one-dimensional.  A co-group whose cocycle is not a coboundary
//! (needing a two-dimensional projective irrep) and a non-abelian co-group
//! return fewer solutions and are left to the higher-dimensional batch: this
//! module never guesses a catalogue.
//!
//! The caller must build the constants with the **same** exact point the
//! character is evaluated at.  A reconstructed little-group operation is not a
//! lattice translate of its representative (its translation can carry quarters),
//! so taking the constants at the raw folded point and evaluating at the reduced
//! one changes the phase: that mix produced `1 +- i` multiplicities on 30 probes
//! (ordinal 14090, SG 226 W5 -> #98) before it was fixed.

use crate::mathfunc::Mat3I;

use super::super::{ExactSeitz, Lattice, Rat, SubductionError, Vec3R, exact_primitive_basis};
use super::super::{bloch_phase, strict_sg_hall_ops};
use super::StarError;

/// The identity rotation, as stored in every Hall operation table.
const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// Upper bound on the search grid used per generator.
///
/// The co-groups behind the R4 gap groups have order at most four here, and one
/// or two generators, so the bound is never reached in practice; it keeps a
/// pathological co-group from turning the solver into a hot loop.  Exceeding it
/// is a fail-closed empty catalogue, never a partial one.
const MAX_GRID: i128 = 200_000;

/// Largest little co-group the exact solver searches.
///
/// R4 batch 2a needs the one-dimensional solution set for `|P_q| <= 4`; batch 2b
/// also needs it for the six-element `D_3` family, where it supplies the gauge
/// (its two solutions are `Hom(D_3, U(1))`) rather than a catalogue.  Anything
/// bigger is out of scope: the search below would grow with the order, so it
/// returns an empty catalogue and keeps the missing-data behaviour.
const MAX_ORDER: usize = 6;

/// Upper bound on `combinations x |P_q|^2`, the solver's work per co-group.
const MAX_WORK: i128 = 4_000_000;

/// The little co-group of one exact child point.
#[derive(Debug, Clone)]
pub(super) struct LittleCoGroup {
    /// One child operation per distinct rotation fixing `q`, identity first.
    representatives: Vec<ExactSeitz>,
    /// `turns[i][j]` in `[0, 1)`: the cocycle phase `q . L_ij` of the product
    /// `s_i s_j` measured against the representative of the same rotation.
    turns: Vec<Vec<Rat>>,
    /// The exact point the co-group belongs to.
    q: Vec3R,
}

impl LittleCoGroup {
    /// Number of distinct rotations fixing `q`.
    pub(super) fn order(&self) -> usize {
        self.representatives.len()
    }

    /// Position of one rotation's representative.
    fn position(&self, rotation: Mat3I) -> Option<usize> {
        self.representatives
            .iter()
            .position(|operation| operation.rotation() == rotation)
    }

    /// Position of the product rotation `R_i R_j`, or a fail-closed error.
    fn product_position(&self, i: usize, j: usize) -> Result<usize, StarError> {
        let product = self.representatives[i].compose(&self.representatives[j])?;
        self.position(product.rotation())
            .ok_or_else(|| StarError::LittleCoGroupNotClosed {
                q: [self.q.get(0), self.q.get(1), self.q.get(2)],
            })
    }
}

/// One exact little co-group at the child point `q`.
pub(super) fn little_co_group(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<LittleCoGroup, StarError> {
    let child_cell = Lattice::new(exact_primitive_basis(child_sg)?)?;
    let hall = strict_sg_hall_ops(child_sg)?;
    let mut representatives: Vec<ExactSeitz> = Vec::new();
    for operation in &hall.operations {
        if !child_reciprocal.preserves(operation.rotation(), q)? {
            continue;
        }
        if representatives
            .iter()
            .any(|kept| kept.rotation() == operation.rotation())
        {
            continue;
        }
        representatives.push(*operation);
    }
    // The identity rotation always fixes `q`.  Its representative must be the
    // *zero-translation* identity: the gauge fixes `psi(identity) = 0`, so a
    // centring representative here would silently drop its Bloch phase.  Any
    // other identity-rotation operation is reached again through the phase of
    // its own translation.
    let identity = representatives
        .iter()
        .position(|operation| {
            operation.rotation() == IDENTITY_ROTATION && operation.translation().is_zero()
        })
        .ok_or_else(|| StarError::LittleCoGroupNotClosed {
            q: [q.get(0), q.get(1), q.get(2)],
        })?;
    representatives.swap(0, identity);
    let order = representatives.len();
    let mut turns = vec![vec![Rat::ZERO; order]; order];
    for i in 0..order {
        for j in 0..order {
            let product = representatives[i].compose(&representatives[j])?;
            let position = representatives
                .iter()
                .position(|kept| kept.rotation() == product.rotation())
                .ok_or_else(|| StarError::LittleCoGroupNotClosed {
                    q: [q.get(0), q.get(1), q.get(2)],
                })?;
            // `s_i s_j` and the representative of the same rotation are the same
            // element of the little group, so they differ by a lattice vector.
            let defect = product
                .translation()
                .checked_sub(representatives[position].translation())?;
            if !child_cell.contains(&defect)? {
                return Err(StarError::LittleCoGroupNotClosed {
                    q: [q.get(0), q.get(1), q.get(2)],
                });
            }
            // The unreduced product defect carries the phase; reducing it first
            // would silently drop it.
            let mut value = Rat::ZERO;
            for axis in 0..3 {
                value = value.checked_add(q.get(axis).checked_mul(defect.get(axis))?)?;
            }
            turns[i][j] = fractional(value)?;
        }
    }
    Ok(LittleCoGroup {
        representatives,
        turns,
        q: *q,
    })
}

/// Every one-dimensional projective character of the co-group.
///
/// The result is the complete solution set of
/// `psi_i + psi_j - psi_k == turns[i][j] (mod 1)`, with `psi` in `[0, 1)`.
/// Every returned vector is verified against every equation, so a returned
/// catalogue is never wrong; a search that cannot finish (too many generators
/// or a grid beyond [`MAX_GRID`]) returns an empty catalogue instead of a
/// partial one.
pub(super) fn one_dimensional_characters(
    co_group: &LittleCoGroup,
) -> Result<Vec<Vec<Rat>>, StarError> {
    let order = co_group.order();
    if order > MAX_ORDER {
        return Ok(Vec::new());
    }
    let identity = 0usize; // `little_co_group` puts the identity first.
    let generators = co_group.generators()?;
    let mut modulus = 1i128;
    for row in &co_group.turns {
        for value in row {
            modulus = lcm(modulus, value.denominator());
        }
    }
    // A gauge can be finer than the cocycle (a C2 cocycle of 1/3 is gauged by a
    // character of 1/6), so the grid is scaled by the group order.
    modulus = modulus
        .checked_mul(i128::try_from(order).map_err(|_| {
            StarError::Subduction(SubductionError::RationalOverflow {
                operation: "little co-group order",
            })
        })?)
        .ok_or(StarError::Subduction(SubductionError::RationalOverflow {
            operation: "cocycle modulus",
        }))?;
    if generators.is_empty() || modulus <= 0 {
        return Ok(Vec::new());
    }
    let mut combinations = 1i128;
    for _ in &generators {
        combinations = combinations.saturating_mul(modulus);
        if combinations > MAX_GRID
            || combinations.saturating_mul((order * order) as i128) > MAX_WORK
        {
            return Ok(Vec::new());
        }
    }
    let mut solutions = Vec::new();
    let mut odometer = vec![0i128; generators.len()];
    loop {
        if let Some(psi) = co_group.complete(&generators, &odometer, modulus, identity)?
            && !solutions.contains(&psi)
        {
            solutions.push(psi);
        }
        // Advance the odometer.
        let mut position = 0usize;
        loop {
            if position == odometer.len() {
                return Ok(solutions);
            }
            odometer[position] += 1;
            if odometer[position] < modulus {
                break;
            }
            odometer[position] = 0;
            position += 1;
        }
    }
}

impl LittleCoGroup {
    /// A generating set of positions, greedily closed under products.
    fn generators(&self) -> Result<Vec<usize>, StarError> {
        let order = self.order();
        let mut generators: Vec<usize> = Vec::new();
        let mut span = vec![0usize]; // the identity is position zero
        for candidate in 0..order {
            if span.contains(&candidate) {
                continue;
            }
            generators.push(candidate);
            span = self.close_span(&span, &generators)?;
            if span.len() == order {
                break;
            }
        }
        if span.len() != order {
            return Err(StarError::LittleCoGroupNotClosed {
                q: [self.q.get(0), self.q.get(1), self.q.get(2)],
            });
        }
        Ok(generators)
    }

    /// The subgroup generated by `generators` together with `span`, as positions.
    fn close_span(&self, span: &[usize], generators: &[usize]) -> Result<Vec<usize>, StarError> {
        let mut closed = span.to_vec();
        let mut changed = true;
        while changed {
            changed = false;
            for left in closed.clone() {
                for right in generators.iter().copied().chain(std::iter::once(left)) {
                    let product = self.product_position(left, right)?;
                    if !closed.contains(&product) {
                        closed.push(product);
                        changed = true;
                    }
                }
            }
        }
        closed.sort_unstable();
        Ok(closed)
    }

    /// Complete one generator assignment into a full character, verified on
    /// every equation; `None` when the assignment is inconsistent.
    fn complete(
        &self,
        generators: &[usize],
        assignment: &[i128],
        modulus: i128,
        identity: usize,
    ) -> Result<Option<Vec<Rat>>, StarError> {
        let order = self.order();
        let mut psi: Vec<Option<Rat>> = vec![None; order];
        psi[identity] = Some(Rat::ZERO);
        for (generator, value) in generators.iter().zip(assignment) {
            psi[*generator] = Some(Rat::new(value.rem_euclid(modulus), modulus)?);
        }
        loop {
            let mut changed = false;
            for i in 0..order {
                for j in 0..order {
                    let k = self.product_position(i, j)?;
                    let (left, right, target) = (psi[i], psi[j], psi[k]);
                    let filled = match (left, right, target) {
                        (Some(a), Some(b), None) => Some((k, a.checked_add(b)?.checked_sub(self.turns[i][j])?)),
                        (Some(a), None, Some(c)) => Some((j, c.checked_add(self.turns[i][j])?.checked_sub(a)?)),
                        (None, Some(b), Some(c)) => Some((i, c.checked_add(self.turns[i][j])?.checked_sub(b)?)),
                        _ => None,
                    };
                    if let Some((slot, value)) = filled {
                        let value = fractional(value)?;
                        match psi[slot] {
                            Some(existing) if existing != value => return Ok(None),
                            Some(_) => {}
                            None => {
                                psi[slot] = Some(value);
                                changed = true;
                            }
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }
        let mut values = Vec::with_capacity(order);
        for value in psi {
            match value {
                Some(value) => values.push(value),
                None => return Ok(None),
            }
        }
        // Verify every equation of the completed character.
        for i in 0..order {
            for j in 0..order {
                let k = self.product_position(i, j)?;
                let residual = fractional(
                    values[i]
                        .checked_add(values[j])?
                        .checked_sub(values[k])?
                        .checked_sub(self.turns[i][j])?,
                )?;
                if !residual.is_zero() {
                    return Ok(None);
                }
            }
        }
        Ok(Some(values))
    }
}

/// One constructed target of a **higher-dimensional** little co-group.
pub(super) struct ProjectiveTarget {
    /// Complex dimension of the little-group irrep.
    pub dimension: u8,
    /// `(rotation, chi_rep(R) * exp(-2 pi i q.t_R))` per little co-group
    /// rotation, so that `chi(R, T) = constant(R) * exp(+2 pi i q.T)` evaluates
    /// the little-group character on any reconstructed operation exactly as the
    /// one-dimensional catalogue does.
    pub constants: Vec<(Mat3I, num_complex::Complex64)>,
}

/// The projective character table of a small little co-group.
///
/// Two families are covered, each behind its own structural gate, and anything
/// else returns an empty table (the caller keeps `MissingChildStarData`):
///
/// * `|P_q| = 4`, abelian, every non-identity element of order two, and a
///   **non-degenerate** commutator pairing `beta_ij = phi_ij - phi_ji`.  Then
///   the twisted group algebra is `M_2(C)`, so there is exactly one irreducible
///   representation: dimension two and character `(2, 0, 0, 0)`.  The trace
///   argument does **not** use `omega(g,g) = +-1`: `g^2 = e` holds only for the
///   rotation part, and the cocycle stores the phase of the lattice translation
///   the representative product picks up, so `omega(g,g)` is an arbitrary root
///   of unity (the archived groups really show `1/6`, `1/4`, `1/3`, `2/3`,
///   `3/4` and `5/6`).  What is used instead is `u_g^2 = omega(g,g) * 1` with
///   `omega(g,g) != 0`: the eigenvalues of `u_g` are `+-sqrt(omega(g,g))`, which
///   are distinct, and `u_g` cannot be scalar -- a scalar `u_g` would force
///   `omega(g,h) = omega(h,g)` for every `h`, making `g` orthogonal to the whole
///   group and contradicting the pairing's non-degeneracy -- so both eigenspaces
///   are one-dimensional and the trace vanishes.
/// * `|P_q| = 6`, non-abelian (hence `D_3`), with a **coboundary** cocycle: the
///   gauge `psi` exists and has exactly two solutions (the size of
///   `Hom(D_3, U(1))`), and the projective irreps are the gauged ordinary ones,
///   with dimensions `{2, 1, 1}`.  The *set* of targets does not depend on which
///   solution is used, because the two solutions differ by the sign character
///   and tensoring with it permutes the ordinary irreps.
pub(super) fn projective_targets(
    child_sg: u8,
    q: &Vec3R,
    child_reciprocal: &Lattice,
) -> Result<Vec<ProjectiveTarget>, StarError> {
    let co_group = little_co_group(child_sg, q, child_reciprocal)?;
    let order = co_group.order();
    let q = co_group.q;
    let identity = 0usize; // `little_co_group` puts the identity first.
    if order == 4 {
        return abelian_four_targets(&co_group, identity, &q);
    }
    if order == 6 {
        return dihedral_six_targets(&co_group, identity, &q);
    }
    Ok(Vec::new())
}

/// `(rotation, chi_rep(R) * exp(-2 pi i q.t_R))` for one co-group character row.
fn constants_of(
    co_group: &LittleCoGroup,
    q: &Vec3R,
    traces: &[num_complex::Complex64],
) -> Result<Vec<(Mat3I, num_complex::Complex64)>, StarError> {
    let mut out = Vec::with_capacity(co_group.order());
    for (position, operation) in co_group.representatives.iter().enumerate() {
        let phase = bloch_phase(q, operation.translation())?.conj();
        out.push((operation.rotation(), traces[position] * phase));
    }
    Ok(out)
}

/// The one-class, four-element abelian family (see [`projective_targets`]).
fn abelian_four_targets(
    co_group: &LittleCoGroup,
    identity: usize,
    q: &Vec3R,
) -> Result<Vec<ProjectiveTarget>, StarError> {
    let order = co_group.order();
    // Abelian, and every non-identity element of order two (so the group is
    // C2 x C2, not C4).
    for i in 0..order {
        for j in 0..order {
            let product = co_group.product_position(i, j)?;
            let other = co_group.product_position(j, i)?;
            if product != other {
                return Ok(Vec::new());
            }
        }
    }
    for i in 0..order {
        if i == identity {
            continue;
        }
        if co_group.product_position(i, i)? != identity {
            return Ok(Vec::new());
        }
    }
    // Non-degenerate commutator pairing: only the identity is orthogonal to
    // everything, which is what makes the twisted algebra M_2(C).
    for i in 0..order {
        let orthogonal = (0..order).all(|j| {
            let beta = co_group.turns[i][j].checked_sub(co_group.turns[j][i]);
            matches!(beta, Ok(value) if fractional(value).is_ok_and(|v| v.is_zero()))
        });
        if orthogonal != (i == identity) {
            return Ok(Vec::new());
        }
    }
    let mut traces = vec![num_complex::Complex64::new(0.0, 0.0); order];
    traces[identity] = num_complex::Complex64::new(2.0, 0.0);
    Ok(vec![ProjectiveTarget {
        dimension: 2,
        constants: constants_of(co_group, q, &traces)?,
    }])
}

/// The six-element dihedral family (see [`projective_targets`]).
fn dihedral_six_targets(
    co_group: &LittleCoGroup,
    identity: usize,
    q: &Vec3R,
) -> Result<Vec<ProjectiveTarget>, StarError> {
    let order = co_group.order();
    // Element orders identify the two three-fold and three two-fold elements.
    let mut orders = Vec::with_capacity(order);
    for i in 0..order {
        let mut power = i;
        let mut value = 1usize;
        while power != identity {
            power = co_group.product_position(power, i)?;
            value += 1;
            if value > order {
                return Ok(Vec::new());
            }
        }
        orders.push(value);
    }
    if orders[identity] != 1
        || orders.iter().filter(|value| **value == 3).count() != 2
        || orders.iter().filter(|value| **value == 2).count() != 3
    {
        return Ok(Vec::new());
    }
    // Coboundary: the one-dimensional characters must exist, and there must be
    // exactly |Hom(D3, U(1))| = 2 of them.
    let gauge = one_dimensional_characters(co_group)?;
    if gauge.len() != 2 {
        return Ok(Vec::new());
    }
    let psi = &gauge[0];
    let sign: Vec<f64> = orders
        .iter()
        .map(|value| if *value == 2 { -1.0 } else { 1.0 })
        .collect();
    let standard: Vec<f64> = orders
        .iter()
        .map(|value| match value {
            1 => 2.0,
            3 => -1.0,
            _ => 0.0,
        })
        .collect();
    let mut targets = Vec::with_capacity(3);
    for (dimension, values) in [(1u8, vec![1.0; order]), (1, sign), (2, standard)] {
        let traces: Vec<num_complex::Complex64> = psi
            .iter()
            .zip(&values)
            .map(|(psi, value)| {
                num_complex::Complex64::from_polar(1.0, std::f64::consts::TAU * psi.to_f64())
                    * value
            })
            .collect();
        // Gate: every projective character has norm one, and the dimensions
        // exhaust the group.
        let norm: f64 = traces.iter().map(|value| value.norm_sqr()).sum::<f64>()
            / f64::from(u32::try_from(order).map_err(|_| {
                StarError::Subduction(SubductionError::RationalOverflow {
                    operation: "little co-group order",
                })
            })?);
        if (norm - 1.0).abs() > 1e-9 {
            return Ok(Vec::new());
        }
        targets.push(ProjectiveTarget {
            dimension,
            constants: constants_of(co_group, q, &traces)?,
        });
    }
    let total: u32 = targets
        .iter()
        .map(|target| u32::from(target.dimension) * u32::from(target.dimension))
        .sum();
    if total != u32::try_from(order).unwrap_or(0) {
        return Ok(Vec::new());
    }
    Ok(targets)
}

/// The per-rotation constants `psi_R - q . t_R` of one catalogue character.
///
/// With them the character of **any** little-group operation `(R, T)` is
/// `exp(2 pi i (constant(R) + q . T))`: the representative itself contributes
/// `exp(2 pi i psi_R)`, and every other operation of the same rotation differs
/// from it by a lattice translation, whose Bloch phase is exactly the extra
/// `q . (T - t_R)`.
pub(super) fn constants(
    co_group: &LittleCoGroup,
    psi: &[Rat],
) -> Result<Vec<(Mat3I, Rat)>, StarError> {
    let mut out = Vec::with_capacity(co_group.order());
    for (position, operation) in co_group.representatives.iter().enumerate() {
        let mut constant = psi[position];
        for axis in 0..3 {
            constant = constant.checked_sub(
                co_group
                    .q
                    .get(axis)
                    .checked_mul(operation.translation().get(axis))?,
            )?;
        }
        out.push((operation.rotation(), fractional(constant)?));
    }
    Ok(out)
}

/// The character of one little-group operation under one catalogue entry.
pub(super) fn character_value(
    constants: &[(Mat3I, Rat)],
    q: &Vec3R,
    operation: &ExactSeitz,
) -> Result<num_complex::Complex64, StarError> {
    let rotation = operation.rotation();
    let constant = constants
        .iter()
        .find(|(known, _)| *known == rotation)
        .map(|(_, value)| *value)
        .ok_or(StarError::ConstructedRotationNotCovered {
            q: [q.get(0), q.get(1), q.get(2)],
        })?;
    let mut turns = constant;
    for axis in 0..3 {
        turns = turns.checked_add(
            q.get(axis)
                .checked_mul(operation.translation().get(axis))?,
        )?;
    }
    Ok(num_complex::Complex64::from_polar(
        1.0,
        std::f64::consts::TAU * turns.to_f64(),
    ))
}

/// The character of one little-group operation under one table entry:
/// `constant(R) * exp(2 pi i q.T)`, with `constant(R)` carrying the reduced
/// point's representative phase.
pub(super) fn table_character_value(
    constants: &[(Mat3I, num_complex::Complex64)],
    q: &Vec3R,
    operation: &ExactSeitz,
) -> Result<num_complex::Complex64, StarError> {
    let rotation = operation.rotation();
    let constant = constants
        .iter()
        .find(|(known, _)| *known == rotation)
        .map(|(_, value)| *value)
        .ok_or(StarError::ConstructedRotationNotCovered {
            q: [q.get(0), q.get(1), q.get(2)],
        })?;
    Ok(constant * bloch_phase(q, operation.translation())?)
}

/// The fractional part in `[0, 1)` of an exact rational.
pub(super) fn fractional(value: Rat) -> Result<Rat, SubductionError> {
    let denominator = value.denominator();
    let remainder = value.numerator().rem_euclid(denominator);
    Rat::new(remainder, denominator)
}

/// Least common multiple of two positive integers.
fn lcm(left: i128, right: i128) -> i128 {
    if left == 0 || right == 0 {
        return 0;
    }
    let mut a = left.abs();
    let mut b = right.abs();
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    (left.abs() / a).saturating_mul(right.abs())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::irrep::query;
    use crate::irrep::subduction::{SUBDUCTION_TOLERANCE, character_of, inline_k_vector};
    use num_complex::Complex64;

    /// A cyclic co-group of order two with a prescribed generator cocycle.
    fn cyclic_two(group_turn: Rat) -> LittleCoGroup {
        LittleCoGroup {
            representatives: vec![
                ExactSeitz::identity(),
                ExactSeitz::new([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], Vec3R::zero()),
            ],
            turns: vec![vec![Rat::ZERO, Rat::ZERO], vec![Rat::ZERO, group_turn]],
            q: Vec3R::zero(),
        }
    }

    /// `2 psi_g == phi_gg` has exactly two solutions on the unit circle, so a
    /// two-fold co-group always has two one-dimensional characters.
    #[test]
    fn a_two_fold_co_group_has_two_characters() {
        let characters = one_dimensional_characters(&cyclic_two(Rat::ZERO)).unwrap();
        assert_eq!(characters.len(), 2);
        assert!(characters.contains(&vec![Rat::ZERO, Rat::ZERO]));
        assert!(characters.contains(&vec![Rat::ZERO, Rat::new(1, 2).unwrap()]));

        // A half-turn cocycle moves them to the two imaginary characters.
        let characters =
            one_dimensional_characters(&cyclic_two(Rat::new(1, 2).unwrap())).unwrap();
        assert_eq!(characters.len(), 2);
        assert!(characters.contains(&vec![Rat::ZERO, Rat::new(1, 4).unwrap()]));
        assert!(characters.contains(&vec![Rat::ZERO, Rat::new(3, 4).unwrap()]));

        // A one-third turn is gauged by a sixth: the grid must be finer than the
        // cocycle itself, or the character would be missed.
        let characters =
            one_dimensional_characters(&cyclic_two(Rat::new(1, 3).unwrap())).unwrap();
        assert_eq!(characters.len(), 2);
        assert!(characters.contains(&vec![Rat::ZERO, Rat::new(1, 6).unwrap()]));
        assert!(characters.contains(&vec![Rat::ZERO, Rat::new(2, 3).unwrap()]));
    }

    /// The character values of one target on every co-group representative, in
    /// the co-group's canonical order.
    fn character_vector(target: &ProjectiveTarget, co_group: &LittleCoGroup) -> Vec<Complex64> {
        co_group
            .representatives
            .iter()
            .map(|operation| {
                table_character_value(&target.constants, &co_group.q, operation).unwrap()
            })
            .collect()
    }

    fn close(left: &[Complex64], right: &[Complex64]) -> bool {
        left.len() == right.len()
            && left
                .iter()
                .zip(right)
                .all(|(a, b)| (a - b).norm() <= SUBDUCTION_TOLERANCE)
    }

    /// The two batch-2b families, described by their structure rather than by a
    /// pinned row: a non-degenerate four-element co-group has exactly one
    /// two-dimensional irrep with character `(2, 0, 0, 0)`, and `D_3` with a
    /// coboundary cocycle has the gauged ordinary irreps `{1, 1, 2}`.
    #[test]
    fn the_projective_tables_cover_only_the_two_gated_families() {
        let cell = Lattice::new(exact_primitive_basis(43).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        // Ordinal 3988 folds onto child #43 with a four-element C2 x C2
        // co-group whose commutator pairing is non-degenerate.
        let four_q = Vec3R::new([
            Rat::new(0, 1).unwrap(),
            Rat::new(1, 1).unwrap(),
            Rat::new(1, 2).unwrap(),
        ]);
        let four = projective_targets(43, &four_q, &reciprocal).unwrap();
        assert_eq!(four.len(), 1, "one two-dimensional irrep");
        assert_eq!(four[0].dimension, 2);
        // The whole character row, not only its identity entry: two on the
        // identity and zero on every other co-group element.  The trace
        // argument covers arbitrary `omega(g,g)` (see `projective_targets`), so
        // the pinned zeros are the load-bearing part.
        let four_group = little_co_group(43, &four_q, &reciprocal).unwrap();
        assert_eq!(four_group.order(), 4);
        assert!(
            close(
                &character_vector(&four[0], &four_group),
                &[2.0.into(), 0.0.into(), 0.0.into(), 0.0.into()]
            ),
            "{:?}",
            character_vector(&four[0], &four_group)
        );

        // A six-element non-abelian co-group: child #160 at (0, 0, 3/4).
        let cell = Lattice::new(exact_primitive_basis(160).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        let six_q = Vec3R::new([
            Rat::new(0, 1).unwrap(),
            Rat::new(0, 1).unwrap(),
            Rat::new(3, 4).unwrap(),
        ]);
        let six = projective_targets(160, &six_q, &reciprocal).unwrap();
        let mut dimensions: Vec<u8> = six.iter().map(|target| target.dimension).collect();
        dimensions.sort_unstable();
        assert_eq!(dimensions, vec![1, 1, 2], "gauged ordinary D3 irreps");
        let total: u32 = six
            .iter()
            .map(|target| u32::from(target.dimension) * u32::from(target.dimension))
            .sum();
        assert_eq!(total, 6, "the irreps exhaust the co-group");
        // The full table in the co-group's canonical order: the two one
        // dimensional targets are the trivial and the sign representation, and
        // the two dimensional one is the standard representation -- `(2, -1, -1,
        // 0, 0, 0)` is exactly the vector reviewer B asked to pin (the order is
        // identity, the two three-fold rotations, the three two-fold ones).
        let six_group = little_co_group(160, &six_q, &reciprocal).unwrap();
        assert_eq!(six_group.order(), 6);
        let mut rows: Vec<(u8, Vec<Complex64>)> = six
            .iter()
            .map(|target| (target.dimension, character_vector(target, &six_group)))
            .collect();
        rows.sort_by_key(|row| {
            let negative = row.1.iter().filter(|value| value.re < 0.0).count();
            (row.0, negative)
        });
        let one = Complex64::new(1.0, 0.0);
        let minus_one = Complex64::new(-1.0, 0.0);
        let zero = Complex64::new(0.0, 0.0);
        let two = Complex64::new(2.0, 0.0);
        assert_eq!(rows.len(), 3);
        assert!(close(&rows[0].1, &[one, one, one, one, one, one]), "{:?}", rows[0]);
        assert!(
            close(&rows[1].1, &[one, one, one, minus_one, minus_one, minus_one]),
            "{:?}",
            rows[1]
        );
        assert!(
            close(&rows[2].1, &[two, minus_one, minus_one, zero, zero, zero]),
            "{:?}",
            rows[2]
        );
        assert_eq!(rows.iter().map(|row| row.0).collect::<Vec<_>>(), vec![1, 1, 2]);

        // Out of scope: a cubic Gamma point (order 48) has no table.
        let cell = Lattice::new(exact_primitive_basis(221).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        assert!(
            projective_targets(221, &Vec3R::new([Rat::new(0, 1).unwrap(); 3]), &reciprocal)
                .unwrap()
                .is_empty()
        );
    }

    /// The one-dimensional solver's negative direction, which no real data
    /// exercises: it is only ever *used* when it returns a full catalogue, so a
    /// partial or empty answer must be reachable and must stay empty.
    #[test]
    fn the_one_dimensional_solver_returns_nothing_instead_of_a_subset() {
        // A four-element co-group with a non-degenerate pairing has no
        // one-dimensional projective representation at all: the solutions of
        // `psi_i + psi_j - psi_k == turns_ij` number zero, not "a few".
        let cell = Lattice::new(exact_primitive_basis(43).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        let q = Vec3R::new([
            Rat::new(0, 1).unwrap(),
            Rat::new(1, 1).unwrap(),
            Rat::new(1, 2).unwrap(),
        ]);
        let co_group = little_co_group(43, &q, &reciprocal).unwrap();
        assert_eq!(co_group.order(), 4);
        let characters = one_dimensional_characters(&co_group).unwrap();
        assert!(
            characters.len() < co_group.order(),
            "a non-coboundary cocycle must not yield a full one-dimensional set"
        );
        assert!(characters.is_empty(), "{characters:?}");

        // `MAX_ORDER` is a real gate: the order-48 cubic Gamma point returns an
        // empty set instead of a partial one.
        let cell = Lattice::new(exact_primitive_basis(221).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        let gamma = Vec3R::new([Rat::new(0, 1).unwrap(); 3]);
        let cubic = little_co_group(221, &gamma, &reciprocal).unwrap();
        assert_eq!(cubic.order(), 48);
        assert!(one_dimensional_characters(&cubic).unwrap().is_empty());

        // The grid/work caps are the second fail-closed layer.  The archive's
        // cocycles live on a 1/12 grid, so no real co-group comes near the cap
        // (the largest real search is `12 x 6 = 72` per generator, well under
        // `MAX_GRID`); this synthetic Klein four with a 1/512 cocycle reaches it
        // and must return nothing rather than a subset.
        let scaled = LittleCoGroup {
            representatives: vec![
                ExactSeitz::identity(),
                ExactSeitz::new([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], Vec3R::zero()),
                ExactSeitz::new([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], Vec3R::zero()),
                ExactSeitz::new([[1, 0, 0], [0, -1, 0], [0, 0, -1]], Vec3R::zero()),
            ],
            turns: vec![
                vec![Rat::ZERO; 4],
                vec![
                    Rat::ZERO,
                    Rat::ZERO,
                    Rat::new(1, 512).unwrap(),
                    Rat::new(1, 512).unwrap(),
                ],
                vec![
                    Rat::ZERO,
                    Rat::new(511, 512).unwrap(),
                    Rat::ZERO,
                    Rat::new(1, 512).unwrap(),
                ],
                vec![
                    Rat::ZERO,
                    Rat::new(511, 512).unwrap(),
                    Rat::new(511, 512).unwrap(),
                    Rat::ZERO,
                ],
            ],
            q: Vec3R::zero(),
        };
        let mut modulus = 1i128;
        for row in &scaled.turns {
            for value in row {
                modulus = lcm(modulus, value.denominator());
            }
        }
        assert!(
            modulus.saturating_mul(modulus) > MAX_GRID,
            "the witness must exceed MAX_GRID, modulus {modulus}"
        );
        assert!(one_dimensional_characters(&scaled).unwrap().is_empty());
    }

    /// The decisive evidence for the catalogue: at every pinned child `k` whose
    /// little co-group is non-trivial and whose irreps are all one-dimensional,
    /// the engine's own validated little-group character row must be **one of**
    /// the catalogue characters, on every little-group operation.
    ///
    /// This is also the *only* external check the constructed targets have.  A
    /// constructed target carries no frozen CIR source, so two of the five
    /// production counters are vacuous for it, and the child Gamma star -- the
    /// only folded star the pinned identity-frequency table constrains -- always
    /// finds a stored row (every space group has a trivial Gamma irrep), so a
    /// constructed star never feeds that comparison.  Keep the counts below
    /// pinned, and see the R5 report's independence table.
    #[test]
    fn the_catalogue_reproduces_pinned_little_group_characters() {
        let mut records_used = 0usize;
        let mut operations_compared = 0usize;
        let mut projective = 0usize;
        let mut deferred = 0usize;
        for sg in 1u8..=230 {
            let cell = Lattice::new(exact_primitive_basis(sg).unwrap()).unwrap();
            let reciprocal = cell.reciprocal().unwrap();
            let hall = strict_sg_hall_ops(sg).unwrap();
            for record in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
                let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() else {
                    continue;
                };
                if row.dimension() != 1 {
                    continue;
                }
                let k = inline_k_vector(record).unwrap();
                let co_group = little_co_group(sg, &k, &reciprocal).unwrap();
                if co_group.order() <= 1 {
                    continue;
                }
                if co_group.order() > MAX_ORDER {
                    deferred += 1;
                    continue;
                }
                let characters = one_dimensional_characters(&co_group).unwrap();
                let operations: Vec<ExactSeitz> = hall
                    .operations
                    .iter()
                    .copied()
                    .filter(|operation| reciprocal.preserves(operation.rotation(), &k).unwrap())
                    .collect();
                let one_dimensional = characters.len() == co_group.order();
                let matched = if one_dimensional {
                    characters.iter().any(|psi| {
                        let constants = constants(&co_group, psi).unwrap();
                        operations.iter().enumerate().all(|(index, operation)| {
                            let pinned =
                                character_of(&row, record.ml, operation, &cell, &k, index).unwrap();
                            let expected =
                                character_value(&constants, &co_group.q, operation).unwrap();
                            (pinned - expected).norm() <= SUBDUCTION_TOLERANCE
                        })
                    })
                } else {
                    // R4 batch 2b: a higher-dimensional little irrep, compared
                    // the same way against the gated projective table.
                    let targets = projective_targets(sg, &k, &reciprocal).unwrap();
                    if targets.is_empty() {
                        deferred += 1;
                        continue;
                    }
                    let used = targets.iter().any(|target| {
                        operations.iter().enumerate().all(|(index, operation)| {
                            let pinned =
                                character_of(&row, record.ml, operation, &cell, &k, index).unwrap();
                            let expected = table_character_value(
                                &target.constants,
                                &co_group.q,
                                operation,
                            )
                            .unwrap();
                            (pinned - expected).norm() <= SUBDUCTION_TOLERANCE
                        })
                    });
                    projective += 1;
                    used
                };
                assert!(
                    matched,
                    "SG {sg} {} at k = {:?}: the pinned little-group character is not in the \
                     one-dimensional catalogue",
                    record.ml, k
                );
                records_used += 1;
                operations_compared += operations.len();
            }
        }
        // The counts are pinned, not bounded: this test is the **only** external
        // evidence the constructed targets have (they carry no frozen CIR source
        // and the pinned identity-frequency table never sees them), so a silent
        // shrinkage of the compared set must fail here rather than pass a
        // `>= 100` threshold.  The numbers are the ones reported in
        // `docs/subduction-r4-batches.md` and `docs/subduction-audit.md`.
        assert_eq!(records_used, 1328, "pinned rows compared");
        assert_eq!(operations_compared, 7578, "little-group operations compared");
        assert_eq!(projective, 94, "rows answered through the higher-dimensional table");
        assert_eq!(deferred, 1660, "rows deferred to the higher-dimensional batch");
        println!(
            "catalogue cross-check: {records_used} pinned rows, {operations_compared} operations, \
             {projective} of them through the higher-dimensional table, {deferred} deferred"
        );
    }
}
