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
//! solutions, which is the size of a coset of `Hom(P_q, U(1))` and therefore
//! exactly the case where **every** irreducible projective representation of
//! `P_q` is one-dimensional.  A co-group whose cocycle is not a coboundary
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
use super::super::strict_sg_hall_ops;
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

/// Largest little co-group the one-dimensional batch solves.
///
/// The R4 batch-2a groups have `|P_q| <= 4` (order 2: 452, order 3: 4, order 4
/// with four omega-regular classes: 170).  A bigger co-group is out of scope --
/// its irreps may need higher dimensions, and the search below would grow with
/// its order -- so it returns an empty catalogue and keeps the missing-data
/// behaviour.
const MAX_ORDER: usize = 4;

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

    /// The decisive evidence for the catalogue: at every pinned child `k` whose
    /// little co-group is non-trivial and whose irreps are all one-dimensional,
    /// the engine's own validated little-group character row must be **one of**
    /// the catalogue characters, on every little-group operation.
    #[test]
    fn the_catalogue_reproduces_pinned_little_group_characters() {
        let mut records_used = 0usize;
        let mut operations_compared = 0usize;
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
                if characters.len() != co_group.order() {
                    deferred += 1;
                    continue;
                }
                let operations: Vec<ExactSeitz> = hall
                    .operations
                    .iter()
                    .copied()
                    .filter(|operation| reciprocal.preserves(operation.rotation(), &k).unwrap())
                    .collect();
                let matched = characters.iter().any(|psi| {
                    let constants = constants(&co_group, psi).unwrap();
                    operations.iter().enumerate().all(|(index, operation)| {
                        let pinned =
                            character_of(&row, record.ml, operation, &cell, &k, index).unwrap();
                        let expected = character_value(&constants, &co_group.q, operation).unwrap();
                        (pinned - expected).norm() <= SUBDUCTION_TOLERANCE
                    })
                });
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
        assert!(
            records_used >= 100,
            "expected a large sample of non-trivial little co-groups, got {records_used}"
        );
        assert!(
            operations_compared >= 300,
            "expected hundreds of compared operations, got {operations_compared}"
        );
        println!(
            "catalogue cross-check: {records_used} pinned rows, {operations_compared} operations, \
             {deferred} deferred to the higher-dimensional batch"
        );
    }
}
