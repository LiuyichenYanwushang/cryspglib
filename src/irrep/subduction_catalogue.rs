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
//! The caller uses this one-dimensional catalogue only when the solver returns
//! exactly `|P_q|` solutions.  The solution set is either empty or a coset of
//! `Hom(P_q, U(1))`, whose size is `|P_q / [P_q, P_q]|`; requiring that size to
//! be `|P_q|` therefore means `P_q` is abelian *and* the cocycle is a
//! coboundary, which is exactly the case where **every** irreducible projective
//! representation of `P_q` is one-dimensional.  A co-group whose cocycle is not
//! a coboundary, and a non-abelian co-group, return fewer solutions; the small
//! higher-dimensional families those cover are built separately by
//! [`projective_targets`], each behind its own finite structural gates.  A
//! co-group outside both paths returns an empty catalogue: this module never
//! guesses one.
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

/// Largest little co-group the exact solver searches.
///
/// R4 batch 2a needs the one-dimensional solution set for `|P_q| <= 4`; batch 2b
/// also needs it for the six-element `D_3` family, where it supplies the gauge
/// (its two solutions are `Hom(D_3, U(1))`) rather than a catalogue.  This batch
/// adds the eight-element `D4` family: the solver enumerates its four gauges and
/// [`projective_targets`] turns them into the four gauge-twisted one-dimensional
/// characters plus the gauge-twisted ordinary two-dimensional character.
/// Raising this bound also lets the generic all-one-dimensional path handle
/// order-eight abelian co-groups whenever their cocycle is a coboundary.
/// Anything bigger -- for example the order-48 cubic `Gamma` point, where the
/// generator candidates alone are far more numerous -- is out of scope and
/// returns an empty catalogue, so the caller keeps `MissingChildStarData`.
const MAX_ORDER: usize = 8;

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
/// Candidate phases are derived from each generator's finite order, so search
/// size depends on the co-group rather than the cocycle denominator. Every
/// returned vector is verified against every equation; unsupported group
/// orders and inconsistent systems return no catalogue rather than a subset.
pub(super) fn one_dimensional_characters(
    co_group: &LittleCoGroup,
) -> Result<Vec<Vec<Rat>>, StarError> {
    let order = co_group.order();
    if order > MAX_ORDER {
        return Ok(Vec::new());
    }
    let identity = 0usize; // `little_co_group` puts the identity first.
    let generators = co_group.generators()?;
    let choices = generators
        .iter()
        .map(|generator| co_group.generator_values(*generator, identity))
        .collect::<Result<Vec<_>, _>>()?;
    let mut solutions = Vec::new();
    let mut odometer = vec![0usize; generators.len()];
    loop {
        let assignment: Vec<Rat> = choices
            .iter()
            .zip(&odometer)
            .map(|(values, index)| values[*index])
            .collect();
        if let Some(psi) = co_group.complete(&generators, &assignment, identity)?
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
            if odometer[position] < choices[position].len() {
                break;
            }
            odometer[position] = 0;
            position += 1;
        }
    }
}

impl LittleCoGroup {
    /// All possible values of a projective character on one generator.
    ///
    /// If `g` has order `m`, summing the character equations for
    /// `e, g, ..., g^(m-1)` multiplied by `g` gives
    /// `m * psi(g) = sum_r omega(g^r, g) (mod 1)`. Thus there are exactly `m`
    /// candidate roots, independent of the cocycle denominator.
    fn generator_values(
        &self,
        generator: usize,
        identity: usize,
    ) -> Result<Vec<Rat>, StarError> {
        let mut power = identity;
        let mut order = 0usize;
        let mut phase = Rat::ZERO;
        loop {
            phase = fractional(phase.checked_add(self.turns[power][generator])?)?;
            power = self.product_position(power, generator)?;
            order += 1;
            if power == identity {
                break;
            }
            if order >= self.order() {
                return Err(StarError::LittleCoGroupNotClosed {
                    q: [self.q.get(0), self.q.get(1), self.q.get(2)],
                });
            }
        }
        let generator_order = order;
        let order = i128::try_from(generator_order).map_err(|_| {
            StarError::Subduction(SubductionError::RationalOverflow {
                operation: "little co-group generator order",
            })
        })?;
        let denominator = Rat::from_integer(order);
        let mut values = Vec::with_capacity(generator_order);
        for branch in 0..generator_order {
            let branch = i128::try_from(branch).map_err(|_| {
                StarError::Subduction(SubductionError::RationalOverflow {
                    operation: "little co-group generator root",
                })
            })?;
            let value = phase
                .checked_add(Rat::from_integer(branch))?
                .checked_div(denominator)?;
            let value = fractional(value)?;
            if !values.contains(&value) {
                values.push(value);
            }
        }
        Ok(values)
    }

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
        assignment: &[Rat],
        identity: usize,
    ) -> Result<Option<Vec<Rat>>, StarError> {
        let order = self.order();
        let mut psi: Vec<Option<Rat>> = vec![None; order];
        psi[identity] = Some(Rat::ZERO);
        for (generator, value) in generators.iter().zip(assignment) {
            psi[*generator] = Some(*value);
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
/// Three higher-dimensional families are covered, each behind its own
/// structural gate, and anything else returns an empty table (the caller keeps
/// `MissingChildStarData`):
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
/// * `|P_q| = 8`, an exact `D4` multiplication table with a **coboundary**
///   cocycle: the four gauges are `Hom(D4, U(1)) = D4/[D4,D4] = C2 x C2`, and
///   the projective irreps are the gauged ordinary ones, with dimensions
///   `{1, 1, 1, 1, 2}`.  Recognition and the character construction are finite
///   gates on the multiplication table; see [`d4_central_involution`] and
///   [`dihedral_eight_targets`].  A non-coboundary `D4` factor system -- for
///   example the section cocycle of the non-split extension
///   `1 -> C2 -> D16 -> D4 -> 1` -- has no one-dimensional solution and is
///   rejected rather than approximated by the coboundary table.
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
    if order == 8 {
        return dihedral_eight_targets(&co_group, identity, &q);
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

/// Recognize an exact `D4` multiplication table and return its central
/// involution, or `None` when the table is not `D4`.
///
/// Both gates are finite and structural -- no labels, no characters:
///
/// * the element-order profile is exactly `(1, 5, 2)`: one identity, five
///   involutions and two elements of order four.  Among the groups of order
///   eight this profile is unique to `D4`: `C2^3` has seven involutions,
///   `C4 x C2` three, `Q8` one, and `C8` has four elements of order eight.
/// * a generating relation of the standard presentation actually holds: some
///   order-four `r` and involution `s` satisfy `s r s = r^-1`, and `{r, s}`
///   spans the whole table, so the group really is
///   `<r, s | r^4 = s^2 = 1, s r s = r^-1>`.
///
/// The central involution is `r^2`; it is also checked to be the **unique**
/// involution commuting with every element, which is the element the ordinary
/// two-dimensional character is pinned on.
fn d4_central_involution(
    co_group: &LittleCoGroup,
    identity: usize,
) -> Result<Option<usize>, StarError> {
    let order = co_group.order();
    if order != 8 {
        return Ok(None);
    }
    let mut orders = Vec::with_capacity(order);
    for i in 0..order {
        let mut power = i;
        let mut value = 1usize;
        while power != identity {
            power = co_group.product_position(power, i)?;
            value += 1;
            if value > order {
                return Ok(None);
            }
        }
        orders.push(value);
    }
    if orders[identity] != 1
        || orders.iter().filter(|value| **value == 2).count() != 5
        || orders.iter().filter(|value| **value == 4).count() != 2
    {
        return Ok(None);
    }
    let mut central = None;
    'rotation: for rotation in 0..order {
        if orders[rotation] != 4 {
            continue;
        }
        let squared = co_group.product_position(rotation, rotation)?;
        let inverse = co_group.product_position(squared, rotation)?;
        for (reflection, &reflection_order) in orders.iter().enumerate() {
            if reflection_order != 2 {
                continue;
            }
            let conjugate = co_group
                .product_position(co_group.product_position(reflection, rotation)?, reflection)?;
            if conjugate == inverse
                && co_group
                    .close_span(&[identity], &[rotation, reflection])?
                    .len()
                    == order
            {
                central = Some(squared);
                break 'rotation;
            }
        }
    }
    let Some(central) = central else {
        return Ok(None);
    };
    if orders[central] != 2 {
        return Ok(None);
    }
    for (i, &element_order) in orders.iter().enumerate() {
        if element_order != 2 {
            continue;
        }
        let mut is_central = true;
        for j in 0..order {
            if co_group.product_position(i, j)? != co_group.product_position(j, i)? {
                is_central = false;
                break;
            }
        }
        if is_central != (i == central) {
            return Ok(None);
        }
    }
    Ok(Some(central))
}

/// The eight-element dihedral family (see [`projective_targets`]).
///
/// **Finite proof gates**, all on the multiplication table:
///
/// * recognition: [`d4_central_involution`] proves the table is `D4` and names
///   its unique central involution;
/// * coboundary: `one_dimensional_characters` returns exactly
///   `|Hom(D4, U(1))| = |D4/[D4,D4]| = 4` solutions.  A non-coboundary factor
///   system has no one-dimensional solution at all, and a co-group with the
///   `D4` profile but only some gauges cannot exist, so `!= 4` is the whole
///   coboundary gate.  The `D16` section-cocycle test drives the rejection.
/// * characters: the four one-dimensional irreps are `exp(2 pi i psi)` for the
///   four gauges, and the two-dimensional irrep is `exp(2 pi i psi_0)` times
///   the ordinary standard character.  In the standard realization on `R^2`,
///   `r -> [[0,-1],[1,0]]` and `s -> diag(1,-1)`, the ordinary traces are `2`
///   at the identity, `-2` at the central involution `r^2` and `0` on the other
///   six elements.  The gauge factors out because `omega = delta psi` makes
///   `g -> exp(-2 pi i psi(g)) u_g` an ordinary representation, hence
///   `tr(u_g) = exp(2 pi i psi(g)) tr(v_g)`.  The four gauges differ by the
///   ordinary linear characters (`C2 x C2`), which are trivial on the centre, so
///   the two-dimensional trace at the central involution `z` is the
///   gauge-independent `-2 exp(2 pi i psi(z))`, the same square root of
///   `omega(z,z)` in every gauge.  The phase is retained without assuming a
///   particular value for `omega(z,z)`.  The scalar equation
///   `chi(g) chi(h) = omega(g,h) chi(gh)` is deliberately never used: it is
///   false for the trace of a two-dimensional representation.
/// * completeness and orthonormality: the five rows are pairwise orthonormal
///   under `(1/|P|) sum_g chi_a(g) conj(chi_b(g))`, and their squared dimensions
///   sum to `|P| = 8`.  This is a **trace-level Gram matrix**; it is the finite
///   gate that replaces any use of the false scalar equation above, and it also
///   shows the five rows exhaust the projective irreps.
fn dihedral_eight_targets(
    co_group: &LittleCoGroup,
    identity: usize,
    q: &Vec3R,
) -> Result<Vec<ProjectiveTarget>, StarError> {
    let order = co_group.order();
    let Some(central) = d4_central_involution(co_group, identity)? else {
        return Ok(Vec::new());
    };
    // Coboundary gate: the gauge coset has exactly |Hom(D4, U(1))| = 4 members.
    let gauge = one_dimensional_characters(co_group)?;
    if gauge.len() != 4 {
        return Ok(Vec::new());
    }
    let mut rows: Vec<Vec<num_complex::Complex64>> = Vec::with_capacity(5);
    for psi in &gauge {
        rows.push(
            psi.iter()
                .map(|value| {
                    num_complex::Complex64::from_polar(1.0, std::f64::consts::TAU * value.to_f64())
                })
                .collect(),
        );
    }
    rows.push(
        gauge[0]
            .iter()
            .enumerate()
            .map(|(position, psi)| {
                let standard = if position == identity {
                    2.0
                } else if position == central {
                    -2.0
                } else {
                    0.0
                };
                num_complex::Complex64::from_polar(1.0, std::f64::consts::TAU * psi.to_f64())
                    * standard
            })
            .collect(),
    );
    // Trace-level orthonormality of the whole five-row table.
    let scale = 1.0
        / f64::from(u32::try_from(order).map_err(|_| {
            StarError::Subduction(SubductionError::RationalOverflow {
                operation: "little co-group order",
            })
        })?);
    for (a, left) in rows.iter().enumerate() {
        for (b, right) in rows.iter().enumerate() {
            let inner: num_complex::Complex64 = left
                .iter()
                .zip(right)
                .map(|(x, y)| x * y.conj())
                .sum::<num_complex::Complex64>()
                * scale;
            let expected = if a == b { 1.0 } else { 0.0 };
            if (inner - num_complex::Complex64::new(expected, 0.0)).norm() > 1e-9 {
                return Ok(Vec::new());
            }
        }
    }
    let mut targets = Vec::with_capacity(5);
    for (index, row) in rows.iter().enumerate() {
        let dimension = if index + 1 == rows.len() { 2 } else { 1 };
        targets.push(ProjectiveTarget {
            dimension,
            constants: constants_of(co_group, q, row)?,
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

    /// The C2 x C2 and D3 examples from batch 2b, described by their structure
    /// rather than by a pinned row. The separate D4 regression below pins the
    /// order-eight coboundary family.
    #[test]
    fn the_projective_tables_cover_only_the_gated_families() {
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

    /// Multiplicative orders of the co-group elements, indexed by position; the
    /// identity is at position zero for every co-group built here.
    fn element_orders(co_group: &LittleCoGroup) -> Vec<usize> {
        (0..co_group.order())
            .map(|element| {
                let mut power = element;
                let mut value = 1usize;
                while power != 0 {
                    power = co_group.product_position(power, element).unwrap();
                    value += 1;
                }
                value
            })
            .collect()
    }

    /// The classification's order-eight `D4` witness: child #136 (`P4_2/mnm`)
    /// at `q = (0, 0, 1/3)`, a point on the four-fold axis whose little
    /// co-group is `4mm = D4`.  The cocycle is a coboundary there, so the
    /// projective catalogue is the four gauged ordinary one-dimensional
    /// characters plus the gauged ordinary two-dimensional one.
    #[test]
    fn the_order_eight_d4_catalogue_is_four_gauges_and_the_standard_rep() {
        let cell = Lattice::new(exact_primitive_basis(136).unwrap()).unwrap();
        let reciprocal = cell.reciprocal().unwrap();
        let q = Vec3R::new([Rat::ZERO, Rat::ZERO, Rat::new(1, 3).unwrap()]);
        let co_group = little_co_group(136, &q, &reciprocal).unwrap();
        assert_eq!(
            co_group.order(),
            8,
            "the 4mm little co-group has order eight"
        );
        let orders = element_orders(&co_group);
        let profile: Vec<usize> = (1..=4)
            .map(|value| orders.iter().filter(|order| **order == value).count())
            .collect();
        assert_eq!(
            profile,
            vec![1, 5, 0, 2],
            "the D4 element-order profile: one identity, five involutions, no element of \
             order three, two elements of order four"
        );
        let central = d4_central_involution(&co_group, 0)
            .unwrap()
            .expect("the (1, 5, 2) table passes the D4 generating-relation gate");
        assert_eq!(
            orders[central], 2,
            "the central involution is an involution"
        );

        let gauge = one_dimensional_characters(&co_group).unwrap();
        assert_eq!(gauge.len(), 4, "|Hom(D4, U(1))| = |D4/[D4,D4]| = 4");

        let targets = projective_targets(136, &q, &reciprocal).unwrap();
        let dimensions: Vec<u8> = targets.iter().map(|target| target.dimension).collect();
        assert_eq!(dimensions, vec![1, 1, 1, 1, 2]);
        let total: u32 = targets
            .iter()
            .map(|target| u32::from(target.dimension) * u32::from(target.dimension))
            .sum();
        assert_eq!(total, 8, "sum of squared dimensions = |P_q|");

        let rows: Vec<Vec<Complex64>> = targets
            .iter()
            .map(|target| character_vector(target, &co_group))
            .collect();
        for row in &rows[..4] {
            for value in row {
                assert!(
                    (value.norm() - 1.0).abs() <= SUBDUCTION_TOLERANCE,
                    "a one-dimensional gauge row must have unit modulus: {value}"
                );
            }
        }
        for (a, left) in rows.iter().enumerate() {
            for (b, right) in rows.iter().enumerate() {
                let inner: Complex64 = left
                    .iter()
                    .zip(right)
                    .map(|(x, y)| x * y.conj())
                    .sum::<Complex64>()
                    / 8.0;
                let expected = if a == b { 1.0 } else { 0.0 };
                assert!(
                    (inner - Complex64::new(expected, 0.0)).norm() <= SUBDUCTION_TOLERANCE,
                    "rows {a} and {b} are not orthonormal: {inner}"
                );
            }
        }
        // The two-dimensional row is the gauge-twisted ordinary standard
        // character.  Its ordinary traces are fixed without a character table by
        // the standard realization on `R^2`, `r -> [[0,-1],[1,0]]` and
        // `s -> diag(1,-1)`: trace `2` at the identity, `-2` at `r^2` (the
        // unique central involution) and `0` on the other six elements.  The
        // gauge multiplies every trace by a unit phase (see
        // `dihedral_eight_targets`); at this witness the cocycle is trivial on
        // the centre and the solver's first gauge takes the value `0` there, so
        // the returned row is exactly the ordinary character.  (The scalar
        // equation `chi(g) chi(h) = omega(g,h) chi(gh)` is never used: it is
        // false for a two-dimensional trace.)
        for (position, value) in rows[4].iter().enumerate() {
            let expected = if position == 0 {
                Complex64::new(2.0, 0.0)
            } else if position == central {
                Complex64::new(-2.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            assert!(
                (value - expected).norm() <= SUBDUCTION_TOLERANCE,
                "the standard D4 row at position {position}: {value} != {expected}"
            );
        }
    }

    /// The larger one-dimensional search bound also admits abelian order-eight
    /// co-groups. They stay on the generic scalar path and are not misidentified
    /// as D4 by the high-dimensional catalogue gate.
    #[test]
    fn an_abelian_order_eight_cogroup_has_eight_scalar_characters() {
        let zero = Vec3R::zero();
        let quarter_turn =
            ExactSeitz::new([[0, -1, 0], [1, 0, 0], [0, 0, 1]], zero);
        let inversion = ExactSeitz::new([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], zero);
        let mut representatives = Vec::with_capacity(8);
        let mut power = ExactSeitz::identity();
        for _ in 0..4 {
            representatives.push(power);
            representatives.push(power.compose(&inversion).unwrap());
            power = power.compose(&quarter_turn).unwrap();
        }
        let co_group = LittleCoGroup {
            turns: vec![vec![Rat::ZERO; 8]; 8],
            representatives,
            q: zero,
        };
        assert_eq!(co_group.order(), 8);
        assert_eq!(
            one_dimensional_characters(&co_group).unwrap().len(),
            8,
            "C4 x C2 has eight ordinary scalar irreps"
        );
        assert!(
            d4_central_involution(&co_group, 0).unwrap().is_none(),
            "the C4 x C2 element-order profile is not D4"
        );
    }

    /// The `D16` section cocycle: a genuine **non-coboundary** `D4` factor
    /// system.  The test refuses to accept the coboundary `D4` table for it.
    ///
    /// `D16 = <R, S | R^8 = S^2 = 1, S R S = R^-1>` has centre `<R^4>`, and its
    /// quotient by that centre is
    /// `D4 = <r, s | r^4 = s^2 = 1, s r s = r^-1>` (`r = R Z`, `s = S Z`).  The
    /// set-theoretic section `r^a s^b -> R^a S^b` (`a` in `0..4`, `b` in
    /// `0..2`) is not a homomorphism, because `R^4 != 1`; its defect
    /// `sigma(g) sigma(h) sigma(gh)^-1` lands in the centre, so the factor
    /// system takes only the values `0` and `1/2`.  A gauge `psi` with
    /// `delta psi = omega` would make `g -> exp(-2 pi i psi(g)) sigma(g)` a
    /// homomorphism and would split the extension; `D16` is not `C2 x D4` (its
    /// centre has order two, while `C2 x D4` has three central involutions), so
    /// no gauge exists.  The cocycle identity is additionally verified by brute
    /// force below, so the negative case is a valid factor system and not a
    /// single-entry mutation.
    #[test]
    fn a_non_coboundary_d4_factor_system_returns_no_catalogue() {
        let co_group = d16_section_cocycle_on_d4();
        assert_eq!(co_group.order(), 8);
        // The recognition gates pass: the table really is D4, so the emptiness
        // asserted below can only come from the coboundary gate.
        assert_eq!(
            element_orders(&co_group)
                .iter()
                .filter(|value| **value == 2)
                .count(),
            5
        );
        assert!(
            d4_central_involution(&co_group, 0).unwrap().is_some(),
            "the D16 section cocycle lives on an exact D4 multiplication table"
        );
        // Valid factor system: the cocycle identity
        // `omega(i,j) omega(ij,k) == omega(j,k) omega(i,jk)` on every triple.
        for i in 0..8 {
            for j in 0..8 {
                for k in 0..8 {
                    let ij = co_group.product_position(i, j).unwrap();
                    let jk = co_group.product_position(j, k).unwrap();
                    let left = fractional(
                        co_group.turns[i][j]
                            .checked_add(co_group.turns[ij][k])
                            .unwrap(),
                    )
                    .unwrap();
                    let right = fractional(
                        co_group.turns[j][k]
                            .checked_add(co_group.turns[i][jk])
                            .unwrap(),
                    )
                    .unwrap();
                    assert_eq!(left, right, "cocycle identity at ({i}, {j}, {k})");
                }
            }
        }
        // Non-coboundary: no one-dimensional projective character, because a
        // gauge would split `D16 -> D4`.  The D4 branch must therefore return
        // nothing rather than the coboundary table.
        let gauge = one_dimensional_characters(&co_group).unwrap();
        assert!(
            gauge.is_empty(),
            "a non-coboundary D4 factor system has no gauge: {gauge:?}"
        );
        assert!(
            dihedral_eight_targets(&co_group, 0, &Vec3R::zero())
                .unwrap()
                .is_empty(),
            "the non-coboundary factor system must not return the D4 coboundary table"
        );
    }

    /// The eight symmetries of the square as `3 x 3` integer matrices, in the
    /// order `r^a s^b` (`a` in `0..4`, `b` in `0..2`), so position `2a + b` is
    /// the group element `r^a s^b` and position zero is the identity.
    fn d4_representatives() -> Vec<ExactSeitz> {
        let r = ExactSeitz::new([[0, -1, 0], [1, 0, 0], [0, 0, 1]], Vec3R::zero());
        let s = ExactSeitz::new([[1, 0, 0], [0, -1, 0], [0, 0, 1]], Vec3R::zero());
        let mut representatives = Vec::with_capacity(8);
        let mut power = ExactSeitz::identity();
        for _ in 0..4 {
            representatives.push(power);
            representatives.push(power.compose(&s).unwrap());
            power = power.compose(&r).unwrap();
        }
        representatives
    }

    /// The `D16 -> D4` section cocycle as a hand-built [`LittleCoGroup`].
    ///
    /// `D16` is modelled exactly by pairs `(a, b)` with `a` in `Z/8`, `b` in
    /// `Z/2` and `(a,b)(c,d) = (a + (-1)^b c, b + d)`; the section lifts the D4
    /// normal form `r^a s^b` (`a` in `0..4`) to `(a, b)`, and the defect of a
    /// product is either `(0, 0)` or the centre `(4, 0)`.
    fn d16_section_cocycle_on_d4() -> LittleCoGroup {
        let representatives = d4_representatives();
        let order = representatives.len();
        let mut turns = Vec::with_capacity(order);
        for i in 0..order {
            let mut row = Vec::with_capacity(order);
            for j in 0..order {
                let (a, b) = ((i / 2) as i32, (i % 2) as i32);
                let (c, d) = ((j / 2) as i32, (j % 2) as i32);
                let signed_c = if b == 0 { c } else { -c };
                // The D4 product and the D16 product of the two lifts.
                let product_d4 = ((a + signed_c).rem_euclid(4), (b + d) % 2);
                let lift = |x: i32, y: i32| (x.rem_euclid(8), y.rem_euclid(2));
                let (left_a, left_b) = {
                    let (a0, b0) = lift(a, b);
                    let (c0, d0) = lift(c, d);
                    let signed = if b0 == 0 { c0 } else { -c0 };
                    ((a0 + signed).rem_euclid(8), (b0 + d0).rem_euclid(2))
                };
                let (right_a, right_b) = lift(product_d4.0, product_d4.1);
                let defect = (
                    (left_a - right_a).rem_euclid(8),
                    (left_b - right_b).rem_euclid(2),
                );
                row.push(match defect {
                    (0, 0) => Rat::ZERO,
                    (4, 0) => Rat::new(1, 2).unwrap(),
                    _ => panic!("the D16 section defect left the centre: {defect:?}"),
                });
            }
            turns.push(row);
        }
        LittleCoGroup {
            representatives,
            turns,
            q: Vec3R::zero(),
        }
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
    }

    /// Generator-order enumeration is independent of the cocycle denominator.
    #[test]
    fn generator_order_search_handles_large_cocycle_denominators() {
        let mut scaled = LittleCoGroup {
            representatives: vec![
                ExactSeitz::identity(),
                ExactSeitz::new([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], Vec3R::zero()),
                ExactSeitz::new([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], Vec3R::zero()),
                ExactSeitz::new([[1, 0, 0], [0, -1, 0], [0, 0, -1]], Vec3R::zero()),
            ],
            turns: vec![vec![Rat::ZERO; 4]; 4],
            q: Vec3R::zero(),
        };

        // Build a genuine coboundary from a 1-cochain with denominator 512.
        // Its four characters are this cochain plus Hom(C2 x C2, U(1)).
        let expected = [
            Rat::ZERO,
            Rat::new(1, 512).unwrap(),
            Rat::new(2, 512).unwrap(),
            Rat::new(3, 512).unwrap(),
        ];
        for left in 0..4 {
            for right in 0..4 {
                let product = scaled.product_position(left, right).unwrap();
                scaled.turns[left][right] = fractional(
                    expected[left]
                        .checked_add(expected[right])
                        .unwrap()
                        .checked_sub(expected[product])
                        .unwrap(),
                )
                .unwrap();
            }
        }

        let characters = one_dimensional_characters(&scaled).unwrap();
        assert_eq!(characters.len(), 4, "the full Hom(C2 x C2,U(1)) coset");
        assert!(
            characters
                .iter()
                .any(|character| character.as_slice() == expected.as_slice()),
            "the cochain used to build the cocycle is one of its gauges"
        );
        for character in &characters {
            for left in 0..4 {
                for right in 0..4 {
                    let product = scaled.product_position(left, right).unwrap();
                    let residual = fractional(
                        character[left]
                            .checked_add(character[right])
                            .unwrap()
                            .checked_sub(character[product])
                            .unwrap()
                            .checked_sub(scaled.turns[left][right])
                            .unwrap(),
                    )
                    .unwrap();
                    assert!(residual.is_zero(), "{character:?} at ({left}, {right})");
                }
            }
        }
    }

    #[test]
    fn the_trivial_co_group_has_its_unique_projective_character() {
        let trivial = LittleCoGroup {
            representatives: vec![ExactSeitz::identity()],
            turns: vec![vec![Rat::ZERO]],
            q: Vec3R::zero(),
        };
        assert_eq!(
            one_dimensional_characters(&trivial).unwrap(),
            vec![vec![Rat::ZERO]]
        );
    }

    /// The decisive evidence for the catalogue: at every pinned child `k` whose
    /// little co-group is non-trivial and whose irreps are all one-dimensional,
    /// the engine's own validated little-group character row must be **one of**
    /// the catalogue characters, on every little-group operation.
    ///
    /// Rows of dimension two are compared too: a two-dimensional pinned little
    /// irrep at a co-group inside the gated families must be one of the
    /// constructed two-dimensional targets.  This is where the `D4` family is
    /// cross-checked against real data -- the `W5`/`X5`/`M5` rows of the
    /// tetragonal and cubic space groups are exactly the relevant dimension-2
    /// `D4` source rows.  (Two-dimensional rows whose co-group is a
    /// **non-coboundary** `D4`, e.g. SG 212/213 `X1`/`X2` or SG 227 `W1`/`W2`,
    /// are not matchable by construction: they are counted as deferred, never
    /// as matched, and the negative test below drives the rejection directly.)
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
        let mut d4_one_dimensional = 0usize;
        let mut two_dimensional = 0usize;
        let mut d4_two_dimensional = 0usize;
        for sg in 1u8..=230 {
            let cell = Lattice::new(exact_primitive_basis(sg).unwrap()).unwrap();
            let reciprocal = cell.reciprocal().unwrap();
            let hall = strict_sg_hall_ops(sg).unwrap();
            for record in query::irreps_of(sg).iter().filter(|record| !record.spinor) {
                let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() else {
                    continue;
                };
                let dimension = row.dimension();
                if dimension != 1 && dimension != 2 {
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
                    // Every projective irrep is one-dimensional here, so a
                    // dimension-two pinned row would contradict the pinned
                    // dimension; fail loudly instead of silently skipping.
                    assert_eq!(
                        dimension, 1,
                        "SG {sg} {} at k = {k:?}: a two-dimensional row at a co-group with only \
                         one-dimensional projective irreps",
                        record.ml
                    );
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
                    // R4 batch 2b and the `D4` batch: a higher-dimensional little
                    // irrep, compared the same way against the gated projective
                    // table.
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
                    // The full D4 table carries four scalar source rows and
                    // one two-dimensional source row per character set.
                    if co_group.order() == 8 && dimension == 1 {
                        d4_one_dimensional += 1;
                    }
                    if dimension == 2 {
                        two_dimensional += 1;
                        if co_group.order() == 8 {
                            // A non-empty order-eight table is the `D4` one (the
                            // recognition and coboundary gates decide it), so
                            // these are exactly the dimension-two `D4` source
                            // rows the batch was written for.
                            d4_two_dimensional += 1;
                        }
                    }
                    used
                };
                assert!(
                    matched,
                    "SG {sg} {} at k = {:?}: the pinned little-group character is not in the \
                     constructed projective catalogue",
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
        // `>= 100` threshold.
        //
        // Before the `D4` batch (one-dimensional solver gate `MAX_ORDER = 6`,
        // dimension-one rows only) the four counters were
        // `records_used = 1328`, `operations_compared = 7578`, `projective = 94`,
        // `deferred = 1660`.  After raising the gate to 8 and extending the
        // corpus to dimension-two rows they are the values pinned below:
        // `2661 / 22302 / 719 / 293 / 121 / 1282`, where the fifth counter is the
        // dimension-two `D4` subset of the fourth.  Order-eight rows moved out of
        // `deferred` (1660 -> 1282), while the newly considered dimension-two
        // rows with a co-group of order > 8 moved in; two-dimensional rows with a
        // non-coboundary `D4` co-group stay deferred because their table is
        // intentionally not constructed.
        println!(
            "catalogue cross-check: {records_used} pinned rows, {operations_compared} operations, \
             {projective} of them through the higher-dimensional table ({d4_one_dimensional} \
             D4 dimension-one, {two_dimensional} dimension-two, {d4_two_dimensional} of them D4), \
             {deferred} deferred"
        );
        assert_eq!(records_used, 2661, "pinned rows compared");
        assert_eq!(
            operations_compared, 22302,
            "little-group operations compared"
        );
        assert_eq!(
            projective, 719,
            "rows answered through the higher-dimensional table"
        );
        assert_eq!(d4_one_dimensional, 332, "one-dimensional D4 source rows");
        assert_eq!(two_dimensional, 293, "dimension-two rows so compared");
        assert_eq!(d4_two_dimensional, 121, "dimension-two D4 rows so compared");
        assert_eq!(
            deferred, 1282,
            "rows deferred to the higher-dimensional batch"
        );
    }
}
