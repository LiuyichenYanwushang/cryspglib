//! Round 24 (R9 stage 2b, F1): the **complete** shift solution set of one fixed
//! candidate map, exactly -- no step grid, no sampling of the continuous kernel.
//!
//! The shipped `search_embeddings` solved the anchor equation of **one** class
//! exactly and then walked the affine set `particular + ker(I - R_p)` with the
//! step set `{+-1, +-1/2}` (round 15m).  Round 23's review measured two
//! consequences on the public API: the walk cannot reach a solution that differs
//! from the anchor's particular solution by a **lattice-shifted branch** (SG 56
//! UNI 119 needs `delta_1 = 3/4` where the exact equation gives `1/4`), and the
//! step set cannot close a component whose denominator is a quarter (SG 89
//! UNI 111).  Both are properties of the walk, not of the mathematics.  For a
//! fixed map the solution set of
//!
//! ```text
//! (I - R_p(c)) delta  ==  tau_p(c) - A tau_c     (mod L_P),  every class c
//! ```
//!
//! is a finite union of cosets of the common real kernel
//!
//! ```text
//! V = { v in R^3 : (I - R_p(c)) v = 0 for every class c },
//! ```
//!
//! and this module computes that union **exactly and completely**:
//!
//! * the parent lattice `L_P = Z^3 + span_Z(generators)` is made integral,
//!   `K = D L_P` with the `Z`-basis `B` of [`parent_lattice_basis`], and the
//!   substitution `delta = C y / D` with `C = B^T` (so `delta in L_P` iff
//!   `y in Z^3`) turns the class equations into integer congruences
//!   `P_c y == q_c (mod e Z^3)`, `e = det C`, `P_c = adj(C) M_c C`,
//!   `q_c = D adj(C) t_c`;
//! * the stacked system is decided by integer Euclid elimination plus a
//!   Smith-style diagonalization -- **unimodular** row and column operations
//!   only, so the modulus lattice `e Z^m` and the domain `R^3/Z^3` are preserved
//!   exactly;
//! * the solution set is reported as one representative per discrete component
//!   plus a basis of `V`: `S = union over branches of (branch + V)`, in the
//!   parent frame, reduced modulo `L_P`.  The continuous directions are never
//!   enumerated, and the branches are exactly the components (the module docs of
//!   [`solve_joint_congruence`] prove the mixed-radix count).
//!
//! # Contract
//!
//! * [`solve_shift_components`] returns `Ok` with an **empty** representative
//!   list exactly when the congruence system has no solution: a proof, not "the
//!   search did not find one".  The kernel may still be reported in that case
//!   (it is a property of the coefficient matrices alone).
//! * `Err(Arithmetic)` (an `i128` limit), `Err(Capacity)` ([`MAX_COMPONENTS`] or
//!   [`MAX_STEPS`]) and `Err(NotWellDefined)` (a synthetic caller whose
//!   congruence does not descend to `R^3/Z^3`) mean the answer is **unknown**; a
//!   caller that treats them as "no solution" would be wrong, and the production
//!   search counts them as unsettled instead (round 23 review F2's rule).  The
//!   production reduction always satisfies the precondition -- see
//!   [`solve_joint_congruence`].
//! * The component *set* is canonical.  The representative of a component is the
//!   elimination's own point, reduced modulo `L_P` by [`canonical_translation_checked`]
//!   (the range-free lexicographic minimum); its class modulo `V + L_P` is what
//!   the mathematics pins, and the list is sorted so the output is a
//!   deterministic function of the input.

use super::{
    MagneticContractError, Mat3I, Rat, ShiftComponents, canonical_translation_checked, gcd,
    parent_lattice_basis,
};

/// Largest number of discrete components [`solve_joint_congruence`] will
/// materialize.  Each component contributes one exact point; the cap turns a
/// pathological diagonal (say `d = (100000, 0, 0)`) into an explicit unknown
/// answer instead of an allocation by the million.
const MAX_COMPONENTS: i128 = 1 << 16;

/// Largest number of elementary reduction steps in one call.  Every loop here
/// provably terminates (each pass strictly lowers the smallest non-zero
/// magnitude of the current submatrix); the cap only turns a hypothetical
/// implementation bug into an error instead of a hang.
const MAX_STEPS: u64 = 1_000_000;

/// Why the complete shift solution could not be computed.  Infeasibility is
/// **not** reported here: it is `Ok` with no representatives.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ShiftEnumerationError {
    /// An `i128` product or sum left the representable range.
    Arithmetic,
    /// [`MAX_COMPONENTS`] or [`MAX_STEPS`] was exceeded.
    Capacity,
    /// The row count and the number of right-hand sides differ.  Reported
    /// instead of silently truncating the system with `zip`: the two slices are
    /// built together by the caller, so this is an internal contract violation.
    LengthMismatch { rows: usize, rhs: usize },
    /// The congruence does not descend to `R^3/Z^3`: some coefficient row has an
    /// entry that `modulus` does not divide, so `rows * Z^3` is not inside
    /// `modulus Z^3` and the solution set is not a union of classes modulo the
    /// cell translations.  The production caller never hits this -- its rows are
    /// `adj(C) M C` with `M = I - R` and `R` a parent rotation, which preserves
    /// the parent lattice, so every entry is divisible by `det C` -- but a
    /// synthetic caller can, and then the honest answer is "unknown" rather than
    /// a branch list that is not a solution set.
    NotWellDefined,
}

/// The solution of one integer congruence system `rows x == rhs (mod modulus)`
/// for `x in R^3/Z^3`, as a finite description.
#[derive(Clone, Debug, PartialEq, Eq)]
struct JointSolution {
    /// One point per discrete component, reduced modulo `Z^3`.
    branches: Vec<[Rat; 3]>,
    /// Integer basis of the common real kernel; every branch is a coset of its
    /// span modulo `Z^3`.
    kernel: Vec<[i128; 3]>,
}

/// Solve the class equations of one fixed map completely.
///
/// `equations` holds `(M_c, t_c)` for every class, with `M_c = I - R_p(c)`
/// integer and `t_c = tau_p(c) - A tau_c` rational.  `parent_lattice` lists the
/// parent's non-zero pure translations; the cell translations `Z^3` are implicit
/// and are added by [`parent_lattice_basis`], exactly as in
/// [`super::canonical_translation`].
pub(crate) fn solve_shift_components(
    equations: &[(Mat3I, [Rat; 3])],
    parent_lattice: &[[Rat; 3]],
) -> Result<ShiftComponents, ShiftEnumerationError> {
    let (denominator, basis) =
        parent_lattice_basis(parent_lattice).map_err(|_| ShiftEnumerationError::Arithmetic)?;
    let change = transpose(&basis);
    let modulus = determinant(&change)?
        .checked_abs()
        .ok_or(ShiftEnumerationError::Arithmetic)?;
    if modulus == 0 {
        return Err(ShiftEnumerationError::Arithmetic);
    }
    let adjugate = adjugate(&change)?;
    let mut rows: Vec<[i128; 3]> = Vec::with_capacity(3 * equations.len());
    let mut rhs: Vec<Rat> = Vec::with_capacity(3 * equations.len());
    for (matrix, target) in equations {
        let matrix = ints_to_wide(matrix);
        let product = multiply(&multiply(&adjugate, &matrix)?, &change)?;
        let moved = mul_vector(&adjugate, target)?;
        for row in 0..3 {
            rows.push(product[row]);
            rhs.push(
                moved[row]
                    .checked_mul(Rat::from_integer(denominator))
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?,
            );
        }
    }
    let solution = match solve_joint_congruence(&rows, &rhs, modulus)? {
        Some(solution) => solution,
        None => {
            return Ok(ShiftComponents {
                representatives: Vec::new(),
                kernel: Vec::new(),
            });
        }
    };
    let mut representatives: Vec<[Rat; 3]> = Vec::with_capacity(solution.branches.len());
    for branch in &solution.branches {
        let mut point = [Rat::ZERO; 3];
        for row in 0..3 {
            let mut sum = Rat::ZERO;
            for (column, value) in branch.iter().enumerate() {
                let term = value
                    .checked_mul(Rat::from_integer(change[row][column]))
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
                sum = sum
                    .checked_add(term)
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
            }
            point[row] = sum;
        }
        for value in point.iter_mut() {
            *value = value
                .checked_div(Rat::from_integer(denominator))
                .map_err(|_| ShiftEnumerationError::Arithmetic)?;
        }
        representatives.push(to_parent_frame(point, parent_lattice)?);
    }
    // `Rat` is not `Ord` (it is exhaustive over a sign-normalized fraction), so
    // the deterministic order is by the reduced numerator/denominator pairs.
    representatives
        .sort_by_key(|point| point.map(|value| (value.numerator(), value.denominator())));
    representatives.dedup();
    let mut kernel: Vec<[Rat; 3]> = Vec::with_capacity(solution.kernel.len());
    for direction in &solution.kernel {
        let mut point = [Rat::ZERO; 3];
        for row in 0..3 {
            let mut sum = Rat::ZERO;
            for (column, value) in direction.iter().enumerate() {
                let term = Rat::from_integer(change[row][column])
                    .checked_mul(Rat::from_integer(*value))
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
                sum = sum
                    .checked_add(term)
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
            }
            point[row] = sum;
        }
        for value in point.iter_mut() {
            *value = value
                .checked_div(Rat::from_integer(denominator))
                .map_err(|_| ShiftEnumerationError::Arithmetic)?;
        }
        kernel.push(point);
    }
    Ok(ShiftComponents {
        representatives,
        kernel,
    })
}

/// Reduce a shift into the parent frame's canonical class representative.
fn to_parent_frame(
    point: [Rat; 3],
    parent_lattice: &[[Rat; 3]],
) -> Result<[Rat; 3], ShiftEnumerationError> {
    canonical_translation_checked(point, parent_lattice).map_err(|error| match error {
        MagneticContractError::ReductionCapacity { .. } => ShiftEnumerationError::Capacity,
        _ => ShiftEnumerationError::Arithmetic,
    })
}
/// Solve `rows x == rhs (mod modulus)` for `x in R^3/Z^3`, exactly and
/// completely.
///
/// **Preconditions.**  `rows` and `rhs` must have the same length
/// ([`ShiftEnumerationError::LengthMismatch`] otherwise, instead of a silent
/// `zip` truncation), and the congruence must descend to the torus: `modulus` has
/// to divide every entry of `rows`, so that `rows * Z^3 <= modulus Z^m` and
/// adding a cell translation to `x` cannot change the left-hand side.  A
/// violated precondition is [`ShiftEnumerationError::NotWellDefined`], never a
/// branch list: without it the "solution set modulo `Z^3`" is not a solution set at all
/// (the round-24 test sweep caught exactly this with `modulus = 2` and odd
/// coefficients).  [`solve_shift_components`] satisfies it by construction: its
/// rows are `adj(C) M C` with `M = I - R` and `R` a rotation of the parent, and
/// `R` preserves the parent lattice `C Z^3 / D`, so `M C Z^3 <= C Z^3` and
/// `adj(C) M C Z^3 <= adj(C) C Z^3 = det(C) Z^3`.
///
/// `Ok(None)` means **proved infeasible**.  `Ok(Some(solution))` describes the
/// whole solution set `union over branches of (branch + span(kernel))` modulo
/// `Z^3`, and the correctness argument is:
///
/// 1. **Row reduction.**  Swapping rows, negating a row, or `row_i += q row_j`
///    (`q in Z`) is multiplication by a unimodular integer matrix; because
///    `U (modulus Z^m) = modulus Z^m`, the system is preserved exactly, and the
///    zero rows that remain read `0 == rhs_i (mod modulus)`.  If any of those
///    fails, the system is inconsistent -- and since every later operation is
///    also invertible over `Q`, that is a proof of infeasibility.
/// 2. **Diagonalization.**  Unimodular row and column operations on the `3 x 3`
///    pivot block (padding rows are zero) give `D = U P V` diagonal, the row
///    operations carried into the right-hand side and the column change
///    accumulated in `basis` (`P V = U^-1 D`, and `x = basis z`, so the system
///    becomes `D z == U rhs (mod modulus Z^3)`).
/// 3. **Decoupling.**  A row with `d_i != 0` reads `d_i z_i == beta_i
///    (mod modulus)`; `z_i = (beta_i + modulus t) / d_i` is a real solution for
///    **every** integer `t`, so a non-zero diagonal is never infeasible (round
///    24 corrected the first draft, which demanded `gcd(d_i, modulus) | beta_i`
///    -- that is the condition for an *integer* solution and wrongly refuted
///    SG 56 UNI 119).  Two choices of `t` give the same class modulo `Z` exactly
///    when `d_i | modulus (t - t')`, i.e. when `t` agrees modulo
///    `d_i / gcd(d_i, modulus)`: that quotient -- call it `c_i` -- is the number
///    of distinct classes, and the residues `t = 0, .., c_i - 1` enumerate them.
/// 4. **Free directions.**  A row with `d_i == 0` reads `0 == beta_i
///    (mod modulus)`, already checked as a zero row, and `z_i` is free; the
///    corresponding column of `basis` spans the common real kernel, because
///    `P basis = U^-1 D` annihilates exactly the zero-diagonal columns of `D`.
/// 5. **Completeness.**  Every solution is `basis z` with `z_i` in the coset of
///    step 3 (pivot rows) or real (zero rows), which is exactly one of the
///    enumerated branches plus a kernel element.  The branches are pairwise
///    distinct modulo `Z^3`: two of them differ by `sum_i (modulus Δt_i/d_i) v_i`
///    with `0 < |Δt_i| < c_i`, and `modulus Δt_i / d_i` is not an integer there
///    (its denominator in lowest terms is `d_i/gcd(d_i, modulus)` divided by
///    `gcd(Δt_i, ...)`, which is `> 1` because `gcd(d_i/g, Δt_i) < d_i/g`), so
///    the difference is not in `basis Z^3 = Z^3`.
fn solve_joint_congruence(
    rows: &[[i128; 3]],
    rhs: &[Rat],
    modulus: i128,
) -> Result<Option<JointSolution>, ShiftEnumerationError> {
    if rows.len() != rhs.len() {
        return Err(ShiftEnumerationError::LengthMismatch {
            rows: rows.len(),
            rhs: rhs.len(),
        });
    }
    let modulus = modulus
        .checked_abs()
        .ok_or(ShiftEnumerationError::Arithmetic)?;
    if modulus == 0 {
        return Err(ShiftEnumerationError::Arithmetic);
    }
    // The domain is `R^3/Z^3`, so the congruence has to be well defined there:
    // adding a cell translation to `x` must not change the left-hand side
    // modulo `modulus Z^m`, which is exactly `rows * Z^3 <= modulus Z^m`, i.e.
    // `modulus` divides every entry.
    for row in rows {
        for entry in row {
            if entry % modulus != 0 {
                return Err(ShiftEnumerationError::NotWellDefined);
            }
        }
    }
    let mut budget = Budget::default();
    let mut work: Vec<([i128; 3], Rat)> = rows.iter().copied().zip(rhs.iter().copied()).collect();
    let rank = row_echelon(&mut work, &mut budget)?;
    for (_, value) in work.iter().skip(rank) {
        if !is_multiple_of_integer(value, modulus) {
            return Ok(None);
        }
    }
    let mut matrix = [[0i128; 3]; 3];
    let mut target = [Rat::ZERO; 3];
    for index in 0..rank {
        matrix[index] = work[index].0;
        target[index] = work[index].1;
    }
    let mut basis = [[1i128, 0, 0], [0, 1, 0], [0, 0, 1]];
    diagonalize(&mut matrix, &mut target, &mut basis, &mut budget)?;
    let mut sizes = [0i128; 3];
    let mut total: i128 = 1;
    for index in 0..3 {
        let value = matrix[index][index];
        if value == 0 {
            continue;
        }
        // `i128::MIN` has no positive counterpart: the magnitude itself is
        // needed here (it multiplies into the component count), so the honest
        // answer is an arithmetic limit.
        let size = value
            .checked_abs()
            .ok_or(ShiftEnumerationError::Arithmetic)?;
        sizes[index] = size / gcd(size, modulus);
        total = total
            .checked_mul(sizes[index])
            .ok_or(ShiftEnumerationError::Arithmetic)?;
        if total > MAX_COMPONENTS {
            return Err(ShiftEnumerationError::Capacity);
        }
    }
    let capacity = usize::try_from(total).map_err(|_| ShiftEnumerationError::Capacity)?;
    let mut branches = Vec::with_capacity(capacity);
    let mut residues = [0i128; 3];
    for _ in 0..total {
        let mut free = [Rat::ZERO; 3];
        for index in 0..3 {
            if sizes[index] == 0 {
                continue;
            }
            let shift = residues[index]
                .checked_mul(modulus)
                .ok_or(ShiftEnumerationError::Arithmetic)?;
            let numerator = target[index]
                .checked_add(Rat::from_integer(shift))
                .map_err(|_| ShiftEnumerationError::Arithmetic)?;
            free[index] = numerator
                .checked_div(Rat::from_integer(matrix[index][index]))
                .map_err(|_| ShiftEnumerationError::Arithmetic)?;
        }
        let mut image = [Rat::ZERO; 3];
        for row in 0..3 {
            let mut sum = Rat::ZERO;
            for column in 0..3 {
                if basis[row][column] == 0 {
                    continue;
                }
                let term = free[column]
                    .checked_mul(Rat::from_integer(basis[row][column]))
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
                sum = sum
                    .checked_add(term)
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
            }
            image[row] = fraction_part(sum)?;
        }
        branches.push(image);
        let mut index = 0;
        while index < 3 {
            if sizes[index] == 0 {
                index += 1;
                continue;
            }
            residues[index] += 1;
            if residues[index] < sizes[index] {
                break;
            }
            residues[index] = 0;
            index += 1;
        }
    }
    let mut kernel = Vec::new();
    for index in 0..3 {
        if sizes[index] != 0 {
            continue;
        }
        let mut vector = [0i128; 3];
        for row in 0..3 {
            vector[row] = basis[row][index];
        }
        kernel.push(vector);
    }
    Ok(Some(JointSolution { branches, kernel }))
}

/// Step counter for the two documented caps.
#[derive(Default)]
struct Budget {
    steps: u64,
}

impl Budget {
    fn step(&mut self) -> Result<(), ShiftEnumerationError> {
        self.steps += 1;
        if self.steps > MAX_STEPS {
            return Err(ShiftEnumerationError::Capacity);
        }
        Ok(())
    }
}

/// Whether `value` is an integer multiple of `modulus`.
fn is_multiple_of_integer(value: &Rat, modulus: i128) -> bool {
    if value.is_zero() {
        return true;
    }
    match value.checked_div(Rat::from_integer(modulus)) {
        Ok(quotient) => quotient.is_integer(),
        Err(_) => false,
    }
}

/// `value - factor * other`, checked.
fn subtract_product(value: i128, factor: i128, other: i128) -> Result<i128, ShiftEnumerationError> {
    factor
        .checked_mul(other)
        .and_then(|product| value.checked_sub(product))
        .ok_or(ShiftEnumerationError::Arithmetic)
}

/// `value` reduced into `[0, 1)` as an exact rational.
fn fraction_part(value: Rat) -> Result<Rat, ShiftEnumerationError> {
    Rat::new(
        value.numerator().rem_euclid(value.denominator()),
        value.denominator(),
    )
    .map_err(|_| ShiftEnumerationError::Arithmetic)
}

/// Integer row echelon of `[rows | rhs]` over the three columns; returns
/// `rank rows`.
///
/// `rows[..rank]` is an upper-triangular integer basis of the row module with
/// the transformed right-hand sides; `rows[rank..]` are zero rows.  Pivot
/// selection takes the smallest non-zero magnitude of the column, which is what
/// makes each pass strictly decrease it and the loop terminate.
fn row_echelon(
    rows: &mut [([i128; 3], Rat)],
    budget: &mut Budget,
) -> Result<usize, ShiftEnumerationError> {
    let mut pivot = 0usize;
    for column in 0..3 {
        if pivot == rows.len() {
            break;
        }
        loop {
            budget.step()?;
            // `unsigned_abs` keeps `i128::MIN` usable as a magnitude (its
            // `abs` would panic in debug and wrap in release).
            let mut best: Option<(usize, u128)> = None;
            for (row, entry) in rows.iter().enumerate().skip(pivot) {
                let value = entry.0[column];
                if value == 0 {
                    continue;
                }
                let magnitude = value.unsigned_abs();
                if best.is_none_or(|(_, smallest)| magnitude < smallest) {
                    best = Some((row, magnitude));
                }
            }
            let Some((found, _)) = best else { break };
            rows.swap(pivot, found);
            let pivot_row = rows[pivot];
            let mut remainder = false;
            for entry in rows.iter_mut().skip(pivot + 1) {
                if entry.0[column] == 0 {
                    continue;
                }
                budget.step()?;
                let factor = entry.0[column]
                    .checked_div(pivot_row.0[column])
                    .ok_or(ShiftEnumerationError::Arithmetic)?;
                for (axis, value) in entry.0.iter_mut().enumerate() {
                    *value = subtract_product(*value, factor, pivot_row.0[axis])?;
                }
                entry.1 = entry
                    .1
                    .checked_sub(
                        pivot_row
                            .1
                            .checked_mul(Rat::from_integer(factor))
                            .map_err(|_| ShiftEnumerationError::Arithmetic)?,
                    )
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
                if entry.0[column] != 0 {
                    remainder = true;
                }
            }
            if !remainder {
                break;
            }
        }
        if rows[pivot].0[column] != 0 {
            pivot += 1;
        }
    }
    Ok(pivot)
}

/// Diagonalize the padded `3 x 3` matrix by unimodular row and column
/// operations, applying the row operations to `target` and accumulating the
/// column change in `basis` (so `matrix = U P basis` at every step, starting
/// from `basis = I`).
///
/// Each pass moves the smallest non-zero magnitude of the current submatrix to
/// the diagonal and clears its row and column with exact division remainders; a
/// non-zero remainder is strictly smaller, so the loop terminates.  A zero
/// remainder everywhere means the diagonal position is clean.
fn diagonalize(
    matrix: &mut [[i128; 3]; 3],
    target: &mut [Rat; 3],
    basis: &mut [[i128; 3]; 3],
    budget: &mut Budget,
) -> Result<(), ShiftEnumerationError> {
    for pivot in 0..3 {
        loop {
            budget.step()?;
            let mut best: Option<(usize, usize, u128)> = None;
            for (row, values) in matrix.iter().enumerate().skip(pivot) {
                for (column, value) in values.iter().enumerate().skip(pivot) {
                    if *value == 0 {
                        continue;
                    }
                    let magnitude = value.unsigned_abs();
                    if best.is_none_or(|(_, _, smallest)| magnitude < smallest) {
                        best = Some((row, column, magnitude));
                    }
                }
            }
            let Some((row, column, _)) = best else { break };
            if row != pivot {
                matrix.swap(pivot, row);
                target.swap(pivot, row);
            }
            if column != pivot {
                for line in 0..3 {
                    matrix[line].swap(pivot, column);
                    basis[line].swap(pivot, column);
                }
            }
            let pivot_value = matrix[pivot][pivot];
            let pivot_row = matrix[pivot];
            let pivot_target = target[pivot];
            let mut dirty = false;
            for other in (pivot + 1)..3 {
                if matrix[other][pivot] == 0 {
                    continue;
                }
                let factor = matrix[other][pivot]
                    .checked_div(pivot_value)
                    .ok_or(ShiftEnumerationError::Arithmetic)?;
                for axis in 0..3 {
                    matrix[other][axis] =
                        subtract_product(matrix[other][axis], factor, pivot_row[axis])?;
                }
                target[other] = target[other]
                    .checked_sub(
                        pivot_target
                            .checked_mul(Rat::from_integer(factor))
                            .map_err(|_| ShiftEnumerationError::Arithmetic)?,
                    )
                    .map_err(|_| ShiftEnumerationError::Arithmetic)?;
                if matrix[other][pivot] != 0 {
                    dirty = true;
                }
            }
            for other in (pivot + 1)..3 {
                if matrix[pivot][other] == 0 {
                    continue;
                }
                let factor = matrix[pivot][other]
                    .checked_div(pivot_value)
                    .ok_or(ShiftEnumerationError::Arithmetic)?;
                for values in matrix.iter_mut() {
                    values[other] = subtract_product(values[other], factor, values[pivot])?;
                }
                for values in basis.iter_mut() {
                    values[other] = subtract_product(values[other], factor, values[pivot])?;
                }
                if matrix[pivot][other] != 0 {
                    dirty = true;
                }
            }
            if !dirty {
                break;
            }
        }
    }
    Ok(())
}

/// Widen an integer matrix to `i128`.
fn ints_to_wide(matrix: &Mat3I) -> [[i128; 3]; 3] {
    let mut out = [[0i128; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            out[row][column] = i128::from(matrix[row][column]);
        }
    }
    out
}

/// Transpose a `3 x 3` matrix.
fn transpose(matrix: &[[i128; 3]; 3]) -> [[i128; 3]; 3] {
    let mut out = [[0i128; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            out[row][column] = matrix[column][row];
        }
    }
    out
}

/// Exact `3 x 3` product with checked arithmetic.
fn multiply(
    left: &[[i128; 3]; 3],
    right: &[[i128; 3]; 3],
) -> Result<[[i128; 3]; 3], ShiftEnumerationError> {
    let mut out = [[0i128; 3]; 3];
    for (row, values) in out.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            let mut sum = 0i128;
            for index in 0..3 {
                sum = sum
                    .checked_add(
                        left[row][index]
                            .checked_mul(right[index][column])
                            .ok_or(ShiftEnumerationError::Arithmetic)?,
                    )
                    .ok_or(ShiftEnumerationError::Arithmetic)?;
            }
            *value = sum;
        }
    }
    Ok(out)
}

/// Exact matrix-vector product with a rational vector, checked.
fn mul_vector(
    matrix: &[[i128; 3]; 3],
    vector: &[Rat; 3],
) -> Result<[Rat; 3], ShiftEnumerationError> {
    let mut out = [Rat::ZERO; 3];
    for (row, value) in out.iter_mut().enumerate() {
        let mut sum = Rat::ZERO;
        for (column, entry) in vector.iter().enumerate() {
            if matrix[row][column] == 0 {
                continue;
            }
            let term = entry
                .checked_mul(Rat::from_integer(matrix[row][column]))
                .map_err(|_| ShiftEnumerationError::Arithmetic)?;
            sum = sum
                .checked_add(term)
                .map_err(|_| ShiftEnumerationError::Arithmetic)?;
        }
        *value = sum;
    }
    Ok(out)
}

/// Exact `3 x 3` determinant with checked arithmetic (six-term Leibniz
/// expansion, the same shape the schoolbook uses).
fn determinant(matrix: &[[i128; 3]; 3]) -> Result<i128, ShiftEnumerationError> {
    // (row, column) triples of the six permutations, with their signs.
    let terms: [[i128; 7]; 6] = [
        [0, 0, 1, 1, 2, 2, 1],
        [0, 1, 1, 2, 2, 0, 1],
        [0, 2, 1, 0, 2, 1, 1],
        [0, 2, 1, 1, 2, 0, -1],
        [0, 0, 1, 2, 2, 1, -1],
        [0, 1, 1, 0, 2, 2, -1],
    ];
    let mut value = 0i128;
    for term in terms {
        let addend = matrix[term[0] as usize][term[1] as usize]
            .checked_mul(matrix[term[2] as usize][term[3] as usize])
            .and_then(|product| product.checked_mul(matrix[term[4] as usize][term[5] as usize]))
            .and_then(|product| product.checked_mul(term[6]))
            .ok_or(ShiftEnumerationError::Arithmetic)?;
        value = value
            .checked_add(addend)
            .ok_or(ShiftEnumerationError::Arithmetic)?;
    }
    Ok(value)
}

/// The two row/column indices left after deleting one, in **increasing** order.
///
/// The order matters: taking the remaining indices cyclically instead flips the
/// sign of the minor whenever the deleted index is `1`, and that mistake is
/// invisible for `adj(I) = I` and for `M = c I` (where `adj(C) c I C = c det(C)
/// I` holds whatever the adjugate's signs are).  Round 24 hit it at SG 167
/// UNI 1333, where the wrong adjugate broke the divisibility `P = adj(C) M C in
/// e Z^3` -- the precondition check reported `NotWellDefined` instead of
/// silently deciding on a wrong system.
const REMAINING: [[usize; 2]; 3] = [[1, 2], [0, 2], [0, 1]];

/// Cofactor `(-1)^(row+column)` times the minor of `matrix` at
/// `(row, column)`.
fn cofactor(
    matrix: &[[i128; 3]; 3],
    row: usize,
    column: usize,
) -> Result<i128, ShiftEnumerationError> {
    let [first_row, second_row] = REMAINING[row];
    let [first_column, second_column] = REMAINING[column];
    let minor = matrix[first_row][first_column]
        .checked_mul(matrix[second_row][second_column])
        .and_then(|value| {
            value.checked_sub(
                matrix[first_row][second_column].checked_mul(matrix[second_row][first_column])?,
            )
        })
        .ok_or(ShiftEnumerationError::Arithmetic)?;
    if (row + column).is_multiple_of(2) {
        Ok(minor)
    } else {
        minor.checked_neg().ok_or(ShiftEnumerationError::Arithmetic)
    }
}

/// Adjugate `adj(A)` with `adj(A)[row][column] = cofactor(A, column, row)`.
fn adjugate(matrix: &[[i128; 3]; 3]) -> Result<[[i128; 3]; 3], ShiftEnumerationError> {
    let mut out = [[0i128; 3]; 3];
    for (row, values) in out.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = cofactor(matrix, column, row)?;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn ratio(numerator: i128, denominator: i128) -> Rat {
        Rat::new(numerator, denominator).expect("small exact ratio")
    }

    /// The definition of the system, evaluated directly and **without** the
    /// solver: `rows x - rhs in modulus Z^m`.
    fn defines(rows: &[[i128; 3]], rhs: &[Rat], modulus: i128, point: &[Rat; 3]) -> bool {
        rows.iter().zip(rhs.iter()).all(|(row, value)| {
            let mut sum = Rat::ZERO;
            for (axis, entry) in point.iter().enumerate() {
                let term = entry
                    .checked_mul(Rat::from_integer(row[axis]))
                    .expect("small product");
                sum = sum.checked_add(term).expect("small sum");
            }
            let difference = sum.checked_sub(*value).expect("small difference");
            is_multiple_of_integer(&difference, modulus)
        })
    }

    /// Membership in the parent lattice, through the production reduction (a
    /// different code path from the solver's).
    fn in_parent_frame(point: [Rat; 3], lattice: &[[Rat; 3]]) -> bool {
        canonical_translation_checked(point, lattice)
            .unwrap()
            .iter()
            .all(|value| value.is_zero())
    }

    /// One equation of the complete solver, checked against the definition.
    fn satisfies(equation: &(Mat3I, [Rat; 3]), point: &[Rat; 3], lattice: &[[Rat; 3]]) -> bool {
        let mut difference = [Rat::ZERO; 3];
        for (row, value) in difference.iter_mut().enumerate() {
            let mut sum = Rat::ZERO;
            for (column, entry) in point.iter().enumerate() {
                let term = entry
                    .checked_mul(Rat::from_integer(i128::from(equation.0[row][column])))
                    .expect("small product");
                sum = sum.checked_add(term).expect("small sum");
            }
            *value = sum.checked_sub(equation.1[row]).expect("small difference");
        }
        in_parent_frame(difference, lattice)
    }

    fn key(point: &[Rat; 3]) -> [(i128, i128); 3] {
        point.map(|value| (value.numerator(), value.denominator()))
    }

    fn grid_solutions(
        rows: &[[i128; 3]],
        rhs: &[Rat],
        modulus: i128,
        resolution: i128,
    ) -> Vec<[Rat; 3]> {
        let mut found: BTreeSet<[(i128, i128); 3]> = BTreeSet::new();
        for x in 0..resolution {
            for y in 0..resolution {
                for z in 0..resolution {
                    let point = [
                        Rat::new(x, resolution).expect("grid point"),
                        Rat::new(y, resolution).expect("grid point"),
                        Rat::new(z, resolution).expect("grid point"),
                    ];
                    if defines(rows, rhs, modulus, &point) {
                        found.insert(key(&point));
                    }
                }
            }
        }
        found
            .into_iter()
            .map(|values| values.map(|(numerator, denominator)| ratio(numerator, denominator)))
            .collect()
    }

    /// Round 24 F1's SG 56 UNI 119 class 1, isolated: the **exact** equation
    /// `2 delta_1 = -1/2` has the single solution `1/4`, but the congruence
    /// `2 delta_1 == -1/2 (mod 1)` also has `3/4`, which differs from `1/4` by a
    /// half cell.  That lattice-shifted branch is exactly what the withdrawn
    /// `{+-1, +-1/2}` walk could not reach from the anchor's particular solution;
    /// the complete set must contain it, and the unconstrained second coordinate
    /// must come back as the kernel rather than as sampled points.
    #[test]
    fn the_shift_solution_set_contains_the_lattice_shifted_branch() {
        let equation = (
            [[2, 0, 0], [0, 0, 0], [0, 0, 2]],
            [ratio(-1, 2), Rat::ZERO, ratio(1, 2)],
        );
        let solution = solve_shift_components(&[equation], &[]).expect("exact");
        assert_eq!(
            solution.kernel,
            vec![[Rat::ZERO, Rat::ONE, Rat::ZERO]],
            "the unconstrained coordinate is the continuous kernel, not a sample"
        );
        let expected: Vec<[Rat; 3]> = vec![
            [ratio(1, 4), Rat::ZERO, ratio(1, 4)],
            [ratio(1, 4), Rat::ZERO, ratio(3, 4)],
            [ratio(3, 4), Rat::ZERO, ratio(1, 4)],
            [ratio(3, 4), Rat::ZERO, ratio(3, 4)],
        ];
        assert_eq!(solution.representatives, expected);
        for representative in &solution.representatives {
            assert!(
                satisfies(&equation, representative, &[]),
                "the representative {representative:?} must satisfy the definition"
            );
        }
        // The same equation through the independent grid oracle.  Because the
        // second coordinate is the continuous kernel, the quarter grid samples
        // the coset: every branch must be a grid point, and every grid point
        // must share a branch's pinned coordinates -- the kernel can move only
        // the second one.
        let rows = equation.0.map(|row| row.map(i128::from));
        let rhs = [equation.1[0], equation.1[1], equation.1[2]];
        let grid = grid_solutions(&rows, &rhs, 1, 4);
        for branch in &solution.representatives {
            assert!(
                grid.contains(branch),
                "branch {branch:?} is not a grid point"
            );
        }
        for point in &grid {
            assert!(
                solution
                    .representatives
                    .iter()
                    .any(|branch| branch[0] == point[0] && branch[2] == point[2]),
                "grid point {point:?} is in no component of the solver"
            );
        }
        assert_eq!(
            grid.len(),
            4 * 4,
            "the quarter grid samples the kernel coset"
        );
    }

    /// The determinant of a 3x3 integer matrix, six-term expansion -- an
    /// independent implementation, used for the count oracle below.
    fn determinant_of(rows: &[[i128; 3]; 3]) -> i128 {
        let mut value = 0i128;
        for (a, b, c, sign) in [
            (0usize, 1usize, 2usize, 1i128),
            (1, 2, 0, 1),
            (2, 0, 1, 1),
            (2, 1, 0, -1),
            (0, 2, 1, -1),
            (1, 0, 2, -1),
        ] {
            value += sign * rows[0][a] * rows[1][b] * rows[2][c];
        }
        value
    }

    /// The rank of the rows, by minors -- an independent computation.
    fn rank_of(rows: &[[i128; 3]]) -> usize {
        let any3 = if rows.len() == 3 {
            determinant_of(&[rows[0], rows[1], rows[2]]) != 0
        } else {
            false
        };
        if any3 {
            return 3;
        }
        for left in 0..rows.len() {
            for right in (left + 1)..rows.len() {
                for a in 0..3 {
                    for b in (a + 1)..3 {
                        if rows[left][a] * rows[right][b] - rows[left][b] * rows[right][a] != 0 {
                            return 2;
                        }
                    }
                }
            }
        }
        if rows.iter().any(|row| row.iter().any(|entry| *entry != 0)) {
            return 1;
        }
        0
    }

    /// The adjugate and the determinant against their definitions: `A adj(A) =
    /// det(A) I` on a deterministic list that includes matrices where deleting
    /// index `1` matters (the round-24 sign bug), singular matrices, and the
    /// basis that exposed it.
    #[test]
    fn the_adjugate_and_determinant_agree_with_their_definitions() {
        let matrices: [[[i128; 3]; 3]; 8] = [
            [[1, 2, 0], [3, 4, 0], [0, 0, 5]],
            [[0, 1, 0], [-1, 1, 0], [0, 0, -1]],
            [[1, 0, 0], [2, 3, 0], [2, 0, 3]],
            [[2, 0, 0], [0, 3, 0], [0, 0, 5]],
            [[0, 1, 2], [3, 4, 5], [6, 7, 8]],
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            [[-2, 5, 1], [3, -1, 4], [0, 6, -3]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        ];
        for matrix in matrices {
            let adjugate = adjugate(&matrix).expect("exact");
            let determinant = determinant(&matrix).expect("exact");
            let product = multiply(&matrix, &adjugate).expect("exact");
            let expected = [
                [determinant, 0, 0],
                [0, determinant, 0],
                [0, 0, determinant],
            ];
            assert_eq!(product, expected, "matrix {matrix:?}");
            let other = multiply(&adjugate, &matrix).expect("exact");
            assert_eq!(other, expected, "adj(A) A for {matrix:?}");
            // The six-term Leibniz expansion, written out again.
            let mut leibniz = 0i128;
            for (a, b, c, sign) in [
                (0usize, 1usize, 2usize, 1i128),
                (1, 2, 0, 1),
                (2, 0, 1, 1),
                (2, 1, 0, -1),
                (0, 2, 1, -1),
                (1, 0, 2, -1),
            ] {
                leibniz += sign * matrix[0][a] * matrix[1][b] * matrix[2][c];
            }
            assert_eq!(determinant, leibniz, "matrix {matrix:?}");
        }
    }

    /// The solver against the bare definition and against independent count and
    /// dimension facts, on a deterministic linear-congruential sweep of small
    /// systems.
    ///
    /// The oracle has four parts, and which of them applies is decided by the
    /// system, never by the solver's answer (except for the resolution cap,
    /// which only *narrows* the grid comparison):
    ///
    /// * **soundness**: every branch must satisfy `rows x - rhs in modulus Z^m`;
    /// * **continuous kernel**: every kernel vector must annihilate every row,
    ///   and there must be exactly `3 - rank(rows)` of them (the rank comes from
    ///   the minors of the rows);
    /// * **discrete completeness**: for a system whose kernel is trivial or
    ///   axis-aligned, a grid at the branches' own resolution -- refined once to
    ///   show it is saturated, and capped at 12 -- must contain exactly the
    ///   solver's components (the free coordinates are collapsed);
    /// * **count**: the solution set of a 3x3 system with a non-zero
    ///   determinant is one coset of `rows^-1(modulus Z^3)`, whose index over
    ///   `Z^3` is `|det rows| / modulus^3` -- that is how many classes modulo
    ///   `Z^3` the answer must have (the well-definedness precondition makes the
    ///   ratio an integer; a value below 1 would mean no solution, which the
    ///   sweep records as `empty_count`).
    #[test]
    fn the_joint_congruence_agrees_with_an_independent_oracle() {
        let fractions = [
            Rat::ZERO,
            ratio(1, 2),
            ratio(-1, 2),
            ratio(1, 3),
            ratio(-1, 4),
            ratio(1, 6),
        ];
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = |bound: u64| -> u64 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % bound
        };
        let mut grid_compared = 0usize;
        let mut structure_only = 0usize;
        let mut counted = 0usize;
        let mut empty_counted = 0usize;
        let mut infeasible = 0usize;
        for round in 0..400u64 {
            let row_count = 1 + next(3) as usize;
            let modulus = 1 + next(3) as i128;
            // Half the sweep is deliberately sparse (`{-1, 0, 1}` up to the
            // modulus): its kernels are axis-aligned often enough for the grid
            // comparison to apply, while the dense half exercises the count
            // oracle and the non-axis-aligned structure.
            let sparse = round % 2 == 0;
            // The solver's precondition: the congruence descends to `R^3/Z^3`
            // exactly when the modulus divides every coefficient.
            let mut rows = Vec::with_capacity(row_count);
            for _ in 0..row_count {
                let mut row = [0i128; 3];
                for entry in row.iter_mut() {
                    *entry = if sparse {
                        modulus * (next(3) as i128 - 1)
                    } else {
                        modulus * (next(7) as i128 - 3)
                    };
                }
                rows.push(row);
            }
            let rhs: Vec<Rat> = (0..row_count)
                .map(|_| fractions[next(fractions.len() as u64) as usize])
                .collect();
            let solution = solve_joint_congruence(&rows, &rhs, modulus).expect("decided");
            let Some(solution) = solution else {
                infeasible += 1;
                // Proved infeasible: no grid point may satisfy the definition.
                assert!(
                    grid_solutions(&rows, &rhs, modulus, 12).is_empty(),
                    "rows {rows:?} rhs {rhs:?} modulus {modulus}: proved infeasible, but a \
                     point of the grid satisfies the definition"
                );
                continue;
            };
            for branch in &solution.branches {
                assert!(
                    defines(&rows, &rhs, modulus, branch),
                    "rows {rows:?} rhs {rhs:?} modulus {modulus}: the branch {branch:?} does \
                     not satisfy the definition"
                );
            }
            for direction in &solution.kernel {
                for row in &rows {
                    let value: i128 = (0..3).map(|axis| row[axis] * direction[axis]).sum();
                    assert_eq!(
                        value, 0,
                        "rows {rows:?}: the kernel vector {direction:?} does not annihilate {row:?}"
                    );
                }
            }
            assert_eq!(
                solution.kernel.len(),
                3 - rank_of(&rows),
                "rows {rows:?}: the kernel must have dimension 3 - rank"
            );
            if rows.len() == 3 && rank_of(&rows) == 3 {
                // Count oracle: the solution set is one coset of
                // `Lambda = rows^-1(modulus Z^3)`, whose index over `Z^3` is
                // `|det rows| / modulus^3` (the well-definedness precondition
                // makes that an integer).
                let determinant = determinant_of(&[rows[0], rows[1], rows[2]]).abs();
                let total = modulus * modulus * modulus;
                if determinant % total == 0 {
                    assert_eq!(
                        solution.branches.len() as i128,
                        determinant / total,
                        "rows {rows:?} rhs {rhs:?} modulus {modulus}: the count oracle"
                    );
                    counted += 1;
                } else {
                    assert!(
                        solution.branches.is_empty(),
                        "rows {rows:?} rhs {rhs:?} modulus {modulus}: the count oracle says \
                         no solution"
                    );
                    empty_counted += 1;
                }
                continue;
            }
            // Discrete completeness on a saturated grid, restricted to kernels
            // that are trivial or axis-aligned (then the free coordinates are a
            // subset of the axes and collapsing is exact).
            let mut resolution = modulus;
            for branch in &solution.branches {
                for value in branch {
                    resolution = lcm(resolution, value.denominator());
                }
            }
            let axis_aligned = !solution.kernel.is_empty()
                && solution
                    .kernel
                    .iter()
                    .all(|direction| direction.iter().filter(|entry| **entry != 0).count() == 1)
                && solution
                    .kernel
                    .iter()
                    .map(|direction| (0..3).find(|axis| direction[*axis] != 0).expect("non-zero"))
                    .collect::<BTreeSet<_>>()
                    .len()
                    == solution.kernel.len();
            let trivial = solution.kernel.is_empty();
            if (trivial || axis_aligned) && resolution <= 24 {
                let free: BTreeSet<usize> = solution
                    .kernel
                    .iter()
                    .map(|direction| (0..3).find(|axis| direction[*axis] != 0).expect("non-zero"))
                    .collect();
                let shrunk = |point: &[Rat; 3]| -> Vec<(i128, i128)> {
                    (0..3)
                        .filter(|axis| !free.contains(axis))
                        .map(|axis| (point[axis].numerator(), point[axis].denominator()))
                        .collect()
                };
                let coarse = grid_solutions(&rows, &rhs, modulus, resolution);
                let fine = grid_solutions(&rows, &rhs, modulus, resolution * 2);
                let coarse_classes: BTreeSet<Vec<(i128, i128)>> =
                    coarse.iter().map(shrunk).collect();
                let fine_classes: BTreeSet<Vec<(i128, i128)>> = fine.iter().map(shrunk).collect();
                // With the free coordinates collapsed, refining the grid must add
                // nothing: the pinned coordinates are already exact.
                assert_eq!(
                    coarse_classes, fine_classes,
                    "rows {rows:?} rhs {rhs:?} modulus {modulus}: the collapsed grid is not \
                     saturated at resolution {resolution}"
                );
                let solver_classes: BTreeSet<Vec<(i128, i128)>> =
                    solution.branches.iter().map(shrunk).collect();
                assert_eq!(
                    coarse_classes, solver_classes,
                    "rows {rows:?} rhs {rhs:?} modulus {modulus}: solver components differ \
                     from the definition on the grid"
                );
                for branch in &solution.branches {
                    assert!(
                        coarse.contains(branch),
                        "rows {rows:?}: the branch {branch:?} is not a grid point"
                    );
                }
                grid_compared += 1;
            } else {
                structure_only += 1;
            }
        }
        println!(
            "sweep grid={grid_compared} structure={structure_only} count={counted} \
             empty_count={empty_counted} infeasible={infeasible}"
        );
        assert!(grid_compared > 20, "grid-compared systems: {grid_compared}");
        assert!(
            counted + empty_counted > 50,
            "count-oracle systems: {counted}"
        );
        assert!(
            structure_only > 20,
            "structure-only systems: {structure_only}"
        );
        assert!(infeasible > 0, "infeasible systems: {infeasible}");
    }

    /// A rank-deficient system: the pinned constants, the continuous kernel and
    /// the branch/coset structure all have to come out of the solver.
    #[test]
    fn the_joint_congruence_reports_the_continuous_kernel() {
        let rows = [[1i128, 0, 0]];
        let rhs = [ratio(1, 3)];
        let solution = solve_joint_congruence(&rows, &rhs, 1)
            .expect("decided")
            .expect("feasible");
        assert_eq!(solution.branches, vec![[ratio(1, 3), Rat::ZERO, Rat::ZERO]]);
        assert_eq!(
            solution.kernel,
            vec![[0, 1, 0], [0, 0, 1]],
            "x_2 and x_3 are unconstrained"
        );
        // Every point of the coset satisfies the definition, and a point off it
        // does not.
        for first in 0..4 {
            for second in 0..3 {
                let point = [ratio(1 + 3 * first, 3), ratio(second, 2), ratio(-1, 5)];
                assert!(defines(&rows, &rhs, 1, &point), "{point:?}");
            }
        }
        assert!(!defines(
            &rows,
            &rhs,
            1,
            &[ratio(2, 3), Rat::ZERO, Rat::ZERO]
        ));
        // The kernel is exactly `ker_R rows`: the solver's basis spans it, so a
        // zero coefficient matrix must return all of R^3 as kernel and a single
        // component.
        let zero = solve_joint_congruence(&[[0, 0, 0]], &[Rat::ZERO], 1)
            .expect("decided")
            .expect("feasible");
        assert_eq!(zero.branches, vec![[Rat::ZERO; 3]]);
        assert_eq!(zero.kernel.len(), 3);
    }

    /// `i128::MIN` has no positive counterpart, so every magnitude in the
    /// solver is taken with `unsigned_abs` (pivot selection) or `checked_abs`
    /// (the values that enter the arithmetic).  The boundary must be an explicit
    /// error, never a panic, a wrap-around or a verdict.
    #[test]
    fn the_minimum_boundary_is_an_error_not_an_answer() {
        // A modulus whose absolute value does not fit.
        assert_eq!(
            solve_joint_congruence(&[[1, 0, 0]], &[Rat::ZERO], i128::MIN),
            Err(ShiftEnumerationError::Arithmetic)
        );
        // A diagonal entry whose absolute value does not fit.
        assert_eq!(
            solve_joint_congruence(&[[i128::MIN, 0, 0]], &[Rat::ZERO], 1),
            Err(ShiftEnumerationError::Arithmetic)
        );
        // `i128::MIN` as a pivot magnitude must not panic: the entry is selected
        // and the answer is the same explicit error.
        assert!(
            solve_joint_congruence(
                &[[i128::MIN, i128::MIN, 0], [0, i128::MIN, 0]],
                &[Rat::ZERO, Rat::ZERO],
                1,
            )
            .is_err()
        );
        // A row/right-hand-side mismatch is an internal contract violation.
        assert_eq!(
            solve_joint_congruence(&[[1, 0, 0]], &[], 1),
            Err(ShiftEnumerationError::LengthMismatch { rows: 1, rhs: 0 })
        );
    }

    /// The boundaries: an inconsistency is an `Ok(None)` proof, while a
    /// component explosion and an `i128` overflow are **errors**, never a
    /// refutation (round 23 review F2's rule).
    #[test]
    fn the_joint_congruence_proves_inconsistency_but_never_capacity() {
        // 2 x_1 == 1/2 and x_1 == 0 have no common solution, although each row
        // alone is solvable.
        let inconsistent =
            solve_joint_congruence(&[[1, 0, 0], [2, 0, 0]], &[Rat::ZERO, ratio(1, 2)], 1)
                .expect("decided");
        assert_eq!(inconsistent, None, "the two rows contradict each other");
        // 100000 components exceed the documented cap.
        let capacity = solve_joint_congruence(&[[100_000, 0, 0]], &[Rat::ZERO], 1);
        assert_eq!(capacity, Err(ShiftEnumerationError::Capacity));
        // A product that leaves i128: the two large entries meet in the
        // elimination.
        let huge = 30_000_000_000_000_000_000i128;
        let overflow = solve_joint_congruence(&[[1, huge, 0], [huge, 0, 0]], &[Rat::ZERO; 2], 1);
        assert_eq!(overflow, Err(ShiftEnumerationError::Arithmetic));
    }

    /// A centred parent lattice (`L_P = Z^3 + (1/2,1/2,0)Z`): the reduction goes
    /// through `K = D L_P`, the modulus is the lattice index, and the answer is
    /// compared against a grid oracle in the **parent frame** that uses the
    /// production `canonical_translation` membership -- a different code path
    /// from the solver's integer elimination.
    #[test]
    fn the_centred_parent_lattice_is_handled_exactly() {
        let lattice = vec![[ratio(1, 2), ratio(1, 2), Rat::ZERO]];
        // `M = 2 I` is `I - R` for `R = -I`, an inversion, which preserves any
        // lattice -- the precondition of the reduction.  (The first draft used
        // `diag(2, 1, 2)`, which is not `I - R` for an invertible `R`: it does
        // not preserve the centred lattice, and the solver correctly refused the
        // answer instead of returning a branch list.)
        let equations = vec![(
            [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
            [ratio(1, 2), Rat::ZERO, ratio(1, 2)],
        )];
        let solution = solve_shift_components(&equations, &lattice).expect("exact");
        assert!(
            solution.kernel.is_empty(),
            "the coefficient matrix is invertible, so the solution set is discrete"
        );
        assert!(!solution.representatives.is_empty());
        for representative in &solution.representatives {
            assert!(
                satisfies(&equations[0], representative, &lattice),
                "representative {representative:?} must land in the parent lattice"
            );
        }
        // Independent oracle: enumerate the parent-frame grid and compare.
        let resolution = 8i128;
        let mut grid: BTreeSet<[(i128, i128); 3]> = BTreeSet::new();
        for x in 0..resolution {
            for y in 0..resolution {
                for z in 0..resolution {
                    let point = [
                        ratio(x, resolution),
                        ratio(y, resolution),
                        ratio(z, resolution),
                    ];
                    if satisfies(&equations[0], &point, &lattice) {
                        let canonical = to_parent_frame(point, &lattice).unwrap();
                        grid.insert(key(&canonical));
                    }
                }
            }
        }
        let grid: Vec<[Rat; 3]> = grid
            .into_iter()
            .map(|values| values.map(|(numerator, denominator)| ratio(numerator, denominator)))
            .collect();
        assert_eq!(
            solution.representatives, grid,
            "the centred-lattice solution set must agree with the grid oracle"
        );
        // The finer grid must not add anything: the oracle is saturated.
        let refined = 2 * resolution;
        let mut finer: BTreeSet<[(i128, i128); 3]> = BTreeSet::new();
        for x in 0..refined {
            for y in 0..refined {
                for z in 0..refined {
                    let point = [ratio(x, refined), ratio(y, refined), ratio(z, refined)];
                    if satisfies(&equations[0], &point, &lattice) {
                        finer.insert(key(&to_parent_frame(point, &lattice).unwrap()));
                    }
                }
            }
        }
        assert_eq!(finer.len(), grid.len(), "the grid oracle is saturated");
    }
    fn lcm(left: i128, right: i128) -> i128 {
        left / gcd(left, right) * right
    }
}
