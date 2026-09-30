//! R9 stage 1: the magnetic operation set of one UNI, exactly, as a **group**.
//!
//! Milestone `docs/subduction-next-milestones.md` R9 ("磁输入契约与真实嵌入")
//! starts by fixing the input contract: an ordinary space-group irrep label does
//! not determine the antiunitary action, so a magnetic request has to name a
//! magnetic space group (UNI number) whose operations are then *verified*, never
//! guessed.  This module is the first layer of that:
//!
//! * [`magnetic_operations`] turns one UNI's operations into **exact** rational
//!   data ([`MagneticOperation`]) with the time-reversal flag preserved.  The
//!   translations of the bundled magnetic database live on denominators up to
//!   6 (measured over all 1,421 UNIs referenced by the isotropy tables); a
//!   translation that is not a small fraction is refused rather than rounded.
//! * The database lists operations **modulo the group's own translation
//!   lattice**: UNI 3 contains the primed operation `(E | 1/2)` whose square is
//!   `(E | 1)`, and `1` is a lattice vector of that group, so the square *is*
//!   the identity class.  [`translation_lattice`] extracts that lattice (the
//!   pure translations of the set) and [`canonical_translation`] reduces any
//!   translation modulo it, so every structural claim below is a claim about
//!   classes, which is what the data means.
//! * [`verify_group`] **proves** the set is a group: it closes a small
//!   generating set breadth-first (composition followed by
//!   [`canonical_translation`]) and requires the closure to cover every
//!   operation class, so identity, closure and inverses all follow without a
//!   quadratic scan.  It additionally proves that the unitary operations form a
//!   subgroup and that the antiunitary operations form a single coset of it.
//!
//! The database lists a **centred** group both ways: UNI 20 contains `(E | 0)`
//! and `(E | 1/2,1/2,0)`, which are the same element modulo that group's
//! centring lattice, so the list has more operations than the quotient
//! `G / L` has classes.  [`MagneticOperationSet`] therefore exposes both views
//! — `operations` (as listed, exact, in database order) and `classes` (one
//! representative per coset of the translation lattice) — and every structural
//! claim is a claim about `classes`.
//!
//! What this stage deliberately does **not** do (that is stage 2): express the
//! operations in the parent's frame through the isotropy record's `basis` and
//! `origin`.  Building the embedding before the frame convention is pinned
//! against an independent anchor would be exactly the guessing R9 forbids.

use std::collections::HashMap;

use crate::irrep::corep;
use crate::irrep::subduction::{Mat3I, Mat3R, Rat, Vec3R};

/// Largest denominator [`magnetic_operations`] accepts when it converts the
/// database's floating-point translations to exact fractions.
const MAX_TRANSLATION_DENOMINATOR: i128 = 96;

/// Lattice coefficients tried when a translation is reduced to its canonical
/// class representative.  The database's pure translations have components at
/// most 1, so `-2..=2` covers every reduction needed for a translation in
/// `[0, 1)`.
const REDUCTION_RANGE: i32 = 2;

/// One magnetic symmetry operation in the setting of its UNI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MagneticOperation {
    /// Integer rotation matrix.
    pub rotation: Mat3I,
    /// Fractional translation, exact.  Meaningful modulo the group's own
    /// translation lattice; see [`canonical_translation`].
    pub translation: [Rat; 3],
    /// `false` = unitary (unprimed), `true` = antiunitary (primed).
    pub time_reversal: bool,
}

/// The operation set of one magnetic space group (UNI number), exactly.
#[derive(Debug, Clone)]
pub struct MagneticOperationSet {
    /// Magnetic space group UNI number (1–1651).
    pub uni: usize,
    /// Operations in the database's order, exactly as listed.
    pub operations: Vec<MagneticOperation>,
    /// One representative per class modulo the translation lattice
    /// ([`translation_lattice`]), in first-seen order.
    pub classes: Vec<MagneticOperation>,
    /// How many listed operations are unitary (unprimed).
    pub unitary_count: usize,
    /// How many listed operations are antiunitary (primed).
    pub antiunitary_count: usize,
    /// How many classes are unitary.
    pub unitary_classes: usize,
    /// How many classes are antiunitary.
    pub antiunitary_classes: usize,
}

impl MagneticOperationSet {
    /// Total number of listed operations (may exceed the class count for a
    /// centred group; see the module docs).
    pub fn len(&self) -> usize {
        self.operations.len()
    }

    /// Number of classes modulo the translation lattice, i.e. the order of the
    /// quotient the group structure is about.
    pub fn class_count(&self) -> usize {
        self.classes.len()
    }

    /// Whether the set is empty (never for a known UNI).
    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }

    /// Index of the unitary subgroup among the **classes**, which is `1` or
    /// `2` for every magnetic space group.
    pub fn unitary_index(&self) -> usize {
        self.class_count()
            .checked_div(self.unitary_classes)
            .unwrap_or(0)
    }
}

/// Why a magnetic operation set could not be built or verified.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MagneticContractError {
    /// The UNI number is outside 1–1651 or has no operations.
    UnknownUni { uni: usize },
    /// A translation is not a fraction with a small denominator, so the exact
    /// contract cannot represent it (never rounded silently).
    TranslationOffGrid {
        uni: usize,
        index: usize,
        component: usize,
        value: String,
    },
    /// The operations do not generate a set equal to themselves: identity,
    /// closure or inverses fail.
    NotAGroup {
        uni: usize,
        generated: usize,
        total: usize,
    },
    /// The unitary operations alone are not closed.
    UnitaryNotClosed {
        uni: usize,
        generated: usize,
        unitary: usize,
    },
    /// The antiunitary operations are not a single coset of the unitary
    /// subgroup.
    AntiunitaryNotCoset {
        uni: usize,
        coset: usize,
        antiunitary: usize,
    },
    /// No operation at all (empty set).
    Empty { uni: usize },
    /// The set does not match a fresh read of the database it came from.
    DatabaseMismatch {
        uni: usize,
        index: usize,
        field: &'static str,
    },
    /// The set's own bookkeeping (classes and counts) contradicts its
    /// operations, so it was mutated after construction.
    MetadataMismatch {
        uni: usize,
        field: &'static str,
        recorded: usize,
        derived: usize,
    },
    /// Exact arithmetic failed (overflow or a non-integral rotation product).
    Arithmetic,
}

impl std::fmt::Display for MagneticContractError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownUni { uni } => write!(f, "no magnetic operations for UNI {uni}"),
            Self::TranslationOffGrid {
                uni,
                index,
                component,
                value,
            } => write!(
                f,
                "UNI {uni} operation {index} translation component {component} = {value} is not \
                 a fraction with denominator <= {MAX_TRANSLATION_DENOMINATOR}"
            ),
            Self::NotAGroup {
                uni,
                generated,
                total,
            } => write!(
                f,
                "UNI {uni}: the operations generate only {generated} of {total} classes, so the \
                 set is not closed"
            ),
            Self::UnitaryNotClosed {
                uni,
                generated,
                unitary,
            } => write!(
                f,
                "UNI {uni}: the unitary operations generate {generated} but the subgroup has \
                 {unitary} elements"
            ),
            Self::AntiunitaryNotCoset {
                uni,
                coset,
                antiunitary,
            } => write!(
                f,
                "UNI {uni}: one antiunitary coset has {coset} classes but {antiunitary} \
                 antiunitary operations exist"
            ),
            Self::Empty { uni } => write!(f, "UNI {uni} has no operations"),
            Self::DatabaseMismatch { uni, index, field } => write!(
                f,
                "UNI {uni}: operation {index} differs from the database in {field}"
            ),
            Self::MetadataMismatch {
                uni,
                field,
                recorded,
                derived,
            } => write!(
                f,
                "UNI {uni}: {field} is recorded as {recorded} but the operations give {derived}"
            ),
            Self::Arithmetic => write!(f, "exact arithmetic failed while composing operations"),
        }
    }
}

impl std::error::Error for MagneticContractError {}

/// Exact conversion of one database translation component.
fn exact_component(
    uni: usize,
    index: usize,
    component: usize,
    value: f64,
) -> Result<Rat, MagneticContractError> {
    for denominator in 1..=MAX_TRANSLATION_DENOMINATOR {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < 1e-9 {
            return Rat::new(scaled.round() as i128, denominator).map_err(|_| {
                MagneticContractError::TranslationOffGrid {
                    uni,
                    index,
                    component,
                    value: format!("{value}"),
                }
            });
        }
    }
    Err(MagneticContractError::TranslationOffGrid {
        uni,
        index,
        component,
        value: format!("{value}"),
    })
}

/// The exact operations of one magnetic space group (UNI number).
///
/// The operations keep the database's order and the time-reversal flag; the
/// translations are converted to [`Rat`] and refused when off the small grid.
pub fn magnetic_operations(uni: usize) -> Result<MagneticOperationSet, MagneticContractError> {
    let ops =
        corep::get_magnetic_operations(uni).ok_or(MagneticContractError::UnknownUni { uni })?;
    if ops.is_empty() {
        return Err(MagneticContractError::Empty { uni });
    }
    let mut operations = Vec::with_capacity(ops.len());
    let mut unitary_count = 0usize;
    for (index, op) in ops.iter().enumerate() {
        let translation = [
            exact_component(uni, index, 0, op.translation[0])?,
            exact_component(uni, index, 1, op.translation[1])?,
            exact_component(uni, index, 2, op.translation[2])?,
        ];
        if !op.time_reversal {
            unitary_count += 1;
        }
        operations.push(MagneticOperation {
            rotation: op.rotation,
            translation,
            time_reversal: op.time_reversal,
        });
    }
    let antiunitary_count = operations.len() - unitary_count;
    let mut set = MagneticOperationSet {
        uni,
        operations,
        classes: Vec::new(),
        unitary_count,
        antiunitary_count,
        unitary_classes: 0,
        antiunitary_classes: 0,
    };
    set.classes = classify(&set);
    set.unitary_classes = set
        .classes
        .iter()
        .filter(|operation| !operation.time_reversal)
        .count();
    set.antiunitary_classes = set.classes.len() - set.unitary_classes;
    Ok(set)
}

/// One representative per class modulo the set's translation lattice, in
/// first-seen order.
fn classify(set: &MagneticOperationSet) -> Vec<MagneticOperation> {
    let lattice = translation_lattice(set);
    let mut representatives: Vec<MagneticOperation> = Vec::new();
    let mut keys: Vec<MagneticKey> = Vec::new();
    for operation in &set.operations {
        let key = MagneticKey::of(operation, &lattice);
        if keys.contains(&key) {
            continue;
        }
        keys.push(key);
        representatives.push(*operation);
    }
    representatives
}

/// The translation lattice of a set: the non-zero pure translations
/// `(E | τ, unprimed)` it contains, exact.
///
/// These generate the lattice modulo which the database's operations are given;
/// a magnetic group's pure translations are always unitary (an antiunitary pure
/// translation would flip every spin without moving anything).
pub fn translation_lattice(set: &MagneticOperationSet) -> Vec<[Rat; 3]> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let mut lattice: Vec<[Rat; 3]> = Vec::new();
    for operation in &set.operations {
        if operation.rotation != identity || operation.time_reversal {
            continue;
        }
        if operation.translation.iter().all(|value| value.is_zero()) {
            continue;
        }
        if !lattice.contains(&operation.translation) {
            lattice.push(operation.translation);
        }
    }
    lattice
}

/// Component in `[0, 1)`: the fractional part of an exact rational.
fn fractional(value: Rat) -> Rat {
    let floor = value.numerator().div_euclid(value.denominator());
    Rat::new(
        value.numerator() - floor * value.denominator(),
        value.denominator(),
    )
    .expect("a non-zero denominator stays non-zero")
}

fn less(left: &[Rat; 3], right: &[Rat; 3]) -> bool {
    for axis in 0..3 {
        let a = left[axis];
        let b = right[axis];
        let lhs = a.numerator() * b.denominator();
        let rhs = b.numerator() * a.denominator();
        if lhs != rhs {
            return lhs < rhs;
        }
    }
    false
}

/// The canonical class representative of a translation modulo `lattice`, in
/// `[0, 1)` per component.
///
/// Every combination of the lattice generators with coefficients in
/// `-2..=2` is tried, reduced into `[0, 1)^3` and the exact lexicographic
/// minimum is returned, so the result depends only on the class.
pub fn canonical_translation(translation: [Rat; 3], lattice: &[[Rat; 3]]) -> [Rat; 3] {
    let mut best = [
        fractional(translation[0]),
        fractional(translation[1]),
        fractional(translation[2]),
    ];
    if lattice.is_empty() {
        return best;
    }
    let mut coefficients = vec![0i32; lattice.len()];
    loop {
        let mut candidate = translation;
        for (generator, coefficient) in lattice.iter().zip(coefficients.iter()) {
            for axis in 0..3 {
                let shift = generator[axis]
                    .checked_mul(Rat::from_integer(i128::from(*coefficient)))
                    .expect("small exact product");
                candidate[axis] = candidate[axis]
                    .checked_sub(shift)
                    .expect("small exact sum");
            }
        }
        let candidate = [
            fractional(candidate[0]),
            fractional(candidate[1]),
            fractional(candidate[2]),
        ];
        if less(&candidate, &best) {
            best = candidate;
        }
        // Odometer over -REDUCTION_RANGE..=REDUCTION_RANGE.
        let mut position = 0;
        loop {
            if position == coefficients.len() {
                return best;
            }
            coefficients[position] += 1;
            if coefficients[position] <= REDUCTION_RANGE {
                break;
            }
            coefficients[position] = -REDUCTION_RANGE;
            position += 1;
        }
    }
}

/// Group product `left * right` in the magnetic convention: rotations multiply,
/// translations add with the left rotation applied to the right translation,
/// and the time-reversal flags add modulo two.  The result's translation is a
/// raw sum; compare it through [`canonical_translation`].
pub fn compose(
    left: &MagneticOperation,
    right: &MagneticOperation,
) -> Result<MagneticOperation, MagneticContractError> {
    let left_rotation = Mat3R::from_ints(left.rotation);
    let right_rotation = Mat3R::from_ints(right.rotation);
    let rotation = left_rotation
        .checked_mul(&right_rotation)
        .and_then(|product| product.to_int_matrix())
        .map_err(|_| MagneticContractError::Arithmetic)?;
    let rotated = left_rotation
        .checked_mul_vector(&Vec3R::new(right.translation))
        .map_err(|_| MagneticContractError::Arithmetic)?;
    let translation = Vec3R::new(left.translation)
        .checked_add(&rotated)
        .map_err(|_| MagneticContractError::Arithmetic)?;
    Ok(MagneticOperation {
        rotation,
        translation: *translation.as_array(),
        time_reversal: left.time_reversal ^ right.time_reversal,
    })
}

/// Hashable exact key of an operation **class** modulo a fixed lattice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct MagneticKey {
    rotation: Mat3I,
    translation: [(i128, i128); 3],
    time_reversal: bool,
}

impl MagneticKey {
    fn of(operation: &MagneticOperation, lattice: &[[Rat; 3]]) -> Self {
        let translation = canonical_translation(operation.translation, lattice);
        Self {
            rotation: operation.rotation,
            translation: [
                (translation[0].numerator(), translation[0].denominator()),
                (translation[1].numerator(), translation[1].denominator()),
                (translation[2].numerator(), translation[2].denominator()),
            ],
            time_reversal: operation.time_reversal,
        }
    }
}

/// Breadth-first closure of `seed` inside `universe` (all keys modulo
/// `lattice`); returns the number of distinct classes reached, or `None` when a
/// product leaves `universe`.
fn closure_size(
    universe: &HashMap<MagneticKey, usize>,
    lattice: &[[Rat; 3]],
    seed: &[MagneticOperation],
) -> Option<usize> {
    let mut seen: HashMap<MagneticKey, ()> = HashMap::new();
    let mut frontier: Vec<MagneticOperation> = Vec::new();
    for operation in seed {
        let key = MagneticKey::of(operation, lattice);
        if !universe.contains_key(&key) {
            return None;
        }
        if seen.insert(key, ()).is_none() {
            frontier.push(*operation);
        }
    }
    while let Some(current) = frontier.pop() {
        for other in seed {
            for (left, right) in [(current, *other), (*other, current)] {
                let product = compose(&left, &right).ok()?;
                let key = MagneticKey::of(&product, lattice);
                if !universe.contains_key(&key) {
                    return None;
                }
                if seen.insert(key, ()).is_none() {
                    frontier.push(product);
                }
            }
        }
    }
    Some(seen.len())
}

/// Greedily generate `universe` from `candidates`, returning how many classes
/// were reached; closure under the chosen seeds is what "is a group" means here.
fn generate(
    universe: &HashMap<MagneticKey, usize>,
    lattice: &[[Rat; 3]],
    candidates: &[MagneticOperation],
) -> Result<usize, MagneticContractError> {
    let mut seeds: Vec<MagneticOperation> = Vec::new();
    let mut reached = 0usize;
    for operation in candidates {
        let key = MagneticKey::of(operation, lattice);
        if seeds
            .iter()
            .any(|seed| MagneticKey::of(seed, lattice) == key)
        {
            continue;
        }
        seeds.push(*operation);
        reached =
            closure_size(universe, lattice, &seeds).ok_or(MagneticContractError::NotAGroup {
                uni: 0,
                generated: seeds.len(),
                total: universe.len(),
            })?;
        if reached == universe.len() {
            break;
        }
    }
    Ok(reached)
}

/// Prove that `set` is a group (modulo its translation lattice), that its
/// unitary part is a subgroup and that its antiunitary part is one coset of that
/// subgroup.
pub fn verify_group(set: &MagneticOperationSet) -> Result<(), MagneticContractError> {
    // Re-derive everything from `operations` first: a set whose bookkeeping was
    // edited after construction must not be verified through the stale copy.
    let lattice = translation_lattice(set);
    let mut checks: [(&'static str, usize, usize); 5] = [("", 0, 0); 5];
    let derived_classes = classify(set);
    let derived_unitary = set
        .operations
        .iter()
        .filter(|operation| !operation.time_reversal)
        .count();
    let derived_unitary_classes = derived_classes
        .iter()
        .filter(|operation| !operation.time_reversal)
        .count();
    checks[0] = ("class count", set.class_count(), derived_classes.len());
    checks[1] = ("unitary operations", set.unitary_count, derived_unitary);
    checks[2] = (
        "antiunitary operations",
        set.antiunitary_count,
        set.operations.len() - derived_unitary,
    );
    checks[3] = (
        "unitary classes",
        set.unitary_classes,
        derived_unitary_classes,
    );
    checks[4] = (
        "antiunitary classes",
        set.antiunitary_classes,
        derived_classes.len() - derived_unitary_classes,
    );
    for (field, recorded, derived) in checks {
        if recorded != derived {
            return Err(MagneticContractError::MetadataMismatch {
                uni: set.uni,
                field,
                recorded,
                derived,
            });
        }
    }
    let classes = derived_classes;

    let universe: HashMap<MagneticKey, usize> = classes
        .iter()
        .enumerate()
        .map(|(index, operation)| (MagneticKey::of(operation, &lattice), index))
        .collect();
    if universe.len() != classes.len() {
        return Err(MagneticContractError::NotAGroup {
            uni: set.uni,
            generated: universe.len(),
            total: classes.len(),
        });
    }

    let reached = generate(&universe, &lattice, &classes)?;
    if reached != classes.len() {
        return Err(MagneticContractError::NotAGroup {
            uni: set.uni,
            generated: reached,
            total: classes.len(),
        });
    }

    // The unitary classes are a subgroup on their own.
    let unitary: Vec<MagneticOperation> = classes
        .iter()
        .copied()
        .filter(|operation| !operation.time_reversal)
        .collect();
    if unitary.len() != set.unitary_classes {
        return Err(MagneticContractError::UnitaryNotClosed {
            uni: set.uni,
            generated: unitary.len(),
            unitary: set.unitary_classes,
        });
    }
    let unitary_keys: HashMap<MagneticKey, usize> = unitary
        .iter()
        .enumerate()
        .map(|(index, operation)| (MagneticKey::of(operation, &lattice), index))
        .collect();
    if unitary_keys.len() != unitary.len() {
        return Err(MagneticContractError::UnitaryNotClosed {
            uni: set.uni,
            generated: unitary_keys.len(),
            unitary: unitary.len(),
        });
    }
    let reached = generate(&unitary_keys, &lattice, &unitary).map_err(|_| {
        MagneticContractError::UnitaryNotClosed {
            uni: set.uni,
            generated: 0,
            unitary: unitary.len(),
        }
    })?;
    if reached != unitary.len() {
        return Err(MagneticContractError::UnitaryNotClosed {
            uni: set.uni,
            generated: reached,
            unitary: unitary.len(),
        });
    }

    // The antiunitary classes are a single coset `a * H`.
    if set.antiunitary_classes > 0 {
        let representative = classes
            .iter()
            .find(|operation| operation.time_reversal)
            .expect("antiunitary classes exist");
        let mut coset: Vec<MagneticKey> = Vec::new();
        for unitary_operation in &unitary {
            let product = compose(representative, unitary_operation).map_err(|_| {
                MagneticContractError::AntiunitaryNotCoset {
                    uni: set.uni,
                    coset: 0,
                    antiunitary: set.antiunitary_classes,
                }
            })?;
            if !product.time_reversal {
                return Err(MagneticContractError::AntiunitaryNotCoset {
                    uni: set.uni,
                    coset: coset.len(),
                    antiunitary: set.antiunitary_classes,
                });
            }
            let key = MagneticKey::of(&product, &lattice);
            if !coset.contains(&key) {
                coset.push(key);
            }
        }
        let inside = coset.iter().all(|key| universe.contains_key(key));
        if !inside || coset.len() != set.antiunitary_classes {
            return Err(MagneticContractError::AntiunitaryNotCoset {
                uni: set.uni,
                coset: coset.len(),
                antiunitary: set.antiunitary_classes,
            });
        }
    }
    Ok(())
}

/// Compare a set with a **fresh read of the database** it came from.
///
/// This is the anchor the group axioms cannot provide.  A lattice rescaling
/// (editing a pure translation so that the set generates a different, equally
/// consistent magnetic group) is invisible to [`verify_group`]; it is caught
/// here, exactly as a flipped time-reversal flag, a moved coset representative
/// or a swapped pair of records is.  Stage 2's parent-frame embedding will use
/// the same anchor for `basis`/`origin`.
pub fn verify_against_database(set: &MagneticOperationSet) -> Result<(), MagneticContractError> {
    let fresh = magnetic_operations(set.uni)?;
    if fresh.operations.len() != set.operations.len() {
        return Err(MagneticContractError::DatabaseMismatch {
            uni: set.uni,
            index: set.operations.len().min(fresh.operations.len()),
            field: "operation count",
        });
    }
    for (index, (left, right)) in set.operations.iter().zip(fresh.operations.iter()).enumerate() {
        if left.rotation != right.rotation {
            return Err(MagneticContractError::DatabaseMismatch {
                uni: set.uni,
                index,
                field: "rotation",
            });
        }
        if left.time_reversal != right.time_reversal {
            return Err(MagneticContractError::DatabaseMismatch {
                uni: set.uni,
                index,
                field: "time reversal",
            });
        }
        if left.translation != right.translation {
            return Err(MagneticContractError::DatabaseMismatch {
                uni: set.uni,
                index,
                field: "translation",
            });
        }
    }
    Ok(())
}
