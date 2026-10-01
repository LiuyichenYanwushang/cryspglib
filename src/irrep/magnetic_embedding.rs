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
//! * [`canonical_translation`] reduces a translation modulo the lattice by
//!   walking the **whole coset** (breadth-first over the generators, in the
//!   finite torus of the common denominator), so the representative depends
//!   only on the class and no coefficient range has to be guessed.
//! * [`verify_group`] **proves** the set is a group with a direct check: every
//!   pair product must land back in the class set, the identity class must be
//!   present and every class must reach the identity by repeated multiplication
//!   (which also produces its inverse).  It additionally proves that the unitary
//!   operations form a subgroup and that the antiunitary operations form a
//!   single coset of it, and that the translation lattice is **rotation
//!   invariant** -- without that last check the quotient would not be a group
//!   even when the pairwise products happen to close (round-r9a finding: the
//!   generation argument silently assumed it).
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
//!
//! # What the group axioms cannot see (round r9a)
//!
//! A **self-consistent rescaling** of the ambient lattice -- editing a pure
//! translation so that every representative moves with it -- yields a different
//! but equally consistent magnetic group, and no axiom can reject it; a
//! *non*-consistent edit is caught ([`MagneticContractError::MetadataMismatch`]).
//! A group-wide flip of every time-reversal flag is likewise invisible for some
//! groups (the reviewer measured 517 of the 1,421 referenced UNIs passing it
//! once the bookkeeping is refreshed; the only literal `{E, T·E}` set in the
//! database is UNI 2, which the isotropy tables do not reference).  Both are why
//! [`verify_against_database`] exists and why R9's acceptance is "align the
//! operation set item by item with an independent source".

use std::collections::HashMap;

use crate::irrep::corep;
use crate::irrep::subduction::{Mat3I, Mat3R, Rat, Vec3R};

/// Largest denominator [`magnetic_operations`] accepts when it converts the
/// database's floating-point translations to exact fractions.
const MAX_TRANSLATION_DENOMINATOR: i128 = 96;

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
    /// The translation lattice is not invariant under the set's rotations, so
    /// the quotient cannot be a group (a conjugating rotation would leave the
    /// lattice, and with it the set).
    LatticeNotRotationInvariant { uni: usize },
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
            Self::LatticeNotRotationInvariant { uni } => write!(
                f,
                "UNI {uni}: the translation lattice is not invariant under the group's rotations"
            ),
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

/// Exact conversion of one floating-point fraction with a small denominator.
fn exact_fraction(value: f64) -> Option<Rat> {
    for denominator in 1..=MAX_TRANSLATION_DENOMINATOR {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < 1e-9 {
            return Rat::new(scaled.round() as i128, denominator).ok();
        }
    }
    None
}

/// Exact conversion of one database translation component.
fn exact_component(
    uni: usize,
    index: usize,
    component: usize,
    value: f64,
) -> Result<Rat, MagneticContractError> {
    exact_fraction(value).ok_or(MagneticContractError::TranslationOffGrid {
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
    lattice_of_operations(&set.operations)
}

/// The non-zero pure translations of an operation list, exact.
pub fn lattice_of_operations(operations: &[MagneticOperation]) -> Vec<[Rat; 3]> {
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let mut lattice: Vec<[Rat; 3]> = Vec::new();
    for operation in operations {
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

/// Greatest common divisor.
fn gcd(mut left: i128, mut right: i128) -> i128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.abs()
}

/// The canonical class representative of a translation modulo `lattice`, in
/// `[0, 1)` per component.
///
/// The coset `translation + lattice` is walked **exhaustively** in the finite
/// torus of the common denominator of all entries (breadth-first over the
/// lattice generators and their negatives), and the exact lexicographic minimum
/// is returned.  No coefficient range has to be chosen: round r9a showed that a
/// fixed range can miss the canonical representative (e.g. `1/2` modulo `1/6`
/// needs coefficient three), so the reduction is range-free by construction.
pub fn canonical_translation(translation: [Rat; 3], lattice: &[[Rat; 3]]) -> [Rat; 3] {
    let mut denominator = 1i128;
    for value in translation.iter().chain(lattice.iter().flatten()) {
        denominator = denominator / gcd(denominator, value.denominator()) * value.denominator();
    }
    let scale = |value: Rat| -> i128 {
        value.numerator() * (denominator / value.denominator())
    };
    let reduce = |value: i128| -> i128 { value.rem_euclid(denominator) };
    let start = [
        reduce(scale(translation[0])),
        reduce(scale(translation[1])),
        reduce(scale(translation[2])),
    ];
    let generators: Vec<[i128; 3]> = lattice
        .iter()
        .map(|generator| {
            [
                reduce(scale(generator[0])),
                reduce(scale(generator[1])),
                reduce(scale(generator[2])),
            ]
        })
        .filter(|generator| *generator != [0, 0, 0])
        .collect();

    let mut best = start;
    let mut seen: std::collections::HashSet<[i128; 3]> = std::collections::HashSet::new();
    let mut queue: Vec<[i128; 3]> = vec![start];
    seen.insert(start);
    while let Some(current) = queue.pop() {
        if current < best {
            best = current;
        }
        for generator in &generators {
            for sign in [1i128, -1] {
                let candidate = [
                    reduce(current[0] + sign * generator[0]),
                    reduce(current[1] + sign * generator[1]),
                    reduce(current[2] + sign * generator[2]),
                ];
                if seen.insert(candidate) {
                    queue.push(candidate);
                }
            }
        }
    }
    [
        Rat::new(best[0], denominator).expect("non-zero denominator"),
        Rat::new(best[1], denominator).expect("non-zero denominator"),
        Rat::new(best[2], denominator).expect("non-zero denominator"),
    ]
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

/// Whether `vector` lies in the lattice (mod the integer translations), judged
/// by whether its canonical class representative is the zero class.
fn in_lattice(vector: [Rat; 3], lattice: &[[Rat; 3]]) -> bool {
    canonical_translation(vector, lattice)
        .iter()
        .all(|value| value.is_zero())
}

/// Prove that `set` is a group (modulo its translation lattice), that its
/// unitary part is a subgroup and that its antiunitary part is one coset of that
/// subgroup.
///
/// The group claim is a direct check, not a generation shortcut: every pair of
/// classes must multiply into the class set, the identity class must be present,
/// and every class must return to the identity under repeated multiplication
/// (which exhibits its inverse).  The translation lattice must additionally be
/// invariant under every rotation in the set -- otherwise conjugating a
/// translation by that rotation leaves the lattice and the quotient is not a
/// group even if the products of the listed representatives happen to close.
pub fn verify_group(set: &MagneticOperationSet) -> Result<(), MagneticContractError> {
    // Re-derive everything from `operations` first: a set whose bookkeeping was
    // edited after construction must not be verified through the stale copy.
    let lattice = translation_lattice(set);
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
    let checks: [(&'static str, usize, usize); 5] = [
        ("class count", set.class_count(), derived_classes.len()),
        ("unitary operations", set.unitary_count, derived_unitary),
        (
            "antiunitary operations",
            set.antiunitary_count,
            set.operations.len() - derived_unitary,
        ),
        (
            "unitary classes",
            set.unitary_classes,
            derived_unitary_classes,
        ),
        (
            "antiunitary classes",
            set.antiunitary_classes,
            derived_classes.len() - derived_unitary_classes,
        ),
    ];
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

    // The lattice must be invariant under every rotation of the set.
    for operation in &classes {
        let rotation = Mat3R::from_ints(operation.rotation);
        for generator in &lattice {
            let image = rotation
                .checked_mul_vector(&Vec3R::new(*generator))
                .map_err(|_| MagneticContractError::Arithmetic)?;
            if !in_lattice(*image.as_array(), &lattice) {
                return Err(MagneticContractError::LatticeNotRotationInvariant { uni: set.uni });
            }
        }
    }

    // Full pairwise closure: every product of classes lands in the class set.
    for left in &classes {
        for right in &classes {
            let product = compose(left, right)?;
            if !universe.contains_key(&MagneticKey::of(&product, &lattice)) {
                return Err(MagneticContractError::NotAGroup {
                    uni: set.uni,
                    generated: 0,
                    total: classes.len(),
                });
            }
        }
    }

    // The identity class is present, and every class reaches it by repeated
    // multiplication (so every class has an inverse).
    let identity: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let is_identity = |operation: &MagneticOperation| {
        operation.rotation == identity
            && !operation.time_reversal
            && canonical_translation(operation.translation, &lattice)
                .iter()
                .all(|value| value.is_zero())
    };
    if !classes.iter().any(is_identity) {
        return Err(MagneticContractError::NotAGroup {
            uni: set.uni,
            generated: 0,
            total: classes.len(),
        });
    }
    for operation in &classes {
        let mut power = *operation;
        let mut reached = false;
        for _ in 0..classes.len() {
            if is_identity(&power) {
                reached = true;
                break;
            }
            power = compose(&power, operation)?;
        }
        if !reached {
            return Err(MagneticContractError::NotAGroup {
                uni: set.uni,
                generated: 0,
                total: classes.len(),
            });
        }
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
    for left in &unitary {
        for right in &unitary {
            let product = compose(left, right)?;
            if !unitary_keys.contains_key(&MagneticKey::of(&product, &lattice)) {
                return Err(MagneticContractError::UnitaryNotClosed {
                    uni: set.uni,
                    generated: 0,
                    unitary: unitary.len(),
                });
            }
        }
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

// ---------------------------------------------------------------------------
// R9 stage 2a: the record's geometry, exactly, and how far it embeds.
// ---------------------------------------------------------------------------

/// Exact geometry of one magnetic isotropy record.
///
/// The stored `basis` and `origin` are integers in the parent's **primitive**
/// frame (see [`crate::irrep::types::MagneticIsotropyRecord`]); this type keeps
/// them exactly, together with the lattice index `|det(basis)|`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagneticGeometry {
    /// Parent space-group number (1–230).
    pub parent_sg: u8,
    /// Magnetic space group UNI number.
    pub uni: usize,
    /// Subgroup lattice vectors in the parent's primitive frame (rows), exact.
    pub basis: [[Rat; 3]; 3],
    /// Origin shift in the parent's primitive frame, exact.
    pub origin: [Rat; 3],
    /// Absolute determinant of the basis, i.e. the lattice index.
    pub index: i128,
}

impl MagneticGeometry {
    /// The coordinate map from the subgroup's primitive frame into the parent's
    /// primitive frame: `M = basis^T` (the rows of `basis` are lattice vectors).
    pub fn map_matrix(&self) -> Mat3R {
        let mut out = [[Rat::ZERO; 3]; 3];
        for (row, values) in out.iter_mut().enumerate() {
            for (column, value) in values.iter_mut().enumerate() {
                *value = self.basis[column][row];
            }
        }
        Mat3R::new(out)
    }

    /// Whether another integer basis spans the lattice this geometry describes.
    ///
    /// Both directions of the change of basis must be integral (equivalently,
    /// the change is unimodular), which is the setting-free statement of "the
    /// two rows describe the same subgroup lattice".  This is the stronger
    /// cross-table comparison: the magnetic and ordinary tables often keep the
    /// same lattice in different settings, which an entry-by-entry comparison
    /// reports as a difference.
    pub fn spans_same_lattice_as(&self, other: [[i32; 3]; 3]) -> bool {
        let mut other_matrix = [[Rat::ZERO; 3]; 3];
        for (row, values) in other.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                other_matrix[row][column] = Rat::from_integer(i128::from(*value));
            }
        }
        let self_matrix = Mat3R::new(self.basis);
        let other_matrix = Mat3R::new(other_matrix);
        let Ok(self_inverse) = self_matrix.inverse() else {
            return false;
        };
        let Ok(other_inverse) = other_matrix.inverse() else {
            return false;
        };
        // `basis` rows are the lattice vectors, so the change of basis acts on
        // the RIGHT: `other = M * self` with `M = other * self^{-1}`.  Writing
        // this the other way round tests the column lattices and is wrong
        // whenever |det| > 1 (round r9b, item 1: 2,346 rows were rejected that
        // do span the same lattice).
        let forward = other_matrix.checked_mul(&self_inverse);
        let backward = self_matrix.checked_mul(&other_inverse);
        matches!(forward, Ok(matrix) if matrix.is_integral())
            && matches!(backward, Ok(matrix) if matrix.is_integral())
    }

    /// Whether the record keeps the parent's basis and origin.
    pub fn is_identity(&self) -> bool {
        self.basis == [[Rat::ONE, Rat::ZERO, Rat::ZERO], [Rat::ZERO, Rat::ONE, Rat::ZERO], [Rat::ZERO, Rat::ZERO, Rat::ONE]]
            && self.origin.iter().all(|value| value.is_zero())
    }
}

/// Why a record's geometry could not be read or does not embed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeometryError {
    /// The parent space-group number is outside 1–230.
    UnknownParentSpaceGroup { sg: u8 },
    /// The origin encoding `[x, y, z, d]` has a non-positive denominator.
    BadOriginEncoding { uni: usize, denominator: i32 },
    /// The basis is singular, so it does not describe a lattice.
    DegenerateBasis { uni: usize },
    /// The parent's operation table could not be read.
    UnreadableParent { sg: u8 },
    /// The parent's translations are not on the small exact grid.
    ParentTranslationOffGrid { sg: u8 },
}

impl std::fmt::Display for GeometryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownParentSpaceGroup { sg } => write!(f, "parent space group {sg} is unknown"),
            Self::BadOriginEncoding { uni, denominator } => write!(
                f,
                "UNI {uni}: origin denominator {denominator} is not positive"
            ),
            Self::DegenerateBasis { uni } => {
                write!(f, "UNI {uni}: the subgroup basis is singular")
            }
            Self::UnreadableParent { sg } => {
                write!(f, "the operations of parent space group {sg} are unavailable")
            }
            Self::ParentTranslationOffGrid { sg } => write!(
                f,
                "parent space group {sg} has a translation off the exact grid"
            ),
        }
    }
}

impl std::error::Error for GeometryError {}

/// Read one magnetic isotropy record's geometry exactly.
pub fn geometry_of(
    parent_sg: u8,
    record: &crate::irrep::types::MagneticIsotropyRecord,
) -> Result<MagneticGeometry, GeometryError> {
    if parent_sg == 0 || parent_sg > 230 {
        return Err(GeometryError::UnknownParentSpaceGroup { sg: parent_sg });
    }
    let denominator = i128::from(record.origin[3]);
    if denominator <= 0 {
        return Err(GeometryError::BadOriginEncoding {
            uni: record.mag_sg,
            denominator: record.origin[3],
        });
    }
    let mut basis = [[Rat::ZERO; 3]; 3];
    for (row, values) in record.basis.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            basis[row][column] = Rat::from_integer(i128::from(*value));
        }
    }
    let determinant = Mat3R::new(basis)
        .determinant()
        .map_err(|_| GeometryError::DegenerateBasis { uni: record.mag_sg })?;
    let index = determinant
        .to_integer()
        .map_err(|_| GeometryError::DegenerateBasis { uni: record.mag_sg })?
        .abs();
    if index == 0 {
        return Err(GeometryError::DegenerateBasis { uni: record.mag_sg });
    }
    let origin = [
        Rat::new(i128::from(record.origin[0]), denominator).expect("positive denominator"),
        Rat::new(i128::from(record.origin[1]), denominator).expect("positive denominator"),
        Rat::new(i128::from(record.origin[2]), denominator).expect("positive denominator"),
    ];
    Ok(MagneticGeometry {
        parent_sg,
        uni: record.mag_sg,
        basis,
        origin,
        index,
    })
}

/// How a record's geometry compares with the ordinary isotropy table's row for
/// the same parent irrep and direction label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometryAgreement {
    /// The ordinary row has the same basis and origin.
    Agree,
    /// The ordinary row differs in the basis only.
    BasisDiffers,
    /// The ordinary row differs in the origin only.
    OriginDiffers,
    /// The ordinary row differs in both.
    BothDiffer,
    /// The ordinary table has no row with this direction label.
    NoOrdinaryRow,
    /// The parent irrep or its isotropy table is unknown.
    NoOrdinaryIrrep,
}

/// Compare a magnetic record's geometry with the ordinary isotropy table.
///
/// Both tables are generated from the same program but by different code paths,
/// so agreement is a cross-source check of the stored geometry (and a tampered
/// basis or origin breaks it).  Disagreement is *not* an error: the two tables
/// legitimately pick different representatives for the same direction (measured:
/// 15,239 of 16,721 rows have a labelled ordinary row, 8,287 of those agree in
/// both basis and origin).
pub fn compare_with_ordinary_geometry(
    parent_sg: u8,
    ml: &str,
    direction_label: &str,
    geometry: &MagneticGeometry,
) -> GeometryAgreement {
    let Some(irrep) = crate::irrep::query::irreps_of(parent_sg)
        .iter()
        .find(|irrep| irrep.ml == ml)
    else {
        return GeometryAgreement::NoOrdinaryIrrep;
    };
    let mut saw_row = false;
    let mut basis_differs = false;
    let mut origin_differs = false;
    for ordinary in irrep.subgroups() {
        if ordinary.direction_label != direction_label {
            continue;
        }
        saw_row = true;
        let mut basis_same = true;
        for row in 0..3 {
            for column in 0..3 {
                if geometry.basis[row][column]
                    != Rat::from_integer(i128::from(ordinary.basis[row][column]))
                {
                    basis_same = false;
                }
            }
        }
        if ordinary.origin[3] <= 0 {
            origin_differs = true;
            continue;
        }
        let origin_same = (0..3).all(|axis| {
            geometry.origin[axis]
                == Rat::new(
                    i128::from(ordinary.origin[axis]),
                    i128::from(ordinary.origin[3]),
                )
                .expect("positive denominator")
        });
        if basis_same && origin_same {
            return GeometryAgreement::Agree;
        }
        if !basis_same {
            basis_differs = true;
        }
        if !origin_same {
            origin_differs = true;
        }
    }
    if !saw_row {
        return GeometryAgreement::NoOrdinaryRow;
    }
    match (basis_differs, origin_differs) {
        (true, true) => GeometryAgreement::BothDiffer,
        (true, false) => GeometryAgreement::BasisDiffers,
        (false, true) => GeometryAgreement::OriginDiffers,
        (false, false) => GeometryAgreement::Agree,
    }
}

/// How far one record's operation set embeds into its parent space group.
///
/// Every count is measured, never assumed.  `unitary_in_parent` counts the
/// unitary classes whose image (through [`MagneticGeometry::map_matrix`], in the
/// parent's primitive frame) is an operation of the parent space group modulo
/// the parent lattice; `antiunitary_rotations_in_parent` counts antiunitary
/// classes whose image rotation appears among the parent's rotations (so the
/// magnetic group sits inside the parent's grey group).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParentContainment {
    /// Classes of the magnetic operation set.
    pub classes: usize,
    /// Unitary classes of the magnetic operation set.
    pub unitary_classes: usize,
    /// Images whose rotation stayed integral.
    pub integral_rotations: usize,
    /// Unitary classes that are parent operations.
    pub unitary_in_parent: usize,
    /// Antiunitary classes whose rotation is a parent rotation.
    pub antiunitary_rotations_in_parent: usize,
}

/// Measure the parent-frame image of one record's operations.
///
/// The map used here is the record's own `basis`/`origin` in the parent's
/// **primitive** frame.  For a primitive parent that is the whole story; for a
/// centred parent the magnetic group's tabulated setting still has to be aligned
/// with the record's lattice, which is why this function reports counts instead
/// of asserting containment (see the R9 ledger entry).
pub fn measure_parent_containment(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Result<ParentContainment, GeometryError> {
    let parent = crate::irrep::query::symmetry_operations_of(geometry.parent_sg)
        .map_err(|_| GeometryError::UnreadableParent {
            sg: geometry.parent_sg,
        })?;
    let mut parent_operations = Vec::with_capacity(parent.operations.len());
    for operation in parent.operations.iter() {
        let translation = [
            exact_fraction(operation.translation[0]),
            exact_fraction(operation.translation[1]),
            exact_fraction(operation.translation[2]),
        ];
        let (Some(x), Some(y), Some(z)) = (translation[0], translation[1], translation[2]) else {
            return Err(GeometryError::ParentTranslationOffGrid {
                sg: geometry.parent_sg,
            });
        };
        parent_operations.push(MagneticOperation {
            rotation: operation.rotation,
            translation: [x, y, z],
            time_reversal: false,
        });
    }
    let parent_lattice = lattice_of_operations(&parent_operations);
    let parent_rotations: Vec<Mat3I> = parent_operations
        .iter()
        .map(|operation| operation.rotation)
        .collect();

    let map = geometry.map_matrix();
    let inverse = map.inverse().expect("non-singular basis");

    let mut result = ParentContainment {
        classes: set.class_count(),
        unitary_classes: set.unitary_classes,
        integral_rotations: 0,
        unitary_in_parent: 0,
        antiunitary_rotations_in_parent: 0,
    };
    for operation in &set.classes {
        let rotation = Mat3R::from_ints(operation.rotation);
        let conjugated = map
            .checked_mul(&rotation)
            .and_then(|product| product.checked_mul(&inverse))
            .expect("exact product");
        let Ok(mapped_rotation) = conjugated.to_int_matrix() else {
            continue;
        };
        result.integral_rotations += 1;
        // (I - R_p) * origin, with R_p integral.
        let mut shift = [Rat::ZERO; 3];
        for (row, cell) in shift.iter_mut().enumerate() {
            let mut sum = Rat::ZERO;
            for (column, value) in geometry.origin.iter().enumerate() {
                let entry = i128::from(if row == column { 1 } else { 0 })
                    - i128::from(mapped_rotation[row][column]);
                sum = sum
                    .checked_add(
                        value
                            .checked_mul(Rat::from_integer(entry))
                            .expect("exact product"),
                    )
                    .expect("exact sum");
            }
            *cell = sum;
        }
        let moved = map
            .checked_mul_vector(&Vec3R::new(operation.translation))
            .expect("exact product");
        let mut mapped_translation = [Rat::ZERO; 3];
        for axis in 0..3 {
            mapped_translation[axis] = moved
                .get(axis)
                .checked_add(shift[axis])
                .expect("exact sum");
        }
        if operation.time_reversal {
            if parent_rotations.contains(&mapped_rotation) {
                result.antiunitary_rotations_in_parent += 1;
            }
            continue;
        }
        let found = parent_operations.iter().any(|candidate| {
            if candidate.rotation != mapped_rotation {
                return false;
            }
            let mut difference = [Rat::ZERO; 3];
            for axis in 0..3 {
                difference[axis] = mapped_translation[axis]
                    .checked_sub(candidate.translation[axis])
                    .expect("exact difference");
            }
            canonical_translation(difference, &parent_lattice)
                .iter()
                .all(|value| value.is_zero())
        });
        if found {
            result.unitary_in_parent += 1;
        }
    }
    Ok(result)
}

/// The 48 signed permutation matrices, i.e. the finite setting ambiguity of a
/// cubic-axis convention.
fn signed_permutations() -> Vec<Mat3I> {
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut out = Vec::new();
    for permutation in permutations {
        for signs in 0..8 {
            let mut matrix = [[0i32; 3]; 3];
            for (row, column) in permutation.iter().enumerate() {
                let sign = if signs & (1 << row) == 0 { 1 } else { -1 };
                matrix[row][*column] = sign;
            }
            out.push(matrix);
        }
    }
    out
}

/// What the validated candidate-setting search found for one record.
///
/// This is an **instrument, not a decision procedure**: the candidate family is
/// the 48 signed permutations of the record's own basis, and no record has a
/// unique survivor, so a "the setting is pinned" claim cannot be made from it.
/// The survivor count is also **not** containment evidence: round 15k showed
/// that 436 records with zero survivors are already fully contained by the
/// naive map, and that UNI 1221/1333 contain under the identity map.  Read
/// [`measure_parent_containment`] for containment; read this only for how much
/// the setting family discriminates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SettingSearchReport {
    /// How many candidates were tried.
    pub candidates: usize,
    /// How many survived every check.
    pub survivors: usize,
    /// The naive map's containment measurement, for comparison.
    pub naive: ParentContainment,
}

/// Run the candidate-setting search for one record (see [`SettingSearchReport`]).
pub fn search_parent_setting(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Result<SettingSearchReport, GeometryError> {
    let parent = crate::irrep::query::symmetry_operations_of(geometry.parent_sg)
        .map_err(|_| GeometryError::UnreadableParent {
            sg: geometry.parent_sg,
        })?;
    let mut parent_operations = Vec::with_capacity(parent.operations.len());
    for operation in parent.operations.iter() {
        let translation = [
            exact_fraction(operation.translation[0]),
            exact_fraction(operation.translation[1]),
            exact_fraction(operation.translation[2]),
        ];
        let (Some(x), Some(y), Some(z)) = (translation[0], translation[1], translation[2]) else {
            return Err(GeometryError::ParentTranslationOffGrid {
                sg: geometry.parent_sg,
            });
        };
        parent_operations.push(MagneticOperation {
            rotation: operation.rotation,
            translation: [x, y, z],
            time_reversal: false,
        });
    }
    let parent_lattice = lattice_of_operations(&parent_operations);
    let parent_rotations: Vec<Mat3I> = parent_operations
        .iter()
        .map(|operation| operation.rotation)
        .collect();

    // NOTE (round 15k): the base map is the record's own basis, transposed.
    // Multiplying by the group's centring matrix here (as an earlier revision
    // did) is wrong: for UNI 1221 and UNI 1333 the *identity* map already maps
    // every listed operation onto a parent operation, yet that revision
    // reported zero survivors.  The centring is a property of the tabulated
    // cell, not an extra coordinate change to apply on top of the record.
    let convention = Mat3R::new(geometry.basis);
    let mut columns = [[Rat::ZERO; 3]; 3];
    for (row, values) in columns.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = convention.row(column)[row];
        }
    }
    let base = Mat3R::new(columns);

    let permutations = signed_permutations();
    let mut survivors = 0usize;
    for permutation in &permutations {
        let candidate = match Mat3R::from_ints(*permutation).checked_mul(&base) {
            Ok(matrix) => matrix,
            Err(_) => continue,
        };
        let Ok(inverse) = candidate.inverse() else {
            continue;
        };
        let mut accepted = true;
        for operation in &set.classes {
            let rotation = Mat3R::from_ints(operation.rotation);
            let Ok(conjugated) = candidate
                .checked_mul(&rotation)
                .and_then(|product| product.checked_mul(&inverse))
            else {
                accepted = false;
                break;
            };
            let Ok(mapped_rotation) = conjugated.to_int_matrix() else {
                accepted = false;
                break;
            };
            if operation.time_reversal {
                if !parent_rotations.contains(&mapped_rotation) {
                    accepted = false;
                    break;
                }
                continue;
            }
            let mut shift = [Rat::ZERO; 3];
            for (index, cell) in shift.iter_mut().enumerate() {
                let mut sum = Rat::ZERO;
                for (column, value) in geometry.origin.iter().enumerate() {
                    let entry = i128::from(if index == column { 1 } else { 0 })
                        - i128::from(mapped_rotation[index][column]);
                    sum = sum
                        .checked_add(
                            value
                                .checked_mul(Rat::from_integer(entry))
                                .expect("exact product"),
                        )
                        .expect("exact sum");
                }
                *cell = sum;
            }
            let moved = candidate
                .checked_mul_vector(&Vec3R::new(operation.translation))
                .map_err(|_| GeometryError::DegenerateBasis { uni: geometry.uni })?;
            let mut mapped_translation = [Rat::ZERO; 3];
            for axis in 0..3 {
                mapped_translation[axis] = moved
                    .get(axis)
                    .checked_add(shift[axis])
                    .expect("exact sum");
            }
            let found = parent_operations.iter().any(|candidate| {
                if candidate.rotation != mapped_rotation {
                    return false;
                }
                let mut difference = [Rat::ZERO; 3];
                for axis in 0..3 {
                    difference[axis] = mapped_translation[axis]
                        .checked_sub(candidate.translation[axis])
                        .expect("exact difference");
                }
                canonical_translation(difference, &parent_lattice)
                    .iter()
                    .all(|value| value.is_zero())
            });
            if !found {
                accepted = false;
                break;
            }
        }
        if accepted {
            survivors += 1;
        }
    }
    Ok(SettingSearchReport {
        candidates: permutations.len(),
        survivors,
        naive: measure_parent_containment(geometry, set)?,
    })
}

/// The image of the set's translation lattice under the record's map.
pub fn embedded_translation_lattice(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Vec<[Rat; 3]> {
    let map = geometry.map_matrix();
    let mut image: Vec<[Rat; 3]> = Vec::new();
    for generator in &translation_lattice(set) {
        let moved = map
            .checked_mul_vector(&Vec3R::new(*generator))
            .expect("exact product");
        let values = *moved.as_array();
        if !image.contains(&values) {
            image.push(values);
        }
    }
    image
}
