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

use std::collections::{HashMap, HashSet, VecDeque};

use crate::irrep::corep;
use crate::irrep::subduction::{Mat3I, Mat3R, Rat, SubductionError, Vec3R};

/// Round 24 (R9 stage 2b, F1): the exact **complete** shift solution set of one
/// fixed candidate map.  Kept in its own file because the elimination there is
/// independent of the embedding conventions here; see its module docs for the
/// contract and the completeness argument.
pub(crate) mod complete_shift;

/// The complete set of origin shifts that make one fixed candidate map an
/// embedding -- round 24's replacement for the finite `{+-1, +-1/2}` walk.
///
/// The solution set of a fixed map is a finite union of cosets of the common
/// real kernel
///
/// ```text
/// V = { v : (I - R_p) v = 0 for every class whose mapped rotation is R_p },
/// ```
///
/// so `representatives` lists **one** exact point per discrete component and
/// `kernel` describes the continuous directions instead of sampling them: every
/// shift that works is `representative + v + lambda` for some listed
/// representative, some `v` in `span_R(kernel)` and some `lambda` in the
/// parent's translation lattice.  Two shifts that differ by an element of `V`
/// produce the *same* image operation set, which is why the representatives --
/// not all rational points of the coset -- are the canonical finite answer.
///
/// An empty `representatives` list is a **proof** that no shift makes this map
/// an embedding: either a map-only gate fails (the image of the source
/// translation lattice leaves the parent's lattice, or a conjugated class
/// rotation is not a parent rotation) or the class congruences are jointly
/// inconsistent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShiftComponents {
    /// One representative per discrete component, canonical modulo the parent's
    /// translation lattice and sorted deterministically.  Empty means there is
    /// no solution.
    pub representatives: Vec<[Rat; 3]>,
    /// Basis of the common continuous kernel, in the parent's frame.  Empty when
    /// the shift coset is a single point modulo the lattice.
    pub kernel: Vec<[Rat; 3]>,
}

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
    /// Translation reduction exceeded its explicit coset budget; undecided.
    ReductionCapacity { limit: usize },
    /// The parent's primitive basis could not be read.
    ParentPrimitiveBasisUnavailable { sg: u8 },
    /// The parent's primitive basis is not on the exact small grid.
    ParentBasisOffGrid { sg: u8 },
    /// The record's basis is singular in the parent's conventional frame.
    SingularConventionalBasis { uni: usize },
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
            Self::Arithmetic => write!(
                f,
                "exact arithmetic failed while composing or reducing operations"
            ),
            Self::ReductionCapacity { limit } => write!(
                f,
                "magnetic translation reduction exceeded {limit} cosets; result unknown"
            ),
            Self::ParentPrimitiveBasisUnavailable { sg } => {
                write!(
                    f,
                    "the primitive basis of parent space group {sg} is unavailable"
                )
            }
            Self::ParentBasisOffGrid { sg } => write!(
                f,
                "the primitive basis of parent space group {sg} is off the exact grid"
            ),
            Self::SingularConventionalBasis { uni } => write!(
                f,
                "UNI {uni}: the record's basis is singular in the parent's conventional frame"
            ),
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
    set.classes = classify(&set)?;
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
fn classify(set: &MagneticOperationSet) -> Result<Vec<MagneticOperation>, MagneticContractError> {
    let lattice = translation_lattice(set);
    let mut representatives: Vec<MagneticOperation> = Vec::new();
    let mut keys: Vec<MagneticKey> = Vec::new();
    for operation in &set.operations {
        let key = MagneticKey::of(operation, &lattice)?;
        if keys.contains(&key) {
            continue;
        }
        keys.push(key);
        representatives.push(*operation);
    }
    Ok(representatives)
}

/// The non-zero pure translations `(E | τ, unprimed)` the set lists, exact.
///
/// These are the generators the database writes down; they are **not** the whole
/// translation lattice.  The missing part is the cell translations `Z^3` (the
/// operations are tabulated modulo the group's translation lattice, and
/// [`canonical_translation`] always reduces each component modulo one), which is
/// what [`full_translation_lattice`] adds.
///
/// Antiunitary pure translations `(E | τ, T)` need **no** extra generator, and
/// `full_translation_lattice` deliberately does not add one: such an operation
/// squares to the unitary lattice translation `(E | 2τ)`, which the quotient
/// already forces into `Z^3 +` this list.  Round 20 measured that directly (all
/// 1,651 UNIs, referenced and unreferenced): every product of listed operations
/// that is a pure translation, and every doubled antiunitary translation, lands
/// in `Z^3 +` this list -- 0 escapes.  The earlier doc claim that an antiunitary
/// pure translation "would flip every spin without moving anything" was wrong:
/// the corpus has them (788 with `τ ≠ 0` plus 230 with `τ = 0`).
pub fn translation_lattice(set: &MagneticOperationSet) -> Vec<[Rat; 3]> {
    lattice_of_operations(&set.operations)
}

/// Generators of the **full** translation lattice `L_M` of a set, exact.
///
/// `L_M = Z^3 + span_Z(listed non-zero unprimed pure translations)`, in the
/// frame the operations are tabulated in.  This is the lattice a containment
/// test has to push through the map: checking one representative per class says
/// nothing about it.  Round 20's counterexample is SG 143 `M1/S1` UNI 1, whose
/// rational dual map sends `e_2` to `(0, 1/2, 0)`: the class check passes (the
/// group has a single class) but the image of `Z^3` leaves the parent lattice.
pub fn full_translation_lattice(set: &MagneticOperationSet) -> Vec<[Rat; 3]> {
    let mut generators: Vec<[Rat; 3]> = Vec::new();
    for axis in 0..3 {
        let mut unit = [Rat::ZERO; 3];
        unit[axis] = Rat::ONE;
        generators.push(unit);
    }
    for generator in translation_lattice(set) {
        if !generators.contains(&generator) {
            generators.push(generator);
        }
    }
    generators
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
/// product of the three per-coordinate denominator grids (breadth-first over the
/// lattice generators and their negatives), and the exact lexicographic minimum
/// is returned.  No coefficient range has to be chosen: round r9a showed that a
/// fixed range can miss the canonical representative (e.g. `1/2` modulo `1/6`
/// needs coefficient three). The complete coset is explored within an explicit
/// state budget rather than an arbitrary coefficient box.
///
/// # Panics
/// Panics if a coordinate denominator exceeds i128 or the coset exceeds 65,536
/// states. Use [canonical_translation_checked] to receive these limits as errors.
pub fn canonical_translation(translation: [Rat; 3], lattice: &[[Rat; 3]]) -> [Rat; 3] {
    canonical_translation_checked(translation, lattice)
        .expect("exact magnetic translation reduction exceeded its arithmetic or capacity limit")
}

/// Checked lexicographic reduction modulo the listed lattice generators and
/// the implicit integer translations.
///
/// Each coordinate has its own common denominator. Unrelated denominators in
/// different coordinates must not be multiplied together. Integer parts are
/// removed *before* scaling, and modular addition never forms an overflowing
/// sum. Thus any returned representative is exact in both debug and release.
///
/// An arithmetic/capacity error is an **unknown** result, never evidence that a
/// vector is or is not in the lattice. The legacy [canonical_translation]
/// wrapper panics on these errors; validation paths use this checked entry.
///
/// ponytail: bounded coset enumeration (65,536 states); use a checked lattice
/// normal form if larger synthetic translation quotients need to be supported.
pub fn canonical_translation_checked(
    translation: [Rat; 3],
    lattice: &[[Rat; 3]],
) -> Result<[Rat; 3], MagneticContractError> {
    const MAX_STATES: usize = 65_536;
    let mut denominators = [1i128; 3];
    for axis in 0..3 {
        let mut denominator = translation[axis].denominator();
        for generator in lattice {
            let other = generator[axis].denominator();
            denominator = (denominator / gcd(denominator, other))
                .checked_mul(other)
                .ok_or(MagneticContractError::Arithmetic)?;
        }
        denominators[axis] = denominator;
    }
    let scale = |point: &[Rat; 3]| -> [i128; 3] {
        std::array::from_fn(|axis| {
            // 0 <= remainder < d and d divides D, hence remainder * (D/d) < D.
            let value = point[axis];
            value.numerator().rem_euclid(value.denominator())
                * (denominators[axis] / value.denominator())
        })
    };
    let start = scale(&translation);
    let generators: Vec<_> = lattice.iter().map(scale).filter(|g| *g != [0; 3]).collect();
    let mut seen = HashSet::from([start]);
    let mut queue = VecDeque::from([start]);
    let mut best = start;
    while let Some(point) = queue.pop_front() {
        best = best.min(point);
        for generator in &generators {
            for add in [true, false] {
                let next = std::array::from_fn(|axis| {
                    let (x, g, d) = (point[axis], generator[axis], denominators[axis]);
                    if add {
                        if x >= d - g { x - (d - g) } else { x + g }
                    } else if x < g {
                        d - (g - x)
                    } else {
                        x - g
                    }
                });
                if seen.insert(next) {
                    if seen.len() > MAX_STATES {
                        return Err(MagneticContractError::ReductionCapacity { limit: MAX_STATES });
                    }
                    queue.push_back(next);
                }
            }
        }
    }
    let mut reduced = [Rat::ZERO; 3];
    for axis in 0..3 {
        reduced[axis] = Rat::new(best[axis], denominators[axis])
            .map_err(|_| MagneticContractError::Arithmetic)?;
    }
    Ok(reduced)
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
    fn of(
        operation: &MagneticOperation,
        lattice: &[[Rat; 3]],
    ) -> Result<Self, MagneticContractError> {
        let translation = canonical_translation_checked(operation.translation, lattice)?;
        Ok(Self {
            rotation: operation.rotation,
            translation: [
                (translation[0].numerator(), translation[0].denominator()),
                (translation[1].numerator(), translation[1].denominator()),
                (translation[2].numerator(), translation[2].denominator()),
            ],
            time_reversal: operation.time_reversal,
        })
    }
}

/// Whether `vector` lies in the lattice (mod the integer translations), judged
/// by whether its canonical class representative is the zero class.
fn in_lattice(vector: [Rat; 3], lattice: &[[Rat; 3]]) -> Result<bool, MagneticContractError> {
    Ok(canonical_translation_checked(vector, lattice)?
        .iter()
        .all(|value| value.is_zero()))
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
    let derived_classes = classify(set)?;
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
        .map(|(index, operation)| Ok((MagneticKey::of(operation, &lattice)?, index)))
        .collect::<Result<_, MagneticContractError>>()?;
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
            if !in_lattice(*image.as_array(), &lattice)? {
                return Err(MagneticContractError::LatticeNotRotationInvariant { uni: set.uni });
            }
        }
    }

    // Full pairwise closure: every product of classes lands in the class set.
    for left in &classes {
        for right in &classes {
            let product = compose(left, right)?;
            if !universe.contains_key(&MagneticKey::of(&product, &lattice)?) {
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
    let identity_key = MagneticKey {
        rotation: identity,
        translation: [(0, 1); 3],
        time_reversal: false,
    };
    if !universe.contains_key(&identity_key) {
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
            if MagneticKey::of(&power, &lattice)? == identity_key {
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
        .map(|(index, operation)| Ok((MagneticKey::of(operation, &lattice)?, index)))
        .collect::<Result<_, MagneticContractError>>()?;
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
            if !unitary_keys.contains_key(&MagneticKey::of(&product, &lattice)?) {
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
            let product = compose(representative, unitary_operation)?;
            if !product.time_reversal {
                return Err(MagneticContractError::AntiunitaryNotCoset {
                    uni: set.uni,
                    coset: coset.len(),
                    antiunitary: set.antiunitary_classes,
                });
            }
            let key = MagneticKey::of(&product, &lattice)?;
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
    for (index, (left, right)) in set
        .operations
        .iter()
        .zip(fresh.operations.iter())
        .enumerate()
    {
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
        self.basis
            == [
                [Rat::ONE, Rat::ZERO, Rat::ZERO],
                [Rat::ZERO, Rat::ONE, Rat::ZERO],
                [Rat::ZERO, Rat::ZERO, Rat::ONE],
            ]
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
    /// Exact arithmetic failed, or bounded lattice reduction was undecided.
    /// Neither case is a proof of non-containment.
    Arithmetic,
    /// The record's basis is singular in the parent's conventional frame.
    SingularConventionalBasis { uni: usize },
    /// A claimed embedding does not land in the parent's operation set: the
    /// class with this rotation has no parent operation at its image.
    NotInParent { uni: usize, rotation: Mat3I },
    /// The parent's primitive basis could not be read.
    ParentPrimitiveBasisUnavailable { sg: u8 },
    /// The parent's primitive basis is not on the exact small grid.
    ParentBasisOffGrid { sg: u8 },
    /// The image of the group's full translation lattice leaves the parent's
    /// lattice.  The per-class checks cannot see this, because they only place
    /// one representative per class and [`canonical_translation`] reduces
    /// modulo the cell translations.
    SourceLatticeNotInParent { uni: usize },
    /// Round 24: the complete shift solution set could not be computed -- an
    /// exact-arithmetic limit or the documented component cap.  This is an
    /// *unknown* answer, never "no shift": a caller must not turn it into a
    /// refutation.
    ShiftEnumerationUnavailable { uni: usize },
    /// The supplied ordinary embedding cannot supply the unitary-subgroup
    /// bridge (mismatched identity, invalid magnetic data, or unavailable frame).
    /// This is not a proof that no magnetic embedding exists.
    UnitaryBridgeUnavailable { uni: usize, reason: String },
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
                write!(
                    f,
                    "the operations of parent space group {sg} are unavailable"
                )
            }
            Self::ParentTranslationOffGrid { sg } => write!(
                f,
                "parent space group {sg} has a translation off the exact grid"
            ),
            Self::Arithmetic => write!(
                f,
                "exact arithmetic or lattice reduction limit reached while converting frames"
            ),
            Self::SingularConventionalBasis { uni } => write!(
                f,
                "UNI {uni}: the record's basis is singular in the parent's conventional frame"
            ),
            Self::NotInParent { uni, rotation } => write!(
                f,
                "UNI {uni}: the claimed embedding does not map the class with rotation \
                 {rotation:?} into the parent"
            ),
            Self::ParentPrimitiveBasisUnavailable { sg } => {
                write!(
                    f,
                    "the primitive basis of parent space group {sg} is unavailable"
                )
            }
            Self::ParentBasisOffGrid { sg } => write!(
                f,
                "the primitive basis of parent space group {sg} is off the exact grid"
            ),
            Self::SourceLatticeNotInParent { uni } => write!(
                f,
                "UNI {uni}: the image of the group's full translation lattice leaves the \
                 parent's lattice"
            ),
            Self::UnitaryBridgeUnavailable { uni, reason } => {
                write!(
                    f,
                    "UNI {uni}: unitary-subgroup bridge unavailable: {reason}"
                )
            }
            Self::ShiftEnumerationUnavailable { uni } => write!(
                f,
                "UNI {uni}: the complete shift solution set could not be computed \
                 (exact-arithmetic limit or component cap); the answer is unknown"
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
/// Every count is measured, never assumed.  These are **diagnostic counts of the
/// naive raw map**, not a containment verdict: [`measure_parent_containment`]
/// applies `geometry.map_matrix()` -- the record's `basis` transposed, stored in
/// the parent's **primitive** frame -- to the parent's operations, which are
/// tabulated in the parent's **conventional** frame, without converting between
/// the two.  The counts reproduce that arithmetic exactly (an independent
/// checker rebuilt it on all 16,721 rows), which is what makes them a useful
/// instrument; it does not turn them into a physical embedding statement.
///
/// `unitary_in_parent` counts the unitary classes whose image is an operation of
/// the parent space group modulo the parent lattice;
/// `antiunitary_rotations_in_parent` counts antiunitary classes whose image
/// **rotation** appears among the parent's rotations -- that alone does **not**
/// place the class in the parent, so it is not evidence that the group sits
/// inside the parent's grey group; read `antiunitary_classes_in_parent` for that
/// (round 15p, L1).  Neither count checks the image of the full source
/// translation lattice (round 20); [`verify_embedding`] does.
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
    /// Antiunitary classes whose conjugated rotation is a parent rotation
    /// **only** -- the weaker reading, kept for continuity; it does not place
    /// the class in the parent (round 15o).
    pub antiunitary_rotations_in_parent: usize,
    /// Antiunitary classes fully placed in the parent: rotation is a parent
    /// rotation **and** the translation equation holds.  This is the class-level
    /// notion the older search uses; since round 20 [`verify_embedding`] adds the
    /// full-translation-lattice gate on top of it.
    pub antiunitary_classes_in_parent: usize,
}

/// Measure the naive map's class-level placement for one record.
///
/// The map is the record's own `basis` transposed with its `origin`, both stored
/// in the parent's **primitive** frame, while the parent operations come from
/// `symmetry_operations_of` in the parent's **conventional** frame.  This
/// function intentionally reproduces that raw arithmetic -- an independent
/// checker can and did rebuild it row by row -- so read the counts as a
/// diagnostic instrument, not as a frame-anchored embedding:
///
/// * it places one representative per class, so it says nothing about the image
///   of the group's full translation lattice (round 20; use
///   [`verify_embedding`] for that);
/// * it never converts primitive to conventional, so for a centred parent the
///   counts mix two frames.
///
/// [`basis_in_parent_conventional`] and [`origin_in_parent_conventional`] do the
/// conversion for callers that want the parent's conventional frame.
pub fn measure_parent_containment(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Result<ParentContainment, GeometryError> {
    let parent = crate::irrep::query::symmetry_operations_of(geometry.parent_sg).map_err(|_| {
        GeometryError::UnreadableParent {
            sg: geometry.parent_sg,
        }
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
    let inverse = map.inverse().map_err(|error| match error {
        SubductionError::SingularMatrix => GeometryError::DegenerateBasis { uni: geometry.uni },
        _ => GeometryError::Arithmetic,
    })?;

    let mut result = ParentContainment {
        classes: set.class_count(),
        unitary_classes: set.unitary_classes,
        integral_rotations: 0,
        unitary_in_parent: 0,
        antiunitary_rotations_in_parent: 0,
        antiunitary_classes_in_parent: 0,
    };
    for operation in &set.classes {
        let rotation = Mat3R::from_ints(operation.rotation);
        let conjugated = map
            .checked_mul(&rotation)
            .and_then(|product| product.checked_mul(&inverse))
            .map_err(|_| GeometryError::Arithmetic)?;
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
                            .map_err(|_| GeometryError::Arithmetic)?,
                    )
                    .map_err(|_| GeometryError::Arithmetic)?;
            }
            *cell = sum;
        }
        let moved = map
            .checked_mul_vector(&Vec3R::new(operation.translation))
            .map_err(|_| GeometryError::Arithmetic)?;
        let mut mapped_translation = [Rat::ZERO; 3];
        for axis in 0..3 {
            mapped_translation[axis] = moved
                .get(axis)
                .checked_add(shift[axis])
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        let found = if let Some(candidate) = parent_operations
            .iter()
            .find(|candidate| candidate.rotation == mapped_rotation)
        {
            let mut difference = [Rat::ZERO; 3];
            for axis in 0..3 {
                difference[axis] = mapped_translation[axis]
                    .checked_sub(candidate.translation[axis])
                    .map_err(|_| GeometryError::Arithmetic)?;
            }
            in_lattice(difference, &parent_lattice).map_err(|_| GeometryError::Arithmetic)?
        } else {
            false
        };
        if operation.time_reversal {
            if parent_rotations.contains(&mapped_rotation) {
                result.antiunitary_rotations_in_parent += 1;
            }
            if found {
                result.antiunitary_classes_in_parent += 1;
            }
            continue;
        }
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
/// Its acceptance criterion is the **class-level** one (rotation integral plus
/// the translation equation for one representative per class); it does not check
/// the image of the full source translation lattice and does not convert frames
/// (round 20), so a survivor is not a containment proof.
///
/// The survivor count is also **not** containment evidence: round 15k showed
/// that records with zero survivors can already be fully placed by the naive map
/// (413 of them under the strict all-classes diagnostic reading in force since
/// round 15p; 436 under that round's weaker, rotation-only reading for the
/// primed classes), and that UNI 1221/1333 place under the identity map.
/// For a verdict read [`verify_embedding`]; for the naive diagnostic read
/// [`measure_parent_containment`]; read this only for how much the setting
/// family discriminates.
///
/// Round 22 note: the *verdict* now comes from [`parent_embedding_verdict`],
/// which searches the wider `{S . B}` family and certifies an empty result
/// exactly.  This report stays as the narrow diagnostic and is kept for its
/// pinned histogram.
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
    let parent = crate::irrep::query::symmetry_operations_of(geometry.parent_sg).map_err(|_| {
        GeometryError::UnreadableParent {
            sg: geometry.parent_sg,
        }
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
            // Round 15o: the antiunitary classes get the same treatment as the
            // unitary ones -- time reversal does not change the translation
            // equation, so a rotation-only check would under-count.
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
                                .map_err(|_| GeometryError::Arithmetic)?,
                        )
                        .map_err(|_| GeometryError::Arithmetic)?;
                }
                *cell = sum;
            }
            let moved = candidate
                .checked_mul_vector(&Vec3R::new(operation.translation))
                .map_err(|_| GeometryError::Arithmetic)?;
            let mut mapped_translation = [Rat::ZERO; 3];
            for axis in 0..3 {
                mapped_translation[axis] = moved
                    .get(axis)
                    .checked_add(shift[axis])
                    .map_err(|_| GeometryError::Arithmetic)?;
            }
            let found = if let Some(candidate) = parent_operations
                .iter()
                .find(|candidate| candidate.rotation == mapped_rotation)
            {
                let mut difference = [Rat::ZERO; 3];
                for axis in 0..3 {
                    difference[axis] = mapped_translation[axis]
                        .checked_sub(candidate.translation[axis])
                        .map_err(|_| GeometryError::Arithmetic)?;
                }
                in_lattice(difference, &parent_lattice).map_err(|_| GeometryError::Arithmetic)?
            } else {
                false
            };
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

/// The record's basis expressed in the parent's **conventional** frame.
///
/// The stored basis is given in the parent's *primitive* frame (see
/// [`crate::irrep::types::MagneticIsotropyRecord`]), while the tabulated
/// operations live in the conventional frame.  Any comparison with those
/// operations must therefore go through this conversion, which is the same one
/// [`crate::irrep::isotropy::basis_in_parent_conventional`] performs -- exact
/// here, and refusing an off-grid parent basis instead of rounding.
pub fn basis_in_parent_conventional(
    geometry: &MagneticGeometry,
) -> Result<[[Rat; 3]; 3], GeometryError> {
    let primitive = parent_primitive_basis_exact(geometry.parent_sg)?;
    let mut converted = [[Rat::ZERO; 3]; 3];
    for (row, values) in converted.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            let mut sum = Rat::ZERO;
            for (inner, generator) in primitive.iter().enumerate() {
                sum = sum
                    .checked_add(
                        geometry.basis[row][inner]
                            .checked_mul(generator[column])
                            .map_err(|_| GeometryError::Arithmetic)?,
                    )
                    .map_err(|_| GeometryError::Arithmetic)?;
            }
            *value = sum;
        }
    }
    Ok(converted)
}

/// The record's origin expressed in the parent's conventional frame.
pub fn origin_in_parent_conventional(
    geometry: &MagneticGeometry,
) -> Result<[Rat; 3], GeometryError> {
    let primitive = parent_primitive_basis_exact(geometry.parent_sg)?;
    let mut converted = [Rat::ZERO; 3];
    for (column, value) in converted.iter_mut().enumerate() {
        let mut sum = Rat::ZERO;
        for (inner, generator) in primitive.iter().enumerate() {
            sum = sum
                .checked_add(
                    geometry.origin[inner]
                        .checked_mul(generator[column])
                        .map_err(|_| GeometryError::Arithmetic)?,
                )
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        *value = sum;
    }
    Ok(converted)
}

/// `|det|` of the record's basis in the parent's conventional frame.
///
/// This is the covolume ratio between the subgroup lattice and the parent's
/// conventional cell, so it is a **rational** number: for a centred parent the
/// conventional cell is larger than the primitive one the basis is stored in
/// (`|det| = |det(basis)| / index(parent centring)`).
pub fn conventional_lattice_index(geometry: &MagneticGeometry) -> Result<Rat, GeometryError> {
    let converted = basis_in_parent_conventional(geometry)?;
    let determinant = Mat3R::new(converted)
        .determinant()
        .map_err(|_| GeometryError::SingularConventionalBasis { uni: geometry.uni })?;
    Ok(if determinant.numerator() < 0 {
        determinant
            .checked_neg()
            .map_err(|_| GeometryError::Arithmetic)?
    } else {
        determinant
    })
}

/// The parent's primitive basis in its conventional frame, exactly.
fn parent_primitive_basis_exact(sg: u8) -> Result<[[Rat; 3]; 3], GeometryError> {
    let primitive = crate::irrep::isotropy::parent_primitive_basis(sg)
        .map_err(|_| GeometryError::ParentPrimitiveBasisUnavailable { sg })?;
    let mut exact = [[Rat::ZERO; 3]; 3];
    for (row, values) in primitive.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            exact[row][column] =
                exact_fraction(*value).ok_or(GeometryError::ParentBasisOffGrid { sg })?;
        }
    }
    Ok(exact)
}

/// Whether the image of the set's **full** translation lattice lies in the
/// parent's lattice -- the containment condition the per-class checks cannot
/// see, shared by [`verify_embedding`] and [`embed_in_parent_conventional`].
///
/// Round 20: without it a rational map whose image of `Z^3` leaves the parent
/// lattice was returned and accepted (SG 143 `M1/S1` UNI 1).  The check is a
/// function of the map, so the search runs it once per candidate, not once per
/// shift.
/// The round-20 containment gate in checked form: `Ok(false)` is a proof that
/// the map cannot embed (the image of the full source translation lattice leaves
/// the parent lattice), `Err` an exact-arithmetic limit -- an **unknown** answer,
/// not a refutation.
///
/// Round 24: both users of this gate (the search's `excluded_by_source_lattice`
/// screen and [`verify_embedding`]) take the checked form, so an arithmetic limit
/// can never be counted as a refutation nor reported as
/// [`GeometryError::SourceLatticeNotInParent`].
fn source_lattice_lands_in_parent_checked(
    map: &Mat3R,
    set: &MagneticOperationSet,
    parent_lattice: &[[Rat; 3]],
) -> Result<bool, GeometryError> {
    for generator in full_translation_lattice(set) {
        let image = map
            .checked_mul_vector(&Vec3R::new(generator))
            .map_err(|_| GeometryError::Arithmetic)?;
        if !in_lattice(*image.as_array(), parent_lattice).map_err(|_| GeometryError::Arithmetic)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether the map realises the record's lattice **exactly**: `A L_M = L_record`,
/// with an exact-arithmetic limit reported as an error instead of as `false`.
///
/// Both inclusions are needed.  Forward: every generator of the full source
/// lattice ([`full_translation_lattice`]) maps into the record's row lattice
/// (`(M^-1)^T v` integral, the row convention).  Reverse: every row of the
/// record's lattice pulls back through `A^-1` into `L_M` (membership modulo the
/// cell translations, via [`in_lattice`]).  Round 20 replaced the earlier test,
/// which looked only at the listed generators in the forward direction and was
/// therefore vacuously true whenever the set lists no pure translation.
///
/// This is the predicate behind [`ParentEmbedding::realises_record_lattice`],
/// exposed so that a caller who selected an embedding elsewhere can recompute the
/// flag instead of trusting it.  It derives the record's conventional-frame
/// lattice and the `(M^-1)^T` membership matrix once, in one place, because two
/// earlier versions of this test had the transpose or the frame wrong.
///
/// `Err` is an exact-arithmetic limit or a singular record basis; a singular
/// **map** is `Ok(false)`, because no shift makes a singular map a change of
/// basis.
pub fn realises_record_lattice_checked(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
    map: &Mat3R,
) -> Result<bool, GeometryError> {
    let record_lattice = basis_in_parent_conventional(geometry)?;
    let record_membership = match Mat3R::new(record_lattice).inverse() {
        Ok(inverse) => inverse.transpose(),
        // A singular conventional basis is a property of the record; anything
        // else is an exact-arithmetic limit and must not be reported as one.
        Err(SubductionError::SingularMatrix) => {
            return Err(GeometryError::SingularConventionalBasis { uni: geometry.uni });
        }
        Err(_) => return Err(GeometryError::Arithmetic),
    };
    for generator in full_translation_lattice(set) {
        let image = map
            .checked_mul_vector(&Vec3R::new(generator))
            .map_err(|_| GeometryError::Arithmetic)?;
        let in_record = record_membership
            .checked_mul_vector(&image)
            .map_err(|_| GeometryError::Arithmetic)?;
        if !in_record.as_array().iter().all(|value| value.is_integer()) {
            return Ok(false);
        }
    }
    let map_inverse = match map.inverse() {
        Ok(inverse) => inverse,
        Err(SubductionError::SingularMatrix) => return Ok(false),
        Err(_) => return Err(GeometryError::Arithmetic),
    };
    for row in &record_lattice {
        let pullback = map_inverse
            .checked_mul_vector(&Vec3R::new(*row))
            .map_err(|_| GeometryError::Arithmetic)?;
        if !in_lattice(*pullback.as_array(), &translation_lattice(set))
            .map_err(|_| GeometryError::Arithmetic)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether the map realises the record's lattice **exactly**: `A L_M = L_record`.
///
/// Both inclusions are needed.  Forward: every generator of the full source
/// lattice ([`full_translation_lattice`]) maps into the record's row lattice
/// (`(M^-1)^T v` integral, the row convention).  Reverse: every row of the
/// record's lattice pulls back through `A^-1` into `L_M` (membership modulo the
/// cell translations, via [`in_lattice`]).  Round 20 replaced the earlier test,
/// which looked only at the listed generators in the forward direction and was
/// therefore vacuously true whenever the set lists no pure translation.
fn realises_record_lattice(
    map: &Mat3R,
    source_generators: &[[Rat; 3]],
    source_lattice: &[[Rat; 3]],
    record_lattice: &[[Rat; 3]; 3],
    record_inverse_transposed: &Mat3R,
) -> bool {
    let in_record_lattice = |vector: [Rat; 3]| -> bool {
        record_inverse_transposed
            .checked_mul_vector(&Vec3R::new(vector))
            .map(|image| image.as_array().iter().all(|value| value.is_integer()))
            .unwrap_or(false)
    };
    let forward = source_generators.iter().all(|generator| {
        match map.checked_mul_vector(&Vec3R::new(*generator)) {
            Ok(image) => in_record_lattice(*image.as_array()),
            Err(_) => false,
        }
    });
    if !forward {
        return false;
    }
    let Ok(inverse) = map.inverse() else {
        return false;
    };
    record_lattice
        .iter()
        .all(|row| match inverse.checked_mul_vector(&Vec3R::new(*row)) {
            Ok(pullback) => in_lattice(*pullback.as_array(), source_lattice).unwrap_or(false),
            Err(_) => false,
        })
}

/// The parent's operations, rotation set and translation-lattice generators,
/// all exact -- the shared prologue of every embedding query.
///
/// The parent's operation table lists only the **non-zero** pure translations;
/// the cell translations `Z^3` stay implicit everywhere, exactly as in
/// [`canonical_translation`] and [`in_lattice`].
struct ParentFrame {
    /// `(rotation, translation)` in database order.
    operations: Vec<(Mat3I, [Rat; 3])>,
    /// The distinct rotations, in first-seen order.
    rotations: Vec<Mat3I>,
    /// The non-zero pure translations (generators of `L_P` over `Z^3`).
    lattice: Vec<[Rat; 3]>,
}

fn parent_frame(sg: u8) -> Result<ParentFrame, GeometryError> {
    let parent = crate::irrep::query::symmetry_operations_of(sg)
        .map_err(|_| GeometryError::UnreadableParent { sg })?;
    let mut operations: Vec<(Mat3I, [Rat; 3])> = Vec::with_capacity(parent.operations.len());
    let mut rotations: Vec<Mat3I> = Vec::new();
    let mut lattice: Vec<[Rat; 3]> = Vec::new();
    for operation in &parent.operations {
        let translation = [
            exact_fraction(operation.translation[0])
                .ok_or(GeometryError::ParentTranslationOffGrid { sg })?,
            exact_fraction(operation.translation[1])
                .ok_or(GeometryError::ParentTranslationOffGrid { sg })?,
            exact_fraction(operation.translation[2])
                .ok_or(GeometryError::ParentTranslationOffGrid { sg })?,
        ];
        if !rotations.contains(&operation.rotation) {
            rotations.push(operation.rotation);
        }
        if operation.rotation == IDENTITY_ROTATION
            && translation.iter().any(|value| !value.is_zero())
            && !lattice.contains(&translation)
        {
            lattice.push(translation);
        }
        operations.push((operation.rotation, translation));
    }
    Ok(ParentFrame {
        operations,
        rotations,
        lattice,
    })
}

/// The conjugated rotation of every class, or `None` when the map is refuted
/// for a map-only reason: a conjugated rotation that is not integral, or not a
/// parent rotation, can never be placed by any shift.
///
/// All parent operations with a given rotation differ by a lattice vector --
/// `(R|t_1)(R|t_2)^-1 = (E | t_1 - t_2)` -- so one representative per rotation
/// is the whole coset modulo `L_P`.
/// The rotation screen in checked form.
///
/// Three answers: `Ok(Some(..))` when every conjugated rotation is integral and
/// a parent rotation; `Ok(None)` when some conjugated rotation is **not
/// integral** (a proof: no shift can change a rotation) or lies outside the
/// parent's rotation set (the same proof); `Err` when the exact arithmetic
/// cannot decide it -- an unknown answer that must not be read as a refutation.
fn mapped_class_rotations_checked(
    map: &Mat3R,
    inverse: &Mat3R,
    set: &MagneticOperationSet,
    parent_rotations: &[Mat3I],
) -> Result<Option<Vec<Mat3I>>, GeometryError> {
    let mut mapped = Vec::with_capacity(set.classes.len());
    for operation in &set.classes {
        let rotation = Mat3R::from_ints(operation.rotation);
        let product = map
            .checked_mul(&rotation)
            .and_then(|product| product.checked_mul(inverse))
            .map_err(|_| GeometryError::Arithmetic)?;
        let integer = match product.to_int_matrix() {
            Ok(integer) => integer,
            // Not an integer matrix: no shift can make `A R A^-1` a parent
            // rotation, so this is a proof, not an arithmetic limit.
            Err(SubductionError::NotAnInteger { .. }) => return Ok(None),
            Err(_) => return Err(GeometryError::Arithmetic),
        };
        if !parent_rotations.contains(&integer) {
            return Ok(None);
        }
        mapped.push(integer);
    }
    Ok(Some(mapped))
}

/// The class equations of one map: `(I - R_p(c)) delta == tau_p(c) - A tau_c`,
/// one pair per class, antiunitary classes included (time reversal does not
/// change the translation equation).
fn class_equations(
    map: &Mat3R,
    set: &MagneticOperationSet,
    mapped_rotations: &[Mat3I],
    parent_operations: &[(Mat3I, [Rat; 3])],
) -> Result<Vec<(Mat3I, [Rat; 3])>, GeometryError> {
    let mut equations = Vec::with_capacity(set.classes.len());
    for (index, operation) in set.classes.iter().enumerate() {
        let rotation = mapped_rotations[index];
        let Some((_, parent_translation)) = parent_operations
            .iter()
            .find(|(candidate, _)| *candidate == rotation)
        else {
            // Unreachable: the rotation screen accepted every mapped rotation.
            return Err(GeometryError::Arithmetic);
        };
        let moved = map
            .checked_mul_vector(&Vec3R::new(operation.translation))
            .map_err(|_| GeometryError::Arithmetic)?;
        let mut matrix = [[0i32; 3]; 3];
        let mut target = [Rat::ZERO; 3];
        for row in 0..3 {
            for column in 0..3 {
                let identity = if row == column { 1 } else { 0 };
                matrix[row][column] = identity - rotation[row][column];
            }
            target[row] = parent_translation[row]
                .checked_sub(moved.get(row))
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        equations.push((matrix, target));
    }
    Ok(equations)
}

/// The complete shift solution set of one fixed candidate map (round 24).
///
/// For a fixed map the set of shifts that place **every** class on a parent
/// operation is a finite union of cosets of the common continuous kernel; this
/// returns one exact representative per component plus a basis of that kernel,
/// instead of the finite `{+-1, +-1/2}` walk the search used to run.  The
/// returned shift set is exactly the set of shifts for which
/// [`verify_embedding`] accepts `(map, shift)`:
///
/// * a map whose source translation lattice does not land in the parent's
///   lattice, whose conjugated class rotation is not a parent rotation or is not
///   integral, or whose matrix is singular, has an **empty** solution set -- a
///   proof, because no shift can change any of those;
/// * otherwise the class congruences are decided exactly, and an empty result
///   means they are jointly inconsistent (which can happen even when every
///   class is placeable on its own);
/// * otherwise every listed representative is a genuine solution, and every
///   solution is `representative + v + lambda` with `v` in the span of `kernel`
///   and `lambda` in the parent's translation lattice.
///
/// `Err(GeometryError::ShiftEnumerationUnavailable)` and
/// `Err(GeometryError::Arithmetic)` are **unknown** answers (exact-arithmetic
/// limits or the component cap), never "no shift": every checked step of the
/// prologue distinguishes a proof from an arithmetic limit instead of folding
/// the latter into the empty set.
pub fn embedding_shift_components(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
    map: Mat3R,
) -> Result<ShiftComponents, GeometryError> {
    let parent = parent_frame(geometry.parent_sg)?;
    let empty = ShiftComponents {
        representatives: Vec::new(),
        kernel: Vec::new(),
    };
    // Every step below distinguishes a **proof** (the map is refuted) from an
    // arithmetic **limit** (the answer is unknown): an empty `ShiftComponents`
    // is documented as "no shift makes this map an embedding", so it must never
    // absorb an unknown answer.  The search uses the `bool` wrappers for the same
    // conditions because its screen counters are pinned measurements; this door
    // is the conservative reading.
    if !source_lattice_lands_in_parent_checked(&map, set, &parent.lattice)? {
        return Ok(empty);
    }
    let inverse = match map.inverse() {
        Ok(inverse) => inverse,
        // A singular map is not a change of basis: no shift can make it an
        // embedding.  Any other failure is an arithmetic limit.
        Err(SubductionError::SingularMatrix) => return Ok(empty),
        Err(_) => return Err(GeometryError::Arithmetic),
    };
    let Some(mapped_rotations) =
        mapped_class_rotations_checked(&map, &inverse, set, &parent.rotations)?
    else {
        return Ok(empty);
    };
    let equations = class_equations(&map, set, &mapped_rotations, &parent.operations)?;
    let solution = complete_shift::solve_shift_components(&equations, &parent.lattice)
        .map_err(|_| GeometryError::ShiftEnumerationUnavailable { uni: geometry.uni })?;
    Ok(ShiftComponents {
        representatives: solution.representatives,
        kernel: solution.kernel,
    })
}

/// Check one claimed embedding `(map, shift)` against the parent's operations.
///
/// This is the same condition the search uses, exposed so that a caller can
/// verify a map it obtained elsewhere -- for example a small explicit embedding
/// found by an independent search, which the family above may not contain.
///
/// The check has two halves, and both are needed:
///
/// * the image of the group's **full** translation lattice `L_M`
///   ([`full_translation_lattice`]) must lie in the parent's lattice.  Per-class
///   placement does not imply this, and round 20 measured the gap;
/// * every class representative (antiunitary ones included) must land on a
///   parent operation.
///
/// Success therefore means the tabulated group maps into the parent's grey
/// group, not merely that its representatives do.
pub fn verify_embedding(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
    map: Mat3R,
    shift: [Rat; 3],
) -> Result<(), GeometryError> {
    let parent = crate::irrep::query::symmetry_operations_of(geometry.parent_sg).map_err(|_| {
        GeometryError::UnreadableParent {
            sg: geometry.parent_sg,
        }
    })?;
    let mut parent_operations: Vec<(Mat3I, [Rat; 3])> = Vec::new();
    let mut parent_rotations: Vec<Mat3I> = Vec::new();
    let mut parent_lattice: Vec<[Rat; 3]> = Vec::new();
    for operation in &parent.operations {
        let translation = [
            exact_fraction(operation.translation[0]).ok_or(
                GeometryError::ParentTranslationOffGrid {
                    sg: geometry.parent_sg,
                },
            )?,
            exact_fraction(operation.translation[1]).ok_or(
                GeometryError::ParentTranslationOffGrid {
                    sg: geometry.parent_sg,
                },
            )?,
            exact_fraction(operation.translation[2]).ok_or(
                GeometryError::ParentTranslationOffGrid {
                    sg: geometry.parent_sg,
                },
            )?,
        ];
        if !parent_rotations.contains(&operation.rotation) {
            parent_rotations.push(operation.rotation);
        }
        if operation.rotation == IDENTITY_ROTATION
            && translation.iter().any(|value| !value.is_zero())
            && !parent_lattice.contains(&translation)
        {
            parent_lattice.push(translation);
        }
        parent_operations.push((operation.rotation, translation));
    }
    // Complete containment first: the source lattice is a property of the map
    // alone, and no per-class equation below can replace this gate.
    if !source_lattice_lands_in_parent_checked(&map, set, &parent_lattice)? {
        return Err(GeometryError::SourceLatticeNotInParent { uni: geometry.uni });
    }
    let inverse = map
        .inverse()
        .map_err(|_| GeometryError::SingularConventionalBasis { uni: geometry.uni })?;
    let mut mapped_rotations = Vec::with_capacity(set.classes.len());
    for operation in &set.classes {
        let rotation = Mat3R::from_ints(operation.rotation);
        let conjugated = map
            .checked_mul(&rotation)
            .and_then(|product| product.checked_mul(&inverse))
            .map_err(|_| GeometryError::Arithmetic)?;
        let mapped = conjugated
            .to_int_matrix()
            .map_err(|_| GeometryError::Arithmetic)?;
        if !parent_rotations.contains(&mapped) {
            return Err(GeometryError::NotInParent {
                uni: geometry.uni,
                rotation: operation.rotation,
            });
        }
        mapped_rotations.push(mapped);
    }
    if !embedding_is_consistent(
        &map,
        &shift,
        set,
        &mapped_rotations,
        &parent_operations,
        &parent_lattice,
    )? {
        return Err(GeometryError::NotInParent {
            uni: geometry.uni,
            rotation: IDENTITY_ROTATION,
        });
    }
    Ok(())
}

/// How a candidate map was built.
///
/// The names below state the matrix the family actually multiplies by, which
/// round 20 found was not what the older docs said: the "dual" variants are the
/// **inverse** of the frame basis, because the doubles are taken from
/// `B = W^T` as `(B^T)^-T = B^-1`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EmbeddingKind {
    /// A signed permutation of the tabulated coordinates, no basis conversion:
    /// `S . I`.
    Plain,
    /// A signed permutation composed with the transpose of the record's
    /// primitive-frame basis: `S . W_prim^T`.
    PrimitiveBasis,
    /// A signed permutation composed with the transpose of the record's
    /// conventional-frame basis: `S . W_conv^T`.
    ConventionalBasis,
    /// A signed permutation composed with the **inverse** of the primitive-frame
    /// basis: `S . W_prim^-1` (the variant name is historical; see the enum
    /// doc, and note the pinned counts belong to this matrix).
    PrimitiveDual,
    /// A signed permutation composed with the **inverse** of the
    /// conventional-frame basis: `S . W_conv^-1`.
    ConventionalDual,
    /// The ordinary H-to-parent map composed with the identified M-to-H frame.
    /// Produced only by [`embed_via_unitary_subgroup`], not the R24 map family.
    UnitarySubgroup,
}

/// One way to place the tabulated magnetic group inside the parent.
#[derive(Clone, Debug, PartialEq)]
pub struct ParentEmbedding {
    /// The map applied to the magnetic group's coordinates.
    ///
    /// Rational in general: the dual bases carry denominators (the record's
    /// conventional basis has index 1/3 for a rhombohedral parent, say).
    pub map: Mat3R,
    /// The exact origin shift, in the parent's conventional coordinates,
    /// reduced modulo the parent's translation lattice.
    ///
    /// Round 24: this is the representative of one **discrete component** of the
    /// map's complete solution set ([`embedding_shift_components`]): shifting it
    /// by any element of the component's common kernel `V` gives the same image
    /// operation set, and by any parent lattice vector the same class.  The
    /// components are the finitely many cosets of `V + L_P`, so the entry list is
    /// finite and normative; the walk it replaced sampled arbitrary points of the
    /// same cosets instead.
    pub shift: [Rat; 3],
    /// How the map was built.
    pub kind: EmbeddingKind,
    /// Whether the map realises the record's lattice **exactly**:
    /// `A . L_M = L_record`, where `L_M` is the group's full translation lattice
    /// ([`full_translation_lattice`]: the cell translations `Z^3` plus the
    /// listed non-zero unprimed pure translations) and `L_record` is the
    /// record's `basis` in the parent's conventional frame.
    ///
    /// Both inclusions are tested: the source generators map into the record's
    /// row lattice (`(M^-1)^T v` integral), and the record's lattice rows pull
    /// back through `A^-1` into `L_M`.  Round 20: the earlier test looked only
    /// at the listed generators in the forward direction, so it was vacuously
    /// true for a set that lists no pure translation -- SG 1 `Z1/P1` UNI 3 was
    /// reported as realising its lattice while `A e_1 = e_3` is not even in it.
    /// The forward-only test was also wrong in the other direction: a map whose
    /// image lattice is a proper sublattice of the record's (synthetic SG 5
    /// UNI 20 with a record lattice doubled, `W = 2I`, `A = I`) was reported as
    /// failing to realise it although the whole group is contained.
    ///
    /// This flag is separate from [`verify_embedding`]: a map can
    /// place the group in the parent (containment holds) without realising the
    /// record's lattice, and the search reports both readings instead of
    /// collapsing them.
    pub realises_record_lattice: bool,
}

/// Propose an exact magnetic embedding through an existing ordinary embedding of H.
///
/// If `T` maps the identified unitary subgroup H into the parent and `Q` is
/// the M-to-H data-Hall frame, the candidate linear map is `A = T Q`.
/// The floating-point frame is only a *proposal*: its entries must lie on the
/// existing small rational grid, and every returned map is checked with exact
/// arithmetic against the magnetic group's complete operation set and lattice.
///
/// The ordinary embedding must have the same parent and the identified H space
/// group. The caller selects it (for isotropy records, normally by the same
/// condensing irrep and direction label); this function does not assert that the
/// two order-parameter directions are identical. It certifies grey-parent
/// containment and `A L_M = L_record`, not an order-parameter stabilizer.
///
/// Origins are solved anew, simultaneously for **all** magnetic classes using
/// [`embedding_shift_components`]. In particular, neither an ordinary origin
/// nor a floating-point frame shift is blindly reused. The returned list has one
/// verified witness per discrete origin component; no basis uniqueness is claimed.
///
/// An empty list only rejects this one derived map. An error means the bridge
/// was unavailable or verification was undecided; neither disproves other maps.
/// This does not change [`parent_embedding_verdict`] or its restricted family
/// absence certificates. No representation or non-Γ phase transport is performed.
pub fn embed_via_unitary_subgroup(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
    unitary_embedding: &crate::irrep::subduction::SubgroupEmbedding,
) -> Result<Vec<ParentEmbedding>, GeometryError> {
    let unavailable = |reason: String| GeometryError::UnitaryBridgeUnavailable {
        uni: geometry.uni,
        reason,
    };
    if geometry.uni != set.uni || geometry.parent_sg != unitary_embedding.parent_sg() {
        return Err(unavailable(
            "magnetic UNI or parent identity differs".into(),
        ));
    }
    verify_against_database(set)
        .and_then(|()| verify_group(set))
        .map_err(|error| unavailable(format!("invalid magnetic operation set: {error}")))?;
    let h = crate::irrep::corep::cached_unitary_subgroup_with_hall(set.uni)
        .ok_or_else(|| unavailable("unitary subgroup identification unavailable".into()))?;
    if h.sg != usize::from(unitary_embedding.subgroup_sg()) {
        return Err(unavailable(format!(
            "identified H is SG {}, supplied ordinary subgroup is SG {}",
            h.sg,
            unitary_embedding.subgroup_sg()
        )));
    }
    let frame = h
        .msg_to_data
        .as_ref()
        .ok_or_else(|| unavailable("M-to-H data-Hall frame unavailable".into()))?;
    let mut q = [[Rat::ZERO; 3]; 3];
    for (row, values) in q.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = exact_fraction(frame.basis[row][column])
                .ok_or_else(|| unavailable("M-to-H linear frame off the rational grid".into()))?;
        }
    }
    let map = unitary_embedding
        .transform()
        .matrix()
        .checked_mul(&Mat3R::new(q))
        .map_err(|_| GeometryError::Arithmetic)?;
    if !realises_record_lattice_checked(geometry, set, &map)? {
        return Ok(Vec::new());
    }
    let components = embedding_shift_components(geometry, set, map)?;
    let mut out = Vec::with_capacity(components.representatives.len());
    for shift in components.representatives {
        verify_embedding(geometry, set, map, shift)?;
        out.push(ParentEmbedding {
            map,
            shift,
            kind: EmbeddingKind::UnitarySubgroup,
            realises_record_lattice: true,
        });
    }
    Ok(out)
}

/// What the embedding search decided for one record (round 22, extended in round
/// 23a and round 24).
#[derive(Clone, Debug)]
pub struct ParentSettingVerdict {
    /// Every verified embedding of the tabulated group into the parent: one
    /// entry per `(kind, map, discrete shift component)` triple.  Round 24
    /// replaced the `{+-1, +-1/2}` shift walk with the **complete** solution set
    /// of each map ([`embedding_shift_components`]), so the entry count is the
    /// number of components -- shifts that differ by an element of the common
    /// continuous kernel produce the same image operation set and are represented
    /// once, instead of being sampled at whatever points the walk happened to
    /// reach.  The same `(map, component)` reached through several bases of the
    /// family counts once per `kind`; `distinct_record_map_shift` in the census
    /// is the key-collapsed reading.
    pub embeddings: Vec<ParentEmbedding>,
    /// `true` when `embeddings` is empty **and** every candidate map was excluded
    /// by an exact condition of an embedding -- see the counters below.  All of
    /// them are consequences of what an embedding is:
    ///
    /// * the image of the full source translation lattice must lie in the
    ///   parent's lattice (round 20's gate);
    /// * every class rotation must be a parent rotation;
    /// * every class must be **placeable**: it must land on *some* parent
    ///   operation with the mapped rotation, so the anchor congruence
    ///   `(I - R_p) delta == tau_p - A tau (mod L_parent)` has to be solvable for
    ///   one of them (round 23a's exact joint test, antiunitary classes
    ///   included);
    /// * the class congruences must be **jointly** consistent (round 24): a map
    ///   that can place every class on its own can still have no single shift
    ///   that places them all at once, and the complete solution set decides
    ///   exactly that.
    ///
    /// So `true` means "**no member of this candidate family embeds the group**",
    /// which is exactly the stage-2b acceptance criterion.  The residual is
    /// `candidates_without_witness`: maps whose solution set could not be decided
    /// (an arithmetic limit or the component cap of the complete solver), or
    /// whose constructed representative failed the consistency check.  Round 24's
    /// complete solver makes that residual **0 in production** -- every other
    /// map is refuted by one of the four exact conditions or yields at least one
    /// verified entry -- so an empty result **with `candidates_without_witness
    /// == 0`** (which is what this flag asserts) is a decision for the family
    /// rather than a statement about a finite walk.  Round 24 extended that rule
    /// to the screens themselves: an arithmetic limit in the source-lattice or
    /// rotation gate leaves a candidate undecided and increments the residual,
    /// instead of being credited as a refutation.
    ///
    /// Maps that could not even be built (an arithmetic limit while inverting a
    /// family basis, or while multiplying a signed permutation by it) count in
    /// the residual as well.  A **singular** candidate map is the one skip that
    /// is genuinely refuted -- no shift makes a singular map a change of basis,
    /// and [`verify_embedding`] rejects it -- so it contributes to no counter and
    /// does not make this flag unsound.
    ///
    /// Scope fences.  The candidate family is still a family (`{S . B}`), so
    /// `true` never claims that the record is not a subgroup of the parent --
    /// only that no member of this family embeds it, while some map *outside* the
    /// family may still embed it (SG 3 UNI 24 is the standing witness).  `false`
    /// means "not settled here": either there is a witness, or some candidate is
    /// undecided (`candidates_without_witness`); it never means "no member
    /// embeds".
    /// Round 22's withdrawn attempt is history: it read the anchor congruence per
    /// normal vector -- a necessary condition only -- and, worse, never ran on its
    /// certified rows.  The replacement is the joint congruence: the round-22
    /// review measured that it lifts the certified zero rows from 5,649 of 8,294
    /// to at least 8,106, and round 24 lifts the residual to the whole corpus.
    pub absence_is_certified: bool,
    /// Maps whose image of the full source translation lattice leaves the
    /// parent's lattice (round 20's gate) -- impossible for an embedding.
    pub excluded_by_source_lattice: usize,
    /// Maps that conjugate some class rotation out of the parent's rotation set.
    pub excluded_by_rotation: usize,
    /// Maps refuted by per-class placeability (round 23a): some class cannot be
    /// placed on any parent operation with its mapped rotation, because the exact
    /// joint congruence has no solution there.  Counting is exclusive in the
    /// order the conditions are checked, so this counts only maps that
    /// passed both screens above.
    pub excluded_by_anchor_congruence: usize,
    /// Maps refuted by the **joint** congruence (round 24): they pass the three
    /// screens above, but the complete solution set of all class equations at
    /// once is empty.  Counting is exclusive, so this counts only maps that
    /// passed the three screens.  (The per-class screen is necessary, not
    /// sufficient: a class can be placeable only through a lattice-shifted
    /// right-hand side while no single shift places all of them.)
    pub excluded_by_joint_congruence: usize,
    /// Maps that passed all four exact conditions above and still produced no
    /// verified embedding.  Round 24: this is the *undecided* residual -- an
    /// arithmetic limit or component cap inside the complete solver, or a
    /// constructed representative that failed the consistency check -- and it is
    /// 0 on the production corpus.  It is never a refutation and never a witness.
    pub candidates_without_witness: usize,
}

/// Solve for every map in the candidate family that embeds the tabulated group
/// into the parent, with the origin shift **solved exactly** rather than taken
/// from the record, and report whether an empty result is a proof.
///
/// The family is `{S . B}` for the 48 signed permutations `S` and the five bases
/// `B` in `{I, W_prim^T, W_conv^T, W_prim^-1, W_conv^-1}` (the last two are the
/// `EmbeddingKind` "dual" variants; see the enum doc for why they are inverses
/// and not inverse-transposes).  The identity is included (it is the
/// `S = I`, `B = I` member), because round 17 measured that some records are
/// tabulated directly in the parent's frame and a family built only from the
/// record's basis misses them entirely.
///
/// Acceptance is the same condition [`verify_embedding`] uses: the image of the
/// **full** source translation lattice must lie in the parent's lattice (round
/// 20), and every class representative must land on a parent operation with the
/// solved shift.  Round 20's counterexample -- SG 143 `M1/S1` UNI 1 with
/// `A = W^-1` -- is no longer returned.  Before the complete solution set, every
/// candidate map is also screened for per-class placeability (round 23a).
///
/// Round 22 turned the empty result from "this walk found nothing" into a
/// decidable statement for the family; round 23a added the third exact
/// exclusion -- per-class placeability, decided by the joint congruence rather
/// than by round 22's withdrawn per-normal test.  Round 24 removed the finite
/// `{+-1, +-1/2}` shift walk altogether: each surviving map gets its **complete**
/// shift solution set ([`embedding_shift_components`]), so a map that passes the
/// three screens and has no single shift is now refuted exactly
/// (`excluded_by_joint_congruence`) instead of leaving an unsettled residual.
/// An empty result therefore no longer depends on a walk at all; see
/// [`ParentSettingVerdict`].  Use [`embed_in_parent_conventional`] when only the
/// witnesses are wanted.
pub fn parent_embedding_verdict(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Result<ParentSettingVerdict, GeometryError> {
    search_embeddings(geometry, set)
}

/// The verified embeddings of `set` in the parent, in the order the family
/// generates them: one entry per `(kind, map, discrete shift component)` triple,
/// so the same `(map, component)` reached through several bases of the family
/// appears once per `kind` (the `distinct_record_map_shift` census key collapses
/// those).  The components come from the **complete** solution set of each map
/// (round 24, [`embedding_shift_components`]), not from a finite shift walk;
/// shifts that differ by an element of a map's common continuous kernel give the
/// same image operation set and are represented once.
///
/// An empty result is **not by itself a proof**: read
/// [`ParentSettingVerdict::absence_is_certified`] on the verdict instead.  Empty
/// plus `candidates_without_witness == 0` (equivalently
/// `absence_is_certified`) means every candidate of the family was refuted by an
/// exact condition; an arithmetic limit or the component cap inside the complete
/// solver leaves a candidate undecided, and then the empty list is "not settled
/// here".  Either way it is not a proof that the record is not a subgroup of the
/// parent: the candidate family is `{S . B}`.
///
/// The `realises_record_lattice` flag separates "contained" from "the record's
/// lattice is realised exactly".
pub fn embed_in_parent_conventional(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Result<Vec<ParentEmbedding>, GeometryError> {
    Ok(parent_embedding_verdict(geometry, set)?.embeddings)
}

/// The search itself: [`parent_embedding_verdict`] is the public door.
fn search_embeddings(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Result<ParentSettingVerdict, GeometryError> {
    let parent = parent_frame(geometry.parent_sg)?;
    let parent_operations: &[(Mat3I, [Rat; 3])] = &parent.operations;
    let parent_rotations: &[Mat3I] = &parent.rotations;
    let parent_lattice: &[[Rat; 3]] = &parent.lattice;
    // `realises_record_lattice` asks for **equality** `A L_M = L_record` (round
    // 20).  The forward half maps the full source translation lattice
    // (`Z^3` + the listed pure translations); the backward half pulls the
    // record's lattice rows back through `A^-1` and asks for membership in
    // `L_M`.  The earlier version tested only the listed generators in the
    // forward direction, which is vacuously true whenever the set lists none,
    // and never tested the reverse.
    let source_generators: Vec<[Rat; 3]> = full_translation_lattice(set);
    let source_lattice: Vec<[Rat; 3]> = translation_lattice(set);
    let record_lattice = basis_in_parent_conventional(geometry)?;
    // The record's `basis` rows are the lattice vectors, so membership is
    // `(M^-1)^T v` integral -- NOT `M^-1 v`, which tests the column lattice
    // (round 15i fixed exactly this confusion in `spans_same_lattice_as`, and
    // round 15m caught it again here).
    let record_inverse_transposed = match Mat3R::new(record_lattice).inverse() {
        Ok(inverse) => inverse.transpose(),
        // A singular conventional basis is a property of the record; anything
        // else is an exact-arithmetic limit and must not be reported as one.
        Err(SubductionError::SingularMatrix) => {
            return Err(GeometryError::SingularConventionalBasis { uni: geometry.uni });
        }
        Err(_) => return Err(GeometryError::Arithmetic),
    };

    let mut bases: Vec<(EmbeddingKind, Mat3R)> = Vec::new();
    bases.push((EmbeddingKind::Plain, Mat3R::identity()));
    let transpose_of = |matrix: &Mat3R| -> Mat3R {
        let mut columns = [[Rat::ZERO; 3]; 3];
        for (row, values) in columns.iter_mut().enumerate() {
            for (column, value) in values.iter_mut().enumerate() {
                *value = matrix.row(column)[row];
            }
        }
        Mat3R::new(columns)
    };
    bases.push((
        EmbeddingKind::PrimitiveBasis,
        transpose_of(&Mat3R::new(geometry.basis)),
    ));
    bases.push((
        EmbeddingKind::ConventionalBasis,
        transpose_of(&Mat3R::new(record_lattice)),
    ));
    // The "dual" members: `base` is already `W^T`, so `transpose(inverse(base))`
    // is `(W^T)^-T = W^-1`, i.e. the **inverse** of the frame basis.  The variant
    // names and the older docs called this the inverse-transpose (round 20
    // flagged the mismatch); the pinned counts belong to the matrix this code
    // actually builds, so the doc was corrected rather than the family.
    // A family member that cannot even be built is an **undecided** candidate,
    // never a silent omission: the family is part of the claim (round 24).
    let mut unbuilt_family_members = 0usize;
    let mut duals: Vec<(EmbeddingKind, Mat3R)> = Vec::new();
    for (kind, base) in &bases {
        if *kind == EmbeddingKind::Plain {
            continue;
        }
        let Ok(inverse) = base.inverse() else {
            unbuilt_family_members += 1;
            continue;
        };
        let dual = transpose_of(&inverse);
        let dual_kind = match kind {
            EmbeddingKind::PrimitiveBasis => EmbeddingKind::PrimitiveDual,
            EmbeddingKind::ConventionalBasis => EmbeddingKind::ConventionalDual,
            EmbeddingKind::Plain => continue,
            EmbeddingKind::PrimitiveDual
            | EmbeddingKind::ConventionalDual
            | EmbeddingKind::UnitarySubgroup => continue,
        };
        duals.push((dual_kind, dual));
    }
    bases.extend(duals);

    let mut found: Vec<ParentEmbedding> = Vec::new();
    // Round 22/22b: an empty result is a proof when **every** candidate map was
    // excluded by an exact necessary condition.  The exclusions are counted
    // separately because the round-22 review found that the first version
    // credited the anchor congruence for proofs that the screens had made.
    //
    // Round 24 (F1): the fourth counter is the **joint** refutation.  The
    // `{+-1, +-1/2}` shift walk is gone; each surviving map is handed to
    // [`embedding_shift_components`], which decides its whole shift solution set
    // exactly.  An empty solution set can only mean the class congruences are
    // jointly inconsistent -- which can happen even when every class is
    // placeable on its own -- and that is a proof for this map, so it is
    // counted.  `candidates_without_witness` therefore stops meaning "the walk
    // found nothing" and becomes the honest residual: maps whose solution set
    // could not be decided (arithmetic limit or component cap) or whose
    // constructed representatives failed the consistency check.  In production
    // it is 0, and `absence_is_certified` is then a decision for the family.
    let mut excluded_by_source_lattice = 0usize;
    let mut excluded_by_rotation = 0usize;
    let mut excluded_by_anchor_congruence = 0usize;
    let mut excluded_by_joint_congruence = 0usize;
    // The residual starts with the family members that could not be built (an
    // arithmetic limit while inverting a basis is not a refutation of anything).
    let mut candidates_without_witness = unbuilt_family_members;
    for (kind, base) in &bases {
        for permutation in signed_permutations() {
            let Ok(map) = Mat3R::from_ints(permutation).checked_mul(base) else {
                // The map could not be built: an arithmetic limit, not a
                // refutation (round 24: an unbuilt candidate must not vanish, or
                // a row could be certified on an incomplete family).
                candidates_without_witness += 1;
                continue;
            };
            let inverse = match map.inverse() {
                Ok(inverse) => inverse,
                // A singular map is not a change of basis, and the public
                // verifier rejects it, so no shift can make it an embedding:
                // refuted, not undecided.
                Err(SubductionError::SingularMatrix) => continue,
                // Anything else is an exact-arithmetic limit.
                Err(_) => {
                    candidates_without_witness += 1;
                    continue;
                }
            };
            // Complete containment (round 20): the image of the full source
            // translation lattice must lie in the parent's lattice.  The class
            // loop below cannot see this, and the condition depends on the map
            // only, so it runs once here rather than once per shift.
            match source_lattice_lands_in_parent_checked(&map, set, parent_lattice) {
                // A proof: the image leaves the parent lattice.
                Ok(false) => {
                    excluded_by_source_lattice += 1;
                    continue;
                }
                // An arithmetic limit is **unknown**, never a refutation
                // (round 24: the invariant "arithmetic failure is not no
                // solution" applies to this production verdict too).
                Err(_) => {
                    candidates_without_witness += 1;
                    continue;
                }
                Ok(true) => {}
            }
            // Rotations must land in the parent's rotation set.
            let mapped_rotations =
                match mapped_class_rotations_checked(&map, &inverse, set, parent_rotations) {
                    // A proof: a conjugated rotation that is not integral or not
                    // a parent rotation can never be placed by any shift.
                    Ok(None) => {
                        excluded_by_rotation += 1;
                        continue;
                    }
                    Err(_) => {
                        candidates_without_witness += 1;
                        continue;
                    }
                    Ok(Some(mapped_rotations)) => mapped_rotations,
                };
            // Round 23a: per-class placeability, refuted **before** the shift
            // solution is computed.  An embedding has to place every class --
            // antiunitary classes included, whose translation equation is the
            // same one -- on *some* parent operation with the mapped rotation, so
            // a class that fits on none of them refutes this map for good.
            // `anchor_congruence_is_feasible` decides each of those equations
            // exactly; this is the joint congruence, not round 22's per-normal
            // necessary condition.
            //
            // It is deliberately **not** used to filter the complete solver
            // below: a class can be placeable only through a **lattice-shifted**
            // right-hand side, and the per-class screen is only a necessary
            // condition.  The screen is kept because it is an independent exact
            // refutation with its own pinned counter.
            let mut placeable = true;
            for (index, anchor) in set.classes.iter().enumerate() {
                let rotation = mapped_rotations[index];
                let Ok(moved) = map.checked_mul_vector(&Vec3R::new(anchor.translation)) else {
                    // An arithmetic failure is not a refutation: leave the class
                    // to the complete solver, which fails the same way and
                    // reports no witness.
                    continue;
                };
                let moved = *moved.as_array();
                let mut difference = [[0i32; 3]; 3];
                for row in 0..3 {
                    for column in 0..3 {
                        let identity = if row == column { 1 } else { 0 };
                        difference[row][column] = identity - rotation[row][column];
                    }
                }
                let mut any = false;
                for (parent_rotation, parent_translation) in parent_operations {
                    if *parent_rotation != rotation {
                        continue;
                    }
                    let mut target = [Rat::ZERO; 3];
                    let mut exact = true;
                    for axis in 0..3 {
                        match parent_translation[axis].checked_sub(moved[axis]) {
                            Ok(value) => target[axis] = value,
                            Err(_) => {
                                exact = false;
                                break;
                            }
                        }
                    }
                    if !exact {
                        // Round 24: an arithmetic limit while forming this
                        // right-hand side is **not** a refutation of the map.
                        // Treat the class as placeable and let
                        // `class_equations` report the limit (or the complete
                        // solver book it as undecided), exactly as the `moved`
                        // failure above does.  Continuing to the next matching
                        // parent operation instead would leave `any == false`
                        // and credit `excluded_by_anchor_congruence` for a
                        // computation that never finished.
                        any = true;
                        break;
                    }
                    // An arithmetic limit inside the decider is **not** a
                    // refutation (round-23 review F2): treat the class as
                    // placeable and let the complete solver decide, exactly as
                    // the `moved` failure above does.  Only an exact
                    // `Ok(false)` refutes the map.
                    match anchor_congruence_is_feasible(&difference, &target, parent_lattice) {
                        Ok(true) => {
                            any = true;
                            break;
                        }
                        Ok(false) => {}
                        Err(_) => {
                            any = true;
                            break;
                        }
                    }
                }
                if !any {
                    placeable = false;
                    break;
                }
            }
            if !placeable {
                excluded_by_anchor_congruence += 1;
                continue;
            }
            // Round 24 (F1): the **complete** shift solution set of this map.
            // Every class equation is decided at once -- the walk that sampled
            // one anchor's affine set with the step set `{+-1, +-1/2}` could not
            // reach a lattice-shifted branch (SG 56 UNI 119) or a quarter-denominator
            // component (SG 89 UNI 111).
            //
            // The solution set is a finite union of cosets of the common kernel,
            // so one entry per component is the canonical finite answer; the
            // continuous directions are described by the kernel, never sampled.
            // Each representative is re-checked through `embedding_is_consistent`
            // -- the same condition the public verifier applies -- so a returned
            // entry can never be a claim the verifier rejects.
            let equations = class_equations(&map, set, &mapped_rotations, parent_operations)?;
            match complete_shift::solve_shift_components(&equations, parent_lattice) {
                // An arithmetic limit or the component cap: the answer is
                // **unknown**, so the map is unsettled, never refuted (round-23
                // review F2).
                Err(_) => candidates_without_witness += 1,
                // Passing all three screens and still having no single shift is a
                // proof for this map: the class congruences are jointly
                // inconsistent.
                Ok(solution) if solution.representatives.is_empty() => {
                    excluded_by_joint_congruence += 1;
                }
                Ok(solution) => {
                    let mut witnessed = false;
                    for shift in solution.representatives {
                        if !embedding_is_consistent(
                            &map,
                            &shift,
                            set,
                            &mapped_rotations,
                            parent_operations,
                            parent_lattice,
                        )
                        .unwrap_or(false)
                        {
                            // A failed recheck cannot certify absence: with no
                            // witness this map remains unsettled below.
                            continue;
                        }
                        // The flag is a function of the MAP (and the record),
                        // not of the candidate's provenance: the map is always
                        // expressed in the parent's conventional frame, so the
                        // record's conventional-frame lattice is the target for
                        // every kind.
                        let entry = ParentEmbedding {
                            map,
                            shift,
                            kind: *kind,
                            realises_record_lattice: realises_record_lattice(
                                &map,
                                &source_generators,
                                &source_lattice,
                                &record_lattice,
                                &record_inverse_transposed,
                            ),
                        };
                        if !found.contains(&entry) {
                            found.push(entry);
                        }
                        witnessed = true;
                    }
                    if !witnessed {
                        candidates_without_witness += 1;
                    }
                }
            }
        }
    }
    Ok(ParentSettingVerdict {
        absence_is_certified: candidates_without_witness == 0 && found.is_empty(),
        excluded_by_source_lattice,
        excluded_by_rotation,
        excluded_by_anchor_congruence,
        excluded_by_joint_congruence,
        candidates_without_witness,
        embeddings: found,
    })
}

// ---------------------------------------------------------------------------
// Round 24 (F1): the `{+-1, +-1/2}` shift walk is **withdrawn**.
//
// `shift_coset(particular, kernel, parent_lattice)` used to walk the affine set
// `particular + ker(I - R_p)` of **one** anchor class with the step set
// `{+-1, +-1/2}`, and the search kept the points the consistency check
// accepted.  Round 23's review (F1) measured what that cannot do, on the public
// API: the walk cannot reach a solution that differs from the anchor's exact
// particular solution by a lattice-shifted branch (SG 56 UNI 119 needs
// `delta_1 = 3/4` where `2 delta_1 == -1/2` is solved exactly by `1/4`), and the
// step set cannot close a component with a quarter denominator (SG 89 UNI 111).
// Both witnesses pass `verify_embedding` while the walk returned nothing, so the
// walk under-reported real embeddings -- and, being a finite sample of a
// continuum, it could never prove that no shift was missed.
//
// The replacement is [`complete_shift::solve_shift_components`]: the stacked
// class congruences are decided by exact integer elimination, the solution set
// is returned as one representative per discrete component plus a basis of the
// common real kernel, and no step grid exists.  `canonical_shift` (the
// one-component reduction the walk used) is gone with it; the
// `embedding_shift_components` door and the search both use the complete set.
// ---------------------------------------------------------------------------

/// Check every class of the set under one `(map, shift)` pair -- the
/// **antiunitary** classes included: time reversal does not change the
/// translation equation, so a primed class `(R|tau)T` needs a parent operation
/// with rotation `A R A^-1` and matching translation, exactly like a unitary
/// one.  Round 15n measured that skipping them let 46% of the returned pairs
/// through.
fn embedding_is_consistent(
    map: &Mat3R,
    shift: &[Rat; 3],
    set: &MagneticOperationSet,
    mapped_rotations: &[Mat3I],
    parent_operations: &[(Mat3I, [Rat; 3])],
    parent_lattice: &[[Rat; 3]],
) -> Result<bool, GeometryError> {
    if mapped_rotations.len() != set.classes.len() {
        return Err(GeometryError::Arithmetic);
    }
    for (operation, rotation) in set.classes.iter().zip(mapped_rotations) {
        // Parent operations with the same rotation differ by a parent lattice
        // vector, so one representative suffices. Do not turn arithmetic errors
        // into a false non-containment verdict through an any/all predicate.
        let Some((_, translation)) = parent_operations.iter().find(|(r, _)| r == rotation) else {
            return Ok(false);
        };
        let moved = map
            .checked_mul_vector(&Vec3R::new(operation.translation))
            .map_err(|_| GeometryError::Arithmetic)?;
        let rotated_shift = Mat3R::from_ints(*rotation)
            .checked_mul_vector(&Vec3R::new(*shift))
            .map_err(|_| GeometryError::Arithmetic)?;
        let mut difference = [Rat::ZERO; 3];
        for axis in 0..3 {
            difference[axis] = moved
                .get(axis)
                .checked_add(shift[axis])
                .and_then(|v| v.checked_sub(rotated_shift.get(axis)))
                .and_then(|v| v.checked_sub(translation[axis]))
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        if !in_lattice(difference, parent_lattice).map_err(|_| GeometryError::Arithmetic)? {
            return Ok(false);
        }
    }
    Ok(true)
}
/// Exact Gaussian elimination on a 3x4 augmented system; `None` when
/// inconsistent, otherwise a particular solution and a kernel basis.
/// Solve `system` over the rationals.
///
/// `Ok(None)` means **inconsistent** (no solution); `Err(Arithmetic)` means the
/// elimination hit an arithmetic limit and the system was *not* decided.  The
/// two must not be conflated: the anchor-congruence refutation treats the first
/// as a proof and must treat the second as undecided (round-23 review F2).
/// A particular solution plus a basis of the kernel.
type LinearSolution = ([Rat; 3], Vec<[Rat; 3]>);

fn solve_linear_system(mut system: [[Rat; 4]; 3]) -> Result<Option<LinearSolution>, GeometryError> {
    let mut pivots: Vec<usize> = Vec::new();
    let mut row = 0usize;
    for column in 0..3 {
        let mut found = None;
        for (candidate, values) in system.iter().enumerate().skip(row) {
            if !values[column].is_zero() {
                found = Some(candidate);
                break;
            }
        }
        let Some(found) = found else { continue };
        system.swap(row, found);
        let pivot = system[row][column];
        for entry in system[row].iter_mut() {
            *entry = entry
                .checked_div(pivot)
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        let pivot_row = system[row];
        for (other, values) in system.iter_mut().enumerate() {
            if other == row || values[column].is_zero() {
                continue;
            }
            let factor = values[column];
            for (entry, value) in values.iter_mut().enumerate() {
                let term = pivot_row[entry]
                    .checked_mul(factor)
                    .map_err(|_| GeometryError::Arithmetic)?;
                *value = value
                    .checked_sub(term)
                    .map_err(|_| GeometryError::Arithmetic)?;
            }
        }
        pivots.push(column);
        row += 1;
    }
    for values in system.iter().skip(row) {
        if !values[3].is_zero() {
            return Ok(None);
        }
    }
    let mut particular = [Rat::ZERO; 3];
    for (index, column) in pivots.iter().enumerate() {
        particular[*column] = system[index][3];
    }
    let mut kernel = Vec::new();
    for column in 0..3 {
        if pivots.contains(&column) {
            continue;
        }
        let mut vector = [Rat::ZERO; 3];
        vector[column] = Rat::ONE;
        for (index, pivot) in pivots.iter().enumerate() {
            vector[*pivot] = system[index][column]
                .checked_neg()
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        kernel.push(vector);
    }
    Ok(Some((particular, kernel)))
}

// ---------------------------------------------------------------------------
// Round 23a: the exact joint congruence behind one class equation.
//
// `search_embeddings` places a class by solving the anchor equation
//
//     (I - R_p) delta  ==  tau_p - A tau        (mod L_P)                (*)
//
// for the origin shift `delta`, where `R_p` is the parent operation the class
// lands on, `L_P = Z^3 + span_Z(parent's listed pure translations)` is the
// parent's translation lattice, `A` the candidate map and `tau` the class
// translation.  Round 22 refuted candidate maps through a **per-normal** reading
// of (*), which is only a necessary condition; round 22b withdrew it and measured
// that the correct joint reading lifts the certified zero rows from 5,649 to at
// least 8,106 of 8,294.  The functions below are that joint reading.
//
// With `M = I - R_p` (integer) and `t = tau_p - A tau` (rational), (*) has a
// solution `delta in R^3` exactly when
//
//     t  in  im_Q(M) + L_P.                                            (**)
//
// Three exact facts turn (**) into a finite integer computation:
//
// * `im_Q(M) = (ker_Q(M^T))^perp`, so (**) holds iff `Y (t + lambda) = 0` for
//   some `lambda in L_P`, where the rows of `Y` are a `Z`-basis of the
//   **integer** kernel `{y in Z^3 : M^T y = 0}`.  The integer kernel is not
//   optional: it is a *saturated* sublattice of `Z^3`, and a rational basis --
//   for instance one row-reduced RREF vector -- does not generate it.  Round 22
//   measured what a test built on such a basis answers on SG 178 / UNI 112.
// * With `D` the common denominator of the parent's pure translations, `K = D
//   L_P` is an integer lattice, so every `lambda in L_P` is `B^T z / D` for a
//   unique `z in Z^3` and a `Z`-basis `B` of `K`.  Substituting turns (*) into
//   the plain integer system
//
//       A z = b,   A = Y B^T  (integer 3x3),   b = -D (Y t).            (***)
//
// * `A z = b` over `z in Z^3` is decided by [`integer_system_is_solvable`],
//   which never truncates a coefficient range and whose completeness argument is
//   written at the function.
//
// `Z^3` is implicit in `L_P` throughout, as in [`canonical_translation`] and
// [`in_lattice`]: the parent's operation table lists only the **non-zero** pure
// translations, and forgetting the cell translations makes every integer target
// unplaceable (round 22's first version had exactly that bug).
// ---------------------------------------------------------------------------

/// Dot product of two integer vectors.
fn dot_product(left: [i128; 3], right: [i128; 3]) -> i128 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

/// Cross product of two integer vectors.
fn cross_product(left: [i128; 3], right: [i128; 3]) -> [i128; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

/// `vector / gcd(vector)` for a non-zero integer vector, `None` for the zero one.
fn primitive_vector(vector: [i128; 3]) -> Option<[i128; 3]> {
    let divisor = gcd(gcd(vector[0], vector[1]), vector[2]);
    if divisor == 0 {
        return None;
    }
    Some([
        vector[0] / divisor,
        vector[1] / divisor,
        vector[2] / divisor,
    ])
}

/// `(g, s, t)` with `g = gcd(left, right) >= 0` and `s left + t right = g`.
fn extended_gcd(left: i128, right: i128) -> (i128, i128, i128) {
    let (mut old_remainder, mut remainder) = (left, right);
    let (mut old_left, mut left_coefficient) = (1i128, 0i128);
    let (mut old_right, mut right_coefficient) = (0i128, 1i128);
    while remainder != 0 {
        let quotient = old_remainder / remainder;
        let next = old_remainder - quotient * remainder;
        old_remainder = remainder;
        remainder = next;
        let next = old_left - quotient * left_coefficient;
        old_left = left_coefficient;
        left_coefficient = next;
        let next = old_right - quotient * right_coefficient;
        old_right = right_coefficient;
        right_coefficient = next;
    }
    if old_remainder < 0 {
        (-old_remainder, -old_left, -old_right)
    } else {
        (old_remainder, old_left, old_right)
    }
}

/// A `Z`-basis `{u, w}` with `u x w = +-normal` of the saturated plane
/// `{y in Z^3 : y . normal = 0}` for a **primitive** non-zero `normal`.
///
/// `u` is `(0, -c/d, b/d)` and `w` is `(d, -a s, -a t)` with `d = gcd(b, c)` and
/// Bézout coefficients `s b + t c = d`.  Both are orthogonal to `normal`, and
/// `u x w = normal`, so the pair generates the whole plane: a rank-two sublattice
/// of `Z^3` is saturated exactly when its two generators have a primitive cross
/// product.  A rational basis of the same plane does not have that property,
/// which is why round 22's RREF vector is not usable here.
fn plane_basis(normal: [i128; 3]) -> Result<[[i128; 3]; 2], GeometryError> {
    let [a, b, c] = normal;
    if b == 0 && c == 0 {
        // `normal` is primitive, so it is `+-e_1` and the plane is `y_1 = 0`.
        return Ok([[0, 1, 0], [0, 0, 1]]);
    }
    let (divisor, s, t) = extended_gcd(b, c);
    if divisor == 0 {
        return Err(GeometryError::Arithmetic);
    }
    Ok([[0, -c / divisor, b / divisor], [divisor, -a * s, -a * t]])
}

/// A `Z`-basis of the integer left kernel `{y in Z^3 : M^T y = 0}` of a 3x3
/// integer matrix -- the orthogonal complement of its column space inside `Z^3`.
///
/// The kernel of an integral matrix is saturated, so its rank determines it
/// completely in dimension three: no non-zero column means `M = 0` and the kernel
/// is all of `Z^3`; one independent column `v` means the plane `v^perp`
/// ([`plane_basis`] of the primitive `v`); two independent columns `v, w` mean
/// the line spanned by the primitive `v x w`; three independent columns mean the
/// kernel is `{0}`.
fn integer_left_kernel(matrix: &Mat3I) -> Result<Vec<[i128; 3]>, GeometryError> {
    let column = |index: usize| -> [i128; 3] {
        [
            i128::from(matrix[0][index]),
            i128::from(matrix[1][index]),
            i128::from(matrix[2][index]),
        ]
    };
    let columns = [column(0), column(1), column(2)];
    let Some(first) = columns.iter().copied().find(|value| *value != [0, 0, 0]) else {
        return Ok(vec![[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    };
    let independent = columns
        .iter()
        .copied()
        .find(|value| cross_product(first, *value) != [0, 0, 0]);
    let Some(second) = independent else {
        let generator = primitive_vector(first).ok_or(GeometryError::Arithmetic)?;
        return Ok(plane_basis(generator)?.to_vec());
    };
    let normal = primitive_vector(cross_product(first, second)).ok_or(GeometryError::Arithmetic)?;
    for value in &columns {
        if dot_product(*value, normal) != 0 {
            return Ok(Vec::new());
        }
    }
    Ok(vec![normal])
}

/// A `Z`-basis of `Z^3 ∩ span_Q(rows)` -- the saturation of the row lattice.
///
/// [`integer_system_is_solvable`] quantifies over exactly this lattice: it is the
/// annihilator of `Z^3 + ker_Q(A)` inside `R^3`, and it is a subgroup, so a basis
/// is a complete test set.  As in [`integer_left_kernel`] the rank decides the
/// answer in dimension three.
fn saturated_row_basis(rows: &[[i128; 3]]) -> Result<Vec<[i128; 3]>, GeometryError> {
    let Some(first) = rows.iter().copied().find(|row| *row != [0, 0, 0]) else {
        return Ok(Vec::new());
    };
    let independent = rows
        .iter()
        .copied()
        .find(|row| cross_product(first, *row) != [0, 0, 0]);
    let Some(second) = independent else {
        return Ok(vec![
            primitive_vector(first).ok_or(GeometryError::Arithmetic)?,
        ]);
    };
    let normal = primitive_vector(cross_product(first, second)).ok_or(GeometryError::Arithmetic)?;
    for row in rows {
        if dot_product(*row, normal) != 0 {
            return Ok(vec![[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        }
    }
    Ok(plane_basis(normal)?.to_vec())
}

/// Decide `A z = b` for `z in Z^3` exactly, with `rows <= 3` equations.
///
/// Let `W = ker_Q(A)` and let `z_0` be a rational particular solution (exact
/// Gaussian elimination; no solution over `Q` means none over `Z`).  The rational
/// solutions are `z_0 + W`, so an integer solution exists exactly when `z_0` lies
/// in the closed subgroup `Z^3 + W` of `R^3`.  That subgroup's annihilator is
/// `(Z^3)^perp ∩ W^perp = Z^3 ∩ span_Q(rows of A)`, the saturation returned by
/// [`saturated_row_basis`], and characters separate points on the compact torus
/// `R^3/(Z^3 + W)`, so
///
/// ```text
/// z_0 in Z^3 + W   iff   y . z_0 in Z for every y in that saturation.
/// ```
///
/// The set of `y` that pass is a subgroup, so a basis is enough, and any rational
/// particular solution gives the same verdict: two of them differ by an element
/// of `W`, which every `y in W^perp` annihilates.
///
/// The saturation is what makes this complete rather than merely consistent: with
/// `A = [[2, 0, 0]]` and `b = 1` the rational solution is `z_0 = (1/2, 0, 0)`,
/// and only the primitive `y = (1, 0, 0)` rejects it -- the row `(2, 0, 0)`
/// itself would accept it.
fn integer_system_is_solvable(
    a: &[[i128; 3]; 3],
    b: &[Rat; 3],
    rows: usize,
) -> Result<bool, GeometryError> {
    let mut system = [[Rat::ZERO; 4]; 3];
    for row in 0..rows {
        for column in 0..3 {
            system[row][column] = Rat::from_integer(a[row][column]);
        }
        system[row][3] = b[row];
    }
    let particular = match solve_linear_system(system) {
        Ok(Some((particular, _kernel))) => particular,
        // Inconsistent: no integer (or real) solution exists.  A proof.
        Ok(None) => return Ok(false),
        // An arithmetic limit is *not* a refutation: report it as "undecided"
        // and let the caller refuse to refute this map (round-23 review F2).
        Err(error) => return Err(error),
    };
    for y in saturated_row_basis(&a[..rows])? {
        let mut value = Rat::ZERO;
        for axis in 0..3 {
            let term = particular[axis]
                .checked_mul(Rat::from_integer(y[axis]))
                .map_err(|_| GeometryError::Arithmetic)?;
            value = value
                .checked_add(term)
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        if !value.is_integer() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// A `Z`-basis of the lattice generated by `rows`, in elimination form.
///
/// Only row swaps and `row_i -= q row_j` with `q in Z` are used, so the row
/// lattice is preserved exactly.  The elimination is the integer Euclid loop: as
/// long as the current column still has a non-zero entry at or below the pivot
/// row, the row with the smallest absolute value is swapped up and the multiples
/// of the pivot are subtracted from the rows below, which strictly lowers the
/// smallest non-zero absolute value in that column -- so the loop terminates with
/// the column clear below the pivot, and the three pivot rows are a basis.
///
/// `Err(GeometryError::Arithmetic)` when the rows span less than `Z^3`.  In
/// production they never do: [`parent_lattice_basis`] always feeds the cell
/// translations `D e_i`, so the lattice is full rank, and a caller that dropped
/// them is told instead of silently deciding on a smaller lattice.
fn hermite_basis(rows: &[[i128; 3]]) -> Result<[[i128; 3]; 3], GeometryError> {
    let mut reduced: Vec<[i128; 3]> = rows
        .iter()
        .copied()
        .filter(|row| *row != [0, 0, 0])
        .collect();
    let mut pivot_row = 0usize;
    for column in 0..3 {
        if pivot_row >= reduced.len() {
            break;
        }
        loop {
            let mut chosen: Option<usize> = None;
            for index in pivot_row..reduced.len() {
                let value = reduced[index][column];
                if value == 0 {
                    continue;
                }
                if chosen.is_none_or(|best| value.abs() < reduced[best][column].abs()) {
                    chosen = Some(index);
                }
            }
            let Some(chosen) = chosen else { break };
            reduced.swap(pivot_row, chosen);
            let pivot = reduced[pivot_row][column];
            let pivot_values = reduced[pivot_row];
            let mut changed = false;
            for row in reduced.iter_mut().skip(pivot_row + 1) {
                let value = row[column];
                if value == 0 {
                    continue;
                }
                let quotient = value / pivot;
                for (entry, pivot_entry) in row.iter_mut().zip(pivot_values.iter()) {
                    *entry -= quotient * pivot_entry;
                }
                changed = true;
            }
            if !changed {
                break;
            }
        }
        if reduced[pivot_row][column] != 0 {
            pivot_row += 1;
        }
    }
    if pivot_row != 3 {
        return Err(GeometryError::Arithmetic);
    }
    Ok([reduced[0], reduced[1], reduced[2]])
}

/// A basis of the parent's translation lattice, scaled to integers.
///
/// Returns `(D, B)` with `D` a common denominator of the listed generators and
/// the rows of `B` a `Z`-basis of `K = D L_P = D Z^3 + span_Z(D g)`, where
/// `L_P = Z^3 + span_Z(g)` is the parent's lattice.  The cell translations
/// `D e_i` are part of the generator list on purpose: the parent's operation
/// table lists only the non-zero pure translations, and the implicit `Z^3` is
/// what makes an integer target placeable at all.  Scaling by `D` keeps `K`
/// inside `Z^3`, so all later arithmetic is integer arithmetic.
fn parent_lattice_basis(generators: &[[Rat; 3]]) -> Result<(i128, [[i128; 3]; 3]), GeometryError> {
    let mut denominator = 1i128;
    for value in generators.iter().flatten() {
        denominator = denominator / gcd(denominator, value.denominator()) * value.denominator();
    }
    let mut rows: Vec<[i128; 3]> = Vec::new();
    for axis in 0..3 {
        let mut unit = [0i128; 3];
        unit[axis] = denominator;
        rows.push(unit);
    }
    for generator in generators {
        let mut scaled = [0i128; 3];
        for axis in 0..3 {
            scaled[axis] =
                generator[axis].numerator() * (denominator / generator[axis].denominator());
        }
        if !rows.contains(&scaled) {
            rows.push(scaled);
        }
    }
    Ok((denominator, hermite_basis(&rows)?))
}

/// The `(M, L_P)`-dependent part of the anchor congruence, built once per mapped
/// rotation: the saturated integer kernel of `M^T` and the scaled coefficient
/// matrix `A = Y B^T` of `(***)`.
///
/// `AnchorSolver::new` returns `None` when `M` is invertible over `Q`: the kernel
/// is trivial and the congruence is solvable for **every** target.
struct AnchorSolver {
    /// Rows `y_j` of a `Z`-basis of `{y in Z^3 : M^T y = 0}`.
    kernel: Vec<[i128; 3]>,
    /// `A[j][i] = y_j . B_i`, the integer system of `(***)`.
    scaled: [[i128; 3]; 3],
    /// `D`, the denominator that makes the parent lattice integral.
    denominator: i128,
}

impl AnchorSolver {
    fn new(
        matrix: &Mat3I,
        denominator: i128,
        basis: &[[i128; 3]; 3],
    ) -> Result<Option<Self>, GeometryError> {
        let kernel = integer_left_kernel(matrix)?;
        if kernel.is_empty() {
            return Ok(None);
        }
        let mut scaled = [[0i128; 3]; 3];
        for (index, y) in kernel.iter().enumerate() {
            for column in 0..3 {
                scaled[index][column] = dot_product(*y, basis[column]);
            }
        }
        Ok(Some(Self {
            kernel,
            scaled,
            denominator,
        }))
    }

    /// Is `M delta == target (mod L_P)` solvable?  `false` is a proof.
    fn is_feasible(&self, target: &[Rat; 3]) -> Result<bool, GeometryError> {
        let rows = self.kernel.len();
        let mut b = [Rat::ZERO; 3];
        for (index, y) in self.kernel.iter().enumerate() {
            let mut value = Rat::ZERO;
            for axis in 0..3 {
                let term = target[axis]
                    .checked_mul(Rat::from_integer(y[axis]))
                    .map_err(|_| GeometryError::Arithmetic)?;
                value = value
                    .checked_add(term)
                    .map_err(|_| GeometryError::Arithmetic)?;
            }
            b[index] = value
                .checked_mul(Rat::from_integer(-self.denominator))
                .map_err(|_| GeometryError::Arithmetic)?;
        }
        integer_system_is_solvable(&self.scaled, &b, rows)
    }
}

/// Exact decision of the anchor congruence `M delta == target (mod L_parent)`.
///
/// `parent_generators` lists the parent's non-zero pure translations; the cell
/// translations `Z^3` are implicit, exactly as in [`canonical_translation`] and
/// [`in_lattice`].  `Ok(false)` is a proof that no `delta in R^3` solves the
/// congruence, so a class whose every matching parent operation returns `false`
/// can never be placed -- that is the refutation [`search_embeddings`] runs
/// before its shift walk.  See the block comment above for the two exact readings
/// of `target in im_Q(M) + L_parent` this is built from.
fn anchor_congruence_is_feasible(
    matrix: &Mat3I,
    target: &[Rat; 3],
    parent_generators: &[[Rat; 3]],
) -> Result<bool, GeometryError> {
    let (denominator, basis) = parent_lattice_basis(parent_generators)?;
    match AnchorSolver::new(matrix, denominator, &basis)? {
        Some(solver) => solver.is_feasible(target),
        None => Ok(true),
    }
}

const IDENTITY_ROTATION: Mat3I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// The image of the set's **full** translation lattice under the record's map.
///
/// The generators are [`full_translation_lattice`]: the cell translations `Z^3`
/// plus the listed non-zero unprimed pure translations.  Round 20: the earlier
/// version mapped only the listed ones, so for a primitive group it returned an
/// empty image although the lattice is `Z^3` -- the name promised the lattice,
/// the body returned a subset of its generators.
pub fn embedded_translation_lattice(
    geometry: &MagneticGeometry,
    set: &MagneticOperationSet,
) -> Vec<[Rat; 3]> {
    let map = geometry.map_matrix();
    let mut image: Vec<[Rat; 3]> = Vec::new();
    for generator in &full_translation_lattice(set) {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: [[i32; 3]; 3]) -> Mat3R {
        let mut exact = [[Rat::ZERO; 3]; 3];
        for (row, values) in rows.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                exact[row][column] = Rat::from_integer(i128::from(*value));
            }
        }
        Mat3R::new(exact)
    }

    fn unit_generators() -> Vec<[Rat; 3]> {
        let mut generators = Vec::new();
        for axis in 0..3 {
            let mut unit = [Rat::ZERO; 3];
            unit[axis] = Rat::ONE;
            generators.push(unit);
        }
        generators
    }

    fn cell_rows() -> [[Rat; 3]; 3] {
        let mut rows = [[Rat::ZERO; 3]; 3];
        for (axis, row) in rows.iter_mut().enumerate() {
            row[axis] = Rat::ONE;
        }
        rows
    }

    /// `A = 2I` maps `Z^3` into `Z^3` (forward holds) but `A^-1 Z^3 = (1/2)Z^3`
    /// is not inside `Z^3`, so it does **not** realise the cell lattice.  A
    /// forward-only test would call it a realisation.
    #[test]
    fn the_lattice_realisation_needs_both_inclusions() {
        let generators = unit_generators();
        let rows = cell_rows();
        let identity = Mat3R::identity();
        assert!(realises_record_lattice(
            &identity,
            &generators,
            &[],
            &rows,
            &identity
        ));
        let doubled = matrix([[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
        // Sanity: the forward half alone is satisfied.
        assert!(generators.iter().all(|generator| {
            doubled
                .checked_mul_vector(&Vec3R::new(*generator))
                .map(|image| image.as_array().iter().all(|value| value.is_integer()))
                .unwrap_or(false)
        }));
        assert!(
            !realises_record_lattice(&doubled, &generators, &[], &rows, &identity),
            "A = 2I does not realise the cell lattice"
        );
    }

    /// The record's membership test uses the **row** lattice.  With
    /// `M = [[2,1,0],[0,1,0],[0,0,1]]` the point `(1,1,0)` is not in the row
    /// lattice `{(2a, a+b, c)}` (`2a = 1` has no integer solution), while the
    /// column test `M^-1 v` accepts it; the shear
    /// `A = [[1,0,0],[1,1,0],[0,0,1]]` sends `e_1` exactly there.
    ///
    /// Scope (round-20 review, F4): this pins the helper given the already
    /// transposed matrix.  The production call site that builds
    /// `record_inverse_transposed` is pinned by the census re-derivation in
    /// `the_embedding_census_on_the_stratified_subset_is_pinned`, not here.
    #[test]
    fn the_record_membership_uses_the_row_lattice() {
        let generators = unit_generators();
        let rows = {
            let mut rows = cell_rows();
            rows[0][0] = Rat::from_integer(2);
            rows[0][1] = Rat::ONE;
            rows
        };
        let inverse_transposed = Mat3R::new(rows)
            .inverse()
            .expect("non-singular")
            .transpose();
        let inverse = Mat3R::new(rows).inverse().expect("non-singular");
        let shear = matrix([[1, 0, 0], [1, 1, 0], [0, 0, 1]]);
        let image = shear
            .checked_mul_vector(&Vec3R::new([Rat::ONE, Rat::ZERO, Rat::ZERO]))
            .expect("product");
        let row_coordinates = inverse_transposed
            .checked_mul_vector(&image)
            .expect("product");
        let column_coordinates = inverse.checked_mul_vector(&image).expect("product");
        assert!(
            !row_coordinates
                .as_array()
                .iter()
                .all(|value| value.is_integer()),
            "the row convention must reject (1,1,0)"
        );
        assert!(
            column_coordinates
                .as_array()
                .iter()
                .all(|value| value.is_integer()),
            "the column convention is the wrong test and accepts it"
        );
        assert!(!realises_record_lattice(
            &shear,
            &generators,
            &[],
            &rows,
            &inverse_transposed
        ));
    }

    /// Membership in the one-dimensional subgroup of `R` generated by `values`,
    /// i.e. in `y . (Z^3 + span_Z(pure translations))` for a normal vector `y`
    /// (the cell translations contribute `y_i`).  This is the per-normal reading
    /// round 22 shipped; the hand cases use it to show what that reading answers
    /// where the joint congruence is the decision.
    fn one_dimensional_membership(value: Rat, values: &[Rat]) -> bool {
        let mut denominator = 1i128;
        for entry in values {
            denominator = denominator / gcd(denominator, entry.denominator()) * entry.denominator();
        }
        let mut divisor = 0i128;
        for entry in values {
            divisor = gcd(
                divisor,
                entry.numerator() * (denominator / entry.denominator()),
            );
        }
        divisor != 0 && value.numerator() * (denominator / value.denominator()) % divisor == 0
    }

    /// Round 23a: the exact joint congruence `M delta == target (mod L)` on
    /// hand-written cases.  Every `false` below is a proof that no `delta in R^3`
    /// solves the congruence, and each case pins one piece of the argument:
    ///
    /// * `M = 0` (the parent rotation is the identity) reads `0 == target (mod L)`,
    ///   so the target has to be a lattice vector -- and the implicit `Z^3` is
    ///   what makes an integer target one when the parent lists no pure
    ///   translation at all;
    /// * the mirror `M = diag(0, 0, 2)` has a two-dimensional kernel.  On the
    ///   C-centred D2 case of the round-23 card the per-normal test over the
    ///   kernel basis `e_1, e_2` passes on **both** normals while the joint
    ///   congruence has no solution -- the reason a per-normal test cannot
    ///   certify;
    /// * the same target becomes solvable as soon as the parent carries the
    ///   C-centring `(1/2, 1/2, 0)`: the lattice, not the normals, decides;
    /// * an invertible `M` is solvable for every target (the kernel is empty);
    /// * the saturation of the *coefficient* system matters: with the two-fold
    ///   `M = diag(0, 2, 2)` and the parent lattice `Z^3 + (0, 0, 1/2)` the scaled
    ///   row is `(2, 0, 0)`, and only its primitive `(1, 0, 0)` rejects
    ///   `target = (1/2, 0, 0)`;
    /// * the measured SG 178 / UNI 112 witness of round 22 stays solvable: the
    ///   mapped rotation there gives `M = [[2,0,0],[1,0,0],[0,0,2]]` and
    ///   `target = (-1, -1, 1/6)` with `L = Z^3`, whose integer kernel is the line
    ///   spanned by the primitive `(1, -2, 0)` -- a rational RREF vector is not a
    ///   generator of that lattice.
    #[test]
    fn the_anchor_congruence_handles_the_hand_cases() {
        let half = |value: i128| Rat::new(value, 2).expect("a half has a non-zero denominator");
        let fraction =
            |value: i128, denominator: i128| Rat::new(value, denominator).expect("non-zero");

        // (1) `M = 0`: the target must be a lattice vector, and the lattice is
        // `Z^3` plus the listed pure translations.
        let zero = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        assert!(anchor_congruence_is_feasible(&zero, &[Rat::ZERO; 3], &[]).expect("exact"));
        assert!(
            anchor_congruence_is_feasible(&zero, &[Rat::ONE, Rat::ZERO, Rat::ZERO], &[])
                .expect("exact"),
            "an integer target is placeable through the implicit Z^3"
        );
        assert!(
            !anchor_congruence_is_feasible(&zero, &[half(1), Rat::ZERO, Rat::ZERO], &[])
                .expect("exact"),
            "a half-integer target is not a cell translation"
        );
        assert!(
            anchor_congruence_is_feasible(
                &zero,
                &[half(1), Rat::ZERO, Rat::ZERO],
                &[[half(1), Rat::ZERO, Rat::ZERO]]
            )
            .expect("exact"),
            "the listed pure translation places the half-integer target"
        );

        // (2) the mirror, the C-centred D2 counterexample and its per-normal
        // shadow (the kernel basis `e_1, e_2` generates the kernel lattice, and
        // `e_1 . L = (1/2)Z`, `e_2 . L = (1/2)Z` both accept their own value).
        let mirror = [[0, 0, 0], [0, 0, 0], [0, 0, 2]];
        let centered = [[half(1), half(1), Rat::ZERO]];
        let d2 = [half(1), Rat::ZERO, Rat::ZERO];
        assert!(
            !anchor_congruence_is_feasible(&mirror, &d2, &centered).expect("exact"),
            "the joint congruence is unsolvable although both normals pass"
        );
        assert!(
            one_dimensional_membership(half(1), &[Rat::ONE, Rat::ZERO, Rat::ZERO, half(1)]),
            "e_1 . target lies in e_1 . L"
        );
        assert!(
            one_dimensional_membership(Rat::ZERO, &[Rat::ZERO, Rat::ONE, Rat::ZERO, half(1)]),
            "e_2 . target lies in e_2 . L"
        );

        // (3) the C-centring turns the other half-integer target solvable.
        let centered_target = [half(1), half(1), Rat::ZERO];
        assert!(
            !anchor_congruence_is_feasible(&mirror, &centered_target, &[]).expect("exact"),
            "without the centring the first coordinate is pinned to Z"
        );
        assert!(
            anchor_congruence_is_feasible(&mirror, &centered_target, &centered).expect("exact"),
            "the C-centring places (1/2, 1/2, 0)"
        );

        // (4) an invertible `M` is solvable for every target.
        let invertible = [[2, 0, 0], [0, 1, 0], [0, 0, 1]];
        assert!(
            anchor_congruence_is_feasible(
                &invertible,
                &[half(1), fraction(1, 3), fraction(1, 6)],
                &[]
            )
            .expect("exact")
        );
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        assert!(
            anchor_congruence_is_feasible(&identity, &[fraction(5, 7), Rat::ZERO, Rat::ZERO], &[])
                .expect("exact")
        );

        // (5) the saturation of the coefficient system: the two-fold about `e_1`
        // with the parent lattice `Z^3 + (0, 0, 1/2)` gives the non-primitive row
        // `(2, 0, 0)`, and `1/2` can never be placed in the first coordinate.
        let twofold = [[0, 0, 0], [0, 2, 0], [0, 0, 2]];
        let c_glide = [[Rat::ZERO, Rat::ZERO, half(1)]];
        assert!(
            !anchor_congruence_is_feasible(&twofold, &[half(1), Rat::ZERO, Rat::ZERO], &c_glide)
                .expect("exact"),
            "the first coordinate of M delta is always 0 and Z^3 + (0,0,1/2) is integral there"
        );
        assert!(
            anchor_congruence_is_feasible(&twofold, &[Rat::ONE, Rat::ZERO, Rat::ZERO], &c_glide)
                .expect("exact"),
            "an integer target in the first coordinate is placeable"
        );

        // (6) the measured SG 178 / UNI 112 witness: `lambda = (-1, 0, 0)` solves
        // it, with the integer kernel spanned by the primitive `(1, -2, 0)`.
        let mapped = [[-1, 0, 0], [-1, 1, 0], [0, 0, -1]];
        let mut difference = [[0i32; 3]; 3];
        for row in 0..3 {
            for column in 0..3 {
                difference[row][column] = (if row == column { 1 } else { 0 }) - mapped[row][column];
            }
        }
        assert_eq!(difference, [[2, 0, 0], [1, 0, 0], [0, 0, 2]]);
        assert!(
            anchor_congruence_is_feasible(
                &difference,
                &[Rat::from_integer(-1), Rat::from_integer(-1), fraction(1, 6)],
                &[]
            )
            .expect("exact"),
            "round 22 measured this witness as solvable"
        );

        // (7) the same saturation on the *plane* branch of the coefficient
        // system: the rows `(-1, -2, 0)` and `(1, 0, 0)` have the non-primitive
        // cross product `(0, 0, 2)`, and the rational solution has
        // `z_0 = (1/2, 0, 0, ...)` -- only the primitive normal `(0, 0, 1)` (i.e.
        // the plane `z = 0`) rejects `target = (-1/4, -1/4, 0)`.
        assert!(
            !anchor_congruence_is_feasible(
                &mirror,
                &[fraction(-1, 4), fraction(-1, 4), Rat::ZERO],
                &centered
            )
            .expect("exact"),
            "(-1/4, -1/4, 0) is not placeable modulo Z^3 + (1/2, 1/2, 0)"
        );
    }

    /// Round 23a: the exact congruence against a **different** code path.
    ///
    /// Every (matrix, target, lattice) case enumerates `delta` on the `1/24` grid
    /// and asks the production lattice membership ([`in_lattice`], i.e. the
    /// audited `canonical_translation` BFS, which shares no code with the integer
    /// system above) whether `M delta - target` is a parent lattice vector.  The
    /// two directions are asserted where they are conclusive:
    ///
    /// * a grid hit must be accepted by the decider.  A `false` from the decider
    ///   claims that **no** `delta in R^3` works, so this half is the soundness
    ///   direction (no false refutations) and it is conclusive whatever the grid
    ///   misses;
    /// * a `true` from the decider must have a grid hit.  The grid alone cannot
    ///   prove infeasibility, so this half is only as strong as the case choice:
    ///   the targets have denominators dividing `12`, the lattices dividing `6`,
    ///   and the matrices have determinant `0`, `+-1` or `2`, so every solvable
    ///   case has a solution on the `1/24` grid (halving a `1/24` component is a
    ///   `1/48` component, which is why the determinants stay at `2`).  A case
    ///   the grid cannot decide fails here instead of passing silently.
    ///
    /// The enumeration is memoised on the exact difference vector and stops at the
    /// first hit: for a singular `M` only a few dozen difference values occur, so
    /// the expensive membership call runs once per distinct value rather than once
    /// per grid point.  The D2 counterexample of the round-23 card is one of the
    /// cases and is asserted separately, together with both answer classes.
    #[test]
    fn the_anchor_congruence_agrees_with_a_brute_force_grid() {
        let grid = 24i128;
        let half = |value: i128| Rat::new(value, 2).expect("a half has a non-zero denominator");
        let fraction =
            |value: i128, denominator: i128| Rat::new(value, denominator).expect("non-zero");
        let matrices: [[[i32; 3]; 3]; 12] = [
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 2]],
            [[0, 0, 0], [0, 2, 0], [0, 0, 2]],
            [[2, 0, 0], [1, 0, 0], [0, 0, 2]],
            [[1, -1, 0], [0, 1, 0], [0, 0, 0]],
            [[1, 1, 0], [0, 1, 1], [0, 0, 1]],
            [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
            [[1, 0, 0], [0, 1, 0], [0, 0, 2]],
            [[2, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        ];
        // The first lattice is the plain cell lattice (no listed pure
        // translation); the rest cover the centrings and half translations the
        // corpus has.
        let lattices: Vec<Vec<[Rat; 3]>> = vec![
            vec![],
            vec![[half(1), half(1), Rat::ZERO]],
            vec![[half(1), Rat::ZERO, Rat::ZERO]],
            vec![[Rat::ZERO, Rat::ZERO, half(1)]],
            vec![[fraction(1, 3), fraction(1, 3), fraction(1, 3)]],
            vec![[half(1), half(1), half(1)]],
        ];
        let targets: Vec<[Rat; 3]> = vec![
            [Rat::ZERO, Rat::ZERO, Rat::ZERO],
            [Rat::ONE, Rat::ZERO, Rat::ZERO],
            [half(1), Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, half(1), Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, half(1)],
            [half(1), half(1), Rat::ZERO],
            [half(1), half(1), half(1)],
            [fraction(1, 3), fraction(1, 3), fraction(1, 3)],
            [fraction(1, 4), Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, fraction(1, 4)],
            [fraction(1, 12), Rat::ZERO, Rat::ZERO],
            [Rat::ZERO, fraction(1, 12), Rat::ZERO],
            [Rat::ZERO, Rat::ZERO, fraction(1, 12)],
            [half(1), fraction(1, 3), fraction(1, 6)],
            [fraction(1, 6), fraction(-1, 6), Rat::ZERO],
            [fraction(-1, 4), fraction(-1, 4), Rat::ZERO],
            [fraction(1, 3), Rat::ZERO, fraction(1, 3)],
            [fraction(3, 4), half(1), Rat::ZERO],
            [fraction(5, 6), fraction(1, 4), fraction(1, 12)],
        ];

        let mut cases = 0usize;
        let mut solvable = 0usize;
        let mut unsolvable = 0usize;
        for matrix in &matrices {
            for lattice in &lattices {
                for target in &targets {
                    cases += 1;
                    let decided = anchor_congruence_is_feasible(matrix, target, lattice)
                        .expect("the congruence is decided exactly");
                    // `24 * target` is integral: every target denominator divides
                    // the grid.
                    let mut base = [0i128; 3];
                    for axis in 0..3 {
                        base[axis] = target[axis].numerator() * (grid / target[axis].denominator());
                    }
                    let mut seen: std::collections::HashSet<[i128; 3]> =
                        std::collections::HashSet::new();
                    let mut witness: Option<[i128; 3]> = None;
                    'grid: for first in 0..grid {
                        for second in 0..grid {
                            for third in 0..grid {
                                let point = [first, second, third];
                                let mut scaled = [0i128; 3];
                                for row in 0..3 {
                                    let mut sum = 0i128;
                                    for column in 0..3 {
                                        sum += i128::from(matrix[row][column]) * point[column];
                                    }
                                    scaled[row] = sum - base[row];
                                }
                                if !seen.insert(scaled) {
                                    continue;
                                }
                                let mut difference = [Rat::ZERO; 3];
                                for axis in 0..3 {
                                    difference[axis] = Rat::new(scaled[axis], grid)
                                        .expect("the grid denominator is non-zero");
                                }
                                if in_lattice(difference, lattice).unwrap() {
                                    witness = Some(point);
                                    break 'grid;
                                }
                            }
                        }
                    }
                    assert!(
                        decided || witness.is_none(),
                        "the decider refuted a solvable case: M = {matrix:?}, L = {lattice:?}, \
                         target = {target:?} is solved at delta = {witness:?}/{grid}"
                    );
                    assert!(
                        !decided || witness.is_some(),
                        "the decider accepted a case the {grid}-grid does not solve: \
                         M = {matrix:?}, L = {lattice:?}, target = {target:?}"
                    );
                    if decided {
                        solvable += 1;
                    } else {
                        unsolvable += 1;
                    }
                }
            }
        }
        assert!(
            cases >= 1_000,
            "the grid oracle must cover at least 1,000 cases, got {cases}"
        );
        assert!(
            solvable > 0 && unsolvable > 0,
            "both answers must occur: {solvable} solvable, {unsolvable} unsolvable"
        );

        // The D2 counterexample of the round-23 card, asserted on its own.
        let mirror = [[0, 0, 0], [0, 0, 0], [0, 0, 2]];
        let centered = vec![[half(1), half(1), Rat::ZERO]];
        let d2 = [half(1), Rat::ZERO, Rat::ZERO];
        assert!(!anchor_congruence_is_feasible(&mirror, &d2, &centered).expect("exact"));
        assert!(
            one_dimensional_membership(half(1), &[Rat::ONE, Rat::ZERO, Rat::ZERO, half(1)]),
            "the first normal passes the per-normal test"
        );
        assert!(
            one_dimensional_membership(Rat::ZERO, &[Rat::ZERO, Rat::ONE, Rat::ZERO, half(1)]),
            "the second normal passes the per-normal test"
        );
        println!(
            "anchor congruence grid: {cases} cases, {solvable} solvable, {unsolvable} unsolvable"
        );
    }
}
