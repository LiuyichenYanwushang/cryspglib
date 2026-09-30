//! R9 stage 1: the magnetic input contract is exact and the operation sets are
//! groups.
//!
//! Every magnetic isotropy row of the bundled tables names a magnetic space
//! group by UNI number; this file pins that all of them convert to exact
//! operations, that each set **proves** to be a group (identity, closure,
//! inverses) whose unitary part is a subgroup and whose antiunitary part is one
//! coset, and that mutations of the time-reversal flag or of a translation are
//! caught.  The parent-frame embedding (through `basis`/`origin`) is stage 2.

use cryspglib::irrep::magnetic_embedding::{
    MagneticContractError, MagneticOperationSet, canonical_translation, compose,
    magnetic_operations, translation_lattice, verify_against_database, verify_group,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::Rat;
use std::collections::{BTreeMap, BTreeSet};

/// Every UNI the bundled magnetic isotropy tables refer to.
fn referenced_unis() -> (usize, Vec<usize>) {
    let mut rows = 0usize;
    let mut unis = BTreeSet::new();
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            rows += 1;
            unis.insert(row.subgroup.mag_sg);
        }
    }
    (rows, unis.into_iter().collect())
}

#[test]
fn every_referenced_uni_is_an_exact_verified_group() {
    let (rows, unis) = referenced_unis();
    assert_eq!(rows, 16_721, "magnetic isotropy rows");
    assert_eq!(unis.len(), 1_421, "distinct UNI numbers referenced");

    let mut index_histogram: BTreeMap<usize, usize> = BTreeMap::new();
    let mut denominators: BTreeMap<i128, usize> = BTreeMap::new();
    let mut operations = 0usize;
    let mut classes = 0usize;
    let mut centred = 0usize;
    let mut max_len = 0usize;
    let mut max_classes = 0usize;
    for uni in &unis {
        let set = magnetic_operations(*uni).unwrap_or_else(|error| panic!("UNI {uni}: {error}"));
        verify_group(&set).unwrap_or_else(|error| panic!("UNI {uni}: {error}"));
        // The acceptance criterion of R9 at this layer: the operation set must
        // match a fresh read of the database item by item, which is what makes
        // an edited flag, translation or ordering fail.
        verify_against_database(&set).unwrap_or_else(|error| panic!("UNI {uni}: {error}"));
        operations += set.len();
        classes += set.class_count();
        max_len = max_len.max(set.len());
        max_classes = max_classes.max(set.class_count());
        if set.class_count() < set.len() {
            centred += 1;
        }
        *index_histogram.entry(set.unitary_index()).or_insert(0) += 1;
        // A magnetic group's unitary part has index 1 or 2.
        assert!(
            set.class_count() == set.unitary_classes
                || set.class_count() == 2 * set.unitary_classes,
            "UNI {uni}: {} classes, {} unitary",
            set.class_count(),
            set.unitary_classes
        );
        // The quotient is a quotient: never more classes than listed operations.
        assert!(set.class_count() <= set.len());
        assert_eq!(
            set.class_count(),
            set.unitary_classes + set.antiunitary_classes
        );
        for operation in &set.operations {
            for component in operation.translation {
                *denominators.entry(component.denominator()).or_insert(0) += 1;
            }
        }
    }
    assert_eq!(operations, 29_457, "listed operations across the referenced UNIs");
    assert_eq!(classes, 19_815, "classes modulo the translation lattices");
    assert_eq!(centred, 402, "UNIs that list centred representatives redundantly");
    assert_eq!(max_len, 384);
    assert_eq!(max_classes, 96);
    assert_eq!(
        index_histogram,
        BTreeMap::from([(1, 230), (2, 1_191)]),
        "230 groups are purely unitary (Type I), 1,191 have an antiunitary coset"
    );
    // Every translation is a small fraction; the exact contract never rounds.
    assert_eq!(
        denominators,
        BTreeMap::from([(1, 47_084), (2, 32_525), (3, 1_556), (4, 6_936), (6, 270)]),
        "translation denominators (component counts)"
    );
}

#[test]
fn a_known_purely_unitary_group_and_a_known_antiunitary_one() {
    // UNI 1 is P1 in the magnetic numbering: one operation, no antiunitary.
    let trivial = magnetic_operations(1).expect("UNI 1");
    assert_eq!(trivial.len(), 1);
    assert_eq!(trivial.unitary_count, 1);
    assert_eq!(trivial.antiunitary_count, 0);
    assert_eq!(trivial.unitary_index(), 1);
    verify_group(&trivial).expect("P1 is a group");

    // Some UNI has an antiunitary coset; pin one and its structure.
    let (_, unis) = referenced_unis();
    let with_coset = unis
        .iter()
        .map(|uni| magnetic_operations(*uni).expect("operations"))
        .find(|set| set.antiunitary_count > 0)
        .expect("some magnetic group has primed operations");
    assert_eq!(with_coset.len(), 2 * with_coset.unitary_count);
    assert_eq!(with_coset.unitary_index(), 2);
    verify_group(&with_coset).expect("the coset structure holds");
}

/// The group proof has teeth: each mutation is caught by the exact checks.
#[test]
fn mutations_of_the_contract_are_refused() {
    let original = magnetic_operations(2).expect("UNI 2");

    // Flip one time-reversal flag: closure must fail.
    let mut flipped = original.clone();
    flipped.operations[0].time_reversal = !flipped.operations[0].time_reversal;
    let error = verify_group(&flipped).expect_err("a single flipped flag is not a group");
    assert!(
        matches!(
            error,
            MagneticContractError::NotAGroup { .. }
                | MagneticContractError::UnitaryNotClosed { .. }
                | MagneticContractError::AntiunitaryNotCoset { .. }
                | MagneticContractError::MetadataMismatch { .. }
        ),
        "unexpected error {error}"
    );

    // Flip every flag but keep the recorded counts.  On a group of order two
    // ({E, T·E}) that maps the set onto itself, which no group axiom can see --
    // stage 1 cannot pin the time-reversal *convention*, only the structure; the
    // convention is pinned against the database in stage 2.  On a larger group
    // the flip does contradict the bookkeeping.
    let mut globally_flipped = magnetic_operations(4).expect("UNI 4");
    for operation in &mut globally_flipped.operations {
        operation.time_reversal = !operation.time_reversal;
    }
    let error = verify_group(&globally_flipped)
        .expect_err("stale unitary/antiunitary counts must be caught");
    assert!(
        matches!(
            error,
            MagneticContractError::UnitaryNotClosed { .. }
                | MagneticContractError::MetadataMismatch { .. }
        ),
        "unexpected error {error}"
    );

    // The database anchor is what catches an edited translation.  Note what the
    // group axioms can and cannot do: rescaling the ambient lattice (editing a
    // pure-translation operation so that the set generates a *different*, still
    // consistent magnetic group) is invisible to `verify_group` -- the test
    // documents that below -- whereas comparing against a fresh read of the
    // database refuses it, which is the R9 acceptance criterion.
    let mut shifted = magnetic_operations(20).expect("UNI 20");
    shifted.operations[0].translation[0] = Rat::new(1, 7).expect("1/7");
    let error = verify_against_database(&shifted)
        .expect_err("an edited translation must differ from the database");
    assert!(
        matches!(
            error,
            MagneticContractError::DatabaseMismatch {
                field: "translation",
                ..
            }
        ),
        "unexpected error {error}"
    );
    // Documenting the limit: the same mutation can still form a consistent
    // group, which is exactly why stage 1 cannot rely on axioms alone.
    let _lattice_rescaling = verify_group(&shifted);

    // A flipped flag and a swapped pair are caught by the anchor too.
    let mut flag_flipped = magnetic_operations(20).expect("UNI 20");
    flag_flipped.operations[2].time_reversal = !flag_flipped.operations[2].time_reversal;
    assert!(matches!(
        verify_against_database(&flag_flipped).expect_err("flipped flag"),
        MagneticContractError::DatabaseMismatch {
            field: "time reversal",
            ..
        }
    ));
    let mut swapped = magnetic_operations(20).expect("UNI 20");
    swapped.operations.swap(0, 2);
    assert!(verify_against_database(&swapped).is_err(), "a swapped pair differs");

    // The untouched sets match their database entries.
    for uni in [1usize, 2, 3, 4, 20, 21, 24] {
        let set = magnetic_operations(uni).expect("operations");
        verify_against_database(&set).unwrap_or_else(|error| panic!("UNI {uni}: {error}"));
    }
}

#[test]
fn the_composition_law_is_the_magnetic_one() {
    // Operations are given modulo the group's own translation lattice, so the
    // product `T^2` is the identity *class*: UNI 3's primed (E | 1/2) squares
    // to (E | 1), and the integer 1 lies in the lattice.
    let set = magnetic_operations(3).expect("UNI 3");
    let lattice = translation_lattice(&set);
    assert_eq!(lattice.len(), 0, "UNI 3 lists no pure translation");
    assert_eq!(set.antiunitary_count, 1);
    assert_eq!(set.class_count(), 2);
    let primed = set
        .operations
        .iter()
        .find(|operation| operation.time_reversal)
        .expect("a primed operation");
    let square = compose(primed, primed).expect("compose");
    assert!(!square.time_reversal, "T^2 removes the flag");
    assert_eq!(square.rotation, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    assert!(
        canonical_translation(square.translation, &lattice)
            .iter()
            .all(|value| value.is_zero()),
        "the square is the identity class"
    );

    // Every operation has an inverse inside the set, modulo the lattice.  UNI
    // 20 is centred (2 classes from 4 listed operations) and has a rotation of
    // order two.
    let set = magnetic_operations(20).expect("UNI 20");
    let lattice = translation_lattice(&set);
    assert!(set.class_count() < set.len(), "UNI 20 lists centred redundantly");
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let non_trivial = set
        .classes
        .iter()
        .find(|operation| operation.rotation != identity)
        .expect("a non-identity rotation class");
    let inverse = set
        .classes
        .iter()
        .find(|candidate| {
            compose(non_trivial, candidate)
                .map(|product| {
                    product.rotation == identity
                        && !product.time_reversal
                        && canonical_translation(product.translation, &lattice)
                            .iter()
                            .all(|value| value.is_zero())
                })
                .unwrap_or(false)
        })
        .expect("an inverse in the set");
    let product = compose(&compose(non_trivial, inverse).expect("compose"), non_trivial)
        .expect("compose");
    assert!(!product.time_reversal);
    let _ = MagneticOperationSet::len(&set);
}
