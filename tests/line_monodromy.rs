//! Monodromy of the frozen line sources under a reciprocal-lattice shift.
//!
//! The contract, its derivation and the fingerprint mechanism are documented in
//! [`cryspglib::irrep::line_monodromy`].  The tests here pin four separate
//! things:
//!
//! 1. the algebraic identities of the map itself (`M_0 = id`, `M_-K = M_K^-1`,
//!    `M_{K1+K2} = M_K2 . M_K1`, `C . M_K = M_-K . C`) on every line parent;
//! 2. the **pinned** table obeys the physical consequence of the contract -- the
//!    trivial content is gauge invariant, so it is constant on a monodromy orbit
//!    (`freq(alpha) == freq(M_v(alpha))`); this is a statement about the frozen
//!    data, not about the engine;
//! 3. the engine transports a parameter step by that map: the decomposition at
//!    `t + 1` for label `alpha` is the decomposition at `t` for
//!    `M_v(alpha)` -- *not* the decomposition of `alpha` again, which is the
//!    over-strong `M == 1` reading R6.1 shipped and R6.2 revoked.
//! 4. reciprocal shifts are validated against the exact parent lattice and each
//!    source little group; primitive along-line shifts transport every pinned row.
use cryspglib::irrep::line_monodromy::{
    LabelImage, LineMonodromy, complex_conjugation, line_direction, line_parents, line_sources,
    line_table, monodromy,
};
use cryspglib::irrep::subduction::star::decompose::{
    LineSubduction, official_line_parameter, subduce_line_at_parameter,
};
use cryspglib::irrep::subduction::{
    ExactSeitz, Lattice, Mat3R, Rat, SubductionError, SubgroupEmbedding, Vec3R,
};
use cryspglib::irrep::{LabelConvention, isotropy, query};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};

/// One step of a source's own line: `k(t + 1) = k(t) + v`.
fn one_step(table: &cryspglib::irrep::w_little_characters_data::LittleCharacterTable) -> Vec3R {
    line_direction(table).expect("the frozen direction parses")
}

/// A parameter value `official + steps`.
fn stepped(steps: i128) -> Rat {
    let official = official_line_parameter().expect("the official parameter is valid");
    Rat::new(
        official.numerator() + steps * official.denominator(),
        official.denominator(),
    )
    .expect("the stepped parameter is valid")
}

/// The uniquely determined part of a map.
fn unique_part(map: &LineMonodromy) -> BTreeMap<&'static str, &'static str> {
    map.images()
        .filter_map(|(from, to)| to.unique().map(|target| (from, target)))
        .collect()
}

fn checked_monodromy(parent: u8, shift: &Vec3R) -> LineMonodromy {
    monodromy(parent, shift).unwrap_or_else(|error| {
        panic!("SG {parent}: expected reciprocal-lattice shift {shift}: {error}")
    })
}

fn parent_direct_lattice(parent: u8) -> Lattice {
    let primitive = isotropy::parent_primitive_basis(parent).expect("valid parent SG");
    let direct = Mat3R::new(primitive.map(|row| {
        row.map(|value| Rat::from_grid(value, 12).expect("ITA basis is on the 1/12 grid"))
    }));
    Lattice::new(direct).expect("primitive basis is nonsingular")
}

fn parent_reciprocal_lattice(parent: u8) -> Lattice {
    parent_direct_lattice(parent)
        .reciprocal()
        .expect("reciprocal basis is nonsingular")
}

fn gcd(mut left: i128, mut right: i128) -> i128 {
    left = left.abs();
    right = right.abs();
    while right != 0 {
        (left, right) = (right, left % right);
    }
    left
}

fn primitive_reciprocal_step(
    table: &cryspglib::irrep::w_little_characters_data::LittleCharacterTable,
) -> (Vec3R, usize) {
    let direction = line_direction(table).expect("frozen direction parses");
    let reciprocal = parent_reciprocal_lattice(table.space_group);
    let coordinates = reciprocal
        .coordinates(&direction)
        .expect("direction has reciprocal coordinates");
    assert!(
        coordinates
            .as_array()
            .iter()
            .all(|value| value.is_integer())
    );
    let divisor = coordinates
        .as_array()
        .iter()
        .fold(0, |divisor, value| gcd(divisor, value.numerator()));
    assert!(divisor > 0, "frozen line direction is nonzero");
    let scale = Rat::new(1, divisor).expect("positive reciprocal step");
    let primitive = Vec3R::new(
        direction
            .as_array()
            .map(|value| value.checked_mul(scale).expect("step component")),
    );
    (
        primitive,
        usize::try_from(divisor).expect("step index fits usize"),
    )
}

fn twist_is_character(
    table: &cryspglib::irrep::w_little_characters_data::LittleCharacterTable,
    shift: &Vec3R,
) -> bool {
    table.operations.iter().all(|operation| {
        let rotation = operation.rotation.map(|row| row.map(i32::from));
        let residual = Mat3R::from_ints(rotation)
            .transpose()
            .checked_mul_vector(shift)
            .and_then(|image| shift.checked_sub(&image));
        residual.is_ok_and(|residual| {
            table.operations.iter().all(|other| {
                let Some(translation) =
                    cryspglib::irrep::line_monodromy::operation_translation(other)
                else {
                    return false;
                };
                let mut pairing = Rat::ZERO;
                for axis in 0..3 {
                    let Ok(term) = residual.get(axis).checked_mul(translation.get(axis)) else {
                        return false;
                    };
                    let Ok(sum) = pairing.checked_add(term) else {
                        return false;
                    };
                    pairing = sum;
                }
                pairing.is_integer()
            })
        })
    })
}

fn exact_operation(
    operation: &cryspglib::irrep::w_little_characters_data::LittleOperation,
) -> ExactSeitz {
    ExactSeitz::new(
        operation.rotation.map(|row| row.map(i32::from)),
        cryspglib::irrep::line_monodromy::operation_translation(operation)
            .expect("frozen translation parses"),
    )
}

fn dot(left: &Vec3R, right: &Vec3R) -> Rat {
    let mut product = Rat::ZERO;
    for axis in 0..3 {
        product = product
            .checked_add(
                left.get(axis)
                    .checked_mul(right.get(axis))
                    .expect("small frozen rational product"),
            )
            .expect("small frozen rational sum");
    }
    product
}

/// Independently verify the twist by multiplying the exact affine operations,
/// reducing their product modulo the parent lattice, and checking phase
/// multiplicativity against that representative.
fn twist_matches_seitz_products(
    table: &cryspglib::irrep::w_little_characters_data::LittleCharacterTable,
    shift: &Vec3R,
    parent_lattice: &Lattice,
) -> bool {
    for left in table.operations {
        for right in table.operations {
            let product = exact_operation(left)
                .compose(&exact_operation(right))
                .expect("frozen Seitz product is exact");
            let representative = table.operations.iter().find(|candidate| {
                let candidate = exact_operation(candidate);
                candidate.rotation() == product.rotation()
                    && product
                        .translation()
                        .checked_sub(candidate.translation())
                        .and_then(|difference| parent_lattice.contains(&difference))
                        .is_ok_and(|is_lattice| is_lattice)
            });
            let Some(representative) = representative else {
                panic!(
                    "SG {} {}: frozen little group is not closed modulo its parent lattice",
                    table.space_group, table.label
                );
            };
            let multiplicativity_residual = dot(shift, exact_operation(left).translation())
                .checked_add(dot(shift, exact_operation(right).translation()))
                .and_then(|value| {
                    value.checked_sub(dot(shift, exact_operation(representative).translation()))
                })
                .expect("frozen phase sum is exact");
            if !multiplicativity_residual.is_integer() {
                return false;
            }
        }
    }
    true
}

fn shift_is_fixed_by_reciprocal_action(
    table: &cryspglib::irrep::w_little_characters_data::LittleCharacterTable,
    shift: &Vec3R,
) -> bool {
    table.operations.iter().all(|operation| {
        let rotation = operation.rotation.map(|row| row.map(i32::from));
        Mat3R::from_ints(rotation)
            .inverse()
            .and_then(|inverse| inverse.transpose().checked_mul_vector(shift))
            .is_ok_and(|image| image == *shift)
    })
}

fn inverse(map: &BTreeMap<&'static str, &'static str>) -> BTreeMap<&'static str, &'static str> {
    let mut out = BTreeMap::new();
    for (from, to) in map {
        assert!(
            out.insert(*to, *from).is_none(),
            "monodromy map is not injective at {to}"
        );
    }
    out
}

fn compose(
    first: &BTreeMap<&'static str, &'static str>,
    second: &BTreeMap<&'static str, &'static str>,
) -> BTreeMap<&'static str, &'static str> {
    first
        .iter()
        .map(|(from, middle)| (*from, *second.get(middle).expect("composable")))
        .collect()
}

/// The four algebraic identities for every distinct frozen direction of each
/// parent. Each shift is applied to every source of that parent whose twist is
/// a character, so this also checks valid cross-line twists.
///
/// The identities are asserted on the **uniquely determined** part of each map;
/// the frozen little-group tables do not always determine it (two sources can
/// share a fingerprint), and those labels are reported so the boundary is
/// visible instead of being guessed.
#[test]
fn the_monodromy_contract_identities_hold() {
    let zero = Vec3R::zero();
    let mut non_trivial: BTreeSet<(&'static str, &'static str)> = BTreeSet::new();
    let mut identity_undetermined = 0usize;
    let mut identity_checked = 0usize;
    let mut inverse_checked = 0usize;
    let mut doubled_checked = 0usize;
    let mut conjugation_checked = 0usize;
    for parent in line_parents() {
        let mut directions: Vec<Vec3R> = Vec::new();
        for source in line_sources(parent) {
            let v = one_step(source);
            if !directions.contains(&v) {
                directions.push(v);
            }
        }
        let identity = checked_monodromy(parent, &zero);
        for (label, target) in unique_part(&identity) {
            assert_eq!(
                target, label,
                "SG {parent}: M_0 must be the identity ({label})"
            );
            identity_checked += 1;
        }
        identity_undetermined += identity.undetermined().len();
        for v in directions {
            let forward = checked_monodromy(parent, &v);
            let backward = checked_monodromy(parent, &v.checked_neg().expect("negated direction"));
            let doubled = checked_monodromy(parent, &v.checked_add(&v).expect("doubled direction"));
            let conjugate = complex_conjugation(parent);
            let (forward, backward, doubled) = (
                unique_part(&forward),
                unique_part(&backward),
                unique_part(&doubled),
            );
            let conjugate = unique_part(&conjugate);
            for (label, target) in &forward {
                if *target != *label {
                    non_trivial.insert((*label, *target));
                }
            }
            // M_-K = M_K^-1, on the labels both directions determine.
            let forward_inverse = inverse(&forward);
            for (label, target) in &backward {
                if let Some(expected) = forward_inverse.get(label) {
                    assert_eq!(
                        target, expected,
                        "SG {parent}: M_-K must invert M_K at {label}"
                    );
                    inverse_checked += 1;
                }
            }
            // M_{2K} = M_K . M_K
            let composed = compose(&forward, &forward);
            for (label, target) in &doubled {
                if let Some(expected) = composed.get(label) {
                    assert_eq!(
                        target, expected,
                        "SG {parent}: M_{{K+K}} = M_K . M_K at {label}"
                    );
                    doubled_checked += 1;
                }
            }
            // C . M_K = M_-K . C
            for (label, image) in &forward {
                let (Some(left), Some(middle)) =
                    (conjugate.get(image).copied(), conjugate.get(label).copied())
                else {
                    continue;
                };
                if let Some(expected) = backward.get(middle) {
                    assert_eq!(
                        left, *expected,
                        "SG {parent}: conjugation must commute with the monodromy as \
                         C . M_K = M_-K . C (label {label})"
                    );
                    conjugation_checked += 1;
                }
            }
        }
    }
    assert_eq!(identity_checked, 73, "M_0 covers all frozen sources");
    assert_eq!(identity_undetermined, 0, "M_0 is uniquely determined");
    assert_eq!(
        (inverse_checked, doubled_checked, conjugation_checked),
        (136, 136, 136)
    );
    assert!(
        !non_trivial.is_empty(),
        "no parent has a non-trivial monodromy: the corpus would not exercise the contract"
    );
    println!(
        "monodromy contract: identity_checks={identity_checked} inverse_checks={inverse_checked} \
         composition_checks={doubled_checked} conjugation_checks={conjugation_checked} \
         identity_undetermined={identity_undetermined}"
    );
    println!("non-trivial monodromy moves: {non_trivial:?}");
}

/// Check all three parent reciprocal-basis generators against every source's
/// little group. These generic twists need not be steps along each source's
/// line; compare the exact criterion with Seitz products and check composition.
#[test]
fn reciprocal_basis_shifts_are_checked_per_little_group() {
    let mut supported = 0usize;
    let mut unsupported = 0usize;
    let mut missing = 0usize;
    let mut unique_images = 0usize;
    let mut ambiguous_images = 0usize;
    let mut compositions = 0usize;

    for parent in line_parents() {
        let sources = line_sources(parent);
        let reciprocal = parent_reciprocal_lattice(parent);
        let direct = parent_direct_lattice(parent);
        let shifts: Vec<_> = (0..3)
            .map(|axis| Vec3R::new(*reciprocal.rows().row(axis)))
            .collect();
        let maps: Vec<_> = shifts
            .iter()
            .map(|shift| checked_monodromy(parent, shift))
            .collect();

        for (axis, map) in maps.iter().enumerate() {
            for source in &sources {
                let preserves = twist_is_character(source, &shifts[axis]);
                assert_eq!(
                    preserves,
                    twist_matches_seitz_products(source, &shifts[axis], &direct),
                    "SG {parent} {} reciprocal generator {axis}: Seitz product gate",
                    source.label
                );
                let image = map.image(source.label).expect("map covers every source");
                assert_eq!(
                    matches!(image, LabelImage::UnsupportedShift),
                    !preserves,
                    "SG {parent} {} reciprocal generator {axis}: applicability",
                    source.label
                );
                if !preserves {
                    unsupported += 1;
                    continue;
                }
                supported += 1;
                match image {
                    LabelImage::Unique(_) => unique_images += 1,
                    LabelImage::Ambiguous(_) => ambiguous_images += 1,
                    LabelImage::Missing => missing += 1,
                    LabelImage::UnsupportedShift => unreachable!("checked above"),
                }
            }
        }

        for left in 0..3 {
            for right in left + 1..3 {
                let sum = shifts[left]
                    .checked_add(&shifts[right])
                    .expect("reciprocal basis sum");
                let combined = checked_monodromy(parent, &sum);
                for source in &sources {
                    if !twist_is_character(source, &shifts[left])
                        || !twist_is_character(source, &shifts[right])
                    {
                        continue;
                    }
                    if let Some(middle) = maps[left].unique_image(source.label)
                        && let Some(target) = maps[right].unique_image(middle)
                    {
                        assert_eq!(
                            combined.unique_image(source.label),
                            Some(target),
                            "SG {parent} {} mixed reciprocal shift",
                            source.label
                        );
                        compositions += 1;
                    }
                }
            }
        }
    }

    assert_eq!(
        (
            supported,
            unsupported,
            unique_images,
            ambiguous_images,
            missing,
            compositions
        ),
        (159, 60, 159, 0, 0, 159)
    );
    println!(
        "reciprocal-shift domain: supported={supported} unsupported={unsupported} \
         source_missing={missing} composition_images={compositions}"
    );
}

#[test]
fn exact_phase_matching_is_stable_for_large_reciprocal_shifts() {
    // The extra (4N,0,0) contributes an integer phase on every SG 203 frozen
    // translation. This also catches implementations that feed huge angles to
    // floating-point sin/cos before reducing the rational phase modulo one.
    let base = checked_monodromy(203, &Vec3R::from_ints([2, 0, 0]));
    let large = checked_monodromy(
        203,
        &Vec3R::new([Rat::from_integer(4_000_000_000_002), Rat::ZERO, Rat::ZERO]),
    );
    for (label, image) in base.images() {
        assert_eq!(large.image(label), Some(image), "label {label}");
    }
    assert_eq!(large.unique_image("DT1"), Some("DT3"));
    assert_eq!(large.unique_image("SM1"), Some("SM2"));
}

#[test]
fn non_reciprocal_shifts_are_rejected_exactly() {
    // SG 5 is C-centred, so (1,0,0) is not in its reciprocal lattice.
    assert!(matches!(
        monodromy(5, &Vec3R::from_ints([1, 0, 0])),
        Err(SubductionError::NonReciprocalShift { sg: 5, .. })
    ));
    assert!(monodromy(5, &Vec3R::from_ints([1, 1, 0])).is_ok());
    // SG 203 is F-centred: all-even (2,0,0) is reciprocal, mixed parity is not.
    assert!(matches!(
        monodromy(203, &Vec3R::from_ints([1, 0, 0])),
        Err(SubductionError::NonReciprocalShift { sg: 203, .. })
    ));
    assert!(monodromy(203, &Vec3R::from_ints([2, 0, 0])).is_ok());
}

#[test]
fn reciprocal_shift_is_checked_against_each_source_little_group() {
    // This is a cross-line twist, not a step along either source's frozen line.
    // (2,0,0) is reciprocal for F-centred SG 203. It preserves the SM mirror
    // little group and swaps SM1/SM2; it also twists DT1 to DT3.
    let map = checked_monodromy(203, &Vec3R::from_ints([2, 0, 0]));
    assert_eq!(map.unique_image("SM1"), Some("SM2"));
    let dt = line_table(203, "DT1").expect("SG 203 DT1 source");
    let shift = Vec3R::from_ints([2, 0, 0]);
    assert!(!shift_is_fixed_by_reciprocal_action(dt, &shift));
    assert!(twist_is_character(dt, &shift));
    assert_eq!(map.unique_image("DT1"), Some("DT3"));
}

#[test]
fn along_line_images_match_the_frozen_parent_labels() {
    const MOVES: &[(u8, &str, &str)] = &[
        (203, "DT1", "DT2"),
        (203, "DT2", "DT1"),
        (203, "DT3", "DT4"),
        (203, "DT4", "DT3"),
        (210, "DT1", "DT2"),
        (210, "DT2", "DT1"),
        (210, "DT3", "DT4"),
        (210, "DT4", "DT3"),
        (227, "DT1", "DT3"),
        (227, "DT3", "DT1"),
        (227, "DT2", "DT4"),
        (227, "DT4", "DT2"),
        (228, "DT1", "DT3"),
        (228, "DT3", "DT1"),
        (228, "DT2", "DT4"),
        (228, "DT4", "DT2"),
    ];

    for parent in line_parents() {
        for source in line_sources(parent) {
            let expected = MOVES
                .iter()
                .find(|(sg, label, _)| *sg == parent && *label == source.label)
                .map_or(source.label, |(_, _, image)| *image);
            let map = checked_monodromy(parent, &one_step(source));
            assert_eq!(
                map.unique_image(source.label),
                Some(expected),
                "SG {parent} {} along its own line",
                source.label
            );
        }
    }
}

#[test]
fn primitive_reciprocal_steps_cover_every_frozen_line_shift() {
    let mut unique = 0usize;
    let mut ambiguous = 0usize;
    let mut missing = 0usize;
    let mut unsupported = 0usize;
    let mut parameter_step_matches = 0usize;
    let mut step_orders = BTreeMap::new();

    for parent in line_parents() {
        for source in line_sources(parent) {
            let (step, parameter_steps) = primitive_reciprocal_step(source);
            *step_orders.entry(parameter_steps).or_insert(0usize) += 1;
            let map = checked_monodromy(parent, &step);
            match map.image(source.label).expect("source has a map entry") {
                LabelImage::Unique(_) => unique += 1,
                LabelImage::Ambiguous(_) => ambiguous += 1,
                LabelImage::Missing => missing += 1,
                LabelImage::UnsupportedShift => unsupported += 1,
            }

            let direction = line_direction(source).expect("frozen direction parses");
            let full_step = checked_monodromy(parent, &direction);
            let orbit = map.orbit(source.label, parameter_steps);
            if let (Some(orbit), Some(expected)) = (orbit, full_step.unique_image(source.label)) {
                assert_eq!(
                    orbit[parameter_steps], expected,
                    "SG {parent} {}: the primitive reciprocal step raised to {parameter_steps} \
                     must equal the frozen parameter step",
                    source.label
                );
                parameter_step_matches += 1;
            }
        }
    }
    println!(
        "primitive line reciprocal steps: unique={unique} ambiguous={ambiguous} \
         missing={missing} unsupported={unsupported} parameter_step_matches={parameter_step_matches}"
    );
    assert_eq!((unique, ambiguous, missing, unsupported), (73, 0, 0, 0));
    assert_eq!(parameter_step_matches, 73);
    assert_eq!(step_orders, BTreeMap::from([(1, 73)]));
    println!("primitive line step orders: {step_orders:?}");
    assert_eq!(
        unsupported, 0,
        "a shift along the line preserves its little group"
    );
}

/// Every isotropy record of one parent space group.
fn subgroups_of(sg: u8) -> Vec<isotropy::IsotropySubgroup> {
    let mut out = Vec::new();
    for record in query::irreps_of(sg) {
        if record.spinor || record.subgroups().is_empty() {
            continue;
        }
        let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
            .unwrap_or_else(|error| panic!("SG {sg} {}: {error}", record.ml));
        out.extend(subgroups);
    }
    out
}

/// Every pinned parametric-k row once: `(subgroup, row)`.
fn pinned_rows() -> Vec<(
    isotropy::IsotropySubgroup,
    isotropy::OtherWaveVectorSubduction,
)> {
    let mut out = Vec::new();
    for parent in line_parents() {
        for subgroup in subgroups_of(parent) {
            let Ok(rows) = subgroup.other_wave_vector_subduction() else {
                continue;
            };
            for row in rows {
                out.push((subgroup, row));
            }
        }
    }
    out
}

#[derive(Clone, Copy)]
struct FrozenLineTransport {
    primitive_parameter_steps: usize,
    primitive_image:
        Option<&'static cryspglib::irrep::w_little_characters_data::LittleCharacterTable>,
    full_step_image:
        Option<&'static cryspglib::irrep::w_little_characters_data::LittleCharacterTable>,
}

/// Precompute source-level shifts and label images once, rather than repeating
/// them for every isotropy row carrying the same parent line source.
fn frozen_line_transports() -> BTreeMap<(u8, &'static str), FrozenLineTransport> {
    let mut transports = BTreeMap::new();
    for parent in line_parents() {
        for source in line_sources(parent) {
            let (primitive_shift, primitive_parameter_steps) = primitive_reciprocal_step(source);
            let primitive_map = checked_monodromy(parent, &primitive_shift);
            let primitive_image = primitive_map
                .unique_image(source.label)
                .and_then(|label| line_table(parent, label));
            let full_map = checked_monodromy(parent, &one_step(source));
            let full_step_image = full_map
                .unique_image(source.label)
                .and_then(|label| line_table(parent, label));
            assert!(
                transports
                    .insert(
                        (parent, source.label),
                        FrozenLineTransport {
                            primitive_parameter_steps,
                            primitive_image,
                            full_step_image,
                        },
                    )
                    .is_none(),
                "SG {parent} {} appears twice in the frozen source list",
                source.label
            );
        }
    }
    transports
}

/// The trivial content of the pinned table, as a map from row label to the set
/// of frequencies the pinned rows carry for this subgroup.
fn pinned_frequencies(subgroup: &isotropy::IsotropySubgroup) -> BTreeMap<&'static str, u16> {
    let mut out = BTreeMap::new();
    if let Ok(rows) = subgroup.other_wave_vector_subduction() {
        for row in rows {
            out.insert(row.parent_ml, row.frequency);
        }
    }
    out
}

/// The physical half of the contract, checked on the **pinned** table: the
/// trivial content is a gauge-invariant number, so a label and its monodromy
/// image must carry the same pinned frequency.
///
/// This is the statement that makes `M != 1` compatible with the pinned data; it
/// is independent of the engine, and it fails loudly if a future frozen table
/// disagrees with the map.
#[test]
fn the_pinned_frequencies_are_constant_on_monodromy_orbits() {
    let mut checked = 0usize;
    let mut moved_pairs: BTreeSet<(&'static str, &'static str)> = BTreeSet::new();
    let transports = frozen_line_transports();
    for (subgroup, row) in pinned_rows() {
        let Some(image_table) = transports
            .get(&(subgroup.parent_sg, row.parent_ml))
            .and_then(|source| source.full_step_image)
        else {
            continue;
        };
        if image_table.label == row.parent_ml {
            continue;
        }
        let frequencies = pinned_frequencies(&subgroup);
        let Some(other) = frequencies.get(image_table.label) else {
            // The subgroup's row list does not carry the image: the transport
            // claim is about the parent's frozen table, not about this row list,
            // so there is nothing to compare here.
            continue;
        };
        assert_eq!(
            other, &row.frequency,
            "ordinal {} {} and its monodromy image {} must carry the same pinned frequency",
            subgroup.ordinal, row.parent_ml, image_table.label
        );
        checked += 1;
        moved_pairs.insert((row.parent_ml, image_table.label));
    }
    assert_eq!(
        checked, 816,
        "pin every non-trivial monodromy frequency pair"
    );
    println!("monodromy-invariant pinned frequencies: {checked} ({moved_pairs:?})");
}

/// The decomposition of a line subduction as a comparable key: the block
/// geometry plus every target's multiplicity, dimension, label, source number
/// and component, and the stored child `k` the block was matched against.
fn key(result: &LineSubduction) -> Vec<String> {
    result
        .blocks()
        .iter()
        .map(|block| {
            let targets: Vec<String> = block
                .targets()
                .iter()
                .map(|target| {
                    format!(
                        "{}x{} {:?} {:?}",
                        target.multiplicity, target.dimension, target.ml, target.irnumber
                    )
                })
                .collect();
            format!(
                "k=({},{},{}) star={} arms={} dim={} little={} [{}]",
                block.stored_k().get(0),
                block.stored_k().get(1),
                block.stored_k().get(2),
                block.star_size(),
                block.arm_count(),
                block.block_dimension(),
                block.little_dimension(),
                targets.join("; ")
            )
        })
        .collect()
}

/// The engine transport contract, on every pinned parametric-k row, for the
/// primitive reciprocal-lattice step along that row's line.
///
/// ```text
/// decompose(alpha, t + 1/d)  ==  decompose(M_{v/d}(alpha), t)
/// ```
///
/// `d` is the largest integer for which `v/d` remains in the parent's
/// reciprocal lattice. Thus these calls describe the same wave vector modulo
/// the full lattice, including cases where `d > 1`; the label is what moves.
#[test]
fn the_engine_transports_every_primitive_line_reciprocal_step() {
    let official = official_line_parameter().expect("the official parameter is valid");
    let rows = pinned_rows();
    assert_eq!(rows.len(), 5_756, "the pinned line corpus has 5,756 rows");
    let transports = frozen_line_transports();
    let outcomes = rows
        .par_iter()
        .map(|(subgroup, row)| {
            let source = transports.get(&(subgroup.parent_sg, row.parent_ml))?;
            let image_table = source.primitive_image?;
            let parameter = official
                .checked_add(
                    Rat::new(1, source.primitive_parameter_steps as i128)
                        .expect("positive step index"),
                )
                .expect("shifted parameter");
            let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(subgroup) else {
                return None;
            };
            // The transported query uses the same frozen table at the first
            // reciprocal-equivalent parameter. Its character twist must equal
            // the monodromy image at the official parameter.
            let transported = subduce_line_at_parameter(
                subgroup,
                &embedding,
                line_table(subgroup.parent_sg, row.parent_ml).expect("frozen source exists"),
                parameter,
            )
            .unwrap_or_else(|error| {
                panic!(
                    "ordinal {} {} at t={parameter}: {error}",
                    subgroup.ordinal, row.parent_ml,
                )
            });
            let reference = subduce_line_at_parameter(subgroup, &embedding, image_table, official)
                .unwrap_or_else(|error| {
                    panic!(
                        "ordinal {} {} (monodromy image of {}) at t=1/4: {error}",
                        subgroup.ordinal, image_table.label, row.parent_ml
                    )
                });
            assert_eq!(
                key(&transported),
                key(&reference),
                "ordinal {} {}: the decomposition at t={parameter} must equal the one of its \
                 monodromy image {} at t=1/4",
                subgroup.ordinal,
                row.parent_ml,
                image_table.label
            );
            Some(image_table.label != row.parent_ml)
        })
        .collect::<Vec<_>>();
    let checked = outcomes.iter().filter(|outcome| outcome.is_some()).count();
    let moved = outcomes
        .iter()
        .filter(|outcome| **outcome == Some(true))
        .count();
    let undetermined = outcomes.len() - checked;
    assert_eq!(checked, 5_756, "every pinned row must be compared");
    assert_eq!(
        undetermined, 0,
        "transport may not silently skip pinned rows"
    );
    assert_eq!(
        moved, 816,
        "pin the number of transported non-identity rows"
    );
    println!(
        "engine transport: checked={checked} non-trivial_images={moved} \
         undetermined_images={undetermined}"
    );
}

/// The revoked reading, kept as a named counterexample: the engine must **not**
/// return the same decomposition for the same label one parameter step later
/// whenever the monodromy is non-trivial.
///
/// Without this, a regression that re-introduced the canonical wave vector (and
/// with it `M == 1`) would pass every other test in this file.
#[test]
fn a_parameter_step_is_not_the_identity_on_a_nontrivial_monodromy() {
    let shifted = stepped(1);
    let transports = frozen_line_transports();
    let official = official_line_parameter().expect("official");
    let witnesses = pinned_rows()
        .par_iter()
        .filter_map(|(subgroup, row)| {
            let source = transports.get(&(subgroup.parent_sg, row.parent_ml))?;
            let image_table = source.full_step_image?;
            if image_table.label == row.parent_ml {
                return None;
            }
            let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).ok()?;
            let source_table = line_table(subgroup.parent_sg, row.parent_ml)?;
            let reference =
                subduce_line_at_parameter(subgroup, &embedding, source_table, official).ok()?;
            let same_label =
                subduce_line_at_parameter(subgroup, &embedding, source_table, shifted).ok()?;
            (key(&same_label) != key(&reference)).then_some((subgroup.ordinal, row.parent_ml))
        })
        .collect::<Vec<_>>();
    assert_eq!(witnesses.len(), 40, "pin every same-label counterexample");
    println!(
        "same-label shift witnesses: {} ({:?})",
        witnesses.len(),
        &witnesses[..witnesses.len().min(5)]
    );
}

/// The unused-image check for the map itself: every label of every line parent
/// has an image, and the labels whose image the frozen data does not determine
/// are reported (never silently defaulted to the identity).
#[test]
fn every_frozen_label_has_a_reportable_image() {
    let mut total = 0usize;
    let mut ambiguous = 0usize;
    let mut missing = 0usize;
    for parent in line_parents() {
        for source in line_sources(parent) {
            let map = checked_monodromy(parent, &one_step(source));
            let image = map
                .image(source.label)
                .unwrap_or_else(|| panic!("SG {parent} {} has no image", source.label));
            total += 1;
            match image {
                LabelImage::Unique(_) => {}
                LabelImage::Ambiguous(_) => ambiguous += 1,
                LabelImage::Missing => missing += 1,
                LabelImage::UnsupportedShift => panic!(
                    "SG {parent} {} must support its own frozen line direction",
                    source.label
                ),
            }
            let orbit = map
                .orbit(source.label, 2)
                .unwrap_or_else(|| panic!("SG {parent} {} orbit", source.label));
            assert_eq!(orbit[0], source.label);
        }
    }
    assert_eq!(total, 73, "the frozen corpus has 73 sources");
    println!("monodromy images: total={total} ambiguous={ambiguous} missing={missing}");
}
