//! Monodromy of the frozen line sources under a reciprocal-lattice shift.
//!
//! The contract, its derivation and the fingerprint mechanism are documented in
//! [`cryspglib::irrep::line_monodromy`].  The tests here pin three separate
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
use cryspglib::irrep::line_monodromy::{
    LabelImage, LineMonodromy, complex_conjugation, line_direction, line_parents, line_sources,
    line_table, monodromy,
};
use cryspglib::irrep::subduction::star::decompose::{
    LineSubduction, official_line_parameter, subduce_line_at_parameter,
};
use cryspglib::irrep::subduction::{Rat, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::{LabelConvention, isotropy, query};
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

/// The four algebraic identities of the monodromy contract, on every parent that
/// carries line sources, with the shift taken along each source's own direction
/// (where the twist is a genuine little-group character).
///
/// The identities are asserted on the **uniquely determined** part of each map;
/// the frozen little-group tables do not always determine it (two sources can
/// share a fingerprint), and those labels are reported so the boundary is
/// visible instead of being guessed.
#[test]
fn the_monodromy_contract_identities_hold() {
    let zero = Vec3R::zero();
    let mut non_trivial: BTreeSet<(&'static str, &'static str)> = BTreeSet::new();
    let mut undetermined_total = 0usize;
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
        let identity = monodromy(parent, &zero);
        for (label, target) in unique_part(&identity) {
            assert_eq!(
                target, label,
                "SG {parent}: M_0 must be the identity ({label})"
            );
            identity_checked += 1;
        }
        undetermined_total += identity.undetermined().len();
        for v in directions {
            let forward = monodromy(parent, &v);
            let backward = monodromy(
                parent,
                &v.checked_neg().expect("negated direction"),
            );
            let doubled = monodromy(
                parent,
                &v.checked_add(&v).expect("doubled direction"),
            );
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
    assert!(
        identity_checked > 0 && inverse_checked > 0 && doubled_checked > 0,
        "the contract identities were never exercised"
    );
    assert!(
        !non_trivial.is_empty(),
        "no parent has a non-trivial monodromy: the corpus would not exercise the contract"
    );
    println!(
        "monodromy contract: identity_checks={identity_checked} inverse_checks={inverse_checked} \
         composition_checks={doubled_checked} conjugation_checks={conjugation_checked} \
         undetermined_labels={undetermined_total}"
    );
    println!("non-trivial monodromy moves: {non_trivial:?}");
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
fn pinned_rows() -> Vec<(isotropy::IsotropySubgroup, isotropy::OtherWaveVectorSubduction)> {
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

/// The trivial content of the pinned table, as a map from row label to the set
/// of frequencies the pinned rows carry for this subgroup.
fn pinned_frequencies(
    subgroup: &isotropy::IsotropySubgroup,
) -> BTreeMap<&'static str, u16> {
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
    for (subgroup, row) in pinned_rows() {
        let Some(table) = line_table(subgroup.parent_sg, row.parent_ml) else {
            continue;
        };
        let map = monodromy(subgroup.parent_sg, &one_step(table));
        let Some(image) = map.unique_image(row.parent_ml) else {
            continue;
        };
        if image == row.parent_ml {
            continue;
        }
        let frequencies = pinned_frequencies(&subgroup);
        let Some(other) = frequencies.get(image) else {
            // The subgroup's row list does not carry the image: the transport
            // claim is about the parent's frozen table, not about this row list,
            // so there is nothing to compare here.
            continue;
        };
        assert_eq!(
            other, &row.frequency,
            "ordinal {} {} and its monodromy image {} must carry the same pinned frequency",
            subgroup.ordinal, row.parent_ml, image
        );
        checked += 1;
        moved_pairs.insert((row.parent_ml, image));
    }
    assert!(
        checked > 0,
        "no pinned row exercised a non-trivial monodromy image"
    );
    println!(
        "monodromy-invariant pinned frequencies: {checked} ({moved_pairs:?})"
    );
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

/// The engine transport contract, on every pinned parametric-k row:
///
/// ```text
/// decompose(alpha, t + 1)  ==  decompose(M_v(alpha), t)
/// ```
///
/// Every frozen direction is a parent reciprocal lattice vector, so the two
/// calls describe the same band at the same point of the Brillouin zone; the
/// label is what moves, and `M_v` is that movement.  Comparing `alpha` with
/// itself at both parameters instead is the over-strong `M == 1` reading the
/// R6.1 gauge gate shipped and R6.2 revoked.
#[test]
fn the_engine_transports_a_parameter_step_by_the_monodromy_map() {
    let official = official_line_parameter().expect("the official parameter is valid");
    let shifted = stepped(1);
    let mut checked = 0usize;
    let mut moved = 0usize;
    let mut undetermined = 0usize;
    for (subgroup, row) in pinned_rows() {
        let Some(table) = line_table(subgroup.parent_sg, row.parent_ml) else {
            continue;
        };
        let map = monodromy(subgroup.parent_sg, &one_step(table));
        let Some(image) = map.unique_image(row.parent_ml) else {
            undetermined += 1;
            continue;
        };
        let Some(image_table) = line_table(subgroup.parent_sg, image) else {
            undetermined += 1;
            continue;
        };
        let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&subgroup) else {
            continue;
        };
        // The transported query: the *same* frozen table, one parameter step
        // later.  Its characters are the twisted ones, so its answer must be the
        // one the monodromy image carries at the official parameter.
        let transported = subduce_line_at_parameter(&subgroup, &embedding, table, shifted)
            .unwrap_or_else(|error| {
                panic!(
                    "ordinal {} {} at t=5/4: {error}",
                    subgroup.ordinal, row.parent_ml
                )
            });
        let reference = subduce_line_at_parameter(&subgroup, &embedding, image_table, official)
            .unwrap_or_else(|error| {
                panic!(
                    "ordinal {} {} (monodromy image of {}) at t=1/4: {error}",
                    subgroup.ordinal, image, row.parent_ml
                )
            });
        assert_eq!(
            key(&transported),
            key(&reference),
            "ordinal {} {}: the decomposition of {} at t=5/4 must equal the one of its \
             monodromy image {} at t=1/4",
            subgroup.ordinal,
            row.parent_ml,
            row.parent_ml,
            image
        );
        checked += 1;
        if image != row.parent_ml {
            moved += 1;
        }
    }
    assert!(
        checked > 1_000,
        "the transport must be checked on the pinned corpus, not on a sample"
    );
    assert!(
        moved > 0,
        "no checked row used a non-trivial monodromy image"
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
    let mut witnesses: Vec<(usize, &'static str)> = Vec::new();
    for (subgroup, row) in pinned_rows() {
        let Some(table) = line_table(subgroup.parent_sg, row.parent_ml) else {
            continue;
        };
        let map = monodromy(subgroup.parent_sg, &one_step(table));
        let Some(image) = map.unique_image(row.parent_ml) else {
            continue;
        };
        if image == row.parent_ml {
            continue;
        }
        let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&subgroup) else {
            continue;
        };
        let (Ok(reference), Ok(same_label)) = (
            subduce_line_at_parameter(&subgroup, &embedding, table, official_line_parameter().expect("official")),
            subduce_line_at_parameter(&subgroup, &embedding, table, shifted),
        ) else {
            continue;
        };
        if key(&same_label) != key(&reference) {
            witnesses.push((subgroup.ordinal, row.parent_ml));
        }
    }
    assert!(
        !witnesses.is_empty(),
        "the corpus must contain a row where the same-label shift is visibly wrong"
    );
    println!("same-label shift witnesses: {} ({:?})", witnesses.len(), &witnesses[..witnesses.len().min(5)]);
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
            let map = monodromy(parent, &one_step(source));
            let image = map
                .image(source.label)
                .unwrap_or_else(|| panic!("SG {parent} {} has no image", source.label));
            total += 1;
            match image {
                LabelImage::Unique(_) => {}
                LabelImage::Ambiguous(_) => ambiguous += 1,
                LabelImage::Missing => missing += 1,
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
