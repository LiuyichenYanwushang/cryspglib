//! Full-corpus census of the R9 magnetic embedding search (round 21).
//!
//! `tests/magnetic_geometry.rs` pins the search on a deterministic 851-row
//! stratified subset.  The R9 acceptance text asks for the whole corpus
//! ("首批 Type-I/III/IV 固定例通过后，再扫描全部 16,721 条记录，单独统计覆盖"),
//! and this example is that pass: it runs the production search
//! ([`embed_in_parent_conventional`]) on every magnetic isotropy record and
//! re-checks every returned pair.
//!
//! What the re-check adds over a plain re-run:
//!
//! * membership in `L = Z^3 + span_Z(generators)` -- used for the parent lattice
//!   `L_P` and for the source lattice `L_M` -- is decided by **enumerating the
//!   whole finite residue set** on the common denominator, not by the production
//!   lexicographic-minimum BFS in `canonical_translation`;
//! * the record's row lattice is rebuilt from [`basis_in_parent_conventional`]
//!   through [`Lattice`] rather than from the search's own local variables;
//! * every returned pair must additionally pass the public [`verify_embedding`].
//!
//! A returned pair is **contained** when the image of the full source
//! translation lattice lies in `L_P`, and **realising** when in addition
//! `A L_M = L_record`.  The census reports how many records have 0 / 1 / many of
//! each; `--gate` turns every invariant -- including the pinned coverage numbers
//! -- into an exit code.
//!
//! Usage:
//!
//! ```text
//! magnetic_embedding_census             # print the census
//! magnetic_embedding_census --gate      # exit 1 unless every invariant holds
//! ```
//!
//! `RAYON_NUM_THREADS` limits the parallelism (one record per task).

use cryspglib::irrep::magnetic_embedding::{
    MagneticOperation, basis_in_parent_conventional, embed_in_parent_conventional,
    full_translation_lattice, geometry_of, lattice_of_operations, magnetic_operations,
    translation_lattice, verify_embedding,
};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, Vec3R};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashSet};
use std::process::ExitCode;

/// One magnetic isotropy record, with its parent space group.
type Row = (u8, cryspglib::irrep::types::MagneticIsotropyRecord);

/// A map's three rows, as a hashable key.
type MapRows = ([Rat; 3], [Rat; 3], [Rat; 3]);

/// Per-record census result.
#[derive(Default, Clone)]
struct RecordStat {
    entries: usize,
    realising: usize,
    distinct_pairs: usize,
    distinct_maps: usize,
    /// Returned pairs whose image of `L_M` leaves the parent lattice: must be 0.
    containment_violations: usize,
    /// Returned pairs whose flag disagrees with the recomputed equality: must be 0.
    flag_mismatches: usize,
    /// Returned pairs the public verifier rejects: must be 0.
    verifier_failures: usize,
    /// The record could not be read or searched: must be 0.
    errors: usize,
}

/// Greatest common divisor of positive integers.
fn gcd(mut left: i128, mut right: i128) -> i128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.abs()
}

/// Exact rational from a floating value on the small crystallographic grid.
fn exact(value: f64) -> Option<Rat> {
    for denominator in 1..=12i128 {
        let scaled = value * denominator as f64;
        if (scaled - scaled.round()).abs() < 1e-9 {
            return Rat::new(scaled.round() as i128, denominator).ok();
        }
    }
    None
}

/// `Z^3 + span_Z(generators)` as an explicit set of fractional parts.
///
/// A vector is in the lattice exactly when its fractional part (componentwise,
/// modulo one) equals the fractional part of some integer combination of the
/// generators.  The coefficients only matter modulo the generator's own
/// denominator `d_i` (because `d_i * g_i` is an integer vector), so the whole
/// reachable set is enumerated once, in exact rational arithmetic.  This is a
/// different algorithm from the production lexicographic-minimum BFS in
/// `canonical_translation`; the two must agree, and
/// `the_residue_lattice_agrees_with_the_production_membership` pins that on a
/// grid of lattices and vectors.
struct ResidueLattice {
    reachable: HashSet<[Rat; 3]>,
}

impl ResidueLattice {
    fn new(generators: &[[Rat; 3]]) -> Self {
        let mut reachable = HashSet::new();
        enumerate_fractions(generators, 0, [Rat::ZERO; 3], &mut |sum| {
            reachable.insert(sum);
        });
        Self { reachable }
    }

    fn contains(&self, vector: [Rat; 3]) -> bool {
        self.reachable.contains(&fractional_part(vector))
    }
}

/// Componentwise fractional part in `[0, 1)`.
fn fractional_part(vector: [Rat; 3]) -> [Rat; 3] {
    vector.map(fractional_part_value)
}

/// One value's fractional part in `[0, 1)`.
fn fractional_part_value(value: Rat) -> Rat {
    let numerator = value.numerator();
    let denominator = value.denominator();
    let whole = numerator.div_euclid(denominator);
    Rat::new(numerator - whole * denominator, denominator).expect("in [0, 1)")
}

/// Visit every reachable fractional part of an integer combination.
fn enumerate_fractions(
    generators: &[[Rat; 3]],
    index: usize,
    accumulator: [Rat; 3],
    visit: &mut impl FnMut([Rat; 3]),
) {
    if index == generators.len() {
        visit(accumulator);
        return;
    }
    let generator = generators[index];
    let mut period = 1i128;
    for value in generator {
        period = period / gcd(period, value.denominator()) * value.denominator();
    }
    for coefficient in 0..period {
        let mut next = accumulator;
        for axis in 0..3 {
            let term = generator[axis]
                .checked_mul(Rat::from_integer(coefficient))
                .expect("exact product");
            next[axis] = fractional_part_value(
                accumulator[axis].checked_add(term).expect("exact sum"),
            );
        }
        enumerate_fractions(generators, index + 1, next, visit);
    }
}

/// The parent's listed non-zero pure translations, exactly.
fn parent_centring(sg: u8) -> Option<Vec<[Rat; 3]>> {
    let operations = query::symmetry_operations_of(sg).ok()?;
    let mut converted = Vec::with_capacity(operations.operations.len());
    for operation in &operations.operations {
        converted.push(MagneticOperation {
            rotation: operation.rotation,
            translation: [
                exact(operation.translation[0])?,
                exact(operation.translation[1])?,
                exact(operation.translation[2])?,
            ],
            time_reversal: false,
        });
    }
    Some(lattice_of_operations(&converted))
}

/// Census one record with the independent checker.
fn census(row: &Row, parent_lattice: &ResidueLattice) -> RecordStat {
    let (sg, record) = row;
    let mut stat = RecordStat::default();
    let Ok(geometry) = geometry_of(*sg, record) else {
        stat.errors += 1;
        return stat;
    };
    let Ok(set) = magnetic_operations(record.mag_sg) else {
        stat.errors += 1;
        return stat;
    };
    let Ok(embeddings) = embed_in_parent_conventional(&geometry, &set) else {
        stat.errors += 1;
        return stat;
    };
    stat.entries = embeddings.len();

    let source_generators = full_translation_lattice(&set);
    let source_lattice = ResidueLattice::new(&translation_lattice(&set));
    let Ok(record_rows) = basis_in_parent_conventional(&geometry) else {
        stat.errors += 1;
        return stat;
    };
    let Ok(record_lattice) = Lattice::new(Mat3R::new(record_rows)) else {
        stat.errors += 1;
        return stat;
    };
    let in_record = |vector: [Rat; 3]| -> bool {
        record_lattice
            .coordinates(&Vec3R::new(vector))
            .map(|coordinates| coordinates.as_array().iter().all(|value| value.is_integer()))
            .unwrap_or(false)
    };

    let mut pairs: HashSet<(MapRows, [Rat; 3])> = HashSet::new();
    let mut maps: HashSet<MapRows> = HashSet::new();
    for entry in &embeddings {
        let image = |generator: [Rat; 3]| -> Option<[Rat; 3]> {
            entry
                .map
                .checked_mul_vector(&Vec3R::new(generator))
                .ok()
                .map(|vector| *vector.as_array())
        };
        // Containment: the image of every generator of the full source lattice
        // lies in the parent lattice.
        let contained = source_generators
            .iter()
            .all(|generator| image(*generator).is_some_and(|v| parent_lattice.contains(v)));
        if !contained {
            stat.containment_violations += 1;
        }
        // Equality `A L_M = L_record`, recomputed from the record again.
        let forward = contained
            && source_generators
                .iter()
                .all(|generator| image(*generator).is_some_and(in_record));
        let reverse = entry
            .map
            .inverse()
            .map(|inverse| {
                record_rows.iter().all(|row| {
                    inverse
                        .checked_mul_vector(&Vec3R::new(*row))
                        .map(|pullback| source_lattice.contains(*pullback.as_array()))
                        .unwrap_or(false)
                })
            })
            .unwrap_or(false);
        if entry.realises_record_lattice != (forward && reverse) {
            stat.flag_mismatches += 1;
        }
        if entry.realises_record_lattice {
            stat.realising += 1;
        }
        if verify_embedding(&geometry, &set, entry.map, entry.shift).is_err() {
            stat.verifier_failures += 1;
        }
        let map_key = (*entry.map.row(0), *entry.map.row(1), *entry.map.row(2));
        pairs.insert((map_key, entry.shift));
        maps.insert(map_key);
    }
    stat.distinct_pairs = pairs.len();
    stat.distinct_maps = maps.len();
    stat
}

/// Every magnetic isotropy record of the corpus, in database order.
fn rows() -> Vec<Row> {
    let mut rows = Vec::new();
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            rows.push((sg, row.subgroup));
        }
    }
    rows
}

/// Pinned coverage of the full corpus (round 21, HEAD with the round-20
/// containment gate and the equality flag).  `--gate` fails when any of these
/// moves, so the numbers in the ledger cannot drift away from the code.
const PINNED_ROWS: usize = 16_721;
const PINNED_ENTRIES: usize = 2_668_100;
const PINNED_DISTINCT_PAIRS: usize = 1_432_672;
const PINNED_DISTINCT_MAPS: usize = 333_500;
const PINNED_WITHOUT_EMBEDDING: usize = 8_294;
const PINNED_WITHOUT_REALISING: usize = 8_946;
const PINNED_SINGLE_REALISING: usize = 0;

fn main() -> ExitCode {
    let gate = std::env::args().any(|argument| argument == "--gate");
    let started = std::time::Instant::now();
    let rows = rows();
    let parent_lattices: Vec<ResidueLattice> = (1..=230u8)
        .map(|sg| {
            ResidueLattice::new(
                &parent_centring(sg).unwrap_or_else(|| panic!("parent {sg} is unavailable")),
            )
        })
        .collect();
    let stats: Vec<RecordStat> = rows
        .par_iter()
        .map(|row| census(row, &parent_lattices[usize::from(row.0) - 1]))
        .collect();

    let mut entries = 0usize;
    let mut distinct_pairs = 0usize;
    let mut distinct_maps = 0usize;
    let mut no_embedding = 0usize;
    let mut no_realising = 0usize;
    let mut unique_realising = 0usize;
    let mut with_realising = 0usize;
    let mut entries_histogram: BTreeMap<usize, usize> = BTreeMap::new();
    let mut realising_histogram: BTreeMap<usize, usize> = BTreeMap::new();
    let mut containment_violations = 0usize;
    let mut flag_mismatches = 0usize;
    let mut verifier_failures = 0usize;
    let mut errors = 0usize;
    for stat in &stats {
        entries += stat.entries;
        distinct_pairs += stat.distinct_pairs;
        distinct_maps += stat.distinct_maps;
        if stat.entries == 0 {
            no_embedding += 1;
        }
        if stat.realising == 0 {
            no_realising += 1;
        } else {
            with_realising += 1;
        }
        if stat.realising == 1 {
            unique_realising += 1;
        }
        *entries_histogram.entry(stat.entries).or_insert(0) += 1;
        *realising_histogram.entry(stat.realising).or_insert(0) += 1;
        containment_violations += stat.containment_violations;
        flag_mismatches += stat.flag_mismatches;
        verifier_failures += stat.verifier_failures;
        errors += stat.errors;
    }

    println!("rows {}", rows.len());
    println!("entries {entries}");
    println!("distinct_record_map_shift {distinct_pairs}");
    println!("distinct_record_map {distinct_maps}");
    println!("rows_without_embedding {no_embedding}");
    println!("rows_without_realising {no_realising}");
    println!("rows_with_realising {with_realising}");
    println!("rows_with_single_realising {unique_realising}");
    println!("containment_violations {containment_violations}");
    println!("flag_mismatches {flag_mismatches}");
    println!("verifier_failures {verifier_failures}");
    println!("errors {errors}");
    println!("entries_histogram {entries_histogram:?}");
    println!("realising_histogram {realising_histogram:?}");
    println!("elapsed {:.1}s", started.elapsed().as_secs_f64());

    if gate {
        let mut violations = Vec::new();
        let pinned: [(&str, usize, usize); 7] = [
            ("rows", rows.len(), PINNED_ROWS),
            ("entries", entries, PINNED_ENTRIES),
            ("distinct (record, map, shift)", distinct_pairs, PINNED_DISTINCT_PAIRS),
            ("distinct (record, map)", distinct_maps, PINNED_DISTINCT_MAPS),
            ("rows without embedding", no_embedding, PINNED_WITHOUT_EMBEDDING),
            ("rows without realising", no_realising, PINNED_WITHOUT_REALISING),
            ("rows with a single realising map", unique_realising, PINNED_SINGLE_REALISING),
        ];
        for (label, actual, expected) in pinned {
            if actual != expected {
                violations.push(format!("{label}: {actual} != {expected}"));
            }
        }
        if containment_violations != 0 {
            violations.push(format!("containment violations {containment_violations}"));
        }
        if flag_mismatches != 0 {
            violations.push(format!("flag mismatches {flag_mismatches}"));
        }
        if verifier_failures != 0 {
            violations.push(format!("verifier failures {verifier_failures}"));
        }
        if errors != 0 {
            violations.push(format!("record errors {errors}"));
        }
        if no_realising < no_embedding {
            violations.push("a row with no embedding cannot have a realising one".to_string());
        }
        if with_realising + no_realising != rows.len() {
            violations.push("coverage accounting does not add up".to_string());
        }
        if !violations.is_empty() {
            for violation in &violations {
                eprintln!("magnetic_embedding_census: {violation}");
            }
            return ExitCode::from(1);
        }
        println!("VERDICT complete scope=magnetic-embedding-corpus gate=--gate");
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use cryspglib::irrep::magnetic_embedding::canonical_translation;

    fn rational(numerator: i128, denominator: i128) -> Rat {
        Rat::new(numerator, denominator).expect("fraction")
    }

    /// The independent residue enumeration must decide exactly what the
    /// production membership primitive decides.  This is the self-check that
    /// caught the round-21 scaling bug (`ResidueLattice` used to scale the query
    /// vector by an integer division that ignored its own denominator, so every
    /// vector looked like a lattice vector when the generators were integral).
    #[test]
    fn the_residue_lattice_agrees_with_the_production_membership() {
        let lattices: Vec<Vec<[Rat; 3]>> = vec![
            Vec::new(),                                          // Z^3
            vec![[rational(1, 2), rational(1, 2), Rat::ZERO]],   // C centring
            vec![[rational(1, 2), rational(1, 2), rational(1, 2)]], // I centring
            vec![
                // R centring (both conventional vectors)
                [rational(2, 3), rational(1, 3), rational(1, 3)],
                [rational(1, 3), rational(2, 3), rational(2, 3)],
            ],
            vec![
                [rational(1, 2), Rat::ZERO, Rat::ZERO],
                [Rat::ZERO, rational(1, 3), Rat::ZERO],
                [Rat::ZERO, Rat::ZERO, rational(1, 6)],
            ],
            vec![[rational(1, 4), rational(1, 4), Rat::ZERO]],
        ];
        let values = [
            rational(-5, 6),
            rational(-2, 3),
            rational(-1, 2),
            rational(-1, 4),
            Rat::ZERO,
            rational(1, 4),
            rational(1, 3),
            rational(1, 2),
            rational(2, 3),
            Rat::ONE,
        ];
        let mut checked = 0usize;
        for generators in &lattices {
            let lattice = ResidueLattice::new(generators);
            for x in values {
                for y in values {
                    for z in values {
                        let vector = [x, y, z];
                        let mine = lattice.contains(vector);
                        let production = canonical_translation(vector, generators)
                            .iter()
                            .all(|value| value.is_zero());
                        assert_eq!(
                            mine, production,
                            "generators {generators:?} vector {vector:?}"
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert_eq!(checked, 6 * 1_000, "vectors compared");
    }

    /// The parent lattices of the corpus are read exactly: every listed centring
    /// translation is on the small grid and the generated lattice agrees with the
    /// production membership on the centring vector itself.
    #[test]
    fn every_parent_lattice_is_readable_and_contains_its_centring() {
        for sg in 1..=230u8 {
            let centring = parent_centring(sg).expect("parent operations");
            let lattice = ResidueLattice::new(&centring);
            for generator in &centring {
                assert!(
                    lattice.contains(*generator),
                    "sg {sg}: the generator {generator:?} must be in its own lattice"
                );
            }
            assert!(
                lattice.contains([Rat::ONE, Rat::ZERO, Rat::ZERO]),
                "sg {sg}: Z^3 is part of every parent lattice"
            );
        }
    }

    /// The residue enumeration must also agree with the production primitive on
    /// the **actual source lattices of the corpus**, not only on textbook ones.
    #[test]
    fn the_residue_enumeration_agrees_on_every_corpus_source_lattice() {
        let fractions = [
            Rat::ZERO,
            rational(1, 6),
            rational(1, 4),
            rational(1, 3),
            rational(1, 2),
            rational(2, 3),
            rational(3, 4),
            Rat::ONE,
        ];
        let mut checked = 0usize;
        let mut referenced = std::collections::BTreeSet::new();
        for sg in 1..=230u8 {
            for row in query::magnetic_isotropy_subgroups_of(sg) {
                referenced.insert(row.subgroup.mag_sg);
            }
        }
        for uni in referenced {
            let set = magnetic_operations(uni).expect("operations");
            let listed = translation_lattice(&set);
            let lattice = ResidueLattice::new(&listed);
            // Test every lattice vector, its negatives, and the grid of
            // fractional vectors around the generators' own denominators.
            let mut probes: Vec<[Rat; 3]> = vec![
                [Rat::ZERO; 3],
                [Rat::ONE, Rat::ZERO, Rat::ZERO],
                [Rat::ZERO, Rat::ONE, Rat::ZERO],
                [Rat::ZERO, Rat::ZERO, Rat::ONE],
            ];
            for generator in &listed {
                probes.push(*generator);
                probes.push(generator.map(|value| value.checked_neg().expect("negate")));
            }
            for x in fractions {
                for y in fractions {
                    for z in fractions {
                        probes.push([x, y, z]);
                    }
                }
            }
            for probe in probes {
                let mine = lattice.contains(probe);
                let production = canonical_translation(probe, &listed)
                    .iter()
                    .all(|value| value.is_zero());
                assert_eq!(mine, production, "UNI {uni}: {probe:?} with {listed:?}");
                checked += 1;
            }
        }
        assert!(checked > 500_000, "vectors compared: {checked}");
    }
}
