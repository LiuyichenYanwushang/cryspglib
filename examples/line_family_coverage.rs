//! R6.2: the parameter-family coverage of the 5,756 pinned parametric-k rows.
//!
//! This is the executable evidence behind `docs/subduction-r6-coverage.md`.  It
//! computes, for every pinned line row and without using any engine internal:
//!
//! * the **exact critical set** `S = { t : some arm folds onto the child Gamma }`
//!   from the frozen direction, the parent's data-Hall rotations and the child
//!   lattice.  Arm `a` folds onto the child Gamma exactly when `t * (a . T l)` is
//!   an integer for the three child lattice basis vectors `l`, so `S` is a union
//!   of rational subgroups and can be enumerated exactly;
//! * the engine's trivial content at the four critical parameters
//!   `t = 0, 1/4, 1/2, 3/4` (with explicit errors counted separately);
//! * a dense **off-grid** sample (`t = j/24`, `j` not a multiple of six) where no
//!   arm may fold onto the child Gamma and the content must stay zero;
//! * the **conjugation oracle**: `t = -1/4` describes the conjugate of the
//!   parent irrep at `t = 1/4`, and conjugation preserves the multiplicity of
//!   the child's trivial representation, so the content there must equal the
//!   pinned one of the conjugate source.  `3/4` and `7/4` are one and two
//!   reciprocal lattice steps further along the *same* line, i.e. the monodromy
//!   images `M_v(C(alpha))`, `M_v^2(C(alpha))` at the anchor
//!   (`cryspglib::irrep::line_monodromy`, `docs/subduction-conventions.md` §16);
//!   because the pinned frequency is constant on monodromy orbits, the partner's
//!   pinned frequency remains the oracle for all three.
//!
//! Usage:
//!
//! ```text
//! line_family_coverage                 # summary
//! line_family_coverage --gate          # exit 1 unless the documented numbers reproduce
//! line_family_coverage --output out.tsv  # per-row evidence table
//! ```
use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::line_monodromy::{line_table, monodromy, operation_translation};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_line_at_parameter;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::w_little_characters_data::{LittleCharacterTable, W_LITTLE_CHARACTERS};
use cryspglib::{HallNumber, SymmetryOps};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::sync::Arc;

/// The conjugate-parameter gap on **real** little-group tables
/// (`t = -1/4`, `3/4`, `7/4`): rows whose trivial content differs from the
/// oracle, and rows where the engine fails closed with a non-integral
/// multiplicity.
///
/// History, three stages: R6.1 measured `187 / 58` against a *canonicalized* wave
/// vector and recorded it here as a known limitation; R6.2 stage 1 (this round's
/// first half) kept the numbers while it re-derived what the parameters mean; the
/// monodromy transport (`line_wave_vector` evaluating the frozen table at the
/// **raw** `t v`, plus the image label for a step) is the second half, and it
/// makes both counts **zero**: the frozen table is no longer paired with a wave
/// vector it does not describe, so the conjugate coset transports exactly and no
/// multiplicity is non-integral.
const CONJUGATE_GAP_CONTENT: usize = 0;
const CONJUGATE_GAP_ERRORS: usize = 0;

/// The character-level **transport** comparison: rows where the restricted parent
/// character at `t = 5/4` is not, operation by operation, the character the
/// monodromy image `M_v(label)` has at `t = 1/4` (tolerance 1e-9, all tables, real
/// and complex).
///
/// This is the derived relation (one lattice step along the line multiplies each
/// arm's term by `exp(2 pi i (R_i^-T v).T_h)`, which is exactly the twist that
/// turns table `alpha` into table `M_v(alpha)`), and it is strictly stronger than
/// the multiplicity count: a row can carry the right trivial content and still
/// have a wrong character.
///
/// History: R6.1 compared `chi(3/4)` with `conj(chi(1/4))` *without* the twist and
/// measured `223 / 5510`; that form is not the relation the convention implies
/// (`3/4 = -1/4 + 1` conjugates **and** shifts by one lattice step, and the
/// elementwise factor is arm-dependent), so it is replaced here by the transport
/// form at `5/4`, which is exact on the whole corpus.
const REAL_CHARACTER_GAP: usize = 0;
const REAL_CHARACTER_COMPARED: usize = 5756;

/// The same counts for the **complex** tables (SG 209/210 `DT3`/`DT4`).  Here
/// conjugating swaps the source, so the oracle is the pinned frequency of the
/// conjugate *partner* in the same isotropy record; the engine satisfies it
/// exactly, which is why the gap above is a real-table statement.
const COMPLEX_CONJUGATE_GAP_CONTENT: usize = 0;
const COMPLEX_CONJUGATE_GAP_ERRORS: usize = 0;

/// The four critical parameters of the frozen corpus.
const CRITICAL: [(&str, (i128, i128)); 4] = [
    ("0", (0, 1)),
    ("1/4", (1, 4)),
    ("1/2", (1, 2)),
    ("3/4", (3, 4)),
];

/// Parameters on the conjugate coset, `k(t) = -k(1/4) + n v`: `-1/4` is the
/// conjugate point itself and `3/4`, `7/4` are its first two monodromy images, so
/// the pinned frequency of the conjugate source is an oracle for all three
/// (constant on monodromy orbits; see `cryspglib::irrep::line_monodromy`).
const CONJUGATE: [(&str, (i128, i128)); 3] = [("-1/4", (-1, 4)), ("3/4", (3, 4)), ("7/4", (7, 4))];

fn direction_of(table: &LittleCharacterTable) -> Vec3R {
    let mut values = [Rat::ZERO; 3];
    for (axis, text) in table.direction.iter().enumerate() {
        let text = text.trim();
        let (numerator, denominator) = text.split_once('/').unwrap_or((text, "1"));
        values[axis] = Rat::new(numerator.parse().unwrap(), denominator.parse().unwrap()).unwrap();
    }
    Vec3R::new(values)
}

/// The primitive lattice of one space group, exact on the pinned twelfth grid.
fn lattice_of(sg: u8) -> Lattice {
    let basis = isotropy::parent_primitive_basis(sg).expect("pinned basis");
    let mut rows = [[Rat::ZERO; 3]; 3];
    for (i, row) in basis.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            rows[i][j] = Rat::from_grid(*value, 12).expect("twelfth grid");
        }
    }
    Lattice::new(Mat3R::new(rows)).expect("basis matrix")
}

fn dot(left: &Vec3R, right: &Vec3R) -> Rat {
    let mut total = Rat::ZERO;
    for axis in 0..3 {
        total = total
            .checked_add(
                left.get(axis)
                    .checked_mul(right.get(axis))
                    .expect("product"),
            )
            .expect("sum");
    }
    total
}

/// The positive generator of `{t : t * c` is an integer`}` for `c = p/q` in
/// lowest terms, or `None` when the condition holds for every `t` (`c = 0`).
fn subgroup_generator(c: &Rat) -> Option<Rat> {
    if c.is_zero() {
        return None;
    }
    Rat::new(c.denominator().abs(), c.numerator().abs()).ok()
}

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

/// `lcm` of two positive rationals, or the non-zero one when the other is zero.
fn rational_lcm(a: &Rat, b: &Rat) -> Rat {
    if a.is_zero() {
        return *b;
    }
    if b.is_zero() {
        return *a;
    }
    let numerator =
        (a.numerator().abs() / gcd(a.numerator().abs(), b.numerator().abs())) * b.numerator().abs();
    let denominator = gcd(a.denominator().abs(), b.denominator().abs());
    Rat::new(numerator, denominator).expect("lcm")
}

/// The exact critical set of one row: the union over arms of the points where
/// that arm folds onto the child Gamma.
///
/// For one arm the condition is an **intersection** over the three child basis
/// vectors (`t` must satisfy all three divisibility conditions), so the arm's
/// own set is the subgroup generated by the divisibility-`lcm` of the three
/// rational generators, which is `{k * n/d mod 1}` with `d` points.  Across arms
/// the support is a **union**, so the point sets are merged, not intersected:
/// merging the generators with an `lcm` would compute the intersection and
/// wrongly drop `t = 1/4` for the rows whose arms have different generators.
fn arm_critical_points(
    arm: &Vec3R,
    embedding: &SubgroupEmbedding,
    child_lattice: &Lattice,
) -> Vec<(i128, i128)> {
    let mut generator: Option<Rat> = None;
    for basis_row in 0..3 {
        let l = Vec3R::new(*child_lattice.rows().row(basis_row));
        let t_l = embedding
            .transform()
            .matrix()
            .checked_mul_vector(&l)
            .expect("T l");
        if let Some(value) = subgroup_generator(&dot(arm, &t_l)) {
            generator = Some(match generator {
                Some(existing) => rational_lcm(&existing, &value),
                None => value,
            });
        }
    }
    match generator {
        None => vec![(0, 1)],
        Some(value) => (0..value.denominator().abs())
            .map(|k| {
                let scaled = Rat::new(k * value.numerator().abs(), value.denominator().abs())
                    .expect("point");
                let reduced = Rat::new(
                    scaled.numerator().rem_euclid(scaled.denominator()),
                    scaled.denominator(),
                )
                .expect("reduced point");
                (reduced.numerator(), reduced.denominator())
            })
            .collect(),
    }
}

/// One row's four critical-parameter contents, errors kept distinct.
type CriticalArity = (Option<u32>, Option<u32>, Option<u32>, Option<u32>);

struct Row {
    ordinal: usize,
    parent: u8,
    child: u8,
    label: &'static str,
    arms: usize,
    pinned: u32,
    critical_denominator: i128,
    /// Content at `t = 0, 1/4, 1/2, 3/4`; `None` is an explicit engine error.
    critical: [Option<u32>; 4],
    /// Off-grid sample: arms folding onto the child Gamma, and the engine's
    /// content when it was asked (a sample of the rows).
    off_grid_arms: usize,
    off_grid_content: Option<Option<u32>>,
    /// The conjugation oracle at `t = -1/4, 3/4, 7/4`.
    conjugate: [Option<u32>; 3],
    /// Whether the frozen little-group table is real (`D* = D`), i.e. whether the
    /// oracle is the *own* pinned frequency rather than the partner's.
    real: bool,
    /// The pinned frequency of the exact conjugate source in the same record,
    /// when that record lists it. The corpus test pins a unique source partner
    /// for every frozen table, so `None` means that this record has no partner row.
    partner_pinned: Option<u32>,
    partner: Option<&'static str>,
    /// Character-level **transport** check: `Some((equal, max deviation))` when
    /// the row decomposes at `t = 5/4`, at `t = 1/4` for its monodromy image, and
    /// the image is uniquely determined; `None` otherwise.
    character_transport: Option<(bool, f64)>,
    shift_content: Option<u32>,
}

struct RowContext {
    subgroup: IsotropySubgroup,
    embedding: SubgroupEmbedding,
    child: u8,
    child_lattice: Lattice,
    child_reciprocal: Lattice,
    rotations: Vec<Mat3R>,
}

struct RowJob {
    context: Arc<RowContext>,
    table: &'static LittleCharacterTable,
    pinned: u32,
    partner: Option<&'static str>,
    partner_pinned: Option<u32>,
}

fn collect_jobs() -> Vec<RowJob> {
    let mut jobs = Vec::new();
    for sg in 1..=230u8 {
        for record in query::irreps_of(sg) {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            let Ok(subgroups) = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
            else {
                continue;
            };
            for subgroup in subgroups {
                let Ok(pinned_rows) = subgroup.other_wave_vector_subduction() else {
                    continue;
                };
                if pinned_rows.is_empty() {
                    continue;
                }
                let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&subgroup) else {
                    continue;
                };
                let parent = subgroup.parent_sg;
                let child = embedding.subgroup_sg();
                let child_lattice = lattice_of(child);
                let child_reciprocal = child_lattice.reciprocal().expect("child reciprocal");
                let hall = cryspglib::irrep::generated_data::SG_DATA_HALL[usize::from(parent)];
                let ops =
                    SymmetryOps::from_hall_number(HallNumber::try_from(usize::from(hall)).unwrap())
                        .expect("data-Hall operations");
                let rotations: Vec<Mat3R> = ops
                    .iter()
                    .map(|operation| Mat3R::from_ints(operation.rotation))
                    .collect();
                let context = Arc::new(RowContext {
                    subgroup,
                    embedding,
                    child,
                    child_lattice,
                    child_reciprocal,
                    rotations,
                });
                for row in &pinned_rows {
                    let Some(table) = W_LITTLE_CHARACTERS.iter().find(|table| {
                        usize::from(table.space_group) == usize::from(parent)
                            && table.label == row.parent_ml
                    }) else {
                        continue;
                    };
                    let partner = conjugate_partner(table);
                    let partner_pinned = partner.and_then(|partner| {
                        pinned_rows
                            .iter()
                            .find(|candidate| candidate.parent_ml == partner)
                            .map(|candidate| u32::from(candidate.frequency))
                    });
                    jobs.push(RowJob {
                        context: Arc::clone(&context),
                        table,
                        pinned: u32::from(row.frequency),
                        partner,
                        partner_pinned,
                    });
                }
            }
        }
    }
    jobs
}

fn compute_row(job: &RowJob, index: usize) -> Row {
    let context = &job.context;
    let subgroup = &context.subgroup;
    let embedding = &context.embedding;
    let table = job.table;
    // The arms: images of the frozen direction under the parent rotations,
    // deduplicated by exact equality.
    let direction = direction_of(table);
    let mut arms: Vec<Vec3R> = Vec::new();
    for rotation in &context.rotations {
        let action = rotation.inverse().expect("inverse").transpose();
        let image = action.checked_mul_vector(&direction).expect("arm");
        if !arms.contains(&image) {
            arms.push(image);
        }
    }
    // The exact critical set: one rational subgroup generator per arm, from
    // the three child lattice basis vectors.
    let mut critical_points: std::collections::BTreeSet<(i128, i128)> = Default::default();
    for arm in &arms {
        critical_points.extend(arm_critical_points(arm, embedding, &context.child_lattice));
    }
    let critical_denominator = critical_points.len() as i128;
    // Content at the four critical parameters.
    let mut critical = [None; 4];
    for (parameter_index, (_, (numerator, denominator))) in CRITICAL.iter().enumerate() {
        critical[parameter_index] = content(
            subgroup,
            embedding,
            table,
            Rat::new(*numerator, *denominator).expect("parameter"),
        );
    }
    // Off-grid sample: arms folding onto the child Gamma, plus the engine's
    // answer on one in every 25 rows. The row index is assigned before workers
    // start, so the sampling pattern remains deterministic under parallelism.
    let mut off_grid_arms = 0usize;
    let sampled = index.is_multiple_of(25);
    let mut off_grid_content = None;
    for j in 1..24i128 {
        if j % 6 == 0 {
            continue;
        }
        let t = Rat::new(j, 24).expect("parameter");
        let mut hits = 0usize;
        for arm in &arms {
            let mut scaled = [Rat::ZERO; 3];
            for (axis, value) in scaled.iter_mut().enumerate() {
                *value = t.checked_mul(arm.get(axis)).expect("scaled arm");
            }
            let q = embedding
                .transform()
                .matrix()
                .transpose()
                .checked_mul_vector(&Vec3R::new(scaled))
                .expect("folded arm");
            if context.child_reciprocal.contains(&q).expect("gamma test") {
                hits += 1;
            }
        }
        off_grid_arms += hits;
        if sampled && j == 1 {
            off_grid_content = Some(content(subgroup, embedding, table, t));
        }
    }
    // The conjugation oracle.
    let mut conjugate = [None; 3];
    for (parameter_index, (_, (numerator, denominator))) in CONJUGATE.iter().enumerate() {
        conjugate[parameter_index] = content(
            subgroup,
            embedding,
            table,
            Rat::new(*numerator, *denominator).expect("parameter"),
        );
    }
    let shift_content = content(
        subgroup,
        embedding,
        table,
        Rat::new(5, 4).expect("parameter"),
    );
    Row {
        ordinal: subgroup.ordinal,
        parent: subgroup.parent_sg,
        child: context.child,
        label: table.label,
        arms: arms.len(),
        pinned: job.pinned,
        critical_denominator,
        critical,
        off_grid_arms,
        off_grid_content,
        conjugate,
        real: is_real(table),
        partner_pinned: job.partner_pinned,
        partner: job.partner,
        character_transport: character_transport(subgroup, embedding, table),
        shift_content,
    }
}

fn write_rows(writer: &mut dyn std::io::Write, rows: &[Row]) {
    writeln!(
        writer,
        "ordinal\tparent\tchild\tlabel\treal\tpartner\tpartner_pinned\tarms\tpinned\tcritical_denominator\tc0\tc1_4\tc1_2\tc3_4\toff_grid_arms\toff_grid_content\tc_m1_4\tc3_4_conj\tc7_4\tshift_1_4"
    )
    .expect("write header");
    let text = |value: Option<u32>| match value {
        Some(value) => value.to_string(),
        None => "error".to_string(),
    };
    for entry in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            entry.ordinal,
            entry.parent,
            entry.child,
            entry.label,
            entry.real,
            entry.partner.unwrap_or("-"),
            entry
                .partner_pinned
                .map_or_else(|| "-".to_string(), |value| value.to_string()),
            entry.arms,
            entry.pinned,
            entry.critical_denominator,
            text(entry.critical[0]),
            text(entry.critical[1]),
            text(entry.critical[2]),
            text(entry.critical[3]),
            entry.off_grid_arms,
            entry.off_grid_content.map_or_else(|| "-".to_string(), text),
            text(entry.conjugate[0]),
            text(entry.conjugate[1]),
            text(entry.conjugate[2]),
            text(entry.shift_content),
        )
        .expect("write row");
    }
}

fn collect(output: Option<&mut dyn std::io::Write>) -> Vec<Row> {
    let jobs = collect_jobs();
    // Vec's indexed parallel iterator preserves input order on collect. Rows
    // are emitted below serially, so reports stay byte-for-byte deterministic.
    let rows = jobs
        .par_iter()
        .enumerate()
        .map(|(index, job)| compute_row(job, index))
        .collect::<Vec<_>>();
    if let Some(writer) = output {
        write_rows(writer, &rows);
    }
    rows
}

fn content(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
    t: Rat,
) -> Option<u32> {
    subduce_line_at_parameter(subgroup, embedding, table, t)
        .ok()
        .and_then(|result| result.trivial_content().ok())
}

/// Compare the restricted parent character at `t = 5/4`, operation by operation,
/// with the one the monodromy image of the label has at `t = 1/4`, through the
/// public reconstruction arrays.  `None` when either parameter fails or the image
/// is not determined.
fn character_transport(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    table: &'static LittleCharacterTable,
) -> Option<(bool, f64)> {
    let image = monodromy(table.space_group, &direction_of(table))
        .ok()?
        .unique_image(table.label)
        .and_then(|image| line_table(table.space_group, image))?;
    let anchor =
        subduce_line_at_parameter(subgroup, embedding, image, Rat::new(1, 4).ok()?).ok()?;
    let shifted =
        subduce_line_at_parameter(subgroup, embedding, table, Rat::new(5, 4).ok()?).ok()?;
    let (anchor_characters, _) = anchor.reconstruction();
    let (shifted_characters, _) = shifted.reconstruction();
    if anchor_characters.len() != shifted_characters.len() {
        return None;
    }
    let mut deviation = 0.0f64;
    for (left, right) in anchor_characters.iter().zip(shifted_characters) {
        deviation = deviation.max((left - right).norm());
    }
    Some((deviation <= 1e-9, deviation))
}

/// Whether the frozen little-group table is real, i.e. `D* = D`.
fn is_real(table: &LittleCharacterTable) -> bool {
    table
        .operations
        .iter()
        .all(|operation| operation.character[1] == 0)
}

type ExactCharacterKey = ([[i32; 3]; 3], [(i128, i128); 3], [i32; 2]);

/// Sorted exact operation/character keys, retaining translations as well as
/// rotations so different Seitz representatives cannot be conflated.
fn character_keys(table: &LittleCharacterTable, conjugate: bool) -> Option<Vec<ExactCharacterKey>> {
    let mut keys = table
        .operations
        .iter()
        .map(|operation| {
            let translation = operation_translation(operation)?;
            let imaginary = if conjugate {
                operation.character[1].checked_neg()?
            } else {
                operation.character[1]
            };
            Some((
                operation.rotation.map(|row| row.map(i32::from)),
                std::array::from_fn(|axis| {
                    let value = translation.get(axis);
                    (value.numerator(), value.denominator())
                }),
                [operation.character[0], imaginary],
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    keys.sort_unstable();
    Some(keys)
}

/// The conjugate partner label, if a same-parent source has the exact
/// conjugated character on the same operations. Missing data stays missing.
fn conjugate_partner(table: &LittleCharacterTable) -> Option<&'static str> {
    let wanted = character_keys(table, true)?;
    let mut found = None;
    for candidate in W_LITTLE_CHARACTERS {
        if candidate.space_group != table.space_group {
            continue;
        }
        let Some(actual) = character_keys(candidate, false) else {
            continue;
        };
        if actual == wanted {
            if found.is_some() {
                return None;
            }
            found = Some(candidate.label);
        }
    }
    found
}

fn main() -> std::process::ExitCode {
    let mut gate = false;
    let mut output: Option<String> = None;
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--gate" => gate = true,
            "--output" => output = arguments.next(),
            other => {
                eprintln!("unknown argument {other}");
                return std::process::ExitCode::from(2);
            }
        }
    }
    let mut file = output.map(|path| {
        std::fs::File::create(&path).unwrap_or_else(|error| panic!("create {path}: {error}"))
    });
    let rows = collect(file.as_mut().map(|file| file as &mut dyn std::io::Write));
    let mut critical_sizes: BTreeMap<i128, usize> = Default::default();
    let mut arity_histogram: BTreeMap<CriticalArity, usize> = Default::default();
    let mut pinned_mismatch = 0usize;
    let mut pinned_error = 0usize;
    let mut off_grid_arms = 0usize;
    let mut off_grid_content_bad = 0usize;
    let mut off_grid_error = 0usize;
    let mut off_grid_sampled = 0usize;
    let mut shift_mismatch = 0usize;
    let mut shift_error = 0usize;
    // Split by table type: a real table has `D* = D` (own pinned is the oracle),
    // a complex one has the conjugate *partner* as its oracle.
    let mut conjugate: BTreeMap<(&str, bool), (usize, usize, usize)> = Default::default();
    let mut complex_own_difference = 0usize;
    let mut character_compared = 0usize;
    let mut character_gap = 0usize;
    let mut character_witnesses = Vec::new();
    let mut conjugate_witnesses = Vec::new();
    for row in &rows {
        *critical_sizes.entry(row.critical_denominator).or_default() += 1;
        *arity_histogram
            .entry((
                row.critical[0],
                row.critical[1],
                row.critical[2],
                row.critical[3],
            ))
            .or_default() += 1;
        match row.critical[1] {
            Some(value) if value == row.pinned => {}
            Some(_) => pinned_mismatch += 1,
            None => pinned_error += 1,
        }
        off_grid_arms += row.off_grid_arms;
        if let Some(result) = row.off_grid_content {
            off_grid_sampled += 1;
            match result {
                Some(0) => {}
                Some(_) => off_grid_content_bad += 1,
                None => off_grid_error += 1,
            }
        }
        match row.shift_content {
            Some(value) if value == row.pinned => {}
            Some(_) => shift_mismatch += 1,
            None => shift_error += 1,
        }
        if !row.real && row.partner_pinned != Some(row.pinned) {
            complex_own_difference += 1;
        }
        if let Some((equal, deviation)) = row.character_transport {
            character_compared += 1;
            if !equal {
                character_gap += 1;
                if character_witnesses.len() < 5 {
                    character_witnesses.push(format!(
                        "ordinal {} {}: max |chi(1/4)* - chi(3/4)| = {deviation:.3e}",
                        row.ordinal, row.label
                    ));
                }
            }
        }
        for (index, value) in row.conjugate.iter().enumerate() {
            let entry = conjugate.entry((CONJUGATE[index].0, row.real)).or_default();
            let Some(expected) = row.partner_pinned else {
                entry.2 += 1;
                continue;
            };
            match value {
                Some(found) if *found == expected => {}
                Some(value) => {
                    entry.0 += 1;
                    if conjugate_witnesses.len() < 5 {
                        conjugate_witnesses.push(format!(
                            "ordinal {} {} (real={}) t={}: content {value} != oracle pinned {expected}",
                            row.ordinal, row.label, row.real, CONJUGATE[index].0
                        ));
                    }
                }
                None => entry.1 += 1,
            }
        }
    }
    println!("rows={}", rows.len());
    println!("exact critical set: generator denominator histogram {critical_sizes:?}");
    println!("content at (t=0, 1/4, 1/2, 3/4): arity histogram {arity_histogram:?}");
    println!(
        "pinned oracle at t = 1/4: mismatched={pinned_mismatch} errors={pinned_error} (of {} rows)",
        rows.len()
    );
    println!(
        "off-grid sample t = j/24 (j not a multiple of 6): arms folding onto the child Gamma {off_grid_arms} over {} points; engine checked on {off_grid_sampled} rows, content != 0: {off_grid_content_bad}, errors: {off_grid_error}",
        rows.len() * 18
    );
    println!(
        "equivalent parameter t = 5/4: content mismatched={shift_mismatch} errors={shift_error} (complete-decomposition comparison lives in audit_irrep_subduction)"
    );
    for (name, _) in CONJUGATE.iter() {
        for real in [true, false] {
            let (mismatch, error, missing) =
                conjugate.get(&(*name, real)).copied().unwrap_or((0, 0, 0));
            println!(
                "conjugation oracle t = {name} real_table={real}: content mismatched={mismatch} errors={error} partner_row_missing={missing}"
            );
        }
    }
    println!(
        "complex tables where the own pinned differs from the partner's (expected, not a gap): {complex_own_difference}"
    );
    println!(
        "character-level transport (t = 5/4 against the monodromy image at t = 1/4): compared={character_compared} not_transported={character_gap}"
    );
    for witness in &character_witnesses {
        println!("  {witness}");
    }
    for witness in &conjugate_witnesses {
        println!("  {witness}");
    }
    if gate {
        let real_gap = CONJUGATE.iter().all(|(name, _)| {
            conjugate
                .get(&(*name, true))
                .copied()
                .map(|(wrong, errors, _)| (wrong, errors))
                == Some((CONJUGATE_GAP_CONTENT, CONJUGATE_GAP_ERRORS))
        });
        let complex_gap = CONJUGATE.iter().all(|(name, _)| {
            conjugate
                .get(&(*name, false))
                .copied()
                .map(|(wrong, errors, _)| (wrong, errors))
                == Some((COMPLEX_CONJUGATE_GAP_CONTENT, COMPLEX_CONJUGATE_GAP_ERRORS))
        });
        let known_gap = real_gap && complex_gap;
        let character_documented =
            character_compared == REAL_CHARACTER_COMPARED && character_gap == REAL_CHARACTER_GAP;
        let ok = pinned_mismatch == 0
            && pinned_error == 0
            && off_grid_arms == 0
            && off_grid_content_bad == 0
            && shift_mismatch == 0
            && shift_error == 0
            && character_documented
            && known_gap
            && critical_sizes.len() == 1
            && critical_sizes.get(&4).copied() == Some(rows.len());
        if ok {
            println!("gate: ok (documented numbers reproduce)");
            return std::process::ExitCode::SUCCESS;
        }
        println!("gate: FAILED");
        return std::process::ExitCode::from(1);
    }
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_frozen_source_has_one_exact_conjugate_partner() {
        for source in W_LITTLE_CHARACTERS {
            assert!(
                conjugate_partner(source).is_some(),
                "SG {} {} has no unique exact conjugate source",
                source.space_group,
                source.label
            );
        }
    }

    /// The whole corpus: the exact critical set is the quarter grid, the pinned
    /// parameter matches, nothing folds onto the child Gamma off the grid, and
    /// the conjugate-parameter gap is **zero** on both the real and the complex
    /// tables — measured against the **correct** oracle, which is the source's
    /// own pinned frequency only when its little-group table is real and the
    /// conjugate partner's otherwise.
    ///
    /// The two halves of that zero: `t = -1/4` is the conjugate point itself
    /// (`chi(-k) = chi(k)*` for a real table, so the own pinned frequency is the
    /// oracle), and `3/4`, `7/4` are one and two reciprocal steps along the line,
    /// i.e. the monodromy images of the conjugate label.  The pinned frequency is
    /// constant on monodromy orbits (`tests/line_monodromy.rs`), so the partner's
    /// pinned frequency remains the oracle there; what changed in R6.2 is that
    /// the engine now evaluates the frozen table at the raw `t v` instead of a
    /// canonicalized representative, which is what made the 187 content
    /// mismatches and the 58 non-integral multiplicities disappear.
    #[test]
    fn the_documented_family_coverage_reproduces() {
        let rows = collect(None);
        assert_eq!(rows.len(), 5756, "the pinned line corpus");
        assert!(
            rows.iter().all(|row| row.critical_denominator == 4),
            "every row's arm set folds onto the child Gamma exactly on (1/4)Z"
        );
        assert!(
            rows.iter().all(|row| row.critical[1] == Some(row.pinned)),
            "t = 1/4 must reproduce every pinned frequency"
        );
        assert!(
            rows.iter().all(|row| row.off_grid_arms == 0),
            "no arm may reach the child Gamma off the quarter grid"
        );
        assert!(
            rows.iter().all(|row| match row.off_grid_content {
                Some(Some(0)) | None => true,
                Some(Some(_)) => false,
                Some(None) => true,
            }),
            "the engine's off-grid content is zero wherever it answers"
        );
        assert!(
            rows.iter().all(|row| row.shift_content == Some(row.pinned)),
            "t = 5/4 preserves the pinned trivial content after label transport"
        );
        let compared = rows
            .iter()
            .filter(|row| row.character_transport.is_some())
            .count();
        let gap = rows
            .iter()
            .filter(|row| matches!(row.character_transport, Some((false, _))))
            .count();
        assert_eq!(
            (compared, gap),
            (REAL_CHARACTER_COMPARED, REAL_CHARACTER_GAP),
            "character-level transport gap: chi(5/4) must equal the monodromy image's \
             chi(1/4) operation by operation",
        );
        for (index, (name, _)) in CONJUGATE.iter().enumerate() {
            for (real, expected) in [
                (true, (CONJUGATE_GAP_CONTENT, CONJUGATE_GAP_ERRORS)),
                (
                    false,
                    (COMPLEX_CONJUGATE_GAP_CONTENT, COMPLEX_CONJUGATE_GAP_ERRORS),
                ),
            ] {
                let (wrong, errors, missing) =
                    rows.iter()
                        .fold((0, 0, 0), |(wrong, errors, missing), row| {
                            if row.real != real {
                                return (wrong, errors, missing);
                            }
                            let Some(expected) = row.partner_pinned else {
                                return (wrong, errors, missing + 1);
                            };
                            match row.conjugate[index] {
                                Some(value) if value == expected => (wrong, errors, missing),
                                Some(_) => (wrong + 1, errors, missing),
                                None => (wrong, errors + 1, missing),
                            }
                        });
                assert_eq!(
                    (wrong, errors),
                    expected,
                    "conjugate-parameter gap at t = {name} on real_table={real} — \
                     if the transport fix landed, update the CONJUGATE_GAP_* constants",
                );
                // The four rows whose record does not list the conjugate partner
                // have no oracle and are reported, not counted as a gap.
                assert_eq!(
                    missing,
                    if real { 0 } else { 4 },
                    "partner rows missing at t = {name}",
                );
            }
        }
    }
}
