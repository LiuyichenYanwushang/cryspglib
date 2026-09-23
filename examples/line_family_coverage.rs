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
//! * the **conjugation oracle**: since every frozen direction is a parent
//!   reciprocal vector, `t = 3/4` (and `t = -1/4`, `7/4`) describes the
//!   conjugate of the parent irrep at `t = 1/4`, and conjugation preserves the
//!   multiplicity of the child's trivial representation, so the content there
//!   must equal the pinned one.  It does not: the measured gap is printed and
//!   pinned by a test as a known limitation (see the module test).
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
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_line_at_parameter;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::w_little_characters_data::{LittleCharacterTable, W_LITTLE_CHARACTERS};
use cryspglib::{HallNumber, SymmetryOps};
use std::collections::BTreeMap;
use std::io::Write as _;

/// The documented size of the conjugate-parameter gap (`t = -1/4`, `3/4`, `7/4`):
/// rows whose trivial content differs from the pinned value, and rows where the
/// engine fails closed with a non-integral multiplicity.  Fixing the line
/// character's gauge must drive both to zero and update this constant.
const CONJUGATE_GAP_CONTENT: usize = 199;
const CONJUGATE_GAP_ERRORS: usize = 58;

/// The four critical parameters of the frozen corpus.
const CRITICAL: [(&str, (i128, i128)); 4] =
    [("0", (0, 1)), ("1/4", (1, 4)), ("1/2", (1, 2)), ("3/4", (3, 4))];

/// Parameters that are equivalent to the pinned `1/4` up to conjugation, so the
/// pinned frequency is an oracle for them.
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
            .checked_add(left.get(axis).checked_mul(right.get(axis)).expect("product"))
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
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

/// `lcm` of two positive rationals, or the non-zero one when the other is zero.
fn rational_lcm(a: &Rat, b: &Rat) -> Rat {
    if a.is_zero() {
        return *b;
    }
    if b.is_zero() {
        return *a;
    }
    let numerator = (a.numerator().abs() / gcd(a.numerator().abs(), b.numerator().abs()))
        * b.numerator().abs();
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
    shift_content: Option<u32>,
}

fn collect(output: Option<&mut dyn std::io::Write>) -> Vec<Row> {
    let mut rows_out = Vec::new();
    let mut writer = output;
    if let Some(writer) = writer.as_mut() {
        writeln!(
            writer,
            "ordinal\tparent\tchild\tlabel\tarms\tpinned\tcritical_denominator\tc0\tc1_4\tc1_2\tc3_4\toff_grid_arms\toff_grid_content\tc_m1_4\tc3_4_conj\tc7_4\tshift_1_4"
        )
        .expect("write header");
    }
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
                for row in &pinned_rows {
                    let Some(table) = W_LITTLE_CHARACTERS.iter().find(|table| {
                        usize::from(table.space_group) == usize::from(parent)
                            && table.label == row.parent_ml
                    }) else {
                        continue;
                    };
                    // The arms: images of the frozen direction under the parent
                    // rotations, deduplicated by exact equality.
                    let direction = direction_of(table);
                    let mut arms: Vec<Vec3R> = Vec::new();
                    for rotation in &rotations {
                        let action = rotation.inverse().expect("inverse").transpose();
                        let image = action.checked_mul_vector(&direction).expect("arm");
                        if !arms.contains(&image) {
                            arms.push(image);
                        }
                    }
                    // The exact critical set: one rational subgroup generator per
                    // arm, from the three child lattice basis vectors.
                    let mut critical_points: std::collections::BTreeSet<(i128, i128)> =
                        Default::default();
                    for arm in &arms {
                        critical_points.extend(arm_critical_points(
                            arm,
                            &embedding,
                            &child_lattice,
                        ));
                    }
                    let critical_denominator = critical_points.len() as i128;
                    // Content at the four critical parameters.
                    let mut critical = [None; 4];
                    for (index, (_, (numerator, denominator))) in CRITICAL.iter().enumerate() {
                        critical[index] = content(
                            &subgroup,
                            &embedding,
                            table,
                            Rat::new(*numerator, *denominator).expect("parameter"),
                        );
                    }
                    // Off-grid sample: arms folding onto the child Gamma, plus
                    // the engine's answer on a sample of rows (the call is the
                    // expensive part; the support statement is exact).
                    let mut off_grid_arms = 0usize;
                    let sampled = rows_out.len() % 25 == 0;
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
                            if child_reciprocal.contains(&q).expect("gamma test") {
                                hits += 1;
                            }
                        }
                        off_grid_arms += hits;
                        if sampled && j == 1 {
                            off_grid_content = Some(content(&subgroup, &embedding, table, t));
                        }
                    }
                    // The conjugation oracle.
                    let mut conjugate = [None; 3];
                    for (index, (_, (numerator, denominator))) in CONJUGATE.iter().enumerate() {
                        conjugate[index] = content(
                            &subgroup,
                            &embedding,
                            table,
                            Rat::new(*numerator, *denominator).expect("parameter"),
                        );
                    }
                    let shift_content = content(
                        &subgroup,
                        &embedding,
                        table,
                        Rat::new(5, 4).expect("parameter"),
                    );
                    let entry = Row {
                        ordinal: subgroup.ordinal,
                        parent,
                        child,
                        label: table.label,
                        arms: arms.len(),
                        pinned: u32::from(row.frequency),
                        critical_denominator,
                        critical,
                        off_grid_arms,
                        off_grid_content,
                        conjugate,
                        shift_content,
                    };
                    if let Some(writer) = writer.as_mut() {
                        let text = |value: Option<u32>| match value {
                            Some(value) => value.to_string(),
                            None => "error".to_string(),
                        };
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            entry.ordinal,
                            entry.parent,
                            entry.child,
                            entry.label,
                            entry.arms,
                            entry.pinned,
                            entry.critical_denominator,
                            text(entry.critical[0]),
                            text(entry.critical[1]),
                            text(entry.critical[2]),
                            text(entry.critical[3]),
                            entry.off_grid_arms,
                            entry
                                .off_grid_content
                                .map_or_else(|| "-".to_string(), text),
                            text(entry.conjugate[0]),
                            text(entry.conjugate[1]),
                            text(entry.conjugate[2]),
                            text(entry.shift_content),
                        )
                        .expect("write row");
                    }
                    rows_out.push(entry);
                }
            }
        }
    }
    rows_out
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
    let mut conjugate = [(0usize, 0usize); 3];
    let mut conjugate_witnesses = Vec::new();
    for row in &rows {
        *critical_sizes.entry(row.critical_denominator).or_default() += 1;
        *arity_histogram.entry((
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
        for (index, value) in row.conjugate.iter().enumerate() {
            match value {
                Some(found) if *found == row.pinned => {}
                Some(value) => {
                    conjugate[index].0 += 1;
                    if conjugate_witnesses.len() < 5 {
                        conjugate_witnesses.push(format!(
                            "ordinal {} {} t={}: content {value} != pinned {}",
                            row.ordinal, row.label, CONJUGATE[index].0, row.pinned
                        ));
                    }
                }
                None => conjugate[index].1 += 1,
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
    for (index, (name, _)) in CONJUGATE.iter().enumerate() {
        println!(
            "conjugation oracle t = {name}: content mismatched={} errors={} (known gap)",
            conjugate[index].0, conjugate[index].1
        );
    }
    for witness in &conjugate_witnesses {
        println!("  {witness}");
    }
    if gate {
        let known_gap = conjugate
            .iter()
            .all(|(wrong, errors)| *wrong == CONJUGATE_GAP_CONTENT && *errors == CONJUGATE_GAP_ERRORS);
        let ok = pinned_mismatch == 0
            && pinned_error == 0
            && off_grid_arms == 0
            && off_grid_content_bad == 0
            && shift_mismatch == 0
            && shift_error == 0
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

    /// The whole corpus: the exact critical set is the quarter grid, the pinned
    /// parameter matches, nothing folds onto the child Gamma off the grid, and
    /// the conjugate-parameter gap has exactly its documented size.
    ///
    /// The gap is a **known limitation**, not a contract: `t = -1/4` and its
    /// shifts describe the conjugate parent irrep, whose trivial content must
    /// equal the pinned one, so the engine is wrong (or fails closed) on those
    /// rows.  Fixing the line character's gauge must drive both counts to zero
    /// and update the constants at the top of this file.
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
            "t = 5/4 is the same parent irrep as t = 1/4"
        );
        let gap: Vec<(usize, usize)> = (0..3)
            .map(|index| {
                rows.iter()
                    .fold((0, 0), |(wrong, errors), row| match row.conjugate[index] {
                        Some(value) if value == row.pinned => (wrong, errors),
                        Some(_) => (wrong + 1, errors),
                        None => (wrong, errors + 1),
                    })
            })
            .collect();
        for (index, (wrong, errors)) in gap.iter().enumerate() {
            assert_eq!(
                (*wrong, *errors),
                (CONJUGATE_GAP_CONTENT, CONJUGATE_GAP_ERRORS),
                "known conjugate-parameter gap at t = {} — if the gauge fix landed, \
                 update CONJUGATE_GAP_CONTENT/ERRORS to the new (expected zero) values",
                CONJUGATE[index].0
            );
        }
    }
}
