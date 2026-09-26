//! Exact parameter-domain census of the frozen parametric-k line sources.
//!
//! R6.5 replaces denominator sampling with an exact partition of the parameter
//! line.  For every frozen source the module
//! [`cryspglib::irrep::subduction::star::line_domain`] computes, in exact
//! rational arithmetic, the parameters `t in [0, 1)` at which some parent
//! operation joins the little group of `k(t) = t . v`; the same computation is
//! repeated after folding `v` through each isotropy record's embedding, which is
//! where the *child* little co-group of the folded star grows.
//!
//! This example joins the two halves with the engine itself: at every parameter
//! of the partition (plus one generic sample per source) it runs the production
//! decomposition and records
//!
//! * the class of the answer -- `stored` (pinned child rows answered it),
//!   `constructed` (the projective/bloch catalogue built the targets),
//!   `unsupported` (`MissingChildStarData`), or `error`;
//! * the child little co-group order at that parameter from the census;
//! * the production `ParameterKind`, which must equal the census's notion of an
//!   enhanced parameter;
//! * the trivial content, so the anchor rows stay comparable with the pinned
//!   table.
//!
//! Usage:
//!
//! ```text
//! line_domain_census                 # summary + parameter partition
//! line_domain_census --gate          # exit 1 unless every invariant holds
//! line_domain_census --output out.tsv
//! line_domain_census --sequential    # one thread (default: rayon over records)
//! ```
//!
//! The gate deliberately accepts `unsupported > 0`: the census exists to *name*
//! the unsupported domains, not to pretend they are closed.  `--require-covered`
//! turns a non-empty unsupported set into a failure, for the round that closes
//! them all.

use cryspglib::irrep::line_monodromy::{line_direction, line_table};
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, ParameterKind, official_line_parameter, subduce_line_at_parameter,
};
use cryspglib::irrep::subduction::star::line_domain::{
    ParentDomain, child_exceptional_parameters, little_co_group_order, minimal_parameter_step,
    minimal_parameter_step_via_coordinates, parent_domain, reciprocal_lattice,
    require_reciprocal_direction, rotation_set, verify_against_grid,
};
use cryspglib::irrep::subduction::{Rat, SubgroupEmbedding, Vec3R, fold_wave_vector};
use cryspglib::irrep::w_little_characters_data::{LittleCharacterTable, W_LITTLE_CHARACTERS};
use cryspglib::irrep::{LabelConvention, isotropy, query};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufWriter, Write as _};
use std::process::ExitCode;

const USAGE: &str = "\
line_domain_census [--gate] [--require-covered] [--sequential] [--output <path>]

  --gate               exit 1 unless the census invariants hold
  --require-covered    additionally fail when any (record, parameter) is unsupported
  --sequential         one thread (default parallel over isotropy records)
  --output <path>      write the per-probe table as TSV
";

/// The generic parameter used to prove that a non-exceptional point is answered
/// from the frozen line little group.  It is off the quarter grid and not one of
/// the exceptional parameters of any frozen source (`{0, 1/2}`).
const GENERIC_SAMPLE: (i128, i128) = (1, 7);

/// How the engine answered one (record, label, parameter) probe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum TargetClass {
    /// Pinned child rows carried the folded star.
    Stored,
    /// The constructed catalogue (Bloch phase or projective targets) answered.
    Constructed,
    /// `MissingChildStarData`: the census must explain it as an enhanced child
    /// little co-group.
    Unsupported,
    /// Any other engine error: a hard failure of the census run.
    Error,
}

impl TargetClass {
    const fn label(self) -> &'static str {
        match self {
            Self::Stored => "stored",
            Self::Constructed => "constructed",
            Self::Unsupported => "unsupported",
            Self::Error => "error",
        }
    }
}

/// One probe of the census.
struct Probe {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    label: &'static str,
    parameter: Rat,
    parent_order: usize,
    parent_added: usize,
    child_order: usize,
    class: TargetClass,
    parameter_kind: Option<ParameterKind>,
    content: Option<u32>,
    detail: String,
}

/// One isotropy record that carries at least one parametric-k row.
struct Record {
    ordinal: usize,
    parent_sg: u8,
    child_sg: u8,
    subgroup: isotropy::IsotropySubgroup,
    /// Frozen source label and the pinned frequency of its row, so the anchor
    /// probe can be compared with the pinned table instead of only printed.
    labels: Vec<(&'static str, u16)>,
}

/// What one record contributed to the census.
struct RecordReport {
    probes: Vec<Probe>,
    /// `(checks, mismatches)` of the child-side grid cross-check.
    child_grid: (usize, usize),
    /// `(checks, mismatches)` of the child-side two-algorithm step comparison.
    child_algorithms: (usize, usize),
    /// The record's candidate parameters (parent union child).
    child_parameters: Vec<Rat>,
    failures: Vec<String>,
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => code,
        Err(message) => {
            eprintln!("line_domain_census: {message}");
            ExitCode::from(3)
        }
    }
}

fn run() -> Result<ExitCode, String> {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let gate = arguments.iter().any(|argument| argument == "--gate");
    let require_covered = arguments.iter().any(|argument| argument == "--require-covered");
    let sequential = arguments.iter().any(|argument| argument == "--sequential");
    if arguments.iter().any(|argument| argument == "--help" || argument == "-h") {
        print!("{USAGE}");
        return Ok(ExitCode::SUCCESS);
    }
    let output = arguments
        .iter()
        .position(|argument| argument == "--output")
        .map(|index| {
            arguments
                .get(index + 1)
                .cloned()
                .ok_or_else(|| "--output needs a path".to_string())
        })
        .transpose()?;

    let domains = source_domains()?;
    let records = records()?;
    let official = official_line_parameter().map_err(|error| error.to_string())?;
    let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1).map_err(|e| e.to_string())?;

    let per_record: Vec<RecordReport> = if sequential {
        records.iter().map(|record| probe_record(record, &domains, official, generic)).collect()
    } else {
        records
            .par_iter()
            .map(|record| probe_record(record, &domains, official, generic))
            .collect()
    };
    let mut probes: Vec<Probe> = Vec::new();
    let mut probe_errors: Vec<String> = Vec::new();
    let mut child_grid = (0usize, 0usize);
    let mut child_algorithms = (0usize, 0usize);
    let mut child_union: Vec<Rat> = Vec::new();
    for (record, report) in records.iter().zip(per_record) {
        if report.probes.is_empty() {
            probe_errors.push(format!("ordinal {} produced no probe", record.ordinal));
        }
        probe_errors.extend(report.failures);
        child_grid.0 += report.child_grid.0;
        child_grid.1 += report.child_grid.1;
        child_algorithms.0 += report.child_algorithms.0;
        child_algorithms.1 += report.child_algorithms.1;
        for parameter in report.child_parameters {
            if !child_union.contains(&parameter) {
                child_union.push(parameter);
            }
        }
        probes.extend(report.probes);
    }
    child_union = sorted_parameters(child_union);

    if let Some(path) = &output {
        let mut writer = BufWriter::new(
            std::fs::File::create(path).map_err(|error| format!("cannot create {path}: {error}"))?,
        );
        writeln!(
            writer,
            "ordinal\tparent_sg\tchild_sg\tlabel\tparameter\tparent_order\tparent_added\t\
             child_order\tclass\tparameter_kind\tcontent\tdetail"
        )
        .map_err(|error| error.to_string())?;
        for probe in &probes {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                probe.ordinal,
                probe.parent_sg,
                probe.child_sg,
                probe.label,
                probe.parameter,
                probe.parent_order,
                probe.parent_added,
                probe.child_order,
                probe.class.label(),
                match probe.parameter_kind {
                    Some(ParameterKind::LineIrrep) => "line",
                    Some(ParameterKind::Formal) => "formal",
                    None => "-",
                },
                probe
                    .content
                    .map_or_else(|| "-".to_string(), |value| value.to_string()),
                probe.detail,
            )
            .map_err(|error| error.to_string())?;
        }
    }

    let (algorithm_checks, algorithm_mismatches, algorithm_failures) =
        step_algorithm_cross_check(&domains);
    let evidence = CensusEvidence {
        child_grid,
        child_algorithms,
        child_union: &child_union,
        algorithm_checks,
        algorithm_mismatches,
    };
    report(&domains, &records, &probes, &probe_errors, &evidence);

    let mut violations: Vec<String> = probe_errors;
    violations.extend(algorithm_failures);
    check_invariants(&domains, &records, &probes, &evidence, &mut violations);

    let unsupported = probes
        .iter()
        .filter(|probe| probe.class == TargetClass::Unsupported)
        .count();
    let errors = probes.iter().filter(|probe| probe.class == TargetClass::Error).count();
    println!(
        "probes={} stored={} constructed={} unsupported={} errors={}",
        probes.len(),
        probes.iter().filter(|p| p.class == TargetClass::Stored).count(),
        probes.iter().filter(|p| p.class == TargetClass::Constructed).count(),
        unsupported,
        errors
    );
    for violation in &violations {
        eprintln!("census violation: {violation}");
    }
    if gate || require_covered {
        if !violations.is_empty() {
            println!("gate: FAILED ({} violation(s))", violations.len());
            return Ok(ExitCode::from(1));
        }
        if require_covered && unsupported > 0 {
            println!("gate: FAILED ({unsupported} unsupported probe(s), --require-covered)");
            return Ok(ExitCode::from(1));
        }
        println!("gate: ok");
    }
    Ok(ExitCode::SUCCESS)
}

/// The exact domain of every frozen source, checked against its own table.
fn source_domains() -> Result<BTreeMap<(u8, &'static str), ParentDomain>, String> {
    let mut domains = BTreeMap::new();
    for table in W_LITTLE_CHARACTERS {
        let domain = parent_domain(table).map_err(|error| {
            format!("SG {} {}: parameter domain: {error}", table.space_group, table.label)
        })?;
        if domain.generic_order != table.operations.len() {
            return Err(format!(
                "SG {} {}: the generic stabiliser has order {} but the frozen table lists {} \
                 operation(s)",
                table.space_group,
                table.label,
                domain.generic_order,
                table.operations.len()
            ));
        }
        domains.insert((table.space_group, table.label), domain);
    }
    Ok(domains)
}

/// Every isotropy record that carries parametric-k rows, with its labels.
fn records() -> Result<Vec<Record>, String> {
    let mut out = Vec::new();
    for sg in 1..=230u8 {
        for record in query::irreps_of(sg) {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
                .map_err(|error| format!("SG {sg} {}: {error}", record.ml))?;
            for subgroup in subgroups {
                let Ok(rows) = subgroup.other_wave_vector_subduction() else {
                    continue;
                };
                let mut labels: Vec<(&'static str, u16)> = Vec::new();
                for row in rows {
                    match labels.iter().find(|(label, _)| *label == row.parent_ml) {
                        Some((_, frequency)) => {
                            if *frequency != row.frequency {
                                return Err(format!(
                                    "ordinal {}: label {} carries two pinned frequencies \
                                     ({frequency} and {})",
                                    subgroup.ordinal, row.parent_ml, row.frequency
                                ));
                            }
                        }
                        None => labels.push((row.parent_ml, row.frequency)),
                    }
                }
                if labels.is_empty() {
                    continue;
                }
                out.push(Record {
                    ordinal: subgroup.ordinal,
                    parent_sg: sg,
                    child_sg: u8::try_from(subgroup.record.sg).unwrap_or(0),
                    subgroup,
                    labels,
                });
            }
        }
    }
    Ok(out)
}

fn gcd(left: i128, right: i128) -> i128 {
    let (mut a, mut b) = (left.abs(), right.abs());
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    a
}

/// Exact ordering key for the census grid: `Rat` deliberately has no `Ord`, so
/// the parameters are compared through a common denominator.  The census grid is
/// bounded by the source directions (denominators at most 12 in the frozen
/// corpus), so the multiplication cannot overflow; a violation would be a
/// finding, not a silent wrap.
fn sorted_parameters(mut values: Vec<Rat>) -> Vec<Rat> {
    let mut common: i128 = 1;
    for value in &values {
        let divisor = gcd(common, value.denominator());
        common = (common / divisor)
            .checked_mul(value.denominator())
            .expect("the census parameter grid stays small");
    }
    values.sort_by_key(|value| value.numerator() * (common / value.denominator()));
    values
}

/// The parameters of one record: `(probed parameters with their parent order, the
/// child's own candidate parameters)`.
type RecordParameters = (Vec<(Rat, usize)>, Vec<Rat>);

/// The partition of one record: the parent's exceptional parameters, the folded
/// child's exceptional parameters, the official anchor and one generic sample.
fn record_parameters(
    record: &Record,
    domain: &ParentDomain,
    table: &'static LittleCharacterTable,
    official: Rat,
    generic: Rat,
) -> Result<RecordParameters, String> {
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup)
        .map_err(|error| format!("ordinal {}: {error}", record.ordinal))?;
    let mut parameters: Vec<(Rat, usize)> = Vec::new();
    let push = |parameter: Rat, order: usize, parameters: &mut Vec<(Rat, usize)>| {
        match parameters.iter_mut().find(|(value, _)| *value == parameter) {
            Some((_, slot)) => *slot = (*slot).max(order),
            None => parameters.push((parameter, order)),
        }
    };
    for entry in &domain.exceptional {
        push(entry.parameter, 0, &mut parameters);
    }
    let folded = child_exceptional_parameters(&embedding, table)
        .map_err(|error| format!("ordinal {} child domain: {error}", record.ordinal))?;
    let child_candidates: Vec<Rat> = folded.iter().map(|(parameter, _)| *parameter).collect();
    for (parameter, _) in &folded {
        push(*parameter, 0, &mut parameters);
    }
    for parameter in [official, generic] {
        push(parameter, 0, &mut parameters);
    }
    let mut ordered: Vec<(Rat, usize)> = Vec::with_capacity(parameters.len());
    let sorted = sorted_parameters(parameters.iter().map(|(value, _)| *value).collect());
    for value in sorted {
        let order = parameters
            .iter()
            .find(|(parameter, _)| *parameter == value)
            .map_or(0, |(_, order)| *order);
        ordered.push((value, order));
    }
    Ok((ordered, child_candidates))
}

/// Probe one record at every parameter of its partition.
///
/// The child little co-group order is recomputed at **every** probed parameter
/// with [`little_co_group_order`], not only at the parameters the child census
/// lists, so no probe ever carries an unknown order.  The child-side enumeration
/// is also cross-checked against the group on a uniform grid right here.
fn probe_record(
    record: &Record,
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    official: Rat,
    generic: Rat,
) -> RecordReport {
    let mut out = Vec::new();
    let mut failures = Vec::new();
    let mut child_grid = (0usize, 0usize);
    let mut child_algorithms = (0usize, 0usize);
    let mut child_parameters: Vec<Rat> = Vec::new();
    let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup) else {
        failures.push(format!("ordinal {}: embedding rejected", record.ordinal));
        return RecordReport {
            probes: out,
            child_grid,
            child_algorithms: (0, 0),
            child_parameters,
            failures,
        };
    };
    let Ok(child_reciprocal) = reciprocal_lattice(record.child_sg) else {
        failures.push(format!(
            "ordinal {}: child #{} has no reciprocal lattice",
            record.ordinal, record.child_sg
        ));
        return RecordReport {
            probes: out,
            child_grid,
            child_algorithms: (0, 0),
            child_parameters,
            failures,
        };
    };
    let Ok(child_rotations) = rotation_set(record.child_sg) else {
        failures.push(format!(
            "ordinal {}: child #{} has no rotation set",
            record.ordinal, record.child_sg
        ));
        return RecordReport {
            probes: out,
            child_grid,
            child_algorithms: (0, 0),
            child_parameters,
            failures,
        };
    };
    // The child frame must be the child's own: the lattice and the rotation set
    // are rebuilt from the embedding's subgroup, the lattice must be invariant
    // under every rotation used, and the record's subgroup number must agree.
    // Neither the grid check nor the two-algorithm comparison can see a wrong
    // frame (both are fed the same lattice), so this is the only control that
    // does; the mutations of the third review round are what it exists for.
    // Limitation, measured in the fourth round: a *sibling* space group with the
    // same reciprocal lattice and the same rotation set (109 of the 116 child
    // groups have one) is indistinguishable here.  That is immaterial because the
    // frame is consumed only through those two values; what the assertion rules
    // out is a frame that would change them.
    if record.child_sg != embedding.subgroup_sg() {
        failures.push(format!(
            "ordinal {}: record child #{} != embedding subgroup #{}",
            record.ordinal,
            record.child_sg,
            embedding.subgroup_sg()
        ));
    }
    match reciprocal_lattice(embedding.subgroup_sg()) {
        Ok(expected) => {
            for (name, left, right) in [
                ("rebuilt", &child_reciprocal, &expected),
                ("expected", &expected, &child_reciprocal),
            ] {
                for row in 0..3 {
                    let vector = Vec3R::new(*right.rows().row(row));
                    if !left.contains(&vector).unwrap_or(false) {
                        failures.push(format!(
                            "ordinal {}: the {name} child lattice is not the embedding's \
                             reciprocal lattice (row {row})",
                            record.ordinal
                        ));
                    }
                }
            }
        }
        Err(error) => failures.push(format!(
            "ordinal {}: no reciprocal lattice for child #{}: {error}",
            record.ordinal,
            embedding.subgroup_sg()
        )),
    }
    match rotation_set(embedding.subgroup_sg()) {
        Ok(expected) => {
            if child_rotations != expected {
                failures.push(format!(
                    "ordinal {}: the child rotation set is not the embedding subgroup's",
                    record.ordinal
                ));
            }
        }
        Err(error) => failures.push(format!(
            "ordinal {}: no rotation set for child #{}: {error}",
            record.ordinal,
            embedding.subgroup_sg()
        )),
    }
    for rotation in &child_rotations {
        let Ok(action) = cryspglib::irrep::subduction::Mat3R::from_ints(*rotation).inverse()
        else {
            failures.push(format!("ordinal {}: child rotation is singular", record.ordinal));
            continue;
        };
        let Ok(action) = action.transpose().inverse() else {
            failures.push(format!("ordinal {}: child rotation is singular", record.ordinal));
            continue;
        };
        for row in 0..3 {
            let basis = Vec3R::new(*child_reciprocal.rows().row(row));
            match action.checked_mul_vector(&basis) {
                Ok(image) => match child_reciprocal.contains(&image) {
                    Ok(true) => {}
                    Ok(false) => failures.push(format!(
                        "ordinal {}: child rotation {rotation:?} does not preserve the child \
                         lattice (row {row})",
                        record.ordinal
                    )),
                    Err(error) => failures.push(format!(
                        "ordinal {}: child lattice test failed: {error}",
                        record.ordinal
                    )),
                },
                Err(error) => failures.push(format!(
                    "ordinal {}: child rotation action failed: {error}",
                    record.ordinal
                )),
            }
        }
    }
    for (label, pinned) in &record.labels {
        let Some(table) = line_table(record.parent_sg, label) else {
            failures.push(format!(
                "ordinal {}: no frozen table for SG {} {label}",
                record.ordinal, record.parent_sg
            ));
            continue;
        };
        let Some(domain) = domains.get(&(record.parent_sg, label)) else {
            failures.push(format!(
                "ordinal {}: no census domain for SG {} {label}",
                record.ordinal, record.parent_sg
            ));
            continue;
        };
        let Ok(direction) = line_direction(table).ok_or(()) else {
            failures.push(format!(
                "ordinal {}: SG {} {label} has no parsable direction",
                record.ordinal, record.parent_sg
            ));
            continue;
        };
        // The folded direction is `q(t) = t . q1`; both the grid cross-check and
        // the per-parameter order use `q1`.
        let Ok(folded) = fold_wave_vector(embedding.transform(), &direction) else {
            failures.push(format!("ordinal {}: folding failed", record.ordinal));
            continue;
        };
        if !folded.is_zero()
            && let Err(error) =
                require_reciprocal_direction(&child_reciprocal, &folded, record.child_sg, "")
        {
            failures.push(format!(
                "ordinal {}: the folded direction is out of scope: {error}",
                record.ordinal
            ));
            continue;
        }
        // The child frame gets the same arithmetic-independent control as the
        // parent: the centring scan and the lattice coordinate map must produce
        // the same step for every child rotation.  This is what makes the
        // child-side claim independent rather than a second reading of the same
        // algebra (a wrong *frame* is caught by the frame checks, a wrong step
        // by this comparison).
        if !folded.is_zero() {
            for rotation in &child_rotations {
                let image = match cryspglib::irrep::subduction::Mat3R::from_ints(*rotation)
                    .inverse()
                    .and_then(|matrix| matrix.transpose().checked_mul_vector(&folded))
                {
                    Ok(image) => image,
                    Err(error) => {
                        failures.push(format!(
                            "ordinal {}: child rotation image: {error}",
                            record.ordinal
                        ));
                        continue;
                    }
                };
                let Ok(w) = image.checked_sub(&folded) else {
                    failures.push(format!(
                        "ordinal {}: child rotation difference failed",
                        record.ordinal
                    ));
                    continue;
                };
                match (
                    minimal_parameter_step(&child_reciprocal, &w),
                    minimal_parameter_step_via_coordinates(&child_reciprocal, &w),
                ) {
                    (Ok(scan), Ok(coordinates)) => {
                        child_algorithms.0 += 1;
                        if scan != coordinates {
                            child_algorithms.1 += 1;
                            failures.push(format!(
                                "ordinal {} child rotation {rotation:?}: centring scan {:?} != \
                                 coordinate map {:?}",
                                record.ordinal,
                                scan.map(|value| value.to_string()),
                                coordinates.map(|value| value.to_string())
                            ));
                        }
                    }
                    (scan, coordinates) => failures.push(format!(
                        "ordinal {} child rotation {rotation:?}: step failed ({:?} / {:?})",
                        record.ordinal,
                        scan.err().map(|error| error.to_string()),
                        coordinates.err().map(|error| error.to_string())
                    )),
                }
            }
        }
        for denominator in [24i128, 120] {
            match verify_against_grid(&child_reciprocal, &folded, &child_rotations, denominator) {
                Ok(check) => {
                    child_grid.0 += check.checks;
                    child_grid.1 += check.mismatches;
                }
                Err(error) => failures.push(format!(
                    "ordinal {} child grid 1/{denominator}: {error}",
                    record.ordinal
                )),
            }
        }
        let Ok((parameters, child_candidates)) =
            record_parameters(record, domain, table, official, generic)
        else {
            failures.push(format!(
                "ordinal {}: parameter partition failed for {label}",
                record.ordinal
            ));
            continue;
        };
        for parameter in child_candidates {
            if !child_parameters.contains(&parameter) {
                child_parameters.push(parameter);
            }
        }
        for (parameter, _) in parameters {
            let child_order = if folded.is_zero() {
                child_rotations.len()
            } else {
                match little_co_group_order(&child_reciprocal, &folded, &child_rotations, parameter)
                {
                    Ok(order) => order,
                    Err(error) => {
                        failures.push(format!(
                            "ordinal {} {label} t={parameter}: child order: {error}",
                            record.ordinal
                        ));
                        continue;
                    }
                }
            };
            let parent = domain
                .exceptional
                .iter()
                .find(|entry| entry.parameter == parameter);
            let result = subduce_line_at_parameter(
                &record.subgroup,
                &embedding,
                table,
                parameter,
            );
            let (class, parameter_kind, content, detail) = match result {
                Ok(result) => {
                    let targets: Vec<_> = result
                        .blocks()
                        .iter()
                        .flat_map(|block| block.targets())
                        .collect();
                    let stored = targets.iter().all(|target| target.irnumber.is_some());
                    let content = result.trivial_content().map_err(|error| error.to_string());
                    (
                        if stored { TargetClass::Stored } else { TargetClass::Constructed },
                        Some(result.parameter_kind()),
                        content.as_ref().ok().copied(),
                        format!(
                            "blocks={} targets={} {}",
                            result.blocks().len(),
                            targets.len(),
                            match &content {
                                Ok(value) => format!("content={value}"),
                                Err(error) => format!("content error: {error}"),
                            }
                        ),
                    )
                }
                Err(FullStarError::MissingChildStarData { sg, points, .. }) => (
                    TargetClass::Unsupported,
                    None,
                    None,
                    format!("missing child data: sg={sg} points={points}"),
                ),
                Err(error) => (TargetClass::Error, None, None, error.to_string()),
            };
            let anchor_mismatch = if parameter == official {
                match (content, *pinned) {
                    (Some(value), expected) if value != u32::from(expected) => Some(format!(
                        "anchor content {value} != pinned {expected}"
                    )),
                    (None, _) if class != TargetClass::Error => {
                        Some("anchor has no trivial content".to_string())
                    }
                    _ => None,
                }
            } else {
                None
            };
            out.push(Probe {
                ordinal: record.ordinal,
                parent_sg: record.parent_sg,
                child_sg: record.child_sg,
                label,
                parameter,
                parent_order: parent.map_or(domain.generic_order, |entry| entry.order),
                parent_added: parent.map_or(0, |entry| entry.added),
                child_order,
                class,
                parameter_kind,
                content,
                detail: match anchor_mismatch {
                    Some(note) => format!("{note}; {}", detail),
                    None => detail,
                },
            });
        }
    }
    child_parameters.sort_by_key(|value| value.numerator() * 10_080 / value.denominator());
    child_parameters.dedup();
    RecordReport {
        probes: out,
        child_grid,
        child_algorithms,
        child_parameters,
        failures,
    }
}

/// Human-readable summary: the partition, the classes and the open boundary.
/// The census evidence that does not live in a single probe.
struct CensusEvidence<'a> {
    /// `(checks, mismatches)` of the child-side grid cross-check over all records.
    child_grid: (usize, usize),
    /// `(checks, mismatches)` of the child-side two-algorithm comparison.
    child_algorithms: (usize, usize),
    /// The union of the records' child candidate parameters.
    child_union: &'a [Rat],
    /// How often the centring scan and the coordinate map route were compared,
    /// and how often they disagreed.
    algorithm_checks: usize,
    algorithm_mismatches: usize,
}

fn report(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    records: &[Record],
    probes: &[Probe],
    errors: &[String],
    evidence: &CensusEvidence,
) {
    println!(
        "sources={} records={} probes={} probe_errors={}",
        domains.len(),
        records.len(),
        probes.len(),
        errors.len()
    );
    println!(
        "child-side grid cross-check: {} predicates, {} mismatch(es); candidate parameters {{{}}}",
        evidence.child_grid.0,
        evidence.child_grid.1,
        evidence
            .child_union
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
    println!(
        "probes without a child order: {}",
        probes.iter().filter(|probe| probe.child_order == 0).count()
    );
    println!(
        "step algorithms: parent {} rotation comparisons / {} disagreement(s); \
         child {} / {}",
        evidence.algorithm_checks,
        evidence.algorithm_mismatches,
        evidence.child_algorithms.0,
        evidence.child_algorithms.1
    );
    let mut parent_shapes: BTreeMap<Vec<(String, usize, usize)>, usize> = BTreeMap::new();
    for domain in domains.values() {
        let shape: Vec<(String, usize, usize)> = domain
            .exceptional
            .iter()
            .map(|entry| (entry.parameter.to_string(), entry.order, entry.added))
            .collect();
        *parent_shapes.entry(shape).or_insert(0) += 1;
    }
    println!("parent-exceptional parameter shapes (source count):");
    for (shape, count) in &parent_shapes {
        let text: Vec<String> = shape
            .iter()
            .map(|(parameter, order, added)| format!("t={parameter} order={order} +{added}"))
            .collect();
        println!("  {count:>3} source(s): {}", text.join(", "));
    }
    let mut per_source: BTreeMap<(&'static str, u8), Vec<Rat>> = BTreeMap::new();
    for probe in probes {
        let entry = per_source.entry((probe.label, probe.parent_sg)).or_default();
        if !entry.contains(&probe.parameter) {
            entry.push(probe.parameter);
        }
    }
    for parameters in per_source.values_mut() {
        *parameters = sorted_parameters(std::mem::take(parameters));
    }
    let mut partition_shapes: BTreeMap<Vec<String>, usize> = BTreeMap::new();
    for parameters in per_source.values() {
        let shape: Vec<String> = parameters.iter().map(|value| value.to_string()).collect();
        *partition_shapes.entry(shape).or_insert(0) += 1;
    }
    println!("probed parameter sets (source count):");
    for (shape, count) in partition_shapes.iter().take(12) {
        println!("  {count:>3} source(s): {{{}}}", shape.join(", "));
    }
    if partition_shapes.len() > 12 {
        println!("  ... {} more distinct sets", partition_shapes.len() - 12);
    }
    let mut classes: BTreeMap<(&'static str, String), usize> = BTreeMap::new();
    for probe in probes {
        *classes
            .entry((probe.class.label(), probe.parameter.to_string()))
            .or_insert(0) += 1;
    }
    println!("class by parameter:");
    for ((class, parameter), count) in &classes {
        println!("  t={parameter:<6} {class:<12} {count}");
    }
    let unsupported: Vec<&Probe> = probes
        .iter()
        .filter(|probe| probe.class == TargetClass::Unsupported)
        .collect();
    if !unsupported.is_empty() {
        let mut by_child: BTreeMap<(u8, usize), usize> = BTreeMap::new();
        for probe in &unsupported {
            *by_child.entry((probe.child_sg, probe.child_order)).or_insert(0) += 1;
        }
        println!("unsupported probes by (child space group, child co-group order):");
        for ((child, order), count) in &by_child {
            println!("  child #{child} order={order}: {count}");
        }
        for probe in unsupported.iter().take(5) {
            println!(
                "  witness ordinal {} {} t={} child #{}: {}",
                probe.ordinal, probe.label, probe.parameter, probe.child_sg, probe.detail
            );
        }
    }
}

/// Compare the centring scan with the coordinate-map route for every source
/// rotation: an arithmetic-independent control of the step, since the two
/// algorithms share no code beyond the rational helpers.
fn step_algorithm_cross_check(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
) -> (usize, usize, Vec<String>) {
    let mut checks = 0usize;
    let mut mismatches = 0usize;
    let mut failures = Vec::new();
    for ((sg, label), domain) in domains {
        let Ok(lattice) = reciprocal_lattice(*sg) else {
            failures.push(format!("SG {sg} {label}: no reciprocal lattice"));
            continue;
        };
        let Ok(rotations) = rotation_set(*sg) else {
            failures.push(format!("SG {sg} {label}: no rotation set"));
            continue;
        };
        for rotation in &rotations {
            let image = match cryspglib::irrep::subduction::Mat3R::from_ints(*rotation)
                .inverse()
                .and_then(|matrix| matrix.transpose().checked_mul_vector(&domain.direction))
            {
                Ok(image) => image,
                Err(error) => {
                    failures.push(format!("SG {sg} {label}: rotation image: {error}"));
                    continue;
                }
            };
            let Ok(w) = image.checked_sub(&domain.direction) else {
                failures.push(format!("SG {sg} {label}: rotation difference failed"));
                continue;
            };
            match (
                minimal_parameter_step(&lattice, &w),
                minimal_parameter_step_via_coordinates(&lattice, &w),
            ) {
                (Ok(scan), Ok(coordinates)) => {
                    checks += 1;
                    if scan != coordinates {
                        mismatches += 1;
                        failures.push(format!(
                            "SG {sg} {label} rotation {rotation:?}: centring scan {:?} != \
                             coordinate map {:?}",
                            scan.map(|value| value.to_string()),
                            coordinates.map(|value| value.to_string())
                        ));
                    }
                }
                (scan, coordinates) => failures.push(format!(
                    "SG {sg} {label} rotation {rotation:?}: step failed ({:?} / {:?})",
                    scan.err().map(|error| error.to_string()),
                    coordinates.err().map(|error| error.to_string())
                )),
            }
        }
    }
    (checks, mismatches, failures)
}

/// Every invariant the census claims, checked against the probes.
fn check_invariants(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    records: &[Record],
    probes: &[Probe],
    evidence: &CensusEvidence,
    violations: &mut Vec<String>,
) {
    // 1. The frozen tables are the generic stabiliser (already enforced while
    //    building the domains; repeated here so a future refactor cannot drop it).
    for ((sg, label), domain) in domains {
        let Some(table) = line_table(*sg, label) else {
            violations.push(format!("SG {sg} {label}: no frozen table"));
            continue;
        };
        if domain.generic_order != table.operations.len() {
            violations.push(format!(
                "SG {sg} {label}: generic order {} != frozen operations {}",
                domain.generic_order,
                table.operations.len()
            ));
        }
        if domain.exceptional.is_empty() {
            violations.push(format!("SG {sg} {label}: no exceptional parameter at all"));
        }
        if !domain
            .exceptional
            .iter()
            .any(|entry| entry.parameter.is_zero())
        {
            violations.push(format!("SG {sg} {label}: t = 0 is not in the domain"));
        }
        for entry in &domain.exceptional {
            if entry.order < domain.generic_order {
                violations.push(format!(
                    "SG {sg} {label}: parameter {} has order {} below the generic {}",
                    entry.parameter, entry.order, domain.generic_order
                ));
            }
        }
    }
    // 2. The production `ParameterKind` agrees with the census's notion of an
    //    enhanced parameter.
    for probe in probes {
        let Some(kind) = probe.parameter_kind else {
            continue;
        };
        let enhanced = domains
            .get(&(probe.parent_sg, probe.label))
            .is_some_and(|domain| domain.is_exceptional(&probe.parameter));
        let formal = kind == ParameterKind::Formal;
        if formal != enhanced {
            violations.push(format!(
                "ordinal {} {} t={}: ParameterKind {:?} but census enhanced={enhanced}",
                probe.ordinal, probe.label, probe.parameter, kind
            ));
        }
    }
    // 3. Every probe carries a real child little co-group order (never the
    //    sentinel zero), a trivial co-group is always answerable, and an
    //    unsupported probe must be explained by an enhanced child co-group.
    for probe in probes {
        if probe.child_order == 0 {
            violations.push(format!(
                "ordinal {} {} t={}: no child little co-group order was computed",
                probe.ordinal, probe.label, probe.parameter
            ));
        }
        if probe.class == TargetClass::Unsupported && probe.child_order <= 1 {
            violations.push(format!(
                "ordinal {} {} t={}: unsupported although the child co-group has order {}",
                probe.ordinal, probe.label, probe.parameter, probe.child_order
            ));
        }
    }
    // 3b. An engine error is a hard failure of the census run, under every flag:
    //     `errors` is not allowed to be a pure diagnostic counter.
    for probe in probes.iter().filter(|probe| probe.class == TargetClass::Error) {
        violations.push(format!(
            "ordinal {} {} t={}: engine error: {}",
            probe.ordinal, probe.label, probe.parameter, probe.detail
        ));
    }
    // 3c. The child-side census is cross-checked against the child group on a
    //     uniform grid, and its corpus-wide candidate set is exactly the eighth
    //     grid measured for this corpus.  Both are assertions, not report lines.
    if evidence.child_grid.1 > 0 {
        violations.push(format!(
            "the child enumeration disagrees with the child group in {} of {} grid predicates",
            evidence.child_grid.1, evidence.child_grid.0
        ));
    }
    if evidence.child_grid.0 == 0 {
        violations.push("the child grid cross-check never ran".to_string());
    }
    let expected_child: Vec<Rat> = [
        (0i128, 1i128),
        (1, 8),
        (1, 4),
        (3, 8),
        (1, 2),
        (5, 8),
        (3, 4),
        (7, 8),
    ]
    .iter()
    .map(|(numerator, denominator)| Rat::new(*numerator, *denominator).expect("rational"))
    .collect();
    if evidence.child_union != expected_child.as_slice() {
        violations.push(format!(
            "the child candidate parameters are {:?}, expected the eighth grid {:?}",
            evidence
                .child_union
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>(),
            expected_child
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
        ));
    }
    if evidence.algorithm_mismatches > 0 {
        violations.push(format!(
            "the centring scan and the coordinate map disagree in {} of {} rotations",
            evidence.algorithm_mismatches, evidence.algorithm_checks
        ));
    }
    if evidence.algorithm_checks < 2_000 {
        violations.push(format!(
            "the two step algorithms were compared only {} times",
            evidence.algorithm_checks
        ));
    }
    if evidence.child_algorithms.1 > 0 {
        violations.push(format!(
            "the child step algorithms disagree in {} of {} comparisons",
            evidence.child_algorithms.1, evidence.child_algorithms.0
        ));
    }
    if evidence.child_algorithms.0 == 0 {
        violations.push("the child two-algorithm comparison never ran".to_string());
    }
    // 3e. The frozen directions are in scope for the residue representation.
    for ((sg, label), domain) in domains {
        match reciprocal_lattice(*sg).and_then(|lattice| {
            require_reciprocal_direction(&lattice, &domain.direction, *sg, label)
        }) {
            Ok(()) => {}
            Err(error) => violations.push(format!("SG {sg} {label}: {error}")),
        }
    }
    // 4. The official anchor is answered for every pinned row, and its trivial
    //    content equals the pinned frequency (compared while probing; the label
    //    detail carries the mismatch).
    let official = official_line_parameter().ok();
    let mut anchor = 0usize;
    for probe in probes {
        if Some(probe.parameter) != official {
            continue;
        }
        anchor += 1;
        if probe.class == TargetClass::Unsupported {
            violations.push(format!(
                "ordinal {} {}: the official anchor is unsupported",
                probe.ordinal, probe.label
            ));
        }
        if probe.detail.contains("!= pinned") {
            violations.push(format!(
                "ordinal {} {}: {detail}",
                probe.ordinal,
                probe.label,
                detail = probe.detail
            ));
        }
    }
    if anchor != 5_756 {
        violations.push(format!("the anchor probes {anchor} rows, expected 5756"));
    }
    // 5. Every probe that the engine answered carries a trivial content.
    for probe in probes {
        if probe.class != TargetClass::Error
            && probe.class != TargetClass::Unsupported
            && probe.content.is_none()
        {
            violations.push(format!(
                "ordinal {} {} t={}: answered without a trivial content",
                probe.ordinal, probe.label, probe.parameter
            ));
        }
    }
    // 6. No probe was dropped: every record contributed at least one.
    if records.is_empty() {
        violations.push("no isotropy record carries parametric-k rows".to_string());
    }
    // 7b. The frozen table is the generic stabiliser as a *set*, not only in
    //     order: compare its rotations with those that fix the direction exactly.
    for ((sg, label), domain) in domains {
        let Ok(rotations) = rotation_set(*sg) else {
            violations.push(format!("SG {sg} {label}: no rotation set"));
            continue;
        };
        let mut generic_flat: Vec<Vec<i32>> = Vec::new();
        for rotation in &rotations {
            let Ok(image) = cryspglib::irrep::subduction::Mat3R::from_ints(*rotation)
                .inverse()
                .and_then(|matrix| matrix.transpose().checked_mul_vector(&domain.direction))
            else {
                violations.push(format!("SG {sg} {label}: rotation image failed"));
                continue;
            };
            if image == domain.direction {
                generic_flat.push(rotation.iter().flatten().copied().collect());
            }
        }
        let Some(table) = line_table(*sg, label) else {
            violations.push(format!("SG {sg} {label}: no frozen table"));
            continue;
        };
        let mut frozen: Vec<Vec<i32>> = table
            .operations
            .iter()
            .map(|operation| {
                operation
                    .rotation
                    .iter()
                    .flatten()
                    .map(|value| i32::from(*value))
                    .collect()
            })
            .collect();
        generic_flat.sort();
        frozen.sort();
        if generic_flat != frozen {
            violations.push(format!(
                "SG {sg} {label}: the frozen rotation set is not the generic stabiliser \
                 (generic {} vs frozen {} rotations)",
                generic_flat.len(),
                frozen.len()
            ));
        }
    }
    // 8. The exact enumeration is cross-checked against the group itself on a
    //    uniform grid.  The two methods share only the lattice membership test.
    let mut grid_checks = 0usize;
    let mut grid_mismatches = 0usize;
    for ((sg, label), domain) in domains {
        let Ok(lattice) = reciprocal_lattice(*sg) else {
            violations.push(format!("SG {sg} {label}: no parent reciprocal lattice"));
            continue;
        };
        let Ok(rotations) = rotation_set(*sg) else {
            violations.push(format!("SG {sg} {label}: no rotation set"));
            continue;
        };
        for denominator in [24i128, 120] {
            match verify_against_grid(&lattice, &domain.direction, &rotations, denominator) {
                Ok(check) => {
                    grid_checks += check.checks;
                    grid_mismatches += check.mismatches;
                }
                Err(error) => violations.push(format!(
                    "SG {sg} {label}: grid check 1/{denominator}: {error}"
                )),
            }
        }
    }
    if grid_mismatches > 0 {
        violations.push(format!(
            "the enumeration disagrees with the group in {grid_mismatches} of \
             {grid_checks} grid predicates"
        ));
    }
    if grid_checks == 0 {
        violations.push("the parent grid cross-check never ran".to_string());
    }
    println!("grid cross-check: {grid_checks} predicates, {grid_mismatches} mismatch(es)");
    // 7. The generic sample is never a formal parameter and never unsupported.
    let generic = Rat::new(GENERIC_SAMPLE.0, GENERIC_SAMPLE.1);
    if let Ok(generic) = generic {
        for probe in probes.iter().filter(|probe| probe.parameter == generic) {
            if probe.parameter_kind == Some(ParameterKind::Formal) {
                violations.push(format!(
                    "ordinal {} {}: the generic sample is labelled Formal",
                    probe.ordinal, probe.label
                ));
            }
            if probe.class == TargetClass::Unsupported {
                violations.push(format!(
                    "ordinal {} {}: the generic sample is unsupported",
                    probe.ordinal, probe.label
                ));
            }
        }
    }
}
