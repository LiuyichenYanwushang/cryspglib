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
    ParentDomain, child_exceptional_parameters, parent_domain, parent_reciprocal, rotation_set,
    verify_against_grid,
};
use cryspglib::irrep::subduction::{Rat, SubgroupEmbedding};
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
    labels: Vec<&'static str>,
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

    let per_record: Vec<Vec<Probe>> = if sequential {
        records.iter().map(|record| probe_record(record, &domains, official, generic)).collect()
    } else {
        records
            .par_iter()
            .map(|record| probe_record(record, &domains, official, generic))
            .collect()
    };
    let mut probes: Vec<Probe> = Vec::new();
    let mut probe_errors: Vec<String> = Vec::new();
    for (record, rows) in records.iter().zip(per_record) {
        if rows.is_empty() {
            probe_errors.push(format!("ordinal {} produced no probe", record.ordinal));
        }
        probes.extend(rows);
    }

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

    report(&domains, &records, &probes, &probe_errors);

    let mut violations: Vec<String> = probe_errors;
    check_invariants(&domains, &records, &probes, &mut violations);

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
                let mut labels: Vec<&'static str> = Vec::new();
                for row in rows {
                    if !labels.contains(&row.parent_ml) {
                        labels.push(row.parent_ml);
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

/// The partition of one record: the parent's exceptional parameters, the folded
/// child's exceptional parameters, the official anchor and one generic sample.
fn record_parameters(
    record: &Record,
    domain: &ParentDomain,
    table: &'static LittleCharacterTable,
    official: Rat,
    generic: Rat,
) -> Result<Vec<(Rat, usize)>, String> {
    let embedding = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup)
        .map_err(|error| format!("ordinal {}: {error}", record.ordinal))?;
    let direction = line_direction(table)
        .ok_or_else(|| format!("SG {} {}: unparsable direction", record.parent_sg, table.label))?;
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
    let folded = child_exceptional_parameters(&embedding, &direction)
        .map_err(|error| format!("ordinal {} child domain: {error}", record.ordinal))?;
    for (parameter, order) in folded {
        push(parameter, order, &mut parameters);
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
    Ok(ordered)
}

/// Probe one record at every parameter of its partition.
fn probe_record(
    record: &Record,
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    official: Rat,
    generic: Rat,
) -> Vec<Probe> {
    let mut out = Vec::new();
    for label in &record.labels {
        let Some(table) = line_table(record.parent_sg, label) else {
            continue;
        };
        let Some(domain) = domains.get(&(record.parent_sg, label)) else {
            continue;
        };
        let Ok(embedding) = SubgroupEmbedding::from_isotropy_subgroup(&record.subgroup) else {
            continue;
        };
        let Ok(parameters) = record_parameters(record, domain, table, official, generic) else {
            continue;
        };
        for (parameter, child_order) in parameters {
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
                detail,
            });
        }
    }
    out
}

/// Human-readable summary: the partition, the classes and the open boundary.
fn report(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    records: &[Record],
    probes: &[Probe],
    errors: &[String],
) {
    println!(
        "sources={} records={} probes={} probe_errors={}",
        domains.len(),
        records.len(),
        probes.len(),
        errors.len()
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

/// Every invariant the census claims, checked against the probes.
fn check_invariants(
    domains: &BTreeMap<(u8, &'static str), ParentDomain>,
    records: &[Record],
    probes: &[Probe],
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
    // 3. A trivial child little co-group is always answerable, and an
    //    unsupported probe must be explained by an enhanced child co-group.
    for probe in probes {
        if probe.child_order == 1 && probe.class == TargetClass::Unsupported {
            violations.push(format!(
                "ordinal {} {} t={}: unsupported although the child co-group is trivial",
                probe.ordinal, probe.label, probe.parameter
            ));
        }
        if probe.class == TargetClass::Unsupported && probe.child_order <= 1 {
            violations.push(format!(
                "ordinal {} {} t={}: unsupported with child co-group order {}",
                probe.ordinal, probe.label, probe.parameter, probe.child_order
            ));
        }
    }
    // 4. The official anchor is answered for every pinned row.
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
    // 8. The exact enumeration is cross-checked against the group itself on a
    //    uniform grid.  The two methods share only the lattice membership test.
    let mut grid_checks = 0usize;
    let mut grid_mismatches = 0usize;
    for ((sg, label), domain) in domains {
        let Ok(lattice) = parent_reciprocal(*sg) else {
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
            "the enumeration disagrees with the group in {grid_mismatches} of {grid_checks}              grid predicates"
        ));
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
