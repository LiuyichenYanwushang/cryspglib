//! Validate explicit `(setting, child_shift)` embedding conventions.
//!
//! This is a task-9 provenance tool, not part of the public API story: the
//! full-subduction engine freezes one convention per isotropy record, and the
//! convention has to be *derived* from the official ISOTROPY program before it
//! is frozen.  This example runs the engine's own validation on a convention
//! handed in from the outside, so a derivation script never has to re-implement
//! (and never has to be trusted to reproduce) that check.
//!
//! Input is a whitespace-separated table on standard input, one candidate per
//! line, `#` comments allowed:
//!
//! ```text
//! id ordinal u11 u12 u13 u21 u22 u23 u31 u32 u33 ud sx sy sz sd
//! ```
//!
//! `id` is an opaque token echoed back so a caller can tell candidates apart
//! (several candidates may name the same ordinal), `u` is the numerator of the
//! change of basis `U = u / ud`, and `(sx, sy, sz, sd)` the child origin shift
//! `(sx/sd, sy/sd, sz/sd)` -- exactly the two conventions a frozen row carries.
//!
//! Output is one line per candidate: `id OK <representatives> trivial=<n>` when
//! the engine accepts the convention, otherwise `id ERR <error>`.  `trivial` is
//! the multiplicity of the child's trivial representation in the **record's own
//! condensing irrep**, which is the property that makes the subgroup the
//! isotropy subgroup of that direction: a convention that validates but reports
//! zero has placed the subgroup on a different (conjugate) branch, and that is
//! invisible to the embedding check alone.  Three states are kept apart:
//! `trivial=<n>` is a computed multiplicity, `trivial=?` means the child's
//! trivial row or the condensing probe was not found, and `trivial=error` means
//! the decomposition itself failed -- the reason is printed on standard error
//! and never folded into a zero (a zero is a *result* here, so confusing the two
//! would misreport an engine gap as a conjugate branch).  A `--profile` line
//! whose subduction failed prints `profile=error` rather than a shorter profile
//! that could accidentally match a target.  Exit status is 0 when every
//! candidate was evaluated and 1 on malformed input.
//!
//! ```text
//! cargo run --release -p cryspglib --example probe_subduction_settings < candidates.tsv
//! ```

use std::collections::HashMap;
use std::io::{self, BufRead, Write};

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::decompose::subduce_full_star_with_embedding;
use cryspglib::irrep::types::IrrepRecord;
use cryspglib::irrep::{LabelConvention, subduction::SubgroupEmbedding};

/// The child's trivial irrep at Gamma: one-dimensional with an all-ones row.
fn trivial_child_label(sg: u8) -> Option<&'static str> {
    let mut found = None;
    for record in query::irreps_of(sg) {
        if record.spinor || record.dim != 1 {
            continue;
        }
        let k = record.k_vector();
        if k.numerators != [0, 0, 0] {
            continue;
        }
        let Ok(row) = record.ordinary_scalar_selected_arm_block_trace() else {
            continue;
        };
        if row.dimension() == 1 && row.values().iter().all(|value| (value - 1.0).norm() < 1e-9) {
            if found.is_some() {
                return None;
            }
            found = Some(record.ml);
        }
    }
    found
}

/// Trivial content of one parent irrep under one convention.
///
/// A failure is **not** a zero: zero is a real answer here (it is how a
/// convention on a conjugate branch shows up), so the error is propagated with
/// the context needed to diagnose it.
fn trivial_content(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
    trivial: &str,
) -> Result<u32, String> {
    let result = subduce_full_star_with_embedding(subgroup, embedding, probe).map_err(|error| {
        format!(
            "SG {} {} under ordinal {}: {error}",
            subgroup.parent_sg, probe.ml, subgroup.ordinal
        )
    })?;
    Ok(result
        .blocks()
        .iter()
        .map(|block| block.multiplicity(trivial))
        .sum())
}

/// Trivial content of the record's own condensing irrep under one convention.
///
/// `Ok(None)` means the answer is structurally unavailable (no unique trivial
/// child row, or the condensing probe is missing); `Err` means the engine was
/// asked and failed.
fn condensing_trivial(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
) -> Result<Option<u32>, String> {
    let Some(probe) = query::irreps_of(subgroup.parent_sg)
        .iter()
        .find(|record| record.ml == subgroup.irrep_ml)
    else {
        return Ok(None);
    };
    let Some(trivial) = trivial_child_label(embedding.subgroup_sg()) else {
        return Ok(None);
    };
    trivial_content(subgroup, embedding, probe, trivial).map(Some)
}

/// The printed token for one `condensing_trivial` outcome, plus the reason to
/// log when the answer is unavailable or failed.
fn render_trivial(outcome: Result<Option<u32>, String>) -> (String, Option<String>) {
    match outcome {
        Ok(Some(value)) => (value.to_string(), None),
        Ok(None) => ("?".to_string(), None),
        Err(reason) => ("error".to_string(), Some(reason)),
    }
}

/// `label=mult` for every parent irrep whose subduction contains the trivial
/// representation, which is exactly the table the audit compares against.
///
/// One failed subduction fails the whole profile: a profile with an entry
/// silently missing could equal a target profile for the wrong reason.
fn trivial_profile(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
) -> Result<String, String> {
    let Some(trivial) = trivial_child_label(embedding.subgroup_sg()) else {
        return Ok(String::from("?"));
    };
    let mut parts = Vec::new();
    for record in query::irreps_of(subgroup.parent_sg) {
        if record.spinor {
            continue;
        }
        let value = trivial_content(subgroup, embedding, record, trivial)?;
        if value > 0 {
            parts.push(format!("{}={}", record.ml, value));
        }
    }
    Ok(parts.join(","))
}

struct Candidate {
    id: String,
    ordinal: usize,
    setting: [[i32; 3]; 3],
    setting_denominator: i32,
    shift: [i32; 4],
}

fn parse_candidate(line: &str, line_number: usize) -> Result<Option<Candidate>, String> {
    let trimmed = line.split('#').next().unwrap_or("").trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let fields: Vec<&str> = trimmed.split_whitespace().collect();
    if fields.len() != 16 {
        return Err(format!(
            "line {line_number}: expected 16 fields (id, ordinal, U, denominator, shift), found {}",
            fields.len()
        ));
    }
    let id = fields[0].to_string();
    let mut values = [0i32; 15];
    for (slot, field) in values.iter_mut().zip(&fields[1..]) {
        *slot = field
            .parse()
            .map_err(|_| format!("line {line_number}: {field:?} is not an integer"))?;
    }
    if values[0] < 0 {
        return Err(format!("line {line_number}: negative ordinal {}", values[0]));
    }
    let mut setting = [[0i32; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            setting[row][column] = values[1 + row * 3 + column];
        }
    }
    let mut shift = [0i32; 4];
    shift.copy_from_slice(&values[11..15]);
    if shift[3] <= 0 {
        return Err(format!(
            "line {line_number}: shift denominator {} must be positive",
            shift[3]
        ));
    }
    if values[10] <= 0 {
        return Err(format!(
            "line {line_number}: setting denominator {} must be positive",
            values[10]
        ));
    }
    Ok(Some(Candidate {
        id,
        ordinal: values[0] as usize,
        setting,
        setting_denominator: values[10],
        shift,
    }))
}

/// Every ordinary isotropy record of the pinned table, keyed by its ordinal.
fn load_subgroups() -> Result<HashMap<usize, IsotropySubgroup>, String> {
    let mut table = HashMap::new();
    for sg in 1u8..=230 {
        let records = query::irreps_of(sg);
        for record in records.iter() {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            let subgroups = isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)
                .map_err(|error| format!("space group {sg} irrep {}: {error}", record.ml))?;
            for subgroup in subgroups {
                if table.insert(subgroup.ordinal, subgroup).is_some() {
                    return Err(format!("ordinal {} appears twice", subgroup.ordinal));
                }
            }
        }
    }
    Ok(table)
}

fn main() -> Result<(), String> {
    let profile = std::env::args().any(|argument| argument == "--profile");
    let table = load_subgroups()?;
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    let mut unknown = 0usize;
    let mut rejected = 0usize;
    let mut unavailable = 0usize;
    for (index, line) in stdin.lock().lines().enumerate() {
        let line = line.map_err(|error| format!("stdin: {error}"))?;
        let Some(candidate) = parse_candidate(&line, index + 1)? else {
            continue;
        };
        let Some(subgroup) = table.get(&candidate.ordinal) else {
            writeln!(out, "{} ERR unknown_ordinal", candidate.id)
                .map_err(|error| error.to_string())?;
            unknown += 1;
            continue;
        };
        match SubgroupEmbedding::probe_embedding(
            subgroup,
            candidate.setting,
            candidate.setting_denominator,
            candidate.shift,
        ) {
            Ok(embedding) => {
                let (trivial, reason) = render_trivial(condensing_trivial(subgroup, &embedding));
                if let Some(reason) = reason {
                    eprintln!("{}: trivial content unavailable: {reason}", candidate.id);
                    unavailable += 1;
                }
                if profile {
                    let text = match trivial_profile(subgroup, &embedding) {
                        Ok(text) => text,
                        Err(reason) => {
                            eprintln!("{}: profile unavailable: {reason}", candidate.id);
                            unavailable += 1;
                            String::from("error")
                        }
                    };
                    writeln!(
                        out,
                        "{} OK {} trivial={} profile={}",
                        candidate.id,
                        embedding.representatives().len(),
                        trivial,
                        text
                    )
                    .map_err(|error| error.to_string())?;
                } else {
                    writeln!(
                        out,
                        "{} OK {} trivial={}",
                        candidate.id,
                        embedding.representatives().len(),
                        trivial
                    )
                    .map_err(|error| error.to_string())?;
                }
            }
            Err(error) => {
                writeln!(out, "{} ERR {error}", candidate.id)
                    .map_err(|error| error.to_string())?;
                rejected += 1;
            }
        }
    }
    out.flush().map_err(|error| error.to_string())?;
    eprintln!(
        "probed candidates, {rejected} rejected, {unknown} unknown ordinals, \
         {unavailable} with an unavailable trivial content"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three states stay apart: a computed zero is a result, an unavailable
    /// answer is `?`, and a failure is `error` with its reason kept.
    #[test]
    fn the_trivial_token_distinguishes_zero_unknown_and_failure() {
        assert_eq!(render_trivial(Ok(Some(0))), ("0".to_string(), None));
        assert_eq!(render_trivial(Ok(Some(3))), ("3".to_string(), None));
        assert_eq!(render_trivial(Ok(None)), ("?".to_string(), None));
        let (token, reason) = render_trivial(Err("missing child data".to_string()));
        assert_eq!(token, "error");
        assert_eq!(reason.as_deref(), Some("missing child data"));
    }

    /// A real convention reports a computed multiplicity, so the wiring is not
    /// only exercised through the failure path.
    #[test]
    fn an_accepted_convention_reports_a_computed_multiplicity() {
        let table = load_subgroups().expect("pinned records");
        let subgroup = table.get(&26).expect("ordinal 26");
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).expect("embedding");
        let (token, reason) = render_trivial(condensing_trivial(subgroup, &embedding));
        assert!(reason.is_none(), "{reason:?}");
        assert!(
            token.parse::<u32>().is_ok(),
            "a successful convention must report a number, got {token}"
        );
        assert!(trivial_profile(subgroup, &embedding).is_ok());
    }
}
