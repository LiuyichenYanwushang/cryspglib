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
//! Output is one line per candidate: `id OK <representatives>` when the engine
//! accepts the convention, otherwise `id ERR <error>`.  Exit status is 0 when
//! every candidate was evaluated and 1 on malformed input.
//!
//! ```text
//! cargo run --release -p cryspglib --example probe_subduction_settings < candidates.tsv
//! ```

use std::collections::HashMap;
use std::io::{self, BufRead, Write};

use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::{LabelConvention, subduction::SubgroupEmbedding};

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
    let table = load_subgroups()?;
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    let mut unknown = 0usize;
    let mut rejected = 0usize;
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
                writeln!(out, "{} OK {}", candidate.id, embedding.representatives().len())
                    .map_err(|error| error.to_string())?;
            }
            Err(error) => {
                writeln!(out, "{} ERR {error}", candidate.id)
                    .map_err(|error| error.to_string())?;
                rejected += 1;
            }
        }
    }
    out.flush().map_err(|error| error.to_string())?;
    eprintln!("probed candidates, {rejected} rejected, {unknown} unknown ordinals");
    Ok(())
}
