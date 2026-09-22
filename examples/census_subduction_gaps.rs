//! Replay identity-only audit rows and inventory *all* folded child stars.
//!
//! Usage: census_subduction_gaps <audit.tsv> > gaps.tsv
//! A stored scalar k match means data is reachable, not that decomposition is
//! proven complete. No characters, multiplicities, or missing labels are invented.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::FoldedStar;
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, subduce_full_star_with_embedding,
};
use cryspglib::irrep::subduction::star::scalar_star::ScalarStar;
use cryspglib::irrep::subduction::{Lattice, Mat3R, Rat, SubgroupEmbedding, Vec3R};
use cryspglib::irrep::types::{CompoundCharacterSemantics, IrrepRecord};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
type Requests = BTreeMap<usize, BTreeSet<(u8, u8, String)>>;

fn read_requests(input: impl BufRead) -> Result<Requests> {
    let mut lines = input.lines();
    let header = lines.next().ok_or("empty audit")??;
    if header
        != "kind\tordinal\tparent_sg\tsubgroup_sg\tdirection\tdomain\tarms\tsize\tprobe_ml\tprobe_sg\tprobe_k\tprobe_src\tstored\tcomputed\tstatus\tdetail"
    {
        return Err("unexpected audit header".into());
    }
    let mut requests = Requests::new();
    for (line_number, line) in lines.enumerate() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 16 {
            return Err(format!("audit line {}: expected 16 fields", line_number + 2).into());
        }
        if !fields[15].contains("identity_only ") {
            continue;
        }
        if !matches!(fields[0], "identity" | "absent") || fields[2] != fields[9] {
            return Err(format!("unexpected identity-only row: {line}").into());
        }
        let ordinal = fields[1].parse()?;
        let row = (fields[2].parse()?, fields[3].parse()?, fields[8].to_owned());
        let rows = requests.entry(ordinal).or_default();
        if rows.iter().any(|(_, _, ml)| ml == &row.2) {
            return Err(format!("duplicate probe {ordinal} {}", row.2).into());
        }
        rows.insert(row);
    }
    if requests.is_empty() {
        return Err("no identity-only probes in audit".into());
    }
    Ok(requests)
}

fn subgroups() -> Result<BTreeMap<usize, IsotropySubgroup>> {
    let mut out = BTreeMap::new();
    for sg in 1..=230 {
        for record in query::irreps_of(sg) {
            if record.spinor || record.subgroups().is_empty() {
                continue;
            }
            for subgroup in isotropy::isotropy_subgroups(sg, record.ml, LabelConvention::Cdml)? {
                let ordinal = subgroup.ordinal;
                if out.insert(ordinal, subgroup).is_some() {
                    return Err(format!("duplicate subgroup ordinal {ordinal}").into());
                }
            }
        }
    }
    Ok(out)
}

fn child_reciprocal(sg: u8) -> Result<Lattice> {
    let mut rows = [[Rat::ZERO; 3]; 3];
    for (i, row) in isotropy::parent_primitive_basis(sg)?.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            rows[i][j] = Rat::from_grid(*value, 12)?;
        }
    }
    Ok(Lattice::new(Mat3R::new(rows))?.reciprocal()?)
}

/// Mirror only effective-k reachability in decompose::child_components_at.
/// Conjugate realifications also provide -k; distinct sums use the stored k.
fn scalar_k_points(sg: u8) -> Result<Vec<(&'static str, Vec3R)>> {
    let mut out = Vec::new();
    for record in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
        let k = record.k_vector();
        let mut values = [Rat::ZERO; 3];
        for (value, numerator) in values.iter_mut().zip(k.numerators) {
            *value = Rat::new(i128::from(numerator), i128::from(k.denominator))?;
        }
        let k = Vec3R::new(values);
        out.push((record.ml, k));
        if record.compound_character_semantics()
            == Some(CompoundCharacterSemantics::ConjugateRealification)
        {
            out.push((record.ml, k.checked_neg()?));
        }
    }
    Ok(out)
}

fn matching_labels(
    star: &FoldedStar,
    reciprocal: &Lattice,
    stored: &[(&'static str, Vec3R)],
) -> Result<BTreeSet<&'static str>> {
    let mut labels = BTreeSet::new();
    for point in star.points() {
        for (label, k) in stored {
            if reciprocal.same_mod(point.q(), k)? {
                labels.insert(*label);
            }
        }
    }
    Ok(labels)
}

fn format_q(q: &Vec3R) -> String {
    (0..3)
        .map(|i| q.get(i).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn replay_missing(
    subgroup: &IsotropySubgroup,
    embedding: &SubgroupEmbedding,
    probe: &'static IrrepRecord,
) -> Result<(Vec3R, usize)> {
    match subduce_full_star_with_embedding(subgroup, embedding, probe) {
        Err(FullStarError::MissingChildStarData { sg, q, points })
            if usize::from(sg) == subgroup.record.sg =>
        {
            Ok((Vec3R::new(q), points))
        }
        other => Err(format!(
            "ordinal {} {} no longer fails with missing child data: {other:?}",
            subgroup.ordinal, probe.ml
        )
        .into()),
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 1 {
        return Err("usage: census_subduction_gaps <audit.tsv> > gaps.tsv".into());
    }
    let requests = read_requests(std::io::BufReader::new(std::fs::File::open(&args[0])?))?;
    let contexts = subgroups()?;
    let mut cache = BTreeMap::new();
    let (mut probes, mut stars, mut missing) = (0usize, 0usize, 0usize);
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "ordinal\tparent_sg\tcondensing_cdml\tdirection\tprobe_cdml\tchild_sg\tstar_index\tstar_size\tarm_count\tblock_dimension\tparent_dimension\tgamma\tstatus\tmatching_cdml\tq\tcanonical_q\tsetting_numerator\tsetting_denominator\tchild_shift"
    )?;
    for (ordinal, rows) in &requests {
        let subgroup = contexts.get(ordinal).ok_or("unknown ordinal")?;
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup)?;
        let child = embedding.subgroup_sg();
        let reciprocal = child_reciprocal(child)?;
        let stored = scalar_k_points(child)?;
        for (parent, requested_child, ml) in rows {
            if *parent != subgroup.parent_sg || *requested_child != child {
                return Err(format!("stale audit context for ordinal {ordinal}").into());
            }
            let probe = query::irreps_of(*parent)
                .iter()
                .find(|r| r.ml == ml)
                .ok_or("unknown probe")?;
            let (first_q, first_points) = replay_missing(subgroup, &embedding, probe)?;
            if let std::collections::btree_map::Entry::Vacant(entry) =
                cache.entry((*parent, probe.ml))
            {
                entry.insert(ScalarStar::new(probe)?);
            }
            let star = &cache[&(*parent, probe.ml)];
            let folded = star.folded_stars(&embedding)?;
            if folded.iter().map(FoldedStar::block_dimension).sum::<u32>() != star.dimension() {
                return Err(format!("block dimensions do not sum for {ordinal} {ml}").into());
            }
            let mut first_missing = None;
            for (index, folded) in folded.iter().enumerate() {
                let labels = matching_labels(folded, &reciprocal, &stored)?;
                let mut raw = Vec::new();
                let mut reduced = Vec::new();
                let mut gamma = false;
                for point in folded.points() {
                    raw.push(format_q(point.q()));
                    reduced.push(format_q(&reciprocal.reduce(point.q())?.representative));
                    gamma |= reciprocal.contains(point.q())?;
                }
                reduced.sort();
                if folded.block_dimension() == 0 || (labels.is_empty() && gamma) {
                    return Err(
                        format!("unexpected empty block or missing Gamma: {ordinal} {ml}").into(),
                    );
                }
                let status = if labels.is_empty() {
                    first_missing.get_or_insert((*folded.points()[0].q(), folded.star_size()));
                    missing += 1;
                    "missing_discrete_scalar_data"
                } else {
                    "stored_k_reachable"
                };
                writeln!(
                    out,
                    "{ordinal}\t{parent}\t{}\t{}\t{ml}\t{child}\t{index}\t{}\t{}\t{}\t{}\t{gamma}\t{status}\t{}\t{}\t{}\t{:?}\t{}\t{}",
                    subgroup.irrep_ml,
                    subgroup.record.direction_label,
                    folded.star_size(),
                    folded.arm_count(),
                    folded.block_dimension(),
                    star.dimension(),
                    labels.into_iter().collect::<Vec<_>>().join(";"),
                    raw.join(";"),
                    reduced.join(";"),
                    embedding.setting(),
                    embedding.setting_denominator(),
                    format_q(embedding.child_shift())
                )?;
                stars += 1;
            }
            if first_missing != Some((first_q, first_points)) {
                return Err(
                    format!("reachability disagrees with engine error for {ordinal} {ml}").into(),
                );
            }
            probes += 1;
        }
    }
    out.flush()?;
    eprintln!(
        "records={} probes={probes} stars={stars} missing_stars={missing} reachable_stars={} replay_errors=0 dimension_errors=0",
        requests.len(),
        stars - missing
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn w1_has_three_missing_stars_not_only_the_first_reported_error() {
        let contexts = subgroups().unwrap();
        let subgroup = &contexts[&13345];
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).unwrap();
        let probe = query::irreps_of(225).iter().find(|r| r.ml == "W1").unwrap();
        let (q, points) = replay_missing(subgroup, &embedding, probe).unwrap();
        assert_eq!(format_q(&q), "-3/4,1/4,-3/4");
        assert_eq!(points, 2);
        let star = ScalarStar::new(probe).unwrap();
        let folded = star.folded_stars(&embedding).unwrap();
        assert_eq!(folded.len(), 3);
        let reciprocal = child_reciprocal(8).unwrap();
        let stored = scalar_k_points(8).unwrap();
        for block in &folded {
            assert!(
                matching_labels(block, &reciprocal, &stored)
                    .unwrap()
                    .is_empty()
            );
            assert_eq!(block.block_dimension(), 2);
        }
        assert_eq!(star.dimension(), 6);
    }
}
