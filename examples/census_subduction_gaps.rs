//! Replay identity-only audit rows and inventory *all* folded child stars.
//!
//! Usage: census_subduction_gaps <audit.tsv> > gaps.tsv
//! A stored scalar k match means data is reachable, not that decomposition is
//! proven complete. No characters, multiplicities, or missing labels are invented.
//!
//! A star with no pinned row is reported as `constructed_target` when the engine
//! can construct the child irrep there (the Bloch phase of a trivial little
//! co-group, or a complete one-dimensional projective catalogue), and as
//! `missing_discrete_scalar_data` otherwise.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};

use cryspglib::irrep::LabelConvention;
use cryspglib::irrep::isotropy::{self, IsotropySubgroup};
use cryspglib::irrep::query;
use cryspglib::irrep::subduction::star::FoldedStar;
use cryspglib::irrep::subduction::star::decompose::{
    FullStarError, constructed_targets_available, subduce_full_star_with_embedding,
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

/// Whether the production entry point can answer a star with no pinned row.
///
/// This asks the engine's own [`constructed_targets_available`] at the star's
/// canonical first point -- the point `decompose::select_representative`
/// constructs at -- instead of re-deriving the criterion here, so the census
/// cannot drift from the engine when a batch adds a new constructed source.
fn star_is_constructed(child: u8, star: &FoldedStar, reciprocal: &Lattice) -> Result<bool> {
    let point = star
        .points()
        .first()
        .ok_or("a folded star must have at least one point")?;
    Ok(constructed_targets_available(child, point.q(), reciprocal)?)
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
    let (mut probes, mut stars, mut missing, mut constructed) = (0usize, 0usize, 0usize, 0usize);
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
                let status = if !labels.is_empty() {
                    "stored_k_reachable"
                } else if star_is_constructed(child, folded, &reciprocal)? {
                    // No pinned row, but the engine constructs the child irrep at
                    // this star: the Bloch phase of a trivial co-group (R4 batch 1)
                    // or a complete one-dimensional projective catalogue (batch 2a).
                    constructed += 1;
                    "constructed_target"
                } else {
                    first_missing.get_or_insert((*folded.points()[0].q(), folded.star_size()));
                    missing += 1;
                    "missing_discrete_scalar_data"
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
        "records={} probes={probes} stars={stars} missing_stars={missing} \
         constructed_stars={constructed} reachable_stars={} replay_errors=0 dimension_errors=0",
        requests.len(),
        stars - missing - constructed
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Count the per-star statuses exactly as the manifest does.
    fn status_counts(ordinal: usize, ml: &str) -> (usize, usize, usize) {
        let contexts = subgroups().unwrap();
        let subgroup = &contexts[&ordinal];
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).unwrap();
        let child = embedding.subgroup_sg();
        let probe = query::irreps_of(subgroup.parent_sg)
            .iter()
            .find(|r| r.ml == ml)
            .unwrap();
        let star = ScalarStar::new(probe).unwrap();
        let folded = star.folded_stars(&embedding).unwrap();
        let reciprocal = child_reciprocal(child).unwrap();
        let stored = scalar_k_points(child).unwrap();
        let (mut reachable, mut constructed, mut missing) = (0, 0, 0);
        for block in &folded {
            if !matching_labels(block, &reciprocal, &stored)
                .unwrap()
                .is_empty()
            {
                reachable += 1;
            } else if star_is_constructed(child, block, &reciprocal).unwrap() {
                constructed += 1;
            } else {
                missing += 1;
            }
        }
        let total = star.dimension();
        assert_eq!(
            folded
                .iter()
                .map(FoldedStar::block_dimension)
                .sum::<u32>(),
            total
        );
        (reachable, constructed, missing)
    }

    /// R4 batch 1: the three folded stars of ordinal 13345 `W1` have no pinned
    /// row anywhere, but every one of them has a trivial little co-group, so the
    /// census reports them as constructed instead of missing and the engine
    /// answers the probe completely.
    #[test]
    fn a_trivial_co_group_star_is_constructed_not_missing() {
        let contexts = subgroups().unwrap();
        let subgroup = &contexts[&13345];
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).unwrap();
        let probe = query::irreps_of(225).iter().find(|r| r.ml == "W1").unwrap();
        assert!(
            subduce_full_star_with_embedding(subgroup, &embedding, probe).is_ok(),
            "batch 1 answers every folded star of this context"
        );
        assert_eq!(status_counts(13345, "W1"), (0, 3, 0));
    }

    /// Ordinal 13346 `W1` folds onto five stars without pinned rows, four of
    /// them on the parametric `B` line with a two-fold little co-group whose
    /// cocycle is a coboundary.  R4 batch 2a answers all five with the
    /// one-dimensional projective catalogue, so the probe is complete.
    #[test]
    fn a_parametric_star_is_answered_by_the_one_dimensional_catalogue() {
        let contexts = subgroups().unwrap();
        let subgroup = &contexts[&13346];
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).unwrap();
        let probe = query::irreps_of(225).iter().find(|r| r.ml == "W1").unwrap();
        assert!(
            subduce_full_star_with_embedding(subgroup, &embedding, probe).is_ok(),
            "the one-dimensional catalogue answers every folded star here"
        );
        assert_eq!(status_counts(13346, "W1"), (0, 5, 0));
    }

    /// The batch boundary after 2a: ordinal 3988 folds onto child #43 stars
    /// whose four-element little co-group has a **single** omega-regular class,
    /// so its irreps are two-dimensional.  That is the higher-dimensional batch,
    /// the star stays `missing_discrete_scalar_data`, and the probe is answered
    /// only by the identity-only entry point.
    #[test]
    fn a_two_dimensional_co_group_keeps_its_missing_status() {
        let contexts = subgroups().unwrap();
        let subgroup = &contexts[&3988];
        let embedding = SubgroupEmbedding::from_isotropy_subgroup(subgroup).unwrap();
        let probe = query::irreps_of(109).iter().find(|r| r.ml == "P1").unwrap();
        let (q, points) = replay_missing(subgroup, &embedding, probe).unwrap();
        assert_eq!(format_q(&q), "0,1,1/2");
        assert_eq!(points, 1);
        // Both folded stars of this probe need the higher-dimensional batch.
        assert_eq!(status_counts(3988, "P1"), (0, 0, 2));
    }
}
