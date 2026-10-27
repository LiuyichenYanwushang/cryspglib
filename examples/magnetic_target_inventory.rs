//! Catalog audit, not a source/record decomposition coverage claim.
//! Run with an optional output filename prefix; pair with audit_magnetic_target_inventory.py.
use cryspglib::irrep::{
    magnetic_summary::{SourceIrrepSummary, magnetic_irrep_summary_by_uni_partial},
    query,
    types::{CharacterRow, CompoundSelectedArmCharacter, IrrepRecord},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufWriter, Write},
};
fn clean(s: impl ToString) -> String {
    s.to_string().replace(['\t', '\n', '\r'], " ")
}
fn output(prefix: &str, suffix: &str, header: &str) -> BufWriter<File> {
    let mut w = BufWriter::new(File::create(format!("{prefix}{suffix}.tsv")).unwrap());
    writeln!(w, "{header}").unwrap();
    w
}
fn raw(owner: &IrrepRecord, row: Option<&CharacterRow>) -> String {
    let m = owner.matrices();
    let n = owner.characters().len();
    let d = owner.dim as usize;
    let shape = m.len() == n * d * d && !m.is_empty();
    let finite = m.iter().all(|x| x.is_finite());
    let trace = row.map(|r| {
        shape
            && r.len() == n
            && r.dimension() == d
            && (0..n).all(|i| {
                let z = r.values()[i];
                z.im.abs() < 1e-9
                    && (z.re - (0..d).map(|j| m[i * d * d + j * d + j]).sum::<f64>()).abs() < 1e-9
            })
    });
    format!(
        "raw_elements={},raw_ops={n},raw_dim={d},raw_shape={shape},finite={finite},typed_trace={trace:?}",
        m.len()
    )
}
fn row_info(r: &CharacterRow) -> String {
    format!(
        "typed_dim={},typed_ops={},complex={}",
        r.dimension(),
        r.len(),
        r.values().iter().any(|z| z.im.abs() > 1e-9)
    )
}
fn source_key(s: &SourceIrrepSummary) -> String {
    format!("H{}:{}:d{}:spin{}", s.sg, s.ml, s.dim, u8::from(s.spinor))
}
fn inspect(s: &SourceIrrepSummary) -> (usize, usize, usize, String, String) {
    let mut owners = 0;
    let mut authoritative = 0;
    let mut standalone = 0;
    let mut kinds = Vec::new();
    let mut details = Vec::new();
    for o in query::irreps_of(s.sg) {
        let meta = o.compound_metadata();
        if o.spinor != s.spinor
            || !(o.ml == s.ml || meta.is_some_and(|m| m.cir_labels.contains(&s.ml)))
        {
            continue;
        }
        owners += 1;
        let id = format!(
            "owner={}#{}[{:?}]",
            o.ml,
            o.id().index(),
            o.source_identity()
        );
        if s.spinor {
            match o.spinor_selected_arm_view() {
                Ok(r) => {
                    authoritative += 1;
                    standalone += 1;
                    kinds.push("spinor-character".to_string());
                    details.push(format!(
                        "{id}:typed_dim={},typed_ops={},complex={};{}",
                        r.dimension(),
                        r.len(),
                        r.values().iter().any(|z| z.im.abs() > 1e-9),
                        raw(o, None)
                    ));
                }
                Err(e) => details.push(format!("{id}:ERROR={e:?}")),
            };
        } else if meta.is_some() {
            match o.compound_selected_arm_view() {
                Ok(view) => {
                    let block = row_info(view.block_trace());
                    let raw_info = raw(o, Some(view.block_trace()));
                    let (sem, parts) = match view {
                        CompoundSelectedArmCharacter::ConjugateRealification { seed, .. } => {
                            ("realification", vec![seed])
                        }
                        CompoundSelectedArmCharacter::DistinctComponentSum {
                            first,
                            second,
                            ..
                        } => ("distinct", vec![first, second]),
                    };
                    let mut matched = false;
                    for p in parts {
                        if p.label == s.ml {
                            matched = true;
                            authoritative += 1;
                            kinds.push(format!("constituent-{sem}"));
                            details.push(format!(
                                "{id}:component={},irnumber={}:{};owner_block={block};{raw_info}",
                                p.component,
                                p.irnumber,
                                row_info(&p.row)
                            ));
                        }
                    }
                    if !matched {
                        details.push(format!(
                            "{id}:NO_EXPOSED_MATCH;sem={sem};owner_block={block};{raw_info}"
                        ));
                    }
                }
                Err(e) => details.push(format!("{id}:ERROR={e:?}")),
            };
        } else {
            match o.ordinary_scalar_selected_arm_block_trace() {
                Ok(r) => {
                    authoritative += 1;
                    standalone += 1;
                    kinds.push(
                        if s.dim == 1 {
                            "standalone-scalar"
                        } else {
                            "standalone-higher"
                        }
                        .to_string(),
                    );
                    details.push(format!("{id}:{};{}", row_info(&r), raw(o, Some(&r))));
                }
                Err(e) => details.push(format!("{id}:ERROR={e:?}")),
            };
        }
    }
    let status = if authoritative == 0 {
        "missing".to_string()
    } else if authoritative > 1 {
        "ambiguous".to_string()
    } else {
        kinds[0].clone()
    };
    (
        owners,
        authoritative,
        standalone,
        status,
        details.join(" | "),
    )
}
fn main() {
    let prefix = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/magnetic-target-".into());
    let mut uniout = output(
        &prefix,
        "uni",
        "uni\tbns\tmagnetic_type\tparent_sg\th_sg\tstatus\tgamma_points\tgamma_targets\tother_targets\tgamma_unresolved\tother_unresolved\tisotropy_records\tisotropy_witness\tdetail",
    );
    let mut targets = output(
        &prefix,
        "gamma",
        "uni\tindex\tbns\tmagnetic_type\th_sg\tanti_order\tlabel\traw_type\tfamily\tspinor\tdimension\th_source_count\th_dimensions\th_keys\tbackend\tcompleteness\tunitary_norm\tisotropy_candidates\tordinary_candidates\tmagnetic_candidates",
    );
    let mut sources = output(
        &prefix,
        "sources",
        "key\th_sg\tlabel\tdimension\tspinor\tnamed_owners\tauthoritative_rows\tstandalone_rows\tstatus\tdetail",
    );
    let mut unresolved = output(
        &prefix,
        "unresolved",
        "uni\tbns\tmagnetic_type\th_sg\tk_label\tis_gamma\tsource_irrep\tspinor\tminimum_dimension\tclassified_type\tbackend\tclassified_dimension\treason",
    );
    let mut refs = BTreeMap::<usize, Vec<String>>::new();
    for sg in 1..=230 {
        for ir in query::irreps_of(sg).iter().filter(|r| !r.spinor) {
            for (i, r) in ir.magnetic_subgroups().iter().enumerate() {
                refs.entry(r.mag_sg)
                    .or_default()
                    .push(format!("SG{sg}:{}:row{i}:{}", ir.ml, r.direction));
            }
        }
    }
    let mut seen = BTreeSet::new();
    let mut total = 0;
    let mut fatal = 0;
    for uni in 1..=1651 {
        let rs = refs.get(&uni);
        let rc = rs.map_or(0, Vec::len);
        let rw = rs.and_then(|r| r.first()).cloned().unwrap_or_default();
        match magnetic_irrep_summary_by_uni_partial(uni) {
            Err(e) => {
                fatal += 1;
                writeln!(
                    uniout,
                    "{uni}\t\t\t\t\tFATAL\t0\t0\t0\t0\t0\t{rc}\t{rw}\t{}",
                    clean(format!("{e:?}"))
                )
                .unwrap();
            }
            Ok(p) => {
                let s = &p.summary;
                let mt = format!("{:?}", s.magnetic_type);
                let gammas = s
                    .kpoints
                    .iter()
                    .filter(|k| k.coords.0 == 0 && k.coords.1 == 0 && k.coords.2 == 0)
                    .collect::<Vec<_>>();
                let glabels = gammas
                    .iter()
                    .map(|g| g.label.as_str())
                    .collect::<BTreeSet<_>>();
                let gu = p
                    .unresolved_coreps
                    .iter()
                    .filter(|u| glabels.contains(u.k_label.as_str()))
                    .count();
                let nt = gammas.iter().map(|k| k.coreps.len()).sum::<usize>();
                total += nt;
                writeln!(
                    uniout,
                    "{uni}\t{}\t{mt}\t{}\t{}\tOK\t{}\t{nt}\t{}\t{gu}\t{}\t{rc}\t{rw}\t",
                    s.bns_label,
                    s.parent_sg,
                    s.unitary_sg,
                    gammas.len(),
                    s.kpoints.iter().map(|k| k.coreps.len()).sum::<usize>() - nt,
                    p.unresolved_coreps.len() - gu
                )
                .unwrap();
                for u in &p.unresolved_coreps {
                    writeln!(
                        unresolved,
                        "{uni}\t{}\t{mt}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:?}\t{:?}\t{:?}\t{}",
                        s.bns_label,
                        u.sg,
                        u.k_label,
                        glabels.contains(u.k_label.as_str()),
                        u.source_irrep,
                        u.spinor,
                        u.minimum_dimension,
                        u.classified_type,
                        u.wigner_source,
                        u.classified_dimension,
                        clean(&u.reason)
                    )
                    .unwrap();
                }
                for g in gammas {
                    for (index, c) in g.coreps.iter().enumerate() {
                        assert!(!c.source_irreps.is_empty());
                        let spin = c.source_irreps[0].spinor;
                        assert!(c.source_irreps.iter().all(|x| x.spinor == spin));
                        let family = if g.antiunitary_order == 0 {
                            if c.source_irreps.len() > 1 {
                                "unitary-packet".into()
                            } else {
                                "unitary-single".into()
                            }
                        } else {
                            format!("magnetic-{:?}", c.corep_type)
                        };
                        let mut norm = 0.;
                        for (i, op) in g.operations.iter().enumerate() {
                            if !op.time_reversal {
                                norm += c.characters[i].expect("unitary missing").norm_sqr();
                            }
                        }
                        norm /= g.unitary_order as f64;
                        for h in &c.source_irreps {
                            let key = source_key(h);
                            if seen.insert(key.clone()) {
                                let (o, a, st, status, d) = inspect(h);
                                writeln!(
                                    sources,
                                    "{key}\t{}\t{}\t{}\t{}\t{o}\t{a}\t{st}\t{status}\t{}",
                                    h.sg,
                                    h.ml,
                                    h.dim,
                                    h.spinor,
                                    clean(d)
                                )
                                .unwrap();
                            }
                        }
                        writeln!(targets,"{uni}\t{index}\t{}\t{mt}\t{}\t{}\t{}\t{:?}\t{family}\t{spin}\t{}\t{}\t{}\t{}\t{:?}\t{:?}\t{norm:.9}\t{}\t{}\t{}",s.bns_label,s.unitary_sg,g.antiunitary_order,clean(&c.label),c.corep_type,c.dim,c.source_irreps.len(),c.source_irreps.iter().map(|h|h.dim.to_string()).collect::<Vec<_>>().join(","),c.source_irreps.iter().map(source_key).collect::<Vec<_>>().join(";"),c.source,c.completeness,c.isotropy_candidates.len(),c.isotropy_candidates.iter().map(|i|i.ordinary.len()).sum::<usize>(),c.isotropy_candidates.iter().map(|i|i.magnetic.len()).sum::<usize>()).unwrap();
                    }
                }
            }
        }
        if uni % 100 == 0 {
            eprintln!("UNI {uni}/1651 gamma_targets={total} fatal={fatal}");
        }
        uniout.flush().unwrap();
        targets.flush().unwrap();
        sources.flush().unwrap();
        unresolved.flush().unwrap();
    }
    println!(
        "uni=1651 gamma_targets={total} unique_H_sources={} fatal={fatal} isotropy_records={} isotropy_UNIs={}",
        seen.len(),
        refs.values().map(Vec::len).sum::<usize>(),
        refs.len()
    );
}
