//! Stored-record/typed-view inventory; successful construction is not subduction certification.
use cryspglib::irrep::subduction::star::scalar_star::ScalarStar;
use cryspglib::irrep::{
    query, types::CompoundSelectedArmCharacter, w_little_characters_data::W_LITTLE_CHARACTERS,
};
use std::{
    collections::BTreeMap,
    io::{self, Write},
};
fn main() {
    let mut out = io::BufWriter::new(io::stdout().lock());
    writeln!(out, "id\tsg\tml\tspinor\tdim\tkx\tky\tkz\tkd\tis_point\tfamily\tselected_dim\tcharacter_count\tmatrix_count\tview\tstar_components\tstar_arms\tstar_dim\tiso_count\tmag_iso_count\tsource_identity\tcir_sources").unwrap();
    let mut stats = BTreeMap::<String, usize>::new();
    for sg in 1..=230 {
        for r in query::irreps_of(sg) {
            let family = if r.spinor {
                "spinor"
            } else if r.compound_metadata().is_some() {
                "compound"
            } else {
                "ordinary"
            };
            let view = if r.spinor {
                r.spinor_selected_arm_view()
                    .map(|v| (v.dimension(), v.len()))
            } else if family == "compound" {
                r.compound_selected_arm_view()
                    .map(|v| (v.block_trace().dimension(), v.block_trace().len()))
            } else {
                r.ordinary_scalar_selected_arm_block_trace()
                    .map(|v| (v.dimension(), v.len()))
            };
            let (sd, status) = match view {
                Ok((d, _)) => (d, "ok".to_string()),
                Err(e) => (0, format!("{e:?}")),
            };
            let (comps, arms, dim) = if !r.spinor {
                match ScalarStar::new(r) {
                    Ok(s) => (s.components().len(), s.arm_count(), s.dimension()),
                    Err(e) => {
                        eprintln!("STAR_ERROR SG{} {} {:?}", sg, r.ml, e);
                        (0, 0, 0)
                    }
                }
            } else {
                (0, 0, 0)
            };
            let (iso, mag) = if r.spinor {
                (0, 0)
            } else {
                (r.subgroups().len(), r.magnetic_subgroups().len())
            };
            let region = if [r.kx, r.ky, r.kz] == [0; 3] {
                "gamma"
            } else {
                "non_gamma"
            };
            *stats.entry(format!("{family}/{region}")).or_default() += 1;
            *stats.entry(format!("{family}/view_{status}")).or_default() += 1;
            if !r.spinor {
                *stats
                    .entry(format!("{family}/star_components_{comps}"))
                    .or_default() += 1;
                *stats
                    .entry(format!("{family}/star_arms_{arms}"))
                    .or_default() += 1;
            }
            if r.spinor && r.spin_character_imag().iter().any(|v| v.abs() > 1e-10) {
                *stats
                    .entry("spinor/nonreal_character".to_string())
                    .or_default() += 1;
            }
            if let Ok(v) = r.compound_selected_arm_view() {
                *stats
                    .entry(
                        match v {
                            CompoundSelectedArmCharacter::ConjugateRealification { .. } => {
                                "compound/conjugate_realification"
                            }
                            CompoundSelectedArmCharacter::DistinctComponentSum { .. } => {
                                "compound/distinct_component_sum"
                            }
                        }
                        .to_string(),
                    )
                    .or_default() += 1;
            }
            let cir_sources = if let Some(m) = r.compound_metadata() {
                format!("{},{}", m.cir_irnumbers[0], m.cir_irnumbers[1])
            } else if let cryspglib::irrep::types::IrrepSourceIdentity::OrdinaryScalar {
                cir_irnumber,
            } = r.source_identity()
            {
                cir_irnumber.to_string()
            } else {
                String::new()
            };
            writeln!(out,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}",r.id().index(),sg,r.ml,r.spinor,r.dim,r.kx,r.ky,r.kz,r.kd,r.is_point(),family,sd,r.characters().len(),r.matrices().len(),status,comps,arms,dim,iso,mag,r.source_identity(),cir_sources).unwrap();
        }
    }
    eprintln!("TOTAL {}", query::total_irreps());
    for (k, v) in stats {
        eprintln!("{k}\t{v}");
    }
    eprintln!("W_LINES {}", W_LITTLE_CHARACTERS.len());
    for t in W_LITTLE_CHARACTERS {
        eprintln!(
            "W\t{}\t{}\t{}\t{:?}\t{}\t{}",
            t.space_group,
            t.label,
            t.k_label,
            t.direction,
            t.dimension,
            t.operations.len()
        );
    }
    let mut n = 0;
    let mut zero = vec![];
    for sg in 1..=230 {
        let (r, _, _) = cryspglib::irrep::types::IrrepRecord::spin_ops_for_sg(sg);
        n += r.len() / 9;
        if r.is_empty() {
            zero.push(sg);
        }
    }
    eprintln!("SPIN_SG_OPERATIONS {} missing {:?}", n, zero);
}
