use cryspglib::irrep::magnetic_embedding::{embed_in_parent_conventional, geometry_of, magnetic_operations};
use cryspglib::irrep::query;

fn main() {
    let stride: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(20);
    let witnesses: [(u8, usize); 6] = [(167, 1333), (142, 1221), (88, 741), (227, 1630), (3, 24), (1, 1)];
    let mut histogram: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    let mut realising: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    let mut unique_realising = 0usize;
    let mut rows = 0usize;
    let mut witness_counts: std::collections::BTreeMap<(u8, usize), (usize, usize)> =
        std::collections::BTreeMap::new();
    for sg in 1..=230u8 {
        for (index, row) in query::magnetic_isotropy_subgroups_of(sg).into_iter().enumerate() {
            let record = row.subgroup;
            let is_witness = witnesses.iter().any(|(w, u)| *w == sg && *u == record.mag_sg);
            if !is_witness && (usize::from(sg) + index) % stride != 0 {
                continue;
            }
            rows += 1;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("ops");
            let embeddings = embed_in_parent_conventional(&geometry, &set).expect("search");
            *histogram.entry(embeddings.len()).or_insert(0) += 1;
            let good = embeddings.iter().filter(|e| e.realises_record_lattice).count();
            *realising.entry(good).or_insert(0) += 1;
            if good == 1 {
                unique_realising += 1;
            }
            if is_witness {
                witness_counts.entry((sg, record.mag_sg)).or_insert((embeddings.len(), good));
            }
        }
    }
    println!("rows {rows}");
    println!("zero {}/{} realising_zero {}/{} unique {unique_realising}",
        histogram.get(&0).copied().unwrap_or(0), rows,
        realising.get(&0).copied().unwrap_or(0), rows);
    println!("witnesses {witness_counts:?}");
}
