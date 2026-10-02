use cryspglib::irrep::magnetic_embedding::{geometry_of, magnetic_operations, search_parent_setting};
use cryspglib::irrep::query;

fn main() {
    let mut histogram: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    let mut unique = 0usize;
    let mut permissive = 0usize;
    let mut contained_but_zero = 0usize;
    let mut rows = 0usize;
    for sg in 1..=230u8 {
        for row in query::magnetic_isotropy_subgroups_of(sg) {
            rows += 1;
            let record = row.subgroup;
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(record.mag_sg).expect("ops");
            let report = search_parent_setting(&geometry, &set).expect("search");
            *histogram.entry(report.survivors).or_insert(0) += 1;
            if report.survivors == 1 {
                unique += 1;
            }
            let full = report.naive.unitary_in_parent == report.naive.unitary_classes
                && report.naive.antiunitary_classes_in_parent == set.antiunitary_classes;
            if full {
                permissive += 1;
                if report.survivors == 0 {
                    contained_but_zero += 1;
                }
            }
        }
    }
    println!("rows {rows} unique {unique} permissive {permissive} contained_but_zero {contained_but_zero}");
    let mut pairs: Vec<(usize, usize)> = histogram.into_iter().collect();
    pairs.sort();
    println!("histogram {pairs:?}");
    for (sg, uni) in [(167u8, 1333usize), (142, 1221), (88, 741), (227, 1630), (3, 24), (1, 1)] {
        if let Some(record) = query::magnetic_isotropy_subgroups_of(sg)
            .into_iter()
            .find(|row| row.subgroup.mag_sg == uni)
            .map(|row| row.subgroup)
        {
            let geometry = geometry_of(sg, &record).expect("geometry");
            let set = magnetic_operations(uni).expect("ops");
            let report = search_parent_setting(&geometry, &set).expect("search");
            println!("witness ({sg},{uni}) survivors {}", report.survivors);
        }
    }
}
