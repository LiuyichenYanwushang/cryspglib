//! Emit a canonical binary snapshot of the magnetic database/runtime APIs.
//!
//! This example is a parity oracle for the typed Python provenance loader.  It
//! intentionally obtains database values from the public raw tables and
//! public runtime accessors; it does not read the committed JSON artifact.

use cryspglib::{MagneticType, msg_database, spg_database};
use std::io::{self, Write};

const MAGIC: &[u8; 8] = b"CRYMDBP\0";
const VERSION: u32 = 1;
const SECTION_COUNT: u32 = 13;
const SPG_HALL_COUNT: usize = 531;
const SPG_HALL_SETTINGS: usize = 530;
const SPG_OPERATION_COUNT: usize = 8147;
const SPG_STANDARD_OPERATION_END: usize = 7389;
const MSG_UNI_COUNT: usize = 1652;
const MSG_HALL_SLOTS: usize = 18;
const MSG_OPERATION_COUNT: usize = 76683;
const MSG_ACTIVE_SPAN_COUNT: usize = 4479;
const ALT_VALUE_COUNT: usize = 536;
const ROTATION_PAYLOAD: i32 = 3_i32.pow(9);
const TRANSLATION_DENOMINATOR: f64 = 12.0;
const SPACE_OPERATION_SCALE: i32 = ROTATION_PAYLOAD * 12_i32.pow(3);

type Result<T> = std::result::Result<T, String>;

fn error<T>(message: impl Into<String>) -> Result<T> {
    Err(message.into())
}

fn append_i32(payload: &mut Vec<u8>, value: i32) {
    payload.extend_from_slice(&value.to_le_bytes());
}

fn append_u8(payload: &mut Vec<u8>, value: u8) {
    payload.push(value);
}

fn append_u16(payload: &mut Vec<u8>, value: u16) {
    payload.extend_from_slice(&value.to_le_bytes());
}

fn append_u32(payload: &mut Vec<u8>, value: u32) {
    payload.extend_from_slice(&value.to_le_bytes());
}

fn append_u64(payload: &mut Vec<u8>, value: u64) {
    payload.extend_from_slice(&value.to_le_bytes());
}

fn append_string(payload: &mut Vec<u8>, value: &str, label: &str) -> Result<()> {
    if value.as_bytes().contains(&0) {
        return error(format!("{label} contains NUL"));
    }
    let length = u32::try_from(value.len())
        .map_err(|_| format!("{label} is too long for u32"))?;
    append_u32(payload, length);
    payload.extend_from_slice(value.as_bytes());
    Ok(())
}

fn checked_u16(value: usize, label: &str) -> Result<u16> {
    u16::try_from(value).map_err(|_| format!("{label} does not fit u16"))
}

fn checked_u32(value: usize, label: &str) -> Result<u32> {
    u32::try_from(value).map_err(|_| format!("{label} does not fit u32"))
}

fn exact_translation(value: f64, label: &str) -> Result<u8> {
    if !value.is_finite() {
        return error(format!("{label} is not finite"));
    }
    if value == 0.0 && value.to_bits() == (-0.0_f64).to_bits() {
        return error(format!("{label} is negative zero"));
    }
    let scaled = value * TRANSLATION_DENOMINATOR;
    let rounded = scaled.round();
    if scaled.to_bits() != rounded.to_bits() {
        return error(format!("{label} is not an exact denominator-12 value"));
    }
    if !(0.0..=11.0).contains(&rounded) {
        return error(format!("{label} numerator is out of range"));
    }
    let numerator = rounded as u8;
    let canonical = f64::from(numerator) / TRANSLATION_DENOMINATOR;
    if value.to_bits() != canonical.to_bits() {
        return error(format!("{label} is not canonical denominator-12 f64"));
    }
    Ok(numerator)
}

fn append_rotation(payload: &mut Vec<u8>, rotation: [[i32; 3]; 3], label: &str) -> Result<()> {
    for (row_index, row) in rotation.into_iter().enumerate() {
        for (column_index, value) in row.into_iter().enumerate() {
            if !(-1..=1).contains(&value) {
                return error(format!(
                    "{label} rotation[{row_index}][{column_index}] is outside -1..1"
                ));
            }
            let value = i8::try_from(value)
                .map_err(|_| format!("{label} rotation does not fit i8"))?;
            payload.push(value as u8);
        }
    }
    Ok(())
}

fn append_spg_operation(
    payload: &mut Vec<u8>,
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
    label: &str,
) -> Result<()> {
    append_rotation(payload, rotation, label)?;
    for (index, value) in translation.into_iter().enumerate() {
        append_u8(payload, exact_translation(value, &format!("{label} translation[{index}]"))?);
    }
    Ok(())
}

fn append_magnetic_operation(
    payload: &mut Vec<u8>,
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
    time_reversal: bool,
    label: &str,
) -> Result<()> {
    append_spg_operation(payload, rotation, translation, label)?;
    append_u8(payload, u8::from(time_reversal));
    Ok(())
}

fn append_section(frame: &mut Vec<u8>, tag: &[u8; 4], count: usize, payload: Vec<u8>) -> Result<()> {
    let count = u64::try_from(count).map_err(|_| format!("{} count overflow", String::from_utf8_lossy(tag)))?;
    let payload_length = u64::try_from(payload.len())
        .map_err(|_| format!("{} payload length overflow", String::from_utf8_lossy(tag)))?;
    frame.extend_from_slice(tag);
    append_u64(frame, count);
    append_u64(frame, payload_length);
    frame.extend_from_slice(&payload);
    Ok(())
}

fn build_sgno() -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    for hall in 0..SPG_HALL_COUNT {
        let raw = spg_database::SPACEGROUP_TYPES[hall].number;
        let api = spg_database::get_spacegroup_type(hall).number;
        if raw != api {
            return error(format!("SGNO Hall {hall} raw/API mismatch"));
        }
        append_i32(&mut payload, i32::try_from(raw).map_err(|_| format!("SGNO[{hall}] out of i32"))?);
    }
    Ok(payload)
}

fn build_sgix() -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    let sentinel = spg_database::SYMMETRY_OPERATION_INDEX[0];
    if sentinel != [0, 0] || spg_database::get_operation_index(0) != (0, 0) {
        return error("SGIX sentinel mismatch");
    }
    append_i32(&mut payload, sentinel[0]);
    append_i32(&mut payload, sentinel[1]);
    let mut previous_end = 1usize;
    for hall in 1..SPG_HALL_COUNT {
        let raw = spg_database::SYMMETRY_OPERATION_INDEX[hall];
        let api = spg_database::get_operation_index(hall);
        if raw[0] < 0 || raw[1] < 0 || (raw[0] as usize, raw[1] as usize) != api {
            return error(format!("SGIX Hall {hall} raw/API mismatch"));
        }
        let order = raw[0] as usize;
        let offset = raw[1] as usize;
        if order == 0 || offset != previous_end || offset + order > SPG_STANDARD_OPERATION_END {
            return error(format!("SGIX Hall {hall} has invalid standard span"));
        }
        previous_end = offset + order;
        append_i32(&mut payload, raw[0]);
        append_i32(&mut payload, raw[1]);
    }
    if previous_end != SPG_STANDARD_OPERATION_END {
        return error("SGIX standard boundary mismatch");
    }
    Ok(payload)
}

fn build_sgrw() -> Result<Vec<u8>> {
    if spg_database::SYMMETRY_OPERATIONS.len() != SPG_OPERATION_COUNT
        || spg_database::SYMMETRY_OPERATIONS[0] != 0
    {
        return error("SGRW census or sentinel mismatch");
    }
    let mut payload = Vec::with_capacity(SPG_OPERATION_COUNT * 4);
    for (index, &code) in spg_database::SYMMETRY_OPERATIONS.iter().enumerate() {
        if index > 0 && !(0..SPACE_OPERATION_SCALE).contains(&code) {
            return error(format!("SGRW[{index}] encoding out of range"));
        }
        append_i32(&mut payload, code);
    }
    Ok(payload)
}

fn build_mtyp() -> Result<Vec<u8>> {
    if msg_database::MAGNETIC_SPACEGROUP_TYPES.len() != MSG_UNI_COUNT {
        return error("MTYP census mismatch");
    }
    let mut payload = Vec::new();
    for uni in 0..MSG_UNI_COUNT {
        let raw = &msg_database::MAGNETIC_SPACEGROUP_TYPES[uni];
        let api = msg_database::get_magnetic_spacegroup_type(uni);
        if raw.uni_number != api.uni_number
            || raw.litvin_number != api.litvin_number
            || raw.bns_number != api.bns_number
            || raw.og_number != api.og_number
            || raw.number != api.number
            || raw.type_ != api.type_
        {
            return error(format!("MTYP UNI {uni} raw/API mismatch"));
        }
        for (value, name) in [
            (raw.uni_number, "uni"),
            (raw.litvin_number, "litvin"),
            (raw.number, "parent spacegroup"),
        ] {
            append_i32(
                &mut payload,
                i32::try_from(value).map_err(|_| format!("MTYP[{uni}] {name} out of i32"))?,
            );
        }
        let magnetic_type = raw.type_ as i32;
        if uni == 0 && magnetic_type != MagneticType::NonMagnetic as i32 {
            return error("MTYP sentinel type mismatch");
        }
        append_i32(&mut payload, magnetic_type);
        append_string(&mut payload, raw.bns_number, &format!("MTYP[{uni}].bns"))?;
        append_string(&mut payload, raw.og_number, &format!("MTYP[{uni}].og"))?;
    }
    Ok(payload)
}

fn build_muni() -> Result<Vec<u8>> {
    if msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING.len() != MSG_UNI_COUNT {
        return error("MUNI census mismatch");
    }
    let mut payload = Vec::new();
    for (uni, raw) in msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING.iter().enumerate() {
        let count = raw[0];
        let first = raw[1];
        if count < 0 || first < 0 {
            return error(format!("MUNI[{uni}] is negative"));
        }
        if uni == 0 {
            if *raw != [0, 0] {
                return error("MUNI sentinel mismatch");
            }
        } else if count == 0
            || count as usize > MSG_HALL_SLOTS
            || first == 0
            || first as usize > SPG_HALL_SETTINGS
            || first as usize + count as usize - 1 > SPG_HALL_SETTINGS
        {
            return error(format!("MUNI[{uni}] range mismatch"));
        }
        append_i32(&mut payload, count);
        append_i32(&mut payload, first);
    }
    Ok(payload)
}

fn build_mhll() -> Result<Vec<u8>> {
    if msg_database::MAGNETIC_SPACEGROUP_HALL_MAPPING.len() != SPG_HALL_COUNT {
        return error("MHLL census mismatch");
    }
    let mut payload = Vec::new();
    for (hall, raw) in msg_database::MAGNETIC_SPACEGROUP_HALL_MAPPING.iter().enumerate() {
        let api = msg_database::get_uni_candidates(hall)
            .ok_or_else(|| format!("MHLL Hall {hall} has no API mapping"))?;
        if raw[0] < 0 || raw[1] < 0 || [raw[0] as usize, raw[1] as usize] != api {
            return error(format!("MHLL Hall {hall} raw/API mismatch"));
        }
        if hall == 0 {
            if *raw != [0, 0] {
                return error("MHLL sentinel mismatch");
            }
        } else if raw[0] < 1 || raw[1] < raw[0] || raw[1] as usize >= MSG_UNI_COUNT {
            return error(format!("MHLL Hall {hall} range mismatch"));
        }
        append_i32(&mut payload, raw[0]);
        append_i32(&mut payload, raw[1]);
    }
    Ok(payload)
}

fn build_midx() -> Result<Vec<u8>> {
    if msg_database::MAGNETIC_SPACEGROUP_OPERATION_INDEX.len() != MSG_UNI_COUNT {
        return error("MIDX UNI census mismatch");
    }
    let mut payload = Vec::new();
    for (uni, row) in msg_database::MAGNETIC_SPACEGROUP_OPERATION_INDEX.iter().enumerate() {
        let active_count = if uni == 0 {
            0
        } else {
            msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni][0] as usize
        };
        for (slot, raw) in row.iter().enumerate() {
            if raw[0] < 0 || raw[1] < 0 {
                return error(format!("MIDX[{uni}][{slot}] is negative"));
            }
            let order = raw[0] as usize;
            let offset = raw[1] as usize;
            if slot < active_count {
                if order == 0 || offset == 0 || offset + order > MSG_OPERATION_COUNT {
                    return error(format!("MIDX[{uni}][{slot}] active span invalid"));
                }
            } else if *raw != [0, 0] {
                return error(format!("MIDX[{uni}][{slot}] inactive span nonzero"));
            }
            append_i32(&mut payload, raw[0]);
            append_i32(&mut payload, raw[1]);
        }
    }
    Ok(payload)
}

fn build_mraw() -> Result<Vec<u8>> {
    if msg_database::MAGNETIC_SYMMETRY_OPERATIONS.len() != MSG_OPERATION_COUNT
        || msg_database::MAGNETIC_SYMMETRY_OPERATIONS[0] != 0
    {
        return error("MRAW census or sentinel mismatch");
    }
    let mut payload = Vec::with_capacity(MSG_OPERATION_COUNT * 4);
    for (index, &code) in msg_database::MAGNETIC_SYMMETRY_OPERATIONS.iter().enumerate() {
        if index > 0 && !(0..(2 * SPACE_OPERATION_SCALE)).contains(&code) {
            return error(format!("MRAW[{index}] encoding out of range"));
        }
        append_i32(&mut payload, code);
    }
    Ok(payload)
}

fn build_malt() -> Result<Vec<u8>> {
    if msg_database::ALTERNATIVE_TRANSFORMATIONS.len() != MSG_UNI_COUNT {
        return error("MALT UNI census mismatch");
    }
    let mut payload = Vec::new();
    let mut occurrences = 0usize;
    for (uni, row) in msg_database::ALTERNATIVE_TRANSFORMATIONS.iter().enumerate() {
        let active_count = if uni == 0 {
            0
        } else {
            msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni][0] as usize
        };
        for (slot, values) in row.iter().enumerate() {
            let first_zero = values.iter().position(|&value| value == 0);
            if slot >= active_count {
                if values.iter().any(|&value| value != 0) {
                    return error(format!("MALT[{uni}][{slot}] inactive tail nonzero"));
                }
            } else {
                let first_zero = first_zero
                    .ok_or_else(|| format!("MALT[{uni}][{slot}] lacks terminator"))?;
                if first_zero > 6 || values[first_zero..].iter().any(|&value| value != 0) {
                    return error(format!("MALT[{uni}][{slot}] terminator/tail invalid"));
                }
                for (index, &value) in values[..first_zero].iter().enumerate() {
                    if !(0..SPACE_OPERATION_SCALE).contains(&value) || value == 0 {
                        return error(format!("MALT[{uni}][{slot}][{index}] encoding invalid"));
                    }
                }
                occurrences += first_zero;
            }
            for &value in values {
                append_i32(&mut payload, value);
            }
        }
    }
    if occurrences != ALT_VALUE_COUNT {
        return error(format!("MALT occurrence census {occurrences} != {ALT_VALUE_COUNT}"));
    }
    Ok(payload)
}

fn build_sdec() -> Result<Vec<u8>> {
    let mut payload = Vec::with_capacity((SPG_OPERATION_COUNT - 1) * 20);
    for index in 1..SPG_OPERATION_COUNT {
        let encoded = spg_database::SYMMETRY_OPERATIONS[index];
        if !(0..SPACE_OPERATION_SCALE).contains(&encoded) || encoded == 0 {
            return error(format!("SDEC[{index}] encoding invalid"));
        }
        let (rotation, translation) = spg_database::decode_symmetry(encoded);
        append_u32(&mut payload, checked_u32(index, &format!("SDEC[{index}] index"))?);
        append_i32(&mut payload, encoded);
        append_spg_operation(&mut payload, rotation, translation, &format!("SDEC[{index}]"))?;
    }
    Ok(payload)
}

fn build_sapi() -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    for hall in 1..=SPG_HALL_SETTINGS {
        let operations = spg_database::get_spacegroup_operations(hall)
            .ok_or_else(|| format!("SAPI Hall {hall} runtime lookup failed"))?;
        let raw = spg_database::SYMMETRY_OPERATION_INDEX[hall];
        if operations.len() != raw[0] as usize {
            return error(format!("SAPI Hall {hall} span/API order mismatch"));
        }
        if operations.rot.len() != operations.trans.len()
            || operations.rot.len() != raw[0] as usize
        {
            return error(format!("SAPI Hall {hall} runtime vector lengths mismatch"));
        }
        append_u16(&mut payload, checked_u16(hall, &format!("SAPI[{hall}] Hall"))?);
        append_u16(&mut payload, checked_u16(operations.len(), &format!("SAPI[{hall}] count"))?);
        for (index, (&rotation, &translation)) in operations
            .rot
            .iter()
            .zip(operations.trans.iter())
            .enumerate()
        {
            append_spg_operation(&mut payload, rotation, translation, &format!("SAPI[{hall}][{index}]"))?;
        }
    }
    Ok(payload)
}

fn check_magnetic_runtime_operation(
    encoded: i32,
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
    time_reversal: bool,
    label: &str,
) -> Result<()> {
    if !(0..(2 * SPACE_OPERATION_SCALE)).contains(&encoded) || encoded == 0 {
        return error(format!("{label} raw encoding invalid"));
    }
    let time = encoded / SPACE_OPERATION_SCALE;
    if time != i32::from(time_reversal) {
        return error(format!("{label} time reversal/raw mismatch"));
    }
    let (expected_rotation, expected_translation) =
        spg_database::decode_symmetry(encoded % SPACE_OPERATION_SCALE);
    if rotation != expected_rotation {
        return error(format!("{label} runtime rotation/raw mismatch"));
    }
    for (index, value) in translation.into_iter().enumerate() {
        let expected = exact_translation(
            expected_translation[index],
            &format!("{label} expected translation[{index}]"),
        )?;
        let actual = exact_translation(value, &format!("{label} translation[{index}]"))?;
        if actual != expected {
            return error(format!("{label} runtime translation/raw mismatch"));
        }
    }
    Ok(())
}

fn compare_magnetic_alias(
    uni: usize,
    first_hall: usize,
    alias: &cryspglib::symmetry::MagneticSymmetry,
    real: &cryspglib::symmetry::MagneticSymmetry,
) -> Result<()> {
    if alias.rot.len() != alias.trans.len()
        || alias.rot.len() != alias.timerev.len()
        || real.rot.len() != real.trans.len()
        || real.rot.len() != real.timerev.len()
        || alias.rot.len() != real.rot.len()
    {
        return error(format!(
            "MAPI UNI {uni} Hall=0/{first_hall} alias vector lengths mismatch"
        ));
    }
    for index in 0..alias.rot.len() {
        if alias.rot[index] != real.rot[index]
            || alias.timerev[index] != real.timerev[index]
            || alias.trans[index]
                .iter()
                .zip(real.trans[index].iter())
                .any(|(left, right)| left.to_bits() != right.to_bits())
        {
            return error(format!(
                "MAPI UNI {uni} Hall=0/{first_hall} alias differs at operation {index}"
            ));
        }
    }
    Ok(())
}

fn compare_transformation_alias(
    uni: usize,
    first_hall: usize,
    alias: &cryspglib::symmetry::Symmetry,
    real: &cryspglib::symmetry::Symmetry,
) -> Result<()> {
    if alias.rot.len() != alias.trans.len()
        || real.rot.len() != real.trans.len()
        || alias.rot.len() != real.rot.len()
    {
        return error(format!(
            "TAPI UNI {uni} Hall=0/{first_hall} alias vector lengths mismatch"
        ));
    }
    for index in 0..alias.rot.len() {
        if alias.rot[index] != real.rot[index]
            || alias.trans[index]
                .iter()
                .zip(real.trans[index].iter())
                .any(|(left, right)| left.to_bits() != right.to_bits())
        {
            return error(format!(
                "TAPI UNI {uni} Hall=0/{first_hall} alias differs at transformation {index}"
            ));
        }
    }
    Ok(())
}

fn validate_zero_alias_contract() -> Result<()> {
    if spg_database::get_operation(0).is_some()
        || spg_database::get_spacegroup_operations(0).is_some()
    {
        return error("SPG Hall=0 unexpectedly accepted as a real query");
    }
    if msg_database::get_spacegroup_operations(0, 0).is_some()
        || msg_database::get_spacegroup_operations(0, 1).is_some()
        || msg_database::get_std_transformations(0, 0).is_some()
        || msg_database::get_std_transformations(0, 1).is_some()
    {
        return error("MSG UNI=0 unexpectedly accepted as a real query");
    }
    Ok(())
}

fn validate_uni_zero_alias(uni: usize, first_hall: usize) -> Result<()> {
    if !(1..=SPG_HALL_SETTINGS).contains(&first_hall) {
        return error(format!("MAPI UNI {uni} first Hall is out of range"));
    }
    let alias_operations = msg_database::get_spacegroup_operations(uni, 0)
        .ok_or_else(|| format!("MAPI UNI {uni} Hall=0 alias lookup failed"))?;
    let real_operations = msg_database::get_spacegroup_operations(uni, first_hall)
        .ok_or_else(|| format!("MAPI UNI {uni} first Hall {first_hall} lookup failed"))?;
    compare_magnetic_alias(uni, first_hall, &alias_operations, &real_operations)?;

    let alias_transformations = msg_database::get_std_transformations(uni, 0)
        .ok_or_else(|| format!("TAPI UNI {uni} Hall=0 alias lookup failed"))?;
    let real_transformations = msg_database::get_std_transformations(uni, first_hall)
        .ok_or_else(|| format!("TAPI UNI {uni} first Hall {first_hall} lookup failed"))?;
    compare_transformation_alias(
        uni,
        first_hall,
        &alias_transformations,
        &real_transformations,
    )?;
    Ok(())
}

fn build_mapi() -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    let mut visited = vec![false; MSG_OPERATION_COUNT];
    let mut span_count = 0usize;
    for uni in 1..MSG_UNI_COUNT {
        let mapping = msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni];
        let count = usize::try_from(mapping[0]).map_err(|_| format!("MAPI UNI {uni} count negative"))?;
        let first = usize::try_from(mapping[1]).map_err(|_| format!("MAPI UNI {uni} first negative"))?;
        validate_uni_zero_alias(uni, first)?;
        for slot in 0..count {
            let hall = first + slot;
            if !(1..=SPG_HALL_SETTINGS).contains(&hall) {
                return error(format!("MAPI UNI {uni} Hall {hall} out of range"));
            }
            let raw_span = msg_database::MAGNETIC_SPACEGROUP_OPERATION_INDEX[uni][slot];
            let order = usize::try_from(raw_span[0]).map_err(|_| format!("MAPI {uni}/{hall} order negative"))?;
            let offset = usize::try_from(raw_span[1]).map_err(|_| format!("MAPI {uni}/{hall} offset negative"))?;
            if order == 0 || offset == 0 || offset + order > MSG_OPERATION_COUNT {
                return error(format!("MAPI {uni}/{hall} span invalid"));
            }
            for index in offset..offset + order {
                if visited[index] {
                    return error(format!("MAPI raw index {index} visited twice"));
                }
                visited[index] = true;
            }
            let operations = msg_database::get_spacegroup_operations(uni, hall)
                .ok_or_else(|| format!("MAPI UNI {uni} Hall {hall} runtime lookup failed"))?;
            if operations.len() != order
                || operations.rot.len() != operations.trans.len()
                || operations.rot.len() != operations.timerev.len()
                || operations.rot.len() != order
            {
                return error(format!("MAPI {uni}/{hall} span/API order mismatch"));
            }
            append_u16(&mut payload, checked_u16(uni, &format!("MAPI {uni} UNI"))?);
            append_u16(&mut payload, checked_u16(hall, &format!("MAPI {uni}/{hall} Hall"))?);
            append_u16(&mut payload, checked_u16(order, &format!("MAPI {uni}/{hall} count"))?);
            for (index, ((&rotation, &translation), &time_reversal)) in operations
                .rot
                .iter()
                .zip(operations.trans.iter())
                .zip(operations.timerev.iter())
                .enumerate()
            {
                let raw_index = offset + index;
                let encoded = msg_database::MAGNETIC_SYMMETRY_OPERATIONS[raw_index];
                check_magnetic_runtime_operation(
                    encoded,
                    rotation,
                    translation,
                    time_reversal,
                    &format!("MAPI {uni}/{hall}[{index}]"),
                )?;
                append_magnetic_operation(
                    &mut payload,
                    rotation,
                    translation,
                    time_reversal,
                    &format!("MAPI {uni}/{hall}[{index}]"),
                )?;
            }
            span_count += 1;
        }
    }
    if span_count != MSG_ACTIVE_SPAN_COUNT {
        return error(format!("MAPI span census {span_count} != {MSG_ACTIVE_SPAN_COUNT}"));
    }
    let visited_count = visited.iter().skip(1).filter(|&&value| value).count();
    if visited[0] || visited_count != MSG_OPERATION_COUNT - 1
        || visited.iter().skip(1).any(|&value| !value)
    {
        return error(format!("MAPI visited raw indices {visited_count}, expected {}", MSG_OPERATION_COUNT - 1));
    }
    Ok(payload)
}

fn check_identity(rotation: [[i32; 3]; 3], translation: [f64; 3], label: &str) -> Result<()> {
    if rotation != [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
        return error(format!("{label} is not identity rotation"));
    }
    for (index, value) in translation.into_iter().enumerate() {
        if exact_translation(value, &format!("{label} translation[{index}]"))? != 0 {
            return error(format!("{label} is not zero translation"));
        }
    }
    Ok(())
}

fn build_tapi() -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    let mut occurrences = 0usize;
    let mut record_count = 0usize;
    for uni in 1..MSG_UNI_COUNT {
        let mapping = msg_database::MAGNETIC_SPACEGROUP_UNI_MAPPING[uni];
        let count = usize::try_from(mapping[0]).map_err(|_| format!("TAPI UNI {uni} count negative"))?;
        let first = usize::try_from(mapping[1]).map_err(|_| format!("TAPI UNI {uni} first negative"))?;
        for slot in 0..count {
            let hall = first + slot;
            let raw_values = msg_database::ALTERNATIVE_TRANSFORMATIONS[uni][slot];
            let first_zero = raw_values
                .iter()
                .position(|&value| value == 0)
                .ok_or_else(|| format!("TAPI {uni}/{hall} lacks alternative terminator"))?;
            if first_zero > 6 || raw_values[first_zero..].iter().any(|&value| value != 0) {
                return error(format!("TAPI {uni}/{hall} alternative tail invalid"));
            }
            let transformations = msg_database::get_std_transformations(uni, hall)
                .ok_or_else(|| format!("TAPI UNI {uni} Hall {hall} runtime lookup failed"))?;
            if transformations.len() != first_zero + 1 {
                return error(format!("TAPI {uni}/{hall} count mismatch"));
            }
            if transformations.rot.len() != transformations.trans.len()
                || transformations.rot.len() != first_zero + 1
            {
                return error(format!("TAPI {uni}/{hall} runtime vector lengths mismatch"));
            }
            let (identity_rotation, identity_translation) =
                (transformations.rot[0], transformations.trans[0]);
            check_identity(identity_rotation, identity_translation, &format!("TAPI {uni}/{hall}[0]"))?;
            append_u16(&mut payload, checked_u16(uni, &format!("TAPI {uni} UNI"))?);
            append_u16(&mut payload, checked_u16(hall, &format!("TAPI {uni}/{hall} Hall"))?);
            append_u16(&mut payload, checked_u16(transformations.len(), &format!("TAPI {uni}/{hall} count"))?);
            for index in 0..transformations.len() {
                let rotation = transformations.rot[index];
                let translation = transformations.trans[index];
                if index > 0 {
                    let encoded = raw_values[index - 1];
                    let (expected_rotation, expected_translation) = spg_database::decode_symmetry(encoded);
                    if rotation != expected_rotation {
                        return error(format!("TAPI {uni}/{hall}[{index}] rotation mismatch"));
                    }
                    for component in 0..3 {
                        let actual = exact_translation(translation[component], &format!("TAPI {uni}/{hall}[{index}] translation[{component}]"))?;
                        let expected = exact_translation(expected_translation[component], &format!("TAPI {uni}/{hall}[{index}] expected translation[{component}]"))?;
                        if actual != expected {
                            return error(format!("TAPI {uni}/{hall}[{index}] translation mismatch"));
                        }
                    }
                }
                append_spg_operation(&mut payload, rotation, translation, &format!("TAPI {uni}/{hall}[{index}]"))?;
            }
            occurrences += first_zero;
            record_count += 1;
        }
    }
    if record_count != MSG_ACTIVE_SPAN_COUNT {
        return error(format!("TAPI record census {record_count} != {MSG_ACTIVE_SPAN_COUNT}"));
    }
    if occurrences != ALT_VALUE_COUNT {
        return error(format!("TAPI raw occurrence census {occurrences} != {ALT_VALUE_COUNT}"));
    }
    Ok(payload)
}

fn build_frame() -> Result<Vec<u8>> {
    validate_zero_alias_contract()?;
    let sgno = build_sgno()?;
    let sgix = build_sgix()?;
    let sgrw = build_sgrw()?;
    let mtyp = build_mtyp()?;
    let muni = build_muni()?;
    let mhll = build_mhll()?;
    let midx = build_midx()?;
    let mraw = build_mraw()?;
    let malt = build_malt()?;
    let sdec = build_sdec()?;
    let sapi = build_sapi()?;
    let mapi = build_mapi()?;
    let tapi = build_tapi()?;

    let mut frame = Vec::new();
    frame.extend_from_slice(MAGIC);
    append_u32(&mut frame, VERSION);
    append_u32(&mut frame, SECTION_COUNT);
    append_section(&mut frame, b"SGNO", SPG_HALL_COUNT, sgno)?;
    append_section(&mut frame, b"SGIX", SPG_HALL_COUNT, sgix)?;
    append_section(&mut frame, b"SGRW", SPG_OPERATION_COUNT, sgrw)?;
    append_section(&mut frame, b"MTYP", MSG_UNI_COUNT, mtyp)?;
    append_section(&mut frame, b"MUNI", MSG_UNI_COUNT, muni)?;
    append_section(&mut frame, b"MHLL", SPG_HALL_COUNT, mhll)?;
    append_section(&mut frame, b"MIDX", MSG_UNI_COUNT * MSG_HALL_SLOTS, midx)?;
    append_section(&mut frame, b"MRAW", MSG_OPERATION_COUNT, mraw)?;
    append_section(&mut frame, b"MALT", MSG_UNI_COUNT * MSG_HALL_SLOTS, malt)?;
    append_section(&mut frame, b"SDEC", SPG_OPERATION_COUNT - 1, sdec)?;
    append_section(&mut frame, b"SAPI", SPG_HALL_SETTINGS, sapi)?;
    append_section(&mut frame, b"MAPI", MSG_ACTIVE_SPAN_COUNT, mapi)?;
    append_section(&mut frame, b"TAPI", MSG_ACTIVE_SPAN_COUNT, tapi)?;
    Ok(frame)
}

fn run() -> Result<()> {
    let frame = build_frame()?;
    let stdout = io::stdout();
    let mut output = io::BufWriter::new(stdout.lock());
    output.write_all(&frame).map_err(|error| format!("write frame: {error}"))?;
    output.flush().map_err(|error| format!("flush frame: {error}"))?;
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("magnetic database parity dump failed: {error}");
        std::process::exit(1);
    }
}
