//! Identification of a detected collinear spin group against the pinned
//! collinear spin-space-group database.
//!
//! The identification is deliberately restricted to the **standard-cell
//! exact** case: the input crystal must already be described in the database
//! record's basis (`G0` or `L0`) and conventional cell.  No basis, origin or
//! supercell search is performed here.  A non-standard description therefore
//! produces [`SpinSpaceGroupIdentificationError::NoExactMatch`] rather than a
//! silently guessed number.
//!
//! The collinear order parameter is a signed scalar.  Its finite quotient is
//! fixed by the spatial operation and the channel sign; the continuous spin
//! rotations about the magnetic axis and the universal spin-only stabilizer are
//! the canonical kernel and are not part of the signature.  Consequently no
//! oriented magnetic axis is needed and a global spin rotation of the magnetic
//! structure does not change the identified number.

use std::collections::{BTreeMap, HashMap};
use std::sync::OnceLock;

use super::collinear::CollinearAction;
use super::database::{
    SpinSpaceGroupNumber, SpinSpaceGroupOperation, SpinSpaceGroupRecord,
    collinear_spin_space_groups,
};
use super::{CollinearSpinGroup, Mat3I, Vec3};

const TRANSLATION_TOLERANCE: f64 = 1e-8;
const TRANSLATION_QUANTUM: f64 = 1e9;

/// Result of a successful collinear spin-space-group identification.
#[derive(Debug, Clone, Copy)]
pub struct CollinearSsgIdentification {
    /// Canonical Chen-Liu identifier.
    pub number: SpinSpaceGroupNumber,
    /// Pinned database record.
    pub record: &'static SpinSpaceGroupRecord,
}

impl CollinearSsgIdentification {
    /// Canonical identifier string, e.g. `221.221.1.1.L`.
    pub fn identifier(self) -> String {
        self.number.to_string()
    }

    /// Plain-text International notation of the identified record.
    pub fn international_symbol(self) -> &'static str {
        self.record.international_symbol
    }
}

/// Failure modes of standard-cell exact collinear spin-space-group
/// identification.
#[derive(Debug, thiserror::Error)]
pub enum SpinSpaceGroupIdentificationError {
    /// The detected group has no operations.
    #[error("cannot identify an empty spin group")]
    EmptyGroup,
    /// Two operations with the same spatial part carry different channel
    /// signs.  This cannot happen for a validated group and indicates an
    /// inconsistent or corrupted certificate.
    #[error("inconsistent channel for spatial rotation {rotation:?}")]
    InconsistentChannel { rotation: Mat3I },
    /// No pinned record matches the signature in the input cell/frame.
    ///
    /// This is the explicit refusal for non-standard cells, settings,
    /// supercells and structures not covered by the pinned database.
    #[error(
        "no pinned collinear spin-space group matches this signature in the database standard cell; \
         the input cell/setting may be non-standard or the group is not in the pinned release"
    )]
    NoExactMatch,
    /// More than one pinned record matches.  The current release has unique
    /// quotient signatures, so this is a defensive error.
    #[error("multiple pinned records match the collinear signature: {identifiers:?}")]
    Ambiguous { identifiers: Vec<&'static str> },
}

/// One operation of the identification signature.
#[derive(Debug, Clone, Copy, PartialEq)]
struct SignatureOperation {
    rotation: Mat3I,
    translation: Vec3,
    channel: CollinearAction,
}

/// Canonical quotient signature used for database matching.
#[derive(Debug, Clone)]
struct Signature {
    operations: Vec<SignatureOperation>,
}

fn normalize_translation(value: f64) -> f64 {
    let mut value = value.rem_euclid(1.0);
    if (value - 1.0).abs() < TRANSLATION_TOLERANCE || value.abs() < TRANSLATION_TOLERANCE {
        value = 0.0;
    }
    value
}

fn quantize_translation(value: f64) -> i64 {
    (normalize_translation(value) * TRANSLATION_QUANTUM).round() as i64
}

fn flat_rotation(rotation: Mat3I) -> [i32; 9] {
    let mut flat = [0i32; 9];
    for row in 0..3 {
        for column in 0..3 {
            flat[3 * row + column] = rotation[row][column];
        }
    }
    flat
}

fn operation_key(rotation: Mat3I, translation: Vec3) -> ([i32; 9], [i64; 3]) {
    (
        flat_rotation(rotation),
        [
            quantize_translation(translation[0]),
            quantize_translation(translation[1]),
            quantize_translation(translation[2]),
        ],
    )
}

fn translations_close(left: Vec3, right: Vec3) -> bool {
    (0..3).all(|axis| {
        let difference = (normalize_translation(left[axis]) - normalize_translation(right[axis]))
            .rem_euclid(1.0);
        difference < TRANSLATION_TOLERANCE || (1.0 - difference) < TRANSLATION_TOLERANCE
    })
}

fn channel_sign(channel: CollinearAction) -> i64 {
    match channel {
        CollinearAction::Preserve => 1,
        CollinearAction::Swap => -1,
    }
}

impl Signature {
    fn from_group(group: &CollinearSpinGroup) -> Result<Self, SpinSpaceGroupIdentificationError> {
        Self::new(
            group
                .operations()
                .iter()
                .map(|operation| (operation.rotation, operation.translation, operation.action)),
        )
    }

    fn new(
        operations: impl IntoIterator<Item = (Mat3I, Vec3, CollinearAction)>,
    ) -> Result<Self, SpinSpaceGroupIdentificationError> {
        let mut table: BTreeMap<([i32; 9], [i64; 3]), SignatureOperation> = BTreeMap::new();
        for (rotation, translation, channel) in operations {
            let normalized = [
                normalize_translation(translation[0]),
                normalize_translation(translation[1]),
                normalize_translation(translation[2]),
            ];
            let key = operation_key(rotation, normalized);
            match table.entry(key) {
                std::collections::btree_map::Entry::Occupied(entry) => {
                    if entry.get().channel != channel {
                        return Err(SpinSpaceGroupIdentificationError::InconsistentChannel {
                            rotation,
                        });
                    }
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(SignatureOperation {
                        rotation,
                        translation: normalized,
                        channel,
                    });
                }
            }
        }
        if table.is_empty() {
            return Err(SpinSpaceGroupIdentificationError::EmptyGroup);
        }
        Ok(Self {
            operations: table.into_values().collect(),
        })
    }

    fn exact_key(&self) -> Vec<i64> {
        let mut key = Vec::with_capacity(13 * self.operations.len());
        for operation in &self.operations {
            let (rotation, translation) = operation_key(operation.rotation, operation.translation);
            key.extend(rotation.map(i64::from));
            key.extend(translation);
            key.push(channel_sign(operation.channel));
        }
        key
    }

    fn matches_record(&self, record: &SpinSpaceGroupRecord) -> bool {
        if record.operations.len() != self.operations.len() {
            return false;
        }
        let mut used = vec![false; record.operations.len()];
        for operation in &self.operations {
            let mut found = false;
            for (index, candidate) in record.operations.iter().enumerate() {
                if used[index] {
                    continue;
                }
                if candidate.rotation == operation.rotation
                    && translations_close(candidate.translation, operation.translation)
                {
                    if candidate.channel != operation.channel {
                        return false;
                    }
                    used[index] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
        }
        true
    }
}

type ExactIndex = HashMap<Vec<i64>, Vec<&'static SpinSpaceGroupRecord>>;

fn exact_index() -> &'static ExactIndex {
    static INDEX: OnceLock<ExactIndex> = OnceLock::new();
    INDEX.get_or_init(|| {
        let mut index: ExactIndex = HashMap::new();
        for record in collinear_spin_space_groups() {
            let mut operations = Vec::with_capacity(record.operations.len());
            for operation in record.operations {
                let normalized = [
                    normalize_translation(operation.translation[0]),
                    normalize_translation(operation.translation[1]),
                    normalize_translation(operation.translation[2]),
                ];
                operations.push(SignatureOperation {
                    rotation: operation.rotation,
                    translation: normalized,
                    channel: operation.channel,
                });
            }
            operations
                .sort_by_key(|operation| operation_key(operation.rotation, operation.translation));
            let signature = Signature { operations };
            index.entry(signature.exact_key()).or_default().push(record);
        }
        index
    })
}

fn identify_signature(
    signature: &Signature,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    if let Some(candidates) = exact_index().get(&signature.exact_key()) {
        let matches: Vec<_> = candidates
            .iter()
            .copied()
            .filter(|record| signature.matches_record(record))
            .collect();
        return resolve_matches(matches);
    }

    // Fallback for signatures whose fractional translations differ from the
    // pinned table only by sub-tolerance floating-point noise.
    let matches: Vec<_> = collinear_spin_space_groups()
        .filter(|record| signature.matches_record(record))
        .collect();
    resolve_matches(matches)
}

fn resolve_matches(
    matches: Vec<&'static SpinSpaceGroupRecord>,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    match matches.as_slice() {
        [] => Err(SpinSpaceGroupIdentificationError::NoExactMatch),
        [record] => Ok(CollinearSsgIdentification {
            number: record.number,
            record,
        }),
        _ => Err(SpinSpaceGroupIdentificationError::Ambiguous {
            identifiers: matches.iter().map(|record| record.identifier).collect(),
        }),
    }
}

/// Identify a detected collinear spin group against the pinned database.
///
/// The input group must describe a crystal already expressed in the database
/// record's standard cell and basis; otherwise [`SpinSpaceGroupIdentificationError::NoExactMatch`]
/// is returned.  The anti-unitary kind is not part of the quotient signature:
/// for the same spatial operation both kinds carry the same channel sign.
pub fn identify_collinear_spin_space_group(
    group: &CollinearSpinGroup,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    let signature = Signature::from_group(group)?;
    identify_signature(&signature)
}

/// Identify a quotient signature directly.
///
/// This is the lower-level entry point for consumers that already hold a
/// spatial-operation/channel signature (for example a pinned database record)
/// rather than a detected [`CollinearSpinGroup`].  The standard-cell exact
/// restriction and the error semantics are the same as
/// [`identify_collinear_spin_space_group`].
pub fn identify_collinear_signature(
    operations: &[SpinSpaceGroupOperation],
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    let signature = Signature::new(
        operations
            .iter()
            .map(|operation| (operation.rotation, operation.translation, operation.channel)),
    )?;
    identify_signature(&signature)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_pinned_signature_identifies_itself() {
        for record in collinear_spin_space_groups() {
            let signature =
                Signature::new(record.operations.iter().map(|operation| {
                    (operation.rotation, operation.translation, operation.channel)
                }))
                .expect("pinned signature is valid");
            let identification =
                identify_signature(&signature).expect("pinned signature must identify itself");
            assert_eq!(identification.number, record.number);
            assert_eq!(identification.record.identifier, record.identifier);
        }
    }

    #[test]
    fn empty_signature_is_rejected() {
        assert!(matches!(
            Signature::new(Vec::new()),
            Err(SpinSpaceGroupIdentificationError::EmptyGroup)
        ));
    }

    #[test]
    fn inconsistent_channel_is_rejected() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let result = Signature::new([
            (identity, [0.0; 3], CollinearAction::Preserve),
            (identity, [0.0; 3], CollinearAction::Swap),
        ]);
        assert!(matches!(
            result,
            Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
        ));
    }

    #[test]
    fn wrong_channel_does_not_match() {
        let record = collinear_spin_space_group_by_identifier_for_test("221.200.1.1.L");
        let mut operations: Vec<_> = record
            .operations
            .iter()
            .map(|operation| (operation.rotation, operation.translation, operation.channel))
            .collect();
        operations[0].2 = match operations[0].2 {
            CollinearAction::Preserve => CollinearAction::Swap,
            CollinearAction::Swap => CollinearAction::Preserve,
        };
        let signature = Signature::new(operations).expect("valid signature");
        assert!(matches!(
            identify_signature(&signature),
            Err(SpinSpaceGroupIdentificationError::NoExactMatch)
        ));
    }

    #[test]
    fn shifted_cell_does_not_match() {
        let record = collinear_spin_space_group_by_identifier_for_test("221.221.1.1.L");
        let operations = record.operations.iter().map(|operation| {
            (
                operation.rotation,
                [
                    operation.translation[0] + 0.1,
                    operation.translation[1] + 0.1,
                    operation.translation[2] + 0.1,
                ],
                operation.channel,
            )
        });
        let signature = Signature::new(operations).expect("valid signature");
        assert!(matches!(
            identify_signature(&signature),
            Err(SpinSpaceGroupIdentificationError::NoExactMatch)
        ));
    }

    fn collinear_spin_space_group_by_identifier_for_test(
        identifier: &str,
    ) -> &'static SpinSpaceGroupRecord {
        super::super::database::collinear_spin_space_group_by_identifier(identifier)
            .expect("test record exists")
    }
}
