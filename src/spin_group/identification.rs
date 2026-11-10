//! Identification of a detected collinear spin group against the pinned
//! collinear spin-space-group database.
//!
//! The identification is deliberately restricted to the **standard-cell
//! exact** case: the input operation data must already be expressed in the
//! database record's basis (`G0` or `L0`) and conventional cell.  No basis or
//! supercell search is performed here.  A frame change that genuinely changes
//! the operation data produces
//! [`SpinSpaceGroupIdentificationError::NoExactMatch`] rather than a silently
//! guessed number.  Origin choices are a gauge: if they leave the operation
//! data unchanged, they still identify the same record.
//!
//! The collinear order parameter is a signed scalar.  Its finite quotient is
//! fixed by the spatial operation and the channel sign; the full spin-only
//! stabilizer of the collinear order parameter (continuous rotations about the
//! magnetic axis, axis-preserving mirrors, and the anti-unitary spin-only
//! partner) is the canonical kernel and is not part of the signature.
//! Consequently no oriented magnetic axis is needed and a global spin rotation
//! of the magnetic structure does not change the identified number.
//!
//! Fractional translations are compared with [`DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE`]
//! (or an explicit tolerance through the `*_with_tolerance` entry points);
//! rotations and channel signs are compared exactly.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::collinear::CollinearAction;
use super::database::{
    SpinSpaceGroupNumber, SpinSpaceGroupOperation, SpinSpaceGroupRecord,
    collinear_spin_space_groups,
};
use super::{CollinearSpinGroup, Mat3I, Vec3};

/// Default fractional-translation tolerance for standard-cell exact
/// identification.
pub const DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE: f64 = 1e-8;

/// Largest accepted explicit fractional-translation tolerance.
///
/// The explicit tolerance is meant for sub-tolerance floating-point noise, not
/// for merging genuinely distinct fractional translations.  The pinned release
/// has a smallest same-rotation separation of `1/3`; `1e-2` stays far below
/// that and is exactly the bound validated by the independent final review.
pub const MAX_COLLINEAR_IDENTIFICATION_TOLERANCE: f64 = 1e-2;

const TRANSLATION_QUANTUM: f64 = 1e9;
const NORMALIZATION_EPSILON: f64 = 1e-12;

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
    /// A fractional translation is NaN or infinite.
    #[error("non-finite fractional translation {value} on axis {axis}")]
    NonFiniteTranslation { axis: usize, value: f64 },
    /// The requested matching tolerance is not a finite value in
    /// `(0, MAX_COLLINEAR_IDENTIFICATION_TOLERANCE]`.
    #[error("invalid identification tolerance {tolerance}")]
    InvalidTolerance { tolerance: f64 },
    /// No pinned record matches the signature in the input frame.
    ///
    /// This is the explicit refusal for genuinely non-standard cells,
    /// settings, supercells and structures not covered by the pinned database.
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

fn validate_tolerance(tolerance: f64) -> Result<(), SpinSpaceGroupIdentificationError> {
    if !tolerance.is_finite()
        || tolerance <= 0.0
        || tolerance > MAX_COLLINEAR_IDENTIFICATION_TOLERANCE
    {
        return Err(SpinSpaceGroupIdentificationError::InvalidTolerance { tolerance });
    }
    Ok(())
}

fn normalize_translation(value: f64) -> f64 {
    let mut value = value.rem_euclid(1.0);
    if (value - 1.0).abs() < NORMALIZATION_EPSILON || value.abs() < NORMALIZATION_EPSILON {
        value = 0.0;
    }
    value
}

fn quantize_translation(value: f64) -> i64 {
    (normalize_translation(value) * TRANSLATION_QUANTUM).round() as i64
}

fn flat_rotation(rotation: Mat3I) -> [i32; 9] {
    let mut flat = [0i32; 9];
    for (row, rotation_row) in rotation.iter().enumerate() {
        for (column, value) in rotation_row.iter().enumerate() {
            flat[3 * row + column] = *value;
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

fn translations_close(left: Vec3, right: Vec3, tolerance: f64) -> bool {
    (0..3).all(|axis| {
        let difference = (normalize_translation(left[axis]) - normalize_translation(right[axis]))
            .rem_euclid(1.0);
        difference < tolerance || (1.0 - difference) < tolerance
    })
}

fn channel_sign(channel: CollinearAction) -> i64 {
    match channel {
        CollinearAction::Preserve => 1,
        CollinearAction::Swap => -1,
    }
}

impl Signature {
    fn from_group(
        group: &CollinearSpinGroup,
        tolerance: f64,
    ) -> Result<Self, SpinSpaceGroupIdentificationError> {
        Self::new(
            group
                .operations()
                .iter()
                .map(|operation| (operation.rotation, operation.translation, operation.action)),
            tolerance,
        )
    }

    fn new(
        operations: impl IntoIterator<Item = (Mat3I, Vec3, CollinearAction)>,
        tolerance: f64,
    ) -> Result<Self, SpinSpaceGroupIdentificationError> {
        validate_tolerance(tolerance)?;
        let mut signature: Vec<SignatureOperation> = Vec::new();
        for (rotation, translation, channel) in operations {
            for (axis, value) in translation.into_iter().enumerate() {
                if !value.is_finite() {
                    return Err(SpinSpaceGroupIdentificationError::NonFiniteTranslation {
                        axis,
                        value,
                    });
                }
            }
            let normalized = [
                normalize_translation(translation[0]),
                normalize_translation(translation[1]),
                normalize_translation(translation[2]),
            ];
            if let Some(existing) = signature.iter_mut().find(|existing| {
                existing.rotation == rotation
                    && translations_close(existing.translation, normalized, tolerance)
            }) {
                if existing.channel != channel {
                    return Err(SpinSpaceGroupIdentificationError::InconsistentChannel {
                        rotation,
                    });
                }
                continue;
            }
            signature.push(SignatureOperation {
                rotation,
                translation: normalized,
                channel,
            });
        }
        if signature.is_empty() {
            return Err(SpinSpaceGroupIdentificationError::EmptyGroup);
        }
        signature.sort_by_key(|operation| operation_key(operation.rotation, operation.translation));
        Ok(Self {
            operations: signature,
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

    fn matches_record(&self, record: &SpinSpaceGroupRecord, tolerance: f64) -> bool {
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
                    && translations_close(candidate.translation, operation.translation, tolerance)
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

fn record_signature(record: &SpinSpaceGroupRecord) -> Signature {
    let mut operations: Vec<_> = record
        .operations
        .iter()
        .map(|operation| SignatureOperation {
            rotation: operation.rotation,
            translation: operation.translation.map(normalize_translation),
            channel: operation.channel,
        })
        .collect();
    operations.sort_by_key(|operation| operation_key(operation.rotation, operation.translation));
    Signature { operations }
}

fn exact_index() -> &'static ExactIndex {
    static INDEX: OnceLock<ExactIndex> = OnceLock::new();
    INDEX.get_or_init(|| {
        let mut index: ExactIndex = HashMap::new();
        for record in collinear_spin_space_groups() {
            let signature = record_signature(record);
            index.entry(signature.exact_key()).or_default().push(record);
        }
        index
    })
}

fn identify_signature(
    signature: &Signature,
    tolerance: f64,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    if let Some(candidates) = exact_index().get(&signature.exact_key()) {
        let matches: Vec<_> = candidates
            .iter()
            .copied()
            .filter(|record| signature.matches_record(record, tolerance))
            .collect();
        return resolve_matches(matches);
    }

    // Fallback for signatures whose fractional translations differ from the
    // pinned table only by sub-tolerance floating-point noise.
    let matches: Vec<_> = collinear_spin_space_groups()
        .filter(|record| signature.matches_record(record, tolerance))
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
/// record's standard cell and basis; otherwise
/// [`SpinSpaceGroupIdentificationError::NoExactMatch`] is returned.  The
/// anti-unitary kind is not part of the quotient signature: for the same
/// spatial operation both kinds carry the same channel sign.  Fractional
/// translations are compared with
/// [`DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE`]; use
/// [`identify_collinear_spin_space_group_with_tolerance`] to override it.
pub fn identify_collinear_spin_space_group(
    group: &CollinearSpinGroup,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    identify_collinear_spin_space_group_with_tolerance(
        group,
        DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE,
    )
}

/// Identify a detected collinear spin group with an explicit fractional
/// translation tolerance.
pub fn identify_collinear_spin_space_group_with_tolerance(
    group: &CollinearSpinGroup,
    tolerance: f64,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    validate_tolerance(tolerance)?;
    let signature = Signature::from_group(group, tolerance)?;
    identify_signature(&signature, tolerance)
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
    identify_collinear_signature_with_tolerance(
        operations,
        DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE,
    )
}

/// Identify a quotient signature with an explicit fractional-translation
/// tolerance.
pub fn identify_collinear_signature_with_tolerance(
    operations: &[SpinSpaceGroupOperation],
    tolerance: f64,
) -> Result<CollinearSsgIdentification, SpinSpaceGroupIdentificationError> {
    validate_tolerance(tolerance)?;
    let signature = Signature::new(
        operations
            .iter()
            .map(|operation| (operation.rotation, operation.translation, operation.channel)),
        tolerance,
    )?;
    identify_signature(&signature, tolerance)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOLERANCE: f64 = DEFAULT_COLLINEAR_IDENTIFICATION_TOLERANCE;

    #[test]
    fn every_pinned_signature_identifies_itself() {
        for record in collinear_spin_space_groups() {
            let signature = Signature::new(
                record.operations.iter().map(|operation| {
                    (operation.rotation, operation.translation, operation.channel)
                }),
                TOLERANCE,
            )
            .expect("pinned signature is valid");
            let identification = identify_signature(&signature, TOLERANCE)
                .expect("pinned signature must identify itself");
            assert_eq!(identification.number, record.number);
            assert_eq!(identification.record.identifier, record.identifier);
        }
    }

    #[test]
    fn empty_signature_is_rejected() {
        assert!(matches!(
            Signature::new(Vec::new(), TOLERANCE),
            Err(SpinSpaceGroupIdentificationError::EmptyGroup)
        ));
    }

    #[test]
    fn inconsistent_channel_is_rejected() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let result = Signature::new(
            [
                (identity, [0.0; 3], CollinearAction::Preserve),
                (identity, [0.0; 3], CollinearAction::Swap),
            ],
            TOLERANCE,
        );
        assert!(matches!(
            result,
            Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
        ));
    }

    #[test]
    fn near_duplicate_conflicting_channel_is_rejected() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let result = Signature::new(
            [
                (identity, [0.0; 3], CollinearAction::Preserve),
                (identity, [5e-9, 0.0, 0.0], CollinearAction::Swap),
            ],
            TOLERANCE,
        );
        assert!(matches!(
            result,
            Err(SpinSpaceGroupIdentificationError::InconsistentChannel { .. })
        ));
    }

    #[test]
    fn same_channel_near_duplicate_collapses() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let signature = Signature::new(
            [
                (identity, [0.0; 3], CollinearAction::Preserve),
                (identity, [5e-9, 0.0, 0.0], CollinearAction::Preserve),
            ],
            TOLERANCE,
        )
        .expect("near duplicate should collapse");
        assert_eq!(signature.operations.len(), 1);
        let identification = identify_signature(&signature, TOLERANCE).unwrap();
        assert_eq!(identification.identifier(), "1.1.1.1.L");
    }

    #[test]
    fn non_finite_translation_is_rejected() {
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let result = Signature::new(
                [(identity, [0.0, value, 0.0], CollinearAction::Preserve)],
                TOLERANCE,
            );
            assert!(matches!(
                result,
                Err(SpinSpaceGroupIdentificationError::NonFiniteTranslation { .. })
            ));
        }
    }

    #[test]
    fn invalid_tolerance_is_rejected() {
        for tolerance in [0.0, -1e-8, f64::NAN, f64::INFINITY, 1.1e-2, 0.5] {
            let result = Signature::new(Vec::new(), tolerance);
            assert!(matches!(
                result,
                Err(SpinSpaceGroupIdentificationError::InvalidTolerance { .. })
            ));
        }
    }

    #[test]
    fn explicit_tolerance_is_capped_below_operation_separation() {
        let record = collinear_spin_space_group_by_identifier_for_test("146.146.1.1.L");
        let signature = Signature::new(
            record
                .operations
                .iter()
                .map(|operation| (operation.rotation, operation.translation, operation.channel)),
            MAX_COLLINEAR_IDENTIFICATION_TOLERANCE,
        )
        .expect("pinned signature is valid");
        let identification = identify_signature(&signature, MAX_COLLINEAR_IDENTIFICATION_TOLERANCE)
            .expect("the capped tolerance must stay safe");
        assert_eq!(identification.number, record.number);
        assert!(matches!(
            Signature::new(Vec::new(), 1.1e-2),
            Err(SpinSpaceGroupIdentificationError::InvalidTolerance { .. })
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
        let signature = Signature::new(operations, TOLERANCE).expect("valid signature");
        assert!(matches!(
            identify_signature(&signature, TOLERANCE),
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
        let signature = Signature::new(operations, TOLERANCE).expect("valid signature");
        assert!(matches!(
            identify_signature(&signature, TOLERANCE),
            Err(SpinSpaceGroupIdentificationError::NoExactMatch)
        ));
    }

    #[test]
    fn tolerance_can_absorb_small_noise() {
        let record = collinear_spin_space_group_by_identifier_for_test("221.221.1.1.L");
        let noise = 5e-8;
        let operations: Vec<_> = record
            .operations
            .iter()
            .map(|operation| {
                (
                    operation.rotation,
                    [
                        operation.translation[0] + noise,
                        operation.translation[1] + noise,
                        operation.translation[2] + noise,
                    ],
                    operation.channel,
                )
            })
            .collect();
        let default = Signature::new(operations.clone(), TOLERANCE).expect("valid signature");
        assert!(matches!(
            identify_signature(&default, TOLERANCE),
            Err(SpinSpaceGroupIdentificationError::NoExactMatch)
        ));
        let relaxed = Signature::new(operations, 1e-7).expect("valid signature");
        let identification = identify_signature(&relaxed, 1e-7).expect("noise is within tolerance");
        assert_eq!(identification.identifier(), "221.221.1.1.L");
    }

    fn collinear_spin_space_group_by_identifier_for_test(
        identifier: &str,
    ) -> &'static SpinSpaceGroupRecord {
        super::super::database::collinear_spin_space_group_by_identifier(identifier)
            .expect("test record exists")
    }
}
