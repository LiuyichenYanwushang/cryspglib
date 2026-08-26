//! Validated crystallographic database identifiers.
//!
//! Space-group, Hall, and UNI numbers are large finite domains. They are
//! represented as bounded newtypes rather than enums: an enum with hundreds or
//! thousands of variants would duplicate the generated database while making
//! lookup and regeneration less convenient.

use std::fmt;
use std::num::{NonZeroU8, NonZeroU16};

/// Error returned when a raw integer is outside an identifier's database
/// range.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("{kind} number {value} is outside {minimum}..={maximum}")]
pub struct InvalidIdentifier {
    kind: &'static str,
    value: usize,
    minimum: usize,
    maximum: usize,
}

impl InvalidIdentifier {
    /// Name of the rejected identifier domain.
    pub const fn kind(self) -> &'static str {
        self.kind
    }

    /// Rejected raw value.
    pub const fn value(self) -> usize {
        self.value
    }

    /// Inclusive valid range.
    pub const fn valid_range(self) -> (usize, usize) {
        (self.minimum, self.maximum)
    }
}

/// International space-group number in `1..=230`.
///
/// The non-zero representation also keeps `Option<SpaceGroupNumber>` compact.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpaceGroupNumber(NonZeroU8);

impl SpaceGroupNumber {
    pub const MIN: usize = 1;
    pub const MAX: usize = 230;

    /// Construct from a raw number, returning `None` outside `1..=230`.
    pub const fn new(value: u8) -> Option<Self> {
        if value == 0 || value as usize > Self::MAX {
            None
        } else {
            match NonZeroU8::new(value) {
                Some(value) => Some(Self(value)),
                None => None,
            }
        }
    }

    /// Return the database index as `usize`.
    pub const fn get(self) -> usize {
        self.0.get() as usize
    }
}

impl TryFrom<usize> for SpaceGroupNumber {
    type Error = InvalidIdentifier;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        u8::try_from(value)
            .ok()
            .and_then(Self::new)
            .ok_or(InvalidIdentifier {
                kind: "space-group",
                value,
                minimum: Self::MIN,
                maximum: Self::MAX,
            })
    }
}

impl From<SpaceGroupNumber> for usize {
    fn from(value: SpaceGroupNumber) -> Self {
        value.get()
    }
}

impl fmt::Display for SpaceGroupNumber {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.get().fmt(formatter)
    }
}

macro_rules! bounded_nonzero_u16 {
    ($name:ident, $kind:literal, $maximum:literal, $docs:literal) => {
        #[doc = $docs]
        #[repr(transparent)]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(NonZeroU16);

        impl $name {
            pub const MIN: usize = 1;
            pub const MAX: usize = $maximum;

            /// Construct from a raw number, returning `None` outside the
            /// database range.
            pub const fn new(value: u16) -> Option<Self> {
                if value == 0 || value as usize > Self::MAX {
                    None
                } else {
                    match NonZeroU16::new(value) {
                        Some(value) => Some(Self(value)),
                        None => None,
                    }
                }
            }

            /// Return the database index as `usize`.
            pub const fn get(self) -> usize {
                self.0.get() as usize
            }
        }

        impl TryFrom<usize> for $name {
            type Error = InvalidIdentifier;

            fn try_from(value: usize) -> Result<Self, Self::Error> {
                u16::try_from(value)
                    .ok()
                    .and_then(Self::new)
                    .ok_or(InvalidIdentifier {
                        kind: $kind,
                        value,
                        minimum: Self::MIN,
                        maximum: Self::MAX,
                    })
            }
        }

        impl From<$name> for usize {
            fn from(value: $name) -> Self {
                value.get()
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.get().fmt(formatter)
            }
        }
    };
}

bounded_nonzero_u16!(HallNumber, "Hall", 530, "Hall setting number in `1..=530`.");
bounded_nonzero_u16!(
    UniNumber,
    "UNI",
    1651,
    "UNI magnetic-space-group number in `1..=1651`."
);

#[cfg(test)]
mod tests {
    use super::{HallNumber, SpaceGroupNumber, UniNumber};
    use std::mem::size_of;

    #[test]
    fn boundaries_are_validated_once() {
        assert_eq!(SpaceGroupNumber::try_from(1).unwrap().get(), 1);
        assert_eq!(SpaceGroupNumber::try_from(230).unwrap().get(), 230);
        assert!(SpaceGroupNumber::try_from(0).is_err());
        assert!(SpaceGroupNumber::try_from(231).is_err());

        assert_eq!(HallNumber::try_from(530).unwrap().get(), 530);
        assert!(HallNumber::try_from(0).is_err());
        assert!(HallNumber::try_from(531).is_err());

        assert_eq!(UniNumber::try_from(1651).unwrap().get(), 1651);
        assert!(UniNumber::try_from(0).is_err());
        assert!(UniNumber::try_from(1652).is_err());
        assert!(UniNumber::try_from(usize::MAX).is_err());
    }

    #[test]
    fn representations_and_options_stay_compact() {
        assert_eq!(size_of::<SpaceGroupNumber>(), 1);
        assert_eq!(size_of::<Option<SpaceGroupNumber>>(), 1);
        assert_eq!(size_of::<HallNumber>(), 2);
        assert_eq!(size_of::<Option<HallNumber>>(), 2);
        assert_eq!(size_of::<UniNumber>(), 2);
        assert_eq!(size_of::<Option<UniNumber>>(), 2);
    }
}
