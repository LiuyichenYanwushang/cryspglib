#!/usr/bin/env python3
"""Freeze the audited direct-source data--Hall authority as JSON.

This module is a serializer and validator for the already-derived authority;
it is not a second Hall-selection implementation.  ``build_artifact`` calls
only the public, zero-argument data--Hall derivation and the public exact
source loader used to obtain the two audit anchors.  The frame convention is
the direct source frame ``P=I, p=0``.  A ``source_to_hall`` shift means
``source = hall + shift`` (the unwrapped integer ``L``), while a
``hall_to_source`` shift means ``hall = source + shift`` (the unwrapped
centering/normalization ``C``); every shift is stored as an integer numerator
over 12.

The committed sidecar binds the external PIR/CIR ZIPs and the committed
spglib artifact/manifest by their fixed path, byte length, and SHA-256.  It
does not bind Python source hashes: the code is versioned by Git, and keeping
those hashes here would make the freezer self-referential.  This stage has no
runtime loader.  The final artifact/manifest byte commitments below are used
only as an integrity gate for this frozen pair, never as Hall-search input.
"""

from __future__ import annotations

from collections import Counter
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Optional


SCHEMA = "cryspglib-iso-irrep-data-hall-v1"
MANIFEST_SCHEMA = "cryspglib-iso-irrep-data-hall-manifest-v1"
GENERATOR_SCHEMA_VERSION = "1"
TRANSLATION_DENOMINATOR = 12
FRAME_SEMANTICS = "direct-source-frame-P=I-p=0"
MAPPING_SEMANTICS = {
    "source_to_hall": "source=hall+shift",
    "hall_to_source": "hall=source+shift",
}

_MODULE_DIR = Path(__file__).resolve().parent
_REPOSITORY_ROOT = _MODULE_DIR.parent

_CENTERING_RESIDUES = {
    "P": ((0, 0, 0),),
    "A": ((0, 0, 0), (0, 6, 6)),
    "B": ((0, 0, 0), (6, 0, 6)),
    "C": ((0, 0, 0), (6, 6, 0)),
    "F": ((0, 0, 0), (0, 6, 6), (6, 0, 6), (6, 6, 0)),
    "I": ((0, 0, 0), (6, 6, 6)),
    "R": ((0, 0, 0), (4, 8, 8), (8, 4, 4)),
}
_CENTERING_ORDER = ("P", "A", "B", "C", "F", "I", "R")
_IDENTITY_BASIS = [1, 0, 0, 0, 1, 0, 0, 0, 1]
_ZERO_ORIGIN = [0, 0, 0]

# These are the only production commitments consumed by this serializer.
# They intentionally duplicate the already-pinned values so the pre/post byte
# check below detects any input replacement before a sidecar can be emitted.
_INPUT_SPECS = (
    (
        "pir_zip", "isotropy_subgroup/PIR_data.zip", 1_235_319,
        "e909a4f0121688b0590ccaec10b0276171bc24619cf7eb562ba441268c01e121",
    ),
    (
        "cir_zip", "isotropy_subgroup/CIR_data.zip", 1_555_153,
        "f4edcb2852b83a86d1b58f29fb862d9124a227cfc90f9e1ae17d2c97585264e6",
    ),
    (
        "spglib_artifact", "scripts/data/spglib_magnetic_provenance_v1.json",
        1_537_875,
        "933a52a6696e7f6a1a2e426825ad92c377c6e96330e18c5c045d659798d740b9",
    ),
    (
        "spglib_manifest",
        "scripts/data/spglib_magnetic_provenance_v1.manifest.json", 570,
        "6a9e1b64c190c30a556d63e51e5b896b967d33e8821714beb745ae699fab84bf",
    ),
)

_EXPECTED_CENSUS = {
    "pir_records": 10_294,
    "cir_records": 11_202,
    "source_representatives": 2_609,
    "raw_unique": 220,
    "raw_ambiguous": 10,
    "raw_missing": 0,
    "raw_ambiguous_spacegroups": [5, 8, 9, 12, 15, 21, 38, 39, 65, 67],
    "filtered_unique": 230,
    "filtered_ambiguous": 0,
    "filtered_missing": 0,
    "selected_hall_operations": 4_425,
    "source_to_hall": 2_609,
    "source_to_hall_nonzero": 0,
    "hall_to_source": 4_425,
    "hall_to_source_nonzero": 1_816,
    "expanded_normalization_nonzero": 410,
    "centering_counts": [
        ["P", 149], ["A", 4], ["B", 0], ["C", 16],
        ["F", 16], ["I", 38], ["R", 7],
    ],
}

# These commitments identify the one sidecar pair produced by this schema
# version.  They are an integrity boundary, not a derivation or Hall search
# input.  A changed derivation must deliberately introduce a new schema/frozen
# pair rather than silently becoming accepted by this validator.
FINAL_ARTIFACT_BYTE_LENGTH = 697_730
FINAL_ARTIFACT_SHA256 = (
    "35bcb00958021eb6fc5a330f8dbf85a80be78ccec324f441e6138cdba4b617e0"
)
FINAL_MANIFEST_BYTE_LENGTH = 869
FINAL_MANIFEST_SHA256 = (
    "bc6aa7a94d698f2193e7cb623b16dded2dd8e0307d502cf28b68c554f364d7e2"
)


class FreezeError(ValueError):
    """Base class for typed sidecar construction and validation failures."""


class FreezeIntegrityError(FreezeError):
    """A committed input, byte stream, or atomic output is not trustworthy."""


class FreezeSchemaError(FreezeError):
    """A JSON value does not obey the exact sidecar schema."""


class FreezeInvariantError(FreezeError):
    """A structurally valid sidecar violates a semantic invariant."""


def _error(context: str, message: str):
    raise FreezeSchemaError(f"{context}: {message}")


def _validate_json_tree(value, context: str = "$"):
    """Require the exact built-in JSON value types used by this sidecar."""

    if type(value) is dict:
        for key, child in value.items():
            if type(key) is not str:
                _error(context, "object keys must be exact strings")
            if any(ord(char) < 0x20 or ord(char) > 0x7E for char in key):
                _error(context, "object keys must be printable ASCII")
            _validate_json_tree(child, f"{context}.{key}")
        return
    if type(value) is list:
        for index, child in enumerate(value):
            _validate_json_tree(child, f"{context}[{index}]")
        return
    if type(value) is str:
        if any(ord(char) < 0x20 or ord(char) > 0x7E for char in value):
            _error(context, "strings must be printable ASCII")
        return
    if type(value) is int:
        return
    _error(context, "null, bool, float, tuple, and non-built-in values are forbidden")


def canonical_json(value) -> bytes:
    """Encode one exact printable-ASCII JSON tree with one final LF."""

    try:
        _validate_json_tree(value)
        encoded = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except RecursionError as error:
        raise FreezeSchemaError("JSON value is too deeply nested") from error
    except (TypeError, ValueError, UnicodeError) as error:
        raise FreezeSchemaError("unable to encode canonical JSON") from error
    return encoded + b"\n"


def _pairs_without_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise FreezeSchemaError(f"duplicate JSON object key {key!r}")
        result[key] = value
    return result


def _parse_canonical_json(data: bytes, context: str):
    if type(data) is not bytes:
        raise FreezeSchemaError(f"{context} must be exact bytes")
    if not data.endswith(b"\n"):
        raise FreezeSchemaError(f"{context} must have exactly one final LF")
    if any(byte >= 0x80 for byte in data):
        raise FreezeSchemaError(f"{context} must be ASCII UTF-8")
    try:
        text = data[:-1].decode("utf-8", errors="strict")
        value = json.loads(
            text,
            object_pairs_hook=_pairs_without_duplicates,
            parse_float=lambda value: (_ for _ in ()).throw(
                FreezeSchemaError(f"{context} contains a float")
            ),
            parse_constant=lambda value: (_ for _ in ()).throw(
                FreezeSchemaError(f"{context} contains {value}")
            ),
        )
    except FreezeError:
        raise
    except RecursionError as error:
        raise FreezeSchemaError(f"{context} is too deeply nested") from error
    except (UnicodeError, ValueError, TypeError) as error:
        raise FreezeSchemaError(f"{context} is not valid JSON") from error
    try:
        _validate_json_tree(value, context)
        if canonical_json(value) != data:
            raise FreezeSchemaError(f"{context} is not canonical JSON")
    except FreezeError:
        raise
    except RecursionError as error:
        raise FreezeSchemaError(f"{context} is too deeply nested") from error
    return value


def _require_keys(value, keys, context: str):
    if type(value) is not dict:
        raise FreezeSchemaError(f"{context} must be an object")
    if set(value) != set(keys):
        raise FreezeSchemaError(f"{context} has missing or extra keys")


def _require_int(value, context: str, *, minimum: Optional[int] = None):
    if type(value) is not int:
        raise FreezeSchemaError(f"{context} must be an exact integer")
    if minimum is not None and value < minimum:
        raise FreezeInvariantError(f"{context} is below {minimum}")
    return value


def _require_str(value, context: str):
    if type(value) is not str:
        raise FreezeSchemaError(f"{context} must be an exact string")
    return value


def _input_descriptors():
    return {
        name: {"path": path, "bytes": byte_count, "sha256": digest}
        for name, path, byte_count, digest in _INPUT_SPECS
    }


def _read_external_inputs():
    """Read and verify each external input exactly once per pass."""

    observed = []
    for name, relative_path, expected_bytes, expected_sha in _INPUT_SPECS:
        path = _REPOSITORY_ROOT / relative_path
        try:
            payload = path.read_bytes()
        except OSError as error:
            raise FreezeIntegrityError(f"unable to read {name}: {path}") from error
        if len(payload) != expected_bytes:
            raise FreezeIntegrityError(f"{name} byte length changed")
        digest = hashlib.sha256(payload).hexdigest()
        if digest != expected_sha:
            raise FreezeIntegrityError(f"{name} SHA-256 changed")
        observed.append((name, payload))
    return tuple(observed)


def _load_authority_modules():
    """Import local authority modules only after the pre-read gate."""

    try:
        from . import derive_iso_irrep_data_hall as derive
        from . import iso_irrep_exact as exact
    except ImportError:  # pragma: no cover - direct scripts/ invocation
        import derive_iso_irrep_data_hall as derive
        import iso_irrep_exact as exact
    return derive, exact


def _fraction_numerator_over_12(value, context: str) -> int:
    scaled = value * TRANSLATION_DENOMINATOR
    if scaled.denominator != 1:
        raise FreezeInvariantError(f"{context} is not exact over 12")
    return scaled.numerator


def _frame_to_json(frame, universe):
    try:
        pir_anchor = universe.pir_irnumbers[0]
        cir_anchor = universe.cir_irnumbers[0]
    except (AttributeError, IndexError, TypeError) as error:
        raise FreezeInvariantError(
            f"SG{frame.spacegroup} lacks a PIR/CIR anchor"
        ) from error
    if type(pir_anchor) is not int or pir_anchor <= 0:
        raise FreezeInvariantError(f"SG{frame.spacegroup} PIR anchor is invalid")
    if type(cir_anchor) is not int or cir_anchor <= 0:
        raise FreezeInvariantError(f"SG{frame.spacegroup} CIR anchor is invalid")
    if universe.space_group_symbol != frame.source_symbol:
        raise FreezeInvariantError(f"SG{frame.spacegroup} source symbol changed")
    if universe.centering.value != frame.centering:
        raise FreezeInvariantError(f"SG{frame.spacegroup} centering changed")
    return {
        "spacegroup": frame.spacegroup,
        "source_symbol": frame.source_symbol,
        "centering": frame.centering,
        "pir_anchor_irnumber": pir_anchor,
        "cir_anchor_irnumber": cir_anchor,
        "raw_candidate_halls": list(frame.raw_candidate_halls),
        "data_hall": frame.data_hall,
        "basis": list(_IDENTITY_BASIS),
        "origin_numerator": list(_ZERO_ORIGIN),
        "source_operation_count": frame.source_operation_count,
        "hall_operation_count": frame.hall_operation_count,
        "source_to_hall": [
            {
                "source_operation_index": mapping.source_operation_index,
                "hall_operation_index": mapping.hall_operation_index,
                "lattice_shift_numerator": list(mapping.shift_numerator),
            }
            for mapping in frame.source_to_hall
        ],
        "hall_to_source": [
            {
                "hall_operation_index": mapping.hall_operation_index,
                "source_operation_index": mapping.source_operation_index,
                "lattice_shift_numerator": list(mapping.shift_numerator),
            }
            for mapping in frame.hall_to_source
        ],
    }


def _census_to_json(census):
    def distribution(value):
        return [[list(key), count] for key, count in value]

    return {
        "pir_records": census.pir_records,
        "cir_records": census.cir_records,
        "source_representatives": census.source_representatives,
        "raw_unique": census.raw_unique,
        "raw_ambiguous": census.raw_ambiguous,
        "raw_missing": census.raw_missing,
        "raw_ambiguous_spacegroups": list(census.raw_ambiguous_spacegroups),
        "filtered_unique": census.filtered_unique,
        "filtered_ambiguous": census.filtered_ambiguous,
        "filtered_missing": census.filtered_missing,
        "selected_hall_operations": census.selected_hall_operations,
        "source_to_hall": census.source_to_hall,
        "source_to_hall_nonzero": census.source_to_hall_nonzero,
        "hall_to_source": census.hall_to_source,
        "hall_to_source_nonzero": census.hall_to_source_nonzero,
        "hall_to_source_shifts": distribution(census.hall_to_source_shifts),
        "hall_to_source_cosets": distribution(census.hall_to_source_cosets),
        "expanded_normalization_nonzero": census.expanded_normalization_nonzero,
        "expanded_normalization_shifts": distribution(
            census.expanded_normalization_shifts
        ),
        "centering_counts": [list(row) for row in census.centering_counts],
    }


def build_artifact() -> dict:
    """Build one validated artifact from the pinned authority chain."""

    pre = _read_external_inputs()
    try:
        derive, exact = _load_authority_modules()
        authority = derive.derive_data_hall_frames()
        if type(authority) is not derive.ExactDataHallDatabase:
            raise FreezeInvariantError("public derivation returned a wrong type")
        frames, census = tuple(authority)
        source_database = exact.load_exact_iso_irrep_sources()
        spacegroups = []
        for spacegroup in range(1, 231):
            frame = frames[spacegroup - 1]
            universe = source_database.source_universe(spacegroup)
            spacegroups.append(
                _frame_to_json(frame, universe)
            )
        artifact = {
            "schema": SCHEMA,
            "translation_denominator": TRANSLATION_DENOMINATOR,
            "frame_semantics": FRAME_SEMANTICS,
            "mapping_semantics": dict(MAPPING_SEMANTICS),
            "inputs": _input_descriptors(),
            "census": _census_to_json(census),
            "spacegroups": spacegroups,
        }
        _validate_artifact(artifact)
        # Serialization is part of the pre/post mutation window.  The caller
        # may serialize again; canonical_json is deterministic.
        canonical_json(artifact)
    except FreezeError:
        raise
    except Exception as error:
        raise FreezeInvariantError("authority derivation failed") from error
    post = _read_external_inputs()
    if post != pre:
        raise FreezeIntegrityError("external input bytes changed during derivation")
    return artifact


def _validate_distribution(value, context: str):
    if type(value) is not list:
        raise FreezeSchemaError(f"{context} must be an array")
    rows = []
    for index, row in enumerate(value):
        if type(row) is not list or len(row) != 2:
            raise FreezeSchemaError(f"{context}[{index}] is malformed")
        key, count = row
        if type(key) is not list or len(key) != 3:
            raise FreezeSchemaError(f"{context}[{index}] key is malformed")
        key_tuple = tuple(_require_int(item, f"{context}[{index}] key") for item in key)
        count = _require_int(count, f"{context}[{index}] count", minimum=1)
        rows.append((key_tuple, count))
    if rows != sorted(rows) or len({key for key, _ in rows}) != len(rows):
        raise FreezeInvariantError(f"{context} is not canonical")
    return rows


def _validate_mapping(mapping, expected_keys, context: str):
    _require_keys(mapping, expected_keys, context)
    first = _require_int(mapping[expected_keys[0]], f"{context}.{expected_keys[0]}", minimum=0)
    second = _require_int(mapping[expected_keys[1]], f"{context}.{expected_keys[1]}", minimum=0)
    shift = mapping["lattice_shift_numerator"]
    if type(shift) is not list or len(shift) != 3:
        raise FreezeSchemaError(f"{context}.lattice_shift_numerator is malformed")
    shift = tuple(_require_int(value, f"{context}.lattice_shift_numerator") for value in shift)
    return first, second, shift


def _validate_record(record):
    keys = (
        "spacegroup", "source_symbol", "centering", "pir_anchor_irnumber",
        "cir_anchor_irnumber", "raw_candidate_halls", "data_hall", "basis",
        "origin_numerator", "source_operation_count", "hall_operation_count",
        "source_to_hall", "hall_to_source",
    )
    _require_keys(record, keys, "spacegroup record")
    spacegroup = _require_int(record["spacegroup"], "record.spacegroup", minimum=1)
    if spacegroup > 230:
        raise FreezeInvariantError("record.spacegroup is outside 1..230")
    symbol = _require_str(record["source_symbol"], "record.source_symbol")
    centering = _require_str(record["centering"], "record.centering")
    if centering not in _CENTERING_RESIDUES:
        raise FreezeInvariantError("record.centering is unknown")
    if not symbol or symbol[0] != centering:
        raise FreezeInvariantError("record symbol/centering mismatch")
    for name in ("pir_anchor_irnumber", "cir_anchor_irnumber"):
        _require_int(record[name], f"record.{name}", minimum=1)
    halls = record["raw_candidate_halls"]
    if type(halls) is not list or not halls:
        raise FreezeSchemaError("record.raw_candidate_halls is malformed")
    halls_tuple = tuple(_require_int(value, "record.raw_candidate_halls", minimum=1) for value in halls)
    if any(value > 530 for value in halls_tuple) or halls_tuple != tuple(sorted(set(halls_tuple))):
        raise FreezeInvariantError("record.raw_candidate_halls is not ordered")
    data_hall = _require_int(record["data_hall"], "record.data_hall", minimum=1)
    if data_hall > 530 or data_hall not in halls_tuple:
        raise FreezeInvariantError("record.data_hall is not a raw candidate")
    if record["basis"] != _IDENTITY_BASIS:
        raise FreezeInvariantError("record basis is not identity")
    if record["origin_numerator"] != _ZERO_ORIGIN:
        raise FreezeInvariantError("record origin is not zero")
    source_count = _require_int(record["source_operation_count"], "record.source_operation_count", minimum=1)
    hall_count = _require_int(record["hall_operation_count"], "record.hall_operation_count", minimum=1)
    residues = _CENTERING_RESIDUES[centering]
    if hall_count != source_count * len(residues):
        raise FreezeInvariantError("record operation counts disagree")

    source_maps = record["source_to_hall"]
    hall_maps = record["hall_to_source"]
    if type(source_maps) is not list or len(source_maps) != source_count:
        raise FreezeSchemaError("record.source_to_hall cardinality mismatch")
    if type(hall_maps) is not list or len(hall_maps) != hall_count:
        raise FreezeSchemaError("record.hall_to_source cardinality mismatch")
    source_to_hall = []
    for index, mapping in enumerate(source_maps):
        source_index, hall_index, shift = _validate_mapping(
            mapping,
            ("source_operation_index", "hall_operation_index", "lattice_shift_numerator"),
            f"record.source_to_hall[{index}]",
        )
        if source_index != index:
            raise FreezeInvariantError("source_to_hall is not source ordered")
        if any(value % TRANSLATION_DENOMINATOR for value in shift):
            raise FreezeInvariantError("source_to_hall shift is not integral")
        source_to_hall.append((source_index, hall_index, shift))
    if len({hall_index for _, hall_index, _ in source_to_hall}) != source_count:
        raise FreezeInvariantError("source_to_hall Hall indices are not unique")
    if any(hall_index >= hall_count for _, hall_index, _ in source_to_hall):
        raise FreezeInvariantError("source_to_hall Hall index is out of range")

    hall_to_source = []
    residues_by_source = [[] for _ in range(source_count)]
    for index, mapping in enumerate(hall_maps):
        hall_index, source_index, shift = _validate_mapping(
            mapping,
            ("hall_operation_index", "source_operation_index", "lattice_shift_numerator"),
            f"record.hall_to_source[{index}]",
        )
        if hall_index != index or source_index >= source_count:
            raise FreezeInvariantError("hall_to_source index is out of range/order")
        residue = tuple(value % TRANSLATION_DENOMINATOR for value in shift)
        residues_by_source[source_index].append(residue)
        hall_to_source.append((hall_index, source_index, shift))
    expected_residues = set(residues)
    for source_index, source_residues in enumerate(residues_by_source):
        if len(source_residues) != len(expected_residues) or set(source_residues) != expected_residues:
            raise FreezeInvariantError(
                f"hall_to_source residues incomplete for source {source_index}"
            )
    for source_index, hall_index, shift in source_to_hall:
        inverse = hall_to_source[hall_index]
        if inverse[1] != source_index or inverse[2] != tuple(-value for value in shift):
            raise FreezeInvariantError("source/Hall mapping directions disagree")
    return source_count, hall_count, halls_tuple, source_to_hall, hall_to_source


def _aggregate_census(spacegroups, pir_records, cir_records):
    raw_unique = raw_ambiguous = raw_missing = 0
    raw_ambiguous_sgs = []
    centering_counts = Counter()
    source_representatives = selected_hall_operations = 0
    source_to_hall = source_to_hall_nonzero = 0
    hall_to_source = hall_to_source_nonzero = 0
    hall_shifts = Counter()
    hall_cosets = Counter()
    expanded_shifts = Counter()
    for record in spacegroups:
        source_count, hall_count, halls, source_maps, hall_maps = _validate_record(record)
        if len(halls) == 0:
            raw_missing += 1
        elif len(halls) == 1:
            raw_unique += 1
        else:
            raw_ambiguous += 1
            raw_ambiguous_sgs.append(record["spacegroup"])
        centering_counts[record["centering"]] += 1
        source_representatives += source_count
        selected_hall_operations += hall_count
        source_to_hall += len(source_maps)
        source_to_hall_nonzero += sum(shift != (0, 0, 0) for _, _, shift in source_maps)
        hall_to_source += len(hall_maps)
        for _, _, shift in hall_maps:
            hall_shifts[shift] += 1
            if shift != (0, 0, 0):
                hall_to_source_nonzero += 1
            residue = tuple(value % TRANSLATION_DENOMINATOR for value in shift)
            hall_cosets[residue] += 1
            expanded_shifts[tuple(
                residue[index] - shift[index] for index in range(3)
            )] += 1

    def distribution(counter):
        return [[list(key), count] for key, count in sorted(counter.items())]

    return {
        "pir_records": pir_records,
        "cir_records": cir_records,
        "source_representatives": source_representatives,
        "raw_unique": raw_unique,
        "raw_ambiguous": raw_ambiguous,
        "raw_missing": raw_missing,
        "raw_ambiguous_spacegroups": raw_ambiguous_sgs,
        "filtered_unique": 230,
        "filtered_ambiguous": 0,
        "filtered_missing": 0,
        "selected_hall_operations": selected_hall_operations,
        "source_to_hall": source_to_hall,
        "source_to_hall_nonzero": source_to_hall_nonzero,
        "hall_to_source": hall_to_source,
        "hall_to_source_nonzero": hall_to_source_nonzero,
        "hall_to_source_shifts": distribution(hall_shifts),
        "hall_to_source_cosets": distribution(hall_cosets),
        "expanded_normalization_nonzero": sum(
            count for shift, count in expanded_shifts.items()
            if shift != (0, 0, 0)
        ),
        "expanded_normalization_shifts": distribution(expanded_shifts),
        "centering_counts": [
            [name, centering_counts.get(name, 0)] for name in _CENTERING_ORDER
        ],
    }


def _validate_census(census, spacegroups):
    keys = (
        "pir_records", "cir_records", "source_representatives", "raw_unique",
        "raw_ambiguous", "raw_missing", "raw_ambiguous_spacegroups",
        "filtered_unique", "filtered_ambiguous", "filtered_missing",
        "selected_hall_operations", "source_to_hall", "source_to_hall_nonzero",
        "hall_to_source", "hall_to_source_nonzero", "hall_to_source_shifts",
        "hall_to_source_cosets", "expanded_normalization_nonzero",
        "expanded_normalization_shifts", "centering_counts",
    )
    _require_keys(census, keys, "census")
    for key in keys:
        if key not in ("raw_ambiguous_spacegroups", "hall_to_source_shifts", "hall_to_source_cosets", "expanded_normalization_shifts", "centering_counts"):
            _require_int(census[key], f"census.{key}", minimum=0)
    ambiguous = census["raw_ambiguous_spacegroups"]
    if type(ambiguous) is not list or any(type(value) is not int for value in ambiguous):
        raise FreezeSchemaError("census.raw_ambiguous_spacegroups is malformed")
    if ambiguous != sorted(set(ambiguous)):
        raise FreezeInvariantError("census raw ambiguity list is not ordered")
    for key in ("hall_to_source_shifts", "hall_to_source_cosets", "expanded_normalization_shifts"):
        _validate_distribution(census[key], f"census.{key}")
    rows = census["centering_counts"]
    if type(rows) is not list or len(rows) != len(_CENTERING_ORDER):
        raise FreezeSchemaError("census.centering_counts is malformed")
    for index, row in enumerate(rows):
        if type(row) is not list or len(row) != 2 or type(row[0]) is not str:
            raise FreezeSchemaError("census.centering_counts row is malformed")
        if row[0] != _CENTERING_ORDER[index]:
            raise FreezeInvariantError("census centering order changed")
        _require_int(row[1], "census centering count", minimum=0)
    expected = _aggregate_census(
        spacegroups, census["pir_records"], census["cir_records"]
    )
    if census != expected:
        raise FreezeInvariantError("census disagrees with spacegroup mappings")
    for key, value in _EXPECTED_CENSUS.items():
        if census[key] != value:
            raise FreezeInvariantError(f"known census mismatch in {key}")


def _validate_inputs(inputs):
    _require_keys(inputs, tuple(name for name, *_ in _INPUT_SPECS), "inputs")
    expected = _input_descriptors()
    for name, descriptor in inputs.items():
        _require_keys(descriptor, ("path", "bytes", "sha256"), f"inputs.{name}")
        if descriptor != expected[name]:
            raise FreezeIntegrityError(f"inputs.{name} does not match pinned commitment")


def _validate_artifact(artifact):
    keys = (
        "schema", "translation_denominator", "frame_semantics",
        "mapping_semantics", "inputs", "census", "spacegroups",
    )
    _require_keys(artifact, keys, "artifact")
    if artifact["schema"] != SCHEMA:
        raise FreezeSchemaError("artifact schema mismatch")
    if artifact["translation_denominator"] != TRANSLATION_DENOMINATOR:
        raise FreezeInvariantError("translation denominator mismatch")
    if artifact["frame_semantics"] != FRAME_SEMANTICS:
        raise FreezeInvariantError("frame semantics mismatch")
    _require_keys(artifact["mapping_semantics"], ("source_to_hall", "hall_to_source"), "mapping_semantics")
    if artifact["mapping_semantics"] != MAPPING_SEMANTICS:
        raise FreezeInvariantError("mapping semantics mismatch")
    _validate_inputs(artifact["inputs"])
    spacegroups = artifact["spacegroups"]
    if type(spacegroups) is not list or len(spacegroups) != 230:
        raise FreezeSchemaError("spacegroups must be a 230-entry array")
    for index, record in enumerate(spacegroups, 1):
        if type(record) is not dict or record.get("spacegroup") != index:
            raise FreezeInvariantError("spacegroups are not ordered 1..230")
        _validate_record(record)
    _validate_census(artifact["census"], spacegroups)


def build_manifest(artifact_bytes: bytes) -> dict:
    """Build the independent manifest for one already-canonical artifact."""

    artifact = _parse_canonical_json(artifact_bytes, "artifact")
    _validate_artifact(artifact)
    return {
        "schema": MANIFEST_SCHEMA,
        "generator_schema_version": GENERATOR_SCHEMA_VERSION,
        "inputs": _input_descriptors(),
        "artifact": {
            "path": "iso_irrep_data_hall_v1.json",
            "bytes": len(artifact_bytes),
            "sha256": hashlib.sha256(artifact_bytes).hexdigest(),
        },
    }


def _parse_and_validate_uncommitted_pair(artifact_bytes: bytes, manifest_bytes: bytes):
    """Validate a canonical pair structurally for internal build/test use.

    This seam deliberately does not accept the frozen authority commitment;
    it is never an authority acceptance boundary.
    """

    artifact = _parse_canonical_json(artifact_bytes, "artifact")
    manifest = _parse_canonical_json(manifest_bytes, "manifest")
    _validate_artifact(artifact)
    _require_keys(
        manifest,
        ("schema", "generator_schema_version", "inputs", "artifact"),
        "manifest",
    )
    if manifest["schema"] != MANIFEST_SCHEMA:
        raise FreezeSchemaError("manifest schema mismatch")
    if manifest["generator_schema_version"] != GENERATOR_SCHEMA_VERSION:
        raise FreezeSchemaError("manifest generator version mismatch")
    _validate_inputs(manifest["inputs"])
    if manifest["inputs"] != artifact["inputs"]:
        raise FreezeInvariantError("manifest and artifact inputs differ")
    _require_keys(manifest["artifact"], ("path", "bytes", "sha256"), "manifest.artifact")
    artifact_entry = manifest["artifact"]
    if artifact_entry["path"] != "iso_irrep_data_hall_v1.json":
        raise FreezeSchemaError("manifest artifact path mismatch")
    if artifact_entry["bytes"] != len(artifact_bytes):
        raise FreezeIntegrityError("manifest artifact byte length mismatch")
    if artifact_entry["sha256"] != hashlib.sha256(artifact_bytes).hexdigest():
        raise FreezeIntegrityError("manifest artifact SHA-256 mismatch")
    return artifact, manifest


def _verify_final_bytes(data, expected_length: int, expected_sha256: str,
                        context: str):
    if type(data) is not bytes:
        raise FreezeSchemaError(f"{context} must be exact bytes")
    if len(data) != expected_length:
        raise FreezeIntegrityError(f"{context} byte length is not frozen")
    if hashlib.sha256(data).hexdigest() != expected_sha256:
        raise FreezeIntegrityError(f"{context} SHA-256 is not frozen")


def parse_and_validate_pair(artifact_bytes: bytes, manifest_bytes: bytes):
    """Accept only the committed frozen pair, then run full validation.

    The final byte/SHA gate is an integrity check for this sidecar pair.  It
    is not a Hall-selection input and does not participate in derivation.
    General structural checking for build/test seams is private.
    """

    _verify_final_bytes(
        artifact_bytes, FINAL_ARTIFACT_BYTE_LENGTH, FINAL_ARTIFACT_SHA256,
        "artifact",
    )
    _verify_final_bytes(
        manifest_bytes, FINAL_MANIFEST_BYTE_LENGTH, FINAL_MANIFEST_SHA256,
        "manifest",
    )
    return _parse_and_validate_uncommitted_pair(artifact_bytes, manifest_bytes)


def _path_argument(value, context: str) -> Path:
    try:
        raw = os.fspath(value)
    except TypeError as error:
        raise FreezeSchemaError(f"{context} is not a path") from error
    if type(raw) is not str:
        raise FreezeSchemaError(f"{context} is not a native path")
    if "\x00" in raw:
        raise FreezeSchemaError(f"{context} contains a NUL byte")
    try:
        os.fsencode(raw)
        return Path(raw).absolute()
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeSchemaError(f"{context} is invalid") from error


def _reject_same_target(output: Path, manifest: Path):
    if output == manifest:
        raise FreezeInvariantError("artifact and manifest paths must differ")
    try:
        if output.exists() and manifest.exists() and os.path.samefile(output, manifest):
            raise FreezeInvariantError("artifact and manifest paths share an inode")
    except FreezeError:
        raise
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeIntegrityError("unable to inspect output inodes") from error


def _replacement_identity(path: Path, context: str) -> Path:
    """Resolve the directory used by ``os.replace`` before any build work."""

    try:
        parent = path.parent.resolve(strict=True)
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeIntegrityError(
            f"{context} parent cannot be resolved"
        ) from error
    try:
        is_directory = parent.is_dir()
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeIntegrityError(
            f"{context} parent cannot be inspected"
        ) from error
    if not is_directory:
        raise FreezeIntegrityError(f"{context} parent is not a directory")
    return parent / path.name


def _reject_same_replacement_identity(output: Path, manifest: Path):
    output_identity = _replacement_identity(output, "output")
    manifest_identity = _replacement_identity(manifest, "manifest")
    if output_identity == manifest_identity:
        raise FreezeInvariantError(
            "artifact and manifest replacement identities must differ"
        )
    _reject_same_target(output_identity, manifest_identity)
    return output_identity, manifest_identity


def _make_temp(path: Path, payload: bytes) -> Path:
    try:
        parent_exists = path.parent.exists()
        parent_is_directory = path.parent.is_dir()
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeIntegrityError(
            f"unable to inspect output directory: {path.parent}"
        ) from error
    if not parent_exists or not parent_is_directory:
        raise FreezeIntegrityError(f"output directory does not exist: {path.parent}")
    descriptor = None
    temporary = None
    try:
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
        )
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = None
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        return Path(temporary)
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except (OSError, ValueError):
                pass
        if temporary is not None:
            try:
                os.unlink(temporary)
            except (OSError, ValueError, UnicodeError):
                pass
        raise FreezeIntegrityError(f"unable to stage {path}") from error


def _fsync_directory(path: Path):
    descriptor = None
    try:
        descriptor = os.open(str(path), os.O_RDONLY)
        os.fsync(descriptor)
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeIntegrityError(f"unable to fsync output directory {path}") from error
    finally:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except (OSError, ValueError):
                pass


def write_outputs(output, manifest):
    """Build and atomically replace two files in commit-marker order.

    Each file is replaced atomically, with the manifest written last as the
    commit marker.  The pair is not transactional: if the second replacement
    fails, a new artifact may coexist with a missing or old manifest, and a
    consumer will fail closed on the pair's frozen commitment mismatch.
    """

    output_path = _path_argument(output, "output")
    manifest_path = _path_argument(manifest, "manifest")
    output_path, manifest_path = _reject_same_replacement_identity(
        output_path, manifest_path
    )
    artifact = build_artifact()
    artifact_bytes = canonical_json(artifact)
    manifest_bytes = canonical_json(build_manifest(artifact_bytes))
    parse_and_validate_pair(artifact_bytes, manifest_bytes)
    temporary_paths = []
    try:
        temporary_paths.append(_make_temp(output_path, artifact_bytes))
        temporary_paths.append(_make_temp(manifest_path, manifest_bytes))
        os.replace(str(temporary_paths[0]), str(output_path))
        temporary_paths[0] = None
        os.replace(str(temporary_paths[1]), str(manifest_path))
        temporary_paths[1] = None
        _fsync_directory(output_path.parent)
        if manifest_path.parent != output_path.parent:
            _fsync_directory(manifest_path.parent)
    except FreezeError:
        raise
    except (OSError, RuntimeError, ValueError, UnicodeError) as error:
        raise FreezeIntegrityError("atomic output replacement failed") from error
    finally:
        for temporary in temporary_paths:
            if temporary is not None:
                try:
                    os.unlink(str(temporary))
                except (OSError, ValueError, UnicodeError):
                    pass
    return artifact_bytes, manifest_bytes


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args(argv)
    try:
        write_outputs(args.output, args.manifest)
    except FreezeError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = [
    "FreezeError", "FreezeIntegrityError", "FreezeInvariantError",
    "FreezeSchemaError", "MANIFEST_SCHEMA", "SCHEMA",
    "canonical_json", "build_artifact", "build_manifest",
    "parse_and_validate_pair", "write_outputs",
]
