"""Focused tests for the pinned generator archive boundary."""

import os
import sys
import struct
import tempfile
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
import zipfile
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))
import generate_irrep_data as generator
import parse_spinor_data


_EXPECTED_HALL_OPERATIONS_BYTE_LENGTH = 481408
_EXPECTED_HALL_OPERATIONS_SHA256 = (
    "ebd1cf36668fb8c0efd633b2d7728c51ca1b404a3cc02ed871ece47b46a0d1c8"
)


class PinnedArchiveBoundaryTests(unittest.TestCase):
    def setUp(self):
        self._iso_dir = generator.ISO_DIR
        self._hashes = generator.PINNED_ARCHIVE_SHA256

    def tearDown(self):
        generator.ISO_DIR = self._iso_dir
        generator.PINNED_ARCHIVE_SHA256 = self._hashes

    @staticmethod
    def _write_zip(directory, members):
        path = os.path.join(directory, "iso.zip")
        with zipfile.ZipFile(path, "w") as archive:
            for name, contents in members.items():
                archive.writestr(name, contents)
        return path

    @staticmethod
    def _hall_tree():
        path = os.path.join(generator.SCRIPT_DIR, "hall_operations.json")
        with open(path, "rb") as stream:
            return json.loads(stream.read().decode("utf-8"))

    @staticmethod
    def _hall_tree_payload(tree):
        return json.dumps(tree, separators=(",", ":")).encode("utf-8")

    def test_valid_zip_is_read_without_extracted_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            self._write_zip(directory, {"nested/data.txt": "from zip\n"})
            generator.ISO_DIR = directory
            with generator._open_zip_path("iso.zip", "data.txt") as stream:
                self.assertEqual(stream.readlines(), ["from zip\n"])

    def test_altered_extracted_file_cannot_override_zip(self):
        with tempfile.TemporaryDirectory() as directory:
            self._write_zip(directory, {"nested/data.txt": "from zip\n"})
            extracted = os.path.join(directory, "iso")
            os.makedirs(extracted)
            with open(os.path.join(extracted, "data.txt"), "w") as stream:
                stream.write("altered extracted data\n")
            generator.ISO_DIR = directory
            with generator._open_zip_path("iso.zip", "data.txt") as stream:
                self.assertEqual(stream.read(), "from zip\n")

    def test_missing_zip_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            generator.ISO_DIR = directory
            with self.assertRaisesRegex(FileNotFoundError, "required ZIP archive"):
                generator._open_zip_path("iso.zip", "data.txt")

    def test_bad_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            for archive in generator.PINNED_ARCHIVE_SHA256:
                with open(os.path.join(directory, archive), "wb") as stream:
                    stream.write(b"archive bytes")
            generator.ISO_DIR = directory
            expected = dict(generator.PINNED_ARCHIVE_SHA256)
            expected["PIR_data.zip"] = "0" * 64
            generator.PINNED_ARCHIVE_SHA256 = expected
            with self.assertRaisesRegex(ValueError, "pinned archive hash mismatch"):
                generator._verify_pinned_archives()

    def test_irreptables_record_provenance_is_pinned(self):
        expected = parse_spinor_data.IRREPTABLES_RECORD_SHA256
        try:
            parse_spinor_data.IRREPTABLES_RECORD_SHA256 = "0" * 64
            with self.assertRaisesRegex(ValueError, "RECORD hash mismatch"):
                parse_spinor_data._verified_irreptables_distribution()
        finally:
            parse_spinor_data.IRREPTABLES_RECORD_SHA256 = expected

    def test_irreptables_version_is_pinned(self):
        expected = parse_spinor_data.IRREPTABLES_VERSION
        try:
            parse_spinor_data.IRREPTABLES_VERSION = "0.0.0"
            with self.assertRaisesRegex(FileNotFoundError, "exactly one"):
                parse_spinor_data._verified_irreptables_distribution()
        finally:
            parse_spinor_data.IRREPTABLES_VERSION = expected

    def test_spin_source_manifest_rejects_missing_and_corrupt_files(self):
        with tempfile.TemporaryDirectory() as directory:
            package_root = os.path.dirname(directory)
            manifest = {}
            with self.assertRaisesRegex(ValueError, "expected 1 pinned spin source files"):
                parse_spinor_data._verify_spin_source_files(
                    directory, package_root, manifest, expected_count=1
                )
            source = os.path.join(directory, "irreps-SG=3-spin.dat")
            with open(source, "wb") as stream:
                stream.write(b"pinned source")
            relative = os.path.relpath(source, package_root).replace(os.sep, "/")
            manifest[relative] = "0" * 64
            with self.assertRaisesRegex(ValueError, "spin source file hash mismatch"):
                parse_spinor_data._verify_spin_source_files(
                    directory, package_root, manifest, expected_count=1
                )

    def test_spin_source_rows_keep_raw_file_order_ordinals(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix="-spin.dat") as stream:
            stream.write(
                "SG=3\n"
                "symmetries=\n"
                "kpoint GM : 0 0 0 : 1\n"
                "-GM1 1 1.0\n"
                "-GM2 1 1.0\n"
            )
            stream.flush()
            _sg, _ops, irreps = parse_spinor_data.parse_spinor_file(stream.name)
            self.assertEqual([row["source_row_ordinal"] for row in irreps], [0, 1])

    def test_spin_source_sg_duplicates_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "each SG exactly once"):
            parse_spinor_data._validate_spin_source_sgs(
                list(range(1, 230)) + [229]
            )

    def test_ambiguous_suffix_member_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            self._write_zip(directory, {
                "first/data.txt": "first",
                "second/data.txt": "second",
            })
            generator.ISO_DIR = directory
            with self.assertRaisesRegex(ValueError, "ambiguous archive member"):
                generator._open_zip_path("iso.zip", "data.txt")

    def test_legacy_hall_table_has_independent_pinned_commitment(self):
        path = os.path.join(generator.SCRIPT_DIR, "hall_operations.json")
        with open(path, "rb") as stream:
            payload = stream.read()
        self.assertEqual(len(payload), _EXPECTED_HALL_OPERATIONS_BYTE_LENGTH)
        self.assertEqual(
            hashlib.sha256(payload).hexdigest(), _EXPECTED_HALL_OPERATIONS_SHA256
        )
        self.assertEqual(
            generator.HALL_OPERATIONS_BYTE_LENGTH,
            _EXPECTED_HALL_OPERATIONS_BYTE_LENGTH,
        )
        self.assertEqual(
            generator.HALL_OPERATIONS_SHA256, _EXPECTED_HALL_OPERATIONS_SHA256
        )

    def test_hall_payload_parser_accepts_full_aggregate(self):
        path = os.path.join(generator.SCRIPT_DIR, "hall_operations.json")
        with open(path, "rb") as stream:
            payload = stream.read()
        parsed = generator._parse_hall_operations_payload(payload)
        self.assertEqual(sorted(parsed), list(range(1, 231)))
        self.assertEqual(
            sum(len(entry[1]) for entries in parsed.values() for entry in entries),
            7388,
        )

    def test_hall_payload_parser_rejects_missing_hall4(self):
        tree = self._hall_tree()
        del tree["4"]
        with self.assertRaisesRegex(ValueError, "root keys"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_hall_payload_parser_rejects_hall4_extra_field(self):
        tree = self._hall_tree()
        tree["4"]["extra"] = 0
        with self.assertRaisesRegex(ValueError, "Hall4 entry keys"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_hall_payload_parser_rejects_duplicate_root_key1(self):
        tree = self._hall_tree()
        entries = []
        encoded_first = json.dumps(tree["1"], separators=(",", ":"))
        entries.extend(["\"1\":" + encoded_first,
                        "\"1\":" + encoded_first])
        for hall_num in range(2, 531):
            hall_key = str(hall_num)
            entries.append(
                json.dumps(hall_key) + ":" +
                json.dumps(tree[hall_key], separators=(",", ":"))
            )
        payload = ("{" + ",".join(entries) + "}").encode("utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate JSON object key"):
            generator._parse_hall_operations_payload(payload)

    def test_hall_payload_parser_rejects_hall4_missing_trans(self):
        tree = self._hall_tree()
        del tree["4"]["trans"]
        with self.assertRaisesRegex(ValueError, "Hall4 entry keys"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_hall_payload_parser_rejects_nonfinite_translation(self):
        tree = self._hall_tree()
        tree["4"]["trans"][0][0] = float("nan")
        with self.assertRaisesRegex(ValueError, "non-finite"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_hall_payload_parser_rejects_parallel_length_mismatch(self):
        tree = self._hall_tree()
        tree["4"]["trans"].pop()
        with self.assertRaisesRegex(ValueError, "Hall4 rots/trans length"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_hall_payload_parser_rejects_duplicate_seitz_row(self):
        tree = self._hall_tree()
        hall_key = next(
            key for key in tree if len(tree[key]["rots"]) > 1)
        tree[hall_key]["rots"][1] = list(tree[hall_key]["rots"][0])
        tree[hall_key]["trans"][1] = list(tree[hall_key]["trans"][0])
        with self.assertRaisesRegex(ValueError, "duplicates a Seitz row"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_hall_payload_parser_rejects_nonexact_component_types(self):
        tree = self._hall_tree()
        tree["4"]["sg"] = True
        with self.assertRaisesRegex(ValueError, "Hall4 sg"):
            generator._parse_hall_operations_payload(self._hall_tree_payload(tree))

    def test_pir_parallel_offsets_are_strictly_linked(self):
        args = dict(
            sg=[4],
            ml=["GM1"],
            char_starts=[0],
            char_counts=[1],
            chars_flat=[1.0],
            pir_rots_flat=[1, 0, 0, 0, 1, 0, 0, 0, 1],
            pir_rot_starts=[0],
            pir_trans_flat=[0.0, 0.0, 0.0],
            pir_trans_starts=[0],
            little_chars_real=[1.0],
            little_chars_imag=[0.0],
            little_chars_valid=[1],
        )
        generator._validate_pir_storage_alignment(**args)
        for field, value in (
                ("pir_rot_starts", [1]),
                ("pir_trans_starts", [3]),
                ("char_starts", [1])):
            malformed = dict(args)
            malformed[field] = value
            with self.assertRaisesRegex(ValueError, "offset"):
                generator._validate_pir_storage_alignment(**malformed)


class DataHallSelectionTests(unittest.TestCase):
    @staticmethod
    def _sg5_source_and_target():
        database = generator.load_committed_data_hall_provenance()
        frame = database.frames[4]
        sg_halls = generator._load_hall_operations()
        selected = [entry for entry in sg_halls[5] if entry[0] == 9]
        if len(selected) != 1:
            raise AssertionError("synthetic SG5 witness needs unique Hall9")
        _hall_number, hall_rots, hall_trans = selected[0]
        return database, frame, hall_rots, hall_trans, sg_halls

    def test_sidecar_sg5_mapping_ignores_other_raw_candidates(self):
        _database, frame, hall_rots, hall_trans, sg_halls = (
            self._sg5_source_and_target())
        source_rots = hall_rots[:2]
        source_trans = hall_trans[:2]
        self.assertEqual(
            generator._sidecar_source_hall_mapping(
                frame, 5, "GM1", source_rots, source_trans,
                hall_rots, hall_trans),
            [0, 1, 0, 1],
        )

        permuted = {spacegroup: list(entries)
                    for spacegroup, entries in sg_halls.items()}
        permuted[5].reverse()
        choices = generator._prepare_sidecar_hall_choices(
            generator.load_committed_data_hall_provenance(), permuted)
        self.assertEqual(choices[5][0], 9)
        self.assertEqual(len(choices[5][2]), frame.hall_operation_count)
        self.assertIn(10, [entry[0] for entry in permuted[5]])

    def test_selected_hall_missing_or_duplicated_fails_closed(self):
        database, frame, _hall_rots, _hall_trans, sg_halls = (
            self._sg5_source_and_target())
        missing_frames = list(database.frames)
        missing_frames[4] = SimpleNamespace(
            spacegroup=5, data_hall=999, hall_operation_count=frame.hall_operation_count)
        missing_database = SimpleNamespace(frames=tuple(missing_frames))
        with self.assertRaisesRegex(ValueError, "legacy table matches"):
            generator._prepare_sidecar_hall_choices(missing_database, sg_halls)

        duplicated = {spacegroup: list(entries)
                      for spacegroup, entries in sg_halls.items()}
        selected = next(entry for entry in duplicated[5] if entry[0] == 9)
        duplicated[5].append(selected)
        with self.assertRaisesRegex(ValueError, "legacy table matches"):
            generator._prepare_sidecar_hall_choices(database, duplicated)

    def test_selected_rotation_and_sidecar_shift_mismatch_fail_closed(self):
        _database, frame, hall_rots, hall_trans, _sg_halls = (
            self._sg5_source_and_target())
        source_rots = hall_rots[:2]
        source_trans = hall_trans[:2]
        bad_rots = [list(rotation) for rotation in hall_rots]
        bad_rots[0][0] = 0
        with self.assertRaisesRegex(ValueError, "rotation mismatch"):
            generator._sidecar_source_hall_mapping(
                frame, 5, "GM1", source_rots, source_trans,
                bad_rots, hall_trans)

        bad_bindings = list(frame.hall_to_source)
        bad_bindings[0] = SimpleNamespace(
            hall_operation_index=0, source_operation_index=0,
            shift_numerator=(12, 0, 0))
        bad_frame = SimpleNamespace(
            source_operation_count=frame.source_operation_count,
            hall_operation_count=frame.hall_operation_count,
            hall_to_source=tuple(bad_bindings),
        )
        exact_target = generator._ExactScalarHallTarget(
            5, 9,
            tuple(binding.source_operation_index
                  for binding in frame.hall_to_source),
            tuple(tuple(binding.shift_numerator)
                  for binding in frame.hall_to_source),
            tuple(tuple(rotation) for rotation in hall_rots),
            tuple((0, 0, 0) for _ in hall_rots),
            tuple((0.0, 0.0, 0.0) for _ in hall_rots),
        )
        with self.assertRaisesRegex(ValueError, "exact Hall mapping mismatch"):
            generator._sidecar_source_hall_mapping(
                bad_frame, 5, "GM1", source_rots, source_trans,
                hall_rots, hall_trans, exact_target=exact_target)

    def test_compound_padding_uses_sidecar_source_to_hall(self):
        _database, frame, hall_rots, hall_trans, sg_halls = (
            self._sg5_source_and_target())
        choices = generator._prepare_sidecar_hall_choices(
            generator.load_committed_data_hall_provenance(), sg_halls)
        source_rots = hall_rots[:2]
        cir_rots = [value for rotation in source_rots for value in rotation]
        with mock.patch.object(
                generator, "_load_hall_operations",
                side_effect=AssertionError("padding must not search Hall candidates")):
            plans = generator._build_padding_plans(
                [5], ["compound"], [0], [1], [2], cir_rots, [None],
                sg_hall_choice=choices)
        self.assertEqual(plans, [(0, frame.hall_operation_count, [0, 1])])

    def test_scalar_selection_has_no_legacy_candidate_score(self):
        source = Path(generator.__file__).read_text(encoding="utf-8")
        self.assertNotIn("best_exact_count", source)
        self.assertNotIn("rot_cache", source)

    def test_decimal_phase_regression_does_not_use_exact_sidecar_shift(self):
        phased = generator._phase_character(
            1.0 + 0.0j,
            [0.3333333333, 0.0, 0.0],
            [1.0 / 3.0, 0.0, 0.0],
            (1, 0, 0, 2),
        )
        self.assertEqual(phased.imag, -1.0471976378421116e-10)


class Radical4CodebookTests(unittest.TestCase):
    def test_codebook_has_exact_pinned_spellings_and_tuples(self):
        self.assertEqual(
            frozenset(generator._PIR_RADICAL4_CODEBOOK),
            generator.PIR_MATRIX_TOKEN_SPELLINGS,
        )
        self.assertEqual(len(generator._PIR_RADICAL4_CODEBOOK), 25)
        self.assertEqual(
            generator._decode_pir_matrix_token("-0.96593"),
            generator.Radical4(0, -1, 0, -1),
        )
        self.assertEqual(
            generator._decode_pir_matrix_token("0.70711"),
            generator.Radical4(0, 2, 0, 0),
        )
        for spelling in ("+1", "1.00000", "0.70710", "unknown"):
            with self.assertRaises(ValueError):
                generator._decode_pir_matrix_token(spelling)

    def test_radical4_signs_and_exact_trace_cancellation(self):
        value = generator.Radical4(1, -2, 3, -4)
        self.assertTrue((-value + value).is_zero())
        self.assertEqual(
            generator._decode_pir_matrix_token("0.96593")
            + generator._decode_pir_matrix_token("0.25882"),
            generator.Radical4(0, 0, 0, 2),
        )
        self.assertAlmostEqual(
            (
                generator._decode_pir_matrix_token("0.96593")
                + generator._decode_pir_matrix_token("0.25882")
            ).materialize(),
            2**0.5 * 3**0.5 / 2,
        )
        self.assertEqual(
            generator._decode_pir_matrix_token("0.68301")
            + generator._decode_pir_matrix_token("0.18301"),
            generator.Radical4(0, 0, 2, 0),
        )
        self.assertTrue(
            (
                generator._decode_pir_matrix_token("0.25000")
                + generator._decode_pir_matrix_token("-0.25000")
            ).is_zero()
        )

    def test_scalar_formatter_is_shortest_ieee_roundtrip(self):
        for value in (1.047197638e-10, 1e-300, 0.0, 1.0, -1.0):
            literal = generator._format_scalar_roundtrip_f64(value)
            self.assertEqual(struct.pack(">d", float(literal)), struct.pack(">d", value))
        self.assertEqual(generator._format_scalar_roundtrip_f64(1.0), "1.0")
        with self.assertRaises(ValueError):
            generator._format_scalar_roundtrip_f64(-0.0)


class PirStructureTests(unittest.TestCase):
    @staticmethod
    def _header(
        irnumber=1,
        sg=2,
        space_group_symbol="P 1",
        ir_label="GM1",
        dim=3,
        irtype=1,
        kcount=1,
        pmkcount=1,
        opcount=1,
    ):
        return (
            f'  {irnumber}  {sg}  "{space_group_symbol}"  "{ir_label}"  '
            f"{dim}  {irtype}  {kcount}  {pmkcount}  {opcount}  "
        )

    @staticmethod
    def _operation_row():
        return "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"

    def test_exact_scalar_operation_decoder_normalizes_to_denominator12(self):
        row = [2, 0, 0, 1, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2]
        decoded = generator._decode_exact_scalar_operation(row, "PIR SG1 op0")
        self.assertEqual(
            decoded,
            generator._ExactScalarOperation(
                (1, 0, 0, 0, 1, 0, 0, 0, 1), (6, 0, 0)
            ),
        )

    def test_exact_scalar_operation_decoder_rejects_bad_rows(self):
        valid = [int(token) for token in self._operation_row().split()]
        cases = []
        cases.append((valid[:-1], "15 integers, expected 16"))
        nonexact = list(valid)
        nonexact[0] = True
        cases.append((nonexact, "non-exact integer"))
        bad_denominator = list(valid)
        bad_denominator[15] = 5
        cases.append((bad_denominator, "does not divide"))
        bad_rotation_division = list(valid)
        bad_rotation_division[0] = 1
        bad_rotation_division[15] = 2
        cases.append((bad_rotation_division, "not divisible"))
        bad_bottom = list(valid)
        bad_bottom[12] = 1
        cases.append((bad_bottom, "homogeneous bottom"))
        bad_determinant = list(valid)
        bad_determinant[10] = 0
        cases.append((bad_determinant, "determinant"))
        bad_domain = list(valid)
        bad_domain[0] = 2
        cases.append((bad_domain, "integer domain"))
        for row, message in cases:
            with self.assertRaisesRegex(ValueError, message):
                generator._decode_exact_scalar_operation(row, "synthetic op")

    @staticmethod
    def _kvector(nonzero=False):
        values = [0] * 16
        values[3] = 2
        if nonzero:
            values[4] = 1
        return values

    def test_kspecial_comes_from_augmented_kvector_not_label(self):
        self.assertTrue(generator._pir_kvector_is_special(self._kvector()))
        self.assertFalse(generator._pir_kvector_is_special(self._kvector(True)))

    def test_all_augmented_kvector_components_control_kspecial(self):
        for offset in (4, 5, 6, 8, 9, 10, 12, 13, 14):
            values = self._kvector()
            values[offset] = 1
            self.assertFalse(
                generator._pir_kvector_is_special(values),
                f"offset {offset} was not part of the kspecial test",
            )

    def test_multi_arm_kvector_block_exact_boundary(self):
        values = list(range(32))
        parsed, next_line = generator._read_exact_pir_int_block(
            [" ".join(map(str, values[:16])), " ".join(map(str, values[16:]))],
            0,
            32,
            "synthetic multi-arm record",
        )
        self.assertEqual(parsed, values)
        self.assertEqual(next_line, 2)

    def test_header_requires_exact_official_nine_fields(self):
        parsed = generator._parse_pir_header(self._header(), 17)
        self.assertEqual(parsed, (1, 2, "GM1", 3, 1, 1))
        malformed = (
            "EXTRA " + self._header(),
            self._header() + " EXTRA",
            '1 2 "P 1" "GM1" "THIRD" 3 1 1 1 1',
            '1 2 "P 1" "GM1" 3 1 1 1',
        )
        for line in malformed:
            with self.assertRaisesRegex(ValueError, "malformed PIR header"):
                generator._parse_pir_header(line, 23)

    def test_header_validates_every_integer_field(self):
        fields = (
            "irnumber",
            "sg",
            "dim",
            "irtype",
            "kcount",
            "pmkcount",
            "opcount",
        )
        for field in fields:
            for invalid in ("+1", "-1"):
                with self.assertRaisesRegex(ValueError, "malformed PIR header"):
                    generator._parse_pir_header(
                        self._header(**{field: invalid}), 29
                    )

    def test_header_semantic_ranges_and_kcount_relation(self):
        invalid_headers = (
            self._header(sg=0),
            self._header(sg=231),
            self._header(space_group_symbol=""),
            self._header(ir_label=""),
            self._header(dim=0),
            self._header(irtype=0),
            self._header(irtype=4),
            self._header(kcount=0),
            self._header(pmkcount=0),
            self._header(kcount=3, pmkcount=2),
            self._header(opcount=0),
        )
        for line in invalid_headers:
            with self.assertRaisesRegex(ValueError, "invalid PIR header field"):
                generator._parse_pir_header(line, 31)

        self.assertEqual(
            generator._parse_pir_header(
                self._header(space_group_symbol="P-1", ir_label="GM1+"), 37
            ),
            (1, 2, "GM1+", 3, 1, 1),
        )
        self.assertEqual(
            generator._parse_pir_header(
                self._header(space_group_symbol="P 1", ir_label="GM1-", kcount=2), 41
            ),
            (1, 2, "GM1-", 3, 1, 1),
        )

    def test_pir_irnumber_must_be_global_and_consecutive(self):
        for sequence in ((1, 3), (1, 1), (1, 2, 4)):
            expected = 1
            for line_number, actual in enumerate(sequence, start=1):
                if actual != expected:
                    with self.assertRaisesRegex(ValueError, "unexpected PIR irnumber"):
                        generator._require_pir_irnumber(
                            actual, expected, line_number, 2, "GM1"
                        )
                    break
                expected += 1
            else:
                self.fail(f"synthetic sequence {sequence!r} unexpectedly passed")

    def test_matrix_tokens_require_exact_archive_spellings(self):
        for token in ("2", "nan", "+1", "1.00000", "-0.5000", "malformed"):
            with self.assertRaisesRegex(ValueError, "unknown PIR matrix token"):
                generator._read_exact_pir_float_block(
                    [token], 0, 1, "synthetic matrix token"
                )

    def test_special_payload_has_no_irtranslation(self):
        payload = [self._operation_row(), "1"]
        _op, irtranslation, matrix, spellings, next_line = (
            generator._read_pir_operation_payload(
                payload, 0, 1, True, "label-that-does-not-matter"
            )
        )
        self.assertIsNone(irtranslation)
        self.assertEqual(matrix, [generator.Radical4(4, 0, 0, 0)])
        self.assertEqual(spellings, ["1"])
        self.assertEqual(next_line, 2)

    def test_nonspecial_payload_requires_irtranslation(self):
        payload = [self._operation_row(), "0 0 0 1", "1"]
        _op, irtranslation, matrix, _spellings, next_line = (
            generator._read_pir_operation_payload(
                payload, 0, 1, False, "GM-label-with-nonspecial-kvector"
            )
        )
        self.assertEqual(irtranslation, [0, 0, 0, 1])
        self.assertEqual(matrix, [generator.Radical4(4, 0, 0, 0)])
        self.assertEqual(next_line, 3)

    def test_malformed_pir_structure_is_rejected(self):
        operation = self._operation_row()
        cases = (
            [operation, "1"],  # missing irtranslation
            [operation, "0 0 0 1 9", "1"],  # extra irtranslation token
            [operation, "0 0 0 1"],  # truncated matrix
        )
        for payload in cases:
            with self.assertRaises(ValueError):
                generator._read_pir_operation_payload(
                    payload, 0, 1, False, "synthetic malformed record"
                )

    def test_archive_structural_census(self):
        parsed = generator._parse_pir_characters()
        source_records = parsed[-2]
        census = parsed[-1]
        self.assertEqual(len(source_records), 230)
        self.assertEqual(sum(len(record.operations) for record in source_records), 2609)
        self.assertEqual(census["records"], 10294)
        self.assertEqual(census["irtranslation_rows"], 64588)
        self.assertEqual(census["matrix_scalar_tokens"], 8977752)
        self.assertEqual(len(census["matrix_token_spellings"]), 25)
        self.assertEqual(
            census["matrix_token_spellings"], generator.PIR_MATRIX_TOKEN_SPELLINGS
        )


class CirStructureTests(unittest.TestCase):
    @staticmethod
    def _header(
        irnumber=1,
        sg=1,
        symbol="P1",
        label="GM1",
        dim=1,
        irtype=1,
        kcount=1,
        pmkcount=1,
        opcount=1,
    ):
        return (
            f' {irnumber} {sg} "{symbol}" "{label}" {dim} {irtype} '
            f"{kcount} {pmkcount} {opcount}"
        )

    @staticmethod
    def _kvector(special=True):
        values = [0] * 16
        values[3] = 1
        if not special:
            values[4] = 1
            values[7] = 2
        return " ".join(map(str, values))

    @staticmethod
    def _operation(denominator=1):
        values = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, denominator]
        return " ".join(map(str, values))

    @staticmethod
    def _mini_record(
        special=True,
        matrix="(1,0)",
        irnumber=1,
        label="GM1",
        dim=1,
        kcount=1,
        pmkcount=1,
        opcount=1,
    ):
        lines = ["title 1", "title 2", "title 3"]
        lines.append(CirStructureTests._header(
            irnumber=irnumber,
            label=label,
            dim=dim,
            kcount=kcount,
            pmkcount=pmkcount,
            opcount=opcount,
        ))
        lines.append(CirStructureTests._kvector(special=special))
        lines.append(CirStructureTests._operation())
        if not special:
            lines.append("0 0 0 1")
        lines.append(matrix)
        return lines

    def test_cir_header_requires_exact_fields_and_ascii_integers(self):
        valid = self._header(symbol="P-1", label="GM1+")
        self.assertEqual(generator._parse_cir_header(valid, 4)["label"], "GM1+")
        malformed = (
            "EXTRA " + valid,
            valid + " EXTRA",
            '1 1 "P1" "GM1" "THIRD" 1 1 1 1 1',
            '1 1 "P1" "GM1" 1 1 1 1',
            self._header(irnumber="+1"),
            self._header(sg="-1"),
            self._header(dim="+1"),
            self._header(irtype="-1"),
            self._header(kcount="+1"),
            self._header(pmkcount="-1"),
            self._header(opcount="+1"),
        )
        for line in malformed:
            with self.assertRaisesRegex(ValueError, "malformed CIR header"):
                generator._parse_cir_header(line, 4)

    def test_cir_header_validates_ranges_and_count_relation(self):
        invalid = (
            self._header(irnumber=0),
            self._header(sg=0),
            self._header(sg=231),
            self._header(symbol=""),
            self._header(label=""),
            self._header(dim=0),
            self._header(dim=49),
            self._header(irtype=0),
            self._header(irtype=4),
            self._header(kcount=0),
            self._header(kcount=49),
            self._header(pmkcount=0),
            self._header(opcount=0),
            self._header(kcount=3, pmkcount=2),
        )
        for line in invalid:
            with self.assertRaisesRegex(ValueError, "invalid CIR"):
                generator._parse_cir_header(line, 4)
        self.assertEqual(
            generator._parse_cir_header(
                self._header(symbol="P-1", label="GM1-", kcount=2), 4
            )["kcount"],
            2,
        )

    def test_cir_kvector_multiline_and_exact_boundary(self):
        first = " ".join(str(value) for value in range(16))
        second = " ".join(str(value) for value in range(16, 32))
        values, next_line = generator._read_exact_cir_block(
            [first, second], 0, 32, "synthetic CIR k-vector", generator._parse_cir_integer
        )
        self.assertEqual(values, list(range(32)))
        self.assertEqual(next_line, 2)
        with self.assertRaisesRegex(ValueError, "extra CIR tokens"):
            generator._read_exact_cir_block(
                [first + " 16"], 0, 16, "synthetic CIR k-vector", generator._parse_cir_integer
            )

    def test_cir_integer_tokens_are_canonical_ascii(self):
        self.assertEqual(generator._parse_cir_integer("-1", 1, "test"), -1)
        for token in ("+1", "01", "-0", "１"):
            with self.assertRaisesRegex(ValueError, "non-integer CIR token"):
                generator._parse_cir_integer(token, 1, "test")

        for line_index in (4, 5, 6):
            for bad_token in ("+1", "１"):
                lines = self._mini_record(special=False)
                values = lines[line_index].split()
                values[0] = bad_token
                lines[line_index] = " ".join(values)
                with self.assertRaisesRegex(ValueError, "non-integer CIR token"):
                    generator._parse_cir_lines(lines, validate_census=False)

        negative_kvector = self._mini_record()
        values = negative_kvector[4].split()
        values[0] = "-1"
        negative_kvector[4] = " ".join(values)
        _chars, _matrices, _source_records, _census = generator._parse_cir_lines(
            negative_kvector, validate_census=False
        )

    def test_cir_mini_records_are_structurally_consumed(self):
        chars, matrices, _source_records, census = generator._parse_cir_lines(
            self._mini_record(special=True), validate_census=False
        )
        self.assertEqual(chars[(1, "GM1")]["chars"], [(1.0, 0.0, 1.0)])
        self.assertEqual(matrices[(1, "GM1")], [(1.0, 0.0)])
        self.assertEqual(census["cursor_eof"], 7)
        chars, _matrices, _source_records, census = generator._parse_cir_lines(
            self._mini_record(special=False), validate_census=False
        )
        self.assertEqual(chars[(1, "GM1")]["chars"][0][0], 1.0)
        self.assertEqual(census["irtranslation_rows"], 1)

    def test_cir_rejects_payload_truncation_and_extra_tokens(self):
        operation = self._operation()
        cases = (
            self._mini_record(special=False)[:-2] + ["0 0 0 1"],
            self._mini_record(special=True) + ["EXTRA"],
            self._mini_record(special=True)[:-1] + ["(1,0) (0,0)"],
            self._mini_record(special=False)[:5] + ["0 0 0 1 9", "(1,0)"],
            self._mini_record(special=True)[:5] + [operation, "(1,0)junk"],
        )
        for lines in cases:
            with self.assertRaises(ValueError):
                generator._parse_cir_lines(lines, validate_census=False)

    def test_cir_rejects_bad_operation_and_translation_rows(self):
        for operation in (
            "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0",
            self._operation() + " 0",
            "not-an-operation",
            self._operation(denominator=0),
        ):
            lines = self._mini_record()[:5] + [operation, "(1,0)"]
            with self.assertRaises(ValueError):
                generator._parse_cir_lines(lines, validate_census=False)
        lines = self._mini_record(special=False)
        lines[6] = "0 0 0 0"
        with self.assertRaises(ValueError):
            generator._parse_cir_lines(lines, validate_census=False)

        bottom_bad = self._mini_record()
        bottom = bottom_bad[5].split()
        bottom[12] = "1"
        bottom_bad[5] = " ".join(bottom)
        with self.assertRaisesRegex(ValueError, "bottom row"):
            generator._parse_cir_lines(bottom_bad, validate_census=False)

        nondivisible = self._mini_record()
        operation = nondivisible[5].split()
        operation[0] = "1"
        operation[15] = "2"
        nondivisible[5] = " ".join(operation)
        with self.assertRaisesRegex(ValueError, "not divisible"):
            generator._parse_cir_lines(nondivisible, validate_census=False)

    def test_cir_kvector_denominator_rules_are_explicit(self):
        values = [0] * 16
        values[3] = 1
        values[4] = 1
        values[7] = 0
        lines = self._mini_record(special=True)
        lines[4] = " ".join(map(str, values))
        with self.assertRaisesRegex(ValueError, "zero CIR parameter denominator"):
            generator._parse_cir_lines(lines, validate_census=False)

    def test_cir_nontrivial_little_dim_and_needed_labels_consume_all_input(self):
        dim_two = self._mini_record(
            irnumber=2,
            label="GM2",
            dim=2,
            matrix="(1,0) (0,0) (0,0) (1,0)",
        )
        lines = self._mini_record(label="GM1")[3:] + dim_two[3:]
        chars, matrices, _source_records, census = generator._parse_cir_lines(
            ["title 1", "title 2", "title 3"] + lines,
            needed_labels={(1, "GM2")},
            validate_census=False,
        )
        self.assertEqual(chars[(1, "GM2")]["little_dim"], 2)
        self.assertNotIn((1, "GM1"), matrices)
        self.assertIn((1, "GM2"), matrices)
        self.assertEqual(census["cursor_eof"], len(lines) + 3)

    def test_cir_trailing_blank_or_nonheader_is_rejected(self):
        for trailing in ("", "trailing garbage"):
            with self.assertRaises(ValueError):
                generator._parse_cir_lines(
                    self._mini_record() + [trailing], validate_census=False
                )

    def test_cir_rejects_invalid_complex_spellings_without_silent_zero(self):
        for token in ("CORRUPTED", "1", "(nan,0)", "(1,)", "(1,0)junk", "(1.00000,0)"):
            with self.assertRaises(ValueError):
                generator._parse_cir_lines(
                    self._mini_record(matrix=token), validate_census=False
                )

    def test_cir_complex_components_use_the_pir_radical_decoder(self):
        self.assertEqual(len(generator.CIR_COMPLEX_TOKEN_SPELLINGS), 65)
        real, imag = generator._parse_complex("(0.96593,0.25882)")
        self.assertEqual(real, generator.Radical4(0, 1, 0, 1))
        self.assertEqual(imag, generator.Radical4(0, -1, 0, 1))
        self.assertEqual(
            real + imag,
            generator.Radical4(0, 0, 0, 2),
        )

    def test_cir_irnumber_sequence_is_global(self):
        for sequence in ((1, 3), (1, 1), (1, 2, 4)):
            lines = ["title 1", "title 2", "title 3"]
            for number in sequence:
                lines.extend(self._mini_record(irnumber=number, label=f"GM{number}")[3:])
            with self.assertRaisesRegex(ValueError, "unexpected CIR irnumber"):
                generator._parse_cir_lines(lines, validate_census=False)

    def test_cir_archive_structural_census(self):
        _chars, _matrices, source_records, census = generator._parse_cir_lines(
            generator._read_cir_lines())
        self.assertEqual(len(source_records), 230)
        self.assertEqual(sum(len(record.operations) for record in source_records), 2609)
        self.assertEqual(census, {
            "records": 11202,
            "kvector_ints": 555920,
            "operation_rows": 133246,
            "irtranslation_rows": 68612,
            "complex_tokens": 7121956,
            "complex_token_spellings": generator.CIR_COMPLEX_TOKEN_SPELLINGS,
            "cursor_eof": 877084,
            "irtype_counts": {1: 7796, 2: 155, 3: 3251},
            "kcount_ratio_counts": {1: 6252, 2: 4950},
        })


class ExactScalarProvenanceTests(unittest.TestCase):
    @staticmethod
    def _operation(rotation, translation=(0, 0, 0)):
        return generator._ExactScalarOperation(tuple(rotation), tuple(translation))

    @classmethod
    def _two_operations(cls):
        identity = cls._operation((1, 0, 0, 0, 1, 0, 0, 0, 1))
        inversion = cls._operation((-1, 0, 0, 0, 1, 0, 0, 0, -1))
        return identity, inversion

    def test_archive_snapshot_rejects_order_subset_and_duplicate_tables(self):
        identity, inversion = self._two_operations()
        for archive in ("PIR", "CIR"):
            source_operations = {}
            source_anchors = {}
            generator._record_exact_scalar_archive_operations(
                archive, source_operations, source_anchors, 1, 7,
                (identity, inversion))
            for operations in ((inversion, identity), (identity,)):
                with self.assertRaisesRegex(ValueError, "order/table differs"):
                    generator._record_exact_scalar_archive_operations(
                        archive, source_operations, source_anchors, 1, 8, operations)
            with self.assertRaisesRegex(ValueError, "duplicates"):
                generator._record_exact_scalar_archive_operations(
                    archive, {}, {}, 1, 8, (identity, identity))

    def test_pir_cir_snapshot_mismatch_is_rejected(self):
        identity, inversion = self._two_operations()
        pir = tuple(
            generator._ExactScalarArchiveRecord(sg, sg, (identity,))
            for sg in range(1, 231))
        cir_records = list(pir)
        cir_records[4] = generator._ExactScalarArchiveRecord(5, 5, (inversion,))
        with self.assertRaisesRegex(ValueError, "PIR/CIR SG5 source operation order"):
            generator._merge_exact_scalar_source_frames(pir, tuple(cir_records))

    def test_exact_target_bridge_is_bitwise_and_has_no_modulo(self):
        identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        target = generator._ExactScalarHallTarget(
            1, 1, (0,), ((0, 0, 0),), (identity,), ((4, 0, 0),),
            ((float(1) / 3.0, 0.0, 0.0),))
        rounded = generator._round_exact_translation_to_10_decimal(4)
        self.assertTrue(generator._same_f64_bits(rounded, float("0.3333333333")))
        generator._validate_exact_legacy_hall_bridge(
            target, [identity], [[rounded, 0.0, 0.0]], "synthetic")
        with self.assertRaisesRegex(ValueError, "fixed ten-decimal"):
            generator._validate_exact_legacy_hall_bridge(
                target, [identity], [[0.3333333334, 0.0, 0.0]], "synthetic")
        with self.assertRaisesRegex(ValueError, "0..11"):
            generator._round_exact_translation_to_10_decimal(12)
        with self.assertRaisesRegex(ValueError, "0..11"):
            generator._round_exact_translation_to_10_decimal(-1)

    def test_sg5_exact_mapping_direction_and_source_float_bridge(self):
        database = generator.load_committed_data_hall_provenance()
        frame = database.frames[4]
        hall_rots, hall_trans = next(
            (entry[1], entry[2]) for entry in generator._load_hall_operations()[5]
            if entry[0] == 9)
        identity, inversion = self._two_operations()
        source_frame = generator._ExactScalarSourceFrame(
            5, frame.pir_anchor_irnumber, frame.cir_anchor_irnumber,
            (identity, inversion))
        h2s = frame.hall_to_source
        exact_target = generator._ExactScalarHallTarget(
            5, 9,
            tuple(binding.source_operation_index for binding in h2s),
            tuple(tuple(binding.shift_numerator) for binding in h2s),
            tuple(tuple(rotation) for rotation in hall_rots),
            tuple(tuple(binding.shift_numerator) for binding in h2s),
            tuple(tuple(float(value) / 12.0 for value in binding.shift_numerator)
                  for binding in h2s),
        )
        self.assertEqual(exact_target.hall_to_source[:4], (0, 1, 0, 1))
        self.assertEqual(exact_target.shift_numerators[2], (6, 6, 0))
        source_rots = [hall_rots[0], hall_rots[1]]
        source_trans = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        mapping = generator._sidecar_source_hall_mapping(
            frame, 5, "GM1", source_rots, source_trans,
            hall_rots, hall_trans, source_frame, exact_target)
        self.assertEqual(mapping, [0, 1, 0, 1])
        bad_source_trans = [[0.0, 0.0, 0.0], [1e-12, 0.0, 0.0]]
        with self.assertRaisesRegex(ValueError, "exact /12 binary64 bridge"):
            generator._sidecar_source_hall_mapping(
                frame, 5, "GM1", source_rots, bad_source_trans,
                hall_rots, hall_trans, source_frame, exact_target)


if __name__ == "__main__":
    unittest.main()
