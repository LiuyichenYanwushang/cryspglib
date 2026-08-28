"""Focused tests for the pinned generator archive boundary."""

import os
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(__file__))
import generate_irrep_data as generator


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

    def test_ambiguous_suffix_member_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            self._write_zip(directory, {
                "first/data.txt": "first",
                "second/data.txt": "second",
            })
            generator.ISO_DIR = directory
            with self.assertRaisesRegex(ValueError, "ambiguous archive member"):
                generator._open_zip_path("iso.zip", "data.txt")

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
        self.assertEqual(matrix, [1.0])
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
        self.assertEqual(matrix, [1.0])
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
        census = parsed[-1]
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
        _chars, _matrices, _census = generator._parse_cir_lines(
            negative_kvector, validate_census=False
        )

    def test_cir_mini_records_are_structurally_consumed(self):
        chars, matrices, census = generator._parse_cir_lines(
            self._mini_record(special=True), validate_census=False
        )
        self.assertEqual(chars[(1, "GM1")]["chars"], [(1.0, 0.0, 1.0)])
        self.assertEqual(matrices[(1, "GM1")], [(1.0, 0.0)])
        self.assertEqual(census["cursor_eof"], 7)
        chars, _matrices, census = generator._parse_cir_lines(
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
        chars, matrices, census = generator._parse_cir_lines(
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

    def test_cir_irnumber_sequence_is_global(self):
        for sequence in ((1, 3), (1, 1), (1, 2, 4)):
            lines = ["title 1", "title 2", "title 3"]
            for number in sequence:
                lines.extend(self._mini_record(irnumber=number, label=f"GM{number}")[3:])
            with self.assertRaisesRegex(ValueError, "unexpected CIR irnumber"):
                generator._parse_cir_lines(lines, validate_census=False)

    def test_cir_archive_structural_census(self):
        _chars, _matrices, census = generator._parse_cir_lines(generator._read_cir_lines())
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


if __name__ == "__main__":
    unittest.main()
