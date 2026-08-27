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
        self.assertEqual(census["unmatched_structural_tokens"], 0)


if __name__ == "__main__":
    unittest.main()
