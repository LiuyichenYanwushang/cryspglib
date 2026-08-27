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


if __name__ == "__main__":
    unittest.main()
