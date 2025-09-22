"""Unit tests for module 'primir_quantification.py'."""

import argparse
from pathlib import Path
import sys

import pytest

from ..primir_quantification import main, parse_arguments
from ..validate_bedtools_intersect import FileFormatError


@pytest.fixture
def empty_file():
    """Import path to empty file."""
    empty_file = Path("files/empty_file")

    return empty_file


@pytest.fixture
def intersect_file():
    """Import path to test files with full content."""
    in_intersect = Path("files/in_intersection.intersect")
    out_table = Path("files/primir_quantification")

    return in_intersect, out_table


@pytest.fixture
def intersect_no_extension_files():
    """Import path to test files with no extension on features names."""
    in_intersect = Path("files/in_intersection_no_extension.intersect")
    out_table = Path("files/no_extension_primir_quantification")

    return in_intersect, out_table


@pytest.fixture
def intersect_extension_id_files():
    """Import path to test files with extension and reads id."""
    in_intersect = Path("files/in_intersection.intersect")
    out_table = Path("files/extension_id_primir_quantification")

    return in_intersect, out_table


@pytest.fixture
def intersect_id_files():
    """Import path to test files with reads id."""
    in_intersect = Path("files/in_intersection_no_extension.intersect")
    out_table = Path("files/id_primir_quantification")

    return in_intersect, out_table


@pytest.fixture
def intersect_some_extension_files():
    """Import path to test files with some extension."""
    in_intersect = Path("files/in_intersection_some_extension.intersect")
    out_table = Path("files/some_extension_primir_quantification")

    return in_intersect, out_table


@pytest.fixture
def intersect_collapsed_file():
    """Import path to test files with full content."""
    in_intersect = Path("files/in_intersection_collapsed.intersect")
    out_table = Path("files/collapsed_primir_quantification")

    return in_intersect, out_table


@pytest.fixture
def intersect_nh_file():
    """Import path to test files with full content."""
    in_intersect = Path("files/in_intersection.intersect")
    out_table = Path("files/nh_primir_quantification")

    return in_intersect, out_table


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without input file."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "primir_quantification",
                ],
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, intersect_file):
        """Call with the correct input file."""
        in_intersect, out_table = intersect_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_too_many_input_files(self, monkeypatch, intersect_file):
        """Call with too many input file."""
        in_intersect, out_table = intersect_file

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "primir_quantification",
                    str(in_intersect),
                    str(in_intersect),
                ],
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_all_input(self, monkeypatch, intersect_file):
        """Call with all the options."""
        in_intersect, out_table = intersect_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--id",
                "name",
                "--feat-extension",
                "--read-ids",
                "--collapsed",
                "--nh",
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_intersect_file(self, monkeypatch, capsys, empty_file):
        """Test main function with an empty intersect file."""
        empty_file = empty_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(empty_file),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_file, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_no_extension(
        self, monkeypatch, capsys, intersect_no_extension_files
    ):
        """Test main function with no extension in features names."""
        in_intersect, expected_out = intersect_no_extension_files

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_id_extension(
        self, monkeypatch, capsys, intersect_extension_id_files
    ):
        """Test main function with extension in feature name and read names."""
        in_intersect, expected_out = intersect_extension_id_files

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--feat-extension",
                "--read-ids",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_id(self, monkeypatch, capsys, intersect_id_files):
        """Test main function with read names."""
        in_intersect, expected_out = intersect_id_files

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--read-ids",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_some_extension_file(
        self, monkeypatch, capsys, intersect_some_extension_files
    ):
        """Test main function with read names."""
        in_intersect, expected_out = intersect_some_extension_files

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--feat-extension",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_collpased_nh_file(self, monkeypatch, capsys, intersect_file):
        """Test main function with collapsed alignments and nh value."""
        in_intersect, expected_out = intersect_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--collapsed",
                "--nh",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_collpased_file(
        self, monkeypatch, capsys, intersect_collapsed_file
    ):
        """Test main function with collapsed alignments."""
        in_intersect, expected_out = intersect_collapsed_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--collapsed",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_nh_file(self, monkeypatch, capsys, intersect_nh_file):
        """Test main function with nh value."""
        in_intersect, expected_out = intersect_nh_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_intersect),
                "--nh",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_invalid_file_validation(self, monkeypatch, capsys, tmp_path):
        """Test main function with invalid file."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(FileFormatError, match=r".*16 columns, found 12."):
            main(args)

    def test_main_malformed_read_nh_collapsed_no_delim(
            self,
            monkeypatch,
            capsys,
            tmp_path
    ):
        """Test main function with malformed read names."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\t"
            "read:without:delimiter\t255\t+\t35\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
                "--collapsed",
                "--nh"
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(Exception, match=r".* READ-#reads_NH."):
            main(args)

    def test_main_malformed_read_nh_collapsed_delim(
            self,
            monkeypatch,
            capsys,
            tmp_path
    ):
        """Test main function with malformed read names."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\t"
            "read-with_delimiter\t255\t+\t35\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
                "--collapsed",
                "--nh"
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(Exception, match=r".* READ-#reads_NH."):
            main(args)

    def test_main_malformed_read_nh_no_delim(
            self,
            monkeypatch,
            capsys,
            tmp_path
    ):
        """Test main function with malformed read names."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\t"
            "read:without:delimiter\t255\t+\t35\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
                "--nh"
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(Exception, match=r".* READ_NH."):
            main(args)

    def test_main_malformed_read_nh_delim(
            self,
            monkeypatch,
            capsys,
            tmp_path
    ):
        """Test main function with malformed read names."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\t"
            "read_delimiter\t255\t+\t35\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
                "--nh",
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(Exception, match=r".* READ_NH."):
            main(args)

    def test_main_malformed_read_collapsed_no_delim(
            self,
            monkeypatch,
            capsys,
            tmp_path
    ):
        """Test main function with malformed read names."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\t"
            "read_delimiter\t255\t+\t35\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
                "--collapsed",
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(Exception, match=r".* READ-#reads."):
            main(args)

    def test_main_malformed_read_collapsed_delim(
            self,
            monkeypatch,
            capsys,
            tmp_path
    ):
        """Test main function with malformed read names."""
        malformed_content = (
            "chr1\t.\tmiRNA_primary_transcript\t100\t300\t.\t+\t.\t"
            "ID=MI0003786;Name=hsa-mir-1323_-0_+0\tchr1\t130\t165\t"
            "read-delimiter\t255\t+\t35\n"
        )
        in_file = tmp_path / " malformed.intersect"
        in_file.write_text(malformed_content)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "primir_quantification",
                str(in_file),
                "--collapsed",
            ],
        )

        args = parse_arguments().parse_args()

        with pytest.raises(Exception, match=r".* READ-#reads."):
            main(args)

