"""Unit tests for module 'iso_name_tagging.py'."""

import argparse
from pathlib import Path
import sys

import pytest

from ..iso_name_tagging import main, parse_arguments


@pytest.fixture
def empty_files():
    """Import path to empty files."""
    empty_intersect = Path("files/empty_file")
    empty_sam = Path("files/header_only.sam")

    return empty_intersect, empty_sam


@pytest.fixture
def intersect_sam():
    """Import path to INTERSECT and SAM files."""
    intersect_file = Path("files/in_intersection_mirna.bed")
    sam_file = Path("files/in_alignments_mirna.sam")
    output_file = Path("files/mirna_tag.sam")

    return intersect_file, sam_file, output_file


@pytest.fixture
def intersect_sam_extension():
    """Import path to INTERSECT and SAM files with miRNA extension."""
    intersect_file = Path("files/in_intersection_extended_mirna.bed")
    sam_file = Path("files/in_alignments_mirna.sam")
    output_file = Path("files/mirna_tag_extension.sam")

    return intersect_file, sam_file, output_file


@pytest.fixture
def intersect_sam_id():
    """Import path to INTERSECT and SAM files with miRNA IDs in the output."""
    intersect_file = Path("files/in_intersection_mirna.bed")
    sam_file = Path("files/in_alignments_mirna.sam")
    output_file = Path("files/mirna_tag_id.sam")

    return intersect_file, sam_file, output_file


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_intersect(self, monkeypatch, intersect_sam):
        """Call without intersect file."""
        in_intersect, in_sam, output = intersect_sam

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "iso_name_tagging",
                    "--sam",
                    str(in_sam),
                ],
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_no_sam(self, monkeypatch, intersect_sam):
        """Call without intersect file."""
        in_intersect, in_sam, output = intersect_sam

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "iso_name_tagging",
                    "--intersect",
                    str(in_intersect),
                ],
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, intersect_sam):
        """Call with the correct input files."""
        in_intersect, in_sam, output = intersect_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(in_intersect),
                "--sam",
                str(in_sam),
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_all_input(self, monkeypatch, intersect_sam):
        """Call with all the options."""
        in_intersect, in_sam, output = intersect_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(in_intersect),
                "--sam",
                str(in_sam),
                "--id",
                "alias",
                "--extension",
                "6",
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_intersect_file(
        self, monkeypatch, capsys, empty_files, intersect_sam
    ):
        """Test main function with an empty intersect file."""
        empty_intersect, empty_sam = empty_files

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(empty_intersect),
                "--sam",
                str(empty_sam),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_sam, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_empty_sam_file(
        self, monkeypatch, capsys, empty_files, intersect_sam
    ):
        """Test main function with an empty sam file."""
        empty_intersect, empty_sam = empty_files
        in_intersect, in_sam, output = intersect_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(in_intersect),
                "--sam",
                str(empty_sam),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_sam, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_intersect_sam_file(self, monkeypatch, capsys, intersect_sam):
        """Test main function without options."""
        in_intersect, in_sam, output = intersect_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(in_intersect),
                "--sam",
                str(in_sam),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_intersect_sam_extension_file(
        self, monkeypatch, capsys, intersect_sam_extension
    ):
        """Test main function with extension equals 6."""
        in_intersect, in_sam, output = intersect_sam_extension

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(in_intersect),
                "--sam",
                str(in_sam),
                "--extension",
                "6",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, "r") as out_file:
            assert captured.out == out_file.read()

    def test_main_intersect_sam_file_id(
        self, monkeypatch, capsys, intersect_sam_id
    ):
        """Test main function with id equals id."""
        in_intersect, in_sam, output = intersect_sam_id

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--intersect",
                str(in_intersect),
                "--sam",
                str(in_sam),
                "--id",
                "id",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, "r") as out_file:
            assert captured.out == out_file.read()
