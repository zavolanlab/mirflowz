"""Unit tests for module 'nh_filter.py'."""

import argparse
from pathlib import Path
import sys

import pytest

from ..nh_filter import (
    main,
    parse_arguments,
)


@pytest.fixture
def empty_file():
    """Import path to empty file."""
    empty_file = Path("files/header_only.sam")

    return empty_file


@pytest.fixture
def sam_file():
    """Import path to test files with full content."""
    in_sam = Path("files/in_aln_tag.sam")
    out_sam = Path("files/aln_nh_2_filtered.sam")

    return in_sam, out_sam


@pytest.fixture
def nh_missing_sam_file():
    """Import path to test files with alignments missing the NH tag."""
    sam_file = Path("files/in_aln_tag_missing_nh.sam")

    return sam_file


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_in_sam(self, monkeypatch):
        """Call without the positional argument."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "nh_filter",
                    "--out_file",
                    str(sam_file),
                ],
            )
            parse_arguments().parse_args()
            assert sysex.value.code == 2

    def test_no_out_sam(self, monkeypatch, empty_file):
        """Call without the output file path."""
        sam_file = empty_file

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "nh_filter",
                    str(sam_file),
                ],
            )
            parse_arguments().parse_args()
            assert sysex.value.code == 2

    def test_correct_args(self, monkeypatch, empty_file):
        """Call with all the correct arguments."""
        sam_file = empty_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "nh_filter",
                str(sam_file),
                "--out_file",
                str(sam_file),
                "--max_nh",
                "100",
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_sam_file(self, monkeypatch, tmp_path, empty_file):
        """Test main function with an empty SAM file."""
        empty_file = empty_file
        output = tmp_path / "output.sam"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "nh_filter",
                str(empty_file),
                "--out_file",
                str(output),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(empty_file, "r") as expected, open(output, "r") as out_file:
            assert out_file.read() == expected.read()

    def test_main_correct_sam_file(
        self, monkeypatch, tmp_path, sam_file
    ):
        """Test main function with correct SAM file and maximum NH 2."""
        in_sam, out_sam = sam_file
        output = tmp_path / "output.sam"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "nh_filter",
                str(in_sam),
                "--out_file",
                str(output),
                "--max_nh",
                "2",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_sam, "r") as expected, open(output, "r") as out_file:
            assert out_file.read() == expected.read()

    def test_main_nh_missing_sam_file(
        self, monkeypatch, nh_missing_sam_file
    ):
        """Test main function with some missing NH tag in SAM file."""
        infile = nh_missing_sam_file

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "nh_filter",
                str(infile),
                "--out_file",
                str(infile),
            ],
        )
        args = parse_arguments().parse_args()
        with pytest.raises(KeyError, match=r".* associated NH .*"):
            main(args)
