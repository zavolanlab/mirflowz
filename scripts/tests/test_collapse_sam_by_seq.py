"""Unit tests for module 'collapse_sam_by_seq.py'."""

import argparse
import sys

from pathlib import Path
import pysam
import pytest

sys.path.append("../../")

from scripts.collapse_sam_by_seq import (
    collapse_alignments,
    main,
    parse_arguments,
)


@pytest.fixture
def sam_empty_file():
    empty_file = Path("files/header_only.sam")
    
    return empty_file


@pytest.fixture
def sam_files():
    in_sam = Path("files/in_alignments.sam")
    out_sam = Path("files/collapsed_alignments.sam")

    return in_sam, out_sam


@pytest.fixture
def sam_no_collapse():
    in_sam = Path("files/in_alignments_missing.sam")
    out_sam = Path("files/no_collapsed_alignments.sam")

    return in_sam, out_sam

@pytest.fixture
def sam_no_NH():
    in_sam = Path("files/in_alignments_no_NH.sam")
    out_sam = Path("files/no_NH_collapsed_alignments.sam")

    return in_sam, out_sam


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without input file."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['collapse_sam_by_seq']
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, sam_empty_file):
        """Call with a single input file."""
        sam = sam_empty_file
        monkeypatch.setattr(
            sys, 'argv',
            ['collapse_sam_by_seq',
             str(sam),
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_too_many_inputs(self, monkeypatch, sam_empty_file):
        """Call with too many input files."""
        sam = sam_empty_file
        monkeypatch.setattr(
            sys, 'argv',
            ['collapse_sam_by_seq',
             str(sam), str(sam),
             ]
        )
        with pytest.raises(SystemExit) as sysex:
            parse_arguments().parse_args()
        assert sysex.value.code == 2


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_input(self, monkeypatch, capsys, sam_empty_file):
        """Test main function with empty sam file."""
        in_sam = sam_empty_file

        monkeypatch.setattr(
            sys, 'argv',
            ['collapse_sam_by_seq',
             str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args.infile)
        captured = capsys.readouterr()

        with open(in_sam, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_complete_file(self, monkeypatch, capsys, sam_files):
        """Test main function with no extension in features names."""
        in_sam, expected_out = sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['collapse_sam_by_seq',
             str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args.infile)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_no_collapse(self, monkeypatch, capsys, sam_no_collapse):
        """Test main function with no alignments to collapse."""
        in_sam, expected_out = sam_no_collapse

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args.infile)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_no_NH(self, monkeypatch, capsys, sam_no_NH):
        """Test main function with no NH tag in some alignments."""
        in_sam, expected_out = sam_no_NH

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args.infile)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()
