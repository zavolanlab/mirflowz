"""Unit tests for module 'iso_name_tagging.py'."""

import argparse
from pathlib import Path
import sys

import pytest

sys.path.append("../../")

from scripts.iso_name_tagging import (
    main,
    parse_arguments
)


@pytest.fixture
def empty_files():
    """Import path to empty files."""
    empty_bed = Path("files/empty_file")
    empty_sam = Path("files/header_only.sam")

    return empty_bed, empty_sam


@pytest.fixture
def bed_sam():
    """Import path to BED and SAM files."""
    bed_file = Path("files/in_intersection_mirna.bed")
    sam_file = Path("files/in_alignments_mirna.sam")
    output_file = Path("files/mirna_tag.sam")

    return bed_file, sam_file, output_file


@pytest.fixture
def bed_sam_extension():
    """Import path to BED and SAM files with miRNA extension."""
    bed_file = Path("files/in_intersection_mirna.bed")
    sam_file = Path("files/in_alignments_mirna.sam")
    output_file = Path("files/mirna_tag_extension.sam")

    return bed_file, sam_file, output_file


@pytest.fixture
def bed_sam_id():
    """Import path to BED and SAM files with miRNA IDs in the output."""
    bed_file = Path("files/in_intersection_mirna.bed")
    sam_file = Path("files/in_alignments_mirna.sam")
    output_file = Path("files/mirna_tag_id.sam")

    return bed_file, sam_file, output_file


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_bed(self, monkeypatch, bed_sam):
        """Call without bed file."""
        in_bed, in_sam, output = bed_sam

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['iso_name_tagging',
                 '--sam', str(in_sam),
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2
    
    def test_no_sam(self, monkeypatch, bed_sam):
        """Call without bed file."""
        in_bed, in_sam, output = bed_sam

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['iso_name_tagging',
                 '--bed', str(in_bed),
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, bed_sam):
        """Call with the correct input files."""
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(in_bed),
             '--sam', str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)
    
    def test_all_input(self, monkeypatch, bed_sam):
        """Call with all the options."""
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(in_bed),
             '--sam', str(in_sam),
             '--id', "alias",
             '--extension', '6',
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_bed_file(self, monkeypatch, capsys, empty_files, bed_sam):
        """Test main function with an empty bed file."""
        empty_bed, empty_sam = empty_files

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(empty_bed),
             '--sam', str(empty_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_sam, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_empty_sam_file(self, monkeypatch, capsys, empty_files, bed_sam):
        """Test main function with an empty sam file."""
        empty_bed, empty_sam = empty_files
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(in_bed),
             '--sam', str(empty_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_sam, 'r') as out_file:
            assert captured.out == out_file.read()
    

    def test_main_bed_sam_file(self, monkeypatch, capsys, bed_sam):
        """Test main function without options."""
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(in_bed),
             '--sam', str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, 'r') as out_file:
            assert captured.out == out_file.read()
    
    def test_main_bed_sam_extension_file(self, monkeypatch, capsys, bed_sam_extension):
        """Test main function with extension equals 6."""
        in_bed, in_sam, output = bed_sam_extension

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(in_bed),
             '--sam', str(in_sam),
             '--extension', '6',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_bed_sam_file(self, monkeypatch, capsys, bed_sam_id):
        """Test main function with id equals id."""
        in_bed, in_sam, output = bed_sam_id

        monkeypatch.setattr(
            sys, 'argv',
            ['iso_name_tagging',
             '--bed', str(in_bed),
             '--sam', str(in_sam),
             '--id', 'id'
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, 'r') as out_file:
            assert captured.out == out_file.read()