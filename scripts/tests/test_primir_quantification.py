"""Unit tests for module 'primir_quantification.py'."""

import argparse
from pathlib import Path
import sys

import pytest

sys.path.append("../../")

from scripts.primir_quantification import (
    main,
    parse_arguments
)


@pytest.fixture
def empty_file():
    """Import path to empty file."""
    empty_file = Path("files/empty_file")
    sam_empty_file = Path("files/header_only.sam")

    return empty_file, sam_empty_file


@pytest.fixture
def bed_sam_files():
    """Import path to test BAM/SAM files with full content."""
    in_bed = Path("files/in_intersection.bed")
    in_sam = Path("files/in_alignments.sam")
    out_table = Path("files/sam_primir_quantification")

    return in_bed, in_sam, out_table


@pytest.fixture
def bed_sam_missing_aln_files():
    """Import path to test files with missing alignments."""
    in_bed = Path("files/in_intersection.bed")
    in_sam = Path("files/in_alignments_missing.sam")
    out_table = Path("files/sam_missing_alns_primir_quantification")

    return in_bed, in_sam, out_table


@pytest.fixture
def bed_sam_no_extension_files():
    """Import path to test files with no extension on features names."""
    in_bed = Path("files/in_intersection_no_extension.bed")
    in_sam = Path("files/in_alignments.sam")
    out_table = Path("files/no_extension_primir_quantification")

    return in_bed, in_sam, out_table


@pytest.fixture
def bed_extension_id_files():
    """Import path to test files with extension and reads id."""
    in_bed = Path("files/in_intersection.bed")
    in_sam = Path("files/in_alignments.sam")
    out_table = Path("files/extension_id_primir_quantification")

    return in_bed, in_sam, out_table


@pytest.fixture
def bed_id_files():
    """Import path to test files with reads id."""
    in_bed = Path("files/in_intersection_no_extension.bed")
    in_sam = Path("files/in_alignments.sam")
    out_table = Path("files/id_primir_quantification")

    return in_bed, in_sam, out_table


@pytest.fixture
def bed_id_missing_aln_files():
    """Import path to test files with missing alignments and reads id."""
    in_bed = Path("files/in_intersection_no_extension.bed")
    in_sam = Path("files/in_alignments_missing.sam")
    out_table = Path("files/missing_aln_id_primir_quantification")

    return in_bed, in_sam, out_table


@pytest.fixture
def bed_no_sam_files():
    """Import path to test files with no sam file."""
    in_bed = Path("files/in_intersection_no_extension.bed")
    out_table = Path("files/no_sam_primir_quantification")

    return in_bed, out_table
    

@pytest.fixture
def bed_some_extension_files():
    """Import path to test files with some extension."""
    in_bed = Path("files/in_intersection_some_extension.bed")
    out_table = Path("files/some_extension_primir_quantification")

    return in_bed, out_table


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without input file."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, bed_sam_files):
        """Call with the correct input file."""
        in_bed, in_sam, out_table = bed_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)
    
    def test_too_many_input_files(self, monkeypatch, bed_sam_files):
        """Call with too many input file."""
        in_bed, in_sam, out_table = bed_sam_files
        
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                 str(in_bed), str(in_bed),
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2
    
    def test_all_input(self, monkeypatch, bed_sam_files):
        """Call with the correct input file."""
        in_bed, in_sam, out_table = bed_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             '--id', "name",
             '--feat-extension',
             '--read-ids',
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_bed_file(self, monkeypatch, capsys, empty_file, bed_sam_files):
        """Test main function with an empty bed file."""
        empty_file, empty_sam = empty_file
        in_bed, in_sam, out = bed_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(empty_file),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_file, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_empty_sam_file(self, monkeypatch, capsys, empty_file, bed_sam_files):
        """Test main function with an empty sam file."""
        empty_file, empty_sam = empty_file
        in_bed, in_sam, out = bed_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(empty_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_file, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_complete_sam_input(self, monkeypatch, capsys, bed_sam_files):
        """Test main function with the correct input files."""
        in_bed, in_sam, expected_out = bed_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             '--feat-extension',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_missing_alignment(self, monkeypatch, capsys, bed_sam_missing_aln_files):
        """Test main function with missing alignments in sam file."""
        in_bed, in_sam, expected_out = bed_sam_missing_aln_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             '--feat-extension',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_no_extension(self, monkeypatch, capsys, bed_sam_no_extension_files):
        """Test main function with no extension in features names."""
        in_bed, in_sam, expected_out = bed_sam_no_extension_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_id_extension(self, monkeypatch, capsys, bed_extension_id_files):
        """Test main function with extension in feature name and read names."""
        in_bed, in_sam, expected_out = bed_extension_id_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             '--feat-extension',
             '--read-ids',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_id(self, monkeypatch, capsys, bed_id_files):
        """Test main function with read names."""
        in_bed, in_sam, expected_out = bed_id_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             '--read-ids',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_id_missing_alns(self, monkeypatch, capsys, bed_id_missing_aln_files):
        """Test main function with read names and missing alignments."""
        in_bed, in_sam, expected_out = bed_id_missing_aln_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--sam', str(in_sam),
             '--read-ids',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_no_sam_file(self, monkeypatch, capsys, bed_no_sam_files):
        """Test main function with read names."""
        in_bed, expected_out = bed_no_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_some_extension_file(self, monkeypatch, capsys, bed_some_extension_files):
        """Test main function with read names."""
        in_bed, expected_out = bed_some_extension_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--feat-extension',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()
