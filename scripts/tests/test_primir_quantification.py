"""Unit tests for module 'primir_quantification.py'"""

import argparse
from pathlib import Path
import sys

import pysam
import pytest

sys.path.append("../../")

from scripts.primir_quantification import (
    main,
    parse_arguments
)

@pytest.fixture
def bed_sam_empty_files():
    """Import path to empty bed test file."""
    in_bed_empty = Path("files/in_bed_empty.bed")
    in_sam_empty = Path("files/in_sam_empty.sam")
    out_empty_table = Path("files/out_empty_primir_quantification")

    return in_bed_empty, in_sam_empty, out_empty_table

@pytest.fixture
def bed_sam_files():
    """Import path to test files with correct content."""
    in_bed = Path("files/in_intersection.bed")
    in_sam = Path("files/in_alignments.sam")
    out_table = Path("files/out_primir_quantification")

    return in_bed, in_sam, out_table

@pytest.fixture
def bed_sam_missing_aln_files():
    """Import path to test files with missing alignments."""
    in_bed = Path("files/in_intersection.bed")
    in_sam = Path("files/in_alignments_missing.sam")
    out_table = Path("files/out_missing_alns_primir_quantification")

    return in_bed, in_sam, out_table

class TestParseArguments:
    """"Test 'parse_arguments()' function."""

    def test_no_sam_input(self, monkeypatch, bed_sam_files):
        """Call without sam input file."""
        in_bed, in_sam, out_table = bed_sam_files
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                '--bed', str(in_bed),
                '--outfile', str(out_table),
                ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2
    
    def test_no_bed_input(self, monkeypatch, bed_sam_files):
        """Call without bed input file."""
        in_bed, in_sam, out_table = bed_sam_files

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                '--sam', str(in_sam),
                '--outfile', str(out_table),
                ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_no_output(self, monkeypatch, bed_sam_files):
        """Call without output file."""
        in_bed, in_sam, out_table = bed_sam_files

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                '--bed', str(in_bed),
                '--sam', str(in_sam),
                ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2
    
    def test_correct_input(self, monkeypatch, bed_sam_files):
        """Call with the correct input files."""
        in_bed, in_sam, out_table = bed_sam_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
            '--bed', str(in_bed),
            '--sam', str(in_sam),
            '--outfile', str(out_table),
            ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

        
class TestMain:
    """Test 'main()' function."""

    def test_main_empty_bed_file(self, monkeypatch, tmp_path, bed_sam_empty_files, bed_sam_files):
        """Test main function with an empty bed file."""
        in_empty_bed, in_empty_sam, expected_out = bed_sam_empty_files
        in_bed, in_sam, out = bed_sam_files

        out_table = tmp_path/"primir_quantification"

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
            '--bed', str(in_empty_bed),
            '--sam', str(in_sam),
            '--outfile', str(out_table),
            ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(expected_out, 'r') as expected, open(out_table, 'r') as output:
            assert output.read() == expected.read()

    
    def test_main_empty_sam_file(self, monkeypatch, tmp_path, bed_sam_empty_files, bed_sam_files):
        """Test main function with an empty sam file."""
        in_empty_bed, in_empty_sam, expected_out = bed_sam_empty_files
        in_bed, in_sam, out = bed_sam_files

        out_table = tmp_path/"primir_quantification"

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
            '--bed', str(in_bed),
            '--sam', str(in_empty_sam),
            '--outfile', str(out_table),
            ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(expected_out, 'r') as expected, open(out_table, 'r') as output:
            assert output.read() == expected.read()

    def test_main_correct_input(self, monkeypatch, tmp_path, bed_sam_files):
        """Test main function with the correct input files."""
        in_bed, in_sam, expected_out = bed_sam_files

        out_table = tmp_path/"primir_quantification"

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
            '--bed', str(in_bed),
            '--sam', str(in_sam),
            '--outfile', str(out_table),
            ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(expected_out, 'r') as expected, open(out_table, 'r') as output:
            assert output.read() == expected.read()

    def test_main_missing_alignment(self, monkeypatch, tmp_path, bed_sam_missing_aln_files):
        """Test main function with missing alignments in sam file."""
        in_bed, in_sam, expected_out = bed_sam_missing_aln_files

        out_table = tmp_path/"primir_quantification"

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
            '--bed', str(in_bed),
            '--sam', str(in_sam),
            '--outfile', str(out_table),
            ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(expected_out, 'r') as expected, open(out_table, 'r') as output:
            assert output.read() == expected.read()
        