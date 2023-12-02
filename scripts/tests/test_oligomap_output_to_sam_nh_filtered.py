"""Unit tests for module 'oligomap_output_to_sam_nh_filtered.py'."""

import argparse
from pathlib import Path
import sys

import pytest

from ..oligomap_output_to_sam_nh_filtered import (
    eval_aln,
    Fields,
    get_cigar_md,
    get_sam_fields,
    main,
    parse_arguments
)


@pytest.fixture
def empty_file():
    """Import path to empty file."""
    empty_file = Path("files/empty_file")

    return empty_file


@pytest.fixture
def genome_nh_2():
    """Import path to test files with maximum NH set to 2."""
    oligo_in = Path("files/in_oligomap_output.oligomap")
    oligo_out = Path("files/oligomap_genome_2_nh.sam")

    return oligo_in, oligo_out


@pytest.fixture
def transcriptome_no_nh():
    """Import path to test files with reference set to transcriptome."""
    oligo_in = Path("files/in_oligomap_output.oligomap")
    oligo_out = Path("files/oligomap_transcriptome_no_nh.sam")

    return oligo_in, oligo_out


@pytest.fixture
def single_read():
    """Import path to test file with a single read."""
    oligo_out = Path("files/oligomap_single_read.sam")

    return oligo_out


@pytest.fixture
def aln_fields():
    """Create sample alignment as a Fields class NamedTuple."""
    # Perfect alignment
    field_1 = Fields("read_1", "0", "19", "44377", "255", "19M", '*', '0', '0',
                     "CTACAAAGGGAAGCACTTT", '*', "NM:i:0", "MD:Z:19")

    # Alignment with a mismatch in the first position
    field_2 = Fields("read_1", '0', "19", "53471", "255", "19M", '*', '0', '0',
                     "CTACAAAGGGAAGCACTTT", '*', "NM:i:1", "MD:Z:G18")

    # Alignment with a mismatch in the last position
    field_3 = Fields("read_1", "0", "19", "44278", "255", "19M", '*', '0', '0',
                     "CTACAAAGGGAAGCACTTT", '*', "NM:i:1", "MD:Z:18C")

    # Alignment with a mismatch in the middle of the read sequence
    field_4 = Fields("read_1", "0", "19", "50971", "255", "19M", '*', '0', '0',
                     "CTACAAAGGGAAGCACTTT", '*', "NM:i:1", "MD:Z:14C4")

    # Alignment with an insertion at read's first position
    field_5 = Fields("read_2", "16", "19", "7627", "255", "1I22M", '*', '0',
                     '0', "AAAGCACCTCCAGAGCTTGAAGC", '*', "NM:i:1", "MD:Z:23")

    # Alignment with an insertion in the middle of the read sequence
    field_6 = Fields("read_2", "16", "19", "7886", "255", "9M1I12M", '*', '0',
                     '0', "AAAGCACCTCCAGAGCTTGAAGC", '*', "NM:i:1", "MD:Z:23")

    return [field_1, field_2, field_3, field_4, field_5, field_6]


@pytest.fixture
def alns():
    """Create sample alignment lists to extract the CIGAR and MD strings.

    The list contain the number of errors, the read's sequence, the alignment
    representation in bars and the reference sequence as strings in that
    order.
    """
    # Perfect alignment
    aln1 = []
    aln1.append("0")
    aln1.append("CTACAAAGGGAAGCACTTT")
    aln1.append("|||||||||||||||||||")
    aln1.append("CTACAAAGGGAAGCACTTT")

    # Alignment with a mismatch in the first position
    aln2 = []
    aln2.append("1")
    aln2.append("CTACAAAGGGAAGCACTTT")
    aln2.append(" ||||||||||||||||||")
    aln2.append("GTACAAAGGGAAGCACTTT")

    # Alignment with a mismatch in the last position
    aln3 = []
    aln3.append("1")
    aln3.append("CTACAAAGGGAAGCACTTT")
    aln3.append("|||||||||||||||||| ")
    aln3.append("CTACAAAGGGAAGCACTTC")

    # Alignment with a mismatch in the middle of the sequence
    aln4 = []
    aln4.append("1")
    aln4.append("CTACAAAGGGAAGCACTTT")
    aln4.append("|||||||||||||| ||||")
    aln4.append("CTACAAAGGGAAGCCCTTT")

    # Alignment with an insertion at read's first position
    aln5 = []
    aln5.append("1")
    aln5.append("AAAGCACCTCCAGAGCTTGAAGC")
    aln5.append(" ||||||||||||||||||||||")
    aln5.append("-AAGCACCTCCAGAGCTTGAAGC")

    # Alignment with an insertion at read's last position
    aln6 = []
    aln6.append("1")
    aln6.append("AAAGCACCTCCAGAGCTTGAAGC")
    aln6.append("|||||||||||||||||||||| ")
    aln6.append("AAAGCACCTCCAGAGCTTGAAG-")

    # Alignment with an insertion in the middle of the read
    aln7 = []
    aln7.append("1")
    aln7.append("AAAGCACCTCCAGAGCTTGAAGC")
    aln7.append("||||||||| |||||||||||||")
    aln7.append("AAAGCACCT-CAGAGCTTGAAGC")

    # Alignment with a deletion at read's first position
    aln8 = []
    aln8.append("1")
    aln8.append("-GAAGGCGCTTCACCTTTGGAGT")
    aln8.append(" ||||||||||||||||||||||")
    aln8.append("TGAAGGCGCTTCACCTTTGGAGT")

    # Alignment with a deletion at read's last position
    aln9 = []
    aln9.append("1")
    aln9.append("GAAGGCGCTTCACCTTTGGAGT-")
    aln9.append("|||||||||||||||||||||| ")
    aln9.append("GAAGGCGCTTCACCTTTGGAGTA")

    # Alignment with a deletion in the middle of the read
    aln10 = []
    aln10.append("1")
    aln10.append("GAAGGCGCTTC-CCTTTGGAGT")
    aln10.append("||||||||||| ||||||||||")
    aln10.append("GAAGGCGCTTCACCTTTGGAGT")

    return [aln1, aln2, aln3, aln4, aln5, aln6, aln7, aln8, aln9, aln10]


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_files(self, monkeypatch):
        """Call without input nor output files."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['oligomap_output_to_sam_nh_filtered']
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_in_files(self, monkeypatch, empty_file):
        """Call with in file."""
        empty_in = empty_file

        monkeypatch.setattr(
            sys, 'argv',
            ['oligomap_output_to_sam_nh_filtered',
             str(empty_in),
             ]
        )

        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_all_arguments(self, monkeypatch, genome_nh_2):
        """Call with all the arguments."""
        fa_in, sam_out = genome_nh_2

        monkeypatch.setattr(
            sys, 'argv',
            ['oligomap_output_to_sam_nh_filtered',
             str(fa_in),
             '-n', '100',
             ]
        )

        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestGetCigarMd:
    """Test 'get_cigar_md()' function."""

    def test_perfect_aln(self, alns):
        """Test perfect alignment."""
        result = ("19M", "MD:Z:19")

        assert get_cigar_md(alns[0][0], alns[0][1],
                            alns[0][2], alns[0][3]) == result

    def test_mm_first_pos_aln(self, alns):
        """Test mismatch at read's first position."""
        result = ("19M", "MD:Z:G18")

        assert get_cigar_md(alns[1][0], alns[1][1],
                            alns[1][2], alns[1][3]) == result

    def test_mm_last_pos_aln(self, alns):
        """Test mismatch at read's last position."""
        result = ("19M", "MD:Z:18C")

        assert get_cigar_md(alns[2][0], alns[2][1],
                            alns[2][2], alns[2][3]) == result

    def test_mm_middle_aln(self, alns):
        """Test mismatch in the middle of the read."""
        result = ("19M", "MD:Z:14C4")

        assert get_cigar_md(alns[3][0], alns[3][1],
                            alns[3][2], alns[3][3]) == result

    def test_in_first_pos_aln(self, alns):
        """Test insertion at read's first position."""
        result = ("1I22M", "MD:Z:23")

        assert get_cigar_md(alns[4][0], alns[4][1],
                            alns[4][2], alns[4][3]) == result

    def test_in_last_pos_aln(self, alns):
        """Test insertion at read's last position."""
        result = ("22M1I", "MD:Z:23")

        assert get_cigar_md(alns[5][0], alns[5][1],
                            alns[5][2], alns[5][3]) == result

    def test_in_middle_aln(self, alns):
        """Test insertion in the middle of the read."""
        result = ("9M1I13M", "MD:Z:23")

        assert get_cigar_md(alns[6][0], alns[6][1],
                            alns[6][2], alns[6][3]) == result

    def test_del_first_pos_aln(self, alns):
        """Test deletion at read's first position."""
        result = ("1D22M", "MD:Z:^T22")

        assert get_cigar_md(alns[7][0], alns[7][1],
                            alns[7][2], alns[7][3]) == result

    def test_del_last_pos_aln(self, alns):
        """Test deletion at read's last position."""
        result = ("22M1D", "MD:Z:22^A0")

        assert get_cigar_md(alns[8][0], alns[8][1],
                            alns[8][2], alns[8][3]) == result

    def test_del_middle_aln(self, alns):
        """Test deletion in the middle of the read."""
        result = ("11M1D10M", "MD:Z:11^A10")

        assert get_cigar_md(alns[9][0], alns[9][1],
                            alns[9][2], alns[9][3]) == result


class TestGetSAMFields():
    """Test 'get_sam_fields()' function."""

    def test_pos_strand_no_err(self, alns, aln_fields):
        """Test perfect alignment in the positive strand."""
        line1 = "read_1 (19 nc) 1...19 19 44377...44395"
        line2 = "19"
        line3 = "errors: 0 orientation: +"

        assert get_sam_fields([line1, line2, line3, alns[0][1],
                              alns[0][2], alns[0][3]]) == aln_fields[0]

    def test_neg_strand_one_err(self, alns, aln_fields):
        """Test alignment with an insertion in the negative strand."""
        line1 = "read_2 (23 nc) 1...23 19 7886...7908"
        line2 = "19"
        line3 = "errors: 1 orientation: -"

        assert get_sam_fields([line1, line2, line3, alns[6][1],
                              alns[6][2], alns[6][3]]) == aln_fields[5]


class TestEvalAln:
    """Test ''eval_aln()' function."""

    def test_eval_empty_dict_new_read(self, aln_fields):
        """Test evaluation with a new read and an empty dictionary."""
        d = dict()
        minerr_nh = {"read_0": ['0', 1]}
        aln = aln_fields[0]
        nhfilter = None

        eval_aln(nhfilter, d, minerr_nh, aln)

        assert list(d.keys())[0] == aln.read_name
        assert minerr_nh[aln.read_name] == ['0', 1]

    def test_eval_empty_dict_smaller_error(self, aln_fields):
        """Test evaluation with a smaller error and an empty dictionary."""
        d = dict()
        minerr_nh = {"read_1": ['1', 1]}
        aln = aln_fields[0]
        nhfilter = None

        eval_aln(nhfilter, d, minerr_nh, aln)

        assert list(d.keys())[0] == aln.read_name
        assert minerr_nh[aln.read_name] == ['0', 1]

    def test_increase_nh_no_filter(self, aln_fields):
        """Test evaluation when increasing NH without a maximum value."""
        d = {"read_1": [aln_fields[1], aln_fields[2]]}
        minerr_nh = {"read_1": ['1', 2]}
        aln = aln_fields[3]
        nhfilter = None

        eval_aln(nhfilter, d, minerr_nh, aln)

        assert len(d[aln.read_name]) == 3
        assert minerr_nh[aln.read_name] == ['1', 3]

    def test_exceed_nh_filter_2(self, capsys, aln_fields):
        """Test evaluation when exceeding the maximum NH set to 2."""
        d = {"read_1": [aln_fields[1], aln_fields[2]]}
        minerr_nh = {"read_1": ['1', 2]}
        aln = aln_fields[3]
        nhfilter = 2

        eval_aln(nhfilter, d, minerr_nh, aln)
        captured = capsys.readouterr()

        assert len(d) == 0
        assert minerr_nh[aln.read_name] == ['1', 3]
        assert captured.err == "Filtered by NH | Read read_1 | Errors = 1\n"

    def test_no_exceed_nh_filter_2(self, aln_fields):
        """Test evaluation when increasing NH with maximum value of 2."""
        d = {"read_1": [aln_fields[1]]}
        minerr_nh = {"read_1": ['1', 1]}
        aln = aln_fields[2]
        nhfilter = 2

        eval_aln(nhfilter, d, minerr_nh, aln)

        assert len(d[aln.read_name]) == 2
        assert minerr_nh[aln.read_name] == ['1', 2]

    def test_smaller_min_error(self, capsys, aln_fields):
        """Test evaluation when having a smaller minimumm error."""
        d = {"read_1": [aln_fields[1], aln_fields[2]]}
        minerr_nh = {"read_1": ['1', 2]}
        aln = aln_fields[0]
        nhfilter = None

        eval_aln(nhfilter, d, minerr_nh, aln)
        captured = capsys.readouterr()

        assert len(d[aln.read_name]) == 1
        assert minerr_nh[aln.read_name] == ['0', 1]
        assert captured.err == "Filtered by ERROR | Read read_1 | Errors = 1\n"

    def test_different_read(self, capsys, tmp_path, aln_fields, single_read):
        """Test evaluation when having to write due to a different read."""
        out_file = single_read

        d = {"read_1": [aln_fields[1], aln_fields[2]]}
        minerr_nh = {"read_1": ['1', 2]}
        aln = aln_fields[4]
        nhfilter = None

        eval_aln(nhfilter, d, minerr_nh, aln)
        captured = capsys.readouterr()

        assert list(d.keys())[0] == aln.read_name
        assert len(minerr_nh) == 1
        assert captured.err == "Written read read_1 | Errors = 1 | NH = 2\n"

        with open(out_file, 'r') as expected:
            assert captured.out == expected.read()


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_file(self, monkeypatch, capsys, empty_file):
        """Test main function with an empty file."""
        empty_in = empty_file

        monkeypatch.setattr(
            sys, 'argv',
            ['oligomap_output_to_sam_nh_filtered',
             str(empty_in),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_in, 'r') as expected:
            assert captured.out == expected.read()

    def test_main_max_nh_2(self, monkeypatch, capsys, genome_nh_2):
        """Test main function with NH set to 2."""
        in_file, out_file = genome_nh_2

        monkeypatch.setattr(
            sys, 'argv',
            ['oligomap_output_to_sam_nh_filtered',
             str(in_file),
             '-n', '2',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(out_file, 'r') as expected:
            assert captured.out == expected.read()

    def test_main_no_nh_transcriptome(self, monkeypatch, capsys,
                                      transcriptome_no_nh):
        """Test main function with no NH set for transcriptome mappings."""
        in_file, out_file = transcriptome_no_nh

        monkeypatch.setattr(
            sys, 'argv',
            ['oligomap_output_to_sam_nh_filtered',
             str(in_file),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(out_file, 'r') as expected:
            assert captured.out == expected.read()
