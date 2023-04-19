"""Unit tests for module 'filter_multimappers.py'."""

import argparse
import sys

from pathlib import Path
import pysam
import pytest

sys.path.append("../../")

from scripts.filter_multimappers import (
    count_indels,
    find_best_alignments,
    main,
    parse_arguments,
    write_output
)


@pytest.fixture
def sam_empty_files():
    """Import path to empty test files."""
    in_empty = Path("scripts/tests/files/in_sam_empty.sam")
    out_empty = Path("scripts/tests/files/out_sam_empty.sam")

    return in_empty, out_empty


@pytest.fixture
def sam_multimappers_files():
    """Import path to test files with multimappers."""
    in_multimappers = Path("scripts/tests/files/in_sam_multimappers.sam")
    out_multimappers = Path("scripts/tests/files/out_sam_multimappers.sam")

    return in_multimappers, out_multimappers


@pytest.fixture
def sam_no_multimappers_file():
    """Import path to test files with no multimappers."""
    no_multi = Path("scripts/tests/files/sam_no_multimappers.sam")

    return no_multi

@pytest.fixture
def sam_unique_diff_multimappers_files():
    """Import path to test files with a single multimapper."""
    in_diff_multi = Path("scripts/tests/files/in_sam_diff_multimappers.sam")
    out_diff_multi = Path("scripts/tests/files/out_sam_diff_multimappers.sam")

    return in_diff_multi, out_diff_multi

@pytest.fixture
def sam_unique_equal_multimapper_files():
    """Import path to the test file with a single multimapper."""
    in_sam = Path("scripts/tests/files/in_sam_equal_multimappers.sam")
    out_sam = Path("scripts/tests/files/out_sam_equal_multimappers.sam")

    return in_sam, out_sam


@pytest.fixture
def alns():
    """Create sample AlignedSegment objects."""
    # Alignment with 10 matches
    aln1 = pysam.AlignedSegment()
    aln1.cigartuples = [(0, 10)]
    aln1.query_name = "read1"
    aln1.query_sequence = "TAAAGCGCTT"
    aln1.flag = 0x2
    aln1.reference_id = 19
    aln1.reference_start = 63250
    aln1.set_tag("NH", 2)
    aln1.set_tag("HI", 1)

    # Alignment with 3 insertions
    aln2 = pysam.AlignedSegment()
    aln2.cigartuples = [(0, 5), (1, 3), (0, 5)]
    aln2.query_name = "read1"
    aln2.query_sequence = "GGGGCGTTTT"
    aln2.flag = 0x2
    aln2.reference_id = 19
    aln2.reference_start = 7589
    aln2.set_tag("NH", 2)
    aln2.set_tag("HI", 2)

    # Alignment with 3 deletions
    aln3 = pysam.AlignedSegment()
    aln3.cigartuples = [(0, 10), (2, 3), (0, 5)]
    aln3.query_name = "read2"
    aln3.query_sequence = "GCCAGGTGGCGTTTT"
    aln3.flag = 0x2
    aln3.reference_id = 19
    aln3.reference_start = 142777
    aln3.set_tag("NH", 3)
    aln3.set_tag("HI", 1)

    # Alignments with 1 insertion and 1 deletion
    aln4 = pysam.AlignedSegment()
    aln4.cigartuples = [(0, 10), (1, 1), (0, 3), (2, 1), (0, 2)]
    aln4.query_name = "read2"
    aln4.query_sequence = "AAGCCTCCCACCTAG"
    aln4.flag = 0x2
    aln4.reference_id = 19
    aln4.reference_start = 63251
    aln4.set_tag("NH", 3)
    aln4.set_tag("HI", 2)

    aln5 = pysam.AlignedSegment()
    aln5.cigartuples = [(0, 10), (2, 1), (0, 2), (1, 1), (0, 3)]
    aln5.query_name = "read2"
    aln5.query_sequence = "AGGTGGCGTTTTTCT"
    aln5.flag = 0x2
    aln5.reference_id = 19
    aln5.reference_start = 77595
    aln5.set_tag("NH", 3)
    aln5.set_tag("HI", 2)

    return [aln1, aln2, aln3, aln4, aln5]


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without input file."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['filter_multimappers']
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, sam_no_multimappers_file):
        """Call with a single input file."""
        sam_1 = sam_no_multimappers_file
        monkeypatch.setattr(
            sys, 'argv',
            ['filter_multimappers',
             str(sam_1),
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_too_many_inputs(self, monkeypatch, sam_multimappers_files):
        """Call with too many input files."""
        sam_1, sam_2 = sam_multimappers_files
        monkeypatch.setattr(
            sys, 'argv',
            ['filter_multimappers',
             str(sam_1), str(sam_2),
             ]
        )
        with pytest.raises(SystemExit) as sysex:
            parse_arguments().parse_args()
        assert sysex.value.code == 2


class TestCountIndels:
    """Test 'count_indels()' function."""

    def test_no_indels(self, alns):
        """Test CIGAR string with 10 matches."""
        assert count_indels(alns[0]) == 0

    def test_three_insertions(self, alns):
        """Test CIGAR string with 3 insertions."""
        assert count_indels(alns[1]) == 3

    def test_three_deletions(self, alns):
        """Test CIGAR string with 3 deletions."""
        assert count_indels(alns[2]) == 3

    def test_one_in_one_del(self, alns):
        """Test CIGAR string with one insertion and one deletion."""
        assert count_indels(alns[3]) == 2


class TestFindBestAlignments:
    """Test 'find_best_alignments()' function."""

    def test_find_best_alignments_multimappers(self, alns):
        """Test function with multimappers with different indel count."""
        output = find_best_alignments([alns[0], alns[1]])

        assert len(output) == 1
        assert output[0] == alns[0]

        assert output[0].get_tag("NH") == 1
        assert output[0].get_tag("HI") == 1

    def test_find_best_alignments_equal_multimappers(self, alns):
        """Test function with multimappers with same indel count."""
        output = find_best_alignments([alns[3], alns[4]])

        assert len(output) == 2
        assert output[0] == alns[3]
        assert output[1] == alns[4]

        assert output[0].get_tag("NH") == 2
        assert output[1].get_tag("NH") == 2
        assert output[0].get_tag("HI") == 1
        assert output[1].get_tag("HI") == 2


class TestWriteOutout:
    """Test 'write_output()' function."""

    def test_write_output_one_alignment(self, capsys, sam_multimappers_files):
        """Test funciton with a single alignment."""
        in_sam, out_sam = sam_multimappers_files

        with pysam.AlignmentFile(in_sam, 'r') as in_file:
            alignment = next(in_file)
        
        write_output([alignment])
        captured = capsys.readouterr()

        with pysam.AlignmentFile(out_sam, 'r') as out_file:
            out_alignment = next(out_file)
        
        assert captured.out == out_alignment.to_string() + '\n'

    def test_write_output_multiple_alignments_diff_indels(self, capsys, sam_unique_diff_multimappers_files):
        """Test function with multimappers with different amount of indels."""
        in_sam, out_sam = sam_unique_diff_multimappers_files

        with pysam.AlignmentFile(in_sam, 'r') as in_file:
            alignments_in = [aln for aln in in_file]

        write_output(alignments_in)
        captured = capsys.readouterr()

        with open(out_sam, 'r') as out_file:
            expected_output = out_file.read()

        assert captured.out == expected_output

    def test_write_output_multiple_alignments_equal_indels(self, capsys, sam_unique_equal_multimapper_files):
        """Test function with equal multimappers."""
        in_sam, out_sam = sam_unique_equal_multimapper_files

        with pysam.AlignmentFile(in_sam, 'r') as in_file:
            alignments_in = [aln for aln in in_file]
        
        write_output(alignments_in)
        captured = capsys.readouterr()

        with open(out_sam, 'r') as out_file:
            expected_output = out_file.read()

        assert captured.out == expected_output


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_file(self, capsys, sam_empty_files):
        """Test main function with an empty file."""
        in_sam, out_sam = sam_empty_files

        main(in_sam)
        captured = capsys.readouterr()

        with open(out_sam, 'r') as out_file:
            expected_output = out_file.read()

        assert captured.out == expected_output

    def test_main_multimappers(self, capsys, sam_multimappers_files):
        """Test main function with multimappers."""
        in_sam, out_sam = sam_multimappers_files

        main(in_sam)
        captured = capsys.readouterr()

        with open(out_sam, 'r') as out_file:
            expected_output = out_file.read()

        assert captured.out == expected_output

    def test_main_no_multimappers(self, capsys, sam_no_multimappers_file):
        """Test main function with no multimappers."""
        sam_file = sam_no_multimappers_file

        main(sam_file)
        captured = capsys.readouterr()

        with open(sam_file, 'r') as out_file:
            expected_output = out_file.read()

        assert captured.out == expected_output
