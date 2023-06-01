"""Unit tests for module 'mirna_quantification.py'."""
import argparse
from pathlib import Path
import sys

import pysam
import pytest

sys.path.append("../../")

from scripts.mirna_quantification import (
    collapsed_nh_contribution,
    collapsed_contribution,
    nh_contribution,
    contribution,
    get_name,
    main,
    parse_arguments
)


@pytest.fixture
def empty_file():
    """Import path to empty file."""
    empty_in = Path("files/header_only.sam")
    empty_out = Path("files/empty_file")

    return empty_in, empty_out


@pytest.fixture
def sam_file():
    """Import path to test files with full content."""
    sam_file = Path("files/in_aln_tag.sam")
    out_table = Path("files/iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def iso_mirna_sam_file():
    """Import path to test files with full content."""
    sam_file = Path("files/in_aln_tag.sam")
    iso_out_table = Path("files/iso_quantification")
    mirna_out_table = Path("files/mirna_quantification")

    return sam_file, iso_out_table, mirna_out_table


@pytest.fixture
def xn_tag_sam_file():
    """Import path to test files with feat names in XN tag."""
    sam_file = Path("files/in_aln_tag_XN.sam")
    out_table = Path("files/xn_tag_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def nh_missing_sam_file():
    """Import path to test files with missing NH tag in some alignments."""
    sam_file = Path("files/in_aln_tag_missing_nh.sam")
    out_table = Path("files/missing_nh_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def len_sam_file():
    """Import path to test files with feature length in the output table."""
    sam_file = Path("files/in_aln_tag.sam")
    out_table = Path("files/len_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def read_sam_file():
    """Import path to test files with read IDs in the output table."""
    sam_file = Path("files/in_aln_tag.sam")
    out_table = Path("files/read_ids_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def read_len_sam_file():
    """Import path to test files with read IDs and feature len in the output table."""
    sam_file = Path("files/in_aln_tag.sam")
    out_table = Path("files/len_ids_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def uncollapsed_sam_file():
    """Import path to uncollapsed test files."""
    sam_file = Path("files/in_aln_tag.sam")
    out_table = Path("files/uncollpased_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def uncollapsed_missing_nh_sam_file():
    """Import path to uncollapsed test files."""
    sam_file = Path("files/in_aln_tag_missing_nh.sam")
    out_table = Path("files/uncollpased_missing_nh_iso_mirna_quantification")

    return sam_file, out_table


@pytest.fixture
def alns():
    """Create sample AlignedSegment objects."""
    # Collapsed lignment with NH tag in the name
    aln1 = pysam.AlignedSegment()
    aln1.query_name = "read1-2_3"
    aln1.query_sequence = "TAAAGCGCTT"
    aln1.reference_id = 19
    aln1.reference_start = 63250

    # Uncollapsed alignment with NH tag in the name
    aln2 = pysam.AlignedSegment()
    aln2.query_name = "read2_4"
    aln2.query_sequence = "GGGGCGTTTT"
    aln2.reference_id = 19
    aln2.reference_start = 7589

    # Collpased alignment without NH tag in the name
    aln3 = pysam.AlignedSegment()
    aln3.query_name = "read3-6"
    aln3.query_sequence = "GCCAGGTGGCGTTTT"
    aln3.reference_id = 19
    aln3.reference_start = 142777
    aln3.set_tag("NH", 3)

    # Uncollapsed alignment without NH tag in the name
    aln4 = pysam.AlignedSegment()
    aln4.query_name = "read4"
    aln4.query_sequence = "AAGCCTCCCACCTAG"
    aln4.reference_id = 19
    aln4.reference_start = 63251
    aln4.set_tag("NH", 8)

    # Uncollapsed alignment missing the NH tag
    aln5 = pysam.AlignedSegment()
    aln5.query_name = "read5"
    aln5.query_sequence = "AGGTGGCGTTTTTCT"
    aln5.reference_id = 19
    aln5.reference_start = 77595

    return [aln1, aln2, aln3, aln4, aln5]


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without input file."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['mirna_quantification',
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, empty_file):
        """Call with the correct input file."""
        in_sam = empty_file

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna__quantification',
             str(in_sam),
             '--lib', 'test_lib',
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_too_many_input_files(self, monkeypatch, empty_file):
        """Call with too many input file."""
        in_sam = empty_file

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['mirna_quantification',
                 str(in_sam), str(in_sam),
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_all_input(self, monkeypatch, empty_file):
        """Call with all the options."""
        in_sam = empty_file

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(in_sam),
             '--len',
             '--read-ids',
             '--collapsed',
             '--nh',
             '--tag', 'YW',
             '--outdir', 'Path.cwd()',
             '--lib', 'test_lib',
             '--mir-list', '[mirna, isomir]'
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestGetContribution():
    """Test 'get_contribution()' function."""

    def test_collapsed_nh(self, alns):
        """Test collapsed alignment with NH in the name."""
        assert collapsed_nh_contribution(alns[0]) == 2/3

    def test_uncollpased_nh(self, alns):
        """Test uncollapsed alignment with NH in the name."""
        assert nh_contribution(alns[1]) == 1/4

    def test_collapsed_no_nh(self, alns):
        """Test collapsed alignment without NH in the name."""
        assert collapsed_contribution(alns[2]) == 2

    def test_uncollpased_no_nh(self, alns):
        """Test uncollapsed alignment without NH in the name."""
        assert contribution(alns[3]) == 1/8

    def test_uncollpased_missing_nh(self, alns):
        """Test uncollapsed alignment with missing NH value."""
        assert contribution(alns[4]) == 1


class TestGetName:
    """Test 'get_name()' function."""

    def test_canonical(self):
        """Test canonical miRNA name."""
        name = ["mirna", "hsa-miR-1323"]
        assert get_name("hsa-miR-1323|0|0|22M|22") == name

    def test_iso_0_shift(self):
        """Test isoform with 0 shift."""
        name = ["isomir", "hsa-miR-1323|0|0|18M3I4M|22"]
        assert get_name("hsa-miR-1323|0|0|18M3I4M|22") == name

    def test_iso(self):
        """Test isoform with shift."""
        name = ["isomir", "hsa-miR-1323|2|0|18M3I4M|22"]
        assert get_name("hsa-miR-1323|2|0|18M3I4M|22") == name


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_sam_file(self, monkeypatch, tmp_path, empty_file):
        """Test main function with an empty SAM file."""
        empty_in, empty_out = empty_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(empty_in),
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(empty_out, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_isomir_mirna_sam_file(self, monkeypatch, tmp_path, sam_file):
        """Test main function with complete SAM file."""
        infile, out_table = sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_iso_sam_file(self, monkeypatch, tmp_path, iso_mirna_sam_file):
        """Test main function tabulating only isomiRs."""
        infile, iso_out_table, mirna_out_table = iso_mirna_sam_file
        mirna_output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--outdir', str(tmp_path),
             '--mir-list', "isomir",
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(iso_out_table, 'r') as expected, open(mirna_output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_mirna_sam_file(self, monkeypatch, tmp_path, iso_mirna_sam_file):
        """Test main function tabulating only canonical miRNA."""
        infile, iso_out_table, mirna_out_table = iso_mirna_sam_file
        mirna_output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--outdir', str(tmp_path),
             '--mir-list',  "mirna"
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(mirna_out_table, 'r') as expected, open(mirna_output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_xn_tag_sam_file(self, monkeypatch, tmp_path, xn_tag_sam_file):
        """Test main function with feature name in the XN tag."""
        infile, out_table = xn_tag_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--tag', 'XN',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_nh_missing_sam_file(self, monkeypatch, tmp_path, nh_missing_sam_file):
        """Test main function with some missing NH tag in SAM file."""
        infile, out_table = nh_missing_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_seq_len_sam_file(self, monkeypatch, tmp_path, len_sam_file):
        """Test main function with read lenght in output table."""
        infile, out_table = len_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--len',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_read_sam_file(self, monkeypatch, tmp_path, read_sam_file):
        """Test main function with intersecting read IDs in the output."""
        infile, out_table = read_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--read-ids',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_read_len_sam_file(self, monkeypatch, tmp_path, read_len_sam_file):
        """Test main function with read IDs and feature length im output."""
        infile, out_table = read_len_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--collapsed',
             '--nh',
             '--len',
             '--read-ids',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()


    def test_main_uncollpased_sam_file(self, monkeypatch, tmp_path, uncollapsed_sam_file):
        """Test main function with uncollapsed SAM file."""
        infile, out_table = uncollapsed_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--nh',
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()

    def test_main_uncollpased_missing_nh_sam_file(self, monkeypatch, tmp_path, uncollapsed_missing_nh_sam_file):
        """Test main function with uncollapsed SAM file and missing NH tags."""
        infile, out_table = uncollapsed_missing_nh_sam_file
        output = tmp_path/"mirna_counts_lib"

        monkeypatch.setattr(
            sys, 'argv',
            ['mirna_quantification',
             str(infile),
             '--outdir', str(tmp_path),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(out_table, 'r') as expected, open(output, 'r') as out_file:
            assert out_file.read() == expected.read()
