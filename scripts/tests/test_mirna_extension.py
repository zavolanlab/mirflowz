"""Unit tests for module 'mirna_extension.py'."""

import argparse
import sys

from pathlib import Path
import gffutils
import pytest

from ..mirna_extension import main, MirnaExtension, parse_arguments


@pytest.fixture
def gff_empty():
    """Import path to empty test file."""
    empty = Path("files/empty_file")

    return empty


@pytest.fixture
def gff_no_extremes():
    """Import path to miRNA annotation files."""
    in_no_extreme = Path("files/in_mirna_anno.gff3")
    out_primir = Path("files/primir_anno.gff3")
    out_mir = Path("files/mir_anno.gff3")

    return in_no_extreme, out_primir, out_mir


@pytest.fixture
def gff_extremes():
    """Import path to miRNA annotation files with extreme miRNA coords."""
    in_extremes = Path("files/in_mirna_extreme_mirs.gff3")
    out_primir = Path("files/extreme_primir_anno.gff3")
    out_mir = Path("files/extreme_mir_anno.gff3")

    return in_extremes, out_primir, out_mir


@pytest.fixture
def gff_extremes_chr():
    """Import path to miRNA annotation files."""
    in_chr_extremes = Path("files/in_mirna_extreme_chr_mirs.gff3")
    out_primir = Path("files/extreme_chr_primir_anno.gff3")
    out_mir = Path("files/extreme_chr_mir_anno.gff3")

    return in_chr_extremes, out_primir, out_mir


@pytest.fixture
def seq_len_tbl():
    """Import path to sequence lengths tables."""
    correct_seq_len = Path("files/seq_len.tsv")
    incorrect_seq_len = Path("files/seq_len_wrong.tsv")

    return correct_seq_len, incorrect_seq_len


class TestSetDb:
    """Test for the 'set_db' method."""

    def test_set_db_file(self, gff_no_extremes):
        """Test setting local db from file."""
        in_file, pre_exp, mir_exp = gff_no_extremes

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        assert miR_obj is not None
        assert isinstance(miR_obj.db, gffutils.FeatureDB)
        assert (
            len(
                list(
                    miR_obj.db.features_of_type("miRNA_primary_transcript")
                )
            )
            == 2
        )
        assert len(list(miR_obj.db.features_of_type("miRNA"))) == 4

    def test_set_db_stdin(self, monkeypatch, gff_no_extremes):
        """Test setting local db from standard input."""
        in_file, pre_exp, mir_exp = gff_no_extremes
        monkeypatch.setattr(sys, "stdin", str(in_file))

        miR_obj = MirnaExtension()
        miR_obj.set_db()

        assert miR_obj is not None
        assert isinstance(miR_obj.db, gffutils.FeatureDB)
        assert (
            len(
                list(
                    miR_obj.db.features_of_type("miRNA_primary_transcript")
                )
            )
            == 2
        )
        assert len(list(miR_obj.db.features_of_type("miRNA"))) == 4

    def test_set_db_empty_file(self, monkeypatch, gff_empty):
        """Test setting local db from empty file."""
        in_file = gff_empty

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        assert miR_obj.db is None
        assert isinstance(miR_obj.db, gffutils.FeatureDB)


class TestSetSeqLengths:
    """Test for the 'set_seq_lengths' method."""

    def test_set_lengths_no_tbl(self, ):
        """Test create sequence lengths dictionary from GFF3 file."""
        in_file, pre_exp, mir_exp = gff_no_extremes

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)
        miR_obj.set_seq_lengths()

        # Assert correct lengths?
        # Assert not empty?

    def test_set_lengths_wrong_tbl(self, ):
        """Test create sequence lengths dictionary from wrong table."""
        # Add tables

        miR_obj = MirnaExtension()
        miR_obj.set_seq_lengths()

        # Assert valueError?
        # Assert not empty?

    def test_set_lengths_correct_tbl(self, ):
        """Test create sequence lengths dictionary from correct table."""
        # Add tables

        miR_obj = MirnaExtension()
        miR_obj.set_seq_lengths()

        # Assert correct lengths?
        # Assert not empty?


class TestProcessPrecursor:
    """Test for the 'process_precursor' method."""

    def test_process_prec_diff_strand_same_coords(self, ):
        """Test processing precursor on different strands, similar coords."""
        # Add files and tables

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        precursor = miR_obj.process_precursor(precursor="name")

        # Assert correct precursor new name
        # Assert correct mature extension ?

    def test_process_prec_miR_out_seq_init_boundaries(self, ):
        """Test processing precursor for outside seq boundaries miRNAs."""
        # Add files and tables

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        precursor = miR_obj.process_precursor(precursor="name")

        # Assert ValueError

    def test_process_prec_unknown_seq(self, ):
        """Test processing precursor with annotated seq ID not in len dict."""
        # Add files and tables

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        precursor = miR_obj.process_precursor(precursor="name")

        # Assert KeyError

    def test_process_prec_no_precursor_extension(self, ):
        """Test processing precursor without precursor extension."""
        # Add files and tables

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        precursor = miR_obj.process_precursor(precursor="name")

        # Assert precursor name
        # Assert correct mature miR extension

    def test_process_prec_extended_miR_out_seq_boundaries(self, ):
        """Test processing precursor with extended miRs out seq boundaries."""
        # Add files

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        precursor = miR_obj.process_precursor(precursor="name")

        # Assert one miR starts at 0
        # Assert one miR ends and chr boundaries

    def test_process_prec_extended_miR_in_seq_boundaries(self, ):
        """Test processing precursor with extended miRs in seq boundaries."""
        # Add files and tables

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        precursor = miR_obj.process_precursor(precursor="name")

        # Assert no miR start = 0
        # Assert no miR ends at chr boundaries


class TestExtendMirnas:
    """Test for the 'extend_mirnas' method."""

    def test_extend_mirnas_no_extension(self, ):
        """Test not extending miRNA coordinates."""
        # Add files

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        miR_obj.extend_mirnas(n=0)

        # Assert db == out_db

    def test_extend_mirnas_extend_6_nts(self, ):
        """Test extending miRNA coordinates 6 positions."""
        # Add files

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        miR_obj.extend_mirnas()

        # Assert correct extension


class TestWriteGFF:
    """Test for the 'write_gff' method."""

    def test_write_precursor_file(self, ):
        """Test writing only precursor GFF3."""
        # Add files

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        miR_obj.extend_mirnas()
        miR_obj.write_gff(feature_type="miRNA_primary_transcript")

        # Assert only precursor?

    def test_write_mature_mir_file(self, ):
        """Test writing only mature miRNA GFF3."""
        # Add files

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        miR_obj.extend_mirnas()
        miR_obj.write_gff(feature_type="miRNA")

        # Assert only miRNA?

    def test_write_both_features_file(self, ):
        """Test writing both precursor and mature miRNA GFF3."""
        # Add files

        miR_obj = MirnaExtension()
        miR_obj.set_db()
        miR_obj.set_seq_lengths()
        miR_obj.extend_mirnas()
        miR_obj.write_gff()

        # Assert both features?

class TestExtendMirnas:
    """Test for the 'extend_mirnas' method."""

    def test_extend_mirnas_no_extreme_coords(self, tmp_path, gff_no_extremes):
        """Test miRNA extension with no extreme coordinates."""
        in_file, pre_exp, mir_exp = gff_no_extremes

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        miR_obj = MirnaExtension()
        miR_obj.load_gff_file(str(in_file))
        miR_obj.extend_mirnas(primir_out=primir_out, mir_out=mir_out)

        with open(primir_out, "r") as output, open(pre_exp, "r") as expected:
            assert output.read() == expected.read()

        with open(mir_out, "r") as output, open(mir_exp, "r") as expected:
            assert output.read() == expected.read()

    def test_extend_mirnas_extreme_coords(self, tmp_path, gff_extremes):
        """Test miRNA extension with miRNAs having extreme coordinates."""
        in_file, pre_exp, mir_exp = gff_extremes

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        miR_obj = MirnaExtension()
        miR_obj.load_gff_file(str(in_file))
        miR_obj.extend_mirnas(primir_out=primir_out, mir_out=mir_out)

        with open(primir_out, "r") as output, open(pre_exp, "r") as expected:
            assert output.read() == expected.read()

        with open(mir_out, "r") as output, open(mir_exp, "r") as expected:
            assert output.read() == expected.read()

    def test_extend_mirnas_extreme_coords_chr_boundaries(
        self, tmp_path, gff_extremes_chr
    ):
        """Test miRNA extension with extreme coordinates and chr boundaries."""
        chr_size, in_file, pre_exp, mir_exp = gff_extremes_chr

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        len_dict = {}
        with open(chr_size, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                len_dict[line[0]] = int(line[1])

        miR_obj = MirnaExtension()
        miR_obj.load_gff_file(str(in_file))
        miR_obj.extend_mirnas(
            primir_out=primir_out, mir_out=mir_out, seq_lengths=len_dict
        )

        with open(primir_out, "r") as output, open(pre_exp, "r") as expected:
            assert output.read() == expected.read()

        with open(mir_out, "r") as output, open(mir_exp, "r") as expected:
            assert output.read() == expected.read()


## REINDENT
class TestParseArguments:
    """Test 'parse_arguments()' function."""
    def test_no_files(self, monkeypatch):
        """Call without input nor output files."""
            with pytest.raises(SystemExit) as sysex:
                monkeypatch.setattr(sys, "argv", ["mirna_extension"])
            parse_arguments().parse_args()
            assert sysex.value.code == 2

            def test_in_files(self, monkeypatch, gff_empty, tmp_path):
                """Call with in and output files."""
            gff_in = gff_empty

            monkeypatch.setattr(
                    sys,
                    "argv",
                    [
                        "mirna_extension",
                        str(gff_in),
                        "--outdir",
                        str(tmp_path),
                        ],
                    )

            args = parse_arguments().parse_args()
            assert isinstance(args, argparse.Namespace)

            def test_all_arguments(self, monkeypatch, gff_extremes_chr, tmp_path):
                """Call with all the arguments."""
            chr_size, gff_in, gff_pre_out, gff_mir_out = gff_extremes_chr

            monkeypatch.setattr(
                    sys,
                    "argv",
                    [
                        "mirna_extension",
                        str(gff_in),
                        "--outdir",
                        str(tmp_path),
                        "--chr",
                        str(chr_size),
                        "--extension",
                        "6",
                        ],
                    )

            args = parse_arguments().parse_args()
            assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_file(self, monkeypatch, gff_empty, tmp_path):
        """Test main function with an empty file."""
        gff_empty = gff_empty

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        monkeypatch.setattr(
                sys,
                "argv",
                [
                    "mirna_extension",
                    str(gff_empty),
                    "--outdir",
                    str(tmp_path),
                    ],
                )
        args = parse_arguments().parse_args()
        main(args)

        with open(gff_empty, "r") as expected, open(primir_out, "r") as output:
            assert output.read() == expected.read()

        with open(gff_empty, "r") as expected, open(mir_out, "r") as output:
            assert output.read() == expected.read()

    def test_main_no_extreme_coords(
            self, monkeypatch, tmp_path, gff_no_extremes
            ):
        """Test main function with no extreme coords."""
        in_gff, pre_gff, mir_gff = gff_no_extremes

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        monkeypatch.setattr(
                sys,
                "argv",
                ["mirna_extension", str(in_gff), "--outdir", str(tmp_path)],
                )
        args = parse_arguments().parse_args()
        main(args)

        with open(pre_gff, "r") as expected, open(primir_out, "r") as output:
            assert output.read() == expected.read()

        with open(mir_gff, "r") as expected, open(mir_out, "r") as output:
            assert output.read() == expected.read()

    def test_main_extreme_coords(self, monkeypatch, tmp_path, gff_extremes):
        """Test main function with extreme coords."""
        in_gff, pre_gff, mir_gff = gff_extremes

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        monkeypatch.setattr(
                sys,
                "argv",
                ["mirna_extension", str(in_gff), "--outdir", str(tmp_path)],
                )
        args = parse_arguments().parse_args()
        main(args)

        with open(pre_gff, "r") as expected, open(primir_out, "r") as output:
            assert output.read() == expected.read()

        with open(mir_gff, "r") as expected, open(mir_out, "r") as output:
            assert output.read() == expected.read()

    def test_main_extreme_coords_limit_size(
            self, monkeypatch, tmp_path, gff_extremes_chr
            ):
        """Test main function with extreme coords and limited by chr size."""
        chr_size, in_gff, pre_gff, mir_gff = gff_extremes_chr

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        monkeypatch.setattr(
                sys,
                "argv",
                [
                    "mirna_extension",
                    str(in_gff),
                    "--outdir",
                    str(tmp_path),
                    "--chr",
                    str(chr_size),
                    ],
                )
        args = parse_arguments().parse_args()
        main(args)

        with open(pre_gff, "r") as expected, open(primir_out, "r") as output:
            assert output.read() == expected.read()

        with open(mir_gff, "r") as expected, open(mir_out, "r") as output:
            assert output.read() == expected.read()
