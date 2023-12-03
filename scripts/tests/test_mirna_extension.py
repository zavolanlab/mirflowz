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
    """Import path to mirna annotation files with extreme miRNA coords."""
    in_extremes = Path("files/in_mirna_extreme_mirs.gff3")
    out_primir = Path("files/extreme_primir_anno.gff3")
    out_mir = Path("files/extreme_mir_anno.gff3")

    return in_extremes, out_primir, out_mir


@pytest.fixture
def gff_extremes_chr():
    """Import path to mirna annotation files and chr size."""
    chr_size = Path("files/chr_size.txt")
    in_chr_extremes = Path("files/in_mirna_extreme_chr_mirs.gff3")
    out_primir = Path("files/extreme_chr_primir_anno.gff3")
    out_mir = Path("files/extreme_chr_mir_anno.gff3")

    return chr_size, in_chr_extremes, out_primir, out_mir


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


class TestLoadGffFile:
    """Test for the 'load_gff_file' method."""

    def test_load_gff_file(self, gff_no_extremes):
        """Test input loading from file."""
        in_file, pre_exp, mir_exp = gff_no_extremes

        mirnaObject = MirnaExtension()
        mirnaObject.load_gff_file(str(in_file))

        assert mirnaObject is not None
        assert isinstance(mirnaObject.db, gffutils.FeatureDB)
        assert (
            len(
                list(
                    mirnaObject.db.features_of_type("miRNA_primary_transcript")
                )
            )
            == 2
        )
        assert len(list(mirnaObject.db.features_of_type("miRNA"))) == 3

    def test_load_gff_file_std(self, monkeypatch, gff_no_extremes):
        """Test input loading from standard input."""
        in_file, pre_exp, mir_exp = gff_no_extremes
        monkeypatch.setattr(sys, "stdin", str(in_file))

        mirnaObject = MirnaExtension()
        mirnaObject.load_gff_file()

        assert mirnaObject is not None
        assert isinstance(mirnaObject.db, gffutils.FeatureDB)
        assert (
            len(
                list(
                    mirnaObject.db.features_of_type("miRNA_primary_transcript")
                )
            )
            == 2
        )
        assert len(list(mirnaObject.db.features_of_type("miRNA"))) == 3


class TestExtendMirnas:
    """Test for the 'extend_mirnas' method."""

    def test_extend_mirnas_no_extreme_coords(self, tmp_path, gff_no_extremes):
        """Test miRNA extension with no extreme coordinates."""
        in_file, pre_exp, mir_exp = gff_no_extremes

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        mirnaObject = MirnaExtension()
        mirnaObject.load_gff_file(str(in_file))
        mirnaObject.extend_mirnas(primir_out=primir_out, mir_out=mir_out)

        with open(primir_out, "r") as output, open(pre_exp, "r") as expected:
            assert output.read() == expected.read()

        with open(mir_out, "r") as output, open(mir_exp, "r") as expected:
            assert output.read() == expected.read()

    def test_extend_mirnas_extreme_coords(self, tmp_path, gff_extremes):
        """Test miRNA extension with miRNAs having extreme coordinates."""
        in_file, pre_exp, mir_exp = gff_extremes

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        mirnaObject = MirnaExtension()
        mirnaObject.load_gff_file(str(in_file))
        mirnaObject.extend_mirnas(primir_out=primir_out, mir_out=mir_out)

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

        mirnaObject = MirnaExtension()
        mirnaObject.load_gff_file(str(in_file))
        mirnaObject.extend_mirnas(
            primir_out=primir_out, mir_out=mir_out, seq_lengths=len_dict
        )

        with open(primir_out, "r") as output, open(pre_exp, "r") as expected:
            assert output.read() == expected.read()

        with open(mir_out, "r") as output, open(mir_exp, "r") as expected:
            assert output.read() == expected.read()
