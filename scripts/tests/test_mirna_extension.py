"""Unit tests for module 'mirna_extension.py'."""

import argparse
import sys

from pathlib import Path
import gffutils  # type: ignore
import pytest

from ..mirna_extension import main, MirnaExtension, parse_arguments


@pytest.fixture
def gff_empty():
    """Import path to empty test file."""
    empty = Path("files/empty_file")

    return empty


@pytest.fixture
def gff_replicas():
    """Import path to miRNA annotation files with replicas."""
    in_replica = Path("files/in_replica_mirna_anno.gff3")
    out_replica = Path("files/replica_mirna_anno.gff3")

    return in_replica, out_replica


@pytest.fixture
def gff_diff_strand():
    """Import path to miRNA annotation files."""
    in_diff_strand = Path("files/in_mirna_anno.gff3")
    out_mir = Path("files/mir_anno.gff3")

    return in_diff_strand, out_mir


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
    """Import path to sequence lengths table."""
    correct_tbl = Path("files/chr_size.tsv")

    return correct_tbl


class TestSetDb:
    """Test for the 'set_db' method."""

    def test_set_db_file(self, gff_diff_strand):
        """Test setting local db from file."""
        in_file, exp_out = gff_diff_strand

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        assert miR_obj is not None
        assert isinstance(miR_obj.db, gffutils.FeatureDB)
        assert (
            len(list(miR_obj.db.features_of_type("miRNA_primary_transcript")))
            == 2
        )
        assert len(list(miR_obj.db.features_of_type("miRNA"))) == 4

    def test_set_db_empty_file(self, gff_empty):
        """Test setting local db from empty file."""
        in_file = gff_empty

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        assert miR_obj.db is None


class TestSetSeqLengths:
    """Test for the 'set_seq_lengths' method."""

    def test_set_lengths_no_tbl(self, gff_diff_strand):
        """Test create sequence lengths dictionary from GFF3 file."""
        in_file, exp_out = gff_diff_strand

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)
        miR_obj.set_seq_lengths()

        assert len(miR_obj.seq_lengths.keys()) == 1
        assert miR_obj.seq_lengths["19"] == 121102

    def test_set_lengths_wrong_type_tbl(
        self, tmp_path, gff_extremes, seq_len_tbl
    ):
        """Test create sequence lengths dictionary from wrong table."""
        in_file, pre_exp, mir_exp = gff_extremes

        table = tmp_path / "files/tbl.tsv"
        table.parent.mkdir()
        table.touch()
        table.write_text("19\t600000bp")

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        with pytest.raises(ValueError, match=r".* integer .*"):
            miR_obj.set_seq_lengths(path=table)

    def test_set_lengths_wrong_len_tbl(
        self, tmp_path, gff_extremes_chr, seq_len_tbl
    ):
        """Test create sequence lengths dictionary from wrong table."""
        in_file, pre_exp, mir_exp = gff_extremes_chr

        table = tmp_path / "files/tbl.tsv"
        table.parent.mkdir()
        table.touch()
        table.write_text("19\t315700")

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        with pytest.raises(Exception, match=r".* exceeds .*"):
            miR_obj.set_seq_lengths(path=table)

    def test_set_lengths_correct_tbl(self, gff_extremes_chr, seq_len_tbl):
        """Test create sequence lengths dictionary from correct table."""
        in_file, pre_exp, mir_exp = gff_extremes_chr
        correct_tbl = seq_len_tbl

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)
        miR_obj.set_seq_lengths(path=correct_tbl)

        assert len(miR_obj.seq_lengths.keys()) == 1


class TestAdjustNames:
    """Test for the 'adjust_names' method."""

    def test_adjust_names_no_replicas(self, gff_replicas):
        """Test adjusting names when no replicas are present."""
        in_file, out_file = gff_replicas
        out_mir = ["hsa-miR-10401-3p", "hsa-miR-10401-5p"]

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        precursor = miR_obj.db["MI0033425"]
        matures = [miR_obj.db["MIMAT0041633"], miR_obj.db["MIMAT0041634"]]

        miR_obj.adjust_names(precursor=precursor, matures=matures)

        assert precursor.attributes["Name"][0] == "hsa-mir-10401"
        assert matures[0].attributes["Name"][0] in out_mir
        assert matures[1].attributes["Name"][0] in out_mir

    def test_adjust_names_id_replica(self, gff_replicas):
        """Test adjusting names when replica integer is in the ID."""
        in_file, out_file = gff_replicas
        out_mir = ["hsa-miR-10401-2-3p", "hsa-miR-10401-2-5p"]

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        precursor = miR_obj.db["MI0033425_2"]
        matures = [miR_obj.db["MIMAT0041633_1"], miR_obj.db["MIMAT0041634_1"]]

        miR_obj.adjust_names(precursor=precursor, matures=matures)

        assert precursor.attributes["Name"][0] == "hsa-mir-10401-2"
        assert matures[0].attributes["Name"][0] in out_mir
        assert matures[1].attributes["Name"][0] in out_mir

    def test_adjust_names_name_replica(self, gff_replicas):
        """Test adjusting names when replica integer is in the Name."""
        in_file, out_file = gff_replicas
        out_mir = ["hsa-miR-16-1-3p", "hsa-miR-16-1-5p"]

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        precursor = miR_obj.db["MI0000070"]
        matures = [miR_obj.db["MIMAT0000069"], miR_obj.db["MIMAT0004489"]]

        miR_obj.adjust_names(precursor=precursor, matures=matures)

        assert precursor.attributes["Name"][0] == "hsa-mir-16-1"
        assert matures[0].attributes["Name"][0] in out_mir
        assert matures[1].attributes["Name"][0] in out_mir

    def test_adjust_names_single(self, gff_replicas):
        """Test adjusting names when a single mature miR is present."""
        in_file, out_file = gff_replicas

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        precursor = miR_obj.db["MI0005764"]
        matures = [miR_obj.db["MIMAT0004984_1"]]

        miR_obj.adjust_names(precursor=precursor, matures=matures)

        assert precursor.attributes["Name"][0] == "hsa-mir-941-2"
        assert matures[0].attributes["Name"][0] == "hsa-miR-941-2"

    def test_adjust_names_replace_replica(self, gff_replicas):
        """Test adjusting names when mature miRs have the replica integer."""
        in_file, out_file = gff_replicas
        out_mir = ["hsa-miR-16-2-3p", "hsa-miR-16-2-5p"]

        miR_obj = MirnaExtension()
        miR_obj.set_db(path=in_file)

        precursor = miR_obj.db["MI0000115"]
        matures = [miR_obj.db["MIMAT0000069"], miR_obj.db["MIMAT0004518"]]

        miR_obj.adjust_names(precursor=precursor, matures=matures)

        assert precursor.attributes["Name"][0] == "hsa-mir-16-2"
        assert matures[0].attributes["Name"][0] in out_mir
        assert matures[1].attributes["Name"][0] in out_mir


class TestProcessPrecursor:
    """Test for the 'process_precursor' method."""

    def test_process_prec_diff_strand_same_coords(
        self,
        gff_diff_strand,
        seq_len_tbl,
    ):
        """Test processing precursor on different strands, similar coords."""
        in_file, out_file = gff_diff_strand
        correct_tbl = seq_len_tbl

        exp_mir_obj = MirnaExtension()
        exp_mir_obj.set_db(out_file)

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(correct_tbl)

        out_pre_1 = miR_obj.process_precursor(
            precursor=miR_obj.db["MI0000779"], n=6
        )
        out_pre_2 = miR_obj.process_precursor(
            precursor=miR_obj.db["MI0017393"], n=6
        )

        assert out_pre_1[0] == exp_mir_obj.db["MI0000779"]
        assert out_pre_2[0] == exp_mir_obj.db["MI0017393"]

        for mir in out_pre_1[1:]:
            assert mir in list(
                exp_mir_obj.db.region(
                    seqid="19",
                    strand="+",
                    start=121034,
                    end=121104,
                    featuretype="miRNA",
                )
            )
        for mir in out_pre_2[1:]:
            assert mir in list(
                exp_mir_obj.db.region(
                    seqid="19",
                    strand="-",
                    start=121032,
                    end=121102,
                    featuretype="miRNA",
                )
            )

    def test_process_prec_miR_out_seq_boundaries(
        self,
        gff_extremes_chr,
        seq_len_tbl,
    ):
        """Test processing precursor for outside seq boundaries miRNAs."""
        in_file, pre_out, mir_out = gff_extremes_chr
        correct_tbl = seq_len_tbl

        exp_pre_obj = MirnaExtension()
        exp_pre_obj.set_db(pre_out)

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(correct_tbl)

        out_pre_1 = miR_obj.process_precursor(
            precursor=miR_obj.db["MI0005757"], n=6
        )
        out_pre_2 = miR_obj.process_precursor(
            precursor=miR_obj.db["MI0003140"], n=6
        )

        assert out_pre_1[0] == exp_pre_obj.db["MI0005757"]
        assert out_pre_2[0] == exp_pre_obj.db["MI0003140"]
        assert out_pre_1[1].end == miR_obj.seq_lengths["19"]
        assert out_pre_2[1].start == 1

    def test_process_prec_unknown_seq(self, gff_extremes, tmp_path):
        """Test processing precursor with annotated seq ID not in len dict."""
        in_file, pre_out, mir_out = gff_extremes

        table = tmp_path / "files/tbl.tsv"
        table.parent.mkdir()
        table.touch()
        table.write_text("13\t600000")

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(table)

        with pytest.raises(KeyError, match=r".* not available .*"):
            miR_obj.process_precursor(precursor=miR_obj.db["MI0005757"])

    def test_process_prec_no_precursor_extension(self, gff_extremes):
        """Test processing precursor without precursor extension."""
        in_file, pre_out, mir_out = gff_extremes

        exp_pre_obj = MirnaExtension()
        exp_mir_obj = MirnaExtension()
        exp_pre_obj.set_db(pre_out)
        exp_mir_obj.set_db(mir_out)

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths()

        out_pre = miR_obj.process_precursor(
            precursor=miR_obj.db["MI0003140"], n=2
        )

        assert out_pre[0] == exp_pre_obj.db["MI0003140"]
        assert out_pre[1:] == list(
            exp_mir_obj.db.region(
                seqid="19", start=6, end=124, featuretype="miRNA"
            )
        )


class TestUpdateDb:
    """Test for the 'update_db' method."""

    def test_update_db_no_extension_no_names(self, gff_diff_strand):
        """Test not extending miRNA coordinates nor changing names."""
        in_file, out_file = gff_diff_strand

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths()
        miR_obj.update_db(n=0)

        for mir in miR_obj.db.features_of_type("miRNA"):
            exp_mir = miR_obj.db_out[mir.id]
            assert mir.attributes["Name"][0] == exp_mir.attributes["Name"][0]

    def test_update_db_no_extension_unique_names(self, gff_replicas):
        """Test not extending miRNA coordinates but changing names."""
        in_file, out_file = gff_replicas

        exp_mir_obj = MirnaExtension()
        exp_mir_obj.set_db(out_file)

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths()
        miR_obj.update_db(n=0)

        for mir in miR_obj.db_out.features_of_type("miRNA"):
            exp_mir = exp_mir_obj.db[mir.id]
            assert mir.attributes["Name"][0] == exp_mir.attributes["Name"][0]

    def test_update_db_extend_6_nts_no_names(
        self, gff_diff_strand, seq_len_tbl
    ):
        """Test extending miRNA coordinates 6 nts but not changing names."""
        in_file, out_file = gff_diff_strand
        correct_tbl = seq_len_tbl

        exp_mir_obj = MirnaExtension()
        exp_mir_obj.set_db(out_file)

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(correct_tbl)
        miR_obj.update_db(n=6)

        for mir in miR_obj.db_out.all_features():
            exp_mir = exp_mir_obj.db[mir.id]
            assert mir.astuple() == exp_mir.astuple()

    def test_update_db_extend_6_nts_unique_names(self, gff_replicas):
        """Test extending miRNA coordinates 6 nts and changing names."""
        in_file, out_file = gff_replicas

        exp_mir_obj = MirnaExtension()
        exp_mir_obj.set_db(out_file)

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths()
        miR_obj.update_db(n=6)

        for mir in miR_obj.db_out.features_of_type("miRNA"):
            exp_mir = exp_mir_obj.db[mir.id]
            assert mir.attributes["Name"][0] == exp_mir.attributes["Name"][0]


class TestWriteGFF:
    """Test for the 'write_gff' method."""

    def test_write_empty_file(self, gff_empty, tmp_path):
        """Test writing an empty file."""
        empty_file = gff_empty

        empty_out = tmp_path / "empty.gff3"

        miR_obj = MirnaExtension()
        miR_obj.set_db(empty_file)
        miR_obj.set_seq_lengths()
        miR_obj.update_db()
        miR_obj.write_gff(path=empty_out)

        with open(empty_out, encoding="utf-8") as output, open(
            empty_file, encoding="utf-8"
        ) as expected:
            assert output.read() == expected.read()

    def test_write_precursor_file(self, gff_extremes, seq_len_tbl, tmp_path):
        """Test writing only precursor GFF3."""
        in_file, pre_exp, mir_exp = gff_extremes
        correct_tbl = seq_len_tbl

        pre_out = tmp_path / "extended_premir_annotation_2_nt.gff3"

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(correct_tbl)
        miR_obj.update_db(n=2)
        miR_obj.write_gff(
            path=pre_out, feature_type="miRNA_primary_transcript"
        )

        with open(pre_out, encoding="utf-8") as output, open(
            pre_exp, encoding="utf-8"
        ) as expected:
            assert output.read() == expected.read()

    def test_write_mature_mir_file(self, gff_extremes, seq_len_tbl, tmp_path):
        """Test writing only mature miRNA GFF3."""
        in_file, pre_exp, mir_exp = gff_extremes
        correct_tbl = seq_len_tbl

        mir_out = tmp_path / "extended_mir_annotation_2_nt.gff3"

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(correct_tbl)
        miR_obj.update_db(n=2)
        miR_obj.write_gff(path=mir_out, feature_type="miRNA")

        with open(mir_out, encoding="utf-8") as output, open(
            mir_exp, encoding="utf-8"
        ) as expected:
            assert output.read() == expected.read()

    def test_write_both_features_file(
        self, tmp_path, gff_diff_strand, seq_len_tbl
    ):
        """Test writing both precursor and mature miRNA GFF3."""
        in_file, exp_out = gff_diff_strand
        correct_tbl = seq_len_tbl

        out_file = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        miR_obj = MirnaExtension()
        miR_obj.set_db(in_file)
        miR_obj.set_seq_lengths(correct_tbl)
        miR_obj.update_db()
        miR_obj.write_gff(path=out_file)

        with open(out_file, encoding="utf-8") as output, open(
            exp_out, encoding="utf-8"
        ) as expected:
            assert output.read() == expected.read()


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

    def test_all_arguments(
        self, monkeypatch, gff_extremes_chr, seq_len_tbl, tmp_path
    ):
        """Call with all the arguments."""
        gff_in, gff_pre_out, gff_mir_out = gff_extremes_chr
        correct_tbl = seq_len_tbl

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "mirna_extension",
                str(gff_in),
                "--outdir",
                str(tmp_path),
                "--chr",
                str(correct_tbl),
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
        in_empty = gff_empty

        primir_out = tmp_path / "extended_primir_annotation_6_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_6_nt.gff3"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "mirna_extension",
                str(in_empty),
                "--outdir",
                str(tmp_path),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(in_empty, encoding="utf-8") as expected, open(
            primir_out, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

        with open(in_empty, encoding="utf-8") as expected, open(
            mir_out, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_no_extreme_coords(self, monkeypatch, tmp_path, gff_extremes):
        """Test main function with no extreme coords."""
        in_gff, pre_gff, mir_gff = gff_extremes

        primir_out = tmp_path / "extended_primir_annotation_2_nt.gff3"
        mir_out = tmp_path / "extended_mirna_annotation_2_nt.gff3"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "mirna_extension",
                str(in_gff),
                "--outdir",
                str(tmp_path),
                "--extension",
                "2",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(pre_gff, encoding="utf-8") as expected, open(
            primir_out, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

        with open(mir_gff, encoding="utf-8") as expected, open(
            mir_out, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_extreme_coords_limit_size(
        self, monkeypatch, tmp_path, gff_extremes_chr, seq_len_tbl
    ):
        """Test main function with extreme coords and limited by chr size."""
        in_gff, pre_gff, mir_gff = gff_extremes_chr
        correct_tbl = seq_len_tbl

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
                str(correct_tbl),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)

        with open(pre_gff, encoding="utf-8") as expected, open(
            primir_out, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

        with open(mir_gff, encoding="utf-8") as expected, open(
            mir_out, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()
