"""Unit tests for module 'iso_name_tagging.py'."""

import argparse
from pathlib import Path
import sys

import pysam
import pytest

from ..iso_name_tagging import (
    attributes_dictionary,
    get_tags,
    main,
    parse_arguments,
    parse_intersect_output,
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
    bed_file = Path("files/in_intersection_extended_mirna.bed")
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


# pylint: disable=line-too-long
@pytest.fixture
def attributes_gtf_gff():
    """Attributes strings in GFF3 and GTF format."""
    gff_attributes = "ID=MIMAT0019022;Alias=MIMAT0019022;Name=hsa-miR-4488;Derives_from=MI0016849"  # noqa: E501
    gtf_attributes = 'ID "MIMAT0019022";Alias "MIMAT0019022";Name "hsa-miR-4488";Derives_from "MI0016849"'  # noqa: E501
    dict_attributes = {
        "id": "MIMAT0019022",
        "alias": "MIMAT0019022",
        "name": "hsa-miR-4488",
        "derives_from": "MI0016849",
    }

    return gff_attributes, gtf_attributes, dict_attributes


# pylint: enable=line-too-long


@pytest.fixture
def alns():
    """Create sample AlignedSegment objects."""
    aln1 = pysam.AlignedSegment()
    aln1.cigarstring = "21M"
    aln1.query_name = "8-2_1"
    aln1.query_sequence = "GAAGGCGCTTCCCTTTGGAGT"
    aln1.reference_id = 19
    aln1.reference_start = 44414
    aln1.set_tag("MD", 21)
    aln1.set_tag("NH", 2)
    aln1.set_tag("HI", 1)

    aln2 = pysam.AlignedSegment()
    aln2.cigarstring = "21M"
    aln2.query_name = "8-2_1"
    aln2.query_sequence = "GAAGGCGCTTCCCTTTGGAGT"
    aln2.reference_id = 19
    aln2.reference_start = 44409
    aln2.set_tag("MD", 21)
    aln2.set_tag("NH", 2)
    aln2.set_tag("HI", 1)

    return aln1, aln2


@pytest.fixture
def intersect_dicts():
    """Dictionaries of the intersect output BED files."""
    inter_name = {
        "22-1_1": [("hsa-miR-1323", 5338, 5359)],
        "41-1_1": [("hsa-miR-1323", 5338, 5359)],
        "44-1_1": [("hsa-miR-1323", 5338, 5359)],
        "46-1_1": [("hsa-miR-1323", 5338, 5359)],
        "48-1_1": [("hsa-miR-1323", 5338, 5359)],
        "66-1_1": [("hsa-miR-1323", 5338, 5359)],
        "92-1_1": [("hsa-miR-1323", 5338, 5359)],
        "14-1_1": [("hsa-miR-498-5p", 7590, 7612)],
        "15-1_1": [("hsa-miR-498-5p", 7590, 7612)],
        "20-1_1": [("hsa-miR-498-5p", 7590, 7612)],
        "34-1_1": [("hsa-miR-498-5p", 7590, 7612)],
        "39-1_1": [("hsa-miR-498-5p", 7590, 7612)],
        "54-1_1": [("hsa-miR-498-5p", 7590, 7612)],
        "26-1_1": [("hsa-miR-498-3p", 7626, 7648)],
        "31-1_1": [("hsa-miR-498-3p", 7626, 7648)],
        "32-1_1": [("hsa-miR-498-3p", 7626, 7648)],
        "49-1_1": [("hsa-miR-498-3p", 7626, 7648)],
        "57-1_1": [("hsa-miR-498-3p", 7626, 7648)],
        "71-1_1": [("hsa-miR-498-3p", 7626, 7648)],
        "13-1_1": [("hsa-miR-524-5p", 44377, 44398)],
        "33-1_1": [("hsa-miR-524-5p", 44377, 44398)],
        "56-1_1": [("hsa-miR-524-5p", 44377, 44398)],
        "59-1_1": [("hsa-miR-524-5p", 44377, 44398)],
        "72-1_1": [("hsa-miR-524-5p", 44377, 44398)],
        "94-1_1": [("hsa-miR-524-5p", 44377, 44398)],
        "8-2_1": [("hsa-miR-524-3p", 44414, 44434)],
        "24-1_1": [("hsa-miR-524-3p", 44414, 44434)],
        "58-1_1": [("hsa-miR-524-3p", 44414, 44434)],
        "90-1_1": [("hsa-miR-524-3p", 44414, 44434)],
        "95-1_1": [("hsa-miR-524-3p", 44414, 44434)],
    }
    inter_id_alias = {
        "22-1_1": [("MIMAT0005795", 5338, 5359)],
        "41-1_1": [("MIMAT0005795", 5338, 5359)],
        "44-1_1": [("MIMAT0005795", 5338, 5359)],
        "46-1_1": [("MIMAT0005795", 5338, 5359)],
        "48-1_1": [("MIMAT0005795", 5338, 5359)],
        "66-1_1": [("MIMAT0005795", 5338, 5359)],
        "92-1_1": [("MIMAT0005795", 5338, 5359)],
        "14-1_1": [("MIMAT0002824", 7590, 7612)],
        "15-1_1": [("MIMAT0002824", 7590, 7612)],
        "20-1_1": [("MIMAT0002824", 7590, 7612)],
        "34-1_1": [("MIMAT0002824", 7590, 7612)],
        "39-1_1": [("MIMAT0002824", 7590, 7612)],
        "54-1_1": [("MIMAT0002824", 7590, 7612)],
        "26-1_1": [("MIMAT0037323", 7626, 7648)],
        "31-1_1": [("MIMAT0037323", 7626, 7648)],
        "32-1_1": [("MIMAT0037323", 7626, 7648)],
        "49-1_1": [("MIMAT0037323", 7626, 7648)],
        "57-1_1": [("MIMAT0037323", 7626, 7648)],
        "71-1_1": [("MIMAT0037323", 7626, 7648)],
        "13-1_1": [("MIMAT0002849", 44377, 44398)],
        "33-1_1": [("MIMAT0002849", 44377, 44398)],
        "56-1_1": [("MIMAT0002849", 44377, 44398)],
        "59-1_1": [("MIMAT0002849", 44377, 44398)],
        "72-1_1": [("MIMAT0002849", 44377, 44398)],
        "94-1_1": [("MIMAT0002849", 44377, 44398)],
        "8-2_1": [("MIMAT0002850", 44414, 44434)],
        "24-1_1": [("MIMAT0002850", 44414, 44434)],
        "58-1_1": [("MIMAT0002850", 44414, 44434)],
        "90-1_1": [("MIMAT0002850", 44414, 44434)],
        "95-1_1": [("MIMAT0002850", 44414, 44434)],
    }
    inter_extension = {
        "22-1_1": [("hsa-miR-1323", 5341, 5356)],
        "41-1_1": [("hsa-miR-1323", 5341, 5356)],
        "44-1_1": [("hsa-miR-1323", 5341, 5356)],
        "46-1_1": [("hsa-miR-1323", 5341, 5356)],
        "48-1_1": [("hsa-miR-1323", 5341, 5356)],
        "66-1_1": [("hsa-miR-1323", 5341, 5356)],
        "92-1_1": [("hsa-miR-1323", 5341, 5356)],
        "14-1_1": [("hsa-miR-498-5p", 7593, 7609)],
        "15-1_1": [("hsa-miR-498-5p", 7593, 7609)],
        "20-1_1": [("hsa-miR-498-5p", 7593, 7609)],
        "34-1_1": [("hsa-miR-498-5p", 7593, 7609)],
        "39-1_1": [("hsa-miR-498-5p", 7593, 7609)],
        "54-1_1": [("hsa-miR-498-5p", 7593, 7609)],
        "26-1_1": [("hsa-miR-498-3p", 7629, 7645)],
        "31-1_1": [("hsa-miR-498-3p", 7629, 7645)],
        "32-1_1": [("hsa-miR-498-3p", 7629, 7645)],
        "49-1_1": [("hsa-miR-498-3p", 7629, 7645)],
        "57-1_1": [("hsa-miR-498-3p", 7629, 7645)],
        "71-1_1": [("hsa-miR-498-3p", 7629, 7645)],
        "13-1_1": [("hsa-miR-524-5p", 44380, 44395)],
        "33-1_1": [("hsa-miR-524-5p", 44380, 44395)],
        "56-1_1": [("hsa-miR-524-5p", 44380, 44395)],
        "59-1_1": [("hsa-miR-524-5p", 44380, 44395)],
        "72-1_1": [("hsa-miR-524-5p", 44380, 44395)],
        "94-1_1": [("hsa-miR-524-5p", 44380, 44395)],
        "8-2_1": [("hsa-miR-524-3p", 44417, 44431)],
        "24-1_1": [("hsa-miR-524-3p", 44417, 44431)],
        "58-1_1": [("hsa-miR-524-3p", 44417, 44431)],
        "90-1_1": [("hsa-miR-524-3p", 44417, 44431)],
        "95-1_1": [("hsa-miR-524-3p", 44417, 44431)],
    }

    return inter_name, inter_id_alias, inter_extension


class TestAttributesDictionary:
    """Test 'attributes_dictionary()' function."""

    def test_attributes_dictionary_gff(self, attributes_gtf_gff):
        """Test get attributes dictionary from GFF3 attributes string."""
        gff_attr, gtf_attr, dict_attr = attributes_gtf_gff

        assert dict_attr == attributes_dictionary(gff_attr)

    def test_attributes_dictionary_gtf(self, attributes_gtf_gff):
        """Test get attributes dictionary from GTF attributes string."""
        gff_attr, gtf_attr, dict_attr = attributes_gtf_gff

        assert dict_attr == attributes_dictionary(gtf_attr)


class TestParseIntersectOutput:
    """Test 'parse_intersect_output()' function."""

    def test_parse_intersect_output_empty_bed(self, empty_files):
        """Test parse intersect empty output."""
        empty_bed, empty_sam = empty_files

        assert parse_intersect_output(intersect_file=empty_bed) is None

    def test_parse_intersect_output_name_id(self, bed_sam, intersect_dicts):
        """Test parse intersect output using 'name' as ID."""
        in_bed, in_sam, out_f = bed_sam
        out_name, out_id_alias, out_ext_3 = intersect_dicts

        out_inter = parse_intersect_output(intersect_file=in_bed)

        assert sorted(out_inter.items()) == sorted(out_name.items())

    def test_parse_intersect_output_alias_id(self, bed_sam, intersect_dicts):
        """Test parse intersect output using 'alias' as ID."""
        in_bed, in_sam, out_f = bed_sam
        out_name, out_id_alias, out_ext_3 = intersect_dicts

        out_inter = parse_intersect_output(intersect_file=in_bed, ID="alias")

        assert sorted(out_inter.items()) == sorted(out_id_alias.items())

    def test_parse_intersect_output_id_id(self, bed_sam, intersect_dicts):
        """Test parse intersect output using 'id' as ID."""
        in_bed, in_sam, out_f = bed_sam
        out_name, out_id_alias, out_ext_3 = intersect_dicts

        out_inter = parse_intersect_output(intersect_file=in_bed, ID="id")

        assert sorted(out_inter.items()) == sorted(out_id_alias.items())

    def test_parse_intersect_output_extend_3(self, bed_sam, intersect_dicts):
        """Test parse intersect output with an extension of 3 nucleotides."""
        in_bed, in_sam, out_f = bed_sam
        out_name, out_id_alias, out_ext_3 = intersect_dicts

        out_inter = parse_intersect_output(intersect_file=in_bed, extension=3)

        assert sorted(out_inter.items()) == sorted(out_ext_3.items())


class TestGetTags:
    """Test 'get_tags()' function."""

    def test_get_tags_within_shift_range_1_extn(self, alns, intersect_dicts):
        """Test get tags within shift range +/-1."""
        aln_1, aln_2 = alns
        out_name, out_id_alias, out_ext_3 = intersect_dicts

        tags = get_tags(
            intersecting_feat=out_name[aln_1.query_name],
            alignment=aln_1,
            shift=1,
        )
        expected = "hsa-miR-524-3p|1|1|21M|21|GAAGGCGCTTCCCTTTGGAGT"

        assert next(iter(tags)) == expected

    def test_get_tags_out_shift_range(self, alns, intersect_dicts):
        """Test get tags outside shift range +/-0."""
        aln_1, aln_2 = alns
        out_name, out_id_alias, out_ext_3 = intersect_dicts

        tags = get_tags(
            intersecting_feat=out_name[aln_1.query_name],
            alignment=aln_2,
            shift=0,
        )
        assert tags == set()


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_bed(self, monkeypatch, bed_sam):
        """Call without bed file."""
        in_bed, in_sam, output = bed_sam

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "iso_name_tagging",
                    "--sam",
                    str(in_sam),
                ],
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_no_sam(self, monkeypatch, bed_sam):
        """Call without bed file."""
        in_bed, in_sam, output = bed_sam

        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    "iso_name_tagging",
                    "--bed",
                    str(in_bed),
                ],
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, bed_sam):
        """Call with the correct input files."""
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(in_bed),
                "--sam",
                str(in_sam),
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)

    def test_all_input(self, monkeypatch, bed_sam):
        """Call with all the options."""
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(in_bed),
                "--sam",
                str(in_sam),
                "--id",
                "alias",
                "--extension",
                "6",
                "--shift",
                "6",
            ],
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_bed_file(
        self, monkeypatch, capsys, empty_files, bed_sam
    ):
        """Test main function with an empty bed file."""
        empty_bed, empty_sam = empty_files

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(empty_bed),
                "--sam",
                str(empty_sam),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_sam, encoding="utf-8") as out_file:
            assert captured.out == out_file.read()

    def test_main_empty_sam_file(
        self, monkeypatch, capsys, empty_files, bed_sam
    ):
        """Test main function with an empty sam file."""
        empty_bed, empty_sam = empty_files
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(in_bed),
                "--sam",
                str(empty_sam),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_sam, encoding="utf-8") as out_file:
            assert captured.out == out_file.read()

    def test_main_bed_sam_file(self, monkeypatch, capsys, bed_sam):
        """Test main function without options."""
        in_bed, in_sam, output = bed_sam

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(in_bed),
                "--sam",
                str(in_sam),
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, encoding="utf-8") as out_file:
            assert captured.out == out_file.read()

    def test_main_bed_sam_extension_file(
        self, monkeypatch, capsys, bed_sam_extension
    ):
        """Test main function with extension and allowed shit equal to 6."""
        in_bed, in_sam, output = bed_sam_extension

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(in_bed),
                "--sam",
                str(in_sam),
                "--extension",
                "6",
                "--shift",
                "6",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, encoding="utf-8") as out_file:
            assert captured.out == out_file.read()

    def test_main_bed_sam_file_id(self, monkeypatch, capsys, bed_sam_id):
        """Test main function with id equals id."""
        in_bed, in_sam, output = bed_sam_id

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "iso_name_tagging",
                "--bed",
                str(in_bed),
                "--sam",
                str(in_sam),
                "--id",
                "id",
            ],
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(output, encoding="utf-8") as out_file:
            assert captured.out == out_file.read()
