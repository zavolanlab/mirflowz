"""Unit tests for module 'validate_bedtools_intersect.py'."""

from pathlib import Path
from itertools import islice

import pytest

from ..validate_bedtools_intersect import (
    FileFormatError,
    Record,
    validate_first_n,
    parse_all,
)


@pytest.fixture
def valid_files():
    """Import path to test files with a valid content."""
    file_8_lines = Path("files/valid_8_lines.intersect")
    file_13_lines = Path("files/valid_13_lines.intersect")

    return file_8_lines, file_13_lines


@pytest.fixture
def invalid_files():
    """Import path to test files with an invalid line."""
    file_8_lines = Path("files/invalid_8_lines.intersect")
    file_13_lines = Path("files/invalid_13_lines.intersect")

    return file_8_lines, file_13_lines


def make_line(sep: str = "\t", **mods) -> str:
    """Build a 'bedtools intersect -wo -s' record, with optional field updates.

    Args:
        sep: Delimiter used to merge all fields into a single string
            (default: TAB).
        **mods: Field -> Value replacements; keys must exist in Record.FIELDS

    Returns: A single str line
    """
    base_rec = {
        "feat_chr": "chr0",
        "feat_source": "X_source",
        "feat_type": "miRNA",
        "feat_start": "1298",
        "feat_end": "1324",
        "feat_score": ".",
        "feat_strand": "+",
        "feat_phase_frame": ".",
        "feat_attrs": "ID=MIR0000001;Alias=MIR0000001;Derives_from=MIR0000000",
        "read_chr": "chr0",
        "read_start": "1299",
        "read_end": "1312",
        "read_name": "read1",
        "read_score": ".",
        "read_strand": "+",
        "overlap_len": "13",
    } | mods

    sorted_parts = [base_rec[field] for field in Record.FIELDS]

    return sep.join(sorted_parts)


class TestFromLine:
    """Test for the 'from_line' class method."""

    def test_from_CSV_line(self):
        """Test parsing a valid line for a CSV file."""
        rec = Record.from_line(make_line(sep=","), sep=",")

        assert isinstance(rec, Record)

        assert rec.feat_chr == "chr0" and rec.read_chr == "chr0"
        assert rec.feat_score == "." and rec.read_score == "."
        assert rec.feat_strand == "+" and rec.read_strand == "+"

        assert rec.feat_source == "X_source"
        assert rec.feat_type == "miRNA"
        assert rec.feat_start == 1298
        assert rec.feat_end == 1324
        assert rec.feat_phase_frame == "."

        assert rec.feat_attrs["id"] == "MIR0000001"
        assert rec.feat_attrs["alias"] == "MIR0000001"
        assert rec.feat_attrs["derives_from"] == "MIR0000000"

        assert rec.read_start == 1299
        assert rec.read_end == 1312
        assert rec.read_name == "read1"

        assert rec.overlap_len == 13

    def test_from_gff3_line(self):
        """Test parsing a valid line with GFF3-style feature attributes."""
        rec = Record.from_line(make_line())

        assert isinstance(rec, Record)

        assert rec.feat_chr == "chr0" and rec.read_chr == "chr0"
        assert rec.feat_score == "." and rec.read_score == "."
        assert rec.feat_strand == "+" and rec.read_strand == "+"

        assert rec.feat_source == "X_source"
        assert rec.feat_type == "miRNA"
        assert rec.feat_start == 1298
        assert rec.feat_end == 1324
        assert rec.feat_phase_frame == "."

        assert rec.feat_attrs["id"] == "MIR0000001"
        assert rec.feat_attrs["alias"] == "MIR0000001"
        assert rec.feat_attrs["derives_from"] == "MIR0000000"

        assert rec.read_start == 1299
        assert rec.read_end == 1312
        assert rec.read_name == "read1"

        assert rec.overlap_len == 13

    def test_from_gtf_line(self):
        """Test parsing a valid line with GTF-style feature attributes."""
        rec = Record.from_line(
            make_line(feat_attrs='gene_id "feature"; gene_name "feat_gtf";')
        )

        assert isinstance(rec, Record)

        assert rec.feat_chr == "chr0" and rec.read_chr == "chr0"
        assert rec.feat_score == "." and rec.read_score == "."
        assert rec.feat_strand == "+" and rec.read_strand == "+"

        assert rec.feat_source == "X_source"
        assert rec.feat_type == "miRNA"
        assert rec.feat_start == 1298
        assert rec.feat_end == 1324
        assert rec.feat_phase_frame == "."

        assert rec.feat_attrs["gene_id"] == "feature"
        assert rec.feat_attrs["gene_name"] == "feat_gtf"

        assert rec.read_start == 1299
        assert rec.read_end == 1312
        assert rec.read_name == "read1"

        assert rec.overlap_len == 13

    def test_from_line_less_fields(self):
        """Test line with less fields than required."""
        line = make_line()
        invalid_line = "\t".join(line.split("\t")[:-1])

        with pytest.raises(
            FileFormatError, match=r"expected 16 columns, found 15"
        ):
            Record.from_line(invalid_line)

    def test_from_line_more_fields(self):
        """Test line with more fields than required."""
        line = make_line()
        invalid_line = line + "\t" + line

        with pytest.raises(
            FileFormatError, match=r"expected 16 columns, found 32"
        ):
            Record.from_line(invalid_line)


class TestCrossCheck:
    """Test for cross_checks method at Record's construction."""

    def test_cross_check_diff_chr(self):
        """Test with feature and read having different chromosomes."""
        with pytest.raises(FileFormatError, match=r"chromosome mismatch: .*"):
            Record(
                feat_chr="chr0",
                feat_source="X_source",
                feat_type="miRNA",
                feat_start=1298,
                feat_end=1324,
                feat_score=".",
                feat_strand="+",
                feat_phase_frame=0,
                feat_attrs={"id": "MIMAT1", "alias": "MIMAT1"},
                read_chr="chr1",
                read_start=1299,
                read_end=1312,
                read_name="read0",
                read_score=".",
                read_strand="+",
                overlap_len=13,
            )

    def test_cross_check_diff_strand(self):
        """Test with feature and read having different strands."""
        with pytest.raises(FileFormatError, match=r"strand mismatch: .*"):
            Record(
                feat_chr="chr0",
                feat_source="X_source",
                feat_type="miRNA",
                feat_start=1298,
                feat_end=1324,
                feat_score=".",
                feat_strand="-",
                feat_phase_frame=0,
                feat_attrs={"id": "MIMAT1", "alias": "MIMAT1"},
                read_chr="chr0",
                read_start=1299,
                read_end=1312,
                read_name="read0",
                read_score=".",
                read_strand="+",
                overlap_len=13,
            )

    def test_cross_check_invalid_feature_range(self):
        """Test with an invalid feature range."""
        with pytest.raises(FileFormatError, match=r"feature range is .*"):
            Record(
                feat_chr="chr0",
                feat_source="X_source",
                feat_type="miRNA",
                feat_start=1324,
                feat_end=1289,
                feat_score=".",
                feat_strand="+",
                feat_phase_frame=0,
                feat_attrs={"id": "MIMAT1", "alias": "MIMAT1"},
                read_chr="chr0",
                read_start=1299,
                read_end=1312,
                read_name="read0",
                read_score=".",
                read_strand="+",
                overlap_len=13,
            )

    def test_cross_check_invalid_read_range(self):
        """Test with an invalid read range."""
        with pytest.raises(FileFormatError, match=r"read range is invalid:.*"):
            Record(
                feat_chr="chr0",
                feat_source="X_source",
                feat_type="miRNA",
                feat_start=1298,
                feat_end=1324,
                feat_score=".",
                feat_strand="+",
                feat_phase_frame=0,
                feat_attrs={"id": "MIMAT1", "alias": "MIMAT1"},
                read_chr="chr0",
                read_start=1312,
                read_end=1299,
                read_name="read0",
                read_score=".",
                read_strand="+",
                overlap_len=13,
            )

    def test_cross_check_no_overlap(self):
        """Test with feature and read not overlapping."""
        with pytest.raises(FileFormatError, match=r"computed overlap is 0 .*"):
            Record(
                feat_chr="chr0",
                feat_source="X_source",
                feat_type="miRNA",
                feat_start=1298,
                feat_end=1324,
                feat_score=".",
                feat_strand="+",
                feat_phase_frame=0,
                feat_attrs={"id": "MIMAT1", "alias": "MIMAT1"},
                read_chr="chr0",
                read_start=12,
                read_end=13,
                read_name="read0",
                read_score=".",
                read_strand="+",
                overlap_len=13,
            )

    def test_cross_check_diff_overlap(self):
        """Test with an overlap different from the annotated one."""
        with pytest.raises(FileFormatError, match=r"overlap length .*"):
            Record(
                feat_chr="chr0",
                feat_source="X_source",
                feat_type="miRNA",
                feat_start=1298,
                feat_end=1324,
                feat_score=".",
                feat_strand="+",
                feat_phase_frame=0,
                feat_attrs={"id": "MIMAT1", "alias": "MIMAT1"},
                read_chr="chr0",
                read_start=1299,
                read_end=1312,
                read_name="read0",
                read_score=".",
                read_strand="+",
                overlap_len=14,
            )


class TestStaticParsers:
    """Test for Record's static parsing methods."""

    def test_parse_text(self):
        """Test parsing text field."""
        assert Record._text("text") == "text"

    def test_parse_valid_as_int(self):
        """Test parsing a correct integer field."""
        assert Record._as_int("1", "field") == 1

    def test_parse_invalid_as_int(self):
        """Test parsing a invalid integer field."""
        with pytest.raises(FileFormatError, match=r".* must be an integer.*"):
            Record._as_int("invalid int", "field")

    def test_parse_valid_as_float_or_dot(self):
        """Test parsing a float as 'float or dot' field."""
        assert Record._as_float_or_dot("3.14", "field") == 3.14

    def test_parse_invalid_as_float_or_dot(self):
        """Test parsing an invalid float or dot field."""
        with pytest.raises(FileFormatError, match=r".* must be a float .*"):
            Record._as_float_or_dot("text", "field")

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_valid_strand_as_strand(self, strand):
        """Test parsing a valid strand field."""
        assert Record._as_strand(strand, "field") == strand

    def test_invalid_strand(self):
        """Test parsing an invalid strand field."""
        with pytest.raises(FileFormatError, match=r".* strand must be .*"):
            Record._as_strand("text", "field")

    @pytest.mark.parametrize("i,o", [(".", "."), ("0", 0), ("1", 1), ("2", 2)])
    def test_valid_as_phase_frame(self, i, o):
        """Test parsing a valid phase/frame field."""
        assert Record._as_phase_frame(i, "field") == o

    @pytest.mark.parametrize("phase_frame", ["3", "text"])
    def test_invalid_phase_frame(self, phase_frame):
        """Test parsing an invalid phase/frame field."""
        with pytest.raises(FileFormatError, match=r".*phase/frame must be .*"):
            Record._as_phase_frame(phase_frame, "field")

    def test_parse_gff3_feat_attrs(self):
        """Test parsing feature GFF3-style attributes."""
        attr_gff = "ID=feature;Alias=feat_gff3;Derives_from=parent"
        exp_attr = {
            "id": "feature",
            "alias": "feat_gff3",
            "derives_from": "parent",
        }

        assert Record._parse_feat_attrs(attr_gff, "field") == exp_attr

    def test_parse_gft_feat_attrs(self):
        """Test parsing feature GFT-style attributes."""
        attr_gtf = 'gene_id "feature"; gene_name "feat_gtf"'
        exp_attr = {"gene_id": "feature", "gene_name": "feat_gtf"}

        assert Record._parse_feat_attrs(attr_gtf, "field") == exp_attr

    def test_parse_invalid_feat_attrs(self):
        """Test parsing invalid attributes style."""
        attr_inv = "This is; an invalid; attributes line"

        with pytest.raises(FileFormatError, match=r".* must be GFF3 .*"):
            Record._parse_feat_attrs(attr_inv, "field")

    def test_parse_duplicate_keys_feat_attrs(self):
        """Test parsing attributes with duplicate keys."""
        attr_dup = 'gene_id "feat1"; gene_id "feat2"; gene_name "feat_gtf"'
        # Assuming last value wins for duplicate keys
        exp_attr = {"gene_id": "feat2", "gene_name": "feat_gtf"}
        assert Record._parse_feat_attrs(attr_dup, "field") == exp_attr


class TestComputeOverlap:
    """Test for the '_compute_overlap' static method."""

    def test_compute_negative_overlap(self):
        """Test for non-overlapping coordinates with result < 0."""
        overlap = Record._compute_overlap(
            f_start=30, f_end=40, r_start=10, r_end=20
        )

        assert overlap == 0

    def test_compute_no_overlap(self):
        """Test for non-overlapping coordinates with result == 0."""
        overlap = Record._compute_overlap(
            f_start=30, f_end=40, r_start=40, r_end=50
        )

        assert overlap == 0

    def test_compute_positive_overlap(self):
        """Test for overlapping coordinates."""
        overlap = Record._compute_overlap(
            f_start=30, f_end=60, r_start=40, r_end=50
        )

        assert overlap == 10


class TestValidateFirstN:
    """Test for the 'validate_first_n()' function."""

    def test_validate_first_n_valid_file(self, valid_files):
        """Test validate first 12 lines with correct file content."""
        short_file, long_file = valid_files

        validate_first_n(long_file, n=12)

    def test_validate_first_n_short_file(self, valid_files):
        """Test validate first 10 lines with shorter valid file content."""
        short_file, long_file = valid_files

        validate_first_n(short_file)

    def test_validate_first_n_invalid_file(self, invalid_files):
        """Test validate first 10 lines with line 3 being invalid."""
        short_file, long_file = invalid_files

        with pytest.raises(
            FileFormatError,
            match=r"Invalid format in line 3: strand mismatch: .*",
        ):
            validate_first_n(short_file)


class TestParseAll:
    """Test for the parse_all()' function."""

    def test_parse_correct_file(self, valid_files):
        """Test parse a whole file with valid content."""
        short_file, long_file = valid_files

        out = list(parse_all(short_file))

        assert [n for n, _ in out] == list(range(1, 9))
        assert all(isinstance(r, Record) for _, r in out)

    def test_parse_incorrect_file(self, invalid_files):
        """Test parse a whole file with an invalid line."""
        short_file, long_file = invalid_files

        it_out = parse_all(long_file)

        nums, recs = zip(*islice(it_out, 11))
        assert nums == tuple(range(1, 12))
        assert all(isinstance(r, Record) for r in recs)

        with pytest.raises(
            FileFormatError,
            match=r"Invalid format in line 12:.*chromosome mismatch: .*"
        ):
            next(it_out)
