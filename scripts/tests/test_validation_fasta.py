"""Unit tests for module 'validation_fasta.py'."""

import argparse
import gzip
from pathlib import Path
import sys
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from ..validation_fasta import (
    compile_trim_pattern,
    main,
    open_fasta,
    parse_and_validate_arguments,
    trim_id,
    write_id_file,
)


@pytest.fixture
def empty_fasta():
    """Import path to an empty FASTA file."""
    return Path("files/empty.fasta")


@pytest.fixture
def fasta_filter_keep():
    """Import paths to FASTA in/out files and to-keep ID list."""
    id_list = Path("files/fasta_id_to_keep.txt")
    in_fasta = Path("files/in_fasta.fa")
    out_fasta = Path("files/filtered_ids.fa")

    return id_list, in_fasta, out_fasta


@pytest.fixture
def fasta_filter_discard():
    """Import paths to FASTA in/out files and to-discard ID list."""
    id_list = Path("files/fasta_id_to_discard.txt")
    in_fasta = Path("files/in_fasta.fa")
    out_fasta = Path("files/filtered_ids.fa")

    return id_list, in_fasta, out_fasta


@pytest.fixture
def fasta_25_len():
    """Import paths to FASTA in/out files with sequence length < 25."""
    in_fasta = Path("files/in_fasta.fa")
    out_fasta = Path("files/filtered_ids.fa")

    return in_fasta, out_fasta


@pytest.fixture
def fasta_id_list():
    """Import paths to FASTA and ID list in/out files."""
    id_list = Path("files/fasta_all_ids.txt")
    in_fasta = Path("files/in_fasta.fa")
    out_fasta = Path("files/all_trimmed_records.fa")

    return id_list, in_fasta, out_fasta


@pytest.fixture
def fasta_trim_default():
    """Import paths to FASTA in/out files with default ID trimming."""
    in_fasta = Path("files/in_fasta.fa")
    out_fasta = Path("files/all_trimmed_records.fa")

    return in_fasta, out_fasta


@pytest.fixture
def fasta_trim_dot():
    """Import paths to FASTA in/out files with ID trimming at dot."""
    in_fasta = Path("files/in_fasta.fa")
    out_fasta = Path("files/all_trimmed_dot_records.fa")

    return in_fasta, out_fasta


@pytest.fixture
def seq_record_dict():
    """Create a dictionary with sample SeqRecord objects."""
    seq_dict = dict()

    # Initial SeqRecord
    seq1 = SeqRecord(
        Seq("AAACCCGGGTTT"),
        id="first.seq",
        description="This is meant to be trimmed"
    )

    # Space-trimmed SeqRecord
    seq2 = SeqRecord(
        Seq("AAACCCGGGTTT"),
        id="first.seq",
        description=""
    )

    # dot-trimmed SeqRecord
    seq3 = SeqRecord(
        Seq("AAACCCGGGTTT"),
        id="first",
        description=""
    )

    seq_dict["full_rec"] = seq1
    seq_dict["space_trimmed"] = seq2
    seq_dict["dot_trimmed"] = seq3

    return seq_dict


class TestOpenFasta:
    """Test 'open_fasta()' function."""

    def test_open_fa(self, tmp_path):
        """Open FASTA file wih '.fa' extension."""
        in_fa = tmp_path/"sample.fa"
        in_fa.write_text(">sample_seq\nACTG\n")

        with open_fasta(in_file=in_fa) as in_f:
            assert ">sample_seq" in [line.strip() for line in in_f.readlines()]

    def test_open_fasta(self, tmp_path):
        """Open FASTA file wih '.fasta' extension."""
        in_fa = tmp_path/"sample.fasta"
        in_fa.write_text(">sample_seq\nACTG\n")

        with open_fasta(in_file=in_fa) as in_f:
            assert "ACTG" in [line.strip() for line in in_f.readlines()]

    def test_open_compressed_fasta(self, tmp_path):
        """Open compressed FASTA."""
        in_comp_fa = tmp_path/"sample.fa.gz"

        with gzip.open(in_comp_fa, "wt") as in_f:
            in_f.write(">sample_seq\nACTG\n")

        with open_fasta(in_file=in_comp_fa) as in_f:
            assert ">sample_seq" in [line.strip() for line in in_f.readlines()]

    def test_open_invalid_extension(self, tmp_path):
        """Open file with an invalid extension."""
        in_txt = tmp_path/"invalid.txt"
        in_txt.write_text("Lorem Ipsum Dolor Sit Amet")

        with pytest.raises(ValueError, match=r".* not a valid extension: ."):
            open_fasta(in_file=in_txt)


class TestCompileTrimPattern:
    """Test 'compile_trim_pattern()' function."""

    def test_compile_default_trim(self):
        """Compile pattern with default trim character (white space)."""
        string = "Lorem Ipsum"

        pattern = compile_trim_pattern(trim_str="")
        match = pattern.match(string)

        assert isinstance(pattern, re.Pattern)
        assert match.group(1) == "Lorem"
        assert match.group(2) == " Ipsum"

    def test_compile_single_trim_char(self):
        """Compile pattern with a single character."""
        string = "Lorem_Ipsum"

        pattern = compile_trim_pattern(trim_str="_")
        match = pattern.match(string)

        assert isinstance(pattern, re.Pattern)
        assert match.group(1) == "Lorem"
        assert match.group(2) == "_Ipsum"

    def test_compile_special_char(self):
        """Compile pattern with a special character."""
        string = "Lorem.Ipsum"

        pattern = compile_trim_pattern(trim_str=".")
        match = pattern.match(string)

        assert isinstance(pattern, re.Pattern)
        assert match.group(1) == "Lorem"
        assert match.group(2) == ".Ipsum"


class TestTrimId:
    """Test 'trim_id()' function."""

    def test_default_trim(self, seq_record_dict):
        """Trim by default character (white space)."""
        rec_dict = seq_record_dict

        i_record = rec_dict["full_rec"]
        o_record = rec_dict["space_trimmed"]

        trim_pat = compile_trim_pattern(trim_str="")
        trimmed_rec = trim_id(seq_rec=i_record, _pattern=trim_pat)

        assert trimmed_rec.id == o_record.id
        assert trimmed_rec.description == o_record.description

    def test_special_char_trim(self, seq_record_dict):
        """Trim by special character (dot)."""
        rec_dict = seq_record_dict

        i_record = rec_dict["full_rec"]
        o_record = rec_dict["dot_trimmed"]

        trim_pat = compile_trim_pattern(trim_str=".")
        trimmed_rec = trim_id(seq_rec=i_record, _pattern=trim_pat)

        assert trimmed_rec.id == o_record.id
        assert trimmed_rec.description == o_record.description

    def test_no_pattern_match(self, seq_record_dict):
        """Trim when no pattern match is found."""
        rec_dict = seq_record_dict

        i_record = rec_dict["full_rec"]

        trim_pat = compile_trim_pattern(trim_str="#")
        trimmed_rec = trim_id(seq_rec=i_record, _pattern=trim_pat)

        assert trimmed_rec.id == i_record.id
        assert trimmed_rec.description == i_record.description


class TestWriteIdFile:
    """Test 'write_id_file()' function."""

    def test_write_empty_id_list(self, tmp_path):
        """Test write an empty list."""
        o_file = tmp_path/"empty_ids.txt"

        write_id_file(out_file=o_file, id_list=[])

        with open(o_file, encoding="utf-8") as output:
            assert output.read() == ""

    def test_write_id_list(self, tmp_path):
        """Test write ID list."""
        o_file = tmp_path/"fasta_id.txt"

        ids = ["first.sequence", "second.sequence"]

        write_id_file(out_file=o_file, id_list=ids)

        with open(o_file, encoding="utf-8") as output:
            assert [line.strip() for line in output.readlines()] == ids

    def test_overwrite_id_list(self, tmp_path):
        """Test overwrite file with ID list."""
        o_file = tmp_path/"fasta_id.txt"
        o_file.write_text("Lorem Ipsum")

        ids = ["first.sequence", "second.sequence"]

        write_id_file(out_file=o_file, id_list=ids)

        with open(o_file, encoding="utf-8") as output:
            assert [line.strip() for line in output.readlines()] == ids


class TestParseAndValidateArguments:
    """Test 'parse_and_validate_arguments()' function."""

    def test_no_files(self, monkeypatch):
        """Call without input nor output files."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(sys, "argv", ["validation_fasta"])
            parse_and_validate_arguments()

        assert sysex.value.code == 2

    def test_required_args(self, monkeypatch, empty_fasta, tmp_path):
        """Call with the required arguments (in/out paths)."""
        in_fasta = empty_fasta

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(in_fasta),
            ],
        )

        args = parse_and_validate_arguments()
        assert isinstance(args, argparse.Namespace)

    def test_filter_no_mode(self, monkeypatch, fasta_id_list, capfd):
        """Call with filter argument and no mode argument."""
        id_list, in_fasta, out_fasta = fasta_id_list

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_fasta),
                "--filter",
                str(id_list),
            ],
        )

        with pytest.raises(SystemExit) as sysex:
            parse_and_validate_arguments()

        out, err = capfd.readouterr()

        assert "Mode argument" in err
        assert sysex.value.code == 2

    def test_mode_no_filter(self, monkeypatch, fasta_id_list, capfd):
        """Call with mode argument and no filter argument."""
        id_list, in_fasta, out_fasta = fasta_id_list

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_fasta),
                "--mode",
                "k",
            ],
        )

        with pytest.raises(SystemExit) as sysex:
            parse_and_validate_arguments()

        out, err = capfd.readouterr()

        assert "Filter argument" in err
        assert sysex.value.code == 2

    def test_invalid_mode_value(self, monkeypatch, fasta_id_list, capfd):
        """Call with invalid mode argument value."""
        id_list, in_fasta, out_fasta = fasta_id_list

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_fasta),
                "--filter",
                str(id_list),
                "--mode",
                "a"
            ],
        )

        with pytest.raises(SystemExit) as sysex:
            parse_and_validate_arguments()

        out, err = capfd.readouterr()

        assert "invalid choice" in err
        assert sysex.value.code == 2

    def test_all_arguments(self, monkeypatch, fasta_id_list, tmp_path):
        """Call with all the arguments."""
        id_list, in_fasta, out_fasta = fasta_id_list

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_fasta),
                "--trim",
                ".",
                "--idlist",
                str(tmp_path),
                "--filter",
                str(id_list),
                "--mode",
                "k",
                "--remove",
                "20",
            ],
        )


class TestMain:
    """Test 'main()' function."""

    def test_main_default(self, monkeypatch, fasta_trim_default, tmp_path):
        """Test main function with default values."""
        in_fasta, out_fasta = fasta_trim_default

        out_tmp = tmp_path/"output.fa"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_tmp),
            ],
        )
        args = parse_and_validate_arguments()
        main(args)

        with open(out_fasta, encoding="utf-8") as expected, open(
            out_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_dot_trim(self, monkeypatch, fasta_trim_dot, tmp_path):
        """Test main function trimming by special character (dot)."""
        in_fasta, out_fasta = fasta_trim_dot

        out_tmp = tmp_path/"output.fa"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_tmp),
                "--trim",
                ".",
            ],
        )
        args = parse_and_validate_arguments()
        main(args)

        with open(out_fasta, encoding="utf-8") as expected, open(
            out_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_filter_keep(self, monkeypatch, fasta_filter_keep, tmp_path):
        """Test main function when filtering by keeping IDs."""
        keep_ids, in_fasta, out_fasta = fasta_filter_keep

        out_tmp = tmp_path/"output.fa"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_tmp),
                "--filter",
                str(keep_ids),
                "--mode",
                "k",
            ],
        )
        args = parse_and_validate_arguments()
        main(args)

        with open(out_fasta, encoding="utf-8") as expected, open(
            out_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_filt_disc(self, monkeypatch, fasta_filter_discard, tmp_path):
        """Test main function when filtering by discarding IDs."""
        discard_ids, in_fasta, out_fasta = fasta_filter_discard

        out_tmp = tmp_path/"output.fa"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_tmp),
                "--filter",
                str(discard_ids),
                "--mode",
                "d",
            ],
        )
        args = parse_and_validate_arguments()
        main(args)

        with open(out_fasta, encoding="utf-8") as expected, open(
            out_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_filt_25_len(self, monkeypatch, fasta_25_len, tmp_path):
        """Test main function when filtering by length < 25nt."""
        in_fasta, out_fasta = fasta_25_len

        out_tmp = tmp_path/"output.fa"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_tmp),
                "--remove",
                "25",
            ],
        )
        args = parse_and_validate_arguments()
        main(args)

        with open(out_fasta, encoding="utf-8") as expected, open(
            out_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

    def test_main_id_list(self, monkeypatch, fasta_id_list, tmp_path):
        """Test main function when output the IDs list."""
        out_ids, in_fasta, out_fasta = fasta_id_list

        out_tmp = tmp_path/"output.fa"
        list_tmp = tmp_path/"id_list.txt"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "validation_fasta",
                str(in_fasta),
                "--output",
                str(out_tmp),
                "--idlist",
                str(list_tmp),
            ],
        )
        args = parse_and_validate_arguments()
        main(args)

        with open(out_fasta, encoding="utf-8") as expected, open(
            out_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()

        with open(out_ids, encoding="utf-8") as expected, open(
            list_tmp, encoding="utf-8"
        ) as output:
            assert output.read() == expected.read()
