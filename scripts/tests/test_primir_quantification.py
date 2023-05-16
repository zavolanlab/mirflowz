"""Unit tests for module 'primir_quantification.py'."""

import argparse
from pathlib import Path
import sys

import pytest

sys.path.append("../../")

from scripts.primir_quantification import (
    main,
    parse_arguments
)


@pytest.fixture
def empty_file():
    """Import path to empty file."""
    empty_file = Path("files/empty_file")

    return empty_file


@pytest.fixture
def bed_file():
    """Import path to test files with full content."""
    in_bed = Path("files/in_intersection.bed")
    out_table = Path("files/primir_quantification")

    return in_bed, out_table


@pytest.fixture
def bed_no_extension_files():
    """Import path to test files with no extension on features names."""
    in_bed = Path("files/in_intersection_no_extension.bed")
    out_table = Path("files/no_extension_primir_quantification")

    return in_bed, out_table


@pytest.fixture
def bed_extension_id_files():
    """Import path to test files with extension and reads id."""
    in_bed = Path("files/in_intersection.bed")
    out_table = Path("files/extension_id_primir_quantification")

    return in_bed, out_table


@pytest.fixture
def bed_id_files():
    """Import path to test files with reads id."""
    in_bed = Path("files/in_intersection_no_extension.bed")
    out_table = Path("files/id_primir_quantification")

    return in_bed, out_table
    

@pytest.fixture
def bed_some_extension_files():
    """Import path to test files with some extension."""
    in_bed = Path("files/in_intersection_some_extension.bed")
    out_table = Path("files/some_extension_primir_quantification")

    return in_bed, out_table


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without input file."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2

    def test_correct_input(self, monkeypatch, bed_file):
        """Call with the correct input file."""
        in_bed, out_table = bed_file

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)
    
    def test_too_many_input_files(self, monkeypatch, bed_file):
        """Call with too many input file."""
        in_bed, out_table = bed_file
        
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(
                sys, 'argv',
                ['primir_quantification',
                 str(in_bed), str(in_bed),
                 ]
            )
            parse_arguments().parse_args()
        assert sysex.value.code == 2
    
    def test_all_input(self, monkeypatch, bed_file):
        """Call with all the options."""
        in_bed, out_table = bed_file

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--id', "name",
             '--feat-extension',
             '--read-ids',
             '--collapsed',
             ]
        )
        args = parse_arguments().parse_args()
        assert isinstance(args, argparse.Namespace)


class TestMain:
    """Test 'main()' function."""

    def test_main_empty_bed_file(self, monkeypatch, capsys, empty_file):
        """Test main function with an empty bed file."""
        empty_file = empty_file

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(empty_file),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(empty_file, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_no_extension(self, monkeypatch, capsys, bed_no_extension_files):
        """Test main function with no extension in features names."""
        in_bed, expected_out = bed_no_extension_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_id_extension(self, monkeypatch, capsys, bed_extension_id_files):
        """Test main function with extension in feature name and read names."""
        in_bed, expected_out = bed_extension_id_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--feat-extension',
             '--read-ids',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_id(self, monkeypatch, capsys, bed_id_files):
        """Test main function with read names."""
        in_bed, expected_out = bed_id_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--read-ids',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_some_extension_file(self, monkeypatch, capsys, bed_some_extension_files):
        """Test main function with read names."""
        in_bed, expected_out = bed_some_extension_files

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--feat-extension',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()

    def test_main_collpased_file(self, monkeypatch, capsys, bed_file):
        """Test main function with collapsed alignments."""
        in_bed, expected_out = bed_file

        monkeypatch.setattr(
            sys, 'argv',
            ['primir_quantification',
             str(in_bed),
             '--collapsed',
             ]
        )
        args = parse_arguments().parse_args()
        main(args)
        captured = capsys.readouterr()

        with open(expected_out, 'r') as out_file:
            assert captured.out == out_file.read()
