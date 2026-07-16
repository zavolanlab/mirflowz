"""Unit tests for the module 'copper.py'."""

import argparse
from pathlib import Path
import sys

import pytest

from ..copper import (
    ColorPileup,
    ColorScheme,
    ColorValues,
    ParsedPileupName,
    RowRender,
    create_dir_pileups,
    main,
    parse_arguments,
)


@pytest.fixture
def scheme():
    """Create a sample color scheme."""
    return ColorScheme(
        adenine="green",
        cytosine="orange",
        guanine="light purple",
        thymine="light blue",
        gap="light gray",
        generic="white",
    )


@pytest.fixture
def pileups():
    """Create dictionary with import paths for test input pileup."""
    pileups_dict = {}

    pileups_dict["empty"] = Path(
        "files/pileups/no_shift/test-lib.hsa-miR-516a-3p.0-shift.tab"
    )

    pileups_dict["arm_0_shift"] = Path(
        "files/pileups/no_shift/test-lib.hsa-miR-1323.0-shift.tab"
    )
    pileups_dict["pre_0_shift"] = Path(
        "files/pileups/no_shift/test-lib.hsa-mir-520a.0-shift.tab"
    )

    pileups_dict["arm_5_shift_canon"] = Path(
        "files/pileups/shift/test-lib.hsa-miR-520a-3p.5-shift.tab"
    )
    pileups_dict["pre_3_shift_canon"] = Path(
        "files/pileups/shift/test-lib.hsa-mir-520a.3-shift.tab"
    )

    return pileups_dict


class TestHelperClasses:
    """Test helper record classes."""

    def test_color_values(self):
        """Create a 'ColorValues' object."""
        colors = ColorValues("#FFFFFF", "#000000")

        assert colors.bg == "#FFFFFF"
        assert colors.text == "#000000"

    def test_color_scheme(self, scheme):
        """Create a 'ColorScheme' object."""
        colors = scheme

        assert colors.adenine == "green"
        assert colors.cytosine == "orange"
        assert colors.guanine == "light purple"
        assert colors.thymine == "light blue"
        assert colors.gap == "light gray"
        assert colors.generic == "white"

    def test_parsed_pileup_name(self):
        """Create a 'ParsedPileupName' object."""
        parsed = ParsedPileupName(
            parts=["LIB", "mir", "2-shift"], title="mir", shift=2
        )

        assert parsed.parts[0] == "LIB"
        assert parsed.title == "mir"
        assert parsed.shift == 2

    def test_row_render(self):
        """Create a 'RowRender' object."""
        rendered = RowRender(html="<div></div>", seq_len=4)

        assert rendered.html == "<div></div>"
        assert rendered.seq_len == 4


class TestColorPileup:
    """Test 'ColorPileup' methods."""

    def test_get_color_scheme(self, scheme):
        """Populate the pileup character-to-color mapping."""
        pileup = ColorPileup()

        pileup.get_color_scheme(scheme)

        assert pileup.color_code["-"] == pileup.COLOR_NAME_TO_HEX[scheme.gap]
        assert (
            pileup.color_code["\\."]
            == pileup.COLOR_NAME_TO_HEX[scheme.generic]
        )
        assert (
            pileup.color_code["\\>"]
            == pileup.COLOR_NAME_TO_HEX[scheme.generic]
        )
        assert (
            pileup.color_code["\\<"]
            == pileup.COLOR_NAME_TO_HEX[scheme.generic]
        )
        assert (
            pileup.color_code["A"] == pileup.COLOR_NAME_TO_HEX[scheme.adenine]
        )
        assert (
            pileup.color_code["C"] == pileup.COLOR_NAME_TO_HEX[scheme.cytosine]
        )
        assert (
            pileup.color_code["G"] == pileup.COLOR_NAME_TO_HEX[scheme.guanine]
        )
        assert (
            pileup.color_code["T"] == pileup.COLOR_NAME_TO_HEX[scheme.thymine]
        )

    def test_get_char_seq_no_shift(self, scheme):
        """Render a sequence with no positions flanking the canonical one."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)

        out = pileup.get_char_seq("GA-CAT-A.", shift=0)

        assert out.count('"char-box char-A">A<') == 3

        assert '"char-box char-T">T<' in out
        assert '"char-box char-.">.<' in out
        assert '"char-box char--">-<' in out

    def test_get_char_seq_shift(self, scheme):
        """Render flanking bases as insertions when shift is allowed."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)

        out = pileup.get_char_seq("GA-CAT-A.", shift=3)

        assert out.count('"char-box char--">-<') == 2
        assert out.count('"char-box char--">A<') == 2
        assert out.count('"char-box char--">G<') == 1
        assert out.count('"char-box char-.">.<') == 1

        assert '"char-box char-C">C<' in out
        assert '"char-box char-A">A<' in out
        assert '"char-box char-T">T<' in out

    def test_create_css(self, tmp_path, scheme):
        """Write a CSS file with character-specific styles."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)

        pileup.css = tmp_path / "pileup_style.css"
        pileup.create_css()

        with open(pileup.css, encoding="utf-8") as css_file:
            css = css_file.read()

        assert "word-wrap: break-word;" in css
        assert ".char-A {" in css
        assert "background-color: #D891EF;" in css
        assert ".char-\\. {" in css

    def test_parse_input_name_shift(self, tmp_path):
        """Parse file-name metadata with non-zero shift."""
        pileup = ColorPileup()
        pileup.input = tmp_path / "LIB.hsa-miR-1-5p.3-shift.tab"

        parsed = pileup._parse_input_name()

        assert parsed.parts == ["LIB", "hsa-miR-1-5p", "3-shift"]
        assert parsed.title == "hsa-miR-1-5p &plusmn; 3nt"
        assert parsed.shift == 3

    def test_parse_input_name_no_shift(self, tmp_path):
        """Parse file-name metadata with zero shift."""
        pileup = ColorPileup()
        pileup.input = tmp_path / "LIB.hsa-miR-2-3p.0-shift.tab"

        parsed = pileup._parse_input_name()

        assert parsed.parts == ["LIB", "hsa-miR-2-3p", "0-shift"]
        assert parsed.title == "hsa-miR-2-3p"
        assert parsed.shift == 0

    def test_format_one_field_row_sequence(self, scheme):
        """Render a one-field sequence-only row."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)
        pileup.counts = "Counts"

        rendered = pileup._format_one_field_row(
            field="GACATA", seq_len=0, shift=0
        )

        assert isinstance(rendered, RowRender)
        assert rendered.seq_len == 6
        assert "<strong></strong>" in rendered.html

    def test_format_one_field_row_counts(self, scheme):
        """Render a one-field counts label row."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)
        pileup.counts = "Counts"

        rendered = pileup._format_one_field_row(
            field="Counts", seq_len=3, shift=0
        )

        assert isinstance(rendered, RowRender)
        assert rendered.seq_len == 3
        assert 'char- "> </span>' in rendered.html

    def test_format_two_field_row_keep_info_false(self, scheme):
        """Suppress a non-numeric value when `keep_info` is `False`."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)
        pileup.keep_info = False

        rendered = pileup._format_two_field_row(
            seq="GATACA", vals="name:1-3", shift=0
        )

        assert rendered.seq_len == 6
        assert "name:1-3" not in rendered.html
        assert '<span class="line-data"></span>' in rendered.html

    def test_format_two_field_row_keep_info_true(self, scheme):
        """Keep a non-numeric value when `keep_info` is `True`."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)
        pileup.keep_info = True

        rendered = pileup._format_two_field_row(
            seq="GATACA", vals="name:1-3", shift=0
        )

        assert rendered.seq_len == 6
        assert '<span class="line-data">name:1-3</span>' in rendered.html

    def test_format_three_field_row(self, scheme):
        """Render a three-field row."""
        pileup = ColorPileup()
        pileup.get_color_scheme(scheme)

        rendered = pileup._format_three_field_row(
            seq="GACATA", vals="13", info="canonical", shift=0
        )

        assert rendered.seq_len == 6
        assert "13  canonical" in rendered.html
        assert "<strong>" in rendered.html

    def test_create_html(self, tmp_path, pileups, scheme):
        """Create an HTML file from a pileup input file."""
        pileup = ColorPileup()

        pileup.input = pileups["arm_5_shift_canon"]
        pileup.html = tmp_path / "output.html"
        pileup.css = tmp_path / "pileup_style.css"

        pileup.group_id = "Library"
        pileup.counts = "Counts"
        pileup.keep_info = False

        pileup.get_color_scheme(scheme)
        pileup.create_html()

        with open(pileup.html, encoding="utf-8") as html_file:
            html = html_file.read()

        assert "<title>hsa-miR-520a-3p &plusmn; 5nt</title>" in html
        assert "<h3>Library: test-lib</h3>" in html
        assert "<h2>hsa-miR-520a-3p &plusmn; 5nt</h2>" in html
        assert "Counts" in html
        assert "< hsa-miR-520a-3p" in html
        assert "19:24287-24320:+" not in html
        assert "pileup_style.css" in html


class TestParseArguments:
    """Test 'parse_arguments()' function."""

    def test_no_input(self, monkeypatch):
        """Call without the required input path."""
        with pytest.raises(SystemExit) as sysex:
            monkeypatch.setattr(sys, "argv", ["copper"])

            parse_arguments().parse_args()

        assert sysex.value.code == 2

    def test_only_required_input(self, monkeypatch, pileups):
        """Call with only the required input positional argument."""
        monkeypatch.setattr(sys, "argv", ["copper", str(pileups["empty"])])

        args = parse_arguments().parse_args()

        assert isinstance(args, argparse.Namespace)

        assert args.input == pileups["empty"]
        assert args.group_id == "Library"
        assert args.counts_id == "Counts"
        assert args.keep_info is False

        assert args.adenine == "green"
        assert args.generic == "white"

    def test_all_arguments(self, monkeypatch, tmp_path):
        """Call with all the supported command-line arguments."""
        in_dir = Path("files/pileups")
        out_dir = tmp_path / "colored_pileups"

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "copper",
                str(in_dir),
                "--outdir",
                str(out_dir),
                "--group_id",
                "Transfection",
                "--counts_id",
                "Mean±CPM",
                "--keep_info",
                "-a",
                "pale green",
                "-c",
                "pale orange",
                "-g",
                "pale red",
                "-t",
                "pale blue",
                "-G",
                "gray",
                "--generic",
                "black",
            ],
        )
        args = parse_arguments().parse_args()

        assert isinstance(args, argparse.Namespace)

        assert args.outdir == out_dir
        assert args.group_id == "Transfection"
        assert args.counts_id == "Mean±CPM"
        assert args.keep_info is True

        assert args.adenine == "pale green"
        assert args.cytosine == "pale orange"
        assert args.guanine == "pale red"
        assert args.thymine == "pale blue"
        assert args.gap == "gray"
        assert args.generic == "black"

    def test_invalid_color_choice(self, monkeypatch, pileups, capfd):
        """Reject a color name not present in the predefined palette."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["copper", str(pileups["empty"]), "-a", "not-a-color"],
        )

        with pytest.raises(SystemExit) as sysex:
            parse_arguments().parse_args()

        out, err = capfd.readouterr()

        assert sysex.value.code == 2
        assert "invalid choice" in err


class TestCreateDirPileups:
    """Test 'create_dir_pileups()' function."""

    def test_create_dir_pileups(self, tmp_path, scheme):
        """Create HTML pileups for all tab files in a directory."""
        out_dir = tmp_path / "out_pileups"
        out_dir.mkdir()

        pileup = ColorPileup()
        pileup.css = out_dir / "pileup_style.css"
        pileup.group_id = "Library"
        pileup.counts = "Counts"
        pileup.keep_info = False
        pileup.get_color_scheme(scheme)

        create_dir_pileups(
            in_dir=Path("files/pileups/no_shift"), out_dir=out_dir, obj=pileup
        )

        assert (out_dir / "pileup_style.css").is_file()
        assert (out_dir / "test-lib_hsa-miR-1323_0-shift.html").is_file()
        assert (out_dir / "test-lib_hsa-miR-516a-3p_0-shift.html").is_file()
        assert (out_dir / "test-lib_hsa-mir-520a_0-shift.html").is_file()


class TestMain:
    """Test 'main()' function."""

    def test_main_single_file(self, monkeypatch, pileups, tmp_path):
        """Process a single pileup file."""
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "copper",
                str(pileups["arm_5_shift_canon"]),
                "--outdir",
                str(tmp_path),
                "--group_id",
                "Library",
                "--counts_id",
                "Counts",
            ],
        )

        args = parse_arguments().parse_args()

        main(args)

        assert (tmp_path / "pileup_style.css").is_file()
        assert (tmp_path / "test-lib_hsa-miR-520a-3p_5-shift.html").is_file()

    def test_main_single_directory(self, monkeypatch, tmp_path):
        """Process pileup files from an input directory."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["copper", "files/pileups/no_shift", "--outdir", str(tmp_path)],
        )

        args = parse_arguments().parse_args()

        main(args)

        assert (tmp_path / "pileup_style.css").is_file()
        assert (tmp_path / "test-lib_hsa-miR-516a-3p_0-shift.html").is_file()
        assert (tmp_path / "test-lib_hsa-miR-1323_0-shift.html").is_file()
        assert (tmp_path / "test-lib_hsa-mir-520a_0-shift.html").is_file()

    def test_main_subdirectory(self, monkeypatch, tmp_path):
        """Process pileup files from each input subdirectory independently."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["copper", "files/pileups", "--outdir", str(tmp_path)],
        )

        args = parse_arguments().parse_args()

        main(args)

        assert (tmp_path / "no_shift" / "pileup_style.css").is_file()
        assert (
            tmp_path / "no_shift" / "test-lib_hsa-miR-516a-3p_0-shift.html"
        ).is_file()

        assert (
            tmp_path / "no_shift" / "test-lib_hsa-miR-1323_0-shift.html"
        ).is_file()

        assert (
            tmp_path / "no_shift" / "test-lib_hsa-mir-520a_0-shift.html"
        ).is_file()

        assert (tmp_path / "shift" / "pileup_style.css").is_file()
        assert (
            tmp_path / "shift" / "test-lib_hsa-miR-520a-3p_5-shift.html"
        ).is_file()

        assert (
            tmp_path / "shift" / "test-lib_hsa-mir-520a_3-shift.html"
        ).is_file()
