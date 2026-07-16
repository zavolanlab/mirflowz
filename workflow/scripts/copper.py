#!/usr/bin/env python
"""COlor-coded PileuP gEneratoR."""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, NamedTuple


class ColorValues(NamedTuple):
    """Background and text colors for a pileup symbol."""

    bg: str
    text: str


@dataclass(frozen=True)
class ColorScheme:
    """Background color names used for pileup symbols."""

    adenine: str
    cytosine: str
    guanine: str
    thymine: str
    gap: str
    generic: str


class ParsedPileupName(NamedTuple):
    """Metadata derived from the pileup input file name."""

    parts: list[str]
    title: str
    shift: int


class RowRender(NamedTuple):
    """Rendered HTML fragment and current sequence length."""

    html: str
    seq_len: int


class ColorPileup:
    """Class for code-coloring pileups.

    Attributes:
        input: Path to the ASCII-style alignment pileup file.
        html: Path to the HTML output color-coded pileup file.
        css: Path to the CSS output file.
        group_id: String to categorize samples based on experimental context.
        counts: String used as the counts column title.
        keep_info: Whether to keep the sequence representation name and
            genomic coordinates on the final display.
        color_code: Dictionary mapping pileup characters to `ColorValues`.
    """

    COLOR_NAME_TO_HEX: Dict[str, ColorValues] = {
        "white": ColorValues("#FFFFFF", "#000000"),
        "light blue": ColorValues("#89CFF0", "#000000"),
        "light gray": ColorValues("#E8E8E8", "#000000"),
        "light green": ColorValues("#CAFF70", "#000000"),
        "light green2": ColorValues("#C0FF3E", "#000000"),
        "light orange": ColorValues("#FFA700", "#000000"),
        "light pink": ColorValues("#FF83F4", "#000000"),
        "light pink2": ColorValues("#FFBBFF", "#000000"),
        "light purple": ColorValues("#D891EF", "#000000"),
        "light red": ColorValues("#FF5349", "#FFFFFF"),
        "light yellow": ColorValues("#FFFF99", "#000000"),
        "pale blue": ColorValues("#B0D6D6", "#000000"),
        "pale green": ColorValues("#BDE1B4", "#000000"),
        "pale orange": ColorValues("#E5CE9A", "#000000"),
        "pale pink": ColorValues("#D6B0C3", "#000000"),
        "pale purple": ColorValues("#C3B0D6", "#000000"),
        "pale red": ColorValues("#D6B0B0", "#000000"),
        "pale yellow": ColorValues("#E5E49A", "#000000"),
        "blue": ColorValues("#00BFFF", "#FFFFFF"),
        "blue2": ColorValues("#1CA9C9", "#FFFFFF"),
        "gray": ColorValues("#ADADAD", "#000000"),
        "green": ColorValues("#76EE00", "#000000"),
        "green2": ColorValues("#A6F432", "#000000"),
        "orange": ColorValues("#FFBF00", "#000000"),
        "orange2": ColorValues("#FFAE42", "#000000"),
        "pink": ColorValues("#FF34B3", "#FFFFFF"),
        "pink2": ColorValues("#FCA4DD", "#000000"),
        "purple": ColorValues("#AB82FF", "#000000"),
        "red": ColorValues("#FF0800", "#FFFFFF"),
        "red2": ColorValues("#FF3030", "#FFFFFF"),
        "red3": ColorValues("#E2062C", "#FFFFFF"),
        "yellow": ColorValues("#FFE135", "#000000"),
        "dark blue": ColorValues("#1C39BB", "#FFFFFF"),
        "dark gray": ColorValues("#454545", "#FFFFFF"),
        "dark green": ColorValues("#458B00", "#FFFFFF"),
        "dark green2": ColorValues("#00DD00", "#FFFFFF"),
        "dark orange": ColorValues("#EE7600", "#FFFFFF"),
        "dark pink": ColorValues("#F400A1", "#FFFFFF"),
        "dark pink2": ColorValues("#E4007C", "#FFFFFF"),
        "dark purple": ColorValues("#9400D3", "#FFFFFF"),
        "dark red": ColorValues("#B31B1B", "#FFFFFF"),
        "dark yellow": ColorValues("#FFDF00", "#000000"),
        "black": ColorValues("#000000", "#FFFFFF"),
    }

    def __init__(self) -> None:
        """Initialize class."""
        self.input: Path
        self.html: Path
        self.css: Path

        self.group_id: str
        self.counts: str
        self.keep_info: bool

        self.color_code: Dict[str, ColorValues] = {}

    def get_color_scheme(self, scheme: ColorScheme) -> None:
        """Create the pileup character-to-color mapping.

        Populate `self.color_code` by mapping pileup characters to their
        background and text hexadecimal color codes using the provided color
        scheme.

        Args:
            scheme: Background color names to use for nucleotide, gap, and
                generic pileup symbols.
        """
        self.color_code["-"] = self.COLOR_NAME_TO_HEX[scheme.gap]
        self.color_code["A"] = self.COLOR_NAME_TO_HEX[scheme.adenine]
        self.color_code["C"] = self.COLOR_NAME_TO_HEX[scheme.cytosine]
        self.color_code["G"] = self.COLOR_NAME_TO_HEX[scheme.guanine]
        self.color_code["T"] = self.COLOR_NAME_TO_HEX[scheme.thymine]
        self.color_code["\\."] = self.COLOR_NAME_TO_HEX[scheme.generic]
        self.color_code["\\>"] = self.COLOR_NAME_TO_HEX[scheme.generic]
        self.color_code["\\<"] = self.COLOR_NAME_TO_HEX[scheme.generic]

    def get_char_seq(self, seq: str, shift: int) -> str:
        """Create HTML string for the aligned sequence.

        Bases outside the canonical representation are treated as insertions
        and colored as such.

        Args:
            seq: Pileup character string.
            shift: Number of nucleotides flanking the canonical sequence.

        Returns:
            HTML string with color-coded characters.
        """
        char_seq = ""

        for i, nt in enumerate(seq):

            # Bases preceding the start position
            if i < shift and nt not in [".", " ", ">"]:
                char = f'<span class="char-box char--">{nt}</span>'

            # Bases beyond the end position
            elif i > len(seq) - shift - 1 and nt not in [".", " ", ">"]:
                char = f'<span class="char-box char--">{nt}</span>'

            # Within sequence bases
            else:
                char = f'<span class="char-box char-{nt}">{nt}</span>'

            char_seq += char

        return char_seq

    def create_css(self) -> None:
        """Write the CSS stylesheet used by all generated pileup HTML files.

        This method creates the stylesheet referenced by the HTML pileup
        output. The stylesheet is composed of three sections:

        1. Global layout rules for page headers and container elements.
        2. Per-line formatting rules for each pileup row and its associated
           text.
        3. Per-character styling rules for nucleotide, gap, and generic
           symbols.

        The character-specific CSS classes are generated dynamically from
        `self.color_code`. Each entry defines the background and text color
        for a pileup character, allowing the HTML output to render color-coded
        sequence positions consistently.

        The final stylesheet is written to `self.css`.
        """
        # Global layout and container styles used by every generated HTML page
        css_body = """
h2 {
   font-size: 16px;
   text-align: center;
   font-family: Arial;
   font-weight: bold;
}

h3 {
   font-size: 16px;
   text-align: center;
   font-family: Arial;
   font-weight: bold;
}

.content {
   width: fit-content;
   margin: 0 auto;
   padding: 10px;
   word-wrap: break-word;
   background: #FFFFFF;
}

.pileup {
   padding: 10px;
   margin-left: 80px;
   margin-right: 80px;
   background: #FFFFFF;
   word-wrap: break-word;
   display: block;
}
        """
        # Styles for each rendered pileup row and its associated text values
        css_line = """
.line {
   font-size: 12px;
   line-height: 15px;
}
.line-data {
   margin-left: 5px;
   word-wrap: break-word;
   font-family: Arial;
   font-size: 12px;
}
.line-sequence {
   display: inline-block;
   word-wrap: break-word;
}
        """
        # Base style shared by all individual character boxes
        css_chars = """
.char-box {
   display: inline-block;
   width: 15px;
   height: 15px;
   margin-right: 1px;
   line-height: 15px;
   text-align: center;
   font-family: serif, Courier New;
}
        """
        # Add one CSS class per pileup character using the active color scheme
        for char, colors in self.color_code.items():
            css_chars += f"""
.char-{char} {{
  background-color: {colors.bg};
  color: {colors.text};
}}
        """

        with open(self.css, "w", encoding="UTF-8") as css_f:
            css_f.write("\n".join([css_body, css_line, css_chars]))

    def _parse_input_name(self) -> ParsedPileupName:
        """Extract file-name components, display title, and shift.

        The input file name must be `LIB.MIRNA.#-shift.tab` where:
            - `LIB` is the library name,
            - `MIRNA` is the mature miRNA the pileup is made for
            - `#` is the maximum allowed shift on both ends of the canonical
               sequence
        """
        parts = str(self.input).rsplit("/", maxsplit=1)[-1].split(".")[0:-1]

        title = parts[1]
        shift = parts[-1].split("-")[0]

        if shift != "0":
            title = "".join([title, " &plusmn; ", shift, "nt"])

        return ParsedPileupName(parts=parts, title=title, shift=int(shift))

    def _format_one_field_row(
        self, field: str, seq_len: int, shift: int
    ) -> RowRender:
        """Render a one-field pileup row and return updated sequence length.

        A one-field row is interpreted either as a sequence-only row or as the
        counts label row, depending on whether the field matches `self.counts`.
        """
        if field != self.counts:
            seq_len = len(field)
            seq, vals, info = field, "", ""

        else:
            seq, vals, info = " " * seq_len, " " * seq_len, field

        html_pileup = f"""
            <div class="line">
                <div class="line-sequence">
                    {self.get_char_seq(seq, shift)}
                <span class="line-data">{vals}<strong>{info}</strong></span>
                </div>
            </div>\n
        """

        return RowRender(html=html_pileup, seq_len=seq_len)

    def _format_two_field_row(
        self, seq: str, vals: str, shift: int
    ) -> RowRender:
        """Render a two-field pileup row and return updated sequence length."""
        # Suppress non-numeric secondary values unless extra
        # information should be preserved in the output
        if not self.keep_info:
            try:
                int(vals)
            except ValueError:
                vals = ""

        html_pileup = f"""
            <div class="line">
                <div class="line-sequence">
                    {self.get_char_seq(seq, shift)}
                <span class="line-data">{vals}</span>
                </div>
            </div>\n
        """

        return RowRender(html=html_pileup, seq_len=len(seq))

    def _format_three_field_row(
        self, seq: str, vals: str, info: str, shift: int
    ) -> RowRender:
        """Render a three-field pileup row."""
        html_pileup = f"""
            <div class="line">
                <div class="line-sequence"><strong>
                    {self.get_char_seq(seq, shift)}
                <span class="line-data">{vals}  {info}</span>
                </strong></div>
            </div>\n
        """

        return RowRender(html=html_pileup, seq_len=len(seq))

    def create_html(self) -> None:
        """Create the color-coded pileup HTML file.

        This method reads a single ASCII-style pileup file from `self.input`
        and converts it into an HTML representation written to `self.html`. The
        HTML document includes a header with the sample group and pileup title,
        and a body where each pileup line is rendered as a sequence of
        color-coded character boxes followed by its associated values or
        annotations.

        The pileup title and allowed sequence shift are inferred from the input
        file name. The stylesheet referenced in the HTML output is taken from
        `self.css`. Sequence characters are converted to HTML using
        `self.get_char_seq()`.

        Input lines are interpreted according to their number of tab-delimited
        fields:
            - 1 field: sequence-only lines or the counts label line
            - 2 fields: sequence and value lines
            - 3 fields: sequence, value, and additional annotation lines

        The final HTML document is written to the path stored in `self.html`.
        """
        parsed = self._parse_input_name()
        css_ref = str(self.css).rsplit("/", maxsplit=1)[-1]

        html_header = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{parsed.title}</title>
  <link rel="stylesheet" type="text/css" href="{css_ref}">
</head>
<body>
        """

        html_pileup = f"""
    <div class="content">

        <div class="pileup">
            <h3>{self.group_id}: {parsed.parts[0]}</h3>
            <h2>{parsed.title}</h2>
        """
        with open(self.input, encoding="UTF-8") as in_p:
            seq_len = 0
            for line in in_p:
                fields = line.strip().split("\t")

                # A single-field line is either a sequence-only row or the row
                # used to display the counts label
                if len(fields) == 1:
                    row_rendered = self._format_one_field_row(
                        field=fields[0], seq_len=seq_len, shift=parsed.shift
                    )

                # A two-field line contains a sequence and an associated value
                elif len(fields) == 2:
                    row_rendered = self._format_two_field_row(
                        seq=fields[0], vals=fields[1], shift=parsed.shift
                    )

                # A three-field line contains a sequence, a value, and an
                # additional annotation that should be emphasized in the output
                else:
                    row_rendered = self._format_three_field_row(
                        seq=fields[0],
                        vals=fields[1],
                        info=fields[2],
                        shift=parsed.shift,
                    )

                html_pileup += row_rendered.html
                seq_len = row_rendered.seq_len

        with open(self.html, "w", encoding="UTF-8") as html_f:
            html_f.write("\n".join([html_header, html_pileup]))
            html_f.write("</div>\n</div>\n</body>\n</html>")


def parse_arguments():
    """Parse command-line arguments."""
    description = """Color-code ASCII-style alignment pileups.

    Generate color-coded pileups in HTML and their corresponding CSS style file
    either for a whole directory or a single file.

    If a directory is provided as input, an HTML file is generated for each
    file with just one CSS file.

    For a proper HTML creation, the counts column title must be specified in
    the CLI argument `--counts_id` (see the "Constraints" section for a more
    detailed explanation on the input format)

    The final color-coded pileup has a two-line header:
        - The first line starts with a label used to categorize samples based
            on experimental context specified in `--group_id`, followed by the
            sample name.
        - The second line contains the mature miRNA the pileup is made for,
            followed by the +/- nucleotide shift allowed at either end of the
            miRNA sequence.

    All the nucleotides outside the canonical sequence use the `--gap`
    background. Therefore, only the positions aligning with the canonical
    sequence are colored.

    If the input pileup contains the canonical sequence, the output pileup
    will display it in bold.


    Constraints:
    - The tabulated file has to follow the format specified in
        https://gist.github.com/deliaBlue/a27f78ab9e80c54b9a021e724384c6e6

    The following assumptions are made:
    - The input file name must be `LIB.MIRNA.#-shift.tab` where `LIB` is the
        library name, `MIRNA` is the mature miRNA the pileup is made for, and
        `#` is the maximum allowed shift on both ends of the canonical
        sequence.
    """
    valid_colors = sorted(ColorPileup.COLOR_NAME_TO_HEX)
    parser = argparse.ArgumentParser(
        description=description,
        epilog=(
            "Valid color names: "
            f'{", ".join(valid_colors)}\n\n'
            "Quote values that contain spaces, for example: -a 'light green'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help=(
            "Path to the input directory containing the ASCII-style pileups"
            " or to a single ASCII-style pileup file."
        ),
        type=Path,
        metavar="DIR",
    )
    parser.add_argument(
        "--outdir",
        help="Path where output files are written to. Default: %(default)s.",
        default=Path.cwd(),
        type=Path,
        metavar="DIR",
    )
    parser.add_argument(
        "--group_id",
        help=(
            "String to categorize samples based on experimental context."
            " Default: %(default)s"
        ),
        default="Library",
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "--counts_id",
        help="String used as the counts column title. Default: %(default)s",
        default="Counts",
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "--keep_info",
        help=(
            "Keep genomic coordinates and feature names on the final display."
            " Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-a",
        "--adenine",
        help=(
            "Adenine background color. Must be one of the predefined color "
            "names. Default: %(default)s"
        ),
        default="green",
        choices=valid_colors,
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "-c",
        "--cytosine",
        help=(
            "Cytosine background color. Must be one of the predefined color "
            "names. Default: %(default)s"
        ),
        default="orange",
        choices=valid_colors,
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "-g",
        "--guanine",
        help=(
            "Guanine background color. Must be one of the predefined color "
            "names. Default: %(default)s"
        ),
        default="light purple",
        choices=valid_colors,
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "-t",
        "--thymine",
        help=(
            "Thymine background color. Must be one of the predefined color "
            "names. Default: %(default)s"
        ),
        default="light blue",
        choices=valid_colors,
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "-G",
        "--gap",
        help=(
            "Gap background color. Must be one of the predefined color names. "
            "Default: %(default)s"
        ),
        default="light gray",
        choices=valid_colors,
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "--generic",
        help=(
            "Non-sequence characters background color. Must be one of the "
            "predefined color names. Default: %(default)s"
        ),
        default="white",
        choices=valid_colors,
        type=str,
        metavar="STR",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 2.0",
        help="show program's version and exit.",
    )
    return parser


def create_dir_pileups(in_dir: Path, out_dir: Path, obj: ColorPileup) -> None:
    """Create color-coded pileups for each file in a directory.

    This function assumes that ASCII-style alignment pileup files have the
    extension ".tab".

    Args:
        in_dir: Path to the input directory where pileup files are stored.
        out_dir: Path to the directory where color-coded pileups are written.
        obj: Initialized `ColorPileup` class with the active color scheme.
    """
    obj.create_css()

    files = [_f for _f in in_dir.iterdir() if str(_f).endswith(".tab")]

    for _f in files:
        sample = str(_f).rsplit("/", maxsplit=1)[-1].split(".")[0:-1]

        obj.input = _f
        obj.html = out_dir / f"{'_'.join(sample)}.html"

        obj.create_html()


def main(args):
    """Generate color-coded HTML alignment pileups."""
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Initialize the pileup renderer
    color_pileup = ColorPileup()

    color_pileup.get_color_scheme(
        ColorScheme(
            gap=args.gap,
            adenine=args.adenine,
            cytosine=args.cytosine,
            guanine=args.guanine,
            thymine=args.thymine,
            generic=args.generic,
        )
    )

    color_pileup.keep_info = args.keep_info
    color_pileup.group_id = args.group_id
    color_pileup.counts = args.counts_id

    if args.input.is_dir():
        sub_dirs = [_sdir for _sdir in args.input.iterdir() if _sdir.is_dir()]

        # If the input directory contains subdirectories, process each one
        # independently and create a matching output subdirectory
        if len(sub_dirs) != 0:
            for _sdir in sub_dirs:
                new_out = args.outdir / str(_sdir).rsplit("/", maxsplit=1)[-1]
                new_out.mkdir(parents=True, exist_ok=True)

                color_pileup.css = new_out / "pileup_style.css"

                create_dir_pileups(
                    in_dir=_sdir, out_dir=new_out, obj=color_pileup
                )
        # Otherwise, process all pileup files directly from the input
        # directory and write them to the output directory
        else:
            color_pileup.css = args.outdir / "pileup_style.css"

            create_dir_pileups(
                in_dir=args.input, out_dir=args.outdir, obj=color_pileup
            )
    # Process a single pileup file
    else:
        sample = str(args.input).rsplit("/", maxsplit=1)[-1].split(".")[0:-1]

        color_pileup.input = args.input
        color_pileup.css = args.outdir / "pileup_style.css"
        color_pileup.html = args.outdir / f"{'_'.join(sample)}.html"

        color_pileup.create_css()
        color_pileup.create_html()


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
