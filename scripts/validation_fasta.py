#!/usr/bin/env python
"""Filter FASTA files."""

import argparse
import gzip
from pathlib import Path
import re
from typing import List, Pattern, TextIO

from Bio import SeqIO, SeqRecord

# ------------------------------------------------------------- #
#   Created: Mar 5, 2019                                        #
#   Author: Paula Iborra                                        #
#   Company: Zavolan Group, Biozentrum, University of Basel     #
# ------------------------------------------------------------- #

def parse_and_validate_arguments():
    """Parse and validate command-line arguments."""
    description = """Process FASTA files.

Process both uncompressed and 'gzip'-compressed FASTA files by trimming,
filtering, and validating sequence records based on user-defined criteria.

Sequence IDs are trimmed at the first occurrence of any specified characters
in '--trim' to standardize naming conventions. If not character string is
provided, the first white space is used.

To filter the FASTA file by sequence IDs, a text file, with one (trimmed) ID
per line, has to be passed to `--filter'. Wheather to keep ('--mode k') or
discard ('--mode d') the sequences with those IDs must be specified.

Sequences exceeding a given length threshold ('--remove') are excluded.

If a path is provided to '--idlist', the resulting sequence IDs are written
one per line in a separate text file.
"""
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "infile",
        help="Path to the input FASTA file",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the output FASTA file",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--trim",
        help=(
            "Character(s) used to trim the ID. Remove anything that follows "
            "the character(s). Default: first white space"
        ),
        type=str,
        nargs="?",
        default="",
    )
    parser.add_argument(
        "--idlist",
        help=(
            "Path to the text file with the final sequences IDs. "
            "One ID per line."
        ),
        type=Path,
    )
    parser.add_argument(
        "-f",
        "--filter",
        help=(
            "Path to the input ID list. Filter FASTA IDs from inpuft ile with "
            "the selected mode. Filter file must contain ONE ID per line."
        ),
        type=Path,
    )
    parser.add_argument(
        "-m",
        "--mode",
        help=(
            "Type of ID filtering for fasta file: keep (k) or discard (d) IDs "
            "contained in the ID list file."
        ),
        choices=("k", "d"),
    )
    parser.add_argument(
        "-r",
        "--remove",
        help="Remove sequences from FASTA file longer than specified length.",
        type=int,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.1.0",
        help="Show program's version number and exit",
    )

    args = parser.parse_args()

    if args.filter and not args.mode:
        parser.error(
            "Mode argument ('--mode', '-m' is required when using the filter"
            " argument ('--filter', '-f'). See '--help' for more information."
        )

    if args.mode and not args.filter:
        parser.error(
            "Filter argument ('--filter', '-f' is required when using the mode"
            " argument ('--mode', '-m'). See '--help' for more information."
        )

    return args


def open_fasta(in_file: Path) -> TextIO:
    """Open a FASTA or FASTA.GZ for text‐mode reading."""
    valid_extension = [".gz", ".fa", ".fasta"]
    suffix = in_file.suffix

    try:
        valid_extension.index(suffix)
    except ValueError as exc:
        raise ValueError(
            f'The provided input file has not a valid extension: "{suffix}".\n'
            "For an uncompressed FASTA file, the valid extensions are '.fa'"
            " and '.fasta'. For a compressed FASTA file, '.gz'."
        ) from exc

    if suffix == ".gz":
        return gzip.open(in_file, "rt")

    return in_file.open("rt", encoding="utf-8")


def write_id_file(out_file: Path, id_list: List[str]) -> None:
    """Write the final sequence IDs, one per line.

    Args:
        out_file: Path to the file where to write the IDs.
        id_list: FASTA IDs to be written.
    """
    with open(out_file, "w", encoding="utf-8") as id_file:
        for seq_id in id_list:
            id_file.write(seq_id + "\n")


def compile_trim_pattern(trim_str: str) -> Pattern[str]:
    """Get a compiled regex pattern to trim at character first occurrence.

    If "trim_str" is empty, white space is used as the default delimiter.

    Args:
        trim_str: Characters used to determine where trimming occurs.

    Returns:
        A compiled regex pattern that captures (1) everything up to the first
        match and (2) the rest of the string.
    """
    if not trim_str:
        new_trim_str = r"\s"

    else:
        new_trim_str = re.escape(trim_str)

    return re.compile(rf"^([^{new_trim_str}]*)(.*)$")


def trim_id(*, seq_rec: SeqRecord, _pattern: Pattern[str]) -> SeqRecord:
    """Trim a FASTA ID using the first-occurrence of any character in _pattern.

    All parameters must be passed by keyword.

    Args:
        seq_rec: A Bio.SeqRecord.SeqRecord to be trimmed in place.
        _pattern: (internal) a pre-compiled regex from "get_trim_pattern".

    Returns:
        The same SeqRecord, with .id and .description possibly updated.
    """
    pattern_match = _pattern.match(seq_rec.id)

    if pattern_match:
        new_id, rest_id = pattern_match.groups()
        seq_rec.id = new_id
        seq_rec.description = ""

    return seq_rec


def main(arguments) -> None:
    """Filter and process a FASTA file."""
    filter_set = None
    out_id_lst: List[str] = []

    trim_pattern = compile_trim_pattern(arguments.trim)

    if arguments.filter:
        with open(arguments.filter, "r", encoding="utf-8") as filt_file:
            filter_set = {line.strip("\n") for line in filt_file}

    with open_fasta(arguments.infile) as in_handle, \
         open(arguments.output, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):

            if trim_pattern:
                record = trim_id(_pattern=trim_pattern, seq_rec=record)

            if filter_set:
                is_in_list = record.id in filter_set

                if (
                    (arguments.mode == "k" and not is_in_list) or
                    (arguments.mode == "d" and is_in_list)
                ):
                    continue

            if arguments.remove and len(record.seq) > arguments.remove:
                continue

            SeqIO.write(record, out_handle, "fasta")
            out_id_lst.append(record.id)

    if arguments.idlist:
        write_id_file(out_file=arguments.idlist, id_list=out_id_lst)


if __name__ == "__main__":
    args = parse_and_validate_arguments()  # pragma:no cover
    main(args)  # pragma: no cover
