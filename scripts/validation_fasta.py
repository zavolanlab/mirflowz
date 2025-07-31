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


with f:
    record = []
    nrec = -1
    inseq = 0

    sys.stdout.write("Parsing FASTA file...")

    for line in f:
        if re.match(r"^>", line):
            nrec += 1
            record.append(Seq())

            # define id of the record
            if not args.trim:
                mobj = re.match(r"^>(\S*)(.*)", line)
            else:
                mobj = re.match(f"^>([^{args.trim}]*)(.*)", line)

            # add id and features
            if mobj:
                record[nrec].id = mobj.group(1)
                record[nrec].features = mobj.group(2)
            inseq = 0
        else:
            if inseq == 0:
                inseq = 1
                record[nrec].seq = line
            else:
                cstring = record[nrec].seq + line
                record[nrec].seq = cstring

    sys.stdout.write("DONE\n")

# ID FILTER LIST #
id_filter = []
if args.filter:
    sys.stdout.write("Filtering FASTA file...")

    with open(args.filter, encoding="utf-8") as filter_file:
        id_filter = [line.rstrip("\n") for line in filter_file]

    sys.stdout.write("DONE\n")


# OUTPUT FASTA FILE #

if args.output:
    sys.stdout.write("Writing FASTA file...")

    with open(args.output, "w", encoding="utf-8") as output:
        if args.filter and args.mode == "k":
            if args.remove:
                for x in range(0, nrec + 1):
                    if (
                        record[x].id in id_filter
                        and len(record[x].seq) - 1 <= args.remove
                    ):
                        output.write(f">{record[x].id}\n{record[x].seq}")
            else:
                for x in range(0, nrec + 1):
                    if record[x].id in id_filter:
                        output.write(f">{record[x].id}\n{record[x].seq}")

        elif args.filter and args.mode == "d":
            if args.remove:
                for x in range(0, nrec + 1):
                    if (
                        record[x].id not in id_filter
                        and len(record[x].seq) - 1 <= args.remove
                    ):
                        output.write(f">{record[x].id}\n{record[x].seq}")
            else:
                for x in range(0, nrec + 1):
                    if record[x].id not in id_filter:
                        output.write(f">{record[x].id}\n{record[x].seq}")

        else:
            if args.remove:
                for x in range(0, nrec + 1):
                    if len(record[x].seq) - 1 <= args.remove:
                        output.write(f">{record[x].id}\n{record[x].seq}")
            else:
                for x in range(0, nrec + 1):
                    output.write(f">{record[x].id}\n{record[x].seq}")
    sys.stdout.write("DONE\n")
