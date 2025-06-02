#!/usr/bin/env python
"""Filter alignments in a SAM file by NH tag."""

import argparse
from pathlib import Path

import pysam


def parse_arguments():
    """Parse command-line arguments."""
    description = """Filter alignments in a SAM file by its NH tag.

For each alignment, check its NH tag, and if the value is higher than the one
specified in `--max_nh`, the aligned read is removed.
"""
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "samfile",
        help="Path to the input SAM file. Rquired!",
        type=Path,
    )
    parser.add_argument(
        "--outfile",
        help="Path to the output SAM file. Required!",
        type=Path,
    )
    parser.add_argument(
        "--max_nh",
        help="Maximum value the NH tag can have for an alignment to be kept.\
              Default: %(default)d.",
        default=100,
        type=int,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.1.0",
        help="Show program's version number and exit",
    )
    return parser


def main(arguments) -> None:
    """Filter alignments by its NH tag value."""
    with (
        pysam.AlignmentFile(arguments.samfile, "r", check_sq=False) as in_sam,
        pysam.AlignmentFile(arguments.outfile, "w", template=in_sam) as o_sam,
    ):
        for alignment in in_sam:
            try:
                nh = alignment.get_tag("NH")

            except KeyError as keyerr:
                raise KeyError(
                    "Missing NH tag: Some alignments do not have the NH"
                    " tag. Please, check that all entries have an"
                    " associated NH tag"
                ) from keyerr

            if nh > arguments.max_nh:
                pass

            o_sam.write(alignment)


if __name__ == "__main__":
    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
