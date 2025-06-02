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


def main():
    """Filter alignments by NH tag."""
    sys.stdout.write(
        f"Removing reads aligned more than {sys.argv[2]} times... \n"
    )

    infile = pysam.Samfile(sys.argv[1], "r", check_sq=False)
    out = pysam.Samfile(sys.argv[3], "w", template=infile)

    keep = True

    for DNAread in infile.fetch():
        intags = DNAread.tags

        for entry in intags:
            if "NH" in entry and entry[1] > int(sys.argv[2]):
                keep = False
        if keep:
            out.write(DNAread)

        keep = True

    out.close()
    sys.stdout.write("DONE!\n")


if __name__ == "__main__":
    main()
