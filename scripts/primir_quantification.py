#!/usr/bin/env python

"""Pri-miR quantification.

This script creates the miRNA primary transcript counting table build upon four
columns: the primary transcript name, the count, the 5' end relative extension
and the 3' end relative extension. First, the SAM file will be traversed to
create a dictionary that will have as keys the alignment names and their NH tag
as value. Thereafter, and under the assumption that the BED file is sorted by
the primary transcript name, the counting will be made. For each pri-miRNA, the
counting will consist on the sum of each of its intersecting alignment; the
contribution is computed as 1/NH. As for the relative extension fields, if the
extensions are not provided in the primary transcript name separated by an
underscore, a '-0' and a '+0' will be add for the 5' and 3' end relative
extension respectively.

Usage:
    primir_quantification.py --bed BED --sam SAM
    primir_quantification.py --bed BED --sam SAM > outfile
"""

import argparse
from pathlib import Path
import sys

import pysam


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Script to quantify pri-miRs."
        )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s 1.0',
        help="Show program's version number and exit"
    )
    parser.add_argument(
        '-b', '--bed',
        help="Path to the pri-miR BED file. This file must be the \
            output of the call `bedtools intersect -wb -s -F 1 -sorted`.",
        type=Path,
        required=True
    )
    parser.add_argument(
        '-s', '--sam',
        help="Path to the SAM file containing the intersecting alignments.",
        type=Path,
        required=True
    )

    return parser

def main(args) -> None:
    """Create pri-miRs counting table."""

    with open(args.bed, 'r') as bedfile:
        if len(bedfile.read()) == 0:
            return
    
    with pysam.AlignmentFile(args.sam, "r") as samfile:
        align_nh =  {}

        for alignment in samfile:
            align_nh[alignment.query_name] = alignment.get_tag("NH")
    
        if len(align_nh.items()) == 0:
            return
    
    with open(args.bed, "r") as bedfile:

        count = 0
        current_name = None

        for line in bedfile:

            line = line.strip().split("\t")
            name = line[9].split(";")[2].split("=")[1]

            if current_name is None:
                current_name = name
                pri_data = name.split("_")
            
            try:
                nh_value = align_nh[line[13]]
            except KeyError:
                continue

            if current_name == name:
                count += (1/nh_value)
            else:
                pri_data.insert(1, str(count))

                if len(pri_data) < 4:
                    pri_data.extend(['-0', '+0'])

                sys.stdout.write("\t".join(pri_data) + "\n")

                current_name = name
                pri_data = name.split("_")
                count = (1/nh_value)
 
        pri_data.insert(1, str(count))

        if len(pri_data) < 4:
            pri_data.extend(['-0', '+0'])

        sys.stdout.write("\t".join(pri_data) + "\n")

if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
