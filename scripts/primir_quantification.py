#!/usr/bin/env python

"""Pri-miR quantification.

This script creates the miRNA primary transcript counting table build upon four
columns: the primary transcript name, the count, the 5' end relative extension
and the 3' end relative extension. First, the SAM file will be traversed to
create a dictionary that will have as keys the alignment names and their NH tag
as value. Thereafter, and under the assumption that the BED file is sorted by
the primary transcript name, the counting will be made. For each pri-miRNA, the
counting will consist on the sum of each of its intersecting alignment; the
contribution is computed as 1/NH. Finally, the pri-miRNA name in the BED file
will be split to retrieve the actual pri-miRNA name and the relative 
extensions.

Usage:
    primir_quantification.py --bed BED --sam SAM --outfile outputfile
"""

import argparse
import os
from pathlib import Path

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
        help="Path to the extended pri-miR intersection BED file.",
        type=Path,
        required=True
    )
    parser.add_argument(
        '-s', '--sam',
        help="Path to the SAM file containing the intersecting alignments.",
        type=Path,
        required=True
    )
    parser.add_argument(
        '--outfile',
        help="Path to the output file.",
        type=Path,
        required=True
    )

    return parser

def main(args) -> None:
    """Create pri-miRs counting table."""
    with open(args.outfile, 'w') as outfile:
        header = ["Name", "5p_extension", "3p_extension", "Count"]
        outfile.write("\t".join(header) + "\n")

    if os.path.getsize(args.bed) == 0:
        return
    
    with pysam.AlignmentFile(args.sam, "r") as samfile:
        align_nh =  {}

        for alignment in samfile:
            align_nh[alignment.query_name] = alignment.get_tag("NH")
    
        if len(align_nh.items()) == 0:
            return
    
    with open(args.bed, "r") as bedfile, open(args.outfile, "a") as outfile:

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
                pri_data.append(str(count))
                outfile.write("\t".join(pri_data) + "\n")

                current_name = name
                pri_data = name.split("_")
                count = (1/nh_value)
 
        pri_data.append(str(count))
        outfile.write("\t".join(pri_data) + "\n")

if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
