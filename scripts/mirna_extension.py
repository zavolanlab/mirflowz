#!/usr/bin/env python

"""Extend miRNAs start and end coordinates by n nucleotides.

This script uses the class MirnaExtension to extend the start and end
coordiantes of the mature miRNAs. Moreover, it will also extend the
corresponding primary transcript start/end coordinates if these are
surpassed by the new start/end coordinates.

Input:
    input_gff_file: 
        Path to the GFF3 annotation file. If not provided, the input will be 
        read from the standard input.
    premir_out_file:
        Path to the primary transcript GFF3 annotation file.
    mir_out_file:
        Path to the mature miRNA GFF3 annotation file.
    extension:
        Number of nucleotides to extend the coordinates. If not provided, the
        extension will be of 6 nucleotides.
    chromosome_size:
        Path to the tabulated file containing the chromosome and its length.
        If not provided, the length will be extracted form the input file.

Output:
    premir_out_file:
        GFF3 file containing only the annotations for the primary transcripts.
    mir_out_file:
        GFF3 file containing only the annotations for the mature miRNAs.

Functions:
    parse_arguments():
        Command-line arguments parser.
    main():
        Extend miRNAs start and end coordinates.

Usage:
    mirna_extension.py [-i GFF3] --premir GFF3 --mir GFF3 [-e int] [--chr Path]
"""

import argparse
import os
import sys

sys.path.append("../")

from scripts.mirnaExtension_class import MirnaExtension


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Script to extend pre-miRNs overhang."
        )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s 1.0',
        help="Show program's version number and exit"
    )
    parser.add_argument(
        '-i', '--input',
        help="Path to the gff3 input file",
        type=str
    )
    parser.add_argument(
        '--premir',
        help="Path to the gff3 pre-miR output file",
        type=str,
        required=True
    )
    parser.add_argument(
        '--mir',
        help="Path to the gff3 miRNA output file.",
        type=str,
        required=True
    )
    parser.add_argument(
        '-e', '--extension',
        help="Number of nucleotides to extend miRs coordinates. Default=6.",
        default=6,
        type=int
    )
    parser.add_argument(
        '--chr',
        help="Path to the tabulated file with the chromosomes and its length.",
        type=str,
        default=None
    )

    return parser


def main(args):
    """Extend miRNAs start/end coordinates."""
    if os.path.getsize(args.input) == 0:
        print("Error: Input file is empty")
        sys.exit(1)

    # Create dictionary with the ref. sequence length
    if args.chr:
        seq_lengths = {}
        with open(args.chr, 'r') as f:
            for line in f:
                ref_seq, length = line.strip().split("\t")
                seq_lengths[ref_seq] = int(length)
    else:
        seq_lengths = None

    m = MirnaExtension(
        gff_file=args.input,
        premir_out=args.premir,
        mir_out=args.mir,
        n=args.extension,
        seq_lengths=seq_lengths
        )
    m.load_gff_file()
    m.extend_mirnas()


if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
