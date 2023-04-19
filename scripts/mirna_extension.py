#!/usr/bin/env python

import argparse
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
        help="Path to the tabulated file with the chromosome and its length.",
        type=str,
        default=None
    )

    return parser


def main():
    """Extend pre-miRNAs overhang."""
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

    args = parse_arguments().parse_args()
    main()
