#!/usr/bin/env python

from mirnaExtension_class import MirnaExtension
import argparse


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
        '-e', '--extension',
        help="Minimum number of base pairs each overhang must be of. Maximum value = 10. Default value = 10.",
        default=10,
        type=int,
        choices=range(11)
        )
    parser.add_argument(
        '--chr',
        help="Path to the file containing a tabulated file with the chromosome and its length.",
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

    m = MirnaExtension(gff_file=args.input, n=args.extension, seq_lengths=seq_lengths)
    m.load_gff_file()
    m.extend_primary_mirnas()


if __name__ == "__main__":

    args = parse_arguments().parse_args()
    main()
