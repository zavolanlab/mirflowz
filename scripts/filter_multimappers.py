#!/usr/bin/env python

"""Filter multimappers in a SAM file by indel count.

This script uses the pysam library to open the input SAM file and iterates
over each alignment in it. It is assumed that the file is ordered by query
name. For each alignment, it checks if it is a secondary or supplementary
alignment, if it is, it skips to the next alignment. Otherwise, it adds the
alignment to a list storing all alignments with the same query name.
If the current query changes, it selects the best alignment(s) from the list
of current alingments and writes them to the standard output. Current query
and current alignments list are reset to the new query and alignment
respectively.
Finally, if there are any remaining alignments in the current alignments
list and selects the best alignment(s) for those queries, writing them to
standard output.

Input:
    input_sam_file: Path to the input SAM file, ordered by query name.

Output:
    standard output: Filtered SAM file

Functions:
    parse_arguments():
        Command-line arguments parser
    count_indels(pysam.AlignedSegment):
        Count the number of indels in an alignment based on its CIGAR string.
    find_best_alignments(list[pysam.AlignedSegment]):
        Find alignments with less indels
    write_output(list[pysam.AlignedSegment]):
        Write the output to the standard output (stdout).
    main(Path):
        Filter multimappers by indels count.

Usage:
    filter_multimappers.py SAM
    filter_multimappers.py SAM > SAM
"""
import argparse
from pathlib import Path
import sys
from typing import List

import pysam


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Script to filter multimappers by indel counts."
        )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s 1.0',
        help="Show program's version number and exit"
        )

    parser.add_argument(
        'infile',
        help="Path to the SAM input file, sorted by query name.",
        type=Path
        )

    return parser


def count_indels(aln: pysam.AlignedSegment) -> int:
    """Count the number of indels in an alignment based on its CIGAR string.

    This function counts the number of indels in the alignment based on the
    alignment CIGAR string returned by .cigartuples. Insertions are encoded as
    1s whereas the deletions are encoded as 2.
    Refer to the .cigartuples documentation for more information at:
    https://pysam.readthedocs.io/en/v0.15.0/api.html#pysam.AlignedSegment.cigartuples

    Args:
        alignments: The read to count insertions for.


    Returns:
        int: The sum of insertions and deletions in the alignment.
    """
    return sum([op[1] for op in aln.cigartuples if op[0] == 1 or op[0] == 2])


def find_best_alignments(alns: List[pysam.AlignedSegment]) -> List[pysam.AlignedSegment]:
    """Find alignments with less indels.

    This function creates a list of tuples with the alignment object and its
    number of indels. Then, computes the minimum number of indels and returns
    a list with the alignments that have no more than that minimum number of
    indels. In addition, it updates the tag 'NH' and 'HI' to match the final
    number of alignments kept and its identifier respectively.

    Args:
        alignments: alignments with the same query name

    Retrns:
        best_alignments: alignments with the less indels
    """
    aln_indels = [(aln, count_indels(alignment=aln)) for aln in alns]
    min_indels = min(aln_indels, key=lambda x: x[1])[1]
    best_alignments = [aln for i, (aln, indels) in enumerate(aln_indels) if indels == min_indels]

    for i in range(len(best_alignments)):
        best_alignments[i].set_tag('NH', len(best_alignments))
        best_alignments[i].set_tag('HI', i + 1)

    return best_alignments


def write_output(alignments: List[pysam.AlignedSegment]) -> None:
    """Write the output to the standard output (stdout).

    Args:
        alignments: alignments with the same query name
    """
    if len(alignments) == 1:
        sys.stdout.write(alignments[0].to_string() + '\n')
    else:
        best_alignments = find_best_alignments(alignments=alignments)
        for alignment in best_alignments:
            sys.stdout.write(alignment.to_string() + '\n')


def main(sam_file: Path) -> None:
    """Filter multimappers by indels count.

    Args:
        sam_file: Path to the input SAM file.
    """
    with pysam.AlignmentFile(sam_file, "r") as samfile:

        sys.stdout.write(str(samfile.header))

        current_query = None
        current_alignments: list[pysam.AlignedSegment] = []

        for alignment in samfile:

            if alignment.is_secondary or alignment.is_supplementary:
                continue

            if current_query is None:
                current_query = alignment.query_name

            if current_query == alignment.query_name:
                current_alignments.append(alignment)

            else:
                write_output(alignments=current_alignments)

                current_query = alignment.query_name
                current_alignments = [alignment]

        write_output(alignments=current_alignments)


if __name__ == "__main__":
    args = parse_arguments().parse_args()
    main(sam_file=args.infile)
