#!/usr/bin/env python

"""Collapse alignments in a SAM file by query sequence and starting position.

Collapse alignments in a SAM file by query sequence and starting position, and
write the collapsed alignment to the standard output. The process is divided
into two steps: first, alignments with the same sequence and starting position
are grouped into a list. Second, the alignments in each list are combined into
a single alignment. The new alignment is given a name starting with "set_"
followed by the number of alignments collapsed, and the "HI" tag is reset to 1.
Two new tags are also added: the "XC" tag is the count of collapsed alignments
divided by the "NH" tag, and the "ZS" tag contains a list of all the alignment
IDs that are combined to create the new one.
"""

import argparse
from pathlib import Path
import sys

import pysam


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description=__doc__
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


def collapse_alignments(alns: list[pysam.AlignedSegment], aln: int) -> pysam.AlignedSegment:
    """Collapse alignments list into a single alignment.

    Combine the alignments in the list into a single alignment. The new
    alignment is given a name starting with "set_" followed by the number of
    alignments collapsed, and the "HI" tag is reset to 1. Two new tags are also
    added: the "XC" tag is the count of collapsed alignments divided by the
    "NH" tag. If the "NH" tag is missing, the value 1 is used. The "ZS" tag
    contains a list of all the alignment IDs that are combined to create
    the new one. If the list contains a single alignment, only the 
    aggregation of the new tags is made.

    Args:
        alns: alignments with the same query sequence
        aln: integer to define the new alignment

    Retrns:
        new_alignment: collapsed alignment
    """
    new_alignment = alns[0]

    try:
        nh_tag = new_alignment.get_tag("NH")
    except KeyError:
        nh_tag = 1

    if len(alns) == 1:
        xc_count = 1 / nh_tag
        new_alignment.set_tags([("NH", nh_tag,),
                                ("XC", xc_count), 
                                ("ZS", new_alignment.query_name)])

    else:
        xc_count = len(alns) / nh_tag
        read_ids = sorted({aln.query_name for aln in alns})
        new_alignment.query_name = f'set_{len(read_ids)}.{aln}'
        new_alignment.set_tags([("NH", nh_tag),
                                ("XC", xc_count), 
                                ("HI", 1), 
                                ("ZS", "_".join(read_ids))])

    return new_alignment


def main(sam_file: Path) -> None:
    """Collapse alignments by query sequence and starting position.

    Args:
        sam_file: Path to the input SAM file.
    """
    with pysam.AlignmentFile(sam_file, 'r') as samfile:

        sys.stdout.write(str(samfile.header))

        try:
            aln = next(samfile)
            current_seq = f'{aln.query_sequence}_{aln.reference_start}'
            current_alignments = [aln]
            count = 1

        except StopIteration:
            return

        for aln in samfile:

            if current_seq == f'{aln.query_sequence}_{aln.reference_start}':
                current_alignments.append(aln)
            
            else:
                collapsed_aln = collapse_alignments(alns=current_alignments, aln=count)
                sys.stdout.write(collapsed_aln.to_string() + '\n')
                
                current_seq = f'{aln.query_sequence}_{aln.reference_start}'
                current_alignments = [aln]
                count += 1

        collapsed_aln = collapse_alignments(alns=current_alignments, aln=count)
        sys.stdout.write(collapsed_aln.to_string() + '\n')
        
if __name__ == "__main__":
    args = parse_arguments().parse_args()
    main(sam_file=args.infile)
