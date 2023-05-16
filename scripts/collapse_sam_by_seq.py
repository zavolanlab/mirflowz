#!/usr/bin/env python

"""Collapse reads in a SAM file by alignment and query sequence.

Collapse alignments in a SAM file by query sequence and starting position, and
write the collapsed alignment to the standard output. The process is divided
into two steps: first, alignments with the same sequence and starting position
are grouped into a list. Second, the alignments in each list are combined into
a single alignment. The new alignment is given a name build of the number of
alignments collapsed and the indicated integer, and the "HI" tag is reset to 1. 
Two new tags may be also add: the "XC" tag is the count of collapsed alignments
divided by the "NH" tag. If the "NH" tag is missing, the value 1 is used. If
`--read-ids` is set, the "ZS" tag will be add. This new tag contains a string
with all the alignment IDs that are combined to create the new one separted by
an underscore.
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
    ),
    parser.add_argument(
        '--read-ids',
        help="Include read IDs of the collpased alignments in the new tag ZS. \
            Default: %(default)s.",
        action='store_true',
        default= False
    )

    return parser


def collapse_alignments(alns: list[pysam.AlignedSegment], aln: int, reads: bool=False) -> pysam.AlignedSegment:
    """Collapse an alignments list into a single alignment.

    Combine the alignments in the list into a single alignment. The new
    alignment is given a name build of the number of alignments collapsed and
    the indicated integer, and the "HI" tag is reset to 1. Two new tags may be
    also add: the "XC" tag is the count of collapsed alignments divided by the
    "NH" tag. If the "NH" tag is missing, the value 1 is used. If `reads` is
    set to True, the "ZS" tag will be add. This new tag contains a string with
    all the alignment IDs that are combined to create the new one separted by
    an underscore.

    Args:
        alns: alignments with the same query sequence
        aln: integer to define the new alignment
        reads: flag to do or do not add the new tag

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
        new_alignment.set_tags([("NH", nh_tag),
                                ("MD", new_alignment.get_tag("MD")),
                                ("XC", xc_count)])
        if reads:
            new_alignment.set_tag("ZS", new_alignment.query_name)

    else:
        xc_count = len(alns) / nh_tag
        read_ids = sorted({aln.query_name for aln in alns})
        new_alignment.set_tags([("NH", nh_tag),
                                ("MD", new_alignment.get_tag("MD")),
                                ("XC", xc_count), 
                                ("HI", 1)])
        if reads:
            new_alignment.set_tag("ZS", "_".join(read_ids))
        
    new_alignment.query_name = f'{len(alns)}.{aln}'

    return new_alignment


def main(args) -> None:
    """Collapse reads by alignment and query sequence."""
    with pysam.AlignmentFile(args.infile, 'r') as samfile:

        sys.stdout.write(str(samfile.header))

        try:
            aln = next(samfile)
            current_seq = f'{aln.query_sequence}_{aln.reference_start}_{aln.cigar}_{aln.get_tag("MD")}_{aln.flag & 0x04}'
            current_alignments = [aln]
            count = 1

        except StopIteration:
            return

        for aln in samfile:

            if current_seq == f'{aln.query_sequence}_{aln.reference_start}_{aln.cigar}_{aln.get_tag("MD")}_{aln.flag & 0x04}':
                current_alignments.append(aln)
            
            else:
                collapsed_aln = collapse_alignments(alns=current_alignments, aln=count, reads=args.read_ids)
                sys.stdout.write(collapsed_aln.to_string() + '\n')
                
                current_seq = f'{aln.query_sequence}_{aln.reference_start}_{aln.cigar}_{aln.get_tag("MD")}_{aln.flag & 0x04}'
                current_alignments = [aln]
                count += 1

        collapsed_aln = collapse_alignments(alns=current_alignments, aln=count, reads=args.read_ids)
        sys.stdout.write(collapsed_aln.to_string() + '\n')
        
if __name__ == "__main__":
    args = parse_arguments().parse_args() # pragma: no cover
    main(args) # pragma: no cover
