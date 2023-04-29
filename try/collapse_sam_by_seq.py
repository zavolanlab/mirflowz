#!/usr/bin/env python

"""Collapse alignments in a SAM file by query sequence.

This script uses the pysam library to open the input SAM file and iterates
over each alignment in it to create a dictionary which has the sequence as key
and the different alignments as value. Thereafter, if and entry on the 
dictionary contains more than one alignment, it will be collapsed into a single
alignment. The alignment collapsing implies the merging of the different
queries IDs into a single one, and the update of the "HI" tag. In addition,
all the entries in the dictionary will get a new tag, "XC", containing the 
division between to amount of collapsed alignments in that entry and the 
"NH" tag. 

Input:
    input_sam_file: Path to the input SAM file, ordered by query name.

Output:
    standard output: Collapsed SAM file

Functions:
    parse_arguments():
        Command-line arguments parser
    collapse_alignments():
        Collapse alignments by query sequence.
    write_output():
        Write the output to the standard output (stdout).
    main():
        Collapse alignments by query sequence.

Usage:
    collapse_sam_by_seq.py SAM
    collapse_sam_by_seq.py SAM > SAM
"""

import argparse
from pathlib import Path
import sys
from typing import List

import pysam

def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Script to collapse alignments by query sequence."
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

def collapse_alignments(alignments: List[pysam.AlignedSegment]) -> pysam.AlignedSegment:
    """Collapse alignments by query sequence.

    This function collapse the different alignments with the same query
    sequence into a single. For this purpose the query ID of the new alignment
    will consist on the all the IDs merged into one, and the new tag "XC"
    will be created to store the total count. This new value is defined as the
    division of the number of collapsed alignments by the "NH" tag.
    
    Args:
        alignments: alignments with the same query sequence

    Retrns:
        new_alignment: collapsed alignment
    """
    new_alignment = alignments[0]
            
    read_ids = [a.query_name for a in alignments]
    new_id = "_".join(read_ids)
    new_alignment.query_name = new_id
            
    xc_count = len(read_ids) / new_alignment.get_tag("NH")
    new_alignment.set_tag("XC", xc_count)
    new_alignment.set_tag("HI", 1)

    return new_alignment

def write_output(alignments: List[pysam.AlignedSegment]) -> None:
    """Write the output to the standard output (stdout).

    Args:
        alignments: alignments with the same sequence
    """
    if len(alignments) == 1:
        xc_count = 1 / alignments[0].get_tag("NH")
        alignments[0].set_tag("XC", xc_count)
        sys.stdout.write(alignments[0].to_string() + '\n')

    else:
        new_alignment = collapse_alignments(alignments)
        sys.stdout.write(new_alignment.to_string() + '\n')

def main(sam_file: Path) -> None:
    """Collapse alignments by query sequence.

    Args:
        sam_file: Path to the input SAM file.
    """
    
    with pysam.AlignmentFile(sam_file, "r") as samfile:

        sys.stdout.write(str(samfile.header))
    
        alignments = {}
    
        for aln in samfile:

            seq = aln.query_sequence
            
            if seq in alignments:
                alignments[seq].append(aln)
            else:
                alignments[seq] = [aln]
    
        for seq, align_list in alignments.items():
            write_output(alignments=align_list)

if __name__ == "__main__":
    args = parse_arguments().parse_args()
    main(sam_file=args.infile)
