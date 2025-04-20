#!/usr/bin/env python

# pylint: disable=line-too-long

"""Filter miRNA reads mapped to multiple locations by indel count.

From all the alignments with the same query name and edit distance, keep the
ones with the higher amount of indels. Supplementary alignments are dismissed
as they are part of a chimeric alignment, which is composed of multiple linear
alignments with minimal overlap.

Additionally, the 'NH' and 'HI' tags are updated to match the new amount of
alignments and keep their identifier within the new set respectively. If the
CLI flag `--nh` is set, query names are appended the suffix '_#' were '#' is
the updated alignment's NH tag.

The following assumptions are made:
    - The input SAM file is sorted by query name.


Examples
---------
Example 1: Different number of InDels 
    IN SAM records:
        read-1	0	19	77595	255	8M1D14M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:3G1T2^A14	NH:i:2	NM:i:3	XA:Z:Q	XI:i:1
        read-1	0	19	330456	255	4M1D1M1I3M1D13M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:4^G4^A13	NH:i:2	NM:i:3	XA:Z:Q	XI:i:0

    Alignments:
        CTGACATC-AGTGATTCTCCTGC
        ||| | || |||||||||||||| (1 InDel, 2 mismatches, discarded)
        CTGGCTTCAAGTGATTCTCCTGC

        CTGA-CATCA-GTGATTCTCCTGC
        |||| | ||| ||||||||||||| (3 InDels, 0 mismatches, retained)
        CTGAGC-TCAAGTGATTCTCCTGC

    Command:
        filter_multimappers.py SAM > out_SAM

    OUT SAM record:
        read-1	0	19	330456	255	4M1D1M1I3M1D13M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:4^G4^A13	NH:i:1	HI:i:1  NM:i:3	XA:Z:Q	XI:i:0

Example 2: Equal number of InDels
    IN SAM records:
        read-2	0	19	142777	255	5M1I15M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	NM:i:3	XA:Z:Q	XI:i:0
        read-2	0	19	270081	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14G0G4	NH:i:3	NM:i:3	XA:Z:Q	XI:i:2
        read-2	0	19	545543	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	NM:i:3	XA:Z:Q	XI:i:1

    Alignments:
        GCTTCAAGCCTCCCACCTAGC
        ||||| |||||||||  |||| (1 Indel, 2 mismatches, retained)
        GCTTC-AGCCTCCCAAGTAGC

        GCTTCAAGCCTCCCACCTAGC
        |||||| ||||||||  |||| (1 Indel, 2 mismatches, retained)
        GCTTCA-GCCTCCCAGGTAGC

        GCTTCAAGCCTCCCACCTAGC
        |||||| ||||||||  |||| (1 Indel, 2 mismatches, retained)
        GCTTCA-GCCTCCCAAGTAGC

    Command:
        filter_multimappers.py SAM --nh > out_SAM

    OUT SAM record:
        read-2_3	0	19	142777	255	5M1D15M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	HI:i:1  NM:i:3	XA:Z:Q	XI:i:0
        read-2_3	0	19	270081	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14G0G4	NH:i:3	HI:i:2  NM:i:3	XA:Z:Q	XI:i:2
        read-2_3	0	19	545543	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	HI:i:3  NM:i:3	XA:Z:Q	XI:i:1
"""  # noqa: E501
# pylint: enable=line-too-long

import argparse
from pathlib import Path
import sys
from typing import List

import pysam


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.1.0",
        help="Show program's version number and exit",
    )

    parser.add_argument(
        "infile",
        help="Path to the SAM input file, sorted by query name.",
        type=Path,
    )

    parser.add_argument(
        "--nh",
        help=(
            "If set, the NH tag will be included in the alignment name after"
            " an underscore. Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )

    return parser


def count_indels(aln: pysam.libcalignedsegment.AlignedSegment) -> int:
    """Count the number of indels in an alignment based on its CIGAR string.

    This function counts the number of indels in the alignment based on the
    alignment CIGAR string returned by .cigartuples. Insertions are encoded as
    1s whereas the deletions are encoded as 2.
    Refer to the .cigartuples documentation for more information at:
    https://pysam.readthedocs.io/en/v0.15.0/api.html#pysam.AlignedSegment.cigartuples

    Args:
        aln: alignment to count insertions and deletions for.

    Returns:
        int: sum of insertions and deletions in that alignment.
    """
    cigar = aln.cigartuples
    assert isinstance(cigar, list)

    return sum(op[1] for op in cigar if op[0] == 1 or op[0] == 2)


def find_best_alignments(
    alns: List[pysam.AlignedSegment], nh: bool = False
) -> List[pysam.AlignedSegment]:
    """Find alignments with less indels.

    Using the function `count_indels` each alignment is assigned with its
    number of indels. Then, only those alignments with the higher amount of
    indels are returned. In addition, the 'NH' and 'HI' tags are updated to
    match the new amount of alignments and keep its identifier within the set
    respectively. If `nh` is set to `True`, all query names are appended a
    suffix `_#` were '`#` is the new alignment's NH tag.

    Args:
        alignments: alignments with the same query name

    Returns:
        best_alignments: alignments with the more indels
    """
    if len(alns) == 1:
        if nh:
            name = f'{alns[0].query_name}_{alns[0].get_tag("NH")}'
            alns[0].query_name = name

        return alns

    aln_indels = [(aln, count_indels(aln=aln)) for aln in alns]
    max_indels = max(aln_indels, key=lambda x: x[1])[1]
    best_alignments = [
        aln
        for i, (aln, indels) in enumerate(aln_indels)
        if indels == max_indels
    ]

    for i, best_aln in enumerate(best_alignments):

        if nh:
            name = f"{best_aln.query_name}_{len(best_alignments)}"
            best_aln.query_name = name

        best_aln.set_tag("NH", len(best_alignments))
        best_aln.set_tag("HI", i + 1)

    return best_alignments


def write_output(alns: List[pysam.AlignedSegment]) -> None:
    """Write the output to the standard output (stdout).

    Args:
        alignments: alignments with the same query name
    """
    for alignment in alns:
        sys.stdout.write(alignment.to_string() + "\n")


def main(arguments) -> None:
    """Filter multimappers by indels count."""
    with pysam.AlignmentFile(arguments.infile, "r") as samfile:

        sys.stdout.write(str(samfile.header))

        current_alignments: list[pysam.AlignedSegment] = []
        current_query = None

        for alignment in samfile:

            if alignment.is_supplementary:
                continue

            if current_query is None:
                current_query = alignment.query_name

            if current_query == alignment.query_name:
                current_alignments.append(alignment)

            else:
                current_alignments = find_best_alignments(
                    current_alignments, arguments.nh
                )
                write_output(alns=current_alignments)

                current_query = alignment.query_name
                current_alignments = [alignment]

        if len(current_alignments) > 0:
            current_alignments = find_best_alignments(
                current_alignments, arguments.nh
            )
            write_output(alns=current_alignments)


if __name__ == "__main__":
    args = parse_arguments().parse_args()  # pragma:no cover
    main(args)  # pragma: no cover
