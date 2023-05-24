#!/usr/bin/env python

"""Add intersecting feature(s) into a SAM file as a tag.

Add intersecting feature names from a BED file as a tag to alignments in a SAM
file using the format: feat-name|5p-shift|3p-shift|CIGAR|MD. Multiple
intersecting feature names are separated by a semi-colon. The 5p-shift and the
3p-shift are calculated as the difference between the alignment and the feature
start/end positions. If --extension is provided, it adjusts the feature start
position by adding the given value and subtracts it from the end position. The
--id string specifies the feature name to be used within the attributes column
in the BED file. If either the BED or the SAM file is empty, only the SAM file
header will be returned.
"""

import argparse
from collections import defaultdict, namedtuple
from pathlib import Path
import sys
from typing import Dict, Optional

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
        '-b', '--bed',
        help=(
            "Path to the BED file. This file must be the output of "
            "a bedtools intersect call with -a being a GFF3 file and"
            "-b a BAM file."
        ),
        type=Path,
        required=True
    )
    parser.add_argument(
        '-s', '--sam',
        help="Path to the SAM input file.",
        type=Path,
        required=True
    )
    parser.add_argument(
        '-e', '--extension',
        help=(
            "Number of nucleotides the start and end coordinates of the"
            "annotated feature has been extended. Default: %(default)d."
        ),
        default=0,
        type=int
    )
    parser.add_argument(
        '--id',
        help=(
            "ID used to identify the feature in the name add as tag."
            "The ID must be in lowercase. Default: %(default)s."
        ),
        default="name",
        type=str
    )

    return parser


def attributes_dictionary(attr: str) -> Optional[Dict[str, str]]:
    """Create attributes dicctionary."""
    pairs = attr.split(';')

    if len(pairs[0].split('=')) == 2:
        attr_dict = {p.split('=')[0].lower(): p.split('=')[1] for p in pairs}
    else:
        attr_dict = {p.split('"')[0].strip().lower(): p.split('"')[1] for p in pairs}

    return attr_dict


def parse_intersect_output(intersect_file: Path, id: str = "name", extension: int = 0) -> defaultdict(list):
    """Parse intersect BED file.

    Given a BED file generated by intersecting a GFF file (-a) with a BAM file
    (-b) using bedtools intersect, create a dictionary where the alignment
    names are the keys. The values are lists containing the feature name,
    start position, and end position. The id argument specifies the feature
    name to use, and the extension argument adjusts the feature coordinates by
    adding the given value and subtracts it from the end position. If the BED
    file is empty, `None` is returned.

    Args:
        intersect_file:
            Path to the intersect BED file.
        id:
            ID used to identify the feature. Defaults to "name".
        extension:
            Number of nucleotides the start and end coordinates have to be
            adjusted. Defaults to 0.
    """
    intersect_data = defaultdict(list)
    Fields = namedtuple('Fields', ("feat_chr", "source", "feat_type",
                                   "feat_start", "feat_end", "feat_score",
                                   "strand", "phase", "feat_attributes",
                                   "read_chr", "read_start", "read_end",
                                   "read_name", "read_score", "read_strand",
                                   "overlap_len"))

    with open(intersect_file, 'r') as bedfile:
        for line in bedfile:

            fields = Fields(*line.strip().split('\t'))

            miRNA_name = attributes_dictionary(fields.feat_attributes)[id]
            miRNA_start = int(fields.feat_start) + extension
            miRNA_end = int(fields.feat_end) - extension

            intersect_data[fields.read_name].append((miRNA_name,
                                                     miRNA_start,
                                                     miRNA_end))

    if not intersect_data:
        return None
    else:
        return intersect_data


def get_tags(intersecting_mirna: list, alignment: pysam.AlignedSegment, extension: int = 0) -> set:
    """Get tag for alignment.

    Given an alignment and a list containing the feature name, start position,
    and end position, create a list of strings to be added as a new tag to that
    alignment. The string has the format: feat-name|5p-shift|3p-shift|CIGAR|MD.
    The 5p-shift and 3p-shift are calculated as a difference between the
    feature start/end position and the alignment start/end position. If the
    start and end position of the alignment differs at most by the extension
    argument value to the feature start and end positions respectively, the
    name will be add to the final list.

    Args:
        intersecting_mirna:
            list with the feature name, start and end positions
        alignment:
            alignment to create the tag for
        extension:
            maximum number of nucleotides the alignment start and end positions
            can differ from the feature to count it as an intersecting feature

    Returns:
        tags:
            set of strings containing the new tag
    """
    cigar = alignment.cigarstring
    md = alignment.get_tag('MD')

    tags = []

    for miRNA_name, miRNA_start, miRNA_end in intersecting_mirna:
        shift_5p = alignment.reference_start - miRNA_start + 1
        shift_3p = alignment.reference_end - miRNA_end
        limit = extension + 1

        if -limit < shift_5p < limit and -limit < shift_3p < limit:
            tags.append(f'{miRNA_name}|{shift_5p}|{shift_3p}|{cigar}|{md}')

    return set(tags)


def main(args) -> None:
    """Add intersecting feature(s) into a SAM file as a tag."""
    intersect_data = parse_intersect_output(args.bed, args.id, args.extension)

    with pysam.AlignmentFile(args.sam, 'r') as samfile:

        sys.stdout.write(str(samfile.header))

        if intersect_data is None:
            return

        for alignment in samfile:
            alignment_id = alignment.query_name
            intersecting_miRNAs = intersect_data.get(alignment_id, [])

            tags = get_tags(intersecting_mirna=intersecting_miRNAs,
                            alignment=alignment,
                            extension=args.extension)

            alignment.set_tag('YW', ';'.join(tags))
            sys.stdout.write(alignment.to_string() + '\n')


if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover