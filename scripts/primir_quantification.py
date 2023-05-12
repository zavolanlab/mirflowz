#!/usr/bin/env python

"""Tabulate bedtools intersect output BED file.

Read the input BED file, calculate the sum of the intersecting contributions
for each feature and print the result in tab-delimited format. The name of the
feature is determined by the value of the --id argument, which must match one
of the fields in the attributes column of the original GFF/GTF file used in
the bedtools intersect command. If a SAM file is provided, the contribution
of each alignment is computed as 1/NH.If an alignment is not found, its
contribution is assumed to be 0. If no SAM file is provided, the contribution
of each alignment is assumed to be 1. If either the BED or SAM file is empty,
no output is produced. The output columns are also determined by the flags --read-ids
and --feat-extension. If --read-ids is set, a semicolon-separated list of all
alignment IDs that overlap with the feature will always be in the last column.
If --feat-extension is set, two additional columns will be added to the output,
containing the 5' and 3' end shifts of the feature (if found in the feature
name, separated by an underscore). If --feat-extension is set but --id is set
to something other than "name", no extra columns will be added.
"""

import argparse
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
        'bedfile',
        help="Path to the BED file. This file must be the output of a \
            bedtools intersect call.",
        type=Path
    )
    parser.add_argument(
        '-s', '--sam',
        help="Path to the SAM file containing the intersecting alignments.",
        type=Path
    )
    parser.add_argument(
        '--id',
        help="ID used to identify the feature in the output table. \
            The ID must be in lowercase. Default: %(default)s.",
        default="name",
        type=str
    )
    parser.add_argument(
        '--read-ids',
        help="Include read IDs of the alignments intersecting a feature in \
            the output table. Default: %(default)s.",
        action='store_true',
        default= False
    )
    parser.add_argument(
        '--feat-extension',
        help="If any of the feature's coordinates had been extended, include \
            the extension in the output table. It is assumed that the \
            extensions are found within the feature id 'name' and separated \
            by an underscore. Default: %(default)s.",
        action='store_true',
        default=False
    )

    return parser

def nh_dictionary(samfile: Path) -> Optional[Dict[str, int]]:
    """Create dictionary from SAM file.
    
    Create a dictionary form a SAM file where keys are the alignments IDs
    and values are the NH tag.
    """
    with pysam.AlignmentFile(samfile, 'r') as samfile:
        align_nh =  {}

        for alignment in samfile:
            align_nh[alignment.query_name] = alignment.get_tag('NH')

    return align_nh

def attributes_dictionary(attr: str) -> Optional[Dict[str, str]]:
    """Create attributes dicctionary."""
    pairs = attr.split(';')

    if len(pairs[0].split('=')) == 2:
        attr_dict = {p.split('=')[0].lower(): p.split('=')[1] for p in pairs}
    else:
        attr_dict = {p.split('"')[0].strip().lower(): p.split('"')[1] for p in pairs}

    return attr_dict

def main(args) -> None:
    """Tabulate a bedtools intersect BED file."""

    with open(args.bedfile, 'r') as bedfile:
        if len(bedfile.read()) == 0:
            return

    align_nh: Optional[Dict] = None 
    if args.sam:
        align_nh = nh_dictionary(args.sam)
    
        if len(align_nh.items()) == 0:
            return
    
    with open(args.bedfile, 'r') as bedfile:

        count = 0
        current_name = None
        read_ID = []

        for line in bedfile:

            line = line.strip().split('\t')
            name = attributes_dictionary(line[9])[args.id]
            
            if align_nh:
                try:
                    contribution = 1/align_nh[line[13]]
                except KeyError:
                    continue
            else:
                contribution = 1

            if current_name is None:
                current_name = name
                if args.feat_extension:
                    feat_data = name.split('_')
                    if len(feat_data) == 1:
                        feat_data.extend(['NA', 'NA'])
                else:
                    feat_data = [name]

            if current_name == name:
                count += contribution
                read_ID.append(line[13])

            else:
                feat_data.insert(1, str(count))

                if args.read_ids:
                    feat_data.append(';'.join(sorted(read_ID)))
                    
                sys.stdout.write('\t'.join(feat_data) + '\n')

                if args.feat_extension:
                    feat_data = name.split('_')
                    if len(feat_data) == 1:
                        feat_data.extend(['NA', 'NA'])
                else:
                    feat_data = [name]

                current_name = name
                count = contribution
                read_ID = [line[13]]
 
        feat_data.insert(1, str(count))

        if args.read_ids:
            feat_data.append(';'.join(sorted(read_ID)))

        sys.stdout.write('\t'.join(feat_data) + '\n')

if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
