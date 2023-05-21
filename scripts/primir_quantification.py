#!/usr/bin/env python

"""Tabulate bedtools intersect output BED file.

Read the input BED file, calculate the sum of the intersecting contributions
for each feature and print the result in tab-delimited format. The name of the
feature is determined by the value of the --id argument, which must match one
of the fields in the attributes column of the original GFF/GTF file used in
the bedtools intersect command. The appropriate contribution is based on the
--collapsed and the --nh flags. If --collapsed and --nh are set, the
contribution of each alignment is computed as # of reads/NH. If only 
--collapsed is set, the contribution is # of reads/1. If only --nh is set,
the contribution is 1/NH. Otherwise, the contribution is 1. If the BED file is
empty, no output is produced. The output columns are also determined by the
flags --read-ids and --feat-extension. If --read-ids is set, a
semicolon-separated list of all alignment IDs that overlap with the feature
will always be in the last column. If --feat-extension is set, two additional
columns will be added to the output, containing the 5' and 3' end shifts of the
feature (if found in the feature name, separated by an underscore). If
--feat-extension is set but --id is set to something other than "name", no
extra columns will be added.
"""

import argparse
from collections import namedtuple
from pathlib import Path
import sys
from typing import Dict, Optional


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
        help=(
            "Path to the BED file. This file must be the output of "
            "a bedtools intersect call with -a being a BED file and"
            "-b a BAM file."
        ),
        type=Path
    )
    parser.add_argument(
        '--collapsed',
        help=(
            "Indicate that the file used in bedtools intersect has the"
            "reads collapsed by sequence and alignment. The collapsed name"
            "must be build by the alignment name followed by a '-' and the"
            "number of collpased alignments, i.e 1-4. Default %(default)s."
        ),
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--nh',
        help=(
            "Indicate that the file used in bedtools intersect has the"
            "NH tag in the read query name. The name must be build by the"
            "alignment name followed by an underscore and the NH value,"
            "i.e 1-2_4. Default %(default)s."
        ),
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--id',
        help=(
            "ID used to identify the feature in the output table."
            "The ID must be in lowercase. Default: %(default)s."
        ),
        default="name",
        type=str
    )
    parser.add_argument(
        '--read-ids',
        help=(
            "Include read IDs of the alignments intersecting a feature in"
            "the output table. Default: %(default)s."
        ),
        action='store_true',
        default= False
    )
    parser.add_argument(
        '--feat-extension',
        help=(
            "If any of the feature's coordinates had been extended, include"
            "the extension in the output table. It is assumed that the"
            "extensions are found within the feature id 'name' and separated"
            "by an underscore. Default: %(default)s."
        ),
        action='store_true',
        default=False
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

def get_contribution(query_id:str, collapsed: bool = False, nh: bool = False) -> float:
    """Get contribution of an alignment to the overall count."""
    if collapsed and nh:
        num_reads = int(query_id.split('-')[1].split('_')[0])
        nh_value = int(query_id.split('-')[1].split('_')[1])

    elif not collapsed and nh :
        num_reads = 1
        nh_value = int(query_id.split('_')[1])

    elif collapsed and not nh:
        num_reads = int(query_id.split('-')[1])
        nh_value = 1

    else:
        num_reads = 1
        nh_value = 1
    
    return num_reads/nh_value

def main(args) -> None:
    """Tabulate a bedtools intersect BED file."""

    with open(args.bedfile, 'r') as bedfile:

        Fields = namedtuple('Fields',("feat_chr", "source", "feat_type",
                                      "feat_start", "feat_end", "feat_score", 
                                      "strand", "phase", "feat_attributes", 
                                      "read_chr", "read_start", "read_end", 
                                      "read_name", "read_score", "read_strand"))
        count = 0
        current_name = None
        read_ID = []

        for line in bedfile:

            fields = Fields(*line.strip().split('\t'))

            name = attributes_dictionary(fields.feat_attributes)[args.id]
            contribution = get_contribution(fields.read_name, args.collapsed, args.nh)
            

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
                read_ID.append(fields.read_name)

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
                read_ID = [fields.read_name]
 
        if current_name is not None:
            feat_data.insert(1, str(count))

            if args.read_ids:
                feat_data.append(';'.join(sorted(read_ID)))

            sys.stdout.write('\t'.join(feat_data) + '\n')
        else:
            return

if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
