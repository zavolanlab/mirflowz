#!/usr/bin/env python

# pylint: disable=line-too-long

"""Tabulate 'bedtools intersect -wo -s' output file.

For each intersecting feature from an INTERSECT file, calculate the sum of its
contributions and print the result in tab-delimited format.

The contribution computation is based on the '--collapsed' and '--nh' flags:
    * If no flag is set, the contribution is 1.
    * If only '--nh' is set, the contribution is 1/NH.
    * If only '--collapsed' is set, the contribution is # of reads/1.
    * If '--collapsed' and '--nh' are set, the contribution of each alignment
        is computed as # of reads/NH.
The values for the amount of collapsed reads (#) and the NH value are inferred
from the sequence name which must follow the format 'read-#_NH'.

EXPECTED INPUT FILE
The expected INTERSECT file must be the output of the call 'bedtools intersect
-a GFF33/GTF -b BAM -wo -s'. To ensure the file follows the expected format,
the first 10 lines are going to be validated.
If the INTERSECT file is empty, no output is produced.

OUTPUT TABLE FORMAT
The basic output table has nor header and two columns:
    * The first one contains the intersecting feature name (determined by the
        value of the '--id' CLI argument, which must match one of the fields
        in the attributes column of the original GFF3/GTF file used in the
        bedtools intersect command.
    * The feature contribution.

Three additional columns can be added when using the flags '--read-ids' and/or
'--feat-extension':
    * '--read-ids' always adds as the last column a semicolon-separated list
        of all the alignment IDs that overlap with that feature.
    * '--feat-extension' adds two additional columns holding the 5' and 3' end
        shifts of the feature if, and only if '--id' is set to "name", and
        "name" contains these shifts separated by an underscore.

Examples
--------
Example 1: Contribution when using '--collapsed'
    use case:
        A single feature with several intersecting reads.
        The flag '--collapsed' is used, so contribution equals the # of reads
        per alignment.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2_1	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1_1	255	+	21

    alignments:
        Read ID: 8-2_1
        Number of collapsed reads: 2
        Contribution: 2

        Read ID: 24-1_1
        Number of collapsed reads: 1
        Contribution: 1

    OUT table:
        hsa-mir-524_-0_+0      3


Example 2: Contribution when using '--nh'
    use case:
        A single feature with several intersecting reads.
        The flag '--nh' is used, so contribution equals 1/NH.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2_1	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1_1	255	+	21

    alignments:
        Read ID: 8-2_1
        Number of mapped genomic loci: 1
        Contribution: 1

        Read ID: 24-1_1
        Number of mapped genomic loci: 1
        Contribution: 1

    OUT table:
        hsa-mir-524_-0_+0      2


Example 3: Contribution when using '--collapsed' and '--nh'
    use case:
        A single feature with several intersecting reads.
        The flags '--nh' and '--contribution' is used, so contribution equals
        # of reads/NH.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2_1	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1_1	255	+	21

    alignments:
        Read ID: 8-2_1
        Number of collapsed reads: 2
        Number of mapped genomic loci: 1
        Contribution: 2/1 = 2

        Read ID: 23-1_1
        Number of collapsed reads: 1
        Number of mapped genomic loci: 1
        Contribution: 1/1 = 1

    OUT table:
        hsa-mir-524_-0_+0      3

Example 4: Column with intersecting reads; using '--read-ids'
    use case:
        A single feature with several intersecting reads.
        Each read contributes with 1.
        Read IDs intersecting the feature are appended as the last column.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2_1	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1_1	255	+	21

    alignments:
        Read ID: 8-2_1
        Contribution: 1

        Read ID: 24-1_1
        Contribution: 1

    OUT table:
        hsa-mir-524_-0_+0      2       8-2_1;24-1_1

Example 5: Columns with feature shifts; using '--feat-extension'
    use case:
        A single feature with several intersecting reads.
        Each read contributes with 1.
        Feature start and end coordinates shift are appended as the third and
        fourth column respectively.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-1_+3	19	44413	44434	8-2_1	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-1_+3	19	44413	44434	24-1_1	255	+	21

    alignments:
        Read ID: 8-2_1
        Contribution: 1

        Read ID: 24-1_1
        Contribution: 1

    feature:
        Feature name: hsa-mir-524_-1_+3
        5' shift: -1
        3' shift: +3

    OUT table:
        hsa-mir-524_-1_+3      2       -1       +3
"""  # noqa: E501
# pylint: enable=line-too-long

import argparse
from pathlib import Path
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .validate_bedtools_intersect import (
        FileFormatError,
        parse_all,
        validate_first_n,
    )
else:
    try:
        from .validate_bedtools_intersect import (
            FileFormatError,
            parse_all,
            validate_first_n,
        )
    except ImportError:  # pragma: no cover
        from validate_bedtools_intersect import (
            FileFormatError,
            parse_all,
            validate_first_n,
        )


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.1.0",
        help="Show program's version number and exit",
    )
    parser.add_argument(
        "intersect",
        help=(
            "Path to the INTERSECT file. This file must be the output of"
            " a bedtools intersect call with -a being a BED file and"
            " -b a BAM file."
        ),
        type=Path,
    )
    parser.add_argument(
        "--collapsed",
        help=(
            "Indicate that the file used in bedtools intersect has the"
            " reads collapsed by sequence and alignment. The collapsed name"
            " must be build by the alignment name followed by a '-' and the"
            " number of collpased alignments, i.e 1-4. Default %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--nh",
        help=(
            "Indicate that the file used in bedtools intersect has the"
            " NH tag in the read query name. The name must be build by the"
            " alignment name followed by an underscore and the NH value,"
            " i.e 1-2_4. Default %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--id",
        help=(
            "ID used to identify the feature in the output table."
            " The ID must be in lowercase. Default: %(default)s."
        ),
        default="name",
        type=str,
    )
    parser.add_argument(
        "--read-ids",
        help=(
            "Include read IDs of the alignments intersecting a feature in"
            " the output table. Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--feat-extension",
        help=(
            "If any of the feature's coordinates had been extended, include"
            " the extension in the output table. It is assumed that the"
            " extensions are found within the feature id 'name' and separated"
            " by an underscore. Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )

    return parser


def get_contribution(
    query_id: str, collapsed: bool = False, nh: bool = False
) -> float:
    """Get contribution of an alignment to the overall count."""
    if collapsed and nh:
        num_reads = int(query_id.split("-")[1].split("_")[0])
        nh_value = int(query_id.split("-")[1].split("_")[1])

    elif not collapsed and nh:
        num_reads = 1
        nh_value = int(query_id.split("_")[1])

    elif collapsed and not nh:
        num_reads = int(query_id.split("-")[1])
        nh_value = 1

    else:
        num_reads = 1
        nh_value = 1

    return num_reads / nh_value


def get_initial_data(name: str, feat_extension: bool) -> list[str]:
    """Get the feature name and its extension.

    Args:
        name:
            string with the feature name that can or not include its
            annotation extension in the format name_5-extension_3-extension
        feat_extension:
            specify whether the feature annotation extension has to be a field
            on the final output

    Returns:
        list with the feature name to be found in the final table and the
        number of extended positions (if asked for)
    """
    if feat_extension:
        feat_data = name.split("_")

        if len(feat_data) == 1:
            feat_data.extend(["NA", "NA"])
    else:
        feat_data = [name]

    return feat_data


def main(args) -> None:
    """Tabulate a bedtools intersect file."""
    try:
        validate_first_n(intersect_file=args.intersect, n=10)
    except FileFormatError as err:
        raise err

    count = 0.0
    current_name = None
    read_id = []

    for _, rec in parse_all(intersect_file=args.intersect):

        name = rec.feat_attrs[args.id]
        contribution = get_contribution(rec.read_name, args.collapsed, args.nh)

        if current_name is None:
            current_name = name
            feat_data = get_initial_data(name, args.feat_extension)

        if current_name == name:
            count += contribution
            read_id.append(rec.read_name)

        else:
            feat_data.insert(1, str(count))

            if args.read_ids:
                feat_data.append(";".join(sorted(read_id)))

            sys.stdout.write("\t".join(feat_data) + "\n")

            feat_data = get_initial_data(name, args.feat_extension)

            current_name = name
            count = contribution
            read_id = [rec.read_name]

    if current_name is not None:
        feat_data.insert(1, str(count))

        if args.read_ids:
            feat_data.append(";".join(sorted(read_id)))

        sys.stdout.write("\t".join(feat_data) + "\n")


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
