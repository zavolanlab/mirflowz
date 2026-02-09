#!/usr/bin/env python

# pylint: disable=line-too-long

"""Tabulate 'bedtools intersect -wo -s' output file.

For each intersecting feature from an INTERSECT file, calculate the sum of its
contributions and print the result in tab-delimited format.

Contribution logic (controlled by `--collapsed` and `--nh`):
    * No flags: contribution = 1
    * Only `--nh`: contribution = 1 / NH
    * Only `--collapsed`: contribution = #reads / 1
    * Both `--collapsed` and `--nh`: contribution = #reads / NH

The number of collapsed reads (#reads) and the NH value are inferred from the
read name, which must follow one of these formats:
    * plain: READ
    * NH only: READ_NH
    * collapsed only: READ-#reads
    * collapsed + NH: READ-#reads_NH

EXPECTED INPUT FILE
The expected INTERSECT file must be the output of:
    bedtools intersect -wo -s -a GFF3/GTF -b BAM
The first 10 lines are validated for format. If the INTERSECT file is empty,
no output is produced.

OUTPUT TABLE FORMAT
The basic output table has no header and two columns:
    1) Feature identifier (taken from the attribute named by `--id`,
       which must match a key in the attributes column of the original
       GFF3/GTF features used in the bedtools call)
    2) Feature contribution (float)

Optional columns (controlled by flags):
    * `--read-ids` appends a semicolon-separated list of read IDs that
      overlap the feature (always added as the last column).
    * `--feat-extension` appends two columns with 5′ and 3′ extension sizes.
      These are inferred by splitting the FEATURE-ID value (selected via `--id`)
      on underscores, assuming the convention: NAME_5EXT_3EXT. If that pattern
      isn't present in the selected attribute, the extension columns will be 0.


Examples
--------
Example 1: Contribution when using '--collapsed'
    use case:
        A single feature with several intersecting reads.
        The flag '--collapsed' is used, so contribution equals the # of reads
        per alignment.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1	255	+	21

    alignments:
        Read ID: 8-2
        Number of collapsed reads: 2
        Contribution: 2

        Read ID: 24-1
        Number of collapsed reads: 1
        Contribution: 1

    OUT table:
        hsa-mir-524_-0_+0      3


Example 2: Contribution when using '--nh'
    use case:
        A single feature with several intersecting reads.
        The flag '--nh' is used, so contribution equals 1/NH.

    IN INTERSECT records:
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8_1	255	+	21
        19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24_1	255	+	21

    alignments:
        Read ID: 8_1
        Number of mapped genomic loci: 1
        Contribution: 1

        Read ID: 24_1
        Number of mapped genomic loci: 1
        Contribution: 1

    OUT table:
        hsa-mir-524_-0_+0      2


Example 3: Contribution when using '--collapsed' and '--nh'
    use case:
        A single feature with several intersecting reads.
        The flags '--nh' and '--contribution' are used, so contribution equals
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
            " the call 'bedtools intersect -wo -s -a GFF3/GTF -b BAM'."
        ),
        type=Path,
    )
    parser.add_argument(
        "--collapsed",
        help=(
            "Indicate that reads were collapsed by sequence/alignment before "
            "running bedtools. Read names must be built as READ-#reads (or "
            "READ-#reads_NH when using --nh). Example: abc123-4 or abc123-4_2."
            " Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--nh",
        help=(
            "Indicate that the NH tag is encoded in the read name. Read names "
            "must be built as READ_NH (or READ-#reads_NH when using "
            "--collapsed). Example: abc13_2 or abc13-4_2. "
            "Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--id",
        help=(
            "Attribute key (lowercase) used to identify the feature in the "
            "ouput table. This key must exist in the attributes column of the "
            "GFF3/GTF features. Default: %(default)s."
        ),
        default="name",
        type=str,
    )
    parser.add_argument(
        "--read-ids",
        help=(
            "Append a semicolon-separated list of intersecting read IDs as the"
            " last output column. Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--feat-extension",
        help=(
            "Append two columns with 5′ and 3′ extension sizes of the feature."
            " Assumes the selected attribute (`--id`) encodes NAME_5EXT_3EXT. "
            "If not present, zeros are reported. Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )

    return parser


def get_contribution(
    query_id: str, collapsed: bool = False, nh: bool = False
) -> float:
    """Return the contribution of a single alignment.

    The read name format depends on the flags:
        * collapsed + nh: READ-#reads_NH
        * nh only:        READ_NH
        * collapsed only: READ-#reads
        * neither:        READ

    Args:
        query_id: Read/query name from the INTERSECT record.
        collapsed: If True, parse #reads from the read name.
        nh: If True, parse NH from the read name.

    Returns:
        Contribution as (#reads / NH) following the rules above.
    """
    if collapsed and nh:
        try:
            num_reads = int(query_id.split("-")[1].split("_")[0])
            nh_value = int(query_id.split("-")[1].split("_")[1])

        except IndexError as err:
            raise IndexError(
                f'Malformed read name (no delimiters): "{query_id}".\n'
                "The flags '--collapsed' and '--nh' had been used and the "
                "expected read format is READ-#reads_NH."
            ) from err

        except ValueError as err:
            raise ValueError(
                f'Malformed read name (no int values): "{query_id}".\n'
                "The flags '--collapsed' and '--nh' had been used and the "
                "expected read format is READ-#reads_NH (str-int_int)."
            ) from err

    elif not collapsed and nh:

        num_reads = 1

        try:
            nh_value = int(query_id.split("_")[1])

        except IndexError as err:
            raise IndexError(
                f'Malformed read name (no delimiters): "{query_id}".\n'
                "The flag '--nh' has been used and the expected read format "
                "is READ_NH."
            ) from err

        except ValueError as err:
            raise ValueError(
                f'Malformed read name (no int value): "{query_id}".\n'
                "The flag '--nh' has been used and the expected read format "
                "is READ_NH (str_int)."
            ) from err

    elif collapsed and not nh:
        try:
            num_reads = int(query_id.split("-")[1])

        except IndexError as err:
            raise IndexError(
                f'Malformed read name (no delimiters): "{query_id}".\n'
                "The flag '--collapsed' has been used and the expected read "
                "format is READ-#reads."
            ) from err

        except ValueError as err:
            raise ValueError(
                f'Malformed read name: "{query_id}".\n'
                "The flag '--collapsed' has been used and the expected read "
                "format is READ-#reads (str-int)."
            ) from err

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
    """Tabulate 'bedtools intersect -wo -s' output file."""
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
