#!/usr/bin/env python

# pylint: disable=line-too-long

"""Add intersecting feature(s) into a SAM file as a tag.

Build new names for the intersecting features from an INTERSECT file and add
them as a tag to alignments in a SAM file using the format
FEATURE_ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ. If either the INTERSECT or the
SAM file is empty, only the SAM file header is returned.

EXPECTED INPUT FILES
The expected INTERSECT file must be the output of the call 'bedtools intersect
-a GFF3/GTF -b BAM -wo -s'. To ensure the file follows the expected format, the
first 10 lines are going to be validated. If the GFF3/GTF used in the call  has
the feature start and end coordinates extended, the number of additional
nucleotides can be specified using the CLI option `--extension`.
The SAM file must contain only the reads that have an intersecting feature.

NAME CREATION and TAG ADDITION
For each alignment, the name of the intersecting feature will follow the
format FEATURE_ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ. The CLI option `--id`
specifies the feature identifier to be used as FEATURE_ID from within the
attributes column in the INTERSECT file. The 5p-shift and the 3-p shift values
are the difference between the feature start and end coordinates and the
alignment start and end coordinates. If `--extension` is provided, the feature
start position are adjusted by adding the given value and subtracting it from
the end position. If both, the 5p-shift and the 3p-shift, are within the range
+/- extension + 1 the feature name is added to the alignment as the new tag
"YW". Multiple intersecting feature names are separated by a semi-colon.


Examples
--------
Example 1: Feature intersects alignment; coordinates adjustment and shift allowed
    use case:
        Prior to checking if the feature is intersecting the alignment, its
        coordinates are adjusted by the value specified in `--extension`. In
        addition, the same value is used to specify the +/- shift allowed
        between the feature and the read alignment start and end coordinates.

    command:
        iso_name_tagging.py -i INTERSECT -s SAM --extension 5

    in INTERSECT record:
        19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_1	255	+	21

    in SAM record:
        read_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

    intersection before coordinates adjustment visualization:

        ---|===============================|--- (feature)
        ---------|===================|--------- (read)

    intersection after coordinates adjustment visualization:

        --------|=====================|-------- (feature)
        ---------|===================|--------- (read)

    out SAM record:
        read_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:hsa-miR-1323|1|-1|21M|21|TCAAAACTGAGGGGCATTTTC

    description:
        The feature start and end coordinates after the adjustment are 5337
        and 5360 respectively.
        The read alignment starts at position 5338. As the read has length 21,
        its end position is 5359.
        The feature intersects the read alignment with an overhang within the
        specified shift range (+/- 5) so it is added as a new tag in the output
        SAM record.


Example 2: Feature intersects alignment; no coordinates adjustment or shift allowed
    use case:
        The feature and read alignment coordinates must perfectly match.

    command:
        iso_name_tagging.py -i INTERSECT -s SAM

    in INTERSECT record:
        19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_2	255	+	21

    in SAM record:
        read_2	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

    intersection visualization:

        ---------|===================|--------- (feature)
        ---------|===================|--------- (read)

    out SAM record:
        read_2	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

    description:
        The feature start and end coordinates are 5338 and 5359 respectively.
        The read alignment starts at position 5338. As the read has length 21,
        its end position is 5359.
        The feature perfectly intersects the read alignment so it is added as
        a new tag in the output SAM record.


Example 3: Non-intersecting feature; shift filter not passed
    use case:
        Prior to checking if the feature is intersecting the alignment, its
        coordinates are adjusted by the value specified in `--extension`. In
        addition, the same value is used to specify the +/- shift allowed
        between the feature and the read alignment start and end coordinates.

    command:
        iso_name_tagging.py -i INTERSECT -s SAM --extension 1

    in INTERSECT record:
        19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_3	255	+	21

    in SAM record:
        read_3	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

    intersection before coordinates adjustment visualization:

        ---|===============================|--- (feature)
        ---------|===================|--------- (read)

    intersection after coordinates adjustment visualization:

        ----|=============================|---- (feature)
        ---------|===================|--------- (read)

    out SAM record:
        read_3	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:

    description:
        The feature start and end coordinates after the adjustment are 5333
        and 5364 respectively.
        The read alignment starts at position 5338. As the read has length 21,
        its end position is 5359.
        There is a 5-nucleotide overhang on both ends. Thus, the feature is
        not considered to intersect the read alignment and the tag is an empty
        string.


Example 4: Feature intersects alignment; using feature's "Alias"
    use case:
        The feature and read alignment coordinates must perfectly match.

    command:
        iso_name_tagging.py -i INTERSECT -s SAM --id alias

    in INTERSECT record:
        19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_4	255	+	21

    in SAM record:
        read_4	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

    intersection visualization:

        ---------|===================|--------- (feature)
        ---------|===================|--------- (read)

    out SAM record:
        read_4	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:MIMAT0005795|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

    description:
        The feature start and end coordinates are 5338 and 5359 respectively.
        The read alignment starts at position 5338. As the read has length 21,
        its end position is 5359.
        The feature perfectly intersects the read alignment so it is added as
        a new tag in the output SAM record. In this case, instead of using the
        feature `Name` (default), the `Alias` is used.
"""  # noqa: E501
# pylint: enable=line-too-long

import argparse
from collections import defaultdict
from pathlib import Path
import sys
from typing import Dict, Optional

import pysam

from validate_bedtools_intersect import (
        FileFormatError,
        Record,
        validate_first_n,
        parse_all,
)


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
        version="%(prog)s 1.2.0",
        help="Show program's version number and exit",
    )
    parser.add_argument(
        "-i",
        "--intersect",
        help=(
            "Path to the INTERSECT file. This file must be the output of"
            " the call 'bedtools intersect -wo -s -a GFF3/GTF -b BAM'."
        ),
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-s",
        "--sam",
        help="Path to the SAM input file.",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-e",
        "--extension",
        help=(
            "Number of nucleotides the start and end coordinates of the"
            " annotated features had been extended. Default: %(default)d."
        ),
        default=0,
        type=int,
    )
    parser.add_argument(
        "--id",
        help=(
            "ID used to identify the feature in the name that is added as tag."
            " The ID must be in lowercase. Default: %(default)s."
        ),
        default="name",
        type=str,
    )

    return parser


def parse_intersect_output(
    intersect_file: Path, feat_id: str = "name", extension: int = 0
) -> Optional[Dict[Optional[str], list]]:
    """Parse 'bedtools intersect -wo -s' output file.

    Given an INTERSECT file generated by intersecting a GFF3/GTF file (-a) with
    a BAM file (-b) using bedtools intersect, create a dictionary where the
    alignment names are the keys. The values are lists containing the feature
    name, start and end positions. The 'feat_id' argument specifies
    the feature name to use, and the extension argument adjusts the feature
    coordinates by adding the given value and subtracts it from the end
    position. If the INTERSECT file is empty, `None` is returned.

    Args:
        intersect_file:
            Path to the INTERSECT file.
        feat_id:
            ID used to identify the feature. Defaults to "name".
        extension:
            Number of nucleotides the start and end coordinates have to be
            adjusted. Defaults to 0.
    """
    intersect_data = defaultdict(list)

    for _, rec in parse_all(intersect_file=intersect_file):
        miRNA_name = rec.feat_attrs[feat_id]
        miRNA_start = rec.feat_start + extension
        miRNA_end = rec.feat_end - extension

        intersect_data[rec.read_name].append(
            (miRNA_name, miRNA_start, miRNA_end)
        )

    if not intersect_data:
        return None

    return intersect_data


def get_tags(
    intersecting_mirna: list, alignment: pysam.AlignedSegment, extend: int = 0
) -> set:
    """Get tag for alignment.

    Given an alignment and a list containing the feature name, start position,
    and end position, create a list of strings to be added as a new tag to that
    alignment. The string has the format:
        feature-id|5p-shift|3p-shift|CIGAR|MD|READ_SEQ
    The 5p-shift and 3p-shift are calculated as a difference between the
    feature start/end position and the alignment start/end position. If the
    start and end position of the alignment differs at most by the extension
    argument value to the feature start and end positions respectively,
    the name will be add to the final list.

    Args:
        intersecting_mirna:
            list with the miRNA species name, start and end positions
        alignment:
            alignment to create the tag for
        extend:
            value that sets the range in which both shifts have to be in to
            create a new string for a particular miRNA species

    Returns:
        tags: set of strings containing the new tag
    """
    cigar = alignment.cigarstring
    seq = alignment.query_sequence
    md = alignment.get_tag("MD")

    limit = extend + 1
    tags = []

    for miRNA_name, miRNA_start, miRNA_end in intersecting_mirna:
        shift_5p = alignment.reference_start - miRNA_start + 1
        shift_3p = alignment.reference_end - miRNA_end

        if -limit < shift_5p < limit and -limit < shift_3p < limit:
            tags.append(
                f"{miRNA_name}|{shift_5p}|{shift_3p}|{cigar}|{md}|{seq}"
            )

    return set(tags)


def main(args) -> None:
    """Add intersecting feature(s) into a SAM file as a tag."""
    try:
        validate_first_n(intersect_file=args.intersect, n=10)
    except FileFormatError as err:
        raise err

    intersect_data = parse_intersect_output(
        args.intersect, args.id, args.extension
    )

    with pysam.AlignmentFile(args.sam, "r") as samfile:
        sys.stdout.write(str(samfile.header))

        if intersect_data is None:
            return

        for alignment in samfile:
            alignment_id = alignment.query_name
            intersecting_miRNAs = intersect_data[alignment_id]

            tags = get_tags(
                intersecting_mirna=intersecting_miRNAs,
                alignment=alignment,
                extend=args.extension,
            )

            alignment.set_tag("YW", ";".join(tags))
            sys.stdout.write(alignment.to_string() + "\n")


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
