#!/usr/bin/env python

# pylint: disable=line-too-long

"""Annotate SAM alignments with their intersecting feature(s).

Add a custom tag ("YW") to each alignment in the SAM file if an intersecting
feature is found in the INTERSECT file. The INTERSECT file must be the result
of the call 'bedtools intersect -wo -s -a GFF3/GTF -b BAM'. If either the
INTERSECT or the SAM file is empty, only the SAM file header is returned.

Each matching feature is used to build a tag with the following format:

    FEATURE_ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ

Where:
    - FEATURE_ID: Extracted from the specified attribute in the INTERSECT file
      (default: "name")
    - 5p-shift: Difference between the (possibly adjusted) feature start and
      the alignment start
    - 3p-shift: Difference between the alignment end and the (possibly
      adjusted) feature end
    - CIGAR and MD: The alignment's CIGAR string and MD tag respectively
    - READ_SEQ: The read sequence from the alignment

Optional adjustment:
--extension: Adjust the feature's start and end coordinates by the given value
    (start is increased and end decreased). In addition, requires both the 5'
    and 3' shift values to be within +/- this value to include the feature tag
    in the alignment

If an alignment has multiple intersecting features, the tag values are
concatenated using a semicolon as the separator. If there are no intersecting
features or none of the intersecting features pass the shift filter, the
alignment is skipped.

Examples
--------
Example 1: Feature intersects alignment; coordinates adjustment and shift allowed
    use case:
        Prior to checking if the feature is intersecting the alignment, its
        coordinates are adjusted by the value specified in `--extension`. In
        addition, the same value is used to specify the +/- shift allowed
        between the feature and the read alignment start and end coordinates.

    command:
        annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM --extension 5

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
        annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM

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
        annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM --extension 1

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
        not considered to intersect the read alignment and the alignment is
        not written in the output file.


Example 4: Feature intersects alignment; using feature's "Alias"
    use case:
        The feature and read alignment coordinates must perfectly match.

    command:
        annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM --id alias

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
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import pysam

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
            "Number of nucleotides to adjust the feature coordinates: add to"
            " the start and subtract from the end. Also, maximum allowed"
            " difference between the alignment and feature coordinates at"
            " both ends. The tag is added only if both shifts are within +/-"
            " this value. Its value has to be either 0 or positive integer."
            " Default: %(default)d."
        ),
        default=0,
        choices=range(10**6),
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
) -> Optional[Dict[str, List[Tuple[str, int, int]]]]:
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
        FEATURE-ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ
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

    try:
        md = alignment.get_tag("MD")
    except KeyError as keyerr:
        raise KeyError(
            f'SAM record "{alignment.query_name}" is missing required MD tag'
        ) from keyerr

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
    """Add intersecting feature(s) into a SAM record as a tag."""
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

            assert alignment_id is not None
            intersecting_miRNAs = intersect_data[alignment_id]

            tags = get_tags(
                intersecting_mirna=intersecting_miRNAs,
                alignment=alignment,
                extend=args.extension,
            )

            if len(tags) == 0:
                continue

            alignment.set_tag("YW", ";".join(tags))
            sys.stdout.write(alignment.to_string() + "\n")


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
