#!/usr/bin/env python

# pylint: disable=line-too-long

"""Add intersecting feature(s) into a SAM file as a tag.

Build new names for the intersecting features from a BED file and add them as a
tag to alignments in a SAM file using the format
FEATURE_ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ. If either the BED or the SAM
file is empty, only the SAM file header is returned.

EXPECTED INPUT FILES
The BED file must be the output a bedtools intersect call with `-a` being a
GFF3 file and `-b` a BAM file. If the GFF3 used in the bedtools intersect call
has the features start and end coordinates extended, the number of additional
nucleotides can be specified using the CLI option `--extension`. The SAM file
must contain only the reads that have an intersecting feature.

NAME CREATION and TAG ADDITION
For each alignment, the name of the intersecting feature will follow the
format FEATURE_ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ. The CLI option `--id`
specifies the feature identifier to be used as FEATURE_ID from within the
attributes column in the BED file. The 5p-shift and the 3-p shift values are
the difference between the feature start and end coordinates and the alignment
start and end coordinates. If `--extension` is provided, the feature start
position are adjusted by adding the given value and subtracting it from the
end position. If both, the 5p-shift and the 3p-shift, are within the range
+/- extension + 1 the feature name is added to the alignment as the new tag
"YW". Multiple intersecting feature names are separated by a semi-colon.

EXAMPLES
    Example 1
    in BED record:
        19	.	miRNA	44377	44398	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160	19	44376	44398	13-1_1	1	+	22
    in SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0
    command:
        iso_name_tagging.py -b BED -s SAM
    new name:
        hsa-miR-524-5p|0|0|11M3I11M|22
    out SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0	YW:Z:hsa-miR-524-5p|0|0|11M3I11M|22

    Example 2
    in BED record:
        19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    in SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    command:
        iso_name_tagging.py -b BED -s SAM
    new name:
        ""
    out SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0	YW:Z:

    Example 3
    in BED record:
        19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    in SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    command:
        iso_name_tagging.py -b BED -s SAM --extension 6
    new name:
        hsa-miR-1323|0|-1|21M|21
    out SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0	YW:Z:hsa-miR-1323|0|-1|21M|21

    Example 4
    in BED record:
        19	.	miRNA	44377	44398	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160	19	44376	44398	13-1_1	1	+	22
    in SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0
    command:
        iso_name_tagging.py -b BED -s SAM --id id
    new name:
        MIMAT0002849|0|0|11M3I11M|22
    out SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0	YW:Z:MIMAT0002849|0|0|11M3I11M|22
"""  # noqa: E501
# pylint: enable=line-too-long

import argparse
from collections import defaultdict, namedtuple
from pathlib import Path
import sys
from typing import Dict, Optional

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
        "-b",
        "--bed",
        help=(
            "Path to the BED file. This file must be the output of "
            " a bedtools intersect call with -a being a GFF3 file and"
            " -b a BAM file."
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


def attributes_dictionary(attr: str) -> Dict[str, str]:
    """Create attributes dicctionary."""
    pairs = attr.split(";")

    if len(pairs[0].split("=")) == 2:
        attr_dict = {p.split("=")[0].lower(): p.split("=")[1] for p in pairs}
    else:
        attr_dict = {
            p.split('"')[0].strip().lower(): p.split('"')[1] for p in pairs
        }

    return attr_dict


def parse_intersect_output(
    intersect_file: Path, ID: str = "name", extension: int = 0
) -> Optional[Dict[Optional[str], list]]:
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
    Fields = namedtuple(
        "Fields",
        (
            "feat_chr",
            "source",
            "feat_type",
            "feat_start",
            "feat_end",
            "feat_score",
            "strand",
            "phase",
            "feat_attributes",
            "read_chr",
            "read_start",
            "read_end",
            "read_name",
            "read_score",
            "read_strand",
            "overlap_len",
        ),
    )

    with open(intersect_file, "r", encoding="utf-8") as bedfile:
        for line in bedfile:
            fields = Fields(*line.strip().split("\t"))

            miRNA_name = attributes_dictionary(fields.feat_attributes)[ID]
            miRNA_start = int(fields.feat_start) + extension
            miRNA_end = int(fields.feat_end) - extension

            intersect_data[fields.read_name].append(
                (miRNA_name, miRNA_start, miRNA_end)
            )

    if not intersect_data:
        return None

    return intersect_data


def get_tags(
    intersecting_mirna: list, alignment: pysam.AlignedSegment, extension: int
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
            list with the feature name, start and end positions
        alignment:
            alignment to create the tag for
        extension:
            maximum number of nucleotides the alignment start and end positions
            can differ from the feature to count it as an intersecting feature

    Returns:
        tags: set of strings containing the new tag
    """
    cigar = alignment.cigarstring
    seq = alignment.query_sequence
    md = alignment.get_tag("MD")

    limit = extension + 1
    tags = []

    for miRNA_name, miRNA_start, miRNA_end in intersecting_mirna:
        shift_5p = alignment.reference_start - miRNA_start + 1
        shift_3p = alignment.reference_end - miRNA_end

        if -limit < shift_5p < limit and -limit < shift_3p < limit:
            tags.append(
                f"{miRNA_name}|{shift_5p}|{shift_3p}|{cigar}|{md}|{seq}"
            )

    return set(tags)


def main(arguments) -> None:
    """Add intersecting feature(s) into a SAM file as a tag."""
    intersect_data = parse_intersect_output(
        arguments.bed, arguments.id, arguments.extension
    )

    with pysam.AlignmentFile(arguments.sam, "r") as samfile:
        sys.stdout.write(str(samfile.header))

        if intersect_data is None:
            return

        for alignment in samfile:
            alignment_id = alignment.query_name
            intersecting_miRNAs = intersect_data[alignment_id]

            tags = get_tags(
                intersecting_mirna=intersecting_miRNAs,
                alignment=alignment,
                extension=arguments.extension,
            )

            alignment.set_tag("YW", ";".join(tags))
            sys.stdout.write(alignment.to_string() + "\n")


if __name__ == "__main__":
    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
