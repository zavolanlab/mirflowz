#!/usr/bin/env python

# pylint: disable=line-too-long

"""Add intersecting feature(s) into a SAM file as a tag.

Build new names for the intersecting features from a BED file and add them as
a tag to alignments in a SAM file using the format
FEATURE_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ.

EXPECTED INPUT FILES
The BED file must be the output of a bedtools intersect call with `-a` being a
GFF3 file and `-b` a BAM file. The SAM file must only
contain alignments with an intersecting feature. If either the BED or the SAM
file is empty, only the SAM file header is returned.

NAME CREATION and TAG ADDITION
For each alignment, the name of the intersecting feature follows the
format FEATURE_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ. The CLI option `--id`
specifies the feature identifier to be used as FEATURE_ID from the attributes
column in the BED file. The 5' and 3' shift values are the difference between
the alignment and its intersecting feature(s) start and end coordinates
respectively. If `--extension` is provided, features start and end positions
are adjusted by adding and subtracting respectively the given value. If
`--shift` is provided, and both, the 5' and 3'-end shifts, are within the
range +/- `--shift` the feature name is added to the alignment as the new tag
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
        ""
    out SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0	YW:Z:
    explanation:
        The aligned read has a length of 22 whereas the feature has a length of
        21. As both have the same starting position, there is a 1bp shift at
        the 3' end. Given that no extension nor shift are provided in the
        script call, no adjustments to the annotated coordinates are made and
        no shifts are allowed. Hence, no tag is added.

    Example 2
    in BED record:
        19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    in SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    command:
        iso_name_tagging.py -b BED -s SAM
    new name:
        hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC
    out SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0	YW:Z:hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC
    explanation:
        The aligned read and the annotated featrue have the same start and end
        positions. Given that no extension are provided in the script call, no
        coordinates adjustments are made. And there is no shift on ether end,
        the new tag is added.

    Example 3
    in BED record:
        19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    in SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    command:
        iso_name_tagging.py -b BED -s SAM --extension 6 --shift 7
    new name:
        hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC
    out SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0	YW:Z:hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC
    explanation:
        The feature's start and end coordinates are adjusted by amount of bp
        specified by the CLI argument `--extension`. The alignment start and
        end positions differ from the adjusted miRNA ones by 6 bp. As there is
        a +/- 7 bp shift allowed and both shifts are within this range, the new
        tag is added.

    Example 4
    in BED record:
        19	.	miRNA	44377	44404	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160	19	44376	44401	13-1_1	1	+	22
    in SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:25	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0
    command:
        iso_name_tagging.py -b BED -s SAM --extension 6 --shift 7 --id id
    new name:
        MIMAT0002849|6|4|11M3I11M|25|CTACAAAGGGAGGTAGCACTTTCTC
    out SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:25	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0	YW:Z:MIMAT0002849|6|4|11M3I11M|25|CTACAAAGGGAGGTAGCACTTTCTC
    explanation:
        The feature's start and end coordinates are adjusted by amount of bp
        specified by the CLI argument `--extension`. The alignment start and
        end positions differ from the adjusted miRNA ones by 6 bp. As there is
        a +/- 7 bp shift allowed and both shifts are within this range, the new
        tag is added. This time using the feature ID instead of the Name.
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
        version="%(prog)s 1.0.0",
        help="Show program's version number and exit",
    )
    parser.add_argument(
        "-b",
        "--bed",
        help=(
            "Path to the BED file. This file must be the output of "
            " a bedtools intersect call with `-a` being a GFF3 file and"
            " `-b` a BAM file."
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
        "-S",
        "--shift",
        help=(
            "Absoulte difference in nucleotides allowed between the alignment"
            " and its intersecting features' start and end coordinates."
            " Default: %(default)d."
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
    """Create attributes dictionary."""
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

    Given a BED file generated by intersecting a GFF3 file (-a) with a BAM file
    (-b) using bedtools intersect, create a dictionary where the alignment
    names are the keys, and the values are lists containing the feature name,
    the start and the end positions. The `ID` argument specifies the feature
    name to use, and the `extension` argument adjusts the feature coordinates
    by adding and subtracting the given value to the start and end positions
    respectively. If the BED file is empty, `None` is returned.

    Args:
        intersect_file:
            Path to the intersect BED file.
        ID:
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

            feat_name = attributes_dictionary(fields.feat_attributes)[ID]
            feat_start = int(fields.feat_start) + extension
            feat_end = int(fields.feat_end) - extension

            intersect_data[fields.read_name].append(
                (feat_name, feat_start, feat_end)
            )

    if not intersect_data:
        return None

    return intersect_data


def get_tags(
    intersecting_feat: list, alignment: pysam.AlignedSegment, shift: int = 0
) -> set:
    """Get alignment's new tag.

    Given an alignment and a list containing the name, the (extended) start and
    end positions of all the intersecting features, create a set of strings
    to be used as a new custom tag on that alignment.
    Each string follows the format:
        FEATURE_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ
    The 5' end shift is computed as the difference between the feature and the
    alignment starting coordinate. Similarly, the 3' end shift is the result of
    subtracting the feature's end position to the alignment's one. If both
    shifts values are within the range +/- `shift`, the name is added to the
    final list. Note that `pysam` assumes 0-base index therefore, the addition
    of one base is required to compute the 5' end shift.

    Args:
        intersecting_feat:
            list with the miRNA species name, start and end positions
        alignment:
            alignment to create the tag for
        shift:
            value that sets the range in which both shifts have to be in to
            create a new string for a particular feature

    Returns:
        tags:
            set of strings for each valid intersecting feature
    """
    cigar = alignment.cigarstring
    seq = alignment.query_sequence
    md = alignment.get_tag("MD")

    tags = []

    for feat_name, feat_start, feat_end in intersecting_feat:
        shift_5 = alignment.reference_start - feat_start + 1
        shift_3 = alignment.reference_start + alignment.query_length - feat_end
        if -shift <= shift_5 <= shift and -shift <= shift_3 <= shift:
            tags.append(f"{feat_name}|{shift_5}|{shift_3}|{cigar}|{md}|{seq}")

    return set(tags)


def main(args) -> None:
    """Add intersecting feature(s) into a SAM file as a tag."""
    intersect_data = parse_intersect_output(
        intersect_file=args.bed,
        ID=args.id,
        extension=args.extension,
    )

    with pysam.AlignmentFile(args.sam, "r") as samfile:
        sys.stdout.write(str(samfile.header))

        if intersect_data is None:
            return

        for alignment in samfile:
            alignment_id = alignment.query_name
            intersecting_feats = intersect_data[alignment_id]

            tags = get_tags(
                intersecting_feat=intersecting_feats,
                alignment=alignment,
                shift=args.shift,
            )

            alignment.set_tag("YW", ";".join(tags))
            sys.stdout.write(alignment.to_string() + "\n")


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
