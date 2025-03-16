#!/usr/bin/env python

# pylint: disable=line-too-long

"""Annotate SAM alignments with their intersecting feature(s).

Add a custom tag ("YW") to each alignment in the SAM file if an intersecting
feature is found in the BED file. The BED file must be the result of a
bedtools intersect call where `-a` is a GFF3 file and `-b` a BAM file. If
either the BED or the SAM file is empty, only the SAM file header is returned.

Each matching feature is used to build a tag with the following format:

    FEATURE_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ.

Where:
    - FEATURE_ID: Extracted from the specified attribute in the BED file
    (default: "name").
    - 5'-shift: Difference between the (possibly adjusted) feature start and
    the alignment start.
    - 3'-shift: Difference between the alignment end and the (possibly
    adjusted) feature end.
    - CIGAR and MD: The alignment's CIGAR string and MD tag respectively.
    - READ_SEQ: The read sequence from the alignment.

Optional adjustments:
--extension: Adjust the feature's start and end coordinates by the given value
    (start is increased and end decreased).
--shift: Requires both the 5' and 3' shift values to be within +/- this value
    to include the feature tag in the alignment.

If an alignment has multiple intersecting features, the tag values are
concatenated using a semicolon as the separator. If no valid intersecting
feature is found, the alignment is skipped.

Examples
--------

Example 1: No tag added due to shift mismatch
    BED record:
        19	.	miRNA	44377	44398	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160	19	44376	44398	13-1_1	1	+	22
    SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0
    Command:
        annotate_sam_with_bed_features.py -b BED -s SAM
    Outcome:
        No tag is added because the alignment's 3' end shift exceeds the
        allowed range.

Example 2: Tag added without coordinate adjustments
    BED record:
        19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    Command:
        annotate_sam_with_bed_features.py -b BED -s SAM
    Outcome:
        The tag "hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC" is appended
        because the alignment perfectly matches the feature coordinates.

Example 3: Tag added with coordinate adjustments and allowed shift
    BED record:
        19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    Command:
        annotate_sam_with_bed_features.py -b BED -s SAM --extension 6 --shift 7
    Outcome:
        After adjusting the feature coordinates, both the 5' and 3' shifts are
        within the allowed range, so the tag
        "hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC" is added.

Example 4: Tag added using a different feature identifier
    BED record:
        19	.	miRNA	44377	44404	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160	19	44376	44401	13-1_1	1	+	22
    SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:25	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0
    Command:
        annotate_sam_with_bed_features.py -b BED -s SAM --extension 6 --shift 7 --id id
    Outcome:
        The feature identifier is taken from the `id` attribute (in
        lowercase), and the tag is constructed accordingly
        ("MIMAT0002849|6|4|11M3I11M|25|CTACAAAGGGAGGTAGCACTTTCTC").
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
            " a bedtools intersect call with `-a` as a GFF3 file and `-b as"
            " a BAM file."
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
            " the start and subtract from the end. Default: %(default)d."
        ),
        default=0,
        type=int,
    )
    parser.add_argument(
        "-S",
        "--shift",
        help=(
            "Maximum allowed absoulte difference (in nucleotides) between the"
            " alignment and feature coordinates at both ends. The tag is"
            " added only if both 5' and 3' shifts are within this range"
            " Default: %(default)d."
        ),
        default=0,
        type=int,
    )
    parser.add_argument(
        "--id",
        help=(
            "Feature attribute key to extract the indentifier (must be"
            " lowercase). Default: %(default)s."
        ),
        default="name",
        type=str,
    )

    return parser


def attributes_dictionary(attr: str) -> Dict[str, str]:
    """Convert an attribute string to a dictionary.

    The attribute string is expected to be a semicolon-separated list of
    key-value pairs. The function supports both '='-delimited and quoted
    formats.

    Args:
        attr: The attribute string from a BED record.

    Returns:
        A dictionary with attribute keys in lowercase.
    """
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
    """Parse the BED file produced by a bedtools intersect command.

    Read the BED file and create a dictionary mapping each alignment name to
    a list of intersecting features. Each feature is represented as a tuple
    containing:
        feature identifier, adjusted feature start and adjusted feature end

    The feature identifier is extracted from the attribute column using the
    key specified by `ID`.

    The `extension` parameter is applied to adjust the start and end
    coordinates (start is increased, end is decreased).

    Args:
        intersect_file: Path to the intersect BED file.
        ID: Attribute key to extract the feature identifier (default: "name").
        extension: Nucleotides to adjust the feature coordinates (default: 0).

    Returns:
        A dictionary mapping alignment names to lists of feature tuples. If
        the BED file is empty, returns None.
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
    """Construct a custom tag for an alignment based on intersecting features.

    For each intersecting feature (given as a tuple of identifier, start, and
    end), compute the 5'-shift and 3'-shift values relative to the alignment:
        - 5'-shift = (alignment start - feature start + 1)
        - 3'-shift = (alignment end - feature end)

    The tag for a feature is generated in the format:

        FEATURE_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ

    Only features where both shift values are within +/- `shift` are included.
    If multiple features are valid, the resulting tags are collected in a set.

    Args:
        intersecting_feat: list of tuples (feature identifier, start and end).
        alignment: A pysam.AlignedSegment object
        shift: Maximum allowed absolute shift for both ends (default: 0).

    Returns:
        A set of tag strings constructed from valid intersecting features.
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
    """Annotate alignments in a SAM file with intersecting feature tags."""
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

            if len(intersecting_feats) == 0:
                continue

            tags = get_tags(
                intersecting_feat=intersecting_feats,
                alignment=alignment,
                shift=args.shift,
            )

            if len(tags) == 0:
                continue

            alignment.set_tag("YW", ";".join(tags))
            sys.stdout.write(alignment.to_string() + "\n")


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
