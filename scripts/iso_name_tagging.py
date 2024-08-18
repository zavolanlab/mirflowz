#!/usr/bin/env python

# pylint: disable=line-too-long

"""Add intersecting feature(s) into a SAM file as a tag.

Build new names for the intersecting features from a BED file and add them as
a tag to alignments in a SAM file using the format
miRNA_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ. If either the BED or the SAM
file is empty, only the SAM file header is returned.

EXPECTED INPUT FILES
The BED file must be the output of a bedtools intersect call with `-a` being a
GFF3 file and `-b` a BAM file. If the GFF3 used in the bedtools intersect call
has the features start and end coordinates extended, the number of additional
nucleotides must be specified using the CLI option `--extension`. The SAM file
must contain only the reads that have an intersecting feature.

NAME CREATION and TAG ADDITION
For each alignment, the name of the intersecting feature follows the
format miRNA_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ. The CLI option `--id`
specifies the feature identifier to be used as miRNA_ID from the attributes
column in the BED file. The 5' and 3' shift values are the difference between
the feature (extended) start and end coordinates and the alignment ones. If
`--extension` is provided, the feature start and end positions are adjusted by
adding and subtracting respectively the given value. If both, the 5' and
3'-end shifts, are within the range +/- extension (or equal 0 if no value is
provided) the feature name is added to the alignment as the new tag "YW".
Multiple intersecting feature names are separated by a semi-colon.

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

    Example 3
    in BED record:
        19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	48-1_1	255	+	21
    in SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
    command:
        iso_name_tagging.py -b BED -s SAM --extension 6
    new name:
        hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC
    out SAM record:
        48-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0	YW:Z:hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

    Example 4
    in BED record:
        19	.	miRNA	44377	44404	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160	19	44376	44401	13-1_1	1	+	22
    in SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:25	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0
    command:
        iso_name_tagging.py -b BED -s SAM --extension 6 --id id
    new name:
        MIMAT0002849|6|4|11M3I11M|25|CTACAAAGGGAGGTAGCACTTTCTC
    out SAM record:
        13-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:25	NH:i:1	NM:i:3	RG:Z:A1	YZ:Z:0	YW:Z:MIMAT0002849|6|4|11M3I11M|25|CTACAAAGGGAGGTAGCACTTTCTC
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
    intersecting_mirna: list, alignment: pysam.AlignedSegment, extend: int
) -> set:
    """Get tag for alignment.

    Given an alignment and a list containing the name, start and end
    (extended) positions of all the intersecting miRNA species, create the set
    of strings used as a new custom tag to that alignment.
    Each string follows the format:
        miRNA_ID|5'-shift|3'-shift|CIGAR|MD|READ_SEQ
    The 5' end shift is computed as the difference between the feature and the
    alignment positions. Similarly, the 3'-end shift is the result of
    subtracting the feature end position to the alignment's one. If both shifts
    values are within the range +/-extension, the name is added to the final
    list. Note that `pysam` assumes 0-base index therefore, the addition of
    one base is required to compute the 5' end shift.

    Args:
        intersecting_mirna:
            list with the miRNA species name, start and end positions
        alignment:
            alignment to create the tag for
        extend:
            value that sets the range in which both shifts have to be in to
            create a new string for a particular miRNA species

    Returns:
        tags:
            set of strings for each valid intersecting miRNA
    """
    cigar = alignment.cigarstring
    seq = alignment.query_sequence
    md = alignment.get_tag("MD")

    tags = []

    for miR_name, miR_start, miR_end in intersecting_mirna:
        shift_5 = alignment.reference_start - miR_start + 1
        shift_3 = alignment.reference_start + alignment.query_length - miR_end
        if -extend <= shift_5 <= extend and -extend <= shift_3 <= extend:
            tags.append(f"{miR_name}|{shift_5}|{shift_3}|{cigar}|{md}|{seq}")

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
