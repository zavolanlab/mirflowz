#!/usr/bin/env python

# pylint: disable=line-too-long

"""Quantify miRNAs and corresponding isomiRs.

Read the input SAM file, calculate the contribution sum of the intersecting
alignments for each feature and write the result in tab-delimited format. The
following sections will cover what the SAM file must be like, how the counting
is done, and how the output can look like. The final section contain examples
covering the main modes.

EXPECTED INPUT
The SAM file must contain only the alignments that intersect with canonical
miRNAs and/or isomiRs. For each alignment, the intersecting feature name(s)
must be stored in a tag. If the tag determined by the CLI argument `--tag` is
empty, the alignment is ignored. In addition, the SAM file must be sorted by
the tag storing the feature names. Optionally, the SAM file can contain
collapsed reads, as, e.g., produced by the 'fastx_collapser' tool of the
FASTX-Toolkit: http://hannonlab.cshl.edu/fastx_toolkit/". Finally, read names
can also contain the alignment NH vale (see 'READ NAME FORMAT' section)

FEATURE NAME FORMAT
The name of the intersecting feature(s) has to follow the format:
FEAT_NAME|5p-SHIFT|3p-SHIFT|CIGAR|MD|READ_SEQ.
For isomiRs, the whole name is used. For canonical miRNAs, the `feat_name`.
A feature is classified as canonical if the 5p and 3p shifts are both 0 and the
CIGAR and MD strings are the same. If there are different features names
separated by a semi-colon, they are treated as a unique isomiR name. For
example, the feature name 'feat_name_1|0|0|23M|23|read_seq' would represent
a canonical miRNA and the feature names 'feat_name_2|0|0|22M|0A0A2T17|read_seq'
and 'feat_name_3|0|0|22M|22|read_seq;feat_name_4|1|0|23M|17A3T0C0|read_seq'
would represent isomiRs.

READ NAME FORMAT
The read name must have one of the following formats:
NAME-COUNT_NH, NAME-COUNT, NAME_NH or NAME where NAME is an arbitrary unique
name, COUNT is the number of identical reads that were collapsed and NH is the
alignment NH tag. For example, the query name 'my_query_name-13_2' would
indicate that the read represented by this alignment occurred 13 times before
collapsing and the alignment's NH value is 2. Only if COUNT is in the name the
CLI flag `--collapsed` can be set. Likewise, the CLI flag `--nh` can only be
set if NH is in the name.

COUNT TABULATION
The count for each feature is the sum of all the intersecting alignment's
contribution. The appropriate contribution is the ratio between the number of
reads and the alignment's NH. The CLI flags `--collapsed` and `--nh` determine
if these values are taken from the alignment's name. If `--collapsed` is not
set, the number of reads is considered to be 1. If `--nh` is not set, the NH
value is obtained from the NH tag. If the NH tag is missing, its value is
considered to be 1.

TABULATED FEATURES
Which kind of features will appear in the output file is determined by the CLI
argument `--mir-list`. The three possible options are 'mirna' to only include
canonical miRNAs, 'isomir' to only include isomiRs, and 'mirna' and 'isomir' to
include both canonical miRNAs and isomiRs.

OUTPUT FORMAT
The output is tab-delimited file named 'mir_counts_LIB', where LIB is
the read's library name set in the CLI option `--lib`. If the SAM file is
empty, an empty file is produced. By default each row contains the feature
name and its partial count. Three extra columns can be added by using the CLI
flags `--count`, `--len` and `--read-ids`. If `--count` is set, the output
table will contain the number of best alignments for each feature. If `--len`
is set, the output table will contain the read's length. If both `--count` and
`--len` are set, the count will always be followed by the read's length. If
`--read-ids` is set, a semicolon-separated list of all alignment IDs that
overlap with the feature will always be in the last column.

EXAMPLES
    Example 1
    input: SAM with intersecting feature names in the tag XN
    command: mirna_quantification.py SAM --tag XN
    output: hsa-miR-516b-5p	3.0
            hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG	0.6000000000000001

    Example 2
    input: SAM meeting the minimal characteristics
    command: mirna_quantification.py SAM --mir-list mirna
    output: hsa-miR-517b-3p	1.3333333333333333
            hsa-miR-518b	1.0

    Example 3
    input: SAM meeting the minimal characteristics
    command: mirna_quantification.py SAM --mir-list isomir
    output: hsa-miR-512-3p|0|1|23M|22C|AAGTGCTGTCATAGCTGAGGTAA	1.6666666666666665
            hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG	0.6000000000000001

    Example 4
    input: SAM with read names following the format NAME-COUNT_NH
    command: mirna_quantification.py SAM --collpased --nh
    output: hsa-miR-516b-5p	3.0
            hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG	0.6000000000000001

    Example 5
    input: SAM with read names following the format NAME-COUNT
    command: mirna_quantification.py SAM --collapsed
    output: hsa-miR-516b-5p	3.0
            hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG	0.6000000000000001

    Example 6
    input: SAM with read names following the format NAME_NH
    command: mirna_quantification.py SAM --nh
    output: hsa-miR-516b-5p	3.0
            hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG	0.6000000000000001

    Example 7
    input: SAM with read names following the format NAME_NH
    command: mirna_quantification.py SAM --read-ids --nh --mir-list mirna
    output: hsa-miR-498-5p	4.333333333333333   270396_3
            hsa-miR-1323    12.0    673650_2;906983_4

    Example 8
    input: SAM meeting the minimal characteristics
    command: mirna_quantification.py SAM --count --len --mir-list isomir
    output: hsa-miR-512-3p|0|1|23M|22C0|AAGTGCTGTCATAGCTGAGGTAA	4.333333333333333  6 23
            hsa-miR-512-3p|0|1|23M|3T18C0|AAGGGCTGTCATAGCTGAGGTAA 12.0  8 19
"""  # noqa: E501
# pylint: enable=line-too-long

import argparse
from pathlib import Path
import sys
import re

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
        help="Show program's version number and exit.",
    )
    parser.add_argument(
        "samfile",
        help=(
            "Path to the SAM file containing the intersecting miRNA name(s)."
        ),
        type=Path,
    )
    parser.add_argument(
        "--outdir",
        help="Path to the output directory. Default: %(default)s.",
        default=Path.cwd(),
        type=Path,
    )
    parser.add_argument(
        "--mir-list",
        help=(
            "List of miRNA types to have in the output table."
            " Default: %(default)s."
        ),
        nargs="*",
        default=["isomir", "mirna"],
        type=str,
    )
    parser.add_argument(
        "--lib",
        help=(
            "Library to which the alignments belong to. Default: %(default)s."
        ),
        type=str,
        default="lib",
    )
    parser.add_argument(
        "-t",
        "--tag",
        help=(
            "Indicate the tag storing the intersecting miRNA name."
            " Default: %(default)s."
        ),
        default="YW",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--collapsed",
        help=(
            "Indicate that the SAM file has the reads collapsed by sequence."
            " In that case, the SAM query names are expected to follow the"
            " format NAME-COUNT where NAME is an aribitrary unique name and"
            " COUNT is  is the number of identical reads that were collapsed"
            " i.e my_query_name-13. The required naming format is produced,"
            " e.g., by the 'fastx_collapser' tool of the"
            " FASTX-Toolkit: http://hannonlab.cshl.edu/fastx_toolkit/"
            " Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--nh",
        help=(
            "Indicate that the SAM file has the NH value at the end of the"
            " read query name. In that case, SAM query names are expected to"
            " follow the format NAME_NH where NAME is an aribitrary unique"
            " name and NH is the NH value, i.e my_query_name_4. Default:"
            " %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--count",
        help=(
            "If set, the amount of best alignments for each miRNA is included"
            " in the output table. Default: %(default)s"
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--len",
        help=(
            "If set, the miRNA length is included in the output table."
            " Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--read-ids",
        help=(
            "If set, the read IDs that belong to each miRNA are included in"
            " the output table separated by a semi-colon."
            " Default: %(default)s."
        ),
        action="store_true",
        default=False,
    )

    return parser


def collapsed_nh_contribution(aln: pysam.AlignedSegment) -> float:
    """Get the contribution of the alignment to the overall count.

    The contribution is computed as the ratio of the number of reads collapsed
    in the alignment and the NH value. It is assumed that the alignment query
    name contains the number of collapsed reads as well as the NH value in the
    format NAME-COUNT_NH.

    Args:
        aln:
            Alignment to which the overall contribution is calculated

    Returns:
        Contribution of alignment to overall count
    """
    name = str(aln.query_name)
    values = []
    try:
        if val := re.search(r"\d+_\d+$", name):
            values = val.group().split("_")

        return float(values[0]) / float(values[1])

    except AttributeError:
        sys.stdout.write(
            f"Invalid query name: '{aln.query_name}'.\n"
            + "Cannot calculate contribution.\n"
            + "Check the SAM file validity and CLI options"
            + " --collapsed and --nh.\n"
        )
        raise


def collapsed_contribution(aln: pysam.AlignedSegment) -> float:
    """Get the contribution of the alignment to the overall count.

    The contribution is computed as the ratio of the number of reads collapsed
    in the alignment and the value stored in the NH tag. If the tag is missing,
    the NH value is 1. It is assumed that the alignment query name contains
    the number of collapsed reads in the format NAME-COUNT.

    Args:
        aln:
            Alignment to which the overall contribution is calculated

    Returns:
        Contribution of alignment to overall count
    """
    name = str(aln.query_name)
    collapsed = 0.0
    try:
        if coll := re.search(r"\d+$", name):
            collapsed = float(coll.group())

    except AttributeError as atterr:
        raise AttributeError(
            f'Invalid query name: "{aln.query_name}".\n'
            f"Option --collapsed specified but query name does not include"
            f" the number of collapsed sequences.\nCheck SAM file consistency"
            f" and CLI options --collapsed and --nh."
        ) from atterr

    try:
        nh_value = float(aln.get_tag("NH"))
        return collapsed / nh_value

    except KeyError:
        return collapsed


def nh_contribution(aln: pysam.AlignedSegment) -> float:
    """Get the contribution of the alignment to the overall count.

    The contribution is computed as the ratio of the number of reads collapsed
    in the alignment and the value stored in the NH tag. If the tag is missing,
    the NH value is 1. It is assumed that the alignment query name contains the
    NH value in the format NAME_NH.

    Args:
        aln:
            Alignment to which the overall contribution is calculated

    Returns:
        Contribution of alignment to overall count
    """
    name = str(aln.query_name)
    nh_val = 0.0
    try:
        if cont := re.search(r"\d+$", name):
            nh_val = float(cont.group())

        return 1 / nh_val

    except AttributeError as atterr:
        raise AttributeError(
            f'Invalid query name: "{aln.query_name}".\n'
            f"Option --nh specified but query name does not include NH.\n"
            f"Check SAM file consistency and CLI options --collapsed and --nh."
        ) from atterr


def contribution(aln: pysam.AlignedSegment) -> float:
    """Get the contribution of the alignment to the overall count.

    The contribution is computed as the ratio of the number of reads collapsed
    in the alignment and the value stored in the NH tag. If the tag is missing,
    the overall contribution is 1.

    Args:
        aln:
            Alignment to which the overall contribution is calculated

    Returns:
        Contribution of alignment to overall count
    """
    try:
        return 1 / float(aln.get_tag("NH"))

    except KeyError:
        return 1.0


def get_name(pre_name: str) -> list[str]:
    """Get the final name for the species name.

    Take a string and processes it to obtain the final name for the species
    and the type of miRNA the string belongs to. Only the feat_name is
    returned if the 3p-shift and 5p-shift are 0 and the CIGAR and MD are the
    same. Otherwise, the whole input string is returned.

    Args:
        pre_name:
            string with the format
            FEAT_NAME|5p-SHIFT|3p-SHIFT|CIGAR|MD|READ_SEQ

    Returns:
        list with the species name to be found in the final table and its type
    """
    data_name = pre_name.split("|")
    cigar = re.sub(r"[^0-9]", "", data_name[3])
    md = re.sub(r"[^0-9]", "", data_name[4])

    if data_name[1] == "0" and data_name[2] == "0" and cigar == md:
        return ["mirna", data_name[0]]

    return ["isomir", pre_name]


def write_output(
    name: str, species: list[str], mir_list: list[str], mirna_out: Path
) -> None:
    """Write to the output the correct miRNA type."""
    with open(mirna_out, "a", encoding="utf-8") as mirna:
        if name in mir_list:
            mirna.write("\t".join(species) + "\n")
        else:
            mirna.write("")


def main(args) -> None:
    """Quantify miRNAs and corresponding isomiRs."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"mirna_counts_{args.lib}"

    contribution_type = {
        (True, True): collapsed_nh_contribution,
        (True, False): collapsed_contribution,
        (False, True): nh_contribution,
        (False, False): contribution,
    }

    get_contribution = contribution_type[args.collapsed, args.nh]

    with pysam.AlignmentFile(args.samfile, "r") as samfile:
        try:
            alignment = next(samfile)
            current_species = alignment.get_tag(args.tag)
            if current_species:
                read_ID = [str(alignment.query_name)]
                count = get_contribution(alignment)
                alns_count = 1

        except StopIteration:
            write_output(
                name="",
                species=[],
                mir_list=args.mir_list,
                mirna_out=outfile,
            )

            return

        for alignment in samfile:
            if current_species == "":
                current_species = alignment.get_tag(args.tag)
                count = get_contribution(alignment)
                alns_count = 1
                read_ID = [str(alignment.query_name)]

                continue

            if current_species == alignment.get_tag(args.tag):
                count += get_contribution(alignment)
                alns_count += 1
                read_ID.append(str(alignment.query_name))

            else:
                name = get_name(str(current_species))
                species = [name[1], str(count)]

                if args.count:
                    species.append(str(alns_count))
                if args.len:
                    species.append(str(alignment.query_alignment_length))
                if args.read_ids:
                    species.append(";".join(read_ID))

                write_output(
                    name=name[0],
                    species=species,
                    mir_list=args.mir_list,
                    mirna_out=outfile,
                )

                current_species = alignment.get_tag(args.tag)
                count = get_contribution(alignment)
                alns_count = 1
                read_ID = [str(alignment.query_name)]

        name = get_name(str(current_species))
        species = [name[1], str(count)]

        if args.count:
            species.append(str(alns_count))
        if args.len:
            species.append(str(alignment.query_alignment_length))
        if args.read_ids:
            species.append(";".join(read_ID))

        write_output(
            name=name[0],
            species=species,
            mir_list=args.mir_list,
            mirna_out=outfile,
        )


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
