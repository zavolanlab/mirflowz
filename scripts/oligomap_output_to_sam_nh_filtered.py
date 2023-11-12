#!/usr/bin/env python

"""Transform oligomap output FASTA file to SAM keeping the best alignments.

Read the input file, pick the best alignment(s) for each read (i.e. all
of the alignments have either 0 or 1 error). In addition, if the `--nh-filter`
CLI argument is set, filter the reads with more hits than the number provided.
The following sections will cover what the input file must be like and how the
output can look like.

EXPECTED INPUT
The FASTA file must be the output of mapping your library with the tool
oligomap: "https://github.com/zavolanlab/oligomap". Refer to the "Output
format" seccion in its main README.md for more information.
In addition, alignments have to be sorted by the read name. An example of how
an entry would look like in the EXAMPLES section.

OUTPUT FORMAT
The output consist on the filtered set of alignments from the input file in
SAM format. Only the best alignments per read (i.e. either all the alignments
have 0 or 1 error) are written to the standard output. Moreover, if the
`--nh-filter` CLI argument is given, reads with more htis than the provided
value are discarded. If none of the alignments meet the criteria, nothing is
returned. The fields in the SAM entry are:

field | value
--------------
QNAME | Read's name
FLAG  | Set to 0 if the reference sequence is found in the poisitive strand,
      | and to 16 otherwise.
RNAME | Refernce sequence name
POS   | Alignment's first position in the reference sequence
MAPQ  | Value set to 255 (mapping quality not available)
CIGAR | Alignment's CIGAR string
RNEXT | Value set to * (unavailable information)
PNEXT | Value set to 0 (unavailable information)
TLEN  | Value set to 0 (unavailable information)
SEQ   | Mapped sequence
QUAL  | Value set to * (Phred-scaled base quality not stored)
NM:i: | Alignment's edit distance (Optional field)
MD:Z: | Alignment's MD string (Optional field)
NH:i: | Alignment's NH (Optional field)

EXAMPLES
    input format:
        82-1 (23 nc) 1..23	19	44377..44398
        19
        errors: 1 orientation: +
        CTACAAAGGGAAGCACTTGTCTC
        |||||||||||||||||| ||||
        CTACAAAGGGAAGCACTT-TCTC


        97-1 (22 nc) 1..22	19	5338..5359
        19
        errors: 1 orientation: +
        TCAAAACTGAGGTGCATTTTCT
        |||||||||||| |||||||||
        TCAAAACTGAGGGGCATTTTCT

    output format:
        82-1	0	19	44377	255	18M1I4M	*	0	0	CTACAAAGGGAAGCACTTGTCTC	*	NM:i:1	MD:Z:23 	NH:i:1
        97-1	0	19	5338	255	22M 	*	0	0	TCAAAACTGAGGTGCATTTTCT	*	NM:i:1	MD:Z:12G9	NH:i:1


Paula Iborra. Zavolan Lab.
Adapted version of Alessandro Crippa script.
Refactored and documented by Iris Mestres.
"""  # noqa: E501
# pylint: enable=line-too-long

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
import re
import sys
from typing import Dict, NamedTuple


class Fields(NamedTuple):
    """Class to store an alignment in its different SAM fields."""

    read_name: str
    flag: str
    ref_seq_name: str
    position_in_ref: str
    map_q: str
    cigar_str: str
    r_next: str
    p_next: str
    t_len: str
    sequence: str
    qual: str
    edit_dist: str
    md_str: str


def parse_arguments():
    """Command-line arguments parser."""
    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-v', '--version',
        action="version",
        version="%(prog)s 1.1.0",
        help="Show program's version number and exit"
    )
    parser.add_argument(
        'infile',
        help="Path to the FASTA file resulting from the oligomap mapping.",
        type=Path,
    )
    parser.add_argument(
        '-n', '--nh-filter',
        help=(
            "Add NH tag to output and remove reads that contain more "
            "aligments than the provided NH value (with min error)."
        ),
        type=int,
    )

    return parser


def get_cigar_md(errors: str, sequence: str, bars_line: str,
                 ref_seq: str) -> tuple[str, str]:
    """Get the CIGAR and MD strings.

    Given the read and target sequences, the number of errors and the
    alignment's bar representation, this function first checks if there is any
    error in the mapping. If there is no errors, the CIGAR and MD strings are
    the read's length followed by and "M" and the read's length respectively.
    If there is an error, the function checks if the read has a deletion (the
    read sequence contains a dash), an insertion (the reference sequence
    contains a dash), or a single-point mutation. Finally, it checks the
    error's position in the read sequence and builds the CIGAR and MD strings.

    Args:
        errors:
            number of mapping errors (either 0 or 1)
        sequence:
            read sequence
        bars_line:
            alignment's bar representation
        ref_seq:
            reference sequence

    Returns:
        cigarStr:
            alignment's CIGAR string
        matchingString:
            alignment's MD string with its corresponding tag

    Examples:
       read sequence ->               GAAGGCGCTTCACCTTTGGAGT
       alignment representation ->    ||||||||||||||||||||||
       reference sequence ->          GAAGGCGCTTCACCTTTGGAGT
       errors -> 0

       CIGAR -> 21M
       MD    -> MD:Z:21

       ------------------------------------------------------

       read sequence ->               GAAGGCGCTTCACCTTTGGAGT
       alignment representation ->    ||||||||||| ||||||||||
       reference sequence ->          GAAGGCGCTTC-CCTTTGGAGT
       errors -> 1

       CIGAR -> 11M1I10M
       MD    -> MD:Z:21

       ------------------------------------------------------

       read sequence ->               GAAGGCGCTTC-CCTTTGGAGT
       alignment representation ->    ||||||||||| ||||||||||
       reference sequence ->          GAAGGCGCTTCCCCTTTGGAGT
       errors -> 1

       CIGAR -> 11M1D10M
       MD    -> MD:Z:11C10

       ------------------------------------------------------

       read sequence ->               GAAGGCGCTTCACCTTTGGAGT
       alignment representation ->    ||||||||||| ||||||||||
       reference sequence ->          GAAGGCGCTTCCCCTTTGGAGT
       errors -> 1

       CIGAR -> 21M
       MD    -> MD:Z:11C10
    """
    seq_len = len(sequence)

    if errors == '0':
        return f"{seq_len}M", f"MD:Z:{seq_len}"

    # CASE 1: deletion in the read
    if '-' in sequence:
        indelerr = "1D"

        if bars_line[0] == ' ':
            cigarStr = f"{indelerr}{seq_len - 1}M"
            matchingString = f"MD:Z:^{ref_seq[0]}{seq_len - 1}"

        elif bars_line[-1] == " ":
            cigarStr = f"{seq_len - 1}M{indelerr}"
            matchingString = f"MD:Z:{seq_len - 1}^{ref_seq[seq_len - 1]}0"

        else:
            idx = bars_line.index(' ')
            cigarStr = f"{idx}M{indelerr}{bars_line.count('|') - idx}M"
            matchingString = f"MD:Z:{idx}^{ref_seq[idx]}{seq_len - idx -1}"

        return cigarStr, matchingString

    # CASE 2: insertion in the read
    if '-' in ref_seq:
        indelerr = "1I"

        if bars_line[0] == " ":
            cigarStr = f"{indelerr}{seq_len - 1}M"

        elif bars_line[-1] == " ":
            cigarStr = f"{seq_len - 1}M{indelerr}"

        else:
            idx = bars_line.index(' ')
            cigarStr = f"{idx}M{indelerr}{bars_line.count('|') - idx}M"

        return cigarStr, f"MD:Z:{seq_len}"

    # CASE 3: single point mutation
    if bars_line[0] == " ":
        matchingString = f"MD:Z:{ref_seq[0]}{seq_len - 1}"

    elif bars_line[-1] == " ":
        matchingString = f"MD:Z:{seq_len - 1}{ref_seq[seq_len - 1]}"

    else:
        idx = bars_line.index(' ')
        matchingString = f"MD:Z:{idx}{ref_seq[idx]}{seq_len - idx - 1}"

    return f"{seq_len}M", matchingString


def get_sam_fields(aln: list[str]) -> Fields:
    """Create the read's alignment in SAM format.

    Given the set of lines oligomap provides as an alignment
    representation, this function returns a Fields class object with the
    mandatory fields in a SAM file and three optional ones.

    field | value
    --------------
    QNAME | Read's name
    FLAG  | Set to 0 if the reference sequence is found in the poisitive
          | strand, and to 16 otherwise.
    RNAME | Refernce sequence name
    POS   | Alignment's first position in the reference sequence
    MAPQ  | Value set to 255 (mapping quality not available)
    CIGAR | Alignment's CIGAR string
    RNEXT | Value set to * (unavailable information)
    PNEXT | Value set to 0 (unavailable information)
    TLEN  | Value set to 0 (unavailable information)
    SEQ   | Mapped sequence
    QUAL  | Value set to * (Phred-scaled base quality not stored)
    NM:i: | Alignment's edit distance (Optional field)
    MD:Z: | Alignment's MD string (Optional field)
    NH:i: | Alignment's NH (Optional field)

    Args:
        aln:
            list made up of the lines that form the output alignment after
            mapping with oligomap.

    Returns:
        fields:
            Fields class NamedTuple with the alignment's data as SAM fields.
    """
    seq_name_pos = aln[0].split()
    errors = aln[2].split()[1]
    seq = aln[3].strip()

    cigar, md = get_cigar_md(errors, seq, aln[4][:-1], aln[5].strip())

    fields = Fields(seq_name_pos[0],
                    '0' if aln[2].split()[3] == '+' else "16",
                    aln[1].strip(),
                    seq_name_pos[5].split('.')[0],
                    "255", cigar, '*', '0', '0',
                    re.sub('-', '', seq), '*',
                    "NM:i:0" if errors == '0' else "NM:i:1", md)

    return fields


def eval_aln(nhfilter: int, d: Dict[str, list], minErr_nh: Dict[str, list],
             fields: Fields) -> None:
    """Evaluate an alignment to add, discard or write it to the STDOUT.

    Given a read's alignment, this function first checks if the dictionary
    storing the read alignments is empty. If it is, the function checks if the
    read has an entry in the dictionary storing the current minium error and
    the NH. If there is not, or the minimum number of errors is higher than
    the alignment under evaluation, the entry for that read in both
    dictionaries is set to the data of the alignment under evaluation.

    If the alignments dictionary is not empty, the function checks if the
    read under evaluation is the first in the dictionary. If it is, and the
    number of errors in the current alignment is the same as the minimum error
    for that read, the number of hits is increased by 1. If the number of hits
    does not exceed the maximum NH or the NH is not provided, the alignment
    under evaluation is appended to the read's entry. If a maximum NH is set,
    and the read exceeds it, it is removed. If the minimum number of errors is
    higher than the number of errors of the alignment under evaluation, the
    entry for that read in both dictionaries is set to the data of the
    current alignment.

    If the read is not the first in the dictionary, the function writes all
    the alignments for the read in the first position to the standard output,
    and both dictionaries are reset to only include the data for the
    alignments under evaluation.

    Args:
        nhfilter:
            maximum number of hits an read can have to be kept
        d:
            dictionary with read names as keys and a list with Fields class
            NamedTuples as values
        minErr_nh:
            dictionary with read names as keys and a list with the current
            minimum error and the number of hits in this order as values
        fields:
            read's alignment as a Fields class NamedTuple
    """
    seq_name = fields.read_name
    errors = fields.edit_dist[-1]

    if len(d) == 0:
        if (seq_name not in list(minErr_nh.keys()) or
           errors < minErr_nh[seq_name][0]):

            minErr_nh[seq_name] = [errors, 1]
            d[seq_name] = [fields]
    else:
        if seq_name == list(d.keys())[0]:
            if errors == minErr_nh[seq_name][0]:
                minErr_nh[seq_name][1] += 1

                if nhfilter:

                    if minErr_nh[seq_name][1] <= nhfilter:
                        d[seq_name].append(fields)

                    else:
                        d.clear()
                        sys.stderr.write(f"Filtered by NH | Read {seq_name}" +
                                         f" | Errors = {errors}\n")

                else:
                    d[seq_name].append(fields)

            elif errors < minErr_nh[seq_name][0]:
                sys.stderr.write(f"Filtered by ERROR | Read {seq_name}" +
                                 f" | Errors = {minErr_nh[seq_name][0]}\n")

                minErr_nh[seq_name] = [min(errors, minErr_nh[seq_name][0]), 1]
                d[seq_name] = [fields]
        else:
            for seq, aln in d.items():
                sys.stderr.write(f"Written read {seq} | " +
                                 f"Errors = {minErr_nh[seq][0]} | " +
                                 f"NH = {minErr_nh[seq][1]}\n")
                for field in aln:
                    sys.stdout.write('\t'.join(field) +
                                     f"\tNH:i:{minErr_nh[seq][1]}" + '\n')

            d.clear()
            minErr_nh.clear()

            d[seq_name] = [fields]
            minErr_nh[seq_name] = [errors, 1]


def main(arguments) -> None:
    """Convert the alignments in the oligomap output file to its SAM format."""
    read_seqs: Dict[str, list] = {}
    seq_minError_nh: Dict[str, list] = {}

    with open(arguments.infile, 'r', encoding="utf-8") as in_file:

        sys.stderr.write("##############\nSTART READING...\n##############\n")

        lines = [in_file.readline() for _ in range(6)]
        i = 1

        while lines[0] != "":

            fields = get_sam_fields(lines)

            sys.stderr.write(f"Record:{i} | Sequence:{fields.read_name}\n")

            eval_aln(arguments.nh_filter, read_seqs, seq_minError_nh,
                     fields)
            i += 1

            in_file.readline()
            lines = [in_file.readline() for _ in range(6)]

    if len(read_seqs) > 0:

        for read_name, alignments in read_seqs.items():
            sys.stderr.write(f"Printed read {read_name} | Errors = " +
                             f"{seq_minError_nh[read_name][0]} | " +
                             f"NH = {seq_minError_nh[read_name][1]}\n")

            for aln in alignments:
                sys.stdout.write('\t'.join(aln) +
                                 f"\tNH:i:{seq_minError_nh[read_name][1]}" +
                                 '\n')

    sys.stderr.write("SUCCESSFULLY FINISHED.")


if __name__ == "__main__":
    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
