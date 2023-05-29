#!/usr/bin/env python

"""Classify and tabulate SAM file with the intersecting features name in a tag.

Read the input SAM file, calculate the sum of the intersecting contributions
for each feature and print the result in tab-delimited format. The name of the
intersecting feature(s) should be stored in the specified tag `--tag` with
the format feat_name|5p-shift|3p-shift|CIGAR|MD. The final name of the feature
is determined by the 5p and 3p shifts as well as the CIGAR and MD strings.
If both shifts equal 0 and the CIGAR (excluding the last character) and
the MD strings are the same, the read is considered the canonical feature and
only the `feat_name` is returned. Otherwise the entire name is returned. If
there are different features in the tag separated by a semi-colon, they are
treated as a unique name. The appropriate contribution of each
alignment is computed as the number of reads divided by NH. The flags
`--collapsed` and `--nh` determine if these values are taken from the alignment
name. If both flags are set, the contribution of each alignment is
computed as # of reads/NH and the values are extracted from the name. If only
`--collapsed` is set, the contribution is # of reads/NH_tag. If only `--nh` is
set, the contribution is 1/NH. Otherwise, the contribution is 1/NH_tag. If
`--nh` is not set and the NH tag is not found, the value is set to 1. If the
SAM file is empty, an empty file is produced. The output columns are also
determined by the flags `--read-ids` and `--len`. If `--len` is set, the
length of the species sequence is add as a new column in the output table.
If `--read-ids` is set, a semicolon-separated list of all alignment IDs that
overlap with the feature will always be in the last column. The output file(s)
are determined by the arguments `--mir-list` and `--lib`. The `--mir-list`
determines which tables are created and how the classification is done. If the
value of `--mir-list` is 'iso_mirna', a single table is created that includes
both isomiR and canonical miRNAs. If the value is either 'isomir' or 'mirna',
a single table is created containing only the specified feature type. If the
value is 'isomir' and 'mirna', two separate tables are created for each feature
type. The `--lib` argument determines the suffix of the output table. If it is
not specified, the suffix 'lib' is used.
"""

import argparse
from pathlib import Path
from typing import Optional

import pysam


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description=__doc__
        )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s 1.0',
        help="Show program's version number and exit"
    )
    parser.add_argument(
        'samfile',
        help=(
            "Path to the SAM file sorted by a the tag containing read"
            "intersecting species name."
        ),
        type=Path
    )
    parser.add_argument(
        '--outdir',
        help="Path to the output directory. Default: %(default)s.",
        default=Path.cwd(),
        type=Path
    )
    parser.add_argument(
        '--mir-list',
        help=(
            "List of miRNA species to have in the different output tables."
            "Default: %(default)s"
        ),
        nargs='*',
        default=['iso_mirna'],
        type=str
    )
    parser.add_argument(
        '--lib',
        help=(
            "Library to which the alignments belong to. Default : %(deafult)s"
        ),
        type=str,
        default="lib",
    )
    parser.add_argument(
        '-t', '--tag',
        help=(
            "Indicate the tag storing the read species name in uppercase."
            "Default %(default)s."
        ),
        default='YW',
        type=str
    )
    parser.add_argument(
        '-c', '--collapsed',
        help=(
            "Indicate that the SAM file has the reads collapsed by sequence"
            "and alignment. The collapsed name must be build by the alignment"
            "name followed by a '-' and the number of collpased alignments,"
            "i.e 1-4. Default %(default)s."
        ),
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--nh',
        help=(
            "Indicate that the SAM file has the NH tag in the read query name."
            "The name must be build by the alignment name followed by an"
            "underscore and the NH value, i.e 1-2_4. Default %(default)s."
        ),
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--len',
        help=(
            "Include the species length in the output table."
            "Default: %(default)s."
        ),
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--read-ids',
        help=(
            "Include read IDs of the alignments that are the same species in"
            "the output table. Default: %(default)s."
        ),
        action='store_true',
        default=False
    )

    return parser


def get_out_dirs(mir_list: list[str], library: str, outdir: Path) -> Optional[Path]:
    """Get output directories.
    
    Create directories for each miRNA type specified in a list and return a
    the new directories. If a miRNA type is not found, its corresponding value
    in the returned list is None.

    Args:
        mir_list:
            list with the desired miRNA types tables
        library:
            library where the species to be classified belong to
        outdir:
            path to the output directory
    Returns:
        the new directories and/or None
    """
    mir_out = None
    isomir_out = None
    iso_mirna_out = None

    if 'mirna' in mir_list:
        mir_out = outdir/f'mirna_counts_{library}'
    if 'isomir' in mir_list:
        isomir_out = outdir/f'isomir_counts_{library}'
    if 'iso_mirna' in mir_list:
        iso_mirna_out = outdir/f'iso_mirna_counts_{library}'

    return mir_out, isomir_out, iso_mirna_out


def get_contribution(aln: pysam.AlignedSegment, collapsed: bool = False, nh: bool = False) -> float:
    """Get contribution of an alignment to the overall count.

    Calculate the contribution of the alignment as the number of reads divided
    by NH. The arguments `collapsed` and `nh` determine if these values are
    taken from the alignment name. If both arguments are set to `True`, the
    contribution of each alignment is computed as the number of reads/NH and
    the values are extracted from the name. If only `collapsed` is `True`, the
    contribution is the number of reads/NH tag. If only `nh` is `True`, the
    contribution is 1/NH. Otherwise, the contribution is 1/NH tag. If `nh` is
    `False` and the alignment misses the NH tag, the value is set to 1.

    Args:
        aln:
            alignment to calculate the contribution from
        collapsed:
            indicates if the alignment is collapsed or not
        nh:
            indicates if the NH vale is found on the alignment name
    Returns:
        the contribution of the alignment to the overall total
    """
    query_id = aln.query_name

    if collapsed and nh:
        num_reads = int(query_id.split('-')[1].split('_')[0])
        nh_value = int(query_id.split('-')[1].split('_')[1])

    elif not collapsed and nh:
        num_reads = 1
        nh_value = int(query_id.split('_')[1])

    elif collapsed and not nh:
        num_reads = int(query_id.split('-')[1])
        try:
            nh_value = aln.get_tag("NH")
        except KeyError:
            nh_value = 1

    else:
        num_reads = 1
        try:
            nh_value = aln.get_tag("NH")
        except KeyError:
            nh_value = 1

    return num_reads/nh_value


def get_name(pre_name: str) -> list[str]:
    """Get the final name for the spieces name.

    Take a string and processes it to obtain the final name for the species
    and the type of miRNA the string belongs to. Only the feat_name is
    returned if the 3p-shift and 5p-shift are 0 and the CIGAR and MD are the
    same - excluding the last character in the CIGAR string. Otherwise, the
    whole input string is returned.

    Args:
        pre_name:
            string with the format feat_name|5p-shift|3p-shift|CIGAR|MD

    Returns:
        list with the species name to be found in the final table and its type
    """
    data_name = pre_name.split("|")

    if data_name[1] == '0' and data_name[2] == '0' and data_name[3][:2] == data_name[4]:
        return ['canonical', data_name[0]]

    return ['isomir', pre_name]


def write_output(name: str, species: list[str], mir_out: Path, isomir_out: Path, iso_mirna_out: Path) -> None:
    """Write the output to the correct file."""
    if mir_out:
        with open(mir_out, 'a') as mir:
            if name == 'canonical':
                mir.write('\t'.join(species) + '\n')
            else:
                mir.write('')

    if isomir_out:
        with open(isomir_out, 'a') as isomir:
            if name == 'isomir':
                isomir.write('\t'.join(species) + '\n')
            else:
                isomir.write('')

    if iso_mirna_out:
        with open(iso_mirna_out, 'a') as isomirna:
            if name == '':
                isomirna.write('')
            else:
                isomirna.write('\t'.join(species) + '\n')


def main(args) -> None:
    """Classify and tabulate a SAM file."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mir_out, isomir_out, iso_mirna_out = get_out_dirs(mir_list=args.mir_list,
                                                      library=args.lib,
                                                      outdir=outdir)

    read_ID = []
    count = 0

    with pysam.AlignmentFile(args.samfile, 'r') as samfile:

        try:
            first_aln = next(samfile)
            current_species = first_aln.get_tag(args.tag)
            if current_species:
                read_ID.append(first_aln.query_name)
                count = get_contribution(first_aln, args.collapsed, args.nh)

        except StopIteration:
            write_output(name="", species=[],
                         mir_out=mir_out,
                         isomir_out=isomir_out,
                         iso_mirna_out=iso_mirna_out)
            return

        for alignment in samfile:

            if current_species == '':
                current_species = alignment.get_tag(args.tag)
                continue

            if current_species == alignment.get_tag(args.tag):
                count += get_contribution(alignment, args.collapsed, args.nh)
                if args.read_ids:
                    read_ID.append(alignment.query_name)

            else:
                name = get_name(current_species)
                species = [name[1], str(count)]
                if args.len:
                    species.append(str(alignment.query_alignment_length))
                if args.read_ids:
                    species.append(';'.join(read_ID))

                write_output(name[0], species, mir_out=mir_out,
                             isomir_out=isomir_out,
                             iso_mirna_out=iso_mirna_out)

                current_species = alignment.get_tag(args.tag)
                count = get_contribution(alignment, args.collapsed, args.nh)
                read_ID = [alignment.query_name]

        name = get_name(current_species)
        species = [name[1], str(count)]
        if args.len:
            species.append(str(alignment.query_alignment_length))
        if args.read_ids:
            species.append(';'.join(read_ID))

        write_output(name[0], species, mir_out=mir_out,
                     isomir_out=isomir_out,
                     iso_mirna_out=iso_mirna_out)


if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
