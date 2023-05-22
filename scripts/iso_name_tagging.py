import argparse
from collections import defaultdict, namedtuple
from pathlib import Path
import sys
from typing import Dict, Optional

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
        '-b', '--bed',
        help=(
            "Path to the BED file. This file must be the output of "
            "a bedtools intersect call with -a being a GFF3 file and"
            "-b a BAM file."
        ),
        type=Path,
        required=True
    )
    parser.add_argument(
        '-s', '--sam',
        help="Path to the SAM input file.",
        type=Path,
        required=True
    )
    parser.add_argument(
        '-e', '--extension',
        help=(
            "Number of nucleotides the start and end coordinates of the"
            "annotated feature has been extended. Default: %(default)d.")
            ,
        default=0,
        type=int
    )
    parser.add_argument(
        '--id',
        help=(
            "ID used to identify the feature in the output table."
            "The ID must be in lowercase. Default: %(default)s."
        ),
        default="name",
        type=str
    )

    return parser

def attributes_dictionary(attr: str) -> Optional[Dict[str, str]]:
    """Create attributes dicctionary."""
    pairs = attr.split(';')

    if len(pairs[0].split('=')) == 2:
        attr_dict = {p.split('=')[0].lower(): p.split('=')[1] for p in pairs}
    else:
        attr_dict = {p.split('"')[0].strip().lower(): p.split('"')[1] for p in pairs}

    return attr_dict

def parse_intersect_output(intersect_file: Path, id: str = "name", extension: int = 0) -> defaultdict(list):
    intersect_data = defaultdict(list)
    Fields = namedtuple('Fields',("feat_chr", "source", "feat_type",
                                    "feat_start", "feat_end", "feat_score", 
                                    "strand", "phase", "feat_attributes", 
                                    "read_chr", "read_start", "read_end", 
                                    "read_name", "read_score", "read_strand"))
    
    with open(intersect_file, 'r') as bedfile:
        for line in bedfile:

            fields = Fields(*line.strip().split('\t'))

            miRNA_name = attributes_dictionary(fields.feat_attributes)[id]
            miRNA_start = int(fields.feat_start) + extension
            miRNA_end = int(fields.feat_end) - extension

            intersect_data[fields.read_name].append((miRNA_name, miRNA_start, miRNA_end))
    
    if not intersect_data:
        return None  
    else:
        return intersect_data

def get_tags(intersecting_mirna, alignment):
    cigar = alignment.cigarstring
    md = alignment.get_tag('MD')

    tags = []

    for miRNA_name, miRNA_start, miRNA_end in intersecting_mirna:
        shift_5p = alignment.reference_start - miRNA_start
        shift_3p = miRNA_end - alignment.reference_end
        tags.append(f'{miRNA_name}|{shift_5p}|{shift_3p}|{cigar}|{md}')

    return tags

def main(args) -> None:
    intersect_data = parse_intersect_output(args.bed, args.id, args.extension)

    with pysam.AlignmentFile(args.sam, 'r') as samfile:
        
        sys.stdout.write(str(samfile.header))

        if intersect_data is None:
            return


        for alignment in samfile:
            alignment_id = alignment.query_name
            intersecting_miRNAs = intersect_data.get(alignment_id, [])

            tags = get_tags(intersecting_mirna=intersecting_miRNAs, alignment=alignment)

            alignment.set_tag('YW', ';'.join(tags))
            sys.stdout.write(alignment.to_string() + '\n')

if __name__ == "__main__":

    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
