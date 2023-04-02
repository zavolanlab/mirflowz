#!/usr/bin/env python

import argparse
import pysam
import sys

def parse_arguments():
    '''
    Command-line arguments parser
    '''
    parser = argparse.ArgumentParser(
        description="Script to extend pre-miRNs overhang."
        )

    parser.add_argument(
        '-v','--version', 
        action='version', 
        version='%(prog)s 1.0',
        help="Show program's version number and exit"
        )

    parser.add_argument(
        '-i', '--input', 
        help="Path to the gff3 input file",
        type=str,
        default=None,
        )

    parser.add_argument(
        '-o', '--output', 
        help="Path for the gff3 final file.",
        type=str,
        default=None,
        )
    return parser


def count_matches(read):
    '''
    Counts the number of matches in a read based on its CIGAR string.

    Args:
        read (pysam.AlignedSegment): 
            The read to count matches for.

    Returns:
        int: The number of matches in the read.
    '''
    
    return sum([op[1] for op in read.cigar if op[0] == 0])

def count_indels(read):
    '''
    Counts the number of insertions and/or deletions in a read based on its 
    CIGAR string.

    Args:
        read (pysam.AlignedSegment): The read to count insertions for.

    Returns:
        int: The number of insertions in the read.
    '''
    return sum([op[1] for op in read.cigar if op[0] == 1 or op[0] == 2])

def filter_multimappers(in_sam_file, out_sam_file):
    '''
    Filters multimappers from a SAM file and writes the filtered reads to a new
    SAM file.

    This function scans a SAM file recording how many alignments per read there
    are. Then, for the reads with multiple alignments (multimappers), it counts
    the number of matches and indels it contain. Finally, it writes to a new SAM
    file the alignments with the more matches and the less indels. If different 
    alignments obtain the same score, all of them will be kept.

    Args:
        in_sam_file (str): 
            Path to the input SAM file.
        out_sam_file (str): 
            Path to the output SAM file.

    Returns:
        None
    '''

    # Open input SAM file and create output SAM file
    samfile = pysam.AlignmentFile(in_sam_file, "r")
    out_sam = pysam.AlignmentFile(out_sam_file, "w", template=samfile)

    # Create dictionary to group reads by query name
    read_dict = {}
    for read in samfile:
        # Skip secondary and supplementary alignments
        if read.is_secondary or read.is_supplementary:
            continue
        
        if read.query_name in read_dict:
            read_dict[read.query_name].append(read)
        else:
            read_dict[read.query_name] = [read]

    # Iterate over reads with the same query name
    for read_name in read_dict:
        reads = read_dict[read_name]
        if len(reads) == 1:  
            out_sam.write(reads[0])
        else:  
            # Create list of tuples with read object, number of matches and indels
            read_counts = [(r, count_matches(r), count_indels(r)) for r in reads]
            # Find the maximum counts for matches and indels
            max_counts = max(read_counts, key=lambda x: (x[1], -x[2]))[1:]
            # Create list of all reads with the same maximum counts
            best_reads = [mr[0] for mr in read_counts if mr[1:] == max_counts]
            
            if len(best_reads) == 1:
                out_sam.write(best_reads[0])
            else:
                for best_read in best_reads:
                    out_sam.write(best_read)

    # Close input and output files
    samfile.close()
    out_sam.close()



def main():
    filter_multimappers(args.input, args.output)


if __name__ == "__main__":
    args = parse_arguments().parse_args()
    main()