#!/usr/bin/env python

import argparse
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
        type=argparse.FileType('r'),
        default=sys.stdin,
        nargs='?'
        )
    parser.add_argument(
        '-o', '--output', 
        help="Path for the gff3 final file.",
        type=argparse.FileType('w'),
        default=sys.stdout,
        nargs='?'
        )
    parser.add_argument(
        '-e','--extension',
        help="Minimum number of base pairs each overhang must be of. Maximum value = 10. Default value = 10.",
        default=10,
        type=int,
        choices=range(11)
        )
    parser.add_argument(
        '--chr',
        help="Path to the file containing a tabulated file with the chromosome and its length.",
        type=str
        )
        
    return parser


def extend_overhang(pre_mirna:list, mirna_start:int, mirna_end:int, extension:int, chr_size:int):
    '''Extending pre-miRNAs overhang with respect to its canonical forms.

    This function modifies if necessary the start and/or end positions of 
    the pre-miRNA annotation to have  an overhang of minimum the value
    passed in "extension".

    Args:
        pre (list): pre-miRNA annotation entry
        mirna_start (int): smallest position of the mature form(s)
        mirna_end (int): highest position of the mature form(s)
        extension (int): Minimum number of basepairs to have as overhang 
        chr_size (int): Size of the sequences chromosome

    Returns:
        pre (list): pre-miRNA modified entry 
    '''

    start_overhang = int(mirna_start) - int(pre_mirna[3])
    end_overhang = int(pre_mirna[4]) - int(mirna_end)
    
    if start_overhang < extension and int(pre_mirna[3]) > extension - start_overhang:
        pre_mirna[3] = int(pre_mirna[3]) - extension - start_overhang
    elif start_overhang < extension and int(pre_mirna[3]) < extension - start_overhang:
        pre_mirna[3] = 0
    
    if end_overhang < extension and int(pre_mirna[4]) + extension - end_overhang < int(chr_size):
        pre_mirna[4] = int(pre_mirna[4]) + extension - end_overhang
    elif end_overhang < extension and int(pre_mirna[4]) + extension - end_overhang >= chr_size:
        pre_mirna[4] = chr_size
    
    return pre_mirna



########################### 
###     MAIN PROGRAM    ###
###########################

def main():

    '''Extend pre-miRNA overhang main program
     
    This program modifies when necessary the start and/or end positions of 
    the pre-miRNA in the annotation file to have an overhang of a minimum 
    value passed in "extension" in relation to its mature form(s).
    Annotation entries located at the extrems will only be extended the allowed
    number of bases by the chromosome start and end positions.
     
    Args:
        input (file): pre-miRNA annotation file in gff3 format
        extension (int): Minimum number of basepairs to have as overhang 
        chr (file): Tabulated file with chromosome names and their length.
    
    Returns:
        output (file): A gff3 file with all the entries having a minimum overhang 
        of the number specified in the "extension" parameter.
    '''

    ## Creating chromosome size dictionary
    chr_size_dict = {}
    with open(args.chr, "r") as f:
        for line in f:
            (chr, chr_length) = line.strip().split()
            chr_size_dict[chr] = chr_length

    extension =  args.extension

    ## Read gff3 input file
    with args.input as f:
        records = [line.strip().split() for line in f.readlines() if line[0] != "#"]

    ## Setting values for previous and next entries
    prev_rec = None
    pos_rec = None

    ## Iterating over all records except for the lasts 2
    for idx in range(len(records) - 2):
        
        chr_size = chr_size_dict[records[idx][0]]

        # if not first entry and not a primary transcript
        if prev_rec != None and records[idx][2] != "miRNA_primary_transcript":
            # if having a single mature form
            if prev_rec == pos_rec:
                records[idx -1] = extend_overhang(records[idx-1], records[idx][3], records[idx][4], extension, chr_size)
            
            # if having two mature forms
            elif pos_rec == "miRNA":
                if records[idx][3] < records[idx + 1][4]:
                    records[idx -1] = extend_overhang(records[idx -1], records[idx][3], records[idx + 1][4], extension, chr_size)
                elif records[idx][3] > records[idx + 1][4]:
                    records[idx -1] = extend_overhang(records[idx -1], records[idx + 1][3], records[idx][4], extension, chr_size)
        # Updating previous and next records
        prev_rec = records[idx][2]
        pos_rec = records[idx + 2][2]

    ## Checking last two entries
    # if having a single mature form
    if records[-2][2] == "miRNA_primary_transcript":
        records[-2] = extend_overhang(records[-2], records[-1][3], records[-1][4], extension, chr_size)

    # if having two mature forms
    else:
        if records[-2][3] < records[-1][4]:
            records[-3] = extend_overhang(records[-3], records[-2][3], records[-1][4], extension, chr_size)
        elif records[-2][3] > records[-1][4]:
            records[-3] = extend_overhang(records[-3], records[-1][3], records[-2][4], extension, chr_size)

    ## Write final gff3 file
    with args.output  as gff3_file:
        for record in records:
            seqid, source, feature, start, end, score, strand, phase, attributes = record
            gff3_line = f"{seqid}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
            gff3_file.write(gff3_line)


if __name__ == "__main__":

    args = parse_arguments().parse_args()
    main()