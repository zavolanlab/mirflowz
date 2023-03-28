#!/usr/bin/env python

### EXTEND pre-miRNA OVERHANG ###
# 
# This program modifies if necessary the start and/or end positions of 
# the pre-miRNA in the annotation file to have an overhang of a minimum 
# value passed in "extension" in relation to its mature form(s).
# Annotation entries located at the extrems will only be extended the allowed
# number of bases by the chromosome start and end positions.
# 
# INPUT
# ---------
# i: pre-miRNA annotation file in gff3 format
# extension: Minimum number of basepairs to have as overhang 
# chr: Chromosome size
#
# OUTPUT
# ---------
# o: A gff3 file with all the entries having a minimum overhang of the number
#   specified in the "extension" parameter. 
#

from argparse import ArgumentParser


### ARGUMENTS ###

parser = ArgumentParser(
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
    type=str
    )
parser.add_argument(
    '-o', '--output', 
    help="Path for the gff3 final file.",
    type=str
    )
parser.add_argument(
    '-e','--extension',
    help="Minimum number of base pairs each overhang must be of. Maximum value = 10. Default value = 10.",
    default=10,
    type=int,
    choices=range(0,11)
    )
parser.add_argument(
    '--chr',
    help="Path to the file containing the chromosome and its size.",
    type=str
    )
    
args = parser.parse_args()

with open(args.chr, "r") as f:
    chr_size = int(f.readline().split()[1])

extension =  args.extension


def extend_overhang(pre, mirna_start, mirna_end, extension, chr_size):
    '''
    EXTEND OVERHANG FUNCTION

    This function modifies if necessary the start and/or end positions of 
    the pre-miRNA annotation to have  an overhang of minimum the value
    passed in "extension".

    INPUT
    ---------
    pre: pre-miRNA annotation entry
    mirna_start: smallest position of the mature form(s)
    mirna_end: highest position of the mature form(s)
    extension: Minimum number of basepairs to have as overhang 
    chr_size: Size of the sequences chromosome

    OUTPUT
    ---------
    pre: pre-miRNA modified entry 
    '''

    over_1 = int(mirna_start) - int(pre[3])
    over_2 = int(pre[4]) - int(mirna_end)
    
    if over_1 < extension and int(pre[3]) > extension - over_1:
        pre[3] = int(pre[3]) - extension - over_1
    elif over_1 < extension and int(pre[3]) < extension - over_1:
        pre[3] = 0
    
    if over_2 < extension and int(pre[4]) + extension - over_2 < chr_size:
        pre[4] = int(pre[4]) + extension - over_2
    elif over_2 < extension and int(pre[4]) + extension - over_2 >= chr_size:
        pre[4] = chr_size
    
    return pre



########################### 
###     MAIN PROGRAM    ###
###########################


## Read gff3 input file
with open(args.input, "r") as f:
    records = [line.strip().split() for line in f.readlines() if line[0] != "#"]

## Setting values for previous and next entries
prev_rec = None
pos_rec = None

## Iterating over all records except for the lasts 2
for idx in range(len(records) - 2):

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
with open(args.output, 'w') as gff3_file:
    for record in records:
        seqid, source, feature, start, end, score, strand, phase, attributes = record
        gff3_line = f"{seqid}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
        gff3_file.write(gff3_line)
