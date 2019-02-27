#!/usr/bin/env python

import sys
import re
import HTSeq
from argparse import ArgumentParser
import pandas as pd



### ARGUMENTS ###

parser = ArgumentParser(
    description='Script to filter GTF files'
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help='Show program\'s version number and exit'
    )
parser.add_argument(
    '--idfield', 
    help='Name of field that contains ID of interest. Default: transcript_id',
    default='transcript_id',
    choices=('transcript_id','gene_id', 'exon_id'),
    )
parser.add_argument(
    '--idlist', 
    help='Generate text file with one ID per line'
    )
parser.add_argument(
    '-f','--filter',
    help='Remove IDs from GTF file missing in filter file. Filter file must contain ONE ID per line'
    )
parser.add_argument(
    '-i','--input', 
    required=True, 
    help='Input GTF file', 
    type=str
    )
parser.add_argument(
    '-o','--output',
    required=True,
    help='Output GTF file'
    )
    
args = parser.parse_args()

#list filtered IDs
ids = []

### PARSE GTF FILE ###

sys.stdout.write('Parsing GTF file...')

# parse gtf file
gtf_file = HTSeq.GFF_Reader(args.input)
sys.stdout.write('DONE\n') 

### FILTER GTF ###
# keep everyline not containing the id field, and the ones that contains the id field if id in id_filter list. 

if (args.filter):
    sys.stdout.write('Filtering GTF file - Step 1...')
    
    #open filter file and make it a list
    id_filter=[line.rstrip('\n') for line in open(args.filter)]

    w = open('pregtf.gtf', 'w')

    for gtf_line in gtf_file:
        if gtf_line.type in ['gene', 'transcript', 'exon']:
        ## filter by gene_id
            if args.idfield == 'gene_id':
                if gtf_line.attr['gene_id'] in id_filter:
                    w.write(gtf_line.get_gff_line())

        ## filter by anything else
            else:
                if args.idfield not in gtf_line.attr.keys():
                    w.write(gtf_line.get_gff_line())
                
                else:
                    if gtf_line.attr[args.idfield] in id_filter:
                        w.write(gtf_line.get_gff_line())
    
    sys.stdout.write('DONE\n')

# remove those lines with no childs (empty genes / transcripts)
    
    sys.stdout.write('Filtering GTF file - Step 2...')
    pregtf = open('pregtf.gtf', 'r')
    w = open(args.output, 'w')
    previous = []

    if args.idfield == 'gene_id':
        for line in pregtf:
            w.write(line)

    elif args.idfield == 'transcript_id':
        # store first line
        if pregtf:
            first_line = pregtf.readline().split('\t')
            previous.append(first_line)

        for line in pregtf:
            line=line.split('\t')

            if line[2] == previous[-1][2]:
                if line[2] != 'gene':
                    w.write('\t'.join(previous[-1]))
                    del previous[-1]
                    previous.append(line)
                else:
                    del previous[-1]
                    previous.append(line)
            else:
                w.write('\t'.join(previous[-1]))
                del previous[-1]
                previous.append(line)

        #for last line
        if previous[-1][2] == 'exon':
            w.write('\t'.join(previous[-1]))

    elif args.idfield == 'exon_id':
        for line in reversed(open('pregtf.gtf','r').readlines()):
            line = line.split('\t')
            if line[2] == 'exon':
                previous.append(line)
            else:
                if len(previous) == 0:
                    continue
                else:
                    if line[2] == 'transcript':
                        if previous[-1][2] == 'exon':
                            previous.append(line)
                    elif line[2] == 'gene':
                        if previous[-1][2] == 'transcript':
                            previous.append(line)
                    else:
                        print('FILTER ERROR')
                        break
        for i in previous[::-1]:
            w.write('\t'.join(i))

    sys.stdout.write('DONE\n')


### OUTPUT GTF FILE ### 

# sys.stdout.write('Writing GTF file...')
# w = open(args.output, 'w')
# for gtf_line in gtf_file:
#     if gtf_line.type == 'gene':
#         if (gtf_line.attr['gene_id'] in gene_trx):
#             w.write(gtf_line.get_gff_line())

#     elif gtf_line.type == 'transcript':
#         if (gtf_line.attr['transcript_id'] in trx):
#             w.write(gtf_line.get_gff_line())

#     elif gtf_line.type == 'exon':
#         if (
#             gtf_line.attr['gene_id'] in gene_trx and
#             gtf_line.attr['transcript_id'] in trx
#             ):
#                 w.write(gtf_line.get_gff_line())

#     else:
#         continue

# w.close()
# sys.stdout.write('DONE\n')

### OUTPUT IDs GTF FILE ###

# if (args.idlist):
#     sys.stdout.write('Creating IDs list from GTF file...')
#     with open(args.idlist, 'w') as id_list:
#         id_list.write('\n'.join(trx))
#         id_list.close()
#         sys.stdout.write('DONE\n')


