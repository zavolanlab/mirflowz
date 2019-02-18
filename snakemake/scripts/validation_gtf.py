#!/usr/bin/env python

import gffutils
import sys
import re
import HTSeq
from argparse import ArgumentParser


### ARGUMENTS ###

parser = ArgumentParser(
    description="Script to filter GTF files"
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '--idfield', 
    help="Name of field that contains ID of interest. Default: transcript_id",
    type=str,
    default="transcript_id",
    choices=["transcript_id","gene_id"]
    )
parser.add_argument(
    '--idlist', 
    help="Generate text file with one ID per line"
    )
parser.add_argument(
    '-f','--filter',
    help="Remove IDs from GTF file missing in filter file. Filter file must contain ONE ID per line"
    )
parser.add_argument(
    '-i','--input', 
    required=True, 
    help="Input GTF file", 
    type=str
    )
parser.add_argument(
    '-o','--output',
    required=True,
    help="Output GTF file"
    )
    
args = parser.parse_args()

# dictionary of genes --> transcripts
gene_trx = dict()

# dictionary of min start of accepted genes ids
genes_dict_start = dict()

# dictionary of max end of accepted gene ids
genes_dict_end = dict()

#list of transcripts
trx = []

### PARSE GTF FILE ###

sys.stdout.write('Parsing GTF file...')

# parse gtf file
gtf_file = HTSeq.GFF_Reader(args.input)

for gtf_line in gtf_file:

    #add gene id to dict
    if gtf_line.type == 'gene':
        gene_id = gtf_line.attr['gene_id']
        gene_trx[gene_id] = []

    #add trx id to its gene in dict 
    if gtf_line.type == 'transcript':
        if gtf_line.attr['gene_id'] in gene_trx.keys():
            trx_id = gtf_line.attr['transcript_id']
            gene_trx[gene_id].append(trx_id)
            #add trx id to a list of all transcripts
            trx.append(trx_id)
    else:
        continue

sys.stdout.write('DONE\n') 

### FILTER GTF ###

if (args.filter):
    sys.stdout.write('Filtering GTF file...')
    
    #open filter file and make it a list
    id_filter=[line.rstrip('\n') for line in open(args.filter)]

    if args.idfield == 'transcript_id':
        for gene_id, trxs in gene_trx.items():
            print(gene_id, trxs)
            for trx_id in trxs:
                print(trx_id)
                if trx_id not in id_filter:
                    gene_trx[gene_id].remove(trx_id)
                    trx.remove(trx_id)

    elif args.infield == 'gene_id':
        for gene_id in gene_trx:
            if gene_id not in id_filter:
                gene_trx.pop(x)

    gene_trx = {k: v for k, v in gene_trx.items() if (len(v)>0)}
    sys.stdout.write('DONE\n') 

else:
    gene_trx = {k: v for k, v in gene_trx.items() if (len(v)>0)}

### OUTPUT GTF FILE ###

sys.stdout.write('Writing GTF file...')
w = open(args.output, 'w')
for gtf_line in gtf_file:
    if gtf_line.type == 'gene':
        if (gtf_line.attr['gene_id'] in gene_trx):
            w.write(gtf_line.get_gff_line())

    elif gtf_line.type == 'transcript':
        if (gtf_line.attr['transcript_id'] in trx):
            w.write(gtf_line.get_gff_line())

    elif gtf_line.type == 'exon':
        if (
            gtf_line.attr['gene_id'] in gene_trx and
            gtf_line.attr['transcript_id'] in trx
            ):

                w.write(gtf_line.get_gff_line())

    else:
        continue

w.close()
sys.stdout.write('DONE\n')

### OUTPUT IDs GTF FILE ###

if (args.idlist):
    sys.stdout.write('Creating IDs list from GTF file...')
    with open(args.idlist, 'w') as id_list:
        id_list.write('\n'.join(trx))
        id_list.close()
        sys.stdout.write('DONE\n')


