#!/usr/bin/env python
import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter
import pysam

### Created: Apr 23, 2019
### Author: Paula Iborra 
### Company: Zavolan Group, Biozentrum, University of Basel


### ARGUMENTS ###

parser = ArgumentParser(
    description="Script to quantify miRNA aligments. Default output: one single quantification table."
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '-u','--uniq',
    type=str,
    help="Get the counting table of specific miR type. Name output: prefix + mirna type",
    )
parser.add_argument(
    '-s','--split',
    action='store_true',
    help="Split the counting tables based on miR type. Name outputs: prefix + mirna type",
    )
parser.add_argument(
    '-c','--column',
    help="Column of miR type in bed file.  Default: 8",
    type=int, 
    default='8'
    )
parser.add_argument(
    '--rid',
    help="Read ID column in intersection bed file. Default: 14",
    type=int,
    default='14'
    )
parser.add_argument(
    '--mid',
    help="miRNA ID;alias;name column in intersection bed file. Default: 10 ",
    type=int,
    default='10'
    )
parser.add_argument(
    '--nh',
    help="read count column in intersection bed file. Default: 15 ",
    type=int,
    default='15'
    )
parser.add_argument(
    '-i','--intersection', 
    required=True, 
    help="miRNA intersection bed file" 
    )
parser.add_argument(
    '-p','--prefix',
    type=str,
    default='',
    help="Prefix output name of quantification table."
    )
    
args = parser.parse_args()



def count_dictionary(file):
    total_counts = {}
    u = 0
    
    sys.stderr.write("##### CREATING COUNTING DICTIONARY #####\n")

    for line in file:
        line = line.split('\t')
        name = (line[(args.mid)-1].split(';'))[2].split('=')[1]
        mirnaid = (line[(args.mid)-1].split(';'))[1].split('=')[1]
        readid = line[(args.rid)-1]
        s = line[(args.column)-1]
                
        if (args.uniq): #create the dictionary just for the specific mirna type
            if s == (args.uniq):
                sys.stderr.write("%s\t%s\t%s\t%s\n"%(readid,name,mirnaid,s))
                if u == 0:
                    total_counts[s] = {}
                    u += 1
                if name not in total_counts[s].keys():
                    if readid != '.':
                        total_counts[s][name] = float(1/int(line[(args.nh)-1]))
                    else:
                        total_counts[s][name] = 0
                else:
                    if readid != '.':
                        total_counts[s][name] += float(1/int(line[(args.nh)-1]))
                    else:
                        total_counts[s][name] += 0

        else: # create the dictionary for all mirna types
            sys.stderr.write("%s\t%s\t%s\t%s\n"%(readid,name,mirnaid,s))
            if s not in total_counts.keys(): # create the dictionary for a new mirna type
                total_counts[s] = {}
                if name not in total_counts[s].keys():
                    if readid != '.':
                        total_counts[s][name] = float(1/int(line[(args.nh)-1]))
                    else:
                        total_counts[s][name] = 0
                else:
                    if readid != '.':
                        total_counts[s][name] += float(1/int(line[(args.nh)-1]))
                    else:
                        total_counts[s][name] += 0

            else: # when mirna type is already in the total_count dictionary
                if name not in total_counts[s].keys():
                    if readid != '.':
                        total_counts[s][name] = float(1/int(line[(args.nh)-1]))
                    else:
                        total_counts[s][name] = 0
                else:
                    if readid != '.':
                        total_counts[s][name] += float(1/int(line[(args.nh)-1]))
                    else:
                        total_counts[s][name] += 0

    return(total_counts)

sys.stderr.write("##### READING INTERSECTION FILE #####\n")
f = open(args.intersection, 'r')
counts = count_dictionary(f)

### UNIQ. TO CREATE JUST ONE SPECIFIC MIRNA TYPE COUNT TABLE

sys.stderr.write("##### WRITING COUNTING TABLES #####\n")

if (args.uniq):
    i = (args.uniq)
    out = open(args.prefix, 'w')
    for x in counts[i].keys():
        cnt = counts[i][x]
        if cnt > 0:
            out.write(x + '\t' + str(cnt) + '\n') 

### SPLIT. CREATE A SEPARE COUNT TABLES FILES FOR EACH MIRNA TYPE
if (args.split):
    for i in counts.keys():
        out = open(args.prefix+i, 'w')
        for x in counts[i].keys():
            cnt = counts[i][x]
            if cnt > 0:
                out.write(x + '\t' + str(cnt) + '\n')

### CREATE ONE TABLE WITH ALL MIRNA TYPE COUNTS
if not (args.uniq) and not (args.split):
    out = open(args.prefix, 'w')
    for i in counts.keys():
        for x in counts[i].keys():
            cnt = counts[i][x]
            if cnt > 0:
                out.write(x + '\t' + cnt + '\n')

sys.stderr.write("##### DONE! #####\n")
