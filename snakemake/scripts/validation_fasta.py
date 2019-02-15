#!/usr/bin/env python
import sys
import re
import argparse
import gzip

### Created: Feb 12, 2019
### Author: Paula Iborra 
### Company: Zavolan Group, Biozentrum, University of Basel

### TEST FILES ###
'''
out="/scicore/home/zavolan/devagy74/projects/mir-prepare-annotation/snakemake/results/fasta_out.fa"
fasta="/scicore/home/zavolan/devagy74/projects/mir-prepare-annotation/snakemake/fastatest.fa"
filter="/scicore/home/zavolan/devagy74/projects/mir-prepare-annotation/snakemake/filterlist.txt"
idlist="/scicore/home/zavolan/devagy74/projects/mir-prepare-annotation/snakemake/results/idlist.txt"
fastagz="/scicore/home/zavolan/devagy74/projects/mir-prepare-annotation/snakemake/fastatest.fa.gz"

i.e usage:
python ~devagy74/projects/mir-prepare-annotation/snakemake/scripts/validation_fasta.py --trim="." --idlist $idlist -f $filter -i $fasta -o $out
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Script to filter FASTA files"
    )
    parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0',
        help="Show program's version number and exit")
    parser.add_argument('--trim', 
        help="Character/s used to trim the ID. Remove anything that follows the character/s. \
        Default: first white space",
        type=str, nargs='?', default="")
    parser.add_argument('--idlist', help="Generate text file with one ID per line")
    parser.add_argument('-f','--filter', help="Remove IDs from FASTA file missing in filter file.\
        Filter file must contain ONE ID per line")
    parser.add_argument('-i','--infile', required=True, help="input FASTA file", type=str)
    parser.add_argument('-o','--output', help="Output FASTA file")
    args = parser.parse_args()

### PARSE FASTA FILE ###

class Seq:
    def __init__(self):
        self.id=""
    def __init__(self):
        self.seq=""
    def __init__(self):
        self.features=""
      
record=[]
nrec=-1
inseq=0

if args.infile.endswith('.gz'):
    f = gzip.open(args.infile, 'rt')
else:
    f = open(args.infile)

for line in f:
    if re.match(r'^>' , line):
        nrec+=1
        record.append(Seq())
        if args.trim == None:
            mobj=re.match ( r'^>(\S*)(.*)', line)
        elif args.trim == ".":
            mobj=re.match(r'^>([^\.]*)(.*)' , line)  
        else:
            mobj=re.match(r'^>([^%s]*)(.*)'%args.trim , line)

        if (mobj):
            record[nrec].id=mobj.group(1)
            record[nrec].features=mobj.group(2)
        inseq=0
    else :
        if inseq==0 :
            inseq=1
            record[nrec].seq = line   
        else:
            cstring=record[nrec].seq+line
            record[nrec].seq = cstring


## ID FILTER LIST ##

if (args.filter):
    id_filter=[line.rstrip('\n') for line in open(args.filter)]

## OUTPUT FASTA FILE ##

if (args.output):
    with open(args.output, 'w') as output:
        if (args.filter):
            for x in range(0,nrec+1):
                if record[x].id in id_filter:
                    output.write(">%s\n%s\n"%(re.sub(r'\n', "",record[x].id), re.sub(r'\n', "", record[x].seq)))
                else:
                    continue
        else:
            for x in range(0,nrec+1):
                output.write(">%s\n%s\n"%(re.sub(r'\n', "",record[x].id), re.sub(r'\n', "", record[x].seq)))
    output.close()

## OUPUT LIST IDs ##    

idlist=[]
if (args.idlist):
    with open(args.idlist, 'w') as id_list:
        if (args.filter):
            for x in range(0,nrec+1):
                if record[x].id in id_filter:
                    idlist.append(record[x].id)
                else:
                    continue
        else:
            for x in range(0,nrec+1):
                idlist.append(re.sub(r'\n', "",record[x].id))
        idlist.sort()
        id_list.write('\n'.join(idlist))
    id_list.close()

