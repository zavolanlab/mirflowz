#!/usr/bin/env python
import sys
import re
import gzip
from argparse import ArgumentParser, RawTextHelpFormatter


### Created: Mar 5, 2019
### Author: Paula Iborra 
### Company: Zavolan Group, Biozentrum, University of Basel


### ARGUMENTS ###

parser = ArgumentParser(
    description="Script to filter FASTA files"
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '--trim', 
    help="Character\'s used to trim the ID. Remove anything that follows the character/s. Write \\ infront of \'.\' and \'-\'(i.e trim=\"$\\.\\-|_\").  Default: first white space",
    type=str, 
    nargs='?', 
    default=""
    )
parser.add_argument(
    '--idlist', 
    help="Generate text file with one ID per line"
    )
parser.add_argument(
    '-f','--filter',
    help="Remove IDs from FASTA file missing in filter file. Filter file must contain ONE ID per line"
    )
parser.add_argument(
    '-i','--input', 
    required=True, 
    help="Input FASTA file", 
    type=str
    )
parser.add_argument(
    '-o','--output', 
    help="Output FASTA file"
    )
    
args = parser.parse_args()

### PARSE FASTA FILE ###

class Seq:
    def __init__(self):
        self.id=""
    def __init__(self):
        self.seq=""
    def __init__(self):
        self.features=""

     
record=[] #list of records 
nrec=-1
inseq=0


# open files 
if args.input.endswith('.gz'):
    f = gzip.open(args.input, 'rt')
else:
    f = open(args.input)

# parse fasta file
sys.stdout.write("Parsing FASTA file...")
for line in f:
    if re.match(r'^>' , line):
        nrec+=1
        record.append(Seq())

        # define id of the record
        if not args.trim:
            mobj=re.match ( r'^>(\S*)(.*)', line) 
        else:
            mobj=re.match(r'^>([^%s]*)(.*)'%args.trim , line)

        # add id and features
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
sys.stdout.write("DONE\n")

## ID FILTER LIST ##

if (args.filter):
    sys.stdout.write("Filtering FASTA file...")
    id_filter=[line.rstrip('\n') for line in open(args.filter)]
    sys.stdout.write("DONE\n")


## OUTPUT FASTA FILE ##

if (args.output):
    sys.stdout.write("Writing FASTA file...")
    with open(args.output, 'w') as output:
        if (args.filter):
            for x in range(0,nrec+1):
                if record[x].id in id_filter:
                    output.write(">%s\n%s"%(record[x].id, record[x].seq))
        else:
            for x in range(0,nrec+1):
                output.write(">%s\n%s"%(record[x].id, record[x].seq))
    output.close()
    sys.stdout.write("DONE\n")

## OUPUT LIST IDs ##    

idlist=[]
if (args.idlist):
    sys.stdout.write("Creating IDs list from FASTA file...")
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
    sys.stdout.write("DONE\n")


