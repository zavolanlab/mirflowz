#!/usr/bin/env python
import sys
import re#!/usr/bin/env python
import sys
import re

### Usage: python2 fastaformat.py [FASTA]

class Seq:
    def __init__(self):
        self.id=""
    def __init__(self):
        self.seq=""
      
record=[]
nrec=-1
inseq=0
with open(sys.argv[1]) as f:
    for line in f:
        if re.match(r'^>' , line):
            nrec+=1
            record.append(Seq())
            mobj=re.match(r'>(\w+)(.*)' , line)

            if (mobj):
                record[nrec].id=mobj.group(1)
            inseq=0
        else :
            if inseq==0 :
                inseq=1
                record[nrec].seq = line   
            else:
                cstring=record[nrec].seq+line
                record[nrec].seq = cstring
           
for x in range (0,nrec+1):
    print ">%s\n%s"%(record[x].id, re.sub(r'\n', "", record[x].seq)) 
    
