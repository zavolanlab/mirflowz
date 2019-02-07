#!/usr/bin/env python
import sys
import re#!/usr/bin/env python
import sys
import re

class Seq:
    def __init__(self):
        self.id=""
    def __init__(self):
        self.seq=""
    def __init__(self):
        self.features=""
        

#read lines one at a time
#do not store the whole file in memory if not needed
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
                record[nrec].features=mobj.group(2) 

            inseq=0
        else :
          #  line = re.sub(r'\n', "", line)
            if inseq==0 :
                inseq=1
                record[nrec].seq = line   

            else:
                cstring=record[nrec].seq+line
                record[nrec].seq = cstring
                
           
for x in range (0,nrec+1):
    
    print ">%s\n%s"%(record[x].id, re.sub(r'\n', "", record[x].seq)) 
    #substitution more efficient here
    #this is because the #range < #lines, so we are entering less times the functions. 
