#!/usr/bin/env python

import sys
import argparse
import getopt
import pysam
import os

if sys.argv[1] in ['--help', '-h', '-help']:
    sys.exit("\nDescription: Checks for NH tag to remove reads that aligned more than max_NH value.\nUsage: filter_nh.py [SAM file] [max_NH] [OUTPUT file]\n")
elif len(sys.argv) < 4 or len(sys.argv) > 4:
    sys.exit("\n Arguments ERROR. See [nh_filter.py --help]\n")

def main():
    
    sys.stdout.write("Removing reads aligned more than %s times... \n"%(sys.argv[2]))

    infile = pysam.Samfile(sys.argv[1], "r", check_sq=False)
    out = pysam.Samfile(sys.argv[3] , "w", template = infile )
    
    keep = True

    for DNAread in infile.fetch():
        intags = DNAread.tags
         
        for entry in intags:
            if 'NH' in entry and entry[1] > int(sys.argv[2]):
                keep = False          
        if keep:
            out.write(DNAread)

        keep = True

    out.close()
    sys.stdout.write("DONE!\n")


if __name__ == '__main__':
    main()
