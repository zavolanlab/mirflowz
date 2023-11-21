#!/usr/bin/env python

"""Filter alignments in a SAM file by NH tag.

This script uses the pysam library to open the input SAM file and iterates
over each alignment in it. If the NH tag is higher than the provided `max_NH`
value, the aligned read is removed.

Usage: filter_nh.py [SAM file] [max_NH] [OUTPUT file]
"""

import sys
import pysam

if sys.argv[1] in ['--help', '-h', '-help']:
    sys.exit("\nDescription: Checks for NH tag to remove reads that aligned "
             "more than max_NH value.\nUsage: filter_nh.py [SAM file] [max_NH]"
             "[OUTPUT file]\n")
elif len(sys.argv) < 4 or len(sys.argv) > 4:
    sys.exit("\n Arguments ERROR. See [nh_filter.py --help]\n")


def main():
    """Filter alignments by NH tag."""
    sys.stdout.write(
        f"Removing reads aligned more than {sys.argv[2]} times... \n"
    )

    infile = pysam.Samfile(sys.argv[1], "r", check_sq=False)
    out = pysam.Samfile(sys.argv[3], "w", template=infile)

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
