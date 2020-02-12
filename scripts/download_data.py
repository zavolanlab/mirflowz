#!/usr/bin/env python

import sys
import os
import re
import wget
import tempfile
import shutil
from argparse import ArgumentParser, RawTextHelpFormatter



### Created: Mar 26, 2019
### Author: Paula Iborra 
### Company: Zavolan Group, Biozentrum, University of Basel

### ARGUMENTS ###

parser = ArgumentParser(
    description="Script to prepare necessary input files. Get files from local path / download them. "
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '-i','--input', 
    help="Input text file with one file per line (url or local path) to download",
    required=True, 
    )
parser.add_argument(
    '-o','--output', 
    help="Output path where to download files, with folder name to be created. i.e. your_local_path/new_folder_name",
    required=True, 
    )
    

args = parser.parse_args()


def get_data(files, tmpdir, output_folder):
    typ = '' # define file type (local / url)
    nfiles = [] # list of files to be concatenated
    name = '' #name for single file
    out = '' #name for the files concatenated
    for file in files:
        if file.startswith("#"): #ignore comments in the data list file
            pass
        elif file.startswith('*'):
            out = (re.match(r'([^*]*)',file)).group()
        elif file.startswith('>'):
            typ = file
        elif file == '\n':
            if out != '':
                with open(os.path.join(output_folder,out), 'w') as outfile:
                    for fname in filenames:
                        with open(os.path.join(tmpdir,fname)) as infile:
                            outfile.write(infile.read())
                out = ''
            else:
                if len(nfiles) != 0: 
                    with open(os.path.join(outpur_folder,name), 'w') as outfile:
                        for fname in filenames:
                            with open(os.path.join(tmpdir,fname)) as infile:
                                out.write(infile.read())
                else:
                    continue
        else:
            name = (re.match(r'([^\/]*)',file[::-1])).group()[::-1]
            nfiles.append(name)
            if typ == '>DOWNLOAD':
                download(file, tmpdir, name)
            else:
                local_path(file, tmpdir, name)



def download(file, tmpdir, name):
    wget.download(file, out=os.path.join(tmpdir,name))


def local_path(file, tmpdir, name):
    shutil.copy(file, os.path.join(tmpdir,name))

tmpdir = tempfile.mkdtemp()
output_folder = args.output

files=[line.rstrip('\n') for line in open(args.input, 'r')]
with open(args.input, 'r') as data:
    get_data(data, tmpdir, output_folder)

os.removedirs(tmpdir)  
