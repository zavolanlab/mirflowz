#!/bin/bash

# Setting to strict mode
set -euo pipefail
IFS=$'\n\t'

#### FUNCTIONS

usage()
{
    echo "usage: filter_annotation.sh [[[-f file ] [-o outname]] | [-h]]"
}



#### MAIN
# Test whether all required inputs are present
if [[ $1 == -h ]] || [[ $# != 4 ]]
then
	usage
	exit
fi

# Get parameters
while [ $# -gt 0 ]; do
    case $1 in
        -f | --file )           shift
                                filename=$1
                                ;;
        -o | --out)             shift
                                outname=$1
				;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

printf "\nRemoves miRs with ID _1 in \"%s\" and write output to %s.\n" \
	${filename}\
	${outname}

#Remove lines with _1 ID
awk -F ';' '{print $0}' ${filename} | awk -F '=' '{print $0}'| grep -v "_1" > ${outname}

printf "\nDONE!\n"
