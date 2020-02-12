#!/bin/bash

# Setting to strict mode
set -euo pipefail
IFS=$'\n\t'

#### FUNCTIONS

usage()
{
    echo "usage: split_file_by_colum_value.sh [[[-f file] [-c column] [-s separator]] | [-h]]"
}



#### MAIN
# Test whether all required inputs are present
if [[ $1 == -h ]] || [[ $# != 6 ]]
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
        -c | --column )   	shift
				column=$1
                                ;;
        -s | --separator )   	shift
				separator=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

printf "\nFile splited by colum %i values.\n" \
	${column}\

# Split file by values in requested column
awk -v col="${column}" -v sep="${separator}" -F'$sep' '{print > $col".txt"}' ${filename}

printf "\nDONE!\n"
