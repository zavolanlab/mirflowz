#!/bin/bash

# define function to print usage
display_usage() {
    echo "Sort oligomap aligments based on their numerical names."
    echo ""
    echo "Usage: $0 input_file number_of_threads output_file"
    echo ""
    echo "Args:"
    echo "    input_file: oligomap aligments"
    echo "    number_of_threads: number of threads to run the sorting with"
    echo "    output_file: path to sorted output file"
}

# show usage if user supplied less than two arguments
if [  $# -ne 3 ]
then
    display_usage
    exit 1
fi

# show usage if user has supplied -h or --help
if [[ ( $# == "--help") ||  $# == "-h" ]] 
then
    display_usage
    exit 0
fi

cat $1 | awk -v RS="" '{ gsub("\n", "*"); print }' | sort -n --parallel=$2 | awk -v ORS="\n\n" '{ gsub("\*", "\n"); print }' > $3
