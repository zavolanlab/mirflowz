#!/bin/bash

#########################################################
### Paula Iborra, Biozentrum, University of Basel     ###
### paula.iborra@alumni.esci.upf.edu                  ###
### JAN-2020                                          ###
#########################################################

#####################
###  DESCRIPTION  ###
#####################

# Process genome sequences fasta file.

####################
###  PARAMETERS  ###
####################

# Prefix for filenames
output_dir="$1"
log_dir="$2"   

# Paths (DO NOT CHANGE!)
root="$PWD"
fileDir="${root}/test_files"
resDir="${root}/${output_dir}"
logDir="${root}/${log_dir}"

# Genome File
genomeSeqFile="$3"  

########################
###  PRE-REQUISITES  ###
########################

# Shell options
set -e
set -u
set -o pipefail

# Create directories
mkdir --parents "$resDir"
mkdir --parents "$logDir"

# Create log file
logFile="${logDir}"
rm -fr "$logFile"; touch "$logFile"
>&2 echo "Log written to '$logFile'..."

##############
###  MAIN  ###
##############

## TRIM GENOME SEQUENCES IDs
echo "Triming genome sequenes IDs..." >> "$logFile"

genomeSeqTrim="${resDir}/genome.processed.fa"
cp "${genomeSeqFile}" "${output_dir}/genome.fa.gz"
genomeSeq="${output_dir}/genome.fa.gz"
gunzip "$genomeSeq"
genomeSeq="${output_dir}/genome.fa"
awk '{if ($1 ~ /^>/) {print $1} else {print $0}}' "$genomeSeq" > "$genomeSeqTrim"
rm "${output_dir}/genome.fa"

#############
###  END  ###
#############

echo "Original data in: $fileDir" >> "$logFile"
echo "Processed data in: $resDir" >> "$logFile"
echo "Done. No errors." >> "$logFile"
>&2 echo "Done. No errors."