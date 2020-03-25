#!/bin/bash

#########################################################
### Paula Iborra, Biozentrum, University of Basel     ###
### paula.iborra@alumni.esci.upf.edu                  ###
### JAN-2020                                          ###
#########################################################

#####################
###  DESCRIPTION  ###
#####################

# Download and process genome sequences fasta files.

####################
###  PARAMETERS  ###
####################

# Prefix for filenames
fileNamePrefix="$1"
organism="$2"   

# Paths (DO NOT CHANGE!)
root="$PWD"
#root="$(cd "$(dirname "$0" )" && pwd)"
resDir="${root}/results/${organism}/${fileNamePrefix}/"
rawDir="${resDir}/raw"
logDir="${root}/logs/local/${organism}/${fileNamePrefix}/"

# URLs
# ------
# All URLs variables represent Bash arrays, so that multiple URLs can be provided; in that case, 
# files are concatenated after download
genomeSeqURLs="$3"   #Modified by Iborra P

########################
###  PRE-REQUISITES  ###
########################

# Shell options
set -e
set -u
set -o pipefail

# Create directories
mkdir --parents "$resDir"
mkdir --parents "$rawDir"
mkdir --parents "$logDir"

# Create log file
logFile="${logDir}/$(basename $0 ".sh").log"
rm -f "$logFile"; touch "$logFile"
>&2 echo "Log written to '$logFile'..."

##############
###  MAIN  ###
##############

## GET & FILTER GENE ANNOTATIONS

# Get genome sequences fasta files
echo "Downloading genome sequences files..." >> "$logFile"
# wget -i "${genomeSeqURLs}" --output-document "${rawDir}/${fileNamePrefix}.genome.fa.gz"
# genomeSeq="${resDir}/${fileNamePrefix}.genome.fa.gz"

for url in "${genomeSeqURLs[@]}"; do
    wget "$url" --output-document "${rawDir}/$(basename "$url")" &> /dev/null
done

# Concatenate genome sequences fasta files
echo "Concatenating genome sequences files..." >> "$logFile"
genomeSeq="${resDir}/genome.fa.gz"
for url in "${genomeSeqURLs[@]}"; do
    cat "${rawDir}/$(basename "$url")" >> "$genomeSeq"
done

# Trim genome sequences IDs
echo "Triming genome sequenes IDs..." >> "$logFile"
genomeSeqTrim="${resDir}/genome.processed.fa"
gunzip "$genomeSeq"
genomeSeq="${resDir}/genome.fa"
awk '{if ($1 ~ /^>/) {print $1} else {print $0}}' "$genomeSeq" > "$genomeSeqTrim"
rm "${resDir}/genome.fa"

#############
###  END  ###
#############

echo "Original data in: $rawDir" >> "$logFile"
echo "Processed data in: $resDir" >> "$logFile"
echo "Done. No errors." >> "$logFile"
>&2 echo "Done. No errors."
