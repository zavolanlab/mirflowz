#!/bin/bash

#########################################################
### Alexander Kanitz, Biozentrum, University of Basel ###
### alexander.kanitz@unibas.ch                        ###
### 27-APR-2016                                       ###
#########################################################


#####################
###  DESCRIPTION  ###
#####################

# Obtains and filters genome, gene annotations and transcriptome.


####################
###  PARAMETERS  ###
####################

output_dir="$1"
log_dir="$2" 

# Paths (DO NOT CHANGE!)  #Modified by Iborra P
root="$PWD"
resDir="${root}/${output_dir}"
fileDir="${root}/test_files"
tmpDir="${root}/.tmp"
logDir="${root}/${log_dir}"


# Annotation file 
# ----
# - It is assumed that the specified transcriptome files contain sequences for all transcripts in 
# the (filtered) gene annotations
geneAnnoFile="$3"  

# Filters
# -------
# - All filters are positive filters, i.e. entries meeting the specified filters are kept
# - Separate multiple entries by a single space and quote the whole string
# - Set to empty string '""' if no filtering is desired
# - Transcriptome sequences are filtered according to the transcript annotations remaining after 
# applying gene annotation filters
# - To avoid unexpected results, ensure that the chromosome filter for gene annotations is equal to 
# or a subset of the one for the genome
# - Warnings are issued if sequences for annotated transcripts are absent in the transcriptome
# Gene annotation / transcriptome filters
geneAnnoFilterChromosomes=""                    # e.g. "1 2 3 X"
geneAnnoFilterGeneBiotypes=""                   # e.g. "protein_coding lincRNA"
geneAnnoFilterTranscriptBiotypes=""             # e.g. "protein_coding processed_transcript"
geneAnnoFilterTranscriptSupportLevels="1 2 3"   # e.g. "1 2 3"


########################
###  PRE-REQUISITES  ###
########################

# Shell options
set -e
set -u
set -o pipefail

# Create directories
mkdir --parents "$resDir"
mkdir --parents "$tmpDir"

# Create log file
logFile="${logDir}"
rm -fr "$logFile"; touch "$logFile"
>&2 echo "Log written to '$logFile'..."


##############
###  MAIN  ###
##############

## FILTER GENE ANNOTATIONS

# Filter gene annotations
geneAnnoFilt="${resDir}/gene_annotations.filtered.gtf.gz"
geneAnnoOut="${resDir}/gene_annotations.filtered.gtf"
geneAnnoFiltTmp="${tmpDir}/gene_annotations.filtered.gtf.gz.tmp"
geneAnno="${resDir}/gene_annotations.gtf.gz"
cp "$geneAnnoFile" "$geneAnno"
cp "$geneAnno" "$geneAnnoFiltTmp"

    # Filter requested chromosomes
    # ----------------------------
    # - If filter provided, filters comments and matching chromosomes
    if [ "$geneAnnoFilterChromosomes" != ""  ]; then
        echo "Filtering gene annotations by chromosomes..." >> "$logFile"
        perl -ane 'if(!@ARGV){if(/^#\!/){print}else{$keep=$chr{$F[0]}}}$keep?print:chomp;$chr{$_}=1 if @ARGV' <(echo "$geneAnnoFilterChromosomes" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Filter requested gene biotypes
    # ------------------------------
    # - If filter provided, filters comments and matching gene biotypes
    if [ "$geneAnnoFilterGeneBiotypes" != ""  ]; then
        echo "Filtering gene annotations by gene biotypes..." >> "$logFile"
        perl -ne 'if(/^#\!/){print;$keep=0}elsif(/gene_biotype\s\"(\S+)\"/){$keep=$type{$1}}else{$keep=0}$keep?print:chomp;$type{$_}=1 if @ARGV' <(echo "$geneAnnoFilterGeneBiotypes" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Filter requested transcript biotypes
    # ------------------------------------
    # - If filter provided, filters 'gene' entries, commentss and matching transcript biotypes
    if [ "$geneAnnoFilterTranscriptBiotypes" != ""  ]; then
        echo "Filtering annotations by transcript biotypes..." >> "$logFile"
        perl -ane 'if(/^#\!/||$F[2] eq "gene"){print;$keep=0}elsif(/transcript_biotype\s\"(\S+)\"/){$keep=$type{$1}}else{$keep=0}$keep?print:chomp;$type{$_}=1 if @ARGV' <(echo "$geneAnnoFilterTranscriptBiotypes" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Filter requested transcript support levels
    # ------------------------------------------
    # - If filter provided, filters 'gene' entries, comments and matching transcript support levels
    if [ "$geneAnnoFilterTranscriptSupportLevels" != ""  ]; then
        echo "Filtering annotations by transcript support levels..." >> "$logFile"
        perl -ane 'if(/^#\!/||$F[2] eq "gene"){print;$keep=0}elsif(/transcript_support_level\s\"(\S+?)\"?/){$keep=$level{$1}}else{$keep=0}$keep?print:chomp;$level{$_}=1 if @ARGV' <(echo "$geneAnnoFilterTranscriptSupportLevels" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Remove orphan 'genes' (i.e. 'genes' with all child entries removed) & temporary file
    echo "Removing 'orphan' genes..." >> "$logFile"
    perl -ane 'if ($F[2] eq "gene"){$prev=$_}else{print $prev,$_; $prev=""}' <(zcat $geneAnnoFiltTmp) > "$geneAnnoOut"
    rm "$geneAnnoFiltTmp"

rm "${resDir}/gene_annotations.filtered.gtf.gz"
rm "${resDir}/gene_annotations.gtf.gz"

#############
###  END  ###
#############

echo "Original data in: $fileDir" >> "$logFile"
echo "Processed data in: $resDir" >> "$logFile"
echo "Done. No errors." >> "$logFile"
>&2 echo "Done. No errors."