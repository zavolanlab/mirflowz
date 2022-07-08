#!/usr/bin/env Rscript

# (c) 2019 Paula Iborra, Biozentrum, University of Basel

#################
###  IMPORTS  ###
#################

# Import required packages
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("dplyr"))) == FALSE ) { stop("[ERROR] Package 'dplyr' required! Aborted.") }


#######################
###  PARSE OPTIONS  ###
#######################

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Build description message
description <- "Merge miRNAs quantification tables.\n"
author <- "Author: Paula Iborra, Biozentrum, University of Basel"
version <- "Version: 1.0.0 (JUN-2019)"
requirements <- "Requires: optparse"
msg <- paste(description, author, version, requirements, sep="\n")

# Define list of arguments
option_list <- list(
  make_option(
    "--input_dir",
    action="store",
    type="character",
    default=getwd(),
    help="Absolute path from where input files shall be readed. Required!",
    metavar="directory"
  ),
  make_option(
    "--output_file",
    action="store",
    type="character",
    default=file.path(getwd(), "counts.tab"),
    help="Table output file path. Default: $PWD/counts.tab",
    metavar="directory"
  ),
  make_option(
    c("--prefix"),
    action="store_true",
    type="character",
    default=NULL,
    help="Prefix for reading input files. Default: NULL.",
    metavar="file"
  ),
  make_option(
    c("-h", "--help"),
    action="store_true",
    default=FALSE,
    help="Show this information and die."
  ),
  make_option(
    c("-u", "--usage"),
    action="store_true",
    default=FALSE,
    dest="help",
    help="Show this information and die."
  ),
  make_option(
    c("-v", "--verbose"),
    action="store_true",
    default=FALSE,
    help="Print log messages to STDOUT."
  )
)

# Parse command-line arguments
opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] --input_dir <path/to/input/files>\n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

# Re-assign variables
in.dir <- opt$`input_dir`
prefix <- opt$`prefix`
out.file <- opt$`output_file`
verb <- opt$`verbose`

# Validate required arguments
if ( is.null(in.dir) ) {
  print_help(opt_parser)
  stop("[ERROR] Required argument missing! Aborted.")
}

######################
###    FUNCTIONS   ###
######################

merge_tables <- function(cwd, prefix){
  dataFiles <- dir(cwd, prefix, full.names=TRUE)
  mat <- NULL
  if (length(dataFiles)) {
    mat <- read.table(dataFiles[1], sep='\t')
    sample <- gsub(prefix, "", dataFiles[1])
    colnames(mat)[2] <- basename(sample)
    for (i in seq_len(length(dataFiles)-1)) {
       mat <- full_join(mat, read.table(dataFiles[i+1], sep = "\t"), by='V1')
       sample <- gsub(prefix, "", dataFiles[i+1])
       colnames(mat)[i + 2] <- basename(sample)
    }
    colnames(mat)[1] <- "ID"
  }
  return(mat)
}

######################
###      MAIN      ###
######################
# Write log
if ( verb ) cat("Creating output directory...\n", sep="")

# Create output directories
dir.create(dirname(out.file), recursive=TRUE, showWarnings=FALSE)

# Write log
if ( verb ) cat("Creating table...\n", sep="")

# Create table from input directory files
myTable <- merge_tables(cwd=in.dir, prefix=prefix)

# Write log
if ( verb ) cat(paste("Writing table: ", out.file, "\n", sep=""), sep="")

# Writing table
write.table(myTable, out.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Write log
if ( verb ) cat("Done.\n", sep="")
