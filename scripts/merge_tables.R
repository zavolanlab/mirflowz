#!/usr/bin/env Rscript

# (c) 2019 Paula Iborra, Biozentrum, University of Basel

#################
###  IMPORTS  ###
#################

# Import required packages
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) {
  stop("[ERROR] Package 'optparse' required! Aborted.")
}
if (suppressWarnings(suppressPackageStartupMessages(require("dplyr"))) == FALSE) {
  stop("[ERROR] Package 'dplyr' required! Aborted.")
}


#######################
###  PARSE OPTIONS  ###
#######################

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Build description message
description <- "Merge miRNAs quantification tables.\n"
author <- "Author: Paula Iborra, Biozentrum, University of Basel"
mantainer <- "Refactor and documentation: Iris Mestres"
version <- "Version: 1.1.0 (FEB-2024)"
requirements <- "Requires: optparse, dplyr"
msg <- paste(description, author, mantainer, version, requirements, sep="\n")

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
opt_parser <-
  OptionParser(
    usage = paste(
      "Usage:",
      script,
      "[OPTIONS] --input_dir <path/to/input/files>\n",
      sep = " "
    ),
    option_list = option_list,
    add_help_option = FALSE,
    description = msg
  )
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

#' Read and process input table
#'
#' `get_table()` uses `tryCatch()` to read the file in `tbl_pth`. If the table
#' is empty and an error is raised, a data frame is created.
#'
#' @param tbl_pth Path to the input table.
#' @param prefix String to be removed from the input file name. It must be
#'   present in all the tables to be merged.
#'
#' @returns `get_table()` returns a data frame containing the miRNA species to
#'   be counted in first column, named `ID`, and their counts in that file in
#'   the second one. The name of the second column in the data frame is obtained
#'   by removing the `prefix` from the input file name. If no `prefix` is given,
#'   the whole file name is used.
#'
#'   If the input file is empty, the returned data frame will consist on one row
#'   with a `NA` in both fields.
#'
#' @seealso [tryCatch()] which this function uses.
get_table <- function(tbl_pth, prefix) {
  sample <- gsub(prefix, "", tbl_pth)
  fields <- c("ID", basename(sample))
  
  tryCatch(
    expr = {
      table <- read.table(tbl_pth, sep = '\t', col.names = fields)
      return(table)
    },
    error = function(e) {
      table <- data.frame(matrix(NA, ncol = 2, nrow = 1))
      colnames(table) <- fields
      return(table)
    }
  )
}

#' Merge tables with the same prefix
#'
#' `merge_tables()` takes all the files in `cwd` that start with `prefix` and
#' merge them keeping all the miRNA species present in each of the tables.
#'
#' @details The function `get_table()` is used to make sure that even if an
#'   empty input file is given, the merge can still be done. Thus, before
#'   returning the merged table, the row with a `NA` in the `ID` field, if any,
#'   is removed.
#'
#'   The function `dplyr::full_join()` method is used for the merge, therefore,
#'   if a miRNA species in `ID` is missing in any of the tables being joined,
#'   its value is set to `NA` in that column.
#'
#' @param cwd Path to the input tables directory.
#' @param prefix String used in all the tables to be selected for the merge. If
#'    not provided, all the files in `cwd` are used.
#'
#' @returns `merge_tables()` returns a single data frame, `mat`, with all the
#'   miRNA species present in the input tables in the first column, `ID`, and
#'   their counts. Each input file has it own column.
#'
#'   If all the input tables are empty, the output consists only on the table's
#'   header, and if no files starting with `prefix` are found, nothing is
#'   returned.
#'
#' @seealso [get_table()], [dplyr::full_join()] which this function uses.
merge_tables <- function(cwd, prefix) {
  dataFiles <- dir(cwd, prefix, full.names = TRUE)
  mat <- NULL
  
  if (length(dataFiles)) {
    mat <- get_table(dataFiles[1], prefix)
    
    for (i in seq_len(length(dataFiles) - 1)) {
      mat <- full_join(mat, get_table(dataFiles[i + 1], prefix), by = "ID")
    }
    mat <- filter(mat, !is.na(ID))
  }
  return(mat)
}

######################
###      MAIN      ###
######################
# Write log
if (verb)
  cat("Creating output directory...\n", sep = "")

# Create output directories
dir.create(dirname(out.file),
           recursive = TRUE,
           showWarnings = FALSE)

# Write log
if (verb)
  cat("Creating table...\n", sep = "")

# Create table from input directory files
myTable <- merge_tables(cwd = in.dir, prefix = prefix)

# Write log
if (verb)
  cat(paste("Writing table: ", out.file, "\n", sep = ""), sep = "")

# Writing table
write.table(
  myTable,
  out.file,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

# Write log
if (verb)
  cat("Done.\n", sep = "")
