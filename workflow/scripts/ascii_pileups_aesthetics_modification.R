#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Enhance ASCII-style alignment pileups
#
# Filter each ASCII-style alignment pileup in the provided input directory
# retaining those alignments with at least  `--min-count` amount of reads that
# appear in the first `--max-sequence` lines. If `--canonical` is set, mark the
# read(s) matching the canonical sequence(s) inferred from the reference row in
# the output.
#
# Optionally, if `--split-arms` is set, split the precursor pileup into
# arm-specific pileup(s) with an allowed overhang of &pm; `--overhang`
# nucleotides on either side of the mature arm, and the genomic coordinates are
# adjusted to the final representation.
#
# If `--keep-all` is set, the ASCII-style alignment pileup is written even if
# it has no aligned sequences.
#
# The expected input is a tab-separated pileup files with two columns without
# header:
#   1. Alignment string or reference sequence
#   2. Feature name, genomic coordinates, or read count
#
# Expected pileup formats:
#   1. One arm with precursor
#   2. Two arms with precursor
#
# (c) 2026 Zavolan Lab, Biozentrum, University of Basel
# ------------------------------------------------------------------------------


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD REQUIRED LIBRARIES <---#
if (
  suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE
) {
  stop("[ERROR] Package 'optparse' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("dplyr"))) == FALSE
) {
  stop("[ERROR] Package 'dplyr' required! Aborted.")
}

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly = FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Enhance raw ASCII-style alignment pileups.\n"
author <- "Author: Iris Mestres-Pascual <zavolab-biozentrum@unibas.ch>"
maintainer <- "Maintainer: Iris Mestres-Pascual <zavolab-biozentrum@unibas.ch>"
version <- "Version: 1.0.0 (ABR-2026)"
requirements <- "Requires: dplyr, optparse"
msg <- paste(description, author, maintainer, version, requirements, sep = "\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option.list <- list(
  make_option(
    "--in-dir",
    action = "store",
    type = "character",
    default = NULL,
    help = "Absolute path from where input files shall be read. Required!",
    metavar = "dir"
  ),
  make_option(
    "--out-dir",
    action = "store",
    type = "character",
    default = getwd(),
    help = "Absolute path to where output files shall be written.
    [default \"%default\"]",
    metavar = "dir"
  ),
  make_option(
    "--split-arms",
    action = "store_true",
    default = FALSE,
    help = "Split precursor pileups into one mature-arm pileup per arm. All
      subsequent filtering and formatting steps are then applied independently
      to each arm."
  ),
  make_option(
    "--min-count",
    action = "store",
    type = "numeric",
    default = 1,
    help = "Minimum count for a sequence to be kept. [default %default]",
    metavar = "int"
  ),
  make_option(
    "--max-sequences",
    action = "store",
    type = "numeric",
    default = 30,
    help = "Maximum number of top sequences to be displayed. It is assumed that
      the input ASCII-style alignment pileups are already sorted in the desired
      order. [default %default]",
    metavar = "int"
  ),
  make_option(
    "--overhang",
    action = "store",
    type = "numeric",
    default = NULL,
    help = "If `--split-arms` is set, number of extra positions to retain on
      each side of the mature arm span. Reads extending beyond this interval are
      removed. If omitted, only reads fully contained within the exact arm span
      are kept. [default %default]",
    metavar = "int"
  ),
  make_option(
    "--canonical",
    action = "store_true",
    default = FALSE,
    help = "Mark the aligned read(s) corresponding to the canonical mature
      sequence(s)."
  ),
  make_option(
    "--keep-all",
    action = "store_true",
    default = FALSE,
    help = "Write the ASCII-style alignment pileup even if it has no aligned
      sequences."
  ),
  make_option(
      "--prefix",
      action = "store",
      type = "character",
      default = NULL,
      help = "Prefix to be used in the output file name(s). Required!",
      metavar = "string"
  ),
  make_option(
    c("-h", "--help"),
    action = "store_true",
    default = FALSE,
    help = "Show this information and die."
  ),
  make_option(
    c("-u", "--usage"),
    action = "store_true",
    default = FALSE,
    dest = "help",
    help = "Show usage information and die."
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print log messages to STDOUT."
  )
)

## Parse options
opt.parser <-
  OptionParser(
    usage = paste(
      "Usage:",
      script,
      "--in-dir <path/to/input/pileups>",
      "--prefix=<prefix>",
      "[OPTIONS]\n",
      sep = " "
    ),
    option_list = option.list,
    add_help_option = FALSE,
    description = msg
  )
opt <- parse_args(opt.parser)

# Parse command-line arguments
in.dir <- opt$`in-dir`
out.dir <- opt$`out-dir`
keep.all <- opt$`keep-all`
overhang <- opt$`overhang`
canonical <- opt$`canonical`
min.count <- opt$`min-count`
split.arms <- opt$`split-arms`
max.seq <- opt$`max-sequences`
prefix <- opt$`prefix`
verb <- opt$`verbose`

# Dies if required arguments are missing
if (is.null( in.dir ) | is.null( prefix )) {
  print_help( opt.parser )
  stop("[ERROR] argument missing! Aborted.")
}

# Ensure CLI arguments correction
if ( !is.null(overhang) && !split.arms ) {
  print_help( opt.parser )
  stop(
    "The `--overhang` argument can only be used when the flag `--split-arms` is set"
  )
}
#==========================#
#    PRE-REQUISITES END    #
#==========================#

#=======================#
#    FUNCTIONS START    #
#=======================#

#' Get the span of the non-padding portion of a pileup string.
#'
#' `GetSpan()` returns the first and last positions in a string that are not
#' `"."`. In the context of ASCII-style pileups, this identifies the coordinates
#' occupied by the represented element, such as a mature-arm marker, an aligned
#' read, or a reference sequence.
#'
#' @details The function is intended for pileup strings in which `"."` marks
#' padding or empty positions. It can therefore be used on:
#' \enumerate{
#'   \item mature-feature rows represented with `">"` characters,
#'   \item aligned read strings,
#'   \item genomic/reference sequence rows.
#' }
#'
#' In all cases, the function returns the interval covered by the non-dot
#' characters. If the string contains only `"."` characters, or no characters at
#' all, the function returns `NA` for both boundaries.
#'
#' @param str A character string from an ASCII-style pileup.
#'
#' @returns A named integer vector of length two with elements:
#'   \describe{
#'     \item{start}{The first position in `str` that is not `"."`.}
#'     \item{end}{The last position in `str` that is not `"."`.}
#'   }
#'   If no such position exists, both values are `NA_integer_`.
#'
#' @examples
#' GetSpan("....>>>>>>>>....")
#' GetSpan("....ACGT....")
#' GetSpan("ACTGAGGTCCTCAAAACTGAGG")
#' GetSpan("..........")
GetSpan <- function( str ) {
  # Split string into a character-based vector
  chars <- strsplit( str, "" )[[1]]
  # Retrieve non-dot positions
  pos <- which( chars != "." )

  if ( length( pos ) == 0) {
    return( c( start = NA_integer_, end = NA_integer_ ))
  }

  c( start = min( pos ), end = max( pos ))
}


#' Test whether an aligned substring lies within a reference span.
#'
#' `WithinSpan()` checks whether the aligned portion of a string falls fully
#' inside a reference span, optionally extended by a user-defined overhang.
#'
#' @param str A character string representing an aligned sequence padded with
#'   `"."`.
#' @param ref.span A named vector with elements `start` and `end` giving the
#'   reference span.
#' @param overhang A non-negative integer specifying how many positions outside
#'   the reference span are tolerated on each side.
#'
#' @returns `TRUE` if the aligned substring in `str` lies fully within the
#'   interval `[ref.span["start"] - overhang, ref.span["end"] + overhang]`;
#'   otherwise `FALSE`.
#'
#' @details
#' If `str` contains no aligned substring, the function returns `FALSE`.
#'
#' @examples
#' ref.span <- c(start = 5, end = 10)
#' WithinSpan("....ACGT....", ref.span, overhang = 0)
#' WithinSpan("..ACGT......", ref.span, overhang = 2)
WithinSpan <- function( str, ref.span, overhang ) {
  # Get reference string span
  span <- GetSpan( str = str)

  if ( any( is.na( span ))) {
    return( FALSE )
  }
  # Assess if string lies within reference span +/- overhang
  span[["start"]] >= ref.span[["start"]] - overhang &&
    span[["end"]] <= ref.span[["end"]]  + overhang
}


#' Trim a string to a reference span plus optional overhang.
#'
#' `CutToSpan()` extracts the substring defined by a reference span extended by
#' an optional overhang. The extracted interval is clipped to the actual string
#' boundaries.
#'
#' @param str A character string to trim.
#' @param span A named vector with elements `start` and `end`.
#' @param overhang A non-negative integer specifying how many extra positions to
#'   retain on each side of the span.
#'
#' @returns A character string corresponding to the requested interval.
#'
#' @examples
#' span <- c(start = 5, end = 10)
#' CutToSpan("ABCDEFGHIJKLMN", span, overhang = 2)
CutToSpan <- function( str, span, overhang ) {
  # Get substring start position
  start <- max( 1, span[[ "start" ]] - overhang )
  # Get substring end position
  end <- min( nchar( str ), span[[ "end" ]] + overhang )

  # Get substring
  substr(str, start, end)
}


#' Adjust genomic coordinates after trimming the reference sequence.
#'
#' `AdjustRefCoords()` recalculates the genomic coordinates of a reference row
#' after trimming the displayed reference sequence to an arm span plus optional
#' overhang.
#'
#' @param coord A character string with genomic coordinates in the format
#'   `"chr:start-end:strand"`.
#' @param span A named vector with elements `start` and `end` defining the arm
#'   span in the reference sequence.
#' @param overhang A non-negative integer specifying how many extra positions to
#'   keep on each side of the arm span.
#' @param ref.seq The full reference sequence string before trimming.
#'
#' @returns A character string containing the adjusted genomic coordinates in
#'   the same format as the input.
#'
#' @details
#' The function assumes that character positions in `ref.seq` correspond
#' directly to increasing genomic coordinates in `coord`.
#'
#' @examples
#' AdjustRefCoords(
#'   coord = "19:24242-24325:+",
#'   span = c(start = 14, end = 35),
#'   overhang = 4,
#'   ref.seq = "TCAGGCTGTGACCCTCCAGAGGGAAGTACTTTCTGTTG"
#' )
AdjustRefCoords <- function(coord, span, overhang, ref.seq) {
  # Retrieve coordinates' string parts
  coord.parts <- regexec("^(.+):(\\d+)-(\\d+):([+-])$", coord)
  coord.match <- regmatches(coord, coord.parts)[[1]]

  if (length(coord.match) == 0) {
    stop("Invalid coordinate format: ", coord)
  }

  # Assign string parts
  chrom <- coord.match[2]
  ref.start <- as.integer(coord.match[3])
  strand <- coord.match[5]

  # Get region window
  window.start <- max(1, span[["start"]] - overhang)
  window.end <- min(nchar(ref.seq), span[["end"]] + overhang)

  # Adjust coordinates
  new.start <- ref.start + window.start - 1
  new.end <- ref.start + window.end - 1

  # Construct new coordinate's string
  paste0(chrom, ":", new.start, "-", new.end, ":", strand)
}


#' Split a precursor pileup into arm-specific pileups.
#'
#' `SplitArms()` splits a precursor ASCII-style alignment pileup into one pileup
#' per mature arm. Each output pileup contains the arm annotation row, the
#' trimmed reference row, and only those reads whose aligned span falls within
#' the arm span plus the specified overhang.
#'
#' @param in.pileup A data frame representing one input pileup. It must contain
#'   the columns `seq` and `counts`.
#' @param head.lines An integer giving the number of header lines in the input
#'   pileup. Expected values are:
#'   \describe{
#'     \item{3}{single arm with precursor}
#'     \item{4}{two arms with precursor}
#'   }
#' @param overhang A non-negative integer specifying the allowed extension
#'   beyond the arm span on each side. If `NULL`, it is treated as `0`.
#'
#' @returns A list of data frames, one per mature arm.
#'
#' @details
#' For each mature arm, the function:
#' \enumerate{
#'   \item identifies the arm span from the arm annotation row,
#'   \item trims the arm row and reference row to that span plus overhang,
#'   \item adjusts the genomic coordinates of the reference row accordingly,
#'   \item retains only reads that lie fully within that interval,
#'   \item trims retained reads to the same interval.
#' }
#'
#' @examples
#' # SplitArms(pileup, head.lines = 4, overhang = 4)
SplitArms <- function(in.pileup, head.lines, overhang) {

  # Store read rows
  read.rows <- in.pileup %>%
    dplyr::slice(-seq_len(head.lines))

  # Store index of the arm's representation row(s)
  arm.rows <- if (head.lines == 3) 2 else 2:3

  #---> SPLIT PILEUP PER ARM <---#
  lapply(arm.rows, function(idx) {
    # Store arm span
    ref.span <- GetSpan(in.pileup$seq[idx])
    # Store reference genomic sequence and coordinates
    ref.seq <- in.pileup$seq[head.lines]
    ref.coord <- in.pileup$counts[head.lines]

    #---> MODIFY HEADER <---#
    header <- in.pileup %>%
      # Retrieve header lines
      dplyr::slice(c(idx, head.lines)) %>%
      # Cut header strings to defined window
      dplyr::mutate(
        seq = vapply(
          seq, CutToSpan, character(1), span = ref.span, overhang = overhang
        ),
        # Adjust coordinates to new window
        counts = c(
          counts[1],
          AdjustRefCoords(
            coord = ref.coord,
            span = ref.span,
            overhang = overhang,
            ref.seq = ref.seq
          )
        )
      )

    #---> MODIFY ALIGNED READS <---#
    arm.pileup <- read.rows %>%
      # Keep aligned reads that lie within defined window
      dplyr::filter(
        vapply(
          seq, WithinSpan, logical(1), ref.span = ref.span, overhang = overhang
        )
      ) %>%
      # Cut aligned reads to defined window
      dplyr::mutate(
        seq = vapply(
          seq, CutToSpan, character(1), span = ref.span, overhang = overhang
        )
      )
    # Bind header and aligned reads
    dplyr::bind_rows(header, arm.pileup)
  })
}


#' Filter a pileup by read count and maximum number of displayed sequences.
#'
#' `FilterPileup()` keeps only read rows with counts greater than or equal to
#' `min.count`, then keeps at most the first `max.seq` of those rows. The header
#' rows are preserved, and a separator row labeled `"Counts"` is inserted before
#' the retained reads.
#'
#' @param pileup A data frame representing one pileup, with columns `seq` and
#'   `counts`.
#' @param min.count Minimum count required for a read to be retained.
#' @param max.seq Maximum number of read rows to retain.
#'
#' @returns A filtered pileup data frame.
#'
#' @details
#' The function assumes that rows above the genomic reference row are header
#' rows and that rows below it correspond to aligned reads.
#'
#' @examples
#' # FilterPileup(pileup, min.count = 2, max.seq = 20)
FilterPileup <- function( pileup, min.count, max.seq ) {

  # Store last header line index
  head.lines <- which(grepl("^.*:", pileup[, 2]))

  # Get aligned reads rows
  read.rows <- pileup %>%
    dplyr::slice( -seq_len( head.lines ))

  read.rows <- read.rows %>%
    # Keep alignments with at least `min.count`
    dplyr::filter( as.numeric(counts) >= min.count ) %>%
    # Keep top `max.seq` rows
    dplyr::slice( 1: max.seq )

  pileup[1: head.lines, ] %>%
    # Create counts header row
    dplyr::add_row( seq = "", counts = "Counts") %>%
    # Bind final read alignments subset
    rbind( read.rows )
}


#' Mark the canonical aligned read in a pileup.
#'
#' `AddCanonical()` adds a column named `feat` to a pileup and marks the row or
#' rows corresponding to the canonical mature sequence for each detected arm.
#'
#' @param pileup A data frame representing one pileup, with columns `seq` and
#'   `counts`.
#'
#' @returns The input pileup with an additional column `feat`. Rows matching the
#'   canonical sequence are marked with `"< <feature_name>"` in that column,
#'   while all other rows contain an empty string.
#'
#' @details
#' The canonical sequence is inferred directly from the reference row using the
#' mature arm span:
#' \enumerate{
#'   \item identify the genomic reference row,
#'   \item identify mature arm rows located above the reference row,
#'   \item extract the arm span from each mature arm row,
#'   \item extract the corresponding substring from the reference row,
#'   \item match this sequence against read rows after removing `"."` padding,
#'   \item exclude any read containing `"-"` from being considered canonical.
#' }
#'
#' The function assumes that the row immediately below the reference row is the
#' `"Counts"` separator row and therefore does not consider it a read.
#'
#' @examples
#' # AddCanonical(pileup)
AddCanonical <- function( pileup ) {
  # Create empty field
  pileup$feat <- ""

  # Retrieve row with reference genomic coordinates
  ref.row <- which( grepl( "^.+:\\d+-\\d+:[+-]$", pileup$counts ))[1]
  if ( is.na( ref.row )) {
    warning( "Could not detect the reference row." )
    return( pileup )
  }

  # Retrieve row with arm representation
  arm.rows <- which( seq_len( nrow( pileup )) < ref.row
                     & grepl( "miR", pileup$counts) )
  if ( length( arm.rows ) == 0 ) {
    warning( "Could not detect any mature-arm rows." )
    return( pileup )
  }

  # Retrieve alignment rows index
  is.read <- seq_len( nrow( pileup )) > ( ref.row + 1 )
  # Set which alignments do not have deletions
  has.no.deletion <- !grepl( "-", pileup$seq, fixed = TRUE )
  # Remove flanking dots
  seq.no.dots <- gsub( "\\.", "", pileup$seq )

  for ( arm.row in arm.rows ) {
    # Get canonical sequence name
    feat.name <- pileup$counts[arm.row]
    # Get canonical sequence span
    ref.span <- GetSpan( pileup$seq[arm.row] )

    # Get canonical sequence
    canon.seq <- substr(
      pileup$seq[ref.row],
      ref.span[["start"]],
      ref.span[["end"]]
    )
    # Retrieve alignments within canonical sequence span
    in.arm <- rep( FALSE, nrow( pileup ))
    in.arm[is.read] <- vapply(
      pileup$seq[is.read],
      WithinSpan,
      logical(1),
      ref.span = ref.span,
      overhang = 0
    )

    # Retrieve canonical sequence index
    read.idx <- which( is.read &
                         in.arm &
                         has.no.deletion &
                         seq.no.dots == canon.seq
                       )

    # Indicate canonical sequence
    if ( length( read.idx ) > 0 ) {
      pileup$feat[read.idx] <- paste0( "< ", feat.name )
    } else {
      warning("Canonical sequence not found in pileup for feature: ", feat.name)
    }
  }

  pileup
}
#=====================#
#    FUNCTIONS END    #
#=====================#

#================#
#   MAIN START   #
#================#
#
# For each input pileup file:
#   1. Read the pileup into a two-column data frame.
#   2. Detect the number of header lines from the genomic coordinate row.
#   3. Optionally split precursor pileups into arm-specific pileups.
#   4. Filter each pileup by minimum count and maximum number of displayed reads.
#   5. Optionally annotate canonical read sequence(s).
#   6. Write the resulting pileup(s) to the output directory.
# ------------------------------------------------------------------------------

#---> START MESSAGE <---#
if ( verb ) cat( "Refromating ASCII-style pileup aesthetics...\n", sep = "" )

# If no overhang defined, no overhang allowed
if (is.null(overhang)) {
  overhang <- 0
}


#---> IMPORT FILES <---#
dir.files <- dir( in.dir, full.names = TRUE )
dir.pileups <- dir.files[grepl( "\\.tab$", dir.files )]

for ( file.pileup in dir.pileups ) {

  # Read ASCII-style pileup
  pileup <- read.csv( file.pileup,
                      stringsAsFactors = FALSE,
                      sep = "\t",
                      header = F,
                      col.names = c( "seq", "counts" ))

  # Store index of the last header line
  last.header.row <- which(grepl("^.*:", pileup[, 2]))

  # Skip empty pileups unless `--keep-all` is set
  if ( !keep.all && nrow( pileup ) <= last.header.row ) {
    next
  }

  #---> SPLIT ARMS  <---#
  if ( split.arms ) {
    # Print status message
    if ( verb ) cat( "Splitting arms...\n", sep = "" )

    list.pileup <- SplitArms(
      in.pileup = pileup,
      head.lines = last.header.row,
      overhang = overhang
    )
  } else {
    list.pileup <- list( pileup )
  }


  #---> FILTER PILEUPS  <---#
  # Print status message
  if ( verb ) cat( "Filtering pileups...\n", sep = "" )

  # Filter pileup
  list.pileup <- lapply(
    list.pileup,
    FilterPileup,
    min.count = min.count,
    max.seq = max.seq
)

  #---> ADD CANONICAL SEQUENCE <---#
  if ( canonical ) {
    # Print status message
    if ( verb ) cat( "Adding name to the canonical sequence...\n", sep = "" )

    list.pileup <- lapply( list.pileup, AddCanonical )
  }

  #---> WRITE OUTPUT  <---#
  # Print status message
  if ( verb ) cat( "Writing final pileups...\n", sep = "" )

  # Create output directory if it does not exist
  dir.create( out.dir, recursive = TRUE, showWarnings = FALSE )

  # Write tab-separated output
  lapply(
    list.pileup,
    function( x ){

      out.name <- file.path(
        out.dir,
        paste0(
          prefix,
          ".",
          x$counts[1],
          ".",
          overhang,
          "-shift.tab"
        )
      )

      write.table(
        x,
        file = out.name,
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
    }
  )
}

#---> END MESSAGE <---#
if ( verb ) cat( "Done.\n", sep = "" )
#================#
#    MAIN END    #
#================#
