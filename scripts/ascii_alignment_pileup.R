#!/usr/bin/env Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Sep 29, 2019
### Author: Alexander Kanitz
### Affiliation: Zavolan Group, Biozentrum, University of Basel
### Email: alexander.kanitz@alumni.ethz.ch
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTION PARSER <---#
if (suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
args.all <- commandArgs(trailingOnly = FALSE)
name.field <- "--file="
script <- basename(sub(name.field, "", args.all[grep(name.field, args.all)]))

#---> DESCRIPTION <---#
description <- "Generates an ASCII-style pileup of read alignments in one or more BAM files
against one or more regions specified in a BED file.\n"
author <- "Author: Alexander Kanitz"
affiliation <- "Affiliation: Biozentrum, University of Basel"
email <- "Email: alexander.kanitz@alumni.ethz.ch"
version <- "1.2.1"
version_formatted <- paste("Version:", version, sep=" ")
requirements <- c("optparse", "rtracklayer", "GenomicAlignments", "tools")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, affiliation, email, version_formatted, requirements_txt, sep="\n")
notes <- "Notes:
- For the input queries, consider the `--maximum-region-width` parameter, which
  is provided for safety. While it is possible to increase it, wide regions may
  require excessive memory.
"

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
        make_option(
            "--reference",
            action="store",
            type="character",
            default=NULL,
            help="Reference genome sequence in FASTA format. The file *MUST* be compressed with
            BGZIP. If supplied, the reference sequence for the query region(s) will be added to
            the output. Note that on the first run with a specific reference genome file, an FAI
            index is generated which will take some time.",
            metavar="file"
        ),
        make_option(
            "--annotations",
            action="store",
            type="character",
            default=NULL,
            help="Annotation file in GFF/GTF format used to annotate sequences. If supplied,
            features overlapping the query region(s) will be visualized in the output. Ensure that
            the argument to option `annotation-name-field` corresponds to a field in the
            annotations, otherwise the script will fail.",
            metavar="file"
        ),
        make_option(
            "--output-directory",
            action="store",
            type="character",
            default=getwd(),
            help="Output directory. One output file will be created for each region in `--bed` and
            the filenames will be generated from the basenames of the supplied BAM file(s) and the
            name field (4th column) of the BED file. [default \"%default\"]",
            metavar="dir"
        ),
        make_option(
            "--maximum-region-width",
            action="store",
            type="integer",
            default=200,
            help="Maximum input region width. Use with care as wide regions will use excessive
            resources. [default %default]",
            metavar="int"
        ),
        make_option(
            "--do-not-collapse-alignments",
            action="store_true",
            type="logical",
            default=FALSE,
            help="Show alignments of reads with identical sequences individually."
        ),
        make_option(
            "--minimum-count",
            action="store",
            type="integer",
            default=1,
            help="Alignments of reads with less copies than the specified number will not be
            printed. Option is not considered if `do-not-collapse-alignments` is set.
            [default %default]",
            metavar="int"
        ),
        make_option(
            "--annotation-name-field",
            action="store",
            type="character",
            default="Name",
            help="Annotation field used to populate the `name` column in the output.
            [default \"%default\"]",
            metavar="str"
        ),
        make_option(
            "--padding-character",
            action="store",
            default=".",
            help="Character used for padding alignments. [default \"%default\"]",
            metavar="char"
        ),
        make_option(
            "--indel-character",
            action="store",
            default="-",
            help="Character to denote insertions and deletions in alignments. [default \"%default\"]",
            metavar="char"
        ),
        make_option(
            "--prefix",
            action="store",
            type="character",
            default=NULL,
            help="Prefix to be used in the output file name(s). If not provided
            the input BAM file(s) name will be used instead,",
            metavar="string"
        ),
        make_option(
            c("-h", "--help"),
            action="store_true",
            default=FALSE,
            help="Show this information and die."
        ),
        make_option(
            c("-v", "--verbose"),
            action="store_true",
            default=FALSE,
            help="Print log messages to STDOUT."
        )
)

## Parse options
opt_parser <- OptionParser(
    usage=paste(script, "[--help] [--verbose] [OPTIONS] BED BAM [BAM2 ...]\n"),
    option_list=option_list,
    add_help_option=FALSE,
    description=msg,
    epilogue=notes
)
cli <- parse_args(opt_parser, positional_arguments=c(2, Inf))

# Re-assign CLI arguments
fl.query <- cli$args[[1]]
fl.bam <- cli$args[2:length(cli$args)]
fl.ref <- cli$options[["reference"]]
fl.anno <- cli$options[["annotations"]]
dir.out <- cli$options[["output-directory"]]
prefix.out <- cli$options[["prefix"]]
width.max <- cli$options[["maximum-region-width"]]
collapse <- ! cli$options[["do-not-collapse-alignments"]]
count.min <- cli$options[["minimum-count"]]
char.pad <- cli$options[["padding-character"]]
char.indel <- cli$options[["indel-character"]]
field.name.anno <- cli$options[["annotation-name-field"]]
verb <- cli$options[["verbose"]]
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if (verb) cat("Starting '", script, "'...\n", sep="")

#---> LOAD REQUIRED LIBRARIES <---#
if (verb) cat("Loading libraries...\n", sep="")
for (req in requirements) {
    if ( suppressWarnings(suppressPackageStartupMessages(require(req, character.only=TRUE))) == FALSE ) { stop("Package '", req, "' required!") }
}

#---> IMPORT FILES <---#
# Print status message
if (verb) cat("Importing input files...\n", sep="")
# Load files
bed <- import(con=fl.query)
if (! is.null(fl.ref)) {ref <- FaFile(fl.ref)}
if (! is.null(fl.anno)) {anno <- import(con=fl.anno)}

# Get file prefix from BAM files or from CLI argument
fl.prefix <- if (! is.null(prefix.out)) prefix.out else paste(basename(file_path_sans_ext(fl.bam)), collapse=".")

#--->   <---#
# Print status message
if (verb) cat("Iterating over regions in BED file...\n")
# Iterate over input regions
for(index in seq_along(bed)) {

    #---> PREPARE STACKING OF ALIGNMENTS  <---#
    # Assign current region
    region <- bed[index]
    # Print status message
    if (verb) cat("Processing region '", mcols(region)[["name"]], "'...\n", sep="")
    # Exit with error if region is too wide
    if (width(region) > width.max) {stop("Supplied region too large. Consider increasing the `width.max` parameter, but note that for very large regions, the memory footprint may be excessive.")}
    # Initialize DNAStringSet container object
    seq.out <- DNAStringSet()

    #---> ADD READ ALIGNMENTS  <---#
    # Print status message
    if (verb) cat("Iterating over BAM files...\n")
    # Iterate over BAM files
    for (bam in fl.bam) {
        # Print status message
        if (verb) cat("Adding read alignments for file '", bam,"'...\n", sep="")
        # Get alignments for current BAM file
        seq.stack <- stackStringsFromBam(
            bam,
            param=region,
            D.letter=char.indel,
            N.letter=char.indel,
            Lpadding.letter=char.pad,
            Rpadding.letter=char.pad
        )
        # Add to container
        seq.out <- append(seq.out, seq.stack)
    }
    # Get complement if on minus strand
    if (as.character(strand(region))[[1]] == "-") {
        seq.out <- complement(seq.out)
    }
    # Convert to character
    seq.out <- as.character(seq.out)
    # Convert to dataframe for writing output
    df <- data.frame(seq=seq.out, count=rep(NA, length(seq.out)), stringsAsFactors=FALSE)

    #---> COLLAPSE READ ALIGNMENTS <---#
    if (collapse) {
        # Print status message
        if (verb) cat("Collapsing identical reads/alignments...\n")
        # Get unique alignments
        df <- data.frame(table(df[["seq"]]), stringsAsFactors=FALSE)
        if (! ncol(df) == 2) df <- data.frame(matrix(ncol = 2, nrow = 0))
        colnames(df) <- c("seq", "count")
        df[["seq"]] <- as.character(df[["seq"]])
        # Filter out any alignments that do not make the specified minimum count cutoff
        df <- df[df[["count"]] >= count.min, ]
    }

    #---> SORTING ALIGNMENTS  <---#
    # Print status message
    if (verb) cat("Sorting alignments...\n")
    # Reverse sequence if on minus strand
    if (as.character(strand(region))[[1]] == "-") {
        df[["seq"]] <- reverse(df[["seq"]])
    }
    # Sort by position of first nucleotide, count and position of last nucleotide
    if (nrow(df)) {
        last_char <- nchar(df[["seq"]][[1]])
        pos.nuc.first <- regexpr(paste0("[^", char.pad, "\\.]"), df[["seq"]])
        pos.nuc.last <- last_char - regexpr(paste0("[^", char.pad, "\\.]"), unlist(lapply(df[["seq"]], reverse))) + 1
        df <- df[order(
            pos.nuc.first, df[["count"]], pos.nuc.last,
            decreasing=c(FALSE, TRUE, FALSE)
        ), ]
    }
    # Reverse sequence again if on minus strand
    if (as.character(strand(region))[[1]] == "-") {
        df[["seq"]] <- reverse(df[["seq"]])
    }

    #---> ADD REFERENCE SEQUENCE  <---#
    if (! is.null(fl.ref)) {
        # Print status message
        if (verb) cat("Adding reference sequence...\n")
        # Generate name for reference sequence
        name.ref <- paste(
            seqnames(region),
            paste(start(region), end(region), sep="-"),
            strand(region), sep=":"
        )
        # Get sequence
        row.ref <- getSeq(ref, region)  # will create .fai index if not present
        # Get reverse if on minus strand
        if (as.character(strand(region))[[1]] == "-") {
            row.ref <- reverse(row.ref)
        }
        # Compile row and add to dataframe
        df.ref <- data.frame(seq=as.character(row.ref), count=name.ref, stringsAsFactors=FALSE)
        df <- rbind(df.ref, df)
    }

    #---> ADD ANNOTATIONS FOR REGION <---#
    # Add annotations
    if (! is.null(fl.anno)) {
        # Print status message
        if (verb) cat("Adding overlapping features...\n")
        # Find overlapping features
        features <- anno[to(findOverlaps(region, anno))]
        # Set order of addition from most upstream to most downstream
        if (as.character(strand(region)) == "+") {
            order.features <- order(start(features))
        } else {
            order.features <- rev(order(end(features)))
        }
        # Order features
        features <- features[order.features]
        # Iterate over features
        seqs <- sapply(seq_along(features), function(index) {
            feat <- features[index]
            # Set markup character depending on strand
            char.anno <- ifelse(strand(region) == "+", ">", "<")
            # Calculate lengths of feature (in region) & left/right padding
            diff.start <- start(feat) - start(region)
            n.padl <- max(0, diff.start)
            n.anno <- min(min(0, diff.start) + width(feat), width(region) - n.padl)
            n.padr <- width(region) - n.padl - n.anno
            # Generate string encompassing feature
            str.padl <- paste(rep(char.pad, n.padl), collapse="")
            str.anno <- paste(rep(char.anno, n.anno), collapse="")
            str.padr <- paste(rep(char.pad, n.padr), collapse="")
            str.final <- paste0(str.padl, str.anno, str.padr)
        })
        df.anno <- data.frame(seq=seqs, count=mcols(features)[[field.name.anno]], stringsAsFactors=FALSE)
        df <- rbind(df.anno, df)
    }

    #---> WRITE OUTPUT  <---#
    if (collapse) {
        name.file <- paste(fl.prefix, mcols(region)[["name"]], "min", as.character(count.min), "pileup", "tab", sep=".")
    } else {
        name.file <- paste(fl.prefix, mcols(region)[["name"]], "uncollapsed", "pileup", "tab", sep=".")
    }
    fl.out <- file.path(dir.out, name.file)
    # Print status message
    if (verb) cat("Writing output to file '", fl.out, "'...\n", sep="")
    # Write tab-separated output
    write.table(df, fl.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

#---> END MESSAGE <---#
if (verb) cat("Done.\n\nSession info:\n")
if (verb) print(sessionInfo())
#================#
#    MAIN END    #
#================#
