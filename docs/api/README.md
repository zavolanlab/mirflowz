# API Overview

## Bash

- [`blocksort.sh`](./blocksort.sh.md): Sort oligomap alignments based on their
numerical names.
- [`get_lines_w_pattern.sh`](./get_lines_w_pattern.sh.md): Retrieve all lines
  with the matching pattern `--pattern` within the requested column
  (`--column`).
- [`trim_id_fasta.sh`](./trim_id_fasta.sh.md): Trim the sequence ID of a FASTA
  file to the first white-space.


## Perl

- [`map_chromosomes.pl`](./map_chromosomes.pl.md): Map/rename chromosome
  identifiers in a delimited text file using a tab-delimited mapping table.
- [`sam_remove_duplicates_inferior_alignments_multimappers.pl`](./sam_remove_duplicates_inferior_alignments_multimappers.pl.md):
  Removes duplicate records, then all `QNAME` duplicates except for the one(s)
  with the shortest edit distance. Optionally, multimappers (alignments of
  queries with the same edit distance, but different coordinates) are
  discarded.
- [`sam_trx_to_sam_gen.pl`](./sam_trx_to_sam_gen.pl.md): Re-maps a SAM file
  resulting from aligning a library of sequencing reads against a transcriptome
  to genomic coordinates.
- [`sam_uncollapse.pl`](./sam_uncollapse.pl.md): Reverses the collapsing of
  reads with identical sequences as done with `fastx_collapser`
  ([FASTX Toolkit][docs-fastx]) or similar.


## Python

### Modules

- [`annotate_sam_with_intersecting_features.py`](./annotate_sam_with_intersecting_features.py.md#module-annotate_sam_with_intersecting_featurespy):
  Annotate SAM alignments with their intersecting feature(s).
- [`filter_multimappers.py`](./filter_multimappers.py.md#module-filter_multimapperspy):
  Filter miRNA reads mapped to multiple locations by indel count.
- [`mirna_extension.py`](./mirna_extension.py.md#module-mirna_extensionpy):
  Extend miRNA start and end coordinates and ensure name uniqueness.
- [`mirna_quantification.py`](./mirna_quantification.py.md#module-mirna_quantificationpy):
  Quantify miRNAs and corresponding isomiRs.
- [`nh_filter.py`](./nh_filter.py.md#module-nh_filterpy): Filter alignments in
  a SAM file by NH tag.
- [`oligomap_output_to_sam_nh_filtered.py`](./oligomap_output_to_sam_nh_filtered.py.md#module-oligomap_output_to_sam_nh_filteredpy):
  Transform oligomap output FASTA file to SAM keeping the best alignments.
- [`primir_quantification.py`](./primir_quantification.py.md#module-primir_quantificationpy):
  Tabulate `bedtools intersect -wo -s` output file.
- [`validate_bedtools_intersect.py`](./validate_bedtools_intersect.py.md#module-validate_bedtools_intersectpy):
  Validation utilities for `bedtools intersect -wo -s` output.
- [`validation_fasta.py`](./validation_fasta.py.md#module-validation_fastapy):
  Filter FASTA files.

### Classes

- [`py.AnnotationException`](./mirna_extension.py.md#class-annotationexception):
  A custom exception class for `MirnaExtension` class.
- [`py.MirnaExtension`](./mirna_extension.py.md#class-mirnaextension): Class for
  updating miRNA annotated coordinates and names.
- [`py.Fields`](./oligomap_output_to_sam_nh_filtered.py.md#class-fields): Class
  to store an alignment in its different SAM fields.
- [`py.FileFormatError`](./validate_bedtools_intersect.py.md#class-fileformaterror):
  Raised when a line or field does not match the expected format.
- [`py.Record`](./validate_bedtools_intersect.py.md#class-record): A single
  validated `bedtools intersect -wo -s` record.

### Functions

- [`py.get_tags`](./annotate_sam_with_intersecting_features.py.md#function-get_tags):
  Construct a custom tag for an alignment based on intersecting features.
- [`py.main`](./annotate_sam_with_intersecting_features.py.md#function-main):
  Annotate alignments in a SAM file with intersecting feature tags.
- [`py.parse_arguments`](./annotate_sam_with_intersecting_features.py.md#function-parse_arguments):
  Command-line arguments parser.
- [`py.parse_intersect_output`](./annotate_sam_with_intersecting_features.py.md#function-parse_intersect_output):
  Parse `bedtools intersect -wo -s` output file.
- [`py.count_indels`](./filter_multimappers.py.md#function-count_indels): Count
  the number of indels in an alignment based on its CIGAR string.
- [`py.find_best_alignments`](./filter_multimappers.py.md#function-find_best_alignments):
  Find alignments with more indels.
- [`py.main`](./filter_multimappers.py.md#function-main): Filter multimappers
  by indels count.
- [`py.parse_arguments`](./filter_multimappers.py.md#function-parse_arguments):
  Command-line arguments parser.
- [`py.write_output`](./filter_multimappers.py.md#function-write_output):
  Write the output to the standard output (`STDOUT`).
- [`py.MirnaExtension.__init__`](./mirna_extension.py.md#function-__init__):
  Initialize class.
- [`py.MirnaExtension.adjust_names`](./mirna_extension.py.md#function-adjust_names):
  Adjust miRNA attributes for uniqueness and consistency.
- [`py.MirnaExtension.process_precursor`](./mirna_extension.py.md#function-process_precursor):
  Extend miRNAs start and end coordinates and ensure name uniqueness.
- [`py.MirnaExtension.set_db`](./mirna_extension.py.md#function-set_db):
  Load GFF3 file into `gffutils.FeatureDB`.
- [`py.MirnaExtension.set_seq_lengths`](./mirna_extension.py.md#function-set_seq_lengths):
  Set the reference sequence lengths.
- [`py.MirnaExtension.update_db`](./mirna_extension.py.md#function-update_db):
  Update miRNA annotations in the local database.
- [`py.MirnaExtension.write_gff`](./mirna_extension.py.md#function-write_gff):
  Write features to a GFF3 file.
- [`py.main`](./mirna_extension.py.md#function-main):
  Extend miRNAs start/end coordinates.
- [`py.parse_arguments`](./mirna_extension.py.md#function-parse_arguments):
  Parse command-line arguments.
- [`py.collapsed_contribution`](./mirna_quantification.py.md#function-collapsed_contribution):
  Get the contribution of the alignment to the overall count.
- [`py.collapsed_nh_contribution`](./mirna_quantification.py.md#function-collapsed_nh_contribution):
  Get the contribution of the alignment to the overall count.
- [`py.contribution`](./mirna_quantification.py.md#function-contribution): Get
  the contribution of the alignment to the overall count.
- [`py.get_name`](./mirna_quantification.py.md#function-get_name): Get the final
  name for the species name.
- [`py.main`](./mirna_quantification.py.md#function-main): Quantify miRNAs and
  corresponding isomiRs.
- [`py.nh_contribution`](./mirna_quantification.py.md#function-nh_contribution):
  Get the contribution of the alignment to the overall count.
- [`py.parse_arguments`](./mirna_quantification.py.md#function-parse_arguments):
  Command-line arguments parser.
- [`py.write_output`](./mirna_quantification.py.md#function-write_output):
  Write to the output the correct miRNA type.
- [`py.main`](./nh_filter.py.md#function-main): Filter alignments by their NH
  tag value.
- [`py.parse_arguments`](./nh_filter.py.md#function-parse_arguments): Parse
  command-line arguments.
- [`py.eval_aln`](./oligomap_output_to_sam_nh_filtered.py.md#function-eval_aln):
  Evaluate an alignment to add, discard or write it to the STDOUT.
- [`py.get_cigar_md`](./oligomap_output_to_sam_nh_filtered.py.md#function-get_cigar_md):
  Get the CIGAR and MD strings.
- [`py.get_sam_fields`](./oligomap_output_to_sam_nh_filtered.py.md#function-get_sam_fields):
  Create the read's alignment in SAM format.
- [`py.main`](./oligomap_output_to_sam_nh_filtered.py.md#function-main):
  Convert the alignments in the oligomap output file to SAM format.
- [`py.parse_arguments`](./oligomap_output_to_sam_nh_filtered.py.md#function-parse_arguments):
  Command-line arguments parser.
- [`py.get_contribution`](./primir_quantification.py.md#function-get_contribution):
  Return the contribution of a single alignment.
- [`py.get_initial_data`](./primir_quantification.py.md#function-get_initial_data):
  Get the feature name and its extension.
- [`py.main`](./primir_quantification.py.md#function-main): Tabulate
  `bedtools intersect -wo -s` output file.
- [`py.parse_arguments`](./primir_quantification.py.md#function-parse_arguments):
  Command-line arguments parser.
- [`py.Record.__init__`](./validate_bedtools_intersect.py.md#function-__init__):
  Initialize a validated record.
- [`py.Record.cross_checks`](./validate_bedtools_intersect.py.md#function-cross_checks):
  Validate relationships across fields.
- [`py.parse_all`](./validate_bedtools_intersect.py.md#function-parse_all):
  Stream a file, yielding `(line_number, Record)` for each data line.
- [`py.validate_first_n`](./validate_bedtools_intersect.py.md#function-validate_first_n):
 Validate the first `n` lines of a `bedtools intersect -wo -s` out file.
- [`py.compile_trim_pattern`](./validation_fasta.py.md#function-compile_trim_pattern):
  Get a compiled regex pattern to trim at a character's first occurrence.
- [`py.main`](./validation_fasta.py.md#function-main): Filter and process a
  FASTA file.
- [`py.open_fasta`](./validation_fasta.py.md#function-open_fasta): Open a FASTA
  or FASTA.GZ for text‐mode reading.
- [`py.parse_and_validate_arguments`](./validation_fasta.py.md#function-parse_and_validate_arguments):
  Parse and validate command-line arguments.
- [`py.trim_id`](./validation_fasta.py.md#function-trim_id): Trim a FASTA ID
  using the first-occurrence of any character in `_pattern`.
- [`py.write_id_file`](./validation_fasta.py.md#function-write_id_file): Write
  the final sequence IDs, one per line.


## R

- [`ascii_alignment_pileup.R`](./ascii_alignment_pileup.R.md): Generates an
  ASCII-style pileup of read alignments in one or more BAM files against one
  or more regions specified in a BED file.
- [`gtf_exons_bed.1.1.2.R`](./gtf_exons_bed.R.md): Converts the exon entries
  of a GTF file to a BED file with one line per exon.
- [`merge_tables.R`](./merge_tables.R.md): Merge miRNAs quantification tables.

---
