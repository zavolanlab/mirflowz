# MIRFLOWZ: workflow documentation

This document describes the individual steps of the workflow. For instructions
on installation and usage please see [here](README.md).

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Third-party software used](#third-party-software-used)
- [Description of workflow steps](#description-of-workflow-steps)
  - [Rule graph](#rule-graph)
  - [Preparatory](#preparatory)
    - [Read sample table](#read-sample-table)
    - [Configuration file](#configuration-file)
  - [Snakefile](#snakefile)
    - [`finish`](#finish)
  - [Prepare workflow](#prepare-workflow)
    - [`finish_prepare`](#finish_prepare)
    - [`trim_genome_seq_ids`](#trim_genome_seq_ids)
    - [`extract_transcriptome_seqs`](#extract_transcriptome_seqs)
    - [`trim_transcriptome_seq_ids`](#trim_transcriptome_seq_ids)
    - [`generate_segemehl_index_transcriptome`](#generate_segemehl_index_transcriptome)
    - [`generate_segemehl_index_genome`](#generate_segemehl_index_genome)
    - [`get_exons_gtf`](#get_exons_gtf)
    - [`convert_exons_gtf_to_bed`](#convert_exons_gtf_to_bed)
    - [`create_genome_header`](#create_genome_header)
    - [`map_chr_names`](#create_chr_names)
    - [`create_index_genome_fasta`](#create_index_genome_fasta)
    - [`extract_chr_len`](#extract_chr_len)
    - [`extend_mirs_annotations`](#extend_mirs_annotations)
  - [Map workflow](#map-workflow)
    - [`finish_map`](#finish_map)
    - [`start`](#start)
    - [`fastq_quality_filter`](#fastq_quality_filter)
    - [`fastq_to_fasta`](#fastq_to_fasta)
    - [`format_fasta`](#format_fasta)
    - [`remove_adapters`](#remove_adapters)
    - [`collapse_indentical_reads`](#collapse_indentical_reads)
    - [`map_genome_segemehl`](#map_genome_segemehl)
    - [`map_transcriptome_segemehl`](#map_transcriptome_segemehl)
    - [`filter_fasta_for_oligomap`](#filter_fasta_for_oligomap)
    - [`map_genome_oligomap`](#map_genome_oligomap)
    - [`sort_genome_oligomap`](#sort_genome_oligomap)
    - [`convert_genome_to_sam_oligomap`](#convert_genome_to_sam_oligomap)
    - [`map_transcriptome_oligomap`](#map_transcriptome_oligomap)
    - [`sort_transcriptome_oligomap`](#sort_transcriptome_oligomap)
    - [`convert_transcriptome_to_sam_oligomap`](#convert_transcriptome_to_sam_oligomap)
    - [`merge_genome_maps`](#merge_genome_maps)
    - [`merge_transcriptome_maps`](#merge_transcriptome_maps)
    - [`filter_genome_by_nh`](#filter_genome_by_nh)
    - [`filter_transcriptome_by_nh`](#filter_transcriptome_by_nh)
    - [`remove_header_genome_mappings`](#remove_header_genome_mappings)
    - [`remove_header_transcriptome_mappings`](#remove_header_transcriptome_mappings)
    - [`transcriptome_to_genome_maps`](#transcriptome_to_genome_maps)
    - [`merge_all_maps`](#merge_all_maps)
    - [`add_header_all_maps`](#add_header_all_maps)
    - [`sort_maps_by_id`](#sort_maps_by_id)
    - [`remove_inferiors`](#remove_inferiors)
    - [`filter_by_indels`](#filter_by_indels)
    - [`convert_all_alns_sam_to_bam`](#convert_all_alns_sam_to_bam)
    - [`sort_all_alns_bam_by_position`](#sort_all_alns_bam_by_position)
    - [`index_all_alns_bam`](#index_all_alns_bam)
  - [Quantify workflow](#quantify-workflow)
    - [`finish_quantify`](#finish_quantify)
    - [`intersect_extended_primir`](#intersect_extended_primir)
    - [`filter_sam_by_intersecting_primir`](#filter_sam_by_intersecting_primir)
    - [`convert_intersecting_primir_sam_to_bam`](#convert_intersecting_primir_sam_to_bam)
    - [`sort_intersecting_primir_bam_by_position`](#sort_intersecting_primir_bam_by_position)
    - [`index_intersecting_primir_bam`](#index_intersecting_primir_bam)
    - [`intersect_extended_mirna`](#intersect_extended_mirna)
    - [`filter_sam_by_intersecting_mirna`](#filter_sam_by_intersecting_mirna)
    - [`add_intersecting_mirna_tag`](#add_intersecting_mirna_tag)
    - [`sort_intersecting_mirna_by_feat_tag`](#sort_intersecting_mirna_by_feat_tag)
    - [`quantify_mirna`](#quantify_mirna)
    - [`quantify_primir`](#quantify_primir)
    - [`merge_tables`](#merge_tables)
    - [`uncollapse_reads`](#uncollapse_reads)
    - [`convert_uncollapse_reads_sam_to_bam`](#convert_uncollapse_reads_sam_to_bam)
    - [`sort_uncollapse_reads_bam_by_position`](#sort_uncollapse_reads_bam_by_position)
    - [`index_uncollapse_reads_bam`](#index_uncollapse_reads_bam)
  - [Pileup workflow](#pileup-workflow)
    - [`finish_pileup`](#finish_pileup)
    - [`create_empty_bed`](#create_empty_bed)
    - [`compress_reference_genome`](#compress_reference_genome)
    - [`create_per_library_ascii_pileups`](#create_per_library_ascii_pileups)
    - [`create_per_run_ascii_pileups`](#create_per_run_ascii_pileups)
    - [`create_per_condition_ascii_pileups`](#create_per_condition_ascii_pileups)



## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
| **ASCII-style alignment pileups** | [Apache 2.0][license-apache2] | _"Generates an ASCII-style pileup of read alignments in one or more BAM files against one or more regions specified in a BED file."_ | [code][code-ascii] | 
| **BEDTools** | [GPLv2][license-gpl2] | _"[...] intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][manual-bedtools] / [publication][pub-bedtools] |
| **cufflinks** | [BSL-1.0][license-bsl1] | _"[...] assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples"_ | [code][code-cufflinks] / [manual][docs-cufflinks] / [publication][pub-cufflinks] |
| **cutadapt** | [MIT][license-mit] | _"[...] finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads"_ | [code][code-cutadapt] / [manual][docs-cutadapt] / [publication][pub-cutadapt] |
| **FASTX-Toolkit** | [AGPL-3.0][license-agpl3] | _"[...] collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing"_ | [code][code-fastx] / [manual][docs-fastx] |
| **GFFUtils** | [AFL-3][license-afl3] | _"[...]  a small set of utility programs for working with GFF and GTF files"_ | [code][code-gffutils] / [manual][docs-gffutils] |
| **Oligomap** | [GPLv3][license-gpl3] | _"[...] fast identification of nearly-perfect matches of small RNAs in sequence databases. It allows to exhaustively identify all the perfect and 1-error (where an error is defined to be a mismatch, insertion or deletion) matches of large sets of small RNAs to target sequences"_ | [code][code-oligomap] / [publication][pub-oligomap] |
| **SAMtools** | [MIT][license-mit] | _"[...] suite of programs for interacting with high-throughput sequencing data"_ | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools] |
| **segemehl** | [GPLv3][license-gpl3] | _"[...]  map short sequencer reads to reference genomes"_ | [manual][docs-segemehl] / [publication][pub-segemehl] |

## Description of workflow steps

> The workflow consists of four Snakemake files: A main `Snakefile` and an
individual Snakemake file for each step in the workflow (the file preparation,
the reads mapping, the miRNA quantification and the ASCII-style alignment
pileups generation). The main `Snakefile` contains the configuration file
validation along with the inclusion of the sub-workflows. Individual steps of
the workflow are described briefly, and links to the respective software
manuals are given. Parameters that can be modified by the user (via the samples
table and the configuration file) are also described.

### Rule graph

![rule_graph][rule-graph]

Visual representation of the workflow. Automatically prepared with
[Snakemake][docs-snakemake].

### Preparatory

#### Read sample table

##### Requirements

- tab-separated values (`.tsv`) file
- First row has to contain parameter names as in 
[`samples_table.tsv`](test/test_files/samples_table.tsv) 
- First column used as sample identifiers

Parameter name | Description | Data type(s)
 --- | --- | --- 
sample | Descriptive sample name. | `str`
sample_file | Path of the library file in either `.fa.gz` or `.fastq.gz` format. | `str`
adapter | Required for [Cutadapt](#third-party-software-used). Use a value such as `XXXXXXXXXX` if no adapter is present or if no trimming is desired | `str`
format | There are two allowed values, `fa` and `fastq` according to the library format | `str`


#### Configuration file

Some parameters within the workflow can be modified. Refer to the
[configuration template](#config/config_template.yaml) for a detailed
explanation of each option.

### Snakefile

#### `finish`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - pri-miR intersections file (`.bed`); from
  [**intersect_extended_primir**](#intersect_extended_primir)
  - miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect_extended_mirna)
  - Alignments file, sorted and tagged (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort_intersecting_mirna_by_feat_tag)
  - (iso)miR and/or pri-miR counts table (`.tab`); from
  [**merge_tables**](#merge_tables)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Empty text file (`.txt`); from
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


### Prepare workflow

#### `finish_prepare`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - segemehl genome index (`idx`); from
  [**generate_segemehl_index_genome**](#generate_segemehl_index_genome)
  - segemehl transcriptome index (`idx`); from
  [**generate_segemehl_index_transcriptome**](#generate_segemehl_index_transcriptome)
  - Exon annotations (`.bed`); from
  [**convert_exons_gtf_to_bed**](#convert_exons_gtf_to_bed)
  - Genome header (`.sam`); from
  [**create_genome_header**](#create_genome_header)
  - Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`);
  from [**extract_chr_len**](#extract_chr_len)
  - Primary miRNA transcript (pri-miR) extended annotation (`.gff3`);
  from [**extend_mirs_annotations**](#extend_mirs_annotations)
  - Mature miRNA (miRNA) extended annotation (`.gff3`);
  from [**extend_mirs_annotations**](#extend_mirs_annotations)


#### `trim_genome_seq_ids`

Trim genome sequence IDs with a [**custom script**][custom-script-trim-id].

- **Input**
  - Genome sequence file (`.fasta`)
- **Output**
  - Genome sequence file, trimmed IDs (`.fasta`); used in
  [**extract_transcriptome_seqs**](#extract_transcriptome_seqs),
  [**create_genome_header**](#create_genome_header),
  [**create_index_genome_fasta**](#create_index_genome),
  [**generate_segemehl_index_genome**](#generate_segemehl_index_genome),
  [**mapping_genome_segemehl**](#mapping_genome_segemehl),
  [**mapping_genome_oligomap**](#mapping_genome_oligomap) and
  [**compress_reference_genome**](#compress_reference_genome)


#### `extract_transcriptome_seqs`

Create transcriptome from genome and genome annotations with 
[**cufflinks**](#third-party-software-used).

- **Input**
  - Genome sequence file (`.fasta`)
  - Genome annotation file (.`gtf`)
- **Output**
  - Transcriptome sequence file (`.fasta`); used in 
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)


#### `trim_transcriptome_seq_ids`

Trim transcriptome sequence IDs with a [**custom script**][custom-script-trim-id].

- **Input**
  - Transcriptome sequence file (`.fasta`)
- **Output**
  - Transcriptome sequence, trimmed IDs (`.fasta`); used in
  [**generate_segemehl_index_transcriptome**](#generate_segemehl_index_transcriptome),
  [**mapping_transcriptome_segemehl**](#mapping_transcriptome_segemehl) and
  [**mapping_transcriptome_oligomap**](#mapping_transcriptome_oligomap)


#### `generate_segemehl_index_transcriptome`

Generate transcriptome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The transcriptome index only needs to be generated once for each combination
of genome and annotations and sample sets.

- **Input**
  - Transcriptome sequence file, trimmed IDs (`.fasta`); from
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)
- **Output**
  - segemehl transcriptome index (`.idx`); used in 
  [**mapping_genome_segemehl**](#mapping_genome_segemehl)


#### `generate_segemehl_index_genome`

Generate genome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The genome index only needs to be generated once for each combination
of annotations and sample sets.

- **Input**
  - Genome sequence file with trim IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - segemehl genome index (`.idx`); used in 
  [**mapping_genome_segemehl**](#mapping_genome_segemehl)


#### `get_exons_gtf`

Retrieve exon annotations from genome annotations with a
[**custom script**][custom-script-get-lines].

- **Input**
  - Genomic annotations (`.gtf`)
- **Output**
  - Exon annotations (`.gtf`); used in 
  [**convert_exons_gtf_to_bed**](#convert_exons_gtf_to_bed)

#### `convert_exons_gtf_to_bed`

Convert exon annotations `.gtf` to `.bed` with a
[**custom script**][custom-script-gtf-bed].

- **Input**
  - Exon annotations (`.gtf`); from [**get_exons_gtf**](#get_exons_gtf)
- **Output**
  - Exon annotations (`.bed`); used in
  [**transcriptome_to_genome_maps**](#transcriptome_to_genome_maps)


#### `create_genome_header`

Create `SAM` header for the genome with 
[**SAMtools**](#third-party-software-used).

> Required by [SAMtools](#third-party-software-used) to work with the alignment
file.

- **Input**
  - Genome sequence file, trimmed IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - Genome header (`.sam`); used in [add_header_all_maps](#add_header_all_maps)


#### `map_chr_names`

Map UCSC-like chromosome names with Ensembl-like ones in miRNA annotation
with a [**custom script**][custom-script-map-chr].

> Required by [BEDTools](#third-party-software) to intersect alignments with
miRNA annotations. Several mapping tables are available [here][chr-maps].

- **Input**
  - miRNA annotations (`.gff3`)
  - Tab-separated mappings table (`.tsv`)
- **Output**
  - miRNA annotations with mapped genes(`.gff3`); used in 
  [**extend_mirs_annotations**](#extend_mirs_annotations),
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `create_index_genome_fasta`

Create a `FASTA` index for the genome with 
[**SAMtools**](#third-party-software-used).

- **Input**
  - Genome sequence file, trimmed IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - `FASTA` genome index (`.fa.fai`); used in
  [**extract_chr_len**](#extract_chr_len)


#### `extract_chr_len`

Extract chromosome(s) length from the genome sequence.

- **Input**
  - `FASTA` genome index (`.fa.fai`); from
  [**create_index_genome_fasta**](#create_index_genome_fasta)
- **Output**
  - Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`); used
  in [**extend_mirs_annotations**](#extend_mirs_annotations)


#### `extend_mirs_annotations`

Extend miRNA annotations and split the file by feature with a
[**custom script**][custom-script-mir-ext].

> Mature miRNA regions are extended on both sides to account for isomiR species
with shifted start and/or end positions. If required, pri-miR loci are also 
extended to accommodate the new miRNA coordinates. 

- **Input**
  - miRNA annotations with mapped chromosomes(`.gff3`); from
  [**map_chr_names**](#map_chr_names)
- **Parameters**
  - **config_template.yaml**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended (default 6)
- **Output**
  - Primary miRNA transcript (pri-miR) extended annotation (`.gff3`); used in
  [**intersect_extended_primir**](#intersect_extended_primir)
  - Mature miRNA (miRNA) extended annotation (`.gff3`); used in
  [**intersect_extended_mirna**](#intersect_extended_mirna)

### Map workflow

#### `finish_map`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - `BAM` index file (`.bam.bai`); from
  [**index_all_alns_bam**](#index_all_alns_bam)



#### `start`

Copy and rename read files.

> Local rule.
Depending on the read files format, the output files undergo a quality filter
(`.fastq`) or are directly formatted (`.fa`).

- **Input**
  - Reads file (`.fa.gz`, `.fastq.gz`)
- **Output**
  - Reads file, copied, renamed (`.fa`, `.fastq`); used in
  [**fastq_quality_filter**](#fastq_quality_filter) or 
  [**format_fasta**](#format_fasta)


#### `fastq_quality_filter`

Conduct quality control for reads library with 
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file, copied, renamed (`.fastq`); from [**start**](#start)
- **Parameters**
  - **config_template.yaml**
    - `q_value`: Minimum Q (Phred) score to keep (default 10)
    - `p_value`: Minimum % of bases that must have a Q (Phred) quality
    (default 50)
- **Output**
  - Reads file, filtered (`.fastq`); used in
  [**fastq_to_fasta**](#fastq_to_fasta)


#### `fastq_to_fasta`

Convert reads file from `.fastq` to `.fa` with 
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file, filtered (`.fastq`); from
  [**fastq_quality_filter**](#fastq-quality-filter)
- **Output**
  - Reads file (`.fa`); used in [**format_fasta**](#format_fasta)

#### `format_fasta`

Format reads to appear on a single line with
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**start**](#start) or 
  [**fastq_to_fasta**](#fastq_to_fasta)
- **Output**
  - Reads file, formatted (`.fa`); used in
  [**remove_adapters**](#remove_adapters)


#### `remove_adapters`

Trim adapters and `N` bases at either end. Filter reads by minimum length and
number of inner `N` bases with [**cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**format_fasta**](#format_fasta)
- **Parameters**
  - **samples.tsv**
    - Adapter to be removed; specify in sample table column `adapter`
  - **config_template.yaml**
    - `error_rate`: Fraction of allowed errors in the matched adapters
    (default 0.1)
    - `overlap`: Minimum overlap length between adapter and read to trim the
    bases (default 3)
    - `minimum_length`: Minimum length for a processed read to be kept
    (default 15)
    - `max_n`:  Maximum number of `N` bases for a processed read to be kept
    (default 0)
- **Output**
  - Reads file (`.fa`; used in 
  [**collapse_identical_reads**](#collapse_identical_reads)


#### `collapse_identical_reads`

Collapse and rename identical reads
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**remove_adapters**](#remove_adapters)
- **Output**
  - Reads file, collapsed, rename; used in
  [**filter_fasta_for_oligomap**](#filter_fasta_for_oligomap),
  [**map_genome_segemehl**](#map_genome_segemehl) and
  [**map_transcriptome_segemehl**](#map_transcriptome_segemehl)


#### `map_genome_segemehl`

Align short reads to reference genome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**collapse_identical_reads**](#collapse_identical_reads)
  - Genome sequence (`.fa`); from 
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
  - segemehl genome index (`idx`); from
  [**generate_segemehl_index_genome**](#generate_segemehl_index_genome)
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_genome_maps**](#merge_genome_maps)


#### `map_transcriptome_segemehl`

Align short reads to reference transcriptome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**collapse_identical_reads**](#collapse_identical_reads)
  - Transcriptome sequence (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)
  - segemehl transcriptome index (`idx`); from
  [**generate_segemehl_index_transcriptome**](#generate_segemehl_index_transcriptome)
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_transcriptome_maps**](#merge_transcriptome_maps)


#### `filter_fasta_for_oligomap`

Filter reads by length with a [**custom script**][custom-script-validation].

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**collapse_identical_reads**](#collapse_identical_reads)
- **Parameters**
  - **config_template.yaml**
    - `max_length_reads`: Maximum length of processed reads to map with
    [**oligomap**](#third-party-software-used)
- **Output**
  - Reads file (`.fa`); used in [**map_genome_oligomap**](#map_genome_oligomap)
  and [**map_transcriptome_oligomap**](#map_transcriptome_oligomap)


#### `map_genome_oligomap`

Align short reads to reference genome with
[**oligomap**](#third-party-software-used).

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**filter_fasta_for_oligomap**](#filter_fasta_for_oligomap)
  - Genome sequence (`.fa`); from 
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - Alignments file (`.fa`); used in
  [**sort_genome_oligomap**](#sort_genome_oligomap)
  - Alignment report (`.txt`); used in
  [**sort_genome_oligomap**](#sort_genome_oligomap)


#### `sort_genome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

- **Input**
  - Alignments file (`.fa`); from
  [**map_genome_oligomap**](#map_genome_oligomap)
  - Alignment report (`.txt`); from
  [**map_genome_oligomap**](#map_genome_oligomap)
- **Output**
  - Alignments file, sorted (`.fa`); used in
  [**convert_genome_to_sam_oligomap**](#convert_genome_to_sam_oligomap)
  - Alignment report, sorted (`.txt`); used in
  [**convert_genome_to_sam_oligomap**](#convert_genome_to_sam_oligomap)


#### `convert_genome_to_sam_oligomap`

Convert aligned reads `.fa` to `.sam` and filter alignments by number of hits
with a [**custom script**][custom-script-oligo-sam].

- **Input**
  - Alignments file, sorted (`.fa`); from
  [**sort_genome_oligomap**](#sort_genome_oligomap)
  - Alignment report, sorted (`.txt`); from
  [**sort_genome_oligomap**](#sort_genome_oligomap)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of hits an alignment can have to be kept
- **Output**
  - Alignments file (`.sam`); used in [**merge_genome_maps**](#merge_genome_maps)


#### `map_transcriptome_oligomap`

Align short reads to reference transcriptome with
[**oligomap**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from
  [**filter_fasta_for_oligomap**](#filter_fasta_for_oligomap)
  - Transcriptome sequence (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)
- **Output**
  - Alignments file (`.fa`); used in
  [**sort_transcriptome_oligomap**](#sort_transcriptome_oligomap)
  - Alignment report (`.txt`); used in
  [**sort_transcriptome_oligomap**](#sort_transcriptome_oligomap)


#### `sort_transcriptome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

- **Input**
  - Alignments file (`.fa`); from
  [**map_transcriptome_oligomap**](#map_transcriptome_oligomap)
  - Alignment report (`.txt`); from
  [**map_transcriptome_oligomap**](#map_transcriptome_oligomap)
- **Output**
  - Aligned file, sorted (`.fa`); used in
  [**convert_transcriptome_to_sam_oligomap**](#convert_transcriptome_to_sam_oligomap)
  - Alignment report, sorted (`.txt`); used in
  [**convert_transcriptome_to_sam_oligomap**](#convert_transcriptme_to_sam_oligomap)


#### `convert_transcriptome_to_sam_oligomap`

Convert aligned reads `.fa` to `.sam` and filter alignments by number of hits
with a [**custom script**][custom-script-oligo-sam].

- **Input**
  - Alignments file, sorted (`.fa`); from
  [**sort_transcriptome_oligomap**](#sort_transcriptome_oligomap)
  - Alignment report, sorted (`.txt`); from
  [**sort_transcriptome_oligomap**](#sort_transcriptome_oligomap)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of hits an alignment can have to be kept
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_transcriptome_maps**](#merge_transcriptome_maps)


#### `merge_genome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) genome alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**map_genome_segemehl**](#map_genome_segemehl)
  - Alignments file (`.sam`); from
  [**convert_genome_to_sam_oligomap**](#convert_genome_to_sam_oligomap)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_genome_by_nh**](#filter_genome_by_nh)


#### `merge_transcriptome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) transcriptome alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**map_transcriptome_segemehl**](#map_transcriptome_segemehl)
  - Alignments file (`.sam`); from
  [**convert_transcriptome_to_sam_oligomap**](#convert_transcriptome_to_sam_oligomap)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_transcriptome_by_nh**](#filter_transcriptome_by_nh)


#### `filter_genome_by_nh`

Filter merged genome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

- **Input**
  - Alignments file (`.sam`); from
  [**merge_genome_maps**](#merge_genome_maps)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**remove_header_genome_mappings**](#remove_header_genome_mappings)

#### `filter_transcriptome_by_nh`

Filter merged transcriptome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

- **Input**
  - Alignments file (`.sam`); from
  [**merge_transcriptome_maps**](#merge_transcriptme_maps)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**remove_header_transcriptome_mappings**](#remove_header_transcriptome_mappings)


#### `remove_header_genome_mappings`

Remove the `SAM` header of the genome alignments file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_genome_by_nh**](#filter_genome_by_nh)
- **Output**
  - Alignments file (`.sam`); used in [**merge_all_maps**](#merge_all_maps)


#### `remove_header_transcriptome_mappings`

Remove the `SAM` header of the transcriptome alignments file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_transcriptome_by_nh**](#filter_transcriptome_by_nh)
- **Output**
  - Alignments file (`.sam`); used in 
  [**transcriptome_to_genome_maps**](#transcriptome_to_genome_maps)


#### `transcriptome_to_genome_maps`

Convert the alignments transcriptome coordinates to genomic ones with a
[**custom script**][custom-script-sam-trx].

- **Input**
  - Alignments file (`.sam`); from 
  [**remove_header_transcriptome_mappings**](#remove_header_transcriptome_mappings)
  - Exon annotations (`.bed`); from
  [**convert_exons_gtf_to_bed**](#convert_exons_gtf_to_bed)
- **Output**
  - Alignments file (`.sam`); used in [**merge_all_maps**](#merge_all_maps)


#### `merge_all_maps`

Concatenate the four alignments files into a single file.

- **Input**
  - Alignments file (`.sam`); from
  [**remove_header_genome_mappings**](#remove_header_genome_mappings) and
  [**transcriptome_to_genome_maps**](#transcriptome_to_genome_maps)
- **Output**
  - Alignments file (`.sam`); used in
  [**add_header_all_maps**](#add_header_all_maps)


#### `add_header_all_maps`

Add the `SAM` header to the aligned reads with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from [**merge_all_maps**](#merge_all_maps)
- **Output**
  - Alignments file (`.sam`); used in [**sort_maps_by_id**](#sort_maps_by_id)


#### `sort_maps_by_id`

Sort alignments by reads ID with [**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**add_header_all_maps**](#add_header_all_maps)
- **Output**
  - Alignments file, sorted (`.sam`); used in
  [**remove_inferiors**](#remove_inferiors)


#### `remove_inferiors`

Remove duplicate and inferior alignments with a
[**custom script**][custom-script-remove-dup].

> Alignments are considered to be duplicates if having identical entries for
the fields `QNAME`, `FLAG`, `RNAME`, `POS` and `CIGAR`. 
Alignments are considered to be inferiors if having the same `QNAME` and
a bigger edit distance than the minimum one within the group.


- **Input**
  - Alignments file, sorted (`.sam`); from [**sort_maps_by_id**](#sort_maps_by_id)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_by_indels**](#filter_by_indels)


#### `filter_by_indels`

Remove multimappers favoring mismatches over indels with a
[**custom script**][custom-script-filter-mm].

- **Input**
  - Alignments file, sorted (`.sam`); from
  [**remove_inferiors**](#remove_inferiors)
- **Output**
  - Alignments file (`.sam`); used in
  [**convert_all_alns_sam_to_bam**](#convert_all_alns_sam_to_bam) and
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primr)


#### `convert_all_alns_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with pri-miR annotations.

- **Input**
  - Alignments file (`.sam`); from [**filter_by_indels**](#filter_by_indels)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_all_alns_bam_by_position**](#sort_all_alns_bam_by_position)


#### `sort_all_alns_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with pri-miR annotations more efficiently.

- **Input**
  - Alignments file (`.bam`); from
  [**convert_all_alns_sam_to_bam**](#convert_all_alns_sam_to_bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_all_alns_bam**](#index_all_alns_bam) and
  [**intersect_extended_primir**](#intersect_extended_primir)


#### `index_all_alns_bam`

Create index `BAM` file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
alignments in a genomic region of interest.

- **Input**
  - Alignments file (`.bam`); from
  [**sort_all_alns_bam_by_position**](#sort_all_alns_bam_by_position)
- **Output**
  - `BAM` index file (`.bam.bai`); used in
  [**intersect_extended_primir**](#intersect_extended_primir)


### Quantify workflow

#### `finish_quantify`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir)
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  - Alignments file, sorted, tagged (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort_intersecting_mirna_by_feat_tag)
  - (iso)miR and/or pri-miR counts table (`.tab`); from
  [**merge_tables**](#merge_tables)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)

#### `intersect_extendend_primir`

Intersect the aligned reads with the extended pri-miR annotations with
[**BEDTools**](#third-party-software-used).

- **Input**
  - Alignments file (`.bam`); from
  [**sort_all_alns_bam_by_position**](#sort_all_alns_bam_by_position)
  - pri-miR extended annotations (`.gff3`); from
  [**extend_mirs_annotations**](#extend_mirs_annotations)
- **Output**
  - pri-miR intersections file (`.bed`); used in
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir)
  and [**quantify_primir**](#quantify_primir)


#### `filter_sam_by_intersecting_primir`

Remove alignments that do not intersect with any pri-miR with
[**SAMtools**](#third-party-software-used).

> Required to only intersect alignments within a pri-miR locus.

- **Input**
  - Alignments file (`.sam`); from [**filter_by_indels**](#filter_by_indels) 
  - pri-miR intersections file (`.bed`); from
  [**intersect_extended_primir**](#intersect_extended_primir)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**convert_intersecting_primir_sam_to_bam**](#convert_intersecting_primir_sam_to_bam)
  and [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)


#### `convert_intersecting_primir_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with miRNA annotations.

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_intersecting_primir_bam_by_position**](#sort_intersecting_primir_bam_by_position)


#### `sort_intersecting_primir_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with miRNA annotations more efficiently.

- **Input**
  - Alignments file (`.bam`); from
  [**convert_intersecting_primir_sam_to_bam**](#convert_intersecting_primir_sam_to_bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_intersecting_primir_bam**](#index_intersecting_primir_bam) and
  [**intersect_extended_mirna**](#intersect_extended_mirna)


#### `index_intersecting_primir_bam`

Create index `BAM` file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
alignments in a genomic region of interest.

- **Input**
  - Alignments file (`.bam`); from
  [**sort_intersecting_primir_bam_by_position**](#sort_intersecting_primir_bam_by_position)
- **Output**
  - `BAM` index file (`.bam.bai`); used in
  [**intersect_extended_mirna**](#intersect_extended_mirna)


#### `intersect_extended_mirna`

Intersect the aligned reads with the extended miRNA annotations with
[**BEDTools**](#third-party-software-used).

- **Input**
  - Alignments file (`.bam`); from
  [**sort_intersecting_primir_bam_by_position**](#sort_intersecting_primir_bam_by_position)
  - miRNA extended annotations (`.gff3`); from
  [**extend_mirs_annotations**](#extend_mirs_annotations)
- **Output**
  - miRNA intersections file (`.bed`); used in
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  and [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag)


#### `filter_sam_by_intersecting_mirna`

Remove alignments that do not intersect with any miRNA with
[**SAMtools**](#third-party-software-used).

> Required to efficiently classify the alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir) 
  - miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect_extended_mirna)
- **Output**
  - Alignments file, filtered (`.sam`); used in 
  [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag) and
  [**uncollapse_reads**](#uncollapse_reads)


#### `add_intersecting_mirna_tag`

Classify and add the intersecting (iso)miR to each alignment as a tag
with a [**custom script**][custom-script-iso-tag].

> To classify the reads, miRNA annotations are turned into the original ones.
The alignment start and end shifts relative to the annotations are computed
and append along with the `CIGAR` and `MD` strings to the intersecting miRNA.
If the read is classified as a canonical miRNA, the name will not include
the star and end shift, nor the `CIGAR` and `MD` strings.

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  - miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect_extended_mirna)
- **Parameters**
  - **config_template.yaml**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended (default 6)
- **Output**
  - Alignments file, tagged (`.sam`); used in
  [**sort_intersecting_mirna_by_feat_tag**](#sort_intersecting_mirna_by_feat_tag)


#### `sort_intersecting_mirna_by_feat_tag`

Sort the alignments by the tag containing the classified intersecting miRNA
with [**SAMtools**](#third-party-software-used).

> Required for an efficient quantification.

- **Input**
  - Alignments file, tagged (`.sam`); from
  [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag)
- **Output**
  - Alignments file, tagged, sorted (`.sam`); used in
  [**quantify_mirna**](#quantify_mirna)


#### `quantify_mirna`

Tabulate alignments according to its tag with a
[**custom script**][custom-script-mir-quant].

> Quantification is done with partial counts (_i.e._ each alignment contributes
by the number of collapsed reads divided by the number of hits).

- **Input**
  - Alignments file, sorted, tagged (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort_intersecting_mirna_by_feat_tag)
- **Parameters**
  - **samples.tsv**
    - Library name; specify in sample table column `sample`
  - **config_template.yaml**
    - `mir_list`: miRNA features to be quantified (default isomir, mirna
    pri-miR)
- **Output**
  - (iso)miR counts tab-delimited file; used in
  [**merge_tables**](#merge_tables)


#### `quantify_primir`

Tabulate alignments according to its intersecting pri-miR with a
[**custom script**][custom-script-pri-quant]

> Quantification is done with partial counts (_i.e._ each alignment contributes
by the number of collapsed reads divided by the number of hits).

- **Input**
  - pri-miR intersections file (`.bed`); from
  [**intersecti_extended_primir**](#intersect_extended_primir)
- **Output**
  - pri-miR counts tab-delimited file; used in
  [**merge_tables**](#merge_tables)


#### `merge_tables`

Merge all the tables from the different libraries into a single one with a
[**custom script**][custom-script-merge-tab].

- **Input**
  - counts tab-delimited file; from [**quantify_mirna**](#quantify_mirna)
  and/or [**quantify_primir**](#quantify_primir)
- **Parameters**
  - **cluster_schema.json**
    - `mir_list`: miRNA features to be quantified (default isomir, mirna
    pri-mir)
- **Output**
  - (iso)miR and/or pri-miR counts table (`.tab`)


#### `uncollapse_reads`

Reverse the collapsing of reads with identical sequences as done with
[**FASTX-Toolkit**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
- **Output**
  - Uncollapsed aligned reads (`.sam`); used in
  [**convert_uncollapsed_reads_sam_to_bam**](#convert_uncollapsed_reads_sam_to_bam)


#### `convert_uncollapsed_reads_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)


#### `sort_uncollapsed_reads_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.bam`); from
  [**convert_uncollapsed_reads_sam_to_bam**](#convert_uncollapsed_reads_sam_to_bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam),
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `index_uncollapsed_reads_bam`

Create index `BAM` file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
alignments in a genomic region of interest.

- **Input**
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
- **Output**
  - `BAM` index file (`.bam.bai`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


### Pileup workflow

#### `finish_pileup`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Empty text files (`.txt`); from
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups) and
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups)


#### `create_empty_bed`

Create an empty BED file if the user has not provided one.

> **OPTIONAL RULE.** This rule will be executed if, and only if, the user has
> not provided a BED file in the [configuration file](#configuration-file)
> with the regions the ASCII-style alignment pileups must be performed on.

- **Condition**
  - **config_template.yaml**
    - `bed_file`: BED6 file with all the desired annotation regions to perform
    the ASCII-style alignment pileups on. (Default: None)
- **Output**
  - BED empty file (`.bed`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `compress_reference_genome`

Compress the processed genome with trimmed IDs using `bgzip`. 

> Required to perform the ASCII-style alignment pileups.
> In order to be able to use the `bgzip` command when running the workflow
> using Singularity, [**SAMtools**](#third-party-software-used) is used.

- **Input**
  - Genome sequence file, trimmed IDs (`.fa`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - Genome sequence file, trimmed IDs, `bgzip`ed (`.fa.bz`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `create_per_library_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions across
libraries with [**ASCII-style alignment pilueups**](#third-party-software-used).

> A directory containing the ASCII-style pileups is created for each
> library. If no BED file is provided, the pileups' output directories will
> only contain an empty file.

- **Input**
  - Genome sequence file, trimmed IDs, `bgzip`ed (`.fa.bz`); from
  [**compress_reference_genome**](#compress_reference_genome)
  - miRNA annotations with mapped genes(`.gff3`); from
  [**map_chr_names**](#map_chr_names)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Annotated genomic regions (`.bed`); from workflow input files or
  [**create_empty_bed**](#create_empty_bed)
- **Output**
  - Empty text file (`.txt`)


#### `create_per_run_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions for the whole
run with [**ASCII-style alignment pilueups**](#third-party-software-used). 

> If no BED file is provided, the pileups' output directory will only contain
> an empty file.

- **Input**
  - Genome sequence file, trimmed IDs, `bgzip`ed (`.fa.bz`); from
  [**compress_reference_genome**](#compress_reference_genome)
  - miRNA annotations with mapped genes(`.gff3`); from
  [**map_chr_names**](#map_chr_names)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Annotated genomic regions (`.bed`); from workflow input files or
  [**create_empty_bed**](#create_empty_bed)
- **Output**
  - Empty text file (`.txt`)


#### `create_per_condition_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions across the
different library subsets if provided with
[**ASCII-style alignment pilueups**](#third-party-software-used).

> **OPTIONAL RULE.** The ASCII-style pileups for each annotated region are
> made if, and only if, at least one library subset is specified in the
> [configuration file](#configuration-file). Otherwise, this rule will not be
> executed, and no output will be generated.

- **Input**
  - Genome sequence file, trimmed IDs, `bgzip`ed (`.fa.bz`); from
  [**compress_reference_genome**](#compress_reference_genome)
  - miRNA annotations with mapped genes(`.gff3`); from
  [**map_chr_names**](#map_chr_names)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Annotated genomic regions (`.bed`); from workflow input files or
  [**create_empty_bed**](#create_empty_bed)
- **Parameters**
  - **config_template.yaml**
    - `lib_dict`: Subset(s) of library name(s), as specified in the samples'
    table column `sample` and the subset identifier stored in a dictionary.
    (default: None)
- **Output**
  - Empty text file (`.txt`)


[chr-maps]: <https://github.com/dpryan79/ChromosomeMappings>
[custom-script-blocksort]: scripts/blocksort.sh
[custom-script-filter-mm]: scripts/filter_multimappers.py
[custom-script-get-lines]: scripts/get_lines_w_pattern.sh
[custom-script-gtf-bed]: scripts/gtf_exons_bed.1.1.2.R
[custom-script-iso-tag]: scripts/iso_name_tagging.py
[custom-script-map-chr]: scripts/map_chromosomes.pl
[custom-script-merge-tab]: scripts/merge_tables.R
[custom-script-mir-ext]: scripts/mirna_extension.py
[custom-script-mir-quant]: scripts/mirna_quantification.py
[custom-script-nh-filter]: scripts/nh_filter.py
[custom-script-oligo-sam]: scripts/oligomap_output_to_sam_nh_filtered.py
[custom-script-pri-quant]: scripts/primir_quantification.py
[custom-script-remove-dup]: scripts/sam_remove_duplicates_inferior_alignments_multimappers.pl
[custom-script-sam-trx]: scripts/sam_trx_to_sam_gen.pl
[custom-script-trim-id]: scripts/trim_id_fasta.sh
[custom-script-uncollapse]: scripts/sam_uncollapse.pl
[custom-script-validation]: scripts/validation_fasta.py
[code-ascii]: <https://git.scicore.unibas.ch/zavolan_group/tools/ascii-alignment-pileup>
[code-bedtools]: <https://github.com/arq5x/bedtools2>
[code-cufflinks]: <https://github.com/cole-trapnell-lab/cufflinks>
[code-cutadapt]: <https://github.com/marcelm/cutadapt>
[code-fastx]: <https://github.com/agordon/fastx_toolkit>
[code-gffutils]: <https://github.com/fls-bioinformatics-core/GFFUtils>
[code-oligomap]: <https://github.com/zavolanlab/oligomap>
[code-samtools]: <https://github.com/samtools/samtools>
[docs-bedtools]: <https://bedtools.readthedocs.io/en/latest/>
[docs-cufflinks]: <http://cole-trapnell-lab.github.io/cufflinks/manual/>
[docs-cutadapt]: <https://cutadapt.readthedocs.io/en/stable/>
[docs-fastx]: <http://hannonlab.cshl.edu/fastx_toolkit/commandline.html>
[docs-gffutils]: <https://gffutils.readthedocs.io/en/latest/>
[docs-samtools]: <http://www.htslib.org/doc/samtools.html>
[docs-segemehl]: <http://www.bioinf.uni-leipzig.de/Software/segemehl/>
[docs-snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[license-afl3]: <https://opensource.org/license/afl-3-0-php/>
[license-agpl3]: <https://opensource.org/license/agpl-v3/>
[license-apache2]: <https://opensource.org/license/apache-2-0/>
[license-bsl1]: <https://opensource.org/license/bsl-1-0/>
[license-gpl2]: <https://opensource.org/licenses/GPL-2.0>
[license-gpl3]: <https://opensource.org/license/gpl-3-0/>
[license-mit]: <https://opensource.org/licenses/MIT>
[pub-bedtools]: <https://academic.oup.com/bioinformatics/article/26/6/841/244688>
[pub-cufflinks]: <https://doi.org/10.1038/nprot.2012.016>
[pub-cutadapt]: <https://doi.org/10.14806/ej.17.1.200>
[pub-oligomap]: <https://doi.org/10.1016/j.ymeth.2007.10.002>
[pub-samtools]: <https://doi.org/10.1093/bioinformatics/btp352>
[pub-segemehl]: <https://doi.org/10.1371/journal.pcbi.1000502>
[rule-graph]: images/rule_graph.svg

