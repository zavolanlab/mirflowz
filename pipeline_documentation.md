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
    - [`finish_prepare`](#finish-prepare)
    - [`trim_genome_seq_ids`](#trim-genome-seq-ids)
    - [`extract_transcriptome_seqs`](#extract-transcriptome-seqs)
    - [`trim_transcriptome_seq_ids`](#trim-transcriptome-seq-ids)
    - [`generate_segemehl_index_transcriptome`](#generate-segemehl-index-transcriptome)
    - [`generate_segemehl_index_genome`](#generate-segemehl-index-genome)
    - [`get_exons_gtf`](#get-exons-gtf)
    - [`convert_exons_gtf_to_bed`](#convert-exons-gtf-to-bed)
    - [`create_genome_header`](#create-genome-header)
    - [`map_chr_names`](#create-chr-names)
    - [`create_index_genome_fasta`](#create-index-genome-fasta)
    - [`extract_chr_len`](#extract-chr-len)
    - [`extend_mirs_annotations`](#extend-mirs-annotations)
  - [Map workflow](#map-workflow)
    - [`finish_map`](#finish-map)
    - [`start`](#start)
    - [`fastq_quality_filter`](#fastq-quality-filter)
    - [`fastq_to_fasta`](#fastq-to-fasta)
    - [`format_fasta`](#format-fasta)
    - [`remove_adapters`](#remove-adapters)
    - [`collapse_indentical_reads`](#collapse-indentical-reads)
    - [`map_genome_segemehl`](#map-genome-segemehl)
    - [`map_transcriptome_segemehl`](#map-transcriptome-segemehl)
    - [`filter_fasta_for_oligomap`](#filter-fasta-for-oligomap)
    - [`map_genome_oligomap`](#map-genome-oligomap)
    - [`sort_genome_oligomap`](#sort-genome-oligomap)
    - [`convert_genome_to_sam_oligomap`](#convert-genome-to-sam-oligomap)
    - [`map_transcriptome_oligomap`](#map-transcriptome-oligomap)
    - [`sort_transcriptome_oligomap`](#sort-transcriptome-oligomap)
    - [`convert_transcriptome_to_sam_oligomap`](#convert-transcriptome-to-sam-oligomap)
    - [`merge_genome_maps`](#merge-genome-maps)
    - [`merge_transcriptome_maps`](#merge-transcriptome-maps)
    - [`filter_genome_by_nh`](#filter-genome-by-nh)
    - [`filter_transcriptome_by_nh`](#filter-transcriptome-by-nh)
    - [`remove_header_genome_mappings`](#remove-header-genome-mappings)
    - [`remove_header_transcriptome_mappings`](#remove-header-transcriptome-mappings)
    - [`transcriptome_to_genome_maps`](#transcriptome-to-genome-maps)
    - [`merge_all_maps`](#merge-all-maps)
    - [`add_header_all_maps`](#add-header-all-maps)
    - [`sort_maps_by_id`](#sort-maps-by-id)
    - [`remove_inferiors`](#remove-inferiors)
    - [`filter_by_indels`](#filter-by-indels)
    - [`convert_all_alns_sam_to_bam`](#convert-all-alns-sam-to-bam)
    - [`sort_all_alns_bam_by_position`](#sort-all-alns-bam-by-position)
    - [`index_all_alns_bam`](#index-all-alns-bam)
  - [Quantify workflow](#quantify-workflow)
    - [`finish_quantify`](#finish-quantify)
    - [`intersect_extended_primir`](#intersect-extended-primir)
    - [`filter_sam_by_intersecting_primir`](#filter-sam-by-intersecting-primir)
    - [`convert_intersecting_primir_sam_to_bam`](#convert-intersecting-primir-sam-to-bam)
    - [`sort_intersecting_primir_bam_by_position`](#sort-intersecting-primir-bam-by-position)
    - [`index_intersecting_primir_bam`](#index-intersecting-primir-bam)
    - [`intersect_extended_mirna`](#intersect-extended-mirna)
    - [`filter_sam_by_intersecting_mirna`](#filter-sam-by-intersecting-mirna)
    - [`add_intersecting_mirna_tag`](#add-intersecting-mirna-tag)
    - [`sort_intersecting_mirna_by_feat_tag`](#sort-intersecting-mirna-by-feat-tag)
    - [`quantify_mirna`](#quantify-mirna)
    - [`quantify_primir`](#quantify-primir)
    - [`merge_tables`](#merge-tables)
    - [`uncollapse_reads`](#uncollapse-reads)
    - [`convert_uncollapse_reads_sam_to_bam`](#convert-uncollapse-reads-sam-to-bam)
    - [`sort_uncollapse_reads_bam_by_position`](#sort-uncollapse-reads-bam-by-position)
    - [`index_uncollapse_reads_bam`](#index-uncollapse-reads-bam)



## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
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
the reads mapping and the miRNA quantification). The main `Snakefile` contains
the configuration file validation along with the inclusion of the
sub-workflows. Individual steps of the workflow are described briefly, and
links to the respective software manuals are given. Parameters that can be
modified by the user (via the samples table and the configuration file) are
also described.

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
  [**intersect_extended_primir**](#intersect-extended-primir)
  - miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect-extended-mirna)
  - Alignments file, sorted and tagged (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort-intersecting-mirna-by-feat-tag)
  - (iso)miR and/or pri-miR counts table (`.tab`); from
  [**merge_tables**](#merge-tables)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort-uncollapsed-reads-bam-by-position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index-uncollapsed-reads-bam)


### Prepare workflow

#### `finish_prepare`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - segemehl genome index (`idx`); from
  [**generate_segemehl_index_genome**](#generate-segemehl-index-genome)
  - segemehl transcriptome index (`idx`); from
  [**generate_segemehl_index_transcriptome**](#generate-segemehl-index-transcriptome)
  - Exon annotations (`.bed`); from
  [**convert_exons_gtf_to_bed**](#convert-exons-gtf-to-bed)
  - Genome header (`.sam`); from
  [**create_genome_header**](#create-genome-header)
  - Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`);
  from [**extract_chr_len**](#extract-chr-len)
  - Primary miRNA transcript (pri-miR) extended annotation (`.gff3`);
  from [**extend_mirs_annotations**](#extend-mirs-annotations)
  - Mature miRNA (miRNA) extended annotation (`.gff3`);
  from [**extend_mirs_annotations**](#extend-mirs-annotations)


#### `trim_genome_seq_ids`

Trim genome sequence IDs with a [**custom script**][custom-script-trim-id].

- **Input**
  - Genome sequence file (`.fasta`)
- **Output**
  - Genome sequence file, trimmed IDs (`.fasta`); used in
  [**extract_transcriptome_seqs**](#extract-transcriptome-seqs),
  [**create_genome_header**](#create-genome-header),
  [**create_index_genome_fasta**](#create-index-genome),
  [**generate_segemehl_index_genome**](#generate-segemehl-index-genome),
  [**mapping_genome_segemehl**](#mapping-genome-segemehl) and
  [**mapping_genome_oligomap**](#mapping-genome-oligomap)


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
  [**generate_segemehl_index_transcriptome**](#generate-segemehl-index-transcriptome),
  [**mapping_transcriptome_segemehl**](#mapping-transcriptome-segemehl) and
  [**mapping_transcriptome_oligomap**](#mapping-transcriptome-oligomap)


#### `generate_segemehl_index_transcriptome`

Generate transcriptome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The transcriptome index only needs to be generated once for each combination
of genome and annotations and sample sets.

- **Input**
  - Transcriptome sequence file, trimmed IDs (`.fasta`); from
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)
- **Output**
  - segemehl transcriptome index (`.idx`); used in 
  [**mapping_genome_segemehl**](#mapping-genome-segemehl)


#### `generate_segemehl_index_genome`

Generate genome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The genome index only needs to be generated once for each combination
of annotations and sample sets.

- **Input**
  - Genome sequence file with trim IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - segemehl genome index (`.idx`); used in 
  [**mapping_genome_segemehl**](#mapping-genome-segemehl)


#### `get_exons_gtf`

Retrieve exon annotations from genome annotations with a
[**custom script**][custom-script-get-lines].

- **Input**
  - Genomic annotations (`.gtf`)
- **Output**
  - Exon annotations (`.gtf`); used in 
  [**convert_exons_gtf_to_bed**](#convert-exons-gtf-to-bed)

#### `convert_exons_gtf_to_bed`

Convert exon annotations `.gtf` to `.bed` with a
[**custom script**][custom-script-gtf-bed].

- **Input**
  - Exon annotations (`.gtf`); from [**get_exons_gtf**](#get-exons-gtf)
- **Output**
  - Exon annotations (`.bed`); used in
  [**transcriptome_to_genome_maps**](#transcriptome-to-genome-maps)


#### `create_genome_header`

Create `SAM` header for the genome with 
[**SAMtools**](#third-party-software-used).

> Required by [SAMtools](#third-party-software-used) to work with the alignment
file.

- **Input**
  - Genome sequence file, trimmed IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - Genome header (`.sam`); used in [add_header_all_maps](#add-header-all-maps)


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
  [**extend_mirs_annotations**](#extend-mirs-annotations)


#### `create_index_genome_fasta`

Create a `FASTA` index for the genome with 
[**SAMtools**](#third-party-software-used).

- **Input**
  - Genome sequence file, trimmed IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - `FASTA` genome index (`.fa.fai`); used in
  [**extract_chr_len**](#extract-chr-len)


#### `extract_chr_len`

Extract chromosome(s) length from the genome sequence.

- **Input**
  - `FASTA` genome index (`.fa.fai`); from
  [**create_index_genome_fasta**](#create-index-genome-fasta)
- **Output**
  - Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`); used
  in [**extend_mirs_annotations**](#extend-mirs-annotations)


#### `extend_mirs_annotations`

Extend miRNA annotations and split the file by feature with a
[**custom script**][custom-script-mir-ext].

> Mature miRNA regions are extended on both sides to account for isomiR species
with shifted start and/or end positions. If required, pri-miR loci are also 
extended to accommodate the new miRNA coordinates. 

- **Input**
  - miRNA annotations with mapped chromosomes(`.gff3`); from
  [**map_chr_names**](#map-chr-names)
- **Parameters**
  - **config_template.yaml**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended (default 6)
- **Output**
  - Primary miRNA transcript (pri-miR) extended annotation (`.gff3`); used in
  [**intersect_extended_primir**](#intersect-extended-primir)
  - Mature miRNA (miRNA) extended annotation (`.gff3`); used in
  [**intersect_extended_mirna**](#intersect-extended-mirna)

### Map workflow

#### `finish_map`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - `BAM` index file (`.bam.bai`); from
  [**index_all_alns_bam**](#index-all-alns-bam)



#### `start`

Copy and rename read files.

> Local rule.
Depending on the read files format, the output files undergo a quality filter
(`.fastq`) or are directly formatted (`.fa`).

- **Input**
  - Reads file (`.fa.gz`, `.fastq.gz`)
- **Output**
  - Reads file, copied, renamed (`.fa`, `.fastq`); used in
  [**fastq_quality_filter**](#fastq-quality-filter) or 
  [**format_fasta**](#format-fasta)


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
  [**fastq_to_fasta**](#fastq-to-fasta)


#### `fastq_to_fasta`

Convert reads file from `.fastq` to `.fa` with 
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file, filtered (`.fastq`); from
  [**fastq_quality_filter**](#fastq-quality-filter)
- **Output**
  - Reads file (`.fa`); used in [**format_fasta**](#format-fasta)

#### `format_fasta`

Format reads to appear on a single line with
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**start**](#start) or 
  [**fastq_to_fasta**](#fastq-to-fasta)
- **Output**
  - Reads file, formatted (`.fa`); used in
  [**remove_adapters**](#remove-adapters)


#### `remove_adapters`

Trim adapters and `N` bases at either end. Filter reads by minimum length and
number of inner `N` bases with [**cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**format_fasta**](#format-fasta)
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
  [**collapse_identical_reads**](#collapse-identical-reads)


#### `collapse_identical_reads`

Collapse and rename identical reads
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**remove_adapters**](#remove-adapters)
- **Output**
  - Reads file, collapsed, rename; used in
  [**filter_fasta_for_oligomap**](#filter-fasta-for-oligomap),
  [**map_genome_segemehl**](#map-genome-segemehl) and
  [**map_transcriptome_segemehl**](#map-transcriptome-segemehl)


#### `map_genome_segemehl`

Align short reads to reference genome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**collapse_identical_reads**](#collapse-identical-reads)
  - Genome sequence (`.fa`); from 
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
  - segemehl genome index (`idx`); from
  [**generate_segemehl_index_genome**](#generate-segemehl-index-genome)
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_genome_maps**](#merge-genome-maps)


#### `map_transcriptome_segemehl`

Align short reads to reference transcriptome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**collapse_identical_reads**](#collapse-identical-reads)
  - Transcriptome sequence (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)
  - segemehl transcriptome index (`idx`); from
  [**generate_segemehl_index_transcriptome**](#generate-segemehl-index-transcriptome)
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_transcriptome_maps**](#merge-transcriptome-maps)


#### `filter_fasta_for_oligomap`

Filter reads by length with a [**custom script**][custom-script-validation].

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**collapse_identical_reads**](#collapse-identical-reads)
- **Parameters**
  - **config_template.yaml**
    - `max_length_reads`: Maximum length of processed reads to map with
    [**oligomap**](#third-party-software-used)
- **Output**
  - Reads file (`.fa`); used in [**map_genome_oligomap**](#map-genome-oligomap)
  and [**map_transcriptome_oligomap**](#map-transcriptome-oligomap)


#### `map_genome_oligomap`

Align short reads to reference genome with
[**oligomap**](#third-party-software.used).

- **Input**
  - Reads file, collapsed (`.fa`); from
  [**filter_fasta_for_oligomap**](#filter-fasta-for-oligomap)
  - Genome sequence (`.fa`); from 
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - Alignments file (`.fa`); used in
  [**sort_genome_oligomap**](#sort-genome-oligomap)
  - Alignment report (`.txt`); used in
  [**sort_genome_oligomap**](#sort-genome-oligomap)


#### `sort_genome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

- **Input**
  - Alignments file (`.fa`); from
  [**map_genome_oligomap**](#map-genome-oligomap)
  - Alignment report (`.txt`); from
  [**map_genome_oligomap**](#map-genome-oligomap)
- **Output**
  - Alignments file, sorted (`.fa`); used in
  [**convert_genome_to_sam_oligomap**](#convert-genome-to-sam-oligomap)
  - Alignment report, sorted (`.txt`); used in
  [**convert_genome_to_sam_oligomap**](#convert-genome-to-sam-oligomap)


#### `convert_genome_to_sam_oligomap`

Convert aligned reads `.fa` to `.sam` and filter alignments by number of hits
with a [**custom script**][custom-script-oligo-sam].

- **Input**
  - Alignments file, sorted (`.fa`); from
  [**sort_genome_oligomap**](#sort-genome-oligomap)
  - Alignment report, sorted (`.txt`); from
  [**sort_genome_oligomap**](#sort-genome-oligomap)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of hits an alignment can have to be kept
- **Output**
  - Alignments file (`.sam`); used in [**merge_genome_maps**](#merge-genome-maps)


#### `map_transcriptome_oligomap`

Align short reads to reference transcriptome with
[**oligomap**](#third-party-software.used).

- **Input**
  - Reads file (`.fa`); from
  [**filter_fasta_for_oligomap**](#filter-fasta-for-oligomap)
  - Transcriptome sequence (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)
- **Output**
  - Alignments file (`.fa`); used in
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)
  - Alignment report (`.txt`); used in
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)


#### `sort_transcriptome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

- **Input**
  - Alignments file (`.fa`); from
  [**map_transcriptome_oligomap**](#map-transcriptome-oligomap)
  - Alignment report (`.txt`); from
  [**map_transcriptome_oligomap**](#map-transcriptome-oligomap)
- **Output**
  - Aligned file, sorted (`.fa`); used in
  [**convert_transcriptome_to_sam_oligomap**](#convert-transcriptome-to-sam-oligomap)
  - Alignment report, sorted (`.txt`); used in
  [**convert_transcriptome_to_sam_oligomap**](#convert-transcriptme-to-sam-oligomap)


#### `convert_transcriptome_to_sam_oligomap`

Convert aligned reads `.fa` to `.sam` and filter alignments by number of hits
with a [**custom script**][custom-script-oligo-sam].

- **Input**
  - Alignments file, sorted (`.fa`); from
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)
  - Alignment report, sorted (`.txt`); from
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of hits an alignment can have to be kept
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_transcriptome_maps**](#merge-transcriptome-maps)


#### `merge_genome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) genome alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**map_genome_segemehl**](#map-genome-segemehl)
  - Alignments file (`.sam`); from
  [**convert_genome_to_sam_oligomap**](#convert-genome-to-sam-oligomap)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_genome_by_nh**](#filter-genome-by-nh)


#### `merge_transcriptome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) transcriptome alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**map_transcriptome_segemehl**](#map-transcriptome-segemehl)
  - Alignments file (`.sam`); from
  [**convert_transcriptome_to_sam_oligomap**](#convert-transcriptome-to-sam-oligomap)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_transcriptome_by_nh**](#filter-transcriptome-by-nh)


#### `filter_genome_by_nh`

Filter merged genome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

- **Input**
  - Alignments file (`.sam`); from
  [**merge_genome_maps**](#merge-genome-maps)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**remove_header_genome_mappings**](#remove-header-genome-mappings)

#### `filter_transcriptome_by_nh`

Filter merged transcriptome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

- **Input**
  - Alignments file (`.sam`); from
  [**merge_transcriptome_maps**](#merge-transcriptme-maps)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**remove_header_transcriptome_mappings**](#remove-header-transcriptome-mappings)


#### `remove_header_genome_mappings`

Remove the `SAM` header of the genome alignments file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_genome_by_nh**](#filter-genome-by-nh)
- **Output**
  - Alignments file (`.sam`); used in [**merge_all_maps**](#merge-all-maps)


#### `remove_header_transcriptome_mappings`

Remove the `SAM` header of the transcriptome alignments file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_transcriptome_by_nh**](#filter-transcriptome-by-nh)
- **Output**
  - Alignments file (`.sam`); used in 
  [**transcriptome_to_genome_maps**](#transcriptome-to-genome--maps)


#### `transcriptome_to_genome_maps`

Convert the alignments transcriptome coordinates to genomic ones with a
[**custom script**][custom-script-sam-trx].

- **Input**
  - Alignments file (`.sam`); from 
  [**remove_header_transcriptome_mappings**](#remove-header-transcriptome-mappings)
  - Exon annotations (`.bed`); from
  [**convert_exons_gtf_to_bed**](#convert-exons-gtf-to-bed)
- **Output**
  - Alignments file (`.sam`); used in [**merge_all_maps**](#merge-all-maps)


#### `merge_all_maps`

Concatenate the four alignments files into a single file.

- **Input**
  - Alignments file (`.sam`); from
  [**remove_header_genome_mappings**](#remove-header-genome-mappings) and
  [**transcriptome_to_genome_maps**](#transcriptome-to-genome-maps)
- **Output**
  - Alignments file (`.sam`); used in
  [**add_header_all_maps**](#add-header-all-maps)


#### `add_header_all_maps`

Add the `SAM` header to the aligned reads with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from [**merge_all_maps**](#merge-all-maps)
- **Output**
  - Alignments file (`.sam`); used in [**sort_maps_by_id**](#sort-maps-by-id)


#### `sort_maps_by_id`

Sort alignments by reads ID with [**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**add_header_all_maps**](#add-header-all-maps)
- **Output**
  - Alignments file, sorted (`.sam`); used in
  [**remove_inferiors**](#remove-inferiors)


#### `remove_inferiors`

Remove duplicate and inferior alignments with a
[**custom script**][custom-script-remove-dup].

> Alignments are considered to be duplicates if having identical entries for
the fields `QNAME`, `FLAG`, `RNAME`, `POS` and `CIGAR`. 
Alignments are considered to be inferiors if having the same `QNAME` and
a bigger edit distance than the minimum one within the group.


- **Input**
  - Alignments file, sorted (`.sam`); from [**sort_maps_by_id**](#sort-maps-by-id)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_by_indels**](#filter-by-indels)


#### `filter_by_indels`

Remove multimappers favoring mismatches over indels with a
[**custom script**][custom-script-filter-mm].

- **Input**
  - Alignments file, sorted (`.sam`); from
  [**remove_inferiors**](#remove-inferiors)
- **Output**
  - Alignments file (`.sam`); used in
  [**convert_all_alns_sam_to_bam**](#convert-all-alns-sam-to-bam) and
  [**filter_sam_by_intersecting_primir**](#filter-sam-by-intersecting-primr)


#### `convert_all_alns_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with pri-miR annotations.

- **Input**
  - Alignments file (`.sam`); from [**filter_by_indels**](#filter-by-indels)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_all_alns_bam_by_position**](#sort-all-alns-bam-by-position)


#### `sort_all_alns_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with pri-miR annotations more efficiently.

- **Input**
  - Alignments file (`.bam`); from
  [**convert_all_alns_sam_to_bam**](#convert-all-alns-sam-to-bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_all_alns_bam**](#index-all-alns-bam) and
  [**intersect_extended_primir**](#intersect-extended-primir)


#### `index_all_alns_bam`

Create index `BAM` file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
alignments in a genomic region of interest.

- **Input**
  - Alignments file (`.bam`); from
  [**sort_all_alns_bam_by_position**](#sort-all-alns-bam-by-position)
- **Output**
  - `BAM` index file (`.bam.bai`); used in
  [**intersect_extended_primir**](#intersect-extended-primir)


### Quantify workflow

#### `finish_quantify`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter-sam-by-intersecting-primir)
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter-sam-by-intersecting-mirna)
  - Alignments file, sorted, tagged (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort-intersecting-mirna-by-feat-tag)
  - (iso)miR and/or pri-miR counts table (`.tab`); from
  [**merge_tables**](#merge-tables)
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort-uncollapsed-reads-bam-by-position)
  - `BAM` index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index-uncollapsed-reads-bam)

#### `intersect_extendend_primir`

Intersect the aligned reads with the extended pri-miR annotations with
[**BEDTools**](#third-party-software-used).

- **Input**
  - Alignments file (`.bam`); from
  [**sort_all_alns_bam_by_position**](#sort-all-alns-bam-by-position)
  - pri-miR extended annotations (`.gff3`); from
  [**extend_mirs_annotations**](#extend-mirs-annotations)
- **Output**
  - pri-miR intersections file (`.bed`); used in
  [**filter_sam_by_intersecting_primir**](#filter-sam-by-intersecting-primir)
  and [**quantify_primir**](#quantify-primir)


#### `filter_sam_by_intersecting_primir`

Remove alignments that do not intersect with any pri-miR with
[**SAMtools**](#third-party-software-used).

> Required to only intersect alignments within a pri-miR locus.

- **Input**
  - Alignments file (`.sam`); from [**filter_by_indels**](#filter-by-indels) 
  - pri-miR intersections file (`.bed`); from
  [**intersect_extended_primir**](#intersect-extended-primir)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**convert_intersecting_primir_sam_to_bam**](#convert-intersecting-primir-sam-to-bam)
  and [**filter_sam_by_intersecting_mirna**](#filter-sam-by-intersecting-mirna)


#### `convert_intersecting_primir_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with miRNA annotations.

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter-sam-by-intersecting-primir)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_intersecting_primir_bam_by_position**](#sort-intersecting-primir-bam-by-position)


#### `sort_intersecting_primir_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
with miRNA annotations more efficiently.

- **Input**
  - Alignments file (`.bam`); from
  [**convert_intersecting_primir_sam_to_bam**](#convert-intersecting-primir-sam-to-bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_intersecting_primir_bam**](#index-intersecting-primir-bam) and
  [**intersect_extended_mirna**](#intersect-extended-mirna)


#### `index_intersecting_primir_bam`

Create index `BAM` file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
alignments in a genomic region of interest.

- **Input**
  - Alignments file (`.bam`); from
  [**sort_intersecting_primir_bam_by_position**](#sort-intersecting-primir-bam-by-position)
- **Output**
  - `BAM` index file (`.bam.bai`); used in
  [**intersect_extended_mirna**](#intersect-extended-mirna)


#### `intersect_extended_mirna`

Intersect the aligned reads with the extended miRNA annotations with
[**BEDTools**](#third-party-software-used).

- **Input**
  - Alignments file (`.bam`); from
  [**sort_intersecting_primir_bam_by_position**](#sort-intersecting-primir-bam-by-position)
  - miRNA extended annotations (`.gff3`); from
  [**extend_mirs_annotations**](#extend-mirs-annotations)
- **Output**
  - miRNA intersections file (`.bed`); used in
  [**filter_sam_by_intersecting_mirna**](#filter-sam-by-intersecting-mirna)
  and [**add_intersecting_mirna_tag**](#add-intersecting-mirna-tag)


#### `filter_sam_by_intersecting_mirna`

Remove alignments that do not intersect with any miRNA with
[**SAMtools**](#third-party-software-used).

> Required to efficiently classify the alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter-sam-by-intersecting-primir) 
  - miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect-extended-mirna)
- **Output**
  - Alignments file, filtered (`.sam`); used in 
  [**add_intersecting_mirna_tag**](#add-intersecting-mirna-tag) and
  [**uncollapse_reads**](#uncollapse-reads)


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
  [**filter_sam_by_intersecting_mirna**](#filter-sam-by-intersecting-mirna)
  - miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect-extended-mirna)
- **Parameters**
  - **config_template.yaml**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended (default 6)
- **Output**
  - Alignments file, tagged (`.sam`); used in
  [**sort_intersecting_mirna_by_feat_tag**](#sort-intersecting-mirna-by-feat-tag)


#### `sort_intersecting_mirna_by_feat_tag`

Sort the alignments by the tag containing the classified intersecting miRNA
with [**SAMtools**](#third-party-software-used).

> Required for an efficient quantification.

- **Input**
  - Alignments file, tagged (`.sam`); from
  [**add_intersecting_mirna_tag**](#add-intersecting-mirna-tag)
- **Output**
  - Alignments file, tagged, sorted (`.sam`); used in
  [**quantify_mirna**](#quantify-mirna)


#### `quantify_mirna`

Tabulate alignments according to its tag with a
[**custom script**][custom-script-mir-quant].

> Quantification is done with partial counts (_i.e._ each alignment contributes
by the number of collapsed reads divided by the number of hits).

- **Input**
  - Alignments file, sorted, tagged (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort-intersecting-mirna-by-feat-tag)
- **Parameters**
  - **samples.tsv**
    - Library name; specify in sample table column `sample`
  - **config_template.yaml**
    - `mir_list`: miRNA features to be quantified (default isomir, mirna
    pri-miR)
- **Output**
  - (iso)miR counts tab-delimited file; used in
  [**merge_tables**](#merge-tables)


#### `quantify_primir`

Tabulate alignments according to its intersecting pri-miR with a
[**custom script**][custom-script-pri-quant]

> Quantification is done with partial counts (_i.e._ each alignment contributes
by the number of collapsed reads divided by the number of hits).

- **Input**
  - pri-miR intersections file (`.bed`); from
  [**intersecti_extended_primir**](#intersect-extended-primir)
- **Output**
  - pri-miR counts tab-delimited file; used in
  [**merge_tables**](#merge-tables)


#### `merge_tables`

Merge all the tables from the different libraries into a single one with a
[**custom script**][custom-script-merge-tab].

- **Input**
  - counts tab-delimited file; from [**quantify_mirna**](#quantify-mirna)
  and/or [**quantify_primir**](#quantify-primir)
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
  [**filter_sam_by_intersecting_mirna**](#filter-sam-by-intersecting-mirna)
- **Output**
  - Uncollapsed aligned reads (`.sam`); used in
  [**convert_uncollapsed_reads_sam_to_bam**](#convert-uncollapsed-reads-sam-to-bam)


#### `convert_uncollapsed_reads_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter-sam-by-intersecting-mirna)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_uncollapsed_reads_bam_by_position**](#sort-uncollapsed-reads-bam-by-position)


#### `sort_uncollapsed_reads_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.bam`); from
  [**convert_uncollapsed_reads_sam_to_bam**](#convert-uncollapsed_reads-sam-to-bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_uncollapsed_reads_bam**](#index-uncollapsed-reads-bam)


#### `index_uncollapsed_reads_bam`

Create index `BAM` file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
alignments in a genomic region of interest.

- **Input**
  - Alignments file (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort-uncollapsed-reads-bam-by-position)
- **Output**
  - `BAM` index file (`.bam.bai`)

[chr-maps]: <https://github.com/dpryan79/ChromosomeMappings>
[custom-script-blocksort]: <https://github.com/zavolanlab/mirflowz/scripts/blocksort.sh>
[custom-script-filter-mm]: <https://github.com/zavolanlab/mirflowz/scripts/filter_multimappers.py>
[custom-script-get-lines]: <https://github.com/zavolanlab/mirflowz/scripts/get_lines_w_pattern.sh>
[custom-script-gtf-bed]: <https://github.com/zavolanlab/mirflowz/scripts/gtf_exons_bed.1.1.2.R>
[custom-script-iso-tag]: <https://github.com/zavolanlab/mirflowz/scripts/iso_name_tagging.py>
[custom-script-map-chr]: <https://github.com/zavolanlab/mirflowz/scripts/map_chromosomes.pl>
[custom-script-merge-tab]: <https://github.com/zavolanlab/mirflowz/scripts/merge_tables.R>
[custom-script-mir-ext]: <https://github.com/zavolanlab/mirflowz/scripts/mirna_extension.py>
[custom-script-mir-quant]: <https://github.com/zavolanlab/mirflowz/scripts/mirna_quantification.py>
[custom-script-nh-filter]: <https://github.com/zavolanlab/mirflowz/scripts/nh_filter.py>
[custom-script-oligo-sam]: <https://github.com/zavolanlab/mirflowz/scripts/oligomapOutputToSam_nhfiltered.py>
[custom-script-pri-quant]: <https://github.com/zavolanlab/mirflowz/scripts/primir_quantification.py>
[custom-script-remove-dup]: <https://github.com/zavolanlab/mirflowz/scripts/sam_remove_duplicates_inferior_alignments_multimappers.pl>
[custom-script-sam-trx]: <https://github.com/zavolanlab/mirflowz/scripts/sam_trx_to_sam_gen.pl>
[custom-script-trim-id]: <https://github.com/zavolanlab/mirflowz/scripts/trim_id_fasta.sh>
[custom-script-uncollapse]: <https://github.com/zavolanlab/mirflowz/scripts/sam_uncollapse.pl>
[custom-script-validation]: <https://github.com/zavolanlab/mirflowz/scripts/validation_fasta.py>
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
[license-bsl1]: <https://opensource.org/license/bsl-1-0/>
[license-gpl2]: <https://opensource.org/licenses/GPL-2.0>
[license-gpl3]: <https://opensource.org/license/gpl-3-0/>
[license-mit]: <https://opensource.org/licenses/MIT>
[pub-bedtools]: <https://academic.oup.com/bioinformatics/article/26/6/841/244688>
[pub-cufflinks]: <https://doi.org/10.1038/nprot.2012.016>
[pub-cutadapt]: <https://doi.org/10.14806/ej.17.1.200>
[pub-oligomap]: <https://doi.org/10.1016/j.ymeth.2007.10.002 >
[pub-samtools]: <https://doi.org/10.1093/bioinformatics/btp352>
[pub-segemehl]: <https://doi.org/10.1371/journal.pcbi.1000502>
[rule-graph]: images/rule_graph.svg
