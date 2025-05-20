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
| **ASCII-style alignment pileups** | [Apache 2.0][license-apache2] | _"Generates ASCII-style pileups of read alignments in one or more BAM files for one or more genomic regions."_ | [code][code-ascii] | 
| **BEDTools** | [GPLv2][license-gpl2] | _"[...] intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][docs-bedtools] / [publication][pub-bedtools] |
| **cufflinks** | [BSL-1.0][license-bsl1] | _"[...] assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples"_ | [code][code-cufflinks] / [manual][docs-cufflinks] / [publication][pub-cufflinks] |
| **cutadapt** | [MIT][license-mit] | _"[...] finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads"_ | [code][code-cutadapt] / [manual][docs-cutadapt] / [publication][pub-cutadapt] |
| **FASTX-Toolkit** | [AGPL-3.0][license-agpl3] | _"[...] collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing"_ | [code][code-fastx] / [manual][docs-fastx] |
| **GFFUtils** | [AFL-3][license-afl3] | _"[...] a small set of utility programs for working with GFF and GTF files"_ | [code][code-gffutils] / [manual][docs-gffutils] |
| **Oligomap** | [GPLv3][license-gpl3] | _"[...] fast identification of nearly-perfect matches of small RNAs in sequence databases. It allows to exhaustively identify all the perfect and 1-error (where an error is defined to be a mismatch, insertion or deletion) matches of large sets of small RNAs to target sequences"_ | [code][code-oligomap] / [publication][pub-oligomap] |
| **SAMtools** | [MIT][license-mit] | _"[...] suite of programs for interacting with high-throughput sequencing data"_ | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools] |
| **segemehl** | [GPLv3][license-gpl3] | _"[...] map short sequencer reads to reference genomes"_ | [manual][docs-segemehl] / [publication][pub-segemehl] |

## Description of workflow steps

> The workflow consists of five Snakemake files: A main `Snakefile` and an
> individual Snakemake file each for the genome resources preparation, the
> reads mapping, the miRNA quantification and the ASCII-style pileups
> generation. The main `Snakefile` contains the configuration file validation
> and imports the various functional modules described below. Individual steps
> of the workflow are described briefly along with some examples, and links to
> the respective software manuals are given. Parameters that can be modified by
> the user (via the samples table and the configuration file) are also
> described.

### Rule graph

![rule_graph][rule-graph]

Visual representation of the workflow. Automatically prepared with
[Snakemake][docs-snakemake].

### Preparatory

#### Read sample table

##### Requirements

- Tab-separated values (`.tsv`) file
- First row has to contain parameter names as in 
[`samples_table.tsv`](test/test_files/samples_table.tsv) 
- First column used as sample identifiers

Parameter name | Description | Data type(s)
 --- | --- | --- 
sample | Arbitrary name for the miRNA sequence library. | `str`
sample_file | Path to the `gzip`ped miRNA sequencing library file. The path must be relative to the directory where the workflow will be run. | `str`
adapter | Sequence of the 3'-end adapter used during library preparation. Required for [Cutadapt](#third-party-software-used). Use a value such as `XXXXXXXXXX` if no adapter is present or if no trimming is desired. | `str`
format | One of `fa`/`fasta` or `fq`/`fastq`, if the library file is in FASTA or FASTQ format, respectively. | `str`


#### Configuration file

Some parameters within the workflow can be modified. Refer to the
[configuration template](#config/config_template.yaml) for a detailed
explanation of each option.


### Snakefile

#### `finish`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - (**Workflow output**) SAM file with the pri-miR intersecting alignments;
  from [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir)
  - (**Workflow output**) SAM file with the mature miRNA intersecting alignments; from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  - (**Workflow output**) (iso)miR and/or pri-miR counts table (`.tab`); from
  [**merge_tables**](#merge_tables)
  - (**Workflow output**) BAM file with the contributing alignments, sorted;
  from [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - (**Workflow output**) BAM index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - (**Workflow output**) Empty text file (`.txt`)
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
  - SAM header (`.sam`); from
  [**create_genome_header**](#create_genome_header)
  - Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`);
  from [**extract_chr_len**](#extract_chr_len)
  - Primary miRNA transcript (pri-miR) extended annotations (`.gff3`);
  from [**extend_mirs_annotations**](#extend_mirs_annotations)
  - Mature miRNA (miRNA) extended annotations (`.gff3`);
  from [**extend_mirs_annotations**](#extend_mirs_annotations)


#### `trim_genome_seq_ids`

Trim genome sequence IDs with a [**custom script**][custom-script-trim-id].

- **Input**
  - (**Workflow input**) Genome sequence, `gzip`ed (`.fa.gz`/`.fasta.gz`)
- **Output**
  - Genome sequence, trimmed IDs (`.fa`); used in
  [**extract_transcriptome_seqs**](#extract_transcriptome_seqs),
  [**create_genome_header**](#create_genome_header),
  [**create_index_genome_fasta**](#create_index_genome),
  [**generate_segemehl_index_genome**](#generate_segemehl_index_genome),
  [**mapping_genome_segemehl**](#mapping_genome_segemehl),
  [**mapping_genome_oligomap**](#mapping_genome_oligomap) and
  [**compress_reference_genome**](#compress_reference_genome)


#### `extract_transcriptome_seqs`

Create transcriptome from genomic sequence and annotations with 
[**cufflinks**](#third-party-software-used).

- **Input**
  - (**Workflow input**) Genome annotations, `gzip`ed (`.gtf.gz`)
  - Genome sequence, trimmed IDs (`.fa`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - Transcriptome sequence (`.fa`); used in 
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)


#### `trim_transcriptome_seq_ids`

Trim transcriptome sequence IDs with a
[**custom script**][custom-script-trim-id].

- **Input**
  - Transcriptome sequence (`.fa`); from
  [**extract_transcriptome_seqs**](#extract_transcriptome_seqs)
- **Output**
  - Transcriptome sequence, trimmed IDs (`.fa`); used in
  [**generate_segemehl_index_transcriptome**](#generate_segemehl_index_transcriptome),
  [**mapping_transcriptome_segemehl**](#mapping_transcriptome_segemehl) and
  [**mapping_transcriptome_oligomap**](#mapping_transcriptome_oligomap)


#### `generate_segemehl_index_transcriptome`

Generate transcriptome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The transcriptome index only needs to be generated once for each combination
> of transcriptome sequence and annotations, and sample sets.

- **Input**
  - Transcriptome sequence, trimmed IDs (`.fa`); from
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)
- **Output**
  - segemehl transcriptome index (`.idx`); used in 
  [**mapping_transcriptome_segemehl**](#mapping_transcriptome_segemehl)


#### `generate_segemehl_index_genome`

Generate genome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The genome index only needs to be generated once for each combination
> of annotations and sample sets.

- **Input**
  - Genome sequence, trimmed IDs (`.fa`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - segemehl genome index (`.idx`); used in 
  [**mapping_genome_segemehl**](#mapping_genome_segemehl)


#### `get_exons_gtf`

Retrieve exon annotations from genome annotations with a
[**custom script**][custom-script-get-lines].

- **Input**
  - (**Workflow input**) Genome annotations, `gzip`ed (`.gtf.gz`)
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

Create SAM header for the genome with 
[**SAMtools**](#third-party-software-used).

> Required by [SAMtools](#third-party-software-used) to work with the
> alignment files.

- **Input**
  - Genome sequence, trimmed IDs (`.fa`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - SAM genome header (`.sam`); used in
  [add_header_all_maps](#add_header_all_maps)


#### `map_chr_names`

Map UCSC-like chromosome names with Ensembl-like ones in miRNA annotations
with a [**custom script**][custom-script-map-chr].

> Required by [BEDTools](#third-party-software) to intersect alignments with
> miRNA annotations. Several mapping tables are available [here][chr-maps].

- **Input**
  - (**Workflow input**) miRNA annotations (`.gff3`)
  - (**Workflow input**) Tab-separated chromosome name mappings table (`.tsv`)
- **Output**
  - miRNA annotations, mapped chromosome name(s) (`.gff3`); used in 
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `create_index_genome_fasta`

Create a FASTA index for the genome with 
[**SAMtools**](#third-party-software-used).

- **Input**
  - Genome sequence, trimmed IDs (`.fa`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - FASTA genome index (`.fa.fai`); used in
  [**extract_chr_len**](#extract_chr_len)


#### `extract_chr_len`

Extract chromosome(s) length from the genome sequence.

> Required to ensure that the extended annotations in generated in the
> [**extend_mirs_annotations**](#extend_mirs_annotations) rule do not exceed
> the chromosome(s) boundaries.

- **Input**
  - FASTA genome index (`.fa.fai`); from
  [**create_index_genome_fasta**](#create_index_genome_fasta)
- **Output**
  - Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`); used
  in [**extend_mirs_annotations**](#extend_mirs_annotations)


#### `extend_mirs_annotations`

Extend miRNA annotations, ensure feature names uniqueness and split the file
by feature with a [**custom script**][custom-script-mir-ext].

> Adjust miRNAs' 'Name' attribute to account for the different genomic
> locations the miRNA sequence is annotated on and ensure their uniqueness.
> The name format is `SPECIES-mir-NAME-#` for pri-miRs, and 
> `SPECIES-miR-NAME-#-ARM` or `SPECIES-miR-NAME-#` for mature miRNA with both
> or just one arm respectively, where `#` is the replica integer. If a pri-miR
> has a replica but its number is set in the 'ID' attribute, the first instance
> does not has a suffix but the other one(s) do. If a precursor has no other
> occurrences, no further modifications are made. On the other hand,
> mature miRNA regions are extended on both sides to account for isomiR species
> with shifted start and/or end positions without exceeding chromosome(s)
> boundaries. If required, pri-miR loci are also extended to accommodate the
> new miRNA coordinates. In addition, pri-miR names are modified to record the
> final positions by appending `_-y` and `_+x` to them, where `y` is the 5'
> shift and `x` the 3' shift.

- **Input**
  - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
  [**map_chr_names**](#map_chr_names)
- **Parameters**
  - **config_template.yaml**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended at most (default: 6)
- **Output**
  - Primary miRNA transcript (pri-miR) extended annotations (`.gff3`); used in
  [**intersect_extended_primir**](#intersect_extended_primir)
  - Mature miRNA (miRNA) extended annotations (`.gff3`); used in
  [**intersect_extended_mirna**](#intersect_extended_mirna)
- **Examples**

```console
Example 1 | Extension | Mature miRNA extension

IN:
    pri-miR entry:
        19	.	miRNA_primary_transcript	2517	2614	.	+	.	ID=MI0003141;Alias=MI0003141;Name=hsa-mir-512-2
    mature miRNA entry:
        19	.	miRNA	2536	2558	.	+	.	ID=MIMAT0002822_1;Alias=MIMAT0002822;Name=hsa-miR-512-5p;Derives_from=MI0003141
    extension:
        6
OUT:
    pri-miR entry:
        19	.	miRNA_primary_transcript	2517	2614	.	+	.	ID=MI0003141;Alias=MI0003141;Name=hsa-mir-512-2_-0_+0
    mature miRNA entry:
        19	.	miRNA	2530	2564	.	+	.	ID=MIMAT0002822_1;Alias=MIMAT0002822;Name=hsa-miR-512-2-5p;Derives_from=MI0003141


Example 2 | Extension | Mature miRNA and pri-miR extension

IN:
    pri-miR entry:
        19	.	miRNA_primary_transcript	9	122	.	+	.	ID=MI0003140;Alias=MI0003140;Name=hsa-mir-512-1
    mature miRNA entry:
        19	.	miRNA	12	74	.	+	.	ID=MIMAT0002822;Alias=MIMAT0002822;Name=hsa-miR-512-5p;Derives_from=MI0003140
    extension:
        6
OUT:
    pri-miR entry:
        19	.	miRNA_primary_transcript	6	122	.	+	.	ID=MI0003140;Alias=MI0003140;Name=hsa-mir-512-1_-3_+0
    mature miRNA entry:
        19	.	miRNA	6	80	.	+	.	ID=MIMAT0002822;Alias=MIMAT0002822;Name=hsa-miR-512-1-5p;Derives_from=MI0003140


Example 3 | Extension | Matrue miRNA exceeding chromosome boundaries extension

IN:
    pri-miR entry:
        19	.	miRNA_primary_transcript	2	122	.	+	.	ID=MI0003140;Alias=MI0003140;Name=hsa-mir-512-1
    mature miRNA entry:
        19	.	miRNA	3	74	.	+	.	ID=MIMAT0002822;Alias=MIMAT0002822;Name=hsa-miR-512-5p;Derives_from=MI0003140
    extension:
        6
OUT:
    pri-miR entry:
        19	.	miRNA_primary_transcript	1	122	.	+	.	ID=MI0003140;Alias=MI0003140;Name=hsa-mir-512-1_-1_+0
    mature miRNA entry:
        19	.	miRNA	1	80	.	+	.	ID=MIMAT0002822;Alias=MIMAT0002822;Name=hsa-miR-512-1-5p;Derives_from=MI0003140


Example 4 | Name uniqueness | Replica number in the ID

IN:
    pri-miR entries:
        chr21	.	miRNA_primary_transcript	8206563	8206618	.	+	.	ID=MI0033425;Alias=MI0033425;Name=hsa-mir-10401
        chr21	.	miRNA_primary_transcript	8250772	8250827	.	+	.	ID=MI0033425_2;Alias=MI0033425;Name=hsa-mir-10401
    mature miRNA entries:
        chr21	.	miRNA	8206563	8206582	.	+	.	ID=MIMAT0041633;Alias=MIMAT0041633;Name=hsa-miR-10401-5p;Derives_from=MI0033425
        chr21	.	miRNA	8206598	8206618	.	+	.	ID=MIMAT0041634;Alias=MIMAT0041634;Name=hsa-miR-10401-3p;Derives_from=MI0033425
        chr21	.	miRNA	8250772	8250791	.	+	.	ID=MIMAT0041633_1;Alias=MIMAT0041633;Name=hsa-miR-10401-5p;Derives_from=MI0033425
        chr21	.	miRNA	8250807	8250827	.	+	.	ID=MIMAT0041634_1;Alias=MIMAT0041634;Name=hsa-miR-10401-3p;Derives_from=MI0033425
OUT:
    pri-miR entries:
        chr21	.	miRNA_primary_transcript	8206563	8206618	.	+	.	ID=MI0033425;Alias=MI0033425;Name=hsa-mir-10401
        chr21	.	miRNA_primary_transcript	8250772	8250827	.	+	.	ID=MI0033425_2;Alias=MI0033425;Name=hsa-mir-10401-2
    mature miRNA entries:
        chr21	.	miRNA	8206563	8206582	.	+	.	ID=MIMAT0041633;Alias=MIMAT0041633;Name=hsa-miR-10401-5p;Derives_from=MI0033425
        chr21	.	miRNA	8206598	8206618	.	+	.	ID=MIMAT0041634;Alias=MIMAT0041634;Name=hsa-miR-10401-3p;Derives_from=MI0033425
        chr21	.	miRNA	8250772	8250791	.	+	.	ID=MIMAT0041633_1;Alias=MIMAT0041633;Name=hsa-miR-10401-2-5p;Derives_from=MI0033425
        chr21	.	miRNA	8250807	8250827	.	+	.	ID=MIMAT0041634_1;Alias=MIMAT0041634;Name=hsa-miR-10401-2-3p;Derives_from=MI0033425


Example 5 | Name uniqueness | Replica number in the Name; single mature arm

IN:
    pri-miR entries:
        chr21	.	miRNA_primary_transcript	8205315	8205406	.	+	.	ID=MI0022559;Alias=MI0022559;Name=hsa-mir-6724-1
        chr21	.	miRNA_primary_transcript	8249505	8249596	.	+	.	ID=MI0031516;Alias=MI0031516;Name=hsa-mir-6724-2
    mature miRNA entries:
        chr21	.	miRNA	8205325	8205347	.	+	.	ID=MIMAT0025856;Alias=MIMAT0025856;Name=hsa-miR-6724-5p;Derives_from=MI0022559
        chr21	.	miRNA	8249515	8249537	.	+	.	ID=MIMAT0025856_1;Alias=MIMAT0025856;Name=hsa-miR-6724-5p;Derives_from=MI0031516
OUT:
    pri-miR entries:
        chr21	.	miRNA_primary_transcript	8205315	8205406	.	+	.	ID=MI0022559;Alias=MI0022559;Name=hsa-mir-6724-1
        chr21	.	miRNA_primary_transcript	8249505	8249596	.	+	.	ID=MI0031516;Alias=MI0031516;Name=hsa-mir-6724-2
    mature miRNA entries:
        chr21	.	miRNA	8205325	8205347	.	+	.	ID=MIMAT0025856;Alias=MIMAT0025856;Name=hsa-miR-6724-1-5p;Derives_from=MI0022559
        chr21	.	miRNA	8249515	8249537	.	+	.	ID=MIMAT0025856_1;Alias=MIMAT0025856;Name=hsa-miR-6724-2-5p;Derives_from=MI0031516


Example 6 | Name uniqueness | Both mature miRNA arms but just one with the replica number

IN:
    pri-miR entries:
        chr2	.	miRNA_primary_transcript	135665397	135665478	.	+	.	ID=MI0000447;Alias=MI0000447;Name=hsa-mir-128-1
        chr3	.	miRNA_primary_transcript	35744476	35744559	.	+	.	ID=MI0000727;Alias=MI0000727;Name=hsa-mir-128-2
    mature miRNA entries:
        chr2	.	miRNA	135665446	135665466	.	+	.	ID=MIMAT0000424;Alias=MIMAT0000424;Name=hsa-miR-128-3p;Derives_from=MI0000447
        chr2	.	miRNA	135665411	135665433	.	+	.	ID=MIMAT0026477;Alias=MIMAT0026477;Name=hsa-miR-128-1-5p;Derives_from=MI0000447
        chr3	.	miRNA	35744527	35744547	.	+	.	ID=MIMAT0000424_1;Alias=MIMAT0000424;Name=hsa-miR-128-3p;Derives_from=MI0000727
        chr3	.	miRNA	35744490	35744512	.	+	.	ID=MIMAT0031095;Alias=MIMAT0031095;Name=hsa-miR-128-2-5p;Derives_from=MI0000727
OUT:
    pri-miR entries:
        chr2	.	miRNA_primary_transcript	135665397	135665478	.	+	.	ID=MI0000447;Alias=MI0000447;Name=hsa-mir-128-1
        chr3	.	miRNA_primary_transcript	35744476	35744559	.	+	.	ID=MI0000727;Alias=MI0000727;Name=hsa-mir-128-2
    mature miRNA entries:
        chr2	.	miRNA	135665446	135665466	.	+	.	ID=MIMAT0000424;Alias=MIMAT0000424;Name=hsa-miR-128-1-3p;Derives_from=MI0000447
        chr2	.	miRNA	135665411	135665433	.	+	.	ID=MIMAT0026477;Alias=MIMAT0026477;Name=hsa-miR-128-1-5p;Derives_from=MI0000447
        chr3	.	miRNA	35744527	35744547	.	+	.	ID=MIMAT0000424_1;Alias=MIMAT0000424;Name=hsa-miR-128-2-3p;Derives_from=MI0000727
        chr3	.	miRNA	35744490	35744512	.	+	.	ID=MIMAT0031095;Alias=MIMAT0031095;Name=hsa-miR-128-2-5p;Derives_from=MI0000727
```

### Map workflow

#### `finish_map`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - BAM index file (`.bam.bai`); from
  [**index_all_alns_bam**](#index_all_alns_bam)


#### `start`

Copy and rename read files.

> Local rule.
> Depending on the library file format, the output file undergoes a quality
> filter (`fa`/`.fastq`) or is directly formatted (`.fa`/`.fasta`).

- **Input**
  - (**Workflow input**) miRNA sequencing library, `gzip`ed
  (`.fa.gz`/`.fasta.gz` or `.fq.gz`/`.fastq.gz`)
- **Output**
  - miRNA sequencing library, copied, renamed (`.fa`, `.fastq`); used in
  [**fastq_quality_filter**](#fastq_quality_filter) and/or 
  [**format_fasta**](#format_fasta)


#### `fastq_quality_filter`

Conduct quality control for reads library with
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - miRNA sequencing library, copied, renamed (`.fastq`); from
  [**start**](#start)
- **Parameters**
  - **config_template.yaml**
    - `q_value`: Minimum Q (Phred) score to keep (default: 10)
    - `p_value`: Minimum % of bases that must have a Q (Phred) quality
    (default: 50)
- **Output**
  - miRNA sequencing library, filtered (`.fastq`); used in
  [**fastq_to_fasta**](#fastq_to_fasta)


#### `fastq_to_fasta`

Convert reads file from FASTQ to FASTA with 
[**fastx_toolkit**](#third-party-software-used).

> Sequence identifiers are renamed to numbers.

- **Input**
  - miRNA sequencing library, filtered (`.fastq`); from
  [**fastq_quality_filter**](#fastq-quality-filter)
- **Output**
  - miRNA sequencing library (`.fa`); used in [**format_fasta**](#format_fasta)


#### `format_fasta`

Format read's sequences to appear on a single line with
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - miRNA sequencing library (`.fa`); from [**start**](#start) or 
  [**fastq_to_fasta**](#fastq_to_fasta)
- **Output**
  - miRNA sequencing library, formatted (`.fasta`); used in
  [**remove_adapters**](#remove_adapters)


#### `remove_adapters`

Trim 3' adapters and `N` bases at either end. Filter reads by minimum length
and number of inner `N` bases with [**cutadapt**](#third-party-software-used).

- **Input**
  - miRNA sequencing library, formatted (`.fasta`); from
  [**format_fasta**](#format_fasta)
- **Parameters**
  - **samples.tsv**
    - Adapter to be removed; specified in the sample's table column `adapter`
  - **config_template.yaml**
    - `error_rate`: Fraction of allowed errors in the matched adapters
    (default: 0.1)
    - `overlap`: Minimum overlap length between adapter and read to trim the
    bases (default: 3)
    - `minimum_length`: Minimum length for a processed read to be kept
    (default: 15)
    - `max_n`: Maximum number of inner `N` bases for a processed read to be
    kept (default: 0)
- **Output**
  - miRNA sequencing library, filtered, without adapters (`.fasta`); used in 
  [**collapse_identical_reads**](#collapse_identical_reads)


#### `collapse_identical_reads`

Collapse and rename identical reads
[**fastx_toolkit**](#third-party-software-used).

> Sequences are renamed in the format `R-N`, where `R` is the assigned number
> to the unique entry, and `N` is the amount of identical sequences within the 
> library collapsed in it.

- **Input**
  - miRNA sequencing library, filtered, without adapters (`.fasta`); from 
  [**remove_adapters**](#remove_adapters)
- **Output**
  - miRNA sequencing library, collapsed, renamed (`.fasta`); used in
  [**filter_fasta_for_oligomap**](#filter_fasta_for_oligomap),
  [**map_genome_segemehl**](#map_genome_segemehl) and
  [**map_transcriptome_segemehl**](#map_transcriptome_segemehl)


#### `map_genome_segemehl`

Align short reads to reference genome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - miRNA sequencing library, collapsed, renamed (`.fasta`); from
  [**collapse_identical_reads**](#collapse_identical_reads)
  - Genome sequence, trimmed IDs (`.fa`); from 
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
  - segemehl genome index (`.idx`); from
  [**generate_segemehl_index_genome**](#generate_segemehl_index_genome)
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_genome_maps**](#merge_genome_maps)


#### `map_transcriptome_segemehl`

Align short reads to reference transcriptome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - miRNA sequencing library, collapsed, renamed (`.fasta`); from
  [**collapse_identical_reads**](#collapse_identical_reads)
  - Transcriptome sequence, trimmed IDs (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)
  - segemehl transcriptome index (`.idx`); from
  [**generate_segemehl_index_transcriptome**](#generate_segemehl_index_transcriptome)
- **Output**
  - Alignments file (`.sam`); used in
  [**merge_transcriptome_maps**](#merge_transcriptome_maps)


#### `filter_fasta_for_oligomap`

Filter reads by length with a [**custom script**][custom-script-validation].

> Required for an optimal mapping speed. Oligomap is specifically written for
> short reads. Therefore, reads with more bases than the default maximum (30
> nts) makes the mapping slower.

- **Input**
  - miRNA sequencing library, collapsed, renamed (`.fasta`); from
  [**collapse_identical_reads**](#collapse_identical_reads)
- **Parameters**
  - **config_template.yaml**
    - `max_length_reads`: Maximum length of processed reads to be mapped with
    [**oligomap**](#third-party-software-used) (default: 30)
- **Output**
  - miRNA sequencing library, collapsed, filtered (`.fasta`); used in 
  [**map_genome_oligomap**](#map_genome_oligomap) and
  [**map_transcriptome_oligomap**](#map_transcriptome_oligomap)


#### `map_genome_oligomap`

Align short reads to reference genome with
[**oligomap**](#third-party-software-used).

> Refer to Oligomap's [**Output format section**][oligomap-out] for a specific
> explanation and examples on the output format.

- **Input**
  - miRNA sequencing library, collapsed, filtered (`.fasta`); from
  [**filter_fasta_for_oligomap**](#filter_fasta_for_oligomap)
  - Genome sequence, trimmed IDs (`.fa`); from 
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - Alignments file (`.oligomap`); used in
  [**sort_genome_oligomap**](#sort_genome_oligomap)
  - Alignments report (`.txt`)


#### `sort_genome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

- **Input**
  - Alignments file (`.oligomap`); from
  [**map_genome_oligomap**](#map_genome_oligomap)
- **Output**
  - Alignments file, sorted (`.oligomap`); used in
  [**convert_genome_to_sam_oligomap**](#convert_genome_to_sam_oligomap)


#### `convert_genome_to_sam_oligomap`

Convert aligned reads `.oligomap` to `.sam` and filter alignments by number of
hits with a [**custom script**][custom-script-oligo-sam].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

- **Input**
  - Alignments file, sorted (`.oligomap`); from
  [**sort_genome_oligomap**](#sort_genome_oligomap)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default: 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**merge_genome_maps**](#merge_genome_maps)


#### `map_transcriptome_oligomap`

Align short reads to reference transcriptome with
[**oligomap**](#third-party-software-used).

> Refer to Oligomap's [Output format section][oligomap-out] for a specific
> explanation and examples on the output format.

- **Input**
  - miRNA sequencing library, collapsed, filtered (`.fasta`); from
  [**filter_fasta_for_oligomap**](#filter_fasta_for_oligomap)
  - Transcriptome sequence, trimmed IDs (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim_transcriptome_seq_ids)
- **Output**
  - Alignments file (`.oligomap`); used in
  [**sort_transcriptome_oligomap**](#sort_transcriptome_oligomap)
  - Alignments report (`.txt`)


#### `sort_transcriptome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

- **Input**
  - Alignments file (`.oligomap`); from
  [**map_transcriptome_oligomap**](#map_transcriptome_oligomap)
- **Output**
  - Alignments file, sorted (`.oligomap`); used in
  [**convert_transcriptome_to_sam_oligomap**](#convert_transcriptome_to_sam_oligomap)


#### `convert_transcriptome_to_sam_oligomap`

Convert aligned reads `.oligomap` to `.sam` and filter alignments by number of
hits with a [**custom script**][custom-script-oligo-sam].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

- **Input**
  - Alignments file, sorted (`.oligomap`); from
  [**sort_transcriptome_oligomap**](#sort_transcriptome_oligomap)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default: 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**merge_transcriptome_maps**](#merge_transcriptome_maps)


#### `merge_genome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) genome alignments.

- **Input**
  - Alignments file (`.sam`); from
  [**map_genome_segemehl**](#map_genome_segemehl)
  - Alignments file, filtered (`.sam`); from
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
  - Alignments file, filtered (`.sam`); from
  [**convert_transcriptome_to_sam_oligomap**](#convert_transcriptome_to_sam_oligomap)
- **Output**
  - Alignments file (`.sam`); used in
  [**filter_transcriptome_by_nh**](#filter_transcriptome_by_nh)


#### `filter_genome_by_nh`

Filter merged genome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

- **Input**
  - Alignments file (`.sam`); from
  [**merge_genome_maps**](#merge_genome_maps)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default: 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**remove_header_genome_mappings**](#remove_header_genome_mappings)


#### `filter_transcriptome_by_nh`

Filter merged transcriptome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

- **Input**
  - Alignments file (`.sam`); from
  [**merge_transcriptome_maps**](#merge_transcriptme_maps)
- **Parameters**
  - **config_template.yaml**
    - `nh`: Maximum number of mappings per read to be kept (default: 100)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**remove_header_transcriptome_mappings**](#remove_header_transcriptome_mappings)


#### `remove_header_genome_mappings`

Remove the SAM header of the genome alignments file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_genome_by_nh**](#filter_genome_by_nh)
- **Output**
  - Alignments file, without SAM header (`.sam`); used in
  [**merge_all_maps**](#merge_all_maps)


#### `remove_header_transcriptome_mappings`

Remove the SAM header of the transcriptome alignments file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file (`.sam`); from
  [**filter_transcriptome_by_nh**](#filter_transcriptome_by_nh)
- **Output**
  - Alignments file, without SAM header (`.sam`); used in 
  [**transcriptome_to_genome_maps**](#transcriptome_to_genome_maps)


#### `transcriptome_to_genome_maps`

Convert the alignments' transcriptome coordinates to genomic ones with a
[**custom script**][custom-script-sam-trx].

- **Input**
  - Alignments file, without SAM header (`.sam`); from 
  [**remove_header_transcriptome_mappings**](#remove_header_transcriptome_mappings)
  - Exon annotations (`.bed`); from
  [**convert_exons_gtf_to_bed**](#convert_exons_gtf_to_bed)
- **Output**
  - Alignments file, without SAM header (`.sam`); used in
  [**merge_all_maps**](#merge_all_maps)


#### `merge_all_maps`

Concatenate the four alignment files into a single one.

- **Input**
  - Alignments files, without SAM header (`.sam`); from
  [**remove_header_genome_mappings**](#remove_header_genome_mappings) and
  [**transcriptome_to_genome_maps**](#transcriptome_to_genome_maps)
- **Output**
  - Alignments file, without SAM header (`.sam`); used in
  [**add_header_all_maps**](#add_header_all_maps)


#### `add_header_all_maps`

Add the SAM header to the aligned reads merged file with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file, without SAM header (`.sam`); from
  [**merge_all_maps**](#merge_all_maps)
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
> the fields `QNAME`, `FLAG`, `RNAME`, `POS` and `CIGAR`. 
> Alignments are considered to be inferiors if having the same `QNAME` and
> a bigger edit distance than the smaller one within the group. The tags `NH`
> (number of hits) and `HI` (query hit index) are updated accordingly.

- **Input**
  - Alignments file, sorted (`.sam`); from
  [**sort_maps_by_id**](#sort_maps_by_id)
- **Output**
  - Alignments file, filtered (`.sam`); used in
  [**filter_by_indels**](#filter_by_indels)
- **Examples**

```console
Example 1 | Remove duplicates

IN:
    1-2	0	19	44414	1	21M	*	0	0	GAAGGCGCTTCCCTTTGGAGT	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:21	RG:Z:A1	YZ:Z:0
    1-2	0	19	44414	255	21M	*	0	0	GAAGGCGCTTCCCTTTGGAGT	*	NM:i:0	MD:Z:21	NH:i:1
OUT:
    1-2	0	19	44414	255	21M	*	0	0	GAAGGCGCTTCCCTTTGGAGT	*	MD:Z:21	NH:i:1	NM:i:0


Example 2 | Remove inferiors single alignment

IN:
    1-704	16	19	207362	1	18M	*	0	0	CCCGGGCCCGGCGCGCCG	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:18	RG:Z:A1	YZ:Z:0
    1-704	272	19	471264	1	16M1I1M	*	0	0	CCCGGGCCCGGCGCGCCG	*	HI:i:1	NH:i:2	NM:i:2	MD:Z:11G5	RG:Z:A1	YZ:Z:0
OUT:
    1-704	16	19	207362	1	18M	*	0	0	CCCGGGCCCGGCGCGCCG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:18	RG:Z:A1	YZ:Z:0


Example 3 | Remove inferiors multiple alignments

IN:

    1-1197	0	19	56327	1	15M	*	0	0	TATGGCACTGGTAGA	*	HI:i:0	NH:i:4	NM:i:2	MD:Z:1C11T1	RG:Z:A1	YZ:Z:0
    1-1197	256	19	68983	1	15M	*	0	0	TATGGCACTGGTAGA	*	HI:i:1	NH:i:4	NM:i:3	MD:Z:1C10AT1	RG:Z:A1	YZ:Z:0
    1-1197	256	19	76967	1	15M	*	0	0	TATGGCACTGGTAGA	*	HI:i:2	NH:i:4	NM:i:2	MD:Z:1C11T1	RG:Z:A1	YZ:Z:0
    1-1197	256	19	92363	1	15M	*	0	0	TATGGCACTGGTAGA	*	HI:i:4	NH:i:4	NM:i:3	MD:Z:1C11TA	RG:Z:A1	YZ:Z:0

OUT:
    1-1197	0	19	56327	1	15M	*	0	0	TATGGCACTGGTAGA	*	HI:i:0	NH:i:2	NM:i:2	MD:Z:1C11T1	RG:Z:A1	YZ:Z:0
    1-1197	256	19	76967	1	15M	*	0	0	TATGGCACTGGTAGA	*	HI:i:1	NH:i:2	NM:i:2	MD:Z:1C11T1	RG:Z:A1	YZ:Z:0
```


#### `filter_by_indels`

Filter multimappers favoring InDels over mismatches with a
[**custom script**][custom-script-filter-mm].

> Given that InDels are more frequent in miRNAs than mismatches, as
> demonstrated by [Saunders et al. (2017)][cite_saunders],
> [Neilsen et al. (2012)][cite_neilsen] and
> [Schumauch et al. (2024)][cite_schumauch], only those "multimappers" (defined
> here as alignments of the same read mapping to different genomic loci with
> the same edit distance) that contain a higher or equal number of InDels
> compared to mismatches are retained.

- **Input**
  - Alignments file, sorted, filtered (`.sam`); from
  [**remove_inferiors**](#remove_inferiors)
- **Output**
  - Alignments file, sorted, filtered (`.sam`); used in
  [**convert_all_alns_sam_to_bam**](#convert_all_alns_sam_to_bam) and
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primr)
- **Examples**

```console
Example 1 | Different number of InDels

IN:
    1-1	0	19	77595	255	8M1D14M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:3G1T2^A14	NH:i:2	NM:i:3	XA:Z:Q	XI:i:1
    1-1	0	19	330456	255	4M1D1M1I3M1D13M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:4^G4^A13	NH:i:2	NM:i:3	XA:Z:Q	XI:i:0

Alignments:
    CTGACATC-AGTGATTCTCCTGC
    ||| | || |||||||||||||| (1 InDel, 2 mismatches, discarded)
    CTGGCTTCAAGTGATTCTCCTGC

    CTGA-CATCA-GTGATTCTCCTGC
    |||| | ||| ||||||||||||| (3 InDels, 0 mismatches, retained)
    CTGAGC-TCAAGTGATTCTCCTGC

OUT:
    1-1	0	19	330456	255	4M1D1M1I3M1D13M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:4^G4^A13	NH:i:1	HI:i:1  NM:i:3	XA:Z:Q	XI:i:0


Example 2 | Equal number of InDels

IN:
    1-1	0	19	142777	255	5M1I15M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	NM:i:3	XA:Z:Q	XI:i:0
    1-1	0	19	270081	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14G0G4	NH:i:3	NM:i:3	XA:Z:Q	XI:i:2
    1-1	0	19	545543	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	NM:i:3	XA:Z:Q	XI:i:1

Alignments:
    GCTTCAAGCCTCCCACCTAGC
    ||||| |||||||||  |||| (1 Indel, 2 mismatches, retained)
    GCTTC-AGCCTCCCAAGTAGC

    GCTTCAAGCCTCCCACCTAGC
    |||||| ||||||||  |||| (1 Indel, 2 mismatches, retained)
    GCTTCA-GCCTCCCAGGTAGC

    GCTTCAAGCCTCCCACCTAGC
    |||||| ||||||||  |||| (1 Indel, 2 mismatches, retained)
    GCTTCA-GCCTCCCAAGTAGC

OUT:
    1-1	0	19	142777	255	5M1D15M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	HI:i:1  NM:i:3	XA:Z:Q	XI:i:0
    1-1	0	19	270081	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14G0G4	NH:i:3	HI:i:2  NM:i:3	XA:Z:Q	XI:i:2
    1-1	0	19	545543	255	6M1I14M	*	0	0	GCTTCAAGCCTCCCACCTAGC	*	MD:Z:14A0G4	NH:i:3	HI:i:3  NM:i:3	XA:Z:Q	XI:i:1
```


#### `convert_all_alns_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
> with pri-miR annotations.

- **Input**
  - Alignments file, filtered (`.sam`); from [**filter_by_indels**](#filter_by_indels)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_all_alns_bam_by_position**](#sort_all_alns_bam_by_position)


#### `sort_all_alns_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
> with pri-miR annotations more efficiently.

- **Input**
  - Alignments file (`.bam`); from
  [**convert_all_alns_sam_to_bam**](#convert_all_alns_sam_to_bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_all_alns_bam**](#index_all_alns_bam) and
  [**intersect_extended_primir**](#intersect_extended_primir)


#### `index_all_alns_bam`

Create index BAM file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
> alignments in a genomic region of interest.

- **Input**
  - Alignments file, sorted (`.bam`); from
  [**sort_all_alns_bam_by_position**](#sort_all_alns_bam_by_position)
- **Output**
  - BAM index file (`.bam.bai`); used in
  [**intersect_extended_primir**](#intersect_extended_primir)


### Quantify workflow

#### `finish_quantify`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir)
  and [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  - (iso)miR and/or pri-miR counts table (`.tab`); from
  [**merge_tables**](#merge_tables)
  - Alignments file, uncollapsed, sorted (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - BAM index file (`.bam.bai`); from
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)


#### `intersect_extendend_primir`

Intersect the aligned reads with the extended pri-miR annotations with
[**BEDTools**](#third-party-software-used).

> Only those alignments fully intersecting a (possibly extended) pri-miR
> annotated region are kept.

- **Input**
  - Alignments file, sorted (`.bam`); from
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

> Required to only intersect alignments within a (possibly extended) pri-miR
> locus.

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_by_indels**](#filter_by_indels) 
  - pri-miR intersections file (`.bed`); from
  [**intersect_extended_primir**](#intersect_extended_primir)
- **Output**
  - (**Workflow output**) Alignments file, filtered (`.sam`); used in
  [**convert_intersecting_primir_sam_to_bam**](#convert_intersecting_primir_sam_to_bam)
  and [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)


#### `convert_intersecting_primir_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
> with miRNA annotations.

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir)
- **Output**
  - Alignments file (`.bam`); used in
  [**sort_intersecting_primir_bam_by_position**](#sort_intersecting_primir_bam_by_position)


#### `sort_intersecting_primir_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

> Required by [BEDTools](#third-party-software-used) to intersect alignments
> with miRNA annotations more efficiently.

- **Input**
  - Alignments file (`.bam`); from
  [**convert_intersecting_primir_sam_to_bam**](#convert_intersecting_primir_sam_to_bam)
- **Output**
  - Alignments file, sorted (`.bam`); used in
  [**index_intersecting_primir_bam**](#index_intersecting_primir_bam) and
  [**intersect_extended_mirna**](#intersect_extended_mirna)


#### `index_intersecting_primir_bam`

Create index BAM file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
> alignments in a genomic region of interest.

- **Input**
  - Alignments file, sorted (`.bam`); from
  [**sort_intersecting_primir_bam_by_position**](#sort_intersecting_primir_bam_by_position)
- **Output**
  - BAM index file (`.bam.bai`); used in
  [**intersect_extended_mirna**](#intersect_extended_mirna)


#### `intersect_extended_mirna`

Intersect the aligned reads with the extended miRNA annotations with
[**BEDTools**](#third-party-software-used).

> Only those alignments fully intersecting an extended mature miRNA annotated
> region are kept.

- **Input**
  - Alignments file, sorted (`.bam`); from
  [**sort_intersecting_primir_bam_by_position**](#sort_intersecting_primir_bam_by_position)
  - Mature miRNA extended annotations (`.gff3`); from
  [**extend_mirs_annotations**](#extend_mirs_annotations)
- **Output**
  - Mature miRNA intersections file (`.bed`); used in
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  and [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag)


#### `filter_sam_by_intersecting_mirna`

Remove alignments that do not intersect with any miRNA with
[**SAMtools**](#third-party-software-used).

> Required to efficiently classify the alignments.

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_primir**](#filter_sam_by_intersecting_primir) 
  - Mature miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect_extended_mirna)
- **Output**
  - (**Workflow output**) Alignments file, filtered (`.sam`); used in 
  [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag) and
  [**uncollapse_reads**](#uncollapse_reads)


#### `add_intersecting_mirna_tag`

Classify and add the intersecting (iso)miR to each alignment as a tag
with a [**custom script**][custom-script-iso-tag].

> In this step, the mature miRNA annotated regions are used instead of the
> extended ones. Each alignment gets an extra tag (`YW:Z`) with either the
> (iso)miR(s) it is considered to really intersect with or an empty string
> otherwise. The format of the intersecting mature miRNA species is
> `miRNA_name|5p-shift|3p-shift|CIGAR|MD|read_name`, where `5p-shift` and
> `3p-shift` are the difference between the miRNA start and end coordinates and
> the alignment's ones respectively.

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
  - Mature miRNA intersections file (`.bed`); from
  [**intersect_extended_mirna**](#intersect_extended_mirna)
- **Parameters**
  - **config_template.yaml**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended (default: 6)
- **Output**
  - Alignments file, tagged (`.sam`); used in
  [**sort_intersecting_mirna_by_feat_tag**](#sort_intersecting_mirna_by_feat_tag)
- **Examples**

```console
Example 1 | Intersecting a canoncial mature miRNA

IN miRNA annotations:
    chr19	.	miRNA	44377	44398	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160
IN SAM record:
    1-1_1	0	19	44377	255	22M	*	0	0	CTACAAAGGGAAGCACTTTCTC	*	MD:Z:22	NH:i:1	NM:i:0
NEW TAG:
	YW:Z:hsa-miR-524-5p|0|0|22M|22|CTACAAAGGGAAGCACTTTCTC


Example 2 | Intersecting an isomiR (no shifts)

IN miRNA annotations:
    chr19	.	miRNA	44377	44398	.	+	.	ID=MIMAT0002849;Alias=MIMAT0002849;Name=hsa-miR-524-5p;Derives_from=MI0003160
IN SAM record:
    1-1_1	0	19	44377	1	11M3I11M	*	0	0	CTACAAAGGGAGGTAGCACTTTCTC	*	HI:i:0	MD:Z:22	NH:i:1	NM:i:3
NEW TAG:
    YW:Z:hsa-miR-524-5p|0|0|11M3I11M|22|CTACAAAGGGAGGTAGCACTTTCTC


Example 3 | Intersecting an isomiR (no InDels nor mismatches)

IN miRNA annotations:
    chr19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786
IN SAM record:
    1-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
NEW TAG:
    YW:Z:hsa-miR-1323|0|-1|21M|21|TCAAAACTGAGGGGCATTTTC


Example 4 | Not intersecting an (iso)miR

IN miRNA annotations:
    chr19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786
IN SAM record:
    1-1_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0
NEW TAG:
    YW:Z:
```

#### `sort_intersecting_mirna_by_feat_tag`

Sort the alignments by the tag containing the classified intersecting (iso)miR
with [**SAMtools**](#third-party-software-used).

> Required for an efficient quantification.

- **Input**
  - Alignments file, tagged (`.sam`); from
  [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag)
- **Output**
  - Alignments file, tagged, sorted (`.sam`); used in
  [**quantify_mirna**](#quantify_mirna)


#### `quantify_mirna`

Tabulate alignments according to its new tag (`YW:Z`) with a
[**custom script**][custom-script-mir-quant].

> Each alignment contributes to the miRNA species in its `YW:Z` tag by `R/N`,
> where `R` is the number of collapsed reads in that alignment, and `N` is the
> number of genomic and/or transcriptomic loci it aligns to. The values of
> both, `R` and `N` are inferred from the sequence name which follows the
> format `ID-R_N`. The resulting table has a row for each mature miRNA species
> (isomiR, canonical miRNA or both) with the name format set in
> [**add_intersecting_mirna_tag**](#add_intersecting_mirna_tag) unless
> considered a canonical miRNA, which only keeps the annotated mature miRNA
> name. A miRNA species is considered to be canonical if there are no shifts
> between its start and end positions and the aligned read ones, and there
> are no mismatches nor InDels.

- **Input** 
  - Alignments file, tagged, sorted (`.sam`); from
  [**sort_intersecting_mirna_by_feat_tag**](#sort_intersecting_mirna_by_feat_tag)
- **Parameters**
  - **samples.tsv**
    - Library name; specified in the sample's table column `sample`
  - **config_template.yaml**
    - `mir_list`: miRNA features to be quantified (default: isomir, mirna
    pri-miR)
- **Output**
  - (iso)miR counts tab-delimited file; used in
  [**merge_tables**](#merge_tables)
- **Examples**

```console
Example 1 | Canonical miRNA and isomiR

IN SAM record:
    10-4_2	0	19	34627	255	21M	*	0	0	AAAGTGCTTCCTTTTAGAGGG	*	MD:Z:21	NM:i:0	NH:i:2	HI:i:1	YW:Z:hsa-miR-520b-3p|0|0|21M|21
    10-4_2	0	19	40866	255	21M	*	0	0	AAAGTGCTTCCTTTTAGAGGG	*	MD:Z:21	NM:i:0	NH:i:2	HI:i:2	YW:Z:hsa-miR-520c-3p|0|-1|21M|21

Data:
    Alignment:
        Read ID: 10
        Number of collapsed reads: 4
        Number of mapped genomic loci: 2
        Contribution: 4/2 = 2

    miRNA species:
        Tag name: hsa-miR-520b-3p|0|0|21M|21|AAAGTGCTTCCTTTTAGAGGG
        Type: Canonical
        Table name: hsa-miR-520b-3p
        Total count: 2

        Tag name: hsa-miR-520c-3p|0|-1|21M|21|AAAGTGCTTCCTTTTAGAGGG
        Type: isomiR
        Table name: hsa-miR-520c-3p|0|-1|21M|21|AAAGTGCTTCCTTTTAGAGGG
        Total count: 2

OUT table:
    ID	                                                lib_name
    hsa-miR-520b-3p         	                        2
    hsa-miR-520c-3p|0|-1|21M|21|AAAGTGCTTCCTTTTAGAGGG	2


Example 2 | Different isomiRs

IN SAM record:
    599-1_3	0	19	27804	255	20M	*	0	0	AAAGTGCTTCCTTTTAGAGG	*	MD:Z:20	NM:i:0	NH:i:3	HI:i:1	YW:Z:hsa-miR-526b-3p|1|-1|20M|20
    599-1_3	0	19	34627	255	20M	*	0	0	AAAGTGCTTCCTTTTAGAGG	*	MD:Z:20	NM:i:0	NH:i:3	HI:i:2	YW:Z:hsa-miR-520b-3p|0|-1|20M|20
    599-1_3	0	19	40866	255	20M	*	0	0	AAAGTGCTTCCTTTTAGAGG	*	MD:Z:20	NM:i:0	NH:i:3	HI:i:3	YW:Z:hsa-miR-520c-3p|0|-2|20M|20

Data:
    Alignment:
        Read ID: 599
        Number of collapsed reads: 1
        Number of mapped genomic loci: 3
        Contribution: 1/3 = 0.33

    miRNA species:
        Tag name: hsa-miR-526b-3p|1|-1|20M|20|AAAGTGCTTCCTTTTAGAGG
        Type: isomiR
        Table name: hsa-miR-526b-3p|1|-1|20M|20|AAAGTGCTTCCTTTTAGAGG
        Total count: 0.33

        Tag name: hsa-miR-520b-3p|0|-1|20M|20|AAAGTGCTTCCTTTTAGAGG 
        Type: isomiR
        Table name: hsa-miR-520b-3p|0|-1|20M|20|AAAGTGCTTCCTTTTAGAGG
        Total count: 0.33

        Tag name: hsa-miR-520c-3p|0|-2|20M|20|AAAGTGCTTCCTTTTAGAGG
        Type: isomiR
        Table name: hsa-miR-520c-3p|0|-2|20M|20|AAAGTGCTTCCTTTTAGAGG 
        Total count: 0.33
                      
OUT table:            
    ID	                                                lib_name
    hsa-miR-520b-3p|0|-1|20M|20|AAAGTGCTTCCTTTTAGAGG	0.33
    hsa-miR-520c-3p|0|-2|20M|20|AAAGTGCTTCCTTTTAGAGG    0.33
    hsa-miR-526b-3p|1|-1|20M|20|AAAGTGCTTCCTTTTAGAGG    0.33
```


#### `quantify_primir`

Tabulate alignments according to its intersecting pri-miR with a
[**custom script**][custom-script-pri-quant]

> Each alignment contributes to the pri-miR it intersects with by `R/N`, where
> `R` is the number of collapsed reads in that alignment, and `N` is the
> number of genomic and/or transcriptomic loci it aligns to. The values of
> both, `R` and `N` are inferred from the sequence name which follows the
> format `ID-R_N`. The resulting table has a row for each pri-miR with the
> name format set in [**mirna_extension**](#mirna_extension).

- **Input**
  - pri-miR intersections file (`.bed`); from
  [**intersect_extended_primir**](#intersect_extended_primir)
- **Output**
  - pri-miR counts tab-delimited file; used in
  [**merge_tables**](#merge_tables)
- **Examples**

```console
Example 1 | One single pri-miR with different alignments

IN BED records:
    19	.	miRNA_primary_transcript	27766	27788	.	+	.	ID=MI0003150;Alias=MI0003150;Name=hsa-mir-526b_-0_+0	19	27765	27788	68-2_1	255	+
    19	.	miRNA_primary_transcript	27766	27787	.	+	.	ID=MI0003150;Alias=MI0003150;Name=hsa-mir-526b_-0_+0	19	27765	27787	316-1_7	1	+
    19	.	miRNA_primary_transcript	27804	27823	.	+	.	ID=MI0003150;Alias=MI0003150;Name=hsa-mir-526b_-0_+0	19	27803	27823	599-1_3	255	+
    19	.	miRNA_primary_transcript	27805	27822	.	+	.	ID=MI0003150;Alias=MI0003150;Name=hsa-mir-526b_-0_+0	19	27804	27822	226-1_4	1	+

Alignments:
    Read ID: 68
    Number of collapsed reads: 2
    Number of mapped genomic loci: 1
    Contribution: 2/1 = 2

    Read ID: 316
    Number of collapsed reads: 1
    Number of mapped genomic loci: 7
    Contribution: 1/7 = 0.143

    Read ID: 599
    Number of collapsed reads: 1
    Number of mapped genomic loci: 3
    Contribution: 1/3 = 0.33

    Read ID: 226
    Number of collapsed reads: 1
    Number of mapped genomic loci: 4
    Contribution: 1/4 = 0.25

OUT table:            
    ID	                lib_name
    hsa-mir-526b_-0_+0	2.723


Example 2 | Different pri-miRs for a single read

IN BED records:
    19	.	miRNA_primary_transcript	40866	40886	.	+	.	ID=MI0003158;Alias=MI0003158;Name=hsa-mir-520c_-0_+0	19	40865	40886	10-4_2	255	+
    19	.	miRNA_primary_transcript	34627	34647	.	+	.	ID=MI0003155;Alias=MI0003155;Name=hsa-mir-520b_-5_+6	19	34626	34647	10-4_2	255	+

Alignment:
    Read ID: 10
    Number of collapsed reads: 4
    Number of mapped genomic loci: 2
    Contribution: 4/2 = 2

OUT table:            
    ID	                lib_name
    hsa-mir-520c_-0_+0	2
    hsa-mir-520b_-5_+6	2
```


#### `merge_tables`

Merge all the tables from the different libraries into a single one with a
[**custom script**][custom-script-merge-tab].

> The final table(s) containing the counting data from all libraries for the
> (iso)miRs and/or pri-miRs have a row per miRNA species and a column per
> sample library. If a miRNA species is not found in a certain library, its
> value is set to `NA`.

- **Input**
  - Counts tab-delimited file; from [**quantify_mirna**](#quantify_mirna)
  and/or [**quantify_primir**](#quantify_primir)
- **Parameters**
  - **cluster_schema.json**
    - `mir_list`: miRNA features to be quantified (default: isomir, mirna
    pri-mir)
- **Output**
  - (**Workflow output**) (iso)miR and/or pri-miR counts table (`.tab`)
- **Example**

```console
IN library 1
    ID                                                       lib_1
    hsa-miR-524-5p          	                             1	 
    hsa-miR-524-5p|0|0|22M|9G12|CTACAAAGGTAAGCACTTTCTC       1    
    hsa-miR-524-5p|0|0|22M|9G9C2|CTACAAAGGTAAGCACTTTATC      1    

IN library 2
    ID                                                      lib_2
    hsa-miR-524-5p          	                            1
    hsa-miR-1283	                                        1
    hsa-miR-498-3p|0|0|23M|8T14|AAAGCACCACCAGAGCTTGAAGC	    1

IN library 3
    ID   lib_3

OUT table
    ID                                                     lib_1  lib_2  lib_3
    hsa-miR-524-5p                                   	    1	    1       NA
    hsa-miR-524-5p|0|0|22M|9G12|CTACAAAGGTAAGCACTTTCTC	    1   	NA      NA
    hsa-miR-524-5p|0|0|22M|9G9C2|CTACAAAGGTAAGCACTTTATC     1   	NA      NA
    hsa-miR-1283	                                        NA	    1       NA
    hsa-miR-498-3p|0|0|23M|8T14|AAAGCACCACCAGAGCTTGAAGC	    NA	    1       NA
```


#### `uncollapse_reads`

Reverse the collapsing of reads with identical sequences as done with
[**FASTX-Toolkit**](#third-party-software-used) with a
[**custom script**][custom-script-uncollapse].

- **Input**
  - Alignments file, filtered (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
- **Output**
  - Uncollapsed aligned reads (`.sam`); used in
  [**convert_uncollapsed_reads_sam_to_bam**](#convert_uncollapsed_reads_sam_to_bam)


#### `convert_uncollapsed_reads_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file, uncollapsed (`.sam`); from
  [**filter_sam_by_intersecting_mirna**](#filter_sam_by_intersecting_mirna)
- **Output**
  - Alignments file, uncollapsed (`.bam`); used in
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)


#### `sort_uncollapsed_reads_bam_by_position`

Sort alignments by position with [**SAMtools**](#third-party-software-used).

- **Input**
  - Alignments file, uncollapsed (`.bam`); from
  [**convert_uncollapsed_reads_sam_to_bam**](#convert_uncollapsed_reads_sam_to_bam)
- **Output**
  - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `index_uncollapsed_reads_bam`

Create index BAM file with [**SAMtools**](#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
> alignments in a genomic region of interest.

- **Input**
  - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
- **Output**
  - (**Workflow output**) BAM index file (`.bam.bai`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


### Pileup workflow

#### `finish_pileup`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - (**Workflow output**) Empty text file (`.txt`)
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
  - Empty BED file (`.bed`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `compress_reference_genome`

Compress the processed genome with trimmed IDs using `bgzip` with
[**SAMtools**](#third-party-software-used). 

> Required to perform the ASCII-style alignment pileups.

- **Input**
  - Genome sequence, trimmed IDs (`.fa`); from
  [**trim_genome_seq_ids**](#trim_genome_seq_ids)
- **Output**
  - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); used in
  [**create_per_library_ascii_pileups**](#create_per_library_ascii_pileups),
  [**create_per_run_ascii_pileups**](#create_per_run_ascii_pileups) and/or
  [**create_per_condition_ascii_pileups**](#create_per_condition_ascii_pileups)


#### `create_per_library_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions across
libraries with [**ASCII-style alignment pileups**](#third-party-software-used).

> A directory containing the ASCII-style pileups is created for each
> library. If no BED file is provided, the pileups' output directories will
> only contain an empty file.

- **Input**
  - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); from
  [**compress_reference_genome**](#compress_reference_genome)
  - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
  [**map_chr_names**](#map_chr_names)
  - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - (**Workflow output**) BAM index file (`.bam.bai`); used in
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Annotated genomic regions (`.bed`); from workflow input files or
  [**create_empty_bed**](#create_empty_bed)
- **Output**
  - (**Workflow output**) Empty text file (`.txt`)


#### `create_per_run_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions for the whole
run with [**ASCII-style alignment pileups**](#third-party-software-used). 

> If no BED file is provided, the pileups' output directory will only contain
> an empty file.

- **Input**
  - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); from
  [**compress_reference_genome**](#compress_reference_genome)
  - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
  [**map_chr_names**](#map_chr_names)
  - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - (**Workflow output**) BAM index file (`.bam.bai`); used in
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Annotated genomic regions (`.bed`); from workflow input files or
  [**create_empty_bed**](#create_empty_bed)
- **Output**
  - (**Workflow output**) Empty text file (`.txt`)


#### `create_per_condition_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions across the
different library subsets if provided with
[**ASCII-style alignment pileups**](#third-party-software-used).

> **OPTIONAL RULE.** The ASCII-style pileups for each annotated region are
> made if, and only if, at least one library subset is specified in the
> [configuration file](#configuration-file). Otherwise, this rule will not be
> executed, and no output will be generated.

- **Input**
  - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); from
  [**compress_reference_genome**](#compress_reference_genome)
  - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
  [**map_chr_names**](#map_chr_names)
  - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
  [**sort_uncollapsed_reads_bam_by_position**](#sort_uncollapsed_reads_bam_by_position)
  - (**Workflow output**) BAM index file (`.bam.bai`); used in
  [**index_uncollapsed_reads_bam**](#index_uncollapsed_reads_bam)
  - Annotated genomic regions (`.bed`); from workflow input files or
  [**create_empty_bed**](#create_empty_bed)
- **Parameters**
  - **config_template.yaml**
    - `lib_dict`: Dictionary of arbitrary condition names (keys) and library
    names to aggregate alignment pileups for (values; MUST correspond to names
    in samples table) (default: None)
- **Output**
  - Empty text file (`.txt`)


[chr-maps]: <https://github.com/dpryan79/ChromosomeMappings>
[cite_neilsen]:<https://www.sciencedirect.com/science/article/pii/S0168952512001126>
[cite_saunders]: <https://pubmed.ncbi.nlm.nih.gov/17360642/>
[cite_schumauch]: <https://www.biorxiv.org/content/10.1101/2024.03.28.587190v1>
[custom-script-blocksort]: scripts/blocksort.sh
[custom-script-filter-mm]: scripts/filter_multimappers.py
[custom-script-get-lines]: scripts/get_lines_w_pattern.sh
[custom-script-gtf-bed]: scripts/gtf_exons_bed.1.1.2.R
[custom-script-iso-tag]: scripts/annotate_sam_with_bed_features.py
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
[oligomap-out]: <https://github.com/zavolanlab/oligomap#output-format>
[pub-bedtools]: <https://academic.oup.com/bioinformatics/article/26/6/841/244688>
[pub-cufflinks]: <https://doi.org/10.1038/nprot.2012.016>
[pub-cutadapt]: <https://doi.org/10.14806/ej.17.1.200>
[pub-oligomap]: <https://doi.org/10.1016/j.ymeth.2007.10.002>
[pub-samtools]: <https://doi.org/10.1093/bioinformatics/btp352>
[pub-segemehl]: <https://doi.org/10.1371/journal.pcbi.1000502>
[rule-graph]: images/rule_graph.svg

