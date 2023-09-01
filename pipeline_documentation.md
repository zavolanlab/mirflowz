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
  - [Prepare workflow](#prepare-workflow)
  - [Map workflow](#map-workflow)
  - [Quantify workflow](#quantify-workflow)


## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
| **bedtools** | [GPLv2][license-gpl2] | _"[...] intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][code-bedtools] |
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

- tab-separated values (`.tsv` and `.csv`) file
- First row has to contain parameter names as in 
[`samples_table.csv`](test/test_files/samples_table.csv) 
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

### Prepare workflow

#### `trim_genome_seq_ids`

Trim genome sequence IDs.

- **Input**
  - Genome sequence file (`.fasta`)
- **Output**
  - Genome sequence file with trim IDs (`.fasta`); used in
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

Trim transcriptome sequence IDs.

- **Input**
  - Transcriptome sequence file (`.fasta`)
- **Output**
  - Transcriptome sequence file with trim IDs (`.fasta`); used in
  [**generate_segemehl_index_transcriptome**](#generate-segemehl-index-transcriptome),
  [**mapping_transcriptome_segemehl**](#mapping-transcriptome-segemehl) and
  [**mapping_transcriptome_oligomap**](#mapping-transcriptome-oligomap)


#### `generate_segemehl_index_transcriptome`

Generate transcriptome index for [**segemehl**](#third-party-software-used)
short read aligner.

> The transcriptome index only needs to be generated once for each combination
of genome and annotations and sample sets.

- **Input**
  - Transcriptome sequence file with trim IDs (`.fasta`); from
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)
- **Output**
  - segemehl index (`.idx`); used in 
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
  - segemehl index (`.idx`); used in 
  [**mapping_genome_segemehl**](#mapping-genome-segemehl)


#### `get_exons_gtf`

Retrieve exon annotations from genome annotations.

- **Input**
  - Genomic annotations (`.gtf`)
- **Output**
  - Exon annotations (`.gtf`); used in 
  [**convert_exons_gtf_to_bed**](#convert-exons-gtf-to-bed)

#### `convert_exons_gtf_to_bed`

Convert exon annotations `.gtf` to `.bed`.

- **Input**
  - Exon annotations (`.gtf`); from [**get_exons_gtf**](#get-exons-gtf)
- **Output**
  - Exon annotations (`.bed`); used in
  [**transcriptome_to_genome_maps**](#transcriptome-to-genome-maps)


#### `create_genome_header`

Create `SAM` header for the genome with 
[**samtools**](#third-party-software-used).

> Required by [samtools](#third-party-software-used) to work with the alignment
file.

- **Input**
  - Genome sequence file with trim IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - Genome header (`.sam`); used in [add_header_all_maps](#add-header-all-maps)


#### `map_chr_names`

Map UCSC-like chromosome names with Ensembl-like ones in miRNA annotation.

> Required by [bedtools](#third-party-software) to intersect alignments with
miRNA annotations. Several mapping tables are available [here][chr-maps].

- **Input**
  - miRNA annotations (`.gff3`)
  - Tab-separated mappings table (`.txt`)
- **Output**
  - miRNA annotations with mapped genes(`.gff3`); used in 
  [**extend_mirs_annotations**](#extend-mirs-annotations)


#### `create_index_genome_fasta`

Create a `FASTA` index for the genome with 
[**samtools**](#third-party-software-used).

- **Input**
  - Genome sequence file with trim IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - `FASTA` genome index (`.fa.fai`); used in
  [**extract_chr_len**](#extract-chr-len)


#### `extract-chr-len`

Extract chromosome(s) length from the genome sequence.

- **Input**
  - `FASTA` genome index (`.fa.fai`); from
  [**create_index_genome_fasta**](#create-index-genome-fasta)
- **Output**
  - Tab-separated table mapping chromosome name(s) and length(s) (`.txt`); used
  in [**extend_mirs_annotations**](#extend-mirs-annotations)


#### `extend-mirs-annotations`

Extend miRNA annotations and split the file by feature.

> Mature miRNA regions are extended on both sides to account for isomiR species
with shifted start and/or end positions. If required, pri-miR loci are also 
extended to accommodate the new miRNA coordinates. 

- **Input**
  - miRNA annotations with mapped genes(`.gff3`); from
  [**map_chr_names**](#map-chr-names)
- **Parameters**
  - **config_schema.json**
    - `extension`: Number of nucleotides by which mature miRNA annotated
    regions are extended (default: 6)
- **Output**
  - Primary miRNA transcript (pri-mir) extended annotation (`.gff3`); used in
  [**intersect_extended_primir**](#intersect-extended-primir)
  - Mature miRNA (miRNA) extended annotation (`.gff3`); used in
  [**intersect_extended_mirna**](#intersect-extended-mirna)

### Map workflow

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
  - Reads file (`.fastq`); from [**start**](#start)
- **Parameters**
  - **config_schema.json**
    - `q_value`: Minimum Q (Phred) score to keep (default 10)
    - `p_value`: Minimum % of bases that must have a Q (Phred) quality
    (default 50)
- **Output**
  - Reads file filtered (`.fastq`); used in
  [**fastq_to_fasta**](#fastq-to-fasta)


#### `fastq_to_fasta`

Convert reads file from `.fastq` to `.fa` with 
[**fastx_toolkit**](#third-party-software-used).

- **Input**
  - Reads file (`.fastq`); from
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
  - Reads file (`.fa`); used in
  [**remove_adapters**](#remove-adapters)


#### `remove_adapters`

Trim adapters and `N` bases at either end. Filter reads by minimum length and
number of inner `N` bases with [**cutadapt**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from [**format_fasta**](#format-fasta)
- **Parameters**
  - **samples.csv**
    - Adapter to be removed; specify in sample table column `adapter`
  - **config_schema.json**
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
  - Collapsed and rename reads file; used in
  [**filter_fasta_for_oligomap**](#filter-fasta-for-oligomap),
  [**map_genome_segemehl**](#map-genome-segemehl) and
  [**map_transcriptome_segemehl**](#map-transcriptome-segemehl)


#### `map_genome_segemehl`

Align short reads to reference genome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from
  [**collapse_identical_reads**](#collapse-identical-reads)
  - Genome sequence (`.fa`); from 
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
  - segemehl index (`idx`); from
  [**generate_segemehl_index_genome**](#generate-segemehl-index-genome)
- **Output**
  - Aligned reads file (`.sam`); used in
  [**merge_genome_maps**](#merge-genome-maps)


#### `map_transcriptome_segemehl`

Align short reads to reference transcriptome with 
[**segemehl**](#third-party-software-used).

- **Input**
  - Reads file (`.fa`); from
  [**collapse_identical_reads**](#collapse-identical-reads)
  - Transcriptome sequence (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)
  - segemehl index (`idx`); from
  [**generate_segemehl_index_transcriptome**](#generate-segemehl-index-transcriptome)
- **Output**
  - Aligned reads file (`.sam`); used in
  [**merge_transcriptome_maps**](#merge-transcriptome-maps)


#### `filter_fasta_for_oligomap`

Filter reads by length.

- **Input**
  - Reads file (`.fa`); from
  [**collapse_identical_reads**](#collapse-identical-reads)
- **Parameters**
  - **config_schema.json**
    - `max_length_reads`: Maximum length of processed reads to map with
    [**oligomap**](#third-party-software-used)
- **Output**
  - Reads file (`.fa`); used in [**map_genome_oligomap**](#map-genome-oligomap)
  and [**map_transcriptome_oligomap**](#map-transcriptome-oligomap)


#### `map_genome_oligomap`

Align short reads to reference genome with
[**oligomap**](#third-party-software.used).

- **Input**
  - Reads file (`.fa`); from
  [**filter_fasta_for_oligomap**](#filter-fasta-for-oligomap)
  - Genome sequence (`.fa`); from 
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
- **Output**
  - Aligned reads file (`.fa`); used in
  [**sort_genome_oligomap**](#sort-genome-oligomap)
  - Alignment report (`.txt`); used in
  [**sort_genome_oligomap**](#sort-genome-oligomap)


#### `sort_genome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name.

- **Input**
  - Aligned reads file (`.fa`); from
  [**map_genome_oligomap**](#map-genome-oligomap)
  - Alignment report (`.txt`); from
  [**map_genome_oligomap**](#map-genome-oligomap)
- **Output**
  - Aligned sorted reads file (`.fa`); used in
  [**convert_genome_to_sam_oligomap](#convert-genome-to-sam-oligomap)
  - Alignment sorted report (`.txt`); used in
  [**convert_genome_to_sam_oligomap](#convert-genome-to-sam-oligomap)


#### `convert_genome_to_sam_oligomap`

Convert aligned reads `.fa` to `.sam` and filter alignments by number of hits.

- **Input**
  - Aligned reads file (`.fa`); from
  [**sort_genome_oligomap**](#sort-genome-oligomap)
  - Alignment report (`.txt`); from
  [**sort_genome_oligomap**](#sort-genome-oligomap)
- **Parameters**
  - **config_schema.json**
    - `nh`: Maximum number of hits an alignment can have to be kept
- **Output**
  - Aligned reads (`.sam`); used in [**merge_genome_maps**](#merge-genome-maps)


#### `map_transcriptome_oligomap`

Align short reads to reference transcriptome with
[**oligomap**](#third-party-software.used).

- **Input**
  - Reads file (`.fa`); from
  [**filter_fasta_for_oligomap**](#filter-fasta-for-oligomap)
  - Transcriptome sequence (`.fa`); from 
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)
- **Output**
  - Aligned reads file (`.fa`); used in
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)
  - Alignment report (`.txt`); used in
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)


#### `sort_transcriptome_oligomap`

Sort [**oligomap**](#third-party-software-used) alignments by query name.

- **Input**
  - Aligned reads file (`.fa`); from
  [**map_transcriptome_oligomap**](#map-transcriptome-oligomap)
  - Alignment report (`.txt`); from
  [**map_transcriptome_oligomap**](#map-transcriptome-oligomap)
- **Output**
  - Aligned sorted reads file (`.fa`); used in
  [**convert_transcriptome_to_sam_oligomap](#convert-transcriptome-to-sam-oligomap)
  - Alignment sorted report (`.txt`); used in
  [**convert_transcriptome_to_sam_oligomap](#convert-transcriptme-to-sam-oligomap)


#### `convert_transcriptome_to_sam_oligomap`

Convert aligned reads `.fa` to `.sam` and filter alignments by number of hits.

- **Input**
  - Aligned reads file (`.fa`); from
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)
  - Alignment report (`.txt`); from
  [**sort_transcriptome_oligomap**](#sort-transcriptome-oligomap)
- **Parameters**
  - **config_schema.json**
    - `nh`: Maximum number of hits an alignment can have to be kept
- **Output**
  - Aligned reads (`.sam`); used in
  [**merge_transcriptome_maps**](#merge-transcriptome-maps)


#### `merge_genome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) genome alignments.

- **Input**
  - Aligned reads (`.sam`); from
  [**map_genome_segemehl**](#map-genome-segemehl)
  - Aligned reads (`.sam`); from
  [**convert_genome_to_sam_oligomap**](#convert-genome-to-sam-oligomap)
- **Output**
  - Aligned reads (`.sam`); used in
  [**filter_genome_by_nh**](#filter-genome-by-nh)


#### `merge_transcriptome_maps`

Concatenate [**segemehl**](#third-party-software-used) and
[**oligomap**](#third-party-software-used) transcriptome alignments.

- **Input**
  - Aligned reads (`.sam`); from
  [**map_transcriptome_segemehl**](#map-transcriptome-segemehl)
  - Aligned reads (`.sam`); from
  [**convert_transcriptome_to_sam_oligomap**](#convert-transcriptome-to-sam-oligomap)
- **Output**
  - Aligned reads (`.sam`); used in
  [**filter_transcriptome_by_nh**](#filter-transcriptome-by-nh)

EXPANDING SOON.

### Quantify workflow

COMING SOON.

[chr-maps]: <https://github.com/dpryan79/ChromosomeMappings>
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
[pub-cufflinks]: <https://doi.org/10.1038/nprot.2012.016>
[pub-cutadapt]: <https://doi.org/10.14806/ej.17.1.200>
[pub-oligomap]: <https://doi.org/10.1016/j.ymeth.2007.10.002 >
[pub-samtools]: <https://doi.org/10.1093/bioinformatics/btp352>
[pub-segemehl]: <https://doi.org/10.1371/journal.pcbi.1000502>
[rule-graph]: images/rule_graph.svg
