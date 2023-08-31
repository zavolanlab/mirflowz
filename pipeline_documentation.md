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
  [**exons_gtf_to_bed**](#exons-gtf-to-bed)

#### `exons_gtf_to_bed`

Convert the 

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

#### `finish_prepare`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Genome sequence file with trim IDs (`.fasta`); from
  [**trim_genome_seq_ids**](#trim-genome-seq-ids)
  - Transcriptome sequence file with trim IDs (`.fasta`); from
  [**trim_transcriptome_seq_ids**](#trim-transcriptome-seq-ids)

### Map workflow


### Quantify workflow



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
