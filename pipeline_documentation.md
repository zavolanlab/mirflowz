# MIRFLOWZ: workflow documentation

This document describes the individual steps of the workflow. For instructions
on installation and usage please see [here](README.md).

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Third-party software used](#third-party-software-used)
- [Description of workflow steps](#description-of-workflow-steps)
  - [Rule graph](#rule-graph)

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
