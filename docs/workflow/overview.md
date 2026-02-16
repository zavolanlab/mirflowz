# Workflow Overview

_MIRFLOWZ_ consists of a main `Snakefile` and four functional modules. In the
`Snakefile`, the configuration file is validated and various modules are
imported. In addition, a handler for both, a successful and a failed run are
set. If the workflow finishes without any errors, all the intermediate files
are removed, otherwise, a log file is created.

<div align="center">
    <img src=../../images/workflow_overview.png>
</div>

The modules [(1)](#prepare-module) process the genome resources,
[(2)](#map-module) map and [(3)](#quantify-module), and
[(4)](#ascii-style-alignment-pileups-module) are describe in detail below.

## Prepare module

_MIRFLOWZ_ initially processes and indexes the genome resources provided by
the user. The regions corresponding to mature miRNAs are extended by a fixed
but user-adjustable number of nucleotides on both sides to accommodate isomiR
species with shifted start and/or end positions. If necessary, pri-miR loci
are extended to adjust to the new miRNA coordinates.

<div align="center">
    <img width="90%" src=../../images/mirna_extension.png>
</div>

In addition, to account for the different genomic locations a miRNA sequence
can be annotated, the name of these sequences are modified to have the format
`SPECIES-mir-NUMBER[LETTER]-#` for pri-miRs, and
`SPECIES-miR-NUMBER[LETTER]-#-ARM` or `SPECIES-miR-NUMBER[LETTER]-#` for mature
miRNAs with both or just one arm respectively, where:

- `#` is the paralog number (replica/locus index), included when multiple loci
  express the same or similar miRNAs
- `LETTER` denotes a sequence variant of mature miRNA (paralogous variant with
  similar but not identical sequences)


??? tip "Wanna see a more technical step by step explanation?"

    For a detailed description of this module, with examples and explanations
    for each rule, see the [prepare module page](./modules/prepare.md).

## Map module

The user-provided short-read small RNA-seq libraries undergo quality filtering
(skipped if libraries are provided in FASTA rather than FASTQ), followed by
adapter removal. The resulting reads are independently mapped to both the
genome and transcriptome using two distinct aligners: [Segemehl][segemehl] and
our in-house tool [Oligomap][oligomap].

Segemehl implements a fast heuristic strategy that returns the alignment(s)
with the smallest edit distance. Oligomap, on the other hand, implements a
slower and more restricted approach that reports all the alignments with an
edit distance of at most 1. The combination of the fast and flexible results
and the strict selection ensures results with a higher fidelity than if only
one of the tools was to be used.

Two merging steps are done in order to have all the alignments in a single
file.

1. The transcriptome and the genome mappings from both aligners are fused and
   only those alignments with a smaller NH than the one provided are kept.

2. Transcriptomic coordinates are turned into genomic ones and alignments are
   combined into a single file.

<div align="center">
    <img src=../../images/merging_strategy.png>
</div>

Duplicate alignments resulting from the partially redundant mapping strategy
are discarded and only the best alignments for each read are retained (_i.e._
the ones with the smallest edit distance. In addition, and due to the
alignment's aggregation, a second filtering according to the new NH is
performed: if a read has been aligned beyond a specified threshold, it is
removed due to (1) performance reasons as the file size can rapidly increase,
and (2) the fact that each read contributes to each count `1/n` where `n` is
the number of genomic loci it aligns to and a large `n` makes the contribution
negligible.

A final filter is made to further increase the classification accuracy and
reduce the amount of multimappers (defined here as alignments of the same read
aligning to different genomic loci with the same edit distance). Given that
isomiRs are known to contain more indels than mismatches when compared to the
canonical sequence they come from, as demonstrated by
[Neilsen et al. (2012)][cite_neilsen], [Saunders et al. (2017)][cite_saunders],
and [Schmauch et al. (2024)][cite_schmauch], only those multimappers that
contain a higher or equal amount of indels compared to mismatches are retained.
Note that some multimappers might still be present if the number of indels is
the same across alignments.

??? tip "Wanna see a more technical step by step explanation?"

    For a detailed description of this module, with examples and explanations
    for each rule, see the [map module page](./modules/map.md).

## Quantify module

!!! info "isomiRs notation"

    A sequence is considered to be an isomiR if it has a shift on either end,
    an indel and/or a mismtach on its sequence when compared to the canonical
    miRNA it maps and intersects with.

    _MIRFLOWZ_ employs an unambiguous notation to classify isomiRs using the
    format `miRNA_name|5p-shift|3p-shift|CIGAR|MD|READ_SEQ`, where `5p-shift`
    and `3p-shift` represent the differences between the annotated mature miRNA
    start and end positions and those of the read alignment, respectively.

The filtered alignments are subsequently intersected with the user-provided,
pre-processed miRNA annotation files using [bedtools][bedtools]. Each
alignment is classified according to the miRNA species it fully intersects with
in order to do the counts.

<div align="center">
    <img src=../../images/mirna_intersection.png>
</div>

Counts are tabulated separately for reads consistent with either miRNA
precursors, mature miRNA and/or isomiRs, and all library counts are fused into
a single table. Note that an alignment is only counted towards a given miRNA
(or isomiR) species if one of its alignments fully falls within the (previously
extended) locus annotated for that miRNA. Specifically, reads contribute with
`1/n` for each miRNA for which that is the case, where `n` is the total number
of genomic loci the read aligns to. Under this criterion, the precursor counts
contain reads that intersect with its mature arm(s), its hairpin sequence
and/ir the whole precursor itself

<div align="center">
    <img src=../../images/read_contribution.png>
</div>


??? tip "Wanna see a more technical step by step explanation?"

    For a detailed description of this module, with examples and explanations
    for each rule, see the [quantify module page](./modules/quantify.md).

## ASCII-style alignment pileups module

Finally, to visualize the distribution of read alignments around miRNA loci,
ASCII-style alignment pileups are optionally generated for user-defined regions
of interest.

```console
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	test-mir
....>>>>>>>>>>>>>>>>>>>>>>.....................................................	test-mir-5p
.......................................................>>>>>>>>>>>>>>>>>>>>>...	test-mir-3p
GGGATGAGGTAGTAGGTTGTATAGTTTTAGGGTCACACCCACCACTGGGAGATAACTATACAATCTACTGTCTTTCCTA	test_ref:3618-3696:+
ACCATGAGGTAGTAGGTTGTATAGTT.....................................................	1
..CATGAGGTAGTAGGTTGTATAGTT.....................................................	10
..GA-GAGGTAGTAGGTTGTATAGTT.....................................................	2
...A-GAGGTAGTAGGTTGTATAGTT.....................................................	19
....TGAGGTAGTAGGTTGTATAGTT.....................................................	17
.....GAGGTAGTAGGTTGTATAGTT.....................................................	33
......AGGTAGTAGGTTGTATAGTT.....................................................	9
......AGGTAGTAGGTTGTATAGTTT....................................................	2
.......GGTAGTAGGTTGTATAGTT.....................................................	7
..................................................GATAACTATACAATCTACTGTCTT.....	1
.....................................................AACTATACAATCTACT..........	1
.......................................................CTATACAATCTACTGTCTTTCT..	28
.......................................................CTATACAATCTACTGTCTTTC-T.	22
.......................................................CTATACAATCTACTGTCTTTCC..	19
.......................................................CTATACAATCTACTGTCTTTC...	12
.......................................................CTATACAATCTACTGTCTTTCTT.	2
.......................................................CTATACAATCTACTGTC.......	1
.......................................................CTATACAATCTACTGTCTT.....	1
.......................................................CTATACAATCTACTGTCTTTCG..	1
........................................................TATACAATCTACTGTCTTTCT..	4
........................................................TATACAATCTACTGTCTTTC-T.	4
........................................................TATACAATCTACTGTCTTTC...	2
........................................................TATACAATCTACTGTCTTTCC..	1
........................................................TATACAATCTACTGTCTTTCCT.	1
```

??? tip "Wanna see a more technical step by step explanation?"

    For a detailed description of this module, with examples and explanations
    for each rule, see the [pileups module page](./modules/pileups.md).
