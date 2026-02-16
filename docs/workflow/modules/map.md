# Map Module

<!-- Change? -->

Welcome to the technical documentation for the map module. This page aims
to be a detailed documentation of each rule within the module by stating its
inputs, outputs (and how they relate to other rules), configurable parameters,
and the software used. Moreover, when needed, there will be explanations and
examples of what that particular rule does.

The schema below is a visual representation of the individual module steps
and how they are related.

<div align="center">
    <img src=../../../images/map.png>
</div>

## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

|              **Name**             |          **License**          |                                                                                                                                       **Tag line**                                                                                                                                      |                                   **More info**                                  |
|:---------------------------------:|:-----------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------:|
| **cutadapt**                      | [MIT][license-mit]            | _"[...] finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads"_                                                                                                                                    | [code][code-cutadapt] / [manual][docs-cutadapt] / [publication][pub-cutadapt]    |
| **FASTX-Toolkit**                 | [AGPL-3.0][license-agpl3]     | _"[...] collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing"_                                                                                                                                                                                              | [code][code-fastx] / [manual][docs-fastx]                                        |
| **Oligomap**                      | [GPLv3][license-gpl3]         | _"[...] fast identification of nearly-perfect matches of small RNAs in sequence databases. It allows to exhaustively identify all the perfect and 1-error (where an error is defined to be a mismatch, insertion or deletion) matches of large sets of small RNAs to target sequences"_ | [code][code-oligomap] / [publication][pub-oligomap]                              |
| **SAMtools**                      | [MIT][license-mit]            | _"[...] suite of programs for interacting with high-throughput sequencing data"_                                                                                                                                                                                                        | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools]    |
| **segemehl**                      | [GPLv3][license-gpl3]         | _"[...] map short sequencer reads to reference genomes"_                                                                                                                                                                                                                                | [manual][docs-segemehl] / [publication][pub-segemehl]                            |


## Configuration file

Some parameters within the workflow can be modified. Refer to the
[configuration template][mirflowz_conf] for a detailed explanation of each
option.

## Map Workflow

### `finish_map`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

=== "Input"

    BAM index file (`.bam.bai`); from
    [**index_all_alns_bam**](map.md#index_all_alns_bam)


### `start`

=== "Input"

    (**Workflow input**) miRNA sequencing library, `gzip`ed
    (`.fa.gz`/`.fasta.gz` or `fq.gz`/`.fastaq.gz`)

=== "Output"

    miRNA sequencing library, copied, renamed (`.fa`, `.fastq`); used in
    [**fasta_quality_filter**](map.md#fastq_quality_filter) or
    [**format_fasta**](map.md#format_fasta)


### `fastq_quality_filter`

Conduct quality control for reads library with
[**fastx_toolkit**](map.md#third-party-software-used).

=== "Input"

    miRNA sequencing library, copied, renamed (`.fastq`); from
    [**start**](map.md#start)

=== "Parameters"

    - **config_template.yaml**
        - `q_value`: Minimum Q (Phred) score to keep (default: 10)
        - `p_value`: Minimum % of bases that must have a Q (Phred) quality
          (default: 50)

=== "Output"

    miRNA sequencing library, filtered (`.fastq`); used in
    [**fastq_to_fasta**](map.md#fastq_to_fasta)


### `fastq_to_fasta`

Convert reads file from FASTQ to FASTA with
[**fastx_toolkit**](map.md#third-party-software-used).

> Sequence identifiers are renamed to numbers.

=== "Input"

    miRNA sequencºing library, filtered (`.fastq`); from
    [**fastq_quality_filter**](map.md#fastq_quality_filter)

=== "Output"

    miRNA sequencing library (`.fa`); used in
    [**format_fasta**](map.md#format_fasta)


### `format_fasta`

Format read's sequences to appear on a single line with
[**fastx_toolkit**](map.md#third-party-software-used).

=== "Input"

    miRNA sequencing library (`.fa`); from [**start**](map.md#start) or
    [**fastq_to_fasta**](map.md#fastq_to_fasta)

=== "Output"

    miRNA sequencing library, formatted (`.fasta`); used in
    [**remove_adapters**](map.md#remove_adapters)


### `remove_adapters`

Trim 3' adapters and `N` bases at either end. Filter reads by minimum length
and number of inner `N` bases with [**cutadapt**](map.md#third-party-software-used).

=== "Input"

    miRNA sequencing library, formatted (`.fasta`); from
    [**format_fasta**](map.md#format_fasta)

=== "Parameters"

    - **samples.tsv**
        - Adapter to be removed; specified in the sample's table column
        `adapter`
    - **config_template.yaml**
        - `error_rate`: Fraction of allowed errors in the matched adapters
        (default: 0.1)
        - `overlap`: Minimum overlap length between adapter and read to trim
        the bases (default: 3)
        - `minimum_length`: Minimum length for a processed read to be kept
        (default: 15)
        - `max_n`: Maximum number of inner `N` bases for a processed read to
        be kept (default: 0)

=== "Output"

    miRNA sequencing library, filtered, without adapters (`.fasta`); used in
    [**collapse_identical_reads**](map.md#collapse_identical_reads)


### `collapse_identical_reads`

Collapse and rename identical reads
[**fastx_toolkit**](map.md#third-party-software-used).

> Sequences are renamed in the format `R-N`, where `R` is the assigned number
> to the unique entry, and `N` is the amount of identical sequences within the
> library collapsed in it.

=== "Input"

    miRNA sequencing library, filtered, without adapters (`.fasta`); from
    [**remove_adapters**](map.md#remove_adapters)

=== "Output"

    miRNA sequencing library, collapsed, renamed (`.fasta`); used in
    [**filter_fasta_for_oligomap**](map.md#filter_fasta_for_oligomap),
    [**map_genome_segemehl**](map.md#map_genome_segemehl) and
    [**map_transcriptome_segemehl**](map.md#map_transcriptome_segemehl)


### `map_genome_segemehl`

Align short reads to reference genome with
[**segemehl**](map.md#third-party-software-used).

=== "Input"

    - miRNA sequencing library, collapsed, renamed (`.fasta`); from
    [**collapse_identical_reads**](map.md#collapse_identical_reads)
    - Genome sequence, trimmed IDs (`.fa`); from
    [**trim_genome_seq_ids**](prepare.md#trim_genome_seq_ids)
    - segemehl genome index (`.idx`); from
    [**generate_segemehl_index_genome**](prepare.md#generate_segemehl_index_genome)

=== "Output"

    Alignments file (`.sam`); used in
    [**merge_genome_maps**](map.md#merge_genome_maps)


### `map_transcriptome_segemehl`

Align short reads to reference transcriptome with
[**segemehl**](map.md#third-party-software-used).

=== "Input"

    - miRNA sequencing library, collapsed, renamed (`.fasta`); from
    [**collapse_identical_reads**](map.md#collapse_identical_reads)
    - Transcriptome sequence, trimmed IDs (`.fa`); from
    [**trim_transcriptome_seq_ids**](prepare.md#trim_transcriptome_seq_ids)
    - segemehl transcriptome index (`.idx`); from
    [**generate_segemehl_index_transcriptome**](prepare.md#generate_segemehl_index_transcriptome)

=== "Output"

    Alignments file (`.sam`); used in
    [**merge_transcriptome_maps**](map.md#merge_transcriptome_maps)


### `filter_fasta_for_oligomap`

Filter reads by length with a [**custom script**][custom-script-validation].

> Required for an optimal mapping speed. Oligomap is specifically written for
> short reads. Therefore, reads with more bases than the default maximum (30
> nts) makes the mapping slower.

=== "Input"

    miRNA sequencing library, collapsed, renamed (`.fasta`); from
    [**collapse_identical_reads**](map.md#collapse_identical_reads)

=== "Parameters"

    - **config_template.yaml**
        - `max_length_reads`: Maximum length of processed reads to be mapped
           with [**oligomap**](map.md#third-party-software-used) (default: 30)

=== "Output"

    miRNA sequencing library, collapsed, filtered (`.fasta`); used in
    [**map_genome_oligomap**](map.md#map_genome_oligomap) and
    [**map_transcriptome_oligomap**](map.md#map_transcriptome_oligomap)


### `map_genome_oligomap`

Align short reads to reference genome with
[**oligomap**](map.md#third-party-software-used).

> Refer to Oligomap's [**Output format section**][oligomap-out] for a specific
> explanation and examples on the output format.

=== "Input"

    - miRNA sequencing library, collapsed, filtered (`.fasta`); from
    [**filter_fasta_for_oligomap**](map.md#filter_fasta_for_oligomap)
    - Genome sequence, trimmed IDs (`.fa`); from
    [**trim_genome_seq_ids**](prepare.md#trim_genome_seq_ids)

=== "Output"

    - Alignments file (`.oligomap`); used in
    [**sort_genome_oligomap**](map.md#sort_genome_oligomap)
    - Alignments report (`.txt`)


### `sort_genome_oligomap`

Sort [**oligomap**](map.md#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

=== "Input"

    Alignments file (`.oligomap`); from
    [**map_genome_oligomap**](map.md#map_genome_oligomap)

=== "Output"

    Alignments file, sorted (`.oligomap`); used in
    [**convert_genome_to_sam_oligomap**](map.md#convert_genome_to_sam_oligomap)


### `convert_genome_to_sam_oligomap`

Convert aligned reads `.oligomap` to `.sam` and filter alignments by number of
hits with a [**custom script**][custom-script-oligo-sam].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

=== "Input"

    Alignments file, sorted (`.oligomap`); from
    [**sort_genome_oligomap**](map.md#sort_genome_oligomap)

=== "Parameters"

    - **config_template.yaml**
        - `nh`: Maximum number of mappings per read to be kept (default: 100)

=== "Output"

    Alignments file, filtered (`.sam`); used in
    [**merge_genome_maps**](map.md#merge_genome_maps)


### `map_transcriptome_oligomap`

Align short reads to reference transcriptome with
[**oligomap**](map.md#third-party-software-used).

> Refer to Oligomap's [Output format section][oligomap-out] for a specific
> explanation and examples on the output format.

=== "Input"

    - miRNA sequencing library, collapsed, filtered (`.fasta`); from
    [**filter_fasta_for_oligomap**](map.md#filter_fasta_for_oligomap)
    - Transcriptome sequence, trimmed IDs (`.fa`); from
    [**trim_transcriptome_seq_ids**](prepare.md#trim_transcriptome_seq_ids)

=== "Output"

    - Alignments file (`.oligomap`); used in
    [**sort_transcriptome_oligomap**](map.md#sort_transcriptome_oligomap)
    - Alignments report (`.txt`)


### `sort_transcriptome_oligomap`

Sort [**oligomap**](map.md#third-party-software-used) alignments by query name
with a [**custom script**][custom-script-blocksort].

=== "Input"

    Alignments file (`.oligomap`); from
    [**map_transcriptome_oligomap**](map.md#map_transcriptome_oligomap)

=== "Output"

    Alignments file, sorted (`.oligomap`); used in
    [**convert_transcriptome_to_sam_oligomap**](map.md#convert_transcriptome_to_sam_oligomap)


### `convert_transcriptome_to_sam_oligomap`

Convert aligned reads `.oligomap` to `.sam` and filter alignments by number of
hits with a [**custom script**][custom-script-oligo-sam].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

=== "Input"

    Alignments file, sorted (`.oligomap`); from
    [**sort_transcriptome_oligomap**](map.md#sort_transcriptome_oligomap)

=== "Parameters"

    - **config_template.yaml**
        - `nh`: Maximum number of mappings per read to be kept (default: 100)

=== "Output"

    Alignments file, filtered (`.sam`); used in
    [**merge_transcriptome_maps**](map.md#merge_transcriptome_maps)


### `merge_genome_maps`

Concatenate [**segemehl**](map.md#third-party-software-used) and
[**oligomap**](map.md#third-party-software-used) genome alignments.

=== "Input"

    - Alignments file (`.sam`); from
    [**map_genome_segemehl**](map.md#map_genome_segemehl)
    - Alignments file, filtered (`.sam`); from
    [**convert_genome_to_sam_oligomap**](map.md#convert_genome_to_sam_oligomap)

=== "Output"

    Alignments file (`.sam`); used in
    [**filter_genome_by_nh**](map.md#filter_genome_by_nh)


### `merge_transcriptome_maps`

Concatenate [**segemehl**](map.md#third-party-software-used) and
[**oligomap**](map.md#third-party-software-used) transcriptome alignments.

=== "Input"

    - Alignments file (`.sam`); from
    [**map_transcriptome_segemehl**](map.md#map_transcriptome_segemehl)
    - Alignments file, filtered (`.sam`); from
    [**convert_transcriptome_to_sam_oligomap**](map.md#convert_transcriptome_to_sam_oligomap)

=== "Output"

    Alignments file (`.sam`); used in
    [**filter_transcriptome_by_nh**](map.md#filter_transcriptome_by_nh)


### `filter_genome_by_nh`

Filter merged genome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

=== "Input"

    Alignments file (`.sam`); from
    [**merge_genome_maps**](map.md#merge_genome_maps)

=== "Parameters"

    - **config_template.yaml**
        - `nh`: Maximum number of mappings per read to be kept (default: 100)

=== "Output"

    Alignments file, filtered (`.sam`); used in
    [**remove_header_genome_mappings**](map.md#remove_header_genome_mappings)


### `filter_transcriptome_by_nh`

Filter merged transcriptome alignments by the number of hits with a
[**custom script**][custom-script-nh-filter].

> If a read has been aligned beyond a specified threshold, it is removed due
> to (1) performance reasons as the file size can rapidly increase, and (2)
> the fact that each read contributes to each count `1/N` where `N` is the
> number of genomic loci it aligns to and a large `N` makes the contribution
> negligible.

=== "Input"

    Alignments file (`.sam`); from
    [**merge_transcriptome_maps**](map.md#merge_transcriptome_maps)

=== "Parameters"

    - **config_template.yaml**
        - `nh`: Maximum number of mappings per read to be kept (default: 100)

=== "Output"

    Alignments file, filtered (`.sam`); used in
    [**remove_header_transcriptome_mappings**](map.md#remove_header_transcriptome_mappings)


### `remove_header_genome_mappings`

Remove the SAM header of the genome alignments file with
[**SAMtools**](map.md#third-party-software-used).

=== "Input"

    Alignments file (`.sam`); from
    [**filter_genome_by_nh**](map.md#filter_genome_by_nh)

=== "Output"

    Alignments file, without SAM header (`.sam`); used in
    [**merge_all_maps**](map.md#merge_all_maps)


### `remove_header_transcriptome_mappings`

Remove the SAM header of the transcriptome alignments file with
[**SAMtools**](map.md#third-party-software-used).

=== "Input"

    Alignments file (`.sam`); from
    [**filter_transcriptome_by_nh**](map.md#filter_transcriptome_by_nh)

=== "Output"

    Alignments file, without SAM header (`.sam`); used in
    [**transcriptome_to_genome_maps**](map.md#transcriptome_to_genome_maps)


### `transcriptome_to_genome_maps`

Convert the alignments' transcriptome coordinates to genomic ones with a
[**custom script**][custom-script-sam-trx].

=== "Input"

    - Alignments file, without SAM header (`.sam`); from
    [**remove_header_transcriptome_mappings**](map.md#remove_header_transcriptome_mappings)
    - Exon annotations (`.bed`); from
    [**convert_exons_gtf_to_bed**](prepare.md#convert_exons_gtf_to_bed)

=== "Output"

    Alignments file, without SAM header (`.sam`); used in
    [**merge_all_maps**](map.md#merge_all_maps)


### `merge_all_maps`

Concatenate the four alignment files into a single one.

=== "Input"

    Alignments files, without SAM header (`.sam`); from
    [**remove_header_genome_mappings**](map.md#remove_header_genome_mappings) and
    [**transcriptome_to_genome_maps**](map.md#transcriptome_to_genome_maps)

=== "Output"

    Alignments file, without SAM header (`.sam`); used in
    [**add_header_all_maps**](map.md#add_header_all_maps)


### `add_header_all_maps`

Add the SAM header to the aligned reads merged file with
[**SAMtools**](map.md#third-party-software-used).

=== "Input"

    Alignments file, without SAM header (`.sam`); from
    [**merge_all_maps**](map.md#merge_all_maps)

=== "Output"

    Alignments file (`.sam`); used in [**sort_maps_by_id**](map.md#sort_maps_by_id)


### `sort_maps_by_id`

Sort alignments by reads ID with [**SAMtools**](map.md#third-party-software-used).

=== "Input"

    Alignments file (`.sam`); from
    [**add_header_all_maps**](map.md#add_header_all_maps)

=== "Output"

    Alignments file, sorted (`.sam`); used in
    [**remove_inferiors**](map.md#remove_inferiors)


### `remove_inferiors`

Remove duplicate and inferior alignments with a
[**custom script**][custom-script-remove-dup].

> Alignments are considered to be duplicates if having identical entries for
> the fields `QNAME`, `FLAG`, `RNAME`, `POS` and `CIGAR`.
> Alignments are considered to be inferiors if having the same `QNAME` and
> a bigger edit distance than the smaller one within the group. The tags `NH`
> (number of hits) and `HI` (query hit index) are updated accordingly.

=== "Input"

    Alignments file, sorted (`.sam`); from
    [**sort_maps_by_id**](map.md#sort_maps_by_id)

=== "Output"

    Alignments file, filtered (`.sam`); used in
    [**filter_by_indels**](map.md#filter_by_indels)

=== "Remove Duplicates Example"

    ```console
    Example 1 | Remove duplicates

    IN:
        1-2	0	19	44414	1	21M	*	0	0	GAAGGCGCTTCCCTTTGGAGT	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:21	RG:Z:A1	YZ:Z:0
        1-2	0	19	44414	255	21M	*	0	0	GAAGGCGCTTCCCTTTGGAGT	*	NM:i:0	MD:Z:21	NH:i:1
    OUT:
        1-2	0	19	44414	255	21M	*	0	0	GAAGGCGCTTCCCTTTGGAGT	*	MD:Z:21	NH:i:1	NM:i:0
    ```


=== "Remove Inferiors Examples"

    ```console
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


### `filter_by_indels`

Filter multimappers favoring indels over mismatches with a
[**custom script**][custom-script-filter-mm].

> Given that indels are more frequent in miRNAs than mismatches, as
> demonstrated by [Saunders et al. (2017)][cite_saunders],
> [Neilsen et al. (2012)][cite_neilsen] and
> [Schmauch et al. (2024)][cite_schmauch], only those "multimappers" (defined
> here as alignments of the same read mapping to different genomic loci with
> the same edit distance) that contain a higher or equal number of indels
> compared to mismatches are retained.

=== "Input"

    Alignments file, sorted, filtered (`.sam`); from
    [**remove_inferiors**](map.md#remove_inferiors)

=== "Output"

    Alignments file, sorted, filtered (`.sam`); used in
    [**convert_all_alns_sam_to_bam**](map.md#convert_all_alns_sam_to_bam) and
    [**filter_sam_by_intersecting_primir**](quantify.md#filter_sam_by_intersecting_primir)

=== "Examples"

    ```console
    Example 1 | Different number of indels

    IN:
        1-1	0	19	77595	255	8M1D14M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:3G1T2^A14	NH:i:2	NM:i:3	XA:Z:Q	XI:i:1
        1-1	0	19	330456	255	4M1D1M1I3M1D13M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:4^G4^A13	NH:i:2	NM:i:3	XA:Z:Q	XI:i:0

    Alignments:
        CTGACATC-AGTGATTCTCCTGC
        ||| | || |||||||||||||| (1 indel, 2 mismatches, discarded)
        CTGGCTTCAAGTGATTCTCCTGC

        CTGA-CATCA-GTGATTCTCCTGC
        |||| | ||| ||||||||||||| (3 indels, 0 mismatches, retained)
        CTGAGC-TCAAGTGATTCTCCTGC

    OUT:
        1-1	0	19	330456	255	4M1D1M1I3M1D13M	*	0	0	CTGACATCAGTGATTCTCCTGC	*	MD:Z:4^G4^A13	NH:i:1	HI:i:1  NM:i:3	XA:Z:Q	XI:i:0


    Example 2 | Equal number of indels

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
[**SAMtools**](map.md#third-party-software-used).

> Required by [BEDTools](map.md#third-party-software-used) to intersect alignments
> with pri-miR annotations.

=== "Input"

    Alignments file, filtered (`.sam`); from
    [**filter_by_indels**](map.md#filter_by_indels)

=== "Output"

    Alignments file (`.bam`); used in
    [**sort_all_alns_bam_by_position**](map.md#sort_all_alns_bam_by_position)


### `sort_all_alns_bam_by_position`

Sort alignments by position with [**SAMtools**](map.md#third-party-software-used).

> Required by [BEDTools](map.md#third-party-software-used) to intersect alignments
> with pri-miR annotations more efficiently.

=== "Input"

    Alignments file (`.bam`); from
    [**convert_all_alns_sam_to_bam**](map.md#convert_all_alns_sam_to_bam)

=== "Output"

    Alignments file, sorted (`.bam`); used in
    [**index_all_alns_bam**](map.md#index_all_alns_bam) and
    [**intersect_extended_primir**](quantify.md#intersect_extended_primir)


### `index_all_alns_bam`

Create index BAM file with [**SAMtools**](map.md#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
> alignments in a genomic region of interest.

=== "Input"

    Alignments file, sorted (`.bam`); from
    [**sort_all_alns_bam_by_position**](map.md#sort_all_alns_bam_by_position)

=== "Output"

    BAM index file (`.bam.bai`); used in
    [**intersect_extended_primir**](quantify.md#intersect_extended_primir)
