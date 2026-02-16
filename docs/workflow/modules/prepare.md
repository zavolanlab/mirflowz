# Prepare Module
<!-- Change? -->

Welcome to the technical documentation for the prepare module. This page aims
to be a detailed documentation of each rule within the module by stating its
inputs, outputs (and how they relate to other rules), configurable parameters,
and the software used. Moreover, when needed, there will be explanations and
examples of what that particular rule does.

The schema below is a visual representation of the individual module steps
and how they are related.

<div align="center">
    <img src=../../../images/prepare.png>
</div>


## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

|      Name     |         License         |                                                               Tag Line                                                               |                                     More Info                                    |
|:-------------:|:-----------------------:|:------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------:|
| **cufflinks** | [BSL-1.0][license-bsl1] | _"[...] assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples"_ | [code][code-cufflinks] / [manual][docs-cufflinks] / [publication][pub-cufflinks] |
| **GFFUtils**  | [AFL-3][license-afl3]   | _"[...] a small set of utility programs for working with GFF and GTF files"_                                                         | [code][code-gffutils] / [manual][docs-gffutils]                                  |
| **SAMtools**  | [MIT][license-mit]      | _"[...] suite of programs for interacting with high-throughput sequencing data"_                                                     | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools]    |
| **segemehl**  | [GPLv3][license-gpl3]   | _"[...] map short sequencer reads to reference genomes"_                                                                             | [manual][docs-segemehl] / [publication][pub-segemehl]                            |


## Configuration file

Some parameters within the workflow can be modified. Refer to the
[configuration template][mirflowz_conf] for a detailed explanation of each
option.

## Prepare Workflow

### `finish_prepare`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

=== "Input"

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


### `trim_genome_seq_ids`

Trim genome sequence IDs with a [**custom script**][custom-script-trim-id].

=== "Input"

    (**Workflow input**) Genome sequence, `gzip`ed (`.fa.gz`/`.fasta.gz`)

=== "Output"

    Genome sequence, trimmed IDs (`.fa`); used in
    [`extract_transcriptome_seqs`](prepare.md#extract_transcriptome_seqs),
    [`create_genome_header`](prepare.md#create_genome_header),
    [`create_index_genome_fasta`](prepare.md#create_index_genome_fasta),
    [`generate_segemehl_index_genome`](prepare.md#generate_segemehl_index_genome),
    [`map_genome_segemehl`](map.md#map_genome_segemehl),
    [`map_genome_oligomap`](map.md#map_genome_oligomap) and
    [`compress_reference_genome`](pileups.md#compress_reference_genome)

### `extract_transcriptome_seqs`

Create transcriptome from genomic sequence and annotations with
[**cufflinks**](prepare.md#third-party-software-used).

=== "Input"

    - (**Workflow input**) Genome annotations, `gzip`ed (`.gtf.gz`)
    - Genome sequence, trimmed IDs (`.fa`); from
      [`trim_genome_seq_ids`](prepare.md#trim_genome_seq_ids)


=== "Output"

    Transcriptome sequence (`.fa`); used in
    [`trim_transcriptome_seq_ids`](prepare.md#trim_transcriptome_seq_ids)


### `trim_transcriptome_seq_ids`

Trim transcriptome sequence IDs with a
[**custom script**][custom-script-trim-id].

=== "Input"

    Transcriptome sequence (`.fa`); from
    [`extract_transcriptome_seqs`](prepare.md#extract_transcriptome_seqs)

=== "Output"

    Transcriptome sequence, trimmed IDs (`.fa`); used in
    [`generate_segemehl_index_transcriptome`](prepare.md#generate_segemehl_index_transcriptome),
    [`map_transcriptome_oligomap`](map.md#map_transcriptome_oligomap) and
    [`map_transcriptome_segemehl`](map.md#map_transcriptome_segemehl)

### `generate_segemehl_index_transcriptome`

Generate transcriptome index for [**segemehl**](prepare.md#third-party-software-used)
short read aligner.

> The transcriptome index only needs to be generated once for each combination
> of transcriptome sequence and annotations, and sample sets.

=== "Input"

    Transcriptome sequence, trimmed IDs (`.fa`); from
    [**trim_transcriptome_seq_ids**](prepare.md#trim_transcriptome_seq_ids)

=== "Output"

    segemehl transcriptome index (`.idx`); used in
    [**map_transcriptome_segemehl**](map.md#map_transcriptome_segemehl)

### `generate_segemehl_index_genome`

Generate genome index for [**segemehl**](prepare.md#third-party-software-used) short read
aligner.

> The genome index only needs to be generated once for each combination
> of annotations, and sample sets.

=== "Input"

    Genome sequence, trimmed IDs (`.fa`); from
    [**trim_genome_seq_ids**](prepare.md#trim_genome_seq_ids)

=== "Output"

    segemehl genome index (`.idx`); used in
    [**map_genome_segemehl**](map.md#map_genome_segemehl)

### `get_exons_gtf`

Retrieve exon annotations from genome annotations with a
[**custom script**][custom-script-get-lines].

=== "Input"

    (**Workflow input**) Genome annotations, `gzip`ed (`.gtf.gz`)

=== "Output"

    Exon annotations (`.gtf`); used in
    [**convert_exons_gtf_to_bed**](prepare.md#convert_exons_gtf_to_bed)


### `convert_exons_gtf_to_bed`

Convert exon annotations `.gtf` to `.bed` with a
[**custom script**][custom-script-gtf-bed].

=== "Input"

    Exon annotations (`.gtf`); from
    [**get_exons_gtf**](prepare.md#get_exons_gtf)

=== "Output"

    Exon annotations (`.bed`); used in
    [**transcriptome_to_genome_maps**](map.md#transcriptome_to_genome_maps)


### `create_genome_header`

Create SAM header for the genome with
[**SAMtools**](prepare.md#third-party-software-used).

> Required by [SAMtools](prepare.md#third-party-software-used) to work with the alignment
> files.

=== "Input"

    Genome sequence, trimmed IDs (`.fa`); from
    [**trim_genome_seq_ids**](prepare.md#trim_genome_seq_ids)

=== "Output"

    SAM genome header (`.sam`); used in
    [**add_header_all_maps**](map.md#add_header_all_maps)



### `map_chr_names`

Map UCSC-like chromosome names with Ensembl-like ones in miRNA annotations
with a [**custom script**][custom-script-map-chr].

> Required by [BEDTools](prepare.md#third-party-software-used) to intersect alignments with
> miRNA annotations. Several mapping tables are available [here][chr_map].

=== "Input"

    - (**Workflow input**) miRNA annotations (`.gff3`)
    - (**Workflow input**) Tab-separated chromosome name mappings table
      (`.tsv`)

=== "Output"

    miRNA annotations, mapped chromosome name(s) (`.gff3`); used in
    [**extend_mirs_annotations**](prepare.md#extend_mirs_annotations),
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups) and/or
    [**create_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups)


### `create_index_genome_fasta`

Create a FASTA index for the genome with
[**SAMtools**](prepare.md#third-party-software-used).

=== "Input"

    Genome sequence, trimmed IDs (`.fa`); from
    [**trim_genome_seq_ids**](prepare.md#trim_genome_seq_ids)

=== "Output"

    FASTA genome index (`.fa.fai`); used in
    [**extract_chr_len**](prepare.md#extract_chr_len)


### `extract_chr_len`

Extract chromosome(s) length from the genome sequence.

> Required to ensure that the extended annotations in generated in the
> [**extend_mirs_annotations**](prepare.md#extend_mirs_annotations) rule do not
> exceed the chromosome(s) boundaries.


=== "Input"

    FASTA genome index (`.fa.fai`); from
    [**create_index_genome_fasta**](prepare.md#create_index_genome_fasta)

=== "Output"

    Tab-separated table mapping chromosome name(s) and length(s) (`.tsv`); used
    in [**extend_mirs_annotations**](prepare.md#extend_mirs_annotations)

### `extend_mirs_annotations`

Extend miRNA annotations, ensure feature names uniqueness and split the file
by feature type with a [**custom script**][custom-script-mir-ext].

=== "Input"

    miRNA annotations, mapped chromosome name(s) (`.gff3`); from
    [**map_chr_names**](prepare.md#map_chr_names)

=== "Parameters"

    - **config_template.yaml**
        - `extension`: Number of nucleotides by which mature miRNA annotated
           regions are extended at most (default: 6)

=== "Output"

    - Primary miRNA transcript (pri-miR) extended annotations (`.gff3`); used
    in [**intersect_extended_primir**](quantify.md#intersect_extended_primir)
    - Mature miRNA (miRNA) extended annotations (`.gff3`); used in
    [**intersect_extended_mirna**](quantify.md#intersect_extended_mirna)


=== "Extension Examples"

    ```console
    Example 1 | Mature miRNA extension

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


    Example 2 | Mature miRNA and pri-miR extension

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


    Example 3 | Matrue miRNA exceeding chromosome boundaries extension

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
    ```

=== "Name Uniqueness Examples"

    ```console
    Example 4 | Replica number in the ID; 'Derives_from' update

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
            chr21	.	miRNA	8250772	8250791	.	+	.	ID=MIMAT0041633_1;Alias=MIMAT0041633;Name=hsa-miR-10401-2-5p;Derives_from=MI0033425_2
            chr21	.	miRNA	8250807	8250827	.	+	.	ID=MIMAT0041634_1;Alias=MIMAT0041634;Name=hsa-miR-10401-2-3p;Derives_from=MI0033425_2


    Example 5 | Replica number in the Name; single mature arm

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


    Example 6 | Both mature miRNA arms but just one with the replica number

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


    Example 7 | Different precursor and mature miRNA "NAME" in Name

    IN:
        pri-miR entries:
            chr19	.	miRNA_primary_transcript	45628	45714	.	+	.	ID=MI0003161;Alias=MI0003161;Name=hsa-mir-517a
            chr19	.	miRNA_primary_transcript	54436	54502	.	+	.	ID=MI0003165;Alias=MI0003165;Name=hsa-mir-517b
        mature miRNA entries:
            chr19	.	miRNA	45642	45663	.	+	.	ID=MIMAT0002851;Alias=MIMAT0002851;Name=hsa-miR-517-5p;Derives_from=MI0003161
            chr19	.	miRNA	45681	45702	.	+	.	ID=MIMAT0002852;Alias=MIMAT0002852;Name=hsa-miR-517a-3p;Derives_from=MI0003161
            chr19	.	miRNA	54441	54462	.	+	.	ID=MIMAT0002851_1;Alias=MIMAT0002851;Name=hsa-miR-517-5p;Derives_from=MI0003165
            chr19	.	miRNA	54478	54499	.	+	.	ID=MIMAT0002857;Alias=MIMAT0002857;Name=hsa-miR-517b-3p;Derives_from=MI0003165
    OUT:
        pri-miR entries:
            chr19	.	miRNA_primary_transcript	45628	45714	.	+	.	ID=MI0003161;Alias=MI0003161;Name=hsa-mir-517a
            chr19	.	miRNA_primary_transcript	54436	54502	.	+	.	ID=MI0003165;Alias=MI0003165;Name=hsa-mir-517b
        mature miRNA entries:
            chr19	.	miRNA	45642	45663	.	+	.	ID=MIMAT0002851;Alias=MIMAT0002851;Name=hsa-miR-517a-5p;Derives_from=MI0003161
            chr19	.	miRNA	45681	45702	.	+	.	ID=MIMAT0002852;Alias=MIMAT0002852;Name=hsa-miR-517a-3p;Derives_from=MI0003161
            chr19	.	miRNA	54441	54462	.	+	.	ID=MIMAT0002851_1;Alias=MIMAT0002851;Name=hsa-miR-517b-5p;Derives_from=MI0003165
            chr19	.	miRNA	54478	54499	.	+	.	ID=MIMAT0002857;Alias=MIMAT0002857;Name=hsa-miR-517b-3p;Derives_from=MI0003165
    ```
