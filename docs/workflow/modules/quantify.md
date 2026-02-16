# Quantify Module
<!-- Change? -->

Welcome to the technical documentation for the quantify module. This page aims
to be a detailed documentation of each rule within the module by stating its
inputs, outputs (and how they relate to other rules), configurable parameters,
and the software used. Moreover, when needed, there will be explanations and
examples of what that particular rule does.

The schema below is a visual representation of the individual module steps
and how they are related.

<div align="center">
    <img src=../../../images/quantify.png>
</div>


## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
| **BEDTools** | [GPLv2][license-gpl2] | _"[...] intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][docs-bedtools] / [publication][pub-bedtools] |
| **FASTX-Toolkit** | [AGPL-3.0][license-agpl3] | _"[...] collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing"_ | [code][code-fastx] / [manual][docs-fastx] |
| **SAMtools** | [MIT][license-mit] | _"[...] suite of programs for interacting with high-throughput sequencing data"_ | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools] |


## Configuration file

Some parameters within the workflow can be modified. Refer to the
[configuration template][mirflowz_conf] for a detailed explanation of each
option.

## Quantify Workflow

### `finish_quantify`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

=== "Input"

      - Alignments file, filtered (`.sam`); from
    [**filter_sam_by_intersecting_primir**](quantify.md#filter_sam_by_intersecting_primir)
    and [**filter_sam_by_intersecting_mirna**](quantify.md#filter_sam_by_intersecting_mirna)
    - (iso)miR and/or pri-miR counts table (`.tab`); from
    [**merge_tables**](quantify.md#merge_tables)
    - Alignments file, uncollapsed, sorted (`.bam`); from
    [**sort_uncollapsed_reads_bam_by_position**](quantify.md#sort_uncollapsed_reads_bam_by_position)
    - BAM index file (`.bam.bai`); from
    [**index_uncollapsed_reads_bam**](quantify.md#index_uncollapsed_reads_bam)


### `intersect_extended_primir`

Intersect the aligned reads with the extended pri-miR annotations with
[**BEDTools**](quantify.md#third-party-software-used).

> Only those alignments fully intersecting a (possibly extended) pri-miR
> annotated region are kept.

=== "Input"

    - Alignments file, sorted (`.bam`); from
    [**sort_all_alns_bam_by_position**](map.md#sort_all_alns_bam_by_position)
    - pri-miR extended annotations (`.gff3`); from
    [**extend_mirs_annotations**](prepare.md#extend_mirs_annotations)

=== "Output"

    pri-miR intersections file (`.intersect`); used in
    [**filter_sam_by_intersecting_primir**](quantify.md#filter_sam_by_intersecting_primir)
    and [**quantify_primir**](quantify.md#quantify_primir)


### `filter_sam_by_intersecting_primir`

Remove alignments that do not intersect with any pri-miR with
[**SAMtools**](quantify.md#third-party-software-used).

> Required to only intersect alignments within a (possibly extended) pri-miR
> locus.

=== "Input"

    - Alignments file, filtered (`.sam`); from
    [**filter_by_indels**](map.md#filter_by_indels)
    - pri-miR intersections file (`.intersect`); from
    [**intersect_extended_primir**](quantify.md#intersect_extended_primir)

=== "Output"

    (**Workflow output**) Alignments file, filtered (`.sam`); used in
    [**convert_intersecting_primir_sam_to_bam**](quantify.md#convert_intersecting_primir_sam_to_bam)
    and [**filter_sam_by_intersecting_mirna**](quantify.md#filter_sam_by_intersecting_mirna)


### `convert_intersecting_primir_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](quantify.md#third-party-software-used).

> Required by [BEDTools](quantify.md#third-party-software-used) to intersect alignments
> with miRNA annotations.

=== "Input"

    Alignments file, filtered (`.sam`); from
    [**filter_sam_by_intersecting_primir**](quantify.md#filter_sam_by_intersecting_primir)

=== "Output"

    Alignments file (`.bam`); used in
    [**sort_intersecting_primir_bam_by_position**](quantify.md#sort_intersecting_primir_bam_by_position)


### `sort_intersecting_primir_bam_by_position`

Sort alignments by position with [**SAMtools**](quantify.md#third-party-software-used).

> Required by [BEDTools](quantify.md#third-party-software-used) to intersect alignments
> with miRNA annotations more efficiently.

=== "Input"

    Alignments file (`.bam`); from
    [**convert_intersecting_primir_sam_to_bam**](quantify.md#convert_intersecting_primir_sam_to_bam)

=== "Output"

    Alignments file, sorted (`.bam`); used in
    [**index_intersecting_primir_bam**](quantify.md#index_intersecting_primir_bam) and
    [**intersect_extended_mirna**](quantify.md#intersect_extended_mirna)


### `index_intersecting_primir_bam`

Create index BAM file with [**SAMtools**](quantify.md#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
> alignments in a genomic region of interest.

=== "Input"

    Alignments file, sorted (`.bam`); from
    [**sort_intersecting_primir_bam_by_position**](quantify.md#sort_intersecting_primir_bam_by_position)

=== "Output"

    BAM index file (`.bam.bai`); used in
    [**intersect_extended_mirna**](quantify.md#intersect_extended_mirna)


### `intersect_extended_mirna`

Intersect the aligned reads with the extended miRNA annotations with
[**BEDTools**](quantify.md#third-party-software-used).

> Only those alignments fully intersecting an extended mature miRNA annotated
> region are kept.

=== "Input"

    - Alignments file, sorted (`.bam`); from
    [**sort_intersecting_primir_bam_by_position**](quantify.md#sort_intersecting_primir_bam_by_position)
    - Mature miRNA extended annotations (`.gff3`); from
    [**extend_mirs_annotations**](prepare.md#extend_mirs_annotations)

=== "Output"

    Mature miRNA intersections file (`.intersect`); used in
    [**filter_sam_by_intersecting_mirna**](quantify.md#filter_sam_by_intersecting_mirna)
    and [**add_intersecting_mirna_tag**](quantify.md#add_intersecting_mirna_tag)


### `filter_sam_by_intersecting_mirna`

Remove alignments that do not intersect with any miRNA with
[**SAMtools**](quantify.md#third-party-software-used).

> Required to efficiently classify the alignments.

=== "Input"

    - Alignments file, filtered (`.sam`); from
    [**filter_sam_by_intersecting_primir**](quantify.md#filter_sam_by_intersecting_primir)
    - Mature miRNA intersections file (`.intersect`); from
    [**intersect_extended_mirna**](quantify.md#intersect_extended_mirna)

=== "Output"

    (**Workflow output**) Alignments file, filtered (`.sam`); used in
    [**add_intersecting_mirna_tag**](quantify.md#add_intersecting_mirna_tag) and
    [**uncollapse_reads**](quantify.md#uncollapse_reads)


### `add_intersecting_mirna_tag`

Classify and add the intersecting (iso)miR to each alignment as a tag
with a [**custom script**][custom-script-iso-tag].

> In this step, the mature miRNA annotated regions are used instead of the
> extended ones. Each alignment gets an extra tag (`YW:Z`) with the (iso)miR(s)
> it is considered to really intersect with using the format:
> `miRNA_name|5p-shift|3p-shift|CIGAR|MD|READ_SEQ`, where `5p-shift` and
> `3p-shift` are the difference between the miRNA start and end coordinates
> and the alignment's ones respectively.

=== "Input"

    - Alignments file, filtered (`.sam`); from
    [**filter_sam_by_intersecting_mirna**](quantify.md#filter_sam_by_intersecting_mirna)
    - Mature miRNA intersections file (`.intersect`); from
    [**intersect_extended_mirna**](quantify.md#intersect_extended_mirna)

=== "Parameters"

      - **config_template.yaml**
        - `extension`: Number of nucleotides by which mature miRNA annotated
        regions are extended (default: 6)

=== "Output"

    Alignments file, tagged (`.sam`); used in
    [**sort_intersecting_mirna_by_feat_tag**](quantify.md#sort_intersecting_mirna_by_feat_tag)

=== "Examples"

    ```console
    Example 1 | Feature intersects alignment | coordinates adjustment and shift allowed
        use case:
            Prior to checking if the feature is intersecting the alignment, its
            coordinates are adjusted by the value specified in `--extension`. In
            addition, the same value is used to specify the +/- shift allowed
            between the feature and the read alignment start and end coordinates.

        command:
            annotate_sam_with_intersecting_features.py -b INTERSECT -s SAM --extension 5

        in INTERSECT record:
            19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_1	255	+	21

        in SAM record:
            read_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

        intersection before coordinates adjustment visualization:

            ---|===============================|--- (feature)
            ---------|===================|--------- (read)

        intersection after coordinates adjustment visualization:

            --------|=====================|-------- (feature)
            ---------|===================|--------- (read)

        out SAM record:
            read_1	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:hsa-miR-1323|1|-1|21M|21|TCAAAACTGAGGGGCATTTTC

        description:
            The feature start and end coordinates after the adjustment are 5337
            and 5360 respectively.
            The read alignment starts at position 5338. As the read has length 21,
            its end position is 5359.
            The feature intersects the read alignment with an overhang within the
            specified shift range (+/- 5) so it is added as a new tag in the output
            SAM record.


    Example 2 | Feature intersects alignment | no coordinates adjustment or shift allowed
        use case:
            The feature and read alignment coordinates must perfectly match.

        command:
            annotate_sam_with_intersecting_features.py -b INTERSECT -s SAM

        in INTERSECT record:
            19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_2	255	+	21

        in SAM record:
            read_2	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

        intersection visualization:

            ---------|===================|--------- (feature)
            ---------|===================|--------- (read)

        out SAM record:
            read_2	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

        description:
            The feature start and end coordinates are 5338 and 5359 respectively.
            The read alignment starts at position 5338. As the read has length 21,
            its end position is 5359.
            The feature perfectly intersects the read alignment so it is added as
            a new tag in the output SAM record.


    Example 3 | Non-intersecting feature | shift filter not passed
        use case:
            Prior to checking if the feature is intersecting the alignment, its
            coordinates are adjusted by the value specified in `--extension`. In
            addition, the same value is used to specify the +/- shift allowed
            between the feature and the read alignment start and end coordinates.

        command:
            annotate_sam_with_intersecting_features.py -b INTERSECT -s SAM --extension 1

        in INTERSECT record:
            19	.	miRNA	5332	5365	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_3	255	+	21

        in SAM record:
            read_3	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

        intersection before coordinates adjustment visualization:

            ---|===============================|--- (feature)
            ---------|===================|--------- (read)

        intersection after coordinates adjustment visualization:

            ----|=============================|---- (feature)
            ---------|===================|--------- (read)

        out SAM record:
            read_3	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:

        description:
            The feature start and end coordinates after the adjustment are 5333
            and 5364 respectively.
            The read alignment starts at position 5338. As the read has length 21,
            its end position is 5359.
            There is a 5-nucleotide overhang on both ends. Thus, the feature is
            not considered to intersect the read alignment and the alignment is
            not written in the output file.


    Example 4 | Feature intersects alignment | using feature's "Alias"
        use case:
            The feature and read alignment coordinates must perfectly match.

        command:
            annotate_sam_with_intersecting_features.py -b INTERSECT -s SAM --id alias

        in INTERSECT record:
            19	.	miRNA	5338	5359	.	+	.	ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786	19	5337	5358	read_4	255	+	21

        in SAM record:
            read_4	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0

        intersection visualization:

            ---------|===================|--------- (feature)
            ---------|===================|--------- (read)

        out SAM record:
            read_4	0	19	5338	255	21M	*	0	0	TCAAAACTGAGGGGCATTTTC	*	MD:Z:21	NH:i:1	NM:i:0  YW:Z:MIMAT0005795|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

        description:
            The feature start and end coordinates are 5338 and 5359 respectively.
            The read alignment starts at position 5338. As the read has length 21,
            its end position is 5359.
            The feature perfectly intersects the read alignment so it is added as
            a new tag in the output SAM record. In this case, instead of using the
            feature `Name` (default), the `Alias` is used.
    ```


### `sort_intersecting_mirna_by_feat_tag`

Sort the alignments by the tag containing the classified intersecting (iso)miR
with [**SAMtools**](quantify.md#third-party-software-used).

> Required for an efficient quantification.

=== "Input"

    Alignments file, tagged (`.sam`); from
    [**add_intersecting_mirna_tag**](quantify.md#add_intersecting_mirna_tag)

=== "Output"

    Alignments file, tagged, sorted (`.sam`); used in
    [**quantify_mirna**](quantify.md#quantify_mirna)


### `quantify_mirna`

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
> are no mismatches nor indels.

=== "Input"

    Alignments file, tagged, sorted (`.sam`); from
    [**sort_intersecting_mirna_by_feat_tag**](quantify.md#sort_intersecting_mirna_by_feat_tag)

=== "Parameters"

    - **samples.tsv**
        - Library name; specified in the sample's table column `sample`
    - **config_template.yaml**
        - `mir_list`: miRNA features to be quantified (default: isomir, mirna
        pri-miR)

=== "Output"

    (iso)miR counts tab-delimited file; used in
    [**merge_tables**](quantify.md#merge_tables)

=== "Examples"

    ```console
    Example 1 | Canonical miRNA and isomiR

    IN SAM record:
        10-4_2	0	19	34627	255	21M	*	0	0	AAAGTGCTTCCTTTTAGAGGG	*	MD:Z:21	NM:i:0	NH:i:2	HI:i:1	YW:Z:hsa-miR-520b-3p|0|0|21M|21|AAAGTGCTTCCTTTTAGAGGG
        10-4_2	0	19	40866	255	21M	*	0	0	AAAGTGCTTCCTTTTAGAGGG	*	MD:Z:21	NM:i:0	NH:i:2	HI:i:2	YW:Z:hsa-miR-520c-3p|0|-1|21M|21|AAAGTGCTTCCTTTTAGAGGG

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
        599-1_3	0	19	27804	255	20M	*	0	0	AAAGTGCTTCCTTTTAGAGG	*	MD:Z:20	NM:i:0	NH:i:3	HI:i:1	YW:Z:hsa-miR-526b-3p|1|-1|20M|20|AAAGTGCTTCCTTTTAGAGG
        599-1_3	0	19	34627	255	20M	*	0	0	AAAGTGCTTCCTTTTAGAGG	*	MD:Z:20	NM:i:0	NH:i:3	HI:i:2	YW:Z:hsa-miR-520b-3p|0|-1|20M|20|AAAGTGCTTCCTTTTAGAGG
        599-1_3	0	19	40866	255	20M	*	0	0	AAAGTGCTTCCTTTTAGAGG	*	MD:Z:20	NM:i:0	NH:i:3	HI:i:3	YW:Z:hsa-miR-520c-3p|0|-2|20M|20|AAAGTGCTTCCTTTTAGAGG

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
        hsa-miR-520c-3p|0|-2|20M|20|AAAGTGCTTCCTTTTAGAGG	0.33
        hsa-miR-526b-3p|1|-1|20M|20|AAAGTGCTTCCTTTTAGAGG	0.33
    ```


### `quantify_primir`

Tabulate alignments according to its intersecting pri-miR with a
[**custom script**][custom-script-pri-quant]

> Each alignment contributes to the pri-miR it intersects with by `R/N`, where
> `R` is the number of collapsed reads in that alignment, and `N` is the
> number of genomic and/or transcriptomic loci it aligns to. The values of
> both, `R` and `N` are inferred from the sequence name which follows the
> format `ID-R_N`. The resulting table has a row for each pri-miR with the
> name format set in
> [**extend_mirs_annotations**](prepare.md#extend_mirs_annotations).

=== "Input"

    pri-miR intersections file (`.intersect`); from
    [**intersect_extended_primir**](quantify.md#intersect_extended_primir)

=== "Output"

    pri-miR counts tab-delimited file; used in
    [**merge_tables**](quantify.md#merge_tables)

=== "Examples"

    ```console
    Example 1 | Contribution when using '--collapsed'
        use case:
            A single feature with several intersecting reads.
            The flag '--collapsed' is used, so contribution equals the # of reads
            per alignment.

        IN INTERSECT records:
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2	255	+	21
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1	255	+	21

        alignments:
            Read ID: 8-2
            Number of collapsed reads: 2
            Contribution: 2

            Read ID: 24-1
            Number of collapsed reads: 1
            Contribution: 1

        OUT table:
            hsa-mir-524_-0_+0      3


    Example 2 | Contribution when using '--nh'
        use case:
            A single feature with several intersecting reads.
            The flag '--nh' is used, so contribution equals 1/NH.

        IN INTERSECT records:
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8_1	255	+	21
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24_1	255	+	21

        alignments:
            Read ID: 8_1
            Number of mapped genomic loci: 1
            Contribution: 1

            Read ID: 24_1
            Number of mapped genomic loci: 1
            Contribution: 1

        OUT table:
            hsa-mir-524_-0_+0      2


    Example 3 | Contribution when using '--collapsed' and '--nh'
        use case:
            A single feature with several intersecting reads.
            The flags '--nh' and '--contribution' are used, so contribution equals
            # of reads/NH.

        IN INTERSECT records:
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2_1	255	+	21
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1_1	255	+	21

        alignments:
            Read ID: 8-2_1
            Number of collapsed reads: 2
            Number of mapped genomic loci: 1
            Contribution: 2/1 = 2

            Read ID: 23-1_1
            Number of collapsed reads: 1
            Number of mapped genomic loci: 1
            Contribution: 1/1 = 1

        OUT table:
            hsa-mir-524_-0_+0      3

    Example 4 | Column with intersecting reads; using '--read-ids'
        use case:
            A single feature with several intersecting reads.
            Each read contributes with 1.
            Read IDs intersecting the feature are appended as the last column.

        IN INTERSECT records:
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	8-2_1	255	+	21
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0	19	44413	44434	24-1_1	255	+	21

        alignments:
            Read ID: 8-2_1
            Contribution: 1

            Read ID: 24-1_1
            Contribution: 1

        OUT table:
            hsa-mir-524_-0_+0      2       8-2_1;24-1_1

    Example 5 | Columns with feature shifts; using '--feat-extension'
        use case:
            A single feature with several intersecting reads.
            Each read contributes with 1.
            Feature start and end coordinates shift are appended as the third and
            fourth column respectively.

        IN INTERSECT records:
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-1_+3	19	44413	44434	8-2_1	255	+	21
            19	.	miRNA_primary_transcript	44362	44448	.	+	.	ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-1_+3	19	44413	44434	24-1_1	255	+	21

        alignments:
            Read ID: 8-2_1
            Contribution: 1

            Read ID: 24-1_1
            Contribution: 1

        feature:
            Feature name: hsa-mir-524_-1_+3
            5' shift: -1
            3' shift: +3

        OUT table:
            hsa-mir-524_-1_+3      2       -1       +3
    ```


### `merge_tables`

Merge all the tables from the different libraries into a single one with a
[**custom script**][custom-script-merge-tab].

> The final table(s) containing the counting data from all libraries for the
> (iso)miRs and/or pri-miRs have a row per miRNA species and a column per
> sample library. If a miRNA species is not found in a certain library, its
> value is set to `NA`.

=== "Input"

    Counts tab-delimited file; from [**quantify_mirna**](quantify.md#quantify_mirna)
    and/or [**quantify_primir**](quantify.md#quantify_primir)

=== "Parameters"

    - **cluster_schema.json**
        - `mir_list`: miRNA features to be quantified (default: isomir, mirna
        pri-mir)

=== "Output"

    (**Workflow output**) (iso)miR and/or pri-miR counts table (`.tab`)

=== "Example"

    ```console
    IN library 1
        ID                              lib_1
        hsa-miR-524-5p          	    1
        hsa-miR-524-5p|0|0|22M|9G12	    1
        hsa-miR-524-5p|0|0|22M|9G9C2	1

    IN library 2
        ID                              lib_2
        hsa-miR-524-5p          	    1
        hsa-miR-1283	                1
        hsa-miR-1283|-1|-2|21M|21	    1

    IN library 3
        ID                              lib_3

    OUT table
        ID                              lib_1  lib_2  lib_3
        hsa-miR-524-5p          	    1	    1       NA
        hsa-miR-524-5p|0|0|22M|9G12	    1   	NA      NA
        hsa-miR-524-5p|0|0|22M|9G9C2	1   	NA      NA
        hsa-miR-1283	                NA	    1       NA
        hsa-miR-1283|-1|-2|21M|21	    NA	    1       NA
    ```


### `uncollapse_reads`

Reverse the collapsing of reads with identical sequences as done with
[**FASTX-Toolkit**](quantify.md#third-party-software-used) with a
[**custom script**][custom-script-uncollapse].

=== "Input"

    Alignments file, filtered (`.sam`); from
    [**filter_sam_by_intersecting_mirna**](quantify.md#filter_sam_by_intersecting_mirna)

=== "Output"

    Uncollapsed aligned reads (`.sam`); used in
    [**convert_uncollapsed_reads_sam_to_bam**](quantify.md#convert_uncollapsed_reads_sam_to_bam)


### `convert_uncollapsed_reads_sam_to_bam`

Convert alignments `.sam` file to `.bam` with
[**SAMtools**](quantify.md#third-party-software-used).

=== "Input"

    Alignments file, uncollapsed (`.sam`); from
    [**filter_sam_by_intersecting_mirna**](quantify.md#filter_sam_by_intersecting_mirna)

=== "Output"

    Alignments file, uncollapsed (`.bam`); used in
    [**sort_uncollapsed_reads_bam_by_position**](quantify.md#sort_uncollapsed_reads_bam_by_position)


### `sort_uncollapsed_reads_bam_by_position`

Sort alignments by position with [**SAMtools**](quantify.md#third-party-software-used).

=== "Input"

    Alignments file, uncollapsed (`.bam`); from
    [**convert_uncollapsed_reads_sam_to_bam**](quantify.md#convert_uncollapsed_reads_sam_to_bam)

=== "Output"

    (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); used in
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups) and/or
    [**create_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups)


### `index_uncollapsed_reads_bam`

Create index BAM file with [**SAMtools**](quantify.md#third-party-software-used).

> Indexing is required by genome viewers such as IGV to quickly display
> alignments in a genomic region of interest.

=== "Input"

    (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
    [**sort_uncollapsed_reads_bam_by_position**](quantify.md#sort_uncollapsed_reads_bam_by_position)

=== "Output"

    (**Workflow output**) BAM index file (`.bam.bai`); used in
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups) and/or
    [**create_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups)
