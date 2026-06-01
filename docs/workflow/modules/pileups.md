# Pileup Module
<!-- Change? -->

Welcome to the technical documentation for the pileup module. This page aims
to be a detailed documentation of each rule within the module by stating its
inputs, outputs (and how they relate to other rules), configurable parameters,
and the software used. Moreover, when needed, there will be explanations and
examples of what that particular rule does.

The schema below is a visual representation of the individual module steps
and how they are related.

<div align="center">
    <img src=../../../images/pileups.png>
</div>

## Third-party software used

> Tag lines were taken from the developers' websites (code repository or manual)

| Name | License | Tag line | More info |
| --- | --- | --- | --- |
| **ASCII-style alignment pileups** | [Apache 2.0][license-apache2] | _"Generates ASCII-style pileups of read alignments in one or more BAM files for one or more genomic regions."_ | [code][code-ascii] |
| **SAMtools** | [MIT][license-mit] | _"[...] suite of programs for interacting with high-throughput sequencing data"_ | [code][code-samtools] / [manual][docs-samtools] / [publication][pub-samtools] |


## Configuration file

Some parameters within the workflow can be modified. Refer to the
[configuration template][mirflowz_conf] for a detailed explanation of each
option.

## Pileup Workflow

### `finish_pileup`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

=== "Input"

    (**Workflow output**) Empty text file (`.txt`)
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups),
    [**modify_per_library_ascii_pileups**](pileups.md#modify_per_library_ascii_pileups),
    [**modify_per_run_ascii_pileups**](pileups.md#modify_per_run_ascii_pileups),
    and [**color_code_ascii_pileups**](pileups.md#color_code_ascii_pileups)


### `create_empty_bed`

Create an empty BED file if the user has not provided one.

> **OPTIONAL RULE.** This rule will be executed if, and only if, the user has
> not provided a BED file in the [configuration file][mirflowz_conf]
> with the regions the ASCII-style alignment pileups must be performed on.

=== "Condition"

    - **config_template.yaml**
        - `bed_file`: BED6 file with all the desired annotation regions to perform
        the ASCII-style alignment pileups on. (Default: None)

=== "Output"

    Empty BED file (`.bed`); used in
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups) and/or
    [**create_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups)


### `compress_reference_genome`

Compress the processed genome with trimmed IDs using `bgzip` with
[**SAMtools**](pileups.md#third-party-software-used).

> Required to perform the ASCII-style alignment pileups.

=== "Input"

    Genome sequence, trimmed IDs (`.fa`); from
    [**trim_genome_seq_ids**](prepare.md#trim_genome_seq_ids)

=== "Output"

    Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); used in
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups) and/or
    [**create_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups)


### `create_per_library_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions across
libraries with [**ASCII-style alignment pileups**](pileups.md#third-party-software-used).

> A directory containing the ASCII-style pileups is created for each
> library. If no BED file is provided, the pileups' output directories will
> only contain an empty file.

=== "Input"

    - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); from
    [**compress_reference_genome**](pileups.md#compress_reference_genome)
    - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
    [**map_chr_names**](prepare.md#map_chr_names)
    - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
    [**sort_uncollapsed_reads_bam_by_position**](quantify.md#sort_uncollapsed_reads_bam_by_position)
    - (**Workflow output**) BAM index file (`.bam.bai`); used in
    [**index_uncollapsed_reads_bam**](quantify.md#index_uncollapsed_reads_bam)
    - Annotated genomic regions (`.bed`); from workflow input files or
    [**create_empty_bed**](pileups.md#create_empty_bed)

=== "Parameters"

    - **config_template.yaml**
        - `sort_by`: ASCII-style alignment pileups can be sorted by the first
        nucleotide's position from left-to-right ('position') or by counts in
        descending order ('counts') (default: 'position')

=== "Output"

    (**Workflow output**) Empty text file (`.txt`)


### `create_per_run_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions for the whole
run with [**ASCII-style alignment pileups**](pileups.md#third-party-software-used).

> If no BED file is provided, the pileups' output directory will only contain
> an empty file.

=== "Input"

    - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); from
    [**compress_reference_genome**](pileups.md#compress_reference_genome)
    - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
    [**map_chr_names**](prepare.md#map_chr_names)
    - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
    [**sort_uncollapsed_reads_bam_by_position**](quantify.md#sort_uncollapsed_reads_bam_by_position)
    - (**Workflow output**) BAM index file (`.bam.bai`); used in
    [**index_uncollapsed_reads_bam**](quantify.md#index_uncollapsed_reads_bam)
    - Annotated genomic regions (`.bed`); from workflow input files or
    [**create_empty_bed**](pileups.md#create_empty_bed)

=== "Parameters"

    - **config_template.yaml**
        - `sort_by`: ASCII-style alignment pileups can be sorted by the first
        nucleotide's position from left-to-right ('position') or by counts in
        descending order ('counts') (default: 'position')

=== "Output"

    (**Workflow output**) Empty text file (`.txt`)


### `create_per_condition_ascii_pileups`

Create ASCII-style pileups for all the desired annotated regions across the
different library subsets if provided with
[**ASCII-style alignment pileups**](pileups.md#third-party-software-used).

> **OPTIONAL RULE.** The ASCII-style pileups for each annotated region are
> made if, and only if, at least one library subset is specified in the
> [configuration file](#configuration-file). Otherwise, this rule will not be
> executed, and no output will be generated.

=== "Condition"

    - **config_template.yaml**
        - `lib_dict`: Dictionary of arbitrary condition names (keys) and library
        names to aggregate alignment pileups for (values; MUST correspond to names
        in samples table) (default: None)

=== "Input"

    - Genome sequence, trimmed IDs, `bgzip`ed (`.fa.bz`); from
    [**compress_reference_genome**](pileups.md#compress_reference_genome)
    - miRNA annotations, mapped chromosome name(s) (`.gff3`); from
    [**map_chr_names**](prepare.md#map_chr_names)
    - (**Workflow output**) Alignments file, uncollapsed, sorted (`.bam`); from
    [**sort_uncollapsed_reads_bam_by_position**](quantify.md#sort_uncollapsed_reads_bam_by_position)
    - (**Workflow output**) BAM index file (`.bam.bai`); used in
    [**index_uncollapsed_reads_bam**](quantify.md#index_uncollapsed_reads_bam)
    - Annotated genomic regions (`.bed`); from workflow input files or
    [**create_empty_bed**](pileups.md#create_empty_bed)

=== "Parameters"

    - **config_template.yaml**
        - `lib_dict`: Dictionary of arbitrary condition names (keys) and library
        names to aggregate alignment pileups for (values; MUST correspond to names
        in samples table) (default: None)
        - `sort_by`: ASCII-style alignment pileups can be sorted by the first
        nucleotide's position from left-to-right ('position') or by counts in
        descending order ('counts') (default: 'position')

=== "Output"

    Empty text file (`.txt`)


### `modify_per_library_ascii_pileups`

Modify the generated ASCII-style pileups for all the desired annotated regions
across libraries with a [**custom script**][custom-script-ascii-mod].

> A directory containing the modified ASCII-style pileups is created for each
> library. If no ASCII-style alignment pileups were created, the modified
> pileups' output directories will only contain an empty file.

=== "Input"

    (**Workflow output**) Empty text file (`.txt`); from
    [**create_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups)

=== "Parameters"

    - **config_template.yaml**
        - `sort_by`: ASCII-style alignment pileups can be sorted by the first
          nucleotide's position from left-to-right ('position') or by counts
          in descending order ('counts') (default: 'position')
        - `split`: Split precursor pileups into one mature-arm pileup per arm
          (default: 'true')
        - `canonical`: Mark the aligned read corresponding to the canonical
          sequence (default: 'true')
        - `keep_all`: Write the pileup even if it has no aligned sequences
          (default: 'true')
        - `min_count_dict`: Dictionary with the minimum count for a sequence to
          be kept for each pileup group (default: 'min_count_dict["lib"] = 1')
        - `max_seq`: Maximum number of top sequences to display (default: 15)
        - `extension`: Extension of the mature miRNA start and end coordinates
          in bp (default: 6)

=== "Output"

    (**Workflow output**) Empty text file (`.txt`)


### `modify_per_run_ascii_pileups`

Modify the generated ASCII-style pileups for all the desired annotated regions
for the whole run with a [**custom script**][custom-script-ascii-mod].

> A directory containing the modified ASCII-style pileups is created for the
> whole run. If no ASCII-style alignment pileups were created, the modified
> pileups' output directories will only contain an empty file.

=== "Input"

    (**Workflow output**) Empty text file (`.txt`); from
    [**create_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups)

=== "Parameters"

    - **config_template.yaml**
        - `sort_by`: ASCII-style alignment pileups can be sorted by the first
          nucleotide's position from left-to-right ('position') or by counts
          in descending order ('counts') (default: 'position')
        - `split`: Split precursor pileups into one mature-arm pileup per arm
          (default: 'true')
        - `canonical`: Mark the aligned read corresponding to the canonical
          sequence (default: 'true')
        - `keep_all`: Write the pileup even if it has no aligned sequences
          (default: 'true')
        - `min_count_dict`: Dictionary with the minimum count for a sequence to
          be kept for each pileup group (default: 'min_count_dict["run"] = 1')
        - `max_seq`: Maximum number of top sequences to display (default: 15)
        - `extension`: Extension of the mature miRNA start and end coordinates
          in bp (default: 6)

=== "Output"

    (**Workflow output**) Empty text file (`.txt`)


### `modify_per_condition_ascii_pileups`

Modify the generated ASCII-style pileups for all the desired annotated regions
across the different library subsets if provided with a
[**custom script**][custom-script-ascii-mod].

> **OPTIONAL RULE.** The ASCII-style pileups for each annotated region are
> modified if, and only if, at least one library subset is specified in the
> [configuration file](#configuration-file). Otherwise, this rule will not be
> executed, and no output will be generated.

=== "Condition"

    - **config_template.yaml**
        - `lib_dict`: Dictionary of arbitrary condition names (keys) and library
        names to aggregate alignment pileups for (values; MUST correspond to names
        in samples table) (default: None)

=== "Input"

    (**Workflow output**) Empty text file (`.txt`); from
    [**create_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups)

=== "Parameters"

    - **config_template.yaml**
        - `sort_by`: ASCII-style alignment pileups can be sorted by the first
          nucleotide's position from left-to-right ('position') or by counts
          in descending order ('counts') (default: 'position')
        - `split`: Split precursor pileups into one mature-arm pileup per arm
          (default: 'true')
        - `canonical`: Mark the aligned read corresponding to the canonical
          sequence (default: 'true')
        - `keep_all`: Write the pileup even if it has no aligned sequences
          (default: 'true')
        - `min_count_dict`: Dictionary with the minimum count for a sequence to
          be kept for each pileup group (default:
          'min_count_dict["condition"] = 1')
        - `max_seq`: Maximum number of top sequences to display (default: 15)
        - `extension`: Extension of the mature miRNA start and end coordinates
          in bp (default: 6)

=== "Output"

    (**Workflow output**) Empty text file (`.txt`)


### `color_code_ascii_pileups`

Color-code all the modified ASCII-style pileups for all the desired annotated
regions with a [**custom script**][custom-script-copper].

> A directory containing the color-coded ASCII-style pileups is created. If no
> ASCII-style alignment pileups were created, the color-coded pileups' output
> directories will only contain an empty file and their corresponding CSS style
> file.

===  "Input"

    (**Workflow output**) Empty text file (`.txt`); from
    [**modify_per_library_ascii_pileups**](pileups.md#create_per_library_ascii_pileups),
    [**modify_per_run_ascii_pileups**](pileups.md#create_per_run_ascii_pileups), and
    [**modify_per_condition_ascii_pileups**](pileups.md#create_per_condition_ascii_pileups),

=== "Parameters"

      - **config_template.yaml**
          - `keep_info`: keep genomic coordinates and feature names on the
            final file (default: 'true')
          - `color_dict`: Dictionary with the character-to-color mapping. See
            available colors in the
            [module overview](../overview.md#ascii-style-alignmnet-pileups-module).
            (default: 'adenine'='green', 'cytosine'='orange',
            'guanine'='"light purple"', 'thymine'='"light blue"',
            'gap'='"light gray"', and 'generic'='white')

=== "Output"

    (**Workflow output**) Empty text file (`.txt`)
