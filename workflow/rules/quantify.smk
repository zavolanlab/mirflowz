###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Pipeline to quantify miRNAs, including isomiRs, from miRNA-seq alignments.
###############################################################################

import os
import pandas as pd

from pathlib import Path

###############################################################################
### Reading samples table
###############################################################################

samples_table = pd.read_csv(
    config["samples"],
    header=0,
    index_col=0,
    comment="#",
    engine="python",
    sep="\t",
)


###############################################################################
### Global configuration
###############################################################################


localrules:
    finish_quantify,


###############################################################################
### Finish rule
###############################################################################


rule finish_quantify:
    input:
        primir_intersect_sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir.sam",
        ),
        mirna_intersect_sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna.sam",
        ),
        intersect_sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_sorted_tag.sam",
        ),
        table=os.path.join(
            config["output_dir"],
            "TABLES",
            "all_{mir}_counts.tab",
        ),
        uncollapsed_bam=expand(
            os.path.join(
                config["output_dir"],
                "{sample}",
                "alignments_intersecting_mirna_uncollapsed_sorted.bam",
            ),
            sample=pd.unique(samples_table.index.values),
        ),
        uncollapsed_bai=expand(
            os.path.join(
                config["output_dir"],
                "{sample}",
                "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
            ),
            sample=pd.unique(samples_table.index.values),
        ),


###############################################################################
### Intersection with extended pri-miRs
###############################################################################


rule intersect_extended_primir:
    input:
        alignment=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_all_sorted_{sample}.bam",
        ),
        primir=expand(
            os.path.join(
                config["output_dir"],
                "extended_primir_annotation_{extension}_nt.gff3",
            ),
            extension=config["extension"],
        ),
    output:
        intersect=os.path.join(
            config["output_dir"],
            "{sample}",
            "intersected_extended_primir.bed",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "intersect_extended_primir_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "intersect_extended_primir_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    shell:
        "(bedtools intersect \
        -wb \
        -s \
        -F 1 \
        -sorted \
        -b {input.alignment} \
        -a {input.primir} \
        -bed \
        > {output.intersect} \
        ) &> {log}"


###############################################################################
### Filter SAM file with intersecting alignments (pri-miR)
###############################################################################


rule filter_sam_by_intersecting_primir:
    input:
        alignments=os.path.join(
            config["output_dir"], "{sample}", "alignments_all.sam"
        ),
        intersect=os.path.join(
            config["output_dir"],
            "{sample}",
            "intersected_extended_primir.bed",
        ),
    output:
        sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "filter_sam_by_intersecting_primir_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "filter_sam_by_intersecting_primir_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "((samtools view \
        -H {input.alignments}; \
        awk 'NR==FNR {{bed[$13]=1; next}} $1 in bed' \
        {input.intersect} {input.alignments} \
        ) > {output.sam} \
        ) &> {log}"


###############################################################################
### Convert SAM to BAM
###############################################################################


rule convert_intersecting_primir_sam_to_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir.sam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "convert_intersecting_primir_sam_to_bam_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "convert_intersecting_primir_sam_to_bam_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_intersecting_primir_bam_by_position:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir_sorted.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "sort_intersecting_primir_bam_by_position_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "sort_intersecting_primir_bam_by_position_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort -n {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_intersecting_primir_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir_sorted.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir_sorted.bam.bai",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "index_intersecting_primir_bam_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "index_intersecting_primir_bam_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Intersection with extended miRNAs
###############################################################################


rule intersect_extended_mirna:
    input:
        alignment=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir_sorted.bam",
        ),
        mirna=expand(
            os.path.join(
                config["output_dir"],
                "extended_mirna_annotation_{extension}_nt.gff3",
            ),
            extension=config["extension"],
        ),
    output:
        intersect=os.path.join(
            config["output_dir"], "{sample}", "intersected_extended_mirna.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "intersect_extended_mirna_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "intersect_extended_mirna_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    shell:
        "(bedtools intersect \
        -wo \
        -s \
        -F 1 \
        -b {input.alignment} \
        -a {input.mirna} \
        -bed \
        > {output.intersect} \
        ) &> {log}"


###############################################################################
### Filter SAM file with intersecting alignments (miRNAs)
###############################################################################


rule filter_sam_by_intersecting_mirna:
    input:
        alignments=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_primir.sam",
        ),
        intersect=os.path.join(
            config["output_dir"], "{sample}", "intersected_extended_mirna.bed"
        ),
    output:
        sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "filter_sam_by__intersecting_mirna_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "filter_sam_by_intersecting_mirna_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "((samtools view \
        -H {input.alignments}; \
        awk 'NR==FNR {{bed[$13]=1; next}} $1 in bed' \
        {input.intersect} {input.alignments} \
        ) > {output.sam} \
        ) &> {log}"


###############################################################################
### Add tag with intersecting miRNAs
###############################################################################


rule add_intersecting_mirna_tag:
    input:
        alignments=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna.sam",
        ),
        intersect=os.path.join(
            config["output_dir"], "{sample}", "intersected_extended_mirna.bed"
        ),
        script=os.path.join(config["scripts_dir"], "iso_name_tagging.py"),
    output:
        sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_tag.sam",
        ),
    params:
        extension=config["extension"],
        cluster_log=os.path.join(
            config["cluster_log"], "add_intersecting_mirna_tag_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "add_intersecting_mirna_tag_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        --bed {input.intersect} \
        --sam {input.alignments} \
        --extension {params.extension} \
        > {output.sam} \
        ) &> {log}"


###############################################################################
### Sort by feature tag
###############################################################################


rule sort_intersecting_mirna_by_feat_tag:
    input:
        sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_tag.sam",
        ),
    output:
        sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_sorted_tag.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "sort_intersecting_mirna_by_feat_tag_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "sort_intersecting_mirna_by_feat_tag_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort -t YW -O SAM {input.sam} > {output.sam}) &> {log}"


###############################################################################
### miRNAs counting table - miRNA
###############################################################################


rule quantify_mirna:
    input:
        alignments=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_sorted_tag.sam",
        ),
        script=os.path.join(config["scripts_dir"], "mirna_quantification.py"),
    output:
        table=os.path.join(
            config["output_dir"], "TABLES", "mirna_counts_{sample}"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "quantify_mirna_{sample}.log"
        ),
        mir_list=config["mir_list"],
        library="{sample}",
        out_dir=lambda wildcards, output: Path(output[0]).parent,
    log:
        os.path.join(config["local_log"], "quantify_mirna_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0"
    shell:
        "(python \
        {input.script} \
        {input.alignments} \
        --collapsed \
        --nh \
        --lib {params.library} \
        --outdir {params.out_dir} \
        --mir-list {params.mir_list} \
        ) &> {log}"


################################################################################
#### miRNAs counting table - miRNA_primary
################################################################################


rule quantify_primir:
    input:
        intersect=os.path.join(
            config["output_dir"],
            "{sample}",
            "intersected_extended_primir.bed",
        ),
        script=os.path.join(config["scripts_dir"], "primir_quantification.py"),
    output:
        table=os.path.join(
            config["output_dir"],
            "TABLES",
            "pri-mir_counts_{sample}",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "quantify_primir_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "quantify_primir_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0"
    shell:
        "(python \
        {input.script} \
        {input.intersect} \
        --collapsed \
        --nh \
        > {output.table} \
        ) &> {log}"


################################################################################
#### Merge counting tables for all samples by mature/primary/isomirs forms.
################################################################################


rule merge_tables:
    input:
        table=expand(
            os.path.join(
                config["output_dir"], "TABLES", "{mir}_counts_{sample}"
            ),
            sample=pd.unique(samples_table.index.values),
            mir=[mir for mir in config["mir_list"] if mir != "isomir"],
        ),
        script=os.path.join(config["scripts_dir"], "merge_tables.R"),
    output:
        table=os.path.join(
            config["output_dir"], "TABLES", "all_{mirna}_counts.tab"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "merge_tables_{mirna}.log"
        ),
        prefix="{mirna}_counts_",
        input_dir=lambda wildcards, input: Path(input[0]).parent,
    log:
        os.path.join(config["local_log"], "merge_tables_{mirna}.log"),
    container:
        "docker://zavolab/r-tidyverse:3.5.3"
    shell:
        "(Rscript \
        {input.script} \
        --input_dir {params.input_dir} \
        --output_file {output.table} \
        --prefix {params.prefix} \
        --verbose \
        ) &> {log}"


###############################################################################
### Uncollapse reads
###############################################################################


rule uncollapse_reads:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna.sam",
        ),
        script=os.path.join(config["scripts_dir"], "sam_uncollapse.pl"),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "uncollapse_reads_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "uncollapse_reads_{sample}.log"),
    container:
        "docker://perl:5.37.10"
    shell:
        "(perl {input.script} \
        --suffix \
        --in {input.maps} \
        --out {output.maps} \
        ) &> {log}"


###############################################################################
### Convert SAM to BAM
###############################################################################


rule convert_uncollpased_reads_sam_to_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed.sam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "convert_uncollapsed_reads_sam_to_bam_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "convert_uncollapsed_reads_sam_to_bam_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_uncollpased_reads_bam_by_position:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed_sorted.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "sort_uncollapsed_reads_bam_by_position_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "sort_uncollapsed_reads_bam_by_position_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_uncollapsed_reads_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed_sorted.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "index_uncollapsed_reads_bam_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "index_uncollapsed_reads_bam_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"
