###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Pipeline to quantify miRNAs, including isomiRs, from miRNA-seq alignments.
###############################################################################

import pandas as pd
from snakemake.utils import validate

from pathlib import Path


###############################################################################
### Configuration validation
###############################################################################

validate(config, Path("../../config/config_schema.json"))


###############################################################################
### Paths configuration
###############################################################################


ENV_DIR = Path(f"{workflow.basedir}/envs")
OUT_DIR = Path(config["output_dir"])
SCRIPTS_DIR = Path(config["scripts_dir"])

CLUSTER_LOG = Path(config["cluster_log"])
LOCAL_LOG = Path(config["local_log"])


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
        primir_intersect_sam=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_primir.sam",
        mirna_intersect_sam=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna.sam",
        intersect_sam=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_sorted_tag.sam",
        table=OUT_DIR / "TABLES" / "all_{mir}_counts.tab",
        uncollapsed_bam=expand(
            OUT_DIR
            / "{sample}"
            / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
            sample=pd.unique(samples_table.index.values),
        ),
        uncollapsed_bai=expand(
            OUT_DIR
            / "{sample}"
            / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
            sample=pd.unique(samples_table.index.values),
        ),


###############################################################################
### Intersection with extended pri-miRs
###############################################################################


rule intersect_extended_primir:
    input:
        alignment=OUT_DIR / "{sample}" / "alignments_all_sorted_{sample}.bam",
        primir=expand(
            OUT_DIR / "extended_primir_annotation_{extension}_nt.gff3",
            extension=config["extension"],
        ),
    output:
        intersect=OUT_DIR / "{sample}" / "intersected_extended_primir.bed",
    params:
        cluster_log=CLUSTER_LOG / "intersect_extended_primir_{sample}.log",
    log:
        LOCAL_LOG / "intersect_extended_primir_{sample}.log",
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    conda:
        ENV_DIR / "bedtools.yaml"
    shell:
        "(bedtools intersect \
        -wb \
        -s \
        -F 1 \
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
        alignments=OUT_DIR / "{sample}" / "alignments_all.sam",
        intersect=OUT_DIR / "{sample}" / "intersected_extended_primir.bed",
    output:
        sam=OUT_DIR / "{sample}" / "alignments_intersecting_primir.sam",
    params:
        cluster_log=CLUSTER_LOG
        / "filter_sam_by_intersecting_primir_{sample}.log",
    log:
        LOCAL_LOG / "filter_sam_by_intersecting_primir_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
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
        maps=OUT_DIR / "{sample}" / "alignments_intersecting_primir.sam",
    output:
        maps=temp(OUT_DIR / "{sample}" / "alignments_intersecting_primir.bam"),
    params:
        cluster_log=CLUSTER_LOG
        / "convert_intersecting_primir_sam_to_bam_{sample}.log",
    log:
        LOCAL_LOG / "convert_intersecting_primir_sam_to_bam_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_intersecting_primir_bam_by_position:
    input:
        maps=OUT_DIR / "{sample}" / "alignments_intersecting_primir.bam",
    output:
        maps=temp(OUT_DIR / "{sample}" / "alignments_intersecting_primir_sorted.bam"),
    params:
        cluster_log=CLUSTER_LOG
        / "sort_intersecting_primir_bam_by_position_{sample}.log",
    log:
        LOCAL_LOG / "sort_intersecting_primir_bam_by_position_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools sort -n {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_intersecting_primir_bam:
    input:
        maps=OUT_DIR / "{sample}" / "alignments_intersecting_primir_sorted.bam",
    output:
        maps=temp(OUT_DIR
        / "{sample}"
        / "alignments_intersecting_primir_sorted.bam.bai"),
    params:
        cluster_log=CLUSTER_LOG / "index_intersecting_primir_bam_{sample}.log",
    log:
        LOCAL_LOG / "index_intersecting_primir_bam_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Intersection with extended miRNAs
###############################################################################


rule intersect_extended_mirna:
    input:
        alignment=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_primir_sorted.bam",
        mirna=expand(
            OUT_DIR / "extended_mirna_annotation_{extension}_nt.gff3",
            extension=config["extension"],
        ),
    output:
        intersect=OUT_DIR / "{sample}" / "intersected_extended_mirna.bed",
    params:
        cluster_log=CLUSTER_LOG / "intersect_extended_mirna_{sample}.log",
    log:
        LOCAL_LOG / "intersect_extended_mirna_{sample}.log",
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    conda:
        ENV_DIR / "bedtools.yaml"
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
        alignments=OUT_DIR / "{sample}" / "alignments_intersecting_primir.sam",
        intersect=OUT_DIR / "{sample}" / "intersected_extended_mirna.bed",
    output:
        sam=OUT_DIR / "{sample}" / "alignments_intersecting_mirna.sam",
    params:
        cluster_log=CLUSTER_LOG
        / "filter_sam_by__intersecting_mirna_{sample}.log",
    log:
        LOCAL_LOG / "filter_sam_by_intersecting_mirna_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
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
        alignments=OUT_DIR / "{sample}" / "alignments_intersecting_mirna.sam",
        intersect=OUT_DIR / "{sample}" / "intersected_extended_mirna.bed",
        script=SCRIPTS_DIR / "iso_name_tagging.py",
    output:
        sam=OUT_DIR / "{sample}" / "alignments_intersecting_mirna_tag.sam",
    params:
        extension=config["extension"],
        cluster_log=CLUSTER_LOG / "add_intersecting_mirna_tag_{sample}.log",
    log:
        LOCAL_LOG / "add_intersecting_mirna_tag_{sample}.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    conda:
        ENV_DIR / "pysam.yaml"
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
        sam=OUT_DIR / "{sample}" / "alignments_intersecting_mirna_tag.sam",
    output:
        sam=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_sorted_tag.sam",
    params:
        cluster_log=CLUSTER_LOG
        / "sort_intersecting_mirna_by_feat_tag_{sample}.log",
    log:
        LOCAL_LOG / "sort_intersecting_mirna_by_feat_tag_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools sort -t YW -O SAM {input.sam} > {output.sam}) &> {log}"


###############################################################################
### miRNAs counting table - miRNA
###############################################################################


rule quantify_mirna:
    input:
        alignments=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_sorted_tag.sam",
        script=SCRIPTS_DIR / "mirna_quantification.py",
    output:
        table=temp(OUT_DIR / "TABLES" / "mirna_counts_{sample}"),
    params:
        cluster_log=CLUSTER_LOG / "quantify_mirna_{sample}.log",
        mir_list=config["mir_list"],
        library="{sample}",
        out_dir=OUT_DIR / "TABLES",
    log:
        LOCAL_LOG / "quantify_mirna_{sample}.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0"
    conda:
        ENV_DIR / "pysam.yaml"
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
        intersect=OUT_DIR / "{sample}" / "intersected_extended_primir.bed",
        script=SCRIPTS_DIR / "primir_quantification.py",
    output:
        table=temp(OUT_DIR / "TABLES" / "pri-mir_counts_{sample}"),
    params:
        cluster_log=CLUSTER_LOG / "quantify_primir_{sample}.log",
    log:
        LOCAL_LOG / "quantify_primir_{sample}.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0"
    conda:
        ENV_DIR / "pysam.yaml"
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
            OUT_DIR / "TABLES" / "{mir}_counts_{sample}",
            sample=pd.unique(samples_table.index.values),
            mir=[mir for mir in config["mir_list"] if mir != "isomir"],
        ),
        script=SCRIPTS_DIR / "merge_tables.R",
    output:
        table=OUT_DIR / "TABLES" / "all_{mirna}_counts.tab",
    params:
        cluster_log=CLUSTER_LOG / "merge_tables_{mirna}.log",
        prefix="{mirna}_counts_",
        input_dir=OUT_DIR / "TABLES",
    log:
        LOCAL_LOG / "merge_tables_{mirna}.log",
    container:
        "docker://zavolab/r-tidyverse:3.5.3"
    conda:
        ENV_DIR / "r.yaml"
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
        maps=OUT_DIR / "{sample}" / "alignments_intersecting_mirna.sam",
        script=SCRIPTS_DIR / "sam_uncollapse.pl",
    output:
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed.sam",
    params:
        cluster_log=CLUSTER_LOG / "uncollapse_reads_{sample}.log",
    log:
        LOCAL_LOG / "uncollapse_reads_{sample}.log",
    container:
        "docker://perl:5.37.10"
    conda:
        ENV_DIR / "perl.yaml"
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
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed.sam",
    output:
        maps=temp(OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed.bam"),
    params:
        cluster_log=CLUSTER_LOG
        / "convert_uncollapsed_reads_sam_to_bam_{sample}.log",
    log:
        LOCAL_LOG / "convert_uncollapsed_reads_sam_to_bam_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_uncollpased_reads_bam_by_position:
    input:
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed.bam",
    output:
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
    params:
        cluster_log=CLUSTER_LOG
        / "sort_uncollapsed_reads_bam_by_position_{sample}.log",
    log:
        LOCAL_LOG / "sort_uncollapsed_reads_bam_by_position_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_uncollapsed_reads_bam:
    input:
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
    output:
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
    params:
        cluster_log=CLUSTER_LOG / "index_uncollapsed_reads_bam_{sample}.log",
    log:
        LOCAL_LOG / "index_uncollapsed_reads_bam_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"
