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
        intersect_sam=os.path.join(
            config["output_dir"], "{sample}", "intersectedAlignments.sam"
        ),
        table=expand(
           os.path.join(
               config["output_dir"],
               "TABLES",
               "counts.{mir}.tab",
           ),
           mir=config["mir_list"],
        ),


###############################################################################
### BAM to BED
###############################################################################


rule bamtobed:
    input:
        alignment=os.path.join(
            config["output_dir"],
            "{sample}",
            "convertedSortedMappings_{sample}.bam",
        ),
    output:
        alignment=os.path.join(
            config["output_dir"], "{sample}", "alignments.bed12"
        ),
    params:
        cluster_log=os.path.join(config["cluster_log"], "bamtobed_{sample}.log"),
    log:
        os.path.join(config["local_log"], "bamtobed_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    shell:
        "(bedtools bamtobed \
        -bed12 \
        -tag NH \
        -i {input.alignment} \
        > {output.alignment} \
        ) &> {log}"


###############################################################################
### Sort alignments
###############################################################################


rule sort_alignments:
    input:
        alignment=os.path.join(
            config["output_dir"], "{sample}", "alignments.bed12"
        ),
    output:
        alignment=os.path.join(
            config["output_dir"], "{sample}", "sorted.alignments.bed12"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sortalignment_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "sortalignment_{sample}.log"),
    resources:
        mem=4,
        threads=8,
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(sort \
        -k1,1 \
        -k2,2n \
        {input.alignment} \
        > {output.alignment} \
        ) &> {log}"


###############################################################################
### Intersection with extended pri-miRs
###############################################################################


rule intersect_extended_premir:
    input:
        alignment=os.path.join(
            config["output_dir"], "{sample}", "sorted.alignments.bed12"
        ),
        premir=os.path.join(
            config["output_dir"], "premir_extended_annotations.bed"
        ),
    output:
        intersect=os.path.join(
            config["output_dir"], "{sample}", 
            "intersection_extended_premir.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], 
            "intersection_extended_premir_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], 
        "intersection_extended_premir_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    shell:
        "(bedtools intersect \
        -wb \
        -s \
        -F 1 \
        -sorted \
        -b {input.alignment} \
        -a {input.premir} \
        > {output.intersect} \
        ) &> {log}"


###############################################################################
### Intersection with extended pri-miRs
###############################################################################


rule intersection_sam_file:
    input:
        alignments=os.path.join(
            config["output_dir"], "{sample}", "uncollapsedMappings.sam"
        ),
        intersect=os.path.join(
            config["output_dir"], "{sample}", 
            "intersection_extended_premir.bed"
        ),
    output:
        sam=os.path.join(
            config["output_dir"], "{sample}", "intersectedAlignments.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], 
            "intersection_sam_file_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], 
        "intersection_sam_file_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "((samtools view \
        -H {input.alignments}; \
        awk 'NR==FNR {{bed[$14]=1; next}} $1 in bed' \
        {input.intersect} {input.alignments} \
        ) > {output.sam} \
        ) &> {log}"


###############################################################################
### miRNAs counting table - miRNA
###############################################################################
#
#
#rule quant_mirna:
#   input:
#       intersect=os.path.join(
#           config["output_dir"], "{sample}", "intersect_mirna.bed"
#       ),
#       script=os.path.join(config["scripts_dir"], "mirna_quantification.py"),
#   output:
#       table=os.path.join(
#           config["output_dir"], "TABLES", "miRNA_counts_{sample}"
#       ),
#   params:
#       cluster_log=os.path.join(
#           config["cluster_log"], "quant_mirna_miRNA_{sample}.log"
#       ),
#       prefix=lambda wildcards, output: output[0],
#   log:
#       os.path.join(config["local_log"], "quant_mirna_miRNA_{sample}.log"),
#   container:
#       "docker://quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0"
#   shell:
#       "(python \
#       {input.script} \
#       -i {input.intersect} \
#       --uniq=miRNA \
#       -p={params.prefix} \
#       ) &> {log}"
#
#
################################################################################
#### miRNAs counting table - miRNA_primary
################################################################################


rule quant_mirna_pri:
    input:
        intersect=os.path.join(
            config["output_dir"], "{sample}", 
            "intersection_extended_premir.bed"
        ),
        script=os.path.join(config["scripts_dir"], "mirna_quantification.py"),
    output:
        table=os.path.join(
           config["output_dir"],
           "TABLES",
           "miRNA_primary_transcript_counts_{sample}",
        ),
    params:
        cluster_log=os.path.join(
           config["cluster_log"],
           "quant_mirna_miRNA_primary_transcript_{sample}.log",
        ),
        prefix=lambda wildcards, output: output[0],
    log:
        os.path.join(
           config["local_log"],
           "quant_mirna_miRNA_primary_transcript_{sample}.log",
        ),
    container:
       "docker://quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0"
    shell:
       "(python \
       {input.script} \
       -i {input.intersect} \
       --uniq=miRNA_primary_transcript \
       -p={params.prefix} \
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
           mir=config["mir_list"],
       ),
       script=os.path.join(config["scripts_dir"], "merge_tables.R"),
   output:
       table=os.path.join(config["output_dir"], "TABLES", "counts.{mir}.tab"),
   params:
       cluster_log=os.path.join(
           config["cluster_log"], "merge_tables_{mir}.log"
       ),
       prefix="{mir}_counts_",
       input_dir=lambda wildcards, input: Path(input[0]).parent,
   log:
       os.path.join(config["local_log"], "merge_tables_{mir}.log"),
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

