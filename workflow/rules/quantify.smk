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
        uncollapsed_bam=expand(
            os.path.join(
                config["output_dir"],
                "{sample}",
                "uncollapsedSortedMappings_{sample}.bam",
            ),
            sample=pd.unique(samples_table.index.values),
        ),
        uncollapsed_bai=expand(
            os.path.join(
                config["output_dir"],
                "{sample}",
                "uncollapsedSortedMappings_{sample}.bam.bai",
            ),
            sample=pd.unique(samples_table.index.values),
        )


###############################################################################
### Intersection with extended pri-miRs
###############################################################################


rule intersect_extended_premir:
    input:
        alignment=os.path.join(
            config["output_dir"],
            "{sample}",
            "convertedSortedMappings_{sample}.bam",
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
        -bed \
        > {output.intersect} \
        ) &> {log}"


###############################################################################
### Filter SAM file with intersected alignments
###############################################################################


rule intersection_sam_file:
    input:
        alignments=os.path.join(
            config["output_dir"], "{sample}", "removeMultimappers.sam"
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
### Convert SAM to BAM
###############################################################################


rule convert_intersected_premir_to_bam:
    input:
        maps=os.path.join(
            config["output_dir"], 
            "{sample}", 
            "intersectedAlignments.sam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"], 
            "{sample}", 
            "intersected_premirAlignments.bam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], 
            "convert_intersected_premir_to_bam_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], 
            "convert_intersected_premir_to_bam_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_intersected_premir_by_position:
    input:
        maps=os.path.join(
            config["output_dir"], 
            "{sample}", 
            "intersected_premirAlignments.bam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"], 
            "{sample}", 
            "sorted_intersected_premirAlignments.bam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "sort_intersected_premir_by_position_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"],
            "sort__intersected_premir_by_position_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_intersected_premir_bam:
    input:
        maps=os.path.join(
            config["output_dir"], 
            "{sample}", 
            "sorted_intersected_premirAlignments.bam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"], 
            "{sample}", 
            "sorted_intersected_premirAlignments.bam.bai"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "index_bam_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "index_bam_{sample}.log"),
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
            "sorted_intersected_premirAlignments.bam"
        ),
        premir=os.path.join(
            config["output_dir"], "mirna_extended_annotations.bed"
        ),
    output:
        intersect=os.path.join(
            config["output_dir"], "{sample}", 
            "intersection_extended_mirna.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], 
            "intersection_extended_mirna_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], 
        "intersection_extended_mirna_{sample}.log"),
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
        -bed \
        > {output.intersect} \
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
        script=os.path.join(config["scripts_dir"], "primir_quantification.py"),
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


###############################################################################
### Uncollapse reads
###############################################################################


rule uncollapse_reads:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "intersectedAlignments.sam"
        ),
        script=os.path.join(config["scripts_dir"], "sam_uncollapse.pl"),
    output:
        maps=os.path.join(
            config["output_dir"], "{sample}", "uncollapsedMappings.sam"
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


rule uncollpased_convert_to_bam:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "uncollapsedMappings.sam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"], "{sample}", "uncollapsedMappings.bam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "convert_uncollpased_to_bam_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], 
            "convert_uncollapsed_to_bam_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_uncollpased_by_position:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "uncollapsedMappings.bam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "uncollapsedSortedMappings_{sample}.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sort_uncollapsed_by_position_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], 
            "sort_uncollapsed_by_position_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule uncollapsed_index_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "uncollapsedSortedMappings_{sample}.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "uncollapsedSortedMappings_{sample}.bam.bai",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "uncollapsed_index_bam_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], 
            "uncollapsed_index_bam_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"

