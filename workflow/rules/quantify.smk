###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Pipeline to quantify miRNAs, including isomiRs, from miRNA-seq alignments.
###############################################################################
#
# USAGE (from the file's directory):
#
# snakemake \
#    --snakefile="quanitfy.smk" \
#    --cores 4 \
#    --use-singularity \
#    --singularity-args "--bind $PWD/../" \ 
#    --printshellcmds \
#    --rerun-incomplete \
#    --verbose
#
# IMPORTANT when executing this file alone:
## * You must modify the config.yaml.
## * Uncomment the configfile line.
################################################################################

import os
import pandas as pd

###############################################################################
### Reading samples' table
###############################################################################

samples_table = pd.read_csv(
    config["samples"],
    header = 0,
    index_col = 0,
    comment = "#",
    engine = "python",
    sep = "\t",
)

# Rules that require internet connection for downloading files are included
# in the localrules
localrules:
    finish_quantify,


###############################################################################
### Finish rule
###############################################################################


rule finish_quantify:
    input:
        table1=expand(
            os.path.join(
                config["output_dir"], 
                "TABLES", 
                "counts.{mir}.tab",
            ),
            mir=config["mir_list"],
        ),
        table2=os.path.join(
            config["output_dir"], 
            "TABLES", 
            "counts.isomirs.tab",
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
        cluster_log=os.path.join(
            config["cluster_log"], "bamtobed_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "bamtobed_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.27.1--hd03093a_6"
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
    singularity:
        "docker://ubuntu:bionic-20221215"
    shell:
        "(sort \
        -k1,1 \
        -k2,2n \
        {input.alignment} \
        > {output.alignment} \
        ) &> {log}"


###############################################################################
### miRNAs intersection
###############################################################################


rule intersect_mirna:
    input:
        alignment=os.path.join(
            config["output_dir"], "{sample}", "sorted.alignments.bed12"
        ),
        mirna=os.path.join(
            config["output_dir"], "mirna_annotations.bed"
        ),
    output:
        intersect=os.path.join(
            config["output_dir"], "{sample}", "intersect_mirna.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "intersection_mirna_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "intersection_mirna_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.27.1--hd03093a_6"
    shell:
        "(bedtools intersect \
        -wao \
        -s \
        -F 1 \
        -sorted \
        -b {input.alignment} \
        -a {input.mirna} \
        > {output.intersect} \
        ) &> {log}"


###############################################################################
### isomiRs intersection
###############################################################################


# rule intersect_isomirs:
#     input:
#         alignment=os.path.join(
#             config["output_dir"], "{sample}", "sorted.alignments.bed12"
#         ),
#         isomirs=os.path.join(
#            config["output_dir"], "isomirs_annotation.bed",
#         ),
#     output:
#         intersect=os.path.join(
#             config["output_dir"], "{sample}", "intersect_isomirs.bed"
#         ),
#     params:
#         cluster_log=os.path.join(
#             config["cluster_log"], "intersection_isomirs_{sample}.log"
#         ),
#     log:
#         os.path.join(config["local_log"], "intersection_isomirs_{sample}.log"),
#     singularity:
#         "docker://quay.io/biocontainers/bedtools:2.27.1--hd03093a_6"
#     shell:
#         "(bedtools intersect \
#         -wao \
#         -s \
#         -F 1 \
#         -sorted \
#         -b {input.alignment} \
#         -a {input.isomirs} \
#         > {output.intersect} \
#         ) &> {log}"


###############################################################################
### miRNAs counting table - miRNA
###############################################################################


rule quant_mirna:
    input:
        intersect=os.path.join(
            config["output_dir"], "{sample}", "intersect_mirna.bed"
        ),
        script=os.path.join(config["scripts_dir"], "mirna_quantification.py"),
    output:
        table=os.path.join(
            config["output_dir"], "TABLES", "miRNA_counts_{sample}"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "quant_mirna_miRNA_{sample}.log"
        ),
        prefix=os.path.join(
            config["output_dir"], "TABLES", "miRNA_counts_{sample}"
        ),
    log:
        os.path.join(config["local_log"], "quant_mirna_miRNA_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python \
        {input.script} \
        -i {input.intersect} \
        --uniq=miRNA \
        -p={params.prefix} \
        ) &> {log}"


###############################################################################
### miRNAs counting table - miRNA_primary
###############################################################################


rule quant_mirna_pri:
    input:
        intersect=os.path.join(
            config["output_dir"], "{sample}", "intersect_mirna.bed"
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
        prefix=os.path.join(
            config["output_dir"],
            "TABLES",
            "miRNA_primary_transcript_counts_{sample}",
        ),
    log:
        os.path.join(
            config["local_log"],
            "quant_mirna_miRNA_primary_transcript_{sample}.log",
        ),
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python \
        {input.script} \
        -i {input.intersect} \
        --uniq=miRNA_primary_transcript \
        -p={params.prefix} \
        ) &> {log}"


###############################################################################
### isomiRs counting table
###############################################################################


# rule quant_isomirs:
#     input:
#         intersect=os.path.join(
#             config["output_dir"], "{sample}", "intersect_isomirs.bed"
#         ),
#         script=os.path.join(config["scripts_dir"], "mirna_quantification.py"),
#     output:
#         table=os.path.join(
#             config["output_dir"], "TABLES", "isomirs_counts_{sample}"
#         ),
#     params:
#         cluster_log=os.path.join(
#             config["cluster_log"], "quant_isomirs_{sample}.log"
#         ),
#         prefix=os.path.join(
#             config["output_dir"], "TABLES", "isomirs_counts_{sample}"
#         ),
#     log:
#         os.path.join(config["local_log"], "quant_isomirs_{sample}.log"),
#     singularity:
#         "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
#     shell:
#         "(python \
#         {input.script} \
#         -i {input.intersect} \
#         --uniq=miRNA \
#         -p={params.prefix} \
#         ) &> {log}"


###############################################################################
### Merge counting tables for all samples by mature/primary/isomirs forms.
###############################################################################


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
        input_dir=os.path.join(config["output_dir"], "TABLES"),
    log:
        os.path.join(config["local_log"], "merge_tables_{mir}.log"),
    singularity:
        "docker://zavolab/r-tidyverse:3.5.3"
    shell:
        "(Rscript \
        {input.script} \
        --input_dir {params.input_dir} \
        --output_file {output.table} \
        --prefix {params.prefix} \
        --verbose \
        ) &> {log}"