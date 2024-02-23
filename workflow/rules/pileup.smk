###############################################################################
# (c) 2024 Iris Mestres, Zavolan Lab, Biozentrum, University of Basel
# (@) iris.mestres@alumn.esci.upf.edu
#
# Workflow to create ASCII-style pileups of read alignments.
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
INTERMEDIATES_DIR = Path(config["intermediates_dir"])
OUT_DIR = Path(config["output_dir"])
PILEUP_DIR = Path(config["pileups_dir"])
SCRIPTS_DIR = Path(config["scripts_dir"])

CLUSTER_LOG = Path(config["cluster_log"])
LOCAL_LOG = Path(config["local_log"])


###############################################################################
### Including functions
###############################################################################


include: "common.smk"


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
    finish_pileup,


###############################################################################
### Finish rule
###############################################################################


rule finish_pileup:
    input:
        piles_run=PILEUP_DIR / "all/check_file.txt",
        piles_lib=expand(
            PILEUP_DIR / "{sample}" / "check_file.txt",
            sample=pd.unique(samples_table.index.values),
        ),


###############################################################################
### Generate empty BED file if not user-provided
###############################################################################


if config["bed_file"] == "":

    rule create_empty_bed:
        output:
            create_empty_bed_file(config, INTERMEDIATES_DIR),
        params:
            cluster_log=CLUSTER_LOG / "create_empty_bed.log",
        log:
            LOCAL_LOG / "create_empty_bed.log",
        container:
            "docker://ubuntu:lunar-20221207"
        shell:
            "(touch {output})"


###############################################################################
### Compress reference genome with trimmed IDs
###############################################################################


rule compress_reference_genome:
    input:
        genome=INTERMEDIATES_DIR / "genome_processed.fa",
    output:
        genome=INTERMEDIATES_DIR / "genome_processed.fa.bz",
    params:
        cluster_log=CLUSTER_LOG / "compress_reference_genome.log",
    log:
        LOCAL_LOG / "compress_reference_genome.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(bgzip < {input.genome} > {output.genome}) &> {log}"


###############################################################################
### Generate ASCII-style pileups (per library)
###############################################################################


rule create_per_library_ascii_pileups:
    input:
        annotations=INTERMEDIATES_DIR / "mirna_annotations.gff3",
        maps=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
        maps_index=OUT_DIR
        / "{sample}"
        / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
        reference=INTERMEDIATES_DIR / "genome_processed.fa.bz",
        regions=config["bed_file"],
        script=SCRIPTS_DIR / "ascii_alignment_pileup.R",
    output:
        piles=PILEUP_DIR / "{sample}" / "check_file.txt",
    params:
        cluster_log=CLUSTER_LOG / "pileups_{sample}.log",
        out_dir=lambda wildcards: expand(
            PILEUP_DIR / "{sample}", sample=pd.unique(samples_table.index.values)
        ),
        prefix="{sample}",
    log:
        LOCAL_LOG / "pileups_{sample}.log",
    container:
        "docker://zavolab/ascii-alignment-pileup:1.1.1"
    conda:
        ENV_DIR / "r.yaml"
    shell:
        "(touch {output.piles} && Rscript {input.script} \
        --verbose \
        --annotations={input.annotations} \
        --reference={input.reference} \
        --prefix={params.prefix} \
        --output-directory {params.out_dir} \
        {input.regions} \
        {input.maps} \
        ) &> {log}"


###############################################################################
### Generate ASCII-style pileups (per run)
###############################################################################


rule create_per_run_ascii_pileups:
    input:
        annotations=INTERMEDIATES_DIR / "mirna_annotations.gff3",
        maps=expand(
            OUT_DIR
            / "{sample}"
            / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
            sample=pd.unique(samples_table.index.values),
        ),
        maps_index=expand(
            OUT_DIR
            / "{sample}"
            / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
            sample=pd.unique(samples_table.index.values),
        ),
        reference=INTERMEDIATES_DIR / "genome_processed.fa.bz",
        regions=config["bed_file"],
        script=SCRIPTS_DIR / "ascii_alignment_pileup.R",
    output:
        piles=PILEUP_DIR / "all/check_file.txt",
    params:
        cluster_log=CLUSTER_LOG / "pileups_whole_run.log",
        out_dir=PILEUP_DIR / "all",
        prefix="all_samples",
    log:
        LOCAL_LOG / "pileups_whole_run.log",
    container:
        "docker://zavolab/ascii-alignment-pileup:1.1.1"
    conda:
        ENV_DIR / "r.yaml"
    shell:
        "(touch {output.piles} && Rscript {input.script} \
        --verbose \
        --annotations={input.annotations} \
        --reference={input.reference} \
        --prefix={params.prefix} \
        --output-directory {params.out_dir} \
        {input.regions} \
        {input.maps} \
        ) &> {log}"


###############################################################################
### Generate ASCII-style pileups (per experiment design)
###############################################################################

if config["lib_dict"] != None:
    condition = list(config["lib_dict"].keys())

    rule create_per_condition_ascii_pileups:
        input:
            annotations=INTERMEDIATES_DIR / "mirna_annotations.gff3",
            maps=lambda wildcards: expand(
                OUT_DIR
                / "{group}"
                / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
                group=config["lib_dict"][wildcards.condition],
            ),
            maps_index=lambda wildcards: expand(
                OUT_DIR
                / "{group}"
                / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
                group=config["lib_dict"][wildcards.condition],
            ),
            reference=INTERMEDIATES_DIR / "genome_processed.fa.bz",
            regions=config["bed_file"],
            script=SCRIPTS_DIR / "ascii_alignment_pileup.R",
        output:
            piles=PILEUP_DIR / "{condition}" / "check_file_{condition}.txt",
        params:
            cluster_log=CLUSTER_LOG / "pileups_condition_{condition}.log",
            out_dir=lambda wildcards: expand(
                PILEUP_DIR / "{condition}", condition=wildcards.condition
            ),
            prefix="{condition}",
        log:
            LOCAL_LOG / "pileups_condition_{condition}.log",
        container:
            "docker://zavolab/ascii-alignment-pileup:1.1.1"
        conda:
            ENV_DIR / "r.yaml"
        shell:
            "(touch {output.piles} && Rscript {input.script} \
            --verbose \
            --annotations={input.annotations} \
            --reference={input.reference} \
            --prefix={params.prefix} \
            --output-directory {params.out_dir} \
            {input.regions} \
            {input.maps} \
            ) &> {log}"
