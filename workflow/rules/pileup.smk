###############################################################################
# (c) 2024 Iris Mestres, Zavolan Lab, Biozentrum, University of Basel
# (@) zavolab-biozentrum@unibas.ch
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
        piles_raw_run=PILEUP_DIR / "raw/all/check_file.txt",
        piles_raw_lib=expand(
            PILEUP_DIR / "raw" / "{sample}" / "check_file.txt",
            sample=pd.unique(samples_table.index.values),
        ),
        piles_mod_run=PILEUP_DIR / "mod/all/check_file.txt",
        piles_mod_lib=expand(
            PILEUP_DIR / "mod" / "{sample}" / "check_file.txt",
            sample=pd.unique(samples_table.index.values),
        ),


###############################################################################
### Generate empty BED file if not user-provided
###############################################################################


if config["bed_file"] == "":

    rule create_empty_bed:
        output:
            create_empty_bed_file(config, INTERMEDIATES_DIR),
        log:
            LOCAL_LOG / "create_empty_bed.log",
        container:
            "docker://ubuntu:lunar-20221207"
        params:
            cluster_log=CLUSTER_LOG / "create_empty_bed.log",
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
    log:
        LOCAL_LOG / "compress_reference_genome.log",
    conda:
        ENV_DIR / "samtools.yaml"
    container:
        "docker://quay.io/biocontainers/samtools:1.21--h96c455f_1"
    params:
        cluster_log=CLUSTER_LOG / "compress_reference_genome.log",
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
        piles=PILEUP_DIR / "raw" / "{sample}" / "check_file.txt",
    log:
        LOCAL_LOG / "pileups_raw_{sample}.log",
    conda:
        ENV_DIR / "r.yaml"
    container:
        "docker://zavolab/ascii-alignment-pileup:1.1.1"
    params:
        cluster_log=CLUSTER_LOG / "pileups_raw_{sample}.log",
        out_dir=lambda wildcards: expand(
            PILEUP_DIR / "raw" / "{sample}", sample=[wildcards.sample]
        ),
        prefix="{sample}",
        sort=config["sort_by"],
    shell:
        "(touch {output.piles} && Rscript {input.script} \
        --verbose \
        --annotations={input.annotations} \
        --reference={input.reference} \
        --sort-by={params.sort} \
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
        piles=PILEUP_DIR / "raw/all/check_file.txt",
    log:
        LOCAL_LOG / "pileups_raw_whole_run.log",
    conda:
        ENV_DIR / "r.yaml"
    container:
        "docker://zavolab/ascii-alignment-pileup:1.1.1"
    resources:
        mem=16,
    params:
        cluster_log=CLUSTER_LOG / "pileups_raw_whole_run.log",
        out_dir=PILEUP_DIR / "raw" / "all",
        prefix="all_samples",
        sort=config["sort_by"],
    shell:
        "(touch {output.piles} && Rscript {input.script} \
        --verbose \
        --annotations={input.annotations} \
        --reference={input.reference} \
        --prefix={params.prefix} \
        --sort-by={params.sort} \
        --output-directory {params.out_dir} \
        {input.regions} \
        {input.maps} \
        ) &> {log}"


###############################################################################
### Generate ASCII-style pileups (per experiment design)
###############################################################################

if config["lib_dict"] != None:
    cond = list(config["lib_dict"].keys())

    rule create_per_condition_ascii_pileups:
        input:
            annotations=INTERMEDIATES_DIR / "mirna_annotations.gff3",
            maps=lambda wildcards: expand(
                OUT_DIR
                / "{group}"
                / "alignments_intersecting_mirna_uncollapsed_sorted.bam",
                group=config["lib_dict"][wildcards.cond],
            ),
            maps_index=lambda wildcards: expand(
                OUT_DIR
                / "{group}"
                / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
                group=config["lib_dict"][wildcards.cond],
            ),
            reference=INTERMEDIATES_DIR / "genome_processed.fa.bz",
            regions=config["bed_file"],
            script=SCRIPTS_DIR / "ascii_alignment_pileup.R",
        output:
            piles=PILEUP_DIR / "raw" / "{cond}" / "check_file_{cond}.txt",
        log:
            LOCAL_LOG / "pileups_raw_condition_{cond}.log",
        conda:
            ENV_DIR / "r.yaml"
        container:
            "docker://zavolab/ascii-alignment-pileup:1.1.1"
        params:
            cluster_log=CLUSTER_LOG / "pileups_raw_condition_{cond}.log",
            out_dir=lambda wildcards: expand(
                PILEUP_DIR / "raw" / "{cond}", cond=wildcards.cond
            ),
            prefix="{cond}",
            sort=config["sort_by"],
        shell:
            "(touch {output.piles} && Rscript {input.script} \
            --verbose \
            --annotations={input.annotations} \
            --reference={input.reference} \
            --prefix={params.prefix} \
            --sort-by={params.sort} \
            --output-directory {params.out_dir} \
            {input.regions} \
            {input.maps} \
            ) &> {log}"


###############################################################################
### Modify ASCII-style pileups (per library)
###############################################################################


rule modify_per_library_ascii_pileups:
    input:
        piles=PILEUP_DIR / "raw" / "{sample}" / "check_file.txt",
        script=SCRIPTS_DIR / "ascii_pileups_aesthetics_modification.R",
    output:
        piles=PILEUP_DIR / "mod" / "{sample}" / "check_file.txt",
    log:
        LOCAL_LOG / "pileups_mod_{sample}.log",
    conda:
        ENV_DIR / "r.yaml"
    container:
        "docker://zavolab/r-tidyverse:3.5.3"
    params:
        cluster_log=CLUSTER_LOG / "pileups_mod_{sample}.log",
        in_dir=lambda wildcards: expand(
            PILEUP_DIR / "raw" / "{sample}", sample=[wildcards.sample]
        ),
        out_dir=lambda wildcards: expand(
            PILEUP_DIR / "mod" / "{sample}", sample=[wildcards.sample]
        ),
        prefix="{sample}",
        split_arms=lambda wc: "--split-arms" if config["split"] else "",
        canonical=lambda wc: "--canonical" if config["canonical"] else "",
        keep_all=lambda wc: "--keep-all" if config["keep_all"] else "",
        min_count=config["min_count_dict"]["lib"],
        max_seq=config["max_seq"],
        overhang=config["extension"],
    shell:
        "(touch {output.piles} && Rscript {input.script} \
        --verbose \
        --in-dir={params.in_dir} \
        --prefix={params.prefix} \
        --out-dir {params.out_dir} \
        --min-count {params.min_count} \
        --max-sequences {params.max_seq} \
        --overhang {params.overhang} \
        {params.split_arms} {params.canonical} {params.keep_all} \
        ) &> {log}"


###############################################################################
### Modify ASCII-style pileups (per run)
###############################################################################


rule modify_per_run_ascii_pileups:
    input:
        piles=PILEUP_DIR / "raw/all/check_file.txt",
        script=SCRIPTS_DIR / "ascii_pileups_aesthetics_modification.R",
    output:
        piles=PILEUP_DIR / "mod/all/check_file.txt",
    log:
        LOCAL_LOG / "pileups_mod_whole_run.log",
    conda:
        ENV_DIR / "r.yaml"
    container:
        "docker://zavolab/r-tidyverse:3.5.3"
    resources:
        mem=16,
    params:
        cluster_log=CLUSTER_LOG / "pileups_mod_whole_run.log",
        in_dir=PILEUP_DIR / "raw" / "all",
        out_dir=PILEUP_DIR / "mod" / "all",
        prefix="all-samples",
        split_arms=lambda wc: "--split-arms" if config["split"] else "",
        canonical=lambda wc: "--canonical" if config["canonical"] else "",
        keep_all=lambda wc: "--keep-all" if config["keep_all"] else "",
        min_count=config["min_count_dict"]["run"],
        max_seq=config["max_seq"],
        overhang=config["extension"],
    shell:
        "(touch {output.piles} && Rscript {input.script} \
        --verbose \
        --in-dir={params.in_dir} \
        --prefix={params.prefix} \
        --out-dir {params.out_dir} \
        --min-count {params.min_count} \
        --max-sequences {params.max_seq} \
        --overhang {params.overhang} \
        {params.split_arms} {params.canonical} {params.keep_all} \
        ) &> {log}"


###############################################################################
### Modify ASCII-style pileups (per experiment design)
###############################################################################

if config["lib_dict"] != None:
    cond = list(config["lib_dict"].keys())

    rule modify_per_condition_ascii_pileups:
        input:
            piles=PILEUP_DIR / "raw" / "{cond}" / "check_file_{cond}.txt",
            script=SCRIPTS_DIR / "ascii_pileups_aesthetics_modification.R",
        output:
            piles=PILEUP_DIR / "mod" / "{cond}" / "check_file_{cond}.txt",
        log:
            LOCAL_LOG / "pileups_mod_condition_{cond}.log",
        conda:
            ENV_DIR / "r.yaml"
        container:
            "docker://zavolab/r-tidyverse:3.5.3"
        params:
            cluster_log=CLUSTER_LOG / "pileups_mod_condition_{cond}.log",
            in_dir=lambda wildcards: expand(
                PILEUP_DIR / "raw" / "{cond}", cond=wildcards.cond
            ),
            out_dir=lambda wildcards: expand(
                PILEUP_DIR / "mod" / "{cond}", cond=wildcards.cond
            ),
            prefix="{cond}",
            split_arms=lambda wc: "--split-arms" if config["split"] else "",
            canonical=lambda wc: "--canonical" if config["canonical"] else "",
            keep_all=lambda wc: "--keep-all" if config["keep_all"] else "",
            min_count=config["min_count_dict"]["condition"],
            max_seq=config["max_seq"],
            overhang=config["extension"],
        shell:
            "(touch {output.piles} && Rscript {input.script} \
            --verbose \
            --in-dir={params.in_dir} \
            --prefix={params.prefix} \
            --out-dir {params.out_dir} \
            --min-count {params.min_count} \
            --max-sequences {params.max_seq} \
            --overhang {params.overhang} \
            {params.split_arms} {params.canonical} {params.keep_all} \
            ) &> {log}"
