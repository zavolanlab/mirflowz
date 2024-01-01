###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Snakemake workflow to download and prepare the necessary files
# for smallRNA-seq related workflows.
#
###############################################################################

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
SCRIPTS_DIR = Path(config["scripts_dir"])

CLUSTER_LOG = Path(config["cluster_log"])
LOCAL_LOG = Path(config["local_log"])


###############################################################################
### Global configuration
###############################################################################


localrules:
    finish_prepare,


###############################################################################
### Finish rule
###############################################################################


rule finish_prepare:
    input:
        idx_transcriptome=INTERMEDIATES_DIR / "segemehl_transcriptome_index.idx",
        idx_genome=INTERMEDIATES_DIR / "segemehl_genome_index.idx",
        exons=INTERMEDIATES_DIR / "exons.bed",
        header=INTERMEDIATES_DIR / "genome_header.sam",
        chrsize=INTERMEDIATES_DIR / "chr_size.tsv",
        extended_mir=expand(
            INTERMEDIATES_DIR / "extended_mirna_annotation_{extension}_nt.gff3",
            extension=config["extension"],
        ),
        extended_primir=expand(
            INTERMEDIATES_DIR / "extended_primir_annotation_{extension}_nt.gff3",
            extension=config["extension"],
        ),


###############################################################################
### Trim genome IDs
###############################################################################


rule trim_genome_seq_ids:
    input:
        genome=config["genome_file"],
        script=SCRIPTS_DIR / "trim_id_fasta.sh",
    output:
        genome=INTERMEDIATES_DIR / "genome_processed.fa",
    params:
        cluster_log=CLUSTER_LOG / "genome_process.log",
    log:
        LOCAL_LOG / "genome_process.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(zcat {input.genome} | {input.script} > {output.genome}) &> {log}"


###############################################################################
### Extract transcriptome sequences in FASTA from genome.
###############################################################################


rule extract_transcriptome_seqs:
    input:
        genome=INTERMEDIATES_DIR / "genome_processed.fa",
        gtf=config["gtf_file"],
    output:
        fasta=INTERMEDIATES_DIR / "transcriptome.fa",
    params:
        cluster_log=CLUSTER_LOG / "extract_transcriptome_seqs.log",
    log:
        LOCAL_LOG / "extract_transcriptome_seqs.log",
    container:
        "docker://quay.io/biocontainers/cufflinks:2.2.1--py27_2"
    conda:
        ENV_DIR / "cufflinks.yaml"
    shell:
        "(zcat {input.gtf} | gffread -w {output.fasta} -g {input.genome}) &> {log}"


###############################################################################
### Trim transcript IDs from FASTA file
###############################################################################


rule trim_transcriptome_seq_ids:
    input:
        fasta=INTERMEDIATES_DIR / "transcriptome.fa",
        script=SCRIPTS_DIR / "trim_id_fasta.sh",
    output:
        fasta=INTERMEDIATES_DIR / "transcriptome_trimmed_id.fa",
    params:
        cluster_log=CLUSTER_LOG / "trim_transcriptome.log",
    log:
        LOCAL_LOG / "trim_transcriptome.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.fasta} | {input.script} > {output.fasta}) &> {log}"


###############################################################################
### Generate segemehl index for transcripts
###############################################################################


rule generate_segemehl_index_transcriptome:
    input:
        fasta=INTERMEDIATES_DIR / "transcriptome_trimmed_id.fa",
    output:
        idx=INTERMEDIATES_DIR / "segemehl_transcriptome_index.idx",
    params:
        cluster_log=CLUSTER_LOG / "generate_segemehl_index_transcriptome.log",
    log:
        LOCAL_LOG / "generate_segemehl_index_transcriptome.log",
    resources:
        mem=10,
        threads=8,
        time=6,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    conda:
        ENV_DIR / "segemehl.yaml"
    shell:
        "(segemehl.x -x {output.idx} -d {input.fasta}) &> {log}"


###############################################################################
### Generate segemehl index for genome
###############################################################################


rule generate_segemehl_index_genome:
    input:
        genome=INTERMEDIATES_DIR / "genome_processed.fa",
    output:
        idx=INTERMEDIATES_DIR / "segemehl_genome_index.idx",
    params:
        cluster_log=CLUSTER_LOG / "generate_segemehl_index_genome.log",
    log:
        LOCAL_LOG / "generate_segemehl_index_genome.log",
    resources:
        mem=60,
        threads=8,
        time=6,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    conda:
        ENV_DIR / "segemehl.yaml"
    shell:
        "(segemehl.x -x {output.idx} -d {input.genome}) &> {log}"


###############################################################################
### GTF file of exons (genomic coordinates)
###############################################################################


rule get_exons_gtf:
    input:
        gtf=config["gtf_file"],
        script=SCRIPTS_DIR / "get_lines_w_pattern.sh",
    output:
        exons=INTERMEDIATES_DIR / "exons.gtf",
    params:
        cluster_log=CLUSTER_LOG / "get_exons_gtf.log",
    log:
        LOCAL_LOG / "get_exons_gtf.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(bash \
        {input.script} \
        -f {input.gtf} \
        -c 3 \
        -p exon \
        -o {output.exons} \
        ) &> {log}"


###############################################################################
### Convert GTF file of exons to BED file
###############################################################################


rule convert_exons_gtf_to_bed:
    input:
        exons=INTERMEDIATES_DIR / "exons.gtf",
        script=SCRIPTS_DIR / "gtf_exons_bed.1.1.2.R",
    output:
        exons=INTERMEDIATES_DIR / "exons.bed",
    params:
        cluster_log=CLUSTER_LOG / "exons_gtf_to_bed.log",
    log:
        LOCAL_LOG / "exons_gtf_to_bed.log",
    container:
        "docker://zavolab/r-zavolab:3.5.1"
    conda:
        ENV_DIR / "r.yaml"
    shell:
        "(Rscript \
        {input.script} \
        --gtf {input.exons} \
        -o {output.exons} \
        ) &> {log}"


###############################################################################
### Create header for SAM file
###############################################################################


rule create_genome_header:
    input:
        genome=INTERMEDIATES_DIR / "genome_processed.fa",
    output:
        header=INTERMEDIATES_DIR / "genome_header.sam",
    params:
        cluster_log=CLUSTER_LOG / "create_genome_header.log",
    log:
        LOCAL_LOG / "create_genome_header.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools dict -o {output.header} --uri=NA {input.genome}) &> {log}"


###############################################################################
### Mapping chromosomes names, UCSC <-> ENSEMBL
###############################################################################


rule map_chr_names:
    input:
        anno=config["mirna_file"],
        script=SCRIPTS_DIR / "map_chromosomes.pl",
        map_chr=config["map_chr_file"],
    output:
        gff=INTERMEDIATES_DIR / "mirna_annotations.gff3",
    params:
        cluster_log=CLUSTER_LOG / "map_chr_names.log",
        column="1",
        delimiter="TAB",
    log:
        LOCAL_LOG / "map_chr_names.log",
    container:
        "docker://perl:5.37.10"
    conda:
        ENV_DIR / "perl.yaml"
    shell:
        "(perl {input.script} \
        {input.anno} \
        {params.column} \
        {params.delimiter} \
        {input.map_chr} \
        {output.gff} \
        ) &> {log}"


###############################################################################
### Index genome fasta file
###############################################################################


rule create_index_genome_fasta:
    input:
        genome=INTERMEDIATES_DIR / "genome_processed.fa",
    output:
        genome=INTERMEDIATES_DIR / "genome_processed.fa.fai",
    params:
        cluster_log=CLUSTER_LOG / "create_index_genome_fasta.log",
    log:
        LOCAL_LOG / "create_index_genome_fasta.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools faidx {input.genome}) &> {log}"


###############################################################################
### Extract chromosome length
###############################################################################


rule extract_chr_len:
    input:
        genome=INTERMEDIATES_DIR / "genome_processed.fa.fai",
    output:
        chrsize=INTERMEDIATES_DIR / "chr_size.tsv",
    params:
        cluster_log=CLUSTER_LOG / "extract_chr_len.log",
    log:
        LOCAL_LOG / "extract_chr_len.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cut -f1,2 {input.genome} > {output.chrsize}) &> {log}"


###############################################################################
### Extend miRNAs annotations
###############################################################################


rule extend_mirs_annotations:
    input:
        gff3=INTERMEDIATES_DIR / "mirna_annotations.gff3",
        chrsize=INTERMEDIATES_DIR / "chr_size.tsv",
        script=SCRIPTS_DIR / "mirna_extension.py",
    output:
        extended_mir=expand(
            INTERMEDIATES_DIR / "extended_mirna_annotation_{extension}_nt.gff3",
            extension=config["extension"],
        ),
        extended_primir=expand(
            INTERMEDIATES_DIR / "extended_primir_annotation_{extension}_nt.gff3",
            extension=config["extension"],
        ),
    params:
        cluster_log=CLUSTER_LOG / "extend_mirs_annotations.log",
        out_dir=INTERMEDIATES_DIR,
        extension=config["extension"],
    log:
        LOCAL_LOG / "extend_mirs_annotations.log",
    container:
        "docker://quay.io/biocontainers/gffutils:0.11.1--pyh7cba7a3_0"
    conda:
        ENV_DIR / "gffutils.yaml"
    shell:
        "(python {input.script} \
        {input.gff3} \
        --chr {input.chrsize} \
        --extension {params.extension} \
        --outdir {params.out_dir} \
        ) &> {log}"
