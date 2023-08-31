###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Snakemake workflow to download and prepare the necessary files
# for smallRNA-seq related workflows.
#
###############################################################################

import os

from pathlib import Path

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
        idx_transcriptome=os.path.join(
            config["output_dir"],
            "segemehl_transcriptome_index.idx",
        ),
        idx_genome=os.path.join(
            config["output_dir"],
            "segemehl_genome_index.idx",
        ),
        exons=os.path.join(
            config["output_dir"],
            "exons.bed",
        ),
        header=os.path.join(
            config["output_dir"],
            "genome_header.sam",
        ),
        chrsize=os.path.join(config["output_dir"], "chr_size.txt"),
        extended_mir=expand(
            os.path.join(
                config["output_dir"],
                "extended_mirna_annotation_{extension}_nt.gff3",
            ),
            extension=config["extension"],
        ),
        extended_primir=expand(
            os.path.join(
                config["output_dir"],
                "extended_primir_annotation_{extension}_nt.gff3",
            ),
            extension=config["extension"],
        ),


###############################################################################
### Trim genome IDs
###############################################################################


rule trim_genome_seq_id:
    input:
        genome=config["genome_file"],
        script=os.path.join(config["scripts_dir"], "trim_id_fasta.sh"),
    output:
        genome=os.path.join(config["output_dir"], "genome_processed.fa"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "genome_process.log",
        ),
    log:
        os.path.join(config["local_log"], "genome_process.log"),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(zcat {input.genome} | {input.script} > {output.genome}) &> {log}"


###############################################################################
### Extract transcriptome sequences in FASTA from genome.
###############################################################################


rule extract_transcriptome_seqs:
    input:
        genome=os.path.join(config["output_dir"], "genome_processed.fa"),
        gtf=config["gtf_file"],
    output:
        fasta=os.path.join(config["output_dir"], "transcriptome.fa"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "extract_transcriptome_seqs.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "extract_transcriptome_seqs.log",
        ),
    container:
        "docker://quay.io/biocontainers/cufflinks:2.2.1--py27_2"
    shell:
        "(zcat {input.gtf} | gffread -w {output.fasta} -g {input.genome}) &> {log}"


###############################################################################
### Trim transcript IDs from FASTA file
###############################################################################


rule trim_fasta_seq_ids:
    input:
        fasta=os.path.join(config["output_dir"], "transcriptome.fa"),
        script=os.path.join(config["scripts_dir"], "trim_id_fasta.sh"),
    output:
        fasta=os.path.join(config["output_dir"], "transcriptome_trimmed_id.fa"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "trim_transcriptome.log"
        ),
    log:
        os.path.join(config["local_log"], "trim_transcriptome.log"),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.fasta} | {input.script} > {output.fasta}) &> {log}"


###############################################################################
### Generate segemehl index for transcripts
###############################################################################


rule generate_segemehl_index_transcriptome:
    input:
        fasta=os.path.join(config["output_dir"], "transcriptome_trimmed_id.fa"),
    output:
        idx=os.path.join(
            config["output_dir"],
            "segemehl_transcriptome_index.idx",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "generate_segemehl_index_transcriptome.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "generate_segemehl_index_transcriptome.log",
        ),
    resources:
        mem=10,
        threads=8,
        time=6,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    shell:
        "(segemehl.x -x {output.idx} -d {input.fasta}) &> {log}"


###############################################################################
### Generate segemehl index for genome
###############################################################################


rule generate_segemehl_index_genome:
    input:
        genome=os.path.join(config["output_dir"], "genome_processed.fa"),
    output:
        idx=os.path.join(config["output_dir"], "segemehl_genome_index.idx"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "generate_segemehl_index_genome.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "generate_segemehl_index_genome.log",
        ),
    resources:
        mem=60,
        threads=8,
        time=6,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    shell:
        "(segemehl.x -x {output.idx} -d {input.genome}) &> {log}"


###############################################################################
### GTF file of exons (genomic coordinates)
###############################################################################


rule get_exons_gtf:
    input:
        gtf=config["gtf_file"],
        script=os.path.join(config["scripts_dir"], "get_lines_w_pattern.sh"),
    output:
        exons=os.path.join(config["output_dir"], "exons.gtf"),
    params:
        cluster_log=os.path.join(config["cluster_log"], "get_exons_gtf.log"),
    log:
        os.path.join(config["local_log"], "get_exons_gtf.log"),
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


rule exons_gtf_to_bed:
    input:
        exons=os.path.join(config["output_dir"], "exons.gtf"),
        script=os.path.join(config["scripts_dir"], "gtf_exons_bed.1.1.2.R"),
    output:
        exons=os.path.join(config["output_dir"], "exons.bed"),
    params:
        cluster_log=os.path.join(config["cluster_log"], "exons_gtf_to_bed.log"),
    log:
        os.path.join(config["local_log"], "exons_gtf_to_bed.log"),
    container:
        "docker://zavolab/r-zavolab:3.5.1"
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
        genome=os.path.join(config["output_dir"], "genome_processed.fa"),
    output:
        header=os.path.join(config["output_dir"], "genome_header.sam"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "create_genome_header.log"
        ),
    log:
        os.path.join(config["local_log"], "create_genome_header.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools dict -o {output.header} --uri=NA {input.genome}) &> {log}"


###############################################################################
### Mapping chromosomes names, UCSC <-> ENSEMBL
###############################################################################


rule map_chr_names:
    input:
        anno=config["mirna_file"],
        script=os.path.join(config["scripts_dir"], "map_chromosomes.pl"),
        map_chr=config["map_chr_file"],
    output:
        gff=os.path.join(config["output_dir"], "mirna_annotations.gff3"),
    params:
        cluster_log=os.path.join(config["cluster_log"], "map_chr_names.log"),
        column="1",
        delimiter="TAB",
    log:
        os.path.join(config["local_log"], "map_chr_names.log"),
    container:
        "docker://perl:5.37.10"
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
        genome=os.path.join(config["output_dir"], "genome_processed.fa"),
    output:
        genome=os.path.join(config["output_dir"], "genome_processed.fa.fai"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "create_index_genome_fasta.log"
        ),
    log:
        os.path.join(config["local_log"], "create_index_genome_fasta.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools faidx {input.genome}) &> {log}"


###############################################################################
### Extract chromosome length
###############################################################################


rule extract_chr_len:
    input:
        genome=os.path.join(config["output_dir"], "genome_processed.fa.fai"),
    output:
        chrsize=os.path.join(config["output_dir"], "chr_size.txt"),
    params:
        cluster_log=os.path.join(config["cluster_log"], "extract_chr_len.log"),
    log:
        os.path.join(config["local_log"], "extract_chr_len.log"),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cut -f1,2 {input.genome} > {output.chrsize}) &> {log}"


###############################################################################
### Extend miRNAs annotations
###############################################################################


rule extend_mirs_annotations:
    input:
        gff3=os.path.join(config["output_dir"], "mirna_annotations.gff3"),
        chrsize=os.path.join(config["output_dir"], "chr_size.txt"),
        script=os.path.join(config["scripts_dir"], "mirna_extension.py"),
    output:
        extended_mir=expand(
            os.path.join(
                config["output_dir"],
                "extended_mirna_annotation_{extension}_nt.gff3",
            ),
            extension=config["extension"],
        ),
        extended_primir=expand(
            os.path.join(
                config["output_dir"],
                "extended_primir_annotation_{extension}_nt.gff3",
            ),
            extension=config["extension"],
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "extend_mirs_annotations.log"
        ),
        out_dir=lambda wildcards, input: Path(input[0]).parent,
        extension=config["extension"],
    log:
        os.path.join(config["local_log"], "extend_mirs_annotations.log"),
    container:
        "docker://quay.io/biocontainers/gffutils:0.11.1--pyh7cba7a3_0"
    shell:
        "(python {input.script} \
        {input.gff3} \
        --chr {input.chrsize} \
        --extension {params.extension} \
        --outdir {params.out_dir} \
        ) &> {log}"
