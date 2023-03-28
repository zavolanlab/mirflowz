###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Snakemake workflow to download and prepare the necessary files 
# for smallRNA-seq related workflows.
#
###############################################################################
#
# USAGE (from the file's directory):
#
# snakemake \
#    --snakefile="prepare.smk" \
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

localrules:
    finish_prepare,

###############################################################################
### Finish rule
###############################################################################
rule finish_prepare:
    input:
        idx_transcriptome=os.path.join(
                config["output_dir"],
                "transcriptome_index_segemehl.idx",

        ),
        idx_genome=os.path.join(
                config["output_dir"],
                "genome_index_segemehl.idx",
        ),
        exons=os.path.join(
                config["output_dir"], 
                "exons.bed",
        ),
        header=os.path.join(
                config["output_dir"],
                "headerOfCollapsedFasta.sam",
        ),
        chrsize=os.path.join(
            config["output_dir"], "chr_size.txt"
        ),
        mirnafilt=os.path.join(
                config["output_dir"], 
                "extended_mirna_filtered.bed",
        ),


###############################################################################
### Trim genome IDs
###############################################################################

rule trim_genome_seq_id:
    input:
        genome=config["genome_file"],
    output:
        genome=os.path.join(config["output_dir"], "genome.processed.fa"),
    params:
        dir_out=config["output_dir"],
        cluster_log=os.path.join(
            config["cluster_log"], "genome_process.log",
        ),
    log:
        os.path.join(config["local_log"], "genome_process.log"),
    singularity:
        "docker://ubuntu:bionic-20221215"
    shell:
        """(zcat {input.genome} | 
        awk \
        -F" " \
        "/^>/ {{print \$1; next}} 1" \
        > {output.genome} \
        ) &> {log}"""


###############################################################################
### Extract transcriptome sequences in FASTA from genome.
###############################################################################


rule extract_transcriptome_seqs:
    input:
        genome=os.path.join(
            config["output_dir"], "genome.processed.fa"
        ),
        gtf=config["gtf_file"],
    output:
        fasta=os.path.join(
            config["output_dir"], "transcriptome.fa"
        ),
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
    singularity:
        "docker://quay.io/biocontainers/cufflinks:2.2.1--py27_2"
    shell:
        "(zcat {input.gtf} | gffread -w {output.fasta} -g {input.genome}) &> {log}"


###############################################################################
### Trim transcript IDs from FASTA file
###############################################################################


rule trim_fasta:
    input:
        fasta=os.path.join(
            config["output_dir"], "transcriptome.fa"
        ),
        script=os.path.join(config["scripts_dir"], "validation_fasta.py"),
    output:
        fasta=os.path.join(
            config["output_dir"], "transcriptome_idtrim.fa"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "trim_fasta.log"
        ),
    log:
        os.path.join(config["local_log"], "trim_fasta.log"),
    singularity:
        "docker://ubuntu:bionic-20221215"
    shell:
        """(awk \
        -F" " \
        "/^>/ {{print \$1; next}} 1" \
        {input.fasta} \
        > {output.fasta} \
        ) &> {log}"""


###############################################################################
### Generate segemehl index for transcripts
###############################################################################


rule generate_segemehl_index_transcriptome:
    input:
        fasta=os.path.join(
            config["output_dir"], "transcriptome_idtrim.fa"
        ),
    output:
        idx=os.path.join(
            config["output_dir"],
            "transcriptome_index_segemehl.idx",
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
    singularity:
        "docker://quay.io/biocontainers/segemehl:0.2.0--hfb9b9cc_7"
    shell:
        "(segemehl.x -x {output.idx} -d {input.fasta}) &> {log}"


###############################################################################
### Generate segemehl index for genome
###############################################################################


rule generate_segemehl_index_genome:
    input:
        genome=os.path.join(
            config["output_dir"], "genome.processed.fa"
        ),
    output:
        idx=os.path.join(
            config["output_dir"], "genome_index_segemehl.idx"
        ),
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
    singularity:
        "docker://quay.io/biocontainers/segemehl:0.2.0--hfb9b9cc_7"
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
        cluster_log=os.path.join(
            config["cluster_log"], "get_exons_gtf.log"
        ),
    log:
        os.path.join(config["local_log"], "get_exons_gtf.log"),
    singularity:
        "docker://ubuntu:bionic-20221215"
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


rule gtftobed:
    input:
        exons=os.path.join(config["output_dir"], "exons.gtf"),
        script=os.path.join(config["scripts_dir"], "gtf_exons_bed.1.1.2.R"),
    output:
        exons=os.path.join(config["output_dir"], "exons.bed"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "gtftobed.log"
        ),
    log:
        os.path.join(config["local_log"], "gtftobed.log"),
    singularity:
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


rule create_header_genome:
    input:
        genome=os.path.join(
            config["output_dir"], "genome.processed.fa"
        ),
    output:
        header=os.path.join(
            config["output_dir"], "headerOfCollapsedFasta.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "create_header_genome.log"
        ),
    log:
        os.path.join(
            config["local_log"], "create_header_genome.log"
        ),
    singularity:
        "docker://quay.io/biocontainers/samtools:1.8--2"
    shell:
        "(samtools dict -o {output.header} --uri=NA {input.genome}) &> {log}"

###############################################################################
### Index genome fasta file
###############################################################################


rule create_index_fasta:
    input:
        genome=os.path.join(
            config["output_dir"], "genome.processed.fa"
        ),
    output:
        genome=os.path.join(
            config["output_dir"], "genome.processed.fa.fai"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "create_index_fasta.log"
        ),
    log:
        os.path.join(
            config["local_log"], "create_index_fasta.log"
        ),
    singularity:
        "docker://quay.io/biocontainers/samtools:1.8--2"
    shell:
        "(samtools faidx {input.genome}) &> {log}"


###############################################################################
### Extract chromosome length
###############################################################################


rule extract_chr_len:
    input:
        genome=os.path.join(
            config["output_dir"], "genome.processed.fa.fai"
        ),
    output:
        chrsize=os.path.join(
            config["output_dir"], "chr_size.txt"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "extract_chr_len.log"
        ),
    log:
        os.path.join(config["local_log"], "extract_chr_len.log"),
    singularity:
        "docker://ubuntu:bionic-20221215"
    shell:
        "(cut -f1,2 {input.genome} > {output.chrsize}) &> {log}"


###############################################################################
### Mapping chromosomes names, UCSC <-> ENSEMBL
###############################################################################


rule map_chr_names:
    input:
        anno=config["mirna_file"],
        script=os.path.join(config["scripts_dir"], "map_chromosomes.pl"),
        map_chr=config["map_chr_file"],
    output:
        gff=os.path.join(
            config["output_dir"], "mirna_annotations.gff3"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "map_chr_names.log"
        ),
        column="1",
        delimiter="TAB",
    log:
        os.path.join(config["local_log"], "map_chr_names.log"),
    singularity:
        "docker://perl:5.28"
    shell:
        "(perl {input.script} \
        {input.anno} \
        {params.column} \
        {params.delimiter} \
        {input.map_chr} \
        {output.gff} \
        ) &> {log}"


###############################################################################
### Etending pre-miRNA overhang
###############################################################################

rule extend_premir:
    input:
        gff=os.path.join(
            config["output_dir"], "mirna_annotations.gff3"
        ),
        script=os.path.join(
            config["scripts_dir"], "extend_overhang.py"
        ),
        chrsize=os.path.join(
            config["output_dir"], "chr_size.txt"
        ),
    output:
        gff=os.path.join(
            config["output_dir"], "extended_mirna_annotations.gff3"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "extend_overhang.log"
        ),
        extension=config["overhang"],
    log:
        os.path.join(config["local_log"], "extend_overhang.log"),
    singularity:
        "docker://python:3.9.16"
    shell:
        "(python {input.script} \
        -i {input.gff} \
        -e {params.extension} \
        --chr {input.chrsize} \
        -o {output.gff} \
        )&>{log}"


###############################################################################
### GFF to BED (improve intersect memory efficient allowing to use -sorted)
###############################################################################


rule gfftobed:
    input:
        gff=os.path.join(
            config["output_dir"], "extended_mirna_annotations.gff3"
        ),
    output:
        bed=os.path.join(
            config["output_dir"], "extended_mirna_annotations.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "gfftobed2.log"
        ),
        out_dir=config["output_dir"]
    log:
        os.path.join(config["local_log"], "gfftobed2.log"),
    singularity:
        "docker://quay.io/biocontainers/bedops:2.4.35--h6bb024c_2"
    shell:
        "(convert2bed -i gff < {input.gff} \
        --sort-tmpdir={params.out_dir} \
        > {output.bed} \
        ) &> {log}"


###############################################################################
### Extract mature miRNA
###############################################################################


rule filter_mature_mirs:
    input:
        bed=os.path.join(
            config["output_dir"], "extended_mirna_annotations.bed"
        ),
    output:
        bed=os.path.join(
            config["output_dir"], "mirna_mature_filtered.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "filter_mature_mirs.log"
        ),
        precursor="miRNA_primary_transcript",
    log:
        os.path.join(
            config["local_log"], "filter_mature_mirs.log"
        ),
    singularity:
        "docker://ubuntu:bionic-20221215"
    shell:
        "(grep -v {params.precursor} {input.bed} > {output.bed}) &> {log}"