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
        mirnafilt=os.path.join(
                config["output_dir"], 
                "mirna_filtered.bed",
        ),
        isomirs=os.path.join(
                config["output_dir"], 
                "isomirs_annotation.bed",
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
        "docker://ubuntu:focal-20210416"
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
        "docker://ubuntu:focal-20210416"
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
        "docker://ubuntu:focal-20210416"
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
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
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
        gff=os.path.join(
            config["output_dir"], "mirna_chr_mapped.gff3"
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
        "docker://quay.io/biocontainers/perl:5.26.2"
    shell:
        "(perl {input.script} \
        {input.anno} \
        {params.column} \
        {params.delimiter} \
        {input.map_chr} \
        {output.gff} \
        ) &> {log}"


###############################################################################
### Filtering _1 miR IDs
###############################################################################


rule filter_mir_1_anno:
    input:
        gff=os.path.join(
            config["output_dir"], "mirna_chr_mapped.gff3"
        ),
    output:
        gff=os.path.join(
            config["output_dir"], "mirna_filtered.gff3"
        ),
    params:
        script=os.path.join(config["scripts_dir"], "filter_mir_1_anno.sh"),
        cluster_log=os.path.join(
            config["cluster_log"], "filter_mir_1_anno.log"
        ),
    log:
        os.path.join(
            config["local_log"], "filter_mir_1_anno.log"
        ),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(bash {params.script} -f {input.gff} -o {output.gff}) &> {log}"


###############################################################################
### GFF to BED (improve intersect memory efficient allowing to use -sorted)
###############################################################################


rule gfftobed:
    input:
        gff=os.path.join(
            config["output_dir"], "mirna_filtered.gff3"
        ),
    output:
        bed=os.path.join(
            config["output_dir"], "mirna_filtered.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "gfftobed.log"
        ),
        out_dir=config["output_dir"]
    log:
        os.path.join(config["local_log"], "gfftobed.log"),
    singularity:
        "docker://quay.io/biocontainers/bedops:2.4.35--h6bb024c_2"
    shell:
        "(convert2bed -i gff < {input.gff} \
        --sort-tmpdir={params.out_dir} \
        > {output.bed} \
        ) &> {log}"


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
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
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
        "docker://ubuntu:focal-20210416"
    shell:
        "(cut -f1,2 {input.genome} > {output.chrsize}) &> {log}"


###############################################################################
### Extract mature miRNA
###############################################################################


rule filter_mature_mirs:
    input:
        bed=os.path.join(
            config["output_dir"], "mirna_filtered.bed"
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
        "docker://ubuntu:focal-20210416"
    shell:
        "(grep -v {params.precursor} {input.bed} > {output.bed}) &> {log}"


###############################################################################
### Create isomirs annotation file from mature miRNA 
###############################################################################


rule iso_anno:
    input:
        bed=os.path.join(
            config["output_dir"], "mirna_mature_filtered.bed"
        ),
        chrsize=os.path.join(
            config["output_dir"], "chr_size.txt"
        ),
    output:
        bed=os.path.join(
            config["output_dir"],
            "iso_anno_5p{bp_5p}_3p{bp_3p}.bed",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "iso_anno_5p{bp_5p}_3p{bp_3p}.log",
        ),
        bp_5p=lambda wildcards: wildcards.bp_5p,
        bp_3p=lambda wildcards: wildcards.bp_3p,
    log:
        os.path.join(
            config["local_log"],
            "iso_anno_5p{bp_5p}_3p{bp_3p}.log",
        ),
    singularity:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        "(bedtools slop \
        -i {input.bed} \
        -g {input.chrsize} \
        -l {params.bp_5p} \
        -r {params.bp_3p} \
        > {output.bed} \
        ) &> {log}"


###############################################################################
### Change miRNA names to isomirs names
###############################################################################


rule iso_anno_rename:
    input:
        bed=os.path.join(
            config["output_dir"],
            "iso_anno_5p{bp_5p}_3p{bp_3p}.bed",
        ),
    output:
        bed=os.path.join(
            config["output_dir"],
            "iso_anno_rename_5p{bp_5p}_3p{bp_3p}.bed",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "iso_anno_rename_5p{bp_5p}_3p{bp_3p}.log",
        ),
        bp_5p=lambda wildcards: wildcards.bp_5p,
        bp_3p=lambda wildcards: wildcards.bp_3p,
    log:
        os.path.join(
            config["local_log"],
            "iso_anno_rename_5p{bp_5p}_3p{bp_3p}.log",
        ),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(sed \
        's/;Derives/_5p{params.bp_5p}_3p{params.bp_3p};Derives/' \
        {input.bed} \
        > {output.bed} \
        ) &> {log}"


###############################################################################
### Concatenate all isomirs annotation files
###############################################################################


rule iso_anno_concat:
    input:
        bed=lambda wildcards: expand(
            os.path.join(
                config["output_dir"],
                "iso_anno_rename_5p{bp_5p}_3p{bp_3p}.bed",
            ),
            bp_3p=config["bp_3p"],
            bp_5p=config["bp_5p"],
        ),
    output:
        bed=os.path.join(
            config["output_dir"], "iso_anno_concat.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "iso_anno_concat.log"
        ),
        prefix=os.path.join(
            config["output_dir"], "iso_anno_rename"
        ),
    log:
        os.path.join(config["local_log"], "iso_anno_concat.log"),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(cat {params.prefix}* > {output.bed}) &> {log}"


###############################################################################
### Remove non changing isomirs (5p0_3p0)
###############################################################################


rule iso_anno_final:
    input:
        bed=os.path.join(
            config["output_dir"], "iso_anno_concat.bed"
        ),
    output:
        bed=os.path.join(
            config["output_dir"], "isomirs_annotation.bed"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "iso_anno_final.log"
        ),
        pattern="5p0_3p0",
    log:
        os.path.join(config["local_log"], "iso_anno_final.log"),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(grep -v '{params.pattern}' {input.bed} > {output.bed}) &> {log}"