###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Workflow to map small RNA-seq reads (e.g. from miRNA sequencing libraries).
###############################################################################
#
# USAGE (from the file's directory):
#
# snakemake \
#    --snakefile="map.smk" \
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
### Reading sample and resources tables
###############################################################################

samples_table = pd.read_csv(
    config["samples"],
    header = 0,
    index_col = 0,
    comment = "#",
    engine = "python",
    sep = "\t",
)

###############################################################################
### Funcitons get_sample and get_resource
###############################################################################
# Function to get relevant per sample information from samples and resources table

def get_sample(column_id, sample_id = None):
    if sample_id:
        return str(
            samples_table[column_id][samples_table.index == sample_id][0]
        )
    else:
        return str(
            samples_table[column_id][0]
        )

###############################################################################
### Global configuration
###############################################################################
# Rules that require internet connection for downloading files are included
# in the localrules
localrules:
    start,
    finish_map,

###############################################################################
### Finish rule
###############################################################################


rule finish_map:
    input:
        maps=expand(
            os.path.join(
                config["output_dir"],
                "{sample}",
                "convertedSortedMappings_{sample}.bam.bai",
            ),
            sample=pd.unique(samples_table.index.values),
        ),



###############################################################################
### Start rule (get samples)
###############################################################################

rule start:
    input:
        reads=lambda wildcards: expand(
            pd.Series(samples_table.loc[wildcards.sample, "sample_file"]).values,
            format=get_sample("format"),
        ),
    output:
        reads=os.path.join(
            config["output_dir"],
            "{sample}",
            "{format}",
            "reads.{format}",
        )
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "uncompress_zipped_files_{sample}_{format}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "uncompress_zipped_files_{sample}_{format}.log",
        ),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(zcat {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Quality filter
###############################################################################


rule fastq_quality_filter:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "fastq", "reads.fastq"
        ),
    output:
        reads=os.path.join(
            config["output_dir"], "{sample}", "fastq", "filtered_reads.fastq"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "fastq_quality_filter_{sample}.log"
        ),
        p=config["p_value"],
        q=config["q_value"],
    log:
        os.path.join(config["local_log"], "fastq_quality_filter_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--2"
    shell:
        "(fastq_quality_filter \
        -v \
        -q {params.q} \
        -p {params.p} \
        -i {input.reads} \
        > {output.reads} \
        ) &> {log}"


###############################################################################
### Convert fastq to fasta
###############################################################################


rule fastq_to_fasta:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "fastq", "filtered_reads.fastq"
        ),
    output:
        reads=os.path.join(
            config["output_dir"], "{sample}", "fastq", "reads.fa"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "fastq_to_fasta_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "fastq_to_fasta_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--2"
    shell:
        "(fastq_to_fasta -r -n -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Format fasta file
###############################################################################


rule fasta_formatter:
    input:
        reads=lambda wildcards: os.path.join(
            config["output_dir"],
            wildcards.sample,
            get_sample("format", wildcards.sample),
            "reads.fa",
        ),
    output:
        reads=os.path.join(config["output_dir"], "{sample}", "formatted.fasta"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "fasta_formatter_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "fasta_formatter_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--2"
    shell:
        "(fasta_formatter -w 0 -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Remove adapters
###############################################################################


rule cutadapt:
    input:
        reads=os.path.join(config["output_dir"], "{sample}", "formatted.fasta"),
    output:
        reads=os.path.join(config["output_dir"], "{sample}", "cut.fasta"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "cutadapt_{sample}.log"
        ),
        adapter=lambda wildcards: get_sample("adapter", wildcards.sample),
        error_rate=config["error_rate"],
        minimum_length=config["minimum_length"],
        overlap=config["overlap"],
        max_n=config["max_n"],
    log:
        os.path.join(config["local_log"], "cutadapt_{sample}.log"),
    resources:
        threads=8,
    singularity:
        "docker://quay.io/biocontainers/cutadapt:1.16--py35_2"
    shell:
        "(cutadapt \
        -a {params.adapter} \
        --error-rate {params.error_rate} \
        --minimum-length {params.minimum_length} \
        --overlap {params.overlap} \
        --trim-n \
        --max-n {params.max_n} \
        --cores {resources.threads} \
        -o {output.reads} {input.reads}) &> {log}"


###############################################################################
### Collapse identical reads
###############################################################################


rule fastx_collapser:
    input:
        reads=os.path.join(config["output_dir"], "{sample}", "cut.fasta"),
    output:
        reads=os.path.join(config["output_dir"], "{sample}", "collapsed.fasta"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "fastx_collapser_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "fastx_collapser_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--2"
    shell:
        "(fastx_collapser -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Segemehl genome mapping
###############################################################################


rule mapping_genome_segemehl:
    input:
        reads=os.path.join(config["output_dir"], "{sample}", "collapsed.fasta"),
        genome=os.path.join(config["output_dir"], "genome.processed.fa"),
        genome_index_segemehl=os.path.join(config["output_dir"], "genome_index_segemehl.idx"),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "segemehlGenome_map.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "mapping_genome_segemehl_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "mapping_genome_segemehl_{sample}.log"
        ),
    resources:
        mem=50,
        time=12,
        threads=8,
    singularity:
        "docker://quay.io/biocontainers/segemehl:0.2.0--hfb9b9cc_7"
    shell:
        "(segemehl.x \
        -i {input.genome_index_segemehl} \
        -d {input.genome} \
        -t {threads} \
        -q {input.reads} \
        -outfile {output.gmap} \
        ) &> {log}"


###############################################################################
### Segemehl transcriptome mapping
###############################################################################


rule mapping_transcriptome_segemehl:
    input:
        reads=os.path.join(config["output_dir"], "{sample}", "collapsed.fasta"),
        transcriptome=os.path.join(config["output_dir"], "transcriptome_idtrim.fa"),
        transcriptome_index_segemehl=os.path.join(config["output_dir"], "transcriptome_index_segemehl.idx"),
    output:
        tmap=os.path.join(
            config["output_dir"], "{sample}", "segemehlTranscriptome_map.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "mapping_transcriptome_segemehl_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"], "mapping_transcriptome_segemehl_{sample}.log"
        ),
    resources:
        mem=10,
        time=12,
        threads=8,
    singularity:
        "docker://quay.io/biocontainers/segemehl:0.2.0--hfb9b9cc_7"
    shell:
        "(segemehl.x \
        -i {input.transcriptome_index_segemehl} \
        -d {input.transcriptome} \
        -t {threads} \
        -q {input.reads} \
        -outfile {output.tmap} \
        ) &> {log}"


###############################################################################
### Filter fasta for oligomap mapping
###############################################################################


rule filter_fasta_for_oligomap:
    input:
        reads=os.path.join(config["output_dir"], "{sample}", "collapsed.fasta"),
        script=os.path.join(config["scripts_dir"], "validation_fasta.py"),
    output:
        reads=os.path.join(
            config["output_dir"], "{sample}", "filtered_for_oligomap.fasta"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "filter_fasta_for_oligomap_{sample}.log"
        ),
        max_length_reads=config["max_length_reads"],
    log:
        os.path.join(
            config["local_log"], "filter_fasta_for_oligomap_{sample}.log"
        ),
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        -r {params.max_length_reads} \
        -i {input.reads} \
        -o {output.reads} \
        ) &> {log}"


###############################################################################
### Oligomap genome mapping
###############################################################################


rule mapping_genome_oligomap:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "filtered_for_oligomap.fasta"
        ),
        target=os.path.join(config["output_dir"], "genome.processed.fa"),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_map.fa"
        ),
        report=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_report.txt"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "mapping_genome_oligomap_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "mapping_genome_oligomap_{sample}.log"
        ),
    resources:
        mem=50,
        time=6,
        threads=8,
    singularity:
        "docker://zavolab/oligomap:1.0"
    shell:
        "(oligomap \
        {input.target} \
        {input.reads} \
        -r {output.report} \
        > {output.gmap} \
        ) &> {log}"


###############################################################################
### Oligomap genome sorting
###############################################################################


rule sort_genome_oligomap:
    input:
        tmap=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_map.fa"
        ),
        report=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_report.txt"
        ),
        script=os.path.join(config["scripts_dir"], "blocksort.sh"),
    output:
        sort=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_sorted.fa"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sorting_genome_oligomap_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "sorting_genome_oligomap_{sample}.log"
        ),
    resources:
        threads=8,
        time=6,
    shell:
        "(bash {input.script} \
        {input.tmap} \
        {resources.threads} \
        {output.sort} \
        ) &> {log}"


###############################################################################
### Oligomap genome mapping output to SAM
###############################################################################


rule oligomap_genome_toSAM:
    input:
        report=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_report.txt"
        ),
        sort=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_sorted.fa"
        ),
        script=os.path.join(
            config["scripts_dir"], "oligomapOutputToSam_nhfiltered.py"
        ),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_converted.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "oligomap_genome_toSAM_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(config["local_log"], "oligomap_genome_toSAM_{sample}.log"),
    resources:
        time=1,
        queue=1,
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        -i {input.sort} \
        -n {params.nh} \
        > {output.gmap}) &> {log}"


###############################################################################
### Oligomap trancriptome mapping
###############################################################################


rule mapping_transcriptome_oligomap:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "filtered_for_oligomap.fasta"
        ),
        target=os.path.join(config["output_dir"], "transcriptome_idtrim.fa"),
    output:
        tmap=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_map.fa"
        ),
        report=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_report.txt"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "mapping_transcriptome_oligomap_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"], "mapping_transcriptome_oligomap_{sample}.log"
        ),
    resources:
        mem=10,
        time=6,
        threads=8,
    singularity:
        "docker://zavolab/oligomap:1.0"
    shell:
        "(oligomap \
        {input.target} \
        {input.reads} \
        -s \
        -r {output.report} \
        > {output.tmap} \
        ) &> {log}"


###############################################################################
### Oligomap trancriptome sorting
###############################################################################


rule sort_transcriptome_oligomap:
    input:
        tmap=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_map.fa"
        ),
        report=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_report.txt"
        ),
        script=os.path.join(config["scripts_dir"], "blocksort.sh"),
    output:
        sort=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_sorted.fa"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "sorting_transcriptome_oligomap_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"], "sorting_transcriptome_oligomap_{sample}.log"
        ),
    resources:
        threads=8,
    shell:
        "(bash {input.script} \
        {input.tmap} \
        {resources.threads} \
        {output.sort} \
        ) &> {log}"


###############################################################################
### Oligomap transcriptome mapping ouput to SAM
###############################################################################


rule oligomap_transcriptome_toSAM:
    input:
        report=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_report.txt"
        ),
        sort=os.path.join(
            config["output_dir"], "{sample}", "oligoTranscriptome_sorted.fa"
        ),
        script=os.path.join(
            config["scripts_dir"], "oligomapOutputToSam_nhfiltered.py"
        ),
    output:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligoTranscriptome_converted.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "oligomap_transcriptome_toSAM_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(
            config["local_log"], "oligomap_transcriptome_toSAM_{sample}.log"
        ),
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        -i {input.sort} \
        -n {params.nh} \
        > {output.tmap} \
        ) &> {log}"


###############################################################################
### Merge genome mappings
###############################################################################


rule merge_genome_maps:
    input:
        gmap1=os.path.join(
            config["output_dir"], "{sample}", "segemehlGenome_map.sam"
        ),
        gmap2=os.path.join(
            config["output_dir"], "{sample}", "oligoGenome_converted.sam"
        ),
    output:
        gmaps=os.path.join(
            config["output_dir"], "{sample}", "GenomeMappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "merge_genome_maps_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "merge_genome_maps_{sample}.log"),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(cat {input.gmap1} {input.gmap2} > {output.gmaps}) &> {log}"


###############################################################################
### Merge trancriptome mappings
###############################################################################


rule merge_transcriptome_maps:
    input:
        tmap1=os.path.join(
            config["output_dir"], "{sample}", "segemehlTranscriptome_map.sam"
        ),
        tmap2=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligoTranscriptome_converted.sam",
        ),
    output:
        tmaps=os.path.join(
            config["output_dir"], "{sample}", "TranscriptomeMappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "merge_transcriptome_maps_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "merge_transcriptome_maps_{sample}.log"
        ),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(cat {input.tmap1} {input.tmap2} > {output.tmaps}) &> {log}"


###############################################################################
### Filter NH genome
###############################################################################


rule nh_filter_genome:
    input:
        gmaps=os.path.join(
            config["output_dir"], "{sample}", "GenomeMappings.sam"
        ),
        script=os.path.join(config["scripts_dir"], "nh_filter.py"),
    output:
        gmaps=os.path.join(
            config["output_dir"], "{sample}", "nhfiltered_GenomeMappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "nh_filter_genome_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(config["local_log"], "nh_filter_genome_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        {input.gmaps} \
        {params.nh} \
        {output.gmaps} \
        ) &> {log}"


###############################################################################
### Filter NH transcriptome
###############################################################################


rule filter_nh_transcriptome:
    input:
        tmaps=os.path.join(
            config["output_dir"], "{sample}", "TranscriptomeMappings.sam"
        ),
        script=os.path.join(config["scripts_dir"], "nh_filter.py"),
    output:
        tmaps=os.path.join(
            config["output_dir"],
            "{sample}",
            "nhfiltered_TranscriptomeMappings.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "filter_nh_transcriptome_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(
            config["local_log"], "filter_nh_transcriptome_{sample}.log"
        ),
    singularity:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        {input.tmaps} \
        {params.nh} \
        {output.tmaps} \
        ) &> {log}"


###############################################################################
### Remove header genome mappings
###############################################################################


rule remove_headers_genome:
    input:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "nhfiltered_GenomeMappings.sam"
        ),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "noheader_GenomeMappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "remove_headers_genome_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "remove_headers_genome_{sample}.log"),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "samtools view {input.gmap} > {output.gmap}"


###############################################################################
### Remove header transcriptome mappings
###############################################################################


rule remove_headers_transcriptome:
    input:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "nhfiltered_TranscriptomeMappings.sam",
        ),
    output:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "noheader_TranscriptomeMappings.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "remove_headers_transcriptome_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "remove_headers_transcriptome_{sample}.log"
        ),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "samtools view {input.tmap} > {output.tmap}"


###############################################################################
### Transcriptome to genome coordinates
###############################################################################


rule trans_to_gen:
    input:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "noheader_TranscriptomeMappings.sam",
        ),
        script=os.path.join(config["scripts_dir"], "sam_trx_to_sam_gen.pl"),
        exons=os.path.join(config["output_dir"], "exons.bed"),
    output:
        genout=os.path.join(config["output_dir"], "{sample}", "TransToGen.sam"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "trans_to_gen_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "trans_to_gen_{sample}.log"),
    singularity:
        "docker://quay.io/biocontainers/perl:5.26.2"
    shell:
        "(perl {input.script} \
        --in {input.tmap} \
        --exons {input.exons} \
        --out {output.genout} \
        ) &> {log}"


###############################################################################
### Concatenate genome and trancriptome mappings
###############################################################################


rule cat_mapping:
    input:
        gmap1=os.path.join(config["output_dir"], "{sample}", "TransToGen.sam"),
        gmap2=os.path.join(
            config["output_dir"], "{sample}", "noheader_GenomeMappings.sam"
        ),
    output:
        catmaps=os.path.join(
            config["output_dir"], "{sample}", "catMappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "cat_mapping_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "cat_mapping_{sample}.log"),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(cat {input.gmap1} {input.gmap2} > {output.catmaps}) &> {log}"


###############################################################################
### Add header
###############################################################################


rule add_header:
    input:
        header=os.path.join(config["output_dir"], "headerOfCollapsedFasta.sam"),
        catmaps=os.path.join(
            config["output_dir"], "{sample}", "catMappings.sam"
        ),
    output:
        concatenate=os.path.join(
            config["output_dir"],
            "{sample}",
            "concatenated_header_catMappings.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "add_header_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "add_header_{sample}.log"),
    singularity:
        "docker://ubuntu:focal-20210416"
    shell:
        "(cat {input.header} {input.catmaps} > {output.concatenate}) &> {log}"


###############################################################################
### Sort mapped file by IDs
###############################################################################


rule sort_id:
    input:
        concatenate=os.path.join(
            config["output_dir"],
            "{sample}",
            "concatenated_header_catMappings.sam",
        ),
    output:
        sort=os.path.join(
            config["output_dir"], "{sample}", "header_sorted_catMappings.sam"
        ),
    params:
        cluster_log=os.path.join(config["cluster_log"], "sort_id_{sample}.log"),
    log:
        os.path.join(config["local_log"], "sort_id_{sample}.log"),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "(samtools sort -n -o {output.sort} {input.concatenate}) &> {log}"


###############################################################################
### Remove inferior mappings (keeping multimappers)
###############################################################################


rule remove_inferiors:
    input:
        sort=os.path.join(
            config["output_dir"], "{sample}", "header_sorted_catMappings.sam"
        ),
        script=os.path.join(
            config["scripts_dir"],
            "sam_remove_duplicates_inferior_alignments_multimappers.1_5.pl",
        ),
    output:
        remove_inf=os.path.join(
            config["output_dir"], "{sample}", "removeInferiors.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "remove_inferiors_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "remove_inferiors_{sample}.log"),
    resources:
        mem=15,
        threads=4,
    singularity:
        "docker://quay.io/biocontainers/perl:5.26.2"
    shell:
        "(perl {input.script} \
        --print-header \
        --keep-mm \
        --in {input.sort} \
        --out {output.remove_inf} \
        ) &> {log}"


###############################################################################
### Uncollapse reads
###############################################################################


rule uncollapse_reads:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "removeInferiors.sam"
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
    singularity:
        "docker://quay.io/biocontainers/perl:5.26.2"
    shell:
        "(perl {input.script} \
        --suffix \
        --in {input.maps} \
        --out {output.maps} \
        ) &> {log}"


###############################################################################
### Convert SAM to BAM
###############################################################################


rule convert_to_bam:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "uncollapsedMappings.sam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"], "{sample}", "mappingsConverted.bam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "convert_to_bam_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "convert_to_bam_{sample}.log"),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_by_position:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "mappingsConverted.bam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "convertedSortedMappings_{sample}.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sort_by_position_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "sort_by_position_{sample}.log"),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "convertedSortedMappings_{sample}.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "convertedSortedMappings_{sample}.bam.bai",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "index_bam_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "index_bam_{sample}.log"),
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"
