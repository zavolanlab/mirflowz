###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Workflow to map small RNA-seq reads (e.g. from miRNA sequencing libraries).
###############################################################################

import os
import pandas as pd


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
                "alignments_all_sorted_{sample}.bam.bai",
            ),
            sample=pd.unique(samples_table.index.values),
        ),


###############################################################################
### Start rule (get samples)
###############################################################################


rule start:
    input:
        reads=lambda wildcards: expand(
            pd.Series(
                samples_table.loc[wildcards.sample, "sample_file"]
            ).values,
            format=get_sample("format"),
        ),
    output:
        reads=os.path.join(
            config["output_dir"],
            "{sample}",
            "{format}",
            "reads.{format}",
        ),
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
    container:
        "docker://ubuntu:lunar-20221207"
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
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--H87F3376_10"
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
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--H87F3376_10"
    shell:
        "(fastq_to_fasta -r -n -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Format fasta file
###############################################################################


rule format_fasta:
    input:
        reads=lambda wildcards: os.path.join(
            config["output_dir"],
            wildcards.sample,
            get_sample("format", wildcards.sample),
            "reads.fa",
        ),
    output:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_formatted.fasta"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "format_fasta_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "format_fasta_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--H87F3376_10"
    shell:
        "(fasta_formatter -w 0 -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Remove adapters
###############################################################################


rule remove_adapters:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_formatted.fasta"
        ),
    output:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_trimmed_adapters.fasta"
        ),
    params:
        adapter=lambda wildcards: get_sample("adapter", wildcards.sample),
        error_rate=config["error_rate"],
        minimum_length=config["minimum_length"],
        overlap=config["overlap"],
        max_n=config["max_n"],
        cluster_log=os.path.join(
            config["cluster_log"], "remove_adapters_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "remove_adapters_{sample}.log"),
    resources:
        threads=8,
    container:
        "docker://quay.io/biocontainers/cutadapt:4.3--py310h1425a21_0"
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


rule collapse_identical_reads:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_trimmed_adapters.fasta"
        ),
    output:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_collapsed.fasta"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "collapse__identical_reads_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "collapse_identical_reads_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--H87F3376_10"
    shell:
        "(fastx_collapser -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Segemehl genome mapping
###############################################################################


rule map_genome_segemehl:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_collapsed.fasta"
        ),
        genome=os.path.join(config["output_dir"], "genome_processed.fa"),
        genome_index_segemehl=os.path.join(
            config["output_dir"], "segemehl_genome_index.idx"
        ),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "segemehl_genome_mappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "map_genome_segemehl_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "map_genome_segemehl_{sample}.log"),
    resources:
        mem=50,
        time=12,
        threads=8,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    shell:
        "(segemehl.x \
        -i {input.genome_index_segemehl} \
        -d {input.genome} \
        -t {threads} \
        -e \
        -q {input.reads} \
        -outfile {output.gmap} \
        ) &> {log}"


###############################################################################
### Segemehl transcriptome mapping
###############################################################################


rule map_transcriptome_segemehl:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_collapsed.fasta"
        ),
        transcriptome=os.path.join(
            config["output_dir"], "transcriptome_trimmed_id.fa"
        ),
        transcriptome_index_segemehl=os.path.join(
            config["output_dir"], "segemehl_transcriptome_index.idx"
        ),
    output:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "segemehl_transcriptome_mappings.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "map_transcriptome_segemehl_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"], "map_transcriptome_segemehl_{sample}.log"
        ),
    resources:
        mem=10,
        time=12,
        threads=8,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    shell:
        "(segemehl.x \
        -i {input.transcriptome_index_segemehl} \
        -d {input.transcriptome} \
        -t {threads} \
        -q {input.reads} \
        -e \
        -outfile {output.tmap} \
        ) &> {log}"


###############################################################################
### Filter fasta for oligomap mapping
###############################################################################


rule filter_fasta_for_oligomap:
    input:
        reads=os.path.join(
            config["output_dir"], "{sample}", "reads_collapsed.fasta"
        ),
        script=os.path.join(config["scripts_dir"], "validation_fasta.py"),
    output:
        reads=os.path.join(
            config["output_dir"],
            "{sample}",
            "reads_filtered_for_oligomap.fasta",
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
    container:
        "docker://python:3.9.16"
    shell:
        "(python {input.script} \
        -r {params.max_length_reads} \
        -i {input.reads} \
        -o {output.reads} \
        ) &> {log}"


###############################################################################
### Oligomap genome mapping
###############################################################################


rule map_genome_oligomap:
    input:
        reads=os.path.join(
            config["output_dir"],
            "{sample}",
            "reads_filtered_for_oligomap.fasta",
        ),
        target=os.path.join(config["output_dir"], "genome_processed.fa"),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_mappings.fasta"
        ),
        report=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_report.txt"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "map_genome_oligomap_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "map_genome_oligomap_{sample}.log"),
    resources:
        mem=50,
        time=6,
        threads=8,
    container:
        "docker://quay.io/biocontainers/oligomap:1.0.1--hdcf5f25_0"
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
            config["output_dir"], "{sample}", "oligomap_genome_mappings.fasta"
        ),
        report=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_report.txt"
        ),
        script=os.path.join(config["scripts_dir"], "blocksort.sh"),
    output:
        sort=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_sorted.fasta"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sort_genome_oligomap_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "sort_genome_oligomap_{sample}.log"),
    resources:
        threads=8,
        time=6,
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(bash {input.script} \
        {input.tmap} \
        {resources.threads} \
        {output.sort} \
        ) &> {log}"


###############################################################################
### Oligomap genome mapping output to SAM
###############################################################################


rule oligomap_genome_to_sam:
    input:
        report=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_report.txt"
        ),
        sort=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_sorted.fasta"
        ),
        script=os.path.join(
            config["scripts_dir"], "oligomapOutputToSam_nhfiltered.py"
        ),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_mappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "oligomap_genome_to_sam_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(config["local_log"], "oligomap_genome_to_sam_{sample}.log"),
    resources:
        time=1,
        queue=1,
    container:
        "docker://python:3.9.16"
    shell:
        "(python {input.script} \
        -i {input.sort} \
        -n {params.nh} \
        > {output.gmap}) &> {log}"


###############################################################################
### Oligomap transcriptome mapping
###############################################################################


rule map_transcriptome_oligomap:
    input:
        reads=os.path.join(
            config["output_dir"],
            "{sample}",
            "reads_filtered_for_oligomap.fasta",
        ),
        target=os.path.join(config["output_dir"], "transcriptome_trimmed_id.fa"),
    output:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_mappings.fasta",
        ),
        report=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_report.txt",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "map_transcriptome_oligomap_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"], "map_transcriptome_oligomap_{sample}.log"
        ),
    resources:
        mem=10,
        time=6,
        threads=8,
    container:
        "docker://quay.io/biocontainers/oligomap:1.0.1--hdcf5f25_0"
    shell:
        "(oligomap \
        {input.target} \
        {input.reads} \
        -s \
        -r {output.report} \
        > {output.tmap} \
        ) &> {log}"


###############################################################################
### Oligomap transcriptome sorting
###############################################################################


rule sort_transcriptome_oligomap:
    input:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_mappings.fasta",
        ),
        report=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_report.txt",
        ),
        script=os.path.join(config["scripts_dir"], "blocksort.sh"),
    output:
        sort=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_sorted.fasta",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "sort_transcriptome_oligomap_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"], "sort_transcriptome_oligomap_{sample}.log"
        ),
    resources:
        threads=8,
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(bash {input.script} \
        {input.tmap} \
        {resources.threads} \
        {output.sort} \
        ) &> {log}"


###############################################################################
### Oligomap transcriptome mapping output to SAM
###############################################################################


rule oligomap_transcriptome_to_sam:
    input:
        report=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_report.txt",
        ),
        sort=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_sorted.fasta",
        ),
        script=os.path.join(
            config["scripts_dir"], "oligomapOutputToSam_nhfiltered.py"
        ),
    output:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_mappings.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "oligomap_transcriptome_to_sam_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(
            config["local_log"], "oligomap_transcriptome_to_sam_{sample}.log"
        ),
    container:
        "docker://python:3.9.16"
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
            config["output_dir"], "{sample}", "segemehl_genome_mappings.sam"
        ),
        gmap2=os.path.join(
            config["output_dir"], "{sample}", "oligomap_genome_mappings.sam"
        ),
    output:
        gmaps=os.path.join(
            config["output_dir"], "{sample}", "genome_mappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "merge_genome_maps_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "merge_genome_maps_{sample}.log"),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.gmap1} {input.gmap2} > {output.gmaps}) &> {log}"


###############################################################################
### Merge transcriptome mappings
###############################################################################


rule merge_transcriptome_maps:
    input:
        tmap1=os.path.join(
            config["output_dir"],
            "{sample}",
            "segemehl_transcriptome_mappings.sam",
        ),
        tmap2=os.path.join(
            config["output_dir"],
            "{sample}",
            "oligomap_transcriptome_mappings.sam",
        ),
    output:
        tmaps=os.path.join(
            config["output_dir"], "{sample}", "transcriptome_mappings.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "merge_transcriptome_maps_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "merge_transcriptome_maps_{sample}.log"
        ),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.tmap1} {input.tmap2} > {output.tmaps}) &> {log}"


###############################################################################
### Filter NH genome
###############################################################################


rule filter_genome_by_nh:
    input:
        gmaps=os.path.join(
            config["output_dir"], "{sample}", "genome_mappings.sam"
        ),
        script=os.path.join(config["scripts_dir"], "nh_filter.py"),
    output:
        gmaps=os.path.join(
            config["output_dir"], "{sample}", "genome_mappings_filtered_nh.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "filter_genome_by_nh_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(config["local_log"], "filter_genome_by_nh_{sample}.log"),
    container:
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


rule filter_transcriptome_by_nh:
    input:
        tmaps=os.path.join(
            config["output_dir"], "{sample}", "transcriptome_mappings.sam"
        ),
        script=os.path.join(config["scripts_dir"], "nh_filter.py"),
    output:
        tmaps=os.path.join(
            config["output_dir"],
            "{sample}",
            "transcriptome_mappings_filtered_nh.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "filter_transcriptome_by_nh_{sample}.log"
        ),
        nh=config["nh"],
    log:
        os.path.join(
            config["local_log"], "filter_transcriptome_by_nh_{sample}.log"
        ),
    container:
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


rule remove_header_genome_mappings:
    input:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "genome_mappings_filtered_nh.sam"
        ),
    output:
        gmap=os.path.join(
            config["output_dir"], "{sample}", "genome_mappings_no_header.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "remove_header_genome_mappings_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "remove_header_genome_mappings_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "samtools view {input.gmap} > {output.gmap}"


###############################################################################
### Remove header transcriptome mappings
###############################################################################


rule remove_header_transcriptome_mappings:
    input:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "transcriptome_mappings_filtered_nh.sam",
        ),
    output:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "transcriptome_mappings_no_header.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"],
            "remove_header_transcriptome_mappings_{sample}.log",
        ),
    log:
        os.path.join(
            config["local_log"],
            "remove_header_transcriptome_mappings_{sample}.log",
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "samtools view {input.tmap} > {output.tmap}"


###############################################################################
### Transcriptome to genome coordinates
###############################################################################


rule transcriptome_to_genome_maps:
    input:
        tmap=os.path.join(
            config["output_dir"],
            "{sample}",
            "transcriptome_mappings_no_header.sam",
        ),
        script=os.path.join(config["scripts_dir"], "sam_trx_to_sam_gen.pl"),
        exons=os.path.join(config["output_dir"], "exons.bed"),
    output:
        genout=os.path.join(
            config["output_dir"],
            "{sample}",
            "transcriptome_mappings_to_genome.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "transcriptome_to_genome_maps_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "transcriptome_to_genome_maps_{sample}.log"
        ),
    container:
        "docker://perl:5.37.10"
    shell:
        "(perl {input.script} \
        --in {input.tmap} \
        --exons {input.exons} \
        --out {output.genout} \
        ) &> {log}"


###############################################################################
### Concatenate genome and transcriptome mappings
###############################################################################


rule merge_all_maps:
    input:
        gmap1=os.path.join(
            config["output_dir"],
            "{sample}",
            "transcriptome_mappings_to_genome.sam",
        ),
        gmap2=os.path.join(
            config["output_dir"], "{sample}", "genome_mappings_no_header.sam"
        ),
    output:
        catmaps=os.path.join(
            config["output_dir"], "{sample}", "mappings_all_no_header.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "merge_all_mappings_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "merge_all_mappings_{sample}.log"),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.gmap1} {input.gmap2} > {output.catmaps}) &> {log}"


###############################################################################
### Add header
###############################################################################


rule add_header_all_maps:
    input:
        header=os.path.join(config["output_dir"], "genome_header.sam"),
        catmaps=os.path.join(
            config["output_dir"], "{sample}", "mappings_all_no_header.sam"
        ),
    output:
        concatenate=os.path.join(
            config["output_dir"],
            "{sample}",
            "mappings_all.sam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "add_header_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "add_header_{sample}.log"),
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.header} {input.catmaps} > {output.concatenate}) &> {log}"


###############################################################################
### Sort mapped file by IDs
###############################################################################


rule sort_maps_by_id:
    input:
        concatenate=os.path.join(
            config["output_dir"],
            "{sample}",
            "mappings_all.sam",
        ),
    output:
        sort=os.path.join(
            config["output_dir"], "{sample}", "mappings_all_sorted_by_id.sam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sort_maps_by_id_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "sort_maps_by_id_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort -n -o {output.sort} {input.concatenate}) &> {log}"


###############################################################################
### Remove inferior mappings (keeping multimappers)
###############################################################################


rule remove_inferiors:
    input:
        sort=os.path.join(
            config["output_dir"], "{sample}", "mappings_all_sorted_by_id.sam"
        ),
        script=os.path.join(
            config["scripts_dir"],
            "sam_remove_duplicates_inferior_alignments_multimappers.pl",
        ),
    output:
        remove_inf=os.path.join(
            config["output_dir"],
            "{sample}",
            "mappings_all_removed_inferiors.sam",
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
    container:
        "docker://perl:5.37.10"
    shell:
        "(perl {input.script} \
        --print-header \
        --keep-mm \
        --in {input.sort} \
        --out {output.remove_inf} \
        ) &> {log}"


###############################################################################
### Filter multimappers by indels
###############################################################################


rule filter_by_indels:
    input:
        sam=os.path.join(
            config["output_dir"],
            "{sample}",
            "mappings_all_removed_inferiors.sam",
        ),
        script=os.path.join(
            config["scripts_dir"],
            "filter_multimappers.py",
        ),
    output:
        sam=os.path.join(config["output_dir"], "{sample}", "alignments_all.sam"),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "remove_multimappers_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "remove_multimappers_{sample}.log"),
    resources:
        mem=15,
        threads=4,
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    shell:
        "(python {input.script} \
        {input.sam} \
        --nh \
        > {output.sam} \
        ) &> {log}"


###############################################################################
### Convert SAM to BAM
###############################################################################


rule convert_all_alns_sam_to_bam:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "alignments_all.sam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"], "{sample}", "alignments_all.bam"
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "convert_all_alns_sam_to_bam_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "convert_all_alns_sam_to_bam_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_all_alns_bam_by_position:
    input:
        maps=os.path.join(
            config["output_dir"], "{sample}", "alignments_all.bam"
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_all_sorted_{sample}.bam",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sort_all_alns_bam_by_position_{sample}.log"
        ),
    log:
        os.path.join(
            config["local_log"], "sort_all_alns_bam_by_position_{sample}.log"
        ),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_all_alns_bam:
    input:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_all_sorted.bam",
        ),
    output:
        maps=os.path.join(
            config["output_dir"],
            "{sample}",
            "alignments_all_sorted.bam.bai",
        ),
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "index_all_alns_bam_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "index_all_alns_bam_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"
