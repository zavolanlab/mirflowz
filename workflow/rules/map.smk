###############################################################################
# (c) 2020 Paula Iborra, Zavolan Lab, Biozentrum, University of Basel
# (@) paula.iborradetoledo@unibas.ch / paula.iborra@alumni.esci.upf.edu
#
# Workflow to map small RNA-seq reads (e.g. from miRNA sequencing libraries).
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
OUT_DIR = Path(config["output_dir"])
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
    start,
    finish_map,


###############################################################################
### Finish rule
###############################################################################


rule finish_map:
    input:
        maps=expand(
            OUT_DIR / "{sample}" / "alignments_all_sorted_{sample}.bam.bai",
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
            format=convert_lib_format(get_sample("format")),
        ),
    output:
        reads=OUT_DIR / "{sample}" / "{format}" / "reads.{format}",
    params:
        cluster_log=CLUSTER_LOG
        / "uncompress_zipped_files_{sample}_{format}.log",
    log:
        LOCAL_LOG / "uncompress_zipped_files_{sample}_{format}.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(zcat {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Quality filter
###############################################################################


rule fastq_quality_filter:
    input:
        reads=OUT_DIR / "{sample}" / "fastq" / "reads.fastq",
    output:
        reads=OUT_DIR / "{sample}" / "fastq" / "filtered_reads.fastq",
    params:
        cluster_log=CLUSTER_LOG / "fastq_quality_filter_{sample}.log",
        p=config["p_value"],
        q=config["q_value"],
    log:
        LOCAL_LOG / "fastq_quality_filter_{sample}.log",
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--h87f3376_10"
    conda:
        ENV_DIR / "fastx_toolkit.yaml"
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
        reads=OUT_DIR / "{sample}" / "fastq" / "filtered_reads.fastq",
    output:
        reads=OUT_DIR / "{sample}" / "fastq" / "reads.fa",
    params:
        cluster_log=CLUSTER_LOG / "fastq_to_fasta_{sample}.log",
    log:
        LOCAL_LOG / "fastq_to_fasta_{sample}.log",
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--h87f3376_10"
    conda:
        ENV_DIR / "fastx_toolkit.yaml"
    shell:
        "(fastq_to_fasta -r -n -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Format fasta file
###############################################################################


rule format_fasta:
    input:
        reads=lambda wildcards: OUT_DIR
        / wildcards.sample
        / convert_lib_format(get_sample("format", wildcards.sample))
        / "reads.fa",
    output:
        reads=OUT_DIR / "{sample}" / "reads_formatted.fasta",
    params:
        cluster_log=CLUSTER_LOG / "format_fasta_{sample}.log",
    log:
        LOCAL_LOG / "format_fasta_{sample}.log",
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--h87f3376_10"
    conda:
        ENV_DIR / "fastx_toolkit.yaml"
    shell:
        "(fasta_formatter -w 0 -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Remove adapters
###############################################################################


rule remove_adapters:
    input:
        reads=OUT_DIR / "{sample}" / "reads_formatted.fasta",
    output:
        reads=OUT_DIR / "{sample}" / "reads_trimmed_adapters.fasta",
    params:
        adapter=lambda wildcards: get_sample(
            "adapter", wildcards.sample
        ).upper(),
        error_rate=config["error_rate"],
        minimum_length=config["minimum_length"],
        overlap=config["overlap"],
        max_n=config["max_n"],
        cluster_log=CLUSTER_LOG / "remove_adapters_{sample}.log",
    log:
        LOCAL_LOG / "remove_adapters_{sample}.log",
    resources:
        threads=8,
    container:
        "docker://quay.io/biocontainers/cutadapt:4.3--py310h1425a21_0"
    conda:
        ENV_DIR / "cutadapt.yaml"
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
        reads=OUT_DIR / "{sample}" / "reads_trimmed_adapters.fasta",
    output:
        reads=OUT_DIR / "{sample}" / "reads_collapsed.fasta",
    params:
        cluster_log=CLUSTER_LOG / "collapse_identical_reads_{sample}.log",
    log:
        LOCAL_LOG / "collapse_identical_reads_{sample}.log",
    container:
        "docker://quay.io/biocontainers/fastx_toolkit:0.0.14--h87f3376_10"
    conda:
        ENV_DIR / "fastx_toolkit.yaml"
    shell:
        "(fastx_collapser -i {input.reads} > {output.reads}) &> {log}"


###############################################################################
### Segemehl genome mapping
###############################################################################


rule map_genome_segemehl:
    input:
        reads=OUT_DIR / "{sample}" / "reads_collapsed.fasta",
        genome=OUT_DIR / "genome_processed.fa",
        genome_index_segemehl=OUT_DIR / "segemehl_genome_index.idx",
    output:
        gmap=OUT_DIR / "{sample}" / "segemehl_genome_mappings.sam",
    params:
        cluster_log=CLUSTER_LOG / "map_genome_segemehl_{sample}.log",
    log:
        LOCAL_LOG / "map_genome_segemehl_{sample}.log",
    resources:
        mem=50,
        time=12,
        threads=8,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    conda:
        ENV_DIR / "segemehl.yaml"
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
        reads=OUT_DIR / "{sample}" / "reads_collapsed.fasta",
        transcriptome=OUT_DIR / "transcriptome_trimmed_id.fa",
        transcriptome_index_segemehl=OUT_DIR
        / "segemehl_transcriptome_index.idx",
    output:
        tmap=OUT_DIR / "{sample}" / "segemehl_transcriptome_mappings.sam",
    params:
        cluster_log=CLUSTER_LOG / "map_transcriptome_segemehl_{sample}.log",
    log:
        LOCAL_LOG / "map_transcriptome_segemehl_{sample}.log",
    resources:
        mem=10,
        time=12,
        threads=8,
    container:
        "docker://quay.io/biocontainers/segemehl:0.3.4--hf7d323f_8"
    conda:
        ENV_DIR / "segemehl.yaml"
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
        reads=OUT_DIR / "{sample}" / "reads_collapsed.fasta",
        script=SCRIPTS_DIR / "validation_fasta.py",
    output:
        reads=OUT_DIR / "{sample}" / "reads_filtered_for_oligomap.fasta",
    params:
        cluster_log=CLUSTER_LOG / "filter_fasta_for_oligomap_{sample}.log",
        max_length_reads=config["max_length_reads"],
    log:
        LOCAL_LOG / "filter_fasta_for_oligomap_{sample}.log",
    container:
        "docker://python:3.9.16"
    conda:
        ENV_DIR / "python.yaml"
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
        reads=OUT_DIR / "{sample}" / "reads_filtered_for_oligomap.fasta",
        target=OUT_DIR / "genome_processed.fa",
    output:
        gmap=OUT_DIR / "{sample}" / "oligomap_genome_mappings.fasta",
        report=OUT_DIR / "{sample}" / "oligomap_genome_report.txt",
    params:
        cluster_log=CLUSTER_LOG / "map_genome_oligomap_{sample}.log",
    log:
        LOCAL_LOG / "map_genome_oligomap_{sample}.log",
    resources:
        mem=50,
        time=6,
        threads=8,
    container:
        "docker://quay.io/biocontainers/oligomap:1.0.1--hdcf5f25_0"
    conda:
        ENV_DIR / "oligomap.yaml"
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
        tmap=OUT_DIR / "{sample}" / "oligomap_genome_mappings.fasta",
        report=OUT_DIR / "{sample}" / "oligomap_genome_report.txt",
        script=SCRIPTS_DIR / "blocksort.sh",
    output:
        sort=OUT_DIR / "{sample}" / "oligomap_genome_sorted.fasta",
    params:
        cluster_log=CLUSTER_LOG / "sort_genome_oligomap_{sample}.log",
    log:
        LOCAL_LOG / "sort_genome_oligomap_{sample}.log",
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


rule convert_genome_to_sam_oligomap:
    input:
        sort=OUT_DIR / "{sample}" / "oligomap_genome_sorted.fasta",
        script=SCRIPTS_DIR / "oligomap_output_to_sam_nh_filtered.py",
    output:
        gmap=OUT_DIR / "{sample}" / "oligomap_genome_mappings.sam",
    params:
        cluster_log=CLUSTER_LOG / "oligomap_genome_to_sam_{sample}.log",
        nh=config["nh"],
    log:
        LOCAL_LOG / "oligomap_genome_to_sam_{sample}.log",
    resources:
        time=1,
        queue=1,
    container:
        "docker://python:3.9.16"
    conda:
        ENV_DIR / "python.yaml"
    shell:
        "(python {input.script} \
        {input.sort} \
        -n {params.nh} \
        > {output.gmap}) &> {log}"


###############################################################################
### Oligomap transcriptome mapping
###############################################################################


rule map_transcriptome_oligomap:
    input:
        reads=OUT_DIR / "{sample}" / "reads_filtered_for_oligomap.fasta",
        target=OUT_DIR / "transcriptome_trimmed_id.fa",
    output:
        tmap=OUT_DIR / "{sample}" / "oligomap_transcriptome_mappings.fasta",
        report=OUT_DIR / "{sample}" / "oligomap_transcriptome_report.txt",
    params:
        cluster_log=CLUSTER_LOG / "map_transcriptome_oligomap_{sample}.log",
    log:
        LOCAL_LOG / "map_transcriptome_oligomap_{sample}.log",
    resources:
        mem=10,
        time=6,
        threads=8,
    container:
        "docker://quay.io/biocontainers/oligomap:1.0.1--hdcf5f25_0"
    conda:
        ENV_DIR / "oligomap.yaml"
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
        tmap=OUT_DIR / "{sample}" / "oligomap_transcriptome_mappings.fasta",
        report=OUT_DIR / "{sample}" / "oligomap_transcriptome_report.txt",
        script=SCRIPTS_DIR / "blocksort.sh",
    output:
        sort=OUT_DIR / "{sample}" / "oligomap_transcriptome_sorted.fasta",
    params:
        cluster_log=CLUSTER_LOG / "sort_transcriptome_oligomap_{sample}.log",
    log:
        LOCAL_LOG / "sort_transcriptome_oligomap_{sample}.log",
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


rule convert_transcriptome_to_sam_oligomap:
    input:
        sort=OUT_DIR / "{sample}" / "oligomap_transcriptome_sorted.fasta",
        script=SCRIPTS_DIR / "oligomap_output_to_sam_nh_filtered.py",
    output:
        tmap=OUT_DIR / "{sample}" / "oligomap_transcriptome_mappings.sam",
    params:
        cluster_log=CLUSTER_LOG / "oligomap_transcriptome_to_sam_{sample}.log",
        nh=config["nh"],
    log:
        LOCAL_LOG / "oligomap_transcriptome_to_sam_{sample}.log",
    container:
        "docker://python:3.9.16"
    conda:
        ENV_DIR / "python.yaml"
    shell:
        "(python {input.script} \
        {input.sort} \
        -n {params.nh} \
        > {output.tmap}) &> {log}"



###############################################################################
### Merge genome mappings
###############################################################################


rule merge_genome_maps:
    input:
        gmap1=OUT_DIR / "{sample}" / "segemehl_genome_mappings.sam",
        gmap2=OUT_DIR / "{sample}" / "oligomap_genome_mappings.sam",
    output:
        gmaps=OUT_DIR / "{sample}" / "genome_mappings.sam",
    params:
        cluster_log=CLUSTER_LOG / "merge_genome_maps_{sample}.log",
    log:
        LOCAL_LOG / "merge_genome_maps_{sample}.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.gmap1} {input.gmap2} > {output.gmaps}) &> {log}"


###############################################################################
### Merge transcriptome mappings
###############################################################################


rule merge_transcriptome_maps:
    input:
        tmap1=OUT_DIR / "{sample}" / "segemehl_transcriptome_mappings.sam",
        tmap2=OUT_DIR / "{sample}" / "oligomap_transcriptome_mappings.sam",
    output:
        tmaps=OUT_DIR / "{sample}" / "transcriptome_mappings.sam",
    params:
        cluster_log=CLUSTER_LOG / "merge_transcriptome_maps_{sample}.log",
    log:
        LOCAL_LOG / "merge_transcriptome_maps_{sample}.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.tmap1} {input.tmap2} > {output.tmaps}) &> {log}"


###############################################################################
### Filter NH genome
###############################################################################


rule filter_genome_by_nh:
    input:
        gmaps=OUT_DIR / "{sample}" / "genome_mappings.sam",
        script=SCRIPTS_DIR / "nh_filter.py",
    output:
        gmaps=OUT_DIR / "{sample}" / "genome_mappings_filtered_nh.sam",
    params:
        cluster_log=CLUSTER_LOG / "filter_genome_by_nh_{sample}.log",
        nh=config["nh"],
    log:
        LOCAL_LOG / "filter_genome_by_nh_{sample}.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    conda:
        ENV_DIR / "pysam.yaml"
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
        tmaps=OUT_DIR / "{sample}" / "transcriptome_mappings.sam",
        script=SCRIPTS_DIR / "nh_filter.py",
    output:
        tmaps=OUT_DIR / "{sample}" / "transcriptome_mappings_filtered_nh.sam",
    params:
        cluster_log=CLUSTER_LOG / "filter_transcriptome_by_nh_{sample}.log",
        nh=config["nh"],
    log:
        LOCAL_LOG / "filter_transcriptome_by_nh_{sample}.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    conda:
        ENV_DIR / "pysam.yaml"
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
        gmap=OUT_DIR / "{sample}" / "genome_mappings_filtered_nh.sam",
    output:
        gmap=OUT_DIR / "{sample}" / "genome_mappings_no_header.sam",
    params:
        cluster_log=CLUSTER_LOG / "remove_header_genome_mappings_{sample}.log",
    log:
        LOCAL_LOG / "remove_header_genome_mappings_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "samtools view {input.gmap} > {output.gmap}"


###############################################################################
### Remove header transcriptome mappings
###############################################################################


rule remove_header_transcriptome_mappings:
    input:
        tmap=OUT_DIR / "{sample}" / "transcriptome_mappings_filtered_nh.sam",
    output:
        tmap=OUT_DIR / "{sample}" / "transcriptome_mappings_no_header.sam",
    params:
        cluster_log=CLUSTER_LOG
        / "remove_header_transcriptome_mappings_{sample}.log",
    log:
        LOCAL_LOG / "remove_header_transcriptome_mappings_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "samtools view {input.tmap} > {output.tmap}"


###############################################################################
### Transcriptome to genome coordinates
###############################################################################


rule transcriptome_to_genome_maps:
    input:
        tmap=OUT_DIR / "{sample}" / "transcriptome_mappings_no_header.sam",
        script=SCRIPTS_DIR / "sam_trx_to_sam_gen.pl",
        exons=OUT_DIR / "exons.bed",
    output:
        genout=OUT_DIR / "{sample}" / "transcriptome_mappings_to_genome.sam",
    params:
        cluster_log=CLUSTER_LOG / "transcriptome_to_genome_maps_{sample}.log",
    log:
        LOCAL_LOG / "transcriptome_to_genome_maps_{sample}.log",
    container:
        "docker://perl:5.37.10"
    conda:
        ENV_DIR / "perl.yaml"
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
        gmap1=OUT_DIR / "{sample}" / "transcriptome_mappings_to_genome.sam",
        gmap2=OUT_DIR / "{sample}" / "genome_mappings_no_header.sam",
    output:
        catmaps=OUT_DIR / "{sample}" / "mappings_all_no_header.sam",
    params:
        cluster_log=CLUSTER_LOG / "merge_all_mappings_{sample}.log",
    log:
        LOCAL_LOG / "merge_all_mappings_{sample}.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.gmap1} {input.gmap2} > {output.catmaps}) &> {log}"


###############################################################################
### Add header
###############################################################################


rule add_header_all_maps:
    input:
        header=OUT_DIR / "genome_header.sam",
        catmaps=OUT_DIR / "{sample}" / "mappings_all_no_header.sam",
    output:
        concatenate=OUT_DIR / "{sample}" / "mappings_all.sam",
    params:
        cluster_log=CLUSTER_LOG / "add_header_{sample}.log",
    log:
        LOCAL_LOG / "add_header_{sample}.log",
    container:
        "docker://ubuntu:lunar-20221207"
    shell:
        "(cat {input.header} {input.catmaps} > {output.concatenate}) &> {log}"


###############################################################################
### Sort mapped file by IDs
###############################################################################


rule sort_maps_by_id:
    input:
        concatenate=OUT_DIR / "{sample}" / "mappings_all.sam",
    output:
        sort=OUT_DIR / "{sample}" / "mappings_all_sorted_by_id.sam",
    params:
        cluster_log=CLUSTER_LOG / "sort_maps_by_id_{sample}.log",
    log:
        LOCAL_LOG / "sort_maps_by_id_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools sort -n -o {output.sort} {input.concatenate}) &> {log}"


###############################################################################
### Remove inferior mappings (keeping multimappers)
###############################################################################


rule remove_inferiors:
    input:
        sort=OUT_DIR / "{sample}" / "mappings_all_sorted_by_id.sam",
        script=SCRIPTS_DIR
        / "sam_remove_duplicates_inferior_alignments_multimappers.pl",
    output:
        remove_inf=OUT_DIR / "{sample}" / "mappings_all_removed_inferiors.sam",
    params:
        cluster_log=CLUSTER_LOG / "remove_inferiors_{sample}.log",
    log:
        LOCAL_LOG / "remove_inferiors_{sample}.log",
    resources:
        mem=15,
        threads=4,
    container:
        "docker://perl:5.37.10"
    conda:
        ENV_DIR / "perl.yaml"
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
        sam=OUT_DIR / "{sample}" / "mappings_all_removed_inferiors.sam",
        script=SCRIPTS_DIR / "filter_multimappers.py",
    output:
        sam=OUT_DIR / "{sample}" / "alignments_all.sam",
    params:
        cluster_log=CLUSTER_LOG / "remove_multimappers_{sample}.log",
    log:
        LOCAL_LOG / "remove_multimappers_{sample}.log",
    resources:
        mem=15,
        threads=4,
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    conda:
        ENV_DIR / "pysam.yaml"
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
        maps=OUT_DIR / "{sample}" / "alignments_all.sam",
    output:
        maps=OUT_DIR / "{sample}" / "alignments_all.bam",
    params:
        cluster_log=CLUSTER_LOG / "convert_all_alns_sam_to_bam_{sample}.log",
    log:
        LOCAL_LOG / "convert_all_alns_sam_to_bam_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools view -b {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_all_alns_bam_by_position:
    input:
        maps=OUT_DIR / "{sample}" / "alignments_all.bam",
    output:
        maps=OUT_DIR / "{sample}" / "alignments_all_sorted_{sample}.bam",
    params:
        cluster_log=CLUSTER_LOG / "sort_all_alns_bam_by_position_{sample}.log",
    log:
        LOCAL_LOG / "sort_all_alns_bam_by_position_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools sort {input.maps} > {output.maps}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_all_alns_bam:
    input:
        maps=OUT_DIR / "{sample}" / "alignments_all_sorted.bam",
    output:
        maps=OUT_DIR / "{sample}" / "alignments_all_sorted.bam.bai",
    params:
        cluster_log=CLUSTER_LOG / "index_all_alns_bam_{sample}.log",
    log:
        LOCAL_LOG / "index_all_alns_bam_{sample}.log",
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    conda:
        ENV_DIR / "samtools.yaml"
    shell:
        "(samtools index -b {input.maps} > {output.maps}) &> {log}"
