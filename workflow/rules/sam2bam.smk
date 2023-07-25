###############################################################################
# SAM to BAM process
#
# The following set of rules takes a SAM file, converts it into BAM file,
# sorts it by starting position and generates the corresponding BAI file.
###############################################################################

import os

###############################################################################
### Convert SAM to BAM
###############################################################################


rule convert_to_bam:
    input:
        sam="input.sam",
    output:
        bam="output.bam",
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "convert_to_bam_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "convert_to_bam_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools view -b {input.sam} > {output.bam}) &> {log}"


###############################################################################
### Sort by coordinate position
###############################################################################


rule sort_by_position:
    input:
        bam="output.bam",
    output:
        bamS="sorted.bam",
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "sort_by_position_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "sort_by_position_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools sort {input.bam} > {output.bamS}) &> {log}"


###############################################################################
### Create bam index
###############################################################################


rule index_bam:
    input:
        bamS="sorted.bam",
    output:
        bai="sorted.bam.bai",
    params:
        cluster_log=os.path.join(
            config["cluster_log"], "index_bam_{sample}.log"
        ),
    log:
        os.path.join(config["local_log"], "index_bam_{sample}.log"),
    container:
        "docker://quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"
    shell:
        "(samtools index -b {input.bamS} > {output.bai}) &> {log}"

