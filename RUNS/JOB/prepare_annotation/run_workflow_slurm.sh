#!/bin/bash

# Tear down environment
cleanup () {
    rc=$?
    cd $user_dir
    echo "Exit status: $rc"
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands

#### CHANGE PATHS WITH YOUR ORGANISM AND PREFIX_NAME ####
mkdir -p logs/cluster/ORGANISM/PREFIX_NAME
mkdir -p logs/local/ORGANISM/PREFIX_NAME
mkdir -p results/ORGANISM/PREFIX_NAME
#########################################################

user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir

# Run workflow
snakemake \
    --snakefile="../../../workflow/prepare_annotation/Snakefile" \
    --configfile="config.yaml" \
    --cluster-config="../../../workflow/prepare_annotation/cluster.json" \
    --cluster "sbatch \
        --cpus-per-task={cluster.threads} \
        --mem={cluster.mem} \
        --qos={cluster.queue} \
        --time={cluster.time} \
        --export=JOB_NAME={rule} \
        -o {params.cluster_log} \
        -p scicore \
        --open-mode=append" \
    --use-singularity \
    --singularity-args="--no-home --bind ${PWD}/../../../" \
    --jobscript="../../../jobscript.sh" \
    --cores=256 \
    --printshellcmds \
    --rerun-incomplete \
    --verbose

# Snakemake report
snakemake \
    --snakefile="../../../workflow/prepare_annotation/Snakefile" \
    --configfile="config.yaml" \
    --report="snakemake_report.html"
