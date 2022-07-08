#!/bin/bash

# Tear down environment
cleanup () {
    rc=$?
    rm $(cat intermediate_files.txt)
    cd $user_dir
    echo "Exit status: $rc"
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir

# Have to match directories indicated in config.yaml
mkdir -p logs/cluster
mkdir -p logs/local
mkdir -p results


# Run workflow
snakemake \
    --snakefile="../../../workflow/map/Snakefile" \
    --configfile="config.yaml" \
    --cluster-config="cluster.json" \
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
    --snakefile="../../../workflow/map/Snakefile" \
    --configfile="config.yaml" \
    --report="snakemake_report.html"
