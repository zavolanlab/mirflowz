#!/bin/bash

# Tear down test environment
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
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir

# Have to match directories indicated in config.yaml files
mkdir -p logs/cluster/{homo_sapiens/chrY,test_lib}
mkdir -p logs/local/{homo_sapiens/chrY,test_lib}
mkdir -p results/{homo_sapiens/chrY,test_lib}

# Run test
snakemake \
    --snakefile="../../workflow/Snakefile" \
    --cores=256 \
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
    --jobscript="../../jobscript.sh" \
    --jobs=20 \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../" \
    --printshellcmds \
    --rerun-incomplete \
    --verbose

# Snakemake report
snakemake \
    --snakefile="../../workflow/Snakefile" \
    --report="snakemake_report.html"
