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

# Run workflow
snakemake \
    --printshellcmds \
    --snakefile="../snakemake/Snakefile" \
    --use-singularity \
    --singularity-args "--bind ${PWD}/../" \
    --cores=4 \
    --rerun-incomplete \
    --configfile="config.yaml" \
    --verbose

# Snakemake report
snakemake \
    --snakefile="../snakemake/Snakefile" \
    --configfile="config.yaml" \
    --report="snakemake_report.html"
