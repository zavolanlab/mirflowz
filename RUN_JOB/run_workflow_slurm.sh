#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    #rm -rf .snakemake/
    #rm -rf logs/
    #rm -rf results/
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

# Run tests
snakemake \
    --snakefile="../snakemake/Snakefile" \
    --configfile="config.yaml" \
    --cluster-config="../cluster.json" \
    --cores=256 \
    --jobscript="../jobscript.sh" \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args="--no-home --bind ${PWD},/scicore/home/zavolan/devagy74/projects" \
    --cluster "sbatch \
        --cpus-per-task={cluster.threads} \
        --mem={cluster.mem} \
        --qos={cluster.queue} \
        --time={cluster.time} \
        --export=JOB_NAME={rule} \
        -o {params.cluster_log} \
        -p scicore \
        --open-mode=append" \
    --verbose

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
md5sum --check "expected_output.md5"

# Checksum file generated with
# find results/ \
#     -type f \
#     -name \*\.gz \
#     -exec gunzip '{}' \;
#     > expected_output.files
# md5sum $(cat expected_output.files) > expected_output.md5