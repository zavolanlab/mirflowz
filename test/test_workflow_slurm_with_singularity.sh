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

# Create directories defined as in the configuration files
mkdir -p logs/cluster
mkdir -p logs/local
mkdir -p results

# Run test
snakemake \
    --snakefile="../workflow/Snakefile" \
    --cores=256 \
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
    --jobscript="../jobscript.sh" \
    --jobs=20 \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../" \
    --printshellcmds \
    --rerun-incomplete \
    --latency-wait 60 \
    --no-hooks \
    --verbose

# Snakemake report
snakemake \
    --snakefile="../workflow/Snakefile" \
    --configfile="config.yaml" \
    --report="snakemake_report.html"

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output.md5"

# Generate checksum files
# (run only when using new test data and after verifying results!)
# md5sum $(find results/ -type f) > expected_output.md5
