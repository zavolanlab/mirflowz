#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf .snakemake/
    rm -rf .tmp/
    rm -rf logs/
    rm -rf results/
    rm -rf snakemake_report_*.html
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

# Run tests
snakemake \
    --snakefile="../workflow/prepare_annotation/Snakefile" \
    --configfile="config_prepare_annotation.yaml" \
    --use-singularity \
    --singularity-args "--bind ${PWD}/../" \
    --cores=4 \
    --printshellcmds \
    --rerun-incomplete \
    --verbose

# Snakemake report
snakemake \
    --snakefile="../workflow/prepare_annotation/Snakefile" \
    --configfile="config_prepare_annotation.yaml" \
    --report="snakemake_report_prepare_annotation.html"

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output.md5"

# Checksum file generated with
# find results/ \
#     -type f \
#     > expected_output.files;
# md5sum $(cat expected_output.files) > expected_output.md5
