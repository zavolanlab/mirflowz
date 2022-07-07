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

# Run test: prepare workflow
snakemake \
    --snakefile="../workflow/prepare/Snakefile" \
    --configfile="config_prepare.yaml" \
    --use-singularity \
    --singularity-args "--bind ${PWD}/../" \
    --cores=4 \
    --printshellcmds \
    --rerun-incomplete \
    --verbose

# Run test: map workflow
snakemake \
    --snakefile="../workflow/map/Snakefile" \
    --configfile="config_map.yaml" \
    --use-singularity \
    --singularity-args "--bind ${PWD}/../" \
    --cores=4 \
    --printshellcmds \
    --rerun-incomplete \
    --verbose

# Snakemake report: prepare workflow
snakemake \
    --snakefile="../workflow/prepare/Snakefile" \
    --configfile="config_prepare.yaml" \
    --report="snakemake_report_prepare.html"

# Snakemake report: map workflow
snakemake \
    --snakefile="../workflow/map/Snakefile" \
    --configfile="config_map.yaml" \
    --report="snakemake_report_map.html"

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "expected_output.md5"

# Generate checksum files
# (run only when using new test data and after verifying results!)
# find results/ -type f > expected_output.files;
# md5sum $(cat expected_output.files) > expected_output.md5
