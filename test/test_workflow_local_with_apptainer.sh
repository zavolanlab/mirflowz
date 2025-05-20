#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    cd $user_dir
    echo "Exit status: $rc"
    rm -rf .snakemake
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir

# Run test
snakemake \
    --snakefile="../workflow/Snakefile" \
    --cores 4  \
    --configfile="config.yaml" \
    --software-deployment-method apptainer \
    --apptainer-args "--bind ${PWD}/../" \
    --printshellcmds \
    --rerun-incomplete \
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
