#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    cd $PWD
    echo "Exit status: $rc"
    rm -rf .snakemake
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands

# Store root and test directories
ROOT="$(git rev-parse --show-toplevel)"
TEST="${ROOT}/test"

cd $ROOT

# Run test
snakemake \
    --snakefile="$ROOT/workflow/Snakefile" \
    --cores 4  \
    --configfile="$TEST/test_files/config.yaml" \
    --software-deployment-method conda \
    --printshellcmds \
    --show-failed-logs \
    --rerun-incomplete \
    --no-hooks \
    --verbose

# Snakemake report
snakemake \
    --snakefile="$ROOT/workflow/Snakefile" \
    --configfile="$TEST/test_files/config.yaml" \
    --report="snakemake_report.html"

# Check md5 sum of some output files
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
find results/ -type f -name \*\.zip -exec sh -c 'unzip -o {} -d $(dirname {})' \;
md5sum --check "$TEST/test_integration_workflow/expected_output.md5"

# Generate checksum files
# (run only when using new test data and after verifying results!)
# md5sum $(find results/ -type f) > "$TEST/test_integration_workflow/expected_output.md5"
