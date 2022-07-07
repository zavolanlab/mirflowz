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
    --rulegraph \
    --printshellcmds \
    --dryrun \
    --verbose \
    | dot -Tsvg > "../images/rule_graph_prepare.svg"

# Run test: map workflow
snakemake \
    --snakefile="../workflow/map/Snakefile" \
    --configfile="config_map.yaml" \
    --rulegraph \
    --printshellcmds \
    --dryrun \
    --verbose \
    | dot -Tsvg > "../images/rule_graph_map.svg"
