#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    cd $PWD
    echo "Exit status: $rc"
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands

# Store root and test directories
ROOT="$(git rev-parse --show-toplevel)"

cd $ROOT

# Run tests
snakefmt  --check "$ROOT/workflow"
