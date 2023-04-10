#!/bin/bash

# Setting to strict mode
set -e
set -u
set -o pipefail

# Trimming fasta ID
awk -F" " "/^>/ {{print \$1; next}} 1"