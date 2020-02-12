#!/bin/bash

# Tear down test environment
trap 'rm -rf .snakemake/ && cd $user_dir' EXIT  # quotes command is exected after script exits, regardless of exit status

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir
mkdir -p logs/cluster_log
mkdir -p logs/local_log

# Run tests
snakemake \
-p \
--snakefile="../snakemake/Snakefile" \
--use-singularity \
--singularity-args "--no-home --bind ${PWD},/scicore/home/zavolan/devagy74/projects" \
--cores 256 \
--local-cores 10 \
--rerun-incomplete \
--configfile config.yaml \
--jobscript ../jobscript.sh \
--cluster-config ../cluster.json \
--cluster "sbatch \
	--cpus-per-task={cluster.threads} \
	--mem={cluster.mem} \
	--qos={cluster.queue} \
	--time={cluster.time} \
	-o {params.cluster_log} \
	-p scicore \
	--export=JOB_NAME={rule} \
	--open-mode=append"
find results/ -type f -name \*\.gz -exec gunzip '{}' \;
md5sum --check "test.md5"

# Checksum file generated with
# find results/ \
#     -type f \
#     -exec md5sum '{}' \; \
#     > expected_output.md5
