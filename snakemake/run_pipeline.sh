# set -e

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
-p \
-s Snakefile \
--use-singularity \
--singularity-args "--no-home --bind ${PWD},/scicore/home/zavolan/GROUP/test_datasets/annotation/" \
--cores 256 \
--local-cores 10 \
--rerun-incomplete \
--configfile config.yaml \
--jobscript jobscript.sh \
--cluster-config cluster.json \
--cluster "sbatch \
	--cpus-per-task={cluster.threads} \
	--mem={cluster.mem} \
	--qos={cluster.queue} \
	--time={cluster.time} \
	-o {params.cluster_log} \
	-p scicore \
	--export=JOB_NAME={rule} \
	--open-mode=append"