# mir-prepare-annotation

[Snakemake] workflow to download and prepare the necessary files for
smallRNA-seq related pipelines
[mir-map](https://git.scicore.unibas.ch/zavolan_group/pipelines/mir-map)
and
[mir-quant](https://git.scicore.unibas.ch/zavolan_group/pipelines/mir-quant).

The scheme below is a visual representation of an example run of the
workflow:

> ![workflow_dag](images/rule_graph.svg)

## Installation

### Cloning the repository

Traverse to the desired path on your file system, then clone the repository and
move into it with:

```bash
git clone ssh://git@git.scicore.unibas.ch:2222/zavolan_group/pipelines/mir-prepare-annotation.git
cd mir-prepare-annotation
```

### Setting up a virtual environment

Workflow dependencies can be conveniently installed with the [Conda][conda]
package manager. We recommend that you install
[Miniconda][miniconda-installation] for your system.

For improved reproducibility and reusability of the workflow, as well as an
easy means to run it on a high performance computing (HPC) cluster managed,
e.g., by [Slurm][slurm], all steps of the workflow run in their own container.
As a consequence, running this workflow has very few individual dependencies. It
does, however, require that the container engine [Singularity][singularity] is
installed.

Create and activate the environment with necessary dependencies with conda:

```bash
conda env create -f environment.yml
conda activate mir-pipelines
```

> **NOTE:** If you have root permissions for your system and you do not have
> `singularity` installed globally on your system, you can use Conda to install
> it. In that case, replace `environment.yml` with `environment.root.yml` in
> the first command above.

### Testing your installation

Several tests are prepared to check the integrity of the workflow.

Change into the test directory:

```bash
cd test/
```

#### DAG and rule graph

Execute the following commands to generete DAG and rule graph images. Outputs
will be found in `images/` folder.

```bash
./test_dag.sh
./test_rule_graph.sh
```

#### Run workflow on local machine

Execute the following command to run the test workflow on your local machine:

```bash
./test_workflow_local.sh
```

#### Run workflow via Slurm

Execute the following command to run the test workflow on a Slurm-managed
high-performance computing (HPC) cluster:

```bash
./test_workflow_slurm.sh
```

## Usage

Assuming that you are currently inside the repository's root directory, change
to the execution directory:

```bash
cd RUN_JOB
```

Before running the pipeline adjust the parameters in file `config.yaml`. We
recommend that you create a copy of it for each run.

To start pipeline execution locally:

```bash
./run_workflow_local.sh
```

To start pipeline execution via Slurm:

```bash
./run_workflow_slurm.sh
```

> *This is strongly recommended due to excessive resource needs of some tools!*

After succesful execution of the workflow, results and logs will be found in
`results/` and `logs/` directories, respectively.

> **NOTE:** Depending on the configuration of your Slurm installation or if
> using a different batch scheduling system, you may need to adapt file
> `cluster.json` (located in root directory) and the arguments to options
> `--config` and `--cores` in file `run_workflow_slurm.sh`, located in the
> `RUN_JOB` directory. Consult the manual of your batch scheduling system, as
> well as the section of the Snakemake manual dealing with [cluster execution].

[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[cluster execution]: <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[rule-graph]: images/rule_graph.svg
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[singularity]: <https://sylabs.io/singularity/>
[slurm]: <https://slurm.schedmd.com/documentation.html>

