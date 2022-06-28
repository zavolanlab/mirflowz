# mir-prepare-annotation

[Snakemake][snakemake] workflow to download and prepare the necessary files for
smallRNA-seq related pipelines [mir-map][mir-map] and [mir-quant][mir-quant].

The scheme below is a visual representation of an example run of the workflow:

> ![rule-graph-prep-anno][rule-graph-prep-anno]

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
package manager. We recommend that you install [Miniconda][miniconda] for your
system.

For improved reproducibility and reusability of the workflow, as well as an
easy means to run it on a high performance computing (HPC) cluster managed,
e.g., by [Slurm][slurm], all steps of the workflow run in their own container.
As a consequence, running this workflow has very few individual dependencies.
It does, however, require that the container engine [Singularity][singularity]
is installed.

Create and activate the environment with necessary dependencies with Conda:

```bash
conda env create -f environment.yml
conda activate mir-pipelines
```

> **NOTE:** For faster creation of the environment (and Conda environments in
> general), you can also install [Mamba][mamba] on top of Conda.
>  
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

Execute the following commands to generate DAG and rule graph images. Outputs
will be found in the `images/` directory.

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

> **NOTE:** This was set up to run on the developer's Slurm cluster. Several
> files may need to be modified on other systems, including `jobscript.sh`,
> `workflow/prepare_annotation/cluster.json` and `test/test_workflow_slurm.sh`
> itself (all relative to the repository's root directory). Consult the manual
> of your batch scheduling system, as well as the section of the Snakemake
> manual dealing with [cluster execution].

## Usage

Assuming that you are currently inside the repository's root directory, change
to the run root directory:

```bash
cd RUNS
```

Now make a clean copy of the `JOB` directory and name it what you want, e.g.,
`MY_ANALYSIS`:

```bash
cp -r JOB MY_ANALYSIS
```

Now traverse to the directory from where you will actually execute the pipeline
with:

```bash
cd MY_ANALYSIS/prepare_annotation
```

Before running the pipeline adjust the parameters in file
`config_prepare_annotation.yaml`:

```yaml
```

> **Note:** We expect the genome and gene annotations to be formatted according
> the style used by Ensembl. Other formats are very likely to lead to problems,
> if not in this pipeline, then further down the road in the mapping or
> annotation pipelines. The miRNA annotation file is expected to originate from
> miRBase, or follow their exact layout.

To start pipeline execution locally:

```bash
./run_workflow_local.sh
```

To start pipeline execution via Slurm:

```bash
./run_workflow_slurm.sh
```

> *This is strongly recommended due to excessive resource needs of some tools!*

After successful execution of the workflow, results and logs will be found in
`results/` and `logs/` directories, respectively.

> **Note:** See the note in the installation section for configuring workflow
> runs on your individual batch scheduling system, Slurm or otherwise.

[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[cluster execution]: <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution>
[mamba]: <https://github.com/mamba-org/mamba>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[mir-map]: <https://git.scicore.unibas.ch/zavolan_group/pipelines/mir-map>
[mir-quant]: <https://git.scicore.unibas.ch/zavolan_group/pipelines/mir-quant>
[rule-graph-prep-anno]: images/rule_graph_prepare_annotation.svg
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[singularity]: <https://sylabs.io/singularity/>
[slurm]: <https://slurm.schedmd.com/documentation.html>
