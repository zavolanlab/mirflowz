# _MIRFLOWZ_

_MIRFLOWZ_ is a [Snakemake][snakemake] workflow for mapping and 
quantifying of miRNA and isomiR quantification.




## Installation

The workflow lives inside this repository and will be available for you to run
after following the installation instructions layed out in this section.

### Cloning the repository

Traverse to the desired path on your file system, then clone the repository and
change into it with:

```bash
git clone https://github.com/zavolanlab/mirflowz.git
cd mirflowz
```

### Dependencies

For improved reproducibility and reusability of the workflow, as well as an
easy means to run it on a high performance computing (HPC) cluster managed,
e.g., by [Slurm][slurm], all steps of the workflow run inside their own
containers. As a consequence, running this workflow has very few individual
dependencies. It does, however, require the package manager [Conda][conda] and
the container engine [Singularity][singularity] to be installed before you
proceed.


### Setting up a virtual environment

If you do not already have [Conda][conda] installed globally on your system,
we recommend that you install [Miniconda][miniconda-installation]. For faster creation of
the environment (and Conda environments in general), you can also install
[Mamba][mamba] on top of Conda. In that case, replace `conda` with `mamba` in
the commands below (particularly in `conda env create`).

Create and activate the environment with necessary dependencies with Conda:

```bash
conda env create -f environment.yml
conda activate mirflowz
```

If you have root permissions for your system and you do not already have
`singularity` installed globally on your system, you must update the Conda
environment using the `environment.root.yml` with command below. Mind that the
you must have the environment activated to updated it.

```bash
conda env update -f environment.root.yml
```

### Testing your installation

Several tests are provided to check the integrity of the installation. Follow
the instructions in this section to make sure the workflow is ready to use.

#### Run test workflow on local machine

Execute the following command to run the test workflow on your local machine:

```bash
bash test/test_workflow_local.sh
```

#### Run test workflow via Slurm

Execute the following command to run the test workflow on a Slurm-managed
high-performance computing (HPC) cluster:

```bash
bash test/test_workflow_slurm.sh
```

> **NOTE:** The Slurm tests were configured to run on the developer's cluster.
> Several files may need to be modified if you would like to run tests (and
> the actual workflow) on other systems. These may possibly include the
> following (relative to the repository root directory), but potentially others
> as well:
>  
> * `jobscript.sh`
> * `test/cluster.json`
> * `test/test_workflow_slurm.sh`
>  
> Consult the manual of
> your batch scheduling system, as well as the section of the Snakemake manual
> dealing with [cluster execution].

#### Rule graph

Execute the following command to generate a rule graph image for the workflow.
The output will be found in the `images/` directory in the repository root.

> **NOTE:** It is essential that you run the rule graph test only _after_ 
> running the test workflow. This is because it requires files to be available
> that will only be created when running the workflow.

```bash
bash test/test_rule_graph.sh
```

#### Clean up test results

After successfully running the tests above, you can run the following command
to remove all artifacts generated by the test runs:

```bash
bash test/test_cleanup.sh
```


## Usage

Now that your virtual environment is set up and the workflow is deployed and
tested, you can go ahead and run the workflow on your samples.

### Running the workflow

Assuming that you are currently inside the repository's root directory, 
create a directory from which you will run your workflow and name it whatever 
you want e.g., `MY_ANALYSIS` and head to it.

```bash
mkdir MY_ANALYSIS
cd MY_ANALYSIS
```
Place your library file(s) here. In addition, create a sample table. Fill it 
with the correct entries. The `sample_table.csv` is a tab-separated  file that must have hte following columns:  

- `sample`. This column contains the library name.  
- `sample_file`. In this column, you must provide the path to the library file.
The path must be relative to the worflow's directory.  
- `adapter`.  This field must contain the adapter sequence in capital letters.  
- `format`. In this field you mast state the library format. It can either be 
`fa` if providing a FASTA file or `fastq` if the library is a FASTQ file.  

You can look at the `test/test_files/sample_table.csv` to know what this file 
must look like, or tu use it as a template.

```bash
touch sample_table.csv
```

> **Note:** We expect the genome and gene annotations to be formatted according
> the style used by Ensembl. Other formats are very likely to lead to problems.
> The miRNA annotation file is expected to originate from miRBase, or follow 
> their exact layout.

Copy the `config_template.yaml` to this directory.

```bash
cp ..config/config_template.yaml ./config.yaml
```
Then, using your editor of choice, adjust the parameters of the `config.yaml`.
The file explains what each of the parameters means and how you can meaningfully
fill them in. 

Accordingly to how you want to run the workflow you can either copy the script 
to run it **locally** with

```bash
cp ../test/test_workflow_local.sh ./run_workflow_local.sh
```
or copy the script to run the workflow on a **cluster via Slurm** along with 
the `cluster.json` file with

```bash
cp ../test/test_workflow_slurm.sh ./run_workflow_slurm.sh
cp ../test/cluster.json ./cluster.json
```

In both cases, the final 9 lines must be removed; this can be done with:

```bash
head -n 9 run_workflow_*.sh
```
Finally, you can optionally copy the script to remove all artifacts generated
by the run:

```bash
cp ../test/test_cleanup.sh ./run_cleanup.sh
```

To start workflow execution, run:

```bash
./run_workflow_slurm.sh
```

> **NOTE:** Check back in the installation section to find more information on
> how to run the workflow on your HPC system. Although we do provide a workflow
> runner to execute the workflow locally (`run_workflow_local.sh`) on your 
> laptop or desktop machine, we recommend against that for real-world data, as 
> the resources requirements for running the workflow are very high(can be >50 
Gigs of memory!).

After successful execution of the workflow, results and logs will be found in
`results/` and `logs/` directories, respectively.

## Workflow description

The _MIRFLOWZ_ workflow starts by processing the user-provided "genome 
resources" while preparing indexes and other contigent resources that will be
used in later steps. Afterwards, the user-provided short read smallRNA-seq
library(ies) will be aligned against the references previously generated. For
increased fidelity, it uses two separate aligning tools, [Segemehl][segemehl] 
and our in-house tool [Oligomap][oligomap]. In both cases, reads are aligned
separately to the genome and the transcriptome. Thereafter, alignments are
merged in a way that only the best alignment (or alignments) of each read is 
(are) kept. Finally, the quantification of miRNA expression is done by
intersecting the alignments with the generated annotation files. Intersections 
are computed with [`bedtools`][bedtools] for one or multiple of mature, primary
transcripts and isomiRs. Reads consistent with each miRNA are counted and 
tabulated.


The schema below is a visual representation of the rules building up the 
workflow and how they are related:

> ![rule-graph][rule-graph]

[bedtools]: <https://github.com/arq5x/bedtools2>
[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[cluster execution]: <https://snakemake.readthedocs.io/en/stable/executing/cluster.html>
[mamba]: <https://github.com/mamba-org/mamba>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[oligomap]: <https://bio.tools/oligomap>
[rule-graph]: images/rule_graph.svg
[segemehl]: <https://www.bioinf.uni-leipzig.de/Software/segemehl/>
[singularity]: <https://sylabs.io/singularity/>
[slurm]: <https://slurm.schedmd.com/documentation.html>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>