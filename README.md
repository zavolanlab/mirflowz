# _MIRFLOWZ_

_MIRFLOWZ_ is a [Snakemake][snakemake] workflow for mapping miRNAs and isomiRs.

## Table of Contents

1. [Installation](#installation)
    - [Cloning the repository](#cloning-the-repository)
    - [Dependencies](#dependencies)
    - [Setting up the virtual environment](#setting-up-the-virtual-environment)
    - [Testing your installation](#testing-your-installation)
2. [Usage](#usage)
    - [Preparing inputs](#preparing-inputs)
    - [Running the workflow](#running-the-workflow)
    - [Creating a Snakemake report](#creating-a-snakemake-report)
3. [Workflow description](#workflow-description)
4. [Contributing](#contributing)
5. [License](#license)
6. [Contact](#contact)

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
e.g., by [Slurm][slurm], all steps of the workflow run inside isolated
environments ([Singularity][singularity] containers or [Conda][conda]
environments). As a consequence, running this workflow has only a few individual
dependencies. These are managed by the package manager Conda, which 
needs to be installed on your system before proceeding.

If you do not already have Conda installed globally on your system,
we recommend that you install [Miniconda][miniconda-installation]. For faster
creation of the environment (and Conda environments in general), you can also
install [Mamba][mamba] on top of Conda. In that case, replace `conda` with
`mamba` in the commands below (particularly in `conda env create`).

### Setting up the virtual environment

Create and activate the environment with necessary dependencies with Conda:

```bash
conda env create -f environment.yml
conda activate mirflowz
```

If you plan to run _MIRFLOWZ_ via Conda, we recommend to use the following
command for a faster environment creation specially if it you will run it on a
HPC cluster.

```bash
conda config --set channel_priority strict
```

If you plan to run _MIRFLOWZ_ via Singularity and do not already
have it installed globally on your system, you must further update the Conda
environment using the `environment.root.yml` with the command below.
Mind that you must have the environment activated to update it.

```bash
conda env update -f environment.root.yml
```

> Note that you will need to have root permissions on your system to be able
> to install Singularity. If you want to run _MIRFLOWZ_ on an HPC cluster
> (recommended in almost all cases), ask your systems administrator about
> Singularity.

If you would like to contribute to _MIRFLOWZ_ development, you may find it 
useful to further update your environment with the development dependencies:

```bash
conda env update -f environment.dev.yml
```

### Testing your installation

Several tests are provided to check the integrity of the installation. Follow
the instructions in this section to make sure the workflow is ready to use.

#### Run test workflow on local machine

Execute one of the following commands to run the test workflow on your local
machine:

- Test workflow on local machine with **Singularity**:

```bash
bash test/test_workflow_local_with_singularity.sh
```

- Test workflow on local machine with **Conda**:

```bash
bash test/test_workflow_local_with_conda.sh
```

#### Run test workflow on a cluster via SLURM

Execute one of the following commands to run the test workflow on a
slurm-managed high-performance computing (HPC) cluster:

- Test workflow with **Singularity**:

```bash
bash test/test_workflow_slurm_with_singularity.sh
```


- Test workflow with **Conda**:

```bash
bash test/test_workflow_slurm_with_conda.sh
```

#### Rule graph

Execute the following command to generate a rule graph image for the workflow.
The output will be found in the `images/` directory in the repository root. 

```bash
bash test/test_rule_graph.sh
```

You can see the rule graph below in the
[workflow description](#workflow-description) section.

#### Clean up test results

After successfully running the tests above, you can run the following command
to remove all artifacts generated by the test runs:

```bash
bash test/test_cleanup.sh
```

## Usage

Now that your virtual environment is set up and the workflow is deployed and
tested, you can go ahead and run the workflow on your samples.

### Preparing inputs 

It is suggested to have all the input files for a given run (or hard links 
pointing to them) inside a dedicated directory, for instance under the 
_MIRFLOWZ_ root directory. This way it is easier to keep the data together, 
reproduce an analysis and set up Singularity access to them.  

#### 1. Prepare a sample table

Refer to `test/test_files/sample_table.tsv` to know what this file 
must look like, or use it as a template. 

```bash
touch path/to/your/sample/table.tsv
```
> Fill the sample table according to the following requirements:  
>
> - `sample`. This column contains the library name.  
> - `sample_file`. In this column, you must provide the path to the library file.
> The path must be relative to the working directory.  
> - `adapter`.  This field must contain the adapter sequence in capital letters.  
> - `format`. In this field you must state the library format. It can either be 
> `fa` if providing a FASTA file or `fastq` if the library is a FASTQ file.  
> 

#### 2. Prepare genome resources

There are 4 files you must provide: 

1. A **`gzip`ped FASTA** file containing **reference sequences**, typically the
   genome of the source/organism from which the library was extracted.

2. A **`gzip`ped GTF** file with matching **gene annotations** for the
   reference sequences above.

> _MIRFLOWZ_ expects both the reference sequence and gene annotation files to
> follow [Ensembl][ensembl] style/formatting. If you obtained these files from
> a source other than Ensembl, you may first need to convert them to the
> expected style to avoid issues!

3. An **uncompressed GFF3** file with **microRNA annotations** for the reference
   sequences above.

> _MIRFLOWZ_ expects the miRNA annotations to follow [miRBase][mirbase]
> style/formatting. If you obtained this file from a source other than miRBase,
> you may first need to convert it to the expected style to avoid issues!

4. An **uncompressed tab-separated file** with a **mapping between the
   reference names** used in the miRNA annotation file (column 1; "UCSC style")
   and in the gene annotations and reference sequence files (column 2; "Ensembl
   style"). Values in column 1 are expected to be unique, no header is
   expected, and any additional columns will be ignored. [This
   resource][chrMap] provides such files for various organisms, and in the
   expected format.

> General note: If you want to process the genome resources before use (e.g.,
> filtering), you can do that, but make sure the formats of any modified
> resource files meet the formatting expectations outlined above!

#### 3. Prepare a configuration file

We recommend creating a copy of the
[configuration file template](config/config_template.yaml):

```bash
cp  config/config_template.yaml  path/to/config.yaml
```

Open the new copy in your editor of choice and adjust the configuration
parameters to your liking. The template explains what each of the
parameters means and how you can meaningfully adjust them. 

### Running the workflow

With all the required files in place, you can now run the workflow locally
via Singularity with the following command:  

```bash
snakemake \
    --snakefile="path/to/Snakefile" \
    --cores 4  \
    --configfile="path/to/config.yaml" \
    --use-singularity \
    --singularity-args "--bind ${PWD}/../" \
    --printshellcmds \
    --rerun-incomplete \
    --verbose
```

> **NOTE:** Depending on your working directory, you do not need to use the 
> parameters  `--snakefile` and `--configfile`. For instance, if the `Snakefile`
> is in the same directory or the `workflow/` directory is beneath the current
> working directory, there's no need for the `--snakefile` directory. Refer to 
> the [Snakemake documentation][snakemakeDocu] for more information.

After successful execution of the workflow, results and logs will be found in
the `results/` and `logs/` directories, respectively.

### Creating a Snakemake report

Snakemake provides the option to generate a detailed HTML report on runtime
statistics, workflow topology and results. If you want to create a Snakemake
report, you must run the following command:

```bash
snakemake \
    --snakefile="path/to/Snakefile" \
    --configfile="path/to/config.yaml" \
    --report="snakemake_report.html"
```

> **NOTE:** The report creation must be done after running the workflow in
> order to have the runtime statistics and the results. 

## Workflow description

The _MIRFLOWZ_ workflow first processes and indexes the user-provided genome 
resources. Afterwards, the user-provided short read small-RNA-seq libraries will
be aligned separately against the genome and transcriptome. For increased 
fidelity, two separated aligners, [Segemehl][segemehl] and our in-house tool 
[Oligomap][oligomap], are used. All the resulting alignments are merged such 
that only the best alignments of each read are kept (smallest edit distance).
Finally, alignments are intersected with the user-provided, pre-processed
miRNA annotation file using [BEDTools][bedtools]. Counts are tabulated 
separately for reads consistent with either miRNA precursors, mature miRNA
and/or isomiRs.

> **NOTE:** For a detailed description of each rule, please, refer to the
> [workflow documentation](pipeline_documentation.md)

The schema below is a visual representation of the individual workflow steps
and how they are related:

> ![rule-graph][rule-graph]

## Contributing

_MIRFLOWZ_ is an open-source project which relies on community contributions.
You are welcome to participate by submitting bug reports or feature requests,
taking part in discussions, or proposing fixes and other code changes. Please
refer to the [contributing guidelines](CONTRIBUTING.md) if you are interested in
contribute.

## License

This project is covered by the [MIT License](LICENSE).

## Contact

For questions or suggestions regarding the code, please use the [issue tracker][issue-tracker]. Do not hesitate on contacting us via [email][email] for any other inquiries.

&copy; 2023 [Zavolab, Biozentrum, University of Basel][zavolab]

[bedtools]: <https://github.com/arq5x/bedtools2>
[chrMap]: <https://github.com/dpryan79/ChromosomeMappings>
[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[cluster execution]: <https://snakemake.readthedocs.io/en/stable/executing/cluster.html>
[email]: <zavolab-biozentrum@unibas.ch>
[ensembl]: <https://ensembl.org/>
[issue-tracker]: <https://github.com/zavolanlab/mirflowz/issues>
[mamba]: <https://github.com/mamba-org/mamba>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[mirbase]: <https://mirbase.org/>
[oligomap]: <https://bio.tools/oligomap>
[rule-graph]: images/rule_graph.svg
[segemehl]: <https://www.bioinf.uni-leipzig.de/Software/segemehl/>
[singularity]: <https://apptainer.org/admin-docs/3.8/index.html>
[slurm]: <https://slurm.schedmd.com/documentation.html>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[snakemakeDocu]: <https://snakemake.readthedocs.io/en/stable/executing/cli.html>
[zavolab]: <https://www.biozentrum.unibas.ch/research/researchgroups/overview/unit/zavolan/research-group-mihaela-zavolan/>
