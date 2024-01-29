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
    - [Expected output files](#expected-output-files)
    - [Creating a Snakemake report](#creating-a-snakemake-report)
3. [Workflow description](#workflow-description)
4. [Contributing](#contributing)
5. [License](#license)
6. [Contact](#contact)

## Installation

The workflow lives inside this repository and will be available for you to run
after following the installation instructions laid out in this section.

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

If you plan to run _MIRFLOWZ_ via Conda, we recommend using the following
command for a faster environment creation, specially if you will run it on an
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
_MIRFLOWZ_ root directory. This way, it is easier to keep the data together,
reproduce analysis and set up Singularity access to them.  

#### 1. Prepare a sample table

Refer to `test/test_files/sample_table.tsv` to know what this file 
must look like, or use it as a template. 

```bash
touch path/to/your/sample/table.tsv
```
> Fill the sample table according to the following requirements:  
>
> - `sample`. Arbitrary name for the miRNA sequencing library.
> - `sample_file`. Path to the miRNA sequencing library file. The path must be
> relative to the directory where the workflow will be run.
> - `adapter`. Sequence of the 3'-end adapter used during library preparation.
> - `format`. One of `fa`/`fasta` or `fq`/`fastq`, if the library file is in
> FASTA or FASTQ format, respectively.

#### 2. Prepare genome resources

There are 4 files you must provide: 

1. A **`gzip`ped FASTA** file containing **reference sequences**, typically the
   genome of the source/organism from which the library was extracted.

2. A **`gzip`ped GTF** file with matching **gene annotations** for the
   reference sequences above.

> _MIRFLOWZ_ expects both the reference sequence and gene annotation files to
> follow [Ensembl][ensembl] style/formatting. If you obtained these files from
> a source other than Ensembl, you must ensure that they adhere to the
> expected format by converting them, if necessary.

3. An **uncompressed GFF3** file with **microRNA annotations** for the reference
   sequences above.

> _MIRFLOWZ_ expects the miRNA annotations to follow [miRBase][mirbase]
> style/formatting. If you obtained this file from a source other than miRBase,
> you must ensure that it adheres to the expected format by converting it, if
> necessary.

4. An **uncompressed tab-separated file** with a **mapping between the
   reference names** used in the miRNA annotation file (column 1; "UCSC style")
   and in the gene annotations and reference sequence files (column 2; "Ensembl
   style"). Values in column 1 are expected to be unique, no header is
   expected, and any additional columns will be ignored. [This
   resource][chrMap] provides such files for various organisms, and in the
   expected format.

5. **OPTIONAL**: A **BED6** file with regions for which to produce
   [ASCII-style alignment pileups][ascii-pileups]. If not provided, no pileups
   will be generated. See [here][bed-format] for the expected format.

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
parameters mean and how you can meaningfully adjust them. 

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
> parameters `--snakefile` and `--configfile`. For instance, if the `Snakefile`
> is in the same directory or the `workflow/` directory is beneath the current
> working directory, there's no need for the `--snakefile` directory. Refer to 
> the [Snakemake documentation][snakemakeDocu] for more information.

After successful execution of the workflow, results and logs will be found in
the `results/` and `logs/` directories, respectively.

### Expected output files

Upon successful execution of _MIRFLOWZ_, the tool automatically removes all
intermediate files generated during the process. The final outputs comprise:

1. A SAM file containing alignments intersecting a pri-miR locus. These
alignments intersect with extended start and/or end positions specified in the
provided pri-miR annotations. Please note that they may not contribute to the
final counting and may not appear in the final table.

2. A SAM file containing alignments intersecting a mature miRNA locus. Similar
to the previous file, these alignments intersect with extended start and/or end
positions specified in the provided miRNA annotations. They may not contribute
to the final counting and might be absent from the final table.

3. A BAM file containing the set of alignments contributing to the final
counting and its corresponding index file (`.bam.bai`).

4. Table(s) containing the counting data from all libraries for (iso)miRs
and/or pri-miRs. Each row corresponds to a miRNA species, and each column
represents a sample library. Each read is counted towards all the annotated
miRNA species it aligns to, with 1/n, where n is the number of genomic and/or
transcriptomic loci that read aligns to.

5. **OPTIONAL**. ASCII-style pileups of read alignments produced for individual
libraries, combinations of libraries and/or all libraries of a given run. The
exact number and nature of the outputs depends on the workflow
inputs/parameters. See the
[pileups section](pipeline_documentation.md/#pileup-workflow) for a detailed
description.

To retain all intermediate files, include `--no-hooks` in the workflow call.

```bash
snakemake \
    --snakefile="path/to/Snakefile" \
    --cores 4  \
    --configfile="path/to/config.yaml" \
    --use-conda \
    --printshellcmds \
    --rerun-incomplete \
    --no-hooks \
    --verbose
```

After successful execution of the workflow, the intermediate files will be
found in the `results/intermediates` directory.

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

The _MIRFLOWZ_ workflow initially processes and indexes the genome resources
provided by the user. The regions corresponding to mature miRNAs are extended
on both sides to accommodate isomiR species with shifted start and/or end
positions. If necessary, pri-miR loci are similarly extended to adjust to the
new miRNA coordinates.

Subsequently, the user-provided short-read small RNA-seq libraries undergo
quality filtering if a FASTQ file is provided. Alternatively, adapters are
directly removed. The resulting reads are independently mapped to both the
genome and the transcriptome using two distinct aligners: [Segemehl][segemehl]
and our in-house tool [Oligomap][oligomap]. After the mapping, only the best
alignments for each read, determined by the smallest edit distance, are
retained by merging and filtering the resulting alignments into a single file.

The collection of resulting alignments is then reduced to contain only unique
entries. Due to the short length of the reads and the sequence similarity among
miRNAs, the number of alignments can be high. Therefore, reads aligned beyond a
specified threshold are discarded. To address multimapping, alignments with the
fewest indels are preserved. These alignments are subsequently intersected with
the user-provided, pre-processed miRNA annotation files using
[BEDTools][bedtools]. Note that an alignment will not contribute to the final
count if its start and/or end positions differ significantly from the provided
miRNA annotations, beyond the extension applied to the mature miRNA start
and/or end positions, or by 1 if no extension was applied. Conversely, a
retained read contributes 1/n to all the annotated miRNA species it aligns
with, where `n` is the number of genomic and/or transcriptomic loci it aligns
to.

_MIRFLOWZ_ employs an unambiguous notation to classify isomiRs using the format
`miRNA_name|5p-shift|3p-shift|CIGAR|MD`, where `5p-shift` and `3p-shift`
represent the differences between the annotated mature miRNA start and end
positions and those of the alignment, respectively.

Counts are tabulated separately for reads consistent with either
miRNA precursors, mature miRNA and/or isomiRs and all library counts are
fused into a single table. Finally, ASCII-style alignment pileups are
optionally generated for user-defined regions of interest.

> **NOTE:** For a detailed description of each rule along with some examples,
> please, refer to the [workflow documentation](pipeline_documentation.md).

The schema below is a visual representation of the individual workflow steps
and how they are related:

> ![rule-graph][rule-graph]

## Contributing

_MIRFLOWZ_ is an open-source project which relies on community contributions.
You are welcome to participate by submitting bug reports or feature requests,
taking part in discussions, or proposing fixes and other code changes. Please
refer to the [contributing guidelines](CONTRIBUTING.md) if you are interested
in contribute.

## License

This project is covered by the [MIT License](LICENSE).

## Contact

For questions or suggestions regarding the code, please use the
[issue tracker][issue-tracker]. Do not hesitate on contacting us via
[email][email] for any other inquiries.

&copy; 2023 [Zavolab, Biozentrum, University of Basel][zavolab]

[ascii-pileups]: <https://git.scicore.unibas.ch/zavolan_group/tools/ascii-alignment-pileup>
[bed-format]: <https://gist.github.com/deliaBlue/19ad3740c95937378bd9281bd9d1bc72>
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
