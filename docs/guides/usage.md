# Usage

Learn here how to run _MIRFLOWZ_.

!!! info "Prerequisites"

    This usage example assumes that you have already
    [installed](./installation.md) _MIRFLOWZ_.

## How to analyze short-read small RNA-Seq samples?

Assuming that your current directory is the workflow repository's root
directory, create a directory for your workflow run and traverse into it with:

```sh
mkdir my_run
cd my_run
```

### Preparing input files

It is suggested to have all the input files for a given run (or hard links
pointing to them) inside a dedicated directory, for instance under the
`input_files/` subdirectoy in `my_run/` directory. This way, it is easier to
keep the data together, set up Apptainer access to them and reproduce analyses.

Create this directory and traverse into it with:

```sh
mkdir input_files
cd input_files
```

#### 1. Prepare a sample table

Create an empty sample table:

```sh
touch samples.tsv
```

Use your editor of choice to populate the sample table first with a header,
and then according to the following requirements:

- `sample` Arbitrary name for the miRNA sequencing library
- `sample_file` Path to the miRNA sequencing library file. The path must be
  relative to the directory where the workflow is going to be run
- `adapter` 3'-end adapter sequence used during library preparation
- `format` One of `fa`/`fasta` or `fq`/`fastq`, if the library file is in
  FASTA or FASTQ format respectively.

??? tip "How can I be sure the samples table has the correct format?"

    You may refer to the [test sample table][sample_tsv] to know what the
    samples table must look like, or use it as a template.

#### 2. Prepare the genome resources

There are 4 files you must provide:

1. A `gzip`**ed FASTA** file containing **reference sequences**, typically the
   genome of the source/organism from which the library was extracted.

2. A `gzip`**ed GTF** file with matching **gene annotations** for the
   reference sequences above.

> _MIRFLOWZ_ expects both, the reference sequence and gene annotation files to
> follow [Ensembl][ensembl] style/formatting. If you obtained these files from
> a source other than Ensembl, you must ensure that they adhere to the
> expected format by converting them, if necessary.

3. An **uncompressed GFF3** file with **microRNA annotations** for the
   reference sequences above.

> _MIRFLOWZ_ expects the miRNA annotations to follow [miRBase][mirbase]
> style/formatting. If you obtained this file from a source other than miRBase,
> you must ensure that it adheres to the expected format by converting it, if
> necessary.

4. An **uncompressed tab-separated file** with a **mapping between the**
   **reference names** used in the miRNA annotation file (column 1; "UCSC
   style") and in the gene annotations and reference sequence files (column
   2; "Ensembl style"). Values in column 1 are expected to be unique, no header
   is expected, and any additional columns will be ignored.
   [This resrouce][chr_map] provides such files for various organisms, and in
   the expected format.

5. **OPTTIONAL**: A **BED6** file with the regions for which to produce
   [ASCII-style pileups][pileups]. If not provided, no pileups are generated.
   See [here][bed-format] for the expected format.

??? tip "Can I process the genome resources before use?"

    Yes, any input file can be processed (_e.g._, filtering) before use as
    long as you make sure the format of any modified resource file meet the
    formatting expectations outlines above!

#### 3. Prepare a configuration file

Return to your working directory (`my_run/`) to create the configuration file.
We recommend creating a copy of the
[configuration file template][mirflowz_conf]:

```sh
cp ../config/config_template.yaml config.yaml
```

Open the new copy in your editor of choice and adjust the configuration
parameters to your liking. The template explains what each of the parameters
mean and how you can meaningfully adjust them.

!!! info "Cluster configuration"

    As done with the configuration file, the [cluster JSON file][cluster_json]
    can be copied into your working directory. The default values have to be
    modified to meet your cluster specifications.

### Running the workflow

With all the required files in place, you can now run the workflow locally
within an activated `mirflowz` environment with the following call:

```sh
snakemake \
    --snakefile="path/to/Snakefile" \
    --cores 4  \
    --configfile="path/to/config.yaml" \
    --software-deployment-method conda \
    --printshellcmds \
    --rerun-incomplete \
    --verbose
```

??? question "Want to use Apptainer instead?"

    Change the argument to `--software-deployment-method` from `conda` to
    `apptainer` and add `--apptainer-args "--bind ${PWD}/../"` to execute the
    workflow steps within Apptainer containers.

After successful execution of the workflow, results and logs will be found in
the `results/` and `logs/` directories, respectively.
