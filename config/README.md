# Dependencies installation

Create and activate the virtual environment with the required dependencies
with Conda:

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

For a faster creation of the environment (and Conda environments in general),
you can also install [Mamba][mamba] on top of Conda. In that case, replace
`conda` with `mamba` in the commands above (particularly in 
`conda env create`).

## Running _MIRFLOWZ_ with Singularity

If you want to run _MIRFLOWZ_ via Singularity and do not already
have it installed globally on your system, you must further update the Conda
environment with:

```bash
conda env update -f environment.root.yml
```

> Mind that you must have the environment activated and root permissions on
> your system to install Singularity. If you want to run _MIRFLOWZ_ on an HPC
> cluster (recommended in almost all cases), ask your system administrator
> about Singularity.

# Run the workflow on your own samples

In order to run _MIRFLOWZ_ on your own samples, we recommend having all the
input files inside a dedicated directory. This way, it is easier to keep the
data together and reproduce an analysis. Assuming that your current directory
is the repository's root directory, create a directory to store all your data
and traverse to it with:

```bash
mkdir path/to/your_run
cd path/to/your_run
```

## 1. Prepare the sample table

Create an empty sample table. Refer to the
[sample.tsv](../test/test_files/samples_table.tsv) test file to see what the
table must look like or use it as a template.

```bash
touch samples.tsv
```

> Fill the sample table according to the following requirements:  
>
> - `sample`. Arbitrary name for the miRNA sequencing library.
> - `sample_file`. Path to the miRNA sequencing library file. The path must be
> relative to the directory where the workflow will be run.
> - `adapter`. Sequence of the 3'-end adapter used during library preparation.
> - `format`. One of `fa`/`fasta` or `fq`/`fastq`, if the library file is in
> FASTA or FASTQ format, respectively.

## 2. Prepare the genome resources

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
   [ASCII-style alignment pileups][ascii-pileups] for. If not provided, no
   pileups will be generated. Refer [here][bed-format] to see the expected
   format.

> General note: If you want to process the genome resources before use (e.g.,
> filtering), you can do that, but make sure the formats of any modified
> resource files meet the formatting expectations outlined above!

## 3. Prepare the configuration file

We recommend creating a copy of the
[configuration file template](config_template.yaml).

```bash
cp ../config/config_template.yaml config.yaml

```

Open the new copy in your editor of choice and adjust the configuration
parameters to your liking. The template explains what each of the parameters
mean and how you can meaningfully adjust them.

[ascii-pileups]: <https://git.scicore.unibas.ch/zavolan_group/tools/ascii-alignment-pileup>
[bed-format]: <https://gist.github.com/deliaBlue/19ad3740c95937378bd9281bd9d1bc72>
[chrMap]: <https://github.com/dpryan79/ChromosomeMappings>
[ensembl]: <https://ensembl.org/>
[ensembl-bed]: <https://www.ensembl.org/info/website/upload/bed.html>
[mamba]: <https://github.com/mamba-org/mamba>
[mirbase]: <https://mirbase.org/>
