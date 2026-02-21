<div align="center">
    <img src=images/mirflowz_logo.png>
</div>

# _MIRFLOWZ_

**Welcome to the _MIRFLOWZ_ documentation pages!**

_MIRFLOWZ_ is a generic small RNA-Seq workflow that allows users to classify
and quantify miRNAs and/or isomiRs.

The workflow is developed in [Snakemake][docs-snakemake], a widely used
workflow management system in the bioinformatic community. _MIRFLOWZ_ will
pre-process, align, classify and quantify your bulk small RNA-seq sequencing
libraries with publicly available state-of-the-art bioinformatics tools.


!!! info "_MIRFLOWZ_ notation"

    _MIRFLOWZ_ uses the notation provided by [miRBase][mirbase] (_i.e._,
    "miRNA primary transcript" for precursors and "miRNA" for the canonical
    mature miRNA). This implies that precursors are named "pri-miRs" across the
    workflow instead of pre-miR. This decision is made upon the lack of
    guarantee that "miRNA primary transcripts" are full pre-miR (and pre-miR
    only) sequences.

## How does it work?

_MIRFLOWZ_ requires [Conda][conda] to install the basic dependencies. Each
individual step of the workflow runs either in its own [Apptainer][apptainer]
container or in its own [Conda][conda] virtual environment.

Once the installation is complete, you fill in a [`config.yaml`][mirflowz_conf]
file with parameters and a [`samples.tsv`][sample_tsv] file with
sample-specific information. You can easily trigger _MIRFLOWZ_ by making a call
to Snakemake with the appropriate parameters.

_MIRFLOWZ_ can be executed in different systems or High Performance Computing
(HPC) cluster.

## How to cite


If you use _MIRFLOWZ_ in your work, please kindly cite the following Zenodo
entry:

**zavolanlab/mirflowz: v0.10.0 (v0.10.0).**
_Iris Mestres-Pascual, Alex Kanitz, & Mihalea Zavolan._
(2025). Zenodo.
[https://doi.org/10.5281/zenodo.17120595](https://doi.org/10.5281/zenodo.17120595)

## Reach out

- For _MIRFLOWZ_ usage questions, please use the
  [_MIRFLOWZ_ Q&A forum][mirflowz_discussions] (requires
  [GitHub registration][github_reg]).

- For feature suggestions and bug reports, please use the
  [_MIRFLOWZ_ issue tracker][mirflowz_issues] (requires
  [GitHub registration][github_reg]).

- For any other requests, please reach out to us via [email][email].

## Contributing

We always welcome and duly acknowledge open source contributors, for _MIRFLOWZ_
or any other of [our projects][zavolab_projects]. Simply follow our
[onboarding instructions][mirflowz_contribute] and please mind our
[Code of Conduct][mirflowz_coc]. If you have any questions, do not hesitate to
shoot us an [email][email].

## Acknowledgements

[![Zavolab](images/zavolab_logo.200px.png)](https://www.biozentrum.unibas.ch/research/research-groups/research-groups-a-z/overview/unit/research-group-mihaela-zavolan)
[![Biozentrum, University of Basel](images/biozentrum_logo.200px.png)](https://www.biozentrum.unibas.ch/)
