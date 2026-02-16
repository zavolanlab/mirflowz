# Installation

Find here how to install _MIRFLOWZ_ on your system.

## Requirements/Dependencies

For improved reproducibility and reusability, as well as an easy means to run
the workflow on a high performance computing (HPC) cluster managed, _e.g._, by
[Slurm][slurm], each individual workflow step runs inside an isolated
[Conda][conda] environment or an [Apptainer][apptainer] container. As a
consequence, running _MIRFLOWZ_ only has a few individual dependencies, that
need to be available on your system:

* Linux (tested with `Ubuntu 24.04` and `Arch Linux`)
* [Conda][conda] (tested with `conda 25.3.0`)
* **(Optional)** [Mamba][mamba] (tested with `mamba 2.1.1`)
* **(Optional)** [Apptainer][apptainer] (tested with `apptainer 1.4.0`)

A few additional dependencies are installed via Conda as described further
below.

!!! warning "Other versions, especially older ones, are not guaranteed to work!"

??? question "How do I install Apptainer?"

    Please follow the [official documentation][apptainer-docs] to install
    Apptainer (formerly Singularity) globally and configure its permissions.

??? tip "You want to run _MIRFLOWZ_ on an HPC?"

    When running _MIRFLOWZ_ on a High-Performance Computing cluster, you will
    need to make sure that compatible versions of at least one of Conda (when
    using Snakemake's `--use-conda` option) and Apptainer (when using the
    `--use-apptainer` option) are installed and properly configured on each
    machine of the cluster, including the head node. Reach out to your systems
    administrator to set this up for you.

## Installation Steps

### 1. Clone _MIRFLOWZ_

Clone the [_MIRFLOWZ_ workflow repository][mirflowz] with:

```sh
git clone https://github.com/zavolanlab/mirflowz.git
```

### 2. Ensure Conda is available

To check if you already have Conda installed on your system, type:

```sh
conda --version
```

If Conda is available, you should see an output similar to this:

```sh
conda 25.3.0
```

If it is not installed, you will instead see
<code style="color: red;">command not found: conda</code>.

??? question "What's the best way to install Conda?"

    Conda can be installed in multiple ways. We strongly recommend using the
    [Miniforge][miniforge] distribution, as it is built around the
    community-supplied `conda-forge` channel that is heavily used by _MIRFLOWZ_.
    It also reduces the risk of accidentally violating Conda's licensing
    restrictions by using the `defaults` channel.

    Please refer to the [official documentation][miniforge] for up-to-date
    installation instructions.

??? question "Do I need to install Mamba if I already have Conda?"

    Short answer: no, you don't need to install Mamba on top of Conda but it
    is recommended.

    Longer answer: You can run _MIRFLOWZ_ only with Conda. But, for faster
    creation of the environment (and Conda environments in general), you can
    install Mamba. Note that if you had installed Conda using
    Miniforge you already have Mamba installed.

    If you decide to use Mamba, replace `conda` with `mamba` in the commands
    below (particularly in `conda env create`).

??? tip "My Conda version is not compatible"

    After completing Conda setup, you can install a specific Conda version with
    the following command:

    ```sh
    conda install conda=25.3.0
    ```

??? tip "I do not want to change my Conda version"

    If you already have a specific Conda version on your system that is not
    compatible with _MIRFLOWZ_, and you do not want to change it, no worries.
    Just indicate a different directory when using the interactive Miniforge
    installer. Then source the appropriate `conda.sh` file to switch between
    versions, e.g.:

    ```sh
    source $HOME/miniconda3/etc/profile.d/conda.sh  # OR
    source $HOME/miniconda3_alt/etc/profile.d/conda.sh
    ```

### 3. Set up your _MIRFLOWZ_ environment

In this step, you need to install the remaining workflow dependencies. For
this purpose, there exist two different "environment files". Use this decision
matrix to pick the most suitable one for you:

| I want to run pre-packaged tests &emsp;&emsp;&emsp;&emsp;| Environment file to use &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; |
|:---:| --- |
| | `environment.yml` |
| :check_mark: | `environment.dev.yml` |


To set up the environment, execute the call below, but do not forget to replace
the placeholder `ENVIRONMENT` with the appropriate file from the table above:

```sh
conda env create -f ENVIRONEMTN
```

### 4. Activate the _MIRFLOWZ_ environment

Activate the Conda environment with:

```sh
conda activate mirflowz
```

??? tip "Planning to run _MIRFLOWZ_ via Conda?"

    If you plan to run _MIRFLOWZ_ via Conda, we recommend using the following
    command for a faster environmetn creation, specially if you will run it
    on an HPC cluster.

    ```sh
    conda config --set channel_priority strict
    ```

## Running Installation Tests

We have prepared several tests to check the integrity of the pipeline and its
components. These can be found in the `test/` directory. The most critical of
these tests enable you to execute the entire workflow on a set of small
example input files. Note that for some of these and other tests to complete
successfully, additional
[development dependencies](#3-set-up-your-mirflowz-environment) need to be
installed.

Execute the command below to run the workflow test on your local machine or
on your Slurm-managed HPC cluster.

??? tip "Failing tests do not necessarily indicate a problem!"

    Our tests were developed to guard against code regression over time or as
    a result of proposed changes and are therefore very rigorous. In
    particular, even the minutest changes in outputs, even individual pixels
    or metadata in the produced output images, will cause a test suite to
    fail. However, we cannot rule out such minor changes across systems, and
    therefore, it is quite possible that your test runs may fail, even though
    the workflow is properly installed and functional. Check the logs to see
    whether Snakemake completed successfully. If so, it is very likely that
    everything is fine.

### Running tests on your local machine

Execute one of the following commands to run the test workflow on your local
machine:

- Test workflow running each step inside an **Apptainer container**:

```sh
bash test/test_integration_workflow/test_workflow_local_with_apptainer.sh
```

- Test workflow running each step inside a dedicated **Conda environment**:

```sh
bash test/test_integration_workflow/test_workflow_local_with_conda.sh
```

### Running tests on your Slurm cluster

Execute one of the following commands to run the test workflow on a
Slurm-managed HPC cluster:

- Test workflow running each step inside an **Apptainer container**:

```sh
bash test/test_integration_workflow/test_workflow_slurm_with_apptainer.sh
```

- Test workflow running each step inside a dedicated **Conda environment**:

```sh
bash test/test_integration_workflow/test_workflow_slurm_with_conda.sh
```

??? tip "The Slurm tests are failing for me!"

    Depending on the configuration of your Slurm installation you may need to
    adapt the `--cluster`argument's option `-p` in the Slurm test files and the
    options `queue`, `time`, `threads` and `mem` in the
    [cluster.json][cluster_json] file. Consult the manual of your workload
    manager.

### Clean up test results

After successfully running the tests above, you can run the following command
to remove all artifacts generated by the test runs:

```sh
bash test/test_cleanup.sh
```
