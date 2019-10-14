# mir-prepare-annotation

Downloads and prepares the necessary files for smallRNA-seq related pipelines.

## Installation

Download and install Miniconda 3

Linux users
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Mac users
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

Clone repository
```bash
git clone -b dev https://git.scicore.unibas.ch/AnnotationPipelines/mir-prepare-annotation.git
```

Create virtual environment
```bash
conda create --name mir -c bioconda -c conda-forge snakemake=5.4.0
```

Activate virtual environment
```bash
conda activate mir
```

Deactivate virtual environment
```bash
conda deactivate 
```

## Usage

* Main configuration file: `config.yaml`. Contains all the paths needed for successful snakemake workflow execution as well as all the parameters for a proper set up for the posterior miRNA sequences analysis. Please go carefully through all the entries and adjust all of them according to your desires.

For workflow execution:
```bash
./run_pipeline.sh
```

For creating DAG (Directed Acyclic Graph, collection of all the tasks you want to run, organized in a way that reflects their relationships and dependencies):
```bash
./create_dag.sh
```