# FixCore

[![Snakemake](https://img.shields.io/badge/Snakemake-≥8.20-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Test workflow](https://github.com/SeviJordi/FixCore/actions/workflows/test.yaml/badge.svg)
 [![License (AGPL version 3)](https://img.shields.io/badge/license-GNU%20AGPL%20version%203-green.svg)](COPYING)
 
Snakemake workflow to polish core genome alignments

---

## Installation

You have the option to download the zipped version of this repository (find the button on the right-hand side of the screen). Also it is possible to clone this repository:

```
git clone git@github.com:SeviJordi/FixCore.git
cd FixCore

```

In case your system lacks Snakemake version 8.20 or newer, it is necessary to [download and install it](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Usage

To use the workflow with default settings modify the [target configuration file](/config/target.yaml). You have to provide the folder in which the fasta files for the individual genes are (aligned or not), the extension of the files and a prefix to use in the outputs. Then you can run the workflow with:

```
snakemake --use-conda -c 8  # To run with 8 threads
```

We recomend to run the workflow with at least 8 threads, but any natural number can be provided. Workflow parameters can be modified in the [workflow configuration file](/config/config.yaml).


## Outputs

Once the workflow is ended, it will have produced a concatenated alignment of the curated alignments of each gene and a maximum likelihood phylogeny for that alignment. Also intermedate files are provided to allow different analyses. 
