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

In case your system lacks Snakemake version 8.20.5, it is necessary to [download and install it](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

If ypu have a conda distribution installed, you can use it to install snakemake as follows:

```
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.20.5
conda activate snakemake
```

## Usage

To use the workflow with default settings modify the [target configuration file](/config/target.yaml) to set the input. You can run the pieline on the output of a core gene analysis by providing the path to a folder with a set of multifasta files for gene families with **.fasta** extension.  Also you can provide the path to a folder containing assemblies with **.fasta** extension and set the core tool to "panacota", "roary" or "panaroo" in the [config file](/config/config.yaml). In that case the pipeline will previusly run a pangenome analysis on the genomes to get the core genes and then aply the coorection to those genes. All the parameters used for finding the core genome can be modifyed in the [config file](/config/config.yaml). Then you can run the workflow with:

```
snakemake --use-conda -c 8  # To run with 8 threads
```

We recomend to run the workflow with at least 8 threads, but any natural number can be provided. Workflow parameters can be modified in the [workflow configuration file](/config/config.yaml).


## Outputs

Once the workflow is ended, it will have produced a concatenated alignment of the curated alignments of each gene and a maximum likelihood phylogeny for that alignment. Also intermedate files are provided to allow different analyses. 


## Contributing
Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature-name`).
3. Commit your changes (`git commit -m 'Add new feature'`).
4. Push to the branch (`git push origin feature-name`).
5. Open a Pull Request.

---

## License
FixCore is licensed under the **GNU AGPL v3**. See the [LICENSE](LICENSE) file for details.

## Contact
For issues or questions, open an [issue](https://github.com/SeviJordi/FixCore/issues) on GitHub.

