[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/tadrep/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tadrep.svg)
![PyPI - Status](https://img.shields.io/pypi/status/tadrep.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/tadrep.svg)
[![PyPI](https://img.shields.io/pypi/v/tadrep.svg)](https://pypi.org/project/tadrep)
[![Conda](https://img.shields.io/conda/v/bioconda/tadrep.svg)](https://bioconda.github.io/recipes/tadrep/README.html)

# TaDReP: Targeted Detection and Reconstruction of Plasmids

TaDReP is a tool for the rapid and targeted detection and reconstruction of plasmids within bacterial draft assemblies via contig alignments.

- [Description](#description)
- [Installation](#installation)
- [Examples](#examples)
- [Input & Output](#input-and-output)
- [Usage](#usage)
- [Database](#database)
- [Issues & Feature Requests](#issues)

## Description

TaDReP facilitates the rapid screening of target plasmids within single genomes or genome cohorts.

It detects and reconstructs reference plasmids within bacterial draft assemblies via contig alignments. Therefore, contigs from draft assemblies are aligned against reference plasmid sequences via BLAST+ and rigourously filtered for contig-wise coverage and identity thresholds. Reference plasmids are finally detected and reconstructed upon strict plasmid-wise coverage and identity thresholds.

## Installation

TaDRep can be installed via BioConda and Pip. However, we encourage to use [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to automatically install all required 3rd party dependencies.

### BioConda

```bash
$ conda install -c conda-forge -c bioconda -c defaults tadrep
```

### Pip

```bash
$ python3 -m pip install --user tadrep
```

TaDRep requires [Blast+](https://blast.ncbi.nlm.nih.gov) which must be installed & executable.

## Examples

Simple:

```bash
$ tadrep --genome test/data/draft.fna --plasmids test/data/plasmids.fna
```

Expert: verbose output writing results to *results* directory using 8 threads, a min per contig sequence identity of 75% and a gap sequence length (multiple `N`) of 100.:

```bash
$ tadrep --genome test/data/draft.fna --db refseq -v --min-contig-identity 75 --gap-sequence-length 100 --output results --threads THREADS 8
```

## Input and Output

### Input

TadReP accepts bacterial draft genome assemblies in (zipped) fasta format. Complete reference plasmid sequences are provided either as a multi Fasta file or a custom database. For the latter, please read the [database](#database) section below.

### Output

Per draft genome
For each draft genome TaDReP writes a TSV summary file providing all detected reference plasmids and aligned genome contigs. For each reference plasmid that was detected in a draft assembly, ordered and rearranged contigs are exported as mere contigs and as a single pseudomolecule sequence combining contigs separated by multiple `N`.
Furthermore, for each reconstructed plasmid, the reference plasmid backbone and all contig alignments are visualized as `PDF`.

- `<genome>-summary.tsv`: detailed per contig alignment summary
- `<genome>-<plasmid>-contigs.fna`: ordered and rearranged contigs of the reconstructed plasmid
- `<genome>-<plasmid>-pseudo.fna`: pseudomolecule sequence of the reconstructed plasmid
- `<genome>-<plasmid>.pdf`: visualization of aligned contigs against the detected reference plasmid

If multiple genomes were provided, TaDReP also provides a presence/absence matrix of all detected plasmids for convenient cohort analyses.

## Usage

```bash
usage: TaDReP [--genome GENOME [GENOME ...]] [--plasmids PLASMIDS] [--db DB]
                [--output OUTPUT] [--prefix PREFIX]
                [--min-contig-coverage [1-100]] [--min-contig-identity [1-100]]
                [--min-plasmid-coverage [1-100]] [--min-plasmid-identity [1-100]]
                [--gap-sequence-length GAP_SEQUENCE_LENGTH]
                [--help] [--verbose] [--threads THREADS] [--tmp-dir TMP_DIR] [--version]

Targeted Detection and Reconstruction of Plasmids

Input / Output:
  --genome GENOME [GENOME ...], -g GENOME [GENOME ...]
                        Draft genome path
  --plasmids PLASMIDS, -p PLASMIDS
                        Plasmids path
  --db DB               Directory which contains blast database
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)
  --prefix PREFIX       Prefix for all output files (default = None)

Annotation:
  --min-contig-coverage [1-100]
                        Minimal contig coverage (default = 90%)
  --min-contig-identity [1-100]
                        Maximal contig identity (default = 90%)
  --min-plasmid-coverage [1-100]
                        Minimal plasmid coverage (default = 80%)
  --min-plasmid-identity [1-100]
                        Minimal plasmid identity (default = 90%)
  --gap-sequence-length GAP_SEQUENCE_LENGTH
                        Gap sequence N length (default = 10)

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of available CPUs)
  --tmp-dir TMP_DIR     Temporary directory to store blast hits
  --version             show program's version number and exit
```

## Database


## Issues & Feature Requests

TaDReP is brand new and like in every software, expect some bugs lurking around. So, if you run into any issues with TaDReP, we'd be happy to hear about it.
Therefore, please, execute it in verbose mode (`-v`) and do not hesitate to file an issue including as much information as possible:

- a detailed description of the issue
- command line output
- log file (`tadrep.log`)
- a reproducible example of the issue with an input file that you can share _if possible_
