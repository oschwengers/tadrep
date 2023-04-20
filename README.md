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
  * [Setup](#setup)
  * [Extract](#extract)
  * [Characterize](#characterize)
  * [Cluster](#cluster)
  * [Detect](#detect)
  * [Visualize](#visualize)
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

TaDReP accepts bacterial draft genome assemblies in (zipped) fasta format. Complete reference plasmid sequences are provided either as a multi Fasta file or a custom database. For the latter, please read the [database](#database) section below.

### Output

Per draft genome
For each draft genome TaDReP writes a TSV summary file providing all detected reference plasmids and aligned genome contigs. For each reference plasmid that was detected in a draft assembly, ordered and rearranged contigs are exported as mere contigs and as a single pseudomolecule sequence combining contigs separated by multiple `N`.
Furthermore, for each reconstructed plasmid, the reference plasmid backbone and all contig alignments are visualized as `PDF`.

- `<genome>-summary.tsv`: detailed per contig alignment summary
- `<genome>-<plasmid>-contigs.fna`: ordered and rearranged contigs of the reconstructed plasmid
- `<genome>-<plasmid>-pseudo.fna`: pseudomolecule sequence of the reconstructed plasmid
- `<genome>-<plasmid>.pdf`: visualization of aligned contigs against the detected reference plasmid

If multiple genomes were provided, TaDReP also provides a presence/absence matrix of all detected plasmids for convenient cohort analyses.

---
# Usage

TaDReP is split up into six different submodules to provide easier usage.

```bash
usage: TaDReP [--help] [--verbose] [--threads THREADS] [--tmp-dir TMP_DIR] [--version] [--output OUTPUT] [--prefix PREFIX] {setup,extract,characterize,cluster,detect,visualize} ...

Targeted Detection and Reconstruction of Plasmids

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of available CPUs)
  --tmp-dir TMP_DIR     Temporary directory to store blast hits
  --version             "show program's version number and exit"

General Input / Output:
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)
  --prefix PREFIX       Prefix for all output files (default = None)

Modules:
  {setup,extract,characterize,cluster,detect,visualize}
    setup               Download and prepare inc-types
    extract             Extract unique plasmid sequences
    characterize        Identify plasmids with GC content, Inc types, conjugation genes
    cluster             Cluster related plasmids
    detect              Detect and reconstruct plasmids in draft genomes
    visualize           Visualize plasmid coverage of contigs

Citation:
Schwengers et al. (2022)
TaDReP: Targeted Detection and Reconstruction of Plasmids.
GitHub https://github.com/oschwengers/tadrep
```

---
## Setup
---

The setup module ist currently used to download and modify plasmid incompatibility groups, required to characterize plasmids.

### Example
Verbosely download and write inc-types into folder 'inc-types':
```bash
tadrep -v -o inc-types setup
```

---
## Extract
---

Extract reference plasmid sequences from complete genomes, draft genomes or plasmid files in fasta format.
Type 'genome' extracts all but the longest sequence, which can be adjusted by setting '--discard_longest'.
Type 'draft' in combination with '--header' provides a possibility to extract only sequences with specific header contents.
Type 'plasmid' extracts all sequences from a given file without any filtering.

```bash
usage: TaDReP extract [-h] [--type {genome,plasmid,draft}] [--header HEADER] [--files FILES [FILES ...]] [--discard-longest DISCARD_LONGEST]

optional arguments:
  -h, --help            show this help message and exit

Input:
  --type {genome,plasmid,draft}, -t {genome,plasmid,draft}
                        Type of input files
  --header HEADER       Template for header description inside input files: e.g.: header: ">pl1234" --> --header "pl"
  --files FILES [FILES ...], -f FILES [FILES ...]
                        File path
  --discard-longest DISCARD_LONGEST, -d DISCARD_LONGEST
                        Discard n longest sequences in output
```

### Examples

Extract all sequences from file `plasmids.fna` into folder `showcase` but ignore the two longest.
```bash
tadrep -v -o showcase extract --type genome --discard-longest 2 --files plasmids.fna
```

Extract all sequences from file `plasmids.fna` where `header` contains `pl` into folder `showcase`.
```bash
tadrep -v -o showcase extract --type draft --header "pl" --files plasmids.fna
```

Extract all sequences from file `plasmids.fna` into folder `showcase`.
```bash
tadrep -v -o showcase extract --type plasmid --files plasmids.fna
```

---
## Characterize
---

All plasmids are characterized by following features:

- Length
- GC - content
- Incompatibility types
- coding sequences

This module provides means to copy reference databases from different folders and to import incompatibility types downloaded by [tadrep setup](#setup).

If you downloaded a reference database this is the step to start with.

```bash
usage: TaDReP characterize [-h] [--db DATABASE] [--inc-types INC_TYPES]

optional arguments:
  -h, --help            show this help message and exit

Input:
  --db DATABASE         Import json file from a given database path into working directory
  --inc-types INC_TYPES
                        Import inc-types from given path into working directory
```

### Examples

Characterize plasmids in working directory `showcase` and import inc-types from `/inc-types` folder:
```bash
tadrep -v -o showcase characterize --inc-types /inc-types/inc-types.fasta
```

If inc-types is already present inside the working directory, the parameter `--inc-types` can be omitted:
```bash
tadrep -v -o showcase characterize
```

If you set up a database you can import it into the working directory `showcase` with the `--db` parameter:
```bash
tadrep -v -o showcase characterize --db database/plsdb/db.json --inc-types /inc-types/inc-types.fasta
```

---
## Cluster
---

This module aims to group plasmids with similar features together.
This is planned for a future release, currently the only option is to skip clustering, where each plasmid is put in its own individual group.

```bash
usage: TaDReP cluster [-h] [--skip]

optional arguments:
  -h, --help  show this help message and exit

Parameter:
  --skip, -s  Skips clustering, one group for each plasmid
```

### Example
```bash
tadrep -v -o showcase cluster --skip
```

---
## Detect
---
```bash
usage: TaDReP detect [-h] [--genome GENOME [GENOME ...]] [--min-contig-coverage [1-100]] [--min-contig-identity [1-100]] [--min-plasmid-coverage [1-100]] [--min-plasmid-identity [1-100]]
                     [--gap-sequence-length GAP_SEQUENCE_LENGTH]

optional arguments:
  -h, --help            show this help message and exit

Input / Output:
  --genome GENOME [GENOME ...], -g GENOME [GENOME ...]
                        Draft genome path

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
```

### Examples

---
## Visualize
---
```bash
usage: TaDReP visualize [-h] [--plotstyle {bigarrow,arrow,bigbox,box,bigrbox,rbox}] [--labelcolor LABELCOLOR] [--linewidth LINEWIDTH] [--arrow-shaft-ratio ARROW_SHAFT_RATIO] [--size-ratio SIZE_RATIO]
                        [--labelsize LABELSIZE] [--labelrotation LABELROTATION] [--labelhpos {left,center,right}] [--labelha {left,center,right}] [--interval-start [0-100]] [--number-of-intervals [1-100]]
                        [--omit_ratio [0-100]]

optional arguments:
  -h, --help            show this help message and exit

Style:
  --plotstyle {bigarrow,arrow,bigbox,box,bigrbox,rbox}
                        Contig representation in plot
  --labelcolor LABELCOLOR
                        Contig label color
  --linewidth LINEWIDTH
                        Contig edge linewidth
  --arrow-shaft-ratio ARROW_SHAFT_RATIO
                        Size ratio between arrow head and shaft
  --size-ratio SIZE_RATIO
                        Contig size ratio to track

Label:
  --labelsize LABELSIZE
                        Contig label size
  --labelrotation LABELROTATION
                        Contig label rotation
  --labelhpos {left,center,right}
                        Contig label horizontal position
  --labelha {left,center,right}
                        Contig label horizontal alignment

Gradient:
  --interval-start [0-100]
                        Percentage where gradient should stop
  --number-of-intervals [1-100]
                        Number of gradient intervals

Omit:
  --omit_ratio [0-100]  Omit contigs shorter than X percent of plasmid length from plot
```

### Examples

---
---
## Issues & Feature Requests

TaDReP is brand new and like in every software, expect some bugs lurking around. So, if you run into any issues with TaDReP, we'd be happy to hear about it.
Therefore, please, execute it in verbose mode (`-v`) and do not hesitate to file an issue including as much information as possible:

- a detailed description of the issue
- command line output
- log file (`tadrep.log`)
- a reproducible example of the issue with an input file that you can share _if possible_
