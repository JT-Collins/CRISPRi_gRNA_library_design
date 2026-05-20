# CRISPRi Spacer Design and Oligo Construction Pipeline


- [Overview](#overview)
- [Environment Setup](#environment-setup)
  - [Install WSL](#install-wsl)
  - [Install micromamba](#install-micromamba)
  - [Create environment](#create-environment)
  - [Install ViennaRNA](#install-viennarna)
- [PART I — Spacer Design Engine](#part-i--spacer-design-engine)
  - [Files](#files)
  - [What the Spacer Designer Does](#what-the-spacer-designer-does)
    - [1. PAM Detection](#1-pam-detection)
    - [2. Guide Extraction](#2-guide-extraction)
    - [3. GC Filtering](#3-gc-filtering)
    - [4. ORF Positional Filtering](#4-orf-positional-filtering)
    - [5. Off-target Filtering](#5-off-target-filtering)
    - [6. Guide Selection](#6-guide-selection)
  - [Run](#run)
  - [Config Example](#config-example)
  - [Output Contains](#output-contains)
- [PART II — Oligo Construction](#part-ii--oligo-construction)
  - [File](#file)
  - [What the Oligo Script Does](#what-the-oligo-script-does)
  - [Singleplex Design Logic](#singleplex-design-logic)
  - [Multiplex Array Construction](#multiplex-array-construction)
  - [Junction Sequences](#junction-sequences)
  - [Spacer Order Optimization](#spacer-order-optimization)
  - [RNAfold Scoring](#rnafold-scoring)
- [running the script](#running-the-script)
  - [Run Singleplex](#run-singleplex)
  - [Run Multiplex](#run-multiplex)
  - [Run Optimized Multiplex](#run-optimized-multiplex)
- [Example Input](#example-input)
- [Design Principles](#design-principles)
  - [Spacer Design](#spacer-design)
  - [Multiplex Design](#multiplex-design)
- [Troubleshooting](#troubleshooting)

# Overview

This pipeline consists of two coordinated Python systems:

1.  **Spacer design engine** (`design_library.py`, `crispri_design.py`)
2.  **Oligo construction engine** (`oligos.py`)

It enables: - Genome-scale CRISPRi guide design (Cas9 / Cas12a) -
Off-target filtered spacer selection - Construction of cloning-ready
oligos - Multiplex array optimization using RNA secondary structure

Most of the functionality works on Windows via PowerShell as long as you
have Python 3.x installed. If you want to use the CRISPRi multiplex
functionality that uses viennaRNA
(<https://github.com/ViennaRNA/ViennaRNA>) to determine RNA hairpin
structure (which you probably should) then you will need to run the
script via Ubuntu on WSL. I expect it works on mac too…

# Environment Setup

## Install WSL

The Windows Subsystem for Linux (WSL) allows you to install and run a
Linux distribution directly within Windows. The default distribution is
Ubuntu.

You can install everything needed for WSL with a single command. Open
PowerShell in administrator mode (right-click → Run as administrator),
run the command below, then restart your machine. You can also install
WSL via the Microsoft Store.

``` bash
wsl --install
```

## Install micromamba

`micromamba` is a lightweight version of the `mamba` package manager. It
allows you to create isolated environments that contain their own
programs and dependencies. For example, you could have:

- one environment for a genome annotation pipeline
- another for RNA-seq analysis

Why not just install everything globally? Because many tools depend on
specific (and often incompatible) versions of software. Managing
multiple versions manually is difficult and error-prone. Using
environments solves this problem by allowing you to run, for example,
Python 2.x in one environment and Python 3.8 in another—without
conflicts.

Read about it here <https://mamba.readthedocs.io/en/latest/> and if the
installation below doesn’t work for you there are other ways documented
in the link.

``` bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
mkdir -p ~/micromamba
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
```

## Create environment

``` bash
micromamba create -n crispr python=3.10 pandas
micromamba activate crispr
```

## Install ViennaRNA

``` bash
micromamba install -c bioconda viennarna
RNAfold --version
```

------------------------------------------------------------------------

# PART I — Spacer Design Engine

## Files

- `design_library.py` → this is the script you call on the cmd line
- `crispri_design.py` → just needs to be in the same folder
- `config.txt` → adjust the parameters here to your likeing. The
  defaults can probably remain the same, just change out
  `target_fasta` - a multi-fasta file with the CDS you want targets for.
  it could be a single gene or every gene in the genome and `genome`
  this is the full genome of your bug and is used to check for
  off-target sites.
  
  These scripts were influenced by Ban H., *et al*. Cross-Strain Transferability  of CRISPRi Systems and Design Rules from Laboratory to Clinical *Escherichia coli* Strains. ACS Synth Biol. 2026;15(5):1993-2010.

------------------------------------------------------------------------

## What the Spacer Designer Does

### 1. PAM Detection

Each CRISPR system defines a PAM motif:

- Cas9: `NGG`
- Cas12a: `TTTN`

The script scans every input sequence using regex and identifies valid
PAM positions.

------------------------------------------------------------------------

### 2. Guide Extraction

For each PAM: - Extract adjacent protospacer - Respect orientation (Cas9
vs Cas12a differ) - Generate candidate guides

------------------------------------------------------------------------

### 3. GC Filtering

Each guide is filtered:

``` text
GC% = (G + C) / length
```

Only guides within acceptable GC bounds are retained. This can be
adjusted in `config.txt`.

------------------------------------------------------------------------

### 4. ORF Positional Filtering

CRISPRi is most effective near transcription start sites.

``` text
position <= orf_length * orf_cutoff
```

Example: - `orf_cutoff = 0.6` - Only first 60% of gene retained

------------------------------------------------------------------------

### 5. Off-target Filtering

The genome is indexed using k-mers:

- Fast lookup of potential matches
- Seed region evaluated first

A guide is rejected if:

- ≤2 mismatches overall AND
- ≤1 mismatch in seed region

This approximates genome-wide specificity efficiently. Essentially it
checks the whole genome for things with high homology to the spacer and
if they occur (other than in your gene) that spacer is thrown out. This
limits off targets.

------------------------------------------------------------------------

### 6. Guide Selection

Final guides are selected based on:

- Passing all filters
- Ranking by position and quality
- Limiting to `grna_number` per gene

------------------------------------------------------------------------

## Run

Place all of the files in the same folder and navigate there on the
Ubuntu cmd line.

To run make any adjustments to the config file (this is a simple text
file) and run:

``` bash
python design_library.py config.txt
```

------------------------------------------------------------------------

## Config Example

``` ini
[design]
crispri_system = fn-dcas12a
target_fasta = targets.fasta
grna_number = 6
orf_cutoff = 0.6
min_gc_content = 30
max_gc_content = 65
genome = genome.fasta
prefix = experiment1
```

------------------------------------------------------------------------

## Output Contains

- gene name
- guide sequence
- GC content
- strand
- position

------------------------------------------------------------------------

# PART II — Oligo Construction

Now that you have a list of CRISPR targets (spacers) you will want to
order oligos to put into your plasmid. This is specific to Cas12 - I
have not added the ability to use Cas9 as yet.

## File

- `oligos.py`

------------------------------------------------------------------------

## What the Oligo Script Does

Transforms spacers into:

1.  **Single cloning oligos**
2.  **Multiplex CRISPR arrays**
3.  **Optimized spacer orders using RNAfold**

------------------------------------------------------------------------

## Singleplex Design Logic

This is pretty trivial and doesn’t really need a script, but…

Each spacer becomes:

Top strand:

``` text
5 prime to 3 prime
tagat + spacer + a
```

Bottom strand:

``` text
5 prime to 3 prime
ttaaa + reverse_complement(spacer) + a
```

------------------------------------------------------------------------

## Multiplex Array Construction

This is based upon Liao C, *et al*. Modular one-pot assembly of CRISPR
arrays enables library generation and reveals factors influencing crRNA
biogenesis. Nat Commun. 2019;10(1):2948.PubMed PMID: 31270316.

It has been modified to fit our setup where the first and last repeat is
already on the plasmid.

Arrays follow:

``` text
overhang - spacer1 - repeat - spacer2 - repeat - spacer3 - overhang
```

Additional engineered junctions are inserted to:

- minimize recombination
- improve cloning

------------------------------------------------------------------------

## Junction Sequences

As outlined in the paper it is possible to add unique junction sequences
at the 3’ end of the spacer that enable simple one pot assembly of
CRISPR arrays

``` text
GCTC, GAGT, AACG, TACA, TTCT, AGAA, TGGC, CCCT
```

------------------------------------------------------------------------

## Spacer Order Optimization

Multiple spacer sequences in CRISPR arrays can form hairpin loops
causing some or all of the targets to fail to be repressed. We try to
mitigate this by shuffling the spacers in all possible permutations and
using ViennaRNA to examine their structure.

All spacer permutations are generated using:

``` python
itertools.permutations(spacers)
```

You don’t need to run this it is done automatically upon flagging
`--optimize`

Each is then evaluated using ViennaRNA.

------------------------------------------------------------------------

## RNAfold Scoring

``` bash
RNAfold
```

Metric:

- Minimum Free Energy (MFE)

Interpretation:

- Less negative → less structured → better processing

------------------------------------------------------------------------

# running the script

We can do this in two ways:

## Run Singleplex

``` bash
python oligos.py <input .fasta> --output <output.csv> --mode single
```

## Run Multiplex

``` bash
python oligos.py <input .fasta> --mode multiplex --output <output.csv>
```

## Run Optimized Multiplex

``` bash
python oligos.py <input .fasta> --mode multiplex --optimize --output <output.csv>
```

------------------------------------------------------------------------

# Example Input

Input should be a fasta file

``` fasta
>geneA
TTTATGGGACAAAACCCTACA
>geneB
GCTCGTGAAGAAATGAGTACA
>geneC
TTTATGGGACAAAACCCTACA
```

------------------------------------------------------------------------

# Design Principles

## Spacer Design

- Favor 5′ gene region
- GC 30–65%
- Avoid repetitive motifs

## Multiplex Design

- Keep arrays ≤6 spacers
- Always optimize order
- Avoid stable RNA structures

------------------------------------------------------------------------

# Troubleshooting

| Issue                    | Cause             | Fix                       |
|--------------------------|-------------------|---------------------------|
| No guides                | strict filters    | relax GC or ORF cutoff    |
| Many off-target failures | genome similarity | adjust mismatch tolerance |
| RNAfold missing          | install issue     | reinstall ViennaRNA       |

------------------------------------------------------------------------
