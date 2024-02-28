# Genome assembly summary
A snakemake workflow for generating and combining genome assembly stats, quality estimates, and taxonomy info. Inputs are fasta files of genome assemblies (currently only non-compressed, let me know if adding accepting gzipped as input would be helpful).

## Overview

It uses:

  - [bit](https://github.com/AstrobioMike/bit#bioinformatics-tools-bit) for generating assembly summary stats
  - [checkm2](https://github.com/chklovski/CheckM2#checkm2) for estimating quality of bacteria/archaea
  - [GTDB-tk](https://github.com/Ecogenomics/GTDBTk#gtdb-tk) for assigning taxonomy of bacteria/archaea
  - [eukcc](https://github.com/Finn-Lab/EukCC#eukcc) for estimating quality of eukarya
  - [CAT](https://github.com/dutilh/CAT#cat-and-bat) with the NCBI nr database for assigning taxonomy of eukarya

Before running it, you first need to set some variables in the config.yaml file (there are notes in there). It cannot currently run on a mix of input bacteria/archaea genomes and eukaryotic genomes, it can only run on bacteria/archaea by themselves, or eukarya by themselves (as set by a parameter in the config.yaml file).

In the config.yaml file, you mostly just need to point to where the input fasta files are, specify what their extensions are, and run the snakemake workflow as exemplified below. It will summarize the assemblies, estimate quality, and assign taxonomy via the programs listed above, and ultimately produce an output table like this:

```bash
Assembly        Total contigs  Total length  Ambiguous characters  GC content  Maximum contig length  Minimum contig length  N50        L50  Est. Completeness (%)  Est. Redundancy (%)  Domain    Phylum          Class                Order             Family             Genus           Species
input_genome_A  1              5,276,633     0                     61.17       5,276,633              5,276,633              5,276,633  1    99.99                  0.96                 Bacteria  Proteobacteria  Gammaproteobacteria  Pseudomonadales   Pseudomonadaceae   Pseudomonas_E   Pseudomonas_E fulva
input_MAG_B     4              2,702,105     0                     33.15       2,601,030              30,881                 2,601,030  1    90.22                  2.78                 Bacteria  Firmicutes      Bacilli              Staphylococcales  Staphylococcaceae  Staphylococcus  Staphylococcus saprophyticus
```
All required databases will be setup by the workflow if they don't exist already whenever they are used for the first time.

The workflow can be retrieved programmatically with [bit](https://github.com/AstrobioMike/bit) using `bit-get-genome-summarize-wf`.

If not using `bit`, but just pulling from here, just be sure to install snakemake, e.g.:

```bash
conda install -c conda-forge -c bioconda -c defaults snakemake
```

## Example running
After variables are set in the config.yaml, here's an example of how it could be run:

```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 4 -p
```

- `--use-conda` – this specifies to use the conda environments included in the workflow
- `--conda-prefix` – this allows us to point to where the needed conda environments should be stored. Including this means if we use the workflow on a different dataset somewhere else in the future, it will re-use the same conda environments rather than make new ones. The value listed here, `${CONDA_PREFIX}/envs`, is the default location for conda environments (the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
- `-j` – this lets us set how many jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
- `-p` – specifies to print out each command being run to the screen

See `snakemake -h` for more options and details.

## Version info
A version for the workflow as a whole is denoted at the top of the Snakefile. All versions of programs used can be found in their corresponding conda yaml file in the envs/ directory. 
