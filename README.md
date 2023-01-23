# Genome assembly summary
Snakemake workflow for generating and combining genome assembly stats and taxonomy info. Inputs are fasta files of genome assemblies (currently only non-compressed).

It uses:

  - [bit](https://github.com/AstrobioMike/bit#bioinformatics-tools-bit) for generating assembly summary stats
  - [checkm2](https://github.com/chklovski/CheckM2#checkm2) for estimating quality of bacteria/archaea
  - [GTDB-tk](https://github.com/Ecogenomics/GTDBTk#gtdb-tk) for assigning taxonomy of bacteria/archaea
  - [eukcc](https://github.com/Finn-Lab/EukCC#eukcc) for estimating quality of eukarya
  - [CAT](https://github.com/dutilh/CAT#cat-and-bat) with NCBI database for assigning taxonomy of eukarya

Set needed in config.yaml file, all required databases should be setup by the workflow if they don't exist already whenever they are first used.

The workflow can be retrieved programmatically with [bit](https://github.com/AstrobioMike/bit) using `bit-get-genome-summarize-wf`.
