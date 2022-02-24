# genome-assembly-summary
Snakemake workflow for generating and combine genome assembly stats and taxonomy info. Internally used with GeneLab.

It uses:

  - [bit](https://github.com/AstrobioMike/bit#bioinformatics-tools-bit) for generating assembly summary stats
  - [checkm](https://github.com/Ecogenomics/CheckM#checkm) for estimating quality of bacteria/archaea
  - [eukcc](https://github.com/Finn-Lab/EukCC#eukcc) for estimating quality of eukarya
  - [GTDBtk](https://github.com/Ecogenomics/GTDBTk#gtdb-tk) for assigning taxonomy of bacteria/archaea
  - [CAT](https://github.com/dutilh/CAT#cat-and-bat) with NCBI database for assigning taxonomy of eukarya


As currently writtent, it does not automatically handle databases set up, so they need to be set up and linked to their appropriate conda environments manually when first being used.
