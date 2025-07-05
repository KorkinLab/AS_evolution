# Alternative Splicing Evolution

## Summary
This repository contains the code for our paper:
<b>The role of alternative splicing in driving yet another phase transition in genomic complexity [doi] </b>

## Code
1. `step1_get_data.py`
  - Download genome annotation and protein sequence for all species on Ensembl
  - Calculate descriptive statistics for gene length, exon count, and protein length for each species
2. `step2_merge_all_data.py`
  - Helper script to aggregate result files when step1 got interrupted.
3. `step3_avg_variance_exon.py`
  - Script for Fig.1A
4. `step4_mean_plots_gene_length.py`
  - Script for Fig.1B
5. `step5_mean_plots_protein_length.py`
  - Script for Fig.1B inset

## Data
### Genome Annotation
Protein-coding gene annotations in `GFF3` format are obtained from Ensembl and EnsemblGenomes FTP server. Archaea and Bacteria are not included in this study for the lack of alternative splicing in general. 

### Protein sequences
Protein sequences in `FASTA` format are obtained from Ensembl and EnsemblGenomes FTP server
|Database|Release|Genome Annotation|Protein Sequence|
|--------|-------|-------|-------|
|Ensembl|114|[Link](https://ftp.ensembl.org/pub/release-114/fasta)|[Link](https://ftp.ensembl.org/pub/release-114/fasta)|
|EnsemblMetazoa|61|[Link](http://ftp.ensemblgenomes.org/pub/metazoa/release-61/)|[Link](http://ftp.ensemblgenomes.org/pub/metazoa/release-61/)|
|EnsemblPlants|61|[Link](http://ftp.ensemblgenomes.org/pub/plants/release-61/)|[Link](http://ftp.ensemblgenomes.org/pub/plants/release-61/)|
|EnsemblFungi|61|[Link](http://ftp.ensemblgenomes.org/pub/fungi/release-61/)|[Link](http://ftp.ensemblgenomes.org/pub/fungi/release-61/)|
|EnsemblProtists|61|[Link](http://ftp.ensemblgenomes.org/pub/protists/release-61/)|[Link](http://ftp.ensemblgenomes.org/pub/protists/release-61/)|

## Results
After data acquisition and aggregation, `panel_a_data.tsv` is generated and serves as the primary dataset for statistical analysis and visualization.



