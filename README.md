# Alternative Splicing Evolution

## Summary
This repository contains the code for our paper:
<b>The role of alternative splicing in driving yet another phase transition in genomic complexity [doi] </b>

## Code
1. `step1_get_data.py`
  - iterates through different biological kingdoms and downloads the primary gene annotation file and corresponding protein sequence file
  - calculates descriptive statistics for gene length, exon count, and protein length for each species
2. `step2_merge_all_data.py`
  - aggregates result files
3. `step3_avg_variance_exon.py`
  - generates a plot between the average and standard deviation of exon count and the best linear regression curve
4. `step4_mean_plots_gene_length.py`
  - generates a log-log plot between the mean gene length and mean exon count
5. `step5_mean_plots_protein_length.py`
  - generates a plot between mean exon count and mean protein length

## Data
### Protein coding genes
Protein coding gene annotations in `GFF3` format are obtained from Ensembl and EnsemblGenomes FTP server following the official group categorization. Archaea and Bacteria are not included in this study for not exhibiting alternative splicing the way eukaryotes do in general. 
|Divisions|Release|
|--------|-------|
|[EnsemblVertebrates](https://ftp.ensembl.org/pub/current/))|114|
|[EnsemblMetazoa](http://ftp.ensemblgenomes.org/pub/metazoa/current/)|61|
|[EnsemblPlants](http://ftp.ensemblgenomes.org/pub/plants/current/)|61|
|[EnsemblFungi](http://ftp.ensemblgenomes.org/pub/fungi/current/)|61|
|[EnsemblProtists](http://ftp.ensemblgenomes.org/pub/protists/current/)|61|

### Protein sequences
Protein sequences in `FASTA` format are obtained from Ensembl and EnsemblGenomes FTP server
  - EnsemblVertebrates: `/pub/current/fasta/${species}/pep`
  - Other Divisions: `/pub/current/${kingdom}/fasta/${species}/pep`

## Results
After data acquisation and aggregation, `panel_a_data.tsv` is generated which serves as the primary dataset for statistical analysis and visualization.



