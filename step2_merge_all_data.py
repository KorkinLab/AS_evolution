import pandas as pd
import glob

data_all = pd.DataFrame(columns=['database', 'mean_gene_length', 'variance_gene_length', 'mean_exon_count', 'variance_exon_count',  'mean_protein_length', 'variance_protein_length'])
tsv_files = glob.glob('panel_a_data*.tsv')
for tsv_file in tsv_files:
    data = pd.read_csv(tsv_file, sep='\t', index_col='gtf_name')
    data_all = pd.concat([data_all, data])

# Check for duplicate indices
if data_all.index.duplicated().any():
    print("Warning: Duplicate indices found!")
else:
    print("No duplicate indices.")

data_all.to_csv('panel_a_data.tsv', sep='\t', index_label='gtf_name')