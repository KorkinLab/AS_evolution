import pandas as pd
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np

mean_df = pd.read_csv('panel_a_data.tsv', sep='\t', index_col='gtf_name')

outlier_species = ['symbiodinium_microadriaticum_gca_001939145']
mean_df = mean_df[~mean_df.index.isin(outlier_species)]

feature = 'mean_exon_count'
x = []
y = []
y_label = feature.replace('_',' ').title()

for species, row in mean_df.iterrows():
    mean_exon_count = row.loc['mean_exon_count']
    if mean_exon_count == 1:
        continue
    gene_length = row.loc['mean_protein_length']
    feature_count = row.loc[feature]

    database = row['database']
    if database == 'metazoa':
        color='#b77df0'
    elif database == 'plants':
        color='#34a853'
    elif database == 'vertebrate':
        color='#ea4335'
    elif database == 'fungi':
        color='#3c83f7'
    elif database == 'protists':
        color='#f4c20d'
    else:
        continue

    plt.scatter(feature_count, gene_length, color=color, marker='.', s=8)
    if species == 'homo_sapiens':
        human_gene_length = gene_length
        human_feature_count = feature_count
        

    x.append(gene_length)
    y.append(feature_count)
plt.scatter(human_feature_count, human_gene_length, color='black', marker='o', facecolors='none', s=80, linewidths=2)




# Set labels
plt.ylabel('Mean protein length (AA)', fontsize=16)
plt.xlabel('Mean exon count',fontsize=16)

# # Add legend for color codes (only once per plot)
# legend_elements = [
#     Line2D([0], [0], marker='.', color='w', label='Fungi', markerfacecolor='#3c83f7', markersize=6, linestyle='None'),
#     Line2D([0], [0], marker='.', color='w', label='Protists', markerfacecolor='#f4c20d', markersize=6, linestyle='None'),
#     Line2D([0], [0], marker='.', color='w', label='Plants', markerfacecolor='#34a853', markersize=6, linestyle='None'),
#     Line2D([0], [0], marker='.', color='w', label='Invertebrate', markerfacecolor='#b77df0', markersize=6, linestyle='None'),
#     Line2D([0], [0], marker='.', color='w', label='Vertebrate', markerfacecolor='#ea4335', markersize=6, linestyle='None'),
#     Line2D([0], [0], marker='o', color='black', label='Homo sapiens', markerfacecolor='none', markersize=5, linestyle='None')
# ]
# plt.legend(
#     handles=legend_elements,
#     loc='upper left',
#     borderaxespad=0.1,
#     fontsize=8,
#     ncol=2
# )

# # Add a point to show where protein length stops growing
# plt.scatter(1.5, 0, color='red', marker='|', facecolors='red', s=80)
# plt.text(5, 0.4, '1.5', va='center', ha='left', fontsize=8, color = 'red')
# plt.annotate(
#     '', 
#     xy=(1.5, 0.1), 
#     xytext=(5, 0.4), 
#     arrowprops=dict(arrowstyle='->', color='black', lw=1)
# )
# plt.ylim(bottom=0)
plt.tight_layout()
plt.axhline(y=500, color='#06c0ea', linestyle=(0, (6, 2)), linewidth=1.5)
# plt.yscale('log')
# plt.xscale('log')
plt.yticks([200, 300, 400, 500, 600, 700, 800], ['200', '300', '400', '500', '600', '700', '800'], fontsize=16)
plt.xticks(fontsize=16)
plt.savefig(f'protein_length_exon_plot_only_AS.png', dpi=600)

plt.close()