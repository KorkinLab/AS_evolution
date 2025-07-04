import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mean_df = pd.read_csv('panel_a_data.tsv', sep='\t', index_col='gtf_name')

outlier_species = ['symbiodinium_microadriaticum_gca_001939145']
mean_df = mean_df[~mean_df.index.isin(outlier_species)]

feature ='exon_count'

# Plot variance of feature against average feature
x = []
y = []

for species, row in mean_df.iterrows():
    
    avg_feature = row.loc['mean_' + feature]
    if avg_feature == 1:
        continue
    variance_feature = np.sqrt(row.loc['variance_' + feature])

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
        continue  # Skip if database is not recognized

    plt.scatter(avg_feature, variance_feature, color=color, marker='.', s=5)
    if species == 'homo_sapiens':
        human_avg_feature = avg_feature
        human_variance_feature = variance_feature
        
    x.append(avg_feature)
    y.append(variance_feature)

plt.scatter(human_avg_feature, human_variance_feature, color='black', marker='o', facecolors='none', s=50, linewidths=1.5)


# Fit 1
x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
mask = ~np.isnan(x) & ~np.isnan(y)
slope, intercept, r_value, _, _ = linregress(x[mask], y[mask])
plt.plot(x[mask], slope * x[mask] + intercept, color='black', linestyle='-',linewidth=0.5)
plt.text(0.94, 0.765, f'$r^2$ = {r_value**2:.3f}', transform=plt.gca().transAxes, fontsize=12, va='top', ha='right')
plt.text(0.80, 0.85, f'$\sigma_E$ = {slope:.3f}$\mu_E$ - {abs(intercept):.3f}', transform=plt.gca().transAxes, fontsize=12, va='top', ha='right')


plt.xlabel(f"Mean exon count, $\mu_E$", fontsize=12)
plt.ylabel(f"Standard deviation of exon count, $\sigma_E$", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Add legend for color codes (only once per plot)
legend_elements = [
    Line2D([0], [0], marker='.', color='w', label='Fungi', markerfacecolor='#3c83f7', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='.', color='w', label='Protists', markerfacecolor='#f4c20d', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='.', color='w', label='Plants', markerfacecolor='#34a853', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='.', color='w', label='Invertebrates', markerfacecolor='#b77df0', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='.', color='w', label='Vertebrates', markerfacecolor='#ea4335', markersize=10, linestyle='None'),
    Line2D([0], [0], marker='o', color='black', label='Homo sapiens', markerfacecolor='none', markersize=8,linestyle='None')
]
plt.legend(handles=legend_elements, loc='lower right', borderaxespad=0.1, fontsize=12)

plt.tight_layout()
plt.savefig(f'avg_std_{feature}_plot_only_AS.png', dpi=600)
plt.close()