import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
os.chdir(r"pathway")
metadata = pd.read_csv('sample_information.csv')
abundance = pd.read_csv('Phylum_aggregated_data_relative.csv')
abundance = abundance.set_index('taxonomy').T.reset_index().rename(columns={'index': 'name'})
merged_data = pd.merge(metadata, abundance, on='name')
long_data = pd.melt(
    merged_data,
    id_vars=['name', 'Age_group'],
    var_name='Phylum',
    value_name='Abundance',
    value_vars=abundance.columns[1:]
)
pivot_data = long_data.pivot(index='name', columns='Phylum', values='Abundance').fillna(0)
mean_abundance = pivot_data.mean().sort_values(ascending=False)
top_phyla = mean_abundance.head(10).index.tolist()
pivot_data['Other'] = pivot_data.drop(columns=top_phyla).sum(axis=1)
pivot_data = pivot_data[top_phyla + ['Other']]
desired_order = ["Young", "Middle-aged", "Old"]
metadata['Age_group'] = pd.Categorical(metadata['Age_group'], categories=desired_order, ordered=True)
sample_order = metadata[['name', 'Age_group']].set_index('name').sort_values(by='Age_group')
pivot_data = pivot_data.loc[sample_order.index]
grouped_samples = (
    sample_order
    .join(pivot_data[['d__Bacteria; p__Proteobacteria']])
    .sort_values(by=['Age_group', 'd__Bacteria; p__Proteobacteria'])
)
sorted_indices = grouped_samples.index
pivot_data = pivot_data.loc[sorted_indices]
sample_order = sample_order.loc[sorted_indices]
plt.style.use('ggplot')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.grid'] = False
plt.rcParams['font.family'] = 'Times New Roman'
fig, ax = plt.subplots(figsize=(12, 9.6))
bottom = np.zeros(len(pivot_data))
soft_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
               '#ffd92f', '#e5c494', '#b3b3b3', '#a6cee3', '#b2df8a', '#fb9a99']
for i, phylum in enumerate(pivot_data.columns):
    ax.bar(pivot_data.index, pivot_data[phylum], bottom=bottom, label=phylum, color=soft_colors[i])
    bottom += pivot_data[phylum]
plt.xticks([])
age_groups = sample_order['Age_group'].value_counts().loc[desired_order].cumsum()
desired_order = ["Young", "Middle-aged", "Old"]
age_labels = list(sample_order['Age_group'].unique())
colors = ['#EFBF47', '#3199D3', '#4BAF71']
for i, (age, pos) in enumerate(zip(age_groups.index, age_groups)):
    start_geo = 0 if i == 0 else age_groups.iloc[i - 1]
    end_geo = pos
    label_width = end_geo - start_geo
    ax.add_patch(
        plt.Rectangle((start_geo - 0.5, -0.07), label_width, 0.07, facecolor=colors[i], edgecolor='black', lw=1))
    ax.text((start_geo + end_geo) / 2, -0.035, age_labels[i], ha='center', va='center', fontsize=10, weight='bold',
            color='black')
    ax.add_patch(plt.Rectangle((start_geo - 0.5, 0), label_width, 1, fill=False, edgecolor='black', lw=1.5))
ax.set_ylim(-0.07, 1)
ax.set_xlim(-0.5, len(pivot_data) - 0.5)
legend = plt.legend(
    bbox_to_anchor=(0.5, -0.02),
    loc='upper center',
    ncol=4,
    title="Phylum",
    fontsize=10.22,
    title_fontsize=10.22,
    borderpad=0.5,
    frameon=True,
    framealpha=0.5,
    edgecolor='black'
)
legend.get_frame().set_linewidth(1)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.subplots_adjust(bottom=0.25)


plt.savefig('Phylum.jpg', format='jpg', dpi=300)
plt.show()
merged_output = pd.concat([sample_order, pivot_data], axis=1)
merged_output.to_csv('sample_phylum_abundance.csv', index=True)
