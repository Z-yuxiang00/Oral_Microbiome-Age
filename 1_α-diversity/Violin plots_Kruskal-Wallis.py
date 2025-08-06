import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scikit_posthocs import posthoc_dunn

os.chdir(r"pathway")
data = pd.read_csv(r"combined_data.csv", header=0)

plt.style.use('ggplot')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['font.family'] = 'Times New Roman'
my_pal_gender = {"Young": '#F9766E', "Middle-aged": '#01BA38', "Old": '#619DFF'}
fig, ax = plt.subplots(figsize=(3.5, 4.5))

sns.violinplot(ax=ax, x=data["Age_group"],
               y=data["Simpson"],
               bw_method=0.3, linewidth=1, width=0.8, alpha=1.0,
               palette=my_pal_gender, hue=None, density_norm='area', inner=None, saturation=0.9)

sns.boxplot(
    ax=ax,
    x=data["Age_group"],
    y=data["Simpson"],
    data=data,
    color="black",
    width=0.10,
    zorder=10,
    showcaps=True,
    boxprops={'facecolor': 'white', "zorder": 10, 'linewidth': 0.8, 'zorder': 10},
    medianprops={'color': 'black', 'linewidth': 1},
    showfliers=True,
    whiskerprops={'linewidth': 0.5, "zorder": 10},
    saturation=0.3,
    orient="v"
)

kw_stat, kw_p_val = stats.kruskal(
    data[data["Age_group"] == "Young"]["Simpson"],
    data[data["Age_group"] == "Middle-aged"]["Simpson"],
    data[data["Age_group"] == "Old"]["Simpson"]
)
print(f"Kruskal-Wallis H-statistic: {kw_stat}, p-value: {kw_p_val}")

dunn_results = posthoc_dunn(
    data,
    val_col="Simpson",
    group_col="Age_group",
    p_adjust="bonferroni"
)
print(dunn_results)

def add_significance(ax, x1, x2, y, significance, y_offset_multiplier=1.0):
    y_offset = y * 0.05 * y_offset_multiplier
    line_color = 'black'
    line_width = 1.0
    ax.plot([x1, x1, x2, x2], [y, y + y_offset, y + y_offset, y], color=line_color, linewidth=line_width)
    ax.text((x1 + x2) / 2, y + y_offset, significance, ha='center', va='bottom', fontsize=12, color='black')

groups_list = list(my_pal_gender.keys())
y_base = data["Simpson"].max() + 0.01

for i, group1 in enumerate(groups_list):
    for j, group2 in enumerate(groups_list):
        if i < j and dunn_results.loc[group1, group2] < 0.05:
            significance = "***" if dunn_results.loc[group1, group2] < 0.001 else "**" if dunn_results.loc[group1, group2] < 0.01 else "*"
            add_significance(ax, i, j, y_base, significance)
            y_base += 0.1

ax.tick_params(axis='x', labelsize=6, labelrotation=0, direction='out', length=4, width=1, colors='black')
ax.tick_params(axis='y', labelsize=6, labelrotation=0, direction='out', length=4, width=1, colors='black')
ax.spines['left'].set_edgecolor('black')
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_edgecolor('black')
ax.spines['bottom'].set_linewidth(1.2)
ax.set_xlabel('')
ax.set_ylabel('Simpson', fontweight='bold', color='black', fontsize=6)

p_value = kw_p_val
ax.text(1, 1.03, f'Kruskal-Wallis p-value: {p_value:.4f}', ha='center', va='bottom', fontsize=12,  color='black')

plt.savefig('Violin_plot_Simpson.jpg', format='jpg', dpi=300)
plt.show()
