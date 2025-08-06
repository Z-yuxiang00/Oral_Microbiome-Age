import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.grid'] = False
plt.rcParams['axes.grid.axis'] = 'both'
plt.rcParams['axes.grid.which'] = 'major'
plt.rcParams['font.family'] = 'Times New Roman'
import os
# os.chdir(r"pathway")
data=pd.read_excel(r"Age bin traversal across four data types.xlsx",header=0)
data_long = data.melt(id_vars=['Datatype', 'Algorithm'], var_name='Age bin', value_name='MAE') # 将数据转化为长格式
plt.figure(figsize=(8, 6))
color_palette = {
    'XGBoost': '#7150B4',
    'Neural network': '#3C78BA',
    'RF': '#1D4F7B',
    'Naïve Bayes': '#359787',
    'KNN': '#E69F00',
    'SVM': '#F06292',
    'LDA': '#C63D18',
    'LR': 'red'
}
g = sns.FacetGrid(data_long, col="Datatype", hue="Algorithm", height=5, aspect=1.2, palette=color_palette, col_wrap=2)# 创建 FacetGrid，按数据类型划分子图
g.map(sns.lineplot, 'Age bin', 'MAE', marker="o", markersize=2.5, linewidth=1.5) # 绘制每个子图
g.set_axis_labels('Age bin (years)', 'MAE (years)')
g.add_legend(loc='upper left', bbox_to_anchor=(0.9, 0.95),title='', frameon=False) #frameon=False图例不显示边框
# 保存图形为TIFF格式
plt.savefig(r"Fig. S2.tiff", dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()