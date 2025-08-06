import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patches as mpatches
import os

# 设置样式
plt.style.use('ggplot')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.grid'] = False

# 设置工作目录和读入数据
# os.chdir(r"pathway")
df = pd.read_excel("Selected taxa from four data processing methods via LASSO.xlsx", sheet_name="Sheet1")

# 提取 Phylum 和 Genus 名称
df['Phyla'] = df['Phyla'].str.extract(r'p__(\w+)')
df['Genera'] = df['Genera'].str.extract(r'g__(\w+)')
df['Genera'] = df['Genera'].fillna('Unknown')

# --- 处理 Phyla ---
Phyla_counts = df.groupby(['Methods', 'Phyla']).size().unstack(fill_value=0)

Phyla_order = [p for p in ['Firmicutes', 'Proteobacteria'] if p in Phyla_counts.columns] + \
              [p for p in Phyla_counts.sum().sort_values(ascending=False).index if p not in ['Firmicutes', 'Proteobacteria']]
Phyla_counts = Phyla_counts[Phyla_order]

# --- 处理 Genera ---
required_genera = ['Neisseria', 'Treponema', 'Centipeda']

#5个属
genus_counts_total = df[~df['Genera'].isin(required_genera)]['Genera'].value_counts()
top5_genera = genus_counts_total.head(5).index.tolist()

# 最终展示的 Genera
displayed_genera = required_genera + top5_genera

# 分组计数
Genera_counts = df.groupby(['Methods', 'Genera']).size().unstack(fill_value=0)

# 保证列完整
for genus in displayed_genera:
    if genus not in Genera_counts.columns:
        Genera_counts[genus] = 0

# Others 列
others_genera = [g for g in Genera_counts.columns if g not in displayed_genera]
Genera_counts['Others'] = Genera_counts[others_genera].sum(axis=1)

# 仅保留需要的列
Genera_counts = Genera_counts[displayed_genera + ['Others']]

# --- 方法顺序固定 ---
methods_order = ['binarized data', 'CLR-transformed data', 'log2-transformed data', 'relative abundance data']
Phyla_counts = Phyla_counts.reindex(methods_order)
Genera_counts = Genera_counts.reindex(methods_order)

# 修改 index 区分 Phyla 和 Genera
Phyla_counts.index = [idx + '-P' for idx in Phyla_counts.index]
Genera_counts.index = [idx + '-G' for idx in Genera_counts.index]

# 合并
combined = pd.concat([Phyla_counts, Genera_counts], axis=0).fillna(0)
final_columns = Phyla_order + displayed_genera + ['Others']
combined = combined[final_columns]

# --- 配色 ---
Phyla_colors = ['#4ca69a', '#c1c15e', '#9a7fc9', '#e6735f', '#6d9eda',
                '#d79e47', '#87b840', '#ce7fb4', '#999999'][:len(Phyla_order)]

Genera_colors = ['#7dd0b3', '#f1a976', '#a19fce', '#ec92c7', '#a7c97b',
                 '#f5d849', '#b89ed0', '#8bb3e1'][:len(displayed_genera)]
colors = Phyla_colors + Genera_colors + ['#cccccc']  # Others

# --- 绘图 ---
fig, ax = plt.subplots(figsize=(8, 6))
combined.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.6)

x_labels = [label[:-2] if label.endswith(('-P', '-G')) else label for label in combined.index]
ax.set_xticks(range(len(combined)))
ax.set_xticklabels(x_labels, rotation=0, ha='center', color='black', fontsize=4.9)
ax.xaxis.set_tick_params(pad=-33)
ax.tick_params(axis='x', length=0)

# 虚线
ax.axhline(0, color='black', linewidth=1)
ax.set_xlim(-0.5, len(combined) - 0.5)
ax.set_ylim(-20, 130)

# 灰边框
for bar in ax.patches:
    bar.set_edgecolor('#cccccc')
    bar.set_linewidth(0.5)

# 标签框
Phyla_rect = patches.FancyBboxPatch((-0.5, -9), 4.5, 4,
                                    boxstyle="round,pad=0.4", facecolor="#e6f2ff", edgecolor='none')
Genera_rect = patches.FancyBboxPatch((3.9, -9), 4.5, 4,
                                     boxstyle="round,pad=0.4", facecolor="#ffe6e6", edgecolor='none')
ax.add_patch(Phyla_rect)
ax.add_patch(Genera_rect)

ax.text(1.5, -5.5, "Phyla", ha='center', va='top', fontsize=7, weight='bold')
ax.text(5.5, -5.5, "Genera", ha='center', va='top', fontsize=7, weight='bold')
ax.set_ylabel("ASV Count", color='black')
# 图例
handles, labels = ax.get_legend_handles_labels()
handle_dict = dict(zip(labels, handles))

Phyla_handles = [handle_dict[l] for l in Phyla_order if l in handle_dict]
Genera_handles = [handle_dict[g] for g in displayed_genera if g in handle_dict]

# Others 图例
if 'Others' in combined.columns:
    Genera_handles.append(mpatches.Patch(color='#cccccc', label='Others'))

legend1 = ax.legend(Phyla_handles, Phyla_order, title="Phyla",
                    bbox_to_anchor=(1, 0.95), loc='upper left', fontsize=7, title_fontsize=8)
legend2 = ax.legend(Genera_handles, displayed_genera + ['Others'], title="Genera",
                    bbox_to_anchor=(1, 0.5), loc='upper left', fontsize=7, title_fontsize=8)

for legend in [legend1, legend2]:
    legend.get_frame().set_linewidth(0)
    legend.get_frame().set_edgecolor('none')
ax.add_artist(legend1)

# 美化边框
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(False)
ax.plot([-0.5, -0.5], [0, 120], color='black', linewidth=1)

ax.set_yticks(range(0, 121, 20),)
ax.set_yticklabels([str(y) for y in range(0, 121, 20)], color='black')
ax.tick_params(axis='y', colors='black')

plt.subplots_adjust(bottom=0.28, right=0.75)
plt.savefig("Fig. 3b.tiff", dpi=300)
plt.show()
