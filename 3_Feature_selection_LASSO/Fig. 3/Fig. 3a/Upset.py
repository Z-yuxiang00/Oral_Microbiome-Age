import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
from matplotlib import cm

# 设置字体
plt.rcParams["font.family"] = "Times New Roman"

# 1. 读取 Excel 文件
file_path = "Selected ASVs from four data processing methods via LASSO.xlsx"
df_raw = pd.read_excel(file_path)

# 2. 将每列中的 ASV 收集为集合
sets = {col: set(df_raw[col].dropna()) for col in df_raw.columns}

# 3. 所有唯一 ASV
all_asvs = set().union(*sets.values())

# 4. 构建布尔矩阵
bool_df = pd.DataFrame(index=list(all_asvs))
for method, asv_set in sets.items():
    bool_df[method] = bool_df.index.isin(asv_set)

# 5. 构建 UpSet 数据
upset_data = from_indicators(bool_df.columns.tolist(), bool_df)

# 6. 创建颜色映射（每个交集一种颜色）
n_subsets = 2 ** len(bool_df.columns) - 1
colors = cm.get_cmap('tab10', n_subsets)  # tab10 是明亮对比色系

# 7. 绘图
plt.figure(figsize=(8, 6))
upset = UpSet(
    upset_data,
    show_counts=True,
    subset_size='count',
    facecolor='tab:blue'  # 应用颜色
)
upset.plot()

# 8. 保存图像
plt.savefig("Fig. 3a.tiff", dpi=300, format='tiff')
plt.tight_layout()
plt.show()
