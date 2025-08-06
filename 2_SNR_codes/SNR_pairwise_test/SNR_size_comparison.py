import pandas as pd
from itertools import combinations

# 读取数据
df = pd.read_excel("SNR_all.xlsx")

# 需要比较的列
cols = ["binarized data", "relative abundance", "CLR", "log2"]

# 保存结果
summary = []

# 两两比较
for col1, col2 in combinations(cols, 2):
    greater = (df[col1] > df[col2]).sum()
    less = (df[col1] < df[col2]).sum()
    equal = (df[col1] == df[col2]).sum()

    summary.append({
        "Comparison": f"{col1} vs {col2}",
        f"{col1} > {col2}": greater,
        f"{col1} < {col2}": less,
        f"{col1} = {col2}": equal
    })

# 转为DataFrame
summary_df = pd.DataFrame(summary)

# 显示结果
print(summary_df)

# 可选：保存结果为Excel
summary_df.to_excel("SNR_size_comparison_result.xlsx", index=False)
