import pandas as pd
from scipy.stats import ttest_rel, wilcoxon, normaltest

# 读取Excel文件
file_path = "SNR_all.xlsx"
df = pd.read_excel(file_path)

# 需要比较的四列
columns_to_compare = ["binarized data", "relative abundance", "CLR", "log2"]

# 保存检验结果
results = []

# 所有两两组合
for i in range(len(columns_to_compare)):
    for j in range(i + 1, len(columns_to_compare)):
        col1 = columns_to_compare[i]
        col2 = columns_to_compare[j]

        # 提取非缺失值行（确保成对）
        paired_data = df[[col1, col2]].dropna()

        # 计算差值
        diff = paired_data[col1] - paired_data[col2]

        # 使用 normaltest 检验正态性（适用于 N > 20）
        normal_stat, normal_p = normaltest(diff)

        if normal_p >= 0.05:
            # 若差值近似正态分布，使用配对t检验
            stat, p_val = ttest_rel(paired_data[col1], paired_data[col2])
            test_type = "Paired t-test"
        else:
            # 若非正态，使用Wilcoxon配对秩检验
            stat, p_val = wilcoxon(paired_data[col1], paired_data[col2])
            test_type = "Wilcoxon signed-rank"

        note = "✓ p < 0.05" if p_val < 0.05 else "ns (p ≥ 0.05)"

        results.append({
            "Comparison": f"{col1} vs {col2}",
            "Test": test_type,
            "Normality_p": round(normal_p, 4),
            "Stat": round(stat, 4),
            "p-value": f"{p_val:.4e}",
            "Note": note
        })

# 输出结果表格
results_df = pd.DataFrame(results)
print(results_df)

# 保存为 Excel 文件
results_df.to_excel("SNR_pairwise_test_results.xlsx", index=False)