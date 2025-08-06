import pandas as pd
import numpy as np

# 示例数据（请用你自己的数据替换）
# 假设有两列：SampleID 和 Age
df = pd.read_csv('sample_information.csv')

# 第一步：每个年龄中，若样本数 ≥2，则随机取1个
filtered = (
    df.groupby('age')
    .filter(lambda x: len(x) >= 2)
    .groupby('age')
    .apply(lambda x: x.sample(1, random_state=8))
    .reset_index(drop=True)
)

# 第二步：添加年龄段信息
def get_age_group(age):
    if 20 <= age <= 29:
        return '20-29'
    elif 30 <= age <= 39:
        return '30-39'
    elif 40 <= age <= 49:
        return '40-49'
    elif 50 <= age <= 59:
        return '50-59'
    elif age >= 60:
        return '60+'
    else:
        return 'Other'

filtered['AgeGroup'] = filtered['age'].apply(get_age_group)

# 仅保留五个指定年龄段
filtered = filtered[filtered['AgeGroup'].isin(['20-29', '30-39', '40-49', '50-59', '60+'])]

# 第二轮筛选：每个年龄段中随机取5个（不足5个就取所有）
final_sample = (
    filtered.groupby('AgeGroup')
    .apply(lambda x: x.sample(min(3, len(x)), random_state=8))
    .reset_index(drop=True)
)

# 显示最终筛选结果
print(final_sample)

final_sample.to_excel('Split_135-15.xlsx', index=False)
