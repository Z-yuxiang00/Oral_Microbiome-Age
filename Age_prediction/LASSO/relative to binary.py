import pandas as pd

# 读取数据
data = pd.read_csv('150个样_未抽平数据.csv')

data.set_index(data.columns[0], inplace = True)

# 将非零值变为1，保持0不变
data = data.applymap(lambda x: 1 if x != 0 else 0)

data.to_csv('150个样_未抽平数据_二分类变量.csv', index=True)








