import os
import random
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelEncoder
import warnings
warnings.filterwarnings("ignore")

# 固定随机种子
random.seed(8)
np.random.seed(8)

# 设置工作路径
# os.chdir(r"/home/smy/ZYX/调参-服务器")

# 读取数据
train_data = pd.read_excel("年龄预测_135_Training_二分类_阈值0.5.xlsx")
test_data = pd.read_excel("年龄预测_15_Testing_二分类_阈值0.5.xlsx")

X_train = train_data.drop('age', axis=1)
y_train = train_data['age']
X_test = test_data.drop('age', axis=1)
y_test = test_data['age']

# 分组函数
def age_groups(y, groups):
    bins = np.concatenate(([1], groups, [100]))
    return np.digitize(y, bins) - 1

# ✅ 固定参数
params = {
    'n_estimators': 200,
    'max_depth': 3,
    'learning_rate': 0.1,
    'subsample': 1.0,
    'colsample_bytree': 1.0,
    'gamma': 0,
    'min_child_weight': 1,
    'reg_alpha': 0,
    'reg_lambda': 5
}

# ✅ 固定 N
N = 32
classifiers = []
for i in range(1, N + 1):
    if i == 1:
        groups = list(range(N + 1, 100, N))
    else:
        groups = list(range(i, 101, N))
    classifiers.append(groups)

# 拟合模型
final_models = []
for groups in classifiers:
    y_train_grouped = age_groups(y_train, groups)
    encoder = LabelEncoder()
    y_train_encoded = encoder.fit_transform(y_train_grouped)

    clf = XGBClassifier(
        use_label_encoder=False,
        eval_metric='logloss',
        random_state=8,
        **params
    )
    clf.fit(X_train, y_train_encoded)
    final_models.append((clf, groups, encoder))

# 测试集预测
predicted_ages_test = []
for i in range(len(X_test)):
    x_single = X_test.iloc[[i]]
    votes = np.zeros(100, dtype=int)
    for clf, groups, encoder in final_models:
        pred_group = clf.predict(x_single)
        pred_group = encoder.inverse_transform(pred_group)
        bins = np.concatenate(([1], groups, [100]))
        vote_range = range(bins[pred_group[0]], bins[pred_group[0] + 1])
        for age in vote_range:
            votes[age - 1] += 1
    final_pred = np.argmax(votes) + 1
    predicted_ages_test.append(final_pred)

# 计算评估指标
mae_test = mean_absolute_error(y_test, predicted_ages_test)
rmse_test = mean_squared_error(y_test, predicted_ages_test, squared=False)  # RMSE
r2_test = r2_score(y_test, predicted_ages_test)  # R²

print(f"N={N} 测试集 MAE: {mae_test:.4f}")
print(f"N={N} 测试集 RMSE: {rmse_test:.4f}")
print(f"N={N} 测试集 R²: {r2_test:.4f}")

# 输出预测结果
test_result_df = pd.DataFrame({
    'Actual Age': y_test.values,
    'Predicted Age': predicted_ages_test
})

# 添加评估指标行
test_result_df.loc[len(test_result_df.index)] = ["MAE", mae_test]
test_result_df.loc[len(test_result_df.index)] = ["RMSE", rmse_test]
test_result_df.loc[len(test_result_df.index)] = ["R2", r2_test]

# 保存结果
test_result_df.to_excel("预测结果_age range_1-100.xlsx", index=False)