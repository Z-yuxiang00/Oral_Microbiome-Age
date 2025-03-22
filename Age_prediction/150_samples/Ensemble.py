import os
import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_absolute_error
from sklearn.neural_network import MLPClassifier

os.chdir(r"pathway")
data = pd.read_excel(
    r"Input_150samples+105ASV.xlsx",
    header=0
)
X = data.drop('age', axis=1)
y = data['age']

def age_groups(y, groups):
    bins = np.concatenate(([1], groups, [80]))
    return np.digitize(y, bins) - 1

loo = LeaveOneOut()
N_values = [48]
classifiers = []

for N in N_values:
    for i in range(1, N + 1):
        if i == 1:
            groups = list(range(N + 1, 80, N))
        else:
            groups = list(range(i, 81, N))
        classifiers.append(groups)

predicted_ages = []
actual_ages = []

for train_index, test_index in loo.split(X):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

    models = []
    for groups in classifiers:
        y_train_groups = age_groups(y_train, groups)

        clf = MLPClassifier(hidden_layer_sizes=(30, 15), max_iter=1000, random_state=8)
        clf.fit(X_train, y_train_groups)
        models.append((clf, groups))

    votes_test = np.zeros((80,), dtype=int)
    for clf, groups in models:
        y_pred_test_groups = clf.predict(X_test)
        bins = np.concatenate(([1], groups, [80]))
        for pred_group in y_pred_test_groups:
            vote_range = range(bins[pred_group], bins[pred_group + 1])
            for age in vote_range:
                votes_test[age - 1] += 1

    final_prediction = np.argmax(votes_test) + 1
    predicted_ages.append(final_prediction)
    actual_ages.append(y_test.values[0])

mae = mean_absolute_error(actual_ages, predicted_ages)
print(f"Mean Absolute Error: {mae}")

results_combined = pd.DataFrame({'Actual Age': actual_ages, 'Predicted Age': predicted_ages})
results_combined.to_excel(
    r"Result.xlsx",
    index=False
)

print("预测完成。合并后的结果已保存。")
