"""Random forest 选择器，输入为 D7 标签和文库 × 基因矩阵。

输出为与基因列对应的分数，以及该方法的全样本得分和文库扰动频率表、图。
"""

from sklearn.ensemble import RandomForestClassifier

from ..config import FOREST_MAX_FEATURES, FOREST_TREES, SEED

METHOD = '随机森林'


USES_SCALING = False


def fit(values, labels):
    """以原 VST 输入和 D7 标签拟合随机森林，返回已拟合模型。"""
    return RandomForestClassifier(
        n_estimators=FOREST_TREES,
        max_features=FOREST_MAX_FEATURES,
        class_weight='balanced',
        random_state=SEED,
    ).fit(values, labels)


def scores(model):
    """返回与输入基因列同序的 Gini 不纯度下降重要性。"""
    return model.feature_importances_


def export(ranking, output):
    """从联合排序表写出本方法的结果矩阵和扰动频率图。"""
    from ..common.method_output import write_method_outputs

    write_method_outputs(ranking, output, METHOD, 'm04_random_forest')
