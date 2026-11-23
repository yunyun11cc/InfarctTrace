"""Elastic Net 选择器，输入为 D7 标签和文库 × 基因矩阵。

输出为与基因列对应的分数，以及该方法的全样本得分和文库扰动频率表、图。
"""

from sklearn.linear_model import LogisticRegression

from ..config import ELASTIC_NET_L1_RATIO, LOGISTIC_C, MAX_ITER, SEED

METHOD = 'ElasticNet'


USES_SCALING = True


def fit(values, labels):
    """用已标准化的训练文库拟合 Elastic Net 逻辑回归，返回已拟合模型。"""
    model = LogisticRegression(
        solver='saga',
        l1_ratio=ELASTIC_NET_L1_RATIO,
        C=LOGISTIC_C,
        class_weight='balanced',
        random_state=SEED,
        max_iter=MAX_ITER,
    ).fit(values, labels)
    if model.n_iter_.max() >= MAX_ITER:
        raise RuntimeError('ElasticNet 逻辑回归未收敛')
    return model


def scores(model):
    """返回与输入基因列同序的带符号系数。"""
    return model.coef_[0]


def export(ranking, output):
    """从联合排序表写出本方法的结果矩阵和扰动频率图。"""
    from ..common.method_output import write_method_outputs

    write_method_outputs(ranking, output, METHOD, 'm03_elastic_net')
