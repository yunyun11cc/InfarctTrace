"""按 D7 标签生成重复分层 CV 与 Bootstrap 的文库下标。

输入为与表达矩阵行对齐的标签和重试次数，输出训练、验证或 Bootstrap 下标。
"""

import numpy as np
from sklearn.model_selection import RepeatedStratifiedKFold

from ..config import FOLDS, SEED


def cv_splits(labels, repeats):
    """按标签生成固定种子的重复分层训练和验证下标。"""
    splitter = RepeatedStratifiedKFold(n_splits=FOLDS, n_repeats=repeats, random_state=SEED)
    return splitter.split(np.empty(len(labels)), labels)


def bootstrap_draws(labels, repeats):
    """在各类内有放回抽样，返回每轮排序后的文库下标。"""
    rng = np.random.default_rng(SEED)
    classes = [np.flatnonzero(labels == group) for group in (False, True)]
    draws = np.column_stack(
        [rng.choice(indices, size=(repeats, len(indices)), replace=True) for indices in classes]
    )
    return np.sort(draws, axis=1)
