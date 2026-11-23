"""VST mean difference 选择器，输入为 D7 标签和文库 × 基因矩阵。

输出为与基因列对应的分数，以及该方法的全样本得分和文库扰动频率表、图。
"""

METHOD = '均值差'


def score(vst, labels):
    """返回每个基因的 D7 与非 D7 文库 VST 均值之差。

    vst 为文库 × 基因数组，labels 为同序的 D7 布尔掩码。
    """
    d7_mean = vst[labels].mean(axis=0)
    other_mean = vst[~labels].mean(axis=0)
    return d7_mean - other_mean


def export(ranking, output):
    """从联合排序表写出本方法的结果矩阵和扰动频率图。"""
    from ..common.method_output import write_method_outputs

    write_method_outputs(ranking, output, METHOD, 'm01_vst_mean_diff')
