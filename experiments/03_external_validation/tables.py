"""读写外部评估、单细胞和汇总步骤所用的 TSV。"""

import pandas as pd


def read(path, **kwargs):
    """读取一个验证 TSV，解析参数由调用方指定。"""
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, name, output, **kwargs):
    """以原验证表格式保存 DataFrame，保留原缺失值标记与列序。"""
    frame.to_csv(output / 'tables' / f'{name}.tsv', sep='\t', index=False, na_rep='NA', **kwargs)
