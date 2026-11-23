"""调度 GSE252027 的质控、标准化、差异、候选窗口、富集与绘图步骤。

输入为发现队列原始文件及已选步骤的上游产物，输出写入阶段一配置的目录。
"""

import argparse
from pathlib import Path
import sys


if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = 'experiments.01_preprocessing_differential'

from . import (
    step01_load_and_qc,
    step02_normalization_vst,
    step03_differential_expr,
    step04_window_selection,
    step05_enrichment,
    step06_visualize,
)


def run_model_steps(steps):
    """用同一次主拟合执行所选的标准化和差异步骤。"""
    meta, quant, annotation, inference = step02_normalization_vst.read_inputs()
    samples = meta[meta['主分析纳入'].eq('是')]
    print(f'拟合 {len(samples)} 个文库的 Python 主模型。', flush=True)
    dds = step02_normalization_vst.fit_time_model(quant, samples, inference)

    # 原流程先检验，再从同一个对象导出 VST，最后单独拟合敏感性模型。
    if 3 in steps:
        primary = step03_differential_expr.write_primary(dds, samples, annotation, inference)
    if 2 in steps:
        step02_normalization_vst.export_normalization(dds)
    if 3 in steps:
        sensitivity = meta[meta['敏感性分析纳入'].eq('是')]
        step03_differential_expr.write_sensitivity(
            dds,
            samples,
            sensitivity,
            quant,
            annotation,
            inference,
            primary,
        )


def main():
    """解析所选步骤并按依赖顺序执行，默认运行完整发现流程。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--steps',
        nargs='+',
        type=int,
        choices=range(1, 7),
        default=range(1, 7),
        help='仅运行指定步骤（1–6），复用已有的上游产物。',
    )
    steps = set(parser.parse_args().steps)

    if 1 in steps:
        step01_load_and_qc.run()
    if steps.intersection((2, 3)):
        run_model_steps(steps)
    for number, module in ((4, step04_window_selection), (5, step05_enrichment), (6, step06_visualize)):
        if number in steps:
            module.run()


if __name__ == '__main__':
    main()
