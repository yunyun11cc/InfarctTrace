"""调度两外部队列的清洗、冻结候选评估与绘图。

默认写入阶段三的清洗和评估子目录。单细胞定位与跨队列汇总可单独运行。
"""

import argparse
from pathlib import Path
import sys


if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = 'experiments.03_external_validation'

from .config import OUT


BULK_STEPS = ('preprocess', 'evaluate', 'visualize')
STEPS = (*BULK_STEPS, 'single-cell', 'summary')


def main():
    """解析步骤与输出目录，从已有产物或源数据运行外部验证。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', nargs='+', choices=STEPS, default=BULK_STEPS)
    parser.add_argument('--output-dir', type=Path, default=OUT, help='阶段三的输出根目录。')
    parser.add_argument('--preprocessing-dir', type=Path, help='队列清洗输出及评估输入目录。')
    args = parser.parse_args()

    preprocessing = args.preprocessing_dir or args.output_dir / 'preprocessing'
    run_steps(set(args.steps), preprocessing, args.output_dir)


def run_steps(steps, preprocessing, output):
    """依赖顺序执行所选步骤，仅在请求时加载 Scanpy。"""
    evaluation = output / 'evaluation'
    single_cell_output = output / 'single_cell'
    if 'preprocess' in steps:
        from . import step01_preprocess_cohorts

        step01_preprocess_cohorts.run(preprocessing)
    if 'evaluate' in steps:
        from . import step02_blind_evaluation

        step02_blind_evaluation.run(preprocessing, evaluation)
    if 'visualize' in steps:
        from . import step03_visualize

        step03_visualize.run(preprocessing, evaluation)
    if 'single-cell' in steps:
        from .supplementary import single_cell

        single_cell.run(single_cell_output)
    if 'summary' in steps:
        from .supplementary import summarize

        summarize.run(evaluation, single_cell_output, output / 'summary')


if __name__ == '__main__':
    main()
