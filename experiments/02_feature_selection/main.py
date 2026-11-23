"""调度阶段二筛选、重采样、EFS 聚合、相关分析与出图。

读取阶段一产物和固定通路资产，按筛选模式写入结果目录。
默认运行 22 文库并集的完整流程；已有上游产物时可以选择单独步骤。
"""

import argparse
from pathlib import Path
import sys


if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = 'experiments.02_feature_selection'

from .config import (
    AGGREGATION_JOBS,
    AGGREGATION_DIR,
    BIO_VALIDATION_DIR,
    BOOTSTRAPS,
    CV_REPEATS,
    MODES,
    RESAMPLING_JOBS,
    STAGE_TWO,
)


STEPS = (
    'selection',
    'methods',
    'resampling',
    'efs',
    'ssgsea',
    'markers',
    'table',
    'figures',
    'analysis-figures',
)
INTERSECTION_STEPS = ('selection', 'methods', 'figures')


def main():
    """解析筛选模式与步骤，依次执行所需模块并复用既有上游结果。"""
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--mode', choices=MODES, default='union_22')
    modes.add_argument('--union', dest='mode', action='store_const', const='union_22')
    modes.add_argument('--d3-d7-only', dest='mode', action='store_const', const='intersection_6')
    parser.add_argument('--steps', nargs='+', choices=STEPS)
    parser.add_argument('--output-dir', type=Path, help='筛选模式的输出目录。')
    parser.add_argument('--cv-repeats', type=int, default=CV_REPEATS)
    parser.add_argument('--bootstrap', type=int, default=BOOTSTRAPS)
    parser.add_argument('--jobs', type=int, help='覆盖重采样和 EFS 各自的默认并行进程数。')
    args = parser.parse_args()
    if min(args.cv_repeats, args.bootstrap) < 1 or (args.jobs is not None and args.jobs < 1):
        parser.error('重复次数与并行进程数必须为正数。')

    use_union = args.mode == 'union_22'
    allowed = STEPS if use_union else INTERSECTION_STEPS
    selected = set(args.steps or allowed)
    if not selected.issubset(allowed):
        parser.error('交集模式只支持 selection、methods 和 figures。')

    output = args.output_dir or STAGE_TWO / args.mode
    run_steps(args, selected, output, use_union)


def run_steps(args, selected, output, use_union):
    """按依赖顺序执行选定任务，将各工作流产物写入同一筛选模式。"""
    import matplotlib.pyplot as plt

    from .bio_validation import marker_coexpression, ssgsea_pathway
    from .common import figures, resampling, selection
    from .common.data_loader import read_ranking
    from .methods import m05_efs_rank_agg as efs

    resampling_dir = output / 'resampling'
    aggregation = output / AGGREGATION_DIR
    aggregation_tables = aggregation / 'tables'
    biology = output / BIO_VALIDATION_DIR
    biology_tables = biology / 'tables'
    if 'selection' in selected:
        selection.run_selection(output, d3_d7_only=args.mode == 'intersection_6', use_union=use_union)
    if 'methods' in selected:
        ranking = read_ranking(output)
        method_dir = output / 'methods'
        for method in selection.RANKERS:
            method.export(ranking, method_dir)
    if 'resampling' in selected:
        resampling.run_analysis(
            output,
            resampling_dir,
            args.cv_repeats,
            args.bootstrap,
            args.jobs or RESAMPLING_JOBS,
        )
    if 'efs' in selected:
        efs.run_analysis(output, aggregation_tables, args.cv_repeats, args.jobs or AGGREGATION_JOBS)
    if 'ssgsea' in selected:
        with plt.rc_context():
            ssgsea_pathway.run_top_candidates(output, biology)
        ssgsea_pathway.run_frozen_candidates(biology_tables)
    if 'markers' in selected:
        marker_coexpression.run(biology_tables)
    if 'table' in selected:
        efs.write_candidate_frequencies(output, aggregation_tables)
        efs.write_table(output, aggregation_tables, biology_tables)
    if 'figures' in selected:
        with plt.rc_context():
            figures.plot_selection(output)
        if use_union:
            with plt.rc_context():
                figures.plot_resampling(resampling_dir)
    if 'analysis-figures' in selected:
        with plt.rc_context():
            efs.plot_results(output, aggregation)
        with plt.rc_context():
            figures.plot_biology(biology)


if __name__ == '__main__':
    main()
