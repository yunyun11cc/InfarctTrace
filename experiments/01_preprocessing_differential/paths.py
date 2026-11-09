"""阶段 01 的路径常量。

所有阶段 01 脚本都从这里取目录，不在各模块里重复拼接相对路径。
路径以本文件位置为基准、从阶段目录向上两级定位项目根目录，因此无论从项目根、
`experiments/` 还是阶段目录启动，读写位置都一致。

四类输出分别对应：
- ``matrices/`` 表达与计数矩阵
- ``tables/`` 样本、质控、差异和富集统计表
- ``figures/`` 科研图
- ``reports/`` 运行摘要、版本与校验记录
"""

from pathlib import Path

# 项目根目录：experiments/01_preprocessing_differential/paths.py 的上两级
ROOT = Path(__file__).resolve().parents[2]

# 阶段 01 输出根目录，与 experiments/ 下的阶段目录同名
OUT = ROOT / 'outputs/01_preprocessing_differential'

# 四类输出目录，按用途分开存放，便于索引和对照
MATRICES = OUT / 'matrices'
TABLES = OUT / 'tables'
FIGURES = OUT / 'figures'
REPORTS = OUT / 'reports'
