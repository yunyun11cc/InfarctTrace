"""阶段 01 的项目路径与输出分类。"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/01_preprocessing_differential'
MATRICES = OUT / 'matrices'
TABLES = OUT / 'tables'
FIGURES = OUT / 'figures'
REPORTS = OUT / 'reports'

