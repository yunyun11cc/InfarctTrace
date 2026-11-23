"""阶段二的输入、输出、固定选择器参数与重采样设置。"""

from importlib import import_module
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE_ONE = import_module('experiments.01_preprocessing_differential.config').OUT
STAGE_TWO = ROOT / 'outputs/02_feature_selection'
AGGREGATION_DIR = 'aggregation'
BIO_VALIDATION_DIR = 'bio_validation'
FROZEN_CANDIDATES_FILE = 'frozen_candidates.tsv'
FROZEN_CANDIDATES = STAGE_TWO / 'union_22' / AGGREGATION_DIR / 'tables' / FROZEN_CANDIDATES_FILE
ASSETS = Path(__file__).resolve().parent / 'assets'
PATHWAY_MEMBERSHIP = ASSETS / 'ssgsea_pathway_membership.tsv'
PATHWAY_SOURCES = ASSETS / 'ssgsea_pathway_sources.json'

MODES = ('intersection_6', 'intersection_22', 'union_22')
MAIN_LIBRARY_COUNTS = {0: 4, 1: 4, 3: 2, 7: 4, 14: 4, 28: 4}
D3_D7_LIBRARY_COUNTS = {3: 2, 7: 4}
UNION_SHAPE = (22, 3907)
D7_POSITIVES = 4
TOP_K = 20
SEED = 252027
FOLDS = 4
CV_REPEATS = 50
BOOTSTRAPS = 1_000
RESAMPLING_JOBS = 8
AGGREGATION_JOBS = 6
METHODS = ('均值差', 'LASSO', 'ElasticNet', '随机森林')
CLASSIFIERS = METHODS[1:]
LOGISTIC_C = 1
MAX_ITER = 20_000
LASSO_L1_RATIO = 1.0
ELASTIC_NET_L1_RATIO = 0.5
FOREST_TREES = 512
FOREST_MAX_FEATURES = 'sqrt'
EFS_FEATURE_COUNTS = range(5, 31)
EFS_MAX_ITER = 1_000

CORE_GENES = ('C1qtnf3', 'Cthrc1', 'Postn', 'Wisp1', 'Cemip', 'Cxadr', 'Kif26b', 'Col24a1', 'Mrc1')
MARKERS = {
    'Nlrp3': 'Inflammation',
    'Casp1': 'Inflammation',
    'Il1b': 'Inflammation',
    'Tgfb1': 'Fibrosis',
    'Smad3': 'Fibrosis',
    'Col1a1': 'Fibrosis',
    'Acta2': 'Fibrosis',
}
SSGSEA_MIN_SIZE = 5
SSGSEA_MAX_SIZE = 500
SSGSEA_WEIGHT = 0.25
CORRELATION_PERMUTATIONS = 19_999
TOP_PATHWAY_GENES = 8
