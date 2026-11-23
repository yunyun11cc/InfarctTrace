"""发现队列数据、注释快照、计算参数和输出目录。"""

from os import environ
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DISCOVERY = ROOT / 'datasets/GSE252027'
FB_SAMPLES = DISCOVERY / 'metadata/samples_FB.csv'
FB_ARCHIVE = DISCOVERY / 'original/GSE252027_RAW.tar'
GTF = ROOT / 'datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz'
CACHE = ROOT / '.tools'
KEGG_DIRECTORY = ROOT / 'datasets/annotation/kegg_mmu'
GO_SOURCE = 'https://bioconductor.org/packages/3.22/data/annotation/src/contrib/'
GO_PACKAGES = {'org.Mm.eg.db': 'org.Mm.eg.sqlite', 'GO.db': 'GO.sqlite'}
GO_VERSION = '3.22.0'
KEGG_SOURCES = {
    'pathway_names.tsv': 'https://rest.kegg.jp/list/pathway/mmu',
    'pathway_genes.tsv': 'https://rest.kegg.jp/link/pathway/mmu',
    'release.txt': 'https://rest.kegg.jp/info/pathway',
}

OUT = Path(
    environ.get(
        'INFARCTTRACE_STAGE1_OUTPUT',
        ROOT / 'outputs/01_preprocessing_differential',
    )
)
MATRICES = OUT / 'matrices'
TABLES = OUT / 'tables'
FIGURES = OUT / 'figures'
REPORTS = OUT / 'reports'

DAYS = [0, 1, 3, 7, 14, 28]
D3_EXCLUSIONS = {'GSM7992335', 'GSM7992336'}
D28_SENSITIVITY_EXCLUSIONS = {'GSM7992349'}
MIN_COUNT = 10
MIN_LIBRARIES = 2
INFERENCE_CPUS = 2
ALPHA = 0.05
MIN_ABS_LOG2FC = 1
LRT_BISECTION_STEPS = 48
ORA_SIZE_RANGE = (10, 500)
GSEA_SIZE_RANGE = (15, 500)
GSEA_MIN_RANKED_GENES = 1000
GSEA_PERMUTATIONS = 1000
GSEA_SAMPLE_SIZE = 101
GSEA_THREADS = 2
GSEA_EPS = 1e-10
GSEA_SEED = 252028
