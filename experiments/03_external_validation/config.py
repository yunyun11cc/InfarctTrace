"""外部验证的队列、参考资产、冻结面板与输出目录。"""

from importlib import import_module
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ASSETS = Path(__file__).resolve().parent / 'assets'
SOURCES = ASSETS / 'sources.json'
OUT = ROOT / 'outputs/03_external_validation'
STAGE_ONE = import_module('experiments.01_preprocessing_differential.config')
STAGE_TWO = import_module('experiments.02_feature_selection.config')
DISCOVERY_MAPPING = STAGE_ONE.TABLES / 'enrichment_gene_mapping.tsv'
MOUSE_ANNOTATION = STAGE_ONE.CACHE / 'org.Mm.eg.sqlite'
FROZEN = STAGE_TWO.FROZEN_CANDIDATES
PATHWAY_MEMBERSHIP = STAGE_TWO.PATHWAY_MEMBERSHIP
GENES = list(STAGE_TWO.CORE_GENES)
FROZEN_LF_SHA256 = '7f7254627a79c5db2c664cf073122089d1ee46c9c1b16eec4921b036245fdb5b'

COHORTS = [('GSE111059', 'uninjured'), ('GSE236374', 'sham')]
CHIP_GROUPS = {'uninjured': 3, 'MI_D3': 3, 'MI_D7': 3, 'MI_D14': 2, 'MI_D28': 3}
RNA_GROUPS = {'sham': 3, 'MI_D7': 3, 'MI_D28': 3}
COHORT_SCALES = {'GSE111059': 'SST-RMA log2 signal', 'GSE236374': 'log2(CPM+1)'}
PRIMARY_GROUP = 'MI_D7'
MARKERS = ['Col1a1', 'Col3a1', 'Acta2', 'Tgfb1', 'Nlrp3', 'Il1b']
PATHWAYS = ['NLRP3 inflammasome', 'Extracellular matrix organization', 'Collagen formation']
SQLITE_MEMBER = 'clariomsmousetranscriptcluster.db/inst/extdata/clariomsmousetranscriptcluster.sqlite'
SEED = 252027
AUC_BOOTSTRAPS = 10_000
AUC_QUANTILES = [0.025, 0.975]
WELCH_GATE_ALPHA = 0.05
LOOCV_C = 1
SSGSEA_MIN_SIZE = 5
SSGSEA_MAX_SIZE = 500
SSGSEA_WEIGHT = 0.25

SINGLE_CELL_ARCHIVE = ROOT / 'datasets/GSE163465/original/GSE163465_RAW.tar'
SINGLE_CELL_ARCHIVE_SHA256 = '4ce64a86848b8af3c339a7ab645928d07fcc7344b24178b7f554f8d020243139'
SINGLE_CELL_LABELS = ASSETS / 'gse163465_cluster_labels.tsv'
SINGLE_CELL_LIBRARIES = [
    ('GSM4985022', 'Sham', 'Sham'),
    ('GSM4985023', 'Day_3', 'D3'),
    ('GSM4985024', 'Day_7', 'D7'),
]
SC_MIN_GENES = 300
SC_MAX_GENES = 6500
SC_MAX_MT_PCT = 20
SC_LIBRARY_SCALE = 10_000
SC_HIGHLY_VARIABLE_GENES = 2200
SC_N_PCS = 30
SC_N_NEIGHBORS = 15
SC_NEIGHBOR_JOBS = 4
SC_LEIDEN_RESOLUTION = 0.65
SC_LEIDEN_ITERATIONS = 2
SC_MARKER_SCORE_FLOOR = 0.1
SC_LABEL_THRESHOLD = 0.3

SC_MARKERS = {
    'Macrophage': ['Adgre1', 'Csf1r', 'C1qa', 'C1qb', 'Cd68', 'Mertk'],
    'Monocyte': ['Ly6c2', 'Ccr2', 'Plac8', 'Sell', 'Fcgr3'],
    'Neutrophil': ['Csf3r', 'Cxcr2', 'Retnlg', 'Mmp9', 'S100a8', 'Ngp'],
    'T cell': ['Cd3d', 'Cd3e', 'Trac', 'Cd247'],
    'NK cell': ['Ncr1', 'Klrb1c', 'Klrd1', 'Prf1', 'Nkg7'],
    'B cell': ['Cd79a', 'Cd79b', 'Ms4a1', 'Cd19'],
    'Dendritic cell': ['Flt3', 'Zbtb46', 'Xcr1', 'Clec9a', 'Itgax'],
    'Mast cell': ['Kit', 'Cpa3', 'Tpsb2', 'Mcpt4'],
    'Fibroblast-like': ['Pdgfra', 'Dcn', 'Lum', 'Col1a1', 'Dpt'],
}
SC_AUDIT_MARKERS = ['Ptprc', 'Acta2', 'Tagln', 'Col1a2', 'Tcf21', 'H2-Ab1']
