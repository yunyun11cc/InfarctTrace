"""从 GEO 固定源文件分别整理 GSE111059 芯片和 GSE236374 RNA-seq。

读取 assets/sources.json、发现队列 ID 映射及阶段二冻结表。
输出清洗矩阵、质控、样本纳入、探针映射、候选覆盖与来源报告。
"""

from contextlib import closing
import csv
import gzip
import hashlib
from io import StringIO
import json
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tarfile
from tempfile import NamedTemporaryFile, TemporaryDirectory

import numpy as np
import pandas as pd
from pydeseq2.preprocessing import deseq2_norm

from .config import CHIP_GROUPS, DISCOVERY_MAPPING, FROZEN, ROOT, RNA_GROUPS, SOURCES, SQLITE_MEMBER


def sha256(path):
    """返回指定文件的 SHA-256，用于源文件与上游产物溯源。"""
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def restore_source(record):
    """按来源清单恢复本地源文件，校验 SHA-256 后返回路径。"""
    path = ROOT / 'datasets' / record['dataset'] / 'original' / record['filename']
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        with NamedTemporaryFile(dir=path.parent, suffix='.part', delete=False) as temp:
            staged = Path(temp.name)
        try:
            subprocess.run(
                [
                    'curl',
                    '--fail',
                    '--location',
                    '--silent',
                    '--show-error',
                    '--retry',
                    '2',
                    '--max-time',
                    '120',
                    '--output',
                    str(staged),
                    record['url'],
                ],
                check=True,
            )
            if sha256(staged) != record['sha256']:
                raise ValueError(f'Download digest mismatch: {path.name}')
            staged.replace(path)
        finally:
            staged.unlink(missing_ok=True)
    if sha256(path) != record['sha256']:
        raise ValueError(f'Input digest mismatch: {path}')
    return path


def read_series(path):
    """解析 GEO Series Matrix 的样本字段和数据表，返回字典与文本流。"""
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        lines = handle.read().splitlines()
    start = lines.index('!series_matrix_table_begin')
    stop = lines.index('!series_matrix_table_end')
    fields = {}
    for key, *values in csv.reader(
        (line for line in lines[:start] if line.startswith('!Sample_')), delimiter='\t'
    ):
        fields.setdefault(key.removeprefix('!Sample_'), []).append(values)
    return fields, StringIO('\n'.join(lines[start + 1 : stop]))


def series_samples(fields, platform, tissue, groups, title_parser):
    """用平台、样本描述和预期分组生成按 GSM 索引的样本表。"""
    (gsm,) = fields['geo_accession']
    (titles,) = fields['title']
    (organisms,) = fields['organism_ch1']
    (platforms,) = fields['platform_id']
    samples = pd.DataFrame(
        [
            title_parser(title, accession, features)
            for accession, title, organism, chip, *features in zip(
                gsm, titles, organisms, platforms, *fields['characteristics_ch1'], strict=True
            )
            if organism == 'Mus musculus' and chip == platform
        ]
    ).set_index('GSM')
    if len(samples) != len(gsm) or samples['Group'].value_counts().to_dict() != groups:
        raise ValueError(f'Unexpected samples or groups for {platform}')
    if samples.duplicated(['Group', 'Replicate']).any():
        raise ValueError(f'Duplicate group/replicate for {platform}')
    return samples.assign(Platform=platform, Tissue=tissue, Included=True)


def chip_sample(title, gsm, characteristics):
    """从芯片标题和特征字段解析分组、采样日及组内生物重复。"""
    matched = re.fullmatch(
        r'CF (uninjured|3D post-MI|7D post-MI|2W post-MI|4W post-MI), biological rep([1-3])', title
    )
    if (
        not matched
        or 'cell type: sorted TCF21 lineage-traced cardiac fibroblasts' not in characteristics
    ):
        raise ValueError(f'Unrecognized fibroblast sample: {gsm}, {title}')
    time, replicate = matched.groups()
    day = {'uninjured': 0, '3D post-MI': 3, '7D post-MI': 7, '2W post-MI': 14, '4W post-MI': 28}[time]
    return {
        'GSM': gsm,
        'Title': title,
        'Group': f'MI_D{day}' if day else 'uninjured',
        'Condition': 'MI' if day else 'uninjured',
        'Day': day,
        'Replicate': int(replicate),
        'Source_column': gsm,
    }


def rna_sample(title, gsm, characteristics):
    """交叉核对 RNA-seq 标题和 GEO 特征，返回分组及原计数表列名。"""
    matched = re.fullmatch(r'Heart, (Sham|Myocardial infarction, (7|28)days), rep([1-3])', title)
    if not matched:
        raise ValueError(f'Unrecognized RNA-seq sample: {gsm}, {title}')
    _, day, replicate = matched.groups()
    condition = 'MI' if day else 'Sham'
    time = f'{day} days after MI' if day else 'Sham'
    if not {'tissue: Heart', f'treatment: {condition}', f'time: {time}'}.issubset(characteristics):
        raise ValueError(f'GEO characteristics conflict with title: {gsm}, {title}')
    return {
        'GSM': gsm,
        'Title': title,
        'Group': f'MI_D{day}' if day else 'sham',
        'Condition': condition,
        'Day': int(day) if day else pd.NA,
        'Replicate': int(replicate),
        'Source_column': f'MI_{day + "d" if day else "sham"}_{replicate}',
    }


def align_samples(matrix, samples):
    """将基因或探针 × 文库矩阵按样本表重排列，并统一用 GSM 命名。"""
    source = samples['Source_column']
    if (
        not samples.index.is_unique
        or not source.is_unique
        or not matrix.columns.is_unique
        or set(source) != set(matrix.columns)
    ):
        raise ValueError('Matrix columns must match unique sample identifiers exactly')
    if matrix.empty or not matrix.index.is_unique or matrix.index.hasnans:
        raise ValueError('Matrix requires unique nonmissing row identifiers')
    return matrix.loc[:, source].set_axis(samples.index.to_list(), axis=1)


def read_chip_annotation(path):
    """只读查询固定芯片注释包中每个探针的候选 Entrez 和歧义标记。"""
    with TemporaryDirectory() as tmp, tarfile.open(path) as archive:
        database = Path(tmp) / 'annotation.sqlite'
        with archive.extractfile(SQLITE_MEMBER) as member, database.open('wb') as destination:
            shutil.copyfileobj(member, destination)
        with closing(sqlite3.connect(f'{database.as_uri()}?mode=ro', uri=True)) as connection:
            return pd.read_sql_query('SELECT probe_id, gene_id, is_multiple FROM probes', connection)


def aggregate_probes(matrix, annotation):
    """移除控制、歧义和未映射探针，返回 Entrez 中位表达与探针审计表。"""
    annotated = annotation.groupby('probe_id').agg(
        ids=('gene_id', lambda values: sorted(set(values.dropna()))),
        flagged=('is_multiple', 'max'),
    )
    audit = annotated.reindex(matrix.index)
    ids = audit['ids'].map(lambda value: value if isinstance(value, list) else [])
    counts = ids.str.len()
    control = matrix.index.str.startswith(('AFFX', 'control'))
    audit = pd.DataFrame(
        {
            'Entrez_candidates': ids.str.join('|'),
            'Status': np.select(
                [control, counts.gt(1) | audit['flagged'].eq(1), counts.eq(1)],
                ['control', 'ambiguous', 'mapped'],
                default='unmapped',
            ),
        },
        index=matrix.index,
    )
    selected = audit['Status'].eq('mapped')
    genes = matrix.loc[selected].groupby(audit.loc[selected, 'Entrez_candidates']).median()
    return genes.rename_axis('EntrezID'), audit


def validate_counts(matrix):
    """检查计数为有限的非负整数，返回 int64 矩阵。"""
    values = matrix.to_numpy(dtype=float)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.equal(values, np.floor(values)).all()
    ):
        raise ValueError('RNA-seq counts must be finite nonnegative integers')
    return matrix.astype('int64')


def write_table(frame, path):
    """按原 TSV 格式写出清洗表，gzip 表使用固定压缩时间戳和数值精度。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None
    frame.to_csv(path, sep='\t', compression=compression, float_format='%.12g')


def clean_chip(paths, output):
    """读取并对齐芯片信号，写出探针层和 Entrez 层表达、质控和纳入表。"""
    fields, table = read_series(paths['chip_matrix'])
    samples = series_samples(
        fields, 'GPL23038', 'sorted TCF21-lineage cardiac fibroblasts', CHIP_GROUPS, chip_sample
    )
    probe = align_samples(pd.read_csv(table, sep='\t', index_col='ID_REF'), samples)
    if not np.isfinite(probe.to_numpy()).all():
        raise ValueError('GSE111059 has nonfinite probe expression')
    genes, audit = aggregate_probes(probe, read_chip_annotation(paths['chip_annotation']))
    qc = pd.DataFrame(
        {'Probe_median': probe.median(), 'Probe_IQR': probe.quantile(0.75) - probe.quantile(0.25)}
    )
    write_table(probe, output / 'matrices/probe_sst_rma.tsv.gz')
    write_table(genes, output / 'matrices/gene_sst_rma.tsv.gz')
    write_table(samples, output / 'tables/sample_inclusion.tsv')
    write_table(qc.rename_axis('GSM'), output / 'tables/sample_qc.tsv')
    write_table(audit, output / 'tables/probe_mapping.tsv.gz')
    return genes, audit


def clean_rna(paths, output):
    """读取并校验 FeatureCounts，写出原始计数、归一化表达及质控。"""
    fields, _ = read_series(paths['rna_metadata'])
    samples = series_samples(fields, 'GPL24247', 'bulk heart tissue', RNA_GROUPS, rna_sample)
    counts = validate_counts(
        align_samples(
            pd.read_csv(paths['rna_counts'], sep='\t', index_col='Genesymbol'),
            samples,
        )
    ).rename_axis('Gene_symbol')
    positive = counts.gt(0)
    if not positive.all(axis=1).any():
        raise ValueError('Median of ratios needs genes with nonzero counts in all libraries')
    # 归一化函数以文库为行，写出基因 × 文库格式前再转置。
    normalized, factors = deseq2_norm(counts.T)
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError('Invalid GSE236374 size factors')
    qc = pd.DataFrame(
        {'Total_counts': counts.sum(), 'Genes_detected': positive.sum(), 'Size_factor': factors}
    )
    write_table(counts, output / 'matrices/raw_counts.tsv.gz')
    write_table(np.log2(normalized.T + 1), output / 'matrices/log2_normalized_counts.tsv.gz')
    write_table(samples, output / 'tables/sample_inclusion.tsv')
    write_table(qc.rename_axis('GSM'), output / 'tables/sample_qc.tsv')
    return counts


def candidate_coverage(genes, audit, counts):
    """按冻结 Ensembl ID 对齐注释，汇总两外部队列的覆盖情况。"""
    candidates = pd.read_csv(FROZEN, sep='\t', usecols=['基因ID', '基因符号']).set_index('基因ID')
    mapping = pd.read_csv(DISCOVERY_MAPPING, sep='\t', dtype={'EntrezID': str}).set_index('基因ID')
    coverage = candidates.join(mapping[['EntrezID']], validate='one_to_one')
    selected = audit.loc[audit['Status'].eq('mapped')].reset_index()
    probe_ids = selected.groupby('Entrez_candidates')['ID_REF'].agg('|'.join)
    coverage['GSE111059_probe_ids'] = coverage['EntrezID'].map(probe_ids)
    coverage['GSE111059_present'] = coverage['EntrezID'].isin(genes.index)
    coverage['GSE236374_present'] = coverage['基因符号'].isin(counts.index)
    return coverage.rename_axis('EnsemblID').rename(columns={'基因符号': 'Gene_symbol'})


def source_path(path):
    """记录来源路径，仓库内文件使用相对路径。"""
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


def analysis_table_sha256(path):
    """以 LF 换行计算上游分析表的 SHA-256，匹配 Git 中的内容校验。"""
    content = path.read_bytes().replace(b'\r\n', b'\n')
    return hashlib.sha256(content).hexdigest()


def run(output):
    """校验四份源文件，分别清洗两队列并写出候选覆盖和数据来源报告。"""
    records = json.loads(SOURCES.read_text(encoding='utf-8'))
    paths = {key: restore_source(record) for key, record in records.items()}
    chip, audit = clean_chip(paths, output / 'GSE111059')
    counts = clean_rna(paths, output / 'GSE236374')
    coverage = candidate_coverage(chip, audit, counts)
    write_table(coverage, output / 'tables/candidate_coverage.tsv')
    report = {
        'GSE111059': {
            'groups': CHIP_GROUPS,
            'probe_rows': len(audit),
            'gene_rows': len(chip),
            'probe_status': audit['Status'].value_counts().to_dict(),
            'scale': 'GEO author SST-RMA; median per Entrez ID',
            'candidate_genes_present': int(coverage['GSE111059_present'].sum()),
        },
        'GSE236374': {
            'groups': RNA_GROUPS,
            'count_rows': len(counts),
            'all_zero_genes': int(counts.eq(0).all(axis=1).sum()),
            'size_factor_reference_genes': int(counts.gt(0).all(axis=1).sum()),
            'scale': 'raw FeatureCounts and log2(median-ratio normalized counts + 1)',
            'candidate_genes_present': int(coverage['GSE236374_present'].sum()),
        },
        'sources': {
            key: dict(records[key], bytes=path.stat().st_size, path=source_path(path))
            for key, path in paths.items()
        },
        'candidate_inputs': {
            source_path(path): analysis_table_sha256(path) for path in (FROZEN, DISCOVERY_MAPPING)
        },
    }
    summary = output / 'reports/cleaning_summary.json'
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'GSE111059: {chip.shape}, GSE236374: {counts.shape}; {summary}')
