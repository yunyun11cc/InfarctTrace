"""导入 GSE252027，整理基因矩阵与样本纳入表。"""
import gzip
import re
import tarfile

import numpy as np
import pandas as pd

from paths import ROOT, MATRICES, TABLES


def save_gene_table(frame, path):
    # 导出时统一用基因 ID 做索引，保留 GSM 作为矩阵列名。
    frame.to_csv(path, sep='\t', index_label='基因ID', encoding='utf-8-sig')


def annotation():
    """从与定量参考一致的注释版本读取基因名称和类型。"""
    genes = {}
    source = ROOT / 'datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz'
    with gzip.open(source, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            gid = attrs['gene_id'].split('.')[0]
            genes[gid] = dict(symbol=attrs.get('gene_name', gid),
                              biotype=attrs.get('gene_biotype', 'unknown'))
    gene_info = pd.DataFrame.from_dict(genes, orient='index')
    save_gene_table(gene_info.rename(columns={'symbol': '基因符号', 'biotype': '基因类型'}),
                    MATRICES / 'gene_annotation.tsv.gz')
    return gene_info


def load_fb():
    """按 GSM 读取 FB 文库，检查数值并按基因 ID 对齐。"""
    meta = pd.read_csv(ROOT / 'datasets/GSE252027/metadata/samples_FB.csv').set_index('gsm')
    matrices = {name: {} for name in ('TPM', 'NumReads', 'EffectiveLength')}
    qc, excluded = [], None
    with tarfile.open(ROOT / 'datasets/GSE252027/original/GSE252027_RAW.tar') as archive:
        entries = {m.name.split('_')[0]: m for m in archive.getmembers() if m.isfile()}
        for gsm, row in meta.iterrows():
            with gzip.open(archive.extractfile(entries[gsm]), 'rt') as stream:
                frame = pd.read_csv(stream, sep='\t').set_index('Name')
            if not np.isfinite(frame.to_numpy()).all() or (frame.to_numpy() < 0).any():
                raise ValueError(f'Invalid quantification: {gsm}')
            gene_rows = frame.index.str.match(r'^ENSMUSG\d+(\.\d+)?$')
            if excluded is None:
                excluded = frame.loc[~gene_rows].copy()
                labels = {'Length': '长度', 'EffectiveLength': '有效长度', 'NumReads': '估计片段数'}
                for example, name in [(excluded, 'excluded_non_gene_rows_example.tsv'),
                                      (frame.head(12), 'raw_quant_first12rows.tsv')]:
                    example.rename(columns=labels).to_csv(
                        TABLES / name, sep='\t', index_label='特征编号', encoding='utf-8-sig')
            # 文件混有转录本行，基因级分析仅纳入 ENSMUSG。
            selected = frame.loc[gene_rows].copy()
            selected.index = selected.index.str.replace(r'\.\d+$', '', regex=True)
            if not selected.index.is_unique:
                raise ValueError('Version stripping produced duplicate gene IDs')
            if (selected.EffectiveLength <= 0).any():
                raise ValueError('Nonpositive effective length')
            for name in matrices:
                matrices[name][gsm] = selected[name]
            qc.append(dict(gsm=gsm, day=int(row.day), replicate=int(row.replicate),
                           feature_rows=len(frame), gene_rows=int(gene_rows.sum()),
                           excluded_non_gene_rows=int((~gene_rows).sum()),
                           estimated_reads_all=float(frame.NumReads.sum()),
                           estimated_reads_retained=float(selected.NumReads.sum()),
                           tpm_sum_all=float(frame.TPM.sum()),
                           genes_tpm_ge1=int((selected.TPM >= 1).sum())))
    matrices = {key: pd.DataFrame(value) for key, value in matrices.items()}
    if any(frame.isna().any().any() for frame in matrices.values()):
        raise ValueError('Gene sets differ between FB libraries')
    sample_labels = {'title': '标题', 'source': '样本来源', 'organism': '物种', 'cell_type': '细胞类型',
                     'treatment': '原始处理说明', 'supplementary_files': '定量文件地址', 'day': '时间点（天）',
                     'replicate': '组内重复编号', 'cell_type_abbreviation': '细胞类型简称',
                     'pooled_mice_per_library': '每文库记录的小鼠数', 'paper_qc_status': '原文质控状态'}
    sample_table = meta.rename(columns=sample_labels).copy()
    sample_table['原文质控状态'] = sample_table['原文质控状态'].replace({
        'not_assessed_here': '本表未评估', 'pending_exact_author_exclusions': '待核对作者的样本排除记录'})
    sample_table.to_csv(TABLES / 'discovery_samples.csv', index_label='GSM编号', encoding='utf-8-sig')
    qc_labels = {'gsm': 'GSM编号', 'day': '时间点（天）', 'replicate': '组内重复编号',
                 'feature_rows': '全部特征数', 'gene_rows': '基因行数', 'excluded_non_gene_rows': '排除的非基因行数',
                 'estimated_reads_all': '全部估计片段数', 'estimated_reads_retained': '保留基因的估计片段数',
                 'tpm_sum_all': '全表TPM总和', 'genes_tpm_ge1': 'TPM大于等于1的基因数'}
    pd.DataFrame(qc).rename(columns=qc_labels).to_csv(TABLES / 'discovery_qc.csv', index=False, encoding='utf-8-sig')
    for name, frame in matrices.items():
        save_gene_table(frame, MATRICES / f'FB_{name}.tsv.gz')
    return meta, matrices, len(excluded)


def length_scaled_counts(tpm, effective_length, estimated_reads):
    """按 tximport 的 lengthScaledTPM 算法生成计数尺度输入。"""
    # 每个基因用所纳入文库的平均有效长度，再恢复每列原估计片段总量。
    unscaled = tpm.mul(effective_length.mean(axis=1), axis=0)
    counts = unscaled.mul(estimated_reads.sum(axis=0) / unscaled.sum(axis=0), axis=1)
    if not np.allclose(counts.sum(), estimated_reads.sum(), rtol=1e-10):
        raise ValueError('Length scaling changed library totals')
    return counts.round().astype('int64')


def write_inclusion(meta):
    # 论文记录 D3 排除两个 FB 异常重复，作者脚本保留 D3r1 与 D3r4。
    d3_excluded = (meta.day == 3) & meta.replicate.isin([2, 3])
    d28_uncertain = (meta.day == 28) & (meta.replicate == 4)
    if set(meta.index[d3_excluded]) != {'GSM7992335', 'GSM7992336'} or set(meta.index[d28_uncertain]) != {'GSM7992349'}:
        raise ValueError('样本编号与事先定义的质控规则不符。')
    reason = np.select([d3_excluded, d28_uncertain],
                       ['论文报告 D3 排除两份 FB 文库，作者脚本保留 r1 和 r4',
                        '论文未明确排除该 FB 文库，主分析保留并另做排除敏感性分析'],
                       default='纳入主分析及敏感性分析')
    decision = pd.DataFrame({'GSM编号': meta.index, '时间点（天）': meta.day.to_numpy(),
                             '组内重复编号': meta.replicate.to_numpy(),
                             '主分析纳入': np.where(d3_excluded, '否', '是'),
                             '敏感性分析纳入': np.where(d3_excluded | d28_uncertain, '否', '是'),
                             '依据': reason})
    decision.to_csv(TABLES / 'sample_inclusion.csv', index=False, encoding='utf-8-sig')


def main():
    for directory in (MATRICES, TABLES):
        directory.mkdir(parents=True, exist_ok=True)
    annotation()
    meta, _, _ = load_fb()
    write_inclusion(meta)
    print('已整理 24 个提交文库。主分析纳入 22 个，D28 排除敏感性分析纳入 21 个。')


if __name__ == '__main__':
    main()
