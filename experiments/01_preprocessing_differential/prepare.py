"""导入 GSE252027，整理基因矩阵与样本纳入表。

本模块是阶段 01 的第一步，只做数据整理，不拟合任何模型。它读取三类输入：

1. ``datasets/GSE252027/original/GSE252027_RAW.tar``：GEO 原始定量归档，
   每个 GSM 对应一个 ``.tsv.gz``，来自 Salmon/tximport 流程。
2. ``datasets/GSE252027/metadata/samples_FB.csv``：24 个成纤维细胞（FB）文库
   的样本表，含时间点 ``day`` 与组内重复编号 ``replicate``。
3. ``datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz``：与定量参考一致的
   Ensembl GRCm38 release 87 注释，提供基因符号和基因类型。

输出的四类文件：

- ``matrices/FB_{TPM,NumReads,EffectiveLength}.tsv.gz``：按基因 ID 对齐的
  定量矩阵，列名是 GSM。
- ``matrices/gene_annotation.tsv.gz``：基因 ID → 基因符号/基因类型。
- ``tables/discovery_samples.csv``、``tables/discovery_qc.csv``：样本说明与
  逐文库质控统计。
- ``tables/sample_inclusion.csv``：主分析与敏感性分析各自的纳入决定。
- ``tables/excluded_non_gene_rows_example.tsv``、``raw_quant_first12rows.tsv``：
  展示原始定量文件里非基因行长什么样，便于复核过滤规则。
"""

import gzip
import re
import tarfile

import numpy as np
import pandas as pd

from paths import ROOT, MATRICES, TABLES


# ---------------------------------------------------------------------------
# 1. 通用写出函数
# ---------------------------------------------------------------------------

def save_gene_table(frame, path):
    """按基因 ID 作行索引写出 TSV。

    导出时统一用基因 ID 做索引，保留 GSM 作为矩阵列名。UTF-8 BOM 让 Excel
    直接打开时中文表头不乱码；``.gz`` 后缀由 pandas 自动识别压缩。
    """
    frame.to_csv(path, sep='\t', index_label='基因ID', encoding='utf-8-sig')


# ---------------------------------------------------------------------------
# 2. 基因注释
# ---------------------------------------------------------------------------

def annotation():
    """从与定量参考一致的注释版本读取基因名称和类型。

    GTF 是逐行制表符分隔的文本，第 3 列（``fields[2]``）为特征类型，第 9 列
    （``fields[8]``）是形如 ``gene_id "ENSMUSG..."`` 的属性串。

    返回的 DataFrame 保留英文列名 ``symbol``/``biotype``，供下游 join 使用；
    写入磁盘的副本才改成中文表头。
    """
    genes = {}
    source = ROOT / 'datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz'

    with gzip.open(source, 'rt') as stream:
        for line in stream:
            # 以 # 开头的行是 GTF 头部说明，不含注释记录
            if line.startswith('#'):
                continue

            fields = line.rstrip().split('\t')
            # 只保留 gene 级特征，跳过 exon/transcript 等重复出现的行
            if fields[2] != 'gene':
                continue

            # 把 key "value" 形式的属性解析成字典
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            # 去掉 Ensembl ID 的版本后缀（如 .87），与定量矩阵对齐
            gid = attrs['gene_id'].split('.')[0]

            genes[gid] = dict(symbol=attrs.get('gene_name', gid),
                              biotype=attrs.get('gene_biotype', 'unknown'))

    # orient='index' 让基因 ID 成为行索引
    gene_info = pd.DataFrame.from_dict(genes, orient='index')

    save_gene_table(gene_info.rename(columns={'symbol': '基因符号', 'biotype': '基因类型'}),
                    MATRICES / 'gene_annotation.tsv.gz')

    return gene_info


# ---------------------------------------------------------------------------
# 3. 读取 FB 文库定量结果
# ---------------------------------------------------------------------------

def load_fb():
    """按 GSM 读取 FB 文库，检查数值并按基因 ID 对齐。

    返回 ``(meta, matrices, excluded_count)``：

    - ``meta``：样本表，索引为 GSM，含 ``day`` 和 ``replicate``。
    - ``matrices``：``{'TPM': DataFrame, 'NumReads': ..., 'EffectiveLength': ...}``，
      行是去掉版本号的基因 ID，列是 GSM。
    - ``excluded_count``：第一个文库中被排除的非基因特征行数。

    读入时会做三类硬性检查，不满足就直接报错，避免脏数据流入模型：
    数值必须是有限非负数、去掉版本号后基因 ID 不得重复、有效长度必须为正。
    """
    meta = pd.read_csv(ROOT / 'datasets/GSE252027/metadata/samples_FB.csv').set_index('gsm')

    # 三个字段分别收集，最后统一转成 DataFrame；空字典的键顺序即列顺序
    matrices = {name: {} for name in ('TPM', 'NumReads', 'EffectiveLength')}
    qc, excluded = [], None

    with tarfile.open(ROOT / 'datasets/GSE252027/original/GSE252027_RAW.tar') as archive:
        # 归档内文件名形如 GSM7992326_xxx.tsv.gz，用下划线前的 GSM 建索引
        entries = {m.name.split('_')[0]: m for m in archive.getmembers() if m.isfile()}

        for gsm, row in meta.iterrows():
            # 直接从归档流读取，不落地解压；第一列 Name 是 Salmon 的特征编号
            with gzip.open(archive.extractfile(entries[gsm]), 'rt') as stream:
                frame = pd.read_csv(stream, sep='\t').set_index('Name')

            # 计数、TPM、有效长度都应为有限非负数
            if not np.isfinite(frame.to_numpy()).all() or (frame.to_numpy() < 0).any():
                raise ValueError(f'Invalid quantification: {gsm}')

            # 每个文件混有转录本行（ENSMUST...），基因级分析只保留 ENSMUSG 行
            gene_rows = frame.index.str.match(r'^ENSMUSG\d+(\.\d+)?$')

            # 只导出第一个文库的示例，用来说明被排除的行与原始文件结构
            if excluded is None:
                excluded = frame.loc[~gene_rows].copy()
                labels = {'Length': '长度', 'EffectiveLength': '有效长度', 'NumReads': '估计片段数'}
                for example, name in [(excluded, 'excluded_non_gene_rows_example.tsv'),
                                      (frame.head(12), 'raw_quant_first12rows.tsv')]:
                    example.rename(columns=labels).to_csv(
                        TABLES / name, sep='\t', index_label='特征编号', encoding='utf-8-sig')

            # 文件混有转录本行，基因级分析仅纳入 ENSMUSG
            selected = frame.loc[gene_rows].copy()

            # GTF 里的基因 ID 不带版本号，这里同步去掉，保证三张矩阵与注释可对齐
            selected.index = selected.index.str.replace(r'\.\d+$', '', regex=True)
            if not selected.index.is_unique:
                raise ValueError('Version stripping produced duplicate gene IDs')

            # 有效长度为 0 会让 TPM 反推计数时除零，提前拦截
            if (selected.EffectiveLength <= 0).any():
                raise ValueError('Nonpositive effective length')

            for name in matrices:
                matrices[name][gsm] = selected[name]

            # 逐文库质控：既能核对 Salmon 汇总规模，也能比较过滤前后的保留量
            qc.append(dict(gsm=gsm, day=int(row.day), replicate=int(row.replicate),
                           feature_rows=len(frame), gene_rows=int(gene_rows.sum()),
                           excluded_non_gene_rows=int((~gene_rows).sum()),
                           estimated_reads_all=float(frame.NumReads.sum()),
                           estimated_reads_retained=float(selected.NumReads.sum()),
                           tpm_sum_all=float(frame.TPM.sum()),
                           genes_tpm_ge1=int((selected.TPM >= 1).sum())))

    # 每个字典按键（GSM）顺序对齐成“基因 × 文库”矩阵
    matrices = {key: pd.DataFrame(value) for key, value in matrices.items()}

    # 任一 NaN 都意味着某个文库缺少某些基因行，与“同一套基因”假设冲突
    if any(frame.isna().any().any() for frame in matrices.values()):
        raise ValueError('Gene sets differ between FB libraries')

    # 样本表转中文表头后再写出，同时把状态码翻译成便于阅读的说明
    sample_labels = {'title': '标题', 'source': '样本来源', 'organism': '物种', 'cell_type': '细胞类型',
                     'treatment': '原始处理说明', 'supplementary_files': '定量文件地址', 'day': '时间点（天）',
                     'replicate': '组内重复编号', 'cell_type_abbreviation': '细胞类型简称',
                     'pooled_mice_per_library': '每文库记录的小鼠数', 'paper_qc_status': '原文质控状态'}
    sample_table = meta.rename(columns=sample_labels).copy()
    sample_table['原文质控状态'] = sample_table['原文质控状态'].replace({
        'not_assessed_here': '本表未评估', 'pending_exact_author_exclusions': '待核对作者的样本排除记录'})
    sample_table.to_csv(TABLES / 'discovery_samples.csv', index_label='GSM编号', encoding='utf-8-sig')

    # 质控表同样使用中文表头，方便与结果解读文档逐列对照
    qc_labels = {'gsm': 'GSM编号', 'day': '时间点（天）', 'replicate': '组内重复编号',
                 'feature_rows': '全部特征数', 'gene_rows': '基因行数', 'excluded_non_gene_rows': '排除的非基因行数',
                 'estimated_reads_all': '全部估计片段数', 'estimated_reads_retained': '保留基因的估计片段数',
                 'tpm_sum_all': '全表TPM总和', 'genes_tpm_ge1': 'TPM大于等于1的基因数'}
    pd.DataFrame(qc).rename(columns=qc_labels).to_csv(TABLES / 'discovery_qc.csv', index=False, encoding='utf-8-sig')

    for name, frame in matrices.items():
        save_gene_table(frame, MATRICES / f'FB_{name}.tsv.gz')

    return meta, matrices, len(excluded)


# ---------------------------------------------------------------------------
# 4. TPM 转整数计数（tximport 的 lengthScaledTPM）
# ---------------------------------------------------------------------------

def length_scaled_counts(tpm, effective_length, estimated_reads):
    """按 tximport 的 lengthScaledTPM 算法生成计数尺度输入。

    记基因为 i、文库为 j。先用本次纳入文库的平均有效长度 L̄_i 把 TPM 换回
    长度校正后的丰度 U_ij = TPM_ij × L̄_i，再按每列“估计片段总数 / U 的列和”
    缩放，使每列总量与保留基因行的原始 NumReads 总量一致，最后四舍五入。

    这样得到的整数计数既保留了 TPM 的相对关系，又具备负二项模型所需的
    文库总量尺度。
    """
    # 每个基因用所纳入文库的平均有效长度，再恢复每列原估计片段总量
    unscaled = tpm.mul(effective_length.mean(axis=1), axis=0)
    counts = unscaled.mul(estimated_reads.sum(axis=0) / unscaled.sum(axis=0), axis=1)

    # 取整前做一次守恒检查：各列缩放后的总量应等于原始估计片段总量
    if not np.allclose(counts.sum(), estimated_reads.sum(), rtol=1e-10):
        raise ValueError('Length scaling changed library totals')

    return counts.round().astype('int64')


# ---------------------------------------------------------------------------
# 5. 样本纳入决定
# ---------------------------------------------------------------------------

def write_inclusion(meta):
    """写出主分析与敏感性分析各自的文库纳入表。

    纳入规则来自论文与作者脚本的对照（见 protocol.md）：

    - D3 的两个 FB 异常重复（r2、r3）不进入任何分析；作者脚本保留 r1 与 r4。
    - D28r4 在公开材料中取舍不一致，主分析保留，并在敏感性分析中排除。

    这里用 ``assert`` 式的编号校验锁死对应的 GSM，防止样本表更新后规则悄悄
    错位。
    """
    # 论文记录 D3 排除两个 FB 异常重复，作者脚本保留 D3r1 与 D3r4
    d3_excluded = (meta.day == 3) & meta.replicate.isin([2, 3])
    d28_uncertain = (meta.day == 28) & (meta.replicate == 4)

    if set(meta.index[d3_excluded]) != {'GSM7992335', 'GSM7992336'} or set(meta.index[d28_uncertain]) != {'GSM7992349'}:
        raise ValueError('样本编号与事先定义的质控规则不符。')

    # 每个文库的排除依据，按规则优先级给出可读说明
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


# ---------------------------------------------------------------------------
# 6. 入口
# ---------------------------------------------------------------------------

def main():
    # 只创建本模块会写入的两类目录，图片和报告留给后续模块
    for directory in (MATRICES, TABLES):
        directory.mkdir(parents=True, exist_ok=True)

    annotation()
    meta, _, _ = load_fb()
    write_inclusion(meta)

    print('已整理 24 个提交文库。主分析纳入 22 个，D28 排除敏感性分析纳入 21 个。')


if __name__ == '__main__':
    main()
