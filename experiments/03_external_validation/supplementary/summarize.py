"""从发现队列、外部 Bulk 和单细胞结果生成冻结候选交接表。

读取阶段一 VST 和 Wald、验证统计与单细胞定位表。
写出 cross_cohort_validation.tsv 和验证 manifest。
"""

from hashlib import sha256
import json

import pandas as pd
from sklearn.metrics import roc_auc_score

from ..config import FROZEN, FROZEN_LF_SHA256, GENES, STAGE_ONE, WELCH_GATE_ALPHA
from ..tables import read, save


def run(evaluation, single_cell, output):
    """校验冻结表内容和顺序后汇总各队列证据，并写出交接文件。"""
    source = FROZEN.read_bytes().replace(b'\r\n', b'\n')
    digest = sha256(source).hexdigest()
    assert digest == FROZEN_LF_SHA256
    frozen = read(FROZEN).set_index('基因符号')
    assert frozen.index.tolist() == GENES
    stage = STAGE_ONE.OUT
    wald = read(stage / 'tables/timecourse_wald.tsv.gz')
    main_d7 = wald[wald['比较'].eq('D7_vs_D0')].set_index('基因ID')
    vst = read(stage / 'matrices/timecourse_vst.tsv.gz', index_col=0)
    included = pd.read_csv(stage / 'tables/sample_inclusion.csv').set_index('GSM编号').loc[vst.columns]
    pair = included[included['时间点（天）'].isin([0, 7])]['时间点（天）']
    tests = read(evaluation / 'tables/bulk_expression_tests.tsv')
    primary = tests[tests.Primary].set_index(['Cohort', 'Gene'])
    roc = read(evaluation / 'tables/bulk_roc_metrics.tsv').set_index(['Cohort', 'Score'])
    loc = read(single_cell / 'tables/sc_expression_by_group.tsv')
    deltas = read(single_cell / 'tables/sc_time_deltas.tsv')
    rows = []
    for gene in GENES:
        gid = frozen.loc[gene, '基因ID']
        d7 = main_d7.loc[gid]
        row = {
            '基因': gene,
            'EnsemblID': gid,
            '发现集_D7方向': '上调' if d7['log2表达倍数变化'] > 0 else '下调',
            '发现集_D7_log2FC': d7['log2表达倍数变化'],
            '发现集_D7_P': d7['P值'],
            '发现集_D7_q': d7['校正后P值'],
            '发现集_D7_AUROC_同队列描述': roc_auc_score(pair.eq(7), vst.loc[gid, pair.index]),
        }
        passed = []
        for cohort in ('GSE111059', 'GSE236374'):
            t = primary.loc[cohort, gene]
            a = roc.loc[cohort, gene]
            prefix = f'{cohort}_D7_'
            row.update(
                {
                    prefix + key: value
                    for key, value in [
                        ('方向', {'up': '上调', 'down': '下调', 'flat': '持平'}[t.Direction]),
                        ('log2尺度均值差', t.Log2_scale_difference),
                        ('P', t.P),
                        ('q', t.BH_q),
                        ('AUROC', a.AUROC),
                        ('AUROC_CI下限', a.CI_low),
                        ('AUROC_CI上限', a.CI_high),
                    ]
                }
            )
            passed.append(bool(t.Pass))
        row['双队列晋级_原始P'] = '是' if all(passed) else '否'
        row['双队列同时满足BH_q'] = (
            '是'
            if all(primary.loc[(c, gene), 'BH_q'] < WELCH_GATE_ALPHA for c in ('GSE111059', 'GSE236374'))
            and all(passed)
            else '否'
        )
        for time in ('Sham', 'D3', 'D7'):
            subset = loc[loc.Gene.eq(gene) & loc.Time.eq(time)]
            if len(subset) and subset.Pct_expressing.max() > 0:
                top = subset.loc[subset.Pct_expressing.idxmax()]
                row[f'单细胞_{time}_检出比例最高群'] = top.Cell_type
                row[f'单细胞_{time}_最高检出比例_pct'] = top.Pct_expressing
                row[f'单细胞_{time}_最高检出群细胞数'] = top.Cells
            else:
                row[f'单细胞_{time}_检出比例最高群'] = '未检出'
        for time in ('D3', 'D7'):
            match = deltas[
                deltas.Gene.eq(gene) & deltas.Time.eq(time) & deltas.Cell_type.eq('Macrophage')
            ]
            if len(match) == 1:
                delta = match.iloc[0]
                row[f'巨噬细胞_{time}_检出率变化百分点'] = delta.Delta_detected_pp
                effect = delta.Delta_mean_log1p_CP10K
                row[f'巨噬细胞_{time}_平均表达变化'] = effect
                row[f'巨噬细胞_{time}_描述方向'] = (
                    ('上升' if effect > 0 else '下降' if effect < 0 else '持平')
                    if pd.notna(effect)
                    else 'NA'
                )
        row['单细胞_P'] = row['单细胞_AUROC'] = 'NA（每时点一份公开文库）'
        rows.append(row)
    summary = pd.DataFrame(rows)
    (output / 'tables').mkdir(parents=True, exist_ok=True)
    (output / 'reports').mkdir(parents=True, exist_ok=True)
    save(summary, 'cross_cohort_validation', output, float_format='%.8g')
    passed = summary.loc[summary['双队列晋级_原始P'].eq('是'), '基因'].tolist()
    (output / 'reports/validation_manifest.json').write_text(
        json.dumps(
            dict(
                stage02_candidate_table_sha256_LF=digest,
                frozen_genes=GENES,
                primary_contrasts=['GSE111059: MI_D7 vs uninjured', 'GSE236374: MI_D7 vs sham'],
                gate='two-sided Welch raw P < 0.05 and positive log2-scale mean difference in both',
                promoted=passed,
            ),
            indent=2,
        )
        + '\n',
        encoding='utf-8',
    )
    print(summary[['基因', '双队列晋级_原始P', '双队列同时满足BH_q']].to_string(index=False))
