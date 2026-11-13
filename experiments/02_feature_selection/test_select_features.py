import gzip
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd
from PIL import Image

from select_features import intersect_candidates, paired_omissions, top_membership
from visualize import CHINESE_FONTS, font_manager


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class SelectionTest(unittest.TestCase):
    def test_intersection_uses_ids_and_retains_direction_reversals(self):
        rows = pd.DataFrame([
            ('D7_vs_D0', 'c', 'shared', 1.0, 0.04, '是'),
            ('D7_vs_D0', 'b', 'shared', 2.0, 0.05, '否'),
            ('D3_vs_D0', 'a', 'other', 1.0, 0.01, '是'),
            ('D3_vs_D0', 'c', 'shared', -2.0, 0.01, '是'),
            ('D3_vs_D0', 'b', 'shared', 2.0, 0.01, '是'),
            ('D7_vs_D0', 'a', 'other', 2.0, 0.01, '是'),
        ], columns=['比较', '基因ID', '基因符号', 'log2表达倍数变化', '校正后P值', '通过差异筛选'])
        result = intersect_candidates(rows)
        self.assertEqual(list(result.index), ['a', 'c'])
        self.assertEqual(list(result['相对D0方向']), ['均上调', 'D3下调_D7上调'])

    def test_each_pair_removes_one_library_per_group(self):
        days = pd.Series([3, 7, 3, 7, 7, 7], index=list('abcdef'))
        pairs = list(paired_omissions(days))
        self.assertEqual(len(set(pairs)), 8)
        for d3, d7 in pairs:
            self.assertEqual((days[d3], days[d7]), (3, 7))
            self.assertEqual(days.drop([d3, d7]).value_counts().to_dict(), {7: 3, 3: 1})

    def test_top_membership_keeps_boundary_ties_and_ignores_zero(self):
        np.testing.assert_array_equal(top_membership(np.array([1, 0, 3, 1]), 2), [True, False, True, True])
        np.testing.assert_array_equal(top_membership(np.zeros(3), 20), [False, False, False])

    def test_cli_matches_stage_one_and_reports_all_perturbations(self):
        source = ROOT / 'outputs/01_preprocessing_differential/tables/timecourse_wald.tsv.gz'
        with gzip.open(source, 'rt', encoding='utf-8-sig') as stream:
            wald = pd.read_csv(stream, sep='\t')
        expected = [set(wald.loc[(wald['比较'] == f'D{day}_vs_D0') & (wald['通过差异筛选'] == '是'), '基因ID'])
                    for day in (3, 7)]
        with TemporaryDirectory() as folder:
            output = Path(folder)
            command = [sys.executable, str(HERE / 'select_features.py'), '--output-dir', folder]
            subprocess.run(command, check=True, capture_output=True, text=True)
            result = pd.read_csv(output / 'tables/feature_ranking.tsv', sep='\t')
            report = json.loads((output / 'reports/selection_summary.json').read_text())
            detail = pd.read_csv(output / 'tables/perturbation_selections.tsv', sep='\t')

        self.assertEqual(len(result), 1141)
        self.assertEqual(set(result['基因ID']), expected[0] & expected[1])
        self.assertEqual(result['相对D0方向'].value_counts().to_dict(),
                         {'均上调': 960, '均下调': 178, 'D3下调_D7上调': 2, 'D3上调_D7下调': 1})
        self.assertEqual(report['重抽样次数'], 8)
        for method in ('均值差', 'LASSO', 'ElasticNet', '随机森林'):
            counts = detail.loc[detail['方法'] == method, '基因ID'].value_counts()
            frequencies = result['基因ID'].map(counts).fillna(0) / 8
            np.testing.assert_allclose(result[f'{method}入选频率'], frequencies)

    def test_visualization_cli_writes_three_readable_figures(self):
        if not set(CHINESE_FONTS).intersection(font.name for font in font_manager.fontManager.ttflist):
            self.skipTest('绘图需要中文字体')
        with TemporaryDirectory() as folder:
            for script in ('select_features.py', 'visualize.py'):
                subprocess.run([sys.executable, str(HERE / script), '--output-dir', folder],
                               check=True, capture_output=True, text=True)
            figures = sorted((Path(folder) / 'figures').glob('*.png'))
            self.assertEqual([path.name for path in figures], [
                '01_candidate_intersection.png', '02_candidate_expression.png', '03_selection_stability.png',
            ])
            for path in figures:
                with Image.open(path) as figure:
                    self.assertGreaterEqual(figure.width, 1400)
                    figure.verify()


if __name__ == '__main__':
    unittest.main()
