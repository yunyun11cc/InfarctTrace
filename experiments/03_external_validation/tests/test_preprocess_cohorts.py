"""用小型表达矩阵验证步骤 01 的样本对齐、探针聚合与计数约束。"""

from importlib import import_module

import unittest

import numpy as np
import pandas as pd

preprocess = import_module('experiments.03_external_validation.step01_preprocess_cohorts')


class ValidationInputTest(unittest.TestCase):
    """回归覆盖步骤 01 的三类输入校验。"""

    def test_align_samples_by_explicit_identifiers(self):
        """确认显式样本 ID 映射和列顺序，拒绝不匹配和重复样本。"""
        values = pd.DataFrame({'MI_7d_1': [7], 'MI_sham_1': [2]}, index=['Gene'])
        samples = pd.DataFrame({'Source_column': ['MI_sham_1', 'MI_7d_1']}, index=['GSM1', 'GSM2'])
        expected = pd.DataFrame({'GSM1': [2], 'GSM2': [7]}, index=['Gene'])
        pd.testing.assert_frame_equal(preprocess.align_samples(values, samples), expected)
        with self.assertRaises(ValueError):
            preprocess.align_samples(values.rename(columns={'MI_7d_1': 'other'}), samples)
        with self.assertRaises(ValueError):
            preprocess.align_samples(values, samples.set_axis(['GSM1', 'GSM1']))

    def test_aggregate_probes_excludes_ambiguous_and_control_rows(self):
        """确认只聚合唯一映射探针并报告排除状态。"""
        values = pd.DataFrame(
            {'GSM1': [2.5, 3.5, 9.5, 8, 12, 10]}, index=['p1', 'p2', 'p3', 'p4', 'AFFX-a', 'p5']
        )
        annotation = pd.DataFrame(
            {
                'probe_id': ['p1', 'p2', 'p3', 'p3', 'p4', 'AFFX-a', 'p5'],
                'gene_id': ['1', '1', '1', '2', '3', '4', None],
                'is_multiple': [0, 0, 0, 0, 1, 0, 0],
            }
        )
        genes, audit = preprocess.aggregate_probes(values, annotation)
        self.assertEqual(genes.index.tolist(), ['1'])
        self.assertEqual(genes.at['1', 'GSM1'], 3)
        self.assertEqual(
            audit['Status'].to_dict(),
            {
                'p1': 'mapped',
                'p2': 'mapped',
                'p3': 'ambiguous',
                'p4': 'ambiguous',
                'AFFX-a': 'control',
                'p5': 'unmapped',
            },
        )

    def test_rejects_invalid_counts(self):
        """确认负值、小数、缺失和无穷值会被拒绝。"""
        for value in (-1, 2.5, np.nan, np.inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                preprocess.validate_counts(pd.DataFrame({'GSM1': [value]}, index=['Gene']))
        values = pd.DataFrame({'GSM1': [0, 12.0]}, index=['Gene1', 'Gene2'])
        self.assertEqual(preprocess.validate_counts(values)['GSM1'].tolist(), [0, 12])
