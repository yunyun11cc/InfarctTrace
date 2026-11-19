import unittest

import numpy as np
import pandas as pd

from prepare import aggregate_probes, align_samples, validate_counts


class ValidationInputTest(unittest.TestCase):
    def test_align_samples_by_explicit_identifiers(self):
        values = pd.DataFrame({'MI_7d_1': [7], 'MI_sham_1': [2]}, index=['Gene'])
        samples = pd.DataFrame({'Source_column': ['MI_sham_1', 'MI_7d_1']}, index=['GSM1', 'GSM2'])
        expected = pd.DataFrame({'GSM1': [2], 'GSM2': [7]}, index=['Gene'])
        pd.testing.assert_frame_equal(align_samples(values, samples), expected)
        with self.assertRaises(ValueError):
            align_samples(values.rename(columns={'MI_7d_1': 'other'}), samples)
        with self.assertRaises(ValueError):
            align_samples(values, samples.set_axis(['GSM1', 'GSM1']))

    def test_aggregate_probes_excludes_ambiguous_and_control_rows(self):
        values = pd.DataFrame({'GSM1': [2.5, 3.5, 9.5, 8, 12, 10]},
                              index=['p1', 'p2', 'p3', 'p4', 'AFFX-a', 'p5'])
        annotation = pd.DataFrame({
            'probe_id': ['p1', 'p2', 'p3', 'p3', 'p4', 'AFFX-a', 'p5'],
            'gene_id': ['1', '1', '1', '2', '3', '4', None],
            'is_multiple': [0, 0, 0, 0, 1, 0, 0],
        })
        genes, audit = aggregate_probes(values, annotation)
        self.assertEqual(genes.index.tolist(), ['1'])
        self.assertEqual(genes.at['1', 'GSM1'], 3)
        self.assertEqual(audit['Status'].to_dict(), {
            'p1': 'mapped', 'p2': 'mapped', 'p3': 'ambiguous', 'p4': 'ambiguous',
            'AFFX-a': 'control', 'p5': 'unmapped',
        })

    def test_rejects_invalid_counts(self):
        for value in (-1, 2.5, np.nan, np.inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_counts(pd.DataFrame({'GSM1': [value]}, index=['Gene']))
        values = pd.DataFrame({'GSM1': [0, 12.0]}, index=['Gene1', 'Gene2'])
        self.assertEqual(validate_counts(values)['GSM1'].tolist(), [0, 12])
