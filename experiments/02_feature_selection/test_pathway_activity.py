from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr

from pathway_activity import spearman_permutation


class PathwayActivityTest(unittest.TestCase):
    def test_rank_correlation_matches_scipy_with_ties(self):
        genes = pd.DataFrame({'a': [0, 0, 2, 3, 5, 5], 'b': [3, 2, 0, 1, 4, 5]})
        scores = pd.DataFrame({'p': [6, 5, 2, 2, 1, 0], 'q': [0, 1, 3, 1, 4, 5]})
        rho, p = spearman_permutation(genes, scores, permutations=199)
        expected = [[spearmanr(genes[g], scores[s]).statistic for s in scores] for g in genes]
        np.testing.assert_allclose(rho, expected)
        self.assertTrue(((0 < p) & (p <= 1)).all())

    def test_real_data_writes_scores_correlations_and_heatmap(self):
        script = Path(__file__).with_name('pathway_activity.py')
        with TemporaryDirectory() as folder:
            subprocess.run([sys.executable, str(script), '--output-dir', folder], check=True, capture_output=True)
            output = Path(folder)
            scores = pd.read_csv(output / 'tables/ssgsea_scores.tsv', sep='\t', index_col=0)
            matrix = pd.read_csv(output / 'tables/ssgsea_spearman_rho.tsv', sep='\t', index_col=0)
            results = pd.read_csv(output / 'tables/ssgsea_correlations.tsv', sep='\t')
            with Image.open(output / 'figures/04_ssgsea_correlation.png') as figure:
                figure.verify()
        self.assertEqual(scores.shape, (22, 8))  # 七个得分加采样日。
        self.assertEqual(matrix.shape, (8, 7))
        self.assertEqual(len(results), 56)
        self.assertTrue(np.isfinite(matrix.to_numpy()).all())
        self.assertTrue(results['q_BH'].ge(results['P']).all())
        self.assertEqual(sum(results['Gene_in_pathway']), 2)


if __name__ == '__main__':
    unittest.main()
