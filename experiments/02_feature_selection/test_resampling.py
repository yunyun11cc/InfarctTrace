import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from PIL import Image

from resampling import bootstrap_draws, cv_splits


class ResamplingTest(unittest.TestCase):
    def setUp(self):
        self.labels = np.array([True] * 4 + [False] * 18)

    def test_repeated_cv_keeps_one_d7_per_validation_fold(self):
        splits = list(cv_splits(self.labels, 50))
        self.assertEqual(len(splits), 200)
        for repeat in range(50):
            validation = []
            for train, valid in splits[4 * repeat:4 * (repeat + 1)]:
                self.assertFalse(set(train) & set(valid))
                self.assertEqual(sum(self.labels[train]), 3)
                self.assertEqual(sum(self.labels[valid]), 1)
                self.assertIn(len(train), (16, 17))
                validation.extend(valid)
            self.assertEqual(sorted(validation), list(range(22)))

    def test_bootstrap_keeps_four_d7_and_eighteen_non_d7_with_replacement(self):
        draws = bootstrap_draws(self.labels, 1000)
        self.assertEqual(draws.shape, (1000, 22))
        self.assertTrue(np.all(self.labels[draws].sum(axis=1) == 4))
        self.assertTrue(any(len(set(draw)) < 22 for draw in draws))
        np.testing.assert_array_equal(draws, bootstrap_draws(self.labels, 1000))

    def test_cli_writes_fold_metrics_bootstrap_probabilities_and_figure(self):
        here = Path(__file__).resolve().parent
        with TemporaryDirectory() as folder:
            subprocess.run(
                [sys.executable, str(here / 'resampling.py'), '--output-dir', folder,
                 '--cv-repeats', '1', '--bootstrap', '2', '--jobs', '1'],
                check=True, capture_output=True, text=True,
            )
            subprocess.run(
                [sys.executable, str(here / 'visualize_resampling.py'), '--output-dir', folder],
                check=True, capture_output=True, text=True,
            )
            output = Path(folder)
            metrics = pd.read_csv(output / 'tables/cv_fold_metrics.tsv', sep='\t')
            overlaps = pd.read_csv(output / 'tables/bootstrap_jaccard.tsv', sep='\t')
            probabilities = pd.read_csv(output / 'tables/bootstrap_selection_probability.tsv', sep='\t')
            with Image.open(output / 'figures/03_selection_stability.png') as figure:
                figure.verify()

        self.assertEqual(metrics.groupby('Method')['Fold'].nunique().to_list(), [4, 4, 4])
        self.assertTrue(metrics['D7_valid'].eq(1).all())
        self.assertEqual(overlaps.groupby('Method')['Draw'].nunique().to_list(), [2, 2, 2, 2])
        self.assertTrue(overlaps['Jaccard'].between(0, 1).all())
        self.assertEqual(probabilities['基因ID'].nunique(), 3907)
        self.assertTrue(probabilities['Mean'].between(0, 1).all())


if __name__ == '__main__':
    unittest.main()
