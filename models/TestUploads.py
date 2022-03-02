import sys
import unittest
import traceback

import numpy as np


class TestGridUpload(unittest.TestCase):
    def __init__(self, gdf):
        super().__init__()
        self.gdf = gdf

    def test_columns(self):
        self.assertIn('clus_gmm', self.gdf.columns)

    def test_cells(self):
        self.assertIsNotNone(self.gdf['clus_gmm'].unique()[0])
        cluster_nulls = self.gdf['clus_gmm'].isna().sum()
        self.assertEqual(cluster_nulls, 0, f"{cluster_nulls} null values found in 'clus_gmm' column")

    def test_all(self):
        self.test_columns()
        self.test_cells()
