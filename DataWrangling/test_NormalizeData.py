import unittest
import numpy as np

from Utils.DataManager import DataManager
from DataWrangling.DataTransformer import DataTransformer


class MyTestCase(unittest.TestCase):
    def test_encode(self):
        data_transformer = DataTransformer(peptide_type='neopep', objective='ml', dataset='NCI',
                                           normalizer=DataTransformer.get_normalizer('ml'))
        data_p = DataManager.load_filter_data(peptide_type='neopep', patient='3775')
        data_p, X_p, y_p = data_transformer.apply(data_p)


if __name__ == '__main__':
    unittest.main()
