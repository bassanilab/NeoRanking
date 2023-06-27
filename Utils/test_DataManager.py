from unittest import TestCase
import time
import pandas as pd

from Utils.DataManager import DataManager
from Utils.GlobalParameters import GlobalParameters
from DataWrangling.DataTransformer import DataTransformer
from Utils.Util_fct import *


class TestDataManager(TestCase):

    def test_get_original_data1(self):
        peptide_type = 'neopep'
        patient = 'Patient1'
        start = time.time()
        data = DataManager.load_filter_data(peptide_type, patient=patient)
        print("{0} data loaded for patient {1} in {2:.3f} secs".format(peptide_type, patient, time.time()-start))

        start = time.time()
        data = DataManager.load_filter_data(peptide_type, patient=patient)
        print("{0} data loaded for patient {1} in {2:.3f} secs".format(peptide_type, patient, time.time()-start))

        dataset = 'TESLA'
        start = time.time()
        data = DataManager.load_filter_data(peptide_type, dataset=dataset)
        print("{0} data loaded for dataset {1} in {2:.3f} secs".format(peptide_type, dataset, time.time()-start))

    def test_has_immunogenic_peptides(self):
        peptide_type = 'neopep'
        self.assertTrue(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='Patient1'))
        self.assertTrue(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='TESLA3'))
        self.assertTrue(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='TESLA3'))
        self.assertFalse(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='4217'))
        self.assertFalse(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='4232'))

        peptide_type = 'mutation'
        self.assertTrue(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='Patient1'))
        self.assertTrue(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='TESLA3'))
        self.assertTrue(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='TESLA3'))
        self.assertFalse(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='4217'))
        self.assertFalse(DataManager.has_immunogenic_peptides(peptide_type=peptide_type, patient='4232'))

    def test_get_processed_data(self):
        data, X, y = DataManager.filter_processed_data(peptide_type='mutation', objective='ml')
        self.assertTrue(all(data.loc[y == 1, 'response_type'] == 'CD8'))
        self.assertEqual(len(y), X.shape[0])
        self.assertEqual(len(y), data.shape[0])

        data, X, y = DataManager.filter_processed_data(peptide_type='neopep', objective='ml')
        self.assertTrue(all(data.loc[y == 1, 'response_type'] == 'CD8'))
        self.assertEqual(len(y), X.shape[0])
        self.assertEqual(len(y), data.shape[0])
        self.assertEqual(GlobalParameters.nr_non_immuno_neopeps, sum(data['response_type'] != 'CD8'))

    def test_get_processed_data2(self):
        data, X, y = DataManager.filter_processed_data(peptide_type='mutation', objective='ml', dataset='NCI',
                                                       response_types=['CD8', 'negative'])
        self.assertTrue(all(data.loc[y == 1, 'response_type'] == 'CD8'))
        self.assertTrue(all(data.dataset == 'NCI'))
        self.assertTrue(all((data.response_type == 'CD8') | (data.response_type == 'negative')))
        self.assertEqual(len(y), X.shape[0])
        self.assertEqual(len(y), data.shape[0])

    def test_peptide_count(self):
        peptide_type = 'mutation'
        data, X, y = DataManager.filter_processed_data(peptide_type, objective='plot', sample=False)
        print(data.bestWTPeptideCount_I.describe())
        print(X.bestWTPeptideCount_I.describe())

    def test_read(self):
        peptide_type = 'mutation'
        objective = 'ml'
        ml_sel_data_file_name, norm_data_file_name = DataManager.get_processed_data_files(peptide_type, objective)
        print(norm_data_file_name)
        X = pd.read_csv(norm_data_file_name, sep="\t", header=0, dtype=get_processed_types(peptide_type, objective),
                        engine='pyarrow', dtype_backend='pyarrow')

        print(X.bestWTPeptideCount_I.describe())

    def test_patient(self):
        data = DataManager.load_original_data(peptide_type='neopep', ml_row_selection=False)
        data_p = DataManager.load_filter_data(peptide_type='neopep', patient='3775')
        data_transformer = DataTransformer('NCI', 'neopep', DataTransformer.get_normalizer('ml'), 'ml')
        data_p, X_p, y_p = data_transformer.apply(data_p)
