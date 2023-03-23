import os
from unittest import TestCase
from Utils.Util_fct import *


class TestUtilFct(TestCase):

    def test_find_exe(self):

        exe = find_exe(path='/home/localadmin/Programmes/netchop-3.1/Linux_x86_64', exe='netChop')

        self.assertTrue(os.path.isfile(exe))
        self.assertTrue(is_binary_file(exe))

    def test_find_exe(self):
        parameters = Parameters()
        self.assertTrue(os.path.isfile(find_exe(parameters.get_exe_dir(), 'netMHCpan')))

    def test_get_processed_file(self):
        parameters = Parameters()
        file = get_processed_file(parameters.get_result_dir(), "14MH", 'rt_netmhc', 'long')
        self.assertTrue(os.path.isfile(file))

    def test_get_patients(self):
        patients = get_valid_patients('TESLA/HiTIDE', 'short')
        self.assertEqual(19, len(patients))

    def test_get_patients2(self):
        patients = get_valid_patients(['TESLA1', 'TESLA2', 'TESLA3'], 'short')
        self.assertEqual(3, len(patients))

    def test_get_patients3(self):
        patients = get_valid_patients(['HiTIDE', 'TESLA2', 'TESLA3'], 'short')
        patients_htide = get_valid_patients('HiTIDE', 'short')

        self.assertEqual(len(patients_htide)+2, len(patients))

    def test_get_patients4(self):
        patients = get_valid_patients(['Rosenberg'], 'long')
        self.assertEqual(112, len(patients))

        patients = get_valid_patients(['Rosenberg'], 'short')
        self.assertEqual(80, len(patients))

    def test_get_tested_patients(self):
        patients_test = get_valid_patients(dataset='HiTIDE', peptide_type='long')

        mgr = DataManager()
        patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients('long')))
        self.assertEqual(10, len(patients_test))

    def test_get_normalizer(self):
        self.assertTrue(type(get_normalizer('q')) is QuantileTransformer)
        self.assertTrue(type(get_normalizer('p')) is PowerTransformer)
        self.assertTrue(type(get_normalizer('z')) is StandardScaler)
        self.assertTrue(type(get_normalizer('l')) is FunctionTransformer)
        self.assertTrue(get_normalizer('n') is None)
        self.assertTrue(type(get_normalizer("{'a':'n', 'b':'l'}")) is dict)
        d = get_normalizer("{'a':'q', 'b':'l', 'c':'n'}")
        self.assertTrue(type(d['a']) is QuantileTransformer)
        self.assertTrue(type(d['b']) is FunctionTransformer)
        self.assertTrue(d['c'] is None)
