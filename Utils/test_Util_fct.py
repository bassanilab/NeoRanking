import os
from unittest import TestCase
from Utils.Util_fct import *


class TestUtilFct(TestCase):

    def test_find_exe(self):

        exe = find_exe(path='/home/localadmin/Programmes/netchop-3.1/Linux_x86_64', exe='netChop')

        self.assertTrue(os.path.isfile(exe))
        self.assertTrue(is_binay_file(exe))

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
