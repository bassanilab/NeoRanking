from unittest import TestCase
from pandas.errors import DtypeWarning
from DataManager import *


class TestDataManager(TestCase):

    def test_constructor(self):
        mgr = DataManager()
        patients = mgr.get_valid_patients()
        for p in patients:
            if len(mgr.get_classI_allotypes(p)) == 0:
                print("Patient {0}".format(p))
            self.assertTrue(len(mgr.get_classI_allotypes(p)) > 0)

    def test_constructor2(self):
        mgr = DataManager(immunogenity=True)

        data = mgr.get_processed_data('4253', 'rt', 'short')
        if data is not None:
            self.assertTrue(any(data.apply(lambda row: row['response_type'] in ['CD8', 'CD4/CD8'], axis=1)))

    def test_constructor4(self):
        mgr = DataManager(immunogenity=True)

        for pt in ['long', 'short']:
            patients = mgr.get_immunogenic_patients(pt)
            for p in patients:
                data = mgr.get_processed_data(p, "rt", pt)
                print(f"{p},{pt}")
                if data is not None:
                    self.assertTrue(any(data.apply(lambda row: row['response_type'] in ['CD8', 'CD4/CD8'], axis=1)))

    def test_get_original_data(self):
        mgr = DataManager()
        import warnings
        warnings.filterwarnings(action='ignore', category=ResourceWarning)
        warnings.filterwarnings(action='ignore', category=DtypeWarning)

        for p in mgr.get_valid_patients():
            print("Test loading patient {0} excel long peptide data".format(p))
            data = mgr.get_original_data(p, peptide_type='long')
            self.assertTrue(data is not None and data.shape[0] > 0 and data.shape[1] > 10)
            self.assertEqual(data.shape, mgr.get_original_data(p, peptide_type='long').shape)

        for p in mgr.get_valid_patients():
            print("Test  loading patient {0} excel short peptide data".format(p))
            data = mgr.get_original_data(p, peptide_type='short')
            self.assertTrue(data is not None and data.shape[0] > 0 and data.shape[1] > 10)
            self.assertEqual(data.shape, mgr.get_original_data(p, peptide_type='short').shape)

    # def test_get_processed_data(self):
    #     mgr = DataManager()
    #     file_tag = "netMHC_res"
    #
    #     for p in mgr.get_patients():
    #         data = mgr.get_processed_data(p, file_tag)

    def test_get_classI_allotypes(self):
        mgr = DataManager()

        allotypes = mgr.get_classI_allotypes("TESLA1")
        self.assertEqual(6, len(allotypes))
        self.assertIn('A*02:01', allotypes)
        self.assertIn('A*68:01', allotypes)
        self.assertIn('B*15:07', allotypes)
        self.assertIn('B*44:02', allotypes)
        self.assertIn('C*03:03', allotypes)
        self.assertIn('C*07:04', allotypes)

        allotypes = mgr.get_classI_allotypes("TESLA8")
        self.assertEqual(6, len(allotypes))

    def test_get_processed_file(self):
        mgr = DataManager()

        file = mgr.get_processed_file("0YM1", "hickadoo")
        self.assertEqual(None, file)

    def test_get_original_file(self):
        mgr = DataManager()

        data = mgr.get_original_data("4148")
        self.assertEqual(None, data)
