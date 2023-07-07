from unittest import TestCase
import os.path
from GlobalParameters import *


class TestGlobalParameters(TestCase):

    def test_files_exist(self):
        self.assertTrue(os.path.isdir(GlobalParameters.base_dir))
        self.assertTrue(os.path.isdir(GlobalParameters.data_dir))
        self.assertTrue(os.path.isdir(GlobalParameters.plot_dir))
        self.assertTrue(os.path.isdir(GlobalParameters.classifier_result_dir))

        self.assertTrue(os.path.isfile(GlobalParameters.neopep_data_file))
        self.assertTrue(os.path.isfile(GlobalParameters.mutation_data_file))
        self.assertTrue(os.path.isfile(GlobalParameters.allo_file))



