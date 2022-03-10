from unittest import TestCase
import os.path
from Parameters import *


class TestParameters(TestCase):

    def test_test_file_exist(self):
        self.assertTrue(os.path.isdir(Parameters().get_data_dir()))
        self.assertTrue(os.path.isdir(Parameters().get_result_dir()))
        self.assertTrue(os.path.isdir(Parameters().get_pickle_dir()))
        self.assertTrue(os.path.isdir(Parameters().get_exe_dir()))
        self.assertTrue(os.path.isdir(Parameters().get_plot_dir()))

        self.assertTrue(os.path.isfile(Parameters().get_allotype_file()))
        self.assertTrue(os.path.isfile(Parameters().get_htide_immuno_info_file()))
        self.assertTrue(os.path.isfile(Parameters().get_parkhurst_info_file()))
        self.assertTrue(os.path.isfile(Parameters().get_gartner_info_files()[0]))
        self.assertTrue(os.path.isfile(Parameters().get_gartner_info_files()[1]))
        self.assertTrue(os.path.isfile(Parameters().get_tesla_info_files()[0]))
        self.assertTrue(os.path.isfile(Parameters().get_tesla_info_files()[1]))
        self.assertTrue(os.path.isfile(Parameters().get_human_allele_score_dist_file()))
        self.assertTrue(os.path.isfile(Parameters().get_virus_allele_score_dist_file()))
        self.assertTrue(os.path.isfile(Parameters().get_allotype_file()))



