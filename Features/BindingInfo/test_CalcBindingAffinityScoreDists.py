from unittest import TestCase
from Features.BindingInfo.CalcBindingAffinityScoreDists import *
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *


class TestCalcAllelePropensity(TestCase):
    def test_constructor(self):

        calcAllelePropensity = CalcAllelePropensity('/home/localadmin/Priorization/Plots/allele_score_data_human.txt',
                                                    '0YM1', 'SNV', 'CD8,CD4/CD8', 'CD8,CD4/CD8,negative',
                                                    'rt_netmhc_stab_chop_tap_mbp_tcr')

        self.assertEqual(4, len(calcAllelePropensity.propensities))
        self.assertEqual(6, len(calcAllelePropensity.immunogenicity_ratios))
        self.assertEqual(89, len(calcAllelePropensity.propensities[9]))

        self.assertTrue('HLA-A*01:01' in calcAllelePropensity.propensities[9])
        self.assertTrue('HLA-B*35:01' in calcAllelePropensity.propensities[9])
        self.assertTrue('HLA-B*35:02' in calcAllelePropensity.propensities[9])
        self.assertTrue('HLA-B*35:03' in calcAllelePropensity.propensities[9])

    def test_constructor2(self):

        calcAllelePropensity = CalcAllelePropensity('/home/localadmin/Priorization/Plots/allele_score_data_human.txt',
                                                    '1913', 'SNV', 'CD8,CD4/CD8', 'CD8,CD4/CD8,negative',
                                                    'rt_netmhc_stab_chop_tap_mbp_tcr')

        self.assertEqual(0, len(calcAllelePropensity.propensities))
        self.assertEqual(0, len(calcAllelePropensity.immunogenicity_ratios))

    def test_add_features_long(self):
        mgr = DataManager()
        patients = mgr.get_valid_patients()
        calcAllelePropensity = CalcAllelePropensity('/home/localadmin/Priorization/Plots/allele_score_data_human.txt',
                                                    patients, 'SNV', 'CD8,CD4/CD8', 'CD8,CD4/CD8,negative',
                                                    'rt_netmhc_stab_chop_tap_mbp_tcr')

        data = calcAllelePropensity.add_features('0YM1', 'rt_netmhc_stab_chop_tap_mbp_tcr', 'prop', 'long', write_res=False)

        self.assertTrue('mut_allele_propensity_0' in data.columns)
        self.assertTrue('wt_allele_propensity_0' in data.columns)

    def test_add_features_short(self):
        mgr = DataManager()
        annotator = RosenbergImmunogenicityAnnotatorLong(mgr)
        calcAllelePropensity = CalcAllelePropensity(Parameters().get_human_allele_score_dist_file(),
                                                    annotator.get_patients('parkhurst'), 'SNV', 'CD8,CD4/CD8',
                                                    'CD8,CD4/CD8,negative', 'rt_netmhc', 'long')

        data = calcAllelePropensity.add_features('0YM1', 'rt_stab_chop_tap_mbp_tcr', 'prop', 'short', write_res=False)

        self.assertTrue('mut_allele_propensity' in data.columns)
        self.assertTrue(len(data['mut_allele_propensity'].unique()) > 1)


