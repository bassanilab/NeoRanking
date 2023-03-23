from unittest import TestCase
from DataWrangling.DataLoader import *
from sklearn.preprocessing import *

features = ['rnaseq_ref_support', 'Sample_Tissue_expression_GTEx', 'MIN_MUT_RANK_CI_PRIME',
           'CSCAPE_score', 'GTEx_all_tissues_expression_median',
           'bestWTMatchScore_I', 'bestWTPeptideCount_I', 'rnaseq_gene_expression_quartile',
           'mut_Rank_BA_0', 'mut_Rank_EL_0', 'mut_Rank_BA_1', 'mut_Rank_EL_1', 'mut_nr_weak_binding_alleles_0',
           'mut_nr_weak_binders_0', 'mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_netchop_Ct_score_0']


class TestDataLoader(TestCase):
    def test_load_patients(self):

        data_loader = DataLoader(features=features)

        df1, X1, y1 = data_loader.load_patients(['TESLA1', 'TESLA2'], file_tag='rt')
        df2, X2, y2 = data_loader.load_patients(['TESLA3', 'TESLA4'], file_tag='rt')
        df3, X3, y3 = data_loader.load_patients(['TESLA1', 'TESLA2', 'TESLA3', 'TESLA4'], file_tag='rt')

        self.assertEqual(X1.shape[0]+X2.shape[0], X3.shape[0])
        self.assertEqual(df1.shape[0]+df2.shape[0], df3.shape[0])
        self.assertEqual(len(y1)+len(y2), len(y3))

    def test_load_patients2(self):

        data_loader = DataLoader(features=features)

        df, X, y = data_loader.load_patients(['0YM1'], file_tag='rt')
        self.assertTrue(y is not None)
        self.assertEqual(2, len(np.where(y == 1)[0]))

        df, X, y = data_loader.load_patients(['4262'], file_tag='rt')
        self.assertTrue(y is not None)
        self.assertEqual(1, len(np.where(y == 1)[0]))

        df, X, y = data_loader.load_patients('4235', file_tag='rt')
        self.assertTrue(y is None)

    def test_load_patients3(self):

        data_loader = DataLoader(features=features)

        df, X, y = data_loader.load_patients('1913', file_tag='rt', peptide_type='short')
        self.assertTrue(y is None)

        df, X, y = data_loader.load_patients('14MH', file_tag='rt', peptide_type='short')
        self.assertTrue(df.shape[0] > 0)

    def test_load_patients4(self):

        data_loader = DataLoader(features=Parameters().get_features())

        df, X, y = data_loader.load_patients('0YM1', file_tag='rt', peptide_type='short', nr_non_immuno_rows=500)
        self.assertEqual(94, sum(df.apply(lambda row: row['response_type'] == 'negative', axis=1)))
        self.assertEqual(4, sum(df.apply(lambda row: row['response_type'] in ['CD8', 'CD4/CD8'], axis=1)))

        df, X, y = data_loader.load_patients('3703', file_tag='rt', peptide_type='short', nr_non_immuno_rows=500)
        self.assertEqual(500, sum(df.apply(lambda row: row['response_type'] == 'negative', axis=1)))
        self.assertEqual(4, sum(df.apply(lambda row: row['response_type'] in ['CD8', 'CD4/CD8'], axis=1)))

    def test_load_patients5(self):

        max_netmhc_rank = 0.5
        data_loader = DataLoader(features=Parameters().get_features(), max_netmhc_rank=max_netmhc_rank)

        df, X, y = data_loader.load_patients('0YM1', file_tag='rt', peptide_type='short')
        self.assertTrue(all(df.apply(lambda row: row['mutant_rank_netMHCpan'] <= max_netmhc_rank, axis=1)))

        df, X, y = data_loader.load_patients('3703', file_tag='rt', peptide_type='short')
        self.assertTrue(all(df.apply(lambda row: row['mutant_rank_netMHCpan'] <= max_netmhc_rank, axis=1)))

    def test_load_patients6(self):
        mgr = DataManager()
        patients_with_data = mgr.get_valid_patients('short')
        data_loader = DataLoader(features=Parameters().get_features(), max_netmhc_rank=2)
        for p in patients_with_data:
            df, X, y = data_loader.load_patients(p, file_tag='rt', peptide_type='short')
            if df is not None and 'nb_same_mutation_Intogen' in df.columns:
                mv_cnt = df['nb_same_mutation_Intogen'].isna().sum()
                mv_cnt = X['nb_same_mutation_Intogen'].isna().sum()
                print(f"Patient {p} has {mv_cnt} of {df.shape[0]} missing values")
            elif df is None:
                print(f"Patient {p} has empty dataframe")
            else:
                print(f"Patient {p} has no column nb_same_mutation_Intogen")

    def test_load_patients7(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=Parameters().get_features(),
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=2)
        df, X, y = data_loader.load_patients('3309', file_tag='rt', peptide_type='short')
        mv_cnt = X['nb_same_mutation_Intogen'].isna().sum()
        print(f"Patient 3678 has {mv_cnt} of {df.shape[0]} missing values")

    def test_load_patients8(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=Parameters().get_features(),
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=10000)
        df, X, y = data_loader.load_patients('2098', file_tag='rt_netmhc', peptide_type='long')
        mv_cnt = X['nb_same_mutation_Intogen'].isna().sum()
        print(f"Patient 2098 has {mv_cnt} of {df.shape[0]} missing values")

    def test_load_patients9(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=Parameters().get_features(),
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=10000)
        df, X, y = data_loader.load_patients('2098', file_tag='rt', peptide_type='short', nr_non_immuno_rows=25)
        self.assertEqual(sum(y == 0), 25)
        self.assertEqual(df.shape[0] - sum(y == 1), 25)
        self.assertEqual(X.shape[0] - sum(y == 1), 25)

    def test_load_patients10(self):
        excl_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1",
                      "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "TRBV3", "TRBV5", "TRBV6", "TRBV6-1", "TRBV10", "TRBV10-1",
                      "TRBV11", "TRAV12", "KRT1", "PRSS3"]

        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=['bestWTMatchScore_I', 'bestWTMatchOverlap_I', 'bestMutationScore_I'],
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0)
        df, X, y = data_loader.load_patients('1913', file_tag='rt', peptide_type='long')
        self.assertEqual(df.shape[0], 2895)

        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=['bestWTMatchScore_I', 'bestWTMatchOverlap_I', 'bestMutationScore_I'],
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, excluded_genes=excl_genes)
        df, X, y = data_loader.load_patients('1913', file_tag='rt', peptide_type='long')
        self.assertEqual(df.shape[0], 2893)

    def test_load_patients9(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=Parameters().get_features(),
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=10000)
        df, X, y = data_loader.load_patients('2556', file_tag='rt_stab_chop_tap_mbp_tcr_prop', peptide_type='short', nr_non_immuno_rows=25)
        self.assertTrue(all([mbp in [True, False, np.nan] for mbp in df.mut_is_binding_pos.unique()]))
        self.assertTrue(all([mbp in [True, False, np.nan] for mbp in X.mut_is_binding_pos.unique()]))

    def test_load_patients10(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=Parameters().get_features(),
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=10000)
        df, X, y = data_loader.load_patients('1913', file_tag='rt_stab_chop_tap_mbp_tcr_prop', peptide_type='short')
        self.assertTrue(all([mbp in [True, False, np.nan] for mbp in df.mut_is_binding_pos.unique()]))
        self.assertTrue(all([mbp in [True, False, np.nan] for mbp in X.mut_is_binding_pos.unique()]))

    def test_load_patients11(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=features,
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=10000)
        df, X, y = data_loader.load_patients('14MH', file_tag='rt_stab_chop_tap_mbp_tcr_prop', peptide_type='short')
        self.assertEqual(1, sum(df["response_type"] == "CD8"))

        df, X, y = data_loader.load_patients('14MH', file_tag='rt_stab_chop_tap_mbp_tcr_prop', peptide_type='long')
        self.assertEqual(2, sum(df["response_type"] == "CD8"))

    def test_load_patients12(self):
        data_loader = DataLoader(transformer=DataTransformer(), normalizer=QuantileTransformer(),
                                 features=features,
                                 mutation_types=['SNV'], response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                 immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0, max_netmhc_rank=10000)
        df, X, y = data_loader.load_patients('1HU3', file_tag='rt_stab_chop_tap_mbp_tcr_prop', peptide_type='short')
        self.assertEqual(0, sum(df["response_type"] == "CD8"))

        df, X, y = data_loader.load_patients('1HU3', file_tag='rt_stab_chop_tap_mbp_tcr_prop', peptide_type='long')
        self.assertEqual(2, sum(df["response_type"] == "CD8"))
