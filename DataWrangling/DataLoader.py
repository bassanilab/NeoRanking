import random
import numpy as np
import pandas as pd
import warnings

from pandas.core.common import SettingWithCopyWarning
from pandas.core.dtypes.concat import union_categoricals
from pandas.errors import DtypeWarning
from scipy.stats import rankdata

from Utils.DataManager import DataManager
from DataWrangling.Transform_Data import DataTransformer
from Utils.Parameters import Parameters
from Utils.Util_fct import *


class DataLoader:

    def __init__(self, features=[], transformer=None, normalizer=None, mutation_types=['SNV'],
                 response_types=['CD8', 'CD4/CD8', 'negative'], immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=1,
                 cat_type='float', max_netmhc_rank=-1, cat_encoders=None, excluded_genes=None):
        self.transformer = transformer
        self.normalizer = normalizer
        self.features = features
        self.data_manager = DataManager()
        self.mutation_types = mutation_types
        self.response_types = response_types
        self.immunogenic = immunogenic
        self.min_nr_immuno = min_nr_immuno
        self.cat_type = cat_type
        self.cat_dims = None
        self.cat_idx = None
        self.cat_encoders = cat_encoders
        self.max_netmhc_rank = max_netmhc_rank
        self.excluded_genes = excluded_genes

        return

    def load_patients(self, patients, file_tag, peptide_type='long', verbose=True, required_columns=None,
                      nr_non_immuno_rows=-1):

        warnings.filterwarnings(action='ignore', category=DtypeWarning)
        warnings.filterwarnings(action='ignore', category=UserWarning)
        warnings.filterwarnings(action='ignore', category=SettingWithCopyWarning)

        if isinstance(patients, str):
            patients = [patients]

        combined_df = None
        combined_X = None
        combined_y = None
        for i, p in enumerate(patients):
            if verbose:
                print("Loading patient {0}/{1} ({2} of {3})...".format(p, peptide_type, i+1, len(patients)))
            if combined_df is None:
                df, X, y = self.load_patient(p, file_tag, peptide_type, required_columns)
                if df is not None and np.sum(y == 1) >= self.min_nr_immuno:
                    combined_df = df
                    combined_X = X
                    combined_y = y
            else:
                df, X, y = self.load_patient(p, file_tag, peptide_type)
                if df is not None and np.sum(y == 1) >= self.min_nr_immuno:
                    combined_df, df = self.combine_categories(combined_df, df)
                    combined_df = combined_df.append(df, ignore_index=True)
                    combined_X, X = self.combine_categories(combined_X, X)
                    combined_X = combined_X.append(X, ignore_index=True)
                    combined_y = np.append(combined_y, y)

        if nr_non_immuno_rows > 0:
            combined_df, combined_X, combined_y = \
                self.sample_rows(combined_df, combined_X, combined_y, nr_non_immuno_rows)

        return combined_df, combined_X, combined_y

    def combine_categories(self, df1, df2):
        if self.cat_type not in ['float', 'int']:
            cat_features = \
                [c for c in df1.columns if c in Parameters().get_categorical_features() and c in self.cat_encoders]
            for c in cat_features:
                uc = union_categoricals([df1[c], df2[c]])
                df1.loc[:, c] = pd.Categorical(df1[c], categories=uc.categories)
                df2.loc[:, c] = pd.Categorical(df2[c], categories=uc.categories)

        return df1, df2

    @staticmethod
    def sample_rows(data, X, y, nr_non_immuno_rows):
        if sum(y == 0) < nr_non_immuno_rows:
            return data, X, y

        idx = random.sample(range(sum(y == 0)), nr_non_immuno_rows)
        X_1 = X.loc[y == 1, :]
        X_0 = X.loc[y == 0, :]
        if X_0.shape[0] > nr_non_immuno_rows:
            X_0 = X_0.iloc[idx, :]
        X_s = pd.concat([X_1, X_0])

        X_1 = data.loc[y == 1, :]
        X_0 = data.loc[y == 0, :]
        if X_0.shape[0] > nr_non_immuno_rows:
            X_0 = X_0.iloc[idx, :]
        data_s = pd.concat([X_1, X_0])

        y_0 = y[y == 0]
        y_1 = y[y == 1]
        y_0 = y_0[idx]
        y_s = np.append(y_1, y_0)

        return data_s, X_s, y_s

    def load_patient(self, patient, file_tag, peptide_type='long', required_columns=None):
        if peptide_type == 'long':
            return self.load_patient_long(patient, file_tag, required_columns)
        elif peptide_type == 'short':
            return self.load_patient_short(patient, file_tag, required_columns)
        else:
            return None, None, None

    def load_patient_long(self, patient, file_tag, required_columns=None):
        df = self.data_manager.get_processed_data(patient, file_tag, 'long', required_columns)

        if df is None:
            return None, None, None

        df = self.filter_rows_long(df)

        if df.shape[0] == 0:
            return None, None, None

        if 'response_type' in df.columns:
            y = df.apply(lambda row: int(row['response_type'] in self.immunogenic), axis=1)
#            y = df.apply(lambda row: 2*int(row['response_type'] in self.immunogenic)-1, axis=1)
            y = np.array(y, dtype=int)
        else:
            y = None
            print("No response_type column in dataframe of patient {0}. Skip. ".format(patient))

        df = self.process_columns_long(df, patient)

        if self.transformer is not None and not df.empty:
            df = self.transformer.fill_missing_values(df)

        if self.features is not None and len(self.features) > 0:
            features_sel = [f for f in df.columns if f in self.features]
            df_sel = df.loc[:, features_sel]
        else:
            df_sel = df

        if self.normalizer is not None and not df.empty:
            X = DataTransformer().normalize(df_sel, self.normalizer)
        else:
            X = df_sel.copy()

        if y is not None:
            if len(y) > 0:
                df.insert(0, "response", y)
            else:
                df['response'] = []

        if 'Nb_Samples' in df.columns:
            df.loc[:, 'Nb_Samples'] = df.loc[:, 'Nb_Samples']/df.loc[:, 'Nb_Samples'].max()

        X = self.encode_cat_features(X)

        return df, X, y

    def filter_rows_long(self, df):
        df = df.loc[df.apply(lambda row: row['mutation_type'] in self.mutation_types, axis=1)]

        if df.shape[0] > 0 and 'response_type' in df.columns:
            df = df.loc[df.apply(lambda row: row['response_type'] in self.response_types, axis=1)]

        # Only load data with at least one valid netmhc prediction
        if df.shape[0] > 0 and 'mut_peptide_pos_0' in df.columns:
            df = df.loc[df['mut_peptide_pos_0'].notna()]

        # Only load data with RNA-seq
        if df.shape[0] > 0 and 'rnaseq_gene_expression_quartile' in df.columns:
            df = df.loc[df['rnaseq_gene_expression_quartile'].notna()]

        if df.shape[0] > 0 and self.max_netmhc_rank > 0:
            df = df.loc[df.apply(lambda row: row['mut_Rank_EL_0'] < self.max_netmhc_rank, axis=1)]

        if df.shape[0] > 0 and self.excluded_genes is not None and len(self.excluded_genes) > 0:
            df = df.loc[df.apply(lambda row: row['gene'] not in self.excluded_genes, axis=1)]

        df.reset_index(inplace=True, drop=True)

        return df

    def load_patient_short(self, patient, file_tag, required_columns=None):
        df = self.data_manager.get_processed_data(patient, file_tag, 'short', required_columns)

        if df is None:
            return None, None, None

        df = self.filter_rows_short(df)

        if df.shape[0] == 0:
            return None, None, None

        if 'response_type' in df.columns:
            y = df.apply(lambda row: int(row['response_type'] in self.immunogenic), axis=1)
#            y = df.apply(lambda row: 2*int(row['response_type'] in self.immunogenic)-1, axis=1)
            y = np.array(y, dtype=int)
        else:
            y = None
            print("No response_type column in dataframe of patient {0}. Skip. ".format(patient))

        df = self.process_columns_short(df, patient)

        if self.transformer is not None and not df.empty:
            df = self.transformer.fill_missing_values(df)

        if self.features is not None and len(self.features) > 0:
            features_sel = [f for f in df.columns if f in self.features]
            df_sel = df.loc[:, features_sel]
        else:
            df_sel = df

        if self.normalizer is not None and not df.empty:
            X = DataTransformer().normalize(df_sel, self.normalizer)
        else:
            X = df_sel.copy()

        if y is not None:
            if len(y) > 0:
                df.insert(0, "response", y)
            else:
                df.loc[:, 'response'] = []

        X = self.encode_cat_features(X)

        return df, X, y

    def filter_rows_short(self, df):
        df = df.loc[df.apply(lambda row: row['mutation_type'] in self.mutation_types, axis=1)]

        if df.shape[0] > 0 and 'mut_is_binding_pos' in df.columns:
            df = df.loc[df['mut_is_binding_pos'].notna()]

        if df.shape[0] > 0 and 'response_type' in df.columns:
            df = df.loc[df.apply(lambda row: row['response_type'] in self.response_types, axis=1)]

        if df.shape[0] > 0 and 'Final_Peptide_Class' in df.columns:
            df = df.loc[df.apply(lambda row: row['Final_Peptide_Class'] in ['HLA_I'], axis=1)]

        # Only load data with at least one valid netmhc prediction
        if df.shape[0] > 0 and 'mutant_rank_netMHCpan' in df.columns:
            df = df.loc[df['mutant_rank_netMHCpan'].notna()]

        # Only load data with RNA-seq
        if df.shape[0] > 0 and 'rnaseq_gene_expression_quartile' in df.columns:
            df = df.loc[df['rnaseq_gene_expression_quartile'].notna()]

        # Only load data with at least one valid netmhc prediction
        if df.shape[0] > 0:
            df = df.loc[df.apply(lambda row: row['pep_mut_start'] <= len(row['mutant_seq']), axis=1)]

        if df.shape[0] > 0 and self.max_netmhc_rank > 0:
            df = df.loc[df.apply(lambda row: row['mutant_rank_netMHCpan'] < self.max_netmhc_rank, axis=1)]

        if df.shape[0] > 0 and self.excluded_genes is not None and len(self.excluded_genes) > 0:
            df = df.loc[df.apply(lambda row: row['gene'] not in self.excluded_genes, axis=1)]

        df.reset_index(inplace=True, drop=True)

        return df

    def encode_cat_features(self, df):
        if self.cat_encoders:
            cat_features = \
                [c for c in df.columns if c in Parameters().get_categorical_features() and c in self.cat_encoders]
            self.cat_dims = {}
            for f in cat_features:
                encoder = self.cat_encoders[f]
                self.cat_dims[f] = encoder.get_nr_classes()
                if self.cat_type == 'float':
                    df.loc[:, f] = encoder.transform(df[f].values, 'float')
                    df = df.astype({f: float})
                elif self.cat_type == 'int':
                    df.loc[:, f] = encoder.transform(df[f].values, 'int')
                    df = df.astype({f: int})
                else:
                    df = df.astype({f: 'category'})
        else:
            cat_features = [c for c in df.columns if c in Parameters().get_categorical_features()]
            df = df.astype(dict.fromkeys(cat_features, 'category'))
            self.cat_dims = None

        return df

    def process_columns_long(self, df, patient):

        if '%rnaseq_alt_support' in df.columns:
            df = df.rename(columns={"%rnaseq_alt_support": "rnaseq_alt_support",
                                    "%rnaseq_ref_support": "rnaseq_ref_support"})

        allowed_cols = Parameters().get_features()
        keep_cols = [c for c in df.columns if c in allowed_cols]
        df = df[keep_cols]

        df.loc[:, 'patient'] = np.full(df.shape[0], patient)
        df.loc[:, 'dataset'] = np.full(df.shape[0], get_patient_group(patient))
        df.loc[:, 'train_test'] = np.full(df.shape[0], get_ml_group(patient, 'long'))

        netMHCpan_ranks = [int(c[c.rfind('_')+1:]) for c in df.columns if 'mut_peptide_pos_' in c]
        if len(netMHCpan_ranks) > 1:
            next_best_EL_mut_ranks = np.zeros(df.shape[0])
            next_best_BA_mut_ranks = np.zeros(df.shape[0])
            for i in df.index:
                allele = df.loc[i, 'mut_allele_0']
                mut_EL_rank = 0.0
                mut_BA_rank = 0.0
                for j in netMHCpan_ranks[1:]:
                    if df.loc[i, 'mut_allele_'+str(j)] != allele:
                        mut_EL_rank = df.loc[i, 'mut_Rank_EL_'+str(j)]
                        mut_BA_rank = df.loc[i, 'mut_Rank_BA_'+str(j)]
                        break

                next_best_EL_mut_ranks[i] = mut_EL_rank
                next_best_BA_mut_ranks[i] = mut_BA_rank

            df['next_best_EL_mut_ranks'] = next_best_EL_mut_ranks
            df['next_best_BA_mut_ranks'] = next_best_BA_mut_ranks

            for j in netMHCpan_ranks[0:]:

                if {'mut_Rank_EL_'+str(j), 'wt_Rank_EL_'+str(j)}.issubset(df.columns):
                    df['DAI_'+str(j)] = \
                        df.apply(lambda row: self.calc_dai(row['mut_Rank_EL_'+str(j)], row['wt_Rank_EL_'+str(j)]), axis=1)
                else:
                    df['DAI_'+str(j)] = 0

        return df

    def process_columns_short(self, df, patient):

        if '%rnaseq_alt_support' in df.columns:
            df = df.rename(columns={"%rnaseq_alt_support": "rnaseq_alt_support",
                                    "%rnaseq_ref_support": "rnaseq_ref_support"})

        allowed_cols = Parameters().get_features()
        keep_cols = [c for c in df.columns if c in allowed_cols]
        df = df[keep_cols]

        df.loc[:, 'patient'] = np.full(df.shape[0], patient)
        df.loc[:, 'dataset'] = np.full(df.shape[0], get_patient_group(patient))
        df.loc[:, 'train_test'] = np.full(df.shape[0], get_ml_group(patient, 'short'))

        df.loc[:, 'mut_seqid'] = \
            df.apply(lambda row:
                     row['mut_seqid'] if row['mut_seqid'].startswith(patient)
                     else patient+":"+row['mut_seqid'], axis=1)

        df.loc[:, 'seq_len'] = \
            df.apply(lambda row: len(row['mutant_seq']), axis=1)

        df.loc[:, 'DAI'] = \
            df.apply(lambda row: self.calc_dai(row['mutant_rank_netMHCpan'], row['wt_best_rank_netMHCpan']), axis=1)

        df.loc[:, 'DAI_NetMHC'] = \
            df.apply(lambda row: self.calc_dai(row['mutant_rank_netMHCpan'], row['wt_best_rank_netMHCpan']), axis=1)

        df.loc[:, 'DAI_MixMHC'] = \
            df.apply(lambda row: self.calc_dai(row['mutant_rank'], row['wt_best_rank']), axis=1)

        if 'mut_aa_coeff' in df.columns and 'wt_aa_coeff' in df.columns:
            df.loc[:, 'DAI_aa_coeff'] = \
                df.apply(lambda row: row['mut_aa_coeff'] - row['wt_aa_coeff'], axis=1)

        if 'mut_Rank_Stab' in df.columns and 'wt_Rank_Stab' in df.columns:
            df.loc[:, 'DAI_NetStab'] = \
                df.apply(lambda row: self.calc_dai(row['mut_Rank_Stab'], row['wt_Rank_Stab']), axis=1)

        df.loc[:, 'mutant_other_significant_alleles'] = \
            df.apply(lambda row: self.count_alleles(row['mutant_other_significant_alleles']), axis=1)

        if 'mut_is_binding_pos' in df.columns and 'DAI_MixMHC' in df.columns:
            df.loc[:, 'DAI_MixMHC_mbp'] = \
                df.apply(lambda row: self.calc_dai_mbp(row['DAI_MixMHC'], row['mut_is_binding_pos']), axis=1)

        if 'bestWTMatchType_I' in df.columns and 'bestWTMatchScore_I' in df.columns:
            df.loc[(df['bestWTMatchType_I'] == 'NONE') & (df['bestWTMatchScore_I'] > 0) &
                   (df['bestWTMatchOverlap_I'] == 1), 'bestWTMatchType_I'] = 'INCLUDED'
            df.loc[(df['bestWTMatchType_I'] == 'NONE') & (df['bestWTMatchScore_I'] > 0) &
                   (df['bestWTMatchOverlap_I'] < 1) & (df['bestMutationScore_I'] > 0), 'bestWTMatchType_I'] = \
                'PARTIAL_MUT'
            df.loc[(df['bestWTMatchType_I'] == 'NONE') & (df['bestWTMatchScore_I'] > 0) &
                   (df['bestWTMatchOverlap_I'] < 1) & (df['bestMutationScore_I'] == 0), 'bestWTMatchType_I'] = \
                'PARTIAL'

        return df

    def calc_dai(self, mut_value, wt_value):
        if mut_value > 0 and wt_value > 0:
            return np.log(mut_value) - np.log(wt_value)
        else:
            return 0

    def calc_dai_mbp(self, dai_value, mbp_value):
        if mbp_value:
            return dai_value
        else:
            return 0

    def count_alleles(self, alleles):
        if isinstance(alleles, str):
            return len(alleles.split(','))
        else:
            return 0

    def get_categorical_dim(self):
        return self.cat_dims

    def get_categorical_idx(self, x):
        cat_idx = None
        if self.cat_dims is not None:
            cat_idx = []
            for f in self.cat_dims:
                cat_idx.append(x.columns.get_loc(f))

        return cat_idx

    def get_cat_encoders(self):
        return self.cat_encoders

    # def neodisc_flags_short(self, row):
    #     flags = []
    #     if row['mutant_rank'] < 1:
    #         flags = flags + ['LOW_PREDICTION']
    #     if row['mutant_rank'] < 1:
    #         flags = flags + ['LOW_PREDICTION']
    #
    # def add_neodisc_flags_short(self, df):
    #     flags = df.apply(lambda row: )
