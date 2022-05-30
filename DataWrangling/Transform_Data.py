import numpy as np
import pandas as pd
from numpy import unique
from collections import Counter

from sklearn.preprocessing import FunctionTransformer
import warnings

from Utils.Parameters import *

warnings.filterwarnings(action='ignore', category=UserWarning)


class Encoder:

    def __init__(self, feature):
        self.feature = feature
        self.encoding = {}
        self.unknown_code = -1
        self.nr_classes = 0
        return

    def fit(self, x, y):
        pos_cnt = Counter([x_sel for x_sel, y_sel in zip(x, y) if y_sel == 1])
        tot_cnt = Counter(x)

        self.encoding = {}
        if len(tot_cnt.keys()) == 2 and (True in tot_cnt.keys()) and (False in tot_cnt.keys()):
            self.encoding[True] = 1.0
            self.encoding[False] = 0.0
            self.unknown_code = 0.5
        else:
            for l in tot_cnt:
                if l in pos_cnt:
                    self.encoding[l] = pos_cnt[l]/tot_cnt[l]
                else:
                    self.encoding[l] = 0.0

            self.unknown_code = 0.0

        self.nr_classes = len(self.encoding.keys())

        return self

    def encode(self, label):
        if label in self.encoding:
            return self.encoding[label]
        else:
            return self.unknown_code

    def transform(self, values):
        return list(map(self.encode, values))

    def get_nr_classes(self):
        return self.nr_classes

    def append_to_file(self, encoding_file):
        for c in self.encoding.keys():
            encoding_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.feature, self.nr_classes, c, self.encoding[c]))
        encoding_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.feature, self.nr_classes, 'unknown', self.unknown_code))

    def read_from_file(self, encoding_df):
        df = encoding_df[encoding_df['Feature'] == self.feature]
        self.encoding = {}
        for idx in df.index:
            cat = df.loc[idx, 'Category']
            if cat == 'unknown':
                self.unknown_code = df.loc[idx, 'Value']
            elif cat == 'False':
                self.encoding[False] = df.loc[idx, 'Value']
            elif cat == 'True':
                self.encoding[True] = df.loc[idx, 'Value']
            elif cat == 'nan':
                self.encoding[np.nan] = df.loc[idx, 'Value']
            else:
                self.encoding[cat] = df.loc[idx, 'Value']

        self.nr_classes = len(self.encoding.keys())

    @staticmethod
    def get_file_header():
        return "Feature\tNr_Categories\tCategory\tValue"


class DataTransformer:

    def __init__(self):
        return

    @staticmethod
    def fill_missing_values(df):

        df = DataTransformer.impute_rnaseq_cov(df)
        df['rnaseq_TPM'].fillna(0, inplace=True)

        parameters = Parameters()
        value_dict = {}
        all_features = parameters.get_numerical_features()+parameters.get_ordinal_features()
        num_features = [f for f in df.columns if f in all_features]

        for f in num_features:
            order_rel = parameters.get_order_relation(f)
            min_v = np.min(df[f])
            max_v = np.max(df[f])
#            dv = (max_v-min_v)/3
            dv = 0
            if order_rel == '>':
                fill_val = min_v-dv
            elif order_rel == '<':
                fill_val = max_v+dv
            elif order_rel == '=':
                if f in parameters.get_ordinal_features():
                    counter = Counter(df[f])
                    fill_val = counter.most_common(1)[0][0]
                else:
                    fill_val = np.mean(df[f])

            value_dict[f] = fill_val

        df.fillna(value=value_dict, inplace=True)

        return df

    @staticmethod
    def impute_rnaseq_cov(df):
        for i in df.index:
            if np.isnan(df.loc[i, 'rnaseq_ref_support']):
                if df.loc[i, 'rnaseq_gene_expression_quartile'] == 0:
                    df.loc[i, 'rnaseq_ref_support'] = 100.0
                elif df.loc[i, 'rnaseq_gene_expression_quartile'] == 1:
                    df.loc[i, 'rnaseq_ref_support'] = 100.0
                elif df.loc[i, 'rnaseq_gene_expression_quartile'] == 2:
                    df.loc[i, 'rnaseq_ref_support'] = 100.0
                elif df.loc[i, 'rnaseq_gene_expression_quartile'] == 3:
                    df.loc[i, 'rnaseq_ref_support'] = 89.0
                elif df.loc[i, 'rnaseq_gene_expression_quartile'] == 4:
                    df.loc[i, 'rnaseq_ref_support'] = 67.0

                df.loc[i, 'rnaseq_alt_support'] = 100.0 - df.loc[i, 'rnaseq_ref_support']

        return df

    @staticmethod
    def same_allele(alleles_wt, alleles_mut):

        if str(alleles_wt) == 'nan' or str(alleles_mut) == 'nan':
            return False
        la_wt = str(alleles_wt).split(",")
        for a in la_wt:
            if a in alleles_mut: return True

        return False

    @staticmethod
    def normalize(df, normalizer):
        num_cols = [c for c in df.columns
                    if c in Parameters().get_numerical_features() or c in Parameters().get_ordinal_features()]

        X = df.copy()
        for c in num_cols:
            if type(normalizer) is dict:
                if c in normalizer:
                    norm_transform = normalizer[c]
                else:
                    norm_transform = None
            else:
                norm_transform = normalizer

            if norm_transform is not None:
                x = df[c].to_numpy().reshape(-1, 1)
                if type(norm_transform) is FunctionTransformer and norm_transform.func.__name__ == 'log10':
                    x[x <= 0] = min(x[x > 0])/10
                X.loc[:, c] = norm_transform.fit_transform(x)

        return X
