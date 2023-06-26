import numpy as np
import pandas as pd

from collections import Counter

from sklearn.preprocessing import QuantileTransformer, StandardScaler, PowerTransformer, MinMaxScaler, \
    FunctionTransformer

from Utils.GlobalParameters import GlobalParameters
from Utils.Util_fct import *
from DataWrangling.CatEncoder import CatEncoder


class DataTransformer:

    def __init__(self, peptide_type: str, objective: str, dataset: str = None, normalizer = None ):
        self.cat_dims: dict = {}
        self.cat_idx: list = []
        if dataset:
            self.cat_encoders = CatEncoder.read_cat_encodings(dataset=dataset, peptide_type=peptide_type)
        self.peptide_type: str = peptide_type
        self.normalizer = normalizer
        self.objective = objective

        return

    def apply(self, data: pd.DataFrame) -> list:
        if data is None or data.shape[0] == 0:
            return None, None, None

        if self.peptide_type == 'mutation':
            return self.load_patient_mutation(data)
        elif self.peptide_type == 'neopep':
            return self.load_patient_neopep(data)
        else:
            return None, None, None

    def load_patient_mutation(self, df):
        if self.objective == 'sel':
            df = DataTransformer.filter_rows_mutation(df)
            y = np.array(df.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
            return df, None, y

        if df.shape[0] == 0:
            return None, None, None

        y = np.array(df.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
        if self.objective == 'ml':
            X = df.loc[:, GlobalParameters.ml_features_mutation]
        else:
            X = df.copy()

        if self.objective == 'ml':
            X = self.fill_missing_values(X)

        X = self.normalize(X)

        if self.objective == 'ml':
            X = self.encode_cat_features(X)

        X = X.astype(get_processed_types(self.peptide_type, self.objective))

        return df, X, y

    @staticmethod
    def filter_rows_mutation(df):
        if df.shape[0] > 0:
            df = df.loc[df.mutation_type.apply(lambda r: r == 'SNV')]

        if df.shape[0] > 0 and GlobalParameters.max_netmhc_rank > 0:
            df = df.loc[df.mut_Rank_EL_0.apply(lambda r: r < GlobalParameters.max_netmhc_rank)]

        if df.shape[0] > 0 and len(GlobalParameters.excluded_genes) > 0:
            df = df.loc[df.gene.apply(lambda g: g not in GlobalParameters.excluded_genes)]

        df.reset_index(inplace=True, drop=True)

        return df

    def load_patient_neopep(self, df):
        if self.objective == 'sel':
            df = DataTransformer.filter_rows_neopep(df)
            y = np.array(df.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
            return df, None, y

        if df.shape[0] == 0:
            return None, None, None

        y = np.array(df.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
        if self.objective == 'ml':
            X = df.loc[:, GlobalParameters.ml_features_neopep]
        else:
            X = df.copy()

        if self.objective == 'ml':
            X = self.fill_missing_values(X)

        X = self.normalize(X)

        if self.objective == 'ml':
            X = self.encode_cat_features(X)

        X = X.astype(get_processed_types(self.peptide_type, self.objective))

        return df, X, y

    @staticmethod
    def filter_rows_neopep(df):
        if df.shape[0] > 0:
            df = df.loc[df.mutation_type.apply(lambda r: r == 'SNV')]

        if df.shape[0] > 0 and GlobalParameters.max_netmhc_rank > 0:
            df = df.loc[df.mutant_rank_netMHCpan.apply(lambda r: r < GlobalParameters.max_netmhc_rank)]

        if df.shape[0] > 0 and len(GlobalParameters.excluded_genes) > 0:
            df = df.loc[df.gene.apply(lambda g: g not in GlobalParameters.excluded_genes)]

        df.reset_index(inplace=True, drop=True)

        return df

    def encode_cat_features(self, x_):
        self.cat_dims = {}
        self.cat_idx = []
        types = GlobalParameters.feature_types_neopep if self.peptide_type == 'neopep' \
            else GlobalParameters.feature_types_mutation
        for i, c in enumerate(x_.columns):
            if is_cat_type(types[c]) or types[c] == 'bool':
            # should be if is_cat_type(x_[c].dtype.name) or is_discrete_ordered_type(x_[c].dtype):
                encoder = self.cat_encoders[c]
                self.cat_dims[c] = encoder.get_nr_classes()
                self.cat_idx.append(i)
                if GlobalParameters.cat_type == 'float':
                    x_.loc[:, c] = encoder.transform(x_[c].values, 'float')
                    x_ = x_.astype({c: 'float64'})
                elif GlobalParameters.cat_type == 'int':
                    x_.loc[:, c] = encoder.transform(x_[c].values, 'int')
                    x_ = x_.astype({c: 'int32'})

        return x_

    def get_categorical_dim(self):
        return self.cat_dims

    def get_categorical_idx(self):
        return self.cat_idx

    def get_cat_encoders(self):
        return self.cat_encoders

    def fill_missing_values(self, x_):

        value_dict = {}

        types = GlobalParameters.feature_types_neopep if self.peptide_type == 'neopep' \
            else GlobalParameters.feature_types_mutation
        order_rels = GlobalParameters.ml_feature_mv_neopep if self.peptide_type == 'neopep' \
            else GlobalParameters.ml_feature_mv_mutation
        for i, c in enumerate(x_.columns):
            if not is_cat_type(types[c]):
                order_rel = order_rels[c]
                if order_rel == 'min':
                    fill_val = np.min(x_[c])
                elif order_rel == 'max':
                    fill_val = np.max(x_[c])
                elif order_rel == 'cnt':
                    counter = Counter(x_[c])
                    fill_val = counter.most_common(1)[0][0]
                elif type(order_rel) == float:
                    fill_val = order_rel
                else:
                    fill_val = np.mean(x_[c])

                value_dict[c] = fill_val

        x_.fillna(value=value_dict, inplace=True)

        return x_

    def normalize(self, x_):
        if self.normalizer is None:
            return x_

        for i, c in enumerate(x_.columns):
            if is_cont_type(x_[c].dtype.name):
                if type(self.normalizer) is dict:
                    if c in self.normalizer:
                        norm_transform = self.normalizer[c]
                    else:
                        norm_transform = None
                else:
                    norm_transform = self.normalizer

                if norm_transform:
                    v = x_[c].to_numpy().reshape(-1, 1)
                    if type(norm_transform) is FunctionTransformer and norm_transform.func.__name__ == 'log10':
                        v[v <= 0] = min(v[v > 0])/10
                    x_.loc[:, c] = norm_transform.fit_transform(v)

        return x_

    @staticmethod
    def get_normalizer_(normalizer_tag):
        if normalizer_tag == 'q':
            return QuantileTransformer()
        elif normalizer_tag == 'z':
            return StandardScaler()
        elif normalizer_tag == 'p':
            return PowerTransformer()
        elif normalizer_tag == 'i':
            return MinMaxScaler()
        elif normalizer_tag == 'l':
            return FunctionTransformer(np.log10, inverse_func=lambda x: np.power(10, x), validate=False, check_inverse=True)
        elif normalizer_tag == 'a':
            return FunctionTransformer(np.arcsinh, inverse_func=np.sinh, validate=False, check_inverse=True)
        elif normalizer_tag == 'n':
            return None
        elif normalizer_tag == 'plot':
            d = {}
            for k, v in GlobalParameters.plot_normalization.items():
                d[k] = DataTransformer.get_normalizer_(v)
            return d

    @staticmethod
    def get_normalizer(objective):
        if objective == 'ml':
            return DataTransformer.get_normalizer_(GlobalParameters.normalizer)
        elif objective == 'plot':
            return DataTransformer.get_normalizer_('plot')
        else:
            return None

    @staticmethod
    def get_normalizer_name(normalizer):
        if type(normalizer) is QuantileTransformer:
            return "QuantileTransformer"
        elif type(normalizer) is PowerTransformer:
            return "PowerTransformer"
        elif type(normalizer) is StandardScaler:
            return "StandardScaler"
        elif type(normalizer) is MinMaxScaler:
            return "MinMaxScaler"
        elif type(normalizer) is FunctionTransformer and normalizer.func.__name__ == 'log10':
            return "log10"
        elif type(normalizer) is FunctionTransformer and normalizer.func.__name__ == 'arcsinh':
            return "arcsinh"
        elif normalizer is None:
            return "None"

