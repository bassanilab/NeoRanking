import numpy as np
import pandas as pd
from numpy import unique
from collections import Counter

from sklearn.preprocessing import FunctionTransformer
import warnings

from Utils.GlobalParameters import *

warnings.filterwarnings(action='ignore', category=UserWarning)


class CatEncoder:

    def __init__(self, feature):
        self.feature = feature
        self.float_encoding = {}
        self.int_encoding = {}
        self.unknown_float_code = -1
        self.unknown_int_code = -1
        self.nr_classes = 0
        return

    def fit(self, x, y):
        pos_cnt = Counter([x_sel for x_sel, y_sel in zip(x, y) if y_sel == 1])
        tot_cnt = Counter(x)

        self.float_encoding = {}
        if len(tot_cnt.keys()) == 2 and (True in tot_cnt.keys()) and (False in tot_cnt.keys()):
            self.float_encoding[True] = 1.0
            self.float_encoding[False] = 0.0
            self.unknown_float_code = 0.0
        else:
            s = 0.0
            for l in tot_cnt:
                if l in pos_cnt:
                    self.float_encoding[l] = pos_cnt[l] / tot_cnt[l]
                    s += self.float_encoding[l]
                else:
                    self.float_encoding[l] = 0.0

            for l in tot_cnt:
                self.float_encoding[l] /= s

            self.unknown_float_code = 0.0

        self.float_encoding = dict(sorted(self.float_encoding.items(), key=lambda item: item[1]))
        self.nr_classes = len(self.float_encoding.keys())

        self.encode_int_from_float()

        return self

    def encode_int_from_float(self):
        self.int_encoding = {}
        idx = 1
        for label in self.float_encoding:
            self.int_encoding[label] = idx
            idx += 1
        self.unknown_int_code = 0

    def encode_float(self, label):
        if label in self.float_encoding:
            return self.float_encoding[label]
        else:
            return self.unknown_float_code

    def encode_int(self, label):
        if label in self.int_encoding:
            return self.int_encoding[label]
        else:
            return self.unknown_int_code

    def transform(self, values, cat_type='float'):
        if cat_type == 'float':
            return list(map(self.encode_float, values))
        elif cat_type == 'int':
            return list(map(self.encode_int, values))
        else:
            return values

    def get_nr_classes(self):
        return self.nr_classes

    def get_unknown(self, cat_type='float'):
        if cat_type == 'float':
            return self.unknown_float_code
        elif cat_type == 'int':
            return self.unknown_int_code
        else:
            return np.nan

    def append_to_file(self, encoding_file):
        for c in self.float_encoding.keys():
            encoding_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.feature, self.nr_classes, c, self.float_encoding[c]))
        encoding_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.feature, self.nr_classes, 'unknown', self.unknown_float_code))

    def get_encoding(self, encoding_df):
        df = encoding_df[encoding_df['Feature'] == self.feature]
        self.float_encoding = {}
        for idx in df.index:
            cat = df.loc[idx, 'Category']
            if cat == 'unknown':
                self.unknown_float_code = df.loc[idx, 'Value']
            elif cat == 'False':
                self.float_encoding[False] = df.loc[idx, 'Value']
            elif cat == 'True':
                self.float_encoding[True] = df.loc[idx, 'Value']
            elif cat == 'nan':
                self.float_encoding[np.nan] = df.loc[idx, 'Value']
            else:
                self.float_encoding[cat] = df.loc[idx, 'Value']

        self.nr_classes = len(self.float_encoding.keys())
        self.encode_int_from_float()

    @staticmethod
    def read_cat_encodings(dataset: str, peptide_type: str):
        encoding_file = GlobalParameters().get_cat_to_num_info_file(dataset, peptide_type)
        encoding_df = pd.read_csv(encoding_file, header=0, sep="\t", comment='#')
        features = encoding_df['Feature'].unique()
        encoders = {}
        for f in features:
            encoder = CatEncoder(f)
            encoder.get_encoding(encoding_df)
            encoders[f] = encoder

        return encoders

    @staticmethod
    def get_file_header():
        return "Feature\tNr_Categories\tCategory\tValue"


