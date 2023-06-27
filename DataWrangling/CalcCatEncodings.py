"""
Python code to create encoding of categorical features. Encodings can be trained for mutation and neopep
data and on any training dataset. The training dataset should be different from datasets used for testing
the ML algorithms.
"""

import argparse
import sys

from DataWrangling.CatEncoder import CatEncoder
from Utils.Util_fct import *
from Utils.GlobalParameters import *
from Utils.DataManager import DataManager

parser = argparse.ArgumentParser(description='Preprocess of neopep and mutation data')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-ds', '--dataset', type=str, choices=GlobalParameters.datasets_encoding,
                    help='Dataset used for encoding (NCI_train or NCI)')

if __name__ == "__main__":
    args = parser.parse_args()
    # GlobalParameters has predefined filenames for NCI and NCI_train datasets
    with open(GlobalParameters.get_cat_to_num_info_file(args.dataset, args.peptide_type), mode='w') \
            as encoding_file:

        for arg in vars(args):
            encoding_file.write(f"#{arg}={getattr(args, arg)}\n")
            print(f"{arg}={getattr(args, arg)}")

        data_train = DataManager.filter_data(peptide_type=args.peptide_type, dataset=args.dataset,
                                             response_types=['CD8', 'negative'])

        encoding_file.write(CatEncoder.get_file_header() + "\n")

        ml_features = GlobalParameters.ml_features_neopep if args.peptide_type == 'neopep' \
            else GlobalParameters.ml_features_mutation
        y = np.array(data_train.response_type.apply(lambda row: int(row == 'CD8')), dtype=int)

        for f in data_train.columns:
            if f in ml_features and (data_train[f].dtype.name == 'category' or data_train[f].dtype == bool):
                print("Encoding feature {0} ...".format(f))
                l_enc = CatEncoder(f)
                l_enc.fit(data_train[f].values, y)
                l_enc.append_to_file(encoding_file)
