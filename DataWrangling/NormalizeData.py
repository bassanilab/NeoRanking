"""
Python code to normalize mutation or neopep data files for machine learning (ML) or plotting.
"""

from Utils.Util_fct import *
from Utils.DataManager import DataManager
from Utils.GlobalParameters import *

parser = argparse.ArgumentParser(description='Preprocess of neopep and mutation data')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-ds', '--dataset', type=str, choices=GlobalParameters.datasets_encoding,
                    help='Dataset used for encoding (NCI_train or NCI)')
parser.add_argument('-o', '--objective', type=str, choices=GlobalParameters.objectives,
                    help='Objective for normalization (ml or plot)')

if __name__ == "__main__":
    args = parser.parse_args()

    DataManager.transform_data(args.peptide_type, args.dataset, args.objective)


