"""
Python code to select rows and columns of mutation or neopep data files relevant for machine learning (ML) or plotting.
"""

from Utils.Util_fct import *
from Utils.DataManager import DataManager
from Utils.GlobalParameters import *

parser = argparse.ArgumentParser(description='Select of neopep and mutation data for machine learning or plotting')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')

if __name__ == "__main__":
    args = parser.parse_args()

    DataManager.select_ml_data(args.peptide_type)


