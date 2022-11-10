import argparse
from collections import Counter

import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from os.path import exists

from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-fp', '--file_prefix', type=str, default="Feature", help='Output files prefix')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

warnings.filterwarnings("ignore")

plot_names = []
for f in args.features:
    plot_names.append(feature_dict[f])

df = pd.DataFrame({'Feature name in table': args.features, 'Feature name in plot': plot_names})

result_file = os.path.join(Parameters().get_plot_dir(),
                           "{0}_Features_{1}.txt".format(args.file_prefix, args.peptide_type))

df.to_csv(result_file, sep='\t', header=True, index=False, mode='a')
