import argparse

import numpy as np

from DataWrangling.DataLoader import *
from scipy import stats
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from Visualization.PCAClassifyPeptideBrowser import *
from collections import Counter
from statsmodels.stats.multitest import multipletests
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *


parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')


args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

if not args.features or len(args.features) == 0:
    features = Parameters().get_ml_features()
else:
    features = args.features

data_loader = DataLoader(transformer=None, normalizer=None, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0, max_netmhc_rank=2)

# perform leave one out on training set
dataManager = DataManager()
patients_with_data = dataManager.get_valid_patients(args.peptide_type)
if not args.patients or len(args.patients) == 0:
    patients = patients_with_data
elif args.patients[0] == 'Rosenberg':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('all')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Gartner_train':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('gartner_train')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Gartner_test':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('gartner_test')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Gartner':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('gartner')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Parkhurst':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('parkhurst')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'HiTIDE':
    annotator = NeoDiscImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients()
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'TESLA':
    annotator = TESLAImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients()
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'TESLA/HiTIDE':
    annotator = TESLAImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients()
    annotator = NeoDiscImmunogenicityAnnotatorLong(dataManager)
    patients = patients + annotator.get_patients()
    patients = [p for p in patients if p in patients_with_data]
else:
    patients = np.array(args.patients)

data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)

df_immuno = data_train.loc[data_train.apply(lambda row: row['response_type'] in args.immunogenic, axis=1)]
df_negative = data_train.loc[data_train.apply(lambda row: row['response_type'] == 'negative', axis=1)]
mv_counts_CD8 = {}
mv_counts_negative = {}

for f in data_train.columns:
    mv_counts_CD8[f] = 100*sum(df_immuno[f].isna())/df_immuno.shape[0]
    mv_counts_negative[f] = 100*sum(df_negative[f].isna())/df_negative.shape[0]


with PdfPages(args.pdf) as pp:
    mv_counts = dict(sorted(mv_counts_CD8.items(), key=lambda item: item[1], reverse=True))
    fig, ax = plt.subplots(figsize=(30, 8))
    x = np.arange(len(mv_counts))
    ax.bar(x, mv_counts.values())
    ax.set_ylabel("Missing values in % of immunogenic kmers", fontsize=15)
    ax.set_xticks(x)
    ax.set_xticklabels(mv_counts.keys())
    ax.xaxis.set_tick_params(labelsize=8, labelrotation=90)
    # plt.title("All feature missing value counts sorted", fontsize=12)
    fig.tight_layout()
    pp.savefig(fig)

    mv_counts = dict(sorted(mv_counts_negative.items(), key=lambda item: item[1], reverse=True))
    fig, ax = plt.subplots(figsize=(30, 8))
    x = np.arange(len(mv_counts))
    ax.bar(x, mv_counts.values())
    ax.set_ylabel("Missing values in % of negative kmers", fontsize=15)
    ax.set_xticks(x)
    ax.set_xticklabels(mv_counts.keys())
    ax.xaxis.set_tick_params(labelsize=8, labelrotation=90)
    # plt.title("All feature missing value counts sorted", fontsize=12)
    fig.tight_layout()
    pp.savefig(fig)
