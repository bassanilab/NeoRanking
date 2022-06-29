import argparse
import ast
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import parallel_coordinates
from matplotlib import pyplot as plt
import seaborn as sns
import glob

from Utils.Util_fct import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorShort import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorShort import *

parser = argparse.ArgumentParser(description='Plot score on training set versus score on test set')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-re', '--clf_result_files_re', type=str, help='Comma separated list of clf result files')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class ClassifierResults:

    def __init__(self, lines, name):
        self.config = {}
        self.parse_config(lines)
        self.hyperopt_results = {}
        self.parse_hyperopt_results(lines)
        self.results = pd.DataFrame()
        self.parse_clf_results(lines)
        self.sum_results = pd.Series(dtype=float)
        self.parse_tot_clf_result(lines)
        self.name = name

    def parse_config(self, lines):
        for l in lines:
            if l.startswith("Hyperopt"):
                break
            fields = l.split("=")
            if len(fields) > 1:
                self.config[fields[0]] = fields[1]

    def parse_hyperopt_results(self, lines):
        for l in lines:
            if l.startswith("Patient"):
                break
            if l.startswith("Hyperopt"):
                fields = l.replace("Hyperopt: ", "").split("; ")
                for f in fields:
                    k, v = f.split('=')
                    if k == 'Params':
                        self.hyperopt_results[k] = ast.literal_eval(v)
                    else:
                        self.hyperopt_results[k] = float(v)

    def parse_clf_results(self, lines):
        result_value_list = []
        header = None
        for l in lines:
            if l.startswith("Patient"):
                header = l.split("\t")
                continue
            if l.startswith("nr_patients"):
                break

            if header is not None:
                values = np.array(l.split("\t"))
                result_value_list.append(pd.Series(values))

        self.results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
        self.results.columns = header

    def parse_tot_clf_result(self, lines):
        header = None
        for l in lines:
            if l.startswith("nr_patients"):
                header = l.split("\t")
                continue

            if header is not None:
                values = np.array(l.split("\t"))
                self.sum_results = pd.Series(values, index=header)

    def get_name(self):
        return self.name

    def get_config(self):
        return self.config

    def get_results_data(self):
        return self.results

    def add_to_plot_df(self, plot_df):
        return plot_df.append(pd.DataFrame([[self.hyperopt_results['Score'], self.sum_results['score_train']]],
                                           columns=['Train score', 'Test score']), ignore_index=True)

    def get_patients(self):
        return set(self.results['Patient'])


clf_result_files = glob.glob(os.path.join(Parameters().get_pickle_dir(), args.clf_result_files_re))
plot_df = pd.DataFrame({'Train score': [], 'Test score': []})
patients = set()
i = 0
for result_file in clf_result_files:
    with open(result_file) as file:
        fields = os.path.basename(result_file).split('_')
        name = "{0}_{1}".format(fields[0], i)
        clf_results = ClassifierResults([line.rstrip() for line in file], name)
        plot_df = clf_results.add_to_plot_df(plot_df)


with PdfPages(args.pdf) as pp:

    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.scatterplot(x='Train score', y='Test score', data=plot_df)
    plt.xlabel("Train score", size=20)
    plt.ylabel("Test score", size=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

