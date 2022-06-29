import argparse

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import parallel_coordinates
from matplotlib import pyplot as plt
import seaborn as sns

from Utils.Util_fct import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorShort import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorShort import *

parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
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
        self.results = {}
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

        results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
        results.columns = header

        for index, row in results.iterrows():
            rank_strs = row['CD8_ranks'].split(',')
            peptide_id_strs = row['CD8_peptide_idx'].split(',')
            peptide_strs = row['CD8_mut_seqs'].split(',')
            gene_strs = row['CD8_genes'].split(',')
            nr_immuno = int(row['Nr_immunogenic'])
            ranks = []
            peptide_ids = []
            peptides = []
            genes = []
            for i in range(0, nr_immuno):
                if i < len(rank_strs) and len(rank_strs[i]) > 0:
                    ranks = np.append(ranks, int(rank_strs[i]))
                    peptide_ids = np.append(peptide_ids, peptide_id_strs[i])
                    peptides = np.append(peptides, peptide_strs[i])
                    genes = np.append(genes, gene_strs[i])
                else:
                    ranks = np.append(ranks, int(row['Nr_peptides']))
                    peptide_ids = np.append(peptide_ids, "")
                    peptides = np.append(peptides, "")
                    genes = np.append(genes, "")

            ranks = np.array(ranks, dtype='int32')
            df = pd.DataFrame({self.name: ranks, 'Peptide_id': peptide_ids, 'Mutant_seq': peptides, 'Gene': genes})

            if row['Patient'] not in plot_dfs:
                self.results[row['Patient']] = df

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


plot_df = pd.DataFrame({'Train score': [], 'Test score': [], 'Name': [], 'top_20': [], 'top_50': [], 'top_100': []})
for re in args.clf_result_files_re:
    clf_result_files = glob.glob(os.path.join(Parameters().get_pickle_dir(), re))
    patients = set()
    i = 0
    for result_file in clf_result_files:
        with open(result_file) as file:
            fields = os.path.basename(result_file).split('_')
            name = "{0}".format(fields[0])
            clf_results = ClassifierResults([line.rstrip() for line in file], name)
            plot_df = clf_results.add_to_plot_df(plot_df)


plot_dfs = {}
patients = set()
for result_file, name in zip(args.clf_result_files, args.names):
    with open(result_file) as file:
        print(file)
        clf_results = ClassifierResults([line.rstrip() for line in file], name)
        clf_results.add_to_plot_dfs(plot_dfs)
        patients = set.union(patients, clf_results.get_patients())

classifiers = args.names

if args.neodisc:
    neodisc_results = NeoDiscResults()
    col_name = neodisc_results.get_col_name()
    classifiers = np.append(classifiers, col_name)
    for p in patients:
        neodisc_results.add_to_plot_dfs(plot_dfs, p, args.peptide_type)

plot_df = []
mgr = DataManager()
for patient in plot_dfs.keys():
    df = plot_dfs[patient]
    rank_cols = [c for c in df.columns if c.startswith('Ranks_')]
    for c in rank_cols:
        max_rank = df[c].max(skipna=True)
        df.loc[df[c].isna(), c] = max_rank
        df.loc[:, c] = df[c].astype('int32')
    df['Patient'] = patient
    patient_group = get_patient_group(patient)
    df['Patient_group'] = patient_group
    plot_dfs[patient] = df

plot_df = pd.concat(plot_dfs.values())

if args.gartner:
    gartner_ranks = {3703: [2, 42, 70, 2522], 3881: [6], 3995: [2], 4007: [1], 4014: [84], 4032: [12, 32, 48], 4166: [9],
                     4242: [242], 4268: [2, 5], 4271: [15, 114], 4283: [21], 4310: [1, 8], 4323: [8], 4324: [24392],
                     4350: [2, 40, 1641], 4359: [3, 49]}

    plot_df.insert(len(classifiers), 'Gartner Ranking', -1)
    for p in plot_df['Patient'].unique():
        plot_df.loc[plot_df['Patient'] == p, 'Gartner Ranking'] = gartner_ranks[int(p)]

    classifiers = np.append(classifiers, 'Gartner Ranking')

plot_df_num = plot_df.loc[:, classifiers]


top_20 = plot_df_num.apply(lambda c: sum(c <= 20), axis=0)
top_50 = plot_df_num.apply(lambda c: sum(c <= 50), axis=0)
top_100 = plot_df_num.apply(lambda c: sum(c <= 100), axis=0)
med_rank = plot_df_num.apply(lambda c: c.median(), axis=0)
mean_rank = plot_df_num.apply(lambda c: c.mean(), axis=0)
exp_score_df = plot_df_num.transform(lambda v: np.exp(np.multiply(-0.02, v)), axis=1)
exp_scores = exp_score_df.apply(lambda c: sum(c), axis=0)

sum_plot_df = pd.concat([top_20, top_50, top_100], axis=1).transpose()
sum_plot_df.insert(0, 'Top N', [20, 50, 100])
sum_plot_df = pd.melt(sum_plot_df, id_vars=['Top N'], value_vars=classifiers, var_name='Method',
                      value_name='# CD8+ 8-12mers')

with PdfPages(args.pdf) as pp:
    patients = plot_df['Patient'].unique()
    for p in patients:
        df = plot_df.loc[plot_df.Patient == p, np.append(classifiers, 'Patient')]
        fig = plt.figure(figsize=(10, 6))
        fig.clf()
        g = parallel_coordinates(df, 'Patient', color=['b'])
        plt.yscale('log')
        plt.ylabel("Rank", size=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=15)
        g.figure.tight_layout()
        pp.savefig()
        plt.close()

    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.barplot(x='Method', y='# CD8+ 8-12mers', hue='Top N', data=sum_plot_df)
    plt.ylabel("CD8+ 8-12 mers in top N", size=20)
    plt.xlabel("", size=15)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.barplot(x=classifiers, y=exp_scores)
    plt.ylabel("Sum of exp(-rank/50)", size=15)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()
