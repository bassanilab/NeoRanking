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
parser.add_argument('-r', '--clf_result_files', type=str, nargs='+', help='Comma separated list of clf result files')
parser.add_argument('-n', '--names', type=str, nargs='+', help='Comma separated list of clf test names')
parser.add_argument('-nd', '--neodisc', dest='neodisc', action='store_true', help='Include neodisc prioritization')
parser.add_argument('-ga', '--gartner', dest='gartner', action='store_true', help='Include Gartner et al. prioritization')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class NeoDiscResults:

    def __init__(self):
        self.mgr = DataManager()
        self.col_name = 'NeoDisc Ranking'
        self.annotators = {}

    def get_annotator(self, peptide_type, patient_group):
        if peptide_type == 'long':
            if patient_group == 'Rosenberg':
                annotator = RosenbergImmunogenicityAnnotatorLong(self.mgr)
            elif patient_group == 'TESLA':
                annotator = TESLAImmunogenicityAnnotatorLong(self.mgr)
            else:
                annotator = NeoDiscImmunogenicityAnnotatorLong(self.mgr)
        else:
            if patient_group == 'Rosenberg':
                annotator = RosenbergImmunogenicityAnnotatorShort(self.mgr)
            elif patient_group == 'TESLA':
                annotator = TESLAImmunogenicityAnnotatorShort(self.mgr)
            else:
                annotator = NeoDiscImmunogenicityAnnotatorShort(self.mgr)

        return annotator

    def get_CD8_ranks(self, patient, peptide_data, peptide_type='long'):
        neodisc_data = self.mgr.get_processed_data(patient, 'rt', peptide_type)
        ranks = []
        for index, row in peptide_data.iterrows():
            idx = np.where(neodisc_data['peptide_id'] == row['Peptide_id'])
            if len(idx) > 0:
                ranks = np.append(ranks, idx[0][0] + 1)
            else:
                ranks = np.append(ranks, np.nan)

        return pd.Series(ranks, name=self.col_name, dtype=int)

    def get_col_name(self):
        return self.col_name

    def add_to_plot_dfs(self, plot_dfs, patient, peptide_type='long'):
        if patient not in plot_dfs:
            plot_dfs[patient] = self.get_CD8_ranks(patient, plot_dfs[patient], peptide_type)
        else:
            plot_dfs[patient] = \
                pd.concat([plot_dfs[patient], self.get_CD8_ranks(patient, plot_dfs[patient], peptide_type)], axis=1)


class ClassifierResults:

    def __init__(self, lines, name):
        self.config = {}
        self.parse_config(lines)
        self.results = self.parse_clf_results(lines)
        self.name = name

    def parse_config(self, lines):
        for l in lines:
            if l.startswith("Patient"):
                break
            fields = l.split("=")
            if len(fields) > 1:
                self.config[fields[0]] = fields[1]

    def parse_clf_results(self, lines):
        result_value_list = []
        header = None
        for l in lines:
            if l.startswith("Patient"):
                header = l.split("\t")
                #header.append("Ranking_score")
                continue
            if l.startswith("nr_patients"):
                break

            if header is not None:
                values = np.array(l.split("\t"))
                result_value_list.append(pd.Series(values))

        results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
        results.columns = header
        return results

    def get_name(self):
        return self.name

    def get_config(self):
        return self.config

    def get_results_data(self):
        return self.results

    def add_to_plot_dfs(self, plot_dfs):
        for index, row in self.results.iterrows():
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
                plot_dfs[row['Patient']] = df
            else:
                plot_dfs[row['Patient']] = \
                    plot_dfs[row['Patient']].merge(df, left_on='Peptide_id', right_on='Peptide_id',
                                                   suffixes=('', '_right')).filter(regex="^(?!.*_right)", axis=1)

    def get_patients(self):
        return set(self.results['Patient'])


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
