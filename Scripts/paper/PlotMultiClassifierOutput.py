import argparse

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import parallel_coordinates
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-g', '--patient_groups', type=str, nargs='+', help='Patient groups displayed separately')
parser.add_argument('-t', '--tt', type=str, nargs='+', help='Patient groups displayed separately')
parser.add_argument('-a', '--alpha', type=float, default=0.02, help='Coefficient alpha in score function')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class ClassifierResults:

    def __init__(self, lines, name, index, groups=None, alpha=0.02):
        self.config = {}
        self.parse_config(lines)
        self.hyperopt_results = {}
        self.parse_hyperopt_results(lines)
        self.name = name
        self.replicate_index = index
        self.groups = groups
        self.results = None
        self.alpha = alpha
        self.parse_clf_results(lines)
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result(lines)

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
            scores = self.calc_scores(ranks)
            group = get_patient_group(row['Patient'])
            df = pd.DataFrame({'Patient': row['Patient'], 'Patient_group': group, "Ranks": ranks,
                               'Scores': scores, 'Peptide_id': peptide_ids, 'Mutant_seq': peptides, 'Gene': genes,
                               'Classifier': self.name})

            if self.results is None:
                self.results = df
            else:
                self.results = self.results.append(df, ignore_index=True)

    def calc_tot_clf_result(self, lines):
        self.sum_results['Total_score'] = sum(self.results['Scores'])
        if self.groups is not None:
            for group in self.groups:
                self.sum_results[group+"_score"] = \
                    sum(self.results.loc[self.results['Patient_group'] == group, 'Scores'])

    def get_name(self):
        return self.name

    def get_config(self):
        return self.config

    def get_results_data(self):
        return self.results

    def add_to_sum_df(self, rank_df):
        d = {'Classifier': self.name, 'Score': [self.sum_results['Total_score']]}
        for group in self.groups:
            d[group+"_score"] = [self.sum_results[group+'_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return rank_df.append(df, ignore_index=True)

    def add_to_vector_df(self, vector_df):
        if vector_df is None:
            df = self.results.loc[:, ['Patient', 'Peptide_id', 'Scores']]
            df.rename(columns={'Scores': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            return df
        else:
            df = self.results.loc[:, ['Patient', 'Peptide_id', 'Scores']]
            df.rename(columns={'Scores': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            df = vector_df.merge(df, how='outer', left_on='Peptide_id', right_on='Peptide_id', suffixes=("", "_right"))
            df = df.fillna(0.0)
            return df.drop(columns=[c for c in df.columns if c.endswith("_right")])

    def get_patients(self):
        return set(self.results['Patient'])

    def calc_scores(self, ranks):
        return [np.exp(np.multiply(-self.alpha, r)) for r in ranks]


plot_df = None
vector_df = None
patients = set()
for re in args.clf_result_files_re:
    clf_result_files = glob.glob(os.path.join(Parameters().get_pickle_dir(), re))

    for i, result_file in enumerate(clf_result_files):
        if os.path.getsize(result_file) > 0:
            with open(result_file) as file:
                fields = os.path.basename(result_file).split('_')
                name = "{0}".format(fields[0])
                clf_results = ClassifierResults([line.rstrip() for line in file], name, i,
                                                args.patient_groups, alpha=args.alpha)
                plot_df = clf_results.add_to_sum_df(plot_df)
                vector_df = clf_results.add_to_vector_df(vector_df)

plot_df = plot_df.astype({'Classifier': str, 'Score': float})
for group in args.patient_groups:
    plot_df.astype({group+"_score": float})

pdf = args.pdf.replace(".pdf", "_alpha_{0}.pdf".format(args.alpha))
with PdfPages(pdf) as pp:

    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.boxplot(x="Classifier", y="Score", data=plot_df)
    sns.swarmplot(x="Classifier", y="Score", data=plot_df, color=".25")
    plt.ylabel("sum(exp(-{0:.3f}*rank))".format(args.alpha), size=20)
    plt.xlabel("Classifier", size=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

    for group in args.patient_groups:
        fig = plt.figure(figsize=(10, 6))
        fig.clf()
        g = sns.boxplot(x="Classifier", y=group+"_score", data=plot_df)
        sns.swarmplot(x="Classifier", y=group+"_score", data=plot_df, color=".25")
        plt.ylabel("sum(exp(-{0:.3f}*rank))".format(args.alpha), size=20)
        plt.xlabel("", size=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=15)
        plt.title(group, fontsize=20)
        g.figure.tight_layout()
        pp.savefig()
        plt.close()

    fig = plt.figure(figsize=(10, 10))
    fig.clf()
    pca = PCA(n_components=2)

    X = vector_df.loc[:, [c for c in vector_df.columns if c not in ['Patient', 'Peptide_id']]].to_numpy().transpose()
    x_pca = pca.fit_transform(X)
    classifiers = [c.split("_")[0] for c in vector_df.columns if c not in ['Patient', 'Peptide_id']]
    pca_df = pd.DataFrame({'PCA_1': x_pca[:, 0], 'PCA_2': x_pca[:, 1], 'Classifier': classifiers})
    variance = pca.explained_variance_ratio_

    fg = sns.scatterplot(data=pca_df, x='PCA_1', y='PCA_2', hue='Classifier')
    plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100), size=20)
    plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100), size=20)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

