import argparse

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import parallel_coordinates
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import re


from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
parser.add_argument('-d', '--data_dir', type=str, help='Directory containing clf results')
parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-g', '--patient_group_plot_names', type=str, default='', help='Patient groups plot names')
parser.add_argument('-c', '--classifier_plot_names', type=str, default='', help='Classifier plot names')
parser.add_argument('-a', '--alpha', type=float, default=0.02, help='Coefficient alpha in score function')
parser.add_argument('-nd', '--neodisc', dest='neodisc', action='store_true', help='Include neodisc comparison')
parser.add_argument('-ga', '--gartner', dest='gartner', action='store_true', help='Include Gartner et al. comparison')
parser.add_argument('-r', '--rotation', type=float, default=0.0, help='x-axis label rotation')
parser.add_argument('-ls', '--label_size', type=float, default=30.0, help='Axis label size')
parser.add_argument('-ts', '--tick_size', type=float, default=30.0, help='Axis tick size')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class ClassifierResults:

    def __init__(self, base_name, index, plot_name, alpha=0.02, included_patients=None):
        self.config = {}
        self.name = plot_name
        self.hyperopt_results = {}
        self.parse_config(base_name)
        self.replicate_index = index
        self.groups = set()
        self.included_patients = included_patients
        self.results = None
        self.alpha = alpha
        self.parse_clf_results(base_name)
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()

    def parse_config(self, base_name):
        train_file = base_name+"_train.txt"
        with open(train_file) as file:
            lines = [line for line in file]
            for l in lines:
                if l.startswith("Hyperopt"):
                    break
                fields = l.split("=")
                if len(fields) > 1:
                    self.config[fields[0]] = fields[1]

            self.parse_hyperopt_results(lines)

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

    def parse_clf_results(self, base_name):
        test_file = base_name+"_test.txt"
        with open(test_file) as file:
            lines = [line for line in file]
            result_value_list = []
            header = None
            for l in lines:
                l = l.rstrip()
                if l.startswith("Patient"):
                    header = l.split("\t")
                    continue

                if header is not None:
                    values = np.array(l.split("\t"))
                    result_value_list.append(pd.Series(values))

            results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
            results.columns = header

            for index, row in results.iterrows():
                if self.included_patients and row['Patient'] not in self.included_patients:
                    continue
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

                max_rank = int(row['Nr_peptides'])
                ranks = np.array(ranks, dtype='int32')
                scores = self.calc_scores(ranks)
                group = get_patient_group(row['Patient'])
                self.groups.add(group)
                df = pd.DataFrame({'Patient': row['Patient'], 'Patient_group': group, "Rank": ranks,
                                   'Score': scores, 'Peptide_id': peptide_ids, 'Mutant_seq': peptides, 'Gene': genes,
                                   'Classifier': self.name, 'Max_rank': max_rank})

                if self.results is None:
                    self.results = df
                else:
                    self.results = self.results.append(df, ignore_index=True)

    def calc_tot_clf_result(self):
        self.sum_results['Total_score'] = sum(self.results['Score'])
        if self.groups is not None:
            for group in self.groups:
                self.sum_results[group+"_score"] = \
                    sum(self.results.loc[self.results['Patient_group'] == group, 'Score'])

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
            df = self.results.loc[:, ['Patient', 'Peptide_id', 'Score']]
            df.rename(columns={'Score': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            return df
        else:
            df = self.results.loc[:, ['Patient', 'Peptide_id', 'Score']]
            df.rename(columns={'Score': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            df = vector_df.merge(df, how='outer', left_on='Peptide_id', right_on='Peptide_id', suffixes=("", "_right"))
            df = df.fillna(0.0)
            return df.drop(columns=[c for c in df.columns if c.endswith("_right")])

    def add_to_topN_df(self, top_n_df, patient_group):
        d = {'Classifier': self.name, 'Replicate': self.replicate_index,
             'Top N': [20, 50, 100],
             'Neo_pep_imm count': [self.get_topN_count(max_rank=20, patient_group=patient_group),
                                   self.get_topN_count(max_rank=50, patient_group=patient_group),
                                   self.get_topN_count(max_rank=100, patient_group=patient_group)]
             }
        df = pd.DataFrame(d)
        if top_n_df is None:
            return df
        else:
            return top_n_df.append(df, ignore_index=True)

    def add_to_patient_groups(self, patient_groups):
        return patient_groups.union(self.groups)

    def calc_scores(self, ranks):
        return [np.exp(np.multiply(-self.alpha, r)) for r in ranks]

    def get_topN_count(self, max_rank, patient_group):
        if patient_group == 'all':
            return self.results.apply(lambda row: row['Rank'] <= max_rank, axis=1).sum()
        else:
            return self.results.loc[self.results['Patient_group'] == patient_group, :]\
                .apply(lambda row: row['Rank'] <= max_rank, axis=1).sum()


class NeoDiscResults(ClassifierResults):

    def __init__(self, peptide_data, plot_name='NeoDisc', alpha=0.02, peptide_type='short'):
        self.name = plot_name
        self.groups = set()
        self.peptide_data = peptide_data
        self.peptide_type = peptide_type
        self.replicate_index = 0
        self.results = None
        self.alpha = alpha
        self.parse_clf_results()
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()

    def parse_clf_results(self):
        data_loader = DataLoader(response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                             immunogenic=['CD8', 'CD4/CD8'], mutation_types=['SNV'], min_nr_immuno=0)
        patient = None
        for index, row in self.peptide_data.iterrows():
            if row['Patient'] != patient:  # reload when new patient
                patient = row['Patient']
                data, X, y = data_loader.load_patients(patient, 'rt', 'short', verbose=True)

            idx = np.where(data['peptide_id'] == row['Peptide_id'])[0]
            if len(idx) > 0:
                rank = min(idx) + 1
            else:
                rank = row['Max_rank']

            score = self.calc_scores([rank])
            df = pd.DataFrame({'Patient': row['Patient'], 'Patient_group': row['Patient_group'], "Rank": rank,
                               'Score': score, 'Peptide_id': row['Peptide_id'], 'Mutant_seq': row['Mutant_seq'],
                               'Gene': row['Gene'], 'Classifier': self.name, 'Max_rank': row['Max_rank']})

            self.groups.add(row['Patient_group'])
            if self.results is None:
                self.results = df
            else:
                self.results = self.results.append(df, ignore_index=True)


class GartnerResults(ClassifierResults):

    gartner_ranks = \
        {'3703': [2, 42, 70, 2522], '3881': [6], '3995': [2], '4007': [1], '4014': [84], '4032': [12, 32, 48],
         '4166': [9], '4242': [242], '4268': [2, 5], '4271': [15, 114], '4283': [21], '4310': [1, 8], '4323': [8],
         '4324': [24392], '4350': [2, 40, 1641], '4359': [3, 49]}

    def __init__(self, peptide_data, plot_name='NCI ranking', alpha=0.02, peptide_type='short'):
        self.name = plot_name
        self.groups = set('NCI')
        self.peptide_data = peptide_data
        self.peptide_type = peptide_type
        self.replicate_index = 0
        self.results = None
        self.alpha = alpha
        self.parse_clf_results()
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()

    @staticmethod
    def get_gartner_test_patients():
        return set(GartnerResults.gartner_ranks)

    def parse_clf_results(self):
        patients = GartnerResults.get_gartner_test_patients().intersection(self.peptide_data['Patient'].unique())
        for p in patients:
            ranks = GartnerResults.gartner_ranks[p]
            scores = self.calc_scores(ranks)
            max_rank = self.peptide_data.loc[self.peptide_data['Patient'] == p, 'Max_rank'].unique()[0]
            df = pd.DataFrame({'Patient': p, 'Patient_group': 'NCI', "Rank": ranks,
                               'Score': scores, 'Peptide_id': 'unknown', 'Mutant_seq': 'unknown',
                               'Gene': 'unknown', 'Classifier': self.name, 'Max_rank': max_rank})

            if self.results is None:
                self.results = df
            else:
                self.results = self.results.append(df, ignore_index=True)


included_patients = None
if args.gartner:
    included_patients = GartnerResults.get_gartner_test_patients()

parse_neodisc = True
parse_gartner = True
plot_df = None
vector_df = None
topN_dfs = None
patient_groups = set()
plt_name_dict = ast.literal_eval(args.classifier_plot_names)
for j, regexp in enumerate(args.clf_result_files_re):
    clf_result_regex = glob.glob(os.path.join(args.data_dir, regexp))

    for i, result_file in enumerate(clf_result_regex):
        if os.path.getsize(result_file) > 0:
            base_file_name = re.sub("_test.txt$", "", result_file)
            clf_results = ClassifierResults(base_file_name, i, plt_name_dict[j], alpha=args.alpha,
                                            included_patients=included_patients)
            if topN_dfs is None:
                topN_dfs = {'all': None}
                for g in clf_results.groups:
                    topN_dfs[g] = None
            plot_df = clf_results.add_to_sum_df(plot_df)
            vector_df = clf_results.add_to_vector_df(vector_df)
            patient_groups = clf_results.add_to_patient_groups(patient_groups)
            for ds in topN_dfs:
                topN_dfs[ds] = clf_results.add_to_topN_df(topN_dfs[ds], ds)

            if args.neodisc and parse_neodisc:
                neodisc_results = NeoDiscResults(clf_results.results, alpha=args.alpha)
                plot_df = neodisc_results.add_to_sum_df(plot_df)
                vector_df = neodisc_results.add_to_vector_df(vector_df)
                for ds in topN_dfs:
                    topN_dfs[ds] = neodisc_results.add_to_topN_df(topN_dfs[ds], ds)
                parse_neodisc = False

            if args.gartner and parse_gartner:
                gartner_results = GartnerResults(clf_results.results, alpha=args.alpha)
                plot_df = gartner_results.add_to_sum_df(plot_df)
                for ds in topN_dfs:
                    topN_dfs[ds] = gartner_results.add_to_topN_df(topN_dfs[ds], ds)
                parse_gartner = False


plot_df = plot_df.astype({'Classifier': str, 'Score': float})
for group in patient_groups:
    plot_df.astype({group+"_score": float})

patient_group_dict = ast.literal_eval(args.patient_group_plot_names)

fig_ids = ['a', 'b', 'c', 'd', 'e']
fig_idx = 0

plt.rcParams['text.usetex'] = True
fig = plt.figure(figsize=(10, 6))
g = sns.boxplot(x="Classifier", y="Score", data=plot_df)
sns.swarmplot(x="Classifier", y="Score", data=plot_df, color=".25")
lbs = g.get_xticklabels()
if plt.rcParams['text.usetex']:
    lbs = list(map(lambda l: re.sub(r"([_()])", r"\\\1", l.get_text()), lbs))
g.set_xticklabels(lbs, rotation=args.rotation)
plt.xlabel("")
plt.ylabel('\\textit{rank\_score}', fontsize=args.label_size)
plt.xticks(fontsize=args.label_size)
plt.yticks(fontsize=args.tick_size)
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}.png".format(args.png_prefix, "all_datasets"))
plt.savefig(png_file, bbox_inches='tight')
fig_idx += 1
plt.close()
print('All: '+png_file)

plt.rcParams['text.usetex'] = False
fig = plt.figure(figsize=(10, 6))
g = sns.barplot(x='Classifier', y='Neo_pep_imm count', hue='Top N', data=topN_dfs['all'], estimator=np.mean, ci=95)
lbs = g.get_xticklabels()
g.set_xticklabels(lbs, rotation=args.rotation, fontsize=args.label_size)
plt.ylabel(r'Neo_pep_imm count', size=args.label_size)
plt.xlabel("")
plt.xticks(fontsize=args.label_size)
plt.yticks(fontsize=args.tick_size)
plt.legend(fontsize=args.tick_size)
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_TopN_counts.png".format(args.png_prefix, "all"))
plt.savefig(png_file, bbox_inches='tight')
plt.close()


for group in patient_groups:
    plt.rcParams['text.usetex'] = True
    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.boxplot(x="Classifier", y=group+"_score", data=plot_df)
    sns.swarmplot(x="Classifier", y=group+"_score", data=plot_df, color=".25")
    lbs = g.get_xticklabels()
    if plt.rcParams['text.usetex']:
        lbs = list(map(lambda l: re.sub(r"([_()])", r"\\\1", l.get_text()), lbs))
    g.set_xticklabels(lbs, rotation=args.rotation)
    plt.ylabel(r'\textit{rank\_score}', fontsize=args.label_size)
    plt.xlabel("")
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
#    plt.title(patient_dict[group], fontsize=args.tick_size)
    g.figure.tight_layout()
    png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}.png".format(args.png_prefix, group))
    plt.savefig(png_file, bbox_inches='tight')
    fig_idx += 1
    plt.close()
    print(patient_group_dict[group] + ': ' + png_file)

    plt.rcParams['text.usetex'] = False
    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.barplot(x='Classifier', y='Neo_pep_imm count', hue='Top N', data=topN_dfs[group], estimator=np.mean, ci=95)
    lbs = g.get_xticklabels()
    g.set_xticklabels(lbs, rotation=args.rotation, fontsize=args.label_size)
    plt.ylabel(r'Neo_pep_imm count', size=args.label_size)
    plt.xlabel("")
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
    plt.legend(fontsize=args.tick_size)
    plt.title(patient_group_dict[group], fontsize=args.tick_size)
    g.figure.tight_layout()
    png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_TopN_counts.png".format(args.png_prefix, group))
    plt.savefig(png_file, bbox_inches='tight')
    plt.close()


plt.rcParams['text.usetex'] = False
fig = plt.figure(figsize=(10, 10))
fig.clf()
pca = PCA(n_components=2)

X = vector_df.loc[:, [c for c in vector_df.columns if c not in ['Patient', 'Peptide_id']]].to_numpy().transpose()
x_pca = pca.fit_transform(X)
classifiers = [c.split("_")[0] for c in vector_df.columns if c not in ['Patient', 'Peptide_id']]
pca_df = pd.DataFrame({'PCA_1': x_pca[:, 0], 'PCA_2': x_pca[:, 1], 'Classifier': classifiers})
variance = pca.explained_variance_ratio_

g = sns.scatterplot(data=pca_df, x='PCA_1', y='PCA_2', hue='Classifier', alpha=0.8, s=100)
plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100), size=args.label_size)
plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100), size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
plt.legend(loc="upper left", frameon=True, fontsize=args.tick_size)
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}.png".format(args.png_prefix, "clf_pca"))
plt.savefig(png_file, bbox_inches='tight')
fig_idx += 1
plt.close()

