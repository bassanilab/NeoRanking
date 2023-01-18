import argparse

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import patches
import seaborn as sns
from sklearn.decomposition import PCA
import re

from DataWrangling.DataLoader import DataLoader
from Utils.Util_fct import *
from Utils.Parameters import *

parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
parser.add_argument('-d', '--data_dir', type=str, help='Directory containing clf results')
parser.add_argument('-fp', '--file_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="svg", help='File type for plot (png, svg or pdf')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-g', '--patient_group_plot_names', type=str, default='', help='Patient groups plot names')
parser.add_argument('-cln', '--classifier_plot_names', type=str, default='', help='Classifier order in plots')
parser.add_argument('-o', '--plot_order', type=str, help='Order of classifier plot names')
parser.add_argument('-a', '--alpha', type=float, default=0.02, help='Coefficient alpha in score function')
parser.add_argument('-nd', '--neodisc', dest='neodisc', action='store_true', help='Include neodisc comparison')
parser.add_argument('-ga', '--gartner', dest='gartner', action='store_true', help='Include Gartner et al. comparison')
parser.add_argument('-sr', '--simple_ranking', dest='simple_ranking', action='store_true', help='Include simple ranking comparison')
parser.add_argument('-rot', '--rotation', type=float, default=30.0, help='x-axis label rotation')
parser.add_argument('-las', '--label_size', type=float, default=25.0, help='Axis label size')
parser.add_argument('-xl', '--xlabel', type=str, default="", help='x-label')
parser.add_argument('-tis', '--tick_size', type=float, default=20.0, help='Axis tick size')
parser.add_argument('-tts', '--title_size', type=float, default=20.0, help='title size')
parser.add_argument('-res', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=6.00, help='Figure height in inches')
parser.add_argument('-les', '--legend_size', type=float, default=15, help='Legend size in float')
parser.add_argument('-hy', '--hyperopt', dest='hyperopt', action='store_true', help='Include hyperopt training score')
parser.add_argument('-oc', '--one_color', type=str, default=None, help='all boxes in same color')
parser.add_argument('-bar', '--bar_plot', dest='bar_plot', action='store_true', help='bars instead of boxes')
parser.add_argument('-ttp', '--title_prefix', type=str, default='', help='prefix for plot title')
parser.add_argument('-ylim', '--rank_score_lim', type=str, default='', help='plot limits for rankscore')
parser.add_argument('-cm', '--color_map', type=str, default='', help='color map for classifiers')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class ClassifierResults:

    def __init__(self, base_name, index, plot_name, alpha=0.02, included_patients_=None, peptide_type='short'):
        self.config = {}
        self.name = plot_name
        self.hyperopt_results = {}
        if 'Voting_classifier' not in base_name:
            self.parse_config(base_name)
        self.replicate_index = index
        self.peptide_type = peptide_type
        self.groups = set()
        self.included_patients = included_patients_
        self.results = None
        self.alpha = alpha
        self.max_rank_score = {'all': 0.0}
        self.neo_pep_imm_count = {'all': 0}
        self.parse_clf_results(base_name)
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()

    def parse_config(self, base_name):
        train_file = base_name+"_train.txt"
        with open(train_file) as file_:
            lines = [line for line in file_]
            for line in lines:
                if line.startswith("Hyperopt"):
                    break
                fields = line.split("=")
                if len(fields) > 1:
                    self.config[fields[0]] = fields[1]

            self.parse_hyperopt_results(lines)

    def parse_hyperopt_results(self, lines):
        for line in lines:
            if line.startswith("Patient"):
                break
            if line.startswith("Hyperopt"):
                fields = line.replace("Hyperopt: ", "").split("; ")
                for f in fields:
                    k, v = f.split('=')
                    if k == 'Params':
                        self.hyperopt_results[k] = ast.literal_eval(v)
                    else:
                        self.hyperopt_results[k] = float(v)

    def parse_clf_results(self, base_name):
        test_file = base_name+"_test.txt"
        with open(test_file) as file_:
            lines = [line for line in file_]
            result_value_list = []
            header = None
            for line in lines:
                line = line.rstrip()
                if line.startswith("Patient"):
                    header = line.split("\t")
                    continue

                if header is not None:
                    values = np.array(line.split("\t"))
                    result_value_list.append(pd.Series(values))

            results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
            results.columns = header

            for row_index, row in results.iterrows():
                if self.included_patients and row['Patient'] not in self.included_patients:
                    continue
                if type(row['CD8_ranks']) is not str:
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
                for index in range(0, nr_immuno):
                    if index < len(rank_strs) and len(rank_strs[index]) > 0:
                        ranks = np.append(ranks, int(rank_strs[index]))
                        peptide_ids = np.append(peptide_ids, peptide_id_strs[index])
                        peptides = np.append(peptides, peptide_strs[index])
                        genes = np.append(genes, gene_strs[index])
                    else:
                        ranks = np.append(ranks, int(row['Nr_peptides']))
                        peptide_ids = np.append(peptide_ids, "")
                        peptides = np.append(peptides, "")
                        genes = np.append(genes, "")

                max_rank = int(row['Nr_peptides'])
                ranks = np.array(ranks, dtype='int32')
                scores = self.calc_scores(ranks)
                dataset = get_patient_group(row['Patient'])
                ml_group = get_ml_group(row['Patient'], self.peptide_type)
                if dataset == "NCI":
                    group_ = "{0}_{1}".format(dataset, ml_group)
                else:
                    group_ = dataset
                print("{0}\t{1}\t{2:.3f}".format(row['Patient'], ",".join(rank_strs), sum(scores)))
                self.groups.add(group_)
                if group_ not in self.max_rank_score:
                    self.max_rank_score[group_] = 0.0
                if group_ not in self.neo_pep_imm_count:
                    self.neo_pep_imm_count[group_] = 0
                self.max_rank_score[group_] += sum(self.calc_scores(np.arange(1, len(ranks) + 1)))
                self.neo_pep_imm_count[group_] += len(ranks)
                self.max_rank_score['all'] += sum(self.calc_scores(np.arange(1, len(ranks)+1)))
                self.neo_pep_imm_count['all'] += len(ranks)
                df = pd.DataFrame({'Patient': row['Patient'], 'Patient_group': group_, "Rank": ranks,
                                   'Score': scores, 'Peptide_id': peptide_ids, 'Mutant_seq': peptides, 'Gene': genes,
                                   'Classifier': self.name, 'Max_rank': max_rank})

                if self.results is None:
                    self.results = df
                else:
                    self.results = self.results.append(df, ignore_index=True)

    def calc_tot_clf_result(self):
        self.sum_results['Total_score'] = sum(self.results['Score'])
        if self.groups is not None:
            for group_ in self.groups:
                self.sum_results[group_ + "_score"] = \
                    sum(self.results.loc[self.results['Patient_group'] == group_, 'Score'])

    def get_name(self):
        return self.name

    def get_config(self):
        return self.config

    def get_results_data(self):
        return self.results

    def get_max_score(self):
        return self.max_rank_score

    def get_imm_count(self):
        return self.neo_pep_imm_count

    def get_hyperopt_param(self, param):
        return self.hyperopt_results[param]

    def add_to_sum_df(self, rank_df):
        d = {'Classifier': self.name, 'Score': [self.sum_results['Total_score']]}
        for group_ in self.groups:
            d[group_ + "_score"] = [self.sum_results[group_ + '_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return rank_df.append(df, ignore_index=True)

    def add_to_vector_df(self, vector_df_):
        if vector_df_ is None:
            df = self.results.loc[:, ['Patient', 'Peptide_id', 'Score']]
            df.rename(columns={'Score': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            return df
        else:
            df = self.results.loc[:, ['Patient', 'Peptide_id', 'Score']]
            df.rename(columns={'Score': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            df = vector_df_.merge(df, how='outer', left_on='Peptide_id', right_on='Peptide_id', suffixes=("", "_right"))
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

    def add_to_hyperopt_df(self, hyperopt_df):
        if 'Score' not in self.hyperopt_results:
            return hyperopt_df

        if hyperopt_df is None:
            return pd.DataFrame({'Classifier': self.name, 'rank_score': self.hyperopt_results['Score'],
                                 'running_time': self.hyperopt_results['Time'], 'Replicate': self.replicate_index},
                                index=[0])
        else:
            df = pd.DataFrame({'Classifier': self.name, 'rank_score': self.hyperopt_results['Score'],
                               'running_time': self.hyperopt_results['Time'], 'Replicate': self.replicate_index},
                              index=[0])
            return hyperopt_df.append(df, ignore_index=True)

    def add_to_patient_groups(self, patient_groups_):
        return patient_groups_.union(self.groups)

    def calc_scores(self, ranks):
        # return [np.exp(np.multiply(-self.alpha, r-1)) if r < 200 else 0 for r in ranks]
        return [np.exp(np.multiply(-self.alpha, r-1)) for r in ranks]

    def get_topN_count(self, max_rank, patient_group):
        if patient_group == 'all':
            return self.results.apply(lambda row: row['Rank'] <= max_rank, axis=1).sum()
        else:
            if any(self.results['Patient_group'] == patient_group):
                return self.results.loc[self.results['Patient_group'] == patient_group, :]\
                    .apply(lambda row: row['Rank'] <= max_rank, axis=1).sum()
            else:
                return 0


class SimpleRankingResults(ClassifierResults):

    def __init__(self, base_name, plot_name, alpha=0.02, peptide_type='short'):
        self.config = {}
        self.included_patients = None
        self.name = plot_name
        self.hyperopt_results = {}
        self.replicate_index = 0
        self.peptide_type = peptide_type
        self.groups = set()
        self.results = None
        self.alpha = alpha
        self.max_rank_score = {'all': 0.0}
        self.neo_pep_imm_count = {'all': 0}
        self.parse_clf_results(base_name)
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()


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
                data, X, y = data_loader.load_patients(patient, 'rt', args.peptide_type, verbose=True)

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


class GartnerResultsShort(ClassifierResults):

    gartner_ranks = \
        {'3703': [2, 42, 70, 2522], '3881': [6], '3995': [2], '4007': [1], '4014': [84], '4032': [12, 32, 48],
         '4166': [9], '4242': [242], '4268': [2, 5], '4271': [15, 114], '4283': [21], '4310': [1, 8], '4323': [8],
         '4324': [24392], '4350': [2, 40, 1641], '4359': [3, 49]}

    def __init__(self, peptide_data, plot_name='Gartner et al.\nranking', alpha=0.02):
        self.name = plot_name
        self.groups = set(['NCI_test'])
        self.peptide_data = peptide_data
        self.peptide_type = 'short'
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
        return set(GartnerResultsShort.gartner_ranks)

    def parse_clf_results(self):
        patients = GartnerResultsShort.get_gartner_test_patients().intersection(self.peptide_data['Patient'].unique())
        for p in patients:
            ranks = GartnerResultsShort.gartner_ranks[p]
            scores = self.calc_scores(ranks)
            max_rank = self.peptide_data.loc[self.peptide_data['Patient'] == p, 'Max_rank'].unique()[0]
            df = pd.DataFrame({'Patient': p, 'Patient_group': 'NCI_test', "Rank": ranks,
                               'Score': scores, 'Peptide_id': 'unknown', 'Mutant_seq': 'unknown',
                               'Gene': 'unknown', 'Classifier': self.name, 'Max_rank': max_rank})

            if self.results is None:
                self.results = df
            else:
                self.results = self.results.append(df, ignore_index=True)

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
            return pd.concat([df, top_n_df], ignore_index=True)

    def add_to_sum_df(self, rank_df):
        d = {'Classifier': self.name, 'Score': [self.sum_results['Total_score']]}
        for group_ in self.groups:
            d[group_ + "_score"] = [self.sum_results[group_ + '_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return pd.concat([df, rank_df], ignore_index=True)


class GartnerResultsLong(ClassifierResults):

    gartner_res_file = Parameters().get_gartner_long_result_file()
    gartner_ranks = pd.read_csv(gartner_res_file, sep=',', header=0, index_col=False)
    gartner_ranks = gartner_ranks.astype({'Patient': str})

    def __init__(self, peptide_data, plot_name='Gartner et al.\nranking', alpha=0.02):
        self.name = plot_name
        self.groups = set(['NCI_test'])
        self.peptide_data = peptide_data
        self.peptide_type = 'long'
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
        return set(GartnerResultsLong.gartner_ranks['Patient'])

    def parse_clf_results(self):
        patients = GartnerResultsLong.get_gartner_test_patients().intersection(self.peptide_data['Patient'].unique())
        for p in patients:

            ranks = GartnerResultsLong.gartner_ranks.loc[GartnerResultsLong.gartner_ranks['Patient'] == p,
                                                         'Rank Nmer Model'].to_numpy()
            scores = self.calc_scores(ranks)
            max_rank = self.peptide_data.loc[self.peptide_data['Patient'] == p, 'Max_rank'].unique()[0]
            df = pd.DataFrame({'Patient': p, 'Patient_group': 'NCI_test', "Rank": ranks,
                               'Score': scores, 'Peptide_id': 'unknown', 'Mutant_seq': 'unknown',
                               'Gene': 'unknown', 'Classifier': self.name, 'Max_rank': max_rank})

            if self.results is None:
                self.results = df
            else:
                self.results = self.results.append(df, ignore_index=True)

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
            return pd.concat([df, top_n_df], ignore_index=True)

    def add_to_sum_df(self, rank_df):
        d = {'Classifier': self.name, 'Score': [self.sum_results['Total_score']]}
        for group_ in self.groups:
            d[group_ + "_score"] = [self.sum_results[group_ + '_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return pd.concat([df, rank_df], ignore_index=True)


included_patients = None
if args.gartner:
    if args.peptide_type == 'short':
        included_patients = GartnerResultsShort.get_gartner_test_patients()
    else:
        included_patients = GartnerResultsShort.get_gartner_test_patients()

parse_neodisc = True
parse_gartner = True
parse_simple = True

plot_df = None
vector_df = None
topN_dfs = None
hyperopt_scores = None
patient_groups = set()

color_map = None
ylim_dict = {}
if args.rank_score_lim != '':
    ylim_dict = ast.literal_eval(args.rank_score_lim)
elif args.color_map != '':
    cm = ast.literal_eval(args.color_map)
    for key in cm:
        color_map[key] = sns.color_palette()[cm[key]]

plt_name_dict = ast.literal_eval(args.classifier_plot_names)

topN_palette = sns.color_palette("Greens_d", 3)

max_rank_score = None
tot_imm_count = None
if args.plot_order:
    plot_order = args.plot_order.split(',')
else:
    plot_order = None

for j, regexp in enumerate(args.clf_result_files_re):
    clf_result_file = glob.glob(os.path.join(args.data_dir, regexp))

    for i, result_file in enumerate(clf_result_file):
        if os.path.getsize(result_file) > 0:
            print(result_file)
            base_file_name = re.sub("_test.txt$", "", result_file)
            if 'SimpleRanking' in base_file_name:
                if args.simple_ranking:
                    simple_results = SimpleRankingResults(base_file_name, plt_name_dict[j], alpha=args.alpha)
                    plot_df = simple_results.add_to_sum_df(plot_df)
                    vector_df = simple_results.add_to_vector_df(vector_df)
                    for ds in topN_dfs:
                        topN_dfs[ds] = simple_results.add_to_topN_df(topN_dfs[ds], ds)
            else:
                clf_results = ClassifierResults(base_file_name, i, plt_name_dict[j], alpha=args.alpha,
                                                included_patients_=included_patients)
                if max_rank_score is None:
                    max_rank_score = clf_results.get_max_score()
                if tot_imm_count is None:
                    tot_imm_count = clf_results.get_imm_count()
                if topN_dfs is None:
                    topN_dfs = {'all': None}
                    for g in clf_results.groups:
                        topN_dfs[g] = None
                plot_df = clf_results.add_to_sum_df(plot_df)
                vector_df = clf_results.add_to_vector_df(vector_df)
                patient_groups = clf_results.add_to_patient_groups(patient_groups)
                for ds in topN_dfs:
                    topN_dfs[ds] = clf_results.add_to_topN_df(topN_dfs[ds], ds)
                if args.hyperopt:
                    hyperopt_scores = clf_results.add_to_hyperopt_df(hyperopt_scores)

            if args.neodisc and parse_neodisc:
                neodisc_results = NeoDiscResults(clf_results.results, alpha=args.alpha)
                plot_df = neodisc_results.add_to_sum_df(plot_df)
                vector_df = neodisc_results.add_to_vector_df(vector_df)
                for ds in topN_dfs:
                    topN_dfs[ds] = neodisc_results.add_to_topN_df(topN_dfs[ds], ds)
                parse_neodisc = False

            if args.gartner and parse_gartner:
                if args.peptide_type == 'short':
                    gartner_results = GartnerResultsShort(clf_results.results, alpha=args.alpha)
                    plot_df = gartner_results.add_to_sum_df(plot_df)
                    for ds in topN_dfs:
                        topN_dfs[ds] = gartner_results.add_to_topN_df(topN_dfs[ds], ds)
                    parse_gartner = False
                else:
                    gartner_results = GartnerResultsLong(clf_results.results, alpha=args.alpha)
                    plot_df = gartner_results.add_to_sum_df(plot_df)
                    for ds in topN_dfs:
                        topN_dfs[ds] = gartner_results.add_to_topN_df(topN_dfs[ds], ds)
                    parse_gartner = False


plot_df = plot_df.astype({'Classifier': str, 'Score': float})
for group in patient_groups:
    plot_df.astype({group+"_score": float})

patient_group_dict = ast.literal_eval(args.patient_group_plot_names)

if args.rotation > 0:
    ha = 'center'
else:
    ha = 'center'

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
if args.bar_plot:
    g = sns.barplot(x="Classifier", y="Score", data=plot_df, order=plot_order, color=args.one_color, estimator=np.mean,
                    errorbar=('ci', 95), palette=color_map)
else:
    g = sns.boxplot(x="Classifier", y="Score", data=plot_df, order=plot_order, color=args.one_color, palette=color_map)
sns.swarmplot(x="Classifier", y="Score", data=plot_df, color=".25", order=plot_order)
lbs = g.get_xticklabels()
g.set_xticklabels(lbs, rotation=args.rotation, ha=ha)
if 'all' in ylim_dict:
    g.set(ylim=(ylim_dict['all'], None))
plt.xlabel("")
plt.ylabel('rank_score', fontsize=args.label_size, style='italic')
plt.xticks(fontsize=args.label_size)
plt.yticks(fontsize=args.tick_size)
g.set_title("{0} Maximal rank_score = {1:.3f}".format(args.title_prefix, max_rank_score['all']),
            fontsize=args.title_size)
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}.{2}".format(args.file_prefix, "all_datasets",
                                                                          args.file_type))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()
print('All: '+png_file)

if args.hyperopt:
    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    if args.bar_plot:
        g = sns.barplot(x="Classifier", y="rank_score", data=hyperopt_scores, order=plot_order, color=args.one_color,
                        estimator=np.mean, errorbar=('ci', 95), palette=color_map)
    else:
        g = sns.boxplot(x="Classifier", y="rank_score", data=hyperopt_scores, order=plot_order, color=args.one_color)
    sns.swarmplot(x="Classifier", y="rank_score", data=hyperopt_scores, color=".25", palette=color_map)
    lbs = g.get_xticklabels()
    g.set_xticklabels(lbs, rotation=args.rotation, ha=ha)
    if 'hyperopt' in ylim_dict:
        g.set(ylim=(ylim_dict['hyperopt'], None))
    plt.xlabel(args.xlabel, fontsize=args.label_size)
    plt.ylabel('rank_score', fontsize=args.label_size, style='italic')
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
    png_file = os.path.join(Parameters().get_plot_dir(),
                            "{0}_{1}.{2}".format(args.file_prefix, "hyperopt_score", args.file_type))
    plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()
    print('Hyperopt score: '+png_file)

    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    if args.bar_plot:
        g = sns.barplot(x="Classifier", y="running_time", data=hyperopt_scores, order=plot_order, color=args.one_color,
                        estimator=np.mean, errorbar=('ci', 95), palette=color_map)
    else:
        g = sns.boxplot(x="Classifier", y="running_time", data=hyperopt_scores, order=plot_order, color=args.one_color,
                        palette=color_map)
    sns.swarmplot(x="Classifier", y="running_time", data=hyperopt_scores, color=".25")
    lbs = g.get_xticklabels()
    g.set_xticklabels(lbs, rotation=args.rotation, ha=ha)
    plt.xlabel(args.xlabel, fontsize=args.label_size)
    plt.ylabel('Running time (sec)', fontsize=args.label_size)
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
    png_file = os.path.join(Parameters().get_plot_dir(),
                            "{0}_{1}.{2}".format(args.file_prefix, "hyperopt_time", args.file_type))
    plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()
    print('Hyperopt time: '+png_file)

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.barplot(x='Classifier', y='Neo_pep_imm count', hue='Top N', data=topN_dfs['all'], estimator=np.mean,
                errorbar=('ci', 95), order=plot_order, palette=topN_palette)
lbs = g.get_xticklabels()
g.set_xticklabels(lbs, rotation=args.rotation, fontsize=args.label_size, ha=ha)
peptide_label = 'Neo-pep_imm count' if args.peptide_type == 'short' else 'Mut-seq_imm count'
plt.ylabel(peptide_label, size=args.label_size)
plt.xlabel(args.xlabel, fontsize=args.label_size)
plt.xticks(fontsize=args.label_size)
plt.yticks(fontsize=args.tick_size)
plt.ylim(0, tot_imm_count['all']+3)
plt.axhline(y=tot_imm_count['all'], color="red", linestyle="--", linewidth=2)
handles = [patches.Patch(color=topN_palette[0], label='Top 20'),
           patches.Patch(color=topN_palette[1], label='Top 50'),
           patches.Patch(color=topN_palette[2], label='Top 100')
           ]
sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=3, handles=handles, title="All test sets",
                frameon=False, fontsize=args.legend_size)
png_file = os.path.join(Parameters().get_plot_dir(),
                        "{0}_{1}_TopN_counts.{2}".format(args.file_prefix, "all", args.file_type))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()

txt_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_TopN_counts.txt".format(args.file_prefix, group))
df_agg = topN_dfs['all'].groupby(['Classifier', 'Top N']).agg({'Neo_pep_imm count': 'mean'})
pd.DataFrame(df_agg).to_csv(txt_file, sep='\t', header=True, index=True)
with open(txt_file, 'a') as file:
    file.write("\nTotal count: {0}\n".format(tot_imm_count['all']))


for group in patient_groups:
    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    if args.bar_plot:
        g = sns.barplot(x="Classifier", y=group+"_score", data=plot_df, order=plot_order, color=args.one_color,
                        estimator=np.mean, errorbar=('ci', 95), palette=color_map)
    else:
        g = sns.boxplot(x="Classifier", y=group+"_score", data=plot_df, order=plot_order, color=args.one_color,
                        palette=color_map)
    sns.swarmplot(x="Classifier", y=group+"_score", data=plot_df, color=".25", order=plot_order)
    lbs = g.get_xticklabels()
    g.set_xticklabels(lbs, rotation=args.rotation, ha=ha)
    if group in ylim_dict:
        g.set(ylim=(ylim_dict[group], None))
    plt.ylabel('rank_score', fontsize=args.label_size, style='italic')
    plt.xlabel(args.xlabel, fontsize=args.label_size)
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
#    plt.title(patient_dict[group], fontsize=args.tick_size)
    g.set_title("{0}{1}: maximal rank_score = {2:.3f}".
                format(args.title_prefix, patient_group_dict[group], max_rank_score[group]), fontsize=args.title_size)
    png_file = os.path.join(Parameters().get_plot_dir(),
                            "{0}_{1}.{2}".format(args.file_prefix, group, args.file_type))
    plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()
    print(patient_group_dict[group] + ': ' + png_file)

    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    g = sns.barplot(x='Classifier', y='Neo_pep_imm count', hue='Top N', data=topN_dfs[group], estimator=np.mean,
                    errorbar=('ci', 95), order=plot_order, palette=topN_palette)
    lbs = g.get_xticklabels()
    g.set_xticklabels(lbs, rotation=args.rotation, fontsize=args.label_size, ha=ha)
    peptide_label = 'neo-pep' if args.peptide_type == 'short' else 'mut-seq'
    plt.ylabel("{0}_imm count".format(peptide_label), size=args.label_size)
    plt.xlabel(args.xlabel, fontsize=args.label_size)
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
    plt.ylim(0, tot_imm_count[group]+3)
    plt.axhline(y=tot_imm_count[group], color="red", linestyle="--", linewidth=2)
    handles = [patches.Patch(color=topN_palette[0], label='Top 20'),
               patches.Patch(color=topN_palette[1], label='Top 50'),
               patches.Patch(color=topN_palette[2], label='Top 100')
               ]
    sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.20), ncol=3, handles=handles,
                    frameon=False, fontsize=args.legend_size, title=patient_group_dict[group],
                    title_fontsize=args.legend_size)
    png_file = os.path.join(Parameters().get_plot_dir(),
                            "{0}_{1}_TopN_counts.{2}".format(args.file_prefix, group, args.file_type))
    plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()

    txt_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_TopN_counts.txt".format(args.file_prefix, group))
    df_agg = topN_dfs[group].groupby(['Classifier', 'Top N']).agg({'Neo_pep_imm count': 'mean'})
    pd.DataFrame(df_agg).to_csv(txt_file, sep='\t', header=True, index=True)
    with open(txt_file, 'a') as file:
        file.write("\nTotal count: {0}\n".format(tot_imm_count[group]))


fig = plt.figure()
fig.set_figheight(max(args.figure_height, args.figure_width))
fig.set_figwidth(max(args.figure_height, args.figure_width))
pca = PCA(n_components=2)

X = vector_df.loc[:, [c for c in vector_df.columns if c not in ['Patient', 'Peptide_id']]].to_numpy().transpose()
x_pca = pca.fit_transform(X)
classifiers = [c.rsplit("_", 1)[0] for c in vector_df.columns if c not in ['Patient', 'Peptide_id']]
pca_df = pd.DataFrame({'PCA_1': x_pca[:, 0], 'PCA_2': x_pca[:, 1], 'Classifier': classifiers})
variance = pca.explained_variance_ratio_

g = sns.scatterplot(data=pca_df, x='PCA_1', y='PCA_2', hue='Classifier', alpha=0.8, s=100, palette=color_map)
plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100), size=args.label_size)
plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100), size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
plt.legend(loc="best", frameon=True, fontsize=args.tick_size)
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}.{2}".format(args.file_prefix, "clf_pca", args.file_type))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()

