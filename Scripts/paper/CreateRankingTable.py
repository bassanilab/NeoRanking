import argparse

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import patches
import seaborn as sns
from sklearn.decomposition import PCA
import re

from DataWrangling.DataTransformer import DataTransformer
from Utils.Util_fct import *
from Utils.GlobalParameters import *

parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
parser.add_argument('-d', '--data_dir', type=str, help='Directory containing clf results')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-g', '--patient_group_plot_names', type=str, default='', help='Patient groups plot names')
parser.add_argument('-cln', '--classifier_plot_names', type=str, default='', help='Classifier names in plots')
parser.add_argument('-pn', '--patient_plot_names', type=str, default='', help='Patient names in plots')
parser.add_argument('-o', '--plot_order', type=str, help='Order of classifier plot names')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class ClassifierResults:

    def __init__(self, base_name, index, plot_name, alpha=0.02, peptide_type='short'):
        self.config = {}
        self.name = plot_name
        self.hyperopt_results = {}
        if 'Voting_classifier' not in base_name:
            self.parse_config(base_name)
        self.replicate_index = index
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
                                   'Classifier': self.name, 'Max_rank': max_rank, 'Replicate': self.replicate_index})

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

    def add_to_results(self, ranking_results_df):
        if ranking_results_df is None:
            ranking_results_df = self.results.copy()
        else:
            ranking_results_df = ranking_results_df.append(self.results, ignore_index=True)
        return ranking_results_df

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


parse_neodisc = True
parse_gartner = True
parse_simple = True

plot_df = None
vector_df = None
topN_dfs = None
hyperopt_scores = None
patient_groups = set()

plt_name_dict = ast.literal_eval(args.classifier_plot_names)

max_rank_score = None
tot_imm_count = None
if args.plot_order:
    plot_order = args.plot_order.split(',')
else:
    plot_order = None

results_df = None
for j, regexp in enumerate(args.clf_result_files_re):
    clf_result_file = glob.glob(os.path.join(args.data_dir, regexp))

    for i, result_file in enumerate(clf_result_file):
        if os.path.getsize(result_file) > 0:
            print(result_file)
            base_file_name = re.sub("_test.txt$", "", result_file)
            clf_results = ClassifierResults(base_file_name, i, plt_name_dict[j], alpha=0.005)
            results_df = clf_results.add_to_results(results_df)


patient_group_dict = ast.literal_eval(args.patient_group_plot_names)
patient_dict = ast.literal_eval(args.patient_plot_names)

results_df = results_df.drop(columns=['Score', 'Peptide_id'])

results_df['Patient'] = results_df['Patient'].replace(patient_dict)
results_df['Patient_group'] = results_df['Patient_group'].replace(patient_group_dict)
if args.peptide_type == 'short':
    results_df = results_df.rename(columns={'Patient_group': 'Dataset', 'Mutant_seq': 'Neo-peptide'})
else:
    results_df = results_df.rename(columns={'Patient_group': 'Dataset', 'Mutant_seq': 'Mutation'})

result_file = os.path.join(GlobalParameters().get_plot_dir(), "ML_Classifier_ranking_{0}.txt".format(args.peptide_type))

results_df.to_csv(result_file, sep="\t", index=False, header=True)

