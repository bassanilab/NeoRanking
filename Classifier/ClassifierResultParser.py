import pandas as pd
import numpy as np
from Utils.GlobalParameters import *


class ClassifierResults:

    def __init__(self, peptide_type, clf_result_file, index, plot_name, data, alpha=0.02, included_patients_=None):
        self.config = {}
        self.name = plot_name
        self.replicate_index = index
        self.peptide_type = peptide_type
        self.data = data
        self.datasets = set()
        self.included_patients = included_patients_
        self.results = None
        self.alpha = alpha
        self.max_rank_score = {'all': 0.0}
        self.neo_pep_imm_count = {'all': 0}
        self.parse_clf_results(clf_result_file)
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()

    def parse_clf_results(self, clf_result_file):
        with open(clf_result_file) as file_:
            lines = [line for line in file_]
            result_value_list = []
            header = None
            for line in lines:
                line = line.rstrip()
                if line.startswith("Patient\t") and header is None:
                    header = line.split("\t")
                    continue

                if header is not None and not line.startswith("Patient\t"):
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
                peptide_strs = row['CD8_mut_seqs'].split(',')
                gene_strs = row['CD8_genes'].split(',')
                nr_immuno = int(row['Nr_immunogenic'])
                ranks = []
                peptides = []
                genes = []
                for index in range(0, nr_immuno):
                    if index < len(rank_strs) and len(rank_strs[index]) > 0:
                        ranks = np.append(ranks, int(rank_strs[index]))
                        peptides = np.append(peptides, peptide_strs[index])
                        genes = np.append(genes, gene_strs[index])
                    else:
                        ranks = np.append(ranks, int(row['Nr_peptides']))
                        peptides = np.append(peptides, "")
                        genes = np.append(genes, "")

                max_rank = int(row['Nr_peptides'])
                ranks = np.array(ranks, dtype='int32')
                scores = self.calc_scores(ranks)
                dataset, ml_group = self.get_patient_groups(row['Patient'])
                if dataset == "NCI":
                    group_ = "{0}_{1}".format(dataset, ml_group)
                else:
                    group_ = dataset
                print("{0}\t{1}\t{2:.3f}".format(row['Patient'], ",".join(rank_strs), sum(scores)))
                self.datasets.add(group_)
                if group_ not in self.max_rank_score:
                    self.max_rank_score[group_] = 0.0
                if group_ not in self.neo_pep_imm_count:
                    self.neo_pep_imm_count[group_] = 0
                self.max_rank_score[group_] += sum(self.calc_scores(np.arange(1, len(ranks) + 1)))
                self.neo_pep_imm_count[group_] += len(ranks)
                self.max_rank_score['all'] += sum(self.calc_scores(np.arange(1, len(ranks)+1)))
                self.neo_pep_imm_count['all'] += len(ranks)
                df = pd.DataFrame({'Patient': row['Patient'], 'Patient_group': group_, "Rank": ranks,
                                   'Score': scores, 'Mutant_seq': peptides, 'Gene': genes,
                                   'Classifier': self.name, 'Max_rank': max_rank})

                if self.results is None:
                    self.results = df
                else:
                    self.results = pd.concat([self.results, df], ignore_index=True)

    def calc_tot_clf_result(self):
        self.sum_results['Total_score'] = sum(self.results['Score'])
        if self.datasets is not None:
            for group_ in self.datasets:
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

    def add_to_sum_df(self, rank_df):
        d = {'Classifier': self.name, 'Score': [self.sum_results['Total_score']]}
        for group_ in self.datasets:
            d[group_ + "_score"] = [self.sum_results[group_ + '_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return pd.concat([rank_df, df], ignore_index=True)

    def add_to_vector_df(self, vector_df_):
        if vector_df_ is None:
            df = self.results.loc[:, ['Patient', 'Mutant_seq', 'Score']]
            df.rename(columns={'Score': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            return df
        else:
            df = self.results.loc[:, ['Patient', 'Mutant_seq', 'Score']]
            df.rename(columns={'Score': "{0}_{1}".format(self.name, self.replicate_index)}, inplace=True)
            df['idx'] = df.groupby('Mutant_seq').cumcount()
            vector_df_['idx'] = vector_df_.groupby('Mutant_seq').cumcount()
            df = pd.merge(df, vector_df_, on=["Mutant_seq", 'idx'], how='inner', suffixes=("", "_right")).\
                drop('idx', axis=1)
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
            return pd.concat([top_n_df, df], ignore_index=True)

    def add_to_datasets(self, dataset_):
        return dataset_.union(self.datasets)

    def get_patient_groups(self, patient):
        idx = np.argwhere(self.data['patient'] == patient)[0][0]
        return self.data.iloc[idx, [self.data.columns.get_loc('dataset'), self.data.columns.get_loc('train_test')]]

    def calc_scores(self, ranks):
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

    def __init__(self, peptide_type, clf_result_file, plot_name, data, alpha=0.02, included_patients_=None):
        self.config = {}
        self.included_patients = included_patients_
        self.name = plot_name
        self.replicate_index = 0
        self.peptide_type = peptide_type
        self.datasets = set()
        self.results = None
        self.data = data
        self.alpha = alpha
        self.max_rank_score = {'all': 0.0}
        self.neo_pep_imm_count = {'all': 0}
        self.parse_clf_results(clf_result_file)
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()


class GartnerResultsNeopep(ClassifierResults):

    gartner_ranks = \
        {'3703': [2, 42, 70, 2522], '3881': [6], '3995': [2], '4007': [1], '4014': [84], '4032': [12, 32, 48],
         '4166': [9], '4242': [242], '4268': [2, 5], '4271': [15, 114], '4283': [21], '4310': [1, 8], '4323': [8],
         '4324': [24392], '4350': [2, 40, 1641], '4359': [3, 49]}

    def __init__(self, peptide_data, plot_name='Gartner et al.', alpha=0.02, included_patients_=None):
        self.name = plot_name
        self.included_patients = included_patients_
        self.peptide_data = peptide_data
        self.peptide_type = 'neopep'
        self.replicate_index = 0
        self.datasets = set(['NCI_test'])
        self.results = None
        self.alpha = alpha
        self.parse_clf_results()
        self.sum_results = pd.Series(dtype=float)
        self.group_scores = pd.Series(dtype=float)
        self.group_results = pd.DataFrame()
        self.calc_tot_clf_result()

    @staticmethod
    def get_gartner_test_patients():
        return set(GartnerResultsNeopep.gartner_ranks)

    def parse_clf_results(self):
        patients = GartnerResultsNeopep.get_gartner_test_patients().intersection(self.peptide_data['Patient'].unique())
        for p in patients:
            ranks = GartnerResultsNeopep.gartner_ranks[p]
            scores = self.calc_scores(ranks)
            max_rank = self.peptide_data.loc[self.peptide_data['Patient'] == p, 'Max_rank'].unique()[0]
            df = pd.DataFrame({'Patient': p, 'Patient_group': 'NCI_test', "Rank": ranks,
                               'Score': scores, 'Mutant_seq': 'unknown',
                               'Gene': 'unknown', 'Classifier': self.name, 'Max_rank': max_rank})

            if self.results is None:
                self.results = df
            else:
                self.results = pd.concat([self.results, df], ignore_index=True)

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
        for group_ in self.datasets:
            d[group_ + "_score"] = [self.sum_results[group_ + '_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return pd.concat([df, rank_df], ignore_index=True)


class GartnerResultsMutation(ClassifierResults):

    gartner_res_file = GlobalParameters.gartner_nmer_rank_file
    gartner_ranks = pd.read_csv(gartner_res_file, sep=',', header=0, index_col=False)
    gartner_ranks = gartner_ranks.astype({'Patient': str})

    def __init__(self, peptide_data, plot_name='Gartner et al.', alpha=0.02, included_patients_=None):
        self.name = plot_name
        self.included_patients = included_patients_
        self.datasets = set(['NCI_test'])
        self.peptide_data = peptide_data
        self.peptide_type = 'mutation'
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
        return set(GartnerResultsMutation.gartner_ranks['Patient'])

    def parse_clf_results(self):
        patients = GartnerResultsMutation.get_gartner_test_patients().intersection(self.peptide_data['Patient'].unique())
        for p in patients:

            ranks = GartnerResultsMutation.gartner_ranks.loc[GartnerResultsMutation.gartner_ranks['Patient'] == p,
                                                             'Rank Nmer Model'].to_numpy()
            scores = self.calc_scores(ranks)
            max_rank = self.peptide_data.loc[self.peptide_data['Patient'] == p, 'Max_rank'].unique()[0]
            df = pd.DataFrame({'Patient': p, 'Patient_group': 'NCI_test', "Rank": ranks,
                               'Score': scores, 'Mutant_seq': 'unknown',
                               'Gene': 'unknown', 'Classifier': self.name, 'Max_rank': max_rank})

            if self.results is None:
                self.results = df
            else:
                self.results = pd.concat([self.results, df], ignore_index=True)

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
        for ds in self.datasets:
            d[ds + "_score"] = [self.sum_results[ds + '_score']]
        df = pd.DataFrame(d)
        if rank_df is None:
            return df
        else:
            return pd.concat([df, rank_df], ignore_index=True)

