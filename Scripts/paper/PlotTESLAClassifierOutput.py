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

parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
parser.add_argument('-d', '--data_dir', type=str, help='Directory containing clf results')
parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-tf', '--tesla_result_file', type=str, help='PNG output files prefix')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-g', '--patient_group_plot_names', type=str, default='', help='Patient groups plot names')
parser.add_argument('-c', '--classifier_plot_names', type=str, default='', help='Classifier plot names')
parser.add_argument('-rot', '--rotation', type=float, default=30.0, help='x-axis label rotation')
parser.add_argument('-las', '--label_size', type=float, default=25.0, help='Axis label size')
parser.add_argument('-tis', '--tick_size', type=float, default=15.0, help='Axis tick size')
parser.add_argument('-tts', '--title_size', type=float, default=20.0, help='title size')
parser.add_argument('-res', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=6.00, help='Figure height in inches')
parser.add_argument('-les', '--legend_size', type=float, default=15, help='Legend size in float')
parser.add_argument('-ttp', '--title-prefix', type=str, default='', help='title prefix')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class ClassifierResults:

    def __init__(self, base_name, index, plot_name, peptide_type='short'):
        self.name = plot_name
        self.replicate_index = index
        self.peptide_type = peptide_type
        self.results = None
        self.parse_clf_results(base_name)

    def parse_clf_results(self, base_name):
        test_file = base_name+"_test.txt"
        with open(test_file) as file:
            lines = [line for line in file]
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
                dataset = get_patient_group(row['Patient'])
                ml_group = get_ml_group(row['Patient'], self.peptide_type)
                if dataset == "NCI":
                    group = "{0}_{1}".format(dataset, ml_group)
                else:
                    group = dataset
                df = pd.DataFrame({'Patient': row['Patient'], 'Patient_group': group, "Rank": ranks,
                                   'Peptide_id': peptide_ids, 'Mutant_seq': peptides, 'Gene': genes,
                                   'Classifier': self.name, 'Max_rank': max_rank})

                if self.results is None:
                    self.results = df
                else:
                    self.results = self.results.append(df, ignore_index=True)

    def get_name(self):
        return self.name

    def get_results_data(self):
        return self.results

    def add_to_thresh_clf_res(self, thresh_clf_res):
        patients = []
        dataset = []
        imm_cnts = []
        selected_cnts = []
        for row_index, row in thresh_clf_res.loc[thresh_clf_res['Method'] == 'Threshold Optimization'].iterrows():
            patients.append(row['Patient'])
            dataset.append(row['Dataset'])
            ranks = self.results.loc[self.results['Patient'] == row['Patient'], 'Rank'].to_numpy()
            imm_cnts.append(sum(ranks < row['Selected count']))
            if row['Immunogenic count'] == 0:
                selected_cnts.append(1)
            else:
                selected_cnts.append(ranks[row['Immunogenic count']-1])

        df = pd.DataFrame({'Patient': patients, 'Dataset': dataset, 'Selected count': selected_cnts,
                           'Immunogenic count': imm_cnts, 'Method': self.name})
        return thresh_clf_res.append(df, ignore_index=True)


def get_patient_grp_ml(row, peptide_type):
    dataset = get_patient_group(row['Patient'])
    if dataset == "NCI":
        ml_group = get_ml_group(row['Patient'], peptide_type)
        group = "{0}_{1}".format(dataset, ml_group)
    else:
        group = dataset
    return group


df_tesla = pd.read_csv(args.tesla_result_file, sep='\t', header=0)
df_tesla['Dataset'] = df_tesla.apply(lambda r: get_patient_grp_ml(r, args.peptide_type), axis=1)
df_tesla['Method'] = 'Threshold Optimization'
df_tesla = df_tesla.loc[df_tesla['Immunogenic count'] > 0, :]

plot_df = None
topN_dfs = None
plt_name_dict = ast.literal_eval(args.classifier_plot_names)
max_rank_score = None
tot_imm_count = None
for j, regexp in enumerate(args.clf_result_files_re):
    clf_result_files = glob.glob(os.path.join(args.data_dir, regexp))

    for i, result_file in enumerate(clf_result_files):
        if os.path.getsize(result_file) > 0:
            print(result_file)
            base_file_name = re.sub("_test.txt$", "", result_file)
            clf_results = ClassifierResults(base_file_name, i, plt_name_dict[j])
            df_tesla = clf_results.add_to_thresh_clf_res(df_tesla)

patients = df_tesla.loc[df_tesla['Method'] == 'Threshold Optimization', 'Patient']
values = df_tesla.loc[df_tesla['Method'] == 'Threshold Optimization', 'Selected count']
patients, values = zip(*sorted(zip(patients, values), key=lambda d: d[1], reverse=True))

cnt_thresh = df_tesla.groupby(['Patient', 'Method']).agg({'Selected count': 'mean'})
ratio = []
for p in patients:
    ratio.append(cnt_thresh.loc[(p, 'Threshold Optimization'), 'Selected count'] /
                 cnt_thresh.loc[(p, 'LR'), 'Selected count'])
print('Mean size ratio {0}'.format(np.mean(ratio)))

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.pointplot(data=df_tesla, x="Patient", y="Selected count", hue='Method', order=patients)
lbs = g.get_xticklabels()
g.set_xticklabels(lbs, rotation=90, fontsize=args.tick_size)
g.set(yscale="log")
plt.xlabel("")
plt.ylabel('Neo-peptides', fontsize=args.label_size, style='italic')
png_file = os.path.join(Parameters().get_plot_dir(), "{0}.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()
print('All: '+png_file)

