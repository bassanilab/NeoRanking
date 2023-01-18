import argparse
import glob
import os
import numpy as np
import pandas as pd

from Utils.Parameters import Parameters
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--classifier_dir', type=str, default=Parameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-c', '--classifiers', type=str, nargs='+', help='classifier name')
parser.add_argument('-re', '--classifier_file_re', type=str, nargs='+', help='regexp for classifier files to use')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


def parse_clf_results(test_file, peptide_type, classifier, replicate_index, results_df):
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

        rank_col_name = "{0}_mut-seq_rank_{1}".format(classifier, replicate_index)
        if peptide_type == 'short':
            rank_col_name = "{0}_neo-pep_rank_{1}".format(classifier, replicate_index)

        df_ = pd.DataFrame()
        for row_index, row in results.iterrows():
            rank_strs = row['CD8_ranks'].split(',')
            peptide_id_strs = row['CD8_peptide_idx'].split(',')
            peptide_strs = row['CD8_mut_seqs'].split(',')
            gene_strs = row['CD8_genes'].split(',')
            nr_immuno = int(row['Nr_immunogenic'])
            max_rank = int(row['Nr_peptides'])
            dataset = get_patient_group(row['Patient'])
            ml_group = get_ml_group(row['Patient'], peptide_type)
            if dataset == "NCI":
                group_ = "{0}_{1}".format(dataset, ml_group)
            else:
                group_ = dataset

            for index in range(0, nr_immuno):
                d = {'Patient': row['Patient'], 'Dataset': group_, 'Peptide_id': peptide_id_strs[index],
                     seq_col_name: peptide_strs[index], 'Gene': gene_strs[index], 'MaxRank': max_rank,
                     rank_col_name: rank_strs[index]}
                df_ = df_.append(d, ignore_index=True)

        if results_df is None:
            results_df = df_
        else:
            results_df = results_df.merge(df_, on='Peptide_id', how='outer', suffixes=('', '_drop'))
            drop_cols = [c for c in results_df.columns if 'drop' in c]
            results_df = results_df.drop(columns=drop_cols)

    return results_df


seq_col_name = 'Mut-seq' if args.peptide_type == 'long' else 'Neo-pep'
clf_results = None
for clf, re in zip(args.classifiers, args.classifier_file_re):
    clf_result_files = glob.glob(os.path.join(args.classifier_dir, re))

    for i, file in enumerate(clf_result_files):
        clf_results = parse_clf_results(file, args.peptide_type, clf, i, clf_results)

col_order = ['Patient', 'Dataset', 'Peptide_id', seq_col_name, 'Gene', 'MaxRank']
col_order = col_order + list(np.sort([c for c in clf_results.columns if 'rank' in c]))

clf_results = clf_results[col_order]

prefix = "_".join(args.classifiers)
output_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_ranking.txt".format(prefix, args.peptide_type))
clf_results.to_csv(output_file, header=True, index=False)
