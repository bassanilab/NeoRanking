import argparse
import re

import numpy as np

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Features.MLPrediction.MLPrediction import MLPrediction
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-dl', '--classifier_dir_long', type=str, default=Parameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-cl1', '--classifier1_long_re', type=str, nargs='+', help='classifier files to use for mut-seq')
parser.add_argument('-cl2', '--classifier2_long_re', type=str, nargs='+', help='classifier files to use for mut-seq')
parser.add_argument('-ds', '--classifier_dir_short', type=str, default=Parameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-cs1', '--classifier1_short_re', type=str, nargs='+', help='classifier files to use for neo-pep')
parser.add_argument('-cs2', '--classifier2_short_re', type=str, nargs='+', help='classifier files to use for neo-pep')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, default='NCI-train', help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, default='test', help='patient ids for test set')
parser.add_argument('-il', '--input_file_tag_long', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc long input file (patient)_(input_file_tag).txt')
parser.add_argument('-is', '--input_file_tag_short', type=str, default='stab_chop',
                    help='File tag for neodisc short input file (patient)_(input_file_tag).txt')
parser.add_argument('-fl', '--features_long', type=str, nargs='+', help='Mut-seq features used by classifier')
parser.add_argument('-fs', '--features_short', type=str, nargs='+', help='Neo-pep features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')
parser.add_argument('-cat', '--cat_encoder', type=str, default='float', help='convert categories to numbers')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='genes excluded from prioritization')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')
parser.add_argument('-mrl', '--max_rank_long', type=int, default=50,
                    help='If nr mutations smaller than mrl, no filtering is performed for mutations')
parser.add_argument('-mrs', '--max_rank_short', type=int, default=3, help='Maximal rank of short peptide for a mutation')
parser.add_argument('-rl', '--keep_long_ratio', type=float, default=0.5, help='The ratio of long peptides to keep')

args = parser.parse_args()

patients_test = \
    get_valid_patients(dataset=args.patients_test, peptide_type='short') \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type='short')

mgr = DataManager()
patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients('short')))

normalizer = get_normalizer(args.normalizer)
response_types = ['CD8', 'CD4/CD8', 'negative', 'not_tested']

encodings = read_cat_encodings(args.patients_train, 'long')
data_loader_long = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_long,
                              mutation_types=args.mutation_types, response_types=response_types,
                              immunogenic=args.immunogenic, min_nr_immuno=1, cat_type=args.cat_encoder,
                              cat_encoders=encodings, excluded_genes=args.excluded_genes)

encodings = read_cat_encodings(args.patients_train, 'short')
data_loader_short = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_short,
                               mutation_types=args.mutation_types, response_types=response_types,
                               immunogenic=args.immunogenic, min_nr_immuno=1, cat_type=args.cat_encoder,
                               max_netmhc_rank=20, cat_encoders=encodings, excluded_genes=args.excluded_genes)

classifier1_files_long = []
for wc in args.classifier1_long_re:
    classifier1_files_long = classifier1_files_long + glob.glob(os.path.join(args.classifier_dir_long, wc))

classifier2_files_long = []
if args.classifier2_long_re is not None:
    for wc in args.classifier2_long_re:
        classifier2_files_long = classifier2_files_long + glob.glob(os.path.join(args.classifier_dir_long, wc))

classifier1_files_short = []
for wc in args.classifier1_short_re:
    classifier1_files_short = classifier1_files_short + glob.glob(os.path.join(args.classifier_dir_short, wc))

classifier2_files_short = []
for wc in args.classifier2_short_re:
    classifier2_files_short = classifier2_files_short + glob.glob(os.path.join(args.classifier_dir_short, wc))


def convert_peptide_id_long(row):
    fields = row['peptide_id'].split('|')
    return fields[0]+":"+fields[2]


def add_mutant_id(data_long_):
    mutant_id_long = data_long_.apply(convert_peptide_id_long, axis=1)
    data_long_.loc[:, 'mutant_id'] = mutant_id_long

    return data_long_


def get_learner(classifier_name, data_loader, features, x):
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader.get_categorical_idx(x),
                                            cat_dims=data_loader.get_categorical_dim(),
                                            input_shape=[len(features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


def keep_short_for_best_long(data_l, data_s, x_s, y_s, max_rank_long, max_rank_short, long_perc=0.5):
    nr_long = min(data_l.shape[0], max_rank_long)
    nr_keep = data_l.shape[0] * long_perc
    if nr_keep > nr_long:
        nr_long = round(nr_keep)

    mutant_ids = data_l.loc[data_l['mutation_rank'] <= nr_long, 'mutant_id']

    idx = np.array(data_s.apply(lambda row: any(mutant_ids == row['mut_seqid']), axis=1))
    y_filtered = np.array([y_s[i] for i in range(len(y_s)) if idx[i]])

    data_s = data_s[idx]
    x_s = x_s[idx]

    idx = data_s['rank_in_mutation'] <= max_rank_short
    return data_s[idx], x_s[idx], y_filtered[idx]


result_file_name = os.path.join(Parameters().get_pickle_dir(), "Voting_classifier_filter_{0:.2f}_{1}_{2}_test.txt".
                                format(args.keep_long_ratio, args.max_rank_long, args.max_rank_short))
open(result_file_name, mode='w').close()

with open(result_file_name, mode='a') as result_file:

    for p in patients_test:
        data_long, X_long, y_long = data_loader_long.load_patients(p, args.input_file_tag_long, 'long', verbose=True)
        if data_long is None:
            continue

        voting_clfs = []
        for clf in classifier1_files_long:
            clf_name = os.path.basename(clf).split("_")[0]
            learner = get_learner(clf_name, data_loader_long, args.features_long, X_long)
            classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
            voting_clfs.append((clf_name, classifier, 1.0))

        for clf in classifier2_files_long:
            clf_name = os.path.basename(clf).split("_")[0]
            learner = get_learner(clf_name, data_loader_long, args.features_long, X_long)
            classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
            voting_clfs.append((clf_name, classifier, 1.0))

        ml_prediction = MLPrediction(voting_clfs, data_loader_long)
        data_long = ml_prediction.add_features_long(data_long, X_long)
        data_long = add_mutant_id(data_long)

        data_short, X_short, y_short = data_loader_short.load_patients(p, args.input_file_tag_short, 'short', verbose=True)
        if data_short is None:
            continue

        voting_clfs = []
        for clf in classifier1_files_short:
            clf_name = os.path.basename(clf).split("_")[0]
            learner = get_learner(clf_name, data_loader_short, args.features_short, X_short)
            classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
            voting_clfs.append((clf_name, classifier, 1.0))

        for clf in classifier2_files_short:
            clf_name = os.path.basename(clf).split("_")[0]
            learner = get_learner(clf_name, data_loader_short, args.features_short, X_short)
            classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
            voting_clfs.append((clf_name, classifier, 1.0))

        ml_prediction = MLPrediction(voting_clfs, data_loader_short)
        data_short = ml_prediction.add_features_short(data_short, X_short)

        learner_long = get_learner(clf_name, data_loader_long, args.features_long, X_long)
        data_short, X_short, y_short = \
            learner_long.add_long_prediction_to_short(data_long, data_short, X_short, y_short)

        data_short['mutant_score'] = -data_short['mutant_rank']
        data_short = data_short.sort_values(by=['peptide_score', 'mutant_score'], ascending=False)
        data_short = data_short.reset_index()

        print("=====> Patient {0}, mut_cnt={1}, pep_cnt={2}".format(p, data_long.shape[0], data_short.shape[0]))
        print(data_short.head(min(100, data_short.shape[0]))[['mutant_seq', 'peptide_score', 'rank_in_mutation',
                                                             'mutation_rank', 'mutation_score']].to_string())

        print("-------------------------------------------------------------------------------------------------------")
        print(data_short.loc[data_short['response_type'] == 'CD8',
                             ['mutant_seq', 'peptide_score', 'rank_in_mutation', 'mutation_rank', 'mutation_score']]
              .to_string())

        max_imm_mut_rank = max(data_short.loc[data_short['response_type'] == 'CD8', 'mutation_rank'])
        print("mut_cnt={0}, max_imm_mut_rank={1}, mut_cnt_ratio={2:.2f}".
              format(data_long.shape[0], max_imm_mut_rank, max_imm_mut_rank/data_long.shape[0]))
        print("-------------------------------------------------------------------------------------------------------")

        data_short, X_short, y_short = \
            keep_short_for_best_long(data_long, data_short, X_short, y_short, args.max_rank_long, args.max_rank_short,
                                     args.keep_long_ratio)
        data_short = data_short.reset_index()
        print("-------------------------------------------------------------------------------------------------------")
        print(data_short.loc[data_short['response_type'] == 'CD8',
                             ['mutant_seq', 'peptide_score', 'rank_in_mutation', 'mutation_rank', 'mutation_score']]
              .to_string())

        data_short['in_mutant_score'] = -data_short['rank_in_mutation']
        data_short = data_short.sort_values(by=['peptide_score', 'in_mutant_score', 'mutation_score'], ascending=False)
        data_short = data_short.reset_index(drop=True)
        print("-------------------------------------------------------------------------------------------------------")
        print(data_short.loc[data_short['response_type'] == 'CD8',
                             ['mutant_seq', 'peptide_score', 'rank_in_mutation', 'mutation_rank', 'mutation_score']]
              .to_string())

        data_short['in_mutant_score'] = -data_short['rank_in_mutation']
        data_short = data_short.sort_values(by=['mutation_score', 'peptide_score'], ascending=False)
        data_short = data_short.reset_index(drop=True)
        print("-------------------------------------------------------------------------------------------------------")
        print(data_short.loc[data_short['response_type'] == 'CD8',
                             ['mutant_seq', 'peptide_score', 'rank_in_mutation', 'mutation_rank', 'mutation_score']]
              .to_string())

