import argparse
import time

import numpy as np
import itertools
from scipy.stats import fisher_exact

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-tr', '--patients_train', type=str, default="NCI_train", help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, default="test", help='patient ids for test set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-fp', '--file_prefix', type=str, default='TESLA_ranking', help='Prefix for output file name')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-nt', '--nr_thresholds', type=int, default=30, help='Number of thresholds per feature')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='Response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='Mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='Immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-nn', '--nr_negative', type=int, default=-1, help='Maximal number of non immunogenic samples')

args = parser.parse_args()


def get_quantile_thresholds(x, nr_thresholds=100):
    d = {}
    for f in x.columns:
        if f in Parameters().get_numerical_features():
            d[f] = np.quantile(x[f], np.arange(0.0, 1.0, 1.0/nr_thresholds))
        else:
            d[f] = x[f].unique()

    return d


def find_optimal_thresholds(x, y, nr_thresholds=100):
    thresholds = get_quantile_thresholds(x, nr_thresholds)

    best_pv = 2.0
    best_thresh_comb = ()
    nr_comb = 1
    for f in thresholds:
        nr_comb *= len(thresholds[f])

    j = 0
    for combination in itertools.product(*thresholds.values()):
        if j % 1000 == 0:
            print("{0} of {1} combinations tested. best p-value = {2} ...".format(j, nr_comb, best_pv))
        j += 1
        idx = np.full(len(y), True)
        for i, f in enumerate(thresholds.keys()):
            order_rel = Parameters().get_order_relation(f)
            if order_rel == '>':
                idx = np.logical_and(idx, x[f] >= combination[i])
            else:
                idx = np.logical_and(idx, x[f] < combination[i])

        table = np.array([[sum(y[~idx] == 0), sum(y[idx] == 0)], [sum(y[~idx] == 1), sum(y[idx] == 1)]])
        oddsr, pv = fisher_exact(table, alternative='two-sided')
        if pv < best_pv:
            best_thresh_comb = combination
            best_pv = pv

    return dict(zip(thresholds.keys(), best_thresh_comb)), best_pv


def filter_with_thresholds(x, y, thresholds):

    idx = np.full(len(y), True)
    for i, f in enumerate(thresholds.keys()):
        order_rel = Parameters().get_order_relation(f)
        if order_rel == '>':
            idx = np.logical_and(idx, x[f] >= thresholds[f])
        else:
            idx = np.logical_and(idx, x[f] < thresholds[f])

    return x.loc[idx, :], y[idx]


normalizer = get_normalizer('q')
train_file = "{0}_{1}_{2}_train.txt".format(args.file_prefix, args.patients_train, args.patients_test)
with open(os.path.join(Parameters().get_pickle_dir(), train_file), mode='w') as result_file:
    for arg in vars(args):
        result_file.write(f"{arg}={getattr(args, arg)}\n")
        print(f"{arg}={getattr(args, arg)}")

    encodings = read_cat_encodings(args.patients_train, args.peptide_type)

    data_loader_train = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                                   mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                                   immunogenic=args.immunogenic, min_nr_immuno=0, cat_type='int', cat_encoders=encodings)

    patients_train = \
        get_valid_patients(dataset=args.patients_train, peptide_type=args.peptide_type) \
            if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

    data_train, X_train, y_train = \
        data_loader_train.load_patients(patients_train, args.input_file_tag, args.peptide_type,
                                        nr_non_immuno_rows=args.nr_negative)

    optimal_thresholds, p_value = find_optimal_thresholds(X_train, y_train, args.nr_thresholds)

    print('Classifier = Wells_threshold, best p-value = {0:.5e}'.format(p_value))
    print('Best training threshold rule: ' + str(optimal_thresholds))

    result_file.write('Training patients: {0}\n'.format(','.join(patients_train)))
    result_file.write('Classifier = Wells_threshold, best p-value = {0:.5e}'.format(p_value))
    for f in optimal_thresholds:
        result_file.write('{0}:{1}'.format(f, optimal_thresholds[f]))


train_file = "{0}_{1}_{2}_test.txt".format(args.file_prefix, args.patients_train, args.patients_test)
with open(os.path.join(Parameters().get_pickle_dir(), train_file), mode='w') as result_file:
    result_file.write("Patient\tSelected count\tImmunogenic count\n")
    patients_test = \
        get_valid_patients(dataset=args.patients_test, peptide_type=args.peptide_type) \
            if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)

    data_loader_test = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                                  mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                                  immunogenic=args.immunogenic, min_nr_immuno=1, cat_type='int',
                                  max_netmhc_rank=20, cat_encoders=encodings)

    for p in patients_test:
        data_test, X_test, y_test = data_loader_test.load_patients(p, args.input_file_tag, args.peptide_type, verbose=True)
        if data_test is None:
            continue

        X, y = filter_with_thresholds(X_test, y_test, optimal_thresholds)

        result_file.write("{0}\t{1}\t{2}\n".format(p, len(y), sum(y)))
