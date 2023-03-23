import argparse
import re

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--classifier_dir', type=str, default=Parameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-c1', '--classifier1_file_re', type=str, nargs='+', help='classifier files to use')
parser.add_argument('-c2', '--classifier2_file_re', type=str, nargs='+', help='classifier files to use')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-te', '--dataset_test', type=str, help='dataset ids for testing')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for training')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')
parser.add_argument('-cat', '--cat_encoder', type=str, default='float', help='convert categories to numbers')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='genes excluded from prioritization')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')
parser.add_argument('-nn', '--nr_negative', type=int, default=-1, help='Maximal number of non immunogenic neo-peptides')

args = parser.parse_args()

patients_test = \
    get_valid_patients(dataset=args.dataset_test, peptide_type=args.peptide_type) \
        if args.dataset_test and len(args.dataset_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)

encodings = read_cat_encodings(args.dataset_test, args.peptide_type)

mgr = DataManager()
patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type)))


normalizer = get_normalizer(args.normalizer)

if args.dataset_test == 'NCI_train':
    response_types = ['CD8', 'CD4/CD8', 'negative']
else:
    response_types = ['CD8', 'CD4/CD8', 'negative', 'not_tested']
data_loader_test = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                              mutation_types=args.mutation_types, response_types=response_types,
                              immunogenic=args.immunogenic, min_nr_immuno=1, cat_type=args.cat_encoder,
                              max_netmhc_rank=20, cat_encoders=encodings, excluded_genes=args.excluded_genes)

classifier1_files = []
for wc in args.classifier1_file_re:
    classifier1_files = classifier1_files + glob.glob(os.path.join(args.classifier_dir, wc))

classifier2_files = []
for wc in args.classifier2_file_re:
    classifier2_files = classifier2_files + glob.glob(os.path.join(args.classifier_dir, wc))


def get_learner(classifier_name, x):
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader_test.get_categorical_idx(x),
                                            cat_dims=data_loader_test.get_categorical_dim(),
                                            input_shape=[len(args.features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


for w in np.arange(0.0, 1.1, 0.1):
    vc_result_file = os.path.join(args.classifier_dir, 'Voting_classifier_{0:.2f}_test.txt'.format(w))
    open(vc_result_file, mode='w').close()

voting_clfs = []

data_all, X_all, y_all = \
    data_loader_test.load_patients(patients_test, args.input_file_tag, args.peptide_type, verbose=False)

for p in patients_test:
    data_test, X_test, y_test = data_loader_test.load_patients(p, args.input_file_tag, args.peptide_type, verbose=True)
    if data_test is None:
        continue

    idx = data_all['patient'] != p
    data_train = data_all.loc[idx, :]
    X_train = X_all.loc[idx, :]
    y_train = y_all[idx]
    if args.nr_negative > 0:
        data_train, X_train, y_train = data_loader_test.sample_rows(data_train, X_train, y_train, args.nr_negative)

    for w in np.arange(0.0, 1.1, 0.1):
        weights = []
        for clf in classifier1_files:
            clf_name = os.path.basename(clf).split("_")[0]
            learner = get_learner(clf_name, X_test)
            classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
            voting_clfs.append((clf_name, classifier))
            weights.append(1.0-w)

        for clf in classifier2_files:
            clf_name = os.path.basename(clf).split("_")[0]
            learner = get_learner(clf_name, X_test)
            classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
            if clf_name == 'CatBoost':
                classifier = learner.get_optimization_params().get_classifier(clf_name, clf.get_params())
            classifier = learner.fit_classifier(X_train, y_train, classifier)
            voting_clfs.append((clf_name, classifier))
            weights.append(w)

        vc_result_file = os.path.join(args.classifier_dir, 'Voting_classifier_{0:.2f}_test.txt'.format(w))
        with open(vc_result_file, mode='a') as result_file:
            y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
                learner.test_voting_classifier(voting_clfs, weights, p, data_test, X_test, y_test, max_rank=args.max_rank,
                                               report_file=result_file)

