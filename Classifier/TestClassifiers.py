import argparse
import re

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-c', '--classifier_file_re', type=str, nargs='+', help='classifier files to use')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')
parser.add_argument('-cat', '--cat_encoder', type=str, default='float', help='convert categories to numbers')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='genes excluded from prioritization')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')

args = parser.parse_args()

patients_test = \
    get_valid_patients(patients=args.patients_test, peptide_type=args.peptide_type) \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)

mgr = DataManager()
patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type)))

normalizer = get_normalizer(args.normalizer)
encodings = read_cat_encodings(args.patients_train[0], args.peptide_type)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=1, cat_type=args.cat_encoder,
                         max_netmhc_rank=20, cat_encoders=encodings)

classifier_files = []
for wc in args.classifier_file_re:
    classifier_files = classifier_files + glob.glob(os.path.join(Parameters().get_pickle_dir(), wc))


def get_learner(classifier_name, x):
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader.get_categorical_idx(x),
                                            cat_dims=data_loader.get_categorical_dim(),
                                            input_shape=[len(args.features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


for p in patients_test:
    data_test, X_test, y_test = data_loader.load_patients(p, args.input_file_tag, args.peptide_type, verbose=True)
    if data_test is None:
        continue

    for clf in classifier_files:
        clf_name = os.path.basename(clf).split("_")[0]
        learner = get_learner(clf_name, X_test)
        classifier = learner.load_classifier(clf_name, learner.get_optimization_params(), clf)
        result_file_name = re.sub("\\.\\w+$", "_test.txt", clf)

        with open(result_file_name, mode='a') as result_file:
            y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
                learner.test_classifier(classifier, p, data_test, X_test, y_test, max_rank=args.max_rank,
                                        report_file=result_file)

