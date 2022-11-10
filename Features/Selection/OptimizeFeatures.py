import argparse
import warnings

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *
from Features.Selection.FeatureOptimizer import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-c', '--classifier_tag', type=str, default='SVM', help='classifier to use')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_training', help='Short info for classifier run')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in RandomSearchCV')
parser.add_argument('-nt', '--max_nr_trials', type=int, default=100, help='Number of trials in feature optimization')
parser.add_argument('-nc', '--nr_classifiers', type=int, default=1,
                    help='Number of best classifiers included for voting')
parser.add_argument('-cv', '--nr_cv', type=int, default=5, help='Number of CV layers in RandomSearchCV')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-e', '--nr_epoch', type=int, default=150, help='Number of epochs for DNN training')
parser.add_argument('-ep', '--early_stopping_patience', type=int, default=150,
                    help='Patience for early stopping for DNN training')
parser.add_argument('-b', '--batch_size', type=int, default=150, help='Batch size for DNN training')
parser.add_argument('-cat', '--cat_to_num', dest='cat_to_num', action='store_true',
                    help='convert categories to numbers')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-ct', '--combine_test', dest='combine_test', action='store_true', help='Combine test data')
parser.add_argument('-nn', '--nr_negative', type=int, default=-1, help='Maximal number of non immunogenic samples')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='genes excluded from prioritization')

args = parser.parse_args()

with open(DataManager().get_result_file(args.classifier_tag, args.run_id, args.peptide_type), mode='w') \
        as result_file:
    for arg in vars(args):
        result_file.write(f"{arg}={getattr(args, arg)}\n")
        print(f"{arg}={getattr(args, arg)}")

    normalizer = get_normalizer(args.normalizer)

    data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                             mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                             immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num,
                             max_netmhc_rank=10000)

    patients_train = \
        get_valid_patients(dataset=args.patients_train, peptide_type=args.peptide_type) \
            if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

    data_train, X_train, y_train = data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type,
                                                             nr_non_immuno_rows=args.nr_negative)

    warnings.filterwarnings("ignore")
    featureOptimizer = FeatureOptimizer(args.classifier_tag, X_train, y_train, max_nr_models=3, scorer=args.scorer,
                                        cat_dims=data_loader.get_categorical_dim(), peptide_type=args.peptide_type,
                                        alpha=args.alpha, nr_iter=args.nr_iter, nr_cv=args.nr_cv, shuffle=args.shuffle,
                                        report_file=result_file)

    featureOptimizer.run(args.max_nr_trials)
    featureOptimizer.report_models()
