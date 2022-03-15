import argparse
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from sklearn.preprocessing import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-cfl', '--classifier_long', type=str, default='', help='classifier file for long peptides')
parser.add_argument('-cfs', '--classifier_short', type=str, default='', help='classifier file for short peptides')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_short_peptide_test', help='Short info for classifier run')
parser.add_argument('-fl', '--features_long', type=str, nargs='+', help='Features used by classifier for long peptides')
parser.add_argument('-fs', '--features_short', type=str, nargs='+', help='Features used by classifier for short peptides')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-ct', '--combine_test', dest='combine_test', action='store_true', help='Combine test data')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0, cat_to_num=args.cat_to_num,
                         max_netmhc_rank=10000)

cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
cat_idx = [i for i, f in enumerate(args.features) if f in Parameters().get_categorical_features()]

optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                        cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])

with open(DataManager().get_classifier_result_file(args.classifier_long, 'long'), mode='r') as file:
    classifier_long_tag = os.path.basename(args.classifier_long).split('_')[0]
    classifier_long = \
        PrioritizationLearner.load_classifier(classifier_long_tag, optimizationParams, args.classifier_long)
    file.write('Classifier for long peptides imported from file {0}\n'.format(args.classifier_long))
    learner = PrioritizationLearner(classifier_long_tag, args.scorer, optimizationParams, verbose=args.verbose)

with open(DataManager().get_classifier_result_file(args.classifier_short, 'short'), mode='r') as file:
    classifier_short_tag = os.path.basename(args.classifier_short).split('_')[0]
    classifier_short = \
        PrioritizationLearner.load_classifier(classifier_short_tag, optimizationParams, args.classifier_short)
    file.write('Classifier for short peptides imported from file {0}\n'.format(args.classifier_short))


patients_test = \
    get_valid_patients(patients=args.patients_test, peptide_type=args.peptide_type) \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)
mgr = DataManager()
patients_test = patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type))

tot_negative_test = 0
tot_correct_test = 0
tot_immunogenic_test = 0
tot_score_test = 0
data_loader_long = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_long,
                              mutation_types=args.mutation_types, response_types=args.response_types,
                              immunogenic=args.immunogenic, min_nr_immono=0, cat_to_num=args.cat_to_num,
                              max_netmhc_rank=10000)


if patients_test is not None:
    if not args.combine_test:
        for p in patients_test:
            data_test, X_test, y_test = data_loader.load_patients(p, args.input_file_tag, args.peptide_type)
            y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                learner.test_classifier(classifier_long, p, X_test.to_numpy(), y_test, max_rank=args.max_rank)

            tot_negative_test += len(y_test) - nr_immuno
            tot_correct_test += nr_correct
            tot_immunogenic_test += nr_immuno
            tot_score_test += score
    else:
        data_test, X_test, y_test = data_loader.load_patients(patients_test, args.input_file_tag, args.peptide_type)
        y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
            learner.test_classifier(best_classifier_train, ','.join(patients_test), X_test.to_numpy(), y_test,
                                    max_rank=args.max_rank)

        tot_negative_test += len(y_test) - nr_immuno
        tot_correct_test += nr_correct
        tot_immunogenic_test += nr_immuno
        tot_score_test += score

if args.verbose > 0:
    print('nr_patients\trun_id\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train'.format(args.max_rank))
    print('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}'.format(len(patients_test), classifier_file, tot_correct_test,
                                                    tot_immunogenic_test, tot_negative_test, tot_score_test))
