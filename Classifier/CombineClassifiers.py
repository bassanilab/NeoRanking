import argparse
from sklearn.preprocessing import *
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)


parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='prefix for patient id')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-c', '--classifiers', type=str, nargs='+', help='Classifier type')
parser.add_argument('-clf', '--classifier_files', type=str, nargs='+', help='Classifier binary files')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-r', '--max_rank', type=int, default=20, help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')

args = parser.parse_args()

normalizer = None
if args.normalizer == 'q':
    normalizer = QuantileTransformer()
elif args.normalizer == 'z':
    normalizer = StandardScaler()
elif args.normalizer == 'p':
    normalizer = PowerTransformer()

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0)

cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
cat_idx = [i for i, f in enumerate(args.features) if f in Parameters().get_categorical_features()]

optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                        cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])


patients = DataManager().get_valid_patients()
patients = [p for p in patients if any([ap == p for ap in args.patients])]

data_all, X_all, y_all = data_loader.load_patients(patients, args.input_file_tag)

cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
cat_idx = [i for i, f in enumerate(args.features) if f in Parameters().get_categorical_features()]

optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                        cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])

tags = args.classifiers
files = args.classifier_files
assert len(tags) == len(files)
classifiers = []
for tag, file in zip(tags, files):
    learner = PrioritizationLearner(tag, args.scorer, optimizationParams)
    clf = learner.load(file)
    classifiers.append(clf)


# perform leave one out on training set
tot_negative_test = 0
tot_correct_test = 0
tot_immunogenic_test = 0
tot_score_test = 0
for p in patients:
    data_test, X_test, y_test = data_loader.load_patients(p, args.input_file_tag)
    y_pred, nr_correct, nr_immuno, r, score = \
        learner.test_classifiers(classifiers, p, X_test, y_test, max_rank=args.max_rank)

    tot_negative_test += len(y_test)-nr_immuno
    tot_correct_test += nr_correct
    tot_immunogenic_test += nr_immuno
    tot_score_test += score

print('nr_patients\trun_id\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train'.format(args.max_rank))
print('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}'.format(len(args.patients), args.classifiers, tot_correct_test,
                                                tot_immunogenic_test, tot_negative_test, tot_score_test))

