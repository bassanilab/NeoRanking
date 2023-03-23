import argparse
import time

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-c', '--classifier', type=str, default='SVM', help='classifier to use')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_training', help='Short info for classifier run')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-nt', '--nr_train_patients', type=int, default=-1,
                    help='Number of patients in -tr option considered')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in Hyperopt')
parser.add_argument('-cv', '--nr_cv', type=int, default=5, help='Number of CV layers in RandomSearchCV')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='Response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='Mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='Immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-cat', '--cat_encoder', type=str, default='categorical', help='Convert categories to')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-nn', '--nr_negative', type=int, default=-1, help='Maximal number of non immunogenic samples')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='Genes excluded from prioritization')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')

args = parser.parse_args()

with open(DataManager().get_result_file(args.classifier, args.run_id, args.peptide_type, ext='txt', result_type='clf',
                                        suffix="train"), mode='w') as result_file:
    for arg in vars(args):
        result_file.write(f"{arg}={getattr(args, arg)}\n")
        print(f"{arg}={getattr(args, arg)}")

    normalizer = get_normalizer(args.normalizer)
    encodings = read_cat_encodings(args.patients_train[0], args.peptide_type)

    data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                             mutation_types=args.mutation_types, response_types=args.response_types,
                             immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                             max_netmhc_rank=args.max_rank_netmhc, cat_encoders=encodings)

    patients_train = \
        get_valid_patients(dataset=args.patients_train, peptide_type=args.peptide_type) \
            if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

    if args.nr_train_patients > 1:
        patients_train = patients_train[0:min(args.nr_train_patients, len(patients_train))]

    best_score_train = -np.Inf
    best_param_train = []
    best_classifier_train = None
    tot_correct_train = 0
    tot_immunogenic_train = 0
    tot_score_train = 0
    tot_negative_train = 0

    data_train, X_train, y_train = \
        data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type,
                                  nr_non_immuno_rows=args.nr_negative)

    if args.peptide_type == 'short':
        class_ratio = sum(y_train == 1)/sum(y_train == 0)
    else:
        class_ratio = None

    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader.get_categorical_idx(X_train),
                                            cat_dims=data_loader.get_categorical_dim(),
                                            input_shape=[len(args.features)], class_ratio=class_ratio)

    learner = PrioritizationLearner(args.classifier, args.scorer, optimizationParams, verbose=args.verbose,
                                    nr_iter=args.nr_iter, nr_cv=args.nr_cv, shuffle=args.shuffle)

    cvres, best_classifier, best_score, best_params = \
        learner.optimize_classifier(data_train, X_train, y_train, result_file)

    base_name = os.path.basename(result_file.name).replace("_train.txt", "")
    classifier_file = DataManager().get_classifier_file(args.classifier, base_name)

    # fit best classifier on all data
    PrioritizationLearner.save_classifier(args.classifier, best_classifier, classifier_file)
    if args.verbose > 1:
        print('Classifier = {0:s}, Scorer = {1:s}'.format(args.classifier, args.scorer))
        print('Best training params: ' + str(best_params) + ', ' + args.scorer + ': ' + str(best_score))
        print('Saved to {0:s}'.format(classifier_file))

        result_file.write('Training patients: {0}\n'.format(','.join(patients_train)))
        result_file.write('Saved to {0:s}\n'.format(classifier_file))


