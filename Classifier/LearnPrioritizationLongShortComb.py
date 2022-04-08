import argparse
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from sklearn.preprocessing import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-cfl', '--classifier_long', type=str, default='', help='classifier file for long peptides')
parser.add_argument('-cfs', '--classifier_short', type=str, default='', help='classifier file for short peptides')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-il', '--input_file_tag_long', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file for long peptides (patient)_long_(input_file_tag).txt')
parser.add_argument('-is', '--input_file_tag_short', type=str, default='stab_chop',
                    help='File tag for neodisc input file for short peptides (patient)_short_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_short_peptide_test', help='Short info for classifier run')
parser.add_argument('-fl', '--features_long', type=str, nargs='+', help='Features used by classifier for long peptides')
parser.add_argument('-fs', '--features_short', type=str, nargs='+',
                    help='Features used by classifier for short peptides')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-cat', '--cat_to_num', dest='cat_to_num', action='store_true',
                    help='convert categories to numbers')
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in RandomSearchCV')
parser.add_argument('-cv', '--nr_cv', type=int, default=5, help='Number of CV layers in RandomSearchCV')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-ms', '--max_nr_short', type=int, default=3, help='Maximum nr short peptide scores used')
parser.add_argument('-manl', '--max_nr_long', type=int, default=100,
                    help='Maximal number of ranked long peptides to consider for short peptide ranking')
parser.add_argument('-minl', '--min_nr_long', type=int, default=50,
                    help='If nr mutations smaller than minl, no filtering is performed for mutations')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

response_types = ['CD8', 'CD4/CD8', 'negative']
data_loader_long = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_long,
                              mutation_types=args.mutation_types, response_types=response_types,
                              immunogenic=args.immunogenic, min_nr_immono=0, cat_to_num=args.cat_to_num,
                              max_netmhc_rank=10000)

data_loader_short = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_short,
                               mutation_types=args.mutation_types, response_types=response_types,
                               immunogenic=args.immunogenic, min_nr_immono=0, cat_to_num=args.cat_to_num,
                               max_netmhc_rank=10000)

cat_features = [f for f in args.features_long if f in Parameters().get_categorical_features()]
cat_idx = [i for i, f in enumerate(args.features_long) if f in Parameters().get_categorical_features()]

optimizationParams = \
    OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                       cat_dims=data_loader_long.get_categorical_dim(), input_shape=[len(args.features_long)])

with open(args.classifier_long, mode='r') as result_file:
    classifier_long_tag = os.path.basename(args.classifier_long).split('_')[0]
    classifier_long = \
        PrioritizationLearner.load_classifier(classifier_long_tag, optimizationParams, args.classifier_long)
    learner_long = PrioritizationLearner(classifier_long_tag, args.scorer, optimizationParams, verbose=args.verbose)

with open(args.classifier_short, mode='r') as result_file:
    classifier_short_tag = os.path.basename(args.classifier_short).split('_')[0]
    classifier_short = \
        PrioritizationLearner.load_classifier(classifier_short_tag, optimizationParams, args.classifier_short)
    learner_short = PrioritizationLearner(classifier_short_tag, args.scorer, optimizationParams, verbose=args.verbose)

patients_train_long = \
    get_valid_patients(patients=args.patients_train, peptide_type='long') \
        if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type='long')
patients_train_short = \
    get_valid_patients(patients=args.patients_train, peptide_type='short') \
        if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type='short')
patients_test_long = \
    get_valid_patients(patients=args.patients_test, peptide_type='long') \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type='long')
patients_test_short = \
    get_valid_patients(patients=args.patients_test, peptide_type='short') \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type='short')

mgr = DataManager()
patients_train_long = patients_train_long.intersection(mgr.get_immunogenic_patients('long'))
patients_train_short = patients_train_short.intersection(mgr.get_immunogenic_patients('short'))
patients_train = sorted(patients_train_long.intersection(patients_train_short))
patients_test_long = patients_test_long.intersection(mgr.get_immunogenic_patients('long'))
patients_test_short = patients_test_short.intersection(mgr.get_immunogenic_patients('short'))
patients_test = sorted(patients_test_long.intersection(patients_test_short))

tot_negative_test = 0
tot_correct_test = 0
tot_immunogenic_test = 0
tot_score_test = 0


def train_classifier(data, x, y):
    class_ratio = sum(y == 1) / sum(y == 0)

    optimizationParams = \
        OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                           cat_dims=data_loader_short.get_categorical_dim(), input_shape=[x.shape[1]],
                           class_ratio=class_ratio)

    learner = \
        PrioritizationLearner(classifier_short_tag, args.scorer, optimizationParams, verbose=args.verbose,
                              nr_iter=args.nr_iter, nr_cv=args.nr_cv, shuffle=args.shuffle)

    return learner.optimize_classifier(x.to_numpy(), y)


def keep_short_for_best_long(data_l, data_s, x_s, y_s, factor=0.05):
    min_nr = args.min_nr_long
    if data_l.shape[0] * factor < min_nr:
        nr_long = min(data_l.shape[0], min_nr)
    else:
        nr_long = min(args.max_nr_long, round((data_l.shape[0]) * factor))

    mutant_ids = data_l.sort_values(by='mutation_score', axis=0, ascending=False).\
                     iloc[0:nr_long, np.argwhere(data_l.columns == 'mutant_id')[0][0]]

    idx = np.array(data_s.apply(lambda row: any(mutant_ids == row['mut_seqid']), axis=1))
    y_filtered = np.array([y_s[i] for i in range(len(y_s)) if idx[i]])

    return data_s[idx], x_s[idx], y_filtered


clf_tag = "{0}_{1}".format(classifier_long_tag, classifier_short_tag)
with open(DataManager().get_result_file(clf_tag, args.run_id), mode='w') as result_file:
    for arg in vars(args):
        result_file.write(f"{arg}={getattr(args, arg)}\n")

    if patients_train is not None:
        data_long, X_long, y_long = \
            data_loader_long.load_patients(patients_train, args.input_file_tag_long, 'long', verbose=True)
        data_short, X_short, y_short = \
            data_loader_short.load_patients(patients_train, args.input_file_tag_short, 'short', verbose=True)
        # First iteration ---------------------------------------------------------------------------------
        # add predictions of classifier_long on long peptides to dataframe of short peptides
        data_long_comb, data_short_comb, X_short_comb, y_short_comb = \
            learner_long.add_long_prediction_to_short(classifier_long, data_long, X_long, data_short, X_short, y_short,
                                                      normalizer)
        # train short peptide classifier with long peptide predictions
        cvres, best_classifier_short_1, best_score_short, best_params_short = \
            train_classifier(data_short_comb, X_short_comb, y_short_comb)

        # add predictions of best_classifier_short on short peptides to dataframe of long peptides
        data_long_comb, X_long_comb, y_long_comb, data_pred_comb = \
            learner_short.add_short_prediction_to_long(classifier_short, data_short, X_short, args.max_nr_short,
                                                       data_long_comb, X_long, y_long, normalizer=None)
        # train long peptide classifier with short peptide predictions
        cvres, best_classifier_long_1, best_score_long, best_params_long = \
            train_classifier(data_long_comb, X_long_comb, y_long_comb)

        # # Second iteration ---------------------------------------------------------------------------------
        # # add predictions of classifier_long on long peptides to dataframe of short peptides
        # data_long_comb, data_short_comb, X_short_comb, y_short_comb = \
        #     learner_long.add_long_prediction_to_short(best_classifier_long_1, data_long_comb, X_long_comb,
        #                                               data_short, X_short, y_short, normalizer)
        # # train short peptide classifier with long peptide predictions
        # cvres, best_classifier_short_2, best_score_short, best_params_short = \
        #     train_classifier(data_short_comb, X_short_comb, y_short_comb)
        #
        # # add predictions of best_classifier_short on short peptides to dataframe of long peptides
        # data_long_comb, X_long_comb, y_long_comb, data_pred_comb = \
        #     learner_short.add_short_prediction_to_long(best_classifier_short_2, data_short_comb, X_short_comb,
        #                                                args.max_nr_short, data_long, X_long, y_long, normalizer=None)
        # # train long peptide classifier with short peptide predictions
        # cvres, best_classifier_long_2, best_score_long, best_params_long = \
        #     train_classifier(data_long_comb, X_long_comb, y_long_comb)

    response_types = ['not_tested', 'CD8', 'CD4/CD8', 'negative']
    data_loader_long = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_long,
                                  mutation_types=args.mutation_types, response_types=response_types,
                                  immunogenic=args.immunogenic, min_nr_immono=0, cat_to_num=args.cat_to_num,
                                  max_netmhc_rank=10000)

    data_loader_short = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_short,
                                   mutation_types=args.mutation_types, response_types=response_types,
                                   immunogenic=args.immunogenic, min_nr_immono=0, cat_to_num=args.cat_to_num,
                                   max_netmhc_rank=10000)

    if patients_test is not None:
        for p in patients_test:
            data_long, X_long, y_long = data_loader_long.load_patients(p, args.input_file_tag_long, 'long')
            data_short, X_short, y_short = \
                data_loader_short.load_patients(p, args.input_file_tag_short, 'short', verbose=False)

            # add long prediction probabilities to short dataframe
            data_long_comb, data_short_comb, X_short_comb, y_short_comb = \
                learner_long.add_long_prediction_to_short(classifier_long, data_long, X_long, data_short, X_short,
                                                          y_short, normalizer)

            # add short prediction probabilities to long dataframe
            data_long_comb, X_long_comb, y_long_comb, data_pred_comb = \
                learner_short.add_short_prediction_to_long(classifier_short, data_short, X_short, args.max_nr_short,
                                                           data_long_comb, X_long, y_long, normalizer=None)

            # add long prediction probabilities of retrained classifier to short dataframe
            data_long_comb, data_short_comb, X_short_comb, y_short_comb = \
                learner_long.add_long_prediction_to_short(best_classifier_long_1, data_long_comb, X_long_comb,
                                                          data_short_comb, X_short_comb, y_short_comb, normalizer)

            # # add short prediction probabilities of retrained classifier to long dataframe
            # data_long_comb, X_long_comb, y_long_comb, data_pred_comb = \
            #     learner_short.add_short_prediction_to_long(best_classifier_short_1, data_short_comb, X_short_comb,
            #                                                args.max_nr_short, data_long, X_long, y_long,
            #                                                normalizer=None)
            #
            # # add long prediction probabilities of retrained classifier to short dataframe
            # data_long_comb, data_short_comb, X_short_comb, y_short_comb = \
            #     learner_long.add_long_prediction_to_short(best_classifier_long_2, data_long_comb, X_long_comb,
            #                                               data_short_comb, X_short_comb, y_short_comb, normalizer)
            #
            # # # add short prediction probabilities of retrained classifier to long dataframe
            # data_long_comb, X_long_comb, y_long_comb, data_pred_comb = \
            #     learner_short.add_short_prediction_to_long(best_classifier_short_2, data_short_comb, X_short_comb,
            #                                                args.max_nr_short, data_long, X_long, y_long, normalizer=None)
            #
            data_short_comb, X_short_comb, y_short_comb = \
                keep_short_for_best_long(data_long_comb, data_short_comb, X_short_comb, y_short_comb, factor=0.05)

            y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                learner_short.test_classifier(best_classifier_short_1, p, X_short_comb.to_numpy(), y_short_comb,
                                              max_rank=args.max_rank, report_file=result_file)

            tot_negative_test += len(y_short) - nr_immuno
            tot_correct_test += nr_correct
            tot_immunogenic_test += nr_immuno
            tot_score_test += score

    if args.verbose > 0:
        print('nr_patients\trun_id\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train'.format(args.max_rank))
        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}'.
              format(len(patients_test), args.run_id, tot_correct_test, tot_immunogenic_test,
                     tot_negative_test, tot_score_test))
    result_file.write('nr_patients\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train\n'.
                      format(args.max_rank))
    result_file.write('{0}\t{1}\t{2}\t{3}\t{4:.3f}\n'.
                      format(len(patients_test), tot_correct_test, tot_immunogenic_test, tot_negative_test,
                             tot_score_test))
