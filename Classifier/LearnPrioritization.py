import argparse

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-c', '--classifier', type=str, default='SVM', help='classifier to use')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-loo', '--leave_one_out', dest='leave_one_out', action='store_true',
                    help='leave-one-out CV for training set')
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
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in RandomSearchCV')
parser.add_argument('-nc', '--nr_classifiers', type=int, default=1,
                    help='Number of best classifiers included for voting')
parser.add_argument('-cv', '--nr_cv', type=int, default=5, help='Number of CV layers in RandomSearchCV')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-e', '--nr_epoch', type=int, default=150, help='Number of epochs for DNN training')
parser.add_argument('-ep', '--early_stopping_patience', type=int, default=150,
                    help='Patience for early stopping for DNN training')
parser.add_argument('-b', '--batch_size', type=int, default=150, help='Batch size for DNN training')
parser.add_argument('-cat', '--cat_to_num', dest='cat_to_num', action='store_true',
                    help='convert categories to numbers')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-ct', '--combine_test', dest='combine_test', action='store_true', help='Combine test data')
parser.add_argument('-cf', '--classifier_file', type=str, default='', help='classifier file instead of training')

args = parser.parse_args()

with open(DataManager().get_result_file(args.classifier, args.run_id, args.peptide_type), mode='w') \
        as result_file:
    for arg in vars(args):
        result_file.write(f"{arg}={getattr(args, arg)}\n")
        print(f"{arg}={getattr(args, arg)}")

    normalizer = get_normalizer(args.normalizer)

    data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                             mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                             immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num,
                             max_netmhc_rank=10000)

    cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
    cat_idx = [i for i, f in enumerate(args.features) if f in Parameters().get_categorical_features()]

    if args.classifier_file == '' or not os.path.exists(args.classifier_file):
        patients_train = \
            get_valid_patients(patients=args.patients_train, peptide_type=args.peptide_type) \
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

        if args.leave_one_out and len(patients_train) < 2:
            print("Not enough patients for leave-one-out testing. Need 2 patients at least.")

        if args.leave_one_out:  # perform leave one out on training set
            for p in patients_train:
                data_test, X_test, y_test = data_loader.load_patients(p, args.input_file_tag, args.peptide_type)
                if y_test is None or sum(y_test == 1) == 0:
                    continue

                data_train, X_train, y_train = \
                    data_loader.load_patients(patients_train[patients_train != p], args.input_file_tag, args.peptide_type)

                if args.peptide_type == 'short':
                    class_ratio = sum(y_train == 1)/sum(y_train == 0)
                else:
                    class_ratio = None

                optimizationParams = \
                    OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                       cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)],
                                       class_ratio=class_ratio)

                learner = PrioritizationLearner(args.classifier, args.scorer, optimizationParams, verbose=args.verbose,
                                                nr_iter=args.nr_iter, nr_classifiers=args.nr_classifiers, nr_cv=args.nr_cv,
                                                shuffle=args.shuffle, nr_epochs=args.nr_epoch,
                                                patience=args.early_stopping_patience, batch_size=args.batch_size)

                cvres, best_classifier, best_score, best_params = learner.optimize_classifier(X_train.to_numpy(), y_train)
                y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                    learner.test_classifier(best_classifier, p, X_test, y_test, max_rank=args.max_rank)

                tot_negative_train += len(y_train) - nr_immuno
                tot_correct_train += nr_correct
                tot_immunogenic_train += nr_immuno
                tot_score_train += score

                if best_score > best_score_train:
                    best_score_train = best_score
                    best_param_train = best_params
                    best_classifier_train = best_classifier

            data_train, X_train, y_train = data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type)
            best_classifier_train = learner.fit_classifier(X_train, y_train, classifier=best_classifier_train)
        else:
            data_train, X_train, y_train = data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type)

            if args.peptide_type == 'short':
                class_ratio = sum(y_train == 1)/sum(y_train == 0)
            else:
                class_ratio = None

            optimizationParams = \
                OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                   cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)],
                                   class_ratio=class_ratio)

            learner = PrioritizationLearner(args.classifier, args.scorer, optimizationParams, verbose=args.verbose,
                                            nr_iter=args.nr_iter, nr_classifiers=args.nr_classifiers, nr_cv=args.nr_cv,
                                            shuffle=args.shuffle, nr_epochs=args.nr_epoch,
                                            patience=args.early_stopping_patience, batch_size=args.batch_size)

            cvres, best_classifier, best_score, best_params = learner.optimize_classifier(X_train.to_numpy(), y_train)
            best_score_train = best_score
            best_param_train = best_params
            best_classifier_train = best_classifier

        classifier_file = DataManager().get_classifier_file(args.classifier, args.run_id, args.peptide_type)

        # fit best classifier on all data
        PrioritizationLearner.save_classifier(args.classifier, best_classifier_train, classifier_file)
        if args.verbose > 1:
            print('Classifier = {0:s}, Scorer = {1:s}'.format(args.classifier, args.scorer))
            print('Best training params: ' + str(best_param_train) + ', ' + args.scorer + ': ' + str(best_score_train))
            print('Saved to {0:s}'.format(classifier_file))

            result_file.write('Training patients: {0}\n'.format(','.join(patients_train)))
            result_file.write('Classifier = {0:s}, Scorer = {1:s}\n'.format(args.classifier, args.scorer))
            result_file.write('Best training params: {0}, {1}\n'.format(str(best_param_train), str(best_score_train)))
            result_file.write('Saved to {0:s}\n'.format(classifier_file))

    else:
        optimizationParams = \
            OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                               cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])

        learner = PrioritizationLearner(args.classifier, args.scorer, optimizationParams, verbose=args.verbose,
                                        nr_iter=args.nr_iter, nr_classifiers=args.nr_classifiers, nr_cv=args.nr_cv,
                                        shuffle=args.shuffle, nr_epochs=args.nr_epoch,
                                        patience=args.early_stopping_patience, batch_size=args.batch_size)

        best_classifier_train = \
            PrioritizationLearner.load_classifier(args.classifier, optimizationParams, args.classifier_file)
        result_file.write('Classifier imported from = {0}\n'.format(args.classifier_file))

    patients_test = \
        get_valid_patients(patients=args.patients_test, peptide_type=args.peptide_type) \
            if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)

    mgr = DataManager()
    patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type)))

    tot_negative_test = 0
    tot_correct_test = 0
    tot_immunogenic_test = 0
    tot_score_test = 0
    response_types = ['not_tested', 'negative', 'CD8', 'CD4/CD8']
    data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                             mutation_types=args.mutation_types, response_types=response_types,
                             immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num,
                             max_netmhc_rank=10000)

    if patients_test is not None:
        result_file.write('Test patients: {0}\n'.format(','.join(patients_test)))
        if not args.combine_test:
            patients_test = np.sort(patients_test)
            for p in patients_test:
                data_test, X_test, y_test = \
                    data_loader.load_patients(p, args.input_file_tag, args.peptide_type, verbose=False)
                y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                    learner.test_classifier(best_classifier_train, p, X_test.to_numpy(), y_test,
                                            max_rank=args.max_rank, report_file=result_file)

                tot_negative_test += len(y_test) - nr_immuno
                tot_correct_test += nr_correct
                tot_immunogenic_test += nr_immuno
                tot_score_test += score
        else:
            data_test, X_test, y_test = data_loader.load_patients(patients_test, args.input_file_tag, args.peptide_type)
            y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                learner.test_classifier(best_classifier_train, ','.join(patients_test), X_test.to_numpy(), y_test,
                                        max_rank=args.max_rank, report_file=result_file)

            tot_negative_test += len(y_test) - nr_immuno
            tot_correct_test += nr_correct
            tot_immunogenic_test += nr_immuno
            tot_score_test += score

    if args.verbose > 0:
        print('nr_patients\trun_id\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train'.format(args.max_rank))
        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}'.format(len(patients_test), classifier_file, tot_correct_test,
                                                        tot_immunogenic_test, tot_negative_test, tot_score_test))

    result_file.write('nr_patients\trun_id\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train\n'.
                      format(args.max_rank))
    result_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}\n'.
                      format(len(patients_test), classifier_file, tot_correct_test, tot_immunogenic_test,
                             tot_negative_test, tot_score_test))
