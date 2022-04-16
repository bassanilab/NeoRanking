import argparse
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-cfl', '--classifier_long', type=str, default='', help='classifier file for long peptides')
parser.add_argument('-cfs', '--classifier_short', type=str, default='', help='classifier file for short peptides')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-il', '--input_file_tag_long', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file for long peptides (patient)_long_(input_file_tag).txt')
parser.add_argument('-is', '--input_file_tag_short', type=str, default='stab_chop',
                    help='File tag for neodisc input file for short peptides (patient)_short_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_short_peptide_test', help='Short info for classifier run')
parser.add_argument('-fl', '--features_long', type=str, nargs='+', help='Features used by classifier for long peptides')
parser.add_argument('-fs', '--features_short', type=str, nargs='+', help='Features used by classifier for short peptides')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20,
                    help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-manl', '--max_nr_long', type=int, default=100,
                    help='Maximal number of ranked long peptides to consider for short peptide ranking')
parser.add_argument('-minl', '--min_nr_long', type=int, default=30,
                    help='If nr mutations smaller than minl, no filtering is performed for mutations')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-cat', '--cat_to_num', dest='cat_to_num', action='store_true',
                    help='convert categories to numbers')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

data_loader_long = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_long,
                              mutation_types=args.mutation_types, response_types=args.response_types,
                              immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num,
                              max_netmhc_rank=10000)

cat_features = [f for f in args.features_long if f in Parameters().get_categorical_features()]
cat_idx = [i for i, f in enumerate(args.features_long) if f in Parameters().get_categorical_features()]

optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                        cat_dims=data_loader_long.get_categorical_dim(), input_shape=[len(args.features_long)])

with open(args.classifier_long, mode='r') as result_file:
    classifier_long_tag = os.path.basename(args.classifier_long).split('_')[0]
    classifier_long = \
        PrioritizationLearner.load_classifier(classifier_long_tag, optimizationParams, args.classifier_long)
    learner = PrioritizationLearner(classifier_long_tag, args.scorer, optimizationParams, verbose=args.verbose)

with open(args.classifier_short, mode='r') as result_file:
    classifier_short_tag = os.path.basename(args.classifier_short).split('_')[0]
    classifier_short = \
        PrioritizationLearner.load_classifier(classifier_short_tag, optimizationParams, args.classifier_short)

patients_test_long = \
    get_valid_patients(patients=args.patients_test, peptide_type='long') \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type='long')
patients_test_short = \
    get_valid_patients(patients=args.patients_test, peptide_type='short') \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type='short')

mgr = DataManager()
patients_test_long = patients_test_long.intersection(mgr.get_immunogenic_patients('long'))
patients_test_short = patients_test_short.intersection(mgr.get_immunogenic_patients('short'))
patients_test = sorted(patients_test_long.intersection(patients_test_short))

tot_negative_test = 0
tot_correct_test = 0
tot_immunogenic_test = 0
tot_score_test = 0
data_loader_short = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features_short,
                               mutation_types=args.mutation_types, response_types=args.response_types,
                               immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num,
                               max_netmhc_rank=10000)


clf_tag = "{0}_{1}".format(classifier_long_tag, classifier_short_tag)
with open(DataManager().get_result_file(clf_tag, args.run_id), mode='w') as result_file:
    for arg in vars(args):
        result_file.write(f"{arg}={getattr(args, arg)}\n")

    if patients_test is not None:
        for p in patients_test:
            data_test, X_test, y_test = \
                data_loader_long.load_patients(p, args.input_file_tag_long, 'long', verbose=False)
            min_nr = args.min_nr_long
            if data_test.shape[0]*0.05 < min_nr:
                nr_long = min(data_test.shape[0], min_nr)
            else:
                nr_long = min(args.max_nr_long, round((data_test.shape[0])*0.05))
#            nr_long = args.max_nr_long
            mutant_ids = \
                learner.get_top_n_mutation_ids(classifier_long, data_test, X_test.to_numpy(), max_rank=nr_long)

            data_test, X_test, y_test = \
                data_loader_short.load_patients(p, args.input_file_tag_short, 'short', verbose=False)
            idx = data_test.apply(lambda row: row['mut_seqid'] in mutant_ids, axis=1)
            y_filtered = np.array([y_test[i] for i in range(len(y_test)) if idx[i]])
            y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                learner.test_classifier(classifier_short, p, X_test.loc[idx, :], y_filtered, max_rank=args.max_rank,
                                        report_file=result_file)

            tot_negative_test += len(y_test) - nr_immuno
            tot_correct_test += nr_correct
            tot_immunogenic_test += nr_immuno
            tot_score_test += score

    if args.verbose > 0:
        print('nr_patients\trun_id\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train'.format(args.max_rank))
        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}'.
              format(len(patients_test), args.run_id, tot_correct_test, tot_immunogenic_test, tot_negative_test,
                     tot_score_test))
    result_file.write('nr_patients\tnr_correct_top{0}\tnr_immunogenic\tnr_negative\tscore_train\n'.
                      format(args.max_rank))
    result_file.write('{0}\t{1}\t{2}\t{3}\t{4:.3f}\n'.
                      format(len(patients_test), tot_correct_test, tot_immunogenic_test, tot_negative_test, tot_score_test))
