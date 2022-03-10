import argparse
import pickle
from sklearn.preprocessing import *
from DataWrangling.DataLoader import *
from scipy.stats import rankdata
from Classifier.OptimizationParams import *
import warnings
from Classifier.PrioritizationLearner import *

warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)


parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='prefix for patient id')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-o', '--output_file_tag', type=str, default='ML_SVM',
                    help='File tag for output file (patient)_(input_file_tag).txt')
parser.add_argument('-c', '--classifier', type=str, help='Classifier type')
parser.add_argument('-clf', '--classifier_file', type=str, help='Classifier binary file')
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


clf = PrioritizationLearner.load(args.classifier, optimizationParams, args.classifier_file)

patients = DataManager().get_valid_patients()
patients = [p for p in patients if any([ap == p for ap in args.patients])]
for patient in patients:
    if Parameters().patient_exist(patient):
        print("run predictions for: "+str(patient))
        df, X, y = data_loader.load_patients(patient, args.input_file_tag)

        if df is None:
            continue

        if args.classifier == 'DNN':
            pred_values = clf.predict(X)  # returns array of arrays
            y_pred = np.array([-v[0] for v in pred_values])
        else:
            y_pred = clf.predict_proba(X)[:, 1]
        df['ML_pred'] = y_pred

        sorted_df = df.sort_values(by=['ML_pred'], ascending=False)

        out_file = os.path.join(Parameters().get_result_dir(), patient + "_"+args.input_file_tag + "_" +
                                args.output_file_tag + ".txt")
        sorted_df.to_csv(out_file, sep="\t", header=True, index=False)

        r = rankdata(-y_pred, method='average')[y == 1]
        nr_correct = sum(r <= args.max_rank)
        nr_immuno = sum(y == 1)
        score = OptimizationParams().score(args.scorer, y, y_pred)

        sort_idx = np.argsort(r)
        print("%s\t%d\t%d\t%d\t%d\t%s\t%f" % (patient, nr_correct, nr_immuno, np.min((args.max_rank, len(y))), len(y),
                                              str(r[sort_idx]), score))
    else:
        print("Patient "+str(patient)+" does not have a data file. Skip.")

