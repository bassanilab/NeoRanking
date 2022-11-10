import argparse
import time

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-c', '--classifier', type=str, default='SVM', help='classifier to use')
parser.add_argument('-d', '--classifier_dir', type=str, default=Parameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-tf', '--train_res_file_re', type=str, nargs='+', help='classifier to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in Hyperopt')
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

def parse_train_file(train_file):
    params_ = None
    patients_ = None
    clf_file_ = None
    lines = [line for line in train_file]
    for line in lines:
        if line.startswith("Hyperopt"):
            params_ = parse_hyperopt_results(line)
        elif line.startswith("Training patients"):
            patients_ = line.split(":")[1].strip().split(",")
        elif line.startswith("Saved to"):
            clf_file_ = line.split(" ")[2].strip()

    return params_, patients_, clf_file_


def parse_hyperopt_results(line):
    fields = line.replace("Hyperopt: ", "").split("; ")
    for f in fields:
        k, v = f.split('=')
        if k == 'Params':
            return ast.literal_eval(v)


normalizer = get_normalizer(args.normalizer)
encodings = read_cat_encodings(args.patients_train[0], args.peptide_type)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                         immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                         max_netmhc_rank=args.max_rank_netmhc, cat_encoders=encodings)

patients_train = \
    get_valid_patients(dataset=args.patients_train, peptide_type=args.peptide_type) \
        if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

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

learner = PrioritizationLearner(args.classifier, 'sum_exp_rank', optimizationParams)

train_res_files = []
for wc in args.train_res_file_re:
    train_res_files = train_res_files + glob.glob(os.path.join(args.classifier_dir, wc))

for train_res_file in train_res_files:
    with open(train_res_file, mode='r') as result_file:
        params, patients, classifier_file = parse_train_file(result_file)

    if patients is None:
        # add if patient info was missing
        with open(train_res_file, mode='a') as result_file:
            result_file.write('Training patients: {0}\n'.format(','.join(patients_train)))

    # put classifier binary file in same directory as _train.txt file
    ext = DataManager().get_classifier_ext(args.classifier)
    classifier_file = result_file.name.replace("_train.txt", "."+ext)

    if classifier_file is None:
        # add only if classifier file line was missing
        with open(train_res_file, mode='a') as result_file:
            result_file.write('Saved to {0:s}\n'.format(classifier_file))

    # fit best classifier on all data
    classifier = learner.fit_classifier(X_train, y_train, classifier=None, params=params)
    PrioritizationLearner.save_classifier(args.classifier, classifier, classifier_file)




