import argparse

from DataWrangling.Transform_Data import *
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-p', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-nn', '--nr_negative', type=int, default=-1, help='Maximal number of non immunogenic samples')
parser.add_argument('-eg', '--excluded_genes', type=str, nargs='+', help='genes excluded from prioritization')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')

args = parser.parse_args()


def calc_save_encoding(x, y, file):

    encoding_file.write(Encoder.get_file_header()+"\n")

    cat_features = [f for f in x.columns if f in Parameters().get_categorical_features()]

    for f in cat_features:
        l_enc = Encoder(f)
        l_enc.fit(x[f].values, y)
        l_enc.append_to_file(file)


with open(Parameters().get_cat_to_num_info_file(args.patients_train[0], args.peptide_type), mode='w') as encoding_file:

    for arg in vars(args):
        encoding_file.write(f"#{arg}={getattr(args, arg)}\n")
        print(f"{arg}={getattr(args, arg)}")

    data_loader = DataLoader(transformer=DataTransformer(), features=args.features, mutation_types=args.mutation_types,
                             response_types=['CD8', 'CD4/CD8', 'negative'], immunogenic=args.immunogenic, min_nr_immuno=0,
                             cat_to_num=False, max_netmhc_rank=args.max_rank_netmhc)

    patients_train = \
        get_valid_patients(patients=args.patients_train, peptide_type=args.peptide_type) \
            if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

    df_train, X_train, y_train = data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type)
    calc_save_encoding(X_train, y_train, encoding_file)

