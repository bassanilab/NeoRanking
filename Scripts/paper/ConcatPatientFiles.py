import argparse

from DataWrangling.DataLoader import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Write all patient data into one file')
parser.add_argument('-f', '--data_file', type=str, help='tabular output file')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='Patients')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=0, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included for testing')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-i', '--input_file_tag', type=str, default='rt',
                    help='File tag for neodisc input file (patient)_(peptide_type)_(input_file_tag).txt')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


patients = \
    get_valid_patients(patients=args.patients, peptide_type=args.peptide_type) \
        if args.patients and len(args.patients) > 0 else get_valid_patients(peptide_type=args.peptide_type)

mgr = DataManager()
patients_test = sorted(patients.intersection(mgr.get_immunogenic_patients(args.peptide_type)))

data_loader = DataLoader(mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=args.min_nr_immuno, cat_encoders=None)

data, X, y = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type, verbose=True)
data.to_csv(path_or_buf=args.data_file, sep="\t", index=False, header=True)


