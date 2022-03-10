import argparse
import pickle
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from sklearn.preprocessing import *
import sys


parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_SVM',
                    help='File tag for classifier pickle output file')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')


args = parser.parse_args()

print('Input file tag: '+str(args.input_file_tag))
print('Run ID: '+str(args.run_id))
print('Features: '+str(args.features))
print('Data set: '+str(args.patients))

if args.features is None:
    features = Parameters().get_features()
else:
    features = args.features


data_loader = DataLoader(transformer=None, normalizer=None, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types)

# perform leave one out on training set
patients = np.array(args.patients)
data, X, y = data_loader.load_patients(patients, args.input_file_tag)


out_file = os.path.join(Parameters().get_result_dir(), "all_data_selection.txt")
data.to_csv(out_file, sep="\t", header=True, index=False)

nr_nas = data.isna().sum(axis=0)
print("Feature\tnr_na")
for i in np.arange(data.shape[1]):
    print("{0}\t{1}".format(data.columns[i], nr_nas[i]))
