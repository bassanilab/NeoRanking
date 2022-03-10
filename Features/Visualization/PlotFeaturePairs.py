import argparse
from DataWrangling.DataLoader import *
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-fp', '--feature_pairs', type=str, nargs='+', help='Features pair for pair plots')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-reg', '--regression_line', dest='regression_line', action='store_true',
                    help='draw regression line')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')


args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

normalizer = None
if args.normalizer == 'q':
    normalizer = QuantileTransformer()
    if args.verbose > 0:
        print('Normalizer: QuantileTransformer')
elif args.normalizer == 'z':
    normalizer = StandardScaler()
    if args.verbose > 0:
        print('Normalizer: StandardScaler')
else:
    if args.verbose > 0:
        print('Normalizer: None')

features = []
for fp in args.feature_pairs:
    (f1, f2) = fp.split(',')
    features.append(f1)
    features.append(f2)

features = np.unique(features)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0)

# perform leave one out on training set
dataManager = DataManager()
patients_with_data = dataManager.get_valid_patients(args.peptide_type)
if not args.patients or len(args.patients) == 0:
    patients = patients_with_data
elif args.patients[0] == 'Rosenberg':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('all')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Gartner_train':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('gartner_train')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Gartner_test':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('gartner_test')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Gartner':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('gartner')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'Parkhurst':
    annotator = RosenbergImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients('parkhurst')
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'HiTIDE':
    annotator = NeoDiscImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients()
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'TESLA':
    annotator = TESLAImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients()
    patients = [p for p in patients if p in patients_with_data]
elif args.patients[0] == 'TESLA/HiTIDE':
    annotator = TESLAImmunogenicityAnnotatorLong(dataManager)
    patients = annotator.get_patients()
    annotator = NeoDiscImmunogenicityAnnotatorLong(dataManager)
    patients = patients + annotator.get_patients()
    patients = [p for p in patients if p in patients_with_data]
else:
    patients = np.array(args.patients)

data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)

df = pd.DataFrame(data=X_train, columns=features)
df.insert(0, 'response', y_train)
all_num_cols = Parameters().get_numerical_features()
num_cols = [c for c in df.columns if c in all_num_cols]
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors='coerce')

p_values = {}

font_size = 1.2

if args.verbose > 0 and args.pdf:

    with PdfPages(args.pdf) as pp:

        for fp in args.feature_pairs:

            (f1, f2) = fp.split(',')
            if f1 in Parameters().get_numerical_features() and f2 in Parameters().get_numerical_features():
                g = sns.lmplot(x=f1, y=f2, hue="response", data=df, palette='seismic',
                               line_kws={'color': 'lightgray'}, legend=False, fit_reg=args.regression_line)
                sns.set(font_scale=font_size)
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

            if f1 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
                    f2 in Parameters().get_numerical_features():
                g = sns.catplot(x=f1, y=f2, hue="response", kind="violin", split=True, data=df,
                                palette='seismic', legend=False)
                sns.set(font_scale=font_size)
                plt.legend(loc='upper right')
                xlabels = ["{0:.1f}".format(l) for l in np.unique(df[f1])]
                g.set_xticklabels(xlabels)
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

            if f2 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
                    f1 in Parameters().get_numerical_features():
                g = sns.catplot(x=f2, y=f1, hue="response", kind="violin", split=True, data=df,
                                palette='seismic', legend=False)
                g.figure.tight_layout()
                plt.legend(loc='upper right')
                sns.set(font_scale=font_size)
                xlabels = ["{0:.1f}".format(l) for l in np.unique(df[f2])]
                g.set_xticklabels(xlabels)
                pp.savefig(g.figure)
                g.figure.clf()

            if f2 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
                    f1 in Parameters().get_categorical_features()+Parameters().get_ordinal_features():
                g = sns.catplot(x=f1, col=f2, hue="response", kind="count", data=df, palette='seismic')
                sns.set(font_scale=font_size)
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

