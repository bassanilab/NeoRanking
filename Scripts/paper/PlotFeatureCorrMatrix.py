import argparse
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt

from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')


args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

label_size = 12
figure_size = (20, 20)

normalizer = get_normalizer(args.normalizer)

if not args.features or len(args.features) == 0:
    features = Parameters().get_ml_features()
else:
    features = args.features

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0)

# perform leave one out on training set
patients = get_valid_patients(args.patients)

data, X, y = data_loader.load_patients(patients, args.input_file_tag, peptide_type=args.peptide_type)

p_values = {}

with PdfPages(args.pdf) as pp:

    data = data[features]
    corr = data.corr(method='pearson')
    fig, ax = plt.subplots(figsize=figure_size)
    g = sns.clustermap(corr, cmap="viridis", annot_kws={"size": label_size})
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

