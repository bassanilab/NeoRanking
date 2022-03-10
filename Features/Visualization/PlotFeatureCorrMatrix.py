import argparse

from DataWrangling.DataLoader import *
from scipy import stats
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from Visualization.PCAClassifyPeptideBrowser import *
from collections import Counter


parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')


args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

normalizer = None
if args.normalizer == 'q':
    normalizer = QuantileTransformer()
    print('Normalizer: QuantileTransformer')
elif args.normalizer == 'z':
    normalizer = StandardScaler()
    print('Normalizer: StandardScaler')
else:
    print('Normalizer: None')


data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                         immunogenic=['CD8', 'CD4/CD8'], min_nr_immono=0)

# perform leave one out on training set
patients = np.array(args.patients)

data, X, y = data_loader.load_patients(patients, args.input_file_tag)

if not args.features or len(args.features) == 0:
    features = Parameters().get_features()
else:
    features = args.features

p_values = {}

with PdfPages(args.pdf) as pp:

    data = data[features]
    corr = data.corr(method='pearson')
    fig, ax = plt.subplots(figsize=(30, 30))
    sns.set(font_scale=0.15)
    g = sns.clustermap(corr, cmap="viridis")
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

