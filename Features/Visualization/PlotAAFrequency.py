import argparse

import numpy as np

from DataWrangling.DataLoader import *
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
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')


args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

features = ['aa_mutant', 'mut_is_binding_pos_0', 'mut_allele_0']
data_loader = DataLoader(transformer=DataTransformer(), normalizer=None, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0)

# perform leave one out on training set
patients = np.array(args.patients)

data, X, y = data_loader.load_patients(patients, args.input_file_tag)
data = data.astype({'response': 'int32', 'mut_is_binding_pos_0': 'int32'})

counter_i_b = Counter(data.loc[(data['response'] == 1) & (data['mut_is_binding_pos_0'] == 1), 'aa_mutant'])
counter_i_n = Counter(data.loc[(data['response'] == 1) & (data['mut_is_binding_pos_0'] == 0), 'aa_mutant'])
counter_n_b = Counter(data.loc[(data['response'] == -1) & (data['mut_is_binding_pos_0'] == 1), 'aa_mutant'])
counter_n_n = Counter(data.loc[(data['response'] == -1) & (data['mut_is_binding_pos_0'] == 0), 'aa_mutant'])

H_alleles = data.loc[(data['response'] == 1) & (data['mut_is_binding_pos_0'] == 0) & (data['aa_mutant'] == 'H'),
                     'mut_allele_0']
print(H_alleles)

aas = Parameters().get_AA()

aa_freq = []
response = []
aa_lst = []

for aa in aas:
    aa_lst = np.append(aa_lst, [aa, aa])
    aa_freq = np.append(aa_freq, counter_i_n[aa]/sum(counter_i_n.values()))
    aa_freq = np.append(aa_freq, counter_n_n[aa]/sum(counter_n_n.values()))
    response = np.append(response, [1, 0])

with PdfPages(args.pdf) as pp:

    df = pd.DataFrame({'AminoAcid': aa_lst, 'aa_freq': aa_freq, 'response': response})

    g = sns.barplot(x='AminoAcid', y='aa_freq', hue="response", data=df, palette='pastel')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()
