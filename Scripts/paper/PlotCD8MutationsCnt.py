import argparse
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

from DataWrangling.DataLoader import *
from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

data_loader = DataLoader(response_types=['not_tested', 'CD8', 'CD4/CD8', 'negative'], mutation_types=['SNV'],
                         immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0)

patients = get_valid_patients(None)

mutations_cnt = []
CD8_mutation_cnt = []
patients_group = []
processed_patients = []

for p in patients:
    data, X, y = data_loader.load_patients(p, 'rt', 'long')
    if data is not None:
        mutations_cnt.append(len(y))
        CD8_mutation_cnt.append(sum(y == 1))
        patients_group.append(get_patient_group(p))
        processed_patients.append(p)

mutation_data = pd.DataFrame({'Mutation count': mutations_cnt,
                              'CD8+ immunogenic mutation counts': CD8_mutation_cnt,
                              'Patient group': patients_group, 'Processed patients': processed_patients})

with PdfPages(args.pdf) as pp:

    fig = plt.figure(figsize=(20, 15))
    fig.clf()
    g = sns.lmplot(data=mutation_data, x="Mutation count", y="CD8+ immunogenic mutation counts", hue='Patient group',
                   scatter_kws={"alpha": 0.7, "s": 100}, logx=True, legend=False)
    plt.xlabel("Mutation count", size=15)
    plt.ylabel("CD8+ immunogenic mutation count", size=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    ax = g.axes[0][0]
    ax.set_xscale('log')
    plt.legend(loc='upper left')
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

