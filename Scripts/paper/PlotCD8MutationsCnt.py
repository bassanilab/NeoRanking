import argparse
import seaborn as sns
from matplotlib import pyplot as plt

from DataWrangling.DataLoader import *
from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot correlation between mutation and immunogenic mutation counts')

parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')

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

mutation_data = pd.DataFrame({'Mut_all count': mutations_cnt,
                              'Mut_imm count': CD8_mutation_cnt,
                              'Patient group': patients_group, 'Processed patients': processed_patients})


fig = plt.figure(figsize=(30, 20))
fig.clf()
g = sns.lmplot(data=mutation_data, x="Mut_all count", y="Mut_imm count", hue='Patient group',
               scatter_kws={"alpha": 0.7, "s": 100}, logx=True, legend=False)
plt.xlabel("Mut_all count", size=30)
plt.ylabel("Mut_imm count", size=30)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = g.axes[0][0]
ax.set_xscale('log')
plt.legend(loc='upper left', prop={'size': 15})
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight')
plt.close()


