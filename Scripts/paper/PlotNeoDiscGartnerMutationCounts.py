import argparse
import seaborn as sns
from matplotlib import pyplot as plt

from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot NeoDisc and Gartner mutation counts')

parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-f', '--mutation_info_file', type=str, help='File produced by Compare_NeoDisc_Gartner_counts.py')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


mutation_data = pd.read_csv(args.mutation_info_file, sep='\t', header=0)

mutation_data['Overlap (%)'] = np.round(mutation_data['NeoDisc-Gartner mut_all overlap ratio']*10)*10
mutation_data['SNV Overlap (%)'] = np.round(mutation_data['NeoDisc-Gartner SNV mut overlap ratio']*10)*10

fig = plt.figure(figsize=(10, 10))
fig.clf()
g = sns.scatterplot(data=mutation_data, x="Gartner mut_all", y="NeoDisc mut_all", hue='Overlap (%)',
                    alpha=0.85, size='Overlap (%)', sizes=(40, 300), legend='full', palette='viridis')
plt.xlabel("NCI mut count", size=50)
plt.ylabel("NeoDisc mut count", size=50)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
g.set_xscale('log')
g.set_yscale('log')
plt.legend(loc='upper left', title='%NCI mut in NeoDisc', title_fontsize=25, prop={'size': 20})
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(10, 10))
fig.clf()
g = sns.scatterplot(data=mutation_data, x="Gartner SNV mut", y="NeoDisc SNV mut", hue='SNV Overlap (%)',
                    alpha=0.85, size='SNV Overlap (%)', sizes=(40, 300), legend='full', palette='viridis')
plt.xlabel("NCI SNV mut count", size=30)
plt.ylabel("NeoDisc SNV mut count", size=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
g.set_xscale('log')
g.set_yscale('log')
plt.legend(loc='upper left', title='%NCI SNV mut\n in NeoDisc', title_fontsize=20, prop={'size': 20})
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_SNV.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=(10, 10))
fig.clf()
g = sns.scatterplot(data=mutation_data, x="Gartner InDel mut", y="NeoDisc InDel mut",
                    alpha=0.7, s=300)
plt.xlabel("NCI InDel mut count", size=30)
plt.ylabel("NeoDisc InDel mut count", size=30)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
g.set_xscale('log')
g.set_yscale('log')
# g.set(xlim=(0, 55), ylim=(0, 55))
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_InDel.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(10, 10))
fig.clf()
g = sns.scatterplot(data=mutation_data, x="Gartner FS mut", y="NeoDisc FS mut", alpha=0.7, s=300)
plt.xlabel("NCI FS mut count", size=30)
plt.ylabel("NeoDisc FS mut count", size=30)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
g.set_xscale('log')
g.set_yscale('log')
# g.set(xlim=(0, 55), ylim=(0, 55))
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_FS.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight')
plt.close()


def jitter(values):
    return values + np.random.normal(0.0, 0.05, values.shape)


fig = plt.figure(figsize=(10, 10))
fig.clf()
g = sns.scatterplot(x=jitter(mutation_data["Gartner mut_imm"]), y=jitter(mutation_data["NeoDisc mut_imm"]),
                    alpha=0.6, s=300, edgecolor="gray")
plt.xlabel("NCI mut_imm count", size=30)
plt.ylabel("NeoDisc mut_imm count", size=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
g.set(xlim=(-1, 12), ylim=(-1, 12))
g.figure.tight_layout()
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_CD8.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight')
plt.close()



