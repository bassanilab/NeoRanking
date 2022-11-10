import argparse
import seaborn as sns
from matplotlib import pyplot as plt

from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot NeoDisc and Gartner mutation counts')

parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-f', '--mutation_info_file', type=str, help='File produced by Compare_NeoDisc_Gartner_counts.py')
parser.add_argument('-las', '--label_size', type=str, default='xx-large',
                    help='Axis label size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-tis', '--tick_size', type=str, default='x-large',
                    help='Tick size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-les', '--legend_size', type=str, default='large',
                    help='Legend size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-fiw', '--figure_width', type=float, default=8.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=8.0, help='Figure height in inches')
parser.add_argument('-dpi', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-ps', '--point_size', type=float, default=300, help='Size of points in scatterplot')
parser.add_argument('-cmp', '--color_map', type=str, default='viridis', help='Colormap for overlap')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


def jitter(values):
    return values + np.random.exponential(0.05, values.shape)


mutation_data = pd.read_csv(args.mutation_info_file, sep='\t', header=0)

mutation_data['Overlap (%)'] = np.round(mutation_data['NCI_train_mut_seq_all NeoDisc-Gartner overlap ratio']*10)*10
mutation_data['SNV Overlap (%)'] = np.round(mutation_data['NCI_train_mut_seq NeoDisc-Gartner SNV overlap ratio']*10)*10

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.scatterplot(data=mutation_data, x="NCI_train_mut_seq_all Gartner", y="NCI_train_mut_seq_all NeoDisc",
                    hue='Overlap (%)', alpha=0.9, size='Overlap (%)', sizes=(args.point_size/10, args.point_size),
                    legend='full', palette=args.color_map, edgecolors='0.2', linewidth=0)
plt.xlabel("Gartner mutation\n count", size=args.label_size)
plt.ylabel("NeoDisc mutation\n count", size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
g.set_xscale('log')
g.set_yscale('log')
plt.legend(loc='upper left', title='%Gartner in NeoDisc', title_fontsize=args.legend_size,
           prop={'size': args.legend_size})
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_all.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.scatterplot(data=mutation_data, x="NCI_train_mut_seq Gartner SNV", y="NCI_train_mut_seq NeoDisc SNV",
                    hue='SNV Overlap (%)', alpha=0.9, size='SNV Overlap (%)', edgecolors='0.2', linewidth=0,
                    sizes=(args.point_size/10, args.point_size), legend='full', palette=args.color_map)
plt.xlabel("Gartner SNV count", size=args.label_size)
plt.ylabel("NeoDisc SNV count", size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
g.set_xscale('log')
g.set_yscale('log')
plt.legend(loc='upper left', title='%Gartner in NeoDisc', title_fontsize=args.legend_size,
           prop={'size': args.legend_size})
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_SNV.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.scatterplot(x=jitter(mutation_data["NCI_train_mut_seq Gartner InDel"]),
                    y=jitter(mutation_data["NCI_train_mut_seq NeoDisc InDel"]),
                    alpha=0.7, s=args.point_size, edgecolors='0.1', linewidth=2)
plt.xlabel("Gartner InDel count", size=args.label_size)
plt.ylabel("NeoDisc InDel count", size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
g.set_xscale('symlog')
g.set_yscale('symlog')
mv = max(max(mutation_data["NCI_train_mut_seq Gartner InDel"]),
         max(mutation_data["NCI_train_mut_seq NeoDisc InDel"]))
mv += mv/5
g.set_xlim((-0.5, mv))
g.set_ylim((-0.5, mv))
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_InDel.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.scatterplot(x=jitter(mutation_data["NCI_train_mut_seq Gartner FS"]),
                    y=jitter(mutation_data["NCI_train_mut_seq NeoDisc FS"]),
                    alpha=0.7, s=args.point_size, edgecolors='0.1', linewidth=2)
plt.xlabel("Gartner FS count", size=args.label_size)
plt.ylabel("NeoDisc FS count", size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
g.set_xscale('symlog')
g.set_yscale('symlog')
mv = max(max(mutation_data["NCI_train_mut_seq Gartner FS"]),
         max(mutation_data["NCI_train_mut_seq NeoDisc FS"]))
mv += mv/5
g.set_xlim((-0.5, mv))
g.set_ylim((-0.5, mv))
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_FS.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.scatterplot(x=jitter(mutation_data["NCI_train_mut_seq_imm Gartner SNV"]),
                    y=jitter(mutation_data["NCI_train_mut_seq_imm NeoDisc SNV"]),
                    alpha=0.7, s=args.point_size, edgecolor='0.2', linewidth=2)
plt.xlabel("Gartner immunogenic\n SNV count", size=args.label_size)
plt.ylabel("NeoDisc immunogenic\n SNV count", size=args.label_size)
plt.xticks(fontsize=args.tick_size)
plt.yticks(fontsize=args.tick_size)
mv = max(max(mutation_data["NCI_train_mut_seq_imm Gartner"]), max(mutation_data["NCI_train_mut_seq_imm NeoDisc"]))
mv += mv/5
g.set_xlim((-0.5, mv))
g.set_ylim((-0.5, mv))
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_imm.png".format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()



