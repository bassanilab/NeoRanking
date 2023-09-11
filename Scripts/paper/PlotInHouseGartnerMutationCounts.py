import argparse
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot in-house and Gartner mutation counts')

parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-pt', '--plot_type', type=str, default="SNV_count",
                    choices=['SNV_count', 'InDel_count', 'FS_count', 'SNV_imm_count'],
                    help='Specifies mutation type')
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

if __name__ == "__main__":

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))


    def jitter(values):
        return values + np.random.exponential(0.05, values.shape)

    comp_mutation_data = \
        pd.read_csv(os.path.join(GlobalParameters.data_dir, 'Compare_in-house_Gartner_Counts.txt'), sep='\t', header=0)

    if args.plot_type == 'SNV_count':

        comp_mutation_data['Overlap (%)'] = \
            np.round(comp_mutation_data['NCI_train_mut_seq_all in-house/Gartner overlap ratio'] * 10) * 10
        comp_mutation_data['SNV Overlap (%)'] = \
            np.round(comp_mutation_data['NCI_train_mut_seq in-house/Gartner SNV overlap ratio'] * 10) * 10

        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.scatterplot(data=comp_mutation_data, x="NCI_train_mut_seq Gartner SNV", y="NCI_train_mut_seq in-house SNV",
                            hue='SNV Overlap (%)', alpha=0.9, size='SNV Overlap (%)', edgecolors='0.2', linewidth=0,
                            sizes=(args.point_size/10, args.point_size), legend='full', palette=args.color_map)
        plt.xlabel("Gartner SM SNV count", size=args.label_size)
        plt.ylabel("SM SNV count", size=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        g.set_xscale('log')
        g.set_yscale('log')
        plt.legend(loc='upper left', title='%Gartner SM SNV found', title_fontsize=args.legend_size,
                   prop={'size': args.legend_size})
        plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

    if args.plot_type == 'InDel_count':
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.scatterplot(x=jitter(comp_mutation_data["NCI_train_mut_seq Gartner InDel"]),
                            y=jitter(comp_mutation_data["NCI_train_mut_seq in-house InDel"]),
                            alpha=0.7, s=args.point_size, edgecolors='0.1', linewidth=2)
        plt.xlabel("Gartner SM InDel count", size=args.label_size)
        plt.ylabel("SM InDel count", size=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        g.set_xscale('symlog')
        g.set_yscale('symlog')
        mv = max(max(comp_mutation_data["NCI_train_mut_seq Gartner InDel"]),
                 max(comp_mutation_data["NCI_train_mut_seq in-house InDel"]))
        mv += mv/5
        g.set_xlim((-0.5, mv))
        g.set_ylim((-0.5, mv))
        plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

    if args.plot_type == 'FS_count':
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.scatterplot(x=jitter(comp_mutation_data["NCI_train_mut_seq Gartner FS"]),
                            y=jitter(comp_mutation_data["NCI_train_mut_seq in-house FS"]),
                            alpha=0.7, s=args.point_size, edgecolors='0.1', linewidth=2)
        plt.xlabel("Gartner SM FS count", size=args.label_size)
        plt.ylabel("SM FS count", size=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        g.set_xscale('symlog')
        g.set_yscale('symlog')
        mv = max(max(comp_mutation_data["NCI_train_mut_seq Gartner FS"]),
                 max(comp_mutation_data["NCI_train_mut_seq in-house FS"]))
        mv += mv/5
        g.set_xlim((-0.5, mv))
        g.set_ylim((-0.5, mv))
        plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

    if args.plot_type == 'SNV_imm_count':
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.scatterplot(x=jitter(comp_mutation_data["NCI_train_mut_seq_imm Gartner SNV"]),
                            y=jitter(comp_mutation_data["NCI_train_mut_seq_imm in-house SNV"]),
                            alpha=0.7, s=args.point_size, edgecolor='0.2', linewidth=2)
        plt.xlabel("Gartner immunogenic\n SM SNV count", size=args.label_size)
        plt.ylabel("Immunogenic\n SM SNV count", size=args.label_size)
        plt.xticks(fontsize=args.tick_size)
        plt.yticks(fontsize=args.tick_size)
        mv = max(max(comp_mutation_data["NCI_train_mut_seq_imm Gartner"]), max(comp_mutation_data["NCI_train_mut_seq_imm in-house"]))
        mv += mv/5
        g.set_xlim((-0.5, mv))
        g.set_ylim((-0.5, mv))
        plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()



