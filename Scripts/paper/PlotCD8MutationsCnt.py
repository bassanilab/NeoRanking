import argparse
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

from Utils.Util_fct import *
from Utils.DataManager import DataManager


parser = argparse.ArgumentParser(description='Plot correlation between mutation counts and immunogenic mutation counts')

parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-las', '--label_size', type=str, default='30',
                    help='Axis label size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-tis', '--tick_size', type=str, default='15',
                    help='Tick size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-les', '--legend_size', type=str, default='15',
                    help='Legend size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-fiw', '--figure_width', type=float, default=30.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=20.0, help='Figure height in inches')
parser.add_argument('-dpi', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-ps', '--point_size', type=float, default=100, help='Size of points in scatterplot')
parser.add_argument('-cmp', '--color_map', type=str, default='viridis', help='Colormap for overlap')

if __name__ == "__main__":
    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    mutations_cnt = []
    CD8_mutation_cnt = []
    patients_group = []
    processed_patients = []

    mutation_data = DataManager.filter_original_data(peptide_type='mutation', ml_row_selection=True)

    for patient in mutation_data['patient'].unique():
        data, X, y = DataManager.get_processed_data(peptide_type=args.peptide_type, objective='sel',
                                                    patient=patient, sample=False)
        if data is not None:
            mutations_cnt.append(len(y))
            CD8_mutation_cnt.append(sum(y == 1))
            patients_group.append(data.dataset.unique()[0])
            processed_patients.append(patient)

    mutation_data = pd.DataFrame({'Mut-seq count': mutations_cnt,
                                  'Mut-seq_imm count': CD8_mutation_cnt,
                                  'Patient group': patients_group, 'Processed patients': processed_patients})

    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    g = sns.lmplot(data=mutation_data, x="Mut-seq count", y="Mut-seq_imm count", hue='Patient group',
                   scatter_kws={"alpha": 0.8, "s": args.point_size}, logx=True, legend=False,
                   palette=dict(NCI="blue", TESLA="Orange", HiTIDE="Green"), hue_order=['NCI', 'TESLA', 'HiTIDE'])
    plt.xlabel("Mut-seq count", size=args.label_size)
    plt.ylabel("Mut-seq_imm count", size=args.label_size)
    plt.xticks(fontsize=args.tick_size)
    plt.yticks(fontsize=args.tick_size)
    ax = g.axes[0][0]
    ax.set_xscale('log')
    plt.legend(loc='upper left', prop={'size': args.legend_size})
    plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()


