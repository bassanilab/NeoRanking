import argparse
import seaborn as sns
from matplotlib import pyplot as plt

from DataWrangling.DataTransformer import *
from Utils.Util_fct import *


parser = argparse.ArgumentParser(description='Plot correlation between mutation and immunogenic mutation counts')

parser.add_argument('-fp', '--file_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="svg", help='File type for plot (png, svg or pdf')
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

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

data_loader = DataTransformer(response_types=['not_tested', 'CD8', 'CD4/CD8', 'negative'], mutation_types=['SNV'],
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
plot_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}.{1}".format(args.file_prefix, args.file_type))
plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


