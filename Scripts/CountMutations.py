import argparse
from DataWrangling.DataLoader import *
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')


args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

data_loader = DataLoader(transformer=None, normalizer=None, features=None,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0)

data, X, y = data_loader.load_patients(args.patients, args.input_file_tag)

sub_ttl = "Mutation types = {}".format(args.mutation_types)

if args.pdf:

    with PdfPages(args.pdf) as pp:

        sns.set(font_scale=1.8)
        order = ['not_tested', 'negative', 'CD8', 'CD4']
        g = sns.catplot(x='response_type', col='patient', kind="count", data=data, order=order, col_wrap=5,
                        palette='pastel')
        g.fig.suptitle(sub_ttl)
        g.fig.set_size_inches(37, 20)
        g.fig.subplots_adjust(top=0.9)
        g.set(yscale='log', ylim=(0.9, 5000))

        # iterate through axes
        for ax in g.axes.ravel():
            # add annotations
            for c in ax.containers:
                labels = [f'{(v.get_height()):.0f}' for v in c]
                ax.bar_label(c, labels=labels, label_type='edge')
            ax.margins(y=0.2)

        g.figure.tight_layout()
        pp.savefig(g.figure)
