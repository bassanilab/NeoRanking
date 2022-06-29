import argparse
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from DataWrangling.DataLoader import DataLoader
from DataWrangling.Transform_Data import DataTransformer
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-fp', '--feature_pairs', type=str, nargs='+', help='Features pair for pair plots')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='n',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-reg', '--regression_line', dest='regression_line', action='store_true',
                    help='draw regression line')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-s', '--nr_samples', type=int, default=10000, help='Nr of Sampled rows from dataframe for plots')
parser.add_argument('-lep', '--legend_position', type=str, default='best', help='Legend position in plot')
parser.add_argument('-kdp', '--kd_plot', dest='kd_plot', action='store_true', help='Plot kernel density')
parser.add_argument('-d', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')

color_CD8 = 'darkorange'
color_negative = 'royalblue'
axis_label_size = 12
legend_size = 10
tickmark_size = 8
figure_width = 15
figure_height = 15

args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

features = []
for fp in args.feature_pairs:
    (f1, f2) = fp.split(',')
    features.append(f1)
    features.append(f2)

features = np.unique(features)

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

# perform leave one out on training set
patients = get_valid_patients(args.patients)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0)

data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)

X_1 = X_train.loc[y_train == 1, :]
X_1['response'] = "1"
X_0 = X_train.loc[y_train == 0, :]
X_0['response'] = "0"
if X_0.shape[0] > args.nr_samples:
    X_0 = X_0.sample(n=args.nr_samples)
df = pd.concat([X_1, X_0])
all_num_cols = Parameters().get_numerical_features()
num_cols = [c for c in df.columns if c in all_num_cols]
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors='coerce')
df_small = df.head(1000).sort_values(by=['response'])

p_values = {}

if args.verbose > 0 and args.pdf:

    with PdfPages(args.pdf) as pp:

        for fp in args.feature_pairs:

            fig = plt.figure()
            fig.set_figheight(figure_height)
            fig.set_figwidth(figure_width)
            (f1, f2) = fp.split(',')
            if type(normalizer) == dict:
                norm_f1 = normalizer[f1]
                normalizer_name = get_normalizer_name(norm_f1)
                norm_f2 = normalizer[f2]
                normalizer_name = get_normalizer_name(norm_f2)
            else:
                norm_f1 = normalizer
                normalizer_name = get_normalizer_name(norm_f1)
                norm_f2 = normalizer
                normalizer_name = get_normalizer_name(norm_f2)

            if f1 in Parameters().get_numerical_features() and f2 in Parameters().get_numerical_features():
                g = sns.lmplot(data=df, x=f1, y=f2, hue="response", line_kws={'alpha': 0.7}, hue_order=['0', '1'],
                               scatter_kws={'alpha': 0.5, 's': 1}, legend=False, fit_reg=args.regression_line,
                               palette={'0': color_negative, '1': color_CD8})
                if args.kd_plot:
                    g = g.map_dataframe(sns.kdeplot, x=f1, y=f2, hue="response",
                                        hue_order=['0', '1'], alpha=0.4, fill=True)
                plt.legend(title='Response type', loc=args.legend_position, labels=['negative', 'CD8+'],
                           fontsize=legend_size, title_fontsize=legend_size)

                plt.xlabel(feature_dict[f1], size=axis_label_size)
                plt.ylabel(feature_dict[f2], size=axis_label_size)

                if norm_f1 is not None:
                    x_ticks = g.ax.get_xticks()
                    x_tick_label = \
                        ["{0:.1e}".format(x[0]) for x in norm_f1.inverse_transform(np.array(x_ticks).reshape(-1, 1))]
                    plt.xticks(x_ticks, x_tick_label, fontsize=tickmark_size)
                else:
                    plt.yticks(fontsize=tickmark_size)

                if norm_f2 is not None:
                    y_ticks = g.ax.get_yticks()
                    y_tick_label = \
                        ["{0:.1e}".format(x[0]) for x in norm_f1.inverse_transform(np.array(y_ticks).reshape(-1, 1))]
                    plt.yticks(y_ticks, y_tick_label, fontsize=tickmark_size)
                else:
                    plt.yticks(fontsize=tickmark_size)

                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

            if f1 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
                    f2 in Parameters().get_numerical_features():

                if type(normalizer) == dict:
                    norm_f2 = normalizer[f2]
                    normalizer_name = get_normalizer_name(norm_f2)
                else:
                    norm_f2 = normalizer
                    normalizer_name = get_normalizer_name(norm_f2)

                g = sns.catplot(x=f1, y=f2, hue="response", kind="violin", split=True, data=df,
                                palette={'0': color_negative, '1': color_CD8}, legend=False)
                line1 = mlines.Line2D([], [], color=color_CD8, marker='s', ls='', label='CD8+')
                line0 = mlines.Line2D([], [], color=color_negative, marker='s', ls='', label='negative')
                plt.legend(title='Response type', loc=args.legend_position, handles=[line1, line0],
                           fontsize=legend_size, title_fontsize=legend_size)
                plt.xlabel(feature_dict[f1], size=axis_label_size)
                plt.ylabel(feature_dict[f2], size=axis_label_size)
                plt.xticks(fontsize=tickmark_size)
                plt.yticks(fontsize=tickmark_size)
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

            if f2 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
                    f1 in Parameters().get_numerical_features():
                g = sns.catplot(x=f2, y=f1, hue="response", kind="violin", split=True, data=df,
                                palette={'0': color_negative, '1': color_CD8}, legend=False)
                g.figure.tight_layout()

                line1 = mlines.Line2D([], [], color=color_CD8, marker='s', ls='', label='CD8+')
                line0 = mlines.Line2D([], [], color=color_negative, marker='s', ls='', label='negative')
                plt.legend(title='Response type', loc=args.legend_position, handles=[line1, line0],
                           fontsize=legend_size, title_fontsize=legend_size)
                plt.xlabel(feature_dict[f1], size=axis_label_size)
                plt.ylabel(feature_dict[f2], size=axis_label_size)
                if norm_f2 is not None:
                    y_ticks = g.ax.get_yticks()
                    y_tick_label = \
                        ["{0:.1e}".format(x[0]) for x in norm_f1.inverse_transform(np.array(y_ticks).reshape(-1, 1))]
                    plt.yticks(y_ticks, y_tick_label, fontsize=tickmark_size)
                else:
                    plt.yticks(fontsize=tickmark_size)

                plt.xticks(fontsize=tickmark_size)
                pp.savefig(g.figure)
                g.figure.clf()

            if f2 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
                    f1 in Parameters().get_categorical_features()+Parameters().get_ordinal_features():
                g = sns.catplot(x=f1, col=f2, hue="response", kind="count", data=df, legend=False,
                                palette={'0': color_negative, '1': color_CD8})
                line1 = mlines.Line2D([], [], color=color_CD8, marker='s', ls='', label='CD8+')
                line0 = mlines.Line2D([], [], color=color_negative, marker='s', ls='', label='negative')
                plt.legend(title='Response type', loc=args.legend_position, handles=[line1, line0],
                           fontsize=legend_size, title_fontsize=legend_size)
                plt.xlabel(feature_dict[f1], size=axis_label_size)
                plt.ylabel("Count", size=axis_label_size)
                plt.xticks(fontsize=tickmark_size)
                plt.yticks(fontsize=tickmark_size)
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

