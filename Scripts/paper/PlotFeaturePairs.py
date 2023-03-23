import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import patches

from DataWrangling.DataLoader import DataLoader
from DataWrangling.Transform_Data import DataTransformer
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot correlation between features')

parser.add_argument('-fp', '--file_prefix', type=str, default="Feature_pair", help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", help='File type for plot (png, svg or pdf')
parser.add_argument('-ds', '--dataset', type=str, default='NCI', help='Dataset used to plot feature correlation')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-fpair', '--feature_pairs', type=str, nargs='+', help='Features pair for pair plots')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='n',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None or dictionary)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-reg', '--regression_line', dest='regression_line', action='store_true',
                    help='draw regression line')
parser.add_argument('-kd', '--kd_plot', dest='kd_plot', action='store_true', help='Draw density contour curves')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-s', '--nr_samples', type=int, default=10000, help='Nr of Sampled rows from dataframe for plots')
parser.add_argument('-lep', '--legend_position', type=str, default='best', help='Legend position in plot')
parser.add_argument('-ci', '--color_immunogenic', type=str, default='darkorange', help='Color of immunogenic peptides')
parser.add_argument('-cn', '--color_negative', type=str, default='royalblue', help='Color of negative peptides')
parser.add_argument('-las', '--label_size', type=str, default='x-large',
                    help='Axis label size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-tis', '--tick_size', type=str, default='large',
                    help='Tick size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-les', '--legend_size', type=str, default='large',
                    help='Legend size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-fiw', '--figure_width', type=float, default=25.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=25.0, help='Figure height in inches')
parser.add_argument('-rot', '--rotation', type=float, default=0.0, help='x-axis label rotation')
parser.add_argument('-rf', '--rotate_labels', type=str, nargs='+', help='Features with x-label rotation')
parser.add_argument('-dpi', '--resolution', type=float, default=200, help='Figure resolution in dots per inch')
parser.add_argument('-hl', '--horizontal_line', type=float, default=None, help='Horizontal line at y in figure')
parser.add_argument('-vl', '--vertical_line', type=float, default=None, help='Vertical line at x in figure')
parser.add_argument('-ps', '--point_size', type=float, default=1, help='Size of points in scatterplot')
parser.add_argument('-o', '--cat_order', type=str, nargs='+', help='Order of categorical features')
parser.add_argument('-xr', '--x_range', type=str, default='', help='Range of x value in plot')
parser.add_argument('-yr', '--y_range', type=str, default='', help='Range of y value in plot')

args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

patients = get_valid_patients(args.dataset)

normalizer = get_normalizer(args.normalizer)

if args.peptide_type == 'short':
    imm_label = args.dataset+'_neo-pep_imm'
    neg_label = args.dataset+'_neo-pep_non-imm'
    y_imm_label = 'neo-pep_imm count'
    y_neg_label = 'neo-pep_non-imm count'
else:
    imm_label = args.dataset+'_mut-seq_imm'
    neg_label = args.dataset+'_mut-seq_non-imm'
    y_imm_label = 'mut-seq_imm count'
    y_neg_label = 'mut-seq_non-imm count'

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
patients = get_valid_patients(args.dataset)

result_file = os.path.join(Parameters().get_plot_dir(),
                           "{0}_Hist_{1}_{2}.txt".format(args.file_prefix, args.dataset, args.peptide_type))

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0)

data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)

X_1 = X_train.loc[y_train == 1, :]
X_1['response'] = "1"
X_1['Frequency'] = 1
X_0 = X_train.loc[y_train == 0, :]
X_0['response'] = "0"
X_0['Frequency'] = 1
if X_0.shape[0] > args.nr_samples:
    X_0 = X_0.sample(n=args.nr_samples)
df = pd.concat([X_1, X_0])
all_num_cols = Parameters().get_numerical_features()
num_cols = [c for c in df.columns if c in all_num_cols]
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors='coerce')
df_small = df.head(1000).sort_values(by=['response'])

R_values = {}

for fp in args.feature_pairs:

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

    plot_file = os.path.join(Parameters().get_plot_dir(),
                             "{0}_{1}_{2}_{3}.{4}".format(args.file_prefix, args.dataset, f1, f2, args.file_type))

    if f1 in Parameters().get_numerical_features() and f2 in Parameters().get_numerical_features():
        point_size = df.apply(lambda r: args.point_size*(1+int(r['response'])*2), axis=1).astype(np.float)
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.lmplot(data=df, x=f1, y=f2, hue="response", line_kws={'alpha': 0.7}, hue_order=['0', '1'],
                       scatter_kws={'alpha': 0.5, 's': 1}, legend=True, fit_reg=args.regression_line,
                       palette={'0': args.color_negative, '1': args.color_immunogenic})

        x_min = df[f1].min()-(df[f1].max()-df[f1].min())/20
        y_min = df[f2].min()-(df[f2].max()-df[f2].min())/20
        g.set(xlim=(x_min, None))
        g.set(ylim=(y_min, None))
        plt.scatter(df.loc[df['response'] == '1',  f1], df.loc[df['response'] == '1',  f2], s=10, c=args.color_immunogenic)

        if args.kd_plot:
            g = g.map_dataframe(sns.kdeplot, x=f1, y=f2, hue="response", hue_order=['0', '1'], alpha=0.4, fill=True)

        plt.xlabel(feature_dict[f1], size=args.label_size)
        plt.ylabel(feature_dict[f2], size=args.label_size)

        sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=1, labels=(neg_label, imm_label),
                        title="", frameon=False, fontsize=args.legend_size, markerscale=10.0)

        if norm_f1 is not None:
            x_ticks = g.ax.get_xticks()
            x_ticks = [x for x in x_ticks if x >= x_min]
            x_tick_label = \
                ["{0:.1e}".format(x[0]) for x in norm_f1.inverse_transform(np.array(x_ticks).reshape(-1, 1))]
            plt.xticks(x_ticks, x_tick_label, fontsize=args.tick_size, rotation=args.rotation)
        else:
            plt.xticks(fontsize=args.tick_size, rotation=args.rotation)

        if norm_f2 is not None:
            y_ticks = g.ax.get_yticks()
            y_ticks = [y for y in y_ticks if y >= y_min]
            y_tick_label = \
                ["{0:.1e}".format(y[0]) for y in norm_f2.inverse_transform(np.array(y_ticks).reshape(-1, 1))]
            plt.yticks(y_ticks, y_tick_label, fontsize=args.tick_size)
        else:
            plt.yticks(fontsize=args.tick_size)

        if args.horizontal_line is not None:
            if norm_f2:
                hv = norm_f2.transform(np.array([args.horizontal_line]).reshape(-1, 1))[0]
            else:
                hv = args.horizontal_line
            plt.axhline(y=hv, color="red", linestyle="--", linewidth=2)

        if args.vertical_line is not None:
            if norm_f1:
                vv = norm_f1.transform(np.array([args.vertical_line]).reshape(-1, 1))[0]
            else:
                vv = args.vertical_line
            plt.axvline(x=vv, color="red", linestyle="--", linewidth=2)

        corr = df[[f1, f2]].corr().loc[f1, f2]
        R_values[fp] = corr
        print("Pearson Corr between {0} and {1}: {2:.3f}".format(f1, f2, corr))

        plt.tight_layout()
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

    if f1 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
            f2 in Parameters().get_numerical_features():
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.violinplot(x=f1, y=f2, hue="response", scale="width", split=True, data=df, inner='quartile',
                           order=args.cat_order, palette={'0': args.color_negative, '1': args.color_immunogenic},
                           legend=False)

        if args.horizontal_line is not None:
            plt.plot([df[f2].min(), df[f2].max()], [args.horizontal_line, args.horizontal_line], color='red',
                     linestyle='dashed', linewidth=2)

        line1 = mlines.Line2D([], [], color=args.color_immunogenic, marker='s', ls='', label=imm_label)
        line0 = mlines.Line2D([], [], color=args.color_negative, marker='s', ls='', label=neg_label)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, 2.50), handles=[line1, line0], fontsize=args.legend_size,
                   ncol=1)
        plt.xlabel(feature_dict[f1], size=args.label_size)
        plt.ylabel(feature_dict[f2], size=args.label_size)
        if norm_f2 is not None:
            y_ticks = g.get_yticks()
            y_tick_label = \
                ["{0:.1e}".format(x[0]) for x in norm_f2.inverse_transform(np.array(y_ticks).reshape(-1, 1))]
            plt.yticks(y_ticks, y_tick_label, fontsize=args.tick_size)
        else:
            plt.yticks(fontsize=args.tick_size)

        if f1 in args.rotate_labels:
            rotation = args.rotation
        else:
            rotation = 0.0
        plt.xticks(fontsize=args.tick_size, rotation=rotation)
        plt.tight_layout()
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

    if f2 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
            f1 in Parameters().get_numerical_features():
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        g = sns.violinplot(x=f1, y=f2, hue="response", scale="width", split=True, data=df, inner='quartile',
                           order=args.cat_order, palette={'0': args.color_negative, '1': args.color_immunogenic},
                           legend=False)

        if args.vertical_line is not None:
            plt.plot([args.vertical_line, args.vertical_line], [df[f1].min(), df[f1].max()], color='red',
                     linestyle='dashed', linewidth=2)

        line1 = mlines.Line2D([], [], color=args.color_immunogenic, marker='s', ls='', label=imm_label)
        line0 = mlines.Line2D([], [], color=args.color_negative, marker='s', ls='', label=neg_label)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), handles=[line1, line0], fontsize=args.legend_size,
                  ncol=1)
        plt.xlabel(feature_dict[f1], size=args.label_size)
        plt.ylabel(feature_dict[f2], size=args.label_size)

        if f1 in args.rotate_labels:
            rotation = args.rotation
        else:
            rotation = 0.0

        if norm_f1 is not None:
            x_ticks = g.get_xticks()
            x_tick_label = \
                ["{0:.1e}".format(x[0]) for x in norm_f2.inverse_transform(np.array(x_ticks).reshape(-1, 1))]
            plt.xticks(x_ticks, x_tick_label, fontsize=args.tick_size, rotation=rotation)
        else:
            plt.xticks(fontsize=args.tick_size, rotation=rotation)

        plt.yticks(fontsize=args.tick_size)
        plt.tight_layout()
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

    if f2 in Parameters().get_categorical_features()+Parameters().get_ordinal_features() and \
            f1 in Parameters().get_categorical_features()+Parameters().get_ordinal_features():

        v1_u = df[f1].unique()
        fig = plt.figure()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(len(v1_u) * args.figure_width)

        if f2 in args.rotate_labels:
            rotation = args.rotation
        else:
            rotation = 0.0

        for fig_i, v_u in enumerate(v1_u):
            df_ = df.loc[df[f1] == v_u, ]
            c_f2_imm = df_.loc[df_["response"] == '1', f2].value_counts()
            c_f2_neg = df_.loc[df_["response"] == '0', f2].value_counts()
            v_imm = []
            v_neg = []

            v2_u = df[f2].unique()
            for l_u in v2_u:
                if l_u in c_f2_imm.index:
                    v_imm.append(c_f2_imm[l_u])
                else:
                    v_imm.append(0)
                if l_u in c_f2_neg.index:
                    v_neg.append(c_f2_neg[l_u])
                else:
                    v_neg.append(0)

            if f2 in Parameters().get_ordinal_features():
                v2_u = np.array(v2_u, dtype=int)
            else:
                v2_u = np.array(v2_u, dtype=str)

            if args.cat_order is None or len(args.cat_order) != len(v2_u):
                v2_u, v_imm, v_neg = zip(*sorted(zip(v1_u, v_imm, v_neg), key=lambda triple: triple[0]))
            else:
                wrong_lbls = [lbl for lbl in args.cat_order if lbl not in v2_u]
                if len(wrong_lbls) > 0:
                    print("Wrong labels in cat_order: "+",".join(wrong_lbls)+". No sorting")
                    v2_u, v_imm, v_neg = zip(*sorted(zip(v2_u, v_imm, v_neg), key=lambda triple: triple[0]))
                else:
                    lst = list(zip(v2_u, v_imm, v_neg))
                    lst.sort(key=lambda i: args.cat_order.index(i[0]))
                    v2_u, v_imm, v_neg = zip(*lst)

            ax = plt.subplot2grid((1, len(v1_u)), (0, fig_i))
            lbls_pos = np.arange(len(v2_u))
            ax.bar(x=lbls_pos, height=v_neg, color=args.color_negative, alpha=0.7)
            ax.set_xticks(ticks=lbls_pos)
            ax.set_xticklabels(labels=v2_u, fontsize=args.tick_size, rotation=rotation)
            ax.tick_params(axis='y', which='major', labelsize=args.tick_size)
            ax.set_xlabel(feature_dict[f2], size=args.label_size)
            ax.set_ylabel(y_neg_label, size=args.label_size, color=args.color_negative)

            ax2 = ax.twinx()
            ax2.bar(x=np.array(lbls_pos)+0.1, height=v_imm, color=args.color_immunogenic, alpha=0.7)
            ax2.tick_params(axis='y', which='major', labelsize=args.tick_size)
            ax2.set_ylabel(y_imm_label, size=args.label_size, color=args.color_immunogenic)

            ax.set_title("{0} = {1}".format(feature_dict[f1], v_u), fontsize=args.legend_size)
            handles = [patches.Patch(color=args.color_immunogenic, label=imm_label, alpha=0.7),
                       patches.Patch(color=args.color_negative, label=neg_label, alpha=0.7)]
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.45), handles=handles, fontsize=args.legend_size,
                      ncol=1)
            plt.tight_layout()

            # if args.add_numbers:
            #     max_v = np.max(np.c_[c_f2_imm.to_numpy(), c_f2_neg.to_numpy()], axis=1)
            #     labels = ['{0:.0f}/{1:.0f}'.format(v_x, v_y) for v_x, v_y in zip(x, y)]
            #     for i in range(len(lbls_pos)):
            #         plt.annotate(labels[i], xy=(lbls_pos[i], max_v[i]), xytext=(0, 5), textcoords="offset points",
            #                      ha="center")
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

with open(result_file, 'a') as file:
    for fp in R_values:
        file.write("Pearson Corr between {0}: {1:.3f}".format(fp, R_values[fp]))
