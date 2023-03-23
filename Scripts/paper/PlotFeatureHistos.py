import argparse
from collections import Counter

from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from os.path import exists
from matplotlib import patches

from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-fp', '--file_prefix', type=str, default="Feature", help='Output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", help='File type for plot (png, svg or pdf')
parser.add_argument('-ds', '--datasets', type=str, nargs='+', help='Datasets used to plot feature histograms')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')
parser.add_argument('-n', '--normalizer', type=str, default='n',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None) or feature based dictionary')
parser.add_argument('-rt', '--response_types', type=str, default='', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-an', '--add_numbers', dest='add_numbers', action='store_true',
                    help='Add counts to categorical plots')
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
parser.add_argument('-lep', '--legend_position', type=str, default='best', help='Legend position in plot')
parser.add_argument('-nb', '--number_bins', type=int, default=20, help='Number of bins in histograms')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=10.0, help='Figure height in inches')
parser.add_argument('-rf', '--rotate_labels', type=str, nargs='+', help='Features with x-label rotation')
parser.add_argument('-rot', '--rotation', type=float, default=0.0, help='x-axis label rotation')
parser.add_argument('-dpi', '--resolution', type=float, default=200, help='Figure resolution in dots per inch')
parser.add_argument('-nr', '--nr_plot_rows', type=int, default=-1, help='Number of features (row) per grid plot')
parser.add_argument('-o', '--cat_order', type=str, nargs='+', help='Order of categorical features')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

if not args.features or len(args.features) == 0:
    features = Parameters().get_ml_features()
else:
    features = args.features

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

try:
    rt = ast.literal_eval(args.response_types)
except:
    print('Cannot parse dictionary {}'.format(args.response_types))


warnings.filterwarnings("ignore")

dataset_str = '_'.join(args.datasets)
result_file = os.path.join(Parameters().get_plot_dir(),
                           "{0}_Hist_{1}_{2}.txt".format(args.file_prefix, dataset_str, args.peptide_type))

with open(result_file, "w") as file:
    for arg in vars(args):
        file.write("#{0}={1}\n".format(arg, getattr(args, arg)))

if args.peptide_type == 'short':
    imm_label = 'neo-pep_imm'
    neg_label = 'neo-pep_non-imm'
else:
    imm_label = 'mut-seq_imm'
    neg_label = 'mut-seq_non-imm'

nr_plot_cols = len(args.datasets)
nr_f = len(features)
nr_plot_rows = nr_f if args.nr_plot_rows > nr_f or args.nr_plot_rows < 0 else args.nr_plot_rows

fig = plt.figure()
fig.set_figheight((nr_plot_rows+0.1)*args.figure_height)
fig.set_figwidth((nr_plot_cols+0.1)*args.figure_width)

datasets = {}
p_values = {}
for ds in args.datasets:
    patients = get_valid_patients(ds)
    data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                             mutation_types=args.mutation_types, response_types=rt[ds], immunogenic=args.immunogenic,
                             min_nr_immuno=0, max_netmhc_rank=10000)

    data, X, y = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)
    datasets[ds] = {'data': data, 'X': X, 'y': y}
    p_values[ds] = {}


def plot_feature(f_base_, f_, ds_, i_, j_, data_ds_, x_ds_, y_ds_):
    v_norm = x_ds_[f_base_]
    v = data_ds_[f_base_]
    imm_label_ = "{0}_{1}".format(ds_, imm_label)
    neg_label_ = "{0}_{1}".format(ds_, neg_label)

    print("Feature {}".format(f))
    if type(normalizer) == dict:
        norm_f = normalizer[f_base_]
    else:
        norm_f = normalizer

    if f_base_ in Parameters().get_numerical_features():
        v_norm = np.array(v_norm, dtype=float)
        x = v_norm[y_ds_ == 1]
        y = v_norm[y_ds_ == 0]
#                x = x[np.where(~np.isnan(x))]
#                y = y[np.where(~np.isnan(y))]
        tt = stats.ttest_ind(x, y, nan_policy='omit', equal_var=False, alternative='two-sided')
        p_values[ds_][f_] = tt.pvalue

        bins = np.histogram_bin_edges(v_norm, bins=args.number_bins)

        ax = plt.subplot2grid((nr_plot_rows, nr_plot_cols), (i_ % nr_plot_rows, j_))
        df = pd.DataFrame({f_: y})
        sns.histplot(
            data=df, x=f_base_, color=args.color_negative, fill=True,
            common_norm=False, alpha=.7, linewidth=0, stat="count", legend=False, bins=bins, ax=ax
        )
        ax.set_ylabel("{0} count".format(neg_label), fontsize=args.label_size, color=args.color_negative)
        ax.set_xlabel(feature_dict[f_], fontsize=args.label_size)

        ax2 = ax.twinx()
        df = pd.DataFrame({f_: x})
        sns.histplot(
            data=df, x=f_base_, color=args.color_immunogenic, fill=True, shrink=1.0,
            common_norm=False, alpha=.7, linewidth=0, stat="count", legend=False, bins=bins, ax=ax2
        )
        ax2.set_ylabel("{0} count".format(imm_label), fontsize=args.label_size, color=args.color_immunogenic)

        if norm_f is not None:
            x_ticks = ax.get_xticks()
            x_tick_label = \
                ["{0:.1e}".format(x[0]) for x in norm_f.inverse_transform(np.array(x_ticks).reshape(-1, 1))]
            ax.set_xticks(ticks=x_ticks)
            ax.set_xticklabels(labels=x_tick_label, fontsize=args.tick_size)
            ax.xaxis.set_tick_params(labelrotation=args.rotation)  # rotate anyway independent of feature
            ax2.set_xticks(ticks=x_ticks)
        else:
            ax.tick_params(axis='x', which='major', labelsize=args.tick_size)

        handles = [patches.Patch(color=args.color_immunogenic, label=imm_label_, alpha=0.7),
                   patches.Patch(color=args.color_negative, label=neg_label_, alpha=0.7)]
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), handles=handles, fontsize=args.legend_size, ncol=1)
        ax.set_title("t-test p-value = {0:.2e}".format(tt.pvalue), fontsize=args.legend_size)
        plt.tight_layout()

    if f_base_ in Parameters().get_ordinal_features():
        counts1 = Counter(v[y_ds_ == 1])
        counts0 = Counter(v[y_ds_ == 0])
        x = np.array([])
        y = np.array([])
        for key in v.unique():
            if key in counts1:
                v1 = counts1[key]
            else:
                v1 = 0
            if key in counts0:
                v0 = counts0[key]
            else:
                v0 = 0

            x = np.append(x, v1)
            y = np.append(y, v0)

        idx = (x > 0) | (y > 0)
        cont_table = np.array([x[idx], y[idx]])
        chi2, p, dof, ex = stats.chi2_contingency(cont_table)
        p_values[ds_][f_] = p

        lbls = np.array(v.unique(), dtype=int)

        if args.cat_order is None or len(args.cat_order) != len(lbls):
            lbls, x, y = zip(*sorted(zip(lbls, x, y), key=lambda triple: triple[0]))
        else:
            wrong_lbls = [lbl for lbl in args.cat_order if lbl not in lbls]
            if len(wrong_lbls) > 0:
                print("Wrong labels in cat_order: "+",".join(wrong_lbls)+". No sorting")
                lbls, x, y = zip(*sorted(zip(lbls, x, y), key=lambda triple: triple[0]))
            else:
                lst = zip(lbls, x, y)
                lst.sort(key=lambda i: args.cat_order.index(i[0]))
                lbls, x, y = zip(*lst)

        ax = plt.subplot2grid((nr_plot_rows, nr_plot_cols), (i_ % nr_plot_rows, j_))
        ax.bar(x=lbls, height=y, color=args.color_negative, label=neg_label_, alpha=0.7)
        ax.set_xlabel(feature_dict[f_], size=args.label_size)
        ax.set_ylabel("{0} count".format(neg_label), fontsize=args.label_size, color=args.color_negative)
        ax.tick_params(axis='x', which='major', labelsize=args.tick_size, labelrotation=rotation)
        ax.tick_params(axis='y', which='major', labelsize=args.tick_size)
        ax2 = ax.twinx()
        ax2.bar(x=np.array(lbls)+0.1, height=x, color=args.color_immunogenic, label=imm_label_, alpha=0.7, width=0.8)
        ax2.set_ylabel("{0} count".format(imm_label), fontsize=args.label_size, color=args.color_immunogenic)
        ax2.tick_params(axis='y', which='major', labelsize=args.tick_size)

        handles = [patches.Patch(color=args.color_immunogenic, label=imm_label_, alpha=0.7),
                   patches.Patch(color=args.color_negative, label=neg_label_, alpha=0.7)]
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), handles=handles, fontsize=args.legend_size, ncol=1)
        ax.set_title("Chi2 p-value = {0:.2e}".format(p), fontsize=args.legend_size)
        plt.tight_layout()

    if f_base_ in Parameters().get_categorical_features():
        counts1 = Counter(v[y_ds_ == 1])
        counts0 = Counter(v[y_ds_ == 0])
        x = []
        y = []
        v_u = v.unique()
        for key in v_u:
            if key in counts1:
                v1 = counts1[key]
            else:
                v1 = 0
            if key in counts0:
                v0 = counts0[key]
            else:
                v0 = 0

            x.append(v1)
            y.append(v0)

        cont_table = np.array([x, y])
        chi2, p, dof, ex = stats.chi2_contingency(cont_table)
        p_values[ds_][f_] = p

        lbls_pos = np.array(np.arange(len(v_u)))

        v_u = np.array(v.unique(), dtype=str)
        if args.cat_order is None or len(args.cat_order) != len(v_u):
            v_u, x, y = zip(*sorted(zip(v_u, x, y), key=lambda triple: triple[0]))
        else:
            wrong_lbls = [lbl for lbl in args.cat_order if lbl not in v_u]
            if len(wrong_lbls) > 0:
                print("Wrong labels in cat_order: "+",".join(wrong_lbls)+". No sorting")
                v_u, x, y = zip(*sorted(zip(v_u, x, y), key=lambda triple: triple[0]))
            else:
                lst = list(zip(v_u, x, y))
                lst.sort(key=lambda i: args.cat_order.index(i[0]))
                v_u, x, y = zip(*lst)

        ax = plt.subplot2grid((nr_plot_rows, nr_plot_cols), (i_ % nr_plot_rows, j_))
        ax.bar(x=lbls_pos, height=y, color=args.color_negative, label=neg_label_, alpha=0.7)
        ax.set_xlabel(feature_dict[f_], size=args.label_size)
        ax.set_ylabel("{0} count".format(neg_label), fontsize=args.label_size, color=args.color_negative)
        ax.set_xticks(ticks=lbls_pos)
        ax.set_xticklabels(labels=v_u, fontsize=args.tick_size)
        ax.xaxis.set_tick_params(labelrotation=rotation)
        ax.tick_params(axis='y', which='major', labelsize=args.tick_size)
        ax2 = ax.twinx()
        ax2.bar(x=lbls_pos+0.1, height=x, color=args.color_immunogenic, label=imm_label_, alpha=0.7, width=0.8)
        ax2.set_ylabel("{0} count".format(imm_label), fontsize=args.label_size, color=args.color_immunogenic)
        ax2.tick_params(axis='y', which='major', labelsize=args.tick_size)

        handles = [patches.Patch(color=args.color_immunogenic, label=imm_label_, alpha=0.7),
                   patches.Patch(color=args.color_negative, label=neg_label_, alpha=0.7)]
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.50), handles=handles, fontsize=args.legend_size, ncol=1)
        ax.set_title("Chi2 p-value = {0:.2e}".format(p), fontsize=args.legend_size)
        plt.tight_layout()

        if args.add_numbers:
            max_v = np.max(np.c_[x, y], axis=1)
            labels = ['{0:.0f}/{1:.0f}'.format(v_x, v_y) for v_x, v_y in zip(x, y)]
            for i in range(len(v_u)):
                plt.annotate(labels[i], xy=(i, max_v[i]), xytext=(0, 5), textcoords="offset points", ha="center")


def filter_by_len(feature, data_ds_, X_ds_, y_ds_):
    seq_len = int(feature.split("_")[-1])
    len_idx = data_ds_['seq_len'] == seq_len
    return data_ds_[len_idx], X_ds_[len_idx], y_ds_[len_idx]


plot_idx = 1
for i, f in enumerate(features):

    if f in args.rotate_labels:
        rotation = args.rotation
    else:
        rotation = 0.0

    for j, ds in enumerate(args.datasets):
        data_ds = datasets[ds]['data']
        X_ds = datasets[ds]['X']
        y_ds = datasets[ds]['y']

        if f not in p_values[ds]:
            p_values[ds][f] = np.nan

        if f.startswith('pep_mut_start_'):
            data_ds, X_ds, y_ds = filter_by_len(f, data_ds, X_ds, y_ds)
            f_base = 'pep_mut_start'
        else:
            f_base = f

        if f_base not in X_ds:
            continue

        plot_feature(f_base, f, ds, i, j, data_ds, X_ds, y_ds)

    if i%nr_plot_rows == nr_plot_rows-1 or i == len(features)-1:
        fig.tight_layout()
        if nr_plot_rows == 1:
            file_name = "{0}_PValues_{1}_{2}_{3}_{4}.{5}".\
                format(args.file_prefix, dataset_str, args.peptide_type, f, plot_idx, args.file_type)
        else:
            file_name = "{0}_PValues_{1}_{2}_{3}.{4}".\
                format(args.file_prefix, dataset_str, args.peptide_type, plot_idx, args.file_type)

        plt.tight_layout()
        plot_file = os.path.join(Parameters().get_plot_dir(), file_name)
        plt.savefig(plot_file, bbox_inches='tight')
        plot_idx += 1

        if i < len(features)-1:
            plt.close()
            fig = plt.figure()
            fig.set_figheight((nr_plot_rows+0.1)*args.figure_height)
            fig.set_figwidth((nr_plot_cols+0.1)*args.figure_width)

pv_df = pd.DataFrame(p_values)
pv_df.to_csv(result_file, sep='\t', header=True, index=True, mode='a')
