"""
Plot histograms for the different numerical and categorical features.
"""
from collections import Counter

import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
import pandas as pd

from Utils.Util_fct import *
from Utils.DataManager import DataManager
from DataWrangling.DataTransformer import DataTransformer

parser = argparse.ArgumentParser(description='Plot histograms for immunogenic and non immunogenic feature values')

parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-ds', '--dataset', type=str, action='append', choices=GlobalParameters.datasets,
                    help='Dataset, one of [NCI, NCI_train, NCI_test, TESLA, HiTIDE]')
parser.add_argument('-f', '--feature', type=str,
                    choices=GlobalParameters.features_neopep+GlobalParameters.features_mutation+['pep_mut_start_10'],
                    help='Histograms of the feature\'s immunogenic and non-immunogenic values are plotted')
parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-an', '--add_numbers', dest='add_numbers', action='store_true',
                    help='Add counts to categorical plots')
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
parser.add_argument('-rot', '--rotation', type=float, default=0.0, help='x-axis label rotation')
parser.add_argument('-dpi', '--resolution', type=float, default=200, help='Figure resolution in dots per inch')
parser.add_argument('-nr', '--nr_plot_rows', type=int, default=-1, help='Number of features (row) per grid plot')
parser.add_argument('-o', '--cat_order', type=str, nargs='+', help='Order of categorical features')
parser.add_argument('-lbls', '--x_labels', type=str, nargs='+', help='Rename x-labels with these values')
parser.add_argument('-log', '--log_scale', dest='log_scale', action='store_true', help='Plots counts on log-scale')
parser.add_argument('-nt', '--not_tested', dest='not_tested', action='store_true',
                    help='not tested peptides to be included for TESLA and HiTIDE')
parser.add_argument('-ena', '--exclude_nan', dest='exclude_nan', action='store_true', help='Exclude nan category from plot')
parser.add_argument('-fw', '--frame_width', type=float, default=0.1, help='Width of plot frame')


def plot_feature(peptide_type: str, f_base_, f_, ds_, j_, data_ds_, x_ds_, y_ds_, imm_label_, neg_label_):
    if f_base_ == 'pep_mut_start':
        v_norm = data_ds_[f_base_]
    else:
        v_norm = x_ds_[f_base_]
    v = data_ds_[f_base_]
    imm_label_ds = "{0}_{1}".format(ds_, imm_label_)
    neg_label_ds = "{0}_{1}".format(ds_, neg_label_)
    if args.log_scale:
        units = "log10 count"
    else:
        units = "count"

    pvalue_ = np.nan
    type_ = get_processed_types(peptide_type=peptide_type, objective='plot')[f_base_]
    norm_f = DataTransformer.get_normalizer('plot')[f_base_]

    if is_cont_type(type_):
        v_norm = np.array(v_norm, dtype=float)
        x = v_norm[y_ds_ == 1]
        y = v_norm[y_ds_ == 0]
        x = x[np.where(~np.isnan(x))]
        y = y[np.where(~np.isnan(y))]
        tt = stats.ttest_ind(x, y, nan_policy='omit', equal_var=False, alternative='two-sided')
        pvalue_ = tt.pvalue

        if args.log_scale:
            y = np.add(y, 1)
            x = np.add(x, 1)

        bins = np.histogram_bin_edges(np.append(x, y), bins=args.number_bins)

        ax = plt.subplot2grid((1, nr_plot_cols), (0, j_))
        df = pd.DataFrame({f_: y})
        gl = sns.histplot(
            data=df, x=f_base_, color=GlobalParameters.color_negative, fill=True,
            common_norm=False, alpha=.7, linewidth=0, stat="count", legend=False, bins=bins, ax=ax
        )
        if args.log_scale:
            gl.set_yscale("log")
        ax.set_ylabel("{0} {1}".format(neg_label_, units), fontsize=args.label_size,
                      color=GlobalParameters.color_negative)
        ax.set_xlabel(GlobalParameters.plot_feature_names[f_], fontsize=args.label_size)
        [x.set_linewidth(args.frame_width) for x in ax.spines.values()]

        ax2 = ax.twinx()
        df = pd.DataFrame({f_: x})
        gr = sns.histplot(
            data=df, x=f_base_, color=GlobalParameters.color_immunogenic, fill=True,
            common_norm=False, alpha=.7, linewidth=0, stat="count", legend=False, bins=bins, ax=ax2
        )
        if args.log_scale:
            gr.set_yscale("log")
        ax2.set_ylabel("{0} {1}".format(imm_label_, units), fontsize=args.label_size,
                       color=GlobalParameters.color_immunogenic)
        [x.set_linewidth(args.frame_width) for x in ax2.spines.values()]

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

        handles = [patches.Patch(color=GlobalParameters.color_immunogenic, label=imm_label_ds, alpha=0.7),
                   patches.Patch(color=GlobalParameters.color_negative, label=neg_label_ds, alpha=0.7)]
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), handles=handles, fontsize=args.legend_size, ncol=1)
        ax.set_title("t-test p-value = {0:.2e}".format(tt.pvalue), fontsize=args.legend_size)
        plt.tight_layout()

    elif is_discrete_ordered_type(type_):
        counts1 = Counter(v[y_ds_ == 1])
        counts0 = Counter(v[y_ds_ == 0])
        x = np.array([])
        y = np.array([])
        lbls = np.sort(v.unique())
        if args.exclude_nan:
            lbls = lbls[lbls != np.nan]
        for key in lbls:
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
        pvalue_ = p

        if args.log_scale:
            y = np.log10(np.add(y, 1))
            x = np.log10(np.add(x, 1))

        if args.x_labels and len(args.x_labels) >= len(lbls):
            lbls = args.x_labels[0:len(lbls)]

        lbls_pos = np.array(np.arange(len(lbls)))

        ax = plt.subplot2grid((1, nr_plot_cols), (0, j_))
        ax.bar(x=lbls, height=y, color=GlobalParameters.color_negative, label=neg_label_, alpha=0.7)
        ax.set_xlabel(GlobalParameters.plot_feature_names[f_], size=args.label_size)
        ax.set_ylabel("{0} {1}".format(neg_label_, units), fontsize=args.label_size, color=GlobalParameters.color_negative)
        ax.tick_params(axis='x', which='major', labelsize=args.tick_size, labelrotation=args.rotation)
        ax.tick_params(axis='y', which='major', labelsize=args.tick_size)
        [x.set_linewidth(args.frame_width) for x in ax.spines.values()]
        ax2 = ax.twinx()
        ax2.bar(x=np.array(lbls_pos)+0.1, height=x, color=GlobalParameters.color_immunogenic, label=imm_label_, alpha=0.7, width=0.8)
        ax2.set_ylabel("{0} {1}".format(imm_label_, units), fontsize=args.label_size, color=GlobalParameters.color_immunogenic)
        ax2.tick_params(axis='y', which='major', labelsize=args.tick_size)
        [x.set_linewidth(args.frame_width) for x in ax2.spines.values()]

        handles = [patches.Patch(color=GlobalParameters.color_immunogenic, label=imm_label_ds, alpha=0.7),
                   patches.Patch(color=GlobalParameters.color_negative, label=neg_label_ds, alpha=0.7)]
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), handles=handles, fontsize=args.legend_size, ncol=1)
        ax.set_title("Chi2 p-value = {0:.2e}".format(p), fontsize=args.legend_size)
        plt.tight_layout()

    if is_cat_type(type_):
        counts1 = Counter(v[y_ds_ == 1])
        counts0 = Counter(v[y_ds_ == 0])
        x = []
        y = []
        v_u = v.unique()
        if args.exclude_nan:
            v_u = v_u[~v_u.isna()]
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
        pvalue_ = p

        lbls_pos = np.array(np.arange(len(v_u)))

        v_u = np.array(v_u, dtype=str)
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

        if args.log_scale:
            y = np.log10(np.add(y, 1))
            x = np.log10(np.add(x, 1))

        if args.x_labels and len(args.x_labels) >= len(v_u):
            v_u = args.x_labels[0:len(v_u)]

        ax = plt.subplot2grid((1, nr_plot_cols), (0, j_))
        ax.bar(x=lbls_pos, height=y, color=GlobalParameters.color_negative, label=neg_label_, alpha=0.7)
        ax.set_xlabel(GlobalParameters.plot_feature_names[f_], size=args.label_size)
        ax.set_ylabel("{0} {1}".format(neg_label_, units), fontsize=args.label_size, color=GlobalParameters.color_negative)
        ax.set_xticks(ticks=lbls_pos)
        ax.set_xticklabels(labels=v_u, fontsize=args.tick_size)
        ax.xaxis.set_tick_params(labelrotation=args.rotation)
        ax.tick_params(axis='y', which='major', labelsize=args.tick_size)
        [x.set_linewidth(args.frame_width) for x in ax.spines.values()]
        ax2 = ax.twinx()
        ax2.bar(x=lbls_pos+0.1, height=x, color=GlobalParameters.color_immunogenic, label=imm_label_, alpha=0.7, width=0.8)
        ax2.set_ylabel("{0} {1}".format(imm_label_, units), fontsize=args.label_size, color=GlobalParameters.color_immunogenic)
        ax2.tick_params(axis='y', which='major', labelsize=args.tick_size)
        [x.set_linewidth(args.frame_width) for x in ax2.spines.values()]

        handles = [patches.Patch(color=GlobalParameters.color_immunogenic, label=imm_label_ds, alpha=0.7),
                   patches.Patch(color=GlobalParameters.color_negative, label=neg_label_ds, alpha=0.7)]
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.50), handles=handles, fontsize=args.legend_size, ncol=1)
        ax.set_title("Chi2 p-value = {0:.2e}".format(p), fontsize=args.legend_size)
        plt.tight_layout()

        if args.add_numbers:
            max_v = np.max(np.c_[x, y], axis=1)
            labels = ['{0:.0f}/{1:.0f}'.format(v_x, v_y) for v_x, v_y in zip(x, y)]
            for i in range(len(v_u)):
                plt.annotate(labels[i], xy=(i, max_v[i]), xytext=(0, 5), textcoords="offset points", ha="center")

    return pvalue_


def filter_by_len(feature, data_ds_, X_ds_, y_ds_):
    seq_len = int(feature.split("_")[-1])
    len_idx = data_ds_['seq_len'] == seq_len
    return data_ds_[len_idx], X_ds_[len_idx], y_ds_[len_idx]


if __name__ == "__main__":
    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    if args.peptide_type == 'neopep':
        imm_label = 'neo-pep_imm'
        neg_label = 'neo-pep_non-imm'
    else:
        imm_label = 'mut-seq_imm'
        neg_label = 'mut-seq_non-imm'

    nr_plot_cols = len(args.dataset)

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    fig = plt.figure()
    fig.set_figheight(1.1*args.figure_height)
    fig.set_figwidth((nr_plot_cols+0.1)*args.figure_width)

    dataset_str = "_".join(args.dataset)
    p_values = {}
    for j, ds in enumerate(args.dataset):
        if ds.startswith('NCI'):
            response_types = ['CD8', 'negative']
        else:
            response_types = ['CD8', 'negative', 'not_tested'] if args.not_tested else ['CD8', 'negative']

        data_ds, X_ds, y_ds = DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='plot',
                                                                dataset=ds, response_types=response_types, sample=False)

        if args.feature.startswith('pep_mut_start_'):
            data_ds, X_ds, y_ds = filter_by_len(args.feature, data_ds, X_ds, y_ds)
            f_base = 'pep_mut_start'
        else:
            f_base = args.feature

        p_values[ds] = plot_feature(args.peptide_type, f_base, args.feature, ds, j, data_ds, X_ds, y_ds, imm_label, neg_label)

    plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
    plt.tight_layout()
    plt.savefig(plot_file, bbox_inches='tight', transparent=args.file_type == 'pdf')

    result_file = os.path.join(GlobalParameters.plot_dir, "{0}.txt".format(args.file_name))
    with open(result_file, "w") as file:
        for arg in vars(args):
            file.write("#{0}={1}\n".format(arg, getattr(args, arg)))
        file.write("Dataset\tFeature\tpValue\n")
        for ds in p_values:
            file.write("{0}\t{1}\t{2}\n".format(ds, args.feature, p_values[ds]))

