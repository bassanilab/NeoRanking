"""
PCA scatterplot of rank_score vectors for different classifiers
"""
import argparse
import glob
import ast

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patches
import seaborn as sns
from sklearn.decomposition import PCA
import re

from Utils.GlobalParameters import GlobalParameters
from Utils.Util_fct import *
from Classifier.ClassifierResultParser import ClassifierResults, SimpleRankingResults, GartnerResultsNeopep, \
    GartnerResultsMutation
from Utils.DataManager import DataManager


parser = argparse.ArgumentParser(description='PCA scatterplot of rank_score vectors for different classifiers')
parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory holding classifier result files')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-ds', '--dataset', type=str, choices=GlobalParameters.datasets,
                    help='Dataset, one of [NCI, NCI_train, NCI_test, TESLA, HiTIDE]')
parser.add_argument('-re', '--clf_result_files_re', type=str, nargs='+',
                    help='Comma separated list of clf result file regular expressions')
parser.add_argument('-g', '--patient_group_plot_names', type=str, default='', help='Patient groups plot names')
parser.add_argument('-cln', '--classifier_plot_names', type=str, default='', help='Classifier order in plots')
parser.add_argument('-o', '--plot_order', type=str, help='Order of classifier plot names')
parser.add_argument('-a', '--alpha', type=float, default=0.02, help='Coefficient alpha in score function')
parser.add_argument('-ga', '--gartner', dest='gartner', action='store_true',
                    help='Include ranking from Gartner et al. in comparison (only for NCI_test)')
parser.add_argument('-sr', '--simple_ranking', dest='simple_ranking', action='store_true',
                    help='Include simple ranking with MixMHCpred and NetMHCpan in comparison')
parser.add_argument('-rot', '--rotation', type=float, default=30.0, help='x-axis label rotation')
parser.add_argument('-las', '--label_size', type=float, default=25.0, help='Axis label size')
parser.add_argument('-xl', '--xlabel', type=str, default="", help='x-label')
parser.add_argument('-tis', '--tick_size', type=float, default=20.0, help='Axis tick size')
parser.add_argument('-tts', '--title_size', type=float, default=20.0, help='title size')
parser.add_argument('-res', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=6.00, help='Figure height in inches')
parser.add_argument('-les', '--legend_size', type=float, default=15, help='Legend size in float')
parser.add_argument('-oc', '--one_color', type=str, default=None, help='all boxes in same color')
parser.add_argument('-bar', '--bar_plot', dest='bar_plot', action='store_true', help='bars instead of boxes')
parser.add_argument('-ttp', '--title_prefix', type=str, default='', help='prefix for plot title')
parser.add_argument('-ylim', '--rank_score_lim', type=str, default='', help='plot limits for rankscore')
parser.add_argument('-cm', '--color_map', type=str, default='', help='color map for classifiers')
parser.add_argument('-fw', '--frame_width', type=float, default=0.1, help='Width of plot frame')

if __name__ == "__main__":

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    data = DataManager.filter_selected_data(peptide_type=args.peptide_type, dataset=args.dataset)

    included_patients = set(data['patient'].unique())
    if args.gartner:
        data, X, y = DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='sel',
                                                       dataset='NCI_test', sample=False)
        included_patients = data['patient'].unique().to_numpy()

    parse_gartner = True
    parse_simple = True

    vector_df = None

    color_map = None
    ylim_dict = {}
    if args.rank_score_lim != '':
        ylim_dict = ast.literal_eval(args.rank_score_lim)
    if args.color_map != '':
        color_map = {}
        cm = ast.literal_eval(args.color_map)
        for key in cm:
            color_map[key] = sns.color_palette()[cm[key]]

    plt_name_dict = ast.literal_eval(args.classifier_plot_names)

    topN_palette = sns.color_palette("Greens_d", 3)

    if args.plot_order:
        plot_order = args.plot_order.split(',')
    else:
        plot_order = None

    for j, regexp in enumerate(args.clf_result_files_re):
        clf_result_files = glob.glob(os.path.join(GlobalParameters.classifier_result_dir, args.sub_dir, regexp))

        for i, clf_result_file in enumerate(clf_result_files):
            if os.path.getsize(clf_result_file) > 0:
                print(clf_result_file)
                clf_results = ClassifierResults(args.peptide_type, clf_result_file, i, plt_name_dict[j], data,
                                                alpha=args.alpha, included_patients_=included_patients)
                vector_df = clf_results.add_to_vector_df(vector_df)

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    pca = PCA(n_components=2)

    X = vector_df.loc[:, [c for c in vector_df.columns if c not in ['Patient', 'Mutant_seq']]].to_numpy().transpose()
    x_pca = pca.fit_transform(X)
    classifiers = [c.rsplit("_", 1)[0] for c in vector_df.columns if c not in ['Patient', 'Mutant_seq']]
    pca_df = pd.DataFrame({'PCA_1': x_pca[:, 0], 'PCA_2': x_pca[:, 1], 'Classifier': classifiers})
    variance = pca.explained_variance_ratio_

    g = sns.scatterplot(data=pca_df, x='PCA_1', y='PCA_2', hue='Classifier', alpha=0.8, s=100, palette=color_map)
    plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100), size=args.label_size)
    plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100), size=args.label_size)
    plt.xticks(fontsize=args.tick_size)
    plt.yticks(fontsize=args.tick_size)
    plt.legend(loc="best", frameon=True, fontsize=args.tick_size)
    [x.set_linewidth(args.frame_width) for x in g.axes.spines.values()]
    g.figure.tight_layout()
    plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution, transparent=True)
    plt.close()

