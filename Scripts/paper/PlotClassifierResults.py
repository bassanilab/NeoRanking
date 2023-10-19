import argparse
import glob
import ast

import numpy as np
import pandas as pd
import matplotlib
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import patches
import seaborn as sns
import re

from Utils.GlobalParameters import GlobalParameters
from Utils.Util_fct import *
from Classifier.ClassifierResultParser import ClassifierResults, SimpleRankingResults, GartnerResultsNeopep, \
    GartnerResultsMutation
from Utils.DataManager import DataManager


parser = argparse.ArgumentParser(description='Plot and test difference between classifier ranking')
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
parser.add_argument('-plt', '--plot_type', type=str, choices=['topn_counts', 'rank_score'],
                    help='Plot topN counts or rank_score')
parser.add_argument('-g', '--dataset_plot_names', type=str, default='', help='Patient groups plot names')
parser.add_argument('-cln', '--classifier_plot_names', type=str, default='', help='Classifier order in plots')
parser.add_argument('-o', '--plot_order', type=str, help='Order of classifier plot names')
parser.add_argument('-a', '--alpha', type=float, default=0.02, help='Coefficient alpha in score function')
parser.add_argument('-ga', '--gartner', dest='gartner', action='store_true', help='Include Gartner et al. comparison')
parser.add_argument('-sr', '--simple_ranking', dest='simple_ranking', action='store_true', help='Include simple ranking comparison')
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

    output_file = os.path.join(GlobalParameters.plot_dir, "{0}.txt".format(args.file_name))
    with open(output_file, "w") as file:
        for arg in vars(args):
            file.write("#{0}={1}\n".format(arg, getattr(args, arg)))

    data = DataManager.filter_selected_data(peptide_type=args.peptide_type, dataset=args.dataset)

    included_patients = set(data['patient'].unique())
    if args.gartner:
        assert args.dataset == 'NCI_test', "-ga option only compatible with NCI_test dataset"
        # get included patients from gartner_ranking file !!
        # add patients to simple ranking
        if args.peptide_type == 'neopep':
            included_patients = included_patients.intersection(GartnerResultsNeopep.get_gartner_test_patients())
        else:
            included_patients = included_patients.intersection(GartnerResultsMutation.get_gartner_test_patients())

    plot_df = None
    topN_dfs = None
    datasets = set()

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

    max_rank_score = None
    tot_imm_count = None
    if args.plot_order:
        plot_order = args.plot_order.split(',')
    else:
        plot_order = None

    for j, regexp in enumerate(args.clf_result_files_re):
        clf_result_files = glob.glob(os.path.join(GlobalParameters.classifier_result_dir, args.sub_dir, regexp))

        for i, clf_result_file in enumerate(clf_result_files):
            if os.path.getsize(clf_result_file) > 0:
                print(clf_result_file)
                if 'SimpleRanking' in clf_result_file:
                    if args.simple_ranking:
                        simple_results = SimpleRankingResults(args.peptide_type, clf_result_file, plt_name_dict[j],
                                                              data, alpha=args.alpha, included_patients_=included_patients)
                        plot_df = simple_results.add_to_sum_df(plot_df)
                        for ds in topN_dfs:
                            topN_dfs[ds] = simple_results.add_to_topN_df(topN_dfs[ds], ds)
                else:
                    clf_results = ClassifierResults(args.peptide_type, clf_result_file, i, plt_name_dict[j],
                                                    data, alpha=args.alpha, included_patients_=included_patients)
                    if max_rank_score is None:
                        max_rank_score = clf_results.get_max_score()
                    if tot_imm_count is None:
                        tot_imm_count = clf_results.get_imm_count()
                    if topN_dfs is None:
                        topN_dfs = {'all': None}
                        for g in clf_results.datasets:
                            topN_dfs[g] = None
                    plot_df = clf_results.add_to_sum_df(plot_df)
                    datasets = clf_results.add_to_datasets(datasets)
                    for ds in topN_dfs:
                        topN_dfs[ds] = clf_results.add_to_topN_df(topN_dfs[ds], ds)

                if args.gartner:
                    if args.peptide_type == 'neopep':
                        gartner_results = GartnerResultsNeopep(clf_results.results, alpha=args.alpha,
                                                               included_patients_=included_patients)
                        plot_df = gartner_results.add_to_sum_df(plot_df)
                        for ds in topN_dfs:
                            topN_dfs[ds] = gartner_results.add_to_topN_df(topN_dfs[ds], ds)
                        parse_gartner = False
                    else:
                        gartner_results = GartnerResultsMutation(clf_results.results, alpha=args.alpha,
                                                                 included_patients_=included_patients)
                        plot_df = gartner_results.add_to_sum_df(plot_df)
                        for ds in topN_dfs:
                            topN_dfs[ds] = gartner_results.add_to_topN_df(topN_dfs[ds], ds)
                        parse_gartner = False

    plot_df = plot_df.astype({'Classifier': str, 'Score': float})
    for group in datasets:
        plot_df.astype({group+"_score": float})

    dataset_dict = ast.literal_eval(args.dataset_plot_names)

    if args.rotation > 0:
        ha = 'center'
    else:
        ha = 'center'

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)

    if args.plot_type == 'rank_score':
        if args.bar_plot:
            g = sns.barplot(x="Classifier", y=args.dataset+"_score", data=plot_df, order=plot_order, color=args.one_color,
                            estimator=np.mean, errorbar=('ci', 95), palette=color_map)
        else:
            g = sns.boxplot(x="Classifier", y=args.dataset+"_score", data=plot_df, order=plot_order, color=args.one_color,
                            palette=color_map)
        sns.swarmplot(x="Classifier", y=args.dataset+"_score", data=plot_df, color=".25", order=plot_order)
        lbs = g.get_xticklabels()
        g.set_xticklabels(lbs, rotation=args.rotation, ha=ha)
        if args.dataset in ylim_dict:
            g.set(ylim=(ylim_dict[args.dataset], None))
        plt.ylabel('rank_score', fontsize=args.label_size, style='italic')
        plt.xlabel(args.xlabel, fontsize=args.label_size)
        plt.xticks(fontsize=args.label_size)
        plt.yticks(fontsize=args.tick_size)
        [x.set_linewidth(args.frame_width) for x in g.axes.spines.values()]
    #    plt.title(patient_dict[args.dataset], fontsize=args.tick_size)
        g.set_title("{0}{1}: maximal rank_score = {2:.3f}".
                    format(args.title_prefix, dataset_dict[args.dataset], max_rank_score[args.dataset]), fontsize=args.title_size)
        plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution, transparent=args.file_type == 'pdf')
        plt.close()

        with open(output_file, "a") as file:
            value_dict = {}

            def add_to_dict(grouped_df):
                value_dict[grouped_df.name] = grouped_df[args.dataset+"_score"]

            file.write("{0} t-tests:\n".format(args.dataset + "_rank_score"))
            plot_df.groupby('Classifier').apply(add_to_dict)

            x_values = [*value_dict]
            for i in range(len(x_values) - 1):
                x_1 = x_values[i]
                for j in range(i+1, len(x_values)):
                    x_2 = x_values[j]
                    pv = stats.ttest_ind(value_dict[x_values[i]], value_dict[x_values[j]])
                    file.write("{0} vrs {1}; t-test statistic={2:e}, p-value={3:e}\n".
                               format(x_1, x_2, pv.statistic, pv.pvalue))

    if args.plot_type == 'topn_counts':
        g = sns.barplot(x='Classifier', y='Neo_pep_imm count', hue='Top N', data=topN_dfs[args.dataset], estimator=np.mean,
                        errorbar=('ci', 95), order=plot_order, palette=topN_palette)
        lbs = g.get_xticklabels()
        g.set_xticklabels(lbs, rotation=args.rotation, fontsize=args.label_size, ha=ha)
        peptide_label = 'neo-pep' if args.peptide_type == 'neopep' else 'mut-seq'
        plt.ylabel("{0}_imm count".format(peptide_label), size=args.label_size)
        plt.xlabel(args.xlabel, fontsize=args.label_size)
        plt.xticks(fontsize=args.label_size)
        plt.yticks(fontsize=args.tick_size)
        plt.ylim(0, tot_imm_count[args.dataset]+3)
        plt.axhline(y=tot_imm_count[args.dataset], color="red", linestyle="--", linewidth=2)
        handles = [patches.Patch(color=topN_palette[0], label='Top 20'),
                   patches.Patch(color=topN_palette[1], label='Top 50'),
                   patches.Patch(color=topN_palette[2], label='Top 100')
                   ]
        sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.20), ncol=3, handles=handles,
                        frameon=False, fontsize=args.legend_size, title=dataset_dict[args.dataset],
                        title_fontsize=args.legend_size)
        [x.set_linewidth(args.frame_width) for x in g.axes.spines.values()]

        plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution, transparent=args.file_type == 'pdf')
        plt.close()

        df_agg = topN_dfs[args.dataset].groupby(['Classifier', 'Top N']).agg({'Neo_pep_imm count': 'mean'})
        pd.DataFrame(df_agg).to_csv(output_file, sep='\t', header=True, index=True, mode="a")
        with open(output_file, 'a') as file:
            file.write("\nTotal count: {0}\n".format(tot_imm_count[args.dataset]))

