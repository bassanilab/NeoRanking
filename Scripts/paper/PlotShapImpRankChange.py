"""
Compare Shapley value ranks for mutations and neo-peptide classifiers
"""

import pandas as pd
import numpy as np
import seaborn as sns
import os

import argparse

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from Utils.GlobalParameters import GlobalParameters

parser = argparse.ArgumentParser(description='Compare Shapley value ranks for mutations and neo-peptide classifiers')
parser.add_argument('-fim', '--feat_imp_file_mut', type=str, help='Feature importance tab file for mutation clf')
parser.add_argument('-fin', '--feat_imp_file_neopep', type=str, help='Feature importance tab file for neopep clf')
parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-las', '--label_size', type=float, default=20.0, help='Axis label size')
parser.add_argument('-tis', '--tick_size', type=float, default=18.0, help='Axis tick size')
parser.add_argument('-tts', '--title_size', type=float, default=20.0, help='Title size')
parser.add_argument('-ps', '--point_size', type=float, default=100.0, help='Point size')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=10.00, help='Figure height in inches')
parser.add_argument('-les', '--legend_size', type=float, default=20, help='Legend size in float')

if __name__ == "__main__":

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    df_neopep = pd.read_csv(args.feat_imp_file_neopep, sep="\t", header=0)
    df_mutation = pd.read_csv(args.feat_imp_file_mut, sep="\t", header=0)

    def group_features(f):
        if 'RNAseq' in f or 'TCGA' in f or 'GTEx' in f:
            return 'RNAseq'
        elif 'ipMSDB' in f:
            return 'ipMSDB'
        elif 'Intogen' in f:
            return 'Intogen'
        elif 'Mut' in f and 'Rank' in f:
            return 'Binding'
        elif ('PRIME' in f or 'MixMHCpred' in f or 'NetMHCpan' in f) and 'Rank' in f:
            return 'Binding'
        else:
            return 'Other'


    df_neopep['Feature group'] = df_neopep.Feature.apply(group_features)
    df_mutation['Feature group'] = df_mutation.Feature.apply(group_features)

    feature_map = {
        'Cancer Cell Fraction': 'Cancer Cell Fraction',
        'Clonality': 'Clonality',
        'RNAseq Expression(TPM)': 'RNAseq Expression(TPM)',
        'RNAseq Mutation Coverage': 'RNAseq Mutation Coverage',
        'CSCAPE Score': 'CSCAPE Score',
        'MixMHCpred Rank': 'Minimal Mut MixMHCpred Rank',
        'PRIME Rank': 'Minimal Mut PRIME Rank',
        'NetMHCpan Rank': 'Best Mut EL Rank',
        'GTEx Mean Tissue Expression': 'GTEx Mean Tissue Expression',
        'TCGA Cancer Expression': 'TCGA Cancer Expression',
        'Gene Driver Intogen': 'Gene Driver Intogen',
        'Intogen Same Mutation Count': 'Intogen Same Mutation Count',
        'Intogen Mutation Driver Statement': 'Intogen Mutation Driver Statement',
        'ipMSDB Peptide Score': 'ipMSDB Peptide Score',
        'ipMSDB Peptide Overlap': 'ipMSDB Peptide Overlap',
        'ipMSDB Mutation Score': 'ipMSDB Mutation Score',
        'ipMSDB Peptide Count': 'ipMSDB Peptide Count',
        'NetStab Rank': 'Best Mut Stab Rank',
        'NetChop CT Score': 'Best Mut NetChop Score',
        'NetMHCpan log_Rank DAI': 'BEST EL Rank DAI'}

    features_neo = np.array([])
    features_mut = np.array([])
    neo_scores = np.array([])
    mut_scores = np.array([])
    groups = np.array([])
    markers = np.array([])
    handles = np.array([])

    palette = sns.color_palette(n_colors=len(df_neopep['Feature group'].unique()))
    color_map = {}
    for i, g in enumerate(df_neopep['Feature group'].unique()):
        color_map[g] = palette[i]

    marker_list = Line2D.filled_markers

    for i, (f1, f2) in enumerate(feature_map.items()):
        print(f1, f2)
        g = df_mutation.loc[df_mutation['Feature'] == f2, 'Feature group'].iloc[0]
        if g == 'Other':
            continue
        features_neo = np.append(features_neo, f1)
        features_mut = np.append(features_mut, f2)
        neo_scores = np.append(neo_scores, df_neopep.loc[df_neopep['Feature'] == f1, 'Feature importance'].mean())
        mut_scores = np.append(mut_scores, df_mutation.loc[df_mutation['Feature'] == f2, 'Feature importance'].mean())
        groups = np.append(groups, g)
        m = marker_list[i % (len(marker_list))]
        markers = np.append(markers, m)
        handles = np.append(handles, Line2D([0], [0], color=color_map[g], label=f1, alpha=0.7, marker=m))

    plot_df = pd.DataFrame({'Feature neopep': features_neo,
                            'Feature mutation': features_mut,
                            'Neopep mean importance': neo_scores,
                            'Mutation mean importance': mut_scores,
                            'Feature group': groups})

    fig, ax = plt.subplots()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)

    g = sns.scatterplot(data=plot_df, x="Neopep mean importance", y="Mutation mean importance",
                        hue='Feature group', s=args.point_size)

    plt.title("Dataset: NCI-train", fontsize=args.title_size)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xticks(fontsize=args.tick_size)
    plt.yticks(fontsize=args.tick_size)
    plt.xlabel('Neo-peptide feature importance', fontsize=args.label_size)
    plt.ylabel('Mutation feature importance', fontsize=args.label_size)
    x = np.array(plot_df['Neopep mean importance'])
    plt.plot(x, x, color='r')

    g.legend(fontsize=args.legend_size)

    plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', transparent=True)
    plt.close()

    txt_file = os.path.join(GlobalParameters.plot_dir, "{0}.txt".format(args.file_name))
    plot_df.to_csv(txt_file, sep='\t', header=True, index=False)
    plt.close()
