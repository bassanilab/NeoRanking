import argparse

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import patches

from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot dataset statistics')

parser.add_argument('-d', '--data_dir', type=str, default=GlobalParameters().get_plot_dir(),
                    help='Directory containing Patient_statistics_long/short.txt files')
parser.add_argument('-fp', '--file_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", help='File type for plot (png, svg or pdf')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-rot', '--rotation', type=float, default=30.0, help='x-axis label rotation')
parser.add_argument('-las', '--label_size', type=str, default='x-large',
                    help='Axis label size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-tis', '--tick_size', type=str, default='large',
                    help='Tick size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-les', '--legend_size', type=str, default='large',
                    help='Legend size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-dpi', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=6.0, help='Figure height in inches')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

ds_stats_file = os.path.join(args.data_dir, "Patient_statistics_"+args.peptide_type+".txt")
ds_stats_data = pd.read_csv(ds_stats_file, sep='\t', header=0)

ds_stats_data['Patient group'].replace('Rosenberg', 'NCI', inplace=True)
stat_names = ds_stats_data.columns[2:]

if args.peptide_type == 'long':
    groups = {'Mut-seq count': ['Mut-seq count', 'Mut-seq_tested count', 'Mut-seq_imm count'],
              'RNAseq TPM': ['Mean mut-seq RNAseq TPM', 'Mean mut-seq_tested RNAseq TPM',
                             'Mean mut-seq_imm RNAseq TPM'],
              'RNAseq coverage': ['Mean mut-seq RNAseq coverage', 'Mean mut-seq_tested RNAseq coverage',
                                  'Mean mut-seq_imm RNAseq coverage']
              }
else:
    groups = {'Neo-pep count': ['Neo-pep count', 'Neo-pep_tested count', 'Neo-pep_imm count'],
              'RNAseq TPM': ['Mean mut-seq RNAseq TPM', 'Mean mut-seq_tested RNAseq TPM',
                                       'Mean mut-seq_imm RNAseq TPM'],
              'RNAseq coverage': ['Mean mut-seq RNAseq coverage', 'Mean mut-seq_tested RNAseq coverage',
                               'Mean mut-seq_imm RNAseq coverage'],
              'Allele count': ['Neo-pep allele count', 'Neo-pep_tested allele count', 'Neo-pep_imm allele count'],
              'MixMHCpred rank': ['Mean neo-pep MixMHC rank', 'Mean neo-pep_tested MixMHC rank',
                              'Mean neo-pep_imm MixMHC rank'],
              'Neo-pep_imm count per mut-seq':
                  ['Mean neo-pep_imm count per mut-seq', 'Mean neo-pep_imm count per mut-seq_tested',
                   'Mean neo-pep_imm count per mut-seq_imm']
              }

if args.peptide_type == 'long':
    trafo = {'Mut-seq count': 'log',
             'RNAseq TPM': 'log',
             'RNAseq coverage': 'lin'
             }
    angle = args.rotation
else:
    trafo = {'Neo-pep count': 'log',
             'RNAseq TPM': 'log',
             'RNAseq coverage': 'lin',
             'Allele count': 'lin',
             'MixMHCpred rank': 'log',
             'Neo-pep_imm count per mut-seq': 'log'
             }
    angle = args.rotation

for group in groups:
    if group == 'Neo-pep_imm count per mut-seq':
        continue

    df0 = ds_stats_data[['Patient group', 'Patient', groups[group][0]]]
    df0 = df0.rename(columns={groups[group][0]: group})
    if args.peptide_type == 'long':
        df0['Response type'] = 'Mut-seq'
    else:
        df0['Response type'] = 'Neo-pep'
    df1 = ds_stats_data[['Patient group', 'Patient', groups[group][1]]]
    df1 = df1.rename(columns={groups[group][1]: group})
    if args.peptide_type == 'long':
        df1['Response type'] = 'Mut-seq_tested'
    else:
        df1['Response type'] = 'Neo-pep_tested'

    df = pd.concat([df0, df1], ignore_index=True)
    df2 = ds_stats_data[['Patient group', 'Patient', groups[group][2]]]
    df2 = df2.rename(columns={groups[group][2]: group})
    if args.peptide_type == 'long':
        df2['Response type'] = 'Mut-seq_imm'
    else:
        df2['Response type'] = 'Neo-pep_imm'
    df = pd.concat([df, df2], ignore_index=True)

    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    g = sns.boxplot(hue="Patient group", y=group, x='Response type', data=df, notch=False,
                    palette=dict(NCI="blue", TESLA="Orange", HiTIDE="Green"), hue_order=['NCI', 'TESLA', 'HiTIDE'])
    xlabels = []
    for lbl in g.get_xticklabels():
        if args.peptide_type == 'long':
            xlabels.append(lbl.get_text().replace("Mut-seq_", "Mut-seq\n_"))
        else:
            xlabels.append(lbl.get_text().replace("Neo-pep_", "Neo-pep\n_"))
    g.set_xticklabels(xlabels, rotation=args.rotation)
    plt.xlabel("")
    g.set_ylabel(group, fontsize=args.label_size)
    plt.xticks(fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
    if trafo[group] == 'log':
        g.set_yscale('log')
    ylim = g.get_ylim()
    plt.ylim(ylim[0], ylim[1]+3)
    handles = [patches.Patch(color="blue", label='NCI'),
               patches.Patch(color='Orange', label='TESLA'),
               patches.Patch(color='green', label='HiTIDE')]
    sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.20), ncol=3, handles=handles, title="",
                    frameon=False, fontsize=args.legend_size)
    g.figure.tight_layout()
    tag = group.replace(" ", "_")
    plot_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}_{2}.{3}".
                             format(args.file_prefix, args.peptide_type, tag, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()

if args.peptide_type == 'short':
    group = 'Neo-pep_imm count per mut-seq'
    df0 = ds_stats_data[['Patient group', 'Patient', groups[group][0]]]
    df0 = df0.rename(columns={groups[group][0]: group})
    df0['Response type'] = 'Mut-seq'
    df1 = ds_stats_data[['Patient group', 'Patient', groups[group][1]]]
    df1 = df1.rename(columns={groups[group][1]: group})
    df1['Response type'] = 'Mut-seq_tested'
    df2 = ds_stats_data[['Patient group', 'Patient', groups[group][2]]]
    df2 = df2.rename(columns={groups[group][2]: group})
    df2['Response type'] = 'Mut-seq_imm'

    df = pd.concat([df0, df1, df2], ignore_index=True)

    fig = plt.figure()
    fig.set_figheight(args.figure_height)
    fig.set_figwidth(args.figure_width)
    g = sns.boxplot(hue="Patient group", y=group, x='Response type', data=df, notch=False,
                    palette=dict(NCI="blue", TESLA="Orange", HiTIDE="Green"), hue_order=['NCI', 'TESLA', 'HiTIDE'])
    xlabels = []
    for lbl in g.get_xticklabels():
        if "Neo-pep" in lbl.get_text():
            xlabels.append(lbl.get_text().replace("Neo-pep_", "Neo-pep\n_"))
        elif "Mut-seq" in lbl.get_text():
            xlabels.append(lbl.get_text().replace("Mut-seq_", "Mut-seq\n_"))
    g.set_xticklabels(xlabels, rotation=args.rotation, fontsize=args.label_size)
    plt.yticks(fontsize=args.tick_size)
    plt.xlabel("")
    g.set_ylabel('Neo-pep_imm count\n per mut-seq', fontsize=args.label_size)
    if trafo[group] == 'log':
        g.set_yscale('log')
    ylim = g.get_ylim()
    plt.ylim(ylim[0], ylim[1]+3)
    handles = [patches.Patch(color="blue", label='NCI'),
               patches.Patch(color='Orange', label='TESLA'),
               patches.Patch(color='green', label='HiTIDE')]
    sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.20), ncol=3, handles=handles, title="",
                    frameon=False, fontsize=args.legend_size)
    g.figure.tight_layout()
    tag = group.replace(" ", "_")
    plot_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}_{2}.{3}".
                             format(args.file_prefix, args.peptide_type, tag, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
    plt.close()

