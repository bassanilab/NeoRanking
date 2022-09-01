import argparse

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot dataset statistics')

parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

ds_stats_file = os.path.join(Parameters().get_plot_dir(), "Patient_statistics_"+args.peptide_type+".txt")
ds_stats_data = pd.read_csv(ds_stats_file, sep='\t', header=0)

ds_stats_data['Patient group'].replace('Rosenberg', 'NCI', inplace=True)
stat_names = ds_stats_data.columns[2:]

if args.peptide_type == 'long':
    groups = {'Mut count': ['Mut count', 'Mut_tested count', 'Mut_imm count'],
              'RNA expression (TPM)': ['Mean mut RNAseq TPM', 'Mean mut_tested RNAseq TPM',
                                       'Mean mut_imm RNAseq TPM'],
              'RNA coverage': ['Mean mut RNAseq coverage', 'Mean mut_tested RNAseq coverage',
                               'Mean mut_imm RNAseq coverage']
              }
else:
    groups = {'Neo_pep count': ['Neo_pep count', 'Neo_pep_tested count', 'Neo_pep_imm count'],
              'RNA expression (TPM)': ['Mean mut RNAseq TPM', 'Mean mut_tested RNAseq TPM',
                                       'Mean mut_imm RNAseq TPM'],
              'RNA coverage': ['Mean mut RNAseq coverage', 'Mean mut_tested RNAseq coverage',
                               'Mean mut_imm RNAseq coverage'],
              'Allele count': ['Neo_pep allele count', 'Neo_pep_tested allele count', 'Neo_pep_imm allele count'],
              'MixMHC rank': ['Mean neo_pep MixMHC rank', 'Mean neo_pep_tested MixMHC rank',
                              'Mean neo_pep_imm MixMHC rank'],
              'Neo_pep_imm count per mut':
                  ['Mean neo_pep_imm count per mut', 'Mean neo_pep_imm count per mut_tested'],
              }

if args.peptide_type == 'long':
    trafo = {'Mut count': 'log',
             'RNA expression (TPM)': 'log',
             'RNA coverage': 'lin'
             }
    angle = 0
else:
    trafo = {'Neo_pep count': 'log',
             'RNA expression (TPM)': 'log',
             'RNA coverage': 'lin',
             'Allele count': 'lin',
             'MixMHC rank': 'log',
             'Neo_pep_imm count per mut': 'log'
             }
    angle = 15

for group in groups:
    if group == 'Neo_pep_imm count per mut':
        continue

    df0 = ds_stats_data[['Patient group', 'Patient', groups[group][0]]]
    df0 = df0.rename(columns={groups[group][0]: group})
    if args.peptide_type == 'long':
        df0['Response type'] = 'Mut_all'
    else:
        df0['Response type'] = 'Neo_pep_all'
    df1 = ds_stats_data[['Patient group', 'Patient', groups[group][1]]]
    df1 = df1.rename(columns={groups[group][1]: group})
    if args.peptide_type == 'long':
        df1['Response type'] = 'Mut_tested'
    else:
        df1['Response type'] = 'Neo_pep_tested'

    df = pd.concat([df0, df1], ignore_index=True)
    df2 = ds_stats_data[['Patient group', 'Patient', groups[group][2]]]
    df2 = df2.rename(columns={groups[group][2]: group})
    if args.peptide_type == 'long':
        df2['Response type'] = 'Mut_imm'
    else:
        df2['Response type'] = 'Neo_pep_imm'
    df = pd.concat([df, df2], ignore_index=True)

    fig = plt.figure(figsize=(10, 6))
    g = sns.boxplot(hue="Patient group", y=group, x='Response type', data=df)
    g.set_xticklabels(g.get_xticklabels(), rotation=angle)
    plt.xlabel("", size=20)
    g.set_ylabel(group, fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=20)
    if trafo[group] == 'log':
        g.set_yscale('log')
    if group == 'Neo-peptide count per mutation':
        plt.legend(loc='upper left', frameon=True, fontsize=15)
    else:
        plt.legend(frameon=True, fontsize=15)
    g.figure.tight_layout()
    tag = group.replace(" ", "_")
    png_file = \
        os.path.join(Parameters().get_plot_dir(), "{0}_{1}_{2}.png".format(args.png_prefix, args.peptide_type, tag))
    plt.savefig(png_file, bbox_inches='tight')
    plt.close()

if args.peptide_type == 'short':
    group = 'Neo_pep_imm count per mut'
    df0 = ds_stats_data[['Patient group', 'Patient', groups[group][0]]]
    df0 = df0.rename(columns={groups[group][0]: group})
    df0['Response type'] = 'Mut_all'
    df1 = ds_stats_data[['Patient group', 'Patient', groups[group][1]]]
    df1 = df1.rename(columns={groups[group][1]: group})
    df1['Response type'] = 'Mut_tested'

    df = pd.concat([df0, df1], ignore_index=True)

    fig = plt.figure(figsize=(10, 6))
    g = sns.boxplot(hue="Patient group", y=group, x='Response type', data=df)
    g.set_xticklabels(g.get_xticklabels(), rotation=0)
    plt.xlabel("", size=20)
    g.set_ylabel(group, fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=20)
    if trafo[group] == 'log':
        g.set_yscale('log')
    plt.legend(loc='upper left', frameon=True, fontsize=15)
    g.figure.tight_layout()
    tag = group.replace(" ", "_")
    png_file = \
        os.path.join(Parameters().get_plot_dir(), "{0}_{1}_{2}.png".format(args.png_prefix, args.peptide_type, tag))
    plt.savefig(png_file, bbox_inches='tight')
    plt.close()

