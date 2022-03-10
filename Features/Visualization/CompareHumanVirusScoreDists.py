import argparse

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
from Features.BindingInfo.CalcBindingAffinityScoreDists import *


parser = argparse.ArgumentParser(description='Get alleles that present immunogenic mutation')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-res', '--residual_file', type=str, help='residual output file')
parser.add_argument('-prh', '--pr_human', type=str, help='Promiscuity human data input file')
parser.add_argument('-prv', '--pr_virus', type=str, help='Promiscuity viral data input file')
parser.add_argument('-prp', '--promiscuity_paper', type=str, help='Promiscuity data input file from Manczinger et al')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

parameters = Parameters()

calcHumanAllelePropensity = \
    CalcAllelePropensity(parameters.get_human_allele_score_dist_file(), parameters.get_patients_with_immo(),
                         args.mutation_types, args.immunogenic, args.response_types, args.input_file_tag)

calcVirusAllelePropensity = \
    CalcAllelePropensity(parameters.get_virus_allele_score_dist_file(), parameters.get_patients_with_immo(),
                         args.mutation_types, args.immunogenic, args.response_types, args.input_file_tag)

with PdfPages(args.pdf) as pp:
    ratios = calcHumanAllelePropensity.get_immunogenicity_ratios()
    keys = list(ratios.keys())
    vals = [float(ratios[k]) for k in keys]
    fig, ax = plt.subplots(figsize=(20, 8))
    x = np.arange(len(keys))
    ax.bar(x, vals)
    ax.set_xticks(x)
    ax.set_xticklabels(keys)
    ax.xaxis.set_tick_params(labelsize=8, labelrotation=90)
    fig.tight_layout()
    pp.savefig(fig)
    fig.clf()

    virus_prop = []
    human_prop = []
    for a in keys:
        virus_prop = np.append(virus_prop, calcVirusAllelePropensity.get_propensity(a, 9))
        human_prop = np.append(human_prop, calcHumanAllelePropensity.get_propensity(a, 9))

    slope, intercept, r_value, p_value, std_err = stats.linregress(x=human_prop, y=virus_prop)

    res = virus_prop - (intercept + np.multiply(human_prop, slope))

    fig, ax = plt.subplots(figsize=(10, 8))
    g = sns.regplot(x=human_prop, y=virus_prop)

    sorted_vals = sorted(zip(keys, res, np.arange(len(keys))), key=lambda pair: pair[1])

    for v in sorted_vals[:5]:
        plt.text(human_prop[v[2]], virus_prop[v[2]], v[0], horizontalalignment='left', size='medium',
                 color='black', weight='semibold')
    for v in sorted_vals[-5:]:
        plt.text(human_prop[v[2]], virus_prop[v[2]], v[0], horizontalalignment='left', size='medium',
                 color='black', weight='semibold')

    g.set(xlabel='human propensity', ylabel='viral propensity')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()
