import argparse

import numpy as np

from DataWrangling.DataLoader import *
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy import stats
from sklearn.linear_model import Lasso


parser = argparse.ArgumentParser(description='Get alleles that present immunogenic mutation')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
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

patients = np.array(args.patients)
features = ['response_type', 'mut_allele_0', 'mut_allele_1', 'mut_allele_2', 'mut_weak_binding_alleles_0']

data_loader = DataLoader(transformer=None, normalizer=None, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0)

data, X, y = data_loader.load_patients(patients, args.input_file_tag)

pos_alleles = []
neg_alleles = []
pos_alleles_all = []
neg_alleles_all = []
for i in data.index:
    if data.loc[i, 'response_type'] in args.immunogenic:
        pos_alleles = np.append(pos_alleles, data.loc[i, 'mut_allele_0'])
        pos_alleles_all = np.append(pos_alleles_all, str(data.loc[i, 'mut_weak_binding_alleles_0']).split(","))
    else:
        neg_alleles = np.append(neg_alleles, data.loc[i, 'mut_allele_0'])
        neg_alleles_all = np.append(neg_alleles_all, str(data.loc[i, 'mut_weak_binding_alleles_0']).split(","))

pos_cnt = Counter(pos_alleles)
neg_cnt = Counter(neg_alleles)
tot_cnt = pos_cnt+neg_cnt
ratios = {}
for (allele, cnt) in tot_cnt.items():
    ratios[allele] = pos_cnt[allele]/cnt

ratios = dict(sorted(ratios.items(), key=lambda item: item[1], reverse=True))
print(ratios)

pos_cnt_all = Counter(pos_alleles_all)
neg_cnt_all = Counter(neg_alleles_all)
tot_cnt_all = pos_cnt_all+neg_cnt_all
ratios_all = {}
for (allele, cnt) in tot_cnt_all.items():
    ratios_all[allele] = pos_cnt_all[allele]/cnt

ratios_all = dict(sorted(ratios_all.items(), key=lambda item: item[1], reverse=True))
print(ratios_all)

promiscuity_paper_data = pd.read_csv(args.promiscuity_paper, sep="\t", header=0)
alleles = np.array(promiscuity_paper_data['Allele'], dtype=str)

with PdfPages(args.pdf) as pp:
    x = []
    r = []
    alleles = np.array(promiscuity_paper_data['Allele'], dtype=str)
    for a in ratios.keys():
        a_short = str(a).replace(":", "").replace("*",  "").replace("HLA-", "")
        if a_short in alleles:
            idx = np.argwhere(np.array(promiscuity_paper_data['Allele']) == a_short)[0][0]
            x = np.append(x, promiscuity_paper_data.loc[promiscuity_paper_data.index[idx], 'Promiscuity'])
            r = np.append(r, ratios[a])

    correlation, p_value = stats.pearsonr(x, r)
    ttl = "Paper data and best binders. Correlation={0:.3f}, p-value={1:.3e}".format(correlation, p_value)
    print(ttl)

    g = sns.regplot(x=x, y=r)
    g.set_title(ttl, fontsize=8)
    g.set(xlabel='Promiscuity', ylabel='Best binders ratio')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

    df = pd.DataFrame({'Promiscuity': x, 'ratio': r, 'immunogenic': r > 0})
    t, p_value = stats.ttest_ind(x[r > 0], x[r == 0])
    ttl = "Paper data and best binders. t={0:.3f}, p-value={1:.3e}".format(t, p_value)
    print(ttl)
    g = sns.boxplot(y=x, x="immunogenic", data=df, palette="Set3")
    g.set_title(ttl, fontsize=8)
    g.set(ylabel='Promiscuity')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

    x = []
    r = []
    alleles = np.array(promiscuity_paper_data['Allele'], dtype=str)
    for a in ratios_all.keys():
        a_short = str(a).replace(":", "").replace("*",  "").replace("HLA-", "")
        if a_short in alleles:
            idx = np.argwhere(np.array(promiscuity_paper_data['Allele']) == a_short)[0][0]
            x = np.append(x, promiscuity_paper_data.loc[promiscuity_paper_data.index[idx], 'Promiscuity'])
            r = np.append(r, ratios_all[a])

    correlation, p_value = stats.pearsonr(x, r)
    ttl = "Paper data and weak binders. Correlation={0:.3f}, p-value={1:.3e}".format(correlation, p_value)
    print(ttl)

    g = sns.regplot(x=x, y=r)
    g.set_title(ttl, fontsize=8)
    g.set(xlabel='Promiscuity', ylabel='Weak binders ratio')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

    df = pd.DataFrame({'Promiscuity': x, 'ratio': r, 'immunogenic': r > 0})
    t, p_value = stats.ttest_ind(x[r > 0], x[r == 0])
    ttl = "Paper data and weak binders. t={0:.3f}, p-value={1:.3e}".format(t, p_value)
    print(ttl)
    g = sns.boxplot(y=x, x="immunogenic", data=df, palette="Set3")
    g.set_title(ttl, fontsize=8)
    g.set(ylabel='Promiscuity')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

