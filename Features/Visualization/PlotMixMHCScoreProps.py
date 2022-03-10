import argparse

import numpy as np
import pandas as pd

from DataWrangling.DataLoader import *
from collections import Counter
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from sklearn.linear_model import LogisticRegression


parser = argparse.ArgumentParser(description='Get alleles that present immunogenic mutation')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-prf', '--prf', type=str, help='Promiscuity data input file')
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

data = pd.read_csv(args.prf, sep="\t", header=0)

with PdfPages(args.pdf) as pp:

    for c in data.columns[3:]:

        r_all = []
        v_all = []
        L_all = []
        pred_all = []
        for peptide_length in [9, 10, 11, 12]:
            data_len = data.loc[data['Length'] == peptide_length, :]
            r = []
            v = []
            L = []
            for i in data_len.index:
                a = 'HLA-'+str(data_len.loc[i, 'Allele'])

                if a in ratios:
                    r = np.append(r, ratios[a])
                    v = np.append(v, data_len.loc[i, c])
                    L = np.append(L, peptide_length)

            r_all = np.append(r_all, r)
            v_all = np.append(v_all, v)
            L_all = np.append(L_all, L)

        # correlation, p_value = stats.pearsonr(v, r)
        # ttl = "Best binders. Correlation={0:.3f}, p-value={1:.3e}".format(correlation, p_value)
        # print("feature="+c+", "+ttl)
        #
        # g = sns.regplot(x=v, y=r)
        # g.set_title(ttl, fontsize=10)
        # g.set(xlabel=c, ylabel='Best binder ratio')
        # g.figure.tight_layout()
        # pp.savefig(g.figure)
        # g.figure.clf()
        #
        df = pd.DataFrame({c: v_all, 'ratio': r_all, 'immunogenic': r_all > 0, 'length': L_all})
        t, p_value = stats.ttest_ind(v_all[r_all > 0], v_all[r_all == 0])
        ttl = "Best binders. t={0:.3f}, p-value={1:.3e}".format(t, p_value)
        print("feature="+c+", "+ttl)
        g = sns.boxplot(x="length", y=c, hue="immunogenic", data=df, palette="Set3")
        g.set_title(ttl, fontsize=8)
        g.set(xlabel='peptide length', ylabel=c)
        g.figure.tight_layout()
        pp.savefig(g.figure)
        g.figure.clf()

    r_all = []
    L_all = []
    pred_all = []
    for peptide_length in [9, 10, 11, 12]:
        data_len = data.loc[data['Length'] == peptide_length, :]
        r = []
        L = []
        X = np.empty((0, len(data.columns)-3), int)
        for i in data_len.index:
            a = 'HLA-'+str(data_len.loc[i, 'Allele'])

            if a in ratios:
                r = np.append(r, ratios[a])
                L = np.append(L, peptide_length)
                X = np.append(X, [np.array(data_len.loc[i, :])[3:]], axis=0)

        reg = LogisticRegression(penalty='l2', C=1, class_weight='balanced').fit(X, [1 if rr > 0 else 0 for rr in r])
        p_pred = reg.predict_proba(X)[:, 1]

        r_all = np.append(r_all, r)
        L_all = np.append(L_all, L)
        pred_all = np.append(pred_all, p_pred)

    df = pd.DataFrame({'Predicted': pred_all, 'Ratios': r_all, 'Immunogenic': r_all > 0, 'Length': L_all})
    t, p_value = stats.ttest_ind(pred_all[r_all > 0], pred_all[r_all == 0])
    ttl = "Weak binders. t={0:.3f}, p-value={1:.3e}".format(t, p_value)
    g = sns.boxplot(x="Length", y='Predicted', hue="Immunogenic", data=df, palette="Set3")
    g.set_title(ttl, fontsize=8)
    g.set(xlabel='Peptide length', ylabel='Predicted')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()

    r_all = []
    L_all = []
    X_all = np.empty((0, len(data.columns)-3), int)
    for peptide_length in [9, 10, 11, 12]:
        data_len = data.loc[data['Length'] == peptide_length, :]
        for i in data_len.index:
            a = 'HLA-'+str(data_len.loc[i, 'Allele'])

            if a in ratios:
                r_all = np.append(r_all, ratios[a])
                L_all = np.append(L_all, peptide_length)
                X_all = np.append(X_all, [np.array(data_len.loc[i, :])[3:]], axis=0)

    reg = LogisticRegression(penalty='l2', C=1, class_weight='balanced').\
        fit(X_all, [1 if rr > 0 else 0 for rr in r_all])
    pred_all = reg.predict_proba(X_all)[:, 1]

    df = pd.DataFrame({'Predicted': pred_all, 'Ratios': r_all, 'Immunogenic': r_all > 0, 'Length': L_all})
    t, p_value = stats.ttest_ind(pred_all[r_all > 0], pred_all[r_all == 0])
    ttl = "Weak binders. t={0:.3f}, p-value={1:.3e}".format(t, p_value)
    g = sns.boxplot(x="Length", y='Predicted', hue="Immunogenic", data=df, palette="Set3")
    g.set_title(ttl, fontsize=8)
    g.set(xlabel='Peptide length', ylabel='Predicted')
    g.figure.tight_layout()
    pp.savefig(g.figure)
    g.figure.clf()





