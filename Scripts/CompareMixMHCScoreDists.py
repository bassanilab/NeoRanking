import argparse
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


mgr = DataManager()
fasta_file = Parameters().get_virus_fasta_file()
proteins = set()
promiscuity = CalcBindingAffinityScoreDists(seq_logo_dir=Parameters().get_data_dir(), seq_fasta_file=fasta_file,
                                            proteins=proteins)


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

data_human = pd.read_csv(args.pr_human, sep="\t", header=0)
data_virus = pd.read_csv(args.pr_virus, sep="\t", header=0)

with PdfPages(args.pdf) as pp:
    for c in data_virus.columns[3:]:

        for peptide_length in [9, 10, 11, 12]:
            data_len_hum = data_human.loc[data_human['Length'] == peptide_length, :]
            data_len_vir = data_virus.loc[data_virus['Length'] == peptide_length, :]

            r = []
            v_h = []
            v_v = []
            L = []
            alleles = []
            for i in data_len_hum.index:
                a = 'HLA-'+str(data_len_hum.loc[i, 'Allele'])

                if a in ratios:
                    r = np.append(r, ratios[a])
                    v_h = np.append(v_h, data_len_hum.loc[i, c])
                    v_v = np.append(v_v, data_len_vir.loc[i, c])
                    L = np.append(L, peptide_length)
                    alleles = np.append(alleles, a)

            correlation, p_value = stats.pearsonr(v_h, v_v)
            ttl = "Feature={0}, Length={1}, best binders. Correlation={2:.3f}, p-value={3:.3e}".\
                format(c, peptide_length, correlation, p_value)
            print(ttl)

            g = sns.regplot(x=v_h, y=v_v)
            plt.plot([min(v_h), max(v_h)], [min(v_h), max(v_h)], 'red', linewidth=2)
            g.set_title(ttl, fontsize=10)
            g.set(xlabel='human', ylabel='virus')
            g.figure.tight_layout()
            pp.savefig(g.figure)
            g.figure.clf()

            slope, intercept, r_value, p_value, std_err = stats.linregress(x=v_h, y=v_v)

            res = v_v - (intercept + np.multiply(v_h, slope))

            df = pd.DataFrame({'Alleles': alleles, 'Residuals': res, 'human': v_h, 'virus': v_v, 'immunogenic': r > 0})
            t, p_value = stats.ttest_ind(res[r > 0], res[r == 0])
            ttl = "best binders. t={0:.3f}, p-value={1:.3e}".format(t, p_value)
            g = sns.boxplot(y='Residuals', x="immunogenic", data=df, palette="Set3")
            g.set_title(ttl, fontsize=8)
            g.set(xlabel='', ylabel='Residuals')
            g.figure.tight_layout()
            pp.savefig(g.figure)
            g.figure.clf()

            sorted_pair = sorted(zip(alleles, res), key=lambda pair: pair[1])

            for p in sorted_pair:
                print("{0}\t{1}\t{2}\t{3}".format(c, peptide_length, p[0], p[1]))


