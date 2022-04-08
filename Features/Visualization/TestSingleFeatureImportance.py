import argparse
from DataWrangling.DataLoader import *
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from Visualization.PCAClassifyPeptideBrowser import *
from collections import Counter
from statsmodels.stats.multitest import multipletests
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-ov', '--pdf_overview', type=str, help='PDF output file for all feature p-values')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')


args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

normalizer = None
if args.normalizer == 'q':
    normalizer = QuantileTransformer()
    if args.verbose > 0:
        print('Normalizer: QuantileTransformer')
elif args.normalizer == 'z':
    normalizer = StandardScaler()
    if args.verbose > 0:
        print('Normalizer: StandardScaler')
elif args.normalizer == 'p':
    normalizer = PowerTransformer()
    if args.verbose > 0:
        print('Normalizer: PowerTransformer')
else:
    if args.verbose > 0:
        print('Normalizer: None')


if not args.features or len(args.features) == 0:
    features = Parameters().get_ml_features()
else:
    features = args.features

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immono=0, max_netmhc_rank=2)

# perform leave one out on training set
patients = get_valid_patients(args.patients)

data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type, 'Nb_Samples')
# data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)
p_values = {}

patient_str = args.patients[0] if len(args.patients) == 1 else '_'.join(args.patients)
with open(DataManager().get_result_file('PV', patient_str, args.peptide_type), mode='w') \
        as result_file:

    with PdfPages(args.pdf) as pp:

        firstPage = plt.figure(figsize=(11.69, 8.27))
        firstPage.clf()
        txt = 'This is the title page'
        firstPage.text(0.05, 0.8, ",".join(patients), transform=firstPage.transFigure, size=20, ha="left")
        pp.savefig()
        plt.close()

        for f in features:
            if f not in X_train:
                continue

            print("Feature {}".format(f))

            v_norm = X_train[f]
            v = data_train[f]

            if f in Parameters().get_numerical_features():
                v_norm = np.array(v_norm, dtype=float)
                x = v_norm[y_train == 1]
                y = v_norm[y_train == 0]
#                x = x[np.where(~np.isnan(x))]
#                y = y[np.where(~np.isnan(y))]
                tt = stats.ttest_ind(x, y, nan_policy='omit', equal_var=False, alternative='two-sided')
                p_values[f] = tt.pvalue

                df = pd.DataFrame({f: v_norm, 'response': y_train})
                sns.set(font_scale=2.2)
                g = sns.histplot(
                   data=df, x=f, hue="response",
                   fill=True, common_norm=False, palette="seismic",
                   alpha=.7, linewidth=0, stat="density"
                )

                g.set_title("t-test p-value = {0:.5e}".format(tt.pvalue))
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

            if f in Parameters().get_ordinal_features():
                counts = Counter(v)
                counts1 = Counter(v[y_train == 1])
                counts0 = Counter(v[y_train == 0])
                x = []
                y = []
                for key in counts.keys():
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
                p_values[f] = p

                x = np.divide(x, sum(counts1.values()))
                y = np.divide(y, sum(counts0.values()))

                fig, ax = plt.subplots(figsize=(15, 8))
                lbls = list(counts.keys())
                ax.bar(x=lbls, height=x, color='r', label='CD8+ kmers', alpha=0.7)
                ax.bar(x=lbls, height=y, color='b', label='negative kmers', alpha=0.7)
                ax.xaxis.set_tick_params(labelsize=20)
                plt.title("{0}: chi2-test p-value = {1:.5e}".format(f, p))
                plt.legend()
                fig.tight_layout()
                pp.savefig(fig)
                fig.clf()

            if f in Parameters().get_categorical_features():
                counts = Counter(v)
                counts1 = Counter(v[y_train == 1])
                counts0 = Counter(v[y_train == 0])
                x = []
                y = []
                r = []
                for key in counts.keys():
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
                p_values[f] = p

                x = np.divide(x, sum(counts1.values()))
                y = np.divide(y, sum(counts0.values()))

                fig, ax = plt.subplots(figsize=(15, 8))
                lbls_pos = np.arange(len(counts.values()))
                ax.bar(x=lbls_pos, height=x, color='r', label='CD8+ kmers', alpha=0.7)
                ax.bar(x=lbls_pos, height=y, color='b', label='negative kmers', alpha=0.7)
                ax.set_xticks(lbls_pos)
                ax.set_xticklabels(counts.keys())
                ax.xaxis.set_tick_params(labelsize=20, labelrotation=90)
                plt.title("{0}: chi2-test p-value = {1:.5e}".format(f, p))
                plt.legend()
                fig.tight_layout()
                pp.savefig(fig)
                fig.clf()

    i = 0
    p_values = dict(map(lambda kv: (kv[0], 1 if np.isnan(kv[1]) else kv[1]), p_values.items()))
    p_values = dict(sorted(p_values.items(), key=lambda item: item[1]))
    rejected, pv_corr, _, alpha_corr = \
        multipletests(list(p_values.values()), alpha=0.1, method='bonferroni', is_sorted=True)

    for f in p_values.keys():
        result_file.write("{0}\t{1:.3e}\t{2:b}\n".format(f, p_values[f], rejected[i]))
        print("{0}\t{1:.5f}\t{2:b}".format(f, p_values[f], rejected[i]))
        i += 1

    pass_f = []
    for f in p_values.keys():
        if p_values[f] < 0.01:
            pass_f.append(f)

    result_file.write("{0} of total {1} features pass bonferroni correction at alpha={2}. Corrected alpha={3:.5f}\n".
          format(sum(rejected), len(p_values), 0.1, alpha_corr))
    result_file.write("{0} of a total of {1} features pass p-value threshold of 0.1\n".format(len(pass_f), len(p_values)))
    result_file.write("Features with p-value < 0.1: {0}\n".format(" ".join(pass_f)))

    print("{0} of total {1} features pass bonferroni correction at alpha={2}. Corrected alpha={3:.5f}".
          format(sum(rejected), len(p_values), 0.05, alpha_corr))
    print("{0} of a total of {1} features pass p-value threshold of 0.01".format(len(pass_f), len(p_values)))
    print("Features with p-value < 0.05", " ".join(pass_f))

    with PdfPages(args.pdf_overview) as pp:
        p_values = dict(map(lambda kv: (kv[0], -np.log10(kv[1])), p_values.items()))
#        p_values = dict(sorted(p_values.items(), key=lambda item: item[1], reverse=True))
        fig, ax = plt.subplots(figsize=(30, 8))
        x = np.arange(len(p_values))
        ax.bar(x, p_values.values())
        ax.set_ylabel("-log10(p-value)", fontsize=20)
        ax.set_xticks(x)
        ax.set_xticklabels(p_values.keys())
        ax.xaxis.set_tick_params(labelsize=8, labelrotation=90)
        # plt.title("All feature p-values sorted", fontsize=12)
        fig.tight_layout()
        pp.savefig(fig)
