import argparse
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from os.path import exists

from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-p', '--patients', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-ot', '--output_file_tag', type=str, default='', help='File tag for output files')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features to test (numerical or categorical)')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='n',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-lep', '--legend_position', type=str, default='best', help='Legend position in plot')
parser.add_argument('-an', '--add_numbers', dest='add_numbers', action='store_true',
                    help='Add counts to categorical plots')
parser.add_argument('-clf', '--classifier', type=str, default='', help='classifier file for short peptides')
parser.add_argument('-cf', '--classifier_features', type=str, nargs='+', help='features used to train classifier')
parser.add_argument('-cn', '--classifier_features_normalizer', type=str, default='n', help='normalizer for classifier')
parser.add_argument('-ct', '--classifier_tag', type=str, default='LR', help='Tag for classifier')


args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

color_CD8 = 'darkorange'
color_negative = 'royalblue'
axis_label_size = 15
legend_size = 12
tickmark_size = 12
nr_bins=15
figure_size = (12, 8)

classifier = None
classifier_features = None
if len(args.classifier) > 0 and exists(args.classifier):
    with open(args.classifier, mode='r') as clf_file:
        cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
        cat_idx = [i for i, f in enumerate(args.features) if f in Parameters().get_categorical_features()]

        classifier_tag = os.path.basename(args.classifier).split('_')[0]
        classifier = \
            PrioritizationLearner.load_classifier(classifier_tag, OptimizationParams(), args.classifier)

normalizer = get_normalizer(args.normalizer)

if not args.features or len(args.features) == 0:
    features = Parameters().get_ml_features()
else:
    features = args.features

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0, max_netmhc_rank=10000,
                         classifier=classifier, classifier_features=args.classifier_features,
                         classifier_normalizer=args.classifier_features_normalizer, classifier_tag=args.classifier_tag)

# perform leave one out on training set
patients = get_valid_patients(args.patients)

data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)
p_values = {}

patient_str = args.patients[0] if len(args.patients) == 1 else '_'.join(args.patients)
out_tag = 'PV' if len(args.output_file_tag) == 0 else 'PV_'+args.output_file_tag
with open(DataManager().get_result_file(out_tag, patient_str, args.peptide_type, 'txt', 'plt'), mode='w') \
        as result_file:

    out_tag = 'FeatureHistograms' if len(args.output_file_tag) == 0 else 'FeatureHistograms_'+args.output_file_tag
    pdf_file = DataManager().get_result_file(out_tag, patient_str, args.peptide_type, 'pdf', 'plt')
    with PdfPages(pdf_file) as pp:

        firstPage = plt.figure(figsize=figure_size)
        firstPage.clf()
        txt = 'This is the title page'
        firstPage.text(0.05, 0.8, ",".join(patients), transform=firstPage.transFigure, size=20, ha="left")
        pp.savefig()
        plt.close()

        for f in features:
            if f not in X_train:
                continue

            print("Feature {}".format(f))
            if type(normalizer) == dict:
                normalizer_name = get_normalizer_name(normalizer[f])
            else:
                normalizer_name = get_normalizer_name(normalizer)

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
                fig, ax = plt.subplots(figsize=figure_size)
                g = sns.histplot(
                    data=df, x=f, hue="response", palette={0: color_negative, 1: color_CD8}, fill=True,
                    common_norm=False, alpha=.7, linewidth=0, stat="density", legend=False, bins=nr_bins,
                )

                plt.xlabel(f, size=axis_label_size)
                plt.ylabel("Density", size=axis_label_size)
                plt.xticks(fontsize=tickmark_size)
                plt.yticks(fontsize=tickmark_size)
                g.set_title("Normalizer: {0}, t-test p-value = {1:.5e}".format(normalizer_name, tt.pvalue))
                plt.legend(title='Response type', loc=args.legend_position, labels=['CD8+', 'negative'],
                           fontsize=legend_size, title_fontsize=legend_size)
                g.figure.tight_layout()
                pp.savefig(g.figure)
                g.figure.clf()

            if f in Parameters().get_ordinal_features():
                counts = Counter(v)
                counts1 = Counter(v[y_train == 1])
                counts0 = Counter(v[y_train == 0])
                x = np.array([])
                y = np.array([])
                for key in v.unique():
                    if key in counts1:
                        v1 = counts1[key]
                    else:
                        v1 = 0
                    if key in counts0:
                        v0 = counts0[key]
                    else:
                        v0 = 0

                    x = np.append(x, v1)
                    y = np.append(y, v0)

                idx = (x > 0) | (y > 0)
                cont_table = np.array([x[idx], y[idx]])
                chi2, p, dof, ex = stats.chi2_contingency(cont_table)
                p_values[f] = p

                x = np.divide(x, sum(counts1.values()))
                y = np.divide(y, sum(counts0.values()))

                fig, ax = plt.subplots(figsize=figure_size)
                lbls = list(v.unique())
                ax.bar(x=lbls, height=x, color=color_CD8, label='CD8+', alpha=0.7)
                ax.bar(x=lbls, height=y, color=color_negative, label='negative', alpha=0.7)
                plt.xlabel(f, size=axis_label_size)
                plt.ylabel("Density", size=axis_label_size)
                plt.xticks(fontsize=tickmark_size)
                plt.yticks(fontsize=tickmark_size)
                plt.title("Normalizer: {0}, chi2-test p-value = {1:.5e}".format(normalizer_name, p))
                plt.legend(title='Response type', loc=args.legend_position, labels=['CD8+', 'negative'],
                           fontsize=15, title_fontsize=15)
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
                v_u = v.unique()
                for key in v_u:
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

                x_n = np.divide(x, sum(counts1.values()))
                y_n = np.divide(y, sum(counts0.values()))

                fig, ax = plt.subplots(figsize=figure_size)
                lbls_pos = np.arange(len(v_u))
                bp_cd8 = ax.bar(x=lbls_pos, height=x_n, color=color_CD8, label='CD8+ kmers', alpha=0.7)
                bp_neg = ax.bar(x=lbls_pos, height=y_n, color=color_negative, label='negative kmers', alpha=0.7)
                ax.set_xticks(lbls_pos)
                ax.set_xticklabels(v_u)
                plt.xlabel(f, size=axis_label_size)
                plt.ylabel("Density", size=axis_label_size)
                plt.xticks(fontsize=tickmark_size)
                plt.yticks(fontsize=tickmark_size)
                plt.title("Normalizer: {0}, chi2-test p-value = {1:.5e}".format(normalizer_name, p), )
                plt.legend(title='Response type', loc=args.legend_position, labels=['CD8+', 'negative'],
                           fontsize=15, title_fontsize=15)
                if args.add_numbers:
                    max_v = np.max(np.c_[x_n, y_n], axis=1)
                    labels = ['{0:.0f} CD8+, {1:.0f} neg'.format(v_x, v_y) for v_x, v_y in zip(x, y)]
                    for i in range(len(lbls_pos)):
                        plt.annotate(labels[i], xy=(lbls_pos[i], max_v[i]), xytext=(0, 5), textcoords="offset points",
                                     ha="center")
                plt.margins(y=0.2)
                fig.tight_layout()
                pp.savefig(fig)
                fig.clf()

    i = 0
    p_values = dict(map(lambda kv: (kv[0], 1 if np.isnan(kv[1]) else kv[1]), p_values.items()))
    p_values = dict(sorted(p_values.items(), key=lambda item: item[1]))
    rejected, pv_corr, _, alpha_corr = \
        multipletests(list(p_values.values()), alpha=0.1, method='bonferroni', is_sorted=True)

    for arg in vars(args):
        result_file.write("{0}={1}\n".format(arg, getattr(args, arg)))

    result_file.write("Feature\tpValue\trejected\n")
    for f in p_values.keys():
        result_file.write("{0}\t{1:.5e}\t{2:b}\n".format(f, p_values[f], rejected[i]))
        print("{0}\t{1:.5e}\t{2:b}".format(f, p_values[f], rejected[i]))
        i += 1

    pass_f = []
    for f in p_values.keys():
        if p_values[f] < 0.01:
            pass_f.append(f)

    result_file.write("---------------------------------------------------------------------------------------------")
    result_file.write("{0} of total {1} features pass bonferroni correction at alpha={2}. Corrected alpha={3:.5f}\n".
                      format(sum(rejected), len(p_values), 0.1, alpha_corr))
    result_file.write("{0} of a total of {1} features pass p-value threshold of 0.1\n".format(len(pass_f), len(p_values)))
    result_file.write("Features with p-value < 0.1: {0}\n".format(" ".join(pass_f)))

    print("{0} of total {1} features pass bonferroni correction at alpha={2}. Corrected alpha={3:.5f}".
          format(sum(rejected), len(p_values), 0.05, alpha_corr))
    print("{0} of a total of {1} features pass p-value threshold of 0.01".format(len(pass_f), len(p_values)))
    print("Features with p-value < 0.05", " ".join(pass_f))

    out_tag = 'FeaturePValues' if len(args.output_file_tag) == 0 else 'FeaturePValues_'+args.output_file_tag
    pdf_file = DataManager().get_result_file(out_tag, patient_str, args.peptide_type, 'pdf', 'plt')
    with PdfPages(pdf_file) as pp:
        p_values = dict(map(lambda kv: (kv[0], -np.log10(kv[1])), p_values.items()))
#        p_values = dict(sorted(p_values.items(), key=lambda item: item[1], reverse=True))
        fig, ax = plt.subplots(figsize=(30, 8))
        x = np.arange(len(p_values))
        ax.bar(x, p_values.values())
        ax.set_ylabel("-log10(p-value)", fontsize=axis_label_size)
        ax.set_xticks(x)
        ax.set_xticklabels(p_values.keys())
        ax.xaxis.set_tick_params(labelsize=axis_label_size, labelrotation=70)
        # plt.title("All feature p-values sorted", fontsize=12)
        fig.tight_layout()
        pp.savefig(fig)
