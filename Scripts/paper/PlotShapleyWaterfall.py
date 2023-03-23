import argparse
import os

import pandas as pd
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import shap
from sklearn.decomposition import PCA
import seaborn as sns
import time

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import get_normalizer, get_valid_patients

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--classifier_dir', type=str, default=Parameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-c', '--classifier_wc', type=str, default='', help='classifier to use')
parser.add_argument('-fp', '--file_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", help='File type for plot (png, svg or pdf')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--train_dataset', type=str, help='dataset for training set')
parser.add_argument('-te', '--test_dataset', type=str, help='dataset for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20, help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-cat', '--cat_encoder', type=str, default='categorical', help='Convert categories to')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-nn', '--nr_negative', type=int, default=100000, help='Maximal number of short, non immunogenic samples')
parser.add_argument('-rot', '--rotation', type=float, default=30.0, help='x-axis label rotation')
parser.add_argument('-las', '--label_size', type=float, default=25.0, help='Axis label size')
parser.add_argument('-tis', '--tick_size', type=float, default=20.0, help='Axis tick size')
parser.add_argument('-tts', '--title_size', type=float, default=20.0, help='title size')
parser.add_argument('-res', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=6.00, help='Figure height in inches')
parser.add_argument('-les', '--legend_size', type=float, default=15, help='Legend size in float')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')
parser.add_argument('-cm', '--color_map', type=str, default='', help='color map for classifiers')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

normalizer = get_normalizer(args.normalizer)
encodings = read_cat_encodings(args.train_dataset, args.peptide_type)

data_loader_train = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                               mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                               immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                               cat_encoders=encodings, max_netmhc_rank=10000)

data_loader_test = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                              mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                              immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                              cat_encoders=encodings, max_netmhc_rank=10000)

patients_train = \
    get_valid_patients(dataset=args.train_dataset, peptide_type=args.peptide_type)
patients_test = \
    get_valid_patients(dataset=args.test_dataset, peptide_type=args.peptide_type)

data_train, X_train, y_train = data_loader_train.load_patients(patients_train, args.input_file_tag, args.peptide_type,
                                                               nr_non_immuno_rows=args.nr_negative)

cat_features = [f for f in X_train.columns if f in Parameters().get_categorical_features()]
cat_idx = [X_train.columns.get_loc(col) for col in cat_features]


def get_learner(classifier_name, x):
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader_train.get_categorical_idx(x),
                                            cat_dims=data_loader_train.get_categorical_dim(),
                                            input_shape=[len(args.features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


classifier_files = glob.glob(os.path.join(args.classifier_dir, args.classifier_wc))

classifiers = []
classifier_tags = []
explainers = []
for i, clf_file in enumerate(classifier_files):
    classifier_tag = os.path.basename(clf_file).split('_')[0]
    learner = get_learner(classifier_tag, X_train)
    classifier = learner.load_classifier(classifier_tag, learner.get_optimization_params(), clf_file)

    if classifier_tag in ['XGBoost', 'CatBoost']:
        explainer = shap.TreeExplainer(classifier)
    elif classifier_tag in ['LR']:
        explainer = shap.Explainer(classifier, X_train, feature_names=args.features)
#        explainer = shap.Permutation(classifier.predict_proba, X_train, feature_names=args.features)

    classifiers.append(classifier)
    classifier_tags.append(classifier_tag)
    explainers.append(explainer)

for p in patients_test:
    data_test, X_test, y_test = data_loader_test.load_patients(p, args.input_file_tag, args.peptide_type, verbose=False)

    if y_test is None or sum(y_test) == 0:
        continue

    shap_values_repl = []
    ranks = []
    for i, (classifier_, classifier_tag_, explainer_) in enumerate(zip(classifiers, classifier_tags, explainers)):
        y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
            learner.test_classifier(classifier_tag_, classifier_, p, data_test, X_test, y_test, max_rank=args.max_rank)
        shap_values = explainer_(X_test)
        fn = shap_values.feature_names
        shap_values.feature_names = [feature_dict[f] for f in fn]
        shap_values_repl.append(shap_values)
        ranks.append(r)

    plot_df_top20 = None
    for i in range(min(len(y_test), 20)):
        rank = 0
        for j in range(len(shap_values_repl)):
            if plot_df_top20 is None:
                plot_df_top20 = pd.DataFrame({'Shapley value': shap_values_repl[j][i].values,
                                              'Feature': shap_values_repl[j][i].feature_names,
                                              'Peptide': 'Top 20'})
            else:
                plot_df_top20 = plot_df_top20.append(pd.DataFrame({'Shapley value': shap_values_repl[j][i].values,
                                                                   'Feature': shap_values_repl[j][i].feature_names,
                                                                   'Peptide': 'Top 20'}))

    for i, i_t in enumerate(np.where(y_test == 1)[0]):
        plot_df = None
        rank = 0
        peptide = data_test.loc[i_t, 'mutant_seq']
        for j in range(len(shap_values_repl)):
            if plot_df is None:
                plot_df = pd.DataFrame({'Shapley value': shap_values_repl[j][i_t].values,
                                        'Feature': shap_values_repl[j][i_t].feature_names,
                                        'Peptide': peptide})
            else:
                plot_df = plot_df.append(pd.DataFrame({'Shapley value': shap_values_repl[j][i_t].values,
                                                       'Feature': shap_values_repl[j][i_t].feature_names,
                                                       'Peptide': peptide}))
            rank += ranks[j][i]

            # print("peptide = {0}, sum_shap = {1:.3f}, rank = {2}".
            #       format(peptide, shap_values_repl[j][i_t].base_values+sum(shap_values_repl[j][i_t].values), ranks[j][i]))

        rank /= len(shap_values_repl)
        plot_df = plot_df.append(plot_df_top20)

        fig, ax = plt.subplots()
        fig.set_figheight(3*args.figure_height)
        fig.set_figwidth(args.figure_width)

        mean_val = plot_df.groupby('Feature').agg('mean')
        colors = {}
        for f in mean_val.index:
            colors[f] = 'r' if mean_val.loc[f, 'Shapley value'] >= 0 else 'b'

        g = sns.barplot(data=plot_df, x='Shapley value', y='Feature', hue='Peptide', errorbar='sd')
        plt.title("{0}, {1}\n avg. rank: {2:.2f}".format(p, peptide, rank), fontsize=args.title_size)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=args.tick_size)
        plt.ylabel("")
        plt.xlabel('Shapley value', fontsize=args.label_size)
        g.legend(fontsize=args.legend_size)
        sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=1,
                        title="", frameon=False, fontsize=args.legend_size, markerscale=10.0)
        plt.tight_layout()
        plot_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_{2}_{3}.{4}".
                                 format(args.file_prefix, p, peptide, args.peptide_type, args.file_type))
        plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()
