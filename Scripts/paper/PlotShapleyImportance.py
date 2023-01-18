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
parser.add_argument('-c1', '--classifier1_wc', type=str, default='', help='classifier to use')
parser.add_argument('-c2', '--classifier2_wc', type=str, default='', help='classifier to use')
parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-ds', '--dataset', type=str, help='patient ids for training set')
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
parser.add_argument('-fp', '--feature_pairs', type=str, nargs='+', help='Features pair for pair plots')
parser.add_argument('-rot', '--rotation', type=float, default=30.0, help='x-axis label rotation')
parser.add_argument('-las', '--label_size', type=float, default=25.0, help='Axis label size')
parser.add_argument('-tis', '--tick_size', type=float, default=20.0, help='Axis tick size')
parser.add_argument('-tts', '--title_size', type=float, default=20.0, help='title size')
parser.add_argument('-res', '--resolution', type=float, default=600, help='Figure resolution in dots per inch')
parser.add_argument('-fiw', '--figure_width', type=float, default=10.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=6.00, help='Figure height in inches')
parser.add_argument('-les', '--legend_size', type=float, default=15, help='Legend size in float')
parser.add_argument('-hy', '--hyperopt', dest='hyperopt', action='store_true', help='Include hyperopt training score')
parser.add_argument('-ttp', '--title-prefix', type=str, default='', help='title prefix')
parser.add_argument('-mrn', '--max_rank_netmhc', type=int, default=20000, help='Maximal netmhc rank of short peptide')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

normalizer = get_normalizer(args.normalizer)
encodings = read_cat_encodings(args.dataset, args.peptide_type)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                         immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                         cat_encoders=encodings, max_netmhc_rank=10000)

patients = \
    get_valid_patients(dataset=args.dataset, peptide_type=args.peptide_type)

data, X, y = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type,
                                       nr_non_immuno_rows=args.nr_negative)

cat_features = [f for f in X.columns if f in Parameters().get_categorical_features()]
cat_idx = [X.columns.get_loc(col) for col in cat_features]


def get_learner(classifier_name, x):
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader.get_categorical_idx(x),
                                            cat_dims=data_loader.get_categorical_dim(),
                                            input_shape=[len(args.features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


classifier_files = []
classifier_files = classifier_files + glob.glob(os.path.join(args.classifier_dir, args.classifier1_wc))
classifier_files = classifier_files + glob.glob(os.path.join(args.classifier_dir, args.classifier2_wc))

plot_df = None
tot_importance = {}
for i, clf_file in enumerate(classifier_files):

    classifier_tag = os.path.basename(clf_file).split('_')[0]
    learner = get_learner(classifier_tag, X)
    classifier = learner.load_classifier(classifier_tag, learner.get_optimization_params(), clf_file)

    if classifier_tag in ['XGBoost', 'CatBoost']:
        explainer = shap.TreeExplainer(classifier)
    elif classifier_tag in ['LR']:
        explainer = shap.Explainer(classifier, X, feature_names=args.features)
        classifier = learner.fit_classifier(X, y, classifier)

    shap_values = explainer(X)

    shap_values.feature_names = [feature_dict[f] for f in shap_values.feature_names]
    feature_importance = shap_values.abs.mean(0).values
    feature_importance /= sum(feature_importance)

    for fn, fi in zip(shap_values.feature_names, feature_importance):
        if fn not in tot_importance:
            tot_importance[fn] = 0
        tot_importance[fn] += fi

    df = pd.DataFrame({'Feature': shap_values.feature_names, 'Feature importance': feature_importance,
                       'Classifier': classifier_tag, 'Replicate': clf_file})
    if plot_df is None:
        plot_df = df
    else:
        plot_df = plot_df.append(df)

fig, ax = plt.subplots()
fig.set_figheight(3*args.figure_height)
fig.set_figwidth(args.figure_width)
txt_file = os.path.join(Parameters().get_plot_dir(),
                        "{0}_{1}_{2}_FeatureImp.txt".format(args.png_prefix, args.peptide_type, args.dataset))
plot_df.to_csv(txt_file, sep='\t', index=False, header=True)

sorted_fn = [k for (k,v) in sorted(tot_importance.items(), key=lambda item: item[1], reverse=True)]
g = sns.barplot(data=plot_df, x='Feature importance', y='Feature', hue='Classifier', errorbar='sd',
                hue_order=['LR', 'XGBoost'], palette={'LR': 'red', 'XGBoost': 'blue'}, order=sorted_fn)
plt.title("Dataset: {0}".format(args.dataset), fontsize=args.title_size)
plt.xticks(fontsize=args.label_size)
plt.yticks(fontsize=args.tick_size)
plt.ylabel("")
plt.xlabel('mean(|Shapley value|)', fontsize=args.label_size)
g.legend(fontsize=args.legend_size)
png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_{2}.png".
                        format(args.png_prefix, args.dataset, args.peptide_type))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


