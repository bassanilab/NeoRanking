import argparse
import ast
import os

import pandas as pd
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import shap
from sklearn.decomposition import PCA
import seaborn as sns
import time

from DataWrangling.DataTransformer import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import get_normalizer, get_valid_patients

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d1', '--classifier_dir1', type=str, default=GlobalParameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-d2', '--classifier_dir2', type=str, default=GlobalParameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-c1', '--classifier1_wc', type=str, default='', help='classifier to use')
parser.add_argument('-c2', '--classifier2_wc', type=str, default='', help='classifier to use')
parser.add_argument('-fp', '--file_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", help='File type for plot (png, svg or pdf')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-ds1', '--dataset1', type=str, help='patient ids for training set')
parser.add_argument('-ds2', '--dataset2', type=str, help='patient ids for training set')
parser.add_argument('-i1', '--input_file_tag1', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-i2', '--input_file_tag2', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-ct1', '--clf_tag1', type=str, help='Classifier tag for plot')
parser.add_argument('-ct2', '--clf_tag2', type=str, help='Classifier tag for plot')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20, help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
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
encodings = read_cat_encodings(args.dataset1, args.peptide_type)


def load_data(dataset, input_file_tag, peptide_type):
    if dataset == 'NCI_train':
        response_types = ['CD8', 'CD4/CD8', 'negative']
    else:
        response_types = ['CD8', 'CD4/CD8', 'negative', 'not_tested']

    data_loader = DataTransformer(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                                  mutation_types=args.mutation_types, response_types=response_types,
                                  immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                                  cat_encoders=encodings, max_netmhc_rank=10000)

    patients = \
        get_valid_patients(dataset=dataset, peptide_type=peptide_type)

    data, X, y = data_loader.load_patients(patients, input_file_tag, peptide_type, nr_non_immuno_rows=args.nr_negative)
    return data, X, y


def get_learner(classifier_name, x):
    data_loader = DataTransformer(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                                  mutation_types=args.mutation_types, immunogenic=args.immunogenic, min_nr_immuno=0,
                                  cat_type=args.cat_encoder, cat_encoders=encodings, max_netmhc_rank=10000)
    optimizationParams = OptimizationParams(args.alpha, cat_idx=data_loader.get_categorical_idx(x),
                                            cat_dims=data_loader.get_categorical_dim(),
                                            input_shape=[len(args.features)])

    return PrioritizationLearner(classifier_name, args.scorer, optimizationParams, verbose=args.verbose)


data1, X1, y1 = load_data(args.dataset1, args.input_file_tag1, args.peptide_type)

if args.dataset1 == args.dataset2 and args.input_file_tag1 == args.input_file_tag2:
    data2, X2, y2 = data1, X1, y1
else:
    data2, X2, y2 = load_data(args.dataset2, args.input_file_tag2, args.peptide_type)


cat_features = [f for f in X1.columns if f in GlobalParameters().get_categorical_features()]
cat_idx = [X1.columns.get_loc(col) for col in cat_features]


classifier_files1 = glob.glob(os.path.join(args.classifier_dir1, args.classifier1_wc))
datasets = np.full(len(classifier_files1), args.dataset1)
clf_tags = np.full(len(classifier_files1), args.clf_tag1)
file_tags = np.full(len(classifier_files1), args.input_file_tag1)
classifier_files2 = glob.glob(os.path.join(args.classifier_dir2, args.classifier2_wc))
datasets = np.append(datasets, np.full(len(classifier_files2), args.dataset2))
clf_tags = np.append(clf_tags, np.full(len(classifier_files2), args.clf_tag2))
file_tags = np.append(file_tags, np.full(len(classifier_files2), args.input_file_tag2))
classifier_files = classifier_files1 + classifier_files2

plot_df = None
tot_importance = {}
for i, (clf_file, ds, file_tag, clf_tag) in enumerate(zip(classifier_files, datasets, file_tags, clf_tags)):

    if ds == args.dataset1 and file_tag == args.input_file_tag1:
        data, X, y = data1, X1, y1
    else:
        data, X, y = data2, X2, y2

    classifier_tag = os.path.basename(clf_file).split('_')[0]
    learner = get_learner(classifier_tag, X)
    classifier = learner.load_classifier(classifier_tag, learner.get_optimization_params(), clf_file)

    if classifier_tag in ['XGBoost', 'CatBoost']:
        explainer = shap.TreeExplainer(classifier)
    elif classifier_tag in ['LR']:
        explainer = shap.Explainer(classifier, X, feature_names=X.columns)
#        classifier = learner.fit_classifier(X, y, classifier)

    shap_values = explainer(X)

    shap_values.feature_names = [feature_dict[f] for f in shap_values.feature_names]
    feature_importance = shap_values.abs.mean(0).values
    feature_importance /= sum(feature_importance)

    for fn, fi in zip(shap_values.feature_names, feature_importance):
        if fn not in tot_importance:
            tot_importance[fn] = 0
        tot_importance[fn] += fi

    classifier_tag = classifier_tag +"_" + clf_tag

    df = pd.DataFrame({'Feature': shap_values.feature_names, 'Feature importance': feature_importance,
                       'Classifier': classifier_tag, 'Replicate': clf_file})
    if plot_df is None:
        plot_df = df
    else:
        plot_df = plot_df.append(df)

fig, ax = plt.subplots()
fig.set_figheight(3*args.figure_height)
fig.set_figwidth(args.figure_width)
txt_file = os.path.join(GlobalParameters().get_plot_dir(),
                        "{0}_{1}_FeatureImp.txt".format(args.file_prefix, args.peptide_type))
plot_df.to_csv(txt_file, sep='\t', index=False, header=True)

color_map = None
if args.color_map != '':
    cm = ast.literal_eval(args.color_map)
    color_map = {}
    for key in cm:
        if cm[key].isnumeric():
            color_map[key] = sns.color_palette()[cm[key]]
        else:
            color_map[key] = cm[key]

sorted_fn = [k for (k, v) in sorted(tot_importance.items(), key=lambda item: item[1], reverse=True)]
g = sns.barplot(data=plot_df, x='Feature importance', y='Feature', hue='Classifier', errorbar='sd',
                palette=color_map, order=sorted_fn)
plt.title("Datasets: {0}, {1}".format(args.dataset1, args.dataset2), fontsize=args.title_size)
plt.xticks(fontsize=18)
plt.yticks(fontsize=args.tick_size)
plt.ylabel("")
plt.xlabel('mean(|Shapley value|)', fontsize=args.label_size)
g.legend(fontsize=args.legend_size)
plot_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}.{2}".
                         format(args.file_prefix, args.peptide_type, args.file_type))
plt.savefig(plot_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


