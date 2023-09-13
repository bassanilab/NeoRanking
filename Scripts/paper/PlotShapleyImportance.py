"""
Plot feature importance as mean absolute shapley values
"""

import argparse
import os
import glob
import ast

import pandas as pd
import matplotlib.pyplot as plt
import shap
import seaborn as sns

from DataWrangling.DataTransformer import *
from Classifier.ClassifierManager import *

parser = argparse.ArgumentParser(description='Plot feature importance as mean absolute shapley values')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory holding classifier model files')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-c1', '--classifier1_file_re', type=str, help='classifier files for first batch')
parser.add_argument('-c2', '--classifier2_file_re', type=str, default="", help='classifier files for second batch')
parser.add_argument('-ds1', '--dataset1', type=str, choices=GlobalParameters.datasets,
                    help='Dataset, one of [NCI, NCI_train, NCI_test, TESLA, HiTIDE]')
parser.add_argument('-ds2', '--dataset2', type=str, default="", choices=GlobalParameters.datasets,
                    help='Dataset, one of [NCI, NCI_train, NCI_test, TESLA, HiTIDE]')
parser.add_argument('-ft', '--file_type', type=str, default="pdf", choices=GlobalParameters.plot_file_formats,
                    help='File type for plot (png, svg or pdf)')
parser.add_argument('-fn', '--file_name', type=str, help='Name of plot output file')
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

if __name__ == "__main__":

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))


    def load_data(dataset, peptide_type):
        if dataset == 'NCI_train':
            response_types = ['CD8', 'negative']
        else:
            response_types = ['CD8', 'negative', 'not_tested']

        data, X, y = DataManager.filter_processed_data(peptide_type=peptide_type, objective='ml',
                                                       response_types=response_types,
                                                       dataset=dataset, sample=peptide_type == 'neopep')
        return data, X, y


    def get_clf_mgr(peptide_type: str, dataset_enc: str, classifier_name: str, x: pd.DataFrame, y: list):
        if args.peptide_type == 'neopep':
            class_ratio = sum(y == 1)/sum(y == 0)
        else:
            class_ratio = None

        alpha = GlobalParameters.neopep_alpha if peptide_type == 'neopep' else GlobalParameters.mutation_alpha
        optimizationParams = \
            OptimizationParams(alpha=alpha,
                               cat_idx=DataManager.get_categorical_feature_idx(peptide_type, x),
                               class_ratio=class_ratio)

        return ClassifierManager(classifier_name, 'sum_exp_rank', optimizationParams, verbose=0)


    data1, X1, y1 = load_data(args.dataset1, args.peptide_type)

    if args.dataset1 == args.dataset2:
        data2, X2, y2 = data1, X1, y1
    elif args.dataset2 != "":
        data2, X2, y2 = load_data(args.dataset2, args.peptide_type)
    else:
        data2, X2, y2 = None, None, None

    classifier_files = \
        glob.glob(os.path.join(GlobalParameters.classifier_model_dir, args.sub_dir, args.classifier1_file_re))
    datasets = np.full(len(classifier_files), args.dataset1)
    if args.classifier2_file_re != "":
        classifier_files2 = \
            glob.glob(os.path.join(GlobalParameters.classifier_model_dir, args.sub_dir, args.classifier2_file_re))
        if args.dataset2 != "":
            datasets = np.append(datasets, np.full(len(classifier_files2), args.dataset2))
        else:
            datasets = np.append(datasets, np.full(len(classifier_files2), args.dataset1))
        classifier_files = classifier_files + classifier_files2

    plot_df = None
    tot_importance = {}
    for i, (clf_file, ds) in enumerate(zip(classifier_files, datasets)):

        if ds == args.dataset1:
            data, X, y = data1, X1, y1
        else:
            data, X, y = data2, X2, y2

        classifier_tag = os.path.basename(clf_file).split('_')[0]
        clf_mgr = get_clf_mgr(args.peptide_type, ds, classifier_tag, X, y)
        classifier = clf_mgr.load(clf_file)

        if classifier_tag in ['XGBoost']:
            explainer = shap.TreeExplainer(classifier)
        elif classifier_tag in ['LR']:
            explainer = shap.Explainer(classifier, X, feature_names=X.columns)

        shap_values = explainer(X)

        shap_values.feature_names = [GlobalParameters.plot_feature_names[f] for f in shap_values.feature_names]
        feature_importance = shap_values.abs.mean(0).values
        feature_importance /= sum(feature_importance)

        for fn, fi in zip(shap_values.feature_names, feature_importance):
            if fn not in tot_importance:
                tot_importance[fn] = 0
            tot_importance[fn] += fi

        if args.dataset2 != "" and args.dataset1 != args.dataset2:
            classifier_tag = classifier_tag + "_" + ds

        df = pd.DataFrame({'Feature': shap_values.feature_names, 'Feature importance': feature_importance,
                           'Classifier': classifier_tag, 'Replicate': clf_file})
        if plot_df is None:
            plot_df = df
        else:
            plot_df = pd.concat([plot_df, df], ignore_index=True)

    fig, ax = plt.subplots()
    fig.set_figheight(3*args.figure_height)
    fig.set_figwidth(args.figure_width)
    data_file = os.path.join(GlobalParameters.plot_dir, "{0}.txt".format(args.file_name))
    plot_df.to_csv(data_file, sep='\t', index=False, header=True)

    color_map = None
    if args.color_map != '':
        cm = ast.literal_eval(args.color_map)
        color_map = {}
        for key in cm:
            if type(cm[key]) == int:
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
    plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', transparent=True)
    plt.close()


