import argparse
import os
import sys
import glob
import matplotlib.pyplot as plt
import shap
import seaborn as sns

from DataWrangling.DataTransformer import *
from Classifier.ClassifierManager import *

parser = argparse.ArgumentParser(description='Plot shapley values for a peptide')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory holding classifier model files')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-c', '--classifier_file_re', type=str, help='classifier files regular expression')
parser.add_argument('-tr', '--dataset_train', type=str, choices=GlobalParameters.datasets,
                    help='Dataset used for training, one of [NCI, NCI_train, NCI_test, TESLA, HiTIDE]')
parser.add_argument('-p', '--patient', type=str, help='Patient with immunogenic peptide')
parser.add_argument('-pept', '--peptide', type=str, help='Immunogenic peptide')
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
parser.add_argument('-ttp', '--title-prefix', type=str, default='', help='title prefix')

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

    
    data_train, X_train, y_train = load_data(args.dataset_train, args.peptide_type)

    classifier_files = \
        glob.glob(os.path.join(GlobalParameters.classifier_model_dir, args.sub_dir, args.classifier_file_re))

    classifiers = []
    classifier_tags = []
    explainers = []
    for i, clf_file in enumerate(classifier_files):
        classifier_tag = os.path.basename(clf_file).split('_')[0]
        clf_mgr = ClassifierManager(classifier_tag, 'sum_exp_rank', OptimizationParams(), verbose=0)
        classifier = clf_mgr.load(clf_file)
    
        if classifier_tag in ['XGBoost']:
            explainer = shap.TreeExplainer(classifier)
        elif classifier_tag in ['LR']:
            explainer = shap.Explainer(classifier, X_train, feature_names=X_train.columns)
    
        classifiers.append(classifier)
        classifier_tags.append(classifier_tag)
        explainers.append(explainer)

    data_p, X_p, y_p = \
        DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml', patient=args.patient,
                                          sample=False)

    if not any(data_p.mutant_seq.str.contains(args.peptide)):
        sys.exit("Peptide {0} not found in patient {}. Abort".format(args.peptide, args.patient))

    shap_values_repl = []
    ranks = []
    for i, (classifier_, classifier_tag_, explainer_) in enumerate(zip(classifiers, classifier_tags, explainers)):
        clf_mgr = ClassifierManager(classifier_tag_, 'sum_exp_rank', OptimizationParams(), verbose=0)
        y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
            clf_mgr.test_classifier(classifier=classifier_, peptide_type=args.peptide_type, patient=args.patient,
                                    data=data_p, x=X_p, y=y_p)
        shap_values = explainer_(X_p)
        shap_values.feature_names = [GlobalParameters.plot_feature_names[f] for f in shap_values.feature_names]
        shap_values_repl.append(shap_values)
        ranks.append(np.where(X_sorted.mutant_seq == args.peptide)[0][0]+1)

    plot_df_top20 = None
    for i in range(min(len(y_p), 20)):
        for j in range(len(shap_values_repl)):
            if plot_df_top20 is None:
                plot_df_top20 = pd.DataFrame({'Shapley value': shap_values_repl[j][i].values,
                                              'Feature': shap_values_repl[j][i].feature_names,
                                              'Peptide': 'Top 20'})
            else:
                plot_df_top20 = pd.concat([plot_df_top20,
                                           pd.DataFrame({'Shapley value': shap_values_repl[j][i].values,
                                                         'Feature': shap_values_repl[j][i].feature_names,
                                                         'Peptide': 'Top 20'})],
                                          ignore_index=True)

    plot_df = None
    i_p = np.where(data_p.mutant_seq == args.peptide)[0][0]
    for j in range(len(shap_values_repl)):
        if plot_df is None:
            plot_df = pd.DataFrame({'Shapley value': shap_values_repl[j][i_p].values,
                                    'Feature': shap_values_repl[j][i_p].feature_names,
                                    'Peptide': args.peptide})
        else:
            plot_df = pd.concat([plot_df, pd.DataFrame({'Shapley value': shap_values_repl[j][i_p].values,
                                                        'Feature': shap_values_repl[j][i_p].feature_names,
                                                        'Peptide': args.peptide})],
                                ignore_index=True)

    rank = np.mean(ranks)
    plot_df = pd.concat([plot_df, plot_df_top20], ignore_index=True)

    fig, ax = plt.subplots()
    fig.set_figheight(3*args.figure_height)
    fig.set_figwidth(args.figure_width)

    g = sns.barplot(data=plot_df, x='Shapley value', y='Feature', hue='Peptide')
    plt.title("{0}, {1}\n avg. rank: {2:.2f}".format(args.patient, args.peptide, rank), fontsize=args.title_size)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=args.tick_size)
    plt.ylabel("")
    plt.xlabel('Shapley value', fontsize=args.label_size)
    g.legend(fontsize=args.legend_size)
    sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=1,
                    title="", frameon=False, fontsize=args.legend_size, markerscale=10.0)
    plt.tight_layout()
    plot_file = os.path.join(GlobalParameters.plot_dir, "{0}.{1}".format(args.file_name, args.file_type))
    plt.savefig(plot_file, bbox_inches='tight', transparent=True)
    plt.close()
