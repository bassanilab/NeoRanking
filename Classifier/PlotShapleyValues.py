import argparse
import os
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import shap
from sklearn.decomposition import PCA

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import get_normalizer, get_valid_patients

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-c', '--classifier', type=str, default='', help='classifier to use')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, nargs='+', help='patient ids for test set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_SVM',
                    help='File tag for classifier pickle output file')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-r', '--max_rank', type=int, default=20, help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in HyperOpt search')
parser.add_argument('-cv', '--nr_cv', type=int, default=5, help='Number of CV layers in RandomSearchCV')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.005, help='Coefficient alpha in score function')
parser.add_argument('-sh', '--shuffle', dest='shuffle', action='store_true', help='Shuffle training data')
parser.add_argument('-e', '--nr_epoch', type=int, default=150, help='Number of epochs for DNN training')
parser.add_argument('-ep', '--early_stopping_patience', type=int, default=150,
                    help='Patience for early stopping for DNN training')
parser.add_argument('-b', '--batch_size', type=int, default=150, help='Batch size for DNN training')
parser.add_argument('-cat', '--cat_to_num', dest='cat_to_num', action='store_true', help='convert categories to numbers')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-nn', '--nr_negative', type=int, default=-1, help='Maximal number of short, non immunogenic samples')
parser.add_argument('-fp', '--feature_pairs', type=str, nargs='+', help='Features pair for pair plots')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                         immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num,
                         max_netmhc_rank=10000)

patients_train = \
    get_valid_patients(patients=args.patients_train, peptide_type=args.peptide_type) \
        if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

patients_test = \
    get_valid_patients(patients=args.patients_test, peptide_type=args.peptide_type) \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)
patients_test = patients_test.intersection(DataManager().get_immunogenic_patients(args.peptide_type))

data, X, y = data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type,
                                       nr_non_immuno_rows=args.nr_negative)

cat_features = [f for f in X.columns if f in Parameters().get_categorical_features()]
cat_idx = [X.columns.get_loc(col) for col in cat_features]

optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                        cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])


optimizationParams = \
    OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                       cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])

with open(args.classifier, mode='r') as result_file:
    classifier_tag = os.path.basename(args.classifier).split('_')[0]
    classifier = PrioritizationLearner.load_classifier(classifier_tag, optimizationParams, args.classifier)

if classifier_tag in ['XGBoost', 'CatBoost']:
    explainer = shap.Explainer(classifier)
elif classifier_tag in ['SVM', 'SVM-lin']:
    explainer = shap.KernelExplainer(classifier.predict_proba, X, link="logit")
elif classifier_tag in ['LR']:
    explainer = shap.Explainer(classifier, X, feature_names=args.features)

learner = PrioritizationLearner(classifier_tag, args.scorer, optimizationParams)

shap_values = explainer(X)
#shap_interaction_values = explainer.shap_interaction_values(X_train, columns=args.features)

with PdfPages(args.pdf) as pp:

    for fp in args.feature_pairs:
        (f1, f2) = fp.split(',')
        fig = plt.figure(figsize=(20, 10))
        shap.plots.scatter(shap_values[:, f1], color=shap_values[:, f2], show=False)
        plt.title("Clf: {0}, patients: {1}".format(classifier_tag, args.patients_train))
        fig.tight_layout()
        pp.savefig()  # saves the current figure into a pdf page
        plt.close()

    fig = plt.figure(figsize=(30, 10))
    shap.plots.beeswarm(shap_values, max_display=30, show=False)
    plt.title("Clf: {0}, patients: {1}".format(classifier_tag, args.patients_train))
    fig.tight_layout()
    pp.savefig()  # saves the current figure into a pdf page
    plt.close()

    fig = plt.figure(figsize=(30, 10))
    shap.plots.bar(shap_values, max_display=30, show=False)
    plt.title("Clf: {0}, patients: {1}".format(classifier_tag, args.patients_train))
    fig.tight_layout()
    pp.savefig()  # saves the current figure into a pdf page
    plt.close()

    if args.patients_test is not None:
        for p in patients_test:
            data, X, y = \
                    data_loader.load_patients(p, args.input_file_tag, args.peptide_type, verbose=False)
            y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
                learner.test_classifier(classifier, p, data, X, y, max_rank=args.max_rank)
            shap_values = explainer(X_sorted.loc[X_sorted['response'] == 1, X.columns])
            mut_seq_idx = np.where(X_sorted.columns == 'mutant_seq')[0][0]
            gene_idx = np.where(X_sorted.columns == 'gene')[0][0]
            mut_idx = np.where(X_sorted['response'] == 1)[0]
            for i in range(len(r)):
                fig = plt.figure(figsize=(10, 10))
                shap.plots.waterfall(shap_values[i], max_display=20, show=False)

                ttl = "Clf: {0}, Patient: {1}, mutation: {2}/{3}, rank={4}, score={5:.5f}".\
                    format(classifier_tag, p, X_sorted.iloc[mut_idx[i], mut_seq_idx],
                           X_sorted.iloc[mut_idx[i], gene_idx], r[i], y_pred_sorted.iloc[mut_idx[i]])
                plt.figtext(0.0, 0.99, ttl, fontsize=8)
                fig.tight_layout()
                pp.savefig()  # saves the current figure into a pdf page
                plt.close()

    # shap_values = explainer(pd.DataFrame(X_train, columns=args.features))
    # pca = PCA(n_components=2)
    #
    # x_pca = pca.fit_transform(pd.DataFrame(shap_values, columns=args.features))
    # variance = pca.explained_variance_ratio_
    #
    # x_pca_min = np.min(x_pca[:, 0])
    # x_pca_max = np.max(x_pca[:, 0])
    # y_pca_min = np.min(x_pca[:, 1])
    # y_pca_max = np.max(x_pca[:, 1])
    #
    # x_eps = abs(x_pca_max-x_pca_min)/100
    # y_eps = abs(y_pca_max-y_pca_min)/100
    #
    # x0_pca = np.min(x_pca[:, 0])
    # x1_pca = np.max(x_pca[:, 0])
    # y0_pca = np.min(x_pca[:, 1])
    # y1_pca = np.max(x_pca[:, 1])
    #
    # idx_neg = y_train == 0
    # idx_cd8 = y_train == 1
    #
    # plt.scatter(x_pca[idx_neg, 0], x_pca[idx_neg, 1], c='lightgray', s=5, picker=True, label='negative', alpha=0.3)
    # plt.scatter(x_pca[idx_cd8, 0], x_pca[idx_cd8, 1], c='gray', s=5, picker=True, label='negative top 20', alpha=0.7)
    # plt.title("PCA projection of {0} shap values".format(args.classifier), fontsize=10)
    # plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100))
    # plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100))
    # fig.tight_layout()
    # pp.savefig()  # saves the current figure into a pdf page
    # plt.close()
