import argparse
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import shap
from sklearn.decomposition import PCA


parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-c', '--classifier', type=str, default='SVM', help='classifier to use')
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

args = parser.parse_args()

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


data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0, cat_to_num=args.cat_to_num)


data_train, X_train, y_train = data_loader.load_patients(args.patients_train, args.input_file_tag)

cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
if args.classifier == 'CatBoost':
    # Catboost does not support strings as feature names
    cat_features = list(range(len(cat_features)))

cat_idx = [i for i, f in enumerate(args.features) if f in Parameters().get_categorical_features()]

optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features, cat_idx=cat_idx,
                                        cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])

learner = PrioritizationLearner(args.classifier, args.scorer, optimizationParams, verbose=args.verbose,
                                nr_iter=args.nr_iter, nr_cv=args.nr_cv,
                                shuffle=args.shuffle, nr_epochs=args.nr_epoch, patience=args.early_stopping_patience,
                                batch_size=args.batch_size)
# perform leave one out on training set

best_score_train = -np.Inf
best_param_train = []
best_classifier_train = None
tot_correct_train = 0
tot_immunogenic_train = 0
tot_score_train = 0
tot_negative_train = 0

cvres, best_classifier, best_score, best_params = learner.optimize_classifier(X_train, y_train)

if args.classifier in ['XGBoost', 'CatBoost']:
    explainer = shap.Explainer(best_classifier)
elif args.classifier in ['SVM', 'SVM-lin']:
    explainer = shap.KernelExplainer(best_classifier.predict_proba,
                                     pd.DataFrame(X_train, columns=args.features),
                                     link="logit")
elif args.classifier in ['LR']:
    explainer = shap.Explainer(best_classifier, X_train, feature_names=args.features)

shap_values = explainer(pd.DataFrame(X_train, columns=args.features))
#shap_interaction_values = explainer.shap_interaction_values(X_train, columns=args.features)

with PdfPages(args.pdf) as pp:

    fig = plt.figure(figsize=(20, 10))
    shap.plots.scatter(shap_values[:, "mut_ic50_0"], color=shap_values[:,"GTEx_all_tissues_expression_median"], show=False)
    plt.title("Clf: {0}, training data".format(args.classifier))
    fig.tight_layout()
    pp.savefig()  # saves the current figure into a pdf page
    plt.close()

    fig = plt.figure(figsize=(20, 10))
    shap.plots.scatter(shap_values[:, "bestWTMatchScore_I"], color=shap_values[:, "GTEx_all_tissues_expression_median"], show=False)
    plt.title("Clf: {0}, training data".format(args.classifier))
    fig.tight_layout()
    pp.savefig()  # saves the current figure into a pdf page
    plt.close()
    #plt.savefig("../Plots/dependence_plot.pdf", dpi=400)
    #plt.show()
    fig = plt.figure(figsize=(30, 10))
    shap.plots.beeswarm(shap_values, max_display=30, show=False)
    plt.title(args.classifier)
    fig.tight_layout()
    pp.savefig()  # saves the current figure into a pdf page
    plt.close()

    fig = plt.figure(figsize=(30, 10))
    shap.plots.bar(shap_values, max_display=30, show=False)
    plt.title("Clf: {0}, training data".format(args.classifier))
    fig.tight_layout()
    pp.savefig()  # saves the current figure into a pdf page
    plt.close()

    if args.patients_test is not None:
        for p in args.patients_test:
            data_test, X_test, y_test = data_loader.load_patients(p, args.input_file_tag)
            mut_seq_idx = np.where(data_test.columns == 'mutant_seq')[0][0]
            gene_idx = np.where(data_test.columns == 'gene')[0][0]
            y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
                learner.test_classifier(best_classifier, p, X_test, y_test, max_rank=args.max_rank)
#            r = rankdata(-y_pred, method='average')-1
            shap_values = explainer(pd.DataFrame(X_test[y_test == 1], columns=args.features))
            for i in range(len(r)):
                fig = plt.figure(figsize=(10, 10))
                shap.plots.waterfall(shap_values[i], max_display=20, show=False)
                ttl = "Clf: {0}, Patient: {1}, mutation: {2}/{3}, rank={4}".\
                    format(args.classifier, p, data_test.iloc[mut_idx[i], mut_seq_idx],
                           data_test.iloc[mut_idx[i], gene_idx], r[i])
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
