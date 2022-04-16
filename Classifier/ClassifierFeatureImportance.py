import argparse
from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from sklearn.preprocessing import *
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from Visualization.PCAClassifyPeptideBrowser import *


parser = argparse.ArgumentParser(description='Calculate Feature importance for classifier')
parser.add_argument('-cl', '--classifier', type=str, help='Classifier pickle file')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, nargs='+', help='patient ids for training set')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-id', '--run_id', type=str, default='ML_SVM',
                    help='File tag for classifier pickle output file')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features used by classifier')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='q',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None)')
parser.add_argument('-nt', '--nr_train_patients', type=int, default=-1, help='Number of patients in -tr option considered')
parser.add_argument('-r', '--max_rank', type=int, default=20, help='Maximal rank for predicted immunogenic considered correct')
parser.add_argument('-mi', '--min_nr_immuno', type=int, default=1, help='Minimum nr of immunogenic mutations in sample')
parser.add_argument('-ni', '--nr_iter', type=int, default=30, help='Number of iteration in RandomSearchCV')
parser.add_argument('-nc', '--nr_classifiers', type=int, default=1, help='Number of best classifiers included for voting')
parser.add_argument('-cv', '--nr_cv', type=int, default=5, help='Number of CV layers in RandomSearchCV')
parser.add_argument('-rt', '--response_types', type=str, nargs='+', help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Coefficient alpha in score function')


args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

pickle_file = DataManager().get_classifier_file(args.run_id, clf_tag=args.classifier)

cat_features = [f for f in args.features if f in Parameters().get_categorical_features()]
optimizationParams = OptimizationParams(args.alpha, cat_features=cat_features)

normalizer = None
if args.normalizer == 'q':
    normalizer = QuantileTransformer()
    if args.verbose > 0:
        print('Normalizer: QuantileTransformer')
elif args.normalizer == 'z':
    normalizer = StandardScaler()
    if args.verbose > 0:
        print('Normalizer: StandardScaler')
else:
    if args.verbose > 0:
        print('Normalizer: None')


data_loader = DataLoader(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                         mutation_types=args.mutation_types, response_types=args.response_types,
                         immunogenic=args.immunogenic, min_nr_immuno=0)

# perform leave one out on training set
patients = np.array(args.patients_train)
if args.nr_train_patients > 1:
    patients = patients[0:min(args.nr_train_patients, len(patients))]

# fit best classifier on all data
data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag)
classifier = pickle.load(open(args.classifier, 'rb'))

feature_importance = np.zeros(len(args.features))
if classifier.__class__.__name__ == 'VotingClassifier':
    for clf in classifier.estimators_:
        clf.fit(X_train, y_train)
        if hasattr(clf, 'feature_importances_'):
            feature_importance = feature_importance + np.array(clf.feature_importances_)
        elif hasattr(clf, 'coef_'):
            feature_importance = feature_importance + np.abs(np.array(clf.coef_[0]))
else:
    classifier.fit(X_train, y_train)
    if hasattr(classifier, 'feature_importances_'):
        feature_importance = np.array(classifier.feature_importances_)
    elif hasattr(classifier, 'coef_'):
        feature_importance = np.abs(np.array(classifier.coef_[0]))

for i in np.arange(len(args.features)):
    print("{0}\t{1:.4f}".format(args.features[i], feature_importance[i]))


if args.verbose > 0 and args.pdf:
    pp = PdfPages(args.pdf)

    fig, ax = plt.subplots(figsize=(15, 8))
    x = np.arange(len(feature_importance))
    ax.bar(x, feature_importance)
    ax.set_xticks(x)
    ax.set_xticklabels(args.features)
    ax.xaxis.set_tick_params(labelsize=8, labelrotation=90)
    fig.tight_layout()

    pp.savefig(fig)

    learner = PrioritizationLearner("SVM", args.scorer, optimizationParams, verbose=0,
                                    nr_iter=args.nr_iter, nr_classifiers=args.nr_classifiers, nr_cv=args.nr_cv)

    y_pred, nr_correct, nr_immuno, r, mut_idx, score = \
        learner.test_classifier(classifier, 'All', X_train, y_train, max_rank=args.max_rank, prob=True)
    pca_browser = PCAClassifyPeptideBrowser(pd.DataFrame(data=X_train[:, 3:], columns=args.features[3:]), y_train, y_pred,
                                            "All datasets", True)
    pp.savefig(pca_browser.fig)

    df = pd.DataFrame(data=X_train[:, 3:8],
                      columns=['mut_Rank_BA_0', 'mut_Rank_BA_1', 'mut_Rank_BA_2', 'mut_Rank_BA_3', 'mut_Rank_BA_4'])
    df.insert(0, 'response', y_train)
    df = df.sort_values(by=['response'], ascending=True)
    g = sns.pairplot(df, hue='response')
    g.map_diag(sns.kdeplot, common_norm=False, palette='viridis')
    g.map_lower(sns.kdeplot, levels=5, palette='viridis', alpha=.5, common_norm=False, )
    pp.savefig(g.fig)

    df = pd.DataFrame(data=X_train[:, 8:13],
                      columns=['mut_Rank_Stab_0', 'mut_Rank_Stab_1', 'mut_Rank_Stab_2', 'mut_Rank_Stab_3',
                               'mut_Rank_Stab_4'])
    df.insert(0, 'response', y_train)
    df = df.sort_values(by=['response'], ascending=True)
    g = sns.pairplot(df, hue='response')
    g.map_diag(sns.kdeplot, common_norm=False, palette='viridis')
    g.map_lower(sns.kdeplot, levels=5, palette='viridis', alpha=.5, common_norm=False, )
    pp.savefig(g.fig)

    df = pd.DataFrame(data=X_train[:, 13:18],
                      columns=['mut_netchop_Ct_score_0', 'mut_netchop_Ct_score_1', 'mut_netchop_Ct_score_2',
                               'mut_netchop_Ct_score_3', 'mut_netchop_Ct_score_4'])
    df.insert(0, 'response', y_train)
    df = df.sort_values(by=['response'], ascending=True)
    g = sns.pairplot(df, hue='response')
    g.map_diag(sns.kdeplot, common_norm=False, palette='viridis')
    g.map_lower(sns.kdeplot, levels=5, palette='viridis', alpha=.5, common_norm=False, )
    pp.savefig(g.fig)

    df = pd.DataFrame(data=X_train[:, 18:23],
                      columns=['mut_netchop_Nt_score_0', 'mut_netchop_Nt_score_1', 'mut_netchop_Nt_score_2',
                               'mut_netchop_Nt_score_3', 'mut_netchop_Nt_score_4'])
    df.insert(0, 'response', y_train)
    df = df.sort_values(by=['response'], ascending=True)
    g = sns.pairplot(df, hue='response')
    g.map_diag(sns.kdeplot, common_norm=False, palette='viridis')
    g.map_lower(sns.kdeplot, levels=5, palette='viridis', alpha=.5, common_norm=False, )
    pp.savefig(g.fig)

    df = pd.DataFrame(data=X_train[:, 23:28],
                      columns=['mut_netchop_Int_score_0', 'mut_netchop_Int_score_1', 'mut_netchop_Int_score_2',
                               'mut_netchop_Int_score_3', 'mut_netchop_Int_score_4'])
    df.insert(0, 'response', y_train)
    df = df.sort_values(by=['response'], ascending=True)
    g = sns.pairplot(df, hue='response')
    g.map_diag(sns.kdeplot, common_norm=False, palette='viridis')
    g.map_lower(sns.kdeplot, levels=5, palette='viridis', alpha=.5, common_norm=False, )
    pp.savefig(g.fig)

    pp.close()

