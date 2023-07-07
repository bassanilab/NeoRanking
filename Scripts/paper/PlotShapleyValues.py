import argparse
import os
from sklearn.preprocessing import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import shap
from sklearn.decomposition import PCA
import seaborn as sns
import time

from DataWrangling.DataTransformer import *
from Classifier.ClassifierManager import *
from Utils.Util_fct import get_normalizer, get_valid_patients

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--classifier_dir', type=str, default=GlobalParameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-c', '--classifier', type=str, default='', help='classifier to use')
parser.add_argument('-png', '--png_prefix', type=str, help='PNG output files prefix')
parser.add_argument('-s', '--scorer', type=str, default='sum_exp_rank', help='scorer for RandomSearchCV to use')
parser.add_argument('-tr', '--patients_train', type=str, help='patient ids for training set')
parser.add_argument('-te', '--patients_test', type=str, help='patient ids for test set')
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
encodings = read_cat_encodings(args.patients_train, args.peptide_type)

data_loader = DataTransformer(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                              mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative'],
                              immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                              cat_encoders=encodings, max_netmhc_rank=10000)

patients_train = \
    get_valid_patients(dataset=args.patients_train, peptide_type=args.peptide_type) \
        if args.patients_train and len(args.patients_train) > 0 else get_valid_patients(peptide_type=args.peptide_type)

patients_test = \
    get_valid_patients(dataset=args.patients_test, peptide_type=args.peptide_type) \
        if args.patients_test and len(args.patients_test) > 0 else get_valid_patients(peptide_type=args.peptide_type)
patients_test = patients_test.intersection(DataManager().get_immunogenic_patients(args.peptide_type))

data_train, X_train, y_train = data_loader.load_patients(patients_train, args.input_file_tag, args.peptide_type,
                                                         nr_non_immuno_rows=args.nr_negative)

cat_features = [f for f in X_train.columns if f in GlobalParameters().get_categorical_features()]
cat_idx = [X_train.columns.get_loc(col) for col in cat_features]

optimizationParams = \
    OptimizationParams(args.alpha, cat_idx=cat_idx, cat_dims=data_loader.get_categorical_dim(), input_shape=[len(args.features)])

classifier_tag = os.path.basename(args.classifier).split('_')[0]
classifier = ClassifierManager.load_classifier(classifier_tag, optimizationParams, args.classifier)
classifier.fit(X_train, y_train)

if classifier_tag in ['XGBoost', 'CatBoost']:
    explainer = shap.TreeExplainer(classifier)
elif classifier_tag in ['LR']:
    explainer = shap.Explainer(classifier, X_train, feature_names=args.features)

learner = ClassifierManager(classifier_tag, args.scorer, optimizationParams)

shap_values = explainer(X_train)

start = time.time()
if args.feature_pairs is not None:
    for fp in args.feature_pairs:
        (f1, f2) = fp.split(',')
        fig, ax = plt.subplots()
        fig.set_figheight(args.figure_height)
        fig.set_figwidth(args.figure_width)
        shap.plots.scatter(shap_values[:, f1], color=shap_values[:, f2], ax=ax, show=False)
        plt.title("Clf: {0}, patients: {1}".format(classifier_tag, args.patients_train), fontsize=args.title_size)
        plt.xlabel(feature_dict[f1], size=args.label_size)
        plt.ylabel("Shapley value", size=args.label_size)
        cmp = fig.get_axes()[1]
        cmp.set_ylabel(feature_dict[f2], size=args.label_size)
        png_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}_{2}_{3}.png".
                                format(args.png_prefix, args.patients_train, f1, f2))
        plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
        plt.close()

print("Feature pairs run in {0} sec.".format(time.time()-start))
start = time.time()
fig, ax = plt.subplots()
fig.set_figheight(3*args.figure_height)
fig.set_figwidth(args.figure_width)
fn = shap_values.feature_names
shap_values.feature_names = [feature_dict[f] for f in fn]
shap.plots.beeswarm(shap_values, max_display=len(args.features), show=False)
plt.title("Clf: {0}, patients: {1}".format(classifier_tag, args.patients_train), fontsize=args.title_size)
png_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}_beeswarm.png".
                        format(args.png_prefix, args.patients_train))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()

print("beeswarm run in {0} sec.".format(time.time()-start))
start = time.time()

fig, ax = plt.subplots()
fig.set_figheight(3*args.figure_height)
fig.set_figwidth(args.figure_width)
feature_importance = shap_values.abs.mean(0).values

df = pd.DataFrame({'Feature': shap_values.feature_names, 'Feature importance': feature_importance})
df = df.sort_values(by='Feature importance', ascending=False)
txt_file = os.path.join(GlobalParameters().get_plot_dir(),
                        "{0}_{1}_{2}_FeatureImp.txt".format(args.png_prefix, args.peptide_type, args.patients_train))
df.to_csv(txt_file, sep='\t', index=False, header=True)

sns.barplot(data=df, x='Feature importance', y='Feature', color='gray')
#shap.plots.bar(shap_values, max_display=len(args.features), show=False)
plt.title("Clf: {0}, patients: {1}".format(classifier_tag, args.patients_train), fontsize=args.title_size)
plt.xticks(fontsize=args.label_size)
plt.yticks(fontsize=args.tick_size)
plt.ylabel("")
plt.xlabel('mean(|Shap_value|)', fontsize=args.label_size)
png_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}_shap_mean.png".
                        format(args.png_prefix, args.patients_train))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


print("Feature importance run in {0} sec.".format(time.time()-start))
start = time.time()

mgr = DataManager()
test_ds = ['NCI_test', 'TESLA', 'HiTIDE']

if args.patients_test == 'debug':
    test_ds = ['debug']

# patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type)))
shapley_df = pd.DataFrame()
df_all = pd.DataFrame()

data_loader = DataTransformer(transformer=DataTransformer(), normalizer=normalizer, features=args.features,
                              mutation_types=args.mutation_types, response_types=['CD8', 'CD4/CD8', 'negative', 'not_tested'],
                              immunogenic=args.immunogenic, min_nr_immuno=0, cat_type=args.cat_encoder,
                              cat_encoders=encodings, max_netmhc_rank=10000)


for ds in test_ds:
    patients_test = get_valid_patients(dataset=ds, peptide_type=args.peptide_type)
    patients_test = sorted(patients_test.intersection(mgr.get_immunogenic_patients(args.peptide_type)))

    # data, X, y = data_loader.load_patients(patients_test, args.input_file_tag, args.peptide_type, verbose=False,
    #                                        nr_non_immuno_rows=1000)
    # shap_values = explainer(X.loc[y == 1, X.columns])
    # fn = shap_values.feature_names
    # shap_values.feature_names = [feature_dict[f] for f in fn]
    #
    # fig, ax = plt.subplots()
    # fig.set_figheight(2*args.figure_height)
    # fig.set_figwidth(args.figure_width)
    # shap.plots.heatmap(shap_values, max_display=len(args.features), show=False)
    # plt.title("Clf: {0}, dataset: {1}".format(classifier_tag, ds), fontsize=args.title_size)
    # png_file = os.path.join(Parameters().get_plot_dir(), "{0}_{1}_heatmap.png".format(args.png_prefix, ds))
    # plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
    # plt.close()
    # print("heatmap for {0} run in {1} sec.".format(ds, time.time()-start))
    # start = time.time()
    #
    df = pd.DataFrame()
    for p in patients_test:
        data, X, y = data_loader.load_patients(p, args.input_file_tag, args.peptide_type, verbose=False)

        if y is None or sum(y) == 0:
            continue
        y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
            learner.test_classifier(classifier_tag, classifier, p, data, X, y, max_rank=args.max_rank)
        shap_values = explainer(X)
        shap_values.feature_names = [feature_dict[f] for f in fn]

        print("shap for {0} in {1} run in {2} sec.".format(p, ds, time.time()-start))
        start = time.time()
        for i, i_t in enumerate(np.where(X_sorted['response'] == 1)[0]):
            s = pd.Series(X.iloc[i_t, :], index=X.columns)
            s['Rank'] = r[i]
            s['Size'] = X.shape[0]
            df = df.append(s, ignore_index=True)
            fig = plt.figure(figsize=(10, 10))
            shap.plots.waterfall(shap_values[i_t], max_display=len(args.features), show=False)
            ttl = "Clf: {0}, Patient: {1}, mutation: {2}, rank={3}, score={4:.5f}".\
                format(classifier_tag, p, X_sorted.loc[X_sorted.index[i_t], 'mutant_seq'], r[i], score)
            plt.title(ttl, fontsize=args.title_size)
            png_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_{1}_{2}_{3}_waterfall.png".
                                    format(args.png_prefix, p, ds, i))
            plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
            plt.close()

    df['Dataset'] = ds
    df_all = df_all.append(df, ignore_index=True)

df_all = df_all.sort_values(by=['Rank'], axis=0)
shap_sum = df_all.mean(axis=0)
sorted_f = [f for _, f in sorted(zip(shap_sum, df_all.columns), reverse=True)]
df_all = df_all[sorted_f+['Dataset']]
fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
#df_all = df_all.rename(columns=feature_dict)
y_lab = ["{0} {1}".format(ds, r) for ds, r in zip(df_all['Dataset'], df_all['Rank'])]
df_pca = df_all.drop(columns=['Rank', 'Dataset', 'Size'])
sns.heatmap(df_pca, yticklabels=y_lab, cmap="Spectral")
#plt.title("".format(classifier_tag), fontsize=args.title_size)
#plt.title("{0} shapley values".format(classifier_tag), fontsize=args.title_size)
plt.yticks(fontsize=5)
plt.xticks(rotation=90, fontsize=10)
png_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_shap_heatmap.png". format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()

pca = PCA(n_components=2)

shap_sum = df_pca.mean(axis=0)
df_pca = df_pca.iloc[:, np.where(shap_sum > 0)[0]]
x_pca = pca.fit_transform(df_pca)
variance = pca.explained_variance_ratio_

rank_cat = np.full(df_all.shape[0], "> 100")
rank_cat[df_all['Rank'] < 100] = "Top 100"
rank_cat[df_all['Rank'] < 50] = "Top 50"
rank_cat[df_all['Rank'] < 20] = "Top 20"
pca_df = pd.DataFrame({'PC1': x_pca[:, 0], 'PC2': x_pca[:, 1], 'Rank': df_all['Rank'],
                       'Rank Category': rank_cat, 'Dataset': df_all['Dataset']})

fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
sns.scatterplot(x='PC1', y='PC2',  hue="Rank Category", data=pca_df)
plt.title("PCA of {0} shapley values".format(classifier_tag), fontsize=args.title_size)
plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100), size=args.label_size)
plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100), size=args.label_size)
png_file = os.path.join(GlobalParameters().get_plot_dir(), "{0}_shap_pca.png". format(args.png_prefix))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()
