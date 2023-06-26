import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression
from matplotlib import patches

from DataWrangling.DataTransformer import DataTransformer
from DataWrangling.CatEncoder import DataTransformer
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Plot correlation between features')

parser.add_argument('-fp', '--file_prefix', type=str, default="Feature_pair", help='PNG output files prefix')
parser.add_argument('-ft', '--file_type', type=str, default="svg", help='File type for plot (png, svg or pdf')
parser.add_argument('-ds', '--datasets', type=str, nargs='+', help='Datasets used')
parser.add_argument('-i', '--input_file_tag', type=str, default='netmhc_stab_chop',
                    help='File tag for neodisc input file (patient)_(input_file_tag).txt')
parser.add_argument('-fpair', '--feature_pairs', type=str, help='Features pair for pair plots')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')
parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of reporting')
parser.add_argument('-n', '--normalizer', type=str, default='n',
                    help='Normalizer used by classifier (q: quantile, z: standard, n: None or dictionary)')
parser.add_argument('-rt', '--response_types', type=str, help='response types included')
parser.add_argument('-mt', '--mutation_types', type=str, nargs='+', help='mutation types included')
parser.add_argument('-im', '--immunogenic', type=str, nargs='+', help='immunogenic response_types included')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')
parser.add_argument('-lep', '--legend_position', type=str, default='best', help='Legend position in plot')
parser.add_argument('-ci', '--color_immunogenic', type=str, default='darkorange', help='Color of immunogenic peptides')
parser.add_argument('-cn', '--color_negative', type=str, default='royalblue', help='Color of negative peptides')
parser.add_argument('-las', '--label_size', type=str, default='x-large',
                    help='Axis label size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-tis', '--tick_size', type=str, default='large',
                    help='Tick size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-tts', '--title_size', type=str, default='large', help='Title font size')
parser.add_argument('-les', '--legend_size', type=str, default='large',
                    help='Legend size, either float or one of: xx-small, x-small, small, medium, large, x-large, '
                         'xx-large, larger, or smaller')
parser.add_argument('-fiw', '--figure_width', type=float, default=25.0, help='Figure width in inches')
parser.add_argument('-fih', '--figure_height', type=float, default=25.0, help='Figure height in inches')
parser.add_argument('-rot', '--rotation', type=float, default=0.0, help='x-axis label rotation')
parser.add_argument('-rf', '--rotate_labels', type=str, nargs='+', help='Features with x-label rotation')
parser.add_argument('-dpi', '--resolution', type=float, default=200, help='Figure resolution in dots per inch')

args = parser.parse_args()

if args.verbose > 0:
    for arg in vars(args):
        print(arg, getattr(args, arg))

normalizer = get_normalizer(args.normalizer)

try:
    rt = ast.literal_eval(args.response_types)
except:
    print('Cannot parse normalization dictionary {}'.format(args.response_types))


features = []
(f1, f2) = args.feature_pairs.split(',')
features.append(f1)
features.append(f2)

features = np.unique(features)

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

df_plot = pd.DataFrame()
p_values = {}
txt_coord = {}
ttl = ""
if args.peptide_type == 'short':
    imm_label = 'Neo-pep_imm'
    neg_label = 'Neo-pep_non-imm'
else:
    imm_label = 'Mut-seq_imm'
    neg_label = 'Mut-seq_non-imm'


for ds in args.datasets:
    patients = get_valid_patients(ds)

    data_loader = DataTransformer(transformer=DataTransformer(), normalizer=normalizer, features=features,
                                  mutation_types=args.mutation_types, response_types=rt[ds],
                                  immunogenic=args.immunogenic, min_nr_immuno=0)

    data_train, X_train, y_train = data_loader.load_patients(patients, args.input_file_tag, args.peptide_type)

    X_1 = X_train.loc[y_train == 1, :]
    X_1['response'] = "1"
    X_1['Frequency'] = 1
    X_0 = X_train.loc[y_train == 0, :]
    X_0['response'] = "0"
    X_0['Frequency'] = 1
    df = pd.concat([X_1, X_0])
    all_num_cols = GlobalParameters().get_numerical_features()
    num_cols = [c for c in df.columns if c in all_num_cols]
    df[num_cols] = df[num_cols].apply(pd.to_numeric, errors='coerce')
    df_small = df.head(1000).sort_values(by=['response'])

    (f1, f2) = args.feature_pairs.split(',')
    if type(normalizer) == dict:
        norm_f1 = normalizer[f1]
        normalizer_name = get_normalizer_name(norm_f1)
        norm_f2 = normalizer[f2]
        normalizer_name = get_normalizer_name(norm_f2)
    else:
        norm_f1 = normalizer
        normalizer_name = get_normalizer_name(norm_f1)
        norm_f2 = normalizer
        normalizer_name = get_normalizer_name(norm_f2)

    if f1 in GlobalParameters().get_numerical_features() and f2 in GlobalParameters().get_numerical_features():
        lin_reg = LinearRegression()
        lin_reg.fit(X_0[f1].to_numpy().reshape(-1, 1), X_0[f2])
        res_0 = np.subtract(X_0[f2].to_numpy(), lin_reg.predict(X_0[f1].to_numpy().reshape(-1, 1)))
        res_1 = np.subtract(X_1[f2].to_numpy(), lin_reg.predict(X_1[f1].to_numpy().reshape(-1, 1)))

        mwt = stats.mannwhitneyu(res_0, res_1, alternative='less')
        p_values[ds] = mwt.pvalue
        txt_coord[ds] = max(max(res_1), max(res_0))
        if len(ttl) > 0:
            ttl += ","
        ttl = ttl + "{0}: {1:.3e}".format(ds, mwt.pvalue)

        df_plot = df_plot.append(pd.DataFrame({'Residual': res_0, 'Subset': neg_label, 'Dataset': ds}), ignore_index=True)
        df_plot = df_plot.append(pd.DataFrame({'Residual': res_1, 'Subset': imm_label, 'Dataset': ds}), ignore_index=True)


fig = plt.figure()
fig.set_figheight(args.figure_height)
fig.set_figwidth(args.figure_width)
g = sns.boxplot(data=df_plot, x="Dataset", y="Residual", hue='Subset', hue_order=[neg_label, imm_label],
                palette={neg_label: args.color_negative, imm_label: args.color_immunogenic}, notch=True)
handles = [patches.Patch(color=args.color_immunogenic, label=imm_label, alpha=0.7),
           patches.Patch(color=args.color_negative, label=neg_label, alpha=0.7)]
sns.move_legend(g, loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=1, handles=handles,
                title="", frameon=False, fontsize=args.legend_size, markerscale=10.0)
plt.ylim(df_plot["Residual"].min()*1.05, df_plot["Residual"].max()*1.5)

for i, ds in enumerate(args.datasets):
    plt.text(i-0.2, txt_coord[ds]*1.05, "{0:.3e}".format(p_values[ds]))

plt.xticks(fontsize=args.tick_size)
plt.xlabel('')
plt.ylabel('Regression residuals', size=args.label_size)

png_file = os.path.join(GlobalParameters().get_plot_dir(),
                        "{0}_{1}_{2}_{3}.{4}".format(args.file_prefix, "_".join(args.datasets), f1, f2, args.file_type))
plt.savefig(png_file, bbox_inches='tight', dpi=args.resolution)
plt.close()


