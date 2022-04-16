import argparse

from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import parallel_coordinates
from matplotlib import pyplot as plt
import seaborn as sns

from Utils.Util_fct import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorShort import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorShort import *

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-r', '--clf_result_files', type=str, nargs='+', help='Comma separated list of clf result files')
parser.add_argument('-n', '--names', type=str, nargs='+', help='Comma separated list of clf test names')
parser.add_argument('-nd', '--neodisc', dest='neodisc', action='store_true', help='Include neodisc prioritization')
parser.add_argument('-pt', '--peptide_type', type=str, default='long', help='Peptide type (long or short)')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class NeoDiscResults:

    def __init__(self):
        self.mgr = DataManager()
        self.col_name = 'Neodisc'
        self.annotators = {}

    def get_annotator(self, peptide_type, patient_group):
        if peptide_type == 'long':
            if patient_group == 'Rosenberg':
                annotator = RosenbergImmunogenicityAnnotatorLong(self.mgr)
            elif patient_group == 'TESLA':
                annotator = TESLAImmunogenicityAnnotatorLong(self.mgr)
            else:
                annotator = NeoDiscImmunogenicityAnnotatorLong(self.mgr)
        else:
            if patient_group == 'Rosenberg':
                annotator = RosenbergImmunogenicityAnnotatorShort(self.mgr)
            elif patient_group == 'TESLA':
                annotator = TESLAImmunogenicityAnnotatorShort(self.mgr)
            else:
                annotator = NeoDiscImmunogenicityAnnotatorShort(self.mgr)

        return annotator

    def get_CD8_ranks(self, patient, peptide_type='long'):
        patient_group = get_patient_group(patient)
        if peptide_type not in self.annotators:
            self.annotators[peptide_type] = {}
        if patient_group not in self.annotators[peptide_type]:
            self.annotators[peptide_type][patient_group] = self.get_annotator(peptide_type, patient_group)

        annotator = self.annotators[peptide_type][patient_group]

        data = annotator.annotate_patient(patient)
        ranks = np.where(data['response_type'] == 'CD8')[0]
        return ranks+1

    def get_col_name(self):
        return self.col_name

    def add_to_plot_dfs(self, plot_dfs, patient, peptide_type='long'):
        if patient not in plot_dfs:
            plot_dfs[patient] = pd.Series(self.get_CD8_ranks(patient, peptide_type), name=self.col_name)
        else:
            plot_dfs[patient] = \
                pd.concat([plot_dfs[patient], pd.Series(self.get_CD8_ranks(patient, peptide_type), name='Neodisc')],
                          axis=1)


class ClassifierResults:

    def __init__(self, lines, name):
        self.config = {}
        self.parse_config(lines)
        self.results = self.parse_clf_results(lines)
        self.name = name

    def parse_config(self, lines):
        for l in lines:
            if l.startswith("Patient"):
                break
            fields = l.split("=")
            if len(fields) > 1:
                self.config[fields[0]] = fields[1]

    def parse_clf_results(self, lines):
        result_value_list = []
        header = None
        for l in lines:
            if l.startswith("Patient"):
                header = l.split("\t")
                header.append("Ranking_score")
                continue
            if l.startswith("nr_patients"):
                break

            if header is not None:
                values = np.array(l.split("\t"))
                result_value_list.append(pd.Series(values))

        results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
        results.columns = header
        return results

    def get_name(self):
        return self.name

    def get_config(self):
        return self.config

    def get_results_data(self):
        return self.results

    def add_to_plot_dfs(self, plot_dfs):
        for item in self.results.apply(lambda row: (row['Patient'], row['CD8_ranks']), axis=1):
            values = np.array(item[1].split(','), dtype='int32')
            for i in range(len(values)-1):
                d = values[i+1] - values[i]
                if d > 5:
                    values[i+1] = values[i+1] - int(d/3)
            values = pd.Series(values, name=self.name)

            if item[0] not in plot_dfs:
                plot_dfs[item[0]] = values
            else:
                plot_dfs[item[0]] = pd.concat([plot_dfs[item[0]], values], axis=1)

    def get_patients(self):
        return set(self.results['Patient'])


plot_dfs = {}
patients = set()
for result_file, name in zip(args.clf_result_files, args.names):
    with open(result_file) as file:
        print(file)
        clf_results = ClassifierResults([line.rstrip() for line in file], name)
        clf_results.add_to_plot_dfs(plot_dfs)
        patients = set.union(patients, clf_results.get_patients())

classifiers = args.names

if args.neodisc:
    neodisc_results = NeoDiscResults()
    col_name = neodisc_results.get_col_name()
    classifiers = classifiers.append(col_name)
    for p in patients:
        neodisc_results.add_to_plot_dfs(plot_dfs, p, args.peptide_type)

plot_df = []
mgr = DataManager()
for patient in plot_dfs.keys():
    df = plot_dfs[patient]
    df.fillna(-1, inplace=True)
    df = df.astype('int32')
    max_rank = df.max().max()
    df[df < 0] = max_rank
    df['Patient'] = patient
    patient_group = get_patient_group(patient)
    df['Patient_group'] = patient_group
    plot_dfs[patient] = df

plot_df = pd.concat(plot_dfs.values())
plot_df = plot_df.loc[plot_df['Patient_group'] == 'HiTIDE', args.names + ['Patient']]
plot_df_num = plot_df.loc[:, args.names]

top_20 = plot_df_num.apply(lambda c: sum(c < 20), axis=0)
top_50 = plot_df_num.apply(lambda c: sum(c < 50), axis=0)
top_100 = plot_df_num.apply(lambda c: sum(c < 100), axis=0)
med_rank = plot_df_num.apply(lambda c: c.median(), axis=0)
mean_rank = plot_df_num.apply(lambda c: c.mean(), axis=0)
exp_score_df = plot_df_num.transform(lambda v: np.exp(np.multiply(-0.02, v)), axis=1)
exp_scores = exp_score_df.apply(lambda c: sum(c), axis=0)

sum_plot_df = pd.concat([top_20, top_50, top_100], ignore_index=True)
nr = len(args.names) if args.neodisc else len(args.names)+1
sum_plot_df = pd.concat([sum_plot_df,
                         pd.Series(np.array([np.repeat('Top 20', nr), np.repeat('Top 50', nr),
                                             np.repeat('Top 100', nr)]).flatten())], axis=1, ignore_index=True)
sum_plot_df = pd.concat([sum_plot_df, pd.Series(np.array(np.tile(args.names, 3)))], axis=1, ignore_index=True)

sum_plot_df.columns = ['CD8+ 8-12 mers in top N', 'N', 'Classifier']
with PdfPages(args.pdf) as pp:
    patients = plot_df['Patient'].unique()
    for i in range(len(patients)):
        patient_sel = patients[i:min(i+1, len(patients))]
        df = plot_df[plot_df.Patient.isin(patient_sel)]
        fig = plt.figure(figsize=(10, 6))
        fig.clf()
        g = parallel_coordinates(df, 'Patient', color=['b'])
        plt.yscale('log')
        plt.ylabel("Rank", size=15)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=15)
        g.figure.tight_layout()
        pp.savefig()
        plt.close()

    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.barplot(x='Classifier', y='CD8+ 8-12 mers in top N', hue='N', data=sum_plot_df)
    plt.ylabel("CD8+ 8-12 mers in top N", size=15)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.barplot(x=args.names, y=exp_scores)
    plt.ylabel("Sum of exp(-rank/50)", size=15)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()

    plot_df = pd.melt(plot_df, id_vars='Patient', value_vars=args.names.append(''),
                      var_name='Classifier', value_name='CD8+ 8-12mers ranks')
    fig = plt.figure(figsize=(10, 6))
    fig.clf()
    g = sns.boxplot(x="Classifier", y="CD8+ 8-12mers ranks", data=plot_df, palette="Set3")
    plt.ylabel("CD8+ 8-12mers ranks", size=15)
    plt.yscale('log')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=15)
    g.figure.tight_layout()
    pp.savefig()
    plt.close()
