import pandas as pd
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import parallel_coordinates
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Plot and test difference between immunogenic and non immunogenic feature'
                                             'values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-r', '--clf_result_files', type=str, nargs='+', help='Comma separated list of clf result files')
parser.add_argument('-n', '--names', type=str, nargs='+', help='Comma separated list of clf test names')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


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
            if  len(fields) >  1:
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
                result_value_list.append(pd.Series(l.split("\t")))

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
            if item[0] not in plot_dfs:
                plot_dfs[item[0]] = pd.Series(item[1].split(','), name=self.name)
            else:
                plot_dfs[item[0]] = \
                    pd.concat([plot_dfs[item[0]], pd.Series(item[1].split(','), name=self.name)], axis=1)


plot_dfs = {}
for result_file, name in zip(args.clf_result_files, args.names):
    with open(result_file) as file:
        print(file)
        clf_results = ClassifierResults([line.rstrip() for line in file], name)
        clf_results.add_to_plot_dfs(plot_dfs)

plot_df = []
for patient in plot_dfs.keys():
    df = plot_dfs[patient]
    df.fillna(-1, inplace=True)
    df = df.astype('int32')
    max_rank = df.max().max()
    df[df < 0] = max_rank
    df['Patient'] = patient
    plot_dfs[patient] = df

plot_df = pd.concat(plot_dfs.values())

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
