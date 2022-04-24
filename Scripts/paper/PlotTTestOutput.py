import argparse
import math

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Plot difference between t-test p-values')
parser.add_argument('-pdf', '--pdf', type=str, help='PDF output file')
parser.add_argument('-r', '--tt_result_files', type=str, nargs='+', help='list of t-test result files')
parser.add_argument('-n', '--names', type=str, nargs='+', help='Comma separated list of clf test names')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


class TTestResults:

    def __init__(self, lines, name):
        self.name = name
        self.results = self.parse_tt_results(lines)

    def parse_tt_results(self, lines):
        result_value_list = []
        header = None
        for l in lines:
            if l.startswith("Feature"):
                header = l.split("\t")
                continue
            if l.startswith("----"):
                break

            if header is not None:
                values = np.array(l.split("\t"))
                result_value_list.append(pd.Series(values))

        results = pd.concat(result_value_list, axis=1, ignore_index=True).transpose()
        results.columns = header
        results = results.astype({'pValue': 'float64'})
        results.rename(columns={'pValue': 'pValue'+"_"+self.name, 'rejected': 'rejected'+"_"+self.name},
                       inplace=True)
        return results

    def get_name(self):
        return self.name

    def get_results_data(self):
        return self.results

    def add_to_plot_dfs(self, plot_df):
        if plot_df is None:
            return self.results
        else:
            return plot_df.merge(self.results, on='Feature')


plot_df = None
for result_file, name in zip(args.tt_result_files, args.names):
    with open(result_file) as file:
        print(file)
        clf_results = TTestResults([line.rstrip() for line in file], name)
        plot_df = clf_results.add_to_plot_dfs(plot_df)

for c in plot_df.columns:
    if c.startswith('pValue'):
        plot_df[c].astype('float64')

with PdfPages(args.pdf) as pp:

    for i in range(len(args.names)):
        for j in range(i+1, len(args.names)):

            fig = plt.figure(figsize=(15, 15))
            fig.clf()
            g = sns.scatterplot(data=plot_df, x="pValue_"+args.names[i], y="pValue_"+args.names[j], alpha=0.7, s=300)
            for f in plot_df['Feature']:
                x = plot_df.loc[plot_df['Feature'] == f, "pValue_"+args.names[i]]
                y = plot_df.loc[plot_df['Feature'] == f, "pValue_"+args.names[j]]
                if abs(math.log10(x) - math.log10(y)) > 10:
                    plt.text(x, y, f, fontsize=20, horizontalalignment='center')
            min_pv = max(plot_df["pValue_"+args.names[i]].min(), plot_df["pValue_"+args.names[j]].min())
            g.plot([1, min_pv], [1, min_pv], color='r', linestyle='-', linewidth=2)
            plt.xlabel("pValue_"+args.names[i], size=20)
            plt.ylabel("pValue_"+args.names[j], size=20)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            g.set_xscale('log')
            g.set_yscale('log')
            g.figure.tight_layout()
            pp.savefig()
            plt.close()

