from collections import Counter
import argparse

from Utils.DataManager import *

parser = argparse.ArgumentParser(description='Plot difference between t-test p-values')
parser.add_argument('-o', '--output_file', type=str, help='Text output file with significant features')
parser.add_argument('-r', '--tt_result_file', type=str, help='Comma separated list of t-test result files')
parser.add_argument('-mr', '--max_rank', type=int, default=10000, help='Maximal rank for features sorted by p-value')
parser.add_argument('-pv', '--max_pv', type=float, default=1, help='Maximal p-value for features')
parser.add_argument('-n', '--name', type=str, help='Comma separated list of clf test names')

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
with open(args.tt_result_file) as file:
    clf_results = TTestResults([line.rstrip() for line in file], args.name)
    plot_df = clf_results.add_to_plot_dfs(plot_df)


plot_df = plot_df.nsmallest(args.max_rank, 'pValue_Gartner_train')
features = plot_df.loc[plot_df['pValue_Gartner_train'] < 0.1, 'Feature'].to_numpy()

print(str(len(set(features)))+" features found twice in top {0}".format(args.max_rank))
print(' '.join(set(features)))

with open(args.output_file, mode='w') as result_file:
    result_file.write((' '.join(set(features))))
    result_file.write("\n")
