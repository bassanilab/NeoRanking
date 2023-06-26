import argparse
import ast
import glob
import os

import pandas as pd

from Utils.GlobalParameters import GlobalParameters

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--classifier_dir', type=str, default=GlobalParameters().get_pickle_dir(),
                    help='directory with classifier files')
parser.add_argument('-c', '--classifier_file_re', type=str, default='', help='classifier files to use')
parser.add_argument('-o', '--output_file', type=str, default='', help='Output file with hyperparameters')
parser.add_argument('-cp', '--column_prefix', type=str, default='', help='Prefix for column')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))


def parse_hyperopt_results(train_file):
    with open(train_file) as file_:
        lines = [line for line in file_]
        for line in lines:
            if line.startswith("Hyperopt"):
                print(line)
                fields = line.replace("Hyperopt: ", "").split("; ")
                hyperopt_results_ = {}
                for f in fields:
                    k, v = f.split('=')
                    if k == 'Params':
                        hyperopt_results_[k] = ast.literal_eval(v)
                        if 'class_weight' in hyperopt_results_[k] and type(hyperopt_results_[k]['class_weight']) is dict:
                            hyperopt_results_[k]['class_weight'] = hyperopt_results_[k]['class_weight'][1]
                    else:
                        hyperopt_results_[k] = float(v)

                return hyperopt_results_


clf_result_files = glob.glob(os.path.join(args.classifier_dir, args.classifier_file_re))

hyperp_df = None
for i, file in enumerate(clf_result_files):
    hyperopt_results = parse_hyperopt_results(file)
    if hyperopt_results is None or len(hyperopt_results) == 0:
        continue
    if hyperp_df is None:
        hyperp_df = pd.Series(hyperopt_results['Params'], name="{0}_{1}".format(args.column_prefix, i)).to_frame()
    else:
        hyperp_df = \
            pd.concat([hyperp_df, pd.Series(hyperopt_results['Params'], name="{0}_{1}".format(args.column_prefix, i))],
                      axis=1)

hyperp_df.to_csv(args.output_file, header=True, index=True)
