import argparse

parser = argparse.ArgumentParser(description='Print features and feature names to file')

parser.add_argument('-o', '--output_file', type=str, help='Output file name')
parser.add_argument('-f', '--features', type=str, nargs='+', help='Features')
parser.add_argument('-fd', '--feature_dict', type=str, nargs='+', help='Names of features used in plot')

args = parser.parse_args()

for arg in vars(args):
    print(arg, getattr(args, arg))

feature_dict = {}
for fn in args.feature_dict:
    (f, n) = fn.split(',')
    feature_dict[f] = n

with open(args.output_file, 'w') as file:
    for f in args.features:
        file.write("{0}\t{1}\n".format(f, feature_dict[f]))
