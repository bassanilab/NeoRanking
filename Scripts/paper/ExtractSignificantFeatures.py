from collections import Counter
from Utils.Parameters import *
from os import path
import pandas as pd

gartner_pv_file = path.join(Parameters().get_pickle_dir(), "PV_Gartner_short_03.24.2022-17.54.36_clf_results.txt")
tesla_pv_file = path.join(Parameters().get_pickle_dir(), "PV_TESLA_short_03.24.2022-17.53.08_clf_results.txt")
hitide_pv_file = path.join(Parameters().get_pickle_dir(), "PV_HiTIDE_short_03.24.2022-17.51.34_clf_results.txt")

col_names = ['feature', 'p-value', 'pass_mt']
gartner_pv = pd.read_csv(gartner_pv_file, header=None, names=col_names, sep="\t", skipfooter=3)
tesla_pv = pd.read_csv(tesla_pv_file, header=None, names=col_names, sep="\t", skipfooter=3)
hitide_pv = pd.read_csv(hitide_pv_file, header=None, names=col_names, sep="\t", skipfooter=3)

gartner_features = gartner_pv.loc[gartner_pv['p-value'] < 0.1, 'feature'].to_numpy()
tesla_features = tesla_pv.loc[tesla_pv['p-value'] < 0.1, 'feature'].to_numpy()
hitide_features = hitide_pv.loc[hitide_pv['p-value'] < 0.1, 'feature'].to_numpy()

features = []
for f in gartner_features: features.append(f)
for f in tesla_features: features.append(f)
for f in hitide_features: features.append(f)

counter = Counter(features)

for f in counter:
    print("{0}: {1}".format(f, counter[f]))

print(str(len(set(features)))+" features for short peptides")
print(' '.join(set(features)))

gartner_pv_file = path.join(Parameters().get_pickle_dir(), "PV_Gartner_long_03.26.2022-16.55.55_clf_results.txt")
tesla_pv_file = path.join(Parameters().get_pickle_dir(), "PV_TESLA_long_03.26.2022-16.57.22_clf_results.txt")
hitide_pv_file = path.join(Parameters().get_pickle_dir(), "PV_HiTIDE_long_03.26.2022-16.59.35_clf_results.txt")

col_names = ['feature', 'p-value', 'pass_mt']
gartner_pv = pd.read_csv(gartner_pv_file, header=None, names=col_names, sep="\t", skipfooter=3)
tesla_pv = pd.read_csv(tesla_pv_file, header=None, names=col_names, sep="\t", skipfooter=3)
hitide_pv = pd.read_csv(hitide_pv_file, header=None, names=col_names, sep="\t", skipfooter=3)

gartner_features = gartner_pv.loc[gartner_pv['p-value'] < 0.1, 'feature'].to_numpy()
tesla_features = tesla_pv.loc[tesla_pv['p-value'] < 0.1, 'feature'].to_numpy()
hitide_features = hitide_pv.loc[hitide_pv['p-value'] < 0.1, 'feature'].to_numpy()

features = []
for f in gartner_features: features.append(f)
for f in tesla_features: features.append(f)
for f in hitide_features: features.append(f)

counter = Counter(features)

for f in counter:
    print("{0}: {1}".format(f, counter[f]))

print(str(len(set(features)))+" features for long peptides")
print(' '.join(set(features)))
