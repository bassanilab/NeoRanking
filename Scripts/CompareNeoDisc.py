import numpy as np
import pandas as pd
nd_oym1_long = pd.read_csv("/home/localadmin/Priorization/results/NeoDisc/0YM1_long_netmhc_stab_chop_tap_mbp_tcr.txt",
                           sep="\t", header=0)
ml_oym1_long = pd.read_csv("/home/localadmin/Priorization/results/0YM1_long_rt_netmhc_stab_chop_tap_mbp_tcr.txt",
                           sep="\t", header=0)

print(list(zip(nd_oym1_long.columns, ml_oym1_long.columns)))
idx1 = np.where(nd_oym1_long.columns == 'mut_peptide_pos_0')[0][0]
idx2 = np.where(ml_oym1_long.columns == 'mut_peptide_pos_0')[0][0]

# for i, j in zip(range(idx1, nd_oym1_long.shape[1]), range(idx2, ml_oym1_long.shape[1])):
#     d1 = nd_oym1_long.iloc[:, i].describe()
#     d2 = ml_oym1_long.iloc[:, j].describe()
#     df = pd.DataFrame({'ML': d2, 'ND': d1})
# #    print("****************************")
#     print(nd_oym1_long.columns[i])
#     print(df)

nd_oym1_short = pd.read_csv("/home/localadmin/Priorization/results/NeoDisc/0YM1_short_stab_chop_tap_mbp_tcr.txt",
                            sep="\t", header=0)
ml_oym1_short = pd.read_csv("/home/localadmin/Priorization/results/0YM1_short_rt_stab_chop_tap_mbp_tcr.txt",
                            sep="\t", header=0)

idx = np.where(ml_oym1_short.columns == 'response_type')[0][0]

for i in range(idx+1, nd_oym1_short.shape[1]):
    d1 = nd_oym1_short.iloc[:, i-1].describe()
    d2 = ml_oym1_short.iloc[:, i].describe()
    df = pd.DataFrame({'ML': d2, 'ND': d1})
    print("****************************")
    print(ml_oym1_short.columns[i]+", "+nd_oym1_short.columns[i-1])
    print(df)
