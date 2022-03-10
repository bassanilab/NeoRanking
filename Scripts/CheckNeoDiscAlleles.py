import pandas as pd
from Utils.Parameters import *
from os import path
import glob
import numpy as np


data_train = pd.read_csv(path.join(Parameters().get_data_dir(), "Rosenberg_train_info.txt"), header=0, sep="\t")
data_test = pd.read_csv(path.join(Parameters().get_data_dir(), "Rosenberg_test_info.txt"), header=0, sep="\t")

data = data_train.append(data_test, ignore_index=True)

data_parkhurst = pd.read_csv(path.join(Parameters().get_data_dir(),
                                       "../Rosenberg_data/RosenbergSamplesSummary_Preds_IEDB_GTEx_Cosmic_Intogen_data_ipMSDB.txt"),
                             header=0, sep="\t")


conf_files = glob.glob(path.join(Parameters().get_data_dir(),  "*.config"))

for f in conf_files:
    with open(f) as file:
        hlas_ref = []
        for line in file:
            if line.startswith("PATIENT_ID="):
                patient = line.rstrip().split("=")[1]
                patient = int(patient.replace("\"",""))
                if any(data['ID']==patient):
                    hlas_ref = data.loc[data['ID']==patient,'Patient HLAs']
                    hlas_ref = str(hlas_ref.iloc[0]).split(', ')
            if line.startswith("HLA_I_NETMHCPAN=") and  len(hlas_ref)>0:
                hlas = line.rstrip().split("=")[1].replace("\"","")
                hlas = hlas.split(",")

                cnt_match = 0
                for hla in hlas:
                    if hla in hlas_ref:
                        cnt_match += 1

                print("For patient {0} {1:d} alleles of a total of {2:d} match".format(patient, cnt_match, len(hlas)))
                print("NeoDisc alleles: {0}".format(hlas))
                print("Gartner et al. alleles {0}".format(hlas_ref))



