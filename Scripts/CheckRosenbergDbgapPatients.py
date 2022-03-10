import pandas as pd
from Utils.Parameters import *
from os import path
import glob


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
                    print("Patient " + str(patient) + " present in Gartner et al")
                elif any(data_parkhurst['patient']==patient):
                    print("Patient " + str(patient) + " present in Parkhurst et al")
                else:
                    print("Patient " + str(patient) + " not in Parkhurst or Gartner et al")



