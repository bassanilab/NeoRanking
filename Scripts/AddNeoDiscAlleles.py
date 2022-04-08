# Markus Muller, CHUV, Switzerland
# Script to read HLA class I allotypes from NeoDisc .config files and write them into global file. This file contains
# a patient column and a column containing the class I alleles of a patients as commas separated list

import pandas as pd
from Utils.Parameters import *
from os import path
import glob

conf_files = glob.glob(path.join(Parameters().get_data_dir(), 'config', "*.config"))
hla_allele_file = Parameters().get_allotype_file()

if os.path.isfile(hla_allele_file):
    os.rename(hla_allele_file, hla_allele_file.replace(".txt", "_old.txt"))

hla_alleles = pd.DataFrame(columns=['Patient', 'Alleles'])
for f in conf_files:
    with open(f) as file:
        hlas = ""
        patient = os.path.basename(f).replace(".config", "")
        for line in file:
            if line.startswith("HLA_I="):
                hla_str = line.rstrip().split("=")[1].replace("\"", "")
                hlas = hla_str.split(',')
                hla_ar = []
                for hla in hlas:
                    hla_ar.append(hla[0]+"*"+hla[1:3]+":"+hla[3:5])
                break

        hla_alleles = hla_alleles.append(pd.Series([patient, ",".join(hla_ar)], index=hla_alleles.columns), ignore_index=True)


hla_alleles.to_csv(hla_allele_file, header=True, index=False, sep="\t")
