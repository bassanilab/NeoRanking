import pandas as pd
from Utils.Parameters import *
import re

rosenberg_info_file = Parameters().get_rosenberg_info_file()
data_train = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 1', header=1)
data_test = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 18', header=1)


def add_tag(patient):
    return 'Rosenberg_' + str(patient)


def transform_alleles(alleles):
    return re.sub(r"HLA-([ABC])(\d{2}:\d{2})", r"\1*\2", alleles).replace(", ", ",")


Rosenberg_train_patients = data_train['ID'].transform(add_tag)
Rosenberg_test_patients = data_test['ID'].transform(add_tag)

print("Rosenberg train: [")
for p in Rosenberg_train_patients:
    print('\''+p+'\',')
print(']')

print("Rosenberg test: [")
for p in Rosenberg_test_patients:
    print('\''+p+'\',')
print(']')


df = data_train[['ID', 'Patient HLAs']]
df = df.append(data_test[['ID', 'Patient HLAs']], ignore_index=True)

df = df.drop(df[df.ID == 'Total'].index)

# print(add_tag(1000))
# print(add_tag('1000'))
# print(transform_alleles("HLA-A02:01, HLA-B27:05, HLA-B40:01, HLA-C02:02"))

df['ID'] = df['ID'].transform(add_tag)
df['Patient HLAs'] = df['Patient HLAs'].transform(transform_alleles)

df = df.rename(columns={'ID':  'Patient', 'Patient HLAs': 'Alleles'})

print(df.head())

hla_file = Parameters().get_allotype_file()
hla_data = pd.read_csv(filepath_or_buffer=hla_file, sep="\t", header=0)
hla_data = hla_data.append(df, ignore_index=True)

hla_data = hla_data.groupby('Patient')
hla_data = hla_data.first()
hla_data = hla_data.reset_index()

print(hla_data.head())
new_hla_file = Parameters().get_allotype_file()
new_hla_file = new_hla_file.replace(".txt", "_new.txt")
hla_data.to_csv(new_hla_file, sep='\t', header=True, index=False)
