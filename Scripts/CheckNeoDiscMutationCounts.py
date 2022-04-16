import pandas as pd
from os import path
import glob
from Utils.Parameters import *


data_train = pd.read_csv(path.join(Parameters().get_data_dir(), "Rosenberg_train_info.txt"), header=0, sep="\t")
data_test = pd.read_csv(path.join(Parameters().get_data_dir(), "Rosenberg_test_info.txt"), header=0, sep="\t")

data = data_train.append(data_test, ignore_index=True)

neodisc_files = glob.glob(path.join(Parameters().get_data_dir(),  "*_long.txt"))

for f in neodisc_files:
    with open(f) as file:
        neodisc_data = pd.read_csv(file, header=0, sep="\t")
        peptide_id = neodisc_data.loc[0, 'peptide_id']
        patient = int(peptide_id.split('|')[0])
        if any(data['ID'] == patient):
            mut_cnt = data.loc[data['ID'] == patient, 'Total nmers']
            neodisc_mut_cnt = neodisc_data.shape[0]

            print("Patient {0}. NeoDisc mutation count={1}, Gartner et al. mutation count={2}".
                  format(patient, neodisc_mut_cnt, mut_cnt.iloc[0]))


mutations_train = pd.read_csv(path.join(Parameters().get_data_dir(), "Rosenberg_mutations_train.txt"), header=0, sep="\t")
mutations_test = pd.read_csv(path.join(Parameters().get_data_dir(), "Rosenberg_mutations_test.txt"), header=0, sep="\t")

mutation_data = mutations_train.append(mutations_test, ignore_index=True)

for i in mutation_data.index:
    patient = str(mutation_data.loc[i, 'ID'])
    mut_seq = mutation_data.loc[i, 'Mut Epitope']
    gene = mutation_data.loc[i, 'Gene Name']
    comment = mutation_data.loc[i, 'Comments']

    neodisc_data = None
    for nd_file in neodisc_files:
        if patient in nd_file:
            neodisc_data = pd.read_csv(nd_file, header=0, sep="\t")

    if neodisc_data is None:
        continue

    if not any(mut_seq[1:-1] in seq for seq in neodisc_data['mutant_seq']):
        print("Mutation {0} ({1}) of gene {2} in patient {3} not called by NeoDisc".format(mut_seq, comment, gene, patient))



