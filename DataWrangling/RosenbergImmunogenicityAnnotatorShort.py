import multiprocessing
import pandas as pd
import os

from Utils.DataManager import DataManager
from Utils.Parameters import Parameters


class RosenbergImmunogenicityAnnotatorShort:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.gartner_data_train, self.gartner_data_test = self.read_gartner_info()
        self.gartner_patients_train = set(map(str, set(self.gartner_data_train['ID'])))
        self.gartner_patients_test = set(map(str, set(self.gartner_data_test['ID'])))

    def annotate_response_types(self, patients=None):
        if patients is None:
            patients = set.intersection(set.union(self.gartner_patients_train, self.gartner_patients_test),
                                        self.mgr.get_valid_patients())

        if type(patients) is str:
            patients = [patients]

        for p in patients:
            data = self.annotate_patient(p)

            if data is not None:
                out_file = os.path.join(self.params.get_result_dir(), p + '_short_rt.txt')
                data.to_csv(out_file, sep="\t", header=True, index=False)

    def annotate_patient(self, patient):
        if patient in self.gartner_patients_train:
            return self.annotate(patient, 'train')
        elif patient in self.gartner_patients_test:
            return self.annotate(patient, 'test')
        else:
            return None

    def annotate_row(self, row, mut_id_seq_map, immuno_data, imm_peptides, neg_peptides):
        mut_seq = row['mutant_seq']
        if mut_seq not in imm_peptides and mut_seq not in neg_peptides:
            return 'not_tested'

        if row['mut_seqid'] in mut_id_seq_map:
            mut_seq_long = mut_id_seq_map[row['mut_seqid']]
        else:
            return 'no_mutation_found'

        screening_status = immuno_data.loc[(immuno_data['Peptide Mutant'] == row['mutant_seq']) &
                                           (immuno_data['Mut Epitope'] == mut_seq_long), 'epitope status']

        if any(screening_status == 1.0):
            return 'CD8'
        elif any(screening_status == 0.0):
            return 'negative'
        else:
            return 'not_tested'

    def map_mut_id_mut_seq(self, row, mut_id_seq_map):
        fields = row['peptide_id'].split('|')
        mut_id = fields[0]+":"+fields[2]
        mut_id_seq_map[mut_id] = row['mutant_seq']

    def convert_mit_seqid(self, row, patient):
        return patient+":"+row['mut_seqid']

    def annotate(self, patient, flag):
        data_short = self.mgr.get_original_data(patient, 'short')
        data_long = self.mgr.get_original_data(patient, 'long')

        if data_short is None or data_long is None:
            return None

        id_seq_map = {}
        data_long.apply(self.map_mut_id_mut_seq, args=(id_seq_map,), axis=1)
        data_short.loc[:, 'mut_seqid'] = data_short.apply(self.convert_mit_seqid, args=(patient,), axis=1)

        if flag == 'train':
            immuno_data = self.gartner_data_train.loc[self.gartner_data_train['ID'] == int(patient), ]
        else:
            immuno_data = self.gartner_data_test.loc[self.gartner_data_test['ID'] == int(patient), ]

        imm_peptides = \
            immuno_data.loc[immuno_data['epitope status'] == 1.0, 'Peptide Mutant'].unique()
        neg_peptides = \
            immuno_data.loc[immuno_data['epitope status'] == 0.0, 'Peptide Mutant'].unique()

        data_short['response_type'] = \
            data_short.apply(self.annotate_row, args=(id_seq_map, immuno_data, imm_peptides, neg_peptides), axis=1)

        return data_short

    def read_gartner_info(self):
        train_file, test_file = self.params.get_gartner_info_files(peptide_type='short')
        train_data = pd.read_csv(train_file, header=0, sep="\t")
        test_data = pd.read_csv(test_file, header=0, sep="\t")

        # here we select only the short peptides of mutated 25-mers that were screened with mingenes
        # only those peptides are immunogenic. There are 4 peptides that have 'Screening Status' == 'unscreened'
        # and 'epitope status' == 0.0 (i.e. negative). These are probably wrong annotations and we discard them
        train_data = train_data[train_data['Screening Status'] != 'unscreened']
        test_data = test_data[test_data['Screening Status'] != 'unscreened']

        return train_data, test_data

    def get_patients(self, patient_subset='all'):
        patient_subset = patient_subset.lower()
        if patient_subset == 'all' or patient_subset == 'nci' or patient_subset == 'gartner':
            return set.union(self.gartner_patients_train, self.gartner_patients_test)
        elif patient_subset == 'gartner_train' or patient_subset == 'nci_train':
            return self.gartner_patients_train
        elif patient_subset == 'gartner_test' or patient_subset == 'nci_test':
            return self.gartner_patients_test
        else:
            return None
