import multiprocessing

from Utils.DataManager import *
from Utils.Parameters import *
import pandas as pd


class RosenbergImmunogenicityAnnotatorShort:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.gartner_data_train, self.gartner_data_test = self.read_gartner_info()
        self.gartner_patients_train = set(map(str, set(self.gartner_data_train['ID'])))
        self.gartner_patients_test = set(map(str, set(self.gartner_data_test['ID'])))

    def annotate_response_types(self):

        patients = set.intersection(set.union(self.gartner_patients_train, self.gartner_patients_test),
                                    self.mgr.get_valid_patients())

        for p in patients:
            self.annotate(p, True)
#        with multiprocessing.Pool(len(patients)) as pool:
#            [pool.apply_async(self.annotate, args=(p, True)) for p in patients]

    def annotate_row(self, mutant_seq, immuno_peptides, long_peptides):
        if any(mutant_seq == immuno_peptides):
            return 'CD8'
        elif any([mutant_seq in p for p in long_peptides]):
            return 'negative'
        else:
            return 'not_tested'

    def annotate(self, patient, write_res=False):
        data = self.mgr.get_original_data(patient, 'short')

        if data is None or data.shape[0] == 0:
            return

        if patient in self.gartner_patients_train:
            immuno_peptides = \
                self.gartner_data_train.loc[(self.gartner_data_train['ID'] == int(patient)) &
                                            (self.gartner_data_train['epitope status'] == 1.0), 'Peptide Mutant']
            long_peptides = \
                self.gartner_data_train.loc[self.gartner_data_train['ID'] == int(patient), 'Mut Epitope'].unique()
        elif patient in self.gartner_patients_test:
            immuno_peptides = \
                self.gartner_data_test.loc[(self.gartner_data_test['ID'] == int(patient)) &
                                           (self.gartner_data_test['epitope status'] == 1.0), 'Peptide Mutant']
            long_peptides = \
                self.gartner_data_test.loc[self.gartner_data_test['ID'] == int(patient), 'Mut Epitope'].unique()
        else:
            return None

        data['response_type'] = \
            data.apply(lambda row: self.annotate_row(row['mutant_seq'], immuno_peptides, long_peptides), axis=1)

        if write_res:
            out_file = os.path.join(self.params.get_result_dir(), patient + '_short_rt.txt')
            data.to_csv(out_file, sep="\t", header=True, index=False)

        return data

    def read_gartner_info(self):
        train_file, test_file = self.params.get_gartner_info_files(peptide_type='short')
        train_data = pd.read_csv(train_file, header=0, sep="\t")
        test_data = pd.read_csv(test_file, header=0, sep="\t")

        train_data = train_data[train_data['Screening Status'] != 'unscreened']
        test_data = test_data[test_data['Screening Status'] != 'unscreened']

        return train_data, test_data

    def get_patients(self, patient_subset='all'):
        patient_subset = patient_subset.lower()
        if patient_subset == 'all':
            return set.union(self.gartner_patients_train, self.gartner_patients_test)
        elif patient_subset == 'gartner_train':
            return self.gartner_patients_train
        elif patient_subset == 'gartner_test':
            return self.gartner_patients_test
        else:
            return None
