import pandas as pd

from Utils.DataManager import DataManager
from Utils.Parameters import Parameters


class RosenbergImmunogenicityAnnotatorLong:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.gartner_data_train, self.gartner_data_test = self.read_gartner_info()
        self.parkhurst_data = self.read_parkhurst_info()
        self.gartner_patients_train = set(map(str, set(self.gartner_data_train['ID'])))
        self.gartner_patients_test = set(map(str, set(self.gartner_data_test['ID'])))
        self.parkhurst_patients = set(map(str, set(self.parkhurst_data['ID'])))

    def annotate_response_types(self):
        for p in self.mgr.get_valid_patients():
            data = None
            if p in self.gartner_patients_train:
                data = self.annotate_gartner_train(p)
            elif p in self.gartner_patients_test:
                data = self.annotate_gartner_test(p)
            elif p in self.parkhurst_patients:
                data = self.annotate_parkhurst(p)

            if data is not None:
                out_file = os.path.join(self.params.get_result_dir(), p + '_long_rt.txt')
                data.to_csv(out_file, sep="\t", header=True, index=False)

    def annotate_gartner(self, patient):
        if patient in self.gartner_patients_train:
            return self.annotate_gartner_train(patient)
        elif patient in self.gartner_patients_test:
            return self.annotate_gartner_test(patient)

    def annotate_gartner_train(self, patient):
        data = self.mgr.get_original_data(patient, 'long')
        if data is None:
            return None

        immuno_data = self.gartner_data_train.loc[self.gartner_data_train['ID'] == int(patient), ]

        response_type = []
        for s in data['mutant_seq']:
            idx = RosenbergImmunogenicityAnnotatorLong.intersect(s, immuno_data['Mut Epitope'])
            if idx >= 0:
                if immuno_data.loc[immuno_data.index[idx], 'Screening Status'] == '-':
                    response_type.append('negative')
                elif immuno_data.loc[immuno_data.index[idx], 'Screening Status'] == 'CD8':
                    response_type.append('CD8')
                else:
                    response_type.append('not_tested')
            else:
                response_type.append('not_tested')

        data['response_type'] = response_type

        return data

    def annotate_gartner_test(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='long')
        if data is None:
            return None

        immuno_data = self.gartner_data_test.loc[self.gartner_data_test['ID'] == int(patient), ]

        response_type = []
        for s in data['mutant_seq']:
            idx = RosenbergImmunogenicityAnnotatorLong.intersect(s, immuno_data['Mut Epitope'])
            if idx >= 0:
                if immuno_data.loc[immuno_data.index[idx], 'Screening Status'] == '0':
                    response_type.append('negative')
                elif immuno_data.loc[immuno_data.index[idx], 'Screening Status'] == '1':
                    response_type.append('CD8')
                else:
                    response_type.append('not_tested')
            else:
                response_type.append('not_tested')

        data['response_type'] = response_type

        return data

    def annotate_parkhurst(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='long')
        if data is None:
            return None

        immuno_data = self.parkhurst_data.loc[self.parkhurst_data['ID'] == int(patient), ]

        response_type = []
        for s in data['mutant_seq']:
            idx = RosenbergImmunogenicityAnnotatorLong.intersect(s, immuno_data['Mutant nMER'])
            if idx >= 0:
                if immuno_data.loc[immuno_data.index[idx], 'CD8/CD4 Nmer screening results'] == 'negative':
                    response_type.append('negative')
                elif 'CD8' in immuno_data.loc[immuno_data.index[idx], 'CD8/CD4 Nmer screening results']:
                    response_type.append('CD8')
                elif 'CD4' in immuno_data.loc[immuno_data.index[idx], 'CD8/CD4 Nmer screening results']:
                    response_type.append('negative')
                else:
                    response_type.append('not_tested')
            else:
                response_type.append('not_tested')

        data['response_type'] = response_type

        return data

    def read_gartner_info(self):
        train_file, test_file = self.params.get_gartner_info_files()
        train_data = pd.read_csv(train_file, header=0, sep="\t")
        test_data = pd.read_csv(test_file, header=0, sep="\t")

        return train_data, test_data

    def read_parkhurst_info(self):
        excel_file = self.params.get_parkhurst_info_file()
        return pd.read_excel(open(excel_file, 'rb'), sheet_name='Supplementary 3', header=0)

    @staticmethod
    def match(seq1, seq2, offset):
        if offset <= 0:
            return (seq1 in seq2) | (seq2 in seq1)
        else:
            return (seq1[offset:-offset] in seq2) | (seq2[offset:-offset] in seq1)

    @staticmethod
    def intersect(seq, ref_seqs, offset=0):
        matched = \
            np.argwhere(list(map(lambda s: RosenbergImmunogenicityAnnotatorLong.match(s, seq, offset), ref_seqs)))

        return matched[0][0] if matched.shape[0] > 0 else -1

    def get_patients(self, patient_subset='all'):
        patient_subset = patient_subset.lower()
        if patient_subset == 'all':
            return set.union(self.gartner_patients_train, self.gartner_patients_test, self.parkhurst_patients)
        elif patient_subset == 'gartner_train':
            return self.gartner_patients_train
        elif patient_subset == 'gartner_test':
            return self.gartner_patients_test
        elif patient_subset == 'gartner':
            return set.union(self.gartner_patients_test, self.gartner_patients_train)
        elif patient_subset == 'parkhurst':
            return self.parkhurst_patients
        else:
            return None
