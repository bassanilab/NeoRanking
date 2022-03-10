import numpy as np
from Utils.DataManager import *
from Utils.Parameters import *
import pandas as pd


class TESLAImmunogenicityAnnotatorLong:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.immuno_data = self.read_immuno_info()
        self.tesla_patients = set(self.immuno_data['PATIENT_ID'])

    def annotate_response_types(self):
        for p in self.mgr.get_valid_patients():
            if p.startswith('TESLA'):
                data = self.annotate_patient(p)
                if data is not None:
                    out_file = os.path.join(self.params.get_result_dir(), p + '_long_rt.txt')
                    data.to_csv(out_file, sep="\t", header=True, index=False)

    def annotate_patient(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='long')
        if data is None:
            return None

        patient_id = int(patient.replace("TESLA", ""))
        immuno_data_sel = self.immuno_data.loc[self.immuno_data['PATIENT_ID'] == patient_id, ]

        response_type = []
        for s, p in zip(data['mutant_seq'], data['pep_mut_start']):
            idx = TESLAImmunogenicityAnnotatorLong.intersect(s, p, immuno_data_sel['ALT_EPI_SEQ'])
            if len(idx) > 0:
                if any(immuno_data_sel.loc[immuno_data_sel.index[idx], 'VALIDATED'] == 1):
                    response_type.append('CD8')
                else:
                    response_type.append('negative')
            else:
                response_type.append('not_tested')

        data['response_type'] = response_type

        return data

    def read_immuno_info(self):
        file1, file2 = self.params.get_tesla_info_files()
        with open(file1, 'rb') as file:
            data1 = pd.read_excel(file, sheet_name='master-bindings-selected', header=0)
        with open(file2, 'rb') as file:
            data2 = pd.read_excel(file, sheet_name='SM_Table_S7_PEPTIDE_VALIDATION_', header=0)

        data1 = data1[['PATIENT_ID', 'ALT_EPI_SEQ', 'VALIDATED']]
        data2 = data2[['PATIENT_ID', 'ALT_EPI_SEQ', 'VALIDATED']]

        return pd.concat([data1, data2], ignore_index=True)

    @staticmethod
    def intersect(seq_long, mut_pos, seqs_short):
        def match(s_long, pos, s_short):
            idx = s_long.find(s_short)
            if idx < 0:
                return False
            else:
                return idx < pos <= idx + len(s_short)

        matched = \
            np.array(np.argwhere(list(map(lambda s: match(seq_long, mut_pos, s), seqs_short))))

        return matched.flatten()

    def get_patients(self):
        return set(['TESLA'+str(p) for p in self.tesla_patients])
