import numpy as np
import pandas as pd
import os

from Utils.DataManager import DataManager
from Utils.Parameters import Parameters


class NeoDiscImmunogenicityAnnotatorLongOld:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.immuno_data = self.read_immuno_info()
        self.neodisc_patients = set(self.immuno_data['Patient ID'])

    def annotate_response_types(self):

        for p in self.mgr.get_valid_patients():
            if p in self.neodisc_patients:
                data = self.mgr.get_original_data(p, 'long')
                if data is None:
                    continue
                mutant_seqs = data['mutant_seq']
                p_info = self.get_patient_info(p, ['HLAI', 'HLA-I', 'HLA-I-II'])
                rt_I = NeoDiscImmunogenicityAnnotatorLongOld.get_response_type(mutant_seqs, p_info, 'CD8')
                data['response_type'] = rt_I
                out_file = os.path.join(Parameters().get_result_dir(), p+"_long_rt.txt")
                data.to_csv(out_file, sep="\t", header=True, index=False)

    def read_immuno_info(self):
        excel_file = self.params.get_htide_immuno_info_file()
        with open(excel_file, 'rb') as file:
            data = pd.read_excel(file, sheet_name='Feuil1')
        return data

    def get_patient_info(self, patient, hla):
        mask = (self.immuno_data['Patient ID'] == patient) & \
               (self.immuno_data['Peptide type'].str.startswith('Predicted')) & \
               ((self.immuno_data.apply(lambda row: row['HLA Type'] in hla, axis=1)) |
                (self.immuno_data['HLA Type'].isna()))

        return self.immuno_data.loc[mask, :]

    @staticmethod
    def intersect(seq1, seq2, peptide_type, offset=2):
        if peptide_type == 'Predicted neo (long)':
            return (seq1[offset:-offset] in seq2) | (seq2[offset:-offset] in seq1)
        else:
            return seq1 in seq2

    @staticmethod
    def get_response_type(mutant_seqs, p_info, tag):
        not_tested = []
        is_immuno = []
        is_negative = []
        for s in mutant_seqs:
            bool1 = any(p_info.apply(
                lambda row: NeoDiscImmunogenicityAnnotatorLongOld.intersect(row['Sequence/mutant sequence'], s, row['Peptide type']),
                axis=1))
            not_tested.append(not bool1)
            if bool1:
                bool2 = any(p_info.apply(
                    lambda row: NeoDiscImmunogenicityAnnotatorLongOld.intersect(row['Sequence/mutant sequence'], s, row['Peptide type']) &
                                ((row['Status'] == 'Tested-positive') | (row['Status.1'] == 'Tested-positive')),
                    axis=1))
                is_immuno.append(bool2)
                bool2 = any(p_info.apply(
                    lambda row: NeoDiscImmunogenicityAnnotatorLongOld.intersect(row['Sequence/mutant sequence'], s, row['Peptide type']) &
                                ((row['Status'] == 'Tested-negative') | (row['Status.1'] == 'Tested-negative')),
                    axis=1))
                is_negative.append(bool2)
            else:
                is_immuno.append(False)
                is_negative.append(False)

        rt = np.full(len(mutant_seqs), "not_tested")
        rt[is_negative] = "negative"
        rt[is_immuno] = tag

        return rt

    def get_patients(self):
        return self.neodisc_patients
