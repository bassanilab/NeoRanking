import numpy as np
import pandas as pd
import os

from Utils.DataManager import DataManager
from Utils.Parameters import Parameters


class NeoDiscImmunogenicityAnnotatorLong:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.immuno_data = self.read_immuno_info()
        self.neodisc_patients = set(self.immuno_data['patient_id'])

    def annotate_patient(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='long')
        if data is None:
            return None
        mutant_seqs = data['mutant_seq']
        p_info = self.get_patient_info(patient, ['HLAI', 'HLA_I', 'HLA_I_II', 'HLA-I', 'HLA-I-II'],
                                                ['Predicted neo', 'Predicted Neo', 'Predicted neo (deconvoluted)',
                                                 'Predicted neo (long)', 'Predicted neo (SNV)'])
        rt_I = NeoDiscImmunogenicityAnnotatorLong.get_response_type(mutant_seqs, p_info, 'CD8')
        data['response_type'] = rt_I

        return data

    def annotate_response_types(self, patients=None):
        if patients is None:
            patients = self.neodisc_patients

        for p in self.mgr.get_valid_patients():
            if p in patients:
                data = self.annotate_patient(p)
                if data is None:
                    out_file = os.path.join(Parameters().get_result_dir(), p+"_long_rt.txt")
                    data.to_csv(out_file, sep="\t", header=True, index=False)

    def read_immuno_info(self):
        excel_file = self.params.get_htide_immuno_info_file()
        with open(excel_file, 'rb') as file:
            data = pd.read_excel(file, sheet_name='Feuil1')
        return data

    def get_patient_info(self, patient, hla, peptide_types):
        mask = (self.immuno_data['patient_id'] == patient) & \
               (self.immuno_data.apply(lambda row: row['type'] in peptide_types, axis=1)) & \
               ((self.immuno_data.apply(lambda row: row['hla_class'] in hla, axis=1)) |
                (self.immuno_data['hla_class'].isna()))

        return self.immuno_data.loc[mask, :]

    @staticmethod
    def intersect(seq1, seq2, peptide_type, offset=2):
        if peptide_type == 'Predicted neo (long)':
            return (seq1[offset:-offset] in seq2) | (seq2[offset:-offset] in seq1)
        else:
            return seq1 in seq2

    @staticmethod
    def get_response_type(mutant_seqs, p_info, tag='CD8'):
        response_type = []
        for s in mutant_seqs:
            idx = p_info.apply(
                lambda row: NeoDiscImmunogenicityAnnotatorLong.intersect(row['sequence'], s, row['type']), axis=1)
            if any(idx):
                df = p_info.loc[p_info.index[idx]]
                if any(df.apply(lambda r: 'POSITIVE' in r['reactivity'], axis=1)):
                    response_type.append(tag)
                elif any(df.apply(lambda r: 'NEGATIVE' == r['reactivity'], axis=1)):
                    response_type.append('negative')
                else:
                    response_type.append("not_tested")
            else:
                response_type.append("not_tested")

        return response_type

    # @staticmethod
    # def get_response_type(mutant_seqs, p_info, tag):
    #     not_tested = []
    #     is_immuno = []
    #     is_negative = []
    #     for s in mutant_seqs:
    #         bool1 = any(p_info.apply(
    #             lambda row: NeoDiscImmunogenicityAnnotatorLong.intersect(row['sequence'], s, row['type']), axis=1))
    #         not_tested.append(not bool1)
    #         if bool1:
    #             bool2 = any(p_info.apply(
    #                 lambda row: NeoDiscImmunogenicityAnnotatorLong.intersect(row['sequence'], s, row['type']) &
    #                             ('POSITIVE' in row['reactivity']), axis=1))
    #             is_immuno.append(bool2)
    #             if not bool2:
    #                 bool3 = any(p_info.apply(
    #                     lambda row: NeoDiscImmunogenicityAnnotatorLong.intersect(row['sequence'], s, row['type']) &
    #                                 ('NEGATIVE' == row['reactivity']), axis=1))
    #                 is_negative.append(bool3)
    #             else:
    #                 is_negative.append(False)
    #         else:
    #             is_immuno.append(False)
    #             is_negative.append(False)
    #
    #     rt = np.full(len(mutant_seqs), "not_tested")
    #     rt[is_negative] = "negative"
    #     rt[is_immuno] = tag
    #
    #     return rt

    def get_patients(self):
        return self.neodisc_patients
