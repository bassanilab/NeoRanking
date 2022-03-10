import numpy as np

from Utils.Parameters import *
from Utils.DataManager import *
import pandas as pd


class NeoDiscImmunogenicityAnnotatorShort:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.immuno_data = self.read_immuno_info()
        self.neodisc_patients = set(self.immuno_data['Patient ID'])

    def annotate_response_types(self):

        for p in self.mgr.get_valid_patients():
            if p in self.neodisc_patients:
                data = self.annotate_patient(p)
                if data is not None:
                    out_file = os.path.join(Parameters().get_result_dir(), p+"_short_rt.txt")
                    data.to_csv(out_file, sep="\t", header=True, index=False)

    def read_immuno_info(self):
        excel_file = self.params.get_htide_immuno_info_file()
        return pd.read_excel(open(excel_file, 'rb'), sheet_name='Feuil1')

    def get_patient_info(self, patient, hla, peptide_types):
        mask = (self.immuno_data['Patient ID'] == patient) & \
               (self.immuno_data.apply(lambda row: row['Peptide type'] in peptide_types, axis=1)) & \
               ((self.immuno_data.apply(lambda row: row['HLA Type'] in hla, axis=1)) |
                (self.immuno_data['HLA Type'].isna()))

        return self.immuno_data.loc[mask, :]

    def annotate_patient(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='short')
        if data is None:
            return None
        mutant_seqs = data['mutant_seq']
        p_info = self.get_patient_info(patient, ['HLAI', 'HLA-I', 'HLA-I-II'],
                                       ['Predicted neo', 'Predicted Neo', 'Predicted neo (deconvoluted)'])
        rt_I = NeoDiscImmunogenicityAnnotatorShort.get_response_type(mutant_seqs, p_info, 'CD8')
        data['response_type'] = rt_I

        return data

    @staticmethod
    def get_response_type(mutant_seqs, p_info, tag):
        response_type = []
        for s in mutant_seqs:
            idx = s == p_info['Sequence/mutant sequence']
            if any(idx):
                if any((p_info.loc[p_info.index[idx], 'Status'] == 'Tested-positive') |
                       (p_info.loc[p_info.index[idx], 'Status.1'] == 'Tested-positive')):
                    response_type.append('CD8')
                else:
                    response_type.append('negative')
            else:
                response_type.append("not_tested")

        return response_type

    def get_patients(self):
        return self.neodisc_patients
