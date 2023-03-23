import numpy as np
import pandas as pd
import os

from Utils.DataManager import DataManager
from Utils.Parameters import Parameters


class NeoDiscImmunogenicityAnnotatorShort:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.immuno_data = self.read_immuno_info()
        self.response_cols = [c for c in self.immuno_data.columns if
                              (c.endswith('_actp_facs_type') or c.endswith('_prerep_facs_type'))]
        self.status_cols = [c for c in self.immuno_data.columns if
                            (c.endswith('_actp_status') or c.endswith('_prerep_status'))]
        self.neodisc_patients = set(self.immuno_data['patient_id'])

    def annotate_response_types(self, patients=None):
        if patients is None:
            patients = self.neodisc_patients

        for p in self.mgr.get_valid_patients('short'):
            if p in patients:
                data = self.annotate_patient(p)
                if data is not None:
                    out_file = os.path.join(Parameters().get_result_dir(), p+"_short_rt.txt")
                    data.to_csv(out_file, sep="\t", header=True, index=False)

    def read_immuno_info(self):
        excel_file = self.params.get_htide_immuno_info_file()
        with open(excel_file, 'rb') as file:
            data = pd.read_csv(file, header=0, sep="\t")
        return data

    def get_patient_info(self, patient, peptide_types, reactivity):
        mask = (self.immuno_data['patient_id'] == patient) & \
               (self.immuno_data.apply(lambda row: row['type'] in peptide_types, axis=1)) & \
               (self.immuno_data.apply(lambda row: row['reactivity'] in reactivity, axis=1))

        return self.immuno_data.loc[mask, :]

    def annotate_patient(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='short')
        if data is None:
            return None

        return self.update_response_type(patient, data)

    def update_response_type(self, patient, data):
        mutant_seqs = data['mutant_seq']
        genes = data['gene']
        p_info = self.get_patient_info(patient,
                                       ['MS (NEO-FS)', 'MS (NEO-SNV)', 'Predicted neo', 'Predicted Neo',
                                        'Predicted Neo (SNV)', 'Predicted Neo (FSS)', 'Predicted Neo (DELETION)',
                                        'Predicted Neo (INSERTION)', 'Predicted Neo (DNV)',
                                        'Predicted neo (deconvoluted)'],
                                       ['POSITIVE', 'NEGATIVE', 'POSITIVE;NEGATIVE'])
        rt_I, rt_annot = self.get_response_type(mutant_seqs, genes, p_info, 'CD8')
        data['response_type'] = rt_I
        data['response_annot'] = rt_annot

        return data

    def get_response_annot(self, row):
        row = np.array(row, dtype=str)
        contains_CD8 = 'CD8+' in row
        contains_CD4 = 'CD4+' in row
        annot = []
        if contains_CD8:
            annot.append("CD8+")
        if contains_CD4:
            annot.append("CD4+")
        if not contains_CD8 and not contains_CD4:
            annot.append("ND")

        return annot

    def is_CD8(self, row):
        row = np.array(row, dtype=str)
        contains_CD8 = 'CD8+' in row
        contains_CD4 = 'CD4+' in row
        if contains_CD8:
            return True
        elif contains_CD4:
            return False
        else:
            return True # for no annotation we assume positive test

    def get_response_type(self, mutant_seqs, genes, p_info, tag='CD8'):
        response_type = []
        response_annot = []
        for s, g in zip(mutant_seqs, genes):
            idx = (s == p_info['sequence']) & (g == p_info['gene'])
            if any(idx):
                df = p_info.loc[idx, ['reactivity'] + self.status_cols + self.response_cols]
                if any(df.apply(lambda r: ('POSITIVE' in r['reactivity']) &
                                          any(r[self.status_cols] == 'POSITIVE') &
                                          (self.is_CD8(r[self.response_cols])), axis=1)):
                    response_type.append(tag)
                    annot = []
                    df.apply(lambda r: annot.append(self.get_response_annot(r[self.response_cols]))
                             if ('POSITIVE' in r['reactivity']) & any(r[self.status_cols] == 'POSITIVE') else
                             None, axis=1)
                    annot = set(np.array(annot).flatten())
                    if 'CD8+' in annot:
                        annot = ['CD8+']
                    response_annot.append(",".join(annot))
#                    print("{0} {1} {}".format(s, g))
                else:
#                    print("{0} {1} NEGATIVE".format(s, g))
                    response_type.append("negative")
                    response_annot.append("")
            else:
#                print("{0} {1} not_tested".format(s, g))
                response_type.append("not_tested")
                response_annot.append("")

        return response_type, response_annot

    def get_patients(self):
        return self.neodisc_patients

    def update_immunogenicity_annotation(self, patient, file_tag_input, file_tag_output, mgr=DataManager()):
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, 'short')
        if data is None:
            return

        data = self.update_response_type(patient, data)

        out_file = os.path.join(parameters.get_result_dir(), patient+"_short_"+file_tag_output+".txt")
        data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

    def overwrite_hitide_immunogenicity_annotations(self, file_tag_input):
        for p in [p for p in self.neodisc_patients if p in self.mgr.get_valid_patients('short')]:
            self.update_immunogenicity_annotation(p, file_tag_input, file_tag_input)

