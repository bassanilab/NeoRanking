import numpy as np
import pandas as pd
import os
import difflib

from Utils.DataManager import DataManager
from Utils.Parameters import Parameters


class NeoDiscImmunogenicityAnnotatorLong:

    def __init__(self, mgr=None):
        self.params = Parameters()
        self.mgr = DataManager() if mgr is None else mgr
        self.immuno_data = self.read_immuno_info()
        self.response_cols = [c for c in self.immuno_data.columns if
                              (c.endswith('_actp_facs_type') or c.endswith('_prerep_facs_type'))]
        self.status_cols = [c for c in self.immuno_data.columns if
                            (c.endswith('_actp_status') or c.endswith('_prerep_status'))]
        self.neodisc_patients = set(self.immuno_data['patient_id'])

    def annotate_patient(self, patient):
        data = self.mgr.get_original_data(patient, peptide_type='long')
        if data is None:
            return None

        return self.update_response_type(patient, data)

    def update_response_type(self, patient, data):
        mutant_seqs = data['mutant_seq']
        mut_pos = data['pep_mut_start']
        genes = data['gene']
        p_info = self.get_patient_info(patient,
                                       ['MS (NEO-FS)', 'MS (NEO-SNV)', 'Predicted neo', 'Predicted Neo',
                                        'Predicted Neo (SNV)', 'Predicted Neo (FSS)', 'Predicted Neo (DELETION)',
                                        'Predicted Neo (INSERTION)', 'Predicted Neo (DNV)',
                                        'Predicted neo (deconvoluted)'],
                                       ['POSITIVE', 'NEGATIVE', 'POSITIVE;NEGATIVE'])
        rt_I, rt_annot = self.get_response_type(mutant_seqs, mut_pos, genes, p_info, 'CD8')
        data['response_type'] = rt_I
        data['response_annot'] = rt_annot

        return data

    def annotate_response_types(self, patients=None):
        if patients is None:
            patients = self.neodisc_patients

        for p in self.mgr.get_valid_patients('long'):
            if p in patients:
                data = self.annotate_patient(p)
                if data is None:
                    out_file = os.path.join(Parameters().get_result_dir(), p+"_long_rt.txt")
                    data.to_csv(out_file, sep="\t", header=True, index=False)

    def read_immuno_info(self):
        immo_file = self.params.get_htide_immuno_info_file()
        with open(immo_file, 'rb') as file:
            data = pd.read_csv(file, header=0, sep="\t")
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

    def get_patient_info(self, patient, peptide_types, reactivity):
        mask = (self.immuno_data['patient_id'] == patient) & \
               (self.immuno_data.apply(lambda row: row['type'] in peptide_types, axis=1)) & \
               (self.immuno_data.apply(lambda row: row['reactivity'] in reactivity, axis=1))

        return self.immuno_data.loc[mask, :]

    @staticmethod
    def intersect(agdisc_seq, agdisc_gene, neodisc_seq, neodisc_mut_pos, neodisc_gene):
        if agdisc_gene != neodisc_gene:
            return False
        s = difflib.SequenceMatcher(None, agdisc_seq, neodisc_seq)
        pos1, pos2, length = s.find_longest_match(alo=0, ahi=len(agdisc_seq), blo=0, bhi=len(neodisc_seq))
        if (pos2 <= (neodisc_mut_pos-1) <= pos2+length) and \
                (pos1+length == len(agdisc_seq) or pos2+length == len(neodisc_seq)):
            return True
        else:
            return False

    def get_response_type(self, mutant_seqs, mut_pos, genes, p_info, tag='CD8'):
        response_type = []
        response_annot = []
        for s, p, g in zip(mutant_seqs, mut_pos, genes):
            idx = p_info.apply(
                lambda row: NeoDiscImmunogenicityAnnotatorLong.intersect(row['sequence'], row['gene'], s, p, g), axis=1)
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

    def update_immunogenicity_annotation(self, patient, file_tag_input, file_tag_output, mgr=DataManager()):
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, 'long')
        if data is None:
            return

        data = self.update_response_type(patient, data)

        out_file = os.path.join(parameters.get_result_dir(), patient+"_long_"+file_tag_output+".txt")
        data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

    def overwrite_hitide_immunogenicity_annotations(self, file_tag_input):
        for p in [p for p in self.neodisc_patients if p in self.mgr.get_valid_patients('long')]:
            self.update_immunogenicity_annotation(p, file_tag_input, file_tag_input)

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
