from Utils.DataManager import *


class MutationAAChange:

    def __init__(self):
        """ Importance of non-binding amino acids for TCR recognition"""

        self.aa_prime_regr_coeff = {'A': -0.5, 'C': 0.5, 'D': -0.5, 'E': -4.2, 'F': 3.9, 'G': -1.75, 'H': -1.1,
                                    'I': 2.0, 'K': -3.0, 'L': 2.2, 'M': 1.15, 'N': 1.15, 'P': -1.5, 'Q': -3.15,
                                    'R': -2.85, 'S': -2.2, 'T': 1.85, 'V': 2.35, 'W': 7.9, 'Y': 2.95}

        self.binding_pos_fact = 1.0

        return

    def get_aa_coeff(self, aa, is_binding_pos):
        if aa not in self.aa_prime_regr_coeff:
            return np.nan

        if is_binding_pos:
            return self.binding_pos_fact*self.aa_prime_regr_coeff[aa]
        else:
            return self.aa_prime_regr_coeff[aa]

    def process_peptides(self, mut_aas, wt_aas, is_binding_pos):
        mut_coeffs = []
        wt_coeffs = []
        for i in np.arange(len(mut_aas)):
            mut_coeffs.append(self.get_aa_coeff(mut_aas[i], is_binding_pos[i]))
            wt_coeffs.append(self.get_aa_coeff(wt_aas[i], is_binding_pos[i]))

        return mut_coeffs, wt_coeffs

    @staticmethod
    def merge_with_data_long(data, mut_aa_coeff, wt_aa_coeff, i):
        data['mut_aa_coeff_'+str(i)] = mut_aa_coeff
        data['wt_aa_coeff_'+str(i)] = wt_aa_coeff
        return data

    @staticmethod
    def merge_with_data_short(data, mut_aa_coeff, wt_aa_coeff):
        data['mut_aa_coeff'] = mut_aa_coeff
        data['wt_aa_coeff'] = wt_aa_coeff
        return data

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):
        mgr = DataManager()
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            if 'mut_peptide_0' not in data.columns:
                print('Patient '+patient+' does not contain netmhc predictions. '
                                         'Unable to add MutationAAChange prediction.')
                return None
            data = self.add_features_long(data)
        else:
            if 'mut_is_binding_pos' not in data.columns:
                print('Patient '+patient+' does not contain binding site info. '
                                         'Unable to add MutationAAChange prediction.')
                return None
            data = self.add_features_short(data)

        mgr.put_processed_data(data, patient, file_tag_output, peptide_type)

        if write_res:
            out_file = os.path.join(parameters.get_result_dir(), patient+"_"+peptide_type+"_"+file_tag_output+".txt")
            data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

        return data

    def add_features_long(self, data):
        netMHCpan_ranks = [int(c[c.rfind('_')+1:]) for c in data.columns if 'mut_peptide_pos_' in c]
        for i in netMHCpan_ranks:
            mut_aas, wt_aas, is_binding_pos = MutationAAChange.get_peptide_info_long(data, i)
            mut_aa_coeff, wt_aa_coeff = self.process_peptides(mut_aas, wt_aas, is_binding_pos)
            data = MutationAAChange.merge_with_data_long(data, mut_aa_coeff, wt_aa_coeff, i)

        return data

    def add_features_short(self, data):
        mut_aas, wt_aas, is_binding_pos = MutationAAChange.get_peptide_info_short(data)
        mut_aa_coeff, wt_aa_coeff = self.process_peptides(mut_aas, wt_aas, is_binding_pos)
        data = MutationAAChange.merge_with_data_short(data, mut_aa_coeff, wt_aa_coeff)

        return data

    @staticmethod
    def get_peptide_info_long(data, index):
        mut_aa_column = 'aa_mutant'
        wt_aa_column = 'aa_wt'
        mut_is_binding = 'mut_is_binding_pos_'+str(index)
        if mut_is_binding not in data.columns:
            return None, None, None
        else:
            return np.array(data[mut_aa_column], dtype='str'), np.array(data[wt_aa_column], dtype='str'), \
                   np.array(data[mut_is_binding], dtype='bool')

    @staticmethod
    def get_peptide_info_short(data):
        mut_aa_column = 'aa_mutant'
        wt_aa_column = 'aa_wt'
        mut_is_binding = 'mut_is_binding_pos'
        if mut_is_binding not in data.columns:
            return None, None, None
        else:
            return np.array(data[mut_aa_column], dtype='str'), np.array(data[wt_aa_column], dtype='str'), \
                   np.array(data[mut_is_binding], dtype='bool')
