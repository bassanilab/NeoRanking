from pathlib import Path

import numpy as np

from Utils.DataManager import *


class SequenceLogoMgr:

    def __init__(self, seq_logo_dir, binding_threshold=15):
        """ seq_logo_dir is directory where allele_list.txt file is found and under which the
        sequence logos are stored"""

        self.seq_logo_dir = seq_logo_dir
        self.primary_alleles = []
        self.sequence_logos = {}
        self.allele_map = {}
        self.binding_threshold = binding_threshold

        self.read_primary_alleles()
        self.read_secondary_alleles()

        return

    def read_primary_alleles(self):

        data = pd.read_csv(list(Path(self.seq_logo_dir).rglob("allele_list.txt"))[0], sep="\t", header=0)
        self.primary_alleles = np.unique(data['Allele'])

        for a in self.primary_alleles:
            self.sequence_logos[a] = SequenceLogo(self.seq_logo_dir, a, binding_threshold=self.binding_threshold)

        return

    def read_secondary_alleles(self):

        data = pd.read_csv(list(Path(self.seq_logo_dir).rglob("alleles_mapping.txt"))[0], sep="\t", header=None)

        for i in data.index.values:
            self.allele_map[data.loc[i, 0]] = data.loc[i, 1]
            self.sequence_logos[data.loc[i, 0]] = self.sequence_logos[data.loc[i, 1]]

        return

    def get_sequence_logo(self, allele):
        allele = allele.replace(":", "").replace("*",  "").replace("HLA-", "")
        if allele in self.sequence_logos:
            return self.sequence_logos[allele]
        else:
            return None

    def is_binding_pos(self, allele, length, pos):
        allele = allele.replace(":", "").replace("*",  "").replace("HLA-", "")
        if allele in self.sequence_logos:
            return self.sequence_logos[allele].is_binding_pos(length, pos)
        else:
            return np.nan

    def get_aa_score(self, allele, length, aa, pos):
        allele = allele.replace(":", "").replace("*", "").replace("HLA-", "")
        if allele in self.sequence_logos:
            return self.sequence_logos[allele].get_aa_score(length, aa, pos)
        else:
            return np.nan

    def process_peptides(self, peptides, aas, alleles, mutation_positions):
        is_binding_pos = []
        binding_pos_score = []

        for i in np.arange(len(aas)):
            peptide = str(peptides[i])
            if not (8 <= len(peptide) <= 12 and 1 <= mutation_positions[i] <= len(peptide) and len(aas[i]) == 1):
                is_binding_pos.append(np.nan)
                binding_pos_score.append(np.nan)
                continue

            # print(f"{peptides[i]}, {aas[i]}, {alleles[i]}, {mutation_positions[i]}")
            allele = alleles[i]
            mut_pos = mutation_positions[i]
            pept_len = len(peptide)
            mut_aa = aas[i]

            if ',' in allele:
                allele_list = str(allele).split(',')
                is_binding = np.nan
                best_score = np.nan
                for a in allele_list:
                    score = self.get_aa_score(a, pept_len, mut_aa, mut_pos)
                    if not np.isnan(score) and (np.isnan(best_score) or score > best_score):
                        is_binding = self.is_binding_pos(a, pept_len, mut_pos)
                        best_score = score

                is_binding_pos.append(is_binding)
                binding_pos_score.append(best_score)
            else:
                is_binding_pos.append(self.is_binding_pos(allele, pept_len, mut_pos))
                binding_pos_score.append(self.get_aa_score(allele, pept_len, mut_aa, mut_pos))

        return binding_pos_score, is_binding_pos

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):
        mgr = DataManager()
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            if 'mut_peptide_0' not in data.columns:
                print('Patient '+patient+' does not contain netmhc predictions. Unable to add binding info.')
                return None
            data = self.add_features_long(data)
        else:
            if 'mutant_best_alleles_netMHCpan' not in data.columns:
                print('Patient '+patient+' does not contain netmhc prediction. Unable to add binding info.')
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
            mut_peptides, mut_aas, wt_aas, mut_alleles, mut_positions = SequenceLogoMgr.get_peptide_info_long(data, i)
            mut_binding_score, mut_is_binding_pos = \
                self.process_peptides(mut_peptides, mut_aas, mut_alleles, mut_positions)
            wt_binding_score, mut_is_binding_pos = \
                self.process_peptides(mut_peptides, wt_aas, mut_alleles, mut_positions)
            data = SequenceLogoMgr.merge_with_data_long(data, mut_is_binding_pos, mut_binding_score, wt_binding_score, i)

        return data

    def add_features_short(self, data):
        mut_peptides, mut_aas, mut_alleles, mut_positions = SequenceLogoMgr.get_peptide_info_short(data)
        mut_binding_score, mut_is_binding_pos = \
            self.process_peptides(mut_peptides, mut_aas, mut_alleles, mut_positions)
        data = SequenceLogoMgr.merge_with_data_short(data, mut_is_binding_pos, mut_binding_score)

        return data

    @staticmethod
    def merge_with_data_long(data, mut_is_binding_pos, mut_binding_score, wt_binding_score, i):
        data['mut_is_binding_pos_'+str(i)] = mut_is_binding_pos
        data['wt_binding_score_'+str(i)] = wt_binding_score
        data['mut_binding_score_'+str(i)] = mut_binding_score

        return data

    @staticmethod
    def merge_with_data_short(data, mut_is_binding_pos, mut_binding_score):
        data['mut_is_binding_pos'] = mut_is_binding_pos
        data['mut_binding_score'] = mut_binding_score

        return data

    @staticmethod
    def get_peptide_info_long(data, index):
        mut_peptide = 'mut_peptide_'+str(index)
        mut_aa_column = 'aa_mutant'
        wt_aa_column = 'aa_wt'
        mut_allele_column = 'mut_allele_'+str(index)
        mut_pos_column = 'mut_pos_in_peptide_'+str(index)
        if mut_pos_column not in data.columns:
            return None, None, None, None, None
        else:
            return np.array(data[mut_peptide], dtype='str'), np.array(data[mut_aa_column], dtype='str'), \
                   np.array(data[wt_aa_column], dtype='str'), np.array(data[mut_allele_column], dtype='str'), \
                   np.array(data[mut_pos_column], dtype='int32')

    @staticmethod
    def get_peptide_info_short(data):
        mut_peptide = 'mutant_seq'
        mut_aa_column = 'aa_mutant'
        mut_allele_column = 'mutant_best_alleles_netMHCpan'
        mut_pos_column = 'pep_mut_start'
        if mut_pos_column not in data.columns:
            return None, None, None, None
        else:
            return np.array(data[mut_peptide], dtype='str'), np.array(data[mut_aa_column], dtype='str'), \
                   np.array(data[mut_allele_column], dtype='str'), np.array(data[mut_pos_column], dtype='int32')


class SequenceLogo:

    def __init__(self, seq_logo_dir, allele, binding_threshold=10):

        self.allele = allele
        self.seq_logo_dir = seq_logo_dir
        self.freq_data = {}
        self.freq_data_sums = {}
        self.re = re.compile(r'.*class1_(\d+).*')
        self.threshold = binding_threshold

        self.read_AA_freq_data()

        return

    def read_AA_freq_data(self):

        allele = self.allele.replace(":", "").replace("*", "").replace("HLA-",  "")
        for freq_file in Path(self.seq_logo_dir).rglob(allele + "*"):

            length = int(self.re.match(str(freq_file)).group(1))
            data = pd.read_csv(freq_file, skiprows=6, sep="\t", header=None, names=['AA']+np.arange(length).tolist())
            data.set_index('AA', inplace=True)
            data = data.apply(np.log2)
            if length in self.freq_data:
                self.freq_data[length] = self.freq_data[length].add(data).div(2)
            else:
                self.freq_data[length] = data

        return

    def get_freq_data(self, length):
        return self.freq_data[length]

    def get_freq_data_sum(self, length):
        if length not in self.freq_data_sums:
            self.freq_data_sums[length] = self.freq_data[length].sum(axis=0)
        return self.freq_data_sums[length]

    def is_binding_pos(self, length, position):

        if not 0 < position <= length:
            return False

        if length not in self.freq_data_sums:
            self.freq_data_sums[length] = self.freq_data[length].sum(axis=0)

        return abs(self.freq_data_sums[length][position-1]) > self.threshold

    def get_aa_score(self, length, aa, position):

        if not 0 < position <= length:
            return 0

        data = self.freq_data[length]

        if aa not in data.index:
            return np.nan

        return data.loc[aa, position-1]
