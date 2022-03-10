import io

import pandas as pd

from Utils.Util_fct import *
from Features.NetMHC.NetMHCRunner import *
from Utils.DataManager import *
from Bio import SeqIO


class NetChopRunner(NetMHCRunner):

    def __init__(self, exe_path, score_key='chop_Score', operator='<', model='Cterm-3.0'):

        NetMHCRunner.__init__(self, score_key, operator)
        self.netMHC_path = find_exe(exe_path, exe='netChop')
        self.netmhc_result_skip_rows = 18
        self.model = model
        self.fasta_sequences = {}
        self.read_fasta_seqs()
        return

    def read_fasta_seqs(self):
        with open(Parameters().get_protein_seq_file('37')) as file:
            fasta_sequences = SeqIO.parse(file, 'fasta')

            for fasta_item in fasta_sequences:
                header, sequence = fasta_item.id, str(fasta_item.seq)
                ensembl_prot_name = header.split('|')[0]
                self.fasta_sequences[ensembl_prot_name] = sequence

    def build_cmd(self, peptide_file, peptide_file_type, allele, length, show_cmd=False):
        # netchop [-d dir] [-s] [-t #] [-tdir file] [-v #] [-verbose]
        m = 1 if "20S" in self.model else 0
        cmd = '%s %s -v %d' % (self.netMHC_path, peptide_file, m)
        if show_cmd:
            print(cmd)

        return cmd

    def read_result(self, temp):
        """Read raw results from netMHCstabpan 1.0 output"""
        # lines = temp.decode('ascii').split('\n')
        ignore = ['pos', '-', 'Number', '']

        cols = ['peptide_pos', 'aa', 'cleavage', 'chop_Score', 'Ident']
        res = pd.read_csv(io.BytesIO(temp), comment='#', names=cols, sep=r'\s+',
                          error_bad_lines=False, skiprows=self.netmhc_result_skip_rows).dropna(subset=['aa'])
        res = res[~res.peptide_pos.isin(ignore)]
        res = res.astype({'peptide_pos': 'int32', 'chop_Score': 'float', 'Ident': 'int32'})
        return res

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):

        mgr = DataManager()
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            if 'mut_peptide_0' not in data.columns:
                print('Data file %s does not contain netMHC predictions. Unable to add netChop prediction.')
                return None
            data = self.add_features_long(data)
        else:
            if 'mutant_best_alleles_netMHCpan' not in data.columns:
                print('Data file %s does not contain netMHC predictions. Unable to add netChop prediction.')
                return None
            data = self.add_features_short(data)

        mgr.put_processed_data(data, patient, file_tag_output, peptide_type)

        if write_res:
            out_file = os.path.join(parameters.get_result_dir(), patient+"_"+peptide_type+"_"+file_tag_output+".txt")
            data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

        return data

    def add_features_long(self, data):

        idx, mut_peptides, mut_pos, wt_peptides, protein_ids = NetChopRunner.get_peptide_info_long(data)

        res_list = self.process_peptides(idx, mut_peptides, 'fa', alleles=['a'], show_cmd=True)
        mut_score_dict = self.process_results(pd.concat(res_list, ignore_index=True))

        res = self.combine_netchop_results_long(idx, mut_peptides, mut_pos, mut_score_dict)
        data = NetChopRunner.merge_with_data(data, res)

        return data

    def add_features_short(self, data):

        idx, starts, ends, ext_seq = self.get_peptide_info_short(data)

        res_list = self.process_peptides(idx, ext_seq, 'fa', alleles=['a'], show_cmd=True)
        mut_score_dict = self.process_results(pd.concat(res_list, ignore_index=True))

        res = self.combine_netchop_results_short(idx, mut_score_dict, starts, ends)
        data = NetChopRunner.merge_with_data(data, res)

        return data

    @staticmethod
    def process_results(res_all):
        """ For each peptide return scores for all positions in long peptide. """
        peptide_ids = np.unique(res_all['Ident'])

        res_dict = {}
        for pept_id in peptide_ids:
            df = res_all.loc[res_all['Ident'] == pept_id]
            if df.shape[0] > 0:
                res_dict[pept_id] = np.array(df.loc[:, 'chop_Score'])

        return res_dict

    @staticmethod
    def combine_netchop_results_long(index, mut_peptides, positions, mut_score_dict):
        mut_netchop_scores = []
        for i in range(len(mut_peptides)):
            scores = mut_score_dict[index[i]]
            chop_scores = []
            for j in np.arange(start=max(0, positions[i]-9), stop=min(positions[i]+1, len(mut_peptides[i])-10)):
                nt_chop_score = 0 if j <= 2 else max(scores[max(0, j-8):max(0, j-2)])
                chop_scores.append(3*scores[j+9] + 1.5*nt_chop_score - max(scores[j:j+9]))

            mut_netchop_scores.append(max(chop_scores) if len(chop_scores) > 0 else 0)

        return pd.DataFrame({'mut_netchop_score': mut_netchop_scores}, index=index)

    @staticmethod
    def combine_netchop_results_short(idx, mut_score_dict, starts, ends):

        mut_scores_C = []
        mut_scores_N = []
        mut_scores_I = []
        for i in range(len(starts)):
            if starts[i] < 0:
                mut_scores_C.append(0)
                mut_scores_N.append(0)
                mut_scores_I.append(0)
                continue

            scores = mut_score_dict[idx[i]]
            mut_scores_C.append(scores[ends[i]-1])
            if starts[i] > 3:
                mut_scores_N.append(max(scores[max(0, starts[i]-8):max(0, starts[i]-2)]))
            else:
                mut_scores_N.append(0)
            mut_scores_I.append(max(scores[starts[i]:ends[i]]))

        return pd.DataFrame({'mut_netchop_score_ct': mut_scores_C, 'mut_netchop_score_nt': mut_scores_N,
                             'mut_netchop_score_int': mut_scores_I}, index=idx)

    def get_extended_seq(self, prot_id, wt_seq, mut_seq, protein_coord):

        if pd.isna(wt_seq):
            return -1, -1, ""

        prot_seq = self.fasta_sequences[prot_id]
        for m in re.finditer(wt_seq, prot_seq):
            if m.start() <= protein_coord < m.end():
                start_idx = max(0, m.start()-10)
                end_idx = min(len(prot_seq), m.end()+10)
                ext_seq = prot_seq[start_idx:end_idx]
                ext_seq.replace(wt_seq, mut_seq)
                start_idx = m.start() - start_idx
                end_idx = start_idx + len(wt_seq)
                return start_idx, end_idx, ext_seq

        return -1, -1, ""

    @staticmethod
    def merge_with_data(data, netchop_res):
        return pd.merge(data, netchop_res, how='left', left_index=True, right_index=True)

    @staticmethod
    def get_peptide_info_long(data):
        mut_peptide_column = 'mutant_seq'
        mut_pos_column = 'pep_mut_start'
        wt_peptide_column = 'wt_seq'
        db_column = 'database_entry'
        return data.index, np.array(data[mut_peptide_column], dtype='str'), \
               np.array(data[mut_pos_column], dtype='int'), np.array(data[wt_peptide_column], dtype='str'), \
               np.array(data[db_column], dtype='str')

    def get_peptide_info_short(self, data):
        ext_seq_info = \
            data.apply(lambda r:
                       self.get_extended_seq(r['database_entry'], r['wt_seq'], r['mutant_seq'], r['protein_coord']),
                       axis=1)

        starts = []
        ends = []
        ext_seqs = []
        for (start, end, ext_seq) in ext_seq_info:
            starts.append(start)
            ends.append(end)
            ext_seqs.append(ext_seq)

        return data.index, np.array(starts, dtype='int'), np.array(ends, dtype='int'), np.array(ext_seqs, dtype='str')
