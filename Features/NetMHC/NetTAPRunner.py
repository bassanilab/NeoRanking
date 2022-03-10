import io

from Utils.Util_fct import *
from Features.NetMHC.NetMHCRunner import *
from Utils.DataManager import *


class NetTAPRunner(NetMHCRunner):

    def __init__(self, exe_path, score_key='TAP_score', operator='<'):

        NetMHCRunner.__init__(self, score_key, operator)
        self.netTAP_path = find_exe(exe_path, exe='tapmat_pred_fsa')
        self.nettap_result_skip_rows = 12
        return

    def build_cmd(self, peptide_file, peptide_file_type, allele, length, show_cmd=False):
        # tapmat_pred_fsa -l length -mat tap_matrix_filename
        cmd = '%s -l %d %s' % (self.netTAP_path, length, peptide_file)
        if show_cmd:
            print(cmd)

        return cmd

    def read_result(self, temp):
        """Read raw results from nettap 1.1 output"""
        ignore = ['#', '']

        cols = ['id', 'peptide', 'TAP_score']
        res = pd.read_csv(io.BytesIO(temp), comment='#', names=cols, sep=r'\s+',
                          error_bad_lines=False, skiprows=self.nettap_result_skip_rows).dropna(subset=['TAP_score'])
        res = res[~res.id.isin(ignore)]
        res = res.astype({'id': 'int32', 'peptide': 'str','TAP_score': 'float'})
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

        netMHCpan_ranks = [int(c[c.rfind('_')+1:]) for c in data.columns if 'mut_peptide_pos_' in c]
        for i in netMHCpan_ranks:
            idx, mut_peptides, wt_peptides = NetTAPRunner.get_peptide_info_long(data, i)

            pept_dict = NetMHCRunner.split_length(mut_peptides)
            res_list = []
            for length in pept_dict.keys():
                peptides = pept_dict[length]
                res_list = res_list + self.process_peptides(
                    np.arange(len(peptides)), peptides, 'fa', ['a'], length, show_cmd=True)
            mut_res = self.process_results(pd.concat(res_list, ignore_index=True), mut_peptides)

            pept_dict = NetMHCRunner.split_length(wt_peptides)
            res_list = []
            for length in pept_dict.keys():
                peptides = pept_dict[length]
                res_list = res_list + self.process_peptides(
                    np.arange(len(peptides)), peptides, 'fa', ['a'], length, show_cmd=True)
            wt_res = self.process_results(pd.concat(res_list, ignore_index=True), wt_peptides)

            res = self.combine_results(wt_res, mut_res)
            data = NetTAPRunner.merge_with_data_long(data, res, i)

        return data

    def add_features_short(self, data):

        idx, mut_peptides, wt_peptides = NetTAPRunner.get_peptide_info_short(data)

        pept_dict = NetMHCRunner.split_length(mut_peptides)
        res_list = []
        for length in pept_dict.keys():
            peptides = pept_dict[length]
            res_list = res_list + self.process_peptides(
                np.arange(len(peptides)), peptides, 'fa', ['a'], length, show_cmd=True)
        mut_res = self.process_results(pd.concat(res_list, ignore_index=True), mut_peptides)

        data = NetTAPRunner.merge_with_data_short(data, mut_res)

        return data

    @staticmethod
    def process_results(res_all, peptides):
        """ For each peptide take the best score over the different alleles. """
        res_list = []
        for peptide in peptides:
            df = res_all.loc[res_all['peptide'] == peptide]
            if df.shape[0] > 0:
                df = df.sort_values(by=['TAP_score'], ascending=False, ignore_index=True)
                res_list.append(df.iloc[0, :])  # keep only best hit
            else:
                res_list.append(pd.Series({'id': 0, 'peptide': 'NA', 'TAP_score': 'NA'}))  # add empty row

        res = pd.concat(res_list, axis=1, ignore_index=True).transpose()

        return res

    @staticmethod
    def combine_results(wt_results, mut_results):
        cols = ['TAP_score']
        mut_results = mut_results.loc[:, cols]
        column_names = {}
        for s in mut_results.columns:
            column_names[s] = 'mut_' + s
        mut_results.rename(columns=column_names, inplace=True)

        column_names = {}
        wt_results = wt_results.loc[:, cols]
        for s in wt_results.columns:
            column_names[s] = 'wt_' + s
        wt_results.rename(columns=column_names, inplace=True)

        res = pd.merge(mut_results, wt_results, how='left', left_index=True, right_index=True)
        return res

    @staticmethod
    def merge_with_data_long(data, netchop_res, index):
        column_names = {}
        for c in netchop_res.columns:
            column_names[c] = c + '_' + str(index)
        netchop_res.rename(columns=column_names, inplace=True)

        return pd.merge(data, netchop_res, how='left', left_index=True, right_index=True)

    @staticmethod
    def merge_with_data_short(data, netchop_res):
        return pd.merge(data, netchop_res, how='left', left_index=True, right_index=True)

    @staticmethod
    def get_peptide_info_long(data, index):
        mut_peptide_column = 'mut_peptide_'+str(index)
        wt_peptide_column = 'wt_peptide_'+str(index)
        if mut_peptide_column not in data.columns:
            return None, None, None
        else:
            return data.index, np.array(data[mut_peptide_column], dtype='str'), \
                   np.array(data[wt_peptide_column], dtype='str')

    @staticmethod
    def get_peptide_info_short(data):
        mut_peptide_column = 'mutant_seq'
        wt_peptide_column = 'wt_seq'
        if mut_peptide_column not in data.columns:
            return None, None, None
        else:
            return data.index, np.array(data[mut_peptide_column], dtype='str'), \
                   np.array(data[wt_peptide_column], dtype='str')

