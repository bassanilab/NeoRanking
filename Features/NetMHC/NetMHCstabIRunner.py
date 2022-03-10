import io
from Utils.Util_fct import *
from Features.NetMHC.NetMHCRunner import *
from Utils.DataManager import *


class NetMHCstabIRunner(NetMHCRunner):

    def __init__(self, exe_base_dir, score_key='Rank_Stab', operator='<'):

        NetMHCRunner.__init__(self, score_key, operator)
        self.netMHC_path = find_exe(exe_base_dir, exe='netMHCstabpan')
        self.netmhc_result_skip_rows = 49
        return

    def build_cmd(self, peptide_file, peptide_file_type, allele, length, show_cmd=False):
        input_type = 0 if peptide_file_type == 'fa' else 1
        cmd = '%s -f %s -inptype %d -a %s' % (self.netMHC_path, peptide_file, input_type, allele)
        if show_cmd:
            print(cmd)

        return cmd

    def read_result(self, temp):
        """Read raw results from netMHCstabpan 1.0 output"""

        lines = temp.decode('ascii').split('\n')
        ignore = ['pos', '#', 'Protein', 'PEPLIST', 'HLA', '-']
        lines = [l.lstrip() for l in lines if l != '' and list(filter(l.lstrip().startswith, ignore)) == []]
        cols = ['peptide_pos', 'allele', 'peptide', 'Identity', 'Stab_Score', 'Thalf', 'Rank_Stab']
        data = [l.split()[0:len(cols)] for l in lines]

        res = pd.DataFrame(data=data, columns=cols)
        res = res.astype({'peptide_pos': 'int32', 'Stab_Score': 'float', 'Thalf': 'float', 'Rank_Stab': 'float'})

        if res.shape[0] > 0 and res.loc[res.index[0], 'Identity'] != 'PEPLIST':
            res = res.astype({'Identity': 'int32'})

        return res

    def process_results(self, res_all, peptides):
        """ For each peptide take the best score over the different alleles. """
        res_list = []
        for peptide in peptides:
            df = res_all.loc[res_all['peptide'] == peptide]
            if df.shape[0] > 0:
                df = df.sort_values(by=[self.score_key], ascending=self.operator == '<', ignore_index=True)
                res_list.append(df.iloc[0, :])  # keep only best hit
            else:
                res_list.append(pd.Series({c: '' for c in df.columns}))  # add empty row

        res = pd.concat(res_list, axis=1, ignore_index=True).transpose()

        return res

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):

        mgr = DataManager()
        parameters = Parameters()

        data = mgr.get_processed_data(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            if 'mut_peptide_0' not in data.columns:
                print('Data file %s does not contain netMHC predictions. Unable to add stability prediction.')
                return None
            data = self.add_features_long(data)
        else:
            if 'mutant_best_alleles_netMHCpan' not in data.columns:
                print('Data file %s does not contain netMHC predictions. Unable to add stability prediction.')
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
            idx, mut_peptides, mut_alleles, wt_peptides, wt_alleles = NetMHCstabIRunner.get_peptide_info_long(data, i)

            pept_dict = NetMHCRunner.split_length_allele(mut_peptides, mut_alleles)
            res_list = []
            for (a, l) in pept_dict.keys():
                peptides = pept_dict[(a, l)]
                res_list = res_list + self.process_peptides(
                     np.arange(len(peptides)), peptides, 'pep', [a], l, show_cmd=True)
            mut_res = self.process_results(pd.concat(res_list, ignore_index=True), mut_peptides)

            pept_dict = NetMHCRunner.split_length_allele(wt_peptides, wt_alleles)
            res_list = []
            for (a, l) in pept_dict.keys():
                peptides = pept_dict[(a, l)]
                res_list = res_list + self.process_peptides(
                    np.arange(len(peptides)), peptides, 'pep', [a], l, show_cmd=True)

            wt_res = self.process_results(pd.concat(res_list, ignore_index=True), wt_peptides)

            res = self.combine_netmhc_results(wt_res, mut_res)
            data = NetMHCstabIRunner.merge_with_data_long(data, res, i)

        return data

    def add_features_short(self, data):
        idx, mut_peptides, mut_alleles, wt_peptides, wt_alleles, mut_uniq_pepts, wt_uniq_pepts = \
            NetMHCstabIRunner.get_peptide_info_short(data)

        pept_dict = NetMHCRunner.split_length_allele(mut_peptides, mut_alleles)
        res_list = []
        for (a, l) in pept_dict.keys():
            peptides = pept_dict[(a, l)]
            res_list = res_list + self.process_peptides(
                 np.arange(len(peptides)), peptides, 'pep', [a], l, show_cmd=True)
        mut_res = self.process_results(pd.concat(res_list, ignore_index=True), mut_uniq_pepts)

        pept_dict = NetMHCRunner.split_length_allele(wt_peptides, wt_alleles)
        res_list = []
        for (a, l) in pept_dict.keys():
            peptides = pept_dict[(a, l)]
            res_list = res_list + self.process_peptides(
                np.arange(len(peptides)), peptides, 'pep', [a], l, show_cmd=True)

        wt_res = self.process_results(pd.concat(res_list, ignore_index=True), wt_uniq_pepts)

        res = self.combine_netmhc_results(wt_res, mut_res)
        data = NetMHCstabIRunner.merge_with_data_short(data, res)

        return data

    @staticmethod
    def combine_netmhc_results(netmhc_wt_results, netmhc_mut_results):

        cols = ['Stab_Score', 'Thalf', 'Rank_Stab']
        # cols = ['allele','peptide','Stab_Score', 'Thalf', 'Rank_Stab'] # for debugging
        netmhc_mut_results = netmhc_mut_results.loc[:, cols]
        netmhc_wt_results = netmhc_wt_results.loc[:, cols]

        column_names = {}
        for s in netmhc_mut_results.columns:
            column_names[s] = 'mut_' + s
        netmhc_mut_results.rename(columns=column_names, inplace=True)

        column_names = {}
        for s in netmhc_wt_results.columns:
            column_names[s] = 'wt_' + s
        netmhc_wt_results.rename(columns=column_names, inplace=True)

        res = pd.merge(netmhc_mut_results, netmhc_wt_results, how='left', left_index=True, right_index=True)

        return res

    @staticmethod
    def merge_with_data_long(data, netmhc_res, index):
        column_names = {}
        for c in netmhc_res.columns:
            column_names[c] = c + '_' + str(index)
        netmhc_res.rename(columns=column_names, inplace=True)

        return pd.merge(data, netmhc_res, how='left', left_index=True, right_index=True)

    @staticmethod
    def merge_with_data_short(data, netmhc_res):
        return pd.merge(data, netmhc_res, how='left', left_index=True, right_index=True)

    @staticmethod
    def get_peptide_info_long(data, index):
        mut_peptide_column = 'mut_peptide_'+str(index)
        wt_peptide_column = 'wt_peptide_'+str(index)
        mut_allele_column = 'mut_allele_'+str(index)
        wt_allele_column = 'wt_allele_'+str(index)
        if mut_peptide_column not in data.columns:
            return None, None, None, None, None
        else:
            return np.arange(data.shape[0]), \
                   np.array(data[mut_peptide_column], dtype='str'), \
                   np.array(data[mut_allele_column], dtype='str'), \
                   np.array(data[wt_peptide_column], dtype='str'), \
                   np.array(data[wt_allele_column], dtype='str')

    @staticmethod
    def get_peptide_info_short(data):
        mut_peptide_column = 'mutant_seq'
        wt_peptide_column = 'wt_seq'
        mut_allele_column = 'mutant_best_alleles_netMHCpan'
        wt_allele_column = 'wt_best_alleles_netMHCpan'
        if mut_peptide_column not in data.columns or mut_allele_column not in data.columns:
            return None, None, None, None, None, None, None
        else:
            mut_alleles = []
            mut_peptides = []
            for allele_str, p in zip(data[mut_allele_column], data[mut_peptide_column]):
                alleles = str(allele_str).split(',')
                for a in alleles:
                    mut_alleles.append(a)
                    mut_peptides.append(p)
            wt_alleles = []
            wt_peptides = []
            for allele_str, p in zip(data[wt_allele_column], data[wt_peptide_column]):
                alleles = str(allele_str).split(',')
                for a in alleles:
                    wt_alleles.append(a)
                    wt_peptides.append(p)
            mut_alleles = list(map(lambda a: NetMHCstabIRunner.format_allele_short(a), mut_alleles))
            wt_alleles = list(map(lambda a: NetMHCstabIRunner.format_allele_short(a), wt_alleles))
            return np.arange(data.shape[0]), \
                   np.array(mut_peptides, dtype='str'), \
                   np.array(mut_alleles, dtype='str'), \
                   np.array(wt_peptides, dtype='str'), \
                   np.array(wt_alleles, dtype='str'), \
                   np.array(data[mut_peptide_column], dtype='str'), \
                   np.array(data[wt_peptide_column], dtype='str')

    @staticmethod
    def format_allele_short(allele):
        return 'HLA-'+allele[0:3]+':'+allele[3:]
