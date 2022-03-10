import io

from Utils.Util_fct import *
from Features.NetMHC.NetMHCRunner import *
from Utils.DataManager import *


class NetMHCpanIRunner(NetMHCRunner):

    def __init__(self, exe_base_dir=None, score_key='Rank_EL', operator='<', max_rank_peptide=2, rank_threshold_strong=0.5,
                 rank_threshold_weak=2):

        NetMHCRunner.__init__(self, score_key, operator)
        self.netMHC_path = find_exe(exe_base_dir, exe='netMHCpan')
        self.max_rank_peptide = max_rank_peptide
        self.rank_threshold_strong = rank_threshold_strong
        self.rank_threshold_weak = rank_threshold_weak

        return

    def build_cmd(self, peptide_file, peptide_file_type, allele, length, show_cmd=False):
        input_type = 0 if peptide_file_type == 'fa' else 1
        cmd = '%s -BA -f %s -inptype %d -a %s' % (self.netMHC_path, peptide_file, input_type, allele)
        if show_cmd is True:
            print(cmd)

        return cmd

    def format_allele(self, allele):
        if not allele.startswith('HLA-'):
            allele = 'HLA-'+allele
        return allele.replace('*', '')

    def read_result(self, temp):
        """Read raw results from netMHCpan 4.1b output"""
        lines = temp.decode('ascii').split('\n')
        ignore = ['Pos', '#', 'Protein', 'PEPLIST', 'HLA', '-']
        lines = [l.lstrip() for l in lines if l != '' and list(filter(l.lstrip().startswith, ignore)) == []]
        cols = ['peptide_pos', 'allele', 'peptide', 'core', 'Of', 'Gp', 'Gl', 'Ip', 'Il', 'Icore',
                'Identity', 'Score_EL', 'Rank_EL', 'Score_BA', 'Rank_BA', 'ic50']
        data = [l.split()[0:len(cols)] for l in lines]

        res = pd.DataFrame(data=data, columns=cols)
        res = res.astype({'peptide_pos': 'int32', 'Score_EL': 'float', 'Rank_EL': 'float',
                          'Score_BA': 'float', 'Rank_BA': 'float', 'ic50': 'float'})

        if res.shape[0] > 0 and res.loc[res.index[0], 'Identity'] != 'PEPLIST':
            res = res.astype({'Identity': 'int32'})

        return res

    def process_results_mut(self, netmhc_results, mut_positions):
        """ For each peptide take the best max_rank_peptide alignments over the different alleles.
            Only peptides that contain the mutation site are considered. mut_positions are the
            mutation position in the long peptides (counted from 1)
        """
        peptide_ids = np.unique(netmhc_results['Identity'])

        res_list = []
        for pept_id in peptide_ids:
            mut_pos = mut_positions[pept_id]
            df = netmhc_results.loc[netmhc_results['Identity'] == pept_id]
            mut_in_peptide = [pos <= mut_pos < len(pept)+pos for (pept, pos) in zip(df['peptide'], df['peptide_pos'])]
            df = df.loc[mut_in_peptide]
            if df.shape[0] > 0:
                df = df.sort_values(by=[self.score_key], ascending=self.operator == '<', ignore_index=True)
                cnt_strong_binders = sum(df.Rank_EL < self.rank_threshold_strong)
                cnt_weak_binders = sum(df.Rank_EL < self.rank_threshold_weak)
                binding_alleles = np.unique(df.loc[df.Rank_EL < self.rank_threshold_weak, 'allele'])
                cnt_binding_alleles = len(binding_alleles)
                n = min(self.max_rank_peptide, df.shape[0])
                df = df.iloc[0:n, :]
                df['pos_in_peptide'] = np.subtract(mut_pos+1, df.loc[:, 'peptide_pos'])
                df['nr_strong_binders'] = np.full(n, cnt_strong_binders)
                df['nr_weak_binders'] = np.full(n, cnt_weak_binders)
                df['weak_binding_alleles'] = np.full(n, ",".join(binding_alleles))
                df['nr_weak_binding_alleles'] = np.full(n, cnt_binding_alleles)
                res_list.append(df)

        return pd.concat(res_list, ignore_index=True)

    def process_results_wt(self, wt_res_all, wt_peptides_short):
        """ For each wt peptide take the best alignment over the different alleles.
        """
        res_list = []
        for wt_peptide in wt_peptides_short:
            df = wt_res_all.loc[wt_res_all['peptide'] == wt_peptide]
            if df.shape[0] > 0:
                df = df.sort_values(by=[self.score_key], ascending=self.operator == '<', ignore_index=True)
                res_list.append(df.iloc[0, :])  # keep only best hit
            else:
                res_list.append(pd.DataFrame(columns=df.columns))  # fill with nan

        wt_res = pd.concat(res_list, axis=1, ignore_index=True).transpose()

        return wt_res

    def merge_with_data(self, data, netmhc_res):
        netmhc_res_list = []
        for i in np.arange(self.max_rank_peptide):
            columns = [c+'_'+str(i) for c in netmhc_res.columns]
            netmhc_res_list.append(pd.DataFrame(columns=columns, index=data.index))

        for pept_id in data.index.values:
            df = netmhc_res.loc[netmhc_res['mut_Identity'] == pept_id]
            for i in np.arange(self.max_rank_peptide):
                df_netmhc = netmhc_res_list[i]
                if i < df.shape[0]:
                    df_netmhc.loc[pept_id] = np.array(df.iloc[i, :])

        return pd.concat([data]+netmhc_res_list, axis=1)

    def add_features(self, patient, file_tag_input, file_tag_output, write_res=True):

        mgr = DataManager()

        if not file_tag_input:
            data = mgr.get_original_data(patient)
        else:
            data = mgr.get_processed_data(patient, file_tag_input)

        alleles = mgr.get_classI_allotypes(patient)

        indexes, peptide_id, mut_peptides, wt_peptides, mut_pos = NetMHCpanIRunner.get_peptide_info(data)

        mut_res_list = self.process_peptides(indexes, mut_peptides, 'fa', alleles, show_cmd=True)
        mut_res = self.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_pos)

        wt_peptides_short = NetMHCpanIRunner.get_wt_peptides(mut_res, wt_peptides)
        wt_res_list = self.process_peptides(mut_res.index, wt_peptides_short, 'pep', alleles, show_cmd=True)
        wt_res = self.process_results_wt(pd.concat(wt_res_list, ignore_index=True), wt_peptides_short)

        netmhc_res = NetMHCpanIRunner.combine_netmhc_results(wt_res, mut_res)
        data_netmhc = self.merge_with_data(data, netmhc_res)

        mgr.put_processed_data(data_netmhc, patient, "netmhc", "long")

        if write_res:
            parameters = Parameters()
            out_file = os.path.join(parameters.get_result_dir(), patient+"_long_"+file_tag_output+".txt")
            data_netmhc.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

        return data_netmhc

    @staticmethod
    def combine_netmhc_results(netmhc_wt_results, netmhc_mut_results):

        column_names = {}
        for s in netmhc_mut_results.columns:
            column_names[s] = 'mut_' + s
        netmhc_mut_results.rename(columns=column_names, inplace=True)

        column_names = {}
        cols = ['allele', 'peptide', 'Score_EL', 'Rank_EL', 'Score_BA', 'Rank_BA', 'ic50']
        netmhc_wt_results = netmhc_wt_results.loc[:, cols]
        for s in netmhc_wt_results.columns:
            column_names[s] = 'wt_' + s
        netmhc_wt_results.rename(columns=column_names, inplace=True)

        res = pd.merge(netmhc_mut_results, netmhc_wt_results, how='left', left_index=True, right_index=True)
        return res

    @staticmethod
    def get_wt_peptides(netmhc_res, wt_peptides_long):
        """Returns WT peptides for mutated peptides in netmhc_res. The order is the same as in netmhc_res."""

        wt_peptides = []
        for i in netmhc_res.index:
            pos = netmhc_res.loc[i, 'peptide_pos']-1
            pept_id = netmhc_res.loc[i, 'Identity']
            seq = netmhc_res.loc[i, 'peptide']
            wt_peptides.append(wt_peptides_long[pept_id][pos:(pos+len(seq))])

        return wt_peptides

    @staticmethod
    def get_peptide_info(data):
        return data.index, np.array(data['peptide_id'], dtype='str'), np.array(data['mutant_seq'], dtype='str'), \
               np.array(data['wt_seq'], dtype='str'), np.array(data['pep_mut_start'], dtype='int32')
