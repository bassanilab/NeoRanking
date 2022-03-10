from unittest import TestCase
from NetMHCpanIRunner import *
from Utils.DataManager import *


class TestNetMHCpanIRunner(TestCase):
    def test_process_peptides(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        indexes = [0]
        peptides = ['TPQDLNTMLNTVGGHQAAMQMLKETINEEA']
        alleles = ['HLA-A*03:01']
        mut_res_list = netMHCRunner.process_peptides(indexes, peptides, 'fa', alleles, show_cmd=True)

        self.assertEqual(1, len(mut_res_list))
        self.assertEqual(mut_res_list[0].shape, (86, 17))

    def test_process_results_mut1(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir(), max_rank_peptide=2)

        indexes = [0]
        peptides = ['HQAAMQMLK']
        alleles = ['A*03:01']
        mut_res_list = netMHCRunner.process_peptides(indexes, peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions=[3])

        self.assertEqual((2, 22), mut_res.shape)
        self.assertEqual('HQAAMQMLK', mut_res.loc[0, 'peptide'])

    def test_process_results_mut2(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir(), max_rank_peptide=3)

        indexes = [0]
        peptides = ['TPQDLNTMLNTVGGHQAAMQMLKETINEEA']
        alleles = ['HLA-A*03:01']
        mut_positions = [19]

        mut_res_list = netMHCRunner.process_peptides(indexes, peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True),
                                                   mut_positions=mut_positions)

        self.assertEqual((3, 22), mut_res.shape)
        self.assertEqual('HQAAMQMLK',  mut_res.loc[0, 'peptide'])

    def test_process_results_mut3(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir(), max_rank_peptide=3)

        indexes = [0]
        peptides = ['TPQDLNTMLNTVGGHQAAMQMLKETINEEA']
        alleles = ['HLA-A*03:01']
        mut_positions = [10]

        mut_res_list = netMHCRunner.process_peptides(indexes, peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True),
                                                   mut_positions=mut_positions)

        self.assertEqual((3, 22), mut_res.shape)
        self.assertEqual('TMLNTVGGH',  mut_res.loc[0, 'peptide'])
        self.assertEqual('MLNTVGGHQ',  mut_res.loc[1, 'peptide'])
        self.assertEqual('MLNTVGGH',  mut_res.loc[2, 'peptide'])

    def test_process_results_mut4(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        indexes = [0]
        peptides = ['TPQDLNTMLNTVGGHQAAMQMLKETINEEA']
        alleles = ['A*03:01']
        mut_res_list = netMHCRunner.process_peptides(indexes, peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions=[19])

        self.assertEqual((2, 22), mut_res.shape)
        self.assertEqual('HQAAMQMLK',  mut_res.iloc[0, mut_res.columns.get_loc('peptide')])

    def test_process_results_mut5(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        indexes = np.arange(10)
        peptides = ['LEQVKKQGVDLEQSRIEEEQQMLQ', 'KARRWDEKAVDKSKEKKERLTE', 'TPVLVLMAAVLTVTGPVPVARLHG',
                    'NVIKRLLCAQTFHTRIG', 'PSLWKMMLNIDVSATVFYKAQPV', 'KRSEVAMDFEPERVVAAPQR', 'ITVGQRIGSASFGTVYKGKW',
                    'KLRPPIPPMVILEPYVLSEL', 'NLGTDSDSSRQKSSRD', 'HAERKATMMYPESSTSEQEEAPL']
        mut_positions = [15, 12, 15, 9, 15, 6, 9, 5, 9, 7]
        alleles = ['A*01:01', 'A*03:01', 'B*27:05', 'B*57:01', 'C*01:02', 'C*06:02']

        mut_res_list = netMHCRunner.process_peptides(indexes, peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions)

        self.assertEqual((20, 22), mut_res.shape)

    def test_get_wt_peptides(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        indexes = np.arange(1)
        mut_peptides = ['LEQVKKQGVDLEQSRIEEEQQMLQ']
        wt_peptides =  ['LEQVKKQGVDLEQSRKEEEQQMLQ']

        mut_positions = [15]
        alleles = ['A*01:01']

        mut_res_list = netMHCRunner.process_peptides(indexes, mut_peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions)

        wt_peptides_short = NetMHCpanIRunner.get_wt_peptides(mut_res, wt_peptides)

        self.assertEqual(2, len(wt_peptides_short))
        self.assertEqual(['RKEEEQQML', 'GVDLEQSRK'], wt_peptides_short)

    def test_process_results_wt(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        indexes = np.arange(10)
        mut_peptides = ['LEQVKKQGVDLEQSRIEEEQQMLQ', 'KARRWDEKAVDKSKEKKERLTE', 'TPVLVLMAAVLTVTGPVPVARLHG',
                    'NVIKRLLCAQTFHTRIG', 'PSLWKMMLNIDVSATVFYKAQPV', 'KRSEVAMDFEPERVVAAPQR', 'ITVGQRIGSASFGTVYKGKW',
                    'KLRPPIPPMVILEPYVLSEL', 'NLGTDSDSSRQKSSRD', 'HAERKATMMYPESSTSEQEEAPL']

        wt_peptides = ['LEQVKKQGVDLEQSRKEEEQQMLQ', 'KARRWDEKAVDKLKEKKERLTE', 'TPVLVLMAAVLTVTGAVPVARLHG',
                       'NVIKRLLCARTFHTRIG', 'PSLWKMMLNIDVSATAFYKAQPV', 'KRSEVAKDFEPERVVAAPQR', 'ITVGQRIGSGSFGTVYKGKW',
                       'KLRPPTPPMVILEPYVLSEL', 'NLGTDSDSSPQKSSRD', 'HAERKATRMYPESSTSEQEEAPL']

        mut_positions = [15, 12, 15, 9, 15, 6, 9, 5, 9, 7]
        alleles = ['A*01:01', 'A*03:01', 'B*27:05', 'B*57:01', 'C*01:02', 'C*06:02']

        mut_res_list = netMHCRunner.process_peptides(indexes, mut_peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions)

        wt_peptides_short = NetMHCpanIRunner.get_wt_peptides(mut_res, wt_peptides)
        wt_res_list = netMHCRunner.process_peptides(mut_res.index, wt_peptides_short, 'pep', alleles, show_cmd=True)

        wt_res = netMHCRunner.process_results_wt(pd.concat(wt_res_list, ignore_index=True), wt_peptides_short)

        self.assertEqual(mut_res.shape[0], wt_res.shape[0])
        self.assertTrue((mut_res.index.values == wt_res.index.values).all())

    def test_merge_with_data(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        indexes = np.arange(10)
        mut_peptides = ['LEQVKKQGVDLEQSRIEEEQQMLQ', 'KARRWDEKAVDKSKEKKERLTE', 'TPVLVLMAAVLTVTGPVPVARLHG',
                    'NVIKRLLCAQTFHTRIG', 'PSLWKMMLNIDVSATVFYKAQPV', 'KRSEVAMDFEPERVVAAPQR', 'ITVGQRIGSASFGTVYKGKW',
                    'KLRPPIPPMVILEPYVLSEL', 'NLGTDSDSSRQKSSRD', 'HAERKATMMYPESSTSEQEEAPL']

        wt_peptides = ['LEQVKKQGVDLEQSRKEEEQQMLQ', 'KARRWDEKAVDKLKEKKERLTE', 'TPVLVLMAAVLTVTGAVPVARLHG',
                       'NVIKRLLCARTFHTRIG', 'PSLWKMMLNIDVSATAFYKAQPV', 'KRSEVAKDFEPERVVAAPQR', 'ITVGQRIGSGSFGTVYKGKW',
                       'KLRPPTPPMVILEPYVLSEL', 'NLGTDSDSSPQKSSRD', 'HAERKATRMYPESSTSEQEEAPL']

        mut_positions = [15, 12, 15, 9, 15, 6, 9, 5, 9, 7]
        alleles = ['A*01:01', 'A*03:01', 'B*27:05', 'B*57:01', 'C*01:02', 'C*06:02']

        mut_res_list = netMHCRunner.process_peptides(indexes, mut_peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions)

        wt_peptides_short = NetMHCpanIRunner.get_wt_peptides(mut_res, wt_peptides)
        wt_res_list = netMHCRunner.process_peptides(mut_res.index, wt_peptides_short, 'pos', alleles, show_cmd=True)

        wt_res = netMHCRunner.process_results_wt(pd.concat(wt_res_list, ignore_index=True), wt_peptides_short)

        netmhc_res = NetMHCpanIRunner.combine_netmhc_results(wt_res, mut_res)
        self.assertEqual(mut_res.shape[0], netmhc_res.shape[0])

    def test_combine_netmhc_results(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir())

        patient = '0YM1'
        mgr = DataManager()
        data = mgr.get_original_data(patient)
        alleles = mgr.get_classI_allotypes(patient)

        data = data.head(10)

        indexes, peptide_id, mut_peptides, wt_peptides, mut_positions = NetMHCpanIRunner.get_peptide_info(data)

        mut_res_list = netMHCRunner.process_peptides(indexes, mut_peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions)

        wt_peptides_short = NetMHCpanIRunner.get_wt_peptides(mut_res, wt_peptides)
        wt_res_list = netMHCRunner.process_peptides(mut_res.index, wt_peptides_short, 'pep', alleles, show_cmd=True)

        wt_res = netMHCRunner.process_results_wt(pd.concat(wt_res_list, ignore_index=True), wt_peptides_short)

        netmhc_res = NetMHCpanIRunner.combine_netmhc_results(wt_res, mut_res)
        data_netmhc = netMHCRunner.merge_with_data(data, netmhc_res)

        self.assertEqual(data.shape[0], data_netmhc.shape[0])

        netMHCpan_ranks = [int(c[c.rfind('_')+1:]) for c in data_netmhc.columns if 'mut_peptide_pos_' in c]
        for r in netMHCpan_ranks:
            mut_pos = []
            seq_pairs = zip(data_netmhc['mutant_seq'], data_netmhc['mut_peptide_'+str(r)])
            for seq_long, seq_short in seq_pairs:
                self.assertTrue(seq_short in seq_long)
                mut_pos.append(seq_long.find(seq_short))

            wt_pos = []
            seq_pairs = zip(data_netmhc['wt_seq'], data_netmhc['wt_peptide_'+str(r)])
            for seq_long, seq_short in seq_pairs:
                self.assertTrue(seq_short in seq_long)
                wt_pos.append(seq_long.find(seq_short))

            self.assertTrue(np.max(np.subtract(mut_pos, wt_pos)) < 3)

    def test_combine_netmhc_results2(self):
        parameters = Parameters()
        netMHCRunner = NetMHCpanIRunner(exe_base_dir=parameters.get_exe_dir(), max_rank_peptide=10)

        patient = '4274'
        mgr = DataManager()
        data = mgr.get_original_data(patient)
        alleles = mgr.get_classI_allotypes(patient)

        data = data.head(10)

        indexes, peptide_id, mut_peptides, wt_peptides, mut_positions = NetMHCpanIRunner.get_peptide_info(data)

        mut_res_list = netMHCRunner.process_peptides(indexes, mut_peptides, 'fa', alleles, show_cmd=True)
        mut_res = netMHCRunner.process_results_mut(pd.concat(mut_res_list, ignore_index=True), mut_positions)

        wt_peptides_short = NetMHCpanIRunner.get_wt_peptides(mut_res, wt_peptides)
        wt_res_list = netMHCRunner.process_peptides(mut_res.index, wt_peptides_short, 'pep', alleles, show_cmd=True)

        wt_res = netMHCRunner.process_results_wt(pd.concat(wt_res_list, ignore_index=True), wt_peptides_short)

        netmhc_res = NetMHCpanIRunner.combine_netmhc_results(wt_res, mut_res)
        data_netmhc = netMHCRunner.merge_with_data(data, netmhc_res)

        self.assertEqual(data.shape[0], data_netmhc.shape[0])

        netMHCpan_ranks = [int(c[c.rfind('_')+1:]) for c in data_netmhc.columns if 'mut_peptide_pos_' in c]
        for r in netMHCpan_ranks:
            mut_pos = []
            seq_pairs = zip(data_netmhc['mutant_seq'], data_netmhc['mut_peptide_'+str(r)])
            for seq_long, seq_short in seq_pairs:
                self.assertTrue(seq_short in seq_long)
                mut_pos.append(seq_long.find(seq_short))

            # check that mutant_seq contain mutation
            seq_pairs = zip(data_netmhc['wt_seq'], data_netmhc['mut_peptide_'+str(r)])
            for wt_seq_long, mut_seq_short in seq_pairs:
                self.assertFalse(mut_seq_short in wt_seq_long)

            wt_pos = []
            seq_pairs = zip(data_netmhc['wt_seq'], data_netmhc['wt_peptide_'+str(r)])
            for seq_long, seq_short in seq_pairs:
                self.assertTrue(seq_short in seq_long)
                wt_pos.append(seq_long.find(seq_short))

            self.assertTrue(np.max(np.subtract(mut_pos, wt_pos)) < 3)

            seq_triplet = zip(data_netmhc['mut_peptide_'+str(r)], data_netmhc['wt_peptide_'+str(r)],
                              data_netmhc['mut_pos_in_peptide_'+str(r)])
            for mut_seq, wt_seq, mut_pos in seq_triplet:
                if len(mut_seq) == len(wt_seq):
                    self.assertTrue(mut_seq[mut_pos-1] != wt_seq[mut_pos-1])

