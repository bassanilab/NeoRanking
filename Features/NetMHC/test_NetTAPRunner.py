from unittest import TestCase

from Features.NetMHC.NetTAPRunner import *
from Utils.DataManager import *


class TestNetTAPRunner(TestCase):

    def test_constructor(self):
        nettapRunner = NetTAPRunner(exe_path=os.path.join(Parameters().get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))

        self.assertTrue(os.path.isfile(nettapRunner.netTAP_path))

    def test_build_command(self):
        nettapRunner = NetTAPRunner(exe_path=os.path.join(Parameters().get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))

        cmd = nettapRunner.build_cmd('./data/peptide_file.fa', 'fa', '', 9, show_cmd=False)
        self.assertEqual(str(nettapRunner.netTAP_path)+' -l 9 ./data/peptide_file.fa', cmd)

    def test_process_peptides(self):
        nettapRunner = NetTAPRunner(exe_path=os.path.join(Parameters().get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))

        idx = list(range(5))
        peptides = ['PEPTIDEER', 'THEANTLER', 'GENEVALAC', 'TAMARANAN', 'MARMILLER']
        res_list = nettapRunner.process_peptides(idx, peptides, 'fa', ['a'], length=9, show_cmd=True)

        self.assertEqual(1, len(res_list))
        self.assertEqual((5, 3), res_list[0].shape)
        self.assertTrue(all(p in peptides for p in res_list[0].peptide))

    def test_split_length(self):
        data = DataManager().get_processed_data('0YM1', 'rt_netmhc')
        idx_0, mut_peptides_0, wt_peptides_0 = NetTAPRunner.get_peptide_info_long(data, 0)

        pept_dict = NetMHCRunner.split_length(mut_peptides_0)

        self.assertEqual(4, len(pept_dict))
        for l in pept_dict.keys():
            for p in pept_dict[l]:
                self.assertEqual(l, len(p))

    def test_add_features_long(self):
        nettapRunner = NetTAPRunner(exe_path=os.path.join(Parameters().get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))

        tap_data = nettapRunner.add_features('0YM1', 'rt_netmhc_stab_chop', 'tap', 'long', write_res=False)

        self.assertEqual(tap_data.shape[0], DataManager().get_processed_data('0YM1', 'rt_netmhc', 'long').shape[0])

    def test_add_features_short(self):
        nettapRunner = NetTAPRunner(exe_path=os.path.join(Parameters().get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))

        tap_data = nettapRunner.add_features('0YM1', 'rt', 'tap', 'short', write_res=False)

        self.assertEqual(tap_data.shape[0], DataManager().get_processed_data('0YM1', 'rt', 'short').shape[0])

    def test_add_features_short(self):

        nettapRunner = NetTAPRunner(exe_path=os.path.join(Parameters().get_exe_dir(), 'netCTLpan-1.1', 'Linux_x86_64'))

        mut_peptides = np.array(['LDRILQAGL', 'DRILQAGLD', 'RILQAGLDV', 'ILQAGLDVE', 'LQAGLDVER'], dtype=str)
        pept_dict = NetMHCRunner.split_length(mut_peptides)
        res_list = []
        for length in pept_dict.keys():
            peptides = pept_dict[length]
            res_list = res_list + nettapRunner.process_peptides(
                np.arange(len(peptides)), peptides, 'fa', ['a'], length, show_cmd=True)
        mut_res = nettapRunner.process_results(pd.concat(res_list, ignore_index=True), mut_peptides)

        tap_scores = np.array([0.84000, -1.68600, 0.68800, -1.65400, 1.78400], dtype=float)
        for i in range(len(mut_peptides)):
            self.assertEqual(mut_res.loc[mut_res.index[i], 'peptide'], mut_peptides[i])
            self.assertEqual(mut_res.loc[mut_res.index[i], 'TAP_score'], tap_scores[i])
