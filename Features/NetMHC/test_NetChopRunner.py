from unittest import TestCase

import numpy as np

from Features.NetMHC.NetChopRunner import *
from Utils.DataManager import *


class TestNetChopRunner(TestCase):

    def test_constructor(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        self.assertTrue(os.path.isfile(netChopRunner.netMHC_path))

    def test_fasta(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        self.assertEqual(106300, len(netChopRunner.fasta_sequences))
        self.assertTrue('ENSP00000266383.5' in netChopRunner.fasta_sequences)
        self.assertTrue('LDRILQAGLDVERLSLRNFFHHFHS' in netChopRunner.fasta_sequences['ENSP00000266383.5'])

    def test_get_peptide_info_long(self):
        mgr = DataManager()
        data = mgr.get_processed_data('0YM1', 'rt_netmhc_stab', 'long')

        idx, mut_peptides, mut_pos, wt_peptides, protein_ids = NetChopRunner.get_peptide_info_long(data)

        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')
        self.assertEqual(len(mut_peptides), len(wt_peptides))
        for i in range(len(wt_peptides)):
            print(str(i)+":"+wt_peptides[i]+", "+protein_ids[i])
            self.assertTrue(protein_ids[i] in netChopRunner.fasta_sequences)
            self.assertTrue(wt_peptides[i] in netChopRunner.fasta_sequences[protein_ids[i]])

    def test_process_peptides(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        mgr = DataManager()
        data = mgr.get_processed_data('13LN', 'rt_netmhc_stab')

        idx, mut_peptides, mut_pos, wt_peptides, protein_ids = NetChopRunner.get_peptide_info_long(data)

        res_list = netChopRunner.process_peptides(idx, mut_peptides, 'fa', ['a'], show_cmd=True)

        self.assertEqual(1, len(res_list))
        self.assertTrue(all(id in idx for id in np.unique(res_list[0].Ident)))

    def test_process_results(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        mgr = DataManager()
        data = mgr.get_processed_data('13LN', 'rt_netmhc_stab')

        idx, mut_peptides, mut_pos, wt_peptides, protein_ids = NetChopRunner.get_peptide_info_long(data)

        res_list = netChopRunner.process_peptides(idx, mut_peptides, 'fa', ['a'], show_cmd=True)
        netchop_score_dict = NetChopRunner.process_results(pd.concat(res_list, ignore_index=True))

        for id in idx:
            self.assertTrue(id in netchop_score_dict)

        for id in idx:
            if wt_peptides[id] != 'nan':
                self.assertEqual(len(mut_peptides[id]), len(netchop_score_dict[id]))

    def test_process_results2(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        res_list = netChopRunner.process_peptides(range(0, 1), ['LDRILQAGLDVERLSLRNFFHHFHS'], 'fa', ['a'], show_cmd=True)
        netchop_score_dict = NetChopRunner.process_results(pd.concat(res_list, ignore_index=True))

        self.assertEqual(len('LDRILQAGLDVERLSLRNFFHHFHS'), len(netchop_score_dict[0]))
        ref_scores = [0.800005, 0.027928, 0.262846, 0.158427, 0.754639,
                      0.029949, 0.751662, 0.057636, 0.559149, 0.027707,
                      0.571897, 0.031234, 0.093149, 0.754613, 0.026392,
                      0.836224, 0.295678, 0.054398, 0.128598, 0.385292,
                      0.384771, 0.595091, 0.975838, 0.735233, 0.290312]
        self.assertTrue(all([abs(s1 - s2) < 0.0001 for (s1, s2) in zip(ref_scores, netchop_score_dict[0])]))

    def test_process_results2a(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        mut_peptides = ['LDRILQAGLDVERLSLRNFFHHFHS', 'LDRILQAGLDVERLSLRNFFHHFH', 'RILQAGLDVERLSLRNFFHHFHS',
                        'ILQAGLDVERLSLRNFFHHFHS', 'LDRILQAGLDVERLSLRNFFH']
        res_list = netChopRunner.process_peptides(range(0, 5), mut_peptides, 'fa', ['a'], show_cmd=True)
        netchop_score_dict = NetChopRunner.process_results(pd.concat(res_list, ignore_index=True))

        for i in range(0, 5):
            self.assertEqual(len(mut_peptides[i]), len(netchop_score_dict[i]))

    def test_process_results3(self):
        test_scores = np.full(len('LDRILQAGLDVERLSLRNFFHHFHS'), 0)
        test_scores[12] = 1
        netchop_score_dict = {0: np.full(len('LDRILQAGLDVERLSLRNFFHHFHS'), 1),
                              1: np.arange(len('LDRILQAGLDVERLSLRNFFHHFHS')),
                              2: test_scores}

        idx = [0, 1, 2]
        mut_peptides = ['LDRILQAGLDVERLSLRNFFHHFHS', 'LDRILQAGLDVERLSLRNFFHHFHS', 'LDRILQAGLDVERLSLRNFFHHFHS']
        positions = [12, 12, 12]
        df = NetChopRunner.combine_netchop_results_long(idx, mut_peptides, positions, netchop_score_dict)
        self.assertEqual(3.5, df.loc[0, 'mut_netchop_score'])
        self.assertEqual(56.5, df.loc[1, 'mut_netchop_score'])
        self.assertEqual(3, df.loc[2, 'mut_netchop_score'])

    def test_get_extended_seq(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')
        start, end, ext_seq = netChopRunner.get_extended_seq('ENSP00000323421.3', 'KLKEKKERL', 'KSKEKKERL', 674)

        self.assertEqual(ext_seq.find('KLKEKKERL'), start)
        self.assertEqual(start + len('KLKEKKERL'), end)
        self.assertEqual(29, len(ext_seq))

    def test_process_results4(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        mgr = DataManager()
        data = mgr.get_processed_data('13LN', 'rt_netmhc_stab')

        idx, mut_peptides, mut_pos, wt_peptides, protein_ids = NetChopRunner.get_peptide_info_long(data)

        res_list = netChopRunner.process_peptides(idx, mut_peptides, 'fa', alleles=['a'], show_cmd=True)
        mut_score_dict = netChopRunner.process_results(pd.concat(res_list, ignore_index=True))

        nr_col = data.shape[1]
        res = netChopRunner.combine_netchop_results_long(idx, mut_peptides, mut_pos, mut_score_dict)
        data = NetChopRunner.merge_with_data(data, res)

        self.assertTrue(nr_col+1, data.shape[1])

    def test_add_features(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')
        stab_data = netChopRunner.add_features('0YM1', 'rt_netmhc_stab', 'chop', write_res=False)

    def test_get_peptide_info_short(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        mgr = DataManager()
        data = mgr.get_processed_data('13LN', 'rt', 'short')

        idx, starts, ends, ext_seqs = netChopRunner.get_peptide_info_short(data)

        for i in range(len(starts)):
            if starts[i] >= 0:
                self.assertEqual(ext_seqs[i].find(data.loc[idx[i], 'wt_seq']), starts[i])

    def test_add_features_short(self):
        netChopRunner = NetChopRunner(exe_path='/home/localadmin/Programmes/')

        mgr = DataManager()
        data = mgr.get_processed_data('0YM1', 'rt', 'short')

        data_netchop = netChopRunner.add_features_short(data)

        self.assertTrue(data_netchop.shape[1], data.shape[1] + 3)


