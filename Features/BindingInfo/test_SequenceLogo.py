from unittest import TestCase

from Features.BindingInfo.SequenceLogo import *


class TestSequenceLogoMgr(TestCase):
    def test_get_sequence_logo(self):
        mgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir())
        sequence_logo = mgr.get_sequence_logo("A0101")

        self.assertEqual((20, 9), sequence_logo.get_freq_data(9).shape)
        self.assertEqual((20, 10), sequence_logo.get_freq_data(10).shape)
        self.assertEqual((20, 11), sequence_logo.get_freq_data(11).shape)
        self.assertEqual((20, 12), sequence_logo.get_freq_data(12).shape)
        self.assertEqual((20, 13), sequence_logo.get_freq_data(13).shape)
        self.assertEqual((20, 14), sequence_logo.get_freq_data(14).shape)

    def test_is_binding_pos(self):
        mgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        sequence_logo = mgr.get_sequence_logo("A0101")
        self.assertFalse(sequence_logo.is_binding_pos(9, 1))
        self.assertTrue(sequence_logo.is_binding_pos(9, 2))
        self.assertTrue(sequence_logo.is_binding_pos(9, 3))
        self.assertTrue(sequence_logo.is_binding_pos(9, 9))
        self.assertFalse(sequence_logo.is_binding_pos(9, 4))
        self.assertFalse(sequence_logo.is_binding_pos(9, 8))

        sequence_logo = mgr.get_sequence_logo("A0201")
        self.assertFalse(sequence_logo.is_binding_pos(9, 1))
        self.assertTrue(sequence_logo.is_binding_pos(9, 2))
        self.assertFalse(sequence_logo.is_binding_pos(9, 4))
        self.assertTrue(sequence_logo.is_binding_pos(9, 9))

    def test_process_peptides(self):
        data = DataManager().get_processed_data('0YM1', 'rt_netmhc_stab_chop')
        if 'mut_peptide_pos_0' not in data.columns:
            print('Data file %s does not contain netMHC predictions. Unable to add netchop prediction.')
            return None

        mgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        mut_peptides, mut_aas, wt_aas, mut_alleles, mut_position = SequenceLogoMgr.get_peptide_info_long(data, 0)
        mut_binding_score, mut_is_binding_pos = mgr.process_peptides(mut_peptides, mut_aas, mut_alleles, mut_position)

        self.assertTrue(all(np.array([2, 8]) == np.where(mut_is_binding_pos)[0]))

    def test_process_peptides(self):
        mgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        mut_binding_score, mut_is_binding_pos = mgr.process_peptides(['AKQTFCFIN'], ['I'], ['B3901,C0702'], [8])
        mut_binding_score_0, mut_is_binding_pos_0 = mgr.process_peptides(['AKQTFCFIN'], ['I'], ['B3901'], [8])
        mut_binding_score_1, mut_is_binding_pos_1 = mgr.process_peptides(['AKQTFCFIN'], ['I'], ['C0702'], [8])

        self.assertEqual(max(mut_binding_score_0[0], mut_binding_score_1[0]), mut_binding_score[0])
        self.assertEqual(max(mut_binding_score_0[0], mut_binding_score_1[0]), mut_binding_score[0])

    def test_add_features_long(self):
        sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        data = sequenceLogoMgr.add_features('4207', 'rt_netmhc_stab_chop', 'mbp', write_res=False)

        self.assertTrue('mut_is_binding_pos_0' in data.columns)
        self.assertTrue('wt_binding_score_0' in data.columns)
        self.assertTrue('mut_binding_score_0' in data.columns)

        data = sequenceLogoMgr.add_features('3737', 'rt_netmhc_stab_chop', 'mbp', write_res=False)

        self.assertTrue('mut_is_binding_pos_0' in data.columns)
        self.assertTrue('wt_binding_score_0' in data.columns)
        self.assertTrue('mut_binding_score_0' in data.columns)

        data = sequenceLogoMgr.add_features('4160', 'rt_netmhc_stab_chop', 'mbp', write_res=False)

        self.assertTrue('mut_is_binding_pos_0' in data.columns)
        self.assertTrue('wt_binding_score_0' in data.columns)
        self.assertTrue('mut_binding_score_0' in data.columns)

        data = sequenceLogoMgr.add_features('4107', 'rt_netmhc_stab_chop', 'mbp', write_res=False)

        self.assertTrue('mut_is_binding_pos_0' in data.columns)
        self.assertTrue('wt_binding_score_0' in data.columns)
        self.assertTrue('mut_binding_score_0' in data.columns)

    def test_add_features_short(self):
        sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        data = sequenceLogoMgr.add_features('14MH', 'rt', 'mbp', 'short', write_res=False)

        self.assertTrue('mut_is_binding_pos' in data.columns)
        self.assertTrue('mut_binding_score' in data.columns)
        self.assertTrue(all([mbp in [True, False, np.nan] for mbp in data.mut_is_binding_pos.unique()]))

    def test_add_features_short2(self):
        sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        data = sequenceLogoMgr.add_features('2556', 'rt', 'mbp', 'short', write_res=False)

        self.assertTrue('mut_is_binding_pos' in data.columns)
        self.assertTrue('mut_binding_score' in data.columns)
        self.assertTrue(all([mbp in [True, False, np.nan] for mbp in data.mut_is_binding_pos.unique()]))

    def test_get_aa_score(self):
        sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir(), binding_threshold=20)
        value = sequenceLogoMgr.get_aa_score('C0702', len('PGPSDQPSH'), 'Q', 9)

        self.assertEqual(-2.6500903606166855, value)

    def test_alleles_exist(self):
        mgr = DataManager()
        alleles = mgr.get_all_classI_allotypes()

        sequenceLogoMgr = SequenceLogoMgr(seq_logo_dir=Parameters().get_data_dir())
        for a in alleles:
            print(a)
            sequenceLogo = sequenceLogoMgr.get_sequence_logo(a)
            if not sequenceLogo:
                print("allele {0} is missing.".format(a))




