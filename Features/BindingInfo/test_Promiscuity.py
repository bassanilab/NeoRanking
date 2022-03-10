from CalcBindingAffinityScoreDists import *
from unittest import TestCase
from scipy import stats


class TestPromiscuity(TestCase):
    def test_match_sequence(self):

        fasta_file = path.join(Parameters().get_data_dir(), "tiny_fasta.fa")
        promiscuity = CalcBindingAffinityScoreDists(seq_logo_dir=Parameters().get_data_dir(), seq_fasta_file=fasta_file)
        protein_seq = 'MNGPVDGLCDHSLSEGVFMFTSESVGEGHPDKICDQISDAVLDAHLKQDP'

        sequenceLogo = promiscuity.sequenceLogoMgr.get_sequence_logo('A0101')
        data = sequenceLogo.get_freq_data(9)

        scores = promiscuity.match_sequence(data, 9, protein_seq)

        self.assertEqual(len(protein_seq)-9+1, len(scores))

    def test_match_fasta(self):

        fasta_file = path.join(Parameters().get_data_dir(), "tiny_fasta.fa")
        promiscuity = CalcBindingAffinityScoreDists(seq_logo_dir=Parameters().get_data_dir(), seq_fasta_file=fasta_file)

        scores = promiscuity.match_fasta('A0101', 9)

        self.assertEqual(1503, len(scores))

    def test_get_ipmsdb_proteins(self):

        proteins = DataManager().get_ipMSDB_proteins(100)

        self.assertEqual(1141, len(proteins))

    def test_get_ipmsdb_fasta_proteins(self):

        proteins = DataManager().get_ipMSDB_proteins(100)

        fasta_file = Parameters().get_human_fasta_file()
        promiscuity = CalcBindingAffinityScoreDists(seq_logo_dir=Parameters().get_data_dir(), seq_fasta_file=fasta_file,
                                                    proteins=proteins)

        self.assertEqual(1128, len(promiscuity.fasta_sequences))

    def test_match_fasta2(self):
        proteins = DataManager().get_ipMSDB_proteins(100)

        fasta_file = Parameters().get_human_fasta_file()
        promiscuity = CalcBindingAffinityScoreDists(seq_logo_dir=Parameters().get_data_dir(), seq_fasta_file=fasta_file,
                                                    proteins=proteins)

        scores = promiscuity.match_fasta('A0101', 9)

        print(stats.describe(scores))
