from unittest import TestCase
from TESLAImmunogenicityAnnotatorLong import *
from TESLAImmunogenicityAnnotatorShort import *


class TestTESLAFileConverter(TestCase):

    def test_constructor(self):
        annotator = TESLAImmunogenicityAnnotatorLong()
        self.assertEqual((918, 3), annotator.immuno_data.shape)

    def test_intersect(self):
        seq_short = 'THEAVENGER'
        seq_long = ['PETIDE', 'AVEN', "AVENT", "TESTTEST"]

        idx = TESLAImmunogenicityAnnotatorLong.intersect(seq_short, 6, seq_long)

        self.assertEqual([1], idx)

    def test_annotate_patient_long(self):
        annotator = TESLAImmunogenicityAnnotatorLong()

        data = annotator.annotate_patient('TESLA1')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA1: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(10, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA2')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA2: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(4, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA3')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA3: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(12, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA4')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA4: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(1, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA8')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA8: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(1, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA9')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA9: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(2, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA12')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA12: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA16')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA16: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))

    def test_get_patients_long(self):
        annotator = TESLAImmunogenicityAnnotatorLong()
        self.assertEqual(9, len(annotator.get_patients()))

    def test_annotate_patient_short(self):
        annotator = TESLAImmunogenicityAnnotatorShort()

        data = annotator.annotate_patient('TESLA1')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA1: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(9, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA2')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA2: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(4, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA3')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA3: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(11, sum(data['response_type'] == 'CD8'))
        # 2 peptides KIYTGEKPYK and KIVEMSTSK are missing in neodisc:
        # KIVEMSTSK: mutation nor detected
        # KIYTGEKPYK: shared peptide overlap with other genes

        data = annotator.annotate_patient('TESLA4')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA4: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(1, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA8')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA8: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(1, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA9')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA9: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(2, sum(data['response_type'] == 'CD8'))

        data = annotator.annotate_patient('TESLA12')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA12: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))
        # immunogenic peptide YQANVVWKV missing in neodisc (frameshift)

        data = annotator.annotate_patient('TESLA16')
        idx = data['response_type'] == 'CD8'
        if sum(idx) > 0:
            print('TESLA16: '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))
        # immunogenic peptide YLNEAVFNFV missing in neodisc (no mutation)

    def test_get_patients_short(self):
        annotator = TESLAImmunogenicityAnnotatorLong()
        self.assertEqual(9, len(annotator.get_patients()))
