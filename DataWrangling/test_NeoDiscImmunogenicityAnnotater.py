from unittest import TestCase
from NeoDiscImmunogenicityAnnotatorShortOld import *
from NeoDiscImmunogenicityAnnotatorLong import *
from NeoDiscImmunogenicityAnnotatorShort import *
from NeoDiscImmunogenicityAnnotatorLongOld import *


class TestNeoDiscFileConverter(TestCase):

    def test_constructor(self):
        converter = NeoDiscImmunogenicityAnnotatorLong()
        converter = NeoDiscImmunogenicityAnnotatorShort()

    def test_get_patients(self):
        annotator = NeoDiscImmunogenicityAnnotatorLong()
        self.assertEqual(17, len(annotator.get_patients()))

    def test_annotate_patient_short(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorShort(mgr)

        for p in mgr.get_valid_patients():
            if p in annotator.get_patients():
                data = annotator.annotate_patient(p)
                idx = data['response_type'] == 'CD8'
                if sum(idx) > 0:
                    print(p+': '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
    #            self.assertEqual(9, sum(data['response_type'] == 'CD8'))

    def test_annotate_patient_short_13LN(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorShort(mgr)

        data = annotator.annotate_patient("13LN")

        CD8_cnt = sum(data['response_type'] == 'CD8')
        negative_cnt = sum(data['response_type'] == 'negative')
        not_tested_cnt = sum(data['response_type'] == 'not_tested')
        self.assertEqual(8, CD8_cnt)
        self.assertEqual(159, negative_cnt)
        self.assertEqual(11790, not_tested_cnt)

    def test_annotate_patient_long_13LN(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorLong(mgr)

        data = annotator.annotate_patient("13LN")

        CD8_cnt = sum(data['response_type'] == 'CD8')
        negative_cnt = sum(data['response_type'] == 'negative')
        not_tested_cnt = sum(data['response_type'] == 'not_tested')
        self.assertEqual(4, CD8_cnt)
        self.assertEqual(62, negative_cnt)
        self.assertEqual(168, not_tested_cnt)

    def test_annotate_patient_long_13P4(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorLong(mgr)

        data = annotator.annotate_patient("13P4")

        CD8_cnt = sum(data['response_type'] == 'CD8')
        negative_cnt = sum(data['response_type'] == 'negative')
        not_tested_cnt = sum(data['response_type'] == 'not_tested')
        self.assertEqual(5, CD8_cnt)
        self.assertEqual(30, negative_cnt)
        self.assertEqual(1066, not_tested_cnt)

    def test_constructor2(self):
        converter = NeoDiscImmunogenicityAnnotatorLongOld()
        converter = NeoDiscImmunogenicityAnnotatorShortOld()

    def test_get_patients2(self):
        annotator = NeoDiscImmunogenicityAnnotatorLongOld()
        self.assertEqual(28, len(annotator.get_patients()))

    def test_annotate_patient_short2(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorShortOld(mgr)

        for p in [p for p in mgr.get_valid_patients() if p in annotator.get_patients()]:
            data = annotator.annotate_patient(p)
            idx = data['response_type'] == 'CD8'
            if sum(idx) > 0:
                print(p+': '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type']]))
    #            self.assertEqual(9, sum(data['response_type'] == 'CD8'))

