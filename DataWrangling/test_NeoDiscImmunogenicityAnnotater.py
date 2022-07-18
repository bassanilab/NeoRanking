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

