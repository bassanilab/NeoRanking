from unittest import TestCase
from NeoDiscImmunogenicityAnnotatorShort import *
from NeoDiscImmunogenicityAnnotatorLong import *


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

