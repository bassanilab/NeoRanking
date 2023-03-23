from unittest import TestCase

import numpy as np

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

    def test_get_patient_info(self):
        annotator = NeoDiscImmunogenicityAnnotatorLong()
        p_info = annotator.get_patient_info('058C', ['HLA-I', 'HLA-I/II'],
                                            ['Predicted neo', 'MS (NEO-SNV)', 'Predicted Neo', 'Predicted Neo (SNV)'])

        self.assertEqual(166, p_info.shape[0])

    def test_annotate_patient_short(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorShort(mgr)

        for p in [p for p in mgr.get_valid_patients('short') if p in annotator.get_patients()]:
            data = annotator.annotate_patient(p)
            idx = data['response_type'] == 'CD8'
            if sum(idx) > 0:
                print(p+': '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type', 'response_annot']]))
            else:
                print(p+': no immunogenic peptides')

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

    def test_intersect(self):
        bool = NeoDiscImmunogenicityAnnotatorLong.intersect("LVKTDILAYLKQFKTK", "GPX5", "PVMRWSHRATVSLVKTDILAYLKQF", 13, "GPX5")
        self.assertTrue(bool)

        bool = NeoDiscImmunogenicityAnnotatorLong.intersect("VKTDILAYLKQFKTK", "GPX5", "PVMRWSHRATVSLVKTDILAYLKQF", 13, "GPX5")
        self.assertFalse(bool)

        bool = NeoDiscImmunogenicityAnnotatorLong.intersect("LVKTDILAYLKQFKTK", "GPX4", "PVMRWSHRATVSLVKTDILAYLKQF", 13, "GPX5")
        self.assertFalse(bool)

    def test_annotate_patient_long(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorLong(mgr)

        for p in [p for p in mgr.get_valid_patients('long') if p in annotator.get_patients()]:
            data = annotator.annotate_patient(p)
            idx = data['response_type'] == 'CD8'
            if sum(idx) > 0:
                print(p+': '+str(data.loc[idx, ['gene', 'mutant_seq', 'response_type', 'response_annot']]))
            else:
                print(p+': no immunogenic peptides')

    def test_annotate_patient_long_13LN(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorLong(mgr)

        annotator.update_immunogenicity_annotation("13LN", "rt_netmhc_stab_chop_tap_mbp_tcr",
                                                   "rt_netmhc_stab_chop_tap_mbp_tcr_newRT")

    def test_overwrite_hitide_immunogenicity_annotations_long(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorLong(mgr)

        annotator.overwrite_hitide_immunogenicity_annotations("rt_netmhc_stab_chop_tap_mbp_tcr")

    def test_annotate_patient_short_13LN(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorShort(mgr)

        annotator.update_immunogenicity_annotation("13LN", "rt_stab_chop_tap_mbp_tcr", "rt_stab_chop_tap_mbp_tcr_newRT")

    def test_overwrite_hitide_immunogenicity_annotations_short(self):
        mgr = DataManager()
        annotator = NeoDiscImmunogenicityAnnotatorShort(mgr)

        annotator.overwrite_hitide_immunogenicity_annotations("rt_stab_chop_tap_mbp_tcr")
