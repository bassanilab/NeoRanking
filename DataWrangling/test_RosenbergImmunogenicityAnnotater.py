from unittest import TestCase
from RosenbergImmunogenicityAnnotatorLong import *
from RosenbergImmunogenicityAnnotatorShort import *


class TestRosenbergFileConverter(TestCase):

    @classmethod
    def setUpClass(cls):
        ...

    def test_match(self):
        self.assertEqual(True, RosenbergImmunogenicityAnnotatorLong.match('hans', 'hans', offset=1))
        self.assertEqual(True, RosenbergImmunogenicityAnnotatorLong.match('hans', 'hans', offset=0))
        self.assertEqual(True, RosenbergImmunogenicityAnnotatorLong.match('harald', 'hharal', offset=1))
        self.assertEqual(True, RosenbergImmunogenicityAnnotatorLong.match('ttamara', 'tamaraa', offset=1))
        self.assertEqual(False, RosenbergImmunogenicityAnnotatorLong.match('ttamara', 'tamaraa', offset=0))

    def test_intersect(self):
        lst = ['hans', 'hharal', 'tamaraa', 'george']

        self.assertEqual([0], RosenbergImmunogenicityAnnotatorLong.intersect('hans', lst, offset=1))
        self.assertEqual([0], RosenbergImmunogenicityAnnotatorLong.intersect('hans', lst, offset=0))
        self.assertEqual([1], RosenbergImmunogenicityAnnotatorLong.intersect('harald', lst, offset=1))
        self.assertEqual([2], RosenbergImmunogenicityAnnotatorLong.intersect('ttamara', lst, offset=1))
        self.assertEqual([], RosenbergImmunogenicityAnnotatorLong.intersect('ttamara', lst, offset=0))
        self.assertEqual([], RosenbergImmunogenicityAnnotatorLong.intersect('patricia', lst, offset=1))
        self.assertEqual([], RosenbergImmunogenicityAnnotatorLong.intersect('tamarra', lst, offset=1))

    def test_annotate_gartner_train1(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_train('2098')

        self.assertEqual(3, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_train2(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_train('3309')

        self.assertEqual(2, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_train3(self):

        rosenberg_info_file = Parameters().get_gartner_info_excel_file()
        data_train = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 1', header=1)

        converter = RosenbergImmunogenicityAnnotatorLong()
        for p in converter.mgr.get_valid_patients():
            if p in converter.gartner_patients_train:
                data = converter.annotate_gartner_train(p)
                idx = data['response_type'] == 'CD8'
                s = data_train.loc[data_train['ID'] == int(p), 'CD8+ nmers']
                gartner_mut_cnt = s.iloc[0] if len(s) > 0 else -1
                print("Patient {0} has {1} CD8 mutations. {2} mutations in gartner et al.".
                      format(p, sum(idx), gartner_mut_cnt))
                if sum(idx) > 0:
                    print(data.loc[idx, ['gene', 'mutant_seq', 'response_type']])

    def test_annotate_gartner_train4(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_train('1913')

        self.assertEqual(0, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_test5(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_test('3703')

        self.assertEqual(4, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_test6(self):

        rosenberg_info_file = Parameters().get_gartner_info_excel_file()
        data_test = pd.read_excel(open(rosenberg_info_file, 'rb'), sheet_name='Supplementary Table 18', header=1)

        converter = RosenbergImmunogenicityAnnotatorLong()
        for p in converter.mgr.get_valid_patients():
            if p in converter.gartner_patients_test:
                data = converter.annotate_gartner_test(p)
                idx = data['response_type'] == 'CD8'
                s = data_test.loc[data_test['ID'] == int(p), 'CD8+ nmers']
                gartner_mut_cnt = s.iloc[0] if len(s) > 0 else -1
                print("Patient {0} has {1} CD8 mutations. {2} mutations in gartner et al.".
                      format(p, sum(idx), gartner_mut_cnt))
                if sum(idx) > 0:
                    print(data.loc[idx, ['gene', 'mutant_seq', 'response_type']])

    def test_annotate_parkhurst_test7(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_test('4007')

        self.assertEqual(2, sum(data['response_type'] == 'CD8'))

    def test_annotate_parkhurst_test1(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_test('abcd')

        self.assertEqual(None, data)

    def test_annotate_parkhurst_test8(self):

        parkhurst_info_file = Parameters().get_parkhurst_info_file()
        data_ref = pd.read_excel(open(parkhurst_info_file, 'rb'), sheet_name='Supplementary 3', header=0)

        converter = RosenbergImmunogenicityAnnotatorLong()
        for p in converter.mgr.get_valid_patients():
            if p in converter.parkhurst_patients:
                data = converter.annotate_parkhurst(p)
                idx = data['response_type'] == 'CD8'
                parkhurst_mut_cnt = \
                    sum(['CD8' in res for res in data_ref.loc[data_ref['ID'] == int(p),
                                                              'CD8/CD4 Nmer screening results']])
                print("Patient {0} has {1} CD8 mutations. {2} mutations in parkhurst et al.".
                      format(p, sum(idx), parkhurst_mut_cnt))
                if sum(idx) > 0:
                    print(data.loc[idx, ['gene', 'mutant_seq', 'response_type']])

    def test_get_patients(self):
        annotator = RosenbergImmunogenicityAnnotatorLong()
        self.assertEqual(70, len(annotator.get_patients('gartner_train')))
        self.assertEqual(26, len(annotator.get_patients('gartner_test')))
        self.assertEqual(75, len(annotator.get_patients('parkhurst')))
        self.assertEqual(137, len(annotator.get_patients('all')))

        mgr = DataManager()
        cnt = 0
        for p in mgr.get_valid_patients():
            if p in annotator.get_patients('parkhurst') and p not in annotator.get_patients('gartner_train') and \
                    p not in annotator.get_patients('gartner_test'):
                cnt += 1
        print(cnt)

    def test_annotate_gartner_short1(self):
        converter = RosenbergImmunogenicityAnnotatorShort()
        data = converter.annotate('2098', 'train')

        self.assertEqual(0, sum(data['response_type'] == 'no_mutation_found'))
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))  # 1 mutated peptide not in neodisc file
        self.assertEqual(147, sum(data['response_type'] == 'negative'))

    def test_annotate_gartner_short2(self):
        converter = RosenbergImmunogenicityAnnotatorShort()
        data = converter.annotate('3713', 'train')

        print(data.loc[data['response_type'] == 'CD8', 'mutant_seq'])

        self.assertEqual(0, sum(data['response_type'] == 'no_mutation_found'))
        # this should be 11 CD8 peptides, but one matches to 2 long peptides (1 immuno, one not tested),
        # but is only given for the not_tested one in neodisc output files
        self.assertEqual(10, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_short3(self):
        converter = RosenbergImmunogenicityAnnotatorShort()
        data = converter.annotate('4253', 'test')

        self.assertEqual(0, sum(data['response_type'] == 'no_mutation_found'))
        self.assertEqual(0, sum(data['response_type'] == 'CD8'))
        self.assertEqual(7400, sum(data['response_type'] == 'negative'))

    def test_annotate_gartner_short4(self):
        converter = RosenbergImmunogenicityAnnotatorShort()
        data = converter.annotate('4284', 'train')

        self.assertEqual(0, sum(data['response_type'] == 'no_mutation_found'))
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))  # 1 mutated peptide not in neodisc file
        self.assertEqual(147, sum(data['response_type'] == 'negative'))

    def test_get_patients_long(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        patients = converter.get_patients('all')
        self.assertEqual(137, len(patients))

        mgr = DataManager()
        self.assertEqual(112, len(patients.intersection(mgr.get_valid_patients('short'))))

    def test_get_patients_short(self):
        converter = RosenbergImmunogenicityAnnotatorShort()
        patients = converter.get_patients('all')
        self.assertEqual(96, len(patients))
        mgr = DataManager()
        self.assertEqual(80, len(patients.intersection(mgr.get_valid_patients('short'))))

    def test_annotate_gartner_train9(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_train('4136')
        self.assertEqual(1, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_train10(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_train('3998')
        self.assertEqual(3, sum(data['response_type'] == 'CD8'))

    def test_annotate_gartner_train11(self):
        converter = RosenbergImmunogenicityAnnotatorLong()
        data = converter.annotate_gartner_train('4346')

        rt = data.loc[data['mutant_seq'] == 'YTVMSYDRYLAICYPLRYTSMMSGS', 'response_type']
        self.assertTrue(all(rt == 'negative'))
        self.assertEqual(1, sum(data['response_type'] == 'CD8'))

    def test_intersect2(self):
        lst = ['VLLMPYGYVLNEFQSRQNSSSAQGS', 'SLLPEFVVPYMIYLLAHDPDFTRSQ']
        self.assertEqual([], RosenbergImmunogenicityAnnotatorLong.intersect('VLLMPYGYVLNEFQSCQNSSSAQGS', lst, offset=0))
        self.assertEqual(1, RosenbergImmunogenicityAnnotatorLong.intersect('SLLPEFVVPYMIYLLAHDPDFTRSQ', lst, offset=0))

