import glob
import math
import os
import re
from shutil import which
import ast
from sklearn.preprocessing import QuantileTransformer, StandardScaler, PowerTransformer, \
    FunctionTransformer, MinMaxScaler

from DataWrangling.DataLoader import *
from DataWrangling.Transform_Data import Encoder
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *
from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorShort import *
from DataWrangling.TESLAImmunogenicityAnnotatorShort import *


def is_binary_file(filepathname):
    textchars = bytearray([7, 8, 9, 10, 12, 13, 27]) + bytearray(range(0x20, 0x7f)) + bytearray(range(0x80, 0x100))
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

    with open(filepathname, 'rb') as file:
        if is_binary_string(file.read(1024)):
            return True
        else:
            return False


def find_exe(path=None, exe=None):
    exe_path = None

    if exe is None:
        return None

    if path is None:
        exe_path = which(exe)
    else:
        for exe_file in glob.glob(os.path.join(path, '**', exe), recursive=True):
            exe_path = exe_file
            if is_binary_file(exe_path):
                break
            else:
                exe_path = None

    return exe_path


def get_processed_file(processed_dir, patient, tag, peptide_type='long'):
    return os.path.join(processed_dir, str(patient) + "_" + peptide_type + "_" + tag + '.txt')


def get_peptides_from_headers(fasta_ids, peptides):
    idx = list(map(lambda id: int(id.split(sep="_")[1]), fasta_ids))

    return peptides[idx]


def get_valid_patients(dataset, peptide_type='long'):
    dataManager = DataManager()
    patients_with_data = dataManager.get_valid_patients(peptide_type)
    if dataset is None or len(dataset) == 0:
        return patients_with_data
    if isinstance(dataset, str):
        dataset = [dataset]

    patient_set = set()
    for ds in dataset:
        patient_set = patient_set.union(get_patients_from_group(ds, peptide_type))

    return patient_set.intersection(patients_with_data)


def get_valid_annotated_patients(dataset, peptide_type='long'):
    dataManager = DataManager()
    patients_with_data = dataManager.get_valid_patients(peptide_type)
    if dataset is None or len(dataset) == 0:
        return patients_with_data
    if isinstance(dataset, str):
        dataset = [dataset]

    patient_set = set()
    for ds in dataset:
        patient_set = patient_set.union(get_annotated_patients_from_group(ds, dataManager, peptide_type))

    return patient_set.intersection(patients_with_data)


def get_immunogenic_patients(patients, peptide_type='long'):
    dataManager = DataManager()
    return set.intersection(patients, dataManager.get_immunogenic_patients(peptide_type))


def get_processed_patients(patients, file_tag, peptide_type='long'):
    data_loader = DataLoader(response_types=['not_tested', 'CD8', 'CD4/CD8', 'negative'], mutation_types=['SNV'],
                             immunogenic=['CD8', 'CD4/CD8'], min_nr_immuno=0)
    processed_patients = []
    for p in patients:
        data, X, y = data_loader.load_patients(p, file_tag, peptide_type)
        if data:
            processed_patients.append(p)

    return processed_patients


def get_normalizer(normalizer_tag):
    if normalizer_tag == 'q':
        return QuantileTransformer()
    elif normalizer_tag == 'z':
        return StandardScaler()
    elif normalizer_tag == 'p':
        return PowerTransformer()
    elif normalizer_tag == 'i':
        return MinMaxScaler()
    elif normalizer_tag == 'l':
        return FunctionTransformer(np.log10, inverse_func=lambda x: np.power(10, x), validate=True, check_inverse=True)
    elif normalizer_tag == 'a':
        return FunctionTransformer(np.arcsinh, inverse_func=np.sinh, validate=True, check_inverse=True)
    elif normalizer_tag == 'n':
        return None
    else:
        try:
            d = ast.literal_eval(normalizer_tag)
            for k, v in d.items():
                d[k] = get_normalizer(v)
            return d
        except:
            print('Cannot parse normalization dictionary {}'.format(normalizer_tag))
            return None


def read_cat_encodings(patient_set, peptide_type='long'):
    encoding_file = Parameters().get_cat_to_num_info_file(patient_set, peptide_type)
    encoding_df = pd.read_csv(encoding_file, header=0, sep="\t", comment='#')
    features = encoding_df['Feature'].unique()
    encoders = {}
    for f in features:
        encoder = Encoder(f)
        encoder.read_from_file(encoding_df)
        encoders[f] = encoder

    return encoders


def get_normalizer_name(normalizer):
    if type(normalizer) is QuantileTransformer:
        return "QuantileTransformer"
    elif type(normalizer) is PowerTransformer:
        return "PowerTransformer"
    elif type(normalizer) is StandardScaler:
        return "StandardScaler"
    elif type(normalizer) is MinMaxScaler:
        return "MinMaxScaler"
    elif type(normalizer) is FunctionTransformer and normalizer.func.__name__ == 'log10':
        return "log10"
    elif type(normalizer) is FunctionTransformer and normalizer.func.__name__ == 'arcsinh':
        return "arcsinh"
    elif normalizer is None:
        return "None"


def get_annotated_patients_from_group(patient_group, data_manager, peptide_type='long'):
    patient_tag = patient_group.lower()

    if patient_tag == 'rosenberg' or patient_tag == 'nci':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients('nci')
    elif patient_tag == 'gartner_train' or patient_tag == 'nci_train':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients('nci_train')
    elif patient_tag == 'gartner_test' or patient_tag == 'nci_test':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients(patient_tag)
    elif patient_tag == 'gartner':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients('gartner')
    elif patient_tag == 'parkhurst':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients('parkhurst')
    elif patient_tag == 'hitide':
        annotator = NeoDiscImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else NeoDiscImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients()
    elif patient_tag == 'tesla':
        annotator = TESLAImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else TESLAImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients()
    elif patient_tag == 'tesla/hitide':
        annotator = TESLAImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else TESLAImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients()
        annotator = NeoDiscImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else NeoDiscImmunogenicityAnnotatorShort(data_manager)
        patients = patients.union(annotator.get_patients())
    elif patient_tag == 'test':
        annotator = TESLAImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else TESLAImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients()
        annotator = NeoDiscImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else NeoDiscImmunogenicityAnnotatorShort(data_manager)
        patients = patients.union(annotator.get_patients())
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = patients.union(annotator.get_patients('gartner_test'))
    elif patient_tag == 'debug':
        patients = set(['058C', '0YM1'])
    else:
        patients = set([patient_group])

    return patients


def get_patients_from_group(patient_group, peptide_type='long'):
    patient_tag = patient_group.lower()

    if patient_tag == 'rosenberg' or patient_tag == 'nci':
        if peptide_type == 'long':
            patients = [4258,4285,4232,4234,4112,3978,4257,3737,4196,4072,4271,4237,4223,3309,4274,4180,4239,4230,4141,
                        4077,4246,4000,4220,4202,4071,4345,3881,4081,4268,4110,4270,4281,4007,4213,4189,3716,4231,4115,
                        4214,4060,4326,3971,4160,4298,4151,4253,4145,4323,3942,4236,4265,4272,3812,3948,4186,3998,4348,
                        4346,4240,4255,4126,4134,3784,4238,4217,4301,3995,4262,4129,3879,2556,4095,4278,4252,2098,3713,
                        3703,4259,4247,3678,4069,4284,4359,4200,4283,4316,4158,4207,4324,4014,3795,2591,1913,4263,4166,
                        4317,4245,4171,4136,4287,4244,4107,4275,4090,4242,2369,4350,3775,4155,4329,4338,4037]
        else:
            patients = [4166,4259,2556,4324,4252,4316,4301,4246,4275,4323,3995,3716,4270,2591,4326,4077,4242,4171,3713,
                        4298,4240,3879,3703,4271,4268,4134,1913,4348,3309,4155,4202,2369,3678,4000,3784,3775,4350,3942,
                        4110,4287,4196,4265,4359,4317,3795,3881,4126,4014,4136,4258,4213,4262,4071,4237,2098,4186,4007,
                        4329,4069,4160,4345,4338,4244,4346,4189,4253,4283,4238,3998,4281,4180,4081,4278,4095,4245,4284,
                        4158,4037,4234,4129]
    elif patient_tag == 'gartner_train' or patient_tag == 'nci_train':
        if peptide_type == 'long':
            patients = [4258,4259,4237,3812,4060,4257,4287,4095,4160,4107,4069,3978,4272,4301,4151,4244,4090,2098,4115,
                        4346,4275,2556,3971,4230,4214,4112,4220,4129,4189,4000,3784,4281,4207,4278,4329,4077,4247,4234,
                        3775,4317,4081,4134,1913,4274,4180,3309,3737,4126,4252,4246,4110,4037,4186,4236,4196,4255,4232,
                        4071,4223,2369,3998,3713,4284,3795,3879,4217,4239,4245,4263,4240,3948,4231,2591,3716,4316,4326,
                        4155,4238,3678,4213,4145,4136,4345,4141,4200,4298,4072,4158,4285]
        else:
            patients = [4136,4258,4259,4213,2556,4238,3998,4252,4134,1913,4281,4316,3309,4071,4301,4155,4237,4246,4180,
                        2369,3678,2098,4275,4186,4000,3784,4081,3775,4329,4278,4069,3716,2591,4095,4110,4326,4287,4077,
                        4196,4160,4345,3713,4317,4245,4284,3795,4298,4158,4240,4244,4037,4234,4346,4126,3879,4189,4129]
    elif patient_tag == 'gartner_test' or patient_tag == 'nci_test':
        if peptide_type == 'long':
            patients = [4268,4283,4270,4324,4014,4007,4262,3995,4271,4166,3703,4323,4171,4253,3942,4242,4265,4359,4350,
                        4202,4348,4338,3881]
        else:
            patients = [3703,4166,4283,4271,4262,4324,4268,4348,4202,4007,4323,3995,4350,4270,3942,4265,4242,4359,4171,
                        4338,3881,4253,4014]
    elif patient_tag == 'gartner':
        if peptide_type == 'long':
            patients = [4166,4259,2556,4324,4252,4316,4301,4246,4275,4323,3995,3716,4270,2591,4326,4077,4242,4171,3713,
                        4298,4240,3879,3703,4271,4268,4134,1913,4348,3309,4155,4202,2369,3678,4000,3784,3775,4350,3942,
                        4110,4287,4196,4265,4359,4317,3795,3881,4126,4014,4136,4258,4213,4262,4071,4237,2098,4186,4007,
                        4329,4069,4160,4345,4338,4244,4346,4189,4253,4283,4238,3998,4281,4180,4081,4278,4095,4245,4284,
                        4158,4037,4234,4129]
        else:
            patients = [4166,4259,2556,4324,4252,4316,4301,4246,4275,4323,3995,3716,4270,2591,4326,4077,4242,4171,3713,
                        4298,4240,3879,3703,4271,4268,4134,1913,4348,3309,4155,4202,2369,3678,4000,3784,3775,4350,3942,
                        4110,4287,4196,4265,4359,4317,3795,3881,4126,4014,4136,4258,4213,4262,4071,4237,2098,4186,4007,
                        4329,4069,4160,4345,4338,4244,4346,4189,4253,4283,4238,3998,4281,4180,4081,4278,4095,4245,4284,
                        4158,4037,4234,4129]
    elif patient_tag == 'hitide':
        if peptide_type == 'long':
            patients = ['1EDA', '0ZMN', '0YM1', '16UH', '1IJX', '1IKA', '058C', '14MH', '13P4', '13LN', '1HU3']
        else:
            patients = ['1EDA', '0ZMN', '0YM1', '16UH', '1IJX', '1IKA', '058C', '14MH', '13P4', '13LN', '1HU3']
    elif patient_tag == 'tesla':
        if peptide_type == 'long':
            patients = ['TESLA1', 'TESLA2', 'TESLA9', 'TESLA12', 'TESLA4', 'TESLA3', 'TESLA16', 'TESLA8']
        else:
            patients = ['TESLA1', 'TESLA2', 'TESLA9', 'TESLA12', 'TESLA4', 'TESLA3', 'TESLA16', 'TESLA8']
    elif patient_tag == 'tesla/hitide':
        if peptide_type == 'long':
            patients = ['1EDA', '0ZMN', '0YM1', '16UH', '1IJX', '1IKA', '058C', '14MH', '13P4', '13LN', '1HU3',
                        'TESLA1', 'TESLA2', 'TESLA9', 'TESLA12', 'TESLA4', 'TESLA3', 'TESLA16', 'TESLA8']
        else:
            patients = ['1EDA', '0ZMN', '0YM1', '16UH', '1IJX', '1IKA', '058C', '14MH', '13P4', '13LN', '1HU3',
                        'TESLA1', 'TESLA2', 'TESLA9', 'TESLA12', 'TESLA4', 'TESLA3', 'TESLA16', 'TESLA8']
    elif patient_tag == 'test':
        if peptide_type == 'long':
            patients = ['1EDA', '0ZMN', '0YM1', '16UH', '1IJX', '1IKA', '058C', '14MH', '13P4', '13LN', '1HU3',
                        'TESLA1', 'TESLA2', 'TESLA9', 'TESLA12', 'TESLA4', 'TESLA3', 'TESLA16', 'TESLA8',
                        4268, 4283,4270,4324,4014,4007,4262,3995,4271,4166,3703,4323,4171,4253,3942,4242,4265,4359,4350,
                        4202,4348,4338,3881]
        else:
            patients = ['1EDA', '0ZMN', '0YM1', '16UH', '1IJX', '1IKA', '058C', '14MH', '13P4', '13LN', '1HU3',
                        'TESLA1', 'TESLA2', 'TESLA9', 'TESLA12', 'TESLA4', 'TESLA3', 'TESLA16', 'TESLA8',
                        4268, 4283,4270,4324,4014,4007,4262,3995,4271,4166,3703,4323,4171,4253,3942,4242,4265,4359,4350,
                        4202,4348,4338,3881]
    elif patient_tag == 'debug':
        patients = ['058C', '0YM1']
    else:
        patients = []

    patients = set(np.array(patients, dtype=str))
    return patients


def get_patient_group(patient):
    regex = re.compile("\\d{4}")
    if patient.startswith('TESLA'):
        return 'TESLA'
    if regex.match(patient):
        return 'NCI'
    else:
        return 'HiTIDE'


def get_ml_group(patient, peptide_type):
    regex = re.compile("\\d{4}")
    if patient.startswith('TESLA'):
        return 'test'
    if regex.match(patient):
        if patient in get_patients_from_group('gartner_test', peptide_type):
            return "test"
        else:
            return "train"
    else:
        return 'test'
