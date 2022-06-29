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


def get_valid_patients(patients, peptide_type='long'):
    dataManager = DataManager()
    patients_with_data = dataManager.get_valid_patients(peptide_type)
    if patients is None or len(patients) == 0:
        return patients_with_data
    if isinstance(patients, str):
        patients = [patients]

    patient_set = set()
    for p in patients:
        patient_set = patient_set.union(get_patients_from_group(p, peptide_type))

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


def read_cat_encodings(peptide_type='long'):
    encoding_file = Parameters().get_cat_to_num_info_file(peptide_type)
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


def get_patients_from_group(patient_group, data_manager, peptide_type='long'):
    patient_tag = patient_group.lower()

    if patient_tag == 'rosenberg':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients('all')
    elif patient_tag == 'gartner_train':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager) \
            if peptide_type == 'long' else RosenbergImmunogenicityAnnotatorShort(data_manager)
        patients = annotator.get_patients('gartner_train')
    elif patient_tag == 'gartner_test':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager)
        patients = annotator.get_patients('gartner_test')
    elif patient_tag == 'gartner':
        annotator = RosenbergImmunogenicityAnnotatorLong(data_manager)
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
    else:
        patients = set([patient_group])

    return patients


def get_patient_group(patient):
    regex = re.compile("\\d{4}")
    if patient.startswith('TESLA'):
        return 'TESLA'
    if regex.match(patient):
        return 'Gartner'
    else:
        return 'HiTIDE'
