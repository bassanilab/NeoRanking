import glob
from os import path
import numpy as np
import pandas as pd
from Utils.Parameters import *
import datetime
import re
import warnings
from pandas.errors import DtypeWarning


class DataManager:

    def __init__(self, rna_seq=True, netmhc=True, immunogenity=True):
        self.parameters = Parameters()
        self.original_data_dict = {}
        self.processed_data_dict = {}
        self.allotypes = pd.read_csv(filepath_or_buffer=self.parameters.get_allotype_file(), sep="\t", header=0)
        self.patient_data_file_dict = {}
        self.fill_in_patient_file_dict()
        self.valid_patients = {'long': [], 'short': []}
        self.set_valid_patients(rna_seq, netmhc)
        self.immunogenic_patients = {'long': [], 'short': []}
        if immunogenity:
            self.set_immunogenic_patients()

        return

    def fill_in_patient_file_dict(self):
        self.patient_data_file_dict['long'] = {}
        data_files = glob.glob(os.path.join(self.parameters.get_data_dir(), 'tabdata', '*_long.txt'))
        for f in data_files:
            with open(f) as dafi:
                p = re.split(r'[_-]', os.path.basename(f))[0]
                self.patient_data_file_dict['long'][p] = f

        self.patient_data_file_dict['short'] = {}
        data_files = glob.glob(os.path.join(self.parameters.get_data_dir(), 'tabdata', '*_short.txt'))
        for f in data_files:
            with open(f) as dafi:
                p = re.split(r'[_-]', os.path.basename(f))[0]
                self.patient_data_file_dict['short'][p] = f

    def set_valid_patients(self, rna_seq, netmhc):
        warnings.filterwarnings(action='ignore', category=DtypeWarning)
        data_info_file = os.path.join(self.parameters.get_data_dir(), "Data_validity_info.txt")
        if not os.path.isfile(data_info_file):
            self.check_neodisc_files(data_info_file)

        validity_info = pd.read_csv(data_info_file, sep='\t', header=0)
        validity_info.astype({'Patient': str, 'RNA_seq': bool, 'netmhc_short': bool})
        for i in validity_info.index:
            if (rna_seq and validity_info.loc[i, 'RNA_seq']) or not rna_seq:
                self.valid_patients['long'].append(validity_info.loc[i, 'Patient'])
                if (netmhc and validity_info.loc[i, 'netmhc_short']) or not netmhc:
                    self.valid_patients['short'].append(validity_info.loc[i, 'Patient'])

    def set_immunogenic_patients(self):
        data_info_file = os.path.join(self.parameters.get_data_dir(), "Data_immunogenicity_info.txt")
        if not os.path.isfile(data_info_file):
            self.check_immunogenicity(data_info_file)

        immunogenicity_info = pd.read_csv(data_info_file, sep='\t', header=0)
        immunogenicity_info.astype({'Patient': str, 'CD8_cnt_long': int, 'CD4_cnt_long': int, 'CD8_cnt_short': int,
                                    'CD4_cnt_short': int})
        for i in immunogenicity_info.index:
            if immunogenicity_info.loc[i, 'CD8_cnt_long'] > 0:
                self.immunogenic_patients['long'].append(immunogenicity_info.loc[i, 'Patient'])
            if immunogenicity_info.loc[i, 'CD8_cnt_short'] > 0:
                self.immunogenic_patients['short'].append(immunogenicity_info.loc[i, 'Patient'])

    def get_processed_file(self, patient, tag, peptide_type='long'):
        data_file = glob.glob(
            path.join(self.parameters.get_result_dir(), patient + "_" + peptide_type + "_" + tag + ".txt"))
        if len(data_file) > 0:
            return data_file[0]
        else:
            return None

    def get_original_data(self, patient, peptide_type='long'):
        if peptide_type not in self.original_data_dict:
            self.original_data_dict[peptide_type] = {}

        patient = str(patient)
        if patient in self.original_data_dict[peptide_type]:
            return self.original_data_dict[peptide_type][patient]
        else:
            if patient not in self.patient_data_file_dict[peptide_type]:
                print("Patient {} does not have data file.".format(patient))
                return None
            neoDisc_file = self.patient_data_file_dict[peptide_type][patient]
            data = pd.read_csv(neoDisc_file, header=0, sep="\t")

            if peptide_type == 'short' and 'Peptide_Class' in data.columns:
                data = data.loc[data.apply(lambda row: row['Peptide_Class'] in ['HLA_I', 'HLA_I_II'], axis=1)]

            data.reset_index(inplace=True, drop=True)
            self.original_data_dict[peptide_type][patient] = data

        return data

    def get_processed_data(self, patient, tag, peptide_type='long'):
        patient = str(patient)
        if peptide_type in self.processed_data_dict and tag in self.processed_data_dict[peptide_type] and \
                patient in self.processed_data_dict[peptide_type][tag]:
            return self.processed_data_dict[peptide_type][tag][patient]
        else:
            data_file = self.get_processed_file(patient, tag, peptide_type)
            if data_file is None:
                print("No data file found for peptide type {0}, patient {1} and tag {2}".format(peptide_type, patient,
                                                                                                tag))
                return None

            data = pd.read_csv(data_file, header=0, sep="\t")

            self.put_processed_data(data, patient, tag, peptide_type)

            return data

    def put_processed_data(self, data, patient, tag, peptide_type='long'):
        if peptide_type not in self.processed_data_dict:
            self.processed_data_dict[peptide_type] = {}
        if tag not in self.processed_data_dict[peptide_type]:
            self.processed_data_dict[peptide_type][tag] = {}

        self.processed_data_dict[peptide_type][tag][patient] = data

    def get_valid_patients(self, peptide_type='long'):
        return set(np.sort(self.valid_patients[peptide_type]))

    def get_immunogenic_patients(self, peptide_type='long'):
        return set(np.sort(self.immunogenic_patients[peptide_type]))

    def get_classI_allotypes(self, patient):

        patient = str(patient)
        p = np.array(self.allotypes['Patient'])
        idx, = np.where(p == patient)

        if len(idx) > 0:
            a = str(self.allotypes.loc[idx[0], 'Alleles'])
            return a.split(sep=",")
        else:
            return []

    def get_all_classI_allotypes(self):

        patients = np.array(self.allotypes['Patient'])

        alleles = set()
        for p in patients:
            idx, = np.where(p == patients)
            a = str(self.allotypes.loc[idx[0], 'Alleles'])
            for a in str(self.allotypes.loc[idx[0], 'Alleles']).split(sep=","):
                alleles.add(a)

        return alleles

    def check_neodisc_files(self, data_info_file):
        self.fill_in_patient_file_dict()

        data_validity_info = []
        for p in self.patient_data_file_dict['long'].keys():
            data_long = self.get_original_data(p, 'long')
            rna_seq = data_long is not None and 'rnaseq_TPM' in data_long.columns
            data_short = self.get_original_data(p, 'short')
            netmhc_short = data_short is not None and 'mutant_rank_netMHCpan' in data_short.columns and \
                data_short['mutant_rank_netMHCpan'].unique().shape[0] > 1

            data_validity_info.append([p, rna_seq, netmhc_short])

        pd.DataFrame(data_validity_info, columns=['Patient', 'RNA_seq', 'netmhc_short']). \
            to_csv(path_or_buf=data_info_file, sep="\t", index=False, header=True)

    def check_immunogenicity(self, data_info_file):
        data_validity_info = []
        for p in self.get_valid_patients():
            values = [p]
            for pt in ['long', 'short']:
                data = self.get_processed_data(p, "rt", pt)
                if data is not None:
                    cnt_I = data.apply(lambda row: row['response_type'] in ['CD8', 'CD4/CD8'], axis=1).sum()
                    cnt_II = data.apply(lambda row: row['response_type'] in ['CD4', 'CD4/CD8'], axis=1).sum()
                    values.append(cnt_I)
                    values.append(cnt_II)
                else:
                    values.append(0)
                    values.append(0)
            data_validity_info.append(values)

        pd.DataFrame(data_validity_info,
                     columns=['Patient', 'CD8_cnt_long', 'CD4_cnt_long', 'CD8_cnt_short', 'CD4_cnt_short']). \
            to_csv(path_or_buf=data_info_file, sep="\t", index=False, header=True)

    def get_classifier_file(self, clf_tag, peptide_type):
        if clf_tag == 'DNN':
            ext = 'h5'
        elif clf_tag == 'CatBoost':
            ext = 'cbm'
        elif clf_tag == 'XGBoost':
            ext = 'xgbm'
        elif clf_tag == 'TabNet':
            ext = 'tnm'
        else:
            ext = 'sav'

        date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
        file_name = '{0}_{1}_{2}.{3}'.format(clf_tag, peptide_type, date_time_str, ext)
        classifier_file = path.join(self.parameters.get_pickle_dir(), file_name)
        while os.path.isfile(classifier_file):
            date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
            file_name = '{0}_{1}_{2}.{3}'.format(clf_tag, peptide_type, date_time_str, ext)
            classifier_file = path.join(self.parameters.get_pickle_dir(), file_name)

        return classifier_file

    def get_classifier_result_file(self, clf_tag, peptide_type):

        date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
        file_name = '{0}_{1}_{2}_clf_results.txt'.format(clf_tag, peptide_type, date_time_str)
        classifier_file = path.join(self.parameters.get_pickle_dir(), file_name)
        while os.path.isfile(classifier_file):
            date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
            file_name = '{0}_{1}_{2}_clf_results.txt'.format(clf_tag, peptide_type, date_time_str)
            classifier_file = path.join(self.parameters.get_pickle_dir(), file_name)

        return classifier_file
