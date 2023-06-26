import datetime
import argparse
import os
from os import path
import numpy as np

from Utils.GlobalParameters import GlobalParameters


def is_cont_type(type_):
    '''
    Plot as continuous values, but natural order
    :param type_:
    :return: True if type_ any 'floatxx' or 'intxx' (xx>8)
    '''
    return type_.startswith('float') or type_ in ['int16', 'int32', 'int64']


def is_discrete_ordered_type(type_):
    '''
    Plot as discrete values, but keep int or bool order
    :param type_:
    :return: True if type_ is 'int8' or 'bool'
    '''
    return type_ in ['int8', 'bool']


def is_cat_type(type_):
    '''
    Plot as discrete values (order may be set or not)
    :param type_:
    :return: True if type_ is 'int8', 'bool' or 'category'
    '''
    return type_ in ['category']


def calc_plot_type(peptide_type: str, feature_: str):
    org_type_ = GlobalParameters.feature_types_neopep[feature_] if peptide_type == 'neopep' else \
        GlobalParameters.feature_types_mutation[feature_]
    return 'float64' if is_cont_type(org_type_) else org_type_


def get_processed_types(peptide_type: str, objective: str):
    if objective == 'ml':
        if peptide_type == 'neopep':
            type_dict = \
                {feature: 'float64' for feature in GlobalParameters.ml_features_neopep}
            type_dict['seq_len'] = 'int8'
            return type_dict
        elif peptide_type == 'mutation':
            type_dict = \
                {feature: 'float64' for feature in GlobalParameters.ml_features_mutation}
            return type_dict
        else:
            return None
    elif objective == 'plot':
        if peptide_type == 'neopep':
            type_dict = \
                {feature: calc_plot_type('neopep', feature) for feature in GlobalParameters.feature_types_neopep}
            return type_dict
        elif peptide_type == 'mutation':
            type_dict = \
                {feature: calc_plot_type('mutation', feature) for feature in GlobalParameters.feature_types_mutation}
            return type_dict
        else:
            return None
    else:
        return None


def get_classifier_ext(self, clf_tag):
    if clf_tag == 'DNN':
        return 'h5'
    elif clf_tag == 'CatBoost':
        return 'cbm'
    elif clf_tag == 'XGBoost':
        return 'xgbm'
    elif clf_tag == 'TabNet':
        return 'tnm'
    elif clf_tag == 'Threshold':
        return 'txt'
    else:
        return 'sav'


def get_classifier_file(self, clf_tag, base_name):
    if clf_tag == 'DNN':
        ext = 'h5'
    elif clf_tag == 'CatBoost':
        ext = 'cbm'
    elif clf_tag == 'XGBoost':
        ext = 'xgbm'
    elif clf_tag == 'TabNet':
        ext = 'tnm'
    elif clf_tag == 'Threshold':
        ext = 'txt'
    else:
        ext = 'sav'

    file_name = '{0}.{1}'.format(base_name, ext)
    classifier_file = path.join(self.parameters.get_pickle_dir(), file_name)

    return classifier_file


def get_result_file(self, prefix, run_id, peptide_type='', ext='txt', result_type='clf', suffix="results"):
    if result_type == 'clf':
        file_dir = self.parameters.get_pickle_dir()
    else:
        file_dir = self.parameters.get_plot_dir()

    date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
    if peptide_type != '':
        file_name = '{0}_{1}_{2}_{3}_clf_{4}.{5}'.format(prefix, run_id, peptide_type, date_time_str, suffix, ext)
    else:
        file_name = '{0}_{1}_{2}_clf_{3}.txt'.format(prefix, run_id, date_time_str, suffix)
    result_file = path.join(file_dir, file_name)
    while os.path.isfile(result_file):
        date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
        if peptide_type != '':
            file_name = '{0}_{1}_{2}_{3}_clf_{4}.{5}'. \
                format(prefix, run_id, peptide_type, date_time_str, suffix, ext)
        else:
            file_name = '{0}_{1}_{2}_clf_{3}.txt'.format(prefix, run_id, date_time_str, suffix)
        result_file = path.join(file_dir, file_name)

    return result_file




