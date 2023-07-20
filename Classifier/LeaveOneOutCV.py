"""
Test classifiers by leave-one-out cross validation
"""

import argparse
import re
import glob
from multiprocessing import Pool
import pandas as pd

from catboost import CatBoostClassifier

from Classifier.ClassifierManager import ClassifierManager
from Classifier.OptimizationParams import OptimizationParams
from Utils.GlobalParameters import GlobalParameters
from Utils.Util_fct import *
from Utils.DataManager import DataManager
from sklearn.ensemble import VotingClassifier

parser = argparse.ArgumentParser(description='Test classifier models on train data with leave-one-out CV')
parser.add_argument('-c', '--classifier', type=str, default='LR', choices=GlobalParameters.classifiers,
                    help='classifier to use')
parser.add_argument('-ds', '--dataset', type=str, choices=GlobalParameters.datasets,
                    help='dataset used for training (should be same as used for encoding cat. features)')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory for classifier model files')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-tag', '--run_tag', type=str, help='Tag used in output file')


def get_clf_mgr(peptide_type: str, dataset_enc: str, classifier_name: str, x: pd.DataFrame, y: list):
    if args.peptide_type == 'neopep':
        class_ratio = sum(y == 1)/sum(y == 0)
    else:
        class_ratio = None

    alpha = GlobalParameters.neopep_alpha if peptide_type == 'neopep' else GlobalParameters.mutation_alpha
    features = GlobalParameters.ml_features_neopep \
        if peptide_type == 'neopep' else GlobalParameters.ml_features_mutation
    optimizationParams = \
        OptimizationParams(alpha=alpha,
                           cat_idx=DataManager.get_categorical_feature_idx(peptide_type, x),
                           cat_dims=DataManager.get_category_cnts(dataset_enc, peptide_type, x),
                           input_shape=[len(features)], class_ratio=class_ratio)

    return ClassifierManager(classifier_name, 'sum_exp_rank', optimizationParams, verbose=0)


def run_training(patient):
    global args, data, X, y

    if not DataManager.has_immunogenic_peptides(args.peptide_type, patient):
        return

    loo_result_file = os.path.join(GlobalParameters.classifier_result_dir, args.sub_dir,
                                   'Leave_one_out_{0}_test.txt'.format(args.dataset))

    with open(loo_result_file, mode='w') as result_file:
        idx_train = data['patients'] == patient
        data_train = data.loc[idx_train, :]
        X_train = X.loc[idx_train, :]
        y_train = y[idx_train]

        data_test, X_test, y_test = \
            DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml', patient=patient,
                                              sample=False)
        # hyperopt loop
        cvres, best_classifier, best_score, best_params = \
            clf_mgr.optimize_classifier(data_train, X_train, y_train)

        clf_mgr.test_classifier(args.classifier, best_classifier, args.peptide_type, patient, data_test, X_test,
                                y_test, report_file=result_file)


if __name__ == "__main__":

    args = parser.parse_args()

    data, X, y = \
        DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml',
                                          response_types=['CD8', 'negative'],
                                          dataset=args.dataset_train, sample=args.peptide_type == 'neopep')
    clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, args.classifier, X, y)

    with Pool(processes=len(data.patient.unique())) as pool:
        pool.map(run_training, data.patient.unique())

