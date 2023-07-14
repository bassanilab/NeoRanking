"""
Test classifier models on specified datasets
"""

import argparse
import re
import glob
import pandas as pd

from catboost import CatBoostClassifier

from Classifier.OptimizationParams import OptimizationParams
from Classifier.ClassifierManager import ClassifierManager
from Utils.GlobalParameters import GlobalParameters
from Utils.Util_fct import *
from Utils.DataManager import DataManager
from sklearn.ensemble import VotingClassifier

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory holding classifier model files')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-c', '--classifier_file_re', type=str, nargs='+', help='classifier files to use')
parser.add_argument('-tr', '--dataset_train', type=str, choices=GlobalParameters.datasets_encoding,
                    help='dataset used for encoding')
parser.add_argument('-te', '--datasets_test', type=str, action='append', choices=GlobalParameters.datasets,
                    help='patient ids for test set')

if __name__ == "__main__":

    args = parser.parse_args()

    for arg in vars(args):
        print(arg, getattr(args, arg))

    classifier_model_files = []
    classifier_result_files = []
    for wc in args.classifier_file_re:
        classifier_model_files = \
            classifier_model_files + glob.glob(os.path.join(GlobalParameters.classifier_model_dir, args.sub_dir, wc))
        for classifier_model_file in classifier_model_files:
            result_file_name = re.sub("_clf\\.\\w+$", "_test.txt", os.path.basename(classifier_model_file))
            result_file_name = os.path.join(GlobalParameters.classifier_result_dir, args.sub_dir, result_file_name)
            classifier_result_files.append(result_file_name)
            open(result_file_name, mode='w').close()


    def get_clf_mgr(peptide_type: str, dataset_enc: str, classifier_name: str, x: pd.DataFrame):
        alpha = GlobalParameters.neopep_alpha if peptide_type == 'neopep' else GlobalParameters.mutation_alpha
        features = GlobalParameters.ml_features_neopep \
            if peptide_type == 'neopep' else GlobalParameters.ml_features_mutation
        optimizationParams = \
            OptimizationParams(alpha=alpha,
                               cat_idx=DataManager.get_categorical_feature_idx(peptide_type, x),
                               cat_dims=DataManager.get_category_cnts(dataset_enc, peptide_type, x),
                               input_shape=[len(features)])

        return ClassifierManager(classifier_name, 'sum_exp_rank', optimizationParams, verbose=0)


    for ds in args.datasets_test:
        data, X, y = DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml',
                                                       dataset=ds, sample=False)

        for patient in data.patient.unique():
            if not DataManager.has_immunogenic_peptides(args.peptide_type, patient):
                continue

            data_test, X_test, y_test = \
                DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml', patient=patient,
                                                  sample=False)

            for model_file, result_file in zip(classifier_model_files, classifier_result_files):
                clf_name = os.path.basename(model_file).split("_")[0]
                clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, clf_name, X_test)
                classifier = clf_mgr.load_classifier(clf_name, clf_mgr.get_optimization_params(), model_file)

                with open(result_file, mode='a') as file:
                    y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
                        clf_mgr.test_classifier(clf_name, classifier, args.peptide_type, patient, data_test, X_test, y_test,
                                                report_file=file)

