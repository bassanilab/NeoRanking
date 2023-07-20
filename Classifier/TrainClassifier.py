import argparse
import time
import re
from multiprocessing import Pool

from DataWrangling.DataTransformer import *
from Classifier.ClassifierManager import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Add features to neodisc files')
parser.add_argument('-c', '--classifier', type=str, default='LR', choices=GlobalParameters.classifiers,
                    help='classifier to use')
parser.add_argument('-tr', '--dataset_train', type=str, choices=GlobalParameters.datasets,
                    help='dataset used for training (should be same as used for encoding cat. features)')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory for classifier model files')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-tag', '--run_tag', type=str, help='Tag used in output file')


def run_training(run_index):
    clf_model_file = get_classifier_file(args.classifier, args.sub_dir, args.run_tag, run_index, args.peptide_type)
    clf_param_file = re.sub("_clf\\.\\w+$", "_param.txt", clf_model_file)

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

    with open(clf_param_file, mode='w') as param_file:
        for arg in vars(args):
            param_file.write(f"{arg}={getattr(args, arg)}\n")
            print(f"{arg}={getattr(args, arg)}")

        data_train, X_train, y_train = \
            DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml',
                                              response_types=['CD8', 'negative'],
                                              dataset=args.dataset_train, sample=args.peptide_type == 'neopep')
        clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, args.classifier, X_train, y_train)

        # hyperopt loop
        cvres, best_classifier, best_score, best_params = \
            clf_mgr.optimize_classifier(data_train, X_train, y_train, param_file)

        ClassifierManager.save_classifier(args.classifier, best_classifier, clf_model_file)

        print('Classifier = {0:s}, run index = {1:d}\nBest training params: {2:s}\nsum_exp_rank = {3:.3f}\nSaved to {4:s}'.
              format(args.classifier, run_index, str(best_params), best_score, clf_model_file))

        param_file.write('Training dataset: {0}\n'.format(args.dataset_train))
        param_file.write('Saved to {0:s}\n'.format(clf_model_file))


if __name__ == "__main__":

    args = parser.parse_args()

    with Pool(processes=GlobalParameters.nr_hyperopt_rep) as pool:
        pool.map(run_training, range(GlobalParameters.nr_hyperopt_rep))



