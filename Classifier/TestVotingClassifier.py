import argparse
import re
import glob

from DataWrangling.DataTransformer import *
from Classifier.ClassifierManager import *
from Utils.Util_fct import *

parser = argparse.ArgumentParser(description='Test voting classifier on test data')
parser.add_argument('-d', '--sub_dir', default='', type=str, help='Subdirectory holding classifier model files')
parser.add_argument('-pt', '--peptide_type', type=str, choices=GlobalParameters.peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-c1', '--classifier1_file_re', type=str, nargs='+', help='classifier files for first batch')
parser.add_argument('-c2', '--classifier2_file_re', type=str, nargs='+', help='classifier files for second batch')
parser.add_argument('-w', '--weight', type=float, default=0.5, help='weight of second batch classifier')
parser.add_argument('-tr', '--dataset_train', type=str, choices=GlobalParameters.datasets_encoding,
                    help='dataset used for encoding')
parser.add_argument('-te', '--datasets_test', type=str, action='append', choices=GlobalParameters.datasets,
                    help='patient ids for test set')


if __name__ == "__main__":

    args = parser.parse_args()

    def get_clf_mgr(peptide_type: str, dataset_enc: str, classifier_name: str, x: pd.DataFrame):
        alpha = GlobalParameters.neopep_alpha if peptide_type == 'neopep' else GlobalParameters.mutation_alpha
        optimizationParams = \
            OptimizationParams(alpha=alpha,
                               cat_idx=DataManager.get_categorical_feature_idx(peptide_type, x),
                               cat_dims=DataManager.get_category_cnts(dataset_enc, peptide_type, x))

        return ClassifierManager(classifier_name, 'sum_exp_rank', optimizationParams, verbose=0)


    classifier1_model_files = []
    for wc in args.classifier1_file_re:
        classifier1_model_files = \
            classifier1_model_files + glob.glob(os.path.join(GlobalParameters.classifier_model_dir, args.sub_dir, wc))

    classifier2_model_files = []
    for wc in args.classifier2_file_re:
        classifier2_model_files = \
            classifier2_model_files + glob.glob(os.path.join(GlobalParameters.classifier_model_dir, args.sub_dir, wc))

    vc_result_file = os.path.join(GlobalParameters.classifier_result_dir, args.sub_dir,
                                  'Voting_classifier_{0:.2f}_test.txt'.format(args.weight))
    open(vc_result_file, mode='w').close()
    tesla_score_file = os.path.join(GlobalParameters.classifier_result_dir, args.sub_dir,
                                    'Voting_classifier_{0:.2f}_tesla_scores.txt'.format(args.weight))
    with open(tesla_score_file, mode='w') as score_file:
        score_file.write("Dataset\tPatient\tTTIF\tFR\tAUPRC\n")

    voting_clfs = []

    for ds in args.datasets_test:
        data, X, y = DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml',
                                                       dataset=ds, sample=False)

        for patient in data.patient.unique():
            if not DataManager.has_immunogenic_peptides(args.peptide_type, patient):
                continue

            data_test, X_test, y_test = \
                DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml', patient=patient,
                                                  sample=False)

            weights = []
            for model_file in classifier1_model_files:
                clf_name = os.path.basename(model_file).split("_")[0]
                clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, clf_name, X_test)
                classifier = clf_mgr.load_classifier(clf_name, clf_mgr.get_optimization_params(), model_file)
                voting_clfs.append((clf_name, classifier))
                weights.append(1.0-args.weight)

            for model_file in classifier2_model_files:
                clf_name = os.path.basename(model_file).split("_")[0]
                clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, clf_name, X_test)
                classifier = clf_mgr.load_classifier(clf_name, clf_mgr.get_optimization_params(), model_file)
                voting_clfs.append((clf_name, classifier))
                weights.append(args.weight)

            with open(vc_result_file, mode='a') as result_file:
                y_pred_sorted, X_sorted, nr_correct20, nr_tested20, nr_correct50, nr_tested50, nr_correct100, \
                nr_tested100, nr_immuno, r, score = \
                    clf_mgr.test_voting_classifier(voting_clfs, weights,  args.peptide_type, patient,
                                                   data_test, X_test, y_test, report_file=result_file)

                with open(tesla_score_file, mode='a') as score_file:
                    ClassifierManager.write_tesla_scores(patient, ds, nr_correct20, nr_tested20, nr_correct100,
                                                         nr_immuno, X_sorted, y_pred_sorted, score_file)

