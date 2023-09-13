import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
import pickle
from Classifier.OptimizationParams import *
from sklearn.model_selection import StratifiedKFold
from scipy.stats import rankdata
from hyperopt import hp, fmin, tpe, rand, STATUS_OK, Trials
import time
import os
from filelock import FileLock

from Utils.GlobalParameters import GlobalParameters
from DataWrangling.DataTransformer import DataTransformer
from Utils.DataManager import DataManager

warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)


class ClassifierManager:

    def __init__(self, classifier_tag: str, scorer_name: str, optimization_params: OptimizationParams, verbose: int = 1,
                 shuffle: bool = False):
        """
        Class for training and testing of classifiers
        Args:
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
            scorer_name (str): name of scoring function for hyperopt optimization ('sum_exp_rank',
                               'nr_correct_top100', 'sum_prob_top100')
            optimization_params (OptimizationParams): OptimizationParams object
            verbose (int): 0: (no prints), 1: some prints
            shuffle (bool): shuffle dataframe before training
        """

        self._classifier_tag: str = classifier_tag
        self._optimization_params: OptimizationParams = optimization_params
        # classifier with sklearn interface
        self._classifier = self._optimization_params.get_base_classifier(self._classifier_tag)
        # scorer returning single float value resulting from sklearn.metrics.make_scorer function
        self._classifier_scorer = None
        self._scorer_name: str = scorer_name
        self._verbose: int = verbose
        self._write_header: bool = True
        self._shuffle: bool = shuffle
        self._seed: int = 42
        return

    def optimize_classifier(self, data: pd.DataFrame, x: pd.DataFrame, y: np.ndarray, report_file: str = None) -> tuple:
        """
        Performs Hyperopt search for 'SVM', 'SVM-lin', 'LR', 'XGBoost', and 'CatBoost' classifiers

        Args:
            data (pd.DataFrame): unprocessed dataframe with rows selected for ML
            x (pd.DataFrame): processed dataframe with rows and columns selected for ML
            y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
            report_file (str): file name to write results to. If None nothing is written
        Returns:
            (best, best_classifier, best_loss, best_params)
                best: Hyperopt dictionary holding details of optimization
                best_classifier: the best classifier found in optimization
                best_loss: the best loss
                best_params: parameters of the best classifier
        """

        self._classifier_scorer = self._optimization_params.get_scorer(self._scorer_name, data)
        param_space = self._optimization_params.get_param_space(self._classifier_tag)

        if self._classifier_tag in ['SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost']:

            trials = Trials()  # Initialize an empty trials database for further saving/loading ran iteractions

            start = time.time()

            objective = OptimizationObjective(optimization_params=self._optimization_params,
                                              classifier_tag=self._classifier_tag, x=x, y=y,
                                              metric=self._classifier_scorer)

            best = fmin(objective.score,
                        space=param_space,
                        algo=tpe.suggest,
                        max_evals=GlobalParameters.nr_hyperopt_iter,
                        trials=trials,
                        rstate=np.random.RandomState(self._seed))

            elapsed_time_hopt = time.time() - start

            if self._verbose > 0:
                print("Hyperopt: Score={0:.3f}, Time={1:f}, Params={2:s}".
                      format(((1 - objective.best_loss) * 100), elapsed_time_hopt, str(objective.best_params)))

            if report_file is not None:
                report_file.write("Hyperopt: Score={0:.3f}; Time={1:f}; Params={2:s}\n".
                                  format(((1 - objective.best_loss) * 100), elapsed_time_hopt,
                                         str(objective.best_params)))
                report_file.flush()

            self.fit_classifier(x, y, classifier=objective.best_classifier)

            return best, objective.best_classifier, objective.best_loss, objective.best_params

    def test_classifier(self, classifier, peptide_type: str, patient: str, data: pd.DataFrame, x, y, max_rank=20,
                        report_file=None, sort_columns=[]) -> tuple:
        """
        Tests classifier on given patient's data

        Args:
            classifier (object): classifier to be tested. Classifier object must implement the sklearn predict_proba
                                 method
            data (pd.DataFrame): unprocessed dataframe with rows selected for ML
            x (pd.DataFrame): processed dataframe with rows and columns selected for ML
            y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
            max_rank (int): number of ranks taken into account to calculate topN counts. only needed dor report
                            output, but does not affect testing results
            report_file (str): file name to write results to. If None nothing is written
            sort_columns (list): additional sort columns to resolve ties in predict_proba sorting
        Returns:
            (predicted_proba, x_sorted_annot, nr_correct, nr_immuno, ranks, score)
                predicted_proba: predicted immunogenic probability
                x_sorted_annot: x sorted by predict_proba probability with additional columns:
                                ML_pred: predict_proba values
                                response: response (0/1), Corresponds to sorted vector y
                                mutant_seq: peptide with mutation
                                gene: mutated gene containing mutant_seq
                nr_correct: nr of immunogenic peptides in top max_rank ranks after sorting
                nr_immuno: nr immunogenic peptides for this patient
                ranks: ranks of the immunogenic peptide of this patients
                score: score for the ranking of this patient (as defined by scorer_name in __init__)
        """
        self._classifier_scorer = self._optimization_params.get_scorer(self._scorer_name, data)

        if self._verbose > 1 and self._write_header:
            print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                  "CD8_ranks\tCD8_mut_seqs\tCD8_genes".format(max_rank))

        if report_file is not None:
            lock = FileLock(report_file+".lock")
            with lock:
                if os.path.getsize(report_file) == 0:
                    with open(report_file, mode='w') as file:
                        file.write("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                                   "CD8_ranks\tCD8_mut_seqs\tCD8_genes\n".format(max_rank))

        self._write_header = False

        y_pred = classifier.predict_proba(x)[:, 1]

        X_r = x.copy()
        X_r['ML_pred'] = y_pred
        X_r['response'] = y
        X_r.loc[:, 'gene'] = data.loc[:, 'gene']
        X_r.loc[:, 'mutant_seq'] = data.loc[:, 'mutant_seq']
        for c in sort_columns:
            if peptide_type == 'neopep':
                if GlobalParameters.ml_feature_mv_neopep[c] == 'max':
                    X_r.loc[:, c] = -X_r.loc[:, c]
            elif peptide_type == 'mutation':
                if GlobalParameters.ml_feature_mv_mutation[c] == 'max':
                    X_r.loc[:, c] = -X_r.loc[:, c]
        sort_columns = ['ML_pred'] + sort_columns
        X_r = X_r.sort_values(by=sort_columns, ascending=False)

        r = np.where(X_r['response'] == 1)[0]
        nr_correct = sum(r < max_rank)
        nr_immuno = sum(y == 1)
        score = self._classifier_scorer._score_func(y, y_pred)
        sort_idx = np.argsort(r)
        ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
        mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
        mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
        genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
        gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

        if self._verbose > 1:
            print("%s\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s" %
                  (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), score, ranks_str,
                   mut_seqs_str, gene_str))

        if report_file is not None:
            lock = FileLock(report_file+".lock")
            with lock:
                with open(report_file, mode='a') as file:
                    file.write("{0:s}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\t{5:.5f}\t{6:s}\t{7:s}\t{8:s}\n".
                               format(patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), score,
                                      ranks_str, mut_seqs_str, gene_str))

        return X_r['ML_pred'], X_r, nr_correct, nr_immuno, r, score

    def test_voting_classifier(self, classifiers: list, weights: list, patient: str, data: pd.DataFrame,
                               x: pd.DataFrame, y: np.ndarray, report_file: str = None, sort_columns: list = []) \
            -> tuple:
        """
        Tests classifier on given patient's data

        Args:
            classifiers (list): classifiers to be combined into voting classifier and tested. Classifier objects must
                                implement the sklearn predict_proba method
            weights (list(float)): weights of each classifier in classifiers. weights must be positive.
            data (pd.DataFrame): unprocessed dataframe with rows selected for ML
            x (pd.DataFrame): processed dataframe with rows and columns selected for ML
            y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
            report_file (str): file name to write results to. If None nothing is written
            sort_columns (list): additional sort columns to resolve ties in predict_proba sorting
        Returns:
            (predicted_proba, x_sorted_annot, nr_correct20, nr_tested20, nr_correct50, nr_tested50, nr_correct100,
             nr_tested100, nr_immuno, ranks, score)
                predicted_proba: predicted immunogenic probability
                x_sorted_annot: x sorted by predict_proba probability with additional columns:
                                ML_pred: predict_proba values
                                response: response (0/1), Corresponds to sorted vector y
                                mutant_seq: peptide with mutation
                                gene: mutated gene containing mutant_seq
                nr_correct20: nr of immunogenic peptides in top 20 ranks after sorting
                nr_tested20: nr of tested peptides in top 20 ranks after sorting
                nr_correct50: nr of immunogenic peptides in top 50 ranks after sorting
                nr_tested50: nr of tested peptides in top 50 ranks after sorting
                nr_correct100: nr of immunogenic peptides in top 100 ranks after sorting
                nr_tested100: nr of tested peptides in top 100 ranks after sorting
                nr_immuno: nr immunogenic peptides for this patient
                ranks: ranks of the immunogenic peptide of this patients
                score: score for the ranking of this patient (as defined by scorer_name in __init__)
        """

        self._classifier_scorer = self._optimization_params.get_scorer(self._scorer_name, data)

        if self._verbose > 1 and self._write_header:
            print("Patient\tNr_correct_top20\tNr_tested_top20\tNr_correct_top50\tNr_tested_top50\t"
                  "Nr_correct_top100\tNr_tested_top100\tNr_immunogenic\tNr_peptides\tClf_score\t"
                  "CD8_ranks\tCD8_mut_seqs\tCD8_genes")

        if report_file and os.path.getsize(report_file.name) == 0:
            report_file.write("Patient\tNr_correct_top20\tNr_tested_top20\tNr_correct_top50\tNr_tested_top50\t"
                              "Nr_correct_top100\tNr_tested_top100\tNr_immunogenic\tNr_peptides\tClf_score\t"
                              "CD8_ranks\tCD8_mut_seqs\tCD8_genes\n")

        self._write_header = False

        y_pred = np.full(len(y), 0.0)
        for (w, clf) in zip(weights, classifiers):
            print(clf)
            y_pred = np.add(y_pred, np.array(clf[1].predict_proba(x)[:, 1]) * w)

        X_r = x.copy()
        X_r['ML_pred'] = y_pred
        X_r['response'] = y
        X_r['response_type'] = data['response_type']
        X_r.loc[:, 'gene'] = data.loc[:, 'gene']
        X_r.loc[:, 'mutant_seq'] = data.loc[:, 'mutant_seq']
        for c in sort_columns:
            if GlobalParameters().get_order_relation(c) == '<':
                X_r.loc[:, c] = -X_r.loc[:, c]
        sort_columns = ['ML_pred'] + sort_columns
        X_r = X_r.sort_values(by=sort_columns, ascending=False)

        r = np.where(X_r['response'] == 1)[0]
        rt = np.where(X_r['response_type'] == 'negative')[0]
        nr_correct20 = sum(r < 20)
        nr_tested20 = nr_correct20 + sum(rt < 20)
        nr_correct50 = sum(r < 50)
        nr_tested50 = nr_correct50 + sum(rt < 50)
        nr_correct100 = sum(r < 100)
        nr_tested100 = nr_correct100 + sum(rt < 100)
        nr_immuno = sum(y == 1)
        score = self._classifier_scorer._score_func(y, y_pred)
        sort_idx = np.argsort(r)
        ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
        mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
        mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
        genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
        gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

        if self._verbose > 1:
            print("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s" %
                  (patient, nr_correct20, nr_tested20, nr_correct50, nr_tested50, nr_correct100, nr_tested100,
                   nr_immuno, len(y), score, ranks_str, mut_seqs_str, gene_str))

        if report_file:
            report_file.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s\n" %
                              (patient, nr_correct20, nr_tested20, nr_correct50, nr_tested50,
                               nr_correct100, nr_tested100, nr_immuno, len(y), score, ranks_str,
                               mut_seqs_str, gene_str))
            report_file.flush()

        return X_r['ML_pred'], X_r, nr_correct20, nr_tested20, nr_correct50, nr_tested50, nr_correct100, nr_tested100,\
               nr_immuno, r, score

    def fit_classifier(self, x: pd.DataFrame, y: np.ndarray, classifier=None, params: dict = None) -> object:
        """
        Calls classifier.fit. If params is None, then classifier object is used. Otherwise, a classifier is
        constructed with hyperparameters defined in params
        Args:
            x (pd.DataFrame): processed dataframe with rows and columns selected for ML
            y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
            classifier (object): classifier to be fitted. Classifier object must implement the sklearn fit
                                 method
            params (dict): dictionary with classifiers hyperparameters
        Returns:
            fitted classifier object
        """

        assert classifier is not None or params is not None

        if params is None:
            clf = classifier
        else:
            clf = self._optimization_params.get_classifier(self._classifier_tag, params)

        if self._classifier_tag == 'CatBoost':
            clf.fit(x, y, plot=False)
        else:
            clf.fit(x, y)

        return clf

    def get_optimization_params(self) -> OptimizationParams:
        """
        Getter for self._optimization_params
        Returns:
            OptimizationParams object
        """
        return self._optimization_params

    @staticmethod
    def save_classifier(classifier_tag: str, classifier, classifier_file: str):
        """
        Saves classifier to file.
        Args:
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
            classifier (object): classifier to be fitted. Classifier object must implement the sklearn fit
                                 method
            classifier_file (str): file name of classifier model file
        """

        if classifier_tag in ['CatBoost', 'XGBoost']:
            classifier.save_model(classifier_file)
        else:
            pickle.dump(classifier, open(classifier_file, 'wb'))

    def load(self, classifier_file: str) -> object:
        """
        Loads classifier from file.
        Args:
            classifier_file (str): file name of classifier model file
        Returns:
            classifier object
        """
        self._classifier = ClassifierManager.load_classifier(
            self._classifier_tag, self._optimization_params, classifier_file)
        return self._classifier

    @staticmethod
    def load_classifier(classifier_tag, optimization_params, classifier_file) -> object:
        """
        Saves classifier to file.
        Args:
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
            optimization_params (OptimizationParams): OptimizationParams object
            classifier_file (str): file name of classifier model file
        Returns:
            classifier object
        """

        if classifier_tag in ['CatBoost', 'XGBoost']:
            classifier = optimization_params.get_base_classifier(classifier_tag)
            classifier.load_model(classifier_file)
        else:
            classifier = pickle.load(open(classifier_file, 'rb'))

        return classifier

    @staticmethod
    def write_tesla_scores(patient: str, dataset: str, nr_correct20: int, nr_tested20: int, nr_correct100: int,
                           nr_immuno: int, x_sorted: pd.DataFrame, y_pred_sorted: np.ndarray, report_file: str):
        """
        Calculates TESLA FP, TTIF, and AUPRC scores for a patient and writes them to report file

        Args:
            patient (str): patient id.
            dataset (bool): dataset id. if not provided all datasets are considered
            nr_correct20 (int): nr immunogenic peptides ranked in top 20 for this patient
            nr_tested20 (int): nr tested peptides ranked in top 20 for this patient
            nr_correct100 (int): nr immunogenic peptides ranked in top 100 for this patient
            nr_immuno (int): nr immunogenic peptides for this patient
            x_sorted (pd.DataFrame): ml data matrix sorted by predict_proba
            y_pred_sorted (np.ndarray): sorted response types (1 = immunogenic)
            report_file (str): report file
        """
        idx = x_sorted['response_type'] != 'not_tested'
        y_pred_tesla = y_pred_sorted[idx].to_numpy()
        y_tesla = x_sorted.loc[idx, 'response'].to_numpy()
        ttif = nr_correct20/nr_tested20
        fr = nr_correct100/nr_immuno
        precision, recall, _ = precision_recall_curve(y_tesla, y_pred_tesla)
        auprc = auc(recall, precision)
        report_file.write("{0}\t{1}\t{2:.3f}\t{3:.3f}\t{4:.3}\n".format(dataset, patient, ttif, fr, auprc))


