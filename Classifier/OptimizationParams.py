import os
import numpy as np
import pandas as pd
import sklearn
from sklearn.svm import SVC
from sklearn.metrics import make_scorer
from sklearn.linear_model import LogisticRegression
from catboost import CatBoostClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from hyperopt import hp
from hyperopt.pyll import scope

from Utils.GlobalParameters import GlobalParameters


class OptimizationParams:

    def __init__(self, alpha: float = 0.05, cat_idx: list = None, class_ratio: float = None):
        """
        Class handling hyperparameters for Hyperopt loop

        Args:
            alpha (float): alpha value ro rank_score calculation
            cat_idx (array): list of indexes of columns with categorical features (only used in CatBoost)
            class_ratio (float): #immunogenic neo-peptides/#non-immunogenic neo-peptides
        """
        self.alpha = alpha
        self.cat_idx = cat_idx
        self.class_ratio = class_ratio  # (nr immunogenic peptides)/(nr non-immunogenic peptides)
        return

    def get_class_weights(self):
        if not self.class_ratio or self.class_ratio == 0 or self.class_ratio >= 1:
            return 'balanced'
        else:
            cws = []
            v = int(2.0/self.class_ratio)
            for cw in range(1, int(2.0/self.class_ratio), round(v/20)):
                cws.append({1: cw})
            return hp.choice('class_weight', cws)

    def get_xgb_pos_weights(self):
        if not self.class_ratio or self.class_ratio == 0 or self.class_ratio >= 1:
            return 1
        else:
            v = int(2.0/self.class_ratio)
            return hp.choice('scale_pos_weight', range(1, v, round(v/20)))

    def get_param_space(self, classifier_tag: str):
        """
        Defines parameter space searched durich Hyperopt loop for each classifier

        Args:
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
        Returns:
            parameter_space (dict): dictionary with parameter space for each classifier hyperparameter used in the
            Hyperopt loop
        """

        if classifier_tag == "SVM":
            parameter_space = {
                'C': hp.uniform('C', 0.005, 1.0),
                'gamma': hp.uniform('gamma', 0, 2),
                'class_weight': self.get_class_weights()
            }

        elif classifier_tag == "SVM-lin":
            parameter_space = {
                'C': hp.uniform('C', 0.005, 1.0),
                'class_weight': self.get_class_weights()
            }

        elif classifier_tag == "LR":

            parameter_space = {
                'penalty': hp.choice('penalty', ['l1', 'l2']),
                'C': hp.uniform('C', 0.0, 5.0),
                'class_weight': self.get_class_weights()
            }
        elif classifier_tag == "XGBoost":

            parameter_space = {
                'booster': hp.choice('booster', ['gbtree', 'gblinear']),
                'max_depth': hp.choice('max_depth', [3, 4, 5, 7, 9]),
                'min_child_weight': hp.choice('min_child_weight', np.round(np.arange(0.0,  0.2, 0.01), 5)),
                'learning_rate': hp.loguniform('learning_rate', np.log(0.01), np.log(0.1)),
                'subsample': hp.uniform('subsample', 0.3, 1.0),
                'colsample_bylevel': hp.uniform('colsample_bylevel', 0.4, 1.0),
                'colsample_bytree': hp.uniform('colsample_bytree', 0.4, 1.0),
                'n_estimators': hp.choice('n_estimators', np.arange(50, 1500, 50)),
                'reg_alpha': hp.loguniform('reg_alpha', np.log(0.0001), np.log(1)),
                'gamma':  hp.uniform('gamma', 0.0, 10.0),
                'scale_pos_weight': self.get_xgb_pos_weights(),
            }

        elif classifier_tag == "CatBoost":

            parameter_space = {
                'iterations': hp.choice('iterations', np.round(np.arange(100, 1500, 100))),
                'auto_class_weights': hp.choice('auto_class_weights', ['None', 'Balanced']),
                'subsample': hp.uniform('subsample', 0.3, 1.0),
                'random_strength': scope.int(hp.quniform('random_strength', 1, 20, 1)),
                'learning_rate': hp.loguniform('learning_rate', np.log(0.0001), np.log(1)),
                'l2_leaf_reg': hp.loguniform('l2_leaf_reg', np.log(1), np.log(10)),
                'leaf_estimation_iterations': scope.int(hp.quniform('leaf_estimation_iterations', 1, 20, 1)),
                'depth': scope.int(hp.quniform('depth', 5, 10, 1)),
                'bagging_temperature': hp.uniform('bagging_temperature', 0.0, 1.0)
            }

        return parameter_space

    def get_classifier(self, classifier_tag: str, params: dict):
        """
        Creates classifier object with hyperparameters corresponding to params

        Args:
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
            params (dict): dictionary with parameter values for classifier specified in classifier_tag
        Returns:
            classifier object (Object): classifier object with hyperparameters corresponding to params
        """

        if classifier_tag == "SVM":
            return SVC(probability=True, kernel='rbf', C=params['C'], gamma=params['gamma'],
                       class_weight=params['class_weight'])

        elif classifier_tag == "SVM-lin":
            return SVC(probability=True, kernel='linear', C=params['C'], class_weight=params['class_weight'])

        elif classifier_tag == "LR":
            return LogisticRegression(solver='saga',  penalty=params['penalty'], C=params['C'],
                                      class_weight=params['class_weight'])

        elif classifier_tag == "CatBoost":
            return CatBoostClassifier(
                loss_function='Logloss',
                iterations=params['iterations'],
                subsample=params['subsample'],
                random_strength=params['random_strength'],
                learning_rate=params['learning_rate'],
                l2_leaf_reg=params['l2_leaf_reg'],
                leaf_estimation_iterations=params['leaf_estimation_iterations'],
                depth=params['depth'],
                bagging_temperature=params['bagging_temperature'],
                random_seed=42,
                use_best_model=False,
                cat_features=self.cat_idx,
                auto_class_weights=params['auto_class_weights'],
                silent=True
            )

        elif classifier_tag == "XGBoost":
            return XGBClassifier(
                enable_categorical=True,
                max_depth=params['max_depth'],
                learning_rate=params['learning_rate'],
                n_estimators=params['n_estimators'],
                eval_metric='logloss',
                verbosity=0,
                silent=None,
                objective='binary:logistic',
                booster=params['booster'],
                tree_method='hist',
                n_jobs=int(os.cpu_count()/GlobalParameters.nr_hyperopt_rep),
                nthread=None,
                gamma=0,
                min_child_weight=params['min_child_weight'],
                max_delta_step=0,
                subsample=params['subsample'],
                colsample_bytree=params['colsample_bytree'],
                colsample_bylevel=params['colsample_bylevel'],
                colsample_bynode=1,
                reg_alpha=params['reg_alpha'],
                reg_lambda=1,
                scale_pos_weight=params['scale_pos_weight'],
                base_score=0.5,
                random_state=0,
                seed=None)

    def get_base_classifier(self, classifier_tag):
        """
        Creates classifier object with default hyperparameters

        Args:
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
        Returns:
            classifier object (Object): classifier object with default hyperparameters
        """

        if classifier_tag == "SVM":
            return SVC(probability=True, kernel='rbf')

        elif classifier_tag == "SVM-lin":
            return SVC(probability=True, kernel='linear')

        elif classifier_tag == "LR":
            return LogisticRegression(solver='saga')

        elif classifier_tag == "CatBoost":
            cb = CatBoostClassifier(
                loss_function='Logloss',
                iterations=10,
                learning_rate=0.01,
                random_seed=42,
                logging_level='Silent',
                use_best_model=False,
                cat_features=self.cat_idx,
                auto_class_weights='None',
                silent=True
            )
            return cb

        elif classifier_tag == "XGBoost":
            clf_xgb = XGBClassifier(
                enable_categorical=True,
                max_depth=8,
                learning_rate=0.1,
                n_estimators=1000,
                verbosity=0,
                silent=None,
                objective='binary:logistic',
                booster='gbtree',
                tree_method='hist',
                n_jobs=int(os.cpu_count()/GlobalParameters.nr_hyperopt_rep),
                nthread=None,
                gamma=0,
                min_child_weight=1,
                max_delta_step=0,
                subsample=0.7,
                colsample_bytree=1,
                colsample_bylevel=1,
                colsample_bynode=1,
                reg_alpha=0,
                reg_lambda=1,
                scale_pos_weight=1,
                base_score=0.5,
                random_state=0,
                seed=None)
            return clf_xgb

    def get_scorer(self, scorer_name: str, data: pd.DataFrame):
        """
        Defines a scorer object used as loss score in Hyperopt optimization

        Args:
            scorer_name (str): name of loss score function
            data (pd.DataFrame): dataframe (only used in scorer_name=='sum_exp_rank_pp')

        Returns:
            scorer (Callable): Callable object that returns a scalar score; greater is better.
        """

        if scorer_name == 'sum_exp_rank':
            return make_scorer(self.sum_rank_correct, needs_threshold=True)
        elif scorer_name == 'sum_exp_rank_pp':
            return make_scorer(self.sum_rank_correct_pp, patients=data['patient'].to_numpy(), needs_threshold=True)
        elif scorer_name == 'nr_correct_top100':
            return make_scorer(OptimizationParams.nr_correct_top100, needs_threshold=True)
        elif scorer_name == 'sum_prob_top100':
            return make_scorer(OptimizationParams.sum_prob_correct, needs_threshold=True)
        else:
            print('No scorer with name '+str(scorer_name)+' implemented. Abort')
            return None

    def sum_rank_correct(self, y_true: np.ndarray, y_pred: np.ndarray):
        """
        Rank_score optimization score used in the paper

        Args:
            y_true (np.ndarray): array with true immunogenicity indicators (0: non-immunogenic, 1: immunogenic)
            y_pred (np.ndarray): array with predicted probabilities that peptide is immunogenic

        Returns:
            rank_score (float): sum of rank_scores for all immunogenic peptides
        """
        idx = np.argsort(-y_pred)
        y_true = y_true[idx]
        r = np.where(y_true == 1)[0]

        return np.sum(np.exp(np.multiply(-self.alpha, r)))

    def sum_rank_correct_pp(self, y_true, y_pred, patients):
        score = 0.0
        for p in set(patients):
            idx = p == patients
            y_true_p = y_true[idx]
            y_pred_p = y_pred[idx]
            idx = np.argsort(-y_pred_p)
            y_true_p = y_true_p[idx]
            r = np.where(y_true_p == 1)[0]
            score += np.sum(np.exp(np.multiply(-self.alpha, r)))

        return score

    @staticmethod
    def nr_correct_top100(y_true, y_pred, max_rank=100):
        n = min(len(y_true, max_rank))
        idx = np.argsort(-y_pred)
        return np.sum(y_true[idx][:n] == 1)

    @staticmethod
    def sum_prob_correct(y_true, y_pred, max_rank=100):
        idx = np.argsort(-y_pred)
        n = min(len(y_true, max_rank))
        y_true = y_true[idx][:n]
        y_pred = y_pred[idx][:n]
        n = np.sum(y_true)
        if n == 0:
            s = 0
        else:
            s = np.sum(y_pred[y_true == 1])/n
        return s - np.sum(y_pred[y_true != 1])/(100-n)


class OptimizationObjective:

    def __init__(self, optimization_params: OptimizationParams, classifier_tag: str, x: pd.DataFrame, y: np.ndarray,
                 shuffle: bool = True, metric: str = 'sum_exp_rank'):
        """
        Class defining Hyperopt classifier optimization score function

        Args:
            optimization_params (OptimizationParams): OptimizationParams object
            classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
            x (pd.DataFrame): processed dataframe with rows and columns selected for ML
            y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
            shuffle (bool): shuffle dataframe before training
            metric (str): metric used during classifier parameter optimization (default sum of rank scores)
        """
        self.optimization_params = optimization_params
        self.classifier_tag = classifier_tag
        self.shuffle = shuffle
        self.stratifiedKFold = StratifiedKFold(n_splits=GlobalParameters.nr_hyperopt_cv, shuffle=self.shuffle)
        self.best_loss = np.Inf
        self.best_classifier = None
        self.best_params = None
        self.metric = metric
        self.X = x
        self.y = y

    def score(self, params):
        """
        Function defining loss score used in Hyperopt optimization

        Args:
            params (dict): dictionary with parameter values for classifier
        Returns:
            loss score (float):
        """
        classifier = self.optimization_params.get_classifier(self.classifier_tag, params)

        if self.classifier_tag in ['SVM', 'SVM-lin', 'LR', 'XGBoost']:
            # cross_val_score calls the metric function with arguments self.metric(self.y, classifier.predict(self.X))
            loss = 1 - \
                   cross_val_score(classifier, self.X, self.y, cv=self.stratifiedKFold, scoring=self.metric).mean()

            if loss < self.best_loss:
                self.best_loss = loss
                self.best_classifier = classifier
                self.best_params = params

            return loss

        elif self.classifier_tag == 'CatBoost':
            loss = 0
            for train_index, valid_index in self.stratifiedKFold.split(self.X, self.y):
                X_train, X_valid = self.X.iloc[train_index, :], self.X.iloc[valid_index, :]
                y_train, y_valid = self.y[train_index], self.y[valid_index]

                classifier.fit(X_train, y_train, use_best_model=True, eval_set=(X_valid, y_valid))
                res = classifier.get_best_score()
                loss += res['validation']['Logloss']

            loss = loss/GlobalParameters.nr_hyperopt_cv

            if loss < self.best_loss:
                self.best_loss = loss
                self.best_classifier = classifier
                self.best_params = params

            return loss


