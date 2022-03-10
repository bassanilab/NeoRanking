import numpy as np
import torch
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import make_scorer
from sklearn.linear_model import LogisticRegression
from tensorflow import keras
from catboost import CatBoostClassifier
from xgboost import XGBClassifier
from pytorch_tabnet.tab_model import TabNetClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import log_loss
from hyperopt import hp
from hyperopt.pyll import scope
import tensorflow as tf


class OptimizationObjective:

    def __init__(self, optimization_params, classifier_tag, X, y, nr_cv=5, shuffle=True, metric='binary_crossentropy'):
        self.optimization_params = optimization_params
        self.classifier_tag = classifier_tag
        self.shuffle = shuffle
        self.nr_cv = nr_cv
        self.stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
        self.best_loss = np.Inf
        self.best_classifier = None
        self.best_params = None
        self.metric = metric
        self.X = X
        self.y = y
        # if self.classifier_tag in ['LR']:
        #     self.y = np.array(list(map(lambda v: (v+1)/2, y)), dtype=float)
        # else:
        #     self.y = y

    def score(self, params):

        classifier = self.optimization_params.get_classifier(self.classifier_tag, params)

        if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'LR', 'NNN', 'XGBoost']:
            loss = 1-cross_val_score(classifier, self.X, self.y, cv=self.stratifiedKFold, scoring=self.metric, verbose=False).\
                mean()

            if loss < self.best_loss:
                self.best_loss = loss
                self.best_classifier = classifier
                self.best_params = params

            return loss

        elif self.classifier_tag == 'DNN':
            loss = 0
            for train_index, valid_index in self.stratifiedKFold.split(self.X, self.y):
                X_train, X_valid = self.X[train_index], self.X[valid_index]
                y_train, y_valid = self.y[train_index], self.y[valid_index]

                classifier.fit(X_train, y_train, epochs=params['nr_epochs'], batch_size=params['batch_size'])
                loss += classifier.evaluate(X_valid, y_valid)

            loss = loss/self.nr_cv

            if loss < self.best_loss:
                self.best_loss = loss
                self.best_classifier = classifier
                self.best_params = params

            return loss

        elif self.classifier_tag == 'CatBoost':
            loss = 0
            for train_index, valid_index in self.stratifiedKFold.split(self.X, self.y):
                X_train, X_valid = self.X[train_index], self.X[valid_index]
                y_train, y_valid = self.y[train_index], self.y[valid_index]

                classifier.fit(X_train, y_train, use_best_model=True, eval_set=(X_valid, y_valid))
                res = classifier.get_best_score()
                loss += res['validation']['Logloss']

            loss = loss/self.nr_cv

            if loss < self.best_loss:
                self.best_loss = loss
                self.best_classifier = classifier
                self.best_params = params

            return loss

        elif self.classifier_tag == 'TabNet':
            loss = 0
            for train_index, valid_index in self.stratifiedKFold.split(self.X, self.y):
                X_train, X_valid = self.X[train_index], self.X[valid_index]
                y_train, y_valid = self.y[train_index], self.y[valid_index]

                classifier.fit(
                    X_train=X_train, y_train=y_train,
                    eval_set=[(X_train, y_train), (X_valid, y_valid)],
                    eval_name=['train', 'valid'],
                    eval_metric=['logloss'],
                    max_epochs=1000, patience=20,
                    batch_size=1024, virtual_batch_size=128,
                    num_workers=0,
                    weights=1,
                    drop_last=False
                )
                y_pred = classifier.predict_proba(X_valid)[:, 1]
                loss += log_loss(y_valid, y_pred)

            loss = loss/self.nr_cv

            if loss < self.best_loss:
                self.best_loss = loss
                self.best_classifier = classifier
                self.best_params = params

            return loss


class OptimizationParams:

    def __init__(self, alpha=0.05, cat_features=None, cat_idx=None, cat_dims=None, input_shape=[10]):
        self.alpha = alpha
        self.input_shape = input_shape
        self.cat_features = cat_features
        self.cat_idx = cat_idx
        self.cat_dims = cat_dims
        return

    def get_param_space(self, classifier_tag):

        if classifier_tag == "SVM":
            cws = []
            for cw in range(1, 20, 2):
                cws.append({1: cw})

            parameter_space = {
                'C': hp.uniform('C', 0.005, 1.0),
                'gamma': hp.uniform('gamma', 0, 2),
                'class_weight': hp.choice('class_weight', cws)
            }

        elif classifier_tag == "SVM-lin":
            cws = []
            for cw in range(1, 50, 2):
                cws.append({1: cw})

            parameter_space = {
                'C': hp.uniform('C', 0.005, 1.0),
                'class_weight': hp.choice('class_weight', cws)
            }

        elif classifier_tag == "CART":
            cws = []
            for cw in range(5, 50, 2):
                cws.append({1: cw})

            parameter_space = {
                'max_depth': scope.int(hp.quniform('max_depth', 1, 50, 1)),
                'min_samples_split': scope.int(hp.quniform('min_samples_split', 2, 50, 1)),
                'min_samples_leaf': scope.int(hp.quniform('min_samples_leaf', 1, 100, 1)),
                'class_weight': hp.choice('class_weight', cws)
            }

        elif classifier_tag == "RF":
            cws = []
            for cw in range(5, 50, 2):
                cws.append({1: cw})

            parameter_space = {
                'max_depth': scope.int(hp.quniform('max_depth', 1, 50, 1)),
                'min_samples_split': scope.int(hp.quniform('min_samples_split', 2, 50, 1)),
                'n_estimators': scope.int(hp.quniform('n_estimators', 1, 250, 1)),
                'class_weight': hp.choice('class_weight', cws)
            }

        elif classifier_tag == "ADA":

            parameter_space = {
                'learning_rate': hp.loguniform('learning_rate', np.log(0.001), np.log(10)),
                'n_estimators': scope.int(hp.quniform('n_estimators', 1, 5, 1)),
            }

        elif classifier_tag == "NNN":

            parameter_space = {
                'weights': hp.choice('weights', ['uniform', 'distance']),
                'n_neighbors': scope.int(hp.quniform('n_neighbors', 1, 5, 1)),
                'metric': hp.choice('metric', ['minkowski', 'euclidean'])
            }

        elif classifier_tag == "LR":

            cws = []
            for cw in range(1, 50, 2):
                cws.append({1: cw})

            parameter_space = {
                'penalty': hp.choice('penalty', ['l1', 'l2']),
                'C': hp.uniform('C', 0.0, 2.0),
                'class_weight': hp.choice('class_weight', cws)
            }

            parameter_space = {
                'penalty': hp.choice('penalty', ['l1', 'l2']),
                'C': hp.uniform('C', 0.0, 5.0),
                'class_weight': "balanced"
            }
        elif classifier_tag == "DNN":

            parameter_space = {
                'nr_hidden': 2,
                'nr_neurons': 100,
                'activation': 'tanh',
                'learning_rate': 0.00021136243205598315,
                'batch_size': 32,
                'nr_epochs': 150,
                'scorer_name': 'sum_exp_rank',
                'input_shape': self.input_shape
            }
            # parameter_space = {
            #     'nr_hidden': hp.choice('nr_hidden', [1, 2, 3]),
            #     'nr_neurons': hp.choice('nr_neurons', [20, 30, 50, 100]),
            #     'activation': hp.choice('activation', ['relu', 'selu', 'tanh']),
            #     'learning_rate': hp.loguniform('learning_rate', np.log(0.0001), np.log(0.1)),
            #     'batch_size': 32,
            #     'nr_epochs': 150,
            #     'scorer_name': 'sum_exp_rank',
            #     'input_shape': self.input_shape
            # }

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
                'reg_alpha': hp.loguniform('reg_alpha', np.log(0.0001), np.log(1))
            }

        elif classifier_tag == "CatBoost":

            # parameter_space = {
            #     'iterations': [300, 500],
            #     'learning_rate': [1.0, 0.1, 0.01],
            #     'auto_class_weights': ['None', 'Balanced']
            # }
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

        elif classifier_tag == "TabNet":

            parameter_space = {
                'n_steps': hp.choice('n_steps', np.round(np.arange(1, 8, 1))),
                'n_d': hp.choice('n_d', np.round(np.arange(8, 64, 4))),
                'gamma': hp.uniform('gamma', 1.0, 2.0),
                'learning_rate': hp.loguniform('learning_rate', np.log(0.001), np.log(0.1)),
                'lambda_sparse': hp.loguniform('lambda_sparse', np.log(0.00001), np.log(0.1))
            }

        return parameter_space

    def get_classifier(self, classifier_tag, params):
        if classifier_tag == "SVM":
            return SVC(probability=True, kernel='rbf', C=params['C'], gamma=params['gamma'],
                       class_weight=params['class_weight'])

        elif classifier_tag == "SVM-lin":
            return SVC(probability=True, kernel='linear', C=params['C'], class_weight=params['class_weight'])

        elif classifier_tag == "CART":
            return DecisionTreeClassifier(min_samples_leaf=params['min_samples_leaf'],
                                          max_depth=params['max_depth'],
                                          min_samples_split=params['min_samples_split'],
                                          class_weight=params['class_weight'])
        elif classifier_tag == "RF":
            return RandomForestClassifier(n_estimators=params['n_estimators'],
                                          max_depth=params['max_depth'],
                                          min_samples_split=params['min_samples_split'],
                                          class_weight=params['class_weight'])

        elif classifier_tag == "ADA":
            return AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_depth=1),
                                      n_estimators=params['n_estimators'], learning_rate=params['learning_rate'])

        elif classifier_tag == "NNN":
            return KNeighborsClassifier(n_neighbors=params['n_neighbors'], weights=params['weights'],
                                        metric=params['metric'])

        elif classifier_tag == "LR":
            return LogisticRegression(solver='saga',  penalty=params['penalty'], C=params['C'],
                                      class_weight=params['class_weight'])

        elif classifier_tag == "DNN":

            if params['activation'] == 'relu':
                kernel_initializer = 'he_normal'
            elif params['activation'] == 'selu':
                kernel_initializer = 'lecun_normal'
            elif params['activation'] == 'tanh':
                kernel_initializer = 'glorot_uniform'
            else:
                kernel_initializer = 'glorot_uniform'

            #metrics = [tf.keras.metrics.AUC()]
            metrics = []

            model = keras.models.Sequential()
            model.add(keras.layers.InputLayer(input_shape=params['input_shape']))
            for layer in range(params['nr_hidden']):
                model.add(keras.layers.Dense(params['nr_neurons'], activation=params['activation'],
                                             kernel_initializer=kernel_initializer))

#            model.add(keras.layers.Dense(2, activation='softmax'))
            model.add(keras.layers.Dense(1, activation='linear'))
            optimizer = keras.optimizers.Adam(lr=params['learning_rate'], beta_1=0.9, beta_2=0.999)
            model.compile(loss=keras.losses.logcosh, optimizer=optimizer, metrics=metrics)

            return model

        elif classifier_tag == "CatBoost":
            return CatBoostClassifier(**params, cat_features=self.cat_features)

        elif classifier_tag == "XGBoost":
            return XGBClassifier(
                max_depth=params['max_depth'],
                learning_rate=params['learning_rate'],
                n_estimators=params['n_estimators'],
                eval_metric='logloss',
                verbosity=0,
                silent=None,
                objective='binary:logistic',
                booster=params['booster'],
                tree_method='gpu_hist',
                gpu_id=0,
                n_jobs=-1,
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
                scale_pos_weight=1,
                base_score=0.5,
                random_state=0,
                seed=None)

        elif classifier_tag == "TabNet":
            return TabNetClassifier(
                cat_idxs=self.cat_idx,
                cat_dims=self.cat_dims,
                cat_emb_dim=1,
                n_d=params['n_d'],
                n_a=params['n_d'],  # same value for n_d and n_a
                n_steps=params['n_steps'],
                gamma=params['gamma'],
                lambda_sparse=params['lambda_sparse'],
                optimizer_fn=torch.optim.Adam,
                optimizer_params=dict(lr=params['learning_rate']),
                scheduler_params={"step_size": 20, "gamma": 0.9},
                scheduler_fn=torch.optim.lr_scheduler.StepLR,
                mask_type='entmax'  # "sparsemax"
                )

    def get_base_classifier(self, classifier_tag):
        if classifier_tag == "SVM":
            return SVC(probability=True, kernel='rbf')

        elif classifier_tag == "SVM-lin":
            return SVC(probability=True, kernel='linear')

        elif classifier_tag == "CART":
            return DecisionTreeClassifier()

        elif classifier_tag == "RF":
            return RandomForestClassifier()

        elif classifier_tag == "ADA":
            return AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_depth=1))

        elif classifier_tag == "NNN":
            return KNeighborsClassifier()

        elif classifier_tag == "LR":
            return LogisticRegression(solver='saga')

        elif classifier_tag == "DNN":
            model = keras.models.Sequential()
            model.add(keras.layers.InputLayer(input_shape=[19]))
            model.add(keras.layers.Dense(50, activation='relu', kernel_initializer='he_normal'))
            model.add(keras.layers.Dense(1, activation='tanh'))
            optimizer = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999)
            model.compile(loss=keras.losses.MSE, optimizer=optimizer)

            return model

        elif classifier_tag == "CatBoost":
            cb = CatBoostClassifier(
                loss_function='Logloss',
                iterations=10,
                learning_rate=0.01,
                random_seed=42,
                logging_level='Silent',
                use_best_model=False,
                cat_features=self.cat_features,
                auto_class_weights='None'
            )
            return cb

        elif classifier_tag == "XGBoost":
            clf_xgb = XGBClassifier(
                max_depth=8,
                learning_rate=0.1,
                n_estimators=1000,
                verbosity=0,
                silent=None,
                objective='binary:logistic',
                booster='gbtree',
                tree_method='gpu_hist',
                gpu_id=0,
                n_jobs=-1,
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

        elif classifier_tag == "TabNet":
            return TabNetClassifier(
                cat_idxs=self.cat_idx,
                cat_dims=self.cat_dims,
                cat_emb_dim=1,
                optimizer_fn=torch.optim.Adam,
                optimizer_params=dict(lr=2e-2),
                scheduler_params={"step_size": 50,  "gamma": 0.9},
                scheduler_fn=torch.optim.lr_scheduler.StepLR,
                mask_type='entmax' # "sparsemax"
                )

    def get_scorer(self, scorer_name):
        if scorer_name == 'sum_exp_rank':
            return make_scorer(self.sum_rank_correct, needs_threshold=True)
        elif scorer_name == 'nr_correct_top100':
            return make_scorer(OptimizationParams.nr_correct_top100, needs_threshold=True)
        elif scorer_name == 'sum_prob_top100':
            return make_scorer(OptimizationParams.sum_prob_correct, needs_threshold=True)
        else:
            print('No scorer with name '+str(scorer_name)+' implemented. Abort')
            return None

    def score(self, scorer_name, y_true, y_pred):
        if scorer_name == 'sum_exp_rank':
            return self.sum_rank_correct(y_true, y_pred)
        elif scorer_name == 'nr_correct_top100':
            return OptimizationParams.nr_correct_top100(y_true, y_pred)
        elif scorer_name == 'sum_prob_top100':
            return OptimizationParams.sum_prob_correct(y_true, y_pred)
        else:
            print('No scorer with name '+str(scorer_name)+' implemented. Abort')
            return None

    def sum_rank_correct(self, y_true, y_pred):
        idx = np.argsort(-y_pred)
        y_true = y_true[idx]
        r = np.where(y_true == 1)[0]

        return np.sum(np.exp(np.multiply(-self.alpha, r)))

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
