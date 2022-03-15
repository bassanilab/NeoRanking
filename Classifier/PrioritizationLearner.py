import warnings
from sklearn.exceptions import UndefinedMetricWarning
import pickle
from Classifier.OptimizationParams import *
from sklearn.model_selection import StratifiedKFold
from scipy.stats import rankdata
from hyperopt import hp, fmin, tpe, rand, STATUS_OK, Trials
import time


warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)


class PrioritizationLearner:

    def __init__(self, classifier_tag, scorer_name, optimization_params, verbose=1, nr_iter=100, nr_cv=5,
                 nr_classifiers=1, alpha=0.005, shuffle=False, nr_epochs=150, batch_size=32, patience=10):
        """ Performs GridSearchCV to train best classifier for prioritization"""

        self.classifier_tag = classifier_tag
        self.optimization_params = optimization_params
        self.classifier = self.optimization_params.get_base_classifier(self.classifier_tag)
        self.classifier_scorer = self.optimization_params.get_scorer(scorer_name)

        self.scorer_name = scorer_name
        self.verbose = verbose
        self.nr_iter = nr_iter
        self.nr_cv = nr_cv
        self.nr_classifiers = nr_classifiers
        self.write_header = True
        self.alpha = alpha
        self.shuffle = shuffle
        self.nr_epochs = nr_epochs
        self.batch_size = batch_size
        self.early_stopping_patience = patience
        self.seed = 42
        return

    def optimize_classifier(self, X, y):

        param_space = self.optimization_params.get_param_space(self.classifier_tag)

        if self.classifier_tag == '__CatBoost':
            rnd_search = self.classifier.randomized_search(param_space, n_iter=20, cv=3, X=X, y=y, train_size=0.8,
                                                           refit=True, shuffle=True, stratified=True, plot=False)
            return rnd_search['cv_results'], self.classifier, self.classifier.score(X, y), rnd_search['params']

        elif self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'LR', 'NNN', 'DNN', 'XGBoost', 'CatBoost',
                                     'TabNet']:

            n_iter_hopt = self.nr_iter
            trials = Trials()  # Initialize an empty trials database for further saving/loading ran iteractions

            start = time.time()

            objective = OptimizationObjective(optimization_params=self.optimization_params,
                                              classifier_tag=self.classifier_tag, X=X, y=y, nr_cv=self.nr_cv,
                                              metric=self.classifier_scorer)

            best = fmin(objective.score,
                        space=param_space,
                        algo=tpe.suggest,
                        max_evals=n_iter_hopt,
                        trials=trials,
                        rstate=np.random.RandomState(self.seed))

            elapsed_time_hopt = time.time() - start

            self.fit_classifier(X, y, classifier=objective.best_classifier)
            print("Hyperopt: Score={0:.3f}, Time={1:f}, Params={2:s}".
                  format(((1 - objective.best_loss) * 100), elapsed_time_hopt, str(objective.best_params)))

            return best, objective.best_classifier, objective.best_loss, objective.best_params

    def test_classifier(self, classifier, patient, X, y, max_rank=20, report_file=None):

        if self.verbose > 1 and self.write_header:
            print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tCD8_ranks".format(max_rank))

        if report_file and self.write_header:
            report_file.write("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tCD8_ranks\n".
                              format(max_rank))

        self.write_header = False

        if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
            y_pred = classifier.predict_proba(X)[:, 1]
        else:
            # y_pred = np.array(classifier.predict(X), dtype=float)
            # y_pred = y_pred.flatten()
            y_pred = np.array(classifier.decision_function(X))

        r = rankdata(-y_pred, method='average')[y == 1]
        mut_idx = np.arange(len(y))[y == 1]
        nr_correct = sum(r <= max_rank)
        nr_immuno = sum(y == 1)
        score = self.classifier_scorer._score_func(y, y_pred)

        if self.verbose > 1:
            sort_idx = np.argsort(r)
            ranks_str = ",".join(["{0:.0f}".format(np.floor(r)) for r in r[sort_idx]])
            print("%s\t%d\t%d\t%d\t%d\t%s\t%f" % (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y),
                                                  ranks_str, score))
        if report_file:
            report_file.write("%s\t%d\t%d\t%d\t%d\t%s\t%f\n" %
                              (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), ranks_str, score))

        return y_pred, nr_correct, nr_immuno, r, mut_idx, score

    def test_classifiers(self, classifiers, patient, X, y, max_rank=20):

        if self.verbose > 1 and self.write_header:
            print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tCD8_ranks".format(max_rank))
            self.write_header = False

        y_pred_avg = np.zeros(X.shape[0])
        for classifier in classifiers:
            if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'LR', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
                y_pred = classifier.predict_proba(X)[:, 1]
            else:
                y_pred = np.array(classifier.predict(X), dtype=float)
                y_pred = y_pred.flatten()
    #            y_pred = np.array(classifier.decision_function(X))
            y_pred = np.divide(y_pred, max(y_pred))
            y_pred_avg = np.add(y_pred_avg, y_pred)

        y_pred_avg = np.divide(y_pred_avg, len(classifiers))

        r = rankdata(-y_pred_avg, method='average')[y == 1]
        nr_correct = sum(r <= max_rank)
        nr_immuno = sum(y == 1)
        score = self.classifier_scorer._score_func(y, y_pred_avg)

        if self.verbose > 1:
            sort_idx = np.argsort(r)
            print("%s\t%d\t%d\t%d\t%d\t%s\t%f" % (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y),
                                                  str(r[sort_idx]), score))

        return y_pred, nr_correct, nr_immuno, r, score

    def get_top_n(self, classifier, patient, X, y, max_rank=20):
        if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
            y_pred = classifier.predict_proba(X)[:, 1]
        else:
            # y_pred = np.array(classifier.predict(X), dtype=float)
            # y_pred = y_pred.flatten()
            y_pred = np.array(classifier.decision_function(X))

        r = rankdata(-y_pred, method='average')[y == 1]
        mut_idx = np.arange(len(y))[y == 1]
        nr_correct = sum(r <= max_rank)
        nr_immuno = sum(y == 1)
        score = self.classifier_scorer._score_func(y, y_pred)

        return y_pred, nr_correct, nr_immuno, r, mut_idx, score

    def fit_classifier(self, X, y, classifier=None, params=None):

        assert classifier is not None or params is not None

        if params is None:
            clf = classifier
        else:
            clf = self.optimization_params.get_classifier(self.classifier_tag, params)

        if self.classifier_tag == 'TabNet':
            stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
            for train_index, test_index in stratifiedKFold.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                break
            max_epochs = 1000

            clf.fit(
                X_train=X_train, y_train=y_train,
                eval_set=[(X_train, y_train), (X_test, y_test)],
                eval_name=['train', 'valid'],
                eval_metric=['logloss'],
                max_epochs=max_epochs, patience=self.early_stopping_patience,
                batch_size=1024, virtual_batch_size=128,
                num_workers=0,
                weights=1,
                drop_last=False
            )

        if self.classifier_tag == 'DNN':
            stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
            for train_index, test_index in stratifiedKFold.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                break
            early_stopping_cb = keras.callbacks.EarlyStopping(patience=self.early_stopping_patience,
                                                              restore_best_weights=True)
            clf.fit(X_train, y_train, epochs=self.nr_epochs, batch_size=self.batch_size,
                    validation_data=(X_test, y_test), callbacks=[early_stopping_cb])
            # cl_weight_0 = sum(y)/len(y)
            # rnd_search.fit(X, y, epochs=10, batch_size=32, class_weight={0: cl_weight_0, 1: 1-cl_weight_0})
        elif self.classifier_tag == 'CatBoost':
            stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
            for train_index, test_index in stratifiedKFold.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                break
            clf.fit(X_train, y_train, cat_features=self.optimization_params.cat_features, eval_set=(X_test, y_test),
                    logging_level='Verbose', use_best_model=True, early_stopping_rounds=40, plot=False)

        elif self.classifier_tag == 'XGBoost':
            stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
            for train_index, test_index in stratifiedKFold.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                break

            clf.fit(X_train, y_train, eval_set=[(X_test, y_test)], early_stopping_rounds=40, verbose=10)

        else:
            clf.fit(X, y)

        return clf

    def save(self, classifier_file):
        PrioritizationLearner.save(self.classifier_tag, self.classifier, classifier_file)

    @staticmethod
    def save_classifier(classifier_tag, classifier, classifier_file):

        if classifier_tag == 'DNN':
            classifier.save(classifier_file)
        elif classifier_tag in ['CatBoost', 'XGBoost', 'TabNet']:
            classifier.save_model(classifier_file)
        else:
            pickle.dump(classifier, open(classifier_file, 'wb'))

    def load(self, classifier_file):
        self.classifier = PrioritizationLearner.load_classifier(
            self.classifier_tag, self.optimization_params, classifier_file)
        return self.classifier

    @staticmethod
    def load_classifier(classifier_tag, optimization_params, classifier_file):

        if classifier_tag == 'DNN':
            classifier = keras.models.load_model(classifier_file)
        elif classifier_tag in ['CatBoost', 'XGBoost', 'TabNet']:
            classifier = optimization_params.get_base_classifier(classifier_tag)
            classifier.load_model(classifier_file)
        else:
            classifier = pickle.load(open(classifier_file, 'rb'))

        return classifier

