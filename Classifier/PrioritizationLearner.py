import warnings
import pandas as pd
from sklearn.exceptions import UndefinedMetricWarning
import pickle
from Classifier.OptimizationParams import *
from sklearn.model_selection import StratifiedKFold
from scipy.stats import rankdata
from hyperopt import hp, fmin, tpe, rand, STATUS_OK, Trials
import time

from DataWrangling.DataLoader import *


warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)


class PrioritizationLearner:

    def __init__(self, classifier_tag, scorer_name, optimization_params, verbose=1, nr_iter=100, nr_cv=5,
                 nr_classifiers=1, alpha=0.005, shuffle=False, nr_epochs=150, batch_size=32, patience=10,):
        """ Performs GridSearchCV to train best classifier for prioritization"""

        self.classifier_tag = classifier_tag
        self.optimization_params = optimization_params
        self.classifier = self.optimization_params.get_base_classifier(self.classifier_tag)
        self.classifier_scorer = None
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

    def optimize_classifier(self, data, X, y, report_file=None):

        self.classifier_scorer = self.optimization_params.get_scorer(self.scorer_name, data)
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

            if self.classifier_tag == 'TabNet':
                X = X.to_numpy()

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

            print("Hyperopt: Score={0:.3f}, Time={1:f}, Params={2:s}".
                  format(((1 - objective.best_loss) * 100), elapsed_time_hopt, str(objective.best_params)))

            if report_file:
                report_file.write("Hyperopt: Score={0:.3f}; Time={1:f}; Params={2:s}\n".
                                  format(((1 - objective.best_loss) * 100), elapsed_time_hopt,
                                         str(objective.best_params)))

            self.fit_classifier(X, y, classifier=objective.best_classifier)

            return best, objective.best_classifier, objective.best_loss, objective.best_params

    def test_classifier(self, classifier, patient, data, X, y, max_rank=20, report_file=None, sort_columns=[]):
        self.classifier_scorer = self.optimization_params.get_scorer(self.scorer_name, data)

        if self.verbose > 1 and self.write_header:
            print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                  "CD8_ranks\tCD8_peptide_idx\tCD8_mut_seqs\tCD8_genes".format(max_rank))

        if report_file and os.path.getsize(report_file.name) == 0:
            report_file.write("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                              "CD8_ranks\tCD8_peptide_idx\tCD8_mut_seqs\tCD8_genes\n".format(max_rank))

        self.write_header = False

        if self.classifier_tag in ['LR', 'SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
            y_pred = classifier.predict_proba(X)[:, 1]
        else:
            y_pred = np.array(classifier.decision_function(X))
            # y_pred = np.array(classifier.predict(X), dtype=float)
            # y_pred = y_pred.flatten()
            # y_pred = np.array(classifier.predict(X))

        X_r = X.copy()
        X_r['ML_pred'] = y_pred
        X_r['response'] = y
        X_r.loc[:, 'gene'] = data.loc[:, 'gene']
        X_r.loc[:, 'mutant_seq'] = data.loc[:, 'mutant_seq']
        X_r.loc[:, 'peptide_id'] = data.loc[:, 'peptide_id']
        for c in sort_columns:
            if Parameters().get_order_relation(c) == '<':
                X_r.loc[:, c] = -X_r.loc[:, c]
        sort_columns = ['ML_pred'] + sort_columns
        X_r = X_r.sort_values(by=sort_columns, ascending=False)

        r = np.where(X_r['response'] == 1)[0]
        nr_correct = sum(r < max_rank)
        nr_immuno = sum(y == 1)
        score = self.classifier_scorer._score_func(y, y_pred)
        sort_idx = np.argsort(r)
        ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
        peptide_ids = X_r.loc[X_r['response'] == 1, 'peptide_id'].to_numpy()
        peptide_id_str = ",".join(["{0}".format(s) for s in peptide_ids[sort_idx]])
        mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
        mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
        genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
        gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

        if self.verbose > 1:
            print("%s\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s\t%s" %
                  (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), score, ranks_str, peptide_id_str,
                   mut_seqs_str, gene_str))

        if report_file:
            report_file.write("%s\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s\t%s\n" %
                              (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), score, ranks_str,
                               peptide_id_str, mut_seqs_str, gene_str))

        return X_r['ML_pred'], X_r, nr_correct, nr_immuno, r, score

    def test_classifiers(self, classifiers, patient, X, y, max_rank=20):

        if self.verbose > 1 and self.write_header:
            print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tCD8_ranks".format(max_rank))
            self.write_header = False

        y_pred_avg = np.zeros(X.shape[0])
        for classifier in classifiers:
            if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'LR', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
                y_pred = classifier.predict_proba(X)[:, 1]
            else:
                y_pred = np.array(classifier.decision_function(X))
#                y_pred = np.array(classifier.predict(X), dtype=float)
                y_pred = y_pred.flatten()
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

    def get_top_n_mutation_ids(self, classifier, data, X, max_rank=100):
        if self.classifier_tag in ['LR', 'SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
            y_pred = classifier.predict_proba(X)[:, 1]
        else:
            # y_pred = np.array(classifier.predict(X), dtype=float)
            # y_pred = y_pred.flatten()
            y_pred = np.array(classifier.predict(X))

        mutant_id = data.apply(self.convert_peptide_id_long, axis=1)
        df = pd.DataFrame({'mutant_id': mutant_id, 'prediction_score': y_pred})
        df.sort_values(by=['prediction_score'], ascending=False, ignore_index=True, inplace=True)

        max_rank = min(max_rank, X.shape[0])
        return df.loc[df.index[0:max_rank], 'mutant_id'].to_numpy()

    def convert_peptide_id_long(self, row):
        fields = row['peptide_id'].split('|')
        return fields[0]+":"+fields[2]

    def convert_peptide_id_short(self, row):
        fields = row['peptide_id'].split('|')
        patient = fields[0].split('-')[0]
        return patient+":"+fields[2]

    def get_mutation_score_short(self, row, df_long):
        idx = np.where(df_long['mutant_id'] == row['mut_seqid'])
        if len(idx[0]) > 0:
            return df_long.loc[df_long.index[idx[0][0]], 'mutation_score']
        else:
            return np.nan

    def get_mutation_rank_short(self, row, df_long):
        idx = np.where(df_long['mutant_id'] == row['mut_seqid'])
        if len(idx[0]) > 0:
            return df_long.loc[df_long.index[idx[0][0]], 'mutation_rank']
        else:
            return np.nan

    def get_peptide_score_long(self, row, df_short, peptide_scores_long, max_rank_short):
        idx = df_short['mut_seqid'] == row['mutant_id']
        if sum(idx) > 0:
            scores = df_short.loc[idx, 'peptide_score'].sort_values(ascending=False).head(max_rank_short)
            peptide_scores_long.append(scores.reset_index(drop=True))
        else:
            peptide_scores_long.append(pd.Series(np.full(max_rank_short, np.nan)))

    def add_long_prediction_to_short(self, data_long, data_short, x_short, y_short):
        assert 'mutation_score' in data_long.columns, "No mutation_score in data"
        assert 'mutation_rank' in data_long.columns, "No mutation_rank in data"

        if 'mutant_id' not in data_long.columns:
            mutant_id = data_long.apply(self.convert_peptide_id_long, axis=1)
            data_long.loc[:, 'mutant_id'] = mutant_id

        mutations_scores_short = data_short.apply(self.get_mutation_score_short, args=(data_long,), axis=1)
        mutations_ranks_short = data_short.apply(self.get_mutation_rank_short, args=(data_long,), axis=1)

        idx = np.isfinite(mutations_scores_short)
        data_short = data_short[idx]
        x_short = x_short[idx]
        y_short = y_short[idx]
        x_short.loc[:, 'mutation_score'] = mutations_scores_short[idx]
        data_short.loc[:, 'mutation_score'] = mutations_scores_short[idx]
        x_short.loc[:, 'mutation_rank'] = mutations_ranks_short[idx]
        data_short.loc[:, 'mutation_rank'] = mutations_ranks_short[idx]

        return data_short, x_short, y_short

    def add_short_prediction_to_long(self, data_short, max_rank_short, data_long, x_long, y_long):
        assert 'peptide_score' in data_short.columns, "No mutation_score in data"
        assert 'rank_in_mutation' in data_short.columns, "No mutation_rank in data"

        df_short = data_short.loc[:, ['mut_seqid', 'peptide_score']]
        if 'mutant_id' not in data_long.columns:
            mutant_id_long = data_long.apply(self.convert_peptide_id_long, axis=1)
            data_long.loc[:, 'mutant_id'] = mutant_id_long

        pred_long = \
            pd.DataFrame({'mutant_id': data_long['mutant_id'], 'mutation_score': data_long['mutation_score']},
                         index=data_long.index.copy())

        peptide_scores_long = []
        pred_long.apply(self.get_peptide_score_long, args=(df_short, peptide_scores_long, max_rank_short), axis=1)
        df = pd.concat(peptide_scores_long, axis=1, ignore_index=True).transpose()
        df.reindex_like(data_long)
        names = []
        for i in range(max_rank_short):
            names.append(f"peptide_score_{i}")
        df.columns = names
        idx = df["peptide_score_{0}".format(max_rank_short-1)].notna()
        data_long = data_long.merge(df, how='left', left_index=True, right_index=True)
        data_long = data_long.loc[idx, :]
        x_long = x_long.merge(df, how='left', left_index=True, right_index=True)
        x_long = x_long.loc[idx, :]
        y_long = y_long[idx]

        return data_long, x_long, y_long

    def fit_classifier(self, X, y, classifier=None, params=None):

        assert classifier is not None or params is not None

        if params is None:
            clf = classifier
        else:
            clf = self.optimization_params.get_classifier(self.classifier_tag, params)

        if self.classifier_tag == 'TabNet':
            stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
            for train_index, test_index in stratifiedKFold.split(X, y):
                X_train, X_test = self.X.iloc[train_index, :], self.X.iloc[test_index, :]
                y_train, y_test = self.y[train_index], self.y[test_index]
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
                X_train, X_test = self.X.iloc[train_index, :], self.X.iloc[test_index, :]
                y_train, y_test = self.y[train_index], self.y[test_index]
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
                X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
                y_train, y_test = y[train_index], y[test_index]
                break
            clf.fit(X_train, y_train, cat_features=self.optimization_params.cat_idx, eval_set=(X_test, y_test),
                    logging_level='Verbose', use_best_model=True, early_stopping_rounds=40, plot=False)

        elif self.classifier_tag == 'XGBoost':
            stratifiedKFold = StratifiedKFold(n_splits=self.nr_cv, shuffle=self.shuffle)
            for train_index, test_index in stratifiedKFold.split(X, y):
                X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
                y_train, y_test = y[train_index], y[test_index]
                break

            clf.fit(X_train, y_train, eval_set=[(X_test, y_test)], early_stopping_rounds=40, verbose=10)

        else:
            clf.fit(X, y)

        return clf

    def save(self, classifier_file):
        PrioritizationLearner.save(self.classifier_tag, self.classifier, classifier_file)

    def get_optimization_params(self):
        return self.optimization_params

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

