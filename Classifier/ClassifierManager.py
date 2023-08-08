import warnings
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

    def __init__(self, classifier_tag, scorer_name, optimization_params, verbose=1, shuffle=False):

        self.classifier_tag = classifier_tag
        self.optimization_params = optimization_params
        self.classifier = self.optimization_params.get_base_classifier(self.classifier_tag)
        self.classifier_scorer = None
        self.scorer_name = scorer_name
        self.verbose = verbose
        self.write_header = True
        self.shuffle = shuffle
        self.seed = 42
        return

    def optimize_classifier(self, data, X, y, report_file=None):

        self.classifier_scorer = self.optimization_params.get_scorer(self.scorer_name, data)
        param_space = self.optimization_params.get_param_space(self.classifier_tag)

        if self.classifier_tag == '__CatBoost':
            rnd_search = self.classifier.randomized_search(param_space, n_iter=20, cv=3, X=X, y=y, train_size=0.8,
                                                           refit=True, shuffle=True, stratified=True, plot=False)
            return rnd_search['cv_results'], self.classifier, self.classifier.score(X, y), rnd_search['params']

        elif self.classifier_tag in ['SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost']:

            trials = Trials()  # Initialize an empty trials database for further saving/loading ran iteractions

            start = time.time()

            objective = OptimizationObjective(optimization_params=self.optimization_params,
                                              classifier_tag=self.classifier_tag, X=X, y=y,
                                              metric=self.classifier_scorer)

            best = fmin(objective.score,
                        space=param_space,
                        algo=tpe.suggest,
                        max_evals=GlobalParameters.nr_hyperopt_iter,
                        trials=trials,
                        rstate=np.random.RandomState(self.seed))

            elapsed_time_hopt = time.time() - start

            print("Hyperopt: Score={0:.3f}, Time={1:f}, Params={2:s}".
                  format(((1 - objective.best_loss) * 100), elapsed_time_hopt, str(objective.best_params)))

            if report_file is not None:
                report_file.write("Hyperopt: Score={0:.3f}; Time={1:f}; Params={2:s}\n".
                                  format(((1 - objective.best_loss) * 100), elapsed_time_hopt,
                                         str(objective.best_params)))
                report_file.flush()

            self.fit_classifier(X, y, classifier=objective.best_classifier)

            return best, objective.best_classifier, objective.best_loss, objective.best_params

    def test_classifier(self, classifier, peptide_type, patient, data, X, y, max_rank=20,
                        report_file=None, sort_columns=[]):
        self.classifier_scorer = self.optimization_params.get_scorer(self.scorer_name, data)

        if self.verbose > 1 and self.write_header:
            print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                  "CD8_ranks\tCD8_mut_seqs\tCD8_genes".format(max_rank))

        if report_file is not None:
            lock = FileLock(report_file+".lock")
            with lock:
                if os.path.getsize(report_file) == 0:
                    with open(report_file, mode='w') as file:
                        file.write("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                                   "CD8_ranks\tCD8_mut_seqs\tCD8_genes\n".format(max_rank))

        self.write_header = False

        y_pred = classifier.predict_proba(X)[:, 1]

        X_r = X.copy()
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
        score = self.classifier_scorer._score_func(y, y_pred)
        sort_idx = np.argsort(r)
        ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
        mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
        mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
        genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
        gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

        if self.verbose > 1:
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

    def test_voting_classifier(self, classifiers, weights, peptide_type, patient, data, X, y, report_file=None, sort_columns=[]):
        self.classifier_scorer = self.optimization_params.get_scorer(self.scorer_name, data)

        if self.verbose > 1 and self.write_header:
            print("Patient\tNr_correct_top20\tNr_tested_top20\tNr_correct_top50\tNr_tested_top50\t"
                  "Nr_correct_top100\tNr_tested_top100\tNr_immunogenic\tNr_peptides\tClf_score\t"
                  "CD8_ranks\tCD8_mut_seqs\tCD8_genes")

        if report_file and os.path.getsize(report_file.name) == 0:
            report_file.write("Patient\tNr_correct_top20\tNr_tested_top20\tNr_correct_top50\tNr_tested_top50\t"
                              "Nr_correct_top100\tNr_tested_top100\tNr_immunogenic\tNr_peptides\tClf_score\t"
                              "CD8_ranks\tCD8_mut_seqs\tCD8_genes\n")

        self.write_header = False

        y_pred = np.full(len(y), 0.0)
        for (w, clf) in zip(weights, classifiers):
            print(clf)
            y_pred = np.add(y_pred, np.array(clf[1].predict_proba(X)[:, 1])*w)

        X_r = X.copy()
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
        score = self.classifier_scorer._score_func(y, y_pred)
        sort_idx = np.argsort(r)
        ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
        mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
        mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
        genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
        gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

        if self.verbose > 1:
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

    def get_top_n_mutation_ids(self, classifier, data, X, max_rank=100):
        if self.classifier_tag in ['LR', 'SVM', 'SVM-lin', 'XGBoost', 'CatBoost']:
            y_pred = classifier.predict_proba(X)[:, 1]
        else:
            # y_pred = np.array(classifier.predict(X), dtype=float)
            # y_pred = y_pred.flatten()
            y_pred = np.array(classifier.predict(X))

        mutant_id = data.apply(DataManager.create_mutation_id, axis=1)
        df = pd.DataFrame({'mutant_id': mutant_id, 'prediction_score': y_pred})
        df.sort_values(by=['prediction_score'], ascending=False, ignore_index=True, inplace=True)

        max_rank = min(max_rank, X.shape[0])
        return df.loc[df.index[0:max_rank], 'mutant_id'].to_numpy()

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
            mutant_id = data_long.apply(DataManager.create_mutation_id, axis=1)
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
            mutant_id_long = data_long.apply(DataManager.create_mutation_id, axis=1)
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

        if self.classifier_tag == 'CatBoost':
            clf.fit(X, y, plot=False)
        else:
            clf.fit(X, y)

        return clf

    def get_optimization_params(self):
        return self.optimization_params

    @staticmethod
    def save_classifier(classifier_tag, classifier, classifier_file):

        if classifier_tag in ['CatBoost', 'XGBoost']:
            classifier.save_model(classifier_file)
        else:
            pickle.dump(classifier, open(classifier_file, 'wb'))

    def load(self, classifier_file):
        self.classifier = ClassifierManager.load_classifier(
            self.classifier_tag, self.optimization_params, classifier_file)
        return self.classifier

    @staticmethod
    def load_classifier(classifier_tag, optimization_params, classifier_file):

        if classifier_tag in ['CatBoost', 'XGBoost']:
            classifier = optimization_params.get_base_classifier(classifier_tag)
            classifier.load_model(classifier_file)
        else:
            classifier = pickle.load(open(classifier_file, 'rb'))

        return classifier

    @staticmethod
    def write_tesla_scores(patient, dataset, nr_correct20, nr_tested20, nr_correct100, nr_immuno,
                           x_sorted, y_pred_sorted, report_file):
        idx = x_sorted['response_type'] != 'not_tested'
        y_pred_tesla = y_pred_sorted[idx].to_numpy()
        y_tesla = x_sorted.loc[idx, 'response'].to_numpy()
        ttif = nr_correct20/nr_tested20
        fr = nr_correct100/nr_immuno
        precision, recall, _ = precision_recall_curve(y_tesla, y_pred_tesla)
        auprc = auc(recall, precision)
        report_file.write("{0}\t{1}\t{2:.3f}\t{3:.3f}\t{4:.3}\n".format(dataset, patient, ttif, fr, auprc))


