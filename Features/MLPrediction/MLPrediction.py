from os.path import exists

import pandas as pd

from Utils.DataManager import *
from Classifier.PrioritizationLearner import *


class MLPrediction:

    def __init__(self, classifier, data_loader=None, classifier_tag=None):
        """ Add learned prediction scores"""

        self.classifier_file = None
        self.data_loader = data_loader
        self.classifier_tag = None
        if type(classifier) is str:  # if it's a binary file encoding a classifier
            self.classifier_file = classifier
            self.classifier_tag = os.path.basename(self.classifier_file).split('_')[0]
            self.classifier = \
                PrioritizationLearner.load_classifier(self.classifier_tag, OptimizationParams(), self.classifier_file)
        elif type(classifier) is list:
            self.classifier_file = None
            self.classifier_tag = 'Voting'
            self.classifier = classifier
        else:
            self.classifier_file = None
            self.classifier_tag = classifier_tag
            self.classifier = classifier

        return

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):
        assert self.classifier is not None, "Classifier is None"
        assert self.data_loader is not None, "DataLoader is None"

        mgr = DataManager()

        df, X, y = self.data_loader.load_patients(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            df = self.add_features_long(df, X)
        else:
            df = self.add_features_short(df, X)

        mgr.put_processed_data(df, patient, file_tag_output, peptide_type)

        if write_res:
            data = mgr.get_processed_data(patient, file_tag_input, peptide_type)
            if peptide_type == 'long':
                data = data.merge(df.loc[:, ['mutation_score', 'mutation_rank', 'peptide_id']], how='left',
                                  on='peptide_id')
            else:
                data = data.loc[data['mut_seqid'].notna(), :]
                data = data.merge(df.loc[:, ['peptide_score', 'rank_in_mutation', 'peptide_id']], how='left',
                                  left_on='peptide_id', right_on='peptide_id')

            out_file = os.path.join(Parameters().get_result_dir(), patient+"_"+peptide_type+"_"+file_tag_output+".txt")
            data.to_csv(path_or_buf=out_file, sep="\t", index=False, header=True)

        return data

    @staticmethod
    def add_rank_long(column, df):
        rel = Parameters().get_order_relation(column)
        if rel == '<':
            df.loc[:, 'mutation_rank'] = rankdata(df.loc[:, column], method='average')
        else:
            df.loc[:, 'mutation_rank'] = rankdata(-df.loc[:, column], method='average')

        return df

    def add_features_long(self, df, x):
        if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
            df.loc[:, 'mutation_score'] = self.classifier.predict_proba(x)[:, 1]
        elif self.classifier_tag == 'Voting':
            y_pred = np.full(x.shape[0], 0.0)
            for (tag, clf, w) in self.classifier:
                y_pred = np.add(y_pred, np.array(clf.predict_proba(x)[:, 1])*w)
            df.loc[:, 'mutation_score'] = y_pred
        else:
            df.loc[:, 'mutation_score'] = np.array(self.classifier.decision_function(x))

        return self.add_rank_long('mutation_score', df)

    @staticmethod
    def add_mutation_rank_short(column, df):
        ids = df['mut_seqid'].unique()
        ranks = np.full(df.shape[0], 0)
        rel = Parameters().get_order_relation(column)
        for id in ids:
            idx = df['mut_seqid'] == id
            if rel == '<':
                ranks[idx] = rankdata(df.loc[idx, column], method='average')
            else:
                ranks[idx] = rankdata(-df.loc[idx, column], method='average')

        df.loc[:, 'rank_in_mutation'] = ranks

        return df

    def add_features_short(self, df, x):
        if self.classifier_tag in ['SVM', 'SVM-lin', 'RF', 'CART', 'ADA', 'NNN', 'XGBoost', 'CatBoost', 'TabNet']:
            df.loc[:, 'peptide_score'] = self.classifier.predict_proba(x)[:, 1]
        elif self.classifier_tag == 'Voting':
            y_pred = np.full(x.shape[0], 0.0)
            for (tag, clf, w) in self.classifier:
                y_pred = np.add(y_pred, np.array(clf.predict_proba(x)[:, 1])*w)
            df.loc[:, 'peptide_score'] = y_pred
        else:
            df.loc[:, 'peptide_score'] = np.array(self.classifier.decision_function(x))

        return self.add_mutation_rank_short('peptide_score', df)

