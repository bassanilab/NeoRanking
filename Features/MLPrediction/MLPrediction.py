from os.path import exists

from Utils.DataManager import *
from Classifier.PrioritizationLearner import *


class MLPrediction:

    def __init__(self, classifier_file, data_loader):
        """ Add learned prediction scores"""

        self.classifier_file = None
        self.data_loader = data_loader
        self.classifier_tag = None
        if exists(classifier_file):
            self.classifier_file = classifier_file
            self.classifier_tag = os.path.basename(classifier_file).split('_')[0]
            self.classifier = \
                PrioritizationLearner.load_classifier(self.classifier_tag, OptimizationParams(), classifier_file)

        return

    def add_features(self, patient, file_tag_input, file_tag_output, peptide_type='long', write_res=True):
        assert self.classifier is not None, "Classifier is None"

        mgr = DataManager()

        df, X, y = self.data_loader.load_patients(patient, file_tag_input, peptide_type)

        if peptide_type == 'long':
            data = self.add_features_long(df, X)
        else:
            data = self.add_features_short(df, X)

        mgr.put_processed_data(df, patient, file_tag_output, peptide_type)

        if write_res:
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
        else:
            df.loc[:, 'peptide_score'] = np.array(self.classifier.decision_function(x))

        return self.add_mutation_rank_short('peptide_score', df)

