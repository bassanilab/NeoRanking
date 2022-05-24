from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *


class FeatureOptimizer:

    def __init__(self, classifier_tag, features, data_train, x_train, y_train, max_nr_models=3, normalizer=None,
                 scorer='sum_exp_rank', cat_dims=None, peptide_type='long', alpha=0.05, nr_iter=100, nr_cv=5,
                 shuffle=True):
        self.features = features
        self.data_train = data_train
        self.test_data = {}
        self.X_train = x_train
        self.y_train = y_train
        self.classifier_tag = classifier_tag
        self.max_nr_model = max_nr_models
        self.best_models = {}
        self.normalizer = normalizer
        self.peptide_type = peptide_type
        self.cat_dims = cat_dims
        self.scorer = scorer
        self.nr_iter = nr_iter
        self.nr_cv = nr_cv
        self.shuffle = shuffle

    def add_test_data(self, patient, data, x, y):
        self.test_data[patient] = {'data': data, 'X': x, 'y': y}

    def init_models(self):

        if self.peptide_type == 'short':
            class_ratio = sum(self.y_train == 1)/sum(self.y_train == 0)
        else:
            class_ratio = None

        cat_features = [f for f in self.X_train.columns if f in Parameters().get_categorical_features()]
        cat_idx = [self.X_train.columns.get_loc(col) for col in cat_features]

        optimizationParams = \
            OptimizationParams(self.alpha, cat_features=cat_features, cat_idx=cat_idx, cat_dims=self.cat_dims,
                               input_shape=[len(self.features)], class_ratio=class_ratio)

        learner = PrioritizationLearner(self.classifier_tag, self.scorer, optimizationParams, verbose=0,
                                        nr_iter=self.nr_iter, nr_classifiers=1, nr_cv=self.nr_cv,
                                        shuffle=self.shuffle)

        cvres, best_classifier, best_score, best_params = learner.optimize_classifier(self.X_train, self.y_train)

        for i in range(self.max_nr_model):
            self.best_models[i] = {'classifier':  best_classifier, 'params': best_params, 'score': best_score,
                                   'features': self.features, 'cat_idx': cat_idx, 'nr_trails': 0}

    def change_features_list(self):
        return
