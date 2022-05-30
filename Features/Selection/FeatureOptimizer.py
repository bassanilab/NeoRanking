import numpy as np
import shap

from DataWrangling.DataLoader import *
from Classifier.PrioritizationLearner import *
from Utils.Util_fct import *


class FeatureOptimizer:

    def __init__(self, classifier_tag, x_train, y_train, max_nr_models=3, scorer='sum_exp_rank', cat_dims=None,
                 peptide_type='long', alpha=0.05, nr_iter=100, nr_cv=5, shuffle=True, report_file=None):
        self.classifier_tag = classifier_tag
        self.features = x_train.columns.to_list()
        self.X_train = x_train
        self.y_train = y_train
        self.max_nr_model = max_nr_models
        self.peptide_type = peptide_type
        self.cat_dims = cat_dims
        self.scorer = scorer
        self.nr_iter = nr_iter
        self.nr_cv = nr_cv
        self.shuffle = shuffle
        self.alpha = alpha

        self.class_ratio = None
        self.set_class_ratio()

        self.best_models = {}
        self.report_txt = ""
        self.report_file = report_file
        self.feature_weights = {}
        self.init_models()

    def add_test_data(self, patient, data, x, y):
        self.test_data[patient] = {'data': data, 'X': x, 'y': y}

    def set_class_ratio(self):
        if self.peptide_type == 'short':
            self.class_ratio = sum(self.y_train == 1) / sum(self.y_train == 0)
        else:
            self.class_ratio = None

    def init_learner(self, features):
        cat_idx = [i for i, f in enumerate(features) if f in Parameters().get_categorical_features()]

        optimizationParams = \
            OptimizationParams(self.alpha, cat_idx=cat_idx, cat_dims=self.cat_dims, input_shape=[len(self.features)],
                               class_ratio=self.class_ratio)

        learner = \
            PrioritizationLearner(self.classifier_tag, self.scorer, optimizationParams, verbose=0, alpha=self.alpha,
                                  nr_iter=self.nr_iter, nr_classifiers=1, nr_cv=self.nr_cv, shuffle=self.shuffle)
        return learner, cat_idx

    def init_models(self):
        feature_weights = []
        for i in range(self.max_nr_model):
            self.best_models[i] = self.update_model(self.features)
            feature_weights.append(self.best_models[i]['shap_importance'])

        feature_weights = np.mean(feature_weights, axis=0)
        for i, f in enumerate(self.features):
            self.feature_weights[f] = feature_weights[i]

    def run(self, max_nr_trials=100):
        for i in range(max_nr_trials):
            for j in range(self.max_nr_model):
                self.report_txt = "Round: {0}, Model: {1}".format(i, j)
                model = self.best_models[j]
                new_features = self.change_feature_set(model)
                if new_features is not None:
                    self.best_models[j] = self.update_model(new_features, model)
                    if self.report_file is not None:
                        self.report_file.write(self.report_txt+"\n")
                        print(self.report_txt)
                    else:
                        print(self.report_txt)

        return

    def change_feature_set(self, model):
        features = model['features'].copy()
        if len(features) == 1:
            choices = ['add', 'swap']
            weights = [0.5, 0.5]
        elif len(features) == len(self.features):
            choices = ['remove']
            weights = [1.0]
        else:
            choices = ['add', 'remove', 'swap']
            weights = [0.3, 0.4, 0.3]

        mode = random.choices(choices, weights)[0]
        self.report_txt = "{0}, {1}: ".format(self.report_txt, mode)
        if mode == 'remove':
            for f in self.choose_feature(features, ascending=True, weights=model['shap_importance']):
                features.remove(f)
                self.report_txt = "{0} {1}".format(self.report_txt, f)
            return features
        elif mode == 'add':
            for f in self.choose_feature([nf for nf in self.features if nf not in features]):
                features.append(f)
                self.report_txt = "{0} {1}".format(self.report_txt, f)
            return features
        elif mode == 'swap':
            f1 = self.choose_feature([nf for nf in self.features if nf not in features])
            f2 = self.choose_feature(features, ascending=True, weights=model['shap_importance'])
            for f in f1:
                features.append(f)
                self.report_txt = "{0} +{1}".format(self.report_txt, f)
            for f in f2:
                features.remove(f)
                self.report_txt = "{0} -{1}".format(self.report_txt, f)
            return features
        else:
            return None

    def choose_feature(self, features, k=1, ascending=False, weights=None):
        if weights is None:
            weights = [self.feature_weights[f] for f in features]
        if ascending:
            s = sum(weights)
            weights = [s-w for w in weights]
        return random.choices(features, weights, k=k)

    def update_model(self, features, model=None):
        learner, cat_idx = self.init_learner(features)

        X = self.X_train[features]
        cvres, best_classifier, best_score, best_params = learner.optimize_classifier(X, self.y_train)

        if self.classifier_tag in ['XGBoost', 'CatBoost']:
            explainer = shap.Explainer(best_classifier)
        elif self.classifier_tag in ['SVM', 'SVM-lin']:
            explainer = shap.KernelExplainer(best_classifier.predict_proba, X, link="logit")
        elif self.classifier_tag in ['LR']:
            explainer = shap.Explainer(best_classifier, X, feature_names=self.features)

        shap_values = explainer(X)
        shap_values_mean = np.mean(np.abs(shap_values.values), axis=0)

        if model is None or best_score < model['score']:
            score_hist = [best_score] if model is None else model['score_hist']+[best_score]
            new_model = {'classifier': best_classifier, 'params': best_params, 'score': best_score,
                         'cnt': len(features), 'features': features, 'cat_idx': cat_idx, 'nr_trials': 0,
                         'shap_importance': shap_values_mean, 'score_hist': score_hist}
            self.report_txt = "{0}, Improved, Cnt: {1}, Features: {2}, params: {3}, score: {4}".\
                format(self.report_txt, len(features), features, best_params, best_score)

            return new_model
        else:
            model['nr_trials'] = model['nr_trials']+1
            self.report_txt = "{0}, Unchanged, Cnt: {1}, Features: {2}, params: {3}, score: {4}".\
                format(self.report_txt, len(features), features, best_params, best_score)
            return model

    def get_models(self):
        return self.best_models

    def report_models(self):
        print("")
        print("################################################################")
        for i in range(self.max_nr_model):
            if self.report_file is not None:
                self.report_file.write("Model {0}:".format(i))
                for k in self.best_models[i].keys():
                    self.report_file.write("{0}: {1}\n".format(k, str(self.best_models[i][k])))

            print("Model {0}:".format(i))
            for k in self.best_models[i].keys():
                print("{0}: {1}".format(k, str(self.best_models[i][k])))

