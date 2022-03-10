from sklearn.ensemble import VotingClassifier
from sklearn.metrics import make_scorer
from Classify_functions import *
import warnings
from sklearn.exceptions import UndefinedMetricWarning
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)

global create_plots, verbose

create_plots = False
verbose = 0
nr_clfs = 1
nr_iter = 10

top100_scorer = make_scorer(nr_correct_top100, needs_threshold = True)

data_file = \
    '~/Documents/data/Prediction/Rosenberg_data/RosenbergSamplesSummary_Preds_IEDB_GTEx_Cosmic_Intogen_data_ipMSDB.txt'

df, y = load_data(data_file)

cl_vars = ['l_MUT_RI_MX', 'l_tis_exp_mean', 'FatHMM_CO', 'WT_Ovrl_I', 'rnaseq_exp_q', 'l_rnaseq_alt1']
X = df[cl_vars]
scaler = StandardScaler()
X = scaler.fit_transform(X)

print('Class counter: ' + str(Counter(y)))
patients = np.array(df.patient)

nr_correct_tot = 0
nr_immuno_tot = 0
nr_missed_tot = 0

svm_prob = np.array([])
rf_prob = np.array([])
col = np.array([])
#y = np.random.permutation(y)

feature_importance = np.zeros(len(cl_vars))
splits = train_test_split(patients, y, 0.2)
for train_idx, test_idx in splits:

    X_train = X[train_idx,]
    X_test = X[test_idx,]
    y_train = y[train_idx]
    y_test = y[test_idx]
    p_test = patients[test_idx]

    if np.sum(y_test == 1) == 0:
        continue

    best_rfs, best_scores, best_std = \
        optimize_rf_rdn(X_train, y_train, top100_scorer, create_plots, nr_clfs=nr_clfs, verbose=0,nr_iter=nr_iter)

    estimators = []
    weights = []
    for i in range(len(best_rfs)):
        estimators.append(('rf_%d' %i, best_rfs[i]))
        weights.append(best_scores[i]-best_std[i])
        best_rfs[i].fit(X_train, y_train)
        feature_importance = feature_importance + np.array(best_rfs[i].feature_importances_)

    best_svms, best_scores, best_std = \
        optimize_svm_rdn(X_train, y_train, top100_scorer, create_plots, nr_clfs=nr_clfs, verbose=0,nr_iter=nr_iter)

    for i in range(len(best_svms)):
        estimators.append(('svm_%d' %i, best_svms[i]))
        weights.append(best_scores[i]-best_std[i])


    voting_rf = VotingClassifier(estimators = estimators, voting='soft',weights=None)
    voting_rf.fit(X_train,y_train)

    best_svms[0].fit(X_train,y_train)
    best_rfs[0].fit(X_train,y_train)
    svm_prob = np.append(svm_prob, best_svms[0].predict_proba(X_test)[:,1])
    rf_prob = np.append(rf_prob, best_rfs[0].predict_proba(X_test)[:,1])
    col = np.append(col, y_test)

#    if create_plots: plot_class_results(voting_rf, X_test, y_test, 'predict_proba')
    nr_correct, nr_immuno, nr_missed = test_classifier_on_patients(voting_rf, p_test, X_test, y_test, 20,verbose=1)

    nr_correct_tot += nr_correct
    nr_immuno_tot += nr_immuno
    nr_missed_tot +=  nr_missed

    print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' %(nr_correct, nr_immuno, nr_missed))


s = sum(feature_importance)
for name, score in zip(cl_vars, feature_importance):
    print(name, score/s)

print('True positive cnt: %d, immunogenic cnt: %d, missed patient cnt: %d' % (nr_correct_tot, nr_immuno_tot, nr_missed_tot))



if not create_plots:
    plt.figure(figsize=(8, 8))

    col = plt.cm.rainbow(col/2+0.5)
    plt.scatter(svm_prob, rf_prob, alpha=0.8, s=2, c=col)
    plt.xlabel('SVM class probability')
    plt.ylabel('RF class probability')

    plt.show()


